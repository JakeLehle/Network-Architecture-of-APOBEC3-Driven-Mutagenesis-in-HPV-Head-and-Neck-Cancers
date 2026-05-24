#!/usr/bin/env python3
"""
Step06_Integrated_Neoantigen_Analysis.py
===========================================
Figure 7: Neoantigen Landscape - Integrated Analysis

Pulls together outputs from Steps 03-05 to address the core clinical
questions:

1. ANTIGEN LOSS MUTATIONS: Do CNV-HIGH cells have more mutations that
   destroy existing MHC-I binding (wt binds, mut doesn't)? This is the
   mirror image of neoantigen creation and directly tests immune evasion.

2. EXPRESSION SILENCING: Are SBS2-HIGH neoantigen genes downregulated
   in CNV-HIGH? Transcriptional immune escape.

3. IMMUNE ESCAPE INTEGRATION: Which SBS2-HIGH neoantigens are lost in
   CNV-HIGH through ANY mechanism (antigen loss mutation, expression
   silencing, fusion disruption)?

4. VACCINE TARGET SCORING: Rank all neoantigen genes on multiple axes
   to identify the best therapeutic candidates, considering:
   - Binding affinity and expression in each group
   - Whether the gene is essential (expressed in both groups)
   - Whether CNV-HIGH actively evades this target
   - Group specificity vs shared target potential

Inputs:
    - data/FIG_7/03_mhc_binding/{group}_all_peptide_results.tsv
    - data/FIG_7/03_mhc_binding/{group}_expression_weighted_ranking.tsv
    - data/FIG_7/03_mhc_binding/shared_neoantigen_comparison.tsv
    - data/FIG_7/04_fusion_analysis/cross_group_neoantigen_fusion_overlap.tsv
    - data/FIG_4/00_input/adata_final.h5ad
    - data/FIG_4/01_group_selection/three_group_assignments.tsv

Outputs (to data/FIG_7/05_summary/):
    - antigen_loss_per_group.tsv
    - antigen_loss_gene_detail.tsv
    - expression_silencing_candidates.tsv
    - immune_escape_integration.tsv
    - vaccine_target_scores.tsv
    - step06_integrated_report.txt

Run in NETWORK conda env (needs scanpy).

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from collections import defaultdict
from scipy.stats import mannwhitneyu, fisher_exact
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml")

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

MHC_DIR = config['outputs']['mhc_binding']
FUSION_DIR = config['outputs']['fusion_analysis']
OUTPUT_DIR = config['outputs']['neoantigen_summary']
os.makedirs(OUTPUT_DIR, exist_ok=True)

ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")

MHC_BIND_THRESH = config['parameters']['mhc_binding_threshold']  # 500 nM

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def log_sep(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# =============================================================================
# STEP 0: LOAD DATA
# =============================================================================
log_sep("STEP 0: Load all upstream outputs")

# --- Peptide results (all scored peptides, binders and non-binders) ---
peptide_results = {}
for grp in ['SBS2_HIGH', 'CNV_HIGH']:
    path = os.path.join(MHC_DIR, f"{grp}_all_peptide_results.tsv")
    if os.path.exists(path):
        peptide_results[grp] = pd.read_csv(path, sep='\t')
        log(f"  {grp} all peptide results: {len(peptide_results[grp])} rows")
    else:
        log(f"  {grp}: all_peptide_results.tsv not found")

# --- Expression-weighted rankings ---
rankings = {}
for grp in ['SBS2_HIGH', 'CNV_HIGH']:
    path = os.path.join(MHC_DIR, f"{grp}_expression_weighted_ranking.tsv")
    if os.path.exists(path):
        rankings[grp] = pd.read_csv(path, sep='\t')
        log(f"  {grp} ranking: {len(rankings[grp])} genes")

# --- Shared neoantigen comparison ---
shared_path = os.path.join(MHC_DIR, "shared_neoantigen_comparison.tsv")
shared_df = None
if os.path.exists(shared_path):
    shared_df = pd.read_csv(shared_path, sep='\t')
    log(f"  Shared neoantigen comparison: {len(shared_df)} genes")

# --- Cross-group fusion-neoantigen overlap ---
cross_fusion_path = os.path.join(FUSION_DIR, "cross_group_neoantigen_fusion_overlap.tsv")
cross_fusion_df = None
if os.path.exists(cross_fusion_path):
    cross_fusion_df = pd.read_csv(cross_fusion_path, sep='\t')
    log(f"  Cross-group fusion overlap: {len(cross_fusion_df)} entries")

# --- Within-group fusion-neoantigen crossref ---
within_fusion_path = os.path.join(FUSION_DIR, "neoantigen_fusion_crossref.tsv")
within_fusion_df = None
if os.path.exists(within_fusion_path):
    within_fusion_df = pd.read_csv(within_fusion_path, sep='\t')
    log(f"  Within-group fusion crossref: {len(within_fusion_df)} entries")

# --- Group assignments ---
groups_df = pd.read_csv(GROUP_PATH, sep='\t')
n_sbs2 = (groups_df['group'] == 'SBS2_HIGH').sum()
n_cnv = (groups_df['group'] == 'CNV_HIGH').sum()
log(f"  Groups: SBS2_HIGH={n_sbs2}, CNV_HIGH={n_cnv}")

# =============================================================================
# STEP 1: ANTIGEN LOSS MUTATIONS
# =============================================================================
log_sep("STEP 1: Antigen loss mutations")
log(f"  Antigen loss = wild-type peptide binds MHC-I (IC50 < {MHC_BIND_THRESH}nM)")
log(f"  but somatic mutation DESTROYS binding (IC50 > {MHC_BIND_THRESH}nM)")
log(f"  These are the opposite of neoantigens: mutations that help the cell hide.")

antigen_loss_summary = []
antigen_loss_details = []

for grp in ['SBS2_HIGH', 'CNV_HIGH']:
    if grp not in peptide_results:
        continue

    pep = peptide_results[grp]
    n_total = len(pep)

    # Neoantigens: mut binds, wt doesn't (mutation creates new epitope)
    neo_mask = (pep['mut_ic50'] < MHC_BIND_THRESH) & (pep['wt_ic50'] >= MHC_BIND_THRESH)
    n_neo = neo_mask.sum()

    # Antigen loss: wt binds, mut doesn't (mutation destroys existing epitope)
    loss_mask = (pep['wt_ic50'] < MHC_BIND_THRESH) & (pep['mut_ic50'] >= MHC_BIND_THRESH)
    n_loss = loss_mask.sum()

    # Both bind (mutation changes but doesn't destroy binding)
    both_bind = (pep['mut_ic50'] < MHC_BIND_THRESH) & (pep['wt_ic50'] < MHC_BIND_THRESH)
    n_both = both_bind.sum()

    # Neither binds
    neither = (pep['mut_ic50'] >= MHC_BIND_THRESH) & (pep['wt_ic50'] >= MHC_BIND_THRESH)
    n_neither = neither.sum()

    n_cells = n_sbs2 if grp == 'SBS2_HIGH' else n_cnv

    log(f"\n  {grp} (n={n_cells}):")
    log(f"    Total peptides scored:         {n_total}")
    log(f"    Neoantigen (mut<500, wt>500):  {n_neo} ({n_neo/n_cells:.2f}/cell)")
    log(f"    Antigen loss (wt<500, mut>500): {n_loss} ({n_loss/n_cells:.2f}/cell)")
    log(f"    Both bind (mut<500, wt<500):   {n_both}")
    log(f"    Neither binds:                 {n_neither}")
    log(f"    Neo:Loss ratio:                {n_neo/max(n_loss,1):.2f}")

    antigen_loss_summary.append({
        'group': grp,
        'n_cells': n_cells,
        'total_peptides': n_total,
        'neoantigens': n_neo,
        'neoantigens_per_cell': n_neo / n_cells,
        'antigen_loss': n_loss,
        'antigen_loss_per_cell': n_loss / n_cells,
        'both_bind': n_both,
        'neither_binds': n_neither,
        'neo_to_loss_ratio': n_neo / max(n_loss, 1),
    })

    # Detail: which genes have the most antigen loss events
    if n_loss > 0:
        loss_peps = pep[loss_mask].copy()
        gene_loss_counts = loss_peps['gene'].value_counts()
        log(f"\n    Top genes with antigen loss mutations:")
        for gene, count in gene_loss_counts.head(15).items():
            best_wt = loss_peps[loss_peps['gene'] == gene]['wt_ic50'].min()
            worst_mut = loss_peps[loss_peps['gene'] == gene]['mut_ic50'].min()
            antigen_loss_details.append({
                'group': grp, 'gene': gene,
                'n_loss_peptides': count,
                'best_wt_ic50': best_wt,
                'best_mut_ic50': worst_mut,
            })
            log(f"      {gene:15s}: {count} loss peptides "
                f"(best wt IC50={best_wt:.1f}, mut IC50={worst_mut:.1f})")

# Compare antigen loss rates
if len(antigen_loss_summary) == 2:
    sbs2_s = antigen_loss_summary[0]
    cnv_s = antigen_loss_summary[1]

    log(f"\n  === ANTIGEN LOSS COMPARISON ===")
    log(f"    SBS2_HIGH: {sbs2_s['antigen_loss']:.0f} loss events "
        f"({sbs2_s['antigen_loss_per_cell']:.2f}/cell)")
    log(f"    CNV_HIGH:  {cnv_s['antigen_loss']:.0f} loss events "
        f"({cnv_s['antigen_loss_per_cell']:.2f}/cell)")

    sbs2_ratio = sbs2_s['neo_to_loss_ratio']
    cnv_ratio = cnv_s['neo_to_loss_ratio']
    log(f"    SBS2 neo:loss ratio: {sbs2_ratio:.2f}")
    log(f"    CNV  neo:loss ratio: {cnv_ratio:.2f}")

    if cnv_ratio < sbs2_ratio:
        log(f"    -> CNV_HIGH has a LOWER neo:loss ratio, suggesting")
        log(f"       proportionally more antigen-destroying mutations")
    else:
        log(f"    -> SBS2_HIGH has a lower neo:loss ratio")

# Save
if antigen_loss_summary:
    pd.DataFrame(antigen_loss_summary).to_csv(
        os.path.join(OUTPUT_DIR, "antigen_loss_per_group.tsv"), sep='\t', index=False)
if antigen_loss_details:
    pd.DataFrame(antigen_loss_details).to_csv(
        os.path.join(OUTPUT_DIR, "antigen_loss_gene_detail.tsv"), sep='\t', index=False)

# =============================================================================
# STEP 2: EXPRESSION SILENCING
# =============================================================================
log_sep("STEP 2: Expression silencing of neoantigen genes")
log(f"  Genes that produce neoantigens in SBS2_HIGH but are downregulated")
log(f"  in CNV_HIGH represent transcriptional immune escape.")

if 'SBS2_HIGH' in rankings:
    sbs2_ranking = rankings['SBS2_HIGH']

    # Use cross-group expression columns from Step04
    if 'mean_expr_sbs2' in sbs2_ranking.columns and 'mean_expr_cnv' in sbs2_ranking.columns:
        sbs2_ranking = sbs2_ranking.copy()
        sbs2_ranking['expr_ratio'] = (
            sbs2_ranking['mean_expr_cnv'] / sbs2_ranking['mean_expr_sbs2'].clip(lower=0.01)
        )
        sbs2_ranking['expr_fold_change'] = np.log2(
            sbs2_ranking['mean_expr_cnv'].clip(lower=0.01) /
            sbs2_ranking['mean_expr_sbs2'].clip(lower=0.01)
        )

        # Silenced: expressed in SBS2 but <50% expression in CNV
        silenced = sbs2_ranking[
            (sbs2_ranking['mean_expr_sbs2'] > 1.0) &
            (sbs2_ranking['expr_ratio'] < 0.5)
        ].sort_values('composite_score', ascending=False)

        log(f"\n  SBS2_HIGH neoantigen genes downregulated >2-fold in CNV_HIGH:")
        log(f"  (mean_expr_sbs2 > 1.0 AND mean_expr_cnv/mean_expr_sbs2 < 0.5)")
        log(f"  Found: {len(silenced)} genes")

        if len(silenced) > 0:
            log(f"\n  {'Gene':15s}  {'SBS2_rank':>10s}  {'Composite':>10s}  "
                f"{'ExprSBS2':>9s}  {'ExprCNV':>8s}  {'FoldChg':>8s}  {'IC50':>6s}")
            log(f"  {'----':15s}  {'--------':>10s}  {'--------':>10s}  "
                f"{'-------':>9s}  {'------':>8s}  {'------':>8s}  {'----':>6s}")

            for _, row in silenced.head(25).iterrows():
                log(f"  {row['gene']:15s}  {row['rank']:10d}  "
                    f"{row['composite_score']:10.2f}  "
                    f"{row['mean_expr_sbs2']:9.4f}  {row['mean_expr_cnv']:8.4f}  "
                    f"{row['expr_fold_change']:8.2f}  {row['best_ic50']:6.1f}")

            silenced.to_csv(os.path.join(OUTPUT_DIR, "expression_silencing_candidates.tsv"),
                            sep='\t', index=False)

        # Also report upregulated (for completeness)
        upregulated = sbs2_ranking[
            (sbs2_ranking['mean_expr_sbs2'] > 1.0) &
            (sbs2_ranking['expr_ratio'] > 2.0)
        ].sort_values('composite_score', ascending=False)
        log(f"\n  SBS2 neoantigen genes upregulated >2-fold in CNV_HIGH: {len(upregulated)}")
    else:
        log(f"  Cross-group expression columns not found in ranking")

# =============================================================================
# STEP 3: IMMUNE ESCAPE INTEGRATION
# =============================================================================
log_sep("STEP 3: Immune escape integration")
log(f"  Combining all mechanisms by which CNV-HIGH cells may escape")
log(f"  SBS2-HIGH neoantigens:")
log(f"    - Antigen loss mutations (Step 1)")
log(f"    - Expression silencing (Step 2)")
log(f"    - Fusion disruption (Step05 cross-group overlap)")

if 'SBS2_HIGH' in rankings:
    sbs2_genes = set(rankings['SBS2_HIGH']['gene'])
    escape_evidence = defaultdict(lambda: {
        'antigen_loss': False, 'expression_silenced': False,
        'fusion_disrupted': False, 'mechanisms': [],
    })

    # Antigen loss genes in CNV
    if 'CNV_HIGH' in peptide_results:
        cnv_pep = peptide_results['CNV_HIGH']
        cnv_loss_mask = (cnv_pep['wt_ic50'] < MHC_BIND_THRESH) & \
                        (cnv_pep['mut_ic50'] >= MHC_BIND_THRESH)
        cnv_loss_genes = set(cnv_pep[cnv_loss_mask]['gene'].unique()) if cnv_loss_mask.any() else set()

        # SBS2 neoantigens that have antigen loss in CNV
        sbs2_neo_with_cnv_loss = sbs2_genes & cnv_loss_genes
        for gene in sbs2_neo_with_cnv_loss:
            escape_evidence[gene]['antigen_loss'] = True
            escape_evidence[gene]['mechanisms'].append('antigen_loss')

    # Expression silenced genes
    if 'silenced' in dir() and len(silenced) > 0:
        for gene in silenced['gene']:
            if gene in sbs2_genes:
                escape_evidence[gene]['expression_silenced'] = True
                escape_evidence[gene]['mechanisms'].append('expression_silenced')

    # Fusion disrupted genes (from Step05 cross-group output)
    if cross_fusion_df is not None:
        sbs2_in_cnv = cross_fusion_df[cross_fusion_df['direction'] == 'SBS2_neo_in_CNV_fusion']
        for _, row in sbs2_in_cnv.iterrows():
            gene = row['gene']
            escape_evidence[gene]['fusion_disrupted'] = True
            escape_evidence[gene]['mechanisms'].append('fusion_disrupted')

    # Build integrated table
    escape_rows = []
    for gene, evidence in escape_evidence.items():
        if gene not in sbs2_genes:
            continue
        sbs2_row = rankings['SBS2_HIGH'][rankings['SBS2_HIGH']['gene'] == gene]
        if len(sbs2_row) == 0:
            continue
        sbs2_row = sbs2_row.iloc[0]

        # Check if also neoantigen in CNV
        cnv_genes = set(rankings['CNV_HIGH']['gene']) if 'CNV_HIGH' in rankings else set()
        in_cnv = gene in cnv_genes

        escape_rows.append({
            'gene': gene,
            'sbs2_rank': sbs2_row['rank'],
            'sbs2_composite': sbs2_row['composite_score'],
            'sbs2_ic50': sbs2_row['best_ic50'],
            'sbs2_expr': sbs2_row.get('mean_expr_group', 0),
            'cnv_expr': sbs2_row.get('mean_expr_cnv', 0),
            'also_cnv_neoantigen': in_cnv,
            'antigen_loss_in_cnv': evidence['antigen_loss'],
            'expression_silenced': evidence['expression_silenced'],
            'fusion_disrupted': evidence['fusion_disrupted'],
            'n_escape_mechanisms': len(evidence['mechanisms']),
            'escape_mechanisms': '; '.join(evidence['mechanisms']),
        })

    if escape_rows:
        escape_df = pd.DataFrame(escape_rows).sort_values(
            ['n_escape_mechanisms', 'sbs2_composite'],
            ascending=[False, False]
        )
        escape_df.to_csv(os.path.join(OUTPUT_DIR, "immune_escape_integration.tsv"),
                         sep='\t', index=False)

        log(f"\n  SBS2_HIGH neoantigen genes with evidence of CNV_HIGH immune escape:")
        log(f"  Total genes with any escape evidence: {len(escape_df)}")

        multi = escape_df[escape_df['n_escape_mechanisms'] > 1]
        log(f"  Genes with MULTIPLE escape mechanisms: {len(multi)}")

        log(f"\n  {'Gene':15s}  {'SBS2rk':>7s}  {'Score':>7s}  {'Mechanisms'}")
        log(f"  {'----':15s}  {'------':>7s}  {'-----':>7s}  {'----------'}")
        for _, row in escape_df.head(25).iterrows():
            log(f"  {row['gene']:15s}  {row['sbs2_rank']:7d}  "
                f"{row['sbs2_composite']:7.1f}  {row['escape_mechanisms']}")

# =============================================================================
# STEP 4: VACCINE TARGET SCORING
# =============================================================================
log_sep("STEP 4: Vaccine target scoring")
log(f"  Scoring criteria:")
log(f"    - Binding strength (IC50)")
log(f"    - Expression in target group")
log(f"    - Expression in BOTH groups (essential gene, harder to escape)")
log(f"    - Evidence of immune pressure (escape mechanisms suggest it works)")
log(f"    - Group specificity")

vaccine_rows = []

# Collect all neoantigen genes from both groups
all_neo_genes = set()
for grp in ['SBS2_HIGH', 'CNV_HIGH']:
    if grp in rankings:
        all_neo_genes.update(rankings[grp]['gene'])

for gene in sorted(all_neo_genes):
    row = {'gene': gene}

    # SBS2 data
    if 'SBS2_HIGH' in rankings:
        sbs2_match = rankings['SBS2_HIGH'][rankings['SBS2_HIGH']['gene'] == gene]
        if len(sbs2_match) > 0:
            s = sbs2_match.iloc[0]
            row['sbs2_rank'] = s['rank']
            row['sbs2_composite'] = s['composite_score']
            row['sbs2_ic50'] = s['best_ic50']
            row['sbs2_expr'] = s.get('mean_expr_group', 0)
            row['sbs2_pct_expr'] = s.get('pct_expressing_group', 0)
            row['sbs2_n_neo'] = s.get('n_neoantigen_peptides', 0)
            row['sbs2_n_strong'] = s.get('n_strong_binders', 0)
            row['sbs2_n_diff'] = s.get('n_differential', 0)
        else:
            row['sbs2_rank'] = None

    # CNV data
    if 'CNV_HIGH' in rankings:
        cnv_match = rankings['CNV_HIGH'][rankings['CNV_HIGH']['gene'] == gene]
        if len(cnv_match) > 0:
            c = cnv_match.iloc[0]
            row['cnv_rank'] = c['rank']
            row['cnv_composite'] = c['composite_score']
            row['cnv_ic50'] = c['best_ic50']
            row['cnv_expr'] = c.get('mean_expr_group', 0)
            row['cnv_pct_expr'] = c.get('pct_expressing_group', 0)
            row['cnv_n_neo'] = c.get('n_neoantigen_peptides', 0)
            row['cnv_n_strong'] = c.get('n_strong_binders', 0)
            row['cnv_n_diff'] = c.get('n_differential', 0)
        else:
            row['cnv_rank'] = None

    # Group specificity
    in_sbs2 = row.get('sbs2_rank') is not None
    in_cnv = row.get('cnv_rank') is not None
    if in_sbs2 and in_cnv:
        row['specificity'] = 'shared'
    elif in_sbs2:
        row['specificity'] = 'SBS2_only'
    else:
        row['specificity'] = 'CNV_only'

    # Escape evidence
    if escape_rows:
        escape_match = [e for e in escape_rows if e['gene'] == gene]
        if escape_match:
            row['n_escape_mechanisms'] = escape_match[0]['n_escape_mechanisms']
            row['escape_mechanisms'] = escape_match[0]['escape_mechanisms']
        else:
            row['n_escape_mechanisms'] = 0
            row['escape_mechanisms'] = ''
    else:
        row['n_escape_mechanisms'] = 0
        row['escape_mechanisms'] = ''

    # Fusion disruption (within-group)
    if within_fusion_df is not None:
        gene_fusions = within_fusion_df[within_fusion_df['gene'] == gene]
        row['fusion_groups'] = '; '.join(gene_fusions['group'].unique()) if len(gene_fusions) > 0 else ''
    else:
        row['fusion_groups'] = ''

    # --- VACCINE SCORE ---
    # Higher = better target
    # Components:
    #   binding: best IC50 across groups -> 500/IC50 (1-100 range)
    #   expression: max expression across groups (already log1p, 0-10 range)
    #   breadth: expressed in both groups = harder to escape = bonus
    #   pressure: escape mechanisms suggest it actually works = bonus

    best_ic50 = min(
        row.get('sbs2_ic50', 999), row.get('cnv_ic50', 999)
    )
    binding_score = MHC_BIND_THRESH / max(best_ic50, 1)

    max_expr = max(row.get('sbs2_expr', 0) or 0, row.get('cnv_expr', 0) or 0)
    expr_score = max_expr

    # Breadth bonus: if expressed >1.0 in both groups, 1.5x multiplier
    sbs2_e = row.get('sbs2_expr', 0) or 0
    cnv_e = row.get('cnv_expr', 0) or 0
    breadth_mult = 1.5 if (sbs2_e > 1.0 and cnv_e > 1.0 and in_sbs2 and in_cnv) else 1.0

    # Escape pressure bonus: evidence of active evasion suggests the
    # immune system was recognizing this target
    escape_n = row.get('n_escape_mechanisms', 0)
    pressure_mult = 1.0 + (0.25 * escape_n)  # 25% bonus per mechanism

    row['binding_score'] = round(binding_score, 2)
    row['expr_score'] = round(expr_score, 4)
    row['breadth_multiplier'] = breadth_mult
    row['pressure_multiplier'] = pressure_mult
    row['vaccine_score'] = round(
        expr_score * binding_score * breadth_mult * pressure_mult, 4
    )

    vaccine_rows.append(row)

vaccine_df = pd.DataFrame(vaccine_rows)
vaccine_df = vaccine_df.sort_values('vaccine_score', ascending=False).reset_index(drop=True)
vaccine_df['vaccine_rank'] = range(1, len(vaccine_df) + 1)

vaccine_df.to_csv(os.path.join(OUTPUT_DIR, "vaccine_target_scores.tsv"), sep='\t', index=False)

# --- REPORT TOP CANDIDATES ---
log(f"\n  === TOP 30 VACCINE TARGET CANDIDATES ===")
log(f"\n  {'Vrank':>5s}  {'Gene':15s}  {'VaxScore':>9s}  {'Spec':10s}  "
    f"{'BestIC50':>9s}  {'MaxExpr':>8s}  {'Breadth':>8s}  {'Escape':>7s}  {'Mechanisms'}")
log(f"  {'-----':>5s}  {'----':15s}  {'--------':>9s}  {'----':10s}  "
    f"{'--------':>9s}  {'-------':>8s}  {'-------':>8s}  {'------':>7s}  {'----------'}")

for _, row in vaccine_df.head(30).iterrows():
    best_ic50 = min(row.get('sbs2_ic50', 999) or 999, row.get('cnv_ic50', 999) or 999)
    log(f"  {row['vaccine_rank']:5d}  {row['gene']:15s}  {row['vaccine_score']:9.2f}  "
        f"{row['specificity']:10s}  {best_ic50:9.1f}  "
        f"{row['expr_score']:8.4f}  {row['breadth_multiplier']:8.1f}  "
        f"{row['n_escape_mechanisms']:7d}  {row.get('escape_mechanisms', '')}")

# --- REPORT BY CATEGORY ---
log(f"\n  === TOP 10 SHARED TARGETS (potential dual-population vaccines) ===")
shared_vax = vaccine_df[vaccine_df['specificity'] == 'shared']
for _, row in shared_vax.head(10).iterrows():
    best_ic50 = min(row.get('sbs2_ic50', 999) or 999, row.get('cnv_ic50', 999) or 999)
    log(f"    {row['vaccine_rank']:3d}. {row['gene']:15s} "
        f"(score={row['vaccine_score']:.2f}, IC50={best_ic50:.1f}, "
        f"SBS2_rank={row.get('sbs2_rank','N/A')}, CNV_rank={row.get('cnv_rank','N/A')})")

log(f"\n  === TOP 10 SBS2-UNIQUE TARGETS (maintenance-phase specific) ===")
sbs2_vax = vaccine_df[vaccine_df['specificity'] == 'SBS2_only']
for _, row in sbs2_vax.head(10).iterrows():
    log(f"    {row['vaccine_rank']:3d}. {row['gene']:15s} "
        f"(score={row['vaccine_score']:.2f}, IC50={row.get('sbs2_ic50','N/A')}, "
        f"rank={row.get('sbs2_rank','N/A')}, escape={row['n_escape_mechanisms']})")

log(f"\n  === TOP 10 CNV-UNIQUE TARGETS (productive-phase specific) ===")
cnv_vax = vaccine_df[vaccine_df['specificity'] == 'CNV_only']
for _, row in cnv_vax.head(10).iterrows():
    log(f"    {row['vaccine_rank']:3d}. {row['gene']:15s} "
        f"(score={row['vaccine_score']:.2f}, IC50={row.get('cnv_ic50','N/A')}, "
        f"rank={row.get('cnv_rank','N/A')})")

# --- KEY TARGETS FROM DR. IPPOLITO DISCUSSION ---
log(f"\n  === KEY TARGET CHECK ===")
key_genes = ['ANXA1', 'CD74', 'HLA-A', 'HLA-B', 'HLA-C', 'DSP', 'B2M',
             'KRT5', 'KRT6A', 'MDK', 'SERPINB5']
for gene in key_genes:
    match = vaccine_df[vaccine_df['gene'] == gene]
    if len(match) > 0:
        r = match.iloc[0]
        best_ic50 = min(r.get('sbs2_ic50', 999) or 999, r.get('cnv_ic50', 999) or 999)
        log(f"    {gene:12s}: vaccine_rank={r['vaccine_rank']:3d}, "
            f"score={r['vaccine_score']:.2f}, {r['specificity']}, "
            f"IC50={best_ic50:.1f}, escape={r['n_escape_mechanisms']}")
    else:
        log(f"    {gene:12s}: not in neoantigen gene list")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("STEP 06 COMPLETE")

report_path = os.path.join(OUTPUT_DIR, "step06_integrated_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
log(f"  Vaccine targets: {os.path.join(OUTPUT_DIR, 'vaccine_target_scores.tsv')}")
