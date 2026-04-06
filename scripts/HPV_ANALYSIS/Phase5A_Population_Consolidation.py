#!/usr/bin/env python3
"""
Phase5A_Population_Consolidation.py
======================================
Figure 6/7 Phase 5A:
  Part A: Consolidate two-population model from all axes
  Part B: CNV / stemness deep dive by lifecycle stage
  Part C: SComatic variant profiling per population (APOBEC context, burden)
  Part D: Immune gene expression profiling per population
  Part E: Neoantigen pipeline tool probe (what's available for Phase 5B)

Inputs:
    - data/FIG_6/03_hpv16_genome/hpv16_gene_by_population.tsv
    - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv
    - data/FIG_6/02_populations/population_assignments.tsv
    - data/FIG_4/00_input/adata_final.h5ad
    - SComatic: all_samples.single_cell_genotype.filtered.tsv

Outputs (to data/FIG_6/04_population_profiles/):
    - two_population_assignments.tsv
    - population_concordance.tsv
    - scomatic_per_population.tsv
    - apobec_context_analysis.tsv
    - immune_gene_expression.tsv
    - neoantigen_tool_probe.txt
    - phase5A_diagnostic_report.txt

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
from collections import defaultdict, Counter
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, chi2_contingency
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")

# Population data from Phase 3/4
POP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/02_populations/population_assignments.tsv")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")

# SComatic variants
SCOMATIC_PATH = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"

# Cell Ranger BAM base
BAM_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Immune-related gene panels
IMMUNE_CHECKPOINT_GENES = [
    'CD274',    # PD-L1
    'PDCD1LG2', # PD-L2
    'CD80', 'CD86',  # co-stimulatory
    'CTLA4', 'HAVCR2', 'LAG3', 'TIGIT',  # T cell exhaustion (if in basal)
    'IDO1', 'LGALS9',  # immune suppression
]

ANTIGEN_PRESENTATION_GENES = [
    'HLA-A', 'HLA-B', 'HLA-C',  # MHC-I
    'B2M',   # beta-2 microglobulin (essential for MHC-I)
    'TAP1', 'TAP2',  # peptide transport
    'TAPBP',  # TAP binding protein
    'PSMB8', 'PSMB9', 'PSMB10',  # immunoproteasome
    'CALR', 'CANX',  # ER chaperones
    'ERAP1', 'ERAP2',  # aminopeptidases
]

INTERFERON_RESPONSE_GENES = [
    'IFITM1', 'IFITM2', 'IFITM3',
    'IFI6', 'IFI27', 'IFI44', 'IFI44L',
    'ISG15', 'ISG20',
    'MX1', 'MX2',
    'OAS1', 'OAS2', 'OAS3',
    'STAT1', 'STAT2', 'IRF1', 'IRF7', 'IRF9',
]

CYTOKINE_CHEMOKINE_GENES = [
    'CXCL9', 'CXCL10', 'CXCL11',  # T cell attractants
    'CCL2', 'CCL5',  # monocyte/T cell chemokines
    'IL6', 'IL8', 'IL1B',  # pro-inflammatory
    'TGFB1', 'TGFB2',  # immunosuppressive
    'VEGFA',  # angiogenic/immunosuppressive
]

NEOANTIGEN_RELEVANT_GENES = [
    'MICA', 'MICB',  # NKG2D ligands (stress signals)
    'ULBP1', 'ULBP2',  # NKG2D ligands
    'FAS', 'FASLG',  # death receptor pathway
    'TRAIL',  # TNF-related apoptosis
]

# APOBEC context: TCW motif (W = A or T)
# C>T and C>G mutations at TC[A/T] contexts are APOBEC-signature
APOBEC_CONTEXTS = {'TCA', 'TCT'}

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(msg)

def log_sep(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# =============================================================================
# PART A: CONSOLIDATE TWO-POPULATION MODEL
# =============================================================================
log_sep("PART A: Consolidate two-population model")

# Load population assignments
pop = pd.read_csv(POP_PATH, sep='\t', index_col=0)
log(f"  Population table: {pop.shape}")

# Load HPV16 gene counts
hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
log(f"  HPV16 gene table: {hpv_genes.shape}")

# Merge HPV16 gene data into population table
for col in ['E6', 'E7', 'E1', 'E2', 'E4', 'E5', 'L1', 'L2', 'URR',
            'intergenic', 'total_hpv16_genome_reads', 'early_reads',
            'late_reads', 'early_late_ratio', 'E6E7_reads', 'E6E7_fraction']:
    if col in hpv_genes.columns:
        pop[col] = hpv_genes[col].reindex(pop.index).fillna(0)

log(f"  Merged table: {pop.shape}")

# --- Define two populations using multiple concordant axes ---
# Population 1 (Mutagenic/Immunogenic): high A3A fraction, high early/late, lower CNV
# Population 2 (Stealth/CIN): low A3A fraction, low early/late, higher CNV

# Axis 1: A3A dominance (A3A fraction > median among expressed cells)
a3_expressed = pop[(pop['APOBEC3A'] > 0) | (pop['APOBEC3B'] > 0)]
a3a_frac_median = a3_expressed['A3A_fraction'].median()
pop['axis_A3A_dominant'] = (pop['A3A_fraction'] > a3a_frac_median).astype(int)
log(f"\n  Axis 1: A3A fraction > {a3a_frac_median:.4f} (median among A3-expressing cells)")
log(f"    A3A-dominant: {pop['axis_A3A_dominant'].sum()}")

# Axis 2: HPV16 lifecycle stage (early/late ratio among HPV16+ cells)
hpv_pos = pop[pop['HPV16_status'] == 'HPV16-positive']
el_median = hpv_pos.loc[hpv_pos['early_late_ratio'] > 0, 'early_late_ratio'].median()
pop['axis_early_lifecycle'] = 0
pop.loc[pop['early_late_ratio'] > el_median, 'axis_early_lifecycle'] = 1
log(f"\n  Axis 2: Early/late ratio > {el_median:.2f} (median among HPV16+)")
log(f"    Early lifecycle: {pop['axis_early_lifecycle'].sum()}")

# Axis 3: CNV status (above/below median)
cnv_median = pop['cnv_score'].median()
pop['axis_low_CNV'] = (pop['cnv_score'] <= cnv_median).astype(int)
log(f"\n  Axis 3: CNV score <= {cnv_median:.4f} (median)")
log(f"    Low CNV: {pop['axis_low_CNV'].sum()}")

# Axis 4: Differentiation (low CytoTRACE2 = more differentiated)
cyto_median = pop['CytoTRACE2_Score'].median()
pop['axis_differentiated'] = (pop['CytoTRACE2_Score'] <= cyto_median).astype(int)
log(f"\n  Axis 4: CytoTRACE2 <= {cyto_median:.4f} (median, more differentiated)")
log(f"    Differentiated: {pop['axis_differentiated'].sum()}")

# Composite score: sum of concordant axes
# Pop1 (mutagenic): high A3A, early lifecycle, low CNV, differentiated → score = sum
pop['pop1_score'] = (pop['axis_A3A_dominant'] +
                     pop['axis_early_lifecycle'] +
                     pop['axis_low_CNV'] +
                     pop['axis_differentiated'])

log(f"\n  Composite Pop1 score distribution (0-4):")
for score in range(5):
    n = (pop['pop1_score'] == score).sum()
    log(f"    Score {score}: {n:6d} ({100*n/len(pop):.1f}%)")

# Assign populations: score >= 3 → Pop1 (mutagenic), score <= 1 → Pop2 (stealth)
pop['two_pop'] = 'Intermediate'
pop.loc[pop['pop1_score'] >= 3, 'two_pop'] = 'Pop1_Mutagenic'
pop.loc[pop['pop1_score'] <= 1, 'two_pop'] = 'Pop2_Stealth'

log(f"\n  Two-population assignment:")
for p, n in pop['two_pop'].value_counts().items():
    log(f"    {p}: {n} ({100*n/len(pop):.1f}%)")

# Concordance with SBS2 groups
if 'SBS2_group' in pop.columns:
    log(f"\n  Concordance with SBS2 groups (among assigned cells):")
    assigned = pop.dropna(subset=['SBS2_group'])
    ct = pd.crosstab(assigned['SBS2_group'], assigned['two_pop'])
    log(ct.to_string())

# Profile the two populations
log(f"\n  Population profiles:")
profile_cols = ['APOBEC3A', 'APOBEC3B', 'A3A_fraction', 'cnv_score',
                'CytoTRACE2_Score', 'raw_HPV16', 'early_late_ratio',
                'sig_SBS2', 'sig_SBS13', 'n_genes']
profile_cols = [c for c in profile_cols if c in pop.columns]

log(f"  {'Metric':25s} {'Pop1_Mutagenic':>15s} {'Intermediate':>15s} {'Pop2_Stealth':>15s}")
log(f"  {'-'*25} {'-'*15} {'-'*15} {'-'*15}")
for col in profile_cols:
    vals = []
    for p in ['Pop1_Mutagenic', 'Intermediate', 'Pop2_Stealth']:
        mask = pop['two_pop'] == p
        mean = pop.loc[mask, col].mean()
        vals.append(f"{mean:.4f}")
    log(f"  {col:25s} {vals[0]:>15s} {vals[1]:>15s} {vals[2]:>15s}")

# Statistical tests between Pop1 and Pop2
log(f"\n  Mann-Whitney: Pop1_Mutagenic vs Pop2_Stealth:")
for col in profile_cols:
    v1 = pop.loc[pop['two_pop'] == 'Pop1_Mutagenic', col].dropna()
    v2 = pop.loc[pop['two_pop'] == 'Pop2_Stealth', col].dropna()
    if len(v1) > 5 and len(v2) > 5:
        u, p_val = mannwhitneyu(v1, v2, alternative='two-sided')
        log(f"  {col:25s}: p={p_val:.2e}")

# =============================================================================
# PART B: CNV / STEMNESS DEEP DIVE BY LIFECYCLE STAGE
# =============================================================================
log_sep("PART B: CNV and stemness by HPV lifecycle stage")

# Among HPV16-positive cells, define lifecycle tertiles
hpv_pos = pop[pop['HPV16_status'] == 'HPV16-positive'].copy()
el_vals = hpv_pos.loc[hpv_pos['early_late_ratio'] > 0, 'early_late_ratio']

if len(el_vals) > 30:
    el_33 = np.percentile(el_vals, 33.3)
    el_67 = np.percentile(el_vals, 66.7)
    
    hpv_pos['lifecycle_stage'] = 'no_alignment'
    hpv_pos.loc[hpv_pos['early_late_ratio'] <= el_33, 'lifecycle_stage'] = 'late_dominant'
    hpv_pos.loc[(hpv_pos['early_late_ratio'] > el_33) & (hpv_pos['early_late_ratio'] <= el_67), 'lifecycle_stage'] = 'balanced'
    hpv_pos.loc[hpv_pos['early_late_ratio'] > el_67, 'lifecycle_stage'] = 'early_dominant'
    
    # Transfer to main pop table
    pop['lifecycle_stage'] = 'not_HPV16_positive'
    pop.loc[hpv_pos.index, 'lifecycle_stage'] = hpv_pos['lifecycle_stage']
    
    log(f"  Lifecycle tertile boundaries: early/late ratio = {el_33:.2f}, {el_67:.2f}")
    log(f"\n  {'Stage':18s} {'n':>6s} {'CNV':>7s} {'CytoTR':>7s} {'A3A':>7s} {'A3B':>7s} "
        f"{'SBS2':>7s} {'E6E7%':>6s}")
    log(f"  {'-'*18} {'-'*6} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*6}")
    
    for stage in ['early_dominant', 'balanced', 'late_dominant', 'no_alignment']:
        mask = hpv_pos['lifecycle_stage'] == stage
        n = mask.sum()
        if n == 0:
            continue
        cnv = hpv_pos.loc[mask, 'cnv_score'].mean()
        cyto = hpv_pos.loc[mask, 'CytoTRACE2_Score'].mean()
        a3a = hpv_pos.loc[mask, 'APOBEC3A'].mean()
        a3b = hpv_pos.loc[mask, 'APOBEC3B'].mean()
        sbs2 = hpv_pos.loc[mask, 'sig_SBS2'].mean() if 'sig_SBS2' in hpv_pos.columns else np.nan
        e6e7 = hpv_pos.loc[mask, 'E6E7_fraction'].mean()
        
        log(f"  {stage:18s} {n:6d} {cnv:7.4f} {cyto:7.4f} {a3a:7.3f} {a3b:7.3f} "
            f"{sbs2:7.4f} {e6e7:6.4f}")
    
    # Statistical test: early_dominant vs late_dominant
    early = hpv_pos[hpv_pos['lifecycle_stage'] == 'early_dominant']
    late = hpv_pos[hpv_pos['lifecycle_stage'] == 'late_dominant']
    
    log(f"\n  Mann-Whitney: early_dominant vs late_dominant:")
    for col in ['cnv_score', 'CytoTRACE2_Score', 'APOBEC3A', 'APOBEC3B', 'A3A_fraction']:
        v1 = early[col].dropna()
        v2 = late[col].dropna()
        if len(v1) > 5 and len(v2) > 5:
            u, p_val = mannwhitneyu(v1, v2, alternative='two-sided')
            log(f"    {col:25s}: early={v1.mean():.4f}, late={v2.mean():.4f}, p={p_val:.2e}")

# cnv_leiden clusters colored by lifecycle stage
log(f"\n  cnv_leiden clusters by lifecycle stage (HPV16+ cells):")
top_clusters = hpv_pos['cnv_leiden'].value_counts().head(12).index
log(f"  {'Cluster':10s} {'n':>5s} {'%Early':>7s} {'%Balanced':>10s} {'%Late':>6s} {'%NoAln':>7s}")
log(f"  {'-'*10} {'-'*5} {'-'*7} {'-'*10} {'-'*6} {'-'*7}")
for cluster in top_clusters:
    mask = hpv_pos['cnv_leiden'] == cluster
    n = mask.sum()
    if n == 0:
        continue
    pct_early = 100 * (hpv_pos.loc[mask, 'lifecycle_stage'] == 'early_dominant').sum() / n
    pct_bal = 100 * (hpv_pos.loc[mask, 'lifecycle_stage'] == 'balanced').sum() / n
    pct_late = 100 * (hpv_pos.loc[mask, 'lifecycle_stage'] == 'late_dominant').sum() / n
    pct_no = 100 * (hpv_pos.loc[mask, 'lifecycle_stage'] == 'no_alignment').sum() / n
    log(f"  {str(cluster):10s} {n:5d} {pct_early:6.1f}% {pct_bal:9.1f}% {pct_late:5.1f}% {pct_no:6.1f}%")

# =============================================================================
# PART C: SCOMATIC VARIANT PROFILING PER POPULATION
# =============================================================================
log_sep("PART C: SComatic variant profiling per population")

if os.path.exists(SCOMATIC_PATH):
    log(f"  Loading SComatic variants: {SCOMATIC_PATH}")
    scomatic = pd.read_csv(SCOMATIC_PATH, sep='\t')
    log(f"  Total variant calls: {len(scomatic)}")
    log(f"  Columns: {list(scomatic.columns)[:10]}...")
    
    # Key columns
    log(f"\n  Key column inspection:")
    for col in ['#CHROM', 'Start', 'REF', 'ALT_expected', 'Base_observed', 'CB']:
        if col in scomatic.columns:
            log(f"    {col}: nunique={scomatic[col].nunique()}, first={scomatic[col].iloc[0]}")
    
    # Filter to actual mutations (REF != Base_observed)
    if 'Base_observed' in scomatic.columns and 'REF' in scomatic.columns:
        mutations = scomatic[scomatic['REF'] != scomatic['Base_observed']].copy()
        log(f"  Actual mutations (REF != Base_observed): {len(mutations)}")
    else:
        mutations = scomatic.copy()
    
    # Filter to basal cells (cells in our population table)
    if 'CB' in mutations.columns:
        basal_barcodes = set(pop.index)
        mutations_basal = mutations[mutations['CB'].isin(basal_barcodes)].copy()
        log(f"  Mutations in basal cells: {len(mutations_basal)}")
        
        # Add population assignment
        mutations_basal['two_pop'] = mutations_basal['CB'].map(pop['two_pop'])
        mutations_basal['HPV16_status'] = mutations_basal['CB'].map(pop['HPV16_status'])
        mutations_basal['lifecycle_stage'] = mutations_basal['CB'].map(pop.get('lifecycle_stage', pd.Series()))
        
        # --- C1: Mutation burden per population ---
        log(f"\n  --- C1: Mutation burden per cell by population ---")
        burden = mutations_basal.groupby(['CB', 'two_pop']).size().reset_index(name='n_mutations')
        
        log(f"  {'Population':20s} {'n_cells':>8s} {'mean_mut':>9s} {'median':>7s} {'max':>5s}")
        log(f"  {'-'*20} {'-'*8} {'-'*9} {'-'*7} {'-'*5}")
        for p in ['Pop1_Mutagenic', 'Intermediate', 'Pop2_Stealth']:
            b = burden[burden['two_pop'] == p]
            if len(b) > 0:
                log(f"  {p:20s} {len(b):8d} {b['n_mutations'].mean():9.2f} "
                    f"{b['n_mutations'].median():7.1f} {b['n_mutations'].max():5d}")
        
        # --- C2: APOBEC context analysis ---
        log(f"\n  --- C2: APOBEC trinucleotide context analysis ---")
        
        # Determine mutation type and context
        # We need the trinucleotide context. If not in the data, we categorize by REF>ALT
        mutations_basal['mut_type'] = mutations_basal['REF'] + '>' + mutations_basal['Base_observed']
        
        # Check for context columns
        has_context = False
        for ctx_col in ['Context', 'Trinucleotide', 'trinucleotide', 'context']:
            if ctx_col in mutations_basal.columns:
                has_context = True
                context_col = ctx_col
                break
        
        if has_context:
            log(f"  Using trinucleotide context from column: {context_col}")
            # APOBEC signature: C>T or C>G at TC[A/T] context
            mutations_basal['is_apobec'] = False
            for ctx in APOBEC_CONTEXTS:
                mask = (mutations_basal[context_col].str[0:2] == ctx[:2]) & \
                       (mutations_basal['REF'] == 'C') & \
                       (mutations_basal['Base_observed'].isin(['T', 'G']))
                mutations_basal.loc[mask, 'is_apobec'] = True
        else:
            log(f"  No trinucleotide context column found.")
            log(f"  Using C>T and C>G mutations as proxy (includes non-APOBEC C mutations)")
            mutations_basal['is_apobec'] = (
                (mutations_basal['REF'] == 'C') &
                (mutations_basal['Base_observed'].isin(['T', 'G']))
            )
        
        # APOBEC fraction per population
        log(f"\n  APOBEC-context mutation fraction by population:")
        log(f"  {'Population':20s} {'Total_mut':>10s} {'APOBEC':>8s} {'%APOBEC':>8s}")
        log(f"  {'-'*20} {'-'*10} {'-'*8} {'-'*8}")
        
        apobec_stats = []
        for p in ['Pop1_Mutagenic', 'Intermediate', 'Pop2_Stealth']:
            pmask = mutations_basal['two_pop'] == p
            total = pmask.sum()
            apobec = (pmask & mutations_basal['is_apobec']).sum()
            pct = 100 * apobec / total if total > 0 else 0
            log(f"  {p:20s} {total:10d} {apobec:8d} {pct:7.1f}%")
            apobec_stats.append({'population': p, 'total_mutations': total,
                                 'apobec_mutations': apobec, 'pct_apobec': pct})
        
        # --- C3: Mutation spectrum (6 classes) per population ---
        log(f"\n  --- C3: Mutation spectrum by population ---")
        # Collapse to 6 SBS classes (pyrimidine reference)
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        def to_pyrimidine(ref, alt):
            if ref in ['C', 'T']:
                return f"{ref}>{alt}"
            else:
                return f"{complement[ref]}>{complement[alt]}"
        
        mutations_basal['sbs_class'] = mutations_basal.apply(
            lambda x: to_pyrimidine(x['REF'], x['Base_observed'])
            if x['REF'] in 'ACGT' and x['Base_observed'] in 'ACGT' else 'other', axis=1
        )
        
        sbs_classes = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        pops_3 = ['Pop1_Mutagenic', 'Intermediate', 'Pop2_Stealth']
        
        header = f"  {'Class':>6s}" + "".join(f"  {p:>18s}" for p in pops_3)
        log(header)
        sep_line = f"  {'-'*6}" + "".join(f"  {'-'*18}" for _ in pops_3)
        log(sep_line)
        
        for sbs in sbs_classes:
            line = f"  {sbs:>6s}"
            for p in pops_3:
                pmask = mutations_basal['two_pop'] == p
                total = pmask.sum()
                sbs_count = (pmask & (mutations_basal['sbs_class'] == sbs)).sum()
                pct = 100 * sbs_count / total if total > 0 else 0
                line += f"  {pct:17.1f}%"
            log(line)
        
        # --- C4: Mutation burden per cell (by population, per-cell) ---
        log(f"\n  --- C4: Per-cell APOBEC mutation count by population ---")
        apobec_per_cell = mutations_basal[mutations_basal['is_apobec']].groupby('CB').size()
        total_per_cell = mutations_basal.groupby('CB').size()
        
        for p in ['Pop1_Mutagenic', 'Intermediate', 'Pop2_Stealth']:
            cells_in_pop = set(pop[pop['two_pop'] == p].index)
            cells_with_mut = cells_in_pop & set(total_per_cell.index)
            
            if len(cells_with_mut) > 0:
                apobec_counts = apobec_per_cell.reindex(list(cells_with_mut)).fillna(0)
                total_counts = total_per_cell.reindex(list(cells_with_mut)).fillna(0)
                log(f"  {p}: n_cells_with_mutations={len(cells_with_mut)}, "
                    f"mean_total={total_counts.mean():.2f}, "
                    f"mean_APOBEC={apobec_counts.mean():.2f}")
        
        # Save SComatic analysis
        pd.DataFrame(apobec_stats).to_csv(
            os.path.join(OUTPUT_DIR, "apobec_context_analysis.tsv"), sep='\t', index=False)
    else:
        log(f"  WARNING: No 'CB' column found in SComatic output")
else:
    log(f"  WARNING: SComatic file not found: {SCOMATIC_PATH}")

# =============================================================================
# PART D: IMMUNE GENE EXPRESSION PROFILING
# =============================================================================
log_sep("PART D: Immune gene expression profiling per population")

log(f"  Loading adata for expression extraction...")
adata = sc.read_h5ad(ADATA_PATH)

basal_mask = adata.obs['final_annotation'].str.lower().str.contains('basal')
adata_basal = adata[basal_mask].copy()
log(f"  Basal cells in adata: {adata_basal.n_obs}")

# Extract expression for all immune gene panels
all_immune_genes = {
    'checkpoint': IMMUNE_CHECKPOINT_GENES,
    'antigen_presentation': ANTIGEN_PRESENTATION_GENES,
    'interferon_response': INTERFERON_RESPONSE_GENES,
    'cytokine_chemokine': CYTOKINE_CHEMOKINE_GENES,
    'neoantigen_relevant': NEOANTIGEN_RELEVANT_GENES,
}

immune_results = []

for panel_name, gene_list in all_immune_genes.items():
    log(f"\n  --- {panel_name} panel ---")
    
    for gene in gene_list:
        # Find gene in adata
        gene_found = False
        gene_vals = None
        
        if gene in adata_basal.var_names:
            gene_idx = list(adata_basal.var_names).index(gene)
            raw = adata_basal.X[:, gene_idx]
            if hasattr(raw, 'toarray'):
                gene_vals = raw.toarray().flatten()
            else:
                gene_vals = np.array(raw).flatten()
            gene_found = True
        else:
            # Search by gene_symbol
            for var_col in ['gene_symbol', 'feature_name']:
                if var_col in adata_basal.var.columns:
                    matches = adata_basal.var[adata_basal.var[var_col] == gene]
                    if len(matches) > 0:
                        gene_idx = list(adata_basal.var_names).index(matches.index[0])
                        raw = adata_basal.X[:, gene_idx]
                        if hasattr(raw, 'toarray'):
                            gene_vals = raw.toarray().flatten()
                        else:
                            gene_vals = np.array(raw).flatten()
                        gene_found = True
                        break
        
        if gene_found:
            gene_series = pd.Series(gene_vals, index=adata_basal.obs_names)
            
            # Per-population expression
            for p in ['Pop1_Mutagenic', 'Intermediate', 'Pop2_Stealth']:
                cells = pop[pop['two_pop'] == p].index
                cells_in_adata = list(set(cells) & set(gene_series.index))
                if len(cells_in_adata) > 0:
                    vals = gene_series.loc[cells_in_adata]
                    immune_results.append({
                        'panel': panel_name, 'gene': gene, 'population': p,
                        'n_cells': len(vals), 'mean_expr': vals.mean(),
                        'pct_expressing': 100 * (vals > 0).sum() / len(vals)
                    })
            
            # Print summary
            pop1_mean = gene_series.reindex(pop[pop['two_pop']=='Pop1_Mutagenic'].index).dropna().mean()
            pop2_mean = gene_series.reindex(pop[pop['two_pop']=='Pop2_Stealth'].index).dropna().mean()
            
            # Statistical test
            v1 = gene_series.reindex(pop[pop['two_pop']=='Pop1_Mutagenic'].index).dropna()
            v2 = gene_series.reindex(pop[pop['two_pop']=='Pop2_Stealth'].index).dropna()
            if len(v1) > 10 and len(v2) > 10:
                u, p_val = mannwhitneyu(v1, v2, alternative='two-sided')
                sig = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else ''))
                log(f"    {gene:12s}: Pop1={pop1_mean:.3f}, Pop2={pop2_mean:.3f}, p={p_val:.2e} {sig}")
            else:
                log(f"    {gene:12s}: Pop1={pop1_mean:.3f}, Pop2={pop2_mean:.3f}")
        else:
            log(f"    {gene:12s}: NOT FOUND in adata")

# Save immune gene results
if immune_results:
    immune_df = pd.DataFrame(immune_results)
    immune_df.to_csv(os.path.join(OUTPUT_DIR, "immune_gene_expression.tsv"), sep='\t', index=False)
    log(f"\n  Saved: immune_gene_expression.tsv ({len(immune_df)} rows)")

# Free adata memory
del adata, adata_basal
import gc; gc.collect()

# =============================================================================
# PART E: NEOANTIGEN PIPELINE TOOL PROBE
# =============================================================================
log_sep("PART E: Neoantigen pipeline tool availability")

tools_to_check = {
    # HLA typing
    'arcasHLA': ['arcasHLA', '--help'],
    'optitype': ['OptiTypePipeline.py', '--help'],
    
    # Variant annotation
    'VEP': ['vep', '--help'],
    'SnpEff': ['snpEff', '--help'],
    'bcftools': ['bcftools', '--version'],
    
    # Neoantigen prediction
    'pVACseq': ['pvacseq', '--help'],
    'NetMHCpan': ['netMHCpan', '-h'],
    'MHCflurry': ['mhcflurry-predict', '--help'],
    
    # Fusion detection
    'STAR-Fusion': ['STAR-Fusion', '--help'],
    'arriba': ['arriba', '-h'],
    'FusionCatcher': ['fusioncatcher', '--help'],
    
    # General tools
    'samtools': ['samtools', '--version'],
    'minimap2': ['minimap2', '--version'],
    'STAR': ['STAR', '--version'],
    'bedtools': ['bedtools', '--version'],
    'picard': ['picard', '--version'],
}

available_tools = {}
log(f"  {'Tool':20s} {'Status':>10s} {'Version/Info':>30s}")
log(f"  {'-'*20} {'-'*10} {'-'*30}")

for tool_name, cmd in tools_to_check.items():
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        # Some tools return 0 on --help, some return 1
        version_info = (result.stdout + result.stderr)[:50].strip().replace('\n', ' ')
        available_tools[tool_name] = True
        log(f"  {tool_name:20s} {'AVAILABLE':>10s} {version_info:>30s}")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        available_tools[tool_name] = False
        log(f"  {tool_name:20s} {'MISSING':>10s}")

# Also check for conda-installable tools
log(f"\n  Tools that can be installed via conda (if needed):")
installable = {
    'arcasHLA': 'conda install -c bioconda arcashla',
    'VEP': 'conda install -c bioconda ensembl-vep',
    'MHCflurry': 'pip install mhcflurry && mhcflurry-downloads fetch',
    'pVACseq': 'pip install pvactools',
    'STAR-Fusion': 'conda install -c bioconda star-fusion',
    'arriba': 'conda install -c bioconda arriba',
    'bcftools': 'conda install -c bioconda bcftools',
}

for tool, install_cmd in installable.items():
    if not available_tools.get(tool, False):
        log(f"    {tool}: {install_cmd}")

# Determine viable neoantigen pipeline
log(f"\n  --- Viable neoantigen pipeline ---")

has_hla = available_tools.get('arcasHLA', False) or available_tools.get('optitype', False)
has_vep = available_tools.get('VEP', False) or available_tools.get('SnpEff', False)
has_mhc = available_tools.get('MHCflurry', False) or available_tools.get('NetMHCpan', False) or available_tools.get('pVACseq', False)
has_fusion = available_tools.get('STAR-Fusion', False) or available_tools.get('arriba', False)
has_bam_tools = available_tools.get('samtools', False)

log(f"  HLA typing:          {'YES' if has_hla else 'NEED INSTALL'}")
log(f"  Variant annotation:  {'YES' if has_vep else 'NEED INSTALL'}")
log(f"  MHC binding pred:    {'YES' if has_mhc else 'NEED INSTALL'}")
log(f"  Fusion detection:    {'YES' if has_fusion else 'NEED INSTALL'}")
log(f"  BAM manipulation:    {'YES' if has_bam_tools else 'NEED INSTALL'}")

# Cell Ranger BAM inventory for neoantigen pipeline
log(f"\n  Cell Ranger BAMs available: {len(glob.glob(os.path.join(BAM_BASE, 'SRR*', '*_S1_L001_', 'outs', 'possorted_genome_bam.bam')))}")

log(f"""
  NEOANTIGEN PIPELINE DESIGN (Phase 5B):
  
  Step 1: HLA typing per patient
    Tool: arcasHLA (preferred) or OptiType
    Input: Cell Ranger BAMs
    Output: Per-patient HLA-I alleles (HLA-A, HLA-B, HLA-C)
    
  Step 2: SComatic variants → VCF → protein consequences
    Tool: bcftools + VEP/SnpEff
    Input: SComatic TSV → convert to VCF
    Output: Annotated variants with protein-level consequences
    Filter: missense, frameshift, stop-gain per population
    
  Step 3: Mutant peptide generation
    Tool: Custom script or pVACseq
    Input: Annotated variants + reference proteome
    Output: 8-11mer peptides spanning each mutation
    
  Step 4: MHC-I binding prediction
    Tool: MHCflurry (preferred, no license needed) or NetMHCpan
    Input: Mutant peptides + patient HLA alleles
    Output: Predicted binding affinity per peptide-HLA pair
    Filter: IC50 < 500 nM (binder), < 50 nM (strong binder)
    
  Step 5: Chimeric/fusion read detection
    Tool: Custom split-read analysis from BAMs
    Input: possorted_genome_bam.bam per sample
    Method: Extract reads with SA (supplementary alignment) tag
            that map to different chromosomes or >1Mb apart
            Cross-reference with cell barcodes (CB tag)
    Output: Per-cell chimeric read counts, breakpoint positions
    
  Step 6: Neoantigen burden per population
    Compare: Pop1 vs Pop2 for:
      - Total predicted neoantigens
      - Strong binders (IC50 < 50 nM)
      - Neoantigen per mutation ratio
      - Fusion/chimeric events
      - APOBEC-context neoantigen fraction
""")

# =============================================================================
# SAVE ALL OUTPUTS
# =============================================================================
log_sep("SAVE OUTPUTS")

# Save population assignments with all axes
pop.to_csv(os.path.join(OUTPUT_DIR, "two_population_assignments.tsv"), sep='\t')
log(f"  Saved: two_population_assignments.tsv ({pop.shape})")

# Save report
report_path = os.path.join(OUTPUT_DIR, "phase5A_diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")

log_sep("PHASE 5A COMPLETE")
