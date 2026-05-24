#!/usr/bin/env python3
"""
Diagnostic_Therapeutic_Tiers.py
==================================
Partition neoantigen targets into therapeutic tiers based on:
    - Group specificity (shared, SBS2-only, CNV-only)
    - Escape evidence in CNV (antigen loss, expression silencing, fusion disruption)
    - Dual presence (neoantigens + fusions in both groups)

Reads existing pipeline outputs from Steps 03-06. No recomputation.

Tiers:
    1A. Hot tumor: Shared neoantigens WITH escape evidence (proven immune targets)
    1B. Hot tumor: SBS2-specific neoantigens (maintenance-phase specific)
    2.  Cold tumor: CNV-specific neoantigens (only epitopes in evasive cells)
    3.  Broad coverage: Shared neoantigens, NO escape evidence, high expression
    4.  Dual neoantigen+fusion: present in both groups as neoantigen AND fusion

Run in NETWORK conda env.
"""

import os
import pandas as pd
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MHC_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/03_mhc_binding")
FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis")
SUMMARY_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/05_summary")
OUTPUT_DIR = os.path.join(SUMMARY_DIR, "THERAPEUTIC_TIERS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
# LOAD DATA
# =============================================================================
log_sep("Load pipeline outputs")

# Rankings
sbs2_rank = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_expression_weighted_ranking.tsv"), sep='\t')
cnv_rank = pd.read_csv(os.path.join(MHC_DIR, "CNV_HIGH_expression_weighted_ranking.tsv"), sep='\t')
log(f"  SBS2 ranking: {len(sbs2_rank)} genes")
log(f"  CNV ranking:  {len(cnv_rank)} genes")

sbs2_genes = set(sbs2_rank['gene'])
cnv_genes = set(cnv_rank['gene'])
shared_genes = sbs2_genes & cnv_genes
sbs2_only_genes = sbs2_genes - cnv_genes
cnv_only_genes = cnv_genes - sbs2_genes

log(f"  Shared: {len(shared_genes)}, SBS2-only: {len(sbs2_only_genes)}, CNV-only: {len(cnv_only_genes)}")

# Escape integration
escape_path = os.path.join(SUMMARY_DIR, "immune_escape_integration.tsv")
escape_df = pd.read_csv(escape_path, sep='\t') if os.path.exists(escape_path) else pd.DataFrame()
escape_genes = set(escape_df['gene']) if len(escape_df) > 0 else set()
log(f"  Genes with escape evidence: {len(escape_genes)}")

# Antigen loss detail
loss_path = os.path.join(SUMMARY_DIR, "antigen_loss_gene_detail.tsv")
loss_df = pd.read_csv(loss_path, sep='\t') if os.path.exists(loss_path) else pd.DataFrame()
cnv_loss_genes = set(loss_df[loss_df['group'] == 'CNV_HIGH']['gene']) if len(loss_df) > 0 else set()
sbs2_loss_genes = set(loss_df[loss_df['group'] == 'SBS2_HIGH']['gene']) if len(loss_df) > 0 else set()

# Expression silencing
silenced_path = os.path.join(SUMMARY_DIR, "expression_silencing_candidates.tsv")
silenced_df = pd.read_csv(silenced_path, sep='\t') if os.path.exists(silenced_path) else pd.DataFrame()
silenced_genes = set(silenced_df['gene']) if len(silenced_df) > 0 else set()
log(f"  Expression-silenced in CNV: {len(silenced_genes)}")

# Within-group fusion crossref (neoantigen genes that are also fusion partners in same group)
within_fusion_path = os.path.join(FUSION_DIR, "neoantigen_fusion_crossref.tsv")
within_df = pd.read_csv(within_fusion_path, sep='\t') if os.path.exists(within_fusion_path) else pd.DataFrame()
sbs2_neo_fusion = set(within_df[within_df['group'] == 'SBS2_HIGH']['gene']) if len(within_df) > 0 else set()
cnv_neo_fusion = set(within_df[within_df['group'] == 'CNV_HIGH']['gene']) if len(within_df) > 0 else set()
log(f"  SBS2 neoantigen+fusion genes: {len(sbs2_neo_fusion)}")
log(f"  CNV neoantigen+fusion genes:  {len(cnv_neo_fusion)}")

# Cross-group fusion overlap (SBS2 neoantigens disrupted by CNV fusions)
cross_path = os.path.join(FUSION_DIR, "cross_group_neoantigen_fusion_overlap.tsv")
cross_df = pd.read_csv(cross_path, sep='\t') if os.path.exists(cross_path) else pd.DataFrame()
sbs2_neo_cnv_fusion = set(cross_df[cross_df['direction'] == 'SBS2_neo_in_CNV_fusion']['gene']) if len(cross_df) > 0 else set()
cnv_neo_sbs2_fusion = set(cross_df[cross_df['direction'] == 'CNV_neo_in_SBS2_fusion']['gene']) if len(cross_df) > 0 else set()
log(f"  SBS2 neoantigens in CNV fusions: {len(sbs2_neo_cnv_fusion)}")
log(f"  CNV neoantigens in SBS2 fusions: {len(cnv_neo_sbs2_fusion)}")

# Exclusive fusion gene lists
sbs2_excl_fusions_path = os.path.join(FUSION_DIR, "SBS2_HIGH_exclusive_pairs.tsv")
cnv_excl_fusions_path = os.path.join(FUSION_DIR, "CNV_HIGH_exclusive_pairs.tsv")
sbs2_fusion_genes_all = set()
cnv_fusion_genes_all = set()
if os.path.exists(sbs2_excl_fusions_path):
    sf = pd.read_csv(sbs2_excl_fusions_path, sep='\t')
    sbs2_fusion_genes_all = set(sf['geneA']) | set(sf['geneB'])
if os.path.exists(cnv_excl_fusions_path):
    cf = pd.read_csv(cnv_excl_fusions_path, sep='\t')
    cnv_fusion_genes_all = set(cf['geneA']) | set(cf['geneB'])
log(f"  SBS2 exclusive fusion genes (all): {len(sbs2_fusion_genes_all)}")
log(f"  CNV exclusive fusion genes (all):  {len(cnv_fusion_genes_all)}")

# Vaccine scores
vaccine_path = os.path.join(SUMMARY_DIR, "vaccine_target_scores.tsv")
vaccine_df = pd.read_csv(vaccine_path, sep='\t') if os.path.exists(vaccine_path) else pd.DataFrame()

# =============================================================================
# BUILD COMPREHENSIVE ESCAPE EVIDENCE PER GENE
# =============================================================================
log_sep("Build per-gene escape evidence")

gene_escape = defaultdict(lambda: {
    'antigen_loss_sbs2': False, 'antigen_loss_cnv': False,
    'expression_silenced': False,
    'fusion_in_sbs2': False, 'fusion_in_cnv': False,
    'cross_sbs2_neo_cnv_fusion': False, 'cross_cnv_neo_sbs2_fusion': False,
    'mechanisms': [],
})

for gene in sbs2_loss_genes:
    gene_escape[gene]['antigen_loss_sbs2'] = True
for gene in cnv_loss_genes:
    gene_escape[gene]['antigen_loss_cnv'] = True
for gene in silenced_genes:
    gene_escape[gene]['expression_silenced'] = True
for gene in sbs2_fusion_genes_all:
    gene_escape[gene]['fusion_in_sbs2'] = True
for gene in cnv_fusion_genes_all:
    gene_escape[gene]['fusion_in_cnv'] = True
for gene in sbs2_neo_cnv_fusion:
    gene_escape[gene]['cross_sbs2_neo_cnv_fusion'] = True
for gene in cnv_neo_sbs2_fusion:
    gene_escape[gene]['cross_cnv_neo_sbs2_fusion'] = True

# =============================================================================
# TIER ASSIGNMENT
# =============================================================================
log_sep("Tier assignment")

tier_rows = []

for gene in sorted(sbs2_genes | cnv_genes):
    in_sbs2 = gene in sbs2_genes
    in_cnv = gene in cnv_genes
    esc = gene_escape[gene]

    # Get ranking data
    sbs2_data = sbs2_rank[sbs2_rank['gene'] == gene].iloc[0] if in_sbs2 else None
    cnv_data = cnv_rank[cnv_rank['gene'] == gene].iloc[0] if in_cnv else None

    # CNV-specific escape evidence (relevant for shared and SBS2 genes)
    has_cnv_escape = (
        esc['antigen_loss_cnv'] or
        esc['expression_silenced'] or
        esc['cross_sbs2_neo_cnv_fusion']
    )

    # Any fusion in either group
    has_fusion_sbs2 = esc['fusion_in_sbs2'] or (gene in sbs2_neo_fusion)
    has_fusion_cnv = esc['fusion_in_cnv'] or (gene in cnv_neo_fusion)
    has_fusion_both = has_fusion_sbs2 and has_fusion_cnv

    # Tier assignment
    if in_sbs2 and in_cnv:
        if has_cnv_escape:
            tier = '1A_hot_shared_escaped'
        elif has_fusion_both:
            tier = '4_dual_neo_fusion'
        else:
            tier = '3_broad_coverage'
    elif in_sbs2 and not in_cnv:
        tier = '1B_hot_sbs2_specific'
    elif in_cnv and not in_sbs2:
        tier = '2_cold_cnv_specific'
    else:
        tier = 'unclassified'

    row = {
        'gene': gene,
        'tier': tier,
        'in_sbs2_neo': in_sbs2,
        'in_cnv_neo': in_cnv,
        'sbs2_rank': sbs2_data['rank'] if sbs2_data is not None else None,
        'sbs2_composite': sbs2_data['composite_score'] if sbs2_data is not None else None,
        'sbs2_ic50': sbs2_data['best_ic50'] if sbs2_data is not None else None,
        'sbs2_expr': sbs2_data.get('mean_expr_group', None) if sbs2_data is not None else None,
        'sbs2_pct_expr': sbs2_data.get('pct_expressing_group', None) if sbs2_data is not None else None,
        'cnv_rank': cnv_data['rank'] if cnv_data is not None else None,
        'cnv_composite': cnv_data['composite_score'] if cnv_data is not None else None,
        'cnv_ic50': cnv_data['best_ic50'] if cnv_data is not None else None,
        'cnv_expr': cnv_data.get('mean_expr_group', None) if cnv_data is not None else None,
        'cnv_pct_expr': cnv_data.get('pct_expressing_group', None) if cnv_data is not None else None,
        'antigen_loss_cnv': esc['antigen_loss_cnv'],
        'expression_silenced': esc['expression_silenced'],
        'fusion_in_sbs2': has_fusion_sbs2,
        'fusion_in_cnv': has_fusion_cnv,
        'cross_sbs2_neo_cnv_fusion': esc['cross_sbs2_neo_cnv_fusion'],
        'cross_cnv_neo_sbs2_fusion': esc['cross_cnv_neo_sbs2_fusion'],
        'has_cnv_escape': has_cnv_escape,
    }

    # Best IC50 across groups
    ic50s = [v for v in [row['sbs2_ic50'], row['cnv_ic50']] if v is not None]
    row['best_ic50'] = min(ic50s) if ic50s else None

    # Max expression
    exprs = [v for v in [row['sbs2_expr'], row['cnv_expr']] if v is not None]
    row['max_expr'] = max(exprs) if exprs else 0

    tier_rows.append(row)

tier_df = pd.DataFrame(tier_rows)

# Sort within each tier by composite score
def tier_sort_key(row):
    if row['tier'] in ['1A_hot_shared_escaped', '3_broad_coverage', '4_dual_neo_fusion']:
        return row.get('sbs2_composite', 0) or 0
    elif row['tier'] == '1B_hot_sbs2_specific':
        return row.get('sbs2_composite', 0) or 0
    elif row['tier'] == '2_cold_cnv_specific':
        return row.get('cnv_composite', 0) or 0
    return 0

tier_df['sort_score'] = tier_df.apply(tier_sort_key, axis=1)
tier_df = tier_df.sort_values(['tier', 'sort_score'], ascending=[True, False])
tier_df.to_csv(os.path.join(OUTPUT_DIR, "all_genes_tiered.tsv"), sep='\t', index=False)

# =============================================================================
# REPORT EACH TIER
# =============================================================================

tier_counts = tier_df['tier'].value_counts()
log(f"\n  Tier distribution:")
for tier, count in sorted(tier_counts.items()):
    log(f"    {tier}: {count} genes")

# --- TIER 1A: Hot tumor - shared with escape ---
log_sep("TIER 1A: Hot tumor targets (shared neoantigens WITH CNV escape)")
log(f"  These neoantigens provoked immune recognition during maintenance phase.")
log(f"  CNV-HIGH cells actively evade them. Best for SBS2-dominant hot tumors.")

t1a = tier_df[tier_df['tier'] == '1A_hot_shared_escaped'].head(25)
if len(t1a) > 0:
    log(f"\n  {'Gene':15s}  {'SBS2rk':>7s}  {'SBS2sc':>7s}  {'IC50':>6s}  "
        f"{'SBS2expr':>9s}  {'CNVexpr':>8s}  {'Escape evidence'}")
    log(f"  {'----':15s}  {'------':>7s}  {'------':>7s}  {'----':>6s}  "
        f"{'--------':>9s}  {'-------':>8s}  {'---------------'}")
    for _, r in t1a.iterrows():
        esc_list = []
        if r['antigen_loss_cnv']: esc_list.append('antigen_loss')
        if r['expression_silenced']: esc_list.append('silenced')
        if r['cross_sbs2_neo_cnv_fusion']: esc_list.append('CNV_fusion')
        log(f"  {r['gene']:15s}  {r['sbs2_rank']:7.0f}  {r['sbs2_composite']:7.1f}  "
            f"{r['best_ic50']:6.1f}  {r['sbs2_expr']:9.4f}  "
            f"{r['cnv_expr'] if r['cnv_expr'] else 0:8.4f}  "
            f"{', '.join(esc_list)}")

pd.DataFrame(t1a).to_csv(os.path.join(OUTPUT_DIR, "tier1A_hot_shared_escaped.tsv"), sep='\t', index=False)

# --- TIER 1B: Hot tumor - SBS2 specific ---
log_sep("TIER 1B: Hot tumor targets (SBS2-specific neoantigens)")
log(f"  Only present in maintenance-phase cells. Strongest immune stimulation.")

t1b = tier_df[tier_df['tier'] == '1B_hot_sbs2_specific'].head(25)
if len(t1b) > 0:
    log(f"\n  {'Gene':15s}  {'SBS2rk':>7s}  {'SBS2sc':>7s}  {'IC50':>6s}  "
        f"{'SBS2expr':>9s}  {'%Expr':>6s}  {'#Neo':>5s}  {'#Strong':>7s}")
    log(f"  {'----':15s}  {'------':>7s}  {'------':>7s}  {'----':>6s}  "
        f"{'--------':>9s}  {'-----':>6s}  {'----':>5s}  {'------':>7s}")
    for _, r in t1b.iterrows():
        s = sbs2_rank[sbs2_rank['gene'] == r['gene']]
        n_neo = s.iloc[0].get('n_neoantigen_peptides', '') if len(s) > 0 else ''
        n_strong = s.iloc[0].get('n_strong_binders', '') if len(s) > 0 else ''
        log(f"  {r['gene']:15s}  {r['sbs2_rank']:7.0f}  {r['sbs2_composite']:7.1f}  "
            f"{r['best_ic50']:6.1f}  {r['sbs2_expr']:9.4f}  "
            f"{r['sbs2_pct_expr']:5.1f}%  {n_neo:>5}  {n_strong:>7}")

pd.DataFrame(t1b).to_csv(os.path.join(OUTPUT_DIR, "tier1B_hot_sbs2_specific.tsv"), sep='\t', index=False)

# --- TIER 2: Cold tumor - CNV specific ---
log_sep("TIER 2: Cold tumor targets (CNV-specific neoantigens)")
log(f"  Only epitopes present in productive-phase immune-evasive cells.")
log(f"  Require combination with checkpoint blockade.")

t2 = tier_df[tier_df['tier'] == '2_cold_cnv_specific'].head(25)
if len(t2) > 0:
    log(f"\n  {'Gene':15s}  {'CNVrk':>6s}  {'CNVsc':>6s}  {'IC50':>6s}  "
        f"{'CNVexpr':>8s}  {'%Expr':>6s}  {'#Neo':>5s}  {'#Strong':>7s}")
    log(f"  {'----':15s}  {'-----':>6s}  {'-----':>6s}  {'----':>6s}  "
        f"{'-------':>8s}  {'-----':>6s}  {'----':>5s}  {'------':>7s}")
    for _, r in t2.iterrows():
        c = cnv_rank[cnv_rank['gene'] == r['gene']]
        n_neo = c.iloc[0].get('n_neoantigen_peptides', '') if len(c) > 0 else ''
        n_strong = c.iloc[0].get('n_strong_binders', '') if len(c) > 0 else ''
        log(f"  {r['gene']:15s}  {r['cnv_rank']:6.0f}  {r['cnv_composite']:6.1f}  "
            f"{r['best_ic50']:6.1f}  {r['cnv_expr']:8.4f}  "
            f"{r['cnv_pct_expr']:5.1f}%  {n_neo:>5}  {n_strong:>7}")

pd.DataFrame(t2).to_csv(os.path.join(OUTPUT_DIR, "tier2_cold_cnv_specific.tsv"), sep='\t', index=False)

# --- TIER 3: Broad coverage - shared, no escape ---
log_sep("TIER 3: Broad coverage targets (shared, NO escape evidence)")
log(f"  Present in both populations with no evidence of immune selection.")
log(f"  Likely less immunogenic but provide baseline coverage across all tumors.")
log(f"  Best candidates for combination with Tier 1/2 targets.")

t3 = tier_df[tier_df['tier'] == '3_broad_coverage'].head(25)
if len(t3) > 0:
    log(f"\n  {'Gene':15s}  {'SBS2rk':>7s}  {'CNVrk':>7s}  {'IC50':>6s}  "
        f"{'SBS2expr':>9s}  {'CNVexpr':>8s}  {'SBS2sc':>7s}  {'CNVsc':>7s}")
    log(f"  {'----':15s}  {'------':>7s}  {'------':>7s}  {'----':>6s}  "
        f"{'--------':>9s}  {'-------':>8s}  {'------':>7s}  {'------':>7s}")
    for _, r in t3.iterrows():
        log(f"  {r['gene']:15s}  {r['sbs2_rank']:7.0f}  {r['cnv_rank']:7.0f}  "
            f"{r['best_ic50']:6.1f}  {r['sbs2_expr']:9.4f}  "
            f"{r['cnv_expr'] if r['cnv_expr'] else 0:8.4f}  "
            f"{r['sbs2_composite']:7.1f}  {r['cnv_composite']:7.1f}")

pd.DataFrame(t3).to_csv(os.path.join(OUTPUT_DIR, "tier3_broad_coverage.tsv"), sep='\t', index=False)

# --- TIER 4: Dual neoantigen + fusion ---
log_sep("TIER 4: Dual neoantigen + fusion (both groups)")
log(f"  Genes with neoantigens in both groups AND fusion events in both groups.")
log(f"  Fusion products may create additional novel epitopes.")

t4 = tier_df[tier_df['tier'] == '4_dual_neo_fusion']
if len(t4) > 0:
    log(f"\n  Found {len(t4)} genes:")
    for _, r in t4.iterrows():
        log(f"    {r['gene']:15s}: SBS2_rank={r['sbs2_rank']:.0f}, "
            f"CNV_rank={r['cnv_rank']:.0f}, "
            f"IC50={r['best_ic50']:.1f}, "
            f"SBS2_fusion={r['fusion_in_sbs2']}, CNV_fusion={r['fusion_in_cnv']}")
else:
    log(f"\n  No genes found with neoantigens AND fusions in BOTH groups")

    # Check near-misses: genes with neoantigen in both + fusion in at least one
    near_miss = tier_df[
        (tier_df['in_sbs2_neo']) & (tier_df['in_cnv_neo']) &
        (tier_df['fusion_in_sbs2'] | tier_df['fusion_in_cnv'])
    ]
    if len(near_miss) > 0:
        log(f"\n  Near-misses (neoantigen in both + fusion in one group): {len(near_miss)}")
        for _, r in near_miss.iterrows():
            fusion_where = []
            if r['fusion_in_sbs2']: fusion_where.append('SBS2')
            if r['fusion_in_cnv']: fusion_where.append('CNV')
            log(f"    {r['gene']:15s}: SBS2_rank={r['sbs2_rank']:.0f}, "
                f"CNV_rank={r['cnv_rank']:.0f}, "
                f"fusion in: {', '.join(fusion_where)}")

pd.DataFrame(t4).to_csv(os.path.join(OUTPUT_DIR, "tier4_dual_neo_fusion.tsv"), sep='\t', index=False)

# =============================================================================
# SUMMARY TABLE
# =============================================================================
log_sep("Summary")

log(f"\n  === THERAPEUTIC TIER SUMMARY ===\n")
log(f"  Tier 1A (hot, shared+escaped):  {len(tier_df[tier_df['tier'] == '1A_hot_shared_escaped'])} genes")
log(f"  Tier 1B (hot, SBS2-specific):   {len(tier_df[tier_df['tier'] == '1B_hot_sbs2_specific'])} genes")
log(f"  Tier 2  (cold, CNV-specific):   {len(tier_df[tier_df['tier'] == '2_cold_cnv_specific'])} genes")
log(f"  Tier 3  (broad, no escape):     {len(tier_df[tier_df['tier'] == '3_broad_coverage'])} genes")
log(f"  Tier 4  (dual neo+fusion):      {len(tier_df[tier_df['tier'] == '4_dual_neo_fusion'])} genes")

log(f"\n  === CLINICAL STRATEGY ===\n")
log(f"  Hot tumor (high SBS2:CNV ratio, maintenance-phase dominant):")
log(f"    Primary: Tier 1B targets (ANXA1, S100A8, etc.)")
log(f"    Support: Tier 1A targets (proven immune recognition)")
log(f"    Backbone: Tier 3 targets (broad coverage)")
log(f"")
log(f"  Cold tumor (high CNV:SBS2 ratio, productive-phase dominant):")
log(f"    Primary: Tier 2 targets + checkpoint blockade")
log(f"    Backbone: Tier 3 targets (broad coverage)")
log(f"    Note: Tier 1A/1B targets likely ineffective (already evaded)")
log(f"")
log(f"  Mixed tumor:")
log(f"    Combination: Tier 1A + Tier 2 + Tier 3 backbone")

# Save report
report_path = os.path.join(OUTPUT_DIR, "therapeutic_tier_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")
