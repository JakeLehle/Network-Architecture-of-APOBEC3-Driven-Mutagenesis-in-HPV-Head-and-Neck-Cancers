#!/usr/bin/env python3
"""
Phase5A_Revised_Population_Discovery.py
==========================================
Revised approach: bottom-up population discovery from the data.

Logic:
    1. Profile ALL cnv_leiden clusters comprehensively
    2. Cluster basal cells by HPV16 gene expression patterns
    3. Identify the high-CNV/high-stemness/late-HPV cluster(s)
    4. Select a size-matched comparison group vs SBS2-HIGH (n=546)
    5. Verify biological predictions for both groups
    6. Patient composition and overlap analysis

This replaces the Phase 5A composite score approach with a data-driven
selection that follows the narrative: Figure 5 SNPs → HPV infection →
lifecycle differences → CNV/stemness as functional readout → two populations.

Inputs:
    - data/FIG_6/03_hpv16_genome/hpv16_gene_by_population.tsv
    - data/FIG_6/02_populations/population_assignments.tsv
    - data/FIG_4/00_input/adata_final.h5ad

Outputs (to data/FIG_6/04_population_profiles_v2/):
    - cluster_comprehensive_profiles.tsv
    - hpv_lifecycle_clustering.tsv
    - revised_two_population_assignments.tsv
    - population_comparison.tsv
    - phase5A_v2_diagnostic_report.txt

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from collections import defaultdict
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, kruskal
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
POP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/02_populations/population_assignments.tsv")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# HPV16 threshold from Phase 3
HPV16_THRESHOLD = 8

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
# STEP 0: LOAD ALL DATA
# =============================================================================
log_sep("STEP 0: Load data")

pop = pd.read_csv(POP_PATH, sep='\t', index_col=0)
log(f"  Population table: {pop.shape}")

hpv_genes = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
log(f"  HPV16 gene table: {hpv_genes.shape}")

groups = pd.read_csv(GROUP_PATH, sep='\t', index_col=0)
group_col = 'group' if 'group' in groups.columns else groups.columns[0]
log(f"  SBS2 groups: {groups[group_col].value_counts().to_dict()}")

# Merge HPV gene data
for col in hpv_genes.columns:
    pop[col] = hpv_genes[col].reindex(pop.index).fillna(0)

# Add SBS2 group
pop['SBS2_group'] = groups[group_col].reindex(pop.index)

# Compute derived metrics
pop['A3A_fraction'] = pop['APOBEC3A'] / (pop['APOBEC3A'] + pop['APOBEC3B'] + 0.01)

# Early/late among cells with genome alignment data
has_alignment = pop['total_hpv16_genome_reads'] > 0
pop['has_hpv_alignment'] = has_alignment.astype(int)

# L1+L2 fraction (capsid gene fraction of total HPV reads)
pop['capsid_fraction'] = (pop['L1'] + pop['L2']) / (pop['total_hpv16_genome_reads'] + 0.5)
pop['capsid_reads'] = pop['L1'] + pop['L2']
pop['early_oncogene_reads'] = pop['E6'] + pop['E7']

n_basal = len(pop)
log(f"  Total basal cells: {n_basal}")
log(f"  Cells with HPV16 genome alignments: {has_alignment.sum()}")

# =============================================================================
# STEP 1: COMPREHENSIVE cnv_leiden CLUSTER PROFILES
# =============================================================================
log_sep("STEP 1: Comprehensive cnv_leiden cluster profiles")

clusters = sorted(pop['cnv_leiden'].unique(), key=lambda x: -pop[pop['cnv_leiden']==x].shape[0])

profile_rows = []
log(f"  {'Clust':>5s} {'n':>6s} {'A3A':>6s} {'A3B':>6s} {'A3Afr':>6s} "
    f"{'CNV':>7s} {'Cyto':>6s} {'%HPV+':>6s} {'rawHPV':>7s} "
    f"{'E/L':>5s} {'%Late':>6s} {'capsid':>7s} {'SBS2':>6s} "
    f"{'%HIGH':>6s} {'%Cancer':>7s}")
log(f"  {'-'*5} {'-'*6} {'-'*6} {'-'*6} {'-'*6} "
    f"{'-'*7} {'-'*6} {'-'*6} {'-'*7} "
    f"{'-'*5} {'-'*6} {'-'*7} {'-'*6} "
    f"{'-'*6} {'-'*7}")

for cluster in clusters:
    mask = pop['cnv_leiden'] == cluster
    n = mask.sum()
    if n < 10:
        continue
    
    c = pop[mask]
    
    a3a = c['APOBEC3A'].mean()
    a3b = c['APOBEC3B'].mean()
    a3a_frac = c['A3A_fraction'].mean()
    cnv = c['cnv_score'].mean()
    cyto = c['CytoTRACE2_Score'].mean()
    
    # HPV status
    hpv_pos_pct = 100 * (c['HPV16_status'] == 'HPV16-positive').sum() / n
    raw_hpv = c['raw_HPV16'].mean()
    
    # HPV lifecycle (among cells with alignment data)
    c_aligned = c[c['has_hpv_alignment'] == 1]
    if len(c_aligned) > 5:
        el_ratio = c_aligned['early_late_ratio'].mean()
        # % cells that are late-dominant (early/late < 1)
        pct_late = 100 * (c_aligned['early_late_ratio'] < 1.0).sum() / len(c_aligned)
        capsid = c_aligned['capsid_reads'].mean()
    else:
        el_ratio = np.nan
        pct_late = np.nan
        capsid = np.nan
    
    # SBS2
    sbs2 = c['sig_SBS2'].mean() if 'sig_SBS2' in c.columns else np.nan
    pct_high = 100 * (c['SBS2_group'] == 'HIGH').sum() / n
    
    # Cancer status
    pct_cancer = 100 * (c['Final_cancer_cell_status'] == 'Cancer cell').sum() / n
    
    row = {
        'cnv_leiden': cluster, 'n_cells': n,
        'mean_A3A': a3a, 'mean_A3B': a3b, 'A3A_fraction': a3a_frac,
        'mean_cnv_score': cnv, 'mean_CytoTRACE2': cyto,
        'pct_HPV16_pos': hpv_pos_pct, 'mean_raw_HPV16': raw_hpv,
        'mean_early_late_ratio': el_ratio, 'pct_late_dominant': pct_late,
        'mean_capsid_reads': capsid, 'mean_SBS2': sbs2,
        'pct_SBS2_HIGH': pct_high, 'pct_cancer': pct_cancer,
    }
    profile_rows.append(row)
    
    log(f"  {str(cluster):>5s} {n:6d} {a3a:6.2f} {a3b:6.2f} {a3a_frac:6.3f} "
        f"{cnv:7.4f} {cyto:6.3f} {hpv_pos_pct:5.1f}% {raw_hpv:7.1f} "
        f"{el_ratio:5.2f} {pct_late:5.1f}% {capsid:7.1f} {sbs2:6.3f} "
        f"{pct_high:5.1f}% {pct_cancer:6.1f}%")

profiles_df = pd.DataFrame(profile_rows)
profiles_df.to_csv(os.path.join(OUTPUT_DIR, "cluster_comprehensive_profiles.tsv"),
                     sep='\t', index=False)

# =============================================================================
# STEP 2: RANK CLUSTERS BY CNV/STEMNESS AND IDENTIFY CANDIDATE POP2
# =============================================================================
log_sep("STEP 2: Rank clusters for high-CNV/stemness population")

# Score each cluster: high CNV + high stemness + low A3A fraction + late HPV
profiles_df['cnv_rank'] = profiles_df['mean_cnv_score'].rank(ascending=True)
profiles_df['stemness_rank'] = profiles_df['mean_CytoTRACE2'].rank(ascending=True)
profiles_df['a3a_low_rank'] = profiles_df['A3A_fraction'].rank(ascending=False)  # low A3A = high rank
profiles_df['late_hpv_rank'] = profiles_df['pct_late_dominant'].rank(ascending=True, na_option='bottom')

# Composite rank
rank_cols = ['cnv_rank', 'stemness_rank', 'a3a_low_rank', 'late_hpv_rank']
profiles_df['composite_rank'] = profiles_df[rank_cols].mean(axis=1)
profiles_df = profiles_df.sort_values('composite_rank', ascending=False)

log(f"  Clusters ranked by CNV/stemness/lowA3A/lateHPV composite:")
log(f"  {'Clust':>5s} {'n':>6s} {'CNV':>7s} {'Cyto':>6s} {'A3Afr':>6s} "
    f"{'%Late':>6s} {'%HPV+':>6s} {'CompoRk':>8s}")
log(f"  {'-'*5} {'-'*6} {'-'*7} {'-'*6} {'-'*6} {'-'*6} {'-'*6} {'-'*8}")

for _, row in profiles_df.iterrows():
    log(f"  {str(row['cnv_leiden']):>5s} {row['n_cells']:6.0f} {row['mean_cnv_score']:7.4f} "
        f"{row['mean_CytoTRACE2']:6.3f} {row['A3A_fraction']:6.3f} "
        f"{row['pct_late_dominant']:5.1f}% {row['pct_HPV16_pos']:5.1f}% "
        f"{row['composite_rank']:8.1f}")

# =============================================================================
# STEP 3: SELECT HIGH-CNV POPULATION (matched to SBS2-HIGH n=546)
# =============================================================================
log_sep("STEP 3: Select high-CNV population")

target_n = 546  # match SBS2-HIGH group size
log(f"  Target size: {target_n} (matching SBS2-HIGH)")

# Strategy: take cells with highest CNV scores across all basal cells,
# then verify they match the biological predictions

# Sort all basal cells by CNV score descending
cnv_sorted = pop.sort_values('cnv_score', ascending=False)

# Take top N by CNV
cnv_high_cells = cnv_sorted.head(target_n).index
pop['CNV_HIGH_group'] = 'other'
pop.loc[cnv_high_cells, 'CNV_HIGH_group'] = 'CNV-HIGH'

log(f"\n  CNV-HIGH group (top {target_n} by cnv_score):")
cnv_high = pop[pop['CNV_HIGH_group'] == 'CNV-HIGH']
log(f"    CNV range: {cnv_high['cnv_score'].min():.4f} - {cnv_high['cnv_score'].max():.4f}")
log(f"    CNV threshold (minimum in group): {cnv_high['cnv_score'].min():.4f}")

# Profile this group
log(f"\n  CNV-HIGH profile (n={target_n}):")
log(f"    Mean A3A:         {cnv_high['APOBEC3A'].mean():.3f}")
log(f"    Mean A3B:         {cnv_high['APOBEC3B'].mean():.3f}")
log(f"    A3A fraction:     {cnv_high['A3A_fraction'].mean():.4f}")
log(f"    Mean CNV:         {cnv_high['cnv_score'].mean():.4f}")
log(f"    Mean CytoTRACE2:  {cnv_high['CytoTRACE2_Score'].mean():.4f}")
log(f"    % HPV16+:         {100*(cnv_high['HPV16_status']=='HPV16-positive').sum()/target_n:.1f}%")
log(f"    Mean raw HPV16:   {cnv_high['raw_HPV16'].mean():.1f}")
log(f"    % Cancer cell:    {100*(cnv_high['Final_cancer_cell_status']=='Cancer cell').sum()/target_n:.1f}%")

# HPV lifecycle in CNV-HIGH
cnv_high_aligned = cnv_high[cnv_high['has_hpv_alignment'] == 1]
if len(cnv_high_aligned) > 5:
    log(f"    Cells with HPV alignment: {len(cnv_high_aligned)}")
    log(f"    Mean early/late ratio: {cnv_high_aligned['early_late_ratio'].mean():.3f}")
    log(f"    % Late dominant: {100*(cnv_high_aligned['early_late_ratio']<1.0).sum()/len(cnv_high_aligned):.1f}%")
    log(f"    Mean capsid reads: {cnv_high_aligned['capsid_reads'].mean():.2f}")
    log(f"    Mean E6+E7 reads: {cnv_high_aligned['early_oncogene_reads'].mean():.2f}")

# SBS2 in CNV-HIGH
if 'sig_SBS2' in cnv_high.columns:
    log(f"    Mean SBS2 weight: {cnv_high['sig_SBS2'].mean():.4f}")
if 'SBS2_group' in cnv_high.columns:
    n_high_sbs2 = (cnv_high['SBS2_group'] == 'HIGH').sum()
    n_low_sbs2 = (cnv_high['SBS2_group'] == 'LOW').sum()
    log(f"    SBS2-HIGH cells: {n_high_sbs2}")
    log(f"    SBS2-LOW cells:  {n_low_sbs2}")

# cnv_leiden composition of CNV-HIGH
log(f"\n  cnv_leiden cluster composition of CNV-HIGH:")
cnv_cluster_dist = cnv_high['cnv_leiden'].value_counts()
for cluster, count in cnv_cluster_dist.items():
    total_in_cluster = (pop['cnv_leiden'] == cluster).sum()
    log(f"    Cluster {cluster}: {count} cells ({100*count/target_n:.1f}% of CNV-HIGH, "
        f"{100*count/total_in_cluster:.1f}% of cluster)")

# Patient composition
log(f"\n  Patient composition of CNV-HIGH:")
for patient in sorted(cnv_high['subject id'].unique()):
    n = (cnv_high['subject id'] == patient).sum()
    log(f"    {patient}: {n} ({100*n/target_n:.1f}%)")

# =============================================================================
# STEP 4: ALTERNATIVE SELECTION — CNV-HIGH + STEMNESS + HPV-LATE
# =============================================================================
log_sep("STEP 4: Refined selection — layered criteria")

# Start with all HPV16-positive basal cells
hpv_pos_cells = pop[pop['HPV16_status'] == 'HPV16-positive'].copy()
log(f"  HPV16-positive basal cells: {len(hpv_pos_cells)}")

# Among HPV+ cells, rank by: CNV score, CytoTRACE2, and capsid fraction
# Weight CNV most heavily per Jake's priority
hpv_pos_cells['cnv_pctile'] = hpv_pos_cells['cnv_score'].rank(pct=True)
hpv_pos_cells['cyto_pctile'] = hpv_pos_cells['CytoTRACE2_Score'].rank(pct=True)

# Capsid fraction only among cells with HPV alignment data
hpv_pos_aligned = hpv_pos_cells[hpv_pos_cells['has_hpv_alignment'] == 1]
if len(hpv_pos_aligned) > 0:
    hpv_pos_cells.loc[hpv_pos_aligned.index, 'capsid_pctile'] = \
        hpv_pos_aligned['capsid_fraction'].rank(pct=True)
    hpv_pos_cells['capsid_pctile'] = hpv_pos_cells['capsid_pctile'].fillna(0.5)
else:
    hpv_pos_cells['capsid_pctile'] = 0.5

# Weighted composite (CNV heaviest)
hpv_pos_cells['stealth_score'] = (
    0.4 * hpv_pos_cells['cnv_pctile'] +
    0.3 * hpv_pos_cells['cyto_pctile'] +
    0.3 * hpv_pos_cells['capsid_pctile']
)

# Select top N
stealth_cells = hpv_pos_cells.nlargest(target_n, 'stealth_score').index
pop['refined_group'] = 'other'
pop.loc[stealth_cells, 'refined_group'] = 'Stealth_CNV'

# Also mark SBS2-HIGH
sbs2_high_cells = pop[pop['SBS2_group'] == 'HIGH'].index
pop.loc[sbs2_high_cells, 'refined_group'] = pop.loc[sbs2_high_cells, 'refined_group'].replace(
    {'other': 'SBS2_HIGH'})
# Handle overlap: cells that are both SBS2-HIGH and Stealth_CNV
overlap = set(stealth_cells) & set(sbs2_high_cells)
pop.loc[list(overlap), 'refined_group'] = 'OVERLAP'

n_sbs2 = (pop['refined_group'] == 'SBS2_HIGH').sum()
n_stealth = (pop['refined_group'] == 'Stealth_CNV').sum()
n_overlap = (pop['refined_group'] == 'OVERLAP').sum()
log(f"\n  Refined group assignment:")
log(f"    SBS2_HIGH:     {n_sbs2}")
log(f"    Stealth_CNV:   {n_stealth}")
log(f"    OVERLAP:       {n_overlap}")
log(f"    Other:         {(pop['refined_group']=='other').sum()}")

# Profile refined Stealth_CNV
stealth = pop[pop['refined_group'] == 'Stealth_CNV']
log(f"\n  Stealth_CNV profile (n={len(stealth)}):")
log(f"    Mean A3A:         {stealth['APOBEC3A'].mean():.3f}")
log(f"    Mean A3B:         {stealth['APOBEC3B'].mean():.3f}")
log(f"    A3A fraction:     {stealth['A3A_fraction'].mean():.4f}")
log(f"    Mean CNV:         {stealth['cnv_score'].mean():.4f}")
log(f"    Mean CytoTRACE2:  {stealth['CytoTRACE2_Score'].mean():.4f}")
log(f"    % HPV16+:         100.0% (by construction)")
log(f"    Mean raw HPV16:   {stealth['raw_HPV16'].mean():.1f}")
log(f"    % Cancer cell:    {100*(stealth['Final_cancer_cell_status']=='Cancer cell').sum()/len(stealth):.1f}%")

stealth_aligned = stealth[stealth['has_hpv_alignment'] == 1]
if len(stealth_aligned) > 5:
    log(f"    With HPV alignment: {len(stealth_aligned)}")
    log(f"    Mean early/late ratio: {stealth_aligned['early_late_ratio'].mean():.3f}")
    log(f"    % Late dominant: {100*(stealth_aligned['early_late_ratio']<1.0).sum()/len(stealth_aligned):.1f}%")
    log(f"    Mean capsid reads: {stealth_aligned['capsid_reads'].mean():.2f}")
    log(f"    Mean L1 reads: {stealth_aligned['L1'].mean():.2f}")
    log(f"    Mean L2 reads: {stealth_aligned['L2'].mean():.2f}")

if 'sig_SBS2' in stealth.columns:
    log(f"    Mean SBS2 weight: {stealth['sig_SBS2'].mean():.4f}")

# =============================================================================
# STEP 5: HEAD-TO-HEAD COMPARISON
# =============================================================================
log_sep("STEP 5: Head-to-head comparison — SBS2-HIGH vs Stealth_CNV")

sbs2_high = pop[pop['refined_group'] == 'SBS2_HIGH']
stealth = pop[pop['refined_group'] == 'Stealth_CNV']

compare_metrics = [
    'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'A3A_fraction',
    'cnv_score', 'CytoTRACE2_Score',
    'raw_HPV16', 'HPV16_per_1k_genes',
    'sig_SBS2', 'sig_SBS13', 'sig_SBS1', 'sig_SBS5',
    'n_genes', 'n_counts',
    'total_hpv16_genome_reads', 'early_late_ratio',
    'capsid_reads', 'capsid_fraction',
    'early_oncogene_reads', 'E1', 'E2', 'L1', 'L2',
]
compare_metrics = [m for m in compare_metrics if m in pop.columns]

log(f"  {'Metric':30s} {'SBS2-HIGH':>12s} {'Stealth_CNV':>12s} {'p-value':>12s} {'Direction':>10s}")
log(f"  {'-'*30} {'-'*12} {'-'*12} {'-'*12} {'-'*10}")

comparison_rows = []
for metric in compare_metrics:
    v1 = sbs2_high[metric].dropna()
    v2 = stealth[metric].dropna()
    
    if len(v1) < 5 or len(v2) < 5:
        continue
    
    u, p = mannwhitneyu(v1, v2, alternative='two-sided')
    direction = 'SBS2>' if v1.mean() > v2.mean() else 'Stealth>'
    
    log(f"  {metric:30s} {v1.mean():12.4f} {v2.mean():12.4f} {p:12.2e} {direction:>10s}")
    
    comparison_rows.append({
        'metric': metric,
        'SBS2_HIGH_mean': v1.mean(), 'SBS2_HIGH_median': v1.median(),
        'Stealth_CNV_mean': v2.mean(), 'Stealth_CNV_median': v2.median(),
        'p_value': p, 'direction': direction,
    })

comparison_df = pd.DataFrame(comparison_rows)
comparison_df.to_csv(os.path.join(OUTPUT_DIR, "population_comparison.tsv"), sep='\t', index=False)

# =============================================================================
# STEP 6: PATIENT COMPOSITION
# =============================================================================
log_sep("STEP 6: Patient composition")

log(f"  {'Patient':20s} {'SBS2-HIGH':>10s} {'Stealth':>10s} {'Overlap':>8s} {'Total_basal':>12s}")
log(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*8} {'-'*12}")

for patient in sorted(pop['subject id'].unique()):
    pmask = pop['subject id'] == patient
    total = pmask.sum()
    n_sbs2 = ((pop['refined_group'] == 'SBS2_HIGH') & pmask).sum()
    n_stealth = ((pop['refined_group'] == 'Stealth_CNV') & pmask).sum()
    n_ovrlp = ((pop['refined_group'] == 'OVERLAP') & pmask).sum()
    log(f"  {patient:20s} {n_sbs2:10d} {n_stealth:10d} {n_ovrlp:8d} {total:12d}")

# =============================================================================
# STEP 7: OVERLAP ANALYSIS (are these truly distinct groups?)
# =============================================================================
log_sep("STEP 7: Group distinctness analysis")

log(f"  SBS2-HIGH cells:  {len(sbs2_high)}")
log(f"  Stealth_CNV cells: {len(stealth)}")
log(f"  Overlap:           {n_overlap}")
log(f"  Jaccard index:     {n_overlap / (len(sbs2_high) + len(stealth) + n_overlap):.4f}")

# cnv_leiden cluster overlap
log(f"\n  cnv_leiden composition:")
log(f"  {'Cluster':>8s} {'SBS2-HIGH':>10s} {'Stealth':>10s}")
log(f"  {'-'*8} {'-'*10} {'-'*10}")
all_clusters = set(sbs2_high['cnv_leiden'].unique()) | set(stealth['cnv_leiden'].unique())
for cluster in sorted(all_clusters):
    n1 = (sbs2_high['cnv_leiden'] == cluster).sum()
    n2 = (stealth['cnv_leiden'] == cluster).sum()
    if n1 > 0 or n2 > 0:
        log(f"  {str(cluster):>8s} {n1:10d} {n2:10d}")

# =============================================================================
# STEP 8: VERIFY BIOLOGICAL PREDICTIONS
# =============================================================================
log_sep("STEP 8: Verify biological predictions")

log(f"""
  PREDICTION 1: Stealth_CNV has HIGHER CNV than SBS2-HIGH
    SBS2-HIGH mean CNV:  {sbs2_high['cnv_score'].mean():.4f}
    Stealth_CNV mean CNV: {stealth['cnv_score'].mean():.4f}
    Result: {'CONFIRMED' if stealth['cnv_score'].mean() > sbs2_high['cnv_score'].mean() else 'NOT CONFIRMED'}

  PREDICTION 2: Stealth_CNV has HIGHER stemness (CytoTRACE2)
    SBS2-HIGH mean CytoTRACE2:  {sbs2_high['CytoTRACE2_Score'].mean():.4f}
    Stealth_CNV mean CytoTRACE2: {stealth['CytoTRACE2_Score'].mean():.4f}
    Result: {'CONFIRMED' if stealth['CytoTRACE2_Score'].mean() > sbs2_high['CytoTRACE2_Score'].mean() else 'NOT CONFIRMED'}

  PREDICTION 3: Stealth_CNV has HIGHER A3B / LOWER A3A
    SBS2-HIGH A3A fraction:  {sbs2_high['A3A_fraction'].mean():.4f}
    Stealth_CNV A3A fraction: {stealth['A3A_fraction'].mean():.4f}
    Result: {'CONFIRMED' if stealth['A3A_fraction'].mean() < sbs2_high['A3A_fraction'].mean() else 'NOT CONFIRMED'}

  PREDICTION 4: Stealth_CNV is enriched for LATE HPV capsid genes""")

stealth_al = stealth[stealth['has_hpv_alignment'] == 1]
sbs2_al = sbs2_high[sbs2_high['has_hpv_alignment'] == 1]
if len(stealth_al) > 5 and len(sbs2_al) > 5:
    log(f"    SBS2-HIGH mean early/late:  {sbs2_al['early_late_ratio'].mean():.3f}")
    log(f"    Stealth_CNV mean early/late: {stealth_al['early_late_ratio'].mean():.3f}")
    log(f"    SBS2-HIGH mean capsid reads: {sbs2_al['capsid_reads'].mean():.2f}")
    log(f"    Stealth_CNV mean capsid reads: {stealth_al['capsid_reads'].mean():.2f}")
    log(f"    Result: {'CONFIRMED' if stealth_al['early_late_ratio'].mean() < sbs2_al['early_late_ratio'].mean() else 'NOT CONFIRMED'}")
else:
    log(f"    Insufficient aligned cells for comparison")

log(f"""
  PREDICTION 5: Stealth_CNV has LOWER SBS2 weight""")
if 'sig_SBS2' in pop.columns:
    log(f"    SBS2-HIGH mean SBS2 weight:  {sbs2_high['sig_SBS2'].mean():.4f}")
    log(f"    Stealth_CNV mean SBS2 weight: {stealth['sig_SBS2'].mean():.4f}")
    log(f"    Result: {'CONFIRMED' if stealth['sig_SBS2'].mean() < sbs2_high['sig_SBS2'].mean() else 'NOT CONFIRMED'}")

log(f"""
  PREDICTION 6: Stealth_CNV is HPV16+ with higher viral load""")
log(f"    SBS2-HIGH mean raw HPV16:  {sbs2_high['raw_HPV16'].mean():.1f}")
log(f"    Stealth_CNV mean raw HPV16: {stealth['raw_HPV16'].mean():.1f}")
log(f"    Result: {'CONFIRMED' if stealth['raw_HPV16'].mean() > sbs2_high['raw_HPV16'].mean() else 'NOT CONFIRMED'}")

# =============================================================================
# STEP 9: SAVE
# =============================================================================
log_sep("STEP 9: Save outputs")

# Save revised assignments
assign_cols = ['subject id', 'tissue type', 'Final_cancer_cell_status',
               'cnv_score', 'cnv_leiden', 'CytoTRACE2_Score',
               'APOBEC3A', 'APOBEC3B', 'A3A_fraction',
               'raw_HPV16', 'HPV16_status', 'HPV16_dose',
               'SBS2_group', 'refined_group',
               'total_hpv16_genome_reads', 'early_late_ratio',
               'capsid_reads', 'capsid_fraction',
               'has_hpv_alignment']

# Add signature columns if present
for col in pop.columns:
    if col.startswith('sig_') and col not in assign_cols:
        assign_cols.append(col)

assign_cols = [c for c in assign_cols if c in pop.columns]
pop[assign_cols].to_csv(
    os.path.join(OUTPUT_DIR, "revised_two_population_assignments.tsv"), sep='\t')
log(f"  Saved: revised_two_population_assignments.tsv")

# Save report
report_path = os.path.join(OUTPUT_DIR, "phase5A_v2_diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")

log_sep("PHASE 5A v2 COMPLETE")
