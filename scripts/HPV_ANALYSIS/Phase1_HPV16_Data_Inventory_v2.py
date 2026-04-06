#!/usr/bin/env python3
"""
Phase1_HPV16_Data_Inventory_v2.py
==================================
Figure 6/7 Preparatory Analysis — Data Inventory and HPV16 Integration (CORRECTED)

Fixes from v1:
    - Signature weights transposed (file is signatures x cells, need cells x signatures)
    - HPV16 counts transferred from adata_v_pp.h5ad (not present in adata_final)
    - hpv state metadata confirmed all NaN, removed
    - Master table limited to key columns (not 31k signature columns)

Inputs:
    - data/FIG_4/00_input/adata_final.h5ad
    - data/FIG_6/00_input/adata_v_pp.h5ad           (combined expression + viral, log-normalized)
    - data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv
    - data/FIG_4/00_input/signature_weights_per_cell.txt  (TRANSPOSED: signatures x cells)

Outputs (to data/FIG_6/00_diagnostic_inventory/):
    - basal_cell_master_table.tsv
    - diagnostic_report.txt
    - HPV16_distribution_by_patient.tsv
    - HPV16_x_SBS2_crosstab.tsv
    - HPV16_x_cancer_status_crosstab.tsv
    - HPV16_x_cnv_leiden_crosstab.tsv
    - subgroup_summary_stats.tsv

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import fisher_exact, chi2_contingency, mannwhitneyu, spearmanr
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input paths
ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
ADATA_VIRAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/00_input/adata_v_pp.h5ad")
GROUP_ASSIGNMENTS_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv")
SIGNATURE_WEIGHTS_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")

# Output directory
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/00_diagnostic_inventory")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# A3 family genes
A3_GENES = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H']

# Key signatures to report
KEY_SIGNATURES = ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS13']

# =============================================================================
# LOGGING HELPER
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg)
    report_lines.append(msg)

def log_sep(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# =============================================================================
# STEP 1: LOAD ADATA_FINAL
# =============================================================================
log_sep("STEP 1: Load adata_final.h5ad")
adata = sc.read_h5ad(ADATA_FINAL_PATH)
log(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

# Cell type distribution
ct_counts = adata.obs['final_annotation'].value_counts()
log("  Cell type distribution:")
for ct, count in ct_counts.items():
    log(f"    {ct:40s} {count:6d}  ({100*count/adata.n_obs:.1f}%)")

# Identify basal cells
basal_mask = adata.obs['final_annotation'].str.lower().str.contains('basal')
basal_barcodes = adata.obs_names[basal_mask]
n_basal = len(basal_barcodes)
log(f"\n  Basal cells: {n_basal}")

# =============================================================================
# STEP 2: TRANSFER HPV16 COUNTS FROM VIRAL ADATA
# =============================================================================
log_sep("STEP 2: Transfer HPV16 counts from adata_v_pp.h5ad")

adata_v = sc.read_h5ad(ADATA_VIRAL_PATH)
log(f"  Viral adata: {adata_v.n_obs} cells, {adata_v.n_vars} features")

# Barcode overlap
overlap = set(adata.obs_names) & set(adata_v.obs_names)
log(f"  Barcode overlap with adata_final: {len(overlap)} / {adata_v.n_obs}")

# Extract HPV16 counts
HPV16_FEATURE = 'Human papillomavirus 16'
assert HPV16_FEATURE in adata_v.var_names, f"'{HPV16_FEATURE}' not found in viral adata!"

hpv16_raw = adata_v[:, HPV16_FEATURE].X
if hasattr(hpv16_raw, 'toarray'):
    hpv16_vals = hpv16_raw.toarray().flatten()
else:
    hpv16_vals = np.array(hpv16_raw).flatten()

hpv16_series = pd.Series(hpv16_vals, index=adata_v.obs_names, name='HPV16_reads')
n_pos_all = (hpv16_vals > 0).sum()
log(f"  HPV16+ cells (all types): {n_pos_all} / {adata_v.n_obs} ({100*n_pos_all/adata_v.n_obs:.1f}%)")
log(f"  NOTE: Values are log-normalized (CPM + log1p). Max = {hpv16_vals.max():.2f}")

# Also check a few other notable viruses for context
log("\n  Other notable viral features:")
for virus in ['Human papillomavirus 154', 'Human papillomavirus KC5',
              'human papillomavirus 18', 'Gammapapillomavirus sp.']:
    if virus in adata_v.var_names:
        v = adata_v[:, virus].X
        if hasattr(v, 'toarray'):
            v = v.toarray().flatten()
        else:
            v = np.array(v).flatten()
        n = (v > 0).sum()
        if n > 0:
            log(f"    {virus}: {n} positive cells")

# Free viral adata memory
del adata_v
import gc; gc.collect()

# =============================================================================
# STEP 3: LOAD SIGNATURE WEIGHTS (TRANSPOSED)
# =============================================================================
log_sep("STEP 3: Load signature weights (transposed)")

sig_raw = pd.read_csv(SIGNATURE_WEIGHTS_PATH, sep='\t', index_col=0)
log(f"  Raw shape: {sig_raw.shape} (signatures x cells)")
log(f"  Signature names (index): {list(sig_raw.index)}")

# Transpose: cells x signatures
sig_weights = sig_raw.T
log(f"  Transposed shape: {sig_weights.shape} (cells x signatures)")
log(f"  First 3 cell barcodes: {list(sig_weights.index[:3])}")

# Check overlap with adata
sig_overlap = set(sig_weights.index) & set(adata.obs_names)
log(f"  Barcode overlap with adata_final: {len(sig_overlap)} / {len(sig_weights)}")

# Check overlap with basal cells
sig_basal_overlap = set(sig_weights.index) & set(basal_barcodes)
log(f"  Barcode overlap with basal cells: {len(sig_basal_overlap)}")

# Report key signature stats across ALL cells with weights
log("\n  Signature statistics (all cells with weights):")
for sig in sig_weights.columns:
    vals = sig_weights[sig]
    log(f"    {sig:8s}: mean={vals.mean():.4f}, median={vals.median():.4f}, "
        f"max={vals.max():.4f}, >0: {(vals>0).sum()}/{len(vals)}")

# =============================================================================
# STEP 4: LOAD SBS2 GROUP ASSIGNMENTS
# =============================================================================
log_sep("STEP 4: Load SBS2 group assignments")

groups = pd.read_csv(GROUP_ASSIGNMENTS_PATH, sep='\t', index_col=0)
log(f"  Loaded: {len(groups)} cells")
log(f"  Columns: {list(groups.columns)}")

# Detect group column name
group_col = 'group' if 'group' in groups.columns else groups.columns[0]
log(f"  Using group column: '{group_col}'")
for grp, count in groups[group_col].value_counts().items():
    log(f"    {grp}: {count}")

# =============================================================================
# STEP 5: BUILD BASAL CELL MASTER TABLE
# =============================================================================
log_sep("STEP 5: Build basal cell master table")

master = pd.DataFrame(index=basal_barcodes)
log(f"  Initialized: {len(master)} basal cells")

# --- Patient metadata ---
for col in ['subject id', 'tissue type', 'run_accession', 'Final_cancer_cell_status',
            'cnv_score', 'cnv_leiden', 'CytoTRACE2_Score', 'CytoTRACE2_Potency',
            'n_genes', 'n_counts']:
    if col in adata.obs.columns:
        master[col] = adata.obs.loc[basal_barcodes, col].values

# --- HPV16 reads ---
master['HPV16_reads'] = hpv16_series.reindex(basal_barcodes).fillna(0).values
master['HPV16_positive'] = (master['HPV16_reads'] > 0).astype(int)
n_hpv_basal = master['HPV16_positive'].sum()
log(f"  HPV16+ basal cells: {n_hpv_basal} / {n_basal} ({100*n_hpv_basal/n_basal:.1f}%)")

# --- A3 expression ---
log("\n  A3 family expression in basal cells:")
# Build index array once for efficiency
basal_idx = np.where(basal_mask)[0]
for gene in A3_GENES:
    if gene in adata.var_names:
        gene_pos = list(adata.var_names).index(gene)
        raw = adata.X[basal_idx, gene_pos]
        if hasattr(raw, 'toarray'):
            vals = raw.toarray().flatten()
        else:
            vals = np.array(raw).flatten()
        master[gene] = vals
        n_expr = (vals > 0).sum()
        log(f"    {gene:12s}: {n_expr}/{n_basal} ({100*n_expr/n_basal:.1f}%), mean={vals.mean():.3f}")
    else:
        log(f"    {gene:12s}: NOT FOUND")

# --- SBS2 group ---
master['SBS2_group'] = groups[group_col].reindex(master.index)
n_assigned = master['SBS2_group'].notna().sum()
log(f"\n  SBS2 group assigned: {n_assigned} / {n_basal}")

# --- Signature weights ---
for sig in sig_weights.columns:
    master[f'sig_{sig}'] = sig_weights[sig].reindex(master.index)

n_with_sigs = master['sig_SBS2'].notna().sum()
log(f"  Cells with signature weights: {n_with_sigs}")

# --- Final shape ---
log(f"\n  Master table: {master.shape[0]} rows x {master.shape[1]} columns")
log(f"  Columns: {list(master.columns)}")

# =============================================================================
# STEP 6: HPV16 DISTRIBUTION BY CELL TYPE (all cells, not just basal)
# =============================================================================
log_sep("STEP 6: HPV16 distribution across all cell types")

hpv16_all = hpv16_series.reindex(adata.obs_names).fillna(0)
adata.obs['HPV16_reads'] = hpv16_all.values
adata.obs['HPV16_positive'] = (adata.obs['HPV16_reads'] > 0).astype(int)

log(f"  {'Cell Type':40s} {'Total':>7s} {'HPV16+':>7s} {'%HPV16+':>8s} {'Mean reads':>11s}")
log(f"  {'-'*40} {'-'*7} {'-'*7} {'-'*8} {'-'*11}")
for ct in ct_counts.index:
    ct_mask = adata.obs['final_annotation'] == ct
    total = ct_mask.sum()
    pos = (adata.obs.loc[ct_mask, 'HPV16_positive'] == 1).sum()
    mean_reads = adata.obs.loc[ct_mask, 'HPV16_reads'].mean()
    log(f"  {ct:40s} {total:7d} {pos:7d} {100*pos/total:7.1f}% {mean_reads:10.3f}")

# =============================================================================
# STEP 7: HPV16 x PATIENT (basal cells only)
# =============================================================================
log_sep("STEP 7: HPV16 distribution by patient (basal cells)")

log(f"\n  {'Patient':20s} {'Total':>6s} {'HPV16+':>7s} {'HPV16-':>7s} {'%+':>6s} "
    f"{'Mean reads':>11s} {'Median(+)':>10s} {'Max':>6s}")
log(f"  {'-'*20} {'-'*6} {'-'*7} {'-'*7} {'-'*6} {'-'*11} {'-'*10} {'-'*6}")

patient_hpv_stats = []
for patient in sorted(master['subject id'].unique()):
    pmask = master['subject id'] == patient
    total = pmask.sum()
    pos = (master.loc[pmask, 'HPV16_positive'] == 1).sum()
    neg = total - pos
    mean_r = master.loc[pmask, 'HPV16_reads'].mean()
    
    pos_reads = master.loc[pmask & (master['HPV16_positive'] == 1), 'HPV16_reads']
    med_pos = pos_reads.median() if len(pos_reads) > 0 else 0
    max_r = pos_reads.max() if len(pos_reads) > 0 else 0
    
    log(f"  {patient:20s} {total:6d} {pos:7d} {neg:7d} {100*pos/total:5.1f}% "
        f"{mean_r:10.3f} {med_pos:10.2f} {max_r:5.1f}")
    
    patient_hpv_stats.append({
        'patient': patient, 'total_basal': total,
        'HPV16_pos': pos, 'HPV16_neg': neg,
        'pct_HPV16_pos': 100*pos/total, 'mean_HPV16_reads': mean_r
    })

patient_hpv_df = pd.DataFrame(patient_hpv_stats)
patient_hpv_df.to_csv(os.path.join(OUTPUT_DIR, "HPV16_distribution_by_patient.tsv"),
                       sep='\t', index=False)

# =============================================================================
# STEP 8: HPV16 x SBS2 GROUP
# =============================================================================
log_sep("STEP 8: HPV16 x SBS2 group cross-tabulation")

assigned = master.dropna(subset=['SBS2_group'])
log(f"  Cells with SBS2 assignment: {len(assigned)}")

ct_sbs2 = pd.crosstab(assigned['SBS2_group'], assigned['HPV16_positive'],
                       margins=True)
ct_sbs2.columns = ['HPV16-', 'HPV16+', 'Total']
ct_sbs2.index = [str(i) for i in ct_sbs2.index]
log(f"\n{ct_sbs2.to_string()}")

# Fisher's exact
table_2x2 = ct_sbs2.iloc[:2, :2].values
odds_ratio, fisher_p = fisher_exact(table_2x2)
log(f"\n  Fisher's exact: OR={odds_ratio:.3f}, p={fisher_p:.2e}")

ct_sbs2.to_csv(os.path.join(OUTPUT_DIR, "HPV16_x_SBS2_crosstab.tsv"), sep='\t')

# =============================================================================
# STEP 9: HPV16 x CANCER STATUS (basal cells)
# =============================================================================
log_sep("STEP 9: HPV16 x cancer cell status (basal cells)")

ct_cancer = pd.crosstab(master['Final_cancer_cell_status'], master['HPV16_positive'],
                         margins=True)
ct_cancer.columns = ['HPV16-', 'HPV16+', 'Total']
log(f"\n{ct_cancer.to_string()}")

table_2x2_cancer = ct_cancer.iloc[:2, :2].values
or_c, p_c = fisher_exact(table_2x2_cancer)
log(f"\n  Fisher's exact: OR={or_c:.3f}, p={p_c:.2e}")

ct_cancer.to_csv(os.path.join(OUTPUT_DIR, "HPV16_x_cancer_status_crosstab.tsv"), sep='\t')

# =============================================================================
# STEP 10: HPV16 x TISSUE TYPE (basal cells)
# =============================================================================
log_sep("STEP 10: HPV16 x tissue type (basal cells)")

ct_tissue = pd.crosstab(master['tissue type'], master['HPV16_positive'], margins=True)
ct_tissue.columns = ['HPV16-', 'HPV16+', 'Total']
log(f"\n{ct_tissue.to_string()}")

# =============================================================================
# STEP 11: HPV16 x CNV LEIDEN CLUSTERS (basal cells)
# =============================================================================
log_sep("STEP 11: HPV16 x cnv_leiden clusters (basal cells)")

ct_cnv = pd.crosstab(master['cnv_leiden'], master['HPV16_positive'])
ct_cnv.columns = ['HPV16-', 'HPV16+']
ct_cnv['Total'] = ct_cnv.sum(axis=1)
ct_cnv['%HPV16+'] = (100 * ct_cnv['HPV16+'] / ct_cnv['Total']).round(1)
ct_cnv = ct_cnv.sort_values('Total', ascending=False)
log(f"\n{ct_cnv.to_string()}")

ct_cnv.to_csv(os.path.join(OUTPUT_DIR, "HPV16_x_cnv_leiden_crosstab.tsv"), sep='\t')

# =============================================================================
# STEP 12: QUANTITATIVE COMPARISONS
# =============================================================================
log_sep("STEP 12: Quantitative comparisons across stratifications")

def compare_groups(df, group_col, group_a, group_b, metric_cols, label_a, label_b):
    """Compare two groups on multiple metrics with Mann-Whitney U."""
    log(f"\n  --- {label_a} vs {label_b} ---")
    mask_a = df[group_col] == group_a
    mask_b = df[group_col] == group_b
    n_a, n_b = mask_a.sum(), mask_b.sum()
    log(f"  n({label_a})={n_a}, n({label_b})={n_b}")
    
    log(f"  {'Metric':35s} {'Mean_A':>8s} {'Mean_B':>8s} {'U':>12s} {'p-value':>12s}")
    log(f"  {'-'*35} {'-'*8} {'-'*8} {'-'*12} {'-'*12}")
    
    for col in metric_cols:
        if col not in df.columns:
            continue
        a_vals = df.loc[mask_a, col].dropna()
        b_vals = df.loc[mask_b, col].dropna()
        if len(a_vals) < 5 or len(b_vals) < 5:
            continue
        u, p = mannwhitneyu(a_vals, b_vals, alternative='two-sided')
        log(f"  {col:35s} {a_vals.mean():8.4f} {b_vals.mean():8.4f} {u:12.0f} {p:12.2e}")

metrics = ['HPV16_reads', 'cnv_score', 'CytoTRACE2_Score',
           'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3H',
           'sig_SBS2', 'sig_SBS13', 'sig_SBS1', 'sig_SBS5', 'n_genes', 'n_counts']

# 12A: SBS2-HIGH vs SBS2-LOW (assigned cells only)
compare_groups(master.dropna(subset=['SBS2_group']),
               'SBS2_group', 'HIGH', 'LOW', metrics, 'SBS2-HIGH', 'SBS2-LOW')

# 12B: HPV16+ vs HPV16- (all basal cells)
compare_groups(master, 'HPV16_positive', 1, 0, 
               [m for m in metrics if m != 'HPV16_reads'], 'HPV16+', 'HPV16-')

# 12C: Cancer vs Normal (all basal cells)
compare_groups(master, 'Final_cancer_cell_status', 'Cancer cell', 'Normal cell',
               metrics, 'Cancer', 'Normal')

# =============================================================================
# STEP 13: 2x2 STRATIFICATION (SBS2 group x HPV16 status)
# =============================================================================
log_sep("STEP 13: SBS2 x HPV16 stratification (among assigned cells)")

assigned = master.dropna(subset=['SBS2_group']).copy()

log(f"\n  {'Subgroup':25s} {'n':>5s} {'A3A':>7s} {'A3B':>7s} {'CNV':>7s} "
    f"{'CytoTR':>7s} {'SBS2':>7s} {'SBS13':>7s}")
log(f"  {'-'*25} {'-'*5} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7}")

for grp in ['HIGH', 'LOW']:
    for hpv in [1, 0]:
        mask = (assigned['SBS2_group'] == grp) & (assigned['HPV16_positive'] == hpv)
        n = mask.sum()
        hpv_label = 'HPV+' if hpv == 1 else 'HPV-'
        label = f'{grp}/{hpv_label}'
        
        a3a = assigned.loc[mask, 'APOBEC3A'].mean() if n > 0 else 0
        a3b = assigned.loc[mask, 'APOBEC3B'].mean() if n > 0 else 0
        cnv = assigned.loc[mask, 'cnv_score'].mean() if n > 0 else 0
        cyto = assigned.loc[mask, 'CytoTRACE2_Score'].mean() if n > 0 else 0
        sbs2 = assigned.loc[mask, 'sig_SBS2'].mean() if 'sig_SBS2' in assigned.columns and n > 0 else np.nan
        sbs13 = assigned.loc[mask, 'sig_SBS13'].mean() if 'sig_SBS13' in assigned.columns and n > 0 else np.nan
        
        log(f"  {label:25s} {n:5d} {a3a:7.3f} {a3b:7.3f} {cnv:7.4f} "
            f"{cyto:7.4f} {sbs2:7.4f} {sbs13:7.4f}")

# =============================================================================
# STEP 14: BROADER BASAL LANDSCAPE (beyond the 1092 assigned cells)
# =============================================================================
log_sep("STEP 14: Broader basal cell landscape (all 52k cells)")

log(f"\n  Total basal cells: {n_basal}")
log(f"  In SBS2 HIGH/LOW: {n_assigned} ({100*n_assigned/n_basal:.1f}%)")
log(f"  Unassigned: {n_basal - n_assigned} ({100*(n_basal-n_assigned)/n_basal:.1f}%)")

# HPV16 by cancer status
log("\n  HPV16+ rate by cancer status (all basal cells):")
for status in master['Final_cancer_cell_status'].unique():
    smask = master['Final_cancer_cell_status'] == status
    total = smask.sum()
    pos = (master.loc[smask, 'HPV16_positive'] == 1).sum()
    log(f"    {status}: {pos}/{total} ({100*pos/total:.1f}%) HPV16+")

# HPV16 among unassigned
unassigned = master[master['SBS2_group'].isna()]
log(f"\n  Among unassigned basal cells (n={len(unassigned)}):")
log(f"    HPV16+: {unassigned['HPV16_positive'].sum()} ({100*unassigned['HPV16_positive'].mean():.1f}%)")
log(f"    Cancer: {(unassigned['Final_cancer_cell_status']=='Cancer cell').sum()}")
log(f"    Normal: {(unassigned['Final_cancer_cell_status']=='Normal cell').sum()}")
log(f"    Mean A3A: {unassigned['APOBEC3A'].mean():.3f}")
log(f"    Mean A3B: {unassigned['APOBEC3B'].mean():.3f}")
log(f"    Mean cnv_score: {unassigned['cnv_score'].mean():.4f}")
log(f"    Mean CytoTRACE2: {unassigned['CytoTRACE2_Score'].mean():.4f}")

# SBS2 weight distribution across all basal cells (from sig weights)
if 'sig_SBS2' in master.columns:
    sbs2_all = master['sig_SBS2'].dropna()
    log(f"\n  SBS2 weight distribution (cells with weights, n={len(sbs2_all)}):")
    log(f"    mean={sbs2_all.mean():.4f}, median={sbs2_all.median():.4f}")
    for pct in [25, 50, 75, 90, 95, 99]:
        log(f"    {pct}th percentile: {np.percentile(sbs2_all, pct):.4f}")

# =============================================================================
# STEP 15: CORRELATION ANALYSIS
# =============================================================================
log_sep("STEP 15: Correlations among key variables (all basal cells)")

corr_vars = ['HPV16_reads', 'APOBEC3A', 'APOBEC3B', 'APOBEC3C',
             'cnv_score', 'CytoTRACE2_Score']
corr_vars = [v for v in corr_vars if v in master.columns]

log(f"\n  Spearman correlations (all {n_basal} basal cells):")
log(f"  {'Var1':20s} {'Var2':20s} {'rho':>8s} {'p-value':>12s}")
log(f"  {'-'*20} {'-'*20} {'-'*8} {'-'*12}")

for i in range(len(corr_vars)):
    for j in range(i+1, len(corr_vars)):
        v1 = master[corr_vars[i]].dropna()
        v2 = master[corr_vars[j]].dropna()
        common = v1.index.intersection(v2.index)
        if len(common) > 10:
            rho, p = spearmanr(v1.loc[common], v2.loc[common])
            log(f"  {corr_vars[i]:20s} {corr_vars[j]:20s} {rho:8.4f} {p:12.2e}")

# Also among HPV16+ cells only
hpv_pos_mask = master['HPV16_positive'] == 1
if hpv_pos_mask.sum() > 100:
    log(f"\n  Spearman correlations (HPV16+ basal cells only, n={hpv_pos_mask.sum()}):")
    log(f"  {'Var1':20s} {'Var2':20s} {'rho':>8s} {'p-value':>12s}")
    log(f"  {'-'*20} {'-'*20} {'-'*8} {'-'*12}")
    for i in range(len(corr_vars)):
        for j in range(i+1, len(corr_vars)):
            v1 = master.loc[hpv_pos_mask, corr_vars[i]].dropna()
            v2 = master.loc[hpv_pos_mask, corr_vars[j]].dropna()
            common = v1.index.intersection(v2.index)
            if len(common) > 10:
                rho, p = spearmanr(v1.loc[common], v2.loc[common])
                log(f"  {corr_vars[i]:20s} {corr_vars[j]:20s} {rho:8.4f} {p:12.2e}")

# =============================================================================
# STEP 16: SAVE OUTPUTS
# =============================================================================
log_sep("STEP 16: Save outputs")

# Save master table
master_path = os.path.join(OUTPUT_DIR, "basal_cell_master_table.tsv")
master.to_csv(master_path, sep='\t')
log(f"  Saved: {master_path}")
log(f"    Shape: {master.shape}")

# Save diagnostic report
report_path = os.path.join(OUTPUT_DIR, "diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")

# =============================================================================
# STEP 17: PHASE 1 SUMMARY AND NEXT STEPS
# =============================================================================
log_sep("PHASE 1 SUMMARY")

log(f"""
  DATASET OVERVIEW:
    Total cells:           {adata.n_obs}
    Basal cells:           {n_basal} ({100*n_basal/adata.n_obs:.1f}%)
    HPV16+ basal cells:    {n_hpv_basal} ({100*n_hpv_basal/n_basal:.1f}%)
    In SBS2 HIGH/LOW:      {n_assigned}
    With signature weights: {n_with_sigs}
  
  KEY QUESTIONS FOR PHASE 2:
    1. What threshold separates HPV16-high from HPV16-low?
       (Elbow detection on read distribution among HPV16+ cells)
    2. Do cnv_leiden clusters separate HPV16+/- and SBS2-HIGH/LOW?
    3. What is the A3A vs A3B vs CNV landscape across all basal cells?
    4. Can we define Daiko's two populations from the data?
       (Pop1: high-SBS2/A3A/HPV-active vs Pop2: high-CNV/A3B/stealth)
""")

log("=" * 80)
log("  PHASE 1 v2 COMPLETE")
log("=" * 80)
