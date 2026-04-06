#!/usr/bin/env python3
"""
Phase3_HPV16_Populations_and_Genome.py
========================================
Figure 6/7 Phase 3:
  Step 1: Elbow detection for HPV16 raw count threshold
  Step 2: Define basal cell populations (HPV-neg / HPV-low / HPV-high x cnv_leiden)
  Step 3: A3A vs A3B divergence across populations
  Step 4: Library size correction (partial correlations)
  Step 5: Probe for Kraken2 intermediate files and pilot HPV16 genome alignment

Inputs:
    - data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv
    - Kraken2 intermediate files (unmapped BAMs, kraken output)

Outputs (to data/FIG_6/02_populations/):
    - elbow_detection_report.tsv
    - elbow_detection_plot.pdf
    - population_assignments.tsv
    - population_summary.tsv
    - divergence_tests.tsv
    - partial_correlations.tsv
    
Outputs (to data/FIG_6/03_hpv16_genome/):
    - kraken2_intermediate_probe.txt
    - pilot_hpv16_gene_counts.tsv  (if intermediate files available)

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, kruskal
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
KRAKEN2_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"

MASTER_TABLE_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")

OUTPUT_POP_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/02_populations")
OUTPUT_GENOME_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/03_hpv16_genome")
os.makedirs(OUTPUT_POP_DIR, exist_ok=True)
os.makedirs(OUTPUT_GENOME_DIR, exist_ok=True)

# HPV16 genome annotation (NC_001526.4, 7906 bp)
HPV16_GENES = {
    'URR_5prime': (1, 82),       # Upstream Regulatory Region (wraps around)
    'E6':         (83, 559),
    'E7':         (562, 858),
    'E1':         (865, 2813),
    'E2':         (2755, 3852),
    'E4':         (3332, 3619),   # overlaps E2
    'E5':         (3849, 4100),
    'L2':         (4236, 5657),
    'L1':         (5560, 7155),
    'URR_3prime': (7156, 7906),   # Upstream Regulatory Region
}
HPV16_GENOME_LENGTH = 7906
HPV16_TAXID = '333760'

# =============================================================================
# LOGGING
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
# STEP 0: LOAD DATA
# =============================================================================
log_sep("STEP 0: Load master table")

master = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
log(f"  Master table: {master.shape}")
log(f"  Columns: {list(master.columns)}")

n_basal = len(master)
log(f"  Total basal cells: {n_basal}")
log(f"  HPV16+ (raw > 0): {(master['raw_HPV16'] > 0).sum()}")

# =============================================================================
# STEP 1: ELBOW DETECTION FOR HPV16 THRESHOLD
# =============================================================================
log_sep("STEP 1: Elbow detection for HPV16 raw count threshold")

# Get distribution of raw counts among HPV16+ cells
pos_counts = master.loc[master['raw_HPV16'] > 0, 'raw_HPV16'].values
pos_counts_sorted = np.sort(pos_counts)

log(f"  HPV16+ cells: {len(pos_counts)}")
log(f"  Count range: {pos_counts.min()} to {pos_counts.max()}")

# --- Method 1: Two-line piecewise linear fit (L-method) ---
# Same approach used for SBS2 threshold in Figure 4
# Fit two linear segments to the sorted count distribution and find optimal breakpoint

# Work with log-scale for the long-tailed distribution
log_counts = np.log1p(pos_counts_sorted)
x_vals = np.arange(len(log_counts))

def two_line_residual(break_idx, x, y):
    """Compute total residual for two-line fit at breakpoint."""
    break_idx = int(break_idx)
    if break_idx < 5 or break_idx >= len(x) - 5:
        return 1e10
    
    # Left segment
    x1, y1 = x[:break_idx], y[:break_idx]
    if len(x1) > 1:
        p1 = np.polyfit(x1, y1, 1)
        res1 = np.sum((y1 - np.polyval(p1, x1))**2)
    else:
        res1 = 0
    
    # Right segment
    x2, y2 = x[break_idx:], y[break_idx:]
    if len(x2) > 1:
        p2 = np.polyfit(x2, y2, 1)
        res2 = np.sum((y2 - np.polyval(p2, x2))**2)
    else:
        res2 = 0
    
    return res1 + res2

# Scan breakpoints
n_points = len(log_counts)
scan_range = range(max(10, int(0.01 * n_points)), min(int(0.5 * n_points), n_points - 10), max(1, int(0.001 * n_points)))
residuals = []
for bp in scan_range:
    res = two_line_residual(bp, x_vals, log_counts)
    residuals.append((bp, res))

residuals_df = pd.DataFrame(residuals, columns=['breakpoint_idx', 'residual'])
best_bp = residuals_df.loc[residuals_df['residual'].idxmin(), 'breakpoint_idx']
best_bp = int(best_bp)
threshold_count = pos_counts_sorted[best_bp]

log(f"\n  L-method elbow detection:")
log(f"    Optimal breakpoint index: {best_bp} / {n_points}")
log(f"    Raw count at breakpoint: {threshold_count}")
log(f"    Cells below threshold: {best_bp} ({100*best_bp/n_points:.1f}%)")
log(f"    Cells at or above threshold: {n_points - best_bp} ({100*(n_points-best_bp)/n_points:.1f}%)")

# --- Method 2: Simple fixed thresholds for comparison ---
log(f"\n  Fixed threshold comparison:")
for thresh in [1, 2, 3, 5, 10]:
    n_above = (pos_counts >= thresh).sum()
    n_below = (pos_counts < thresh).sum()
    log(f"    >= {thresh} UMI: {n_above} cells ({100*n_above/len(pos_counts):.1f}% of HPV16+), "
        f"below: {n_below}")

# --- Method 3: Examine count=1 cells more carefully ---
# Are count=1 cells distributed across all patients, or concentrated?
log(f"\n  Cells with exactly 1 UMI read by patient:")
one_read = master[master['raw_HPV16'] == 1]
multi_read = master[master['raw_HPV16'] > 1]
log(f"  {'Patient':20s} {'1-read':>7s} {'>1-read':>8s} {'Ratio':>7s}")
log(f"  {'-'*20} {'-'*7} {'-'*8} {'-'*7}")
for patient in sorted(master['subject id'].unique()):
    n_one = (one_read['subject id'] == patient).sum()
    n_multi = (multi_read['subject id'] == patient).sum()
    ratio = n_one / n_multi if n_multi > 0 else float('inf')
    log(f"  {patient:20s} {n_one:7d} {n_multi:8d} {ratio:7.2f}")

# SC008 (HPV-negative patient) is the key control:
# if SC008 has count=1 cells, those are almost certainly ambient
n_sc008_one = ((master['subject id'] == 'Patient SC008') & (master['raw_HPV16'] == 1)).sum()
n_sc008_any = ((master['subject id'] == 'Patient SC008') & (master['raw_HPV16'] > 0)).sum()
log(f"\n  SC008 (HPV-negative control): {n_sc008_any} cells with any HPV16 reads, {n_sc008_one} with exactly 1")
log(f"  --> If SC008 has count=1 cells, threshold should be >= 2")

# --- Decision ---
# Use the more conservative of L-method or the ambient-informed threshold
if n_sc008_any > 0:
    ambient_thresh = master.loc[master['subject id'] == 'Patient SC008', 'raw_HPV16'].max() + 1
    log(f"\n  Ambient-informed threshold (SC008 max + 1): {ambient_thresh}")
else:
    ambient_thresh = 2  # Default conservative

THRESHOLD = max(threshold_count, ambient_thresh, 2)  # At least 2
log(f"\n  *** SELECTED THRESHOLD: >= {THRESHOLD} UMI reads ***")

# Save elbow data
residuals_df.to_csv(os.path.join(OUTPUT_POP_DIR, "elbow_detection_residuals.tsv"),
                     sep='\t', index=False)

# =============================================================================
# STEP 2: DEFINE POPULATIONS
# =============================================================================
log_sep("STEP 2: Define basal cell populations")

# Primary stratification: HPV16 status based on threshold
master['HPV16_status'] = 'HPV16-negative'
master.loc[master['raw_HPV16'] >= THRESHOLD, 'HPV16_status'] = 'HPV16-positive'
# Cells with 1 to THRESHOLD-1 reads: ambiguous
master.loc[(master['raw_HPV16'] > 0) & (master['raw_HPV16'] < THRESHOLD), 'HPV16_status'] = 'HPV16-ambiguous'

for status, count in master['HPV16_status'].value_counts().items():
    log(f"  {status}: {count} ({100*count/n_basal:.1f}%)")

# Secondary: HPV16-positive split into dose levels (tertiles among positive cells)
hpv_pos_mask = master['HPV16_status'] == 'HPV16-positive'
hpv_pos_vals = master.loc[hpv_pos_mask, 'raw_HPV16']

if len(hpv_pos_vals) > 30:
    tertile_33 = np.percentile(hpv_pos_vals, 33.3)
    tertile_67 = np.percentile(hpv_pos_vals, 66.7)
    
    master['HPV16_dose'] = 'none'
    master.loc[master['HPV16_status'] == 'HPV16-ambiguous', 'HPV16_dose'] = 'ambiguous'
    master.loc[hpv_pos_mask & (master['raw_HPV16'] <= tertile_33), 'HPV16_dose'] = 'low'
    master.loc[hpv_pos_mask & (master['raw_HPV16'] > tertile_33) & (master['raw_HPV16'] <= tertile_67), 'HPV16_dose'] = 'medium'
    master.loc[hpv_pos_mask & (master['raw_HPV16'] > tertile_67), 'HPV16_dose'] = 'high'
    
    log(f"\n  HPV16 dose levels (among HPV16-positive, tertile boundaries: {tertile_33:.0f}, {tertile_67:.0f}):")
    for dose in ['none', 'ambiguous', 'low', 'medium', 'high']:
        n = (master['HPV16_dose'] == dose).sum()
        if n > 0:
            mean_raw = master.loc[master['HPV16_dose'] == dose, 'raw_HPV16'].mean()
            log(f"    {dose:12s}: {n:6d} cells, mean raw HPV16 = {mean_raw:.1f}")

# =============================================================================
# Population matrix: HPV16_status x cnv_leiden
# =============================================================================
log(f"\n  HPV16 status x cnv_leiden clusters:")
ct_pop = pd.crosstab(master['cnv_leiden'], master['HPV16_status'])
ct_pop['Total'] = ct_pop.sum(axis=1)
ct_pop = ct_pop.sort_values('Total', ascending=False)

# Calculate % positive for each cluster
if 'HPV16-positive' in ct_pop.columns:
    ct_pop['%Positive'] = (100 * ct_pop['HPV16-positive'] / ct_pop['Total']).round(1)
log(ct_pop.to_string())

# =============================================================================
# STEP 3: A3A vs A3B DIVERGENCE ACROSS POPULATIONS
# =============================================================================
log_sep("STEP 3: A3A vs A3B divergence across populations")

pop_stats = []
metrics = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'cnv_score', 'CytoTRACE2_Score',
           'raw_HPV16', 'sig_SBS2', 'sig_SBS13', 'n_genes']
metrics = [m for m in metrics if m in master.columns]

# By HPV16 status
log(f"\n  --- By HPV16 status ---")
log(f"  {'Status':18s} {'n':>6s} {'A3A':>7s} {'A3B':>7s} {'CNV':>7s} "
    f"{'CytoTR':>7s} {'SBS2':>7s} {'HPV16':>8s}")
log(f"  {'-'*18} {'-'*6} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*8}")

for status in ['HPV16-positive', 'HPV16-ambiguous', 'HPV16-negative']:
    mask = master['HPV16_status'] == status
    n = mask.sum()
    if n == 0:
        continue
    
    row = {'population': status, 'n_cells': n}
    a3a = master.loc[mask, 'APOBEC3A'].mean()
    a3b = master.loc[mask, 'APOBEC3B'].mean()
    cnv = master.loc[mask, 'cnv_score'].mean()
    cyto = master.loc[mask, 'CytoTRACE2_Score'].mean()
    sbs2 = master.loc[mask, 'sig_SBS2'].mean() if 'sig_SBS2' in master.columns else np.nan
    hpv = master.loc[mask, 'raw_HPV16'].mean()
    
    for m in metrics:
        row[f'{m}_mean'] = master.loc[mask, m].mean()
    pop_stats.append(row)
    
    log(f"  {status:18s} {n:6d} {a3a:7.3f} {a3b:7.3f} {cnv:7.4f} "
        f"{cyto:7.4f} {sbs2:7.4f} {hpv:8.1f}")

# By HPV16 dose (finer resolution)
log(f"\n  --- By HPV16 dose ---")
log(f"  {'Dose':12s} {'n':>6s} {'A3A':>7s} {'A3B':>7s} {'CNV':>7s} "
    f"{'CytoTR':>7s} {'SBS2':>7s} {'HPV16':>8s}")
log(f"  {'-'*12} {'-'*6} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*8}")

for dose in ['none', 'ambiguous', 'low', 'medium', 'high']:
    mask = master['HPV16_dose'] == dose
    n = mask.sum()
    if n == 0:
        continue
    
    a3a = master.loc[mask, 'APOBEC3A'].mean()
    a3b = master.loc[mask, 'APOBEC3B'].mean()
    cnv = master.loc[mask, 'cnv_score'].mean()
    cyto = master.loc[mask, 'CytoTRACE2_Score'].mean()
    sbs2 = master.loc[mask, 'sig_SBS2'].mean() if 'sig_SBS2' in master.columns else np.nan
    hpv = master.loc[mask, 'raw_HPV16'].mean()
    
    log(f"  {dose:12s} {n:6d} {a3a:7.3f} {a3b:7.3f} {cnv:7.4f} "
        f"{cyto:7.4f} {sbs2:7.4f} {hpv:8.1f}")

# Pairwise statistical tests
log(f"\n  --- Kruskal-Wallis across HPV16 dose groups (none/low/med/high) ---")
for metric in metrics:
    groups = []
    for dose in ['none', 'low', 'medium', 'high']:
        vals = master.loc[master['HPV16_dose'] == dose, metric].dropna()
        if len(vals) > 5:
            groups.append(vals)
    if len(groups) >= 3:
        h_stat, kw_p = kruskal(*groups)
        log(f"  {metric:25s}: H={h_stat:.1f}, p={kw_p:.2e}")

# Key pairwise: HPV16-high vs HPV16-negative
log(f"\n  --- Mann-Whitney: HPV16-positive vs HPV16-negative ---")
for metric in metrics:
    if metric == 'raw_HPV16':
        continue
    pos_vals = master.loc[master['HPV16_status'] == 'HPV16-positive', metric].dropna()
    neg_vals = master.loc[master['HPV16_status'] == 'HPV16-negative', metric].dropna()
    if len(pos_vals) > 5 and len(neg_vals) > 5:
        u, p = mannwhitneyu(pos_vals, neg_vals, alternative='two-sided')
        log(f"  {metric:25s}: pos_mean={pos_vals.mean():.4f}, neg_mean={neg_vals.mean():.4f}, p={p:.2e}")

# =============================================================================
# STEP 3B: A3A/A3B RATIO ANALYSIS
# =============================================================================
log_sep("STEP 3B: A3A/A3B expression ratio across populations")

# Calculate A3A / (A3A + A3B) ratio (avoids division by zero)
master['A3A_fraction'] = master['APOBEC3A'] / (master['APOBEC3A'] + master['APOBEC3B'] + 0.01)

log(f"  A3A / (A3A + A3B + 0.01) ratio by HPV16 dose:")
log(f"  {'Dose':12s} {'n':>6s} {'A3A_frac':>9s} {'A3A':>7s} {'A3B':>7s}")
log(f"  {'-'*12} {'-'*6} {'-'*9} {'-'*7} {'-'*7}")
for dose in ['none', 'ambiguous', 'low', 'medium', 'high']:
    mask = master['HPV16_dose'] == dose
    n = mask.sum()
    if n == 0:
        continue
    frac = master.loc[mask, 'A3A_fraction'].mean()
    a3a = master.loc[mask, 'APOBEC3A'].mean()
    a3b = master.loc[mask, 'APOBEC3B'].mean()
    log(f"  {dose:12s} {n:6d} {frac:9.4f} {a3a:7.3f} {a3b:7.3f}")

# By cnv_leiden (top 10 clusters)
log(f"\n  A3A/(A3A+A3B) ratio by cnv_leiden cluster (top 10 by size):")
cnv_sizes = master['cnv_leiden'].value_counts().head(10)
log(f"  {'Cluster':10s} {'n':>6s} {'A3A_frac':>9s} {'A3A':>7s} {'A3B':>7s} "
    f"{'CNV':>7s} {'%HPV16+':>8s}")
log(f"  {'-'*10} {'-'*6} {'-'*9} {'-'*7} {'-'*7} {'-'*7} {'-'*8}")
for cluster in cnv_sizes.index:
    mask = master['cnv_leiden'] == cluster
    n = mask.sum()
    frac = master.loc[mask, 'A3A_fraction'].mean()
    a3a = master.loc[mask, 'APOBEC3A'].mean()
    a3b = master.loc[mask, 'APOBEC3B'].mean()
    cnv = master.loc[mask, 'cnv_score'].mean()
    pct_hpv = 100 * (master.loc[mask, 'HPV16_status'] == 'HPV16-positive').sum() / n
    log(f"  {str(cluster):10s} {n:6d} {frac:9.4f} {a3a:7.3f} {a3b:7.3f} "
        f"{cnv:7.4f} {pct_hpv:7.1f}%")

# =============================================================================
# STEP 4: LIBRARY SIZE CORRECTION (PARTIAL CORRELATIONS)
# =============================================================================
log_sep("STEP 4: Partial correlations controlling for library size")

hpv_pos = master[master['HPV16_status'] == 'HPV16-positive'].copy()
log(f"  HPV16-positive cells: {len(hpv_pos)}")

# Raw HPV16 / n_genes ratio (reads per 1000 genes detected)
hpv_pos['HPV16_per_1k_genes'] = 1000 * hpv_pos['raw_HPV16'] / hpv_pos['n_genes'].clip(lower=1)

log(f"\n  HPV16 per 1000 genes distribution:")
rate_vals = hpv_pos['HPV16_per_1k_genes']
log(f"    Mean:   {rate_vals.mean():.3f}")
log(f"    Median: {rate_vals.median():.3f}")
log(f"    Std:    {rate_vals.std():.3f}")

# Partial correlations: raw_HPV16 ~ target | n_genes
log(f"\n  Partial Spearman correlations (raw_HPV16 ~ X | n_genes):")
log(f"  {'Variable':25s} {'raw rho':>8s} {'partial rho':>12s} {'partial p':>12s}")
log(f"  {'-'*25} {'-'*8} {'-'*12} {'-'*12}")

targets = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'cnv_score', 'CytoTRACE2_Score']
targets = [t for t in targets if t in hpv_pos.columns]

for target in targets:
    # Drop NaN
    subset = hpv_pos[['raw_HPV16', target, 'n_genes']].dropna()
    if len(subset) < 20:
        continue
    
    # Raw correlation
    raw_rho, _ = spearmanr(subset['raw_HPV16'], subset[target])
    
    # Partial correlation via residualization
    # Rank-transform, then regress out n_genes from both variables
    from scipy.stats import rankdata
    r_hpv = rankdata(subset['raw_HPV16'])
    r_target = rankdata(subset[target])
    r_ngenes = rankdata(subset['n_genes'])
    
    # Residualize HPV16 | n_genes
    p_hpv = np.polyfit(r_ngenes, r_hpv, 1)
    res_hpv = r_hpv - np.polyval(p_hpv, r_ngenes)
    
    # Residualize target | n_genes
    p_target = np.polyfit(r_ngenes, r_target, 1)
    res_target = r_target - np.polyval(p_target, r_ngenes)
    
    # Partial correlation
    partial_rho, partial_p = spearmanr(res_hpv, res_target)
    
    log(f"  {target:25s} {raw_rho:8.4f} {partial_rho:12.4f} {partial_p:12.2e}")

# Also: rate-based correlations (HPV16_per_1k_genes)
log(f"\n  Rate-based Spearman correlations (HPV16_per_1k_genes ~ X):")
log(f"  {'Variable':25s} {'rho':>8s} {'p-value':>12s}")
log(f"  {'-'*25} {'-'*8} {'-'*12}")
for target in targets:
    subset = hpv_pos[['HPV16_per_1k_genes', target]].dropna()
    if len(subset) > 20:
        rho, p = spearmanr(subset['HPV16_per_1k_genes'], subset[target])
        log(f"  {target:25s} {rho:8.4f} {p:12.2e}")

# =============================================================================
# STEP 5: HPV16 GENOME MAPPING — PROBE AND PILOT
# =============================================================================
log_sep("STEP 5: HPV16 genome mapping — probe for intermediate files")

# --- 5A: Check what Kraken2 intermediate files exist ---
log("  Checking Kraken2 intermediate files for first 3 samples...")

sample_dirs = sorted(glob.glob(os.path.join(KRAKEN2_BASE, "SRR*", "SRR*_S1_L001_", "outs")))
probe_results = []

for sample_dir in sample_dirs[:5]:  # check first 5
    srr_id = [p for p in sample_dir.split('/') if p.startswith('SRR')][0]
    
    files_to_check = {
        'unmapped_bam': os.path.join(sample_dir, f"{srr_id}_unmapped.bam"),
        'unmapped_fq': os.path.join(sample_dir, f"{srr_id}_unmapped.fq"),
        'kraken_output': os.path.join(sample_dir, f"{srr_id}_kraken.out"),
        'kraken_report': os.path.join(sample_dir, f"{srr_id}_kraken_report.txt"),
    }
    
    # Also check alternative locations/naming
    alt_files = {
        'unmapped_bam_alt': os.path.join(os.path.dirname(sample_dir), f"{srr_id}_unmapped.bam"),
        'kraken_out_alt': os.path.join(os.path.dirname(sample_dir), f"{srr_id}_kraken.out"),
    }
    files_to_check.update(alt_files)
    
    found = {}
    for name, path in files_to_check.items():
        if os.path.exists(path):
            size = os.path.getsize(path)
            found[name] = {'path': path, 'size_MB': size / 1e6}
    
    probe_results.append({'srr_id': srr_id, 'files_found': found})
    
    log(f"\n  [{srr_id}]:")
    if found:
        for name, info in found.items():
            log(f"    [FOUND] {name}: {info['path']} ({info['size_MB']:.1f} MB)")
    else:
        log(f"    No intermediate files found in {sample_dir}")

# Also do a broader search
log(f"\n  Broader search for Kraken2 intermediates...")
for pattern in ['*_unmapped.bam', '*_kraken.out', '*_unmapped.fq', '*_unmapped.fq.gz']:
    found_files = glob.glob(os.path.join(KRAKEN2_BASE, "SRR*", "**", pattern), recursive=True)
    if found_files:
        log(f"  {pattern}: found {len(found_files)} files")
        log(f"    Example: {found_files[0]}")
    else:
        log(f"  {pattern}: NOT FOUND")

# Also check parent directories
for pattern in ['*unmapped*', '*kraken*']:
    found_files = glob.glob(os.path.join(KRAKEN2_BASE, "..", pattern))
    if found_files:
        log(f"  Parent dir {pattern}: {found_files[:3]}")

# --- 5B: Check for Cell Ranger BAM files (fallback) ---
log(f"\n  Checking Cell Ranger BAM files (fallback for re-extraction)...")
bam_pattern = os.path.join(KRAKEN2_BASE, "SRR*", "SRR*_S1_L001_", "outs", "possorted_genome_bam.bam")
bam_files = glob.glob(bam_pattern)
log(f"  Cell Ranger BAMs found: {len(bam_files)}")
if bam_files:
    log(f"    Example: {bam_files[0]}")
    bam_size = os.path.getsize(bam_files[0]) / 1e9
    log(f"    Size: {bam_size:.1f} GB")

# --- 5C: Report HPV16 genome annotation for reference ---
log(f"\n  HPV16 genome annotation (NC_001526.4, {HPV16_GENOME_LENGTH} bp):")
log(f"  {'Gene':12s} {'Start':>6s} {'End':>6s} {'Length':>7s} {'Category':>12s}")
log(f"  {'-'*12} {'-'*6} {'-'*6} {'-'*7} {'-'*12}")

gene_categories = {
    'E6': 'early/oncogene', 'E7': 'early/oncogene',
    'E1': 'early/replication', 'E2': 'early/regulation',
    'E4': 'early/assembly', 'E5': 'early/signaling',
    'L1': 'late/capsid', 'L2': 'late/capsid',
    'URR_5prime': 'regulatory', 'URR_3prime': 'regulatory'
}
for gene, (start, end) in HPV16_GENES.items():
    length = end - start + 1
    cat = gene_categories.get(gene, 'unknown')
    log(f"  {gene:12s} {start:6d} {end:6d} {length:7d} {cat:>12s}")

# --- 5D: Determine approach and create pilot script ---
has_unmapped_bam = any(
    any('unmapped_bam' in k for k in r['files_found']) 
    for r in probe_results
)
has_kraken_output = any(
    any('kraken' in k and 'out' in k for k in r['files_found'])
    for r in probe_results
)
has_cellranger_bam = len(bam_files) > 0

log(f"\n  --- Approach determination ---")
log(f"  Kraken2 unmapped BAMs available: {has_unmapped_bam}")
log(f"  Kraken2 output files available: {has_kraken_output}")
log(f"  Cell Ranger BAMs available: {has_cellranger_bam}")

if has_unmapped_bam and has_kraken_output:
    approach = "KRAKEN_INTERMEDIATE"
    log(f"\n  RECOMMENDED APPROACH: Use Kraken2 intermediates")
    log(f"  Steps:")
    log(f"    1. Read kraken.out, filter for taxid {HPV16_TAXID} (HPV16)")
    log(f"    2. Extract matching reads from unmapped BAM (has CB/UB tags)")
    log(f"    3. Align extracted reads to HPV16 reference (minimap2/bowtie2)")
    log(f"    4. Count per-cell per-gene alignments")
elif has_cellranger_bam:
    approach = "CELLRANGER_BAM"
    log(f"\n  RECOMMENDED APPROACH: Re-extract from Cell Ranger BAMs")
    log(f"  Steps:")
    log(f"    1. Extract unmapped reads from possorted_genome_bam.bam")
    log(f"    2. Align unmapped reads to HPV16 reference genome")
    log(f"    3. Filter for aligned reads, extract CB tags")
    log(f"    4. Count per-cell per-gene alignments")
    log(f"  NOTE: This approach is slower but does not depend on intermediates")
else:
    approach = "NONE"
    log(f"\n  WARNING: No viable approach found")
    log(f"  Need either Kraken2 intermediates or Cell Ranger BAMs")

# =============================================================================
# STEP 6: SAVE OUTPUTS
# =============================================================================
log_sep("STEP 6: Save outputs")

# Population assignments
pop_cols = ['HPV16_status', 'HPV16_dose', 'raw_HPV16', 'HPV16_per_1k_genes'
            ] if 'HPV16_per_1k_genes' in master.columns else ['HPV16_status', 'HPV16_dose', 'raw_HPV16']
# HPV16_per_1k_genes only exists for HPV16+ cells, add to master
if 'HPV16_per_1k_genes' not in master.columns:
    master['HPV16_per_1k_genes'] = 1000 * master['raw_HPV16'] / master['n_genes'].clip(lower=1)

master['A3A_fraction'] = master['APOBEC3A'] / (master['APOBEC3A'] + master['APOBEC3B'] + 0.01)

pop_path = os.path.join(OUTPUT_POP_DIR, "population_assignments.tsv")
master.to_csv(pop_path, sep='\t')
log(f"  Saved: {pop_path}")
log(f"    Shape: {master.shape}")

# Population summary
if pop_stats:
    pop_stats_df = pd.DataFrame(pop_stats)
    pop_stats_df.to_csv(os.path.join(OUTPUT_POP_DIR, "population_summary.tsv"),
                         sep='\t', index=False)
    log(f"  Saved: population_summary.tsv")

# Report
report_path = os.path.join(OUTPUT_POP_DIR, "phase3_diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")

# Genome mapping approach file
approach_path = os.path.join(OUTPUT_GENOME_DIR, "approach_determination.txt")
with open(approach_path, 'w') as f:
    f.write(f"approach={approach}\n")
    f.write(f"has_unmapped_bam={has_unmapped_bam}\n")
    f.write(f"has_kraken_output={has_kraken_output}\n")
    f.write(f"has_cellranger_bam={has_cellranger_bam}\n")
    f.write(f"n_cellranger_bams={len(bam_files)}\n")
    if bam_files:
        f.write(f"example_bam={bam_files[0]}\n")
log(f"  Saved: {approach_path}")

log_sep("PHASE 3 COMPLETE")
log(f"""
  SUMMARY:
    HPV16 threshold: >= {THRESHOLD} UMI reads
    HPV16-positive basal cells: {(master['HPV16_status']=='HPV16-positive').sum()}
    HPV16-negative basal cells: {(master['HPV16_status']=='HPV16-negative').sum()}
    HPV16-ambiguous basal cells: {(master['HPV16_status']=='HPV16-ambiguous').sum()}
    HPV16 genome mapping approach: {approach}
    
  NEXT STEPS:
    1. Generate HPV16 genome alignment script (based on approach: {approach})
    2. Visualize populations on UMAP
    3. Differential expression between populations
    4. Map viral reads to HPV16 genes per population
""")
