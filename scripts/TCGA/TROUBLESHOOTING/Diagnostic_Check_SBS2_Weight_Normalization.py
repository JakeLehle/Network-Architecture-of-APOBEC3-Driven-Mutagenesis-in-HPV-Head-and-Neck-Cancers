#!/usr/bin/env python3
"""
Diagnostic_Check_SBS2_Weight_Normalization.py
================================================
Quick diagnostic to determine whether TCGA SBS mutational signature weights
are normalized to a constant total per sample or represent true absolute
context-specific mutation counts.

Daiko's concern (4/10 meeting): If the signature deconvolution was performed
on a matrix where each sample's 96-context mutation counts were normalized
to the same total (e.g., sum to 1.0, or sum to some fixed constant), then
a sample with 8 total mutations and 4 SBS2-context mutations would get the
same SBS2 weight as a sample with 800 mutations and 400 SBS2-context mutations.
This would make low-mutation-burden samples appear equivalent to high-mutation
samples in SBS2 weight, biasing the analysis.

The values are known to be integer-based, not fractional. So the question is:
do all samples sum to the SAME total (constant normalization), or do the sums
vary across samples (reflecting true absolute mutation counts)?

What this script checks:
  1. Sum all SBS columns per sample
  2. Check whether all sums are identical (constant) or vary
  3. If constant → weights are normalized, need VCF recomputation
  4. If varying → weights are absolute counts, current approach is valid

Input:
  - data/FIG_1/Mutation_Table_Tumors_TCGA.tsv

Output:
  - Console report with verdict
  - data/FIG_1/TROUBLESHOOTING/sbs_weight_normalization_check.tsv
  - data/FIG_1/TROUBLESHOOTING/sbs_weight_normalization_report.txt

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MUTATION_SIGNATURE_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)

# =============================================================================
# MAIN
# =============================================================================
banner("DIAGNOSTIC: SBS Weight Normalization Check")

log(f"  Input: {MUTATION_SIGNATURE_PATH}")

# Load the mutation signature file
log(f"\n  Loading mutation signature table...")
df = pd.read_csv(MUTATION_SIGNATURE_PATH, sep='\t')
log(f"  Shape: {df.shape}")
log(f"  Columns (first 10): {list(df.columns[:10])}")
log(f"  Total columns: {len(df.columns)}")

# Identify SBS columns
sbs_cols = [c for c in df.columns if c.startswith('SBS')]
log(f"\n  SBS signature columns found: {len(sbs_cols)}")
log(f"  SBS columns: {sbs_cols}")

if len(sbs_cols) == 0:
    log("\n  ERROR: No SBS columns found. Check column naming.")
    sys.exit(1)

# Check data types — confirm integer-based
log(f"\n  Data type check (first SBS column: {sbs_cols[0]}):")
log(f"    dtype: {df[sbs_cols[0]].dtype}")
sample_vals = df[sbs_cols[0]].head(10).tolist()
log(f"    First 10 values: {sample_vals}")
all_integer = all(df[c].dropna().apply(lambda x: float(x).is_integer()).all() for c in sbs_cols)
log(f"    All values are integers: {all_integer}")

# =============================================================================
# KEY TEST: Sum all SBS weights per sample
# =============================================================================
banner("TEST: Do all samples sum to the same total?")

df['SBS_total_sum'] = df[sbs_cols].sum(axis=1)

log(f"  Per-sample SBS sum statistics:")
log(f"    n samples:     {len(df)}")
log(f"    Min sum:       {df['SBS_total_sum'].min():.1f}")
log(f"    Max sum:       {df['SBS_total_sum'].max():.1f}")
log(f"    Mean sum:      {df['SBS_total_sum'].mean():.1f}")
log(f"    Median sum:    {df['SBS_total_sum'].median():.1f}")
log(f"    Std sum:       {df['SBS_total_sum'].std():.1f}")
log(f"    Q25:           {df['SBS_total_sum'].quantile(0.25):.1f}")
log(f"    Q75:           {df['SBS_total_sum'].quantile(0.75):.1f}")

n_unique_sums = df['SBS_total_sum'].nunique()
log(f"\n    Unique sum values: {n_unique_sums}")

# Check if all sums are identical (constant normalization)
sum_range = df['SBS_total_sum'].max() - df['SBS_total_sum'].min()
cv = df['SBS_total_sum'].std() / df['SBS_total_sum'].mean() if df['SBS_total_sum'].mean() > 0 else 0
log(f"    Range (max - min): {sum_range:.1f}")
log(f"    Coefficient of variation: {cv:.4f}")

# =============================================================================
# SBS2-SPECIFIC ANALYSIS
# =============================================================================
banner("SBS2 VALUE DISTRIBUTION")

if 'SBS2' in df.columns:
    sbs2 = df['SBS2']
    log(f"  SBS2 statistics:")
    log(f"    Min:              {sbs2.min():.1f}")
    log(f"    Max:              {sbs2.max():.1f}")
    log(f"    Mean:             {sbs2.mean():.1f}")
    log(f"    Median:           {sbs2.median():.1f}")
    log(f"    Samples SBS2 > 0: {(sbs2 > 0).sum()}")
    log(f"    Samples SBS2 > 10:{(sbs2 > 10).sum()}")
    log(f"    Samples SBS2 > 50:{(sbs2 > 50).sum()}")
    log(f"    Samples SBS2 >100:{(sbs2 > 100).sum()}")

    # Show SBS2 as fraction of total per sample
    df['SBS2_fraction'] = df['SBS2'] / df['SBS_total_sum']
    log(f"\n  SBS2 as fraction of total per sample:")
    log(f"    Min fraction:     {df['SBS2_fraction'].min():.4f}")
    log(f"    Max fraction:     {df['SBS2_fraction'].max():.4f}")
    log(f"    Mean fraction:    {df['SBS2_fraction'].mean():.4f}")
    log(f"    Median fraction:  {df['SBS2_fraction'].median():.4f}")

    # Correlation between total mutations and SBS2
    from scipy.stats import spearmanr
    rho, p = spearmanr(df['SBS_total_sum'], df['SBS2'])
    log(f"\n  Spearman correlation (total SBS sum vs SBS2):")
    log(f"    rho = {rho:.4f}, p = {p:.2e}")
    log(f"    (If positive, high-mutation samples tend to have higher SBS2)")
else:
    log(f"  WARNING: 'SBS2' column not found.")
    log(f"  Available SBS columns with '2': {[c for c in sbs_cols if '2' in c]}")

# =============================================================================
# FILTER TO HNSC AND CHECK THERE TOO
# =============================================================================
banner("HNSC-SPECIFIC CHECK")

# Find sample ID column
sample_col = None
for candidate in ['TCGA_Gene_Expression_Entity_ID', 'Entity_ID', 'Sample', 'sample']:
    if candidate in df.columns:
        sample_col = candidate
        break

if sample_col:
    log(f"  Sample column: {sample_col}")

    # Try to identify HNSC samples
    project_col = None
    for candidate in ['Project_ID', 'project', 'cancer_type', 'Cancer_Type']:
        if candidate in df.columns:
            project_col = candidate
            break

    if project_col:
        hnsc_mask = df[project_col].str.contains('HNSC', case=False, na=False)
        log(f"  Project column: {project_col}")
    else:
        hnsc_mask = df[sample_col].str.contains('HNSC', case=False, na=False)
        log(f"  No project column found, matching HNSC in sample IDs")

    n_hnsc = hnsc_mask.sum()
    log(f"  HNSC samples found: {n_hnsc}")

    if n_hnsc > 0:
        hnsc_df = df[hnsc_mask]
        log(f"\n  HNSC SBS sum statistics:")
        log(f"    Min sum:       {hnsc_df['SBS_total_sum'].min():.1f}")
        log(f"    Max sum:       {hnsc_df['SBS_total_sum'].max():.1f}")
        log(f"    Mean sum:      {hnsc_df['SBS_total_sum'].mean():.1f}")
        log(f"    Std sum:       {hnsc_df['SBS_total_sum'].std():.1f}")
        log(f"    Unique sums:   {hnsc_df['SBS_total_sum'].nunique()}")

        hnsc_range = hnsc_df['SBS_total_sum'].max() - hnsc_df['SBS_total_sum'].min()
        hnsc_cv = hnsc_df['SBS_total_sum'].std() / hnsc_df['SBS_total_sum'].mean() if hnsc_df['SBS_total_sum'].mean() > 0 else 0
        log(f"    Range:         {hnsc_range:.1f}")
        log(f"    CV:            {hnsc_cv:.4f}")

        if 'SBS2' in df.columns:
            log(f"\n  HNSC SBS2 statistics:")
            log(f"    Min:           {hnsc_df['SBS2'].min():.1f}")
            log(f"    Max:           {hnsc_df['SBS2'].max():.1f}")
            log(f"    Mean:          {hnsc_df['SBS2'].mean():.1f}")
            log(f"    Samples > 0:   {(hnsc_df['SBS2'] > 0).sum()}")
            log(f"    Samples > 50:  {(hnsc_df['SBS2'] > 50).sum()}")
            log(f"    Samples > 100: {(hnsc_df['SBS2'] > 100).sum()}")
else:
    log(f"  No sample ID column found for HNSC filtering")

# =============================================================================
# EXAMPLE ROWS — show range of variation
# =============================================================================
banner("EXAMPLE ROWS (sorted by total sum)")

cols_to_show = []
if sample_col and sample_col in df.columns:
    cols_to_show.append(sample_col)
if 'SBS2' in df.columns:
    cols_to_show.append('SBS2')
cols_to_show.append('SBS_total_sum')

# Show lowest, median, and highest total sum samples
sorted_df = df.sort_values('SBS_total_sum')
n = len(sorted_df)
examples = pd.concat([
    sorted_df.head(5),
    sorted_df.iloc[n//2 - 2 : n//2 + 3],
    sorted_df.tail(5)
])

log(f"\n  {'Sample':30s}  {'SBS2':>8s}  {'Total_Sum':>10s}")
log(f"  {'-'*30}  {'-'*8}  {'-'*10}")
for _, row in examples.iterrows():
    sid = str(row[sample_col])[:30] if sample_col and sample_col in df.columns else "N/A"
    sbs2_val = f"{row['SBS2']:.0f}" if 'SBS2' in df.columns else "N/A"
    total = row['SBS_total_sum']
    log(f"  {sid:30s}  {sbs2_val:>8s}  {total:>10.1f}")

# =============================================================================
# VERDICT
# =============================================================================
banner("VERDICT")

if n_unique_sums <= 3:
    log(f"  CONSTANT NORMALIZATION DETECTED")
    log(f"  All {len(df)} samples sum to {df['SBS_total_sum'].iloc[0]:.1f} (+/- {sum_range:.1f})")
    log(f"  → The mutation counts were normalized to a constant total per sample")
    log(f"  → A sample with 8 mutations gets the same total as one with 800")
    log(f"  → Daiko's concern IS VALID: SBS2 weight is relative, not absolute")
    log(f"  → ACTION REQUIRED: Download VCFs, recompute from raw mutation counts")
    verdict = "CONSTANT_NORMALIZATION"
elif cv < 0.05:
    log(f"  NEAR-CONSTANT NORMALIZATION (CV = {cv:.4f})")
    log(f"  Sums cluster tightly: {df['SBS_total_sum'].min():.1f} to {df['SBS_total_sum'].max():.1f}")
    log(f"  → Likely normalized with very small sample-to-sample variation")
    log(f"  → Daiko's concern is LIKELY VALID")
    log(f"  → RECOMMEND: Download VCFs to verify and recompute if needed")
    verdict = "NEAR_CONSTANT"
elif cv > 0.3:
    log(f"  ABSOLUTE COUNTS CONFIRMED (CV = {cv:.4f})")
    log(f"  Sums vary widely: {df['SBS_total_sum'].min():.1f} to {df['SBS_total_sum'].max():.1f}")
    log(f"  {n_unique_sums} unique sum values across {len(df)} samples")
    log(f"  → Each sample's SBS weights reflect its actual context-specific mutation counts")
    log(f"  → Samples with more total mutations naturally have higher absolute SBS2 counts")
    log(f"  → Daiko's concern does NOT apply: these are already absolute, not normalized")
    log(f"  → NO VCF DOWNLOAD NEEDED for SBS2 weight recomputation")
    log(f"  → Current SBS2 values are suitable for between-sample comparison")
    verdict = "ABSOLUTE_COUNTS"
else:
    log(f"  AMBIGUOUS (CV = {cv:.4f})")
    log(f"  Sums show moderate variation: {df['SBS_total_sum'].min():.1f} to {df['SBS_total_sum'].max():.1f}")
    log(f"  → Manual inspection recommended")
    log(f"  → Consider downloading a few VCFs to spot-check against the weights")
    verdict = "AMBIGUOUS"

# =============================================================================
# SAVE OUTPUT
# =============================================================================
banner("SAVE OUTPUT")

output_path = os.path.join(OUTPUT_DIR, "sbs_weight_normalization_check.tsv")
save_cols = []
if sample_col and sample_col in df.columns:
    save_cols.append(sample_col)
if 'SBS2' in df.columns:
    save_cols.append('SBS2')
save_cols.append('SBS_total_sum')
if 'SBS2_fraction' in df.columns:
    save_cols.append('SBS2_fraction')

df[save_cols].to_csv(output_path, sep='\t', index=False)
log(f"  Saved: {output_path}")

report_path = os.path.join(OUTPUT_DIR, "sbs_weight_normalization_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")

banner(f"DIAGNOSTIC COMPLETE — VERDICT: {verdict}")
