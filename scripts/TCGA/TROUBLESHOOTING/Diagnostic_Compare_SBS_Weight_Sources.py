#!/usr/bin/env python3
"""
Diagnostic_Compare_SBS_Weight_Sources_v3.py
=============================================
Compare original Mutation_Table_Tumors_TCGA.tsv against SigProfiler outputs.

v3 changes:
  - Uses the EXPLICIT barcode crosswalk already in the original file:
      TCGA_Gene_Expression_Entity_ID = RNA-seq aliquot barcode
      Mutation_Signature__File_Orginal_Entity_ID = WES aliquot barcode
  - Uses TCGA_sample_metadata_final.tsv for Project_ID (cancer type) filtering
  - Also falls back to 16-char truncation and 12-char Case_ID matching
    to maximize overlap and identify any unmatchable samples

Inputs:
  data/FIG_1/Mutation_Table_Tumors_TCGA.tsv
  data/FIG_1/TCGA_sample_metadata_final.tsv
  SHARED/TCGA/VCF/SigProfiler_output/TCGA_SBS_signature_weights.tsv
  SHARED/TCGA/VCF/SigProfiler_output/TCGA_SBS_signature_counts.tsv

Outputs (data/FIG_1/TROUBLESHOOTING/):
  sbs_v3_comparison_report.txt
  sbs_v3_barcode_crosswalk.tsv
  sbs_v3_matched_all.tsv
  sbs_v3_matched_HNSC.tsv

Outputs (data/FIG_1/):
  original_HNSC_SBS_weights.tsv
  new_counts_HNSC_SBS_weights.tsv
  new_fractions_HNSC_SBS_weights.tsv

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

ORIGINAL_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
METADATA_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_sample_metadata_final.tsv")
NEW_WEIGHTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_weights.tsv")
NEW_COUNTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_counts.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TROUBLESHOOTING")
HNSC_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1")
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

def load_sigprofiler_file(path, label):
    """Load SigProfiler file, auto-detecting transposed format."""
    log(f"\n  Loading {label}: {os.path.basename(path)}")
    raw = pd.read_csv(path, sep='\t')
    log(f"    Raw shape: {raw.shape}")

    first_col = raw.columns[0]
    first_col_vals = raw[first_col].astype(str).tolist()

    if any(v.startswith('SBS') for v in first_col_vals):
        log(f"    Transposed format — transposing...")
        raw = raw.set_index(first_col)
        df = raw.T.copy()
        df.index.name = 'Tumor_Sample_Barcode'
        df = df.reset_index()

        # Drop non-sample rows (e.g. Cancer_Type header that got transposed)
        bad = df['Tumor_Sample_Barcode'].str.contains('Cancer_Type|^$', na=True, regex=True)
        if bad.any():
            df = df[~bad].copy()

        sbs_cols = [c for c in df.columns if c.startswith('SBS')]
        for c in sbs_cols:
            df[c] = pd.to_numeric(df[c], errors='coerce')

        log(f"    After transpose: {df.shape} ({len(sbs_cols)} SBS sigs)")
        return df, 'Tumor_Sample_Barcode', sbs_cols
    else:
        sbs_cols = [c for c in raw.columns if c.startswith('SBS')]
        log(f"    Standard format: {raw.shape} ({len(sbs_cols)} SBS sigs)")
        return raw, first_col, sbs_cols

# =============================================================================
# STEP 1: LOAD ALL FILES
# =============================================================================
banner("STEP 1: Load all files")

log(f"\n  Loading ORIGINAL: {os.path.basename(ORIGINAL_PATH)}")
orig = pd.read_csv(ORIGINAL_PATH, sep='\t')
orig_rnaseq_col = 'TCGA_Gene_Expression_Entity_ID'
orig_wes_col = 'Mutation_Signature__File_Orginal_Entity_ID'
orig_sbs_cols = [c for c in orig.columns if c.startswith('SBS')]
log(f"    Shape: {orig.shape}, {len(orig_sbs_cols)} SBS columns")
log(f"    RNA-seq barcode col: {orig_rnaseq_col}")
log(f"    WES barcode col:     {orig_wes_col}")
log(f"    Example RNA-seq: {orig[orig_rnaseq_col].iloc[0]}")
log(f"    Example WES:     {orig[orig_wes_col].iloc[0]}")

log(f"\n  Loading METADATA: {os.path.basename(METADATA_PATH)}")
metadata = pd.read_csv(METADATA_PATH, sep='\t')
log(f"    Shape: {metadata.shape}")
log(f"    Columns: {list(metadata.columns)}")
log(f"    Cancer types: {metadata['Project_ID'].nunique()}")

new_wt, new_wt_sample_col, new_wt_sbs_cols = load_sigprofiler_file(NEW_WEIGHTS_PATH, "NEW WEIGHTS")
new_ct, new_ct_sample_col, new_ct_sbs_cols = load_sigprofiler_file(NEW_COUNTS_PATH, "NEW COUNTS")

# =============================================================================
# STEP 2: BUILD CROSSWALK USING EXPLICIT WES BARCODE
# =============================================================================
banner("STEP 2: Match using explicit WES barcode crosswalk")

# The original file already contains the WES barcode in the 4th column
# New SigProfiler files use WES barcodes as sample IDs
orig_wes_ids = set(orig[orig_wes_col].astype(str).values)
new_ct_ids = set(new_ct[new_ct_sample_col].astype(str).values)
new_wt_ids = set(new_wt[new_wt_sample_col].astype(str).values)

log(f"  Original WES barcodes:     {len(orig_wes_ids)}")
log(f"  New counts barcodes:       {len(new_ct_ids)}")
log(f"  New weights barcodes:      {len(new_wt_ids)}")

# Direct match on WES barcode
direct_overlap = orig_wes_ids & new_ct_ids
log(f"\n  DIRECT MATCH (WES barcode): {len(direct_overlap)}")

# If direct match is low, try fallback strategies
if len(direct_overlap) < len(orig_wes_ids) * 0.5:
    log(f"\n  Direct match < 50%, trying fallback strategies...")

    # Fallback 1: 16-char truncation
    orig_trunc16 = {s[:16] for s in orig_wes_ids}
    ct_trunc16 = {s[:16] for s in new_ct_ids}
    log(f"  Fallback 16-char truncation: {len(orig_trunc16 & ct_trunc16)}")

    # Fallback 2: 12-char Case_ID
    orig_trunc12 = {s[:12] for s in orig_wes_ids}
    ct_trunc12 = {s[:12] for s in new_ct_ids}
    log(f"  Fallback 12-char Case_ID:    {len(orig_trunc12 & ct_trunc12)}")

    # Show example mismatches
    unmatched = list(orig_wes_ids - new_ct_ids)[:5]
    log(f"\n  Example unmatched WES barcodes from original:")
    for bc in unmatched:
        # Find closest match in new
        candidates = [n for n in new_ct_ids if n[:12] == bc[:12]]
        log(f"    {bc}  →  candidates: {candidates[:3]}")

# Use the best matching strategy
if len(direct_overlap) > 100:
    match_method = "direct_wes"
    log(f"\n  Using DIRECT WES barcode matching ({len(direct_overlap)} matches)")
else:
    match_method = "trunc16"
    log(f"\n  Using 16-char truncation matching")

# =============================================================================
# STEP 3: MERGE AND COMPARE
# =============================================================================
banner("STEP 3: Merge and compare SBS2 values")

if match_method == "direct_wes":
    # Merge on WES barcode directly
    m_orig = orig[[orig_rnaseq_col, orig_wes_col, 'SBS2']].rename(
        columns={'SBS2': 'SBS2_original', orig_wes_col: 'wes_barcode', orig_rnaseq_col: 'rnaseq_barcode'})
    m_ct = new_ct[[new_ct_sample_col, 'SBS2']].rename(
        columns={new_ct_sample_col: 'wes_barcode', 'SBS2': 'SBS2_new_count'})
    m_wt = new_wt[[new_wt_sample_col, 'SBS2']].rename(
        columns={new_wt_sample_col: 'wes_barcode', 'SBS2': 'SBS2_new_fraction'})

    merged = m_orig.merge(m_ct, on='wes_barcode', how='inner')
    merged = merged.merge(m_wt, on='wes_barcode', how='inner')
else:
    # 16-char truncation
    orig['trunc16'] = orig[orig_wes_col].astype(str).str[:16]
    new_ct['trunc16'] = new_ct[new_ct_sample_col].astype(str).str[:16]
    new_wt['trunc16'] = new_wt[new_wt_sample_col].astype(str).str[:16]

    m_orig = orig[[orig_rnaseq_col, orig_wes_col, 'trunc16', 'SBS2']].rename(
        columns={'SBS2': 'SBS2_original', orig_wes_col: 'wes_barcode', orig_rnaseq_col: 'rnaseq_barcode'})
    m_ct = new_ct[['trunc16', 'SBS2']].rename(columns={'SBS2': 'SBS2_new_count'})
    m_wt = new_wt[['trunc16', 'SBS2']].rename(columns={'SBS2': 'SBS2_new_fraction'})

    # Deduplicate on trunc16
    m_orig = m_orig.drop_duplicates(subset='trunc16', keep='first')
    m_ct = m_ct.drop_duplicates(subset='trunc16', keep='first')
    m_wt = m_wt.drop_duplicates(subset='trunc16', keep='first')

    merged = m_orig.merge(m_ct, on='trunc16', how='inner')
    merged = merged.merge(m_wt, on='trunc16', how='inner')

log(f"  Matched samples: {len(merged)}")

# Add cancer type from metadata
meta_map = metadata.set_index('Entity_ID')['Project_ID'].to_dict()
merged['cancer_type'] = merged['rnaseq_barcode'].map(meta_map)
n_annotated = merged['cancer_type'].notna().sum()
log(f"  With cancer type (from RNA-seq metadata): {n_annotated}")

# If metadata mapping missed some, try Case_ID mapping
if n_annotated < len(merged):
    missing = merged['cancer_type'].isna()
    log(f"  Missing cancer type: {missing.sum()} — trying Case_ID fallback")
    case_id_map = metadata.drop_duplicates(subset='Case_ID').set_index('Case_ID')['Project_ID'].to_dict()
    merged.loc[missing, 'cancer_type'] = merged.loc[missing, 'rnaseq_barcode'].str[:12].map(case_id_map)
    log(f"  After Case_ID fallback: {merged['cancer_type'].notna().sum()} annotated")

# --- Correlations ---
log(f"\n  Correlations (n={len(merged)}):")
for pair_name, col_a, col_b in [
    ("Original vs New_Count", 'SBS2_original', 'SBS2_new_count'),
    ("Original vs New_Fraction", 'SBS2_original', 'SBS2_new_fraction'),
    ("New_Count vs New_Fraction", 'SBS2_new_count', 'SBS2_new_fraction'),
]:
    a, b = merged[col_a], merged[col_b]
    rho, rho_p = spearmanr(a, b)
    r, r_p = pearsonr(a, b)
    log(f"    {pair_name}:")
    log(f"      Spearman rho = {rho:.4f} (p = {rho_p:.2e})")
    log(f"      Pearson r    = {r:.4f} (p = {r_p:.2e})")

# --- Ratio analysis ---
log(f"\n  Ratio analysis (Original / New_Count):")
both_pos = (merged['SBS2_original'] > 0) & (merged['SBS2_new_count'] > 0)
if both_pos.sum() > 10:
    ratios = merged.loc[both_pos, 'SBS2_original'] / merged.loc[both_pos, 'SBS2_new_count']
    log(f"    n pairs (both > 0): {both_pos.sum()}")
    log(f"    Ratio min:          {ratios.min():.4f}")
    log(f"    Ratio max:          {ratios.max():.4f}")
    log(f"    Ratio mean:         {ratios.mean():.4f}")
    log(f"    Ratio median:       {ratios.median():.4f}")
    log(f"    Ratio std:          {ratios.std():.4f}")
    ratio_cv = ratios.std() / ratios.mean() if ratios.mean() > 0 else float('inf')
    log(f"    Ratio CV:           {ratio_cv:.4f}")

    if ratio_cv < 0.05:
        log(f"    → CONSTANT RATIO: Original ≈ {ratios.median():.2f} × New_Count")
    elif ratio_cv < 0.2:
        log(f"    → NEAR-CONSTANT RATIO")
    else:
        log(f"    → VARIABLE RATIO: different decompositions")

# --- Rounding check ---
log(f"\n  Rounding check:")
rounded = merged['SBS2_new_count'].round(0)
exact = (merged['SBS2_original'] == rounded).sum()
within_1 = ((merged['SBS2_original'] - rounded).abs() <= 1).sum()
within_5 = ((merged['SBS2_original'] - rounded).abs() <= 5).sum()
log(f"    Exact match:   {exact} / {len(merged)} ({100*exact/len(merged):.1f}%)")
log(f"    Within +/-1:   {within_1} / {len(merged)} ({100*within_1/len(merged):.1f}%)")
log(f"    Within +/-5:   {within_5} / {len(merged)} ({100*within_5/len(merged):.1f}%)")

# --- Zero agreement ---
log(f"\n  Zero agreement:")
orig_zero = merged['SBS2_original'] == 0
ct_zero = merged['SBS2_new_count'] == 0
log(f"    Both zero:           {(orig_zero & ct_zero).sum()}")
log(f"    Both non-zero:       {(~orig_zero & ~ct_zero).sum()}")
log(f"    Original zero only:  {(orig_zero & ~ct_zero).sum()}")
log(f"    New count zero only: {(~orig_zero & ct_zero).sum()}")

# --- Rank agreement ---
log(f"\n  Rank agreement:")
rho_ct, _ = spearmanr(merged['SBS2_original'], merged['SBS2_new_count'])
rho_wt, _ = spearmanr(merged['SBS2_original'], merged['SBS2_new_fraction'])
log(f"    Original vs New_Count:    rho = {rho_ct:.4f}")
log(f"    Original vs New_Fraction: rho = {rho_wt:.4f}")

# --- Examples across the range ---
log(f"\n  Example rows (sorted by original SBS2):")
examples = merged.sort_values('SBS2_original')
n = len(examples)
idx = sorted(set([0, 1, 2, n//4, n//2, 3*n//4, n-3, n-2, n-1]))
idx = [i for i in idx if 0 <= i < n]
log(f"    {'WES_Barcode':30s} {'Cancer':>8s}  {'Orig':>8s}  {'NewCt':>8s}  {'NewFrac':>10s}  {'Ratio':>7s}")
log(f"    {'-'*30} {'-'*8}  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*7}")
for i in idx:
    row = examples.iloc[i]
    wes = str(row.get('wes_barcode', ''))[:30]
    ct = str(row.get('cancer_type', ''))[:8]
    ratio = row['SBS2_original'] / row['SBS2_new_count'] if row['SBS2_new_count'] > 0 else float('nan')
    log(f"    {wes:30s} {ct:>8s}  {row['SBS2_original']:8.0f}  {row['SBS2_new_count']:8.0f}  "
        f"{row['SBS2_new_fraction']:10.6f}  {ratio:7.2f}")

# Save
merged.to_csv(os.path.join(OUTPUT_DIR, "sbs_v3_matched_all.tsv"), sep='\t', index=False)
log(f"\n  Saved: sbs_v3_matched_all.tsv ({len(merged)} samples)")

# Save crosswalk
cw_cols = ['rnaseq_barcode', 'wes_barcode', 'cancer_type']
cw_cols = [c for c in cw_cols if c in merged.columns]
merged[cw_cols].to_csv(os.path.join(OUTPUT_DIR, "sbs_v3_barcode_crosswalk.tsv"), sep='\t', index=False)
log(f"  Saved: sbs_v3_barcode_crosswalk.tsv")

# =============================================================================
# STEP 4: FILTER TO HNSC
# =============================================================================
banner("STEP 4: Filter to HNSC")

hnsc = merged[merged['cancer_type'] == 'TCGA-HNSC'].copy()
log(f"  HNSC matched samples: {len(hnsc)}")

if len(hnsc) > 0:
    log(f"\n  HNSC SBS2 comparison:")
    log(f"    Original:  min={hnsc['SBS2_original'].min():.0f}  max={hnsc['SBS2_original'].max():.0f}  "
        f"mean={hnsc['SBS2_original'].mean():.1f}  n>0={(hnsc['SBS2_original']>0).sum()}")
    log(f"    New count: min={hnsc['SBS2_new_count'].min():.0f}  max={hnsc['SBS2_new_count'].max():.0f}  "
        f"mean={hnsc['SBS2_new_count'].mean():.1f}  n>0={(hnsc['SBS2_new_count']>0).sum()}")
    log(f"    New frac:  min={hnsc['SBS2_new_fraction'].min():.6f}  max={hnsc['SBS2_new_fraction'].max():.6f}  "
        f"mean={hnsc['SBS2_new_fraction'].mean():.6f}  n>0={(hnsc['SBS2_new_fraction']>0).sum()}")

    rho_hnsc, p_hnsc = spearmanr(hnsc['SBS2_original'], hnsc['SBS2_new_count'])
    log(f"\n  HNSC Spearman (Original vs New_Count): rho={rho_hnsc:.4f}, p={p_hnsc:.2e}")
    rho_hnsc2, p_hnsc2 = spearmanr(hnsc['SBS2_original'], hnsc['SBS2_new_fraction'])
    log(f"  HNSC Spearman (Original vs New_Frac):  rho={rho_hnsc2:.4f}, p={p_hnsc2:.2e}")

    # Save HNSC comparison
    hnsc.to_csv(os.path.join(OUTPUT_DIR, "sbs_v3_matched_HNSC.tsv"), sep='\t', index=False)
    log(f"\n  Saved: sbs_v3_matched_HNSC.tsv ({len(hnsc)} samples)")

    # Save HNSC-only files for downstream Figure 1 use
    hnsc_rnaseq_ids = set(hnsc['rnaseq_barcode'].values)
    hnsc_wes_ids = set(hnsc['wes_barcode'].values) if 'wes_barcode' in hnsc.columns else set()

    orig_hnsc = orig[orig[orig_rnaseq_col].isin(hnsc_rnaseq_ids)].copy()
    orig_hnsc.to_csv(os.path.join(HNSC_OUTPUT_DIR, "original_HNSC_SBS_weights.tsv"), sep='\t', index=False)
    log(f"  Saved: original_HNSC_SBS_weights.tsv ({len(orig_hnsc)} samples)")

    new_ct_hnsc = new_ct[new_ct[new_ct_sample_col].isin(hnsc_wes_ids)].copy()
    new_ct_hnsc.to_csv(os.path.join(HNSC_OUTPUT_DIR, "new_counts_HNSC_SBS_weights.tsv"), sep='\t', index=False)
    log(f"  Saved: new_counts_HNSC_SBS_weights.tsv ({len(new_ct_hnsc)} samples)")

    new_wt_hnsc = new_wt[new_wt[new_wt_sample_col].isin(hnsc_wes_ids)].copy()
    new_wt_hnsc.to_csv(os.path.join(HNSC_OUTPUT_DIR, "new_fractions_HNSC_SBS_weights.tsv"), sep='\t', index=False)
    log(f"  Saved: new_fractions_HNSC_SBS_weights.tsv ({len(new_wt_hnsc)} samples)")

else:
    log(f"  No HNSC samples found. Check cancer_type values:")
    log(f"    Unique cancer types: {sorted(merged['cancer_type'].dropna().unique())[:20]}")

# =============================================================================
# STEP 5: VERDICT
# =============================================================================
banner("VERDICT")

log(f"""
  MATCHING: Used {'direct WES barcode' if match_method == 'direct_wes' else '16-char truncation'}.
  Matched {len(merged)} samples across original and new SigProfiler files.

  PAN-CANCER SBS2 RANK CORRELATION:
    Original vs New_Count:    rho = {rho_ct:.4f}
    Original vs New_Fraction: rho = {rho_wt:.4f}

  INTERPRETATION:
    Original file: {len(orig_sbs_cols)} signatures (older COSMIC version, integer counts)
    New file:      {len(new_ct_sbs_cols)} signatures (COSMIC v3.4, integer counts)

    Both are absolute mutation count decompositions. They differ because:
      - Different COSMIC signature versions (65 vs 86 signatures)
      - More signatures in v3.4 means mutations get split across more categories
      - This can change per-signature counts but should preserve relative ranking

  FOR FIGURE 1:
    We recommend using the NEW COUNTS file because:
      1. Documented provenance (SigProfilerAssignment v3.4 on your own VCFs)
      2. Current COSMIC version (v3.4)
      3. Absolute counts, not fractions, which is what Daiko wanted
    Note in methods: "SBS mutational signature weights were computed from
    MuTect2-called somatic SNVs using SigProfilerAssignment (COSMIC v3.4,
    exome mode), replacing previously used pre-computed signature weights."

  HNSC FILES: {len(hnsc) if len(hnsc) > 0 else 0} samples saved to data/FIG_1/
""")

# Save report
report_path = os.path.join(OUTPUT_DIR, "sbs_v3_comparison_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")

banner("DIAGNOSTIC COMPLETE")
