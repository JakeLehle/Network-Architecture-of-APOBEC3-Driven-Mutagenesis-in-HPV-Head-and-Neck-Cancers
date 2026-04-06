#!/usr/bin/env python3
"""
Phase2_Raw_HPV16_Counts.py
===========================
Extract raw HPV16 UMI counts from Kraken2 per-sample matrices.

The adata_v_pp.h5ad has CPM+log1p normalized counts that saturate at 13.82,
making HPV16 effectively binary. This script reads the raw 10x-format
Kraken2 matrices to recover actual UMI counts per cell.

Inputs:
    - Per-sample kraken2_filtered_feature_bc_matrix/ directories
    - data/FIG_6/00_diagnostic_inventory/basal_cell_master_table.tsv (from Phase 1)

Outputs (to data/FIG_6/01_raw_hpv16_counts/):
    - raw_HPV16_counts_all_cells.tsv     (barcode, raw_HPV16, sample)
    - raw_HPV16_counts_basal_only.tsv    (merged with master table)
    - hpv16_raw_distribution_report.txt  (diagnostic)

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import mannwhitneyu, spearmanr
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
KRAKEN2_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"

# Master table from Phase 1
MASTER_TABLE_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/00_diagnostic_inventory/basal_cell_master_table.tsv")

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/01_raw_hpv16_counts")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# HPV16 search terms (Kraken2 may use taxid or name)
HPV16_NAMES = ['Human papillomavirus 16', 'Human papillomavirus type 16',
               'human papillomavirus 16', 'HPV16', 'HPV-16']
HPV16_TAXID = '333760'  # NCBI taxonomy ID for HPV16

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
# STEP 1: DISCOVER SAMPLE DIRECTORIES
# =============================================================================
log_sep("STEP 1: Discover Kraken2 sample directories")

# Find all kraken2 matrix directories
matrix_dirs = sorted(glob.glob(os.path.join(KRAKEN2_BASE, "SRR*", "SRR*_S1_L001_", "outs", "kraken2_filtered_feature_bc_matrix")))
log(f"  Found {len(matrix_dirs)} sample directories")

# Extract SRR IDs
sample_info = []
for d in matrix_dirs:
    parts = d.split('/')
    # Find the SRR ID from the path
    for p in parts:
        if p.startswith('SRR'):
            srr_id = p
            break
    sample_info.append({'srr_id': srr_id, 'matrix_dir': d})

srr_ids = [s['srr_id'] for s in sample_info]
log(f"  SRR IDs: {srr_ids[:5]}... ({len(srr_ids)} total)")

# =============================================================================
# STEP 2: PEEK AT FIRST SAMPLE TO UNDERSTAND FORMAT
# =============================================================================
log_sep("STEP 2: Inspect first sample matrix format")

first_dir = sample_info[0]['matrix_dir']
log(f"  Directory: {first_dir}")

# List contents
contents = os.listdir(first_dir)
log(f"  Contents: {contents}")

# Check for standard 10x files
for fname in ['matrix.mtx', 'matrix.mtx.gz', 'barcodes.tsv', 'barcodes.tsv.gz',
              'features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz']:
    fpath = os.path.join(first_dir, fname)
    if os.path.exists(fpath):
        size = os.path.getsize(fpath)
        log(f"  [FOUND] {fname} ({size/1024:.1f} KB)")

# Read features/genes file
features_path = None
for fname in ['features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv']:
    fpath = os.path.join(first_dir, fname)
    if os.path.exists(fpath):
        features_path = fpath
        break

if features_path:
    if features_path.endswith('.gz'):
        features_df = pd.read_csv(features_path, sep='\t', header=None, compression='gzip')
    else:
        features_df = pd.read_csv(features_path, sep='\t', header=None)
    
    log(f"\n  Features file: {features_path}")
    log(f"  Shape: {features_df.shape}")
    log(f"  Columns: {features_df.columns.tolist()}")
    log(f"  First 5 rows:")
    for _, row in features_df.head().iterrows():
        log(f"    {row.tolist()}")
    log(f"  Last 5 rows:")
    for _, row in features_df.tail().iterrows():
        log(f"    {row.tolist()}")
    
    # Search for HPV16
    log(f"\n  Searching for HPV16 in features...")
    for col in features_df.columns:
        if features_df[col].dtype == object:
            hpv_matches = features_df[features_df[col].str.contains('papilloma|HPV|hpv', case=False, na=False)]
            if len(hpv_matches) > 0:
                log(f"  Column {col} HPV matches:")
                for _, row in hpv_matches.iterrows():
                    log(f"    Row {row.name}: {row.tolist()}")
    
    # Also check for taxid match
    for col in features_df.columns:
        taxid_matches = features_df[features_df[col].astype(str) == HPV16_TAXID]
        if len(taxid_matches) > 0:
            log(f"  Taxid {HPV16_TAXID} found in column {col}:")
            for _, row in taxid_matches.iterrows():
                log(f"    Row {row.name}: {row.tolist()}")

# Read barcodes
barcodes_path = None
for fname in ['barcodes.tsv.gz', 'barcodes.tsv']:
    fpath = os.path.join(first_dir, fname)
    if os.path.exists(fpath):
        barcodes_path = fpath
        break

if barcodes_path:
    if barcodes_path.endswith('.gz'):
        barcodes_df = pd.read_csv(barcodes_path, sep='\t', header=None, compression='gzip')
    else:
        barcodes_df = pd.read_csv(barcodes_path, sep='\t', header=None)
    
    log(f"\n  Barcodes file: {barcodes_path}")
    log(f"  Total barcodes: {len(barcodes_df)}")
    log(f"  First 3: {barcodes_df[0].head(3).tolist()}")
    log(f"  Format example: '{barcodes_df[0].iloc[0]}'")

# Read matrix
matrix_path = None
for fname in ['matrix.mtx.gz', 'matrix.mtx']:
    fpath = os.path.join(first_dir, fname)
    if os.path.exists(fpath):
        matrix_path = fpath
        break

if matrix_path:
    mat = mmread(matrix_path)
    log(f"\n  Matrix: {mat.shape} (features x cells), nnz={mat.nnz}")

# =============================================================================
# STEP 3: EXTRACT HPV16 RAW COUNTS FROM ALL SAMPLES
# =============================================================================
log_sep("STEP 3: Extract HPV16 raw counts from all samples")

all_hpv16_counts = []
sample_summaries = []

for sample in sample_info:
    srr_id = sample['srr_id']
    mdir = sample['matrix_dir']
    
    # --- Read features ---
    feat_path = None
    for fname in ['features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv']:
        fpath = os.path.join(mdir, fname)
        if os.path.exists(fpath):
            feat_path = fpath
            break
    
    if feat_path is None:
        log(f"  [{srr_id}] WARNING: No features file found, skipping")
        continue
    
    if feat_path.endswith('.gz'):
        feats = pd.read_csv(feat_path, sep='\t', header=None, compression='gzip')
    else:
        feats = pd.read_csv(feat_path, sep='\t', header=None)
    
    # --- Find HPV16 feature index ---
    hpv16_idx = None
    hpv16_name_found = None
    
    # Search by name in all string columns
    for col in feats.columns:
        if feats[col].dtype == object:
            for name in HPV16_NAMES:
                matches = feats[feats[col] == name]
                if len(matches) > 0:
                    hpv16_idx = matches.index[0]
                    hpv16_name_found = name
                    break
        if hpv16_idx is not None:
            break
    
    # Search by taxid
    if hpv16_idx is None:
        for col in feats.columns:
            matches = feats[feats[col].astype(str) == HPV16_TAXID]
            if len(matches) > 0:
                hpv16_idx = matches.index[0]
                hpv16_name_found = f"taxid:{HPV16_TAXID}"
                break
    
    if hpv16_idx is None:
        log(f"  [{srr_id}] HPV16 not found in features, skipping")
        sample_summaries.append({'srr_id': srr_id, 'status': 'no_HPV16_feature',
                                 'n_cells': 0, 'n_hpv16_pos': 0})
        continue
    
    # --- Read barcodes ---
    bc_path = None
    for fname in ['barcodes.tsv.gz', 'barcodes.tsv']:
        fpath = os.path.join(mdir, fname)
        if os.path.exists(fpath):
            bc_path = fpath
            break
    
    if bc_path.endswith('.gz'):
        barcodes = pd.read_csv(bc_path, sep='\t', header=None, compression='gzip')[0].tolist()
    else:
        barcodes = pd.read_csv(bc_path, sep='\t', header=None)[0].tolist()
    
    # --- Read matrix ---
    mat_path = None
    for fname in ['matrix.mtx.gz', 'matrix.mtx']:
        fpath = os.path.join(mdir, fname)
        if os.path.exists(fpath):
            mat_path = fpath
            break
    
    mat = mmread(mat_path).tocsr()  # features x cells
    
    # --- Extract HPV16 row ---
    hpv16_counts = np.array(mat[hpv16_idx, :].todense()).flatten()
    
    # --- Build barcode mapping (append SRR suffix to match adata_final) ---
    # adata_final format: AAACCTGAGAGCCTAG-1-SRR14340883
    # Kraken2 barcodes might be: AAACCTGAGAGCCTAG-1
    full_barcodes = []
    for bc in barcodes:
        if srr_id not in bc:
            full_barcodes.append(f"{bc}-{srr_id}")
        else:
            full_barcodes.append(bc)
    
    # Store results
    n_pos = (hpv16_counts > 0).sum()
    for bc, count in zip(full_barcodes, hpv16_counts):
        if count > 0:  # only store positive cells to save memory
            all_hpv16_counts.append({'barcode': bc, 'raw_HPV16': int(count), 'sample': srr_id})
    
    # Also store zeros for summary
    summary = {
        'srr_id': srr_id,
        'status': 'ok',
        'hpv16_feature_idx': hpv16_idx,
        'hpv16_feature_name': hpv16_name_found,
        'n_cells': len(barcodes),
        'n_hpv16_pos': n_pos,
        'pct_pos': 100 * n_pos / len(barcodes) if len(barcodes) > 0 else 0,
        'max_raw_count': int(hpv16_counts.max()) if n_pos > 0 else 0,
        'mean_pos_count': float(hpv16_counts[hpv16_counts > 0].mean()) if n_pos > 0 else 0,
        'median_pos_count': float(np.median(hpv16_counts[hpv16_counts > 0])) if n_pos > 0 else 0,
        'total_hpv16_reads': int(hpv16_counts.sum()),
    }
    sample_summaries.append(summary)
    
    log(f"  [{srr_id}] cells={len(barcodes)}, HPV16+={n_pos} ({summary['pct_pos']:.1f}%), "
        f"max={summary['max_raw_count']}, mean(+)={summary['mean_pos_count']:.1f}, "
        f"total_reads={summary['total_hpv16_reads']}")

# =============================================================================
# STEP 4: COMBINE AND REPORT
# =============================================================================
log_sep("STEP 4: Combined raw HPV16 count statistics")

hpv16_df = pd.DataFrame(all_hpv16_counts)
log(f"  Total HPV16+ cells across all samples: {len(hpv16_df)}")

if len(hpv16_df) > 0:
    raw_vals = hpv16_df['raw_HPV16']
    log(f"\n  Raw HPV16 UMI count distribution (among HPV16+ cells):")
    log(f"    Mean:   {raw_vals.mean():.2f}")
    log(f"    Median: {raw_vals.median():.1f}")
    log(f"    Std:    {raw_vals.std():.2f}")
    log(f"    Min:    {raw_vals.min()}")
    log(f"    Max:    {raw_vals.max()}")
    for pct in [10, 25, 50, 75, 90, 95, 99]:
        log(f"    {pct}th:   {np.percentile(raw_vals, pct):.0f}")
    
    # Value counts for low end
    log(f"\n  Raw count value frequencies (first 20):")
    vc = raw_vals.value_counts().sort_index()
    for val, count in vc.head(20).items():
        log(f"    {val:5d} reads: {count:6d} cells ({100*count/len(hpv16_df):.1f}%)")
    if len(vc) > 20:
        log(f"    ... ({len(vc)} unique values total)")

# Per-sample summary table
log(f"\n  Per-sample summary:")
summary_df = pd.DataFrame(sample_summaries)
summary_df = summary_df.sort_values('n_hpv16_pos', ascending=False)
log(f"  {'SRR':15s} {'Cells':>6s} {'HPV16+':>7s} {'%Pos':>6s} {'Max':>5s} "
    f"{'Mean(+)':>8s} {'Med(+)':>7s} {'TotalReads':>11s}")
log(f"  {'-'*15} {'-'*6} {'-'*7} {'-'*6} {'-'*5} {'-'*8} {'-'*7} {'-'*11}")
for _, row in summary_df.iterrows():
    if row['status'] == 'ok':
        log(f"  {row['srr_id']:15s} {row['n_cells']:6d} {row['n_hpv16_pos']:7d} "
            f"{row['pct_pos']:5.1f}% {row['max_raw_count']:5d} "
            f"{row['mean_pos_count']:8.1f} {row['median_pos_count']:7.1f} "
            f"{row['total_hpv16_reads']:11d}")
    else:
        log(f"  {row['srr_id']:15s} {row['status']}")

summary_df.to_csv(os.path.join(OUTPUT_DIR, "per_sample_hpv16_summary.tsv"),
                   sep='\t', index=False)

# =============================================================================
# STEP 5: MAP RAW COUNTS TO BASAL CELL MASTER TABLE
# =============================================================================
log_sep("STEP 5: Map raw HPV16 counts to basal cell master table")

# Load master table
master = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
log(f"  Master table: {master.shape}")

# Create raw count series (default 0 for cells not in hpv16_df)
raw_series = hpv16_df.set_index('barcode')['raw_HPV16']
master['raw_HPV16'] = raw_series.reindex(master.index).fillna(0).astype(int)
master['HPV16_raw_positive'] = (master['raw_HPV16'] > 0).astype(int)

n_pos_basal = master['HPV16_raw_positive'].sum()
log(f"  HPV16+ basal cells (raw): {n_pos_basal} / {len(master)} ({100*n_pos_basal/len(master):.1f}%)")

# Compare with normalized values
if 'HPV16_positive' in master.columns:
    both_pos = ((master['HPV16_raw_positive'] == 1) & (master['HPV16_positive'] == 1)).sum()
    raw_only = ((master['HPV16_raw_positive'] == 1) & (master['HPV16_positive'] == 0)).sum()
    norm_only = ((master['HPV16_raw_positive'] == 0) & (master['HPV16_positive'] == 1)).sum()
    log(f"\n  Concordance with normalized HPV16 calls:")
    log(f"    Both positive:     {both_pos}")
    log(f"    Raw+ / Norm-:      {raw_only}")
    log(f"    Raw- / Norm+:      {norm_only}")

# Distribution of raw counts among basal HPV16+ cells
pos_raw = master.loc[master['HPV16_raw_positive'] == 1, 'raw_HPV16']
if len(pos_raw) > 0:
    log(f"\n  Raw HPV16 UMI distribution (HPV16+ basal cells, n={len(pos_raw)}):")
    log(f"    Mean:   {pos_raw.mean():.2f}")
    log(f"    Median: {pos_raw.median():.1f}")
    log(f"    Std:    {pos_raw.std():.2f}")
    log(f"    Min:    {pos_raw.min()}")
    log(f"    Max:    {pos_raw.max()}")
    for pct in [10, 25, 50, 75, 90, 95, 99]:
        log(f"    {pct}th:   {np.percentile(pos_raw, pct):.0f}")

# =============================================================================
# STEP 6: RAW HPV16 COMPARISONS ACROSS SUBGROUPS
# =============================================================================
log_sep("STEP 6: Raw HPV16 counts across subgroups")

# By SBS2 group
if 'SBS2_group' in master.columns:
    log("\n  --- Raw HPV16 by SBS2 group ---")
    for grp in ['HIGH', 'LOW']:
        mask = master['SBS2_group'] == grp
        vals = master.loc[mask, 'raw_HPV16']
        pos = (vals > 0).sum()
        log(f"  {grp}: n={mask.sum()}, HPV16+={pos} ({100*pos/mask.sum():.1f}%), "
            f"mean(all)={vals.mean():.2f}, mean(+)={vals[vals>0].mean():.2f}" if pos > 0 else
            f"  {grp}: n={mask.sum()}, HPV16+={pos}")
    
    # Mann-Whitney on raw counts
    h_vals = master.loc[master['SBS2_group'] == 'HIGH', 'raw_HPV16'].dropna()
    l_vals = master.loc[master['SBS2_group'] == 'LOW', 'raw_HPV16'].dropna()
    if len(h_vals) > 5 and len(l_vals) > 5:
        u, p = mannwhitneyu(h_vals, l_vals, alternative='two-sided')
        log(f"  Mann-Whitney U (raw HPV16, HIGH vs LOW): U={u:.0f}, p={p:.2e}")

# By cancer status
log("\n  --- Raw HPV16 by cancer status ---")
for status in ['Cancer cell', 'Normal cell']:
    mask = master['Final_cancer_cell_status'] == status
    vals = master.loc[mask, 'raw_HPV16']
    pos = (vals > 0).sum()
    if pos > 0:
        log(f"  {status}: n={mask.sum()}, HPV16+={pos} ({100*pos/mask.sum():.1f}%), "
            f"mean(all)={vals.mean():.2f}, mean(+)={vals[vals>0].mean():.2f}, "
            f"median(+)={vals[vals>0].median():.1f}")

# By patient
log("\n  --- Raw HPV16 by patient (basal cells) ---")
log(f"  {'Patient':20s} {'Total':>6s} {'HPV16+':>7s} {'Mean(all)':>10s} "
    f"{'Mean(+)':>8s} {'Med(+)':>7s} {'Max':>5s}")
log(f"  {'-'*20} {'-'*6} {'-'*7} {'-'*10} {'-'*8} {'-'*7} {'-'*5}")
for patient in sorted(master['subject id'].unique()):
    pmask = master['subject id'] == patient
    vals = master.loc[pmask, 'raw_HPV16']
    pos = (vals > 0).sum()
    mean_all = vals.mean()
    mean_pos = vals[vals > 0].mean() if pos > 0 else 0
    med_pos = vals[vals > 0].median() if pos > 0 else 0
    max_v = vals.max()
    log(f"  {patient:20s} {pmask.sum():6d} {pos:7d} {mean_all:10.2f} "
        f"{mean_pos:8.2f} {med_pos:7.1f} {max_v:5d}")

# =============================================================================
# STEP 7: CORRELATIONS WITH RAW COUNTS
# =============================================================================
log_sep("STEP 7: Correlations using raw HPV16 counts (HPV16+ basal cells)")

hpv_pos = master[master['HPV16_raw_positive'] == 1].copy()
log(f"  HPV16+ basal cells: {len(hpv_pos)}")

corr_targets = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'cnv_score',
                'CytoTRACE2_Score', 'sig_SBS2', 'sig_SBS13', 'n_genes']

log(f"\n  {'Variable':25s} {'rho':>8s} {'p-value':>12s}")
log(f"  {'-'*25} {'-'*8} {'-'*12}")
for var in corr_targets:
    if var in hpv_pos.columns:
        v1 = hpv_pos['raw_HPV16'].dropna()
        v2 = hpv_pos[var].dropna()
        common = v1.index.intersection(v2.index)
        if len(common) > 10:
            rho, p = spearmanr(v1.loc[common], v2.loc[common])
            log(f"  {var:25s} {rho:8.4f} {p:12.2e}")

# =============================================================================
# STEP 8: SAVE OUTPUTS
# =============================================================================
log_sep("STEP 8: Save outputs")

# Save raw counts for all cells (just HPV16+ to save space)
hpv16_df.to_csv(os.path.join(OUTPUT_DIR, "raw_HPV16_counts_all_cells.tsv"),
                sep='\t', index=False)
log(f"  Saved: raw_HPV16_counts_all_cells.tsv ({len(hpv16_df)} rows)")

# Save updated master table with raw counts
master_out = os.path.join(OUTPUT_DIR, "basal_cell_master_table_with_raw_HPV16.tsv")
master.to_csv(master_out, sep='\t')
log(f"  Saved: {master_out}")
log(f"    Shape: {master.shape}")

# Save report
report_path = os.path.join(OUTPUT_DIR, "hpv16_raw_distribution_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")

log_sep("PHASE 2 (Raw HPV16 Counts) COMPLETE")
