#!/usr/bin/env python3
"""
Phase1_5_Index_Diagnostic.py
==============================
Quick diagnostic to:
1. Copy adata_v_pp.h5ad to FIG_6 input directory
2. Compare barcode formats between signature_weights_per_cell.txt and adata_final.h5ad
3. Peek at viral adata contents (organisms detected, cell counts, HPV16 availability)

Run: python scripts/HPV_ANALYSIS/Phase1_5_Index_Diagnostic.py
"""

import os
import re
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# =============================================================================
# PATHS — UPDATE VIRAL ADATA SOURCE PATH AS NEEDED
# =============================================================================
# Source viral adata (update this to wherever you found it)
VIRAL_ADATA_SOURCE = "data/FIG_4/00_input/adata_v_pp.h5ad"
# ^^^ UPDATE THIS PATH if it's somewhere else ^^^

ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
SIG_WEIGHTS_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")

# Destination for viral adata
FIG6_INPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/00_input")
VIRAL_ADATA_DEST = os.path.join(FIG6_INPUT_DIR, "adata_v_pp.h5ad")

# =============================================================================
# STEP 1: Copy viral adata
# =============================================================================
print("=" * 70)
print("STEP 1: Copy viral adata to FIG_6 inputs")
print("=" * 70)

os.makedirs(FIG6_INPUT_DIR, exist_ok=True)

if os.path.exists(VIRAL_ADATA_SOURCE):
    if not os.path.exists(VIRAL_ADATA_DEST):
        print(f"  Copying: {VIRAL_ADATA_SOURCE}")
        print(f"      To:  {VIRAL_ADATA_DEST}")
        shutil.copy2(VIRAL_ADATA_SOURCE, VIRAL_ADATA_DEST)
        print(f"  Done. Size: {os.path.getsize(VIRAL_ADATA_DEST) / 1e9:.2f} GB")
    else:
        print(f"  Already exists: {VIRAL_ADATA_DEST}")
        print(f"  Size: {os.path.getsize(VIRAL_ADATA_DEST) / 1e9:.2f} GB")
else:
    print(f"  WARNING: Source not found at {VIRAL_ADATA_SOURCE}")
    print(f"  Please update VIRAL_ADATA_SOURCE in this script.")
    print(f"  Common locations to check:")
    print(f"    /master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/")
    print(f"    /master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_*/viral/")
    print(f"    Look for: adata_v_pp.h5ad, adata_with_virus.h5ad, viral_counts.h5ad")

# =============================================================================
# STEP 2: Compare barcode index formats
# =============================================================================
print("\n" + "=" * 70)
print("STEP 2: Barcode index format comparison")
print("=" * 70)

# --- adata_final barcodes ---
print("\n--- adata_final.h5ad ---")
adata = sc.read_h5ad(ADATA_FINAL_PATH, backed='r')
adata_barcodes = list(adata.obs_names)
print(f"  Total cells: {len(adata_barcodes)}")
print(f"  First 5 barcodes:")
for bc in adata_barcodes[:5]:
    print(f"    '{bc}'")
print(f"  Last 5 barcodes:")
for bc in adata_barcodes[-5:]:
    print(f"    '{bc}'")

# Analyze barcode structure
sample_bc = adata_barcodes[0]
print(f"\n  Structure analysis of first barcode: '{sample_bc}'")
print(f"    Length: {len(sample_bc)}")
print(f"    Contains '-': {'-' in sample_bc}")
print(f"    Parts split by '-': {sample_bc.split('-')}")

# Check for unique prefixes/suffixes
prefixes = set()
suffixes = set()
for bc in adata_barcodes[:1000]:
    parts = bc.split('-')
    if len(parts) >= 2:
        prefixes.add(parts[0][:8])  # first 8 chars of sequence
        suffixes.add('-'.join(parts[1:]))  # everything after first dash
print(f"  Unique suffix patterns (first 1000): {len(suffixes)}")
if len(suffixes) <= 50:
    print(f"  Suffixes: {sorted(suffixes)[:20]}")

adata.file.close()

# --- signature_weights barcodes ---
print("\n--- signature_weights_per_cell.txt ---")
if os.path.exists(SIG_WEIGHTS_PATH):
    # Just read the header and first few lines
    sig_df = pd.read_csv(SIG_WEIGHTS_PATH, sep='\t', index_col=0, nrows=5)
    print(f"  Columns ({len(sig_df.columns)}): {list(sig_df.columns)[:10]}...")
    print(f"  First 5 index values:")
    for idx in sig_df.index:
        print(f"    '{idx}'")
    
    # Read full index only (skip data)
    sig_full_idx = pd.read_csv(SIG_WEIGHTS_PATH, sep='\t', usecols=[0]).iloc[:, 0]
    # Actually re-read properly
    sig_full = pd.read_csv(SIG_WEIGHTS_PATH, sep='\t', index_col=0)
    sig_barcodes = list(sig_full.index)
    print(f"\n  Total cells in sig weights: {len(sig_barcodes)}")
    print(f"  Total signatures: {len(sig_full.columns)}")
    
    # Structure analysis
    sample_sig = sig_barcodes[0]
    print(f"\n  Structure analysis of first barcode: '{sample_sig}'")
    print(f"    Length: {len(sample_sig)}")
    print(f"    Contains '-': {'-' in sample_sig}")
    print(f"    Parts split by '-': {sample_sig.split('-')}")
    
    # --- Direct overlap ---
    overlap_direct = set(adata_barcodes) & set(sig_barcodes)
    print(f"\n  Direct index overlap: {len(overlap_direct)}")
    
    if len(overlap_direct) == 0:
        print("\n  Attempting barcode matching strategies...")
        
        # Strategy A: sig barcodes might be missing the SRR suffix
        # adata format: SEQUENCE-1-SRR#####
        # sig format might be: SEQUENCE-1  or  SEQUENCE  or  SRR#####_SEQUENCE-1
        
        # Check if adata barcodes contain sig barcodes as substrings
        adata_set_str = adata_barcodes[0]
        sig_set_str = sig_barcodes[0]
        print(f"    adata example: '{adata_set_str}'")
        print(f"    sig example:   '{sig_set_str}'")
        
        # Strategy A: Strip SRR suffix from adata
        adata_stripped = {}
        for bc in adata_barcodes[:100]:
            parts = bc.split('-')
            if len(parts) >= 3:
                stripped = '-'.join(parts[:2])  # keep SEQUENCE-1
                adata_stripped[stripped] = adata_stripped.get(stripped, 0) + 1
        
        overlap_stripped = set(adata_stripped.keys()) & set(sig_barcodes[:100])
        print(f"\n    Strategy A (strip SRR from adata, first 100): {len(overlap_stripped)} matches")
        
        # Strategy B: Strip suffix from sig barcodes
        sig_stripped = {}
        for bc in sig_barcodes[:100]:
            parts = bc.split('-')
            if len(parts) >= 2:
                stripped = parts[0]
                sig_stripped[stripped] = bc
        
        # Strategy C: Check if sig barcodes have a different prefix
        # Try matching the 16-char nucleotide sequence portion
        import re
        adata_seqs = {}
        for bc in adata_barcodes[:100]:
            seq_match = re.match(r'([ACGT]{16})', bc)
            if seq_match:
                adata_seqs[seq_match.group(1)] = bc
        
        sig_seqs = {}
        for bc in sig_barcodes[:100]:
            seq_match = re.search(r'([ACGT]{16})', bc)
            if seq_match:
                sig_seqs[seq_match.group(1)] = bc
        
        overlap_seq = set(adata_seqs.keys()) & set(sig_seqs.keys())
        print(f"    Strategy C (16bp sequence match, first 100): {len(overlap_seq)} matches")
        if len(overlap_seq) > 0:
            # Show example mapping
            example_seq = list(overlap_seq)[0]
            print(f"      Example: adata='{adata_seqs[example_seq]}' <-> sig='{sig_seqs[example_seq]}'")
else:
    print(f"  File not found: {SIG_WEIGHTS_PATH}")

# =============================================================================
# STEP 3: Peek at viral adata
# =============================================================================
print("\n" + "=" * 70)
print("STEP 3: Viral adata contents")
print("=" * 70)

viral_path = VIRAL_ADATA_DEST if os.path.exists(VIRAL_ADATA_DEST) else VIRAL_ADATA_SOURCE

if os.path.exists(viral_path):
    adata_v = sc.read_h5ad(viral_path)
    print(f"  Total cells: {adata_v.n_obs}")
    print(f"  Total features (organisms + genes): {adata_v.n_vars}")
    print(f"  obs columns: {list(adata_v.obs.columns)[:20]}")
    print(f"  var columns: {list(adata_v.var.columns)}")
    
    # First 5 barcodes
    print(f"\n  First 5 barcodes:")
    for bc in adata_v.obs_names[:5]:
        print(f"    '{bc}'")
    
    # Barcode overlap with adata_final
    adata_final_barcodes = set(pd.read_csv(
        ADATA_FINAL_PATH.replace('.h5ad', '_barcodes_would_need_loading.txt'), sep='\t'
    ).index) if False else set()  # placeholder
    
    # Reload adata_final obs_names for overlap check
    adata_reload = sc.read_h5ad(ADATA_FINAL_PATH, backed='r')
    adata_final_bcs = set(adata_reload.obs_names)
    adata_reload.file.close()
    
    viral_bcs = set(adata_v.obs_names)
    overlap_v = adata_final_bcs & viral_bcs
    print(f"\n  Barcode overlap with adata_final: {len(overlap_v)} / {len(viral_bcs)} viral cells")
    
    # Search for HPV16 in var_names
    print(f"\n  Searching for HPV/papillomavirus in var_names...")
    hpv_matches = [v for v in adata_v.var_names if 'papilloma' in v.lower() or 'hpv' in v.lower()]
    print(f"  HPV-related features: {hpv_matches}")
    
    # If found, get basic stats
    for hpv_feat in hpv_matches:
        vals = adata_v[:, hpv_feat].X
        if hasattr(vals, 'toarray'):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()
        n_pos = (vals > 0).sum()
        print(f"\n  '{hpv_feat}':")
        print(f"    Cells with reads: {n_pos} / {adata_v.n_obs} ({100*n_pos/adata_v.n_obs:.1f}%)")
        if n_pos > 0:
            pos_vals = vals[vals > 0]
            print(f"    Among positive cells:")
            print(f"      Mean:   {pos_vals.mean():.2f}")
            print(f"      Median: {np.median(pos_vals):.2f}")
            print(f"      Max:    {pos_vals.max():.2f}")
            print(f"      Std:    {pos_vals.std():.2f}")
    
    # Also show all organisms/features that are NOT standard genes
    # (viral features typically have spaces or non-gene-like names)
    print(f"\n  Non-standard feature names (likely viral/organism):")
    non_gene = [v for v in adata_v.var_names if ' ' in v or not v.replace('-', '').replace('.', '').replace('_', '').isalnum()]
    if len(non_gene) > 0:
        for feat in non_gene[:30]:
            vals = adata_v[:, feat].X
            if hasattr(vals, 'toarray'):
                vals = vals.toarray().flatten()
            else:
                vals = np.array(vals).flatten()
            n_pos = (vals > 0).sum()
            print(f"    '{feat}': {n_pos} cells positive")
    else:
        print(f"    None found with spaces. Total var_names: {adata_v.n_vars}")
        print(f"    First 10: {list(adata_v.var_names[:10])}")
        print(f"    Last 10: {list(adata_v.var_names[-10:])}")
    
    # Check if adata_v is the combined (expression + viral) object
    # or just the viral counts
    n_standard_genes = sum(1 for v in adata_v.var_names if re.match(r'^[A-Z][A-Z0-9]+$', v) or v.startswith('ENSG'))
    print(f"\n  Standard gene-like names: {n_standard_genes} / {adata_v.n_vars}")
    if n_standard_genes > 20000:
        print(f"  --> This is the COMBINED expression + viral object")
    elif n_standard_genes < 1000:
        print(f"  --> This appears to be a VIRAL-ONLY counts object")
    else:
        print(f"  --> Mixed object, check contents")
    
    # Check if expression data is normalized (log1p)
    sample_vals = adata_v.X[:100, :100]
    if hasattr(sample_vals, 'toarray'):
        sample_vals = sample_vals.toarray()
    max_val = sample_vals.max()
    print(f"\n  Max value in first 100x100 block: {max_val:.2f}")
    if max_val < 20:
        print(f"  --> Likely log-normalized (CPM + log1p)")
    else:
        print(f"  --> Likely raw counts")

else:
    print(f"  Viral adata not found at either location.")
    print(f"  Checked: {VIRAL_ADATA_DEST}")
    print(f"  Checked: {VIRAL_ADATA_SOURCE}")

# =============================================================================
# STEP 4: Check hpv state metadata column
# =============================================================================
print("\n" + "=" * 70)
print("STEP 4: Patient-level HPV state metadata")
print("=" * 70)

adata_check = sc.read_h5ad(ADATA_FINAL_PATH, backed='r')
if 'hpv state' in adata_check.obs.columns:
    hpv_state = adata_check.obs[['hpv state', 'subject id']].drop_duplicates()
    print(f"  'hpv state' values by patient:")
    for _, row in hpv_state.iterrows():
        print(f"    {row['subject id']}: {row['hpv state']}")
    
    print(f"\n  Value counts:")
    for val, count in adata_check.obs['hpv state'].value_counts().items():
        print(f"    {val}: {count} cells")
adata_check.file.close()

print("\n" + "=" * 70)
print("DIAGNOSTIC COMPLETE")
print("=" * 70)
