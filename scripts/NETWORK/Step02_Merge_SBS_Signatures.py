#!/usr/bin/env python3
"""
Step02_Merge_SBS_Signatures.py

Merge the cleaned TCGA expression matrix with COSMIC SBS mutation
signature weights. Remove duplicate samples and identify available
cancer types.

Corresponds to original pipeline Steps 10–11.

Input:
    01_cleaned_expression/TCGA_expression_cleaned.parquet
    Mutation_Table_Tumors_TCGA.tsv  (from data/FIG_1/)

Output (-> data/FIG_2/02_merged_with_SBS/):
    TCGA_merged_expression_SBS.parquet  — merged matrix (expression + SBS weights)
    cancer_types_available.txt          — list of cancer types in the dataset
    step02_summary.txt                  — QC summary

Usage:
    python Step02_Merge_SBS_Signatures.py
"""

import os
import json
import numpy as np
import pandas as pd

from network_config import (
    MUTATION_SIGNATURE_PATH, SIGNATURE_SAMPLE_COL,
    BARCODE_PARTS_FOR_MERGE, CLINICAL_COLS,
    A3_GENES, A3_ID_TO_ALIAS,
    DIR_01_CLEANED, DIR_02_MERGED,
    banner, log, ensure_dir
)


# =============================================================================
# LOAD cleaned expression from Step 01
# =============================================================================
banner("[STEP 10] Load cleaned TCGA expression")

expr_path = os.path.join(DIR_01_CLEANED, "TCGA_expression_cleaned.parquet")
log(f"[STEP 10] Reading: {expr_path}")
tcga_df = pd.read_parquet(expr_path)
log(f"[STEP 10] Shape: {tcga_df.shape}")


# =============================================================================
# LOAD SBS signature file
# =============================================================================
banner("[STEP 10] Load SBS mutation signature file")

log(f"[STEP 10] Reading: {MUTATION_SIGNATURE_PATH}")

# Detect format (TSV vs CSV) based on file extension
if MUTATION_SIGNATURE_PATH.endswith(".tsv"):
    sig_df = pd.read_csv(MUTATION_SIGNATURE_PATH, sep="\t")
elif MUTATION_SIGNATURE_PATH.endswith(".csv"):
    sig_df = pd.read_csv(MUTATION_SIGNATURE_PATH)
else:
    # Try TSV first, fall back to CSV
    try:
        sig_df = pd.read_csv(MUTATION_SIGNATURE_PATH, sep="\t")
    except Exception:
        sig_df = pd.read_csv(MUTATION_SIGNATURE_PATH)

log(f"[STEP 10] Signature file shape: {sig_df.shape}")
log(f"[STEP 10] Signature file columns (first 10): {sig_df.columns[:10].tolist()}")

# Verify sample ID column exists
if SIGNATURE_SAMPLE_COL not in sig_df.columns:
    # Try common alternatives
    alt_cols = ["Sample Names", "Sample_Names", "sample_name", "Entity_ID"]
    found = [c for c in alt_cols if c in sig_df.columns]
    if found:
        log(f"[STEP 10] WARNING: '{SIGNATURE_SAMPLE_COL}' not found. Using '{found[0]}' instead.")
        SIGNATURE_SAMPLE_COL_USED = found[0]
    else:
        raise ValueError(
            f"Cannot find sample ID column. Expected '{SIGNATURE_SAMPLE_COL}'. "
            f"Available columns: {sig_df.columns.tolist()}"
        )
else:
    SIGNATURE_SAMPLE_COL_USED = SIGNATURE_SAMPLE_COL


# =============================================================================
# MERGE on shortened TCGA barcode
# =============================================================================
banner("[STEP 10] Merge expression + signatures on shortened barcode")

def shorten_barcode(x, n_parts=BARCODE_PARTS_FOR_MERGE):
    """TCGA-XX-YYYY-ZZZZ-... -> TCGA-XX-YYYY-ZZZZ (first n_parts)"""
    return "-".join(str(x).split("-")[:n_parts])

tcga_df["Entity_ID_short"] = tcga_df["Entity_ID"].apply(shorten_barcode)
sig_df["Entity_ID_short"]  = sig_df[SIGNATURE_SAMPLE_COL_USED].apply(shorten_barcode)

log(f"[STEP 10] Expression unique short IDs: {tcga_df['Entity_ID_short'].nunique()}")
log(f"[STEP 10] Signature unique short IDs:  {sig_df['Entity_ID_short'].nunique()}")

# Inner merge
merged = pd.merge(tcga_df, sig_df, on="Entity_ID_short", how="inner")
log(f"[STEP 10] Merged shape (before dedup): {merged.shape}")


# =============================================================================
# REMOVE duplicated samples (drop BOTH copies)
# =============================================================================
banner("[STEP 10B] Remove duplicated samples")

dup_mask = merged["Entity_ID_short"].duplicated(keep=False)
n_dup_rows = dup_mask.sum()
n_dup_ids  = merged.loc[dup_mask, "Entity_ID_short"].nunique()

log(f"[STEP 10B] Duplicated rows: {n_dup_rows} ({n_dup_ids} unique short IDs)")

if n_dup_rows > 0:
    log("[STEP 10B] Dropping ALL rows with duplicated Entity_ID_short (conservative)")
    merged = merged[~dup_mask].copy()

log(f"[STEP 10B] Shape after dedup: {merged.shape}")


# =============================================================================
# ADD A3 alias columns (convenience copies for scoring/plotting)
# =============================================================================
banner("[STEP 10A] Create A3 alias columns")

for ensg_id, alias in A3_ID_TO_ALIAS.items():
    if ensg_id in merged.columns:
        merged[alias] = pd.to_numeric(merged[ensg_id], errors="coerce")
        log(f"[STEP 10A] Created alias column: {alias} <- {ensg_id}")
    else:
        log(f"[STEP 10A] WARNING: {ensg_id} ({alias}) not found in merged data")


# =============================================================================
# IDENTIFY available cancer types
# =============================================================================
banner("[STEP 11] Identify cancer types")

cancer_types_all = sorted(merged["Project_ID"].dropna().unique().tolist())
log(f"[STEP 11] Cancer types found: {len(cancer_types_all)}")
for ct in cancer_types_all:
    n = (merged["Project_ID"] == ct).sum()
    log(f"  {ct}: {n} samples")


# =============================================================================
# Verify SBS2 column exists
# =============================================================================
if "SBS2" not in merged.columns:
    log("[WARNING] SBS2 column not found in merged data!")
    log(f"  Available SBS columns: {[c for c in merged.columns if c.startswith('SBS')][:20]}")
else:
    sbs2_valid = pd.to_numeric(merged["SBS2"], errors="coerce").dropna()
    log(f"\n[QC] SBS2 column: {len(sbs2_valid)} valid values")
    log(f"  Range: {sbs2_valid.min():.6f} — {sbs2_valid.max():.6f}")
    log(f"  Median: {sbs2_valid.median():.6f}")


# =============================================================================
# SAVE OUTPUTS
# =============================================================================
banner("[SAVE] Writing merged outputs")

out_dir = ensure_dir(DIR_02_MERGED)

# Merged matrix
parquet_path = os.path.join(out_dir, "TCGA_merged_expression_SBS.parquet")
merged.to_parquet(parquet_path, index=False)
log(f"[SAVE] Merged matrix -> {parquet_path}")

# Also save as TSV for inspection (large file — optional)
tsv_path = os.path.join(out_dir, "TCGA_merged_expression_SBS.tsv")
merged.to_csv(tsv_path, sep="\t", index=False)
log(f"[SAVE] Merged matrix (TSV) -> {tsv_path}")

# Cancer types list
ct_path = os.path.join(out_dir, "cancer_types_available.txt")
pd.Series(cancer_types_all).to_csv(ct_path, index=False, header=False)
log(f"[SAVE] Cancer types -> {ct_path}")

# Summary
summary_path = os.path.join(out_dir, "step02_summary.txt")
with open(summary_path, "w") as f:
    f.write("Step 02 — Merge Expression + SBS Signatures\n")
    f.write("=" * 50 + "\n")
    f.write(f"Expression input: {expr_path}\n")
    f.write(f"Signature input: {MUTATION_SIGNATURE_PATH}\n")
    f.write(f"Signature sample column: {SIGNATURE_SAMPLE_COL_USED}\n")
    f.write(f"Merged shape (after dedup): {merged.shape}\n")
    f.write(f"Duplicated rows removed: {n_dup_rows}\n")
    f.write(f"Cancer types: {len(cancer_types_all)}\n")
    for ct in cancer_types_all:
        n = (merged["Project_ID"] == ct).sum()
        f.write(f"  {ct}: {n}\n")
log(f"[SAVE] Summary -> {summary_path}")

banner("[STEP 02 COMPLETE]")
print(f"\n  Output directory: {out_dir}")
print(f"  Merged samples: {len(merged)}")
print(f"  Cancer types: {len(cancer_types_all)}")
