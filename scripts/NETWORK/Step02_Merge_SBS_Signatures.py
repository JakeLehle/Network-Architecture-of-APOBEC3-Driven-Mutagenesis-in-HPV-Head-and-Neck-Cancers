#!/usr/bin/env python3
"""
Step02_Merge_SBS_Signatures.py
==============================

Merge the cleaned TCGA expression matrix with SBS mutation signature
counts from SigProfiler v3.4 (COSMIC v3.4, 86 signatures).

CHANGE FROM PREVIOUS VERSION:
  The old Step02 merged directly with Mutation_Table_Tumors_TCGA.tsv
  (unknown provenance, 65 COSMIC signatures, RNA-barcode-indexed).
  
  This version uses:
    1. SigProfiler v3.4 counts (WES-barcode-indexed, documented provenance)
    2. The old file ONLY as a crosswalk (RNA barcode -> WES barcode)
    3. DIRECT crosswalk matches only (no CASE_ID 12-char fallback)
  
  This matches the v3 Figure 1 pipeline (Step05_Revised_HNSC_A3_vs_SBS2_v3.py)
  and ensures all downstream analyses use the same SBS2 values with
  documented provenance.

Merge strategy:
    Expression table:  Entity_ID (full 28-char RNA barcode)
    Crosswalk:         TCGA_Gene_Expression_Entity_ID -> Mutation_Signature__File_Orginal_Entity_ID
    SigProfiler:       WES barcode (full 28-char)

    RNA barcode --[crosswalk]--> WES barcode --[lookup]--> SBS counts

Corresponds to original pipeline Steps 10-11.

Input:
    01_cleaned_expression/TCGA_expression_cleaned.pkl
    Mutation_Table_Tumors_TCGA.tsv  (crosswalk only: RNA -> WES mapping)
    TCGA_SBS_signature_counts.tsv   (SigProfiler v3.4 output)

Output (-> data/FIG_2/02_merged_with_SBS/):
    TCGA_merged_expression_SBS.pkl  -- merged matrix (expression + SBS counts)
    TCGA_merged_expression_SBS.tsv  -- same, for inspection
    cancer_types_available.txt      -- list of cancer types in the dataset
    step02_summary.txt              -- QC summary

Usage:
    python Step02_Merge_SBS_Signatures.py
"""

import os
import numpy as np
import pandas as pd

from network_config import (
    CLINICAL_COLS, A3_GENES, A3_ID_TO_ALIAS,
    DIR_01_CLEANED, DIR_02_MERGED,
    banner, log, ensure_dir
)

# =============================================================================
# NEW INPUT PATHS
# =============================================================================
# These should be added to network_config.py. Defined here for clarity
# until the config is updated.

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

# SigProfiler v3.4 counts (WES-barcode-indexed, COSMIC v3.4, 86 signatures)
NEW_COUNTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output",
                                "TCGA_SBS_signature_counts.tsv")

# Old file used ONLY as crosswalk (RNA barcode -> WES barcode)
CROSSWALK_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1",
                               "Mutation_Table_Tumors_TCGA.tsv")
CROSSWALK_RNA_COL = "TCGA_Gene_Expression_Entity_ID"
CROSSWALK_WES_COL = "Mutation_Signature__File_Orginal_Entity_ID"


# =============================================================================
# STEP 10 -- Load cleaned expression from Step 01
# =============================================================================
banner("[STEP 10] Load cleaned TCGA expression")

expr_path = os.path.join(DIR_01_CLEANED, "TCGA_expression_cleaned.pkl")
log(f"[STEP 10] Reading: {expr_path}")
tcga_df = pd.read_pickle(expr_path)
log(f"[STEP 10] Shape: {tcga_df.shape}")
log(f"[STEP 10] Entity_ID examples: {tcga_df['Entity_ID'].head(3).tolist()}")
log(f"[STEP 10] Entity_ID barcode length: {tcga_df['Entity_ID'].str.len().median():.0f} chars")


# =============================================================================
# STEP 10A -- Load SigProfiler v3.4 counts (new provenance-tracked source)
# =============================================================================
banner("[STEP 10A] Load SigProfiler v3.4 SBS counts")

log(f"[STEP 10A] Reading: {NEW_COUNTS_PATH}")
raw_ct = pd.read_csv(NEW_COUNTS_PATH, sep="\t")
log(f"[STEP 10A] Raw shape: {raw_ct.shape}")
log(f"[STEP 10A] First column: '{raw_ct.columns[0]}'")
log(f"[STEP 10A] First 5 values of col 0: {raw_ct.iloc[:5, 0].tolist()}")

# The SigProfiler output is TRANSPOSED: signatures as rows, samples as columns.
# Detect this by checking if the first column values look like SBS signature names.
first_col = raw_ct.columns[0]
first_col_vals = raw_ct[first_col].astype(str)

if any(first_col_vals.str.startswith("SBS")):
    log("[STEP 10A] Detected TRANSPOSED format (sigs as rows, samples as cols)")
    log("[STEP 10A] Transposing to samples-as-rows format...")
    raw_ct = raw_ct.set_index(first_col)
    sig_df = raw_ct.T.copy()
    sig_df.index.name = "WES_Barcode"
    sig_df = sig_df.reset_index()
    # Remove any non-sample rows (e.g., Cancer_Type annotations)
    sig_df = sig_df[~sig_df["WES_Barcode"].str.contains(
        "Cancer_Type|^$", na=True, regex=True
    )].copy()
    # Convert SBS columns to numeric
    sbs_cols = [c for c in sig_df.columns if c.startswith("SBS")]
    for c in sbs_cols:
        sig_df[c] = pd.to_numeric(sig_df[c], errors="coerce")
else:
    log("[STEP 10A] Detected standard format (samples as rows)")
    sig_df = raw_ct.rename(columns={first_col: "WES_Barcode"})
    sbs_cols = [c for c in sig_df.columns if c.startswith("SBS")]

log(f"[STEP 10A] After transpose: {sig_df.shape}")
log(f"[STEP 10A] SBS signature columns: {len(sbs_cols)}")
log(f"[STEP 10A] SBS columns (first 10): {sbs_cols[:10]}")
log(f"[STEP 10A] WES_Barcode examples: {sig_df['WES_Barcode'].head(3).tolist()}")
log(f"[STEP 10A] WES_Barcode length: {sig_df['WES_Barcode'].str.len().median():.0f} chars")

if "SBS2" not in sbs_cols:
    raise ValueError("SBS2 column not found in SigProfiler output!")

sbs2_valid = pd.to_numeric(sig_df["SBS2"], errors="coerce").dropna()
log(f"[STEP 10A] SBS2 values: {len(sbs2_valid)} valid")
log(f"[STEP 10A] SBS2 range: {sbs2_valid.min():.0f} -- {sbs2_valid.max():.0f}")
log(f"[STEP 10A] SBS2 median: {sbs2_valid.median():.0f}")
log(f"[STEP 10A] SBS2 > 0: {(sbs2_valid > 0).sum()} ({100*(sbs2_valid > 0).mean():.1f}%)")


# =============================================================================
# STEP 10B -- Load crosswalk (RNA -> WES barcode mapping)
# =============================================================================
banner("[STEP 10B] Load crosswalk (RNA -> WES barcode mapping)")

log(f"[STEP 10B] Reading: {CROSSWALK_PATH}")
log(f"[STEP 10B] RNA column: {CROSSWALK_RNA_COL}")
log(f"[STEP 10B] WES column: {CROSSWALK_WES_COL}")

cw_df = pd.read_csv(CROSSWALK_PATH, sep="\t", usecols=[CROSSWALK_RNA_COL, CROSSWALK_WES_COL])
log(f"[STEP 10B] Crosswalk shape: {cw_df.shape}")

# Build RNA -> WES dictionary (DIRECT crosswalk only)
cw_dict = dict(zip(cw_df[CROSSWALK_RNA_COL].astype(str),
                    cw_df[CROSSWALK_WES_COL].astype(str)))
log(f"[STEP 10B] Crosswalk entries: {len(cw_dict)}")
log(f"[STEP 10B] RNA barcode example: {list(cw_dict.keys())[:2]}")
log(f"[STEP 10B] WES barcode example: {list(cw_dict.values())[:2]}")


# =============================================================================
# STEP 10C -- Crosswalk merge: RNA -> WES -> SBS counts
# =============================================================================
# For each expression sample (RNA barcode), look up the WES barcode
# via the crosswalk, then pull the SBS signature counts from the
# SigProfiler output using the WES barcode. Only DIRECT crosswalk
# matches are kept (no CASE_ID 12-char fallback).
# =============================================================================
banner("[STEP 10C] Crosswalk merge: RNA -> WES -> SBS counts")

# Build WES -> SBS lookup from the SigProfiler table
wes_set = set(sig_df["WES_Barcode"].values)
sig_indexed = sig_df.set_index("WES_Barcode")
log(f"[STEP 10C] WES barcodes in SigProfiler: {len(wes_set)}")

# Match expression Entity_IDs to WES barcodes via crosswalk
expr_ids = tcga_df["Entity_ID"].values
n_total = len(expr_ids)
n_in_crosswalk = 0
n_wes_found = 0
n_matched = 0

# Build mapping: RNA barcode -> row index in sig_df for matched samples
rna_to_wes = {}
for rna_bc in expr_ids:
    rna_str = str(rna_bc)
    if rna_str in cw_dict:
        n_in_crosswalk += 1
        wes_bc = cw_dict[rna_str]
        if wes_bc in wes_set:
            n_wes_found += 1
            rna_to_wes[rna_str] = wes_bc
            n_matched += 1

log(f"[STEP 10C] Expression samples total: {n_total}")
log(f"[STEP 10C] Found in crosswalk: {n_in_crosswalk}")
log(f"[STEP 10C] WES barcode in SigProfiler: {n_wes_found}")
log(f"[STEP 10C] Successfully matched: {n_matched}")

if n_matched == 0:
    raise ValueError(
        "Zero matches! Check barcode formats.\n"
        f"  Expression Entity_ID example: {expr_ids[0]}\n"
        f"  Crosswalk RNA example: {list(cw_dict.keys())[0]}\n"
        f"  Crosswalk WES example: {list(cw_dict.values())[0]}\n"
        f"  SigProfiler WES example: {list(wes_set)[:1]}"
    )

# Build the SBS columns for matched samples
# Strategy: create a DataFrame with RNA barcode as index holding all SBS values,
# then merge with the expression DataFrame.
matched_rows = []
for rna_bc, wes_bc in rna_to_wes.items():
    sbs_row = sig_indexed.loc[wes_bc, sbs_cols].to_dict()
    sbs_row["Entity_ID"] = rna_bc
    sbs_row["WES_Barcode"] = wes_bc
    sbs_row["match_source"] = "DIRECT"
    matched_rows.append(sbs_row)

sbs_matched = pd.DataFrame(matched_rows)
log(f"[STEP 10C] SBS matched table shape: {sbs_matched.shape}")

# Inner merge with expression on Entity_ID
merged = pd.merge(tcga_df, sbs_matched, on="Entity_ID", how="inner")
log(f"[STEP 10C] Merged shape: {merged.shape}")


# =============================================================================
# STEP 10D -- Duplicate safety check
# =============================================================================
# With DIRECT crosswalk matching on full 28-char barcodes, duplicates
# should not occur. If they do, the source files have repeated entries.
# =============================================================================
banner("[STEP 10D] Duplicate safety check")

dup_mask = merged["Entity_ID"].duplicated(keep=False)
n_dup_rows = int(dup_mask.sum())
n_dup_ids = int(merged.loc[dup_mask, "Entity_ID"].nunique()) if n_dup_rows > 0 else 0

if n_dup_rows == 0:
    log("[STEP 10D] No duplicates found (expected with DIRECT crosswalk merge)")
else:
    log(f"[STEP 10D] WARNING: {n_dup_rows} duplicated rows ({n_dup_ids} unique IDs)")
    log("[STEP 10D] Dropping ALL rows with duplicated Entity_ID (conservative)")
    merged = merged[~dup_mask].copy()
    log(f"[STEP 10D] Shape after dedup: {merged.shape}")


# =============================================================================
# STEP 10E -- Create A3 alias columns (convenience copies for scoring/plotting)
# =============================================================================
banner("[STEP 10E] Create A3 alias columns")

for ensg_id, alias in A3_ID_TO_ALIAS.items():
    if ensg_id in merged.columns:
        merged[alias] = pd.to_numeric(merged[ensg_id], errors="coerce")
        log(f"[STEP 10E] Created alias column: {alias} <- {ensg_id}")
    else:
        log(f"[STEP 10E] WARNING: {ensg_id} ({alias}) not found in merged data")


# =============================================================================
# STEP 11 -- Identify available cancer types
# =============================================================================
banner("[STEP 11] Identify cancer types")

cancer_types_all = sorted(merged["Project_ID"].dropna().unique().tolist())
log(f"[STEP 11] Cancer types found: {len(cancer_types_all)}")
for ct in cancer_types_all:
    n = (merged["Project_ID"] == ct).sum()
    log(f"  {ct}: {n} samples")


# =============================================================================
# QC -- Verify SBS2 column in merged data
# =============================================================================
banner("[QC] Verify SBS2 in merged data")

if "SBS2" not in merged.columns:
    log("[QC] CRITICAL: SBS2 column not found in merged data!")
    log(f"[QC] Available SBS columns: {[c for c in merged.columns if c.startswith('SBS')][:20]}")
else:
    sbs2_merged = pd.to_numeric(merged["SBS2"], errors="coerce").dropna()
    log(f"[QC] SBS2 column: {len(sbs2_merged)} valid values")
    log(f"[QC] Range: {sbs2_merged.min():.0f} -- {sbs2_merged.max():.0f}")
    log(f"[QC] Median: {sbs2_merged.median():.0f}")
    log(f"[QC] Mean: {sbs2_merged.mean():.1f}")
    log(f"[QC] SBS2 > 0: {(sbs2_merged > 0).sum()} ({100*(sbs2_merged > 0).mean():.1f}%)")

    # HNSC-specific QC
    hnsc_mask = merged["Project_ID"] == "TCGA-HNSC"
    if hnsc_mask.sum() > 0:
        hnsc_sbs2 = pd.to_numeric(merged.loc[hnsc_mask, "SBS2"], errors="coerce").dropna()
        log(f"\n[QC] HNSC-specific:")
        log(f"  HNSC tumors: {hnsc_mask.sum()}")
        log(f"  SBS2 range: {hnsc_sbs2.min():.0f} -- {hnsc_sbs2.max():.0f}")
        log(f"  SBS2 median: {hnsc_sbs2.median():.0f}")
        log(f"  SBS2 > 0: {(hnsc_sbs2 > 0).sum()} ({100*(hnsc_sbs2 > 0).mean():.1f}%)")


# =============================================================================
# COMPARISON -- Report what changed vs old merge (informational only)
# =============================================================================
banner("[COMPARE] Old vs new merge summary")

log(f"[COMPARE] Previous merge: Mutation_Table_Tumors_TCGA.tsv")
log(f"  Method: direct inner join on RNA barcode (Entity_ID)")
log(f"  Source: unknown COSMIC version, 65 SBS signatures")
log(f"  Previous HNSC count: ~425 tumors")
log(f"")
log(f"[COMPARE] Current merge: SigProfiler v3.4 counts + DIRECT crosswalk")
log(f"  Method: RNA -> WES via crosswalk, then WES -> SBS from SigProfiler")
log(f"  Source: COSMIC v3.4, {len(sbs_cols)} SBS signatures")
log(f"  Current total matched: {len(merged)} samples")
hnsc_n = (merged["Project_ID"] == "TCGA-HNSC").sum() if "Project_ID" in merged.columns else "N/A"
log(f"  Current HNSC count: {hnsc_n} tumors")
log(f"  Match source: DIRECT crosswalk only (no CASE_ID fallback)")


# =============================================================================
# SAVE OUTPUTS
# =============================================================================
banner("[SAVE] Writing merged outputs")

out_dir = ensure_dir(DIR_02_MERGED)

# Pickle (fast I/O between steps)
pickle_path = os.path.join(out_dir, "TCGA_merged_expression_SBS.pkl")
merged.to_pickle(pickle_path)
log(f"[SAVE] Merged matrix (pickle) -> {pickle_path}")

# TSV for inspection
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
    f.write("Step 02 -- Merge Expression + SBS Signatures (UPDATED)\n")
    f.write("=" * 60 + "\n")
    f.write(f"Expression input: {expr_path}\n")
    f.write(f"SigProfiler input: {NEW_COUNTS_PATH}\n")
    f.write(f"Crosswalk input: {CROSSWALK_PATH}\n")
    f.write(f"Crosswalk columns: {CROSSWALK_RNA_COL} -> {CROSSWALK_WES_COL}\n")
    f.write(f"Match method: DIRECT crosswalk only (no CASE_ID fallback)\n")
    f.write(f"SBS signature version: COSMIC v3.4 ({len(sbs_cols)} signatures)\n")
    f.write(f"\n")
    f.write(f"Expression samples: {n_total}\n")
    f.write(f"In crosswalk: {n_in_crosswalk}\n")
    f.write(f"WES in SigProfiler: {n_wes_found}\n")
    f.write(f"Successfully matched: {n_matched}\n")
    f.write(f"Merged shape: {merged.shape}\n")
    f.write(f"Duplicate rows found: {n_dup_rows}\n")
    f.write(f"\n")
    f.write(f"Cancer types: {len(cancer_types_all)}\n")
    for ct in cancer_types_all:
        n = (merged["Project_ID"] == ct).sum()
        f.write(f"  {ct}: {n}\n")
log(f"[SAVE] Summary -> {summary_path}")

banner("[STEP 02 COMPLETE]")
print(f"\n  Output directory: {out_dir}")
print(f"  Merged samples: {len(merged)}")
print(f"  Cancer types: {len(cancer_types_all)}")
print(f"  SBS signatures: {len(sbs_cols)} (COSMIC v3.4)")
