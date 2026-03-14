#!/usr/bin/env python3
"""
Step01_Load_Clean_TCGA.py

Load the pan-cancer TCGA_master_FPKM_UQ.tsv, clean gene identifiers,
remove PAR_Y duplicates, resolve versioned ENSG IDs, filter to tumor
samples, and save a clean expression matrix + ENSG-to-symbol mapping.

Corresponds to original pipeline Steps 01–09.

Input:
    TCGA_master_FPKM_UQ.tsv  (from Figure 1 data download)

Output (-> data/FIG_2/01_cleaned_expression/):
    TCGA_expression_cleaned.pkl       — samples × ENSG expression matrix + clinical cols
    ensg_to_symbol.json               — {ENSG_ID: gene_symbol} mapping
    step01_summary.txt                — QC summary

Usage:
    python Step01_Load_Clean_TCGA.py
"""

import os
import re
import json
import numpy as np
import pandas as pd

# ---- shared config
from network_config import (
    TCGA_EXPRESSION_PATH, N_CLINICAL, CLINICAL_COLS,
    DIR_01_CLEANED, banner, log, ensure_dir
)


# =============================================================================
# STEP 01 — Load raw TSV
# =============================================================================
banner("[STEP 01] Load TCGA gene expression data (raw)")

log(f"[STEP 01] Reading: {TCGA_EXPRESSION_PATH}")
df_raw = pd.read_csv(TCGA_EXPRESSION_PATH, sep="\t", header=None, low_memory=False)
log(f"[STEP 01] Raw shape: {df_raw.shape}")
log("[STEP 01] Rows 0-2 = gene annotations (ENSG IDs, symbols, biotypes)")
log("[STEP 01] Rows 3+  = sample expression data")
log(f"[STEP 01] Columns 0-4 = clinical metadata ({CLINICAL_COLS})")


# =============================================================================
# STEP 02 — Remove _PAR_Y gene columns
# =============================================================================
banner("[STEP 02] Remove _PAR_Y gene columns")

ensg_row0 = df_raw.iloc[0, N_CLINICAL:].astype(str)
par_y_mask = ensg_row0.str.contains("_PAR_Y", na=False)
cols_to_drop = ensg_row0[par_y_mask].index.tolist()

log(f"[STEP 02] Found {len(cols_to_drop)} _PAR_Y columns")

if len(cols_to_drop) > 0:
    # Sanity check: should be all zeros in expression rows
    expr_block = df_raw.iloc[3:, N_CLINICAL:].apply(pd.to_numeric, errors="coerce")
    par_y_sums = expr_block.loc[:, cols_to_drop].sum()
    if (par_y_sums != 0).any():
        log("[STEP 02] WARNING: Some _PAR_Y columns have non-zero sums — review before dropping")
    else:
        log("[STEP 02] OK: All _PAR_Y columns sum to zero — safe to drop")

df = df_raw.drop(columns=cols_to_drop).copy()
log(f"[STEP 02] Shape after dropping _PAR_Y: {df.shape}")


# =============================================================================
# STEP 03 — Assign column headers
# =============================================================================
banner("[STEP 03] Assign column headers from row 0")

clinical_headers = df.iloc[0, :N_CLINICAL].astype(str).tolist()
gene_headers_raw = df.iloc[0, N_CLINICAL:].astype(str).tolist()
df.columns = clinical_headers + gene_headers_raw

log(f"[STEP 03] Total columns: {len(df.columns)}")
log(f"[STEP 03] First 10: {df.columns[:10].tolist()}")


# =============================================================================
# STEP 04 — Extract gene annotations (rows 0, 1, 2)
# =============================================================================
banner("[STEP 04] Extract gene annotations")

ensg_ids_raw  = df.iloc[0, N_CLINICAL:].astype(str)    # versioned ENSG IDs
gene_symbols  = df.iloc[1, N_CLINICAL:].astype(str)    # gene symbols
biotype_flags = df.iloc[2, N_CLINICAL:].astype(str)    # biotype info

log(f"[STEP 04] Example ENSG IDs:   {ensg_ids_raw.iloc[:5].tolist()}")
log(f"[STEP 04] Example symbols:    {gene_symbols.iloc[:5].tolist()}")
log(f"[STEP 04] Example biotypes:   {biotype_flags.iloc[:5].tolist()}")


# =============================================================================
# STEP 05 — Clean ENSG IDs (remove .XX version suffix if present)
# =============================================================================
# NOTE (Jake): This is a SAFETY GUARD, not an active transformation.
# Jake's R pipeline (Master_TCGA_RNA-seq_Counts_Table.R) already produces
# unversioned ENSG IDs, so this step should be a no-op (0 duplicates).
#
# WHY KEEP IT: Downstream pathway tools (g:Profiler, clusterProfiler,
# Enrichr, Reactome) require bare ENSG IDs (no .XX suffix). If the upstream
# R pipeline ever changes or a different input file is used, this catches it.
#
# EARMARKED FOR REMOVAL: If confirmed unnecessary after full pipeline run,
# this block can be safely deleted — no downstream steps depend on the
# cleaning function itself, only on the resulting IDs being unversioned.
# =============================================================================
banner("[STEP 05] Clean ENSG IDs (safety check for version suffixes)")

def clean_ensg(x):
    """ENSG00000141510.18 -> ENSG00000141510"""
    return re.sub(r"\.\d+$", "", str(x))

# Check whether any versioned IDs actually exist before cleaning
n_versioned = ensg_ids_raw.str.contains(r"\.\d+$", regex=True, na=False).sum()
log(f"[STEP 05] ENSG IDs with version suffix (.XX): {n_versioned} / {len(ensg_ids_raw)}")

if n_versioned == 0:
    log("[STEP 05] No versioned IDs found — IDs are already clean (safety check passed)")
    ensg_ids_clean = ensg_ids_raw.copy()
else:
    log(f"[STEP 05] Stripping version suffixes from {n_versioned} IDs...")
    ensg_ids_clean = ensg_ids_raw.map(clean_ensg)

n_total = len(ensg_ids_raw)
n_unique = ensg_ids_clean.nunique()
n_dups = n_total - n_unique

log(f"[STEP 05] Total gene columns: {n_total}")
log(f"[STEP 05] Unique base ENSG IDs: {n_unique}")
log(f"[STEP 05] Duplicates from version removal: {n_dups}")
if n_dups > 0:
    log(f"[STEP 05] WARNING: {n_dups} duplicates will be merged by SUM in Step 08")


# =============================================================================
# STEP 06 — Build expression matrix (rows 3+) and metadata
# =============================================================================
banner("[STEP 06] Build expression + metadata matrices")

meta = df.iloc[3:, :N_CLINICAL].copy().reset_index(drop=True)
meta.columns = CLINICAL_COLS

expr = df.iloc[3:, N_CLINICAL:].copy().reset_index(drop=True)
expr = expr.apply(pd.to_numeric, errors="coerce")

log(f"[STEP 06] meta shape: {meta.shape}")
log(f"[STEP 06] expr shape: {expr.shape}")


# =============================================================================
# STEP 07 — Rename expression columns to cleaned ENSG IDs
# =============================================================================
banner("[STEP 07] Rename columns to cleaned ENSG IDs")

expr.columns = ensg_ids_clean.tolist()
log(f"[STEP 07] First 10 gene columns: {expr.columns[:10].tolist()}")


# =============================================================================
# STEP 08 — Merge duplicate columns (if version removal created duplicates)
# =============================================================================
banner("[STEP 08] Merge duplicate ENSG columns (SUM)")

if expr.columns.duplicated().any():
    n_dup_cols = expr.columns.duplicated().sum()
    log(f"[STEP 08] Merging {n_dup_cols} duplicate columns by SUM...")
    expr = expr.T.groupby(level=0).sum(min_count=1).T
    log(f"[STEP 08] Expression shape after merge: {expr.shape}")
else:
    log("[STEP 08] No duplicate columns — nothing to merge")


# =============================================================================
# STEP 09 — Build ENSG -> gene symbol mapping
# =============================================================================
banner("[STEP 09] Build ENSG -> GeneSymbol mapping")

ensg_to_symbol = {}
for e_clean, sym in zip(ensg_ids_clean.tolist(), gene_symbols.tolist()):
    if e_clean not in ensg_to_symbol:
        ensg_to_symbol[e_clean] = sym

log(f"[STEP 09] Mapping size: {len(ensg_to_symbol)}")
log(f"[STEP 09] Examples: {list(ensg_to_symbol.items())[:5]}")


# =============================================================================
# STEP 09A — Combine meta + expr and filter to tumor samples only
# =============================================================================
# Uses the Tissue_Type metadata column (set to "Tumor" or "Normal" by
# Jake's R pipeline in Master_TCGA_RNA-seq_Counts_Table.R) rather than
# parsing the TCGA barcode. This is more direct and reliable.
# =============================================================================
banner("[STEP 09A] Combine and filter to tumor samples")

tcga_df = pd.concat([meta, expr], axis=1)
log(f"[STEP 09A] Combined shape: {tcga_df.shape}")

# Verify Tissue_Type column exists and inspect values
if "Tissue_Type" not in tcga_df.columns:
    raise ValueError(
        "Tissue_Type column not found in metadata. "
        f"Available columns: {tcga_df.columns[:10].tolist()}"
    )

tissue_counts = tcga_df["Tissue_Type"].value_counts(dropna=False)
log(f"[STEP 09A] Tissue_Type distribution:")
for val, count in tissue_counts.items():
    log(f"  {val}: {count}")

n_tumor  = (tcga_df["Tissue_Type"] == "Tumor").sum()
n_normal = (tcga_df["Tissue_Type"] == "Normal").sum()
n_other  = len(tcga_df) - n_tumor - n_normal

log(f"[STEP 09A] Summary: Tumor={n_tumor}, Normal={n_normal}, Other={n_other}")

tcga_df = tcga_df[tcga_df["Tissue_Type"] == "Tumor"].copy()
log(f"[STEP 09A] After tumor-only filter: {tcga_df.shape}")


# =============================================================================
# SAVE OUTPUTS
# =============================================================================
banner("[SAVE] Writing cleaned outputs")

out_dir = ensure_dir(DIR_01_CLEANED)

# Save expression as pickle (fast I/O between steps, no pyarrow dependency)
pickle_path = os.path.join(out_dir, "TCGA_expression_cleaned.pkl")
tcga_df.to_pickle(pickle_path)
log(f"[SAVE] Expression matrix (pickle) -> {pickle_path}")

# Also save as TSV for inspection in spreadsheet/text editor
tsv_path = os.path.join(out_dir, "TCGA_expression_cleaned.tsv")
tcga_df.to_csv(tsv_path, sep="\t", index=False)
log(f"[SAVE] Expression matrix (TSV) -> {tsv_path}")

# Save symbol mapping
json_path = os.path.join(out_dir, "ensg_to_symbol.json")
with open(json_path, "w") as f:
    json.dump(ensg_to_symbol, f, indent=2)
log(f"[SAVE] ENSG->Symbol mapping -> {json_path}")

# Save biotype mapping (needed by Step03 for protein-coding filter)
ensg_to_biotype = {}
for e_clean, bio in zip(ensg_ids_clean.tolist(), biotype_flags.tolist()):
    if e_clean not in ensg_to_biotype:
        ensg_to_biotype[e_clean] = bio

biotype_path = os.path.join(out_dir, "ensg_to_biotype.json")
with open(biotype_path, "w") as f:
    json.dump(ensg_to_biotype, f, indent=2)
log(f"[SAVE] ENSG->Biotype mapping -> {biotype_path}")

# Report biotype distribution
from collections import Counter
bio_counts = Counter(ensg_to_biotype.values())
log(f"[SAVE] Biotype distribution (top 10):")
for bio, cnt in bio_counts.most_common(10):
    log(f"  {bio}: {cnt}")

# Summary
summary_path = os.path.join(out_dir, "step01_summary.txt")
with open(summary_path, "w") as f:
    f.write("Step 01 — Load & Clean TCGA Expression Data\n")
    f.write("=" * 50 + "\n")
    f.write(f"Input file: {TCGA_EXPRESSION_PATH}\n")
    f.write(f"Raw shape: {df_raw.shape}\n")
    f.write(f"PAR_Y columns removed: {len(cols_to_drop)}\n")
    f.write(f"ENSG IDs with version suffix: {n_versioned} (safety check)\n")
    f.write(f"Duplicate ENSG columns merged: {n_dups}\n")
    f.write(f"Tumor filter: Tissue_Type == 'Tumor'\n")
    f.write(f"Tumor samples retained: {len(tcga_df)}\n")
    f.write(f"Final gene columns: {len(tcga_df.columns) - len(CLINICAL_COLS)}\n")
    f.write(f"ENSG->Symbol mapping entries: {len(ensg_to_symbol)}\n")
log(f"[SAVE] Summary -> {summary_path}")

banner("[STEP 01 COMPLETE]")
print(f"\n  Output directory: {out_dir}")
print(f"  Tumor samples: {len(tcga_df)}")
print(f"  Gene features: {len(tcga_df.columns) - len(CLINICAL_COLS)}")
