#!/usr/bin/env python3
"""
Diagnostic_Figure4_Issues.py
==============================

Quick diagnostic to investigate three issues before the final Figure 4 rerun:

  1. CELL COUNT: Why 707 vs 546? Check adata details and L-method determinism.
  2. GENE SYMBOL OVERLAP: Do the SC and TCGA gene namespaces actually align?
     Check TCGA community genes against the SC expression matrix index.
     Check Harris interactor genes similarly.
  3. A3 GENE NETWORK MEMBERSHIP: Why are A3 genes absent from communities?
     Report their max |delta-rho| values.

This reads existing outputs — no heavy computation needed.

Usage:
  conda run -n NETWORK python Diagnostic_Figure4_Issues.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import json
import pickle
import numpy as np
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIG — hardcoded paths (adjust if needed)
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

# SC data
SC_HIGH_PATH = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv")
SC_LOW_PATH = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv")
SC_GROUPS_PATH = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv")
SC_PARTITION_PATH = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_best_partition.csv")
SC_DIFFEXPR_PATH = os.path.join(BASE_DIR, "data/FIG_4/02_differential_expression/SC_selected_genes_filtered.csv")
SC_CORR_DIFF_PATH = os.path.join(BASE_DIR, "data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl")

# TCGA / Figure 2 data
TCGA_PARTITION_PATH = os.path.join(BASE_DIR, "data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv")
ENSG_TO_SYMBOL_PATH = os.path.join(BASE_DIR, "data/FIG_2/01_cleaned_expression/ensg_to_symbol.json")

# Harris interactors
HARRIS_ALL_PATH = os.path.join(BASE_DIR, "data/FIG_4/00_input/Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(BASE_DIR, "data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt")

# adata (for cell count check)
ADATA_PATH = os.path.join(BASE_DIR, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH = os.path.join(BASE_DIR, "data/FIG_4/00_input/signature_weights_per_cell.txt")

# A3 genes
A3_SYMBOLS = ["APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"]


def banner(title):
    line = "=" * 80
    print(f"\n{line}\n{title}\n{line}", flush=True)


def log(msg):
    print(msg, flush=True)


# =============================================================================
# DIAGNOSTIC 1: CELL COUNT AND ADATA DETAILS
# =============================================================================

def check_cell_counts():
    banner("DIAGNOSTIC 1: CELL COUNT AND ADATA DETAILS")

    # Check current group assignments
    if os.path.exists(SC_GROUPS_PATH):
        groups = pd.read_csv(SC_GROUPS_PATH, sep="\t")
        n_high = (groups["group"] == "HIGH").sum()
        n_low = (groups["group"] == "LOW").sum()
        log(f"Current group assignments: {n_high} HIGH / {n_low} LOW")
    else:
        log("WARNING: Group assignments file not found")

    # Check adata
    if os.path.exists(ADATA_PATH):
        log(f"\nLoading adata to check cell counts...")
        import scanpy as sc
        adata = sc.read_h5ad(ADATA_PATH)
        log(f"  Total cells in adata: {adata.n_obs:,}")
        log(f"  Total genes in adata: {adata.n_vars:,}")

        # Cell type counts
        if 'final_annotation' in adata.obs.columns:
            ct_counts = adata.obs['final_annotation'].value_counts()
            log(f"\n  Cell type counts:")
            for ct, n in ct_counts.items():
                marker = " <<<" if ct == "basal cell" else ""
                log(f"    {ct}: {n:,}{marker}")
        else:
            log("  WARNING: 'final_annotation' not in adata.obs")
            log(f"  Available obs columns: {list(adata.obs.columns)}")

        # Check SBS2 in obs
        if 'SBS2' in adata.obs.columns:
            log(f"\n  SBS2 already in adata.obs: YES")
            basal = adata[adata.obs['final_annotation'] == 'basal cell']
            sbs2 = basal.obs['SBS2']
            log(f"  Basal cells: {len(sbs2):,}")
            log(f"  SBS2 > 0: {(sbs2 > 0).sum():,}")
            log(f"  SBS2 == 0: {(sbs2 == 0).sum():,}")
        else:
            log(f"\n  SBS2 NOT in adata.obs — must be loaded from weights file")

        # Check weights file
        if os.path.exists(WEIGHTS_PATH):
            weights = pd.read_csv(WEIGHTS_PATH, sep="\t", index_col=0)
            log(f"\n  Weights file: {weights.shape[0]} signatures × {weights.shape[1]} cells")
            log(f"  Signatures: {list(weights.index)}")

            # Check overlap with adata cell barcodes
            adata_barcodes = set(adata.obs_names)
            weight_barcodes = set(weights.columns)
            overlap = adata_barcodes & weight_barcodes
            log(f"  Adata barcodes: {len(adata_barcodes):,}")
            log(f"  Weight barcodes: {len(weight_barcodes):,}")
            log(f"  Overlap: {len(overlap):,}")
            log(f"  In adata but not weights: {len(adata_barcodes - weight_barcodes):,}")

        # Gene index format check
        log(f"\n  Gene index (adata.var):")
        log(f"  var columns: {list(adata.var.columns)}")
        log(f"  var index name: {adata.var.index.name}")
        log(f"  First 10 var index values: {list(adata.var.index[:10])}")

        # Check for gene_symbol column
        for col in ['gene_symbol', 'feature_name', 'gene_name', 'symbol']:
            if col in adata.var.columns:
                log(f"  '{col}' column found — first 10: {list(adata.var[col][:10])}")

        del adata  # Free memory
    else:
        log(f"WARNING: adata not found at {ADATA_PATH}")


# =============================================================================
# DIAGNOSTIC 2: GENE SYMBOL OVERLAP
# =============================================================================

def check_gene_overlap():
    banner("DIAGNOSTIC 2: GENE SYMBOL OVERLAP (SC vs TCGA vs Harris)")

    # --- Load SC gene index ---
    log("Loading SC expression matrix gene index...")
    sc_genes_high = set(pd.read_csv(SC_HIGH_PATH, sep="\t", index_col=0, nrows=0).index)
    log(f"  SC HIGH genes: {len(sc_genes_high)}")

    sc_n_ensg = sum(1 for g in sc_genes_high if str(g).startswith("ENSG"))
    sc_n_symbol = len(sc_genes_high) - sc_n_ensg
    log(f"  Format: {sc_n_symbol} symbols + {sc_n_ensg} ENSG IDs")

    # Also get the selected genes (those that entered the network)
    sc_selected = set(pd.read_csv(SC_DIFFEXPR_PATH)["gene"])
    log(f"  SC selected (network-ready) genes: {len(sc_selected)}")

    # SC community genes
    if os.path.exists(SC_PARTITION_PATH):
        sc_comm = set(pd.read_csv(SC_PARTITION_PATH)["gene"])
        log(f"  SC community genes: {len(sc_comm)}")
    else:
        sc_comm = set()

    # --- Load TCGA community genes (raw ENSG) ---
    log(f"\nLoading TCGA community genes...")
    tcga_df = pd.read_csv(TCGA_PARTITION_PATH)
    gene_col = "gene" if "gene" in tcga_df.columns else tcga_df.columns[0]
    tcga_genes_ensg = set(tcga_df[gene_col])
    log(f"  TCGA community genes (ENSG): {len(tcga_genes_ensg)}")

    # --- Load ENSG → symbol mapping ---
    with open(ENSG_TO_SYMBOL_PATH) as f:
        ensg_to_symbol = json.load(f)
    log(f"  ENSG→symbol mapping: {len(ensg_to_symbol)} entries")

    # --- Convert TCGA to symbols ---
    tcga_symbols = set()
    tcga_unmapped = set()
    for g in tcga_genes_ensg:
        sym = ensg_to_symbol.get(g)
        if sym:
            tcga_symbols.add(sym)
        else:
            tcga_unmapped.add(g)
            tcga_symbols.add(g)  # Keep ENSG as fallback

    log(f"  TCGA community genes (symbols): {len(tcga_symbols)}")
    log(f"  TCGA unmapped (kept as ENSG): {len(tcga_unmapped)}")

    # --- Check TCGA symbols against SC gene index ---
    banner("TCGA community genes vs SC expression matrix")

    tcga_in_sc_all = tcga_symbols & sc_genes_high
    tcga_in_sc_selected = tcga_symbols & sc_selected
    tcga_in_sc_comm = tcga_symbols & sc_comm

    log(f"  TCGA community genes:          {len(tcga_symbols)}")
    log(f"  Present in SC expression:      {len(tcga_in_sc_all)} ({100*len(tcga_in_sc_all)/len(tcga_symbols):.1f}%)")
    log(f"  Present in SC selected (DE):   {len(tcga_in_sc_selected)} ({100*len(tcga_in_sc_selected)/len(tcga_symbols):.1f}%)")
    log(f"  Present in SC communities:     {len(tcga_in_sc_comm)} ({100*len(tcga_in_sc_comm)/len(tcga_symbols):.1f}%)")

    tcga_missing_from_sc = tcga_symbols - sc_genes_high
    log(f"\n  TCGA genes MISSING from SC expression matrix: {len(tcga_missing_from_sc)}")
    if len(tcga_missing_from_sc) <= 50:
        for g in sorted(tcga_missing_from_sc):
            log(f"    {g}")
    else:
        log(f"  (showing first 50 of {len(tcga_missing_from_sc)}):")
        for g in sorted(tcga_missing_from_sc)[:50]:
            log(f"    {g}")

    # --- Check for near-matches (case, prefix issues) ---
    banner("Near-match analysis (potential naming mismatches)")

    sc_genes_upper = {g.upper(): g for g in sc_genes_high}
    tcga_missing_upper = {g.upper(): g for g in tcga_missing_from_sc}

    near_matches = []
    for tc_upper, tc_orig in tcga_missing_upper.items():
        if tc_upper in sc_genes_upper:
            sc_orig = sc_genes_upper[tc_upper]
            near_matches.append((tc_orig, sc_orig))

    log(f"  Case-insensitive matches found: {len(near_matches)}")
    for tc, sc in near_matches[:30]:
        log(f"    TCGA: '{tc}' ↔ SC: '{sc}'")

    # --- Harris interactors ---
    banner("Harris interactors vs SC gene index")

    for label, path in [("ALL", HARRIS_ALL_PATH), ("A3B-only", HARRIS_A3B_PATH)]:
        if not os.path.exists(path):
            log(f"  {label}: file not found")
            continue

        harris_genes = set()
        with open(path) as f:
            for line in f:
                gene = line.strip().split("\t")[0].strip()
                if gene and not gene.startswith("#"):
                    harris_genes.add(gene)

        harris_in_sc_all = harris_genes & sc_genes_high
        harris_in_sc_selected = harris_genes & sc_selected
        harris_in_sc_comm = harris_genes & sc_comm

        log(f"\n  {label} ({len(harris_genes)} genes):")
        log(f"    In SC expression matrix: {len(harris_in_sc_all)} ({100*len(harris_in_sc_all)/len(harris_genes):.1f}%)")
        log(f"    In SC selected (DE):     {len(harris_in_sc_selected)} ({100*len(harris_in_sc_selected)/len(harris_genes):.1f}%)")
        log(f"    In SC communities:       {len(harris_in_sc_comm)} ({100*len(harris_in_sc_comm)/len(harris_genes):.1f}%)")

        harris_missing = harris_genes - sc_genes_high
        log(f"    Missing from SC expression: {len(harris_missing)}")
        if harris_missing:
            # Check case-insensitive
            harris_near = []
            for h in harris_missing:
                if h.upper() in sc_genes_upper:
                    harris_near.append((h, sc_genes_upper[h.upper()]))
            if harris_near:
                log(f"    Case-insensitive matches:")
                for h, s in harris_near:
                    log(f"      Harris: '{h}' ↔ SC: '{s}'")

        # Show which Harris genes ARE in SC
        if harris_in_sc_all:
            log(f"    Harris genes present in SC: {sorted(harris_in_sc_all)}")


# =============================================================================
# DIAGNOSTIC 3: A3 GENES IN DIFF NETWORK
# =============================================================================

def check_a3_in_network():
    banner("DIAGNOSTIC 3: A3 GENES IN DIFF CORRELATION MATRIX")

    if not os.path.exists(SC_CORR_DIFF_PATH):
        log("WARNING: DIFF correlation matrix not found")
        return

    log("Loading DIFF correlation matrix (this may take a moment)...")
    with open(SC_CORR_DIFF_PATH, "rb") as f:
        corr_diff = pickle.load(f)
    log(f"  Matrix shape: {corr_diff.shape}")

    genes_in_matrix = set(corr_diff.index)

    for a3 in A3_SYMBOLS:
        if a3 not in genes_in_matrix:
            log(f"\n  {a3}: NOT in correlation matrix")
            continue

        # Get this gene's row in the DIFF matrix
        row = corr_diff.loc[a3].drop(a3)  # Drop self
        abs_row = row.abs()

        max_val = abs_row.max()
        max_partner = abs_row.idxmax()
        n_above_03 = (abs_row >= 0.30).sum()
        n_above_04 = (abs_row >= 0.40).sum()
        n_above_045 = (abs_row >= 0.45).sum()
        n_above_05 = (abs_row >= 0.50).sum()

        log(f"\n  {a3}:")
        log(f"    Max |delta-rho|: {max_val:.4f} (with {max_partner})")
        log(f"    Edges at |delta-rho| >= 0.30: {n_above_03}")
        log(f"    Edges at |delta-rho| >= 0.40: {n_above_04}")
        log(f"    Edges at |delta-rho| >= 0.45: {n_above_045} {'← current threshold' if n_above_045 == 0 else ''}")
        log(f"    Edges at |delta-rho| >= 0.50: {n_above_05}")

        # Top 10 partners
        top10 = abs_row.nlargest(10)
        log(f"    Top 10 DIFF partners:")
        for partner, val in top10.items():
            signed = corr_diff.loc[a3, partner]
            direction = "+" if signed > 0 else "-"
            in_comm = ""
            if os.path.exists(SC_PARTITION_PATH):
                comm_df = pd.read_csv(SC_PARTITION_PATH)
                match = comm_df[comm_df["gene"] == partner]
                if len(match) > 0:
                    in_comm = f" [Community {int(match['community'].iloc[0])}]"
            log(f"      {partner}: {direction}{val:.4f}{in_comm}")

    # Check total A3 family pairwise correlations
    banner("A3 family pairwise DIFF correlations")
    a3_in_matrix = [a3 for a3 in A3_SYMBOLS if a3 in genes_in_matrix]
    if len(a3_in_matrix) >= 2:
        log(f"  A3 genes in matrix: {a3_in_matrix}")
        for i, a3_1 in enumerate(a3_in_matrix):
            for a3_2 in a3_in_matrix[i+1:]:
                val = corr_diff.loc[a3_1, a3_2]
                log(f"  {a3_1} ↔ {a3_2}: delta-rho = {val:+.4f}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    start = datetime.now()
    banner("FIGURE 4 DIAGNOSTIC — INVESTIGATING LOW OVERLAP")
    log(f"Start: {start}")

    check_cell_counts()
    check_gene_overlap()
    check_a3_in_network()

    elapsed = datetime.now() - start
    banner(f"DIAGNOSTIC COMPLETE | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
