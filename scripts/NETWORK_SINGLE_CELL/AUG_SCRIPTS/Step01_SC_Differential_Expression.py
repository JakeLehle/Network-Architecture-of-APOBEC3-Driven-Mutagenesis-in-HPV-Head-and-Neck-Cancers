#!/usr/bin/env python3
"""
Step01_SC_Differential_Expression.py
=====================================

Figure 4 -- Step 01: Differential expression between two cell populations
using scanpy's built-in rank_genes_groups (Wilcoxon method with BH-FDR).

Takes a comparison name as command-line argument:
  SBS2_VS_CNV     -- SBS2_HIGH vs CNV_HIGH
  SBS2_VS_NORMAL  -- SBS2_HIGH vs NORMAL
  CNV_VS_NORMAL   -- CNV_HIGH vs NORMAL

Loads the adata object directly and uses the three-group assignments from
Step00B. Data is already log-normalized so NO additional normalization is
applied.

Gene selection: BH-adjusted p-value (FDR) < 0.05. No force-keep for SC
since A3A and A3B pass FDR naturally in all three comparisons.

Input:
  data/FIG_4/00_input/adata_final.h5ad
  data/FIG_4/01_group_selection/three_group_assignments.tsv

Output (-> data/FIG_4/02_differential_expression/):
  SC_diffexpr_stats.csv          -- per-gene DEG results (from scanpy)
  SC_selected_genes.csv          -- genes passing FDR < 0.05
  SC_selected_genes_filtered.csv -- network-ready genes
  SC_volcano.png                 -- volcano plot
  SC_manhattan.png               -- Manhattan-style plot
  SC_filtering_summary.txt       -- filtering log

Usage:
  conda run -n NETWORK python Step01_SC_Differential_Expression.py SBS2_VS_CNV
  conda run -n NETWORK python Step01_SC_Differential_Expression.py SBS2_VS_NORMAL
  conda run -n NETWORK python Step01_SC_Differential_Expression.py CNV_VS_NORMAL

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime

from network_config_SC import (
    DIR_00_INPUT, DIR_01_GROUPS, DIR_02_DE,
    ADATA_FINAL_PATH,
    A3_GENES_SYMBOLS, A3_SYMBOL_TO_ALIAS, BIOMARKERS_SYMBOLS,
    MIN_CELLS_DETECTED, FDR_THRESHOLD, RAW_P_THRESHOLD, LOGFC_THRESHOLD,
    FORCE_KEEP_A3,
    banner, log, ensure_dir
)

# =============================================================================
# COMPARISON DEFINITIONS
# =============================================================================
# For each comparison: test_group is "HIGH" in the network context,
# ref_group is "LOW". scanpy reports DEGs for test vs reference.
COMPARISONS = {
    "SBS2_VS_CNV": {
        "test_group": "SBS2_HIGH",
        "ref_group": "CNV_HIGH",
        "description": "SBS2-HIGH vs CNV-HIGH (divergent fates)",
    },
    "SBS2_VS_NORMAL": {
        "test_group": "SBS2_HIGH",
        "ref_group": "NORMAL",
        "description": "SBS2-HIGH vs NORMAL (latent A3A-driven)",
    },
    "CNV_VS_NORMAL": {
        "test_group": "CNV_HIGH",
        "ref_group": "NORMAL",
        "description": "CNV-HIGH vs NORMAL (productive A3B-driven)",
    },
}

# Three-group assignments file
THREE_GROUP_PATH = os.path.join(DIR_01_GROUPS, "three_group_assignments.tsv")

# A3 genes to track (symbols)
A3_TRACK = ["APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
            "APOBEC3F", "APOBEC3G", "APOBEC3H"]


# =============================================================================
# MAIN
# =============================================================================

def main():
    start_time = datetime.now()

    # =========================================================================
    # Parse comparison from command line
    # =========================================================================
    if len(sys.argv) < 2:
        print("Usage: python Step01_SC_Differential_Expression.py <COMPARISON>")
        print("  COMPARISON: SBS2_VS_CNV | SBS2_VS_NORMAL | CNV_VS_NORMAL")
        sys.exit(1)

    comparison_name = sys.argv[1].upper()
    if comparison_name not in COMPARISONS:
        print(f"ERROR: Unknown comparison '{comparison_name}'")
        print(f"  Valid options: {', '.join(COMPARISONS.keys())}")
        sys.exit(1)

    comp = COMPARISONS[comparison_name]
    test_group = comp["test_group"]
    ref_group = comp["ref_group"]

    banner(f"[STEP 01] SC Differential Expression -- {comparison_name}")
    log(f"Start time: {start_time}")
    log(f"Comparison: {comp['description']}")
    log(f"Test group (HIGH): {test_group}")
    log(f"Reference group (LOW): {ref_group}")

    ensure_dir(DIR_02_DE)
    filter_log = []

    # =========================================================================
    # 1. Load adata
    # =========================================================================
    banner("[STEP 01.1] Load adata")

    log(f"Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  Full adata: {adata.n_obs} cells x {adata.n_vars} genes")
    filter_log.append(f"Full adata: {adata.n_obs} cells x {adata.n_vars} genes")

    # =========================================================================
    # 2. Load three-group assignments and add to adata
    # =========================================================================
    banner("[STEP 01.2] Load group assignments")

    log(f"Loading: {THREE_GROUP_PATH}")
    groups_df = pd.read_csv(THREE_GROUP_PATH, sep="\t")
    log(f"  Group assignments: {len(groups_df)} cells")
    log(f"  Groups: {groups_df['group'].value_counts().to_dict()}")

    # Map cell barcodes to groups
    group_map = dict(zip(groups_df["cell_barcode"], groups_df["group"]))

    # Add to adata.obs
    adata.obs["three_group"] = adata.obs.index.map(lambda x: group_map.get(x, "OTHER"))

    group_counts = adata.obs["three_group"].value_counts()
    log(f"  Mapped to adata.obs:")
    for g, n in group_counts.items():
        log(f"    {g}: {n} cells")

    # =========================================================================
    # 3. Subset to the two groups for this comparison
    # =========================================================================
    banner("[STEP 01.3] Subset to comparison groups")

    mask = adata.obs["three_group"].isin([test_group, ref_group])
    adata_sub = adata[mask].copy()

    n_test = (adata_sub.obs["three_group"] == test_group).sum()
    n_ref = (adata_sub.obs["three_group"] == ref_group).sum()
    log(f"  {test_group}: {n_test} cells")
    log(f"  {ref_group}: {n_ref} cells")
    log(f"  Total subset: {adata_sub.n_obs} cells x {adata_sub.n_vars} genes")
    filter_log.append(f"Subset: {n_test} {test_group} + {n_ref} {ref_group} = {adata_sub.n_obs} cells")

    # =========================================================================
    # 4. Filter low-detection genes
    # =========================================================================
    banner("[STEP 01.4] Filter low-detection genes")

    n_before = adata_sub.n_vars
    sc.pp.filter_genes(adata_sub, min_cells=MIN_CELLS_DETECTED)
    n_after = adata_sub.n_vars
    n_removed = n_before - n_after

    log(f"  Before filter: {n_before} genes")
    log(f"  Min cells detected: {MIN_CELLS_DETECTED}")
    log(f"  After filter: {n_after} genes (removed {n_removed})")
    filter_log.append(f"Gene filter (>={MIN_CELLS_DETECTED} cells): {n_before} -> {n_after}")

    # =========================================================================
    # 4b. Report A3 gene status
    # =========================================================================
    banner("[STEP 01.4b] A3 gene detection status")

    genes_in_adata = set(adata_sub.var_names)
    for a3 in A3_TRACK:
        alias = A3_SYMBOL_TO_ALIAS.get(a3, a3)
        if a3 in genes_in_adata:
            # Count detection in each group
            test_cells = adata_sub[adata_sub.obs["three_group"] == test_group]
            ref_cells = adata_sub[adata_sub.obs["three_group"] == ref_group]

            if hasattr(test_cells[:, a3].X, 'toarray'):
                n_test_det = int((test_cells[:, a3].X.toarray() > 0).sum())
                n_ref_det = int((ref_cells[:, a3].X.toarray() > 0).sum())
            else:
                n_test_det = int((test_cells[:, a3].X > 0).sum())
                n_ref_det = int((ref_cells[:, a3].X > 0).sum())

            log(f"  {alias} ({a3}): {test_group}={n_test_det} cells, "
                f"{ref_group}={n_ref_det} cells -- IN MATRIX")
        else:
            log(f"  {alias} ({a3}): NOT IN MATRIX (filtered out)")

    # =========================================================================
    # 4c. Biomarker status
    # =========================================================================
    banner("[STEP 01.4c] Biomarker gene status")

    for cell_type, markers in BIOMARKERS_SYMBOLS.items():
        in_matrix = [m for m in markers if m in genes_in_adata]
        if in_matrix:
            log(f"  {cell_type}: {', '.join(in_matrix)} -- IN MATRIX")
        else:
            log(f"  {cell_type}: not in expression matrix")

    # =========================================================================
    # 5. Run scanpy rank_genes_groups (Wilcoxon with BH-FDR)
    # =========================================================================
    banner("[STEP 01.5] Differential expression (scanpy rank_genes_groups)")

    log(f"  Method: wilcoxon")
    log(f"  Groupby: three_group")
    log(f"  Reference: {ref_group}")
    log(f"  Data is already log-normalized (no additional normalization)")

    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="three_group",
        method="wilcoxon",
        reference=ref_group,
        use_raw=False,
    )

    # Extract results for the test group
    results_df = sc.get.rank_genes_groups_df(adata_sub, group=test_group)

    log(f"  DEG results: {len(results_df)} genes")
    log(f"  Columns: {list(results_df.columns)}")

    # =========================================================================
    # 6. Build stats table in our standard format
    # =========================================================================
    banner("[STEP 01.6] Build stats table")

    # Rename columns to match our pipeline format
    # scanpy columns: names, scores, logfoldchanges, pvals, pvals_adj
    stats_df = pd.DataFrame({
        "gene": results_df["names"].values,
        "p_value": results_df["pvals"].values,
        "fdr": results_df["pvals_adj"].values,
        "log2FC": results_df["logfoldchanges"].values,
        "score": results_df["scores"].values,
    })

    # Add detection percentages
    test_cells = adata_sub[adata_sub.obs["three_group"] == test_group]
    ref_cells = adata_sub[adata_sub.obs["three_group"] == ref_group]

    pct_test = {}
    pct_ref = {}
    for gene in stats_df["gene"]:
        if gene in adata_sub.var_names:
            idx = list(adata_sub.var_names).index(gene)
            if hasattr(test_cells.X, 'toarray'):
                pct_test[gene] = float((test_cells.X[:, idx].toarray() > 0).mean())
                pct_ref[gene] = float((ref_cells.X[:, idx].toarray() > 0).mean())
            else:
                pct_test[gene] = float((test_cells.X[:, idx] > 0).mean())
                pct_ref[gene] = float((ref_cells.X[:, idx] > 0).mean())
        else:
            pct_test[gene] = 0.0
            pct_ref[gene] = 0.0

    stats_df["pct_test"] = stats_df["gene"].map(pct_test)
    stats_df["pct_ref"] = stats_df["gene"].map(pct_ref)

    # Sort by p-value
    stats_df = stats_df.sort_values("p_value").reset_index(drop=True)

    # Report counts
    n_sig_raw = (stats_df["p_value"] < RAW_P_THRESHOLD).sum()
    n_sig_fdr = (stats_df["fdr"] < FDR_THRESHOLD).sum()

    log(f"  Genes tested: {len(stats_df)}")
    log(f"  Significant (raw p < {RAW_P_THRESHOLD}): {n_sig_raw}")
    log(f"  Significant (FDR < {FDR_THRESHOLD}): {n_sig_fdr}")
    filter_log.append(f"Tested: {len(stats_df)}, raw p<0.05: {n_sig_raw}, FDR<0.05: {n_sig_fdr}")

    # FDR sanity check
    n_bad = int(np.sum(stats_df["fdr"].values < stats_df["p_value"].values * 0.99))
    if n_bad > 0:
        log(f"  WARNING: {n_bad} genes have FDR < raw p (should not happen)")
    else:
        log(f"  FDR sanity check: PASSED (all FDR >= raw p)")

    # =========================================================================
    # 7. Gene selection (FDR < 0.05)
    # =========================================================================
    banner("[STEP 01.7] Gene selection")

    # Select by FDR
    if LOGFC_THRESHOLD > 0:
        selected_mask = (stats_df["fdr"] < FDR_THRESHOLD) & (stats_df["log2FC"].abs() >= LOGFC_THRESHOLD)
    else:
        selected_mask = stats_df["fdr"] < FDR_THRESHOLD

    log(f"  FDR < {FDR_THRESHOLD}: {(stats_df['fdr'] < FDR_THRESHOLD).sum()} genes")
    log(f"  Selected (FDR + FC filter): {selected_mask.sum()} genes")

    # Report A3 status at FDR threshold
    log(f"\n  A3 gene DE status:")
    for a3 in A3_TRACK:
        alias = A3_SYMBOL_TO_ALIAS.get(a3, a3)
        match = stats_df[stats_df["gene"] == a3]
        if len(match) > 0:
            row = match.iloc[0]
            passes = "SELECTED" if row["fdr"] < FDR_THRESHOLD else "not selected"
            log(f"    {alias}: raw_p={row['p_value']:.2e}, fdr={row['fdr']:.2e}, "
                f"log2FC={row['log2FC']:.3f} -- {passes}")
        else:
            log(f"    {alias}: not in results")

    stats_df["selected"] = selected_mask
    selected_df = stats_df[stats_df["selected"]].copy()
    log(f"  Final selected genes: {len(selected_df)}")
    filter_log.append(f"Selected (FDR < {FDR_THRESHOLD}): {len(selected_df)}")

    # =========================================================================
    # 8. Network-readiness filter
    # =========================================================================
    banner("[STEP 01.8] Network-readiness filter")

    sel_genes = selected_df["gene"].tolist()
    ready_genes = []
    genes_in_var = set(adata_sub.var_names)

    for gene in sel_genes:
        if gene not in genes_in_var:
            continue
        idx = list(adata_sub.var_names).index(gene)
        if hasattr(test_cells.X, 'toarray'):
            h_std = float(test_cells.X[:, idx].toarray().std())
            l_std = float(ref_cells.X[:, idx].toarray().std())
        else:
            h_std = float(test_cells.X[:, idx].std())
            l_std = float(ref_cells.X[:, idx].std())

        if h_std > 0 and l_std > 0:
            ready_genes.append(gene)

    n_dropped = len(sel_genes) - len(ready_genes)
    log(f"  Network-ready genes: {len(ready_genes)} (dropped {n_dropped} constant-in-one-group)")
    filter_log.append(f"Network-ready: {len(ready_genes)} (dropped {n_dropped})")

    # =========================================================================
    # 9. Save results
    # =========================================================================
    banner("[STEP 01.9] Save results")

    # Full stats
    stats_path = os.path.join(DIR_02_DE, "SC_diffexpr_stats.csv")
    stats_df.to_csv(stats_path, index=False)
    log(f"  [SAVE] Full stats ({len(stats_df)} genes) -> {stats_path}")

    # Selected genes
    sel_path = os.path.join(DIR_02_DE, "SC_selected_genes.csv")
    selected_df[["gene", "p_value", "fdr", "log2FC"]].to_csv(sel_path, index=False)
    log(f"  [SAVE] Selected genes ({len(selected_df)}) -> {sel_path}")

    # Network-ready genes
    ready_path = os.path.join(DIR_02_DE, "SC_selected_genes_filtered.csv")
    pd.DataFrame({"gene": ready_genes}).to_csv(ready_path, index=False)
    log(f"  [SAVE] Network-ready genes ({len(ready_genes)}) -> {ready_path}")

    # =========================================================================
    # 10. Volcano plot
    # =========================================================================
    banner("[STEP 01.10] Volcano plot")

    fig, ax = plt.subplots(figsize=(10, 7))
    neg_log_p = -np.log10(stats_df["p_value"].clip(lower=1e-300))

    ns_mask = ~selected_mask
    ax.scatter(stats_df.loc[ns_mask, "log2FC"], neg_log_p[ns_mask],
               c="gray", alpha=0.3, s=8, label="NS", zorder=1)

    up_mask = selected_mask & (stats_df["log2FC"] > 0)
    ax.scatter(stats_df.loc[up_mask, "log2FC"], neg_log_p[up_mask],
               c="#ed6a5a", alpha=0.6, s=14, label=f"Up (n={up_mask.sum()})", zorder=2)

    down_mask = selected_mask & (stats_df["log2FC"] < 0)
    ax.scatter(stats_df.loc[down_mask, "log2FC"], neg_log_p[down_mask],
               c="#7eb0d5", alpha=0.6, s=14, label=f"Down (n={down_mask.sum()})", zorder=2)

    # Label A3 genes
    for a3 in A3_TRACK:
        match = stats_df[stats_df["gene"] == a3]
        if len(match) > 0:
            row = match.iloc[0]
            idx = match.index[0]
            alias = A3_SYMBOL_TO_ALIAS.get(a3, a3)
            ax.annotate(alias, (row["log2FC"], neg_log_p[idx]),
                        fontsize=8, fontweight="bold", color="#2d6a4f",
                        xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{comparison_name} | Volcano (scanpy Wilcoxon, FDR < {FDR_THRESHOLD})")
    ax.legend(fontsize=9, loc="upper right")
    plt.tight_layout()

    volcano_path = os.path.join(DIR_02_DE, "SC_volcano.png")
    plt.savefig(volcano_path, dpi=300)
    plt.close()
    log(f"  [SAVE] Volcano -> {volcano_path}")

    # =========================================================================
    # 11. Manhattan plot
    # =========================================================================
    banner("[STEP 01.11] Manhattan plot")

    sorted_stats = stats_df.sort_values("p_value").reset_index(drop=True)
    y = -np.log10(sorted_stats["p_value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(range(len(y)), y, s=4, c="#7eb0d5", alpha=0.5)
    ax.axhline(-np.log10(FDR_THRESHOLD), ls="--", c="red", lw=0.8,
               label=f"FDR = {FDR_THRESHOLD}")
    ax.set_xlabel("Genes (sorted by p-value)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{comparison_name} | Manhattan (scanpy Wilcoxon)")
    ax.legend()
    plt.tight_layout()

    manhattan_path = os.path.join(DIR_02_DE, "SC_manhattan.png")
    plt.savefig(manhattan_path, dpi=300)
    plt.close()
    log(f"  [SAVE] Manhattan -> {manhattan_path}")

    # =========================================================================
    # 12. Filtering summary
    # =========================================================================
    banner("[STEP 01.12] Filtering summary")

    summary_path = os.path.join(DIR_02_DE, "SC_filtering_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"SC Differential Expression Summary -- {comparison_name}\n")
        f.write(f"{'=' * 60}\n")
        f.write(f"Comparison: {comp['description']}\n")
        f.write(f"Test group: {test_group} ({n_test} cells)\n")
        f.write(f"Reference group: {ref_group} ({n_ref} cells)\n")
        f.write(f"Method: scanpy rank_genes_groups (wilcoxon)\n")
        f.write(f"Selection: FDR < {FDR_THRESHOLD}\n")
        f.write(f"Force-keep A3: {FORCE_KEEP_A3}\n\n")
        for line in filter_log:
            f.write(f"{line}\n")
    log(f"  [SAVE] Summary -> {summary_path}")

    # =========================================================================
    # DONE
    # =========================================================================
    elapsed = datetime.now() - start_time
    banner(f"[STEP 01 COMPLETE] {len(ready_genes)} network-ready genes | "
           f"FDR < {FDR_THRESHOLD} | Elapsed: {elapsed}")

    log(f"  Comparison: {comparison_name}")
    log(f"  Test: {test_group} ({n_test} cells), Ref: {ref_group} ({n_ref} cells)")
    log(f"  Genes tested: {len(stats_df)}")
    log(f"  FDR < {FDR_THRESHOLD}: {n_sig_fdr}")
    log(f"  Network-ready: {len(ready_genes)}")


if __name__ == "__main__":
    main()
