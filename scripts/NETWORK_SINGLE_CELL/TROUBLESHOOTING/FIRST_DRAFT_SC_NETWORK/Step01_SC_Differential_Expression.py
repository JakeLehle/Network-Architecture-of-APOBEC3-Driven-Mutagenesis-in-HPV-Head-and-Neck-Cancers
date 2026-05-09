#!/usr/bin/env python3
"""
Step01_SC_Differential_Expression.py
=====================================

Figure 4 — Step 01: Differential Expression between SBS2-HIGH and SBS2-LOW
basal epithelial cells.

Mirrors Step03_Differential_Expression.py from the Figure 2 (TCGA bulk)
pipeline but simplified: groups are already defined by Step 00 (L-method
elbow detection), so no A3 matching logic is needed here.

NOTE on gene IDs:
  The SC expression matrices use a HYBRID index: gene symbols for annotated
  genes (APOBEC3B, KRT5, etc.) and ENSG IDs only for unannotated transcripts.
  This script uses gene symbols for all A3/biomarker matching. No protein-coding
  filter is applied because Cell Ranger output already contains ~20K protein-
  coding genes (matching the TCGA post-biotype-filter starting point).

Workflow:
  1. Load HIGH and LOW expression matrices (genes × cells)
  2. Find common genes (intersection of row indices)
  3. Filter low-detection genes (detected in >= MIN_CELLS_DETECTED cells)
  4. Log1p transform
  5. Wilcoxon rank-sum test per gene (HIGH vs LOW)
  6. BH-FDR correction
  7. log2FC = mean(log1p(HIGH)) - mean(log1p(LOW))
  8. Select genes: raw p < 0.05
  9. Force-keep A3 family genes
  10. Save results + diagnostic plots

Input:
  data/FIG_4/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv
  data/FIG_4/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv

Output (→ data/FIG_4/02_differential_expression/):
  SC_diffexpr_stats.csv          — per-gene Wilcoxon results
  SC_selected_genes.csv          — genes passing thresholds
  SC_selected_genes_filtered.csv — network-ready genes (non-constant in both)
  SC_volcano.png                 — volcano plot
  SC_manhattan.png               — Manhattan-style plot
  SC_filtering_summary.txt       — step-by-step filtering log

Usage:
  conda run -n NETWORK python Step01_SC_Differential_Expression.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime

from network_config_SC import (
    DIR_01_GROUPS, DIR_02_DE,
    HIGH_EXPR_PATH, LOW_EXPR_PATH,
    A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    MIN_CELLS_DETECTED, RAW_P_THRESHOLD, LOGFC_THRESHOLD,
    FORCE_KEEP_A3,
    banner, log, ensure_dir
)


# =============================================================================
# BH-FDR CORRECTION
# =============================================================================

def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction. Returns array of adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    adjusted = pvals * n / ranked
    # Enforce monotonicity (step-down)
    adjusted = np.minimum.accumulate(adjusted[np.argsort(ranked)[::-1]])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    return adjusted


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    start_time = datetime.now()
    banner("[STEP 01] Single-Cell Differential Expression — HIGH vs LOW")
    log(f"Start time: {start_time}")

    ensure_dir(DIR_02_DE)

    # Track filtering for summary
    filter_log = []

    # =========================================================================
    # 1. Load expression matrices
    # =========================================================================
    banner("[STEP 01.1] Load expression matrices")

    log(f"Reading HIGH matrix: {HIGH_EXPR_PATH}")
    high_df = pd.read_csv(HIGH_EXPR_PATH, sep="\t", index_col=0)
    log(f"  HIGH shape: {high_df.shape[0]} genes × {high_df.shape[1]} cells")

    log(f"Reading LOW matrix: {LOW_EXPR_PATH}")
    low_df = pd.read_csv(LOW_EXPR_PATH, sep="\t", index_col=0)
    log(f"  LOW shape: {low_df.shape[0]} genes × {low_df.shape[1]} cells")

    filter_log.append(f"HIGH input: {high_df.shape[0]} genes × {high_df.shape[1]} cells")
    filter_log.append(f"LOW input: {low_df.shape[0]} genes × {low_df.shape[1]} cells")

    # Report gene ID format
    sample_ids = list(high_df.index[:10])
    n_ensg = sum(1 for g in high_df.index if str(g).startswith("ENSG"))
    n_symbol = high_df.shape[0] - n_ensg
    log(f"  Gene ID format: {n_symbol} symbols + {n_ensg} ENSG IDs (hybrid)")
    log(f"  Sample IDs: {sample_ids}")
    filter_log.append(f"Gene ID format: {n_symbol} symbols + {n_ensg} ENSG IDs")

    # =========================================================================
    # 2. Find common genes
    # =========================================================================
    banner("[STEP 01.2] Find common genes")

    common_genes = sorted(set(high_df.index) & set(low_df.index))
    log(f"  Common genes: {len(common_genes)}")
    filter_log.append(f"Common genes (intersection): {len(common_genes)}")

    high_df = high_df.loc[common_genes]
    low_df = low_df.loc[common_genes]

    # =========================================================================
    # 3. Filter low-detection genes
    # =========================================================================
    banner("[STEP 01.3] Filter low-detection genes")

    # Count cells where gene is detected (expression > 0), across both groups
    high_detected = (high_df > 0).sum(axis=1)
    low_detected = (low_df > 0).sum(axis=1)
    total_detected = high_detected + low_detected

    pass_detection = total_detected >= MIN_CELLS_DETECTED
    n_before = len(common_genes)
    filtered_genes = sorted(pass_detection[pass_detection].index.tolist())

    # Force-keep A3 genes
    if FORCE_KEEP_A3:
        for a3 in A3_GENES:
            if a3 in common_genes and a3 not in filtered_genes:
                filtered_genes.append(a3)
                alias = A3_ID_TO_ALIAS.get(a3, a3)
                log(f"  Force-kept {alias} (detected in {total_detected.get(a3, 0)} cells)")
        filtered_genes = sorted(set(filtered_genes))

    n_after = len(filtered_genes)
    log(f"  Detection filter (>= {MIN_CELLS_DETECTED} cells): {n_before} → {n_after} genes")
    filter_log.append(f"Detection filter (>= {MIN_CELLS_DETECTED} cells): {n_before} → {n_after}")

    high_df = high_df.loc[filtered_genes]
    low_df = low_df.loc[filtered_genes]

    # Report A3 gene status
    banner("[STEP 01.3b] A3 gene detection status")
    for a3 in A3_GENES:
        alias = A3_ID_TO_ALIAS.get(a3, a3)
        if a3 in filtered_genes:
            h_det = int(high_detected.get(a3, 0))
            l_det = int(low_detected.get(a3, 0))
            log(f"  {alias} ({a3}): HIGH={h_det} cells, LOW={l_det} cells — INCLUDED")
        elif a3 in common_genes:
            h_det = int(high_detected.get(a3, 0))
            l_det = int(low_detected.get(a3, 0))
            log(f"  {alias} ({a3}): HIGH={h_det} cells, LOW={l_det} cells — BELOW THRESHOLD")
        else:
            log(f"  {alias} ({a3}): NOT in common gene set")

    # Report biomarker gene status
    banner("[STEP 01.3c] Biomarker gene status")
    for bm in BIOMARKERS:
        if bm in filtered_genes:
            log(f"  {bm}: INCLUDED")
        elif bm in common_genes:
            log(f"  {bm}: below detection threshold")
        else:
            log(f"  {bm}: not in expression matrix")

    # =========================================================================
    # 4. Log1p transform
    # =========================================================================
    banner("[STEP 01.4] Log1p transform")

    high_log = np.log1p(high_df)
    low_log = np.log1p(low_df)
    log(f"  Applied log1p to {len(filtered_genes)} genes × {high_df.shape[1]} + {low_df.shape[1]} cells")

    # =========================================================================
    # 5. Wilcoxon rank-sum test
    # =========================================================================
    banner("[STEP 01.5] Wilcoxon rank-sum test (HIGH vs LOW)")

    results = []
    n_genes = len(filtered_genes)
    checkpoint = max(1, n_genes // 20)

    for idx, gene in enumerate(filtered_genes):
        if (idx + 1) % checkpoint == 0 or (idx + 1) == n_genes:
            log(f"  Testing gene {idx + 1}/{n_genes} ({100*(idx+1)/n_genes:.0f}%)...")

        h_vals = high_log.loc[gene].values
        l_vals = low_log.loc[gene].values

        # Skip if constant in both groups
        if np.std(h_vals) == 0 and np.std(l_vals) == 0:
            results.append({
                "gene": gene,
                "p_value": 1.0,
                "log2FC": 0.0,
                "mean_high": float(np.mean(h_vals)),
                "mean_low": float(np.mean(l_vals)),
                "pct_high": float((high_df.loc[gene] > 0).mean()),
                "pct_low": float((low_df.loc[gene] > 0).mean()),
                "test_status": "constant"
            })
            continue

        try:
            stat, pval = mannwhitneyu(h_vals, l_vals, alternative="two-sided")
        except ValueError:
            pval = 1.0

        log2fc = float(np.mean(h_vals) - np.mean(l_vals))

        results.append({
            "gene": gene,
            "p_value": pval,
            "log2FC": log2fc,
            "mean_high": float(np.mean(h_vals)),
            "mean_low": float(np.mean(l_vals)),
            "pct_high": float((high_df.loc[gene] > 0).mean()),
            "pct_low": float((low_df.loc[gene] > 0).mean()),
            "test_status": "tested"
        })

    stats_df = pd.DataFrame(results)
    n_tested = (stats_df["test_status"] == "tested").sum()
    n_constant = (stats_df["test_status"] == "constant").sum()
    log(f"  Tested: {n_tested} genes, Constant (skipped): {n_constant} genes")
    filter_log.append(f"Wilcoxon tested: {n_tested} (constant: {n_constant})")

    # =========================================================================
    # 6. BH-FDR correction
    # =========================================================================
    banner("[STEP 01.6] BH-FDR correction")

    stats_df["fdr"] = bh_fdr(stats_df["p_value"].values)
    log(f"  FDR < 0.05: {(stats_df['fdr'] < 0.05).sum()} genes")
    log(f"  FDR < 0.10: {(stats_df['fdr'] < 0.10).sum()} genes")
    log(f"  FDR < 0.25: {(stats_df['fdr'] < 0.25).sum()} genes")

    # =========================================================================
    # 7. Gene selection (raw p < 0.05, matching Figure 2 approach)
    # =========================================================================
    banner("[STEP 01.7] Gene selection")

    # Primary selection: raw p-value
    pass_p = stats_df["p_value"] < RAW_P_THRESHOLD
    if LOGFC_THRESHOLD > 0:
        pass_fc = stats_df["log2FC"].abs() >= LOGFC_THRESHOLD
        selected_mask = pass_p & pass_fc
    else:
        selected_mask = pass_p

    n_selected = selected_mask.sum()
    log(f"  Raw p < {RAW_P_THRESHOLD}: {pass_p.sum()} genes")
    log(f"  Selected (p + FC filter): {n_selected} genes")

    # Force-keep A3 genes
    if FORCE_KEEP_A3:
        a3_forced = 0
        for a3 in A3_GENES:
            mask = stats_df["gene"] == a3
            if mask.any() and not selected_mask[mask].values[0]:
                selected_mask[mask] = True
                alias = A3_ID_TO_ALIAS.get(a3, a3)
                p = stats_df.loc[mask, "p_value"].values[0]
                log(f"  Force-kept {alias} (p={p:.4f})")
                a3_forced += 1
        if a3_forced > 0:
            log(f"  Force-kept {a3_forced} A3 genes")

    stats_df["selected"] = selected_mask
    selected_df = stats_df[stats_df["selected"]].copy()
    log(f"  Final selected genes: {len(selected_df)}")
    filter_log.append(f"Selected genes (raw p < {RAW_P_THRESHOLD}): {len(selected_df)}")

    # =========================================================================
    # 8. Network-readiness filter
    # =========================================================================
    banner("[STEP 01.8] Network-readiness filter")
    # Ensure selected genes are non-constant in both groups (needed for Spearman)

    sel_genes = selected_df["gene"].tolist()
    ready_genes = []
    for gene in sel_genes:
        h_std = float(high_df.loc[gene].std())
        l_std = float(low_df.loc[gene].std())
        if h_std > 0 and l_std > 0:
            ready_genes.append(gene)

    n_dropped = len(sel_genes) - len(ready_genes)
    log(f"  Network-ready genes: {len(ready_genes)} (dropped {n_dropped} constant-in-one-group)")
    filter_log.append(f"Network-ready genes: {len(ready_genes)} (dropped {n_dropped})")

    # =========================================================================
    # 9. Save results
    # =========================================================================
    banner("[STEP 01.9] Save results")

    # Full stats
    stats_path = os.path.join(DIR_02_DE, "SC_diffexpr_stats.csv")
    stats_df.to_csv(stats_path, index=False)
    log(f"  [SAVE] Full stats ({len(stats_df)} genes) → {stats_path}")

    # Selected genes
    sel_path = os.path.join(DIR_02_DE, "SC_selected_genes.csv")
    selected_df.to_csv(sel_path, index=False)
    log(f"  [SAVE] Selected genes ({len(selected_df)}) → {sel_path}")

    # Network-ready genes
    ready_df = selected_df[selected_df["gene"].isin(ready_genes)].copy()
    ready_path = os.path.join(DIR_02_DE, "SC_selected_genes_filtered.csv")
    ready_df.to_csv(ready_path, index=False)
    log(f"  [SAVE] Network-ready genes ({len(ready_df)}) → {ready_path}")

    # =========================================================================
    # 10. Volcano plot
    # =========================================================================
    banner("[STEP 01.10] Volcano plot")

    fig, ax = plt.subplots(figsize=(10, 7))

    # Color by significance
    sig = stats_df["p_value"] < RAW_P_THRESHOLD
    nonsig = ~sig

    ax.scatter(stats_df.loc[nonsig, "log2FC"],
               -np.log10(stats_df.loc[nonsig, "p_value"].clip(lower=1e-300)),
               s=4, c="lightgray", alpha=0.5, label=f"NS (n={nonsig.sum()})")

    ax.scatter(stats_df.loc[sig, "log2FC"],
               -np.log10(stats_df.loc[sig, "p_value"].clip(lower=1e-300)),
               s=6, c="steelblue", alpha=0.6, label=f"p<{RAW_P_THRESHOLD} (n={sig.sum()})")

    # Highlight A3 genes
    for a3 in A3_GENES:
        row = stats_df[stats_df["gene"] == a3]
        if len(row) == 1:
            alias = A3_ID_TO_ALIAS.get(a3, a3)
            x = float(row["log2FC"].iloc[0])
            y = -np.log10(float(row["p_value"].clip(lower=1e-300).iloc[0]))
            ax.scatter(x, y, s=60, c="firebrick", edgecolors="black",
                       linewidths=0.5, zorder=5)
            ax.annotate(alias, (x, y), fontsize=7, fontweight="bold",
                        xytext=(5, 5), textcoords="offset points")

    ax.axhline(-np.log10(RAW_P_THRESHOLD), ls="--", c="red", lw=0.8, alpha=0.6)
    ax.set_xlabel("log2FC (mean log1p HIGH − LOW)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"SC Basal Cells | Wilcoxon Rank-Sum (n={high_df.shape[1]} vs {low_df.shape[1]})")
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()

    volcano_path = os.path.join(DIR_02_DE, "SC_volcano.png")
    plt.savefig(volcano_path, dpi=300)
    plt.close()
    log(f"  [SAVE] Volcano → {volcano_path}")

    # =========================================================================
    # 11. Manhattan plot
    # =========================================================================
    banner("[STEP 01.11] Manhattan plot")

    sorted_stats = stats_df.sort_values("p_value").reset_index(drop=True)
    y = -np.log10(sorted_stats["p_value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(range(len(y)), y, s=4, c="steelblue", alpha=0.5)
    ax.axhline(-np.log10(RAW_P_THRESHOLD), ls="--", c="red", lw=0.8)
    ax.set_xlabel("Genes (sorted by p-value)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("SC Basal Cells | Manhattan — Wilcoxon Rank-Sum")
    plt.tight_layout()

    manhattan_path = os.path.join(DIR_02_DE, "SC_manhattan.png")
    plt.savefig(manhattan_path, dpi=300)
    plt.close()
    log(f"  [SAVE] Manhattan → {manhattan_path}")

    # =========================================================================
    # 12. Filtering summary
    # =========================================================================
    banner("[STEP 01.12] Filtering summary")

    summary_path = os.path.join(DIR_02_DE, "SC_filtering_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Step 01 — Single-Cell Differential Expression Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("Filtering steps:\n")
        for line in filter_log:
            f.write(f"  {line}\n")
        f.write(f"\nA3 gene results:\n")
        for a3 in A3_GENES:
            alias = A3_ID_TO_ALIAS.get(a3, a3)
            row = stats_df[stats_df["gene"] == a3]
            if len(row) == 1:
                p = float(row["p_value"].iloc[0])
                fc = float(row["log2FC"].iloc[0])
                sel = "SELECTED" if bool(row["selected"].values[0]) else "not selected"
                f.write(f"  {alias}: p={p:.6f}, log2FC={fc:.4f} [{sel}]\n")
            else:
                f.write(f"  {alias}: NOT FOUND in expression data\n")

    log(f"  [SAVE] Summary → {summary_path}")

    # =========================================================================
    # DONE
    # =========================================================================
    elapsed = datetime.now() - start_time
    banner(f"[STEP 01 COMPLETE] {len(ready_genes)} network-ready genes | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
