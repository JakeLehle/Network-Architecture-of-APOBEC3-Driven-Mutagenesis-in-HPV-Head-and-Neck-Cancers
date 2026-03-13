#!/usr/bin/env python3
"""
Step03_Differential_Expression.py

For each cancer type, compare gene expression between high-SBS2 and
low-SBS2 tumors using Welch's t-test. Select significant genes,
always retaining A3 family members. Produce volcano and Manhattan plots.

Corresponds to original pipeline Step 12.

Input:
    02_merged_with_SBS/TCGA_merged_expression_SBS.pkl

Output (-> data/FIG_2/03_differential_expression/{cancer_type}/):
    {cancer_type}_diffexpr_stats.csv         — per-gene t-test results
    {cancer_type}_selected_genes.csv         — genes passing thresholds
    {cancer_type}_volcano.png                — volcano plot
    {cancer_type}_manhattan.png              — Manhattan-style plot
    {cancer_type}_cancer_df2.pkl         — filtered dataframe (clinical + selected genes + SBS2/A3)

Usage:
    python Step03_Differential_Expression.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from network_config import (
    CANCER_TYPES, CLINICAL_COLS, A3_GENES, BIOMARKERS,
    SBS2_HIGH_PERCENTILE, SBS2_LOW_PERCENTILE, P_VALUE_THRESHOLD,
    FORCE_KEEP_A3,
    DIR_02_MERGED, DIR_03_DIFFEXPR,
    banner, log, ensure_dir
)


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
    # Enforce monotonicity
    adjusted = np.minimum.accumulate(adjusted[np.argsort(ranked)[::-1]])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    return adjusted


# =============================================================================
# LOAD merged data
# =============================================================================
banner("[STEP 12] Load merged expression + SBS data")

merged_path = os.path.join(DIR_02_MERGED, "TCGA_merged_expression_SBS.pkl")
log(f"[STEP 12] Reading: {merged_path}")
merged = pd.read_pickle(merged_path)
log(f"[STEP 12] Shape: {merged.shape}")

# Identify all ENSG gene columns
all_gene_cols = [c for c in merged.columns if isinstance(c, str) and c.startswith("ENSG")]
log(f"[STEP 12] Total ENSG gene columns: {len(all_gene_cols)}")

# Check A3 and biomarker presence
present_a3 = [g for g in A3_GENES if g in merged.columns]
missing_a3 = [g for g in A3_GENES if g not in merged.columns]
present_bm = [g for g in BIOMARKERS if g in merged.columns]
log(f"[STEP 12] A3 genes present: {len(present_a3)} | missing: {missing_a3}")
log(f"[STEP 12] Biomarkers present: {len(present_bm)}")


# =============================================================================
# LOOP over cancer types
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"[CANCER] {cancer_type}", char="=")

    cancer_dir = ensure_dir(os.path.join(DIR_03_DIFFEXPR, cancer_type))

    # ---- Filter to this cancer type
    cancer_df = merged[merged["Project_ID"] == cancer_type].copy()
    log(f"[STEP 12.1] Samples for {cancer_type}: {len(cancer_df)}")

    if len(cancer_df) < 25:
        log("[SKIP] Too few samples (<25)")
        continue

    # ---- Prepare SBS2
    if "SBS2" not in cancer_df.columns:
        log("[SKIP] SBS2 column missing")
        continue

    cancer_df["SBS2"] = pd.to_numeric(cancer_df["SBS2"], errors="coerce")
    cancer_df = cancer_df.dropna(subset=["SBS2"]).reset_index(drop=True)

    if len(cancer_df) < 40:
        log("[SKIP] Too few valid SBS2 samples (<40)")
        continue

    # ---- Define high/low SBS2 groups (for t-test)
    high_thr = cancer_df["SBS2"].quantile(SBS2_HIGH_PERCENTILE)
    low_thr  = cancer_df["SBS2"].quantile(SBS2_LOW_PERCENTILE)

    high_df = cancer_df[cancer_df["SBS2"] >= high_thr]
    low_df  = cancer_df[cancer_df["SBS2"] <= low_thr]

    log(f"[STEP 12.2] SBS2 thresholds: low20={low_thr:.6f}, high80={high_thr:.6f}")
    log(f"[STEP 12.2] High group: {len(high_df)} | Low group: {len(low_df)}")

    if len(high_df) < 10 or len(low_df) < 10:
        log("[SKIP] Too few samples per group (<10)")
        continue

    # ---- Run Welch's t-test for each gene
    banner(f"[STEP 12.3] Welch's t-test ({cancer_type})")

    results = []
    for gene in all_gene_cols:
        vals_high = pd.to_numeric(high_df[gene], errors="coerce").dropna().values
        vals_low  = pd.to_numeric(low_df[gene], errors="coerce").dropna().values

        if len(vals_high) < 3 or len(vals_low) < 3:
            continue

        t_stat, p_val = ttest_ind(vals_high, vals_low, equal_var=False)

        mean_high = np.mean(vals_high)
        mean_low  = np.mean(vals_low)
        log2fc = np.log2((mean_high + 1) / (mean_low + 1))

        results.append({
            "gene": gene,
            "t_stat": t_stat,
            "p_value": p_val,
            "mean_high_SBS2": mean_high,
            "mean_low_SBS2": mean_low,
            "log2FC": log2fc,
        })

    stats_df = pd.DataFrame(results)
    log(f"[STEP 12.3] Genes tested: {len(stats_df)}")

    # ---- FDR correction
    stats_df["p_adj"] = bh_fdr(stats_df["p_value"].values)
    stats_df = stats_df.sort_values("p_value").reset_index(drop=True)

    n_sig_raw = (stats_df["p_value"] < P_VALUE_THRESHOLD).sum()
    n_sig_fdr = (stats_df["p_adj"] < 0.05).sum()
    log(f"[STEP 12.4] Significant (raw p < {P_VALUE_THRESHOLD}): {n_sig_raw}")
    log(f"[STEP 12.4] Significant (FDR < 0.05): {n_sig_fdr}")

    # ---- Select genes
    selected = stats_df[stats_df["p_value"] < P_VALUE_THRESHOLD]["gene"].tolist()

    # Force-keep A3 genes
    if FORCE_KEEP_A3:
        for g in A3_GENES:
            if g in stats_df["gene"].values and g not in selected:
                selected.append(g)
                log(f"[STEP 12.5] Force-kept A3 gene: {g}")

    selected = list(dict.fromkeys(selected))  # deduplicate, preserve order
    log(f"[STEP 12.5] Total selected genes: {len(selected)}")

    # Report biomarker status
    tested_set = set(stats_df["gene"])
    log("[STEP 12.6] Biomarkers tested/selected:")
    for g in BIOMARKERS:
        tested = "TESTED" if g in tested_set else "NOT_TESTED"
        sel = "SELECTED" if g in selected else "NOT_SELECTED"
        log(f"  {g}: {tested} / {sel}")

    # ---- Save stats
    stats_csv = os.path.join(cancer_dir, f"{cancer_type}_diffexpr_stats.csv")
    stats_df.to_csv(stats_csv, index=False)
    log(f"[SAVE] Stats -> {stats_csv}")

    # ---- Save selected gene list
    sel_csv = os.path.join(cancer_dir, f"{cancer_type}_selected_genes.csv")
    pd.DataFrame({"gene": selected}).to_csv(sel_csv, index=False)
    log(f"[SAVE] Selected genes -> {sel_csv}")

    # ---- Volcano plot
    banner(f"[STEP 12.6A] Volcano plot ({cancer_type})")

    fig, ax = plt.subplots(figsize=(10, 7))
    neg_log_p = -np.log10(stats_df["p_value"].clip(lower=1e-300))

    # Non-significant
    ns_mask = stats_df["p_value"] >= P_VALUE_THRESHOLD
    ax.scatter(stats_df.loc[ns_mask, "log2FC"], neg_log_p[ns_mask],
               c="gray", alpha=0.4, s=8, label="NS")

    # Significant
    sig_mask = ~ns_mask
    ax.scatter(stats_df.loc[sig_mask, "log2FC"], neg_log_p[sig_mask],
               c="firebrick", alpha=0.6, s=12, label=f"p < {P_VALUE_THRESHOLD}")

    ax.axhline(-np.log10(P_VALUE_THRESHOLD), ls="--", c="black", lw=0.8)
    ax.set_xlabel("log2(Fold Change: High SBS2 / Low SBS2)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{cancer_type} | Volcano (High vs Low SBS2)")
    ax.legend(fontsize=9)
    plt.tight_layout()

    volcano_path = os.path.join(cancer_dir, f"{cancer_type}_volcano.png")
    plt.savefig(volcano_path, dpi=300)
    plt.close()
    log(f"[SAVE] Volcano -> {volcano_path}")

    # ---- Manhattan-style plot
    banner(f"[STEP 12.6B] Manhattan plot ({cancer_type})")

    sorted_stats = stats_df.sort_values("p_value").reset_index(drop=True)
    y = -np.log10(sorted_stats["p_value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(range(len(y)), y, s=4, c="steelblue", alpha=0.5)
    ax.axhline(-np.log10(P_VALUE_THRESHOLD), ls="--", c="red", lw=0.8)
    ax.set_xlabel("Genes (sorted by p-value)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{cancer_type} | Manhattan (High vs Low SBS2)")
    plt.tight_layout()

    manhattan_path = os.path.join(cancer_dir, f"{cancer_type}_manhattan.png")
    plt.savefig(manhattan_path, dpi=300)
    plt.close()
    log(f"[SAVE] Manhattan -> {manhattan_path}")

    # ---- Build cancer_df2 (clinical + selected genes + signal columns)
    banner(f"[STEP 12.7] Build cancer_df2 ({cancer_type})")

    keep_cols = list(CLINICAL_COLS)

    # Add SBS2 and A3 alias columns for downstream scoring
    for c in ["SBS2", "A3A", "A3B", "A3C", "A3H"]:
        if c in cancer_df.columns and c not in keep_cols:
            keep_cols.append(c)

    # Add selected genes
    keep_cols += [g for g in selected if g in cancer_df.columns and g not in keep_cols]

    cancer_df2 = cancer_df[keep_cols].copy()
    log(f"[STEP 12.7] cancer_df2 shape: {cancer_df2.shape}")

    # Sanity checks
    dup_cols = cancer_df2.columns[cancer_df2.columns.duplicated()].tolist()
    if dup_cols:
        log(f"[WARNING] Duplicate columns in cancer_df2: {dup_cols}")
    else:
        log("[OK] No duplicate columns")

    # Save
    df2_path = os.path.join(cancer_dir, f"{cancer_type}_cancer_df2.pkl")
    cancer_df2.to_pickle(df2_path)
    log(f"[SAVE] cancer_df2 -> {df2_path}")

    log(f"\n[STEP 12 COMPLETE for {cancer_type}]")
    print(f"  Selected genes: {len(selected)}")
    print(f"  Significant (raw): {n_sig_raw}")
    print(f"  Significant (FDR): {n_sig_fdr}")

banner("[STEP 03 COMPLETE — ALL CANCER TYPES]")
