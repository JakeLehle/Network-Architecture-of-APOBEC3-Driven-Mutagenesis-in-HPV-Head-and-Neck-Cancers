#!/usr/bin/env python3
"""
Step04_Define_Groups.py

Compute the APOBEC activity score (A3A + A3B, normalized), restrict
to high-A3 tumors, then split by SBS2 into TOP (high SBS2) and
BOTTOM (low SBS2) groups. These groups are the basis for differential
network analysis.

Produces the SBS2-vs-A3-score selection plot → Figure 2a candidate.

Corresponds to original pipeline Step 13.

Input:
    03_differential_expression/{cancer_type}/{cancer_type}_cancer_df2.parquet
    03_differential_expression/{cancer_type}/{cancer_type}_selected_genes.csv

Output (-> data/FIG_2/04_TOP_BOTTOM_groups/{cancer_type}/):
    {ct}_TOP_group.parquet / .csv         — TOP group (high SBS2, high A3)
    {ct}_BOTTOM_group.parquet / .csv      — BOTTOM group (low SBS2, high A3)
    {ct}_selected_genes_filtered.csv      — genes present in both groups
    {ct}_SBS2_vs_A3score_selection.png    — selection plot (Figure 2a)
    {ct}_group_summary.txt               — group size/threshold summary

Usage:
    python Step04_Define_Groups.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from network_config import (
    CANCER_TYPES,
    A3_SCORE_PERCENTILE, GROUP_SBS2_HIGH_PERCENTILE, GROUP_SBS2_LOW_PERCENTILE,
    MIN_GROUP_SIZE,
    DIR_03_DIFFEXPR, DIR_04_GROUPS,
    banner, log, ensure_dir
)


def compute_a3_score(df, a_col="A3A", b_col="A3B"):
    """
    Compute APOBEC3 activity score = normalized(A3A) + normalized(A3B).
    Each gene is z-scored within the cohort before summing.
    """
    a = pd.to_numeric(df[a_col], errors="coerce")
    b = pd.to_numeric(df[b_col], errors="coerce")

    a_z = (a - a.mean()) / (a.std() + 1e-12)
    b_z = (b - b.mean()) / (b.std() + 1e-12)

    return a_z + b_z


# =============================================================================
# LOOP over cancer types
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"[STEP 13] Define TOP/BOTTOM groups — {cancer_type}", char="=")

    cancer_dir = ensure_dir(os.path.join(DIR_04_GROUPS, cancer_type))

    # ---- Load cancer_df2 from Step 03
    df2_path = os.path.join(DIR_03_DIFFEXPR, cancer_type, f"{cancer_type}_cancer_df2.parquet")
    if not os.path.exists(df2_path):
        log(f"[SKIP] cancer_df2 not found: {df2_path}")
        continue

    cancer_df2 = pd.read_parquet(df2_path)
    log(f"[STEP 13] Loaded cancer_df2: {cancer_df2.shape}")

    # ---- Load selected genes
    sel_path = os.path.join(DIR_03_DIFFEXPR, cancer_type, f"{cancer_type}_selected_genes.csv")
    selected_genes = pd.read_csv(sel_path)["gene"].tolist()
    log(f"[STEP 13] Selected genes from Step 03: {len(selected_genes)}")

    # ---- Check required columns
    needed = ["SBS2", "A3A", "A3B"]
    missing = [c for c in needed if c not in cancer_df2.columns]
    if missing:
        log(f"[SKIP] Missing columns: {missing}")
        continue

    # ---- Clean and compute A3 score
    df_s13 = cancer_df2.copy()
    for c in needed:
        df_s13[c] = pd.to_numeric(df_s13[c], errors="coerce")
    df_s13 = df_s13.dropna(subset=needed)
    log(f"[STEP 13.1] After cleanup: {len(df_s13)} samples")

    if len(df_s13) < 30:
        log("[SKIP] Too few samples (<30)")
        continue

    df_s13["A3_score"] = compute_a3_score(df_s13, "A3A", "A3B")

    # ---- Restrict to high-A3 tumors (above median)
    a3_thr = float(df_s13["A3_score"].quantile(A3_SCORE_PERCENTILE))
    high_a3 = df_s13[df_s13["A3_score"] >= a3_thr].copy()

    log(f"[STEP 13.2] A3_score threshold (median): {a3_thr:.6f}")
    log(f"[STEP 13.2] High-A3 samples: {len(high_a3)} / {len(df_s13)}")

    if len(high_a3) < 20:
        log("[SKIP] Too few high-A3 samples (<20)")
        continue

    # ---- Split by SBS2 within high-A3
    sbs_hi = float(high_a3["SBS2"].quantile(GROUP_SBS2_HIGH_PERCENTILE))
    sbs_lo = float(high_a3["SBS2"].quantile(GROUP_SBS2_LOW_PERCENTILE))

    group_top = high_a3[high_a3["SBS2"] >= sbs_hi].copy()
    group_bot = high_a3[high_a3["SBS2"] <= sbs_lo].copy()

    log(f"[STEP 13.3] SBS2 thresholds (within high-A3): low20={sbs_lo:.6f}, high80={sbs_hi:.6f}")
    log(f"[STEP 13.3] TOP: {len(group_top)} | BOTTOM: {len(group_bot)}")

    if len(group_top) < MIN_GROUP_SIZE or len(group_bot) < MIN_GROUP_SIZE:
        log(f"[SKIP] Groups too small (need >= {MIN_GROUP_SIZE})")
        continue

    # ---- Filter selected genes to those present in both groups
    selected_filtered = [g for g in selected_genes
                         if g in group_top.columns and g in group_bot.columns]
    selected_filtered = list(dict.fromkeys(selected_filtered))

    # Remove genes constant in BOTH groups
    top_nunique = group_top[selected_filtered].apply(pd.to_numeric, errors="coerce").nunique()
    bot_nunique = group_bot[selected_filtered].apply(pd.to_numeric, errors="coerce").nunique()
    const_both = (top_nunique <= 1) & (bot_nunique <= 1)
    n_const = const_both.sum()
    if n_const > 0:
        log(f"[STEP 13.4] Removing {n_const} genes constant in both groups")
        selected_filtered = [g for g in selected_filtered if not const_both.get(g, False)]

    log(f"[STEP 13.4] Genes after filtering: {len(selected_filtered)}")

    # ---- Selection plot (Figure 2a candidate)
    banner(f"[STEP 13.5] Selection plot ({cancer_type})")

    fig, ax = plt.subplots(figsize=(9, 7))

    # All samples (gray)
    ax.scatter(df_s13["A3_score"], df_s13["SBS2"],
               c="lightgray", s=15, alpha=0.5, label="All samples", zorder=1)

    # Low-A3 region
    ax.axvline(a3_thr, ls="--", c="black", lw=0.8, alpha=0.5)

    # SBS2 thresholds within high-A3
    ax.axhline(sbs_hi, ls=":", c="firebrick", lw=0.8, alpha=0.6)
    ax.axhline(sbs_lo, ls=":", c="steelblue", lw=0.8, alpha=0.6)

    # TOP group (red)
    ax.scatter(group_top["A3_score"], group_top["SBS2"],
               c="firebrick", s=30, alpha=0.8, edgecolors="black", linewidths=0.3,
               label=f"TOP (n={len(group_top)})", zorder=3)

    # BOTTOM group (blue)
    ax.scatter(group_bot["A3_score"], group_bot["SBS2"],
               c="steelblue", s=30, alpha=0.8, edgecolors="black", linewidths=0.3,
               label=f"BOTTOM (n={len(group_bot)})", zorder=3)

    ax.set_xlabel("A3 Activity Score (z-scored A3A + A3B)", fontsize=12)
    ax.set_ylabel("SBS2 Weight", fontsize=12)
    ax.set_title(f"{cancer_type} | SBS2 vs A3 Score — TOP/BOTTOM Selection", fontsize=13)
    ax.legend(fontsize=10, loc="upper left")
    plt.tight_layout()

    plot_path = os.path.join(cancer_dir, f"{cancer_type}_SBS2_vs_A3score_selection.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    log(f"[SAVE] Selection plot -> {plot_path}")

    # ---- Save group DataFrames
    top_parquet = os.path.join(cancer_dir, f"{cancer_type}_TOP_group.parquet")
    bot_parquet = os.path.join(cancer_dir, f"{cancer_type}_BOTTOM_group.parquet")
    group_top.to_parquet(top_parquet, index=False)
    group_bot.to_parquet(bot_parquet, index=False)
    log(f"[SAVE] TOP group -> {top_parquet}")
    log(f"[SAVE] BOTTOM group -> {bot_parquet}")

    # Also CSV for easy inspection
    group_top.to_csv(top_parquet.replace(".parquet", ".csv"), index=False)
    group_bot.to_csv(bot_parquet.replace(".parquet", ".csv"), index=False)

    # Save filtered gene list
    genes_path = os.path.join(cancer_dir, f"{cancer_type}_selected_genes_filtered.csv")
    pd.DataFrame({"gene": selected_filtered}).to_csv(genes_path, index=False)
    log(f"[SAVE] Filtered genes -> {genes_path}")

    # ---- Summary
    summary_path = os.path.join(cancer_dir, f"{cancer_type}_group_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Step 04 — TOP/BOTTOM Group Definition ({cancer_type})\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total samples (after cleanup): {len(df_s13)}\n")
        f.write(f"A3 score threshold (median): {a3_thr:.6f}\n")
        f.write(f"High-A3 samples: {len(high_a3)}\n")
        f.write(f"SBS2 threshold (high, 80th): {sbs_hi:.6f}\n")
        f.write(f"SBS2 threshold (low, 20th): {sbs_lo:.6f}\n")
        f.write(f"TOP group size: {len(group_top)}\n")
        f.write(f"BOTTOM group size: {len(group_bot)}\n")
        f.write(f"Selected genes (after filter): {len(selected_filtered)}\n")
    log(f"[SAVE] Summary -> {summary_path}")

    # ---- Console summary
    print(f"\n[STEP 13 SUMMARY — {cancer_type}]")
    print(f"  Clean samples: {len(df_s13)}")
    print(f"  High-A3 samples: {len(high_a3)}")
    print(f"  TOP: {len(group_top)} | BOTTOM: {len(group_bot)}")
    print(f"  Selected genes: {len(selected_filtered)}")

banner("[STEP 04 COMPLETE — ALL CANCER TYPES]")
