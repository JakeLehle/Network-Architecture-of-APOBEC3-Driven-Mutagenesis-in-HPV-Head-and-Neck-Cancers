#!/usr/bin/env python3
"""
Step03_Differential_Expression.py

Gene filtering, group definition, and differential expression analysis
aligned with the single-cell DEG workflow for methodological consistency.

Workflow (mirrors single-cell pipeline):
  1. Filter to protein-coding genes only (using biotype from Step01)
  2. Filter low-expression genes (detected in ≥20 samples)
  3. Filter to cancer type (HNSC)
  4. Define A3-controlled groups:
     a. Sum raw A3A + A3B per tumor
     b. Keep tumors above median of summed A3
     c. Within high-A3: top 25% SBS2 = HIGH, bottom 25% = LOW
  5. Log1p transform expression values
  6. Wilcoxon rank-sum test (matches scanpy method='wilcoxon')
  7. BH-FDR correction
  8. Select genes: FDR < 0.05 AND |log2FC| > 1.5
  9. Force-keep A3 family genes

The HIGH/LOW groups defined here are saved and inherited by Step04
for network construction (single source of truth for group definitions).

Corresponds to original pipeline Step 12.

Input:
    02_merged_with_SBS/TCGA_merged_expression_SBS.pkl
    01_cleaned_expression/ensg_to_biotype.json
    01_cleaned_expression/ensg_to_symbol.json

Output (-> data/FIG_2/03_differential_expression/{cancer_type}/):
    {ct}_diffexpr_stats.csv              — per-gene Wilcoxon results
    {ct}_selected_genes.csv              — genes passing thresholds
    {ct}_volcano.png                     — volcano plot
    {ct}_manhattan.png                   — Manhattan-style plot
    {ct}_cancer_df2.pkl                  — filtered dataframe for downstream
    {ct}_SBS2_HIGH_group.pkl / .csv      — HIGH SBS2 group (for Step04)
    {ct}_SBS2_LOW_group.pkl / .csv       — LOW SBS2 group (for Step04)
    {ct}_group_definition_summary.txt    — group sizes, thresholds
    {ct}_filtering_summary.txt           — gene filtering log

Usage:
    python Step03_Differential_Expression.py
"""

import os
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from network_config import (
    CANCER_TYPES, CLINICAL_COLS, A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    MIN_SAMPLES_DETECTED,
    A3_SUM_PERCENTILE, SBS2_HIGH_PERCENTILE, SBS2_LOW_PERCENTILE,
    FDR_THRESHOLD, RAW_P_THRESHOLD, LOGFC_THRESHOLD,
    FORCE_KEEP_A3,
    DIR_01_CLEANED, DIR_02_MERGED, DIR_03_DIFFEXPR,
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
# LOAD merged data + mappings
# =============================================================================
banner("[STEP 12] Load merged expression + SBS data")

merged_path = os.path.join(DIR_02_MERGED, "TCGA_merged_expression_SBS.pkl")
log(f"[STEP 12] Reading: {merged_path}")
merged = pd.read_pickle(merged_path)
log(f"[STEP 12] Shape: {merged.shape}")

# Load biotype mapping
biotype_path = os.path.join(DIR_01_CLEANED, "ensg_to_biotype.json")
log(f"[STEP 12] Reading biotype mapping: {biotype_path}")
with open(biotype_path) as f:
    ensg_to_biotype = json.load(f)
log(f"[STEP 12] Biotype mapping entries: {len(ensg_to_biotype)}")

# Load symbol mapping (for labeling)
symbol_path = os.path.join(DIR_01_CLEANED, "ensg_to_symbol.json")
with open(symbol_path) as f:
    ensg_to_symbol = json.load(f)

# Identify all ENSG gene columns
all_gene_cols = [c for c in merged.columns if isinstance(c, str) and c.startswith("ENSG")]
log(f"[STEP 12] Total ENSG gene columns: {len(all_gene_cols)}")


# =============================================================================
# STEP 12.1 — Filter to protein-coding genes only
# =============================================================================
# Aligned with single-cell workflow: SC data starts with ~21K protein-coding
# genes from the reference genome. Here we explicitly filter using the biotype
# annotation from Row 3 of the master FPKM-UQ file.
# =============================================================================
banner("[STEP 12.1] Filter to protein-coding genes")

protein_coding_genes = [g for g in all_gene_cols if ensg_to_biotype.get(g) == "protein_coding"]
non_coding_removed = len(all_gene_cols) - len(protein_coding_genes)

log(f"[STEP 12.1] Protein-coding genes: {len(protein_coding_genes)}")
log(f"[STEP 12.1] Non-coding genes removed: {non_coding_removed}")

# Verify A3 genes survived the filter
a3_in_pc = [g for g in A3_GENES if g in protein_coding_genes]
a3_lost  = [g for g in A3_GENES if g not in protein_coding_genes]
log(f"[STEP 12.1] A3 genes in protein-coding set: {len(a3_in_pc)}")
if a3_lost:
    log(f"[STEP 12.1] WARNING: A3 genes NOT in protein-coding: {a3_lost}")


# =============================================================================
# STEP 12.2 — Filter low-expression genes
# =============================================================================
# Analog of sc.pp.filter_genes(min_cells=20):
# Require gene to have FPKM-UQ > 0 in at least MIN_SAMPLES_DETECTED samples.
# Applied within each cancer type (below), but we can pre-screen globally.
# =============================================================================
banner("[STEP 12.2] Filter low-expression genes (global pre-screen)")

# This is a soft pre-screen — the cancer-type-specific filter below is definitive.
# Here we just log the global numbers for awareness.
global_detected = (merged[protein_coding_genes] > 0).sum(axis=0)
globally_sparse = (global_detected < MIN_SAMPLES_DETECTED).sum()
log(f"[STEP 12.2] Genes detected in <{MIN_SAMPLES_DETECTED} samples (globally): {globally_sparse}")
log(f"[STEP 12.2] (Cancer-specific filtering applied below)")


# =============================================================================
# LOOP over cancer types
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"[CANCER] {cancer_type}", char="=")

    cancer_dir = ensure_dir(os.path.join(DIR_03_DIFFEXPR, cancer_type))

    # ---- Filter to this cancer type
    cancer_df = merged[merged["Project_ID"] == cancer_type].copy()
    log(f"[STEP 12.3] Samples for {cancer_type}: {len(cancer_df)}")

    if len(cancer_df) < 25:
        log("[SKIP] Too few samples (<25)")
        continue

    # ---- Verify required columns
    if "SBS2" not in cancer_df.columns:
        log("[SKIP] SBS2 column missing")
        continue

    for col in ["SBS2", "A3A", "A3B"]:
        cancer_df[col] = pd.to_numeric(cancer_df[col], errors="coerce")
    cancer_df = cancer_df.dropna(subset=["SBS2", "A3A", "A3B"]).reset_index(drop=True)
    log(f"[STEP 12.3] Samples after cleanup: {len(cancer_df)}")

    if len(cancer_df) < 40:
        log("[SKIP] Too few valid samples (<40)")
        continue

    # ---- Cancer-specific low-expression filter on protein-coding genes
    banner(f"[STEP 12.4] Low-expression filter ({cancer_type})")

    pc_expr = cancer_df[protein_coding_genes].apply(pd.to_numeric, errors="coerce")
    n_detected = (pc_expr > 0).sum(axis=0)
    genes_passing = n_detected[n_detected >= MIN_SAMPLES_DETECTED].index.tolist()
    genes_removed = len(protein_coding_genes) - len(genes_passing)

    log(f"[STEP 12.4] Protein-coding genes input: {len(protein_coding_genes)}")
    log(f"[STEP 12.4] Genes detected in >={MIN_SAMPLES_DETECTED} samples: {len(genes_passing)}")
    log(f"[STEP 12.4] Low-expression genes removed: {genes_removed}")

    # Ensure A3 genes are retained even if low-expression
    if FORCE_KEEP_A3:
        for g in A3_GENES:
            if g in protein_coding_genes and g not in genes_passing:
                genes_passing.append(g)
                log(f"[STEP 12.4] Force-kept low-expression A3 gene: {g}")

    gene_universe = list(dict.fromkeys(genes_passing))  # deduplicate
    log(f"[STEP 12.4] Final gene universe for DEG: {len(gene_universe)}")

    # ---- Save filtering summary
    filter_summary_path = os.path.join(cancer_dir, f"{cancer_type}_filtering_summary.txt")
    with open(filter_summary_path, "w") as f:
        f.write(f"Gene Filtering Summary — {cancer_type}\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total ENSG columns in data: {len(all_gene_cols)}\n")
        f.write(f"After protein-coding filter: {len(protein_coding_genes)}\n")
        f.write(f"After low-expression filter (>={MIN_SAMPLES_DETECTED} samples): {len(genes_passing)}\n")
        f.write(f"Final gene universe: {len(gene_universe)}\n")

    # =========================================================================
    # STEP 12.5 — Define A3-controlled HIGH/LOW SBS2 groups
    # =========================================================================
    # Logic:
    #   1. Sum raw A3A + A3B FPKM-UQ per tumor
    #   2. Keep only tumors ABOVE median of this sum
    #      → ensures both groups have similar A3 expression
    #   3. Within high-A3 tumors, rank by SBS2
    #   4. Top 25% = SBS2-HIGH, bottom 25% = SBS2-LOW
    #      → equal group sizes, controlled for A3
    #
    # This is the SINGLE SOURCE OF TRUTH for group definitions.
    # Step04 inherits these groups directly.
    # =========================================================================
    banner(f"[STEP 12.5] Define A3-controlled groups ({cancer_type})")

    # Raw sum of A3A + A3B (not z-scored)
    cancer_df["A3_sum"] = cancer_df["A3A"] + cancer_df["A3B"]

    a3_median = float(cancer_df["A3_sum"].quantile(A3_SUM_PERCENTILE))
    high_a3_df = cancer_df[cancer_df["A3_sum"] >= a3_median].copy()

    log(f"[STEP 12.5] A3A + A3B sum — median threshold: {a3_median:.4f}")
    log(f"[STEP 12.5] Tumors above median (high-A3): {len(high_a3_df)} / {len(cancer_df)}")

    if len(high_a3_df) < 20:
        log("[SKIP] Too few high-A3 tumors (<20)")
        continue

    # Within high-A3: rank by SBS2, take bottom 25% first, then match from top
    # This guarantees exactly equal group sizes
    high_a3_ranked = high_a3_df.sort_values("SBS2", ascending=True).reset_index(drop=True)
    n_high_a3 = len(high_a3_ranked)
    n_per_group = int(np.floor(n_high_a3 * 0.25))

    # Bottom 25% (lowest SBS2) — take first n_per_group from ranked list
    group_low = high_a3_ranked.iloc[:n_per_group].copy()

    # Top 25% (highest SBS2) — take last n_per_group from ranked list
    group_high = high_a3_ranked.iloc[-n_per_group:].copy()

    sbs2_low_max  = float(group_low["SBS2"].max())    # upper boundary of LOW group
    sbs2_high_min = float(group_high["SBS2"].min())    # lower boundary of HIGH group

    log(f"[STEP 12.5] High-A3 tumors ranked by SBS2: {n_high_a3}")
    log(f"[STEP 12.5] Tumors per group (25%): {n_per_group}")
    log(f"[STEP 12.5] SBS2-LOW group:  n={len(group_low)}, SBS2 range: "
        f"{group_low['SBS2'].min():.6f} — {sbs2_low_max:.6f}")
    log(f"[STEP 12.5] SBS2-HIGH group: n={len(group_high)}, SBS2 range: "
        f"{sbs2_high_min:.6f} — {group_high['SBS2'].max():.6f}")

    # Verify A3 expression is similar between groups
    a3_mean_high = group_high["A3_sum"].mean()
    a3_mean_low  = group_low["A3_sum"].mean()
    log(f"[STEP 12.5] Mean A3_sum — HIGH: {a3_mean_high:.4f}, LOW: {a3_mean_low:.4f}")

    if len(group_high) < 8 or len(group_low) < 8:
        log("[SKIP] Groups too small (<8)")
        continue

    # Save group definitions (Step04 inherits these)
    group_high.to_pickle(os.path.join(cancer_dir, f"{cancer_type}_SBS2_HIGH_group.pkl"))
    group_low.to_pickle(os.path.join(cancer_dir, f"{cancer_type}_SBS2_LOW_group.pkl"))
    group_high.to_csv(os.path.join(cancer_dir, f"{cancer_type}_SBS2_HIGH_group.csv"), index=False)
    group_low.to_csv(os.path.join(cancer_dir, f"{cancer_type}_SBS2_LOW_group.csv"), index=False)
    log(f"[SAVE] Group definitions saved (for Step04 inheritance)")

    # Group summary
    group_summary_path = os.path.join(cancer_dir, f"{cancer_type}_group_definition_summary.txt")
    with open(group_summary_path, "w") as f:
        f.write(f"Group Definition Summary — {cancer_type}\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total tumors: {len(cancer_df)}\n")
        f.write(f"A3A + A3B sum median threshold: {a3_median:.4f}\n")
        f.write(f"High-A3 tumors (above median): {len(high_a3_df)}\n")
        f.write(f"Tumors per group (25% of high-A3): {n_per_group}\n")
        f.write(f"SBS2-LOW max: {sbs2_low_max:.6f}\n")
        f.write(f"SBS2-HIGH min: {sbs2_high_min:.6f}\n")
        f.write(f"SBS2-HIGH group size: {len(group_high)}\n")
        f.write(f"SBS2-LOW group size: {len(group_low)}\n")
        f.write(f"Mean A3_sum HIGH: {a3_mean_high:.4f}\n")
        f.write(f"Mean A3_sum LOW: {a3_mean_low:.4f}\n")

    # =========================================================================
    # STEP 12.6 — Log1p transform + Wilcoxon rank-sum test
    # =========================================================================
    # Aligned with single-cell: sc.pp.log1p() then sc.tl.rank_genes_groups
    # with method='wilcoxon'.
    #
    # FPKM-UQ is already library-size normalized (analogous to CPM), so we
    # skip the CPM step and go straight to log1p. The Wilcoxon rank-sum test
    # is identical to scipy.stats.mannwhitneyu (one-sided → we use two-sided).
    #
    # log2FC is computed as: mean(log1p(HIGH)) - mean(log1p(LOW))
    # This matches scanpy's logfoldchanges output.
    # =========================================================================
    banner(f"[STEP 12.6] Wilcoxon rank-sum test ({cancer_type})")

    results = []
    n_tested = 0

    for gene in gene_universe:
        vals_high = pd.to_numeric(group_high[gene], errors="coerce").dropna().values
        vals_low  = pd.to_numeric(group_low[gene], errors="coerce").dropna().values

        if len(vals_high) < 3 or len(vals_low) < 3:
            continue

        # Log1p transform (matching sc.pp.log1p)
        log_high = np.log1p(vals_high)
        log_low  = np.log1p(vals_low)

        # Wilcoxon rank-sum test (two-sided) — scipy's mannwhitneyu
        try:
            stat, p_val = mannwhitneyu(log_high, log_low, alternative="two-sided")
        except ValueError:
            # All values identical between groups
            continue

        # Log fold change from log-transformed means
        # (matches scanpy's logfoldchanges = mean_group - mean_rest on log1p data)
        mean_log_high = np.mean(log_high)
        mean_log_low  = np.mean(log_low)
        logfc = mean_log_high - mean_log_low

        # Also store raw means for reference
        mean_raw_high = np.mean(vals_high)
        mean_raw_low  = np.mean(vals_low)

        results.append({
            "gene": gene,
            "gene_symbol": ensg_to_symbol.get(gene, gene),
            "U_stat": stat,
            "p_value": p_val,
            "logFC": logfc,
            "mean_log1p_HIGH": mean_log_high,
            "mean_log1p_LOW": mean_log_low,
            "mean_raw_HIGH": mean_raw_high,
            "mean_raw_LOW": mean_raw_low,
        })
        n_tested += 1

    stats_df = pd.DataFrame(results)
    log(f"[STEP 12.6] Genes tested: {n_tested}")

    # ---- FDR correction (computed for reference, not used for selection)
    stats_df["p_adj"] = bh_fdr(stats_df["p_value"].values)
    stats_df = stats_df.sort_values("p_value").reset_index(drop=True)

    # Report both raw and FDR counts for transparency
    n_sig_raw = (stats_df["p_value"] < RAW_P_THRESHOLD).sum()
    n_sig_fdr = (stats_df["p_adj"] < FDR_THRESHOLD).sum()
    n_sig_sel = ((stats_df["p_value"] < RAW_P_THRESHOLD) & (stats_df["logFC"].abs() > LOGFC_THRESHOLD)).sum()
    n_up      = ((stats_df["p_value"] < RAW_P_THRESHOLD) & (stats_df["logFC"] > LOGFC_THRESHOLD)).sum()
    n_down    = ((stats_df["p_value"] < RAW_P_THRESHOLD) & (stats_df["logFC"] < -LOGFC_THRESHOLD)).sum()

    log(f"[STEP 12.7] Significant (raw p < {RAW_P_THRESHOLD}): {n_sig_raw}")
    log(f"[STEP 12.7] Significant (FDR < {FDR_THRESHOLD}): {n_sig_fdr} (reference only)")
    log(f"[STEP 12.7] Selected (raw p < {RAW_P_THRESHOLD} AND |logFC| > {LOGFC_THRESHOLD}): {n_sig_sel}")
    log(f"[STEP 12.7] Upregulated in HIGH SBS2: {n_up}")
    log(f"[STEP 12.7] Downregulated in HIGH SBS2: {n_down}")

    # ---- Select genes (raw p-value + logFC thresholds)
    sig_mask = (stats_df["p_value"] < RAW_P_THRESHOLD) & (stats_df["logFC"].abs() > LOGFC_THRESHOLD)
    selected = stats_df.loc[sig_mask, "gene"].tolist()

    # Force-keep A3 genes
    if FORCE_KEEP_A3:
        for g in A3_GENES:
            if g in stats_df["gene"].values and g not in selected:
                selected.append(g)
                sym = ensg_to_symbol.get(g, g)
                row = stats_df[stats_df["gene"] == g].iloc[0]
                log(f"[STEP 12.8] Force-kept A3 gene: {sym} "
                    f"(p={row['p_value']:.4f}, logFC={row['logFC']:.3f})")

    selected = list(dict.fromkeys(selected))  # deduplicate, preserve order
    log(f"[STEP 12.8] Total selected genes: {len(selected)}")

    # Report biomarker status
    tested_set = set(stats_df["gene"])
    log("[STEP 12.8] Biomarkers tested/selected:")
    for g in BIOMARKERS:
        tested = "TESTED" if g in tested_set else "NOT_TESTED"
        sel = "SELECTED" if g in selected else "NOT_SELECTED"
        sym = ensg_to_symbol.get(g, g)
        log(f"  {sym} ({g}): {tested} / {sel}")

    # ---- Save stats
    stats_csv = os.path.join(cancer_dir, f"{cancer_type}_diffexpr_stats.csv")
    stats_df.to_csv(stats_csv, index=False)
    log(f"[SAVE] Stats -> {stats_csv}")

    # ---- Save selected gene list
    sel_csv = os.path.join(cancer_dir, f"{cancer_type}_selected_genes.csv")
    sel_df = stats_df[stats_df["gene"].isin(selected)][["gene", "gene_symbol", "logFC", "p_adj"]].copy()
    sel_df.to_csv(sel_csv, index=False)
    log(f"[SAVE] Selected genes -> {sel_csv}")

    # =========================================================================
    # STEP 12.8B — Filter selected genes for network readiness
    # =========================================================================
    # Ensures selected genes are usable in correlation networks (Step 04+):
    #   1. Gene must be present as a column in both HIGH and LOW groups
    #   2. Gene must not be constant (nunique <= 1) in BOTH groups
    # This was previously done in a separate Step04 but now lives here
    # since Step03 is the single source of truth for groups + gene selection.
    # =========================================================================
    banner(f"[STEP 12.8B] Filter genes for network readiness ({cancer_type})")

    # Filter to genes present in both group DataFrames
    selected_filtered = [g for g in selected
                         if g in group_high.columns and g in group_low.columns]
    n_missing = len(selected) - len(selected_filtered)
    if n_missing > 0:
        log(f"[STEP 12.8B] Removed {n_missing} genes not present in both groups")

    # Remove genes constant in BOTH groups (would break correlation computation)
    high_nunique = group_high[selected_filtered].apply(pd.to_numeric, errors="coerce").nunique()
    low_nunique  = group_low[selected_filtered].apply(pd.to_numeric, errors="coerce").nunique()
    const_both = (high_nunique <= 1) & (low_nunique <= 1)
    n_const = int(const_both.sum())
    if n_const > 0:
        log(f"[STEP 12.8B] Removing {n_const} genes constant in both groups")
        selected_filtered = [g for g in selected_filtered if not const_both.get(g, False)]

    selected_filtered = list(dict.fromkeys(selected_filtered))  # deduplicate
    log(f"[STEP 12.8B] Network-ready genes: {len(selected_filtered)}")

    # Save filtered gene list (this is what Step04/correlation networks will use)
    filtered_csv = os.path.join(cancer_dir, f"{cancer_type}_selected_genes_filtered.csv")
    pd.DataFrame({"gene": selected_filtered}).to_csv(filtered_csv, index=False)
    log(f"[SAVE] Filtered genes -> {filtered_csv}")

    # =========================================================================
    # STEP 12.9 — Volcano plot
    # =========================================================================
    banner(f"[STEP 12.9] Volcano plot ({cancer_type})")

    fig, ax = plt.subplots(figsize=(10, 7))
    neg_log_p = -np.log10(stats_df["p_value"].clip(lower=1e-300))

    # Non-significant (gray)
    ns_mask = ~sig_mask
    ax.scatter(stats_df.loc[ns_mask, "logFC"], neg_log_p[ns_mask],
               c="gray", alpha=0.3, s=8, label="NS", zorder=1)

    # Significant UP (red)
    up_mask = sig_mask & (stats_df["logFC"] > 0)
    ax.scatter(stats_df.loc[up_mask, "logFC"], neg_log_p[up_mask],
               c="firebrick", alpha=0.6, s=14, label=f"Up (n={up_mask.sum()})", zorder=2)

    # Significant DOWN (blue)
    down_mask = sig_mask & (stats_df["logFC"] < 0)
    ax.scatter(stats_df.loc[down_mask, "logFC"], neg_log_p[down_mask],
               c="steelblue", alpha=0.6, s=14, label=f"Down (n={down_mask.sum()})", zorder=2)

    # Threshold lines
    ax.axhline(-np.log10(RAW_P_THRESHOLD), ls="--", c="black", lw=0.8, alpha=0.5)
    ax.axvline(LOGFC_THRESHOLD, ls="--", c="black", lw=0.8, alpha=0.5)
    ax.axvline(-LOGFC_THRESHOLD, ls="--", c="black", lw=0.8, alpha=0.5)

    # Label A3 genes
    for g in A3_GENES:
        if g in stats_df["gene"].values:
            row = stats_df[stats_df["gene"] == g].iloc[0]
            idx = stats_df[stats_df["gene"] == g].index[0]
            sym = ensg_to_symbol.get(g, g)
            ax.annotate(sym, (row["logFC"], neg_log_p[idx]),
                        fontsize=7, fontweight="bold", color="darkgreen",
                        xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("log2 Fold Change (HIGH SBS2 / LOW SBS2)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{cancer_type} | Volcano — Wilcoxon (A3-controlled groups)")
    ax.legend(fontsize=9, loc="upper right")
    plt.tight_layout()

    volcano_path = os.path.join(cancer_dir, f"{cancer_type}_volcano.png")
    plt.savefig(volcano_path, dpi=300)
    plt.close()
    log(f"[SAVE] Volcano -> {volcano_path}")

    # =========================================================================
    # STEP 12.10 — Manhattan-style plot
    # =========================================================================
    banner(f"[STEP 12.10] Manhattan plot ({cancer_type})")

    sorted_stats = stats_df.sort_values("p_value").reset_index(drop=True)
    y = -np.log10(sorted_stats["p_value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(range(len(y)), y, s=4, c="steelblue", alpha=0.5)
    ax.axhline(-np.log10(RAW_P_THRESHOLD), ls="--", c="red", lw=0.8)
    ax.set_xlabel("Genes (sorted by p-value)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{cancer_type} | Manhattan — Wilcoxon (A3-controlled groups)")
    plt.tight_layout()

    manhattan_path = os.path.join(cancer_dir, f"{cancer_type}_manhattan.png")
    plt.savefig(manhattan_path, dpi=300)
    plt.close()
    log(f"[SAVE] Manhattan -> {manhattan_path}")

    # =========================================================================
    # STEP 12.11 — Selection plot (SBS2 vs A3 sum with groups highlighted)
    # =========================================================================
    banner(f"[STEP 12.11] Selection plot ({cancer_type})")

    fig, ax = plt.subplots(figsize=(9, 7))

    # All tumors (gray)
    ax.scatter(cancer_df["A3_sum"], cancer_df["SBS2"],
               c="lightgray", s=15, alpha=0.5, label="All tumors", zorder=1)

    # A3 threshold
    ax.axvline(a3_median, ls="--", c="black", lw=0.8, alpha=0.5, label="A3 median")

    # SBS2 thresholds
    ax.axhline(sbs2_high_min, ls=":", c="firebrick", lw=0.8, alpha=0.6)
    ax.axhline(sbs2_low_max, ls=":", c="steelblue", lw=0.8, alpha=0.6)

    # HIGH group (red)
    ax.scatter(group_high["A3_sum"], group_high["SBS2"],
               c="firebrick", s=30, alpha=0.8, edgecolors="black", linewidths=0.3,
               label=f"SBS2-HIGH (n={len(group_high)})", zorder=3)

    # LOW group (blue)
    ax.scatter(group_low["A3_sum"], group_low["SBS2"],
               c="steelblue", s=30, alpha=0.8, edgecolors="black", linewidths=0.3,
               label=f"SBS2-LOW (n={len(group_low)})", zorder=3)

    ax.set_xlabel("A3A + A3B (raw FPKM-UQ sum)", fontsize=12)
    ax.set_ylabel("SBS2 Weight", fontsize=12)
    ax.set_title(f"{cancer_type} | SBS2 vs A3 Sum — Group Selection", fontsize=13)
    ax.legend(fontsize=10, loc="upper left")
    plt.tight_layout()

    selection_path = os.path.join(cancer_dir, f"{cancer_type}_SBS2_vs_A3sum_selection.png")
    plt.savefig(selection_path, dpi=300)
    plt.close()
    log(f"[SAVE] Selection plot -> {selection_path}")

    # =========================================================================
    # STEP 12.12 — Build cancer_df2 (clinical + selected genes + signal columns)
    # =========================================================================
    banner(f"[STEP 12.12] Build cancer_df2 ({cancer_type})")

    keep_cols = list(CLINICAL_COLS)

    # Add SBS2, A3 alias columns, and A3_sum
    for c in ["SBS2", "A3A", "A3B", "A3C", "A3H", "A3_sum"]:
        if c in cancer_df.columns and c not in keep_cols:
            keep_cols.append(c)

    # Add selected genes
    keep_cols += [g for g in selected if g in cancer_df.columns and g not in keep_cols]

    cancer_df2 = cancer_df[keep_cols].copy()
    log(f"[STEP 12.12] cancer_df2 shape: {cancer_df2.shape}")

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

    # ---- Console summary
    log(f"\n[STEP 12 COMPLETE for {cancer_type}]")
    print(f"  Gene universe (protein-coding, expressed): {len(gene_universe)}")
    print(f"  SBS2-HIGH group: {len(group_high)} | SBS2-LOW group: {len(group_low)}")
    print(f"  Genes tested: {n_tested}")
    print(f"  Significant (raw p + logFC): {n_sig_sel}")
    print(f"  Selected (incl. force-kept A3): {len(selected)}")

banner("[STEP 03 COMPLETE — ALL CANCER TYPES]")
