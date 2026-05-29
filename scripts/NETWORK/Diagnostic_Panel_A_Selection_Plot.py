#!/usr/bin/env python3
"""
Diagnostic_Panel_A_Selection_Plot.py
=====================================
Standalone script to regenerate the SBS2 vs A3A+A3B selection plot
at 1/3 original width for tight figure layout.

Everything identical to the supplement selection plot in
Generate_Figure2_Panels.py except figsize width compressed 3x.

Usage:
    conda run -n NETWORK python Diagnostic_Panel_A_Selection_Plot.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

from network_config import (
    CANCER_TYPES, DIR_02_MERGED, DIR_03_DIFFEXPR,
    FIG2_ROOT, banner, log, ensure_dir
)

# =============================================================================
# STYLE (identical to Generate_Figure2_Panels.py)
# =============================================================================
COLOR_HIGH  = "#ed6a5a"
COLOR_LOW   = "#5b8e7d"
COLOR_CREAM = "#f4f1bb"
COLOR_CORAL = "#ed6a5a"
COLOR_GRAY  = "#d0d0d0"

FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22

DPI = 300

FIG_DIR = ensure_dir(os.path.join(FIG2_ROOT, "FIGURE_2_PANELS"))

# =============================================================================
# MAIN
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"Selection Plot (compressed) -- {cancer_type}")
    ct_fig_dir = ensure_dir(os.path.join(FIG_DIR, cancer_type))

    # ---- Load data ----
    merged = pd.read_pickle(os.path.join(DIR_02_MERGED,
                            "TCGA_merged_expression_SBS.pkl"))
    all_hnsc = merged[merged["Project_ID"] == cancer_type].copy()
    for c in ["SBS2", "A3A", "A3B"]:
        all_hnsc[c] = pd.to_numeric(all_hnsc[c], errors="coerce")
    all_hnsc = all_hnsc.dropna(subset=["SBS2", "A3A", "A3B"])
    all_hnsc["A3_sum"] = all_hnsc["A3A"] + all_hnsc["A3B"]
    log(f"All HNSC tumors: {len(all_hnsc)}")

    high_path = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                             f"{cancer_type}_SBS2_HIGH_group.pkl")
    low_path  = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                             f"{cancer_type}_SBS2_LOW_group.pkl")
    if not os.path.exists(high_path):
        log(f"[SKIP] Groups not found"); continue
    group_high = pd.read_pickle(high_path)
    group_low  = pd.read_pickle(low_path)
    log(f"HIGH: {len(group_high)}, LOW: {len(group_low)}")

    # ---- Thresholds ----
    a3_median = float(all_hnsc["A3_sum"].median())
    sbs2_low_max = float(group_low["SBS2"].max())
    sbs2_high_min = float(group_high["SBS2"].min())
    x_max_data = all_hnsc["A3_sum"].max()
    y_max_data = all_hnsc["SBS2"].max()

    # ---- Outlier detection (broken axis logic, unchanged) ----
    a3_sorted = all_hnsc["A3_sum"].sort_values(ascending=False).values
    has_outliers = False
    if len(a3_sorted) >= 3:
        best_gap_ratio = 1.0
        best_gap_idx = -1
        for i in range(min(10, len(a3_sorted) - 1)):
            ratio = (a3_sorted[i] / a3_sorted[i + 1]
                     if a3_sorted[i + 1] > 0 else 1.0)
            if ratio > best_gap_ratio:
                best_gap_ratio = ratio
                best_gap_idx = i
        if best_gap_ratio > 1.5 and best_gap_idx >= 0:
            cluster_max = a3_sorted[best_gap_idx + 1]
            outlier_min = a3_sorted[best_gap_idx]
            x_break_start = cluster_max * 1.10
            x_break_end = outlier_min * 0.90
            has_outliers = True

    # ---- Figure (1/3 width) ----
    # Original: (16, 10) with outliers, (14, 10) without
    # Compressed: (5.33, 10) with outliers, (4.67, 10) without
    if has_outliers:
        fig, (ax_main, ax_break) = plt.subplots(
            1, 2, sharey=True, figsize=(5.33, 10),
            gridspec_kw={"width_ratios": [3, 1], "wspace": 0.03})
        axes_list = [ax_main, ax_break]
    else:
        fig, ax_main = plt.subplots(figsize=(4.67, 10))
        ax_break = None
        axes_list = [ax_main]

    for ax in axes_list:
        y_pad = y_max_data * 0.05
        y_max = y_max_data + y_pad
        if ax == ax_main:
            x_lo = -2
            x_hi = x_break_start if has_outliers else x_max_data * 1.05
        elif ax == ax_break:
            x_lo, x_hi = x_break_end, x_max_data * 1.05

        cream_poly = MplPolygon(
            [(a3_median, 0), (a3_median, sbs2_low_max),
             (x_hi, sbs2_low_max), (x_hi, 0)],
            closed=True, facecolor=COLOR_CREAM, alpha=0.75,
            edgecolor="none", zorder=0)
        ax.add_patch(cream_poly)
        coral_poly = MplPolygon(
            [(a3_median, sbs2_high_min), (a3_median, y_max),
             (x_hi, y_max), (x_hi, sbs2_high_min)],
            closed=True, facecolor=COLOR_CORAL, alpha=0.55,
            edgecolor="none", zorder=0)
        ax.add_patch(coral_poly)
        if sbs2_high_min > sbs2_low_max:
            mid_poly = MplPolygon(
                [(a3_median, sbs2_low_max), (a3_median, sbs2_high_min),
                 (x_hi, sbs2_high_min), (x_hi, sbs2_low_max)],
                closed=True, facecolor="#e8e4c8", alpha=0.55,
                edgecolor="none", zorder=0)
            ax.add_patch(mid_poly)

        ax.axvline(a3_median, ls="--", c="#4D4D4D", lw=0.8, zorder=1)
        ax.axhline(sbs2_high_min, ls=":", c="#4D4D4D", lw=0.8, alpha=0.6)
        ax.axhline(sbs2_low_max, ls=":", c="#4D4D4D", lw=0.8, alpha=0.6)

        ax.scatter(all_hnsc["A3_sum"], all_hnsc["SBS2"],
                   facecolors=COLOR_GRAY, edgecolors="#000000",
                   linewidths=0.3, s=100, alpha=0.5, zorder=2)
        ax.scatter(group_high["A3_sum"], group_high["SBS2"],
                   facecolors=COLOR_HIGH, edgecolors="#000000",
                   linewidths=0.8, s=100, alpha=0.85, zorder=3,
                   label=f"SBS2-HIGH (n={len(group_high)})"
                   if ax == ax_main else "")
        ax.scatter(group_low["A3_sum"], group_low["SBS2"],
                   facecolors=COLOR_LOW, edgecolors="#000000",
                   linewidths=0.8, s=100, alpha=0.85, zorder=3,
                   label=f"SBS2-LOW (n={len(group_low)})"
                   if ax == ax_main else "")

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(-y_pad * 0.5, y_max)
        ax.tick_params(axis="both", labelsize=FONT_TICK)

    ax_main.set_xlabel("A3A + A3B Expression (FPKM-UQ)",
                       fontsize=FONT_AXIS, labelpad=12)
    ax_main.set_ylabel("SBS2 Mutational Weight",
                       fontsize=FONT_AXIS, labelpad=12)
    ax_main.legend(fontsize=FONT_LEGEND - 4, loc="upper left",
                   framealpha=0.9)

    if ax_break is not None:
        ax_break.tick_params(axis="y", left=False, labelleft=False)
        d = 0.015
        for ax_b, side in [(ax_main, "right"), (ax_break, "left")]:
            kwargs = dict(transform=ax_b.transAxes, color="#000000",
                          clip_on=False, lw=1.5)
            if side == "right":
                ax_b.plot((1-d, 1+d), (-d, +d), **kwargs)
                ax_b.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
                ax_b.spines["right"].set_visible(False)
            else:
                ax_b.plot((-d, +d), (-d, +d), **kwargs)
                ax_b.plot((-d, +d), (1-d, 1+d), **kwargs)
                ax_b.spines["left"].set_visible(False)

    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(ct_fig_dir,
                    f"{cancer_type}_Panel_A_selection_compressed.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close()
    log(f"[SAVE] Compressed selection plot -> {ct_fig_dir}")

banner("DONE")
