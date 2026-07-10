#!/usr/bin/env python3
"""
Generate_Selection_Robustness_Figures.py
=========================================

Polished supplemental figures for the selection-robustness diagnostic.

Reads the diagnostic's on-disk output (simple_group_membership.tsv) plus adata,
and produces three standalone figures supporting the supplemental-methods
argument that the Step00B composite selection is robust:

  1. UMAP_selection_comparison       three UMAPs (canonical / 2-term / 3-term),
                                      SBS2 arm coral, CNV arm mustard, over a
                                      faint all-cell background. Simple-method
                                      panels carry their recovery % vs canonical.
  2. Metric_distributions            2x2 boxplots of SBS2 / CNV / A3A fraction /
                                      A3B fraction across all six groups, split
                                      into SBS2 arm (coral) and CNV arm (mustard).
  3. Dominance_split                 stacked bars, % A3A-dominant vs % A3B-dominant
                                      among A3-expressing cells, per group. The
                                      two-term groups saw no A3 input, so their
                                      split is intrinsic.

The figures are downstream of the diagnostic: nothing is re-selected here. Per-cell
metrics are recomputed by the SAME code path Step00B / the diagnostic use, so the
values plotted are identical to those the diagnostic reported.

Style: 28-34 pt text (ticks 26 for the crowded six-group axes), hex colors,
PDF + PNG at 300 DPI.

Usage:
    conda run -n NETWORK python Generate_Selection_Robustness_Figures.py
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Find network_config_SC.py up the tree (script may live in TROUBLESHOOTING/)
_here = os.path.dirname(os.path.abspath(__file__))
_search = _here
for _ in range(5):
    if os.path.exists(os.path.join(_search, "network_config_SC.py")):
        if _search not in sys.path:
            sys.path.insert(0, _search)
        break
    _parent = os.path.dirname(_search)
    if _parent == _search:
        break
    _search = _parent
else:
    raise ModuleNotFoundError(
        f"network_config_SC.py not found in any parent of {_here}.")

from network_config_SC import (
    DIR_01_GROUPS, ADATA_FINAL_PATH, WEIGHTS_PATH,
    TARGET_CELL_TYPE, banner, log, ensure_dir
)

# =============================================================================
# CONFIGURATION
# =============================================================================

DIAG_DIR = os.path.join(DIR_01_GROUPS, "selection_robustness_diagnostic")
MEMBERSHIP_FILE = os.path.join(DIAG_DIR, "simple_group_membership.tsv")
OUT_DIR = ensure_dir(os.path.join(DIAG_DIR, "FIGURES"))

N_CELLS = 546

# Palette (paper-consistent)
COLOR_SBS2 = "#ed6a5a"   # coral  -> SBS2 arm / A3A-dominant
COLOR_CNV  = "#F6D155"   # mustard-> CNV arm  / A3B-dominant
COLOR_BG   = "#f0f0f0"   # faint background

FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 28
FONT_TICK   = 26

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
})

# Column name in membership TSV -> friendly label / arm
GROUP_SPEC = [
    ("canonical_SBS2", "Canonical\nSBS2", "sbs2"),
    ("s2_SBS2",        "2-term\nSBS2",    "sbs2"),
    ("s3_SBS2",        "3-term\nSBS2",    "sbs2"),
    ("canonical_CNV",  "Canonical\nCNV",  "cnv"),
    ("s2_CNV",         "2-term\nCNV",     "cnv"),
    ("s3_CNV",         "3-term\nCNV",     "cnv"),
]


# =============================================================================
# LOAD
# =============================================================================

def load_metrics():
    """Load adata; recompute per-cell metrics by the Step00B/diagnostic path."""
    banner("Load adata and recompute per-cell metrics")

    log(f"  Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  Total cells: {adata.n_obs:,}")

    basal = adata[adata.obs["final_annotation"] == TARGET_CELL_TYPE].copy()
    log(f"  Basal cells: {basal.n_obs:,}")

    weights_df = pd.read_csv(WEIGHTS_PATH, sep="\t", index_col=0)
    if "SBS2" in weights_df.index:
        sbs2_map = weights_df.loc["SBS2"].to_dict()
        basal.obs["SBS2"] = basal.obs_names.map(lambda x: sbs2_map.get(x, 0.0)).astype(float)

    gene_names = list(basal.var_names)
    for gene in ["APOBEC3A", "APOBEC3B"]:
        if gene in gene_names:
            vals = basal.X[:, gene_names.index(gene)]
            vals = vals.toarray().flatten() if hasattr(vals, "toarray") else np.array(vals).flatten()
            basal.obs[gene] = vals
        else:
            basal.obs[gene] = 0.0
    basal.obs["A3_sum"] = basal.obs["APOBEC3A"] + basal.obs["APOBEC3B"]
    basal.obs["A3A_fraction"] = basal.obs["APOBEC3A"] / (basal.obs["APOBEC3A"] + basal.obs["APOBEC3B"] + 0.01)
    basal.obs["A3B_fraction"] = basal.obs["APOBEC3B"] / (basal.obs["APOBEC3A"] + basal.obs["APOBEC3B"] + 0.01)
    for col in ["cnv_score", "CytoTRACE2_Score"]:
        if col in basal.obs.columns:
            basal.obs[col] = basal.obs[col].astype(float)

    m = basal.obs[["SBS2", "cnv_score", "A3A_fraction", "A3B_fraction", "A3_sum"]].copy()
    return adata, m


def load_membership():
    """Read the diagnostic membership TSV into per-group barcode sets."""
    if not os.path.exists(MEMBERSHIP_FILE):
        raise FileNotFoundError(
            f"Membership file not found: {MEMBERSHIP_FILE}\n"
            f"  Run Diagnostic_Selection_Robustness.py first.")
    df = pd.read_csv(MEMBERSHIP_FILE, sep="\t", index_col=0)
    sets = {}
    for col, _, _ in GROUP_SPEC:
        flag = df[col].astype(str).str.strip().str.lower().isin(["true", "1"])
        sets[col] = set(df.index[flag])
        log(f"  {col}: {len(sets[col]):,} cells")
    return sets


def pct_overlap(a, b):
    a, b = set(a), set(b)
    return 100.0 * len(a & b) / len(a) if a else 0.0


# =============================================================================
# FIGURE 1: UMAP COMPARISON
# =============================================================================

def fig_umap(adata, sets):
    banner("FIGURE 1: UMAP selection comparison")
    umap = adata.obsm["X_umap"]
    obs_names = adata.obs_names

    methods = [
        ("Canonical composite", sets["canonical_SBS2"], sets["canonical_CNV"], None),
        ("Two-term (SBS2 + CNV)", sets["s2_SBS2"], sets["s2_CNV"],
         (pct_overlap(sets["s2_SBS2"], sets["canonical_SBS2"]),
          pct_overlap(sets["s2_CNV"], sets["canonical_CNV"]))),
        ("Three-term (+ A3 dominance)", sets["s3_SBS2"], sets["s3_CNV"],
         (pct_overlap(sets["s3_SBS2"], sets["canonical_SBS2"]),
          pct_overlap(sets["s3_CNV"], sets["canonical_CNV"]))),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(39, 13))
    for ax, (title, sbs2_set, cnv_set, recov) in zip(axes, methods):
        ax.scatter(umap[:, 0], umap[:, 1], s=3, c=COLOR_BG, alpha=0.3,
                   edgecolors="none", rasterized=True, zorder=0)
        c_mask = obs_names.isin(cnv_set)
        ax.scatter(umap[c_mask, 0], umap[c_mask, 1], s=22, c=COLOR_CNV, alpha=0.75,
                   edgecolors="#000000", linewidths=0.2, rasterized=True, zorder=2)
        s_mask = obs_names.isin(sbs2_set)
        ax.scatter(umap[s_mask, 0], umap[s_mask, 1], s=22, c=COLOR_SBS2, alpha=0.75,
                   edgecolors="#000000", linewidths=0.2, rasterized=True, zorder=3)

        t = title
        if recov is not None:
            t += f"\nrecovered: SBS2 {recov[0]:.0f}%  /  CNV {recov[1]:.0f}%"
        ax.set_title(t, fontsize=FONT_TITLE, pad=12)
        ax.set_xlabel("UMAP 1", fontsize=FONT_AXIS)
        if ax is axes[0]:
            ax.set_ylabel("UMAP 2", fontsize=FONT_AXIS)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    handles = [
        Patch(fc=COLOR_SBS2, ec="#000000", label="SBS2-HIGH (n=546)"),
        Patch(fc=COLOR_CNV, ec="#000000", label="CNV-HIGH (n=546)"),
        Patch(fc=COLOR_BG, ec="#999999", label="Other cells"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=FONT_LEGEND,
               framealpha=0.9, bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout(rect=(0, 0.04, 1, 1))
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(OUT_DIR, f"UMAP_selection_comparison.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log("  [SAVE] UMAP_selection_comparison")


# =============================================================================
# FIGURE 2: METRIC DISTRIBUTIONS
# =============================================================================

def fig_distributions(m, sets):
    banner("FIGURE 2: Metric distributions")
    metrics = [("SBS2", "SBS2 weight"), ("cnv_score", "CNV score"),
               ("A3A_fraction", "A3A fraction"), ("A3B_fraction", "A3B fraction")]

    fig, axes = plt.subplots(2, 2, figsize=(24, 20))
    for ax, (col, title) in zip(axes.flat, metrics):
        data, colors = [], []
        for key, _, arm in GROUP_SPEC:
            vals = m.loc[m.index.isin(sets[key]), col].dropna().values
            data.append(vals)
            colors.append(COLOR_SBS2 if arm == "sbs2" else COLOR_CNV)

        bp = ax.boxplot(data, patch_artist=True, showmeans=True, showfliers=False,
                        widths=0.6,
                        medianprops=dict(color="#000000", lw=2.5),
                        meanprops=dict(marker="D", markerfacecolor="#000000",
                                       markeredgecolor="#000000", markersize=9),
                        whiskerprops=dict(color="#555555", lw=1.5),
                        capprops=dict(color="#555555", lw=1.5))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.85); patch.set_edgecolor("#000000")

        ax.axvline(3.5, color="#999999", ls=":", lw=2)
        ax.set_xticks(range(1, len(GROUP_SPEC) + 1))
        ax.set_xticklabels([g[1] for g in GROUP_SPEC], fontsize=FONT_TICK)
        ax.tick_params(axis="y", labelsize=FONT_TICK)
        ax.set_title(title, fontsize=FONT_TITLE, pad=12)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    handles = [Patch(fc=COLOR_SBS2, ec="#000000", label="SBS2 arm"),
               Patch(fc=COLOR_CNV, ec="#000000", label="CNV arm")]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=FONT_LEGEND,
               framealpha=0.9, bbox_to_anchor=(0.5, -0.01))
    plt.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(OUT_DIR, f"Metric_distributions.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log("  [SAVE] Metric_distributions")


# =============================================================================
# FIGURE 3: DOMINANCE SPLIT
# =============================================================================

def fig_dominance(m, sets):
    banner("FIGURE 3: A3 dominance split")
    labels, a3a_pct, a3b_pct, n_pos = [], [], [], []
    for key, lab, _ in GROUP_SPEC:
        sub = m.loc[m.index.isin(sets[key])]
        pos = sub[sub["A3_sum"] > 0]
        n = len(pos)
        n_pos.append(n)
        labels.append(lab)
        a = 100.0 * (pos["A3A_fraction"] > 0.5).mean() if n else 0.0
        a3a_pct.append(a)
        a3b_pct.append(100.0 - a if n else 0.0)
        log(f"  {lab.replace(chr(10),' '):16s}: A3+ {n:3d}  "
            f"A3A-dom {a:5.1f}%  A3B-dom {100-a if n else 0:5.1f}%")

    x = np.arange(len(GROUP_SPEC))
    fig, ax = plt.subplots(figsize=(20, 12))
    for xi, a, b in zip(x, a3a_pct, a3b_pct):
        # dominant (larger) fraction on the bottom; color tracks dominance type
        # (A3A coral, A3B mustard), so SBS2 arms read coral-bottom and CNV arms
        # read mustard-bottom.
        if a >= b:
            bot_val, bot_col = a, COLOR_SBS2
            top_val, top_col = b, COLOR_CNV
        else:
            bot_val, bot_col = b, COLOR_CNV
            top_val, top_col = a, COLOR_SBS2
        ax.bar(xi, bot_val, color=bot_col, edgecolor="#000000", lw=1.2)
        ax.bar(xi, top_val, bottom=bot_val, color=top_col, edgecolor="#000000", lw=1.2)
        # percent of the dominant section, centered inside it
        if bot_val > 0:
            ax.text(xi, bot_val / 2.0, f"{bot_val:.1f}%", ha="center", va="center",
                    fontsize=FONT_LEGEND, fontweight="bold", color="#000000")
    for xi, n in zip(x, n_pos):
        ax.text(xi, 101.5, f"n={n}", ha="center", va="bottom", fontsize=FONT_TICK)

    ax.axvline(2.5, color="#999999", ls=":", lw=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=FONT_TICK)
    ax.set_ylabel("% of A3-expressing cells", fontsize=FONT_AXIS)
    ax.set_ylim(0, 108)
    ax.tick_params(axis="y", labelsize=FONT_TICK)
    ax.set_title("A3A / A3B dominance is intrinsic to the populations",
                 fontsize=FONT_TITLE, pad=14)
    legend_handles = [Patch(fc=COLOR_SBS2, ec="#000000", label="A3A-dominant"),
                      Patch(fc=COLOR_CNV, ec="#000000", label="A3B-dominant")]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               fontsize=FONT_LEGEND, framealpha=0.9, bbox_to_anchor=(0.5, -0.02))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    plt.tight_layout(rect=(0, 0.05, 1, 1))
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(OUT_DIR, f"Dominance_split.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log("  [SAVE] Dominance_split")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("Selection Robustness Figures")
    log(f"  Output: {OUT_DIR}")

    sets = load_membership()
    adata, m = load_metrics()

    fig_umap(adata, sets)
    fig_distributions(m, sets)
    fig_dominance(m, sets)

    banner("FIGURES COMPLETE")
    for f in sorted(os.listdir(OUT_DIR)):
        log(f"  {f}")


if __name__ == "__main__":
    main()
