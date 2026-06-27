#!/usr/bin/env python3
"""
Diagnostic_Fig3_Dominance_Boxplots_SBS2_CNV.py
==============================================
2026_NMF_PAPER | Figure 3 diagnostic (Diako review) -- dominance box-and-whisker

Box plots on basal cells expressing EITHER enzyme (A3A + A3B > 0), split by A3
dominance exactly as the slider defines it:
  frac = A3A / (A3A + A3B);  A3A-dominant frac > 0.5,  A3B-dominant frac < 0.5.

  Rows    : top = Tumor, bottom = Normal-adjacent
  Columns : left  = y-axis is SBS2 weight
            right = y-axis is CNV (inferCNV cnv_score)
  Two boxes per panel: A3A-dominant vs A3B-dominant.

TWO FIGURES ARE WRITTEN:
  (1) ..._SBS2_CNV   : all expressing-either basal cells.
  (2) ..._SBS2pos    : same layout, but cells with SBS2 == 0 are dropped, so the
      A3A-dom vs A3B-dom separation on SBS2 is not buried by the zero pile.
      (the whole figure is filtered to SBS2 > 0, so the CNV column there is CNV
      among SBS2-bearing cells.)

PLUMBING matches Generate_Figure3_A3_Dominance_Slider.py exactly:
  basal = final_annotation == "basal cell"
  SBS2  = signature_weights_per_cell.txt row "SBS2" (0.0 where missing)
  CNV   = inferCNV cnv_score
  tissue split on source_name vs NORMAL_SOURCE (GEO typo intentional)

DISPLAY CHOICES:
  - boxes and the Mann-Whitney U test use ALL cells in each group.
  - overlaid jittered points are a random subsample (N_JITTER per group) for
    legibility only; a diamond marks the group MEAN (computed on all cells).
  - y-axis is SHARED DOWN EACH COLUMN (tumor and normal comparable).
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================== CONFIG ==============================
PROJECT_ROOT     = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")
OUTPUT_DIR       = os.path.join(PROJECT_ROOT, "data/FIG_3/figures/TROUBLESHOOTING")

ANNOT_COL  = "final_annotation"
SOURCE_COL = "source_name"
CNV_COL    = "cnv_score"
A3A_GENE   = "APOBEC3A"
A3B_GENE   = "APOBEC3B"

TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"  # GEO typo intentional

DOM_MARGIN = 0.0     # split at frac = 0.5 (matches slider default)

COLOR_A3A  = "#b5341f"   # muted red  (A3A-dominant)
COLOR_A3B  = "#3b6ea5"   # steel blue (A3B-dominant)

N_JITTER   = 1500        # max overlaid points per group (display only)
JIT_W      = 0.16        # horizontal jitter width
PT_ALPHA   = 0.30
Y_PCT      = 99.0        # shared y-limit per column at this percentile (+ headroom)

FS_TICK, FS_LAB, FS_TITLE, FS_ANNOT, FS_LEG = 28, 30, 32, 28, 26
DPI = 300
SEED = 0
# ====================================================================


def log(m=""):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


def load_basal():
    log(f"Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    b = adata[adata.obs[ANNOT_COL].astype(str) == TARGET_CELL_TYPE].copy()
    log(f"  Basal cells: {b.n_obs:,}")

    if "SBS2" not in b.obs.columns:
        w = pd.read_csv(WEIGHTS_PATH, sep="\t", index_col=0)
        sbs2_map = w.loc["SBS2"].to_dict()
        b.obs["SBS2"] = b.obs_names.map(lambda x: sbs2_map.get(x, 0.0)).astype(float)
    b.obs["SBS2"] = b.obs["SBS2"].astype(float)

    genes = list(b.var_names)
    for g in (A3A_GENE, A3B_GENE):
        if g not in genes:
            raise KeyError(f"{g} not found in adata.var_names")
        v = b.X[:, genes.index(g)]
        v = v.toarray().flatten() if hasattr(v, "toarray") else np.asarray(v).flatten()
        b.obs[g] = v
    b.obs[CNV_COL] = b.obs[CNV_COL].astype(float)

    b.obs["tissue_grp"] = np.where(
        b.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, "normal", "tumor")

    df = b.obs[["SBS2", CNV_COL, A3A_GENE, A3B_GENE, "tissue_grp"]].copy()
    df = df[(df[A3A_GENE] + df[A3B_GENE]) > 0].copy()       # expressing either
    df["frac"] = df[A3A_GENE] / (df[A3A_GENE] + df[A3B_GENE])
    df["dom"] = np.where(df["frac"] > 0.5 + DOM_MARGIN, "A3A",
                np.where(df["frac"] < 0.5 - DOM_MARGIN, "A3B", "mid"))
    df = df[df["dom"] != "mid"].copy()
    log(f"  Expressing-either basal: {len(df):,} "
        f"(A3A-dom {int((df['dom']=='A3A').sum()):,}, A3B-dom {int((df['dom']=='A3B').sum()):,})")
    return df


def fmt_p(p):
    if not np.isfinite(p):
        return "p = n/a"
    if p < 1e-300:
        return "p < 1e-300"
    return f"p = {p:.2e}" if p < 0.01 else f"p = {p:.3f}"


def panel(ax, sub, ycol, ylim, rng):
    groups = [("A3A", COLOR_A3A), ("A3B", COLOR_A3B)]
    data = {g: sub.loc[sub["dom"] == g, ycol].values for g, _ in groups}

    bp = ax.boxplot([data["A3A"], data["A3B"]], positions=[0, 1], widths=0.55,
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color="#111111", linewidth=2.6),
                    whiskerprops=dict(color="#111111", linewidth=1.8),
                    capprops=dict(color="#111111", linewidth=1.8),
                    boxprops=dict(linewidth=1.8))
    for patch, (_, color) in zip(bp["boxes"], groups):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
        patch.set_edgecolor("#111111")

    for pos, (g, color) in enumerate(groups):
        vals = data[g]
        if vals.size == 0:
            continue
        show = vals if vals.size <= N_JITTER else rng.choice(vals, N_JITTER, replace=False)
        jx = pos + rng.uniform(-JIT_W, JIT_W, show.size)
        ax.scatter(jx, show, s=10, color=color, alpha=PT_ALPHA,
                   edgecolors="none", rasterized=True, zorder=2)
        ax.plot(pos, float(np.mean(vals)), marker="D", ms=15, mfc="#ffffff",
                mec=color, mew=3.0, zorder=5)

    if data["A3A"].size and data["A3B"].size:
        p = mannwhitneyu(data["A3A"], data["A3B"], alternative="two-sided").pvalue
    else:
        p = np.nan
    ax.text(0.5, 0.97, fmt_p(p), transform=ax.transAxes, ha="center", va="top",
            fontsize=FS_ANNOT, fontweight="bold")

    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["A3A-dom", "A3B-dom"], fontsize=FS_TICK)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return p


def build_figure(df_in, out_stem, suptitle):
    rng = np.random.default_rng(SEED)
    columns = [("SBS2", "SBS2 weight"), (CNV_COL, "CNV (cnv_score)")]
    tissues = [("tumor", "Tumor"), ("normal", "Normal-adjacent")]

    # shared y-limit per column (both tissues), with headroom for the p-value text
    col_ylim = {}
    for col_key, _ in columns:
        lo = float(min(0.0, df_in[col_key].min()))
        hi = float(np.percentile(df_in[col_key].values, Y_PCT)) or 1.0
        span = (hi - lo) or 1.0
        col_ylim[col_key] = (lo - 0.02 * span, hi + 0.18 * span)

    fig, axes = plt.subplots(2, 2, figsize=(18, 18))

    for r, (tkey, tlabel) in enumerate(tissues):
        sub = df_in[df_in["tissue_grp"] == tkey]
        for c, (col_key, col_label) in enumerate(columns):
            ax = axes[r, c]
            p = panel(ax, sub, col_key, col_ylim[col_key], rng)
            ax.set_title(f"{tlabel}  |  {col_label}", fontsize=FS_TITLE, fontweight="bold")
            ax.set_ylabel(col_label, fontsize=FS_LAB)
            n_a = int((sub["dom"] == "A3A").sum())
            n_b = int((sub["dom"] == "A3B").sum())
            log(f"  [{out_stem.split('_')[-1]}] {tlabel:<16} {col_label:<16} | "
                f"A3A-dom n={n_a:,} mean={sub.loc[sub['dom']=='A3A', col_key].mean():.4f} | "
                f"A3B-dom n={n_b:,} mean={sub.loc[sub['dom']=='A3B', col_key].mean():.4f} | "
                f"{fmt_p(p)}")

    leg = [Line2D([0], [0], marker="s", linestyle="none", markerfacecolor=COLOR_A3A,
                  markeredgecolor="#111111", markersize=18, alpha=0.5, label="A3A-dominant"),
           Line2D([0], [0], marker="s", linestyle="none", markerfacecolor=COLOR_A3B,
                  markeredgecolor="#111111", markersize=18, alpha=0.5, label="A3B-dominant"),
           Line2D([0], [0], marker="D", linestyle="none", markerfacecolor="#ffffff",
                  markeredgecolor="#111111", markersize=14, markeredgewidth=2.5, label="group mean")]
    fig.legend(handles=leg, loc="lower center", ncol=3, fontsize=FS_LEG,
               frameon=False, bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(suptitle, fontsize=FS_TITLE + 2, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0.03, 1, 0.98))

    png = os.path.join(OUTPUT_DIR, out_stem + ".png")
    pdf = os.path.join(OUTPUT_DIR, out_stem + ".pdf")
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved:\n  {png}\n  {pdf}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = load_basal()

    # (1) all expressing-either basal cells
    build_figure(
        df,
        "Diagnostic_Fig3_dominance_boxplots_SBS2_CNV",
        "Basal cells (expressing either A3): SBS2 and CNV by A3 dominance",
    )

    # (2) same layout, SBS2 == 0 cells dropped so the SBS2 separation is visible
    df_pos = df[df["SBS2"] > 0].copy()
    log(f"  SBS2 > 0 subset: {len(df_pos):,} cells "
        f"(A3A-dom {int((df_pos['dom']=='A3A').sum()):,}, "
        f"A3B-dom {int((df_pos['dom']=='A3B').sum()):,})")
    build_figure(
        df_pos,
        "Diagnostic_Fig3_dominance_boxplots_SBS2pos",
        "Basal cells (expressing either A3, SBS2 > 0): SBS2 and CNV by A3 dominance",
    )


if __name__ == "__main__":
    main()
