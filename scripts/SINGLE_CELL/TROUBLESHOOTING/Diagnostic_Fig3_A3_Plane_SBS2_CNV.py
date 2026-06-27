#!/usr/bin/env python3
"""
Diagnostic_Fig3_A3_Plane_SBS2_CNV.py
====================================
2026_NMF_PAPER | Figure 3 diagnostic (Diako review) -- A3A-vs-A3B expression plane

Four-panel dot plot on basal cells expressing EITHER enzyme (A3A + A3B > 0).
Every panel: y-axis = APOBEC3A expression, x-axis = APOBEC3B expression.
  Rows    : top = Tumor, bottom = Normal-adjacent
  Columns : left  = dots colored AND sized by SBS2 weight
            right = dots colored AND sized by CNV (inferCNV cnv_score)

Read we are after: in the tumor SBS2 panel, the large bright dots should ride
UP the A3A axis; in the tumor CNV panel they should ride OUT the A3B axis;
normal-adjacent should stay flat in both. This makes the A3A:SBS2 / A3B:CNV
link visible at the single-cell level, which the global Spearman rho (~0.149,
~0.147) understates.

PLUMBING matches Generate_Figure3_A3_Dominance_Slider.py exactly:
  basal = final_annotation == "basal cell"
  SBS2  = signature_weights_per_cell.txt row "SBS2" (0.0 where missing)
  CNV   = inferCNV cnv_score
  tissue split on source_name vs NORMAL_SOURCE (GEO typo intentional)

ENCODING CHOICES:
  - color + size encode the SAME weight value (raw, with a size floor so
    zero-weight cells stay visible as the smallest dots: honest, not hidden).
  - color/size scale is SHARED DOWN EACH COLUMN (tumor and normal comparable).
  - A3A/A3B axis limits are SHARED ACROSS ALL FOUR PANELS.
  - high-weight dots are drawn last (on top of low-weight dots).
  - SBS2 column ramp: light grey -> coral  (#ed6a5a, SBS2 pole)
    CNV  column ramp: light grey -> mustard (#F6D155, CNV  pole)

LAYOUT (v2, requested):
  - colorbars moved off to the right of each column on dedicated axes.
  - wider gap between the two columns.
  - the dot-size key is a single legend in the bottom margin, off the plotting area.
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
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

# encoding
COLOR_SBS2 = "#ed6a5a"   # coral   (SBS2 pole)
COLOR_CNV  = "#F6D155"   # mustard (CNV pole)
GREY_LOW   = "#e8e8e8"   # low-weight end of both ramps (matches slider neutral)
CMAP_SBS2  = LinearSegmentedColormap.from_list("sbs2", [GREY_LOW, COLOR_SBS2])
CMAP_CNV   = LinearSegmentedColormap.from_list("cnv",  [GREY_LOW, COLOR_CNV])

SIZE_MIN, SIZE_MAX = 8, 150     # marker area (pts^2); floor keeps zero-weight cells visible
PT_ALPHA   = 0.80
AXIS_PCT   = 99.5               # expression axis limits clipped at this percentile (display)
SCALE_PCT  = 99.0               # color/size saturate at this percentile of the column

# layout
CBAR_W      = 0.016             # colorbar width (figure fraction)
CBAR_GAP    = 0.016             # gap between a column's right edge and its colorbar
ADJ = dict(left=0.075, right=0.85, top=0.92, bottom=0.15, wspace=0.55, hspace=0.22)

FS_TICK, FS_LAB, FS_TITLE, FS_ANNOT, FS_LEG = 28, 30, 32, 28, 26
DPI = 300
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
    # universe: basal cells expressing EITHER enzyme (drops the double-negative pile)
    df = df[(df[A3A_GENE] + df[A3B_GENE]) > 0].copy()
    log(f"  Expressing-either basal cells: {len(df):,} "
        f"(tumor {int((df['tissue_grp']=='tumor').sum()):,}, "
        f"normal {int((df['tissue_grp']=='normal').sum()):,})")
    return df


def size_for(vals, vmin, vmax):
    frac = np.clip((vals - vmin) / (vmax - vmin) if vmax > vmin else np.zeros_like(vals), 0, 1)
    return SIZE_MIN + (SIZE_MAX - SIZE_MIN) * frac


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = load_basal()

    # shared expression axis limits across all four panels
    xmax = float(np.percentile(df[A3B_GENE].values, AXIS_PCT)) or 1.0
    ymax = float(np.percentile(df[A3A_GENE].values, AXIS_PCT)) or 1.0
    xlim = (-0.03 * xmax, 1.05 * xmax)
    ylim = (-0.03 * ymax, 1.05 * ymax)

    columns = [
        ("SBS2",  "SBS2 weight", CMAP_SBS2),
        (CNV_COL, "CNV (cnv_score)", CMAP_CNV),
    ]
    tissues = [("tumor", "Tumor"), ("normal", "Normal-adjacent")]

    # per-column shared scale (vmin=0, vmax = SCALE_PCT over both tissues)
    col_scale = {}
    for col_key, _, _ in columns:
        vmax = float(np.percentile(df[col_key].values, SCALE_PCT)) or 1.0
        col_scale[col_key] = (0.0, vmax)

    fig, axes = plt.subplots(2, 2, figsize=(22, 18))
    fig.subplots_adjust(**ADJ)

    for r, (tkey, tlabel) in enumerate(tissues):
        sub = df[df["tissue_grp"] == tkey]
        for c, (col_key, col_label, cmap) in enumerate(columns):
            ax = axes[r, c]
            vmin, vmax = col_scale[col_key]
            norm = Normalize(vmin, vmax)

            order = np.argsort(sub[col_key].values)      # low first, high drawn last
            x = sub[A3B_GENE].values[order]
            y = sub[A3A_GENE].values[order]
            w = sub[col_key].values[order]
            sizes = size_for(w, vmin, vmax)

            ax.scatter(x, y, c=w, cmap=cmap, norm=norm, s=sizes, alpha=PT_ALPHA,
                       edgecolors="none", rasterized=True)

            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_title(f"{tlabel}  |  dots = {col_label}", fontsize=FS_TITLE, fontweight="bold")
            ax.tick_params(axis="both", labelsize=FS_TICK)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if c == 0:
                ax.set_ylabel("APOBEC3A expression", fontsize=FS_LAB)
            if r == 1:
                ax.set_xlabel("APOBEC3B expression", fontsize=FS_LAB)

            rho_a3a, _ = spearmanr(sub[A3A_GENE], sub[col_key])
            rho_a3b, _ = spearmanr(sub[A3B_GENE], sub[col_key])
            log(f"  {tlabel:<16} {col_label:<16} n={len(sub):>6,} | "
                f"rho(A3A,{col_key})={rho_a3a:+.3f}  rho(A3B,{col_key})={rho_a3b:+.3f}")

    # one colorbar per column, on a dedicated axis to the right of that column
    for c, (col_key, col_label, cmap) in enumerate(columns):
        pos_top = axes[0, c].get_position()
        pos_bot = axes[1, c].get_position()
        cax = fig.add_axes([pos_top.x1 + CBAR_GAP, pos_bot.y0, CBAR_W, pos_top.y1 - pos_bot.y0])
        sm = ScalarMappable(norm=Normalize(*col_scale[col_key]), cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label(col_label, fontsize=FS_LAB)
        cbar.ax.tick_params(labelsize=FS_TICK)

    # single dot-size key in the bottom margin, off the plotting area.
    # sizing is the same SIZE_MIN..SIZE_MAX mapping for both columns (only the
    # absolute weight differs, which the colorbars show), so the key is relative.
    size_handles = []
    for s, lab in zip([SIZE_MIN, (SIZE_MIN + SIZE_MAX) / 2, SIZE_MAX], ["low", "mid", "high"]):
        size_handles.append(Line2D([0], [0], marker="o", linestyle="none",
                                   markerfacecolor="#9aa0a6", markeredgecolor="none",
                                   markersize=np.sqrt(s), label=lab))
    fig.legend(handles=size_handles, title="dot size and color = weight (per-column scale)",
               loc="lower center", ncol=3, fontsize=FS_LEG, title_fontsize=FS_LEG,
               frameon=False, bbox_to_anchor=(0.5, 0.015), handletextpad=0.6, columnspacing=2.4)

    fig.suptitle("Basal cells (expressing either A3): A3A-vs-A3B plane by mutational weight",
                 fontsize=FS_TITLE + 2, fontweight="bold", y=0.975)

    png = os.path.join(OUTPUT_DIR, "Diagnostic_Fig3_A3_plane_SBS2_CNV.png")
    pdf = os.path.join(OUTPUT_DIR, "Diagnostic_Fig3_A3_plane_SBS2_CNV.pdf")
    fig.savefig(png, dpi=DPI)
    fig.savefig(pdf)
    plt.close(fig)
    log(f"Saved:\n  {png}\n  {pdf}")


if __name__ == "__main__":
    main()
