#!/usr/bin/env python3
"""
Generate_Figure3_A3_SBS2_CNV_Scatter.py
=======================================
2026_NMF_PAPER | NEW Figure 3 (v3 restructure)

Basal-cell SBS2 (x) vs a combined CNV + Stemness phenotype (y), points colored
by A3A or A3B expression, split tumor (left column) vs normal-adjacent (right).

Purpose: introduce the A3A<->SBS2-HIGH / A3B<->CNV-HIGH color key and show that
A3A tracks elevated SBS2 while A3B tracks the chromosomally-unstable / stem-like
state, motivating the group selection used in Figure 4.

WHY A COMBINED Y-AXIS
  The raw inferCNV cnv_score is binned to only 24 discrete levels, which renders
  as horizontal tracks. CytoTRACE2_Score (stemness) is continuous (~155k unique
  values) and positively correlated with cnv_score, and the CNV-HIGH group in the
  methods is already defined on a composite of both. We therefore plot a magnitude
  blend of the two min-max-normalized scores, which is smooth and stays anchored
  to the same A3B-associated cancer pole. This is NOT the old "agreement_score"
  (a QC concordance metric that peaks at both extremes); it is a weighted average
  so that high-CNV/high-stemness cells sit high and low/low cells sit low.

PROTECTED CONSTRAINTS
  - A3A is colored with the SBS2-HIGH coral (#ed6a5a); A3B with the CNV-HIGH
    mustard (#F6D155). Colormap base is gray (#9e9e9e) so zero/low cells stay
    visible. This figure INTRODUCES that key; keep it consistent downstream.
  - Load + field logic mirrors Step00B_Apply_A3CNV_Groups_and_Export.py:
      basal      = obs['final_annotation'] == TARGET_CELL_TYPE
      SBS2       = mapped per-barcode from the weights file -> obs['SBS2']
      A3A/A3B    = pulled from adata.X (var_names 'APOBEC3A' / 'APOBEC3B')
      CNV        = obs['cnv_score'] ; stemness = obs['CytoTRACE2_Score']
      tumor/norm = obs['source_name'] == NORMAL_SOURCE -> normal, else tumor
  - Fonts 28-34, hex colors only, PDF + PNG at 300 DPI, no descriptive titles
    over panels (column/row headers only).

CHANGELOG
  v1  initial build (SBS2 vs raw cnv_score)
  v2  y-axis = magnitude blend of normalized cnv_score + CytoTRACE2_Score
      ("CNV + Stemness phenotype"), to remove 24-level cnv_score banding;
      SBS2 x-axis clipped at SBS2_CLIP_PCTL; bottom-left rho now vs the blend;
      prints Spearman(cnv_score, CytoTRACE2_Score) as a direction check.

CONFIRMED CONSTANTS (match Step00B / GEO):
  TARGET_CELL_TYPE = "basal cell"
  NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ============================== CONFIG ==============================
PROJECT_ROOT     = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")
OUTPUT_DIR       = os.path.join(PROJECT_ROOT, "data/FIG_3/figures")

ANNOT_COL  = "final_annotation"
SOURCE_COL = "source_name"
CNV_COL    = "cnv_score"
CYTO_COL   = "CytoTRACE2_Score"
A3A_GENE   = "APOBEC3A"
A3B_GENE   = "APOBEC3B"

TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"  # GEO typo intentional

# combined y-axis = W_CNV * norm(cnv_score) + W_CYTO * norm(CytoTRACE2_Score)
W_CNV, W_CYTO = 0.5, 0.5      # tunable; raise W_CNV to keep the axis more CNV-anchored
Y_LABEL = "CNV + Stemness phenotype"

# SBS2 x-axis: clip the long upper tail (max ~11.5, p99 ~2) so the bulk isn't crushed
SBS2_CLIP_PCTL = 99.5

# True  = plot every basal cell (no detectable SBS2 sits at x=0)
# False = plot only basal cells with a detectable SBS2 weight (> 0)
PLOT_ALL_BASAL = True

# standard NMF palette (this figure introduces the A3 color key); gray base = visible zeros
COLOR_SBS2_HIGH = "#ed6a5a"   # coral   -> A3A
COLOR_CNV_HIGH  = "#F6D155"   # mustard -> A3B
CMAP_A3A = LinearSegmentedColormap.from_list("a3a", ["#9e9e9e", "#f4a294", COLOR_SBS2_HIGH, "#b5341f"])
CMAP_A3B = LinearSegmentedColormap.from_list("a3b", ["#9e9e9e", "#f6e08f", COLOR_CNV_HIGH, "#a87f12"])

DPI = 300
FS_TICK, FS_LAB, FS_HEAD, FS_CBAR = 28, 32, 34, 28
S_TUMOR, S_NORMAL = 10, 24          # normal-adjacent is sparse, so larger points
# ====================================================================


def log(msg=""):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def minmax(x):
    x = np.asarray(x, dtype=float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    return (x - lo) / (hi - lo) if hi > lo else np.zeros_like(x)


def load_basal():
    """Mirror Step00B: subset basal, map SBS2 from weights, pull A3 from X."""
    if not os.path.exists(ADATA_FINAL_PATH):
        log(f"ERROR: not found: {ADATA_FINAL_PATH}")
        sys.exit(1)
    log(f"Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    if TARGET_CELL_TYPE not in set(adata.obs[ANNOT_COL].astype(str)):
        log(f"ERROR: TARGET_CELL_TYPE '{TARGET_CELL_TYPE}' not in {ANNOT_COL}.")
        log(f"  Available: {sorted(adata.obs[ANNOT_COL].astype(str).unique())}")
        sys.exit(1)
    b = adata[adata.obs[ANNOT_COL].astype(str) == TARGET_CELL_TYPE].copy()
    log(f"  Basal cells: {b.n_obs:,}")

    # SBS2 weights (per-barcode map, exactly as Step00B)
    if "SBS2" not in b.obs.columns:
        w = pd.read_csv(WEIGHTS_PATH, sep="\t", index_col=0)
        sbs2_map = w.loc["SBS2"].to_dict()
        b.obs["SBS2"] = b.obs_names.map(lambda x: sbs2_map.get(x, 0.0)).astype(float)
    b.obs["SBS2"] = b.obs["SBS2"].astype(float)

    # A3A / A3B from expression matrix
    genes = list(b.var_names)
    for g in (A3A_GENE, A3B_GENE):
        if g in genes:
            v = b.X[:, genes.index(g)]
            v = v.toarray().flatten() if hasattr(v, "toarray") else np.asarray(v).flatten()
            b.obs[g] = v
        else:
            log(f"  WARNING: {g} not in var_names; setting 0.0")
            b.obs[g] = 0.0

    # CNV + stemness
    for col in (CNV_COL, CYTO_COL):
        if col not in b.obs.columns:
            numeric = [c for c in b.obs.columns if b.obs[c].dtype.kind in "fi"]
            log(f"ERROR: '{col}' not in obs. Numeric columns available: {numeric}")
            sys.exit(1)
        b.obs[col] = b.obs[col].astype(float)

    # tumor vs normal-adjacent
    if NORMAL_SOURCE not in set(b.obs[SOURCE_COL].astype(str)):
        log(f"ERROR: NORMAL_SOURCE '{NORMAL_SOURCE}' not in {SOURCE_COL}.")
        log(f"  Available: {sorted(b.obs[SOURCE_COL].astype(str).unique())}")
        sys.exit(1)
    b.obs["tissue_grp"] = np.where(
        b.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, "normal", "tumor")
    log(f"  Tumor: {(b.obs['tissue_grp'] == 'tumor').sum():,} | "
        f"Normal-adjacent: {(b.obs['tissue_grp'] == 'normal').sum():,}")

    return b.obs[["SBS2", CNV_COL, CYTO_COL, A3A_GENE, A3B_GENE, "tissue_grp"]].copy()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    obs = load_basal()

    # direction check (should be positive: both up in cancer cells)
    r_all, p_all = spearmanr(obs[CNV_COL], obs[CYTO_COL])
    tmask = obs["tissue_grp"] == "tumor"
    r_tum, p_tum = spearmanr(obs.loc[tmask, CNV_COL], obs.loc[tmask, CYTO_COL])
    log(f"  Spearman cnv_score vs CytoTRACE2: all basal rho={r_all:.3f} (p={p_all:.1e}); "
        f"tumor rho={r_tum:.3f} (p={p_tum:.1e})")
    if r_all < 0:
        log("  WARNING: negative correlation -- blend would cancel; reconsider before trusting.")

    # combined CNV + Stemness phenotype (min-max over ALL basal, then weighted avg)
    obs["cnv_stem"] = W_CNV * minmax(obs[CNV_COL].values) + W_CYTO * minmax(obs[CYTO_COL].values)

    # clip SBS2 upper tail for display (Spearman below uses the raw values)
    cap = float(np.nanpercentile(obs["SBS2"].values, SBS2_CLIP_PCTL))
    cap = max(cap, 1e-6)
    n_pinned = int((obs["SBS2"] > cap).sum())
    obs["SBS2_plot"] = np.minimum(obs["SBS2"].values, cap)
    log(f"  SBS2 x clipped at p{SBS2_CLIP_PCTL} = {cap:.3f} ({n_pinned} cells pinned to edge)")

    def split(grp):
        d = obs[obs["tissue_grp"] == grp]
        return d if PLOT_ALL_BASAL else d[d["SBS2"] > 0]

    tum, nor = split("tumor"), split("normal")
    log(f"  Plotting tumor n={len(tum):,} | normal n={len(nor):,} (PLOT_ALL_BASAL={PLOT_ALL_BASAL})")

    a3a_vmax = float(np.percentile(obs[A3A_GENE].values, 99)) or 1.0
    a3b_vmax = float(np.percentile(obs[A3B_GENE].values, 99)) or 1.0

    fig = plt.figure(figsize=(18, 15))
    L, R, B, T = 0.11, 0.85, 0.09, 0.90
    gs = fig.add_gridspec(2, 2, left=L, right=R, bottom=B, top=T, wspace=0.14, hspace=0.16)
    ax = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(2)]

    def panel(a, d, gene, cmap, vmax, s):
        return a.scatter(d["SBS2_plot"], d["cnv_stem"], c=d[gene], cmap=cmap,
                         vmin=0, vmax=vmax, s=s, alpha=0.80,
                         edgecolors="none", rasterized=True)

    sc_a = panel(ax[0][0], tum, A3A_GENE, CMAP_A3A, a3a_vmax, S_TUMOR)
    panel(ax[0][1], nor, A3A_GENE, CMAP_A3A, a3a_vmax, S_NORMAL)
    sc_b = panel(ax[1][0], tum, A3B_GENE, CMAP_A3B, a3b_vmax, S_TUMOR)
    panel(ax[1][1], nor, A3B_GENE, CMAP_A3B, a3b_vmax, S_NORMAL)

    # rho where each enzyme should track: A3A vs SBS2 (raw), A3B vs combined phenotype
    def rho(a, d, gene, axis_col):
        if len(d) > 10:
            r, _ = spearmanr(d[gene], d[axis_col])
            a.text(0.96, 0.95, f"rho = {r:.2f}", transform=a.transAxes,
                   ha="right", va="top", fontsize=FS_TICK - 4, color="#333333")
    rho(ax[0][0], tum, A3A_GENE, "SBS2")
    rho(ax[1][0], tum, A3B_GENE, "cnv_stem")

    for row in ax:
        for a in row:
            a.set_xlim(-0.02 * cap, cap * 1.02)
            a.set_ylim(-0.02, 1.02)
            a.tick_params(labelsize=FS_TICK, length=6, width=1.3)
            for sp in a.spines.values():
                sp.set_linewidth(1.4)

    ax[1][0].set_xlabel("SBS2 weight", fontsize=FS_LAB)
    ax[1][1].set_xlabel("SBS2 weight", fontsize=FS_LAB)
    ax[0][0].set_ylabel(Y_LABEL, fontsize=FS_LAB)
    ax[1][0].set_ylabel(Y_LABEL, fontsize=FS_LAB)
    ax[0][0].set_title("Tumor", fontsize=FS_HEAD, pad=16)
    ax[0][1].set_title("Normal-adjacent", fontsize=FS_HEAD, pad=16)

    fig.text(0.028, 0.705, "A3A", fontsize=FS_HEAD, rotation=90,
             va="center", ha="center", color="#b5341f")
    fig.text(0.028, 0.285, "A3B", fontsize=FS_HEAD, rotation=90,
             va="center", ha="center", color="#a87f12")

    cax_a = fig.add_axes([0.875, 0.53, 0.020, 0.36])
    cb_a = fig.colorbar(sc_a, cax=cax_a)
    cb_a.set_label("A3A expression", fontsize=FS_CBAR)
    cb_a.ax.tick_params(labelsize=FS_TICK - 4)

    cax_b = fig.add_axes([0.875, 0.10, 0.020, 0.36])
    cb_b = fig.colorbar(sc_b, cax=cax_b)
    cb_b.set_label("A3B expression", fontsize=FS_CBAR)
    cb_b.ax.tick_params(labelsize=FS_TICK - 4)

    png = os.path.join(OUTPUT_DIR, "Figure3_A3_SBS2_CNV_scatter.png")
    pdf = os.path.join(OUTPUT_DIR, "Figure3_A3_SBS2_CNV_scatter.pdf")
    fig.savefig(png, dpi=DPI)
    fig.savefig(pdf)
    plt.close(fig)
    log(f"Saved:\n  {png}\n  {pdf}")


if __name__ == "__main__":
    main()
