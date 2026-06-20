#!/usr/bin/env python3
"""
Generate_Figure3_A3_Dominance_Slider.py
=======================================
2026_NMF_PAPER | NEW Figure 3 (v3 restructure) -- DOMINANCE version

Four-row slider showing where A3A-dominant vs A3B-dominant basal cells sit on a
within-tissue SBS2-vs-CNV axis, in tumor and in normal-adjacent tissue.
  right pole = SBS2-high,  left pole = CNV-high (productive).

WHY THIS DESIGN (breadcrumbs; all numbers from Diagnostic_Fig3_Lean_Decompose.py)
--------------------------------------------------------------------------------
DROPPED, and why:
  (1) STEMNESS (CytoTRACE2) is no longer blended into the productive pole.
      The A3A fraction anti-correlates with stemness in BOTH tumor (rho=-0.50)
      and normal (rho=-0.67), so stemness is NOT tumor-specific and blending it
      in masks the contrast we want. Evidence: frac-vs-CNV = -0.437 (tumor) vs
      +0.070 (normal) is tumor-specific, but frac-vs-(CNV+stem) = -0.538 (tumor)
      vs -0.590 (normal) is not. So the productive pole is CNV ALONE.
  (2) The A3>0 per-enzyme conditioning is gone. Conditioning on expressing cells
      deleted the co-occurrence signal (all-basal A3A-vs-SBS2 = +0.140 collapsed
      to -0.035 among A3A>0 cells). Dominance is instead defined among cells
      expressing EITHER enzyme.
  (3) Absolute A3 expression is no longer the axis driver. SBS2 does not separate
      the two enzymes (A3A-vs-SBS2 +0.140, A3B-vs-SBS2 +0.108 on all basal);
      absolute levels dilute the signal.

ADDED, and why:
  (1) A3 DOMINANCE fraction  frac = A3A / (A3A + A3B), among basal cells with
      A3A+A3B > 0. Cells split into A3A-dominant (frac > 0.5 + MARGIN) and
      A3B-dominant (frac < 0.5 - MARGIN). The fraction isolates which enzyme is
      in play and differentiates strongly on CNV (frac-vs-CNV = -0.437 in tumor)
      where absolute levels did not.
  (2) x-axis = within-tissue  z(SBS2) - z(CNV)  (ranks standardized WITHIN each
      tissue's expressing-either universe, so "center" is that tissue's own
      median and the normal/tumor marginal differences don't push normal off
      center artefactually). CNV is the inferCNV cnv_score (24 discrete levels);
      a small display-only jitter (X_JITTER) breaks the banding. Medians and
      means are computed on the UNJITTERED values.

  Honest caveat for the text: the separation is carried mainly by CNV (low CNV
  on the A3A-dominant side), not by high SBS2, since SBS2 is shared between the
  enzymes. This fits the maintenance(low-CNV)/productive(high-CNV) framing.

CHANGELOG
  v1  scatter SBS2 vs cnv_score                                   [shelved]
  v2  scatter SBS2 vs CNV+Stemness blend                          [shelved]
  v3  lean slider, z_A3*(z_SBS2 - z_cnvstem), expressing cells    [shelved]
  v4  THIS: dominance slider; drop stemness + A3>0 cond. + abs-level driver;
      group by A3A/(A3A+A3B); x = within-tissue z(SBS2)-z(CNV), CNV alone jittered
  v4.1 marker nodes (median=circle, mean=diamond), bold zero line, point alpha 0.80

CONFIRMED CONSTANTS (match Step00B / GEO):
  TARGET_CELL_TYPE = "basal cell"
  NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata, spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

# ============================== CONFIG ==============================
PROJECT_ROOT     = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")
OUTPUT_DIR       = os.path.join(PROJECT_ROOT, "data/FIG_3/figures")

ANNOT_COL  = "final_annotation"
SOURCE_COL = "source_name"
CNV_COL    = "cnv_score"          # inferCNV; 24 discrete levels -> jitter for display
A3A_GENE   = "APOBEC3A"
A3B_GENE   = "APOBEC3B"
# NOTE: CytoTRACE2_Score is intentionally NOT loaded -- stemness dropped (see header).

TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"  # GEO typo intentional

DOM_MARGIN = 0.0   # 0.0 -> split at frac=0.5; set 0.1 for a cleaner 0.6/0.4 split (drops ambiguous middle)
X_JITTER   = 0.07  # display-only horizontal jitter to break the 24-level CNV banding

COLOR_SBS2 = "#ed6a5a"     # coral   -> right pole (SBS2-high)
COLOR_CNV  = "#F6D155"     # mustard -> left pole (CNV-high / productive)
CMAP_X     = LinearSegmentedColormap.from_list("x", [COLOR_CNV, "#e8e8e8", COLOR_SBS2])
POINT_ALPHA = 0.95         # high alpha so sparse normal-adjacent points stay visible
MARK_MEAN   = "#2b50aa"    # mean-diamond accent color (pops against the mustard/coral cloud)
FS_TICK, FS_LAB, FS_POLE, FS_ANNOT = 28, 30, 30, 26
DPI = 300
# ====================================================================


def log(m=""):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


def zrank(x):
    r = rankdata(np.asarray(x, dtype=float))     # average ranks for ties
    s = r.std()
    return (r - r.mean()) / s if s > 0 else np.zeros_like(r)


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
        v = b.X[:, genes.index(g)]
        v = v.toarray().flatten() if hasattr(v, "toarray") else np.asarray(v).flatten()
        b.obs[g] = v
    b.obs[CNV_COL] = b.obs[CNV_COL].astype(float)

    b.obs["tissue_grp"] = np.where(
        b.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, "normal", "tumor")
    log(f"  Tumor: {(b.obs['tissue_grp']=='tumor').sum():,} | "
        f"Normal-adjacent: {(b.obs['tissue_grp']=='normal').sum():,}")
    return b.obs[["SBS2", CNV_COL, A3A_GENE, A3B_GENE, "tissue_grp"]].copy()


def process_tissue(obs, grp):
    """Among basal cells expressing EITHER A3 in this tissue: compute the A3A
    dominance fraction, the within-tissue SBS2-vs-CNV position, and split by dominance."""
    g = obs[obs["tissue_grp"] == grp].copy()
    g = g[(g[A3A_GENE] + g[A3B_GENE]) > 0].copy()
    g["frac"] = g[A3A_GENE] / (g[A3A_GENE] + g[A3B_GENE])
    g["zs"] = zrank(g["SBS2"].values)     # standardized WITHIN this tissue's universe
    g["zc"] = zrank(g[CNV_COL].values)
    g["x"]  = g["zs"] - g["zc"]
    a3a = g[g["frac"] > 0.5 + DOM_MARGIN]
    a3b = g[g["frac"] < 0.5 - DOM_MARGIN]
    rs, _ = spearmanr(g["frac"], g["SBS2"])
    rc, _ = spearmanr(g["frac"], g[CNV_COL])
    log(f"  {grp}: either-expressing n={len(g):,} | frac-vs-SBS2={rs:+.3f} frac-vs-CNV={rc:+.3f}")
    return a3a, a3b


def build_slider(rows, out_png, out_pdf=None, seed=0):
    """rows: list (top->bottom) of dicts with keys 'label' and 'vals' (np.array of x)."""
    rng = np.random.default_rng(seed)
    ycs = list(range(len(rows) - 1, -1, -1))
    allv = np.concatenate([r["vals"] for r in rows if len(r["vals"])])
    lim = float(np.percentile(np.abs(allv), 99)) or 1.0
    norm = Normalize(-lim, lim)

    fig, ax = plt.subplots(figsize=(17, 12))
    BAND = 0.30
    for r, yc in zip(rows, ycs):
        v = r["vals"]
        if v.size == 0:
            ax.text(0, yc, "no cells", fontsize=FS_ANNOT, ha="center", va="center", color="#999")
            continue
        xd = np.clip(v + rng.uniform(-X_JITTER, X_JITTER, v.size), -lim, lim)  # jitter = display only
        yd = yc + rng.uniform(-BAND, BAND, v.size)
        ax.scatter(xd, yd, c=v, cmap=CMAP_X, norm=norm, s=50, alpha=POINT_ALPHA,
                   edgecolors="none", rasterized=True, zorder=2)
        med, mean = float(np.median(v)), float(np.mean(v))   # computed on UNJITTERED values
        # median = solid line + white circle node; mean = dashed line + colored diamond node
        ax.plot([med, med], [yc - 0.30, yc + 0.30], color="#111111", lw=2.4, zorder=5)
        ax.plot([mean, mean], [yc - 0.30, yc + 0.30], color="#111111", lw=2.4, ls="--", zorder=5)
        ax.plot(med, yc, marker="o", ms=14, mfc="white", mec="#111111", mew=2.4, zorder=6)
        ax.plot(mean, yc, marker="D", ms=13, mfc=MARK_MEAN, mec="#111111", mew=1.8, zorder=6)
        ax.text(-lim * 0.99, yc + 0.46, f"n={v.size:,}   median={med:+.2f}   mean={mean:+.2f}",
                fontsize=FS_ANNOT, ha="left", va="center", color="#333333")

    ax.axvline(0, color="#5a5a5a", ls="-", lw=2.2, alpha=0.65, zorder=3)  # bold zero line, above cloud
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.7, len(rows) - 1 + 0.95)
    ax.set_yticks(ycs)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=FS_TICK)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    ax.set_xlabel("per-cell position:  z(SBS2) \u2212 z(CNV)   [within tissue]", fontsize=FS_LAB)

    top = len(rows) - 1 + 0.72
    ax.text(-lim * 0.98, top, "\u2190 CNV-high (productive)", fontsize=FS_POLE,
            fontweight="bold", color="#a87f12", ha="left", va="center")
    ax.text(lim * 0.98, top, "SBS2-high \u2192", fontsize=FS_POLE,
            fontweight="bold", color="#b5341f", ha="right", va="center")
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.4)

    leg = [Line2D([0], [0], color="#111111", lw=2.4, marker="o", mfc="white",
                  mec="#111111", mew=2, ms=12, label="median"),
           Line2D([0], [0], color="#111111", lw=2.4, ls="--", marker="D", mfc=MARK_MEAN,
                  mec="#111111", mew=1.5, ms=11, label="mean")]
    ax.legend(handles=leg, loc="lower right", fontsize=FS_ANNOT, frameon=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    if out_pdf:
        fig.savefig(out_pdf)
    plt.close(fig)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    obs = load_basal()
    a3a_t, a3b_t = process_tissue(obs, "tumor")
    a3a_n, a3b_n = process_tissue(obs, "normal")
    rows = [
        {"label": "A3A-dom \u00b7 Tumor",           "vals": a3a_t["x"].values},
        {"label": "A3B-dom \u00b7 Tumor",           "vals": a3b_t["x"].values},
        {"label": "A3A-dom \u00b7 Normal-adjacent", "vals": a3a_n["x"].values},
        {"label": "A3B-dom \u00b7 Normal-adjacent", "vals": a3b_n["x"].values},
    ]
    for r, d in zip(rows, [a3a_t, a3b_t, a3a_n, a3b_n]):
        if len(d):
            log(f"  {r['label']:<26} n={len(d):>6,}  median(x)={np.median(d['x']):+.3f}  "
                f"mean(x)={np.mean(d['x']):+.3f}  | med z(SBS2)={np.median(d['zs']):+.2f} "
                f"med z(CNV)={np.median(d['zc']):+.2f}")
    png = os.path.join(OUTPUT_DIR, "Figure3_A3_dominance_slider.png")
    pdf = os.path.join(OUTPUT_DIR, "Figure3_A3_dominance_slider.pdf")
    build_slider(rows, png, pdf)
    log(f"Saved:\n  {png}\n  {pdf}")


if __name__ == "__main__":
    main()
