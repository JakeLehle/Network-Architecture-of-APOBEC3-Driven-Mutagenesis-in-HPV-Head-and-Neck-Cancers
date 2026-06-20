#!/usr/bin/env python3
"""
Generate_Figure3_A3_Lean_Slider.py
===================================
2026_NMF_PAPER | NEW Figure 3 (v3 restructure)

A per-cell "lean" slider in four rows (A3A tumor/normal, A3B tumor/normal),
showing which mutational fate each enzyme's expression co-ranks with:
  right pole = SBS2,  left pole = CNV + Stemness phenotype.

PER-CELL LEAN (confirmed metric)
  Within each row's enzyme-EXPRESSING subset (A3X > 0), standardize the ranks
  (z of average ranks) of A3X, SBS2, and the CNV+Stemness composite, then:
      lean(i) = z_A3X(i) * ( z_SBS2(i) - z_CNVstem(i) )
  Sign sets the side (positive -> SBS2/right), magnitude sets distance from center.
  A negative association with one fate amplifies the lean toward the other
  (subtracting a negative). The MEAN of lean over a row equals exactly
  rho(A3X, SBS2) - rho(A3X, CNV+Stem); the MEDIAN is the robust center.
  Both are drawn (mean dashed, median solid) and printed.

WHY EXPRESSING-CELLS-ONLY
  Restricting to A3X > 0 removes the shared-zero mass that would otherwise
  inflate the rank agreement and make tumor vs normal non-comparable.

PROTECTED CONSTRAINTS
  - Right pole / SBS2 = coral (#ed6a5a); left pole / CNV+Stem = mustard (#F6D155);
    points colored by lean on a mustard->gray->coral diverging map.
  - CNV+Stemness composite = 0.5*minmax(cnv_score) + 0.5*minmax(CytoTRACE2_Score)
    over all basal cells (same blend as the scatter version).
  - Load/field logic mirrors Step00B (final_annotation, SBS2 from weights,
    APOBEC3A/B from X, cnv_score, CytoTRACE2_Score, source_name split).
  - Fonts 28-34, hex colors, PDF + PNG at 300 DPI, no titles over the panel.

CHANGELOG
  v1  scatter (SBS2 vs cnv_score)            [shelved]
  v2  scatter (SBS2 vs CNV+Stemness blend)   [shelved]
  v3  per-cell lean slider, 4 rows, expressing cells only, mean + median

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
CNV_COL    = "cnv_score"
CYTO_COL   = "CytoTRACE2_Score"
A3A_GENE   = "APOBEC3A"
A3B_GENE   = "APOBEC3B"

TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"  # GEO typo intentional

W_CNV, W_CYTO = 0.5, 0.5   # CNV+Stemness composite weights

COLOR_SBS2 = "#ed6a5a"     # coral   -> right pole (SBS2)
COLOR_CNV  = "#F6D155"     # mustard -> left pole (CNV + Stemness)
CMAP_LEAN  = LinearSegmentedColormap.from_list("lean", [COLOR_CNV, "#e8e8e8", COLOR_SBS2])
FS_TICK, FS_LAB, FS_POLE, FS_ANNOT = 28, 30, 30, 26
DPI = 300
# ====================================================================


def log(msg=""):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def minmax(x):
    x = np.asarray(x, dtype=float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    return (x - lo) / (hi - lo) if hi > lo else np.zeros_like(x)


def zrank(x):
    r = rankdata(np.asarray(x, dtype=float))     # average ranks for ties
    s = r.std()                                  # ddof=0 -> mean(z_a*z_b) == Spearman rho
    return (r - r.mean()) / s if s > 0 else np.zeros_like(r)


def compute_lean(a3x, sbs2, cnv):
    """lean = z_A3 * (z_SBS2 - z_CNVstem). mean(lean) == rho(A3,SBS2) - rho(A3,CNVstem)."""
    za, zs, zc = zrank(a3x), zrank(sbs2), zrank(cnv)
    return za * (zs - zc)


def build_slider(rows, out_png, out_pdf=None, seed=0):
    """rows: list (top->bottom) of dicts with keys 'label' and 'lean' (np.array)."""
    rng = np.random.default_rng(seed)
    n_rows = len(rows)
    ycs = list(range(n_rows - 1, -1, -1))

    all_lean = np.concatenate([r["lean"] for r in rows if r["lean"].size])
    lim = float(np.percentile(np.abs(all_lean), 99)) or 1.0
    norm = Normalize(-lim, lim)

    fig, ax = plt.subplots(figsize=(17, 12))
    BAND = 0.30
    for r, yc in zip(rows, ycs):
        if r["lean"].size == 0:
            ax.text(0, yc, "no expressing cells", fontsize=FS_ANNOT, ha="center", va="center", color="#999")
            continue
        x = np.clip(r["lean"], -lim, lim)
        y = yc + rng.uniform(-BAND, BAND, size=x.size)
        ax.scatter(x, y, c=r["lean"], cmap=CMAP_LEAN, norm=norm,
                   s=9, alpha=0.45, edgecolors="none", rasterized=True, zorder=2)
        med, mean = float(np.median(r["lean"])), float(np.mean(r["lean"]))
        ax.plot([med, med], [yc - 0.34, yc + 0.34], color="black", lw=3.0, zorder=4)
        ax.plot([mean, mean], [yc - 0.34, yc + 0.34], color="black", lw=3.0, ls="--", zorder=4)
        ax.text(-lim * 0.99, yc + 0.46,
                f"n={r['lean'].size:,}   median={med:+.2f}   mean (\u0394\u03c1)={mean:+.2f}",
                fontsize=FS_ANNOT, ha="left", va="center", color="#333333")

    ax.axvline(0, color="#444444", ls=":", lw=1.6, zorder=1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.7, n_rows - 1 + 0.95)
    ax.set_yticks(ycs)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=FS_TICK)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    ax.set_xlabel("per-cell lean   =   z(A3) \u00d7 ( z(SBS2) \u2212 z(CNV+Stemness) )", fontsize=FS_LAB)

    top = n_rows - 1 + 0.72
    ax.text(-lim * 0.98, top, "\u2190 CNV + Stemness", fontsize=FS_POLE, fontweight="bold",
            color="#a87f12", ha="left", va="center")
    ax.text(lim * 0.98, top, "SBS2 \u2192", fontsize=FS_POLE, fontweight="bold",
            color="#b5341f", ha="right", va="center")
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.4)

    leg = [Line2D([0], [0], color="black", lw=3, label="median"),
           Line2D([0], [0], color="black", lw=3, ls="--", label="mean (\u0394\u03c1)")]
    ax.legend(handles=leg, loc="lower right", fontsize=FS_ANNOT, frameon=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    if out_pdf:
        fig.savefig(out_pdf)
    plt.close(fig)


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
        if g in genes:
            v = b.X[:, genes.index(g)]
            v = v.toarray().flatten() if hasattr(v, "toarray") else np.asarray(v).flatten()
            b.obs[g] = v
        else:
            log(f"  WARNING: {g} not in var_names; setting 0.0")
            b.obs[g] = 0.0

    for col in (CNV_COL, CYTO_COL):
        b.obs[col] = b.obs[col].astype(float)

    # CNV + Stemness composite (over all basal)
    b.obs["cnv_stem"] = W_CNV * minmax(b.obs[CNV_COL].values) + W_CYTO * minmax(b.obs[CYTO_COL].values)

    b.obs["tissue_grp"] = np.where(
        b.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, "normal", "tumor")
    log(f"  Tumor: {(b.obs['tissue_grp']=='tumor').sum():,} | "
        f"Normal-adjacent: {(b.obs['tissue_grp']=='normal').sum():,}")
    return b.obs[["SBS2", "cnv_stem", A3A_GENE, A3B_GENE, "tissue_grp"]].copy()


def row_lean(obs, gene, grp, label):
    sub = obs[(obs["tissue_grp"] == grp) & (obs[gene] > 0)]
    if len(sub) < 3:
        log(f"  {label}: {len(sub)} expressing cells (too few)")
        return {"label": label, "lean": np.array([])}
    lean = compute_lean(sub[gene].values, sub["SBS2"].values, sub["cnv_stem"].values)
    rs, _ = spearmanr(sub[gene], sub["SBS2"])
    rc, _ = spearmanr(sub[gene], sub["cnv_stem"])
    log(f"  {label}: n={len(sub):,}  median={np.median(lean):+.3f}  mean={np.mean(lean):+.3f}  "
        f"| rho(SBS2)={rs:+.3f} rho(CNV+Stem)={rc:+.3f}  check[rho_s-rho_c]={rs-rc:+.3f}")
    return {"label": label, "lean": lean}


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    obs = load_basal()
    rows = [
        row_lean(obs, A3A_GENE, "tumor",  "A3A \u00b7 Tumor"),
        row_lean(obs, A3A_GENE, "normal", "A3A \u00b7 Normal-adjacent"),
        row_lean(obs, A3B_GENE, "tumor",  "A3B \u00b7 Tumor"),
        row_lean(obs, A3B_GENE, "normal", "A3B \u00b7 Normal-adjacent"),
    ]
    png = os.path.join(OUTPUT_DIR, "Figure3_A3_lean_slider.png")
    pdf = os.path.join(OUTPUT_DIR, "Figure3_A3_lean_slider.pdf")
    build_slider(rows, png, pdf)
    log(f"Saved:\n  {png}\n  {pdf}")


if __name__ == "__main__":
    main()
