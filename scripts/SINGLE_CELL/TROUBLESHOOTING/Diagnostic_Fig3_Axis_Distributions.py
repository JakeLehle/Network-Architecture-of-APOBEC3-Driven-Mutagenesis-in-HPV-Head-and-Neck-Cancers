#!/usr/bin/env python3
"""
Diagnostic_Fig3_Axis_Distributions.py
=====================================
2026_NMF_PAPER | Fig 3 axis diagnostic

Profiles the four Fig 3 quantities on basal cells (SBS2, cnv_score, APOBEC3A,
APOBEC3B): range, spread, number of unique values, smallest gap between unique
values (quantization / "track" detection), zero fraction, and the most common
values. Saves histograms. Read-only; writes nothing back to the adata.

Run, then paste the printed report (or share the PNG) so we can choose how to
scale the axes for the figure.

Uses the confirmed constants from Generate_Figure3_A3_SBS2_CNV_Scatter.py.
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================== CONFIG ==============================
PROJECT_ROOT     = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")
OUTPUT_DIR       = os.path.join(PROJECT_ROOT, "data/FIG_3/diagnostics")

ANNOT_COL  = "final_annotation"
SOURCE_COL = "source_name"
CNV_COL    = "cnv_score"
A3A_GENE   = "APOBEC3A"
A3B_GENE   = "APOBEC3B"

TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"  # GEO typo intentional
# ====================================================================


def log(msg=""):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def load_basal():
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    b = adata[adata.obs[ANNOT_COL].astype(str) == TARGET_CELL_TYPE].copy()
    log(f"Basal cells: {b.n_obs:,}")

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
            b.obs[g] = 0.0

    b.obs[CNV_COL] = b.obs[CNV_COL].astype(float)
    b.obs["tissue_grp"] = np.where(
        b.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, "normal", "tumor")
    return b.obs[["SBS2", CNV_COL, A3A_GENE, A3B_GENE, "tissue_grp"]].copy()


def profile(name, vals):
    v = np.asarray(vals, dtype=float)
    v = v[~np.isnan(v)]
    n = v.size
    nz = int((v == 0).sum())
    uniq = np.unique(v)
    diffs = np.diff(uniq)
    min_gap = float(diffs[diffs > 0].min()) if np.any(diffs > 0) else float("nan")
    pct = np.percentile(v, [1, 5, 25, 50, 75, 95, 99]) if n else [np.nan] * 7

    print(f"\n  --- {name} ---")
    print(f"    n              : {n:,}")
    print(f"    zeros          : {nz:,} ({100*nz/max(n,1):.1f}%)")
    print(f"    unique values  : {uniq.size:,}")
    print(f"    min / max      : {v.min():.5g} / {v.max():.5g}")
    print(f"    mean / median  : {v.mean():.5g} / {np.median(v):.5g}")
    print(f"    std            : {v.std():.5g}")
    print(f"    min nonzero gap: {min_gap:.5g}   <- regular spacing => quantized")
    print(f"    pctl 1/5/25/50/75/95/99:")
    print(f"      {', '.join(f'{p:.4g}' for p in pct)}")
    # most common values (rounded for grouping)
    s = pd.Series(v).round(4).value_counts().head(12)
    print(f"    top values (rounded 4dp): count")
    for val, cnt in s.items():
        print(f"      {val:>10}: {cnt:,}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    obs = load_basal()
    tum = obs[obs["tissue_grp"] == "tumor"]
    nor = obs[obs["tissue_grp"] == "normal"]

    cols = ["SBS2", CNV_COL, A3A_GENE, A3B_GENE]

    print("\n" + "=" * 70)
    print("ALL BASAL")
    print("=" * 70)
    for c in cols:
        profile(c, obs[c].values)

    print("\n" + "=" * 70)
    print(f"TUMOR BASAL (n={len(tum):,})")
    print("=" * 70)
    for c in cols:
        profile(c, tum[c].values)

    print("\n" + "=" * 70)
    print(f"NORMAL-ADJACENT BASAL (n={len(nor):,})")
    print("=" * 70)
    for c in cols:
        profile(c, nor[c].values)

    # histograms on tumor basal (the dense group where banding appears)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, c in zip(axes.ravel(), cols):
        ax.hist(tum[c].values, bins=80, color="#4682b4", edgecolor="none")
        ax.set_title(f"{c} (tumor basal)")
        ax.set_xlabel(c)
        ax.set_ylabel("cells")
    fig.tight_layout()
    hist_png = os.path.join(OUTPUT_DIR, "fig3_axis_histograms_tumor.png")
    fig.savefig(hist_png, dpi=200)
    plt.close(fig)

    # CNV unique-value spacing plot (shows quantization directly)
    cnv_u = np.unique(tum[CNV_COL].values)
    fig2, ax2 = plt.subplots(figsize=(12, 4))
    ax2.plot(np.arange(cnv_u.size), cnv_u, marker="o", ms=3, lw=0.8, color="#ed6a5a")
    ax2.set_title(f"cnv_score sorted unique values (n_unique={cnv_u.size})")
    ax2.set_xlabel("rank of unique value")
    ax2.set_ylabel("cnv_score")
    fig2.tight_layout()
    cnv_png = os.path.join(OUTPUT_DIR, "fig3_cnv_unique_spacing.png")
    fig2.savefig(cnv_png, dpi=200)
    plt.close(fig2)

    log(f"\nSaved:\n  {hist_png}\n  {cnv_png}")


if __name__ == "__main__":
    main()
