#!/usr/bin/env python3
"""
Diagnostic_Fig3_Lean_Decompose.py
=================================
Pull apart the lean result. For A3A, A3B, and the A3A dominance fraction
A3A/(A3A+A3B), report Spearman rho against SBS2, cnv_score (alone),
CytoTRACE2 (alone), and the cnv_stem blend, for tumor and normal, computed
both on expressing cells only and on all basal cells.

Goal: see (1) whether stemness, not CNV, drives the strong negative on the
CNV+Stem pole, and (2) how much conditioning on A3>0 removed the co-occurrence
signal, and (3) whether the A3A fraction differentiates the way absolute levels don't.
"""
import os
from datetime import datetime
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr

PROJECT_ROOT     = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")
ANNOT_COL, SOURCE_COL = "final_annotation", "source_name"
CNV_COL, CYTO_COL = "cnv_score", "CytoTRACE2_Score"
A3A_GENE, A3B_GENE = "APOBEC3A", "APOBEC3B"
TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE = "normal tissue adjucent to head and neck squamous cell carcinoma"

def log(m=""): print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)
def mm(x):
    x = np.asarray(x, float); lo, hi = np.nanmin(x), np.nanmax(x)
    return (x-lo)/(hi-lo) if hi > lo else np.zeros_like(x)

def rr(a, b):
    if len(a) < 5: return np.nan
    r, _ = spearmanr(a, b); return r

adata = sc.read_h5ad(ADATA_FINAL_PATH)
b = adata[adata.obs[ANNOT_COL].astype(str) == TARGET_CELL_TYPE].copy()
log(f"Basal: {b.n_obs:,}")
w = pd.read_csv(WEIGHTS_PATH, sep="\t", index_col=0); m = w.loc["SBS2"].to_dict()
b.obs["SBS2"] = b.obs_names.map(lambda x: m.get(x, 0.0)).astype(float)
genes = list(b.var_names)
for g in (A3A_GENE, A3B_GENE):
    v = b.X[:, genes.index(g)]; v = v.toarray().flatten() if hasattr(v, "toarray") else np.asarray(v).flatten()
    b.obs[g] = v
b.obs[CNV_COL] = b.obs[CNV_COL].astype(float); b.obs[CYTO_COL] = b.obs[CYTO_COL].astype(float)
b.obs["cnv_stem"] = 0.5*mm(b.obs[CNV_COL].values) + 0.5*mm(b.obs[CYTO_COL].values)
b.obs["frac"] = b.obs[A3A_GENE] / (b.obs[A3A_GENE] + b.obs[A3B_GENE]).replace(0, np.nan)
b.obs["tissue_grp"] = np.where(b.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, "normal", "tumor")
obs = b.obs

def block(title, df):
    print(f"\n=== {title} (n={len(df):,}) ===")
    print(f"{'var':>6} | {'SBS2':>8} {'cnv':>8} {'CytoT2':>8} {'cnv+stem':>9}")
    for var, col, mask in [("A3A", A3A_GENE, df[A3A_GENE] > 0),
                           ("A3B", A3B_GENE, df[A3B_GENE] > 0),
                           ("frac", "frac", (df[A3A_GENE] + df[A3B_GENE]) > 0)]:
        d = df[mask]
        print(f"{var:>6} | {rr(d[col], d['SBS2']):>8.3f} {rr(d[col], d[CNV_COL]):>8.3f} "
              f"{rr(d[col], d[CYTO_COL]):>8.3f} {rr(d[col], d['cnv_stem']):>9.3f}")

for grp in ["tumor", "normal"]:
    g = obs[obs["tissue_grp"] == grp]
    block(f"{grp.upper()}  (expressing cells per row)", g)

print("\n--- ALL BASAL, no expression conditioning (co-occurrence included) ---")
for grp in ["tumor", "normal"]:
    g = obs[obs["tissue_grp"] == grp]
    print(f"\n=== {grp.upper()} all basal (n={len(g):,}) ===")
    print(f"{'var':>6} | {'SBS2':>8} {'cnv':>8} {'CytoT2':>8} {'cnv+stem':>9}")
    for lab, var in [("A3A", A3A_GENE), ("A3B", A3B_GENE)]:
        print(f"{lab:>6} | {rr(g[var], g['SBS2']):>8.3f} {rr(g[var], g[CNV_COL]):>8.3f} "
              f"{rr(g[var], g[CYTO_COL]):>8.3f} {rr(g[var], g['cnv_stem']):>9.3f}")
