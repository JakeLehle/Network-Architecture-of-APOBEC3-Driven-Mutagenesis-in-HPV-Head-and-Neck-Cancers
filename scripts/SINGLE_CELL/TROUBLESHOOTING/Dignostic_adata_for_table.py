#!/usr/bin/env python
"""
Diagnostic: inventory adata_final.h5ad for the basal-cell export table.

Run in the NETWORK env on titan. Reads only; writes nothing to the object.
Goal: confirm exactly where each requested quantity lives (and under what
column name / layer) BEFORE we write the extraction script.

Requested table columns (per boss):
  SBS2 / SBS5 / SBS40 / SBS? weight, CNV score, stemness score,
  A3A/A3B/A3C/A3D/A3F/A3G/A3H expression, HPV status (8-UMI threshold),
  donor ID, source (tumor vs normal-adjacent), cell type.
Filter: basal cells in annotation_final.
"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc

# ---- EDIT if your path differs -------------------------------------------
ADATA_PATH = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/adata_final.h5ad"
# external HPV table (raw UMIs) referenced in the SC pipeline
HPV_TABLE_GLOB = "basal_cell_master_table_with_raw_HPV"
# --------------------------------------------------------------------------

pd.set_option("display.max_rows", 200)
pd.set_option("display.width", 160)

print("=" * 78)
print("LOADING:", ADATA_PATH)
adata = sc.read_h5ad(ADATA_PATH)
print(adata)
print("\nshape (cells x genes):", adata.shape)

# ---- 1. annotation_final: confirm exact basal label -----------------------
print("\n" + "=" * 78)
print("1. annotation_final value_counts (confirm exact basal label string)")
print("=" * 78)
ann_col = "annotation_final"
if ann_col in adata.obs.columns:
    print(adata.obs[ann_col].value_counts(dropna=False))
else:
    print(f"!! '{ann_col}' NOT in obs. Candidate annotation-like columns:")
    print([c for c in adata.obs.columns if "annot" in c.lower() or "cell_type" in c.lower()
           or "celltype" in c.lower()])

# ---- 2. full obs column inventory with dtype ------------------------------
print("\n" + "=" * 78)
print("2. ALL obs columns (name | dtype | n_unique | example)")
print("=" * 78)
for c in adata.obs.columns:
    s = adata.obs[c]
    try:
        ex = s.dropna().iloc[0]
    except Exception:
        ex = "NA"
    print(f"  {c:42s} | {str(s.dtype):12s} | nuniq={s.nunique():>6} | ex={ex}")

# ---- 3. SBS / signature weight columns ------------------------------------
print("\n" + "=" * 78)
print("3. Signature-weight columns (grep SBS / weight / signature / nmf)")
print("=" * 78)
sig_cols = [c for c in adata.obs.columns
            if any(k in c.lower() for k in ["sbs", "weight", "signature", "nmf"])]
print("matched:", sig_cols)
for c in sig_cols:
    s = pd.to_numeric(adata.obs[c], errors="coerce")
    print(f"  {c:42s} min={np.nanmin(s):.4g} med={np.nanmedian(s):.4g} "
          f"max={np.nanmax(s):.4g} nan={s.isna().sum()}")

# ---- 4. CNV / stemness columns --------------------------------------------
print("\n" + "=" * 78)
print("4. CNV / stemness candidate columns")
print("=" * 78)
for key in ["cnv", "stem"]:
    hits = [c for c in adata.obs.columns if key in c.lower()]
    print(f"  '{key}' ->", hits)
    for c in hits:
        s = pd.to_numeric(adata.obs[c], errors="coerce")
        print(f"      {c}: min={np.nanmin(s):.4g} med={np.nanmedian(s):.4g} max={np.nanmax(s):.4g}")

# ---- 5. donor / source columns --------------------------------------------
print("\n" + "=" * 78)
print("5. donor ID / source (tumor vs normal-adjacent) candidates")
print("=" * 78)
for key in ["donor", "patient", "sample", "subject", "source", "tissue",
            "tumor", "normal", "adjacent", "adjucent", "condition", "status"]:
    hits = [c for c in adata.obs.columns if key in c.lower()]
    if hits:
        print(f"  '{key}' ->", hits)
        for c in hits:
            print(f"      {c} value_counts:")
            print(adata.obs[c].value_counts(dropna=False).to_string().replace("\n", "\n        "))

# ---- 6. HPV: in-object UMI column AND external table check ----------------
print("\n" + "=" * 78)
print("6. HPV: look for raw-UMI column in obs (8-UMI threshold source)")
print("=" * 78)
hpv_cols = [c for c in adata.obs.columns
            if any(k in c.lower() for k in ["hpv", "umi", "viral", "kraken", "e7", "e6"])]
print("  obs HPV/UMI candidates:", hpv_cols)
for c in hpv_cols:
    s = pd.to_numeric(adata.obs[c], errors="coerce")
    if s.notna().any():
        print(f"      {c}: min={np.nanmin(s):.4g} med={np.nanmedian(s):.4g} "
              f"max={np.nanmax(s):.4g} ; cells>=8 = {(s>=8).sum()}")

import glob, os
print("\n  external HPV master table search:")
search_roots = ["/master/jlehle/WORKING/2026_NMF_PAPER",
                os.path.dirname(ADATA_PATH)]
found = []
for root in search_roots:
    found += glob.glob(os.path.join(root, "**", f"*{HPV_TABLE_GLOB}*"), recursive=True)
found = sorted(set(found))
print("   ", found if found else "  (none found under search roots)")
if found:
    t = pd.read_csv(found[0], sep=None, engine="python", nrows=5)
    print("    columns:", list(t.columns))
    # barcode-overlap sanity check
    full = pd.read_csv(found[0], sep=None, engine="python")
    bc_col = next((c for c in full.columns
                   if c.lower() in ("barcode", "cell", "cell_id", "index")), full.columns[0])
    overlap = len(set(full[bc_col]).intersection(set(adata.obs_names)))
    print(f"    barcode col guess: {bc_col} | overlap with adata.obs_names: "
          f"{overlap}/{adata.n_obs}")

# ---- 7. APOBEC3 genes in var ----------------------------------------------
print("\n" + "=" * 78)
print("7. APOBEC3 family genes present in var_names")
print("=" * 78)
a3 = ["APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
      "APOBEC3F", "APOBEC3G", "APOBEC3H"]
vn = set(map(str, adata.var_names))
for g in a3:
    print(f"  {g:10s} present={g in vn}")
print("  any APOBEC* in var:", [g for g in adata.var_names if str(g).upper().startswith("APOBEC")])

# ---- 8. expression units: X state + layers --------------------------------
print("\n" + "=" * 78)
print("8. Expression representation (pick units for A3 columns deliberately)")
print("=" * 78)
print("  layers:", list(adata.layers.keys()))
Xs = adata.X[:50].toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X[:50])
print(f"  X dtype={adata.X.dtype} min={Xs.min():.4g} max={Xs.max():.4g} "
      f"looks_integer={np.allclose(Xs, np.round(Xs))}")
for ln in adata.layers.keys():
    L = adata.layers[ln]
    Ls = L[:50].toarray() if hasattr(L, "toarray") else np.asarray(L[:50])
    print(f"  layer[{ln}] min={Ls.min():.4g} max={Ls.max():.4g} "
          f"looks_integer={np.allclose(Ls, np.round(Ls))}")

# ---- 9. obsm / uns quick peek ---------------------------------------------
print("\n" + "=" * 78)
print("9. obsm / uns keys (in case scores live here)")
print("=" * 78)
print("  obsm:", list(adata.obsm.keys()))
print("  uns :", list(adata.uns.keys()))

print("\n" + "=" * 78)
print("DONE. Paste this whole output back and I'll write the extraction script.")
print("=" * 78)
