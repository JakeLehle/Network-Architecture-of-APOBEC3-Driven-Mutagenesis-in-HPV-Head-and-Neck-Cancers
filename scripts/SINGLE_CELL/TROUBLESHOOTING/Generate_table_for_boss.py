#!/usr/bin/env python
"""
Generate_Basal_Cell_Table_for_Boss.py

Build the per-cell basal table requested by the PI, reshaped to his column layout.

WHY NOT adata_final directly:
  The adata_final object does NOT carry per-cell SBS signature weights, and its
  HPV field ('hpv state') is all-NaN. Both required quantities live in the basal
  master table, which is ALREADY subset to basal cells (n=52,126) and already
  carries cnv_score, CytoTRACE2 stemness, APOBEC3A-H expression, donor, and tissue.
  That table was derived from adata_final, so it is the faithful and complete spine.

DECISIONS surfaced at top (flip if the printed diagnostic block disagrees):
  SBS_QUESTION_COL : which signature fills the boss's open "SBS ?" column.
                     Default sig_SBS13 (C>G APOBEC partner to SBS2's C>T).
  HPV_POS_COL/READS: which positive call is the 8-UMI threshold. Default the
                     "raw" pair; the script VERIFIES min(reads | positive)==8.

A3 expression units are whatever the master table stored; reported in diagnostics
so the units are a deliberate choice, not an assumption.

Run in NETWORK env on titan.
"""

import os
import numpy as np
import pandas as pd

# ---- paths ----------------------------------------------------------------
MASTER = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/"
          "01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
OUTDIR = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/01_raw_hpv16_counts")
OUTBASE = "basal_cell_table_for_PI"

# ---- DECISIONS (confirm against the DIAGNOSTIC BLOCK printed below) --------
SBS_QUESTION_COL = "sig_SBS13"   # the boss's open "SBS ?" column (APOBEC C>G)
# HPV status is the Phase3 L-method call: raw_HPV16 >= 8 UMI (same as Fig 6).
# NOTE: master['HPV16_raw_positive'] is a >0 call (threshold 1), NOT this; do
# not use it. Ambiguous cells (1-7 reads) fold into -1 in the binary column;
# the raw UMI count is exported alongside so they can be re-tiered.
HPV_READS_COL  = "raw_HPV16"
HPV_THRESHOLD  = 8
# ---------------------------------------------------------------------------

pd.set_option("display.width", 160)

df = pd.read_csv(MASTER, sep="\t")
df = df.rename(columns={"Unnamed: 0": "cell_barcode"})


def _as_bool(series):
    """Parse a positive-call column that may be bool / int / 'True'/'False'."""
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])


print("=" * 78)
print("DIAGNOSTIC BLOCK  (read this before handing the table to the PI)")
print("=" * 78)
print("source:", MASTER)
print("rows (basal cells):", len(df))

# basal confirmation
if "Final_cancer_cell_status" in df.columns:
    print("\nFinal_cancer_cell_status:")
    print(df["Final_cancer_cell_status"].value_counts().to_string())

# full SBS menu + the four we will export
sbs_cols = [c for c in df.columns if c.startswith("sig_SBS")]
print("\nAll available SBS signatures:", sbs_cols)
print("Exporting SBS2 / SBS5 / SBS40a / %s :" % SBS_QUESTION_COL.replace("sig_", ""))
for c in ["sig_SBS2", "sig_SBS5", "sig_SBS40a", SBS_QUESTION_COL]:
    if c in df.columns:
        s = pd.to_numeric(df[c], errors="coerce")
        print(f"   {c:12s} min={s.min():.4g}  med={s.median():.4g}  "
              f"max={s.max():.4g}  nonzero={(s > 0).sum()}")
    else:
        print(f"   {c:12s} !! NOT FOUND in table")

# A3 units check
print("\nA3 expression ranges (decide if these units are what the PI wants):")
a3 = ["APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
      "APOBEC3F", "APOBEC3G", "APOBEC3H"]
for g in a3:
    s = pd.to_numeric(df[g], errors="coerce").dropna()
    looks_int = np.allclose(s, np.round(s)) if len(s) else False
    print(f"   {g}: min={s.min():.4g}  max={s.max():.4g}  "
          f"looks_integer(raw_counts)={looks_int}")

# HPV threshold: the correct 8-UMI call is raw_HPV16 >= 8 (Phase3 L-method).
print("\nHPV tiering from raw_HPV16 (correct 8-UMI call, NOT HPV16_raw_positive):")
_raw = pd.to_numeric(df[HPV_READS_COL], errors="coerce").fillna(0)
print(f"   negative (0)    : {int((_raw == 0).sum())}")
print(f"   ambiguous (1-7) : {int(((_raw >= 1) & (_raw < HPV_THRESHOLD)).sum())}"
      f"   <- folded into -1")
print(f"   positive (>=8)  : {int((_raw >= HPV_THRESHOLD).sum())}")
if "HPV16_raw_positive" in df.columns:
    _wrong = int((pd.to_numeric(df["HPV16_raw_positive"], errors="coerce") > 0.5).sum())
    print(f"   [for contrast] HPV16_raw_positive (>0 call) would flag {_wrong} "
          f"-> do NOT use")

# donor / tissue
print("\ntissue type:")
print(df["tissue type"].value_counts().to_string())
print("\ndistinct donors:", df["subject id"].nunique())

# ---- BUILD THE TABLE ------------------------------------------------------
raw_hpv = pd.to_numeric(df[HPV_READS_COL], errors="coerce").fillna(0)
hpv_status = np.where(raw_hpv >= HPV_THRESHOLD, 1, -1)

# clean source label without touching any GEO metadata column
src = df["tissue type"].astype(str).str.lower().map(
    {"tumor": "tumor", "normal": "normal-adjacent"}).fillna(df["tissue type"])

q_label = SBS_QUESTION_COL.replace("sig_", "")  # e.g. SBS13

out = pd.DataFrame({
    "Cell barcode": df["cell_barcode"],
    "SBS2 weight": pd.to_numeric(df["sig_SBS2"], errors="coerce"),
    "SBS5 weight": pd.to_numeric(df["sig_SBS5"], errors="coerce"),
    "SBS40 weight": pd.to_numeric(df["sig_SBS40a"], errors="coerce"),
    f"{q_label} weight": pd.to_numeric(df[SBS_QUESTION_COL], errors="coerce"),
    "CNV score": pd.to_numeric(df["cnv_score"], errors="coerce"),
    "Stemness score (CytoTRACE2)": pd.to_numeric(df["CytoTRACE2_Score"], errors="coerce"),
    "A3A expression": pd.to_numeric(df["APOBEC3A"], errors="coerce"),
    "A3B expression": pd.to_numeric(df["APOBEC3B"], errors="coerce"),
    "A3C expression": pd.to_numeric(df["APOBEC3C"], errors="coerce"),
    "A3D expression": pd.to_numeric(df["APOBEC3D"], errors="coerce"),
    "A3F expression": pd.to_numeric(df["APOBEC3F"], errors="coerce"),
    "A3G expression": pd.to_numeric(df["APOBEC3G"], errors="coerce"),
    "A3H expression": pd.to_numeric(df["APOBEC3H"], errors="coerce"),
    "HPV status (1=pos, -1=neg)": hpv_status,
    "HPV16 raw UMI count": raw_hpv.astype(int).values,
    "Donor ID": df["subject id"],
    "Source (tumor vs normal-adjacent)": src,
    "Cell type": "basal cell",
})

os.makedirs(OUTDIR, exist_ok=True)
tsv = os.path.join(OUTDIR, OUTBASE + ".tsv")
out.to_csv(tsv, sep="\t", index=False)
print("\n" + "=" * 78)
print("WROTE:", tsv, " shape:", out.shape)
try:
    xlsx = os.path.join(OUTDIR, OUTBASE + ".xlsx")
    out.to_excel(xlsx, index=False)
    print("WROTE:", xlsx)
except Exception as e:
    print("xlsx skipped (openpyxl not installed?):", e)
print("=" * 78)
print("\nFirst 5 rows:")
print(out.head().to_string())
