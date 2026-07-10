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

NEW (this revision):
  Adds a "Final canonical group" column (SBS2-HIGH / CNV-HIGH / NORMAL /
  unassigned) sourced from three_group_assignments.tsv, and a CANONICAL GROUP
  MAPPING AUDIT block that verifies the join in both directions before handoff.

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

# canonical group assignments (the same file Step00B / Figure 4 write & use)
GROUPS_FILE = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/"
               "01_group_selection/three_group_assignments.tsv")
GROUP_LABEL_MAP = {"SBS2_HIGH": "SBS2-HIGH", "CNV_HIGH": "CNV-HIGH", "NORMAL": "NORMAL"}
N_EXPECTED = 546  # cells per canonical group

# ---- DECISIONS (confirm against the DIAGNOSTIC BLOCK printed below) --------
SBS_QUESTION_COL = "sig_SBS13"   # the boss's open "SBS ?" column (APOBEC C>G)
# HPV status is the Phase3 L-method call: raw_HPV16 >= 8 UMI (same as Fig 6).
# NOTE: master['HPV16_raw_positive'] is a >0 call (threshold 1), NOT this; do
# not use it. Ambiguous cells (1-7 reads) fold into -1 in the binary column;
# the raw UMI count is exported alongside so they can be re-tiered.
HPV_READS_COL  = "raw_HPV16"
HPV_THRESHOLD  = 8
# Canonical Fig 6 gate is the AND of the UMI threshold and >0 genome reads.
# This master table stores the genome-reads field as 'HPV16_reads'
# (the Fig 6 'total_hpv16_genome_reads' quantity under a shorter name); the
# per-group cross-check below against Fig 6 (197/446/8) confirms it's correct.
HPV_GENOME_COL = "HPV16_reads"
# Manuscript (Fig 6) HPV+ counts per canonical group, for the audit cross-check.
# These are a reference to confirm against the paper, NOT a hard pass/fail.
EXPECTED_HPV_POS = {"SBS2-HIGH": 197, "CNV-HIGH": 446, "NORMAL": 8}
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

# HPV gate: canonical Fig 6 definition is raw_HPV16 >= 8 UMI AND
# total_hpv16_genome_reads > 0. Resolve the genome-reads column, then apply both.
print("\nHPV gating (canonical Fig 6 definition: raw>=8 AND genome>0):")
raw_hpv = pd.to_numeric(df[HPV_READS_COL], errors="coerce").fillna(0)

genome_candidates = [c for c in df.columns if "hpv" in c.lower()
                     and ("genome" in c.lower() or "read" in c.lower())]
if HPV_GENOME_COL in df.columns:
    genome_reads = pd.to_numeric(df[HPV_GENOME_COL], errors="coerce").fillna(0)
    print(f"   genome-reads column: '{HPV_GENOME_COL}' (found)")
else:
    print(f"   !! genome-reads column '{HPV_GENOME_COL}' NOT FOUND")
    print(f"   candidate columns matching hpv+genome/read: {genome_candidates}")
    raise ValueError(
        f"Set HPV_GENOME_COL to the correct column name. Candidates: {genome_candidates}")

hpv_pos_bool = (raw_hpv >= HPV_THRESHOLD) & (genome_reads > 0)
n_raw_only = int(((raw_hpv >= HPV_THRESHOLD) & (genome_reads <= 0)).sum())

print(f"   negative (raw 0)        : {int((raw_hpv == 0).sum())}")
print(f"   ambiguous (raw 1-7)     : {int(((raw_hpv >= 1) & (raw_hpv < HPV_THRESHOLD)).sum())}"
      f"   <- folded into -1")
print(f"   raw >= 8 (UMI only)      : {int((raw_hpv >= HPV_THRESHOLD).sum())}")
print(f"   canonical HPV+ (raw>=8 AND genome>0): {int(hpv_pos_bool.sum())}")
print(f"   demoted (raw>=8 but genome==0)      : {n_raw_only}   <- now -1, not +1")
if "HPV16_raw_positive" in df.columns:
    _wrong = int((pd.to_numeric(df["HPV16_raw_positive"], errors="coerce") > 0.5).sum())
    print(f"   [for contrast] HPV16_raw_positive (>0 call) would flag {_wrong} "
          f"-> do NOT use")

# donor / tissue
print("\ntissue type:")
print(df["tissue type"].value_counts().to_string())
print("\ndistinct donors:", df["subject id"].nunique())

# ---- CANONICAL GROUP LABELS + MAPPING AUDIT -------------------------------
assign = pd.read_csv(GROUPS_FILE, sep="\t")
if "cell_barcode" not in assign.columns:
    assign = assign.rename(columns={assign.columns[0]: "cell_barcode"})

bc_to_label = {bc: GROUP_LABEL_MAP.get(g, g)
               for bc, g in zip(assign["cell_barcode"], assign["group"])}
final_group = df["cell_barcode"].map(bc_to_label).fillna("unassigned")

print("\n" + "=" * 78)
print("CANONICAL GROUP MAPPING AUDIT  (verify before handoff)")
print("=" * 78)
print("source:", GROUPS_FILE)

file_counts = assign["group"].value_counts()
n_dup = int(assign["cell_barcode"].duplicated().sum())
print("\nassignments file:")
for g in ["SBS2_HIGH", "CNV_HIGH", "NORMAL"]:
    print(f"   {g:10s}: {int(file_counts.get(g, 0))}")
print(f"   total rows: {len(assign)}   duplicate barcodes in file: {n_dup}")

mapped = final_group.value_counts()
print("\nmapped onto master table (the new output column):")
for lab in ["SBS2-HIGH", "CNV-HIGH", "NORMAL", "unassigned"]:
    print(f"   {lab:12s}: {int(mapped.get(lab, 0))}")

master_bc = set(df["cell_barcode"])
missing = {}
print("\ncanonical barcodes present in master table:")
for g, pretty in GROUP_LABEL_MAP.items():
    grp_bc = set(assign.loc[assign["group"] == g, "cell_barcode"])
    miss = grp_bc - master_bc
    missing[pretty] = miss
    print(f"   {pretty:12s}: {len(grp_bc) - len(miss)}/{len(grp_bc)} found")
    if miss:
        print(f"       MISSING {len(miss)} e.g. {list(miss)[:5]}")

# explicit pass/fail checks
checks = [("no duplicate barcodes in file", n_dup == 0)]
for g, pretty in GROUP_LABEL_MAP.items():
    fc = int(file_counts.get(g, 0))
    checks.append((f"{pretty}: {N_EXPECTED} rows in file", fc == N_EXPECTED))
    checks.append((f"{pretty}: all barcodes found in master", len(missing[pretty]) == 0))
    checks.append((f"{pretty}: {N_EXPECTED} rows in output column",
                   int(mapped.get(pretty, 0)) == N_EXPECTED))

print("\nRESULT:")
for name, ok in checks:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}")
overall = all(ok for _, ok in checks)
print(f"\n   OVERALL: {'PASS -- mapping verified' if overall else 'FAIL -- see above'}")

# ---- PER-GROUP HPV BREAKDOWN (canonical gate; cross-check against Fig 6) ----
print("\n" + "=" * 78)
print("PER-GROUP HPV BREAKDOWN  (canonical gate; compare HPV+ to Fig 6)")
print("=" * 78)
for pretty in ["SBS2-HIGH", "CNV-HIGH", "NORMAL"]:
    gmask = (final_group == pretty).values
    grp_raw = raw_hpv[gmask]
    grp_pos = int(hpv_pos_bool[gmask].sum())
    neg = int((grp_raw == 0).sum())
    amb = int(((grp_raw >= 1) & (grp_raw < HPV_THRESHOLD)).sum())
    rawpos = int((grp_raw >= HPV_THRESHOLD).sum())
    exp = EXPECTED_HPV_POS.get(pretty)
    flag = ("  [MATCH Fig 6]" if grp_pos == exp
            else f"  [CHECK vs Fig 6={exp}]") if exp is not None else ""
    print(f"   {pretty:10s}: neg(0)={neg:3d}  amb(1-7)={amb:3d}  "
          f"raw>=8={rawpos:3d}  canonical HPV+={grp_pos:3d}{flag}")
print("   Fig 6 reference HPV+ counts: SBS2-HIGH 197 / CNV-HIGH 446 / NORMAL 8 "
      "(confirm against manuscript)")

# ---- BUILD THE TABLE ------------------------------------------------------
# raw_hpv and hpv_pos_bool (canonical gate) were computed in the diagnostic block.
hpv_status = np.where(hpv_pos_bool, 1, -1)

# clean source label without touching any GEO metadata column
src = df["tissue type"].astype(str).str.lower().map(
    {"tumor": "tumor", "normal": "normal-adjacent"}).fillna(df["tissue type"])

q_label = SBS_QUESTION_COL.replace("sig_", "")  # e.g. SBS13

out = pd.DataFrame({
    "Cell barcode": df["cell_barcode"],
    "Final canonical group": final_group,
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

# show a few cells from each canonical group so the PI can spot-check the join
print("\nSample rows per canonical group:")
for pretty in ["SBS2-HIGH", "CNV-HIGH", "NORMAL"]:
    sample = out.loc[out["Final canonical group"] == pretty,
                     ["Cell barcode", "Final canonical group", "SBS2 weight",
                      "CNV score", "HPV status (1=pos, -1=neg)"]].head(3)
    print(f"\n  {pretty}:")
    print(sample.to_string(index=False))
