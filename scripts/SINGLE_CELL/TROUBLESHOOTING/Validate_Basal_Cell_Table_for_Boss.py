#!/usr/bin/env python
"""
Validate_Basal_Cell_Table_for_Boss.py

Sanity-check the basal export table BEFORE it goes to the PI. Three blocks:

  A. INTERNAL CONSISTENCY of the written table (dupes, NaNs, counts, basal-only).
  B. HPV THRESHOLD AUDIT. Proves that master['HPV16_raw_positive'] is a >0 call,
     NOT the 8-UMI call, and quantifies how many cells that mis-call would flip.
     The correct status is raw_HPV16 >= 8 (Phase3 L-method threshold, Fig 6).
  C. BIOLOGY. Reproduces the Figure 3 dominance result the PI will eyeball:
       A3A_fraction = A3A / (A3A + A3B + 0.01), among cells expressing EITHER.
       Expected (repo, Generate_Figure3_A3_Dominance_Slider.py):
         A3A_fraction vs CNV : ~ -0.44 in TUMOR   (high CNV -> A3B-dominant)
                               ~ +0.07 in NORMAL  (tumor-specific)
         A3A_fraction vs stemness : ~ -0.5 in BOTH (NOT tumor-specific)
     and tabulates CNV / stemness by A3A- vs A3B-dominant calls.

Run in NETWORK env on titan.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu

MASTER = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/"
          "01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
TABLE = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/"
         "01_raw_hpv16_counts/basal_cell_table_for_PI.tsv")
HPV_THRESHOLD = 8          # Phase3 L-method + SC008 ambient control
DOM_MARGIN = 0.0           # A3A-dominant: frac > 0.5+margin; A3B: frac < 0.5-margin

pd.set_option("display.width", 170)


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


master = pd.read_csv(MASTER, sep="\t").rename(columns={"Unnamed: 0": "cell_barcode"})
table = pd.read_csv(TABLE, sep="\t") if os.path.exists(TABLE) else None

# =====================================================================
# A. INTERNAL CONSISTENCY OF THE WRITTEN TABLE
# =====================================================================
hr("A. INTERNAL CONSISTENCY (written output table)")
if table is None:
    print("  !! output table not found; run the generator first:", TABLE)
else:
    print("  rows:", len(table), "| columns:", len(table.columns))
    dup = table["Cell barcode"].duplicated().sum()
    print(f"  duplicate barcodes: {dup}  (must be 0)")
    print("  NaNs per column (any > 0 worth explaining to the PI):")
    nn = table.isna().sum()
    for c in table.columns:
        if nn[c] > 0:
            print(f"     {c:38s} {nn[c]}")
    if nn.sum() == 0:
        print("     none")
    print("\n  HPV status as written:")
    print(table["HPV status (1=pos, -1=neg)"].value_counts().to_string())
    print("\n  tissue split:")
    print(table["Source (tumor vs normal-adjacent)"].value_counts().to_string())
    print("  distinct donors:", table["Donor ID"].nunique())
    # basal-only check (the 'Cell type' constant is only honest if true)
    print("\n  master barcode overlap with table:",
          len(set(master["cell_barcode"]).intersection(table["Cell barcode"])),
          "/", len(table))

# =====================================================================
# B. HPV THRESHOLD AUDIT  (the headline catch)
# =====================================================================
hr("B. HPV THRESHOLD AUDIT")
raw = pd.to_numeric(master["raw_HPV16"], errors="coerce").fillna(0)

# what the pre-baked boolean actually is
if "HPV16_raw_positive" in master.columns:
    rawpos = master["HPV16_raw_positive"].astype(float) > 0.5
    print(f"  master['HPV16_raw_positive']: n_pos={int(rawpos.sum())}, "
          f"min raw among 'positive' = {raw[rawpos].min():.0f}  "
          f"(==1 confirms it is a >0 call, NOT 8-UMI)")

neg = int((raw == 0).sum())
amb = int(((raw >= 1) & (raw < HPV_THRESHOLD)).sum())
pos = int((raw >= HPV_THRESHOLD).sum())
print(f"\n  CORRECT tiering from raw_HPV16:")
print(f"     negative (0)      : {neg}")
print(f"     ambiguous (1-7)   : {amb}   <- folded into -1 in the binary table")
print(f"     positive (>=8)    : {pos}")
print(f"\n  Mis-call magnitude if >0 were used instead of >=8:")
print(f"     cells flagged +ve by >0 but NOT by >=8 (the over-call): "
      f"{int((raw >= 1).sum()) - pos}")

# cross-tab correct HPV vs tissue (positives should sit in tumor)
hpv_bin = np.where(raw >= HPV_THRESHOLD, 1, -1)
ct = pd.crosstab(master["tissue type"], pd.Series(hpv_bin, index=master.index,
                                                  name="HPV(>=8)"))
print("\n  CORRECT HPV(>=8) vs tissue (positives should concentrate in tumor):")
print(ct.to_string())

# does the written table match the correct call?
if table is not None:
    written = table.set_index("Cell barcode")["HPV status (1=pos, -1=neg)"]
    correct = pd.Series(hpv_bin, index=master["cell_barcode"])
    common = written.index.intersection(correct.index)
    mismatch = int((written.loc[common] != correct.loc[common]).sum())
    print(f"\n  written HPV column vs correct >=8 call: {mismatch} mismatches "
          f"of {len(common)}  (0 == table already fixed)")

# =====================================================================
# C. BIOLOGY: A3 DOMINANCE vs CNV  (reproduce Fig 3)
# =====================================================================
hr("C. BIOLOGY: A3A/A3B dominance vs CNV (vs Fig 3 expectation)")
m = master.copy()
for col in ["APOBEC3A", "APOBEC3B", "cnv_score", "CytoTRACE2_Score", "sig_SBS2"]:
    m[col] = pd.to_numeric(m[col], errors="coerce")

# units sanity: are A3 columns lognorm (what the figures use) or raw counts?
for g in ["APOBEC3A", "APOBEC3B"]:
    s = m[g].dropna()
    print(f"  {g}: min={s.min():.3g} max={s.max():.3g} "
          f"looks_integer(raw)={np.allclose(s, np.round(s))}")

m["A3A_fraction"] = m["APOBEC3A"] / (m["APOBEC3A"] + m["APOBEC3B"] + 0.01)
m["A3B_fraction"] = m["APOBEC3B"] / (m["APOBEC3A"] + m["APOBEC3B"] + 0.01)
expr_either = m[(m["APOBEC3A"] + m["APOBEC3B"]) > 0].copy()
print(f"\n  cells expressing either A3A or A3B: {len(expr_either)} / {len(m)}")

print("\n  A3A_fraction vs CNV   (EXPECTED ~ -0.44 tumor / ~ +0.07 normal):")
for tis in ["tumor", "normal"]:
    sub = expr_either[expr_either["tissue type"] == tis]
    if len(sub) > 10:
        rho, p = spearmanr(sub["A3A_fraction"], sub["cnv_score"])
        r, _ = pearsonr(sub["A3A_fraction"], sub["cnv_score"])
        print(f"     {tis:7s} n={len(sub):6d}  spearman={rho:+.3f} (p={p:.1e})  pearson={r:+.3f}")

print("\n  A3A_fraction vs stemness  (EXPECTED ~ -0.5 in BOTH; NOT tumor-specific):")
for tis in ["tumor", "normal"]:
    sub = expr_either[expr_either["tissue type"] == tis]
    if len(sub) > 10:
        rho, _ = spearmanr(sub["A3A_fraction"], sub["CytoTRACE2_Score"])
        print(f"     {tis:7s} spearman={rho:+.3f}")

print("\n  SBS2 does NOT separate the enzymes (sanity; expect both small +):")
for g in ["APOBEC3A", "APOBEC3B"]:
    sub = m[(m["tissue type"] == "tumor")].dropna(subset=[g, "sig_SBS2"])
    rho, _ = spearmanr(sub[g], sub["sig_SBS2"])
    print(f"     tumor {g} vs sig_SBS2 spearman={rho:+.3f}")

# A3A- vs A3B-dominant cells: CNV/stemness contrast (the PI's exact phrasing)
hr("C2. A3A-dominant vs A3B-dominant cells (tumor, expressing either)")
tum = expr_either[expr_either["tissue type"] == "tumor"].copy()
a3a_dom = tum[tum["A3A_fraction"] > 0.5 + DOM_MARGIN]
a3b_dom = tum[tum["A3A_fraction"] < 0.5 - DOM_MARGIN]
print(f"  A3A-dominant n={len(a3a_dom)} | A3B-dominant n={len(a3b_dom)}")
print(f"  {'metric':18s} {'A3A-dom med':>12s} {'A3B-dom med':>12s}  expectation")
for col, exp in [("cnv_score", "A3B-dom HIGHER"),
                 ("CytoTRACE2_Score", "A3B-dom HIGHER"),
                 ("sig_SBS2", "A3A-dom >= A3B-dom")]:
    ma, mb = a3a_dom[col].median(), a3b_dom[col].median()
    print(f"  {col:18s} {ma:12.4f} {mb:12.4f}  {exp}")
if len(a3a_dom) > 5 and len(a3b_dom) > 5:
    u, p = mannwhitneyu(a3a_dom["cnv_score"].dropna(),
                        a3b_dom["cnv_score"].dropna(), alternative="two-sided")
    print(f"  Mann-Whitney CNV (A3A-dom vs A3B-dom): p={p:.2e}")

# CNV-quartile view of A3B dominance (clean 'distribution' table for the PI)
hr("C3. A3B dominance across CNV quartiles (tumor) -- monotonic increase expected")
tum["cnv_q"] = pd.qcut(tum["cnv_score"], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"],
                       duplicates="drop")
g = tum.groupby("cnv_q", observed=True).agg(
    n=("A3B_fraction", "size"),
    median_A3A=("APOBEC3A", "median"),
    median_A3B=("APOBEC3B", "median"),
    median_A3B_fraction=("A3B_fraction", "median"),
    pct_A3B_dominant=("A3A_fraction", lambda x: 100 * (x < 0.5).mean()),
)
print(g.round(4).to_string())

print("\n" + "=" * 78)
print("If tumor A3A_fraction-vs-CNV lands near -0.44 and A3B dominance rises across")
print("CNV quartiles, the table reconciles with Fig 3. Confirm the HPV block shows")
print("0 mismatches after regenerating with the >=8 fix before sending to the PI.")
print("=" * 78)
