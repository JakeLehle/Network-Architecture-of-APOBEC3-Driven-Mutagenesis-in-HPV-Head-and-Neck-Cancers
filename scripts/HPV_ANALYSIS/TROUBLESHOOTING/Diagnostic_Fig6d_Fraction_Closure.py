#!/usr/bin/env python3
"""
Diagnostic_Fig6d_Fraction_Closure.py
====================================
2026_NMF_PAPER | Figure 6d audit -- do the HPV16 category fractions close?

WHY THIS EXISTS
---------------
Section 5.3 reports five HPV16 read categories for the two tumour populations:

    SBS2-HIGH   URR 63.5   maintenance 25.9   assembly 10.3   oncogene 0.53
    CNA-HIGH    URR 63.5   maintenance 13.2   assembly 17.8   oncogene 0.79

On a single denominator those must sum to at most 100% with amplification
making up the remainder. SBS2-HIGH sums to 100.23% BEFORE amplification is
counted, yet the same paragraph reports an E5 amplification contrast. One of
the numbers is on a different denominator, a different estimator, or a
different cell set.

Reading the pipeline, the likely cause is an ESTIMATOR MIX:

  - Generate_Figure6_Lifecycle_Panels.py Panel F reports the MEAN OF PER-CELL
    FRACTIONS (gene / total_hpv16_genome_reads, no pseudocount) on the gated
    set. The phase and gene percentages in the prose come from there.
  - Diagnostic_Figure6_HostMarkers_and_IntegrationProxy.py Section 2 reports
    URR POOLED (sum URR / sum total) and notes the prose should cite pooled.

Mean-of-per-cell-fractions and pooled are different estimators of the same
quantity and are not additive with one another. There is also an `intergenic`
read class in the gene-count table that the manuscript never mentions, so
total = URR + ORF_sum + intergenic and the five prose categories were never
going to close even on one estimator.

WHAT THIS SCRIPT DOES
---------------------
  1. Reproduces the Panel F gate exactly (raw_HPV16 >= 8 AND total > 0).
  2. Computes EVERY read class both ways, pooled and mean-of-per-cell, for all
     three populations: URR, intergenic, and the four functional gene groups.
  3. Prints the closure sum for each estimator. A correct estimator closes to
     100.00% by construction; anything else is a bug worth chasing.
  4. Diffs each manuscript number against both estimators and reports which one
     it came from.
  5. Prints the amplification and intergenic values, which the manuscript does
     not currently report and which are the missing mass.
  6. Emits a copy-pasteable sentence for each of the two options.

READ-ONLY. Writes only to the diagnostics directory. Touches no figure input.

Run:
    conda run -n NETWORK python Diagnostic_Fig6d_Fraction_Closure.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

# ============================== CONFIG ==============================
ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

THREE_GROUP_PATH = os.path.join(ROOT, "data/FIG_4/01_group_selection/"
                                      "three_group_assignments.tsv")
HPV_GENE_PATH    = os.path.join(ROOT, "data/FIG_6/03_hpv16_genome/"
                                      "per_cell_hpv16_gene_counts.tsv")
MASTER_HPV_PATH  = os.path.join(ROOT, "data/FIG_6/01_raw_hpv16_counts/"
                                      "basal_cell_master_table_with_raw_HPV16.tsv")
OUT_DIR = os.path.join(ROOT, "data/FIG_6/00_diagnostics")

TOTAL_COL       = "total_hpv16_genome_reads"
RAW_COL         = "raw_HPV16"
HPV16_THRESHOLD = 8

POP_ORDER  = ["SBS2_HIGH", "CNV_HIGH", "NORMAL"]
POP_LABELS = {"SBS2_HIGH": "SBS2-HIGH", "CNV_HIGH": "CNA-HIGH", "NORMAL": "NORMAL"}

# Functional viral gene groups, matching Methods 7.4 and Fig. 6d.
HPV16_GROUPS = {
    "Maintenance":   ["E1", "E2"],
    "Amplification": ["E4", "E5"],
    "Oncogene":      ["E6", "E7"],
    "Packaging":     ["L1", "L2"],
}
ORF_GENES = ["E1", "E2", "E4", "E5", "E6", "E7", "L1", "L2"]

# Section 5.3 as currently written. Edit here if the text changes.
# value is None where the manuscript reports nothing.
MANUSCRIPT = {
    ("SBS2_HIGH", "URR"):           63.5,
    ("CNV_HIGH",  "URR"):           63.5,
    ("NORMAL",    "URR"):           64.9,
    ("SBS2_HIGH", "Maintenance"):   25.9,
    ("CNV_HIGH",  "Maintenance"):   13.2,
    ("SBS2_HIGH", "Packaging"):     10.3,
    ("CNV_HIGH",  "Packaging"):     17.8,
    ("SBS2_HIGH", "Oncogene"):      0.53,
    ("CNV_HIGH",  "Oncogene"):      0.79,
    ("SBS2_HIGH", "Amplification"): None,
    ("CNV_HIGH",  "Amplification"): None,
    ("SBS2_HIGH", "intergenic"):    None,
    ("CNV_HIGH",  "intergenic"):    None,
}

# A manuscript value is called a match if it is within this many percentage
# points of the recomputed value. 0.05 covers 1-decimal rounding.
MATCH_TOL = 0.05
# ====================================================================

RUN_TS = datetime.now().strftime("%Y%m%d_%H%M%S")
_loglines = []


def log(m=""):
    print(m, flush=True)
    _loglines.append(m)


def banner(t):
    log("")
    log("=" * 78)
    log(f"  {t}")
    log("=" * 78)


def die(m):
    log(f"  ERROR: {m}")
    sys.exit(1)


# =============================================================================
# STEP 0: LOAD AND REPRODUCE THE PANEL F GATE
# =============================================================================
banner("STEP 0: INPUTS AND PANEL F GATE")

for p in (THREE_GROUP_PATH, HPV_GENE_PATH, MASTER_HPV_PATH):
    if not os.path.exists(p):
        die(f"not found: {p}")
os.makedirs(OUT_DIR, exist_ok=True)

groups = pd.read_csv(THREE_GROUP_PATH, sep="\t", index_col=0)
hpv    = pd.read_csv(HPV_GENE_PATH,    sep="\t", index_col=0)
master = pd.read_csv(MASTER_HPV_PATH,  sep="\t", index_col=0)
log(f"  three_group_assignments : {len(groups):,} rows")
log(f"  per_cell_hpv16_gene_counts: {len(hpv):,} rows")
log(f"  basal_cell_master_table : {len(master):,} rows")

# group column
gcol = next((c for c in ("group", "population", "assignment")
             if c in groups.columns), None)
if gcol is None:
    die(f"no group column in {THREE_GROUP_PATH}. Columns: {list(groups.columns)}")
log(f"  group column: '{gcol}'  values: "
    f"{sorted(groups[gcol].astype(str).unique())}")

if TOTAL_COL not in hpv.columns:
    die(f"'{TOTAL_COL}' not in gene-count table. Columns: {list(hpv.columns)}")
if RAW_COL not in master.columns:
    die(f"'{RAW_COL}' not in master table. Columns: {list(master.columns)}")

# ---- column inventory: what read classes actually exist -------------------
present_orfs = [g for g in ORF_GENES if g in hpv.columns]
missing_orfs = [g for g in ORF_GENES if g not in hpv.columns]
has_urr = "URR" in hpv.columns
other_cols = [c for c in hpv.columns
              if c not in set(present_orfs) | {TOTAL_COL, "URR"}]
log("")
log(f"  ORFs present : {', '.join(present_orfs)}"
    + (f"   MISSING: {', '.join(missing_orfs)}" if missing_orfs else ""))
log(f"  URR column   : {'yes' if has_urr else 'NO'}")
log(f"  other columns: {other_cols if other_cols else 'none'}")
if not has_urr:
    die("no URR column; this diagnostic cannot run without it.")

# Read classes to audit. Anything numeric that is not an ORF, not URR and not
# the total is a candidate extra class (intergenic being the expected one).
extra_classes = [c for c in other_cols
                 if pd.api.types.is_numeric_dtype(hpv[c])]
log(f"  extra numeric read classes treated as categories: "
    f"{extra_classes if extra_classes else 'none'}")

# ---- gate ----------------------------------------------------------------
hpv = hpv.copy()
hpv["group"] = hpv.index.map(lambda x: groups[gcol].get(x, "other"))
hpv[RAW_COL] = hpv.index.map(master[RAW_COL])

gated = hpv[hpv["group"].isin(POP_ORDER)
            & (hpv[RAW_COL] >= HPV16_THRESHOLD)
            & (hpv[TOTAL_COL] > 0)].copy()

log("")
log(f"  Panel F gate: {RAW_COL} >= {HPV16_THRESHOLD} AND {TOTAL_COL} > 0")
log(f"  gated cells: {len(gated):,}")
for p in POP_ORDER:
    log(f"    {POP_LABELS[p]:<10s} {(gated['group'] == p).sum():>5d}")
log("  Manuscript states 197 / 446 / 8. Any disagreement here means the gate")
log("  or an upstream input has changed and nothing below is comparable.")

if len(gated) == 0:
    die("gate returned no cells")

gated["ORF_sum"] = gated[present_orfs].sum(axis=1)


# =============================================================================
# STEP 1: BOTH ESTIMATORS FOR EVERY READ CLASS
# =============================================================================
banner("STEP 1: EVERY READ CLASS, POOLED AND MEAN-OF-PER-CELL")

log("  pooled   = sum(class reads) / sum(total reads) across the group")
log("  per-cell = mean over cells of (class reads / that cell's total reads)")
log("  These estimate the same quantity and are NOT additive with each other.")

categories = ["URR"] + list(HPV16_GROUPS.keys()) + extra_classes

rows = []
for p in POP_ORDER:
    sub = gated[gated["group"] == p]
    if len(sub) == 0:
        continue
    tot_sum = sub[TOTAL_COL].sum()
    for cat in categories:
        if cat in HPV16_GROUPS:
            genes = [g for g in HPV16_GROUPS[cat] if g in sub.columns]
            reads = sub[genes].sum(axis=1) if genes else pd.Series(0.0, index=sub.index)
            members = ", ".join(genes) if genes else "none present"
        else:
            reads = sub[cat]
            members = cat
        pooled  = 100.0 * reads.sum() / tot_sum if tot_sum > 0 else np.nan
        percell = 100.0 * (reads / sub[TOTAL_COL]).mean()
        rows.append({"group": POP_LABELS[p], "n": len(sub), "category": cat,
                     "members": members,
                     "pooled_pct": pooled, "percell_mean_pct": percell})
res = pd.DataFrame(rows)

for p in POP_ORDER:
    lab = POP_LABELS[p]
    sub = res[res["group"] == lab]
    if sub.empty:
        continue
    n = int(sub["n"].iloc[0])
    log("")
    log(f"  {lab}  (n = {n})")
    log(f"    {'category':<16s} {'members':<12s} {'pooled':>10s} {'per-cell':>10s}")
    log(f"    {'-'*16} {'-'*12} {'-'*10} {'-'*10}")
    for _, r in sub.iterrows():
        log(f"    {r['category']:<16s} {r['members']:<12s} "
            f"{r['pooled_pct']:>9.2f}% {r['percell_mean_pct']:>9.2f}%")
    log(f"    {'-'*16} {'-'*12} {'-'*10} {'-'*10}")
    tp = sub["pooled_pct"].sum()
    tc = sub["percell_mean_pct"].sum()
    log(f"    {'CLOSURE SUM':<16s} {'':<12s} {tp:>9.2f}% {tc:>9.2f}%")
    for name, val in (("pooled", tp), ("per-cell", tc)):
        if abs(val - 100.0) > 0.5:
            log(f"    WARNING: {name} sum is {val:.2f}%, not 100%. Either a read "
                f"class is missing from the category list above, or the class "
                f"columns and {TOTAL_COL} use different denominators.")


# =============================================================================
# STEP 2: WHICH ESTIMATOR DID EACH MANUSCRIPT NUMBER COME FROM
# =============================================================================
banner("STEP 2: MANUSCRIPT NUMBER PROVENANCE")

log(f"  A value is called a match within {MATCH_TOL} percentage points.")
log("")
log(f"  {'group':<10s} {'category':<16s} {'text':>8s} {'pooled':>9s} "
    f"{'per-cell':>9s}  {'verdict'}")
log(f"  {'-'*10} {'-'*16} {'-'*8} {'-'*9} {'-'*9}  {'-'*28}")

prov_rows = []
for (p, cat), stated in MANUSCRIPT.items():
    hit = res[(res["group"] == POP_LABELS[p]) & (res["category"] == cat)]
    if hit.empty:
        log(f"  {POP_LABELS[p]:<10s} {cat:<16s} {'--':>8s} {'--':>9s} "
            f"{'--':>9s}  NOT COMPUTED")
        continue
    pooled  = float(hit["pooled_pct"].iloc[0])
    percell = float(hit["percell_mean_pct"].iloc[0])
    if stated is None:
        verdict = "not reported in text"
        stated_s = "--"
    else:
        dp, dc = abs(stated - pooled), abs(stated - percell)
        mp, mc = dp <= MATCH_TOL, dc <= MATCH_TOL
        if mp and mc:
            verdict = "matches BOTH (identical)"
        elif mp:
            verdict = "from POOLED"
        elif mc:
            verdict = "from PER-CELL"
        else:
            verdict = f"MATCHES NEITHER (d {min(dp, dc):.2f})"
        stated_s = f"{stated:.2f}"
    log(f"  {POP_LABELS[p]:<10s} {cat:<16s} {stated_s:>8s} {pooled:>8.2f}% "
        f"{percell:>8.2f}%  {verdict}")
    prov_rows.append({"group": POP_LABELS[p], "category": cat,
                      "manuscript": stated, "pooled_pct": pooled,
                      "percell_mean_pct": percell, "verdict": verdict})


# =============================================================================
# STEP 3: CLOSURE OF THE FIVE CATEGORIES THE MANUSCRIPT REPORTS
# =============================================================================
banner("STEP 3: DOES THE MANUSCRIPT'S OWN SET CLOSE")

log("  The text reports URR + maintenance + packaging + oncogene and asserts an")
log("  amplification contrast. Below, each estimator is summed over the four")
log("  reported categories, then amplification and any extra class are added.")

for p in ("SBS2_HIGH", "CNV_HIGH"):
    lab = POP_LABELS[p]
    sub = res[res["group"] == lab].set_index("category")
    if sub.empty:
        continue
    log("")
    log(f"  {lab}")
    reported = ["URR", "Maintenance", "Packaging", "Oncogene"]
    for est, col in (("pooled", "pooled_pct"), ("per-cell", "percell_mean_pct")):
        four = sum(float(sub.loc[c, col]) for c in reported if c in sub.index)
        amp  = float(sub.loc["Amplification", col]) if "Amplification" in sub.index else np.nan
        extra = sum(float(sub.loc[c, col]) for c in extra_classes if c in sub.index)
        log(f"    {est:<9s} four reported = {four:7.2f}%   "
            f"+ amplification {amp:6.2f}%  + extra {extra:6.2f}%  "
            f"= {four + (0 if np.isnan(amp) else amp) + extra:7.2f}%")

    # the mixed version the manuscript currently prints
    mixed = 0.0
    parts = []
    for c in reported:
        if c not in sub.index:
            continue
        col = "pooled_pct" if c == "URR" else "percell_mean_pct"
        mixed += float(sub.loc[c, col])
        parts.append(f"{c}={float(sub.loc[c, col]):.2f}")
    log(f"    AS PRINTED (URR pooled + rest per-cell): {mixed:.2f}%   "
        f"[{'; '.join(parts)}]")
    if mixed > 100.0:
        log("    This exceeds 100% before amplification, which is the error in "
            "the current text.")


# =============================================================================
# STEP 4: RECOMMENDED SENTENCES
# =============================================================================
banner("STEP 4: COPY-PASTEABLE REPLACEMENTS")

def val(p, cat, col):
    hit = res[(res["group"] == POP_LABELS[p]) & (res["category"] == cat)]
    return float(hit[col].iloc[0]) if not hit.empty else float("nan")

for est, col, name in (("pooled", "pooled_pct", "ALL-POOLED"),
                       ("per-cell", "percell_mean_pct", "ALL-PER-CELL")):
    log("")
    log(f"  OPTION: {name}")
    log(f"    In all three populations, approximately two-thirds of HPV16 reads")
    log(f"    mapped to the non-coding upstream regulatory region (URR) rather")
    log(f"    than to a viral open reading frame "
        f"({val('SBS2_HIGH','URR',col):.1f}% in SBS2-HIGH, "
        f"{val('CNV_HIGH','URR',col):.1f}% in CNA-HIGH, and "
        f"{val('NORMAL','URR',col):.1f}% in NORMAL cells).")
    log(f"    ... maintenance {val('SBS2_HIGH','Maintenance',col):.1f}% versus "
        f"{val('CNV_HIGH','Maintenance',col):.1f}% ...")
    log(f"    ... packaging {val('SBS2_HIGH','Packaging',col):.1f}% versus "
        f"{val('CNV_HIGH','Packaging',col):.1f}% ...")
    log(f"    ... oncogene {val('SBS2_HIGH','Oncogene',col):.2f}% and "
        f"{val('CNV_HIGH','Oncogene',col):.2f}% ...")
    log(f"    ... amplification {val('SBS2_HIGH','Amplification',col):.2f}% versus "
        f"{val('CNV_HIGH','Amplification',col):.2f}% (NOT currently in the text)")

log("")
log("  Whichever option is taken, Methods 7.4 and the Fig. 6d legend should")
log("  name the estimator explicitly, because pooled and mean-of-per-cell are")
log("  both defensible and the reader cannot tell them apart from the number.")
if extra_classes:
    log(f"  The text should also account for the {', '.join(extra_classes)} read")
    log("  class, or state that the reported categories do not sum to 100%")
    log("  because that class is excluded.")


# =============================================================================
# SAVE
# =============================================================================
banner("SAVE")

res.to_csv(os.path.join(OUT_DIR, f"fig6d_readclass_both_estimators_{RUN_TS}.tsv"),
           sep="\t", index=False)
log(f"  [SAVE] fig6d_readclass_both_estimators_{RUN_TS}.tsv")

if prov_rows:
    pd.DataFrame(prov_rows).to_csv(
        os.path.join(OUT_DIR, f"fig6d_manuscript_provenance_{RUN_TS}.tsv"),
        sep="\t", index=False)
    log(f"  [SAVE] fig6d_manuscript_provenance_{RUN_TS}.tsv")

logp = os.path.join(OUT_DIR, f"fig6d_fraction_closure_{RUN_TS}.txt")
with open(logp, "w") as fh:
    fh.write("\n".join(_loglines) + "\n")
log(f"  [SAVE] {os.path.basename(logp)}")

banner("DONE")
