#!/usr/bin/env python
"""
Diagnostic_A3_negative_fraction.py

Computes the fraction of basal epithelial cells with NO detectable A3A or A3B
expression, both overall and within the SBS2-positive subset. The SBS2-positive
number is the one for the Results 4.1 sentence ("A subset of SBS2-positive basal
epithelial cells (X%) carried no detectable A3A or A3B expression").

"No detectable expression" = log-normalized expression == 0 (i.e. no captured
reads). A3A/A3B columns in the master table are log-normalized, so 0 means the
transcript was not captured in that cell.

Two sanity cross-checks against numbers already in the manuscript:
  A3A-positive %  should be ~25.7%   (Results 4.1)
  A3B-positive %  should be ~33.2%   (Results 4.1)

Run in NETWORK env on titan. Read-only.
"""

import pandas as pd

MASTER = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/"
          "01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")

df = pd.read_csv(MASTER, sep="\t").rename(columns={"Unnamed: 0": "cell_barcode"})

a3a = pd.to_numeric(df["APOBEC3A"], errors="coerce").fillna(0)
a3b = pd.to_numeric(df["APOBEC3B"], errors="coerce").fillna(0)
sbs2 = pd.to_numeric(df["sig_SBS2"], errors="coerce")   # NaN = never refit

double_neg = (a3a <= 0) & (a3b <= 0)
n_basal = len(df)

print("=" * 70)
print("A3-negative fraction in basal epithelial cells")
print("=" * 70)
print(f"basal cells total: {n_basal:,}")

# cross-checks against Results 4.1
print("\nsanity cross-checks (expect ~25.7% A3A+, ~33.2% A3B+):")
print(f"   A3A-positive: {int((a3a > 0).sum()):,} ({100*(a3a > 0).mean():.1f}%)")
print(f"   A3B-positive: {int((a3b > 0).sum()):,} ({100*(a3b > 0).mean():.1f}%)")

# all basal cells
print("\nall basal cells with NO A3A and NO A3B:")
print(f"   {int(double_neg.sum()):,} / {n_basal:,} ({100*double_neg.mean():.1f}%)")

# SBS2-positive subset, broken out by A3A/A3B status
sbs2_pos = sbs2 > 0
n_sbs2 = int(sbs2_pos.sum())
a3a_s = a3a[sbs2_pos]
a3b_s = a3b[sbs2_pos]

def pct(x):
    return 100 * x / n_sbs2 if n_sbs2 else 0.0

n_a3a_zero = int((a3a_s <= 0).sum())
n_a3b_zero = int((a3b_s <= 0).sum())

both_pos = int(((a3a_s > 0) & (a3b_s > 0)).sum())
only_a3a = int(((a3a_s > 0) & (a3b_s <= 0)).sum())
only_a3b = int(((a3a_s <= 0) & (a3b_s > 0)).sum())
neither  = int(((a3a_s <= 0) & (a3b_s <= 0)).sum())

print("\nSBS2-positive basal cells (sig_SBS2 > 0):")
print(f"   total: {n_sbs2:,}   (expect 5,911 from Results 4.1)")

print("\n   marginal (each enzyme independently):")
print(f"     A3A == 0: {n_a3a_zero:,} ({pct(n_a3a_zero):.1f}%)")
print(f"     A3B == 0: {n_a3b_zero:,} ({pct(n_a3b_zero):.1f}%)")

print("\n   four-way partition (mutually exclusive, sums to total):")
print(f"     A3A>0 & A3B>0  (both present): {both_pos:,} ({pct(both_pos):.1f}%)")
print(f"     A3A>0 & A3B==0 (A3A only)    : {only_a3a:,} ({pct(only_a3a):.1f}%)")
print(f"     A3A==0 & A3B>0 (A3B only)    : {only_a3b:,} ({pct(only_a3b):.1f}%)")
print(f"     A3A==0 & A3B==0 (neither)    : {neither:,} ({pct(neither):.1f}%)")
_chk = both_pos + only_a3a + only_a3b + neither
print(f"     partition sum check: {_chk:,} == {n_sbs2:,}  "
      f"{'OK' if _chk == n_sbs2 else 'MISMATCH'}")

print(f"\n   >>> SENTENCE NUMBER (neither A3A nor A3B): {pct(neither):.1f}% <<<")
print("=" * 70)
