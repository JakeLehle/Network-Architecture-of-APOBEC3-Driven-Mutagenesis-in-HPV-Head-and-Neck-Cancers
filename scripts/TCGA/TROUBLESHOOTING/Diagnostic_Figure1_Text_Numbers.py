#!/usr/bin/env python3
"""
Diagnostic_Figure1_Text_Numbers.py
=====================================
Extract all key numbers needed for writing Section 4.1 results text
and Figure 1 caption from the matched HNSC dataset.

Input: data/FIG_1/HNSC_A3_SBS2_matched_v2.tsv
Output: Console report + data/FIG_1/TROUBLESHOOTING/figure1_text_numbers.txt

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
DATA_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "HNSC_A3_SBS2_matched_v2.tsv")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)

report_lines = []
def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)

# =============================================================================
banner("LOAD DATA")
df = pd.read_csv(DATA_PATH, sep='\t')
log(f"  Loaded: {DATA_PATH}")
log(f"  Shape: {df.shape}")
log(f"  Columns: {list(df.columns)}")

# =============================================================================
banner("1. BASIC COUNTS")
n_total = len(df)
n_sbs2_pos = (df['SBS2'] > 0).sum()
n_sbs2_zero = (df['SBS2'] == 0).sum()
median_sbs2 = df['SBS2'].median()
mean_sbs2 = df['SBS2'].mean()

log(f"  Total HNSC tumors: {n_total}")
log(f"  SBS2 > 0: {n_sbs2_pos} ({100*n_sbs2_pos/n_total:.1f}%)")
log(f"  SBS2 = 0: {n_sbs2_zero} ({100*n_sbs2_zero/n_total:.1f}%)")
log(f"  SBS2 range: {df['SBS2'].min():.0f} - {df['SBS2'].max():.0f}")
log(f"  SBS2 median: {median_sbs2:.0f}")
log(f"  SBS2 mean: {mean_sbs2:.1f}")
log(f"  SBS2 Q25: {df['SBS2'].quantile(0.25):.0f}")
log(f"  SBS2 Q75: {df['SBS2'].quantile(0.75):.0f}")

# =============================================================================
banner("2. A3 EXPRESSION SUMMARY")
for gene in ['A3A', 'A3B', 'A3A_plus_A3B']:
    vals = df[gene]
    log(f"  {gene}: min={vals.min():.2f}  max={vals.max():.2f}  "
        f"mean={vals.mean():.2f}  median={vals.median():.2f}")

# Median A3A and A3B for stratification reference
median_a3a = df['A3A'].median()
median_a3b = df['A3B'].median()
log(f"\n  Median A3A: {median_a3a:.2f}")
log(f"  Median A3B: {median_a3b:.2f}")

# =============================================================================
banner("3. REGION COUNTS (from Panel 1a)")
if 'region' in df.columns:
    for region in ['teal', 'coral', 'cream']:
        n = (df['region'] == region).sum()
        log(f"  {region}: n={n} ({100*n/n_total:.1f}%)")
else:
    log("  WARNING: 'region' column not found. Run Step05 first.")

# =============================================================================
banner("4. LOWEST-A3 TUMORS IN HIGH-SBS2 GROUP")
high_sbs2 = df[df['SBS2'] > median_sbs2]
lowest_4 = high_sbs2.nsmallest(4, 'A3A_plus_A3B')
log(f"  Tumors with SBS2 > {median_sbs2:.0f} (median): {len(high_sbs2)}")
log(f"  4 lowest A3A+A3B in this group:")
for _, r in lowest_4.iterrows():
    eid = r['Entity_ID'][:12] if 'Entity_ID' in r.index else 'N/A'
    log(f"    {eid}: A3A={r['A3A']:.2f}, A3B={r['A3B']:.2f}, "
        f"A3A+A3B={r['A3A_plus_A3B']:.2f}, SBS2={r['SBS2']:.0f}")

# Also: how many have A3A+A3B < 1.0 in the high group?
n_low_a3_high_sbs2 = (high_sbs2['A3A_plus_A3B'] < 1.0).sum()
log(f"\n  High-SBS2 tumors with A3A+A3B < 1.0: {n_low_a3_high_sbs2}")
log(f"  High-SBS2 tumors with A3A+A3B < 2.0: {(high_sbs2['A3A_plus_A3B'] < 2.0).sum()}")
log(f"  High-SBS2 tumors with A3A+A3B < 5.0: {(high_sbs2['A3A_plus_A3B'] < 5.0).sum()}")

# =============================================================================
banner("5. A3B BASELINE / A3A AMPLIFICATION (Panel 1b context)")

# Stratify by median A3A and A3B
df['A3A_status'] = np.where(df['A3A'] >= median_a3a, 'HIGH', 'LOW')
df['A3B_status'] = np.where(df['A3B'] >= median_a3b, 'HIGH', 'LOW')
df['quadrant'] = df['A3B_status'] + '_A3B_' + df['A3A_status'] + '_A3A'

log(f"  Median-split SBS2 by A3A/A3B quadrant:")
log(f"  {'Quadrant':30s} {'n':>5s} {'median_SBS2':>12s} {'mean_SBS2':>10s} {'SBS2>0':>8s}")
log(f"  {'-'*30} {'-'*5} {'-'*12} {'-'*10} {'-'*8}")

quadrant_order = ['LOW_A3B_LOW_A3A', 'HIGH_A3B_LOW_A3A', 'LOW_A3B_HIGH_A3A', 'HIGH_A3B_HIGH_A3A']
for q in quadrant_order:
    sub = df[df['quadrant'] == q]
    n = len(sub)
    med = sub['SBS2'].median()
    mn = sub['SBS2'].mean()
    pos = (sub['SBS2'] > 0).sum()
    log(f"  {q:30s} {n:5d} {med:12.1f} {mn:10.1f} {pos:8d}")

# Wilcoxon tests between key comparisons
log(f"\n  Key pairwise comparisons (Wilcoxon rank-sum):")
comparisons = [
    ('LOW_A3B_LOW_A3A', 'HIGH_A3B_LOW_A3A', 'Effect of A3B alone'),
    ('LOW_A3B_LOW_A3A', 'LOW_A3B_HIGH_A3A', 'Effect of A3A alone'),
    ('HIGH_A3B_LOW_A3A', 'HIGH_A3B_HIGH_A3A', 'Effect of adding A3A to A3B'),
    ('LOW_A3B_LOW_A3A', 'HIGH_A3B_HIGH_A3A', 'Both vs neither'),
]
for q1, q2, label in comparisons:
    g1 = df[df['quadrant'] == q1]['SBS2']
    g2 = df[df['quadrant'] == q2]['SBS2']
    if len(g1) > 5 and len(g2) > 5:
        stat, p = mannwhitneyu(g1, g2, alternative='two-sided')
        log(f"    {label}:")
        log(f"      {q1} (median={g1.median():.1f}) vs {q2} (median={g2.median():.1f})")
        log(f"      p = {p:.2e}")

# =============================================================================
banner("6. INDIVIDUAL A3 CORRELATIONS WITH SBS2")
for gene in ['A3A', 'A3B', 'APOBEC3C', 'APOBEC3H']:
    col = gene if gene in df.columns else None
    if col is None:
        continue
    rho, p = spearmanr(df[col], df['SBS2'])
    log(f"  {gene} vs SBS2: Spearman rho = {rho:.4f}, p = {p:.2e}")

# A3A+A3B combined
rho, p = spearmanr(df['A3A_plus_A3B'], df['SBS2'])
log(f"  A3A+A3B vs SBS2: Spearman rho = {rho:.4f}, p = {p:.2e}")

# =============================================================================
banner("7. DATA PROVENANCE NOTE")
log(f"""
  SBS2 weights source: SigProfilerAssignment v3.4 (COSMIC v3.4, exome mode)
  Built from: MuTect2-called PASS-filtered somatic SNVs from GDC
  Expression source: TCGA STAR-Counts FPKM-UQ
  Matching: RNA-seq ↔ WES via direct crosswalk + Case_ID (sample-type constrained)
  
  Old file had: 522 HNSC tumors, 65 SBS signatures
  New file has: {n_total} HNSC tumors, 86 SBS signatures (COSMIC v3.4)
  
  NOTE FOR TEXT: Update all numbers from the old 522/6-panel analysis.
  The simplified 2-panel figure uses A3A+A3B only (not A3A+A3B+A3C+A3H),
  so the x-axis and region assignments differ from the original.
""")

# =============================================================================
banner("8. TEXT-READY NUMBERS SUMMARY")
log(f"""
  FOR SECTION 4.1:
    n tumors = {n_total}
    SBS2 range: 0 to {df['SBS2'].max():.0f}
    Tumors with SBS2 > 0: {n_sbs2_pos}/{n_total} ({100*n_sbs2_pos/n_total:.1f}%)
    Median SBS2: {median_sbs2:.0f}
    
    Region counts (Panel 1a):
      Coral (A3 + HIGH SBS2): {(df['region']=='coral').sum() if 'region' in df.columns else '?'}
      Cream (A3 + LOW SBS2): {(df['region']=='cream').sum() if 'region' in df.columns else '?'}
      Teal (no A3 + HIGH SBS2): {(df['region']=='teal').sum() if 'region' in df.columns else '?'}
    
    Low-A3 outliers in high-SBS2 group: {n_low_a3_high_sbs2} (A3A+A3B < 1.0)
    
  FOR PANEL 1b NARRATIVE:
    A3B-high/A3A-low median SBS2: {df[df['quadrant']=='HIGH_A3B_LOW_A3A']['SBS2'].median():.1f}
    A3B-low/A3A-high median SBS2: {df[df['quadrant']=='LOW_A3B_HIGH_A3A']['SBS2'].median():.1f}
    Both high median SBS2: {df[df['quadrant']=='HIGH_A3B_HIGH_A3A']['SBS2'].median():.1f}
    Neither median SBS2: {df[df['quadrant']=='LOW_A3B_LOW_A3A']['SBS2'].median():.1f}
    
  FOR METHODS:
    A3A vs SBS2 Spearman rho (for methods, NOT figure)
    A3B vs SBS2 Spearman rho (for methods, NOT figure)
    Matching: {len(df)} tumors matched from RNA-seq (FPKM-UQ) and WES (SigProfiler v3.4)
""")

# Save
report_path = os.path.join(OUTPUT_DIR, "figure1_text_numbers.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
banner("COMPLETE")
