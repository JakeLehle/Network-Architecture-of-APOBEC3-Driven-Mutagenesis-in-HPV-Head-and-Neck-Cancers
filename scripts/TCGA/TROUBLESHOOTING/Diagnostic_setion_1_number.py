#!/usr/bin/env python3
"""
Diagnostic_Section1_Numbers.py
================================
Pull every number needed for the ~500-word Section 1 rewrite.
Run from anywhere; all paths are absolute.

Outputs a plain-text report to stdout.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, kruskal

PROJECT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG1 = os.path.join(PROJECT, "data", "FIG_1")

sep = "=" * 70

# ============================================================
# 1. LOAD MATCHED DATASET
# ============================================================
print(f"\n{sep}")
print("SECTION 1 TEXT NUMBERS — FIGURE 1")
print(sep)

fpath = os.path.join(FIG1, "HNSC_A3_SBS2_matched_v3.tsv")
if not os.path.exists(fpath):
    # Try alternate name
    fpath = os.path.join(FIG1, "HNSC_A3_SBS2_matched.tsv")
if not os.path.exists(fpath):
    print(f"ERROR: Cannot find matched data in {FIG1}")
    print("Available files:")
    for f in sorted(os.listdir(FIG1)):
        print(f"  {f}")
    exit(1)

df = pd.read_csv(fpath, sep='\t')
print(f"\nLoaded: {fpath}")
print(f"Shape: {df.shape}")
print(f"Columns: {list(df.columns)}")

n_total = len(df)
print(f"\n{'—'*50}")
print(f"TOTAL TUMORS: {n_total}")
print(f"{'—'*50}")

# ============================================================
# 2. SBS2 BASICS
# ============================================================
print(f"\n--- SBS2 Distribution ---")
print(f"  SBS2 > 0:  {(df['SBS2'] > 0).sum()} / {n_total} ({100*(df['SBS2']>0).sum()/n_total:.1f}%)")
print(f"  SBS2 = 0:  {(df['SBS2'] == 0).sum()} / {n_total}")
print(f"  Median:    {df['SBS2'].median():.1f}")
print(f"  Mean:      {df['SBS2'].mean():.1f}")
print(f"  Max:       {df['SBS2'].max():.0f}")
print(f"  Q25:       {df['SBS2'].quantile(0.25):.1f}")
print(f"  Q75:       {df['SBS2'].quantile(0.75):.1f}")

# ============================================================
# 3. A3 EXPRESSION
# ============================================================
print(f"\n--- A3 Expression ---")
for col in ['A3A', 'A3B']:
    if col in df.columns:
        print(f"  {col}: median={df[col].median():.2f}, mean={df[col].mean():.2f}, max={df[col].max():.2f}")

if 'A3A_plus_A3B' in df.columns:
    a3sum_col = 'A3A_plus_A3B'
elif 'A3_sum' in df.columns:
    a3sum_col = 'A3_sum'
else:
    # Try to compute
    if 'A3A' in df.columns and 'A3B' in df.columns:
        df['A3_sum'] = df['A3A'] + df['A3B']
        a3sum_col = 'A3_sum'
    else:
        a3sum_col = None
        print("  WARNING: No A3 sum column found")

if a3sum_col:
    print(f"  {a3sum_col}: median={df[a3sum_col].median():.2f}, mean={df[a3sum_col].mean():.2f}")

# ============================================================
# 4. REGION COUNTS (Panel 1a)
# ============================================================
print(f"\n--- Region Counts (Panel 1a) ---")
if 'region' in df.columns:
    for region in df['region'].unique():
        n = (df['region'] == region).sum()
        print(f"  {region}: n={n} ({100*n/n_total:.1f}%)")
    
    # Tumors near the boundary: low A3 but some SBS2
    if a3sum_col:
        median_sbs2 = df['SBS2'].median()
        near_boundary = df[(df[a3sum_col] < 1.0) & (df['SBS2'] > median_sbs2)]
        print(f"\n  Near-boundary tumors (A3sum < 1.0 AND SBS2 > median {median_sbs2:.0f}): {len(near_boundary)}")
        if len(near_boundary) > 0:
            for _, row in near_boundary.iterrows():
                print(f"    A3sum={row[a3sum_col]:.2f}, SBS2={row['SBS2']:.0f}")
        
        # Tumors with zero A3 and any SBS2
        zero_a3_any_sbs2 = df[(df[a3sum_col] == 0) & (df['SBS2'] > 0)]
        print(f"  Tumors with A3sum=0 AND SBS2>0: {len(zero_a3_any_sbs2)}")
        
        # Tumors in teal region (if applicable)
        if 'teal' in df['region'].values:
            teal = df[df['region'] == 'teal']
            print(f"  Teal region A3sum range: {teal[a3sum_col].min():.2f} - {teal[a3sum_col].max():.2f}")
            print(f"  Teal region SBS2 range: {teal['SBS2'].min():.0f} - {teal['SBS2'].max():.0f}")
else:
    print("  No 'region' column — listing all available columns:")
    print(f"  {list(df.columns)}")

# ============================================================
# 5. SPEARMAN CORRELATIONS
# ============================================================
print(f"\n--- Spearman Correlations with SBS2 ---")
for col in ['A3A', 'A3B']:
    if col in df.columns:
        rho, p = spearmanr(df[col], df['SBS2'])
        print(f"  {col} vs SBS2: rho={rho:.4f}, p={p:.2e}")

if a3sum_col:
    rho, p = spearmanr(df[a3sum_col], df['SBS2'])
    print(f"  {a3sum_col} vs SBS2: rho={rho:.4f}, p={p:.2e}")

# ============================================================
# 6. QUADRANT ANALYSIS (Panel 1c)
# ============================================================
print(f"\n--- Quadrant Analysis (Panel 1b/1c) ---")
if 'quadrant' in df.columns:
    for q in df['quadrant'].unique():
        sub = df[df['quadrant'] == q]
        print(f"  {q}: n={len(sub)}, median SBS2={sub['SBS2'].median():.1f}")
    
    # Key comparisons
    if 'HIGH_A3A_LOW_A3B' in df['quadrant'].values and 'LOW_A3B_LOW_A3A' in df['quadrant'].values:
        g1 = df[df['quadrant'] == 'HIGH_A3A_LOW_A3B']['SBS2']
        g0 = df[df['quadrant'] == 'LOW_A3B_LOW_A3A']['SBS2']
        stat, p = mannwhitneyu(g1, g0, alternative='two-sided')
        print(f"\n  A3A-high/A3B-low vs both-low: U={stat:.0f}, p={p:.2e}")
    
    # Also try alternate naming
    for q_high in ['HIGH_A3B_HIGH_A3A', 'HIGH_A3A_HIGH_A3B']:
        for q_a3a in ['HIGH_A3A_LOW_A3B', 'LOW_A3B_HIGH_A3A']:
            if q_high in df['quadrant'].values and q_a3a in df['quadrant'].values:
                g_both = df[df['quadrant'] == q_high]['SBS2']
                g_a3a = df[df['quadrant'] == q_a3a]['SBS2']
                stat, p = mannwhitneyu(g_both, g_a3a, alternative='two-sided')
                print(f"  Both-high vs A3A-high: U={stat:.0f}, p={p:.2e}")
    
    # Kruskal-Wallis across all 4
    groups = [df[df['quadrant'] == q]['SBS2'].values for q in df['quadrant'].unique()]
    if len(groups) == 4:
        stat, p = kruskal(*groups)
        print(f"  Kruskal-Wallis (4 groups): H={stat:.1f}, p={p:.2e}")
elif 'A3A' in df.columns and 'A3B' in df.columns:
    # Compute quadrants from medians
    med_a3a = df['A3A'].median()
    med_a3b = df['A3B'].median()
    print(f"  (Computing from medians: A3A={med_a3a:.2f}, A3B={med_a3b:.2f})")
    
    conditions = [
        ('LOW_A3A_LOW_A3B', (df['A3A'] <= med_a3a) & (df['A3B'] <= med_a3b)),
        ('HIGH_A3A_LOW_A3B', (df['A3A'] > med_a3a) & (df['A3B'] <= med_a3b)),
        ('LOW_A3A_HIGH_A3B', (df['A3A'] <= med_a3a) & (df['A3B'] > med_a3b)),
        ('HIGH_A3A_HIGH_A3B', (df['A3A'] > med_a3a) & (df['A3B'] > med_a3b)),
    ]
    for name, mask in conditions:
        sub = df[mask]
        print(f"  {name}: n={len(sub)}, median SBS2={sub['SBS2'].median():.1f}")
    
    # Key p-values
    low_low = df[(df['A3A'] <= med_a3a) & (df['A3B'] <= med_a3b)]['SBS2']
    high_a3a = df[(df['A3A'] > med_a3a) & (df['A3B'] <= med_a3b)]['SBS2']
    high_a3b = df[(df['A3A'] <= med_a3a) & (df['A3B'] > med_a3b)]['SBS2']
    both_high = df[(df['A3A'] > med_a3a) & (df['A3B'] > med_a3b)]['SBS2']
    
    stat, p = mannwhitneyu(high_a3a, low_low, alternative='two-sided')
    print(f"\n  A3A-high vs both-low: p={p:.2e}")
    stat, p = mannwhitneyu(high_a3b, low_low, alternative='two-sided')
    print(f"  A3B-high vs both-low: p={p:.2e}")
    stat, p = mannwhitneyu(both_high, high_a3a, alternative='two-sided')
    print(f"  Both-high vs A3A-high: p={p:.2e}")
    stat, p = mannwhitneyu(both_high, high_a3b, alternative='two-sided')
    print(f"  Both-high vs A3B-high: p={p:.2e}")

# ============================================================
# 7. A3B SATURATION
# ============================================================
print(f"\n--- A3B Saturation Analysis ---")
if 'A3B' in df.columns:
    # Find percentile where A3B-SBS2 correlation breaks down
    percentiles = list(range(50, 96, 5))
    print(f"  A3B threshold sweep (Spearman rho for tumors ABOVE each percentile):")
    for pct in percentiles:
        thresh = np.percentile(df['A3B'], pct)
        above = df[df['A3B'] >= thresh]
        if len(above) >= 10:
            rho, p = spearmanr(above['A3B'], above['SBS2'])
            print(f"    >{pct}th pct (A3B>={thresh:.2f}, n={len(above)}): rho={rho:.3f}, p={p:.2e}")
    
    print(f"\n  A3B threshold sweep (Spearman rho for tumors BELOW each percentile):")
    for pct in percentiles:
        thresh = np.percentile(df['A3B'], pct)
        below = df[df['A3B'] <= thresh]
        if len(below) >= 10:
            rho, p = spearmanr(below['A3B'], below['SBS2'])
            print(f"    <{pct}th pct (A3B<={thresh:.2f}, n={len(below)}): rho={rho:.3f}, p={p:.2e}")

# ============================================================
# 8. GERMLINE / SOMATIC SUMMARY
# ============================================================
print(f"\n--- Germline & Somatic Summary ---")
print("  (These numbers come from separate analyses; listing known values)")
print("  Germline SNPs tested: 6,938")
print("  Nominal p<0.05: 93")
print("  BH-significant (adj p<0.25): 0")
print("  Top KEGG: B cell receptor signaling, Influenza A")
print()
print("  Somatic: HIGH vs LOW matched groups (n=53 each)")
print("  Mutation burden ratio: ~1.5x (HIGH > LOW)")

# Check if group assignment files exist
group_file = os.path.join(FIG1, "network_group_assignments_v3.tsv")
if not os.path.exists(group_file):
    group_file = os.path.join(FIG1, "network_group_assignments.tsv")
if os.path.exists(group_file):
    groups = pd.read_csv(group_file, sep='\t')
    print(f"\n  Group assignment file found: {group_file}")
    print(f"  Columns: {list(groups.columns)}")
    if 'group' in groups.columns:
        for g in groups['group'].unique():
            print(f"    {g}: n={len(groups[groups['group']==g])}")
    elif 'Group' in groups.columns:
        for g in groups['Group'].unique():
            print(f"    {g}: n={len(groups[groups['Group']==g])}")

# ============================================================
# 9. FILE INVENTORY
# ============================================================
print(f"\n--- FIG_1 Directory Contents ---")
for f in sorted(os.listdir(FIG1)):
    fsize = os.path.getsize(os.path.join(FIG1, f))
    print(f"  {f}  ({fsize/1024:.1f} KB)")

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
