#!/usr/bin/env python3
"""
Diagnostic_Section4_PanelD.py
================================
Extract the Panel D validation numbers from existing chain validation files.
Summarizes: activating internal DIFF, inhibiting internal DIFF,
cross-chain DIFF, A3-to-activating, A3-to-inhibiting.
"""

import os
import numpy as np
import pandas as pd

PROJECT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG_4 = os.path.join(PROJECT, "data", "FIG_4")
CHAIN_DIR = os.path.join(FIG_4, "DIAGNOSTIC_CHAIN_VALIDATION")

sep = "=" * 70

print(f"\n{sep}")
print("PANEL D VALIDATION NUMBERS")
print(sep)

# ============================================================
# 1. LOAD CHAIN GENE LIST
# ============================================================
chain_list = pd.read_csv(os.path.join(CHAIN_DIR, "chain_gene_list.tsv"), sep='\t')
activating = chain_list[chain_list['chain_type'] == 'activating']['gene'].tolist()
inhibiting = chain_list[chain_list['chain_type'] == 'inhibiting']['gene'].tolist()
a3_genes = chain_list[chain_list['chain_type'] == 'A3']['gene'].tolist()

print(f"\n  Activating chain: {len(activating)} genes")
for g in activating:
    harris = chain_list[chain_list['gene'] == g]['is_harris'].values[0]
    source = chain_list[chain_list['gene'] == g]['source'].values[0]
    h_tag = " [Harris]" if harris else ""
    print(f"    {g} ({source}){h_tag}")

print(f"\n  Inhibiting chain: {len(inhibiting)} genes")
harris_inh = chain_list[(chain_list['chain_type'] == 'inhibiting') & (chain_list['is_harris'] == True)]
print(f"    Harris interactors in inhibiting: {len(harris_inh)}")
for _, row in harris_inh.iterrows():
    print(f"      {row['gene']} ({row['source']})")

# Count by source
inh_by_source = chain_list[chain_list['chain_type'] == 'inhibiting']['source'].value_counts()
print(f"    By source: {inh_by_source.to_dict()}")

print(f"\n  A3 genes: {a3_genes}")

# ============================================================
# 2. LOAD DIFF MATRIX (SBS2 vs CNV)
# ============================================================
print(f"\n{sep}")
print("DIFF MATRIX (SBS2-HIGH vs CNV-HIGH correlations)")
print(sep)

diff_path = os.path.join(CHAIN_DIR, "DIFF_chain_genes.tsv")
diff = pd.read_csv(diff_path, sep='\t', index_col=0)
print(f"  DIFF matrix shape: {diff.shape}")

# Check which genes are present
all_chain = activating + inhibiting + a3_genes
present = [g for g in all_chain if g in diff.index]
missing = [g for g in all_chain if g not in diff.index]
print(f"  Present in DIFF: {len(present)} / {len(all_chain)}")
if missing:
    print(f"  Missing: {missing}")

act_present = [g for g in activating if g in diff.index]
inh_present = [g for g in inhibiting if g in diff.index]
a3_present = [g for g in a3_genes if g in diff.index]

# ============================================================
# 3. COMPUTE PANEL D SUMMARY STATS
# ============================================================
print(f"\n{sep}")
print("PANEL D COMPARISONS")
print(sep)

def pairwise_diff(genes1, genes2, label):
    """Compute pairwise DIFF values between two gene lists."""
    vals = []
    for g1 in genes1:
        for g2 in genes2:
            if g1 != g2 and g1 in diff.index and g2 in diff.columns:
                vals.append(diff.loc[g1, g2])
    vals = np.array(vals)
    if len(vals) == 0:
        print(f"\n  {label}: NO PAIRS")
        return
    print(f"\n  {label}:")
    print(f"    Pairs: {len(vals)}")
    print(f"    Mean DIFF:   {vals.mean():.4f}")
    print(f"    Median DIFF: {np.median(vals):.4f}")
    print(f"    Std:         {vals.std():.4f}")
    print(f"    Range:       [{vals.min():.4f}, {vals.max():.4f}]")
    print(f"    Positive (stronger in SBS2): {(vals > 0).sum()}/{len(vals)} ({100*(vals>0).sum()/len(vals):.1f}%)")
    print(f"    Negative (stronger in CNV):  {(vals < 0).sum()}/{len(vals)} ({100*(vals<0).sum()/len(vals):.1f}%)")
    return vals

# Comparison 1: Activating internal
print("\n  PREDICTION: Activating genes should co-express MORE tightly in SBS2-HIGH (positive DIFF)")
act_vals = pairwise_diff(act_present, act_present, "Activating-Activating (internal)")

# Comparison 2: Inhibiting internal
print("\n  PREDICTION: Inhibiting genes should co-express MORE tightly in CNV-HIGH (negative DIFF)")
inh_vals = pairwise_diff(inh_present, inh_present, "Inhibiting-Inhibiting (internal)")

# Comparison 3: Cross-chain
print("\n  PREDICTION: Cross-chain should be near zero or mixed (independent programs)")
cross_vals = pairwise_diff(act_present, inh_present, "Activating-Inhibiting (cross-chain)")

# Comparison 4: A3 to activating
print("\n  PREDICTION: A3 to activating - tests recoupling in SBS2 context")
for a3g in a3_present:
    alias = a3g.replace("APOBEC3", "A3")
    a3_act_vals = []
    for g in act_present:
        if g != a3g:
            a3_act_vals.append(diff.loc[a3g, g])
    a3_act_vals = np.array(a3_act_vals)
    if len(a3_act_vals) > 0:
        print(f"\n    {alias} -> Activating chain:")
        print(f"      Mean DIFF: {a3_act_vals.mean():.4f}")
        print(f"      Positive: {(a3_act_vals > 0).sum()}/{len(a3_act_vals)}")
        print(f"      Negative: {(a3_act_vals < 0).sum()}/{len(a3_act_vals)}")
        for g, v in zip(act_present, a3_act_vals):
            print(f"        {alias} <-> {g}: {v:.4f}")

# Comparison 5: A3 to inhibiting
print("\n  PREDICTION: A3 to inhibiting - tests decoupling")
for a3g in a3_present:
    alias = a3g.replace("APOBEC3", "A3")
    a3_inh_vals = []
    for g in inh_present[:20]:  # Sample first 20 for readability
        if g != a3g:
            a3_inh_vals.append(diff.loc[a3g, g])
    a3_inh_vals = np.array(a3_inh_vals)
    if len(a3_inh_vals) > 0:
        # Get all values for summary
        all_inh = [diff.loc[a3g, g] for g in inh_present if g != a3g]
        all_inh = np.array(all_inh)
        print(f"\n    {alias} -> Inhibiting chain (all {len(all_inh)} genes):")
        print(f"      Mean DIFF: {all_inh.mean():.4f}")
        print(f"      Positive: {(all_inh > 0).sum()}/{len(all_inh)}")
        print(f"      Negative: {(all_inh < 0).sum()}/{len(all_inh)}")

# Comparison 6: A3 to A3
print(f"\n  A3-A3 correlation:")
for i, g1 in enumerate(a3_present):
    for g2 in a3_present[i+1:]:
        v = diff.loc[g1, g2]
        a1 = g1.replace("APOBEC3", "A3")
        a2 = g2.replace("APOBEC3", "A3")
        print(f"    {a1} <-> {a2}: {v:.4f}")

# ============================================================
# 4. LOAD CHAIN VALIDATION DETAIL (from earlier diagnostic)
# ============================================================
print(f"\n{sep}")
print("CHAIN VALIDATION DETAIL (per-gene summary)")
print(sep)

detail_path = os.path.join(CHAIN_DIR, "chain_validation_detail.tsv")
if os.path.exists(detail_path):
    detail = pd.read_csv(detail_path, sep='\t')
    print(f"  Columns: {list(detail.columns)}")
    
    # Show activating genes
    print(f"\n  Activating chain genes:")
    act_detail = detail[detail['chain_type'] == 'activating']
    for _, row in act_detail.iterrows():
        h = " [Harris]" if row.get('is_harris', False) else ""
        print(f"    {row['gene']:15s} mean_diff_to_act={row['mean_diff_to_activating']:.4f}  "
              f"mean_diff_to_inh={row['mean_diff_to_inhibiting']:.4f}  "
              f"mean_diff_to_A3={row['mean_diff_to_A3']:.4f}{h}")
    
    # Show A3 genes
    print(f"\n  A3 genes:")
    a3_detail = detail[detail['chain_type'] == 'A3']
    for _, row in a3_detail.iterrows():
        alias = row['gene'].replace("APOBEC3", "A3")
        print(f"    {alias:15s} mean_diff_to_act={row['mean_diff_to_activating']:.4f}  "
              f"mean_diff_to_inh={row['mean_diff_to_inhibiting']:.4f}  "
              f"mean_diff_to_A3={row['mean_diff_to_A3']:.4f}")

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
