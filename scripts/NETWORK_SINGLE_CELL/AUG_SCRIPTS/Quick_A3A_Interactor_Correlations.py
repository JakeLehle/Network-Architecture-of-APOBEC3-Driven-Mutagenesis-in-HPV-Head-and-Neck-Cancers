#!/usr/bin/env python3
"""
Quick diagnostic: Check A3A's DIFF correlations with Known A3-Interactors
and report whether they share direct edges or connect through shared neighbors.
"""
import pickle
import numpy as np
import pandas as pd

DIFF_PATH = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl"
HIGH_PATH = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_HIGH.pkl"
LOW_PATH  = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_LOW.pkl"

INTERACTORS = ["HNRNPA2B1", "HSPD1", "RPL5", "TIMM8B", "RPL3"]
A3_GENES = ["APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H"]
THRESHOLD = 0.40

print("Loading correlation matrices...")
with open(DIFF_PATH, "rb") as f: corr_diff = pickle.load(f)
with open(HIGH_PATH, "rb") as f: corr_high = pickle.load(f)
with open(LOW_PATH, "rb") as f: corr_low = pickle.load(f)
print(f"Matrix shape: {corr_diff.shape}")

# 1. Direct correlations: A3A with each interactor
print("\n" + "="*70)
print("A3A DIRECT CORRELATIONS WITH KNOWN A3-INTERACTORS")
print("="*70)
print(f"{'Gene':<15} {'DIFF':>8} {'HIGH':>8} {'LOW':>8} {'Direction':>10} {'Edge@0.40':>10}")
print("-"*70)
for g in INTERACTORS:
    if g in corr_diff.index and "APOBEC3A" in corr_diff.index:
        d = corr_diff.loc["APOBEC3A", g]
        h = corr_high.loc["APOBEC3A", g]
        l = corr_low.loc["APOBEC3A", g]
        direction = "GAINED" if d > 0 else "LOST"
        edge = "YES" if abs(d) >= THRESHOLD else "no"
        print(f"{g:<15} {d:>+8.4f} {h:>+8.4f} {l:>+8.4f} {direction:>10} {edge:>10}")
    else:
        print(f"{g:<15} NOT IN MATRIX")

# 2. A3B with each interactor (for comparison)
print("\n" + "="*70)
print("A3B DIRECT CORRELATIONS WITH KNOWN A3-INTERACTORS")
print("="*70)
print(f"{'Gene':<15} {'DIFF':>8} {'HIGH':>8} {'LOW':>8}")
print("-"*70)
for g in INTERACTORS:
    if g in corr_diff.index and "APOBEC3B" in corr_diff.index:
        d = corr_diff.loc["APOBEC3B", g]
        h = corr_high.loc["APOBEC3B", g]
        l = corr_low.loc["APOBEC3B", g]
        print(f"{g:<15} {d:>+8.4f} {h:>+8.4f} {l:>+8.4f}")

# 3. Shared neighbors: genes connected to BOTH A3A and each interactor at threshold
print("\n" + "="*70)
print("SHARED NEIGHBORS (genes with |delta-rho| >= 0.40 to BOTH A3A and interactor)")
print("="*70)
a3a_row = corr_diff.loc["APOBEC3A"].abs()
a3a_neighbors = set(a3a_row[a3a_row >= THRESHOLD].index) - {"APOBEC3A"}
print(f"A3A direct neighbors at |delta-rho| >= {THRESHOLD}: {len(a3a_neighbors)}")
if a3a_neighbors:
    print(f"  Neighbors: {', '.join(sorted(a3a_neighbors))}")

for g in INTERACTORS:
    if g not in corr_diff.index:
        continue
    g_row = corr_diff.loc[g].abs()
    g_neighbors = set(g_row[g_row >= THRESHOLD].index) - {g}
    shared = a3a_neighbors & g_neighbors
    print(f"\n{g} neighbors at threshold: {len(g_neighbors)}")
    print(f"  Shared with A3A: {len(shared)}")
    if shared:
        print(f"  Shared genes: {', '.join(sorted(shared))}")

# 4. Top 10 A3A neighbors with their interactor status
print("\n" + "="*70)
print("A3A TOP 20 DIFF PARTNERS (with interactor/community annotation)")
print("="*70)

# Load community assignments
part = pd.read_csv("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/04_communities/SC_best_partition.csv")
g2c = dict(zip(part["gene"], part["community"]))

a3a_abs = corr_diff.loc["APOBEC3A"].abs().drop("APOBEC3A").sort_values(ascending=False)
print(f"{'Rank':<5} {'Gene':<20} {'DIFF':>8} {'HIGH':>8} {'LOW':>8} {'Comm':>5} {'Interactor':>12}")
print("-"*75)
for rank, (gene, val) in enumerate(a3a_abs.head(20).items(), 1):
    d = corr_diff.loc["APOBEC3A", gene]
    h = corr_high.loc["APOBEC3A", gene]
    l = corr_low.loc["APOBEC3A", gene]
    comm = g2c.get(gene, "-")
    is_int = "A3-INTERACTOR" if gene in INTERACTORS else ""
    print(f"{rank:<5} {gene:<20} {d:>+8.4f} {h:>+8.4f} {l:>+8.4f} {str(comm):>5} {is_int:>12}")

# 5. Pairwise correlations among the 5 interactors themselves
print("\n" + "="*70)
print("PAIRWISE DIFF CORRELATIONS AMONG A3-INTERACTORS")
print("="*70)
present = [g for g in INTERACTORS if g in corr_diff.index]
for i, g1 in enumerate(present):
    for g2 in present[i+1:]:
        d = corr_diff.loc[g1, g2]
        h = corr_high.loc[g1, g2]
        print(f"  {g1} <-> {g2}: DIFF={d:+.4f}, HIGH={h:+.4f}")

print("\n" + "="*70)
print("DONE")
print("="*70)
