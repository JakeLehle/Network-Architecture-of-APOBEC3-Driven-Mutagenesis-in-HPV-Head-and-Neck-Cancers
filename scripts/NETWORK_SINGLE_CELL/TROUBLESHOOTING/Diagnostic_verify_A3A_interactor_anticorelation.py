#!/usr/bin/env python3
"""
Diagnostic_Verify_A3A_Interactor_Anticorrelation.py
=====================================================

Quick sanity check: are the 4 C0 interactors (HNRNPA2B1, HSPD1, RPL5, TIMM8B)
truly anti-correlated with APOBEC3A in the DIFF network?

Checks three independent sources:
  1. Raw DIFF correlation matrix (SC_corr_DIFF.pkl)
     -> rho_DIFF = rho_HIGH - rho_LOW for each A3A-interactor pair
  2. HIGH and LOW correlation matrices separately
     -> Shows whether the anti-correlation is from LOST co-expression
        (positive in LOW, near-zero in HIGH) vs GAINED anti-correlation
  3. Graph edge weights in SC_G_comm.gpickle
     -> Confirms what the network plot is actually drawing

Also reports:
  - Whether A3A has a direct edge to each interactor (or indirect)
  - The sign and magnitude of each connection
  - Interactor-to-interactor correlations (they should be positively
    co-expressed with each other per Quick_A3A_Interactor_Correlations.py)

Usage:
  conda run -n NETWORK python Diagnostic_Verify_A3A_Interactor_Anticorrelation.py

Author: Jake Lehle / Claude (2026 NMF Paper, grant QC)
"""

import os, pickle, sys
import numpy as np
import pandas as pd
import networkx as nx

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
CORR_DIR = os.path.join(BASE_DIR, "data/FIG_4/03_correlation_networks/corr_matrices")
SC_GRAPH = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_G_comm.gpickle")
SC_PARTITION = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_best_partition.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "data/FIG_4/TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# The 4 interactors in C0
C0_INTERACTORS = ["HNRNPA2B1", "HSPD1", "RPL5", "TIMM8B"]
A3A = "APOBEC3A"

# =============================================================================
# LOGGING
# =============================================================================
report = []
def log(msg=""):
    print(msg, flush=True)
    report.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 70)
    log(f"  {title}")
    log(char * 70)

# =============================================================================
# 1. LOAD CORRELATION MATRICES
# =============================================================================
banner("1. LOAD CORRELATION MATRICES")

corr_files = {
    "DIFF": os.path.join(CORR_DIR, "SC_corr_DIFF.pkl"),
    "HIGH": os.path.join(CORR_DIR, "SC_corr_HIGH.pkl"),
    "LOW":  os.path.join(CORR_DIR, "SC_corr_LOW.pkl"),
}

corr = {}
for label, path in corr_files.items():
    if os.path.exists(path):
        with open(path, "rb") as f:
            corr[label] = pickle.load(f)
        log(f"  {label}: {corr[label].shape[0]} x {corr[label].shape[1]}")
    else:
        log(f"  WARNING: {label} matrix not found at {path}")

# Check that A3A and all interactors are in the matrices
for label, mat in corr.items():
    genes_in = [A3A] + C0_INTERACTORS
    missing = [g for g in genes_in if g not in mat.index]
    if missing:
        log(f"  WARNING: {label} matrix missing genes: {missing}")
    else:
        log(f"  {label}: all target genes present")

# =============================================================================
# 2. A3A <-> INTERACTOR CORRELATIONS (all three matrices)
# =============================================================================
banner("2. A3A <-> INTERACTOR CORRELATIONS")

log(f"\n  {'Gene':<15s} {'rho_DIFF':>10s} {'rho_HIGH':>10s} {'rho_LOW':>10s} {'Interpretation'}")
log(f"  {'-'*15} {'-'*10} {'-'*10} {'-'*10} {'-'*30}")

for gene in C0_INTERACTORS:
    vals = {}
    for label in ["DIFF", "HIGH", "LOW"]:
        if label in corr and A3A in corr[label].index and gene in corr[label].columns:
            vals[label] = corr[label].loc[A3A, gene]
        else:
            vals[label] = np.nan

    # Interpret the pattern
    if vals["DIFF"] < -0.1:
        if vals["LOW"] > 0.1 and vals["HIGH"] < 0.1:
            interp = "LOST co-expression"
        elif vals["HIGH"] < -0.1:
            interp = "GAINED anti-correlation"
        else:
            interp = "Mixed (DIFF negative)"
    elif vals["DIFF"] > 0.1:
        interp = "GAINED co-expression (!)"
    else:
        interp = "Weak/no DIFF signal"

    log(f"  {gene:<15s} {vals['DIFF']:>+10.4f} {vals['HIGH']:>+10.4f} {vals['LOW']:>+10.4f} {interp}")

# =============================================================================
# 3. GRAPH EDGE WEIGHTS (what the network plot actually draws)
# =============================================================================
banner("3. GRAPH EDGE WEIGHTS")

if os.path.exists(SC_GRAPH):
    with open(SC_GRAPH, "rb") as f:
        G_full = pickle.load(f)
    log(f"  Full graph: {G_full.number_of_nodes()} nodes, {G_full.number_of_edges()} edges")

    # Check partition to confirm these are all in C0
    part_df = pd.read_csv(SC_PARTITION)
    g2c = dict(zip(part_df["gene"], part_df["community"]))

    log(f"\n  Community assignments:")
    for gene in [A3A] + C0_INTERACTORS:
        c = g2c.get(gene, "NOT FOUND")
        in_graph = gene in G_full.nodes()
        log(f"    {gene:<15s}  community={c}  in_graph={in_graph}")

    # Direct edges from A3A to each interactor
    log(f"\n  Direct edges A3A <-> interactors:")
    log(f"  {'Gene':<15s} {'Has Edge':>10s} {'Weight':>10s} {'|Weight|':>10s} {'Sign':>8s}")
    log(f"  {'-'*15} {'-'*10} {'-'*10} {'-'*10} {'-'*8}")

    for gene in C0_INTERACTORS:
        if G_full.has_edge(A3A, gene):
            w = G_full[A3A][gene].get("weight", np.nan)
            sign = "POS" if w > 0 else "NEG"
            log(f"  {gene:<15s} {'YES':>10s} {w:>+10.4f} {abs(w):>10.4f} {sign:>8s}")
        else:
            log(f"  {gene:<15s} {'NO':>10s} {'---':>10s} {'---':>10s} {'---':>8s}")

    # Also check shortest path if no direct edge
    log(f"\n  Shortest paths A3A <-> interactors (if no direct edge):")
    c0_genes = set(part_df[part_df["community"] == 0]["gene"])
    c0_in_graph = [n for n in c0_genes if n in G_full.nodes()]
    G_c0 = G_full.subgraph(c0_in_graph).copy()

    for gene in C0_INTERACTORS:
        if not G_full.has_edge(A3A, gene):
            if A3A in G_c0 and gene in G_c0:
                try:
                    path = nx.shortest_path(G_c0, A3A, gene)
                    path_str = " -> ".join(path)
                    log(f"    {gene}: path length {len(path)-1}: {path_str}")

                    # Report edge signs along the path
                    for i in range(len(path) - 1):
                        u, v = path[i], path[i+1]
                        w = G_c0[u][v].get("weight", np.nan)
                        log(f"      {u} -> {v}: weight={w:+.4f} ({'POS' if w>0 else 'NEG'})")
                except nx.NetworkXNoPath:
                    log(f"    {gene}: NO PATH in C0 subgraph")
            else:
                log(f"    {gene}: not reachable in C0 subgraph")
        else:
            log(f"    {gene}: has direct edge (see above)")
else:
    log(f"  WARNING: Graph not found at {SC_GRAPH}")

# =============================================================================
# 4. INTERACTOR <-> INTERACTOR CORRELATIONS
# =============================================================================
banner("4. INTERACTOR <-> INTERACTOR CORRELATIONS")

log(f"\n  Checking whether interactors are positively co-expressed with each other:")
log(f"  {'Pair':<30s} {'rho_DIFF':>10s} {'rho_HIGH':>10s} {'rho_LOW':>10s}")
log(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*10}")

for i, g1 in enumerate(C0_INTERACTORS):
    for g2 in C0_INTERACTORS[i+1:]:
        vals = {}
        for label in ["DIFF", "HIGH", "LOW"]:
            if label in corr and g1 in corr[label].index and g2 in corr[label].columns:
                vals[label] = corr[label].loc[g1, g2]
            else:
                vals[label] = np.nan
        pair = f"{g1} <-> {g2}"
        log(f"  {pair:<30s} {vals['DIFF']:>+10.4f} {vals['HIGH']:>+10.4f} {vals['LOW']:>+10.4f}")

# =============================================================================
# 5. SUMMARY VERDICT
# =============================================================================
banner("5. SUMMARY VERDICT")

if "DIFF" in corr:
    all_neg = True
    for gene in C0_INTERACTORS:
        if A3A in corr["DIFF"].index and gene in corr["DIFF"].columns:
            rho = corr["DIFF"].loc[A3A, gene]
            if rho >= 0:
                log(f"  *** {gene}: rho_DIFF = {rho:+.4f} -- NOT anti-correlated! ***")
                all_neg = False

    if all_neg:
        log(f"  CONFIRMED: All 4 C0 interactors have negative DIFF correlation with A3A.")
        log(f"  The anti-correlation claim is supported by the data.")
    else:
        log(f"  WARNING: At least one interactor is NOT anti-correlated with A3A in DIFF.")
        log(f"  Review the numbers above carefully before making this claim in the grant.")

    # Special attention to TIMM8B (Jake flagged this one)
    log("")
    if A3A in corr["DIFF"].index and "TIMM8B" in corr["DIFF"].columns:
        rho_diff = corr["DIFF"].loc[A3A, "TIMM8B"]
        rho_high = corr["HIGH"].loc[A3A, "TIMM8B"] if "HIGH" in corr else np.nan
        rho_low  = corr["LOW"].loc[A3A, "TIMM8B"] if "LOW" in corr else np.nan
        log(f"  TIMM8B SPOTLIGHT (flagged by Jake):")
        log(f"    rho_DIFF = {rho_diff:+.4f}")
        log(f"    rho_HIGH = {rho_high:+.4f}")
        log(f"    rho_LOW  = {rho_low:+.4f}")
        if rho_diff < 0:
            log(f"    -> DIFF is negative. TIMM8B qualifies as anti-correlated.")
            if not G_full.has_edge(A3A, "TIMM8B"):
                log(f"    -> But NOTE: no direct graph edge A3A<->TIMM8B.")
                log(f"       TIMM8B may be connected through intermediate genes,")
                log(f"       which explains the distal position in the layout.")
                log(f"       Anti-correlation still holds in the correlation matrix")
                log(f"       even without a direct edge in the thresholded network.")
        else:
            log(f"    -> DIFF is NOT negative. Double-check this interactor.")

# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUTPUT_DIR, "verify_A3A_interactor_anticorrelation.txt")
with open(report_path, "w") as f:
    f.write("\n".join(report))
log(f"\n  Report saved: {report_path}")

banner("DONE")
