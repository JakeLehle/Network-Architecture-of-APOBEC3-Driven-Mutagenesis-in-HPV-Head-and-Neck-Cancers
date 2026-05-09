#!/usr/bin/env python3
"""
Step04_SC_Centrality_Metrics.py
================================

Figure 4 — Step 04: Per-gene centrality metrics for the single-cell
HIGH, LOW, and DIFF co-expression networks.

Mirrors Step06_Centrality_Metrics.py from the Figure 2 (TCGA bulk) pipeline.

Computes: degree, betweenness centrality, closeness centrality,
eigenvector centrality, and strength (sum of |edge weights|) for each gene.

Input:
  data/FIG_4/03_correlation_networks/graph_objects/SC_G_*.gpickle
  data/FIG_4/04_communities/SC_best_partition.csv

Output (→ data/FIG_4/05_centrality_metrics/):
  SC_{HIGH|LOW|DIFF}_metrics.csv  — per-gene centrality tables
  SC_hub_genes_DIFF.csv           — top hub genes by multiple metrics

Usage:
  conda run -n NETWORK python Step04_SC_Centrality_Metrics.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime

from network_config_SC import (
    DIR_03_NETWORKS, DIR_04_COMMUNITIES, DIR_05_CENTRALITY,
    A3_GENES, A3_ID_TO_ALIAS,
    banner, log, ensure_dir, load_ensg_to_symbol
)


# =============================================================================
# CENTRALITY COMPUTATION
# =============================================================================

def compute_centrality(G, label=""):
    """Compute per-node centrality metrics for graph G."""
    log(f"  Computing centrality for {label} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")

    if G.number_of_nodes() == 0:
        return pd.DataFrame()

    # Degree
    degree = dict(G.degree())

    # Strength (sum of |edge weights|)
    strength = {}
    for n in G.nodes():
        s = sum(abs(float(G[n][nb].get("abs_weight", G[n][nb].get("weight", 1.0))))
                for nb in G.neighbors(n))
        strength[n] = s

    # Betweenness centrality
    log(f"    Betweenness...")
    betweenness = nx.betweenness_centrality(G, weight="abs_weight")

    # Closeness centrality
    log(f"    Closeness...")
    closeness = nx.closeness_centrality(G, distance=None)

    # Eigenvector centrality
    log(f"    Eigenvector...")
    try:
        eigenvector = nx.eigenvector_centrality(G, max_iter=1000, weight="abs_weight")
    except nx.PowerIterationFailedConvergence:
        log(f"    WARNING: Eigenvector centrality did not converge, using zeros")
        eigenvector = {n: 0.0 for n in G.nodes()}

    # Combine
    metrics = pd.DataFrame({
        "gene": list(G.nodes()),
        "degree": [degree[n] for n in G.nodes()],
        "strength": [strength[n] for n in G.nodes()],
        "betweenness": [betweenness[n] for n in G.nodes()],
        "closeness": [closeness[n] for n in G.nodes()],
        "eigenvector": [eigenvector[n] for n in G.nodes()],
    })

    # Rank (1 = highest for each metric)
    for col in ["degree", "strength", "betweenness", "closeness", "eigenvector"]:
        metrics[f"{col}_rank"] = metrics[col].rank(ascending=False, method="min").astype(int)

    # Composite hub score: mean of ranks
    rank_cols = [c for c in metrics.columns if c.endswith("_rank")]
    metrics["hub_score"] = metrics[rank_cols].mean(axis=1)
    metrics["hub_rank"] = metrics["hub_score"].rank(ascending=True, method="min").astype(int)

    return metrics.sort_values("hub_rank")


# =============================================================================
# MAIN
# =============================================================================

def main():
    start_time = datetime.now()
    banner("[STEP 04] Single-Cell Centrality Metrics")
    log(f"Start time: {start_time}")

    ensure_dir(DIR_05_CENTRALITY)

    ensg_to_symbol = load_ensg_to_symbol()

    # Load community assignments
    part_path = os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv")
    if os.path.exists(part_path):
        part_df = pd.read_csv(part_path)
        gene_to_comm = dict(zip(part_df["gene"], part_df["community"]))
        log(f"Loaded community assignments: {len(gene_to_comm)} genes")
    else:
        gene_to_comm = {}
        log("WARNING: No community assignments found")

    # Process each network
    graph_dir = os.path.join(DIR_03_NETWORKS, "graph_objects")

    for net_name in ["HIGH", "LOW", "DIFF"]:
        banner(f"[STEP 04] {net_name} network centrality")

        graph_path = os.path.join(graph_dir, f"SC_G_{net_name.lower()}_noiso.gpickle")
        if not os.path.exists(graph_path):
            log(f"  WARNING: {graph_path} not found, skipping")
            continue

        with open(graph_path, "rb") as f:
            G = pickle.load(f)

        metrics = compute_centrality(G, label=net_name)
        if metrics.empty:
            continue

        # Add annotations
        metrics["symbol"] = metrics["gene"].map(
            lambda g: ensg_to_symbol.get(g, A3_ID_TO_ALIAS.get(g, g))
        )
        metrics["community"] = metrics["gene"].map(lambda g: gene_to_comm.get(g, -1))
        metrics["is_A3"] = metrics["gene"].isin(A3_GENES)

        # Save
        out_path = os.path.join(DIR_05_CENTRALITY, f"SC_{net_name}_metrics.csv")
        metrics.to_csv(out_path, index=False)
        log(f"  [SAVE] {net_name} metrics ({len(metrics)} genes) → {out_path}")

        # Report top 20 hubs
        log(f"\n  Top 20 hub genes ({net_name}):")
        log(f"  {'Rank':>5} {'Symbol':>15} {'Deg':>5} {'Between':>8} {'Eigen':>8} {'Comm':>5}")
        for _, row in metrics.head(20).iterrows():
            sym = row["symbol"]
            if row["is_A3"]:
                sym = f"**{sym}**"
            log(f"  {int(row['hub_rank']):>5} {sym:>15} {int(row['degree']):>5} "
                f"{row['betweenness']:>8.4f} {row['eigenvector']:>8.4f} "
                f"{int(row['community']):>5}")

    # =========================================================================
    # Extract top DIFF hub genes specifically
    # =========================================================================
    banner("[STEP 04] Save DIFF hub genes")

    diff_path = os.path.join(DIR_05_CENTRALITY, "SC_DIFF_metrics.csv")
    if os.path.exists(diff_path):
        diff_metrics = pd.read_csv(diff_path)
        # Top 50 hub genes
        hubs = diff_metrics.head(50).copy()
        hub_path = os.path.join(DIR_05_CENTRALITY, "SC_hub_genes_DIFF.csv")
        hubs.to_csv(hub_path, index=False)
        log(f"  [SAVE] Top 50 DIFF hubs → {hub_path}")

        # A3 genes specifically
        a3_metrics = diff_metrics[diff_metrics["is_A3"]].copy()
        if len(a3_metrics) > 0:
            log(f"\n  A3 gene centrality in DIFF network:")
            for _, row in a3_metrics.iterrows():
                log(f"    {row['symbol']}: hub_rank={int(row['hub_rank'])}, "
                    f"degree={int(row['degree'])}, community={int(row['community'])}")

    elapsed = datetime.now() - start_time
    banner(f"[STEP 04 COMPLETE] Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
