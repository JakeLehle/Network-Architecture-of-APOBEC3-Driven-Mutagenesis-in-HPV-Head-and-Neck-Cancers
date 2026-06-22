#!/usr/bin/env python3
"""
Step02.1_SC_Sweep_DIFF_Threshold.py
=====================================

Diagnostic script for detailed DIFF network threshold exploration.
Reads the already-computed HIGH and LOW correlation matrices from Step 02,
then sweeps DIFF thresholds with finer granularity and includes Leiden
community counts at each threshold.

Run AFTER Step 02 (correlation matrices must already exist).

This is the single-cell analog of Step04.1_Sweep_DIFF_Threshold.py from
the Figure 2 pipeline.

Input:
  data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_HIGH.pkl
  data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_LOW.pkl

Output (→ data/FIG_4/03_correlation_networks/):
  sweep_detailed.csv         — full sweep results
  sweep_detailed_plot.png    — multi-panel visualization

Usage:
  conda run -n NETWORK python Step02.1_SC_Sweep_DIFF_Threshold.py

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
    DIR_03_NETWORKS, COMMUNITY_BASE_SEED,
    A3_GENES, A3_ID_TO_ALIAS,
    banner, log, ensure_dir
)


# =============================================================================
# SWEEP PARAMETERS
# =============================================================================
SWEEP_RANGE = np.arange(0.30, 0.95, 0.05)
LEIDEN_RESOLUTIONS = [0.6, 0.8, 1.0]


# =============================================================================
# HELPERS
# =============================================================================

def build_weighted_graph(corr_df, threshold):
    """Build graph from correlation matrix with |corr| >= threshold."""
    genes = list(corr_df.index)
    G = nx.Graph()
    G.add_nodes_from(genes)
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            w = float(corr_df.iat[i, j])
            if abs(w) >= threshold and abs(w) < 0.999999:
                G.add_edge(genes[i], genes[j], weight=w, abs_weight=abs(w))
    return G


def quick_leiden(G, resolution=1.0, seed=42):
    """Run Leiden once and return {node: community} dict."""
    try:
        import leidenalg
        import igraph as ig
    except ImportError:
        log("WARNING: leidenalg not available, skipping community detection")
        return {}

    mapping = {n: i for i, n in enumerate(G.nodes())}
    reverse = {i: n for n, i in mapping.items()}

    ig_edges = []
    ig_weights = []
    for u, v, d in G.edges(data=True):
        ig_edges.append((mapping[u], mapping[v]))
        w = abs(float(d.get("abs_weight", d.get("weight", 1.0))))
        ig_weights.append(max(w, 1e-10))

    g = ig.Graph(n=len(mapping), edges=ig_edges, directed=False)
    g.es["weight"] = ig_weights

    partition = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=resolution, seed=seed
    )

    comm_map = {}
    for comm_id, members in enumerate(partition):
        for node_idx in members:
            comm_map[reverse[node_idx]] = comm_id
    return comm_map


# =============================================================================
# MAIN
# =============================================================================

def main():
    start_time = datetime.now()
    banner("[STEP 02.1] Detailed DIFF Threshold Sweep")
    log(f"Start time: {start_time}")

    # Load correlation matrices
    corr_dir = os.path.join(DIR_03_NETWORKS, "corr_matrices")

    high_path = os.path.join(corr_dir, "SC_corr_HIGH.pkl")
    low_path = os.path.join(corr_dir, "SC_corr_LOW.pkl")

    log(f"Loading HIGH correlation matrix: {high_path}")
    with open(high_path, "rb") as f:
        corr_high = pickle.load(f)

    log(f"Loading LOW correlation matrix: {low_path}")
    with open(low_path, "rb") as f:
        corr_low = pickle.load(f)

    corr_diff = corr_high - corr_low
    log(f"DIFF matrix: {corr_diff.shape}")

    # Sweep
    banner("[SWEEP] Threshold scan with community detection")

    results = []
    header = (f"{'Threshold':>10} {'Nodes':>7} {'Edges':>8} {'AvgDeg':>7} "
              f"{'Density':>10} {'Comp':>5} {'LCC':>6}")

    for res in LEIDEN_RESOLUTIONS:
        header += f" {'C(r='+f'{res:.1f}'+')':>10}"
    print(header)
    print("-" * len(header))

    for dt in SWEEP_RANGE:
        dt = round(dt, 2)
        G = build_weighted_graph(corr_diff, dt)

        # Remove isolates
        isolates = [n for n in G.nodes() if G.degree(n) == 0]
        G_noiso = G.copy()
        G_noiso.remove_nodes_from(isolates)

        n_nodes = G_noiso.number_of_nodes()
        n_edges = G_noiso.number_of_edges()
        avg_deg = 2 * n_edges / n_nodes if n_nodes > 0 else 0
        density = nx.density(G_noiso) if n_nodes > 1 else 0

        if n_nodes > 0:
            components = list(nx.connected_components(G_noiso))
            n_comp = len(components)
            lcc = max(components, key=len)
            lcc_size = len(lcc)
        else:
            n_comp = 0
            lcc_size = 0

        row = {
            "threshold": dt, "nodes": n_nodes, "edges": n_edges,
            "avg_degree": avg_deg, "density": density,
            "components": n_comp, "lcc_size": lcc_size
        }

        line = (f"{dt:>10.2f} {n_nodes:>7} {n_edges:>8} {avg_deg:>7.1f} "
                f"{density:>10.6f} {n_comp:>5} {lcc_size:>6}")

        # Community detection at each resolution (on LCC)
        if n_nodes >= 5 and n_edges >= 3:
            G_lcc = G_noiso.subgraph(lcc).copy() if lcc_size > 0 else G_noiso
            for res in LEIDEN_RESOLUTIONS:
                cm = quick_leiden(G_lcc, resolution=res, seed=COMMUNITY_BASE_SEED)
                n_comms = len(set(cm.values())) if cm else 0
                row[f"comms_r{res}"] = n_comms
                line += f" {n_comms:>10}"
        else:
            for res in LEIDEN_RESOLUTIONS:
                row[f"comms_r{res}"] = 0
                line += f" {'N/A':>10}"

        # Check A3 genes in network
        a3_in = sum(1 for a3 in A3_GENES if a3 in G_noiso.nodes())
        row["a3_genes_in_network"] = a3_in

        results.append(row)
        print(line)

    sweep_df = pd.DataFrame(results)

    # Save
    csv_path = os.path.join(DIR_03_NETWORKS, "sweep_detailed.csv")
    sweep_df.to_csv(csv_path, index=False)
    log(f"\n[SAVE] Sweep CSV → {csv_path}")

    # Plot
    banner("[PLOT] Threshold sweep visualization")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    thresholds = sweep_df["threshold"]

    # Panel 1: Nodes and edges
    ax1 = axes[0, 0]
    ax1.plot(thresholds, sweep_df["nodes"], "o-", c="steelblue", label="Nodes")
    ax1b = ax1.twinx()
    ax1b.plot(thresholds, sweep_df["edges"], "s-", c="firebrick", label="Edges")
    ax1.set_xlabel("|Δρ| threshold")
    ax1.set_ylabel("Nodes", color="steelblue")
    ax1b.set_ylabel("Edges", color="firebrick")
    ax1.set_title("Network Size vs Threshold")

    # Panel 2: Average degree
    ax2 = axes[0, 1]
    ax2.plot(thresholds, sweep_df["avg_degree"], "o-", c="darkgreen")
    ax2.set_xlabel("|Δρ| threshold")
    ax2.set_ylabel("Average Degree")
    ax2.set_title("Average Degree vs Threshold")

    # Panel 3: LCC size
    ax3 = axes[1, 0]
    ax3.plot(thresholds, sweep_df["lcc_size"], "o-", c="purple")
    ax3.set_xlabel("|Δρ| threshold")
    ax3.set_ylabel("LCC Size (nodes)")
    ax3.set_title("Largest Connected Component vs Threshold")

    # Panel 4: Community counts
    ax4 = axes[1, 1]
    for res in LEIDEN_RESOLUTIONS:
        col = f"comms_r{res}"
        if col in sweep_df.columns:
            ax4.plot(thresholds, sweep_df[col], "o-", label=f"res={res}")
    ax4.set_xlabel("|Δρ| threshold")
    ax4.set_ylabel("Communities (Leiden)")
    ax4.set_title("Community Count vs Threshold")
    ax4.legend(fontsize=8)

    plt.tight_layout()
    plot_path = os.path.join(DIR_03_NETWORKS, "sweep_detailed_plot.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    log(f"[SAVE] Sweep plot → {plot_path}")

    elapsed = datetime.now() - start_time
    banner(f"[STEP 02.1 COMPLETE] Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
