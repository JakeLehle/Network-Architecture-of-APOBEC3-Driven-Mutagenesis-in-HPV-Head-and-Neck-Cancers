#!/usr/bin/env python3
"""
Sweep_DIFF_Threshold.py

Diagnostic script to find the optimal DIFF network threshold.
Reads the already-computed TOP and BOTTOM correlation matrices from Step04,
then sweeps multiple DIFF thresholds to show how edge count, node count,
and community structure change.

Run AFTER Step04 (correlation matrices must already exist).

Usage:
    conda run -n NETWORK python Sweep_DIFF_Threshold.py
"""

import os
import json
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from network_config import (
    CANCER_TYPES, A3_GENES, A3_ID_TO_ALIAS,
    CORR_THRESHOLD, COMMUNITY_BASE_SEED,
    DIR_01_CLEANED, DIR_04_NETWORKS, FIG2_ROOT,
    banner, log, ensure_dir
)


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
    import leidenalg
    import igraph as ig

    mapping = {n: i for i, n in enumerate(G.nodes())}
    reverse = {i: n for n, i in mapping.items()}

    ig_edges = [(mapping[u], mapping[v]) for u, v in G.edges()]
    ig_weights = [abs(float(d.get("abs_weight", d.get("weight", 1.0))))
                  for _, _, d in G.edges(data=True)]
    ig_weights = [max(w, 1e-10) for w in ig_weights]

    g = ig.Graph(n=len(mapping), edges=ig_edges, directed=False)
    g.es["weight"] = ig_weights

    partition = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=resolution, seed=seed
    )

    cm = {}
    for cid, members in enumerate(partition):
        for idx in members:
            cm[reverse[idx]] = cid
    return cm


def symbol_or_self(ensg, mapping):
    return mapping.get(str(ensg), str(ensg))


# =============================================================================
# CONFIG
# =============================================================================
DIFF_THRESHOLDS = [0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]
LEIDEN_RESOLUTIONS = [0.6, 0.8, 1.0, 1.2]

# Load symbol mapping
with open(os.path.join(DIR_01_CLEANED, "ensg_to_symbol.json")) as f:
    ensg_to_symbol = json.load(f)


# =============================================================================
# SWEEP
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"DIFF Threshold Sweep — {cancer_type}", char="=")

    # Load correlation matrices
    corr_dir = os.path.join(DIR_04_NETWORKS, cancer_type, "corr_matrices")
    corr_top = pd.read_pickle(os.path.join(corr_dir, f"{cancer_type}_corr_TOP.pkl"))
    corr_bot = pd.read_pickle(os.path.join(corr_dir, f"{cancer_type}_corr_BOTTOM.pkl"))
    corr_diff = corr_top - corr_bot

    log(f"Correlation matrix shape: {corr_diff.shape}")

    # Output directory
    sweep_dir = ensure_dir(os.path.join(FIG2_ROOT, "DIFF_threshold_sweep", cancer_type))

    # ---- Sweep thresholds
    sweep_rows = []

    for dt in DIFF_THRESHOLDS:
        banner(f"DIFF threshold = {dt:.2f}")

        G_diff = build_weighted_graph(corr_diff, dt)

        # Remove isolates
        isolates = [n for n in G_diff.nodes() if G_diff.degree(n) == 0]
        G_diff_noiso = G_diff.copy()
        G_diff_noiso.remove_nodes_from(isolates)

        n_nodes = G_diff_noiso.number_of_nodes()
        n_edges = G_diff_noiso.number_of_edges()
        avg_deg = 2 * n_edges / n_nodes if n_nodes > 0 else 0
        density = nx.density(G_diff_noiso) if n_nodes > 1 else 0

        log(f"  Nodes (non-isolated): {n_nodes}")
        log(f"  Edges: {n_edges}")
        log(f"  Avg degree: {avg_deg:.1f}")
        log(f"  Density: {density:.4f}")

        # LCC
        if n_nodes > 0:
            lcc = max(nx.connected_components(G_diff_noiso), key=len)
            lcc_size = len(lcc)
            n_components = nx.number_connected_components(G_diff_noiso)
        else:
            lcc_size = 0
            n_components = 0

        log(f"  Connected components: {n_components}")
        log(f"  LCC size: {lcc_size}")

        # Quick community detection at multiple resolutions
        if n_nodes >= 5 and n_edges >= 3:
            G_lcc = G_diff_noiso.subgraph(lcc).copy() if lcc_size > 0 else G_diff_noiso

            # Ensure abs_weight
            for u, v, d in G_lcc.edges(data=True):
                d["abs_weight"] = abs(float(d.get("weight", 0.0)))

            for res in LEIDEN_RESOLUTIONS:
                cm = quick_leiden(G_lcc, resolution=res, seed=COMMUNITY_BASE_SEED)
                n_comms = len(set(cm.values()))

                # Community sizes
                comm_sizes = {}
                for n, c in cm.items():
                    comm_sizes[c] = comm_sizes.get(c, 0) + 1
                sizes_sorted = sorted(comm_sizes.values(), reverse=True)

                # Modularity
                from networkx.algorithms.community.quality import modularity as nx_mod
                comm_sets = {}
                for n, c in cm.items():
                    comm_sets.setdefault(c, set()).add(n)
                mod = nx_mod(G_lcc, list(comm_sets.values()), weight="abs_weight")

                log(f"  Leiden (res={res:.1f}): {n_comms} communities, modularity={mod:.4f}")
                log(f"    Sizes: {sizes_sorted[:10]}{'...' if len(sizes_sorted) > 10 else ''}")

                sweep_rows.append({
                    "diff_threshold": dt,
                    "leiden_resolution": res,
                    "n_nodes": n_nodes,
                    "n_edges": n_edges,
                    "avg_degree": avg_deg,
                    "density": density,
                    "n_components": n_components,
                    "lcc_size": lcc_size,
                    "n_communities": n_comms,
                    "modularity": mod,
                    "largest_comm": sizes_sorted[0] if sizes_sorted else 0,
                    "smallest_comm": sizes_sorted[-1] if sizes_sorted else 0,
                    "comm_sizes": str(sizes_sorted[:10]),
                })

                # If we found good structure, identify top hub genes at this config
                if n_comms >= 3 and mod > 0.15:
                    log(f"\n    >>> Good structure at diff={dt:.2f}, res={res:.1f}")
                    log(f"    >>> Top genes per community:")
                    for c_id in sorted(comm_sets.keys()):
                        c_nodes = list(comm_sets[c_id])
                        if len(c_nodes) < 3:
                            continue
                        # Top 5 by degree
                        deg = {n: G_lcc.degree(n) for n in c_nodes}
                        top5 = sorted(deg.items(), key=lambda x: x[1], reverse=True)[:5]
                        names = [f"{symbol_or_self(g, ensg_to_symbol)}({d})" for g, d in top5]
                        log(f"      C{c_id} (n={len(c_nodes)}): {', '.join(names)}")
        else:
            log("  [SKIP] Graph too small for community detection")

    # ---- Save sweep table
    if sweep_rows:
        sweep_df = pd.DataFrame(sweep_rows)
        csv_path = os.path.join(sweep_dir, f"{cancer_type}_DIFF_threshold_sweep.csv")
        sweep_df.to_csv(csv_path, index=False)
        log(f"\n[SAVE] Sweep table -> {csv_path}")

        # ---- Summary plot: edges & communities vs threshold
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Plot 1: Edges vs threshold
        for res in LEIDEN_RESOLUTIONS:
            sub = sweep_df[sweep_df["leiden_resolution"] == res]
            axes[0, 0].plot(sub["diff_threshold"], sub["n_edges"], "o-", label=f"res={res}")
        axes[0, 0].set_xlabel("DIFF Threshold")
        axes[0, 0].set_ylabel("Number of Edges")
        axes[0, 0].set_title("Edge Count vs DIFF Threshold")
        axes[0, 0].set_yscale("log")
        axes[0, 0].legend()

        # Plot 2: Communities vs threshold
        for res in LEIDEN_RESOLUTIONS:
            sub = sweep_df[sweep_df["leiden_resolution"] == res]
            axes[0, 1].plot(sub["diff_threshold"], sub["n_communities"], "o-", label=f"res={res}")
        axes[0, 1].set_xlabel("DIFF Threshold")
        axes[0, 1].set_ylabel("Number of Communities")
        axes[0, 1].set_title("Communities vs DIFF Threshold")
        axes[0, 1].legend()

        # Plot 3: Modularity vs threshold
        for res in LEIDEN_RESOLUTIONS:
            sub = sweep_df[sweep_df["leiden_resolution"] == res]
            axes[1, 0].plot(sub["diff_threshold"], sub["modularity"], "o-", label=f"res={res}")
        axes[1, 0].set_xlabel("DIFF Threshold")
        axes[1, 0].set_ylabel("Modularity")
        axes[1, 0].set_title("Modularity vs DIFF Threshold")
        axes[1, 0].legend()

        # Plot 4: Avg degree vs threshold
        for res in LEIDEN_RESOLUTIONS:
            sub = sweep_df[sweep_df["leiden_resolution"] == res]
            axes[1, 1].plot(sub["diff_threshold"], sub["avg_degree"], "o-", label=f"res={res}")
        axes[1, 1].set_xlabel("DIFF Threshold")
        axes[1, 1].set_ylabel("Average Degree")
        axes[1, 1].set_title("Avg Degree vs DIFF Threshold")
        axes[1, 1].set_yscale("log")
        axes[1, 1].legend()

        plt.suptitle(f"{cancer_type} — DIFF Threshold Sweep", fontsize=14)
        plt.tight_layout()
        plot_path = os.path.join(sweep_dir, f"{cancer_type}_threshold_sweep_summary.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        log(f"[SAVE] Summary plot -> {plot_path}")

    banner("SWEEP COMPLETE")
