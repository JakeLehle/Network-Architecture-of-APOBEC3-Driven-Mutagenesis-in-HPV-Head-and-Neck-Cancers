#!/usr/bin/env python3
"""
Step05_Community_Detection.py

Run Leiden community detection on the DIFF network across multiple
resolutions. Evaluate stability (ARI/NMI across runs), merge small
communities, and save community assignments + network visualizations.

Corresponds to original pipeline Step 15.

Input:
    05_correlation_networks/{cancer_type}/graph_objects/{ct}_G_diff_noiso.gpickle
    01_cleaned_expression/ensg_to_symbol.json

Output (-> data/FIG_2/06_communities/{cancer_type}/):
    sweep/                     — per-resolution community assignments + plots
    {ct}_resolution_sweep.csv  — stability metrics across resolutions
    {ct}_sweep_plots.png       — modularity/ARI/NMI curves
    {ct}_best_partition.csv    — final community assignments (best resolution)
    {ct}_community_gene_lists.csv  — genes per community

Usage:
    python Step05_Community_Detection.py
"""

import os
import json
import pickle
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from network_config import (
    CANCER_TYPES, A3_GENES, BIOMARKERS, A3_ID_TO_ALIAS,
    COMMUNITY_METHOD, COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION,
    COMMUNITY_BASE_SEED, USE_LARGEST_COMPONENT,
    TARGET_BIG_COMMUNITIES, MIN_COMMUNITY_SIZE,
    DIR_01_CLEANED, DIR_04_NETWORKS, DIR_05_COMMUNITIES,
    banner, log, ensure_dir
)


# =============================================================================
# COMMUNITY DETECTION FUNCTIONS
# =============================================================================

def detect_leiden(G, seed=42, resolution=1.0, weight="abs_weight"):
    """Run Leiden community detection. Returns {node: community_id} dict."""
    try:
        import leidenalg
        import igraph as ig
    except ImportError:
        raise ImportError("leidenalg and python-igraph required. Install: pip install leidenalg igraph")

    # Convert networkx -> igraph
    mapping = {n: i for i, n in enumerate(G.nodes())}
    reverse = {i: n for n, i in mapping.items()}

    ig_edges = []
    ig_weights = []
    for u, v, d in G.edges(data=True):
        ig_edges.append((mapping[u], mapping[v]))
        w = abs(float(d.get(weight, d.get("weight", 1.0))))
        ig_weights.append(max(w, 1e-10))

    g = ig.Graph(n=len(mapping), edges=ig_edges, directed=False)
    g.es["weight"] = ig_weights

    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=seed
    )

    comm_map = {}
    for comm_id, members in enumerate(partition):
        for node_idx in members:
            comm_map[reverse[node_idx]] = comm_id

    return comm_map


def partition_to_labels(comm_map, nodes):
    """Convert {node: community} to label array aligned with node list."""
    return np.array([comm_map.get(n, -1) for n in nodes], dtype=int)


def nx9_largest_component(G):
    """Return largest connected component."""
    if G.number_of_edges() == 0:
        return G.copy()
    comp = max(nx.connected_components(G), key=len)
    return G.subgraph(comp).copy()


def nx9_remove_degree0(G):
    """Remove degree-0 nodes."""
    H = G.copy()
    H.remove_nodes_from([n for n in H.nodes() if H.degree(n) == 0])
    return H


def merge_small_communities(comm_map, k_keep=8, min_size=10):
    """Merge small communities into 'Other' bin."""
    sizes = {}
    for n, cid in comm_map.items():
        sizes[cid] = sizes.get(cid, 0) + 1

    top = sorted(sizes.keys(), key=lambda c: sizes[c], reverse=True)[:k_keep]
    top = set(top)
    other_id = max(sizes.keys()) + 1 if sizes else 0

    merged = {}
    for n, cid in comm_map.items():
        if cid in top and sizes[cid] >= min_size:
            merged[n] = cid
        else:
            merged[n] = other_id

    merged_sizes = {}
    for n, cid in merged.items():
        merged_sizes[cid] = merged_sizes.get(cid, 0) + 1

    return merged, sizes, merged_sizes


def symbol_or_self(ensg, mapping):
    return mapping.get(str(ensg), str(ensg))


# =============================================================================
# LOAD symbol mapping
# =============================================================================
symbol_path = os.path.join(DIR_01_CLEANED, "ensg_to_symbol.json")
with open(symbol_path) as f:
    ensg_to_symbol = json.load(f)


# =============================================================================
# LOOP over cancer types
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"[STEP 15] Community Detection — {cancer_type}", char="=")

    cancer_dir = ensure_dir(os.path.join(DIR_05_COMMUNITIES, cancer_type))
    sweep_dir = ensure_dir(os.path.join(cancer_dir, "sweep"))

    # ---- Load DIFF graph
    pkl_path = os.path.join(DIR_04_NETWORKS, cancer_type, "graph_objects",
                            f"{cancer_type}_G_diff_noiso.gpickle")
    if not os.path.exists(pkl_path):
        log(f"[SKIP] DIFF graph not found: {pkl_path}")
        continue

    with open(pkl_path, "rb") as f:
        G_diff_noiso = pickle.load(f)

    log(f"[STEP 15.0] Loaded DIFF graph: nodes={G_diff_noiso.number_of_nodes()}, edges={G_diff_noiso.number_of_edges()}")

    # ---- Prepare graph for community detection
    banner("[STEP 15.1] Prepare graph")

    G_comm = nx9_remove_degree0(G_diff_noiso.copy())

    if USE_LARGEST_COMPONENT and G_comm.number_of_nodes() > 0:
        G_comm_lcc = nx9_largest_component(G_comm)
        log(f"[STEP 15.1] Before LCC: {G_comm.number_of_nodes()} nodes")
        log(f"[STEP 15.1] After LCC:  {G_comm_lcc.number_of_nodes()} nodes")
        G_comm = G_comm_lcc

    # Ensure abs_weight
    for u, v, d in G_comm.edges(data=True):
        d["abs_weight"] = abs(float(d.get("abs_weight", d.get("weight", 0.0))))

    log(f"[STEP 15.1] Community graph: nodes={G_comm.number_of_nodes()}, edges={G_comm.number_of_edges()}")

    if G_comm.number_of_nodes() < 5 or G_comm.number_of_edges() < 3:
        log("[SKIP] Graph too small for community detection")
        continue

    nodes_list = sorted(G_comm.nodes())

    # ---- Resolution sweep
    banner("[STEP 15.2] Resolution sweep")

    sweep_rows = []

    for r in COMMUNITY_RESOLUTIONS:
        res_tag = f"res{r:.2f}".replace(".", "p")
        res_dir = ensure_dir(os.path.join(sweep_dir, res_tag))

        log(f"\n[STEP 15.2] Resolution = {r:.2f}")

        labels_list = []
        ncomms_list = []
        modularities = []
        best_cm = None
        best_mod = -999

        for run in range(RUNS_PER_RESOLUTION):
            seed = COMMUNITY_BASE_SEED + run

            cm = detect_leiden(G_comm, seed=seed, resolution=r, weight="abs_weight")
            labels = partition_to_labels(cm, nodes_list)
            labels_list.append(labels)

            n_comms = len(set(cm.values()))
            ncomms_list.append(n_comms)

            # Compute modularity
            from networkx.algorithms.community.quality import modularity as nx_modularity
            comm_sets = {}
            for n, c in cm.items():
                comm_sets.setdefault(c, set()).add(n)
            mod = nx_modularity(G_comm, list(comm_sets.values()), weight="abs_weight")
            modularities.append(mod)

            if mod > best_mod:
                best_mod = mod
                best_cm = cm.copy()

        # ---- Save best partition (raw)
        raw_csv = os.path.join(res_dir, f"{cancer_type}_partition_raw_{res_tag}.csv")
        pd.DataFrame({
            "gene": list(best_cm.keys()),
            "gene_symbol": [symbol_or_self(g, ensg_to_symbol) for g in best_cm.keys()],
            "community": list(best_cm.values())
        }).to_csv(raw_csv, index=False)

        # ---- Merge small communities
        cm_merged, raw_sizes, merged_sizes = merge_small_communities(
            best_cm, k_keep=TARGET_BIG_COMMUNITIES, min_size=MIN_COMMUNITY_SIZE
        )

        merged_csv = os.path.join(res_dir, f"{cancer_type}_partition_merged_{res_tag}.csv")
        pd.DataFrame({
            "gene": list(cm_merged.keys()),
            "gene_symbol": [symbol_or_self(g, ensg_to_symbol) for g in cm_merged.keys()],
            "community": list(cm_merged.values())
        }).to_csv(merged_csv, index=False)

        # ---- Gene lists per community
        comm_to_genes = {}
        for g, c in cm_merged.items():
            comm_to_genes.setdefault(c, []).append(symbol_or_self(g, ensg_to_symbol))

        gene_list_rows = []
        for c, genes in sorted(comm_to_genes.items(), key=lambda x: len(x[1]), reverse=True):
            gene_list_rows.append({
                "community": c,
                "size": len(genes),
                "genes": ";".join(sorted(genes))
            })
        gene_list_csv = os.path.join(res_dir, f"{cancer_type}_gene_lists_{res_tag}.csv")
        pd.DataFrame(gene_list_rows).to_csv(gene_list_csv, index=False)

        # ---- Community size bar plot
        sizes_sorted = sorted(merged_sizes.items(), key=lambda x: x[1], reverse=True)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar([str(k) for k, _ in sizes_sorted], [v for _, v in sizes_sorted],
               color="steelblue", edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Community ID")
        ax.set_ylabel("Number of genes")
        ax.set_title(f"{cancer_type} | Community sizes | res={r:.2f}")
        plt.tight_layout()
        plt.savefig(os.path.join(res_dir, f"{cancer_type}_comm_sizes_{res_tag}.png"), dpi=300)
        plt.close()

        # ---- Network plot colored by community
        pos = nx.spring_layout(G_comm, seed=COMMUNITY_BASE_SEED, weight="abs_weight")

        unique_comms = sorted(set(cm_merged.values()))
        cmap = plt.cm.get_cmap("tab20", max(len(unique_comms), 1))
        comm_colors = {c: cmap(i) for i, c in enumerate(unique_comms)}

        fig, ax = plt.subplots(figsize=(14, 12))
        nx.draw_networkx_edges(G_comm, pos, ax=ax, alpha=0.1, width=0.3, edge_color="gray")

        for comm_id in unique_comms:
            comm_nodes = [n for n in G_comm.nodes() if cm_merged.get(n) == comm_id]
            nx.draw_networkx_nodes(G_comm, pos, nodelist=comm_nodes, ax=ax,
                                   node_size=30, node_color=[comm_colors[comm_id]],
                                   alpha=0.7, label=f"C{comm_id} (n={len(comm_nodes)})")

        # Highlight A3 genes
        a3_in_graph = [g for g in A3_GENES if g in G_comm.nodes()]
        if a3_in_graph:
            nx.draw_networkx_nodes(G_comm, pos, nodelist=a3_in_graph, ax=ax,
                                   node_size=150, node_color="gold", edgecolors="black", linewidths=1.5)
            labels = {n: symbol_or_self(n, ensg_to_symbol) for n in a3_in_graph}
            nx.draw_networkx_labels(G_comm, pos, labels=labels, ax=ax, font_size=8, font_weight="bold")

        ax.set_title(f"{cancer_type} | DIFF communities (merged) | res={r:.2f}")
        ax.legend(fontsize=7, loc="upper right", ncol=2)
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(res_dir, f"{cancer_type}_DIFF_communities_{res_tag}.png"), dpi=300)
        plt.close()

        log(f"  Saved outputs for res={r:.2f}: comms={np.mean(ncomms_list):.1f}, mod={best_mod:.4f}")

        # ---- Stability metrics
        aris = []
        nmis = []
        for a, b in itertools.combinations(range(RUNS_PER_RESOLUTION), 2):
            aris.append(adjusted_rand_score(labels_list[a], labels_list[b]))
            nmis.append(normalized_mutual_info_score(labels_list[a], labels_list[b]))

        sweep_rows.append({
            "resolution": float(r),
            "n_runs": RUNS_PER_RESOLUTION,
            "ncomms_mean": float(np.mean(ncomms_list)),
            "ncomms_std": float(np.std(ncomms_list)),
            "modularity_mean": float(np.mean(modularities)),
            "modularity_std": float(np.std(modularities)),
            "modularity_best": float(best_mod),
            "ARI_mean": float(np.mean(aris)),
            "ARI_std": float(np.std(aris)),
            "NMI_mean": float(np.mean(nmis)),
            "NMI_std": float(np.std(nmis)),
        })

    # ---- Save sweep table
    banner("[STEP 15.8] Save sweep summary")

    sweep_df = pd.DataFrame(sweep_rows).sort_values("resolution")
    sweep_csv = os.path.join(cancer_dir, f"{cancer_type}_resolution_sweep.csv")
    sweep_df.to_csv(sweep_csv, index=False)
    log(f"[SAVE] Sweep table -> {sweep_csv}")

    # ---- Sweep plots (modularity, ARI, NMI vs resolution)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for ax, metric, label in [
        (axes[0], "modularity_mean", "Modularity"),
        (axes[1], "ARI_mean", "ARI (stability)"),
        (axes[2], "NMI_mean", "NMI (stability)"),
    ]:
        std_col = metric.replace("_mean", "_std")
        ax.errorbar(sweep_df["resolution"], sweep_df[metric],
                     yerr=sweep_df.get(std_col, 0), marker="o", capsize=3)
        ax.set_xlabel("Resolution")
        ax.set_ylabel(label)
        ax.set_title(f"{cancer_type} | {label}")

    plt.tight_layout()
    sweep_plot = os.path.join(cancer_dir, f"{cancer_type}_sweep_plots.png")
    plt.savefig(sweep_plot, dpi=300)
    plt.close()
    log(f"[SAVE] Sweep plots -> {sweep_plot}")

    # ---- Identify best resolution (highest modularity with ARI > 0.7)
    stable = sweep_df[sweep_df["ARI_mean"] > 0.7]
    if len(stable) > 0:
        best_row = stable.loc[stable["modularity_best"].idxmax()]
    else:
        best_row = sweep_df.loc[sweep_df["modularity_best"].idxmax()]

    best_res = best_row["resolution"]
    log(f"\n[STEP 15] Best resolution: {best_res:.2f} "
        f"(modularity={best_row['modularity_best']:.4f}, "
        f"ARI={best_row['ARI_mean']:.3f}, "
        f"comms={best_row['ncomms_mean']:.1f})")

    # Copy best partition as the "official" output
    best_res_tag = f"res{best_res:.2f}".replace(".", "p")
    best_merged_csv = os.path.join(sweep_dir, best_res_tag,
                                    f"{cancer_type}_partition_merged_{best_res_tag}.csv")
    if os.path.exists(best_merged_csv):
        import shutil
        final_part = os.path.join(cancer_dir, f"{cancer_type}_best_partition.csv")
        shutil.copy2(best_merged_csv, final_part)
        log(f"[SAVE] Best partition -> {final_part}")

        best_gene_csv = os.path.join(sweep_dir, best_res_tag,
                                      f"{cancer_type}_gene_lists_{best_res_tag}.csv")
        if os.path.exists(best_gene_csv):
            final_genes = os.path.join(cancer_dir, f"{cancer_type}_community_gene_lists.csv")
            shutil.copy2(best_gene_csv, final_genes)
            log(f"[SAVE] Gene lists -> {final_genes}")

    # Save the community graph itself
    comm_pkl = os.path.join(cancer_dir, f"{cancer_type}_G_comm.gpickle")
    with open(comm_pkl, "wb") as f:
        pickle.dump(G_comm, f)
    log(f"[SAVE] Community graph -> {comm_pkl}")

    log(f"\n[STEP 15 COMPLETE for {cancer_type}]")

banner("[STEP 06 COMPLETE — ALL CANCER TYPES]")
