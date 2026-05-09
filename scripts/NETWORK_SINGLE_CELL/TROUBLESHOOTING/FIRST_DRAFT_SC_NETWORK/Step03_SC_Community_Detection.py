#!/usr/bin/env python3
"""
Step03_SC_Community_Detection.py
=================================

Figure 4 — Step 03: Leiden community detection on the single-cell DIFF
network across multiple resolutions.

Mirrors Step05_Community_Detection.py from the Figure 2 (TCGA bulk) pipeline.

IMPORTANT: This script builds the DIFF graph from the saved correlation
matrix pickle at whatever DIFF_THRESHOLD is set in the config. This means
you can change the threshold and re-run Step 03 without re-running the
expensive Step 02 correlation computation.

Workflow:
  1. Load DIFF correlation matrix from Step 02 pickle
  2. Build graph at DIFF_THRESHOLD, remove isolates
  3. Extract largest connected component (LCC)
  4. Run Leiden at multiple resolutions (15 runs each for stability)
  5. Evaluate partition stability (ARI/NMI across runs)
  6. Select best resolution (highest modularity among stable partitions)
  7. Merge communities smaller than MIN_COMMUNITY_SIZE
  8. Save community assignments, gene lists, and annotated graph

Input:
  data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl

Output (→ data/FIG_4/04_communities/):
  sweep/                          — per-resolution assignments + plots
  SC_resolution_sweep.csv         — stability metrics per resolution
  SC_sweep_plots.png              — modularity/ARI/NMI curves
  SC_best_partition.csv           — final gene-to-community assignments
  SC_community_gene_lists.csv     — genes per community
  SC_community_summary.txt        — community sizes and key genes
  SC_G_comm.gpickle               — annotated graph with community labels

Usage:
  conda run -n NETWORK python Step03_SC_Community_Detection.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import pickle
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from datetime import datetime

from network_config_SC import (
    DIR_03_NETWORKS, DIR_04_COMMUNITIES,
    DIFF_THRESHOLD,
    A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    COMMUNITY_METHOD, COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION,
    COMMUNITY_BASE_SEED, USE_LARGEST_COMPONENT,
    TARGET_BIG_COMMUNITIES, MIN_COMMUNITY_SIZE,
    banner, log, ensure_dir
)


# =============================================================================
# GRAPH CONSTRUCTION
# =============================================================================

def build_weighted_graph(corr_df, threshold):
    """Build undirected weighted graph from correlation matrix with |corr| >= threshold."""
    genes = list(corr_df.index)
    G = nx.Graph()
    G.add_nodes_from(genes)
    n = len(genes)
    for i in range(n):
        for j in range(i + 1, n):
            w = float(corr_df.iat[i, j])
            if abs(w) >= threshold and abs(w) < 0.999999:
                G.add_edge(genes[i], genes[j], weight=w, abs_weight=abs(w))
    return G


def remove_isolated(G):
    """Return copy of G with degree-0 nodes removed."""
    G2 = G.copy()
    isolates = [n for n in G2.nodes() if G2.degree(n) == 0]
    G2.remove_nodes_from(isolates)
    return G2


# =============================================================================
# COMMUNITY DETECTION
# =============================================================================

def detect_leiden(G, seed=42, resolution=1.0, weight="abs_weight"):
    """Run Leiden community detection. Returns {node: community_id} dict."""
    import leidenalg
    import igraph as ig

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
        g, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=resolution, seed=seed
    )

    comm_map = {}
    for comm_id, members in enumerate(partition):
        for node_idx in members:
            comm_map[reverse[node_idx]] = comm_id
    return comm_map


def partition_to_labels(comm_map, nodes):
    """Convert {node: community} to label array aligned with node list."""
    return [comm_map.get(n, -1) for n in nodes]


def compute_modularity(G, comm_map, weight="abs_weight"):
    """Compute modularity of a partition."""
    from networkx.algorithms.community.quality import modularity as nx_mod
    comm_sets = {}
    for n, c in comm_map.items():
        comm_sets.setdefault(c, set()).add(n)
    try:
        return nx_mod(G, list(comm_sets.values()), weight=weight)
    except Exception:
        return 0.0


def merge_small_communities(comm_map, G, min_size, target_big, weight="abs_weight"):
    """
    Merge communities smaller than min_size into nearest large community.

    If no communities meet min_size, the largest community is kept as-is
    and all others are merged into it based on edge weight proximity.
    """
    # Count sizes
    comm_sizes = {}
    for n, c in comm_map.items():
        comm_sizes[c] = comm_sizes.get(c, 0) + 1

    # Sort by size descending
    sorted_comms = sorted(comm_sizes.items(), key=lambda x: -x[1])

    # Identify big communities (>= min_size)
    big_comms = [c for c, s in sorted_comms if s >= min_size][:target_big]

    # If no communities are big enough, keep the largest one as the anchor
    if not big_comms:
        if sorted_comms:
            biggest = sorted_comms[0][0]
            log(f"  No communities >= {min_size} genes. Using largest (C{biggest}, "
                f"n={sorted_comms[0][1]}) as anchor.")
            big_comms = [biggest]
        else:
            log("  WARNING: No communities found at all. Returning original partition.")
            return comm_map

    small_nodes = {n: c for n, c in comm_map.items() if c not in big_comms}

    if not small_nodes:
        return comm_map

    log(f"  Merging {len(small_nodes)} nodes from {len(set(small_nodes.values()))} "
        f"small communities into {len(big_comms)} large communities")

    # For each small node, find which big community it connects to most strongly
    merged = dict(comm_map)
    n_reassigned = 0
    for node in small_nodes:
        if node not in G:
            continue
        # Sum edge weights to each big community
        comm_weights = {c: 0.0 for c in big_comms}
        for neighbor in G.neighbors(node):
            nb_comm = comm_map.get(neighbor)
            if nb_comm in big_comms:
                w = abs(float(G[node][neighbor].get(weight, 1.0)))
                comm_weights[nb_comm] += w

        # Assign to strongest connection, or to largest big community if no connections
        best_weight = max(comm_weights.values())
        if best_weight > 0:
            best = max(comm_weights, key=comm_weights.get)
        else:
            # No edges to any big community — assign to the largest one
            best = big_comms[0]

        merged[node] = best
        n_reassigned += 1

    log(f"  Reassigned {n_reassigned} nodes")
    return merged


def remove_intra_community_isolates(comm_map, G):
    """
    Remove nodes that have zero edges to other members of their own community.
    These are 'orphan' nodes created by merging — they were assigned to a community
    they have no connections to. Matches the TCGA pipeline's nx9_remove_degree0
    behavior applied after community assignment.
    """
    cleaned = {}
    n_removed = 0
    for node, comm in comm_map.items():
        if node not in G:
            n_removed += 1
            continue
        # Check if node has at least one intra-community neighbor
        has_intra = False
        for neighbor in G.neighbors(node):
            if comm_map.get(neighbor) == comm:
                has_intra = True
                break
        if has_intra:
            cleaned[node] = comm
        else:
            n_removed += 1

    if n_removed > 0:
        log(f"  Removed {n_removed} intra-community isolates (no edges within their community)")
    return cleaned


def renumber_communities(comm_map):
    """Renumber communities 0, 1, 2, ... by descending size."""
    comm_sizes = {}
    for c in comm_map.values():
        comm_sizes[c] = comm_sizes.get(c, 0) + 1
    sorted_comms = sorted(comm_sizes, key=lambda c: -comm_sizes[c])
    remap = {old: new for new, old in enumerate(sorted_comms)}
    return {n: remap[c] for n, c in comm_map.items()}


# =============================================================================
# MAIN
# =============================================================================

def main():
    start_time = datetime.now()
    banner("[STEP 03] Single-Cell Community Detection (Leiden)")
    log(f"Start time: {start_time}")

    ensure_dir(DIR_04_COMMUNITIES)
    sweep_dir = ensure_dir(os.path.join(DIR_04_COMMUNITIES, "sweep"))

    # =========================================================================
    # 1. Load DIFF correlation matrix and build graph at config threshold
    # =========================================================================
    banner("[STEP 03.1] Load DIFF correlation matrix and build graph")

    corr_path = os.path.join(DIR_03_NETWORKS, "corr_matrices", "SC_corr_DIFF.pkl")
    log(f"Loading: {corr_path}")
    with open(corr_path, "rb") as f:
        corr_diff = pickle.load(f)
    log(f"  DIFF correlation matrix: {corr_diff.shape}")

    log(f"  Building graph at |delta-rho| >= {DIFF_THRESHOLD}")
    G_diff = build_weighted_graph(corr_diff, DIFF_THRESHOLD)
    G_diff_noiso = remove_isolated(G_diff)

    n_total = G_diff.number_of_nodes()
    n_iso = n_total - G_diff_noiso.number_of_nodes()
    log(f"  Total genes in matrix: {n_total}")
    log(f"  Non-isolated nodes: {G_diff_noiso.number_of_nodes()}")
    log(f"  Edges: {G_diff_noiso.number_of_edges()}")
    log(f"  Isolated nodes removed: {n_iso}")

    # Also save the rebuilt graph for Steps 04 and 05
    graph_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "graph_objects"))
    graph_save_path = os.path.join(graph_dir, "SC_G_diff_noiso.gpickle")
    with open(graph_save_path, "wb") as f:
        pickle.dump(G_diff_noiso, f)
    log(f"  [SAVE] Updated DIFF graph → {graph_save_path}")

    if G_diff_noiso.number_of_nodes() < 5 or G_diff_noiso.number_of_edges() < 3:
        log(f"WARNING: Graph too small for community detection "
            f"({G_diff_noiso.number_of_nodes()} nodes, {G_diff_noiso.number_of_edges()} edges)")
        log(f"  Consider further relaxing DIFF_THRESHOLD (currently {DIFF_THRESHOLD})")
        return

    # =========================================================================
    # 2. Extract LCC
    # =========================================================================
    if USE_LARGEST_COMPONENT and G_diff_noiso.number_of_nodes() > 0:
        banner("[STEP 03.2] Extract largest connected component")
        components = list(nx.connected_components(G_diff_noiso))
        lcc = max(components, key=len)
        G_comm = G_diff_noiso.subgraph(lcc).copy()
        log(f"  Components: {len(components)}")
        log(f"  LCC: {len(lcc)} nodes ({100*len(lcc)/G_diff_noiso.number_of_nodes():.1f}%)")
        log(f"  Non-LCC nodes dropped: {G_diff_noiso.number_of_nodes() - len(lcc)}")
    else:
        G_comm = G_diff_noiso.copy()

    if G_comm.number_of_nodes() < 5 or G_comm.number_of_edges() < 3:
        log("WARNING: LCC too small for community detection. Exiting.")
        return

    # Ensure abs_weight on all edges
    for u, v, d in G_comm.edges(data=True):
        d["abs_weight"] = abs(float(d.get("abs_weight", d.get("weight", 1.0))))

    nodes = sorted(G_comm.nodes())

    # =========================================================================
    # 3. Resolution sweep with stability assessment
    # =========================================================================
    banner("[STEP 03.3] Resolution sweep")

    sweep_results = []

    for res in COMMUNITY_RESOLUTIONS:
        log(f"\n  Resolution {res}:")

        partitions = []
        modularities = []
        for r in range(RUNS_PER_RESOLUTION):
            seed = COMMUNITY_BASE_SEED + r
            cm = detect_leiden(G_comm, seed=seed, resolution=res)
            partitions.append(cm)
            mod = compute_modularity(G_comm, cm)
            modularities.append(mod)

        labels_list = [partition_to_labels(p, nodes) for p in partitions]

        ari_pairs = []
        nmi_pairs = []
        for i, j in itertools.combinations(range(len(labels_list)), 2):
            ari_pairs.append(adjusted_rand_score(labels_list[i], labels_list[j]))
            nmi_pairs.append(normalized_mutual_info_score(labels_list[i], labels_list[j]))

        mean_mod = np.mean(modularities)
        mean_ari = np.mean(ari_pairs)
        mean_nmi = np.mean(nmi_pairs)
        n_comms = len(set(partitions[0].values()))

        log(f"    Communities: {n_comms}")
        log(f"    Modularity: {mean_mod:.4f} (± {np.std(modularities):.4f})")
        log(f"    ARI stability: {mean_ari:.4f}")
        log(f"    NMI stability: {mean_nmi:.4f}")

        sweep_results.append({
            "resolution": res,
            "n_communities": n_comms,
            "modularity_mean": mean_mod,
            "modularity_std": np.std(modularities),
            "ari_mean": mean_ari,
            "ari_std": np.std(ari_pairs),
            "nmi_mean": mean_nmi,
            "nmi_std": np.std(nmi_pairs),
        })

        # Save per-resolution partition (first run)
        res_df = pd.DataFrame([
            {"gene": g, "community": c} for g, c in partitions[0].items()
        ])
        res_path = os.path.join(sweep_dir, f"SC_partition_res{res:.1f}.csv")
        res_df.to_csv(res_path, index=False)

    sweep_df = pd.DataFrame(sweep_results)
    sweep_csv = os.path.join(DIR_04_COMMUNITIES, "SC_resolution_sweep.csv")
    sweep_df.to_csv(sweep_csv, index=False)
    log(f"\n  [SAVE] Resolution sweep → {sweep_csv}")

    # =========================================================================
    # 4. Select best resolution
    # =========================================================================
    banner("[STEP 03.4] Select best resolution")

    # Strategy: highest modularity among reasonably stable partitions (ARI > 0.3)
    stable = sweep_df[sweep_df["ari_mean"] > 0.3]
    if len(stable) == 0:
        log("  WARNING: No stable partitions found. Using highest modularity.")
        stable = sweep_df

    best_idx = stable["modularity_mean"].idxmax()
    best_res = float(stable.loc[best_idx, "resolution"])
    best_mod = float(stable.loc[best_idx, "modularity_mean"])
    best_ari = float(stable.loc[best_idx, "ari_mean"])
    best_n_comms = int(stable.loc[best_idx, "n_communities"])

    log(f"  Best resolution: {best_res}")
    log(f"  Modularity: {best_mod:.4f}")
    log(f"  ARI stability: {best_ari:.4f}")
    log(f"  Communities: {best_n_comms}")

    # Run final partition at best resolution with base seed
    best_partition = detect_leiden(G_comm, seed=COMMUNITY_BASE_SEED, resolution=best_res)

    # =========================================================================
    # 5. Merge small communities
    # =========================================================================
    banner("[STEP 03.5] Merge small communities")

    comm_sizes = {}
    for c in best_partition.values():
        comm_sizes[c] = comm_sizes.get(c, 0) + 1

    log(f"  Before merging: {len(comm_sizes)} communities")
    for c, s in sorted(comm_sizes.items(), key=lambda x: -x[1]):
        log(f"    Community {c}: {s} genes")

    merged = merge_small_communities(
        best_partition, G_comm, MIN_COMMUNITY_SIZE, TARGET_BIG_COMMUNITIES
    )

    # Remove nodes with zero intra-community edges (orphans from merging)
    cleaned = remove_intra_community_isolates(merged, G_comm)

    # Remove small disconnected islands within each community
    # After merging, some communities contain multiple disconnected sub-components
    # (e.g., 2-11 node islands not connected to the main body). Keep only the
    # largest connected component within each community.
    pruned = {}
    n_island_removed = 0
    for comm_id in sorted(set(cleaned.values())):
        comm_nodes = [n for n, c in cleaned.items() if c == comm_id]
        if len(comm_nodes) <= 1:
            for n in comm_nodes:
                pruned[n] = comm_id
            continue
        sub = G_comm.subgraph(comm_nodes)
        components = list(nx.connected_components(sub))
        if len(components) <= 1:
            for n in comm_nodes:
                pruned[n] = comm_id
        else:
            # Keep only the largest connected component
            lcc = max(components, key=len)
            island_nodes = set(comm_nodes) - lcc
            n_island_removed += len(island_nodes)
            if island_nodes:
                sizes = sorted([len(c) for c in components], reverse=True)
                log(f"  Community {comm_id}: {len(components)} components "
                    f"(sizes: {sizes}), keeping LCC ({len(lcc)}), "
                    f"dropping {len(island_nodes)} island nodes")
            for n in lcc:
                pruned[n] = comm_id

    if n_island_removed > 0:
        log(f"  Removed {n_island_removed} island nodes across all communities")

    final_partition = renumber_communities(pruned)

    # Update graph to match final partition (remove nodes not in partition)
    nodes_in_partition = set(final_partition.keys())
    nodes_to_remove = [n for n in G_comm.nodes() if n not in nodes_in_partition]
    if nodes_to_remove:
        G_comm = G_comm.copy()
        G_comm.remove_nodes_from(nodes_to_remove)
        log(f"  Removed {len(nodes_to_remove)} nodes from graph (not in final partition)")

    final_sizes = {}
    for c in final_partition.values():
        final_sizes[c] = final_sizes.get(c, 0) + 1

    log(f"\n  After merging: {len(final_sizes)} communities")
    total_genes = 0
    for c, s in sorted(final_sizes.items()):
        log(f"    Community {c}: {s} genes")
        total_genes += s
    log(f"  Total genes in communities: {total_genes}")

    # =========================================================================
    # 6. Save results
    # =========================================================================
    banner("[STEP 03.6] Save community assignments")

    # Best partition CSV
    part_rows = []
    for gene, comm in sorted(final_partition.items(), key=lambda x: (x[1], x[0])):
        alias = A3_ID_TO_ALIAS.get(gene, "")
        part_rows.append({"gene": gene, "community": comm, "a3_alias": alias})

    part_df = pd.DataFrame(part_rows)
    part_path = os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv")
    part_df.to_csv(part_path, index=False)
    log(f"  [SAVE] Best partition ({len(part_df)} genes) → {part_path}")

    # Community gene lists
    gene_list_rows = []
    for comm in sorted(final_sizes.keys()):
        genes = sorted([g for g, c in final_partition.items() if c == comm])
        gene_list_rows.append({
            "community": comm,
            "size": len(genes),
            "genes": ";".join(genes)
        })

    gene_list_df = pd.DataFrame(gene_list_rows)
    gene_list_path = os.path.join(DIR_04_COMMUNITIES, "SC_community_gene_lists.csv")
    gene_list_df.to_csv(gene_list_path, index=False)
    log(f"  [SAVE] Gene lists → {gene_list_path}")

    # Annotated graph
    for n in G_comm.nodes():
        G_comm.nodes[n]["community"] = final_partition.get(n, -1)

    graph_path = os.path.join(DIR_04_COMMUNITIES, "SC_G_comm.gpickle")
    with open(graph_path, "wb") as f:
        pickle.dump(G_comm, f)
    log(f"  [SAVE] Annotated graph → {graph_path}")

    # =========================================================================
    # 7. Community summary with A3 and biomarker gene locations
    # =========================================================================
    banner("[STEP 03.7] Community summary")

    summary_lines = []
    summary_lines.append(f"Single-Cell Community Detection Summary")
    summary_lines.append(f"=" * 60)
    summary_lines.append(f"DIFF threshold: |delta-rho| >= {DIFF_THRESHOLD}")
    summary_lines.append(f"Resolution: {best_res}")
    summary_lines.append(f"Modularity: {best_mod:.4f}")
    summary_lines.append(f"Total genes in communities: {total_genes}")
    summary_lines.append(f"Communities: {len(final_sizes)}")
    summary_lines.append("")

    for comm in sorted(final_sizes.keys()):
        genes_in_comm = [g for g, c in final_partition.items() if c == comm]
        summary_lines.append(f"Community {comm} ({len(genes_in_comm)} genes):")

        # A3 genes
        for a3 in A3_GENES:
            if a3 in genes_in_comm:
                alias = A3_ID_TO_ALIAS.get(a3, a3)
                summary_lines.append(f"  ** {alias} ({a3}) **")

        # Biomarkers
        for bm in BIOMARKERS:
            if bm in genes_in_comm:
                summary_lines.append(f"  [biomarker] {bm}")

        summary_lines.append("")

    summary_path = os.path.join(DIR_04_COMMUNITIES, "SC_community_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines))
    log(f"  [SAVE] Summary → {summary_path}")

    # Print summary to stdout
    for line in summary_lines:
        print(line)

    # =========================================================================
    # 8. Sweep plots
    # =========================================================================
    banner("[STEP 03.8] Resolution sweep plots")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    axes[0].errorbar(sweep_df["resolution"], sweep_df["modularity_mean"],
                     yerr=sweep_df["modularity_std"], fmt="o-", c="steelblue")
    axes[0].axvline(best_res, ls="--", c="red", alpha=0.5)
    axes[0].set_xlabel("Resolution")
    axes[0].set_ylabel("Modularity")
    axes[0].set_title("Modularity vs Resolution")

    axes[1].errorbar(sweep_df["resolution"], sweep_df["ari_mean"],
                     yerr=sweep_df["ari_std"], fmt="o-", c="darkgreen")
    axes[1].axvline(best_res, ls="--", c="red", alpha=0.5)
    axes[1].set_xlabel("Resolution")
    axes[1].set_ylabel("ARI (stability)")
    axes[1].set_title("Partition Stability (ARI)")

    axes[2].plot(sweep_df["resolution"], sweep_df["n_communities"], "o-", c="purple")
    axes[2].axvline(best_res, ls="--", c="red", alpha=0.5)
    axes[2].set_xlabel("Resolution")
    axes[2].set_ylabel("N Communities")
    axes[2].set_title("Community Count")

    plt.tight_layout()
    plot_path = os.path.join(DIR_04_COMMUNITIES, "SC_sweep_plots.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    log(f"  [SAVE] Sweep plots → {plot_path}")

    # =========================================================================
    # DONE
    # =========================================================================
    elapsed = datetime.now() - start_time
    banner(f"[STEP 03 COMPLETE] {len(final_sizes)} communities, "
           f"{total_genes} genes | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
