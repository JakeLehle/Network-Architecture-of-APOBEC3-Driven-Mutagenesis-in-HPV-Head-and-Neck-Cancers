#!/usr/bin/env python3
"""
Step03_SC_Community_Detection.py
=================================

Figure 4 -- Step 03: Leiden community detection on the single-cell DIFF
network across multiple resolutions.

IMPORTANT: This script builds the DIFF graph from the saved correlation
matrix pickle. The DIFF threshold is auto-selected from Step02.1 sweep
data (component peak + LCC < 300 criterion), or falls back to the config
value if auto-selection is disabled or sweep data is unavailable.

The Leiden resolution is auto-selected via modularity jump detection:
finds the resolution with the largest single-step Q increase, subject
to an ARI stability floor.

Input:
  data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl
  data/FIG_4/03_correlation_networks/sweep_detailed.csv  (from Step02.1)

Output (to data/FIG_4/04_communities/):
  SC_resolution_sweep.csv, SC_best_partition.csv, SC_community_gene_lists.csv,
  SC_community_summary.txt, SC_G_comm.gpickle, SC_selected_parameters.txt,
  SC_sweep_plots.png, sweep/

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
    DIFF_THRESHOLD, DIFF_THRESHOLD_AUTO, DIFF_THRESHOLD_MAX_LCC,
    RESOLUTION_AUTO, RESOLUTION_MIN_ARI,
    A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    COMMUNITY_METHOD, COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION,
    COMMUNITY_BASE_SEED, USE_LARGEST_COMPONENT,
    TARGET_BIG_COMMUNITIES, MIN_COMMUNITY_SIZE,
    banner, log, ensure_dir
)


# =============================================================================
# AUTO-SELECTION: DIFF THRESHOLD
# =============================================================================

def auto_select_threshold(sweep_csv, max_lcc=300, comp_fraction=0.80):
    """
    Select optimal DIFF threshold from sweep data.
    Criterion: peak connected components with LCC < max_lcc.
    Returns selected threshold (float) or None if selection fails.
    """
    if not os.path.exists(sweep_csv):
        log(f"  Sweep data not found: {sweep_csv}")
        return None

    df = pd.read_csv(sweep_csv)

    col_map = {}
    for col in df.columns:
        cl = col.strip().lower()
        if 'threshold' in cl:
            col_map[col] = 'threshold'
        elif cl in ['comp', 'components', 'n_components']:
            col_map[col] = 'components'
        elif cl in ['lcc', 'lcc_size']:
            col_map[col] = 'lcc'
        elif 'node' in cl:
            col_map[col] = 'nodes'
        elif 'edge' in cl:
            col_map[col] = 'edges'
        elif 'avgdeg' in cl or 'avg_deg' in cl or 'avg deg' in cl:
            col_map[col] = 'avg_deg'
    df = df.rename(columns=col_map)

    for col in ['threshold', 'components', 'lcc', 'nodes', 'edges']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.sort_values('threshold').reset_index(drop=True)

    peak_idx = df['components'].idxmax()
    peak_row = df.loc[peak_idx]
    peak_comp = peak_row['components']
    peak_thresh = peak_row['threshold']
    peak_lcc = peak_row['lcc']

    log(f"  Component peak: threshold={peak_thresh}, "
        f"components={int(peak_comp)}, LCC={int(peak_lcc)}")

    if peak_lcc <= max_lcc:
        selected = peak_row
        reason = "Component peak with LCC within limit"
    else:
        candidates = df[(df['threshold'] > peak_thresh) & (df['lcc'] <= max_lcc)]
        if len(candidates) > 0:
            selected = candidates.iloc[0]
            reason = f"LCC at peak ({int(peak_lcc)}) > {max_lcc}, stepped up"
        else:
            above = df[df['threshold'] > peak_thresh]
            if len(above) > 0:
                selected = above.loc[above['lcc'].idxmin()]
                reason = f"No threshold gives LCC < {max_lcc}, using smallest LCC"
            else:
                return None

    sel_comp = selected['components']
    if sel_comp < peak_comp * comp_fraction:
        log(f"  WARNING: Selected has {int(sel_comp)} components "
            f"({sel_comp/peak_comp:.0%} of peak {int(peak_comp)})")

    sel_thresh = float(selected['threshold'])
    log(f"  Auto-selected threshold: {sel_thresh}")
    log(f"    Reason: {reason}")
    log(f"    Nodes: {int(selected['nodes']):,}, Edges: {int(selected['edges']):,}")
    log(f"    Components: {int(sel_comp)}, LCC: {int(selected['lcc'])}")

    return sel_thresh


# =============================================================================
# AUTO-SELECTION: LEIDEN RESOLUTION
# =============================================================================

def auto_select_resolution(sweep_df, min_ari=0.65):
    """
    Select optimal Leiden resolution via modularity jump detection.
    Finds the largest single-step Q increase, subject to ARI floor.
    Returns selected resolution (float).
    """
    df = sweep_df.sort_values('resolution').reset_index(drop=True)
    df['delta_q'] = df['modularity_mean'].diff().fillna(0)

    log(f"\n  Modularity jump analysis:")
    log(f"  {'Res':>5s} {'Comms':>6s} {'Q':>7s} {'dQ':>7s} "
        f"{'ARI':>7s} {'NMI':>7s}")
    log(f"  {'-'*5} {'-'*6} {'-'*7} {'-'*7} {'-'*7} {'-'*7}")
    for _, row in df.iterrows():
        log(f"  {row['resolution']:>5.1f} {int(row['n_communities']):>6} "
            f"{row['modularity_mean']:>7.4f} {row['delta_q']:>+7.4f} "
            f"{row['ari_mean']:>7.4f} {row['nmi_mean']:>7.4f}")

    if len(df) < 2:
        log(f"  Only one resolution tested, using it")
        return float(df.iloc[0]['resolution'])

    jump_candidates = df.iloc[1:].copy()
    best_jump_idx = jump_candidates['delta_q'].idxmax()
    best_jump = df.loc[best_jump_idx]
    best_dq = best_jump['delta_q']

    log(f"\n  Largest modularity jump: dQ={best_dq:+.4f} "
        f"at resolution={best_jump['resolution']}")

    if best_jump['ari_mean'] >= min_ari:
        selected = best_jump
        reason = f"Largest Q jump (dQ={best_dq:+.4f}), ARI passes floor"
    else:
        prev_idx = best_jump_idx - 1
        if prev_idx >= 0:
            prev = df.loc[prev_idx]
            selected = prev
            reason = (f"Jump at r={best_jump['resolution']} has "
                      f"ARI={best_jump['ari_mean']:.3f} < {min_ari}, "
                      f"stepped back to r={prev['resolution']}")
            log(f"  ARI at jump ({best_jump['ari_mean']:.4f}) < floor ({min_ari})")
        else:
            selected = best_jump
            reason = "Largest jump (no earlier resolution available)"

    next_idx = best_jump_idx + 1
    if next_idx < len(df):
        beyond = df.loc[next_idx]
        marginal_q = beyond['modularity_mean'] - selected['modularity_mean']
        marginal_ari = selected['ari_mean'] - beyond['ari_mean']
        if marginal_ari > 0:
            efficiency = marginal_q / marginal_ari
            verdict = "worth exploring" if efficiency > 0.5 else "diminishing returns"
            log(f"  Beyond (r={beyond['resolution']}): "
                f"dQ={marginal_q:+.4f}, ARI cost={marginal_ari:+.4f}, "
                f"efficiency={efficiency:.2f}x ({verdict})")

    sel_res = float(selected['resolution'])
    log(f"\n  Auto-selected resolution: {sel_res}")
    log(f"    Reason: {reason}")
    log(f"    Communities: {int(selected['n_communities'])}")
    log(f"    Modularity: {selected['modularity_mean']:.4f}")
    log(f"    ARI: {selected['ari_mean']:.4f}")

    return sel_res


# =============================================================================
# GRAPH CONSTRUCTION
# =============================================================================

def build_weighted_graph(corr_df, threshold):
    """Build undirected weighted graph from correlation matrix."""
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
    G2 = G.copy()
    isolates = [n for n in G2.nodes() if G2.degree(n) == 0]
    G2.remove_nodes_from(isolates)
    return G2


# =============================================================================
# COMMUNITY DETECTION
# =============================================================================

def detect_leiden(G, seed=42, resolution=1.0, weight="abs_weight"):
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
    return [comm_map.get(n, -1) for n in nodes]


def compute_modularity(G, comm_map, weight="abs_weight"):
    from networkx.algorithms.community.quality import modularity as nx_mod
    comm_sets = {}
    for n, c in comm_map.items():
        comm_sets.setdefault(c, set()).add(n)
    try:
        return nx_mod(G, list(comm_sets.values()), weight=weight)
    except Exception:
        return 0.0


def merge_small_communities(comm_map, G, min_size, target_big, weight="abs_weight"):
    comm_sizes = {}
    for n, c in comm_map.items():
        comm_sizes[c] = comm_sizes.get(c, 0) + 1

    sorted_comms = sorted(comm_sizes.items(), key=lambda x: -x[1])
    big_comms = [c for c, s in sorted_comms if s >= min_size][:target_big]

    if not big_comms:
        if sorted_comms:
            biggest = sorted_comms[0][0]
            log(f"  No communities >= {min_size} genes. Using largest (C{biggest}, "
                f"n={sorted_comms[0][1]}) as anchor.")
            big_comms = [biggest]
        else:
            return comm_map

    small_nodes = {n: c for n, c in comm_map.items() if c not in big_comms}
    if not small_nodes:
        return comm_map

    log(f"  Merging {len(small_nodes)} nodes from {len(set(small_nodes.values()))} "
        f"small communities into {len(big_comms)} large communities")

    merged = dict(comm_map)
    n_reassigned = 0
    for node in small_nodes:
        if node not in G:
            continue
        comm_weights = {c: 0.0 for c in big_comms}
        for neighbor in G.neighbors(node):
            nb_comm = comm_map.get(neighbor)
            if nb_comm in big_comms:
                w = abs(float(G[node][neighbor].get(weight, 1.0)))
                comm_weights[nb_comm] += w

        best_weight = max(comm_weights.values())
        best = max(comm_weights, key=comm_weights.get) if best_weight > 0 else big_comms[0]
        merged[node] = best
        n_reassigned += 1

    log(f"  Reassigned {n_reassigned} nodes")
    return merged


def remove_intra_community_isolates(comm_map, G):
    cleaned = {}
    n_removed = 0
    for node, comm in comm_map.items():
        if node not in G:
            n_removed += 1
            continue
        has_intra = any(comm_map.get(nb) == comm for nb in G.neighbors(node))
        if has_intra:
            cleaned[node] = comm
        else:
            n_removed += 1
    if n_removed > 0:
        log(f"  Removed {n_removed} intra-community isolates")
    return cleaned


def renumber_communities(comm_map):
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
    # 1. Auto-select DIFF threshold
    # =========================================================================
    banner("[STEP 03.1] Select DIFF threshold")

    active_threshold = DIFF_THRESHOLD

    if DIFF_THRESHOLD_AUTO:
        sweep_csv = os.path.join(DIR_03_NETWORKS, "sweep_detailed.csv")
        log(f"  Auto-selection enabled (max LCC = {DIFF_THRESHOLD_MAX_LCC})")
        auto_thresh = auto_select_threshold(sweep_csv, max_lcc=DIFF_THRESHOLD_MAX_LCC)
        if auto_thresh is not None:
            active_threshold = auto_thresh
            log(f"  Using auto-selected threshold: {active_threshold}")
        else:
            log(f"  Auto-selection failed, using config fallback: {active_threshold}")
    else:
        log(f"  Auto-selection disabled, using config value: {active_threshold}")

    # =========================================================================
    # 2. Load DIFF correlation matrix and build graph
    # =========================================================================
    banner("[STEP 03.2] Load DIFF correlation matrix and build graph")

    corr_path = os.path.join(DIR_03_NETWORKS, "corr_matrices", "SC_corr_DIFF.pkl")
    log(f"Loading: {corr_path}")
    with open(corr_path, "rb") as f:
        corr_diff = pickle.load(f)
    log(f"  DIFF correlation matrix: {corr_diff.shape}")

    log(f"  Building graph at |delta-rho| >= {active_threshold}")
    G_diff = build_weighted_graph(corr_diff, active_threshold)
    G_diff_noiso = remove_isolated(G_diff)

    n_total = G_diff.number_of_nodes()
    n_iso = n_total - G_diff_noiso.number_of_nodes()
    log(f"  Total genes in matrix: {n_total}")
    log(f"  Non-isolated nodes: {G_diff_noiso.number_of_nodes()}")
    log(f"  Edges: {G_diff_noiso.number_of_edges()}")
    log(f"  Isolated nodes removed: {n_iso}")

    graph_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "graph_objects"))
    graph_save_path = os.path.join(graph_dir, "SC_G_diff_noiso.gpickle")
    with open(graph_save_path, "wb") as f:
        pickle.dump(G_diff_noiso, f)
    log(f"  [SAVE] Updated DIFF graph -> {graph_save_path}")

    if G_diff_noiso.number_of_nodes() < 5 or G_diff_noiso.number_of_edges() < 3:
        log(f"WARNING: Graph too small ({G_diff_noiso.number_of_nodes()} nodes)")
        return

    # =========================================================================
    # 3. Extract LCC
    # =========================================================================
    if USE_LARGEST_COMPONENT and G_diff_noiso.number_of_nodes() > 0:
        banner("[STEP 03.3] Extract largest connected component")
        components = list(nx.connected_components(G_diff_noiso))
        lcc = max(components, key=len)
        G_comm = G_diff_noiso.subgraph(lcc).copy()
        log(f"  Components: {len(components)}")
        log(f"  LCC: {len(lcc)} nodes "
            f"({100*len(lcc)/G_diff_noiso.number_of_nodes():.1f}%)")
        log(f"  Non-LCC nodes dropped: {G_diff_noiso.number_of_nodes() - len(lcc)}")
    else:
        G_comm = G_diff_noiso.copy()

    if G_comm.number_of_nodes() < 5:
        log("WARNING: LCC too small. Exiting.")
        return

    for u, v, d in G_comm.edges(data=True):
        d["abs_weight"] = abs(float(d.get("abs_weight", d.get("weight", 1.0))))

    nodes = sorted(G_comm.nodes())

    # =========================================================================
    # 4. Resolution sweep
    # =========================================================================
    banner("[STEP 03.4] Resolution sweep")

    sweep_results = []
    all_partitions = {}

    for res in COMMUNITY_RESOLUTIONS:
        log(f"\n  Resolution {res}:")

        partitions = []
        modularities = []
        for r in range(RUNS_PER_RESOLUTION):
            seed = COMMUNITY_BASE_SEED + r
            cm = detect_leiden(G_comm, seed=seed, resolution=res)
            partitions.append(cm)
            modularities.append(compute_modularity(G_comm, cm))

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
        log(f"    Modularity: {mean_mod:.4f} (+/- {np.std(modularities):.4f})")
        log(f"    ARI stability: {mean_ari:.4f}")
        log(f"    NMI stability: {mean_nmi:.4f}")

        sweep_results.append({
            "resolution": res, "n_communities": n_comms,
            "modularity_mean": mean_mod, "modularity_std": np.std(modularities),
            "ari_mean": mean_ari, "ari_std": np.std(ari_pairs),
            "nmi_mean": mean_nmi, "nmi_std": np.std(nmi_pairs),
        })
        all_partitions[res] = partitions[0]

        res_df = pd.DataFrame([{"gene": g, "community": c} for g, c in partitions[0].items()])
        res_df.to_csv(os.path.join(sweep_dir, f"SC_partition_res{res:.1f}.csv"), index=False)

    sweep_df = pd.DataFrame(sweep_results)
    sweep_csv_out = os.path.join(DIR_04_COMMUNITIES, "SC_resolution_sweep.csv")
    sweep_df.to_csv(sweep_csv_out, index=False)
    log(f"\n  [SAVE] Resolution sweep -> {sweep_csv_out}")

    # =========================================================================
    # 5. Select best resolution
    # =========================================================================
    banner("[STEP 03.5] Select best resolution")

    if RESOLUTION_AUTO:
        log(f"  Auto-selection enabled (ARI floor = {RESOLUTION_MIN_ARI})")
        best_res = auto_select_resolution(sweep_df, min_ari=RESOLUTION_MIN_ARI)
    else:
        log(f"  Auto-selection disabled, using max modularity with ARI > 0.3")
        stable = sweep_df[sweep_df["ari_mean"] > 0.3]
        if len(stable) == 0:
            stable = sweep_df
        best_idx = stable["modularity_mean"].idxmax()
        best_res = float(stable.loc[best_idx, "resolution"])

    best_row = sweep_df[sweep_df["resolution"] == best_res].iloc[0]
    best_mod = float(best_row["modularity_mean"])
    best_ari = float(best_row["ari_mean"])
    best_nmi = float(best_row["nmi_mean"])
    best_n_comms = int(best_row["n_communities"])

    log(f"\n  Final selection: resolution={best_res}, Q={best_mod:.4f}, "
        f"ARI={best_ari:.4f}, communities={best_n_comms}")

    best_partition = detect_leiden(G_comm, seed=COMMUNITY_BASE_SEED, resolution=best_res)

    # =========================================================================
    # 6. Merge small communities
    # =========================================================================
    banner("[STEP 03.6] Merge small communities")

    comm_sizes = {}
    for c in best_partition.values():
        comm_sizes[c] = comm_sizes.get(c, 0) + 1

    log(f"  Before merging: {len(comm_sizes)} communities")
    for c, s in sorted(comm_sizes.items(), key=lambda x: -x[1]):
        log(f"    Community {c}: {s} genes")

    merged = merge_small_communities(
        best_partition, G_comm, MIN_COMMUNITY_SIZE, TARGET_BIG_COMMUNITIES)
    cleaned = remove_intra_community_isolates(merged, G_comm)

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
            lcc_sub = max(components, key=len)
            island_nodes = set(comm_nodes) - lcc_sub
            n_island_removed += len(island_nodes)
            if island_nodes:
                sizes = sorted([len(c) for c in components], reverse=True)
                log(f"  Community {comm_id}: {len(components)} components "
                    f"(sizes: {sizes}), keeping LCC ({len(lcc_sub)}), "
                    f"dropping {len(island_nodes)} island nodes")
            for n in lcc_sub:
                pruned[n] = comm_id

    if n_island_removed > 0:
        log(f"  Removed {n_island_removed} island nodes across all communities")

    final_partition = renumber_communities(pruned)

    nodes_in_partition = set(final_partition.keys())
    nodes_to_remove = [n for n in G_comm.nodes() if n not in nodes_in_partition]
    if nodes_to_remove:
        G_comm = G_comm.copy()
        G_comm.remove_nodes_from(nodes_to_remove)
        log(f"  Removed {len(nodes_to_remove)} nodes from graph")

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
    # 7. Save results
    # =========================================================================
    banner("[STEP 03.7] Save community assignments")

    part_rows = [{"gene": g, "community": c, "a3_alias": A3_ID_TO_ALIAS.get(g, "")}
                 for g, c in sorted(final_partition.items(), key=lambda x: (x[1], x[0]))]
    part_df = pd.DataFrame(part_rows)
    part_df.to_csv(os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv"), index=False)
    log(f"  [SAVE] Best partition ({len(part_df)} genes)")

    gene_list_rows = []
    for comm in sorted(final_sizes.keys()):
        genes = sorted([g for g, c in final_partition.items() if c == comm])
        gene_list_rows.append({"community": comm, "size": len(genes), "genes": ";".join(genes)})
    pd.DataFrame(gene_list_rows).to_csv(
        os.path.join(DIR_04_COMMUNITIES, "SC_community_gene_lists.csv"), index=False)
    log(f"  [SAVE] Gene lists")

    for n in G_comm.nodes():
        G_comm.nodes[n]["community"] = final_partition.get(n, -1)
    with open(os.path.join(DIR_04_COMMUNITIES, "SC_G_comm.gpickle"), "wb") as f:
        pickle.dump(G_comm, f)
    log(f"  [SAVE] Annotated graph")

    params_path = os.path.join(DIR_04_COMMUNITIES, "SC_selected_parameters.txt")
    with open(params_path, 'w') as f:
        f.write(f"DIFF_THRESHOLD={active_threshold}\n")
        f.write(f"LEIDEN_RESOLUTION={best_res}\n")
        f.write(f"N_COMMUNITIES={len(final_sizes)}\n")
        f.write(f"N_GENES={total_genes}\n")
        f.write(f"MODULARITY={best_mod:.4f}\n")
        f.write(f"ARI={best_ari:.4f}\n")
        f.write(f"NMI={best_nmi:.4f}\n")
    log(f"  [SAVE] Selected parameters -> {params_path}")

    # =========================================================================
    # 8. Community summary
    # =========================================================================
    banner("[STEP 03.8] Community summary")

    summary_lines = []
    summary_lines.append("Single-Cell Community Detection Summary")
    summary_lines.append("=" * 60)
    summary_lines.append(f"DIFF threshold: |delta-rho| >= {active_threshold}"
                         f" {'(auto-selected)' if DIFF_THRESHOLD_AUTO else '(config)'}")
    summary_lines.append(f"Resolution: {best_res}"
                         f" {'(auto-selected)' if RESOLUTION_AUTO else '(config)'}")
    summary_lines.append(f"Modularity: {best_mod:.4f}")
    summary_lines.append(f"ARI stability: {best_ari:.4f}")
    summary_lines.append(f"Total genes: {total_genes}")
    summary_lines.append(f"Communities: {len(final_sizes)}")
    summary_lines.append("")

    for comm in sorted(final_sizes.keys()):
        genes_in_comm = [g for g, c in final_partition.items() if c == comm]
        summary_lines.append(f"Community {comm} ({len(genes_in_comm)} genes):")
        for a3 in A3_GENES:
            if a3 in genes_in_comm:
                summary_lines.append(f"  ** {A3_ID_TO_ALIAS.get(a3, a3)} ({a3}) **")
        for bm in BIOMARKERS:
            if bm in genes_in_comm:
                summary_lines.append(f"  [biomarker] {bm}")
        summary_lines.append("")

    summary_path = os.path.join(DIR_04_COMMUNITIES, "SC_community_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines))
    log(f"  [SAVE] Summary -> {summary_path}")
    for line in summary_lines:
        print(line)

    # =========================================================================
    # 9. Sweep plots
    # =========================================================================
    banner("[STEP 03.9] Resolution sweep plots")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].errorbar(sweep_df["resolution"], sweep_df["modularity_mean"],
                     yerr=sweep_df["modularity_std"], fmt="o-", c="steelblue")
    axes[0].axvline(best_res, ls="--", c="red", alpha=0.5)
    axes[0].set_xlabel("Resolution"); axes[0].set_ylabel("Modularity")
    axes[0].set_title("Modularity vs Resolution")

    axes[1].errorbar(sweep_df["resolution"], sweep_df["ari_mean"],
                     yerr=sweep_df["ari_std"], fmt="o-", c="darkgreen")
    axes[1].axvline(best_res, ls="--", c="red", alpha=0.5)
    axes[1].set_xlabel("Resolution"); axes[1].set_ylabel("ARI")
    axes[1].set_title("Partition Stability (ARI)")

    axes[2].plot(sweep_df["resolution"], sweep_df["n_communities"], "o-", c="purple")
    axes[2].axvline(best_res, ls="--", c="red", alpha=0.5)
    axes[2].set_xlabel("Resolution"); axes[2].set_ylabel("N Communities")
    axes[2].set_title("Community Count")

    plt.tight_layout()
    plt.savefig(os.path.join(DIR_04_COMMUNITIES, "SC_sweep_plots.png"), dpi=300)
    plt.close()
    log(f"  [SAVE] Sweep plots")

    elapsed = datetime.now() - start_time
    banner(f"[STEP 03 COMPLETE] {len(final_sizes)} communities, "
           f"{total_genes} genes | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
