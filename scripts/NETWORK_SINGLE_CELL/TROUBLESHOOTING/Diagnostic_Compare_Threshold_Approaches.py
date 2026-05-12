#!/usr/bin/env python3
"""
Diagnostic_Compare_Threshold_Approaches.py
=============================================

Side-by-side comparison of threshold selection approaches for the three
Figure 4 single-cell networks:

  v2 (ENRICHMENT): Harris interactor enrichment in A3 2-hop neighborhood.
     Results already computed and stored in each network's THRESHOLD_SELECTION/
     and 04_communities/ directories.

  v3 (SIMPLE MAX): Highest DIFF threshold where any A3 gene that passed DE
     on its own merit (raw p < 0.05, no force-keeping) has degree >= 1 in
     the network. Then standard Leiden resolution selection (modularity jump
     + ARI floor). Matches the TCGA bulk pipeline logic.

Also reports the TCGA bulk network (Fig 2) as a reference point.

For v3, this script actually builds the graph and runs community detection
at the selected threshold so we get real community assignments, A3 placement,
and Harris interactor counts.

Output:
  data/FIG_4/DIAGNOSTIC_AUDIT/approach_comparison_report.txt
  data/FIG_4/DIAGNOSTIC_AUDIT/approach_comparison_summary.csv

Usage:
  conda run -n NETWORK python Diagnostic_Compare_Threshold_Approaches.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import leidenalg as la
import igraph as ig
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_4")
FIG2_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_2")

SC_NETWORKS = {
    "SBS2_VS_CNV": {
        "label": "SBS2-HIGH vs CNV-HIGH",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_CNV"),
    },
    "SBS2_VS_NORMAL": {
        "label": "SBS2-HIGH vs NORMAL",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL"),
    },
    "CNV_VS_NORMAL": {
        "label": "CNV-HIGH vs NORMAL",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL"),
    },
}

A3_GENES = [
    "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
    "APOBEC3F", "APOBEC3G", "APOBEC3H",
]
A3_SHORT = {
    "APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
    "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
    "APOBEC3H": "A3H",
}

HARRIS_ALL_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors_A3B_only.txt")

# Threshold sweep (fine grid)
SWEEP_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                    0.75, 0.80, 0.85, 0.90]

# Resolution sweep for Leiden (matching TCGA pipeline)
RESOLUTION_SWEEP = [0.2, 0.4, 0.6, 0.8, 1.0]
N_LEIDEN_RUNS = 10
ARI_FLOOR = 0.65
RANDOM_SEED = 42
MIN_COMMUNITY_SIZE = 5

RAW_P_THRESHOLD = 0.05

OUTPUT_DIR = os.path.join(FIG4_ROOT, "DIAGNOSTIC_AUDIT")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# LOGGING
# =============================================================================

report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log()
    log(char * 80)
    log(f"  {title}")
    log(char * 80)


# =============================================================================
# HELPERS
# =============================================================================

def load_harris():
    harris_all = set()
    harris_a3b = set()
    if os.path.exists(HARRIS_ALL_PATH):
        try:
            df = pd.read_csv(HARRIS_ALL_PATH, sep='\t')
            harris_all = set(df['gene_symbol'].dropna().values)
        except Exception:
            pass
    if os.path.exists(HARRIS_A3B_PATH):
        try:
            df = pd.read_csv(HARRIS_A3B_PATH, sep='\t')
            harris_a3b = set(df['gene_symbol'].dropna().values)
        except Exception:
            pass
    return harris_all, harris_a3b


def load_diff_matrix(net_dir):
    path = os.path.join(net_dir, "03_correlation_networks",
                        "corr_matrices", "SC_corr_DIFF.pkl")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def load_de_stats(net_dir):
    """Load DE stats to identify which A3 genes passed DE on their own."""
    path = os.path.join(net_dir, "02_differential_expression",
                        "SC_diffexpr_stats.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def build_graph(corr_diff, threshold):
    """Build graph at threshold, remove isolates."""
    genes = list(corr_diff.index)
    abs_diff = np.abs(corr_diff.values)
    np.fill_diagonal(abs_diff, 0)
    mask = np.triu(abs_diff >= threshold, k=1)
    rows, cols = np.where(mask)

    G = nx.Graph()
    G.add_nodes_from(genes)
    for r, c in zip(rows, cols):
        G.add_edge(genes[r], genes[c],
                   weight=float(corr_diff.values[r, c]),
                   abs_weight=float(abs_diff[r, c]))

    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    return G


def run_leiden(G, resolution, n_runs=N_LEIDEN_RUNS, seed=RANDOM_SEED):
    """Run Leiden community detection, return partition and modularity."""
    # Convert NetworkX to igraph
    node_list = list(G.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_list)}

    ig_G = ig.Graph()
    ig_G.add_vertices(len(node_list))
    edges = [(node_to_idx[u], node_to_idx[v]) for u, v in G.edges()]
    ig_G.add_edges(edges)

    weights = [abs(G[u][v].get('weight', G[u][v].get('abs_weight', 1.0)))
               for u, v in G.edges()]
    ig_G.es['weight'] = weights

    best_Q = -1
    best_partition = None
    all_Q = []

    for i in range(n_runs):
        part = la.find_partition(
            ig_G, la.RBConfigurationVertexPartition,
            weights=weights, resolution_parameter=resolution,
            seed=seed + i
        )
        Q = part.modularity
        all_Q.append(Q)
        if Q > best_Q:
            best_Q = Q
            best_partition = part

    # Convert back to gene -> community mapping
    gene_to_comm = {}
    for comm_idx, members in enumerate(best_partition):
        for node_idx in members:
            gene_to_comm[node_list[node_idx]] = comm_idx

    return gene_to_comm, best_Q, np.std(all_Q)


def compute_ari(G, resolution, n_runs=N_LEIDEN_RUNS, seed=RANDOM_SEED):
    """Compute ARI stability across multiple Leiden runs."""
    from sklearn.metrics import adjusted_rand_score

    node_list = list(G.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_list)}

    ig_G = ig.Graph()
    ig_G.add_vertices(len(node_list))
    edges = [(node_to_idx[u], node_to_idx[v]) for u, v in G.edges()]
    ig_G.add_edges(edges)
    weights = [abs(G[u][v].get('weight', G[u][v].get('abs_weight', 1.0)))
               for u, v in G.edges()]
    ig_G.es['weight'] = weights

    partitions = []
    for i in range(n_runs):
        part = la.find_partition(
            ig_G, la.RBConfigurationVertexPartition,
            weights=weights, resolution_parameter=resolution,
            seed=seed + i
        )
        labels = [0] * len(node_list)
        for comm_idx, members in enumerate(part):
            for node_idx in members:
                labels[node_idx] = comm_idx
        partitions.append(labels)

    # Average pairwise ARI
    aris = []
    for i in range(len(partitions)):
        for j in range(i + 1, len(partitions)):
            aris.append(adjusted_rand_score(partitions[i], partitions[j]))

    return np.mean(aris) if aris else 1.0


def merge_small_communities(gene_to_comm, G, min_size=MIN_COMMUNITY_SIZE):
    """Merge communities smaller than min_size into nearest large community."""
    comm_sizes = {}
    for g, c in gene_to_comm.items():
        comm_sizes[c] = comm_sizes.get(c, 0) + 1

    small = {c for c, s in comm_sizes.items() if s < min_size}
    large = {c for c, s in comm_sizes.items() if s >= min_size}

    if not small:
        return gene_to_comm

    merged = dict(gene_to_comm)
    for gene, comm in gene_to_comm.items():
        if comm in small and gene in G:
            # Find which large community this gene connects to most
            neighbor_comms = {}
            for neighbor in G.neighbors(gene):
                nc = gene_to_comm.get(neighbor)
                if nc is not None and nc in large:
                    neighbor_comms[nc] = neighbor_comms.get(nc, 0) + 1
            if neighbor_comms:
                merged[gene] = max(neighbor_comms, key=neighbor_comms.get)

    # Renumber communities sequentially
    unique = sorted(set(merged.values()))
    remap = {old: new for new, old in enumerate(unique)}
    return {g: remap[c] for g, c in merged.items()}


# =============================================================================
# V3: SIMPLE MAX THRESHOLD APPROACH
# =============================================================================

def v3_simple_max(corr_diff, de_stats, net_name, harris_all, harris_a3b):
    """
    V3 approach: highest DIFF threshold keeping any naturally-DE A3 gene
    with degree >= 1, then standard Leiden resolution selection.
    """
    banner(f"V3 SIMPLE MAX: {net_name}")

    genes_in_matrix = set(corr_diff.index)

    # Step 1: Which A3 genes passed DE on their own (raw p < 0.05)?
    log(f"\n  A3 genes in DE stats:")
    a3_natural_de = []
    for a3 in A3_GENES:
        row = de_stats[de_stats['gene'] == a3]
        if len(row) == 0:
            log(f"    {A3_SHORT[a3]}: NOT in DE stats")
            continue
        p = float(row['p_value'].values[0])
        passed = p < RAW_P_THRESHOLD
        status = f"PASSED (p={p:.4f})" if passed else f"force-kept (p={p:.4f})"
        log(f"    {A3_SHORT[a3]}: {status}")
        if passed:
            a3_natural_de.append(a3)

    log(f"\n  A3 genes passing DE naturally: {len(a3_natural_de)} "
        f"({', '.join(A3_SHORT[a] for a in a3_natural_de)})")

    if not a3_natural_de:
        log(f"  No A3 genes passed DE naturally. Cannot apply v3.")
        return None

    # Filter to those present in the DIFF matrix
    a3_in_matrix = [a3 for a3 in a3_natural_de if a3 in genes_in_matrix]
    log(f"  A3 in DIFF matrix: {len(a3_in_matrix)} "
        f"({', '.join(A3_SHORT[a] for a in a3_in_matrix)})")

    if not a3_in_matrix:
        log(f"  No naturally-DE A3 genes in DIFF matrix.")
        return None

    # Step 2: Sweep thresholds, find max where any A3 has degree >= 1
    log(f"\n  Threshold sweep (A3 degree tracking):")
    log(f"  {'Thresh':>7s}  {'Nodes':>6s}  {'Edges':>7s}  "
        f"{'Comp':>5s}  {'LCC':>5s}  {'A3 connected':>40s}")
    log(f"  {'-----':>7s}  {'-----':>6s}  {'-----':>7s}  "
        f"{'----':>5s}  {'---':>5s}  {'-'*40}")

    best_thresh = None
    thresh_details = {}

    for t in SWEEP_THRESHOLDS:
        G = build_graph(corr_diff, t)
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()

        if n_nodes > 0:
            components = list(nx.connected_components(G))
            n_comp = len(components)
            lcc = max(components, key=len)
            n_lcc = len(lcc)
        else:
            n_comp = 0
            n_lcc = 0
            lcc = set()

        # Check A3 degrees
        a3_connected = {}
        for a3 in a3_in_matrix:
            if a3 in G:
                deg = G.degree(a3)
                if deg >= 1:
                    a3_connected[A3_SHORT[a3]] = deg

        # Harris interactors
        harris_present = harris_all & set(G.nodes())

        conn_str = ", ".join(f"{k}={v}" for k, v in a3_connected.items()) if a3_connected else "none"
        log(f"  {t:>7.2f}  {n_nodes:>6d}  {n_edges:>7d}  "
            f"{n_comp:>5d}  {n_lcc:>5d}  {conn_str}")

        thresh_details[t] = {
            'nodes': n_nodes, 'edges': n_edges, 'comp': n_comp,
            'lcc': n_lcc, 'a3_connected': a3_connected,
            'harris_in_net': len(harris_present),
        }

        if a3_connected:
            best_thresh = t  # Keep updating; we want the highest

    if best_thresh is None:
        log(f"\n  No threshold has any A3 gene with degree >= 1.")
        return None

    log(f"\n  Selected threshold: {best_thresh}")
    details = thresh_details[best_thresh]
    log(f"    A3 connected: {details['a3_connected']}")
    log(f"    Network: {details['nodes']} nodes, {details['edges']} edges")
    log(f"    LCC: {details['lcc']}")
    log(f"    Harris in network: {details['harris_in_net']}")

    # Step 3: Build graph at selected threshold, extract LCC
    G = build_graph(corr_diff, best_thresh)
    components = list(nx.connected_components(G))
    lcc = max(components, key=len)
    G_lcc = G.subgraph(lcc).copy()

    log(f"\n  LCC extracted: {G_lcc.number_of_nodes()} nodes, "
        f"{G_lcc.number_of_edges()} edges")

    # Check which A3 genes are in the LCC
    a3_in_lcc = [a3 for a3 in a3_in_matrix if a3 in lcc]
    log(f"  A3 in LCC: {len(a3_in_lcc)} "
        f"({', '.join(A3_SHORT[a] for a in a3_in_lcc)})")

    # Harris in LCC
    harris_in_lcc = harris_all & set(G_lcc.nodes())
    harris_a3b_in_lcc = harris_a3b & set(G_lcc.nodes())
    log(f"  Harris in LCC: {len(harris_in_lcc)}")
    log(f"  Harris A3B-sp in LCC: {len(harris_a3b_in_lcc)}")

    # Step 4: Resolution sweep with standard Leiden
    log(f"\n  Resolution sweep:")
    log(f"  {'Res':>5s}  {'Comms':>6s}  {'Q':>8s}  {'ARI':>8s}")
    log(f"  {'---':>5s}  {'-----':>6s}  {'---':>8s}  {'---':>8s}")

    res_results = []
    for res in RESOLUTION_SWEEP:
        gene_to_comm, Q, Q_std = run_leiden(G_lcc, res)
        n_comms = len(set(gene_to_comm.values()))
        ari = compute_ari(G_lcc, res)
        log(f"  {res:>5.1f}  {n_comms:>6d}  {Q:>8.4f}  {ari:>8.4f}")
        res_results.append({
            'resolution': res, 'communities': n_comms,
            'Q': Q, 'ARI': ari,
        })

    # Select resolution: largest Q jump with ARI >= floor
    res_df = pd.DataFrame(res_results)
    res_df['dQ'] = res_df['Q'].diff().fillna(0)

    log(f"\n  Modularity jump analysis:")
    log(f"  {'Res':>5s}  {'Comms':>6s}  {'Q':>8s}  {'dQ':>8s}  {'ARI':>8s}")
    for _, row in res_df.iterrows():
        log(f"  {row['resolution']:>5.1f}  {row['communities']:>6.0f}  "
            f"{row['Q']:>8.4f}  {row['dQ']:>+8.4f}  {row['ARI']:>8.4f}")

    # Best: largest dQ with ARI >= floor
    passing = res_df[res_df['ARI'] >= ARI_FLOOR].copy()
    if len(passing) > 0:
        best_res_idx = passing['dQ'].idxmax()
        best_res = float(res_df.loc[best_res_idx, 'resolution'])
        best_Q = float(res_df.loc[best_res_idx, 'Q'])
        best_ARI = float(res_df.loc[best_res_idx, 'ARI'])
        best_comms = int(res_df.loc[best_res_idx, 'communities'])
    else:
        # Fallback: highest ARI
        best_res_idx = res_df['ARI'].idxmax()
        best_res = float(res_df.loc[best_res_idx, 'resolution'])
        best_Q = float(res_df.loc[best_res_idx, 'Q'])
        best_ARI = float(res_df.loc[best_res_idx, 'ARI'])
        best_comms = int(res_df.loc[best_res_idx, 'communities'])

    log(f"\n  Selected resolution: {best_res} "
        f"(Q={best_Q:.4f}, ARI={best_ARI:.4f}, comms={best_comms})")

    # Run final partition at selected resolution
    gene_to_comm, Q_final, _ = run_leiden(G_lcc, best_res)
    gene_to_comm = merge_small_communities(gene_to_comm, G_lcc)
    n_final_comms = len(set(gene_to_comm.values()))

    # Community sizes
    comm_sizes = {}
    for g, c in gene_to_comm.items():
        comm_sizes[c] = comm_sizes.get(c, 0) + 1

    log(f"\n  Final communities: {n_final_comms}")
    for c in sorted(comm_sizes):
        log(f"    C{c}: {comm_sizes[c]} genes")

    # A3 placement
    log(f"\n  A3 gene placement:")
    a3_in_communities = {}
    for a3 in A3_GENES:
        if a3 in gene_to_comm:
            c = gene_to_comm[a3]
            deg = G_lcc.degree(a3)
            log(f"    {A3_SHORT[a3]}: C{c} (degree={deg})")
            a3_in_communities[A3_SHORT[a3]] = {'community': c, 'degree': deg}

    # Harris interactor placement
    harris_in_comms = {}
    harris_names_in_comms = []
    for h in harris_all:
        if h in gene_to_comm:
            c = gene_to_comm[h]
            harris_in_comms[c] = harris_in_comms.get(c, 0) + 1
            harris_names_in_comms.append(h)

    log(f"\n  Harris interactors in communities: {len(harris_names_in_comms)}")
    for c in sorted(harris_in_comms):
        log(f"    C{c}: {harris_in_comms[c]}")

    if harris_names_in_comms:
        log(f"  Names: {', '.join(sorted(harris_names_in_comms)[:20])}"
            f"{'...' if len(harris_names_in_comms) > 20 else ''}")

    return {
        'threshold': best_thresh,
        'resolution': best_res,
        'nodes_lcc': G_lcc.number_of_nodes(),
        'edges_lcc': G_lcc.number_of_edges(),
        'communities': n_final_comms,
        'Q': Q_final,
        'ARI': best_ARI,
        'a3_in_comms': len(a3_in_communities),
        'a3_details': a3_in_communities,
        'harris_in_comms': len(harris_names_in_comms),
        'harris_a3b_in_comms': sum(1 for h in harris_names_in_comms if h in harris_a3b),
        'largest_comm': max(comm_sizes.values()) if comm_sizes else 0,
    }


# =============================================================================
# LOAD V2 RESULTS
# =============================================================================

def load_v2_results(net_dir, harris_all, harris_a3b):
    """Load the v2 enrichment results from the existing pipeline output."""
    result = {}

    # Threshold from selection
    param_path = os.path.join(net_dir, "THRESHOLD_SELECTION",
                              "selected_parameters.txt")
    if os.path.exists(param_path):
        with open(param_path) as f:
            for line in f:
                if line.startswith("DIFF_THRESHOLD="):
                    result['threshold'] = float(line.split('=')[1].strip())
                if line.startswith("SELECTION_TIER="):
                    result['tier'] = int(line.split('=')[1].strip())

    # Community assignments
    part_path = os.path.join(net_dir, "04_communities",
                             "SC_best_partition.csv")
    if os.path.exists(part_path):
        part_df = pd.read_csv(part_path)
        gene_to_comm = dict(zip(part_df['gene'], part_df['community']))
        n_comms = len(set(gene_to_comm.values()))
        result['nodes_lcc'] = len(gene_to_comm)
        result['communities'] = n_comms

        # Community sizes
        comm_sizes = {}
        for g, c in gene_to_comm.items():
            comm_sizes[c] = comm_sizes.get(c, 0) + 1
        result['largest_comm'] = max(comm_sizes.values()) if comm_sizes else 0

        # A3 placement
        a3_details = {}
        for a3 in A3_GENES:
            if a3 in gene_to_comm:
                a3_details[A3_SHORT[a3]] = {'community': gene_to_comm[a3]}
        result['a3_in_comms'] = len(a3_details)
        result['a3_details'] = a3_details

        # Harris
        harris_in = [h for h in harris_all if h in gene_to_comm]
        result['harris_in_comms'] = len(harris_in)
        result['harris_a3b_in_comms'] = sum(1 for h in harris_in if h in harris_a3b)

    # Parameters
    param_path2 = os.path.join(net_dir, "04_communities",
                               "SC_selected_parameters.txt")
    if os.path.exists(param_path2):
        with open(param_path2) as f:
            for line in f:
                if 'resolution' in line.lower() and '=' in line:
                    try:
                        result['resolution'] = float(line.split('=')[1].split('#')[0].strip())
                    except ValueError:
                        pass
                if 'modularity' in line.lower() and '=' in line:
                    try:
                        result['Q'] = float(line.split('=')[1].split('#')[0].strip())
                    except ValueError:
                        pass
                if 'ari' in line.lower() and '=' in line:
                    try:
                        result['ARI'] = float(line.split('=')[1].split('#')[0].strip())
                    except ValueError:
                        pass

    return result


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("THRESHOLD APPROACH COMPARISON: v2 (Enrichment) vs v3 (Simple Max)")
    log(f"  Start: {t0}")

    harris_all, harris_a3b = load_harris()
    log(f"  Harris (all): {len(harris_all)}, (A3B-sp): {len(harris_a3b)}")

    comparison_rows = []

    # ---- TCGA bulk reference ----
    banner("REFERENCE: TCGA BULK (Figure 2)")
    log(f"  Threshold: 0.70 (topology-driven)")
    log(f"  Genes: 616, Communities: 14, Q: 0.61")
    log(f"  A3B: C4 (degree=1, passenger)")
    log(f"  Harris interactors: 0")

    comparison_rows.append({
        'network': 'TCGA_BULK',
        'approach': 'TCGA (topology)',
        'threshold': 0.70,
        'nodes': 616,
        'communities': 14,
        'Q': 0.61,
        'ARI': 0.63,
        'a3_in_comms': 1,
        'harris_in_comms': 0,
        'harris_a3b_in_comms': 0,
        'largest_comm': 88,
    })

    # ---- SC networks ----
    for net_name, net_info in SC_NETWORKS.items():
        net_dir = net_info['dir']

        banner(f"NETWORK: {net_name}", char="*")
        log(f"  {net_info['label']}")

        # Load data
        corr_diff = load_diff_matrix(net_dir)
        de_stats = load_de_stats(net_dir)

        if corr_diff is None or de_stats is None:
            log(f"  SKIPPING (missing data)")
            continue

        # ---- V2 results (already computed) ----
        banner(f"V2 ENRICHMENT RESULTS: {net_name}")
        v2 = load_v2_results(net_dir, harris_all, harris_a3b)
        log(f"  Threshold: {v2.get('threshold', '?')}")
        log(f"  Nodes: {v2.get('nodes_lcc', '?')}")
        log(f"  Communities: {v2.get('communities', '?')}")
        log(f"  Q: {v2.get('Q', '?')}")
        log(f"  ARI: {v2.get('ARI', '?')}")
        log(f"  A3 in communities: {v2.get('a3_in_comms', '?')}")
        if v2.get('a3_details'):
            for a3, info in v2['a3_details'].items():
                log(f"    {a3}: C{info['community']}")
        log(f"  Harris in communities: {v2.get('harris_in_comms', '?')}")
        log(f"  Largest community: {v2.get('largest_comm', '?')}")

        v2_row = {
            'network': net_name, 'approach': 'v2 (enrichment)',
            **{k: v2.get(k, 0) for k in ['threshold', 'nodes_lcc',
               'communities', 'Q', 'ARI', 'a3_in_comms',
               'harris_in_comms', 'harris_a3b_in_comms', 'largest_comm']}
        }
        comparison_rows.append(v2_row)

        # ---- V3 results (compute now) ----
        v3 = v3_simple_max(corr_diff, de_stats, net_name, harris_all, harris_a3b)

        if v3 is not None:
            v3_row = {
                'network': net_name, 'approach': 'v3 (simple max)',
                **{k: v3.get(k, 0) for k in ['threshold', 'nodes_lcc',
                   'communities', 'Q', 'ARI', 'a3_in_comms',
                   'harris_in_comms', 'harris_a3b_in_comms', 'largest_comm']}
            }
            comparison_rows.append(v3_row)

    # =================================================================
    # COMPARISON TABLE
    # =================================================================
    banner("SIDE-BY-SIDE COMPARISON")

    comp_df = pd.DataFrame(comparison_rows)
    log(f"\n  {'Network':>18s}  {'Approach':>20s}  {'Thresh':>6s}  "
        f"{'Nodes':>6s}  {'Comms':>5s}  {'Q':>6s}  {'ARI':>6s}  "
        f"{'A3':>3s}  {'Harris':>7s}  {'MaxC':>5s}")
    log(f"  {'-'*18}  {'-'*20}  {'-'*6}  "
        f"{'-'*6}  {'-'*5}  {'-'*6}  {'-'*6}  "
        f"{'-'*3}  {'-'*7}  {'-'*5}")

    for _, row in comp_df.iterrows():
        q_str = f"{row.get('Q', 0):.2f}" if row.get('Q') else "?"
        ari_str = f"{row.get('ARI', 0):.2f}" if row.get('ARI') else "?"
        log(f"  {row['network']:>18s}  {row['approach']:>20s}  "
            f"{row.get('threshold', 0):>6.2f}  "
            f"{row.get('nodes_lcc', 0):>6.0f}  "
            f"{row.get('communities', 0):>5.0f}  "
            f"{q_str:>6s}  {ari_str:>6s}  "
            f"{row.get('a3_in_comms', 0):>3.0f}  "
            f"{row.get('harris_in_comms', 0):>7.0f}  "
            f"{row.get('largest_comm', 0):>5.0f}")

    # Save
    csv_path = os.path.join(OUTPUT_DIR, "approach_comparison_summary.csv")
    comp_df.to_csv(csv_path, index=False)
    log(f"\n  [SAVE] Summary CSV -> {csv_path}")

    report_path = os.path.join(OUTPUT_DIR, "approach_comparison_report.txt")
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    log(f"  [SAVE] Report -> {report_path}")

    banner("COMPARISON COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")


if __name__ == "__main__":
    main()
