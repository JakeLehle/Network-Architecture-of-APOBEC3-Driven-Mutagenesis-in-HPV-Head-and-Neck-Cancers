#!/usr/bin/env python3
"""
Diagnostic_A3_Network_Feature_Sweep.py
=========================================

Comprehensive sweep of network connectivity features for A3 genes and
Harris interactors across all DIFF thresholds and all three Figure 4
single-cell networks.

Four metrics, chosen because they capture indirect/diffuse connectivity
that standard degree-based measures miss:

  1. TOPOLOGICAL OVERLAP (WGCNA-style TOM)
     Measures shared-neighbor structure between A3 genes and Harris
     interactors. High TOM means two genes wire into the same neighborhood
     even without a direct edge. Well-established in gene network literature.

  2. 2-HOP REACHABILITY
     Counts Harris interactors within 2 hops of each A3 gene. Simple,
     interpretable, and directly measures the "funnel" structure where
     A3 genes have few direct edges but those edges connect to high-degree
     hubs that reach the interactor network.

  3. COMMUNICABILITY (matrix exponential)
     Counts ALL walks between two nodes, weighted by walk length. Captures
     connectivity that shortest-path measures miss. Handles disconnected
     components gracefully (communicability -> 0 smoothly). Uses
     scipy.sparse.linalg.expm_multiply for efficiency on large graphs.

  4. K-CORE MEMBERSHIP
     Identifies the densest subgraph each node belongs to. Even with
     degree=2, a node can sit at the edge of a high-k core if its
     neighbors are highly connected. Reveals whether A3 genes are
     structurally embedded in dense regions.

Output:
  data/FIG_4/DIAGNOSTIC_AUDIT/A3_feature_sweep_report.txt
  data/FIG_4/DIAGNOSTIC_AUDIT/A3_feature_sweep_summary.csv
  data/FIG_4/DIAGNOSTIC_AUDIT/A3_feature_sweep_plots.png

Usage:
  conda run -n NETWORK python Diagnostic_A3_Network_Feature_Sweep.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply, eigsh
from scipy.sparse import csr_matrix
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_4")

NETWORKS = {
    "SBS2_VS_CNV": {
        "label": "SBS2-HIGH vs CNV-HIGH",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_CNV"),
        "auto_threshold": 0.65,
    },
    "SBS2_VS_NORMAL": {
        "label": "SBS2-HIGH vs NORMAL",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL"),
        "auto_threshold": 0.65,
    },
    "CNV_VS_NORMAL": {
        "label": "CNV-HIGH vs NORMAL",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL"),
        "auto_threshold": 0.70,
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
# Focus on A3A and A3B for detailed reporting
FOCAL_A3 = ["APOBEC3A", "APOBEC3B"]

HARRIS_ALL_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors_A3B_only.txt")

SWEEP_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                    0.75, 0.80, 0.85, 0.90]

OUTPUT_DIR = os.path.join(FIG4_ROOT, "DIAGNOSTIC_AUDIT")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Communicability: max component size for dense expm (larger uses sparse)
COMM_DENSE_MAX = 3000
# Communicability: normalize adjacency to prevent overflow
COMM_NORMALIZE = True

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
# LOAD HELPERS
# =============================================================================

def load_harris():
    """Load Harris A3 interactors from TSV files."""
    harris_all = set()
    harris_a3b = set()
    if os.path.exists(HARRIS_ALL_PATH):
        try:
            df = pd.read_csv(HARRIS_ALL_PATH, sep='\t')
            harris_all = set(df['gene_symbol'].dropna().values)
        except Exception as e:
            log(f"  WARNING: Could not load Harris interactors: {e}")
    if os.path.exists(HARRIS_A3B_PATH):
        try:
            df = pd.read_csv(HARRIS_A3B_PATH, sep='\t')
            harris_a3b = set(df['gene_symbol'].dropna().values)
        except Exception as e:
            log(f"  WARNING: Could not load A3B interactors: {e}")
    return harris_all, harris_a3b


def load_diff_matrix(net_dir):
    """Load DIFF correlation matrix pickle from a network directory."""
    corr_path = os.path.join(net_dir, "03_correlation_networks",
                             "corr_matrices", "SC_corr_DIFF.pkl")
    if not os.path.exists(corr_path):
        log(f"  WARNING: DIFF matrix not found: {corr_path}")
        return None
    with open(corr_path, "rb") as f:
        return pickle.load(f)


def build_graph_at_threshold(corr_diff, threshold):
    """Build undirected graph from DIFF matrix at given |delta-rho| threshold."""
    genes = list(corr_diff.index)
    n = len(genes)
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    abs_diff = np.abs(corr_diff.values)
    np.fill_diagonal(abs_diff, 0)

    # Upper triangle edges above threshold
    mask = np.triu(abs_diff >= threshold, k=1)
    rows, cols = np.where(mask)

    G = nx.Graph()
    G.add_nodes_from(genes)
    for r, c in zip(rows, cols):
        G.add_edge(genes[r], genes[c],
                   weight=float(corr_diff.values[r, c]),
                   abs_weight=float(abs_diff[r, c]))

    # Remove isolates
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)

    return G, gene_to_idx, genes


# =============================================================================
# METRIC 1: TOPOLOGICAL OVERLAP
# =============================================================================

def compute_tom_pairs(G, a3_genes_in_net, harris_in_net):
    """
    Compute topological overlap between A3 genes and Harris interactors.

    TOM(u, v) = (shared_neighbors(u,v) + a(u,v)) / (min(deg(u), deg(v)) + 1 - a(u,v))

    Returns dict with summary stats and top pairs.
    """
    if not a3_genes_in_net or not harris_in_net:
        return {'max_tom': 0, 'mean_tom': 0, 'n_pairs': 0,
                'top_pairs': [], 'a3a_max': 0, 'a3b_max': 0}

    neighbor_sets = {}
    for node in set(a3_genes_in_net) | set(harris_in_net):
        if node in G:
            neighbor_sets[node] = set(G.neighbors(node))
        else:
            neighbor_sets[node] = set()

    all_toms = []
    a3a_toms = []
    a3b_toms = []
    top_pairs = []

    for a3 in a3_genes_in_net:
        if a3 not in G:
            continue
        deg_a3 = G.degree(a3)
        n_a3 = neighbor_sets.get(a3, set())

        for harris in harris_in_net:
            if harris not in G or harris == a3:
                continue
            deg_h = G.degree(harris)
            n_h = neighbor_sets.get(harris, set())

            shared = len(n_a3 & n_h)
            adjacent = 1 if G.has_edge(a3, harris) else 0
            denom = min(deg_a3, deg_h) + 1 - adjacent
            tom = (shared + adjacent) / denom if denom > 0 else 0

            all_toms.append(tom)
            top_pairs.append((a3, harris, tom, shared, adjacent))

            if a3 == "APOBEC3A":
                a3a_toms.append(tom)
            elif a3 == "APOBEC3B":
                a3b_toms.append(tom)

    if not all_toms:
        return {'max_tom': 0, 'mean_tom': 0, 'n_pairs': 0,
                'top_pairs': [], 'a3a_max': 0, 'a3b_max': 0}

    # Sort top pairs by TOM
    top_pairs.sort(key=lambda x: -x[2])

    return {
        'max_tom': max(all_toms),
        'mean_tom': np.mean(all_toms),
        'n_pairs': len(all_toms),
        'n_nonzero': sum(1 for t in all_toms if t > 0),
        'top_pairs': top_pairs[:15],
        'a3a_max': max(a3a_toms) if a3a_toms else 0,
        'a3b_max': max(a3b_toms) if a3b_toms else 0,
        'a3a_mean': np.mean(a3a_toms) if a3a_toms else 0,
        'a3b_mean': np.mean(a3b_toms) if a3b_toms else 0,
    }


# =============================================================================
# METRIC 2: 2-HOP REACHABILITY
# =============================================================================

def compute_reachability(G, a3_genes_in_net, harris_in_net, harris_a3b_set):
    """
    Count Harris interactors within 2 hops of each A3 gene.
    Returns per-A3 results and aggregate summary.
    """
    results = {}
    all_harris_reached = set()

    for a3 in a3_genes_in_net:
        if a3 not in G:
            continue
        short = A3_SHORT.get(a3, a3)

        # 1-hop
        n1 = set(G.neighbors(a3))
        harris_1hop = n1 & harris_in_net

        # 2-hop
        n2 = set()
        for neighbor in n1:
            n2 |= set(G.neighbors(neighbor))
        n2 -= n1
        n2 -= {a3}
        harris_2hop = n2 & harris_in_net

        harris_total = harris_1hop | harris_2hop
        harris_a3b_reached = harris_total & harris_a3b_set
        all_harris_reached |= harris_total

        results[short] = {
            'degree': G.degree(a3),
            'n1_size': len(n1),
            'n2_size': len(n2),
            'harris_1hop': len(harris_1hop),
            'harris_2hop': len(harris_2hop),
            'harris_total': len(harris_total),
            'harris_a3b_reached': len(harris_a3b_reached),
            'harris_1hop_names': sorted(harris_1hop),
            'harris_2hop_names': sorted(harris_2hop),
        }

    return results, len(all_harris_reached)


# =============================================================================
# METRIC 3: COMMUNICABILITY
# =============================================================================

def compute_communicability(G, a3_genes_in_net, harris_in_net):
    """
    Compute communicability from A3 genes to Harris interactors.

    For small components: dense scipy.linalg.expm
    For large components: scipy.sparse.linalg.expm_multiply

    The adjacency matrix is normalized by spectral radius to prevent overflow.
    """
    from scipy.linalg import expm as dense_expm

    results = {'a3a_sum': 0, 'a3b_sum': 0, 'a3a_max': 0, 'a3b_max': 0,
               'total_sum': 0, 'harris_reached': 0, 'top_pairs': []}

    if not a3_genes_in_net or not harris_in_net:
        return results

    # Find which component each A3 gene is in
    components = list(nx.connected_components(G))
    node_to_comp = {}
    for i, comp in enumerate(components):
        for node in comp:
            node_to_comp[node] = i

    top_pairs = []

    for a3 in a3_genes_in_net:
        if a3 not in G:
            continue
        short = A3_SHORT.get(a3, a3)

        comp_idx = node_to_comp.get(a3)
        if comp_idx is None:
            continue

        comp_nodes = sorted(components[comp_idx])
        comp_size = len(comp_nodes)

        # Harris interactors in same component
        harris_in_comp = [h for h in harris_in_net if h in set(comp_nodes)]
        if not harris_in_comp:
            continue

        node_to_local = {n: i for i, n in enumerate(comp_nodes)}
        a3_local = node_to_local[a3]

        # Build adjacency matrix for this component (unweighted for stability)
        A = nx.adjacency_matrix(G.subgraph(comp_nodes),
                                nodelist=comp_nodes, weight=None)
        A = A.astype(float)

        # Normalize by spectral radius to prevent overflow
        if COMM_NORMALIZE and comp_size > 2:
            try:
                eigenvalues = eigsh(A.astype(float), k=1, which='LM',
                                    return_eigenvectors=False)
                spectral_radius = abs(eigenvalues[0])
                if spectral_radius > 1e-10:
                    A = A / spectral_radius
            except Exception:
                # Fallback: normalize by max degree
                max_deg = max(dict(G.subgraph(comp_nodes).degree()).values())
                if max_deg > 0:
                    A = A / max_deg

        # Compute communicability from A3 gene to all others
        if comp_size <= COMM_DENSE_MAX:
            # Dense: full matrix exponential
            try:
                E = dense_expm(A.toarray())
                comm_vec = E[a3_local, :]
            except Exception as e:
                log(f"    WARNING: dense expm failed for {short}: {e}")
                continue
        else:
            # Sparse: expm_multiply with unit vector
            try:
                e_i = np.zeros(comp_size)
                e_i[a3_local] = 1.0
                comm_vec = expm_multiply(A.tocsc(), e_i)
            except Exception as e:
                log(f"    WARNING: sparse expm_multiply failed for {short}: {e}")
                continue

        # Extract communicability to Harris interactors
        for h in harris_in_comp:
            h_local = node_to_local[h]
            comm_val = comm_vec[h_local]
            top_pairs.append((a3, h, float(comm_val)))

            if a3 == "APOBEC3A":
                results['a3a_sum'] += comm_val
                results['a3a_max'] = max(results['a3a_max'], comm_val)
            elif a3 == "APOBEC3B":
                results['a3b_sum'] += comm_val
                results['a3b_max'] = max(results['a3b_max'], comm_val)
            results['total_sum'] += comm_val

    results['harris_reached'] = len(set(p[1] for p in top_pairs if p[2] > 1e-6))
    top_pairs.sort(key=lambda x: -x[2])
    results['top_pairs'] = top_pairs[:15]

    return results


# =============================================================================
# METRIC 4: K-CORE MEMBERSHIP
# =============================================================================

def compute_kcore(G, a3_genes_in_net, harris_in_net):
    """
    Compute k-core decomposition and report coreness of A3 and Harris genes.
    """
    if G.number_of_nodes() == 0:
        return {'a3_coreness': {}, 'harris_max_coreness': 0,
                'max_coreness': 0, 'a3_in_core_with_harris': False}

    coreness = nx.core_number(G)
    max_k = max(coreness.values()) if coreness else 0

    a3_coreness = {}
    for a3 in a3_genes_in_net:
        if a3 in coreness:
            a3_coreness[A3_SHORT.get(a3, a3)] = coreness[a3]

    harris_coreness = [coreness.get(h, 0) for h in harris_in_net if h in coreness]
    harris_max = max(harris_coreness) if harris_coreness else 0
    harris_mean = np.mean(harris_coreness) if harris_coreness else 0

    # Check if any A3 gene shares a k-core level with Harris interactors
    a3_core_vals = set(a3_coreness.values())
    harris_core_vals = set(coreness.get(h, 0) for h in harris_in_net if h in G)
    shared_cores = a3_core_vals & harris_core_vals - {0}

    return {
        'a3_coreness': a3_coreness,
        'harris_max_coreness': harris_max,
        'harris_mean_coreness': harris_mean,
        'max_coreness': max_k,
        'shared_core_levels': shared_cores,
        'a3_in_core_with_harris': len(shared_cores) > 0,
    }


# =============================================================================
# MAIN SWEEP
# =============================================================================

def sweep_network(net_name, net_info, corr_diff, harris_all, harris_a3b):
    """Run all four metrics across all thresholds for one network."""
    banner(f"FEATURE SWEEP: {net_name}", char="*")
    log(f"  {net_info['label']}")
    log(f"  Auto-selected threshold: {net_info['auto_threshold']}")
    log(f"  DIFF matrix: {corr_diff.shape}")

    genes_in_matrix = set(corr_diff.index)
    harris_present = harris_all & genes_in_matrix
    harris_a3b_present = harris_a3b & genes_in_matrix
    a3_present = [a3 for a3 in A3_GENES if a3 in genes_in_matrix]
    focal_present = [a3 for a3 in FOCAL_A3 if a3 in genes_in_matrix]

    log(f"  A3 genes in matrix: {len(a3_present)} "
        f"({', '.join(A3_SHORT[a] for a in a3_present)})")
    log(f"  Harris interactors in matrix: {len(harris_present)}")

    summary_rows = []
    auto_thresh = net_info['auto_threshold']

    for t in SWEEP_THRESHOLDS:
        is_auto = (t == auto_thresh)
        marker = " <<<AUTO" if is_auto else ""

        banner(f"Threshold {t:.2f}{marker}", char="-")

        # Build graph
        G, gene_to_idx, genes = build_graph_at_threshold(corr_diff, t)
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()

        if n_nodes == 0:
            log(f"  Empty graph at threshold {t}")
            summary_rows.append({
                'network': net_name, 'threshold': t, 'is_auto': is_auto,
                'nodes': 0, 'edges': 0, 'components': 0, 'lcc_size': 0,
                'tom_max': 0, 'tom_a3a_max': 0, 'tom_a3b_max': 0,
                'tom_nonzero_pairs': 0,
                'reach_a3a_harris': 0, 'reach_a3b_harris': 0,
                'reach_total_harris': 0,
                'comm_a3a_sum': 0, 'comm_a3b_sum': 0,
                'comm_harris_reached': 0,
                'kcore_max': 0, 'kcore_a3a': 0, 'kcore_a3b': 0,
                'kcore_harris_max': 0,
            })
            continue

        components = list(nx.connected_components(G))
        lcc = max(components, key=len) if components else set()
        n_comp = len(components)
        n_lcc = len(lcc)

        # Genes present in network at this threshold
        a3_in_net = [a3 for a3 in a3_present if a3 in G]
        harris_in_net = harris_present & set(G.nodes())

        log(f"  Network: {n_nodes} nodes, {n_edges} edges, "
            f"{n_comp} components, LCC={n_lcc}")
        log(f"  A3 in network: {len(a3_in_net)} "
            f"({', '.join(A3_SHORT[a] for a in a3_in_net)})")
        log(f"  Harris in network: {len(harris_in_net)}")

        # ---- METRIC 1: Topological Overlap ----
        log(f"\n  [TOM] Topological Overlap:")
        tom = compute_tom_pairs(G, a3_in_net, harris_in_net)
        log(f"    Pairs evaluated: {tom['n_pairs']}, non-zero: {tom.get('n_nonzero', 0)}")
        log(f"    Max TOM: {tom['max_tom']:.4f}")
        log(f"    A3A max TOM: {tom['a3a_max']:.4f}, "
            f"A3B max TOM: {tom['a3b_max']:.4f}")
        if tom['top_pairs']:
            log(f"    Top A3-Harris TOM pairs:")
            log(f"    {'A3':>8s}  {'Harris':>15s}  {'TOM':>8s}  "
                f"{'Shared':>7s}  {'Direct':>7s}")
            for a3, h, tv, shared, adj in tom['top_pairs'][:10]:
                log(f"    {A3_SHORT.get(a3,a3):>8s}  {h:>15s}  {tv:>8.4f}  "
                    f"{shared:>7d}  {'yes' if adj else 'no':>7s}")

        # ---- METRIC 2: 2-Hop Reachability ----
        log(f"\n  [REACH] 2-Hop Reachability:")
        reach, total_harris_reached = compute_reachability(
            G, a3_in_net, harris_in_net, harris_a3b_present)
        for short, rdata in sorted(reach.items()):
            log(f"    {short}: deg={rdata['degree']}, "
                f"1-hop={rdata['harris_1hop']}, "
                f"2-hop={rdata['harris_2hop']}, "
                f"total Harris={rdata['harris_total']} "
                f"({rdata['harris_a3b_reached']} A3B-sp)")
            if rdata['harris_1hop_names']:
                log(f"      1-hop Harris: {', '.join(rdata['harris_1hop_names'][:8])}"
                    f"{'...' if len(rdata['harris_1hop_names']) > 8 else ''}")
        log(f"    Total unique Harris reached: {total_harris_reached}")

        # ---- METRIC 3: Communicability ----
        log(f"\n  [COMM] Communicability:")
        comm = compute_communicability(G, a3_in_net, harris_in_net)
        log(f"    A3A sum to Harris: {comm['a3a_sum']:.4f}, "
            f"max: {comm['a3a_max']:.4f}")
        log(f"    A3B sum to Harris: {comm['a3b_sum']:.4f}, "
            f"max: {comm['a3b_max']:.4f}")
        log(f"    Harris reached (comm > 1e-6): {comm['harris_reached']}")
        if comm['top_pairs']:
            log(f"    Top communicability pairs:")
            for a3, h, cv in comm['top_pairs'][:8]:
                log(f"      {A3_SHORT.get(a3,a3):>5s} -> {h:>15s}  comm={cv:.6f}")

        # ---- METRIC 4: K-Core ----
        log(f"\n  [KCORE] K-Core Decomposition:")
        kcore = compute_kcore(G, a3_in_net, harris_in_net)
        log(f"    Max coreness in network: {kcore['max_coreness']}")
        if kcore['a3_coreness']:
            core_str = ", ".join(f"{k}={v}" for k, v in
                                sorted(kcore['a3_coreness'].items()))
            log(f"    A3 coreness: {core_str}")
        log(f"    Harris max coreness: {kcore['harris_max_coreness']}")
        log(f"    Harris mean coreness: {kcore['harris_mean_coreness']:.1f}")
        if kcore['shared_core_levels']:
            log(f"    A3 shares core levels with Harris: "
                f"{kcore['shared_core_levels']}")

        # ---- Collect row ----
        reach_a3a = reach.get('A3A', {})
        reach_a3b = reach.get('A3B', {})

        summary_rows.append({
            'network': net_name,
            'threshold': t,
            'is_auto': is_auto,
            'nodes': n_nodes,
            'edges': n_edges,
            'components': n_comp,
            'lcc_size': n_lcc,
            # TOM
            'tom_max': tom['max_tom'],
            'tom_a3a_max': tom['a3a_max'],
            'tom_a3b_max': tom['a3b_max'],
            'tom_nonzero_pairs': tom.get('n_nonzero', 0),
            # Reachability
            'reach_a3a_harris': reach_a3a.get('harris_total', 0),
            'reach_a3b_harris': reach_a3b.get('harris_total', 0),
            'reach_total_harris': total_harris_reached,
            'reach_a3a_1hop': reach_a3a.get('harris_1hop', 0),
            'reach_a3b_1hop': reach_a3b.get('harris_1hop', 0),
            'reach_a3a_2hop': reach_a3a.get('harris_2hop', 0),
            'reach_a3b_2hop': reach_a3b.get('harris_2hop', 0),
            # Communicability
            'comm_a3a_sum': comm['a3a_sum'],
            'comm_a3b_sum': comm['a3b_sum'],
            'comm_a3a_max': comm['a3a_max'],
            'comm_a3b_max': comm['a3b_max'],
            'comm_harris_reached': comm['harris_reached'],
            # K-core
            'kcore_max': kcore['max_coreness'],
            'kcore_a3a': kcore['a3_coreness'].get('A3A', 0),
            'kcore_a3b': kcore['a3_coreness'].get('A3B', 0),
            'kcore_harris_max': kcore['harris_max_coreness'],
            'kcore_harris_mean': kcore['harris_mean_coreness'],
        })

    return summary_rows


# =============================================================================
# PLOTS
# =============================================================================

def generate_plots(all_rows):
    """Generate 4-row x 3-col diagnostic plot (metrics x networks)."""
    banner("GENERATING PLOTS")

    df = pd.DataFrame(all_rows)
    net_names = list(NETWORKS.keys())

    fig, axes = plt.subplots(4, 3, figsize=(24, 26))
    fig.suptitle("A3 Pathway Network Feature Sweep\n"
                 "Four Metrics Across Three Networks",
                 fontsize=16, fontweight='bold')

    metric_labels = [
        "Topological Overlap\n(max A3-Harris TOM)",
        "2-Hop Reachability\n(Harris within 2 hops)",
        "Communicability\n(sum A3->Harris)",
        "K-Core\n(A3 coreness)"
    ]

    for col, net_name in enumerate(net_names):
        ndf = df[df['network'] == net_name].sort_values('threshold')
        auto_t = NETWORKS[net_name]['auto_threshold']

        if ndf.empty:
            continue

        threshs = ndf['threshold'].values

        # Row 0: Topological Overlap
        ax = axes[0, col]
        ax.plot(threshs, ndf['tom_a3a_max'], 'o-', color='#ed6a5a',
                label='A3A max TOM', markersize=5)
        ax.plot(threshs, ndf['tom_a3b_max'], 's-', color='#2c3e50',
                label='A3B max TOM', markersize=5)
        ax.axvline(auto_t, color='red', linestyle='--', alpha=0.5,
                   label=f'Auto={auto_t}')
        ax.set_ylabel(metric_labels[0])
        ax.set_title(f"{net_name}\n({NETWORKS[net_name]['label']})", fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # Row 1: 2-Hop Reachability
        ax = axes[1, col]
        ax.plot(threshs, ndf['reach_a3a_harris'], 'o-', color='#ed6a5a',
                label='A3A', markersize=5)
        ax.plot(threshs, ndf['reach_a3b_harris'], 's-', color='#2c3e50',
                label='A3B', markersize=5)
        ax.plot(threshs, ndf['reach_total_harris'], '^-', color='#f18f01',
                label='Any A3', markersize=5, alpha=0.6)
        ax.axvline(auto_t, color='red', linestyle='--', alpha=0.5)
        ax.set_ylabel(metric_labels[1])
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # Row 2: Communicability
        ax = axes[2, col]
        # Use log scale if values span orders of magnitude
        a3a_comm = ndf['comm_a3a_sum'].values
        a3b_comm = ndf['comm_a3b_sum'].values
        ax.plot(threshs, a3a_comm, 'o-', color='#ed6a5a',
                label='A3A sum', markersize=5)
        ax.plot(threshs, a3b_comm, 's-', color='#2c3e50',
                label='A3B sum', markersize=5)
        ax.axvline(auto_t, color='red', linestyle='--', alpha=0.5)
        ax.set_ylabel(metric_labels[2])
        ax.set_yscale('symlog', linthresh=0.01)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # Row 3: K-Core
        ax = axes[3, col]
        ax.plot(threshs, ndf['kcore_a3a'], 'o-', color='#ed6a5a',
                label='A3A', markersize=5)
        ax.plot(threshs, ndf['kcore_a3b'], 's-', color='#2c3e50',
                label='A3B', markersize=5)
        ax.plot(threshs, ndf['kcore_max'], 'x--', color='gray',
                label='Network max', markersize=4, alpha=0.5)
        ax.plot(threshs, ndf['kcore_harris_mean'], 'd-', color='#f18f01',
                label='Harris mean', markersize=4, alpha=0.6)
        ax.axvline(auto_t, color='red', linestyle='--', alpha=0.5)
        ax.set_ylabel(metric_labels[3])
        ax.set_xlabel("|delta-rho| threshold")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plot_path = os.path.join(OUTPUT_DIR, "A3_feature_sweep_plots.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Plots -> {plot_path}")

    # ---- Additional: threshold recommendation plot ----
    fig2, axes2 = plt.subplots(1, 3, figsize=(24, 7))
    fig2.suptitle("Threshold Selection: Composite A3 Pathway Signal",
                  fontsize=14, fontweight='bold')

    for col, net_name in enumerate(net_names):
        ndf = df[df['network'] == net_name].sort_values('threshold')
        auto_t = NETWORKS[net_name]['auto_threshold']

        if ndf.empty:
            continue

        threshs = ndf['threshold'].values

        # Normalize each metric to 0-1 for comparison
        ax = axes2[col]

        for label, col_name, color, marker in [
            ('TOM (A3B)', 'tom_a3b_max', '#2c3e50', 's'),
            ('2-hop reach', 'reach_total_harris', '#f18f01', '^'),
            ('Comm (A3B)', 'comm_a3b_sum', '#27ae60', 'd'),
            ('K-core (A3B)', 'kcore_a3b', '#8e44ad', 'v'),
        ]:
            vals = ndf[col_name].values.astype(float)
            vmax = vals.max()
            if vmax > 0:
                norm_vals = vals / vmax
            else:
                norm_vals = vals
            ax.plot(threshs, norm_vals, f'{marker}-', color=color,
                    label=label, markersize=5)

        # Also show network size (normalized)
        node_vals = ndf['nodes'].values.astype(float)
        node_max = node_vals.max()
        if node_max > 0:
            ax.plot(threshs, node_vals / node_max, 'x--', color='gray',
                    label='Net size', markersize=4, alpha=0.4)

        ax.axvline(auto_t, color='red', linestyle='--', alpha=0.7,
                   label=f'Auto={auto_t}')
        ax.set_xlabel("|delta-rho| threshold")
        ax.set_ylabel("Normalized metric value")
        ax.set_title(f"{net_name}")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    plot2_path = os.path.join(OUTPUT_DIR, "A3_feature_sweep_composite.png")
    plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Composite plot -> {plot2_path}")


# =============================================================================
# CROSS-NETWORK SUMMARY
# =============================================================================

def print_summary(all_rows):
    """Print a concise cross-network summary table."""
    banner("CROSS-NETWORK SUMMARY")

    df = pd.DataFrame(all_rows)

    for net_name in NETWORKS:
        ndf = df[df['network'] == net_name].sort_values('threshold')
        auto_t = NETWORKS[net_name]['auto_threshold']

        log(f"\n  {net_name} (auto={auto_t}):")
        log(f"  {'Thresh':>7s}  {'Nodes':>6s}  {'LCC':>5s}  "
            f"{'TOM_B':>6s}  {'Rch_A':>6s}  {'Rch_B':>6s}  "
            f"{'Com_A':>8s}  {'Com_B':>8s}  {'Kc_A':>5s}  {'Kc_B':>5s}  "
            f"{'HarNet':>7s}")
        log(f"  {'-----':>7s}  {'-----':>6s}  {'---':>5s}  "
            f"{'-----':>6s}  {'-----':>6s}  {'-----':>6s}  "
            f"{'------':>8s}  {'------':>8s}  {'----':>5s}  {'----':>5s}  "
            f"{'------':>7s}")

        for _, row in ndf.iterrows():
            marker = " <<<" if row['is_auto'] else ""
            log(f"  {row['threshold']:>7.2f}  {row['nodes']:>6.0f}  "
                f"{row['lcc_size']:>5.0f}  "
                f"{row['tom_a3b_max']:>6.3f}  "
                f"{row['reach_a3a_harris']:>6.0f}  "
                f"{row['reach_a3b_harris']:>6.0f}  "
                f"{row['comm_a3a_sum']:>8.3f}  "
                f"{row['comm_a3b_sum']:>8.3f}  "
                f"{row['kcore_a3a']:>5.0f}  "
                f"{row['kcore_a3b']:>5.0f}  "
                f"{row.get('reach_total_harris', 0):>7.0f}{marker}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("A3 PATHWAY NETWORK FEATURE SWEEP")
    log(f"  Start: {t0}")

    harris_all, harris_a3b = load_harris()
    log(f"  Harris (all): {len(harris_all)}, (A3B-sp): {len(harris_a3b)}")

    all_rows = []

    for net_name, net_info in NETWORKS.items():
        corr_diff = load_diff_matrix(net_info['dir'])
        if corr_diff is None:
            continue
        rows = sweep_network(net_name, net_info, corr_diff, harris_all, harris_a3b)
        all_rows.extend(rows)

    # Summary table
    print_summary(all_rows)

    # Save CSV
    summary_df = pd.DataFrame(all_rows)
    csv_path = os.path.join(OUTPUT_DIR, "A3_feature_sweep_summary.csv")
    summary_df.to_csv(csv_path, index=False)
    log(f"\n  [SAVE] Summary CSV -> {csv_path}")

    # Plots
    generate_plots(all_rows)

    # Save report
    report_path = os.path.join(OUTPUT_DIR, "A3_feature_sweep_report.txt")
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    log(f"  [SAVE] Report -> {report_path}")

    banner("FEATURE SWEEP COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")


if __name__ == "__main__":
    main()
