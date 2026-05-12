#!/usr/bin/env python3
"""
Diagnostic_A3_Correlation_Landscape.py
========================================

Examines the DIFF correlation landscape for A3 genes and Harris A3
interactors across all three Figure 4 single-cell networks.

Purpose: Understand why A3 genes and interactors are excluded from the
auto-selected network thresholds, and identify what threshold (if any)
would retain them.

Four diagnostic phases:

  Phase 1 — A3 gene correlation profiles
    For each A3 gene present in each network's DIFF matrix, report the
    full distribution of |delta-rho| values (max, 95th, 90th, 75th, median).
    Shows exactly where A3 genes sit relative to auto-selected thresholds.

  Phase 2 — Harris interactor correlation profiles
    Same analysis aggregated across the 174 Harris interactors. Reports
    what fraction have at least one edge above each threshold. Also
    examines A3-to-Harris and Harris-to-Harris DIFF correlations.

  Phase 3 — Threshold entry/exit tracking
    For A3A, A3B, and Harris interactors, sweep from 0.30 to 0.90 and
    report: how many are connected, how many are in the LCC, and their
    degree distribution at each threshold.

  Phase 4 — A3-centric neighborhood analysis
    At each threshold, what does the local 1-hop and 2-hop neighborhood
    of A3A/A3B look like? Are Harris interactors reachable within 2 hops?

Output:
  data/FIG_4/DIAGNOSTIC_AUDIT/A3_correlation_landscape_report.txt
  data/FIG_4/DIAGNOSTIC_AUDIT/A3_correlation_landscape_plots.png
  data/FIG_4/DIAGNOSTIC_AUDIT/A3_threshold_entry_tracking.csv

Usage:
  conda run -n NETWORK python Diagnostic_A3_Correlation_Landscape.py

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

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_4")

# Three network output directories
NETWORKS = {
    "SBS2_VS_CNV": {
        "label": "SBS2-HIGH vs CNV-HIGH (divergent fates)",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_CNV"),
        "auto_threshold": 0.65,
    },
    "SBS2_VS_NORMAL": {
        "label": "SBS2-HIGH vs NORMAL (mutagenic program entry)",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL"),
        "auto_threshold": 0.65,
    },
    "CNV_VS_NORMAL": {
        "label": "CNV-HIGH vs NORMAL (productive infection entry)",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL"),
        "auto_threshold": 0.70,
    },
}

# A3 genes (symbol format, matching SC expression matrices)
A3_GENES = [
    "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
    "APOBEC3F", "APOBEC3G", "APOBEC3H",
]
A3_SHORT = {
    "APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
    "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
    "APOBEC3H": "A3H",
}

# Harris interactor files (TSV with gene_symbol column)
HARRIS_ALL_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors_A3B_only.txt")

# Threshold sweep values (matching pipeline)
SWEEP_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                    0.75, 0.80, 0.85, 0.90]

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
# LOAD HARRIS INTERACTORS
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


# =============================================================================
# LOAD DIFF CORRELATION MATRIX FOR A NETWORK
# =============================================================================

def load_diff_matrix(net_dir):
    """Load the DIFF correlation matrix pickle from a network directory."""
    corr_path = os.path.join(net_dir, "03_correlation_networks",
                             "corr_matrices", "SC_corr_DIFF.pkl")
    if not os.path.exists(corr_path):
        log(f"  WARNING: DIFF matrix not found: {corr_path}")
        return None

    with open(corr_path, "rb") as f:
        corr_diff = pickle.load(f)
    return corr_diff


# =============================================================================
# PHASE 1: A3 GENE CORRELATION PROFILES
# =============================================================================

def phase1_a3_profiles(corr_diff, net_name, auto_thresh):
    """Report DIFF correlation distribution for each A3 gene."""
    banner(f"PHASE 1: A3 Gene Correlation Profiles — {net_name}")

    genes_in_matrix = set(corr_diff.index)
    log(f"  Matrix genes: {len(genes_in_matrix)}")

    results = []

    for a3 in A3_GENES:
        short = A3_SHORT[a3]
        if a3 not in genes_in_matrix:
            log(f"\n  {short}: NOT in DIFF matrix")
            continue

        # Get absolute DIFF correlations for this gene (exclude self)
        row = corr_diff.loc[a3].drop(a3, errors='ignore')
        abs_row = row.abs()

        # Key percentiles
        pmax = abs_row.max()
        p95 = abs_row.quantile(0.95)
        p90 = abs_row.quantile(0.90)
        p75 = abs_row.quantile(0.75)
        p50 = abs_row.median()

        # Degree at various thresholds
        degrees = {}
        for t in SWEEP_THRESHOLDS:
            degrees[t] = int((abs_row >= t).sum())

        # Top 10 strongest DIFF partners
        top10 = abs_row.nlargest(10)

        log(f"\n  {short} ({a3}):")
        log(f"    Max |delta-rho|: {pmax:.4f}")
        log(f"    95th percentile: {p95:.4f}")
        log(f"    90th percentile: {p90:.4f}")
        log(f"    75th percentile: {p75:.4f}")
        log(f"    Median:          {p50:.4f}")
        log(f"    Auto threshold:  {auto_thresh}")
        log(f"    {'INCLUDED' if pmax >= auto_thresh else 'EXCLUDED'} at auto threshold")
        log(f"")
        log(f"    Degree at each threshold:")
        log(f"    {'Thresh':>8s}  {'Degree':>8s}  {'Status':>10s}")
        log(f"    {'------':>8s}  {'------':>8s}  {'------':>10s}")
        for t in SWEEP_THRESHOLDS:
            status = "AUTO >>>" if t == auto_thresh else ""
            log(f"    {t:>8.2f}  {degrees[t]:>8d}  {status}")

        log(f"")
        log(f"    Top 10 DIFF partners:")
        log(f"    {'Gene':>20s}  {'|delta-rho|':>12s}  {'raw delta':>10s}")
        log(f"    {'-'*20}  {'-'*12}  {'-'*10}")
        for gene, val in top10.items():
            raw = corr_diff.loc[a3, gene]
            log(f"    {gene:>20s}  {val:>12.4f}  {raw:>+10.4f}")

        results.append({
            'gene': a3, 'short': short, 'max': pmax, 'p95': p95,
            'p90': p90, 'p75': p75, 'median': p50,
            'deg_at_auto': degrees.get(auto_thresh, 0),
            'degrees': degrees,
        })

    return results


# =============================================================================
# PHASE 2: HARRIS INTERACTOR CORRELATION PROFILES
# =============================================================================

def phase2_harris_profiles(corr_diff, net_name, auto_thresh,
                           harris_all, harris_a3b):
    """Report DIFF correlation landscape for Harris interactors."""
    banner(f"PHASE 2: Harris Interactor Landscape — {net_name}")

    genes_in_matrix = set(corr_diff.index)
    harris_present = harris_all & genes_in_matrix
    harris_a3b_present = harris_a3b & genes_in_matrix

    log(f"  Harris interactors in DIFF matrix: {len(harris_present)} / {len(harris_all)}")
    log(f"  Harris A3B-specific in DIFF matrix: {len(harris_a3b_present)} / {len(harris_a3b)}")

    if not harris_present:
        log(f"  No Harris interactors found in matrix")
        return {}

    # For each threshold: how many Harris interactors have degree >= 1?
    log(f"\n  Harris interactor connectivity across thresholds:")
    log(f"  {'Thresh':>8s}  {'Connected':>10s}  {'Frac':>8s}  "
        f"{'A3B-sp':>8s}  {'MedDeg':>8s}  {'MaxDeg':>8s}")
    log(f"  {'------':>8s}  {'-'*10}  {'----':>8s}  "
        f"{'------':>8s}  {'------':>8s}  {'------':>8s}")

    threshold_data = {}

    for t in SWEEP_THRESHOLDS:
        connected = 0
        a3b_connected = 0
        degrees_at_t = []

        for h in harris_present:
            row = corr_diff.loc[h].drop(h, errors='ignore')
            deg = int((row.abs() >= t).sum())
            if deg > 0:
                connected += 1
                degrees_at_t.append(deg)
                if h in harris_a3b_present:
                    a3b_connected += 1

        frac = connected / len(harris_present) if harris_present else 0
        med_deg = int(np.median(degrees_at_t)) if degrees_at_t else 0
        max_deg = max(degrees_at_t) if degrees_at_t else 0

        marker = "  AUTO >>>" if t == auto_thresh else ""
        log(f"  {t:>8.2f}  {connected:>10d}  {frac:>8.1%}  "
            f"{a3b_connected:>8d}  {med_deg:>8d}  {max_deg:>8d}{marker}")

        threshold_data[t] = {
            'connected': connected, 'fraction': frac,
            'a3b_connected': a3b_connected,
            'median_degree': med_deg, 'max_degree': max_deg,
        }

    # A3-to-Harris correlations in the DIFF matrix
    log(f"\n  A3-to-Harris DIFF correlations:")
    a3_in_matrix = [a3 for a3 in A3_GENES if a3 in genes_in_matrix]

    for a3 in a3_in_matrix:
        short = A3_SHORT[a3]
        harris_in = sorted(harris_present)
        corr_vals = []
        for h in harris_in:
            if h != a3:
                corr_vals.append(corr_diff.loc[a3, h])

        if not corr_vals:
            continue

        abs_vals = np.abs(corr_vals)
        log(f"\n    {short} -> Harris ({len(corr_vals)} interactors in matrix):")
        log(f"      Max |delta-rho|: {abs_vals.max():.4f}")
        log(f"      95th pctl:       {np.percentile(abs_vals, 95):.4f}")
        log(f"      90th pctl:       {np.percentile(abs_vals, 90):.4f}")
        log(f"      75th pctl:       {np.percentile(abs_vals, 75):.4f}")
        log(f"      Median:          {np.median(abs_vals):.4f}")
        log(f"      Mean:            {np.mean(abs_vals):.4f}")

        # Top 5 Harris interactors by |DIFF correlation| with this A3
        top_idx = np.argsort(abs_vals)[::-1][:5]
        harris_sorted = [harris_in[i] for i in range(len(harris_in))
                         if harris_in[i] != a3]
        abs_sorted = abs_vals.copy()
        top_pairs = sorted(zip(harris_sorted, abs_vals), key=lambda x: -x[1])[:5]
        log(f"      Top 5 interactors:")
        for gene, val in top_pairs:
            raw = corr_diff.loc[a3, gene]
            a3b_flag = " [A3B-sp]" if gene in harris_a3b else ""
            log(f"        {gene:>20s}  |rho|={val:.4f}  raw={raw:+.4f}{a3b_flag}")

    # Harris-to-Harris: average pairwise DIFF correlation
    log(f"\n  Harris-to-Harris DIFF correlations:")
    harris_list = sorted(harris_present)
    if len(harris_list) >= 2:
        pairwise = []
        for i in range(len(harris_list)):
            for j in range(i+1, len(harris_list)):
                pairwise.append(abs(corr_diff.loc[harris_list[i], harris_list[j]]))
        pairwise = np.array(pairwise)
        log(f"    Pairs: {len(pairwise)}")
        log(f"    Max |delta-rho|: {pairwise.max():.4f}")
        log(f"    95th pctl:       {np.percentile(pairwise, 95):.4f}")
        log(f"    90th pctl:       {np.percentile(pairwise, 90):.4f}")
        log(f"    75th pctl:       {np.percentile(pairwise, 75):.4f}")
        log(f"    Median:          {np.median(pairwise):.4f}")
        log(f"    Mean:            {np.mean(pairwise):.4f}")

    return threshold_data


# =============================================================================
# PHASE 3: THRESHOLD ENTRY/EXIT TRACKING
# =============================================================================

def phase3_threshold_tracking(corr_diff, net_name, auto_thresh,
                              harris_all, harris_a3b):
    """Track A3 and Harris interactor inclusion in network and LCC."""
    banner(f"PHASE 3: Threshold Entry/Exit Tracking — {net_name}")

    genes_in_matrix = set(corr_diff.index)
    harris_present = harris_all & genes_in_matrix

    rows = []

    for t in SWEEP_THRESHOLDS:
        # Build graph at this threshold
        abs_diff = corr_diff.abs()
        np.fill_diagonal(abs_diff.values, 0)

        # Find edges above threshold
        mask = abs_diff >= t
        edge_count = int(mask.sum().sum() / 2)

        G = nx.Graph()
        genes = list(corr_diff.index)
        G.add_nodes_from(genes)

        # Add edges (upper triangle only)
        idxs = np.where(np.triu(mask.values, k=1))
        for i, j in zip(idxs[0], idxs[1]):
            G.add_edge(genes[i], genes[j],
                       weight=float(corr_diff.iloc[i, j]))

        # Remove isolates
        G_noiso = G.copy()
        isolates = list(nx.isolates(G_noiso))
        G_noiso.remove_nodes_from(isolates)

        n_nodes = G_noiso.number_of_nodes()
        n_edges = G_noiso.number_of_edges()

        # LCC
        if n_nodes > 0:
            components = list(nx.connected_components(G_noiso))
            lcc = max(components, key=len) if components else set()
            n_lcc = len(lcc)
            n_comp = len(components)
        else:
            lcc = set()
            n_lcc = 0
            n_comp = 0

        # A3 gene status
        a3_in_net = []
        a3_in_lcc = []
        a3_degrees = {}
        for a3 in A3_GENES:
            if a3 in G_noiso.nodes():
                deg = G_noiso.degree(a3)
                a3_degrees[A3_SHORT[a3]] = deg
                a3_in_net.append(A3_SHORT[a3])
                if a3 in lcc:
                    a3_in_lcc.append(A3_SHORT[a3])

        # Harris interactor status
        harris_in_net = sum(1 for h in harris_present if h in G_noiso.nodes())
        harris_in_lcc = sum(1 for h in harris_present if h in lcc)

        marker = " <<<AUTO" if t == auto_thresh else ""
        log(f"\n  Threshold {t:.2f}:{marker}")
        log(f"    Network: {n_nodes} nodes, {n_edges} edges, {n_comp} components, LCC={n_lcc}")
        log(f"    A3 in network: {len(a3_in_net)} ({', '.join(a3_in_net) if a3_in_net else 'none'})")
        log(f"    A3 in LCC:     {len(a3_in_lcc)} ({', '.join(a3_in_lcc) if a3_in_lcc else 'none'})")
        if a3_degrees:
            deg_str = ", ".join(f"{k}={v}" for k, v in a3_degrees.items())
            log(f"    A3 degrees:    {deg_str}")
        log(f"    Harris in network: {harris_in_net} / {len(harris_present)}")
        log(f"    Harris in LCC:     {harris_in_lcc} / {len(harris_present)}")

        rows.append({
            'network': net_name, 'threshold': t,
            'nodes': n_nodes, 'edges': n_edges,
            'components': n_comp, 'lcc_size': n_lcc,
            'a3_in_net': len(a3_in_net), 'a3_in_lcc': len(a3_in_lcc),
            'a3_names': ';'.join(a3_in_net),
            'harris_in_net': harris_in_net, 'harris_in_lcc': harris_in_lcc,
            'a3_degrees': str(a3_degrees),
        })

    return rows


# =============================================================================
# PHASE 4: A3-CENTRIC NEIGHBORHOOD ANALYSIS
# =============================================================================

def phase4_neighborhood(corr_diff, net_name, auto_thresh,
                        harris_all, harris_a3b):
    """Analyze 1-hop and 2-hop neighborhoods of A3A/A3B at key thresholds."""
    banner(f"PHASE 4: A3-Centric Neighborhood — {net_name}")

    genes_in_matrix = set(corr_diff.index)
    harris_present = harris_all & genes_in_matrix

    # Focus on A3A and A3B
    focal_a3 = [a3 for a3 in ["APOBEC3A", "APOBEC3B"] if a3 in genes_in_matrix]
    if not focal_a3:
        log(f"  Neither A3A nor A3B in DIFF matrix for {net_name}")
        return

    # Examine at a few key thresholds: 0.40, 0.50, auto, and one below auto
    key_thresholds = sorted(set([0.35, 0.40, 0.45, 0.50, 0.55,
                                  auto_thresh - 0.05, auto_thresh]))

    for a3 in focal_a3:
        short = A3_SHORT[a3]
        log(f"\n  {short} neighborhood analysis:")

        row = corr_diff.loc[a3].drop(a3, errors='ignore')
        abs_row = row.abs()

        for t in key_thresholds:
            # 1-hop neighbors at this threshold
            neighbors_1 = set(abs_row[abs_row >= t].index)
            n1 = len(neighbors_1)

            # Classify 1-hop neighbors
            harris_1hop = neighbors_1 & harris_present
            a3_1hop = neighbors_1 & set(A3_GENES)

            # 2-hop neighbors: for each 1-hop neighbor, find THEIR neighbors
            neighbors_2 = set()
            harris_2hop_only = set()  # Harris genes reachable in exactly 2 hops
            for n1_gene in neighbors_1:
                n1_row = corr_diff.loc[n1_gene].drop(n1_gene, errors='ignore')
                n1_neighbors = set(n1_row.abs()[n1_row.abs() >= t].index)
                neighbors_2 |= n1_neighbors

            # Remove 1-hop and self from 2-hop
            neighbors_2 -= neighbors_1
            neighbors_2 -= {a3}
            harris_2hop = neighbors_2 & harris_present

            marker = " <<<AUTO" if t == auto_thresh else ""
            log(f"\n    Threshold {t:.2f}:{marker}")
            log(f"      1-hop: {n1} genes ({len(harris_1hop)} Harris, "
                f"{len(a3_1hop)} A3)")
            log(f"      2-hop: {len(neighbors_2)} additional genes "
                f"({len(harris_2hop)} Harris)")

            if harris_1hop:
                log(f"      Harris 1-hop: {', '.join(sorted(harris_1hop)[:10])}"
                    f"{'...' if len(harris_1hop) > 10 else ''}")
            if harris_2hop:
                log(f"      Harris 2-hop: {', '.join(sorted(harris_2hop)[:10])}"
                    f"{'...' if len(harris_2hop) > 10 else ''}")

            # Show top 5 1-hop neighbors with their correlation
            if neighbors_1:
                top5 = abs_row[list(neighbors_1)].nlargest(min(5, len(neighbors_1)))
                log(f"      Top 1-hop partners:")
                for gene, val in top5.items():
                    flags = []
                    if gene in harris_present:
                        flags.append("Harris")
                    if gene in harris_a3b:
                        flags.append("A3B-sp")
                    if gene in set(A3_GENES):
                        flags.append("A3")
                    flag_str = f"  [{', '.join(flags)}]" if flags else ""
                    log(f"        {gene:>20s}  |rho|={val:.4f}{flag_str}")


# =============================================================================
# GENERATE PLOTS
# =============================================================================

def generate_plots(all_a3_results, all_tracking, harris_threshold_data):
    """Generate summary plots across all three networks."""
    banner("GENERATING DIAGNOSTIC PLOTS")

    fig, axes = plt.subplots(3, 3, figsize=(24, 20))
    fig.suptitle("A3 & Harris Interactor Correlation Landscape\n"
                 "Three-Network Diagnostic", fontsize=16, fontweight='bold')

    net_names = list(NETWORKS.keys())

    for col, net_name in enumerate(net_names):
        net_info = NETWORKS[net_name]
        auto_thresh = net_info['auto_threshold']

        # ---- Row 0: A3 gene max |delta-rho| across thresholds ----
        ax = axes[0, col]
        ax.set_title(f"{net_name}\n({net_info['label'].split('(')[0].strip()})",
                     fontsize=11)

        if net_name in all_a3_results:
            for a3_data in all_a3_results[net_name]:
                threshs = sorted(a3_data['degrees'].keys())
                degs = [a3_data['degrees'][t] for t in threshs]
                ax.plot(threshs, degs, 'o-', label=a3_data['short'], markersize=4)

        ax.axvline(auto_thresh, color='red', linestyle='--', alpha=0.7,
                   label=f'Auto={auto_thresh}')
        ax.set_xlabel("|delta-rho| threshold")
        ax.set_ylabel("A3 gene degree")
        ax.legend(fontsize=8)
        ax.set_yscale('symlog', linthresh=1)
        ax.grid(True, alpha=0.3)

        # ---- Row 1: Harris interactor connectivity ----
        ax = axes[1, col]
        if net_name in harris_threshold_data:
            hdata = harris_threshold_data[net_name]
            threshs = sorted(hdata.keys())
            connected = [hdata[t]['connected'] for t in threshs]
            ax.plot(threshs, connected, 's-', color='#f18f01',
                    label='Harris connected', markersize=5)

        ax.axvline(auto_thresh, color='red', linestyle='--', alpha=0.7,
                   label=f'Auto={auto_thresh}')
        ax.set_xlabel("|delta-rho| threshold")
        ax.set_ylabel("Harris interactors connected")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # ---- Row 2: Network size + A3 inclusion ----
        ax = axes[2, col]
        tracking = [r for r in all_tracking if r['network'] == net_name]
        if tracking:
            threshs = [r['threshold'] for r in tracking]
            nodes = [r['nodes'] for r in tracking]
            a3_in = [r['a3_in_net'] for r in tracking]

            ax2 = ax.twinx()
            ax.plot(threshs, nodes, 'o-', color='#2c3e50', label='Network nodes',
                    markersize=4)
            ax2.plot(threshs, a3_in, 's-', color='#ed6a5a', label='A3 in network',
                     markersize=5)
            ax2.set_ylabel("A3 genes in network", color='#ed6a5a')
            ax2.set_ylim(-0.5, 8)

        ax.axvline(auto_thresh, color='red', linestyle='--', alpha=0.7,
                   label=f'Auto={auto_thresh}')
        ax.set_xlabel("|delta-rho| threshold")
        ax.set_ylabel("Network nodes")
        ax.set_yscale('log')
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plot_path = os.path.join(OUTPUT_DIR, "A3_correlation_landscape_plots.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Plots -> {plot_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("DIAGNOSTIC: A3 Correlation Landscape Across Three Networks")
    log(f"  Start: {t0}")
    log(f"  Output: {OUTPUT_DIR}")

    # Load Harris interactors
    harris_all, harris_a3b = load_harris()
    log(f"  Harris interactors (all): {len(harris_all)}")
    log(f"  Harris interactors (A3B-specific): {len(harris_a3b)}")

    # Collect results across networks
    all_a3_results = {}
    all_tracking = []
    harris_threshold_data = {}

    for net_name, net_info in NETWORKS.items():
        net_dir = net_info['dir']
        auto_thresh = net_info['auto_threshold']

        banner(f"NETWORK: {net_name}", char="*")
        log(f"  {net_info['label']}")
        log(f"  Directory: {net_dir}")
        log(f"  Auto-selected threshold: {auto_thresh}")

        # Load DIFF correlation matrix
        corr_diff = load_diff_matrix(net_dir)
        if corr_diff is None:
            log(f"  SKIPPING (no DIFF matrix)")
            continue

        log(f"  DIFF matrix: {corr_diff.shape}")

        # Phase 1: A3 gene profiles
        a3_results = phase1_a3_profiles(corr_diff, net_name, auto_thresh)
        all_a3_results[net_name] = a3_results

        # Phase 2: Harris interactor profiles
        h_data = phase2_harris_profiles(corr_diff, net_name, auto_thresh,
                                         harris_all, harris_a3b)
        harris_threshold_data[net_name] = h_data

        # Phase 3: Threshold entry/exit tracking
        tracking = phase3_threshold_tracking(corr_diff, net_name, auto_thresh,
                                              harris_all, harris_a3b)
        all_tracking.extend(tracking)

        # Phase 4: A3-centric neighborhood
        phase4_neighborhood(corr_diff, net_name, auto_thresh,
                           harris_all, harris_a3b)

    # Generate plots
    generate_plots(all_a3_results, all_tracking, harris_threshold_data)

    # Save tracking CSV
    tracking_df = pd.DataFrame(all_tracking)
    tracking_path = os.path.join(OUTPUT_DIR, "A3_threshold_entry_tracking.csv")
    tracking_df.to_csv(tracking_path, index=False)
    log(f"\n  [SAVE] Tracking CSV -> {tracking_path}")

    # Save report
    report_path = os.path.join(OUTPUT_DIR, "A3_correlation_landscape_report.txt")
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    log(f"  [SAVE] Report -> {report_path}")

    banner("DIAGNOSTIC COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")


if __name__ == "__main__":
    main()
