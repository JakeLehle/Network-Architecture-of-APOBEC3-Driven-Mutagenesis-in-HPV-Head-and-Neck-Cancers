#!/usr/bin/env python3
"""
Auto_Select_Threshold_Tiered.py
==================================

Two-tier DIFF threshold selection for differential co-expression networks.

TIER 1 — Biology-informed (enrichment-based):
  Prerequisites: A3 seed genes AND known interactors (Harris et al.) both
  present in the DE gene list / DIFF correlation matrix.

  Method: At each candidate threshold, compute the 2-hop neighborhood of
  all A3 seed genes. Test whether Harris interactors are enriched in this
  neighborhood using the hypergeometric distribution. Select the threshold
  with the most significant enrichment (lowest BH-corrected p-value).

  Rationale: The threshold where the A3 neighborhood is most specifically
  enriched for known interactors is where the network best resolves the
  A3 pathway biology. At loose thresholds, the neighborhood is the whole
  network (no enrichment). At tight thresholds, the neighborhood is empty.
  The peak enrichment marks the sweet spot.

TIER 2 — Topology-driven (component peak + LCC):
  Fallback when Tier 1 prerequisites are not met (e.g., no Harris
  interactors survive the DE filter, as in the TCGA bulk network).

  Method: Find the threshold where connected components peak. If LCC at
  that threshold < max_lcc, use it. Otherwise step up to the next
  threshold where LCC < max_lcc.

  Rationale: Component peak = natural fragmentation boundaries.
  LCC < 300 = pathway-scale communities.

TIER SELECTION LOGIC:
  1. Count A3 genes in DIFF matrix
  2. Count Harris interactors in DIFF matrix
  3. If >= 1 A3 gene AND >= MIN_HARRIS interactors:
     a. Run Tier 1 enrichment sweep
     b. If best BH-corrected p < ENRICHMENT_ALPHA: use Tier 1 threshold
     c. If not significant: fall back to Tier 2
  4. If prerequisites not met: use Tier 2

OUTPUT:
  selected_parameters.txt     — Key=value file for pipeline consumption
  threshold_selection_report.txt  — Full diagnostic report
  threshold_selection_plots.png   — Enrichment curve + topology overlay

Can be run standalone or imported by Step03 for inline threshold selection.

Usage:
  # For a specific SC network:
  python Auto_Select_Threshold_Tiered.py \
      --corr_pkl path/to/SC_corr_DIFF.pkl \
      --harris path/to/Harris_A3_interactors.txt \
      --output_dir path/to/output/

  # For all three SC networks:
  python Auto_Select_Threshold_Tiered.py --mode three_networks

  # For TCGA bulk:
  python Auto_Select_Threshold_Tiered.py \
      --corr_pkl path/to/TCGA-HNSC_corr_DIFF.pkl \
      --harris path/to/Harris_A3_interactors.txt \
      --output_dir path/to/output/

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import hypergeom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_4")
FIG2_ROOT = os.path.join(PROJECT_ROOT, "data", "FIG_2")

# Three SC networks (default mode)
SC_NETWORKS = {
    "SBS2_VS_CNV": {
        "label": "SBS2-HIGH vs CNV-HIGH (divergent fates)",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_CNV"),
    },
    "SBS2_VS_NORMAL": {
        "label": "SBS2-HIGH vs NORMAL (mutagenic program entry)",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL"),
    },
    "CNV_VS_NORMAL": {
        "label": "CNV-HIGH vs NORMAL (productive infection entry)",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL"),
    },
}

# A3 seed genes (symbol format)
A3_GENES = [
    "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
    "APOBEC3F", "APOBEC3G", "APOBEC3H",
]
A3_SHORT = {
    "APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
    "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
    "APOBEC3H": "A3H",
}

# Harris interactor paths
HARRIS_ALL_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(FIG4_ROOT, "00_input", "Harris_A3_interactors_A3B_only.txt")

# Threshold sweep
SWEEP_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                    0.75, 0.80, 0.85, 0.90]

# Tier 1 parameters
MIN_HARRIS_FOR_TIER1 = 10     # Minimum Harris interactors in matrix for Tier 1
ENRICHMENT_ALPHA = 0.05       # BH-corrected p-value threshold for Tier 1

# Tier 2 parameters
TIER2_MAX_LCC = 300           # Maximum LCC size
TIER2_COMP_FRACTION = 0.80    # Minimum fraction of peak components

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

def load_harris(harris_all_path=HARRIS_ALL_PATH, harris_a3b_path=HARRIS_A3B_PATH):
    """Load Harris A3 interactors from TSV files (gene_symbol column)."""
    harris_all = set()
    harris_a3b = set()

    if os.path.exists(harris_all_path):
        try:
            df = pd.read_csv(harris_all_path, sep='\t')
            harris_all = set(df['gene_symbol'].dropna().values)
        except Exception:
            # Fallback: try plain text (one gene per line)
            try:
                with open(harris_all_path) as f:
                    harris_all = set(line.strip() for line in f
                                     if line.strip() and not line.startswith('#'))
            except Exception as e:
                log(f"  WARNING: Could not load Harris interactors: {e}")

    if os.path.exists(harris_a3b_path):
        try:
            df = pd.read_csv(harris_a3b_path, sep='\t')
            harris_a3b = set(df['gene_symbol'].dropna().values)
        except Exception:
            try:
                with open(harris_a3b_path) as f:
                    harris_a3b = set(line.strip() for line in f
                                     if line.strip() and not line.startswith('#'))
            except Exception:
                pass

    return harris_all, harris_a3b


def load_diff_matrix(corr_pkl_path):
    """Load DIFF correlation matrix from pickle."""
    if not os.path.exists(corr_pkl_path):
        log(f"  ERROR: DIFF matrix not found: {corr_pkl_path}")
        return None
    with open(corr_pkl_path, "rb") as f:
        return pickle.load(f)


def build_graph(corr_diff, threshold):
    """Build undirected graph from DIFF matrix at |delta-rho| >= threshold."""
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

    # Remove isolates
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    return G


def bh_correct(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.array(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    adjusted = pvals * n / ranked
    # Enforce monotonicity
    adjusted = np.minimum.accumulate(adjusted[np.argsort(ranked)[::-1]])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    # Return in original order
    result = np.empty(n)
    result[np.argsort(ranked)[::-1]] = adjusted
    return result


# =============================================================================
# TIER 1: ENRICHMENT-BASED THRESHOLD SELECTION
# =============================================================================

def tier1_enrichment_sweep(corr_diff, a3_present, harris_present, harris_a3b_present):
    """
    At each threshold, test whether Harris interactors are enriched in the
    2-hop neighborhood of A3 genes using the hypergeometric distribution.

    Returns a DataFrame with one row per threshold and the selected threshold.
    """
    banner("TIER 1: Enrichment-Based Threshold Selection")
    log(f"  A3 seed genes in matrix: {len(a3_present)} "
        f"({', '.join(A3_SHORT.get(a, a) for a in a3_present)})")
    log(f"  Harris interactors in matrix: {len(harris_present)}")
    log(f"  Harris A3B-specific in matrix: {len(harris_a3b_present)}")

    rows = []

    for t in SWEEP_THRESHOLDS:
        G = build_graph(corr_diff, t)
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        net_genes = set(G.nodes())

        # Components and LCC
        if n_nodes > 0:
            components = list(nx.connected_components(G))
            lcc = max(components, key=len)
            n_comp = len(components)
            n_lcc = len(lcc)
        else:
            components = []
            lcc = set()
            n_comp = 0
            n_lcc = 0

        # A3 genes in network
        a3_in_net = [a3 for a3 in a3_present if a3 in net_genes]
        harris_in_net = harris_present & net_genes
        harris_a3b_in_net = harris_a3b_present & net_genes

        # Compute 2-hop neighborhood of all A3 genes combined
        a3_neighborhood = set()
        for a3 in a3_in_net:
            # Include the A3 gene itself
            a3_neighborhood.add(a3)
            # 1-hop
            n1 = set(G.neighbors(a3))
            a3_neighborhood |= n1
            # 2-hop
            for neighbor in n1:
                a3_neighborhood |= set(G.neighbors(neighbor))

        # Remove A3 genes from neighborhood for cleaner test
        a3_set = set(a3_present)
        a3_neighborhood -= a3_set

        # Hypergeometric test parameters
        # N = total non-A3 genes in the network
        # K = Harris interactors among those (non-A3)
        # n = non-A3 genes in the 2-hop neighborhood
        # k = Harris interactors in the 2-hop neighborhood (non-A3)
        N = len(net_genes - a3_set)
        K = len(harris_in_net - a3_set)
        n_sample = len(a3_neighborhood)
        k_success = len(a3_neighborhood & (harris_present - a3_set))
        k_a3b_success = len(a3_neighborhood & (harris_a3b_present - a3_set))

        # Expected by chance
        expected = (K * n_sample / N) if N > 0 else 0

        # Hypergeometric p-value (P(X >= k))
        if N > 0 and K > 0 and n_sample > 0 and k_success > 0:
            pval = hypergeom.sf(k_success - 1, N, K, n_sample)
        else:
            pval = 1.0

        # Fold enrichment
        fold = (k_success / expected) if expected > 0 else 0

        # A3 gene details
        a3_degrees = {}
        a3_in_lcc = []
        for a3 in a3_in_net:
            a3_degrees[A3_SHORT.get(a3, a3)] = G.degree(a3)
            if a3 in lcc:
                a3_in_lcc.append(A3_SHORT.get(a3, a3))

        rows.append({
            'threshold': t,
            'nodes': n_nodes,
            'edges': n_edges,
            'components': n_comp,
            'lcc_size': n_lcc,
            'a3_in_net': len(a3_in_net),
            'a3_in_lcc': len(a3_in_lcc),
            'a3_names': ';'.join(A3_SHORT.get(a, a) for a in a3_in_net),
            'harris_in_net': len(harris_in_net),
            'neighborhood_size': n_sample,
            'harris_in_neighborhood': k_success,
            'harris_a3b_in_neighborhood': k_a3b_success,
            'expected': expected,
            'fold_enrichment': fold,
            'pval_raw': pval,
            'a3_degrees': str(a3_degrees),
            'a3_in_lcc_names': ';'.join(a3_in_lcc),
        })

    df = pd.DataFrame(rows)

    # BH correction across thresholds
    df['pval_bh'] = bh_correct(df['pval_raw'].values)
    df['neg_log10_pval'] = -np.log10(df['pval_raw'].replace(0, 1e-300))
    df['neg_log10_bh'] = -np.log10(df['pval_bh'].replace(0, 1e-300))

    # Print summary table
    log(f"\n  {'Thresh':>7s}  {'Nodes':>6s}  {'LCC':>5s}  {'A3net':>5s}  "
        f"{'A3lcc':>5s}  {'HarN':>5s}  {'Nhood':>6s}  {'HarNH':>5s}  "
        f"{'Expct':>6s}  {'Fold':>6s}  {'rawP':>10s}  {'BH_P':>10s}")
    log(f"  {'-----':>7s}  {'-----':>6s}  {'---':>5s}  {'-----':>5s}  "
        f"{'-----':>5s}  {'----':>5s}  {'-----':>6s}  {'-----':>5s}  "
        f"{'-----':>6s}  {'----':>6s}  {'--------':>10s}  {'--------':>10s}")

    for _, row in df.iterrows():
        praw = f"{row['pval_raw']:.2e}" if row['pval_raw'] < 0.01 else f"{row['pval_raw']:.4f}"
        pbh = f"{row['pval_bh']:.2e}" if row['pval_bh'] < 0.01 else f"{row['pval_bh']:.4f}"
        log(f"  {row['threshold']:>7.2f}  {row['nodes']:>6.0f}  {row['lcc_size']:>5.0f}  "
            f"{row['a3_in_net']:>5.0f}  {row['a3_in_lcc']:>5.0f}  "
            f"{row['harris_in_net']:>5.0f}  {row['neighborhood_size']:>6.0f}  "
            f"{row['harris_in_neighborhood']:>5.0f}  "
            f"{row['expected']:>6.1f}  {row['fold_enrichment']:>6.2f}  "
            f"{praw:>10s}  {pbh:>10s}")

    # ---- Selection: rank by fold enrichment, pick first passing raw p ----
    #
    # Rationale: Adjacent thresholds share most edges, so BH correction
    # for independent tests is too conservative. Instead, we rank thresholds
    # by fold enrichment (strongest biological signal first) and select the
    # first that passes a raw p < 0.05 significance cutoff. This asks:
    # "which threshold gives the strongest enrichment that isn't due to chance?"
    #
    # If no threshold passes, Tier 1 fails and we fall back to Tier 2.

    # Must have at least 1 A3 gene and at least 1 Harris interactor in neighborhood
    valid = df[(df['a3_in_net'] >= 1) & (df['harris_in_neighborhood'] >= 1)].copy()

    if len(valid) == 0:
        log(f"\n  Tier 1 FAILED: No threshold has both A3 genes and Harris "
            f"interactors in the 2-hop neighborhood.")
        return df, None, "No valid threshold found"

    # Sort by fold enrichment (descending), then by raw p-value (ascending)
    ranked = valid.sort_values(['fold_enrichment', 'pval_raw'],
                                ascending=[False, True])

    log(f"\n  Ranked by fold enrichment (descending):")
    log(f"  {'Thresh':>7s}  {'Fold':>6s}  {'rawP':>10s}  {'HarNH':>5s}  "
        f"{'Nhood':>6s}  {'A3net':>5s}  {'Pass':>5s}")
    log(f"  {'-----':>7s}  {'----':>6s}  {'--------':>10s}  {'-----':>5s}  "
        f"{'-----':>6s}  {'-----':>5s}  {'----':>5s}")

    for _, row in ranked.iterrows():
        praw = f"{row['pval_raw']:.2e}" if row['pval_raw'] < 0.01 else f"{row['pval_raw']:.4f}"
        passes = "YES" if row['pval_raw'] < ENRICHMENT_ALPHA else "no"
        log(f"  {row['threshold']:>7.2f}  {row['fold_enrichment']:>6.2f}  "
            f"{praw:>10s}  {row['harris_in_neighborhood']:>5.0f}  "
            f"{row['neighborhood_size']:>6.0f}  {row['a3_in_net']:>5.0f}  "
            f"{passes:>5s}")

    # Select: highest fold enrichment that passes raw p < alpha
    passing = ranked[ranked['pval_raw'] < ENRICHMENT_ALPHA]

    if len(passing) == 0:
        # Report the best fold even though it didn't pass
        best_row = ranked.iloc[0]
        log(f"\n  Tier 1 INCONCLUSIVE: No threshold passes raw p < {ENRICHMENT_ALPHA}")
        log(f"    Best fold: {best_row['fold_enrichment']:.2f}x at "
            f"threshold={best_row['threshold']:.2f} "
            f"(p={best_row['pval_raw']:.4f})")
        reason = (f"Best fold={best_row['fold_enrichment']:.2f}x at "
                  f"threshold={best_row['threshold']:.2f}, "
                  f"p={best_row['pval_raw']:.4f} > alpha={ENRICHMENT_ALPHA}")
        return df, None, reason

    # First passing row (highest fold that is significant)
    best_row = passing.iloc[0]
    best_thresh = float(best_row['threshold'])
    best_praw = float(best_row['pval_raw'])
    best_fold = float(best_row['fold_enrichment'])

    log(f"\n  Selected: threshold={best_thresh:.2f} "
        f"(fold={best_fold:.2f}x, p={best_praw:.2e})")
    log(f"    A3 in network: {int(best_row['a3_in_net'])} "
        f"({best_row['a3_names']})")
    log(f"    A3 in LCC: {int(best_row['a3_in_lcc'])} "
        f"({best_row['a3_in_lcc_names']})")
    log(f"    Harris in 2-hop neighborhood: "
        f"{int(best_row['harris_in_neighborhood'])} / "
        f"{int(best_row['harris_in_net'])} in network")
    log(f"    Neighborhood size: {int(best_row['neighborhood_size'])} genes")
    log(f"    Network: {int(best_row['nodes'])} nodes, "
        f"{int(best_row['edges'])} edges, "
        f"LCC={int(best_row['lcc_size'])}")

    reason = (f"Harris enrichment in A3 2-hop neighborhood: "
              f"fold={best_fold:.2f}x, p={best_praw:.2e}, "
              f"{int(best_row['harris_in_neighborhood'])} interactors")
    log(f"\n  Tier 1 SELECTED: threshold={best_thresh:.2f}")
    log(f"    Reason: {reason}")

    return df, best_thresh, reason


# =============================================================================
# TIER 2: TOPOLOGY-DRIVEN THRESHOLD SELECTION
# =============================================================================

def tier2_topology(corr_diff, max_lcc=TIER2_MAX_LCC,
                   comp_fraction=TIER2_COMP_FRACTION):
    """
    Select threshold using component peak + LCC constraint.
    Matches the existing auto_select_threshold logic from the pipeline.
    """
    banner("TIER 2: Topology-Driven Threshold Selection")
    log(f"  Max LCC: {max_lcc}")
    log(f"  Min component fraction: {comp_fraction}")

    rows = []
    for t in SWEEP_THRESHOLDS:
        G = build_graph(corr_diff, t)
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()

        if n_nodes > 0:
            components = list(nx.connected_components(G))
            n_comp = len(components)
            n_lcc = len(max(components, key=len))
        else:
            n_comp = 0
            n_lcc = 0

        rows.append({
            'threshold': t, 'nodes': n_nodes, 'edges': n_edges,
            'components': n_comp, 'lcc': n_lcc,
        })

    df = pd.DataFrame(rows)

    log(f"\n  {'Thresh':>7s}  {'Nodes':>6s}  {'Edges':>7s}  "
        f"{'Comp':>5s}  {'LCC':>5s}")
    log(f"  {'-----':>7s}  {'-----':>6s}  {'-----':>7s}  "
        f"{'----':>5s}  {'---':>5s}")
    for _, row in df.iterrows():
        log(f"  {row['threshold']:>7.2f}  {row['nodes']:>6.0f}  "
            f"{row['edges']:>7.0f}  {row['components']:>5.0f}  "
            f"{row['lcc']:>5.0f}")

    # Find component peak
    peak_idx = df['components'].idxmax()
    peak_row = df.loc[peak_idx]
    peak_thresh = float(peak_row['threshold'])
    peak_comp = int(peak_row['components'])
    peak_lcc = int(peak_row['lcc'])

    log(f"\n  Component peak: threshold={peak_thresh}, "
        f"components={peak_comp}, LCC={peak_lcc}")

    # Apply LCC constraint
    if peak_lcc <= max_lcc:
        selected = peak_row
        reason = f"Component peak (comp={peak_comp}) with LCC={peak_lcc} within limit"
    else:
        candidates = df[(df['threshold'] > peak_thresh) & (df['lcc'] <= max_lcc)]
        if len(candidates) > 0:
            selected = candidates.iloc[0]
            reason = (f"LCC at peak ({peak_lcc}) > {max_lcc}, "
                      f"stepped up from {peak_thresh}")
        else:
            above = df[df['threshold'] > peak_thresh]
            if len(above) > 0:
                selected = above.loc[above['lcc'].idxmin()]
                reason = f"No threshold gives LCC < {max_lcc}, using smallest LCC"
            else:
                selected = peak_row
                reason = "Fallback to peak (no better option)"

    sel_thresh = float(selected['threshold'])
    sel_comp = int(selected['components'])
    sel_lcc = int(selected['lcc'])

    # Verify component count
    if sel_comp < peak_comp * comp_fraction:
        log(f"  WARNING: Selected has {sel_comp} components "
            f"({sel_comp/peak_comp:.0%} of peak {peak_comp})")

    log(f"\n  Tier 2 SELECTED: threshold={sel_thresh:.2f}")
    log(f"    Reason: {reason}")
    log(f"    Nodes: {int(selected['nodes'])}, Edges: {int(selected['edges'])}")
    log(f"    Components: {sel_comp}, LCC: {sel_lcc}")

    return df, sel_thresh, reason


# =============================================================================
# UNIFIED THRESHOLD SELECTOR
# =============================================================================

def select_threshold(corr_diff, harris_all, harris_a3b, network_label=""):
    """
    Unified two-tier threshold selection.

    Returns:
        selected_threshold (float)
        tier_used (int): 1 or 2
        reason (str): Human-readable justification
        tier1_df (pd.DataFrame or None): Enrichment sweep data
        tier2_df (pd.DataFrame or None): Topology sweep data
    """
    banner(f"THRESHOLD SELECTION: {network_label}")

    genes_in_matrix = set(corr_diff.index)
    a3_present = [a3 for a3 in A3_GENES if a3 in genes_in_matrix]
    harris_present = harris_all & genes_in_matrix
    harris_a3b_present = harris_a3b & genes_in_matrix

    log(f"  DIFF matrix: {corr_diff.shape[0]} genes")
    log(f"  A3 genes present: {len(a3_present)} "
        f"({', '.join(A3_SHORT.get(a, a) for a in a3_present)})")
    log(f"  Harris interactors present: {len(harris_present)} / {len(harris_all)}")
    log(f"  Harris A3B-specific present: {len(harris_a3b_present)} / {len(harris_a3b)}")

    # ---- Check Tier 1 prerequisites ----
    tier1_possible = (len(a3_present) >= 1 and
                      len(harris_present) >= MIN_HARRIS_FOR_TIER1)

    if tier1_possible:
        log(f"\n  Tier 1 prerequisites MET "
            f"({len(a3_present)} A3 genes, {len(harris_present)} Harris interactors)")

        tier1_df, tier1_thresh, tier1_reason = tier1_enrichment_sweep(
            corr_diff, a3_present, harris_present, harris_a3b_present)

        if tier1_thresh is not None:
            log(f"\n  FINAL SELECTION: Tier 1, threshold={tier1_thresh:.2f}")
            log(f"    {tier1_reason}")

            # Also run Tier 2 for comparison
            log(f"\n  (Running Tier 2 for comparison only)")
            tier2_df, tier2_thresh, tier2_reason = tier2_topology(corr_diff)
            log(f"  Tier 2 would have selected: {tier2_thresh:.2f} ({tier2_reason})")

            return tier1_thresh, 1, tier1_reason, tier1_df, tier2_df

        else:
            log(f"\n  Tier 1 did not produce a significant result.")
            log(f"    ({tier1_reason})")
            log(f"  Falling back to Tier 2...")
            tier2_df, tier2_thresh, tier2_reason = tier2_topology(corr_diff)

            full_reason = (f"Tier 1 inconclusive ({tier1_reason}), "
                           f"Tier 2: {tier2_reason}")
            log(f"\n  FINAL SELECTION: Tier 2 (fallback), "
                f"threshold={tier2_thresh:.2f}")

            return tier2_thresh, 2, full_reason, tier1_df, tier2_df

    else:
        reasons = []
        if len(a3_present) < 1:
            reasons.append(f"no A3 genes in matrix")
        if len(harris_present) < MIN_HARRIS_FOR_TIER1:
            reasons.append(f"only {len(harris_present)} Harris interactors "
                           f"(need >= {MIN_HARRIS_FOR_TIER1})")
        prereq_reason = "; ".join(reasons)

        log(f"\n  Tier 1 prerequisites NOT MET: {prereq_reason}")
        log(f"  Using Tier 2 directly.")

        tier2_df, tier2_thresh, tier2_reason = tier2_topology(corr_diff)

        full_reason = (f"Tier 1 prerequisites not met ({prereq_reason}). "
                       f"Tier 2: {tier2_reason}")
        log(f"\n  FINAL SELECTION: Tier 2, threshold={tier2_thresh:.2f}")

        return tier2_thresh, 2, full_reason, None, tier2_df


# =============================================================================
# OUTPUT: PARAMETER FILE
# =============================================================================

def write_parameters(output_dir, threshold, tier, reason, network_label=""):
    """Write selected parameters to a key=value file for pipeline consumption."""
    os.makedirs(output_dir, exist_ok=True)
    param_path = os.path.join(output_dir, "selected_parameters.txt")

    with open(param_path, 'w') as f:
        f.write(f"# Auto-selected threshold parameters\n")
        f.write(f"# Generated: {datetime.now().isoformat()}\n")
        f.write(f"# Network: {network_label}\n")
        f.write(f"DIFF_THRESHOLD={threshold:.2f}\n")
        f.write(f"SELECTION_TIER={tier}\n")
        f.write(f"SELECTION_REASON={reason}\n")

    log(f"  [SAVE] Parameters -> {param_path}")
    return param_path


# =============================================================================
# PLOTS
# =============================================================================

def generate_plot(tier1_df, tier2_df, threshold, tier, network_label, output_dir):
    """Generate diagnostic plot showing both tiers."""
    os.makedirs(output_dir, exist_ok=True)

    n_panels = 2 if tier1_df is not None else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 6))
    if n_panels == 1:
        axes = [axes]

    fig.suptitle(f"Threshold Selection: {network_label}\n"
                 f"Selected: {threshold:.2f} (Tier {tier})",
                 fontsize=13, fontweight='bold')

    panel_idx = 0

    # ---- Tier 1 panel (if available) ----
    if tier1_df is not None:
        ax = axes[panel_idx]
        panel_idx += 1

        t1 = tier1_df.sort_values('threshold')
        threshs = t1['threshold'].values

        # Enrichment significance
        ax_left = ax
        ax_right = ax.twinx()

        # -log10(BH p) on left axis
        ax_left.plot(threshs, t1['neg_log10_bh'], 'o-', color='#e74c3c',
                     label='-log10(BH p)', markersize=6, linewidth=2)
        ax_left.axhline(-np.log10(ENRICHMENT_ALPHA), color='#e74c3c',
                        linestyle=':', alpha=0.5,
                        label=f'alpha={ENRICHMENT_ALPHA}')

        # Harris in neighborhood on right axis
        ax_right.plot(threshs, t1['harris_in_neighborhood'], 's--',
                      color='#f18f01', label='Harris in 2-hop',
                      markersize=5, alpha=0.7)
        ax_right.plot(threshs, t1['a3_in_net'], '^--', color='#2c3e50',
                      label='A3 in network', markersize=5, alpha=0.7)

        ax_left.axvline(threshold, color='green', linewidth=2, alpha=0.7,
                        label=f'Selected={threshold:.2f}')

        ax_left.set_xlabel("|delta-rho| threshold", fontsize=11)
        ax_left.set_ylabel("-log10(BH p-value)", color='#e74c3c', fontsize=11)
        ax_right.set_ylabel("Count", color='#f18f01', fontsize=11)
        ax_left.set_title("Tier 1: Harris Enrichment in A3 Neighborhood",
                          fontsize=11)

        # Combined legend
        lines_left, labels_left = ax_left.get_legend_handles_labels()
        lines_right, labels_right = ax_right.get_legend_handles_labels()
        ax_left.legend(lines_left + lines_right,
                       labels_left + labels_right,
                       fontsize=8, loc='upper right')
        ax_left.grid(True, alpha=0.3)

    # ---- Tier 2 panel ----
    if tier2_df is not None:
        ax = axes[panel_idx]
        t2 = tier2_df.sort_values('threshold')
        threshs = t2['threshold'].values

        ax_left = ax
        ax_right = ax.twinx()

        ax_left.plot(threshs, t2['components'], 'o-', color='#3498db',
                     label='Components', markersize=5)
        ax_right.plot(threshs, t2['lcc'], 's-', color='#e67e22',
                      label='LCC size', markersize=5)
        ax_right.axhline(TIER2_MAX_LCC, color='#e67e22', linestyle=':',
                         alpha=0.5, label=f'LCC limit={TIER2_MAX_LCC}')

        ax_left.axvline(threshold, color='green', linewidth=2, alpha=0.7,
                        label=f'Selected={threshold:.2f}')

        ax_left.set_xlabel("|delta-rho| threshold", fontsize=11)
        ax_left.set_ylabel("Connected components", color='#3498db',
                           fontsize=11)
        ax_right.set_ylabel("LCC size", color='#e67e22', fontsize=11)

        tier2_label = ("Tier 2: Topology (reference)"
                       if tier == 1 else "Tier 2: Topology (selected)")
        ax_left.set_title(tier2_label, fontsize=11)

        lines_left, labels_left = ax_left.get_legend_handles_labels()
        lines_right, labels_right = ax_right.get_legend_handles_labels()
        ax_left.legend(lines_left + lines_right,
                       labels_left + labels_right,
                       fontsize=8, loc='upper left')
        ax_left.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plot_path = os.path.join(output_dir, "threshold_selection_plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Plot -> {plot_path}")


# =============================================================================
# RUN FOR ONE NETWORK
# =============================================================================

def run_single(corr_pkl_path, harris_all, harris_a3b, output_dir,
               network_label=""):
    """Run threshold selection for a single network."""
    corr_diff = load_diff_matrix(corr_pkl_path)
    if corr_diff is None:
        return None

    threshold, tier, reason, tier1_df, tier2_df = select_threshold(
        corr_diff, harris_all, harris_a3b, network_label=network_label)

    # Write parameters
    write_parameters(output_dir, threshold, tier, reason, network_label)

    # Save enrichment data if available
    if tier1_df is not None:
        enrich_path = os.path.join(output_dir, "tier1_enrichment_sweep.csv")
        tier1_df.to_csv(enrich_path, index=False)
        log(f"  [SAVE] Enrichment data -> {enrich_path}")

    if tier2_df is not None:
        topo_path = os.path.join(output_dir, "tier2_topology_sweep.csv")
        tier2_df.to_csv(topo_path, index=False)
        log(f"  [SAVE] Topology data -> {topo_path}")

    # Plot
    generate_plot(tier1_df, tier2_df, threshold, tier, network_label, output_dir)

    return threshold, tier, reason


# =============================================================================
# RUN FOR THREE SC NETWORKS
# =============================================================================

def run_three_networks(harris_all, harris_a3b):
    """Run threshold selection for all three Figure 4 SC networks."""
    banner("THREE-NETWORK MODE: Figure 4 Single-Cell Networks")

    results = {}

    for net_name, net_info in SC_NETWORKS.items():
        net_dir = net_info['dir']
        corr_pkl = os.path.join(net_dir, "03_correlation_networks",
                                "corr_matrices", "SC_corr_DIFF.pkl")
        out_dir = os.path.join(net_dir, "THRESHOLD_SELECTION")

        result = run_single(corr_pkl, harris_all, harris_a3b,
                           out_dir, network_label=net_info['label'])
        if result is not None:
            results[net_name] = result

    # Cross-network summary
    banner("CROSS-NETWORK SUMMARY")
    log(f"\n  {'Network':>20s}  {'Threshold':>10s}  {'Tier':>5s}  Reason")
    log(f"  {'-'*20}  {'-'*10}  {'-'*5}  {'-'*40}")
    for net_name, (thresh, tier, reason) in results.items():
        short_reason = reason[:60] + "..." if len(reason) > 60 else reason
        log(f"  {net_name:>20s}  {thresh:>10.2f}  {tier:>5d}  {short_reason}")

    # Unified methods statement
    log(f"\n  Methods statement:")
    tier1_nets = [n for n, (t, tier, r) in results.items() if tier == 1]
    tier2_nets = [n for n, (t, tier, r) in results.items() if tier == 2]

    if tier1_nets and tier2_nets:
        log(f"  Tier 1 (enrichment) applied to: {', '.join(tier1_nets)}")
        log(f"  Tier 2 (topology) applied to: {', '.join(tier2_nets)}")
    elif tier1_nets:
        log(f"  All networks used Tier 1 (enrichment-based selection)")
    else:
        log(f"  All networks used Tier 2 (topology-driven selection)")

    return results


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Two-tier DIFF threshold selection")
    parser.add_argument('--mode', choices=['single', 'three_networks'],
                        default='three_networks',
                        help='Run mode (default: three_networks)')
    parser.add_argument('--corr_pkl',
                        help='Path to DIFF correlation matrix pickle (single mode)')
    parser.add_argument('--harris', default=HARRIS_ALL_PATH,
                        help='Path to Harris interactors file')
    parser.add_argument('--harris_a3b', default=HARRIS_A3B_PATH,
                        help='Path to Harris A3B-specific interactors file')
    parser.add_argument('--output_dir',
                        help='Output directory (single mode)')
    parser.add_argument('--label', default='',
                        help='Network label (single mode)')

    args = parser.parse_args()

    t0 = datetime.now()
    banner("AUTO-SELECT THRESHOLD (TWO-TIER FRAMEWORK)")
    log(f"  Start: {t0}")
    log(f"  Mode: {args.mode}")

    # Load Harris interactors
    harris_all, harris_a3b = load_harris(args.harris, args.harris_a3b)
    log(f"  Harris interactors (all): {len(harris_all)}")
    log(f"  Harris interactors (A3B-specific): {len(harris_a3b)}")

    if args.mode == 'three_networks':
        results = run_three_networks(harris_all, harris_a3b)
    else:
        if not args.corr_pkl or not args.output_dir:
            print("ERROR: --corr_pkl and --output_dir required for single mode")
            sys.exit(1)
        result = run_single(args.corr_pkl, harris_all, harris_a3b,
                           args.output_dir, network_label=args.label)

    # Save full report
    global report_lines
    report_path = os.path.join(
        FIG4_ROOT if args.mode == 'three_networks' else args.output_dir,
        "DIAGNOSTIC_AUDIT",
        "threshold_selection_report.txt"
    )
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    log(f"\n  [SAVE] Full report -> {report_path}")

    banner("THRESHOLD SELECTION COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")


if __name__ == "__main__":
    main()
