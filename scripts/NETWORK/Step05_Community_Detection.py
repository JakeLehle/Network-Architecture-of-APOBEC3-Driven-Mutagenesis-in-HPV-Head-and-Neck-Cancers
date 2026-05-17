#!/usr/bin/env python3
"""
Step05_Community_Detection.py

Run Leiden community detection on the full DIFF network across multiple
resolutions. Evaluate stability (ARI/NMI across runs), merge small
communities (within large components only), and save community assignments
+ network visualizations.

THRESHOLD SELECTION: Maximum Fragmentation Rate
  Sweeps DIFF correlation thresholds and computes the forward difference
  in connected component count between consecutive steps. Among intervals
  where the upper threshold is A3-valid (A3A and A3B have degree >= 1),
  selects the upper threshold of the interval with the largest positive
  delta-component count. Ties broken by the lower (more inclusive) upper
  threshold.

  Rationale: The threshold where connected components increase most
  rapidly is where differential co-expression modules are most actively
  separating. Selecting at this point preserves A3 pathway connectivity
  while capturing the natural resolution of the underlying biology.

COMMUNITY DETECTION: Full-Network Leiden
  Leiden runs on the entire graph (all connected components), not just
  the largest connected component. Large components get subdivided into
  pathway-scale communities. Small satellite components (where A3 genes
  often reside) remain as their own communities.

  Merge logic is component-aware: small communities within large
  connected components are merged (Leiden artifacts), but satellite
  communities that ARE entire small components are preserved as
  biologically meaningful units.

Resolution selection uses composite score = modularity x ARI x evenness
(unchanged from previous version).

Corresponds to original pipeline Step 15.

Input:
    05_correlation_networks/{cancer_type}/corr_matrices/{ct}_corr_DIFF.pkl
    03_differential_expression/{cancer_type}/{ct}_diffexpr_stats.csv
    01_cleaned_expression/ensg_to_symbol.json

Output (-> data/FIG_2/06_communities/{cancer_type}/):
    sweep/                     -- per-resolution community assignments + plots
    {ct}_resolution_sweep.csv  -- stability metrics across resolutions
    {ct}_sweep_diagnostics.png -- 6-panel diagnostic plots
    {ct}_best_partition.csv    -- final community assignments (best resolution)
    {ct}_community_gene_lists.csv  -- genes per community
    {ct}_threshold_sweep.csv   -- threshold selection diagnostics
    {ct}_selected_parameters.txt -- threshold + resolution + metrics
    {ct}_community_summary.txt -- human-readable community report

Usage:
    python Step05_Community_Detection.py
"""

import os
import json
import pickle
import shutil
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from networkx.algorithms.community.quality import modularity as nx_modularity

from network_config import (
    CANCER_TYPES, A3_GENES, BIOMARKERS, A3_ID_TO_ALIAS,
    COMMUNITY_METHOD, COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION,
    COMMUNITY_BASE_SEED, USE_LARGEST_COMPONENT,
    TARGET_BIG_COMMUNITIES, MIN_COMMUNITY_SIZE,
    DIFF_THRESHOLD, SWEEP_THRESHOLDS,
    RAW_P_THRESHOLD, FORCE_KEEP_A3,
    DIR_01_CLEANED, DIR_03_DIFFEXPR, DIR_04_NETWORKS, DIR_05_COMMUNITIES,
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
# THRESHOLD SELECTION: Maximum Fragmentation Rate
# =============================================================================

def auto_select_threshold(corr_diff, a3_seeds, sweep_thresholds):
    """
    Select DIFF threshold at the point of maximum network fragmentation rate.

    Among consecutive threshold pairs where the upper threshold is A3-valid
    (all seed genes have degree >= 1), find the interval with the largest
    positive increase in connected component count. Select the UPPER
    threshold of that interval (where the fragmentation is realized).

    Tiebreaker: if multiple intervals share the same max delta-comp, select
    the one with the lower upper threshold (more inclusive network).

    Fallback: if no positive fragmentation occurs among A3-valid thresholds
    (e.g., network remains monolithic), select the highest A3-valid threshold.

    Parameters
    ----------
    corr_diff : pd.DataFrame
        DIFF correlation matrix (genes x genes).
    a3_seeds : list of str
        Gene IDs (ENSG for TCGA) for A3A and/or A3B.
    sweep_thresholds : list of float
        Thresholds to evaluate.

    Returns
    -------
    selected_threshold : float or None
        The selected threshold, or None if no A3-valid threshold exists.
    sweep_df : pd.DataFrame
        Diagnostic table with per-threshold metrics and delta-comp values.
    """
    matrix_genes = set(corr_diff.index)

    # Verify which A3 seeds are actually in the matrix
    seeds_in_matrix = [g for g in a3_seeds if g in matrix_genes]
    seeds_missing = [g for g in a3_seeds if g not in matrix_genes]

    log(f"  A3 seed genes for threshold selection:")
    for g in seeds_in_matrix:
        alias = A3_ID_TO_ALIAS.get(g, g)
        log(f"    {g} ({alias}) -- IN matrix")
    for g in seeds_missing:
        alias = A3_ID_TO_ALIAS.get(g, g)
        log(f"    {g} ({alias}) -- MISSING from matrix")

    if not seeds_in_matrix:
        log("  WARNING: No A3 seed genes found in DIFF matrix. Cannot auto-select.")
        return None, pd.DataFrame()

    # ---- Sweep thresholds ----
    rows = []
    for t in sorted(sweep_thresholds):
        G = build_weighted_graph(corr_diff, t)
        G_noiso = remove_isolated(G)

        n_nodes = G_noiso.number_of_nodes()
        n_edges = G_noiso.number_of_edges()

        if n_nodes > 0:
            components = list(nx.connected_components(G_noiso))
            n_comp = len(components)
            lcc_size = len(max(components, key=len))
        else:
            n_comp = 0
            lcc_size = 0

        # Check degree of each A3 seed
        seed_degrees = {}
        all_seeds_connected = True
        for g in seeds_in_matrix:
            deg = G_noiso.degree(g) if g in G_noiso else 0
            seed_degrees[g] = deg
            if deg < 1:
                all_seeds_connected = False

        rows.append({
            "threshold": t,
            "nodes": n_nodes,
            "edges": n_edges,
            "components": n_comp,
            "lcc_size": lcc_size,
            "all_seeds_connected": all_seeds_connected,
            **{f"deg_{A3_ID_TO_ALIAS.get(g, g)}": d for g, d in seed_degrees.items()},
        })

    sweep_df = pd.DataFrame(rows).sort_values("threshold").reset_index(drop=True)

    # ---- Log the sweep table ----
    header = (f"\n  {'Thresh':>7s}  {'Nodes':>6s}  {'Edges':>7s}  {'Comp':>5s}  "
              f"{'LCC':>5s}  {'A3 OK':>5s}  ")
    sep_line = (f"  {'-----':>7s}  {'-----':>6s}  {'-----':>7s}  {'----':>5s}  "
                f"{'---':>5s}  {'-----':>5s}  ")
    for g in seeds_in_matrix:
        alias = A3_ID_TO_ALIAS.get(g, g)
        header += f"{'deg_' + alias:>8s}  "
        sep_line += f"{'--------':>8s}  "
    log(header)
    log(sep_line)

    for _, row in sweep_df.iterrows():
        marker = " <--" if row["all_seeds_connected"] else ""
        line = (f"  {row['threshold']:>7.2f}  {row['nodes']:>6.0f}  "
                f"{row['edges']:>7.0f}  {row['components']:>5.0f}  "
                f"{row['lcc_size']:>5.0f}  {'YES' if row['all_seeds_connected'] else 'no':>5s}  ")
        for g in seeds_in_matrix:
            alias = A3_ID_TO_ALIAS.get(g, g)
            col = f"deg_{alias}"
            line += f"{row[col]:>8.0f}  "
        log(line + marker)

    # ---- Compute forward delta-comp between consecutive thresholds ----
    sweep_df["delta_comp"] = sweep_df["components"].diff().shift(-1)
    # delta_comp at row i = comp[i+1] - comp[i]
    # We want to select the UPPER threshold of the best interval,
    # so we track which upper threshold each delta corresponds to.

    log(f"\n  Forward delta-comp (fragmentation rate):")
    log(f"  {'Interval':>15s}  {'Δcomp':>6s}  {'Upper A3-valid':>14s}")
    log(f"  {'--------':>15s}  {'-----':>6s}  {'--------------':>14s}")

    candidates = []
    for i in range(len(sweep_df) - 1):
        t_lower = sweep_df.loc[i, "threshold"]
        t_upper = sweep_df.loc[i + 1, "threshold"]
        delta = int(sweep_df.loc[i + 1, "components"] - sweep_df.loc[i, "components"])
        upper_valid = bool(sweep_df.loc[i + 1, "all_seeds_connected"])

        marker = ""
        if upper_valid and delta > 0:
            candidates.append({
                "t_lower": t_lower,
                "t_upper": t_upper,
                "delta_comp": delta,
            })
            marker = " <-- candidate"

        log(f"  {t_lower:.2f} -> {t_upper:.2f}  {delta:>+6d}  "
            f"{'YES' if upper_valid else 'no':>14s}{marker}")

    # ---- Select threshold ----
    if not candidates:
        # Fallback: highest A3-valid threshold
        valid = sweep_df[sweep_df["all_seeds_connected"]]
        if len(valid) == 0:
            log("\n  WARNING: No A3-valid threshold exists. Cannot auto-select.")
            return None, sweep_df

        selected = float(valid["threshold"].max())
        reason = ("fallback: highest A3-valid threshold "
                  "(no positive fragmentation among A3-valid intervals)")
        log(f"\n  No positive fragmentation candidates found.")
    else:
        # Find maximum delta-comp; tiebreaker = lower upper threshold
        max_delta = max(c["delta_comp"] for c in candidates)
        ties = [c for c in candidates if c["delta_comp"] == max_delta]
        best = min(ties, key=lambda c: c["t_upper"])
        selected = best["t_upper"]
        reason = (f"max fragmentation rate (delta_comp=+{best['delta_comp']}) "
                  f"at interval {best['t_lower']:.2f} -> {best['t_upper']:.2f}")

        if len(ties) > 1:
            log(f"\n  Tiebreaker applied: {len(ties)} intervals with "
                f"delta_comp=+{max_delta}")
            for t in ties:
                log(f"    {t['t_lower']:.2f} -> {t['t_upper']:.2f}")
            log(f"  Selected lower upper threshold: {selected:.2f}")

    sel_row = sweep_df[sweep_df["threshold"] == selected].iloc[0]
    log(f"\n  SELECTED threshold: {selected:.2f}")
    log(f"    Reason: {reason}")
    log(f"    Nodes: {int(sel_row['nodes'])}, Edges: {int(sel_row['edges'])}, "
        f"LCC: {int(sel_row['lcc_size'])}")
    for g in seeds_in_matrix:
        alias = A3_ID_TO_ALIAS.get(g, g)
        col = f"deg_{alias}"
        log(f"    {alias} degree: {int(sel_row[col])}")

    return selected, sweep_df


# =============================================================================
# A3 SEED IDENTIFICATION
# =============================================================================

def identify_a3_seeds(de_dir, cancer_type, force_keep=False):
    """
    Identify which of A3A and A3B should be used as seed genes for threshold
    selection.

    Two modes:
      - force_keep=False: only A3 genes passing DE at raw p < 0.05 are seeds
      - force_keep=True: A3A and A3B are seeds if they exist in the DE results
        (regardless of p-value), because FORCE_KEEP_A3 guarantees they are in
        the gene list and correlation matrix

    For TCGA, genes are ENSG IDs.

    Parameters
    ----------
    de_dir : str
        Path to DE results directory for this cancer type.
    cancer_type : str
        Cancer type label (e.g. "TCGA-HNSC").
    force_keep : bool
        If True, include A3A/A3B as seeds even if they fail DE.

    Returns
    -------
    seeds : list of str
        ENSG IDs for A3A and/or A3B to use as threshold anchors.
    """
    # Look for DE results file
    de_candidates = [
        os.path.join(de_dir, f"{cancer_type}_diffexpr_stats.csv"),
        os.path.join(de_dir, f"{cancer_type}_diffexpr_stats.tsv"),
    ]

    de_path = None
    for p in de_candidates:
        if os.path.exists(p):
            de_path = p
            break

    # A3A and A3B ENSG IDs
    a3a_ensg = "ENSG00000128383"
    a3b_ensg = "ENSG00000179750"

    if de_path is None:
        log(f"  WARNING: No DE results file found in {de_dir}")
        log(f"  Defaulting to both A3A and A3B as seeds")
        return [a3a_ensg, a3b_ensg]

    log(f"  Loading DE results: {de_path}")
    if force_keep:
        log(f"  FORCE_KEEP_A3 is ON: A3A/A3B will be used as seeds "
            f"regardless of DE status")
    sep = "\t" if de_path.endswith(".tsv") else ","
    de_df = pd.read_csv(de_path, sep=sep)

    # Identify gene column
    gene_col = None
    for col in ["gene", "gene_symbol", "Gene", "gene_id"]:
        if col in de_df.columns:
            gene_col = col
            break
    if gene_col is None:
        log(f"  WARNING: Cannot identify gene column. "
            f"Columns: {list(de_df.columns)}")
        return [a3a_ensg, a3b_ensg]

    # Identify p-value column
    pval_col = None
    for col in ["p_value", "pvalue", "pval", "raw_p", "PValue"]:
        if col in de_df.columns:
            pval_col = col
            break
    if pval_col is None:
        log(f"  WARNING: Cannot identify p-value column. "
            f"Columns: {list(de_df.columns)}")
        return [a3a_ensg, a3b_ensg]

    # Check A3A and A3B
    seeds = []
    for ensg, label in [(a3a_ensg, "A3A"), (a3b_ensg, "A3B")]:
        match = de_df[de_df[gene_col] == ensg]
        if len(match) > 0:
            row = match.iloc[0]
            pval = float(row[pval_col])
            if pval < RAW_P_THRESHOLD:
                log(f"    {label} ({ensg}): PASSED DE (p={pval:.2e})")
                seeds.append(ensg)
            elif force_keep:
                log(f"    {label} ({ensg}): FAILED DE (p={pval:.2e}) "
                    f"-- included via FORCE_KEEP")
                seeds.append(ensg)
            else:
                log(f"    {label} ({ensg}): FAILED DE (p={pval:.2e})")
        else:
            log(f"    {label} ({ensg}): not found in DE results")

    if not seeds:
        log(f"  WARNING: Neither A3A nor A3B available as seeds.")

    return seeds


# =============================================================================
# COMMUNITY DETECTION FUNCTIONS
# =============================================================================

def detect_leiden(G, seed=42, resolution=1.0, weight="abs_weight"):
    """Run Leiden community detection. Returns {node: community_id} dict."""
    try:
        import leidenalg
        import igraph as ig
    except ImportError:
        raise ImportError(
            "leidenalg and python-igraph required. "
            "Install: pip install leidenalg igraph"
        )

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


def merge_small_communities(comm_map, G, k_keep=8, min_size=10):
    """
    Component-aware community merge for full-network clustering.

    Within large connected components, small Leiden communities are merged
    into an 'Other' bin (standard behavior). Satellite communities whose
    members all belong to a small connected component are preserved as-is,
    since they represent biologically meaningful co-expression units that
    are naturally disconnected from the main network.

    Parameters
    ----------
    comm_map : dict
        {node: community_id} from Leiden.
    G : nx.Graph
        The full graph (all components) used for community detection.
    k_keep : int
        Maximum number of large communities to retain within the main
        component(s). Smaller communities within large components are
        merged into 'Other'.
    min_size : int
        Minimum community size within large components. Communities
        smaller than this (within large components) are merged.

    Returns
    -------
    merged : dict
        {node: community_id} after merge.
    raw_sizes : dict
        {community_id: size} before merge.
    merged_sizes : dict
        {community_id: size} after merge.
    satellite_comm_ids : set
        Community IDs that were preserved as satellites.
    """
    # Map each node to its connected component
    components = list(nx.connected_components(G))
    node_to_comp_size = {}
    for comp in components:
        sz = len(comp)
        for node in comp:
            node_to_comp_size[node] = sz

    # Compute raw community sizes
    raw_sizes = {}
    for n, cid in comm_map.items():
        raw_sizes[cid] = raw_sizes.get(cid, 0) + 1

    # Classify each community: satellite vs main-component
    # A community is a satellite if ALL its members belong to a connected
    # component smaller than min_size (i.e., the component IS the community)
    comm_members = {}
    for n, cid in comm_map.items():
        comm_members.setdefault(cid, []).append(n)

    satellite_comm_ids = set()
    main_comm_ids = set()
    for cid, members in comm_members.items():
        max_comp_size = max(node_to_comp_size[n] for n in members)
        if max_comp_size < min_size:
            satellite_comm_ids.add(cid)
        else:
            main_comm_ids.add(cid)

    # Among main-component communities, keep top k_keep by size
    main_sizes = {cid: raw_sizes[cid] for cid in main_comm_ids}
    top_main = sorted(main_sizes.keys(),
                      key=lambda c: main_sizes[c], reverse=True)[:k_keep]
    top_main = set(top_main)

    # Assign "Other" ID (higher than any existing community)
    other_id = max(raw_sizes.keys()) + 1 if raw_sizes else 0

    # Build merged map
    merged = {}
    for n, cid in comm_map.items():
        if cid in satellite_comm_ids:
            # Satellite: preserve as-is
            merged[n] = cid
        elif cid in top_main and raw_sizes[cid] >= min_size:
            # Large main-component community: keep
            merged[n] = cid
        else:
            # Small main-component community: merge into Other
            merged[n] = other_id

    # Compute merged sizes
    merged_sizes = {}
    for n, cid in merged.items():
        merged_sizes[cid] = merged_sizes.get(cid, 0) + 1

    return merged, raw_sizes, merged_sizes, satellite_comm_ids


def symbol_or_self(ensg, mapping):
    """Look up gene symbol, return ENSG ID if not found."""
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

    banner(f"[STEP 15] Community Detection -- {cancer_type}", char="=")

    cancer_dir = ensure_dir(os.path.join(DIR_05_COMMUNITIES, cancer_type))
    sweep_dir = ensure_dir(os.path.join(cancer_dir, "sweep"))

    # =========================================================================
    # STEP 15.0: Identify A3 seeds and auto-select DIFF threshold
    # =========================================================================
    banner("[STEP 15.0] DIFF threshold selection (max fragmentation rate)")

    # Identify A3 seeds from DE results
    de_cancer_dir = os.path.join(DIR_03_DIFFEXPR, cancer_type)
    a3_seeds = identify_a3_seeds(de_cancer_dir, cancer_type,
                                 force_keep=FORCE_KEEP_A3)

    corr_diff = None  # Track whether we have the correlation matrix

    if not a3_seeds:
        log("WARNING: No A3 seed genes available. "
            "Falling back to config threshold.")
        active_threshold = DIFF_THRESHOLD
    else:
        log(f"  A3 seeds: {a3_seeds}")

        # Load DIFF correlation matrix
        corr_path = os.path.join(
            DIR_04_NETWORKS, cancer_type, "corr_matrices",
            f"{cancer_type}_corr_DIFF.pkl"
        )
        if not os.path.exists(corr_path):
            log(f"  WARNING: DIFF correlation matrix not found: {corr_path}")
            log(f"  Falling back to config threshold: {DIFF_THRESHOLD}")
            active_threshold = DIFF_THRESHOLD
        else:
            log(f"  Loading DIFF matrix: {corr_path}")
            with open(corr_path, "rb") as f:
                corr_diff = pickle.load(f)
            log(f"  DIFF matrix shape: {corr_diff.shape}")

            selected_threshold, threshold_sweep_df = auto_select_threshold(
                corr_diff, a3_seeds, SWEEP_THRESHOLDS
            )

            if selected_threshold is not None:
                active_threshold = selected_threshold
            else:
                active_threshold = DIFF_THRESHOLD
                log(f"  Using config fallback: {active_threshold}")

            # Save threshold sweep diagnostics
            if len(threshold_sweep_df) > 0:
                thresh_csv = os.path.join(
                    cancer_dir, f"{cancer_type}_threshold_sweep.csv"
                )
                threshold_sweep_df.to_csv(thresh_csv, index=False)
                log(f"  [SAVE] Threshold sweep -> {thresh_csv}")

    log(f"\n  Active DIFF threshold: {active_threshold}")

    # =========================================================================
    # STEP 15.0b: Build DIFF graph at the selected threshold
    # =========================================================================
    banner("[STEP 15.0b] Build DIFF graph at selected threshold")

    if corr_diff is None:
        # Need to load correlation matrix (fallback path)
        corr_path = os.path.join(
            DIR_04_NETWORKS, cancer_type, "corr_matrices",
            f"{cancer_type}_corr_DIFF.pkl"
        )
        if not os.path.exists(corr_path):
            # Last resort: try loading pre-built graph
            pkl_path = os.path.join(
                DIR_04_NETWORKS, cancer_type, "graph_objects",
                f"{cancer_type}_G_diff_noiso.gpickle"
            )
            if not os.path.exists(pkl_path):
                log(f"[SKIP] Neither DIFF matrix nor pre-built graph "
                    f"found for {cancer_type}")
                continue
            log(f"  Loading pre-built graph (legacy fallback): {pkl_path}")
            with open(pkl_path, "rb") as f:
                G_diff_noiso = pickle.load(f)
        else:
            with open(corr_path, "rb") as f:
                corr_diff = pickle.load(f)

    # Build graph from correlation matrix if available
    if corr_diff is not None:
        log(f"  Building graph at |delta-rho| >= {active_threshold}")
        G_diff = build_weighted_graph(corr_diff, active_threshold)
        G_diff_noiso = remove_isolated(G_diff)

        # Save rebuilt graph for downstream steps
        graph_dir = ensure_dir(os.path.join(
            DIR_04_NETWORKS, cancer_type, "graph_objects"
        ))
        graph_save_path = os.path.join(
            graph_dir, f"{cancer_type}_G_diff_noiso.gpickle"
        )
        with open(graph_save_path, "wb") as f:
            pickle.dump(G_diff_noiso, f)
        log(f"  [SAVE] Updated DIFF graph -> {graph_save_path}")

        # Free memory
        del corr_diff

    n_nodes = G_diff_noiso.number_of_nodes()
    n_edges = G_diff_noiso.number_of_edges()
    log(f"  DIFF graph: nodes={n_nodes}, edges={n_edges}")

    if n_nodes < 5 or n_edges < 3:
        log("[SKIP] Graph too small for community detection")
        continue

    # =========================================================================
    # STEP 15.1: Prepare graph for community detection (FULL NETWORK)
    # =========================================================================
    banner("[STEP 15.1] Prepare full network for community detection")

    G_comm = G_diff_noiso.copy()

    # Ensure abs_weight on all edges
    for u, v, d in G_comm.edges(data=True):
        d["abs_weight"] = abs(float(d.get("abs_weight", d.get("weight", 0.0))))

    # Report connected component structure
    components = list(nx.connected_components(G_comm))
    components_sorted = sorted(components, key=len, reverse=True)
    n_components = len(components_sorted)
    lcc_size = len(components_sorted[0])

    log(f"  Full network: {G_comm.number_of_nodes()} nodes, "
        f"{G_comm.number_of_edges()} edges")
    log(f"  Connected components: {n_components}")
    log(f"  LCC size: {lcc_size} nodes "
        f"({100*lcc_size/G_comm.number_of_nodes():.1f}% of network)")

    # Component size distribution
    comp_sizes = [len(c) for c in components_sorted]
    log(f"  Component size distribution:")
    log(f"    Largest 5: {comp_sizes[:5]}")
    if n_components > 5:
        small_count = sum(1 for s in comp_sizes if s < MIN_COMMUNITY_SIZE)
        log(f"    Components < {MIN_COMMUNITY_SIZE} nodes: "
            f"{small_count} ({sum(s for s in comp_sizes if s < MIN_COMMUNITY_SIZE)} genes)")

    # Report A3 seed status
    log(f"\n  A3 seed positions in full network:")
    for g in a3_seeds:
        alias = A3_ID_TO_ALIAS.get(g, g)
        if g in G_comm:
            deg = G_comm.degree(g)
            # Find which component
            for i, comp in enumerate(components_sorted):
                if g in comp:
                    log(f"    {alias}: component {i} "
                        f"(size={len(comp)}), degree={deg}")
                    break
        else:
            log(f"    {alias}: NOT in graph")

    nodes_list = sorted(G_comm.nodes())

    # =========================================================================
    # STEP 15.2: Resolution sweep
    # =========================================================================
    banner("[STEP 15.2] Resolution sweep (full-network Leiden)")

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

            cm = detect_leiden(G_comm, seed=seed, resolution=r,
                               weight="abs_weight")
            labels = partition_to_labels(cm, nodes_list)
            labels_list.append(labels)

            n_comms = len(set(cm.values()))
            ncomms_list.append(n_comms)

            # Compute modularity
            comm_sets = {}
            for n, c in cm.items():
                comm_sets.setdefault(c, set()).add(n)
            mod = nx_modularity(G_comm, list(comm_sets.values()),
                                weight="abs_weight")
            modularities.append(mod)

            if mod > best_mod:
                best_mod = mod
                best_cm = cm.copy()

        # ---- Save best raw partition ----
        raw_csv = os.path.join(
            res_dir, f"{cancer_type}_partition_raw_{res_tag}.csv"
        )
        pd.DataFrame({
            "gene": list(best_cm.keys()),
            "gene_symbol": [symbol_or_self(g, ensg_to_symbol)
                            for g in best_cm.keys()],
            "community": list(best_cm.values())
        }).to_csv(raw_csv, index=False)

        # ---- Component-aware merge ----
        cm_merged, raw_sizes, merged_sizes, satellite_ids = \
            merge_small_communities(
                best_cm, G_comm,
                k_keep=TARGET_BIG_COMMUNITIES,
                min_size=MIN_COMMUNITY_SIZE
            )

        merged_csv = os.path.join(
            res_dir, f"{cancer_type}_partition_merged_{res_tag}.csv"
        )
        pd.DataFrame({
            "gene": list(cm_merged.keys()),
            "gene_symbol": [symbol_or_self(g, ensg_to_symbol)
                            for g in cm_merged.keys()],
            "community": list(cm_merged.values())
        }).to_csv(merged_csv, index=False)

        # ---- Gene lists per community ----
        comm_to_genes = {}
        for g, c in cm_merged.items():
            comm_to_genes.setdefault(c, []).append(
                symbol_or_self(g, ensg_to_symbol)
            )

        gene_list_rows = []
        for c, genes in sorted(comm_to_genes.items(),
                                key=lambda x: len(x[1]), reverse=True):
            is_satellite = c in satellite_ids
            gene_list_rows.append({
                "community": c,
                "size": len(genes),
                "is_satellite": is_satellite,
                "genes": ";".join(sorted(genes))
            })
        gene_list_csv = os.path.join(
            res_dir, f"{cancer_type}_gene_lists_{res_tag}.csv"
        )
        pd.DataFrame(gene_list_rows).to_csv(gene_list_csv, index=False)

        # ---- Community size bar plot ----
        sizes_sorted = sorted(merged_sizes.items(),
                              key=lambda x: x[1], reverse=True)
        fig, ax = plt.subplots(figsize=(10, 5))
        bar_colors = []
        bar_labels = []
        for cid, sz in sizes_sorted:
            if cid in satellite_ids:
                bar_colors.append("#E8913A")  # orange for satellites
                bar_labels.append(f"C{cid}*")
            else:
                bar_colors.append("#4682B4")  # steel blue for main
                bar_labels.append(f"C{cid}")

        ax.bar(bar_labels, [v for _, v in sizes_sorted],
               color=bar_colors, edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Community ID (* = satellite)")
        ax.set_ylabel("Number of genes")
        ax.set_title(f"{cancer_type} | Community sizes | res={r:.2f} | "
                     f"full network")
        plt.tight_layout()
        plt.savefig(os.path.join(
            res_dir, f"{cancer_type}_comm_sizes_{res_tag}.png"), dpi=300)
        plt.close()

        # ---- Network plot colored by community ----
        pos = nx.spring_layout(G_comm, seed=COMMUNITY_BASE_SEED,
                               weight="abs_weight")

        unique_comms = sorted(set(cm_merged.values()))
        cmap = plt.colormaps["tab20"].resampled(max(len(unique_comms), 1))
        comm_colors = {c: cmap(i) for i, c in enumerate(unique_comms)}

        fig, ax = plt.subplots(figsize=(14, 12))
        nx.draw_networkx_edges(G_comm, pos, ax=ax, alpha=0.1, width=0.3,
                               edge_color="gray")

        for comm_id in unique_comms:
            comm_nodes = [n for n in G_comm.nodes()
                          if cm_merged.get(n) == comm_id]
            sat_tag = "*" if comm_id in satellite_ids else ""
            nx.draw_networkx_nodes(
                G_comm, pos, nodelist=comm_nodes, ax=ax,
                node_size=30, node_color=[comm_colors[comm_id]],
                alpha=0.7,
                label=f"C{comm_id}{sat_tag} (n={len(comm_nodes)})"
            )

        # Highlight A3 genes
        a3_in_graph = [g for g in A3_GENES if g in G_comm.nodes()]
        if a3_in_graph:
            nx.draw_networkx_nodes(
                G_comm, pos, nodelist=a3_in_graph, ax=ax,
                node_size=150, node_color="gold",
                edgecolors="black", linewidths=1.5
            )
            labels = {n: symbol_or_self(n, ensg_to_symbol)
                      for n in a3_in_graph}
            nx.draw_networkx_labels(G_comm, pos, labels=labels, ax=ax,
                                    font_size=8, font_weight="bold")

        ax.set_title(f"{cancer_type} | Full-network communities (merged) | "
                     f"res={r:.2f}")
        ax.legend(fontsize=7, loc="upper right", ncol=2)
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(
            res_dir, f"{cancer_type}_DIFF_communities_{res_tag}.png"),
            dpi=300)
        plt.close()

        # ---- Log summary ----
        n_satellite = len(satellite_ids)
        n_main = len(merged_sizes) - n_satellite
        log(f"  Saved outputs for res={r:.2f}: "
            f"comms={np.mean(ncomms_list):.1f} (raw), "
            f"{len(merged_sizes)} (merged: {n_main} main + "
            f"{n_satellite} satellite), mod={best_mod:.4f}")

        # ---- Stability metrics ----
        aris = []
        nmis = []
        for a, b in itertools.combinations(range(RUNS_PER_RESOLUTION), 2):
            aris.append(adjusted_rand_score(labels_list[a], labels_list[b]))
            nmis.append(normalized_mutual_info_score(
                labels_list[a], labels_list[b]))

        # Evenness: based on merged sizes (excluding satellites)
        n_total_genes = len(cm_merged)
        largest_comm_size = max(merged_sizes.values()) if merged_sizes else 1
        n_real_comms = len(merged_sizes)

        sweep_rows.append({
            "resolution": float(r),
            "n_runs": RUNS_PER_RESOLUTION,
            "ncomms_mean": float(np.mean(ncomms_list)),
            "ncomms_std": float(np.std(ncomms_list)),
            "ncomms_merged": len(merged_sizes),
            "n_satellite": n_satellite,
            "modularity_mean": float(np.mean(modularities)),
            "modularity_std": float(np.std(modularities)),
            "modularity_best": float(best_mod),
            "ARI_mean": float(np.mean(aris)),
            "ARI_std": float(np.std(aris)),
            "NMI_mean": float(np.mean(nmis)),
            "NMI_std": float(np.std(nmis)),
            "n_total_genes": n_total_genes,
            "largest_community": largest_comm_size,
            "n_communities_merged": n_real_comms,
        })

    # =========================================================================
    # STEP 15.8: Save sweep summary and select best resolution
    # =========================================================================
    banner("[STEP 15.8] Save sweep summary")

    sweep_df = pd.DataFrame(sweep_rows).sort_values("resolution") \
                                        .reset_index(drop=True)

    # Composite selection metrics
    sweep_df["evenness"] = (
        1.0 - (sweep_df["largest_community"] / sweep_df["n_total_genes"])
    )
    sweep_df["composite_score"] = (
        sweep_df["modularity_best"]
        * sweep_df["ARI_mean"]
        * sweep_df["evenness"]
    )
    sweep_df["delta_modularity"] = sweep_df["modularity_best"].diff().fillna(0)

    max_mod = sweep_df["modularity_best"].max()
    sweep_df["delta_normalized"] = (
        sweep_df["delta_modularity"] / (max_mod + 1e-10)
    )

    sweep_csv = os.path.join(cancer_dir, f"{cancer_type}_resolution_sweep.csv")
    sweep_df.to_csv(sweep_csv, index=False)
    log(f"[SAVE] Sweep table -> {sweep_csv}")

    # Log the full sweep table
    log(f"\n[STEP 15.8] Resolution sweep summary:")
    log(f"{'Res':>5s} | {'Mod':>6s} | {'ARI':>5s} | {'Even':>5s} | "
        f"{'Comp':>6s} | {'dMod':>6s} | {'Comms':>5s} | {'Sat':>3s} | "
        f"{'Largest':>7s}")
    log("-" * 72)
    for _, row in sweep_df.iterrows():
        log(f"{row['resolution']:5.2f} | {row['modularity_best']:6.4f} | "
            f"{row['ARI_mean']:5.3f} | {row['evenness']:5.3f} | "
            f"{row['composite_score']:6.4f} | "
            f"{row['delta_modularity']:+6.4f} | "
            f"{row['ncomms_mean']:5.1f} | {row['n_satellite']:3.0f} | "
            f"{row['largest_community']:7.0f}")

    # =========================================================================
    # Diagnostic sweep plots (6 panels)
    # =========================================================================
    fig, axes = plt.subplots(2, 3, figsize=(22, 12))

    ax = axes[0, 0]
    ax.errorbar(sweep_df["resolution"], sweep_df["modularity_best"],
                yerr=sweep_df["modularity_std"], marker="o", capsize=3,
                color="#4682B4")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Modularity", fontsize=12)
    ax.set_title(f"{cancer_type} | Modularity", fontsize=13)

    ax = axes[0, 1]
    ax.errorbar(sweep_df["resolution"], sweep_df["ARI_mean"],
                yerr=sweep_df["ARI_std"], marker="o", capsize=3,
                color="#B22222")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("ARI (stability)", fontsize=12)
    ax.set_title(f"{cancer_type} | Stability (ARI)", fontsize=13)

    ax = axes[0, 2]
    ax.plot(sweep_df["resolution"], sweep_df["evenness"], "o-",
            color="#228B22")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Evenness", fontsize=12)
    ax.set_title(f"{cancer_type} | Evenness (1 - largest/total)", fontsize=13)
    ax.set_ylim(-0.05, 1.05)

    ax = axes[1, 0]
    ax.plot(sweep_df["resolution"], sweep_df["composite_score"], "o-",
            color="#E8913A", linewidth=2, markersize=8)
    best_idx = sweep_df["composite_score"].idxmax()
    ax.scatter([sweep_df.loc[best_idx, "resolution"]],
               [sweep_df.loc[best_idx, "composite_score"]],
               s=200, c="red", zorder=5, edgecolors="black", linewidths=2)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Composite (Mod x ARI x Even)", fontsize=12)
    ax.set_title(f"{cancer_type} | Composite Score", fontsize=13,
                 fontweight="bold")

    ax = axes[1, 1]
    ax.bar(sweep_df["resolution"], sweep_df["delta_modularity"],
           width=0.08, color="#4682B4", edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Delta Modularity", fontsize=12)
    ax.set_title(f"{cancer_type} | Delta-Modularity (Elbow)", fontsize=13)

    ax = axes[1, 2]
    ax.errorbar(sweep_df["resolution"], sweep_df["ncomms_mean"],
                yerr=sweep_df["ncomms_std"], marker="o", capsize=3,
                color="#800080")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("N Communities (raw Leiden)", fontsize=12)
    ax.set_title(f"{cancer_type} | Community Count", fontsize=13)

    plt.suptitle(
        f"{cancer_type} -- Resolution Selection Diagnostics "
        f"(threshold={active_threshold}, full network)",
        fontsize=16, fontweight="bold"
    )
    plt.tight_layout()
    sweep_plot = os.path.join(
        cancer_dir, f"{cancer_type}_sweep_diagnostics.png"
    )
    plt.savefig(sweep_plot, dpi=300)
    plt.close()
    log(f"[SAVE] Sweep diagnostics -> {sweep_plot}")

    # =========================================================================
    # Select best resolution using composite score
    # =========================================================================
    best_row = sweep_df.loc[sweep_df["composite_score"].idxmax()]
    best_res = best_row["resolution"]

    log(f"\n[STEP 15] SELECTED RESOLUTION: {best_res:.2f}")
    log(f"  Composite score: {best_row['composite_score']:.4f}")
    log(f"  Modularity: {best_row['modularity_best']:.4f}")
    log(f"  ARI: {best_row['ARI_mean']:.3f}")
    log(f"  Evenness: {best_row['evenness']:.3f}")
    log(f"  Communities (raw): {best_row['ncomms_mean']:.1f}")
    log(f"  Communities (merged): {best_row['ncomms_merged']:.0f} "
        f"({best_row['n_satellite']:.0f} satellite)")
    log(f"  Largest community: {best_row['largest_community']:.0f} genes")

    # Copy best partition as the "official" output
    best_res_tag = f"res{best_res:.2f}".replace(".", "p")
    best_merged_csv = os.path.join(
        sweep_dir, best_res_tag,
        f"{cancer_type}_partition_merged_{best_res_tag}.csv"
    )
    if os.path.exists(best_merged_csv):
        final_part = os.path.join(
            cancer_dir, f"{cancer_type}_best_partition.csv"
        )
        shutil.copy2(best_merged_csv, final_part)
        log(f"[SAVE] Best partition -> {final_part}")

        best_gene_csv = os.path.join(
            sweep_dir, best_res_tag,
            f"{cancer_type}_gene_lists_{best_res_tag}.csv"
        )
        if os.path.exists(best_gene_csv):
            final_genes = os.path.join(
                cancer_dir, f"{cancer_type}_community_gene_lists.csv"
            )
            shutil.copy2(best_gene_csv, final_genes)
            log(f"[SAVE] Gene lists -> {final_genes}")

    # Save the community graph (full network, not just LCC)
    comm_pkl = os.path.join(cancer_dir, f"{cancer_type}_G_comm.gpickle")
    with open(comm_pkl, "wb") as f:
        pickle.dump(G_comm, f)
    log(f"[SAVE] Community graph -> {comm_pkl}")

    # =========================================================================
    # Community summary report
    # =========================================================================
    banner("[STEP 15.9] Community summary report")

    # Reload the best merged partition for reporting
    if os.path.exists(os.path.join(
            cancer_dir, f"{cancer_type}_best_partition.csv")):
        best_df = pd.read_csv(os.path.join(
            cancer_dir, f"{cancer_type}_best_partition.csv"
        ))

        summary_path = os.path.join(
            cancer_dir, f"{cancer_type}_community_summary.txt"
        )
        with open(summary_path, "w") as sf:
            sf.write(f"Community Summary: {cancer_type}\n")
            sf.write(f"DIFF threshold: {active_threshold}\n")
            sf.write(f"Leiden resolution: {best_res}\n")
            sf.write(f"Clustering scope: full network "
                     f"(all {n_components} connected components)\n")
            sf.write(f"Total genes: {len(best_df)}\n\n")

            for cid in sorted(best_df["community"].unique()):
                comm_genes = best_df[best_df["community"] == cid]
                symbols = sorted(comm_genes["gene_symbol"].tolist())

                # Check for A3 genes
                a3_here = []
                for g in a3_seeds:
                    alias = A3_ID_TO_ALIAS.get(g, g)
                    if g in comm_genes["gene"].values:
                        a3_here.append(alias)

                a3_tag = f"  [A3: {', '.join(a3_here)}]" if a3_here else ""
                sf.write(f"Community {cid} (n={len(comm_genes)}){a3_tag}\n")
                # Write gene symbols in rows of 10
                for i in range(0, len(symbols), 10):
                    sf.write(f"  {', '.join(symbols[i:i+10])}\n")
                sf.write("\n")

        log(f"[SAVE] Community summary -> {summary_path}")

    # Save selected parameters
    params_path = os.path.join(
        cancer_dir, f"{cancer_type}_selected_parameters.txt"
    )
    with open(params_path, "w") as f:
        f.write(f"DIFF_THRESHOLD={active_threshold}\n")
        f.write(f"THRESHOLD_METHOD=max_fragmentation_rate\n")
        f.write(f"CLUSTERING_SCOPE=full_network\n")
        f.write(f"N_COMPONENTS={n_components}\n")
        f.write(f"LEIDEN_RESOLUTION={best_res}\n")
        f.write(f"RESOLUTION_METHOD=composite_score_mod_x_ari_x_evenness\n")
        f.write(f"MODULARITY={best_row['modularity_best']:.4f}\n")
        f.write(f"ARI={best_row['ARI_mean']:.4f}\n")
        f.write(f"EVENNESS={best_row['evenness']:.4f}\n")
        f.write(f"COMPOSITE_SCORE={best_row['composite_score']:.4f}\n")
        f.write(f"A3_SEEDS={','.join(a3_seeds)}\n")
    log(f"[SAVE] Selected parameters -> {params_path}")

    log(f"\n[STEP 15 COMPLETE for {cancer_type}]")


banner("[STEP 05 COMPLETE -- ALL CANCER TYPES]")
