#!/usr/bin/env python3
"""
Step03_SC_Community_Detection.py
=================================

Figure 4 -- Step 03: Leiden community detection on the full single-cell
DIFF network with unified threshold and resolution selection.

Mirrors Step05_Community_Detection.py from the Figure 2 (TCGA bulk)
pipeline. Uses identical threshold selection and resolution selection
logic so that both TCGA and single-cell networks are built the same way.

THRESHOLD SELECTION: Maximum Fragmentation Rate
  Sweeps DIFF correlation thresholds and computes the forward difference
  in connected component count between consecutive steps. Among intervals
  where the upper threshold is A3-valid (A3A and A3B have degree >= 1),
  selects the upper threshold of the interval with the largest positive
  delta-component count. Ties broken by the lower (more inclusive) upper
  threshold.

COMMUNITY DETECTION: Full-Network Leiden
  Leiden runs on the entire graph (all connected components), not just
  the LCC. Large components get subdivided into pathway-scale communities.
  Small satellite components (where A3 genes often reside) remain as
  their own communities.

  Merge logic is component-aware: small Leiden communities within large
  connected components are merged (Leiden artifacts), but satellite
  communities that ARE entire small components are preserved as
  biologically meaningful units.

Resolution selection uses composite score = modularity x ARI x evenness
(identical to TCGA Step05).

Workflow:
  1. Identify A3 seed genes from DE results
  2. Load DIFF correlation matrix and auto-select threshold
  3. Build graph at selected threshold, remove isolates
  4. Run Leiden on full network at multiple resolutions
  5. Select best resolution via composite score (Mod x ARI x Evenness)
  6. Component-aware merge of small communities
  7. Save community assignments, gene lists, and annotated graph

Input:
  data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl
  data/FIG_4/02_differential_expression/SC_DE_results.csv

Output (to data/FIG_4/04_communities/):
  sweep/                          -- per-resolution assignments
  SC_resolution_sweep.csv         -- stability metrics per resolution
  SC_sweep_diagnostics.png        -- 6-panel diagnostic plots
  SC_best_partition.csv           -- final gene-to-community assignments
  SC_community_gene_lists.csv     -- genes per community
  SC_community_summary.txt        -- community sizes and key genes
  SC_G_comm.gpickle               -- annotated graph with community labels
  SC_selected_parameters.txt      -- threshold + resolution + metrics
  SC_threshold_sweep.csv          -- threshold selection diagnostics

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
from networkx.algorithms.community.quality import modularity as nx_modularity
from datetime import datetime

from network_config_SC import (
    DIR_02_DE, DIR_03_NETWORKS, DIR_04_COMMUNITIES,
    DIFF_THRESHOLD, SWEEP_THRESHOLDS,
    A3_GENES, A3_ID_TO_ALIAS, A3_SYMBOL_TO_ALIAS, BIOMARKERS,
    COMMUNITY_METHOD, COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION,
    COMMUNITY_BASE_SEED, USE_LARGEST_COMPONENT,
    TARGET_BIG_COMMUNITIES, MIN_COMMUNITY_SIZE,
    RAW_P_THRESHOLD, FORCE_KEEP_A3,
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


def get_alias(g):
    """Look up A3 alias for a gene, handling both ENSG and symbol formats."""
    return A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))


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

    Fallback: if no positive fragmentation occurs among A3-valid thresholds,
    select the highest A3-valid threshold.

    Parameters
    ----------
    corr_diff : pd.DataFrame
        DIFF correlation matrix (genes x genes).
    a3_seeds : list of str
        Gene symbols or ENSG IDs for A3A and/or A3B.
    sweep_thresholds : list of float
        Thresholds to evaluate.

    Returns
    -------
    selected_threshold : float or None
        The selected threshold, or None if no A3-valid threshold exists.
    sweep_df : pd.DataFrame
        Diagnostic table with per-threshold metrics.
    """
    matrix_genes = set(corr_diff.index)

    # Verify which A3 seeds are actually in the matrix
    seeds_in_matrix = [g for g in a3_seeds if g in matrix_genes]
    seeds_missing = [g for g in a3_seeds if g not in matrix_genes]

    log(f"  A3 seed genes for threshold selection:")
    for g in seeds_in_matrix:
        log(f"    {g} ({get_alias(g)}) -- IN matrix")
    for g in seeds_missing:
        log(f"    {g} ({get_alias(g)}) -- MISSING from matrix")

    if not seeds_in_matrix:
        log("  WARNING: No A3 seed genes found in DIFF matrix. "
            "Cannot auto-select.")
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
            **{f"deg_{get_alias(g)}": d for g, d in seed_degrees.items()},
        })

    sweep_df = pd.DataFrame(rows).sort_values("threshold").reset_index(drop=True)

    # ---- Log the sweep table ----
    header = (f"\n  {'Thresh':>7s}  {'Nodes':>6s}  {'Edges':>7s}  {'Comp':>5s}  "
              f"{'LCC':>5s}  {'A3 OK':>5s}  ")
    sep_line = (f"  {'-----':>7s}  {'-----':>6s}  {'-----':>7s}  {'----':>5s}  "
                f"{'---':>5s}  {'-----':>5s}  ")
    for g in seeds_in_matrix:
        alias = get_alias(g)
        header += f"{'deg_' + alias:>8s}  "
        sep_line += f"{'--------':>8s}  "
    log(header)
    log(sep_line)

    for _, row in sweep_df.iterrows():
        marker = " <--" if row["all_seeds_connected"] else ""
        line = (f"  {row['threshold']:>7.2f}  {row['nodes']:>6.0f}  "
                f"{row['edges']:>7.0f}  {row['components']:>5.0f}  "
                f"{row['lcc_size']:>5.0f}  "
                f"{'YES' if row['all_seeds_connected'] else 'no':>5s}  ")
        for g in seeds_in_matrix:
            alias = get_alias(g)
            col = f"deg_{alias}"
            line += f"{row[col]:>8.0f}  "
        log(line + marker)

    # ---- Compute forward delta-comp ----
    log(f"\n  Forward delta-comp (fragmentation rate):")
    log(f"  {'Interval':>15s}  {'delta_comp':>6s}  {'Upper A3-valid':>14s}")
    log(f"  {'--------':>15s}  {'-----':>6s}  {'--------------':>14s}")

    candidates = []
    for i in range(len(sweep_df) - 1):
        t_lower = sweep_df.loc[i, "threshold"]
        t_upper = sweep_df.loc[i + 1, "threshold"]
        delta = int(sweep_df.loc[i + 1, "components"]
                    - sweep_df.loc[i, "components"])
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
        valid = sweep_df[sweep_df["all_seeds_connected"]]
        if len(valid) == 0:
            log("\n  WARNING: No A3-valid threshold exists. "
                "Cannot auto-select.")
            return None, sweep_df

        selected = float(valid["threshold"].max())
        reason = ("fallback: highest A3-valid threshold "
                  "(no positive fragmentation among A3-valid intervals)")
        log(f"\n  No positive fragmentation candidates found.")
    else:
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
    log(f"    Nodes: {int(sel_row['nodes'])}, "
        f"Edges: {int(sel_row['edges'])}, "
        f"LCC: {int(sel_row['lcc_size'])}")
    for g in seeds_in_matrix:
        alias = get_alias(g)
        col = f"deg_{alias}"
        log(f"    {alias} degree: {int(sel_row[col])}")

    return selected, sweep_df


# =============================================================================
# A3 SEED IDENTIFICATION
# =============================================================================

def identify_a3_seeds(de_dir, force_keep=False):
    """
    Identify which of A3A and A3B should be used as seed genes for
    threshold selection.

    Two modes:
      - force_keep=False: only A3 genes passing DE are seeds
      - force_keep=True: A3A and A3B are seeds regardless of p-value

    Parameters
    ----------
    de_dir : str
        Path to DE results directory.
    force_keep : bool
        If True, include A3A/A3B as seeds even if they fail DE.

    Returns
    -------
    seeds : list of str
        Gene IDs for A3A and/or A3B to use as threshold anchors.
    """
    de_candidates = [
        os.path.join(de_dir, "SC_diffexpr_stats.csv"),
        os.path.join(de_dir, "SC_diffexpr_stats.tsv"),
        os.path.join(de_dir, "SC_DE_results.csv"),
        os.path.join(de_dir, "SC_DE_results.tsv"),
    ]

    de_path = None
    for p in de_candidates:
        if os.path.exists(p):
            de_path = p
            break

    if de_path is None:
        log(f"  WARNING: No DE results file found in {de_dir}")
        log(f"  Defaulting to both A3A and A3B as seeds")
        return ["APOBEC3A", "APOBEC3B"]

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
        return ["APOBEC3A", "APOBEC3B"]

    # Identify p-value column
    pval_col = None
    for col in ["pvalue", "p_value", "pval", "raw_p", "PValue",
                 "pvals_adj", "pvals"]:
        if col in de_df.columns:
            pval_col = col
            break
    if pval_col is None:
        log(f"  WARNING: Cannot identify p-value column. "
            f"Columns: {list(de_df.columns)}")
        return ["APOBEC3A", "APOBEC3B"]

    # Check A3A and A3B
    a3a_names = ["APOBEC3A", "A3A", "ENSG00000128383"]
    a3b_names = ["APOBEC3B", "A3B", "ENSG00000179750"]

    seeds = []
    for target_names, label in [(a3a_names, "A3A"), (a3b_names, "A3B")]:
        match = de_df[de_df[gene_col].isin(target_names)]
        if len(match) > 0:
            row = match.iloc[0]
            pval = float(row[pval_col])
            gene_id = str(row[gene_col])
            if pval < RAW_P_THRESHOLD:
                log(f"    {label}: PASSED DE (p={pval:.2e}, "
                    f"gene_id={gene_id})")
                seeds.append(gene_id)
            elif force_keep:
                log(f"    {label}: FAILED DE (p={pval:.2e}) "
                    f"-- included via FORCE_KEEP")
                seeds.append(gene_id)
            else:
                log(f"    {label}: FAILED DE (p={pval:.2e})")
        else:
            log(f"    {label}: not found in DE results")

    if not seeds:
        log(f"  WARNING: Neither A3A nor A3B available as seeds.")

    return seeds


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
    return np.array([comm_map.get(n, -1) for n in nodes], dtype=int)


def merge_small_communities(comm_map, G, k_keep=8, min_size=10):
    """
    Component-aware community merge for full-network clustering.

    Within large connected components, small Leiden communities are merged
    into an 'Other' bin. Satellite communities whose members all belong to
    a small connected component are preserved as-is (biologically meaningful
    co-expression units naturally disconnected from the main network).

    Parameters
    ----------
    comm_map : dict
        {node: community_id} from Leiden.
    G : nx.Graph
        The full graph used for community detection.
    k_keep : int
        Max number of large communities to retain within main components.
    min_size : int
        Minimum community size within large components.

    Returns
    -------
    merged : dict
        {node: community_id} after merge.
    raw_sizes : dict
        {community_id: size} before merge.
    merged_sizes : dict
        {community_id: size} after merge.
    satellite_comm_ids : set
        Community IDs preserved as satellites.
    """
    # Map each node to its connected component size
    components = list(nx.connected_components(G))
    node_to_comp_size = {}
    for comp in components:
        sz = len(comp)
        for node in comp:
            node_to_comp_size[node] = sz

    # Raw community sizes
    raw_sizes = {}
    for n, cid in comm_map.items():
        raw_sizes[cid] = raw_sizes.get(cid, 0) + 1

    # Classify: satellite vs main-component
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

    other_id = max(raw_sizes.keys()) + 1 if raw_sizes else 0

    # Build merged map
    merged = {}
    for n, cid in comm_map.items():
        if cid in satellite_comm_ids:
            merged[n] = cid
        elif cid in top_main and raw_sizes[cid] >= min_size:
            merged[n] = cid
        else:
            merged[n] = other_id

    # Merged sizes
    merged_sizes = {}
    for n, cid in merged.items():
        merged_sizes[cid] = merged_sizes.get(cid, 0) + 1

    return merged, raw_sizes, merged_sizes, satellite_comm_ids


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
    banner("[STEP 03] Single-Cell Community Detection "
           "(Full-Network Leiden)")
    log(f"Start time: {start_time}")

    ensure_dir(DIR_04_COMMUNITIES)
    sweep_dir = ensure_dir(os.path.join(DIR_04_COMMUNITIES, "sweep"))

    # =========================================================================
    # 1. Identify A3 seed genes from DE results
    # =========================================================================
    banner("[STEP 03.1] Identify A3 seed genes from DE results")

    a3_seeds = identify_a3_seeds(DIR_02_DE, force_keep=FORCE_KEEP_A3)

    if not a3_seeds:
        log("FATAL: No A3 seed genes available. "
            "Cannot proceed with threshold selection.")
        return

    log(f"\n  A3 seeds for threshold selection: {a3_seeds}")

    # =========================================================================
    # 2. Load DIFF correlation matrix and auto-select threshold
    # =========================================================================
    banner("[STEP 03.2] Load DIFF matrix and auto-select threshold "
           "(max fragmentation rate)")

    corr_path = os.path.join(DIR_03_NETWORKS, "corr_matrices",
                             "SC_corr_DIFF.pkl")
    log(f"Loading: {corr_path}")
    with open(corr_path, "rb") as f:
        corr_diff = pickle.load(f)
    log(f"  DIFF correlation matrix: {corr_diff.shape}")

    selected_threshold, threshold_sweep_df = auto_select_threshold(
        corr_diff, a3_seeds, SWEEP_THRESHOLDS
    )

    if selected_threshold is not None:
        active_threshold = selected_threshold
        log(f"\n  Using auto-selected threshold: {active_threshold}")
    else:
        active_threshold = DIFF_THRESHOLD
        log(f"\n  Auto-selection failed. "
            f"Using config fallback: {active_threshold}")

    # Save threshold sweep diagnostics
    if len(threshold_sweep_df) > 0:
        thresh_csv = os.path.join(DIR_04_COMMUNITIES,
                                  "SC_threshold_sweep.csv")
        threshold_sweep_df.to_csv(thresh_csv, index=False)
        log(f"  [SAVE] Threshold sweep diagnostics -> {thresh_csv}")

    # =========================================================================
    # 3. Build graph at selected threshold
    # =========================================================================
    banner("[STEP 03.3] Build DIFF graph at selected threshold")

    log(f"  Building graph at |delta-rho| >= {active_threshold}")
    G_diff = build_weighted_graph(corr_diff, active_threshold)
    G_diff_noiso = remove_isolated(G_diff)

    n_total = G_diff.number_of_nodes()
    n_iso = n_total - G_diff_noiso.number_of_nodes()
    log(f"  Total genes in matrix: {n_total}")
    log(f"  Non-isolated nodes: {G_diff_noiso.number_of_nodes()}")
    log(f"  Edges: {G_diff_noiso.number_of_edges()}")
    log(f"  Isolated nodes removed: {n_iso}")

    # Free correlation matrix memory
    del corr_diff

    # Save rebuilt graph
    graph_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "graph_objects"))
    graph_save_path = os.path.join(graph_dir, "SC_G_diff_noiso.gpickle")
    with open(graph_save_path, "wb") as f:
        pickle.dump(G_diff_noiso, f)
    log(f"  [SAVE] Updated DIFF graph -> {graph_save_path}")

    if (G_diff_noiso.number_of_nodes() < 5
            or G_diff_noiso.number_of_edges() < 3):
        log(f"WARNING: Graph too small for community detection "
            f"({G_diff_noiso.number_of_nodes()} nodes, "
            f"{G_diff_noiso.number_of_edges()} edges)")
        return

    # =========================================================================
    # 4. Prepare full network for community detection
    # =========================================================================
    banner("[STEP 03.4] Prepare full network for community detection")

    G_comm = G_diff_noiso.copy()

    # Ensure abs_weight on all edges
    for u, v, d in G_comm.edges(data=True):
        d["abs_weight"] = abs(float(
            d.get("abs_weight", d.get("weight", 1.0))
        ))

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
        small_genes = sum(s for s in comp_sizes if s < MIN_COMMUNITY_SIZE)
        log(f"    Components < {MIN_COMMUNITY_SIZE} nodes: "
            f"{small_count} ({small_genes} genes)")

    # Report A3 seed positions
    log(f"\n  A3 seed positions in full network:")
    for g in a3_seeds:
        alias = get_alias(g)
        if g in G_comm:
            deg = G_comm.degree(g)
            for i, comp in enumerate(components_sorted):
                if g in comp:
                    log(f"    {alias}: component {i} "
                        f"(size={len(comp)}), degree={deg}")
                    break
        else:
            log(f"    {alias}: NOT in graph")

    nodes_list = sorted(G_comm.nodes())

    # =========================================================================
    # 5. Resolution sweep with stability assessment
    # =========================================================================
    banner("[STEP 03.5] Resolution sweep (full-network Leiden)")

    sweep_results = []
    all_best_partitions = {}

    for res in COMMUNITY_RESOLUTIONS:
        log(f"\n  Resolution {res}:")

        partitions = []
        modularities = []
        best_cm = None
        best_mod = -999

        for r in range(RUNS_PER_RESOLUTION):
            seed = COMMUNITY_BASE_SEED + r
            cm = detect_leiden(G_comm, seed=seed, resolution=res)
            partitions.append(cm)

            comm_sets = {}
            for n, c in cm.items():
                comm_sets.setdefault(c, set()).add(n)
            mod = nx_modularity(G_comm, list(comm_sets.values()),
                                weight="abs_weight")
            modularities.append(mod)

            if mod > best_mod:
                best_mod = mod
                best_cm = cm.copy()

        labels_list = [partition_to_labels(p, nodes_list)
                       for p in partitions]

        ari_pairs = []
        nmi_pairs = []
        for i, j in itertools.combinations(range(len(labels_list)), 2):
            ari_pairs.append(
                adjusted_rand_score(labels_list[i], labels_list[j]))
            nmi_pairs.append(
                normalized_mutual_info_score(labels_list[i], labels_list[j]))

        # Community count stats
        ncomms_list = [len(set(p.values())) for p in partitions]

        # Largest community size (from best partition)
        comm_sizes_best = {}
        for c in best_cm.values():
            comm_sizes_best[c] = comm_sizes_best.get(c, 0) + 1
        largest_comm = max(comm_sizes_best.values()) if comm_sizes_best else 0
        n_total_genes = len(best_cm)

        log(f"    Communities: {np.mean(ncomms_list):.1f} "
            f"(+/- {np.std(ncomms_list):.1f})")
        log(f"    Modularity: {best_mod:.4f} "
            f"(mean {np.mean(modularities):.4f})")
        log(f"    ARI stability: {np.mean(ari_pairs):.4f}")
        log(f"    Largest community: {largest_comm} genes")

        sweep_results.append({
            "resolution": res,
            "ncomms_mean": float(np.mean(ncomms_list)),
            "ncomms_std": float(np.std(ncomms_list)),
            "modularity_best": float(best_mod),
            "modularity_mean": float(np.mean(modularities)),
            "modularity_std": float(np.std(modularities)),
            "ARI_mean": float(np.mean(ari_pairs)),
            "ARI_std": float(np.std(ari_pairs)),
            "NMI_mean": float(np.mean(nmi_pairs)),
            "NMI_std": float(np.std(nmi_pairs)),
            "n_total_genes": n_total_genes,
            "largest_community": largest_comm,
        })

        all_best_partitions[res] = best_cm

        # Save per-resolution raw partition
        res_df = pd.DataFrame([
            {"gene": g, "community": c} for g, c in best_cm.items()
        ])
        res_path = os.path.join(sweep_dir,
                                f"SC_partition_res{res:.1f}.csv")
        res_df.to_csv(res_path, index=False)

    sweep_df = pd.DataFrame(sweep_results)

    # =========================================================================
    # Composite selection metrics (matches TCGA Step05)
    # =========================================================================
    sweep_df["evenness"] = (
        1.0 - (sweep_df["largest_community"] / sweep_df["n_total_genes"])
    )
    sweep_df["composite_score"] = (
        sweep_df["modularity_best"]
        * sweep_df["ARI_mean"]
        * sweep_df["evenness"]
    )
    sweep_df["delta_modularity"] = (
        sweep_df["modularity_best"].diff().fillna(0)
    )

    sweep_csv = os.path.join(DIR_04_COMMUNITIES,
                             "SC_resolution_sweep.csv")
    sweep_df.to_csv(sweep_csv, index=False)
    log(f"\n  [SAVE] Resolution sweep -> {sweep_csv}")

    # Log full sweep table
    log(f"\n  {'Res':>5s} | {'Mod':>6s} | {'ARI':>5s} | {'Even':>5s} | "
        f"{'Comp':>6s} | {'dMod':>6s} | {'Comms':>5s} | "
        f"{'Largest':>7s}")
    log(f"  " + "-" * 65)
    for _, row in sweep_df.iterrows():
        log(f"  {row['resolution']:5.2f} | "
            f"{row['modularity_best']:6.4f} | "
            f"{row['ARI_mean']:5.3f} | {row['evenness']:5.3f} | "
            f"{row['composite_score']:6.4f} | "
            f"{row['delta_modularity']:+6.4f} | "
            f"{row['ncomms_mean']:5.1f} | "
            f"{row['largest_community']:7.0f}")

    # =========================================================================
    # 6. Select best resolution
    # =========================================================================
    banner("[STEP 03.6] Select best resolution (composite score)")

    best_row = sweep_df.loc[sweep_df["composite_score"].idxmax()]
    best_res = float(best_row["resolution"])

    log(f"  SELECTED RESOLUTION: {best_res:.2f}")
    log(f"    Composite score: {best_row['composite_score']:.4f}")
    log(f"    Modularity: {best_row['modularity_best']:.4f}")
    log(f"    ARI: {best_row['ARI_mean']:.3f}")
    log(f"    Evenness: {best_row['evenness']:.3f}")
    log(f"    Communities: {best_row['ncomms_mean']:.1f}")
    log(f"    Largest community: "
        f"{best_row['largest_community']:.0f} genes")

    best_mod_val = float(best_row["modularity_best"])
    best_ari_val = float(best_row["ARI_mean"])
    best_nmi_val = float(best_row["NMI_mean"])
    best_partition = all_best_partitions[best_res]

    # =========================================================================
    # 7. Component-aware merge
    # =========================================================================
    banner("[STEP 03.7] Component-aware merge")

    # Log raw partition before merge
    raw_sizes = {}
    for c in best_partition.values():
        raw_sizes[c] = raw_sizes.get(c, 0) + 1

    log(f"  Before merging: {len(raw_sizes)} raw Leiden communities")
    for c, s in sorted(raw_sizes.items(), key=lambda x: -x[1])[:15]:
        log(f"    Community {c}: {s} genes")
    if len(raw_sizes) > 15:
        log(f"    ... and {len(raw_sizes) - 15} more")

    merged, pre_sizes, merged_sizes, satellite_ids = \
        merge_small_communities(
            best_partition, G_comm,
            k_keep=TARGET_BIG_COMMUNITIES,
            min_size=MIN_COMMUNITY_SIZE
        )

    # Renumber by descending size for clean output
    final_partition = renumber_communities(merged)

    # Rebuild satellite ID set after renumbering
    # Map old IDs to new IDs
    old_to_new = {}
    for n in merged:
        old_id = merged[n]
        new_id = final_partition[n]
        old_to_new[old_id] = new_id
    final_satellite_ids = set()
    for old_id in satellite_ids:
        if old_id in old_to_new:
            final_satellite_ids.add(old_to_new[old_id])

    final_sizes = {}
    for c in final_partition.values():
        final_sizes[c] = final_sizes.get(c, 0) + 1

    n_satellite = len(final_satellite_ids)
    n_main = len(final_sizes) - n_satellite
    total_genes = sum(final_sizes.values())

    log(f"\n  After merging: {len(final_sizes)} communities "
        f"({n_main} main + {n_satellite} satellite)")
    for c, s in sorted(final_sizes.items()):
        tag = " [satellite]" if c in final_satellite_ids else ""
        log(f"    Community {c}: {s} genes{tag}")
    log(f"  Total genes in communities: {total_genes}")

    # =========================================================================
    # 8. Save community assignments
    # =========================================================================
    banner("[STEP 03.8] Save community assignments")

    # Best partition CSV
    part_rows = []
    for gene, comm in sorted(final_partition.items(),
                              key=lambda x: (x[1], x[0])):
        alias = get_alias(gene)
        is_a3 = alias != gene  # alias differs from gene if it's an A3
        part_rows.append({
            "gene": gene,
            "community": comm,
            "a3_alias": alias if is_a3 else "",
            "is_satellite": comm in final_satellite_ids,
        })

    part_df = pd.DataFrame(part_rows)
    part_path = os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv")
    part_df.to_csv(part_path, index=False)
    log(f"  [SAVE] Best partition ({len(part_df)} genes) -> {part_path}")

    # Community gene lists
    gene_list_rows = []
    for comm in sorted(final_sizes.keys()):
        genes = sorted(
            [g for g, c in final_partition.items() if c == comm]
        )
        gene_list_rows.append({
            "community": comm,
            "size": len(genes),
            "is_satellite": comm in final_satellite_ids,
            "genes": ";".join(genes)
        })

    gene_list_df = pd.DataFrame(gene_list_rows)
    gene_list_path = os.path.join(DIR_04_COMMUNITIES,
                                  "SC_community_gene_lists.csv")
    gene_list_df.to_csv(gene_list_path, index=False)
    log(f"  [SAVE] Gene lists -> {gene_list_path}")

    # Annotated graph
    for n in G_comm.nodes():
        G_comm.nodes[n]["community"] = final_partition.get(n, -1)

    graph_path = os.path.join(DIR_04_COMMUNITIES, "SC_G_comm.gpickle")
    with open(graph_path, "wb") as f:
        pickle.dump(G_comm, f)
    log(f"  [SAVE] Annotated graph -> {graph_path}")

    # Selected parameters
    params_path = os.path.join(DIR_04_COMMUNITIES,
                               "SC_selected_parameters.txt")
    with open(params_path, "w") as f:
        f.write(f"DIFF_THRESHOLD={active_threshold}\n")
        f.write(f"THRESHOLD_METHOD=max_fragmentation_rate\n")
        f.write(f"CLUSTERING_SCOPE=full_network\n")
        f.write(f"N_COMPONENTS={n_components}\n")
        f.write(f"LEIDEN_RESOLUTION={best_res}\n")
        f.write(f"RESOLUTION_METHOD="
                f"composite_score_mod_x_ari_x_evenness\n")
        f.write(f"N_COMMUNITIES={len(final_sizes)}\n")
        f.write(f"N_SATELLITE={n_satellite}\n")
        f.write(f"N_GENES={total_genes}\n")
        f.write(f"MODULARITY={best_mod_val:.4f}\n")
        f.write(f"ARI={best_ari_val:.4f}\n")
        f.write(f"NMI={best_nmi_val:.4f}\n")
        f.write(f"EVENNESS={float(best_row['evenness']):.4f}\n")
        f.write(f"COMPOSITE_SCORE="
                f"{float(best_row['composite_score']):.4f}\n")
        f.write(f"A3_SEEDS={','.join(a3_seeds)}\n")
    log(f"  [SAVE] Selected parameters -> {params_path}")

    # =========================================================================
    # 9. Community summary
    # =========================================================================
    banner("[STEP 03.9] Community summary")

    all_biomarker_genes = set()
    for markers in BIOMARKERS.values():
        all_biomarker_genes.update(markers)

    summary_lines = []
    summary_lines.append(f"Single-Cell Community Detection Summary")
    summary_lines.append(f"=" * 60)
    summary_lines.append(
        f"DIFF threshold: {active_threshold} "
        f"(max fragmentation rate)")
    summary_lines.append(
        f"Leiden resolution: {best_res} "
        f"(composite score: Mod x ARI x Evenness)")
    summary_lines.append(f"Clustering scope: full network "
                         f"(all {n_components} connected components)")
    summary_lines.append(f"Modularity: {best_mod_val:.4f}")
    summary_lines.append(f"ARI stability: {best_ari_val:.4f}")
    summary_lines.append(f"NMI stability: {best_nmi_val:.4f}")
    summary_lines.append(
        f"Evenness: {float(best_row['evenness']):.4f}")
    summary_lines.append(
        f"Composite score: "
        f"{float(best_row['composite_score']):.4f}")
    summary_lines.append(
        f"Communities: {len(final_sizes)} "
        f"({n_main} main + {n_satellite} satellite)")
    summary_lines.append(f"Total genes: {total_genes}")
    summary_lines.append(f"A3 seeds: {', '.join(a3_seeds)}")
    summary_lines.append(f"")

    for c in sorted(final_sizes.keys()):
        genes_in_c = sorted(
            [g for g, cc in final_partition.items() if cc == c]
        )
        sat_tag = " [SATELLITE]" if c in final_satellite_ids else ""
        summary_lines.append(
            f"Community {c}: {len(genes_in_c)} genes{sat_tag}")

        # A3 genes
        a3_in_c = [g for g in genes_in_c
                    if g in A3_ID_TO_ALIAS or g in A3_SYMBOL_TO_ALIAS]
        if a3_in_c:
            aliases = [get_alias(g) for g in a3_in_c]
            summary_lines.append(f"  A3 genes: {', '.join(aliases)}")

        # Biomarkers
        bio_in_c = [g for g in genes_in_c if g in all_biomarker_genes]
        if bio_in_c:
            summary_lines.append(
                f"  Biomarkers: {', '.join(bio_in_c)}")

        # Top hubs by degree (top 5)
        if G_comm.number_of_nodes() > 0:
            degs = [(g, G_comm.degree(g)) for g in genes_in_c
                    if g in G_comm]
            degs.sort(key=lambda x: -x[1])
            top5 = degs[:5]
            summary_lines.append(
                f"  Top hubs: "
                f"{', '.join(f'{g}(d={d})' for g, d in top5)}")

        summary_lines.append("")

    summary_path = os.path.join(DIR_04_COMMUNITIES,
                                "SC_community_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines))
    log(f"  [SAVE] Community summary -> {summary_path}")

    for line in summary_lines:
        log(f"  {line}")

    # =========================================================================
    # 10. Diagnostic sweep plots (6-panel, matches TCGA Step05)
    # =========================================================================
    banner("[STEP 03.10] Resolution sweep diagnostic plots")

    fig, axes = plt.subplots(2, 3, figsize=(22, 12))

    ax = axes[0, 0]
    ax.errorbar(sweep_df["resolution"], sweep_df["modularity_best"],
                yerr=sweep_df["modularity_std"], marker="o", capsize=3,
                color="#4682B4")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Modularity", fontsize=12)
    ax.set_title("Modularity", fontsize=13)

    ax = axes[0, 1]
    ax.errorbar(sweep_df["resolution"], sweep_df["ARI_mean"],
                yerr=sweep_df["ARI_std"], marker="o", capsize=3,
                color="#B22222")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("ARI (stability)", fontsize=12)
    ax.set_title("Stability (ARI)", fontsize=13)

    ax = axes[0, 2]
    ax.plot(sweep_df["resolution"], sweep_df["evenness"], "o-",
            color="#228B22")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Evenness", fontsize=12)
    ax.set_title("Evenness (1 - largest/total)", fontsize=13)
    ax.set_ylim(-0.05, 1.05)

    ax = axes[1, 0]
    ax.plot(sweep_df["resolution"], sweep_df["composite_score"], "o-",
            color="#E8913A", linewidth=2, markersize=8)
    best_idx = sweep_df["composite_score"].idxmax()
    ax.scatter([sweep_df.loc[best_idx, "resolution"]],
               [sweep_df.loc[best_idx, "composite_score"]],
               s=200, c="red", zorder=5, edgecolors="black",
               linewidths=2)
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Composite (Mod x ARI x Even)", fontsize=12)
    ax.set_title("Composite Score", fontsize=13, fontweight="bold")

    ax = axes[1, 1]
    ax.bar(sweep_df["resolution"], sweep_df["delta_modularity"],
           width=0.08, color="#4682B4", edgecolor="black",
           linewidth=0.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Delta Modularity", fontsize=12)
    ax.set_title("Delta-Modularity (Elbow)", fontsize=13)

    ax = axes[1, 2]
    ax.errorbar(sweep_df["resolution"], sweep_df["ncomms_mean"],
                yerr=sweep_df["ncomms_std"], marker="o", capsize=3,
                color="#800080")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("N Communities (raw Leiden)", fontsize=12)
    ax.set_title("Community Count", fontsize=13)

    plt.suptitle(
        f"SC Network -- Resolution Selection Diagnostics "
        f"(threshold={active_threshold}, full network)",
        fontsize=16, fontweight="bold"
    )
    plt.tight_layout()
    sweep_plot = os.path.join(DIR_04_COMMUNITIES,
                              "SC_sweep_diagnostics.png")
    plt.savefig(sweep_plot, dpi=300)
    plt.close()
    log(f"  [SAVE] Sweep diagnostics -> {sweep_plot}")

    # =========================================================================
    # DONE
    # =========================================================================
    elapsed = datetime.now() - start_time
    banner(f"[STEP 03 COMPLETE] {len(final_sizes)} communities "
           f"({n_main} main + {n_satellite} satellite), "
           f"{total_genes} genes | threshold={active_threshold}, "
           f"resolution={best_res} | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
