#!/usr/bin/env python3
"""
Step03_SC_Community_Detection.py
=================================

Figure 4 -- Step 03: Leiden community detection on the single-cell DIFF
network with unified threshold and resolution selection.

Mirrors Step05_Community_Detection.py from the Figure 2 (TCGA bulk) pipeline.
Uses the same threshold selection criterion and resolution selection logic
so that both TCGA and single-cell networks are built identically.

IMPORTANT: This script builds the DIFF graph from the saved correlation
matrix pickle. The DIFF threshold is auto-selected using the unified
criterion: among thresholds where A3A and/or A3B (those passing DE at
raw p < 0.05 without force-keeping) have degree >= 1, select the
threshold with the highest connected component count. Ties broken by
the lower (more inclusive) threshold.

Leiden resolution is selected using composite score = modularity x ARI x
evenness, matching TCGA Step05 exactly.

Workflow:
  1. Load DIFF correlation matrix from Step 02 pickle
  2. Auto-select DIFF threshold (A3 connectivity + component peak)
  3. Build graph at selected threshold, remove isolates
  4. Extract largest connected component (LCC)
  5. Run Leiden at multiple resolutions (15 runs each for stability)
  6. Select best resolution via composite score (Mod x ARI x Evenness)
  7. Merge communities smaller than MIN_COMMUNITY_SIZE
  8. Save community assignments, gene lists, and annotated graph

Input:
  data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl
  data/FIG_4/02_differential_expression/SC_DE_results.csv  (for A3 seed ID)

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
from datetime import datetime

from network_config_SC import (
    DIR_02_DE, DIR_03_NETWORKS, DIR_04_COMMUNITIES,
    DIFF_THRESHOLD, SWEEP_THRESHOLDS,
    A3_GENES, A3_ID_TO_ALIAS, A3_SYMBOL_TO_ALIAS, BIOMARKERS,
    COMMUNITY_METHOD, COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION,
    COMMUNITY_BASE_SEED, USE_LARGEST_COMPONENT,
    TARGET_BIG_COMMUNITIES, MIN_COMMUNITY_SIZE,
    RAW_P_THRESHOLD,
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
# UNIFIED THRESHOLD SELECTION
# =============================================================================

def auto_select_threshold(corr_diff, a3_seeds, sweep_thresholds):
    """
    Unified DIFF threshold selection (identical logic for TCGA and SC).

    Among all thresholds where every A3 seed gene has degree >= 1,
    select the threshold with the highest connected component count.
    Ties broken by the lower (more inclusive) threshold.

    Parameters
    ----------
    corr_diff : pd.DataFrame
        DIFF correlation matrix (genes x genes).
    a3_seeds : list of str
        Gene IDs (symbols or ENSG) that passed DE. Only A3A and/or A3B.
        The caller filters to whichever passed DE at raw p < 0.05.
    sweep_thresholds : list of float
        Thresholds to evaluate (e.g. [0.30, 0.35, ..., 0.90]).

    Returns
    -------
    selected_threshold : float or None
        The selected threshold, or None if no threshold satisfies the criterion.
    sweep_df : pd.DataFrame
        Diagnostic table with per-threshold metrics.
    """
    matrix_genes = set(corr_diff.index)

    # Verify which A3 seeds are actually in the matrix
    seeds_in_matrix = [g for g in a3_seeds if g in matrix_genes]
    seeds_missing = [g for g in a3_seeds if g not in matrix_genes]

    log(f"  A3 seed genes for threshold selection:")
    for g in seeds_in_matrix:
        alias = A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))
        log(f"    {g} ({alias}) -- IN matrix")
    for g in seeds_missing:
        alias = A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))
        log(f"    {g} ({alias}) -- MISSING from matrix")

    if not seeds_in_matrix:
        log("  WARNING: No A3 seed genes found in DIFF matrix. Cannot auto-select.")
        return None, pd.DataFrame()

    # Sweep thresholds
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
            **{f"deg_{A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))}": d
               for g, d in seed_degrees.items()},
        })

    sweep_df = pd.DataFrame(rows)

    # Log the sweep table
    log(f"\n  {'Thresh':>7s}  {'Nodes':>6s}  {'Edges':>7s}  {'Comp':>5s}  "
        f"{'LCC':>5s}  {'A3 OK':>5s}  ", end="")
    for g in seeds_in_matrix:
        alias = A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))
        log(f"{'deg_' + alias:>8s}  ", end="")
    log("")
    log(f"  {'-----':>7s}  {'-----':>6s}  {'-----':>7s}  {'----':>5s}  "
        f"{'---':>5s}  {'-----':>5s}  ", end="")
    for _ in seeds_in_matrix:
        log(f"{'--------':>8s}  ", end="")
    log("")

    for _, row in sweep_df.iterrows():
        marker = " <--" if row["all_seeds_connected"] else ""
        log(f"  {row['threshold']:>7.2f}  {row['nodes']:>6.0f}  "
            f"{row['edges']:>7.0f}  {row['components']:>5.0f}  "
            f"{row['lcc_size']:>5.0f}  {'YES' if row['all_seeds_connected'] else 'no':>5s}  ",
            end="")
        for g in seeds_in_matrix:
            alias = A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))
            col = f"deg_{alias}"
            log(f"{row[col]:>8.0f}  ", end="")
        log(marker)

    # Filter to A3-valid thresholds
    valid = sweep_df[sweep_df["all_seeds_connected"]].copy()

    if len(valid) == 0:
        log("\n  WARNING: No threshold keeps all A3 seeds connected (degree >= 1).")
        log("  Falling back to config DIFF_THRESHOLD.")
        return None, sweep_df

    # Among valid thresholds, find peak component count
    max_comp = valid["components"].max()
    peak_candidates = valid[valid["components"] == max_comp]

    # Ties broken by lower (more inclusive) threshold
    selected = float(peak_candidates["threshold"].min())

    sel_row = sweep_df[sweep_df["threshold"] == selected].iloc[0]
    log(f"\n  SELECTED threshold: {selected:.2f}")
    log(f"    Reason: peak components ({int(sel_row['components'])}) among "
        f"A3-valid thresholds")
    log(f"    Nodes: {int(sel_row['nodes'])}, Edges: {int(sel_row['edges'])}, "
        f"LCC: {int(sel_row['lcc_size'])}")
    for g in seeds_in_matrix:
        alias = A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))
        col = f"deg_{alias}"
        log(f"    {alias} degree: {int(sel_row[col])}")

    return selected, sweep_df


# =============================================================================
# A3 SEED IDENTIFICATION
# =============================================================================

def identify_a3_seeds(de_dir):
    """
    Identify which of A3A and A3B pass DE at raw p < 0.05 without force-keeping.

    Scans DE result files in de_dir for A3A and A3B. Returns list of gene IDs
    (symbols or ENSG depending on what's in the DE file) that passed the
    raw p-value threshold.

    Parameters
    ----------
    de_dir : str
        Path to DE results directory (e.g. data/FIG_4/02_differential_expression/).

    Returns
    -------
    seeds : list of str
        Gene IDs for A3A and/or A3B that passed DE.
    """
    # Look for DE results file
    de_candidates = [
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
    sep = "\t" if de_path.endswith(".tsv") else ","
    de_df = pd.read_csv(de_path, sep=sep)

    # Identify gene column (could be 'gene', 'gene_symbol', 'Gene', etc.)
    gene_col = None
    for col in ["gene", "gene_symbol", "Gene", "gene_id"]:
        if col in de_df.columns:
            gene_col = col
            break
    if gene_col is None:
        log(f"  WARNING: Cannot identify gene column in DE file. "
            f"Columns: {list(de_df.columns)}")
        return ["APOBEC3A", "APOBEC3B"]

    # Identify p-value column
    pval_col = None
    for col in ["pvalue", "p_value", "pval", "raw_p", "PValue"]:
        if col in de_df.columns:
            pval_col = col
            break
    if pval_col is None:
        log(f"  WARNING: Cannot identify p-value column in DE file. "
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
                log(f"    {label}: PASSED DE (p={pval:.2e}, gene_id={gene_id})")
                seeds.append(gene_id)
            else:
                log(f"    {label}: FAILED DE (p={pval:.2e})")
        else:
            log(f"    {label}: not found in DE results")

    if not seeds:
        log(f"  WARNING: Neither A3A nor A3B passed DE. Cannot anchor threshold.")

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
            return comm_map.copy()

    big_set = set(big_comms)
    new_map = {}

    for n, c in comm_map.items():
        if c in big_set:
            new_map[n] = c
        else:
            # Find nearest big community by total edge weight
            best_comm = None
            best_weight = -1
            for neighbor in G.neighbors(n):
                nc = comm_map.get(neighbor)
                if nc in big_set:
                    w = abs(float(G[n][neighbor].get(weight, G[n][neighbor].get("weight", 1.0))))
                    if w > best_weight:
                        best_weight = w
                        best_comm = nc
            if best_comm is not None:
                new_map[n] = best_comm
            else:
                # No direct edge to a big community -- assign to largest
                new_map[n] = big_comms[0]

    n_merged = sum(1 for n in comm_map if comm_map[n] not in big_set and n in new_map)
    log(f"  Merged {n_merged} genes from small communities into {len(big_comms)} large communities")

    return new_map


def remove_intra_community_isolates(comm_map, G):
    """
    Remove nodes that have zero edges to other members of their own community.
    These are 'orphan' nodes created by merging -- they were assigned to a community
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
    # 1. Identify A3 seed genes from DE results
    # =========================================================================
    banner("[STEP 03.1] Identify A3 seed genes from DE results")

    a3_seeds = identify_a3_seeds(DIR_02_DE)

    if not a3_seeds:
        log("FATAL: No A3 seed genes available. Cannot proceed with threshold selection.")
        return

    log(f"\n  A3 seeds for threshold selection: {a3_seeds}")

    # =========================================================================
    # 2. Load DIFF correlation matrix and auto-select threshold
    # =========================================================================
    banner("[STEP 03.2] Load DIFF correlation matrix and auto-select threshold")

    corr_path = os.path.join(DIR_03_NETWORKS, "corr_matrices", "SC_corr_DIFF.pkl")
    log(f"Loading: {corr_path}")
    with open(corr_path, "rb") as f:
        corr_diff = pickle.load(f)
    log(f"  DIFF correlation matrix: {corr_diff.shape}")

    # Auto-select threshold
    selected_threshold, threshold_sweep_df = auto_select_threshold(
        corr_diff, a3_seeds, SWEEP_THRESHOLDS
    )

    if selected_threshold is not None:
        active_threshold = selected_threshold
        log(f"\n  Using auto-selected threshold: {active_threshold}")
    else:
        active_threshold = DIFF_THRESHOLD
        log(f"\n  Auto-selection failed. Using config fallback: {active_threshold}")

    # Save threshold sweep diagnostics
    if len(threshold_sweep_df) > 0:
        thresh_csv = os.path.join(DIR_04_COMMUNITIES, "SC_threshold_sweep.csv")
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

    # Save the rebuilt graph for Steps 04 and 05
    graph_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "graph_objects"))
    graph_save_path = os.path.join(graph_dir, "SC_G_diff_noiso.gpickle")
    with open(graph_save_path, "wb") as f:
        pickle.dump(G_diff_noiso, f)
    log(f"  [SAVE] Updated DIFF graph -> {graph_save_path}")

    if G_diff_noiso.number_of_nodes() < 5 or G_diff_noiso.number_of_edges() < 3:
        log(f"WARNING: Graph too small for community detection "
            f"({G_diff_noiso.number_of_nodes()} nodes, {G_diff_noiso.number_of_edges()} edges)")
        log(f"  Consider relaxing the DIFF threshold (currently {active_threshold})")
        return

    # =========================================================================
    # 4. Extract LCC
    # =========================================================================
    if USE_LARGEST_COMPONENT and G_diff_noiso.number_of_nodes() > 0:
        banner("[STEP 03.4] Extract largest connected component")
        components = list(nx.connected_components(G_diff_noiso))
        lcc = max(components, key=len)
        G_comm = G_diff_noiso.subgraph(lcc).copy()
        log(f"  Components: {len(components)}")
        log(f"  LCC: {len(lcc)} nodes ({100*len(lcc)/G_diff_noiso.number_of_nodes():.1f}%)")
        log(f"  Non-LCC nodes dropped: {G_diff_noiso.number_of_nodes() - len(lcc)}")

        # Report A3 seed status in LCC
        for g in a3_seeds:
            alias = A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g))
            in_lcc = g in lcc
            deg = G_diff_noiso.degree(g) if g in G_diff_noiso else 0
            log(f"    {alias}: {'IN LCC' if in_lcc else 'NOT in LCC'} (degree={deg})")
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
    # 5. Resolution sweep with stability assessment
    # =========================================================================
    banner("[STEP 03.5] Resolution sweep")

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
            mod = compute_modularity(G_comm, cm)
            modularities.append(mod)

            if mod > best_mod:
                best_mod = mod
                best_cm = cm.copy()

        labels_list = [partition_to_labels(p, nodes) for p in partitions]

        ari_pairs = []
        nmi_pairs = []
        for i, j in itertools.combinations(range(len(labels_list)), 2):
            ari_pairs.append(adjusted_rand_score(labels_list[i], labels_list[j]))
            nmi_pairs.append(normalized_mutual_info_score(labels_list[i], labels_list[j]))

        mean_mod = np.mean(modularities)
        mean_ari = np.mean(ari_pairs)
        mean_nmi = np.mean(nmi_pairs)

        # Community count stats
        ncomms_list = [len(set(p.values())) for p in partitions]
        n_comms_mean = np.mean(ncomms_list)
        n_comms_std = np.std(ncomms_list)

        # Largest community size (from best partition)
        comm_sizes_best = {}
        for c in best_cm.values():
            comm_sizes_best[c] = comm_sizes_best.get(c, 0) + 1
        largest_comm = max(comm_sizes_best.values()) if comm_sizes_best else 0
        n_total_genes = len(best_cm)

        log(f"    Communities: {n_comms_mean:.1f} (+/- {n_comms_std:.1f})")
        log(f"    Modularity: {best_mod:.4f} (mean {mean_mod:.4f} +/- {np.std(modularities):.4f})")
        log(f"    ARI stability: {mean_ari:.4f}")
        log(f"    NMI stability: {mean_nmi:.4f}")
        log(f"    Largest community: {largest_comm} genes")

        sweep_results.append({
            "resolution": res,
            "ncomms_mean": n_comms_mean,
            "ncomms_std": n_comms_std,
            "modularity_best": best_mod,
            "modularity_mean": mean_mod,
            "modularity_std": np.std(modularities),
            "ARI_mean": mean_ari,
            "ARI_std": np.std(ari_pairs),
            "NMI_mean": mean_nmi,
            "NMI_std": np.std(nmi_pairs),
            "n_total_genes": n_total_genes,
            "largest_community": largest_comm,
        })

        all_best_partitions[res] = best_cm

        # Save per-resolution partition
        res_df = pd.DataFrame([
            {"gene": g, "community": c} for g, c in best_cm.items()
        ])
        res_path = os.path.join(sweep_dir, f"SC_partition_res{res:.1f}.csv")
        res_df.to_csv(res_path, index=False)

    sweep_df = pd.DataFrame(sweep_results)

    # =========================================================================
    # Compute composite selection metrics (matches TCGA Step05)
    # =========================================================================
    # Evenness: 1 - (largest_community / total_genes)
    # Penalizes degenerate solutions where one community absorbs everything
    sweep_df["evenness"] = 1.0 - (sweep_df["largest_community"] / sweep_df["n_total_genes"])

    # Composite score: modularity x ARI x evenness
    # - modularity rewards partition quality
    # - ARI rewards reproducibility
    # - evenness penalizes one-giant-community solutions
    sweep_df["composite_score"] = (
        sweep_df["modularity_best"] *
        sweep_df["ARI_mean"] *
        sweep_df["evenness"]
    )

    # Delta-modularity: marginal gain at each resolution step
    sweep_df["delta_modularity"] = sweep_df["modularity_best"].diff().fillna(0)

    sweep_csv = os.path.join(DIR_04_COMMUNITIES, "SC_resolution_sweep.csv")
    sweep_df.to_csv(sweep_csv, index=False)
    log(f"\n  [SAVE] Resolution sweep -> {sweep_csv}")

    # Log the full sweep table
    log(f"\n  {'Res':>5s} | {'Mod':>6s} | {'ARI':>5s} | {'Even':>5s} | "
        f"{'Comp':>6s} | {'dMod':>6s} | {'Comms':>5s} | {'Largest':>7s}")
    log(f"  " + "-" * 65)
    for _, row in sweep_df.iterrows():
        log(f"  {row['resolution']:5.2f} | {row['modularity_best']:6.4f} | "
            f"{row['ARI_mean']:5.3f} | {row['evenness']:5.3f} | "
            f"{row['composite_score']:6.4f} | {row['delta_modularity']:+6.4f} | "
            f"{row['ncomms_mean']:5.1f} | {row['largest_community']:7.0f}")

    # =========================================================================
    # 6. Select best resolution using composite score
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
    log(f"    Largest community: {best_row['largest_community']:.0f} genes")

    best_mod = float(best_row["modularity_best"])
    best_ari = float(best_row["ARI_mean"])
    best_nmi = float(best_row["NMI_mean"])

    best_partition = all_best_partitions[best_res]

    # =========================================================================
    # 7. Merge small communities
    # =========================================================================
    banner("[STEP 03.7] Merge small communities")

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
    # After merging, some communities contain multiple disconnected sub-components.
    # Keep only the largest connected component within each community.
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
    # 8. Save community assignments
    # =========================================================================
    banner("[STEP 03.8] Save community assignments")

    # Best partition CSV
    part_rows = []
    for gene, comm in sorted(final_partition.items(), key=lambda x: (x[1], x[0])):
        alias = A3_ID_TO_ALIAS.get(gene, A3_SYMBOL_TO_ALIAS.get(gene, ""))
        part_rows.append({"gene": gene, "community": comm, "a3_alias": alias})

    part_df = pd.DataFrame(part_rows)
    part_path = os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv")
    part_df.to_csv(part_path, index=False)
    log(f"  [SAVE] Best partition ({len(part_df)} genes) -> {part_path}")

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
    log(f"  [SAVE] Gene lists -> {gene_list_path}")

    # Annotated graph
    for n in G_comm.nodes():
        G_comm.nodes[n]["community"] = final_partition.get(n, -1)

    graph_path = os.path.join(DIR_04_COMMUNITIES, "SC_G_comm.gpickle")
    with open(graph_path, "wb") as f:
        pickle.dump(G_comm, f)
    log(f"  [SAVE] Annotated graph -> {graph_path}")

    # Selected parameters
    params_path = os.path.join(DIR_04_COMMUNITIES, "SC_selected_parameters.txt")
    with open(params_path, "w") as f:
        f.write(f"DIFF_THRESHOLD={active_threshold}\n")
        f.write(f"THRESHOLD_METHOD=unified_a3_connectivity_component_peak\n")
        f.write(f"LEIDEN_RESOLUTION={best_res}\n")
        f.write(f"RESOLUTION_METHOD=composite_score_mod_x_ari_x_evenness\n")
        f.write(f"N_COMMUNITIES={len(final_sizes)}\n")
        f.write(f"N_GENES={total_genes}\n")
        f.write(f"MODULARITY={best_mod:.4f}\n")
        f.write(f"ARI={best_ari:.4f}\n")
        f.write(f"NMI={best_nmi:.4f}\n")
        f.write(f"EVENNESS={float(best_row['evenness']):.4f}\n")
        f.write(f"COMPOSITE_SCORE={float(best_row['composite_score']):.4f}\n")
        f.write(f"A3_SEEDS={','.join(a3_seeds)}\n")
    log(f"  [SAVE] Selected parameters -> {params_path}")

    # =========================================================================
    # 9. Community summary
    # =========================================================================
    banner("[STEP 03.9] Community summary")

    summary_lines = []
    summary_lines.append(f"Single-Cell Community Detection Summary")
    summary_lines.append(f"=" * 60)
    summary_lines.append(f"DIFF threshold: {active_threshold} (unified A3 connectivity + component peak)")
    summary_lines.append(f"Leiden resolution: {best_res} (composite score: Mod x ARI x Evenness)")
    summary_lines.append(f"Modularity: {best_mod:.4f}")
    summary_lines.append(f"ARI stability: {best_ari:.4f}")
    summary_lines.append(f"NMI stability: {best_nmi:.4f}")
    summary_lines.append(f"Evenness: {float(best_row['evenness']):.4f}")
    summary_lines.append(f"Composite score: {float(best_row['composite_score']):.4f}")
    summary_lines.append(f"Communities: {len(final_sizes)}")
    summary_lines.append(f"Total genes: {total_genes}")
    summary_lines.append(f"A3 seeds: {', '.join(a3_seeds)}")
    summary_lines.append(f"")

    # Per-community details
    all_biomarker_genes = set()
    for markers in BIOMARKERS.values():
        all_biomarker_genes.update(markers)

    for c in sorted(final_sizes.keys()):
        genes_in_c = sorted([g for g, cc in final_partition.items() if cc == c])
        summary_lines.append(f"Community {c}: {len(genes_in_c)} genes")

        # A3 genes
        a3_in_c = [g for g in genes_in_c if g in A3_ID_TO_ALIAS or g in A3_SYMBOL_TO_ALIAS]
        if a3_in_c:
            aliases = [A3_ID_TO_ALIAS.get(g, A3_SYMBOL_TO_ALIAS.get(g, g)) for g in a3_in_c]
            summary_lines.append(f"  A3 genes: {', '.join(aliases)}")

        # Biomarkers
        bio_in_c = [g for g in genes_in_c if g in all_biomarker_genes]
        if bio_in_c:
            summary_lines.append(f"  Biomarkers: {', '.join(bio_in_c)}")

        # Top genes by degree (top 5)
        if G_comm.number_of_nodes() > 0:
            degs = [(g, G_comm.degree(g)) for g in genes_in_c if g in G_comm]
            degs.sort(key=lambda x: -x[1])
            top5 = degs[:5]
            summary_lines.append(f"  Top hubs: {', '.join(f'{g}(d={d})' for g, d in top5)}")

        summary_lines.append("")

    summary_path = os.path.join(DIR_04_COMMUNITIES, "SC_community_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines))
    log(f"  [SAVE] Community summary -> {summary_path}")

    # Print summary
    for line in summary_lines:
        log(f"  {line}")

    # =========================================================================
    # 10. Diagnostic sweep plots (6-panel, matches TCGA Step05)
    # =========================================================================
    banner("[STEP 03.10] Resolution sweep diagnostic plots")

    fig, axes = plt.subplots(2, 3, figsize=(22, 12))

    # Panel 1: Modularity vs resolution
    ax = axes[0, 0]
    ax.errorbar(sweep_df["resolution"], sweep_df["modularity_best"],
                yerr=sweep_df["modularity_std"], marker="o", capsize=3, color="steelblue")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Modularity", fontsize=12)
    ax.set_title("Modularity", fontsize=13)

    # Panel 2: ARI stability vs resolution
    ax = axes[0, 1]
    ax.errorbar(sweep_df["resolution"], sweep_df["ARI_mean"],
                yerr=sweep_df["ARI_std"], marker="o", capsize=3, color="firebrick")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("ARI (stability)", fontsize=12)
    ax.set_title("Stability (ARI)", fontsize=13)

    # Panel 3: Evenness vs resolution
    ax = axes[0, 2]
    ax.plot(sweep_df["resolution"], sweep_df["evenness"], "o-", color="forestgreen")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Evenness", fontsize=12)
    ax.set_title("Evenness (1 - largest/total)", fontsize=13)
    ax.set_ylim(-0.05, 1.05)

    # Panel 4: Composite score vs resolution (the selection criterion)
    ax = axes[1, 0]
    ax.plot(sweep_df["resolution"], sweep_df["composite_score"], "o-",
            color="darkorange", linewidth=2, markersize=8)
    best_idx = sweep_df["composite_score"].idxmax()
    ax.scatter([sweep_df.loc[best_idx, "resolution"]],
               [sweep_df.loc[best_idx, "composite_score"]],
               s=200, c="red", zorder=5, edgecolors="black", linewidths=2)
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Composite (Mod x ARI x Even)", fontsize=12)
    ax.set_title("Composite Score", fontsize=13, fontweight="bold")

    # Panel 5: Delta-modularity (elbow plot)
    ax = axes[1, 1]
    ax.bar(sweep_df["resolution"], sweep_df["delta_modularity"],
           width=0.08, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Delta Modularity", fontsize=12)
    ax.set_title("Delta-Modularity (Elbow)", fontsize=13)

    # Panel 6: Number of communities vs resolution
    ax = axes[1, 2]
    ax.errorbar(sweep_df["resolution"], sweep_df["ncomms_mean"],
                yerr=sweep_df["ncomms_std"], marker="o", capsize=3, color="purple")
    ax.axvline(best_res, ls="--", c="red", alpha=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("N Communities", fontsize=12)
    ax.set_title("Community Count", fontsize=13)

    plt.suptitle(f"SC Network -- Resolution Selection Diagnostics "
                 f"(threshold={active_threshold})", fontsize=16, fontweight="bold")
    plt.tight_layout()
    sweep_plot = os.path.join(DIR_04_COMMUNITIES, "SC_sweep_diagnostics.png")
    plt.savefig(sweep_plot, dpi=300)
    plt.close()
    log(f"  [SAVE] Sweep diagnostics -> {sweep_plot}")

    # =========================================================================
    # DONE
    # =========================================================================
    elapsed = datetime.now() - start_time
    banner(f"[STEP 03 COMPLETE] {len(final_sizes)} communities, "
           f"{total_genes} genes | threshold={active_threshold}, "
           f"resolution={best_res} | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
