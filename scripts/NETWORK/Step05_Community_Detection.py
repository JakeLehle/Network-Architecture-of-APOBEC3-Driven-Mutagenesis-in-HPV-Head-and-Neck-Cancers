#!/usr/bin/env python3
"""
Step05_Community_Detection.py

Run Leiden community detection on the DIFF network across multiple
resolutions. Evaluate stability (ARI/NMI across runs), merge small
communities, and save community assignments + network visualizations.

UPDATED: Now includes unified DIFF threshold auto-selection, identical
to the single-cell Step03 criterion. Among thresholds where A3A and/or
A3B (those passing DE at raw p < 0.05 without force-keeping) have
degree >= 1, selects the threshold with the highest connected component
count. Ties broken by the lower (more inclusive) threshold.

The DIFF correlation matrix is loaded from Step04's pickle, and the
graph is rebuilt at the auto-selected threshold. This replaces the
previous approach of loading a pre-built graph at a hardcoded threshold.

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
        Gene IDs (ENSG for TCGA, symbols for SC) that passed DE.
        Only A3A and/or A3B.
    sweep_thresholds : list of float
        Thresholds to evaluate.

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
        alias = A3_ID_TO_ALIAS.get(g, g)
        log(f"    {g} ({alias}) -- IN matrix")
    for g in seeds_missing:
        alias = A3_ID_TO_ALIAS.get(g, g)
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
            **{f"deg_{A3_ID_TO_ALIAS.get(g, g)}": d for g, d in seed_degrees.items()},
        })

    sweep_df = pd.DataFrame(rows)

    # Log the sweep table
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

    For TCGA, genes are ENSG IDs. Looks for the DE stats file and checks
    A3A (ENSG00000128383) and A3B (ENSG00000179750).

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
        log(f"  FORCE_KEEP_A3 is ON: A3A/A3B will be used as seeds regardless of DE status")
    sep = "\t" if de_path.endswith(".tsv") else ","
    de_df = pd.read_csv(de_path, sep=sep)

    # Identify gene column
    gene_col = None
    for col in ["gene", "gene_symbol", "Gene", "gene_id"]:
        if col in de_df.columns:
            gene_col = col
            break
    if gene_col is None:
        log(f"  WARNING: Cannot identify gene column. Columns: {list(de_df.columns)}")
        return [a3a_ensg, a3b_ensg]

    # Identify p-value column
    pval_col = None
    for col in ["p_value", "pvalue", "pval", "raw_p", "PValue"]:
        if col in de_df.columns:
            pval_col = col
            break
    if pval_col is None:
        log(f"  WARNING: Cannot identify p-value column. Columns: {list(de_df.columns)}")
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
                log(f"    {label} ({ensg}): FAILED DE (p={pval:.2e}) -- included via FORCE_KEEP")
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

    banner(f"[STEP 15] Community Detection -- {cancer_type}", char="=")

    cancer_dir = ensure_dir(os.path.join(DIR_05_COMMUNITIES, cancer_type))
    sweep_dir = ensure_dir(os.path.join(cancer_dir, "sweep"))

    # =========================================================================
    # NEW: Identify A3 seeds and auto-select DIFF threshold
    # =========================================================================
    banner("[STEP 15.0] Unified DIFF threshold selection")

    # Identify A3 seeds from DE results
    de_cancer_dir = os.path.join(DIR_03_DIFFEXPR, cancer_type)
    a3_seeds = identify_a3_seeds(de_cancer_dir, cancer_type, force_keep=FORCE_KEEP_A3)

    if not a3_seeds:
        log("WARNING: No A3 seed genes available. Falling back to config threshold.")
        active_threshold = DIFF_THRESHOLD
    else:
        log(f"  A3 seeds: {a3_seeds}")

        # Load DIFF correlation matrix
        corr_path = os.path.join(DIR_04_NETWORKS, cancer_type, "corr_matrices",
                                 f"{cancer_type}_corr_DIFF.pkl")
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
                thresh_csv = os.path.join(cancer_dir, f"{cancer_type}_threshold_sweep.csv")
                threshold_sweep_df.to_csv(thresh_csv, index=False)
                log(f"  [SAVE] Threshold sweep -> {thresh_csv}")

    log(f"\n  Active DIFF threshold: {active_threshold}")

    # =========================================================================
    # Build DIFF graph at the selected threshold
    # =========================================================================
    banner("[STEP 15.0b] Build DIFF graph at selected threshold")

    # Check if we already loaded the correlation matrix above
    if 'corr_diff' not in dir():
        # Load it now (fallback path where seeds failed but we still need the graph)
        corr_path = os.path.join(DIR_04_NETWORKS, cancer_type, "corr_matrices",
                                 f"{cancer_type}_corr_DIFF.pkl")
        if not os.path.exists(corr_path):
            # Last resort: try loading pre-built graph at original threshold
            pkl_path = os.path.join(DIR_04_NETWORKS, cancer_type, "graph_objects",
                                    f"{cancer_type}_G_diff_noiso.gpickle")
            if not os.path.exists(pkl_path):
                log(f"[SKIP] Neither DIFF matrix nor pre-built graph found for {cancer_type}")
                continue
            log(f"  Loading pre-built graph (legacy fallback): {pkl_path}")
            with open(pkl_path, "rb") as f:
                G_diff_noiso = pickle.load(f)
        else:
            with open(corr_path, "rb") as f:
                corr_diff = pickle.load(f)

    # Build graph from correlation matrix if we have it
    if 'corr_diff' in dir():
        log(f"  Building graph at |delta-rho| >= {active_threshold}")
        G_diff = build_weighted_graph(corr_diff, active_threshold)
        G_diff_noiso = remove_isolated(G_diff)

        # Save rebuilt graph for downstream steps
        graph_dir = ensure_dir(os.path.join(DIR_04_NETWORKS, cancer_type, "graph_objects"))
        graph_save_path = os.path.join(graph_dir, f"{cancer_type}_G_diff_noiso.gpickle")
        with open(graph_save_path, "wb") as f:
            pickle.dump(G_diff_noiso, f)
        log(f"  [SAVE] Updated DIFF graph -> {graph_save_path}")

        # Clean up to free memory
        del corr_diff

    log(f"  DIFF graph: nodes={G_diff_noiso.number_of_nodes()}, "
        f"edges={G_diff_noiso.number_of_edges()}")

    # ---- Prepare graph for community detection
    banner("[STEP 15.1] Prepare graph")

    G_comm = nx9_remove_degree0(G_diff_noiso.copy())

    if USE_LARGEST_COMPONENT and G_comm.number_of_nodes() > 0:
        G_comm_lcc = nx9_largest_component(G_comm)
        log(f"[STEP 15.1] Before LCC: {G_comm.number_of_nodes()} nodes")
        log(f"[STEP 15.1] After LCC:  {G_comm_lcc.number_of_nodes()} nodes")
        G_comm = G_comm_lcc

        # Report A3 seed status
        for g in a3_seeds:
            alias = A3_ID_TO_ALIAS.get(g, g)
            in_lcc = g in G_comm
            deg = G_diff_noiso.degree(g) if g in G_diff_noiso else 0
            log(f"  {alias}: {'IN LCC' if in_lcc else 'NOT in LCC'} (degree={deg})")

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
        cmap = plt.colormaps["tab20"].resampled(max(len(unique_comms), 1))
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

        # Compute community sizes from merged partition for evenness
        n_total_genes = len(cm_merged)
        largest_comm_size = max(merged_sizes.values()) if merged_sizes else n_total_genes
        n_real_comms = sum(1 for s in merged_sizes.values() if s >= MIN_COMMUNITY_SIZE)

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
            "n_total_genes": n_total_genes,
            "largest_community": largest_comm_size,
            "n_communities_merged": n_real_comms,
        })

    # ---- Save sweep table
    banner("[STEP 15.8] Save sweep summary")

    sweep_df = pd.DataFrame(sweep_rows).sort_values("resolution").reset_index(drop=True)
    sweep_csv = os.path.join(cancer_dir, f"{cancer_type}_resolution_sweep.csv")

    # =========================================================================
    # Compute composite selection metrics
    # =========================================================================
    sweep_df["evenness"] = 1.0 - (sweep_df["largest_community"] / sweep_df["n_total_genes"])

    sweep_df["composite_score"] = (
        sweep_df["modularity_best"] *
        sweep_df["ARI_mean"] *
        sweep_df["evenness"]
    )

    sweep_df["delta_modularity"] = sweep_df["modularity_best"].diff().fillna(0)

    max_mod = sweep_df["modularity_best"].max()
    sweep_df["delta_normalized"] = sweep_df["delta_modularity"] / (max_mod + 1e-10)

    sweep_df.to_csv(sweep_csv, index=False)
    log(f"[SAVE] Sweep table -> {sweep_csv}")

    # ---- Log the full sweep table
    log("\n[STEP 15.8] Resolution sweep summary:")
    log(f"{'Res':>5s} | {'Mod':>6s} | {'ARI':>5s} | {'Even':>5s} | {'Comp':>6s} | {'dMod':>6s} | {'Comms':>5s} | {'Largest':>7s}")
    log("-" * 65)
    for _, row in sweep_df.iterrows():
        log(f"{row['resolution']:5.2f} | {row['modularity_best']:6.4f} | "
            f"{row['ARI_mean']:5.3f} | {row['evenness']:5.3f} | "
            f"{row['composite_score']:6.4f} | {row['delta_modularity']:+6.4f} | "
            f"{row['ncomms_mean']:5.1f} | {row['largest_community']:7.0f}")

    # =========================================================================
    # Diagnostic sweep plots (6 panels)
    # =========================================================================
    fig, axes = plt.subplots(2, 3, figsize=(22, 12))

    # Panel 1: Modularity vs resolution
    ax = axes[0, 0]
    ax.errorbar(sweep_df["resolution"], sweep_df["modularity_best"],
                yerr=sweep_df["modularity_std"], marker="o", capsize=3, color="steelblue")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Modularity", fontsize=12)
    ax.set_title(f"{cancer_type} | Modularity", fontsize=13)

    # Panel 2: ARI stability vs resolution
    ax = axes[0, 1]
    ax.errorbar(sweep_df["resolution"], sweep_df["ARI_mean"],
                yerr=sweep_df["ARI_std"], marker="o", capsize=3, color="firebrick")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("ARI (stability)", fontsize=12)
    ax.set_title(f"{cancer_type} | Stability (ARI)", fontsize=13)

    # Panel 3: Evenness vs resolution
    ax = axes[0, 2]
    ax.plot(sweep_df["resolution"], sweep_df["evenness"], "o-", color="forestgreen")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Evenness", fontsize=12)
    ax.set_title(f"{cancer_type} | Evenness (1 - largest/total)", fontsize=13)
    ax.set_ylim(-0.05, 1.05)

    # Panel 4: Composite score
    ax = axes[1, 0]
    ax.plot(sweep_df["resolution"], sweep_df["composite_score"], "o-",
            color="darkorange", linewidth=2, markersize=8)
    best_idx = sweep_df["composite_score"].idxmax()
    ax.scatter([sweep_df.loc[best_idx, "resolution"]],
              [sweep_df.loc[best_idx, "composite_score"]],
              s=200, c="red", zorder=5, edgecolors="black", linewidths=2)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Composite (Mod x ARI x Even)", fontsize=12)
    ax.set_title(f"{cancer_type} | Composite Score", fontsize=13, fontweight="bold")

    # Panel 5: Delta-modularity
    ax = axes[1, 1]
    ax.bar(sweep_df["resolution"], sweep_df["delta_modularity"],
           width=0.08, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("Delta Modularity", fontsize=12)
    ax.set_title(f"{cancer_type} | Delta-Modularity (Elbow)", fontsize=13)

    # Panel 6: Community count
    ax = axes[1, 2]
    ax.errorbar(sweep_df["resolution"], sweep_df["ncomms_mean"],
                yerr=sweep_df["ncomms_std"], marker="o", capsize=3, color="purple")
    ax.set_xlabel("Resolution", fontsize=12)
    ax.set_ylabel("N Communities", fontsize=12)
    ax.set_title(f"{cancer_type} | Community Count", fontsize=13)

    plt.suptitle(f"{cancer_type} -- Resolution Selection Diagnostics "
                 f"(threshold={active_threshold})", fontsize=16, fontweight="bold")
    plt.tight_layout()
    sweep_plot = os.path.join(cancer_dir, f"{cancer_type}_sweep_diagnostics.png")
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
    log(f"  Communities: {best_row['ncomms_mean']:.1f}")
    log(f"  Largest community: {best_row['largest_community']:.0f} genes")

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

    # Save the community graph
    comm_pkl = os.path.join(cancer_dir, f"{cancer_type}_G_comm.gpickle")
    with open(comm_pkl, "wb") as f:
        pickle.dump(G_comm, f)
    log(f"[SAVE] Community graph -> {comm_pkl}")

    # Save selected parameters (new output)
    params_path = os.path.join(cancer_dir, f"{cancer_type}_selected_parameters.txt")
    with open(params_path, "w") as f:
        f.write(f"DIFF_THRESHOLD={active_threshold}\n")
        f.write(f"THRESHOLD_METHOD=unified_a3_connectivity_component_peak\n")
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
