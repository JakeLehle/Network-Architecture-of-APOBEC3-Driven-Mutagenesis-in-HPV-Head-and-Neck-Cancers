#!/usr/bin/env python3
"""
Tier3A_Leave_One_Patient_Out.py
================================

Leave-one-patient-out (LOPO) network sensitivity analysis.
Mirrors the V4 SC network pipeline (Figure 4) exactly.

For each excluded patient:
  1. Remove that patient's cells from BOTH SBS2-HIGH and NORMAL groups
  2. Re-run the full V4 pipeline:
     a. DE via scanpy rank_genes_groups (Wilcoxon, FDR < 0.05)
     b. Spearman correlation matrices (HIGH / NORMAL / DIFF)
     c. Max fragmentation rate threshold auto-selection
     d. Full-network Leiden with resolution sweep + stability
     e. Component-aware merge (satellites preserved)
  3. Compare to full analysis: A3 placement, Harris interactors,
     concordant chain gene recovery, community structure overlap

The test: if the network is a generalizable program across HNSCC
patients (not driven by one individual), removing any single patient
should preserve the core architecture.

Usage:
  conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC029"
  conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC013"
  conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC001"

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import argparse
import itertools
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.stats import rankdata
import igraph as ig
import leidenalg
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from patient_config import (
    DIR_03_SENSITIVITY, ADATA_PATH, WEIGHTS_PATH, PATIENT_COL,
    CELLTYPE_COL, RANDOM_SEED, DE_FDR_THRESHOLD, FORCE_KEEP_A3,
    DIFF_THRESHOLD, DIFF_THRESHOLD_AUTO, SWEEP_THRESHOLDS,
    COMMUNITY_RESOLUTIONS, RUNS_PER_RESOLUTION, COMMUNITY_BASE_SEED,
    MIN_COMMUNITY_SIZE, A3_GENES_SYMBOLS,
    FULL_ANALYSIS_PARTITION, COMMUNITIES_DIR, NETWORK_DIR,
    banner, log, ensure_dir, load_adata, load_three_groups,
    load_harris_interactors, load_full_analysis_params,
)

# Concordant chain genes from Figure 4 SBS2_VS_NORMAL analysis
ACTIVATING_CHAIN = [
    "RALY", "HNRNPA2B1", "CCL20", "KRT24", "LCN2",
    "LINC00278", "RRAD", "SMOX", "UTY",
]
INHIBITING_CHAIN_ANCHORS = ["SNHG3", "THYN1", "ZNG1A"]


# =============================================================================
# STEP 1: DIFFERENTIAL EXPRESSION (V4: scanpy rank_genes_groups)
# =============================================================================

def run_de(adata, high_cells, normal_cells, run_dir, short):
    """Differential expression using scanpy Wilcoxon with BH-FDR."""
    banner(f"  DE via scanpy (LOPO {short})")

    # Subset adata to basal cells in our two groups
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell']
    all_group_cells = list(high_cells | normal_cells)
    in_adata = [c for c in all_group_cells if c in basal.obs_names]
    adata_sub = basal[in_adata].copy()

    # Assign group labels
    adata_sub.obs['_lopo_group'] = 'NORMAL'
    high_mask = adata_sub.obs_names.isin(high_cells)
    adata_sub.obs.loc[high_mask, '_lopo_group'] = 'SBS2_HIGH'

    n_h = high_mask.sum()
    n_n = (~high_mask).sum()
    log(f"  Cells: {n_h} SBS2-HIGH, {n_n} NORMAL")

    # Run scanpy DE (no additional normalization, data already processed)
    sc.tl.rank_genes_groups(
        adata_sub, groupby='_lopo_group', groups=['SBS2_HIGH'],
        reference='NORMAL', method='wilcoxon', use_raw=False,
        corr_method='benjamini-hochberg',
    )

    # Extract results
    result = adata_sub.uns['rank_genes_groups']
    de_df = pd.DataFrame({
        'gene': result['names']['SBS2_HIGH'],
        'logfoldchange': result['logfoldchanges']['SBS2_HIGH'],
        'pval': result['pvals']['SBS2_HIGH'],
        'pval_adj': result['pvals_adj']['SBS2_HIGH'],
        'score': result['scores']['SBS2_HIGH'],
    })

    # Select DE genes at FDR < 0.05
    de_sig = de_df[de_df['pval_adj'] < DE_FDR_THRESHOLD].copy()
    log(f"  DE genes (FDR < {DE_FDR_THRESHOLD}): {len(de_sig):,}")

    # Check A3 genes
    for a3 in ['APOBEC3A', 'APOBEC3B']:
        row = de_df[de_df['gene'] == a3]
        if len(row) > 0:
            r = row.iloc[0]
            status = "PASS" if r['pval_adj'] < DE_FDR_THRESHOLD else "FAIL"
            log(f"    {a3}: FDR={r['pval_adj']:.2e}, log2FC={r['logfoldchange']:.2f} [{status}]")

    # Save DE results
    de_df.to_csv(os.path.join(run_dir, f"LOPO_{short}_DE_all.tsv"),
                 sep='\t', index=False)
    de_sig.to_csv(os.path.join(run_dir, f"LOPO_{short}_DE_significant.tsv"),
                  sep='\t', index=False)

    de_genes = de_sig['gene'].tolist()
    return adata_sub, de_genes, de_df


# =============================================================================
# STEP 2: SPEARMAN CORRELATION MATRICES
# =============================================================================

def compute_spearman_matrix(X):
    """Compute Spearman correlation matrix via rank transform + Pearson."""
    X_ranked = np.apply_along_axis(rankdata, 0, X)
    rho = np.corrcoef(X_ranked.T)
    np.fill_diagonal(rho, 0)
    return np.nan_to_num(rho)


def run_correlations(adata_sub, high_cells, de_genes, run_dir, short):
    """Compute HIGH, NORMAL, and DIFF Spearman correlation matrices."""
    banner(f"  SPEARMAN CORRELATIONS (LOPO {short})")

    # Get expression for DE genes only
    gene_mask = np.isin(adata_sub.var_names, de_genes)
    gene_idx = np.where(gene_mask)[0]
    gene_names = adata_sub.var_names[gene_idx]

    X = adata_sub.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    X_de = X[:, gene_idx]

    mask_h = adata_sub.obs_names.isin(high_cells)
    X_high = X_de[mask_h]
    X_norm = X_de[~mask_h]

    n_genes = len(gene_names)
    log(f"  Computing {n_genes}x{n_genes} Spearman matrices...")
    log(f"  HIGH: {X_high.shape[0]} cells, NORMAL: {X_norm.shape[0]} cells")

    rho_high = compute_spearman_matrix(X_high)
    rho_norm = compute_spearman_matrix(X_norm)
    diff = rho_high - rho_norm

    log(f"  DIFF range: [{diff.min():.3f}, {diff.max():.3f}]")
    log(f"  DIFF mean abs: {np.abs(diff).mean():.4f}")

    return gene_names, diff


# =============================================================================
# STEP 3a: MAX FRAGMENTATION RATE THRESHOLD SELECTION (V4)
# =============================================================================

def select_threshold_max_frag(gene_names, diff, run_dir, short):
    """Auto-select DIFF threshold via max fragmentation rate.

    Sweeps thresholds, computes connected component count at each.
    Among A3-valid thresholds (A3A and A3B have degree >= 1),
    selects the upper threshold of the interval with the largest
    positive delta-component count.
    """
    banner(f"  THRESHOLD SELECTION (LOPO {short})")

    n_genes = len(gene_names)
    a3a_idx = np.where(gene_names == 'APOBEC3A')[0]
    a3b_idx = np.where(gene_names == 'APOBEC3B')[0]

    sweep_data = []
    for thr in SWEEP_THRESHOLDS:
        # Build adjacency at this threshold
        adj = np.abs(diff) >= thr
        np.fill_diagonal(adj, False)

        # Node degrees
        degrees = adj.sum(axis=0)
        n_nodes = (degrees > 0).sum()
        n_edges = adj.sum() // 2

        # A3 degrees
        a3a_deg = int(degrees[a3a_idx[0]]) if len(a3a_idx) > 0 else 0
        a3b_deg = int(degrees[a3b_idx[0]]) if len(a3b_idx) > 0 else 0
        a3_valid = (a3a_deg >= 1) and (a3b_deg >= 1)

        # Connected components via NetworkX on non-isolate subgraph
        G_tmp = nx.Graph()
        active_genes = gene_names[degrees > 0]
        G_tmp.add_nodes_from(active_genes)
        for i in range(n_genes):
            if degrees[i] == 0:
                continue
            for j in range(i + 1, n_genes):
                if adj[i, j]:
                    G_tmp.add_edge(gene_names[i], gene_names[j])
        n_comp = nx.number_connected_components(G_tmp) if G_tmp.number_of_nodes() > 0 else 0

        sweep_data.append({
            'threshold': thr, 'n_nodes': n_nodes, 'n_edges': n_edges,
            'n_components': n_comp, 'A3A_degree': a3a_deg,
            'A3B_degree': a3b_deg, 'A3_valid': a3_valid,
        })
        log(f"    thr={thr:.2f}: {n_nodes} nodes, {n_edges} edges, "
            f"{n_comp} comp, A3A={a3a_deg}, A3B={a3b_deg} "
            f"{'[A3-valid]' if a3_valid else ''}")

    sweep_df = pd.DataFrame(sweep_data)
    sweep_df.to_csv(os.path.join(run_dir, f"LOPO_{short}_threshold_sweep.tsv"),
                    sep='\t', index=False)

    # Find max fragmentation rate among A3-valid intervals
    selected = DIFF_THRESHOLD  # fallback
    best_delta = -1

    for k in range(len(sweep_data) - 1):
        lower = sweep_data[k]
        upper = sweep_data[k + 1]
        if not upper['A3_valid']:
            continue
        delta_comp = upper['n_components'] - lower['n_components']
        if delta_comp > best_delta:
            best_delta = delta_comp
            selected = upper['threshold']
        elif delta_comp == best_delta and delta_comp > 0:
            # Tiebreak: prefer lower (more inclusive) threshold
            selected = min(selected, upper['threshold'])

    if best_delta <= 0:
        log(f"  WARNING: No positive delta-comp at any A3-valid threshold.")
        log(f"  Falling back to config DIFF_THRESHOLD={DIFF_THRESHOLD}")
        selected = DIFF_THRESHOLD

    log(f"  Selected threshold: {selected:.2f} (delta-comp={best_delta})")
    return selected, sweep_df


# =============================================================================
# STEP 3b: BUILD DIFF NETWORK
# =============================================================================

def build_diff_network(gene_names, diff, threshold):
    """Build DIFF network at the selected threshold. Keep full network."""
    n_genes = len(gene_names)
    G = nx.Graph()
    G.add_nodes_from(gene_names)

    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            if abs(diff[i, j]) >= threshold:
                G.add_edge(gene_names[i], gene_names[j],
                           weight=diff[i, j],
                           abs_weight=abs(diff[i, j]))

    # Remove isolates (nodes with zero edges)
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)

    log(f"  DIFF network: {G.number_of_nodes()} nodes, "
        f"{G.number_of_edges()} edges")
    log(f"  Connected components: {nx.number_connected_components(G)}")

    return G


# =============================================================================
# STEP 3c: FULL-NETWORK LEIDEN WITH RESOLUTION SWEEP (V4)
# =============================================================================

def detect_leiden(G, seed, resolution):
    """Run Leiden on full graph, return node->community dict."""
    node_list = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(node_list)}

    ig_edges = [(node_idx[u], node_idx[v]) for u, v in G.edges()]
    ig_weights = [G[u][v]['abs_weight'] for u, v in G.edges()]

    ig_graph = ig.Graph(n=len(node_list), edges=ig_edges)
    ig_graph.es['weight'] = ig_weights

    partition = leidenalg.find_partition(
        ig_graph, leidenalg.RBConfigurationVertexPartition,
        weights='weight', resolution_parameter=resolution, seed=seed,
    )

    return {node_list[i]: partition.membership[i]
            for i in range(len(node_list))}


def run_leiden_sweep(G, run_dir, short):
    """Resolution sweep with stability assessment. Returns best partition."""
    banner(f"  LEIDEN RESOLUTION SWEEP (LOPO {short})")

    nodes_list = sorted(G.nodes())

    def partition_to_labels(cm):
        return [cm.get(n, -1) for n in nodes_list]

    sweep_results = []
    all_best = {}

    for res in COMMUNITY_RESOLUTIONS:
        partitions = []
        modularities = []
        best_mod = -999
        best_cm = None

        for r in range(RUNS_PER_RESOLUTION):
            seed = COMMUNITY_BASE_SEED + r
            cm = detect_leiden(G, seed=seed, resolution=res)
            partitions.append(cm)

            # Compute modularity
            comm_sets = {}
            for n, c in cm.items():
                comm_sets.setdefault(c, set()).add(n)
            try:
                mod = nx.algorithms.community.quality.modularity(
                    G, list(comm_sets.values()), weight='abs_weight')
            except Exception:
                mod = 0.0
            modularities.append(mod)

            if mod > best_mod:
                best_mod = mod
                best_cm = cm.copy()

        # Stability: pairwise ARI/NMI
        labels_list = [partition_to_labels(p) for p in partitions]
        ari_pairs = []
        nmi_pairs = []
        for i, j in itertools.combinations(range(len(labels_list)), 2):
            ari_pairs.append(adjusted_rand_score(labels_list[i], labels_list[j]))
            nmi_pairs.append(normalized_mutual_info_score(
                labels_list[i], labels_list[j]))

        mean_ari = np.mean(ari_pairs) if ari_pairs else 0
        mean_nmi = np.mean(nmi_pairs) if nmi_pairs else 0
        mean_mod = np.mean(modularities)

        # Evenness: 1 - (largest_comm / total_nodes)
        comm_sizes = {}
        for n, c in best_cm.items():
            comm_sizes[c] = comm_sizes.get(c, 0) + 1
        largest = max(comm_sizes.values()) if comm_sizes else 1
        evenness = 1.0 - (largest / len(nodes_list)) if nodes_list else 0

        # Composite score (V4: modularity x ARI x evenness)
        composite = mean_mod * mean_ari * evenness

        n_comms = len(set(best_cm.values()))
        log(f"  res={res:.1f}: mod={mean_mod:.3f}, ARI={mean_ari:.3f}, "
            f"NMI={mean_nmi:.3f}, even={evenness:.3f}, "
            f"composite={composite:.4f}, comms={n_comms}")

        sweep_results.append({
            'resolution': res, 'modularity': mean_mod,
            'ARI': mean_ari, 'NMI': mean_nmi,
            'evenness': evenness, 'composite': composite,
            'n_communities': n_comms,
        })
        all_best[res] = best_cm

    sweep_df = pd.DataFrame(sweep_results)
    sweep_df.to_csv(os.path.join(run_dir, f"LOPO_{short}_resolution_sweep.tsv"),
                    sep='\t', index=False)

    # Select best resolution by composite score
    best_idx = sweep_df['composite'].idxmax()
    best_res = sweep_df.loc[best_idx, 'resolution']
    log(f"  Selected resolution: {best_res:.1f} "
        f"(composite={sweep_df.loc[best_idx, 'composite']:.4f})")

    return all_best[best_res], best_res, sweep_df


# =============================================================================
# STEP 3d: COMPONENT-AWARE MERGE (V4)
# =============================================================================

def component_aware_merge(G, partition):
    """Merge small Leiden communities within large components,
    but preserve satellites (entire small components) as-is."""
    banner("  COMPONENT-AWARE MERGE")

    # Get connected components
    components = list(nx.connected_components(G))
    comp_sizes = sorted([len(c) for c in components], reverse=True)
    log(f"  Components: {len(components)} "
        f"(sizes: {comp_sizes[:10]}{'...' if len(comp_sizes) > 10 else ''})")

    # Identify which component each node belongs to
    node_to_comp = {}
    for i, comp in enumerate(components):
        for n in comp:
            node_to_comp[n] = i

    # Build community -> set of nodes
    comm_nodes = {}
    for n, c in partition.items():
        comm_nodes.setdefault(c, set()).add(n)

    # Identify large components (size > MIN_COMMUNITY_SIZE)
    large_comp_ids = {i for i, comp in enumerate(components)
                      if len(comp) > MIN_COMMUNITY_SIZE}

    # For each small community, check if it's entirely within a large component
    merged = partition.copy()
    merge_count = 0

    for comm_id, nodes in list(comm_nodes.items()):
        if len(nodes) >= MIN_COMMUNITY_SIZE:
            continue

        # Check if this community is an entire component (satellite)
        comp_ids = {node_to_comp[n] for n in nodes}
        if len(comp_ids) == 1:
            comp_id = comp_ids.pop()
            comp_nodes = components[comp_id] if comp_id < len(components) else set()
            if nodes == comp_nodes or (isinstance(comp_nodes, set) and nodes == comp_nodes):
                # This community IS the entire component: preserve as satellite
                continue

        # Small community within a large component: merge into nearest neighbor
        if not any(node_to_comp.get(n) in large_comp_ids for n in nodes):
            continue

        # Find the neighboring community with most edges to this one
        neighbor_comms = {}
        for n in nodes:
            for nb in G.neighbors(n):
                nb_comm = merged.get(nb)
                if nb_comm is not None and nb_comm != comm_id:
                    neighbor_comms[nb_comm] = neighbor_comms.get(nb_comm, 0) + 1

        if neighbor_comms:
            target = max(neighbor_comms, key=neighbor_comms.get)
            for n in nodes:
                merged[n] = target
            merge_count += 1

    # Renumber communities consecutively
    unique_comms = sorted(set(merged.values()))
    remap = {old: new for new, old in enumerate(unique_comms)}
    merged = {n: remap[c] for n, c in merged.items()}

    # Summary
    final_comm_nodes = {}
    for n, c in merged.items():
        final_comm_nodes.setdefault(c, []).append(n)

    n_main = sum(1 for v in final_comm_nodes.values()
                 if len(v) >= MIN_COMMUNITY_SIZE)
    n_sat = sum(1 for v in final_comm_nodes.values()
                if len(v) < MIN_COMMUNITY_SIZE)
    log(f"  Merged {merge_count} small communities")
    log(f"  Final: {n_main} main + {n_sat} satellite communities")

    return merged, final_comm_nodes


# =============================================================================
# COMPARISON TO FULL ANALYSIS
# =============================================================================

def compare_to_full(G, merged, final_comms, de_genes, run_dir, short):
    """Compare LOPO network to the full SBS2_VS_NORMAL analysis."""
    banner(f"  COMPARISON TO FULL ANALYSIS (LOPO {short})")

    harris_all, harris_a3b = load_harris_interactors()
    full_params = load_full_analysis_params()

    # -- A3 gene status --
    a3_status = {}
    for a3 in ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3G']:
        in_net = a3 in G.nodes()
        deg = G.degree(a3) if in_net else 0
        comm = merged.get(a3, -1)
        comm_size = len(final_comms.get(comm, []))
        a3_status[a3] = {
            'in_network': in_net, 'degree': deg,
            'community': comm, 'community_size': comm_size,
        }
        if in_net:
            log(f"  {a3}: degree={deg}, community={comm} (size={comm_size})")
        else:
            log(f"  {a3}: NOT in network")

    # -- Harris interactors --
    harris_in = {g for g in harris_all if g in G.nodes()}
    log(f"  Harris interactors in network: {len(harris_in)} / {len(harris_all)}")

    # -- Concordant chain gene recovery --
    act_in = [g for g in ACTIVATING_CHAIN if g in G.nodes()]
    inh_in = [g for g in INHIBITING_CHAIN_ANCHORS if g in G.nodes()]
    log(f"  Activating chain genes recovered: {len(act_in)}/{len(ACTIVATING_CHAIN)} "
        f"({', '.join(act_in) if act_in else 'none'})")
    log(f"  Inhibiting chain anchors recovered: {len(inh_in)}/{len(INHIBITING_CHAIN_ANCHORS)} "
        f"({', '.join(inh_in) if inh_in else 'none'})")

    # -- Gene overlap with full analysis partition --
    full_genes = set()
    full_partition_map = {}
    if os.path.exists(FULL_ANALYSIS_PARTITION):
        fp = pd.read_csv(FULL_ANALYSIS_PARTITION)
        gene_col = 'gene' if 'gene' in fp.columns else fp.columns[0]
        comm_col = 'community' if 'community' in fp.columns else fp.columns[1]
        full_genes = set(fp[gene_col].values)
        full_partition_map = dict(zip(fp[gene_col], fp[comm_col]))
        overlap = full_genes & set(G.nodes())
        jaccard = len(overlap) / len(full_genes | set(G.nodes())) if (full_genes | set(G.nodes())) else 0
        log(f"  Gene overlap with full analysis: {len(overlap)} / "
            f"{len(full_genes)} full, {len(set(G.nodes()))} LOPO "
            f"(Jaccard={jaccard:.3f})")

        # ARI between full and LOPO partitions on shared genes
        shared = sorted(overlap)
        if len(shared) > 50:
            full_labels = [full_partition_map.get(g, -1) for g in shared]
            lopo_labels = [merged.get(g, -1) for g in shared]
            ari = adjusted_rand_score(full_labels, lopo_labels)
            log(f"  Community structure ARI (on {len(shared)} shared genes): {ari:.3f}")
        else:
            ari = None
            log(f"  Too few shared genes ({len(shared)}) for ARI comparison")
    else:
        overlap = set()
        jaccard = 0
        ari = None
        log(f"  Full analysis partition not found, skipping overlap comparison")

    # -- A3 wall check (are A3 DIFF edges still all negative?) --
    a3_wall = {}
    for a3 in ['APOBEC3A', 'APOBEC3B']:
        if a3 not in G.nodes():
            a3_wall[a3] = {'n_pos': 0, 'n_neg': 0, 'n_total': 0, 'pct_neg': None}
            continue
        neighbors = list(G.neighbors(a3))
        pos = sum(1 for nb in neighbors if G[a3][nb]['weight'] > 0)
        neg = sum(1 for nb in neighbors if G[a3][nb]['weight'] <= 0)
        pct_neg = neg / len(neighbors) * 100 if neighbors else None
        a3_wall[a3] = {'n_pos': pos, 'n_neg': neg, 'n_total': len(neighbors),
                        'pct_neg': pct_neg}
        log(f"  A3 wall check {a3}: {neg}/{len(neighbors)} negative DIFF edges "
            f"({pct_neg:.0f}%)" if pct_neg is not None else
            f"  A3 wall check {a3}: no edges")

    return {
        'a3_status': a3_status,
        'harris_in_network': len(harris_in),
        'harris_total': len(harris_all),
        'activating_chain_recovered': len(act_in),
        'activating_chain_total': len(ACTIVATING_CHAIN),
        'activating_chain_genes': act_in,
        'inhibiting_anchors_recovered': len(inh_in),
        'inhibiting_anchors_total': len(INHIBITING_CHAIN_ANCHORS),
        'gene_overlap': len(overlap) if overlap else 0,
        'jaccard': jaccard,
        'community_ari': ari,
        'a3_wall': a3_wall,
    }


# =============================================================================
# MAIN LOPO FUNCTION
# =============================================================================

def run_lopo(exclude_patient):
    """Run full V4 network pipeline excluding one patient."""

    short = exclude_patient.replace("Patient ", "")
    banner(f"TIER 3A: LOPO -- EXCLUDING {exclude_patient}")
    run_dir = ensure_dir(os.path.join(DIR_03_SENSITIVITY, f"LOPO_{short}"))

    # -----------------------------------------------------------------
    # Load data and remove patient from BOTH groups
    # -----------------------------------------------------------------
    adata = load_adata()
    sbs2_high, cnv_high, normal = load_three_groups()

    patient_map = adata.obs[PATIENT_COL].to_dict()

    # Remove excluded patient from SBS2-HIGH
    high_orig = len(sbs2_high)
    high_filtered = {c for c in sbs2_high if patient_map.get(c) != exclude_patient}
    high_removed = high_orig - len(high_filtered)

    # Remove excluded patient from NORMAL
    norm_orig = len(normal)
    norm_filtered = {c for c in normal if patient_map.get(c) != exclude_patient}
    norm_removed = norm_orig - len(norm_filtered)

    log(f"  SBS2-HIGH: {high_orig} -> {len(high_filtered)} "
        f"(removed {high_removed} from {exclude_patient})")
    log(f"  NORMAL:    {norm_orig} -> {len(norm_filtered)} "
        f"(removed {norm_removed} from {exclude_patient})")

    if len(high_filtered) < 50:
        log("  ERROR: Too few SBS2-HIGH cells remaining. Skipping.")
        return None

    if len(norm_filtered) < 50:
        log("  ERROR: Too few NORMAL cells remaining. Skipping.")
        return None

    # -----------------------------------------------------------------
    # Step 1: Differential expression
    # -----------------------------------------------------------------
    adata_sub, de_genes, de_df = run_de(
        adata, high_filtered, norm_filtered, run_dir, short)

    if len(de_genes) < 100:
        log(f"  WARNING: Only {len(de_genes)} DE genes. Network may be very sparse.")

    # -----------------------------------------------------------------
    # Step 2: Spearman correlations
    # -----------------------------------------------------------------
    gene_names, diff = run_correlations(
        adata_sub, high_filtered, de_genes, run_dir, short)

    # -----------------------------------------------------------------
    # Step 3a: Threshold selection
    # -----------------------------------------------------------------
    if DIFF_THRESHOLD_AUTO:
        selected_thr, sweep_df = select_threshold_max_frag(
            gene_names, diff, run_dir, short)
    else:
        selected_thr = DIFF_THRESHOLD
        sweep_df = None
        log(f"  Using fixed threshold: {selected_thr}")

    # -----------------------------------------------------------------
    # Step 3b: Build DIFF network
    # -----------------------------------------------------------------
    banner(f"  DIFF NETWORK (LOPO {short})")
    G = build_diff_network(gene_names, diff, selected_thr)

    if G.number_of_nodes() < 10:
        log("  Network too small for community detection. Skipping.")
        result = {
            'excluded': exclude_patient,
            'n_high': len(high_filtered), 'n_normal': len(norm_filtered),
            'high_removed': high_removed, 'normal_removed': norm_removed,
            'n_de_genes': len(de_genes),
            'diff_threshold': selected_thr,
            'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(),
            'n_communities': 0,
            'status': 'NETWORK_TOO_SMALL',
        }
        pd.DataFrame([result]).to_csv(
            os.path.join(run_dir, f"LOPO_{short}_summary.tsv"),
            sep='\t', index=False)
        return result

    # -----------------------------------------------------------------
    # Step 3c: Leiden with resolution sweep
    # -----------------------------------------------------------------
    best_partition, best_res, res_sweep = run_leiden_sweep(G, run_dir, short)

    # -----------------------------------------------------------------
    # Step 3d: Component-aware merge
    # -----------------------------------------------------------------
    merged, final_comms = component_aware_merge(G, best_partition)

    # -----------------------------------------------------------------
    # Save community assignments
    # -----------------------------------------------------------------
    comm_rows = []
    for comm_id, genes_list in sorted(final_comms.items()):
        for g in sorted(genes_list):
            comm_rows.append({
                'gene': g, 'community': comm_id,
                'degree': G.degree(g),
            })
    comm_df = pd.DataFrame(comm_rows)
    comm_df.to_csv(os.path.join(run_dir, f"LOPO_{short}_communities.tsv"),
                   sep='\t', index=False)

    # -----------------------------------------------------------------
    # Comparison to full analysis
    # -----------------------------------------------------------------
    comparison = compare_to_full(G, merged, final_comms, de_genes,
                                  run_dir, short)

    # -----------------------------------------------------------------
    # Build summary result
    # -----------------------------------------------------------------
    a3a = comparison['a3_status'].get('APOBEC3A', {})
    a3b = comparison['a3_status'].get('APOBEC3B', {})

    n_main = sum(1 for v in final_comms.values() if len(v) >= MIN_COMMUNITY_SIZE)
    n_sat = sum(1 for v in final_comms.values() if len(v) < MIN_COMMUNITY_SIZE)

    result = {
        'excluded': exclude_patient,
        'n_high': len(high_filtered),
        'n_normal': len(norm_filtered),
        'high_removed': high_removed,
        'normal_removed': norm_removed,
        'n_de_genes': len(de_genes),
        'diff_threshold': selected_thr,
        'leiden_resolution': best_res,
        'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(),
        'n_components': nx.number_connected_components(G),
        'n_communities_main': n_main,
        'n_communities_satellite': n_sat,
        'n_communities_total': n_main + n_sat,
        'A3A_in_network': a3a.get('in_network', False),
        'A3A_degree': a3a.get('degree', 0),
        'A3A_community': a3a.get('community', -1),
        'A3A_community_size': a3a.get('community_size', 0),
        'A3B_in_network': a3b.get('in_network', False),
        'A3B_degree': a3b.get('degree', 0),
        'A3B_community': a3b.get('community', -1),
        'harris_in_network': comparison['harris_in_network'],
        'harris_total': comparison['harris_total'],
        'activating_chain_recovered': comparison['activating_chain_recovered'],
        'activating_chain_total': comparison['activating_chain_total'],
        'activating_chain_genes': str(comparison['activating_chain_genes']),
        'inhibiting_anchors_recovered': comparison['inhibiting_anchors_recovered'],
        'gene_overlap_with_full': comparison['gene_overlap'],
        'jaccard_with_full': comparison['jaccard'],
        'community_ari_with_full': comparison['community_ari'],
        'a3a_wall_pct_neg': comparison['a3_wall'].get('APOBEC3A', {}).get('pct_neg'),
        'a3b_wall_pct_neg': comparison['a3_wall'].get('APOBEC3B', {}).get('pct_neg'),
        'status': 'COMPLETE',
    }

    # Save summary
    pd.DataFrame([result]).to_csv(
        os.path.join(run_dir, f"LOPO_{short}_summary.tsv"),
        sep='\t', index=False)

    # Save selected parameters (for reference)
    with open(os.path.join(run_dir, f"LOPO_{short}_parameters.txt"), 'w') as f:
        f.write(f"DIFF_THRESHOLD={selected_thr}\n")
        f.write(f"LEIDEN_RESOLUTION={best_res}\n")
        f.write(f"N_COMMUNITIES={n_main + n_sat}\n")
        f.write(f"N_GENES={G.number_of_nodes()}\n")
        f.write(f"N_EDGES={G.number_of_edges()}\n")
        f.write(f"N_HIGH={len(high_filtered)}\n")
        f.write(f"N_NORMAL={len(norm_filtered)}\n")

    log(f"\n  LOPO {short} complete. Saved to {run_dir}")
    return result


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Leave-one-patient-out network sensitivity (V4 pipeline)")
    parser.add_argument('--exclude', type=str, required=True,
                        help='Patient to exclude (e.g. "Patient SC029")')
    args = parser.parse_args()

    result = run_lopo(args.exclude)

    if result:
        banner("RESULT SUMMARY")
        for k, v in result.items():
            log(f"  {k}: {v}")

        # Quick pass/fail assessment
        banner("SENSITIVITY ASSESSMENT")
        if result.get('status') != 'COMPLETE':
            log(f"  FAIL: Pipeline did not complete ({result.get('status')})")
        else:
            checks = []
            if result['A3A_in_network']:
                checks.append("A3A in network: PASS")
            else:
                checks.append("A3A in network: FAIL")

            act_pct = (result['activating_chain_recovered'] /
                       result['activating_chain_total'] * 100)
            if act_pct >= 50:
                checks.append(f"Activating chain: PASS ({act_pct:.0f}% recovered)")
            else:
                checks.append(f"Activating chain: PARTIAL ({act_pct:.0f}% recovered)")

            if result['community_ari_with_full'] is not None:
                if result['community_ari_with_full'] > 0.3:
                    checks.append(f"Community ARI: PASS ({result['community_ari_with_full']:.3f})")
                else:
                    checks.append(f"Community ARI: LOW ({result['community_ari_with_full']:.3f})")

            wall_a3a = result.get('a3a_wall_pct_neg')
            if wall_a3a is not None and wall_a3a >= 90:
                checks.append(f"A3A wall: INTACT ({wall_a3a:.0f}% negative)")
            elif wall_a3a is not None:
                checks.append(f"A3A wall: DEGRADED ({wall_a3a:.0f}% negative)")

            for c in checks:
                log(f"  {c}")


if __name__ == "__main__":
    main()
