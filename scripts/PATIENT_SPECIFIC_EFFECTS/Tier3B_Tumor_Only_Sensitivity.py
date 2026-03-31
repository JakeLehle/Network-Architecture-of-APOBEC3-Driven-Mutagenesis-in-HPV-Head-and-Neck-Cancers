#!/usr/bin/env python3
"""
Tier3B_Tumor_Only_Sensitivity.py (v2 — runs both approaches)
==============================================================

Approach A: Remove 25 normal-tissue HIGH cells, KEEP original 546 LOWs.
            Unbalanced (521 HIGH vs 546 LOW) but preserves control set.
Approach B: Remove 25 normal HIGH, re-match 521 LOWs from tumor-only basal.
            Balanced but different control population.

Comparing both tells us whether the A3A dropout in v1 was due to the
control re-selection or the loss of 25 normal HIGH cells.

Usage: conda run -n NETWORK python Tier3B_Tumor_Only_Sensitivity.py
"""

import os, sys, numpy as np, pandas as pd, scanpy as sc, scipy.sparse
from scipy.stats import ranksums, rankdata
import igraph as ig, leidenalg, networkx as nx
import matplotlib; matplotlib.use('Agg')

from patient_config import *


def run_network_pipeline(X_high, X_low, gene_names, label):
    """
    Run DE -> Spearman -> DIFF -> Leiden -> Centrality.
    Returns result dict.
    """
    n_h, n_l = X_high.shape[0], X_low.shape[0]
    log(f"  Running pipeline: {n_h} HIGH vs {n_l} LOW, {len(gene_names)} genes")

    # DE: Wilcoxon
    banner(f"  DE ({label})")
    stats, pvals = [], []
    for g in range(len(gene_names)):
        if X_high[:,g].std() == 0 and X_low[:,g].std() == 0:
            stats.append(0); pvals.append(1); continue
        s, p = ranksums(X_high[:,g], X_low[:,g])
        stats.append(s); pvals.append(p)
    de_mask = np.array(pvals) < 0.05
    de_genes = gene_names[de_mask]
    de_idx = np.where(de_mask)[0]
    log(f"  DE genes (p<0.05): {len(de_genes)}")

    if len(de_genes) < 50:
        log("  WARNING: Very few DE genes.")
        return {'label': label, 'n_high': n_h, 'n_low': n_l,
                'n_de_genes': len(de_genes), 'n_nodes': 0, 'n_edges': 0,
                'n_communities': 0, 'A3A_in_network': False, 'A3A_degree': 0,
                'A3A_community': -1, 'A3B_in_network': False,
                'interactors_found': '[]'}

    X_de_h = X_high[:, de_idx]
    X_de_l = X_low[:, de_idx]

    # Spearman
    banner(f"  SPEARMAN ({label})")
    n_g = len(de_genes)
    log(f"  Computing {n_g}x{n_g} matrices...")
    Xhr = np.apply_along_axis(rankdata, 0, X_de_h)
    Xlr = np.apply_along_axis(rankdata, 0, X_de_l)
    rho_h = np.corrcoef(Xhr.T); rho_l = np.corrcoef(Xlr.T)
    np.fill_diagonal(rho_h, 0); np.fill_diagonal(rho_l, 0)
    rho_h = np.nan_to_num(rho_h); rho_l = np.nan_to_num(rho_l)
    diff = rho_h - rho_l
    log(f"  DIFF range: [{diff.min():.3f}, {diff.max():.3f}]")

    # Check A3A max DIFF (diagnostic for threshold sensitivity)
    if 'APOBEC3A' in de_genes:
        a3a_idx = list(de_genes).index('APOBEC3A')
        a3a_max_diff = np.max(np.abs(diff[a3a_idx, :]))
        log(f"  A3A max |DIFF|: {a3a_max_diff:.4f} (threshold: {DIFF_THRESHOLD})")

    # DIFF network
    banner(f"  DIFF NETWORK ({label})")
    G = nx.Graph()
    G.add_nodes_from(de_genes)
    for i in range(n_g):
        for j in range(i+1, n_g):
            if abs(diff[i,j]) >= DIFF_THRESHOLD:
                G.add_edge(de_genes[i], de_genes[j], weight=diff[i,j])
    isolates = list(nx.isolates(G)); G.remove_nodes_from(isolates)
    log(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    if G.number_of_nodes() < 10:
        log("  Network too small")
        return {'label': label, 'n_high': n_h, 'n_low': n_l,
                'n_de_genes': len(de_genes), 'n_nodes': G.number_of_nodes(),
                'n_edges': G.number_of_edges(), 'n_communities': 0,
                'A3A_in_network': 'APOBEC3A' in G.nodes(), 'A3A_degree': 0,
                'A3A_community': -1, 'A3B_in_network': 'APOBEC3B' in G.nodes(),
                'interactors_found': '[]'}

    # Leiden
    banner(f"  LEIDEN ({label})")
    node_list = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(node_list)}
    ig_edges = [(node_idx[u], node_idx[v]) for u, v in G.edges()]
    ig_weights = [abs(G[u][v]['weight']) for u, v in G.edges()]
    ig_g = ig.Graph(n=len(node_list), edges=ig_edges)
    ig_g.es['weight'] = ig_weights
    part = leidenalg.find_partition(ig_g, leidenalg.RBConfigurationVertexPartition,
                                    weights='weight', resolution_parameter=LEIDEN_RESOLUTION,
                                    seed=RANDOM_SEED)
    communities = {}
    for ni, c in enumerate(part.membership):
        communities.setdefault(c, []).append(node_list[ni])
    communities = {k: v for k, v in communities.items() if len(v) >= 5}
    log(f"  Communities: {len(communities)}")

    a3a_comm, a3a_deg = -1, 0
    if 'APOBEC3A' in G.nodes():
        a3a_deg = G.degree('APOBEC3A')
        for c, gl in communities.items():
            if 'APOBEC3A' in gl: a3a_comm = c; break
        log(f"  A3A: comm={a3a_comm}, deg={a3a_deg}")
    else:
        log("  A3A NOT in network")

    interactors = [(g, c) for g in C0_INTERACTORS if g in G.nodes()
                   for c, gl in communities.items() if g in gl]
    log(f"  Interactors: {interactors}")

    result = {
        'label': label, 'n_high': n_h, 'n_low': n_l,
        'n_de_genes': len(de_genes), 'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(), 'n_communities': len(communities),
        'total_community_genes': sum(len(v) for v in communities.values()),
        'A3A_in_network': 'APOBEC3A' in G.nodes(),
        'A3A_degree': a3a_deg, 'A3A_community': a3a_comm,
        'A3A_community_size': len(communities.get(a3a_comm, [])),
        'A3B_in_network': 'APOBEC3B' in G.nodes(),
        'interactors_found': str(interactors),
    }
    return result, communities


def main():
    banner("TIER 3B: TUMOR-ONLY SENSITIVITY (v2 — both approaches)")
    run_dir = ensure_dir(os.path.join(DIR_03_SENSITIVITY, "tumor_only"))

    adata = load_adata()
    high_cells, low_cells = load_groups()
    high_in = high_cells & set(adata.obs_names)
    low_in = low_cells & set(adata.obs_names)

    # Identify tumor vs normal HIGH cells
    tissue_map = adata.obs[TISSUE_COL].to_dict()
    high_tumor = {c for c in high_in if tissue_map.get(c) == 'tumor'}
    high_normal = high_in - high_tumor
    log(f"  Original HIGH: {len(high_in)}")
    log(f"  Tumor HIGH: {len(high_tumor)}")
    log(f"  Normal HIGH removed: {len(high_normal)}")

    # Get basal subset for expression extraction
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell']
    genes = basal.var_names.values

    def get_expression(cell_set):
        cells = [c for c in cell_set if c in basal.obs_names]
        sub = basal[cells]
        X = sub.X
        if scipy.sparse.issparse(X): X = X.toarray()
        return X, cells

    # =================================================================
    # APPROACH A: Tumor HIGH + ORIGINAL LOWs (unbalanced)
    # =================================================================
    banner("APPROACH A: TUMOR HIGH + ORIGINAL LOWS (unbalanced)")
    log(f"  521 tumor HIGH vs 546 original LOW")

    X_high_a, cells_h_a = get_expression(high_tumor)
    X_low_a, cells_l_a = get_expression(low_in)
    log(f"  HIGH matrix: {X_high_a.shape}")
    log(f"  LOW matrix: {X_low_a.shape}")

    result_a = run_network_pipeline(X_high_a, X_low_a, genes, "tumor_HIGH_original_LOW")
    if isinstance(result_a, tuple):
        result_a, comms_a = result_a
        comm_rows = [{'gene': g, 'community': c} for c, gl in comms_a.items() for g in gl]
        pd.DataFrame(comm_rows).to_csv(
            os.path.join(run_dir, "approach_A_communities.tsv"), sep='\t', index=False)
    result_a['approach'] = 'A_original_LOW'
    result_a['n_normal_removed'] = len(high_normal)

    # =================================================================
    # APPROACH B: Tumor HIGH + RE-MATCHED TUMOR LOWs (balanced)
    # =================================================================
    banner("APPROACH B: TUMOR HIGH + RE-MATCHED TUMOR LOWS (balanced)")
    log(f"  521 tumor HIGH vs 521 re-matched tumor LOW")

    # Select new LOWs from tumor-only basal cells with SBS2=0
    basal_tumor = basal[basal.obs[TISSUE_COL] == 'tumor']
    low_candidates = [c for c in basal_tumor.obs_names if c not in high_tumor]

    weights = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    sbs2_vals = {c: weights.loc['SBS2', c] if c in weights.columns else 0 for c in low_candidates}
    zero_sbs2 = [c for c, v in sbs2_vals.items() if v == 0]

    np.random.seed(RANDOM_SEED)
    if len(zero_sbs2) >= len(high_tumor):
        low_tumor = set(np.random.choice(zero_sbs2, size=len(high_tumor), replace=False))
    else:
        remaining = sorted([(c, sbs2_vals[c]) for c in low_candidates if c not in zero_sbs2],
                          key=lambda x: x[1])
        needed = len(high_tumor) - len(zero_sbs2)
        low_tumor = set(zero_sbs2) | set([c for c, _ in remaining[:needed]])
    log(f"  Re-matched tumor LOW: {len(low_tumor)}")

    X_high_b, cells_h_b = get_expression(high_tumor)
    X_low_b, cells_l_b = get_expression(low_tumor)
    log(f"  HIGH matrix: {X_high_b.shape}")
    log(f"  LOW matrix: {X_low_b.shape}")

    result_b = run_network_pipeline(X_high_b, X_low_b, genes, "tumor_HIGH_rematched_LOW")
    if isinstance(result_b, tuple):
        result_b, comms_b = result_b
        comm_rows = [{'gene': g, 'community': c} for c, gl in comms_b.items() for g in gl]
        pd.DataFrame(comm_rows).to_csv(
            os.path.join(run_dir, "approach_B_communities.tsv"), sep='\t', index=False)
    result_b['approach'] = 'B_rematched_LOW'
    result_b['n_normal_removed'] = len(high_normal)

    # =================================================================
    # COMPARISON
    # =================================================================
    banner("COMPARISON: APPROACH A vs B")
    comparison = pd.DataFrame([result_a, result_b])
    comp_path = os.path.join(run_dir, "tumor_only_comparison.tsv")
    comparison.to_csv(comp_path, sep='\t', index=False)
    log(f"  Saved: {comp_path}")

    log(f"\n  {'Metric':<30s}  {'A (orig LOW)':>15s}  {'B (rematched)':>15s}")
    log(f"  {'-'*30}  {'-'*15}  {'-'*15}")
    for col in ['n_high', 'n_low', 'n_de_genes', 'n_nodes', 'n_edges',
                'n_communities', 'A3A_in_network', 'A3A_degree', 'A3A_community']:
        va = result_a.get(col, '')
        vb = result_b.get(col, '')
        log(f"  {col:<30s}  {str(va):>15s}  {str(vb):>15s}")

    log(f"\n  A interactors: {result_a.get('interactors_found','')}")
    log(f"  B interactors: {result_b.get('interactors_found','')}")

    if result_a.get('A3A_in_network') and not result_b.get('A3A_in_network'):
        log(f"\n  FINDING: A3A present with original LOWs but absent with re-matched LOWs.")
        log(f"  This confirms the v1 dropout was caused by control re-selection,")
        log(f"  not by removing 25 normal-tissue HIGH cells.")
    elif result_a.get('A3A_in_network') and result_b.get('A3A_in_network'):
        log(f"\n  FINDING: A3A present in both approaches. Network is robust.")
    elif not result_a.get('A3A_in_network') and not result_b.get('A3A_in_network'):
        log(f"\n  FINDING: A3A absent in both approaches. The 25 normal HIGH cells")
        log(f"  may have been contributing disproportionately to A3A's DIFF signal.")

    log("\nTier 3B complete.")


if __name__ == "__main__":
    main()
