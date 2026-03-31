#!/usr/bin/env python3
"""
Tier3A_Leave_One_Patient_Out.py
================================

Leave-one-patient-out (LOPO) network sensitivity analysis.

For each excluded patient (default: SC029, SC013, SC001):
  1. Remove that patient's cells from HIGH group
  2. Subsample LOW to match new HIGH size
  3. Re-run: DE -> Spearman correlations -> DIFF -> Leiden -> Centrality
  4. Compare to full analysis

Usage:
  conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC029"
  conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC013"
  conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC001"
"""

import os, sys, argparse, pickle, numpy as np, pandas as pd, scanpy as sc
import scipy.sparse
from scipy.stats import ranksums, spearmanr
from statsmodels.stats.multitest import multipletests
import igraph as ig
import leidenalg
import networkx as nx
import matplotlib; matplotlib.use('Agg')

from patient_config import *

def run_lopo(exclude_patient):
    """Run full network pipeline excluding one patient."""

    short = exclude_patient.replace("Patient ", "")
    banner(f"TIER 3A: LOPO — EXCLUDING {exclude_patient}")
    run_dir = ensure_dir(os.path.join(DIR_03_SENSITIVITY, f"LOPO_{short}"))

    # Load data
    adata = load_adata()
    high_cells, low_cells = load_groups()

    # Tag cells
    high_in = high_cells & set(adata.obs_names)
    low_in = low_cells & set(adata.obs_names)

    # Remove excluded patient from HIGH group
    patient_map = adata.obs[PATIENT_COL].to_dict()
    high_filtered = {c for c in high_in if patient_map.get(c) != exclude_patient}
    n_removed = len(high_in) - len(high_filtered)
    log(f"  Original HIGH: {len(high_in)}")
    log(f"  Removed ({exclude_patient}): {n_removed}")
    log(f"  Remaining HIGH: {len(high_filtered)}")

    if len(high_filtered) < 50:
        log("  ERROR: Too few HIGH cells remaining. Skipping.")
        return None

    # Subsample LOW to match
    np.random.seed(RANDOM_SEED)
    low_list = list(low_in)
    np.random.shuffle(low_list)
    low_filtered = set(low_list[:len(high_filtered)])
    log(f"  Matched LOW: {len(low_filtered)}")

    # Subset to basal cells in these groups
    all_cells = list(high_filtered | low_filtered)
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell']
    subset_cells = [c for c in all_cells if c in basal.obs_names]
    adata_sub = basal[subset_cells].copy()

    # Tag groups
    adata_sub.obs['_group'] = 'LOW'
    adata_sub.obs.loc[adata_sub.obs_names.isin(high_filtered), '_group'] = 'HIGH'
    n_h = (adata_sub.obs['_group'] == 'HIGH').sum()
    n_l = (adata_sub.obs['_group'] == 'LOW').sum()
    log(f"  Final subset: {n_h} HIGH + {n_l} LOW = {adata_sub.n_obs}")

    # Get expression
    X = adata_sub.X
    if scipy.sparse.issparse(X): X = X.toarray()
    genes = adata_sub.var_names.values
    groups = adata_sub.obs['_group'].values

    # =================================================================
    # DE: Wilcoxon
    # =================================================================
    banner(f"  DE (LOPO {short})")
    mask_h = groups == 'HIGH'
    mask_l = groups == 'LOW'
    stats, pvals = [], []
    for g in range(X.shape[1]):
        if X[mask_h,g].std() == 0 and X[mask_l,g].std() == 0:
            stats.append(0); pvals.append(1); continue
        s, p = ranksums(X[mask_h,g], X[mask_l,g])
        stats.append(s); pvals.append(p)
    de_df = pd.DataFrame({'gene': genes, 'stat': stats, 'pval': pvals})
    de_genes = de_df[de_df['pval'] < 0.05]['gene'].tolist()
    log(f"  DE genes (p<0.05): {len(de_genes)}")

    if len(de_genes) < 100:
        log("  WARNING: Very few DE genes. Network may be sparse.")

    # Subset expression to DE genes
    de_idx = [i for i, g in enumerate(genes) if g in set(de_genes)]
    X_de_h = X[mask_h][:, de_idx]
    X_de_l = X[mask_l][:, de_idx]
    de_gene_names = genes[de_idx]

    # =================================================================
    # Spearman correlations
    # =================================================================
    banner(f"  SPEARMAN CORRELATIONS (LOPO {short})")
    n_genes = len(de_gene_names)
    log(f"  Computing {n_genes}x{n_genes} Spearman matrices...")

    rho_h = np.corrcoef(X_de_h.T)  # Pearson on ranks approx
    rho_l = np.corrcoef(X_de_l.T)

    # Actually use Spearman
    from scipy.stats import rankdata
    X_de_h_ranked = np.apply_along_axis(rankdata, 0, X_de_h)
    X_de_l_ranked = np.apply_along_axis(rankdata, 0, X_de_l)
    rho_h = np.corrcoef(X_de_h_ranked.T)
    rho_l = np.corrcoef(X_de_l_ranked.T)

    np.fill_diagonal(rho_h, 0)
    np.fill_diagonal(rho_l, 0)
    rho_h = np.nan_to_num(rho_h)
    rho_l = np.nan_to_num(rho_l)

    diff = rho_h - rho_l
    log(f"  DIFF range: [{diff.min():.3f}, {diff.max():.3f}]")

    # =================================================================
    # Build DIFF network
    # =================================================================
    banner(f"  DIFF NETWORK (LOPO {short})")
    G = nx.Graph()
    G.add_nodes_from(de_gene_names)
    n_edges = 0
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            if abs(diff[i,j]) >= DIFF_THRESHOLD:
                G.add_edge(de_gene_names[i], de_gene_names[j], weight=diff[i,j])
                n_edges += 1
    log(f"  DIFF network: {G.number_of_nodes()} nodes, {n_edges} edges")

    # Remove isolates
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    log(f"  After removing isolates: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    if G.number_of_nodes() < 10:
        log("  WARNING: Network too small for community detection.")
        result = {
            'excluded': exclude_patient, 'n_high': n_h, 'n_low': n_l,
            'n_de_genes': len(de_genes), 'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(), 'n_communities': 0,
            'A3A_in_network': 'APOBEC3A' in G.nodes(),
            'A3A_degree': 0, 'A3A_community': -1,
        }
        return result

    # =================================================================
    # Leiden community detection
    # =================================================================
    banner(f"  LEIDEN COMMUNITIES (LOPO {short})")

    # Convert to igraph
    node_list = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(node_list)}
    ig_edges = [(node_idx[u], node_idx[v]) for u, v in G.edges()]
    ig_weights = [abs(G[u][v]['weight']) for u, v in G.edges()]

    ig_graph = ig.Graph(n=len(node_list), edges=ig_edges)
    ig_graph.es['weight'] = ig_weights

    partition = leidenalg.find_partition(
        ig_graph, leidenalg.RBConfigurationVertexPartition,
        weights='weight', resolution_parameter=LEIDEN_RESOLUTION, seed=RANDOM_SEED
    )

    membership = partition.membership
    communities = {}
    for node_i, comm in enumerate(membership):
        communities.setdefault(comm, []).append(node_list[node_i])

    # Filter small communities
    communities = {k: v for k, v in communities.items() if len(v) >= 5}
    log(f"  Communities (size >= 5): {len(communities)}")
    total_genes = sum(len(v) for v in communities.values())
    log(f"  Total community genes: {total_genes}")

    # Check A3A
    a3a_comm = -1
    a3a_degree = 0
    if 'APOBEC3A' in G.nodes():
        a3a_degree = G.degree('APOBEC3A')
        for c, genes_list in communities.items():
            if 'APOBEC3A' in genes_list:
                a3a_comm = c
                break
        log(f"  APOBEC3A: community={a3a_comm}, degree={a3a_degree}")
    else:
        log(f"  APOBEC3A: NOT in network")

    # Check interactors
    interactors_found = []
    for gene in C0_INTERACTORS:
        if gene in G.nodes():
            for c, genes_list in communities.items():
                if gene in genes_list:
                    interactors_found.append((gene, c))
                    break
    log(f"  C0 interactors in network: {interactors_found}")

    # Check A3B
    a3b_in = 'APOBEC3B' in G.nodes()
    log(f"  APOBEC3B in network: {a3b_in}")

    result = {
        'excluded': exclude_patient, 'n_high': n_h, 'n_low': n_l,
        'n_de_genes': len(de_genes), 'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(), 'n_communities': len(communities),
        'total_community_genes': total_genes,
        'A3A_in_network': 'APOBEC3A' in G.nodes(),
        'A3A_degree': a3a_degree, 'A3A_community': a3a_comm,
        'A3A_community_size': len(communities.get(a3a_comm, [])),
        'A3B_in_network': a3b_in,
        'interactors_found': str(interactors_found),
    }

    # Save per-run details
    pd.DataFrame([result]).to_csv(
        os.path.join(run_dir, f"LOPO_{short}_summary.tsv"), sep='\t', index=False)

    # Save community assignments
    comm_rows = []
    for c, genes_list in communities.items():
        for g in genes_list:
            comm_rows.append({'gene': g, 'community': c})
    pd.DataFrame(comm_rows).to_csv(
        os.path.join(run_dir, f"LOPO_{short}_communities.tsv"), sep='\t', index=False)

    log(f"\n  LOPO {short} complete. Saved to {run_dir}")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exclude', type=str, required=True,
                       help='Patient to exclude (e.g. "Patient SC029")')
    args = parser.parse_args()

    result = run_lopo(args.exclude)
    if result:
        banner("RESULT SUMMARY")
        for k, v in result.items():
            log(f"  {k}: {v}")

if __name__ == "__main__":
    main()
