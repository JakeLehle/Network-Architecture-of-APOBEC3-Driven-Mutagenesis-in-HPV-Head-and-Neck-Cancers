#!/usr/bin/env python3
"""
Step02_SC_Correlation_Networks.py
==================================

Figure 4 — Step 02: Build Spearman co-expression networks for SBS2-HIGH
and SBS2-LOW basal cell groups, compute the differential network (DIFF),
and apply edge thresholds.

Mirrors Step04_Correlation_Networks.py from Figure 2 (TCGA bulk).

Workflow:
  1. Load expression matrices + selected gene list from Step 01
  2. Subset to network-ready genes
  3. Compute Spearman correlation for HIGH group (genes × genes)
  4. Compute Spearman correlation for LOW group
  5. Compute DIFF = HIGH − LOW
  6. Apply thresholds and build NetworkX graphs
  7. Remove isolated nodes from DIFF
  8. Auto-report: network statistics at multiple thresholds
  9. Save correlation matrices, edge lists, and graph objects

Input:
  data/FIG_4/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv
  data/FIG_4/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv
  data/FIG_4/02_differential_expression/SC_selected_genes_filtered.csv

Output (→ data/FIG_4/03_correlation_networks/):
  corr_matrices/   — HIGH, LOW, DIFF correlation matrices (pickle + CSV)
  heatmaps/        — correlation distribution plots
  edge_lists/      — TSV edge lists for external visualization
  graph_objects/   — NetworkX graph pickles
  threshold_report.txt — multi-threshold auto-report

Usage:
  conda run -n NETWORK python Step02_SC_Correlation_Networks.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime

from network_config_SC import (
    DIR_01_GROUPS, DIR_02_DE, DIR_03_NETWORKS,
    HIGH_EXPR_PATH, LOW_EXPR_PATH,
    A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    CORRELATION_METHOD, CORR_THRESHOLD, DIFF_THRESHOLD,
    SWEEP_THRESHOLDS,
    banner, log, ensure_dir, load_ensg_to_symbol
)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def compute_spearman_matrix(df, gene_list):
    """
    Compute gene-gene Spearman correlation matrix.

    Parameters
    ----------
    df : pd.DataFrame
        Expression matrix (genes × cells), subset to gene_list.
    gene_list : list
        Ordered list of genes (row indices).

    Returns
    -------
    pd.DataFrame
        genes × genes Spearman correlation matrix.
    """
    # spearmanr expects (n_obs, n_vars) — transpose so cells=rows, genes=cols
    mat = df.loc[gene_list].values.T  # (n_cells, n_genes)
    corr, _ = spearmanr(mat)

    # Handle edge case: single gene returns scalar
    if np.isscalar(corr):
        corr = np.array([[1.0]])

    corr_df = pd.DataFrame(corr, index=gene_list, columns=gene_list)
    return corr_df


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


def graph_stats(G, label=""):
    """Compute and log basic graph statistics."""
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    n_iso = sum(1 for n in G.nodes() if G.degree(n) == 0)
    avg_deg = 2 * n_edges / n_nodes if n_nodes > 0 else 0
    density = nx.density(G) if n_nodes > 1 else 0

    if label:
        log(f"  {label}:")
    log(f"    Nodes: {n_nodes} (isolated: {n_iso})")
    log(f"    Edges: {n_edges}")
    log(f"    Avg degree: {avg_deg:.2f}")
    log(f"    Density: {density:.6f}")

    return {
        "nodes": n_nodes, "edges": n_edges, "isolated": n_iso,
        "avg_degree": avg_deg, "density": density
    }


def remove_isolated(G):
    """Return copy of G with degree-0 nodes removed."""
    G2 = G.copy()
    isolates = [n for n in G2.nodes() if G2.degree(n) == 0]
    G2.remove_nodes_from(isolates)
    return G2


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    start_time = datetime.now()
    banner("[STEP 02] Single-Cell Correlation Networks")
    log(f"Start time: {start_time}")

    # Create output subdirectories
    corr_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "corr_matrices"))
    heat_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "heatmaps"))
    edge_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "edge_lists"))
    graph_dir = ensure_dir(os.path.join(DIR_03_NETWORKS, "graph_objects"))

    # =========================================================================
    # 1. Load selected genes
    # =========================================================================
    banner("[STEP 02.1] Load selected genes from Step 01")

    sel_path = os.path.join(DIR_02_DE, "SC_selected_genes_filtered.csv")
    sel_df = pd.read_csv(sel_path)
    gene_list = sorted(sel_df["gene"].tolist())
    log(f"  Network-ready genes: {len(gene_list)}")

    # Report A3 genes in list
    for a3 in A3_GENES:
        alias = A3_ID_TO_ALIAS.get(a3, a3)
        status = "PRESENT" if a3 in gene_list else "ABSENT"
        log(f"  {alias}: {status}")

    # =========================================================================
    # 2. Load expression matrices and subset
    # =========================================================================
    banner("[STEP 02.2] Load expression matrices")

    log(f"  Reading HIGH: {HIGH_EXPR_PATH}")
    high_df = pd.read_csv(HIGH_EXPR_PATH, sep="\t", index_col=0)
    log(f"  Reading LOW: {LOW_EXPR_PATH}")
    low_df = pd.read_csv(LOW_EXPR_PATH, sep="\t", index_col=0)

    # Subset to selected genes
    high_df = high_df.loc[gene_list]
    low_df = low_df.loc[gene_list]
    log(f"  HIGH subset: {high_df.shape[0]} genes × {high_df.shape[1]} cells")
    log(f"  LOW subset: {low_df.shape[0]} genes × {low_df.shape[1]} cells")

    # =========================================================================
    # 3. Compute Spearman correlation matrices
    # =========================================================================
    banner("[STEP 02.3] Compute Spearman correlation — HIGH group")
    log(f"  Computing {len(gene_list)} × {len(gene_list)} matrix ({high_df.shape[1]} cells)...")
    corr_high = compute_spearman_matrix(high_df, gene_list)
    log(f"  HIGH correlation matrix: {corr_high.shape}")

    banner("[STEP 02.4] Compute Spearman correlation — LOW group")
    log(f"  Computing {len(gene_list)} × {len(gene_list)} matrix ({low_df.shape[1]} cells)...")
    corr_low = compute_spearman_matrix(low_df, gene_list)
    log(f"  LOW correlation matrix: {corr_low.shape}")

    # =========================================================================
    # 4. Compute DIFF = HIGH − LOW
    # =========================================================================
    banner("[STEP 02.5] Compute DIFF correlation matrix (HIGH − LOW)")
    corr_diff = corr_high - corr_low
    log(f"  DIFF matrix: {corr_diff.shape}")

    # Distribution statistics
    upper_tri = corr_diff.values[np.triu_indices_from(corr_diff.values, k=1)]
    log(f"  DIFF upper-triangle stats:")
    log(f"    Mean:   {np.mean(upper_tri):.6f}")
    log(f"    Std:    {np.std(upper_tri):.6f}")
    log(f"    Min:    {np.min(upper_tri):.6f}")
    log(f"    Max:    {np.max(upper_tri):.6f}")
    log(f"    |val| > 0.5: {np.sum(np.abs(upper_tri) > 0.5)}")
    log(f"    |val| > 0.7: {np.sum(np.abs(upper_tri) > 0.7)}")
    log(f"    |val| > 0.8: {np.sum(np.abs(upper_tri) > 0.8)}")

    # =========================================================================
    # 5. Save correlation matrices
    # =========================================================================
    banner("[STEP 02.6] Save correlation matrices")

    for name, mat in [("HIGH", corr_high), ("LOW", corr_low), ("DIFF", corr_diff)]:
        pkl_path = os.path.join(corr_dir, f"SC_corr_{name}.pkl")
        csv_path = os.path.join(corr_dir, f"SC_corr_{name}.csv.gz")
        with open(pkl_path, "wb") as f:
            pickle.dump(mat, f)
        mat.to_csv(csv_path, compression="gzip")
        log(f"  [SAVE] {name} → {pkl_path}")

    # =========================================================================
    # 6. Build graphs at primary thresholds
    # =========================================================================
    banner("[STEP 02.7] Build graphs at primary thresholds")

    log(f"\n  HIGH network (|rho| >= {CORR_THRESHOLD}):")
    G_high = build_weighted_graph(corr_high, CORR_THRESHOLD)
    graph_stats(G_high, "G_high")

    log(f"\n  LOW network (|rho| >= {CORR_THRESHOLD}):")
    G_low = build_weighted_graph(corr_low, CORR_THRESHOLD)
    graph_stats(G_low, "G_low")

    log(f"\n  DIFF network (|delta-rho| >= {DIFF_THRESHOLD}):")
    G_diff = build_weighted_graph(corr_diff, DIFF_THRESHOLD)
    graph_stats(G_diff, "G_diff (with isolates)")

    G_diff_noiso = remove_isolated(G_diff)
    graph_stats(G_diff_noiso, "G_diff (no isolates)")

    G_high_noiso = remove_isolated(G_high)
    G_low_noiso = remove_isolated(G_low)

    # LCC stats for DIFF
    if G_diff_noiso.number_of_nodes() > 0:
        components = list(nx.connected_components(G_diff_noiso))
        lcc = max(components, key=len)
        log(f"\n  DIFF connected components: {len(components)}")
        log(f"  DIFF LCC size: {len(lcc)} nodes")
    else:
        log(f"\n  WARNING: DIFF network has no edges at threshold {DIFF_THRESHOLD}")
        log(f"  Consider relaxing DIFF_THRESHOLD (see threshold report below)")

    # =========================================================================
    # 7. Save graphs
    # =========================================================================
    banner("[STEP 02.8] Save graph objects")

    graphs = {
        "SC_G_high": G_high, "SC_G_low": G_low, "SC_G_diff": G_diff,
        "SC_G_high_noiso": G_high_noiso, "SC_G_low_noiso": G_low_noiso,
        "SC_G_diff_noiso": G_diff_noiso,
    }
    for name, G in graphs.items():
        path = os.path.join(graph_dir, f"{name}.gpickle")
        with open(path, "wb") as f:
            pickle.dump(G, f)
        log(f"  [SAVE] {name} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges) → {path}")

    # =========================================================================
    # 8. Save edge lists
    # =========================================================================
    banner("[STEP 02.9] Save edge lists")

    ensg_to_symbol = load_ensg_to_symbol()

    for name, G in [("HIGH", G_high_noiso), ("LOW", G_low_noiso), ("DIFF", G_diff_noiso)]:
        if G.number_of_edges() == 0:
            log(f"  {name}: no edges, skipping edge list")
            continue

        rows = []
        for u, v, d in G.edges(data=True):
            u_sym = ensg_to_symbol.get(u, A3_ID_TO_ALIAS.get(u, u))
            v_sym = ensg_to_symbol.get(v, A3_ID_TO_ALIAS.get(v, v))
            rows.append({
                "source": u, "target": v,
                "source_symbol": u_sym, "target_symbol": v_sym,
                "weight": d.get("weight", 0), "abs_weight": d.get("abs_weight", 0)
            })
        edge_df = pd.DataFrame(rows)
        edge_path = os.path.join(edge_dir, f"SC_edges_{name}.tsv")
        edge_df.to_csv(edge_path, sep="\t", index=False)
        log(f"  [SAVE] {name} edges ({len(edge_df)}) → {edge_path}")

    # =========================================================================
    # 9. Correlation distribution plots
    # =========================================================================
    banner("[STEP 02.10] Correlation distribution plots")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for ax, (name, mat, thresh) in zip(axes, [
        ("HIGH", corr_high, CORR_THRESHOLD),
        ("LOW", corr_low, CORR_THRESHOLD),
        ("DIFF", corr_diff, DIFF_THRESHOLD)
    ]):
        vals = mat.values[np.triu_indices_from(mat.values, k=1)]
        ax.hist(vals, bins=100, color="steelblue", alpha=0.7, edgecolor="none")
        ax.axvline(thresh, ls="--", c="red", lw=1, label=f"+{thresh}")
        ax.axvline(-thresh, ls="--", c="red", lw=1, label=f"-{thresh}")
        ax.set_title(f"{name} (n={len(vals):,})")
        ax.set_xlabel("Correlation" if name != "DIFF" else "Delta-rho")
        ax.legend(fontsize=7)

    plt.tight_layout()
    hist_path = os.path.join(heat_dir, "SC_correlation_distributions.png")
    plt.savefig(hist_path, dpi=300)
    plt.close()
    log(f"  [SAVE] Distributions → {hist_path}")

    # =========================================================================
    # 10. AUTO-REPORT: Multi-threshold sweep
    # =========================================================================
    banner("[STEP 02.11] Multi-threshold diagnostic report")

    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("DIFF NETWORK THRESHOLD SWEEP — AUTO-REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"Genes in network: {len(gene_list)}")
    report_lines.append(f"Cells per group: HIGH={high_df.shape[1]}, LOW={low_df.shape[1]}")
    report_lines.append(f"Primary threshold: |delta-rho| >= {DIFF_THRESHOLD}")
    report_lines.append("")
    report_lines.append(f"{'Threshold':>10} {'Nodes':>8} {'Edges':>8} {'AvgDeg':>8} "
                        f"{'Density':>10} {'Components':>12} {'LCC':>6}")
    report_lines.append("-" * 80)

    sweep_results = []
    for dt in SWEEP_THRESHOLDS:
        G_sweep = build_weighted_graph(corr_diff, dt)
        G_sweep_noiso = remove_isolated(G_sweep)

        n_nodes = G_sweep_noiso.number_of_nodes()
        n_edges = G_sweep_noiso.number_of_edges()
        avg_deg = 2 * n_edges / n_nodes if n_nodes > 0 else 0
        density = nx.density(G_sweep_noiso) if n_nodes > 1 else 0

        if n_nodes > 0:
            components = list(nx.connected_components(G_sweep_noiso))
            n_comp = len(components)
            lcc_size = len(max(components, key=len))
        else:
            n_comp = 0
            lcc_size = 0

        marker = " <<<" if dt == DIFF_THRESHOLD else ""
        report_lines.append(f"{dt:>10.2f} {n_nodes:>8} {n_edges:>8} {avg_deg:>8.1f} "
                            f"{density:>10.6f} {n_comp:>12} {lcc_size:>6}{marker}")

        sweep_results.append({
            "threshold": dt, "nodes": n_nodes, "edges": n_edges,
            "avg_degree": avg_deg, "density": density,
            "components": n_comp, "lcc_size": lcc_size
        })

    report_lines.append("")
    report_lines.append("<<< = primary threshold from config")

    # Print to stdout
    for line in report_lines:
        print(line)

    # Save report
    report_path = os.path.join(DIR_03_NETWORKS, "threshold_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    log(f"\n  [SAVE] Threshold report → {report_path}")

    # Save sweep as CSV
    sweep_df = pd.DataFrame(sweep_results)
    sweep_csv = os.path.join(DIR_03_NETWORKS, "threshold_sweep.csv")
    sweep_df.to_csv(sweep_csv, index=False)
    log(f"  [SAVE] Sweep CSV → {sweep_csv}")

    # =========================================================================
    # DONE
    # =========================================================================
    elapsed = datetime.now() - start_time
    banner(f"[STEP 02 COMPLETE] Elapsed: {elapsed}")

    log(f"\nPrimary DIFF network: {G_diff_noiso.number_of_nodes()} nodes, "
        f"{G_diff_noiso.number_of_edges()} edges")
    log(f"Review threshold_report.txt to decide if thresholds need adjustment")
    log(f"Then run Step03_SC_Community_Detection.py")


if __name__ == "__main__":
    main()
