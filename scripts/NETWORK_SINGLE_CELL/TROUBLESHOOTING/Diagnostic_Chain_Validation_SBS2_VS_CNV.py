#!/usr/bin/env python3
"""
Diagnostic_Chain_Validation_SBS2_VS_CNV.py
=============================================

Confirmatory analysis: do the concordant chain genes identified in the
tumor-vs-normal comparisons show predicted behavior when comparing the
two tumor states directly?

Predictions:
  - Activating chain genes (from SBS2_VS_NORMAL): should show POSITIVE
    DIFF correlations among themselves in SBS2 vs CNV, because their
    coordinated co-expression is specific to SBS2-HIGH.
  - Inhibiting chain genes (from CNV_VS_NORMAL): should show NEGATIVE
    DIFF correlations, because their repressor coordination is specific
    to CNV-HIGH.
  - A3A/A3B: should maintain all-negative DIFF (pulsatile decoupling
    applies in both tumor states).

Method:
  1. Recompute chain membership from the two normal-comparison networks
  2. Load expression data for chain genes + A3A/A3B from adata
  3. Compute Spearman correlations in SBS2-HIGH cells and CNV-HIGH cells
  4. Compute DIFF = corr_SBS2 - corr_CNV
  5. Test predictions by examining within-chain and cross-chain DIFF

Usage:
    conda run -n NETWORK python Diagnostic_Chain_Validation_SBS2_VS_CNV.py
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")
ADATA_PATH = os.path.join(FIG4_ROOT, "00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(FIG4_ROOT,
                          "01_group_selection/three_group_assignments.tsv")
OUTPUT_DIR = os.path.join(FIG4_ROOT, "DIAGNOSTIC_CHAIN_VALIDATION")
os.makedirs(OUTPUT_DIR, exist_ok=True)

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}
A3_TARGETS = ["APOBEC3A", "APOBEC3B"]


def log(msg):
    print(msg, flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# CHAIN IDENTIFICATION (mirrors figure script logic)
# =============================================================================

def identify_chains(net_dir, harris_genes):
    """Identify activator and repressor chain genes from a network."""

    part_df = pd.read_csv(
        os.path.join(net_dir, "04_communities/SC_best_partition.csv"))
    gene_to_comm = dict(zip(part_df["gene"], part_df["community"]))

    with open(os.path.join(net_dir,
              "04_communities/SC_G_comm.gpickle"), "rb") as f:
        G_comm = pickle.load(f)

    de_df = pd.read_csv(
        os.path.join(net_dir,
                     "02_differential_expression/SC_diffexpr_stats.csv"))
    de_log2fc = dict(zip(de_df["gene"], de_df["log2FC"]))

    # Find A3 community
    a3_comm = None
    for g in A3_TARGETS:
        if g in gene_to_comm:
            a3_comm = gene_to_comm[g]
            break
    if a3_comm is None:
        return set(), set(), set()

    comm_nodes = [g for g, c in gene_to_comm.items() if c == a3_comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
    G_sub.remove_nodes_from(isolates)

    a3_set = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)
    harris_in_comm = set(n for n in G_sub.nodes()
                         if n in harris_genes and n not in A3_SYMBOLS)

    # Build concordance subgraphs
    def build_conc(mode):
        a3_nodes = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)
        G_conc = nx.Graph()
        if mode == "activator":
            eligible = set(n for n in G_sub.nodes()
                           if de_log2fc.get(n, 0) > 0) | a3_nodes
            for u, v, d in G_sub.edges(data=True):
                w = d.get("weight", 0)
                if w > 0 and u in eligible and v in eligible:
                    G_conc.add_edge(u, v, weight=w)
        else:
            eligible = set(n for n in G_sub.nodes()
                           if de_log2fc.get(n, 0) <= 0) | a3_nodes
            for u, v, d in G_sub.edges(data=True):
                w = d.get("weight", 0)
                if w < 0 and u in eligible and v in eligible:
                    G_conc.add_edge(u, v, weight=w)
        return G_conc

    G_act = build_conc("activator")
    G_rep = build_conc("repressor")

    # Find chains containing Harris interactors
    activator_nodes = set()
    repressor_nodes = set()

    for gene in harris_in_comm:
        if de_log2fc.get(gene, 0) > 0 and gene in G_act and \
           G_act.degree(gene) > 0:
            comp = set(nx.node_connected_component(G_act, gene))
            activator_nodes |= comp
        if de_log2fc.get(gene, 0) <= 0 and gene in G_rep and \
           G_rep.degree(gene) > 0:
            comp = set(nx.node_connected_component(G_rep, gene))
            repressor_nodes |= comp

    # Add boundary gene chains
    for a3g in (a3_set & set(G_sub.nodes())):
        for nb in G_sub.neighbors(a3g):
            if nb in A3_SYMBOLS:
                continue
            nb_fc = de_log2fc.get(nb, 0)
            for nn in G_sub.neighbors(nb):
                if nn == a3g or nn in A3_SYMBOLS:
                    continue
                nn_w = G_sub[nb][nn].get("weight", 0)
                nn_fc = de_log2fc.get(nn, 0)
                if nb_fc > 0 and nn_fc > 0 and nn_w > 0:
                    if nb in G_act and G_act.degree(nb) > 0:
                        comp = set(nx.node_connected_component(G_act, nb))
                        activator_nodes |= comp
                    break
                elif nb_fc <= 0 and nn_fc <= 0 and nn_w < 0:
                    if nb in G_rep and G_rep.degree(nb) > 0:
                        comp = set(nx.node_connected_component(G_rep, nb))
                        repressor_nodes |= comp
                    break

    activator_nodes -= a3_set
    repressor_nodes -= a3_set

    return activator_nodes, repressor_nodes, harris_in_comm


# =============================================================================
# COMPUTE CORRELATIONS FROM EXPRESSION
# =============================================================================

def compute_spearman_matrix(adata, cells, genes):
    """Compute pairwise Spearman correlation for genes across cells.

    Returns DataFrame (genes x genes).
    """
    # Subset adata to cells and genes
    valid_genes = [g for g in genes if g in adata.var_names]
    missing = set(genes) - set(valid_genes)
    if missing:
        log(f"    Missing from adata: {', '.join(sorted(missing))}")

    adata_sub = adata[cells][:, valid_genes]

    # Get dense expression matrix (cells x genes)
    if hasattr(adata_sub.X, "toarray"):
        X = adata_sub.X.toarray()
    else:
        X = np.array(adata_sub.X)

    n_genes = len(valid_genes)
    corr_matrix = np.zeros((n_genes, n_genes))

    for i in range(n_genes):
        for j in range(i, n_genes):
            if i == j:
                corr_matrix[i, j] = 1.0
            else:
                rho, _ = spearmanr(X[:, i], X[:, j])
                if np.isnan(rho):
                    rho = 0.0
                corr_matrix[i, j] = rho
                corr_matrix[j, i] = rho

    return pd.DataFrame(corr_matrix, index=valid_genes, columns=valid_genes)


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    banner("CHAIN VALIDATION: SBS2 vs CNV")

    import scanpy as sc

    # Load Harris interactors
    harris_df = pd.read_csv(HARRIS_PATH, sep="\t")
    harris_genes = set(harris_df["gene_symbol"].values)
    log(f"Harris A3 interactors: {len(harris_genes)}")

    # -----------------------------------------------------------------
    # Step 1: Identify chain genes from the two normal comparisons
    # -----------------------------------------------------------------
    banner("STEP 1: Identify chain genes")

    sbs2_dir = os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL")
    cnv_dir = os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL")

    act_sbs2, rep_sbs2, harris_sbs2 = identify_chains(sbs2_dir, harris_genes)
    act_cnv, rep_cnv, harris_cnv = identify_chains(cnv_dir, harris_genes)

    log(f"\n  SBS2_VS_NORMAL chains:")
    log(f"    Activating: {len(act_sbs2)} genes "
        f"({', '.join(sorted(act_sbs2)[:10])})")
    log(f"    Inhibiting: {len(rep_sbs2)} genes "
        f"({', '.join(sorted(rep_sbs2)[:10])})")

    log(f"\n  CNV_VS_NORMAL chains:")
    log(f"    Activating: {len(act_cnv)} genes "
        f"({', '.join(sorted(act_cnv)[:10])})")
    log(f"    Inhibiting: {len(rep_cnv)} genes "
        f"({', '.join(sorted(rep_cnv)[:10])})")

    # Combined gene sets
    all_activating = act_sbs2 | act_cnv
    all_inhibiting = rep_sbs2 | rep_cnv
    a3_genes = set(A3_TARGETS)
    overlap = all_activating & all_inhibiting
    if overlap:
        log(f"\n  [NOTE] Genes in both chains: {', '.join(sorted(overlap))}")
        all_inhibiting -= overlap  # prioritize activating classification

    all_chain_genes = all_activating | all_inhibiting | a3_genes

    log(f"\n  Combined for validation:")
    log(f"    Activating: {len(all_activating)}")
    log(f"    Inhibiting: {len(all_inhibiting)}")
    log(f"    A3 targets: {len(a3_genes)}")
    log(f"    Total unique: {len(all_chain_genes)}")

    # Save gene lists
    gene_list_rows = []
    for g in sorted(all_activating):
        gene_list_rows.append({
            "gene": g, "chain_type": "activating",
            "source": "SBS2" if g in act_sbs2 else "CNV",
            "is_harris": g in harris_genes,
        })
    for g in sorted(all_inhibiting):
        gene_list_rows.append({
            "gene": g, "chain_type": "inhibiting",
            "source": "SBS2" if g in rep_sbs2 else "CNV",
            "is_harris": g in harris_genes,
        })
    for g in sorted(a3_genes):
        gene_list_rows.append({
            "gene": g, "chain_type": "A3",
            "source": "both", "is_harris": False,
        })
    pd.DataFrame(gene_list_rows).to_csv(
        os.path.join(OUTPUT_DIR, "chain_gene_list.tsv"),
        sep="\t", index=False)

    # -----------------------------------------------------------------
    # Step 2: Load expression data and compute correlations
    # -----------------------------------------------------------------
    banner("STEP 2: Compute correlations in SBS2-HIGH and CNV-HIGH cells")

    log("Loading adata...")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

    groups_df = pd.read_csv(GROUP_PATH, sep="\t")
    group_map = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))
    adata.obs["three_group"] = adata.obs.index.map(
        lambda x: group_map.get(x, "OTHER"))

    sbs2_cells = adata.obs.index[adata.obs["three_group"] == "SBS2_HIGH"]
    cnv_cells = adata.obs.index[adata.obs["three_group"] == "CNV_HIGH"]
    log(f"  SBS2-HIGH cells: {len(sbs2_cells)}")
    log(f"  CNV-HIGH cells: {len(cnv_cells)}")

    genes_to_test = sorted(all_chain_genes)
    log(f"  Computing Spearman correlations for {len(genes_to_test)} genes...")

    log(f"\n  SBS2-HIGH correlations:")
    corr_sbs2 = compute_spearman_matrix(adata, sbs2_cells, genes_to_test)

    log(f"  CNV-HIGH correlations:")
    corr_cnv = compute_spearman_matrix(adata, cnv_cells, genes_to_test)

    # Align matrices (in case some genes were missing from adata)
    shared_genes = sorted(set(corr_sbs2.index) & set(corr_cnv.index))
    corr_sbs2 = corr_sbs2.loc[shared_genes, shared_genes]
    corr_cnv = corr_cnv.loc[shared_genes, shared_genes]

    # DIFF = corr_SBS2 - corr_CNV
    diff = corr_sbs2 - corr_cnv
    log(f"\n  DIFF matrix: {diff.shape[0]} x {diff.shape[1]}")

    del adata  # free memory

    # Save matrices
    diff.to_csv(os.path.join(OUTPUT_DIR, "DIFF_chain_genes.tsv"), sep="\t")
    corr_sbs2.to_csv(os.path.join(OUTPUT_DIR, "corr_SBS2_chain_genes.tsv"),
                     sep="\t")
    corr_cnv.to_csv(os.path.join(OUTPUT_DIR, "corr_CNV_chain_genes.tsv"),
                    sep="\t")

    # -----------------------------------------------------------------
    # Step 3: Test predictions
    # -----------------------------------------------------------------
    banner("STEP 3: Test predictions")

    # Classify genes present in the matrix
    act_present = sorted(set(shared_genes) & all_activating)
    rep_present = sorted(set(shared_genes) & all_inhibiting)
    a3_present = sorted(set(shared_genes) & a3_genes)

    log(f"  Activating in matrix: {len(act_present)} "
        f"({', '.join(act_present)})")
    log(f"  Inhibiting in matrix: {len(rep_present)} "
        f"({', '.join(rep_present)})")
    log(f"  A3 in matrix: {len(a3_present)} "
        f"({', '.join(a3_present)})")

    # --- Prediction 1: activating chain internal DIFF ---
    if len(act_present) >= 2:
        act_diffs = []
        for i, g1 in enumerate(act_present):
            for g2 in act_present[i+1:]:
                act_diffs.append(diff.loc[g1, g2])
        act_diffs = np.array(act_diffs)
        log(f"\n  PREDICTION 1: Activating chain internal DIFF")
        log(f"    Expected: positive (co-expression stronger in SBS2)")
        log(f"    Pairs: {len(act_diffs)}")
        log(f"    Mean DIFF: {act_diffs.mean():.4f}")
        log(f"    Median DIFF: {np.median(act_diffs):.4f}")
        log(f"    Positive: {(act_diffs > 0).sum()}/{len(act_diffs)} "
            f"({(act_diffs > 0).mean()*100:.1f}%)")
        log(f"    Range: [{act_diffs.min():.4f}, {act_diffs.max():.4f}]")
    else:
        log(f"\n  PREDICTION 1: SKIP (< 2 activating genes)")

    # --- Prediction 2: inhibiting chain internal DIFF ---
    if len(rep_present) >= 2:
        rep_diffs = []
        for i, g1 in enumerate(rep_present):
            for g2 in rep_present[i+1:]:
                rep_diffs.append(diff.loc[g1, g2])
        rep_diffs = np.array(rep_diffs)
        log(f"\n  PREDICTION 2: Inhibiting chain internal DIFF")
        log(f"    Expected: negative (co-expression stronger in CNV)")
        log(f"    Pairs: {len(rep_diffs)}")
        log(f"    Mean DIFF: {rep_diffs.mean():.4f}")
        log(f"    Median DIFF: {np.median(rep_diffs):.4f}")
        log(f"    Negative: {(rep_diffs < 0).sum()}/{len(rep_diffs)} "
            f"({(rep_diffs < 0).mean()*100:.1f}%)")
        log(f"    Range: [{rep_diffs.min():.4f}, {rep_diffs.max():.4f}]")
    else:
        log(f"\n  PREDICTION 2: SKIP (< 2 inhibiting genes)")

    # --- Prediction 3: cross-chain DIFF (should show separation) ---
    if act_present and rep_present:
        cross_diffs = []
        for g1 in act_present:
            for g2 in rep_present:
                cross_diffs.append(diff.loc[g1, g2])
        cross_diffs = np.array(cross_diffs)
        log(f"\n  PREDICTION 3: Cross-chain DIFF (activating-inhibiting)")
        log(f"    Expected: near zero or mixed (different programs)")
        log(f"    Pairs: {len(cross_diffs)}")
        log(f"    Mean DIFF: {cross_diffs.mean():.4f}")
        log(f"    Median DIFF: {np.median(cross_diffs):.4f}")
        log(f"    Positive: {(cross_diffs > 0).sum()}/{len(cross_diffs)}")
        log(f"    Negative: {(cross_diffs < 0).sum()}/{len(cross_diffs)}")

    # --- Prediction 4: A3 edge profiles ---
    for a3g in a3_present:
        alias = A3_ALIAS.get(a3g, a3g)
        log(f"\n  PREDICTION 4: {alias} edges in SBS2 vs CNV DIFF")
        log(f"    Expected: all or mostly negative (decoupled in both)")

        if act_present:
            a3_act = [diff.loc[a3g, g] for g in act_present
                      if g != a3g]
            if a3_act:
                a3_act = np.array(a3_act)
                log(f"    -> Activating genes: mean={a3_act.mean():.4f}, "
                    f"negative={np.sum(a3_act < 0)}/{len(a3_act)}")

        if rep_present:
            a3_rep = [diff.loc[a3g, g] for g in rep_present
                      if g != a3g]
            if a3_rep:
                a3_rep = np.array(a3_rep)
                log(f"    -> Inhibiting genes: mean={a3_rep.mean():.4f}, "
                    f"negative={np.sum(a3_rep < 0)}/{len(a3_rep)}")

        # All other chain genes
        all_others = [g for g in shared_genes
                      if g != a3g and g not in A3_SYMBOLS]
        a3_all = [diff.loc[a3g, g] for g in all_others]
        if a3_all:
            a3_all = np.array(a3_all)
            log(f"    -> All chain genes: mean={a3_all.mean():.4f}, "
                f"negative={np.sum(a3_all < 0)}/{len(a3_all)}, "
                f"positive={np.sum(a3_all > 0)}/{len(a3_all)}")

    # --- Per-gene detail table ---
    detail_rows = []
    for g in shared_genes:
        if g in A3_SYMBOLS:
            chain_type = "A3"
        elif g in all_activating:
            chain_type = "activating"
        elif g in all_inhibiting:
            chain_type = "inhibiting"
        else:
            chain_type = "unknown"

        # Mean DIFF with activating genes
        act_edges = [diff.loc[g, g2] for g2 in act_present
                     if g2 != g]
        rep_edges = [diff.loc[g, g2] for g2 in rep_present
                     if g2 != g]
        a3_edges = [diff.loc[g, g2] for g2 in a3_present
                    if g2 != g]

        detail_rows.append({
            "gene": g,
            "chain_type": chain_type,
            "is_harris": g in harris_genes,
            "mean_diff_to_activating": (
                round(np.mean(act_edges), 4) if act_edges else None),
            "mean_diff_to_inhibiting": (
                round(np.mean(rep_edges), 4) if rep_edges else None),
            "mean_diff_to_A3": (
                round(np.mean(a3_edges), 4) if a3_edges else None),
        })

    df_detail = pd.DataFrame(detail_rows)
    df_detail.to_csv(
        os.path.join(OUTPUT_DIR, "chain_validation_detail.tsv"),
        sep="\t", index=False)
    log(f"\n  [SAVE] Detail table -> chain_validation_detail.tsv")

    # -----------------------------------------------------------------
    # Step 4: Heatmap of DIFF matrix
    # -----------------------------------------------------------------
    banner("STEP 4: DIFF heatmap")

    # Order genes: activating, inhibiting, A3
    ordered = act_present + rep_present + a3_present
    diff_ordered = diff.loc[ordered, ordered]

    # Build color labels for rows
    row_colors = []
    for g in ordered:
        if g in all_activating:
            row_colors.append("#fcad61")
        elif g in all_inhibiting:
            row_colors.append("#aadda4")
        else:
            row_colors.append("#ed6a5a")

    fig, ax = plt.subplots(figsize=(16, 14))

    # Use diverging colormap centered at 0
    vmax = max(abs(diff_ordered.values.min()),
               abs(diff_ordered.values.max()), 0.3)
    im = ax.imshow(diff_ordered.values, cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="auto")

    # Labels
    display_names = [A3_ALIAS.get(g, g) for g in ordered]
    ax.set_xticks(range(len(ordered)))
    ax.set_xticklabels(display_names, rotation=90, fontsize=10)
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels(display_names, fontsize=10)

    # Color-code tick labels
    for i, (tick, color) in enumerate(
            zip(ax.get_yticklabels(), row_colors)):
        tick.set_color(color)
        tick.set_fontweight("bold")
    for i, (tick, color) in enumerate(
            zip(ax.get_xticklabels(), row_colors)):
        tick.set_color(color)
        tick.set_fontweight("bold")

    # Divider lines between groups
    n_act = len(act_present)
    n_rep = len(rep_present)
    if n_act > 0:
        ax.axhline(n_act - 0.5, color="black", linewidth=2)
        ax.axvline(n_act - 0.5, color="black", linewidth=2)
    if n_rep > 0:
        ax.axhline(n_act + n_rep - 0.5, color="black", linewidth=2)
        ax.axvline(n_act + n_rep - 0.5, color="black", linewidth=2)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("DIFF (corr_SBS2 - corr_CNV)", fontsize=14)

    ax.set_title("Chain Gene DIFF Correlations: SBS2-HIGH vs CNV-HIGH\n"
                 "Orange=Activating, Green=Inhibiting, Red=A3",
                 fontsize=16, pad=15)

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(OUTPUT_DIR,
                    f"DIFF_heatmap_chain_genes.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Heatmap -> DIFF_heatmap_chain_genes.pdf/png")

    banner("CHAIN VALIDATION COMPLETE")
    log(f"\nOutput: {OUTPUT_DIR}")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        log(f"  {f}")


if __name__ == "__main__":
    main()
