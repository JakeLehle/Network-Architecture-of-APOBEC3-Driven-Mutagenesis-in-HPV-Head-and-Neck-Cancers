#!/usr/bin/env python3
"""
Step05_Generate_Figure4_Panels.py
==================================

Figure 4 — Step 05: Generate publication-quality panels for Figure 4, including
the community overlap analysis (SC vs TCGA bulk) and Harris A3B interactor mapping.

NOTE on gene ID namespaces:
  SC communities use gene SYMBOLS (APOBEC3B, KRT5, etc.)
  TCGA communities use ENSG IDs (ENSG00000179750, etc.)
  This script converts TCGA community gene lists from ENSG → symbols using
  the ensg_to_symbol.json mapping from Figure 2, so all overlap computations
  happen in the shared symbol namespace.

Proposed Figure 4 panels:
  a — UMAP: HIGH vs LOW basal cell selection (DONE in Step 00)
  b — Community heatmap or network visualization of SC DIFF communities
  c — Overlap matrix: SC communities vs TCGA bulk communities
      (hypergeometric test, -log10 adj.p heatmap, Jaccard + overlap annotations)
  d — Harris A3B interactor enrichment per community

Input:
  data/FIG_4/04_communities/SC_best_partition.csv
  data/FIG_4/04_communities/SC_community_gene_lists.csv
  data/FIG_4/04_communities/SC_G_comm.gpickle
  data/FIG_4/05_centrality_metrics/SC_DIFF_metrics.csv
  data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv
  data/FIG_4/00_input/Harris_A3_interactors.txt
  data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt

Output (→ data/FIG_4/06_overlap_analysis/ and FIGURE_4_PANELS/):
  overlap_matrix.csv              — raw overlap counts
  hypergeom_pvalues.csv           — p-values per SC×TCGA pair
  hypergeom_pvalues_BH.csv        — BH-adjusted p-values
  jaccard_matrix.csv              — Jaccard index per pair
  harris_enrichment.csv           — Harris interactor enrichment per community
  overlap_gene_details.csv        — per-pair overlapping gene lists
  Panel_4b_*.pdf/.png             — community visualization
  Panel_4c_*.pdf/.png             — overlap heatmap
  Panel_4d_*.pdf/.png             — Harris interactor enrichment

Usage:
  conda run -n NETWORK python Step05_Generate_Figure4_Panels.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import hypergeom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from datetime import datetime

from network_config_SC import (
    DIR_04_COMMUNITIES, DIR_05_CENTRALITY, DIR_06_OVERLAP, FIGURE_4_PANELS,
    DIR_03_NETWORKS, FIG2_COMMUNITIES,
    HARRIS_ALL_PATH, HARRIS_A3B_PATH,
    A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    OVERLAP_FDR_METHOD, OVERLAP_P_THRESHOLD, JACCARD_DISPLAY,
    banner, log, ensure_dir, load_ensg_to_symbol, convert_tcga_genes_to_symbols
)


# =============================================================================
# BH-FDR CORRECTION
# =============================================================================

def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    adjusted = pvals * n / ranked
    adjusted = np.minimum.accumulate(adjusted[np.argsort(ranked)[::-1]])[::-1]
    return np.clip(adjusted, 0, 1)


# =============================================================================
# LOAD COMMUNITY DATA
# =============================================================================

def load_sc_communities():
    """Load single-cell community assignments (already in gene symbol space)."""
    path = os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv")
    df = pd.read_csv(path)
    comm_genes = {}
    for comm in sorted(df["community"].unique()):
        genes = df[df["community"] == comm]["gene"].tolist()
        comm_genes[comm] = set(genes)
    return comm_genes


def load_tcga_communities():
    """
    Load TCGA bulk community assignments from Figure 2 and convert ENSG → symbols.

    Figure 2 communities use ENSG IDs. We convert them to gene symbols so they
    match the SC community namespace for overlap computation.
    """
    cancer_type = "TCGA-HNSC"
    path = os.path.join(FIG2_COMMUNITIES, cancer_type, f"{cancer_type}_best_partition.csv")

    if not os.path.exists(path):
        # Try alternative locations
        alt_path = os.path.join(FIG2_COMMUNITIES, f"{cancer_type}_best_partition.csv")
        if os.path.exists(alt_path):
            path = alt_path
        else:
            log(f"WARNING: TCGA community file not found at {path}")
            log(f"  Also tried: {alt_path}")
            return {}

    log(f"Loading TCGA bulk communities: {path}")
    df = pd.read_csv(path)

    # Figure out column names (may vary)
    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    comm_col = "community" if "community" in df.columns else df.columns[-1]

    # Load ENSG → symbol mapping for conversion
    ensg_to_symbol = load_ensg_to_symbol()
    log(f"  ENSG→symbol mapping: {len(ensg_to_symbol)} entries")

    comm_genes_ensg = {}
    for comm in sorted(df[comm_col].unique()):
        genes = df[df[comm_col] == comm][gene_col].tolist()
        comm_genes_ensg[comm] = set(genes)

    # Convert each community from ENSG → symbols
    comm_genes_symbol = {}
    total_converted = 0
    total_genes = 0
    for comm, ensg_set in comm_genes_ensg.items():
        symbol_set = convert_tcga_genes_to_symbols(ensg_set, ensg_to_symbol)
        comm_genes_symbol[comm] = symbol_set
        # Count how many were actually converted (vs kept as ENSG fallback)
        n_converted = sum(1 for g in symbol_set if not g.startswith("ENSG"))
        total_converted += n_converted
        total_genes += len(ensg_set)

    log(f"  Loaded {len(comm_genes_symbol)} TCGA communities")
    log(f"  Converted {total_converted}/{total_genes} genes from ENSG → symbol "
        f"({100*total_converted/total_genes:.1f}%)")

    return comm_genes_symbol


def load_harris_interactors():
    """Load Harris A3 interactor gene lists (gene symbols — direct match to SC)."""
    interactors = {"all": set(), "A3B_only": set()}

    if os.path.exists(HARRIS_ALL_PATH):
        with open(HARRIS_ALL_PATH) as f:
            for line in f:
                gene = line.strip().split("\t")[0].strip()
                if gene and not gene.startswith("#"):
                    interactors["all"].add(gene)
        log(f"  Harris ALL interactors: {len(interactors['all'])}")

    if os.path.exists(HARRIS_A3B_PATH):
        with open(HARRIS_A3B_PATH) as f:
            for line in f:
                gene = line.strip().split("\t")[0].strip()
                if gene and not gene.startswith("#"):
                    interactors["A3B_only"].add(gene)
        log(f"  Harris A3B-only interactors: {len(interactors['A3B_only'])}")

    return interactors


# =============================================================================
# HYPERGEOMETRIC OVERLAP TEST
# =============================================================================

def compute_overlap_matrix(sc_comms, tcga_comms, universe_size):
    """
    Compute overlap statistics between SC and TCGA communities.
    Both are now in gene symbol space.

    For each pair (SC_i, TCGA_j):
      - k = |SC_i ∩ TCGA_j| (overlap count)
      - M = universe_size
      - n = |SC_i|
      - N = |TCGA_j|
      - p = P(X >= k) under Hypergeometric(M, N, n)
      - Jaccard = k / |SC_i ∪ TCGA_j|
    """
    sc_ids = sorted(sc_comms.keys())
    tcga_ids = sorted(tcga_comms.keys())

    overlap_counts = pd.DataFrame(0, index=sc_ids, columns=tcga_ids)
    pvalues = pd.DataFrame(1.0, index=sc_ids, columns=tcga_ids)
    jaccards = pd.DataFrame(0.0, index=sc_ids, columns=tcga_ids)

    overlap_details = []

    for sc_id in sc_ids:
        sc_genes = sc_comms[sc_id]
        for tcga_id in tcga_ids:
            tcga_genes = tcga_comms[tcga_id]
            overlap = sc_genes & tcga_genes
            k = len(overlap)
            n = len(sc_genes)
            N = len(tcga_genes)
            M = universe_size

            overlap_counts.loc[sc_id, tcga_id] = k

            # Hypergeometric p-value: P(X >= k)
            if k > 0:
                pval = hypergeom.sf(k - 1, M, N, n)
            else:
                pval = 1.0
            pvalues.loc[sc_id, tcga_id] = pval

            # Jaccard index
            union = len(sc_genes | tcga_genes)
            jaccards.loc[sc_id, tcga_id] = k / union if union > 0 else 0.0

            if k > 0:
                overlap_details.append({
                    "sc_community": sc_id,
                    "tcga_community": tcga_id,
                    "overlap_count": k,
                    "sc_size": n,
                    "tcga_size": N,
                    "jaccard": k / union if union > 0 else 0,
                    "p_value": pval,
                    "genes": ";".join(sorted(overlap))
                })

    return overlap_counts, pvalues, jaccards, pd.DataFrame(overlap_details)


def apply_bh_to_matrix(pval_df):
    """Apply BH-FDR correction across all tests in a p-value matrix."""
    flat = pval_df.values.flatten()
    adj = bh_fdr(flat)
    adj_df = pd.DataFrame(adj.reshape(pval_df.shape),
                          index=pval_df.index, columns=pval_df.columns)
    return adj_df


# =============================================================================
# HARRIS INTERACTOR ENRICHMENT
# =============================================================================

def harris_enrichment(sc_comms, harris_symbols, universe_size):
    """
    Test enrichment of Harris A3 interactors in each SC community.

    Both SC community genes and Harris interactors are in gene symbol space,
    so no conversion is needed.
    """
    results = []

    # Intersect Harris list with genes actually in any SC community
    all_sc_genes = set()
    for genes in sc_comms.values():
        all_sc_genes |= genes
    harris_in_data = harris_symbols & all_sc_genes
    N_harris = len(harris_in_data)
    log(f"  Harris interactors present in SC network: {N_harris} / {len(harris_symbols)}")

    if N_harris == 0:
        log("  WARNING: No Harris interactors found in SC network genes")

    for comm_id in sorted(sc_comms.keys()):
        comm_genes = sc_comms[comm_id]
        overlap = comm_genes & harris_in_data
        k = len(overlap)
        n = len(comm_genes)

        if k > 0 and N_harris > 0:
            pval = hypergeom.sf(k - 1, universe_size, N_harris, n)
        else:
            pval = 1.0

        results.append({
            "community": comm_id,
            "community_size": n,
            "harris_in_community": k,
            "harris_in_network": N_harris,
            "harris_total": len(harris_symbols),
            "pct_community": 100 * k / n if n > 0 else 0,
            "pct_harris": 100 * k / N_harris if N_harris > 0 else 0,
            "p_value": pval,
            "overlap_genes": ";".join(sorted(overlap)) if overlap else ""
        })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df["fdr"] = bh_fdr(results_df["p_value"].values)
    return results_df


# =============================================================================
# FIGURE PANELS
# =============================================================================

def plot_overlap_heatmap(overlap_counts, adj_pvalues, jaccards, save_dir):
    """
    Panel 4c: Overlap heatmap colored by -log10(adjusted p-value),
    annotated with overlap count and Jaccard index.
    """
    banner("[PANEL 4c] Community Overlap Heatmap")

    # -log10 transform (cap at 10 for visualization)
    log10p = -np.log10(adj_pvalues.clip(lower=1e-10))
    log10p = log10p.clip(upper=10)

    fig, ax = plt.subplots(figsize=(max(8, len(adj_pvalues.columns) * 0.8 + 2),
                                    max(6, len(adj_pvalues.index) * 0.6 + 2)))

    # Heatmap
    im = ax.imshow(log10p.values, cmap="YlOrRd", aspect="auto",
                   vmin=0, vmax=max(3, log10p.values.max()))

    # Annotations
    for i in range(len(adj_pvalues.index)):
        for j in range(len(adj_pvalues.columns)):
            count = int(overlap_counts.iloc[i, j])
            jacc = jaccards.iloc[i, j]
            padj = adj_pvalues.iloc[i, j]

            if count > 0:
                color = "white" if log10p.iloc[i, j] > 2 else "black"
                if JACCARD_DISPLAY:
                    text = f"{count}\n({jacc:.2f})"
                else:
                    text = str(count)
                fontsize = 7 if len(adj_pvalues.columns) > 10 else 8
                ax.text(j, i, text, ha="center", va="center",
                        color=color, fontsize=fontsize,
                        fontweight="bold" if padj < 0.05 else "normal")

    # Labels
    ax.set_xticks(range(len(adj_pvalues.columns)))
    ax.set_xticklabels([f"TCGA C{c}" for c in adj_pvalues.columns],
                       rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(adj_pvalues.index)))
    ax.set_yticklabels([f"SC C{c}" for c in adj_pvalues.index], fontsize=8)

    ax.set_xlabel("TCGA Bulk Communities (Figure 2)")
    ax.set_ylabel("Single-Cell Communities (Figure 4)")
    ax.set_title("Community Overlap: SC vs TCGA Bulk\n"
                 "(color = -log10 adj.p; annotations = count / Jaccard)")

    plt.colorbar(im, ax=ax, shrink=0.6, label="-log10(BH-adjusted p)")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(save_dir, f"Panel_4c_Community_Overlap.{ext}")
        plt.savefig(path, dpi=300, bbox_inches="tight")
        log(f"  [SAVE] Panel 4c → {path}")
    plt.close()


def plot_harris_enrichment(harris_df, save_dir):
    """Panel 4d: Harris A3B interactor enrichment per community."""
    banner("[PANEL 4d] Harris A3B Interactor Enrichment")

    fig, ax = plt.subplots(figsize=(max(6, len(harris_df) * 0.6 + 2), 5))

    x = range(len(harris_df))
    bars = ax.bar(x, harris_df["harris_in_community"],
                  color=["firebrick" if p < 0.05 else "steelblue"
                         for p in harris_df["fdr"]],
                  edgecolor="black", linewidth=0.5)

    # Significance markers
    for i, row in harris_df.iterrows():
        if row["fdr"] < 0.05:
            y = row["harris_in_community"]
            ax.text(i, y + 0.3, "*", ha="center", va="bottom",
                    fontsize=12, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"C{c}" for c in harris_df["community"]], fontsize=8)
    ax.set_xlabel("Single-Cell Community")
    ax.set_ylabel("Harris A3B Interactors in Community")
    ax.set_title("A3B Interactor Enrichment per Community\n(red = FDR < 0.05)")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(save_dir, f"Panel_4d_Harris_Enrichment.{ext}")
        plt.savefig(path, dpi=300, bbox_inches="tight")
        log(f"  [SAVE] Panel 4d → {path}")
    plt.close()


def plot_community_network(G_comm, partition, save_dir):
    """Panel 4b: Network visualization of SC communities."""
    banner("[PANEL 4b] Community Network Visualization")

    if G_comm.number_of_nodes() == 0:
        log("  WARNING: Empty graph, skipping network plot")
        return

    # Get community assignments
    comm_ids = sorted(set(partition.values()))
    n_comms = len(comm_ids)

    # Color map
    cmap = plt.cm.get_cmap("tab20", max(n_comms, 2))
    comm_colors = {c: cmap(i) for i, c in enumerate(comm_ids)}
    node_colors = [comm_colors.get(partition.get(n, -1), "gray") for n in G_comm.nodes()]

    # Layout: spring with community seeding
    log("  Computing layout...")
    pos = nx.spring_layout(G_comm, k=2.0 / np.sqrt(G_comm.number_of_nodes()),
                           iterations=100, seed=42, weight="abs_weight")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Draw edges (thin, low alpha)
    nx.draw_networkx_edges(G_comm, pos, ax=ax, alpha=0.05, width=0.3,
                           edge_color="gray")

    # Draw nodes
    node_sizes = [max(10, 3 * G_comm.degree(n)) for n in G_comm.nodes()]
    nx.draw_networkx_nodes(G_comm, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.8,
                           edgecolors="black", linewidths=0.2)

    # Label A3 genes and top hubs
    labels = {}
    for n in G_comm.nodes():
        if n in A3_GENES:
            labels[n] = A3_ID_TO_ALIAS.get(n, n)
        elif G_comm.degree(n) >= 15:
            labels[n] = n  # Already a symbol in SC data

    if labels:
        nx.draw_networkx_labels(G_comm, pos, labels, ax=ax, font_size=6,
                                font_weight="bold")

    # Legend
    legend_elements = [
        Patch(facecolor=comm_colors[c], edgecolor="black", linewidth=0.5,
              label=f"C{c} (n={sum(1 for v in partition.values() if v == c)})")
        for c in comm_ids[:15]
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=7,
              ncol=2, framealpha=0.8)

    ax.set_title("Single-Cell DIFF Network Communities")
    ax.axis("off")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(save_dir, f"Panel_4b_SC_Network.{ext}")
        plt.savefig(path, dpi=300, bbox_inches="tight")
        log(f"  [SAVE] Panel 4b → {path}")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    start_time = datetime.now()
    banner("[STEP 05] Generate Figure 4 Panels + Overlap Analysis")
    log(f"Start time: {start_time}")

    ensure_dir(DIR_06_OVERLAP)
    ensure_dir(FIGURE_4_PANELS)

    # =========================================================================
    # 1. Load SC communities
    # =========================================================================
    banner("[STEP 05.1] Load single-cell communities")
    sc_comms = load_sc_communities()
    log(f"  SC communities: {len(sc_comms)}")
    for c, genes in sorted(sc_comms.items()):
        # Show a few example genes per community
        sample = sorted(genes)[:5]
        log(f"    C{c}: {len(genes)} genes (e.g. {', '.join(sample)})")

    # =========================================================================
    # 2. Load TCGA bulk communities (converted ENSG → symbols)
    # =========================================================================
    banner("[STEP 05.2] Load TCGA bulk communities (ENSG → symbol conversion)")
    tcga_comms = load_tcga_communities()

    if tcga_comms:
        for c, genes in sorted(tcga_comms.items()):
            # Show a few example genes per community
            sample = sorted(genes)[:5]
            log(f"    C{c}: {len(genes)} genes (e.g. {', '.join(sample)})")
    else:
        log("  WARNING: No TCGA communities loaded. Overlap analysis will be skipped.")

    # =========================================================================
    # 3. Compute overlap (if TCGA communities available)
    # =========================================================================
    details = pd.DataFrame()  # Initialize for later reference

    if tcga_comms:
        banner("[STEP 05.3] Hypergeometric overlap analysis")

        # Gene universe = union of all genes in either analysis (both in symbol space now)
        all_sc_genes = set()
        for genes in sc_comms.values():
            all_sc_genes |= genes

        all_tcga_genes = set()
        for genes in tcga_comms.values():
            all_tcga_genes |= genes

        universe = all_sc_genes | all_tcga_genes
        universe_size = len(universe)
        shared = all_sc_genes & all_tcga_genes

        log(f"  Gene universe: {universe_size}")
        log(f"    SC genes: {len(all_sc_genes)}")
        log(f"    TCGA genes (converted): {len(all_tcga_genes)}")
        log(f"    Shared (same symbol in both): {len(shared)}")
        log(f"    Shared examples: {sorted(shared)[:10]}")

        overlap_counts, pvalues, jaccards, details = compute_overlap_matrix(
            sc_comms, tcga_comms, universe_size
        )

        # BH-FDR correction
        adj_pvalues = apply_bh_to_matrix(pvalues)

        # Save matrices
        overlap_counts.to_csv(os.path.join(DIR_06_OVERLAP, "overlap_matrix.csv"))
        pvalues.to_csv(os.path.join(DIR_06_OVERLAP, "hypergeom_pvalues.csv"))
        adj_pvalues.to_csv(os.path.join(DIR_06_OVERLAP, "hypergeom_pvalues_BH.csv"))
        jaccards.to_csv(os.path.join(DIR_06_OVERLAP, "jaccard_matrix.csv"))
        log(f"  [SAVE] Overlap matrices → {DIR_06_OVERLAP}")

        # Save overlap gene details
        if len(details) > 0:
            details["fdr"] = bh_fdr(details["p_value"].values)
            details.to_csv(os.path.join(DIR_06_OVERLAP, "overlap_gene_details.csv"),
                           index=False)
            log(f"  [SAVE] Overlap details ({len(details)} pairs with overlap)")

        # Report significant overlaps
        sig = details[details["fdr"] < 0.05] if len(details) > 0 else pd.DataFrame()
        log(f"\n  Significant overlaps (FDR < 0.05): {len(sig)}")
        for _, row in sig.iterrows():
            log(f"    SC C{int(row['sc_community'])} × TCGA C{int(row['tcga_community'])}: "
                f"{int(row['overlap_count'])} genes (Jaccard={row['jaccard']:.3f}, "
                f"FDR={row['fdr']:.4f})")
            if row['genes']:
                gene_list = row['genes'].split(";")
                if len(gene_list) <= 10:
                    log(f"      Genes: {row['genes']}")
                else:
                    log(f"      Genes (first 10): {';'.join(gene_list[:10])}...")

        # Panel 4c: Overlap heatmap
        plot_overlap_heatmap(overlap_counts, adj_pvalues, jaccards, FIGURE_4_PANELS)

    # =========================================================================
    # 4. Harris A3B interactor enrichment
    # =========================================================================
    banner("[STEP 05.4] Harris A3B interactor enrichment")

    harris = load_harris_interactors()

    # Universe = all genes in SC communities
    all_sc_genes = set()
    for genes in sc_comms.values():
        all_sc_genes |= genes
    sc_universe = len(all_sc_genes)

    for harris_set_name, harris_symbols in harris.items():
        if not harris_symbols:
            continue

        log(f"\n  Enrichment test: {harris_set_name} ({len(harris_symbols)} genes)")
        harris_df = harris_enrichment(sc_comms, harris_symbols, sc_universe)
        harris_df.to_csv(
            os.path.join(DIR_06_OVERLAP, f"harris_enrichment_{harris_set_name}.csv"),
            index=False
        )
        log(f"  [SAVE] Harris enrichment ({harris_set_name})")

        # Report
        for _, row in harris_df.iterrows():
            if row["harris_in_community"] > 0:
                sig_marker = "***" if row["fdr"] < 0.05 else ""
                log(f"    C{int(row['community'])}: {int(row['harris_in_community'])} / "
                    f"{int(row['community_size'])} ({row['pct_community']:.1f}%) "
                    f"p={row['p_value']:.4f} FDR={row['fdr']:.4f} {sig_marker}")
                if row["overlap_genes"]:
                    log(f"      Genes: {row['overlap_genes']}")

    # Panel 4d: Harris enrichment bar chart (A3B_only)
    a3b_enrich_path = os.path.join(DIR_06_OVERLAP, "harris_enrichment_A3B_only.csv")
    if os.path.exists(a3b_enrich_path):
        a3b_df = pd.read_csv(a3b_enrich_path)
        plot_harris_enrichment(a3b_df, FIGURE_4_PANELS)

    # =========================================================================
    # 5. Panel 4b: Community network visualization
    # =========================================================================
    banner("[STEP 05.5] Community network visualization")

    graph_path = os.path.join(DIR_04_COMMUNITIES, "SC_G_comm.gpickle")
    if os.path.exists(graph_path):
        with open(graph_path, "rb") as f:
            G_comm = pickle.load(f)

        partition = {n: d.get("community", -1) for n, d in G_comm.nodes(data=True)}
        plot_community_network(G_comm, partition, FIGURE_4_PANELS)
    else:
        log("  WARNING: Community graph not found, skipping Panel 4b")

    # =========================================================================
    # 6. Summary report
    # =========================================================================
    banner("[STEP 05.6] Summary report")

    summary_lines = []
    summary_lines.append("Figure 4 — Step 05 Summary")
    summary_lines.append("=" * 60)
    summary_lines.append(f"SC communities: {len(sc_comms)}")
    summary_lines.append(f"TCGA communities: {len(tcga_comms)}")

    if tcga_comms and len(details) > 0:
        n_sig = len(details[details["fdr"] < 0.05]) if "fdr" in details.columns else 0
        summary_lines.append(f"Significant overlaps (FDR < 0.05): {n_sig}")

    summary_lines.append("")
    summary_lines.append("Panels generated:")
    summary_lines.append("  4a — UMAP cell selection (Step 00)")

    for panel, desc in [("4b", "Community network"), ("4c", "Overlap heatmap"),
                        ("4d", "Harris enrichment")]:
        files = [f for f in os.listdir(FIGURE_4_PANELS)
                 if f.startswith(f"Panel_{panel}")] if os.path.exists(FIGURE_4_PANELS) else []
        if files:
            summary_lines.append(f"  {panel} — {desc} ✓")
        else:
            summary_lines.append(f"  {panel} — {desc} (not generated)")

    summary_path = os.path.join(DIR_06_OVERLAP, "step05_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines))

    for line in summary_lines:
        print(line)

    elapsed = datetime.now() - start_time
    banner(f"[STEP 05 COMPLETE] Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
