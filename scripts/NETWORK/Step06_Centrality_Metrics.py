#!/usr/bin/env python3
"""
Step06_Centrality_Metrics.py

Compute per-gene network centrality metrics (degree, betweenness,
closeness, eigenvector, strength) for TOP, BOTTOM, and DIFF networks.
Identify hub genes whose network role changes between conditions.

Corresponds to original pipeline Step 16.

Input:
    05_correlation_networks/{cancer_type}/graph_objects/{ct}_G_*.gpickle
    06_communities/{cancer_type}/{ct}_best_partition.csv
    01_cleaned_expression/ensg_to_symbol.json

Output (-> data/FIG_2/07_centrality_metrics/{cancer_type}/):
    {ct}_{TOP|BOTTOM|DIFF}_metrics.csv     — per-gene centrality tables
    {ct}_hub_genes_DIFF.csv                — top hub genes by multiple metrics
    {ct}_metrics_distributions.png         — metric histograms
    {ct}_centrality_comparison.png         — TOP vs BOTTOM centrality scatter

Usage:
    python Step06_Centrality_Metrics.py
"""

import os
import json
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from network_config import (
    CANCER_TYPES, A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    DIR_01_CLEANED, DIR_04_NETWORKS, DIR_05_COMMUNITIES, DIR_06_CENTRALITY,
    banner, log, ensure_dir
)


# =============================================================================
# METRIC COMPUTATION
# =============================================================================

def compute_metrics(G):
    """Compute per-node centrality metrics. Returns DataFrame."""
    nodes = list(G.nodes())
    if len(nodes) == 0:
        return pd.DataFrame()

    deg = dict(G.degree())
    has_edges = G.number_of_edges() > 0

    btw = nx.betweenness_centrality(G) if has_edges else {n: 0.0 for n in nodes}
    close = nx.closeness_centrality(G) if has_edges else {n: 0.0 for n in nodes}

    # Eigenvector centrality (can fail on disconnected graphs)
    try:
        eigen = nx.eigenvector_centrality(G, max_iter=1000, weight="abs_weight")
    except (nx.PowerIterationFailedConvergence, nx.NetworkXError):
        eigen = {n: 0.0 for n in nodes}

    # Strength = sum of |edge weights|
    strength = {}
    for n in nodes:
        s = 0.0
        for _, d in G[n].items():
            s += abs(float(d.get("abs_weight", d.get("weight", 1.0))))
        strength[n] = s

    return pd.DataFrame({
        "node": nodes,
        "degree": [deg[n] for n in nodes],
        "betweenness": [btw[n] for n in nodes],
        "closeness": [close[n] for n in nodes],
        "eigenvector": [eigen[n] for n in nodes],
        "strength_abs": [strength[n] for n in nodes],
    })


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

    banner(f"[STEP 16] Centrality Metrics — {cancer_type}", char="=")

    cancer_dir = ensure_dir(os.path.join(DIR_06_CENTRALITY, cancer_type))

    # ---- Load graphs
    graph_dir = os.path.join(DIR_04_NETWORKS, cancer_type, "graph_objects")

    graphs = {}
    for name in ["G_top", "G_bot", "G_diff"]:
        pkl_path = os.path.join(graph_dir, f"{cancer_type}_{name}.gpickle")
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                graphs[name] = pickle.load(f)
            log(f"[STEP 16] Loaded {name}: nodes={graphs[name].number_of_nodes()}, edges={graphs[name].number_of_edges()}")
        else:
            log(f"[SKIP] {name} not found")

    if not graphs:
        log("[SKIP] No graphs loaded")
        continue

    # ---- Compute metrics for each graph
    banner("[STEP 16.1] Compute centrality metrics")

    metrics_dfs = {}
    for name, G in graphs.items():
        label = name.replace("G_", "").upper()
        m = compute_metrics(G)
        if len(m) > 0:
            m["gene_symbol"] = m["node"].map(lambda x: symbol_or_self(x, ensg_to_symbol))
            m = m.sort_values("degree", ascending=False)
        metrics_dfs[label] = m

        csv_path = os.path.join(cancer_dir, f"{cancer_type}_{label}_metrics.csv")
        m.to_csv(csv_path, index=False)
        log(f"[SAVE] {label} metrics ({len(m)} genes) -> {csv_path}")

    # ---- Metric distribution plots
    banner("[STEP 16.2] Metric distribution plots")

    for label, m in metrics_dfs.items():
        if len(m) == 0:
            continue

        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes = axes.flatten()
        metric_cols = ["degree", "betweenness", "closeness", "eigenvector", "strength_abs"]

        for i, col in enumerate(metric_cols):
            if col in m.columns:
                vals = pd.to_numeric(m[col], errors="coerce").dropna().values
                axes[i].hist(vals, bins=40, color="steelblue", edgecolor="black", linewidth=0.3)
                axes[i].set_xlabel(col)
                axes[i].set_ylabel("Gene count")
                axes[i].set_title(f"{cancer_type} | {label} | {col}")

        axes[-1].axis("off")
        plt.suptitle(f"{cancer_type} | {label} Metric Distributions", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(cancer_dir, f"{cancer_type}_{label}_metric_distributions.png"), dpi=300)
        plt.close()

    log(f"[SAVE] Distribution plots -> {cancer_dir}")

    # ---- Hub genes from DIFF network
    banner("[STEP 16.3] Identify DIFF hub genes")

    if "DIFF" in metrics_dfs and len(metrics_dfs["DIFF"]) > 0:
        m_diff = metrics_dfs["DIFF"].copy()

        # Composite hub score: rank across multiple metrics, average ranks
        for col in ["degree", "betweenness", "eigenvector", "strength_abs"]:
            m_diff[f"{col}_rank"] = m_diff[col].rank(ascending=False)

        rank_cols = [c for c in m_diff.columns if c.endswith("_rank")]
        m_diff["hub_score_rank"] = m_diff[rank_cols].mean(axis=1)
        m_diff = m_diff.sort_values("hub_score_rank")

        # Top 50 hub genes
        top_hubs = m_diff.head(50)
        hub_path = os.path.join(cancer_dir, f"{cancer_type}_hub_genes_DIFF.csv")
        top_hubs.to_csv(hub_path, index=False)
        log(f"[SAVE] Top 50 DIFF hubs -> {hub_path}")

        log("\n[STEP 16.3] Top 20 DIFF hub genes:")
        for _, row in top_hubs.head(20).iterrows():
            sym = row["gene_symbol"]
            deg = int(row["degree"])
            btw = row["betweenness"]
            eig = row["eigenvector"]
            is_a3 = " *** A3 ***" if row["node"] in A3_GENES else ""
            is_bm = " [biomarker]" if row["node"] in BIOMARKERS else ""
            log(f"  {sym:12s}  deg={deg:4d}  btw={btw:.4f}  eig={eig:.4f}{is_a3}{is_bm}")

    # ---- Centrality comparison: TOP vs BOTTOM
    banner("[STEP 16.4] Centrality comparison (TOP vs BOTTOM)")

    if "TOP" in metrics_dfs and "BOT" in metrics_dfs:
        m_top = metrics_dfs["TOP"].set_index("node")
        m_bot = metrics_dfs["BOT"].set_index("node")

        shared = sorted(set(m_top.index) & set(m_bot.index))

        if len(shared) > 10:
            fig, axes = plt.subplots(1, 3, figsize=(18, 5))

            for ax, metric in zip(axes, ["degree", "betweenness", "eigenvector"]):
                x = m_bot.loc[shared, metric].values.astype(float)
                y = m_top.loc[shared, metric].values.astype(float)

                ax.scatter(x, y, s=10, alpha=0.4, c="steelblue")

                # Highlight A3 genes
                for g in A3_GENES:
                    if g in shared:
                        ax.scatter(m_bot.loc[g, metric], m_top.loc[g, metric],
                                   s=80, c="gold", edgecolors="black", linewidths=1.0, zorder=5)
                        ax.annotate(symbol_or_self(g, ensg_to_symbol),
                                    (m_bot.loc[g, metric], m_top.loc[g, metric]),
                                    fontsize=7, fontweight="bold")

                # Diagonal
                lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
                ax.plot([0, lim], [0, lim], ls="--", c="gray", lw=0.5)

                ax.set_xlabel(f"BOTTOM {metric}")
                ax.set_ylabel(f"TOP {metric}")
                ax.set_title(f"{cancer_type} | {metric}")

            plt.suptitle(f"{cancer_type} | Centrality: TOP vs BOTTOM", fontsize=14)
            plt.tight_layout()
            comp_path = os.path.join(cancer_dir, f"{cancer_type}_centrality_comparison.png")
            plt.savefig(comp_path, dpi=300)
            plt.close()
            log(f"[SAVE] Centrality comparison -> {comp_path}")

    # ---- Annotate hub genes with community membership (if available)
    banner("[STEP 16.5] Annotate hubs with community membership")

    part_path = os.path.join(DIR_05_COMMUNITIES, cancer_type, f"{cancer_type}_best_partition.csv")
    if os.path.exists(part_path) and "DIFF" in metrics_dfs:
        partition = pd.read_csv(part_path)
        gene_to_comm = dict(zip(partition["gene"], partition["community"]))

        m_diff = metrics_dfs["DIFF"].copy()
        m_diff["community"] = m_diff["node"].map(gene_to_comm)

        annotated_path = os.path.join(cancer_dir, f"{cancer_type}_DIFF_metrics_with_communities.csv")
        m_diff.sort_values("hub_score_rank" if "hub_score_rank" in m_diff.columns else "degree",
                          ascending=True if "hub_score_rank" in m_diff.columns else False
                          ).to_csv(annotated_path, index=False)
        log(f"[SAVE] Annotated metrics -> {annotated_path}")

    log(f"\n[STEP 16 COMPLETE for {cancer_type}]")

banner("[STEP 07 COMPLETE — ALL CANCER TYPES]")
