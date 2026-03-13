#!/usr/bin/env python3
"""
Step05_Correlation_Networks.py

Build Spearman co-expression networks for TOP and BOTTOM groups
separately, compute the DIFF network (TOP - BOTTOM), and apply
thresholds to create weighted graphs.

Corresponds to original pipeline Step 14.

Input:
    04_TOP_BOTTOM_groups/{cancer_type}/{ct}_TOP_group.parquet
    04_TOP_BOTTOM_groups/{cancer_type}/{ct}_BOTTOM_group.parquet
    04_TOP_BOTTOM_groups/{cancer_type}/{ct}_selected_genes_filtered.csv
    01_cleaned_expression/ensg_to_symbol.json

Output (-> data/FIG_2/05_correlation_networks/{cancer_type}/):
    corr_matrices/     — TOP, BOTTOM, DIFF correlation matrices (.parquet + .csv)
    heatmaps/          — clustered heatmaps + correlation distribution plots
    edge_lists/        — network edge lists (for Cytoscape/Gephi)
    network_plots/     — TOP, BOTTOM, DIFF network visualizations

Usage:
    python Step05_Correlation_Networks.py
"""

import os
import json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from network_config import (
    CANCER_TYPES, CLINICAL_COLS, A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    CORRELATION_METHOD, CORR_THRESHOLD, DIFF_THRESHOLD,
    DIR_01_CLEANED, DIR_04_GROUPS, DIR_05_NETWORKS,
    banner, log, ensure_dir
)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def build_weighted_graph(corr_df, threshold):
    """Build undirected weighted graph from correlation matrix with |corr| >= threshold."""
    corr_df = corr_df.fillna(0)
    genes = list(corr_df.index)
    G = nx.Graph()
    G.add_nodes_from(genes)

    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            w = float(corr_df.iat[i, j])
            if abs(w) >= threshold and abs(w) < 0.999999:
                G.add_edge(genes[i], genes[j], weight=w, abs_weight=abs(w))
    return G


def save_graph_csv(G, out_dir, prefix):
    """Save edge list and node metrics to CSV."""
    ensure_dir(out_dir)

    # Edges
    edges = []
    for u, v, d in G.edges(data=True):
        edges.append({
            "source": u, "target": v,
            "weight": d.get("weight", 1.0),
            "abs_weight": d.get("abs_weight", abs(d.get("weight", 1.0)))
        })
    pd.DataFrame(edges).to_csv(os.path.join(out_dir, f"{prefix}_edges.csv"), index=False)

    # Nodes
    deg = dict(G.degree())
    btw = nx.betweenness_centrality(G) if G.number_of_edges() > 0 else {n: 0.0 for n in G.nodes()}
    nodes_df = pd.DataFrame({
        "node": list(G.nodes()),
        "degree": [deg[n] for n in G.nodes()],
        "betweenness": [btw[n] for n in G.nodes()],
    })
    nodes_df.to_csv(os.path.join(out_dir, f"{prefix}_nodes.csv"), index=False)


def remove_isolates_in_both(G_top, G_bot, G_union):
    """Remove genes isolated in BOTH top and bottom networks."""
    H = G_union.copy()
    for n in list(H.nodes()):
        if G_top.degree(n) == 0 and G_bot.degree(n) == 0:
            H.remove_node(n)
    return H


def get_highlight_nodes(G):
    """Get A3 genes + biomarkers present in the graph."""
    nodes = set(G.nodes())
    highlights = [g for g in A3_GENES if g in nodes]
    highlights += [g for g in BIOMARKERS if g in nodes]
    return list(dict.fromkeys(highlights))


def symbol_or_self(ensg, mapping):
    """Return gene symbol if available, otherwise ENSG ID."""
    return mapping.get(str(ensg), str(ensg))


# =============================================================================
# LOAD symbol mapping
# =============================================================================
symbol_path = os.path.join(DIR_01_CLEANED, "ensg_to_symbol.json")
with open(symbol_path) as f:
    ensg_to_symbol = json.load(f)
log(f"[INIT] Loaded ENSG->symbol mapping: {len(ensg_to_symbol)} entries")


# =============================================================================
# LOOP over cancer types
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"[STEP 14] Correlation Networks — {cancer_type}", char="=")

    # ---- Setup output directories
    base_dir = ensure_dir(os.path.join(DIR_05_NETWORKS, cancer_type))
    corr_dir = ensure_dir(os.path.join(base_dir, "corr_matrices"))
    heat_dir = ensure_dir(os.path.join(base_dir, "heatmaps"))
    edge_dir = ensure_dir(os.path.join(base_dir, "edge_lists"))
    plot_dir = ensure_dir(os.path.join(base_dir, "network_plots"))

    # ---- Load groups
    top_path = os.path.join(DIR_04_GROUPS, cancer_type, f"{cancer_type}_TOP_group.parquet")
    bot_path = os.path.join(DIR_04_GROUPS, cancer_type, f"{cancer_type}_BOTTOM_group.parquet")
    genes_path = os.path.join(DIR_04_GROUPS, cancer_type, f"{cancer_type}_selected_genes_filtered.csv")

    for p in [top_path, bot_path, genes_path]:
        if not os.path.exists(p):
            log(f"[SKIP] Missing file: {p}")
            continue

    top_df = pd.read_parquet(top_path)
    bot_df = pd.read_parquet(bot_path)
    selected_genes = pd.read_csv(genes_path)["gene"].tolist()

    log(f"[STEP 14.0] TOP samples: {len(top_df)}")
    log(f"[STEP 14.0] BOTTOM samples: {len(bot_df)}")
    log(f"[STEP 14.0] Selected genes: {len(selected_genes)}")

    if len(selected_genes) < 10:
        log("[SKIP] Too few genes (<10)")
        continue

    # ---- Extract numeric expression matrices (samples × genes)
    top_num = top_df[selected_genes].apply(pd.to_numeric, errors="coerce")
    bot_num = bot_df[selected_genes].apply(pd.to_numeric, errors="coerce")

    # ---- Compute Spearman correlations
    banner("[STEP 14.1] Compute Spearman correlations")

    corr_top = top_num.corr(method=CORRELATION_METHOD)
    corr_bot = bot_num.corr(method=CORRELATION_METHOD)
    corr_diff = corr_top - corr_bot

    log(f"[STEP 14.1] Correlation matrix shape: {corr_top.shape}")

    # Save correlation matrices
    corr_top.to_parquet(os.path.join(corr_dir, f"{cancer_type}_corr_TOP.parquet"))
    corr_bot.to_parquet(os.path.join(corr_dir, f"{cancer_type}_corr_BOTTOM.parquet"))
    corr_diff.to_parquet(os.path.join(corr_dir, f"{cancer_type}_corr_DIFF.parquet"))
    log(f"[SAVE] Correlation matrices -> {corr_dir}")

    # ---- Correlation distribution plots
    banner("[STEP 14.2] Correlation distribution plots")

    def upper_tri_vals(mat):
        a = mat.values
        iu = np.triu_indices_from(a, k=1)
        return a[iu]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, vals, label, color in [
        (axes[0], upper_tri_vals(corr_top), "TOP", "firebrick"),
        (axes[1], upper_tri_vals(corr_bot), "BOTTOM", "steelblue"),
        (axes[2], upper_tri_vals(corr_diff), "DIFF (TOP-BOT)", "purple"),
    ]:
        ax.hist(vals, bins=80, color=color, alpha=0.7, edgecolor="black", linewidth=0.3)
        ax.set_xlabel("Spearman rho")
        ax.set_ylabel("Gene pairs")
        ax.set_title(f"{cancer_type} | {label}")
        ax.axvline(0, ls="--", c="black", lw=0.5)
    plt.tight_layout()
    dist_path = os.path.join(heat_dir, f"{cancer_type}_corr_distributions.png")
    plt.savefig(dist_path, dpi=300)
    plt.close()
    log(f"[SAVE] Distribution plot -> {dist_path}")

    # Violin plot
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.violinplot(
        [upper_tri_vals(corr_top), upper_tri_vals(corr_bot), upper_tri_vals(corr_diff)],
        showmeans=False, showmedians=True, showextrema=True
    )
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["TOP", "BOTTOM", "DIFF"])
    ax.set_ylabel("Spearman correlation")
    ax.set_title(f"{cancer_type} | Correlation distributions")
    plt.tight_layout()
    violin_path = os.path.join(heat_dir, f"{cancer_type}_corr_violin.png")
    plt.savefig(violin_path, dpi=300)
    plt.close()
    log(f"[SAVE] Violin -> {violin_path}")

    # ---- Build graphs
    banner("[STEP 14.3] Build weighted graphs")

    G_top  = build_weighted_graph(corr_top, CORR_THRESHOLD)
    G_bot  = build_weighted_graph(corr_bot, CORR_THRESHOLD)
    G_diff = build_weighted_graph(corr_diff, DIFF_THRESHOLD)

    log(f"[STEP 14.3] G_top:  nodes={G_top.number_of_nodes()}, edges={G_top.number_of_edges()}")
    log(f"[STEP 14.3] G_bot:  nodes={G_bot.number_of_nodes()}, edges={G_bot.number_of_edges()}")
    log(f"[STEP 14.3] G_diff: nodes={G_diff.number_of_nodes()}, edges={G_diff.number_of_edges()}")

    # Ensure abs_weight on DIFF edges
    for u, v, d in G_diff.edges(data=True):
        d["abs_weight"] = abs(d.get("weight", 0.0))

    # Remove isolates from DIFF
    G_diff_noiso = G_diff.copy()
    G_diff_noiso.remove_nodes_from([n for n in G_diff_noiso.nodes() if G_diff_noiso.degree(n) == 0])
    log(f"[STEP 14.3] G_diff (no isolates): nodes={G_diff_noiso.number_of_nodes()}, edges={G_diff_noiso.number_of_edges()}")

    # ---- Save edge lists
    save_graph_csv(G_top, edge_dir, f"{cancer_type}_TOP")
    save_graph_csv(G_bot, edge_dir, f"{cancer_type}_BOTTOM")
    save_graph_csv(G_diff, edge_dir, f"{cancer_type}_DIFF")
    log(f"[SAVE] Edge lists -> {edge_dir}")

    # ---- Network visualization
    banner("[STEP 14.4] Network plots")

    # Compute shared layout from the DIFF network
    if G_diff_noiso.number_of_nodes() > 0:
        pos = nx.spring_layout(G_diff_noiso, seed=42, weight="abs_weight",
                               k=2.0/np.sqrt(max(G_diff_noiso.number_of_nodes(), 1)),
                               iterations=200)
    else:
        pos = {}

    for G, label, color_pos, color_neg in [
        (G_top, "TOP", "firebrick", "steelblue"),
        (G_bot, "BOTTOM", "steelblue", "firebrick"),
        (G_diff_noiso, "DIFF", "darkorange", "purple"),
    ]:
        if G.number_of_nodes() == 0:
            log(f"[SKIP] {label} graph is empty")
            continue

        # Use shared pos where available; spring_layout for missing nodes
        local_pos = {n: pos[n] for n in G.nodes() if n in pos}
        missing = [n for n in G.nodes() if n not in local_pos]
        if missing:
            extra = nx.spring_layout(nx.subgraph(G, missing), seed=42)
            local_pos.update(extra)

        fig, ax = plt.subplots(figsize=(14, 12))

        # Draw edges
        edges = list(G.edges(data=True))
        if edges:
            weights = [abs(d.get("weight", 0.5)) for _, _, d in edges]
            max_w = max(weights) if weights else 1
            edge_widths = [0.3 + 2 * w / max_w for w in weights]
            edge_colors = [color_pos if d.get("weight", 0) > 0 else color_neg for _, _, d in edges]
            nx.draw_networkx_edges(G, local_pos, ax=ax,
                                   width=edge_widths, edge_color=edge_colors, alpha=0.3)

        # Draw nodes
        highlights = get_highlight_nodes(G)
        regular = [n for n in G.nodes() if n not in highlights]

        nx.draw_networkx_nodes(G, local_pos, nodelist=regular, ax=ax,
                               node_size=30, node_color="lightgray", alpha=0.6)
        if highlights:
            nx.draw_networkx_nodes(G, local_pos, nodelist=highlights, ax=ax,
                                   node_size=120, node_color="gold", edgecolors="black", linewidths=1.0)

        # Label highlighted nodes with gene symbols
        labels = {n: symbol_or_self(n, ensg_to_symbol) for n in highlights}
        nx.draw_networkx_labels(G, local_pos, labels=labels, ax=ax, font_size=7, font_weight="bold")

        ax.set_title(f"{cancer_type} | {label} network | "
                     f"nodes={G.number_of_nodes()}, edges={G.number_of_edges()}", fontsize=13)
        ax.axis("off")
        plt.tight_layout()

        net_path = os.path.join(plot_dir, f"{cancer_type}_{label}_network.png")
        plt.savefig(net_path, dpi=300)
        plt.close()
        log(f"[SAVE] {label} network -> {net_path}")

    # ---- Save graphs as pickles for downstream steps
    import pickle
    pickle_dir = ensure_dir(os.path.join(base_dir, "graph_objects"))
    for G_obj, name in [(G_top, "G_top"), (G_bot, "G_bot"), (G_diff, "G_diff"), (G_diff_noiso, "G_diff_noiso")]:
        pkl_path = os.path.join(pickle_dir, f"{cancer_type}_{name}.gpickle")
        with open(pkl_path, "wb") as f:
            pickle.dump(G_obj, f)
    log(f"[SAVE] Graph pickles -> {pickle_dir}")

    log(f"\n[STEP 14 COMPLETE for {cancer_type}]")
    print(f"  G_top:  {G_top.number_of_nodes()} nodes, {G_top.number_of_edges()} edges")
    print(f"  G_bot:  {G_bot.number_of_nodes()} nodes, {G_bot.number_of_edges()} edges")
    print(f"  G_diff: {G_diff_noiso.number_of_nodes()} nodes, {G_diff_noiso.number_of_edges()} edges")

banner("[STEP 05 COMPLETE — ALL CANCER TYPES]")
