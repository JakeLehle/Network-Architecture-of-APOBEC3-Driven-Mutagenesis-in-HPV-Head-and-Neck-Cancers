#!/usr/bin/env python3
"""
Plot_C0_Zoom_Expanded.py
==========================

Generate an expanded, high-visibility zoom of SC Community 0 (A3A community).
Uses increased spring constant k to spread the dense core, larger figure,
and highlights the bipolar architecture (A3A epithelial arm vs interactor arm).

Usage:
  conda run -n NETWORK python Plot_C0_Zoom_Expanded.py
"""

import os, pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
SC_PARTITION = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_best_partition.csv")
SC_GRAPH     = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_G_comm.gpickle")
HARRIS_ALL   = os.path.join(BASE_DIR, "data/FIG_4/00_input/Harris_A3_interactors.txt")
CORR_DIFF    = os.path.join(BASE_DIR, "data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl")
OUTPUT_DIR   = os.path.join(BASE_DIR, "data/FIG_4/FIGURE_4_PANELS")

COLOR_A3       = "#ed6a5a"
COLOR_TCGA     = "#f18f01"
COLOR_HARRIS   = "#fed766"
COLOR_A3_TEXT  = "#c0392b"
COLOR_HARRIS_TEXT = "#b8860b"

A3_GENES = ["APOBEC3A","APOBEC3B","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G","APOBEC3H"]
INTERACTORS = {"HNRNPA2B1","HSPD1","RPL5","TIMM8B","RPL3"}

# A3A's top positive and negative DIFF partners (for labeling priority)
A3A_TOP_POS = {"SULT2B1","SPRR2A","SPRR3","EMP1","CYSRT1","RHOD","CLDN1",
               "PRSS22","IL1RN","CEACAM5","HS3ST1","SPRR2D","MUC4","SCEL"}
A3A_TOP_NEG = {"SOD1","HINT1","RPS15","RPS27A","RPL4"}

def main():
    print("Loading data...")

    # Load partition
    part_df = pd.read_csv(SC_PARTITION)
    g2c = dict(zip(part_df["gene"], part_df["community"]))
    c0_genes = set(part_df[part_df["community"]==0]["gene"])
    print(f"  C0 genes: {len(c0_genes)}")

    # Load graph and extract C0 subgraph
    with open(SC_GRAPH, "rb") as f:
        G_full = pickle.load(f)
    c0_in_graph = [n for n in c0_genes if n in G_full.nodes()]
    G = G_full.subgraph(c0_in_graph).copy()
    print(f"  C0 subgraph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Remove any remaining degree-0 nodes
    isolates = [n for n in G.nodes() if G.degree(n) == 0]
    G.remove_nodes_from(isolates)
    print(f"  After isolate removal: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Load Harris interactors
    harris = set()
    if os.path.exists(HARRIS_ALL):
        with open(HARRIS_ALL) as f:
            for line in f:
                g = line.strip().split("\t")[0].strip()
                if g and not g.startswith("#"): harris.add(g)

    # Load DIFF correlations for A3A edge annotation
    with open(CORR_DIFF, "rb") as f:
        corr_diff = pickle.load(f)

    # ---- Layout: spread the dense core ----
    # Key: increase k (ideal distance between nodes) significantly
    # Also increase iterations for better convergence
    print("  Computing expanded layout...")
    pos = nx.spring_layout(
        G, seed=42, weight="abs_weight",
        k=4.0 / np.sqrt(max(G.number_of_nodes(), 1)),  # 2x the default spread
        iterations=500,  # more iterations for convergence
        scale=2.0  # scale up coordinates
    )

    # ---- Classify nodes ----
    deg = dict(G.degree())
    md = max(deg.values()) if deg else 1
    nodes = list(G.nodes())

    a3_nodes = [n for n in nodes if n in A3_GENES]
    harris_nodes = [n for n in nodes if n in harris and n not in A3_GENES]
    a3a_pos_partners = [n for n in nodes if n in A3A_TOP_POS]
    a3a_neg_partners = [n for n in nodes if n in A3A_TOP_NEG or n in INTERACTORS]
    regular_nodes = [n for n in nodes if n not in a3_nodes and n not in harris_nodes
                     and n not in a3a_pos_partners and n not in a3a_neg_partners]

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(24, 22))

    # Edges: color by DIFF sign, scale width by magnitude
    for u, v, d in G.edges(data=True):
        w = d.get("weight", 0)
        color = "firebrick" if w > 0 else "steelblue"
        width = 0.3 + 3.0 * abs(w)
        alpha = min(0.5, 0.15 + 0.3 * abs(w))
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, alpha=alpha, linewidth=width)

    # Node sizes (scaled up)
    base, scale = 200, 1200
    ns_func = lambda n: base + scale * (deg[n] / md)

    # Regular nodes
    if regular_nodes:
        sizes = [ns_func(n) for n in regular_nodes]
        nx.draw_networkx_nodes(G, pos, nodelist=regular_nodes, ax=ax,
                               node_size=sizes, node_color="#9ecae1",  # light blue
                               alpha=0.7, edgecolors="black", linewidths=0.4)

    # A3A positive partners (epithelial arm) — warm color
    if a3a_pos_partners:
        sizes = [ns_func(n) for n in a3a_pos_partners]
        nx.draw_networkx_nodes(G, pos, nodelist=a3a_pos_partners, ax=ax,
                               node_size=sizes, node_color="#fdae6b",  # light orange
                               alpha=0.85, edgecolors="black", linewidths=0.6)

    # A3A negative partners (interactor arm) — cool color
    if a3a_neg_partners:
        sizes = [ns_func(n) for n in a3a_neg_partners]
        nx.draw_networkx_nodes(G, pos, nodelist=a3a_neg_partners, ax=ax,
                               node_size=sizes, node_color="#bcbddc",  # light purple
                               alpha=0.85, edgecolors="black", linewidths=0.6)

    # Harris interactor rings (yellow)
    harris_in = [n for n in nodes if n in harris]
    if harris_in:
        sizes = [2.5 * ns_func(n) for n in harris_in]
        nx.draw_networkx_nodes(G, pos, nodelist=harris_in, ax=ax,
                               node_size=sizes, node_color="none",
                               edgecolors=COLOR_HARRIS, linewidths=5.0, alpha=0.9)

    # A3 gene nodes (red, prominent)
    if a3_nodes:
        sizes = [2.5 * ns_func(n) for n in a3_nodes]
        nx.draw_networkx_nodes(G, pos, nodelist=a3_nodes, ax=ax,
                               node_size=[ns_func(n) * 1.5 for n in a3_nodes],
                               node_color=COLOR_A3, alpha=1.0,
                               edgecolors="black", linewidths=2.0)
        # Red ring
        nx.draw_networkx_nodes(G, pos, nodelist=a3_nodes, ax=ax,
                               node_size=sizes, node_color="none",
                               edgecolors=COLOR_A3, linewidths=5.0, alpha=1.0)

    # ---- Labels with repulsion ----
    # Select which genes to label
    labels = {}
    label_colors = {}

    # Always label: A3 genes, interactors, top A3A partners, top hubs
    for n in a3_nodes:
        labels[n] = n
        label_colors[n] = COLOR_A3_TEXT
    for n in harris_in:
        labels[n] = n
        label_colors[n] = COLOR_HARRIS_TEXT
    for n in a3a_pos_partners:
        labels[n] = n
        label_colors[n] = "#d35400"  # dark orange
    for n in a3a_neg_partners:
        if n not in labels:
            labels[n] = n
            label_colors[n] = "#7b4f9e"  # purple for neg partners

    # Top hubs not already labeled
    top_hubs = sorted(deg.items(), key=lambda x: -x[1])[:25]
    for n, d in top_hubs:
        if n not in labels:
            labels[n] = n
            label_colors[n] = "black"

    # Simple repulsion
    label_pos = {n: np.array(pos[n], dtype=float) + np.array([0, 0.04]) for n in labels}
    data_range = 4.0  # approximate from scale=2.0
    char_w = data_range * 0.006

    for _ in range(80):
        for i, n1 in enumerate(list(label_pos.keys())):
            for n2 in list(label_pos.keys())[i+1:]:
                p1 = label_pos[n1]; p2 = label_pos[n2]
                diff = p1 - p2; dist = np.linalg.norm(diff)
                min_d = char_w * max(len(labels[n1]), len(labels[n2])) * 0.5
                if dist < min_d and dist > 1e-8:
                    push = 0.03 * (min_d - dist) * diff / dist
                    label_pos[n1] += push; label_pos[n2] -= push

    for n in labels:
        col = label_colors.get(n, "black")
        nxy = pos[n]; lxy = label_pos[n]
        offset = np.linalg.norm(lxy - nxy)
        fs = 16 if n in A3_GENES else 13 if n in harris else 11
        fw = "bold"

        if offset > char_w * 2:
            ax.annotate(labels[n], xy=nxy, xytext=lxy,
                        fontsize=fs, fontweight=fw, color=col, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8),
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.4))
        else:
            ax.annotate(labels[n], xy=lxy,
                        fontsize=fs, fontweight=fw, color=col, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))

    # ---- Legend ----
    legend = [
        Patch(fc=COLOR_A3, ec="black", lw=2, label="APOBEC3A"),
        Patch(fc="#fdae6b", ec="black", lw=0.5, label="A3A + co-expressed (epithelial arm)"),
        Patch(fc="#bcbddc", ec="black", lw=0.5, label="A3A anti-correlated (interactor arm)"),
        Patch(fc="none", ec=COLOR_HARRIS, lw=4, label="Known A3-Interactor"),
        Patch(fc="#9ecae1", ec="black", lw=0.5, label="Other C0 genes"),
        Patch(fc="firebrick", ec="none", alpha=0.5, label="Gained co-expression (red edges)"),
        Patch(fc="steelblue", ec="none", alpha=0.5, label="Lost co-expression (blue edges)"),
    ]
    ax.legend(handles=legend, loc="upper left", fontsize=14, framealpha=0.9,
              title="Community 0 — Bipolar Architecture", title_fontsize=16)

    ax.set_title(f"SC Community 0 — {G.number_of_nodes()} genes, {G.number_of_edges()} edges\n"
                 f"A3A + epithelial differentiation program ↔ Known A3 protein interactors",
                 fontsize=20, pad=15)
    ax.axis("off")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(OUTPUT_DIR, f"Panel_4d_C0_zoom_expanded.{ext}")
        plt.savefig(path, dpi=300, bbox_inches="tight")
        print(f"  [SAVE] {path}")
    plt.close()

    print("\nDONE")

if __name__ == "__main__":
    main()
