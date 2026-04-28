#!/usr/bin/env python3
"""
Remake_C0_Network_Grant.py
===========================

Grant revision of the SC Community 0 (A3A community) network zoom.
Based on Plot_C0_Zoom_Expanded.py with the following changes:

  1. Known A3 interactor labels 5x larger (13 -> 65pt)
  2. Node colors changed so the two functional arms pop against dark background
     - Regular nodes: near-black (#1a1a1a) -- recedes, lets arms + edges dominate
     - Epithelial arm: bold orange (#f28c28) -- A3A positively correlated partners
     - Interactor arm: bold emerald (#00a878) -- A3A anti-correlated / known interactors
     - A3 genes: unchanged (red)
     - Harris interactor yellow rings: unchanged
  3. Edges thicker: 1.0 + 5.0*|w| instead of 0.3 + 3.0*|w|
  4. No TCGA shared node orange rings (not needed for C0-only grant figure)
  5. Font sizes 28-34 range for all text per grant style

Reads:
  data/FIG_4/04_communities/SC_best_partition.csv
  data/FIG_4/04_communities/SC_G_comm.gpickle
  data/FIG_4/00_input/Harris_A3_interactors.txt
  data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl

Saves:
  data/FIG_4/FIGURE_4_PANELS/Panel_C0_network_grant.pdf/.png

Usage:
  conda run -n NETWORK python Remake_C0_Network_Grant.py

Author: Jake Lehle / Claude (2026 NMF Paper, grant revision)
"""

import os, pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# =============================================================================
# PATHS  (same as Plot_C0_Zoom_Expanded.py)
# =============================================================================
BASE_DIR     = "/master/jlehle/WORKING/2026_NMF_PAPER"
SC_PARTITION = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_best_partition.csv")
SC_GRAPH     = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_G_comm.gpickle")
HARRIS_ALL   = os.path.join(BASE_DIR, "data/FIG_4/00_input/Harris_A3_interactors.txt")
CORR_DIFF    = os.path.join(BASE_DIR, "data/FIG_4/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl")
OUTPUT_DIR   = os.path.join(BASE_DIR, "data/FIG_4/FIGURE_4_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# COLORS
# =============================================================================
# A3 and highlight colors (unchanged)
COLOR_A3          = "#ed6a5a"
COLOR_HARRIS      = "#fed766"
COLOR_A3_TEXT     = "#c0392b"
COLOR_HARRIS_TEXT = "#b8860b"

# GRANT CHANGE v2: black background nodes, bold saturated arm colors
COLOR_NODE_REGULAR    = "#1a1a1a"   # near-black for all other C0 genes
COLOR_NODE_EPITHELIAL = "#f28c28"   # bold warm orange (epithelial arm, A3A+ partners)
COLOR_NODE_INTERACTOR = "#00a878"   # bold emerald green (interactor arm, A3A- partners)
COLOR_BLACK           = "#000000"

# Edge colors (unchanged)
COLOR_EDGE_POS = "firebrick"   # gained co-expression
COLOR_EDGE_NEG = "steelblue"   # lost co-expression

# =============================================================================
# FONT SIZE (grant style: 28-34 for main text)
# =============================================================================
FONT_SIZE = 28

# =============================================================================
# GENE SETS  (unchanged from Plot_C0_Zoom_Expanded.py)
# =============================================================================
A3_GENES = ["APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
            "APOBEC3F", "APOBEC3G", "APOBEC3H"]
INTERACTORS = {"HNRNPA2B1", "HSPD1", "RPL5", "TIMM8B", "RPL3"}

# Subset of interactors with confirmed A3A anti-correlation (rho_DIFF < -0.1)
# TIMM8B excluded: rho_DIFF = -0.037 (near zero), path to A3A is all positive
# edges (A3A->CYSRT1->KLK13->TIMM8B). It's in C0 via community structure,
# not via A3A anti-correlation. Still gets Harris yellow ring for identification.
INTERACTORS_ANTI_A3A = {"HNRNPA2B1", "HSPD1", "RPL5"}

# A3A's top positive and negative DIFF partners (for labeling priority)
A3A_TOP_POS = {"SULT2B1", "SPRR2A", "SPRR3", "EMP1", "CYSRT1", "RHOD", "CLDN1",
               "PRSS22", "IL1RN", "CEACAM5", "HS3ST1", "SPRR2D", "MUC4", "SCEL"}
A3A_TOP_NEG = {"SOD1", "HINT1", "RPS15", "RPS27A", "RPL4"}


def main():
    print("=" * 70)
    print("  GRANT REVISION: Community 0 Network Plot")
    print("=" * 70)

    # ---- Load partition ----
    print("\nLoading data...")
    part_df = pd.read_csv(SC_PARTITION)
    g2c = dict(zip(part_df["gene"], part_df["community"]))
    c0_genes = set(part_df[part_df["community"] == 0]["gene"])
    print(f"  C0 genes: {len(c0_genes)}")

    # ---- Load graph, extract C0 subgraph ----
    with open(SC_GRAPH, "rb") as f:
        G_full = pickle.load(f)
    c0_in_graph = [n for n in c0_genes if n in G_full.nodes()]
    G = G_full.subgraph(c0_in_graph).copy()
    print(f"  C0 subgraph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Remove degree-0 isolates
    isolates = [n for n in G.nodes() if G.degree(n) == 0]
    G.remove_nodes_from(isolates)
    print(f"  After isolate removal: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # ---- Load Harris interactors ----
    harris = set()
    if os.path.exists(HARRIS_ALL):
        with open(HARRIS_ALL) as f:
            for line in f:
                g = line.strip().split("\t")[0].strip()
                if g and not g.startswith("#"):
                    harris.add(g)
    print(f"  Harris interactors loaded: {len(harris)}")

    # ---- Load DIFF correlation matrix ----
    with open(CORR_DIFF, "rb") as f:
        corr_diff = pickle.load(f)

    # ---- Layout (same as original: expanded spring) ----
    print("  Computing expanded layout...")
    pos = nx.spring_layout(
        G, seed=42, weight="abs_weight",
        k=4.0 / np.sqrt(max(G.number_of_nodes(), 1)),
        iterations=500,
        scale=2.0
    )

    # ---- Classify nodes ----
    deg = dict(G.degree())
    md = max(deg.values()) if deg else 1
    nodes = list(G.nodes())

    a3_nodes        = [n for n in nodes if n in A3_GENES]
    harris_nodes    = [n for n in nodes if n in harris and n not in A3_GENES]
    a3a_pos_partners = [n for n in nodes if n in A3A_TOP_POS]
    a3a_neg_partners = [n for n in nodes if n in A3A_TOP_NEG or n in INTERACTORS_ANTI_A3A]
    regular_nodes   = [n for n in nodes if n not in a3_nodes and n not in harris_nodes
                       and n not in a3a_pos_partners and n not in a3a_neg_partners]

    print(f"  Node classes: A3={len(a3_nodes)}, harris={len(harris_nodes)}, "
          f"pos_partners={len(a3a_pos_partners)}, neg_partners={len(a3a_neg_partners)}, "
          f"regular={len(regular_nodes)}")

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(24, 22))

    # === EDGES ===
    # GRANT CHANGE: thicker edges (was 0.3 + 3.0*|w|, now 1.0 + 5.0*|w|)
    for u, v, d in G.edges(data=True):
        w = d.get("weight", 0)
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        width = 1.0 + 5.0 * abs(w)                      # <-- THICKER
        alpha = min(0.6, 0.20 + 0.35 * abs(w))           # slightly more opaque too
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, alpha=alpha, linewidth=width)

    # === NODE SIZES (unchanged scaling) ===
    base, scale = 200, 1200
    def ns_func(n):
        return base + scale * (deg[n] / md)

    # --- Regular nodes: NEAR-BLACK (was gray) ---
    if regular_nodes:
        sizes = [ns_func(n) for n in regular_nodes]
        nx.draw_networkx_nodes(G, pos, nodelist=regular_nodes, ax=ax,
                               node_size=sizes,
                               node_color=COLOR_NODE_REGULAR,
                               alpha=0.6, edgecolors="#444444", linewidths=0.3)

    # --- Epithelial arm (A3A positive partners): BOLD ORANGE ---
    if a3a_pos_partners:
        sizes = [ns_func(n) for n in a3a_pos_partners]
        nx.draw_networkx_nodes(G, pos, nodelist=a3a_pos_partners, ax=ax,
                               node_size=sizes,
                               node_color=COLOR_NODE_EPITHELIAL,
                               alpha=0.95, edgecolors="black", linewidths=0.8)

    # --- Interactor arm (A3A negative partners): BOLD EMERALD GREEN ---
    if a3a_neg_partners:
        sizes = [ns_func(n) for n in a3a_neg_partners]
        nx.draw_networkx_nodes(G, pos, nodelist=a3a_neg_partners, ax=ax,
                               node_size=sizes,
                               node_color=COLOR_NODE_INTERACTOR,
                               alpha=0.95, edgecolors="black", linewidths=0.8)

    # --- Harris interactor rings (yellow, unchanged) ---
    harris_in = [n for n in nodes if n in harris]
    if harris_in:
        sizes = [2.5 * ns_func(n) for n in harris_in]
        nx.draw_networkx_nodes(G, pos, nodelist=harris_in, ax=ax,
                               node_size=sizes, node_color="none",
                               edgecolors=COLOR_HARRIS, linewidths=5.0, alpha=0.9)

    # --- A3 gene nodes (red, prominent, unchanged) ---
    if a3_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=a3_nodes, ax=ax,
                               node_size=[ns_func(n) * 1.5 for n in a3_nodes],
                               node_color=COLOR_A3, alpha=1.0,
                               edgecolors="black", linewidths=2.0)
        # Red ring
        sizes = [2.5 * ns_func(n) for n in a3_nodes]
        nx.draw_networkx_nodes(G, pos, nodelist=a3_nodes, ax=ax,
                               node_size=sizes, node_color="none",
                               edgecolors=COLOR_A3, linewidths=5.0, alpha=1.0)

    # NOTE: No TCGA shared node orange rings (removed per PI directive)

    # === LABELS ===
    labels = {}
    label_colors = {}

    # A3 genes
    for n in a3_nodes:
        labels[n] = n
        label_colors[n] = COLOR_A3_TEXT

    # Harris interactors
    for n in harris_in:
        labels[n] = n
        label_colors[n] = COLOR_HARRIS_TEXT

    # A3A positive partners (epithelial arm)
    for n in a3a_pos_partners:
        labels[n] = n
        label_colors[n] = "#c46a00"  # dark orange, matches epithelial arm

    # A3A negative partners
    for n in a3a_neg_partners:
        if n not in labels:
            labels[n] = n
            label_colors[n] = "#006e50"  # dark emerald, matches interactor arm

    # Top hubs not already labeled
    top_hubs = sorted(deg.items(), key=lambda x: -x[1])[:25]
    for n, d in top_hubs:
        if n not in labels:
            labels[n] = n
            label_colors[n] = "black"

    # ---- Label positioning with repulsion (same algorithm) ----
    label_pos = {n: np.array(pos[n], dtype=float) + np.array([0, 0.04])
                 for n in labels}
    data_range = 4.0
    char_w = data_range * 0.006

    for _ in range(80):
        for i, n1 in enumerate(list(label_pos.keys())):
            for n2 in list(label_pos.keys())[i + 1:]:
                p1 = label_pos[n1]
                p2 = label_pos[n2]
                diff = p1 - p2
                dist = np.linalg.norm(diff)
                min_d = char_w * max(len(labels[n1]), len(labels[n2])) * 0.5
                if dist < min_d and dist > 1e-8:
                    push = 0.03 * (min_d - dist) * diff / dist
                    label_pos[n1] += push
                    label_pos[n2] -= push

    # ---- Draw labels ----
    for n in labels:
        col = label_colors.get(n, "black")
        nxy = pos[n]
        lxy = label_pos[n]
        offset = np.linalg.norm(lxy - nxy)

        # GRANT CHANGE: Harris interactor font is 5x larger
        if n in harris:
            fs = 65                          # <-- 5x the original 13
        elif n in A3_GENES:
            fs = FONT_SIZE                   # 28 (grant style)
        else:
            fs = 11                          # unchanged for other labels

        fw = "bold"

        if offset > char_w * 2:
            ax.annotate(labels[n], xy=nxy, xytext=lxy,
                        fontsize=fs, fontweight=fw, color=col,
                        ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15",
                                  fc="white", ec="none", alpha=0.8),
                        arrowprops=dict(arrowstyle="-", color="gray",
                                        lw=0.5, alpha=0.4))
        else:
            ax.annotate(labels[n], xy=lxy,
                        fontsize=fs, fontweight=fw, color=col,
                        ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15",
                                  fc="white", ec="none", alpha=0.8))

    # ---- Legend ----
    legend = [
        Patch(fc=COLOR_A3, ec="black", lw=2, label="APOBEC3A"),
        Patch(fc=COLOR_NODE_EPITHELIAL, ec="black", lw=0.5,
              label="A3A + co-expressed (epithelial arm)"),
        Patch(fc=COLOR_NODE_INTERACTOR, ec="black", lw=0.5,
              label="A3A anti-correlated (interactor arm)"),
        Patch(fc="none", ec=COLOR_HARRIS, lw=4,
              label="Known A3-Interactor"),
        Patch(fc=COLOR_NODE_REGULAR, ec="black", lw=0.5,
              label="Other C0 genes"),
        Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.5,
              label="Gained co-expression (red edges)"),
        Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.5,
              label="Lost co-expression (blue edges)"),
    ]
    ax.legend(handles=legend, loc="upper left", fontsize=FONT_SIZE - 8,
              framealpha=0.9,
              title="Community 0", title_fontsize=FONT_SIZE - 4)

    ax.set_title(
        f"SC Community 0 -- {G.number_of_nodes()} genes, {G.number_of_edges()} edges\n"
        f"A3A + epithelial differentiation program <-> Known A3 protein interactors",
        fontsize=FONT_SIZE, pad=15
    )
    ax.axis("off")
    plt.tight_layout()

    # ---- Save ----
    for ext in ["pdf", "png"]:
        path = os.path.join(OUTPUT_DIR, f"Panel_C0_network_grant.{ext}")
        plt.savefig(path, dpi=300, bbox_inches="tight")
        print(f"  [SAVE] {path}")
    plt.close()

    print("\nDONE")


if __name__ == "__main__":
    main()
