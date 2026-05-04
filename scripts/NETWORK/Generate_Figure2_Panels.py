#!/usr/bin/env python3
"""
Generate_Figure2_Panels.py
===========================

Generate publication-quality Figure 2 panels from pipeline outputs.
Reads from Steps 03-06 outputs + node importance scores and produces:

  Panel A  -- SBS2 vs A3 sum selection plot (s=100 for consistency)
  Panel B  -- Full community network (intra_score-based node/label sizing)
  Panel C  -- Community 4 (A3B) zoomed network
  Supplement -- Methods graphical abstract: HIGH -> LOW -> DIFF -> network
  Community zooms -- All communities, intra_score-based sizing

Style: 28-34pt text, all hex color codes, PDF + PNG (300 DPI).

Run AFTER the full pipeline (Steps 01-06) has completed.

Usage:
    conda run -n NETWORK python Generate_Figure2_Panels.py
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
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon as MplPolygon, Patch
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from network_config import (
    CANCER_TYPES, A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    COMMUNITY_BASE_SEED,
    DIR_01_CLEANED, DIR_02_MERGED, DIR_03_DIFFEXPR,
    DIR_04_NETWORKS, DIR_05_COMMUNITIES, DIR_06_CENTRALITY,
    FIG2_ROOT, banner, log, ensure_dir
)


# =============================================================================
# LOAD SHARED DATA
# =============================================================================
banner("[INIT] Load shared data")

# Symbol mapping
with open(os.path.join(DIR_01_CLEANED, "ensg_to_symbol.json")) as f:
    ensg_to_symbol = json.load(f)

def sym(ensg):
    return ensg_to_symbol.get(str(ensg), str(ensg))


# Figure 1 color palette
COLOR_TEAL  = "#9bc1bc"
COLOR_CORAL = "#ed6a5a"
COLOR_CREAM = "#f4f1bb"
COLOR_HIGH  = "#ed6a5a"   # coral for HIGH SBS2
COLOR_LOW   = "#5b8e7d"   # muted teal for LOW SBS2
COLOR_GRAY  = "#d0d0d0"

# Network edge colors (all hex)
COLOR_EDGE_POS = "#b22222"   # firebrick - gained co-expression
COLOR_EDGE_NEG = "#4682b4"   # steelblue - lost co-expression

# A3 gene highlight colors
COLOR_A3      = "#ed6a5a"
COLOR_A3_TEXT = "#c0392b"

# Font sizes (28-34 range for publication)
FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22

# Label background alpha (halved so nodes show through)
LABEL_BG_ALPHA = 0.42

# Output directory
FIG_DIR = ensure_dir(os.path.join(FIG2_ROOT, "FIGURE_2_PANELS"))


# =============================================================================
# HELPER: Label repulsion algorithm (ported from grant script)
# =============================================================================
def repel_labels(pos, labels, data_range=12.0, n_iters=120, push_strength=0.04):
    """
    Iteratively push overlapping labels apart.
    Returns dict of adjusted label positions.
    Labels shifted up from nodes so underlying structure stays visible.
    """
    label_pos = {n: np.array(pos[n], dtype=float) + np.array([0, 0.16])
                 for n in labels}
    char_w = data_range * 0.012

    for _ in range(n_iters):
        keys = list(label_pos.keys())
        for i, n1 in enumerate(keys):
            for n2 in keys[i + 1:]:
                p1 = label_pos[n1]
                p2 = label_pos[n2]
                diff = p1 - p2
                dist = np.linalg.norm(diff)
                min_d = char_w * max(len(str(labels[n1])),
                                     len(str(labels[n2]))) * 0.55
                if dist < min_d and dist > 1e-8:
                    push = push_strength * (min_d - dist) * diff / dist
                    label_pos[n1] = label_pos[n1] + push
                    label_pos[n2] = label_pos[n2] - push

    return label_pos


def draw_labels_with_boxes(ax, pos, label_pos, labels, font_sizes,
                           font_colors, char_w=0.06):
    """
    Draw labels with semi-transparent white background boxes and leader lines.
    """
    for n in labels:
        col = font_colors.get(n, "#000000")
        fs  = font_sizes.get(n, 12)
        nxy = np.array(pos[n])
        lxy = np.array(label_pos[n])
        offset = np.linalg.norm(lxy - nxy)

        kwargs = dict(
            fontsize=fs, fontweight="bold", color=col,
            ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.15",
                      fc="white", ec="none", alpha=LABEL_BG_ALPHA)
        )

        if offset > char_w * 2.5:
            ax.annotate(labels[n], xy=nxy, xytext=lxy, **kwargs,
                        arrowprops=dict(arrowstyle="-", color="#888888",
                                        lw=0.6, alpha=0.5))
        else:
            ax.annotate(labels[n], xy=lxy, **kwargs)


# =============================================================================
# HELPER: Build community-ordered gene list with hierarchical clustering
# =============================================================================
def build_community_ordered_genes(partition_df, ref_matrix, gene_to_comm):
    """
    Order genes by community with hierarchical clustering within and between.
    Returns (ordered_genes, comm_boundaries, comm_ids_sorted).
    """
    comm_genes = [g for g in partition_df["gene"] if g in ref_matrix.index]

    # Group by community
    comm_to_genes = {}
    for g in comm_genes:
        c = gene_to_comm.get(g, -1)
        comm_to_genes.setdefault(c, []).append(g)

    # Cluster within each community
    clustered = {}
    for c, genes in comm_to_genes.items():
        if len(genes) <= 2:
            clustered[c] = genes
            continue
        sub = ref_matrix.loc[genes, genes].values
        sub = np.nan_to_num(sub, nan=0.0)
        try:
            dist = 1.0 - np.abs(sub)
            np.fill_diagonal(dist, 0)
            dist = (dist + dist.T) / 2.0
            dist = np.clip(dist, 0, 2)
            Z = linkage(squareform(dist, checks=False), method="average")
            order = leaves_list(Z)
            clustered[c] = [genes[i] for i in order]
        except Exception:
            clustered[c] = genes

    # Cluster between communities
    comm_ids = list(clustered.keys())
    if len(comm_ids) > 2:
        n_c = len(comm_ids)
        cd = np.zeros((n_c, n_c))
        for i in range(n_c):
            for j in range(i + 1, n_c):
                block = ref_matrix.loc[clustered[comm_ids[i]],
                                       clustered[comm_ids[j]]].values
                d = 1.0 - np.nanmean(np.abs(block))
                cd[i, j] = cd[j, i] = d
        Z_c = linkage(squareform(cd, checks=False), method="average")
        comm_ids_sorted = [comm_ids[i] for i in leaves_list(Z_c)]
    else:
        comm_ids_sorted = sorted(comm_ids,
                                 key=lambda c: len(clustered[c]), reverse=True)

    # Build final ordered list
    ordered_genes = []
    comm_boundaries = []
    for c in comm_ids_sorted:
        start = len(ordered_genes)
        ordered_genes.extend(clustered[c])
        comm_boundaries.append((start, len(ordered_genes), c))

    return ordered_genes, comm_boundaries, comm_ids_sorted


# =============================================================================
# LOOP CANCER TYPES
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"FIGURE 2 -- {cancer_type}", char="=")

    ct_fig_dir = ensure_dir(os.path.join(FIG_DIR, cancer_type))

    # ---- Load merged data (for all tumors)
    merged = pd.read_pickle(os.path.join(DIR_02_MERGED, "TCGA_merged_expression_SBS.pkl"))
    all_hnsc = merged[merged["Project_ID"] == cancer_type].copy()
    for c in ["SBS2", "A3A", "A3B"]:
        all_hnsc[c] = pd.to_numeric(all_hnsc[c], errors="coerce")
    all_hnsc = all_hnsc.dropna(subset=["SBS2", "A3A", "A3B"])
    all_hnsc["A3_sum"] = all_hnsc["A3A"] + all_hnsc["A3B"]
    log(f"All HNSC tumors: {len(all_hnsc)}")

    # ---- Load groups from Step03
    high_path = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                             f"{cancer_type}_SBS2_HIGH_group.pkl")
    low_path  = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                             f"{cancer_type}_SBS2_LOW_group.pkl")
    if not os.path.exists(high_path):
        log(f"[SKIP] Groups not found")
        continue
    group_high = pd.read_pickle(high_path)
    group_low  = pd.read_pickle(low_path)
    log(f"HIGH: {len(group_high)}, LOW: {len(group_low)}")

    # ---- Load correlation matrices (TOP, BOTTOM, DIFF)
    corr_dir = os.path.join(DIR_04_NETWORKS, cancer_type, "corr_matrices")

    corr_diff_path = os.path.join(corr_dir, f"{cancer_type}_corr_DIFF.pkl")
    corr_top_path  = os.path.join(corr_dir, f"{cancer_type}_corr_TOP.pkl")
    corr_bot_path  = os.path.join(corr_dir, f"{cancer_type}_corr_BOTTOM.pkl")

    if not os.path.exists(corr_diff_path):
        log("[SKIP] DIFF correlation matrix not found")
        continue
    corr_diff = pd.read_pickle(corr_diff_path)
    log(f"DIFF corr matrix: {corr_diff.shape}")

    corr_top = pd.read_pickle(corr_top_path) if os.path.exists(corr_top_path) else None
    corr_bot = pd.read_pickle(corr_bot_path) if os.path.exists(corr_bot_path) else None
    if corr_top is not None:
        log(f"TOP corr matrix: {corr_top.shape}")
    if corr_bot is not None:
        log(f"BOTTOM corr matrix: {corr_bot.shape}")

    # ---- Load community assignments
    part_path = os.path.join(DIR_05_COMMUNITIES, cancer_type,
                             f"{cancer_type}_best_partition.csv")
    if not os.path.exists(part_path):
        log("[SKIP] Best partition not found")
        continue
    partition_df = pd.read_csv(part_path)
    gene_to_comm = dict(zip(partition_df["gene"], partition_df["community"]))
    log(f"Partition: {len(gene_to_comm)} genes, "
        f"{partition_df['community'].nunique()} communities")

    # ---- Load community graph
    comm_graph_path = os.path.join(DIR_05_COMMUNITIES, cancer_type,
                                   f"{cancer_type}_G_comm.gpickle")
    G_comm = None
    if os.path.exists(comm_graph_path):
        with open(comm_graph_path, "rb") as f:
            G_comm = pickle.load(f)
        log(f"Community graph: {G_comm.number_of_nodes()} nodes, "
            f"{G_comm.number_of_edges()} edges")

    # ---- Load node importance scores
    sizing_path = os.path.join(FIG2_ROOT, "DIAGNOSTIC_AUDIT",
                               "Figure2_Node_Sizing.tsv")
    node_intra = {}
    node_inter = {}
    if os.path.exists(sizing_path):
        sizing_df = pd.read_csv(sizing_path, sep="\t")
        for _, row in sizing_df.iterrows():
            node_intra[row["ensg_id"]]      = float(row["intra_score"])
            node_intra[row["gene_symbol"]]   = float(row["intra_score"])
            node_inter[row["ensg_id"]]       = float(row["inter_score"])
            node_inter[row["gene_symbol"]]   = float(row["inter_score"])
        log(f"Node importance scores loaded: {len(sizing_df)} genes")
    else:
        log("[WARNING] Figure2_Node_Sizing.tsv not found, falling back to degree")

    # Community color map (used across panels)
    if G_comm is not None:
        comm_to_nodelist = {}
        for n in G_comm.nodes():
            c = gene_to_comm.get(n, -1)
            comm_to_nodelist.setdefault(c, []).append(n)
        unique_comms = sorted(comm_to_nodelist.keys())
        n_comms = len(unique_comms)
        cmap_comms = plt.colormaps["tab20"]
        comm_color_map = {c: mcolors.to_hex(cmap_comms(i))
                          for i, c in enumerate(unique_comms)}

    # =====================================================================
    # PANEL A -- SBS2 vs A3 Sum Selection Plot
    # =====================================================================
    banner("[PANEL A] Selection plot")

    a3_median = float(all_hnsc["A3_sum"].median())
    sbs2_low_max  = float(group_low["SBS2"].max())
    sbs2_high_min = float(group_high["SBS2"].min())
    x_max_data = all_hnsc["A3_sum"].max()
    y_max_data = all_hnsc["SBS2"].max()

    # ---- Data-driven x-axis break ----
    a3_sorted = all_hnsc["A3_sum"].sort_values(ascending=False).values
    log(f"[PANEL A] Top 5 A3_sum values: {a3_sorted[:5].round(2).tolist()}")

    has_outliers = False
    if len(a3_sorted) >= 3:
        best_gap_ratio = 1.0
        best_gap_idx = -1
        n_check = min(10, len(a3_sorted) - 1)
        for i in range(n_check):
            ratio = a3_sorted[i] / a3_sorted[i + 1] if a3_sorted[i + 1] > 0 else 1.0
            if ratio > best_gap_ratio:
                best_gap_ratio = ratio
                best_gap_idx = i

        log(f"[PANEL A] Largest gap: position {best_gap_idx} -> {best_gap_idx+1}, "
            f"values {a3_sorted[best_gap_idx]:.2f} vs {a3_sorted[best_gap_idx+1]:.2f}, "
            f"ratio {best_gap_ratio:.2f}")

        if best_gap_ratio > 1.5 and best_gap_idx >= 0:
            outlier_min   = a3_sorted[best_gap_idx]
            cluster_max   = a3_sorted[best_gap_idx + 1]
            x_break_start = cluster_max * 1.10
            x_break_end   = outlier_min * 0.90
            has_outliers   = True
            n_outliers     = best_gap_idx + 1
            log(f"[PANEL A] {n_outliers} outlier(s), break: "
                f"{x_break_start:.1f} -- {x_break_end:.1f}")

    if has_outliers:
        fig, (ax_main, ax_break) = plt.subplots(
            1, 2, sharey=True, figsize=(16, 10),
            gridspec_kw={"width_ratios": [3, 1], "wspace": 0.03}
        )
        axes_list = [ax_main, ax_break]
    else:
        fig, ax_main = plt.subplots(figsize=(14, 10))
        ax_break = None
        axes_list = [ax_main]

    for ax in axes_list:
        y_pad = y_max_data * 0.05
        y_max = y_max_data + y_pad

        if ax == ax_main:
            x_lo = -2
            x_hi = x_break_start if has_outliers else x_max_data * 1.05
        elif ax == ax_break:
            x_lo, x_hi = x_break_end, x_max_data * 1.05

        # Background regions (right of A3 median)
        cream_poly = MplPolygon(
            [(a3_median, 0), (a3_median, sbs2_low_max),
             (x_hi, sbs2_low_max), (x_hi, 0)],
            closed=True, facecolor=COLOR_CREAM, alpha=0.75,
            edgecolor="none", zorder=0)
        ax.add_patch(cream_poly)

        coral_poly = MplPolygon(
            [(a3_median, sbs2_high_min), (a3_median, y_max),
             (x_hi, y_max), (x_hi, sbs2_high_min)],
            closed=True, facecolor=COLOR_CORAL, alpha=0.55,
            edgecolor="none", zorder=0)
        ax.add_patch(coral_poly)

        if sbs2_high_min > sbs2_low_max:
            mid_poly = MplPolygon(
                [(a3_median, sbs2_low_max), (a3_median, sbs2_high_min),
                 (x_hi, sbs2_high_min), (x_hi, sbs2_low_max)],
                closed=True, facecolor="#e8e4c8", alpha=0.55,
                edgecolor="none", zorder=0)
            ax.add_patch(mid_poly)

        # Boundary lines
        ax.axvline(a3_median, ls="--", c="#4D4D4D", lw=0.8, zorder=1)
        ax.axhline(sbs2_high_min, ls=":", c="#4D4D4D", lw=0.8, alpha=0.6, zorder=1)
        ax.axhline(sbs2_low_max,  ls=":", c="#4D4D4D", lw=0.8, alpha=0.6, zorder=1)

        # Data points -- ALL s=100 for GitHub consistency
        ax.scatter(all_hnsc["A3_sum"], all_hnsc["SBS2"],
                   facecolors=COLOR_GRAY, edgecolors="#000000", linewidths=0.3,
                   s=100, alpha=0.5, zorder=2)

        ax.scatter(group_high["A3_sum"], group_high["SBS2"],
                   facecolors=COLOR_HIGH, edgecolors="#000000", linewidths=0.8,
                   s=100, alpha=0.85, zorder=3,
                   label=f"SBS2-HIGH (n={len(group_high)})" if ax == ax_main else "")

        ax.scatter(group_low["A3_sum"], group_low["SBS2"],
                   facecolors=COLOR_LOW, edgecolors="#000000", linewidths=0.8,
                   s=100, alpha=0.85, zorder=3,
                   label=f"SBS2-LOW (n={len(group_low)})" if ax == ax_main else "")

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(-y_pad * 0.5, y_max)
        ax.tick_params(axis="both", labelsize=FONT_TICK)

    ax_main.set_xlabel("A3A + A3B Expression (FPKM-UQ)", fontsize=FONT_AXIS, labelpad=12)
    ax_main.set_ylabel("SBS2 Mutational Weight", fontsize=FONT_AXIS, labelpad=12)
    ax_main.legend(fontsize=FONT_LEGEND - 4, loc="upper left", framealpha=0.9)

    ax_main.annotate("High A3\n(analysis zone)",
                     xy=(a3_median * 1.3, y_max_data * 0.85),
                     fontsize=20, fontweight="bold", color="#4D4D4D",
                     ha="left", fontstyle="italic")
    ax_main.annotate(f"(n = {len(all_hnsc)})",
                     xy=(0.98, 0.98), xycoords="axes fraction",
                     fontsize=FONT_LEGEND, ha="right", va="top")

    if ax_break is not None:
        ax_break.set_xlabel("")
        ax_break.tick_params(axis="y", left=False, labelleft=False)
        d = 0.015
        for ax_b, side in [(ax_main, "right"), (ax_break, "left")]:
            kwargs = dict(transform=ax_b.transAxes, color="#000000",
                          clip_on=False, lw=1.5)
            if side == "right":
                ax_b.plot((1-d, 1+d), (-d, +d), **kwargs)
                ax_b.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
                ax_b.spines["right"].set_visible(False)
            else:
                ax_b.plot((-d, +d), (-d, +d), **kwargs)
                ax_b.plot((-d, +d), (1-d, 1+d), **kwargs)
                ax_b.spines["left"].set_visible(False)

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(ct_fig_dir,
                    f"{cancer_type}_Panel_A_selection.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"[SAVE] Panel A -> {ct_fig_dir}")

    # =====================================================================
    # PANEL B -- Full Community Network (intra_score-based sizing)
    # =====================================================================
    banner("[PANEL B] Full community network (intra_score sizing)")

    if G_comm is not None and G_comm.number_of_nodes() > 0:

        # ---- Community-aware exploded layout (same algorithm as before) ----
        INTER_SCALE = 6.0
        INTRA_SCALE = 1.0
        ITERS_LOCAL  = 200
        ITERS_GLOBAL = 200

        rng = np.random.default_rng(COMMUNITY_BASE_SEED)

        # Build community super-graph
        CG = nx.Graph()
        for c in unique_comms:
            CG.add_node(c)
        for u, v, d in G_comm.edges(data=True):
            cu = gene_to_comm.get(u, -1)
            cv = gene_to_comm.get(v, -1)
            if cu == cv:
                continue
            w = abs(float(d.get("abs_weight", d.get("weight", 1.0))))
            if CG.has_edge(cu, cv):
                CG[cu][cv]["weight"] += w
            else:
                CG.add_edge(cu, cv, weight=w)

        if CG.number_of_edges() > 0:
            pos_comm = nx.spring_layout(
                CG, seed=COMMUNITY_BASE_SEED, weight="weight",
                k=2.0 / np.sqrt(max(CG.number_of_nodes(), 1)),
                iterations=ITERS_GLOBAL)
        else:
            pos_comm = nx.circular_layout(CG)

        pos_comm = {c: INTER_SCALE * np.asarray(p, dtype=float)
                    for c, p in pos_comm.items()}

        # Layout within each community
        pos = {}
        for c in unique_comms:
            c_nodes = comm_to_nodelist[c]
            sub = G_comm.subgraph(c_nodes).copy()
            if len(c_nodes) == 1:
                pos[c_nodes[0]] = pos_comm.get(c, np.zeros(2))
                continue
            if sub.number_of_edges() > 0:
                k_val = 2.0 / np.sqrt(max(sub.number_of_nodes(), 1))
                pos_sub = nx.spring_layout(sub, seed=COMMUNITY_BASE_SEED,
                                           weight="abs_weight",
                                           k=k_val, iterations=ITERS_LOCAL)
            else:
                pos_sub = {n: rng.normal(0, 0.05, size=2) for n in c_nodes}

            coords = np.array([pos_sub[n] for n in c_nodes], dtype=float)
            coords -= coords.mean(axis=0, keepdims=True)
            coords *= INTRA_SCALE
            center = pos_comm.get(c, np.zeros(2, dtype=float))
            for i, n in enumerate(c_nodes):
                pos[n] = center + coords[i]

        log(f"  Exploded layout: {len(pos)} nodes across {n_comms} communities")

        # ---- Sizing functions using intra_score ----
        def get_intra(n):
            return node_intra.get(n, 0.1)

        def node_size_full(n):
            """Full network node sizing."""
            return 100 + 1160 * get_intra(n)

        def label_fontsize_full(n):
            """Full network label sizing."""
            return 15 + 28 * get_intra(n)

        # ---- Plot ----
        fig, ax = plt.subplots(figsize=(26, 24))

        # === EDGES (thick, grant-style) ===
        edges = list(G_comm.edges(data=True))
        if edges:
            intra_edges = []
            inter_edges = []
            for u, v, d in edges:
                cu = gene_to_comm.get(u, -1)
                cv = gene_to_comm.get(v, -1)
                if cu == cv:
                    intra_edges.append((u, v, d))
                else:
                    inter_edges.append((u, v, d))

            # Intra-community: thick, visible
            for u, v, d in intra_edges:
                w = d.get("weight", 0)
                color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
                width = 0.8 + 4.0 * abs(w)
                alpha = min(0.55, 0.15 + 0.35 * abs(w))
                ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                        color=color, alpha=alpha, linewidth=width, zorder=1)

            # Inter-community: lighter, dashed
            for u, v, d in inter_edges:
                w = d.get("weight", 0)
                width = 0.4 + 1.5 * abs(w)
                ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                        color="#999999", alpha=0.10, linewidth=width,
                        linestyle=(0, (3, 3)), zorder=0)

        # === NODES (colored by community, sized by intra_score) ===
        a3_in_graph = [g for g in A3_GENES if g in G_comm.nodes()]

        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if not c_nodes:
                continue
            sizes = [node_size_full(n) for n in c_nodes]
            nx.draw_networkx_nodes(
                G_comm, pos, nodelist=c_nodes, ax=ax,
                node_size=sizes,
                node_color=[comm_color_map[c]] * len(c_nodes),
                alpha=0.85, linewidths=0.5, edgecolors="#000000",
                label=f"C{c} (n={len(c_nodes)})")

        # A3 gene nodes: prominent with ring
        if a3_in_graph:
            a3_sizes = [node_size_full(n) * 1.8 for n in a3_in_graph]
            nx.draw_networkx_nodes(
                G_comm, pos, nodelist=a3_in_graph, ax=ax,
                node_size=a3_sizes,
                node_color=COLOR_A3, alpha=1.0,
                edgecolors="#000000", linewidths=2.5)
            # Red ring
            ring_sizes = [s * 1.8 for s in a3_sizes]
            nx.draw_networkx_nodes(
                G_comm, pos, nodelist=a3_in_graph, ax=ax,
                node_size=ring_sizes, node_color="none",
                edgecolors=COLOR_A3, linewidths=4.0, alpha=1.0)

        # === LABELS (top 3 per community by intra_score + always A3) ===
        labels = {}
        font_sizes = {}
        font_colors = {}

        # A3 genes always labeled
        for n in a3_in_graph:
            labels[n] = sym(n)
            font_sizes[n] = FONT_AXIS
            font_colors[n] = COLOR_A3_TEXT

        # Top 3 hubs per community by intra_score
        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 5:
                continue
            scored = [(n, get_intra(n)) for n in c_nodes]
            scored.sort(key=lambda x: -x[1])
            for n, sc in scored[:3]:
                if n not in labels:
                    labels[n] = sym(n)
                    font_sizes[n] = label_fontsize_full(n)
                    font_colors[n] = "#000000"

        log(f"  Labeling {len(labels)} nodes (top hubs + A3 genes)")

        # Label repulsion
        label_pos = repel_labels(pos, labels, data_range=14.0,
                                 n_iters=300, push_strength=0.08)

        # Draw labels with background boxes
        draw_labels_with_boxes(ax, pos, label_pos, labels,
                               font_sizes, font_colors, char_w=0.07)

        # Community center labels (Cx annotation near each cluster)
        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 3:
                continue
            cx = np.mean([pos[n][0] for n in c_nodes])
            cy = np.max([pos[n][1] for n in c_nodes]) + 0.25
            ax.annotate(f"C{c}", xy=(cx, cy), fontsize=18, fontweight="bold",
                        color=comm_color_map[c], ha="center", va="bottom",
                        alpha=0.8)

        # Legend
        legend_handles = [
            Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3 gene"),
            Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.5,
                  label="Gained co-expression"),
            Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.5,
                  label="Lost co-expression"),
        ]
        ax.legend(handles=legend_handles, loc="upper left",
                  fontsize=FONT_LEGEND, framealpha=0.9,
                  title="Legend", title_fontsize=FONT_LEGEND + 2)

        ax.set_title(
            f"{cancer_type} | DIFF Network -- {G_comm.number_of_nodes()} genes, "
            f"{G_comm.number_of_edges()} edges, {n_comms} communities",
            fontsize=FONT_TITLE, pad=15)
        ax.axis("off")
        plt.tight_layout()

        for ext in ["pdf", "png"]:
            plt.savefig(os.path.join(ct_fig_dir,
                        f"{cancer_type}_Panel_B_network_full.{ext}"),
                        dpi=300, bbox_inches="tight")
        plt.close()
        log(f"[SAVE] Panel B (network) -> {ct_fig_dir}")

        # =================================================================
        # PANEL C -- Community 4 (A3B) Zoom
        # =================================================================
        banner("[PANEL C] Community 4 (A3B) zoom")

        # Find community containing A3B
        a3b_comm = None
        for n in a3_in_graph:
            if sym(n) == "APOBEC3B" or n == "APOBEC3B":
                a3b_comm = gene_to_comm.get(n, None)
                break
        if a3b_comm is None:
            # Fallback: check by symbol in partition
            for _, row in partition_df.iterrows():
                if sym(row["gene"]) == "APOBEC3B" or row["gene"] == "APOBEC3B":
                    a3b_comm = row["community"]
                    break

        if a3b_comm is not None:
            c4_nodes = comm_to_nodelist.get(a3b_comm, [])
            G_c4 = G_comm.subgraph(c4_nodes).copy()

            # Remove isolates
            isolates = [n for n in G_c4.nodes() if G_c4.degree(n) == 0]
            G_c4.remove_nodes_from(isolates)

            log(f"  A3B community: C{a3b_comm}, {G_c4.number_of_nodes()} nodes, "
                f"{G_c4.number_of_edges()} edges")

            if G_c4.number_of_nodes() > 0 and G_c4.number_of_edges() > 0:

                # Layout with expanded spring (grant style)
                pos_c4 = nx.spring_layout(
                    G_c4, seed=COMMUNITY_BASE_SEED, weight="abs_weight",
                    k=4.0 / np.sqrt(max(G_c4.number_of_nodes(), 1)),
                    iterations=500, scale=2.0)

                def node_size_zoom(n):
                    """Zoom panel node sizing."""
                    return 750 + 5000 * get_intra(n)

                def label_fontsize_zoom(n):
                    """Zoom panel label sizing."""
                    return 18 + 28 * get_intra(n)

                fig, ax = plt.subplots(figsize=(22, 20))

                # Edges (doubled thickness for zoom readability)
                for u, v, d in G_c4.edges(data=True):
                    w = d.get("weight", 0)
                    color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
                    width = 2.0 + 10.0 * abs(w)
                    alpha = min(0.6, 0.20 + 0.35 * abs(w))
                    ax.plot([pos_c4[u][0], pos_c4[v][0]],
                            [pos_c4[u][1], pos_c4[v][1]],
                            color=color, alpha=alpha, linewidth=width)

                # Regular nodes (community color, intra_score sizing)
                c4_regular = [n for n in G_c4.nodes() if n not in a3_in_graph]
                if c4_regular:
                    sizes = [node_size_zoom(n) for n in c4_regular]
                    nx.draw_networkx_nodes(
                        G_c4, pos_c4, nodelist=c4_regular, ax=ax,
                        node_size=sizes,
                        node_color=comm_color_map.get(a3b_comm, "#888888"),
                        alpha=0.90, edgecolors="#000000", linewidths=0.8)

                # A3 gene nodes
                a3_c4 = [n for n in a3_in_graph if n in G_c4.nodes()]
                if a3_c4:
                    a3_sizes = [node_size_zoom(n) * 1.5 for n in a3_c4]
                    nx.draw_networkx_nodes(
                        G_c4, pos_c4, nodelist=a3_c4, ax=ax,
                        node_size=a3_sizes,
                        node_color=COLOR_A3, alpha=1.0,
                        edgecolors="#000000", linewidths=2.5)
                    # Red ring
                    ring_sizes = [s * 2.0 for s in a3_sizes]
                    nx.draw_networkx_nodes(
                        G_c4, pos_c4, nodelist=a3_c4, ax=ax,
                        node_size=ring_sizes, node_color="none",
                        edgecolors=COLOR_A3, linewidths=5.0, alpha=1.0)

                # Labels: all nodes labeled in zoom (community is small)
                labels_c4 = {}
                fsizes_c4 = {}
                fcolors_c4 = {}

                for n in G_c4.nodes():
                    labels_c4[n] = sym(n)
                    if n in a3_c4:
                        fsizes_c4[n] = FONT_AXIS
                        fcolors_c4[n] = COLOR_A3_TEXT
                    else:
                        fsizes_c4[n] = label_fontsize_zoom(n)
                        fcolors_c4[n] = "#000000"

                # Label repulsion
                lp_c4 = repel_labels(pos_c4, labels_c4, data_range=5.0,
                                     n_iters=250, push_strength=0.07)

                draw_labels_with_boxes(ax, pos_c4, lp_c4, labels_c4,
                                       fsizes_c4, fcolors_c4, char_w=0.05)

                # Legend
                legend_c4 = [
                    Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3B"),
                    Patch(fc=comm_color_map.get(a3b_comm, "#888888"),
                          ec="#000000", lw=0.5,
                          label=f"Community {a3b_comm} genes"),
                    Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.5,
                          label="Gained co-expression"),
                    Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.5,
                          label="Lost co-expression"),
                ]
                ax.legend(handles=legend_c4, loc="upper left",
                          fontsize=FONT_LEGEND, framealpha=0.9,
                          title=f"Community {a3b_comm}",
                          title_fontsize=FONT_LEGEND + 2)

                ax.set_title(
                    f"Community {a3b_comm} (A3B) -- {G_c4.number_of_nodes()} genes, "
                    f"{G_c4.number_of_edges()} edges\n"
                    f"Top KEGG: Wnt signaling / Basal cell carcinoma",
                    fontsize=FONT_TITLE, pad=15)
                ax.axis("off")
                plt.tight_layout()

                for ext in ["pdf", "png"]:
                    plt.savefig(os.path.join(ct_fig_dir,
                                f"{cancer_type}_Panel_C_C{a3b_comm}_zoom.{ext}"),
                                dpi=300, bbox_inches="tight")
                plt.close()
                log(f"[SAVE] Panel C (C{a3b_comm} zoom) -> {ct_fig_dir}")
        else:
            log("[WARNING] Could not identify A3B community for Panel C")

        # =================================================================
        # Per-community zooms (supplemental, intra_score sizing)
        # =================================================================
        banner("[SUPPLEMENT] Per-community zoomed plots")

        zoom_dir = ensure_dir(os.path.join(ct_fig_dir, "community_zooms"))

        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 5:
                continue

            G_sub = G_comm.subgraph(c_nodes).copy()
            if G_sub.number_of_edges() == 0:
                continue

            pos_sub = nx.spring_layout(
                G_sub, seed=COMMUNITY_BASE_SEED, weight="abs_weight",
                k=3.0 / np.sqrt(max(G_sub.number_of_nodes(), 1)),
                iterations=300)

            fig, ax = plt.subplots(figsize=(16, 14))

            # Edges (doubled thickness for zoom readability)
            for u, v, d in G_sub.edges(data=True):
                w = d.get("weight", 0)
                color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
                width = 1.6 + 8.0 * abs(w)
                alpha = min(0.5, 0.15 + 0.3 * abs(w))
                ax.plot([pos_sub[u][0], pos_sub[v][0]],
                        [pos_sub[u][1], pos_sub[v][1]],
                        color=color, alpha=alpha, linewidth=width)

            # Nodes (intra_score sizing)
            sizes = [500 + 3500 * get_intra(n) for n in G_sub.nodes()]
            nx.draw_networkx_nodes(
                G_sub, pos_sub, ax=ax, node_size=sizes,
                node_color=[comm_color_map[c]] * G_sub.number_of_nodes(),
                alpha=0.85, edgecolors="#000000", linewidths=0.8)

            # A3 gene highlighting
            a3_here = [g for g in A3_GENES if g in G_sub.nodes()]
            if a3_here:
                a3_s = [(500 + 3500 * get_intra(n)) * 1.5 for n in a3_here]
                nx.draw_networkx_nodes(
                    G_sub, pos_sub, nodelist=a3_here, ax=ax,
                    node_size=a3_s, node_color=COLOR_A3,
                    edgecolors="#000000", linewidths=2.0)

            # Labels (all for small communities, top hubs for large)
            zoom_labels = {}
            zoom_fs = {}
            zoom_fc = {}
            if len(c_nodes) <= 40:
                for n in G_sub.nodes():
                    zoom_labels[n] = sym(n)
                    zoom_fs[n] = 15 + 21 * get_intra(n)
                    zoom_fc[n] = COLOR_A3_TEXT if n in a3_here else "#000000"
            else:
                scored = [(n, get_intra(n)) for n in G_sub.nodes()]
                scored.sort(key=lambda x: -x[1])
                for n, sc in scored[:15]:
                    zoom_labels[n] = sym(n)
                    zoom_fs[n] = 15 + 21 * get_intra(n)
                    zoom_fc[n] = "#000000"
                for n in a3_here:
                    zoom_labels[n] = sym(n)
                    zoom_fs[n] = FONT_AXIS
                    zoom_fc[n] = COLOR_A3_TEXT

            if zoom_labels:
                zl_pos = repel_labels(pos_sub, zoom_labels, data_range=3.0,
                                      n_iters=250, push_strength=0.07)
                draw_labels_with_boxes(ax, pos_sub, zl_pos, zoom_labels,
                                       zoom_fs, zoom_fc, char_w=0.04)

            a3_note = ""
            a3_sym_here = [sym(g) for g in a3_here]
            if a3_sym_here:
                a3_note = f" -- {', '.join(a3_sym_here)}"

            ax.set_title(
                f"Community {c} -- {len(c_nodes)} genes, "
                f"{G_sub.number_of_edges()} edges{a3_note}",
                fontsize=FONT_AXIS, pad=10)
            ax.axis("off")
            plt.tight_layout()

            for ext in ["png", "pdf"]:
                plt.savefig(os.path.join(zoom_dir,
                            f"{cancer_type}_community_{c:02d}.{ext}"),
                            dpi=300, bbox_inches="tight")
            plt.close()
            log(f"  [SAVE] Community {c} ({len(c_nodes)} genes)")

        log(f"[SAVE] All community zooms -> {zoom_dir}")

    else:
        log("[SKIP] No community graph available for Panels B/C")

    # =====================================================================
    # SUPPLEMENT -- Methods Graphical Abstract
    # HIGH -> LOW -> DIFF heatmaps -> mini network
    # =====================================================================
    banner("[SUPPLEMENT] Methods graphical abstract")

    if corr_top is not None and corr_bot is not None and G_comm is not None:

        # Build community-ordered gene list
        ref_matrix = corr_top
        ordered_genes, comm_boundaries, comm_ids_sorted = \
            build_community_ordered_genes(partition_df, ref_matrix, gene_to_comm)
        log(f"  Ordered {len(ordered_genes)} genes across "
            f"{len(comm_ids_sorted)} communities")

        fig = plt.figure(figsize=(52, 13))
        gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.12)

        matrices_info = [
            ("HIGH (SBS2-high)", corr_top,  "plasma"),
            ("LOW (SBS2-low)",   corr_bot,  "plasma"),
            ("DIFF (HIGH - LOW)", corr_diff, "plasma"),
        ]

        for idx, (label, matrix, cmap_name) in enumerate(matrices_info):
            ax = fig.add_subplot(gs[0, idx])

            genes_present = [g for g in ordered_genes if g in matrix.index]
            heatmap_data = matrix.loc[genes_present, genes_present].values

            vmin = np.nanpercentile(heatmap_data, 1)
            vmax = np.nanpercentile(heatmap_data, 99)
            if vmax <= vmin:
                vmax = 1.0

            im = ax.imshow(heatmap_data, aspect="auto",
                           interpolation="nearest",
                           cmap=cmap_name, vmin=vmin, vmax=vmax)

            # Community boundary lines
            for start, end, c in comm_boundaries:
                if start > 0:
                    ax.axhline(start - 0.5, color="white", lw=2.5, alpha=0.9)
                    ax.axvline(start - 0.5, color="white", lw=2.5, alpha=0.9)

            # Community labels
            for start, end, c in comm_boundaries:
                mid = (start + end) / 2
                n_genes = end - start
                if n_genes >= 8:
                    ax.annotate(f"C{c}", xy=(len(genes_present) + 3, mid),
                                fontsize=FONT_LEGEND, fontweight="bold",
                                va="center", ha="left",
                                annotation_clip=False)

            cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.10)
            cbar.set_label(f"Spearman rho", fontsize=FONT_LEGEND, labelpad=8)
            cbar.ax.tick_params(labelsize=FONT_TICK)

            ax.set_title(label, fontsize=FONT_TITLE, pad=12)
            ax.set_xticks([])
            ax.set_yticks([])

            # Operator labels between panels
            if idx == 0:
                ax.annotate("  -  ", xy=(1.06, 0.5),
                            xycoords="axes fraction",
                            fontsize=50, fontweight="bold", color="#333333",
                            ha="center", va="center")
            elif idx == 1:
                ax.annotate("  =  ", xy=(1.06, 0.5),
                            xycoords="axes fraction",
                            fontsize=50, fontweight="bold", color="#333333",
                            ha="center", va="center")
            elif idx == 2:
                ax.annotate("  >>  ", xy=(1.06, 0.5),
                            xycoords="axes fraction",
                            fontsize=40, fontweight="bold", color="#333333",
                            ha="center", va="center")

        # Mini network in 4th panel
        ax_net = fig.add_subplot(gs[0, 3])

        # Use same exploded layout positions, just rescale for this axes
        if pos:
            all_xy = np.array(list(pos.values()))
            xmin, ymin = all_xy.min(axis=0)
            xmax, ymax = all_xy.max(axis=0)

            # Draw edges (simplified)
            for u, v, d in G_comm.edges(data=True):
                cu = gene_to_comm.get(u, -1)
                cv = gene_to_comm.get(v, -1)
                if cu != cv:
                    continue
                w = d.get("weight", 0)
                color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
                ax_net.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                            color=color, alpha=0.15, linewidth=0.5)

            # Draw nodes
            for c in unique_comms:
                c_nodes = comm_to_nodelist.get(c, [])
                if not c_nodes:
                    continue
                xs = [pos[n][0] for n in c_nodes]
                ys = [pos[n][1] for n in c_nodes]
                sizes = [15 + 80 * get_intra(n) for n in c_nodes]
                ax_net.scatter(xs, ys, s=sizes,
                               c=comm_color_map[c], alpha=0.8,
                               edgecolors="#000000", linewidths=0.2)

            # A3 genes
            for n in a3_in_graph:
                ax_net.scatter(pos[n][0], pos[n][1], s=120,
                               c=COLOR_A3, edgecolors="#000000",
                               linewidths=1.5, zorder=5)

            # Community labels
            for c in unique_comms:
                c_nodes = comm_to_nodelist.get(c, [])
                if len(c_nodes) < 5:
                    continue
                cx = np.mean([pos[n][0] for n in c_nodes])
                cy = np.max([pos[n][1] for n in c_nodes]) + 0.15
                ax_net.annotate(f"C{c}", xy=(cx, cy), fontsize=14,
                                fontweight="bold",
                                color=comm_color_map[c],
                                ha="center", va="bottom", alpha=0.9)

        ax_net.set_title("Network", fontsize=FONT_TITLE, pad=12)
        ax_net.axis("off")

        plt.savefig(os.path.join(ct_fig_dir,
                    f"{cancer_type}_Supplement_methods_overview.png"),
                    dpi=300, bbox_inches="tight")
        plt.savefig(os.path.join(ct_fig_dir,
                    f"{cancer_type}_Supplement_methods_overview.pdf"),
                    bbox_inches="tight")
        plt.close()
        log(f"[SAVE] Supplement methods overview -> {ct_fig_dir}")

    else:
        missing = []
        if corr_top is None:
            missing.append("TOP")
        if corr_bot is None:
            missing.append("BOTTOM")
        if G_comm is None:
            missing.append("graph")
        log(f"[SKIP] Supplement: missing {', '.join(missing)}")

    # =====================================================================
    # SUMMARY
    # =====================================================================
    banner(f"FIGURE 2 COMPLETE -- {cancer_type}")
    print(f"\n  Output directory: {ct_fig_dir}")
    print(f"  Panel A: Selection plot (PDF + PNG)")
    print(f"  Panel B: Full network with intra_score sizing (PDF + PNG)")
    print(f"  Panel C: A3B community zoom (PDF + PNG)")
    print(f"  Supplement: Methods graphical abstract (PDF + PNG)")
    print(f"  Community zooms: All communities (PNG + PDF)")

banner("ALL PANELS COMPLETE")
