#!/usr/bin/env python3
"""
Generate_Figure2_Panels.py
===========================

Generate publication-quality Figure 2 panels from V4 pipeline outputs.
Reads from Steps 03-06 outputs + node importance scores and produces:

  Panel A  -- Methods graphical abstract: HIGH - LOW = DIFF -> network
  Panel B  -- Full community network (two-tier layout: LCC + satellites)
  Panel B2 -- LCC-only version (fallback if satellites clutter)
  Panel C  -- Community zoom insets (A3B community, epithelial signaling, etc.)
  Supplement -- SBS2 vs A3 sum selection plot

Style: 28-34pt text, all hex color codes, PDF + PNG (300 DPI).

Run AFTER the full pipeline (Steps 01-08) has completed.

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
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from network_config import (
    CANCER_TYPES, A3_GENES, A3_ID_TO_ALIAS, BIOMARKERS,
    COMMUNITY_BASE_SEED, MIN_COMMUNITY_SIZE,
    DIR_01_CLEANED, DIR_02_MERGED, DIR_03_DIFFEXPR,
    DIR_04_NETWORKS, DIR_05_COMMUNITIES, DIR_06_CENTRALITY,
    FIG2_ROOT, banner, log, ensure_dir
)


# =============================================================================
# SHARED SETTINGS (hard-won, do not change)
# =============================================================================

with open(os.path.join(DIR_01_CLEANED, "ensg_to_symbol.json")) as f:
    ensg_to_symbol = json.load(f)

def sym(ensg):
    return ensg_to_symbol.get(str(ensg), str(ensg))

# Color palette
COLOR_TEAL  = "#9bc1bc"
COLOR_CORAL = "#ed6a5a"
COLOR_CREAM = "#f4f1bb"
COLOR_HIGH  = "#ed6a5a"
COLOR_LOW   = "#5b8e7d"
COLOR_GRAY  = "#d0d0d0"
COLOR_EDGE_POS = "#b22222"
COLOR_EDGE_NEG = "#4682b4"
COLOR_A3       = "#ed6a5a"
COLOR_A3_TEXT  = "#c0392b"
COLOR_SAT      = "#e8913a"   # orange for satellite communities

# Font sizes (28-34 range)
FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22
LABEL_BG_ALPHA = 0.42

FIG_DIR = ensure_dir(os.path.join(FIG2_ROOT, "FIGURE_2_PANELS"))


# =============================================================================
# LABEL REPULSION (ported from grant script, do not change)
# =============================================================================

def repel_labels(pos, labels, data_range=12.0, n_iters=120, push_strength=0.04):
    label_pos = {n: np.array(pos[n], dtype=float) + np.array([0, 0.16])
                 for n in labels}
    char_w = data_range * 0.012
    for _ in range(n_iters):
        keys = list(label_pos.keys())
        for i, n1 in enumerate(keys):
            for n2 in keys[i + 1:]:
                p1, p2 = label_pos[n1], label_pos[n2]
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
                      fc="white", ec="none", alpha=LABEL_BG_ALPHA))
        if offset > char_w * 2.5:
            ax.annotate(labels[n], xy=nxy, xytext=lxy, **kwargs,
                        arrowprops=dict(arrowstyle="-", color="#888888",
                                        lw=0.6, alpha=0.5))
        else:
            ax.annotate(labels[n], xy=lxy, **kwargs)


# =============================================================================
# TWO-TIER LAYOUT: LCC communities spread out, satellites in ring
# =============================================================================

def two_tier_layout(G_comm, gene_to_comm, min_community_size=10,
                    inter_scale=8.0, intra_scale=1.2, sat_ring_radius=14.0,
                    seed=42):
    """
    Build exploded layout with two tiers:
      1. Main communities (in LCC or large components) get spring layout
         with generous spacing
      2. Satellite communities (entire small components) placed in a ring
         around the outside

    Returns dict {node: (x, y)} for all nodes.
    """
    rng = np.random.default_rng(seed)

    # Classify communities
    components = list(nx.connected_components(G_comm))
    node_to_comp_size = {}
    for comp in components:
        sz = len(comp)
        for n in comp:
            node_to_comp_size[n] = sz

    comm_to_nodelist = {}
    for n in G_comm.nodes():
        c = gene_to_comm.get(n, -1)
        comm_to_nodelist.setdefault(c, []).append(n)

    main_comms = []
    sat_comms = []
    for c, nodes in comm_to_nodelist.items():
        max_comp = max(node_to_comp_size.get(n, 1) for n in nodes)
        if max_comp >= min_community_size:
            main_comms.append(c)
        else:
            sat_comms.append(c)

    main_comms.sort(key=lambda c: len(comm_to_nodelist[c]), reverse=True)
    sat_comms.sort(key=lambda c: len(comm_to_nodelist[c]), reverse=True)

    # --- Tier 1: Main communities with spring layout ---
    # Build super-graph of main communities
    CG = nx.Graph()
    for c in main_comms:
        CG.add_node(c)
    for u, v, d in G_comm.edges(data=True):
        cu = gene_to_comm.get(u, -1)
        cv = gene_to_comm.get(v, -1)
        if cu == cv or cu not in main_comms or cv not in main_comms:
            continue
        w = abs(float(d.get("abs_weight", d.get("weight", 1.0))))
        if CG.has_edge(cu, cv):
            CG[cu][cv]["weight"] += w
        else:
            CG.add_edge(cu, cv, weight=w)

    if CG.number_of_edges() > 0:
        pos_comm = nx.spring_layout(
            CG, seed=seed, weight="weight",
            k=3.0 / np.sqrt(max(CG.number_of_nodes(), 1)),
            iterations=300)
    else:
        pos_comm = nx.circular_layout(CG)

    pos_comm = {c: inter_scale * np.asarray(p, dtype=float)
                for c, p in pos_comm.items()}

    # --- Tier 2: Satellites in a ring ---
    if sat_comms:
        n_sat = len(sat_comms)
        for i, c in enumerate(sat_comms):
            angle = 2 * np.pi * i / n_sat
            pos_comm[c] = np.array([
                sat_ring_radius * np.cos(angle),
                sat_ring_radius * np.sin(angle)
            ])

    # --- Layout within each community ---
    pos = {}
    for c, c_nodes in comm_to_nodelist.items():
        center = pos_comm.get(c, np.zeros(2))
        if len(c_nodes) == 1:
            pos[c_nodes[0]] = center
            continue

        sub = G_comm.subgraph(c_nodes).copy()
        if sub.number_of_edges() > 0:
            k_val = 2.5 / np.sqrt(max(sub.number_of_nodes(), 1))
            pos_sub = nx.spring_layout(sub, seed=seed, weight="abs_weight",
                                       k=k_val, iterations=200)
        else:
            pos_sub = {n: rng.normal(0, 0.05, size=2) for n in c_nodes}

        coords = np.array([pos_sub[n] for n in c_nodes], dtype=float)
        coords -= coords.mean(axis=0, keepdims=True)

        # Scale main communities larger, satellites smaller
        if c in main_comms:
            coords *= intra_scale
        else:
            coords *= 0.3  # compact satellites

        for i, n in enumerate(c_nodes):
            pos[n] = center + coords[i]

    return pos, main_comms, sat_comms, comm_to_nodelist


# =============================================================================
# COMMUNITY-ORDERED GENE LIST (for heatmaps)
# =============================================================================

def build_community_ordered_genes(partition_df, ref_matrix, gene_to_comm):
    comm_genes = [g for g in partition_df["gene"] if g in ref_matrix.index]
    comm_to_genes = {}
    for g in comm_genes:
        c = gene_to_comm.get(g, -1)
        comm_to_genes.setdefault(c, []).append(g)

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
        try:
            Z_c = linkage(squareform(cd, checks=False), method="average")
            comm_ids_sorted = [comm_ids[i] for i in leaves_list(Z_c)]
        except Exception:
            comm_ids_sorted = sorted(comm_ids,
                                     key=lambda c: len(clustered[c]),
                                     reverse=True)
    else:
        comm_ids_sorted = sorted(comm_ids,
                                 key=lambda c: len(clustered[c]),
                                 reverse=True)

    ordered_genes = []
    comm_boundaries = []
    for c in comm_ids_sorted:
        start = len(ordered_genes)
        ordered_genes.extend(clustered[c])
        comm_boundaries.append((start, len(ordered_genes), c))

    return ordered_genes, comm_boundaries, comm_ids_sorted


# =============================================================================
# MAIN LOOP
# =============================================================================

for cancer_type in CANCER_TYPES:

    banner(f"FIGURE 2 -- {cancer_type}", char="=")
    ct_fig_dir = ensure_dir(os.path.join(FIG_DIR, cancer_type))

    # ---- Load data ----
    merged = pd.read_pickle(os.path.join(DIR_02_MERGED,
                            "TCGA_merged_expression_SBS.pkl"))
    all_hnsc = merged[merged["Project_ID"] == cancer_type].copy()
    for c in ["SBS2", "A3A", "A3B"]:
        all_hnsc[c] = pd.to_numeric(all_hnsc[c], errors="coerce")
    all_hnsc = all_hnsc.dropna(subset=["SBS2", "A3A", "A3B"])
    all_hnsc["A3_sum"] = all_hnsc["A3A"] + all_hnsc["A3B"]
    log(f"All HNSC tumors: {len(all_hnsc)}")

    # Groups
    high_path = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                             f"{cancer_type}_SBS2_HIGH_group.pkl")
    low_path  = os.path.join(DIR_03_DIFFEXPR, cancer_type,
                             f"{cancer_type}_SBS2_LOW_group.pkl")
    if not os.path.exists(high_path):
        log(f"[SKIP] Groups not found"); continue
    group_high = pd.read_pickle(high_path)
    group_low  = pd.read_pickle(low_path)
    log(f"HIGH: {len(group_high)}, LOW: {len(group_low)}")

    # Correlation matrices
    corr_dir = os.path.join(DIR_04_NETWORKS, cancer_type, "corr_matrices")
    corr_diff = pd.read_pickle(os.path.join(corr_dir,
                               f"{cancer_type}_corr_DIFF.pkl"))
    corr_top_path = os.path.join(corr_dir, f"{cancer_type}_corr_TOP.pkl")
    corr_bot_path = os.path.join(corr_dir, f"{cancer_type}_corr_BOTTOM.pkl")
    corr_top = pd.read_pickle(corr_top_path) if os.path.exists(corr_top_path) else None
    corr_bot = pd.read_pickle(corr_bot_path) if os.path.exists(corr_bot_path) else None
    log(f"DIFF corr matrix: {corr_diff.shape}")

    # Partition
    part_path = os.path.join(DIR_05_COMMUNITIES, cancer_type,
                             f"{cancer_type}_best_partition.csv")
    if not os.path.exists(part_path):
        log("[SKIP] Partition not found"); continue
    partition_df = pd.read_csv(part_path)
    gene_to_comm = dict(zip(partition_df["gene"], partition_df["community"]))

    # Detect satellite communities from partition
    is_satellite = {}
    if "is_satellite" in partition_df.columns:
        for _, row in partition_df.iterrows():
            is_satellite[row["community"]] = bool(row["is_satellite"])
    log(f"Partition: {len(gene_to_comm)} genes, "
        f"{partition_df['community'].nunique()} communities")

    # Community graph
    comm_graph_path = os.path.join(DIR_05_COMMUNITIES, cancer_type,
                                   f"{cancer_type}_G_comm.gpickle")
    if not os.path.exists(comm_graph_path):
        log("[SKIP] Community graph not found"); continue
    with open(comm_graph_path, "rb") as f:
        G_comm = pickle.load(f)
    log(f"Community graph: {G_comm.number_of_nodes()} nodes, "
        f"{G_comm.number_of_edges()} edges")

    # Node importance scores
    sizing_path = os.path.join(FIG2_ROOT, "DIAGNOSTIC_AUDIT",
                               "Figure2_Node_Sizing.tsv")
    node_intra = {}
    if os.path.exists(sizing_path):
        sizing_df = pd.read_csv(sizing_path, sep="\t")
        for _, row in sizing_df.iterrows():
            if "ensg_id" in sizing_df.columns:
                node_intra[row["ensg_id"]] = float(row["intra_score"])
            if "gene_symbol" in sizing_df.columns:
                node_intra[row["gene_symbol"]] = float(row["intra_score"])
        log(f"Node scores loaded: {len(sizing_df)} genes")
    else:
        log("[WARNING] Node_Sizing.tsv not found, using degree fallback")

    def get_intra(n):
        return node_intra.get(n, 0.1)

    # Build layouts
    pos_full, main_comms, sat_comms, comm_to_nodelist = two_tier_layout(
        G_comm, gene_to_comm,
        min_community_size=MIN_COMMUNITY_SIZE,
        inter_scale=8.0, intra_scale=1.2, sat_ring_radius=14.0,
        seed=COMMUNITY_BASE_SEED)

    unique_comms = sorted(comm_to_nodelist.keys())
    n_comms = len(unique_comms)
    cmap_comms = plt.colormaps["tab20"]
    comm_color_map = {}
    for i, c in enumerate(unique_comms):
        if c in sat_comms:
            comm_color_map[c] = COLOR_SAT
        else:
            comm_color_map[c] = mcolors.to_hex(cmap_comms(i % 20))

    log(f"Layout: {len(main_comms)} main + {len(sat_comms)} satellite")

    # =================================================================
    # PANEL A -- Methods Overview (HIGH - LOW = DIFF -> network)
    # =================================================================
    banner("[PANEL A] Methods graphical abstract")

    if corr_top is not None and corr_bot is not None:
        ref_matrix = corr_top
        ordered_genes, comm_boundaries, comm_ids_sorted = \
            build_community_ordered_genes(partition_df, ref_matrix,
                                          gene_to_comm)
        log(f"  Ordered {len(ordered_genes)} genes across "
            f"{len(comm_ids_sorted)} communities")

        fig = plt.figure(figsize=(52, 13))
        gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.12)

        matrices_info = [
            ("HIGH (SBS2-high)", corr_top, "plasma"),
            ("LOW (SBS2-low)", corr_bot, "plasma"),
            ("DIFF (HIGH - LOW)", corr_diff, "plasma"),
        ]

        for idx, (label, matrix, cmap_name) in enumerate(matrices_info):
            ax = fig.add_subplot(gs[0, idx])
            genes_present = [g for g in ordered_genes if g in matrix.index]
            hdata = matrix.loc[genes_present, genes_present].values
            vmin = np.nanpercentile(hdata, 1)
            vmax = np.nanpercentile(hdata, 99)
            if vmax <= vmin:
                vmax = 1.0

            im = ax.imshow(hdata, aspect="auto", interpolation="nearest",
                           cmap=cmap_name, vmin=vmin, vmax=vmax)

            # Community boundaries
            for start, end, c in comm_boundaries:
                if start > 0:
                    ax.axhline(start - 0.5, color="white", lw=2.5, alpha=0.9)
                    ax.axvline(start - 0.5, color="white", lw=2.5, alpha=0.9)

            for start, end, c in comm_boundaries:
                mid = (start + end) / 2
                if end - start >= 8:
                    ax.annotate(f"C{c}", xy=(len(genes_present) + 3, mid),
                                fontsize=FONT_LEGEND, fontweight="bold",
                                va="center", ha="left",
                                annotation_clip=False)

            cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.10)
            cbar.set_label("Spearman rho", fontsize=FONT_LEGEND, labelpad=8)
            cbar.ax.tick_params(labelsize=FONT_TICK)
            ax.set_title(label, fontsize=FONT_TITLE, pad=12)
            ax.set_xticks([]); ax.set_yticks([])

            # Operator labels
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

        # Mini network in 4th panel (main communities only for clarity)
        ax_net = fig.add_subplot(gs[0, 3])
        for u, v, d in G_comm.edges(data=True):
            cu = gene_to_comm.get(u, -1)
            cv = gene_to_comm.get(v, -1)
            if cu != cv:
                continue
            if cu not in main_comms:
                continue
            w = d.get("weight", 0)
            color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
            ax_net.plot([pos_full[u][0], pos_full[v][0]],
                        [pos_full[u][1], pos_full[v][1]],
                        color=color, alpha=0.12, linewidth=0.4)

        for c in main_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if not c_nodes:
                continue
            xs = [pos_full[n][0] for n in c_nodes]
            ys = [pos_full[n][1] for n in c_nodes]
            sizes = [12 + 70 * get_intra(n) for n in c_nodes]
            ax_net.scatter(xs, ys, s=sizes, c=comm_color_map[c],
                           alpha=0.8, edgecolors="#000000", linewidths=0.15)

        a3_in_graph = [g for g in A3_GENES if g in G_comm.nodes()]
        for n in a3_in_graph:
            if n in pos_full:
                ax_net.scatter(pos_full[n][0], pos_full[n][1], s=100,
                               c=COLOR_A3, edgecolors="#000000",
                               linewidths=1.2, zorder=5)

        ax_net.set_title("Network", fontsize=FONT_TITLE, pad=12)
        ax_net.axis("off")

        for ext in ["pdf", "png"]:
            plt.savefig(os.path.join(ct_fig_dir,
                        f"{cancer_type}_Panel_A_methods_overview.{ext}"),
                        dpi=300, bbox_inches="tight")
        plt.close()
        log(f"[SAVE] Panel A -> {ct_fig_dir}")

    # =================================================================
    # PANEL B -- Full Network (two-tier layout, all communities)
    # =================================================================
    banner("[PANEL B] Full network (two-tier layout)")

    def draw_network(ax, G, pos, comms_to_show, label_top_n=3,
                     node_size_fn=None, label_fs_fn=None,
                     show_sat_label=True):
        """Draw network on axes with specified communities."""
        if node_size_fn is None:
            node_size_fn = lambda n: 100 + 1160 * get_intra(n)
        if label_fs_fn is None:
            label_fs_fn = lambda n: 15 + 28 * get_intra(n)

        comms_set = set(comms_to_show)

        # Edges (intra only for clarity, inter as faint dashed)
        for u, v, d in G.edges(data=True):
            cu = gene_to_comm.get(u, -1)
            cv = gene_to_comm.get(v, -1)
            if cu not in comms_set and cv not in comms_set:
                continue
            w = d.get("weight", 0)
            if cu == cv:
                color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
                width = 0.8 + 4.0 * abs(w)
                alpha = min(0.55, 0.15 + 0.35 * abs(w))
                ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                        color=color, alpha=alpha, linewidth=width, zorder=1)
            else:
                if cu in comms_set and cv in comms_set:
                    width = 0.4 + 1.5 * abs(w)
                    ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                            color="#999999", alpha=0.08, linewidth=width,
                            linestyle=(0, (3, 3)), zorder=0)

        # Nodes
        for c in comms_to_show:
            c_nodes = comm_to_nodelist.get(c, [])
            if not c_nodes:
                continue
            sizes = [node_size_fn(n) for n in c_nodes]
            nx.draw_networkx_nodes(
                G, pos, nodelist=c_nodes, ax=ax,
                node_size=sizes,
                node_color=[comm_color_map[c]] * len(c_nodes),
                alpha=0.85, linewidths=0.5, edgecolors="#000000",
                label=f"C{c} (n={len(c_nodes)})")

        # A3 gene highlighting
        a3_here = [g for g in A3_GENES if g in G.nodes()
                    and gene_to_comm.get(g, -1) in comms_set]
        if a3_here:
            a3_sizes = [node_size_fn(n) * 1.8 for n in a3_here]
            nx.draw_networkx_nodes(
                G, pos, nodelist=a3_here, ax=ax,
                node_size=a3_sizes, node_color=COLOR_A3,
                alpha=1.0, edgecolors="#000000", linewidths=2.5)
            ring_sizes = [s * 1.8 for s in a3_sizes]
            nx.draw_networkx_nodes(
                G, pos, nodelist=a3_here, ax=ax,
                node_size=ring_sizes, node_color="none",
                edgecolors=COLOR_A3, linewidths=4.0, alpha=1.0)

        # Labels
        labels = {}
        font_sizes = {}
        font_colors = {}

        for n in a3_here:
            labels[n] = sym(n)
            font_sizes[n] = FONT_AXIS
            font_colors[n] = COLOR_A3_TEXT

        for c in comms_to_show:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 5:
                if show_sat_label and len(c_nodes) >= 2:
                    for n in c_nodes:
                        if n not in labels:
                            labels[n] = sym(n)
                            font_sizes[n] = 16
                            font_colors[n] = "#555555"
                continue
            scored = [(n, get_intra(n)) for n in c_nodes]
            scored.sort(key=lambda x: -x[1])
            for n, sc in scored[:label_top_n]:
                if n not in labels:
                    labels[n] = sym(n)
                    font_sizes[n] = label_fs_fn(n)
                    font_colors[n] = "#000000"

        if labels:
            label_pos = repel_labels(pos, labels, data_range=14.0,
                                     n_iters=300, push_strength=0.08)
            draw_labels_with_boxes(ax, pos, label_pos, labels,
                                   font_sizes, font_colors, char_w=0.07)

        # Community center labels
        for c in comms_to_show:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 3:
                continue
            cx = np.mean([pos[n][0] for n in c_nodes])
            cy = np.max([pos[n][1] for n in c_nodes]) + 0.25
            ax.annotate(f"C{c}", xy=(cx, cy), fontsize=18,
                        fontweight="bold", color=comm_color_map[c],
                        ha="center", va="bottom", alpha=0.8)

    # Full network (all communities)
    fig, ax = plt.subplots(figsize=(28, 26))
    draw_network(ax, G_comm, pos_full, unique_comms,
                 label_top_n=3, show_sat_label=True)

    legend_handles = [
        Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3 gene"),
        Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.5,
              label="Gained co-expression"),
        Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.5,
              label="Lost co-expression"),
        Patch(fc=COLOR_SAT, ec="#000000", lw=0.5,
              label="Satellite community"),
    ]
    ax.legend(handles=legend_handles, loc="upper left",
              fontsize=FONT_LEGEND, framealpha=0.9,
              title="Legend", title_fontsize=FONT_LEGEND + 2)
    ax.set_title(
        f"{cancer_type} | DIFF Network -- {G_comm.number_of_nodes()} genes, "
        f"{G_comm.number_of_edges()} edges\n"
        f"{len(main_comms)} main + {len(sat_comms)} satellite communities",
        fontsize=FONT_TITLE, pad=15)
    ax.axis("off")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(ct_fig_dir,
                    f"{cancer_type}_Panel_B_network_full.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"[SAVE] Panel B (full network) -> {ct_fig_dir}")

    # LCC-only version (backup)
    fig, ax = plt.subplots(figsize=(26, 24))
    draw_network(ax, G_comm, pos_full, main_comms,
                 label_top_n=3, show_sat_label=False)

    legend_handles_lcc = [
        Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3 gene"),
        Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.5,
              label="Gained co-expression"),
        Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.5,
              label="Lost co-expression"),
    ]
    ax.legend(handles=legend_handles_lcc, loc="upper left",
              fontsize=FONT_LEGEND, framealpha=0.9)
    ax.set_title(
        f"{cancer_type} | DIFF Network (main communities only) -- "
        f"{sum(len(comm_to_nodelist[c]) for c in main_comms)} genes",
        fontsize=FONT_TITLE, pad=15)
    ax.axis("off")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(ct_fig_dir,
                    f"{cancer_type}_Panel_B2_network_LCC.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"[SAVE] Panel B2 (LCC only) -> {ct_fig_dir}")

    # =================================================================
    # PANEL C -- Community Zoom Insets
    # =================================================================
    banner("[PANEL C] Community zoom insets")

    # Identify communities to zoom on (mentioned in text)
    zoom_targets = {}

    # A3B community
    for n in A3_GENES:
        if sym(n) == "APOBEC3B" and n in gene_to_comm:
            a3b_comm = gene_to_comm[n]
            zoom_targets[a3b_comm] = f"A3B community (C{a3b_comm})"
            break

    # Epithelial signaling community (C8 in latest run, 16 genes)
    # Find the community with CXCL2/CXCL3
    for g, c in gene_to_comm.items():
        if sym(g) in ["CXCL2", "CXCL3"]:
            if c not in zoom_targets:
                zoom_targets[c] = f"Epithelial signaling (C{c})"
            break

    # ZNF/HSV1 community (C6, 33 genes)
    for g, c in gene_to_comm.items():
        if sym(g) in ["ZNF571", "ZNF175"]:
            if c not in zoom_targets:
                zoom_targets[c] = f"ZNF cluster (C{c})"
            break

    zoom_dir = ensure_dir(os.path.join(ct_fig_dir, "community_zooms"))

    for c, desc in zoom_targets.items():
        c_nodes = comm_to_nodelist.get(c, [])
        G_sub = G_comm.subgraph(c_nodes).copy()
        isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
        G_sub.remove_nodes_from(isolates)

        if G_sub.number_of_nodes() == 0:
            log(f"  [SKIP] {desc}: no connected nodes")
            continue

        log(f"  Zooming: {desc} ({G_sub.number_of_nodes()} nodes, "
            f"{G_sub.number_of_edges()} edges)")

        pos_zoom = nx.spring_layout(
            G_sub, seed=COMMUNITY_BASE_SEED, weight="abs_weight",
            k=4.0 / np.sqrt(max(G_sub.number_of_nodes(), 1)),
            iterations=500, scale=2.0)

        def ns_zoom(n): return 750 + 5000 * get_intra(n)
        def fs_zoom(n): return 18 + 28 * get_intra(n)

        fig, ax = plt.subplots(figsize=(22, 20))

        # Edges
        for u, v, d in G_sub.edges(data=True):
            w = d.get("weight", 0)
            color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
            width = 2.0 + 10.0 * abs(w)
            alpha = min(0.6, 0.20 + 0.35 * abs(w))
            ax.plot([pos_zoom[u][0], pos_zoom[v][0]],
                    [pos_zoom[u][1], pos_zoom[v][1]],
                    color=color, alpha=alpha, linewidth=width)

        # Nodes
        regular = [n for n in G_sub.nodes() if n not in a3_in_graph]
        if regular:
            sizes = [ns_zoom(n) for n in regular]
            nx.draw_networkx_nodes(
                G_sub, pos_zoom, nodelist=regular, ax=ax,
                node_size=sizes,
                node_color=comm_color_map.get(c, "#888888"),
                alpha=0.90, edgecolors="#000000", linewidths=0.8)

        a3_sub = [n for n in a3_in_graph if n in G_sub.nodes()]
        if a3_sub:
            a3_sizes = [ns_zoom(n) * 1.5 for n in a3_sub]
            nx.draw_networkx_nodes(
                G_sub, pos_zoom, nodelist=a3_sub, ax=ax,
                node_size=a3_sizes, node_color=COLOR_A3,
                alpha=1.0, edgecolors="#000000", linewidths=2.5)
            ring_sizes = [s * 2.0 for s in a3_sizes]
            nx.draw_networkx_nodes(
                G_sub, pos_zoom, nodelist=a3_sub, ax=ax,
                node_size=ring_sizes, node_color="none",
                edgecolors=COLOR_A3, linewidths=5.0, alpha=1.0)

        # Labels (all for small communities, top hubs for large)
        z_labels = {}; z_fs = {}; z_fc = {}
        if len(list(G_sub.nodes())) <= 40:
            for n in G_sub.nodes():
                z_labels[n] = sym(n)
                z_fs[n] = fs_zoom(n) if n not in a3_sub else FONT_AXIS
                z_fc[n] = COLOR_A3_TEXT if n in a3_sub else "#000000"
        else:
            scored = [(n, get_intra(n)) for n in G_sub.nodes()]
            scored.sort(key=lambda x: -x[1])
            for n, sc in scored[:15]:
                z_labels[n] = sym(n)
                z_fs[n] = fs_zoom(n)
                z_fc[n] = "#000000"
            for n in a3_sub:
                z_labels[n] = sym(n)
                z_fs[n] = FONT_AXIS
                z_fc[n] = COLOR_A3_TEXT

        if z_labels:
            zl_pos = repel_labels(pos_zoom, z_labels, data_range=5.0,
                                  n_iters=250, push_strength=0.07)
            draw_labels_with_boxes(ax, pos_zoom, zl_pos, z_labels,
                                   z_fs, z_fc, char_w=0.05)

        ax.set_title(f"{desc}\n{G_sub.number_of_nodes()} genes, "
                     f"{G_sub.number_of_edges()} edges",
                     fontsize=FONT_TITLE, pad=15)
        ax.axis("off")
        plt.tight_layout()

        for ext in ["pdf", "png"]:
            plt.savefig(os.path.join(zoom_dir,
                        f"{cancer_type}_zoom_C{c}_{desc.split('(')[0].strip().replace(' ','_')}.{ext}"),
                        dpi=300, bbox_inches="tight")
        plt.close()
        log(f"  [SAVE] {desc}")

    # Also generate ALL community zooms (supplemental)
    banner("[SUPPLEMENT] All community zooms")
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
        for u, v, d in G_sub.edges(data=True):
            w = d.get("weight", 0)
            color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
            width = 1.6 + 8.0 * abs(w)
            alpha = min(0.5, 0.15 + 0.3 * abs(w))
            ax.plot([pos_sub[u][0], pos_sub[v][0]],
                    [pos_sub[u][1], pos_sub[v][1]],
                    color=color, alpha=alpha, linewidth=width)

        sizes = [500 + 3500 * get_intra(n) for n in G_sub.nodes()]
        nx.draw_networkx_nodes(
            G_sub, pos_sub, ax=ax, node_size=sizes,
            node_color=[comm_color_map[c]] * G_sub.number_of_nodes(),
            alpha=0.85, edgecolors="#000000", linewidths=0.8)

        a3_here = [g for g in A3_GENES if g in G_sub.nodes()]
        if a3_here:
            a3_s = [(500 + 3500 * get_intra(n)) * 1.5 for n in a3_here]
            nx.draw_networkx_nodes(
                G_sub, pos_sub, nodelist=a3_here, ax=ax,
                node_size=a3_s, node_color=COLOR_A3,
                edgecolors="#000000", linewidths=2.0)

        z_labels = {}; z_fs = {}; z_fc = {}
        if len(c_nodes) <= 40:
            for n in G_sub.nodes():
                z_labels[n] = sym(n)
                z_fs[n] = 15 + 21 * get_intra(n)
                z_fc[n] = COLOR_A3_TEXT if n in a3_here else "#000000"
        else:
            scored = [(n, get_intra(n)) for n in G_sub.nodes()]
            scored.sort(key=lambda x: -x[1])
            for n, sc in scored[:15]:
                z_labels[n] = sym(n)
                z_fs[n] = 15 + 21 * get_intra(n)
                z_fc[n] = "#000000"
            for n in a3_here:
                z_labels[n] = sym(n)
                z_fs[n] = FONT_AXIS
                z_fc[n] = COLOR_A3_TEXT

        if z_labels:
            zl_pos = repel_labels(pos_sub, z_labels, data_range=3.0,
                                  n_iters=250, push_strength=0.07)
            draw_labels_with_boxes(ax, pos_sub, zl_pos, z_labels,
                                   z_fs, z_fc, char_w=0.04)

        a3_note = ""
        a3_sym = [sym(g) for g in a3_here]
        if a3_sym:
            a3_note = f" -- {', '.join(a3_sym)}"
        sat_tag = " [satellite]" if c in sat_comms else ""

        ax.set_title(f"Community {c} -- {len(c_nodes)} genes, "
                     f"{G_sub.number_of_edges()} edges{a3_note}{sat_tag}",
                     fontsize=FONT_AXIS, pad=10)
        ax.axis("off")
        plt.tight_layout()

        for ext in ["png", "pdf"]:
            plt.savefig(os.path.join(zoom_dir,
                        f"{cancer_type}_community_{c:02d}.{ext}"),
                        dpi=300, bbox_inches="tight")
        plt.close()
        log(f"  [SAVE] Community {c} ({len(c_nodes)} genes)")

    # =================================================================
    # SUPPLEMENT -- Selection Plot
    # =================================================================
    banner("[SUPPLEMENT] SBS2 vs A3 sum selection plot")

    a3_median = float(all_hnsc["A3_sum"].median())
    sbs2_low_max = float(group_low["SBS2"].max())
    sbs2_high_min = float(group_high["SBS2"].min())
    x_max_data = all_hnsc["A3_sum"].max()
    y_max_data = all_hnsc["SBS2"].max()

    a3_sorted = all_hnsc["A3_sum"].sort_values(ascending=False).values
    has_outliers = False
    if len(a3_sorted) >= 3:
        best_gap_ratio = 1.0
        best_gap_idx = -1
        for i in range(min(10, len(a3_sorted) - 1)):
            ratio = (a3_sorted[i] / a3_sorted[i + 1]
                     if a3_sorted[i + 1] > 0 else 1.0)
            if ratio > best_gap_ratio:
                best_gap_ratio = ratio
                best_gap_idx = i
        if best_gap_ratio > 1.5 and best_gap_idx >= 0:
            cluster_max = a3_sorted[best_gap_idx + 1]
            outlier_min = a3_sorted[best_gap_idx]
            x_break_start = cluster_max * 1.10
            x_break_end = outlier_min * 0.90
            has_outliers = True

    if has_outliers:
        fig, (ax_main, ax_break) = plt.subplots(
            1, 2, sharey=True, figsize=(16, 10),
            gridspec_kw={"width_ratios": [3, 1], "wspace": 0.03})
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

        ax.axvline(a3_median, ls="--", c="#4D4D4D", lw=0.8, zorder=1)
        ax.axhline(sbs2_high_min, ls=":", c="#4D4D4D", lw=0.8, alpha=0.6)
        ax.axhline(sbs2_low_max, ls=":", c="#4D4D4D", lw=0.8, alpha=0.6)

        ax.scatter(all_hnsc["A3_sum"], all_hnsc["SBS2"],
                   facecolors=COLOR_GRAY, edgecolors="#000000",
                   linewidths=0.3, s=100, alpha=0.5, zorder=2)
        ax.scatter(group_high["A3_sum"], group_high["SBS2"],
                   facecolors=COLOR_HIGH, edgecolors="#000000",
                   linewidths=0.8, s=100, alpha=0.85, zorder=3,
                   label=f"SBS2-HIGH (n={len(group_high)})"
                   if ax == ax_main else "")
        ax.scatter(group_low["A3_sum"], group_low["SBS2"],
                   facecolors=COLOR_LOW, edgecolors="#000000",
                   linewidths=0.8, s=100, alpha=0.85, zorder=3,
                   label=f"SBS2-LOW (n={len(group_low)})"
                   if ax == ax_main else "")

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(-y_pad * 0.5, y_max)
        ax.tick_params(axis="both", labelsize=FONT_TICK)

    ax_main.set_xlabel("A3A + A3B Expression (FPKM-UQ)",
                       fontsize=FONT_AXIS, labelpad=12)
    ax_main.set_ylabel("SBS2 Mutational Weight",
                       fontsize=FONT_AXIS, labelpad=12)
    ax_main.legend(fontsize=FONT_LEGEND - 4, loc="upper left",
                   framealpha=0.9)

    if ax_break is not None:
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
                    f"{cancer_type}_Supplement_selection.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"[SAVE] Selection plot -> {ct_fig_dir}")

    # =================================================================
    banner(f"FIGURE 2 COMPLETE -- {cancer_type}")
    print(f"\n  Output: {ct_fig_dir}")
    print(f"  Panel A: Methods overview (HIGH - LOW = DIFF -> network)")
    print(f"  Panel B: Full network (two-tier layout)")
    print(f"  Panel B2: LCC-only network (backup)")
    print(f"  Panel C: Community zoom insets")
    print(f"  Supplement: Selection plot")

banner("ALL PANELS COMPLETE")
