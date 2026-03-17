#!/usr/bin/env python3
"""
Generate_Figure2_Panels.py

Generate publication-quality Figure 2 panels from pipeline outputs.
Reads from Steps 03–06 outputs and produces three panels:

  Panel A — SBS2 vs A3 sum selection plot (Figure 1 style)
  Panel B — Community-ordered correlation heatmap (DIFF matrix)
  Panel C — Full community network + per-community zoomed subplots

Style matches Figure 1: colored background regions, broken x-axis,
large fonts (26pt axes), consistent color palette.

Run AFTER the full pipeline (Steps 01–06) has completed.

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
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
import seaborn as sns

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

# Output directory
FIG_DIR = ensure_dir(os.path.join(FIG2_ROOT, "FIGURE_2_PANELS"))

# =============================================================================
# LOOP CANCER TYPES
# =============================================================================
for cancer_type in CANCER_TYPES:

    banner(f"FIGURE 2 — {cancer_type}", char="=")

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
    high_path = os.path.join(DIR_03_DIFFEXPR, cancer_type, f"{cancer_type}_SBS2_HIGH_group.pkl")
    low_path  = os.path.join(DIR_03_DIFFEXPR, cancer_type, f"{cancer_type}_SBS2_LOW_group.pkl")
    if not os.path.exists(high_path):
        log(f"[SKIP] Groups not found")
        continue
    group_high = pd.read_pickle(high_path)
    group_low  = pd.read_pickle(low_path)
    log(f"HIGH: {len(group_high)}, LOW: {len(group_low)}")

    # ---- Load correlation DIFF matrix
    corr_diff_path = os.path.join(DIR_04_NETWORKS, cancer_type, "corr_matrices",
                                   f"{cancer_type}_corr_DIFF.pkl")
    if os.path.exists(corr_diff_path):
        corr_diff = pd.read_pickle(corr_diff_path)
        log(f"DIFF corr matrix: {corr_diff.shape}")
    else:
        log("[SKIP] DIFF correlation matrix not found")
        continue

    # ---- Load community assignments
    part_path = os.path.join(DIR_05_COMMUNITIES, cancer_type, f"{cancer_type}_best_partition.csv")
    if not os.path.exists(part_path):
        log("[SKIP] Best partition not found")
        continue
    partition_df = pd.read_csv(part_path)
    gene_to_comm = dict(zip(partition_df["gene"], partition_df["community"]))
    log(f"Partition: {len(gene_to_comm)} genes, {partition_df['community'].nunique()} communities")

    # ---- Load community graph
    comm_graph_path = os.path.join(DIR_05_COMMUNITIES, cancer_type, f"{cancer_type}_G_comm.gpickle")
    if os.path.exists(comm_graph_path):
        with open(comm_graph_path, "rb") as f:
            G_comm = pickle.load(f)
        log(f"Community graph: {G_comm.number_of_nodes()} nodes, {G_comm.number_of_edges()} edges")
    else:
        G_comm = None

    # =====================================================================
    # PANEL A — SBS2 vs A3 Sum Selection Plot (Figure 1 style)
    # =====================================================================
    banner("[PANEL A] Selection plot")

    a3_median = float(all_hnsc["A3_sum"].median())

    # Determine SBS2 boundary (max SBS2 in LOW group)
    sbs2_low_max  = float(group_low["SBS2"].max())
    sbs2_high_min = float(group_high["SBS2"].min())

    x_max_data = all_hnsc["A3_sum"].max()
    y_max_data = all_hnsc["SBS2"].max()

    # Broken x-axis: main range and outlier range
    x_break_start = 400
    x_break_end   = 600
    has_outliers = x_max_data > x_break_start

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
            x_lo, x_hi = -2, min(x_break_start, x_max_data * 1.05)
        elif ax == ax_break:
            x_lo, x_hi = x_break_end, x_max_data * 1.05

        # ---- Background regions
        # Cream region: below SBS2 boundary, right of A3 median
        cream_poly = MplPolygon(
            [(a3_median, 0), (a3_median, sbs2_low_max), (x_hi, sbs2_low_max), (x_hi, 0)],
            closed=True, facecolor=COLOR_CREAM, alpha=0.75, edgecolor="none", zorder=0
        )
        ax.add_patch(cream_poly)

        # Coral region: above SBS2 HIGH boundary, right of A3 median
        coral_poly = MplPolygon(
            [(a3_median, sbs2_high_min), (a3_median, y_max), (x_hi, y_max), (x_hi, sbs2_high_min)],
            closed=True, facecolor=COLOR_CORAL, alpha=0.55, edgecolor="none", zorder=0
        )
        ax.add_patch(coral_poly)

        # Teal region: left of A3 median (low A3 — not selected for analysis)
        teal_poly = MplPolygon(
            [(x_lo, 0), (x_lo, y_max), (a3_median, y_max), (a3_median, 0)],
            closed=True, facecolor=COLOR_TEAL, alpha=0.55, edgecolor="none", zorder=0
        )
        ax.add_patch(teal_poly)

        # Middle zone (between low_max and high_min) — lighter cream
        if sbs2_high_min > sbs2_low_max:
            mid_poly = MplPolygon(
                [(a3_median, sbs2_low_max), (a3_median, sbs2_high_min),
                 (x_hi, sbs2_high_min), (x_hi, sbs2_low_max)],
                closed=True, facecolor="#e8e4c8", alpha=0.55, edgecolor="none", zorder=0
            )
            ax.add_patch(mid_poly)

        # ---- Boundary lines
        ax.axvline(a3_median, ls="--", c="grey30", lw=0.8, zorder=1)
        ax.axhline(sbs2_high_min, ls=":", c="grey30", lw=0.8, alpha=0.6, zorder=1)
        ax.axhline(sbs2_low_max, ls=":", c="grey30", lw=0.8, alpha=0.6, zorder=1)

        # ---- Data points
        # All tumors (gray)
        ax.scatter(all_hnsc["A3_sum"], all_hnsc["SBS2"],
                   facecolors=COLOR_GRAY, edgecolors="black", linewidths=0.3,
                   s=50, alpha=0.5, zorder=2)

        # HIGH group (coral)
        ax.scatter(group_high["A3_sum"], group_high["SBS2"],
                   facecolors=COLOR_HIGH, edgecolors="black", linewidths=0.8,
                   s=70, alpha=0.85, zorder=3,
                   label=f"SBS2-HIGH (n={len(group_high)})" if ax == ax_main else "")

        # LOW group (teal)
        ax.scatter(group_low["A3_sum"], group_low["SBS2"],
                   facecolors=COLOR_LOW, edgecolors="black", linewidths=0.8,
                   s=70, alpha=0.85, zorder=3,
                   label=f"SBS2-LOW (n={len(group_low)})" if ax == ax_main else "")

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(-y_pad * 0.5, y_max)

        ax.tick_params(axis="both", labelsize=22)

    # ---- Axis labels and styling
    ax_main.set_xlabel("A3A + A3B Expression (FPKM-UQ)", fontsize=26, labelpad=12)
    ax_main.set_ylabel("SBS2 Mutational Weight", fontsize=26, labelpad=12)
    ax_main.legend(fontsize=18, loc="upper left", framealpha=0.9)

    # Region count labels
    n_high_a3 = (all_hnsc["A3_sum"] >= a3_median).sum()
    n_low_a3  = (all_hnsc["A3_sum"] < a3_median).sum()
    ax_main.annotate(f"Low A3\nn = {n_low_a3}",
                     xy=(a3_median * 0.25, y_max_data * 0.85),
                     fontsize=18, fontweight="bold", color="grey30", ha="center")

    # Total count
    ax_main.annotate(f"(n = {len(all_hnsc)})",
                     xy=(0.98, 0.98), xycoords="axes fraction",
                     fontsize=20, ha="right", va="top")

    if ax_break is not None:
        ax_break.set_xlabel("")
        ax_break.tick_params(axis="y", left=False, labelleft=False)

        # Break marks
        d = 0.015
        for ax_b, side in [(ax_main, "right"), (ax_break, "left")]:
            kwargs = dict(transform=ax_b.transAxes, color="k", clip_on=False, lw=1.5)
            if side == "right":
                ax_b.plot((1 - d, 1 + d), (-d, +d), **kwargs)
                ax_b.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
                ax_b.spines["right"].set_visible(False)
            else:
                ax_b.plot((-d, +d), (-d, +d), **kwargs)
                ax_b.plot((-d, +d), (1 - d, 1 + d), **kwargs)
                ax_b.spines["left"].set_visible(False)

    plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_A_selection.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_A_selection.pdf"), bbox_inches="tight")
    plt.close()
    log(f"[SAVE] Panel A -> {ct_fig_dir}")

    # =====================================================================
    # PANEL B — Community-ordered correlation heatmap (DIFF matrix)
    # =====================================================================
    banner("[PANEL B] Community correlation heatmap")

    # Get genes that are in both the DIFF matrix and the partition
    comm_genes = [g for g in partition_df["gene"] if g in corr_diff.index]

    # Order genes by community, then by degree within community
    comm_to_genes = {}
    for g in comm_genes:
        c = gene_to_comm.get(g, -1)
        comm_to_genes.setdefault(c, []).append(g)

    # Sort communities by size (descending), put "Other"/largest-ID last
    comm_ids_sorted = sorted(comm_to_genes.keys(),
                             key=lambda c: len(comm_to_genes[c]), reverse=True)

    # Build ordered gene list
    ordered_genes = []
    comm_boundaries = []  # (start_idx, end_idx, comm_id) for annotation bars
    comm_labels = []

    for c in comm_ids_sorted:
        genes_in_comm = comm_to_genes[c]
        # Sort by degree within community if we have the graph
        if G_comm is not None:
            genes_in_comm = sorted(genes_in_comm,
                                   key=lambda g: G_comm.degree(g) if g in G_comm else 0,
                                   reverse=True)
        start = len(ordered_genes)
        ordered_genes.extend(genes_in_comm)
        end = len(ordered_genes)
        comm_boundaries.append((start, end, c))
        comm_labels.append(f"C{c}\n({len(genes_in_comm)})")

    # Subset and order the DIFF correlation matrix
    ordered_genes_present = [g for g in ordered_genes if g in corr_diff.index]
    heatmap_data = corr_diff.loc[ordered_genes_present, ordered_genes_present].values

    # Determine symmetric color scale
    vmax = np.nanpercentile(np.abs(heatmap_data), 99)
    if vmax == 0:
        vmax = 1.0

    # Plot
    fig, ax = plt.subplots(figsize=(16, 14))

    im = ax.imshow(heatmap_data, aspect="auto", interpolation="nearest",
                   cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    # Community boundary lines
    for start, end, c in comm_boundaries:
        if start > 0:
            ax.axhline(start - 0.5, color="black", lw=1.2, alpha=0.8)
            ax.axvline(start - 0.5, color="black", lw=1.2, alpha=0.8)

    # Community labels on the side
    for start, end, c in comm_boundaries:
        mid = (start + end) / 2
        n_genes = end - start
        if n_genes >= 10:  # Only label substantial communities
            ax.annotate(f"C{c} ({n_genes})", xy=(len(ordered_genes_present) + 5, mid),
                        fontsize=11, fontweight="bold", va="center", ha="left",
                        annotation_clip=False)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.12)
    cbar.set_label("Differential Correlation (TOP − BOTTOM)", fontsize=20, labelpad=12)
    cbar.ax.tick_params(labelsize=18)

    ax.set_xlabel("Genes (ordered by community)", fontsize=26, labelpad=12)
    ax.set_ylabel("Genes (ordered by community)", fontsize=26, labelpad=12)
    ax.set_title(f"{cancer_type} | DIFF Correlation — Community Structure", fontsize=22, pad=15)

    # Remove individual gene tick labels (too many)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_B_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_B_heatmap.pdf"), bbox_inches="tight")
    plt.close()
    log(f"[SAVE] Panel B -> {ct_fig_dir}")

    # =====================================================================
    # PANEL C — Full community network plot
    # =====================================================================
    banner("[PANEL C] Community network plot")

    if G_comm is not None and G_comm.number_of_nodes() > 0:
        # Community-aware layout
        comm_to_nodelist = {}
        for n in G_comm.nodes():
            c = gene_to_comm.get(n, -1)
            comm_to_nodelist.setdefault(c, []).append(n)

        # Community color map
        unique_comms = sorted(comm_to_nodelist.keys())
        n_comms = len(unique_comms)
        cmap_comms = plt.cm.get_cmap("tab20", max(n_comms, 1))
        comm_color_map = {c: cmap_comms(i) for i, c in enumerate(unique_comms)}

        # Spring layout with community seeding
        pos = nx.spring_layout(G_comm, seed=COMMUNITY_BASE_SEED,
                               weight="abs_weight",
                               k=3.0 / np.sqrt(max(G_comm.number_of_nodes(), 1)),
                               iterations=300)

        # ---- Full network plot
        fig, ax = plt.subplots(figsize=(18, 16))

        # Draw edges (light)
        nx.draw_networkx_edges(G_comm, pos, ax=ax, alpha=0.05, width=0.2,
                               edge_color="gray")

        # Draw nodes colored by community
        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) == 0:
                continue
            nx.draw_networkx_nodes(G_comm, pos, nodelist=c_nodes, ax=ax,
                                   node_size=25, node_color=[comm_color_map[c]],
                                   alpha=0.75, linewidths=0.3, edgecolors="black",
                                   label=f"C{c} (n={len(c_nodes)})")

        # Highlight A3 genes
        a3_in_graph = [g for g in A3_GENES if g in G_comm.nodes()]
        if a3_in_graph:
            nx.draw_networkx_nodes(G_comm, pos, nodelist=a3_in_graph, ax=ax,
                                   node_size=200, node_color="gold",
                                   edgecolors="black", linewidths=2.0, zorder=5)
            labels_a3 = {n: sym(n) for n in a3_in_graph}
            nx.draw_networkx_labels(G_comm, pos, labels=labels_a3, ax=ax,
                                    font_size=10, font_weight="bold", font_color="darkred")

        # Highlight top 3 hub genes per community (by degree)
        hub_labels = {}
        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 5:
                continue
            deg = {n: G_comm.degree(n) for n in c_nodes}
            top3 = sorted(deg.items(), key=lambda x: x[1], reverse=True)[:3]
            for g, d in top3:
                if g not in a3_in_graph:
                    hub_labels[g] = sym(g)

        if hub_labels:
            nx.draw_networkx_labels(G_comm, pos, labels=hub_labels, ax=ax,
                                    font_size=7, font_weight="bold", font_color="black",
                                    alpha=0.8)

        ax.set_title(f"{cancer_type} | DIFF Network — Leiden Communities",
                     fontsize=22, pad=15)
        ax.legend(fontsize=9, loc="upper right", ncol=2, framealpha=0.9,
                  markerscale=2)
        ax.axis("off")
        plt.tight_layout()

        plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_C_network_full.png"),
                    dpi=300, bbox_inches="tight")
        plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_C_network_full.pdf"),
                    bbox_inches="tight")
        plt.close()
        log(f"[SAVE] Panel C (full) -> {ct_fig_dir}")

        # ---- Individual community zoomed subplots
        banner("[PANEL C] Per-community zoomed plots")

        zoom_dir = ensure_dir(os.path.join(ct_fig_dir, "community_zooms"))

        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 5:
                continue

            # Extract subgraph
            G_sub = G_comm.subgraph(c_nodes).copy()
            if G_sub.number_of_edges() == 0:
                continue

            # Layout for this community
            pos_sub = nx.spring_layout(G_sub, seed=COMMUNITY_BASE_SEED,
                                       weight="abs_weight",
                                       k=2.0 / np.sqrt(max(G_sub.number_of_nodes(), 1)),
                                       iterations=200)

            fig, ax = plt.subplots(figsize=(12, 10))

            # Edge coloring by weight sign
            edge_colors = []
            edge_widths = []
            for u, v, d in G_sub.edges(data=True):
                w = d.get("weight", 0)
                edge_colors.append("firebrick" if w > 0 else "steelblue")
                edge_widths.append(0.5 + 2.0 * abs(w))

            nx.draw_networkx_edges(G_sub, pos_sub, ax=ax,
                                   edge_color=edge_colors, width=edge_widths,
                                   alpha=0.4)

            # Node sizing by degree
            degrees = dict(G_sub.degree())
            max_deg = max(degrees.values()) if degrees else 1
            node_sizes = [40 + 300 * (degrees[n] / max_deg) for n in G_sub.nodes()]

            nx.draw_networkx_nodes(G_sub, pos_sub, ax=ax,
                                   node_size=node_sizes,
                                   node_color=[comm_color_map[c]],
                                   alpha=0.85, edgecolors="black", linewidths=0.8)

            # Label all genes if small community, top genes if large
            if len(c_nodes) <= 40:
                labels = {n: sym(n) for n in G_sub.nodes()}
                font_size = max(6, 10 - len(c_nodes) // 10)
            else:
                # Top 15 by degree
                top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:15]
                labels = {n: sym(n) for n, d in top_nodes}
                font_size = 8

            # Highlight A3 genes in this community
            a3_here = [g for g in A3_GENES if g in G_sub.nodes()]
            if a3_here:
                nx.draw_networkx_nodes(G_sub, pos_sub, nodelist=a3_here, ax=ax,
                                       node_size=250, node_color="gold",
                                       edgecolors="black", linewidths=2.0, zorder=5)
                for g in a3_here:
                    labels[g] = sym(g)

            nx.draw_networkx_labels(G_sub, pos_sub, labels=labels, ax=ax,
                                    font_size=font_size, font_weight="bold")

            ax.set_title(f"Community {c} — {len(c_nodes)} genes, "
                         f"{G_sub.number_of_edges()} edges",
                         fontsize=18, pad=10)
            ax.axis("off")
            plt.tight_layout()

            zoom_path = os.path.join(zoom_dir, f"{cancer_type}_community_{c:02d}.png")
            plt.savefig(zoom_path, dpi=300, bbox_inches="tight")
            plt.close()

            log(f"  [SAVE] Community {c} ({len(c_nodes)} genes) -> {zoom_path}")

        log(f"[SAVE] All community zooms -> {zoom_dir}")

    else:
        log("[SKIP] No community graph available for Panel C")

    # =====================================================================
    # SUMMARY
    # =====================================================================
    banner(f"FIGURE 2 COMPLETE — {cancer_type}")
    print(f"\n  Output directory: {ct_fig_dir}")
    print(f"  Panel A: Selection plot (PDF + PNG)")
    print(f"  Panel B: Community heatmap (PDF + PNG)")
    print(f"  Panel C: Full network + per-community zooms")

banner("ALL PANELS COMPLETE")
