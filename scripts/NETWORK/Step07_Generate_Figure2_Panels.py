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

    # Load TOP correlation matrix (for Panel B heatmap)
    corr_top_path = os.path.join(DIR_04_NETWORKS, cancer_type, "corr_matrices",
                                  f"{cancer_type}_corr_TOP.pkl")
    if os.path.exists(corr_top_path):
        corr_top = pd.read_pickle(corr_top_path)
        log(f"TOP corr matrix: {corr_top.shape}")
    else:
        corr_top = None
        log("[WARNING] TOP correlation matrix not found — Panel B will use DIFF")

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

    # ---- Data-driven x-axis break
    # Scan from the highest A3 values downward to find the largest gap.
    # The break goes in the biggest gap between consecutive ranked tumors.
    a3_sorted = all_hnsc["A3_sum"].sort_values(ascending=False).values
    log(f"[PANEL A] Top 5 A3_sum values: {a3_sorted[:5].round(2).tolist()}")

    has_outliers = False
    if len(a3_sorted) >= 3:
        # Check gaps between consecutive top values (top1-top2, top2-top3, etc.)
        # Find the largest ratio gap in the top 10 values
        best_gap_ratio = 1.0
        best_gap_idx = -1
        n_check = min(10, len(a3_sorted) - 1)

        for i in range(n_check):
            higher = a3_sorted[i]
            lower  = a3_sorted[i + 1]
            ratio  = higher / lower if lower > 0 else 1.0
            if ratio > best_gap_ratio:
                best_gap_ratio = ratio
                best_gap_idx = i

        log(f"[PANEL A] Largest gap: position {best_gap_idx} -> {best_gap_idx+1}, "
            f"values {a3_sorted[best_gap_idx]:.2f} vs {a3_sorted[best_gap_idx+1]:.2f}, "
            f"ratio {best_gap_ratio:.2f}")

        if best_gap_ratio > 1.5 and best_gap_idx >= 0:
            outlier_min = a3_sorted[best_gap_idx]       # lowest outlier value
            cluster_max = a3_sorted[best_gap_idx + 1]   # highest main-cluster value
            x_break_start = cluster_max * 1.10   # 10% past main cluster
            x_break_end   = outlier_min * 0.90   # 10% before outlier
            has_outliers = True
            n_outliers = best_gap_idx + 1
            log(f"[PANEL A] {n_outliers} outlier(s) above gap")
            log(f"[PANEL A] Breaking x-axis: {x_break_start:.1f} — {x_break_end:.1f}")
        else:
            log("[PANEL A] No substantial gap (ratio < 1.5) — no x-axis break")
    else:
        log("[PANEL A] Too few tumors for break analysis")

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
            if has_outliers:
                x_lo, x_hi = -2, x_break_start
            else:
                x_lo, x_hi = -2, x_max_data * 1.05
        elif ax == ax_break:
            x_lo, x_hi = x_break_end, x_max_data * 1.05

        # ---- Background regions (only right of A3 median — the analysis zone)
        # Cream region: below SBS2 LOW boundary, right of A3 median
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

        # Middle zone (between LOW max and HIGH min) — lighter cream
        if sbs2_high_min > sbs2_low_max:
            mid_poly = MplPolygon(
                [(a3_median, sbs2_low_max), (a3_median, sbs2_high_min),
                 (x_hi, sbs2_high_min), (x_hi, sbs2_low_max)],
                closed=True, facecolor="#e8e4c8", alpha=0.55, edgecolor="none", zorder=0
            )
            ax.add_patch(mid_poly)

        # ---- Boundary lines
        ax.axvline(a3_median, ls="--", c="#4D4D4D", lw=0.8, zorder=1)
        ax.axhline(sbs2_high_min, ls=":", c="#4D4D4D", lw=0.8, alpha=0.6, zorder=1)
        ax.axhline(sbs2_low_max, ls=":", c="#4D4D4D", lw=0.8, alpha=0.6, zorder=1)

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

    # Region label — mark the analysis zone
    ax_main.annotate("High A3\n(analysis zone)",
                     xy=(a3_median * 1.3, y_max_data * 0.85),
                     fontsize=16, fontweight="bold", color="#4D4D4D", ha="left",
                     fontstyle="italic")

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
    # PANEL B — Dual heatmap: TOP + DIFF with hierarchical clustering
    # =====================================================================
    # Shows BOTH the TOP (high-SBS2) and DIFF (TOP-BOTTOM) correlation
    # matrices, ordered by community with hierarchical clustering applied:
    #   1. Within each community: genes are re-ordered by dendrogram
    #   2. Between communities: communities are ordered by linkage
    # This reveals co-expression structure within communities and
    # similarity between communities.
    # =====================================================================

    banner("[PANEL B] Dual heatmap: TOP + DIFF with hierarchical clustering")
 
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform, pdist
 
    # Select matrices to plot
    matrices_to_plot = []
    if corr_top is not None:
        matrices_to_plot.append(("TOP (High SBS2)", corr_top, "plasma"))
    if corr_diff is not None:
        matrices_to_plot.append(("DIFF (TOP − BOTTOM)", corr_diff, "plasma"))
    elif corr_top is None:
        log("[SKIP] No correlation matrices available for Panel B")
 
    if len(matrices_to_plot) > 0:
 
        # ---- Build hierarchically-clustered gene order ----
 
        # Get genes present in both the partition and the matrices
        ref_matrix = matrices_to_plot[0][1]
        comm_genes = [g for g in partition_df["gene"] if g in ref_matrix.index]
 
        # Group genes by community
        comm_to_genes = {}
        for g in comm_genes:
            c = gene_to_comm.get(g, -1)
            comm_to_genes.setdefault(c, []).append(g)
 
        # --- Step 1: Hierarchical clustering WITHIN each community ---
        clustered_comm_to_genes = {}
        for c, genes in comm_to_genes.items():
            if len(genes) <= 2:
                clustered_comm_to_genes[c] = genes
                continue
 
            # Use the TOP matrix (or DIFF if TOP unavailable) for clustering
            sub_matrix = ref_matrix.loc[genes, genes].values
            # Replace NaN with 0
            sub_matrix = np.nan_to_num(sub_matrix, nan=0.0)
 
            try:
                # Convert correlation to distance: 1 - |corr|
                dist_matrix = 1.0 - np.abs(sub_matrix)
                np.fill_diagonal(dist_matrix, 0)
                # Ensure symmetry
                dist_matrix = (dist_matrix + dist_matrix.T) / 2.0
                # Clip to avoid negative distances from floating point
                dist_matrix = np.clip(dist_matrix, 0, 2)
 
                condensed = squareform(dist_matrix, checks=False)
                Z = linkage(condensed, method="average")
                order = leaves_list(Z)
                clustered_comm_to_genes[c] = [genes[i] for i in order]
            except Exception as e:
                log(f"  [WARNING] Clustering failed for community {c}: {e}")
                clustered_comm_to_genes[c] = genes
 
        # --- Step 2: Hierarchical clustering BETWEEN communities ---
        comm_ids = list(clustered_comm_to_genes.keys())
 
        if len(comm_ids) > 2:
            # Compute community-level distance from mean inter-community correlation
            n_comms = len(comm_ids)
            comm_dist = np.zeros((n_comms, n_comms))
 
            for i in range(n_comms):
                for j in range(i + 1, n_comms):
                    genes_i = clustered_comm_to_genes[comm_ids[i]]
                    genes_j = clustered_comm_to_genes[comm_ids[j]]
                    # Mean absolute correlation between communities
                    block = ref_matrix.loc[genes_i, genes_j].values
                    mean_corr = np.nanmean(np.abs(block))
                    dist = 1.0 - mean_corr
                    comm_dist[i, j] = dist
                    comm_dist[j, i] = dist
 
            condensed_comm = squareform(comm_dist, checks=False)
            Z_comm = linkage(condensed_comm, method="average")
            comm_order = leaves_list(Z_comm)
            comm_ids_sorted = [comm_ids[i] for i in comm_order]
        else:
            # Sort by size (descending) if too few communities
            comm_ids_sorted = sorted(comm_ids,
                                     key=lambda c: len(clustered_comm_to_genes[c]),
                                     reverse=True)
 
        # --- Step 3: Build final ordered gene list ---
        ordered_genes = []
        comm_boundaries = []
        for c in comm_ids_sorted:
            genes_in_comm = clustered_comm_to_genes[c]
            start = len(ordered_genes)
            ordered_genes.extend(genes_in_comm)
            end = len(ordered_genes)
            comm_boundaries.append((start, end, c))
 
        log(f"  Ordered {len(ordered_genes)} genes across {len(comm_ids_sorted)} communities")
 
        # --- Step 4: Plot each matrix ---
        n_panels = len(matrices_to_plot)
        fig, axes = plt.subplots(1, n_panels, figsize=(16 * n_panels, 14))
        if n_panels == 1:
            axes = [axes]
 
        for ax, (label, matrix, cmap_name) in zip(axes, matrices_to_plot):
            # Subset and order
            genes_present = [g for g in ordered_genes if g in matrix.index]
            heatmap_data = matrix.loc[genes_present, genes_present].values
 
            vmin = np.nanpercentile(heatmap_data, 1)
            vmax = np.nanpercentile(heatmap_data, 99)
            if vmax <= vmin:
                vmax = 1.0
 
            im = ax.imshow(heatmap_data, aspect="auto", interpolation="nearest",
                           cmap=cmap_name, vmin=vmin, vmax=vmax)
 
            # Community boundary lines
            for start, end, c in comm_boundaries:
                if start > 0:
                    ax.axhline(start - 0.5, color="white", lw=2, alpha=0.8)
                    ax.axvline(start - 0.5, color="white", lw=2, alpha=0.8)
 
            # Community labels on the right side
            for start, end, c in comm_boundaries:
                mid = (start + end) / 2
                n_genes = end - start
                if n_genes >= 5:
                    ax.annotate(f"C{c} ({n_genes})", xy=(len(genes_present) + 5, mid),
                                fontsize=20, fontweight="bold", va="center", ha="left",
                                annotation_clip=False)
 
            cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.12)
            cbar.set_label(f"Spearman rho ({label})", fontsize=24, labelpad=12)
            cbar.ax.tick_params(labelsize=22)
 
            ax.set_xlabel("Genes (clustered within communities)", fontsize=28, labelpad=12)
            ax.set_ylabel("Genes (clustered within communities)", fontsize=28, labelpad=12)
            ax.set_title(f"{cancer_type} | {label}", fontsize=28, pad=15)
            ax.set_xticks([])
            ax.set_yticks([])
 
        plt.tight_layout(w_pad=4)
        plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_B_dual_heatmap.png"),
                    dpi=300, bbox_inches="tight")
        plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_B_dual_heatmap.pdf"),
                    bbox_inches="tight")
        plt.close()
        log(f"[SAVE] Panel B (dual heatmap) -> {ct_fig_dir}")
 
    # =====================================================================
    # PANEL C — Full community network plot (exploded layout)
    # =====================================================================
    banner("[PANEL C] Community network plot (exploded layout)")
 
    if G_comm is not None and G_comm.number_of_nodes() > 0:
 
        # ---- Group nodes by community ----
        comm_to_nodelist = {}
        for n in G_comm.nodes():
            c = gene_to_comm.get(n, -1)
            comm_to_nodelist.setdefault(c, []).append(n)
 
        unique_comms = sorted(comm_to_nodelist.keys())
        n_comms = len(unique_comms)
        cmap_comms = plt.cm.get_cmap("tab20", max(n_comms, 1))
        comm_color_map = {c: cmap_comms(i) for i, c in enumerate(unique_comms)}
 
        # ---- Community-aware exploded layout ----
        # Step 1: Build a super-graph where each community is a node.
        #         Edge weight = sum of |edge weights| between communities.
        # Step 2: Layout the super-graph with large spacing.
        # Step 3: Layout nodes within each community (local spring).
        # Step 4: Combine: pos[node] = inter_scale * center + intra_scale * local
 
        INTER_SCALE = 6.0    # How far apart communities are placed
        INTRA_SCALE = 1.0    # How spread out nodes are within a community
        K_LOCAL = None        # spring k for within-community (auto)
        ITERS_LOCAL = 200
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
 
        # Layout super-graph
        if CG.number_of_edges() > 0:
            pos_comm = nx.spring_layout(CG, seed=COMMUNITY_BASE_SEED,
                                        weight="weight",
                                        k=2.0 / np.sqrt(max(CG.number_of_nodes(), 1)),
                                        iterations=ITERS_GLOBAL)
        else:
            # No inter-community edges: arrange in a circle
            pos_comm = nx.circular_layout(CG)
 
        # Scale community centers
        pos_comm = {c: INTER_SCALE * np.asarray(p, dtype=float)
                    for c, p in pos_comm.items()}
 
        # Layout nodes within each community
        pos = {}
        for c in unique_comms:
            c_nodes = comm_to_nodelist[c]
            sub = G_comm.subgraph(c_nodes).copy()
 
            if len(c_nodes) == 1:
                pos[c_nodes[0]] = pos_comm.get(c, np.zeros(2))
                continue
 
            if sub.number_of_edges() > 0:
                k_val = K_LOCAL or (2.0 / np.sqrt(max(sub.number_of_nodes(), 1)))
                pos_sub = nx.spring_layout(sub, seed=COMMUNITY_BASE_SEED,
                                           weight="abs_weight",
                                           k=k_val, iterations=ITERS_LOCAL)
            else:
                # No internal edges: small jitter cloud
                pos_sub = {n: rng.normal(0, 0.05, size=2) for n in c_nodes}
 
            # Center and scale local coordinates
            coords = np.array([pos_sub[n] for n in c_nodes], dtype=float)
            coords = coords - coords.mean(axis=0, keepdims=True)
            coords = INTRA_SCALE * coords
 
            center = pos_comm.get(c, np.zeros(2, dtype=float))
            for i, n in enumerate(c_nodes):
                pos[n] = center + coords[i]
 
        log(f"  Exploded layout: INTER_SCALE={INTER_SCALE}, INTRA_SCALE={INTRA_SCALE}")
        log(f"  {len(pos)} nodes positioned across {n_comms} communities")
 
        # ---- Full network plot ----
        fig, ax = plt.subplots(figsize=(22, 20))
 
        # Draw edges — VISIBLE (increased alpha and width)
        edges = list(G_comm.edges(data=True))
        if edges:
            # Separate intra-community and inter-community edges
            intra_edges = []
            inter_edges = []
            for u, v, d in edges:
                cu = gene_to_comm.get(u, -1)
                cv = gene_to_comm.get(v, -1)
                if cu == cv:
                    intra_edges.append((u, v, d))
                else:
                    inter_edges.append((u, v, d))
 
            # Intra-community edges: more visible
            if intra_edges:
                intra_colors = []
                intra_widths = []
                for u, v, d in intra_edges:
                    w = d.get("weight", 0)
                    intra_colors.append("firebrick" if w > 0 else "steelblue")
                    intra_widths.append(0.3 + 1.5 * abs(w))
 
                nx.draw_networkx_edges(G_comm, pos,
                                       edgelist=[(u, v) for u, v, d in intra_edges],
                                       ax=ax, alpha=0.25, width=intra_widths,
                                       edge_color=intra_colors)
 
            # Inter-community edges: lighter
            if inter_edges:
                nx.draw_networkx_edges(G_comm, pos,
                                       edgelist=[(u, v) for u, v, d in inter_edges],
                                       ax=ax, alpha=0.08, width=0.3,
                                       edge_color="gray", style="dashed")
 
        # Draw nodes colored by community
        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) == 0:
                continue
 
            # Node sizing by degree
            degrees = {n: G_comm.degree(n) for n in c_nodes}
            max_deg = max(degrees.values()) if degrees else 1
            node_sizes = [30 + 200 * (degrees[n] / max_deg) for n in c_nodes]
 
            nx.draw_networkx_nodes(G_comm, pos, nodelist=c_nodes, ax=ax,
                                   node_size=node_sizes,
                                   node_color=[comm_color_map[c]],
                                   alpha=0.80, linewidths=0.4, edgecolors="black",
                                   label=f"C{c} (n={len(c_nodes)})")
 
        # Highlight A3 genes
        a3_in_graph = [g for g in A3_GENES if g in G_comm.nodes()]
        if a3_in_graph:
            nx.draw_networkx_nodes(G_comm, pos, nodelist=a3_in_graph, ax=ax,
                                   node_size=300, node_color="gold",
                                   edgecolors="black", linewidths=2.5)
            labels_a3 = {n: sym(n) for n in a3_in_graph}
            nx.draw_networkx_labels(G_comm, pos, labels=labels_a3, ax=ax,
                                    font_size=20, font_weight="bold",
                                    font_color="darkred")
 
        # Highlight top 3 hub genes per community (by degree) — DOUBLED text size
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
                                    font_size=14, font_weight="bold",
                                    font_color="black", alpha=0.9)
 
        ax.set_title(f"{cancer_type} | DIFF Network — Leiden Communities (Exploded Layout)",
                     fontsize=26, pad=15)
        ax.legend(fontsize=18, loc="upper right", ncol=2,
                  framealpha=0.9, markerscale=2,
                  title="Communities", title_fontsize=20)
        ax.axis("off")
        plt.tight_layout()
 
        plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_C_network_full.png"),
                    dpi=300, bbox_inches="tight")
        plt.savefig(os.path.join(ct_fig_dir, f"{cancer_type}_Panel_C_network_full.pdf"),
                    bbox_inches="tight")
        plt.close()
        log(f"[SAVE] Panel C (exploded layout) -> {ct_fig_dir}")
 
        # ---- Individual community zoomed subplots (with updated text sizes) ----
        banner("[PANEL C] Per-community zoomed plots")
 
        zoom_dir = ensure_dir(os.path.join(ct_fig_dir, "community_zooms"))
 
        for c in unique_comms:
            c_nodes = comm_to_nodelist.get(c, [])
            if len(c_nodes) < 5:
                continue
 
            G_sub = G_comm.subgraph(c_nodes).copy()
            if G_sub.number_of_edges() == 0:
                continue
 
            pos_sub = nx.spring_layout(G_sub, seed=COMMUNITY_BASE_SEED,
                                       weight="abs_weight",
                                       k=2.0 / np.sqrt(max(G_sub.number_of_nodes(), 1)),
                                       iterations=200)
 
            fig, ax = plt.subplots(figsize=(14, 12))
 
            # Edges — visible
            edge_colors = []
            edge_widths = []
            for u, v, d in G_sub.edges(data=True):
                w = d.get("weight", 0)
                edge_colors.append("firebrick" if w > 0 else "steelblue")
                edge_widths.append(0.5 + 2.0 * abs(w))
 
            nx.draw_networkx_edges(G_sub, pos_sub, ax=ax,
                                   edge_color=edge_colors, width=edge_widths,
                                   alpha=0.35)
 
            # Nodes
            degrees = dict(G_sub.degree())
            max_deg = max(degrees.values()) if degrees else 1
            node_sizes = [60 + 400 * (degrees[n] / max_deg) for n in G_sub.nodes()]
 
            nx.draw_networkx_nodes(G_sub, pos_sub, ax=ax,
                                   node_size=node_sizes,
                                   node_color=[comm_color_map[c]],
                                   alpha=0.85, edgecolors="black", linewidths=0.8)
 
            # Labels — doubled font sizes
            if len(c_nodes) <= 40:
                labels = {n: sym(n) for n in G_sub.nodes()}
                font_size = max(10, 16 - len(c_nodes) // 8)
            else:
                top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:15]
                labels = {n: sym(n) for n, d in top_nodes}
                font_size = 14
 
            # A3 genes
            a3_here = [g for g in A3_GENES if g in G_sub.nodes()]
            if a3_here:
                nx.draw_networkx_nodes(G_sub, pos_sub, nodelist=a3_here, ax=ax,
                                       node_size=350, node_color="gold",
                                       edgecolors="black", linewidths=2.0)
                for g in a3_here:
                    labels[g] = sym(g)
 
            nx.draw_networkx_labels(G_sub, pos_sub, labels=labels, ax=ax,
                                    font_size=font_size, font_weight="bold")
 
            ax.set_title(f"Community {c} — {len(c_nodes)} genes, "
                         f"{G_sub.number_of_edges()} edges",
                         fontsize=22, pad=10)
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
