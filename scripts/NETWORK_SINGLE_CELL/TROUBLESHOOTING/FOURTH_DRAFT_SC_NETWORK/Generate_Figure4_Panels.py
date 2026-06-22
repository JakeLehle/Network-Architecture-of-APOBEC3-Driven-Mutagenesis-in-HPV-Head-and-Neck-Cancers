#!/usr/bin/env python3
"""
Generate_Figure4_Panels.py
===========================

Generate publication-quality Figure 4 panels for the three-network
single-cell analysis.

Layout concept:
  - Central panel: Three-population UMAP (SBS2-HIGH / CNV-HIGH / NORMAL)
  - Two A3 community zoom panels (SBS2_VS_NORMAL, CNV_VS_NORMAL)
  - SBS2_VS_CNV parked for separate diagnostic work
  - Supplemental: Full network plots

Node coloring (A3 zoom panels):
  - Up in tumor (positive log2FC):  #fcad61
  - Up in normal (negative log2FC): #aadda4
  - APOBEC3 family genes:           #ed6a5a (coral)
  - Harris A3 interactor ring:      #F6D155 (mustard)

Edge coloring:
  - Gained co-expression (pos DIFF): #b22222 (firebrick)
  - Lost co-expression (neg DIFF):   #4682b4 (steelblue)

Style: Inherits Figure 2 conventions (28-34pt, hex codes, PDF+PNG 300 DPI).

Usage:
    conda run -n NETWORK python Generate_Figure4_Panels.py
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
FIG_DIR = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS")
os.makedirs(FIG_DIR, exist_ok=True)

ADATA_PATH = os.path.join(FIG4_ROOT, "00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(FIG4_ROOT, "01_group_selection/three_group_assignments.tsv")
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")

NETWORKS = [
    {
        "name": "SBS2_VS_NORMAL",
        "label": "SBS2-HIGH vs NORMAL",
        "short": "Mutagenic entry",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL"),
    },
    {
        "name": "CNV_VS_NORMAL",
        "label": "CNV-HIGH vs NORMAL",
        "short": "Productive infection",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL"),
    },
]

# SBS2_VS_CNV parked for now (separate diagnostic needed)
NETWORK_SBS2_VS_CNV = {
    "name": "SBS2_VS_CNV",
    "label": "SBS2-HIGH vs CNV-HIGH",
    "short": "Divergent fates",
    "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_CNV"),
}

# A3 gene identifiers (symbol format for SC)
A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}

# =============================================================================
# STYLE
# =============================================================================

# UMAP population colors
COLOR_SBS2_HIGH = "#ed6a5a"   # coral
COLOR_CNV_HIGH  = "#F6D155"   # mustard yellow
COLOR_NORMAL    = "#4682b4"   # steelblue

# Network node colors (DE direction)
COLOR_UP_TUMOR  = "#fcad61"   # up in tumor (positive log2FC)
COLOR_UP_NORMAL = "#aadda4"   # up in normal (negative log2FC)
COLOR_A3        = "#ed6a5a"   # APOBEC3 family genes
COLOR_A3_TEXT   = "#c0392b"
COLOR_HARRIS    = "#F6D155"   # Harris interactor ring

# Network edge colors
COLOR_EDGE_POS  = "#b22222"   # firebrick (gained co-expression)
COLOR_EDGE_NEG  = "#4682b4"   # steelblue (lost co-expression)

# Satellite/supplemental
COLOR_SAT       = "#e8913a"   # orange for satellites

FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22
LABEL_BG_ALPHA = 0.42

COMMUNITY_BASE_SEED = 42

# How many top hub genes to label in large communities
TOP_N_LABELS = 20


# =============================================================================
# LABEL REPULSION (from Figure 2, do not change)
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
                    label_pos[n1] += push
                    label_pos[n2] -= push
    return label_pos


def draw_labels_with_boxes(ax, pos, label_pos, labels, font_sizes,
                           font_colors, char_w=0.06):
    for n in labels:
        col = font_colors.get(n, "#000000")
        fs = font_sizes.get(n, 12)
        nxy = np.array(pos[n])
        lxy = np.array(label_pos[n])
        offset = np.linalg.norm(lxy - nxy)
        kwargs = dict(
            fontsize=fs, fontweight="bold", color=col,
            ha="center", va="center", zorder=9,
            bbox=dict(boxstyle="round,pad=0.15",
                      fc="white", ec="none", alpha=LABEL_BG_ALPHA))
        if offset > char_w * 2.5:
            ax.annotate(labels[n], xy=nxy, xytext=lxy, **kwargs,
                        arrowprops=dict(arrowstyle="-", color="#888888",
                                        lw=0.6, alpha=0.5))
        else:
            ax.annotate(labels[n], xy=lxy, **kwargs)


def log(msg):
    print(msg, flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# LOAD HARRIS INTERACTORS
# =============================================================================

def load_harris_interactors():
    """Load Harris A3 interactor gene set."""
    if not os.path.exists(HARRIS_PATH):
        log(f"  [WARNING] Harris interactors not found: {HARRIS_PATH}")
        return set()
    harris_df = pd.read_csv(HARRIS_PATH, sep="\t")
    harris_genes = set(harris_df["gene_symbol"].values)
    log(f"  Harris A3 interactors loaded: {len(harris_genes)} genes")
    return harris_genes


# =============================================================================
# LOAD DATA
# =============================================================================

def load_network_data(net_config):
    """Load partition, graph, node scores, and DE stats for one network."""
    net_dir = net_config["dir"]
    result = {"name": net_config["name"], "label": net_config["label"],
              "short": net_config["short"]}

    # Partition
    part_path = os.path.join(net_dir, "04_communities/SC_best_partition.csv")
    if os.path.exists(part_path):
        part_df = pd.read_csv(part_path)
        result["partition"] = part_df
        result["gene_to_comm"] = dict(zip(part_df["gene"], part_df["community"]))
    else:
        log(f"  [SKIP] Partition not found: {part_path}")
        return None

    # Community graph
    graph_path = os.path.join(net_dir, "04_communities/SC_G_comm.gpickle")
    if os.path.exists(graph_path):
        with open(graph_path, "rb") as f:
            result["G_comm"] = pickle.load(f)
    else:
        result["G_comm"] = None

    # Node scores
    scores_path = os.path.join(net_dir,
                               "DIAGNOSTIC_AUDIT/Figure4_Node_Sizing.tsv")
    result["node_intra"] = {}
    if os.path.exists(scores_path):
        scores_df = pd.read_csv(scores_path, sep="\t")
        for _, row in scores_df.iterrows():
            result["node_intra"][row["gene_symbol"]] = float(row["intra_score"])

    # Parameters
    params_path = os.path.join(net_dir,
                               "04_communities/SC_selected_parameters.txt")
    result["params"] = {}
    if os.path.exists(params_path):
        with open(params_path) as f:
            for line in f:
                if "=" in line:
                    k, v = line.strip().split("=", 1)
                    result["params"][k] = v

    # DE stats (for log2FC-based node coloring)
    de_path = os.path.join(net_dir,
                           "02_differential_expression/SC_diffexpr_stats.csv")
    result["de_log2fc"] = {}
    if os.path.exists(de_path):
        de_df = pd.read_csv(de_path)
        result["de_log2fc"] = dict(zip(de_df["gene"], de_df["log2FC"]))
        n_up = sum(1 for v in result["de_log2fc"].values() if v > 0)
        n_down = sum(1 for v in result["de_log2fc"].values() if v < 0)
        log(f"  DE stats loaded: {len(result['de_log2fc'])} genes "
            f"({n_up} up in tumor, {n_down} up in normal)")
    else:
        log(f"  [WARNING] DE stats not found: {de_path}")

    return result


# =============================================================================
# PANEL 1: THREE-POPULATION UMAP
# =============================================================================

def generate_umap_panel(adata, groups_df):
    """Generate the central three-population UMAP."""
    banner("[PANEL 1] Three-Population UMAP")

    # Map groups to adata
    group_map = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))
    adata.obs["three_group"] = adata.obs.index.map(
        lambda x: group_map.get(x, "OTHER"))

    # Subset to cells with group assignments
    mask = adata.obs["three_group"] != "OTHER"
    adata_sub = adata[mask].copy()
    log(f"  Cells in three groups: {mask.sum()}")

    # Also get all basal cells for background
    basal_mask = adata.obs.get("final_annotation",
                               pd.Series("", index=adata.obs.index))
    if hasattr(basal_mask, "str"):
        basal_mask = basal_mask.str.contains("basal", case=False, na=False)
    else:
        basal_mask = pd.Series(False, index=adata.obs.index)

    fig, ax = plt.subplots(figsize=(16, 14))

    # Background: all basal cells in very light gray
    if basal_mask.sum() > 0:
        adata_basal = adata[basal_mask]
        if "X_umap" in adata_basal.obsm:
            umap = adata_basal.obsm["X_umap"]
            ax.scatter(umap[:, 0], umap[:, 1], s=3, c="#f0f0f0",
                       alpha=0.3, rasterized=True, zorder=0)

    if "X_umap" not in adata_sub.obsm:
        log("  [WARNING] X_umap not in adata_sub. Using full adata UMAP.")
        umap_all = adata.obsm.get("X_umap", None)
        if umap_all is not None:
            idx = [list(adata.obs.index).index(c)
                   for c in adata_sub.obs.index
                   if c in set(adata.obs.index)]
            umap_sub = umap_all[idx]
        else:
            log("  [ERROR] No UMAP coordinates available")
            plt.close()
            return
    else:
        umap_sub = adata_sub.obsm["X_umap"]

    groups_array = adata_sub.obs["three_group"].values

    # Plot order: NORMAL first (background), then conditions on top
    for group, color, label, zorder in [
        ("NORMAL",    COLOR_NORMAL,    "NORMAL (n=546)",    1),
        ("CNV_HIGH",  COLOR_CNV_HIGH,  "CNV-HIGH (n=546)",  2),
        ("SBS2_HIGH", COLOR_SBS2_HIGH, "SBS2-HIGH (n=546)", 3),
    ]:
        mask_g = groups_array == group
        if mask_g.sum() > 0:
            ax.scatter(umap_sub[mask_g, 0], umap_sub[mask_g, 1],
                       s=25, c=color, alpha=0.7, edgecolors="#000000",
                       linewidths=0.2, label=label, zorder=zorder,
                       rasterized=True)

    ax.legend(fontsize=FONT_LEGEND, loc="lower left", framealpha=0.9,
              markerscale=3)
    ax.set_xlabel("UMAP 1", fontsize=FONT_AXIS)
    ax.set_ylabel("UMAP 2", fontsize=FONT_AXIS)
    ax.tick_params(labelsize=FONT_TICK)
    ax.set_title("Three-Population Epithelial Cell Framework",
                 fontsize=FONT_TITLE, pad=15)

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Panel_UMAP_three_pop.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] UMAP panel -> {FIG_DIR}")


# =============================================================================
# A3 COMMUNITY ZOOM (DE-colored, Harris-ringed, density-aware)
# =============================================================================


def generate_a3_zoom(net_data, harris_genes):
    """Generate chain-only subgraph plot (no background community).

    Extracts just the concordant chain nodes + A3 genes, lays them out
    independently so they fill the plot space naturally.
    Full community plot generated separately as supplement.
    """
    name = net_data["name"]
    label = net_data["label"]
    short = net_data["short"]
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    node_intra = net_data["node_intra"]
    de_log2fc = net_data.get("de_log2fc", {})

    banner(f"[A3 CHAINS] {name}: {label}")

    if G_comm is None:
        log("  [SKIP] No community graph"); return

    def get_intra(n):
        return node_intra.get(n, 0.1)

    a3_comm = None
    for g in ["APOBEC3A", "APOBEC3B"]:
        if g in gene_to_comm:
            a3_comm = gene_to_comm[g]
            break
    if a3_comm is None:
        log("  [SKIP] A3A/A3B not in partition"); return

    comm_nodes = [g for g, c in gene_to_comm.items() if c == a3_comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
    G_sub.remove_nodes_from(isolates)
    n_nodes = G_sub.number_of_nodes()
    n_edges = G_sub.number_of_edges()
    log(f"  A3 community: C{a3_comm}, {n_nodes} nodes, {n_edges} edges")
    if n_nodes == 0:
        log("  [SKIP] Empty"); return

    a3_here = [n for n in G_sub.nodes() if n in A3_SYMBOLS]
    a3_set = set(a3_here)
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

    # Find chain nodes
    activator_chain_nodes = set()
    repressor_chain_nodes = set()
    for gene in harris_in_comm:
        if de_log2fc.get(gene, 0) > 0 and gene in G_act and G_act.degree(gene) > 0:
            activator_chain_nodes |= set(nx.node_connected_component(G_act, gene))
        if de_log2fc.get(gene, 0) <= 0 and gene in G_rep and G_rep.degree(gene) > 0:
            repressor_chain_nodes |= set(nx.node_connected_component(G_rep, gene))

    # Boundary gene chains
    for a3g in a3_here:
        for nb in G_sub.neighbors(a3g):
            if nb in A3_SYMBOLS: continue
            nb_fc = de_log2fc.get(nb, 0)
            for nn in G_sub.neighbors(nb):
                if nn == a3g or nn in A3_SYMBOLS: continue
                nn_w = G_sub[nb][nn].get("weight", 0)
                nn_fc = de_log2fc.get(nn, 0)
                if nb_fc > 0 and nn_fc > 0 and nn_w > 0:
                    if nb in G_act and G_act.degree(nb) > 0:
                        activator_chain_nodes |= set(nx.node_connected_component(G_act, nb))
                    break
                elif nb_fc <= 0 and nn_fc <= 0 and nn_w < 0:
                    if nb in G_rep and G_rep.degree(nb) > 0:
                        repressor_chain_nodes |= set(nx.node_connected_component(G_rep, nb))
                    break

    activator_chain_nodes -= a3_set
    repressor_chain_nodes -= a3_set
    all_chain_nodes = activator_chain_nodes | repressor_chain_nodes
    featured = all_chain_nodes | a3_set
    harris_in_chains = harris_in_comm & all_chain_nodes

    log(f"  Activating: {len(activator_chain_nodes)}, Inhibiting: {len(repressor_chain_nodes)}")
    log(f"  Harris in chains: {len(harris_in_chains)}")

    # Wall and bridge edges
    wall_edges = []
    concordant_bridge_edges = []
    for a3g in a3_here:
        for nb in G_sub.neighbors(a3g):
            if nb not in all_chain_nodes: continue
            w = G_sub[a3g][nb].get("weight", 0)
            nb_fc = de_log2fc.get(nb, 0)
            if (nb_fc > 0 and w > 0) or (nb_fc <= 0 and w < 0):
                concordant_bridge_edges.append((nb, a3g, w))
            else:
                wall_edges.append((nb, a3g, w))

    log(f"  Wall edges: {len(wall_edges)}, Bridge edges: {len(concordant_bridge_edges)}")

    # ===== EXTRACT CHAIN-ONLY SUBGRAPH =====
    G_chain = G_sub.subgraph(featured).copy()
    chain_iso = [n for n in G_chain.nodes() if G_chain.degree(n) == 0]
    G_chain.remove_nodes_from(chain_iso)
    n_chain = G_chain.number_of_nodes()
    e_chain = G_chain.number_of_edges()
    log(f"  Chain subgraph: {n_chain} nodes, {e_chain} edges")
    if n_chain == 0:
        log("  [SKIP] No chain nodes with edges"); return

    # Layout on JUST the chain subgraph
    if n_chain <= 10:
        pos = nx.spring_layout(G_chain, seed=42, k=3.0, iterations=500, scale=3.0)
    elif n_chain <= 50:
        pos = nx.spring_layout(G_chain, seed=42, weight="abs_weight",
                               k=5.0/np.sqrt(n_chain), iterations=800, scale=4.0)
    else:
        pos = nx.spring_layout(G_chain, seed=42, weight="abs_weight",
                               k=8.0/np.sqrt(n_chain), iterations=1000, scale=6.0)

    def ns_c(n): return 600 + 4000 * get_intra(n)
    def ns_h(n): return 800 + 5000 * get_intra(n)
    def ns_a(n): return max(2000, 1200 + 6000 * get_intra(n))
    def lfs(n): return 16 + 16 * get_intra(n)

    fig, ax = plt.subplots(figsize=(28, 26))

    wall_set = set((min(u,v), max(u,v)) for u,v,_ in wall_edges)
    bridge_set = set((min(u,v), max(u,v)) for u,v,_ in concordant_bridge_edges)

    # zorder 1: chain edges
    for u, v, d in G_chain.edges(data=True):
        ek = (min(u,v), max(u,v))
        if ek in wall_set or ek in bridge_set: continue
        w = d.get("weight", 0)
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        width = 2.5 + 10.0 * abs(w)
        alpha = min(0.75, 0.30 + 0.40 * abs(w))
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, alpha=alpha, linewidth=width, zorder=1)

    # zorder 2: wall edges (dashed) and bridge edges
    for u, v, w in wall_edges:
        if u not in pos or v not in pos: continue
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        width = 3.5 + 10.0 * abs(w)
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, alpha=0.85, linewidth=width, linestyle="--", zorder=2)
    for u, v, w in concordant_bridge_edges:
        if u not in pos or v not in pos: continue
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        width = 3.5 + 10.0 * abs(w)
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, alpha=0.85, linewidth=width, zorder=2)

    # zorder 3: chain nodes
    act_plot = sorted(activator_chain_nodes & set(G_chain.nodes()) - harris_in_chains)
    if act_plot:
        ax.scatter([pos[n][0] for n in act_plot], [pos[n][1] for n in act_plot],
                   s=[ns_c(n) for n in act_plot], c=COLOR_UP_TUMOR,
                   alpha=0.95, edgecolors="#000000", linewidths=1.0, zorder=3)
    rep_plot = sorted(repressor_chain_nodes & set(G_chain.nodes()) - harris_in_chains)
    if rep_plot:
        ax.scatter([pos[n][0] for n in rep_plot], [pos[n][1] for n in rep_plot],
                   s=[ns_c(n) for n in rep_plot], c=COLOR_UP_NORMAL,
                   alpha=0.95, edgecolors="#000000", linewidths=1.0, zorder=3)

    # zorder 4: Harris interactors
    h_list = sorted(harris_in_chains & set(G_chain.nodes()))
    if h_list:
        hc = [COLOR_UP_TUMOR if n in activator_chain_nodes else COLOR_UP_NORMAL for n in h_list]
        hs = [ns_h(n) for n in h_list]
        rs = [s * 2.2 for s in hs]
        ax.scatter([pos[n][0] for n in h_list], [pos[n][1] for n in h_list],
                   s=rs, facecolors="none", edgecolors=COLOR_HARRIS, linewidths=4.5, alpha=1.0, zorder=4)
        ax.scatter([pos[n][0] for n in h_list], [pos[n][1] for n in h_list],
                   s=hs, c=hc, alpha=0.95, edgecolors="#000000", linewidths=1.2, zorder=4)

    # zorder 5: A3 genes
    a3_c = [n for n in a3_here if n in G_chain.nodes()]
    if a3_c:
        asizes = [ns_a(n) for n in a3_c]
        rsizes = [s * 2.0 for s in asizes]
        ax.scatter([pos[n][0] for n in a3_c], [pos[n][1] for n in a3_c],
                   s=rsizes, facecolors="none", edgecolors=COLOR_A3, linewidths=5.0, alpha=1.0, zorder=5)
        ax.scatter([pos[n][0] for n in a3_c], [pos[n][1] for n in a3_c],
                   s=asizes, c=COLOR_A3, alpha=1.0, edgecolors="#000000", linewidths=2.5, zorder=5)

    # zorder 6: Labels
    labels = {}; fsd = {}; fcd = {}
    for n in a3_c:
        labels[n] = A3_ALIAS.get(n, n); fsd[n] = FONT_AXIS; fcd[n] = COLOR_A3_TEXT
    for n in h_list:
        labels[n] = n; fsd[n] = max(lfs(n), 16); fcd[n] = "#5a4a00"
    for n in set(G_chain.nodes()) - a3_set - harris_in_chains:
        labels[n] = n; fsd[n] = max(lfs(n), 14)
        fcd[n] = "#8b4513" if n in activator_chain_nodes else "#2e6e2e"

    if labels:
        lp = repel_labels(pos, labels, data_range=8.0, n_iters=500, push_strength=0.08)
        draw_labels_with_boxes(ax, pos, lp, labels, fsd, fcd, char_w=0.05)

    n_act = len(activator_chain_nodes & set(G_chain.nodes()))
    n_rep = len(repressor_chain_nodes & set(G_chain.nodes()))
    legend_handles = [
        Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3"),
        Patch(fc=COLOR_UP_TUMOR, ec="#000000", lw=0.5, label=f"A3 Activating ({n_act})"),
        Patch(fc=COLOR_UP_NORMAL, ec="#000000", lw=0.5, label=f"A3 Inhibiting ({n_rep})"),
        Patch(fc="white", ec=COLOR_HARRIS, lw=3, label=f"A3 interactor ({len(h_list)})"),
        Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.6, label="Gained co-expression"),
        Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.6, label="Lost co-expression"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=FONT_LEGEND-2, framealpha=0.9)

    thresh = net_data["params"].get("DIFF_THRESHOLD", "?")
    ax.set_title(f"{label} ({short})\n"
                 f"Concordant Chain Subnetwork: {n_chain} genes, {e_chain} edges | "
                 f"threshold={thresh}\n"
                 f"Activating: {n_act}, Inhibiting: {n_rep}, Wall edges: {len(wall_edges)}",
                 fontsize=FONT_TITLE-2, pad=15)
    ax.axis("off")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Panel_{name}_A3_chains.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] {name} chain panel -> {FIG_DIR}")



# =============================================================================
# PANEL 4: SBS2_VS_CNV CHAIN VALIDATION HEATMAP
# =============================================================================

def identify_chain_genes(net_dir, harris_genes):
    """Identify activator and repressor chain genes from one network.
    Same logic as generate_a3_zoom chain detection."""
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

    a3_comm = None
    for g in ["APOBEC3A", "APOBEC3B"]:
        if g in gene_to_comm:
            a3_comm = gene_to_comm[g]
            break
    if a3_comm is None:
        return set(), set()

    comm_nodes = [g for g, c in gene_to_comm.items() if c == a3_comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
    G_sub.remove_nodes_from(isolates)

    a3_set = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)
    harris_in_comm = set(n for n in G_sub.nodes()
                         if n in harris_genes and n not in A3_SYMBOLS)

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

    activator_nodes = set()
    repressor_nodes = set()

    for gene in harris_in_comm:
        if de_log2fc.get(gene, 0) > 0 and gene in G_act and \
           G_act.degree(gene) > 0:
            activator_nodes |= set(nx.node_connected_component(G_act, gene))
        if de_log2fc.get(gene, 0) <= 0 and gene in G_rep and \
           G_rep.degree(gene) > 0:
            repressor_nodes |= set(nx.node_connected_component(G_rep, gene))

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
                        activator_nodes |= set(
                            nx.node_connected_component(G_act, nb))
                    break
                elif nb_fc <= 0 and nn_fc <= 0 and nn_w < 0:
                    if nb in G_rep and G_rep.degree(nb) > 0:
                        repressor_nodes |= set(
                            nx.node_connected_component(G_rep, nb))
                    break

    activator_nodes -= a3_set
    repressor_nodes -= a3_set
    return activator_nodes, repressor_nodes


# =============================================================================
# COMPOSITE FIGURE (UMAP + 2 chain plots + heatmap)
# =============================================================================

def generate_composite():
    """Assemble composite figure with all four panels."""
    banner("[COMPOSITE] Assembling Figure 4")

    umap_path = os.path.join(FIG_DIR, "Panel_UMAP_three_pop.png")
    if not os.path.exists(umap_path):
        log("  [SKIP] UMAP panel not found"); return

    fig = plt.figure(figsize=(36, 28))
    gs = gridspec.GridSpec(2, 2, hspace=0.10, wspace=0.05,
                           height_ratios=[1, 1.2])

    # Top left: UMAP
    ax_umap = fig.add_subplot(gs[0, 0])
    img = plt.imread(umap_path)
    ax_umap.imshow(img)
    ax_umap.axis("off")
    ax_umap.set_title("a", fontsize=36, fontweight="bold",
                      loc="left", pad=10)

    # Top right: SBS2_VS_NORMAL A3 chains
    sbs2_path = os.path.join(FIG_DIR, "Panel_SBS2_VS_NORMAL_A3_chains.png")
    if os.path.exists(sbs2_path):
        ax = fig.add_subplot(gs[0, 1])
        img = plt.imread(sbs2_path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title("b", fontsize=36, fontweight="bold",
                     loc="left", pad=10)
    else:
        log("  [SKIP] SBS2_VS_NORMAL chains panel not found")

    # Bottom left: CNV_VS_NORMAL chains
    cnv_path = os.path.join(FIG_DIR, "Panel_CNV_VS_NORMAL_A3_chains.png")
    if os.path.exists(cnv_path):
        ax = fig.add_subplot(gs[1, 0])
        img = plt.imread(cnv_path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title("c", fontsize=36, fontweight="bold",
                     loc="left", pad=10)
    else:
        log("  [SKIP] CNV_VS_NORMAL chains panel not found")

    # Bottom right: SBS2_VS_CNV dot plot
    heatmap_path = os.path.join(FIG_DIR,
                                "Panel_D_option_C_dotplot.png")
    if os.path.exists(heatmap_path):
        ax = fig.add_subplot(gs[1, 1])
        img = plt.imread(heatmap_path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title("d", fontsize=36, fontweight="bold",
                     loc="left", pad=10)
    else:
        log("  [SKIP] SBS2_VS_CNV dot plot not found")

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Figure4_composite.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Composite figure -> {FIG_DIR}")


# =============================================================================
# SUPPLEMENTAL: FULL NETWORK PLOTS
# =============================================================================

def generate_full_network(net_data, harris_genes):
    """Generate full network plot for one comparison (supplemental)."""
    name = net_data["name"]
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    node_intra = net_data["node_intra"]

    banner(f"[SUPPLEMENT] Full network: {name}")

    if G_comm is None or G_comm.number_of_nodes() == 0:
        log("  [SKIP] No graph"); return

    def get_intra(n):
        return node_intra.get(n, 0.1)

    # Build community structure
    comm_to_nodelist = {}
    for n in G_comm.nodes():
        c = gene_to_comm.get(n, -1)
        comm_to_nodelist.setdefault(c, []).append(n)

    # Classify satellites
    components = list(nx.connected_components(G_comm))
    node_to_comp_size = {}
    for comp in components:
        for n in comp:
            node_to_comp_size[n] = len(comp)

    main_comms = []
    sat_comms = []
    for c, nodes in comm_to_nodelist.items():
        max_comp = max(node_to_comp_size.get(n, 1) for n in nodes)
        if max_comp >= 10:
            main_comms.append(c)
        else:
            sat_comms.append(c)

    unique_comms = sorted(comm_to_nodelist.keys())
    cmap = plt.colormaps["tab20"]
    comm_color_map = {}
    for i, c in enumerate(unique_comms):
        if c in sat_comms:
            comm_color_map[c] = COLOR_SAT
        else:
            comm_color_map[c] = mcolors.to_hex(cmap(i % 20))

    # Layout
    pos = nx.spring_layout(G_comm, seed=COMMUNITY_BASE_SEED,
                           weight="abs_weight",
                           k=3.0 / np.sqrt(max(G_comm.number_of_nodes(), 1)),
                           iterations=200)

    fig, ax = plt.subplots(figsize=(24, 22))

    # Edges (intra only for readability)
    for u, v, d in G_comm.edges(data=True):
        cu = gene_to_comm.get(u, -1)
        cv = gene_to_comm.get(v, -1)
        if cu != cv:
            continue
        w = d.get("weight", 0)
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, alpha=0.08, linewidth=0.3, zorder=0)

    # Nodes by community
    for c in unique_comms:
        c_nodes = comm_to_nodelist.get(c, [])
        if not c_nodes:
            continue
        sizes = [30 + 200 * get_intra(n) for n in c_nodes]
        xs = [pos[n][0] for n in c_nodes]
        ys = [pos[n][1] for n in c_nodes]
        ax.scatter(xs, ys, s=sizes, c=comm_color_map[c],
                   alpha=0.7, edgecolors="#000000", linewidths=0.2,
                   zorder=1)

    # Harris interactor highlights
    harris_in_graph = [g for g in harris_genes
                       if g in G_comm.nodes() and g not in A3_SYMBOLS]
    for n in harris_in_graph:
        ax.scatter(pos[n][0], pos[n][1], s=120,
                   facecolors="none", edgecolors=COLOR_HARRIS,
                   linewidths=2.0, zorder=4)

    # A3 genes
    a3_in = [g for g in A3_SYMBOLS if g in G_comm.nodes()]
    for n in a3_in:
        ax.scatter(pos[n][0], pos[n][1], s=200,
                   c=COLOR_A3, edgecolors="#000000", linewidths=1.5,
                   zorder=5)
        alias = A3_ALIAS.get(n, n)
        ax.annotate(alias, xy=pos[n],
                    xytext=(pos[n][0] + 0.02, pos[n][1] + 0.02),
                    fontsize=14, fontweight="bold", color=COLOR_A3_TEXT,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white",
                              ec="none", alpha=0.7))

    thresh = net_data["params"].get("DIFF_THRESHOLD", "?")
    n_main = len(main_comms)
    n_sat = len(sat_comms)
    n_harris = len(harris_in_graph) + len(
        [g for g in a3_in if g in harris_genes])
    ax.set_title(f"{net_data['label']} ({net_data['short']})\n"
                 f"{G_comm.number_of_nodes()} genes, "
                 f"{G_comm.number_of_edges()} edges | "
                 f"{n_main} main + {n_sat} satellite | "
                 f"threshold={thresh} | "
                 f"{n_harris} Harris interactors",
                 fontsize=FONT_TITLE - 4, pad=15)
    ax.axis("off")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR,
                    f"Supplement_{name}_full_network.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] {name} full network -> {FIG_DIR}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("GENERATE FIGURE 4 PANELS")

    # Load Harris interactors
    log("Loading Harris A3 interactors...")
    harris_genes = load_harris_interactors()

    # Load adata for UMAP
    log("\nLoading adata...")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Load group assignments
    groups_df = pd.read_csv(GROUP_PATH, sep="\t")
    log(f"  Groups: {len(groups_df)} cells")

    # Panel 1: UMAP
    generate_umap_panel(adata, groups_df)

    # Free adata memory
    del adata

    # Load network data (SBS2_VS_NORMAL and CNV_VS_NORMAL only)
    all_net_data = []
    for net in NETWORKS:
        log(f"\nLoading {net['name']}...")
        data = load_network_data(net)
        if data is not None:
            all_net_data.append(data)
            log(f"  Graph: {data['G_comm'].number_of_nodes()} nodes, "
                f"{data['G_comm'].number_of_edges()} edges")

    # Panel 2+3: A3 community zooms (normal comparisons only)
    for net_data in all_net_data:
        generate_a3_zoom(net_data, harris_genes)

    # Panel 4: SBS2_VS_CNV chain validation heatmap
    # Panel 4: run Generate_Panel_D_Options.py separately

    # Composite figure
    generate_composite()

    # Supplemental: full networks
    for net_data in all_net_data:
        generate_full_network(net_data, harris_genes)

    banner("FIGURE 4 COMPLETE")
    log(f"\nOutput: {FIG_DIR}")
    log(f"Files:")
    for f in sorted(os.listdir(FIG_DIR)):
        if not os.path.isdir(os.path.join(FIG_DIR, f)):
            log(f"  {f}")


if __name__ == "__main__":
    main()
