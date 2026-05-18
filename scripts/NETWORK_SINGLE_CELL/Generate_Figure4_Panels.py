#!/usr/bin/env python3
"""
Generate_Figure4_Panels.py
===========================

Generate publication-quality Figure 4 panels for the three-network
single-cell analysis.

Layout concept:
  - Central panel: Three-population UMAP (SBS2-HIGH / CNV-HIGH / NORMAL)
  - Three surrounding A3 community zoom panels (one per comparison)
  - Supplemental: Full network plots + non-A3 community zooms

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
from matplotlib.patches import Patch, FancyArrowPatch
import matplotlib.gridspec as gridspec

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
FIG_DIR = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS")
os.makedirs(FIG_DIR, exist_ok=True)

ADATA_PATH = os.path.join(FIG4_ROOT, "00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(FIG4_ROOT, "01_group_selection/three_group_assignments.tsv")

NETWORKS = [
    {
        "name": "SBS2_VS_CNV",
        "label": "SBS2-HIGH vs CNV-HIGH",
        "short": "Divergent fates",
        "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_CNV"),
    },
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

# A3 gene identifiers (symbol format for SC)
A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}

# =============================================================================
# STYLE (inherited from Figure 2, do not change)
# =============================================================================

COLOR_SBS2_HIGH = "#ed6a5a"   # coral
COLOR_CNV_HIGH  = "#4682b4"   # steelblue
COLOR_NORMAL    = "#d0d0d0"   # gray
COLOR_EDGE_POS  = "#b22222"   # firebrick
COLOR_EDGE_NEG  = "#4682b4"   # steelblue
COLOR_A3        = "#ed6a5a"   # coral
COLOR_A3_TEXT   = "#c0392b"
COLOR_SAT       = "#e8913a"   # orange for satellites

FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22
LABEL_BG_ALPHA = 0.42

COMMUNITY_BASE_SEED = 42


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
            ha="center", va="center",
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
# LOAD DATA
# =============================================================================

def load_network_data(net_config):
    """Load partition, graph, and node scores for one network."""
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

    # Correlation matrix (for A3 neighborhood coloring)
    corr_path = os.path.join(net_dir,
                             "03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl")
    if os.path.exists(corr_path):
        with open(corr_path, "rb") as f:
            result["corr_diff"] = pickle.load(f)
    else:
        result["corr_diff"] = None

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

    # Subset to basal cells with group assignments
    mask = adata.obs["three_group"] != "OTHER"
    adata_sub = adata[mask].copy()
    log(f"  Cells in three groups: {mask.sum()}")

    # Also get all basal cells for background
    basal_mask = adata.obs.get("final_annotation", pd.Series("", index=adata.obs.index))
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
        log("  [WARNING] X_umap not in adata. Using full adata UMAP coords.")
        # Try to get UMAP from full adata for the subset cells
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

    # Plot order: NORMAL first (background), then the two conditions
    for group, color, label, zorder in [
        ("NORMAL", COLOR_NORMAL, "NORMAL (n=546)", 1),
        ("CNV_HIGH", COLOR_CNV_HIGH, "CNV-HIGH (n=546)", 2),
        ("SBS2_HIGH", COLOR_SBS2_HIGH, "SBS2-HIGH (n=546)", 3),
    ]:
        mask_g = groups_array == group
        if mask_g.sum() > 0:
            ax.scatter(umap_sub[mask_g, 0], umap_sub[mask_g, 1],
                       s=25, c=color, alpha=0.7, edgecolors="#000000",
                       linewidths=0.2, label=label, zorder=zorder,
                       rasterized=True)

    ax.legend(fontsize=FONT_LEGEND, loc="upper right", framealpha=0.9,
              markerscale=3)
    ax.set_xlabel("UMAP 1", fontsize=FONT_AXIS)
    ax.set_ylabel("UMAP 2", fontsize=FONT_AXIS)
    ax.tick_params(labelsize=FONT_TICK)
    ax.set_title("Three-Population Basal Cell Framework",
                 fontsize=FONT_TITLE, pad=15)

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Panel_UMAP_three_pop.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] UMAP panel -> {FIG_DIR}")


# =============================================================================
# A3 COMMUNITY ZOOM
# =============================================================================

def generate_a3_zoom(net_data):
    """Generate A3 community zoom for one network."""
    name = net_data["name"]
    label = net_data["label"]
    short = net_data["short"]
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    node_intra = net_data["node_intra"]
    corr_diff = net_data.get("corr_diff", None)

    banner(f"[A3 ZOOM] {name}: {label}")

    if G_comm is None:
        log("  [SKIP] No community graph"); return

    def get_intra(n):
        return node_intra.get(n, 0.1)

    # Find A3 community
    a3_in_graph = [g for g in A3_SYMBOLS if g in G_comm.nodes()]
    if not a3_in_graph:
        log("  [SKIP] No A3 genes in graph"); return

    # Get the community containing A3A or A3B
    a3_comm = None
    for g in ["APOBEC3A", "APOBEC3B"]:
        if g in gene_to_comm:
            a3_comm = gene_to_comm[g]
            break
    if a3_comm is None:
        log("  [SKIP] A3A/A3B not in partition"); return

    # Get community nodes
    comm_nodes = [g for g, c in gene_to_comm.items() if c == a3_comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()

    # Remove isolates
    isolates = [n for n in G_sub.nodes() if G_sub.degree(n) == 0]
    G_sub.remove_nodes_from(isolates)

    log(f"  A3 community: C{a3_comm}, {G_sub.number_of_nodes()} nodes, "
        f"{G_sub.number_of_edges()} edges")

    if G_sub.number_of_nodes() == 0:
        log("  [SKIP] Empty after isolate removal"); return

    # Layout
    if G_sub.number_of_nodes() <= 5:
        pos = nx.spring_layout(G_sub, seed=COMMUNITY_BASE_SEED,
                               k=2.0, iterations=200, scale=1.0)
    else:
        pos = nx.spring_layout(
            G_sub, seed=COMMUNITY_BASE_SEED, weight="abs_weight",
            k=4.0 / np.sqrt(max(G_sub.number_of_nodes(), 1)),
            iterations=500, scale=2.0)

    # Color nodes by A3 correlation if available
    node_colors = {}
    a3_seed = "APOBEC3A" if "APOBEC3A" in G_sub else "APOBEC3B"

    if corr_diff is not None and a3_seed in corr_diff.index:
        for n in G_sub.nodes():
            if n == a3_seed:
                node_colors[n] = COLOR_A3
            elif n in A3_SYMBOLS:
                node_colors[n] = COLOR_A3
            elif n in corr_diff.index:
                rho = corr_diff.loc[a3_seed, n]
                if rho > 0:
                    node_colors[n] = "#f4a582"  # light coral (positive)
                else:
                    node_colors[n] = "#92c5de"  # light blue (negative)
            else:
                node_colors[n] = "#d0d0d0"
    else:
        for n in G_sub.nodes():
            node_colors[n] = COLOR_A3 if n in A3_SYMBOLS else "#92c5de"

    # Sizing
    def ns(n):
        if G_sub.number_of_nodes() <= 10:
            return 2000 + 6000 * get_intra(n)
        else:
            return 750 + 5000 * get_intra(n)

    def fs(n):
        if G_sub.number_of_nodes() <= 10:
            return 22 + 20 * get_intra(n)
        else:
            return 14 + 24 * get_intra(n)

    fig, ax = plt.subplots(figsize=(20, 18))

    # Edges
    for u, v, d in G_sub.edges(data=True):
        w = d.get("weight", 0)
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        width = 2.0 + 10.0 * abs(w)
        alpha = min(0.6, 0.20 + 0.35 * abs(w))
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, alpha=alpha, linewidth=width, zorder=1)

    # Nodes
    a3_here = [n for n in G_sub.nodes() if n in A3_SYMBOLS]
    regular = [n for n in G_sub.nodes() if n not in A3_SYMBOLS]

    if regular:
        r_colors = [node_colors[n] for n in regular]
        r_sizes = [ns(n) for n in regular]
        nx.draw_networkx_nodes(G_sub, pos, nodelist=regular, ax=ax,
                               node_size=r_sizes, node_color=r_colors,
                               alpha=0.90, edgecolors="#000000",
                               linewidths=0.8, zorder=2)

    if a3_here:
        a3_sizes = [ns(n) * 1.5 for n in a3_here]
        nx.draw_networkx_nodes(G_sub, pos, nodelist=a3_here, ax=ax,
                               node_size=a3_sizes, node_color=COLOR_A3,
                               alpha=1.0, edgecolors="#000000",
                               linewidths=2.5, zorder=3)
        ring_sizes = [s * 2.0 for s in a3_sizes]
        nx.draw_networkx_nodes(G_sub, pos, nodelist=a3_here, ax=ax,
                               node_size=ring_sizes, node_color="none",
                               edgecolors=COLOR_A3, linewidths=5.0,
                               alpha=1.0, zorder=3)

    # Labels
    labels = {}; font_sizes = {}; font_colors = {}

    # Always label A3 genes
    for n in a3_here:
        alias = A3_ALIAS.get(n, n)
        labels[n] = alias
        font_sizes[n] = FONT_AXIS
        font_colors[n] = COLOR_A3_TEXT

    # Label all nodes for small communities, top hubs for large
    if G_sub.number_of_nodes() <= 20:
        for n in G_sub.nodes():
            if n not in labels:
                labels[n] = n
                font_sizes[n] = fs(n)
                font_colors[n] = "#000000"
    else:
        scored = [(n, get_intra(n)) for n in G_sub.nodes()]
        scored.sort(key=lambda x: -x[1])
        for n, sc in scored[:15]:
            if n not in labels:
                labels[n] = n
                font_sizes[n] = fs(n)
                font_colors[n] = "#000000"

    if labels:
        dr = 5.0 if G_sub.number_of_nodes() <= 10 else 5.0
        lp = repel_labels(pos, labels, data_range=dr,
                          n_iters=250, push_strength=0.07)
        draw_labels_with_boxes(ax, pos, lp, labels,
                               font_sizes, font_colors, char_w=0.05)

    # Legend
    legend_handles = [
        Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3"),
        Patch(fc="#f4a582", ec="#000000", lw=0.5,
              label=f"Positive DIFF with {A3_ALIAS.get(a3_seed, a3_seed)}"),
        Patch(fc="#92c5de", ec="#000000", lw=0.5,
              label=f"Negative DIFF with {A3_ALIAS.get(a3_seed, a3_seed)}"),
        Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.5,
              label="Gained co-expression"),
        Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.5,
              label="Lost co-expression"),
    ]
    ax.legend(handles=legend_handles, loc="upper left",
              fontsize=FONT_LEGEND - 2, framealpha=0.9)

    thresh = net_data["params"].get("DIFF_THRESHOLD", "?")
    ax.set_title(f"{label} ({short})\n"
                 f"A3 Community C{a3_comm}: "
                 f"{G_sub.number_of_nodes()} genes, "
                 f"{G_sub.number_of_edges()} edges | "
                 f"threshold={thresh}",
                 fontsize=FONT_TITLE - 2, pad=15)
    ax.axis("off")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR,
                    f"Panel_{name}_A3_zoom.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] {name} A3 zoom -> {FIG_DIR}")


# =============================================================================
# COMPOSITE FIGURE (UMAP + 3 zooms)
# =============================================================================

def generate_composite(zoom_paths):
    """Assemble the composite figure with UMAP center and 3 zoom panels."""
    banner("[COMPOSITE] Assembling Figure 4")

    # Check if all component images exist
    umap_path = os.path.join(FIG_DIR, "Panel_UMAP_three_pop.png")
    if not os.path.exists(umap_path):
        log("  [SKIP] UMAP panel not found"); return

    fig = plt.figure(figsize=(36, 28))
    gs = gridspec.GridSpec(2, 3, hspace=0.15, wspace=0.1,
                           height_ratios=[1.2, 1])

    # Top center: UMAP (spans 2 columns)
    ax_umap = fig.add_subplot(gs[0, 0:2])
    img = plt.imread(umap_path)
    ax_umap.imshow(img)
    ax_umap.axis("off")
    ax_umap.set_title("a", fontsize=36, fontweight="bold",
                      loc="left", pad=10)

    # Top right: SBS2_VS_CNV zoom
    panel_labels = ["b", "c", "d"]
    positions = [(0, 2), (1, 0), (1, 1)]  # grid positions

    for i, net_name in enumerate(["SBS2_VS_CNV", "SBS2_VS_NORMAL",
                                   "CNV_VS_NORMAL"]):
        zoom_path = os.path.join(FIG_DIR,
                                 f"Panel_{net_name}_A3_zoom.png")
        if not os.path.exists(zoom_path):
            log(f"  [SKIP] {net_name} zoom not found")
            continue

        row, col = positions[i]
        ax = fig.add_subplot(gs[row, col])
        img = plt.imread(zoom_path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(panel_labels[i], fontsize=36, fontweight="bold",
                     loc="left", pad=10)

    # Arrow annotations between UMAP and zoom panels would go here
    # (requires fine-tuning of coordinates after seeing initial output)

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Figure4_composite.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Composite figure -> {FIG_DIR}")


# =============================================================================
# SUPPLEMENTAL: FULL NETWORK PLOTS
# =============================================================================

def generate_full_network(net_data):
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

    # Layout (spring, generous spacing for main communities)
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

    # Nodes
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
    ax.set_title(f"{net_data['label']} ({net_data['short']})\n"
                 f"{G_comm.number_of_nodes()} genes, "
                 f"{G_comm.number_of_edges()} edges | "
                 f"{n_main} main + {n_sat} satellite | "
                 f"threshold={thresh}",
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

    # Load adata for UMAP
    log("Loading adata...")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Load group assignments
    groups_df = pd.read_csv(GROUP_PATH, sep="\t")
    log(f"  Groups: {len(groups_df)} cells")

    # Panel 1: UMAP
    generate_umap_panel(adata, groups_df)

    # Free adata memory
    del adata

    # Load all network data
    all_net_data = []
    for net in NETWORKS:
        log(f"\nLoading {net['name']}...")
        data = load_network_data(net)
        if data is not None:
            all_net_data.append(data)
            log(f"  Loaded: {data['G_comm'].number_of_nodes()} nodes")

    # Panel 2: A3 community zooms
    for net_data in all_net_data:
        generate_a3_zoom(net_data)
        # Free correlation matrix after use
        if "corr_diff" in net_data:
            del net_data["corr_diff"]

    # Composite figure
    generate_composite([])

    # Supplemental: full networks
    for net_data in all_net_data:
        generate_full_network(net_data)

    banner("FIGURE 4 COMPLETE")
    log(f"\nOutput: {FIG_DIR}")
    log(f"Files:")
    for f in sorted(os.listdir(FIG_DIR)):
        if not os.path.isdir(os.path.join(FIG_DIR, f)):
            log(f"  {f}")


if __name__ == "__main__":
    main()
