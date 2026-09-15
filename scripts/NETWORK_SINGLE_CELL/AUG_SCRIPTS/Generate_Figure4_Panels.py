#!/usr/bin/env python3
"""
Generate_Figure4_Panels_v2.py
=============================

Figure 4 main panels for the three-network single-cell analysis.

  a        Three-population UMAP (SBS2-HIGH / CNV-HIGH / NORMAL)
  b        SBS2-HIGH vs NORMAL, A3 community C0, concordant chains (both)
  c-top    CNV-HIGH vs NORMAL, A3A community C1, ACTIVATING chain only
  c-bottom CNV-HIGH vs NORMAL, A3A community C1, INHIBITING chain only
  d        CNV-HIGH vs NORMAL, A3B community C0, concordant chains (both)
  composite  a + b (left), c-top + c-bottom (right), d (below)

Chain logic is the prior version's, with three additions:
  focus_gene   which enzyme's community to render.
  mode         "both" / "activator" / "repressor" filter.
  HOP_CUTOFF   concordant chains with a node within this many hops of A3 are
               plotted (default 2), and the one-hop bridge node connecting A3
               to a 2-hop chain is included so A3 itself is drawn connected.
Edge styling: only the discordant negative edges sitting directly on A3 (the
A3 wall) are dashed; connecting bridge edges and chain-internal edges are solid.

Style: 28-34 pt text, hex colors, PDF + PNG at 300 DPI.

Usage:
    conda run -n NETWORK python Generate_Figure4_Panels_v2.py
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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
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
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")

SBS2_DIR = os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL")
CNV_DIR = os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL")

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}

# Chain reach: include concordant chains with a node within this many hops of A3
HOP_CUTOFF = 2

# Palette
COLOR_SBS2_HIGH = "#ed6a5a"
COLOR_CNV_HIGH  = "#F6D155"
COLOR_NORMAL    = "#4682b4"
COLOR_UP_TUMOR  = "#fcad61"
COLOR_UP_NORMAL = "#aadda4"
COLOR_A3        = "#ed6a5a"
COLOR_A3_TEXT   = "#c0392b"
COLOR_HARRIS    = "#F6D155"
COLOR_BRIDGE    = "#cccccc"
COLOR_EDGE_POS  = "#b22222"
COLOR_EDGE_NEG  = "#4682b4"

FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22
LABEL_BG_ALPHA = 0.42


def log(msg):
    print(msg, flush=True)

def banner(t, c="="):
    print(f"\n{c*80}\n  {t}\n{c*80}", flush=True)


# =============================================================================
# LABEL ENGINE (unchanged from Figure 2 / prior version)
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
        nxy = np.array(pos[n]); lxy = np.array(label_pos[n])
        offset = np.linalg.norm(lxy - nxy)
        kwargs = dict(fontsize=fs, fontweight="bold", color=col,
                      ha="center", va="center", zorder=9,
                      bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                ec="none", alpha=LABEL_BG_ALPHA))
        if offset > char_w * 2.5:
            ax.annotate(labels[n], xy=nxy, xytext=lxy, **kwargs,
                        arrowprops=dict(arrowstyle="-", color="#888888",
                                        lw=0.6, alpha=0.5))
        else:
            ax.annotate(labels[n], xy=lxy, **kwargs)


# =============================================================================
# LOADERS
# =============================================================================

def load_harris_interactors():
    if not os.path.exists(HARRIS_PATH):
        log(f"  [WARNING] Harris not found: {HARRIS_PATH}"); return set()
    df = pd.read_csv(HARRIS_PATH, sep="\t")
    genes = set(df["gene_symbol"].values)
    log(f"  Harris A3 interactors: {len(genes)}")
    return genes


def load_network_data(net_dir, name, label, short):
    result = {"name": name, "label": label, "short": short}
    part = pd.read_csv(os.path.join(net_dir, "04_communities/SC_best_partition.csv"))
    result["gene_to_comm"] = dict(zip(part["gene"], part["community"]))
    with open(os.path.join(net_dir, "04_communities/SC_G_comm.gpickle"), "rb") as f:
        result["G_comm"] = pickle.load(f)
    sizing = os.path.join(net_dir, "DIAGNOSTIC_AUDIT/Figure4_Node_Sizing.tsv")
    result["node_intra"] = {}
    if os.path.exists(sizing):
        s = pd.read_csv(sizing, sep="\t")
        result["node_intra"] = {r["gene_symbol"]: float(r["intra_score"])
                                for _, r in s.iterrows()}
    params = os.path.join(net_dir, "04_communities/SC_selected_parameters.txt")
    result["params"] = {}
    if os.path.exists(params):
        with open(params) as f:
            for line in f:
                if "=" in line:
                    k, v = line.strip().split("=", 1); result["params"][k] = v
    de = os.path.join(net_dir, "02_differential_expression/SC_diffexpr_stats.csv")
    de_df = pd.read_csv(de)
    result["de_log2fc"] = dict(zip(de_df["gene"], de_df["log2FC"]))
    return result


# =============================================================================
# PANEL a: UMAP
# =============================================================================

def generate_umap_panel(adata, groups_df):
    banner("[PANEL a] Three-population UMAP")
    gmap = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))
    adata.obs["three_group"] = adata.obs.index.map(lambda x: gmap.get(x, "OTHER"))
    sub = adata[adata.obs["three_group"] != "OTHER"].copy()
    umap = sub.obsm["X_umap"]
    grp = sub.obs["three_group"].values

    fig, ax = plt.subplots(figsize=(16, 14))
    if "X_umap" in adata.obsm:
        bg = adata.obsm["X_umap"]
        ax.scatter(bg[:, 0], bg[:, 1], s=3, c="#f0f0f0", alpha=0.3,
                   rasterized=True, zorder=0)
    for g, color, lab, z in [("NORMAL", COLOR_NORMAL, "NORMAL (n=546)", 1),
                             ("CNV_HIGH", COLOR_CNV_HIGH, "CNV-HIGH (n=546)", 2),
                             ("SBS2_HIGH", COLOR_SBS2_HIGH, "SBS2-HIGH (n=546)", 3)]:
        m = grp == g
        if m.sum():
            ax.scatter(umap[m, 0], umap[m, 1], s=25, c=color, alpha=0.7,
                       edgecolors="#000000", linewidths=0.2, label=lab,
                       zorder=z, rasterized=True)
    ax.legend(fontsize=FONT_LEGEND, loc="lower left", framealpha=0.9, markerscale=3)
    ax.set_xlabel("UMAP 1", fontsize=FONT_AXIS); ax.set_ylabel("UMAP 2", fontsize=FONT_AXIS)
    ax.tick_params(labelsize=FONT_TICK)
    ax.set_title("Three-Population Epithelial Cell Framework", fontsize=FONT_TITLE, pad=15)
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Panel_a_UMAP.{ext}"), dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Panel_a_UMAP")


# =============================================================================
# CHAIN ZOOM (prior logic + 2-hop reach, bridge nodes, A3-only dashed walls)
# =============================================================================

def generate_a3_zoom(net_data, harris_genes, focus_gene, mode, out_tag, title_suffix=""):
    name, label, short = net_data["name"], net_data["label"], net_data["short"]
    G_comm = net_data["G_comm"]
    gene_to_comm = net_data["gene_to_comm"]
    node_intra = net_data["node_intra"]
    de_log2fc = net_data["de_log2fc"]

    banner(f"[CHAIN] {out_tag}: {label} | focus {focus_gene} | mode {mode} | hops {HOP_CUTOFF}")

    def get_intra(n): return node_intra.get(n, 0.1)
    def dirof(g): return "up_tumor" if de_log2fc.get(g, 0) > 0 else "up_normal"

    if focus_gene not in gene_to_comm:
        log(f"  [SKIP] {focus_gene} not in partition"); return
    a3_comm = gene_to_comm[focus_gene]

    comm_nodes = [g for g, c in gene_to_comm.items() if c == a3_comm]
    G_sub = G_comm.subgraph(comm_nodes).copy()
    G_sub.remove_nodes_from([n for n in G_sub.nodes() if G_sub.degree(n) == 0])
    log(f"  Community C{a3_comm}: {G_sub.number_of_nodes()} nodes, {G_sub.number_of_edges()} edges")
    if G_sub.number_of_nodes() == 0:
        log("  [SKIP] empty"); return

    a3_here = [n for n in G_sub.nodes() if n in A3_SYMBOLS]
    a3_set = set(a3_here)
    harris_in_comm = set(n for n in G_sub.nodes()
                         if n in harris_genes and n not in A3_SYMBOLS)

    def build_conc(m):
        a3_nodes = set(n for n in G_sub.nodes() if n in A3_SYMBOLS)
        Gc = nx.Graph()
        if m == "activator":
            elig = set(n for n in G_sub.nodes() if de_log2fc.get(n, 0) > 0) | a3_nodes
            for u, v, d in G_sub.edges(data=True):
                w = d.get("weight", 0)
                if w > 0 and u in elig and v in elig: Gc.add_edge(u, v, weight=w)
        else:
            elig = set(n for n in G_sub.nodes() if de_log2fc.get(n, 0) <= 0) | a3_nodes
            for u, v, d in G_sub.edges(data=True):
                w = d.get("weight", 0)
                if w < 0 and u in elig and v in elig: Gc.add_edge(u, v, weight=w)
        return Gc

    G_act = build_conc("activator")
    G_rep = build_conc("repressor")

    activator_chain_nodes, repressor_chain_nodes = set(), set()
    # (1) Harris-anchored chains (any distance)
    for gene in harris_in_comm:
        if de_log2fc.get(gene, 0) > 0 and gene in G_act and G_act.degree(gene) > 0:
            activator_chain_nodes |= set(nx.node_connected_component(G_act, gene))
        if de_log2fc.get(gene, 0) <= 0 and gene in G_rep and G_rep.degree(gene) > 0:
            repressor_chain_nodes |= set(nx.node_connected_component(G_rep, gene))
    # (2) every concordant chain with a node within HOP_CUTOFF hops of A3
    near_a3 = set()
    for a3g in a3_here:
        near_a3 |= set(nx.single_source_shortest_path_length(
            G_sub, a3g, cutoff=HOP_CUTOFF).keys())
    for node in near_a3 - a3_set:
        if node in G_act and G_act.degree(node) > 0:
            activator_chain_nodes |= set(nx.node_connected_component(G_act, node))
        if node in G_rep and G_rep.degree(node) > 0:
            repressor_chain_nodes |= set(nx.node_connected_component(G_rep, node))

    activator_chain_nodes -= a3_set
    repressor_chain_nodes -= a3_set
    # mode filter
    if mode == "activator":
        repressor_chain_nodes = set()
    elif mode == "repressor":
        activator_chain_nodes = set()

    all_chain_nodes = activator_chain_nodes | repressor_chain_nodes
    # bridge intermediates: one-hop A3 neighbors adjacent to a chain node, so A3
    # plots connected when the chain is two hops out (e.g. the RALY activator).
    bridge_nodes = set()
    for a3g in a3_here:
        for x in G_sub.neighbors(a3g):
            if x in A3_SYMBOLS or x in all_chain_nodes:
                continue
            if any(y in all_chain_nodes for y in G_sub.neighbors(x)):
                bridge_nodes.add(x)
    featured = all_chain_nodes | a3_set | bridge_nodes
    harris_in_chains = harris_in_comm & all_chain_nodes
    log(f"  Activating: {len(activator_chain_nodes)}, Inhibiting: {len(repressor_chain_nodes)}, "
        f"bridges: {len(bridge_nodes)}, Harris in chains: {len(harris_in_chains)}")

    G_chain = G_sub.subgraph(featured).copy()
    G_chain.remove_nodes_from([n for n in G_chain.nodes() if G_chain.degree(n) == 0])
    n_chain, e_chain = G_chain.number_of_nodes(), G_chain.number_of_edges()
    log(f"  Chain subgraph: {n_chain} nodes, {e_chain} edges")
    if n_chain == 0:
        log("  [NOTE] no connected chain nodes; drawing A3 + any fragments")

    layout_nodes = list(G_chain.nodes()) if n_chain > 0 else list(featured & set(G_sub.nodes()))
    G_layout = G_sub.subgraph(layout_nodes).copy()
    nL = G_layout.number_of_nodes()
    if nL == 0:
        log("  [SKIP] nothing to draw"); return
    if nL <= 10:
        pos = nx.spring_layout(G_layout, seed=42, k=3.0, iterations=500, scale=3.0)
    elif nL <= 50:
        pos = nx.spring_layout(G_layout, seed=42, weight="abs_weight",
                               k=5.0 / np.sqrt(nL), iterations=800, scale=4.0)
    else:
        pos = nx.spring_layout(G_layout, seed=42, weight="abs_weight",
                               k=8.0 / np.sqrt(nL), iterations=1000, scale=6.0)

    def ns_c(n): return 600 + 4000 * get_intra(n)
    def ns_h(n): return 800 + 5000 * get_intra(n)
    def ns_a(n): return max(2000, 1200 + 6000 * get_intra(n))
    def lfs(n): return 16 + 16 * get_intra(n)

    fig, ax = plt.subplots(figsize=(28, 26))

    # edges: only discordant negative edges touching A3 are dashed (the A3 wall);
    # concordant A3 edges, bridge-to-chain edges, and chain-internal edges solid.
    n_wall = 0
    for u, v, d in G_layout.edges(data=True):
        w = d.get("weight", 0.0)
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        if (u in A3_SYMBOLS) or (v in A3_SYMBOLS):
            nb = v if u in A3_SYMBOLS else u
            concordant = ((dirof(nb) == "up_tumor" and w > 0) or
                          (dirof(nb) == "up_normal" and w < 0))
            if concordant:
                ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                        color=color, alpha=0.85, linewidth=3.5 + 10.0 * abs(w), zorder=2)
            else:
                n_wall += 1
                ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                        color=color, alpha=0.85, linewidth=3.5 + 10.0 * abs(w),
                        linestyle="--", zorder=2)
        else:
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color=color, alpha=min(0.75, 0.30 + 0.40 * abs(w)),
                    linewidth=2.5 + 10.0 * abs(w), zorder=1)

    drawn = set(G_layout.nodes())
    act_plot = sorted((activator_chain_nodes & drawn) - harris_in_chains)
    rep_plot = sorted((repressor_chain_nodes & drawn) - harris_in_chains)
    bridge_plot = sorted((bridge_nodes & drawn) - harris_in_chains - a3_set
                         - activator_chain_nodes - repressor_chain_nodes)
    if act_plot:
        ax.scatter([pos[n][0] for n in act_plot], [pos[n][1] for n in act_plot],
                   s=[ns_c(n) for n in act_plot], c=COLOR_UP_TUMOR, alpha=0.95,
                   edgecolors="#000000", linewidths=1.0, zorder=3)
    if rep_plot:
        ax.scatter([pos[n][0] for n in rep_plot], [pos[n][1] for n in rep_plot],
                   s=[ns_c(n) for n in rep_plot], c=COLOR_UP_NORMAL, alpha=0.95,
                   edgecolors="#000000", linewidths=1.0, zorder=3)
    if bridge_plot:
        ax.scatter([pos[n][0] for n in bridge_plot], [pos[n][1] for n in bridge_plot],
                   s=[0.55 * ns_c(n) for n in bridge_plot], c=COLOR_BRIDGE, alpha=0.9,
                   edgecolors="#000000", linewidths=0.8, zorder=3)
    h_list = sorted(harris_in_chains & drawn)
    if h_list:
        hc = [COLOR_UP_TUMOR if n in activator_chain_nodes else COLOR_UP_NORMAL for n in h_list]
        hs = [ns_h(n) for n in h_list]
        ax.scatter([pos[n][0] for n in h_list], [pos[n][1] for n in h_list],
                   s=[s * 2.2 for s in hs], facecolors="none", edgecolors=COLOR_HARRIS,
                   linewidths=4.5, zorder=4)
        ax.scatter([pos[n][0] for n in h_list], [pos[n][1] for n in h_list],
                   s=hs, c=hc, alpha=0.95, edgecolors="#000000", linewidths=1.2, zorder=4)
    a3_c = [n for n in a3_here if n in drawn]
    if a3_c:
        asz = [ns_a(n) for n in a3_c]
        ax.scatter([pos[n][0] for n in a3_c], [pos[n][1] for n in a3_c],
                   s=[s * 2.0 for s in asz], facecolors="none", edgecolors=COLOR_A3,
                   linewidths=5.0, zorder=5)
        ax.scatter([pos[n][0] for n in a3_c], [pos[n][1] for n in a3_c],
                   s=asz, c=COLOR_A3, alpha=1.0, edgecolors="#000000", linewidths=2.5, zorder=5)

    labels, fsd, fcd = {}, {}, {}
    for n in a3_c:
        labels[n] = A3_ALIAS.get(n, n); fsd[n] = FONT_AXIS; fcd[n] = COLOR_A3_TEXT
    for n in h_list:
        labels[n] = n; fsd[n] = max(lfs(n), 16); fcd[n] = "#5a4a00"
    for n in drawn - a3_set - harris_in_chains:
        labels[n] = n; fsd[n] = max(lfs(n), 14)
        if n in activator_chain_nodes:
            fcd[n] = "#8b4513"
        elif n in repressor_chain_nodes:
            fcd[n] = "#2e6e2e"
        else:
            fcd[n] = "#555555"   # bridge intermediate
    if labels:
        lp = repel_labels(pos, labels, data_range=8.0, n_iters=500, push_strength=0.08)
        draw_labels_with_boxes(ax, pos, lp, labels, fsd, fcd, char_w=0.05)

    n_act = len(activator_chain_nodes & drawn)
    n_rep = len(repressor_chain_nodes & drawn)
    legend_handles = [Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3")]
    if mode in ("both", "activator"):
        legend_handles.append(Patch(fc=COLOR_UP_TUMOR, ec="#000000", lw=0.5,
                                    label=f"A3 Activating ({n_act})"))
    if mode in ("both", "repressor"):
        legend_handles.append(Patch(fc=COLOR_UP_NORMAL, ec="#000000", lw=0.5,
                                    label=f"A3 Inhibiting ({n_rep})"))
    legend_handles.append(Patch(fc="white", ec=COLOR_HARRIS, lw=3,
                                label=f"A3 interactor ({len(h_list)})"))
    if bridge_plot:
        legend_handles.append(Patch(fc=COLOR_BRIDGE, ec="#000000",
                                    label=f"Bridge to chain ({len(bridge_plot)})"))
    legend_handles += [
        Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.6, label="Gained co-expression"),
        Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.6, label="Lost co-expression"),
        Line2D([0], [0], color=COLOR_EDGE_NEG, lw=3, linestyle="--", label="A3 wall edge"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=FONT_LEGEND - 2, framealpha=0.9)

    thr = net_data["params"].get("DIFF_THRESHOLD", "?")
    suff = f" — {title_suffix}" if title_suffix else ""
    ax.set_title(f"{label} ({short}){suff}\n"
                 f"C{a3_comm}: {n_chain} genes, {e_chain} edges | threshold={thr}\n"
                 f"Activating: {n_act}, Inhibiting: {n_rep}, A3 wall edges: {n_wall}",
                 fontsize=FONT_TITLE - 2, pad=15)
    ax.axis("off"); plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Panel_{out_tag}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Panel_{out_tag}")


# =============================================================================
# COMPOSITE
# =============================================================================

def assemble_composite():
    banner("[COMPOSITE] Figure 4")
    panels = {
        "a": "Panel_a_UMAP.png",
        "b": "Panel_b_SBS2_VS_NORMAL_chains.png",
        "btop": "Panel_b_top_SBS2_A3A_activator.png",
        "bbot": "Panel_b_bottom_SBS2_A3A_inhibitor.png",
        "c": "Panel_c_CNV_VS_NORMAL_chains.png",
        "ctop": "Panel_c_top_CNV_A3A_activator.png",
        "cbot": "Panel_c_bottom_CNV_A3A_inhibitor.png",
        "d": "Panel_d_CNV_A3B_chains.png",
    }
    imgs = {k: os.path.join(FIG_DIR, v) for k, v in panels.items()}
    fig = plt.figure(figsize=(38, 56))
    gs = gridspec.GridSpec(3, 2, hspace=0.06, wspace=0.04,
                           height_ratios=[1, 1, 1.1])
    placements = [
        ("btop", gs[0, 0], "b"), ("ctop", gs[0, 1], "c"),
        ("bbot", gs[1, 0], ""),  ("cbot", gs[1, 1], ""),
        ("a",    gs[2, 0], "a"), ("d",    gs[2, 1], "d"),
    ]
    for key, cell, letter in placements:
        if not os.path.exists(imgs[key]):
            log(f"  [SKIP] missing {panels[key]}"); continue
        ax = fig.add_subplot(cell)
        ax.imshow(plt.imread(imgs[key])); ax.axis("off")
        if letter:
            ax.set_title(letter, fontsize=40, fontweight="bold", loc="left", pad=8)
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Figure4_composite.{ext}"), dpi=200, bbox_inches="tight")
    plt.close()
    log("  [SAVE] Figure4_composite")

# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("GENERATE FIGURE 4 PANELS (v2)")
    harris = load_harris_interactors()

    log("\nLoading adata for UMAP...")
    adata = sc.read_h5ad(ADATA_PATH)
    groups = pd.read_csv(GROUP_PATH, sep="\t")
    generate_umap_panel(adata, groups)
    del adata

    sbs2 = load_network_data(SBS2_DIR, "SBS2_VS_NORMAL", "SBS2-HIGH vs NORMAL", "Mutagenic entry")
    cnv = load_network_data(CNV_DIR, "CNV_VS_NORMAL", "CNV-HIGH vs NORMAL", "Productive infection")

    generate_a3_zoom(sbs2, harris, focus_gene="APOBEC3A", mode="both",
                     out_tag="b_SBS2_VS_NORMAL_chains")
    generate_a3_zoom(sbs2, harris, focus_gene="APOBEC3A", mode="activator",
                     out_tag="b_top_SBS2_A3A_activator", title_suffix="Activating chain")
    generate_a3_zoom(sbs2, harris, focus_gene="APOBEC3A", mode="repressor",
                     out_tag="b_bottom_SBS2_A3A_inhibitor", title_suffix="Inhibiting chain")
    generate_a3_zoom(cnv, harris, focus_gene="APOBEC3A", mode="both",
                     out_tag="c_CNV_VS_NORMAL_chains")    
    generate_a3_zoom(cnv, harris, focus_gene="APOBEC3A", mode="activator",
                     out_tag="c_top_CNV_A3A_activator", title_suffix="Activating chain")
    generate_a3_zoom(cnv, harris, focus_gene="APOBEC3A", mode="repressor",
                     out_tag="c_bottom_CNV_A3A_inhibitor", title_suffix="Inhibiting chain")
    generate_a3_zoom(cnv, harris, focus_gene="APOBEC3B", mode="both",
                     out_tag="d_CNV_A3B_chains")

    assemble_composite()
    banner("FIGURE 4 COMPLETE")
    log(f"Output: {FIG_DIR}")
    for f in sorted(os.listdir(FIG_DIR)):
        if f.endswith(".png") and not os.path.isdir(os.path.join(FIG_DIR, f)):
            log(f"  {f}")


if __name__ == "__main__":
    main()
