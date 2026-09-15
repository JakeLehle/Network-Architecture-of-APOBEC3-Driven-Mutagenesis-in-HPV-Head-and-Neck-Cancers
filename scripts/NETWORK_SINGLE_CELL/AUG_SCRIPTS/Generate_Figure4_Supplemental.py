#!/usr/bin/env python3
"""
Generate_Figure4_Supplemental.py
=================================

Supplemental figure showing the differential-network construction, for both
networks (primarily SBS2-HIGH vs NORMAL).

Per network:
  Supp_<name>_triple_heatmap     HIGH / LOW / DIFF Spearman matrices
                                 (community-ordered; HIGH - LOW = DIFF)
  Supp_<name>_global_network     full differential network, exploded by gene
                                 group so it can be read
  Supp_<name>_community_C<k>      full render of each A3-containing community
                                 (CNV gives C1 = A3A and C0 = A3B; the C0 render
                                 is the full A3B-community view discussed as a
                                 possible panel-d alternative)

Inputs (per NETWORK_<name>/):
  03_correlation_networks/corr_matrices/SC_corr_{HIGH,LOW,DIFF}.pkl
  04_communities/SC_best_partition.csv, SC_G_comm.gpickle, SC_selected_parameters.txt
  DIAGNOSTIC_AUDIT/Figure4_Node_Sizing.tsv
  00_input/Harris_A3_interactors.txt

Style: 28-34 pt text, hex colors, PDF + PNG at 300 DPI.

Usage:
    conda run -n NETWORK python Generate_Figure4_Supplemental.py
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.collections import LineCollection
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
FIG_DIR = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS", "SUPPLEMENT")
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")
os.makedirs(FIG_DIR, exist_ok=True)

NETWORKS = [
    {"name": "SBS2_VS_NORMAL", "label": "SBS2-HIGH vs NORMAL",
     "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL")},
    {"name": "CNV_VS_NORMAL", "label": "CNV-HIGH vs NORMAL",
     "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL")},
]

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}

COLOR_UP_TUMOR  = "#fcad61"
COLOR_UP_NORMAL = "#aadda4"
COLOR_A3        = "#ed6a5a"
COLOR_A3_TEXT   = "#c0392b"
COLOR_HARRIS    = "#F6D155"
COLOR_EDGE_POS  = "#b22222"
COLOR_EDGE_NEG  = "#4682b4"
COLOR_INTER     = "#bbbbbb"

FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22
LABEL_BG_ALPHA = 0.42

INTER_SCALE = 6.0
INTRA_SCALE = 1.0

# Global network only: drop gene groups (satellites) with this many nodes or fewer
GLOBAL_DROP_GROUPS_LE = 6


def log(msg):
    print(msg, flush=True)

def banner(t, c="="):
    print(f"\n{c*80}\n  {t}\n{c*80}", flush=True)


# =============================================================================
# LABEL ENGINE
# =============================================================================

def repel_labels(pos, labels, data_range=12.0, n_iters=200, push_strength=0.05):
    label_pos = {n: np.array(pos[n], dtype=float) + np.array([0, 0.16]) for n in labels}
    char_w = data_range * 0.012
    for _ in range(n_iters):
        keys = list(label_pos.keys())
        for i, n1 in enumerate(keys):
            for n2 in keys[i + 1:]:
                p1, p2 = label_pos[n1], label_pos[n2]
                diff = p1 - p2; dist = np.linalg.norm(diff)
                min_d = char_w * max(len(str(labels[n1])), len(str(labels[n2]))) * 0.55
                if dist < min_d and dist > 1e-8:
                    push = push_strength * (min_d - dist) * diff / dist
                    label_pos[n1] += push; label_pos[n2] -= push
    return label_pos


def draw_labels_with_boxes(ax, pos, label_pos, labels, font_sizes, font_colors, char_w=0.06):
    for n in labels:
        col = font_colors.get(n, "#000000"); fs = font_sizes.get(n, 12)
        nxy = np.array(pos[n]); lxy = np.array(label_pos[n])
        offset = np.linalg.norm(lxy - nxy)
        kwargs = dict(fontsize=fs, fontweight="bold", color=col, ha="center", va="center",
                      zorder=9, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                                          alpha=LABEL_BG_ALPHA))
        if offset > char_w * 2.5:
            ax.annotate(labels[n], xy=nxy, xytext=lxy, **kwargs,
                        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.6, alpha=0.5))
        else:
            ax.annotate(labels[n], xy=lxy, **kwargs)


# =============================================================================
# LOADERS
# =============================================================================

def load_harris():
    if not os.path.exists(HARRIS_PATH):
        return set()
    df = pd.read_csv(HARRIS_PATH, sep="\t")
    col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
    return set(df[col].astype(str))


def load_net(cfg, harris):
    net_dir = cfg["dir"]
    net = {"name": cfg["name"], "label": cfg["label"], "harris": harris}
    part = pd.read_csv(os.path.join(net_dir, "04_communities/SC_best_partition.csv"))
    net["partition_df"] = part
    net["gene_to_comm"] = dict(zip(part["gene"], part["community"]))
    with open(os.path.join(net_dir, "04_communities/SC_G_comm.gpickle"), "rb") as f:
        net["G_comm"] = pickle.load(f)
    de = pd.read_csv(os.path.join(net_dir, "02_differential_expression/SC_diffexpr_stats.csv"))
    net["de_log2fc"] = dict(zip(de["gene"], de["log2FC"]))
    sizing = os.path.join(net_dir, "DIAGNOSTIC_AUDIT/Figure4_Node_Sizing.tsv")
    net["node_intra"] = {}
    if os.path.exists(sizing):
        s = pd.read_csv(sizing, sep="\t")
        net["node_intra"] = {r["gene_symbol"]: float(r["intra_score"]) for _, r in s.iterrows()}
    params = os.path.join(net_dir, "04_communities/SC_selected_parameters.txt")
    net["params"] = {}
    if os.path.exists(params):
        with open(params) as f:
            for line in f:
                if "=" in line:
                    k, v = line.strip().split("=", 1); net["params"][k] = v
    cdir = os.path.join(net_dir, "03_correlation_networks/corr_matrices")
    net["corr"] = {}
    for tag in ["HIGH", "LOW", "DIFF"]:
        p = os.path.join(cdir, f"SC_corr_{tag}.pkl")
        if os.path.exists(p):
            with open(p, "rb") as f:
                net["corr"][tag] = pickle.load(f)
        else:
            log(f"  [WARNING] missing {p}")
    return net


# =============================================================================
# TRIPLE HEATMAP (HIGH / LOW / DIFF)
# =============================================================================

def plot_triple_heatmap(net):
    banner(f"[SUPP] Triple heatmap: {net['label']}")
    corr = net["corr"]
    if not all(k in corr for k in ["HIGH", "LOW", "DIFF"]):
        log("  [SKIP] need HIGH, LOW and DIFF matrices"); return
    high = corr["HIGH"]
    g2c = net["gene_to_comm"]
    genes_in = [g for g in net["partition_df"]["gene"] if g in high.index]
    c2g = {}
    for g in genes_in:
        c2g.setdefault(g2c[g], []).append(g)

    # within-community clustering on HIGH co-expression
    clustered = {}
    for c, genes in c2g.items():
        if len(genes) <= 2:
            clustered[c] = genes; continue
        sub = np.nan_to_num(high.loc[genes, genes].values, nan=0.0)
        try:
            d = np.clip((1.0 - np.abs(sub) + 1.0 - np.abs(sub).T) / 2, 0, None)
            np.fill_diagonal(d, 0)
            Z = linkage(squareform(d, checks=False), method="average")
            clustered[c] = [genes[i] for i in leaves_list(Z)]
        except Exception:
            clustered[c] = genes

    # between-community ordering
    cids = sorted(clustered.keys())
    if len(cids) > 2:
        mc = np.zeros((len(cids), len(cids)))
        for i, ci in enumerate(cids):
            for j, cj in enumerate(cids):
                if i == j: continue
                gi, gj = clustered[ci], clustered[cj]
                if gi and gj:
                    mc[i, j] = np.nanmean(np.abs(high.loc[gi, gj].values))
        try:
            dc = np.clip((1 - mc + 1 - mc.T) / 2, 0, None); np.fill_diagonal(dc, 0)
            Zc = linkage(squareform(dc, checks=False), method="average")
            cids = [cids[i] for i in leaves_list(Zc)]
        except Exception:
            pass

    ordered, bounds = [], []
    for c in cids:
        s = len(ordered); ordered.extend(clustered[c]); bounds.append((s, len(ordered), c))
    log(f"  {len(ordered)} genes across {len(cids)} gene groups")

    fig, axes = plt.subplots(1, 3, figsize=(48, 14))
    panels = [("HIGH", "HIGH (tumor)"), ("LOW", "LOW (normal)"), ("DIFF", "DIFF (HIGH-LOW)")]
    for ax, (tag, plabel) in zip(axes, panels):
        m = corr[tag]
        gp = [g for g in ordered if g in m.index]
        data = m.loc[gp, gp].values
        vmin = np.nanpercentile(data, 1)
        vmax = np.nanpercentile(data, 99)
        if vmax <= vmin: vmax = 1.0
        im = ax.imshow(data, aspect="auto", interpolation="nearest",
                       cmap="plasma", vmin=vmin, vmax=vmax)
        # gene-group boundary lines
        for s, e, c in bounds:
            if s > 0:
                ax.axhline(s - 0.5, color="white", lw=2, alpha=0.8)
                ax.axvline(s - 0.5, color="white", lw=2, alpha=0.8)
        # gene-group size labels on the right
        for s, e, c in bounds:
            if e - s >= 5:
                ax.annotate(f"C{c} ({e - s})", xy=(len(gp) + 5, (s + e) / 2),
                            fontsize=20, fontweight="bold", va="center", ha="left",
                            annotation_clip=False)
        cb = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.12)
        cb.set_label(f"Spearman rho ({plabel})", fontsize=24, labelpad=12)
        cb.ax.tick_params(labelsize=22)
        ax.set_xlabel("Genes (clustered within gene groups)", fontsize=28, labelpad=12)
        ax.set_ylabel("Genes (clustered within gene groups)", fontsize=28, labelpad=12)
        ax.set_title(f"{net['label']} | {plabel}", fontsize=28, pad=15)
        ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout(w_pad=4)
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Supp_{net['name']}_triple_heatmap.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Supp_{net['name']}_triple_heatmap")


# =============================================================================
# GLOBAL EXPLODED NETWORK
# =============================================================================

def exploded_layout(G, g2c, seed=42):
    c2n = {}
    for n in G.nodes():
        c2n.setdefault(g2c.get(n, -1), []).append(n)
    ucomms = sorted(c2n.keys())
    CG = nx.Graph()
    for c in ucomms:
        CG.add_node(c)
    for u, v, d in G.edges(data=True):
        cu, cv = g2c.get(u, -1), g2c.get(v, -1)
        if cu == cv: continue
        w = abs(float(d.get("abs_weight", d.get("weight", 1))))
        if CG.has_edge(cu, cv): CG[cu][cv]["weight"] += w
        else: CG.add_edge(cu, cv, weight=w)
    pos_c = nx.spring_layout(CG, seed=seed, weight="weight", iterations=200)
    for c in pos_c:
        pos_c[c] = INTER_SCALE * np.array(pos_c[c], dtype=float)
    rng = np.random.default_rng(seed); pos = {}
    for c, cn in c2n.items():
        sub = G.subgraph(cn).copy()
        if sub.number_of_nodes() == 1:
            ps = {cn[0]: np.array([0., 0.])}
        elif sub.number_of_edges() > 0:
            ps = nx.spring_layout(sub, seed=seed, weight="abs_weight",
                                  k=2. / np.sqrt(max(len(cn), 1)), iterations=200)
        else:
            ps = {n: rng.normal(0, 0.05, size=2) for n in cn}
        arr = np.array([ps[n] for n in cn], dtype=float)
        arr -= arr.mean(axis=0, keepdims=True); arr *= INTRA_SCALE
        ctr = np.array(pos_c.get(c, [0., 0.]), dtype=float)
        for i, n in enumerate(cn):
            pos[n] = ctr + arr[i]
    return pos, c2n, ucomms


def plot_global_network(net):
    banner(f"[SUPP] Global network: {net['label']}")
    G = net["G_comm"]; g2c = net["gene_to_comm"]; harris = net["harris"]
    node_intra = net["node_intra"]
    if G.number_of_nodes() == 0:
        log("  [SKIP] empty graph"); return
    # drop small gene groups (satellites); always keep any A3-containing group
    sizes = {}
    for n in G.nodes():
        c = g2c.get(n, -1); sizes[c] = sizes.get(c, 0) + 1
    a3_comms = {g2c.get(g) for g in A3_SYMBOLS if g in g2c}
    keep_comms = {c for c, s in sizes.items() if s > GLOBAL_DROP_GROUPS_LE} | (a3_comms & set(sizes))
    keep_nodes = [n for n in G.nodes() if g2c.get(n, -1) in keep_comms]
    log(f"  gene groups: {len(sizes)} total, dropping {len(sizes) - len(keep_comms)} "
        f"with <= {GLOBAL_DROP_GROUPS_LE} nodes, keeping {len(keep_comms)}")
    G = G.subgraph(keep_nodes).copy()
    if G.number_of_nodes() == 0:
        log("  [SKIP] nothing left after filtering"); return
    pos, c2n, ucomms = exploded_layout(G, g2c, seed=42)
    cmap = plt.get_cmap("tab20", max(len(ucomms), 1))
    ccol = {c: cmap(i % 20) for i, c in enumerate(ucomms)}

    fig, ax = plt.subplots(figsize=(26, 24))
    segs_pos, segs_neg, segs_inter = [], [], []
    for u, v, d in G.edges(data=True):
        if g2c.get(u) != g2c.get(v):
            segs_inter.append([pos[u], pos[v]]); continue
        (segs_pos if d.get("weight", 0) > 0 else segs_neg).append([pos[u], pos[v]])
    ax.add_collection(LineCollection(segs_inter, colors=COLOR_INTER, linewidths=0.2,
                                     alpha=0.05, zorder=0))
    ax.add_collection(LineCollection(segs_neg, colors=COLOR_EDGE_NEG, linewidths=0.3,
                                     alpha=0.10, zorder=1))
    ax.add_collection(LineCollection(segs_pos, colors=COLOR_EDGE_POS, linewidths=0.4,
                                     alpha=0.15, zorder=1))

    def gi(n): return node_intra.get(n, 0.1)
    for c in ucomms:
        cn = c2n.get(c, [])
        if not cn: continue
        ax.scatter([pos[n][0] for n in cn], [pos[n][1] for n in cn],
                   s=[30 + 200 * gi(n) for n in cn], color=[ccol[c]] * len(cn),
                   alpha=0.75, edgecolors="#000000", linewidths=0.2, zorder=2, rasterized=True)
    hh = [n for n in G.nodes() if n in harris and n not in A3_SYMBOLS]
    if hh:
        ax.scatter([pos[n][0] for n in hh], [pos[n][1] for n in hh], s=200,
                   facecolors="none", edgecolors=COLOR_HARRIS, linewidths=2.5, zorder=4)
    a3 = [n for n in G.nodes() if n in A3_SYMBOLS]
    for n in a3:
        ax.scatter(pos[n][0], pos[n][1], s=340, c=COLOR_A3, edgecolors="#000000",
                   linewidths=1.5, zorder=5)

    deg = dict(G.degree())
    labels, fsd, fcd = {}, {}, {}
    for n in a3:
        labels[n] = A3_ALIAS.get(n, n); fsd[n] = FONT_AXIS - 6; fcd[n] = COLOR_A3_TEXT
    for c, cn in c2n.items():
        for n in sorted(cn, key=lambda x: -deg.get(x, 0))[:3]:
            if n not in labels:
                labels[n] = n; fsd[n] = 14; fcd[n] = "#333333"
    lp = repel_labels(pos, labels, data_range=20.0, n_iters=200, push_strength=0.05)
    draw_labels_with_boxes(ax, pos, lp, labels, fsd, fcd, char_w=0.08)

    thr = net["params"].get("DIFF_THRESHOLD", "?")
    le = [Patch(fc=ccol[c], ec="#000000", lw=0.4, label=f"C{c} ({len(c2n.get(c, []))})")
          for c in ucomms[:14]]
    le += [Patch(fc=COLOR_A3, ec="#000000", label="APOBEC3"),
           Patch(fc="white", ec=COLOR_HARRIS, lw=2.5, label="A3 interactor"),
           Patch(fc=COLOR_EDGE_POS, ec="none", label="Gained co-expression"),
           Patch(fc=COLOR_EDGE_NEG, ec="none", label="Lost co-expression")]
    ax.legend(handles=le, loc="upper left", fontsize=14, ncol=2, framealpha=0.9,
              title="Gene groups", title_fontsize=16)
    ax.set_title(f"{net['label']}: full differential network (gene groups > {GLOBAL_DROP_GROUPS_LE} nodes)\n"
                 f"{G.number_of_nodes()} genes, {G.number_of_edges()} edges, "
                 f"{len(ucomms)} gene groups shown | threshold={thr}",
                 fontsize=FONT_TITLE - 4, pad=15)
    ax.axis("off"); ax.autoscale(); plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Supp_{net['name']}_global_network.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Supp_{net['name']}_global_network")


# =============================================================================
# FULL A3-COMMUNITY RENDER
# =============================================================================

def plot_a3_community(net, comm):
    banner(f"[SUPP] A3 community C{comm}: {net['label']}")
    G = net["G_comm"]; g2c = net["gene_to_comm"]; de = net["de_log2fc"]
    node_intra = net["node_intra"]; harris = net["harris"]

    nodes = [g for g, c in g2c.items() if c == comm and g in G]
    Gs = G.subgraph(nodes).copy()
    Gs.remove_nodes_from([n for n in Gs.nodes() if Gs.degree(n) == 0])
    n_nodes, n_edges = Gs.number_of_nodes(), Gs.number_of_edges()
    log(f"  C{comm}: {n_nodes} nodes, {n_edges} edges")
    if n_nodes == 0:
        log("  [SKIP] empty"); return

    def is_up(g): return de.get(g, 0) > 0
    wall = sum(1 for u, v, d in Gs.edges(data=True)
               if d.get("weight", 0) < 0 and is_up(u) and is_up(v))
    wall_pct = 100 * wall / n_edges if n_edges else 0

    log("  computing layout...")
    pos = nx.spring_layout(Gs, seed=42, weight="abs_weight",
                           k=2.0 / np.sqrt(n_nodes), iterations=150, scale=6.0)

    fig, ax = plt.subplots(figsize=(24, 22))
    segs_neg, segs_pos = [], []
    for u, v, d in Gs.edges(data=True):
        (segs_pos if d.get("weight", 0) > 0 else segs_neg).append([pos[u], pos[v]])
    ax.add_collection(LineCollection(segs_neg, colors=COLOR_EDGE_NEG, linewidths=0.3,
                                     alpha=0.06, zorder=0))
    ax.add_collection(LineCollection(segs_pos, colors=COLOR_EDGE_POS, linewidths=0.5,
                                     alpha=0.20, zorder=1))

    def ns(g): return 8 + 55 * node_intra.get(g, 0.1)
    up = [n for n in Gs.nodes() if is_up(n) and n not in A3_SYMBOLS]
    dn = [n for n in Gs.nodes() if not is_up(n) and n not in A3_SYMBOLS]
    ax.scatter([pos[n][0] for n in up], [pos[n][1] for n in up], s=[ns(n) for n in up],
               c=COLOR_UP_TUMOR, alpha=0.85, edgecolors="none", zorder=2, rasterized=True)
    ax.scatter([pos[n][0] for n in dn], [pos[n][1] for n in dn], s=[ns(n) for n in dn],
               c=COLOR_UP_NORMAL, alpha=0.85, edgecolors="none", zorder=2, rasterized=True)
    hh = sorted((set(Gs.nodes()) & harris) - A3_SYMBOLS)
    if hh:
        ax.scatter([pos[n][0] for n in hh], [pos[n][1] for n in hh],
                   s=[max(220, 6 * ns(n)) for n in hh], facecolors="none",
                   edgecolors=COLOR_HARRIS, linewidths=3.0, zorder=4)
    a3 = [n for n in Gs.nodes() if n in A3_SYMBOLS]
    if a3:
        ax.scatter([pos[n][0] for n in a3], [pos[n][1] for n in a3], s=2400,
                   facecolors="none", edgecolors=COLOR_A3, linewidths=5.0, zorder=5)
        ax.scatter([pos[n][0] for n in a3], [pos[n][1] for n in a3], s=1200,
                   c=COLOR_A3, edgecolors="#000000", linewidths=2.5, zorder=5)

    deg = dict(Gs.degree())
    top = [g for g, _ in sorted(deg.items(), key=lambda x: -x[1]) if g not in A3_SYMBOLS][:15]
    labels, fsd, fcd = {}, {}, {}
    for n in a3:
        labels[n] = A3_ALIAS.get(n, n); fsd[n] = FONT_AXIS; fcd[n] = COLOR_A3_TEXT
    for n in hh:
        labels[n] = n; fsd[n] = 18; fcd[n] = "#5a4a00"
    for n in top:
        if n not in labels:
            labels[n] = n; fsd[n] = 16; fcd[n] = "#8b4513" if is_up(n) else "#2e6e2e"
    lp = repel_labels(pos, labels, data_range=12.0, n_iters=400, push_strength=0.06)
    draw_labels_with_boxes(ax, pos, lp, labels, fsd, fcd, char_w=0.06)

    legend = [Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3"),
              Patch(fc=COLOR_UP_TUMOR, ec="none", label="Up in tumor"),
              Patch(fc=COLOR_UP_NORMAL, ec="none", label="Up in normal"),
              Patch(fc="white", ec=COLOR_HARRIS, lw=3, label=f"A3 interactor ({len(hh)})"),
              Patch(fc=COLOR_EDGE_POS, ec="none", label="Gained co-expression"),
              Patch(fc=COLOR_EDGE_NEG, ec="none", label="Lost co-expression")]
    ax.legend(handles=legend, loc="upper left", fontsize=FONT_LEGEND - 2, framealpha=0.9)
    a3lbl = ", ".join(A3_ALIAS.get(n, n) for n in a3) or "none"
    ax.set_title(f"{net['label']}: Community {comm} ({a3lbl})\n"
                 f"{n_nodes} genes, {n_edges} edges, {wall_pct:.0f}% wall edges",
                 fontsize=FONT_TITLE - 2, pad=15)
    ax.axis("off"); ax.autoscale(); plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR, f"Supp_{net['name']}_community_C{comm}.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Supp_{net['name']}_community_C{comm}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("GENERATE FIGURE 4 SUPPLEMENT")
    harris = load_harris()
    log(f"Harris interactors: {len(harris)}")
    for cfg in NETWORKS:
        log(f"\nLoading {cfg['name']}...")
        net = load_net(cfg, harris)
        plot_triple_heatmap(net)
        plot_global_network(net)
        a3_comms = sorted({net["gene_to_comm"].get(g) for g in ["APOBEC3A", "APOBEC3B"]
                           if g in net["gene_to_comm"]})
        for c in a3_comms:
            plot_a3_community(net, c)
        del net  # free corr matrices before next network
    banner("SUPPLEMENT COMPLETE")
    log(f"Output: {FIG_DIR}")
    for f in sorted(os.listdir(FIG_DIR)):
        if f.endswith(".png"):
            log(f"  {f}")


if __name__ == "__main__":
    main()
