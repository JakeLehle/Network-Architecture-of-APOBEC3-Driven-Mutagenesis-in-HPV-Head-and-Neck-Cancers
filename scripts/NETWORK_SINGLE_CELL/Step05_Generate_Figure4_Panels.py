#!/usr/bin/env python3
"""
Step05_Generate_Figure4_Panels.py
==================================

Figure 4 — Step 05: Generate publication-quality panels matching Figure 2 style.

Panels:
  4a — UMAP: HIGH vs LOW basal cell selection (larger text) [from Step 00]
  4b — Dual heatmap: HIGH + DIFF correlation (community-sorted, plasma cmap)
  4c — Overlap matrix: SC vs TCGA bulk communities (hypergeometric)
  4d — Exploded community network (full) with gene highlighting
  4d_zooms — Per-community zoomed subplots
  4e — Known A3-Interactor enrichment

Gene highlighting on network plots:
  - A3 genes: full name (APOBEC3A), red circle (#ed6a5a), red text
  - TCGA shared genes: orange circle (#f18f01), orange text
  - Known A3-Interactors: yellow circle (#fed766), dark gold text

Usage:
  conda run -n NETWORK python Step05_Generate_Figure4_Panels.py
"""

import os, sys, json, pickle
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from datetime import datetime

from network_config_SC import (
    DIR_01_GROUPS, DIR_02_DE, DIR_03_NETWORKS,
    DIR_04_COMMUNITIES, DIR_05_CENTRALITY, DIR_06_OVERLAP, FIGURE_4_PANELS,
    FIG2_COMMUNITIES, DIFF_THRESHOLD,
    HARRIS_ALL_PATH, HARRIS_A3B_PATH,
    A3_GENES, A3_ID_TO_ALIAS, A3_SYMBOL_TO_ALIAS, BIOMARKERS,
    OVERLAP_P_THRESHOLD, JACCARD_DISPLAY, COMMUNITY_BASE_SEED,
    banner, log, ensure_dir, load_ensg_to_symbol, convert_tcga_genes_to_symbols
)

# Colors
COLOR_A3          = "#ed6a5a"
COLOR_TCGA        = "#f18f01"
COLOR_HARRIS      = "#fed766"
COLOR_A3_TEXT     = "#c0392b"   # darker red for readability
COLOR_TCGA_TEXT   = "#d35400"   # darker orange for readability
COLOR_HARRIS_TEXT = "#b8860b"   # dark goldenrod for readability on white
INTER_SCALE  = 6.0
INTRA_SCALE  = 1.0

# Node size multipliers (base sizes scaled up to match increased font/edge)
NODE_BASE_GLOBAL = 120      # was 60
NODE_SCALE_GLOBAL = 800     # was 400
NODE_BASE_ZOOM = 160        # was 80
NODE_SCALE_ZOOM = 1000      # was 500

def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float); n = len(pvals)
    if n == 0: return np.array([])
    order = np.argsort(pvals); ranked = np.empty(n)
    ranked[order] = np.arange(1, n+1)
    adj = pvals * n / ranked
    adj = np.minimum.accumulate(adj[np.argsort(ranked)[::-1]])[::-1]
    return np.clip(adj, 0, 1)


# =============================================================================
# LABEL REPULSION — prevent overlapping gene names
# =============================================================================

def draw_labels_with_repulsion(G, pos, labels, ax, font_size=14, font_weight="bold",
                               color_map=None, default_color="black",
                               iterations=50, repel_force=0.02):
    """
    Draw gene labels with simple iterative repulsion to reduce overlap.

    color_map: dict mapping node -> text color (for class-colored labels)
    """
    if not labels:
        return

    # Get label positions (start at node positions)
    label_nodes = list(labels.keys())
    label_pos = {n: np.array(pos[n], dtype=float) for n in label_nodes if n in pos}

    if len(label_pos) < 2:
        # No repulsion needed for 0-1 labels
        for n in label_pos:
            col = color_map.get(n, default_color) if color_map else default_color
            ax.annotate(labels[n], label_pos[n], fontsize=font_size, fontweight=font_weight,
                        color=col, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
        return

    # Estimate label extents for repulsion (rough bbox in data coords)
    # Get the data range to scale font size to data units
    all_x = [pos[n][0] for n in G.nodes() if n in pos]
    all_y = [pos[n][1] for n in G.nodes() if n in pos]
    if not all_x:
        return
    data_range = max(max(all_x) - min(all_x), max(all_y) - min(all_y), 1e-6)
    char_width = data_range * 0.008 * (font_size / 12.0)

    # Offset labels slightly above node
    for n in label_pos:
        label_pos[n] = label_pos[n] + np.array([0, char_width * 1.5])

    # Iterative repulsion
    nodes_list = list(label_pos.keys())
    for _ in range(iterations):
        for i, n1 in enumerate(nodes_list):
            for j, n2 in enumerate(nodes_list):
                if i >= j:
                    continue
                p1 = label_pos[n1]
                p2 = label_pos[n2]
                diff = p1 - p2
                dist = np.linalg.norm(diff)

                # Estimate overlap threshold based on label lengths
                min_dist = char_width * max(len(labels[n1]), len(labels[n2])) * 0.5
                if dist < min_dist and dist > 1e-8:
                    push = repel_force * (min_dist - dist) * diff / dist
                    label_pos[n1] = label_pos[n1] + push
                    label_pos[n2] = label_pos[n2] - push

    # Draw labels with leader lines
    for n in label_nodes:
        if n not in label_pos or n not in pos:
            continue
        col = color_map.get(n, default_color) if color_map else default_color
        node_xy = pos[n]
        label_xy = label_pos[n]

        # Draw thin leader line if label moved away from node
        offset_dist = np.linalg.norm(label_xy - node_xy)
        if offset_dist > char_width * 2:
            ax.annotate(
                labels[n], xy=node_xy, xytext=label_xy,
                fontsize=font_size, fontweight=font_weight, color=col,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7),
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.5)
            )
        else:
            ax.annotate(
                labels[n], xy=label_xy,
                fontsize=font_size, fontweight=font_weight, color=col,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7)
            )


def build_color_map(labels, a3_nodes, tcga_nodes, harris_nodes):
    """Build a dict mapping labeled nodes to their class text color."""
    cmap = {}
    for n in labels:
        if n in a3_nodes:
            cmap[n] = COLOR_A3_TEXT
        elif n in tcga_nodes:
            cmap[n] = COLOR_TCGA_TEXT
        elif n in harris_nodes:
            cmap[n] = COLOR_HARRIS_TEXT
        else:
            cmap[n] = "black"
    return cmap


# ---- Loaders ----
def load_sc_communities():
    df = pd.read_csv(os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv"))
    comm_genes = {}
    for c in sorted(df["community"].unique()):
        comm_genes[c] = set(df[df["community"]==c]["gene"])
    return comm_genes, dict(zip(df["gene"], df["community"]))

def load_tcga_communities():
    ct = "TCGA-HNSC"
    path = os.path.join(FIG2_COMMUNITIES, ct, f"{ct}_best_partition.csv")
    if not os.path.exists(path):
        alt = os.path.join(FIG2_COMMUNITIES, f"{ct}_best_partition.csv")
        path = alt if os.path.exists(alt) else path
    if not os.path.exists(path): return {}
    df = pd.read_csv(path)
    gc = "gene" if "gene" in df.columns else df.columns[0]
    e2s = load_ensg_to_symbol()
    out = {}
    for c in sorted(df["community"].unique()):
        out[c] = convert_tcga_genes_to_symbols(set(df[df["community"]==c][gc]), e2s)
    return out

def load_harris():
    out = {"all": set(), "A3B_only": set()}
    for k, p in [("all", HARRIS_ALL_PATH), ("A3B_only", HARRIS_A3B_PATH)]:
        if os.path.exists(p):
            with open(p) as f:
                for line in f:
                    g = line.strip().split("\t")[0].strip()
                    if g and not g.startswith("#"): out[k].add(g)
    return out

def get_tcga_shared(sc_comms, tcga_comms):
    all_t = set(); [all_t.update(v) for v in tcga_comms.values()]
    all_s = set(); [all_s.update(v) for v in sc_comms.values()]
    return all_s & all_t

# ---- Panel 4b: Dual Heatmap ----
def plot_dual_heatmap(save_dir):
    banner("[PANEL 4b] Dual heatmap: HIGH + DIFF")
    corr_dir = os.path.join(DIR_03_NETWORKS, "corr_matrices")
    with open(os.path.join(corr_dir, "SC_corr_HIGH.pkl"), "rb") as f: corr_high = pickle.load(f)
    with open(os.path.join(corr_dir, "SC_corr_DIFF.pkl"), "rb") as f: corr_diff = pickle.load(f)
    part_df = pd.read_csv(os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv"))
    g2c = dict(zip(part_df["gene"], part_df["community"]))
    comm_genes_list = [g for g in part_df["gene"] if g in corr_high.index]
    c2g = {}
    for g in comm_genes_list: c2g.setdefault(g2c.get(g,-1),[]).append(g)

    clustered = {}
    for c, genes in c2g.items():
        if len(genes) <= 2: clustered[c] = genes; continue
        sub = np.nan_to_num(corr_high.loc[genes,genes].values, nan=0.0)
        try:
            dist = np.clip((1.0-np.abs(sub)+1.0-np.abs(sub).T)/2, 0, None)
            np.fill_diagonal(dist, 0)
            Z = linkage(squareform(dist, checks=False), method="average")
            clustered[c] = [genes[i] for i in leaves_list(Z)]
        except: clustered[c] = genes

    cids = sorted(clustered.keys())
    if len(cids) > 2:
        mc = np.zeros((len(cids),len(cids)))
        for i, ci in enumerate(cids):
            for j, cj in enumerate(cids):
                if i==j: continue
                gi = [g for g in clustered[ci] if g in corr_high.index]
                gj = [g for g in clustered[cj] if g in corr_high.index]
                if gi and gj: mc[i,j] = np.nanmean(np.abs(corr_high.loc[gi,gj].values))
        try:
            dc = np.clip((1-mc+1-mc.T)/2, 0, None); np.fill_diagonal(dc, 0)
            Zc = linkage(squareform(dc, checks=False), method="average")
            cids_sorted = [cids[i] for i in leaves_list(Zc)]
        except: cids_sorted = cids
    else: cids_sorted = cids

    ordered = []; boundaries = []
    for c in cids_sorted:
        s = len(ordered); ordered.extend(clustered[c]); boundaries.append((s, len(ordered), c))
    log(f"  Ordered {len(ordered)} genes across {len(cids_sorted)} communities")

    fig, axes = plt.subplots(1, 2, figsize=(32, 14))
    for ax, (label, mat) in zip(axes, [("HIGH (SBS2-active)", corr_high), ("DIFF (HIGH-LOW)", corr_diff)]):
        gp = [g for g in ordered if g in mat.index]
        data = mat.loc[gp,gp].values
        vmin, vmax = np.nanpercentile(data,1), np.nanpercentile(data,99)
        if vmax <= vmin: vmax = 1.0
        im = ax.imshow(data, aspect="auto", interpolation="nearest", cmap="plasma", vmin=vmin, vmax=vmax)
        for s,e,c in boundaries:
            if s > 0: ax.axhline(s-0.5, color="white", lw=2, alpha=0.8); ax.axvline(s-0.5, color="white", lw=2, alpha=0.8)
        for s,e,c in boundaries:
            if e-s >= 5: ax.annotate(f"C{c} ({e-s})", xy=(len(gp)+5, (s+e)/2), fontsize=26, fontweight="bold", va="center", ha="left", annotation_clip=False)
        cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.12)
        cbar.set_label(f"Spearman rho ({label})", fontsize=30, labelpad=12); cbar.ax.tick_params(labelsize=22)
        ax.set_xlabel("Genes (clustered within communities)", fontsize=34, labelpad=12)
        ax.set_ylabel("Genes (clustered within communities)", fontsize=34, labelpad=12)
        ax.set_title(f"SC Basal | {label}", fontsize=34, pad=15); ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout(w_pad=4)
    for ext in ["pdf","png"]:
        plt.savefig(os.path.join(save_dir, f"Panel_4b_dual_heatmap.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(); log(f"  [SAVE] Panel 4b heatmap")

# ---- Exploded Layout ----
def exploded_layout(G, g2c, seed=42):
    c2n = {}
    for n in G.nodes(): c2n.setdefault(g2c.get(n,-1),[]).append(n)
    ucomms = sorted(c2n.keys())
    CG = nx.Graph()
    for c in ucomms: CG.add_node(c)
    for u,v,d in G.edges(data=True):
        cu,cv = g2c.get(u,-1), g2c.get(v,-1)
        if cu==cv: continue
        w = abs(float(d.get("abs_weight", d.get("weight",1))))
        if CG.has_edge(cu,cv): CG[cu][cv]["weight"] += w
        else: CG.add_edge(cu,cv,weight=w)
    pos_c = nx.spring_layout(CG, seed=seed, weight="weight", iterations=200)
    for c in pos_c: pos_c[c] = INTER_SCALE * np.array(pos_c[c], dtype=float)
    rng = np.random.default_rng(seed); pos = {}
    for c, cn in c2n.items():
        sub = G.subgraph(cn).copy()
        if sub.number_of_nodes()==1: ps = {cn[0]: np.array([0.,0.])}
        elif sub.number_of_edges()>0: ps = nx.spring_layout(sub, seed=seed, weight="abs_weight", k=2./np.sqrt(max(len(cn),1)), iterations=200)
        else: ps = {n: rng.normal(0,0.05,size=2) for n in cn}
        arr = np.array([ps[n] for n in cn], dtype=float)
        arr -= arr.mean(axis=0, keepdims=True); arr *= INTRA_SCALE
        ctr = np.array(pos_c.get(c,[0.,0.]), dtype=float)
        for i,n in enumerate(cn): pos[n] = ctr + arr[i]
    return pos, c2n, ucomms

# ---- Panel 4d: Full Network ----
def plot_network(G, g2c, tcga_shared, harris_genes, save_dir):
    banner("[PANEL 4d] Full community network")
    if G.number_of_nodes()==0: return {},{},[]
    pos, c2n, ucomms = exploded_layout(G, g2c, seed=COMMUNITY_BASE_SEED)
    nc = len(ucomms); cmap = plt.cm.get_cmap("tab20", max(nc,1))
    ccmap = {c: cmap(i) for i,c in enumerate(ucomms)}
    fig, ax = plt.subplots(figsize=(26,24))

    # Edges
    intra = [(u,v,d) for u,v,d in G.edges(data=True) if g2c.get(u)==g2c.get(v)]
    inter = [(u,v,d) for u,v,d in G.edges(data=True) if g2c.get(u)!=g2c.get(v)]
    if intra:
        ec = ["firebrick" if d.get("weight",0)>0 else "steelblue" for u,v,d in intra]
        ew = [0.3+1.5*abs(d.get("weight",0)) for u,v,d in intra]
        nx.draw_networkx_edges(G, pos, edgelist=[(u,v) for u,v,d in intra], ax=ax, alpha=0.3, width=ew, edge_color=ec)
    if inter:
        nx.draw_networkx_edges(G, pos, edgelist=[(u,v) for u,v,d in inter], ax=ax, alpha=0.1, width=0.4, edge_color="gray", style="dashed")

    # Base nodes (scaled up)
    deg = dict(G.degree()); md = max(deg.values()) if deg else 1
    nodes = list(G.nodes())
    ns = [NODE_BASE_GLOBAL + NODE_SCALE_GLOBAL*(deg[n]/md) for n in nodes]
    nc_list = [ccmap.get(g2c.get(n,-1),"gray") for n in nodes]
    nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax, node_size=ns, node_color=nc_list, alpha=0.85, edgecolors="black", linewidths=0.5)

    # Highlight rings
    rs = 2.5
    a3_set = set(A3_GENES)
    harris_nn = [n for n in nodes if n in harris_genes]
    tcga_nn = [n for n in nodes if n in tcga_shared and n not in a3_set]
    a3_nn = [n for n in nodes if n in a3_set]

    if harris_nn:
        nx.draw_networkx_nodes(G, pos, nodelist=harris_nn, ax=ax, node_size=[rs*(NODE_BASE_GLOBAL+NODE_SCALE_GLOBAL*(deg[n]/md)) for n in harris_nn], node_color="none", edgecolors=COLOR_HARRIS, linewidths=5.0, alpha=0.9)
        log(f"  Known A3-Interactors in network: {len(harris_nn)} — {sorted(harris_nn)}")
    if tcga_nn:
        nx.draw_networkx_nodes(G, pos, nodelist=tcga_nn, ax=ax, node_size=[rs*(NODE_BASE_GLOBAL+NODE_SCALE_GLOBAL*(deg[n]/md)) for n in tcga_nn], node_color="none", edgecolors=COLOR_TCGA, linewidths=5.0, alpha=0.9)
        log(f"  TCGA shared in network: {len(tcga_nn)} — {sorted(tcga_nn)}")
    if a3_nn:
        nx.draw_networkx_nodes(G, pos, nodelist=a3_nn, ax=ax, node_size=[rs*(NODE_BASE_GLOBAL+NODE_SCALE_GLOBAL*(deg[n]/md)) for n in a3_nn], node_color="none", edgecolors=COLOR_A3, linewidths=5.0, alpha=1.0)
        log(f"  A3 in network: {[f'{n} (C{g2c.get(n)})' for n in a3_nn]}")

    # Build labels
    labels = {}
    for n in a3_nn: labels[n] = n
    for n in tcga_nn: labels[n] = n
    for n in harris_nn: labels[n] = n
    for c, cn in c2n.items():
        for n in sorted(cn, key=lambda x: deg.get(x,0), reverse=True)[:3]:
            if n not in labels: labels[n] = n

    # Class-colored text with repulsion
    color_map = build_color_map(labels, set(a3_nn), set(tcga_nn), set(harris_nn))
    draw_labels_with_repulsion(G, pos, labels, ax, font_size=32, font_weight="bold",
                               color_map=color_map, iterations=80, repel_force=0.025)

    # Legend
    le = [Patch(facecolor=ccmap[c], edgecolor="black", lw=0.5, label=f"C{c} ({len(c2n.get(c,[]))})") for c in ucomms[:14]]
    le += [Patch(fc="none", ec=COLOR_A3, lw=3, label="APOBEC3"),
           Patch(fc="none", ec=COLOR_TCGA, lw=3, label="TCGA shared"),
           Patch(fc="none", ec=COLOR_HARRIS, lw=3, label="Known A3-Interactor")]
    ax.legend(handles=le, loc="upper left", fontsize=16, ncol=2, framealpha=0.9, title="Communities", title_fontsize=18)
    ax.set_title(f"SC Basal DIFF Network — 14 Communities (|Δρ| ≥ {DIFF_THRESHOLD})", fontsize=26, pad=15)
    ax.axis("off"); plt.tight_layout()
    for ext in ["pdf","png"]: plt.savefig(os.path.join(save_dir, f"Panel_4d_network_full.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(); log(f"  [SAVE] Panel 4d full network")
    return pos, ccmap, ucomms

# ---- Per-community zooms ----
def plot_zooms(G, g2c, ccmap, ucomms, tcga_shared, harris_genes, save_dir):
    banner("[PANEL 4d] Per-community zooms")
    zdir = ensure_dir(os.path.join(save_dir, "community_zooms"))
    c2n = {}
    for n in G.nodes(): c2n.setdefault(g2c.get(n,-1),[]).append(n)
    a3_set = set(A3_GENES)

    for c in ucomms:
        cn = c2n.get(c,[])
        if len(cn) < 3: continue
        Gs = G.subgraph(cn).copy()
        if Gs.number_of_edges()==0: continue
        ps = nx.spring_layout(Gs, seed=COMMUNITY_BASE_SEED, weight="abs_weight", k=2./np.sqrt(max(len(cn),1)), iterations=200)
        fig, ax = plt.subplots(figsize=(16,14))
        ec = ["firebrick" if d.get("weight",0)>0 else "steelblue" for u,v,d in Gs.edges(data=True)]
        ew = [0.5+2.5*abs(d.get("weight",0)) for u,v,d in Gs.edges(data=True)]
        nx.draw_networkx_edges(Gs, ps, ax=ax, edge_color=ec, width=ew, alpha=0.35)

        deg = dict(Gs.degree()); md = max(deg.values()) if deg else 1
        ns = [NODE_BASE_ZOOM + NODE_SCALE_ZOOM*(deg[n]/md) for n in cn]
        nx.draw_networkx_nodes(Gs, ps, ax=ax, node_size=ns, node_color=[ccmap.get(c,"gray")]*len(cn), alpha=0.85, edgecolors="black", linewidths=0.8)

        rs = 2.5; sn = list(Gs.nodes())
        hh = [n for n in sn if n in harris_genes]
        if hh: nx.draw_networkx_nodes(Gs, ps, nodelist=hh, ax=ax, node_size=[rs*(NODE_BASE_ZOOM+NODE_SCALE_ZOOM*(deg[n]/md)) for n in hh], node_color="none", edgecolors=COLOR_HARRIS, linewidths=4.0)
        tt = [n for n in sn if n in tcga_shared and n not in a3_set]
        if tt: nx.draw_networkx_nodes(Gs, ps, nodelist=tt, ax=ax, node_size=[rs*(NODE_BASE_ZOOM+NODE_SCALE_ZOOM*(deg[n]/md)) for n in tt], node_color="none", edgecolors=COLOR_TCGA, linewidths=4.0)
        aa = [n for n in sn if n in a3_set]
        if aa: nx.draw_networkx_nodes(Gs, ps, nodelist=aa, ax=ax, node_size=[rs*(NODE_BASE_ZOOM+NODE_SCALE_ZOOM*(deg[n]/md)) for n in aa], node_color="none", edgecolors=COLOR_A3, linewidths=5.0)

        if len(cn) <= 50:
            labels = {n:n for n in sn}; fs = max(8, 14-len(cn)//10)
        else:
            top = sorted(deg.items(), key=lambda x:x[1], reverse=True)[:20]
            labels = {n:n for n,d in top}
            for n in aa+tt+hh: labels[n] = n
            fs = 30

        color_map = build_color_map(labels, set(aa), set(tt), set(hh))
        draw_labels_with_repulsion(Gs, ps, labels, ax, font_size=fs, font_weight="bold",
                                   color_map=color_map, iterations=60, repel_force=0.03)

        spec = [A3_SYMBOL_TO_ALIAS.get(n,n) for n in aa]
        extra = f" — {', '.join(spec)}" if spec else ""
        ax.set_title(f"Community {c} — {len(cn)} genes, {Gs.number_of_edges()} edges{extra}", fontsize=22, pad=10)
        ax.axis("off"); plt.tight_layout()
        for ext in ["png","pdf"]: plt.savefig(os.path.join(zdir, f"SC_community_{c:02d}.{ext}"), dpi=300, bbox_inches="tight")
        plt.close(); log(f"  [SAVE] Community {c}")
    log(f"  [SAVE] All zooms -> {zdir}")

# ---- Overlap ----
def compute_overlap(sc_c, tc_c, univ):
    si = sorted(sc_c.keys()); ti = sorted(tc_c.keys())
    oc = pd.DataFrame(0, index=si, columns=ti)
    pv = pd.DataFrame(1., index=si, columns=ti)
    jc = pd.DataFrame(0., index=si, columns=ti); det = []
    for s in si:
        for t in ti:
            o = sc_c[s] & tc_c[t]; k = len(o)
            oc.loc[s,t] = k
            pv.loc[s,t] = hypergeom.sf(k-1, univ, len(tc_c[t]), len(sc_c[s])) if k>0 else 1.
            u = len(sc_c[s]|tc_c[t]); jc.loc[s,t] = k/u if u>0 else 0
            if k>0: det.append({"sc_community":s,"tcga_community":t,"overlap_count":k,"jaccard":k/u if u>0 else 0,"p_value":pv.loc[s,t],"genes":";".join(sorted(o))})
    return oc, pv, jc, pd.DataFrame(det)

def plot_overlap(oc, ap, jc, sd):
    banner("[PANEL 4c] Overlap Heatmap")
    lp = -np.log10(ap.clip(lower=1e-10)).clip(upper=10)
    fig, ax = plt.subplots(figsize=(max(10,len(ap.columns)*0.9+3), max(8,len(ap.index)*0.7+3)))
    im = ax.imshow(lp.values, cmap="YlOrRd", aspect="auto", vmin=0, vmax=max(3,lp.values.max()))
    for i in range(len(ap.index)):
        for j in range(len(ap.columns)):
            cnt = int(oc.iloc[i,j])
            if cnt > 0:
                col = "white" if lp.iloc[i,j]>2 else "black"
                txt = f"{cnt}\n({jc.iloc[i,j]:.2f})" if JACCARD_DISPLAY else str(cnt)
                ax.text(j,i,txt, ha="center", va="center", color=col, fontsize=9, fontweight="bold" if ap.iloc[i,j]<0.05 else "normal")
    ax.set_xticks(range(len(ap.columns))); ax.set_xticklabels([f"TCGA C{c}" for c in ap.columns], rotation=45, ha="right", fontsize=12)
    ax.set_yticks(range(len(ap.index))); ax.set_yticklabels([f"SC C{c}" for c in ap.index], fontsize=12)
    ax.set_xlabel("TCGA Bulk Communities (Figure 2)", fontsize=16); ax.set_ylabel("SC Communities (Figure 4)", fontsize=16)
    ax.set_title("Community Overlap: SC vs TCGA Bulk", fontsize=18)
    plt.colorbar(im, ax=ax, shrink=0.6, label="-log10(BH adj.p)"); plt.tight_layout()
    for ext in ["pdf","png"]: plt.savefig(os.path.join(sd, f"Panel_4c_Community_Overlap.{ext}"), dpi=300, bbox_inches="tight")
    plt.close()

# ---- Known A3-Interactor Enrichment ----
def harris_enrich(sc_c, hs, label, sd, od):
    banner(f"[Known A3-Interactor Enrichment] {label}")
    all_sc = set(); [all_sc.update(v) for v in sc_c.values()]
    hid = hs & all_sc; nh = len(hid); univ = len(all_sc)
    log(f"  {label}: {nh}/{len(hs)} in network")
    rows = []
    for c in sorted(sc_c.keys()):
        o = sc_c[c] & hid; k = len(o)
        pv = hypergeom.sf(k-1, univ, nh, len(sc_c[c])) if k>0 and nh>0 else 1.
        rows.append({"community":c,"size":len(sc_c[c]),"known_interactors":k,"p_value":pv,"genes":";".join(sorted(o)) if o else ""})
    df = pd.DataFrame(rows)
    if len(df): df["fdr"] = bh_fdr(df["p_value"].values)
    df.to_csv(os.path.join(od, f"known_A3_interactor_enrichment_{label}.csv"), index=False)
    for _, r in df.iterrows():
        if r["known_interactors"]>0: log(f"    C{int(r['community'])}: {int(r['known_interactors'])} p={r['p_value']:.4f} FDR={r['fdr']:.4f} {r['genes']}")
    if label == "A3B_only":
        fig, ax = plt.subplots(figsize=(max(8,len(df)*0.7),5))
        ax.bar(range(len(df)), df["known_interactors"], color=["firebrick" if p<0.05 else "steelblue" for p in df["fdr"]], edgecolor="black", lw=0.5)
        ax.set_xticks(range(len(df))); ax.set_xticklabels([f"C{c}" for c in df["community"]], fontsize=12)
        ax.set_xlabel("SC Community", fontsize=14); ax.set_ylabel("Known A3-Interactors", fontsize=14)
        ax.set_title("Known A3-Interactor Enrichment (red=FDR<0.05)", fontsize=16); plt.tight_layout()
        for ext in ["pdf","png"]: plt.savefig(os.path.join(sd, f"Panel_4e_A3_Interactor_Enrichment.{ext}"), dpi=300, bbox_inches="tight")
        plt.close()

# ---- Main ----
def main():
    t0 = datetime.now()
    banner("[STEP 05] Generate Figure 4 Panels"); log(f"Start: {t0}")
    ensure_dir(DIR_06_OVERLAP); ensure_dir(FIGURE_4_PANELS)
    sc_c, g2c = load_sc_communities()
    tc_c = load_tcga_communities()
    harris = load_harris(); all_harris = harris["all"] | harris["A3B_only"]
    tcga_sh = get_tcga_shared(sc_c, tc_c)
    log(f"  SC: {len(sc_c)} comms, {sum(len(v) for v in sc_c.values())} genes")
    log(f"  TCGA: {len(tc_c)} comms | Shared: {len(tcga_sh)} | Known A3-Interactors: {len(harris['all'])}")

    plot_dual_heatmap(FIGURE_4_PANELS)

    if tc_c:
        all_s = set(); [all_s.update(v) for v in sc_c.values()]
        all_t = set(); [all_t.update(v) for v in tc_c.values()]
        univ = len(all_s | all_t)
        oc, pv, jc, det = compute_overlap(sc_c, tc_c, univ)
        ap = pd.DataFrame(bh_fdr(pv.values.flatten()).reshape(pv.shape), index=pv.index, columns=pv.columns)
        oc.to_csv(os.path.join(DIR_06_OVERLAP,"overlap_matrix.csv"))
        ap.to_csv(os.path.join(DIR_06_OVERLAP,"hypergeom_pvalues_BH.csv"))
        jc.to_csv(os.path.join(DIR_06_OVERLAP,"jaccard_matrix.csv"))
        if len(det): det["fdr"]=bh_fdr(det["p_value"].values); det.to_csv(os.path.join(DIR_06_OVERLAP,"overlap_gene_details.csv"),index=False)
        sig = det[det["fdr"]<0.05] if len(det) and "fdr" in det.columns else pd.DataFrame()
        log(f"  Significant overlaps: {len(sig)}")
        plot_overlap(oc, ap, jc, FIGURE_4_PANELS)

    gpath = os.path.join(DIR_04_COMMUNITIES, "SC_G_comm.gpickle")
    if os.path.exists(gpath):
        with open(gpath,"rb") as f: Gc = pickle.load(f)
        pos, ccm, uc = plot_network(Gc, g2c, tcga_sh, all_harris, FIGURE_4_PANELS)
        if pos: plot_zooms(Gc, g2c, ccm, uc, tcga_sh, all_harris, FIGURE_4_PANELS)

    for label, genes in harris.items():
        if genes: harris_enrich(sc_c, genes, label, FIGURE_4_PANELS, DIR_06_OVERLAP)

    banner(f"[STEP 05 COMPLETE] Elapsed: {datetime.now()-t0}")

if __name__ == "__main__":
    main()
