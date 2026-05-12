#!/usr/bin/env python3
"""
Generate_Figure4_Panels.py
===========================

Generate publication-quality Figure 4 panels from SC network pipeline outputs.
Supersedes Step05_Generate_Figure4_Panels.py with Figure 2 styling:
  - intra_score-based node/label sizing (from Figure4_Node_Sizing.tsv)
  - Label repulsion algorithm with background boxes
  - Three-class gene highlighting with colored rings
  - Doubled edge thickness on zoom panels
  - Community center annotations

Panels:
  4a -- UMAP: HIGH vs A3+CNV-matched LOW cell selection
  4b -- Dual heatmap: HIGH + DIFF correlation (community-sorted, plasma cmap)
  4c -- SC vs TCGA community overlap matrix (hypergeometric)
  4d -- Exploded community network (full) with three-class highlighting
  4e -- Harris A3 interactor enrichment bar chart
  Supplement -- Methods graphical abstract (HIGH -> LOW -> DIFF -> network)
  Community zooms -- All communities with intra_score sizing + gene class rings

Gene highlighting (ring priority: A3 > TCGA shared > Harris interactors):
  - A3 genes:          coral fill (#ed6a5a), red ring, red text
  - TCGA shared genes: orange ring (#f18f01), orange text
  - Harris interactors: gold ring (#fed766), gold text

Run AFTER the full pipeline (Steps 00B-04) and Compute_Node_Importance_Scores_SC.py.

Usage:
    conda run -n NETWORK python Generate_Figure4_Panels.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os, sys, json, pickle
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
from datetime import datetime

from network_config_SC import (
    DIR_00_INPUT, DIR_01_GROUPS, DIR_02_DE, DIR_03_NETWORKS,
    DIR_04_COMMUNITIES, DIR_05_CENTRALITY, DIR_06_OVERLAP, FIGURE_4_PANELS,
    FIG2_COMMUNITIES, FIG2_CLEANED, FIG4_ROOT,
    HARRIS_ALL_PATH, HARRIS_A3B_PATH,
    ADATA_FINAL_PATH,
    A3_GENES_SYMBOLS, A3_SYMBOL_TO_ALIAS,
    OVERLAP_P_THRESHOLD, COMMUNITY_BASE_SEED,
    TARGET_CELL_TYPE,
    banner, log, ensure_dir, load_ensg_to_symbol, convert_tcga_genes_to_symbols
)


# =============================================================================
# CONSTANTS
# =============================================================================

# Colors (all hex)
COLOR_HIGH     = "#ed6a5a"   # coral -- SBS2-HIGH cells
COLOR_LOW      = "#5b8e7d"   # muted green -- A3+CNV-matched LOW
COLOR_GRAY     = "#e0e0e0"   # other cells

COLOR_A3       = "#ed6a5a"   # coral -- APOBEC3 genes
COLOR_A3_TEXT  = "#c0392b"   # dark red
COLOR_TCGA     = "#f18f01"   # orange -- TCGA shared genes
COLOR_HARRIS   = "#fed766"   # gold -- Harris interactors
COLOR_EDGE_POS = "#b22222"   # firebrick -- gained co-expression
COLOR_EDGE_NEG = "#4682b4"   # steelblue -- lost co-expression

# Layout
INTER_SCALE = 6.0
INTRA_SCALE = 1.0

# Font sizes (28-34 range for publication)
FONT_TITLE  = 32
FONT_AXIS   = 30
FONT_LEGEND = 24
FONT_TICK   = 22

# Label background alpha
LABEL_BG_ALPHA = 0.42

# Output
FIG_DIR = ensure_dir(FIGURE_4_PANELS)

A3_SET = set(A3_GENES_SYMBOLS)


# =============================================================================
# BH FDR
# =============================================================================

def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = (np.arange(n) + 1) * pvals[order] / n
    for i in range(n - 2, -1, -1):
        ranked[order[i]] = min(ranked[order[i]], ranked[order[i + 1]])
    return np.clip(ranked, 0, 1)


# =============================================================================
# LABEL REPULSION (ported from Figure 2)
# =============================================================================

def repel_labels(pos, labels, data_range=12.0, n_iters=120, push_strength=0.04):
    """Iteratively push overlapping labels apart. Returns adjusted positions."""
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
    """Draw labels with semi-transparent white background boxes."""
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
# COMMUNITY-ORDERED GENE LIST
# =============================================================================

def build_community_ordered_genes(partition_df, ref_matrix, gene_to_comm):
    """Order genes by community with hierarchical clustering within and between."""
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
                gi = [g for g in clustered[comm_ids[i]] if g in ref_matrix.index]
                gj = [g for g in clustered[comm_ids[j]] if g in ref_matrix.index]
                if gi and gj:
                    cd[i, j] = cd[j, i] = 1.0 - np.nanmean(
                        np.abs(ref_matrix.loc[gi, gj].values))
        try:
            Z_c = linkage(squareform(cd, checks=False), method="average")
            comm_ids_sorted = [comm_ids[i] for i in leaves_list(Z_c)]
        except Exception:
            comm_ids_sorted = sorted(comm_ids,
                                     key=lambda c: len(clustered[c]), reverse=True)
    else:
        comm_ids_sorted = sorted(comm_ids,
                                 key=lambda c: len(clustered[c]), reverse=True)

    ordered_genes = []
    comm_boundaries = []
    for c in comm_ids_sorted:
        start = len(ordered_genes)
        ordered_genes.extend(clustered[c])
        comm_boundaries.append((start, len(ordered_genes), c))
    return ordered_genes, comm_boundaries, comm_ids_sorted


# =============================================================================
# LOAD ALL DATA
# =============================================================================

def load_all():
    """Load all pipeline outputs needed for panel generation."""
    banner("[INIT] Load all data")
    data = {}

    # Partition
    part_path = os.path.join(DIR_04_COMMUNITIES, "SC_best_partition.csv")
    data['partition_df'] = pd.read_csv(part_path)
    data['gene_to_comm'] = dict(zip(data['partition_df']['gene'],
                                     data['partition_df']['community']))
    log(f"  Partition: {len(data['partition_df'])} genes, "
        f"{data['partition_df']['community'].nunique()} communities")

    # Community graph
    graph_path = os.path.join(DIR_04_COMMUNITIES, "SC_G_comm.gpickle")
    with open(graph_path, "rb") as f:
        data['G_comm'] = pickle.load(f)
    log(f"  Community graph: {data['G_comm'].number_of_nodes()} nodes, "
        f"{data['G_comm'].number_of_edges()} edges")

    # Correlation matrices
    corr_dir = os.path.join(DIR_03_NETWORKS, "corr_matrices")
    for name in ['HIGH', 'LOW', 'DIFF']:
        path = os.path.join(corr_dir, f"SC_corr_{name}.pkl")
        if os.path.exists(path):
            with open(path, "rb") as f:
                data[f'corr_{name.lower()}'] = pickle.load(f)
            log(f"  {name} corr: {data[f'corr_{name.lower()}'].shape}")

    # Node importance scores
    sizing_path = os.path.join(FIG4_ROOT, "DIAGNOSTIC_AUDIT", "Figure4_Node_Sizing.tsv")
    data['node_intra'] = {}
    data['node_inter'] = {}
    if os.path.exists(sizing_path):
        sizing_df = pd.read_csv(sizing_path, sep="\t")
        for _, row in sizing_df.iterrows():
            data['node_intra'][row['gene_symbol']] = float(row['intra_score'])
            data['node_inter'][row['gene_symbol']] = float(row['inter_score'])
        log(f"  Node scores: {len(sizing_df)} genes")
    else:
        log(f"  WARNING: Figure4_Node_Sizing.tsv not found, falling back to degree")

    # Harris interactors
    data['harris_all'] = set()
    data['harris_a3b'] = set()
    if os.path.exists(HARRIS_ALL_PATH):
        with open(HARRIS_ALL_PATH) as f:
            data['harris_all'] = set(line.strip() for line in f if line.strip())
    if os.path.exists(HARRIS_A3B_PATH):
        with open(HARRIS_A3B_PATH) as f:
            data['harris_a3b'] = set(line.strip() for line in f if line.strip())
    log(f"  Harris: {len(data['harris_all'])} all, {len(data['harris_a3b'])} A3B-specific")

    # TCGA shared genes (from overlap analysis)
    data['tcga_shared'] = set()
    overlap_path = os.path.join(DIR_06_OVERLAP, "overlap_gene_details.csv")
    if os.path.exists(overlap_path):
        odf = pd.read_csv(overlap_path)
        col = 'gene' if 'gene' in odf.columns else 'gene_symbol'
        if col in odf.columns:
            data['tcga_shared'] = set(odf[col].values)
    log(f"  TCGA shared: {len(data['tcga_shared'])} genes")

    # TCGA bulk communities (for overlap panel)
    data['tcga_communities'] = {}
    try:
        tcga_comms = convert_tcga_genes_to_symbols(load_tcga_bulk_communities())
        data['tcga_communities'] = tcga_comms
    except Exception:
        try:
            data['tcga_communities'] = load_tcga_bulk_communities_fallback()
        except Exception:
            log(f"  TCGA bulk communities: not loaded")

    # Group assignments
    groups_path = os.path.join(DIR_01_GROUPS, "SC_Basal_group_assignments.tsv")
    groups = pd.read_csv(groups_path, sep='\t')
    data['high_cells'] = groups[groups['group'] == 'HIGH']['cell_barcode'].tolist()
    data['low_cells']  = groups[groups['group'] == 'LOW']['cell_barcode'].tolist()
    log(f"  Groups: {len(data['high_cells'])} HIGH, {len(data['low_cells'])} LOW")

    # Build community data structures
    G = data['G_comm']
    g2c = data['gene_to_comm']
    data['comm_to_nodes'] = {}
    for n in G.nodes():
        c = g2c.get(n, -1)
        data['comm_to_nodes'].setdefault(c, []).append(n)
    data['unique_comms'] = sorted(data['comm_to_nodes'].keys())
    cmap = plt.colormaps["tab20"]
    data['comm_colors'] = {c: mcolors.to_hex(cmap(i))
                           for i, c in enumerate(data['unique_comms'])}

    return data


def load_tcga_bulk_communities():
    """Load TCGA bulk communities for overlap analysis."""
    tcga_comms = {}
    cancer_type = "TCGA-HNSC"
    part_path = os.path.join(FIG2_COMMUNITIES, cancer_type,
                              f"{cancer_type}_best_partition.csv")
    if os.path.exists(part_path):
        df = pd.read_csv(part_path)
        ensg_to_sym = load_ensg_to_symbol()
        for _, row in df.iterrows():
            c = int(row['community'])
            gene = ensg_to_sym.get(row['gene'], row['gene'])
            tcga_comms.setdefault(c, set()).add(gene)
    return tcga_comms


def load_tcga_bulk_communities_fallback():
    """Fallback: try loading from community gene lists."""
    tcga_comms = {}
    cancer_type = "TCGA-HNSC"
    gene_list_path = os.path.join(FIG2_COMMUNITIES, cancer_type,
                                   f"{cancer_type}_community_gene_lists.csv")
    if os.path.exists(gene_list_path):
        df = pd.read_csv(gene_list_path)
        ensg_to_sym = load_ensg_to_symbol()
        for _, row in df.iterrows():
            c = int(row['community'])
            genes = row['genes'].split(';')
            tcga_comms[c] = set(ensg_to_sym.get(g, g) for g in genes)
    return tcga_comms


# =============================================================================
# GENE CLASS HELPERS
# =============================================================================

def get_gene_class(gene, harris_all, tcga_shared):
    """Return gene class with priority: A3 > TCGA > Harris > none."""
    if gene in A3_SET:
        return 'a3'
    if gene in tcga_shared:
        return 'tcga'
    if gene in harris_all:
        return 'harris'
    return 'none'


def get_ring_color(gene_class):
    """Ring color for gene class."""
    return {'a3': COLOR_A3, 'tcga': COLOR_TCGA,
            'harris': COLOR_HARRIS, 'none': None}[gene_class]


def get_text_color(gene_class):
    """Label text color for gene class."""
    return {'a3': COLOR_A3_TEXT, 'tcga': COLOR_TCGA,
            'harris': "#b8860b", 'none': "#000000"}[gene_class]


# =============================================================================
# PANEL 4a: UMAP
# =============================================================================

def panel_4a_umap(data):
    """UMAP with HIGH and A3+CNV-matched LOW highlighted."""
    banner("[PANEL 4a] UMAP cell selection")

    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    umap = adata.obsm['X_umap']

    fig, ax = plt.subplots(figsize=(14, 12))

    # All cells gray
    ax.scatter(umap[:, 0], umap[:, 1], c=COLOR_GRAY, s=2, alpha=0.15,
               edgecolors='none', rasterized=True)

    # LOW (A3+CNV matched)
    low_mask = adata.obs_names.isin(data['low_cells'])
    ax.scatter(umap[low_mask, 0], umap[low_mask, 1],
               c=COLOR_LOW, s=25, alpha=0.8,
               edgecolors='#000000', linewidths=0.3, rasterized=True)

    # HIGH (on top)
    high_mask = adata.obs_names.isin(data['high_cells'])
    ax.scatter(umap[high_mask, 0], umap[high_mask, 1],
               c=COLOR_HIGH, s=25, alpha=0.8,
               edgecolors='#000000', linewidths=0.3, rasterized=True)

    legend = [
        Patch(fc=COLOR_HIGH, ec='#000000', label=f"SBS2-HIGH (n={high_mask.sum():,})"),
        Patch(fc=COLOR_LOW, ec='#000000', label=f"A3+CNV LOW (n={low_mask.sum():,})"),
        Patch(fc=COLOR_GRAY, ec='#999999', label="Other cells"),
    ]
    ax.legend(handles=legend, fontsize=FONT_LEGEND - 4, framealpha=0.9,
              loc='upper right')
    ax.set_xlabel('UMAP 1', fontsize=FONT_AXIS)
    ax.set_ylabel('UMAP 2', fontsize=FONT_AXIS)
    ax.set_title('Basal Cell Selection', fontsize=FONT_TITLE, pad=15)
    ax.set_frame_on(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    plt.tight_layout()

    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(FIG_DIR, f"SC_Panel_4a_UMAP.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    del adata
    log(f"  [SAVE] Panel 4a UMAP")


# =============================================================================
# PANEL 4b: DUAL HEATMAP
# =============================================================================

def panel_4b_heatmap(data):
    """Dual heatmap: HIGH + DIFF correlation, community-sorted."""
    banner("[PANEL 4b] Dual heatmap (HIGH + DIFF)")

    corr_high = data.get('corr_high')
    corr_diff = data.get('corr_diff')
    if corr_high is None or corr_diff is None:
        log("  [SKIP] Missing correlation matrices")
        return

    partition_df = data['partition_df']
    gene_to_comm = data['gene_to_comm']

    ordered_genes, comm_boundaries, _ = build_community_ordered_genes(
        partition_df, corr_high, gene_to_comm)

    fig, axes = plt.subplots(1, 2, figsize=(26, 12))

    for ax, (label, matrix) in zip(axes, [
        ("HIGH (SBS2-high)", corr_high),
        ("DIFF (HIGH - LOW)", corr_diff),
    ]):
        genes_present = [g for g in ordered_genes if g in matrix.index]
        hm = matrix.loc[genes_present, genes_present].values

        vmin = np.nanpercentile(hm, 1)
        vmax = np.nanpercentile(hm, 99)
        if vmax <= vmin:
            vmax = 1.0

        im = ax.imshow(hm, aspect="auto", interpolation="nearest",
                        cmap="plasma", vmin=vmin, vmax=vmax)

        for start, end, c in comm_boundaries:
            if start > 0:
                ax.axhline(start - 0.5, color="white", lw=2.0, alpha=0.9)
                ax.axvline(start - 0.5, color="white", lw=2.0, alpha=0.9)

        for start, end, c in comm_boundaries:
            mid = (start + end) / 2
            if end - start >= 8:
                ax.annotate(f"C{c}", xy=(len(genes_present) + 3, mid),
                            fontsize=FONT_LEGEND - 4, fontweight="bold",
                            va="center", ha="left", annotation_clip=False)

        cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.08)
        cbar.set_label("Spearman rho", fontsize=FONT_LEGEND - 4, labelpad=6)
        cbar.ax.tick_params(labelsize=FONT_TICK - 4)

        ax.set_title(label, fontsize=FONT_TITLE, pad=12)
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(FIG_DIR, f"SC_Panel_4b_dual_heatmap.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Panel 4b dual heatmap")


# =============================================================================
# PANEL 4c: OVERLAP MATRIX
# =============================================================================

def panel_4c_overlap(data):
    """SC vs TCGA community overlap heatmap with hypergeometric test."""
    banner("[PANEL 4c] SC vs TCGA overlap")

    sc_comms = {}
    for _, row in data['partition_df'].iterrows():
        c = int(row['community'])
        sc_comms.setdefault(c, set()).add(row['gene'])

    tcga_comms = data.get('tcga_communities', {})
    if not tcga_comms:
        log("  [SKIP] No TCGA communities loaded")
        return

    all_sc = set()
    for v in sc_comms.values():
        all_sc.update(v)
    all_tcga = set()
    for v in tcga_comms.values():
        all_tcga.update(v)
    universe = len(all_sc | all_tcga)

    sc_ids = sorted(sc_comms.keys())
    tc_ids = sorted(tcga_comms.keys())

    overlap = np.zeros((len(sc_ids), len(tc_ids)))
    pvals   = np.zeros((len(sc_ids), len(tc_ids)))

    for i, si in enumerate(sc_ids):
        for j, tj in enumerate(tc_ids):
            shared = sc_comms[si] & tcga_comms[tj]
            overlap[i, j] = len(shared)
            k = len(shared)
            M = universe
            n = len(sc_comms[si])
            N = len(tcga_comms[tj])
            pvals[i, j] = hypergeom.sf(k - 1, M, N, n) if k > 0 else 1.0

    fdr = bh_fdr(pvals.flatten()).reshape(pvals.shape)

    # Save overlap data
    ensure_dir(DIR_06_OVERLAP)
    pd.DataFrame(overlap, index=[f"SC_C{c}" for c in sc_ids],
                 columns=[f"TCGA_C{c}" for c in tc_ids]).to_csv(
        os.path.join(DIR_06_OVERLAP, "overlap_matrix.csv"))
    pd.DataFrame(fdr, index=[f"SC_C{c}" for c in sc_ids],
                 columns=[f"TCGA_C{c}" for c in tc_ids]).to_csv(
        os.path.join(DIR_06_OVERLAP, "hypergeom_pvalues_BH.csv"))

    # Plot
    fig, ax = plt.subplots(figsize=(14, 12))
    im = ax.imshow(overlap, cmap="YlOrRd", aspect="auto")

    for i in range(len(sc_ids)):
        for j in range(len(tc_ids)):
            val = int(overlap[i, j])
            if val > 0:
                star = "*" if fdr[i, j] < 0.05 else ""
                ax.text(j, i, f"{val}{star}", ha='center', va='center',
                        fontsize=FONT_TICK - 6, fontweight='bold')

    ax.set_xticks(range(len(tc_ids)))
    ax.set_xticklabels([f"C{c}" for c in tc_ids], fontsize=FONT_TICK - 4, rotation=45)
    ax.set_yticks(range(len(sc_ids)))
    ax.set_yticklabels([f"C{c}" for c in sc_ids], fontsize=FONT_TICK - 4)
    ax.set_xlabel("TCGA Bulk Communities", fontsize=FONT_AXIS - 2)
    ax.set_ylabel("SC Communities", fontsize=FONT_AXIS - 2)
    ax.set_title("Cross-Resolution Community Overlap\n(* = FDR < 0.05)",
                 fontsize=FONT_TITLE, pad=15)

    cbar = plt.colorbar(im, ax=ax, shrink=0.7)
    cbar.set_label("Shared Genes", fontsize=FONT_LEGEND, labelpad=8)
    cbar.ax.tick_params(labelsize=FONT_TICK - 4)

    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(FIG_DIR, f"SC_Panel_4c_overlap.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Panel 4c overlap matrix")


# =============================================================================
# PANEL 4d: FULL NETWORK (intra_score sizing, three-class highlighting)
# =============================================================================

def panel_4d_network(data):
    """Exploded community network with intra_score sizing and gene class rings."""
    banner("[PANEL 4d] Full network (intra_score sizing)")

    G = data['G_comm']
    g2c = data['gene_to_comm']
    comm_to_nodes = data['comm_to_nodes']
    unique_comms = data['unique_comms']
    ccmap = data['comm_colors']
    node_intra = data['node_intra']
    harris_all = data['harris_all']
    tcga_shared = data['tcga_shared']

    def get_intra(n):
        return node_intra.get(n, 0.1)

    # ---- Exploded layout ----
    rng = np.random.default_rng(COMMUNITY_BASE_SEED)

    CG = nx.Graph()
    for c in unique_comms:
        CG.add_node(c)
    for u, v, d in G.edges(data=True):
        cu, cv = g2c.get(u, -1), g2c.get(v, -1)
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
            iterations=200)
    else:
        pos_comm = nx.circular_layout(CG)

    pos_comm = {c: INTER_SCALE * np.asarray(p, dtype=float)
                for c, p in pos_comm.items()}

    pos = {}
    for c in unique_comms:
        c_nodes = comm_to_nodes[c]
        sub = G.subgraph(c_nodes).copy()
        if len(c_nodes) == 1:
            pos[c_nodes[0]] = pos_comm.get(c, np.zeros(2))
            continue
        if sub.number_of_edges() > 0:
            k_val = 2.0 / np.sqrt(max(sub.number_of_nodes(), 1))
            pos_sub = nx.spring_layout(sub, seed=COMMUNITY_BASE_SEED,
                                       weight="abs_weight",
                                       k=k_val, iterations=200)
        else:
            pos_sub = {n: rng.normal(0, 0.05, size=2) for n in c_nodes}

        coords = np.array([pos_sub[n] for n in c_nodes], dtype=float)
        coords -= coords.mean(axis=0, keepdims=True)
        coords *= INTRA_SCALE
        center = pos_comm.get(c, np.zeros(2, dtype=float))
        for i, n in enumerate(c_nodes):
            pos[n] = center + coords[i]

    log(f"  Layout: {len(pos)} nodes, {len(unique_comms)} communities")

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(26, 24))

    # Edges
    for u, v, d in G.edges(data=True):
        cu, cv = g2c.get(u, -1), g2c.get(v, -1)
        w = d.get("weight", 0)
        if cu == cv:
            color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
            width = 0.8 + 4.0 * abs(w)
            alpha = min(0.55, 0.15 + 0.35 * abs(w))
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color=color, alpha=alpha, linewidth=width, zorder=1)
        else:
            width = 0.4 + 1.5 * abs(w)
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color="#999999", alpha=0.10, linewidth=width,
                    linestyle=(0, (3, 3)), zorder=0)

    # Nodes (community-colored, intra_score sized)
    a3_in_graph = [g for g in A3_GENES_SYMBOLS if g in G.nodes()]

    for c in unique_comms:
        c_nodes = comm_to_nodes.get(c, [])
        if not c_nodes:
            continue
        sizes = [100 + 1160 * get_intra(n) for n in c_nodes]
        nx.draw_networkx_nodes(
            G, pos, nodelist=c_nodes, ax=ax, node_size=sizes,
            node_color=[ccmap[c]] * len(c_nodes),
            alpha=0.85, linewidths=0.5, edgecolors="#000000",
            label=f"C{c} (n={len(c_nodes)})")

    # Three-class ring highlighting (A3 gets filled + ring, others get ring only)
    # A3 genes: coral fill, enlarged, red ring
    if a3_in_graph:
        a3_sizes = [100 + 1160 * get_intra(n) for n in a3_in_graph]
        a3_sizes_big = [s * 1.8 for s in a3_sizes]
        nx.draw_networkx_nodes(
            G, pos, nodelist=a3_in_graph, ax=ax,
            node_size=a3_sizes_big, node_color=COLOR_A3, alpha=1.0,
            edgecolors="#000000", linewidths=2.5)
        ring_sizes = [s * 1.8 for s in a3_sizes_big]
        nx.draw_networkx_nodes(
            G, pos, nodelist=a3_in_graph, ax=ax,
            node_size=ring_sizes, node_color="none",
            edgecolors=COLOR_A3, linewidths=4.0, alpha=1.0)

    # TCGA shared: orange ring (exclude A3 genes)
    tcga_in_graph = [n for n in G.nodes()
                     if n in tcga_shared and n not in A3_SET]
    if tcga_in_graph:
        ring_sizes = [(100 + 1160 * get_intra(n)) * 2.5 for n in tcga_in_graph]
        nx.draw_networkx_nodes(
            G, pos, nodelist=tcga_in_graph, ax=ax,
            node_size=ring_sizes, node_color="none",
            edgecolors=COLOR_TCGA, linewidths=4.0, alpha=0.9)

    # Harris interactors: gold ring (exclude A3 and TCGA)
    harris_in_graph = [n for n in G.nodes()
                       if n in harris_all and n not in A3_SET
                       and n not in tcga_shared]
    if harris_in_graph:
        ring_sizes = [(100 + 1160 * get_intra(n)) * 2.5 for n in harris_in_graph]
        nx.draw_networkx_nodes(
            G, pos, nodelist=harris_in_graph, ax=ax,
            node_size=ring_sizes, node_color="none",
            edgecolors=COLOR_HARRIS, linewidths=3.5, alpha=0.9)

    # Labels (top 3 hubs per community + all A3/TCGA/Harris genes)
    labels = {}
    font_sizes = {}
    font_colors = {}

    for n in a3_in_graph:
        labels[n] = n
        font_sizes[n] = FONT_AXIS
        font_colors[n] = COLOR_A3_TEXT

    for n in tcga_in_graph:
        if n not in labels:
            labels[n] = n
            font_sizes[n] = 15 + 28 * get_intra(n)
            font_colors[n] = COLOR_TCGA

    for n in harris_in_graph:
        if n not in labels:
            labels[n] = n
            font_sizes[n] = 15 + 28 * get_intra(n)
            font_colors[n] = "#b8860b"

    for c in unique_comms:
        c_nodes = comm_to_nodes.get(c, [])
        if len(c_nodes) < 5:
            continue
        scored = [(n, get_intra(n)) for n in c_nodes]
        scored.sort(key=lambda x: -x[1])
        for n, sc in scored[:3]:
            if n not in labels:
                labels[n] = n
                font_sizes[n] = 15 + 28 * get_intra(n)
                font_colors[n] = "#000000"

    log(f"  Labels: {len(labels)} nodes")

    label_pos = repel_labels(pos, labels, data_range=14.0,
                             n_iters=300, push_strength=0.08)
    draw_labels_with_boxes(ax, pos, label_pos, labels,
                           font_sizes, font_colors, char_w=0.07)

    # Community center annotations
    for c in unique_comms:
        c_nodes = comm_to_nodes.get(c, [])
        if len(c_nodes) < 3:
            continue
        cx = np.mean([pos[n][0] for n in c_nodes])
        cy = np.max([pos[n][1] for n in c_nodes]) + 0.25
        ax.annotate(f"C{c}", xy=(cx, cy), fontsize=18, fontweight="bold",
                    color=ccmap[c], ha="center", va="bottom", alpha=0.8)

    # Legend
    legend_handles = [
        Patch(fc=COLOR_A3, ec="#000000", lw=2, label="APOBEC3 gene"),
        Patch(fc="white", ec=COLOR_TCGA, lw=3, label="TCGA shared gene"),
        Patch(fc="white", ec=COLOR_HARRIS, lw=3, label="Harris A3 interactor"),
        Patch(fc=COLOR_EDGE_POS, ec="none", alpha=0.5, label="Gained co-expression"),
        Patch(fc=COLOR_EDGE_NEG, ec="none", alpha=0.5, label="Lost co-expression"),
    ]
    ax.legend(handles=legend_handles, loc="upper left",
              fontsize=FONT_LEGEND, framealpha=0.9,
              title="Legend", title_fontsize=FONT_LEGEND + 2)

    ax.set_title(
        f"SC DIFF Network -- {G.number_of_nodes()} genes, "
        f"{G.number_of_edges()} edges, {len(unique_comms)} communities",
        fontsize=FONT_TITLE, pad=15)
    ax.axis("off")
    plt.tight_layout()

    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(FIG_DIR, f"SC_Panel_4d_network_full.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Panel 4d full network")

    return pos  # Return layout for supplement panel


# =============================================================================
# PANEL 4e: HARRIS ENRICHMENT
# =============================================================================

def panel_4e_harris(data):
    """Harris A3 interactor enrichment per community."""
    banner("[PANEL 4e] Harris interactor enrichment")

    sc_comms = {}
    for _, row in data['partition_df'].iterrows():
        c = int(row['community'])
        sc_comms.setdefault(c, set()).add(row['gene'])

    harris_all = data['harris_all']
    if not harris_all:
        log("  [SKIP] No Harris interactors")
        return

    all_genes = set()
    for v in sc_comms.values():
        all_genes.update(v)
    universe = len(all_genes)

    results = []
    for cid in sorted(sc_comms.keys()):
        genes = sc_comms[cid]
        overlap = genes & harris_all
        k = len(overlap)
        n = len(genes)
        N = len(harris_all & all_genes)
        p = hypergeom.sf(k - 1, universe, N, n) if k > 0 else 1.0
        results.append({
            'community': cid, 'size': n, 'harris_count': k,
            'harris_genes': ';'.join(sorted(overlap)),
            'p_value': p
        })

    rdf = pd.DataFrame(results)
    rdf['fdr'] = bh_fdr(rdf['p_value'].values)
    rdf.to_csv(os.path.join(DIR_06_OVERLAP, "harris_enrichment_all.csv"), index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(14, 8))
    cids = rdf['community'].values
    counts = rdf['harris_count'].values
    colors = [data['comm_colors'].get(c, '#888888') for c in cids]

    bars = ax.bar(range(len(cids)), counts, color=colors, edgecolor='#000000', linewidth=0.8)

    for i, (c, cnt, fdr) in enumerate(zip(cids, counts, rdf['fdr'])):
        if cnt > 0:
            star = " *" if fdr < 0.05 else ""
            ax.text(i, cnt + 0.2, f"{cnt}{star}", ha='center', va='bottom',
                    fontsize=FONT_TICK - 4, fontweight='bold')

    ax.set_xticks(range(len(cids)))
    ax.set_xticklabels([f"C{c}" for c in cids], fontsize=FONT_TICK - 2)
    ax.set_ylabel("Harris A3 Interactors", fontsize=FONT_AXIS - 2)
    ax.set_xlabel("SC Community", fontsize=FONT_AXIS - 2)
    ax.set_title("Known A3 Protein Interactor Enrichment\n(* = FDR < 0.05)",
                 fontsize=FONT_TITLE, pad=15)
    ax.tick_params(axis='y', labelsize=FONT_TICK - 2)

    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(FIG_DIR, f"SC_Panel_4e_harris_enrichment.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Panel 4e Harris enrichment")


# =============================================================================
# COMMUNITY ZOOMS (intra_score sizing + three-class rings)
# =============================================================================

def plot_community_zooms(data):
    """Per-community zoomed plots with intra_score sizing and gene class rings."""
    banner("[SUPPLEMENT] Per-community zooms")

    G = data['G_comm']
    g2c = data['gene_to_comm']
    comm_to_nodes = data['comm_to_nodes']
    unique_comms = data['unique_comms']
    ccmap = data['comm_colors']
    node_intra = data['node_intra']
    harris_all = data['harris_all']
    tcga_shared = data['tcga_shared']

    zoom_dir = ensure_dir(os.path.join(FIG_DIR, "community_zooms"))

    def get_intra(n):
        return node_intra.get(n, 0.1)

    for c in unique_comms:
        c_nodes = comm_to_nodes.get(c, [])
        if len(c_nodes) < 5:
            continue

        G_sub = G.subgraph(c_nodes).copy()
        if G_sub.number_of_edges() == 0:
            continue

        pos_sub = nx.spring_layout(
            G_sub, seed=COMMUNITY_BASE_SEED, weight="abs_weight",
            k=3.0 / np.sqrt(max(G_sub.number_of_nodes(), 1)),
            iterations=300)

        fig, ax = plt.subplots(figsize=(16, 14))

        # Edges (doubled thickness)
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
            node_color=[ccmap[c]] * G_sub.number_of_nodes(),
            alpha=0.85, edgecolors="#000000", linewidths=0.8)

        # Three-class rings
        a3_here = [n for n in G_sub.nodes() if n in A3_SET]
        tcga_here = [n for n in G_sub.nodes()
                     if n in tcga_shared and n not in A3_SET]
        harris_here = [n for n in G_sub.nodes()
                       if n in harris_all and n not in A3_SET
                       and n not in tcga_shared]

        if a3_here:
            a3_s = [(500 + 3500 * get_intra(n)) * 1.5 for n in a3_here]
            nx.draw_networkx_nodes(
                G_sub, pos_sub, nodelist=a3_here, ax=ax,
                node_size=a3_s, node_color=COLOR_A3,
                edgecolors="#000000", linewidths=2.0)
            ring_s = [s * 2.0 for s in a3_s]
            nx.draw_networkx_nodes(
                G_sub, pos_sub, nodelist=a3_here, ax=ax,
                node_size=ring_s, node_color="none",
                edgecolors=COLOR_A3, linewidths=5.0, alpha=1.0)

        if tcga_here:
            ring_s = [(500 + 3500 * get_intra(n)) * 2.5 for n in tcga_here]
            nx.draw_networkx_nodes(
                G_sub, pos_sub, nodelist=tcga_here, ax=ax,
                node_size=ring_s, node_color="none",
                edgecolors=COLOR_TCGA, linewidths=4.0, alpha=0.9)

        if harris_here:
            ring_s = [(500 + 3500 * get_intra(n)) * 2.5 for n in harris_here]
            nx.draw_networkx_nodes(
                G_sub, pos_sub, nodelist=harris_here, ax=ax,
                node_size=ring_s, node_color="none",
                edgecolors=COLOR_HARRIS, linewidths=3.5, alpha=0.9)

        # Labels
        zoom_labels = {}
        zoom_fs = {}
        zoom_fc = {}
        special = set(a3_here + tcga_here + harris_here)

        if len(c_nodes) <= 40:
            for n in G_sub.nodes():
                zoom_labels[n] = n
                gc = get_gene_class(n, harris_all, tcga_shared)
                zoom_fs[n] = 15 + 21 * get_intra(n)
                zoom_fc[n] = get_text_color(gc)
                if n in A3_SET:
                    zoom_fs[n] = FONT_AXIS
        else:
            scored = [(n, get_intra(n)) for n in G_sub.nodes()]
            scored.sort(key=lambda x: -x[1])
            for n, sc in scored[:15]:
                zoom_labels[n] = n
                gc = get_gene_class(n, harris_all, tcga_shared)
                zoom_fs[n] = 15 + 21 * get_intra(n)
                zoom_fc[n] = get_text_color(gc)
            for n in special:
                zoom_labels[n] = n
                gc = get_gene_class(n, harris_all, tcga_shared)
                zoom_fs[n] = max(zoom_fs.get(n, 0), 15 + 21 * get_intra(n))
                zoom_fc[n] = get_text_color(gc)
                if n in A3_SET:
                    zoom_fs[n] = FONT_AXIS

        if zoom_labels:
            zl_pos = repel_labels(pos_sub, zoom_labels, data_range=3.0,
                                  n_iters=250, push_strength=0.07)
            draw_labels_with_boxes(ax, pos_sub, zl_pos, zoom_labels,
                                   zoom_fs, zoom_fc, char_w=0.04)

        # Title
        a3_sym = [A3_SYMBOL_TO_ALIAS.get(n, n) for n in a3_here]
        tcga_count = len(tcga_here)
        harris_count = len(harris_here)
        notes = []
        if a3_sym:
            notes.append(', '.join(a3_sym))
        if tcga_count:
            notes.append(f"{tcga_count} TCGA shared")
        if harris_count:
            notes.append(f"{harris_count} Harris")
        note_str = f" -- {' | '.join(notes)}" if notes else ""

        ax.set_title(
            f"Community {c} -- {len(c_nodes)} genes, "
            f"{G_sub.number_of_edges()} edges{note_str}",
            fontsize=FONT_AXIS, pad=10)
        ax.axis("off")
        plt.tight_layout()

        for ext in ['png', 'pdf']:
            plt.savefig(os.path.join(zoom_dir, f"SC_community_{c:02d}.{ext}"),
                        dpi=300, bbox_inches='tight')
        plt.close()
        log(f"  Community {c} ({len(c_nodes)} genes)")

    log(f"  [SAVE] All zooms -> {zoom_dir}")


# =============================================================================
# SUPPLEMENT: METHODS OVERVIEW
# =============================================================================

def supplement_methods(data, network_pos):
    """Methods graphical abstract: HIGH -> LOW -> DIFF -> mini network."""
    banner("[SUPPLEMENT] Methods graphical abstract")

    corr_high = data.get('corr_high')
    corr_low  = data.get('corr_low')
    corr_diff = data.get('corr_diff')
    G = data['G_comm']

    if corr_high is None or corr_low is None or corr_diff is None:
        log("  [SKIP] Missing correlation matrices")
        return
    if network_pos is None:
        log("  [SKIP] No network positions available")
        return

    partition_df = data['partition_df']
    gene_to_comm = data['gene_to_comm']
    comm_to_nodes = data['comm_to_nodes']
    unique_comms = data['unique_comms']
    ccmap = data['comm_colors']
    node_intra = data['node_intra']

    ordered_genes, comm_boundaries, _ = build_community_ordered_genes(
        partition_df, corr_high, gene_to_comm)

    fig = plt.figure(figsize=(52, 13))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.12)

    for idx, (label, matrix, cmap_name) in enumerate([
        ("HIGH (SBS2-high)", corr_high, "plasma"),
        ("LOW (A3+CNV control)", corr_low, "plasma"),
        ("DIFF (HIGH - LOW)", corr_diff, "plasma"),
    ]):
        ax = fig.add_subplot(gs[0, idx])
        genes_present = [g for g in ordered_genes if g in matrix.index]
        hm = matrix.loc[genes_present, genes_present].values

        vmin = np.nanpercentile(hm, 1)
        vmax = np.nanpercentile(hm, 99)
        if vmax <= vmin:
            vmax = 1.0

        im = ax.imshow(hm, aspect="auto", interpolation="nearest",
                        cmap=cmap_name, vmin=vmin, vmax=vmax)

        for start, end, c in comm_boundaries:
            if start > 0:
                ax.axhline(start - 0.5, color="white", lw=2.5, alpha=0.9)
                ax.axvline(start - 0.5, color="white", lw=2.5, alpha=0.9)

        for start, end, c in comm_boundaries:
            mid = (start + end) / 2
            if end - start >= 8:
                ax.annotate(f"C{c}", xy=(len(genes_present) + 3, mid),
                            fontsize=FONT_LEGEND, fontweight="bold",
                            va="center", ha="left", annotation_clip=False)

        cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.10)
        cbar.set_label("Spearman rho", fontsize=FONT_LEGEND, labelpad=8)
        cbar.ax.tick_params(labelsize=FONT_TICK)

        ax.set_title(label, fontsize=FONT_TITLE, pad=12)
        ax.set_xticks([])
        ax.set_yticks([])

        if idx == 0:
            ax.annotate("  -  ", xy=(1.06, 0.5), xycoords="axes fraction",
                        fontsize=50, fontweight="bold", color="#333333",
                        ha="center", va="center")
        elif idx == 1:
            ax.annotate("  =  ", xy=(1.06, 0.5), xycoords="axes fraction",
                        fontsize=50, fontweight="bold", color="#333333",
                        ha="center", va="center")
        elif idx == 2:
            ax.annotate("  >>  ", xy=(1.06, 0.5), xycoords="axes fraction",
                        fontsize=40, fontweight="bold", color="#333333",
                        ha="center", va="center")

    # Mini network
    ax_net = fig.add_subplot(gs[0, 3])
    pos = network_pos

    def get_intra(n):
        return node_intra.get(n, 0.1)

    for u, v, d in G.edges(data=True):
        cu, cv = gene_to_comm.get(u, -1), gene_to_comm.get(v, -1)
        if cu != cv:
            continue
        w = d.get("weight", 0)
        color = COLOR_EDGE_POS if w > 0 else COLOR_EDGE_NEG
        ax_net.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color=color, alpha=0.15, linewidth=0.5)

    for c in unique_comms:
        c_nodes = comm_to_nodes.get(c, [])
        if not c_nodes:
            continue
        xs = [pos[n][0] for n in c_nodes]
        ys = [pos[n][1] for n in c_nodes]
        sizes = [15 + 80 * get_intra(n) for n in c_nodes]
        ax_net.scatter(xs, ys, s=sizes, c=ccmap[c], alpha=0.8,
                       edgecolors="#000000", linewidths=0.2)

    for n in A3_GENES_SYMBOLS:
        if n in pos:
            ax_net.scatter(pos[n][0], pos[n][1], s=120,
                           c=COLOR_A3, edgecolors="#000000",
                           linewidths=1.5, zorder=5)

    for c in unique_comms:
        c_nodes = comm_to_nodes.get(c, [])
        if len(c_nodes) < 5:
            continue
        cx = np.mean([pos[n][0] for n in c_nodes])
        cy = np.max([pos[n][1] for n in c_nodes]) + 0.15
        ax_net.annotate(f"C{c}", xy=(cx, cy), fontsize=14,
                        fontweight="bold", color=ccmap[c],
                        ha="center", va="bottom", alpha=0.9)

    ax_net.set_title("Network", fontsize=FONT_TITLE, pad=12)
    ax_net.axis("off")

    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(FIG_DIR, f"SC_Supplement_methods_overview.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  [SAVE] Supplement methods overview")


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("GENERATE FIGURE 4 PANELS")
    log(f"  Start: {t0}")
    log(f"  Output: {FIG_DIR}")

    data = load_all()

    # Panel 4a: UMAP
    panel_4a_umap(data)

    # Panel 4b: Dual heatmap
    panel_4b_heatmap(data)

    # Panel 4c: Overlap matrix
    panel_4c_overlap(data)

    # Panel 4d: Full network
    network_pos = panel_4d_network(data)

    # Panel 4e: Harris enrichment
    panel_4e_harris(data)

    # Community zooms
    plot_community_zooms(data)

    # Supplement: Methods overview
    supplement_methods(data, network_pos)

    banner("FIGURE 4 PANELS COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")
    log(f"  Output: {FIG_DIR}")


if __name__ == "__main__":
    main()
