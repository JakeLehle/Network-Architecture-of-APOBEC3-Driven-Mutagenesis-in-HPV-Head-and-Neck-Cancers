"""
network_utils.py

Utility library for APOBEC gene-network analysis.

Main capabilities:
1) Build gene–gene correlation networks (weighted graphs)
2) Compute APOBEC A3 scores and QC plots
3) Compute network metrics
4) Plot networks with fixed layouts and optional annotations
5) Compare TOP vs BOTTOM networks and analyze DIFF networks
6) Community detection + community-aware layouts
7) Correlation diagnostics and clustered heatmaps
8) Export helpers (CSVs, edge lists, partitions)

Designed to be imported into the main analysis pipeline.

============================================================
TABLE OF CONTENTS
============================================================
1. Imports
2. Logging
3. Graph construction & export
4. APOBEC A3 score & QC plots
5. Network metrics
6. Layout control & node filtering
7. Network plotting (fixed positions)
8. Multiple testing & correlation diagnostics
9. Clustered heatmaps & dendrograms
10. Community detection (Leiden/Louvain) + layouts
11. Community heatmap utilities
12. SBS2 vs A3 selection plot
13. Community ordering helpers + simple heatmaps
14. Modularity + community detection v2
15. Community-aware layout v3 + zoomed community plots
16. Resolution sweep stability
17. Violin correlations
18. Save partition outputs
19. DIFF graph strategies (top-k, soft threshold WGCNA-like)
20. Volcano & Manhattan plots
21. Network saving helpers (STEP 8)
22. Misc utilities (ENS cleanup, banners, id shorten, etc.)
============================================================
"""

# ============================================================
# 1) Imports (tidied: single import block, no duplicates)
# ============================================================
import os
import re
import json
import csv

import numpy as np
import pandas as pd
import networkx as nx

import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.offsetbox import AnchoredText
from matplotlib import colormaps

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


# ============================================================
# 2) Logging
# ============================================================

# Helper: lightweight console logger used throughout the pipeline for optional verbosity.
def log(msg: str, verbose: bool = True):
    """
    Simple console logger.

    Parameters
    ----------
    msg : str
        Message to print.
    verbose : bool
        If False, suppress output.
    """
    if verbose:
        print(msg, flush=True)


# ============================================================
# 3) Graph construction & export
# ============================================================

# Helper: build an undirected weighted gene graph from a correlation matrix with a hard |corr| threshold.
def build_weighted_graph_from_corr(corr_df: pd.DataFrame, threshold: float) -> nx.Graph:
    """
    Build an undirected weighted graph from a correlation matrix.

    - Nodes = genes
    - Edges kept if |correlation| >= threshold
    - Keeps isolated nodes
    - Edge attributes:
        * weight     : signed correlation
        * abs_weight : absolute correlation
    """
    corr_df = corr_df.copy().fillna(0)
    corr_df = corr_df.apply(pd.to_numeric, errors="coerce").fillna(0)

    genes = list(corr_df.index)
    G = nx.Graph()
    G.add_nodes_from(genes)

    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            g1, g2 = genes[i], genes[j]
            w = float(corr_df.iat[i, j])
            if abs(w) >= threshold and abs(w) < 0.999999:
                G.add_edge(g1, g2, weight=w, abs_weight=abs(w))
    return G


# Helper: export edges and basic node centralities to CSV for downstream tools (Cytoscape/Gephi, QC).
def save_graph_csv(G: nx.Graph, out_dir: str, out_prefix: str):
    """
    Save network structure and node metrics to CSV files.

    Outputs:
    - <prefix>_edges.csv : source, target, weight, abs_weight
    - <prefix>_nodes.csv : node, degree, betweenness
    """
    os.makedirs(out_dir, exist_ok=True)

    edges = []
    for u, v, data in G.edges(data=True):
        edges.append({
            "source": u,
            "target": v,
            "weight": float(data.get("weight", 1.0)),
            "abs_weight": float(data.get("abs_weight", abs(data.get("weight", 1.0))))
        })
    pd.DataFrame(edges).to_csv(os.path.join(out_dir, f"{out_prefix}_edges.csv"), index=False)

    deg = dict(G.degree())
    btw = nx.betweenness_centrality(G) if G.number_of_edges() > 0 else {n: 0.0 for n in G.nodes()}

    nodes_df = pd.DataFrame({
        "node": list(G.nodes()),
        "degree": [deg[n] for n in G.nodes()],
        "betweenness": [btw[n] for n in G.nodes()],
    })
    nodes_df.to_csv(os.path.join(out_dir, f"{out_prefix}_nodes.csv"), index=False)


# ============================================================
# 4) APOBEC A3 score & QC plots
# ============================================================

# Helper: compute the APOBEC A3 score using only A3A and A3B, normalized within the cohort.
def compute_a3_score_AB(df_in: pd.DataFrame, a_col="A3A", b_col="A3B") -> pd.Series:
    """
    Compute APOBEC A3 score using ONLY A3A and A3B.

    - Each gene is normalized by its max value within the cancer type
    - A3_score = normalized A3A + normalized A3B
    """
    tmp = df_in[[a_col, b_col]].copy()
    tmp[a_col] = pd.to_numeric(tmp[a_col], errors="coerce")
    tmp[b_col] = pd.to_numeric(tmp[b_col], errors="coerce")

    a_max = tmp[a_col].max()
    b_max = tmp[b_col].max()

    a_norm = 0.0 if (pd.isna(a_max) or a_max == 0) else (tmp[a_col] / a_max)
    b_norm = 0.0 if (pd.isna(b_max) or b_max == 0) else (tmp[b_col] / b_max)

    return (a_norm + b_norm).astype(float)


# Helper: produce quick QC figures for A3 score relationships and distributions.
def plot_a3_scatter_and_distribution(cancer_df: pd.DataFrame, out_dir: str, cancer_type: str):
    """
    Generate QC plots for APOBEC A3 activity.

    Plots:
    1) Scatter: A3_score vs SBS2 (each dot = donor)
    2) Histogram: distribution of A3_score across donors
    """
    os.makedirs(out_dir, exist_ok=True)

    plt.figure(figsize=(7, 6))
    plt.scatter(cancer_df["A3_score"], cancer_df["SBS2"], alpha=0.6)
    plt.xlabel("A3_score (A3A + A3B)")
    plt.ylabel("SBS2")
    plt.title(f"{cancer_type} | donor-level scatter")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{cancer_type}_scatter_A3score_vs_SBS2.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.hist(cancer_df["A3_score"].dropna().values, bins=40)
    plt.xlabel("A3_score")
    plt.ylabel("Donor count")
    plt.title(f"{cancer_type} | A3_score distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{cancer_type}_dist_A3score.png"), dpi=300)
    plt.close()


# ============================================================
# 5) Network metrics
# ============================================================

# Helper: compute per-node graph metrics for downstream analysis/plots.
def compute_graph_metrics_df(G: nx.Graph) -> pd.DataFrame:
    """
    Compute per-gene network metrics.

    Metrics:
    - degree
    - betweenness centrality
    - strength_abs (sum of |edge weights|)
    - closeness centrality
    """
    nodes = list(G.nodes())
    deg = dict(G.degree())
    betw = nx.betweenness_centrality(G) if G.number_of_edges() > 0 else {n: 0.0 for n in nodes}
    close = nx.closeness_centrality(G) if G.number_of_edges() > 0 else {n: 0.0 for n in nodes}

    strength = {}
    for n in nodes:
        s = 0.0
        for _, data in G[n].items():
            s += float(data.get("abs_weight", abs(data.get("weight", 1.0))))
        strength[n] = s

    return pd.DataFrame({
        "node": nodes,
        "degree": [deg[n] for n in nodes],
        "betweenness": [betw[n] for n in nodes],
        "strength_abs": [strength[n] for n in nodes],
        "closeness": [close[n] for n in nodes],
    })


# Helper: plot histograms of metric distributions to help pick thresholds / sanity check networks.
def plot_metric_distributions(metrics_df: pd.DataFrame, out_dir: str, prefix: str):
    """
    Plot distributions of network metrics.
    Helps decide correlation thresholds by visual inspection.
    """
    os.makedirs(out_dir, exist_ok=True)

    for col in ["degree", "betweenness", "strength_abs", "closeness"]:
        plt.figure(figsize=(7, 5))
        vals = pd.to_numeric(metrics_df[col], errors="coerce").dropna().values
        plt.hist(vals, bins=40)
        plt.xlabel(col)
        plt.ylabel("Gene count")
        plt.title(f"{prefix} | {col} distribution")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{prefix}_dist_{col}.png"), dpi=300)
        plt.close()


# ============================================================
# 6) Layout control & node filtering
# ============================================================

# Helper: create a stable layout so multiple graphs can be compared visually using fixed node positions.
def get_constant_layout_from_bottom(G_bottom: nx.Graph, G_union: nx.Graph, seed: int = 42) -> dict:
    """
    Generate a fixed node layout shared across multiple plots.

    Strategy:
    - Prefer layout computed from bottom network
    - Ensure all union nodes receive coordinates
    """
    base = G_bottom if G_bottom.number_of_edges() > 0 else G_union
    pos = nx.spring_layout(base, seed=seed)

    for n in G_union.nodes():
        if n not in pos:
            pos[n] = np.array([0.0, 0.0])
    return pos


# Helper: remove genes that are isolated in BOTH top and bottom networks (to reduce clutter in union plots).
def remove_isolates_in_both(G_top: nx.Graph, G_bottom: nx.Graph, G_any: nx.Graph) -> nx.Graph:
    """
    Remove genes that are isolated in BOTH top and bottom networks.
    """
    keep = []
    for n in G_any.nodes():
        if (G_top.degree(n) > 0) or (G_bottom.degree(n) > 0):
            keep.append(n)
    return G_any.subgraph(keep).copy()


# ============================================================
# 7) Network plotting (fixed positions)
# ============================================================

# Helper: plot a single network using a fixed pos dict; optional legend textbox; optional APOBEC labels.
def plot_global_network_fixedpos(
    G: nx.Graph,
    pos: dict,
    out_path: str,
    a3_genes: list,
    title: str,
    weighted: bool = True,
    dpi: int = 300,
    label_apobec: bool = True,
    legend_text=None,
    legend_loc="upper left"
):
    node_colors = ["red" if n in a3_genes else "gray" for n in G.nodes()]

    if weighted and G.number_of_edges() > 0:
        abs_w = np.array([float(G[u][v].get("abs_weight", 1.0)) for u, v in G.edges()])
        if abs_w.max() > abs_w.min():
            widths = 0.5 + 5.0 * (abs_w - abs_w.min()) / (abs_w.max() - abs_w.min())
        else:
            widths = np.ones_like(abs_w) * 2.0
    else:
        widths = 0.8

    # use explicit fig/ax
    fig, ax = plt.subplots(figsize=(18, 14))

    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=40, node_color=node_colors, alpha=0.9)
    nx.draw_networkx_edges(G, pos, ax=ax, width=widths, alpha=0.30)

    if label_apobec:
        labels = {n: n for n in G.nodes() if n in a3_genes}
        nx.draw_networkx_labels(G, pos, ax=ax, labels=labels, font_size=10)

    ax.set_title(title)
    ax.axis("off")

    # legend textbox
    if legend_text is not None and str(legend_text).strip() != "":
        at = AnchoredText(legend_text, loc=legend_loc, prop={"size": 8}, frameon=True)
        at.patch.set_alpha(0.85)
        ax.add_artist(at)

    # prevents legend cropping
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# Helper: plot union of TOP and BOTTOM networks and color edges by uniqueness/shared status.
def plot_combined_network_fixedpos(G_top, G_bottom, pos, out_path, a3_genes, title, dpi=300):
    """
    Plot combined top + bottom networks using fixed positions.

    Edge coloring:
    - Red  : edges unique to top group
    - Blue : edges unique to bottom group
    - Gray : edges shared between both
    """
    G_union = nx.compose(G_top, G_bottom)

    top_edges = set(tuple(sorted(e)) for e in G_top.edges())
    bot_edges = set(tuple(sorted(e)) for e in G_bottom.edges())
    only_top = top_edges - bot_edges
    only_bot = bot_edges - top_edges
    both = top_edges & bot_edges

    node_colors = ["red" if n in a3_genes else "gray" for n in G_union.nodes()]

    plt.figure(figsize=(18, 14))
    nx.draw_networkx_nodes(G_union, pos, node_size=40, node_color=node_colors, alpha=0.9)
    nx.draw_networkx_edges(G_union, pos, edgelist=list(both), alpha=0.10)
    nx.draw_networkx_edges(G_union, pos, edgelist=list(only_top), edge_color="red", alpha=0.35)
    nx.draw_networkx_edges(G_union, pos, edgelist=list(only_bot), edge_color="blue", alpha=0.35)

    labels = {n: n for n in G_union.nodes() if n in a3_genes}
    nx.draw_networkx_labels(G_union, pos, labels=labels, font_size=10)

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


# ============================================================
# 8) Multiple testing & correlation diagnostics
# ============================================================

# Helper: Benjamini–Hochberg procedure used for q-values (FDR control).
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini–Hochberg FDR correction.
    Returns q-values (same length as pvals).
    """
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


# Helper: pull upper-triangle values from a correlation matrix for distribution diagnostics.
def corr_upper_values(corr_df: pd.DataFrame, drop_diag: bool = True) -> np.ndarray:
    """
    Extract upper-triangle correlation values from a square correlation matrix.
    Returns a 1D numpy array of values (optionally excludes diagonal).
    """
    M = corr_df.values.astype(float)
    n = M.shape[0]
    k = 1 if drop_diag else 0
    vals = M[np.triu_indices(n, k=k)]
    vals = vals[~np.isnan(vals)]
    return vals


# Helper: extract only correlations above a threshold (mirrors what becomes edges in graphs).
def corr_edge_values(corr_df: pd.DataFrame, threshold: float) -> np.ndarray:
    """
    Extract only correlation values whose absolute value >= threshold (excluding diagonal).
    Useful to compare 'all correlations' vs 'edges we keep' distributions.
    """
    M = corr_df.values.astype(float)
    n = M.shape[0]
    vals = M[np.triu_indices(n, k=1)]
    vals = vals[~np.isnan(vals)]
    vals = vals[np.abs(vals) >= threshold]
    return vals


# Helper: generate box+jitter and histogram summaries for TOP/BOTTOM/DIFF correlation distributions.
def plot_corr_distributions(
    corr_top: pd.DataFrame,
    corr_bottom: pd.DataFrame,
    out_dir: str,
    prefix: str,
    threshold: float = None,
    seed: int = 42,
    dpi: int = 300
):
    """
    Creates distribution plots for:
      - all correlations (upper triangle, excluding diagonal)
      - optional: edge-only correlations where |corr| >= threshold
    """
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)

    corr_diff = corr_top - corr_bottom  # signed difference

    # ---- ALL values ----
    vals_top = corr_upper_values(corr_top)
    vals_bot = corr_upper_values(corr_bottom)
    vals_diff = corr_upper_values(corr_diff)

    _box_jitter_plot(
        [vals_top, vals_bot, vals_diff],
        labels=["TOP", "BOTTOM", "DIFF (TOP-BOT)"],
        title=f"{prefix} | Correlation distributions (ALL pairs)",
        out_path=os.path.join(out_dir, f"{prefix}_corr_box_jitter_ALL.png"),
        rng=rng,
        dpi=dpi
    )

    _hist_plot(
        [vals_top, vals_bot, vals_diff],
        labels=["TOP", "BOTTOM", "DIFF (TOP-BOT)"],
        title=f"{prefix} | Correlation histograms (ALL pairs)",
        out_path=os.path.join(out_dir, f"{prefix}_corr_hist_ALL.png"),
        dpi=dpi
    )

    # ---- EDGE-only values ----
    if threshold is not None:
        e_top = corr_edge_values(corr_top, threshold)
        e_bot = corr_edge_values(corr_bottom, threshold)
        e_diff = corr_edge_values(corr_diff, threshold)

        _box_jitter_plot(
            [e_top, e_bot, e_diff],
            labels=[f"TOP | |corr|≥{threshold}", f"BOTTOM | |corr|≥{threshold}", f"DIFF | |diff|≥{threshold}"],
            title=f"{prefix} | Correlation distributions (EDGE-only)",
            out_path=os.path.join(out_dir, f"{prefix}_corr_box_jitter_EDGES.png"),
            rng=rng,
            dpi=dpi
        )

        _hist_plot(
            [e_top, e_bot, e_diff],
            labels=[f"TOP | |corr|≥{threshold}", f"BOTTOM | |corr|≥{threshold}", f"DIFF | |diff|≥{threshold}"],
            title=f"{prefix} | Correlation histograms (EDGE-only)",
            out_path=os.path.join(out_dir, f"{prefix}_corr_hist_EDGES.png"),
            dpi=dpi
        )


# Helper: internal plotting helper for boxplot + jitter overlay (kept internal by convention).
def _box_jitter_plot(arrs, labels, title, out_path, rng, dpi=300):
    """
    Internal helper: box + jitter dots overlay for multiple arrays.
    """
    plt.figure(figsize=(10, 5))
    plt.boxplot(arrs, labels=labels, showfliers=True)

    # jittered dots
    for i, a in enumerate(arrs, start=1):
        if len(a) == 0:
            continue
        # downsample dots if huge
        a_plot = a
        if len(a_plot) > 5000:
            idx = rng.choice(len(a_plot), size=5000, replace=False)
            a_plot = a_plot[idx]
        x = rng.normal(loc=i, scale=0.06, size=len(a_plot))
        plt.scatter(x, a_plot, s=4, alpha=0.25)

    plt.title(title)
    plt.ylabel("Correlation value")
    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


# Helper: internal histogram helper for overlaid distributions (kept internal by convention).
def _hist_plot(arrs, labels, title, out_path, dpi=300):
    """
    Internal helper: overlaid histograms for multiple arrays.
    """
    plt.figure(figsize=(10, 5))
    for a, lab in zip(arrs, labels):
        if len(a) == 0:
            continue
        plt.hist(a, bins=60, alpha=0.35, label=lab)
    plt.title(title)
    plt.xlabel("Correlation value")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


# ============================================================
# 9) Clustered heatmaps & dendrograms
# ============================================================

# Helper: cluster a correlation matrix and save standard + clustered heatmaps and dendrogram.
def clustered_heatmap_and_dendrogram(
    corr_df: pd.DataFrame,
    out_dir: str,
    prefix: str,
    method: str = "average",
    dpi: int = 300
):
    """
    Creates:
      - Standard heatmap (original order)
      - Clustered heatmap (rows/cols reordered by hierarchical clustering)
      - Dendrogram plot

    Clustering distance:
      dist = 1 - corr (works for corr in [-1,1])
    """
    os.makedirs(out_dir, exist_ok=True)

    genes = list(corr_df.index)
    M = corr_df.values.astype(float)
    M = np.nan_to_num(M, nan=0.0)

    # distance matrix for clustering: 1 - corr
    dist = 1.0 - M
    np.fill_diagonal(dist, 0.0)
    dist_condensed = squareform(dist, checks=False)

    Z = linkage(dist_condensed, method=method)

    # dendrogram order
    dend = dendrogram(Z, no_plot=True)
    order = dend["leaves"]
    genes_ord = [genes[i] for i in order]
    M_ord = M[np.ix_(order, order)]

    # 1) standard heatmap
    plt.figure(figsize=(8, 7))
    plt.imshow(M, aspect="auto", interpolation="nearest")
    plt.title(f"{prefix} | Heatmap (original order)")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_heatmap.png"), dpi=dpi)
    plt.close()

    # 2) clustered heatmap
    plt.figure(figsize=(8, 7))
    plt.imshow(M_ord, aspect="auto", interpolation="nearest")
    plt.title(f"{prefix} | Heatmap (clustered)")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_heatmap_clustered.png"), dpi=dpi)
    plt.close()

    # 3) dendrogram
    plt.figure(figsize=(10, 4))
    dendrogram(Z, labels=None, leaf_rotation=90)
    plt.title(f"{prefix} | Dendrogram")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_dendrogram.png"), dpi=dpi)
    plt.close()

    return genes_ord, Z


# Helper: scale matrix to [-1, 1] to make distance computation comparable for DIFF matrices.
def scale_to_unit_corr(M: np.ndarray) -> np.ndarray:
    """
    Scale a matrix to [-1, 1] by dividing by max abs value (safe for 0).
    """
    max_abs = np.max(np.abs(M))
    if max_abs == 0 or np.isnan(max_abs):
        return np.zeros_like(M)
    return M / max_abs


# Helper: clustered heatmap+dendrogram for DIFF matrices (TOP - BOTTOM).
def clustered_heatmap_and_dendrogram_diff(
    corr_top: pd.DataFrame,
    corr_bottom: pd.DataFrame,
    out_dir: str,
    prefix: str,
    method: str = "average",
    dpi: int = 300
):
    """
    Same as clustered_heatmap_and_dendrogram but for signed DIFF matrix:
      corr_diff = corr_top - corr_bottom (signed)
    We scale diff to [-1,1] before clustering distance.
    """
    os.makedirs(out_dir, exist_ok=True)

    corr_diff = (corr_top - corr_bottom).copy()
    genes = list(corr_diff.index)
    M = np.nan_to_num(corr_diff.values.astype(float), nan=0.0)

    M_scaled = scale_to_unit_corr(M)  # now in [-1,1]

    dist = 1.0 - M_scaled
    np.fill_diagonal(dist, 0.0)
    dist_condensed = squareform(dist, checks=False)
    Z = linkage(dist_condensed, method=method)

    dend = dendrogram(Z, no_plot=True)
    order = dend["leaves"]
    genes_ord = [genes[i] for i in order]
    M_ord = M[np.ix_(order, order)]  # plot real signed diff (not scaled)

    # standard diff heatmap
    plt.figure(figsize=(8, 7))
    plt.imshow(M, aspect="auto", interpolation="nearest")
    plt.title(f"{prefix} | DIFF Heatmap (TOP-BOT) original order")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_diff_heatmap.png"), dpi=dpi)
    plt.close()

    # clustered diff heatmap
    plt.figure(figsize=(8, 7))
    plt.imshow(M_ord, aspect="auto", interpolation="nearest")
    plt.title(f"{prefix} | DIFF Heatmap (TOP-BOT) clustered")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_diff_heatmap_clustered.png"), dpi=dpi)
    plt.close()

    # dendrogram
    plt.figure(figsize=(10, 4))
    dendrogram(Z, labels=None, leaf_rotation=90)
    plt.title(f"{prefix} | DIFF Dendrogram")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_diff_dendrogram.png"), dpi=dpi)
    plt.close()

    return genes_ord, Z


# Helper: convert linkage clustering output into gene->cluster_id for labeling/plotting.
def cluster_labels_from_linkage(
    genes_in_order: list,
    Z,
    k: int
) -> dict:
    """
    Given hierarchical linkage Z and number of clusters k,
    returns dict: gene -> cluster_id (1..k).
    """
    labels = fcluster(Z, t=k, criterion="maxclust")
    # IMPORTANT: labels returned correspond to the original gene order (not leaves order)
    # Here we assume genes_in_order is the original order used to compute Z? Not leaves.
    # In our heatmap funcs we used original gene list for Z, so pass that gene list there.
    return {g: int(c) for g, c in zip(genes_in_order, labels)}


# Helper: export gene->cluster assignments for downstream inspection.
def save_gene_cluster_labels(cluster_map: dict, out_path: str):
    """
    Save gene->cluster mapping as CSV.
    """
    df = pd.DataFrame({"gene": list(cluster_map.keys()), "cluster": list(cluster_map.values())})
    df.to_csv(out_path, index=False)


# ============================================================
# 10) Cluster-colored network plotting
# ============================================================

# Helper: plot network with node colors based on a cluster mapping (e.g., hierarchical clusters).
def plot_global_network_clusters_fixedpos(
    G: nx.Graph,
    pos: dict,
    out_path: str,
    cluster_map: dict,
    a3_genes: list,
    title: str,
    weighted: bool = True,
    dpi: int = 300,
    label_apobec: bool = True
):
    """
    Plot network with node colors based on cluster labels.
    APOBEC nodes are still labeled (optional).
    Cluster colors are auto-assigned by matplotlib colormap.
    """
    import matplotlib.cm as cm

    nodes = list(G.nodes())
    clusters = np.array([cluster_map.get(n, 0) for n in nodes], dtype=int)

    # map cluster ids -> colors via colormap
    maxc = int(clusters.max()) if clusters.size else 0
    cmap = cm.get_cmap("tab20", max(1, maxc + 1))
    node_colors = [cmap(c) for c in clusters]

    # edge widths
    if weighted and G.number_of_edges() > 0:
        abs_w = np.array([G[u][v].get("abs_weight", 1.0) for u, v in G.edges()])
        if abs_w.max() > abs_w.min():
            widths = 0.5 + 5.0 * (abs_w - abs_w.min()) / (abs_w.max() - abs_w.min())
        else:
            widths = np.ones_like(abs_w) * 2.0
    else:
        widths = 0.8

    plt.figure(figsize=(18, 14))
    nx.draw_networkx_nodes(G, pos, node_size=45, node_color=node_colors, alpha=0.9)
    nx.draw_networkx_edges(G, pos, width=widths, alpha=0.25)

    if label_apobec:
        labels = {n: n for n in nodes if n in a3_genes}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10)

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


# ============================================================
# 11) DIFF graph builders
# ============================================================

# Helper: build a DIFF graph using corr_top - corr_bottom, keeping edges by |diff| threshold.
def build_diff_graph_from_corr(corr_top: pd.DataFrame,
                               corr_bottom: pd.DataFrame,
                               threshold: float,
                               use_abs_diff: bool = True) -> nx.Graph:
    """
    Build a DIFF graph using corr_diff = corr_top - corr_bottom.
    Keep edges where abs(corr_diff) >= threshold (or corr_diff >= threshold if use_abs_diff=False).
    Store:
      - weight = corr_diff (signed)
      - abs_weight = abs(corr_diff)
    Keeps all nodes.
    """
    corr_top = corr_top.copy().fillna(0)
    corr_bottom = corr_bottom.copy().fillna(0)

    corr_diff = (corr_top - corr_bottom).fillna(0)
    genes = list(corr_diff.index)

    G = nx.Graph()
    G.add_nodes_from(genes)

    M = corr_diff.values.astype(float)
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            w = float(M[i, j])
            keep = (abs(w) >= threshold) if use_abs_diff else (w >= threshold)
            if keep and abs(w) < 0.9999999:
                G.add_edge(genes[i], genes[j], weight=w, abs_weight=abs(w))

    return G


# ============================================================
# 12) Community detection (Leiden/Louvain) + layouts
# ============================================================

# Helper: detect communities on a graph using Leiden (preferred) or Louvain (fallback).
def detect_communities_on_graph(G: nx.Graph,
                                method: str = "leiden",
                                seed: int = 42) -> dict:
    """
    Detect communities/modules in a graph.
    Returns: node -> community_id (0..K-1)

    method:
      - "leiden" (requires igraph + leidenalg)
      - "louvain" (requires python-louvain)
    """
    if G.number_of_nodes() == 0:
        return {}

    # if graph has no edges, every node is its own community
    if G.number_of_edges() == 0:
        return {n: i for i, n in enumerate(G.nodes())}

    if method.lower() == "leiden":
        try:
            import igraph as ig
            import leidenalg

            nodes = list(G.nodes())
            idx = {n: i for i, n in enumerate(nodes)}

            edges = []
            weights = []
            for u, v, d in G.edges(data=True):
                edges.append((idx[u], idx[v]))
                # use abs_weight so communities group by strong changes
                weights.append(float(d.get("abs_weight", abs(d.get("weight", 1.0)))))

            g = ig.Graph(n=len(nodes), edges=edges, directed=False)
            g.es["weight"] = weights

            # Leiden partition (modularity)
            part = leidenalg.find_partition(
                g,
                leidenalg.RBConfigurationVertexPartition,
                weights=g.es["weight"],
                seed=seed
            )

            comm_map = {}
            for cid, members in enumerate(part):
                for vid in members:
                    comm_map[nodes[vid]] = cid
            return comm_map

        except Exception:
            # fall back to louvain if leiden fails
            method = "louvain"

    if method.lower() == "louvain":
        try:
            import community as community_louvain  # python-louvain
            # build weight attribute
            H = G.copy()
            for u, v, d in H.edges(data=True):
                d["weight"] = float(d.get("abs_weight", abs(d.get("weight", 1.0))))
            part = community_louvain.best_partition(H, weight="weight", random_state=seed)
            # python-louvain returns node->community_id
            return {n: int(c) for n, c in part.items()}
        except Exception as e:
            raise RuntimeError("Need leidenalg/igraph or python-louvain for community detection.") from e

    raise ValueError(f"Unknown community method: {method}")


# Helper: create a community-separated layout from a DIFF graph and its community map.
def community_aware_layout_from_diff(G_diff: nx.Graph,
                                     comm_map: dict,
                                     seed: int = 42) -> dict:
    """
    Build node positions based ONLY on DIFF graph and its detected communities.
    Produces a module-separated layout.
    Returns: pos dict usable for TOP/BOTTOM/DIFF plots.
    """
    rng = np.random.default_rng(seed)

    nodes = list(G_diff.nodes())
    if len(nodes) == 0:
        return {}

    # If no comm map or trivial, fallback to spring layout on diff
    if not comm_map or len(set(comm_map.values())) <= 1:
        return nx.spring_layout(G_diff if G_diff.number_of_edges() > 0 else nx.Graph(G_diff), seed=seed)

    # group nodes by community
    comm_to_nodes = {}
    for n in nodes:
        c = comm_map.get(n, -1)
        comm_to_nodes.setdefault(c, []).append(n)

    # build community graph (super-nodes)
    CG = nx.Graph()
    for c in comm_to_nodes.keys():
        CG.add_node(c)

    # edge weight between communities = sum abs_weight between nodes
    for u, v, d in G_diff.edges(data=True):
        cu = comm_map.get(u, -1)
        cv = comm_map.get(v, -1)
        if cu == cv:
            continue
        w = float(d.get("abs_weight", abs(d.get("weight", 1.0))))
        if CG.has_edge(cu, cv):
            CG[cu][cv]["weight"] += w
        else:
            CG.add_edge(cu, cv, weight=w)

    # layout community centroids
    pos_comm = nx.spring_layout(CG, seed=seed, weight="weight")

    # layout within each community and scale it down
    pos = {}
    for c, cnodes in comm_to_nodes.items():
        sub = G_diff.subgraph(cnodes).copy()
        # internal layout
        if sub.number_of_edges() > 0 and sub.number_of_nodes() > 1:
            pos_sub = nx.spring_layout(sub, seed=seed, weight="abs_weight")
        else:
            # isolated module: random tiny jitter
            pos_sub = {n: rng.normal(0, 0.01, size=2) for n in cnodes}

        # scale internal layout small
        arr = np.array(list(pos_sub.values()))
        if len(arr) > 0:
            arr = arr - arr.mean(axis=0)
            scale = 0.20  # module size
            pos_sub = {n: scale * v for n, v in zip(pos_sub.keys(), arr)}

        # shift to community centroid
        centroid = pos_comm.get(c, np.array([0.0, 0.0]))
        for n in cnodes:
            pos[n] = pos_sub[n] + centroid

    return pos


# Helper: plot a network with node colors = community id using a fixed pos layout.
def plot_network_by_community_fixedpos(G: nx.Graph,
                                       pos: dict,
                                       comm_map: dict,
                                       out_path: str,
                                       title: str,
                                       a3_genes: list = None,
                                       weighted: bool = True,
                                       dpi: int = 300):
    """
    Plot network using fixed positions and node colors = community.
    APOBEC genes optionally labeled.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    nodes = list(G.nodes())
    if len(nodes) == 0:
        return

    # community ids
    comm_ids = np.array([comm_map.get(n, -1) for n in nodes], dtype=int)

    uniq = sorted(set(comm_ids.tolist()))
    remap = {c: i for i, c in enumerate(uniq)}
    comm_idx = np.array([remap[c] for c in comm_ids], dtype=int)

    cmap = colormaps.get_cmap("tab20").resampled(max(1, len(uniq)))
    node_colors = [cmap(i) for i in comm_idx]

    # edge widths
    if weighted and G.number_of_edges() > 0:
        abs_w = np.array([G[u][v].get("abs_weight", 1.0) for u, v in G.edges()])
        if abs_w.max() > abs_w.min():
            widths = 0.5 + 5.0 * (abs_w - abs_w.min()) / (abs_w.max() - abs_w.min())
        else:
            widths = np.ones_like(abs_w) * 2.0
    else:
        widths = 0.8

    plt.figure(figsize=(18, 14))
    nx.draw_networkx_nodes(G, pos, node_size=45, node_color=node_colors, alpha=0.90)
    nx.draw_networkx_edges(G, pos, width=widths, alpha=0.25)

    if a3_genes:
        labels = {n: n for n in nodes if n in a3_genes}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10)

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


# Helper: community-aware layout (blob layout) given a partition dict, used to separate modules visually.
def community_layout_from_partition(
    G: nx.Graph,
    partition: dict,
    seed: int = 42,
    scale: float = 10.0,
    inner_scale: float = 1.5
):
    """
    Community-aware layout:
    - partition: dict[node] -> community_id
    - Places communities as separate blobs
    - Then places nodes within each blob
    Returns pos dict usable for TOP/BOTTOM/DIFF with fixed positions.
    """
    rng = np.random.default_rng(seed)

    # 1) Build community -> nodes
    comm_to_nodes = {}
    for n, c in partition.items():
        comm_to_nodes.setdefault(c, []).append(n)

    communities = sorted(comm_to_nodes.keys())

    # 2) Build "community graph" (super-nodes)
    CG = nx.Graph()
    CG.add_nodes_from(communities)

    # Add weighted edges between communities based on number of inter-community edges
    for u, v in G.edges():
        cu, cv = partition.get(u), partition.get(v)
        if cu is None or cv is None or cu == cv:
            continue
        if CG.has_edge(cu, cv):
            CG[cu][cv]["weight"] += 1.0
        else:
            CG.add_edge(cu, cv, weight=1.0)

    # 3) Layout community graph (global positions for each module)
    if CG.number_of_edges() > 0 and CG.number_of_nodes() > 1:
        pos_comm = nx.spring_layout(CG, seed=seed, weight="weight")
    else:
        # fallback: spread communities in a grid-ish way
        pos_comm = {}
        k = int(np.ceil(np.sqrt(len(communities))))
        for i, c in enumerate(communities):
            pos_comm[c] = np.array([i % k, i // k], dtype=float)

    # normalize & scale the community centers
    comm_centers = np.array(list(pos_comm.values()), dtype=float)
    comm_centers -= comm_centers.mean(axis=0, keepdims=True)
    max_abs = np.max(np.abs(comm_centers)) if comm_centers.size else 1.0
    if max_abs == 0:
        max_abs = 1.0
    comm_centers = (comm_centers / max_abs) * scale
    for i, c in enumerate(communities):
        pos_comm[c] = comm_centers[i]

    # 4) Layout nodes inside each community, then translate to the community center
    pos = {}
    for c in communities:
        nodes = comm_to_nodes[c]
        sub = G.subgraph(nodes).copy()

        if sub.number_of_nodes() == 1:
            pos[nodes[0]] = pos_comm[c] + rng.normal(0, 0.01, size=2)
            continue

        # local layout inside community
        if sub.number_of_edges() > 0:
            local = nx.spring_layout(sub, seed=seed, weight="abs_weight")
        else:
            # no edges: just scatter locally
            local = {n: rng.normal(0, 0.2, size=2) for n in sub.nodes()}

        # normalize local coords then scale
        local_coords = np.array(list(local.values()), dtype=float)
        local_coords -= local_coords.mean(axis=0, keepdims=True)
        local_max = np.max(np.abs(local_coords)) if local_coords.size else 1.0
        if local_max == 0:
            local_max = 1.0
        local_coords = (local_coords / local_max) * inner_scale

        for i, n in enumerate(sub.nodes()):
            pos[n] = pos_comm[c] + local_coords[i]

    return pos


# ============================================================
# 13) nx9 helper functions (cleanup, safety, merging)
# ============================================================

# Helper: keep only the largest connected component to focus on the main network structure.
def nx9_largest_component(G: nx.Graph) -> nx.Graph:
    """Return largest connected component; if no edges, return copy."""
    if G.number_of_edges() == 0:
        return G.copy()
    comp = max(nx.connected_components(G), key=len)
    return G.subgraph(comp).copy()


# Helper: remove degree-0 nodes (isolates) to simplify plots/community detection.
def nx9_remove_degree0(G: nx.Graph) -> nx.Graph:
    """Remove nodes with degree 0."""
    H = G.copy()
    H.remove_nodes_from([n for n in H.nodes() if H.degree(n) == 0])
    return H


# Helper: merge small communities into "Other" to avoid many tiny modules in plots.
def nx9_merge_communities_to_topK(comm_map: dict, k_keep: int = 5, min_size: int = 40):
    """
    Merge small communities into an 'Other' bin to avoid thousands of tiny groups.
    - keep top k_keep communities by size (and size >= min_size)
    - everything else -> Other community id (k_keep)
    Returns: merged_map, raw_sizes, merged_sizes
    """
    # community sizes
    sizes = {}
    for n, cid in comm_map.items():
        sizes[cid] = sizes.get(cid, 0) + 1

    # keep only the biggest k_keep
    top = sorted(sizes.keys(), key=lambda c: sizes[c], reverse=True)[:k_keep]
    top = set(top)

    other_id = k_keep
    merged = {}
    for n, cid in comm_map.items():
        if (cid in top) and (sizes[cid] >= min_size):
            merged[n] = cid
        else:
            merged[n] = other_id

    merged_sizes = {}
    for n, cid in merged.items():
        merged_sizes[cid] = merged_sizes.get(cid, 0) + 1

    return merged, sizes, merged_sizes


# Helper: ensure pos contains coordinates for every node (important after subgraphing/unioning).
def nx9_safe_fill_pos(pos: dict, G: nx.Graph, seed: int = 42):
    """Ensure pos has all nodes in G (add missing at origin)."""
    for n in G.nodes():
        if n not in pos:
            pos[n] = np.array([0.0, 0.0])
    return pos


# ============================================================
# 14) Community-aware layout v2 (parameterized)
# ============================================================

# Helper: parameterized community-aware layout for more controllable inter/intra module spacing.
def community_aware_layout_from_diff_v2(
    G, partition, seed=42,
    inter_scale=4.0,
    intra_scale=1.0,
    k_local=None,
    k_global=None,
    iters_local=200,
    iters_global=200,
    edge_weight_attr="abs_weight",
):
    """
    Community-aware layout computed from a DIFF graph G and a partition dict {node -> community_id}.
    Returns pos: dict {node -> np.array([x,y])}

    - Inter-community structure from a super-graph (communities as nodes)
    - Intra-community structure from each community induced subgraph
    - Final: pos[node] = inter_scale * center[comm] + intra_scale * local[node]
    """
    rng = np.random.default_rng(seed)

    nodes = list(G.nodes())
    if len(nodes) == 0:
        return {}

    # If partition is missing or only 1 community, fall back to spring layout on G
    if not partition:
        return nx.spring_layout(G, seed=seed, weight="weight", k=k_global, iterations=iters_global)

    comm_ids = {partition.get(n, -1) for n in nodes}
    if len(comm_ids) <= 1:
        return nx.spring_layout(G, seed=seed, weight="weight", k=k_global, iterations=iters_global)

    # Group nodes by community
    comm_to_nodes = {}
    for n in nodes:
        c = partition.get(n, -1)
        comm_to_nodes.setdefault(c, []).append(n)

    # Build community graph CG
    CG = nx.Graph()
    for c in comm_to_nodes.keys():
        CG.add_node(c)

    # Edge weight between communities = sum of abs(edge_weight) across cross-community edges
    for u, v, d in G.edges(data=True):
        cu = partition.get(u, -1)
        cv = partition.get(v, -1)
        if cu == cv:
            continue

        w = d.get(edge_weight_attr, None)
        if w is None:
            w = d.get("weight", 1.0)
        try:
            w = float(w)
        except Exception:
            w = 1.0
        w = abs(w)

        if CG.has_edge(cu, cv):
            CG[cu][cv]["weight"] += w
        else:
            CG.add_edge(cu, cv, weight=w)

    # If CG has no edges, place comms on random points
    if CG.number_of_edges() == 0:
        pos_comm = {c: rng.normal(0, 1.0, size=2) for c in CG.nodes()}
    else:
        pos_comm = nx.spring_layout(
            CG, seed=seed, weight="weight",
            k=k_global, iterations=iters_global
        )

    # Convert to numpy arrays and apply inter_scale
    pos_comm = {c: inter_scale * np.asarray(p, dtype=float) for c, p in pos_comm.items()}

    # Local layouts per community
    pos = {}

    for c, cnodes in comm_to_nodes.items():
        sub = G.subgraph(cnodes).copy()

        if sub.number_of_nodes() == 1:
            only = cnodes[0]
            local = np.zeros(2, dtype=float)
            pos[only] = pos_comm.get(c, np.zeros(2)) + intra_scale * local
            continue

        if sub.number_of_edges() > 0:
            weight_attr = edge_weight_attr if any(edge_weight_attr in dd for _, _, dd in sub.edges(data=True)) else "weight"
            pos_sub = nx.spring_layout(
                sub, seed=seed, weight=weight_attr,
                k=k_local, iterations=iters_local
            )
            local_coords = np.array([pos_sub[n] for n in cnodes], dtype=float)
        else:
            local_coords = rng.normal(0, 0.05, size=(len(cnodes), 2))

        # Center local coords at (0,0)
        local_coords = local_coords - local_coords.mean(axis=0, keepdims=True)

        # Apply intra_scale
        local_coords = intra_scale * local_coords

        center = pos_comm.get(c, np.zeros(2, dtype=float))

        for i, n in enumerate(cnodes):
            pos[n] = center + local_coords[i]

    return pos


# ============================================================
# 15) Heatmap with community bar
# ============================================================

# Helper: render a matrix heatmap with an aligned community bar and optional separators/count labels.
def plot_heatmap_with_community_bar(mat_df, gene_order, comm_map, out_path, title="",
                                   draw_separators=True, annotate_counts=True):
    # reorder matrix
    genes = [g for g in gene_order if g in mat_df.index]
    M = mat_df.loc[genes, genes].values

    # build community vector aligned to genes
    comms = [comm_map.get(g, "Other") for g in genes]

    # find community blocks (consecutive segments)
    boundaries = []
    blocks = []
    start = 0
    for i in range(1, len(comms) + 1):
        if i == len(comms) or comms[i] != comms[i-1]:
            end = i
            blocks.append((comms[i-1], start, end))
            boundaries.append(end)
            start = i

    # heatmap
    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[20, 1], hspace=0.05)

    ax = fig.add_subplot(gs[0])
    vmax = np.nanmax(np.abs(M)) if M.size else 1.0
    if vmax == 0: vmax = 1.0

    im = ax.imshow(M, aspect="auto", interpolation="nearest", vmin=-vmax, vmax=vmax)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])

    # separators on heatmap
    if draw_separators:
        for _, s, e in blocks[:-1]:
            ax.axhline(e - 0.5, linewidth=1.2)
            ax.axvline(e - 0.5, linewidth=1.2)

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # community bar
    axb = fig.add_subplot(gs[1])
    # encode communities into integers
    uniq = list(dict.fromkeys(comms))
    cmap_idx = {c: i for i, c in enumerate(uniq)}
    bar = np.array([cmap_idx[c] for c in comms])[None, :]
    axb.imshow(bar, aspect="auto", interpolation="nearest")
    axb.set_xticks([])
    axb.set_yticks([])

    # separators on bar + annotate counts
    if draw_separators:
        for _, s, e in blocks[:-1]:
            axb.axvline(e - 0.5, linewidth=1.2)

    if annotate_counts:
        # put labels centered in each block
        for c, s, e in blocks:
            mid = (s + e) / 2
            count = e - s
            axb.text(mid, 0, f"{c} (n={count})", ha="center", va="center", fontsize=8, rotation=0)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


# ============================================================
# 16) SBS2 vs A3 selection plot
# ============================================================

# Helper: visualize the donor selection logic for High-A3 then High/Low SBS2 split.
def plot_sbs_vs_a3_selection(
        df,
        a3_thr,
        sbs_low_thr,
        sbs_high_thr,
        group1_idx,
        group2_idx,
        out_path,
        title=""
):
    """
    Visualize SBS2 vs A3_score with selection cutoffs.

    Top panel:
        SBS2 (Y) vs A3_score (X)
        - vertical line: A3_score cutoff
        - horizontal lines: SBS2 20% and 80% (within High-A3)
        - colored selected groups

    Bottom panel:
        Distribution of A3_score
        - histogram + KDE
        - mean + median lines
        - cutoff line
    """
    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])

    ax_scatter = fig.add_subplot(gs[0])
    ax_hist = fig.add_subplot(gs[1])

    # TOP PANEL — Scatter
    ax_scatter.scatter(
        df["A3_score"], df["SBS2"],
        s=20, alpha=0.3, color="gray", label="All samples"
    )

    # highlight groups
    ax_scatter.scatter(
        df.loc[group1_idx, "A3_score"],
        df.loc[group1_idx, "SBS2"],
        s=40, color="red", label="Group1 (High SBS2 in High-A3)"
    )

    ax_scatter.scatter(
        df.loc[group2_idx, "A3_score"],
        df.loc[group2_idx, "SBS2"],
        s=40, color="blue", label="Group2 (Low SBS2 in High-A3)"
    )

    # Cutoff lines
    ax_scatter.axvline(a3_thr, color="black", linestyle="--", label="A3 median cutoff")
    ax_scatter.axhline(sbs_high_thr, color="red", linestyle="--", label="SBS2 80% (High-A3)")
    ax_scatter.axhline(sbs_low_thr, color="blue", linestyle="--", label="SBS2 20% (High-A3)")

    ax_scatter.set_xlabel("A3_score")
    ax_scatter.set_ylabel("SBS2")
    ax_scatter.set_title(title)
    ax_scatter.legend(loc="best")
    ax_scatter.grid(alpha=0.2)

    # BOTTOM PANEL — A3 Distribution
    sns.histplot(df["A3_score"], bins=40, kde=True, ax=ax_hist, color="gray")

    mean_val = df["A3_score"].mean()
    median_val = df["A3_score"].median()

    ax_hist.axvline(mean_val, color="green", linestyle="-", label=f"Mean = {mean_val:.3f}")
    ax_hist.axvline(median_val, color="black", linestyle="--", label=f"Median = {median_val:.3f}")
    ax_hist.axvline(a3_thr, color="red", linestyle=":", label=f"Cutoff = {a3_thr:.3f}")

    ax_hist.set_xlabel("A3_score")
    ax_hist.set_ylabel("Count")
    ax_hist.legend(loc="best")
    ax_hist.grid(alpha=0.2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


# ============================================================
# 17) Community ordering helpers + simple heatmaps
# ============================================================

# Helper: order genes by community size (descending), alphabetical within community, "Other" last.
def nx9_order_genes_by_community(comm_map, genes=None, other_label="Other"):
    """
    Return a community-ordered gene list.
    - Communities ordered by size descending.
    - Genes within each community sorted alphabetically (stable).
    - "Other" community (if present) placed last.
    """
    if genes is None:
        genes = list(comm_map.keys())
    genes = [g for g in genes if g in comm_map]

    # community -> genes
    comm_to_genes = {}
    for g in genes:
        c = comm_map.get(g, other_label)
        comm_to_genes.setdefault(c, []).append(g)

    # order communities by size (desc), but put "Other" last
    comm_items = list(comm_to_genes.items())
    comm_items.sort(key=lambda kv: len(kv[1]), reverse=True)

    if other_label in comm_to_genes:
        comm_items = [kv for kv in comm_items if kv[0] != other_label] + [(other_label, comm_to_genes[other_label])]

    ordered = []
    for c, glist in comm_items:
        ordered.extend(sorted(glist))  # stable inside module

    return ordered


# Helper: safely subset a matrix to ordered genes (intersection only) to avoid KeyErrors.
def nx9_restrict_and_order_matrix(mat_df, ordered_genes):
    """Safely subset matrix to ordered_genes (intersection only) and return ordered matrix + final gene list."""
    g = [x for x in ordered_genes if x in mat_df.index]
    g = [x for x in g if x in mat_df.columns]
    return mat_df.loc[g, g], g


# Helper: minimal heatmap helper for quick inspection of community-ordered matrices.
def nx9_plot_simple_heatmap(mat_df, out_path, title, vmax=None):
    arr = mat_df.values
    # optional symmetric scale for DIFF
    if vmax is None:
        vmax = np.nanmax(np.abs(arr)) if arr.size else 1.0
        if vmax == 0:
            vmax = 1.0

    plt.figure(figsize=(10, 9))
    plt.imshow(arr, aspect="auto", interpolation="nearest", vmin=-vmax, vmax=vmax)
    plt.colorbar(label="value")
    plt.title(title)
    plt.xlabel("genes (community order)")
    plt.ylabel("genes (community order)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


# ============================================================
# 18) Modularity + community detection v2
# ============================================================

# Helper: convert node->community dict to list-of-sets format for modularity computations.
def nx9_partition_to_communities(comm_map):
    """Convert {node: community_id} to list-of-sets (community sets)."""
    comm_to_nodes = {}
    for n, c in comm_map.items():
        comm_to_nodes.setdefault(c, set()).add(n)
    return list(comm_to_nodes.values())


# Helper: compute modularity for a given partition on a NetworkX graph.
def nx9_modularity_nx(G, comm_map, weight="weight"):
    """Compute modularity using NetworkX quality.modularity."""
    from networkx.algorithms.community.quality import modularity
    comms = nx9_partition_to_communities(comm_map)
    if len(comms) == 0:
        return np.nan
    return float(modularity(G, comms, weight=weight))


# Helper: updated community detection wrapper with resolution + optional modularity output.
def detect_communities_on_graph_v2(
    G,
    method="leiden",
    seed=42,
    resolution=1.0,
    weight="weight",
    return_modularity=False
):
    """
    NEW version (keeps old function untouched).
    Supports resolution and optional modularity output.

    Returns:
      - comm_map: dict {node: community_id}
      - modularity (optional): float
    """
    if G is None or G.number_of_nodes() == 0:
        if return_modularity:
            return {}, np.nan
        return {}

    method = (method or "").lower()

    # 1) LEIDEN (igraph + leidenalg)
    if method == "leiden":
        try:
            import igraph as ig
            import leidenalg as la

            nodes = list(G.nodes())
            idx = {n: i for i, n in enumerate(nodes)}

            edges = []
            weights = []
            for u, v, d in G.edges(data=True):
                edges.append((idx[u], idx[v]))
                w = float(d.get(weight, 1.0))
                weights.append(w)

            igG = ig.Graph(n=len(nodes), edges=edges, directed=False)
            if len(weights) == len(edges) and len(edges) > 0:
                igG.es["weight"] = weights
                part = la.find_partition(
                    igG,
                    la.RBConfigurationVertexPartition,
                    weights="weight",
                    resolution_parameter=float(resolution),
                    seed=int(seed)
                )
            else:
                part = la.find_partition(
                    igG,
                    la.RBConfigurationVertexPartition,
                    resolution_parameter=float(resolution),
                    seed=int(seed)
                )

            membership = part.membership
            comm_map = {nodes[i]: int(membership[i]) for i in range(len(nodes))}

            if return_modularity:
                mod = nx9_modularity_nx(G, comm_map, weight=weight)
                return comm_map, mod
            return comm_map

        except Exception:
            method = "louvain"

    # 2) LOUVAIN (python-louvain)
    if method == "louvain":
        try:
            import community as community_louvain  # python-louvain
            comm_map = community_louvain.best_partition(
                G, weight=weight, resolution=float(resolution), random_state=int(seed)
            )
            comm_map = {k: int(v) for k, v in comm_map.items()}

            if return_modularity:
                mod = nx9_modularity_nx(G, comm_map, weight=weight)
                return comm_map, mod
            return comm_map
        except Exception as e:
            raise RuntimeError(
                "Could not run Louvain. Install python-louvain: pip install python-louvain"
            ) from e

    raise ValueError(f"Unknown method='{method}'. Use 'leiden' or 'louvain'.")


# ============================================================
# 19) Community-aware layout v3 + zoomed community plots
# ============================================================

# Helper: tuned community-aware layout (v3) with tighter inter-module spacing and larger intra spread.
def community_aware_layout_from_diff_v3(
    G, comm_map, seed=42,
    inter_scale=1.8,
    intra_scale=2.5,
    k_local=None,
    k_global=None,
    iters_local=200,
    iters_global=200,
    weight_attr="abs_weight"
):
    rng = np.random.default_rng(seed)

    nodes = list(G.nodes())
    if len(nodes) == 0:
        return {}

    # fallback if comm_map missing / trivial
    if (not comm_map) or len(set(comm_map.values())) <= 1:
        return nx.spring_layout(G if G.number_of_edges() > 0 else nx.Graph(G), seed=seed)

    # group nodes by community
    comm_to_nodes = {}
    for n in nodes:
        c = comm_map.get(n, "Other")
        comm_to_nodes.setdefault(c, []).append(n)

    # build community graph (super-nodes)
    CG = nx.Graph()
    for c in comm_to_nodes:
        CG.add_node(c)

    for u, v, d in G.edges(data=True):
        cu = comm_map.get(u, "Other")
        cv = comm_map.get(v, "Other")
        if cu == cv:
            continue
        w = float(d.get(weight_attr, abs(d.get("weight", 1.0))))
        if CG.has_edge(cu, cv):
            CG[cu][cv]["weight"] += w
        else:
            CG.add_edge(cu, cv, weight=w)

    # community centroid layout (global)
    pos_comm = nx.spring_layout(
        CG,
        seed=seed,
        weight="weight",
        k=k_global,
        iterations=iters_global
    )

    # scale centroids to control inter-community distance
    for c in pos_comm:
        pos_comm[c] = inter_scale * np.array(pos_comm[c], dtype=float)

    # layout within each community
    pos = {}
    for c, cnodes in comm_to_nodes.items():
        sub = G.subgraph(cnodes).copy()

        if sub.number_of_nodes() == 1:
            pos_sub = {cnodes[0]: np.array([0.0, 0.0])}
        elif sub.number_of_edges() > 0:
            pos_sub = nx.spring_layout(
                sub,
                seed=seed,
                weight=weight_attr if any(weight_attr in sub[u][v] for u, v in sub.edges()) else "weight",
                k=k_local,
                iterations=iters_local
            )
        else:
            pos_sub = {n: rng.normal(0, 0.05, size=2) for n in cnodes}

        # center & scale inside community
        arr = np.array([pos_sub[n] for n in cnodes], dtype=float)
        arr = arr - arr.mean(axis=0, keepdims=True)
        arr = intra_scale * arr
        pos_sub = {n: arr[i] for i, n in enumerate(cnodes)}

        centroid = np.array(pos_comm.get(c, [0.0, 0.0]), dtype=float)
        for n in cnodes:
            pos[n] = pos_sub[n] + centroid

    return pos


# Helper: generate per-community zoomed-in plots using a global layout and bounding boxes.
def plot_zoomed_communities(
    G, pos, comm_map, out_dir, title_prefix,
    a3_genes=None, pad=0.25,
    min_nodes=10,
    weighted=True,
    label_apobec=True
):
    os.makedirs(out_dir, exist_ok=True)
    a3_genes = set(a3_genes or [])

    # group nodes by community
    comm_to_nodes = {}
    for n in G.nodes():
        c = comm_map.get(n, "Other")
        comm_to_nodes.setdefault(c, []).append(n)

    # stable order: largest communities first
    comm_order = sorted(comm_to_nodes.keys(), key=lambda c: len(comm_to_nodes[c]), reverse=True)

    for c in comm_order:
        nodes_c = comm_to_nodes[c]
        if len(nodes_c) < min_nodes:
            continue

        sub = G.subgraph(nodes_c).copy()

        # ensure APOBEC nodes are included if they exist in the full graph AND this community contains them
        a3_in = [g for g in nodes_c if g in a3_genes]

        # get bounding box from global pos for those nodes
        pts = np.array([pos[n] for n in nodes_c if n in pos], dtype=float)
        if pts.size == 0:
            continue

        xmin, ymin = pts.min(axis=0)
        xmax, ymax = pts.max(axis=0)

        dx = (xmax - xmin) if xmax > xmin else 1.0
        dy = (ymax - ymin) if ymax > ymin else 1.0

        xmin -= pad * dx; xmax += pad * dx
        ymin -= pad * dy; ymax += pad * dy

        plt.figure(figsize=(8, 8))
        ax = plt.gca()

        # draw edges
        nx.draw_networkx_edges(sub, pos, ax=ax, alpha=0.25)

        # draw nodes
        nx.draw_networkx_nodes(sub, pos, ax=ax, node_size=30, alpha=0.9)

        # highlight APOBEC
        if a3_in:
            nx.draw_networkx_nodes(
                sub, pos, ax=ax,
                nodelist=a3_in, node_size=140, alpha=0.95
            )
            if label_apobec:
                nx.draw_networkx_labels(
                    sub, pos,
                    labels={g: g for g in a3_in},
                    font_size=10
                )

        ax.set_title(f"{title_prefix} | Community={c} | n={len(nodes_c)}")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.axis("off")

        out_path = os.path.join(out_dir, f"{title_prefix}_COMM_{str(c)}_zoom.png")
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()


# ============================================================
# 20) Resolution sweep stability
# ============================================================

# Helper: sweep resolution parameter and assess stability across multiple random seeds.
def nx9_resolution_sweep_stability(
    G,
    out_dir,
    method="leiden",
    resolutions=None,
    n_runs=10,
    base_seed=42,
    weight="abs_weight",
    verbose=True
):
    """
    Runs community detection multiple times per resolution.
    Saves:
      - sweep_summary.csv
      - partitions per (resolution, seed)
      - plots: ncomms vs resolution, modularity vs resolution, ARI vs resolution, NMI vs resolution
    """
    import itertools
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    os.makedirs(out_dir, exist_ok=True)
    part_dir = os.path.join(out_dir, "partitions")
    plot_dir = os.path.join(out_dir, "plots")
    os.makedirs(part_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    if resolutions is None:
        resolutions = np.round(np.linspace(0.2, 2.0, 10), 2)

    nodes_eval = list(G.nodes())

    def partition_to_labels(cm):
        comm_ids = [cm.get(n, -1) for n in nodes_eval]
        uniq = {c: i for i, c in enumerate(sorted(set(comm_ids)))}
        return np.array([uniq[c] for c in comm_ids], dtype=int)

    rows = []

    for r in resolutions:
        labels_list = []
        modularities = []
        ncomms_list = []

        for i in range(n_runs):
            seed = base_seed + i
            cm, mod = detect_communities_on_graph_v2(
                G, method=method, seed=seed, resolution=float(r),
                weight=weight, return_modularity=True
            )

            # save partition
            part_path = os.path.join(part_dir, f"partition__res{r:.2f}__seed{seed}.json")
            with open(part_path, "w") as f:
                json.dump(cm, f)

            labels_list.append(partition_to_labels(cm))
            ncomms_list.append(len(set(cm.values())))
            modularities.append(mod if mod is not None else np.nan)

        # pairwise stability
        aris, nmis = [], []
        for a, b in itertools.combinations(range(n_runs), 2):
            aris.append(adjusted_rand_score(labels_list[a], labels_list[b]))
            nmis.append(normalized_mutual_info_score(labels_list[a], labels_list[b]))

        rows.append({
            "resolution": float(r),
            "n_runs": int(n_runs),
            "ncomms_mean": float(np.mean(ncomms_list)),
            "ncomms_std": float(np.std(ncomms_list)),
            "modularity_mean": float(np.nanmean(modularities)),
            "modularity_std": float(np.nanstd(modularities)),
            "ARI_mean": float(np.mean(aris)),
            "ARI_std": float(np.std(aris)),
            "NMI_mean": float(np.mean(nmis)),
            "NMI_std": float(np.std(nmis)),
        })

    sweep_df = pd.DataFrame(rows).sort_values("resolution")
    sweep_csv = os.path.join(out_dir, "sweep_summary.csv")
    sweep_df.to_csv(sweep_csv, index=False)

    # plot 1: ncomms only
    plt.figure(figsize=(8, 5))
    x = sweep_df["resolution"].values
    y = sweep_df["ncomms_mean"].values
    yerr = sweep_df["ncomms_std"].values
    plt.plot(x, y, marker="o")
    plt.fill_between(x, y - yerr, y + yerr, alpha=0.2)
    plt.xlabel("resolution")
    plt.ylabel("#communities (mean ± std)")
    plt.title("Resolution sweep: #communities")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "ncomms_only.png"), dpi=300)
    plt.close()

    # plot 2: ncomms colored by ARI
    plt.figure(figsize=(8, 5))
    sc = plt.scatter(
        sweep_df["resolution"], sweep_df["ncomms_mean"],
        c=sweep_df["ARI_mean"], s=60
    )
    plt.colorbar(sc, label="ARI_mean")
    plt.xlabel("resolution")
    plt.ylabel("#communities (mean)")
    plt.title("Resolution sweep: #communities colored by stability (ARI)")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "ncomms_colored_by_ARI.png"), dpi=300)
    plt.close()

    for col, fname, ylabel in [
        ("ARI_mean", "ARI_vs_resolution.png", "ARI (mean)"),
        ("NMI_mean", "NMI_vs_resolution.png", "NMI (mean)"),
        ("modularity_mean", "modularity_vs_resolution.png", "Modularity (mean)"),
    ]:
        plt.figure(figsize=(8, 5))
        plt.plot(sweep_df["resolution"], sweep_df[col], marker="o")
        plt.xlabel("resolution")
        plt.ylabel(ylabel)
        plt.title(f"Resolution sweep: {ylabel}")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, fname), dpi=300)
        plt.close()

    return sweep_df, sweep_csv, part_dir, plot_dir


# ============================================================
# 21) Violin correlations
# ============================================================

# Helper: quick violin plot comparing TOP/BOTTOM/DIFF correlation distributions.
def plot_corr_violin(corr_top, corr_bottom, corr_diff, out_path_prefix):
    def upper_vals(mat):
        a = mat.values
        iu = np.triu_indices_from(a, k=1)
        return a[iu]

    top_vals = upper_vals(corr_top)
    bot_vals = upper_vals(corr_bottom)
    diff_vals = upper_vals(corr_diff)

    plt.figure(figsize=(8, 4))
    plt.violinplot([top_vals, bot_vals, diff_vals], showmeans=False, showmedians=True, showextrema=True)
    plt.xticks([1, 2, 3], ["TOP", "BOTTOM", "DIFF (TOP-BOT)"])
    plt.ylabel("Spearman correlation")
    plt.title("Correlation distributions (violin)")
    plt.tight_layout()
    plt.savefig(f"{out_path_prefix}_violin.png", dpi=300)
    plt.close()


# ============================================================
# 22) Save partition outputs
# ============================================================

# Helper: save community partition map and a per-community gene list, plus a colored plot.
def save_partition_outputs(G, pos, comm_map, out_dir, tag):
    os.makedirs(out_dir, exist_ok=True)

    # 1) JSON mapping
    json_path = os.path.join(out_dir, f"{tag}_comm_map.json")
    with open(json_path, "w") as f:
        json.dump(comm_map, f, indent=2)

    # 2) CSV gene list per community
    comm_to_genes = {}
    for g, c in comm_map.items():
        comm_to_genes.setdefault(str(c), []).append(g)

    rows = []
    for c, genes in sorted(comm_to_genes.items(), key=lambda x: -len(x[1])):
        rows.append({
            "community": c,
            "size": len(genes),
            "genes": ";".join(sorted(genes))
        })

    csv_path = os.path.join(out_dir, f"{tag}_community_geneList.csv")
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    # 3) Plot (colored by community)
    # If you already have plot_network_by_community_fixedpos, use it:
    from network_utils import plot_network_by_community_fixedpos
    plot_path = os.path.join(out_dir, f"{tag}_community_plot.png")
    plot_network_by_community_fixedpos(
        G, pos, comm_map,
        out_path=plot_path,
        title=tag,
        a3_genes=["A3A","A3B","A3C","A3D","A3F","A3G","A3H"],
        weighted=True
    )


# ============================================================
# 23) DIFF graph strategies (TOPK, soft WGCNA-like)
# ============================================================

# Helper: build DIFF graph by keeping top-k strongest diff edges per node (no global threshold).
def build_diff_graph_topk(corr_top: pd.DataFrame,
                          corr_bottom: pd.DataFrame,
                          k: int = 10,
                          use_abs_diff: bool = True) -> nx.Graph:
    """
    Build a DIFF graph without a global threshold:
    For each node i, connect it to the TOP-k nodes by |diff(i,j)|.
    """
    genes = list(corr_top.index)
    diff = (corr_top - corr_bottom).fillna(0.0)
    if use_abs_diff:
        score = diff.abs()
    else:
        score = diff

    G = nx.Graph()
    G.add_nodes_from(genes)

    for i, g in enumerate(genes):
        # take topk neighbors excluding self
        s = score.loc[g].copy()
        s[g] = -np.inf
        top_neighbors = s.nlargest(k).index.tolist()

        for h in top_neighbors:
            w_signed = float(diff.loc[g, h])
            w_abs = float(abs(w_signed))
            if not np.isfinite(w_abs) or w_abs <= 0:
                continue
            if G.has_edge(g, h):
                # keep strongest abs
                if w_abs > G[g][h].get("abs_weight", 0):
                    G[g][h]["weight"] = w_signed
                    G[g][h]["abs_weight"] = w_abs
            else:
                G.add_edge(g, h, weight=w_signed, abs_weight=w_abs)

    return G


# Helper: WGCNA-like soft adjacency builder for a single correlation matrix.
def build_soft_wgcna_graph(corr: pd.DataFrame,
                           power: int = 6,
                           keep_topk_per_node: int = 30) -> nx.Graph:
    """
    WGCNA-like adjacency:
      adjacency(i,j) = |corr(i,j)|^power
    """
    genes = list(corr.index)
    A = corr.abs().fillna(0.0) ** power
    np.fill_diagonal(A.values, 0.0)

    G = nx.Graph()
    G.add_nodes_from(genes)

    for g in genes:
        s = A.loc[g]
        top_neighbors = s.nlargest(keep_topk_per_node).index.tolist()
        for h in top_neighbors:
            w = float(A.loc[g, h])
            if w <= 0 or not np.isfinite(w):
                continue
            if g == h:
                continue
            if not G.has_edge(g, h):
                G.add_edge(g, h, weight=w, abs_weight=w)
            else:
                if w > G[g][h].get("abs_weight", 0):
                    G[g][h]["weight"] = w
                    G[g][h]["abs_weight"] = w
    return G


# Helper: build a DIFF graph from soft adjacency matrices (|corr|^power) and keep top-k per node by |diff|.
def build_diff_graph_soft_wgcna(corr_top: pd.DataFrame,
                                corr_bottom: pd.DataFrame,
                                power: int = 6,
                                topk_per_node: int = 30) -> nx.Graph:
    """
    Create soft adjacency for each group, then DIFF as difference in adjacency.
    For community detection, we keep abs(diff) as abs_weight.
    """
    A_top = corr_top.abs().fillna(0.0) ** power
    A_bot = corr_bottom.abs().fillna(0.0) ** power
    D = (A_top - A_bot).fillna(0.0)
    genes = list(D.index)

    G = nx.Graph()
    G.add_nodes_from(genes)
    np.fill_diagonal(D.values, 0.0)

    for g in genes:
        s = D.loc[g].abs()
        s[g] = -np.inf
        neigh = s.nlargest(topk_per_node).index.tolist()
        for h in neigh:
            w_signed = float(D.loc[g, h])
            w_abs = float(abs(w_signed))
            if w_abs <= 0 or not np.isfinite(w_abs):
                continue
            if not G.has_edge(g, h):
                G.add_edge(g, h, weight=w_signed, abs_weight=w_abs)
            else:
                if w_abs > G[g][h].get("abs_weight", 0):
                    G[g][h]["weight"] = w_signed
                    G[g][h]["abs_weight"] = w_abs
    return G


# ============================================================
# 24) Volcano & Manhattan
# ============================================================

# Helper: volcano plot for differential statistics tables (log2FC vs -log10(p)).
def plot_volcano(stats_df, outpath, title=None, q_thr=0.05):
    """Volcano plot: log2FC vs -log10(p-value)."""
    if stats_df.empty:
        return
    df = stats_df.copy()
    df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))

    sig_mask = df["qvalue"] <= q_thr
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df["log2FC_high_vs_low"], df["neglog10p"],
               c="lightgray", s=10, alpha=0.7, label="NS")
    ax.scatter(df.loc[sig_mask, "log2FC_high_vs_low"],
               df.loc[sig_mask, "neglog10p"],
               c="red", s=12, alpha=0.8, label=f"FDR ≤ {q_thr:g}")

    ax.axhline(-np.log10(0.05), color="black", ls="--", lw=1, alpha=0.6)
    ax.axvline(0.0, color="black", ls=":", lw=1, alpha=0.6)

    ax.set_xlabel("log2FC (High SBS2 vs Low SBS2)")
    ax.set_ylabel("-log10(p-value)")
    if title:
        ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# Helper: manhattan-like plot using gene index as pseudo genomic position.
def plot_manhattan_index(stats_df, outpath, title=None, q_thr=0.05):
    """
    Manhattan-like plot using gene index as pseudo-position.
    """
    if stats_df.empty:
        return
    df = stats_df.copy().reset_index(drop=True)
    df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    df["idx"] = np.arange(len(df))

    sig_mask = df["qvalue"] <= q_thr

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.scatter(df["idx"], df["neglog10p"],
               c="gray", s=6, alpha=0.6, label="NS")
    ax.scatter(df.loc[sig_mask, "idx"],
               df.loc[sig_mask, "neglog10p"],
               c="red", s=8, alpha=0.8, label=f"FDR ≤ {q_thr:g}")

    ax.axhline(-np.log10(0.05), color="black", ls="--", lw=1, alpha=0.6)

    ax.set_xlabel("Gene index")
    ax.set_ylabel("-log10(p-value)")
    if title:
        ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# ============================================================
# 25) Network saving helpers (STEP 8)
# ============================================================

# Helper: remove isolates (degree 0 nodes) from a graph before exporting/plotting.
def remove_isolated_nodes(G: nx.Graph) -> nx.Graph:
    H = G.copy()
    isolates = list(nx.isolates(H))
    H.remove_nodes_from(isolates)
    return H


# Helper: export a graph as a CSV-like edge list (weighted or unweighted).
def save_graph_edgelist(G: nx.Graph, outpath: str, weighted: bool = True):
    """Save graph as edge list (CSV-like)."""
    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f)
        if weighted:
            writer.writerow(["source", "target", "weight"])
            for u, v, d in G.edges(data=True):
                w = float(d.get("weight", 0.0))
                writer.writerow([u, v, w])
        else:
            writer.writerow(["source", "target"])
            for u, v in G.edges():
                writer.writerow([u, v])


# ============================================================
# 26) Misc utilities (ENS cleanup, banners, presence checks, etc.)
# ============================================================

# Helper: strip Ensembl version suffix (".N") from ENSG identifiers.
def clean_ensg(x: str) -> str:
    # remove trailing ".digits" version
    return re.sub(r"\.\d+$", "", str(x))


# Helper: print a visible console banner block for pipeline stage separation.
def banner(title: str, char: str = "=", width: int = 100):
    line = char * width
    print("\n" + line)
    print(title)
    print(line)


# Helper: shorten long Entity_ID strings for readability (first 4 dash-separated parts).
def shorten_entity_id(x):
    return "-".join(str(x).split("-")[:4])


# Helper: ensure every edge has abs_weight, useful if some graphs were created without it.
def ensure_abs_weight(G):
    for u, v, d in G.edges(data=True):
        w = float(d.get("weight", 0.0))
        d["abs_weight"] = abs(float(d.get("abs_weight", w)))


# Helper: report which target genes are present/missing in a graph (quick biomarker sanity check).
def report_gene_presence(tag, G, genes):
    present = [g for g in genes if g in G.nodes]
    missing = [g for g in genes if g not in G.nodes]
    print(f"\n[{tag}] present: {len(present)} | missing: {len(missing)}")
    if present:
        print("  present:", present)
    if missing:
        print("  missing:", missing)


# Helper: write both weighted and unweighted edge lists for the same graph (convenience export).
def save_edges_both_formats(G, out_prefix):
    save_graph_edgelist(G, out_prefix + "_weighted.csv", weighted=True)
    save_graph_edgelist(G, out_prefix + "_unweighted.csv", weighted=False)


# Helper: convert node->community map into integer labels aligned to a node list (for ARI/NMI comparisons).
def partition_to_labels(cm: dict, nodes: list):
    """Convert {node:community} -> integer labels aligned to nodes list."""
    ids = [cm.get(n, -1) for n in nodes]
    uniq = {c: i for i, c in enumerate(sorted(set(ids)))}
    return np.array([uniq[c] for c in ids], dtype=int)


# Helper: save a simple curve plot to disk (used in sweeps).
def save_curve(x, y, title, ylabel, out_path):
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, marker="o")
    plt.xlabel("resolution")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


# Helper: export per-community gene lists as CSV (one row per community).
def save_gene_lists_from_map(cm_map: dict, out_csv: str):
    """Save one row per community with gene list."""
    comm_to_genes = {}
    for g, c in cm_map.items():
        comm_to_genes.setdefault(c, []).append(g)

    rows = []
    for c, genes in comm_to_genes.items():
        rows.append({"community": c, "size": len(genes), "genes": ";".join(sorted(genes))})

    pd.DataFrame(rows).sort_values("size", ascending=False).to_csv(out_csv, index=False)




























# """
# network_utils.py
#
# Utility functions for:
# - Building weighted gene-gene networks from correlations
# - Computing APOBEC A3 scores
# - Plotting QC figures and network visualizations
# - Computing and visualizing network metrics
# - Enforcing fixed layouts and filtering isolated nodes
#
# Designed to be imported into the main analysis pipeline.
# """
#
# import os
# import numpy as np
# import pandas as pd
# import networkx as nx
# import matplotlib.pyplot as plt
# from matplotlib.offsetbox import AnchoredText
#
#
# from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
# from scipy.spatial.distance import squareform
#
# # ============================================================
# # Logging
# # ============================================================
# def log(msg: str, verbose: bool = True):
#     """
#     Simple console logger.
#
#     Parameters
#     ----------
#     msg : str
#         Message to print.
#     verbose : bool
#         If False, suppress output.
#     """
#     if verbose:
#         print(msg, flush=True)
#
#
# # ============================================================
# # Graph construction & export
# # ============================================================
# def build_weighted_graph_from_corr(corr_df: pd.DataFrame, threshold: float) -> nx.Graph:
#     """
#     Build an undirected weighted graph from a correlation matrix.
#
#     - Nodes = genes
#     - Edges kept if |correlation| >= threshold
#     - Keeps isolated nodes
#     - Edge attributes:
#         * weight     : signed correlation
#         * abs_weight : absolute correlation
#
#     Parameters
#     ----------
#     corr_df : pd.DataFrame
#         Correlation matrix (genes x genes).
#     threshold : float
#         Minimum absolute correlation to keep an edge.
#
#     Returns
#     -------
#     nx.Graph
#         Weighted gene network.
#     """
#     corr_df = corr_df.copy().fillna(0)
#     corr_df = corr_df.apply(pd.to_numeric, errors="coerce").fillna(0)
#
#     genes = list(corr_df.index)
#     G = nx.Graph()
#     G.add_nodes_from(genes)
#
#     for i in range(len(genes)):
#         for j in range(i + 1, len(genes)):
#             g1, g2 = genes[i], genes[j]
#             w = float(corr_df.iat[i, j])
#             if abs(w) >= threshold and abs(w) < 0.999999:
#                 G.add_edge(g1, g2, weight=w, abs_weight=abs(w))
#     return G
#
#
# def save_graph_csv(G: nx.Graph, out_dir: str, out_prefix: str):
#     """
#     Save network structure and node metrics to CSV files.
#
#     Outputs:
#     - <prefix>_edges.csv : source, target, weight, abs_weight
#     - <prefix>_nodes.csv : node, degree, betweenness
#
#     Parameters
#     ----------
#     G : nx.Graph
#         Network to export.
#     out_dir : str
#         Output directory.
#     out_prefix : str
#         Filename prefix.
#     """
#     os.makedirs(out_dir, exist_ok=True)
#
#     edges = []
#     for u, v, data in G.edges(data=True):
#         edges.append({
#             "source": u,
#             "target": v,
#             "weight": float(data.get("weight", 1.0)),
#             "abs_weight": float(data.get("abs_weight", abs(data.get("weight", 1.0))))
#         })
#     pd.DataFrame(edges).to_csv(os.path.join(out_dir, f"{out_prefix}_edges.csv"), index=False)
#
#     deg = dict(G.degree())
#     btw = nx.betweenness_centrality(G) if G.number_of_edges() > 0 else {n: 0.0 for n in G.nodes()}
#
#     nodes_df = pd.DataFrame({
#         "node": list(G.nodes()),
#         "degree": [deg[n] for n in G.nodes()],
#         "betweenness": [btw[n] for n in G.nodes()],
#     })
#     nodes_df.to_csv(os.path.join(out_dir, f"{out_prefix}_nodes.csv"), index=False)
#
#
# # ============================================================
# # APOBEC A3 score & QC plots
# # ============================================================
# def compute_a3_score_AB(df_in: pd.DataFrame, a_col="A3A", b_col="A3B") -> pd.Series:
#     """
#     Compute APOBEC A3 score using ONLY A3A and A3B.
#
#     - Each gene is normalized by its max value
#       within the cancer type
#     - A3_score = normalized A3A + normalized A3B
#
#     Parameters
#     ----------
#     df_in : pd.DataFrame
#         DataFrame containing A3A and A3B columns.
#
#     Returns
#     -------
#     pd.Series
#         A3_score per sample (donor).
#     """
#     tmp = df_in[[a_col, b_col]].copy()
#     tmp[a_col] = pd.to_numeric(tmp[a_col], errors="coerce")
#     tmp[b_col] = pd.to_numeric(tmp[b_col], errors="coerce")
#
#     a_max = tmp[a_col].max()
#     b_max = tmp[b_col].max()
#
#     a_norm = 0.0 if (pd.isna(a_max) or a_max == 0) else (tmp[a_col] / a_max)
#     b_norm = 0.0 if (pd.isna(b_max) or b_max == 0) else (tmp[b_col] / b_max)
#
#     return (a_norm + b_norm).astype(float)
#
#
# def plot_a3_scatter_and_distribution(cancer_df: pd.DataFrame, out_dir: str, cancer_type: str):
#     """
#     Generate QC plots for APOBEC A3 activity.
#
#     Plots:
#     1) Scatter: A3_score vs SBS2 (each dot = donor)
#     2) Histogram: distribution of A3_score across donors
#
#     Parameters
#     ----------
#     cancer_df : pd.DataFrame
#         DataFrame containing A3_score and SBS2.
#     out_dir : str
#         Directory to save plots.
#     cancer_type : str
#         Cancer name for plot titles.
#     """
#     os.makedirs(out_dir, exist_ok=True)
#
#     plt.figure(figsize=(7, 6))
#     plt.scatter(cancer_df["A3_score"], cancer_df["SBS2"], alpha=0.6)
#     plt.xlabel("A3_score (A3A + A3B)")
#     plt.ylabel("SBS2")
#     plt.title(f"{cancer_type} | donor-level scatter")
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{cancer_type}_scatter_A3score_vs_SBS2.png"), dpi=300)
#     plt.close()
#
#     plt.figure(figsize=(7, 5))
#     plt.hist(cancer_df["A3_score"].dropna().values, bins=40)
#     plt.xlabel("A3_score")
#     plt.ylabel("Donor count")
#     plt.title(f"{cancer_type} | A3_score distribution")
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{cancer_type}_dist_A3score.png"), dpi=300)
#     plt.close()
#
#
# # ============================================================
# # Network metrics
# # ============================================================
# def compute_graph_metrics_df(G: nx.Graph) -> pd.DataFrame:
#     """
#     Compute per-gene network metrics.
#
#     Metrics:
#     - degree
#     - betweenness centrality
#     - strength_abs (sum of |edge weights|)
#     - closeness centrality
#
#     Parameters
#     ----------
#     G : nx.Graph
#         Gene network.
#
#     Returns
#     -------
#     pd.DataFrame
#         One row per gene with metrics.
#     """
#     nodes = list(G.nodes())
#     deg = dict(G.degree())
#     betw = nx.betweenness_centrality(G) if G.number_of_edges() > 0 else {n: 0.0 for n in nodes}
#     close = nx.closeness_centrality(G) if G.number_of_edges() > 0 else {n: 0.0 for n in nodes}
#
#     strength = {}
#     for n in nodes:
#         s = 0.0
#         for _, data in G[n].items():
#             s += float(data.get("abs_weight", abs(data.get("weight", 1.0))))
#         strength[n] = s
#
#     return pd.DataFrame({
#         "node": nodes,
#         "degree": [deg[n] for n in nodes],
#         "betweenness": [betw[n] for n in nodes],
#         "strength_abs": [strength[n] for n in nodes],
#         "closeness": [close[n] for n in nodes],
#     })
#
#
# def plot_metric_distributions(metrics_df: pd.DataFrame, out_dir: str, prefix: str):
#     """
#     Plot distributions of network metrics.
#
#     Helps decide correlation thresholds by visual inspection.
#
#     Parameters
#     ----------
#     metrics_df : pd.DataFrame
#         Output of compute_graph_metrics_df().
#     out_dir : str
#         Directory for plots.
#     prefix : str
#         Filename prefix.
#     """
#     os.makedirs(out_dir, exist_ok=True)
#
#     for col in ["degree", "betweenness", "strength_abs", "closeness"]:
#         plt.figure(figsize=(7, 5))
#         vals = pd.to_numeric(metrics_df[col], errors="coerce").dropna().values
#         plt.hist(vals, bins=40)
#         plt.xlabel(col)
#         plt.ylabel("Gene count")
#         plt.title(f"{prefix} | {col} distribution")
#         plt.tight_layout()
#         plt.savefig(os.path.join(out_dir, f"{prefix}_dist_{col}.png"), dpi=300)
#         plt.close()
#
#
# # ============================================================
# # Layout control & node filtering
# # ============================================================
# def get_constant_layout_from_bottom(G_bottom: nx.Graph, G_union: nx.Graph, seed: int = 42) -> dict:
#     """
#     Generate a fixed node layout shared across multiple plots.
#
#     Strategy:
#     - Prefer layout computed from bottom network
#     - Ensure all union nodes receive coordinates
#
#     Parameters
#     ----------
#     G_bottom : nx.Graph
#         Bottom-group network.
#     G_union : nx.Graph
#         Union of top and bottom networks.
#
#     Returns
#     -------
#     dict
#         Node -> (x, y) layout dictionary.
#     """
#     base = G_bottom if G_bottom.number_of_edges() > 0 else G_union
#     pos = nx.spring_layout(base, seed=seed)
#
#     for n in G_union.nodes():
#         if n not in pos:
#             pos[n] = np.array([0.0, 0.0])
#     return pos
#
#
# def remove_isolates_in_both(G_top: nx.Graph, G_bottom: nx.Graph, G_any: nx.Graph) -> nx.Graph:
#     """
#     Remove genes that are isolated in BOTH top and bottom networks.
#
#     Parameters
#     ----------
#     G_top : nx.Graph
#         Top network.
#     G_bottom : nx.Graph
#         Bottom network.
#     G_any : nx.Graph
#         Graph to filter.
#
#     Returns
#     -------
#     nx.Graph
#         Filtered graph without double-isolated nodes.
#     """
#     keep = []
#     for n in G_any.nodes():
#         if (G_top.degree(n) > 0) or (G_bottom.degree(n) > 0):
#             keep.append(n)
#     return G_any.subgraph(keep).copy()
#
#
# # ============================================================
# # Network plotting (fixed positions)
# from matplotlib.offsetbox import AnchoredText
# import matplotlib.pyplot as plt
# import numpy as np
# import networkx as nx
#
# def plot_global_network_fixedpos(
#     G: nx.Graph,
#     pos: dict,
#     out_path: str,
#     a3_genes: list,
#     title: str,
#     weighted: bool = True,
#     dpi: int = 300,
#     label_apobec: bool = True,
#     legend_text=None,
#     legend_loc="upper left"
# ):
#     node_colors = ["red" if n in a3_genes else "gray" for n in G.nodes()]
#
#     if weighted and G.number_of_edges() > 0:
#         abs_w = np.array([float(G[u][v].get("abs_weight", 1.0)) for u, v in G.edges()])
#         if abs_w.max() > abs_w.min():
#             widths = 0.5 + 5.0 * (abs_w - abs_w.min()) / (abs_w.max() - abs_w.min())
#         else:
#             widths = np.ones_like(abs_w) * 2.0
#     else:
#         widths = 0.8
#
#     # ✅ use explicit fig/ax
#     fig, ax = plt.subplots(figsize=(18, 14))
#
#     nx.draw_networkx_nodes(G, pos, ax=ax, node_size=40, node_color=node_colors, alpha=0.9)
#     nx.draw_networkx_edges(G, pos, ax=ax, width=widths, alpha=0.30)
#
#     if label_apobec:
#         labels = {n: n for n in G.nodes() if n in a3_genes}
#         nx.draw_networkx_labels(G, pos, ax=ax, labels=labels, font_size=10)
#
#     ax.set_title(title)
#     ax.axis("off")
#
#     # ✅ legend textbox
#     if legend_text is not None and str(legend_text).strip() != "":
#         at = AnchoredText(legend_text, loc=legend_loc, prop={"size": 8}, frameon=True)
#         at.patch.set_alpha(0.85)
#         ax.add_artist(at)
#
#     # ✅ IMPORTANT: prevents legend cropping
#     fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
#     plt.close(fig)
#
#
# def plot_combined_network_fixedpos(G_top, G_bottom, pos, out_path, a3_genes, title, dpi=300):
#     """
#     Plot combined top + bottom networks using fixed positions.
#
#     Edge coloring:
#     - Red  : edges unique to top group
#     - Blue : edges unique to bottom group
#     - Gray : edges shared between both
#
#     Parameters
#     ----------
#     G_top, G_bottom : nx.Graph
#         Networks to compare.
#     """
#     G_union = nx.compose(G_top, G_bottom)
#
#     top_edges = set(tuple(sorted(e)) for e in G_top.edges())
#     bot_edges = set(tuple(sorted(e)) for e in G_bottom.edges())
#     only_top = top_edges - bot_edges
#     only_bot = bot_edges - top_edges
#     both = top_edges & bot_edges
#
#     node_colors = ["red" if n in a3_genes else "gray" for n in G_union.nodes()]
#
#     plt.figure(figsize=(18, 14))
#     nx.draw_networkx_nodes(G_union, pos, node_size=40, node_color=node_colors, alpha=0.9)
#     nx.draw_networkx_edges(G_union, pos, edgelist=list(both), alpha=0.10)
#     nx.draw_networkx_edges(G_union, pos, edgelist=list(only_top), edge_color="red", alpha=0.35)
#     nx.draw_networkx_edges(G_union, pos, edgelist=list(only_bot), edge_color="blue", alpha=0.35)
#
#     labels = {n: n for n in G_union.nodes() if n in a3_genes}
#     nx.draw_networkx_labels(G_union, pos, labels=labels, font_size=10)
#
#     plt.title(title)
#     plt.axis("off")
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=dpi)
#     plt.close()
#
#
#
# def bh_fdr(pvals: np.ndarray) -> np.ndarray:
#     """
#     Benjamini–Hochberg FDR correction.
#     Returns q-values (same length as pvals).
#     """
#     p = np.asarray(pvals, dtype=float)
#     n = len(p)
#     order = np.argsort(p)
#     ranked = p[order]
#     q = ranked * n / (np.arange(1, n + 1))
#     q = np.minimum.accumulate(q[::-1])[::-1]
#     q = np.clip(q, 0, 1)
#     out = np.empty_like(q)
#     out[order] = q
#     return out
#
#
#
# def corr_upper_values(corr_df: pd.DataFrame, drop_diag: bool = True) -> np.ndarray:
#     """
#     Extract upper-triangle correlation values from a square correlation matrix.
#     Returns a 1D numpy array of values (optionally excludes diagonal).
#     """
#     M = corr_df.values.astype(float)
#     n = M.shape[0]
#     k = 1 if drop_diag else 0
#     vals = M[np.triu_indices(n, k=k)]
#     vals = vals[~np.isnan(vals)]
#     return vals
#
#
# def corr_edge_values(corr_df: pd.DataFrame, threshold: float) -> np.ndarray:
#     """
#     Extract only correlation values whose absolute value >= threshold (excluding diagonal).
#     Useful to compare 'all correlations' vs 'edges we keep' distributions.
#     """
#     M = corr_df.values.astype(float)
#     n = M.shape[0]
#     vals = M[np.triu_indices(n, k=1)]
#     vals = vals[~np.isnan(vals)]
#     vals = vals[np.abs(vals) >= threshold]
#     return vals
#
#
# def plot_corr_distributions(
#     corr_top: pd.DataFrame,
#     corr_bottom: pd.DataFrame,
#     out_dir: str,
#     prefix: str,
#     threshold: float = None,
#     seed: int = 42,
#     dpi: int = 300
# ):
#     """
#     Creates distribution plots for:
#       - all correlations (upper triangle, excluding diagonal)
#       - optional: edge-only correlations where |corr| >= threshold
#
#     Outputs:
#       - <prefix>_corr_box_jitter_ALL.png
#       - <prefix>_corr_box_jitter_EDGES.png (if threshold is not None)
#       - <prefix>_corr_hist_ALL.png
#       - <prefix>_corr_hist_EDGES.png (if threshold is not None)
#     """
#     os.makedirs(out_dir, exist_ok=True)
#     rng = np.random.default_rng(seed)
#
#     corr_diff = corr_top - corr_bottom  # signed difference
#
#     # ---- ALL values ----
#     vals_top = corr_upper_values(corr_top)
#     vals_bot = corr_upper_values(corr_bottom)
#     vals_diff = corr_upper_values(corr_diff)
#
#     _box_jitter_plot(
#         [vals_top, vals_bot, vals_diff],
#         labels=["TOP", "BOTTOM", "DIFF (TOP-BOT)"],
#         title=f"{prefix} | Correlation distributions (ALL pairs)",
#         out_path=os.path.join(out_dir, f"{prefix}_corr_box_jitter_ALL.png"),
#         rng=rng,
#         dpi=dpi
#     )
#
#     _hist_plot(
#         [vals_top, vals_bot, vals_diff],
#         labels=["TOP", "BOTTOM", "DIFF (TOP-BOT)"],
#         title=f"{prefix} | Correlation histograms (ALL pairs)",
#         out_path=os.path.join(out_dir, f"{prefix}_corr_hist_ALL.png"),
#         dpi=dpi
#     )
#
#     # ---- EDGE-only values ----
#     if threshold is not None:
#         e_top = corr_edge_values(corr_top, threshold)
#         e_bot = corr_edge_values(corr_bottom, threshold)
#         e_diff = corr_edge_values(corr_diff, threshold)  # edges by |diff| threshold (optional meaning)
#         # Note: for diff, edge-only by |diff|>=thr might be too strict; still useful
#
#         _box_jitter_plot(
#             [e_top, e_bot, e_diff],
#             labels=[f"TOP | |corr|≥{threshold}", f"BOTTOM | |corr|≥{threshold}", f"DIFF | |diff|≥{threshold}"],
#             title=f"{prefix} | Correlation distributions (EDGE-only)",
#             out_path=os.path.join(out_dir, f"{prefix}_corr_box_jitter_EDGES.png"),
#             rng=rng,
#             dpi=dpi
#         )
#
#         _hist_plot(
#             [e_top, e_bot, e_diff],
#             labels=[f"TOP | |corr|≥{threshold}", f"BOTTOM | |corr|≥{threshold}", f"DIFF | |diff|≥{threshold}"],
#             title=f"{prefix} | Correlation histograms (EDGE-only)",
#             out_path=os.path.join(out_dir, f"{prefix}_corr_hist_EDGES.png"),
#             dpi=dpi
#         )
#
#
# def _box_jitter_plot(arrs, labels, title, out_path, rng, dpi=300):
#     """
#     Internal helper: box + jitter dots overlay for multiple arrays.
#     """
#     plt.figure(figsize=(10, 5))
#     plt.boxplot(arrs, labels=labels, showfliers=True)
#
#     # jittered dots
#     for i, a in enumerate(arrs, start=1):
#         if len(a) == 0:
#             continue
#         # downsample dots if huge
#         a_plot = a
#         if len(a_plot) > 5000:
#             idx = rng.choice(len(a_plot), size=5000, replace=False)
#             a_plot = a_plot[idx]
#         x = rng.normal(loc=i, scale=0.06, size=len(a_plot))
#         plt.scatter(x, a_plot, s=4, alpha=0.25)
#
#     plt.title(title)
#     plt.ylabel("Correlation value")
#     plt.xticks(rotation=15, ha="right")
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=dpi)
#     plt.close()
#
#
# def _hist_plot(arrs, labels, title, out_path, dpi=300):
#     """
#     Internal helper: overlaid histograms for multiple arrays.
#     """
#     plt.figure(figsize=(10, 5))
#     for a, lab in zip(arrs, labels):
#         if len(a) == 0:
#             continue
#         plt.hist(a, bins=60, alpha=0.35, label=lab)
#     plt.title(title)
#     plt.xlabel("Correlation value")
#     plt.ylabel("Count")
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=dpi)
#     plt.close()
#
#
#
# def clustered_heatmap_and_dendrogram(
#     corr_df: pd.DataFrame,
#     out_dir: str,
#     prefix: str,
#     method: str = "average",
#     dpi: int = 300
# ):
#     """
#     Creates:
#       - Standard heatmap (original order)
#       - Clustered heatmap (rows/cols reordered by hierarchical clustering)
#       - Dendrogram plot
#
#     Clustering distance:
#       dist = 1 - corr (works for corr in [-1,1])
#     """
#     os.makedirs(out_dir, exist_ok=True)
#
#     genes = list(corr_df.index)
#     M = corr_df.values.astype(float)
#     M = np.nan_to_num(M, nan=0.0)
#
#     # distance matrix for clustering: 1 - corr
#     dist = 1.0 - M
#     np.fill_diagonal(dist, 0.0)
#     dist_condensed = squareform(dist, checks=False)
#
#     Z = linkage(dist_condensed, method=method)
#
#     # dendrogram order
#     dend = dendrogram(Z, no_plot=True)
#     order = dend["leaves"]
#     genes_ord = [genes[i] for i in order]
#     M_ord = M[np.ix_(order, order)]
#
#     # 1) standard heatmap
#     plt.figure(figsize=(8, 7))
#     plt.imshow(M, aspect="auto", interpolation="nearest")
#     plt.title(f"{prefix} | Heatmap (original order)")
#     plt.colorbar()
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{prefix}_heatmap.png"), dpi=dpi)
#     plt.close()
#
#     # 2) clustered heatmap
#     plt.figure(figsize=(8, 7))
#     plt.imshow(M_ord, aspect="auto", interpolation="nearest")
#     plt.title(f"{prefix} | Heatmap (clustered)")
#     plt.colorbar()
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{prefix}_heatmap_clustered.png"), dpi=dpi)
#     plt.close()
#
#     # 3) dendrogram
#     plt.figure(figsize=(10, 4))
#     dendrogram(Z, labels=None, leaf_rotation=90)
#     plt.title(f"{prefix} | Dendrogram")
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{prefix}_dendrogram.png"), dpi=dpi)
#     plt.close()
#
#     return genes_ord, Z
#
#
#
#
# def scale_to_unit_corr(M: np.ndarray) -> np.ndarray:
#     """
#     Scale a matrix to [-1, 1] by dividing by max abs value (safe for 0).
#     """
#     max_abs = np.max(np.abs(M))
#     if max_abs == 0 or np.isnan(max_abs):
#         return np.zeros_like(M)
#     return M / max_abs
#
#
# def clustered_heatmap_and_dendrogram_diff(
#     corr_top: pd.DataFrame,
#     corr_bottom: pd.DataFrame,
#     out_dir: str,
#     prefix: str,
#     method: str = "average",
#     dpi: int = 300
# ):
#     """
#     Same as clustered_heatmap_and_dendrogram but for signed DIFF matrix:
#       corr_diff = corr_top - corr_bottom (signed)
#     We scale diff to [-1,1] before clustering distance.
#     """
#     os.makedirs(out_dir, exist_ok=True)
#
#     corr_diff = (corr_top - corr_bottom).copy()
#     genes = list(corr_diff.index)
#     M = np.nan_to_num(corr_diff.values.astype(float), nan=0.0)
#
#     M_scaled = scale_to_unit_corr(M)  # now in [-1,1]
#
#     dist = 1.0 - M_scaled
#     np.fill_diagonal(dist, 0.0)
#     dist_condensed = squareform(dist, checks=False)
#     Z = linkage(dist_condensed, method=method)
#
#     dend = dendrogram(Z, no_plot=True)
#     order = dend["leaves"]
#     genes_ord = [genes[i] for i in order]
#     M_ord = M[np.ix_(order, order)]  # plot real signed diff (not scaled)
#
#     # standard diff heatmap
#     plt.figure(figsize=(8, 7))
#     plt.imshow(M, aspect="auto", interpolation="nearest")
#     plt.title(f"{prefix} | DIFF Heatmap (TOP-BOT) original order")
#     plt.colorbar()
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{prefix}_diff_heatmap.png"), dpi=dpi)
#     plt.close()
#
#     # clustered diff heatmap
#     plt.figure(figsize=(8, 7))
#     plt.imshow(M_ord, aspect="auto", interpolation="nearest")
#     plt.title(f"{prefix} | DIFF Heatmap (TOP-BOT) clustered")
#     plt.colorbar()
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{prefix}_diff_heatmap_clustered.png"), dpi=dpi)
#     plt.close()
#
#     # dendrogram
#     plt.figure(figsize=(10, 4))
#     dendrogram(Z, labels=None, leaf_rotation=90)
#     plt.title(f"{prefix} | DIFF Dendrogram")
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f"{prefix}_diff_dendrogram.png"), dpi=dpi)
#     plt.close()
#
#     return genes_ord, Z
#
#
#
# def cluster_labels_from_linkage(
#     genes_in_order: list,
#     Z,
#     k: int
# ) -> dict:
#     """
#     Given hierarchical linkage Z and number of clusters k,
#     returns dict: gene -> cluster_id (1..k).
#     """
#     labels = fcluster(Z, t=k, criterion="maxclust")
#     # IMPORTANT: labels returned correspond to the original gene order (not leaves order)
#     # Here we assume genes_in_order is the original order used to compute Z? Not leaves.
#     # In our heatmap funcs we used original gene list for Z, so pass that gene list there.
#     return {g: int(c) for g, c in zip(genes_in_order, labels)}
#
#
# def save_gene_cluster_labels(cluster_map: dict, out_path: str):
#     """
#     Save gene->cluster mapping as CSV.
#     """
#     df = pd.DataFrame({"gene": list(cluster_map.keys()), "cluster": list(cluster_map.values())})
#     df.to_csv(out_path, index=False)
#
#
#
#
# def plot_global_network_clusters_fixedpos(
#     G: nx.Graph,
#     pos: dict,
#     out_path: str,
#     cluster_map: dict,
#     a3_genes: list,
#     title: str,
#     weighted: bool = True,
#     dpi: int = 300,
#     label_apobec: bool = True
# ):
#     """
#     Plot network with node colors based on cluster labels.
#     APOBEC nodes are still labeled (optional).
#     Cluster colors are auto-assigned by matplotlib colormap.
#     """
#     import matplotlib.cm as cm
#
#     nodes = list(G.nodes())
#     clusters = np.array([cluster_map.get(n, 0) for n in nodes], dtype=int)
#
#     # map cluster ids -> colors via colormap
#     maxc = int(clusters.max()) if clusters.size else 0
#     cmap = cm.get_cmap("tab20", max(1, maxc + 1))
#     node_colors = [cmap(c) for c in clusters]
#
#     # edge widths
#     if weighted and G.number_of_edges() > 0:
#         abs_w = np.array([G[u][v].get("abs_weight", 1.0) for u, v in G.edges()])
#         if abs_w.max() > abs_w.min():
#             widths = 0.5 + 5.0 * (abs_w - abs_w.min()) / (abs_w.max() - abs_w.min())
#         else:
#             widths = np.ones_like(abs_w) * 2.0
#     else:
#         widths = 0.8
#
#     plt.figure(figsize=(18, 14))
#     nx.draw_networkx_nodes(G, pos, node_size=45, node_color=node_colors, alpha=0.9)
#     nx.draw_networkx_edges(G, pos, width=widths, alpha=0.25)
#
#     if label_apobec:
#         labels = {n: n for n in nodes if n in a3_genes}
#         nx.draw_networkx_labels(G, pos, labels=labels, font_size=10)
#
#     plt.title(title)
#     plt.axis("off")
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=dpi)
#     plt.close()
#
#
#
#
# def build_diff_graph_from_corr(corr_top: pd.DataFrame,
#                                corr_bottom: pd.DataFrame,
#                                threshold: float,
#                                use_abs_diff: bool = True) -> nx.Graph:
#     """
#     Build a DIFF graph using corr_diff = corr_top - corr_bottom.
#     Keep edges where abs(corr_diff) >= threshold (or corr_diff >= threshold if use_abs_diff=False).
#     Store:
#       - weight = corr_diff (signed)
#       - abs_weight = abs(corr_diff)
#     Keeps all nodes.
#     """
#     corr_top = corr_top.copy().fillna(0)
#     corr_bottom = corr_bottom.copy().fillna(0)
#
#     corr_diff = (corr_top - corr_bottom).fillna(0)
#     genes = list(corr_diff.index)
#
#     G = nx.Graph()
#     G.add_nodes_from(genes)
#
#     M = corr_diff.values.astype(float)
#     for i in range(len(genes)):
#         for j in range(i + 1, len(genes)):
#             w = float(M[i, j])
#             keep = (abs(w) >= threshold) if use_abs_diff else (w >= threshold)
#             if keep and abs(w) < 0.9999999:
#                 G.add_edge(genes[i], genes[j], weight=w, abs_weight=abs(w))
#
#     return G
#
#
# def detect_communities_on_graph(G: nx.Graph,
#                                 method: str = "leiden",
#                                 seed: int = 42) -> dict:
#     """
#     Detect communities/modules in a graph.
#     Returns: node -> community_id (0..K-1)
#
#     method:
#       - "leiden" (requires igraph + leidenalg)
#       - "louvain" (requires python-louvain)
#     """
#     if G.number_of_nodes() == 0:
#         return {}
#
#     # if graph has no edges, every node is its own community
#     if G.number_of_edges() == 0:
#         return {n: i for i, n in enumerate(G.nodes())}
#
#     if method.lower() == "leiden":
#         try:
#             import igraph as ig
#             import leidenalg
#
#             nodes = list(G.nodes())
#             idx = {n: i for i, n in enumerate(nodes)}
#
#             edges = []
#             weights = []
#             for u, v, d in G.edges(data=True):
#                 edges.append((idx[u], idx[v]))
#                 # use abs_weight so communities group by strong changes
#                 weights.append(float(d.get("abs_weight", abs(d.get("weight", 1.0)))))
#
#             g = ig.Graph(n=len(nodes), edges=edges, directed=False)
#             g.es["weight"] = weights
#
#             # Leiden partition (modularity)
#             part = leidenalg.find_partition(
#                 g,
#                 leidenalg.RBConfigurationVertexPartition,
#                 weights=g.es["weight"],
#                 seed=seed
#             )
#
#             comm_map = {}
#             for cid, members in enumerate(part):
#                 for vid in members:
#                     comm_map[nodes[vid]] = cid
#             return comm_map
#
#         except Exception:
#             # fall back to louvain if leiden fails
#             method = "louvain"
#
#     if method.lower() == "louvain":
#         try:
#             import community as community_louvain  # python-louvain
#             # build weight attribute
#             H = G.copy()
#             for u, v, d in H.edges(data=True):
#                 d["weight"] = float(d.get("abs_weight", abs(d.get("weight", 1.0))))
#             part = community_louvain.best_partition(H, weight="weight", random_state=seed)
#             # python-louvain returns node->community_id
#             return {n: int(c) for n, c in part.items()}
#         except Exception as e:
#             raise RuntimeError("Need leidenalg/igraph or python-louvain for community detection.") from e
#
#     raise ValueError(f"Unknown community method: {method}")
#
#
#
# def community_aware_layout_from_diff(G_diff: nx.Graph,
#                                      comm_map: dict,
#                                      seed: int = 42) -> dict:
#     """
#     Build node positions based ONLY on DIFF graph and its detected communities.
#     Produces a module-separated layout.
#     Returns: pos dict usable for TOP/BOTTOM/DIFF plots.
#     """
#     rng = np.random.default_rng(seed)
#
#     nodes = list(G_diff.nodes())
#     if len(nodes) == 0:
#         return {}
#
#     # If no comm map or trivial, fallback to spring layout on diff
#     if not comm_map or len(set(comm_map.values())) <= 1:
#         return nx.spring_layout(G_diff if G_diff.number_of_edges() > 0 else nx.Graph(G_diff), seed=seed)
#
#     # group nodes by community
#     comm_to_nodes = {}
#     for n in nodes:
#         c = comm_map.get(n, -1)
#         comm_to_nodes.setdefault(c, []).append(n)
#
#     # build community graph (super-nodes)
#     CG = nx.Graph()
#     for c in comm_to_nodes.keys():
#         CG.add_node(c)
#
#     # edge weight between communities = sum abs_weight between nodes
#     for u, v, d in G_diff.edges(data=True):
#         cu = comm_map.get(u, -1)
#         cv = comm_map.get(v, -1)
#         if cu == cv:
#             continue
#         w = float(d.get("abs_weight", abs(d.get("weight", 1.0))))
#         if CG.has_edge(cu, cv):
#             CG[cu][cv]["weight"] += w
#         else:
#             CG.add_edge(cu, cv, weight=w)
#
#     # layout community centroids
#     pos_comm = nx.spring_layout(CG, seed=seed, weight="weight")
#
#     # layout within each community and scale it down
#     pos = {}
#     for c, cnodes in comm_to_nodes.items():
#         sub = G_diff.subgraph(cnodes).copy()
#         # internal layout
#         if sub.number_of_edges() > 0 and sub.number_of_nodes() > 1:
#             pos_sub = nx.spring_layout(sub, seed=seed, weight="abs_weight")
#         else:
#             # isolated module: random tiny jitter
#             pos_sub = {n: rng.normal(0, 0.01, size=2) for n in cnodes}
#
#         # scale internal layout small
#         arr = np.array(list(pos_sub.values()))
#         if len(arr) > 0:
#             arr = arr - arr.mean(axis=0)
#             scale = 0.20  # module size
#             pos_sub = {n: scale * v for n, v in zip(pos_sub.keys(), arr)}
#
#         # shift to community centroid
#         centroid = pos_comm.get(c, np.array([0.0, 0.0]))
#         for n in cnodes:
#             pos[n] = pos_sub[n] + centroid
#
#     return pos
#
#
#
# def plot_network_by_community_fixedpos(G: nx.Graph,
#                                        pos: dict,
#                                        comm_map: dict,
#                                        out_path: str,
#                                        title: str,
#                                        a3_genes: list = None,
#                                        weighted: bool = True,
#                                        dpi: int = 300):
#     """
#     Plot network using fixed positions and node colors = community.
#     APOBEC genes optionally labeled.
#     """
#     import matplotlib.cm as cm
#
#     os.makedirs(os.path.dirname(out_path), exist_ok=True)
#
#     nodes = list(G.nodes())
#     if len(nodes) == 0:
#         return
#
#     # community ids
#     comm_ids = np.array([comm_map.get(n, -1) for n in nodes], dtype=int)
#     import matplotlib.cm as cm
#     from matplotlib import colormaps
#
#     uniq = sorted(set(comm_ids.tolist()))
#     remap = {c: i for i, c in enumerate(uniq)}
#     comm_idx = np.array([remap[c] for c in comm_ids], dtype=int)
#
#     cmap = colormaps.get_cmap("tab20").resampled(max(1, len(uniq)))
#     node_colors = [cmap(i) for i in comm_idx]
#
#     # edge widths
#     if weighted and G.number_of_edges() > 0:
#         abs_w = np.array([G[u][v].get("abs_weight", 1.0) for u, v in G.edges()])
#         if abs_w.max() > abs_w.min():
#             widths = 0.5 + 5.0 * (abs_w - abs_w.min()) / (abs_w.max() - abs_w.min())
#         else:
#             widths = np.ones_like(abs_w) * 2.0
#     else:
#         widths = 0.8
#
#     plt.figure(figsize=(18, 14))
#     nx.draw_networkx_nodes(G, pos, node_size=45, node_color=node_colors, alpha=0.90)
#     nx.draw_networkx_edges(G, pos, width=widths, alpha=0.25)
#
#     if a3_genes:
#         labels = {n: n for n in nodes if n in a3_genes}
#         nx.draw_networkx_labels(G, pos, labels=labels, font_size=10)
#
#     plt.title(title)
#     plt.axis("off")
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=dpi)
#     plt.close()
#
#
# import numpy as np
# import networkx as nx
#
# def community_layout_from_partition(
#     G: nx.Graph,
#     partition: dict,
#     seed: int = 42,
#     scale: float = 10.0,
#     inner_scale: float = 1.5
# ):
#     """
#     Community-aware layout:
#     - partition: dict[node] -> community_id
#     - Places communities as separate blobs (like your example figures)
#     - Then places nodes within each blob
#     Returns pos dict usable for TOP/BOTTOM/DIFF with fixed positions.
#     """
#
#     rng = np.random.default_rng(seed)
#
#     # 1) Build community -> nodes
#     comm_to_nodes = {}
#     for n, c in partition.items():
#         comm_to_nodes.setdefault(c, []).append(n)
#
#     communities = sorted(comm_to_nodes.keys())
#
#     # 2) Build "community graph" (super-nodes)
#     CG = nx.Graph()
#     CG.add_nodes_from(communities)
#
#     # Add weighted edges between communities based on number of inter-community edges
#     for u, v in G.edges():
#         cu, cv = partition.get(u), partition.get(v)
#         if cu is None or cv is None or cu == cv:
#             continue
#         if CG.has_edge(cu, cv):
#             CG[cu][cv]["weight"] += 1.0
#         else:
#             CG.add_edge(cu, cv, weight=1.0)
#
#     # 3) Layout community graph (global positions for each module)
#     if CG.number_of_edges() > 0 and CG.number_of_nodes() > 1:
#         pos_comm = nx.spring_layout(CG, seed=seed, weight="weight")
#     else:
#         # fallback: spread communities in a grid-ish way
#         pos_comm = {}
#         k = int(np.ceil(np.sqrt(len(communities))))
#         for i, c in enumerate(communities):
#             pos_comm[c] = np.array([i % k, i // k], dtype=float)
#
#     # normalize & scale the community centers
#     comm_centers = np.array(list(pos_comm.values()), dtype=float)
#     comm_centers -= comm_centers.mean(axis=0, keepdims=True)
#     max_abs = np.max(np.abs(comm_centers)) if comm_centers.size else 1.0
#     if max_abs == 0:
#         max_abs = 1.0
#     comm_centers = (comm_centers / max_abs) * scale
#     for i, c in enumerate(communities):
#         pos_comm[c] = comm_centers[i]
#
#     # 4) Layout nodes inside each community, then translate to the community center
#     pos = {}
#     for c in communities:
#         nodes = comm_to_nodes[c]
#         sub = G.subgraph(nodes).copy()
#
#         if sub.number_of_nodes() == 1:
#             pos[nodes[0]] = pos_comm[c] + rng.normal(0, 0.01, size=2)
#             continue
#
#         # local layout inside community
#         if sub.number_of_edges() > 0:
#             local = nx.spring_layout(sub, seed=seed, weight="abs_weight")
#         else:
#             # no edges: just scatter locally
#             local = {n: rng.normal(0, 0.2, size=2) for n in sub.nodes()}
#
#         # normalize local coords then scale
#         local_coords = np.array(list(local.values()), dtype=float)
#         local_coords -= local_coords.mean(axis=0, keepdims=True)
#         local_max = np.max(np.abs(local_coords)) if local_coords.size else 1.0
#         if local_max == 0:
#             local_max = 1.0
#         local_coords = (local_coords / local_max) * inner_scale
#
#         for i, n in enumerate(sub.nodes()):
#             pos[n] = pos_comm[c] + local_coords[i]
#
#     return pos
#
#
#
#
#
# # -----------------------------
# # NEW helper functions (NEW NAMES only)
# # -----------------------------
# def nx9_largest_component(G: nx.Graph) -> nx.Graph:
#     """Return largest connected component; if no edges, return copy."""
#     if G.number_of_edges() == 0:
#         return G.copy()
#     comp = max(nx.connected_components(G), key=len)
#     return G.subgraph(comp).copy()
#
# def nx9_remove_degree0(G: nx.Graph) -> nx.Graph:
#     """Remove nodes with degree 0."""
#     H = G.copy()
#     H.remove_nodes_from([n for n in H.nodes() if H.degree(n) == 0])
#     return H
#
# def nx9_merge_communities_to_topK(comm_map: dict, k_keep: int = 5, min_size: int = 40):
#     """
#     Merge small communities into an 'Other' bin to avoid thousands of tiny groups.
#     - keep top k_keep communities by size (and size >= min_size)
#     - everything else -> Other community id (k_keep)
#     Returns: merged_map, raw_sizes, merged_sizes
#     """
#     # community sizes
#     sizes = {}
#     for n, cid in comm_map.items():
#         sizes[cid] = sizes.get(cid, 0) + 1
#
#     # keep only the biggest k_keep
#     top = sorted(sizes.keys(), key=lambda c: sizes[c], reverse=True)[:k_keep]
#     top = set(top)
#
#     other_id = k_keep
#     merged = {}
#     for n, cid in comm_map.items():
#         if (cid in top) and (sizes[cid] >= min_size):
#             merged[n] = cid
#         else:
#             merged[n] = other_id
#
#     merged_sizes = {}
#     for n, cid in merged.items():
#         merged_sizes[cid] = merged_sizes.get(cid, 0) + 1
#
#     return merged, sizes, merged_sizes
#
# def nx9_safe_fill_pos(pos: dict, G: nx.Graph, seed: int = 42):
#     """Ensure pos has all nodes in G (add missing at origin)."""
#     for n in G.nodes():
#         if n not in pos:
#             pos[n] = np.array([0.0, 0.0])
#     return pos
#
#
#
#
# import numpy as np
# import networkx as nx
#
# def community_aware_layout_from_diff_v2(
#     G, partition, seed=42,
#     inter_scale=4.0,     # bigger => communities farther apart
#     intra_scale=1.0,     # bigger => nodes within community more spread
#     k_local=None,        # spring_layout k inside each community
#     k_global=None,       # spring_layout k for community graph
#     iters_local=200,
#     iters_global=200,
#     edge_weight_attr="abs_weight",   # fallback to "weight" if missing
# ):
#     """
#     Community-aware layout computed from a DIFF graph G and a partition dict {node -> community_id}.
#     Returns pos: dict {node -> np.array([x,y])}
#
#     - Inter-community structure from a super-graph (communities as nodes)
#     - Intra-community structure from each community induced subgraph
#     - Final: pos[node] = inter_scale * center[comm] + intra_scale * local[node]
#     """
#
#     rng = np.random.default_rng(seed)
#
#     nodes = list(G.nodes())
#     if len(nodes) == 0:
#         return {}
#
#     # If partition is missing or only 1 community, fall back to spring layout on G
#     if not partition:
#         return nx.spring_layout(G, seed=seed, weight="weight", k=k_global, iterations=iters_global)
#
#     comm_ids = {partition.get(n, -1) for n in nodes}
#     if len(comm_ids) <= 1:
#         return nx.spring_layout(G, seed=seed, weight="weight", k=k_global, iterations=iters_global)
#
#     # Group nodes by community
#     comm_to_nodes = {}
#     for n in nodes:
#         c = partition.get(n, -1)
#         comm_to_nodes.setdefault(c, []).append(n)
#
#     # -------------------------
#     # Build community graph CG
#     # -------------------------
#     CG = nx.Graph()
#     for c in comm_to_nodes.keys():
#         CG.add_node(c)
#
#     # Edge weight between communities = sum of abs(edge_weight) across cross-community edges
#     for u, v, d in G.edges(data=True):
#         cu = partition.get(u, -1)
#         cv = partition.get(v, -1)
#         if cu == cv:
#             continue
#
#         # Choose a weight to summarize cross-community connection strength
#         w = d.get(edge_weight_attr, None)
#         if w is None:
#             w = d.get("weight", 1.0)
#         try:
#             w = float(w)
#         except Exception:
#             w = 1.0
#         w = abs(w)
#
#         if CG.has_edge(cu, cv):
#             CG[cu][cv]["weight"] += w
#         else:
#             CG.add_edge(cu, cv, weight=w)
#
#     # If CG has no edges (communities disconnected), place comms on random points
#     if CG.number_of_edges() == 0:
#         pos_comm = {c: rng.normal(0, 1.0, size=2) for c in CG.nodes()}
#     else:
#         pos_comm = nx.spring_layout(
#             CG, seed=seed, weight="weight",
#             k=k_global, iterations=iters_global
#         )
#
#     # Convert to numpy arrays and apply inter_scale
#     pos_comm = {c: inter_scale * np.asarray(p, dtype=float) for c, p in pos_comm.items()}
#
#     # -------------------------
#     # Local layouts per community
#     # -------------------------
#     pos = {}
#
#     for c, cnodes in comm_to_nodes.items():
#         sub = G.subgraph(cnodes).copy()
#
#         if sub.number_of_nodes() == 1:
#             # single node community
#             only = cnodes[0]
#             local = np.zeros(2, dtype=float)
#             pos[only] = pos_comm.get(c, np.zeros(2)) + intra_scale * local
#             continue
#
#         if sub.number_of_edges() > 0:
#             # Use local spring layout within the module
#             # Use abs_weight if present; otherwise weight
#             # (spring_layout uses the attr name you pass as `weight`)
#             weight_attr = edge_weight_attr if any(edge_weight_attr in dd for _,_,dd in sub.edges(data=True)) else "weight"
#             pos_sub = nx.spring_layout(
#                 sub, seed=seed, weight=weight_attr,
#                 k=k_local, iterations=iters_local
#             )
#             # numpy
#             local_coords = np.array([pos_sub[n] for n in cnodes], dtype=float)
#         else:
#             # no internal edges => tiny jitter cloud
#             local_coords = rng.normal(0, 0.05, size=(len(cnodes), 2))
#
#         # Center local coords at (0,0)
#         local_coords = local_coords - local_coords.mean(axis=0, keepdims=True)
#
#         # Apply intra_scale
#         local_coords = intra_scale * local_coords
#
#         center = pos_comm.get(c, np.zeros(2, dtype=float))
#
#         for i, n in enumerate(cnodes):
#             pos[n] = center + local_coords[i]
#
#     return pos
#
#
#
#
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import colormaps
# #
# # def plot_heatmap_with_community_bar(
# #     mat_df, gene_order, comm_map, out_path,
# #     title="", vmax=None
# # ):
# #     """
# #     mat_df: square DataFrame (genes x genes), signed diff allowed
# #     gene_order: list of genes in the desired order
# #     comm_map: dict gene -> community label (string/int)
# #     Saves heatmap + top colorbar aligned with gene order.
# #     """
# #     genes = [g for g in gene_order if g in mat_df.index]
# #     M = mat_df.loc[genes, genes].values
# #
# #     # Map communities to 0..K-1 for colors
# #     labels = [comm_map.get(g, "Other") for g in genes]
# #     uniq = list(dict.fromkeys(labels))  # preserve order
# #     lab_to_idx = {lab: i for i, lab in enumerate(uniq)}
# #     idx = np.array([lab_to_idx[x] for x in labels], dtype=int)
# #
# #     import matplotlib.pyplot as plt
# #
# #     base_cmap = plt.get_cmap("tab20")
# #     cmap_cat = base_cmap
# #     colors = [cmap_cat(i % cmap_cat.N) for i in range(len(uniq))]
# #
# #
# #     # Figure layout: top bar + heatmap
# #     fig = plt.figure(figsize=(10, 9))
# #     gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[0.25, 8], hspace=0.05)
# #
# #     ax_bar = fig.add_subplot(gs[0, 0])
# #     ax_hm = fig.add_subplot(gs[1, 0])
# #
# #
# #     # Top community bar
# #     colors = np.asarray(colors)
# #     print(colors.shape)  # should be (num_genes, 4) or (num_genes, 3)
# #
# #     ax_bar.imshow(colors[np.newaxis, :, :], aspect="auto")
# #     ax_bar.set_xticks([])
# #     ax_bar.set_yticks([])
# #     ax_bar.set_title(title, fontsize=12)
# #
# #     # Heatmap
# #     if vmax is None:
# #         vmax = np.nanpercentile(np.abs(M), 99)  # robust scaling
# #     im = ax_hm.imshow(M, aspect="auto", vmin=-vmax, vmax=vmax)
# #     ax_hm.set_xticks([])
# #     ax_hm.set_yticks([])
# #
# #     cbar = fig.colorbar(im, ax=ax_hm, fraction=0.046, pad=0.02)
# #     cbar.set_label("DIFF correlation (top - bottom)")
# #
# #     fig.tight_layout()
# #     fig.savefig(out_path, dpi=300)
# #     plt.close(fig)
#
#
# import numpy as np
# import matplotlib.pyplot as plt
#
# def plot_heatmap_with_community_bar(mat_df, gene_order, comm_map, out_path, title="",
#                                    draw_separators=True, annotate_counts=True):
#     # reorder matrix
#     genes = [g for g in gene_order if g in mat_df.index]
#     M = mat_df.loc[genes, genes].values
#
#     # build community vector aligned to genes
#     comms = [comm_map.get(g, "Other") for g in genes]
#
#     # find community blocks (consecutive segments)
#     boundaries = []
#     blocks = []
#     start = 0
#     for i in range(1, len(comms) + 1):
#         if i == len(comms) or comms[i] != comms[i-1]:
#             end = i
#             blocks.append((comms[i-1], start, end))
#             boundaries.append(end)
#             start = i
#
#     # heatmap
#     fig = plt.figure(figsize=(11, 10))
#     gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[20, 1], hspace=0.05)
#
#     ax = fig.add_subplot(gs[0])
#     vmax = np.nanmax(np.abs(M)) if M.size else 1.0
#     if vmax == 0: vmax = 1.0
#
#     im = ax.imshow(M, aspect="auto", interpolation="nearest", vmin=-vmax, vmax=vmax)
#     ax.set_title(title)
#     ax.set_xticks([])
#     ax.set_yticks([])
#
#     # separators on heatmap
#     if draw_separators:
#         for _, s, e in blocks[:-1]:
#             ax.axhline(e - 0.5, linewidth=1.2)
#             ax.axvline(e - 0.5, linewidth=1.2)
#
#     fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
#
#     # community bar
#     axb = fig.add_subplot(gs[1])
#     # encode communities into integers
#     uniq = list(dict.fromkeys(comms))
#     cmap_idx = {c: i for i, c in enumerate(uniq)}
#     bar = np.array([cmap_idx[c] for c in comms])[None, :]
#     axb.imshow(bar, aspect="auto", interpolation="nearest")
#     axb.set_xticks([])
#     axb.set_yticks([])
#
#     # separators on bar + annotate counts
#     if draw_separators:
#         for _, s, e in blocks[:-1]:
#             axb.axvline(e - 0.5, linewidth=1.2)
#
#     if annotate_counts:
#         # put labels centered in each block
#         for c, s, e in blocks:
#             mid = (s + e) / 2
#             count = e - s
#             axb.text(mid, 0, f"{c} (n={count})", ha="center", va="center", fontsize=8, rotation=0)
#
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=300)
#     plt.close()
#
#
#
#
# import matplotlib.pyplot as plt
# import seaborn as sns
# import numpy as np
#
#
# def plot_sbs_vs_a3_selection(
#         df,
#         a3_thr,
#         sbs_low_thr,
#         sbs_high_thr,
#         group1_idx,
#         group2_idx,
#         out_path,
#         title=""
# ):
#     """
#     Visualize SBS2 vs A3_score with selection cutoffs.
#
#     Top panel:
#         SBS2 (Y) vs A3_score (X)
#         - vertical line: A3_score cutoff
#         - horizontal lines: SBS2 20% and 80% (within High-A3)
#         - colored selected groups
#
#     Bottom panel:
#         Distribution of A3_score
#         - histogram + KDE
#         - mean + median lines
#         - cutoff line
#     """
#
#     fig = plt.figure(figsize=(10, 10))
#     gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])
#
#     ax_scatter = fig.add_subplot(gs[0])
#     ax_hist = fig.add_subplot(gs[1])
#
#     # ===============================
#     # TOP PANEL — Scatter
#     # ===============================
#     ax_scatter.scatter(
#         df["A3_score"], df["SBS2"],
#         s=20, alpha=0.3, color="gray", label="All samples"
#     )
#
#     # highlight groups
#     ax_scatter.scatter(
#         df.loc[group1_idx, "A3_score"],
#         df.loc[group1_idx, "SBS2"],
#         s=40, color="red", label="Group1 (High SBS2 in High-A3)"
#     )
#
#     ax_scatter.scatter(
#         df.loc[group2_idx, "A3_score"],
#         df.loc[group2_idx, "SBS2"],
#         s=40, color="blue", label="Group2 (Low SBS2 in High-A3)"
#     )
#
#     # Cutoff lines
#     ax_scatter.axvline(a3_thr, color="black", linestyle="--", label="A3 median cutoff")
#     ax_scatter.axhline(sbs_high_thr, color="red", linestyle="--", label="SBS2 80% (High-A3)")
#     ax_scatter.axhline(sbs_low_thr, color="blue", linestyle="--", label="SBS2 20% (High-A3)")
#
#     ax_scatter.set_xlabel("A3_score")
#     ax_scatter.set_ylabel("SBS2")
#     ax_scatter.set_title(title)
#     ax_scatter.legend(loc="best")
#     ax_scatter.grid(alpha=0.2)
#
#     # ===============================
#     # BOTTOM PANEL — A3 Distribution
#     # ===============================
#     sns.histplot(df["A3_score"], bins=40, kde=True, ax=ax_hist, color="gray")
#
#     mean_val = df["A3_score"].mean()
#     median_val = df["A3_score"].median()
#
#     ax_hist.axvline(mean_val, color="green", linestyle="-", label=f"Mean = {mean_val:.3f}")
#     ax_hist.axvline(median_val, color="black", linestyle="--", label=f"Median = {median_val:.3f}")
#     ax_hist.axvline(a3_thr, color="red", linestyle=":", label=f"Cutoff = {a3_thr:.3f}")
#
#     ax_hist.set_xlabel("A3_score")
#     ax_hist.set_ylabel("Count")
#     ax_hist.legend(loc="best")
#     ax_hist.grid(alpha=0.2)
#
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=300)
#     plt.close()
#
#
#
#
#
# import numpy as np
# import pandas as pd
# import os
#
# def nx9_order_genes_by_community(comm_map, genes=None, other_label="Other"):
#     """
#     Return a community-ordered gene list.
#     - Communities ordered by size descending.
#     - Genes within each community sorted alphabetically (stable).
#     - "Other" community (if present) placed last.
#     """
#     if genes is None:
#         genes = list(comm_map.keys())
#     genes = [g for g in genes if g in comm_map]
#
#     # community -> genes
#     comm_to_genes = {}
#     for g in genes:
#         c = comm_map.get(g, other_label)
#         comm_to_genes.setdefault(c, []).append(g)
#
#     # order communities by size (desc), but put "Other" last
#     comm_items = list(comm_to_genes.items())
#     comm_items.sort(key=lambda kv: len(kv[1]), reverse=True)
#
#     if other_label in comm_to_genes:
#         comm_items = [kv for kv in comm_items if kv[0] != other_label] + [(other_label, comm_to_genes[other_label])]
#
#     ordered = []
#     for c, glist in comm_items:
#         ordered.extend(sorted(glist))  # stable inside module
#
#     return ordered
#
#
# def nx9_restrict_and_order_matrix(mat_df, ordered_genes):
#     """Safely subset matrix to ordered_genes (intersection only) and return ordered matrix + final gene list."""
#     g = [x for x in ordered_genes if x in mat_df.index]
#     g = [x for x in g if x in mat_df.columns]
#     return mat_df.loc[g, g], g
#
#
# def nx9_plot_simple_heatmap(mat_df, out_path, title, vmax=None):
#     arr = mat_df.values
#     # optional symmetric scale for DIFF
#     if vmax is None:
#         vmax = np.nanmax(np.abs(arr)) if arr.size else 1.0
#         if vmax == 0:
#             vmax = 1.0
#
#     plt.figure(figsize=(10, 9))
#     plt.imshow(arr, aspect="auto", interpolation="nearest", vmin=-vmax, vmax=vmax)
#     plt.colorbar(label="value")
#     plt.title(title)
#     plt.xlabel("genes (community order)")
#     plt.ylabel("genes (community order)")
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=300)
#     plt.close()
#
#
#
#
#
# import numpy as np
# import networkx as nx
#
# def nx9_partition_to_communities(comm_map):
#     """Convert {node: community_id} to list-of-sets (community sets)."""
#     comm_to_nodes = {}
#     for n, c in comm_map.items():
#         comm_to_nodes.setdefault(c, set()).add(n)
#     return list(comm_to_nodes.values())
#
# def nx9_modularity_nx(G, comm_map, weight="weight"):
#     """Compute modularity using NetworkX quality.modularity."""
#     from networkx.algorithms.community.quality import modularity
#     comms = nx9_partition_to_communities(comm_map)
#     if len(comms) == 0:
#         return np.nan
#     return float(modularity(G, comms, weight=weight))
#
# def detect_communities_on_graph_v2(
#     G,
#     method="leiden",
#     seed=42,
#     resolution=1.0,
#     weight="weight",
#     return_modularity=False
# ):
#     """
#     NEW version (keeps old function untouched).
#     Supports resolution and optional modularity output.
#
#     Returns:
#       - comm_map: dict {node: community_id}
#       - modularity (optional): float
#     """
#
#     # ---------
#     # Safety
#     # ---------
#     if G is None or G.number_of_nodes() == 0:
#         if return_modularity:
#             return {}, np.nan
#         return {}
#
#     method = (method or "").lower()
#
#     # ==========================================================
#     # 1) LEIDEN (igraph + leidenalg)
#     # ==========================================================
#     if method == "leiden":
#         try:
#             import igraph as ig
#             import leidenalg as la
#
#             nodes = list(G.nodes())
#             idx = {n: i for i, n in enumerate(nodes)}
#
#             # build igraph edges + weights
#             edges = []
#             weights = []
#             for u, v, d in G.edges(data=True):
#                 edges.append((idx[u], idx[v]))
#                 w = float(d.get(weight, 1.0))
#                 weights.append(w)
#
#             igG = ig.Graph(n=len(nodes), edges=edges, directed=False)
#             if len(weights) == len(edges) and len(edges) > 0:
#                 igG.es["weight"] = weights
#                 part = la.find_partition(
#                     igG,
#                     la.RBConfigurationVertexPartition,
#                     weights="weight",
#                     resolution_parameter=float(resolution),
#                     seed=int(seed)
#                 )
#             else:
#                 part = la.find_partition(
#                     igG,
#                     la.RBConfigurationVertexPartition,
#                     resolution_parameter=float(resolution),
#                     seed=int(seed)
#                 )
#
#             # membership list -> comm_map
#             membership = part.membership
#             comm_map = {nodes[i]: int(membership[i]) for i in range(len(nodes))}
#
#             if return_modularity:
#                 # leidenalg partition has quality() which is modularity-like for RBConfiguration
#                 # but to keep consistent with NetworkX graphs, we compute NX modularity:
#                 mod = nx9_modularity_nx(G, comm_map, weight=weight)
#                 return comm_map, mod
#             return comm_map
#
#         except Exception:
#             # If leiden libs not available, fall through to Louvain
#             method = "louvain"
#
#     # ==========================================================
#     # 2) LOUVAIN (python-louvain: community)
#     # ==========================================================
#     if method == "louvain":
#         try:
#             import community as community_louvain  # python-louvain
#             # best_partition supports resolution
#             comm_map = community_louvain.best_partition(
#                 G, weight=weight, resolution=float(resolution), random_state=int(seed)
#             )
#             comm_map = {k: int(v) for k, v in comm_map.items()}
#
#             if return_modularity:
#                 mod = nx9_modularity_nx(G, comm_map, weight=weight)
#                 return comm_map, mod
#             return comm_map
#         except Exception as e:
#             raise RuntimeError(
#                 "Could not run Louvain. Install python-louvain: pip install python-louvain"
#             ) from e
#
#     raise ValueError(f"Unknown method='{method}'. Use 'leiden' or 'louvain'.")
#
#
#
#
#
# def community_aware_layout_from_diff_v3(
#     G, comm_map, seed=42,
#     inter_scale=1.8,     # SMALLER => communities closer (try 1.0 ~ 2.0)
#     intra_scale=2.5,     # BIGGER  => community internal spread larger (try 1.5 ~ 4.0)
#     k_local=None,
#     k_global=None,
#     iters_local=200,
#     iters_global=200,
#     weight_attr="abs_weight"
# ):
#     import numpy as np
#     import networkx as nx
#
#     rng = np.random.default_rng(seed)
#
#     nodes = list(G.nodes())
#     if len(nodes) == 0:
#         return {}
#
#     # fallback if comm_map missing / trivial
#     if (not comm_map) or len(set(comm_map.values())) <= 1:
#         return nx.spring_layout(G if G.number_of_edges() > 0 else nx.Graph(G), seed=seed)
#
#     # group nodes by community
#     comm_to_nodes = {}
#     for n in nodes:
#         c = comm_map.get(n, "Other")
#         comm_to_nodes.setdefault(c, []).append(n)
#
#     # build community graph (super-nodes)
#     CG = nx.Graph()
#     for c in comm_to_nodes:
#         CG.add_node(c)
#
#     for u, v, d in G.edges(data=True):
#         cu = comm_map.get(u, "Other")
#         cv = comm_map.get(v, "Other")
#         if cu == cv:
#             continue
#         w = float(d.get(weight_attr, abs(d.get("weight", 1.0))))
#         if CG.has_edge(cu, cv):
#             CG[cu][cv]["weight"] += w
#         else:
#             CG.add_edge(cu, cv, weight=w)
#
#     # community centroid layout (global)
#     pos_comm = nx.spring_layout(
#         CG,
#         seed=seed,
#         weight="weight",
#         k=k_global,
#         iterations=iters_global
#     )
#
#     # scale centroids to control inter-community distance
#     for c in pos_comm:
#         pos_comm[c] = inter_scale * np.array(pos_comm[c], dtype=float)
#
#     # layout within each community
#     pos = {}
#     for c, cnodes in comm_to_nodes.items():
#         sub = G.subgraph(cnodes).copy()
#
#         if sub.number_of_nodes() == 1:
#             pos_sub = {cnodes[0]: np.array([0.0, 0.0])}
#         elif sub.number_of_edges() > 0:
#             pos_sub = nx.spring_layout(
#                 sub,
#                 seed=seed,
#                 weight=weight_attr if any(weight_attr in sub[u][v] for u,v in sub.edges()) else "weight",
#                 k=k_local,
#                 iterations=iters_local
#             )
#         else:
#             pos_sub = {n: rng.normal(0, 0.05, size=2) for n in cnodes}
#
#         # center & scale inside community
#         arr = np.array([pos_sub[n] for n in cnodes], dtype=float)
#         arr = arr - arr.mean(axis=0, keepdims=True)
#         arr = intra_scale * arr  # bigger => more spread
#         pos_sub = {n: arr[i] for i, n in enumerate(cnodes)}
#
#         centroid = np.array(pos_comm.get(c, [0.0, 0.0]), dtype=float)
#         for n in cnodes:
#             pos[n] = pos_sub[n] + centroid
#
#     return pos
#
#
#
#
# def plot_zoomed_communities(
#     G, pos, comm_map, out_dir, title_prefix,
#     a3_genes=None, pad=0.25,
#     min_nodes=10,
#     weighted=True,
#     label_apobec=True
# ):
#     import os
#     import numpy as np
#     import matplotlib.pyplot as plt
#     import networkx as nx
#
#     os.makedirs(out_dir, exist_ok=True)
#     a3_genes = set(a3_genes or [])
#
#     # group nodes by community
#     comm_to_nodes = {}
#     for n in G.nodes():
#         c = comm_map.get(n, "Other")
#         comm_to_nodes.setdefault(c, []).append(n)
#
#     # stable order: largest communities first
#     comm_order = sorted(comm_to_nodes.keys(), key=lambda c: len(comm_to_nodes[c]), reverse=True)
#
#     for c in comm_order:
#         nodes_c = comm_to_nodes[c]
#         if len(nodes_c) < min_nodes:
#             continue
#
#         sub = G.subgraph(nodes_c).copy()
#
#         # ensure APOBEC nodes are included if they exist in the full graph AND this community contains them
#         a3_in = [g for g in nodes_c if g in a3_genes]
#
#         # get bounding box from global pos for those nodes
#         pts = np.array([pos[n] for n in nodes_c if n in pos], dtype=float)
#         if pts.size == 0:
#             continue
#
#         xmin, ymin = pts.min(axis=0)
#         xmax, ymax = pts.max(axis=0)
#
#         dx = (xmax - xmin) if xmax > xmin else 1.0
#         dy = (ymax - ymin) if ymax > ymin else 1.0
#
#         xmin -= pad * dx; xmax += pad * dx
#         ymin -= pad * dy; ymax += pad * dy
#
#         plt.figure(figsize=(8, 8))
#         ax = plt.gca()
#
#         # draw edges
#         nx.draw_networkx_edges(sub, pos, ax=ax, alpha=0.25)
#
#         # draw nodes
#         nx.draw_networkx_nodes(sub, pos, ax=ax, node_size=30, alpha=0.9)
#
#         # highlight APOBEC
#         if a3_in:
#             nx.draw_networkx_nodes(
#                 sub, pos, ax=ax,
#                 nodelist=a3_in, node_size=140, alpha=0.95
#             )
#             if label_apobec:
#                 nx.draw_networkx_labels(
#                     sub, pos,
#                     labels={g: g for g in a3_in},
#                     font_size=10
#                 )
#
#         ax.set_title(f"{title_prefix} | Community={c} | n={len(nodes_c)}")
#         ax.set_xlim(xmin, xmax)
#         ax.set_ylim(ymin, ymax)
#         ax.axis("off")
#
#         out_path = os.path.join(out_dir, f"{title_prefix}_COMM_{str(c)}_zoom.png")
#         plt.tight_layout()
#         plt.savefig(out_path, dpi=300)
#         plt.close()
#
#
#
#
# def nx9_resolution_sweep_stability(
#     G,
#     out_dir,
#     method="leiden",
#     resolutions=None,
#     n_runs=10,
#     base_seed=42,
#     weight="abs_weight",
#     verbose=True
# ):
#     """
#     Runs community detection multiple times per resolution.
#     Saves:
#       - sweep_summary.csv
#       - partitions per (resolution, seed)
#       - plots: ncomms vs resolution, modularity vs resolution, ARI vs resolution, NMI vs resolution
#       - ONE-PLOT option: ncomms only (with std shading)
#       - ONE-PLOT enhanced: ncomms colored by ARI (optional)
#     """
#     import os, json, itertools
#     import numpy as np
#     import pandas as pd
#     import matplotlib.pyplot as plt
#     from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
#
#     os.makedirs(out_dir, exist_ok=True)
#     part_dir = os.path.join(out_dir, "partitions")
#     plot_dir = os.path.join(out_dir, "plots")
#     os.makedirs(part_dir, exist_ok=True)
#     os.makedirs(plot_dir, exist_ok=True)
#
#     if resolutions is None:
#         resolutions = np.round(np.linspace(0.2, 2.0, 10), 2)
#
#     nodes_eval = list(G.nodes())
#
#     def partition_to_labels(cm):
#         comm_ids = [cm.get(n, -1) for n in nodes_eval]
#         uniq = {c: i for i, c in enumerate(sorted(set(comm_ids)))}
#         return np.array([uniq[c] for c in comm_ids], dtype=int)
#
#     rows = []
#
#     for r in resolutions:
#         labels_list = []
#         modularities = []
#         ncomms_list = []
#
#         for i in range(n_runs):
#             seed = base_seed + i
#             cm, mod = detect_communities_on_graph_v2(
#                 G, method=method, seed=seed, resolution=float(r),
#                 weight=weight, return_modularity=True
#             )
#
#             # save partition
#             part_path = os.path.join(part_dir, f"partition__res{r:.2f}__seed{seed}.json")
#             with open(part_path, "w") as f:
#                 json.dump(cm, f)
#
#             labels_list.append(partition_to_labels(cm))
#             ncomms_list.append(len(set(cm.values())))
#             modularities.append(mod if mod is not None else np.nan)
#
#         # pairwise stability
#         aris, nmis = [], []
#         for a, b in itertools.combinations(range(n_runs), 2):
#             aris.append(adjusted_rand_score(labels_list[a], labels_list[b]))
#             nmis.append(normalized_mutual_info_score(labels_list[a], labels_list[b]))
#
#         rows.append({
#             "resolution": float(r),
#             "n_runs": int(n_runs),
#             "ncomms_mean": float(np.mean(ncomms_list)),
#             "ncomms_std": float(np.std(ncomms_list)),
#             "modularity_mean": float(np.nanmean(modularities)),
#             "modularity_std": float(np.nanstd(modularities)),
#             "ARI_mean": float(np.mean(aris)),
#             "ARI_std": float(np.std(aris)),
#             "NMI_mean": float(np.mean(nmis)),
#             "NMI_std": float(np.std(nmis)),
#         })
#
#     sweep_df = pd.DataFrame(rows).sort_values("resolution")
#     sweep_csv = os.path.join(out_dir, "sweep_summary.csv")
#     sweep_df.to_csv(sweep_csv, index=False)
#
#     # ---- plot 1: ncomms only (what you asked)
#     plt.figure(figsize=(8, 5))
#     x = sweep_df["resolution"].values
#     y = sweep_df["ncomms_mean"].values
#     yerr = sweep_df["ncomms_std"].values
#     plt.plot(x, y, marker="o")
#     plt.fill_between(x, y - yerr, y + yerr, alpha=0.2)
#     plt.xlabel("resolution")
#     plt.ylabel("#communities (mean ± std)")
#     plt.title("Resolution sweep: #communities")
#     plt.tight_layout()
#     plt.savefig(os.path.join(plot_dir, "ncomms_only.png"), dpi=300)
#     plt.close()
#
#     # ---- plot 2: one-plot enhanced (ncomms colored by ARI)
#     plt.figure(figsize=(8, 5))
#     sc = plt.scatter(
#         sweep_df["resolution"], sweep_df["ncomms_mean"],
#         c=sweep_df["ARI_mean"], s=60
#     )
#     plt.colorbar(sc, label="ARI_mean")
#     plt.xlabel("resolution")
#     plt.ylabel("#communities (mean)")
#     plt.title("Resolution sweep: #communities colored by stability (ARI)")
#     plt.tight_layout()
#     plt.savefig(os.path.join(plot_dir, "ncomms_colored_by_ARI.png"), dpi=300)
#     plt.close()
#
#     # (Optional) keep your other plots too:
#     for col, fname, ylabel in [
#         ("ARI_mean", "ARI_vs_resolution.png", "ARI (mean)"),
#         ("NMI_mean", "NMI_vs_resolution.png", "NMI (mean)"),
#         ("modularity_mean", "modularity_vs_resolution.png", "Modularity (mean)"),
#     ]:
#         plt.figure(figsize=(8, 5))
#         plt.plot(sweep_df["resolution"], sweep_df[col], marker="o")
#         plt.xlabel("resolution")
#         plt.ylabel(ylabel)
#         plt.title(f"Resolution sweep: {ylabel}")
#         plt.tight_layout()
#         plt.savefig(os.path.join(plot_dir, fname), dpi=300)
#         plt.close()
#
#     return sweep_df, sweep_csv, part_dir, plot_dir
#
#
#
#
# import numpy as np
# import matplotlib.pyplot as plt
#
# def plot_corr_violin(corr_top, corr_bottom, corr_diff, out_path_prefix):
#     def upper_vals(mat):
#         a = mat.values
#         iu = np.triu_indices_from(a, k=1)
#         return a[iu]
#
#     top_vals = upper_vals(corr_top)
#     bot_vals = upper_vals(corr_bottom)
#     diff_vals = upper_vals(corr_diff)
#
#     plt.figure(figsize=(8, 4))
#     plt.violinplot([top_vals, bot_vals, diff_vals], showmeans=False, showmedians=True, showextrema=True)
#     plt.xticks([1, 2, 3], ["TOP", "BOTTOM", "DIFF (TOP-BOT)"])
#     plt.ylabel("Spearman correlation")
#     plt.title("Correlation distributions (violin)")
#     plt.tight_layout()
#     plt.savefig(f"{out_path_prefix}_violin.png", dpi=300)
#     plt.close()
#
#
#
#
# import json
# import pandas as pd
# import os
# import matplotlib.pyplot as plt
#
# def save_partition_outputs(G, pos, comm_map, out_dir, tag):
#     os.makedirs(out_dir, exist_ok=True)
#
#     # 1) JSON mapping
#     json_path = os.path.join(out_dir, f"{tag}_comm_map.json")
#     with open(json_path, "w") as f:
#         json.dump(comm_map, f, indent=2)
#
#     # 2) CSV gene list per community
#     comm_to_genes = {}
#     for g, c in comm_map.items():
#         comm_to_genes.setdefault(str(c), []).append(g)
#
#     rows = []
#     for c, genes in sorted(comm_to_genes.items(), key=lambda x: -len(x[1])):
#         rows.append({
#             "community": c,
#             "size": len(genes),
#             "genes": ";".join(sorted(genes))
#         })
#
#     csv_path = os.path.join(out_dir, f"{tag}_community_geneList.csv")
#     pd.DataFrame(rows).to_csv(csv_path, index=False)
#
#     # 3) Plot (colored by community)
#     # If you already have plot_network_by_community_fixedpos, use it:
#     from network_utils import plot_network_by_community_fixedpos
#     plot_path = os.path.join(out_dir, f"{tag}_community_plot.png")
#     plot_network_by_community_fixedpos(
#         G, pos, comm_map,
#         out_path=plot_path,
#         title=tag,
#         a3_genes=["A3A","A3B","A3C","A3D","A3F","A3G","A3H"],
#         weighted=True
#     )
#
#
#
#
#
#
# # =============================================================================
# # Helper 1 — DIFF graph strategy: TOPK edges per node (no hard threshold)
# # =============================================================================
# def build_diff_graph_topk(corr_top: pd.DataFrame,
#                           corr_bottom: pd.DataFrame,
#                           k: int = 10,
#                           use_abs_diff: bool = True) -> nx.Graph:
#     """
#     Build a DIFF graph without a global threshold:
#     For each node i, connect it to the TOP-k nodes by |diff(i,j)|.
#     """
#     genes = list(corr_top.index)
#     diff = (corr_top - corr_bottom).fillna(0.0)
#     if use_abs_diff:
#         score = diff.abs()
#     else:
#         score = diff
#
#     G = nx.Graph()
#     G.add_nodes_from(genes)
#
#     for i, g in enumerate(genes):
#         # take topk neighbors excluding self
#         s = score.loc[g].copy()
#         s[g] = -np.inf
#         top_neighbors = s.nlargest(k).index.tolist()
#
#         for h in top_neighbors:
#             w_signed = float(diff.loc[g, h])
#             w_abs = float(abs(w_signed))
#             if not np.isfinite(w_abs) or w_abs <= 0:
#                 continue
#             if G.has_edge(g, h):
#                 # keep strongest abs
#                 if w_abs > G[g][h].get("abs_weight", 0):
#                     G[g][h]["weight"] = w_signed
#                     G[g][h]["abs_weight"] = w_abs
#             else:
#                 G.add_edge(g, h, weight=w_signed, abs_weight=w_abs)
#
#     return G
#
#
# # =============================================================================
# # Helper 2 — DIFF graph strategy: Soft-threshold (WGCNA-like)
# # =============================================================================
# def build_soft_wgcna_graph(corr: pd.DataFrame,
#                            power: int = 6,
#                            keep_topk_per_node: int = 30) -> nx.Graph:
#     """
#     WGCNA-like adjacency:
#       adjacency(i,j) = |corr(i,j)|^power
#     """
#     genes = list(corr.index)
#     A = corr.abs().fillna(0.0) ** power
#     np.fill_diagonal(A.values, 0.0)
#
#     G = nx.Graph()
#     G.add_nodes_from(genes)
#
#     for g in genes:
#         s = A.loc[g]
#         top_neighbors = s.nlargest(keep_topk_per_node).index.tolist()
#         for h in top_neighbors:
#             w = float(A.loc[g, h])
#             if w <= 0 or not np.isfinite(w):
#                 continue
#             if g == h:
#                 continue
#             if not G.has_edge(g, h):
#                 G.add_edge(g, h, weight=w, abs_weight=w)
#             else:
#                 if w > G[g][h].get("abs_weight", 0):
#                     G[g][h]["weight"] = w
#                     G[g][h]["abs_weight"] = w
#     return G
#
#
# def build_diff_graph_soft_wgcna(corr_top: pd.DataFrame,
#                                 corr_bottom: pd.DataFrame,
#                                 power: int = 6,
#                                 topk_per_node: int = 30) -> nx.Graph:
#     """
#     Create soft adjacency for each group, then DIFF as difference in adjacency.
#     For community detection, we keep abs(diff) as abs_weight.
#     """
#     A_top = corr_top.abs().fillna(0.0) ** power
#     A_bot = corr_bottom.abs().fillna(0.0) ** power
#     D = (A_top - A_bot).fillna(0.0)
#     genes = list(D.index)
#
#     G = nx.Graph()
#     G.add_nodes_from(genes)
#     np.fill_diagonal(D.values, 0.0)
#
#     for g in genes:
#         s = D.loc[g].abs()
#         s[g] = -np.inf
#         neigh = s.nlargest(topk_per_node).index.tolist()
#         for h in neigh:
#             w_signed = float(D.loc[g, h])
#             w_abs = float(abs(w_signed))
#             if w_abs <= 0 or not np.isfinite(w_abs):
#                 continue
#             if not G.has_edge(g, h):
#                 G.add_edge(g, h, weight=w_signed, abs_weight=w_abs)
#             else:
#                 if w_abs > G[g][h].get("abs_weight", 0):
#                     G[g][h]["weight"] = w_signed
#                     G[g][h]["abs_weight"] = w_abs
#     return G
#
#
# # =============================================================================
# # STEP 6 PLOTS — Volcano & Manhattan
# # =============================================================================
# def plot_volcano(stats_df, outpath, title=None, q_thr=0.05):
#     """Volcano plot: log2FC vs -log10(p-value)."""
#     if stats_df.empty:
#         return
#     df = stats_df.copy()
#     df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
#
#     sig_mask = df["qvalue"] <= q_thr
#     fig, ax = plt.subplots(figsize=(7, 5))
#     ax.scatter(df["log2FC_high_vs_low"], df["neglog10p"],
#                c="lightgray", s=10, alpha=0.7, label="NS")
#     ax.scatter(df.loc[sig_mask, "log2FC_high_vs_low"],
#                df.loc[sig_mask, "neglog10p"],
#                c="red", s=12, alpha=0.8, label=f"FDR ≤ {q_thr:g}")
#
#     ax.axhline(-np.log10(0.05), color="black", ls="--", lw=1, alpha=0.6)
#     ax.axvline(0.0, color="black", ls=":", lw=1, alpha=0.6)
#
#     ax.set_xlabel("log2FC (High SBS2 vs Low SBS2)")
#     ax.set_ylabel("-log10(p-value)")
#     if title:
#         ax.set_title(title)
#     ax.legend(frameon=False, fontsize=8)
#     plt.tight_layout()
#     plt.savefig(outpath, dpi=300)
#     plt.close()
#
#
# def plot_manhattan_index(stats_df, outpath, title=None, q_thr=0.05):
#     """
#     Manhattan-like plot using gene index as pseudo-position.
#     """
#     if stats_df.empty:
#         return
#     df = stats_df.copy().reset_index(drop=True)
#     df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
#     df["idx"] = np.arange(len(df))
#
#     sig_mask = df["qvalue"] <= q_thr
#
#     fig, ax = plt.subplots(figsize=(8, 4))
#     ax.scatter(df["idx"], df["neglog10p"],
#                c="gray", s=6, alpha=0.6, label="NS")
#     ax.scatter(df.loc[sig_mask, "idx"],
#                df.loc[sig_mask, "neglog10p"],
#                c="red", s=8, alpha=0.8, label=f"FDR ≤ {q_thr:g}")
#
#     ax.axhline(-np.log10(0.05), color="black", ls="--", lw=1, alpha=0.6)
#
#     ax.set_xlabel("Gene index")
#     ax.set_ylabel("-log10(p-value)")
#     if title:
#         ax.set_title(title)
#     ax.legend(frameon=False, fontsize=8)
#     plt.tight_layout()
#     plt.savefig(outpath, dpi=300)
#     plt.close()
#
#
# # =============================================================================
# # Network saving helpers (STEP 8)
# # =============================================================================
# def remove_isolated_nodes(G: nx.Graph) -> nx.Graph:
#     H = G.copy()
#     isolates = list(nx.isolates(H))
#     H.remove_nodes_from(isolates)
#     return H
#
#
# def save_graph_edgelist(G: nx.Graph, outpath: str, weighted: bool = True):
#     """Save graph as edge list (CSV-like)."""
#     import csv
#     with open(outpath, "w", newline="") as f:
#         writer = csv.writer(f)
#         if weighted:
#             writer.writerow(["source", "target", "weight"])
#             for u, v, d in G.edges(data=True):
#                 w = float(d.get("weight", 0.0))
#                 writer.writerow([u, v, w])
#         else:
#             writer.writerow(["source", "target"])
#             for u, v in G.edges():
#                 writer.writerow([u, v])
#
#
#
#
# def clean_ensg(x: str) -> str:
#     # remove trailing ".digits" version
#     return re.sub(r"\.\d+$", "", str(x))
#
#
#
#
# # ------------------------------------------------------------------
# # Visibility helper
# # ------------------------------------------------------------------
# def banner(title: str, char: str = "=", width: int = 100):
#     line = char * width
#     print("\n" + line)
#     print(title)
#     print(line)
#
#
#
#
# # Create a "short" version of Entity_ID with first 4 parts
# def shorten_entity_id(x):
#     return "-".join(str(x).split("-")[:4])
#
#
#
#
# # ensure abs_weight exists
# def ensure_abs_weight(G):
#     for u, v, d in G.edges(data=True):
#         w = float(d.get("weight", 0.0))
#         d["abs_weight"] = abs(float(d.get("abs_weight", w)))
#
#
#
#
# # biomarker presence in graphs
# def report_gene_presence(tag, G, genes):
#     present = [g for g in genes if g in G.nodes]
#     missing = [g for g in genes if g not in G.nodes]
#     print(f"\n[{tag}] present: {len(present)} | missing: {len(missing)}")
#     if present:
#         print("  present:", present)
#     if missing:
#         print("  missing:", missing)
#
#
#
#
# def save_edges_both_formats(G, out_prefix):
#     save_graph_edgelist(G, out_prefix + "_weighted.csv", weighted=True)
#     save_graph_edgelist(G, out_prefix + "_unweighted.csv", weighted=False)
#
#
#
# def partition_to_labels(cm: dict, nodes: list):
#     """Convert {node:community} -> integer labels aligned to nodes list."""
#     ids = [cm.get(n, -1) for n in nodes]
#     uniq = {c: i for i, c in enumerate(sorted(set(ids)))}
#     return np.array([uniq[c] for c in ids], dtype=int)
#
#
# def save_curve(x, y, title, ylabel, out_path):
#     plt.figure(figsize=(8, 5))
#     plt.plot(x, y, marker="o")
#     plt.xlabel("resolution")
#     plt.ylabel(ylabel)
#     plt.title(title)
#     plt.tight_layout()
#     plt.savefig(out_path, dpi=300)
#     plt.close()
#
#
# def save_gene_lists_from_map(cm_map: dict, out_csv: str):
#     """Save one row per community with gene list."""
#     comm_to_genes = {}
#     for g, c in cm_map.items():
#         comm_to_genes.setdefault(c, []).append(g)
#
#     rows = []
#     for c, genes in comm_to_genes.items():
#         rows.append({"community": c, "size": len(genes), "genes": ";".join(sorted(genes))})
#
#     pd.DataFrame(rows).sort_values("size", ascending=False).to_csv(out_csv, index=False)
#
#
# # A3_ALIASES = ["A3A", "A3B", "A3C", "A3D", "A3G", "A3H"]
# #
# #
# # def nodes_that_exist(G, gene_list):
# #     return [g for g in gene_list if g in G.nodes]
# #
# # def get_a3_nodes_for_graph(G):
# #     # ENSG mode
# #     if any(g in G.nodes for g in A3_GENES):
# #         return [g for g in A3_GENES if g in G.nodes]
# #     # alias mode
# #     if any(a in G.nodes for a in A3_ALIASES):
# #         return [a for a in A3_ALIASES if a in G.nodes]
# #     return []
# #
# # def get_biomarker_nodes_for_graph(G):
# #     # Biomarkers are ENSG in your code. If you ever add alias columns, you can extend here.
# #     return [b for b in BIOMARKERS if b in G.nodes]