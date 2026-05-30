#!/usr/bin/env python3
"""
Generate_Panel_D_v2.py
========================

Simplified Panel D for Figure 4: three-comparison dot plot showing
chain co-expression validation between SBS2-HIGH and CNV-HIGH cells.

Three comparisons only (clear narrative arc):
  1. Activating chain internal   -> stronger in SBS2-HIGH
  2. Inhibiting chain internal   -> gradient (NORM > SBS2 > CNV)
  3. Cross-chain (Act vs Inh)    -> near zero (independent programs)

Caches the DIFF matrix to avoid recomputing correlations on subsequent runs.

Usage:
    conda run -n NETWORK python Generate_Panel_D_v2.py
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")
ADATA_PATH = os.path.join(FIG4_ROOT, "00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(FIG4_ROOT, "01_group_selection/three_group_assignments.tsv")
FIG_DIR = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS")
CACHE_DIR = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS", "CACHE")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}

# Style (inherits Figure 2/4 conventions)
COLOR_ACTIVATING = "#fcad61"   # orange, up in tumor
COLOR_INHIBITING = "#aadda4"   # green, up in normal
COLOR_CROSSCHAIN = "#888888"   # gray, independence
FONT_TITLE = 32
FONT_AXIS  = 30
FONT_TICK  = 28
FONT_LABEL = 26


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# CHAIN IDENTIFICATION (reused from Generate_Figure4_Panels.py)
# =============================================================================

def identify_chain_genes(net_dir, harris_genes):
    """Identify activator and repressor chain genes from one network."""
    part_path = os.path.join(net_dir, "04_communities/SC_best_partition.csv")
    graph_path = os.path.join(net_dir, "04_communities/SC_G_comm.gpickle")
    de_path = os.path.join(net_dir, "02_differential_expression/SC_diffexpr_stats.csv")

    for p in [part_path, graph_path, de_path]:
        if not os.path.exists(p):
            log(f"  [ERROR] Required file missing: {p}")
            raise FileNotFoundError(p)

    part_df = pd.read_csv(part_path)
    gene_to_comm = dict(zip(part_df["gene"], part_df["community"]))

    with open(graph_path, "rb") as f:
        G_comm = pickle.load(f)

    de_df = pd.read_csv(de_path)
    de_log2fc = dict(zip(de_df["gene"], de_df["log2FC"]))

    # Find A3 community
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

    # Harris-anchored chains
    for gene in harris_in_comm:
        if de_log2fc.get(gene, 0) > 0 and gene in G_act and G_act.degree(gene) > 0:
            activator_nodes |= set(nx.node_connected_component(G_act, gene))
        if de_log2fc.get(gene, 0) <= 0 and gene in G_rep and G_rep.degree(gene) > 0:
            repressor_nodes |= set(nx.node_connected_component(G_rep, gene))

    # Boundary gene chains (A3 neighbor -> neighbor's neighbor)
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
                        activator_nodes |= set(nx.node_connected_component(G_act, nb))
                    break
                elif nb_fc <= 0 and nn_fc <= 0 and nn_w < 0:
                    if nb in G_rep and G_rep.degree(nb) > 0:
                        repressor_nodes |= set(nx.node_connected_component(G_rep, nb))
                    break

    activator_nodes -= a3_set
    repressor_nodes -= a3_set
    return activator_nodes, repressor_nodes


# =============================================================================
# CORRELATION COMPUTATION
# =============================================================================

def compute_corr(adata, cells, genes):
    """Compute pairwise Spearman correlation matrix for genes in cells."""
    valid_genes = [g for g in genes if g in adata.var_names]
    missing = set(genes) - set(valid_genes)
    if missing:
        log(f"    Missing from adata: {len(missing)} genes")

    adata_sub = adata[cells][:, valid_genes]
    if hasattr(adata_sub.X, "toarray"):
        X = adata_sub.X.toarray()
    else:
        X = np.array(adata_sub.X)

    n = len(valid_genes)
    corr = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if i == j:
                corr[i, j] = 1.0
            else:
                rho, _ = spearmanr(X[:, i], X[:, j])
                if np.isnan(rho):
                    rho = 0.0
                corr[i, j] = rho
                corr[j, i] = rho

    return pd.DataFrame(corr, index=valid_genes, columns=valid_genes)


# =============================================================================
# DOT PLOT (3 comparisons only)
# =============================================================================

def generate_dotplot(diff, all_activating, all_inhibiting, valid_genes):
    """Three-comparison dot plot: activating, inhibiting, cross-chain."""
    banner("GENERATING PANEL D (3-comparison dot plot)")

    act = sorted(set(valid_genes) & all_activating)
    rep = sorted(set(valid_genes) & all_inhibiting)

    def get_pairs(group1, group2, exclude_self=True):
        vals = []
        for g1 in group1:
            for g2 in group2:
                if exclude_self and g1 == g2:
                    continue
                if g1 in diff.index and g2 in diff.columns:
                    vals.append(diff.loc[g1, g2])
        return np.array(vals) if vals else np.array([0.0])

    # ---- Build exactly 3 comparisons ----
    comparisons = []

    # 1. Activating chain internal
    if len(act) >= 2:
        vals = get_pairs(act, act)
        n_pos = (vals > 0).sum()
        comparisons.append({
            "label": f"Activating chain\n({n_pos}/{len(vals)} positive edges)",
            "mean": np.mean(vals),
            "median": np.median(vals),
            "values": vals,
            "color": COLOR_ACTIVATING,
        })
        log(f"  Activating internal: mean={np.mean(vals):.4f}, "
            f"median={np.median(vals):.4f}, "
            f"{n_pos}/{len(vals)} positive ({n_pos/len(vals)*100:.1f}%)")

    # 2. Inhibiting chain internal
    if len(rep) >= 2:
        vals = get_pairs(rep, rep)
        n_pos = (vals > 0).sum()
        comparisons.append({
            "label": f"Inhibiting chain\n({n_pos}/{len(vals)} positive edges)",
            "mean": np.mean(vals),
            "median": np.median(vals),
            "values": vals,
            "color": COLOR_INHIBITING,
        })
        log(f"  Inhibiting internal: mean={np.mean(vals):.4f}, "
            f"median={np.median(vals):.4f}, "
            f"{n_pos}/{len(vals)} positive ({n_pos/len(vals)*100:.1f}%)")

    # 3. Cross-chain (activating vs inhibiting)
    if act and rep:
        vals = get_pairs(act, rep, exclude_self=False)
        n_pos = (vals > 0).sum()
        comparisons.append({
            "label": f"Cross-chain\n({n_pos}/{len(vals)} positive edges)",
            "mean": np.mean(vals),
            "median": np.median(vals),
            "values": vals,
            "color": COLOR_CROSSCHAIN,
        })
        log(f"  Cross-chain: mean={np.mean(vals):.4f}, "
            f"median={np.median(vals):.4f}, "
            f"{n_pos}/{len(vals)} positive ({n_pos/len(vals)*100:.1f}%)")

    if not comparisons:
        log("  [ERROR] No comparisons could be computed")
        return

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(14, 8))

    n_comp = len(comparisons)
    y_positions = list(range(n_comp))[::-1]  # top to bottom

    for i, (comp, y) in enumerate(zip(comparisons, y_positions)):
        mean_val = comp["mean"]
        median_val = comp["median"]

        # Dot size proportional to |mean DIFF|
        dot_size = 300 + 3000 * abs(mean_val)

        # Jittered individual edges (very faint)
        jitter = np.random.RandomState(42).normal(0, 0.08, len(comp["values"]))
        ax.scatter(comp["values"], [y] * len(comp["values"]) + jitter,
                   s=15, c=comp["color"], alpha=0.15, zorder=1,
                   rasterized=True)

        # Mean dot
        ax.scatter(mean_val, y, s=dot_size, c=comp["color"],
                   edgecolors="#000000", linewidths=2.0, alpha=0.95,
                   zorder=3)

        # Median diamond
        ax.scatter(median_val, y, s=120, c="white", edgecolors="#000000",
                   linewidths=1.5, marker="D", zorder=4)

        # Value annotation
        ax.text(mean_val, y + 0.30, f"mean = {mean_val:+.3f}",
                ha="center", va="bottom", fontsize=FONT_TICK - 4,
                fontweight="bold", color=comp["color"],
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec="none", alpha=0.8))

    # Zero line
    ax.axvline(0, color="#000000", linewidth=1.5, linestyle="-",
               alpha=0.6, zorder=0)

    # Y-axis
    ax.set_yticks(y_positions)
    ax.set_yticklabels([c["label"] for c in comparisons], fontsize=FONT_LABEL)
    for tick, comp in zip(ax.get_yticklabels(), comparisons):
        tick.set_color(comp["color"])
        tick.set_fontweight("bold")

    # X-axis
    ax.set_xlabel("DIFF  (corr$_{SBS2-HIGH}$ \u2212 corr$_{CNV-HIGH}$)",
                  fontsize=FONT_AXIS)
    ax.tick_params(axis="x", labelsize=FONT_TICK)

    # Shaded regions
    ax.axvspan(-1, 0, alpha=0.04, color=COLOR_INHIBITING, zorder=0)
    ax.axvspan(0, 1, alpha=0.04, color=COLOR_ACTIVATING, zorder=0)
    ax.text(-0.02, -0.7, "Stronger in\nCNV-HIGH",
            ha="right", fontsize=FONT_TICK - 6, color="#2e6e2e",
            style="italic")
    ax.text(0.02, -0.7, "Stronger in\nSBS2-HIGH",
            ha="left", fontsize=FONT_TICK - 6, color="#b5651d",
            style="italic")

    # Axis limits
    all_vals = np.concatenate([c["values"] for c in comparisons])
    x_pad = max(abs(all_vals.min()), abs(all_vals.max())) * 1.3
    ax.set_xlim(-x_pad, x_pad)
    ax.set_ylim(-1.0, n_comp - 0.3)

    # Size legend
    for ref_val, ref_label in [(0.05, "0.05"), (0.10, "0.10"),
                                (0.15, "0.15")]:
        ref_size = 300 + 3000 * ref_val
        ax.scatter([], [], s=ref_size, c="#cccccc", edgecolors="#000000",
                   linewidths=1.5, label=f"|mean DIFF| = {ref_label}")
    ax.scatter([], [], s=120, c="white", edgecolors="#000000",
               linewidths=1.5, marker="D", label="Median")
    ax.legend(loc="upper right", fontsize=FONT_TICK - 6, framealpha=0.9,
              title="Dot size", title_fontsize=FONT_TICK - 6)

    ax.set_title("Chain Co-expression Validation\nSBS2-HIGH vs CNV-HIGH",
                 fontsize=FONT_TITLE, pad=15)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    for ext in ["pdf", "png"]:
        out_path = os.path.join(FIG_DIR, f"Panel_D_v2_dotplot.{ext}")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        log(f"  [SAVE] {out_path}")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("PANEL D v2: Three-Comparison Dot Plot")

    # Load Harris interactors
    harris_df = pd.read_csv(HARRIS_PATH, sep="\t")
    harris_genes = set(harris_df["gene_symbol"].values)
    log(f"Harris A3 interactors: {len(harris_genes)}")

    # Identify chain genes from both networks
    sbs2_dir = os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL")
    cnv_dir = os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL")

    log("Identifying chains from SBS2_VS_NORMAL...")
    act_sbs2, rep_sbs2 = identify_chain_genes(sbs2_dir, harris_genes)
    log(f"  Activating: {len(act_sbs2)}, Inhibiting: {len(rep_sbs2)}")

    log("Identifying chains from CNV_VS_NORMAL...")
    act_cnv, rep_cnv = identify_chain_genes(cnv_dir, harris_genes)
    log(f"  Activating: {len(act_cnv)}, Inhibiting: {len(rep_cnv)}")

    all_activating = act_sbs2 | act_cnv
    all_inhibiting = rep_sbs2 | rep_cnv
    # Resolve overlaps: if a gene is in both, keep as activating
    overlap = all_activating & all_inhibiting
    if overlap:
        log(f"  Overlap resolved (kept as activating): {len(overlap)}")
        all_inhibiting -= overlap

    all_chain_genes = sorted(all_activating | all_inhibiting)
    log(f"  Total chain genes: {len(all_chain_genes)} "
        f"(act={len(all_activating)}, inh={len(all_inhibiting)})")

    # Check for cached DIFF matrix
    cache_path = os.path.join(CACHE_DIR, "DIFF_chain_genes.tsv")
    cache_genes_path = os.path.join(CACHE_DIR, "chain_gene_sets.pkl")

    if os.path.exists(cache_path) and os.path.exists(cache_genes_path):
        log(f"\nLoading cached DIFF matrix from {cache_path}")
        diff = pd.read_csv(cache_path, sep="\t", index_col=0)
        with open(cache_genes_path, "rb") as f:
            cached = pickle.load(f)
        # Verify cache matches current chain genes
        if cached.get("all_chain_genes") == set(all_chain_genes):
            log(f"  Cache valid: {diff.shape[0]} x {diff.shape[1]}")
            valid_genes = list(diff.index)
        else:
            log("  Cache stale (chain genes changed), recomputing...")
            diff = None
    else:
        diff = None

    if diff is None:
        # Compute from scratch
        log("\nLoading adata (this may take a moment)...")
        import scanpy as sc
        adata = sc.read_h5ad(ADATA_PATH)
        log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

        groups_df = pd.read_csv(GROUP_PATH, sep="\t")
        group_map = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))
        adata.obs["three_group"] = adata.obs.index.map(
            lambda x: group_map.get(x, "OTHER"))

        sbs2_cells = list(
            adata.obs.index[adata.obs["three_group"] == "SBS2_HIGH"])
        cnv_cells = list(
            adata.obs.index[adata.obs["three_group"] == "CNV_HIGH"])
        log(f"  SBS2-HIGH cells: {len(sbs2_cells)}")
        log(f"  CNV-HIGH cells: {len(cnv_cells)}")

        valid_genes = [g for g in all_chain_genes if g in adata.var_names]
        log(f"  Valid genes in adata: {len(valid_genes)} / {len(all_chain_genes)}")

        log("  Computing SBS2-HIGH Spearman correlations...")
        corr_sbs2 = compute_corr(adata, sbs2_cells, valid_genes)
        log("  Computing CNV-HIGH Spearman correlations...")
        corr_cnv = compute_corr(adata, cnv_cells, valid_genes)

        # Align
        shared = sorted(set(corr_sbs2.index) & set(corr_cnv.index))
        corr_sbs2 = corr_sbs2.loc[shared, shared]
        corr_cnv = corr_cnv.loc[shared, shared]
        diff = corr_sbs2 - corr_cnv
        valid_genes = shared

        log(f"  DIFF matrix: {diff.shape[0]} x {diff.shape[1]}")

        # Cache
        diff.to_csv(cache_path, sep="\t")
        corr_sbs2.to_csv(os.path.join(CACHE_DIR, "corr_SBS2_chain_genes.tsv"),
                         sep="\t")
        corr_cnv.to_csv(os.path.join(CACHE_DIR, "corr_CNV_chain_genes.tsv"),
                        sep="\t")
        with open(cache_genes_path, "wb") as f:
            pickle.dump({"all_chain_genes": set(all_chain_genes),
                         "all_activating": all_activating,
                         "all_inhibiting": all_inhibiting}, f)
        log(f"  Cached to {CACHE_DIR}")

        del adata

    # Generate the plot
    generate_dotplot(diff, all_activating, all_inhibiting, valid_genes)

    banner("PANEL D v2 COMPLETE")
    log(f"\nOutput: {FIG_DIR}")
    for f in sorted(os.listdir(FIG_DIR)):
        if "Panel_D" in f:
            log(f"  {f}")


if __name__ == "__main__":
    main()
