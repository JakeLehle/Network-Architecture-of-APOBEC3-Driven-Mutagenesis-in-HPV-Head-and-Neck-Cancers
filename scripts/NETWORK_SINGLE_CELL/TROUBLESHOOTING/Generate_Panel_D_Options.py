#!/usr/bin/env python3
"""
Generate_Panel_D_Options.py
=============================

Two candidate visualizations for Figure 4 panel d (SBS2 vs CNV
chain validation):

Option B: Focused gene heatmap with hierarchically sorted genes
          within activating/inhibiting/A3 groups.

Option C: Grouped dot plot showing mean DIFF with dot size scaled
          by absolute correlation strength.

Usage:
    conda run -n NETWORK python Generate_Panel_D_Options.py
"""

import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")
ADATA_PATH = os.path.join(FIG4_ROOT, "00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(FIG4_ROOT,
                          "01_group_selection/three_group_assignments.tsv")
FIG_DIR = os.path.join(FIG4_ROOT, "FIGURE_4_PANELS")
os.makedirs(FIG_DIR, exist_ok=True)

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
            "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
            "APOBEC3H": "A3H"}

# Style
COLOR_UP_TUMOR  = "#fcad61"
COLOR_UP_NORMAL = "#aadda4"
COLOR_A3        = "#ed6a5a"
COLOR_HARRIS    = "#F6D155"
FONT_TITLE  = 32
FONT_AXIS   = 28
FONT_LEGEND = 22
FONT_TICK   = 20


def log(msg):
    print(msg, flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# CHAIN IDENTIFICATION (same logic as figure script)
# =============================================================================

def identify_chain_genes(net_dir, harris_genes):
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
# COMPUTE CORRELATIONS
# =============================================================================

def compute_corr(adata, cells, genes):
    sub = adata[cells][:, genes]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.array(sub.X)
    n = len(genes)
    mat = np.zeros((n, n))
    for i in range(n):
        mat[i, i] = 1.0
        for j in range(i + 1, n):
            rho, _ = spearmanr(X[:, i], X[:, j])
            if np.isnan(rho):
                rho = 0.0
            mat[i, j] = rho
            mat[j, i] = rho
    return pd.DataFrame(mat, index=genes, columns=genes)


# =============================================================================
# SELECT FOCUSED GENES FOR OPTION B
# =============================================================================

def select_focused_genes(all_activating, all_inhibiting, harris_genes,
                         diff, valid_genes):
    """Select a focused set of genes for the heatmap.
    Activating: all 9 (small enough already).
    Inhibiting: Harris interactors + boundary genes + top hub genes.
    A3: A3A, A3B.
    """
    act_present = sorted(set(valid_genes) & all_activating)

    # For inhibiting: take Harris interactors first, then genes with
    # strongest mean absolute DIFF to A3
    rep_present = sorted(set(valid_genes) & all_inhibiting)
    a3_genes = ["APOBEC3A", "APOBEC3B"]
    a3_present = [g for g in a3_genes if g in valid_genes]

    # Known boundary genes from the diagnostic
    boundary_genes = {"TIMM8B", "PRSS27", "IL1B", "CA12", "LGALSL",
                      "SPRR1B", "IL1RN", "TMEM45B", "SNHG3", "THYN1",
                      "KRT24", "CCL20", "RRAD"}

    rep_harris = [g for g in rep_present if g in harris_genes]
    rep_boundary = [g for g in rep_present if g in boundary_genes
                    and g not in rep_harris]

    # Top remaining by mean |DIFF| to A3
    rep_used = set(rep_harris) | set(rep_boundary)
    rep_remaining = [g for g in rep_present if g not in rep_used]
    if rep_remaining and a3_present:
        scored = []
        for g in rep_remaining:
            mean_abs = np.mean([abs(diff.loc[g, a]) for a in a3_present
                                if a in diff.index and g in diff.index])
            scored.append((g, mean_abs))
        scored.sort(key=lambda x: -x[1])
        rep_top = [g for g, _ in scored[:8]]
    else:
        rep_top = []

    focused_rep = sorted(set(rep_harris) | set(rep_boundary) | set(rep_top))
    # Cap at ~20 inhibiting genes
    focused_rep = focused_rep[:20]

    return act_present, focused_rep, a3_present


# =============================================================================
# HIERARCHICAL SORT WITHIN GROUPS
# =============================================================================

def hierarchical_sort(genes, diff_matrix):
    """Sort genes within a group by hierarchical clustering on DIFF."""
    if len(genes) <= 2:
        return genes
    sub = diff_matrix.loc[genes, genes].values
    # Convert correlation-like matrix to distance
    dist = 1.0 - sub
    np.fill_diagonal(dist, 0)
    # Make symmetric and fix numerical issues
    dist = (dist + dist.T) / 2
    dist = np.clip(dist, 0, None)
    try:
        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        order = leaves_list(Z)
        return [genes[i] for i in order]
    except Exception as e:
        log(f"    [WARNING] Hierarchical sort failed: {e}")
        return genes


# =============================================================================
# OPTION B: FOCUSED GENE HEATMAP
# =============================================================================

def generate_option_b(diff, all_activating, all_inhibiting, harris_genes,
                      valid_genes):
    """Focused heatmap with hierarchically sorted genes."""
    banner("[OPTION B] Focused Gene Heatmap")

    act_present, rep_focused, a3_present = select_focused_genes(
        all_activating, all_inhibiting, harris_genes, diff, valid_genes)

    log(f"  Activating: {len(act_present)} genes")
    log(f"  Inhibiting (focused): {len(rep_focused)} genes")
    log(f"  A3: {len(a3_present)} genes")

    # Hierarchical sort within each group
    act_sorted = hierarchical_sort(act_present, diff)
    rep_sorted = hierarchical_sort(rep_focused, diff)
    ordered = act_sorted + rep_sorted + a3_present

    diff_ordered = diff.loc[ordered, ordered]

    fig, ax = plt.subplots(figsize=(18, 16))

    vmax = max(abs(diff_ordered.values.min()),
               abs(diff_ordered.values.max()), 0.3)
    im = ax.imshow(diff_ordered.values, cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="auto")

    display_names = []
    for g in ordered:
        name = A3_ALIAS.get(g, g)
        if g in harris_genes:
            name = name + " *"
        display_names.append(name)

    ax.set_xticks(range(len(ordered)))
    ax.set_xticklabels(display_names, rotation=90, fontsize=11)
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels(display_names, fontsize=11)

    # Color-code labels
    for g, ytick, xtick in zip(ordered, ax.get_yticklabels(),
                                ax.get_xticklabels()):
        if g in all_activating:
            color = "#b5651d"  # dark orange
        elif g in all_inhibiting:
            color = "#2e6e2e"  # dark green
        else:
            color = COLOR_A3
        ytick.set_color(color)
        ytick.set_fontweight("bold")
        xtick.set_color(color)
        xtick.set_fontweight("bold")

    # Divider lines
    n_act = len(act_sorted)
    n_rep = len(rep_sorted)
    for pos_line in [n_act - 0.5, n_act + n_rep - 0.5]:
        ax.axhline(pos_line, color="black", linewidth=2.5)
        ax.axvline(pos_line, color="black", linewidth=2.5)

    # Group labels on right side
    if n_act > 0:
        ax.text(len(ordered) + 0.8, n_act / 2 - 0.5, "A3\nActivating",
                fontsize=FONT_LEGEND - 2, fontweight="bold",
                color="#b5651d", ha="left", va="center")
    if n_rep > 0:
        ax.text(len(ordered) + 0.8, n_act + n_rep / 2 - 0.5,
                "A3\nInhibiting",
                fontsize=FONT_LEGEND - 2, fontweight="bold",
                color="#2e6e2e", ha="left", va="center")
    ax.text(len(ordered) + 0.8, n_act + n_rep + len(a3_present) / 2 - 0.5,
            "APOBEC3", fontsize=FONT_LEGEND - 2, fontweight="bold",
            color=COLOR_A3, ha="left", va="center")

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.12)
    cbar.set_label("DIFF (corr$_{SBS2}$ \u2212 corr$_{CNV}$)",
                   fontsize=FONT_LEGEND)
    cbar.ax.tick_params(labelsize=FONT_TICK - 2)

    ax.set_title("SBS2-HIGH vs CNV-HIGH\nChain Gene Co-expression Validation",
                 fontsize=FONT_TITLE, pad=15)

    # Footnote
    ax.text(0.0, -0.06, "* Harris A3 interactor",
            transform=ax.transAxes, fontsize=14, style="italic",
            color="#666666")

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR,
                    f"Panel_D_option_B_heatmap.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Option B -> Panel_D_option_B_heatmap.pdf/png")


# =============================================================================
# OPTION C: GROUPED DOT PLOT
# =============================================================================

def generate_option_c(diff, all_activating, all_inhibiting, valid_genes):
    """Grouped dot plot: mean DIFF for each block pair, dot size = |DIFF|."""
    banner("[OPTION C] Grouped Dot Plot")

    act = sorted(set(valid_genes) & all_activating)
    rep = sorted(set(valid_genes) & all_inhibiting)
    a3a = "APOBEC3A" if "APOBEC3A" in valid_genes else None
    a3b = "APOBEC3B" if "APOBEC3B" in valid_genes else None

    def get_pairs(group1, group2, exclude_self=True):
        vals = []
        for g1 in group1:
            for g2 in group2:
                if exclude_self and g1 == g2:
                    continue
                if g1 in diff.index and g2 in diff.columns:
                    vals.append(diff.loc[g1, g2])
        return np.array(vals) if vals else np.array([0.0])

    # Define the comparisons
    comparisons = []

    # Within-group
    if len(act) >= 2:
        vals = get_pairs(act, act)
        comparisons.append({
            "label": f"Activating\nwithin\n(n={len(vals)})",
            "short": "Act-Act",
            "mean": np.mean(vals), "median": np.median(vals),
            "values": vals,
            "color": COLOR_UP_TUMOR, "group": "within",
            "prediction": "Positive\n(SBS2-specific)",
        })

    if len(rep) >= 2:
        vals = get_pairs(rep, rep)
        comparisons.append({
            "label": f"Inhibiting\nwithin\n(n={len(vals)})",
            "short": "Inh-Inh",
            "mean": np.mean(vals), "median": np.median(vals),
            "values": vals,
            "color": COLOR_UP_NORMAL, "group": "within",
            "prediction": "Gradient\n(NORM>SBS2>CNV)",
        })

    # Cross-group
    if act and rep:
        vals = get_pairs(act, rep, exclude_self=False)
        comparisons.append({
            "label": f"Activating\nvs Inhibiting\n(n={len(vals)})",
            "short": "Act-Inh",
            "mean": np.mean(vals), "median": np.median(vals),
            "values": vals,
            "color": "#888888", "group": "cross",
            "prediction": "Near zero\n(independent)",
        })

    # A3A to groups
    if a3a:
        if act:
            vals = get_pairs([a3a], act, exclude_self=False)
            comparisons.append({
                "label": f"A3A to\nActivating\n(n={len(vals)})",
                "short": "A3A-Act",
                "mean": np.mean(vals), "median": np.median(vals),
                "values": vals,
                "color": COLOR_A3, "group": "A3A",
                "prediction": "Negative\n(decoupled)",
            })
        if rep:
            vals = get_pairs([a3a], rep, exclude_self=False)
            comparisons.append({
                "label": f"A3A to\nInhibiting\n(n={len(vals)})",
                "short": "A3A-Inh",
                "mean": np.mean(vals), "median": np.median(vals),
                "values": vals,
                "color": COLOR_A3, "group": "A3A",
                "prediction": "Negative\n(decoupled)",
            })

    # A3B to groups
    if a3b:
        if act:
            vals = get_pairs([a3b], act, exclude_self=False)
            comparisons.append({
                "label": f"A3B to\nActivating\n(n={len(vals)})",
                "short": "A3B-Act",
                "mean": np.mean(vals), "median": np.median(vals),
                "values": vals,
                "color": "#c0392b", "group": "A3B",
                "prediction": "Positive\n(recoupling)",
            })
        if rep:
            vals = get_pairs([a3b], rep, exclude_self=False)
            comparisons.append({
                "label": f"A3B to\nInhibiting\n(n={len(vals)})",
                "short": "A3B-Inh",
                "mean": np.mean(vals), "median": np.median(vals),
                "values": vals,
                "color": "#c0392b", "group": "A3B",
                "prediction": "Mild positive\n(partial)",
            })

    # Build the plot
    fig, ax = plt.subplots(figsize=(16, 12))

    n_comp = len(comparisons)
    y_positions = list(range(n_comp))[::-1]  # top to bottom

    for i, (comp, y) in enumerate(zip(comparisons, y_positions)):
        mean_val = comp["mean"]
        abs_mean = abs(mean_val)
        vals = comp["values"]

        # Dot size scales with absolute mean DIFF
        dot_size = 200 + 2000 * abs_mean

        # Individual data points (jittered, small, transparent)
        if len(vals) <= 100:
            jitter = np.random.default_rng(42).normal(0, 0.08, len(vals))
            ax.scatter(vals, [y] * len(vals) + jitter,
                       s=15, c=comp["color"], alpha=0.25, zorder=1)
        else:
            # Too many points: show percentile whiskers instead
            p25, p75 = np.percentile(vals, [25, 75])
            p5, p95 = np.percentile(vals, [5, 95])
            ax.plot([p5, p95], [y, y], color=comp["color"],
                    alpha=0.4, linewidth=2, zorder=1)
            ax.plot([p25, p75], [y, y], color=comp["color"],
                    alpha=0.7, linewidth=5, zorder=2)

        # Mean dot (large, prominent)
        ax.scatter([mean_val], [y], s=dot_size, c=comp["color"],
                   edgecolors="#000000", linewidths=2.0, alpha=0.95,
                   zorder=4)

        # Median marker
        ax.scatter([comp["median"]], [y], s=80, c="white",
                   edgecolors="#000000", linewidths=1.0, marker="D",
                   alpha=0.9, zorder=5)

        # Value annotation
        ax.text(mean_val, y + 0.35, f"{mean_val:+.3f}",
                ha="center", va="bottom", fontsize=FONT_TICK,
                fontweight="bold", color=comp["color"],
                bbox=dict(boxstyle="round,pad=0.1", fc="white",
                          ec="none", alpha=0.8),
                zorder=6)

        # Prediction text on right
        ax.text(0.82, y, comp["prediction"],
                ha="left", va="center", fontsize=13,
                color="#555555", style="italic",
                transform=ax.get_yaxis_transform())

    # Zero line
    ax.axvline(0, color="#000000", linewidth=1.5, linestyle="-",
               alpha=0.6, zorder=0)

    # Axis setup
    ax.set_yticks(y_positions)
    ax.set_yticklabels([c["label"] for c in comparisons],
                       fontsize=FONT_TICK)

    # Color the y-axis labels
    for tick, comp in zip(ax.get_yticklabels(), comparisons):
        tick.set_color(comp["color"])
        tick.set_fontweight("bold")

    ax.set_xlabel("DIFF (corr$_{SBS2}$ \u2212 corr$_{CNV}$)",
                  fontsize=FONT_AXIS)
    ax.tick_params(axis="x", labelsize=FONT_TICK)

    # Shade regions
    ax.axvspan(-1, 0, alpha=0.04, color=COLOR_UP_NORMAL, zorder=0)
    ax.axvspan(0, 1, alpha=0.04, color=COLOR_UP_TUMOR, zorder=0)
    ax.text(-0.02, -0.8, "Stronger in CNV-HIGH",
            ha="right", fontsize=14, color="#2e6e2e", style="italic")
    ax.text(0.02, -0.8, "Stronger in SBS2-HIGH",
            ha="left", fontsize=14, color="#b5651d", style="italic")

    ax.set_xlim(-0.5, 0.8)
    ax.set_ylim(-1.2, n_comp - 0.3)

    # Size legend
    for ref_val, ref_label in [(0.05, "0.05"), (0.10, "0.10"),
                                (0.15, "0.15")]:
        ref_size = 200 + 2000 * ref_val
        ax.scatter([], [], s=ref_size, c="#cccccc", edgecolors="#000000",
                   linewidths=1, label=f"|DIFF| = {ref_label}")
    ax.scatter([], [], s=80, c="white", edgecolors="#000000",
               linewidths=1, marker="D", label="Median")
    ax.legend(loc="upper right", fontsize=14, framealpha=0.9,
              title="Dot size = |mean DIFF|", title_fontsize=14)

    ax.set_title("SBS2-HIGH vs CNV-HIGH\n"
                 "Chain Co-expression Validation",
                 fontsize=FONT_TITLE, pad=15)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(FIG_DIR,
                    f"Panel_D_option_C_dotplot.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log(f"  [SAVE] Option C -> Panel_D_option_C_dotplot.pdf/png")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("PANEL D OPTIONS: B (heatmap) and C (dot plot)")

    # Load Harris
    harris_df = pd.read_csv(HARRIS_PATH, sep="\t")
    harris_genes = set(harris_df["gene_symbol"].values)
    log(f"Harris A3 interactors: {len(harris_genes)}")

    # Identify chain genes
    sbs2_dir = os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL")
    cnv_dir = os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL")

    act_sbs2, rep_sbs2 = identify_chain_genes(sbs2_dir, harris_genes)
    act_cnv, rep_cnv = identify_chain_genes(cnv_dir, harris_genes)

    all_activating = act_sbs2 | act_cnv
    all_inhibiting = rep_sbs2 | rep_cnv
    overlap = all_activating & all_inhibiting
    all_inhibiting -= overlap

    a3_genes = {"APOBEC3A", "APOBEC3B"}
    all_chain_genes = sorted(all_activating | all_inhibiting | a3_genes)

    log(f"  Activating: {len(all_activating)}")
    log(f"  Inhibiting: {len(all_inhibiting)}")
    log(f"  Total: {len(all_chain_genes)}")

    # Load expression and compute DIFF
    log("\nLoading adata...")
    adata = sc.read_h5ad(ADATA_PATH)
    groups_df = pd.read_csv(GROUP_PATH, sep="\t")
    group_map = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))
    adata.obs["three_group"] = adata.obs.index.map(
        lambda x: group_map.get(x, "OTHER"))

    sbs2_cells = list(
        adata.obs.index[adata.obs["three_group"] == "SBS2_HIGH"])
    cnv_cells = list(
        adata.obs.index[adata.obs["three_group"] == "CNV_HIGH"])

    valid_genes = [g for g in all_chain_genes if g in adata.var_names]
    log(f"  Valid genes in adata: {len(valid_genes)}")

    log("  Computing SBS2-HIGH correlations...")
    corr_sbs2 = compute_corr(adata, sbs2_cells, valid_genes)
    log("  Computing CNV-HIGH correlations...")
    corr_cnv = compute_corr(adata, cnv_cells, valid_genes)

    diff = corr_sbs2 - corr_cnv
    del adata
    log(f"  DIFF matrix: {diff.shape[0]} x {diff.shape[1]}")

    # Generate both options
    generate_option_b(diff, all_activating, all_inhibiting,
                      harris_genes, valid_genes)
    generate_option_c(diff, all_activating, all_inhibiting, valid_genes)

    banner("PANEL D OPTIONS COMPLETE")
    log(f"\nOutput: {FIG_DIR}")
    for f in sorted(os.listdir(FIG_DIR)):
        if "Panel_D" in f:
            log(f"  {f}")


if __name__ == "__main__":
    main()
