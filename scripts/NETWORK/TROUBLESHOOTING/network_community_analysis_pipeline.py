'''

# APOBEC3 – SBS2 Network Analysis Pipeline

Author: Dr. Mohadeseh Soleimanpour
Email: [msoleimanpour@txbiomed.org]
Institution: Texas Biomedical Research Institute

---

## Pipeline Overview

This pipeline analyzes the relationship between APOBEC activity
(SBS2 mutational signature) and gene expression patterns in TCGA
tumor samples. The workflow integrates mutation signatures,
gene expression, statistical testing, correlation networks,
and community detection to identify gene modules associated
with APOBEC-driven mutagenesis.

The analysis focuses particularly on APOBEC3 family genes and
selected immune and tumor microenvironment biomarkers.

---

## Main Analysis Steps

STEP 01 – Define input data paths
Load TCGA gene expression data and the SBS mutation
signature dataset.

STEP 02 – Load TCGA gene expression matrix
Import raw TCGA expression file and inspect structure.

STEP 03 – Remove _PAR_Y genes
Remove pseudoautosomal Y-linked genes to avoid redundant
duplicated genomic regions.

STEP 04 – Assign column headers
Extract clinical metadata columns and gene expression
columns from the raw dataset.

STEP 05 – Extract gene annotations
Retrieve ENSG identifiers, gene symbols, and gene biotype
information from annotation rows.

STEP 06 – Clean ENSG identifiers
Remove version numbers (e.g., ENSG000001.15 → ENSG000001)
to standardize gene identifiers.

STEP 07 – Build expression matrix
Separate clinical metadata from the numeric expression matrix.

STEP 08 – Resolve duplicated genes
Merge duplicate ENSG columns created by version removal.

STEP 09 – Build ENSG → gene symbol mapping
Create a dictionary used for labeling plots and outputs.

STEP 10 – Merge TCGA expression with SBS mutation signatures
Match samples using TCGA barcode identifiers.

STEP 10B – Remove duplicated samples
Drop samples appearing multiple times after merge to avoid
statistical bias.

STEP 11 – Identify available cancer types
Extract and save cancer types present in the dataset.

STEP 12 – Differential expression analysis
Compare gene expression between:
High SBS2 samples
Low SBS2 samples
using Welch's t-test.

STEP 13 – Define APOBEC activity groups
Compute an APOBEC activity score using A3A and A3B expression,
then define:

```
    TOP group    → high SBS2 within high APOBEC score samples
    BOTTOM group → low SBS2 within high APOBEC score samples
```

STEP 14 – Correlation network construction
Build gene co-expression networks for:

```
    TOP samples
    BOTTOM samples
    DIFF network (difference between TOP and BOTTOM)
```

STEP 15 – Community detection
Identify network modules using Leiden clustering across
multiple resolution parameters and evaluate stability.

STEP 16 – Network metrics
Compute graph-based centrality measures including:

```
    degree
    betweenness
    closeness
    eigenvector centrality
```

These metrics highlight genes whose network roles change
between high and low APOBEC mutational activity.

---

## Outputs

The pipeline produces:

• Differential expression tables
• Volcano and Manhattan plots
• Correlation matrices
• Network edge lists
• Network visualizations
• Community gene modules
• Graph metric tables
• Zoomed ego-network visualizations

All results are saved in a structured directory hierarchy
organized by cancer type and analysis step.

---

## Important Notes

• Gene identifiers are kept as ENSG IDs for analysis consistency.
• Gene symbols are used only for visualization and reporting.
• APOBEC3 genes are always retained in the analysis even if
they are not statistically significant in the initial test.
• Biomarker genes are included when present but are not forced
into the final gene selection.

---

## End of Header

'''



A3_ALIASES = ["A3A", "A3B", "A3C", "A3D","A3F", "A3G", "A3H"]

def symbol_or_self(x: str) -> str:
    return ensg_to_symbol.get(str(x), str(x))

def relabel_graph_and_pos_to_symbols(G, pos):
    """
    Plot-only relabel:
    ENSG nodes -> gene symbols (when available).
    Returns (G_sym, pos_sym, name_map)
    """
    name_map = {n: symbol_or_self(n) for n in G.nodes()}
    G_sym = nx.relabel_nodes(G, name_map, copy=True)
    pos_sym = {name_map.get(n, n): p for n, p in pos.items() if n in name_map}
    return G_sym, pos_sym, name_map


def get_a3_nodes_for_graph(G):
    # ENSG mode
    if any(g in G.nodes for g in A3_GENES):
        return [g for g in A3_GENES if g in G.nodes]
    # alias mode (rare in your pipeline, but safe)
    if any(a in G.nodes for a in A3_ALIASES):
        return [a for a in A3_ALIASES if a in G.nodes]
    return []

def get_biomarker_nodes_for_graph(G):
    return [b for b in BIOMARKERS if b in G.nodes]

def get_highlight_nodes(G):
    # A3 + biomarkers-if-present
    return list(dict.fromkeys(get_a3_nodes_for_graph(G) + get_biomarker_nodes_for_graph(G)))

def add_nodes_if_missing(G, nodes):
    for n in nodes:
        if n not in G:
            G.add_node(n)

def remove_isolates_but_keep(G, keep_nodes):
    """Remove degree-0 nodes EXCEPT the ones we want to keep visible."""
    H = G.copy()
    for n in list(H.nodes()):
        if H.degree(n) == 0 and n not in keep_nodes:
            H.remove_node(n)
    return H



def relabel_graph_pos_and_partition_to_symbols(G, pos, cm, ensg_to_symbol):
    """
    Relabel graph nodes, position dict, and community map from ENSG -> symbol
    for plotting only.

    Returns:
        G_sym, pos_sym, cm_sym, name_map
    """
    name_map = {n: ensg_to_symbol.get(str(n), str(n)) for n in G.nodes()}

    # graph
    G_sym = nx.relabel_nodes(G, name_map, copy=True)

    # positions
    pos_sym = {name_map[n]: p for n, p in pos.items() if n in name_map}

    # partition / community map
    cm_sym = {name_map[n]: c for n, c in cm.items() if n in name_map}

    return G_sym, pos_sym, cm_sym, name_map



def plot_network_by_community_fixedpos_symbols(
    G,
    pos,
    community_map,
    out_path,
    ensg_to_symbol,
    a3_genes=None,
    title=None,
    weighted=True
):
    """
    Wrapper around plot_network_by_community_fixedpos that converts
    ENSG node IDs -> gene symbols for plotting only.

    The underlying graph and analysis remain ENSG-based.
    """

    # -------------------------------
    # Build mapping ENSG -> symbol
    # -------------------------------
    name_map = {
        n: ensg_to_symbol.get(str(n), str(n))
        for n in G.nodes()
    }

    # -------------------------------
    # Relabel graph
    # -------------------------------
    G_sym = nx.relabel_nodes(G, name_map, copy=True)

    # -------------------------------
    # Relabel positions
    # -------------------------------
    pos_sym = {
        name_map[n]: p
        for n, p in pos.items()
        if n in name_map
    }

    # -------------------------------
    # Relabel community map
    # -------------------------------
    cm_sym = {
        name_map[n]: c
        for n, c in community_map.items()
        if n in name_map
    }

    # -------------------------------
    # Relabel highlight genes
    # -------------------------------
    highlight_sym = None
    if a3_genes is not None:
        highlight_sym = [
            name_map[g] for g in a3_genes
            if g in name_map
        ]

    # -------------------------------
    # Call ORIGINAL plotting function
    # -------------------------------
    plot_network_by_community_fixedpos(
        G_sym,
        pos_sym,
        cm_sym,
        out_path=out_path,
        title=title,
        a3_genes=highlight_sym,
        weighted=weighted
    )


def save_gene_lists_from_map_symbols_only(cm_map: dict, out_csv: str, ensg_to_symbol: dict):
    comm_to_genes = {}
    for g, c in cm_map.items():
        comm_to_genes.setdefault(c, []).append(g)

    rows = []
    for c, genes in comm_to_genes.items():
        genes_sorted = sorted(genes)
        symbols_sorted = [ensg_to_symbol.get(g, g) for g in genes_sorted]

        rows.append({
            "community": c,
            "size": len(genes_sorted),
            "genes": ";".join(map(str, symbols_sorted)),
        })

    pd.DataFrame(rows).sort_values("size", ascending=False).to_csv(out_csv, index=False)

def save_gene_lists_from_map_symbols(cm_map: dict, out_csv: str, ensg_to_symbol: dict):
    """
    Save community gene lists using gene symbols.
    """
    comm_to_genes = {}

    for g, c in cm_map.items():
        comm_to_genes.setdefault(c, []).append(g)

    rows = []

    for c, genes in comm_to_genes.items():
        genes_sorted = sorted(genes)
        symbols_sorted = [ensg_to_symbol.get(g, g) for g in genes_sorted]

        rows.append({
            "community": c,
            "size": len(genes_sorted),
            "genes_symbol": ";".join(map(str, symbols_sorted))
        })

    pd.DataFrame(rows).sort_values("size", ascending=False).to_csv(out_csv, index=False)

def relabel_graph_and_pos_to_symbols_safe(G, pos, ensg_to_symbol):
    """
    Plot-only relabel:
    ENSG nodes -> gene symbols when available.
    Keeps original name if symbol missing.
    """
    name_map = {n: ensg_to_symbol.get(str(n), str(n)) for n in G.nodes()}
    G_sym = nx.relabel_nodes(G, name_map, copy=True)
    pos_sym = {name_map[n]: p for n, p in pos.items() if n in name_map}
    return G_sym, pos_sym, name_map


def compute_basic_node_metrics(G):
    """
    Return a dataframe of common node metrics.
    """
    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=[
            "gene", "degree", "degree_centrality",
            "betweenness", "closeness", "eigenvector"
        ])

    degree_dict = dict(G.degree())
    degree_cent = nx.degree_centrality(G)

    if G.number_of_edges() > 0:
        betweenness = nx.betweenness_centrality(G, weight=None)
        closeness = nx.closeness_centrality(G)
        try:
            eigenvector = nx.eigenvector_centrality_numpy(G)
        except Exception:
            eigenvector = {n: np.nan for n in G.nodes()}
    else:
        betweenness = {n: 0.0 for n in G.nodes()}
        closeness = {n: 0.0 for n in G.nodes()}
        eigenvector = {n: np.nan for n in G.nodes()}

    rows = []
    for n in G.nodes():
        rows.append({
            "gene": n,
            "degree": degree_dict.get(n, 0),
            "degree_centrality": degree_cent.get(n, 0.0),
            "betweenness": betweenness.get(n, 0.0),
            "closeness": closeness.get(n, 0.0),
            "eigenvector": eigenvector.get(n, np.nan),
        })

    return pd.DataFrame(rows)


def build_metric_change_table(G_top, G_bot, ensg_to_symbol=None):
    """
    Compare node metrics between TOP and BOTTOM.
    """
    mt = compute_basic_node_metrics(G_top).rename(columns={
        "degree": "degree_top",
        "degree_centrality": "degree_centrality_top",
        "betweenness": "betweenness_top",
        "closeness": "closeness_top",
        "eigenvector": "eigenvector_top",
    })

    mb = compute_basic_node_metrics(G_bot).rename(columns={
        "degree": "degree_bottom",
        "degree_centrality": "degree_centrality_bottom",
        "betweenness": "betweenness_bottom",
        "closeness": "closeness_bottom",
        "eigenvector": "eigenvector_bottom",
    })

    m = pd.merge(mt, mb, on="gene", how="outer").fillna(0)

    m["abs_degree_diff"] = (m["degree_top"] - m["degree_bottom"]).abs()
    m["abs_betweenness_diff"] = (m["betweenness_top"] - m["betweenness_bottom"]).abs()
    m["abs_closeness_diff"] = (m["closeness_top"] - m["closeness_bottom"]).abs()
    m["abs_eigenvector_diff"] = (m["eigenvector_top"] - m["eigenvector_bottom"]).abs()

    m["importance_score"] = (
        m["abs_degree_diff"].fillna(0) +
        m["abs_betweenness_diff"].fillna(0) +
        m["abs_closeness_diff"].fillna(0) +
        m["abs_eigenvector_diff"].fillna(0)
    )

    if ensg_to_symbol is not None:
        m["gene_symbol"] = m["gene"].map(lambda x: ensg_to_symbol.get(str(x), str(x)))

    return m.sort_values("importance_score", ascending=False).reset_index(drop=True)



def plot_zoomed_metric_bars(
    gene,
    metric_df,
    out_path,
    ensg_to_symbol=None
):
    """
    Bar plots comparing TOP vs BOTTOM metrics for one gene.
    """
    row = metric_df.loc[metric_df["gene"] == gene]
    if row.empty:
        return

    row = row.iloc[0]
    gene_label = ensg_to_symbol.get(str(gene), str(gene)) if ensg_to_symbol else str(gene)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    metrics = [
        ("degree_top", "degree_bottom", "Degree"),
        ("betweenness_top", "betweenness_bottom", "Betweenness"),
        ("closeness_top", "closeness_bottom", "Closeness"),
        ("eigenvector_top", "eigenvector_bottom", "Eigenvector"),
    ]

    for ax, (c1, c2, title) in zip(axes.flatten(), metrics):
        ax.bar(["TOP", "BOTTOM"], [row[c1], row[c2]])
        ax.set_title(f"{title} | {gene_label}")
        ax.set_ylabel(title)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()



def relabel_subgraph_to_symbols(G_sub, ensg_to_symbol):
    """
    Convert ENSG node labels in a subgraph to gene symbols for plotting.
    """
    name_map = {n: ensg_to_symbol.get(str(n), str(n)) for n in G_sub.nodes()}
    G_sym = nx.relabel_nodes(G_sub, name_map, copy=True)
    return G_sym, name_map


def plot_zoomed_ego_triptych(
    gene,
    G_top,
    G_bot,
    G_diff,
    out_path,
    ensg_to_symbol=None,
    use_symbols=False
):
    """
    Plot ego networks around one gene in TOP / BOTTOM / DIFF side-by-side
    using ONE shared layout across all three panels.

    Layout priority:
      1) ego graph from DIFF (if gene is present there)
      2) otherwise ego graph from union(TOP, BOTTOM, DIFF)
    """

    # ------------------------------------------------------------
    # Build ego subgraphs
    # ------------------------------------------------------------
    sub_top = nx.ego_graph(G_top, gene, radius=1) if gene in G_top.nodes else nx.Graph()
    sub_bot = nx.ego_graph(G_bot, gene, radius=1) if gene in G_bot.nodes else nx.Graph()
    sub_diff = nx.ego_graph(G_diff, gene, radius=1) if gene in G_diff.nodes else nx.Graph()

    gene_label = ensg_to_symbol.get(str(gene), str(gene)) if ensg_to_symbol else str(gene)

    # ------------------------------------------------------------
    # Build union graph for a shared node set
    # ------------------------------------------------------------
    sub_union = nx.compose(nx.compose(sub_top, sub_bot), sub_diff)

    # nothing to plot
    if sub_union.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_title(f"{gene_label}: not present in TOP/BOTTOM/DIFF")
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()
        return

    # ------------------------------------------------------------
    # Shared layout anchor:
    # prefer DIFF ego graph if available, otherwise union
    # ------------------------------------------------------------
    if sub_diff.number_of_nodes() > 0:
        pos_shared = nx.spring_layout(sub_diff, seed=42)
        # fill missing nodes from union while preserving existing positions
        missing_nodes = [n for n in sub_union.nodes if n not in pos_shared]
        if missing_nodes:
            sub_union_with_pos = sub_union.copy()
            pos_shared = nx.spring_layout(
                sub_union_with_pos,
                seed=42,
                pos=pos_shared,
                fixed=list(pos_shared.keys())
            )
    else:
        pos_shared = nx.spring_layout(sub_union, seed=42)

    # ------------------------------------------------------------
    # Optional symbol relabeling
    # ------------------------------------------------------------
    if use_symbols and ensg_to_symbol is not None:
        def relabel_subgraph(G_sub):
            name_map = {n: ensg_to_symbol.get(str(n), str(n)) for n in G_sub.nodes()}
            G_sub_sym = nx.relabel_nodes(G_sub, name_map, copy=True)
            return G_sub_sym, name_map

        sub_top_plot, map_top = relabel_subgraph(sub_top)
        sub_bot_plot, map_bot = relabel_subgraph(sub_bot)
        sub_diff_plot, map_diff = relabel_subgraph(sub_diff)

        # relabel shared positions into symbol space
        union_name_map = {n: ensg_to_symbol.get(str(n), str(n)) for n in sub_union.nodes()}
        pos_shared_plot = {union_name_map[n]: p for n, p in pos_shared.items() if n in union_name_map}
    else:
        sub_top_plot = sub_top
        sub_bot_plot = sub_bot
        sub_diff_plot = sub_diff
        pos_shared_plot = pos_shared

    # ------------------------------------------------------------
    # Plot three panels with the SAME positions
    # ------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, G_sub, title in [
        (axes[0], sub_top_plot, "TOP"),
        (axes[1], sub_bot_plot, "BOTTOM"),
        (axes[2], sub_diff_plot, "DIFF"),
    ]:
        if G_sub.number_of_nodes() == 0:
            ax.set_title(f"{title}: not present")
            ax.axis("off")
            continue

        pos_sub = {n: pos_shared_plot[n] for n in G_sub.nodes() if n in pos_shared_plot}

        nx.draw(
            G_sub,
            pos_sub,
            ax=ax,
            with_labels=True,
            node_size=350,
            font_size=7
        )

        ax.set_title(f"{title} around {gene_label}")
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()




# =============================================================================
# STEP 00 — Imports + Logger
# =============================================================================
import os
import re
import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from dataclasses import dataclass
import json
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from sklearn.decomposition import NMF
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score



from network_utils import (
    # logging
    log as ulog,

    # selection helpers
    bh_fdr,
    compute_a3_score_AB,
    plot_a3_scatter_and_distribution,
    plot_sbs_vs_a3_selection,

    # correlation/heatmap helpers
    plot_corr_distributions,      # will be slightly extended here
    plot_corr_violin,             # will be slightly extended here
    clustered_heatmap_and_dendrogram,
    clustered_heatmap_and_dendrogram_diff,
    cluster_labels_from_linkage,
    save_gene_cluster_labels,
    plot_heatmap_with_community_bar,

    # graph build + utilities (existing)
    build_weighted_graph_from_corr,
    build_diff_graph_from_corr,
    remove_isolates_in_both,
    compute_graph_metrics_df,
    plot_metric_distributions,

    # communities + layout
    detect_communities_on_graph,        # may NOT support resolution
    detect_communities_on_graph_v2,     # supports resolution + modularity + weight
    nx9_remove_degree0,
    nx9_largest_component,
    nx9_merge_communities_to_topK,
    nx9_safe_fill_pos,
    nx9_order_genes_by_community,
    nx9_restrict_and_order_matrix,
    community_aware_layout_from_diff_v3,

    # plotting networks
    plot_network_by_community_fixedpos,
    plot_global_network_fixedpos,
    plot_combined_network_fixedpos,


    build_diff_graph_topk,
    build_soft_wgcna_graph,
    build_diff_graph_soft_wgcna,
    plot_volcano,
    plot_manhattan_index,
    remove_isolated_nodes,
    save_graph_edgelist,

    clean_ensg,
    banner,
    shorten_entity_id,
    ensure_abs_weight,
    report_gene_presence,
    save_edges_both_formats,
    partition_to_labels,
    save_curve

)


VERBOSE = True

def log(msg: str):
    if VERBOSE:
        print(msg)


# =============================================================================
# STEP 01 — Input/Output paths
# =============================================================================
log("[STEP 01] Setting input paths...")
TCGA_TSV_PATH = (
    "/Users/mohadeshsoleimanour/TX Biomed Dropbox/Mohadeseh Soleimanpour-moghadam/"
    "A3-DECONVOLUTION-CANCER/APOBEC3-gene-expression-profiling-in-cancer_PAPER/"
    "Data-and-Analysis/TCGA_Gene_Expression_Data/TCGA_master_FPKM_UQ.tsv"
)
log(f"[STEP 01] TCGA TSV path: {TCGA_TSV_PATH}")
# signature_file = "InputData/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv"

signature_file = "/Users/mohadeshsoleimanour/TX Biomed Dropbox/Mohadeseh Soleimanpour-moghadam/A3-DECONVOLUTION-CANCER/APOBEC3-gene-expression-profiling-in-cancer_PAPER/Data-and-Analysis/Cosmic SBS Dataset/Mutation_Table_Tumors_TCGA_Mohadeseh.csv"

log(f"[STEP 01] signature_file  path: {signature_file}")

cancer_types = ["TCGA-HNSC"]  # modify as needed

RUN_ROOT = "@@Final_Network_community_code"
os.makedirs(RUN_ROOT, exist_ok=True)

# =============================================================================
# STEP 02 — Load TCGA gene expression TSV (raw)
# =============================================================================
log("[STEP 02] Loading TCGA gene expression data (raw)...")
df_raw = pd.read_csv(TCGA_TSV_PATH, sep="\t", header=None)
log(f"[STEP 02] Raw shape: {df_raw.shape}")
log("[STEP 02] Explanation: rows 0/1/2 contain gene annotations; rows 3+ are samples.")


# =============================================================================
# STEP 03 — Remove _PAR_Y gene columns safely
# =============================================================================
log("[STEP 03] Detecting _PAR_Y gene columns (if any)...")

ensg_row0 = df_raw.iloc[0, 5:].astype(str)  # row 0 holds ENSG IDs (often versioned)
par_y_mask = ensg_row0.str.contains("_PAR_Y", na=False)
cols_to_drop = ensg_row0[par_y_mask].index.tolist()

log(f"[STEP 03] Found {len(cols_to_drop)} _PAR_Y columns.")
log("[STEP 03] Explanation: _PAR_Y genes are pseudoautosomal Y-linked duplicates; commonly removed to avoid redundancy.")

# sanity check: are they all zeros in expression rows?
if len(cols_to_drop) > 0:
    expr_block = df_raw.iloc[3:, 5:].apply(pd.to_numeric, errors="coerce")
    par_y_sums = expr_block.loc[:, cols_to_drop].sum()
    if (par_y_sums != 0).any():
        log("[STEP 03] WARNING: Some _PAR_Y columns have non-zero sums (unexpected). Review before dropping.")
    else:
        log("[STEP 03] OK: All _PAR_Y columns sum to zero. Safe to drop.")

df = df_raw.drop(columns=cols_to_drop).copy()
log(f"[STEP 03] Shape after dropping _PAR_Y: {df.shape}")


# =============================================================================
# STEP 04 — Assign proper column headers (5 clinical + ENSG IDs from row 0)
# =============================================================================
log("[STEP 04] Assigning column headers...")

clinical_headers = df.iloc[0, :5].astype(str).tolist()
gene_headers_raw = df.iloc[0, 5:].astype(str).tolist()
df.columns = clinical_headers + gene_headers_raw

log(f"[STEP 04] Total columns assigned: {len(df.columns)}")
log(f"[STEP 04] First 10 columns: {df.columns[:10].tolist()}")
log("[STEP 04] Explanation: first 5 are meta columns, remaining are gene columns (currently as in the file).")


# =============================================================================
# STEP 05 — Extract gene annotations (row 0/1/2)
# =============================================================================
log("[STEP 05] Extracting gene annotations from rows 0/1/2...")

ensg_ids_raw = df.iloc[0, 5:].astype(str)     # ENSG with versions like ENSG...15
gene_symbols = df.iloc[1, 5:].astype(str)     # gene symbols like TP53, APOBEC3A, ...
biotype_flag = df.iloc[2, 5:].astype(str)     # protein_coding / biotype / flag

log(f"[STEP 05] Example ENSG IDs: {ensg_ids_raw.iloc[:5].tolist()}")
log(f"[STEP 05] Example Gene symbols: {gene_symbols.iloc[:5].tolist()}")
log(f"[STEP 05] Example biotype/flag: {biotype_flag.iloc[:5].tolist()}")
log("[STEP 05] Explanation: These rows describe the gene columns; samples start at row index 3.")


# =============================================================================
# STEP 06 — Clean ENSG IDs (remove .XX versions)
# =============================================================================
log("[STEP 06] Cleaning ENSG IDs (removing version suffix .XX)...")

def clean_ensg(x: str) -> str:
    # remove trailing ".digits" version
    return re.sub(r"\.\d+$", "", str(x))

ensg_ids_clean = ensg_ids_raw.map(clean_ensg)

log(f"[STEP 06] Total gene columns: {len(ensg_ids_raw)}")
log(f"[STEP 06] Unique base ENSG IDs: {ensg_ids_clean.nunique()}")
log("[STEP 06] Explanation: We drop versions because most mappings/lists use base ENSG IDs (no .XX).")


# =============================================================================
# STEP 06A — Detect and resolve duplicates created by version removal (IMPORTANT)
# =============================================================================
log("[STEP 06A] Checking for duplicated base ENSG IDs after cleaning...")

dup_counts = ensg_ids_clean.value_counts()
dup_bases = dup_counts[dup_counts > 1]

if len(dup_bases) == 0:
    log("[STEP 06A] No duplicates created by version removal. Great.")
else:
    log(f"[STEP 06A] WARNING: {len(dup_bases)} base ENSG IDs have multiple version-columns.")
    log("[STEP 06A] We will merge duplicates by SUM (other options: mean/max).")

    # We'll merge later on the expression matrix (STEP 08A) after converting to numeric.


# =============================================================================
# STEP 07 — Build expression matrix (rows 3+), keep meta separate
# =============================================================================
log("[STEP 07] Building meta + expression matrices...")

meta = df.iloc[3:, :5].copy()
expr = df.iloc[3:, 5:].copy().apply(pd.to_numeric, errors="coerce")

log(f"[STEP 07] meta shape: {meta.shape}")
log(f"[STEP 07] expr shape: {expr.shape}")
log("[STEP 07] Explanation: expr is numeric; meta stays as strings/IDs.")


# =============================================================================
# STEP 08 — Apply cleaned ENSG IDs as expr columns
# =============================================================================
log("[STEP 08] Renaming expression columns to cleaned ENSG IDs (no .XX)...")
expr.columns = ensg_ids_clean.tolist()
log(f"[STEP 08] Example expr columns: {expr.columns[:10].tolist()}")


# =============================================================================
# STEP 08A — Merge duplicate columns (only if needed)
# =============================================================================
if expr.columns.duplicated().any():
    log("[STEP 08A] Merging duplicated genes created by version removal (SUM across duplicates)...")

    # groupby over columns (axis=1) merges duplicates with same column name
    expr = expr.groupby(level=0, axis=1).sum(min_count=1)

    log(f"[STEP 08A] New expr shape after merge: {expr.shape}")
    log("[STEP 08A] Explanation: ENSG.1 + ENSG.2 now become one ENSG column.")
else:
    log("[STEP 08A] No duplicate columns to merge.")


# =============================================================================
# STEP 09 — Build mapping base ENSG -> gene symbol (for labeling)
# =============================================================================
log("[STEP 09] Building ENSG->GeneSymbol mapping (base ENSG)...")

# If duplicates exist, mapping can be ambiguous; we take the first symbol encountered.
ensg_to_symbol = {}
for e_clean, sym in zip(ensg_ids_clean.tolist(), gene_symbols.tolist()):
    if e_clean not in ensg_to_symbol:
        ensg_to_symbol[e_clean] = sym

log(f"[STEP 09] Example mapping items: {list(ensg_to_symbol.items())[:5]}")
log("[STEP 09] Explanation: mapping is used for plots/labels; analysis keys stay as ENSG IDs.")


# =============================================================================
# STEP 09A — Create sample-level expression table (meta + expr)
# =============================================================================
log("====================================================")
log("[STEP 09A] Build sample-level table: meta + expr (ready for merges)")
log("====================================================")

tcga_expr_df = pd.concat([meta.reset_index(drop=True), expr.reset_index(drop=True)], axis=1)
log(f"[STEP 09A] tcga_expr_df shape: {tcga_expr_df.shape}")
log(f"[STEP 09A] Columns (first 10): {tcga_expr_df.columns[:10].tolist()}")



# =============================================================================
# STEP 09A — Create sample-level expression table (meta + expr)
# =============================================================================
log("====================================================")
log("[STEP 09A] Sanity check:   [STEP 12A] Duplicated Entity_ID_short inside merged_TCGA: 102 ")
log("====================================================")


def tcga_sample_type(barcode: str) -> str:
    parts = str(barcode).split("-")
    if len(parts) < 4:
        return ""
    return parts[3][:2]   # "01", "11", etc.

tcga_expr_df["sample_type"] = tcga_expr_df["Entity_ID"].apply(tcga_sample_type)

print(tcga_expr_df["sample_type"].value_counts(dropna=False).head(20))
print("Normals (11) rows:", (tcga_expr_df["sample_type"] == "11").sum())
print("Tumors (01) rows:", (tcga_expr_df["sample_type"] == "01").sum())

tcga_expr_df = tcga_expr_df[tcga_expr_df["sample_type"] == "01"].copy()
print("After tumor-only filter:", tcga_expr_df.shape)

# ------------------------------------------------------------------
# Visibility helper
# ------------------------------------------------------------------
def banner(title: str, char: str = "=", width: int = 100):
    line = char * width
    print("\n" + line)
    print(title)
    print(line)


# =============================================================================
# STEP 10 — Load SBS signature file + merge with TCGA expression (Gene IDs)
# =============================================================================
banner("[STEP 10] Load SBS signature file + merge with TCGA expression (Gene IDs)")

log("[STEP 10] Loading signature file...")
signature_df = pd.read_csv(signature_file)
log(f"[STEP 10] signature_df shape: {signature_df.shape}")

def shorten_entity_id(x):
    """Shorten to first 4 parts: TCGA-XX-YYYY-ZZZZ (matches many TCGA tables)."""
    return "-".join(str(x).split("-")[:4])

# IMPORTANT: use tcga_expr_df (not df_TCGA_New)
tcga_merged_base = tcga_expr_df.copy()

if "Entity_ID" not in tcga_merged_base.columns:
    raise ValueError("tcga_expr_df must contain 'Entity_ID' in the first 5 clinical columns.")

tcga_merged_base["Entity_ID_short"] = tcga_merged_base["Entity_ID"].apply(shorten_entity_id)

if "Sample Names" not in signature_df.columns:
    raise ValueError("signature_df must contain a column named 'Sample Names'.")

signature_df["Entity_ID_short"] = signature_df["Sample Names"].apply(shorten_entity_id)

merged_TCGA = pd.merge(
    tcga_merged_base,
    signature_df,
    on="Entity_ID_short",
    how="inner"
)



# =============================================================================
# REMOVE duplicated Entity_ID_short (drop BOTH rows)
# =============================================================================
banner("[STEP 10B] Remove duplicated samples (drop both rows)")

"""
WHY THIS STEP EXISTS
--------------------

During STEP 10 we merge TCGA expression data with the SBS signature table
using a shortened TCGA barcode (Entity_ID_short).

In TCGA datasets it is common that the same short barcode appears multiple times
because of:
    • multiple aliquots of the same tumor sample
    • multiple WES runs for the same sample
    • repeated rows in the signature table
    • metadata differences while referring to the same biological sample

When such rows exist on BOTH sides of the merge, pandas performs a
many-to-many merge and duplicates the resulting rows.

Example:
    expression table contains:
        TCGA-XX-1234-01A

    signature table contains:
        TCGA-XX-1234-01A
        TCGA-XX-1234-01A

The merge produces TWO rows for the same sample.

In our dataset we detected:
    51 duplicated samples → 102 duplicated rows.

To avoid:
    • double counting samples
    • biasing statistical tests
    • inflating correlations or network edges

we REMOVE ALL duplicated keys entirely rather than keeping one.

This ensures that only samples with a clean one-to-one mapping between
expression and mutation signature remain in the dataset.

If future analyses need these samples, a more sophisticated strategy could be:
    • selecting the best aliquot
    • averaging replicates
    • merging by patient-level barcode

For the current analysis we prefer strict one-to-one sample mapping.
"""


before_rows = merged_TCGA.shape[0]

# identify duplicated keys
dup_mask = merged_TCGA["Entity_ID_short"].duplicated(keep=False)

dup_keys = merged_TCGA.loc[dup_mask, "Entity_ID_short"].unique()
log(f"[STEP 10B] Duplicated keys detected: {len(dup_keys)}")

# remove all duplicated rows
merged_TCGA = merged_TCGA.loc[~dup_mask].copy()

after_rows = merged_TCGA.shape[0]

log(f"[STEP 10B] Rows before: {before_rows}")
log(f"[STEP 10B] Rows after removing duplicates: {after_rows}")
log(f"[STEP 10B] Rows removed: {before_rows - after_rows}")



log(f"[STEP 10] merged_TCGA shape: {merged_TCGA.shape}")
log(f"[STEP 10] Example merged Entity_ID_short: {merged_TCGA['Entity_ID_short'].head().tolist()}")

log("====================================================")
log("@@@[STEP 09A] Sanity check:   [STEP 12A] Duplicated Entity_ID_short inside merged_TCGA: 102 ")
log("====================================================")
# keys duplicated in merged_TCGA (based on the merge key you used)
dup_keys = (
    merged_TCGA.loc[merged_TCGA["Entity_ID_short"].duplicated(keep=False), "Entity_ID_short"]
    .astype(str)
    .unique()
)

print("Number of duplicated Entity_ID_short keys:", len(dup_keys))
print("First 20 duplicated keys:", dup_keys[:20])

dup_rows = merged_TCGA[merged_TCGA["Entity_ID_short"].isin(dup_keys)].copy()

# Keep a focused set of columns that show the identity on BOTH sides
cols_wanted = []

# expression full barcode (your table)
if "Entity_ID" in dup_rows.columns:
    cols_wanted.append("Entity_ID")

# signature full barcode column
if "Sample Names" in dup_rows.columns:
    cols_wanted.append("Sample Names")

# also keep merge key for grouping
cols_wanted.append("Entity_ID_short")

# plus a few helpful metadata columns if present
for c in ["Project_ID", "Tissue_Type", "Case_ID", "File_ID"]:
    if c in dup_rows.columns:
        cols_wanted.append(c)

dup_rows_view = dup_rows[cols_wanted].copy()

print("Rows involved in duplicated keys:", dup_rows_view.shape[0])
print(dup_rows_view.head(30))


dup_rows_view["dup_count_for_key"] = dup_rows_view.groupby("Entity_ID_short")["Entity_ID_short"].transform("size")

dup_rows_view = dup_rows_view.sort_values(
    ["dup_count_for_key", "Entity_ID_short", "Entity_ID", "Sample Names"],
    ascending=[False, True, True, True]
)

print("\nTop duplicated keys by multiplicity:")
print(dup_rows_view[["Entity_ID_short", "dup_count_for_key"]].drop_duplicates().head(20))

out_csv = os.path.join(RUN_ROOT, "duplicated_merge_rows_FULL_TCGA_IDS.csv")
dup_rows_view.to_csv(out_csv, index=False)
print("Saved:", out_csv)


# =============================================================================
# SANITY CHECK — Correlation of duplicated samples
# =============================================================================
banner("[SANITY] Correlation check for duplicated samples")

sanity_dir = os.path.join(RUN_ROOT, "sanity_check_duplicates")
os.makedirs(sanity_dir, exist_ok=True)

# gene columns
gene_cols = [c for c in merged_TCGA.columns if isinstance(c, str) and c.startswith("ENSG")]

# find duplicated keys
dup_keys = merged_TCGA.loc[
    merged_TCGA["Entity_ID_short"].duplicated(keep=False),
    "Entity_ID_short"
].unique()

log(f"[SANITY] Duplicated sample keys: {len(dup_keys)}")

corr_rows = []

for key in dup_keys:

    subset = merged_TCGA[merged_TCGA["Entity_ID_short"] == key]

    # should be exactly 2 rows
    if subset.shape[0] != 2:
        log(f"[SANITY] WARNING: {key} has {subset.shape[0]} rows (expected 2)")
        continue

    row1 = subset.iloc[0]
    row2 = subset.iloc[1]

    expr1 = pd.to_numeric(row1[gene_cols], errors="coerce")
    expr2 = pd.to_numeric(row2[gene_cols], errors="coerce")

    corr = expr1.corr(expr2, method="spearman")

    corr_rows.append({
        "Entity_ID_short": key,
        "Entity_ID_1": row1.get("Entity_ID", ""),
        "Entity_ID_2": row2.get("Entity_ID", ""),
        "corr_spearman": corr
    })

# convert to dataframe
corr_df = pd.DataFrame(corr_rows)

# save results
corr_path = os.path.join(sanity_dir, "duplicate_sample_correlations.csv")
corr_df.to_csv(corr_path, index=False)

log(f"[SANITY] Saved correlations -> {corr_path}")




log("====================================================")
log("====================================================")

# =============================================================================
# STEP 10A — Add A3* alias columns using ENSG IDs (do NOT replace gene IDs)
# =============================================================================
banner("[STEP 10A] Create A3 alias columns (A3A/A3B/...) from ENSG IDs")

# ENSG IDs (CLEANED, no .XX) you want to keep all the way
A3_GENES = [
    "ENSG00000128383",  # A3A
    "ENSG00000179750",  # A3B
    "ENSG00000244509",  # A3C
    "ENSG00000243811",  # A3D
    "ENSG00000128394",  # A3F
    "ENSG00000239713",  # A3G
    "ENSG00000100298",  # A3H
]

# Biomarkers: DO NOT force keep; only check if they were tested/selected
BIOMARKERS = [

    # CD4-Positive alpha-beta T cell
    "ENSG00000168685",  # IL7R
    "ENSG00000198851",  # CD3E

    # B Cells
    "ENSG00000007314",  # CD79A
    "ENSG00000156738",  # MS4A1 (CD20)

    # fibroblasts
    "ENSG00000011465",  # DCN
    "ENSG00000108821",  # COL1A1

    # CD8-Positive alpha-beta T cell
    "ENSG00000153563",  # CD8A
    "ENSG00000271503",  # CCL5

    # basal cell
    "ENSG00000186081",  # KRT5
    "ENSG00000166886",  # TACSTD2 (EpCAM)

    # macrophage
    "ENSG00000129226",  # CD68
    "ENSG00000066336",  # SPI1 (PU.1)

    # endothelial cell
    "ENSG00000261371",  # PECAM1
    "ENSG00000110799",  # VWF

    # smooth muscle cell
    "ENSG00000149591",  # TAGLN
    "ENSG00000143248",  # RGS5

    # regulatory T cell
    "ENSG00000181847",  # TIGIT
    "ENSG00000163599",  # CTLA4

    # mast cell
    "ENSG00000172236",  # TPSAB1
    "ENSG00000163736",  # CPA3

    # plasmacytoid dendritic cell
    "ENSG00000185291",  # IL3RA
    "ENSG00000162692",  # LILRA4

    # myeloid dendritic cell
    "ENSG00000110395",  # LAMP3
    "ENSG00000126353",  # CCR7
]
#
# BIOMARKERS = [
# "IL7R","CD3E",
# "CD79A","MS4A1",
# "DCN","COL1A1",
# "CD8A","CCL5",
# "KRT5","TACSTD2",
# "CD68","SPI1",
# "PECAM1","VWF",
# "TAGLN","RGS5",
# "TIGIT","CTLA4",
# "TPSAB1","CPA3",
# "IL3RA","LILRA4",
# "LAMP3","CCR7"
# ]


# Convenience alias columns (copies) for scoring/plots
A3_ID_TO_ALIAS = {
    "ENSG00000128383": "A3A",
    "ENSG00000179750": "A3B",
    "ENSG00000244509": "A3C",
    "ENSG00000243811": "A3D",
    "ENSG00000128394": "A3F",
    "ENSG00000239713": "A3G",
    "ENSG00000100298": "A3H",
}

for ensg_id, alias in A3_ID_TO_ALIAS.items():
    if ensg_id in merged_TCGA.columns:
        merged_TCGA[alias] = pd.to_numeric(merged_TCGA[ensg_id], errors="coerce")
    else:
        log(f"[STEP 10A] WARNING: {ensg_id} not found in merged_TCGA (cannot create {alias}).")

log("[STEP 10A] Explanation: analysis stays on ENSG IDs; A3A/A3B columns are convenience copies for plotting/scoring.")


# =============================================================================
# STEP 11 — Cancer types list from Project_ID + save
# =============================================================================
banner("[STEP 11] Extract cancer types from Project_ID + save")

if "Project_ID" not in merged_TCGA.columns:
    raise ValueError("merged_TCGA must contain 'Project_ID'.")

cancer_types_all = sorted(merged_TCGA["Project_ID"].dropna().unique().tolist())
log(f"[STEP 11] Total cancer types found: {len(cancer_types_all)}")

os.makedirs(RUN_ROOT, exist_ok=True)
pd.Series(cancer_types_all).to_csv(
    f"{RUN_ROOT}/cancer_types_from_TCGA.txt",
    index=False,
    header=False
)
log(f"[STEP 11] Saved cancer types -> {RUN_ROOT}/cancer_types_from_TCGA.txt")

# =============================================================================
# STEP 12 — Loop per cancer type + t-test (Gene IDs) + force-keep A3 ONLY
# =============================================================================
banner("[STEP 12] Per-cancer t-test (High SBS2 vs Low SBS2) using ENSG IDs")

P_THR = 0.025  # raw p-value selection threshold

# all ENSG gene columns in merged_TCGA
all_gene_cols = [c for c in merged_TCGA.columns if isinstance(c, str) and c.startswith("ENSG")]
log(f"[STEP 12] Total ENSG gene columns in merged_TCGA: {len(all_gene_cols)}")

# sanity: duplicated column names?
if pd.Index(merged_TCGA.columns).duplicated().any():
    dup_cols = pd.Index(merged_TCGA.columns)[pd.Index(merged_TCGA.columns).duplicated()].unique().tolist()
    log(f"[STEP 12] WARNING: duplicated column names exist: {dup_cols[:10]} (showing first 10)")
else:
    log("[STEP 12] OK: No duplicated column names in merged_TCGA.")

# A3 genes / biomarkers availability
present_a3 = [g for g in A3_GENES if g in merged_TCGA.columns]
missing_a3 = [g for g in A3_GENES if g not in merged_TCGA.columns]
log(f"[STEP 12] A3 genes present: {len(present_a3)} | missing: {missing_a3}")

present_bm = [g for g in BIOMARKERS if g in merged_TCGA.columns]
missing_bm = [g for g in BIOMARKERS if g not in merged_TCGA.columns]
log(f"[STEP 12] Biomarkers present: {len(present_bm)} | missing: {missing_bm}")

# =============================================================================
# STEP 12A — Merge-key diagnostics (helps catch silent merge issues)
# =============================================================================
banner("[STEP 12A] Merge-key diagnostics (Entity_ID_short)")

if "Entity_ID_short" in merged_TCGA.columns:
    dup_key_total = int(merged_TCGA["Entity_ID_short"].duplicated().sum())
    log(f"[STEP 12A] Duplicated Entity_ID_short inside merged_TCGA: {dup_key_total}")
else:
    log("[STEP 12A] NOTE: merged_TCGA has no Entity_ID_short column (unexpected).")

# =============================================================================
# Main loop
# =============================================================================
for cancer_type in cancer_types:

    banner(f"[CANCER] {cancer_type}", char="=")

    cancer_dir = os.path.join(RUN_ROOT, cancer_type)
    os.makedirs(cancer_dir, exist_ok=True)

    # -------------------------
    # STEP 12.1 — Filter by cancer type
    # -------------------------
    log("[STEP 12.1] Filtering merged_TCGA by cancer type...")
    cancer_df = merged_TCGA[merged_TCGA["Project_ID"] == cancer_type].copy()
    log(f"[STEP 12.1] cancer_df shape: {cancer_df.shape}")

    if len(cancer_df) < 25:
        log("[SKIP] Too few samples for this cancer.")
        continue

    # -------------------------
    # STEP 12.2 — Prepare SBS2 + define groups
    # -------------------------
    if "SBS2" not in cancer_df.columns:
        log("[SKIP] SBS2 column missing.")
        continue

    tmp = cancer_df.copy()
    tmp["SBS2"] = pd.to_numeric(tmp["SBS2"], errors="coerce")
    tmp = tmp.dropna(subset=["SBS2"]).reset_index(drop=True)

    if len(tmp) < 40:
        log("[SKIP] Too few valid SBS2 samples.")
        continue

    high_thr = tmp["SBS2"].quantile(0.80)
    low_thr  = tmp["SBS2"].quantile(0.20)

    high_df = tmp[tmp["SBS2"] >= high_thr].copy()
    low_df  = tmp[tmp["SBS2"] <= low_thr].copy()

    log(f"[STEP 12.2] High group: {len(high_df)} | Low group: {len(low_df)}")
    if len(high_df) < 10 or len(low_df) < 10:
        log("[SKIP] Too few samples per group.")
        continue

    # -------------------------
    # STEP 12.3 — Candidate genes (ENSG IDs)
    # -------------------------
    # start from all ENSG columns present
    candidate_genes = [g for g in all_gene_cols if g in tmp.columns]

    # ensure A3 + biomarkers are at least included in TESTING if present
    candidate_genes = list(dict.fromkeys(
        candidate_genes
        + [g for g in A3_GENES if g in tmp.columns]
        + [g for g in BIOMARKERS if g in tmp.columns]
    ))

    log(f"[STEP 12.3] Candidate genes for t-test: {len(candidate_genes)}")
    log(f"[STEP 12.3] A3 present in tmp: {[g for g in A3_GENES if g in tmp.columns]}")
    log(f"[STEP 12.3] Biomarkers present in tmp: {[g for g in BIOMARKERS if g in tmp.columns]}")

    high_num = high_df[candidate_genes].apply(pd.to_numeric, errors="coerce")
    low_num  = low_df[candidate_genes].apply(pd.to_numeric, errors="coerce")

    # -------------------------
    # STEP 12.4 — Welch t-tests
    # -------------------------
    log("[STEP 12.4] Running Welch t-tests...")

    pvals, log2fc, mean_high, mean_low, tested = [], [], [], [], []
    eps = 1e-9

    for g in candidate_genes:
        x = high_num[g].dropna()
        y = low_num[g].dropna()
        if len(x) < 8 or len(y) < 8:
            continue

        _, p = ttest_ind(x, y, equal_var=False, nan_policy="omit")
        mh = float(x.mean())
        ml = float(y.mean())
        fc = float(np.log2((mh + eps) / (ml + eps)))

        tested.append(g)
        pvals.append(p)
        log2fc.append(fc)
        mean_high.append(mh)
        mean_low.append(ml)

    if len(tested) == 0:
        log("[SKIP] No genes tested (too many NaNs or small groups).")
        continue

    pvals = np.array(pvals)
    qvals = bh_fdr(pvals)

    stats_df = pd.DataFrame({
        "gene": tested,  # ENSG IDs
        "pvalue": pvals,
        "qvalue": qvals,
        "mean_high": mean_high,
        "mean_low": mean_low,
        "log2FC_high_vs_low": log2fc
    }).sort_values("pvalue")

    # -------------------------
    # STEP 12.5 — Save stats
    # -------------------------
    diff_dir = os.path.join(cancer_dir, "12_diff_SBS2")
    os.makedirs(diff_dir, exist_ok=True)

    stats_path = os.path.join(diff_dir, f"{cancer_type}_SBS2_ttest_stats.csv")
    stats_df.to_csv(stats_path, index=False)
    log(f"[STEP 12.5] Saved stats -> {stats_path}")

    # -------------------------
    # STEP 12.6 — Select genes by p-value + FORCE keep A3 ONLY
    # -------------------------
    selected_genes = stats_df.loc[stats_df["pvalue"] <= P_THR, "gene"].tolist()

    # FORCE keep A3 genes only
    for g in A3_GENES:
        if g in tmp.columns and g not in selected_genes:
            selected_genes.insert(0, g)

    selected_genes = list(dict.fromkeys(selected_genes))
    log(f"[STEP 12.6] Selected genes (p ≤ {P_THR} + forced A3): {len(selected_genes)}")

    # Report A3 status
    log("[STEP 12.6] A3 genes kept:")
    for g in A3_GENES:
        log(f"  - {g}: {'YES' if g in selected_genes else 'NO'}")

    # Report biomarker status: tested? selected?
    tested_set = set(stats_df["gene"])
    log("[STEP 12.6] Biomarkers tested + selected (NOT forced):")
    for g in BIOMARKERS:
        tested_flag = "TESTED" if g in tested_set else "NOT_TESTED"
        selected_flag = "SELECTED" if g in selected_genes else "NOT_SELECTED"
        log(f"  - {g}: {tested_flag} / {selected_flag}")

    # -------------------------
    # STEP 12.7 — Build cancer_df2 (clinical + selected genes + SBS2/A3A/A3B)
    # -------------------------
    needed_columns = ["Project_ID", "Tissue_Type", "Case_ID", "File_ID", "Entity_ID"]

    # keep the signal columns needed for Step 13/14
    extra_signal_cols = ["SBS2", "A3A", "A3B"]
    needed_columns += [c for c in extra_signal_cols if c in tmp.columns]

    # add selected ENSG genes
    needed_columns += selected_genes

    # drop missing (safety)
    missing_needed = [c for c in needed_columns if c not in tmp.columns]
    if missing_needed:
        log(f"[STEP 12.7] WARNING: missing columns in tmp (skipping): {missing_needed[:10]}")
    needed_columns = [c for c in needed_columns if c in tmp.columns]

    cancer_df2 = tmp[needed_columns].copy()
    log(f"[STEP 12.7] cancer_df2 shape: {cancer_df2.shape}")
    # =============================================================================
    # STEP 12.8 — Final sanity checks for cancer_df2
    # =============================================================================
    banner("[STEP 12.8] Sanity check: cancer_df2 integrity")

    print("Shape:", cancer_df2.shape)
    print("First 10 columns:", cancer_df2.columns[:10].tolist())

    expected_clinical = ["Project_ID", "Tissue_Type", "Case_ID", "File_ID", "Entity_ID"]
    missing_clinical = [c for c in expected_clinical if c not in cancer_df2.columns]

    print("\nClinical columns present:", [c for c in expected_clinical if c in cancer_df2.columns])
    if missing_clinical:
        print("⚠ Missing clinical columns:", missing_clinical)
    else:
        print("✓ All expected clinical columns present")

    print("\nA3 genes in cancer_df2:")
    for g in A3_GENES:
        print(f"  {g}: {'YES' if g in cancer_df2.columns else 'NO'}")

    print("\nBiomarkers in cancer_df2:")
    for g in BIOMARKERS:
        print(f"  {g}: {'YES' if g in cancer_df2.columns else 'NO'}")

    dup_cols = cancer_df2.columns[cancer_df2.columns.duplicated()].tolist()
    if dup_cols:
        print("\n⚠ Duplicate columns detected:", dup_cols)
    else:
        print("\n✓ No duplicate column names")

    versioned_cols = [c for c in cancer_df2.columns if isinstance(c, str) and "." in c and c.startswith("ENSG")]
    if versioned_cols:
        print("\n⚠ Versioned ENSG IDs detected (unexpected):", versioned_cols[:10])
    else:
        print("\n✓ All ENSG IDs are cleaned (no .XX suffix)")

    gene_only_cols = [c for c in cancer_df2.columns if isinstance(c, str) and c.startswith("ENSG")]
    non_numeric = [c for c in gene_only_cols if not pd.api.types.is_numeric_dtype(cancer_df2[c])]

    if non_numeric:
        print("\n⚠ Non-numeric gene columns detected:", non_numeric[:10])
    else:
        print("\n✓ All gene expression columns are numeric")

    print("\nPreview of expression matrix (first 3 rows, first 5 genes):")
    if len(gene_only_cols) > 0:
        print(cancer_df2[gene_only_cols[:5]].head(3))
    else:
        print("No ENSG gene columns found.")

    # =============================================================================
    # STEP 12.6A — Quick report for plotting
    # =============================================================================
    banner("[STEP 12.6A] Plot report (ttest)")

    print("Top 10 genes by p-value:")
    print(stats_df.head(10)[["gene", "pvalue", "log2FC_high_vs_low"]])

    # show A3 + biomarkers p-values if they were tested
    tested_set = set(stats_df["gene"])
    print("\nA3 stats (if tested):")
    for g in A3_GENES:
        if g in tested_set:
            row = stats_df.loc[stats_df["gene"] == g].iloc[0]
            print(f"  {g} | p={row['pvalue']:.3e} | log2FC={row['log2FC_high_vs_low']:.3f}")
        else:
            print(f"  {g} | NOT_TESTED")

    print("\nBiomarker stats (if tested):")
    for g in BIOMARKERS:
        if g in tested_set:
            row = stats_df.loc[stats_df["gene"] == g].iloc[0]
            print(f"  {g} | p={row['pvalue']:.3e} | log2FC={row['log2FC_high_vs_low']:.3f}")
        else:
            print(f"  {g} | NOT_TESTED")

    # =============================================================================
    # STEP 12.6B — Volcano Plot (with p-value threshold line)
    # =============================================================================
    banner("[STEP 12.6B] Volcano plot")

    y = -np.log10(stats_df["pvalue"].values + 1e-300)
    x = stats_df["log2FC_high_vs_low"].values

    plt.figure(figsize=(7, 6))
    plt.scatter(x, y, s=8, alpha=0.6)

    # threshold line
    plt.axhline(-np.log10(P_THR), linestyle="--", color="red", linewidth=1.5)

    # optional: annotate A3 + biomarkers (only if present in stats_df)
    annotate_genes = list(dict.fromkeys(A3_GENES + BIOMARKERS))
    for g in annotate_genes:
        hit = stats_df[stats_df["gene"] == g]
        if len(hit) == 1:
            gx = float(hit["log2FC_high_vs_low"].iloc[0])
            gy = float(-np.log10(hit["pvalue"].iloc[0] + 1e-300))

            # label: gene symbol if available, else ENSG
            label = ensg_to_symbol.get(g, g) if "ensg_to_symbol" in globals() else g
            plt.text(gx, gy, label, fontsize=8)

    plt.xlabel("log2FC (High SBS2 / Low SBS2)")
    plt.ylabel("-log10(p-value)")
    plt.title(f"{cancer_type} | Volcano (High vs Low SBS2)")
    plt.tight_layout()

    volcano_path = os.path.join(diff_dir, f"{cancer_type}_volcano_p{P_THR}.png")
    plt.savefig(volcano_path, dpi=300)
    plt.close()
    log(f"[STEP 12.6B] Saved volcano plot -> {volcano_path}")

    # =============================================================================
    # STEP 12.6C — Manhattan-like Plot (sorted by p-value)
    # =============================================================================
    banner("[STEP 12.6C] Manhattan-like plot")

    stats_sorted = stats_df.sort_values("pvalue").reset_index(drop=True)

    x = np.arange(len(stats_sorted))
    y = -np.log10(stats_sorted["pvalue"].values + 1e-300)

    plt.figure(figsize=(14, 4))
    plt.scatter(x, y, s=6, alpha=0.6)

    plt.axhline(-np.log10(P_THR), linestyle="--", color="red", linewidth=1.5)

    # optional: mark A3 + biomarkers as vertical ticks
    for g in annotate_genes:
        idx = stats_sorted.index[stats_sorted["gene"] == g].tolist()
        if idx:
            i = idx[0]
            plt.scatter([i], [y[i]], s=40, marker="x")

    plt.xlabel("Genes (sorted by p-value)")
    plt.ylabel("-log10(p-value)")
    plt.title(f"{cancer_type} | Manhattan-like (High vs Low SBS2)")
    plt.tight_layout()

    manhattan_path = os.path.join(diff_dir, f"{cancer_type}_manhattan_p{P_THR}.png")
    plt.savefig(manhattan_path, dpi=300)
    plt.close()
    log(f"[STEP 12.6C] Saved Manhattan plot -> {manhattan_path}")




    print("\n✓ STEP 12 COMPLETE for", cancer_type)

    # =============================================================================
    # STEP 13 — Define TOP and BOTTOM groups using APOBEC score (A3_score)
    # =============================================================================
    banner("[STEP 13] Define TOP/BOTTOM groups using A3_score (High-A3 then SBS2 split)")

    step13_dir = os.path.join(cancer_dir, "13_TOP_BOTTOM_by_A3score")
    os.makedirs(step13_dir, exist_ok=True)

    needed_cols = ["SBS2", "A3A", "A3B"]
    missing_needed_cols = [c for c in needed_cols if c not in cancer_df2.columns]
    if missing_needed_cols:
        log(f"[SKIP][STEP 13] Missing required columns in cancer_df2: {missing_needed_cols}")
        continue

    # numeric conversion
    df_step13 = cancer_df2.copy()
    for c in needed_cols:
        df_step13[c] = pd.to_numeric(df_step13[c], errors="coerce")

    # keep stable index (no reset_index)
    df_step13 = df_step13.dropna(subset=needed_cols)
    log(f"[STEP 13] After SBS2/A3A/A3B cleanup: {df_step13.shape}")

    if len(df_step13) < 30:
        log("[SKIP][STEP 13] Too few samples after cleanup (<30).")
        continue

    # compute A3_score
    df_step13["A3_score"] = compute_a3_score_AB(df_step13, "A3A", "A3B")

    # QC plots
    qc_dir = os.path.join(step13_dir, "13_1_QC_A3score")
    os.makedirs(qc_dir, exist_ok=True)
    plot_a3_scatter_and_distribution(df_step13, qc_dir, cancer_type)
    log(f"[STEP 13.1] Saved A3_score QC plots -> {qc_dir}")

    # High-A3 restriction by median
    a3_thr = float(df_step13["A3_score"].median())
    high_a3_df = df_step13[df_step13["A3_score"] >= a3_thr].copy()

    log(f"[STEP 13.2] A3_score median threshold = {a3_thr:.6f}")
    log(f"[STEP 13.2] High-A3 samples: {len(high_a3_df)} / {len(df_step13)}")

    if len(high_a3_df) < 20:
        log("[SKIP][STEP 13] Too few High-A3 samples (<20).")
        continue

    # Split SBS2 within High-A3: TOP (top 20%), BOTTOM (bottom 20%)
    sbs_hi = float(high_a3_df["SBS2"].quantile(0.80))
    sbs_lo = float(high_a3_df["SBS2"].quantile(0.20))

    group_top = high_a3_df[high_a3_df["SBS2"] >= sbs_hi].copy()
    group_bot = high_a3_df[high_a3_df["SBS2"] <= sbs_lo].copy()

    log(f"[STEP 13.3] SBS2 thresholds within High-A3: low20={sbs_lo:.6f}, high80={sbs_hi:.6f}")
    log(f"[STEP 13.3] TOP size={len(group_top)} | BOTTOM size={len(group_bot)}")

    # Selection plot SBS2 vs A3_score
    selection_plot = os.path.join(step13_dir, f"{cancer_type}_SBS2_vs_A3score_TOP_BOTTOM.png")
    plot_sbs_vs_a3_selection(
        df=df_step13,
        a3_thr=a3_thr,
        sbs_low_thr=sbs_lo,
        sbs_high_thr=sbs_hi,
        group1_idx=group_top.index,
        group2_idx=group_bot.index,
        out_path=selection_plot,
        title=f"{cancer_type} | SBS2 vs A3_score (TOP/BOTTOM inside High-A3)"
    )
    log(f"[STEP 13.4] Saved selection plot -> {selection_plot}")

    # Save groups
    groups_dir = os.path.join(step13_dir, "13_2_groups_csv")
    os.makedirs(groups_dir, exist_ok=True)

    top_path = os.path.join(groups_dir, f"{cancer_type}_TOP_highSBS2_top20_inHighA3.csv")
    bot_path = os.path.join(groups_dir, f"{cancer_type}_BOTTOM_lowSBS2_bottom20_inHighA3.csv")
    group_top.to_csv(top_path, index=False)
    group_bot.to_csv(bot_path, index=False)
    log(f"[STEP 13.5] Saved TOP group -> {top_path}")
    log(f"[STEP 13.5] Saved BOTTOM group -> {bot_path}")

    MIN_GROUP = 8
    if len(group_top) < MIN_GROUP or len(group_bot) < MIN_GROUP:
        log(f"[SKIP][STEP 13] Groups too small for network analysis (need >= {MIN_GROUP}).")
        continue

    top_samples = group_top
    bottom_samples = group_bot

    print("\n[STEP 13 SUMMARY]")
    print(f"  Clean samples: {len(df_step13)}")
    print(f"  High-A3 samples: {len(high_a3_df)}")
    print(f"  TOP samples: {len(top_samples)}")
    print(f"  BOTTOM samples: {len(bottom_samples)}")
    # =============================================================================
    # STEP 14 — Correlations and networks (TOP, BOTTOM, DIFF)
    # =============================================================================
    banner("[STEP 14] Correlations + Networks (TOP/BOTTOM/DIFF)")

    step14_dir = os.path.join(cancer_dir, "14_Corr_Networks")
    os.makedirs(step14_dir, exist_ok=True)

    corr_dir = os.path.join(step14_dir, "14_1_corr_matrices_and_qc")
    heat_dir = os.path.join(step14_dir, "14_2_heatmaps_dendrograms")
    net_dir = os.path.join(step14_dir, "14_3_network_edge_lists")
    plot_dir = os.path.join(step14_dir, "14_4_network_plots")
    for d in [corr_dir, heat_dir, net_dir, plot_dir]:
        os.makedirs(d, exist_ok=True)

    # -----------------------------------------------------------------------------
    # STEP 14.0 — Decide which genes go into correlations (from Step 12 selection)
    # -----------------------------------------------------------------------------
    CORR_THRESHOLD = 0.80
    DIFF_THRESHOLD = 0.40

    selected_genes_corr = [g for g in selected_genes if g in top_samples.columns and g in bottom_samples.columns]
    selected_genes_corr = list(dict.fromkeys(selected_genes_corr))
    log(f"[STEP 14.0] Candidate genes (from Step 12, present in both groups): {len(selected_genes_corr)}")

    if len(selected_genes_corr) < 10:
        log("[SKIP][STEP 14] Too few genes for correlation network (<10).")
        continue

    # numeric matrices (samples x genes)
    top_num = top_samples[selected_genes_corr].apply(pd.to_numeric, errors="coerce")
    bot_num = bottom_samples[selected_genes_corr].apply(pd.to_numeric, errors="coerce")

    # remove genes constant in BOTH groups only (optional but recommended)
    const_top = top_num.nunique(dropna=True) <= 1
    const_bot = bot_num.nunique(dropna=True) <= 1
    const_both = const_top & const_bot

    n_const_both = int(const_both.sum())
    if n_const_both > 0:
        removed = list(top_num.columns[const_both])
        log(f"[STEP 14.0] Removing genes constant in BOTH groups: {n_const_both}")
        top_num = top_num.drop(columns=removed)
        bot_num = bot_num.drop(columns=removed)

    selected_genes_corr = list(top_num.columns)
    log(f"[STEP 14.0] Genes after constant-in-both filter: {len(selected_genes_corr)}")

    if len(selected_genes_corr) < 10:
        log("[SKIP][STEP 14] Too few genes after filtering (<10).")
        continue

    # -----------------------------------------------------------------------------
    # STEP 14.1 — Correlation matrices (Spearman) + QC
    # -----------------------------------------------------------------------------
    log(f"[STEP 14.1] Computing Spearman correlations (thr={CORR_THRESHOLD}) ...")

    corr_top_full = top_num.corr(method="spearman")
    corr_bot_full = bot_num.corr(method="spearman")

    # IMPORTANT: leave NaNs as NaN for now; thresholding will ignore them safely if the builder checks them,
    # otherwise fill with 0 only AFTER you are sure you want that behavior.
    corr_top_full = corr_top_full.fillna(0.0)
    corr_bot_full = corr_bot_full.fillna(0.0)

    corr_diff_full = (corr_top_full - corr_bot_full).fillna(0.0)  # signed diff

    # Save matrices (full, before AND restriction)
    corr_top_path = os.path.join(corr_dir, f"{cancer_type}_corr_TOP_spearman_FULL.csv")
    corr_bot_path = os.path.join(corr_dir, f"{cancer_type}_corr_BOTTOM_spearman_FULL.csv")
    corr_diff_path = os.path.join(corr_dir, f"{cancer_type}_corr_DIFF_signed_FULL.csv")
    corr_top_full.to_csv(corr_top_path)
    corr_bot_full.to_csv(corr_bot_path)
    corr_diff_full.to_csv(corr_diff_path)

    log(f"[STEP 14.1] Saved FULL correlation matrices -> {corr_dir}")

    # QC plots (optional)
    plot_corr_distributions(
        corr_top=corr_top_full,
        corr_bottom=corr_bot_full,
        out_dir=corr_dir,
        prefix=f"{cancer_type}_corrQC_FULL",
        threshold=CORR_THRESHOLD,
        seed=42
    )
    log(f"[STEP 14.1] Saved correlation QC plots -> {corr_dir}")

    print("\n[STEP 14.1 SUMMARY]")
    print("  corr_top_full shape:", corr_top_full.shape)
    print("  corr_bot_full shape:", corr_bot_full.shape)
    print("  corr_diff_full shape:", corr_diff_full.shape)

    # -----------------------------------------------------------------------------
    # STEP 14.2 — Build TOP/BOTTOM graphs, remove isolates, then compute AND gene set
    # -----------------------------------------------------------------------------
    banner("[STEP 14.2] Build graphs TOP/BOTTOM/DIFF (AND logic + isolate removal)")

    log("[STEP 14.2] Building TOP/BOTTOM graphs by |r| >= threshold ...")

    G_top_raw = build_weighted_graph_from_corr(corr_top_full, threshold=CORR_THRESHOLD)
    G_bot_raw = build_weighted_graph_from_corr(corr_bot_full, threshold=CORR_THRESHOLD)

    # remove isolates NOW (this defines your TOP and BOTTOM gene lists)
    G_top = remove_isolated_nodes(G_top_raw)
    G_bot = remove_isolated_nodes(G_bot_raw)

    genes_top = sorted(G_top.nodes())
    genes_bot = sorted(G_bot.nodes())
    genes_and = sorted(set(genes_top) & set(genes_bot))

    print("\n[STEP 14.2 SUMMARY — TOP/BOTTOM]")
    print(f"  TOP raw:    nodes={G_top_raw.number_of_nodes()} edges={G_top_raw.number_of_edges()}")
    print(f"  TOP noiso:  nodes={G_top.number_of_nodes()} edges={G_top.number_of_edges()}")
    print(f"  BOTTOM raw: nodes={G_bot_raw.number_of_nodes()} edges={G_bot_raw.number_of_edges()}")
    print(f"  BOTTOM noiso:nodes={G_bot.number_of_nodes()} edges={G_bot.number_of_edges()}")
    print(f"  AND genes (TOP ∩ BOTTOM): {len(genes_and)}")

    if len(genes_and) < 10:
        log("[SKIP][STEP 14] Too few AND genes for DIFF network (<10).")
        continue

    # -----------------------------------------------------------------------------
    # STEP 14.2B — DIFF network on AND genes only
    # -----------------------------------------------------------------------------
    log(f"[STEP 14.2B] Building DIFF graph on AND genes by |corr_top - corr_bot| >= {DIFF_THRESHOLD} ...")

    # restrict matrices to AND genes
    corr_top = corr_top_full.loc[genes_and, genes_and]
    corr_bot = corr_bot_full.loc[genes_and, genes_and]
    corr_diff = (corr_top - corr_bot).fillna(0.0)

    # Save restricted matrices (AND-only)
    corr_top_and_path = os.path.join(corr_dir, f"{cancer_type}_corr_TOP_spearman_AND.csv")
    corr_bot_and_path = os.path.join(corr_dir, f"{cancer_type}_corr_BOTTOM_spearman_AND.csv")
    corr_diff_and_path = os.path.join(corr_dir, f"{cancer_type}_corr_DIFF_signed_AND.csv")
    corr_top.to_csv(corr_top_and_path)
    corr_bot.to_csv(corr_bot_and_path)
    corr_diff.to_csv(corr_diff_and_path)

    # DIFF edges: abs(diff) >= threshold
    G_diff_raw = build_weighted_graph_from_corr(corr_diff, threshold=DIFF_THRESHOLD)
    G_diff = remove_isolated_nodes(G_diff_raw)

    print("\n[STEP 14.2 SUMMARY — DIFF]")
    print(f"  DIFF raw:   nodes={G_diff_raw.number_of_nodes()} edges={G_diff_raw.number_of_edges()}")
    print(f"  DIFF noiso: nodes={G_diff.number_of_nodes()} edges={G_diff.number_of_edges()}")


    # ensure abs_weight exists for later layout/community steps
    def ensure_abs_weight(G):
        for u, v, d in G.edges(data=True):
            w = float(d.get("weight", 0.0))
            d["abs_weight"] = abs(w)


    ensure_abs_weight(G_top)
    ensure_abs_weight(G_bot)
    ensure_abs_weight(G_diff)

    # Save gene lists
    pd.Series(genes_top, name="gene").to_csv(os.path.join(net_dir, f"{cancer_type}_genes_TOP_noiso.csv"), index=False)
    pd.Series(genes_bot, name="gene").to_csv(os.path.join(net_dir, f"{cancer_type}_genes_BOTTOM_noiso.csv"),
                                             index=False)
    pd.Series(genes_and, name="gene").to_csv(os.path.join(net_dir, f"{cancer_type}_genes_AND.csv"), index=False)
    pd.Series(sorted(G_diff.nodes()), name="gene").to_csv(os.path.join(net_dir, f"{cancer_type}_genes_DIFF_noiso.csv"),
                                                          index=False)

    log("[STEP 14.2] Saved gene lists (TOP/BOTTOM/AND/DIFF) -> net_dir")
    # biomarker presence in graphs
    def report_gene_presence(tag, G, genes):
        present = [g for g in genes if g in G.nodes]
        missing = [g for g in genes if g not in G.nodes]
        print(f"\n[{tag}] present: {len(present)} | missing: {len(missing)}")
        if present:
            print("  present:", present)
        if missing:
            print("  missing:", missing)


    banner("[STEP 14.2A] Biomarker/A3 presence in graphs (WITH isolates)")
    report_gene_presence("A3 in G_top", G_top, A3_GENES)
    report_gene_presence("A3 in G_bot", G_bot, A3_GENES)
    report_gene_presence("A3 in G_diff", G_diff, A3_GENES)
    report_gene_presence("Biomarkers in G_top", G_top, BIOMARKERS)
    report_gene_presence("Biomarkers in G_bot", G_bot, BIOMARKERS)
    report_gene_presence("Biomarkers in G_diff", G_diff, BIOMARKERS)

    # --- Force A3 nodes to exist in graphs (even if isolated)
    must_keep_nodes = [g for g in A3_GENES if g in selected_genes_corr]  # A3 are in selected_genes
    add_nodes_if_missing(G_top, must_keep_nodes)
    add_nodes_if_missing(G_bot, must_keep_nodes)
    add_nodes_if_missing(G_diff, must_keep_nodes)

    # Biomarkers are not mandatory — only add if you want them visible when isolated too:
    bm_present_in_data = [b for b in BIOMARKERS if b in selected_genes_corr]
    # If you want biomarkers to appear only when they have edges, DO NOT add them here.
    # If you want them visible even when isolated, uncomment:
    # add_nodes_if_missing(G_top, bm_present_in_data)
    # add_nodes_if_missing(G_bot, bm_present_in_data)
    # add_nodes_if_missing(G_diff, bm_present_in_data)

    # -----------------------------------------------------------------------------
    # STEP 14.3 — Heatmaps + dendrograms
    # -----------------------------------------------------------------------------
    banner("[STEP 14.3] Heatmaps + dendrograms")

    genes_top_ord, Z_top = clustered_heatmap_and_dendrogram(
        corr_df=corr_top, out_dir=heat_dir, prefix=f"{cancer_type}_TOP", method="average"
    )
    genes_bot_ord, Z_bot = clustered_heatmap_and_dendrogram(
        corr_df=corr_bot, out_dir=heat_dir, prefix=f"{cancer_type}_BOTTOM", method="average"
    )
    genes_diff_ord, Z_diff = clustered_heatmap_and_dendrogram_diff(
        corr_top=corr_top, corr_bottom=corr_bot, out_dir=heat_dir, prefix=f"{cancer_type}_DIFF", method="average"
    )

    # Save DIFF order + clustered DIFF matrix
    order_path = os.path.join(heat_dir, f"{cancer_type}_DIFF_gene_order.csv")
    pd.Series(genes_diff_ord, name="gene").to_csv(order_path, index=False)

    corr_diff_ord = corr_diff.loc[genes_diff_ord, genes_diff_ord]
    mat_path = os.path.join(heat_dir, f"{cancer_type}_DIFF_corrdiff_clustered.csv")
    corr_diff_ord.to_csv(mat_path)

    log(f"[STEP 14.3] Saved DIFF gene order -> {order_path}")
    log(f"[STEP 14.3] Saved clustered DIFF matrix -> {mat_path}")

    # optional dendrogram cut
    k = 6
    cluster_map = cluster_labels_from_linkage(list(corr_top.index), Z_diff, k=k)
    cluster_csv = os.path.join(heat_dir, f"{cancer_type}_DIFF_dendroClusters_k{k}.csv")
    save_gene_cluster_labels(cluster_map, cluster_csv)
    log(f"[STEP 14.3] Saved dendrogram cluster labels -> {cluster_csv}")

    # -----------------------------------------------------------------------------
    # STEP 14.4 — Save edge lists (with isolates + no-isolates)
    # -----------------------------------------------------------------------------
    banner("[STEP 14.4] Save edge lists (with isolates + no-isolates)")


    def save_edges_both_formats(G, out_prefix):
        save_graph_edgelist(G, out_prefix + "_weighted.csv", weighted=True)
        save_graph_edgelist(G, out_prefix + "_unweighted.csv", weighted=False)


    save_edges_both_formats(G_top, os.path.join(net_dir, f"{cancer_type}_TOP"))
    save_edges_both_formats(G_bot, os.path.join(net_dir, f"{cancer_type}_BOTTOM"))
    save_edges_both_formats(G_diff, os.path.join(net_dir, f"{cancer_type}_DIFF"))

    G_top_noiso = G_top.copy()
    G_bot_noiso = G_bot.copy()
    G_diff_noiso = G_diff.copy()

    save_edges_both_formats(G_top_noiso, os.path.join(net_dir, f"{cancer_type}_TOP_NOISO"))
    save_edges_both_formats(G_bot_noiso, os.path.join(net_dir, f"{cancer_type}_BOTTOM_NOISO"))
    save_edges_both_formats(G_diff_noiso, os.path.join(net_dir, f"{cancer_type}_DIFF_NOISO"))

    log(f"[STEP 14.4] Saved edge lists -> {net_dir}")

    banner("[STEP 14.4A] Biomarker/A3 presence in graphs (NO isolates)")
    report_gene_presence("A3 in G_top_noiso", G_top_noiso, A3_GENES)
    report_gene_presence("A3 in G_bot_noiso", G_bot_noiso, A3_GENES)
    report_gene_presence("A3 in G_diff_noiso", G_diff_noiso, A3_GENES)
    report_gene_presence("Biomarkers in G_top_noiso", G_top_noiso, BIOMARKERS)
    report_gene_presence("Biomarkers in G_bot_noiso", G_bot_noiso, BIOMARKERS)
    report_gene_presence("Biomarkers in G_diff_noiso", G_diff_noiso, BIOMARKERS)

    # -----------------------------------------------------------------------------
    # STEP 14.5 — Common layout + network plots
    # -----------------------------------------------------------------------------
    banner("[STEP 14.5] Layout + network plots")

    A3_ALIASES = ["A3A", "A3B", "A3C", "A3D", "A3F", "A3G", "A3H"]


    def get_a3_nodes_for_graph(G):
        if any(g in G.nodes for g in A3_GENES):
            return [g for g in A3_GENES if g in G.nodes]  # ENSG mode
        if any(a in G.nodes for a in A3_ALIASES):
            return [a for a in A3_ALIASES if a in G.nodes]  # alias mode
        return []


    def get_biomarker_nodes_for_graph(G):
        return [b for b in BIOMARKERS if b in G.nodes]


    def get_highlight_nodes(G):
        return list(dict.fromkeys(get_a3_nodes_for_graph(G) + get_biomarker_nodes_for_graph(G)))


    # ------------------------------------------------------------------
    # MASTER highlight list: fixed biology-first set
    # ------------------------------------------------------------------
    MASTER_HIGHLIGHT_NODES = list(dict.fromkeys(
        [g for g in A3_GENES if g in selected_genes_corr] +
        [b for b in BIOMARKERS if b in selected_genes_corr]
    ))

    # force keep highlight nodes in all graphs if they exist in selected genes
    add_nodes_if_missing(G_top, MASTER_HIGHLIGHT_NODES)
    add_nodes_if_missing(G_bot, MASTER_HIGHLIGHT_NODES)
    add_nodes_if_missing(G_diff, MASTER_HIGHLIGHT_NODES)

    # ------------------------------------------------------------------
    # Layout anchored on DIFF no-isolates preferred
    # ------------------------------------------------------------------
    G_layout = G_diff_noiso if G_diff_noiso.number_of_nodes() > 0 else G_diff
    pos = nx.spring_layout(G_layout, seed=42, weight="abs_weight")

    G_union = nx.compose(nx.compose(G_top, G_bot), G_diff)
    pos = nx9_safe_fill_pos(pos, G_union)

    # ------------------------------------------------------------------
    # Keep highlight nodes even in "no isolate" plots
    # ------------------------------------------------------------------
    G_top_keep = remove_isolates_but_keep(G_top, MASTER_HIGHLIGHT_NODES)
    G_bot_keep = remove_isolates_but_keep(G_bot, MASTER_HIGHLIGHT_NODES)
    G_diff_keep = remove_isolates_but_keep(G_diff, MASTER_HIGHLIGHT_NODES)

    # ------------------------------------------------------------------
    # Plot directories
    # ------------------------------------------------------------------
    plots_with_iso_ids = os.path.join(plot_dir, "with_isolates_weighted_geneID")
    plots_with_iso_symbols = os.path.join(plot_dir, "with_isolates_weighted_geneSymbol")
    plots_no_iso_ids = os.path.join(plot_dir, "no_isolates_weighted_geneID")
    plots_no_iso_symbols = os.path.join(plot_dir, "no_isolates_weighted_geneSymbol")
    zoom_dir = os.path.join(plot_dir, "zoomed_in_plots")
    zoom_dir_Gene_Name = os.path.join(plot_dir, "zoomed_in_plots_Gene_Name")

    for d in [
        plots_with_iso_ids,
        plots_with_iso_symbols,
        plots_no_iso_ids,
        plots_no_iso_symbols,
        zoom_dir,
        zoom_dir_Gene_Name
    ]:
        os.makedirs(d, exist_ok=True)

    # ------------------------------------------------------------------
    # SYMBOL versions of graphs
    # ------------------------------------------------------------------
    G_top_sym, pos_top_sym, name_map_top = relabel_graph_and_pos_to_symbols_safe(G_top, pos, ensg_to_symbol)
    G_bot_sym, pos_bot_sym, name_map_bot = relabel_graph_and_pos_to_symbols_safe(G_bot, pos, ensg_to_symbol)
    G_diff_sym, pos_diff_sym, name_map_diff = relabel_graph_and_pos_to_symbols_safe(G_diff, pos, ensg_to_symbol)

    G_top_keep_sym, pos_top_keep_sym, name_map_top_keep = relabel_graph_and_pos_to_symbols_safe(G_top_keep, pos,
                                                                                                ensg_to_symbol)
    G_bot_keep_sym, pos_bot_keep_sym, name_map_bot_keep = relabel_graph_and_pos_to_symbols_safe(G_bot_keep, pos,
                                                                                                ensg_to_symbol)
    G_diff_keep_sym, pos_diff_keep_sym, name_map_diff_keep = relabel_graph_and_pos_to_symbols_safe(G_diff_keep, pos,
                                                                                                   ensg_to_symbol)

    MASTER_HIGHLIGHT_SYMBOLS_TOP = [name_map_top[g] for g in MASTER_HIGHLIGHT_NODES if g in name_map_top]
    MASTER_HIGHLIGHT_SYMBOLS_BOT = [name_map_bot[g] for g in MASTER_HIGHLIGHT_NODES if g in name_map_bot]
    MASTER_HIGHLIGHT_SYMBOLS_DIFF = [name_map_diff[g] for g in MASTER_HIGHLIGHT_NODES if g in name_map_diff]

    MASTER_HIGHLIGHT_SYMBOLS_TOP_KEEP = [name_map_top_keep[g] for g in MASTER_HIGHLIGHT_NODES if g in name_map_top_keep]
    MASTER_HIGHLIGHT_SYMBOLS_BOT_KEEP = [name_map_bot_keep[g] for g in MASTER_HIGHLIGHT_NODES if g in name_map_bot_keep]
    MASTER_HIGHLIGHT_SYMBOLS_DIFF_KEEP = [name_map_diff_keep[g] for g in MASTER_HIGHLIGHT_NODES if
                                          g in name_map_diff_keep]

    # ------------------------------------------------------------------
    # WITH isolates — gene ID plots
    # ------------------------------------------------------------------
    plot_global_network_fixedpos(
        G_top, pos,
        out_path=os.path.join(plots_with_iso_ids, f"{cancer_type}_TOP_weighted_withIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | TOP | weighted | with isolates | gene ID",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_bot, pos,
        out_path=os.path.join(plots_with_iso_ids, f"{cancer_type}_BOTTOM_weighted_withIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | BOTTOM | weighted | with isolates | gene ID",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_diff, pos,
        out_path=os.path.join(plots_with_iso_ids, f"{cancer_type}_DIFF_weighted_withIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | DIFF | weighted | with isolates | gene ID",
        weighted=True,
        label_apobec=True
    )
    plot_combined_network_fixedpos(
        G_top, G_bot, pos,
        out_path=os.path.join(plots_with_iso_ids, f"{cancer_type}_COMBINED_weighted_withIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | Combined | with isolates | gene ID (TOP=red, BOTTOM=blue)"
    )

    # ------------------------------------------------------------------
    # WITH isolates — gene symbol plots
    # ------------------------------------------------------------------
    plot_global_network_fixedpos(
        G_top_sym, pos_top_sym,
        out_path=os.path.join(plots_with_iso_symbols, f"{cancer_type}_TOP_weighted_withIsolates_geneSymbol.png"),
        a3_genes=MASTER_HIGHLIGHT_SYMBOLS_TOP,
        title=f"{cancer_type} | TOP | weighted | with isolates | gene symbol",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_bot_sym, pos_bot_sym,
        out_path=os.path.join(plots_with_iso_symbols, f"{cancer_type}_BOTTOM_weighted_withIsolates_geneSymbol.png"),
        a3_genes=MASTER_HIGHLIGHT_SYMBOLS_BOT,
        title=f"{cancer_type} | BOTTOM | weighted | with isolates | gene symbol",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_diff_sym, pos_diff_sym,
        out_path=os.path.join(plots_with_iso_symbols, f"{cancer_type}_DIFF_weighted_withIsolates_geneSymbol.png"),
        a3_genes=MASTER_HIGHLIGHT_SYMBOLS_DIFF,
        title=f"{cancer_type} | DIFF | weighted | with isolates | gene symbol",
        weighted=True,
        label_apobec=True
    )

    # ------------------------------------------------------------------
    # WITHOUT isolates — gene ID plots
    # ------------------------------------------------------------------
    plot_global_network_fixedpos(
        G_top_keep, pos,
        out_path=os.path.join(plots_no_iso_ids, f"{cancer_type}_TOP_weighted_noIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | TOP | weighted | no isolates except highlights | gene ID",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_bot_keep, pos,
        out_path=os.path.join(plots_no_iso_ids, f"{cancer_type}_BOTTOM_weighted_noIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | BOTTOM | weighted | no isolates except highlights | gene ID",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_diff_keep, pos,
        out_path=os.path.join(plots_no_iso_ids, f"{cancer_type}_DIFF_weighted_noIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | DIFF | weighted | no isolates except highlights | gene ID",
        weighted=True,
        label_apobec=True
    )
    plot_combined_network_fixedpos(
        G_top_keep, G_bot_keep, pos,
        out_path=os.path.join(plots_no_iso_ids, f"{cancer_type}_COMBINED_weighted_noIsolates_geneID.png"),
        a3_genes=MASTER_HIGHLIGHT_NODES,
        title=f"{cancer_type} | Combined | no isolates except highlights | gene ID (TOP=red, BOTTOM=blue)"
    )

    # ------------------------------------------------------------------
    # WITHOUT isolates — gene symbol plots
    # ------------------------------------------------------------------
    plot_global_network_fixedpos(
        G_top_keep_sym, pos_top_keep_sym,
        out_path=os.path.join(plots_no_iso_symbols, f"{cancer_type}_TOP_weighted_noIsolates_geneSymbol.png"),
        a3_genes=MASTER_HIGHLIGHT_SYMBOLS_TOP_KEEP,
        title=f"{cancer_type} | TOP | weighted | no isolates except highlights | gene symbol",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_bot_keep_sym, pos_bot_keep_sym,
        out_path=os.path.join(plots_no_iso_symbols, f"{cancer_type}_BOTTOM_weighted_noIsolates_geneSymbol.png"),
        a3_genes=MASTER_HIGHLIGHT_SYMBOLS_BOT_KEEP,
        title=f"{cancer_type} | BOTTOM | weighted | no isolates except highlights | gene symbol",
        weighted=True,
        label_apobec=True
    )
    plot_global_network_fixedpos(
        G_diff_keep_sym, pos_diff_keep_sym,
        out_path=os.path.join(plots_no_iso_symbols, f"{cancer_type}_DIFF_weighted_noIsolates_geneSymbol.png"),
        a3_genes=MASTER_HIGHLIGHT_SYMBOLS_DIFF_KEEP,
        title=f"{cancer_type} | DIFF | weighted | no isolates except highlights | gene symbol",
        weighted=True,
        label_apobec=True
    )

    # ------------------------------------------------------------------
    # Metric-based zoom target selection
    # ------------------------------------------------------------------
    metric_change_df = build_metric_change_table(G_top, G_bot, ensg_to_symbol=ensg_to_symbol)
    metric_change_csv = os.path.join(zoom_dir, f"{cancer_type}_TOP_vs_BOTTOM_metric_changes.csv")
    metric_change_df.to_csv(metric_change_csv, index=False)

    top_changed_metric_genes = metric_change_df["gene"].head(12).tolist()

    genes_to_plot = list(dict.fromkeys(
        [g for g in A3_GENES if g in G_union.nodes] +
        [b for b in BIOMARKERS if b in G_union.nodes] +
        top_changed_metric_genes
    ))

    zoom_summary_rows = []

    for gene in genes_to_plot:
        gene_symbol = ensg_to_symbol.get(str(gene), str(gene))

        ego_png_symbols = os.path.join(
            zoom_dir,
            f"{cancer_type}_ego_TOP_BOTTOM_DIFF_{gene_symbol}.png"
        )


        # ENSG version
        ego_png_ids = os.path.join(
            zoom_dir,
            f"{cancer_type}_ego_TOP_BOTTOM_DIFF_{gene}.png"
        )
        plot_zoomed_ego_triptych(
            gene=gene,
            G_top=G_top_keep,
            G_bot=G_bot_keep,
            G_diff=G_diff_keep,
            out_path=ego_png_ids,
            ensg_to_symbol=ensg_to_symbol,
            use_symbols=False
        )


        plot_zoomed_ego_triptych(
            gene=gene,
            G_top=G_top_keep,
            G_bot=G_bot_keep,
            G_diff=G_diff_keep,
            out_path=ego_png_symbols,
            ensg_to_symbol=ensg_to_symbol,
            use_symbols=True
        )

        # SYMBOL version
        ego_png_symbols = os.path.join(
            zoom_dir_Gene_Name,
            f"{cancer_type}_ego_TOP_BOTTOM_DIFF_SYMBOL_{gene_symbol}.png"
        )

        plot_zoomed_ego_triptych(
            gene=gene,
            G_top=G_top_keep,
            G_bot=G_bot_keep,
            G_diff=G_diff_keep,
            out_path=ego_png_symbols,
            ensg_to_symbol=ensg_to_symbol,
            use_symbols=True
        )

        bars_png_ids = os.path.join(
            zoom_dir,
            f"{cancer_type}_metrics_TOP_vs_BOTTOM_{gene}.png"
        )

        plot_zoomed_metric_bars(
            gene=gene,
            metric_df=metric_change_df,
            out_path=bars_png_ids,
            ensg_to_symbol=ensg_to_symbol
        )


        bars_png_symbols = os.path.join(
            zoom_dir_Gene_Name,
            f"{cancer_type}_metrics_TOP_vs_BOTTOM_SYMBOL_{gene_symbol}.png"
        )

        plot_zoomed_metric_bars(
            gene=gene,
            metric_df=metric_change_df,
            out_path=bars_png_symbols,
            ensg_to_symbol=ensg_to_symbol
        )

        row = metric_change_df.loc[metric_change_df["gene"] == gene]
        if not row.empty:
            row = row.iloc[0]
            zoom_summary_rows.append({
                "gene_ensg": gene,
                "gene_symbol": gene_symbol,
                "degree_top": row["degree_top"],
                "degree_bottom": row["degree_bottom"],
                "betweenness_top": row["betweenness_top"],
                "betweenness_bottom": row["betweenness_bottom"],
                "closeness_top": row["closeness_top"],
                "closeness_bottom": row["closeness_bottom"],
                "eigenvector_top": row["eigenvector_top"],
                "eigenvector_bottom": row["eigenvector_bottom"],
                "importance_score": row["importance_score"],
                "ego_plot": ego_png_symbols,
                "metric_plot": bars_png_ids,
            })



            zoom_summary_rows.append({
                "gene_ensg": gene,
                "gene_symbol": gene_symbol,
                "degree_top": row["degree_top"],
                "degree_bottom": row["degree_bottom"],
                "betweenness_top": row["betweenness_top"],
                "betweenness_bottom": row["betweenness_bottom"],
                "closeness_top": row["closeness_top"],
                "closeness_bottom": row["closeness_bottom"],
                "eigenvector_top": row["eigenvector_top"],
                "eigenvector_bottom": row["eigenvector_bottom"],
                "importance_score": row["importance_score"],
                "ego_plot": ego_png_symbols,
                "metric_plot": bars_png_symbols,
            })

    zoom_summary_csv = os.path.join(zoom_dir, f"{cancer_type}_zoomed_gene_summary.csv")
    pd.DataFrame(zoom_summary_rows).to_csv(zoom_summary_csv, index=False)

    log(f"[STEP 14.5] Saved network plots -> {plot_dir}")
    log(f"[STEP 14.5] Saved zoomed plots -> {zoom_dir}")
    log("[STEP 14 COMPLETE] (Community sweep is STEP 15)")



    # =============================================================================
    # STEP 15 — Community detection as a resolution loop on DIFF
    # =============================================================================
    banner("[STEP 15] Community detection resolution sweep on DIFF")

    # ---- knobs
    STEP15_METHOD = "leiden"
    STEP15_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
    STEP15_RUNS_PER_RES = 15
    STEP15_BASE_SEED = 42

    STEP15_USE_LCC = True  # use largest connected component
    TARGET_BIG_COMMS = 8  # merge to top K big communities + Other
    MIN_COMM_SIZE = 10  # min size to keep as its own community (else merge)

    # Root directory for step 15 outputs
    step15_root = os.path.join(cancer_dir, "15_COMMUNITIES")
    os.makedirs(step15_root, exist_ok=True)

    # -----------------------------------------------------------------------------
    # STEP 15.1 — Prepare DIFF graph for community detection
    # -----------------------------------------------------------------------------
    banner("[STEP 15.1] Prepare DIFF graph for community detection")

    # ensure we have a no-isolates version
    try:
        G_diff_noiso
    except NameError:
        G_diff_noiso = remove_isolated_nodes(G_diff)

    G_comm = G_diff_noiso.copy()
    G_comm = nx9_remove_degree0(G_comm)

    # optional largest connected component
    if STEP15_USE_LCC and G_comm.number_of_nodes() > 0:
        G_comm = nx9_largest_component(G_comm)

    # ensure abs_weight exists (non-negative) for Leiden
    for u, v, d in G_comm.edges(data=True):
        w = float(d.get("weight", 0.0))
        d["abs_weight"] = abs(float(d.get("abs_weight", w)))

    log(f"[STEP 15.1] DIFF graph for communities: nodes={G_comm.number_of_nodes()} edges={G_comm.number_of_edges()}")

    if G_comm.number_of_nodes() < 5 or G_comm.number_of_edges() < 3:
        log("[SKIP][STEP 15] DIFF graph too small for community detection.")
        continue

    # Presence checks in community graph
    banner("[STEP 15.1A] A3/Biomarker presence in community graph (G_comm)")


    def report_presence_simple(title, nodes, genes):
        present = [g for g in genes if g in nodes]
        missing = [g for g in genes if g not in nodes]
        print(title)
        print("  present:", present if present else "None")
        print("  missing:", missing if missing else "None")


    report_presence_simple("A3 genes:", set(G_comm.nodes()), A3_GENES)
    report_presence_simple("Biomarkers:", set(G_comm.nodes()), BIOMARKERS)

    nodes_eval = list(G_comm.nodes())


    def partition_to_labels(cm: dict, nodes: list):
        """Convert {node:community} -> integer labels aligned to nodes list."""
        ids = [cm.get(n, -1) for n in nodes]
        uniq = {c: i for i, c in enumerate(sorted(set(ids)))}
        return np.array([uniq[c] for c in ids], dtype=int)


    def save_curve(x, y, title, ylabel, out_path):
        plt.figure(figsize=(8, 5))
        plt.plot(x, y, marker="o")
        plt.xlabel("resolution")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()


    def save_gene_lists_from_map(cm_map: dict, out_csv: str):
        """Save one row per community with gene list."""
        comm_to_genes = {}
        for g, c in cm_map.items():
            comm_to_genes.setdefault(c, []).append(g)

        rows = []
        for c, genes in comm_to_genes.items():
            rows.append({"community": c, "size": len(genes), "genes": ";".join(sorted(genes))})

        pd.DataFrame(rows).sort_values("size", ascending=False).to_csv(out_csv, index=False)


    banner("[STEP 15.1B] Community graph tracking")

    print("DIFF graph before Step 15:")
    print(f"  G_diff nodes={G_diff.number_of_nodes()} edges={G_diff.number_of_edges()}")

    print("DIFF noiso copy entering Step 15:")
    print(f"  G_diff_noiso nodes={G_diff_noiso.number_of_nodes()} edges={G_diff_noiso.number_of_edges()}")

    G_comm_before_degree0 = G_diff_noiso.copy()
    print("Before nx9_remove_degree0:")
    print(f"  nodes={G_comm_before_degree0.number_of_nodes()} edges={G_comm_before_degree0.number_of_edges()}")

    G_comm_after_degree0 = nx9_remove_degree0(G_comm_before_degree0.copy())
    print("After nx9_remove_degree0:")
    print(f"  nodes={G_comm_after_degree0.number_of_nodes()} edges={G_comm_after_degree0.number_of_edges()}")

    if STEP15_USE_LCC and G_comm_after_degree0.number_of_nodes() > 0:
        G_comm_lcc = nx9_largest_component(G_comm_after_degree0.copy())
        print("After largest connected component:")
        print(f"  nodes={G_comm_lcc.number_of_nodes()} edges={G_comm_lcc.number_of_edges()}")

        removed_not_in_lcc = sorted(set(G_comm_after_degree0.nodes()) - set(G_comm_lcc.nodes()))
        print(f"  Nodes removed by LCC filter: {len(removed_not_in_lcc)}")
        print(f"  First 20 removed by LCC: {removed_not_in_lcc[:20]}")
    else:
        print("Largest connected component step not applied.")

    print("A3 in G_diff_noiso:", [g for g in A3_GENES if g in G_diff_noiso.nodes()])
    print("A3 in final G_comm:", [g for g in A3_GENES if g in G_comm.nodes()])
    print("Biomarkers in G_diff_noiso:", [g for g in BIOMARKERS if g in G_diff_noiso.nodes()])
    print("Biomarkers in final G_comm:", [g for g in BIOMARKERS if g in G_comm.nodes()])


    # -----------------------------------------------------------------------------
    # STEP 15.2–15.6 — Resolution sweep loop (save partitions + representative plots)
    # -----------------------------------------------------------------------------
    banner("[STEP 15.2] Resolution sweep loop")

    rows = []

    for r in STEP15_RESOLUTIONS:
        res_tag = f"res{r:.2f}".replace(".", "p")
        res_dir = os.path.join(step15_root, f"resolution_{res_tag}")
        os.makedirs(res_dir, exist_ok=True)

        partitions_dir = os.path.join(res_dir, "partitions_json")
        os.makedirs(partitions_dir, exist_ok=True)

        rep_dir = os.path.join(res_dir, "representative_partition")
        os.makedirs(rep_dir, exist_ok=True)

        labels_list = []
        modularities = []
        ncomms_list = []

        log(f"[STEP 15.2] resolution={r:.2f} | runs={STEP15_RUNS_PER_RES}")

        for i in range(STEP15_RUNS_PER_RES):
            seed = STEP15_BASE_SEED + i

            cm, mod = detect_communities_on_graph_v2(
                G_comm,
                method=STEP15_METHOD,
                seed=seed,
                resolution=float(r),
                weight="abs_weight",
                return_modularity=True
            )

            # save per-run partition
            json_path = os.path.join(partitions_dir, f"partition_seed{seed}.json")
            with open(json_path, "w") as f:
                json.dump(cm, f, indent=2)

            lab = partition_to_labels(cm, nodes_eval)
            labels_list.append(lab)
            ncomms_list.append(len(set(cm.values())))
            modularities.append(mod if mod is not None else np.nan)

            # representative partition = first run
            # representative partition = first run
            if i == 0:
                # -------------------------
                # RAW lists (original ENSG)
                # -------------------------
                raw_list_csv = os.path.join(rep_dir, f"{cancer_type}_raw_gene_lists_{res_tag}.csv")
                save_gene_lists_from_map(cm, raw_list_csv)

                # -------------------------
                # RAW lists (NEW SYMBOL version)
                # -------------------------
                raw_list_sym_csv = os.path.join(rep_dir,
                                                f"{cancer_type}_raw_gene_lists_SYMBOLS_{res_tag}.csv"
                                                )
                save_gene_lists_from_map_symbols(cm, raw_list_sym_csv, ensg_to_symbol)

                # -------------------------
                # RAW plot
                # -------------------------
                pos_r = nx.spring_layout(G_comm, seed=STEP15_BASE_SEED, weight="abs_weight")

                raw_png = os.path.join(rep_dir, f"{cancer_type}_DIFF_comm_rep_RAW_{res_tag}.png")

                plot_network_by_community_fixedpos_symbols(
                    G_comm,
                    pos_r,
                    cm,
                    out_path=raw_png,
                    ensg_to_symbol=ensg_to_symbol,
                    title=f"{cancer_type} | DIFF | RAW communities | res={r:.2f}",
                    a3_genes=get_highlight_nodes(G_top),
                    weighted=True
                )

                # -------------------------
                # merge small communities
                # -------------------------
                cm_merged, raw_sizes, merged_sizes = nx9_merge_communities_to_topK(
                    cm, k_keep=TARGET_BIG_COMMS, min_size=MIN_COMM_SIZE
                )

                # -------------------------
                # merged partition map (original ENSG)
                # -------------------------
                merged_map_csv = os.path.join(
                    rep_dir,
                    f"{cancer_type}_rep_partition_merged_top{TARGET_BIG_COMMS}_{res_tag}.csv"
                )

                pd.DataFrame({
                    "gene": list(cm_merged.keys()),
                    "community": list(cm_merged.values())
                }).to_csv(merged_map_csv, index=False)

                # -------------------------
                # merged partition map (NEW SYMBOL version)
                # -------------------------
                merged_map_sym_csv = os.path.join(
                    rep_dir,
                    f"{cancer_type}_rep_partition_merged_SYMBOLS_top{TARGET_BIG_COMMS}_{res_tag}.csv"
                )

                pd.DataFrame({
                    "gene_symbol": [ensg_to_symbol.get(g, g) for g in cm_merged.keys()],
                    "community": list(cm_merged.values())
                }).to_csv(merged_map_sym_csv, index=False)

                # -------------------------
                # merged gene lists (original ENSG)
                # -------------------------
                merged_list_csv = os.path.join(
                    rep_dir,
                    f"{cancer_type}_merged_gene_lists_top{TARGET_BIG_COMMS}_{res_tag}.csv"
                )

                save_gene_lists_from_map(cm_merged, merged_list_csv)

                # -------------------------
                # merged gene lists (NEW SYMBOL version)
                # -------------------------
                merged_list_sym_csv = os.path.join(
                    rep_dir,
                    f"{cancer_type}_merged_gene_lists_SYMBOLS_top{TARGET_BIG_COMMS}_{res_tag}.csv"
                )

                save_gene_lists_from_map_symbols(cm_merged, merged_list_sym_csv, ensg_to_symbol)

                # -------------------------
                # merged community plot
                # -------------------------
                merged_png = os.path.join(
                    rep_dir,
                    f"{cancer_type}_DIFF_comm_rep_MERGED_top{TARGET_BIG_COMMS}_{res_tag}.png"
                )

                plot_network_by_community_fixedpos_symbols(
                    G_comm,
                    pos_r,
                    cm_merged,
                    out_path=merged_png,
                    ensg_to_symbol=ensg_to_symbol,
                    title=f"{cancer_type} | DIFF | merged communities (top{TARGET_BIG_COMMS}+Other) | res={r:.2f}",
                    a3_genes=get_highlight_nodes(G_top),
                    weighted=True
                )

                # -------------------------
                # community size plot
                # -------------------------
                sizes_sorted = sorted(merged_sizes.items(), key=lambda x: x[1], reverse=True)

                labels = [str(k) for k, _ in sizes_sorted]
                values = [v for _, v in sizes_sorted]

                plt.figure(figsize=(8, 5))
                plt.bar(labels, values)
                plt.xlabel("Community ID")
                plt.ylabel("Number of genes")
                plt.title(f"{cancer_type} | Community sizes | res={r:.2f}")

                plt.tight_layout()

                size_plot_path = os.path.join(
                    rep_dir,
                    f"{cancer_type}_community_size_barplot_{res_tag}.png"
                )

                plt.savefig(size_plot_path, dpi=300)
                plt.close()

                log(f"[STEP 15.4] Saved representative RAW plot -> {raw_png}")
                log(f"[STEP 15.4] Saved representative MERGED plot -> {merged_png}")
                log(f"[STEP 15.4] Saved merged gene lists -> {merged_list_csv}")
        # -----------------------------------------------------------------------------
        # STEP 15.7 — Stability metrics (ARI/NMI) per resolution
        # -----------------------------------------------------------------------------
        aris = []
        nmis = []
        for a, b in itertools.combinations(range(STEP15_RUNS_PER_RES), 2):
            aris.append(adjusted_rand_score(labels_list[a], labels_list[b]))
            nmis.append(normalized_mutual_info_score(labels_list[a], labels_list[b]))

        rows.append({
            "resolution": float(r),
            "n_runs": STEP15_RUNS_PER_RES,
            "ncomms_mean": float(np.mean(ncomms_list)),
            "ncomms_std": float(np.std(ncomms_list)),
            "modularity_mean": float(np.nanmean(modularities)),
            "modularity_std": float(np.nanstd(modularities)),
            "ARI_mean": float(np.mean(aris)),
            "ARI_std": float(np.std(aris)),
            "NMI_mean": float(np.mean(nmis)),
            "NMI_std": float(np.std(nmis)),
        })

    # -----------------------------------------------------------------------------
    # STEP 15.8 — Save sweep table + plots
    # -----------------------------------------------------------------------------
    banner("[STEP 15.8] Save sweep table + plots")

    sweep_df = pd.DataFrame(rows).sort_values("resolution")
    sweep_csv = os.path.join(step15_root, f"{cancer_type}_resolution_sweep.csv")
    sweep_df.to_csv(sweep_csv, index=False)
    log(f"[STEP 15.8] Saved resolution sweep CSV -> {sweep_csv}")

    save_curve(
        sweep_df["resolution"], sweep_df["ncomms_mean"],
        title=f"{cancer_type} | #communities vs resolution",
        ylabel="#communities (mean)",
        out_path=os.path.join(step15_root, f"{cancer_type}_ncomms_vs_resolution.png")
    )
    save_curve(
        sweep_df["resolution"], sweep_df["ARI_mean"],
        title=f"{cancer_type} | stability vs resolution (ARI)",
        ylabel="ARI (mean)",
        out_path=os.path.join(step15_root, f"{cancer_type}_ARI_vs_resolution.png")
    )
    save_curve(
        sweep_df["resolution"], sweep_df["NMI_mean"],
        title=f"{cancer_type} | stability vs resolution (NMI)",
        ylabel="NMI (mean)",
        out_path=os.path.join(step15_root, f"{cancer_type}_NMI_vs_resolution.png")
    )
    save_curve(
        sweep_df["resolution"], sweep_df["modularity_mean"],
        title=f"{cancer_type} | modularity vs resolution",
        ylabel="Modularity (mean)",
        out_path=os.path.join(step15_root, f"{cancer_type}_modularity_vs_resolution.png")
    )

    log(f"[STEP 15.8] Saved sweep plots -> {step15_root}")

    # -----------------------------------------------------------------------------
    # STEP 15.9 — Recommend resolution (simple rule)
    # -----------------------------------------------------------------------------
    banner("[STEP 15.9] Recommend a resolution")

    MIN_STABILITY = 0.75
    cand = sweep_df[sweep_df["ARI_mean"] >= MIN_STABILITY]
    if len(cand) > 0:
        best_row = cand.sort_values("modularity_mean", ascending=False).iloc[0]
    else:
        best_row = sweep_df.sort_values(["ARI_mean", "modularity_mean"], ascending=False).iloc[0]

    BEST_RESOLUTION = float(best_row["resolution"])
    log(
        f"[STEP 15.9] Recommended resolution={BEST_RESOLUTION:.2f} "
        f"(ARI={best_row['ARI_mean']:.3f}, NMI={best_row['NMI_mean']:.3f}, "
        f"Mod={best_row['modularity_mean']:.3f}, ncomms~{best_row['ncomms_mean']:.1f})"
    )

    log("[STEP 15 COMPLETE]")

    # =============================================================================
    # STEP 16 — Per-gene network metrics (TOP/BOTTOM/DIFF)
    # =============================================================================
    banner("[STEP 16] Per-gene network metrics (TOP/BOTTOM/DIFF)")

    step16_dir = os.path.join(cancer_dir, "16_network_metrics")
    os.makedirs(step16_dir, exist_ok=True)

    # ---- compute metrics (your util)
    m_top = compute_graph_metrics_df(G_top)
    m_bot = compute_graph_metrics_df(G_bot)
    m_diff = compute_graph_metrics_df(G_diff)

    # ---- save tables
    top_csv = os.path.join(step16_dir, f"{cancer_type}_TOP_metrics.csv")
    bot_csv = os.path.join(step16_dir, f"{cancer_type}_BOTTOM_metrics.csv")
    diff_csv = os.path.join(step16_dir, f"{cancer_type}_DIFF_metrics.csv")

    m_top.to_csv(top_csv, index=False)
    m_bot.to_csv(bot_csv, index=False)
    m_diff.to_csv(diff_csv, index=False)

    log(f"[STEP 16] Saved metrics CSVs -> {step16_dir}")

    # ---- plots (your util)
    plot_metric_distributions(m_top, step16_dir, f"{cancer_type}_TOP")
    plot_metric_distributions(m_bot, step16_dir, f"{cancer_type}_BOTTOM")
    plot_metric_distributions(m_diff, step16_dir, f"{cancer_type}_DIFF")

    log(f"[STEP 16] Saved metric distribution plots -> {step16_dir}")

    # =============================================================================
    # STEP 16A — Check A3 genes + biomarkers exist in metrics outputs
    # =============================================================================
    banner("[STEP 16A] A3/Biomarker presence in metrics tables")


    def metrics_presence_report(df_metrics: pd.DataFrame, label: str, genes: list[str]):
        # different metric functions sometimes name the node column differently
        if "gene" in df_metrics.columns:
            node_col = "gene"
        elif "node" in df_metrics.columns:
            node_col = "node"
        else:
            print(f"[{label}] ⚠ Cannot find a node column (expected 'gene' or 'node').")
            print("Columns:", df_metrics.columns.tolist())
            return

        nodes = set(df_metrics[node_col].astype(str).tolist())
        present = [g for g in genes if g in nodes]
        missing = [g for g in genes if g not in nodes]

        print(f"[{label}] present={len(present)} missing={len(missing)}")
        if present:
            print("  present:", present)
        if missing:
            print("  missing:", missing)


    metrics_presence_report(m_top, "TOP metrics | A3", A3_GENES)
    metrics_presence_report(m_bot, "BOTTOM metrics | A3", A3_GENES)
    metrics_presence_report(m_diff, "DIFF metrics | A3", A3_GENES)

    metrics_presence_report(m_top, "TOP metrics | biomarkers", BIOMARKERS)
    metrics_presence_report(m_bot, "BOTTOM metrics | biomarkers", BIOMARKERS)
    metrics_presence_report(m_diff, "DIFF metrics | biomarkers", BIOMARKERS)

    log(f"[STEP 16 COMPLETE] Finished metrics for {cancer_type}")




# =============================================================================
# FINAL STEP — Save a copy of this script into the main run directory
# =============================================================================
banner("[FINAL STEP] Save a copy of the script")

import shutil
from datetime import datetime

try:
    # path of the currently running script
    current_script_path = os.path.abspath(__file__)

    # optional timestamped backup name
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    script_name = os.path.basename(current_script_path)

    # save one plain copy
    dst_plain = os.path.join(RUN_ROOT, script_name)
    shutil.copy2(current_script_path, dst_plain)

    # save one timestamped backup copy
    dst_backup = os.path.join(
        RUN_ROOT,
        f"{os.path.splitext(script_name)[0]}__backup_{ts}.py"
    )
    shutil.copy2(current_script_path, dst_backup)

    log(f"[FINAL STEP] Saved script copy -> {dst_plain}")
    log(f"[FINAL STEP] Saved timestamped backup -> {dst_backup}")

except NameError:
    log("[FINAL STEP] WARNING: __file__ is not available in this environment.")
    log("[FINAL STEP] This usually happens in notebooks or interactive consoles.")
except Exception as e:
    log(f"[FINAL STEP] ERROR while saving script copy: {e}")