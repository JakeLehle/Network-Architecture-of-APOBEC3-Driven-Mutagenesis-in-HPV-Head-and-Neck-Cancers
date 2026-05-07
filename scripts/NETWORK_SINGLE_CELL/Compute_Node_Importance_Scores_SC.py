#!/usr/bin/env python3
"""
Compute_Node_Importance_Scores_SC.py
======================================

Compute dual importance scores for all genes in the single-cell differential
co-expression network (Figure 4). Adapted from the bulk TCGA version
(Compute_Node_Importance_Scores.py) used for Figure 2.

Scoring:
  1. INTRA score ("local hub"): How central is this gene within its
     own community? Drives node size in Figure 4 panels.

  2. INTER score ("bridge"): How much does this gene connect different
     communities? Drives border/shape annotation in Figure 4 panels.

Components:
  INTRA: intra-community degree, intra-community strength,
         within-community eigenvector centrality
  INTER: inter-community degree, global betweenness centrality,
         number of distinct communities connected to

Key differences from bulk version:
  - SC network uses gene symbols directly (no ENSG-to-symbol mapping)
  - Edge list is TSV (not CSV)
  - No KEGG summary available (skipped in reporting)
  - Output to data/FIG_4/ instead of data/FIG_2/

Outputs:
  - data/FIG_4/DIAGNOSTIC_AUDIT/Supp_Table_SC_Node_Scores.tsv
  - data/FIG_4/DIAGNOSTIC_AUDIT/Figure4_Node_Sizing.tsv
  - data/FIG_4/DIAGNOSTIC_AUDIT/SC_Node_Score_Summary_Report.txt

Usage:
    conda run -n NETWORK python Compute_Node_Importance_Scores_SC.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import numpy as np
import pandas as pd
import networkx as nx
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input paths (SC pipeline outputs)
DIFF_EDGES = os.path.join(
    PROJECT_ROOT, "data/FIG_4/03_correlation_networks/edge_lists/SC_edges_DIFF.tsv"
)
BEST_PARTITION = os.path.join(
    PROJECT_ROOT, "data/FIG_4/04_communities/SC_best_partition.csv"
)
COMM_GRAPH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/04_communities/SC_G_comm.gpickle"
)

# Harris interactors
HARRIS_ALL_PATH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/00_input/Harris_A3_interactors.txt"
)
HARRIS_A3B_PATH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt"
)

# TCGA shared genes (from overlap analysis, if available)
TCGA_OVERLAP_PATH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/06_overlap_analysis/overlap_gene_details.csv"
)

A3_SYMBOLS = {
    'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
    'APOBEC3F', 'APOBEC3G', 'APOBEC3H',
}

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_4/DIAGNOSTIC_AUDIT")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# LOGGING
# =============================================================================

report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)


# =============================================================================
# LOAD DATA
# =============================================================================

def load_data():
    """Load partition, edge list, Harris interactors, and TCGA shared genes."""
    banner("LOAD DATA")

    # Community assignments
    partition = pd.read_csv(BEST_PARTITION)
    partition['community'] = partition['community'].astype(int)
    gene_to_comm = dict(zip(partition['gene'], partition['community']))
    comm_genes = partition.groupby('community')['gene'].apply(set).to_dict()
    log(f"  Partition: {len(partition)} genes in {len(comm_genes)} communities")

    # DIFF edges (TSV for SC pipeline)
    edges_df = pd.read_csv(DIFF_EDGES, sep='\t')
    log(f"  DIFF edges: {len(edges_df)}")
    log(f"  Edge columns: {list(edges_df.columns)}")

    # Build networkx graph from edges (community genes only)
    community_gene_set = set(partition['gene'])

    G = nx.Graph()
    for _, row in edges_df.iterrows():
        s, t = row['source'], row['target']
        if s in community_gene_set and t in community_gene_set:
            G.add_edge(s, t,
                       weight=row.get('weight', 0),
                       abs_weight=row.get('abs_weight', abs(row.get('weight', 0))))

    # Add isolated community genes
    for gene in community_gene_set:
        if gene not in G:
            G.add_node(gene)

    log(f"  Community subgraph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Harris interactors
    harris_all = set()
    harris_a3b = set()
    if os.path.exists(HARRIS_ALL_PATH):
        with open(HARRIS_ALL_PATH) as f:
            harris_all = set(line.strip() for line in f if line.strip())
        log(f"  Harris interactors (all): {len(harris_all)}")
    if os.path.exists(HARRIS_A3B_PATH):
        with open(HARRIS_A3B_PATH) as f:
            harris_a3b = set(line.strip() for line in f if line.strip())
        log(f"  Harris interactors (A3B): {len(harris_a3b)}")

    # TCGA shared genes (optional)
    tcga_shared = set()
    if os.path.exists(TCGA_OVERLAP_PATH):
        overlap_df = pd.read_csv(TCGA_OVERLAP_PATH)
        if 'gene' in overlap_df.columns:
            tcga_shared = set(overlap_df['gene'].values)
        elif 'gene_symbol' in overlap_df.columns:
            tcga_shared = set(overlap_df['gene_symbol'].values)
        log(f"  TCGA shared genes: {len(tcga_shared)}")
    else:
        log(f"  TCGA overlap file not found (skipping)")

    return G, gene_to_comm, comm_genes, community_gene_set, harris_all, harris_a3b, tcga_shared


# =============================================================================
# CLASSIFY EDGES
# =============================================================================

def classify_edges(G, gene_to_comm):
    """Split edges into intra-community and inter-community."""
    banner("CLASSIFY EDGES")

    intra_edges = []
    inter_edges = []

    for u, v, data in G.edges(data=True):
        comm_u = gene_to_comm.get(u, -1)
        comm_v = gene_to_comm.get(v, -1)
        if comm_u == comm_v and comm_u >= 0:
            intra_edges.append((u, v, data))
        else:
            inter_edges.append((u, v, data))

    log(f"  Intra-community edges: {len(intra_edges)}")
    log(f"  Inter-community edges: {len(inter_edges)}")
    if G.number_of_edges() > 0:
        log(f"  Intra fraction: {len(intra_edges)/G.number_of_edges():.1%}")

    return intra_edges, inter_edges


# =============================================================================
# COMPUTE INTRA-COMMUNITY METRICS
# =============================================================================

def compute_intra(G, comm_genes, gene_to_comm, community_gene_set):
    """Compute intra-community degree, strength, and eigenvector centrality."""
    banner("COMPUTE INTRA-COMMUNITY METRICS")

    intra_degree = {}
    intra_strength = {}
    intra_eigenvector = {}

    for cid, genes in sorted(comm_genes.items()):
        subgraph = G.subgraph(genes).copy()

        for gene in genes:
            ideg = 0
            istr = 0.0
            if gene in subgraph:
                for neighbor in subgraph.neighbors(gene):
                    if neighbor in genes:
                        ideg += 1
                        istr += abs(subgraph[gene][neighbor].get('abs_weight', 0))
            intra_degree[gene] = ideg
            intra_strength[gene] = istr

        sub_edges = subgraph.number_of_edges()
        if sub_edges > 0:
            try:
                eig = nx.eigenvector_centrality_numpy(subgraph, weight='abs_weight')
                for gene in genes:
                    intra_eigenvector[gene] = eig.get(gene, 0.0)
            except Exception:
                for gene in genes:
                    intra_eigenvector[gene] = intra_degree[gene] / max(len(genes) - 1, 1)
        else:
            for gene in genes:
                intra_eigenvector[gene] = 0.0

        n_intra = sum(1 for g in genes if intra_degree[g] > 0)
        top_hub = max(genes, key=lambda g: intra_degree[g])
        log(f"  C{cid} ({len(genes)} genes): {sub_edges} intra-edges, "
            f"{n_intra} connected, top hub: {top_hub}({intra_degree[top_hub]})")

    return intra_degree, intra_strength, intra_eigenvector


# =============================================================================
# COMPUTE INTER-COMMUNITY METRICS
# =============================================================================

def compute_inter(G, inter_edges, gene_to_comm, community_gene_set):
    """Compute inter-community degree, strength, betweenness, and bridge count."""
    banner("COMPUTE INTER-COMMUNITY METRICS")

    inter_degree = {g: 0 for g in community_gene_set}
    inter_strength = {g: 0.0 for g in community_gene_set}
    communities_connected = {g: set() for g in community_gene_set}

    for u, v, data in inter_edges:
        comm_u = gene_to_comm.get(u, -1)
        comm_v = gene_to_comm.get(v, -1)

        inter_degree[u] = inter_degree.get(u, 0) + 1
        inter_degree[v] = inter_degree.get(v, 0) + 1

        w = abs(data.get('abs_weight', 0))
        inter_strength[u] = inter_strength.get(u, 0) + w
        inter_strength[v] = inter_strength.get(v, 0) + w

        if comm_v >= 0:
            communities_connected[u].add(comm_v)
        if comm_u >= 0:
            communities_connected[v].add(comm_u)

    n_comms_connected = {g: len(c) for g, c in communities_connected.items()}

    log(f"\n  Computing global betweenness centrality...")
    global_betweenness = nx.betweenness_centrality(G, weight='abs_weight')

    # Report top bridges
    inter_sorted = sorted(community_gene_set, key=lambda g: inter_degree.get(g, 0), reverse=True)
    log(f"\n  Top 15 inter-community bridges:")
    log(f"  {'Gene':15s} {'Comm':>5s} {'InterDeg':>9s} {'InterStr':>9s} "
        f"{'#Comms':>7s} {'GlobBtw':>9s}")
    for gene in inter_sorted[:15]:
        cid = gene_to_comm.get(gene, -1)
        log(f"  {gene:15s} {cid:>5d} {inter_degree[gene]:>9d} "
            f"{inter_strength[gene]:>9.3f} {n_comms_connected[gene]:>7d} "
            f"{global_betweenness.get(gene, 0):>9.5f}")

    return inter_degree, inter_strength, n_comms_connected, global_betweenness


# =============================================================================
# BUILD COMPOSITE SCORES
# =============================================================================

def build_scores(community_gene_set, gene_to_comm,
                 intra_degree, intra_strength, intra_eigenvector,
                 inter_degree, inter_strength, n_comms_connected,
                 global_betweenness, harris_all, harris_a3b, tcga_shared):
    """Build composite intra/inter scores and return DataFrame."""
    banner("BUILD COMPOSITE SCORES")

    rows = []
    for gene in community_gene_set:
        cid = gene_to_comm[gene]

        rows.append({
            'gene_symbol': gene,
            'community': cid,
            # Raw intra metrics
            'intra_degree': intra_degree.get(gene, 0),
            'intra_strength': round(intra_strength.get(gene, 0), 4),
            'intra_eigenvector': round(intra_eigenvector.get(gene, 0), 6),
            # Raw inter metrics
            'inter_degree': inter_degree.get(gene, 0),
            'inter_strength': round(inter_strength.get(gene, 0), 4),
            'n_communities_connected': n_comms_connected.get(gene, 0),
            'global_betweenness': round(global_betweenness.get(gene, 0), 6),
            # Total degree
            'total_degree': intra_degree.get(gene, 0) + inter_degree.get(gene, 0),
            # Flags
            'is_a3_gene': gene in A3_SYMBOLS,
            'is_harris_all': gene in harris_all,
            'is_harris_a3b': gene in harris_a3b,
            'is_tcga_shared': gene in tcga_shared,
        })

    df = pd.DataFrame(rows)

    # Percentile ranks for intra components
    for col in ['intra_degree', 'intra_strength', 'intra_eigenvector']:
        df[f'{col}_pctrank'] = df[col].rank(pct=True)

    df['intra_score'] = df[[
        'intra_degree_pctrank',
        'intra_strength_pctrank',
        'intra_eigenvector_pctrank',
    ]].mean(axis=1).round(4)

    # Percentile ranks for inter components
    for col in ['inter_degree', 'inter_strength', 'n_communities_connected', 'global_betweenness']:
        df[f'{col}_pctrank'] = df[col].rank(pct=True)

    df['inter_score'] = df[[
        'inter_degree_pctrank',
        'inter_strength_pctrank',
        'n_communities_connected_pctrank',
        'global_betweenness_pctrank',
    ]].mean(axis=1).round(4)

    # Drop intermediate rank columns
    rank_cols = [c for c in df.columns if c.endswith('_pctrank')]
    df = df.drop(columns=rank_cols)

    # Sort by intra_score descending
    df = df.sort_values(['intra_score', 'inter_score'], ascending=[False, False])

    return df


# =============================================================================
# REPORTING
# =============================================================================

def report_scores(df, comm_genes, harris_all, harris_a3b, tcga_shared):
    """Print summary tables."""

    banner("TOP 25 GENES BY INTRA SCORE (Local Hubs)")
    log(f"  {'Gene':15s} {'Comm':>5s} {'InDeg':>6s} {'InStr':>7s} {'InEig':>7s} "
        f"{'INTRA':>6s} {'ExDeg':>6s} {'INTER':>6s} {'Flags':>12s}")
    log(f"  {'-'*15} {'-'*5} {'-'*6} {'-'*7} {'-'*7} {'-'*6} {'-'*6} {'-'*6} {'-'*12}")
    for _, row in df.head(25).iterrows():
        flags = []
        if row['is_a3_gene']:      flags.append('A3')
        if row['is_harris_all']:   flags.append('Harris')
        if row['is_tcga_shared']:  flags.append('TCGA')
        flag_str = ','.join(flags) if flags else ''
        log(f"  {row['gene_symbol']:15s} {int(row['community']):>5d} "
            f"{int(row['intra_degree']):>6d} {row['intra_strength']:>7.2f} "
            f"{row['intra_eigenvector']:>7.4f} {row['intra_score']:>6.4f} "
            f"{int(row['inter_degree']):>6d} {row['inter_score']:>6.4f} "
            f"{flag_str:>12s}")

    banner("TOP 25 GENES BY INTER SCORE (Bridges)")
    df_inter = df.sort_values('inter_score', ascending=False)
    log(f"  {'Gene':15s} {'Comm':>5s} {'ExDeg':>6s} {'ExStr':>7s} {'#Comm':>6s} "
        f"{'GloBtw':>8s} {'INTER':>6s} {'INTRA':>6s}")
    log(f"  {'-'*15} {'-'*5} {'-'*6} {'-'*7} {'-'*6} {'-'*8} {'-'*6} {'-'*6}")
    for _, row in df_inter.head(25).iterrows():
        log(f"  {row['gene_symbol']:15s} {int(row['community']):>5d} "
            f"{int(row['inter_degree']):>6d} {row['inter_strength']:>7.2f} "
            f"{int(row['n_communities_connected']):>6d} "
            f"{row['global_betweenness']:>8.5f} {row['inter_score']:>6.4f} "
            f"{row['intra_score']:>6.4f}")

    banner("PER-COMMUNITY TOP 3 LOCAL HUBS")
    for cid in sorted(comm_genes.keys()):
        c_df = df[df['community'] == cid].sort_values('intra_score', ascending=False)
        n_genes = len(c_df)

        top3 = c_df.head(3)
        hub_str = ', '.join(
            f"{r['gene_symbol']}(in={int(r['intra_degree'])},ex={int(r['inter_degree'])})"
            for _, r in top3.iterrows()
        )

        a3_in = [r['gene_symbol'] for _, r in c_df.iterrows() if r['is_a3_gene']]
        harris_in = [r['gene_symbol'] for _, r in c_df.iterrows() if r['is_harris_all']]
        tcga_in = [r['gene_symbol'] for _, r in c_df.iterrows() if r['is_tcga_shared']]

        tags = []
        if a3_in:    tags.append(f"A3: {', '.join(a3_in)}")
        if harris_in: tags.append(f"Harris: {len(harris_in)}")
        if tcga_in:  tags.append(f"TCGA shared: {len(tcga_in)}")
        tag_str = f" | {' | '.join(tags)}" if tags else ""

        log(f"\n  C{cid} ({n_genes} genes){tag_str}")
        log(f"    Top hubs: {hub_str}")

        bridges = c_df[c_df['inter_degree'] > 0].sort_values('inter_score', ascending=False)
        if len(bridges) > 0:
            top_bridge = bridges.iloc[0]
            log(f"    Bridges: {len(bridges)} genes (top: {top_bridge['gene_symbol']}, "
                f"ex_deg={int(top_bridge['inter_degree'])}, "
                f"to {int(top_bridge['n_communities_connected'])} communities)")

    # Special interest genes
    banner("SPECIAL INTEREST GENES")

    # A3 genes in network
    a3_df = df[df['is_a3_gene']].sort_values('intra_score', ascending=False)
    log(f"\n  A3 genes in network: {len(a3_df)}")
    for _, row in a3_df.iterrows():
        log(f"    {row['gene_symbol']:12s} C{int(row['community']):>2d}  "
            f"intra={row['intra_score']:.4f}  inter={row['inter_score']:.4f}  "
            f"deg={int(row['total_degree'])}")

    # Harris interactors in network
    harris_df = df[df['is_harris_all']].sort_values('intra_score', ascending=False)
    log(f"\n  Harris A3 interactors in network: {len(harris_df)}")
    for _, row in harris_df.iterrows():
        a3b_flag = " [A3B-specific]" if row['is_harris_a3b'] else ""
        log(f"    {row['gene_symbol']:12s} C{int(row['community']):>2d}  "
            f"intra={row['intra_score']:.4f}  inter={row['inter_score']:.4f}  "
            f"deg={int(row['total_degree'])}{a3b_flag}")


# =============================================================================
# SAVE
# =============================================================================

def save_outputs(df):
    """Save node scores and figure-ready sizing table."""
    banner("SAVE OUTPUTS")

    # Full scores table
    out_path = os.path.join(OUTPUT_DIR, "Supp_Table_SC_Node_Scores.tsv")
    df_save = df.sort_values(['community', 'intra_score'], ascending=[True, False])
    df_save.to_csv(out_path, sep='\t', index=False)
    log(f"  Node scores table: {out_path}")
    log(f"  Total genes: {len(df_save)}")

    # Figure-ready sizing table (matches Figure2_Node_Sizing.tsv format)
    fig_cols = ['gene_symbol', 'community', 'intra_score', 'inter_score',
                'intra_degree', 'inter_degree', 'total_degree']
    fig_path = os.path.join(OUTPUT_DIR, "Figure4_Node_Sizing.tsv")
    df_save[fig_cols].to_csv(fig_path, sep='\t', index=False)
    log(f"  Figure 4 sizing table: {fig_path}")

    # Report
    report_path = os.path.join(OUTPUT_DIR, "SC_Node_Score_Summary_Report.txt")
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    log(f"  Report: {report_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("COMPUTE NODE IMPORTANCE SCORES -- SINGLE CELL (Figure 4)")
    log(f"  Start: {t0}")

    # Load
    (G, gene_to_comm, comm_genes, community_gene_set,
     harris_all, harris_a3b, tcga_shared) = load_data()

    # Classify edges
    intra_edges, inter_edges = classify_edges(G, gene_to_comm)

    # Compute metrics
    intra_degree, intra_strength, intra_eigenvector = compute_intra(
        G, comm_genes, gene_to_comm, community_gene_set)

    inter_degree, inter_strength, n_comms_connected, global_betweenness = compute_inter(
        G, inter_edges, gene_to_comm, community_gene_set)

    # Build scores
    df = build_scores(
        community_gene_set, gene_to_comm,
        intra_degree, intra_strength, intra_eigenvector,
        inter_degree, inter_strength, n_comms_connected,
        global_betweenness, harris_all, harris_a3b, tcga_shared)

    # Report
    report_scores(df, comm_genes, harris_all, harris_a3b, tcga_shared)

    # Save
    save_outputs(df)

    banner("DONE")
    log(f"  Elapsed: {datetime.now() - t0}")


if __name__ == "__main__":
    main()
