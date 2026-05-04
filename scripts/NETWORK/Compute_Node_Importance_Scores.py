#!/usr/bin/env python3
"""
Compute_Node_Importance_Scores.py
==================================
Compute dual importance scores for all genes in the HNSC differential
co-expression network:

  1. INTRA score ("local hub"): How central is this gene within its
     own community? Drives node size in Figure 2.

  2. INTER score ("bridge"): How much does this gene connect different
     communities? Drives border/shape annotation in Figure 2.

Components:
  INTRA: intra-community degree, intra-community strength,
         within-community eigenvector centrality
  INTER: inter-community degree, global betweenness centrality,
         number of distinct communities connected to

Outputs:
  - Supp_Table_Node_Scores.tsv (publication-ready, all 616 genes)
  - Node_Score_Summary_Report.txt (console report)

Usage:
    conda run -n NETWORK python Compute_Node_Importance_Scores.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os, json
import numpy as np
import pandas as pd
import networkx as nx

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CANCER_TYPE = "TCGA-HNSC"

DIFF_EDGES = os.path.join(
    PROJECT_ROOT, "data/FIG_2/04_correlation_networks", CANCER_TYPE,
    "edge_lists", f"{CANCER_TYPE}_DIFF_edges.csv"
)
BEST_PARTITION = os.path.join(
    PROJECT_ROOT, "data/FIG_2/05_communities", CANCER_TYPE,
    f"{CANCER_TYPE}_best_partition.csv"
)
ENSG_TO_SYMBOL = os.path.join(
    PROJECT_ROOT, "data/FIG_2/01_cleaned_expression/ensg_to_symbol.json"
)
KEGG_SUMMARY = os.path.join(
    PROJECT_ROOT, "data/FIG_2/08_pipeline_summary",
    f"{CANCER_TYPE}_KEGG_enrichment_summary.csv"
)

# Harris interactors
HARRIS_ALL_PATH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/00_input/Harris_A3_interactors.txt"
)
HARRIS_A3B_PATH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt"
)

A3_SYMBOLS = {
    'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
    'APOBEC3F', 'APOBEC3G', 'APOBEC3H',
}

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_2/DIAGNOSTIC_AUDIT")
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
banner("LOAD DATA")

# ENSG to symbol
with open(ENSG_TO_SYMBOL) as f:
    ensg_to_symbol = json.load(f)
log(f"  ENSG-to-symbol mapping: {len(ensg_to_symbol)} entries")

def sym(ensg):
    return ensg_to_symbol.get(str(ensg), str(ensg))

# Community assignments
partition = pd.read_csv(BEST_PARTITION)
partition['community'] = partition['community'].astype(int)
gene_to_comm = dict(zip(partition['gene'], partition['community']))
comm_genes = partition.groupby('community')['gene'].apply(set).to_dict()
log(f"  Partition: {len(partition)} genes in {len(comm_genes)} communities")

# DIFF edges
edges_df = pd.read_csv(DIFF_EDGES)
log(f"  DIFF edges: {len(edges_df)}")

# Build networkx graph (community genes only)
community_gene_set = set(partition['gene'])

G = nx.Graph()
for _, row in edges_df.iterrows():
    s, t = row['source'], row['target']
    # Only include edges where both nodes are in a community
    if s in community_gene_set and t in community_gene_set:
        G.add_edge(s, t, weight=row['weight'], abs_weight=row['abs_weight'])

# Add isolated community genes (degree=0 within community subgraph)
for gene in community_gene_set:
    if gene not in G:
        G.add_node(gene)

log(f"  Community subgraph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Harris interactors
with open(HARRIS_ALL_PATH) as f:
    harris_all = set(line.strip() for line in f if line.strip())
with open(HARRIS_A3B_PATH) as f:
    harris_a3b = set(line.strip() for line in f if line.strip())

# KEGG summary
kegg_summary = pd.read_csv(KEGG_SUMMARY)
comm_kegg = {}
for _, row in kegg_summary.iterrows():
    cid = int(row['community'])
    term = str(row.get('top_kegg_term', 'N/A'))
    pval = row.get('top_kegg_pvalue', None)
    comm_kegg[cid] = (term, pval)


# =============================================================================
# CLASSIFY EDGES: INTRA vs INTER
# =============================================================================
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
log(f"  Intra fraction: {len(intra_edges)/G.number_of_edges():.1%}")


# =============================================================================
# COMPUTE INTRA-COMMUNITY METRICS
# =============================================================================
banner("COMPUTE INTRA-COMMUNITY METRICS")

# Build per-community subgraphs for eigenvector computation
intra_degree = {}
intra_strength = {}
intra_eigenvector = {}

for cid, genes in sorted(comm_genes.items()):
    # Build subgraph for this community
    subgraph = G.subgraph(genes).copy()

    # Intra-degree and intra-strength per node
    for gene in genes:
        # Count only edges to nodes in the same community
        ideg = 0
        istr = 0.0
        if gene in subgraph:
            for neighbor in subgraph.neighbors(gene):
                if neighbor in genes:
                    ideg += 1
                    istr += abs(subgraph[gene][neighbor].get('abs_weight', 0))
        intra_degree[gene] = ideg
        intra_strength[gene] = istr

    # Within-community eigenvector centrality
    # Only compute if subgraph has edges (eigenvector undefined for isolated nodes)
    sub_edges = subgraph.number_of_edges()
    if sub_edges > 0:
        try:
            eig = nx.eigenvector_centrality_numpy(subgraph, weight='abs_weight')
            for gene in genes:
                intra_eigenvector[gene] = eig.get(gene, 0.0)
        except Exception:
            # Fallback: degree centrality within community
            for gene in genes:
                intra_eigenvector[gene] = intra_degree[gene] / max(len(genes) - 1, 1)
    else:
        for gene in genes:
            intra_eigenvector[gene] = 0.0

    n_intra = sum(1 for g in genes if intra_degree[g] > 0)
    top_hub = max(genes, key=lambda g: intra_degree[g])
    log(f"  C{cid} ({len(genes)} genes): {sub_edges} intra-edges, "
        f"{n_intra} connected, top hub: {sym(top_hub)}({intra_degree[top_hub]})")


# =============================================================================
# COMPUTE INTER-COMMUNITY METRICS
# =============================================================================
banner("COMPUTE INTER-COMMUNITY METRICS")

inter_degree = {}
inter_strength = {}
communities_connected = {}  # gene -> set of community IDs

for gene in community_gene_set:
    inter_degree[gene] = 0
    inter_strength[gene] = 0.0
    communities_connected[gene] = set()

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

# Global betweenness (already on the full community subgraph)
log(f"\n  Computing global betweenness centrality...")
global_betweenness = nx.betweenness_centrality(G, weight='abs_weight')

# Report top inter-community bridges
inter_sorted = sorted(community_gene_set, key=lambda g: inter_degree.get(g, 0), reverse=True)
log(f"\n  Top 15 inter-community bridges:")
log(f"  {'Gene':15s} {'Comm':>5s} {'InterDeg':>9s} {'InterStr':>9s} "
    f"{'#Comms':>7s} {'GlobBtw':>9s}")
for gene in inter_sorted[:15]:
    cid = gene_to_comm.get(gene, -1)
    log(f"  {sym(gene):15s} {cid:>5d} {inter_degree[gene]:>9d} "
        f"{inter_strength[gene]:>9.3f} {n_comms_connected[gene]:>7d} "
        f"{global_betweenness.get(gene, 0):>9.5f}")


# =============================================================================
# BUILD COMPOSITE SCORES
# =============================================================================
banner("BUILD COMPOSITE SCORES")

rows = []
for gene in community_gene_set:
    cid = gene_to_comm[gene]
    symbol = sym(gene)

    rows.append({
        'gene_symbol': symbol,
        'ensg_id': gene,
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
        # Total degree (for reference)
        'total_degree': intra_degree.get(gene, 0) + inter_degree.get(gene, 0),
        # Flags
        'is_a3_gene': symbol in A3_SYMBOLS,
        'is_harris_all': symbol in harris_all,
        'is_harris_a3b': symbol in harris_a3b,
    })

df = pd.DataFrame(rows)

# Compute percentile ranks for intra components
for col in ['intra_degree', 'intra_strength', 'intra_eigenvector']:
    df[f'{col}_pctrank'] = df[col].rank(pct=True)

# INTRA composite: mean percentile rank of the three intra metrics
df['intra_score'] = df[[
    'intra_degree_pctrank',
    'intra_strength_pctrank',
    'intra_eigenvector_pctrank',
]].mean(axis=1).round(4)

# Compute percentile ranks for inter components
for col in ['inter_degree', 'inter_strength', 'n_communities_connected', 'global_betweenness']:
    df[f'{col}_pctrank'] = df[col].rank(pct=True)

# INTER composite: mean percentile rank of the four inter metrics
df['inter_score'] = df[[
    'inter_degree_pctrank',
    'inter_strength_pctrank',
    'n_communities_connected_pctrank',
    'global_betweenness_pctrank',
]].mean(axis=1).round(4)

# Drop intermediate rank columns
rank_cols = [c for c in df.columns if c.endswith('_pctrank')]
df = df.drop(columns=rank_cols)

# Sort by intra_score descending (primary), inter_score descending (secondary)
df = df.sort_values(['intra_score', 'inter_score'], ascending=[False, False])


# =============================================================================
# REPORT
# =============================================================================
banner("TOP 25 GENES BY INTRA SCORE (Local Hubs)")

log(f"  {'Gene':15s} {'Comm':>5s} {'InDeg':>6s} {'InStr':>7s} {'InEig':>7s} "
    f"{'INTRA':>6s} {'ExDeg':>6s} {'INTER':>6s}")
log(f"  {'-'*15} {'-'*5} {'-'*6} {'-'*7} {'-'*7} {'-'*6} {'-'*6} {'-'*6}")
for _, row in df.head(25).iterrows():
    log(f"  {row['gene_symbol']:15s} {int(row['community']):>5d} "
        f"{int(row['intra_degree']):>6d} {row['intra_strength']:>7.2f} "
        f"{row['intra_eigenvector']:>7.4f} {row['intra_score']:>6.4f} "
        f"{int(row['inter_degree']):>6d} {row['inter_score']:>6.4f}")


banner("TOP 25 GENES BY INTER SCORE (Bridges)")

df_inter = df.sort_values('inter_score', ascending=False)
log(f"  {'Gene':15s} {'Comm':>5s} {'ExDeg':>6s} {'ExStr':>7s} {'#Comm':>6s} "
    f"{'GloBtw':>8s} {'INTER':>6s} {'InDeg':>6s} {'INTRA':>6s}")
log(f"  {'-'*15} {'-'*5} {'-'*6} {'-'*7} {'-'*6} {'-'*8} {'-'*6} {'-'*6} {'-'*6}")
for _, row in df_inter.head(25).iterrows():
    log(f"  {row['gene_symbol']:15s} {int(row['community']):>5d} "
        f"{int(row['inter_degree']):>6d} {row['inter_strength']:>7.2f} "
        f"{int(row['n_communities_connected']):>6d} "
        f"{row['global_betweenness']:>8.5f} {row['inter_score']:>6.4f} "
        f"{int(row['intra_degree']):>6d} {row['intra_score']:>6.4f}")


banner("PER-COMMUNITY TOP 3 LOCAL HUBS")

for cid in sorted(comm_genes.keys()):
    c_df = df[df['community'] == cid].sort_values('intra_score', ascending=False)
    n_genes = len(c_df)
    kegg_term, kegg_p = comm_kegg.get(cid, ('N/A', None))
    kegg_str = f"{kegg_term}" if kegg_p is None else f"{kegg_term} (p={kegg_p:.2e})"

    top3 = c_df.head(3)
    hub_str = ', '.join(
        f"{r['gene_symbol']}(in={int(r['intra_degree'])},ex={int(r['inter_degree'])})"
        for _, r in top3.iterrows()
    )

    # Check for A3 genes
    a3_in = [r['gene_symbol'] for _, r in c_df.iterrows() if r['is_a3_gene']]
    a3_str = f" | A3: {', '.join(a3_in)}" if a3_in else ""

    log(f"\n  C{cid} ({n_genes} genes){a3_str}")
    log(f"    KEGG: {kegg_str}")
    log(f"    Top hubs: {hub_str}")

    # Any bridges?
    bridges = c_df[c_df['inter_degree'] > 0].sort_values('inter_score', ascending=False)
    if len(bridges) > 0:
        n_bridge = len(bridges)
        top_bridge = bridges.iloc[0]
        log(f"    Bridges: {n_bridge} genes with inter-community edges "
            f"(top: {top_bridge['gene_symbol']}, "
            f"ex_deg={int(top_bridge['inter_degree'])}, "
            f"to {int(top_bridge['n_communities_connected'])} communities)")


# =============================================================================
# SAVE
# =============================================================================
banner("SAVE OUTPUTS")

out_path = os.path.join(OUTPUT_DIR, "Supp_Table_Node_Scores.tsv")
df_save = df.sort_values(['community', 'intra_score'], ascending=[True, False])
df_save.to_csv(out_path, sep='\t', index=False)
log(f"  Node scores table: {out_path}")
log(f"  Total genes: {len(df_save)}")

# Also save a figure-ready version with just the columns needed for node sizing
fig_cols = ['gene_symbol', 'ensg_id', 'community', 'intra_score', 'inter_score',
            'intra_degree', 'inter_degree', 'total_degree']
fig_path = os.path.join(OUTPUT_DIR, "Figure2_Node_Sizing.tsv")
df_save[fig_cols].to_csv(fig_path, sep='\t', index=False)
log(f"  Figure 2 sizing table: {fig_path}")

# Report
report_path = os.path.join(OUTPUT_DIR, "Node_Score_Summary_Report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")

banner("DONE")
