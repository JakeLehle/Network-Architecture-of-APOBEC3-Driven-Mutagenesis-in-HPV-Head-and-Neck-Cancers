#!/usr/bin/env python3
"""
Diagnostic_Network_Community_Audit.py
======================================
Post-rerun audit of the updated HNSC differential co-expression network.
Produces a consolidated report and supplemental-ready tables.

Reads:
  - Community gene lists (05_communities)
  - DIFF metrics with community assignments (06_centrality_metrics)
  - KEGG enrichment summary (08_pipeline_summary)
  - Harris A3 interactor lists (FIG_4/00_input)
  - ENSG-to-symbol mapping (01_cleaned_expression)
  - Pipeline report JSON (08_pipeline_summary)

Outputs:
  - Console report (+ saved as .txt)
  - Supplemental Table S_network (community summary, .tsv)
  - Supplemental Table S_nodes (full node manifest, .tsv)

Usage:
    conda run -n NETWORK python Diagnostic_Network_Community_Audit.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os, json
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CANCER_TYPE = "TCGA-HNSC"

# Input paths
COMMUNITY_GENE_LISTS = os.path.join(
    PROJECT_ROOT, "data/FIG_2/05_communities", CANCER_TYPE,
    f"{CANCER_TYPE}_community_gene_lists.csv"
)
DIFF_METRICS = os.path.join(
    PROJECT_ROOT, "data/FIG_2/06_centrality_metrics", CANCER_TYPE,
    f"{CANCER_TYPE}_DIFF_metrics_with_communities.csv"
)
KEGG_SUMMARY = os.path.join(
    PROJECT_ROOT, "data/FIG_2/08_pipeline_summary",
    f"{CANCER_TYPE}_KEGG_enrichment_summary.csv"
)
PIPELINE_JSON = os.path.join(
    PROJECT_ROOT, "data/FIG_2/08_pipeline_summary",
    f"{CANCER_TYPE}_pipeline_report.json"
)
ENSG_TO_SYMBOL = os.path.join(
    PROJECT_ROOT, "data/FIG_2/01_cleaned_expression/ensg_to_symbol.json"
)
HARRIS_ALL_PATH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/00_input/Harris_A3_interactors.txt"
)
HARRIS_A3B_PATH = os.path.join(
    PROJECT_ROOT, "data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt"
)

# DDR and chromatin gene sets (same as somatic enrichment script)
DDR_GENES = set([
    "OGG1", "MUTYH", "NEIL1", "NEIL2", "NEIL3", "MPG", "NTHL1", "SMUG1",
    "TDG", "UNG", "APEX1", "APEX2", "PNKP", "XRCC1", "PARP1", "PARP2",
    "LIG3", "POLB", "FEN1", "XPA", "XPC", "RAD23B", "CETN2", "DDB1",
    "DDB2", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6", "ERCC8",
    "RPA1", "RPA2", "RPA3", "LIG1", "MLH1", "MLH3", "MSH2", "MSH3",
    "MSH6", "PMS1", "PMS2", "EXO1", "PCNA", "BRCA1", "BRCA2", "RAD51",
    "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3", "RAD54L", "RAD54B",
    "NBN", "MRE11", "RAD50", "PALB2", "BRIP1", "BLM", "WRN", "RECQL",
    "RECQL4", "RECQL5", "ATM", "ATR", "ATRIP", "CHEK1", "CHEK2", "MDC1",
    "RNF8", "RNF168", "TP53BP1", "RBBP8", "XRCC4", "XRCC5", "XRCC6",
    "LIG4", "DCLRE1C", "PRKDC", "NHEJ1", "FANCA", "FANCB", "FANCC",
    "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM",
    "UBE2T", "SLX4", "REV1", "REV3L", "POLH", "POLI", "POLK", "POLQ",
])

CHROMATIN_GENES = set([
    "SMARCA2", "SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1",
    "SMARCD2", "SMARCD3", "SMARCE1", "ARID1A", "ARID1B", "ARID2",
    "PBRM1", "BRD7", "BRD9", "DPF1", "DPF2", "DPF3", "PHF10",
    "KDM6A", "KMT2A", "KMT2C", "KMT2D", "SETD2", "NSD1", "NSD2",
    "EZH2", "SUZ12", "EED", "EP300", "CREBBP", "KAT6A", "KAT6B",
    "DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3",
    "IDH1", "IDH2", "ATRX", "DAXX", "H3-3A", "H3-3B",
    "CHD1", "CHD2", "CHD3", "CHD4", "CHD5", "CHD6", "CHD7", "CHD8",
])

# A3 genes
A3_SYMBOLS = {
    'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
    'APOBEC3F', 'APOBEC3G', 'APOBEC3H',
}

# Output
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

# Community gene lists
comm_df = pd.read_csv(COMMUNITY_GENE_LISTS)
community_genes = {}
for _, row in comm_df.iterrows():
    cid = int(row['community'])
    genes = [g.strip() for g in str(row['genes']).split(';') if g.strip()]
    community_genes[cid] = genes
total_genes = sum(len(g) for g in community_genes.values())
log(f"  Communities: {len(community_genes)}")
log(f"  Total genes in communities: {total_genes}")

# DIFF metrics
metrics = pd.read_csv(DIFF_METRICS)
log(f"  DIFF metrics: {len(metrics)} genes")

# Ensure community column is clean
if 'community' in metrics.columns:
    metrics['community'] = metrics['community'].apply(
        lambda x: int(x) if pd.notna(x) else -1
    )

# Harris interactors
with open(HARRIS_ALL_PATH) as f:
    harris_all = set(line.strip() for line in f if line.strip())
with open(HARRIS_A3B_PATH) as f:
    harris_a3b = set(line.strip() for line in f if line.strip())
log(f"  Harris A3 interactors (all): {len(harris_all)}")
log(f"  Harris A3B-specific: {len(harris_a3b)}")

# KEGG summary
kegg_summary = None
if os.path.exists(KEGG_SUMMARY):
    kegg_summary = pd.read_csv(KEGG_SUMMARY)
    log(f"  KEGG summary: {len(kegg_summary)} communities")

# Pipeline JSON
pipeline_info = None
if os.path.exists(PIPELINE_JSON):
    with open(PIPELINE_JSON) as f:
        pipeline_info = json.load(f)
    log(f"  Pipeline JSON loaded")


# =============================================================================
# SECTION 1: Network-level summary
# =============================================================================
banner("SECTION 1: Network-Level Summary")

# Extract from pipeline JSON if available
if pipeline_info:
    for key in ['n_high', 'n_low', 'sbs2_median_high', 'sbs2_median_low',
                'a3a_median_high', 'a3a_median_low', 'a3b_median_high', 'a3b_median_low',
                'genes_tested', 'genes_significant', 'genes_selected',
                'best_resolution', 'best_modularity', 'best_ari']:
        if key in pipeline_info:
            log(f"  {key}: {pipeline_info[key]}")

# Compute from metrics
in_community = metrics[metrics['community'] >= 0]
log(f"\n  Genes in any community: {len(in_community)}")
log(f"  Genes with degree > 0: {(metrics['degree'] > 0).sum()}")
log(f"  Genes with degree = 0 (isolates): {(metrics['degree'] == 0).sum()}")

# Degree distribution
deg_vals = in_community['degree']
log(f"\n  Degree distribution (community genes only):")
log(f"    Mean: {deg_vals.mean():.1f}")
log(f"    Median: {deg_vals.median():.0f}")
log(f"    Max: {deg_vals.max()}")
log(f"    Genes with degree >= 10: {(deg_vals >= 10).sum()}")
log(f"    Genes with degree >= 5: {(deg_vals >= 5).sum()}")


# =============================================================================
# SECTION 2: Per-Community Summary Table
# =============================================================================
banner("SECTION 2: Per-Community Summary Table")

summary_rows = []

for cid in sorted(community_genes.keys()):
    genes = community_genes[cid]
    n_genes = len(genes)

    # Get metrics for this community
    c_metrics = in_community[in_community['community'] == cid].copy()

    # Top 3 hubs by degree
    top_hubs = c_metrics.nlargest(3, 'degree')
    hub_str = ', '.join(
        f"{row['gene_symbol']}({int(row['degree'])})"
        for _, row in top_hubs.iterrows()
    )

    # A3 genes
    a3_in = [g for g in genes if g in A3_SYMBOLS]
    a3_str = ', '.join(a3_in) if a3_in else 'None'

    # Harris interactors
    harris_in = [g for g in genes if g in harris_all]
    harris_a3b_in = [g for g in genes if g in harris_a3b]

    # DDR genes
    ddr_in = [g for g in genes if g in DDR_GENES]

    # Chromatin genes
    chrom_in = [g for g in genes if g in CHROMATIN_GENES]

    # KEGG
    kegg_term = "N/A"
    kegg_p = "N/A"
    n_kegg_sig = 0
    if kegg_summary is not None:
        k_row = kegg_summary[kegg_summary['community'] == cid]
        if len(k_row) > 0:
            k = k_row.iloc[0]
            kegg_term = str(k.get('top_kegg_term', 'N/A'))
            kegg_p_val = k.get('top_kegg_pvalue', None)
            kegg_p = f"{kegg_p_val:.2e}" if pd.notna(kegg_p_val) else "N/A"
            n_kegg_sig = int(k.get('n_kegg_significant', 0))

    # Mean metrics
    mean_deg = c_metrics['degree'].mean() if len(c_metrics) > 0 else 0
    mean_btw = c_metrics['betweenness'].mean() if len(c_metrics) > 0 else 0
    max_eig = c_metrics['eigenvector'].max() if len(c_metrics) > 0 else 0

    summary_rows.append({
        'community': cid,
        'n_genes': n_genes,
        'mean_degree': round(mean_deg, 1),
        'max_eigenvector': round(max_eig, 4),
        'top_hubs': hub_str,
        'a3_genes': a3_str,
        'n_harris_all': len(harris_in),
        'harris_genes': ';'.join(sorted(harris_in)) if harris_in else '',
        'n_harris_a3b': len(harris_a3b_in),
        'n_ddr': len(ddr_in),
        'ddr_genes': ';'.join(sorted(ddr_in)) if ddr_in else '',
        'n_chromatin': len(chrom_in),
        'chromatin_genes': ';'.join(sorted(chrom_in)) if chrom_in else '',
        'top_kegg_term': kegg_term,
        'kegg_adj_p': kegg_p,
        'n_kegg_significant': n_kegg_sig,
    })

    # Print summary
    log(f"\n  Community {cid} ({n_genes} genes)")
    log(f"    Top hubs: {hub_str}")
    log(f"    A3 genes: {a3_str}")
    log(f"    Harris interactors: {len(harris_in)}"
        f"{' (' + ', '.join(sorted(harris_in)) + ')' if harris_in else ''}")
    if harris_a3b_in:
        log(f"    Harris A3B-specific: {', '.join(sorted(harris_a3b_in))}")
    if ddr_in:
        log(f"    DDR genes: {', '.join(sorted(ddr_in))}")
    if chrom_in:
        log(f"    Chromatin remodelers: {', '.join(sorted(chrom_in))}")
    log(f"    Top KEGG: {kegg_term} (adj p={kegg_p}, {n_kegg_sig} sig terms)")

summary_df = pd.DataFrame(summary_rows)
summary_path = os.path.join(OUTPUT_DIR, "Supp_Table_Network_Community_Summary.tsv")
summary_df.to_csv(summary_path, sep='\t', index=False)
log(f"\n  Saved: {summary_path}")


# =============================================================================
# SECTION 3: Community 4 (A3B) Detailed Gene Table
# =============================================================================
banner("SECTION 3: Community 4 (APOBEC3B) Detailed Audit")

c4_genes = community_genes.get(4, [])
c4_metrics = in_community[in_community['community'] == 4].copy()
c4_metrics = c4_metrics.sort_values('degree', ascending=False)

log(f"  Community 4: {len(c4_genes)} genes")
log(f"  A3B degree in DIFF network: 1 (peripheral, force-kept through DE)")
log(f"  A3B eigenvector centrality: {c4_metrics[c4_metrics['gene_symbol']=='APOBEC3B']['eigenvector'].values[0] if len(c4_metrics[c4_metrics['gene_symbol']=='APOBEC3B']) > 0 else 'N/A':.6f}")

log(f"\n  Wnt pathway genes: WNT5B, WNT7B, WNT16 (3 ligands)")
log(f"  Wnt-related: RAC3 (downstream effector)")
log(f"  Immune-related: NLRC5 (MHC-I transactivator), MALT1 (NF-kB)")

log(f"\n  Full gene table (sorted by degree):")
log(f"  {'Gene':15s} {'Deg':>4s} {'Btw':>8s} {'Eig':>8s} {'Harris':>7s} {'A3B_sp':>7s} {'DDR':>4s} {'Chrom':>6s}")
log(f"  {'-'*15} {'-'*4} {'-'*8} {'-'*8} {'-'*7} {'-'*7} {'-'*4} {'-'*6}")

for _, row in c4_metrics.iterrows():
    gene = row['gene_symbol']
    is_harris = 'Y' if gene in harris_all else ''
    is_a3b_sp = 'Y' if gene in harris_a3b else ''
    is_ddr = 'Y' if gene in DDR_GENES else ''
    is_chrom = 'Y' if gene in CHROMATIN_GENES else ''
    log(f"  {gene:15s} {int(row['degree']):>4d} {row['betweenness']:>8.5f} "
        f"{row['eigenvector']:>8.5f} {is_harris:>7s} {is_a3b_sp:>7s} "
        f"{is_ddr:>4s} {is_chrom:>6s}")


# =============================================================================
# SECTION 4: Community 11 (A3H) Brief Summary
# =============================================================================
banner("SECTION 4: Community 11 (APOBEC3H) Summary")

c11_genes = community_genes.get(11, [])
c11_metrics = in_community[in_community['community'] == 11].copy()
c11_metrics = c11_metrics.sort_values('degree', ascending=False)

log(f"  Community 11: {len(c11_genes)} genes")
a3h_row = c11_metrics[c11_metrics['gene_symbol'] == 'APOBEC3H']
if len(a3h_row) > 0:
    log(f"  A3H degree: {int(a3h_row.iloc[0]['degree'])}")
    log(f"  A3H betweenness: {a3h_row.iloc[0]['betweenness']:.5f}")
    log(f"  A3H eigenvector: {a3h_row.iloc[0]['eigenvector']:.5f}")

log(f"\n  Notable genes: {', '.join(c11_genes[:15])}")

harris_c11 = [g for g in c11_genes if g in harris_all]
if harris_c11:
    log(f"  Harris interactors: {', '.join(sorted(harris_c11))}")
else:
    log(f"  Harris interactors: None")

log(f"\n  Full gene list:")
for _, row in c11_metrics.iterrows():
    gene = row['gene_symbol']
    log(f"    {gene:15s} deg={int(row['degree']):>3d}  btw={row['betweenness']:.5f}"
        f"  eig={row['eigenvector']:.5f}"
        f"{'  [Harris]' if gene in harris_all else ''}")


# =============================================================================
# SECTION 5: Harris Interactor Cross-Reference (All Communities)
# =============================================================================
banner("SECTION 5: Harris A3 Interactors in Network")

all_community_genes = set()
for genes in community_genes.values():
    all_community_genes.update(genes)

harris_in_network = harris_all & all_community_genes
harris_a3b_in_network = harris_a3b & all_community_genes

log(f"  Harris A3 interactors in any community: {len(harris_in_network)} / {len(harris_all)}")
log(f"  Harris A3B-specific in any community: {len(harris_a3b_in_network)} / {len(harris_a3b)}")

if harris_in_network:
    log(f"\n  {'Gene':15s} {'Community':>10s} {'Degree':>7s} {'Btw':>8s} {'A3B_sp':>7s}")
    log(f"  {'-'*15} {'-'*10} {'-'*7} {'-'*8} {'-'*7}")
    for gene in sorted(harris_in_network):
        # Find community
        comm = 'N/A'
        for cid, genes in community_genes.items():
            if gene in genes:
                comm = str(cid)
                break
        # Find metrics
        g_row = metrics[metrics['gene_symbol'] == gene]
        if len(g_row) > 0:
            deg = int(g_row.iloc[0]['degree'])
            btw = g_row.iloc[0]['betweenness']
        else:
            deg = 0
            btw = 0.0
        is_a3b_sp = 'Y' if gene in harris_a3b else ''
        log(f"  {gene:15s} {comm:>10s} {deg:>7d} {btw:>8.5f} {is_a3b_sp:>7s}")
else:
    log(f"\n  No Harris interactors found in any community")

# Also check which Harris interactors are in the DIFF network but not in communities
harris_in_metrics = set()
for _, row in metrics.iterrows():
    if row['gene_symbol'] in harris_all and row['degree'] > 0:
        harris_in_metrics.add(row['gene_symbol'])

harris_diff_only = harris_in_metrics - harris_in_network
if harris_diff_only:
    log(f"\n  Harris interactors in DIFF network but NOT in communities:")
    for gene in sorted(harris_diff_only):
        g_row = metrics[metrics['gene_symbol'] == gene]
        if len(g_row) > 0:
            deg = int(g_row.iloc[0]['degree'])
            log(f"    {gene}: degree={deg}")


# =============================================================================
# SECTION 6: DDR and Chromatin Remodeler Cross-Reference
# =============================================================================
banner("SECTION 6: DDR and Chromatin Genes in Network")

ddr_in_network = DDR_GENES & all_community_genes
chrom_in_network = CHROMATIN_GENES & all_community_genes

log(f"  DDR genes in communities: {len(ddr_in_network)} / {len(DDR_GENES)}")
if ddr_in_network:
    for gene in sorted(ddr_in_network):
        comm = 'N/A'
        for cid, genes in community_genes.items():
            if gene in genes:
                comm = str(cid)
                break
        log(f"    {gene}: Community {comm}")

log(f"\n  Chromatin remodelers in communities: {len(chrom_in_network)} / {len(CHROMATIN_GENES)}")
if chrom_in_network:
    for gene in sorted(chrom_in_network):
        comm = 'N/A'
        for cid, genes in community_genes.items():
            if gene in genes:
                comm = str(cid)
                break
        log(f"    {gene}: Community {comm}")


# =============================================================================
# SECTION 7: Full Node Manifest (Supplemental Table)
# =============================================================================
banner("SECTION 7: Build Full Node Manifest (Supplemental Table)")

# Use only genes in communities
node_rows = []
for _, row in in_community.iterrows():
    gene = row['gene_symbol']
    cid = int(row['community'])

    # Composite score: normalized rank across degree, betweenness, eigenvector
    # (will compute after collecting all rows)
    node_rows.append({
        'gene_symbol': gene,
        'ensg_id': row['node'] if 'node' in row.index else '',
        'community': cid,
        'degree': int(row['degree']),
        'betweenness': round(row['betweenness'], 6),
        'eigenvector': round(row['eigenvector'], 6),
        'strength_abs': round(row['strength_abs'], 4) if 'strength_abs' in row.index else 0,
        'is_a3_gene': gene in A3_SYMBOLS,
        'is_harris_interactor': gene in harris_all,
        'is_harris_a3b_specific': gene in harris_a3b,
        'is_ddr': gene in DDR_GENES,
        'is_chromatin': gene in CHROMATIN_GENES,
    })

nodes_df = pd.DataFrame(node_rows)

# Composite score: mean of percentile ranks for degree, betweenness, eigenvector
for col in ['degree', 'betweenness', 'eigenvector']:
    nodes_df[f'{col}_pctrank'] = nodes_df[col].rank(pct=True)

nodes_df['composite_score'] = nodes_df[
    ['degree_pctrank', 'betweenness_pctrank', 'eigenvector_pctrank']
].mean(axis=1).round(4)

# Drop intermediate rank columns
nodes_df = nodes_df.drop(columns=['degree_pctrank', 'betweenness_pctrank', 'eigenvector_pctrank'])

# Sort by composite score descending
nodes_df = nodes_df.sort_values('composite_score', ascending=False)

nodes_path = os.path.join(OUTPUT_DIR, "Supp_Table_Node_Manifest.tsv")
nodes_df.to_csv(nodes_path, sep='\t', index=False)
log(f"  Node manifest: {len(nodes_df)} genes")
log(f"  Saved: {nodes_path}")

# Top 20 by composite score
log(f"\n  Top 20 genes by composite score:")
log(f"  {'Gene':15s} {'Comm':>5s} {'Deg':>4s} {'Btw':>8s} {'Eig':>8s} {'Comp':>6s} {'Flags':>20s}")
for _, row in nodes_df.head(20).iterrows():
    flags = []
    if row['is_a3_gene']:
        flags.append('A3')
    if row['is_harris_interactor']:
        flags.append('Harris')
    if row['is_ddr']:
        flags.append('DDR')
    if row['is_chromatin']:
        flags.append('Chrom')
    flag_str = ','.join(flags) if flags else ''
    log(f"  {row['gene_symbol']:15s} {int(row['community']):>5d} {int(row['degree']):>4d} "
        f"{row['betweenness']:>8.5f} {row['eigenvector']:>8.5f} "
        f"{row['composite_score']:>6.4f} {flag_str:>20s}")


# =============================================================================
# SECTION 8: Key Narrative Points
# =============================================================================
banner("SECTION 8: Key Narrative Points for Text")

log(f"""
  NETWORK STRUCTURE:
    14 communities, 616 genes, modularity 0.61, resolution 0.80
    Tighter than old run (was 13 communities, 1437 genes, mod 0.53)
    Higher modularity = more internally coherent communities

  A3 GENE PLACEMENT:
    APOBEC3B: Community 4 (50 genes), degree=1, peripheral
      -> force-kept through DE, not a co-expression hub
      -> consistent with SC finding: A3B cofactors are post-transcriptional
    APOBEC3H: Community 11 (25 genes), degree=4
    A3A, A3C, A3D, A3F, A3G: NOT in any community

  COMMUNITY 4 (A3B):
    Dominated by Wnt signaling: WNT5B, WNT7B, WNT16 (3 ligands)
    KEGG: Basal cell carcinoma (adj p=0.017), Wnt signaling (adj p=0.017)
    Biologically coherent: Wnt maintains basal stem cell identity
    A3B is peripheral (degree=1) but community context is relevant

  CELL-TYPE MARKERS:
    0/24 canonical markers in any community (was 2/24 in old run)
    KRT5, TACSTD2, CD68, CD3E, etc. all absent
    Strongest possible case for "bulk cannot resolve cell types"

  COMPARISON TO OLD RUN:
    Old: A3B in C2 with KRT5, A3G in C0 with immune pattern
    New: A3B in C4 with Wnt, A3H in C11, A3G dropped out entirely
    KRT5 lost because it was not DE (p > 0.05) in the updated groups
    A3G lost because it was not DE (p = 0.97)
    The tighter selection (1643 vs 2810 DE genes) pruned noise

  IMPLICATIONS FOR SECTION 4.2:
    Old narrative (A3B+KRT5, A3G+immune) needs complete rewrite
    New narrative: A3B assigned to Wnt/basal-cell-identity community,
      but A3B is peripheral (degree=1), reinforcing that A3B's regulatory
      relationships are not detectable at the co-expression level
    Focus on Community 4's Wnt biology as the cellular context
    0/24 markers strengthens the "need single-cell" argument
""")


# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUTPUT_DIR, "Network_Community_Audit_Report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Full report: {report_path}")

banner("AUDIT COMPLETE")
