#!/usr/bin/env python3
"""
KEGG_A3_Neighborhood_Diagnostic.py
====================================

Comprehensive KEGG enrichment diagnostic for TCGA bulk and all three
single-cell networks. For each network, this script:

  1. Runs KEGG pathway enrichment (via Enrichr) on every community
     with >= 5 genes
  2. Identifies the A3A and A3B neighborhood within their community:
     - Positively correlated genes (co-activated with A3 in DIFF)
     - Negatively correlated genes (anti-correlated with A3 in DIFF)
  3. Runs KEGG on the positive and negative A3 sub-networks separately
  4. Reports the top correlated and anti-correlated genes with A3

Output (to scripts/TROUBLESHOOTING/KEGG_DIAGNOSTIC/):
  {network}/community_KEGG_summary.tsv
  {network}/A3_neighborhood_genes.tsv
  {network}/A3_positive_KEGG.tsv
  {network}/A3_negative_KEGG.tsv
  {network}/diagnostic_report.txt

Usage:
  conda run -n NETWORK python KEGG_A3_Neighborhood_Diagnostic.py

Author: Jake Lehle / Texas Biomedical Research Institute
"""

import os
import sys
import json
import time
import pickle
import numpy as np
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
OUTPUT_ROOT = os.path.join(BASE_DIR, "scripts/TROUBLESHOOTING/KEGG_DIAGNOSTIC")

# Rate limit for Enrichr API (seconds between calls)
ENRICHR_DELAY = 16

# Minimum genes for enrichment
MIN_GENES_ENRICHR = 5

# A3 gene identifiers
A3_ENSG = {
    "ENSG00000128383": "A3A",
    "ENSG00000179750": "A3B",
}
A3_SYMBOLS = {
    "APOBEC3A": "A3A",
    "APOBEC3B": "A3B",
}

# Networks to analyze
NETWORKS = [
    {
        "name": "TCGA-HNSC",
        "partition": os.path.join(BASE_DIR, "data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv"),
        "corr_pkl": os.path.join(BASE_DIR, "data/FIG_2/04_correlation_networks/TCGA-HNSC/corr_matrices/TCGA-HNSC_corr_DIFF.pkl"),
        "ensg_map": os.path.join(BASE_DIR, "data/FIG_2/01_cleaned_expression/ensg_to_symbol.json"),
        "id_type": "ensg",
        "a3_ids": A3_ENSG,
    },
    {
        "name": "NETWORK_SBS2_VS_CNV",
        "partition": os.path.join(BASE_DIR, "data/FIG_4/NETWORK_SBS2_VS_CNV/04_communities/SC_best_partition.csv"),
        "corr_pkl": os.path.join(BASE_DIR, "data/FIG_4/NETWORK_SBS2_VS_CNV/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl"),
        "ensg_map": None,
        "id_type": "symbol",
        "a3_ids": A3_SYMBOLS,
    },
    {
        "name": "NETWORK_SBS2_VS_NORMAL",
        "partition": os.path.join(BASE_DIR, "data/FIG_4/NETWORK_SBS2_VS_NORMAL/04_communities/SC_best_partition.csv"),
        "corr_pkl": os.path.join(BASE_DIR, "data/FIG_4/NETWORK_SBS2_VS_NORMAL/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl"),
        "ensg_map": None,
        "id_type": "symbol",
        "a3_ids": A3_SYMBOLS,
    },
    {
        "name": "NETWORK_CNV_VS_NORMAL",
        "partition": os.path.join(BASE_DIR, "data/FIG_4/NETWORK_CNV_VS_NORMAL/04_communities/SC_best_partition.csv"),
        "corr_pkl": os.path.join(BASE_DIR, "data/FIG_4/NETWORK_CNV_VS_NORMAL/03_correlation_networks/corr_matrices/SC_corr_DIFF.pkl"),
        "ensg_map": None,
        "id_type": "symbol",
        "a3_ids": A3_SYMBOLS,
    },
]


# =============================================================================
# LOGGING
# =============================================================================

report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="=", width=80):
    log("")
    log(char * width)
    log(f"  {title}")
    log(char * width)


# =============================================================================
# ENRICHR WRAPPER
# =============================================================================

def run_enrichr_kegg(gene_list, description="query", delay=ENRICHR_DELAY):
    """
    Run KEGG pathway enrichment via gseapy/Enrichr.

    Returns a DataFrame with enrichment results, or empty DataFrame on failure.
    """
    import gseapy as gp

    if len(gene_list) < MIN_GENES_ENRICHR:
        return pd.DataFrame()

    try:
        time.sleep(delay)
        enr = gp.enrichr(
            gene_list=list(gene_list),
            gene_sets="KEGG_2021_Human",
            organism="human",
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        df = enr.results
        if len(df) > 0:
            df = df.sort_values("Adjusted P-value").reset_index(drop=True)
        return df
    except Exception as e:
        log(f"    Enrichr error: {e}")
        return pd.DataFrame()


def summarize_kegg(df, max_show=5):
    """Return summary string of top KEGG hits."""
    if len(df) == 0:
        return "No results"

    sig = df[df["Adjusted P-value"] < 0.05]
    n_sig = len(sig)

    lines = []
    lines.append(f"{n_sig} significant (FDR<0.05) of {len(df)} tested")
    show = sig.head(max_show) if n_sig > 0 else df.head(3)
    for _, row in show.iterrows():
        term = row.get("Term", "?")
        pval = row.get("Adjusted P-value", 1.0)
        genes = row.get("Genes", "")
        lines.append(f"    {term} (adj.p={pval:.2e}, genes: {genes})")

    return "\n".join(lines)


# =============================================================================
# GENE ID CONVERSION
# =============================================================================

def to_symbols(gene_list, ensg_map):
    """Convert ENSG IDs to gene symbols. Pass-through if already symbols."""
    if ensg_map is None:
        return list(gene_list)
    return [ensg_map.get(g, g) for g in gene_list]


# =============================================================================
# A3 NEIGHBORHOOD ANALYSIS
# =============================================================================

def analyze_a3_neighborhood(a3_gene_id, a3_alias, corr_matrix, community_genes,
                            all_partition_genes, ensg_map, out_dir):
    """
    Analyze the positive and negative correlation neighborhood of an A3 gene.

    Parameters
    ----------
    a3_gene_id : str
        Gene ID (ENSG or symbol) for the A3 gene.
    a3_alias : str
        Short alias (A3A or A3B).
    corr_matrix : pd.DataFrame
        DIFF correlation matrix.
    community_genes : list
        Genes in the same community as the A3 gene.
    all_partition_genes : set
        All genes in the partition (for checking matrix membership).
    ensg_map : dict or None
        ENSG to symbol mapping (None for SC).
    out_dir : str
        Output directory for this network.

    Returns
    -------
    dict with neighborhood stats and KEGG results.
    """
    banner(f"A3 Neighborhood: {a3_alias} ({a3_gene_id})", char="-")

    if a3_gene_id not in corr_matrix.index:
        log(f"  {a3_alias} not in correlation matrix. Skipping.")
        return None

    # Get correlations with all genes in the community
    comm_in_matrix = [g for g in community_genes
                      if g in corr_matrix.index and g != a3_gene_id]

    if len(comm_in_matrix) == 0:
        log(f"  No community genes in correlation matrix. Skipping.")
        return None

    corrs = corr_matrix.loc[a3_gene_id, comm_in_matrix].sort_values(
        ascending=False)

    # Also get correlations with ALL genes in the matrix (not just community)
    all_in_matrix = [g for g in corr_matrix.index if g != a3_gene_id]
    all_corrs = corr_matrix.loc[a3_gene_id, all_in_matrix]

    # Split into positive and negative within community
    pos_genes = corrs[corrs > 0].index.tolist()
    neg_genes = corrs[corrs < 0].index.tolist()

    log(f"  Community size: {len(community_genes)}")
    log(f"  In correlation matrix: {len(comm_in_matrix)}")
    log(f"  Positively correlated: {len(pos_genes)}")
    log(f"  Negatively correlated: {len(neg_genes)}")

    # Top 20 positively correlated
    log(f"\n  Top 20 positively correlated (co-activated):")
    pos_sorted = corrs[corrs > 0].head(20)
    pos_rows = []
    for gene, rho in pos_sorted.items():
        symbol = to_symbols([gene], ensg_map)[0]
        log(f"    {symbol:15s} rho={rho:+.4f}")
        pos_rows.append({"gene_id": gene, "symbol": symbol,
                          "diff_rho": round(float(rho), 4),
                          "direction": "positive"})

    # Top 20 negatively correlated
    log(f"\n  Top 20 negatively correlated (anti-correlated):")
    neg_sorted = corrs[corrs < 0].sort_values().head(20)
    neg_rows = []
    for gene, rho in neg_sorted.items():
        symbol = to_symbols([gene], ensg_map)[0]
        log(f"    {symbol:15s} rho={rho:+.4f}")
        neg_rows.append({"gene_id": gene, "symbol": symbol,
                          "diff_rho": round(float(rho), 4),
                          "direction": "negative"})

    # Save full neighborhood gene list
    all_rows = []
    for gene, rho in corrs.items():
        symbol = to_symbols([gene], ensg_map)[0]
        all_rows.append({
            "gene_id": gene, "symbol": symbol,
            "diff_rho_with_" + a3_alias: round(float(rho), 4),
            "direction": "positive" if rho > 0 else "negative",
        })

    neighborhood_df = pd.DataFrame(all_rows)
    neighborhood_path = os.path.join(
        out_dir, f"{a3_alias}_neighborhood_genes.tsv")
    neighborhood_df.to_csv(neighborhood_path, sep="\t", index=False)
    log(f"\n  [SAVE] {a3_alias} neighborhood -> {neighborhood_path}")

    # ---- KEGG on positive sub-network ----
    pos_symbols = to_symbols(pos_genes, ensg_map)
    # Remove any that are still ENSG (unmapped)
    pos_symbols_clean = [s for s in pos_symbols if not s.startswith("ENSG")]

    log(f"\n  KEGG enrichment: {a3_alias} POSITIVE sub-network "
        f"({len(pos_symbols_clean)} genes)")
    pos_kegg = run_enrichr_kegg(
        pos_symbols_clean,
        description=f"{a3_alias}_positive"
    )
    log(f"  {summarize_kegg(pos_kegg)}")

    if len(pos_kegg) > 0:
        pos_kegg_path = os.path.join(
            out_dir, f"{a3_alias}_positive_KEGG.tsv")
        pos_kegg.to_csv(pos_kegg_path, sep="\t", index=False)
        log(f"  [SAVE] -> {pos_kegg_path}")

    # ---- KEGG on negative sub-network ----
    neg_symbols = to_symbols(neg_genes, ensg_map)
    neg_symbols_clean = [s for s in neg_symbols if not s.startswith("ENSG")]

    log(f"\n  KEGG enrichment: {a3_alias} NEGATIVE sub-network "
        f"({len(neg_symbols_clean)} genes)")
    neg_kegg = run_enrichr_kegg(
        neg_symbols_clean,
        description=f"{a3_alias}_negative"
    )
    log(f"  {summarize_kegg(neg_kegg)}")

    if len(neg_kegg) > 0:
        neg_kegg_path = os.path.join(
            out_dir, f"{a3_alias}_negative_KEGG.tsv")
        neg_kegg.to_csv(neg_kegg_path, sep="\t", index=False)
        log(f"  [SAVE] -> {neg_kegg_path}")

    # ---- Also check: A3 correlations with genes OUTSIDE its community ----
    partition_genes = set(community_genes)
    outside_in_matrix = [g for g in all_in_matrix
                         if g not in partition_genes
                         and g in set(all_partition_genes)]
    if len(outside_in_matrix) > 0:
        outside_corrs = corr_matrix.loc[a3_gene_id, outside_in_matrix]
        top_outside_pos = outside_corrs.sort_values(ascending=False).head(10)
        top_outside_neg = outside_corrs.sort_values(ascending=True).head(10)

        log(f"\n  Top 10 positively correlated OUTSIDE {a3_alias}'s community:")
        for gene, rho in top_outside_pos.items():
            symbol = to_symbols([gene], ensg_map)[0]
            log(f"    {symbol:15s} rho={rho:+.4f}")

        log(f"\n  Top 10 negatively correlated OUTSIDE {a3_alias}'s community:")
        for gene, rho in top_outside_neg.items():
            symbol = to_symbols([gene], ensg_map)[0]
            log(f"    {symbol:15s} rho={rho:+.4f}")

    return {
        "n_positive": len(pos_genes),
        "n_negative": len(neg_genes),
        "pos_kegg_sig": len(pos_kegg[pos_kegg["Adjusted P-value"] < 0.05])
            if len(pos_kegg) > 0 else 0,
        "neg_kegg_sig": len(neg_kegg[neg_kegg["Adjusted P-value"] < 0.05])
            if len(neg_kegg) > 0 else 0,
    }


# =============================================================================
# PROCESS ONE NETWORK
# =============================================================================

def process_network(net_config):
    """Run full KEGG diagnostic for one network."""
    name = net_config["name"]
    banner(f"NETWORK: {name}")

    out_dir = os.path.join(OUTPUT_ROOT, name)
    os.makedirs(out_dir, exist_ok=True)

    # ---- Load partition ----
    if not os.path.exists(net_config["partition"]):
        log(f"  Partition not found: {net_config['partition']}")
        log(f"  SKIPPING {name}")
        return
    part_df = pd.read_csv(net_config["partition"])

    # Identify columns
    gene_col = "gene"
    comm_col = "community"
    for col in ["gene", "gene_symbol", "Gene"]:
        if col in part_df.columns:
            gene_col = col
            break
    for col in ["community", "Community"]:
        if col in part_df.columns:
            comm_col = col
            break

    gene_to_comm = dict(zip(part_df[gene_col], part_df[comm_col]))
    comm_to_genes = {}
    for g, c in gene_to_comm.items():
        comm_to_genes.setdefault(c, []).append(g)

    log(f"  Partition: {len(gene_to_comm)} genes, "
        f"{len(comm_to_genes)} communities")

    # ---- Load ENSG map if needed ----
    ensg_map = None
    if net_config["ensg_map"] and os.path.exists(net_config["ensg_map"]):
        with open(net_config["ensg_map"]) as f:
            ensg_map = json.load(f)
        log(f"  ENSG map: {len(ensg_map)} entries")

    # ---- Load DIFF correlation matrix ----
    if not os.path.exists(net_config["corr_pkl"]):
        log(f"  Correlation matrix not found: {net_config['corr_pkl']}")
        log(f"  SKIPPING {name}")
        return

    log(f"  Loading DIFF correlation matrix...")
    with open(net_config["corr_pkl"], "rb") as f:
        corr_matrix = pickle.load(f)
    log(f"  DIFF matrix shape: {corr_matrix.shape}")

    # =========================================================================
    # PART 1: Community-level KEGG enrichment
    # =========================================================================
    banner(f"COMMUNITY KEGG ENRICHMENT: {name}")

    comm_kegg_rows = []
    sorted_comms = sorted(comm_to_genes.keys(),
                          key=lambda c: len(comm_to_genes[c]), reverse=True)

    for cid in sorted_comms:
        genes = comm_to_genes[cid]
        symbols = to_symbols(genes, ensg_map)
        symbols_clean = [s for s in symbols if not s.startswith("ENSG")]

        # Check if any A3 genes are in this community
        a3_here = []
        for a3_id, a3_alias in net_config["a3_ids"].items():
            if a3_id in genes:
                a3_here.append(a3_alias)

        a3_tag = f" [A3: {', '.join(a3_here)}]" if a3_here else ""

        log(f"\n  Community {cid}: {len(genes)} genes "
            f"({len(symbols_clean)} mapped){a3_tag}")

        if len(symbols_clean) < MIN_GENES_ENRICHR:
            log(f"    [SKIP] Too few genes for enrichment")
            comm_kegg_rows.append({
                "community": cid, "n_genes": len(genes),
                "n_mapped": len(symbols_clean),
                "a3_genes": ", ".join(a3_here),
                "n_sig_kegg": 0,
                "top_term": "N/A", "top_adj_p": 1.0,
                "top_genes": "",
            })
            continue

        kegg = run_enrichr_kegg(symbols_clean,
                                description=f"C{cid}_{name}")
        log(f"    {summarize_kegg(kegg)}")

        # Save per-community KEGG
        if len(kegg) > 0:
            comm_kegg_path = os.path.join(
                out_dir, f"community_{cid}_KEGG.tsv")
            kegg.to_csv(comm_kegg_path, sep="\t", index=False)

            sig = kegg[kegg["Adjusted P-value"] < 0.05]
            top = kegg.iloc[0]
            comm_kegg_rows.append({
                "community": cid, "n_genes": len(genes),
                "n_mapped": len(symbols_clean),
                "a3_genes": ", ".join(a3_here),
                "n_sig_kegg": len(sig),
                "top_term": top.get("Term", "?"),
                "top_adj_p": float(top.get("Adjusted P-value", 1.0)),
                "top_genes": top.get("Genes", ""),
            })
        else:
            comm_kegg_rows.append({
                "community": cid, "n_genes": len(genes),
                "n_mapped": len(symbols_clean),
                "a3_genes": ", ".join(a3_here),
                "n_sig_kegg": 0,
                "top_term": "No results", "top_adj_p": 1.0,
                "top_genes": "",
            })

    # Save community summary
    comm_summary_df = pd.DataFrame(comm_kegg_rows)
    comm_summary_path = os.path.join(out_dir, "community_KEGG_summary.tsv")
    comm_summary_df.to_csv(comm_summary_path, sep="\t", index=False)
    log(f"\n  [SAVE] Community KEGG summary -> {comm_summary_path}")

    # =========================================================================
    # PART 2: A3 neighborhood analysis
    # =========================================================================
    banner(f"A3 NEIGHBORHOOD ANALYSIS: {name}")

    all_partition_genes = set(gene_to_comm.keys())

    for a3_id, a3_alias in net_config["a3_ids"].items():
        if a3_id not in gene_to_comm:
            log(f"\n  {a3_alias} ({a3_id}): NOT in partition. Skipping.")
            continue

        a3_comm = gene_to_comm[a3_id]
        comm_genes = comm_to_genes[a3_comm]
        log(f"\n  {a3_alias}: Community {a3_comm} ({len(comm_genes)} genes)")

        analyze_a3_neighborhood(
            a3_gene_id=a3_id,
            a3_alias=a3_alias,
            corr_matrix=corr_matrix,
            community_genes=comm_genes,
            all_partition_genes=all_partition_genes,
            ensg_map=ensg_map,
            out_dir=out_dir,
        )

    # Free memory
    del corr_matrix

    log(f"\n  {name} COMPLETE")


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("KEGG + A3 NEIGHBORHOOD DIAGNOSTIC")
    log(f"Start: {t0}")
    log(f"Output: {OUTPUT_ROOT}")

    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    # Check if specific network requested via command line
    if len(sys.argv) > 1:
        requested = sys.argv[1]
        matching = [n for n in NETWORKS if n["name"] == requested]
        if matching:
            networks_to_run = matching
            log(f"Running single network: {requested}")
        else:
            log(f"Unknown network: {requested}")
            log(f"Available: {[n['name'] for n in NETWORKS]}")
            return
    else:
        networks_to_run = NETWORKS
        log(f"Running all {len(NETWORKS)} networks")

    for net in networks_to_run:
        try:
            process_network(net)
        except Exception as e:
            log(f"\n  ERROR processing {net['name']}: {e}")
            import traceback
            log(traceback.format_exc())
            continue

    # =========================================================================
    # Cross-network summary
    # =========================================================================
    if len(networks_to_run) > 1:
        banner("CROSS-NETWORK SUMMARY")

        for net in networks_to_run:
            summary_path = os.path.join(
                OUTPUT_ROOT, net["name"], "community_KEGG_summary.tsv")
            if os.path.exists(summary_path):
                df = pd.read_csv(summary_path, sep="\t")
                n_with_sig = len(df[df["n_sig_kegg"] > 0])
                total_sig = df["n_sig_kegg"].sum()
                log(f"\n  {net['name']}:")
                log(f"    Communities analyzed: {len(df)}")
                log(f"    With significant KEGG: {n_with_sig}")
                log(f"    Total significant terms: {total_sig}")

                # A3 community enrichment
                a3_rows = df[df["a3_genes"].str.len() > 0]
                if len(a3_rows) > 0:
                    for _, row in a3_rows.iterrows():
                        log(f"    A3 community (C{int(row['community'])}): "
                            f"{int(row['n_sig_kegg'])} sig KEGG, "
                            f"top: {row['top_term']} "
                            f"(p={row['top_adj_p']:.2e})")

    # Save full report
    elapsed = datetime.now() - t0
    banner(f"DIAGNOSTIC COMPLETE | Elapsed: {elapsed}")

    report_path = os.path.join(OUTPUT_ROOT, "diagnostic_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
    log(f"[SAVE] Full report -> {report_path}")


if __name__ == "__main__":
    main()
