#!/usr/bin/env python3
"""
Diagnostic_HC_Network_HPV_Overlap.py
=====================================

Quick check: do HC-exclusive variant genes that overlap with the SC network
specifically include HPV and antigen presentation pathway genes?

Reads:
  - HC-exclusive KEGG enrichment results (from Analyze_SNP_Tier_Genes.py)
  - SC network partition (from Figure 4 pipeline)
  - Variant-to-gene mapping (from Analyze_SNP_Tier_Genes.py)

Usage:
  conda run -n NETWORK python Diagnostic_HC_Network_HPV_Overlap.py
"""

import os, pandas as pd

from patient_config import (
    DIR_02_SNP, FIGURE_5_PANELS, COMMUNITIES_DIR,
    banner, log, ensure_dir,
)

def main():
    banner("HC-EXCLUSIVE x NETWORK x HPV/ANTIGEN OVERLAP")

    # Paths
    kegg_path = os.path.join(DIR_02_SNP, "gene_analysis",
                             "SNP_tier_HC-exclusive_KEGG.tsv")
    partition_path = os.path.join(COMMUNITIES_DIR, "SC_best_partition.csv")
    mapping_path = os.path.join(DIR_02_SNP, "gene_analysis",
                                "variant_to_gene_mapping.tsv")

    # Load SC network genes
    banner("LOADING SC NETWORK")
    sc_df = pd.read_csv(partition_path)
    col = 'gene' if 'gene' in sc_df.columns else sc_df.columns[0]
    comm_col = 'community' if 'community' in sc_df.columns else sc_df.columns[1]
    sc_network = set(sc_df[col].tolist())
    gene_to_comm = dict(zip(sc_df[col], sc_df[comm_col]))
    log(f"  SC network genes: {len(sc_network)}")

    # Load HC-exclusive variant genes
    banner("LOADING HC-EXCLUSIVE VARIANT GENES")
    mapping = pd.read_csv(mapping_path, sep='\t')
    hc = mapping[mapping['tier'] == 'HC-exclusive']
    hc_genes = set()
    for gs in hc['gene_str'].dropna():
        if gs:
            hc_genes.update(gs.split(','))
    hc_genes.discard('')
    log(f"  HC-exclusive genes: {len(hc_genes)}")

    hc_in_network = hc_genes & sc_network
    log(f"  HC-exclusive genes IN network: {len(hc_in_network)}")

    # Load KEGG results for HC-exclusive
    banner("HC-EXCLUSIVE KEGG PATHWAY GENES")
    if not os.path.exists(kegg_path):
        log(f"  ERROR: {kegg_path} not found")
        return

    kegg = pd.read_csv(kegg_path, sep='\t')
    log(f"  Total KEGG terms: {len(kegg)}")

    # Target pathways: HPV, antigen, immune-related
    target_terms = []
    for _, row in kegg.iterrows():
        term = row['Term'].lower()
        if any(kw in term for kw in ['papilloma', 'hpv', 'antigen',
                                      'graft', 'immune', 'phagosome',
                                      'cell adhesion']):
            target_terms.append(row)

    log(f"\n  Target pathway terms found: {len(target_terms)}")

    all_pathway_genes = set()
    pathway_network_overlap = {}

    for row in target_terms:
        term = row['Term']
        pval = row['Adjusted P-value']
        # Enrichr format: genes are semicolon-separated in 'Genes' column
        genes_raw = row.get('Genes', '')
        if pd.isna(genes_raw) or not genes_raw:
            log(f"\n  {term} (p={pval:.4f}): no gene list available")
            continue

        genes = set(genes_raw.split(';'))
        in_network = genes & sc_network
        all_pathway_genes.update(genes)

        log(f"\n  {term} (p={pval:.4f})")
        log(f"    HC-exclusive genes in pathway: {len(genes)}")
        log(f"    Of those, also in SC network: {len(in_network)}")
        if in_network:
            for g in sorted(in_network):
                comm = gene_to_comm.get(g, '?')
                log(f"      {g} (community {comm})")

        pathway_network_overlap[term] = {
            'pathway_genes': genes,
            'in_network': in_network,
        }

    # Summary: all HPV/immune pathway genes that are also in the network
    banner("SUMMARY: HPV/IMMUNE PATHWAY GENES IN SC NETWORK")
    all_overlap = all_pathway_genes & sc_network
    log(f"  Total unique HPV/immune pathway genes from HC-exclusive: "
        f"{len(all_pathway_genes)}")
    log(f"  Of those in SC network: {len(all_overlap)}")

    if all_overlap:
        log(f"\n  Gene list (with community assignment):")
        for g in sorted(all_overlap):
            comm = gene_to_comm.get(g, '?')
            log(f"    {g}: community {comm}")

    # Also check: what fraction of the 199 HC-in-network genes are
    # in HPV/immune pathways vs other pathways?
    banner("CONTEXT: HC-IN-NETWORK PATHWAY BREAKDOWN")
    hc_network_in_hpv = hc_in_network & all_pathway_genes
    hc_network_not_hpv = hc_in_network - all_pathway_genes
    log(f"  HC genes in network AND HPV/immune pathways: "
        f"{len(hc_network_in_hpv)}")
    log(f"  HC genes in network but NOT HPV/immune: "
        f"{len(hc_network_not_hpv)}")

    # Check ALL significant KEGG terms for the HC-in-network subset
    sig_kegg = kegg[kegg['Adjusted P-value'] < 0.05]
    log(f"\n  All significant KEGG terms ({len(sig_kegg)}):")
    for _, row in sig_kegg.iterrows():
        genes_raw = row.get('Genes', '')
        if pd.isna(genes_raw) or not genes_raw:
            continue
        genes = set(genes_raw.split(';'))
        in_net = genes & sc_network
        if in_net:
            log(f"    {row['Term']} (p={row['Adjusted P-value']:.4f}): "
                f"{len(in_net)}/{len(genes)} pathway genes in network")

    log("\nDiagnostic complete.")


if __name__ == "__main__":
    main()
