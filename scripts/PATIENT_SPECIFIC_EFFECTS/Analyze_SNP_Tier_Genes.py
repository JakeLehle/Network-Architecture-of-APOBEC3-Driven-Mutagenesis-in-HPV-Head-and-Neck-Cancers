#!/usr/bin/env python3
"""
Analyze_SNP_Tier_Genes.py
===========================

For each variant sharing tier (Universal, Broadly shared, HC-exclusive),
map variants to genes and check:
  1. Is the gene in the SC differential co-expression network (Figure 4)?
  2. Is it a known A3 interactor (Harris/McCann/Jang)?
  3. Is it an A3 enzyme?
  4. GSEA (KEGG) on the gene list per tier
  5. Cross-reference with neoantigen databases (if provided)

Requires a GTF file to map variant positions to genes. Set GTF_PATH below.

This script is designed to be reusable for Figure 6 (high-CNV population).

Usage:
  conda run -n NETWORK python Analyze_SNP_Tier_Genes.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os, sys, gzip, numpy as np, pandas as pd
from collections import Counter, defaultdict
import gseapy as gp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

from patient_config import *

# =============================================================================
# CONFIGURATION — SET THESE PATHS
# =============================================================================

# Cell Ranger reference GTF (GRCh38)
# Common locations — update to match your cluster setup
GTF_CANDIDATES = [
    "/master/jlehle/WORKING/SC/ref/GRCh38/genes/genes_unzipped.gtf"
]

# Variant sharing tiers (output from Generate_Supplemental script)
TIER_FILE = os.path.join(FIGURE_5_PANELS, "variant_sharing_tiers.tsv")

# SC network community genes
SC_PARTITION_FILE = os.path.join(COMMUNITIES_DIR, "SC_best_partition.csv")

# Neoantigen databases (optional — set to None if not available)
TSNADB_PATH = None   # e.g., "/path/to/TSNAdb_v2.0.tsv"
NEPDB_PATH  = None   # e.g., "/path/to/NEPdb.tsv"

# Tiers to analyze (skip patient-specific and partially shared)
TIERS_TO_ANALYZE = ['Universal', 'Broadly shared', 'HC-exclusive']


def parse_gtf(gtf_path):
    """
    Parse a GTF file into a dict of chrom -> list of (start, end, gene_name).
    Handles both .gtf and .gtf.gz files.
    """
    banner("PARSING GTF")
    log(f"  File: {gtf_path}")

    gene_coords = defaultdict(list)
    opener = gzip.open if gtf_path.endswith('.gz') else open
    n_genes = 0

    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'gene':
                continue

            chrom = fields[0]
            # Ensure chr prefix
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom

            start = int(fields[3])
            end = int(fields[4])

            # Parse gene_name from attributes
            attrs = fields[8]
            gene_name = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1]
                    break

            if gene_name:
                gene_coords[chrom].append((start, end, gene_name))
                n_genes += 1

    # Sort each chromosome's genes by start position for binary search
    for chrom in gene_coords:
        gene_coords[chrom].sort(key=lambda x: x[0])

    log(f"  Parsed {n_genes:,} gene entries across {len(gene_coords)} chromosomes")
    return gene_coords


def map_variant_to_gene(chrom, pos, gene_coords):
    """Map a single variant position to overlapping gene(s)."""
    if chrom not in gene_coords:
        return []
    genes = []
    for start, end, name in gene_coords[chrom]:
        if start > pos + 10000:  # past our position with buffer
            break
        if start <= pos <= end:
            genes.append(name)
    return genes


def map_all_variants(tier_df, gene_coords):
    """Map all variants to genes. Returns tier_df with added 'genes' column."""
    banner("MAPPING VARIANTS TO GENES")

    gene_lists = []
    n_mapped = 0
    n_total = len(tier_df)

    for i, row in tier_df.iterrows():
        parts = row['variant_id'].split(':')
        if len(parts) < 2:
            gene_lists.append([]); continue
        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            gene_lists.append([]); continue

        genes = map_variant_to_gene(chrom, pos, gene_coords)
        gene_lists.append(genes)
        if genes:
            n_mapped += 1

        if (i + 1) % 10000 == 0:
            log(f"  Processed {i+1:,}/{n_total:,}...")

    tier_df['genes'] = gene_lists
    tier_df['gene_str'] = tier_df['genes'].apply(lambda x: ','.join(x) if x else '')
    log(f"  Mapped {n_mapped:,}/{n_total:,} variants to genes ({100*n_mapped/n_total:.1f}%)")
    return tier_df


def cross_reference_analysis(tier_df, sc_network_genes, harris_all, harris_a3b, a3_genes):
    """Cross-reference per-tier gene lists against network, interactors, A3 enzymes."""
    banner("CROSS-REFERENCE ANALYSIS")

    for tier in TIERS_TO_ANALYZE:
        sub = tier_df[tier_df['tier'] == tier]
        tier_genes = set()
        for gl in sub['genes']:
            tier_genes.update(gl)

        if len(tier_genes) == 0:
            log(f"\n  {tier}: no genes mapped"); continue

        # Cross-reference
        in_network = tier_genes & sc_network_genes
        in_harris = tier_genes & harris_all
        in_harris_a3b = tier_genes & harris_a3b
        in_a3 = tier_genes & a3_genes

        log(f"\n  === {tier} ({len(sub)} variants, {len(tier_genes)} unique genes) ===")
        log(f"    In SC network (Fig 4): {len(in_network)} genes")
        if in_network:
            for g in sorted(in_network)[:20]:
                log(f"      {g}")
            if len(in_network) > 20:
                log(f"      ... and {len(in_network)-20} more")

        log(f"    Known A3 interactors (Harris): {len(in_harris)} genes")
        if in_harris:
            for g in sorted(in_harris): log(f"      {g}")

        log(f"    A3B-specific interactors: {len(in_harris_a3b)} genes")
        if in_harris_a3b:
            for g in sorted(in_harris_a3b): log(f"      {g}")

        log(f"    A3 enzymes: {len(in_a3)} genes")
        if in_a3:
            for g in sorted(in_a3): log(f"      {g}")

    return tier_df


def run_gsea_per_tier(tier_df, out_dir):
    """Run enrichr (KEGG) on gene lists per tier."""
    banner("GSEA PER TIER (KEGG)")

    gsea_results = {}
    for tier in TIERS_TO_ANALYZE:
        sub = tier_df[tier_df['tier'] == tier]
        tier_genes = set()
        for gl in sub['genes']:
            tier_genes.update(gl)

        tier_genes = sorted(tier_genes)
        if len(tier_genes) < 5:
            log(f"  {tier}: only {len(tier_genes)} genes, skipping GSEA"); continue

        log(f"  {tier}: running enrichr on {len(tier_genes)} genes...")

        try:
            enr = gp.enrichr(
                gene_list=tier_genes,
                gene_sets='KEGG_2021_Human',
                no_plot=True,
                verbose=False,
            )
            res = enr.results
            res['tier'] = tier

            sig = res[res['Adjusted P-value'] < 0.05]
            log(f"    {len(res)} pathways tested, {len(sig)} significant (adj.p < 0.05)")

            if len(sig) > 0:
                top = sig.sort_values('Adjusted P-value').head(10)
                for _, r in top.iterrows():
                    log(f"      {r['Term']}: p={r['Adjusted P-value']:.4f}, "
                        f"genes={r['Overlap']}")

            gsea_results[tier] = res
            res.to_csv(os.path.join(out_dir, f"SNP_tier_{tier.replace(' ','_')}_KEGG.tsv"),
                      sep='\t', index=False)

        except Exception as e:
            log(f"    ERROR: {e}")

    return gsea_results


def neoantigen_cross_reference(tier_df, out_dir):
    """Cross-reference variants with neoantigen databases (if available)."""
    banner("NEOANTIGEN DATABASE CROSS-REFERENCE")

    for name, path in [('TSNAdb', TSNADB_PATH), ('NEPdb', NEPDB_PATH)]:
        if path is None or not os.path.exists(path):
            log(f"  {name}: not provided or not found (path={path})")
            log(f"  To enable, set {name.upper()}_PATH in the script config.")
            continue

        log(f"  Loading {name}: {path}")
        try:
            db = pd.read_csv(path, sep='\t')
            log(f"    Shape: {db.shape}")
            log(f"    Columns: {list(db.columns)[:10]}")

            # The cross-reference logic depends on database format
            # Typical fields: gene, mutation, chromosome, position
            # Implement specific matching once database format is confirmed
            log(f"    NOTE: Cross-reference logic needs customization for {name} format.")
            log(f"    Saving gene overlap for manual inspection.")

            # Basic gene-level overlap
            for tier in TIERS_TO_ANALYZE:
                sub = tier_df[tier_df['tier'] == tier]
                tier_genes = set()
                for gl in sub['genes']:
                    tier_genes.update(gl)

                # Try common gene column names
                for gene_col in ['gene', 'Gene', 'gene_name', 'Hugo_Symbol']:
                    if gene_col in db.columns:
                        db_genes = set(db[gene_col].dropna().unique())
                        overlap = tier_genes & db_genes
                        log(f"    {tier}: {len(overlap)} genes overlap with {name}")
                        if overlap:
                            for g in sorted(overlap)[:10]:
                                log(f"      {g}")
                        break
        except Exception as e:
            log(f"    ERROR loading {name}: {e}")


def generate_report(tier_df, gsea_results, out_dir):
    """Save comprehensive report TSV."""
    banner("GENERATING REPORT")

    # Per-tier summary
    summary_rows = []
    for tier in TIERS_TO_ANALYZE:
        sub = tier_df[tier_df['tier'] == tier]
        tier_genes = set()
        for gl in sub['genes']:
            tier_genes.update(gl)

        row = {
            'tier': tier,
            'n_variants': len(sub),
            'n_genes': len(tier_genes),
            'top_genes': ','.join(sorted(tier_genes)[:20]),
        }

        # Add top GSEA pathway if available
        if tier in gsea_results:
            sig = gsea_results[tier][gsea_results[tier]['Adjusted P-value'] < 0.05]
            if len(sig) > 0:
                row['top_kegg'] = sig.iloc[0]['Term']
                row['top_kegg_pval'] = sig.iloc[0]['Adjusted P-value']
                row['n_sig_kegg'] = len(sig)
            else:
                row['top_kegg'] = 'none'; row['top_kegg_pval'] = 1.0; row['n_sig_kegg'] = 0
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    report_path = os.path.join(out_dir, "SNP_tier_gene_analysis_report.tsv")
    summary_df.to_csv(report_path, sep='\t', index=False)
    log(f"  Saved: {report_path}")

    # Full variant-gene mapping
    full_path = os.path.join(out_dir, "variant_to_gene_mapping.tsv")
    tier_df[['variant_id', 'tier', 'n_patients', 'patients', 'gene_str']].to_csv(
        full_path, sep='\t', index=False)
    log(f"  Saved: {full_path}")


def main():
    banner("SNP TIER GENE ANALYSIS")
    out_dir = ensure_dir(os.path.join(DIR_02_SNP, "gene_analysis"))

    # ── Load tier assignments ────────────────────────────────────────────
    if not os.path.exists(TIER_FILE):
        log(f"  ERROR: Tier file not found: {TIER_FILE}")
        log(f"  Run Generate_Supplemental_Patient_Effects_v3.py first.")
        sys.exit(1)

    tier_df = pd.read_csv(TIER_FILE, sep='\t')
    log(f"  Loaded {len(tier_df):,} variants with tier assignments")
    for t in TIERS_TO_ANALYZE:
        log(f"    {t}: {(tier_df['tier']==t).sum():,}")

    # ── Find and parse GTF ───────────────────────────────────────────────
    gtf_path = None
    for candidate in GTF_CANDIDATES:
        if os.path.exists(candidate):
            gtf_path = candidate; break
        # Also check .gz version
        if os.path.exists(candidate + '.gz'):
            gtf_path = candidate + '.gz'; break

    if gtf_path is None:
        log("  WARNING: GTF file not found at any candidate path:")
        for c in GTF_CANDIDATES:
            log(f"    {c}")
        log("  Cannot map variants to genes. Exiting.")
        log("  Update GTF_CANDIDATES in the script with the correct path.")
        sys.exit(1)

    gene_coords = parse_gtf(gtf_path)

    # ── Map variants to genes ────────────────────────────────────────────
    tier_df = map_all_variants(tier_df, gene_coords)

    # ── Load reference gene sets ─────────────────────────────────────────
    banner("LOADING REFERENCE GENE SETS")

    # SC network genes
    sc_network_genes = set()
    if os.path.exists(SC_PARTITION_FILE):
        sc_df = pd.read_csv(SC_PARTITION_FILE)
        col = 'gene' if 'gene' in sc_df.columns else sc_df.columns[0]
        sc_network_genes = set(sc_df[col].tolist())
        log(f"  SC network genes: {len(sc_network_genes)}")
    else:
        log(f"  SC partition file not found: {SC_PARTITION_FILE}")

    # Harris interactors
    harris_all_genes = set()
    if os.path.exists(HARRIS_ALL_PATH):
        with open(HARRIS_ALL_PATH) as f:
            harris_all_genes = {line.strip().split('\t')[0] for line in f
                                if line.strip() and not line.startswith('#')}
        log(f"  Harris ALL interactors: {len(harris_all_genes)}")

    harris_a3b_genes = set()
    if os.path.exists(HARRIS_A3B_PATH):
        with open(HARRIS_A3B_PATH) as f:
            harris_a3b_genes = {line.strip().split('\t')[0] for line in f
                                if line.strip() and not line.startswith('#')}
        log(f"  Harris A3B interactors: {len(harris_a3b_genes)}")

    a3_gene_set = set(A3_GENES_SYMBOLS)
    log(f"  A3 enzymes: {len(a3_gene_set)}")

    # ── Cross-reference ──────────────────────────────────────────────────
    cross_reference_analysis(tier_df, sc_network_genes, harris_all_genes,
                            harris_a3b_genes, a3_gene_set)

    # ── GSEA per tier ────────────────────────────────────────────────────
    gsea_results = run_gsea_per_tier(tier_df, out_dir)

    # ── Neoantigen cross-reference ───────────────────────────────────────
    neoantigen_cross_reference(tier_df, out_dir)

    # ── Report ───────────────────────────────────────────────────────────
    generate_report(tier_df, gsea_results, out_dir)

    log("\nSNP tier gene analysis complete.")


if __name__ == "__main__":
    main()
