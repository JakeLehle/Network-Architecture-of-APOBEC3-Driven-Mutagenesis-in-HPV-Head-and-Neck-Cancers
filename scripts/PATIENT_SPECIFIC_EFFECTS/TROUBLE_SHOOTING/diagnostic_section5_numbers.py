#!/usr/bin/env python3
"""
diagnostic_section5_numbers.py
================================

Diagnostic script to verify all numbers for Section 1.5 (Figure 5):
Patient-Specific Contributions to the SBS2-HIGH Population.

Pulls and verifies numbers for:
  BEAT 1: Patient contribution breakdown (Fig 5A-B)
    - Per-patient total basal cells and SBS2-HIGH cells
    - Chi-square test for over-representation
    - High-contributor identification (SC029, SC013, SC001)
    - Per-patient mean A3A and A3B expression
    - Tissue type breakdown (tumor vs normal in SBS2-HIGH)

  BEAT 2: Shared transcriptional program (Fig 5C) + LOPO (Supp Fig 6)
    - Silhouette score from PCA
    - LOPO: A3 wall integrity, activating chain recovery,
      gene overlap Jaccard, community ARI

  BEAT 3: HC-exclusive variant analysis (Fig 5D-E)
    - Total SComatic variants
    - SNP tier breakdown (universal, broadly shared, etc.)
    - HC-exclusive gene count and KEGG enrichment
    - Network overlap of HC-exclusive HPV/immune genes
    - Confirmation: zero A3 genes, zero activating chain genes
      in HC-exclusive variants

Usage:
    conda run -n NETWORK python diagnostic_section5_numbers.py
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2_contingency
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input data
ADATA_PATH = os.path.join(BASE_DIR, "data/FIG_4/00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/three_group_assignments.tsv")

# SComatic mutations
SCOMATIC_PATH = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"

# Figure 5 output directories
FIG5_ROOT = os.path.join(BASE_DIR, "data/FIG_5")
DIR_00_DIAG = os.path.join(FIG5_ROOT, "00_diagnostics")
DIR_01_EXPR = os.path.join(FIG5_ROOT, "01_patient_expression")
DIR_02_SNP = os.path.join(FIG5_ROOT, "02_snp_haplotype")
DIR_03_SENS = os.path.join(FIG5_ROOT, "03_sensitivity")

# SC network partition (for HC-exclusive overlap check)
PARTITION_PATH = os.path.join(BASE_DIR, "data/FIG_4/NETWORK_SBS2_VS_NORMAL/04_communities/SC_best_partition.csv")

# Column names
PATIENT_COL = "subject id"
TISSUE_COL = "tissue type"
CELLTYPE_COL = "final_annotation"

# Reference gene sets
A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
ACTIVATING_CHAIN = {"RALY", "HNRNPA2B1", "CCL20", "KRT24", "LCN2",
                    "LINC00278", "RRAD", "SMOX", "UTY"}
HIGH_CONTRIBUTORS = ["Patient SC029", "Patient SC013", "Patient SC001"]


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# BEAT 1: PATIENT CONTRIBUTION (Fig 5A-B)
# =============================================================================

def beat1_patient_contributions(adata, sbs2_high_cells):
    """Per-patient basal cell counts, SBS2-HIGH contributions, A3 expression."""
    banner("BEAT 1: PATIENT CONTRIBUTION BREAKDOWN (Fig 5A-B)")

    # Get basal cells
    basal_mask = adata.obs[CELLTYPE_COL].str.contains("basal", case=False, na=False)
    basal = adata[basal_mask].copy()
    log(f"Total basal cells: {basal.n_obs:,}")

    # Map SBS2-HIGH to basal
    high_in_basal = set(sbs2_high_cells) & set(basal.obs_names)
    basal.obs["is_SBS2_HIGH"] = basal.obs_names.isin(high_in_basal)
    log(f"SBS2-HIGH cells in basal: {len(high_in_basal)}")

    # Per-patient breakdown
    patients = sorted(basal.obs[PATIENT_COL].unique())
    log(f"Unique patients: {len(patients)}")

    rows = []
    for p in patients:
        p_mask = basal.obs[PATIENT_COL] == p
        p_basal = basal[p_mask]
        n_total = p_mask.sum()
        n_high = p_basal.obs["is_SBS2_HIGH"].sum()
        n_tumor = (p_basal.obs[TISSUE_COL] == "tumor").sum()
        n_normal = (p_basal.obs[TISSUE_COL] == "normal").sum()

        # A3 expression (across ALL basal cells for this patient)
        a3a_expr = None
        a3b_expr = None
        if "APOBEC3A" in basal.var_names:
            x = p_basal[:, "APOBEC3A"].X
            a3a_expr = float(x.toarray().mean() if hasattr(x, "toarray") else x.mean())
        if "APOBEC3B" in basal.var_names:
            x = p_basal[:, "APOBEC3B"].X
            a3b_expr = float(x.toarray().mean() if hasattr(x, "toarray") else x.mean())

        expected = n_total / basal.n_obs * len(high_in_basal)
        fold = n_high / expected if expected > 0 else 0

        rows.append({
            "patient": p,
            "n_basal_total": n_total,
            "n_SBS2_HIGH": n_high,
            "pct_of_HIGH": n_high / len(high_in_basal) * 100 if high_in_basal else 0,
            "expected": expected,
            "fold_enrichment": fold,
            "n_tumor": n_tumor,
            "n_normal": n_normal,
            "mean_A3A": a3a_expr,
            "mean_A3B": a3b_expr,
        })

    df = pd.DataFrame(rows).sort_values("n_SBS2_HIGH", ascending=False)

    log("\n  Per-patient breakdown (sorted by SBS2-HIGH contribution):")
    log(f"  {'Patient':<20} {'Basal':>8} {'HIGH':>6} {'%HIGH':>7} "
        f"{'Fold':>6} {'A3A':>6} {'A3B':>6} {'Tumor':>7} {'Normal':>7}")
    log(f"  {'-'*80}")
    for _, r in df.iterrows():
        log(f"  {r['patient']:<20} {r['n_basal_total']:>8,} {r['n_SBS2_HIGH']:>6} "
            f"{r['pct_of_HIGH']:>6.1f}% {r['fold_enrichment']:>5.1f}x "
            f"{r['mean_A3A']:>6.2f} {r['mean_A3B']:>6.2f} "
            f"{r['n_tumor']:>7,} {r['n_normal']:>7,}")

    # Top 3 contributors
    top3 = df.head(3)
    top3_pct = top3["pct_of_HIGH"].sum()
    log(f"\n  Top 3 contributors: {', '.join(top3['patient'].values)}")
    log(f"  Combined: {top3_pct:.1f}% of SBS2-HIGH")

    # Zero contributors
    zero = df[df["n_SBS2_HIGH"] == 0]
    if len(zero) > 0:
        log(f"  Zero-contributor patients: {', '.join(zero['patient'].values)}")

    # Tissue type breakdown of SBS2-HIGH
    high_basal = basal[basal.obs["is_SBS2_HIGH"]]
    n_high_tumor = (high_basal.obs[TISSUE_COL] == "tumor").sum()
    n_high_normal = (high_basal.obs[TISSUE_COL] == "normal").sum()
    pct_tumor = n_high_tumor / len(high_in_basal) * 100 if high_in_basal else 0
    log(f"\n  SBS2-HIGH tissue breakdown: {n_high_tumor} tumor ({pct_tumor:.1f}%), "
        f"{n_high_normal} normal ({100-pct_tumor:.1f}%)")

    # Chi-square test for non-uniform patient contribution
    observed = df["n_SBS2_HIGH"].values
    expected_uniform = np.full(len(observed), len(high_in_basal) / len(patients))
    # Use contingency table approach: HIGH vs non-HIGH per patient
    n_non_high = df["n_basal_total"].values - df["n_SBS2_HIGH"].values
    contingency = np.array([df["n_SBS2_HIGH"].values, n_non_high])
    chi2, p_chi, dof, _ = chi2_contingency(contingency)
    log(f"\n  Chi-square test (patient x HIGH/non-HIGH):")
    log(f"    chi2 = {chi2:.2f}, p = {p_chi:.2e}, dof = {dof}")

    # A3 expression ranking context
    log(f"\n  A3 expression context:")
    df_sorted_a3a = df.sort_values("mean_A3A", ascending=False)
    log(f"  Highest mean A3A: {df_sorted_a3a.iloc[0]['patient']} "
        f"(mean={df_sorted_a3a.iloc[0]['mean_A3A']:.2f}, "
        f"HIGH cells={int(df_sorted_a3a.iloc[0]['n_SBS2_HIGH'])}, "
        f"fold={df_sorted_a3a.iloc[0]['fold_enrichment']:.1f}x)")

    # Where do HC patients rank in A3A expression?
    a3a_rank = df_sorted_a3a.reset_index(drop=True)
    for hc in HIGH_CONTRIBUTORS:
        rank_idx = a3a_rank[a3a_rank["patient"] == hc].index
        if len(rank_idx) > 0:
            log(f"  {hc}: A3A rank {rank_idx[0]+1}/{len(patients)} "
                f"(mean={float(a3a_rank.loc[rank_idx[0], 'mean_A3A']):.2f})")

    return df


# =============================================================================
# BEAT 2: SHARED PROGRAM + LOPO (Fig 5C, Supp Fig 6)
# =============================================================================

def beat2_shared_program_and_lopo():
    """Silhouette score and LOPO network sensitivity results."""
    banner("BEAT 2: SHARED PROGRAM (Fig 5C) + LOPO (Supp Fig 6)")

    # Silhouette score
    sil_path = os.path.join(DIR_01_EXPR, "HIGH_cell_silhouette_score.txt")
    if os.path.exists(sil_path):
        with open(sil_path) as f:
            content = f.read().strip()
        log(f"  Silhouette score file contents:\n    {content}")
    else:
        log(f"  [WARNING] Silhouette file not found: {sil_path}")
        # Also check for TSV format
        alt_path = os.path.join(DIR_01_EXPR, "HIGH_cell_similarity_metrics.tsv")
        if os.path.exists(alt_path):
            log(f"  Found alternative: {alt_path}")
            metrics = pd.read_csv(alt_path, sep="\t")
            log(f"  {metrics.to_string()}")

    # LOPO results
    banner("LOPO SENSITIVITY (Supp Fig 6)", char="-")
    lopo_patients = ["SC029", "SC013", "SC001"]

    log(f"  {'Config':<12} {'Chain':>8} {'ARI':>8} {'Jaccard':>8} "
        f"{'A3A Wall':>10} {'A3B Wall':>10} {'HIGH_n':>8} {'NORM_n':>8}")
    log(f"  {'-'*76}")

    # Full analysis baseline
    log(f"  {'Full':<12} {'9/9':>8} {'1.000':>8} {'1.000':>8} "
        f"{'100%':>10} {'100%':>10} {'546':>8} {'546':>8}")

    for patient in lopo_patients:
        summary_path = os.path.join(DIR_03_SENS, f"LOPO_{patient}",
                                     f"LOPO_{patient}_summary.tsv")
        if os.path.exists(summary_path):
            r = pd.read_csv(summary_path, sep="\t").iloc[0]
            chain_rec = int(r.get("activating_chain_recovered", 0))
            chain_tot = int(r.get("activating_chain_total", 9))
            ari = float(r.get("community_ari_with_full", 0))
            jaccard = float(r.get("jaccard_with_full", 0))
            a3a_wall = r.get("a3a_wall_pct_neg", None)
            a3b_wall = r.get("a3b_wall_pct_neg", None)
            n_high = int(r.get("N_HIGH", r.get("high_remaining",
                        r.get("n_high", "?"))))
            n_norm = int(r.get("N_NORMAL", r.get("norm_remaining",
                        r.get("n_normal", "?"))))

            a3a_str = f"{a3a_wall:.0f}%" if pd.notna(a3a_wall) else "N/A"
            a3b_str = f"{a3b_wall:.0f}%" if pd.notna(a3b_wall) else "N/A"

            log(f"  -{patient:<11} {chain_rec}/{chain_tot}:>8 {ari:>8.3f} "
                f"{jaccard:>8.3f} {a3a_str:>10} {a3b_str:>10} "
                f"{n_high:>8} {n_norm:>8}")

            # Also show which chain genes were recovered
            chain_genes = r.get("activating_chain_genes", "")
            if chain_genes and str(chain_genes) != "nan":
                log(f"    Chain genes: {chain_genes}")
        else:
            log(f"  -{patient:<11} [summary file not found]")

        # Also check parameters file for network size
        params_path = os.path.join(DIR_03_SENS, f"LOPO_{patient}",
                                    f"LOPO_{patient}_parameters.txt")
        if os.path.exists(params_path):
            with open(params_path) as f:
                params = f.read().strip()
            log(f"    Parameters: {params}")


# =============================================================================
# BEAT 3: HC-EXCLUSIVE VARIANTS (Fig 5D-E)
# =============================================================================

def beat3_hc_exclusive_variants():
    """SNP tier analysis, HC-exclusive KEGG enrichment, network overlap."""
    banner("BEAT 3: HC-EXCLUSIVE VARIANT ANALYSIS (Fig 5D-E)")

    # --- Total variant count ---
    if os.path.exists(SCOMATIC_PATH):
        log(f"  Loading SComatic variants: {SCOMATIC_PATH}")
        # Read just enough to count
        mut_df = pd.read_csv(SCOMATIC_PATH, sep="\t", low_memory=False)
        log(f"  Total variant rows: {len(mut_df):,}")
        if "FILTER" in mut_df.columns:
            pass_count = (mut_df["FILTER"] == "PASS").sum()
            log(f"  PASS variants: {pass_count:,}")
    else:
        log(f"  [WARNING] SComatic file not found: {SCOMATIC_PATH}")

    # --- SNP tier breakdown ---
    banner("SNP TIER BREAKDOWN", char="-")

    # Try variant_to_gene_mapping.tsv first (has tier column)
    mapping_path = os.path.join(DIR_02_SNP, "gene_analysis",
                                "variant_to_gene_mapping.tsv")
    tier_report_path = os.path.join(DIR_02_SNP, "gene_analysis",
                                     "SNP_tier_gene_analysis_report.tsv")

    if os.path.exists(mapping_path):
        log(f"  Loading variant-to-gene mapping: {mapping_path}")
        mapping = pd.read_csv(mapping_path, sep="\t")
        log(f"  Total mapped variants: {len(mapping):,}")

        if "tier" in mapping.columns:
            tier_counts = mapping["tier"].value_counts()
            log(f"\n  Variant counts by tier:")
            for tier, count in tier_counts.items():
                log(f"    {tier}: {count:,}")

            # HC-exclusive details
            hc = mapping[mapping["tier"] == "HC-exclusive"]
            log(f"\n  HC-exclusive variants: {len(hc):,}")

            if "gene_str" in hc.columns:
                hc_genes = set()
                for gs in hc["gene_str"].dropna():
                    if gs:
                        hc_genes.update(gs.split(","))
                hc_genes.discard("")
                log(f"  HC-exclusive unique genes: {len(hc_genes)}")

                # Check for A3 genes
                a3_in_hc = hc_genes & A3_SYMBOLS
                log(f"\n  A3 genes in HC-exclusive: {len(a3_in_hc)} "
                    f"{'(' + ', '.join(a3_in_hc) + ')' if a3_in_hc else '(NONE)'}")

                # Check for activating chain genes
                chain_in_hc = hc_genes & ACTIVATING_CHAIN
                log(f"  Activating chain genes in HC-exclusive: {len(chain_in_hc)} "
                    f"{'(' + ', '.join(chain_in_hc) + ')' if chain_in_hc else '(NONE)'}")
    else:
        log(f"  [WARNING] Variant mapping not found: {mapping_path}")

    # --- Gene analysis report ---
    if os.path.exists(tier_report_path):
        log(f"\n  Loading tier gene report: {tier_report_path}")
        report = pd.read_csv(tier_report_path, sep="\t")
        log(f"  Tiers in report: {report['tier'].unique().tolist() if 'tier' in report.columns else 'N/A'}")
        if "tier" in report.columns and "n_genes" in report.columns:
            for _, r in report.iterrows():
                log(f"    {r['tier']}: {r.get('n_genes', '?')} genes, "
                    f"{r.get('n_variants', '?')} variants")

    # --- HC-exclusive KEGG enrichment ---
    banner("HC-EXCLUSIVE KEGG ENRICHMENT", char="-")
    kegg_path = os.path.join(DIR_02_SNP, "gene_analysis",
                              "SNP_tier_HC-exclusive_KEGG.tsv")
    if os.path.exists(kegg_path):
        kegg = pd.read_csv(kegg_path, sep="\t")
        log(f"  Total KEGG terms tested: {len(kegg)}")

        sig = kegg[kegg["Adjusted P-value"] < 0.05]
        log(f"  Significant (adj p < 0.05): {len(sig)}")

        # Show all significant + top non-significant
        show = sig if len(sig) > 0 else kegg.head(10)
        log(f"\n  {'Term':<50} {'Adj P':>10} {'Genes':>6}")
        log(f"  {'-'*70}")
        for _, r in show.iterrows():
            genes_raw = r.get("Genes", "")
            n_genes = len(genes_raw.split(";")) if pd.notna(genes_raw) and genes_raw else 0
            sig_marker = "*" if r["Adjusted P-value"] < 0.05 else " "
            log(f"  {sig_marker}{r['Term']:<49} {r['Adjusted P-value']:>10.4f} {n_genes:>6}")

        # Target pathways
        target_keywords = ["papilloma", "hpv", "antigen", "graft",
                           "phagosome", "cell adhesion"]
        target_rows = []
        for _, r in kegg.iterrows():
            term_lower = r["Term"].lower()
            if any(kw in term_lower for kw in target_keywords):
                target_rows.append(r)

        if target_rows:
            log(f"\n  HPV/immune-related pathways:")
            for r in target_rows:
                genes_raw = r.get("Genes", "")
                genes_list = genes_raw.split(";") if pd.notna(genes_raw) and genes_raw else []
                log(f"    {r['Term']}: adj p = {r['Adjusted P-value']:.4f}, "
                    f"{len(genes_list)} genes")
                if genes_list:
                    log(f"      Genes: {', '.join(sorted(genes_list))}")
    else:
        log(f"  [WARNING] HC-exclusive KEGG not found: {kegg_path}")

    # --- Network overlap ---
    banner("HC-EXCLUSIVE x SC NETWORK OVERLAP", char="-")
    if os.path.exists(PARTITION_PATH) and os.path.exists(mapping_path):
        # Load network genes
        part_df = pd.read_csv(PARTITION_PATH)
        col = "gene" if "gene" in part_df.columns else part_df.columns[0]
        comm_col = "community" if "community" in part_df.columns else part_df.columns[1]
        sc_network = set(part_df[col].tolist())
        gene_to_comm = dict(zip(part_df[col], part_df[comm_col]))
        log(f"  SC network genes: {len(sc_network)}")

        # HC-exclusive genes
        mapping = pd.read_csv(mapping_path, sep="\t")
        hc = mapping[mapping["tier"] == "HC-exclusive"]
        hc_genes = set()
        for gs in hc["gene_str"].dropna():
            if gs:
                hc_genes.update(gs.split(","))
        hc_genes.discard("")

        overlap = hc_genes & sc_network
        log(f"  HC-exclusive genes in SC network: {len(overlap)} / {len(hc_genes)}")

        # Check for specific target genes
        target_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "TAP1",
                        "B2M", "CDH1", "CLDN1", "F11R", "MX1", "MDM2"]
        log(f"\n  Target gene check in HC-exclusive AND network:")
        for g in target_genes:
            in_hc = g in hc_genes
            in_net = g in sc_network
            comm = gene_to_comm.get(g, "N/A")
            status = "HC+NET" if (in_hc and in_net) else \
                     "HC only" if in_hc else \
                     "NET only" if in_net else "neither"
            log(f"    {g:<12} {status:<10} community={comm}")

        # If KEGG data available, cross-reference pathway genes with network
        if os.path.exists(kegg_path):
            kegg = pd.read_csv(kegg_path, sep="\t")
            all_pathway_genes = set()
            for _, r in kegg.iterrows():
                term_lower = r["Term"].lower()
                if any(kw in term_lower for kw in target_keywords):
                    genes_raw = r.get("Genes", "")
                    if pd.notna(genes_raw) and genes_raw:
                        all_pathway_genes.update(genes_raw.split(";"))

            pathway_in_network = all_pathway_genes & sc_network
            log(f"\n  HPV/immune pathway genes from HC-exclusive KEGG:")
            log(f"    Total unique pathway genes: {len(all_pathway_genes)}")
            log(f"    In SC network: {len(pathway_in_network)}")
            if pathway_in_network:
                log(f"    Genes (with community):")
                for g in sorted(pathway_in_network):
                    comm = gene_to_comm.get(g, "?")
                    log(f"      {g}: community {comm}")
    else:
        if not os.path.exists(PARTITION_PATH):
            log(f"  [WARNING] Partition not found: {PARTITION_PATH}")
        if not os.path.exists(mapping_path):
            log(f"  [WARNING] Mapping not found: {mapping_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("DIAGNOSTIC: SECTION 1.5 (FIGURE 5) NUMBERS")
    log(f"Base directory: {BASE_DIR}")

    # Load adata
    log("\nLoading adata...")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  adata: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Load group assignments
    groups_df = pd.read_csv(GROUP_PATH, sep="\t")
    group_map = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))
    sbs2_high_cells = set(
        bc for bc, grp in group_map.items() if grp == "SBS2_HIGH")
    log(f"  SBS2-HIGH cells from assignments: {len(sbs2_high_cells)}")

    # BEAT 1: Patient contributions
    patient_df = beat1_patient_contributions(adata, sbs2_high_cells)

    # Free memory
    del adata

    # BEAT 2: Shared program + LOPO
    beat2_shared_program_and_lopo()

    # BEAT 3: HC-exclusive variants
    beat3_hc_exclusive_variants()

    banner("DIAGNOSTIC COMPLETE")


if __name__ == "__main__":
    main()
