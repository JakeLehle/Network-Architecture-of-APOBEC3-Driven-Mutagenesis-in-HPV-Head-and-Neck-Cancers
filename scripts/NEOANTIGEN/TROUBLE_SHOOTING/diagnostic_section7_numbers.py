#!/usr/bin/env python3
"""
diagnostic_section7_numbers.py
================================

Diagnostic script to verify all numbers for Section 1.7 (Figure 7):
Neoantigen Landscape and Therapeutic Target Identification.

Verifies:
  BEAT 1: Neoantigen + fusion burden (Panel A)
    - Total neoantigens per population
    - Strong binders (IC50 < 50nM)
    - Differential neoantigens (mutant binds, WT does not)
    - Per-cell rates (n=546)
    - Fold difference
    - RNA fusion counts per population

  BEAT 2: Venn overlap + tier framework (Panel B)
    - Gene counts: SBS2-only, shared, CNV-only
    - Tier counts: 1A (escaped shared), 1B (SBS2-specific),
      2 (CNV-specific), 3 (broad coverage)
    - Top genes per tier with composite scores
    - Escape mechanisms: silenced genes, multi-escape genes

  BEAT 3: ANXA1 deep dive (Panels C-D)
    - Composite score and rank
    - Best IC50, number of peptides, strong binders
    - Expression in SBS2-HIGH cells (% positive, mean)
    - Hotspot regions (pos 289, pos 321-327)
    - Presence/absence in CNV-HIGH neoantigen list

  BEAT 4: Summary stats for closing
    - Protein-altering variant counts (after germline subtraction)
    - Missense counts
    - Proteome mapping rate

Usage:
    conda run -n NETWORK python diagnostic_section7_numbers.py
"""

import os
import re
import numpy as np
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG7_ROOT = os.path.join(BASE_DIR, "data/FIG_7")

MHC_DIR = os.path.join(FIG7_ROOT, "03_mhc_binding")
FUSION_DIR = os.path.join(FIG7_ROOT, "04_fusion_analysis")
SUMMARY_DIR = os.path.join(FIG7_ROOT, "05_summary")
TIER_DIR = os.path.join(SUMMARY_DIR, "THERAPEUTIC_TIERS")
ANNOTATION_DIR = os.path.join(FIG7_ROOT, "02_snpeff_annotation")
EXPRESSION_DIR = os.path.join(FIG7_ROOT, "05_summary")

N_CELLS = 546


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# BEAT 1: NEOANTIGEN + FUSION BURDEN (Panel A)
# =============================================================================

def beat1_burden():
    """Neoantigen counts, strong binders, per-cell rates, fusions."""
    banner("BEAT 1: NEOANTIGEN + FUSION BURDEN (Panel A)")

    for group in ['SBS2_HIGH', 'CNV_HIGH']:
        neo_path = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
        all_pep_path = os.path.join(MHC_DIR, f"{group}_all_peptide_results.tsv")

        if not os.path.exists(neo_path):
            log(f"  [WARNING] {neo_path} not found")
            continue

        neo = pd.read_csv(neo_path, sep="\t")
        log(f"\n  {group}:")
        log(f"    Total predicted neoantigens: {len(neo):,}")
        log(f"    Per cell: {len(neo)/N_CELLS:.2f}")

        # Strong binders
        if 'mut_ic50' in neo.columns:
            strong = (neo['mut_ic50'] < 50).sum()
            log(f"    Strong binders (IC50 < 50nM): {strong}")
        elif 'best_ic50' in neo.columns:
            strong = (neo['best_ic50'] < 50).sum()
            log(f"    Strong binders (IC50 < 50nM): {strong}")

        # Differential neoantigens (mutant binds, WT does not)
        if 'is_differential' in neo.columns:
            diff = neo['is_differential'].sum()
            log(f"    Differential neoantigens: {diff}")
        elif 'wt_ic50' in neo.columns and 'mut_ic50' in neo.columns:
            diff = ((neo['mut_ic50'] < 500) & (neo['wt_ic50'] >= 500)).sum()
            log(f"    Differential (mut<500, wt>=500): {diff}")

        # Unique genes
        if 'gene' in neo.columns:
            n_genes = neo['gene'].nunique()
            log(f"    Unique neoantigen genes: {n_genes}")

        # All peptide results (for total peptide count)
        if os.path.exists(all_pep_path):
            all_pep = pd.read_csv(all_pep_path, sep="\t")
            log(f"    Total peptides tested: {len(all_pep):,}")
            if 'is_binder' in all_pep.columns:
                log(f"    Binders (IC50 < 500nM): {all_pep['is_binder'].sum():,}")

    # Fold difference
    sbs2_path = os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv")
    cnv_path = os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv")
    if os.path.exists(sbs2_path) and os.path.exists(cnv_path):
        n_sbs2 = len(pd.read_csv(sbs2_path, sep="\t"))
        n_cnv = len(pd.read_csv(cnv_path, sep="\t"))
        fold = n_sbs2 / n_cnv if n_cnv > 0 else 0
        log(f"\n  Fold difference (SBS2/CNV): {fold:.2f}x")
        log(f"  Per-cell: SBS2={n_sbs2/N_CELLS:.2f}, CNV={n_cnv/N_CELLS:.2f}")

    # Fusions
    banner("RNA FUSIONS", char="-")
    fusion_path = os.path.join(FUSION_DIR, "per_group_junction_summary.tsv")
    if os.path.exists(fusion_path):
        fusions = pd.read_csv(fusion_path, sep="\t")
        log(f"  Fusion summary:")
        for _, row in fusions.iterrows():
            log(f"    {row['group']}: {row['total_junctions']} junctions")
    else:
        log(f"  [WARNING] Fusion summary not found: {fusion_path}")

    # Also check for detailed fusion files
    for group in ['SBS2_HIGH', 'CNV_HIGH']:
        fusion_detail = os.path.join(FUSION_DIR, f"{group}_fusions.tsv")
        if os.path.exists(fusion_detail):
            fd = pd.read_csv(fusion_detail, sep="\t")
            log(f"    {group} fusion details: {len(fd)} entries")
            if 'gene_pair' in fd.columns or 'gene1' in fd.columns:
                log(f"    Unique gene pairs: {fd.shape[0]}")


# =============================================================================
# BEAT 2: VENN OVERLAP + TIER FRAMEWORK (Panel B)
# =============================================================================

def beat2_tiers():
    """Gene overlap, tier counts, top genes per tier."""
    banner("BEAT 2: VENN OVERLAP + TIER FRAMEWORK (Panel B)")

    sbs2_path = os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv")
    cnv_path = os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv")

    if not os.path.exists(sbs2_path) or not os.path.exists(cnv_path):
        log("  [WARNING] Neoantigen files not found")
        return

    sbs2_neo = pd.read_csv(sbs2_path, sep="\t")
    cnv_neo = pd.read_csv(cnv_path, sep="\t")

    sbs2_genes = set(sbs2_neo['gene'].dropna().unique())
    cnv_genes = set(cnv_neo['gene'].dropna().unique())
    shared = sbs2_genes & cnv_genes
    sbs2_only = sbs2_genes - cnv_genes
    cnv_only = cnv_genes - sbs2_genes

    log(f"  Venn diagram gene counts:")
    log(f"    SBS2-HIGH total genes: {len(sbs2_genes)}")
    log(f"    CNV-HIGH total genes: {len(cnv_genes)}")
    log(f"    SBS2-only: {len(sbs2_only)}")
    log(f"    Shared: {len(shared)}")
    log(f"    CNV-only: {len(cnv_only)}")

    # Tier breakdown
    banner("THERAPEUTIC TIERS", char="-")
    tier_path = os.path.join(TIER_DIR, "all_genes_tiered.tsv")
    if os.path.exists(tier_path):
        tiers = pd.read_csv(tier_path, sep="\t")
        log(f"  Total tiered genes: {len(tiers)}")

        tier_counts = tiers['tier'].value_counts()
        log(f"\n  Tier counts:")
        for tier, count in sorted(tier_counts.items()):
            log(f"    {tier}: {count}")

        # Top genes per tier (by composite score if available)
        score_col = None
        for candidate in ['composite_score', 'score', 'combined_score',
                          'sbs2_composite', 'cnv_composite']:
            if candidate in tiers.columns:
                score_col = candidate
                break

        if score_col is None:
            # Look for any numeric column that could be a score
            numeric_cols = tiers.select_dtypes(include=[np.number]).columns.tolist()
            log(f"  Available numeric columns: {numeric_cols}")
            if numeric_cols:
                score_col = numeric_cols[0]

        for tier_name in sorted(tiers['tier'].unique()):
            tier_sub = tiers[tiers['tier'] == tier_name]
            if score_col and score_col in tier_sub.columns:
                top5 = tier_sub.nlargest(5, score_col)
                log(f"\n  Top 5 {tier_name} (by {score_col}):")
                for _, row in top5.iterrows():
                    gene = row.get('gene', row.get('gene_symbol', '?'))
                    score = row[score_col]
                    log(f"    {gene}: {score:.1f}")
            else:
                log(f"\n  {tier_name}: {len(tier_sub)} genes (no score column found)")
                if 'gene' in tier_sub.columns:
                    log(f"    First 5: {', '.join(tier_sub['gene'].head(5).tolist())}")

        # Show all columns for reference
        log(f"\n  Tier file columns: {list(tiers.columns)}")
    else:
        log(f"  [WARNING] Tier file not found: {tier_path}")

    # Escape mechanisms
    banner("ESCAPE MECHANISMS", char="-")

    # Expression silencing
    silence_path = os.path.join(TIER_DIR, "expression_silenced_genes.tsv")
    if os.path.exists(silence_path):
        silenced = pd.read_csv(silence_path, sep="\t")
        log(f"  Expression-silenced genes (>2x down in CNV): {len(silenced)}")
        if 'gene' in silenced.columns:
            log(f"    Examples: {', '.join(silenced['gene'].head(10).tolist())}")
    else:
        # Check alternative paths
        for alt in ["silenced_genes.tsv", "escape_silencing.tsv"]:
            alt_path = os.path.join(TIER_DIR, alt)
            if os.path.exists(alt_path):
                silenced = pd.read_csv(alt_path, sep="\t")
                log(f"  Found at {alt}: {len(silenced)} genes")
                break

    # Multi-escape genes
    multi_path = os.path.join(TIER_DIR, "multi_escape_genes.tsv")
    if os.path.exists(multi_path):
        multi = pd.read_csv(multi_path, sep="\t")
        log(f"  Multi-escape genes: {len(multi)}")
        if 'gene' in multi.columns:
            log(f"    Genes: {', '.join(multi['gene'].tolist())}")
    else:
        for alt in ["multi_mechanism_escape.tsv", "convergent_escape.tsv"]:
            alt_path = os.path.join(TIER_DIR, alt)
            if os.path.exists(alt_path):
                multi = pd.read_csv(alt_path, sep="\t")
                log(f"  Found at {alt}: {len(multi)} genes")
                break

    # List all files in TIER_DIR for reference
    if os.path.exists(TIER_DIR):
        log(f"\n  All files in TIER_DIR:")
        for f in sorted(os.listdir(TIER_DIR)):
            fpath = os.path.join(TIER_DIR, f)
            if os.path.isfile(fpath):
                size = os.path.getsize(fpath)
                log(f"    {f} ({size:,} bytes)")


# =============================================================================
# BEAT 3: ANXA1 DEEP DIVE (Panels C-D)
# =============================================================================

def beat3_anxa1():
    """ANXA1 specifics: score, IC50, peptides, hotspots, expression."""
    banner("BEAT 3: ANXA1 DEEP DIVE (Panels C-D)")

    # ANXA1 in neoantigen list
    sbs2_neo = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv"), sep="\t")
    cnv_neo = pd.read_csv(os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv"), sep="\t")

    anxa1_sbs2 = sbs2_neo[sbs2_neo['gene'] == 'ANXA1']
    anxa1_cnv = cnv_neo[cnv_neo['gene'] == 'ANXA1']

    log(f"  ANXA1 in SBS2-HIGH neoantigens: {len(anxa1_sbs2)} peptides")
    log(f"  ANXA1 in CNV-HIGH neoantigens: {len(anxa1_cnv)} peptides")
    log(f"  ANXA1 absent from CNV-HIGH: {len(anxa1_cnv) == 0}")

    # All peptide results for ANXA1
    all_pep_path = os.path.join(MHC_DIR, "SBS2_HIGH_all_peptide_results.tsv")
    if os.path.exists(all_pep_path):
        all_pep = pd.read_csv(all_pep_path, sep="\t")
        anxa1_pep = all_pep[all_pep['gene'] == 'ANXA1']
        log(f"\n  ANXA1 all peptides tested: {len(anxa1_pep)}")

        if 'is_binder' in anxa1_pep.columns:
            n_binders = anxa1_pep['is_binder'].sum()
            log(f"  ANXA1 binders (IC50 < 500nM): {n_binders}")

        if 'mut_ic50' in anxa1_pep.columns:
            strong = (anxa1_pep['mut_ic50'] < 50).sum()
            best_ic50 = anxa1_pep['mut_ic50'].min()
            log(f"  ANXA1 strong binders (IC50 < 50nM): {strong}")
            log(f"  ANXA1 best IC50: {best_ic50:.1f} nM")

        # Hotspot regions
        if 'mut_pos_protein' in anxa1_pep.columns:
            binders = anxa1_pep[anxa1_pep['is_binder'] == True] if 'is_binder' in anxa1_pep.columns else anxa1_pep[anxa1_pep['mut_ic50'] < 500]

            r1 = binders[binders['mut_pos_protein'] == 289]
            r2 = binders[binders['mut_pos_protein'].isin([321, 322, 326, 327])]
            other = binders[~binders['mut_pos_protein'].isin([289, 321, 322, 326, 327])]

            log(f"\n  Hotspot regions (among binders):")
            log(f"    Pos 289: {len(r1)} peptides")
            log(f"    Pos 321-327: {len(r2)} peptides")
            log(f"    Other positions: {len(other)} peptides")

            # Unique positions
            all_mut_pos = binders['mut_pos_protein'].unique()
            log(f"    Unique mutation positions with binders: {sorted(all_mut_pos)}")

    # Composite score
    tier_path = os.path.join(TIER_DIR, "all_genes_tiered.tsv")
    if os.path.exists(tier_path):
        tiers = pd.read_csv(tier_path, sep="\t")
        anxa1_tier = tiers[tiers['gene'] == 'ANXA1'] if 'gene' in tiers.columns else pd.DataFrame()
        if len(anxa1_tier) > 0:
            log(f"\n  ANXA1 tier assignment:")
            for col in anxa1_tier.columns:
                val = anxa1_tier.iloc[0][col]
                log(f"    {col}: {val}")

    # ANXA1 expression in SBS2-HIGH
    banner("ANXA1 EXPRESSION", char="-")

    # Check for expression data in summary
    expr_path = os.path.join(SUMMARY_DIR, "neoantigen_expression_scores.tsv")
    if os.path.exists(expr_path):
        expr = pd.read_csv(expr_path, sep="\t")
        anxa1_expr = expr[expr['gene'] == 'ANXA1'] if 'gene' in expr.columns else pd.DataFrame()
        if len(anxa1_expr) > 0:
            log(f"  ANXA1 expression data:")
            for col in anxa1_expr.columns:
                val = anxa1_expr.iloc[0][col]
                log(f"    {col}: {val}")
    else:
        log(f"  [INFO] Expression scores file not found: {expr_path}")
        # Try alternative
        for alt in ["gene_expression_summary.tsv", "expression_weighted_scores.tsv"]:
            alt_path = os.path.join(SUMMARY_DIR, alt)
            if os.path.exists(alt_path):
                expr = pd.read_csv(alt_path, sep="\t")
                anxa1_expr = expr[expr.iloc[:, 0] == 'ANXA1']
                if len(anxa1_expr) > 0:
                    log(f"  Found in {alt}:")
                    log(f"    {anxa1_expr.to_string()}")
                break

    # ANXA1 somatic variants
    banner("ANXA1 SOMATIC VARIANTS", char="-")
    spa_path = os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv")
    if os.path.exists(spa_path):
        spa = pd.read_csv(spa_path, sep="\t")
        anxa1_var = spa[spa['gene'] == 'ANXA1']
        log(f"  ANXA1 protein-altering variants: {len(anxa1_var)}")

        if 'hgvs_p' in anxa1_var.columns:
            positions = []
            for _, row in anxa1_var.iterrows():
                hgvs = str(row['hgvs_p'])
                match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
                if match:
                    pos = int(match.group(2))
                    positions.append(pos)
                    log(f"    {hgvs} (position {pos})")
            log(f"  Unique positions: {sorted(set(positions))}")


# =============================================================================
# BEAT 4: VARIANT SUMMARY STATS
# =============================================================================

def beat4_summary():
    """Protein-altering variant counts, missense counts, mapping rate."""
    banner("BEAT 4: VARIANT SUMMARY STATS")

    for group in ['SBS2_HIGH', 'CNV_HIGH']:
        spa_path = os.path.join(ANNOTATION_DIR, f"{group}.somatic_protein_altering.tsv")
        if os.path.exists(spa_path):
            spa = pd.read_csv(spa_path, sep="\t")
            log(f"\n  {group}:")
            log(f"    Total protein-altering variants: {len(spa)}")
            log(f"    Per cell: {len(spa)/N_CELLS:.2f}")

            if 'effect' in spa.columns:
                effect_counts = spa['effect'].value_counts()
                log(f"    By effect type:")
                for effect, count in effect_counts.items():
                    log(f"      {effect}: {count}")
            elif 'annotation' in spa.columns:
                ann_counts = spa['annotation'].value_counts()
                log(f"    By annotation:")
                for ann, count in ann_counts.items():
                    log(f"      {ann}: {count}")

            if 'gene' in spa.columns:
                n_genes = spa['gene'].nunique()
                log(f"    Unique genes: {n_genes}")
        else:
            log(f"  [WARNING] {spa_path} not found")

    # Check for germline subtraction stats
    for f in os.listdir(ANNOTATION_DIR):
        if 'germline' in f.lower() or 'subtract' in f.lower():
            log(f"\n  Germline-related file: {f}")

    # List all files in annotation dir for reference
    log(f"\n  All files in ANNOTATION_DIR:")
    if os.path.exists(ANNOTATION_DIR):
        for f in sorted(os.listdir(ANNOTATION_DIR)):
            fpath = os.path.join(ANNOTATION_DIR, f)
            if os.path.isfile(fpath):
                log(f"    {f}")

    # Proteome mapping stats
    banner("PROTEOME MAPPING STATS", char="-")
    map_path = os.path.join(FIG7_ROOT, "03_mhc_binding/proteome_mapping_stats.tsv")
    if os.path.exists(map_path):
        stats = pd.read_csv(map_path, sep="\t")
        log(f"  Proteome mapping stats:")
        log(f"    {stats.to_string()}")
    else:
        # Check alternative locations
        for search_dir in [MHC_DIR, SUMMARY_DIR, FIG7_ROOT]:
            for f in os.listdir(search_dir):
                if 'mapping' in f.lower() or 'proteome' in f.lower():
                    log(f"  Found: {os.path.join(search_dir, f)}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("DIAGNOSTIC: SECTION 1.7 (FIGURE 7) NUMBERS")
    log(f"Base directory: {BASE_DIR}")

    # List FIG_7 directory structure
    log(f"\nFIG_7 directory structure:")
    for d in sorted(os.listdir(FIG7_ROOT)):
        dpath = os.path.join(FIG7_ROOT, d)
        if os.path.isdir(dpath):
            n_files = len([f for f in os.listdir(dpath) if os.path.isfile(os.path.join(dpath, f))])
            log(f"  {d}/ ({n_files} files)")

    beat1_burden()
    beat2_tiers()
    beat3_anxa1()
    beat4_summary()

    banner("DIAGNOSTIC COMPLETE")


if __name__ == "__main__":
    main()
