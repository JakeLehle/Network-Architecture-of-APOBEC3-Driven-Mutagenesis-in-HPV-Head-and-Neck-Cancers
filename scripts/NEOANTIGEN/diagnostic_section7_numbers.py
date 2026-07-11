#!/usr/bin/env python3
"""
diagnostic_section7_numbers.py
================================

Diagnostic script to verify all numbers for Section 4.5 (Figure 7):
Neoantigen Landscape and Therapeutic Target Identification.

This version adds four validation items requested for the manuscript update:
  (A) Per-cell RNA fusion rates (not just totals)
  (B) ANXA1 TCW / APOBEC trinucleotide-context check (ported from
      Diagnostic_ANXA1_TCW_Context.py; SComatic Start is matched directly to
      the annotation pos, i.e. Start == pos after the Step01 Start+1 fix)
  (C) Tier 1A escape-mechanism breakdown (antigen loss / silencing /
      fusion, plus multi-mechanism counts) over the 22 shared+escaped genes
  (D) Explicit mapping of the on-disk 4-tier classifier onto the 3-tier
      manuscript taxonomy, with MDK's tier location made explicit, plus a
      reconciliation check on the small tier4_dual_neo_fusion bucket

BEAT 5 additionally reports TCW-context prevalence under TWO substitution
scopes so both can be kept on hand for the manuscript:
  - C>T and C>G  (SBS2 + SBS13; the current manuscript text)
  - C>T only     (SBS2 arm; if the C>G/SBS13 arm is dropped)

Verifies:
  BEAT 1: Neoantigen + fusion burden (Panel A)
  BEAT 2: Venn overlap + tier framework (Panel B), mapped to 3 text tiers
  BEAT 3: ANXA1 deep dive (Panels C-D), incl. TCW context
  BEAT 4: Summary stats for closing
  BEAT 5: TCW-context prevalence across groups (both substitution scopes)

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

# SComatic master (source of REF_TRI for the TCW check)
SCOMATIC_TSV = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
                "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

N_CELLS = 546

# Classifier-tier -> manuscript-tier map (3-tier text taxonomy)
#   Text Tier 1 (105 shared) = 1A_hot_shared_escaped (22) + 3_broad_coverage (83)
#   Text Tier 2 (276 SBS2)   = 1B_hot_sbs2_specific
#   Text Tier 3 (135 CNV)    = 2_cold_cnv_specific
TEXT_TIER_MAP = {
    '1A_hot_shared_escaped': 'Text Tier 1 (shared) - escaped subset',
    '3_broad_coverage':      'Text Tier 1 (shared) - broad/no-escape subset',
    '1B_hot_sbs2_specific':  'Text Tier 2 (SBS2-specific)',
    '2_cold_cnv_specific':   'Text Tier 3 (CNV-specific)',
}


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

        if 'mut_ic50' in neo.columns:
            strong = (neo['mut_ic50'] < 50).sum()
            log(f"    Strong binders (IC50 < 50nM): {strong}")
        elif 'best_ic50' in neo.columns:
            strong = (neo['best_ic50'] < 50).sum()
            log(f"    Strong binders (IC50 < 50nM): {strong}")

        if 'is_differential' in neo.columns:
            diff = neo['is_differential'].sum()
            log(f"    Differential neoantigens: {diff}")
        elif 'wt_ic50' in neo.columns and 'mut_ic50' in neo.columns:
            diff = ((neo['mut_ic50'] < 500) & (neo['wt_ic50'] >= 500)).sum()
            log(f"    Differential (mut<500, wt>=500): {diff}")

        if 'gene' in neo.columns:
            n_genes = neo['gene'].nunique()
            log(f"    Unique neoantigen genes: {n_genes}")

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

    # -------------------------------------------------------------------------
    # Fusions: totals AND per-cell rates  [ADD A]
    # -------------------------------------------------------------------------
    banner("RNA FUSIONS (totals + per-cell rates)", char="-")
    fusion_path = os.path.join(FUSION_DIR, "per_group_junction_summary.tsv")
    if os.path.exists(fusion_path):
        fusions = pd.read_csv(fusion_path, sep="\t")
        log(f"  Fusion summary (group: total junctions, per-cell rate):")
        # Identify the count column robustly
        count_col = None
        for c in ['total_junctions', 'n_junctions', 'junctions', 'count']:
            if c in fusions.columns:
                count_col = c
                break
        for _, row in fusions.iterrows():
            grp = row['group']
            if count_col is not None:
                total = int(row[count_col])
                # per-cell: NORMAL also has 546 cells in this design
                rate = total / N_CELLS
                log(f"    {grp:12s}: {total:6,d} junctions   ({rate:.1f} per cell)")
            else:
                log(f"    {grp}: {dict(row)}")
        log(f"  (per-cell uses n={N_CELLS} cells per group)")
    else:
        log(f"  [WARNING] Fusion summary not found: {fusion_path}")
        # Fall back to recomputing from all_filtered_junctions.tsv if present
        afj = os.path.join(FUSION_DIR, "all_filtered_junctions.tsv")
        if os.path.exists(afj):
            j = pd.read_csv(afj, sep="\t")
            if 'group' in j.columns:
                log(f"  Recomputed from all_filtered_junctions.tsv:")
                for grp, n in j['group'].value_counts().items():
                    log(f"    {grp:12s}: {n:6,d} junctions   ({n/N_CELLS:.1f} per cell)")


# =============================================================================
# BEAT 2: VENN OVERLAP + TIER FRAMEWORK (Panel B)
# =============================================================================

def beat2_tiers():
    """Gene overlap, tier counts, 3-tier mapping, escape-mechanism breakdown."""
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
    banner("THERAPEUTIC TIERS (classifier, on disk)", char="-")
    tier_path = os.path.join(TIER_DIR, "all_genes_tiered.tsv")
    if not os.path.exists(tier_path):
        log(f"  [WARNING] Tier file not found: {tier_path}")
        return

    tiers = pd.read_csv(tier_path, sep="\t")
    log(f"  Total tiered genes: {len(tiers)}")

    tier_counts = tiers['tier'].value_counts()
    log(f"\n  Classifier tier counts:")
    for tier, count in sorted(tier_counts.items()):
        log(f"    {tier}: {count}")

    # Pick a usable composite score column
    score_col = None
    for candidate in ['sbs2_composite', 'composite_score', 'score',
                      'combined_score', 'cnv_composite', 'sort_score']:
        if candidate in tiers.columns:
            score_col = candidate
            break

    for tier_name in sorted(tiers['tier'].unique()):
        tier_sub = tiers[tiers['tier'] == tier_name]
        if score_col and score_col in tier_sub.columns:
            top5 = tier_sub.nlargest(5, score_col)
            log(f"\n  Top 5 {tier_name} (by {score_col}):")
            for _, row in top5.iterrows():
                gene = row.get('gene', '?')
                score = row[score_col]
                log(f"    {gene}: {score:.1f}" if pd.notna(score) else f"    {gene}: nan")
        else:
            log(f"\n  {tier_name}: {len(tier_sub)} genes (no score column)")

    log(f"\n  Tier file columns: {list(tiers.columns)}")

    # -------------------------------------------------------------------------
    # 3-tier manuscript mapping  [ADD D]
    # -------------------------------------------------------------------------
    banner("MANUSCRIPT 3-TIER MAPPING", char="-")
    log("  Classifier -> manuscript text tier:")
    for ctier, label in TEXT_TIER_MAP.items():
        n = int((tiers['tier'] == ctier).sum())
        log(f"    {ctier:24s} -> {label}  (n={n})")

    n_1a = int((tiers['tier'] == '1A_hot_shared_escaped').sum())
    n_3 = int((tiers['tier'] == '3_broad_coverage').sum())
    n_1b = int((tiers['tier'] == '1B_hot_sbs2_specific').sum())
    n_2 = int((tiers['tier'] == '2_cold_cnv_specific').sum())

    log(f"\n  Collapsed manuscript tiers:")
    log(f"    Text Tier 1 (shared)        = {n_1a} escaped + {n_3} broad = {n_1a + n_3}")
    log(f"    Text Tier 2 (SBS2-specific) = {n_1b}")
    log(f"    Text Tier 3 (CNV-specific)  = {n_2}")
    log(f"    SUM = {n_1a + n_3 + n_1b + n_2}")

    # Where does MDK live?
    if 'gene' in tiers.columns:
        mdk = tiers[tiers['gene'] == 'MDK']
        if len(mdk) > 0:
            mrow = mdk.iloc[0]
            mtier = mrow['tier']
            mscore = mrow[score_col] if score_col in mrow else float('nan')
            log(f"\n  MDK classifier tier: {mtier}  ({TEXT_TIER_MAP.get(mtier, '?')})")
            log(f"  MDK composite: {mscore:.1f}" if pd.notna(mscore) else "  MDK composite: nan")
            # Is MDK the top of the broad (non-escaped shared) subset?
            broad = tiers[tiers['tier'] == '3_broad_coverage']
            if score_col in broad.columns and len(broad) > 0:
                top_broad = broad.nlargest(1, score_col).iloc[0]['gene']
                log(f"  Top of 3_broad_coverage by {score_col}: {top_broad}")
                log(f"  -> Text claim 'Tier 1 led by MDK' refers to the broad/no-escape "
                    f"subset: {'CONFIRMED' if top_broad == 'MDK' else 'CHECK (top is ' + str(top_broad) + ')'}")

    # -------------------------------------------------------------------------
    # Tier 1A escape-mechanism breakdown  [ADD C]
    # -------------------------------------------------------------------------
    banner("TIER 1A ESCAPE-MECHANISM BREAKDOWN (22 shared+escaped genes)", char="-")
    t1a = tiers[tiers['tier'] == '1A_hot_shared_escaped'].copy()
    log(f"  Tier 1A gene count: {len(t1a)}")

    # Mechanism flag columns present in the tier file
    mech_cols = {
        'antigen_loss': 'antigen_loss_cnv',
        'expression_silenced': 'expression_silenced',
        'fusion_in_cnv': 'fusion_in_cnv',
        'cross_sbs2_neo_cnv_fusion': 'cross_sbs2_neo_cnv_fusion',
    }
    present = {k: v for k, v in mech_cols.items() if v in t1a.columns}
    log(f"  Mechanism columns present: {list(present.values())}")

    # Per-mechanism counts within Tier 1A
    for name, col in present.items():
        try:
            n = int(t1a[col].astype(bool).sum())
        except Exception:
            n = int((t1a[col] == True).sum())
        log(f"    {name:30s} ({col}): {n}")

    # Build a per-gene mechanism list and multi-mechanism count
    # Treat antigen loss, expression silencing, and ANY fusion disruption
    # (fusion_in_cnv OR cross_sbs2_neo_cnv_fusion) as the three text mechanisms.
    def gene_mechs(row):
        mechs = []
        if 'antigen_loss_cnv' in row and bool(row['antigen_loss_cnv']):
            mechs.append('antigen_loss')
        if 'expression_silenced' in row and bool(row['expression_silenced']):
            mechs.append('expression_silenced')
        fus = False
        if 'fusion_in_cnv' in row and bool(row['fusion_in_cnv']):
            fus = True
        if 'cross_sbs2_neo_cnv_fusion' in row and bool(row['cross_sbs2_neo_cnv_fusion']):
            fus = True
        if fus:
            mechs.append('fusion_disrupted')
        return mechs

    multi = 0
    log(f"\n  Per-gene mechanisms (Tier 1A):")
    score_for_sort = score_col if score_col in t1a.columns else None
    t1a_sorted = t1a.sort_values(score_for_sort, ascending=False) if score_for_sort else t1a
    for _, row in t1a_sorted.iterrows():
        mechs = gene_mechs(row)
        if len(mechs) >= 2:
            multi += 1
        gene = row.get('gene', '?')
        sc = row[score_for_sort] if score_for_sort else float('nan')
        sc_str = f"{sc:.1f}" if pd.notna(sc) else "nan"
        log(f"    {gene:14s} score={sc_str:>7s}  {'; '.join(mechs) if mechs else '(none)'}")

    log(f"\n  Tier 1A genes with >=2 escape mechanisms: {multi}")
    log(f"  Tier 1A genes with any escape mechanism: {len(t1a)} (by definition of 1A)")

    # -------------------------------------------------------------------------
    # tier4 reconciliation  [ADD D, part 2]
    # -------------------------------------------------------------------------
    banner("TIER4 RECONCILIATION (dual_neo_fusion overlay)", char="-")
    t4_path = os.path.join(TIER_DIR, "tier4_dual_neo_fusion.tsv")
    if os.path.exists(t4_path):
        t4 = pd.read_csv(t4_path, sep="\t")
        log(f"  tier4_dual_neo_fusion rows: {len(t4)}")
        if 'gene' in t4.columns:
            t4_genes = list(t4['gene'])
            log(f"  Genes: {', '.join(map(str, t4_genes))}")
            # Are these genes ALSO present in the main 4 tiers (i.e. overlay, not new)?
            main_genes = set(tiers['gene'])
            overlap = [g for g in t4_genes if g in main_genes]
            log(f"  Of these, already counted in main tiers: {len(overlap)}/{len(t4_genes)}")
            log(f"  -> tier4 is an {'OVERLAY (no double-count risk)' if len(overlap)==len(t4_genes) else 'INDEPENDENT bucket (CHECK accounting)'}")
    else:
        log(f"  tier4_dual_neo_fusion.tsv not found (ok if classifier did not emit it)")

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

# ANXA1 domain boundaries (for context labelling in the TCW table)
ANXA1_DOMAINS = [
    (1, 41, 'N-term'), (42, 109, 'Repeat 1'), (110, 122, 'Linker 1-2'),
    (123, 196, 'Repeat 2'), (197, 220, 'Linker 2-3'),
    (221, 287, 'Repeat 3'), (288, 290, 'Linker 3-4'), (291, 346, 'Repeat 4'),
]
_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def classify_tcw(ref, alt, tri):
    """Strand-aware APOBEC TCW classification.
    Returns (is_tcw, is_c_class) where is_c_class flags C>T/C>G (incl. G>A/G>C
    on the minus strand), i.e. the conventional APOBEC-eligible substitutions.
    """
    if tri is None or len(tri) != 3:
        return False, False
    if ref == 'C' and alt in ('T', 'G'):
        return (tri[0] == 'T' and tri[2] in ('A', 'T')), True
    if ref == 'G' and alt in ('A', 'C'):
        rc = ''.join(_COMP[b] for b in reversed(tri))
        return (rc[0] == 'T' and rc[2] in ('A', 'T')), True
    return False, False


def classify_tcw_ctonly(ref, alt, tri):
    """C>T-only APOBEC TCW classification (SBS2 arm; excludes the C>G/SBS13 arm).
    Returns (is_tcw, is_c_class) where is_c_class flags C>T only
    (incl. G>A on the minus strand).
    """
    if tri is None or len(tri) != 3:
        return False, False
    if ref == 'C' and alt == 'T':
        return (tri[0] == 'T' and tri[2] in ('A', 'T')), True
    if ref == 'G' and alt == 'A':
        rc = ''.join(_COMP[b] for b in reversed(tri))
        return (rc[0] == 'T' and rc[2] in ('A', 'T')), True
    return False, False


def _anxa1_domain(aa_pos):
    for s, e, name in ANXA1_DOMAINS:
        if s <= aa_pos <= e:
            return name
    return 'Unknown'


def beat3_anxa1():
    """ANXA1 specifics: score, IC50, peptides, hotspots, expression, TCW."""
    banner("BEAT 3: ANXA1 DEEP DIVE (Panels C-D)")

    sbs2_neo = pd.read_csv(os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv"), sep="\t")
    cnv_neo = pd.read_csv(os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv"), sep="\t")

    anxa1_sbs2 = sbs2_neo[sbs2_neo['gene'] == 'ANXA1']
    anxa1_cnv = cnv_neo[cnv_neo['gene'] == 'ANXA1']

    log(f"  ANXA1 in SBS2-HIGH neoantigens: {len(anxa1_sbs2)} peptides")
    log(f"  ANXA1 in CNV-HIGH neoantigens: {len(anxa1_cnv)} peptides")
    log(f"  ANXA1 absent from CNV-HIGH: {len(anxa1_cnv) == 0}")

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

        if 'mut_pos_protein' in anxa1_pep.columns:
            binders = (anxa1_pep[anxa1_pep['is_binder'] == True]
                       if 'is_binder' in anxa1_pep.columns
                       else anxa1_pep[anxa1_pep['mut_ic50'] < 500])

            r1 = binders[binders['mut_pos_protein'] == 289]
            r2 = binders[binders['mut_pos_protein'].isin([321, 322, 326, 327])]
            other = binders[~binders['mut_pos_protein'].isin([289, 321, 322, 326, 327])]

            log(f"\n  Hotspot regions (among binders):")
            log(f"    Pos 289: {len(r1)} peptides")
            log(f"    Pos 321-327: {len(r2)} peptides")
            log(f"    Other positions: {len(other)} peptides")
            all_mut_pos = binders['mut_pos_protein'].unique()
            log(f"    Unique mutation positions with binders: {sorted(all_mut_pos)}")

    # Tier assignment
    tier_path = os.path.join(TIER_DIR, "all_genes_tiered.tsv")
    if os.path.exists(tier_path):
        tiers = pd.read_csv(tier_path, sep="\t")
        anxa1_tier = tiers[tiers['gene'] == 'ANXA1'] if 'gene' in tiers.columns else pd.DataFrame()
        if len(anxa1_tier) > 0:
            log(f"\n  ANXA1 tier assignment:")
            for col in anxa1_tier.columns:
                log(f"    {col}: {anxa1_tier.iloc[0][col]}")

    # ANXA1 somatic variants (HGVS list)
    banner("ANXA1 SOMATIC VARIANTS", char="-")
    spa_path = os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv")
    anxa1_var = pd.DataFrame()
    if os.path.exists(spa_path):
        spa = pd.read_csv(spa_path, sep="\t")
        anxa1_var = spa[spa['gene'] == 'ANXA1'].copy()
        log(f"  ANXA1 protein-altering variants: {len(anxa1_var)}")
        positions = []
        for _, row in anxa1_var.iterrows():
            hgvs = str(row.get('hgvs_p', ''))
            m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
            if m:
                pos = int(m.group(2))
                positions.append(pos)
                log(f"    {hgvs} (position {pos}, {_anxa1_domain(pos)})")
        log(f"  Unique positions: {sorted(set(positions))}")
        # Flag positions named vs not named in current text
        named_in_text = {289, 321, 322, 327}
        extra = sorted(set(positions) - named_in_text)
        log(f"  Positions NOT named in current text: {extra}")

    # -------------------------------------------------------------------------
    # ANXA1 TCW / APOBEC context check  [ADD B]
    # Ported from Diagnostic_ANXA1_TCW_Context.py. After the Step01 Start+1 fix
    # the annotation pos equals the SComatic Start directly (Start == pos).
    # -------------------------------------------------------------------------
    banner("ANXA1 TCW / APOBEC CONTEXT CHECK", char="-")
    if len(anxa1_var) == 0:
        log("  [WARNING] No ANXA1 variants loaded; skipping TCW check")
        return

    # Resolve the chrom/pos/ref/alt column names in the somatic file
    cols = {c.lower(): c for c in anxa1_var.columns}
    chrom_c = cols.get('chrom', cols.get('#chrom', 'chrom'))
    pos_c = cols.get('pos', 'pos')
    ref_c = cols.get('ref', 'ref')
    alt_c = cols.get('alt', 'alt')

    # Build the SComatic Start key set from the annotation pos (Start == pos)
    want = {}
    for _, row in anxa1_var.iterrows():
        try:
            chrom = str(row[chrom_c])
            vcf_pos = int(row[pos_c])
        except Exception:
            continue
        scomatic_start = vcf_pos  # annotation pos == SComatic Start (after Step01 Start+1 fix)
        want[(chrom, scomatic_start)] = {
            'vcf_pos': vcf_pos,
            'ref': str(row[ref_c]),
            'alt': str(row[alt_c]),
            'hgvs_p': str(row.get('hgvs_p', '')),
        }

    log(f"  Looking up REF_TRI for {len(want)} ANXA1 variants in SComatic master")
    log(f"  (matching SComatic Start == annotation pos)")

    found = {}
    if os.path.exists(SCOMATIC_TSV):
        with open(SCOMATIC_TSV) as f:
            header = f.readline().rstrip('\n').split('\t')
            ci = header.index('#CHROM')
            si = header.index('Start')
            ri = header.index('REF')
            ti = header.index('REF_TRI')
            for line in f:
                fields = line.rstrip('\n').split('\t')
                key = (fields[ci], int(fields[si]))
                if key in want and key not in found:
                    found[key] = {'ref': fields[ri], 'ref_tri': fields[ti].upper()}
                    if len(found) == len(want):
                        break
        log(f"  Resolved {len(found)}/{len(want)} positions from SComatic")
    else:
        log(f"  [WARNING] SComatic master not found: {SCOMATIC_TSV}")

    # Classify TCW
    results = []
    for key, var in want.items():
        chrom, scstart = key
        ref, alt = var['ref'], var['alt']
        hgvs = var['hgvs_p']
        m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
        aa_pos = int(m.group(2)) if m else 0
        lk = found.get(key)
        if lk is None:
            tri, is_tcw, detail = 'NOT_FOUND', False, 'position not in SComatic'
        else:
            tri = lk['ref_tri']
            is_tcw = False
            detail = ''
            if ref == 'C' and alt in ['T', 'G'] and len(tri) == 3:
                if tri[0] == 'T' and tri[2] in ['A', 'T']:
                    is_tcw = True
                    detail = f"TCW(+): {tri[0]}[C>{alt}]{tri[2]}"
                else:
                    detail = f"non-TCW(+): {tri[0]}[C>{alt}]{tri[2]}"
            elif ref == 'G' and alt in ['A', 'C'] and len(tri) == 3:
                rc = ''.join(_COMP[b] for b in reversed(tri))
                if rc[0] == 'T' and rc[2] in ['A', 'T']:
                    is_tcw = True
                    detail = f"TCW(-): {tri}->rc {rc}"
                else:
                    detail = f"non-TCW(-): {tri}->rc {rc}"
            else:
                detail = f"not C>T/G class: {ref}>{alt}"
        results.append({
            'aa_pos': aa_pos, 'hgvs_p': hgvs, 'vcf_pos': var['vcf_pos'],
            'ref_alt': f"{ref}>{alt}", 'ref_tri': tri, 'is_tcw': is_tcw,
            'label': 'APOBEC' if is_tcw else 'non-APOBEC',
            'domain': _anxa1_domain(aa_pos), 'detail': detail,
        })

    n_apobec = sum(1 for r in results if r['is_tcw'])
    n_non = len(results) - n_apobec
    log(f"\n  Total ANXA1 mutations: {len(results)}")
    log(f"  APOBEC (TCW context):  {n_apobec}")
    log(f"  Non-APOBEC:            {n_non}")
    log(f"\n  {'pos':>5s}  {'hgvs_p':18s}  {'ref>alt':7s}  {'TRI':7s}  {'class':10s}  {'domain':12s}  detail")
    log(f"  {'-'*5}  {'-'*18}  {'-'*7}  {'-'*7}  {'-'*10}  {'-'*12}  {'-'*20}")
    for r in sorted(results, key=lambda x: x['aa_pos']):
        log(f"  {r['aa_pos']:5d}  {r['hgvs_p']:18s}  {r['ref_alt']:7s}  "
            f"{r['ref_tri']:7s}  {r['label']:10s}  {r['domain']:12s}  {r['detail']}")
    log(f"\n  TEXT CLAIM CHECK: 'mutations do not fall in TCW contexts' -> "
        f"{'CONFIRMED (0 APOBEC)' if n_apobec == 0 else 'FALSE: ' + str(n_apobec) + ' are TCW'}")


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
            eff_col = 'effect' if 'effect' in spa.columns else ('annotation' if 'annotation' in spa.columns else None)
            if eff_col:
                log(f"    By {eff_col} type:")
                for effect, count in spa[eff_col].value_counts().items():
                    log(f"      {effect}: {count}")
            if 'gene' in spa.columns:
                log(f"    Unique genes: {spa['gene'].nunique()}")
        else:
            log(f"  [WARNING] {spa_path} not found")

    # Proteome mapping rate from per-group diagnostics
    banner("PROTEOME MAPPING STATS", char="-")
    for group in ['SBS2_HIGH', 'CNV_HIGH']:
        d = os.path.join(MHC_DIR, f"{group}_proteome_mapping_diagnostics.tsv")
        if os.path.exists(d):
            md = pd.read_csv(d, sep="\t")
            log(f"\n  {group} proteome mapping diagnostics: {len(md)} rows")
            # Try to summarise a status column if present
            status_col = None
            for c in ['status', 'result', 'outcome', 'category']:
                if c in md.columns:
                    status_col = c
                    break
            if status_col:
                total = len(md)
                for st, n in md[status_col].value_counts().items():
                    log(f"    {st}: {n} ({100*n/total:.1f}%)")
            else:
                log(f"    columns: {list(md.columns)}")


# =============================================================================
# BEAT 5: TCW-CONTEXT PREVALENCE ACROSS GROUPS  [ADD - population profiling]
# =============================================================================

def _load_scomatic_tri_index(needed_keys):
    """One pass over SComatic master; return {(chrom, 0-based-start): REF_TRI}.
    needed_keys is a set of (chrom, scomatic_start) to keep memory bounded.
    """
    idx = {}
    if not os.path.exists(SCOMATIC_TSV):
        log(f"  [WARNING] SComatic master not found: {SCOMATIC_TSV}")
        return idx
    with open(SCOMATIC_TSV) as f:
        header = f.readline().rstrip('\n').split('\t')
        ci = header.index('#CHROM')
        si = header.index('Start')
        ti = header.index('REF_TRI')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            try:
                key = (fields[ci], int(fields[si]))
            except (ValueError, IndexError):
                continue
            if key in needed_keys and key not in idx:
                idx[key] = fields[ti].upper()
                if len(idx) == len(needed_keys):
                    break
    return idx


def beat5_tcw_prevalence():
    """Compare TCW-context prevalence between SBS2-HIGH and CNV-HIGH across
    two variant sets (protein-altering, neoantigen-producing) and two
    fraction denominators (of C-substitution class, of all variants).

    Runs under TWO substitution scopes so both stay on hand for the manuscript:
      - C>T and C>G  (SBS2 + SBS13; the current manuscript text)
      - C>T only     (SBS2 arm; if the C>G/SBS13 arm is dropped)

    Supports the claim: while ANXA1 itself is non-TCW, SBS2-HIGH may carry a
    higher overall TCW fraction, tying its neoantigen excess to A3 activity.
    """
    banner("BEAT 5: TCW-CONTEXT PREVALENCE (SBS2-HIGH vs CNV-HIGH)")

    try:
        from scipy.stats import fisher_exact
        have_scipy = True
    except Exception:
        have_scipy = False
        log("  [INFO] scipy not available; Fisher tests will be skipped")

    # ---- Step 1: load per-group protein-altering variants (chrom/pos/ref/alt)
    pa = {}
    for group in ['SBS2_HIGH', 'CNV_HIGH']:
        spa_path = os.path.join(ANNOTATION_DIR, f"{group}.somatic_protein_altering.tsv")
        if not os.path.exists(spa_path):
            log(f"  [WARNING] {spa_path} not found")
            continue
        df = pd.read_csv(spa_path, sep="\t")
        cols = {c.lower(): c for c in df.columns}
        chrom_c = cols.get('chrom', cols.get('#chrom'))
        pos_c = cols.get('pos')
        ref_c = cols.get('ref')
        alt_c = cols.get('alt')
        if not all([chrom_c, pos_c, ref_c, alt_c]):
            log(f"  [WARNING] {group}: missing chrom/pos/ref/alt columns "
                f"(have {list(df.columns)})")
            continue
        df = df.rename(columns={chrom_c: 'chrom', pos_c: 'pos',
                                ref_c: 'ref', alt_c: 'alt'})
        # Unique variant per (chrom,pos,ref,alt). Keep hgvs_p + location key.
        df['scstart'] = df['pos'].astype(int)  # annotation pos == SComatic Start (after Step01 Start+1 fix)
        df['location'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str)
        df = df.drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'])
        pa[group] = df
        log(f"  {group}: {len(df)} unique protein-altering variants loaded")

    if not pa:
        log("  [WARNING] No protein-altering variants; aborting BEAT 5")
        return

    # ---- Step 2: identify neoantigen-producing variants (>=1 binder)
    # Map binders back to variants via (location, hgvs_p).
    neo_variant_keys = {}
    for group in ['SBS2_HIGH', 'CNV_HIGH']:
        all_pep_path = os.path.join(MHC_DIR, f"{group}_all_peptide_results.tsv")
        keys = set()
        if os.path.exists(all_pep_path):
            ap = pd.read_csv(all_pep_path, sep="\t")
            if 'is_binder' in ap.columns:
                binders = ap[ap['is_binder'] == True]
            elif 'mut_ic50' in ap.columns:
                binders = ap[ap['mut_ic50'] < 500]
            else:
                binders = ap
            loc_col = 'location' if 'location' in binders.columns else None
            hg_col = 'hgvs_p' if 'hgvs_p' in binders.columns else None
            if loc_col and hg_col:
                keys = set(zip(binders[loc_col].astype(str),
                               binders[hg_col].astype(str)))
        neo_variant_keys[group] = keys
        log(f"  {group}: {len(keys)} unique neoantigen-producing variants (>=1 binder)")

    # ---- Step 3: gather REF_TRI for all needed variants in one SComatic pass
    needed = set()
    for group, df in pa.items():
        for _, r in df.iterrows():
            needed.add((str(r['chrom']), int(r['scstart'])))
    log(f"\n  Resolving REF_TRI for {len(needed)} unique variant positions "
        f"(single SComatic pass)...")
    tri_idx = _load_scomatic_tri_index(needed)
    log(f"  Resolved {len(tri_idx)}/{len(needed)} positions")

    # ---- Step 4: tabulator (substitution scope selectable via classifier)
    def tabulate(df, subset_keys, label, classifier):
        """Return dict of counts for a variant set under a given classifier."""
        n_total = n_cclass = n_tcw = n_resolved = 0
        for _, r in df.iterrows():
            if subset_keys is not None:
                k = (str(r['location']), str(r.get('hgvs_p', '')))
                if k not in subset_keys:
                    continue
            n_total += 1
            tri = tri_idx.get((str(r['chrom']), int(r['scstart'])))
            if tri is None:
                continue
            n_resolved += 1
            is_tcw, is_c = classifier(str(r['ref']), str(r['alt']), tri)
            if is_c:
                n_cclass += 1
            if is_tcw:
                n_tcw += 1
        return {
            'label': label, 'n_total': n_total, 'n_resolved': n_resolved,
            'n_cclass': n_cclass, 'n_tcw': n_tcw,
            'frac_of_c': (n_tcw / n_cclass) if n_cclass else float('nan'),
            'frac_of_all': (n_tcw / n_total) if n_total else float('nan'),
        }

    # ---- Steps 5-6: run, print matrix, and Fisher-test under one definition
    def run_definition(def_label, classifier):
        rows = []
        for group in ['SBS2_HIGH', 'CNV_HIGH']:
            if group not in pa:
                continue
            rows.append((group, 'protein_altering',
                         tabulate(pa[group], None,
                                  f"{group} / protein-altering", classifier)))
            rows.append((group, 'neoantigen_producing',
                         tabulate(pa[group], neo_variant_keys.get(group, set()),
                                  f"{group} / neoantigen-producing", classifier)))

        banner(f"TCW PREVALENCE MATRIX -- {def_label}", char="-")
        log(f"  {'group':10s}  {'set':22s}  {'n':>5s}  {'Cclass':>6s}  "
            f"{'TCW':>4s}  {'TCW/Cclass':>11s}  {'TCW/all':>8s}")
        log(f"  {'-'*10}  {'-'*22}  {'-'*5}  {'-'*6}  {'-'*4}  {'-'*11}  {'-'*8}")
        for group, vset, t in rows:
            log(f"  {group:10s}  {vset:22s}  {t['n_total']:5d}  {t['n_cclass']:6d}  "
                f"{t['n_tcw']:4d}  {t['frac_of_c']:11.3f}  {t['frac_of_all']:8.3f}")

        banner(f"SBS2 vs CNV FISHER TESTS -- {def_label}", char="-")
        by_set = {}
        for group, vset, t in rows:
            by_set.setdefault(vset, {})[group] = t
        for vset, gd in by_set.items():
            if 'SBS2_HIGH' not in gd or 'CNV_HIGH' not in gd:
                continue
            s = gd['SBS2_HIGH']; c = gd['CNV_HIGH']
            a = s['n_tcw']; b = s['n_cclass'] - s['n_tcw']
            cc = c['n_tcw']; d = c['n_cclass'] - c['n_tcw']
            log(f"\n  {vset} (denominator: C-substitution class):")
            log(f"    SBS2: {a} TCW / {s['n_cclass']} C-class ({s['frac_of_c']:.1%})")
            log(f"    CNV:  {cc} TCW / {c['n_cclass']} C-class ({c['frac_of_c']:.1%})")
            if have_scipy and (a + b) > 0 and (cc + d) > 0:
                try:
                    odds, p = fisher_exact([[a, b], [cc, d]])
                    direction = "SBS2 > CNV" if s['frac_of_c'] > c['frac_of_c'] else "CNV >= SBS2"
                    log(f"    Fisher exact: OR={odds:.2f}, p={p:.3g}  ({direction})")
                    log(f"    -> {'higher TCW fraction in SBS2-HIGH, suggesting A3-driven excess' if (s['frac_of_c'] > c['frac_of_c'] and p < 0.05) else 'no significant TCW-fraction difference at p<0.05'}")
                except Exception as e:
                    log(f"    Fisher failed: {e}")
        return rows

    # Both definitions: current text (C>T + C>G) and the SBS2-only (C>T) variant
    run_definition("C>T and C>G  (SBS2 + SBS13; current manuscript text)", classify_tcw)
    run_definition("C>T only  (SBS2 arm)", classify_tcw_ctonly)

    log(f"\n  NOTE: TCW fraction is reported per VARIANT (deduped), not per peptide,")
    log(f"  so a variant producing many binders is counted once.")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("DIAGNOSTIC: SECTION 4.5 (FIGURE 7) NUMBERS")
    log(f"Base directory: {BASE_DIR}")

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
    beat5_tcw_prevalence()

    banner("DIAGNOSTIC COMPLETE")


if __name__ == "__main__":
    main()
