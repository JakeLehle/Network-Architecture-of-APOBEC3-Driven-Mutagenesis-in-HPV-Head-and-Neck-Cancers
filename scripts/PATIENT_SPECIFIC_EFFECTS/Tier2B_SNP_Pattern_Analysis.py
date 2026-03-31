#!/usr/bin/env python3
"""
Tier2B_SNP_Pattern_Analysis.py (v3 — fixed ALT column: Base_observed)
======================================================================

SComatic SNP pattern analysis across patients.

v3 fixes:
  - Uses 'Base_observed' as the ALT allele (the actual observed variant base
    in each cell), not 'ALT_expected' (comma-separated expected alts across
    cell types) or 'REF' (which was the v2 fallback producing REF>REF artifacts)
  - Explicit column mapping for SComatic output format from ClusterCatcher

SComatic column format:
  #CHROM, Start, End, REF, ALT_expected, Cell_type_expected,
  Num_cells_expected, CB, Cell_type_observed, Base_observed,
  Num_reads, Total_depth, REF_TRI, ALT_TRI

Usage: conda run -n NETWORK python Tier2B_SNP_Pattern_Analysis.py
"""

import os, sys, numpy as np, pandas as pd, scanpy as sc
from itertools import combinations
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from patient_config import *


def try_barcode_mapping(mut_barcodes, adata_barcodes):
    """Try multiple strategies to map SComatic barcodes to adata barcodes."""
    adata_set = set(adata_barcodes)
    mut_unique = set(mut_barcodes)

    # Strategy 1: Direct
    direct = mut_unique & adata_set
    log(f"    Strategy 1 (direct): {len(direct):,} / {len(mut_unique):,}")
    if len(direct) > len(mut_unique) * 0.5:
        return {b: b for b in direct}, "direct"

    # Strategy 2: short -> long (strip SRR suffix from adata)
    short_to_long = {}
    for ab in adata_barcodes:
        parts = ab.rsplit('-', 1)
        if len(parts) == 2 and parts[1].startswith('SRR'):
            short_to_long.setdefault(parts[0], []).append(ab)
        else:
            parts2 = ab.split('-')
            if len(parts2) >= 3 and parts2[-1].startswith('SRR'):
                short_to_long.setdefault('-'.join(parts2[:-1]), []).append(ab)
    matched = {mb: short_to_long[mb][0] for mb in mut_unique if mb in short_to_long}
    log(f"    Strategy 2 (short->long): {len(matched):,} / {len(mut_unique):,}")
    if len(matched) > len(mut_unique) * 0.3:
        return matched, "short_to_long"

    # Strategy 3: prefix swap (SRR_barcode -> barcode-SRR)
    prefix_swap = {}
    for mb in mut_unique:
        if '_' in mb:
            parts = mb.split('_', 1)
            swapped = f"{parts[1]}-{parts[0]}"
            if swapped in adata_set:
                prefix_swap[mb] = swapped
    log(f"    Strategy 3 (prefix swap): {len(prefix_swap):,} / {len(mut_unique):,}")
    if len(prefix_swap) > len(mut_unique) * 0.3:
        return prefix_swap, "prefix_swap"

    log("    WARNING: No strategy achieved >30% match")
    best = max([(direct, "direct"), (matched, "short_to_long"),
                (prefix_swap, "prefix_swap")], key=lambda x: len(x[0]))
    if isinstance(best[0], set):
        return {b: b for b in best[0]}, best[1]
    return best[0], best[1]


def main():
    banner("TIER 2B: SNP PATTERN ANALYSIS (v3 — fixed ALT)")
    out_dir = ensure_dir(DIR_02_SNP)

    # =========================================================================
    # STEP 1: LOAD SCOMATIC DATA WITH EXPLICIT COLUMN MAPPING
    # =========================================================================
    banner("STEP 1: LOAD SCOMATIC MUTATION TABLE")
    if not os.path.exists(SCOMATIC_PATH):
        log(f"  ERROR: {SCOMATIC_PATH}"); sys.exit(1)

    log(f"  Loading: {SCOMATIC_PATH}")
    mut_df = pd.read_csv(SCOMATIC_PATH, sep='\t')
    log(f"  Shape: {mut_df.shape}")
    log(f"  Columns: {list(mut_df.columns)}")

    # Explicit column mapping for ClusterCatcher SComatic output
    # #CHROM  Start  End  REF  ALT_expected  Cell_type_expected
    # Num_cells_expected  CB  Cell_type_observed  Base_observed
    # Num_reads  Total_depth  REF_TRI  ALT_TRI
    CHROM = '#CHROM'
    POS   = 'Start'
    REF   = 'REF'
    ALT   = 'Base_observed'   # <<< THIS IS THE FIX: actual observed variant base
    CELL  = 'CB'

    # Verify columns exist
    for name, col in [('CHROM', CHROM), ('POS', POS), ('REF', REF),
                       ('ALT', ALT), ('CELL', CELL)]:
        if col not in mut_df.columns:
            log(f"  ERROR: Expected column '{col}' not found!")
            log(f"  Available: {list(mut_df.columns)}")
            sys.exit(1)

    log(f"  Column mapping: CHROM={CHROM}, POS={POS}, REF={REF}, ALT={ALT}, CELL={CELL}")
    log(f"  Total mutation calls: {len(mut_df):,}")
    log(f"  Unique cells: {mut_df[CELL].nunique():,}")

    # Verify ALT is actually different from REF (sanity check)
    same_as_ref = (mut_df[REF] == mut_df[ALT]).sum()
    diff_from_ref = (mut_df[REF] != mut_df[ALT]).sum()
    log(f"  Base_observed == REF: {same_as_ref:,} ({100*same_as_ref/len(mut_df):.1f}%)")
    log(f"  Base_observed != REF: {diff_from_ref:,} ({100*diff_from_ref/len(mut_df):.1f}%)")

    # Filter to actual variants (Base_observed != REF)
    if same_as_ref > 0:
        log(f"  Filtering to actual variants (Base_observed != REF)...")
        mut_df = mut_df[mut_df[REF] != mut_df[ALT]].copy()
        log(f"  After filtering: {len(mut_df):,} mutation calls")

    # Create variant ID
    mut_df['variant_id'] = (mut_df[CHROM].astype(str) + ':' +
                            mut_df[POS].astype(str) + ':' +
                            mut_df[REF].astype(str) + '>' +
                            mut_df[ALT].astype(str))
    log(f"  Unique variants: {mut_df['variant_id'].nunique():,}")
    log(f"  Example variants: {mut_df['variant_id'].unique()[:5].tolist()}")

    # =========================================================================
    # STEP 2: MAP CELL BARCODES TO PATIENTS
    # =========================================================================
    banner("STEP 2: MAP CELLS TO PATIENTS")
    adata = load_adata()
    high_cells, low_cells = load_groups()

    adata_barcodes = list(adata.obs_names)
    mut_barcodes = mut_df[CELL].unique().tolist()
    log(f"  SComatic barcodes: {len(mut_barcodes):,}")
    log(f"  adata barcodes: {len(adata_barcodes):,}")
    log(f"  SComatic example: {mut_barcodes[0]}")
    log(f"  adata example:    {adata_barcodes[0]}")

    bc_mapping, strategy = try_barcode_mapping(mut_barcodes, adata_barcodes)
    log(f"  Strategy: {strategy} ({len(bc_mapping):,} mapped)")

    if len(bc_mapping) == 0:
        log("  ERROR: No barcodes mapped."); sys.exit(1)

    bc_to_patient  = adata.obs[PATIENT_COL].to_dict()
    bc_to_tissue   = adata.obs[TISSUE_COL].to_dict()
    bc_to_celltype = adata.obs[CELLTYPE_COL].to_dict()

    mut_df['adata_bc'] = mut_df[CELL].map(bc_mapping)
    mut_df['patient']  = mut_df['adata_bc'].map(bc_to_patient)
    mut_df['tissue']   = mut_df['adata_bc'].map(bc_to_tissue)
    mut_df['celltype'] = mut_df['adata_bc'].map(bc_to_celltype)
    mut_df['is_HIGH']  = mut_df['adata_bc'].isin(high_cells)
    mut_df['is_basal'] = mut_df['celltype'] == 'basal cell'

    mapped = mut_df['patient'].notna()
    log(f"  Mapped: {mapped.sum():,} / {len(mut_df):,} ({100*mapped.mean():.1f}%)")
    mut_df = mut_df[mapped].copy()

    log(f"\n  Mutations per patient:")
    for p in sorted(mut_df['patient'].unique()):
        pm = mut_df[mut_df['patient'] == p]
        flag = " <<<" if p in HIGH_CONTRIBUTORS else ""
        log(f"    {p}: {pm['variant_id'].nunique():,} variants, "
            f"{pm[CELL].nunique():,} cells, "
            f"{pm[pm['is_HIGH']]['variant_id'].nunique():,} in HIGH{flag}")

    # =========================================================================
    # STEP 3: PATIENT-LEVEL VARIANT COMPARISON
    # =========================================================================
    banner("STEP 3: PATIENT-LEVEL VARIANT COMPARISON")
    patients = sorted(mut_df['patient'].unique())
    patient_variants = {p: set(mut_df[mut_df['patient']==p]['variant_id'].unique())
                        for p in patients}
    for p in patients:
        log(f"    {p}: {len(patient_variants[p]):,}")

    # Pairwise Jaccard
    n_p = len(patients)
    jaccard_matrix = np.zeros((n_p, n_p))
    shared_matrix = np.zeros((n_p, n_p), dtype=int)
    for i, p1 in enumerate(patients):
        for j, p2 in enumerate(patients):
            s1, s2 = patient_variants[p1], patient_variants[p2]
            inter = len(s1 & s2)
            union = len(s1 | s2)
            jaccard_matrix[i,j] = inter / union if union > 0 else 0
            shared_matrix[i,j] = inter

    # Shared vs patient-specific
    variant_patients = {}
    for p, vs in patient_variants.items():
        for v in vs:
            variant_patients.setdefault(v, set()).add(p)

    shared_variants = {v: ps for v, ps in variant_patients.items() if len(ps) >= 2}
    patient_specific = {v: list(ps)[0] for v, ps in variant_patients.items() if len(ps) == 1}
    log(f"\n  Shared (>=2 patients): {len(shared_variants):,}")
    log(f"  Patient-specific: {len(patient_specific):,}")

    ps_counts = {}
    for v, p in patient_specific.items():
        ps_counts[p] = ps_counts.get(p, 0) + 1
    for p in patients:
        log(f"    {p}: {ps_counts.get(p,0):,} patient-specific")

    # =========================================================================
    # STEP 4: HIGH-CONTRIBUTOR SHARED VARIANTS
    # =========================================================================
    banner("STEP 4: HIGH-CONTRIBUTOR SHARED VARIANTS")
    hc_in = [p for p in HIGH_CONTRIBUTORS if p in patient_variants]
    log(f"  High-contributors with mutations: {hc_in}")

    if len(hc_in) >= 2:
        hc_vars = {p: patient_variants[p] for p in hc_in}
        hc_shared_all = set.intersection(*hc_vars.values())
        log(f"  Shared by ALL high-contributors: {len(hc_shared_all):,}")
        for p1, p2 in combinations(hc_in, 2):
            shared = hc_vars[p1] & hc_vars[p2]
            log(f"    {p1.replace('Patient ','')} & {p2.replace('Patient ','')}: {len(shared):,}")

        non_hc_all = set()
        for p in patients:
            if p not in set(HIGH_CONTRIBUTORS):
                non_hc_all |= patient_variants.get(p, set())
        hc_exclusive = hc_shared_all - non_hc_all
        log(f"  Exclusive to high-contributors: {len(hc_exclusive):,}")
        if len(hc_exclusive) > 0:
            log(f"  Top exclusive variants:")
            for v in sorted(hc_exclusive)[:20]:
                log(f"    {v}")

    # =========================================================================
    # STEP 5: VARIANTS IN SBS2-HIGH CELLS
    # =========================================================================
    banner("STEP 5: VARIANTS IN SBS2-HIGH CELLS")
    high_muts = mut_df[mut_df['is_HIGH']]
    log(f"  Mutations in HIGH cells: {len(high_muts):,}")
    log(f"  Unique variants: {high_muts['variant_id'].nunique():,}")
    high_patient_variants = {}
    for p in patients:
        hp = high_muts[high_muts['patient'] == p]
        hvs = set(hp['variant_id'].unique())
        if len(hvs) > 0:
            high_patient_variants[p] = hvs
            log(f"    {p}: {len(hvs):,}")

    # =========================================================================
    # STEP 6: GERMLINE CANDIDATES
    # =========================================================================
    banner("STEP 6: GERMLINE CANDIDATES")
    basal_muts = mut_df[mut_df['is_basal']]
    log(f"  Basal cell mutations: {len(basal_muts):,}")

    germline_candidates = []
    for p in patients:
        p_basal = basal_muts[basal_muts['patient'] == p]
        if len(p_basal) == 0: continue
        n_bc = p_basal[CELL].nunique()
        if n_bc < 10: continue
        var_cc = p_basal.groupby('variant_id')[CELL].nunique()
        threshold = n_bc * 0.5
        germ = var_cc[var_cc >= threshold]
        if len(germ) > 0:
            log(f"  {p}: {len(germ)} candidates (>{threshold:.0f}/{n_bc} basal cells)")
            for v in germ.index:
                chrom = v.split(':')[0]
                germline_candidates.append({
                    'patient': p, 'variant_id': v, 'chrom': chrom,
                    'n_basal_with_var': int(germ[v]), 'total_basal': n_bc,
                    'fraction': germ[v] / n_bc,
                    'in_HIGH': v in high_patient_variants.get(p, set()),
                    'on_A3_locus': chrom == A3_LOCUS_CHR,
                })

    germ_df = pd.DataFrame(germline_candidates)
    if len(germ_df) > 0:
        germ_path = os.path.join(out_dir, "Tier2B_germline_candidates.tsv")
        germ_df.to_csv(germ_path, sep='\t', index=False)
        log(f"  Saved: {germ_path} ({len(germ_df)} candidates)")
        a3_germ = germ_df[germ_df['on_A3_locus']]
        log(f"  On A3 locus (chr22): {len(a3_germ)}")
        for _, r in a3_germ.iterrows():
            log(f"    {r['patient']}: {r['variant_id']} ({r['n_basal_with_var']}/{r['total_basal']})")
    else:
        log("  No germline candidates found")

    # =========================================================================
    # STEP 7: SAVE AND VISUALIZE
    # =========================================================================
    banner("STEP 7: OUTPUTS")
    summary_rows = []
    for p in patients:
        n_total = len(patient_variants.get(p, set()))
        n_spec = ps_counts.get(p, 0)
        n_hv = len(high_patient_variants.get(p, set()))
        n_germ = len(germ_df[germ_df['patient']==p]) if len(germ_df)>0 else 0
        summary_rows.append({
            'patient': p, 'total_variants': n_total,
            'patient_specific': n_spec, 'shared_with_others': n_total - n_spec,
            'variants_in_HIGH_cells': n_hv, 'germline_candidates': n_germ,
            'high_contributor': p in HIGH_CONTRIBUTORS,
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(out_dir, "Tier2B_patient_variant_summary.tsv"),
                      sep='\t', index=False)
    log(f"  Saved summary")

    if n_p >= 2:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        ax = axes[0]
        labs = [p.replace('Patient ','') for p in patients]
        im = ax.imshow(jaccard_matrix, cmap='YlOrRd', vmin=0,
                       vmax=max(jaccard_matrix.max(), 0.01))
        ax.set_xticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=45, ha='right', fontsize=10)
        ax.set_yticks(range(len(labs)))
        ax.set_yticklabels(labs, fontsize=10)
        for i, p in enumerate(patients):
            if p in HIGH_CONTRIBUTORS:
                ax.get_xticklabels()[i].set_color('#ed6a5a')
                ax.get_xticklabels()[i].set_fontweight('bold')
                ax.get_yticklabels()[i].set_color('#ed6a5a')
                ax.get_yticklabels()[i].set_fontweight('bold')
        plt.colorbar(im, ax=ax, label='Jaccard Similarity', shrink=0.8)
        ax.set_title('Patient Variant Similarity (Jaccard)', fontsize=13, fontweight='bold')

        ax = axes[1]
        ss = summary_df.sort_values('total_variants', ascending=True)
        y = range(len(ss))
        ax.barh(list(y), ss['shared_with_others'].values, color='#5e81ac',
                edgecolor='black', linewidth=0.3, label='Shared')
        ax.barh(list(y), ss['patient_specific'].values,
                left=ss['shared_with_others'].values,
                color='#ed6a5a', edgecolor='black', linewidth=0.3, label='Patient-specific')
        ax.set_yticks(list(y))
        ax.set_yticklabels([p.replace('Patient ','') for p in ss['patient'].values], fontsize=10)
        ax.set_xlabel('Number of Variants', fontsize=12)
        ax.set_title('Variant Composition per Patient', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)

        plt.suptitle('Tier 2B: Patient SNP Pattern Analysis (corrected ALT)',
                     fontsize=15, fontweight='bold', y=1.02)
        plt.tight_layout()
        for ext in ['pdf','png']:
            plt.savefig(os.path.join(out_dir, f"Tier2B_patient_variant_comparison.{ext}"),
                       dpi=300, bbox_inches='tight')
        plt.close(); log("  Saved plots")

    log("\nTier 2B complete.")

if __name__ == "__main__":
    main()
