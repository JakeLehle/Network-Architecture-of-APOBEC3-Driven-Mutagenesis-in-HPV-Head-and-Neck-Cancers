#!/usr/bin/env python3
"""
diagnostic_section6_e6e7_fdr.py
=================================

Quick diagnostic to verify two things:
  1. Pull the figure-wide BH-adjusted q-values from the figure report
  2. Independently recompute the figure-wide BH correction for Panel F
     lifecycle genes to confirm E6/E7 q-values are correct

Also reports the BH q-values for lifecycle PHASE-level comparisons
(maintenance, amplification, oncogene, capsid) so the text can use
FDR-corrected values matching the figure.

Usage:
    conda run -n NETWORK python diagnostic_section6_e6e7_fdr.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG6_PANELS = os.path.join(BASE_DIR, "data/FIG_6/FIGURE_6_PANELS")
REPORT_PATH = os.path.join(FIG6_PANELS, "figure6_lifecycle_report.txt")

HPV_GENE_PATH = os.path.join(BASE_DIR, "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")
GROUP_PATH = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/three_group_assignments.tsv")

HPV16_PHASES = OrderedDict([
    ('Maintenance',   ['E1', 'E2']),
    ('Amplification', ['E4', 'E5']),
    ('Oncogene',      ['E6', 'E7']),
    ('Capsid',        ['L1', 'L2']),
])

POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# PART 1: Check figure report for E6/E7 lines
# =============================================================================

def check_figure_report():
    """Parse the figure report for F_E6 and F_E7 BH q-values."""
    banner("PART 1: Figure report q-values for E6/E7")

    if not os.path.exists(REPORT_PATH):
        log(f"  [WARNING] Report not found: {REPORT_PATH}")
        return

    log(f"  Reading: {REPORT_PATH}")

    # Find the raw vs adjusted table
    with open(REPORT_PATH) as f:
        lines = f.readlines()

    # Look for lines containing F_E6 or F_E7 or the header
    relevant = []
    in_table = False
    for line in lines:
        stripped = line.strip()
        if 'Panel' in stripped and 'Pair' in stripped and 'Raw p' in stripped:
            in_table = True
            relevant.append(stripped)
            continue
        if in_table and '---' in stripped:
            relevant.append(stripped)
            continue
        if in_table and stripped:
            # Check if it's a table row
            if any(label in stripped for label in ['F_E6', 'F_E7', 'F_E1', 'F_E2',
                                                     'F_E4', 'F_E5', 'F_L1', 'F_L2',
                                                     'B_SBS2', 'B_CNV', 'B_CytoTRACE',
                                                     'B_A3A', 'B_A3B', 'D_HPV16']):
                relevant.append(stripped)

    if relevant:
        log(f"\n  All Panel F + B + D entries from figure report:")
        for line in relevant:
            log(f"    {line}")
    else:
        log("  Could not find the BH correction table in the report.")

    # Specifically highlight E6/E7
    log(f"\n  E6 and E7 specific lines:")
    for line in lines:
        if 'F_E6' in line or 'F_E7' in line:
            log(f"    {line.strip()}")


# =============================================================================
# PART 2: Independent recomputation of figure-wide BH
# =============================================================================

def recompute_figure_wide_bh():
    """Recompute BH correction matching the figure script's approach."""
    banner("PART 2: Independent recomputation of figure-wide BH for Panel F")

    # Load data
    groups_df = pd.read_csv(GROUP_PATH, sep="\t")
    groups_map = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))

    hpv_genes = pd.read_csv(HPV_GENE_PATH, sep="\t", index_col=0)
    hpv_genes['group'] = hpv_genes.index.map(lambda x: groups_map.get(x, 'other'))

    total_col = 'total_hpv16_genome_reads'
    hpv_pop = hpv_genes[hpv_genes['group'].isin(POP_ORDER)].copy()
    hpv_pop = hpv_pop[hpv_pop[total_col] > 0].copy()

    log(f"HPV+ cells in three populations: {len(hpv_pop)}")
    for pop in POP_ORDER:
        log(f"  {pop}: {(hpv_pop['group'] == pop).sum()}")

    # Compute per-gene fractions
    all_genes = []
    for phase_genes in HPV16_PHASES.values():
        all_genes.extend(phase_genes)

    for gene in all_genes:
        if gene in hpv_pop.columns:
            hpv_pop[f'{gene}_frac'] = hpv_pop[gene] / hpv_pop[total_col]

    # Compute ALL pairwise raw p-values (3 pairs per gene, 8 genes = 24 tests)
    pairs = [(0, 1), (1, 2), (0, 2)]  # SBS2vCNV, CNVvNORM, SBS2vNORM
    pair_names = ['SBS2 vs CNV', 'CNV vs NORM', 'SBS2 vs NORM']

    all_raw = []
    all_labels = []

    for gene in all_genes:
        frac_col = f'{gene}_frac'
        if frac_col not in hpv_pop.columns:
            continue
        for pi, (i, j) in enumerate(pairs):
            pop_i = POP_ORDER[i]
            pop_j = POP_ORDER[j]
            vals_i = hpv_pop.loc[hpv_pop['group'] == pop_i, frac_col].values
            vals_j = hpv_pop.loc[hpv_pop['group'] == pop_j, frac_col].values

            if len(vals_i) >= 10 and len(vals_j) >= 10:
                _, p = mannwhitneyu(vals_i, vals_j, alternative='two-sided')
                all_raw.append(p)
            else:
                all_raw.append(np.nan)
            all_labels.append(f"F_{gene} | {pair_names[pi]}")

    # BH correction (Panel F only, 24 tests)
    raw_arr = np.array(all_raw, dtype=float)
    valid = ~np.isnan(raw_arr)
    adj_arr = np.full_like(raw_arr, np.nan)
    if valid.sum() > 0:
        _, adj_p, _, _ = multipletests(raw_arr[valid], method='fdr_bh')
        adj_arr[valid] = adj_p

    log(f"\n  Panel F BH correction ({valid.sum()} valid tests out of {len(raw_arr)}):")
    log(f"\n  {'Label':<30} {'Raw p':>14} {'BH q (F only)':>14} {'Stars':>8}")
    log(f"  {'-'*70}")

    for k in range(len(all_labels)):
        raw_p = all_raw[k]
        adj_q = adj_arr[k]
        if np.isnan(raw_p):
            log(f"  {all_labels[k]:<30} {'N.D.':>14} {'N.D.':>14} {'N.D.':>8}")
        else:
            stars = '****' if adj_q < 1e-4 else '***' if adj_q < 1e-3 else '**' if adj_q < 0.01 else '*' if adj_q < 0.05 else 'ns'
            log(f"  {all_labels[k]:<30} {raw_p:>14.4e} {adj_q:>14.4e} {stars:>8}")

    # Now compute PHASE-level BH (4 phases x 3 pairs = 12 tests)
    banner("PHASE-LEVEL BH CORRECTION (for text)", char="-")

    phase_raw = []
    phase_labels = []

    for phase_name, phase_genes in HPV16_PHASES.items():
        frac_cols = [f'{g}_frac' for g in phase_genes if f'{g}_frac' in hpv_pop.columns]
        if not frac_cols:
            continue
        for pi, (i, j) in enumerate(pairs):
            pop_i = POP_ORDER[i]
            pop_j = POP_ORDER[j]
            vals_i = hpv_pop.loc[hpv_pop['group'] == pop_i, frac_cols].sum(axis=1).values
            vals_j = hpv_pop.loc[hpv_pop['group'] == pop_j, frac_cols].sum(axis=1).values

            if len(vals_i) >= 10 and len(vals_j) >= 10:
                _, p = mannwhitneyu(vals_i, vals_j, alternative='two-sided')
                phase_raw.append(p)
            else:
                phase_raw.append(np.nan)
            phase_labels.append(f"{phase_name} | {pair_names[pi]}")

    phase_arr = np.array(phase_raw, dtype=float)
    valid_ph = ~np.isnan(phase_arr)
    adj_phase = np.full_like(phase_arr, np.nan)
    if valid_ph.sum() > 0:
        _, adj_ph_p, _, _ = multipletests(phase_arr[valid_ph], method='fdr_bh')
        adj_phase[valid_ph] = adj_ph_p

    log(f"\n  {'Phase comparison':<35} {'Raw p':>14} {'BH q':>14}")
    log(f"  {'-'*65}")
    for k in range(len(phase_labels)):
        raw_p = phase_raw[k]
        adj_q = adj_phase[k]
        if np.isnan(raw_p):
            log(f"  {phase_labels[k]:<35} {'N.D.':>14} {'N.D.':>14}")
        else:
            log(f"  {phase_labels[k]:<35} {raw_p:>14.4e} {adj_q:>14.4e}")

    # E6/E7 mean comparison with effect size
    banner("E6/E7 BIOLOGICAL EFFECT SIZE", char="-")
    for gene in ['E6', 'E7']:
        frac_col = f'{gene}_frac'
        sbs2_vals = hpv_pop.loc[hpv_pop['group'] == 'SBS2_HIGH', frac_col].values
        cnv_vals = hpv_pop.loc[hpv_pop['group'] == 'CNV_HIGH', frac_col].values

        sbs2_mean = np.mean(sbs2_vals) * 100
        cnv_mean = np.mean(cnv_vals) * 100
        diff = abs(sbs2_mean - cnv_mean)

        sbs2_med = np.median(sbs2_vals) * 100
        cnv_med = np.median(cnv_vals) * 100

        log(f"\n  {gene}:")
        log(f"    SBS2-HIGH: mean={sbs2_mean:.3f}%, median={sbs2_med:.3f}%, "
            f"n={len(sbs2_vals)}")
        log(f"    CNV-HIGH:  mean={cnv_mean:.3f}%, median={cnv_med:.3f}%, "
            f"n={len(cnv_vals)}")
        log(f"    Absolute difference: {diff:.3f} percentage points")
        log(f"    SBS2 zeros: {(sbs2_vals == 0).sum()}/{len(sbs2_vals)} "
            f"({100*(sbs2_vals==0).sum()/len(sbs2_vals):.1f}%)")
        log(f"    CNV zeros:  {(cnv_vals == 0).sum()}/{len(cnv_vals)} "
            f"({100*(cnv_vals==0).sum()/len(cnv_vals):.1f}%)")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("DIAGNOSTIC: E6/E7 FDR VERIFICATION")

    check_figure_report()
    recompute_figure_wide_bh()

    banner("DIAGNOSTIC COMPLETE")


if __name__ == "__main__":
    main()
