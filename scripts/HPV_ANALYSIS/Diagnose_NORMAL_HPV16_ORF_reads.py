#!/usr/bin/env python3
"""
Diagnose_NORMAL_HPV16_ORF_reads.py
==================================
Question: the 8 NORMAL cells with raw_HPV16 >= 8 UMI survive the Panel F
filter (raw_HPV16 >= 8 AND total_hpv16_genome_reads > 0), yet their per-ORF
gene counts (E1..L2) look empty. Where are those reads going?

This is a read-the-tables-only check. It answers:
  1. For each of the 8 NORMAL HPV16+ cells, print every column of the
     gene-count table, plus raw_HPV16 (master), total genome reads,
     the sum across the 8 lifecycle ORFs, and the off-ORF remainder
     (total - ORF_sum).
  2. How many NORMAL HPV16+ cells have ZERO reads in all 8 ORFs.
  3. The same total / ORF_sum / off-ORF breakdown summarized per group
     (SBS2-HIGH, CNV-HIGH, NORMAL), so we can see whether off-ORF reads
     are a NORMAL-only artifact or a universal feature of the alignment.
  4. Whether raw_HPV16 (Cell Ranger feature) and total_hpv16_genome_reads
     (Phase4 genome alignment) agree, since they are separate quantifications.

It does NOT modify any production file or output. Text report only.

Env: NETWORK
Usage: conda run -n NETWORK python Diagnose_NORMAL_HPV16_ORF_reads.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import numpy as np
import pandas as pd
from collections import OrderedDict

# =============================================================================
# CONFIGURATION (paths match the Figure 6 production script)
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

THREE_GROUP_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_4/01_group_selection/three_group_assignments.tsv")
MASTER_TABLE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/01_raw_hpv16_counts/basal_cell_master_table_with_raw_HPV16.tsv")
HPV_GENE_PATH = os.path.join(PROJECT_ROOT,
    "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")

HPV16_THRESHOLD = 8
TOTAL_COL = 'total_hpv16_genome_reads'

POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}

# The 8 lifecycle ORFs scored in Panel F
ORF_GENES = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'L1', 'L2']

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title):
    log("")
    log("=" * 80)
    log(f"  {title}")
    log("=" * 80)

# =============================================================================
# LOAD
# =============================================================================
banner("LOAD")

groups = pd.read_csv(THREE_GROUP_PATH, sep='\t')
cell_to_group = dict(zip(groups['cell_barcode'], groups['group']))
log(f"  Group assignments: {len(groups)} cells")

master = pd.read_csv(MASTER_TABLE_PATH, sep='\t', index_col=0)
master['group'] = master.index.map(lambda x: cell_to_group.get(x, 'other'))
log(f"  Master table: {master.shape}")

hpv = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
log(f"  Gene-count table: {hpv.shape}")
log(f"  Gene-count columns ({len(hpv.columns)}): {list(hpv.columns)}")

# Confirm which ORF columns are actually present
present_orfs = [g for g in ORF_GENES if g in hpv.columns]
missing_orfs = [g for g in ORF_GENES if g not in hpv.columns]
log(f"  ORF columns present: {present_orfs}")
if missing_orfs:
    log(f"  WARNING: ORF columns MISSING from table: {missing_orfs}")
if TOTAL_COL not in hpv.columns:
    log(f"  WARNING: '{TOTAL_COL}' not in gene-count table")

# Any non-ORF, non-total numeric columns? (LCR / intergenic / splice etc.)
known = set(present_orfs) | {TOTAL_COL}
other_cols = [c for c in hpv.columns if c not in known]
log(f"  Other columns in gene-count table (not ORF, not total): {other_cols}")

# Attach group + raw_HPV16 onto the gene-count table by barcode
hpv['group'] = hpv.index.map(lambda x: cell_to_group.get(x, 'other'))
hpv['raw_HPV16'] = hpv.index.map(master['raw_HPV16'])

# Reconstruct the exact Panel F set: raw_HPV16 >= 8 AND total > 0, in the 3 groups
panelF = hpv[hpv['group'].isin(POP_ORDER) &
             (hpv['raw_HPV16'] >= HPV16_THRESHOLD) &
             (hpv[TOTAL_COL] > 0)].copy()
panelF['ORF_sum'] = panelF[present_orfs].sum(axis=1)
panelF['off_ORF'] = panelF[TOTAL_COL] - panelF['ORF_sum']
log(f"  Panel F set reconstructed: {len(panelF)} cells")
for pop in POP_ORDER:
    log(f"    {POP_LABELS[pop]}: {(panelF['group'] == pop).sum()}")

# =============================================================================
# 1. PER-CELL DUMP OF THE 8 NORMAL HPV16+ CELLS
# =============================================================================
banner("1. NORMAL HPV16+ cells (raw_HPV16 >= 8): full per-cell breakdown")

normalF = panelF[panelF['group'] == 'NORMAL'].copy()
log(f"  n = {len(normalF)} NORMAL cells in the Panel F set")
log("")

dump_cols = ['raw_HPV16', TOTAL_COL] + present_orfs + ['ORF_sum', 'off_ORF']
if other_cols:
    dump_cols = ['raw_HPV16', TOTAL_COL] + present_orfs + other_cols + ['ORF_sum', 'off_ORF']

for bc, row in normalF.iterrows():
    log(f"  cell: {bc}")
    log(f"    raw_HPV16 (master feature) : {row['raw_HPV16']:.0f}")
    log(f"    {TOTAL_COL:<26s} : {row[TOTAL_COL]:.0f}")
    orf_str = ", ".join(f"{g}={row[g]:.0f}" for g in present_orfs)
    log(f"    ORF counts                 : {orf_str}")
    if other_cols:
        oth_str = ", ".join(f"{c}={row[c]}" for c in other_cols)
        log(f"    other cols                 : {oth_str}")
    log(f"    ORF_sum                    : {row['ORF_sum']:.0f}")
    log(f"    off-ORF (total - ORF_sum)  : {row['off_ORF']:.0f}")
    log("")

n_all_zero_orf = (normalF['ORF_sum'] == 0).sum()
log(f"  NORMAL cells with ZERO reads across all 8 ORFs: "
    f"{n_all_zero_orf} / {len(normalF)}")
log(f"  NORMAL total genome reads: min={normalF[TOTAL_COL].min():.0f}, "
    f"median={normalF[TOTAL_COL].median():.1f}, max={normalF[TOTAL_COL].max():.0f}")
if (normalF['off_ORF'] < 0).any():
    log("  NOTE: some off-ORF values are NEGATIVE, which means ORF_sum exceeds "
        "total_hpv16_genome_reads. That points to the two counts using different "
        "denominators (e.g. multimapped reads counted per-ORF but not in total).")

# =============================================================================
# 2. PER-GROUP SUMMARY: total vs ORF_sum vs off-ORF
# =============================================================================
banner("2. Per-group breakdown (is off-ORF a NORMAL-only thing?)")

log(f"  {'Group':<12s}  {'n':>5s}  {'mean total':>11s}  {'mean ORF':>9s}  "
    f"{'mean offORF':>12s}  {'%reads offORF':>13s}  {'%cells ORF=0':>13s}")
log(f"  {'-'*12}  {'-'*5}  {'-'*11}  {'-'*9}  {'-'*12}  {'-'*13}  {'-'*13}")
for pop in POP_ORDER:
    sub = panelF[panelF['group'] == pop]
    if len(sub) == 0:
        continue
    mean_total = sub[TOTAL_COL].mean()
    mean_orf = sub['ORF_sum'].mean()
    mean_off = sub['off_ORF'].mean()
    pct_off = 100.0 * sub['off_ORF'].sum() / sub[TOTAL_COL].sum() \
        if sub[TOTAL_COL].sum() > 0 else 0.0
    pct_zero = 100.0 * (sub['ORF_sum'] == 0).sum() / len(sub)
    log(f"  {POP_LABELS[pop]:<12s}  {len(sub):>5d}  {mean_total:>11.2f}  "
        f"{mean_orf:>9.2f}  {mean_off:>12.2f}  {pct_off:>12.1f}%  {pct_zero:>12.1f}%")

log(f"  {'Group':<12s}  {'n':>5s}  {'%reads URR':>11s}  {'%reads interg':>14s}  {'%reads ORF':>11s}")
log(f"  {'-'*12}  {'-'*5}  {'-'*11}  {'-'*14}  {'-'*11}")
for pop in POP_ORDER:
    sub = panelF[panelF['group'] == pop]
    if len(sub) == 0:
        continue
    tot = sub[TOTAL_COL].sum()
    pct_urr = 100.0 * sub['URR'].sum() / tot if tot > 0 else 0.0
    pct_int = 100.0 * sub['intergenic'].sum() / tot if tot > 0 else 0.0
    pct_orf = 100.0 * sub['ORF_sum'].sum() / tot if tot > 0 else 0.0
    log(f"  {POP_LABELS[pop]:<12s}  {len(sub):>5d}  "
        f"{pct_urr:>10.1f}%  {pct_int:>13.1f}%  {pct_orf:>10.1f}%")

# =============================================================================
# 3. PER-ORF detectability across groups (how often is each ORF nonzero?)
# =============================================================================
banner("3. Per-ORF detection rate (% of HPV16+ cells with that ORF > 0)")

log(f"  {'ORF':<6s}  " + "  ".join(f"{POP_LABELS[p]:>10s}" for p in POP_ORDER))
log(f"  {'-'*6}  " + "  ".join("-" * 10 for _ in POP_ORDER))
for g in present_orfs:
    cells = []
    for pop in POP_ORDER:
        sub = panelF[panelF['group'] == pop]
        rate = 100.0 * (sub[g] > 0).sum() / len(sub) if len(sub) > 0 else 0.0
        cells.append(f"{rate:>9.1f}%")
    log(f"  {g:<6s}  " + "  ".join(cells))

# =============================================================================
# 4. raw_HPV16 (feature) vs total_hpv16_genome_reads (alignment) agreement
# =============================================================================
banner("4. raw_HPV16 vs total_hpv16_genome_reads (two separate quantifications)")

log("  These are counted by different steps and need not match. Checking how")
log("  closely they track on the Panel F set:")
log("")
for pop in POP_ORDER:
    sub = panelF[panelF['group'] == pop]
    if len(sub) < 2:
        # still show means for tiny groups, skip correlation
        if len(sub) >= 1:
            log(f"  {POP_LABELS[pop]:<12s}: n={len(sub)}, "
                f"mean raw_HPV16={sub['raw_HPV16'].mean():.1f}, "
                f"mean total_genome={sub[TOTAL_COL].mean():.1f} "
                f"(too few cells for correlation)")
        continue
    r = np.corrcoef(sub['raw_HPV16'].values.astype(float),
                    sub[TOTAL_COL].values.astype(float))[0, 1]
    log(f"  {POP_LABELS[pop]:<12s}: n={len(sub)}, "
        f"mean raw_HPV16={sub['raw_HPV16'].mean():.1f}, "
        f"mean total_genome={sub[TOTAL_COL].mean():.1f}, "
        f"Pearson r={r:.3f}")

# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUTPUT_DIR, "diagnose_normal_hpv16_orf_reads_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")
banner("DIAGNOSTIC COMPLETE")
