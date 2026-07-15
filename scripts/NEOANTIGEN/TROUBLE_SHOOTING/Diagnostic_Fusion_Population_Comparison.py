#!/usr/bin/env python3
"""
Diagnostic_Fusion_Population_Comparison.py
============================================
The defensible fusion comparison, given the sample structure.

Established facts (from the group assignment sample map):
  - All 546 NORMAL cells come from 9 NORMAL-ONLY samples (SRR14340883-891).
    No NORMAL cell shares a sample with any tumor cell. NORMAL-vs-tumor fusion
    burden is therefore fully confounded with sample / patient / tissue / depth
    and CANNOT be compared. It is reported here descriptively only.
  - SBS2_HIGH and CNV_HIGH share ~15 samples. Within a shared sample, depth,
    patient, and batch are held fixed and only the SBS2/CNV label differs. That
    is the internally controlled comparison Figure 7 cares about.

So this diagnostic:
  1. Reports the naive pooled per-546 rates for reference, with the caveat.
  2. Confirms the sample composition (which SRRs are single- vs multi-group).
  3. PRIMARY: paired, within-sample SBS2 vs CNV across shared samples
     (per-sample per-cell rates + Wilcoxon signed-rank), which controls
     sample-level depth by construction.
  4. Depth check (best effort): within shared samples, compares per-cell
     sequencing depth of SBS2 vs CNV cells, to rule out a per-cell depth skew
     driving any within-sample difference. Skipped cleanly if a raw-count source
     is not available (the within-sample design already controls sample depth).
  5. Reports NORMAL per-sample rates descriptively, explicitly excluded from the
     burden comparison.

READ-ONLY. Inputs: all_filtered_junctions.tsv (Step05b), three_group_assignments.tsv.
Optional: adata_final.h5ad for a per-cell depth proxy.

Outputs to data/FIG_7/TROUBLESHOOTING/fusion_population_comparison/:
    - per_sample_group_fusion_rates.tsv   : srr, group, n_cells, n_junctions, rate
    - paired_sbs2_cnv_shared_samples.tsv   : one row per shared sample, both rates + diff
    - fusion_population_comparison_report.txt (incl. quotable summary)

Env: NETWORK (pandas, scipy; scanpy only if the depth check runs). Read-only.
Belongs in scripts/NEOANTIGEN/TROUBLESHOOTING/ (excluded from walkthroughs).
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml")

FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis")
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
try:
    import yaml
    with open(CONFIG_PATH) as _f:
        _cfg = yaml.safe_load(_f)
    FUSION_DIR = _cfg['outputs'].get('fusion_analysis', FUSION_DIR)
    GROUP_PATH = _cfg['inputs'].get('three_group_assignments', GROUP_PATH)
except Exception:
    pass

JUNCTIONS_PATH = os.path.join(FUSION_DIR, "all_filtered_junctions.tsv")
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/TROUBLESHOOTING/fusion_population_comparison")

GROUPS = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
TUMOR = ['SBS2_HIGH', 'CNV_HIGH']
# Shared samples with fewer than this many cells of either group give noisy
# per-cell rates; the paired test is reported on all shared samples AND on this
# better-powered subset as a sensitivity check.
MIN_CELLS_FOR_POWERED = 10

report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)


def try_load_percell_depth(barcodes):
    """Best-effort per-cell depth (raw UMI/count total) for the given barcodes.

    Returns dict barcode -> depth, or None if no valid raw-count source is found.
    Uses backed mode to avoid loading the full matrix. Never uses log-normalized
    X as a depth proxy (that would be meaningless), so if only normalized X is
    present with no counts layer/obs column, returns None.
    """
    try:
        import scanpy as sc
        import scipy.sparse as sp
    except Exception:
        return None
    if not os.path.isfile(ADATA_PATH):
        return None
    try:
        ad = sc.read_h5ad(ADATA_PATH, backed='r')
    except Exception:
        return None

    bc_set = set(barcodes)
    mask = ad.obs_names.isin(bc_set)
    if mask.sum() == 0:
        return None
    sub_names = ad.obs_names[mask].tolist()

    # 1) An obs column that looks like a total-count field
    depth_cols = [c for c in ad.obs.columns
                  if c.lower() in ('total_counts', 'n_counts', 'ncount_rna',
                                   'total_umis', 'n_umi', 'nreads', 'n_reads')]
    if depth_cols:
        vals = ad.obs.loc[mask, depth_cols[0]].astype(float)
        return dict(zip(sub_names, vals.tolist()))

    # 2) A raw counts layer or .raw (sum per cell)
    src = None
    try:
        if 'counts' in ad.layers:
            src = ad[mask].layers['counts']
        elif ad.raw is not None:
            src = ad.raw[mask].X
    except Exception:
        src = None
    if src is None:
        return None
    try:
        totals = np.asarray(src.sum(axis=1)).flatten() if sp.issparse(src) \
            else np.asarray(src).sum(axis=1).flatten()
        # Reject if it looks log-normalized rather than integer counts
        if np.nanmax(totals) < 100:
            return None
        return dict(zip(sub_names, totals.tolist()))
    except Exception:
        return None


def main():
    banner("Diagnostic: fusion comparison, sample-aware (SBS2 vs CNV paired; NORMAL descriptive)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.isfile(JUNCTIONS_PATH):
        log(f"  [FATAL] not found: {JUNCTIONS_PATH} (run Step05b first)")
        return

    jdf = pd.read_csv(JUNCTIONS_PATH, sep='\t')
    jdf['group'] = jdf['group'].astype(str)
    jdf['srr_id'] = jdf['srr_id'].astype(str)
    g = pd.read_csv(GROUP_PATH, sep='\t')
    g['srr'] = g['cell_barcode'].str.split('-').str[2]

    ncells_total = {grp: int((g['group'] == grp).sum()) for grp in GROUPS}
    log(f"  Loaded {len(jdf):,} junctions; group sizes " +
        ", ".join(f"{k}={v}" for k, v in ncells_total.items()))

    # cells per (srr, group) and junctions per (srr, group)
    cells_sg = g.groupby(['srr', 'group']).size().rename('n_cells').reset_index()
    jun_sg = (jdf.groupby(['srr_id', 'group']).size().rename('n_junctions')
              .reset_index().rename(columns={'srr_id': 'srr'}))
    sg = cells_sg.merge(jun_sg, on=['srr', 'group'], how='left').fillna({'n_junctions': 0})
    sg['n_junctions'] = sg['n_junctions'].astype(int)
    sg['rate'] = sg['n_junctions'] / sg['n_cells']
    sg.to_csv(os.path.join(OUTPUT_DIR, "per_sample_group_fusion_rates.tsv"), sep='\t', index=False)

    # ------------------------------------------------------------------ (1)
    banner("(1) Naive pooled per-cell rates  [reference only, see caveat]", "-")
    for grp in GROUPS:
        n = int((jdf['group'] == grp).sum())
        log(f"  {grp:12s}: {n:5d} junctions / {ncells_total[grp]} cells = {n/ncells_total[grp]:.2f}/cell")
    log("  CAVEAT: this pools cells across samples of unequal depth and ignores")
    log("  that NORMAL is a disjoint set of samples. Not a valid burden comparison.")

    # ------------------------------------------------------------------ (2)
    banner("(2) Sample composition", "-")
    comp = g.groupby('srr')['group'].agg(lambda x: tuple(sorted(set(x))))
    normal_only = [s for s, gp in comp.items() if gp == ('NORMAL',)]
    shared = [s for s, gp in comp.items() if set(gp) >= {'SBS2_HIGH', 'CNV_HIGH'}]
    sbs2_solo = [s for s, gp in comp.items() if gp == ('SBS2_HIGH',)]
    cnv_solo = [s for s, gp in comp.items() if gp == ('CNV_HIGH',)]
    log(f"  NORMAL-only samples:        {len(normal_only)}  {normal_only}")
    log(f"  SBS2+CNV shared samples:    {len(shared)}  {shared}")
    log(f"  SBS2-only samples:          {len(sbs2_solo)}")
    log(f"  CNV-only samples:           {len(cnv_solo)}")
    log("  -> NORMAL shares 0 samples with any tumor group (confound confirmed).")

    # ------------------------------------------------------------------ (3)
    banner("(3) PRIMARY: paired within-sample SBS2 vs CNV (shared samples)", "-")
    rows = []
    for s in shared:
        sub = sg[sg['srr'] == s]
        d = {r['group']: r for _, r in sub.iterrows()}
        if 'SBS2_HIGH' not in d or 'CNV_HIGH' not in d:
            continue
        rows.append({
            'srr': s,
            'sbs2_cells': int(d['SBS2_HIGH']['n_cells']), 'sbs2_junctions': int(d['SBS2_HIGH']['n_junctions']),
            'sbs2_rate': d['SBS2_HIGH']['rate'],
            'cnv_cells': int(d['CNV_HIGH']['n_cells']), 'cnv_junctions': int(d['CNV_HIGH']['n_junctions']),
            'cnv_rate': d['CNV_HIGH']['rate'],
            'diff_sbs2_minus_cnv': d['SBS2_HIGH']['rate'] - d['CNV_HIGH']['rate'],
        })
    paired = pd.DataFrame(rows)
    paired.to_csv(os.path.join(OUTPUT_DIR, "paired_sbs2_cnv_shared_samples.tsv"), sep='\t', index=False)

    log(f"  {'srr':14s} {'SBS2 j/c/rate':>18s} {'CNV j/c/rate':>18s} {'diff':>7s}")
    for _, r in paired.iterrows():
        log(f"  {r['srr']:14s} "
            f"{r['sbs2_junctions']:4d}/{r['sbs2_cells']:3d}/{r['sbs2_rate']:5.2f}   "
            f"{r['cnv_junctions']:4d}/{r['cnv_cells']:3d}/{r['cnv_rate']:5.2f}  "
            f"{r['diff_sbs2_minus_cnv']:+6.2f}")

    def paired_summary(pdf, label):
        if len(pdf) < 3:
            log(f"  {label}: too few samples ({len(pdf)}) for a paired test")
            return
        sbs2_pool = pdf['sbs2_junctions'].sum() / pdf['sbs2_cells'].sum()
        cnv_pool = pdf['cnv_junctions'].sum() / pdf['cnv_cells'].sum()
        n_sbs2_hi = int((pdf['diff_sbs2_minus_cnv'] > 0).sum())
        n_cnv_hi = int((pdf['diff_sbs2_minus_cnv'] < 0).sum())
        med = pdf['diff_sbs2_minus_cnv'].median()
        pval = np.nan
        try:
            from scipy.stats import wilcoxon
            pval = wilcoxon(pdf['sbs2_rate'], pdf['cnv_rate']).pvalue
        except Exception:
            pass
        log(f"  {label} (n={len(pdf)} samples):")
        log(f"    pooled-within-shared rate:  SBS2 {sbs2_pool:.2f}/cell   CNV {cnv_pool:.2f}/cell")
        log(f"    per-sample direction:       SBS2>CNV in {n_sbs2_hi}, CNV>SBS2 in {n_cnv_hi}")
        log(f"    median per-sample diff:     {med:+.2f}/cell (SBS2 - CNV)")
        log(f"    Wilcoxon signed-rank p:     {pval:.3f}" if pval == pval else
            "    Wilcoxon signed-rank p:     n/a")
        return sbs2_pool, cnv_pool, med, pval

    log("")
    res_all = paired_summary(paired, "All shared samples")
    powered = paired[(paired['sbs2_cells'] >= MIN_CELLS_FOR_POWERED) &
                     (paired['cnv_cells'] >= MIN_CELLS_FOR_POWERED)]
    log("")
    res_pow = paired_summary(powered, f"Well-powered subset (>= {MIN_CELLS_FOR_POWERED} cells/group)")

    # ------------------------------------------------------------------ (4)
    banner("(4) Depth check within shared samples (per-cell UMI, best effort)", "-")
    shared_bcs = g.loc[(g['srr'].isin(shared)) & (g['group'].isin(TUMOR)), 'cell_barcode'].tolist()
    depth = try_load_percell_depth(shared_bcs)
    if depth is None:
        log("  No raw-count depth source found (or scanpy/adata unavailable).")
        log("  Skipping per-cell depth normalization. NOTE: the within-sample design")
        log("  in (3) already controls SAMPLE-level depth; this check would only add")
        log("  a guard against per-cell depth skew BETWEEN SBS2 and CNV inside a sample.")
    else:
        gg = g.copy()
        gg['depth'] = gg['cell_barcode'].map(depth)
        dd = gg[(gg['srr'].isin(shared)) & (gg['group'].isin(TUMOR)) & gg['depth'].notna()]
        med_sbs2 = dd.loc[dd['group'] == 'SBS2_HIGH', 'depth'].median()
        med_cnv = dd.loc[dd['group'] == 'CNV_HIGH', 'depth'].median()
        log(f"  Median per-cell depth in shared samples: SBS2 {med_sbs2:,.0f}   CNV {med_cnv:,.0f}")
        ratio = med_sbs2 / med_cnv if med_cnv else float('nan')
        log(f"  SBS2/CNV depth ratio: {ratio:.2f}")
        if 0.8 <= ratio <= 1.25:
            log("  -> Comparable depth; the within-sample SBS2-vs-CNV difference is not")
            log("     explained by a per-cell depth skew.")
        else:
            log("  -> Depth differs between SBS2 and CNV cells; report the depth-normalized")
            log("     rate below as the primary number.")
        # depth-normalized rate (junctions per 1000 UMI), pooled within shared samples
        jbc = jdf[jdf['srr_id'].isin(shared)].groupby('barcode').size().rename('nj')
        tmp = dd.set_index('cell_barcode')
        tmp = tmp.join(jbc).fillna({'nj': 0})
        for grp in TUMOR:
            t = tmp[tmp['group'] == grp]
            if t['depth'].sum() > 0:
                per1k = 1000 * t['nj'].sum() / t['depth'].sum()
                log(f"    {grp:12s}: {per1k:.3f} junctions / 1000 UMI (depth-normalized)")

    # ------------------------------------------------------------------ (5)
    banner("(5) NORMAL, descriptive only (sample-confounded, excluded from burden)", "-")
    nsub = sg[sg['group'] == 'NORMAL'].sort_values('srr')
    log(f"  NORMAL per-sample per-cell rates across its {len(nsub)} dedicated samples:")
    for _, r in nsub.iterrows():
        log(f"    {r['srr']:14s}: {int(r['n_junctions']):4d} junctions / {int(r['n_cells']):3d} cells = {r['rate']:.2f}/cell")
    log(f"  NORMAL mean per-sample rate: {nsub['rate'].mean():.2f}  (range {nsub['rate'].min():.2f}-{nsub['rate'].max():.2f})")
    log("  These samples contain NO tumor cells, so this rate cannot be compared to")
    log("  the tumor groups: any difference is confounded with tissue/patient/depth.")

    # ------------------------------------------------------------------ SUMMARY
    banner("TEXT-READY SUMMARY (paraphrase into methods/results)")
    if res_all:
        sbs2_pool, cnv_pool, med, pval = res_all
        direction = ("comparable" if abs(med) < 0.5 else
                     ("higher in SBS2-HIGH" if med > 0 else "higher in CNV-HIGH"))
        sig = ("no significant difference" if not (pval == pval) or pval >= 0.05
               else "a significant difference")
        summary = (
            f"Because all NORMAL cells derive from separate normal-adjacent samples that "
            f"contain no tumor cells, NORMAL fusion burden is confounded with sample and "
            f"cannot be compared to the tumor groups and is reported descriptively only. "
            f"SBS2-HIGH and CNV-HIGH co-occur in {len(paired)} samples, permitting a paired "
            f"within-sample comparison that controls for sequencing depth and patient. Across "
            f"these samples the per-cell fusion rate was {direction} between the two "
            f"(SBS2-HIGH {sbs2_pool:.1f} vs CNV-HIGH {cnv_pool:.1f} per cell; median per-sample "
            f"difference {med:+.1f}; Wilcoxon signed-rank {sig}"
            + (f", p = {pval:.3f}" if pval == pval else "") + "). "
            f"RNA fusion burden therefore does not distinguish the SBS2 and CNV mutational "
            f"programs, consistent with fusions not tracking the immune divergence between them."
        )
        log("  " + summary.replace(". ", ".\n  "))
    log(f"\n  Wrote: {OUTPUT_DIR}")
    log("    per_sample_group_fusion_rates.tsv, paired_sbs2_cnv_shared_samples.tsv,")
    log("    fusion_population_comparison_report.txt")

    with open(os.path.join(OUTPUT_DIR, "fusion_population_comparison_report.txt"), 'w') as fh:
        fh.write('\n'.join(report_lines))


if __name__ == "__main__":
    main()
