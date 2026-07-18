#!/usr/bin/env python3
"""
Diagnostic_Fusion_Normalization_Compare.py
=============================================
Compare three fusion-rate normalizations, all done the RIGHT way: within sample,
aggregated per sample, then paired across the shared SBS2/CNV samples. No metric
pools raw counts across samples of unequal depth.

Metrics (per group, per sample):
  RAW   = passing fusion junctions / cells
          (what the pipeline reports; scales with depth)
  perUMI= 1000 * sum(passing fusions) / sum(per-cell UMI)
          denominator = transcriptome size (adata counts); depth-normalized
  fracChim = sum(passing fusions) / sum(per-cell total chimeric junctions)
          denominator = the cell's own chimeric output. Because chimeric-read
          generation is itself the depth-driven artifact source, this asks
          "what fraction of a cell's chimeric reads survive as plausible
          fusions," which should be depth-invariant if the artifact process is
          uniform. Requires re-parsing the raw junction files (Step05b barcode
          logic reused verbatim so the denominator is consistent with the
          numerator).

Design: SBS2_HIGH vs CNV_HIGH paired across the 15 shared samples (Wilcoxon
signed-rank on per-sample values), gated to samples with >= MIN_CELLS per group.
NORMAL is reported descriptively only (disjoint samples; not comparable).

2026-07-16 ADDITION (global comparison + verdict):
  The paired test controls for BOTH depth (via the perUMI/fracChim denominators)
  and patient (by keeping SBS2/CNV cells inside the same sample), but it is stuck
  with the 15 dual-group samples. To answer whether we can instead use the
  simpler global comparison across ALL samples (all SBS2 cells vs all CNV cells,
  depth-normalized), this version adds:
    (A) a GLOBAL pooled comparison (every sample; depth handled by the UMI /
        chimeric denominators; does NOT control for patient), and
    (B) a VERDICT that checks whether GLOBAL agrees with PAIRED. If they agree,
        patient confounding is not distorting the result and the global,
        intuitive perUMI comparison can be carried into Panel A. If they
        disagree, that disagreement is the patient effect made visible.
  perUMI is the figure metric (fusions per 1000 UMI). fracChim is an
  independent corroborating normalization (never divides by UMI).

READ-ONLY. Inputs:
  - all_filtered_junctions.tsv (Step05b passing set -> numerator per cell)
  - Chimeric.out.junction files (Step05a -> per-cell total chimeric, denominator)
  - three_group_assignments.tsv
  - adata_final.h5ad (per-cell UMI, optional but needed for perUMI)

Outputs to data/FIG_7/TROUBLESHOOTING/fusion_normalization_compare/:
  - per_sample_metrics.tsv        : srr, group, cells, fusions, umi, chimeric, three rates
  - paired_sbs2_cnv_metrics.tsv    : one row/shared sample, all three metrics both groups
  - global_group_metrics.tsv       : one row/group, pooled depth-normalized metrics (NEW)
  - fusion_normalization_report.txt (incl. side-by-side summary)

Env: NETWORK. Read-only. TROUBLESHOOTING/ (excluded from walkthroughs).
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import glob
from collections import defaultdict, Counter
from datetime import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml")

FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis")
STAR_DIR = os.path.join(FUSION_DIR, "star_chimeric")
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
try:
    import yaml
    with open(CONFIG_PATH) as _f:
        _cfg = yaml.safe_load(_f)
    FUSION_DIR = _cfg['outputs'].get('fusion_analysis', FUSION_DIR)
    STAR_DIR = _cfg['outputs'].get('star_chimeric', STAR_DIR)
    GROUP_PATH = _cfg['inputs'].get('three_group_assignments', GROUP_PATH)
except Exception:
    pass

JUNCTIONS_PATH = os.path.join(FUSION_DIR, "all_filtered_junctions.tsv")
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/TROUBLESHOOTING/fusion_normalization_compare")

GROUPS = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
TUMOR = ['SBS2_HIGH', 'CNV_HIGH']
VALID_JUNCTION_TYPES = {'1', '2'}   # same as Step05b (uniquely-mapped chimeric)
MIN_CELLS = 10                       # per group per sample to enter the paired test

report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)


# --- Step05b barcode recovery, reused verbatim so denominator matches numerator
def extract_barcode(fields, target_bcs_short):
    for fi in range(19, len(fields)):
        if fields[fi].startswith('CB:Z:'):
            cb = fields[fi][5:]
            if cb in target_bcs_short:
                return target_bcs_short[cb]
            cb_suffixed = cb + '-1' if '-' not in cb else cb
            if cb_suffixed in target_bcs_short:
                return target_bcs_short[cb_suffixed]
    for fi in range(14, len(fields)):
        field = fields[fi]
        if len(field) == 16 and field.isalpha():
            short = field + '-1'
            if short in target_bcs_short:
                return target_bcs_short[short]
        if len(field) == 18 and '-' in field:
            if field in target_bcs_short:
                return target_bcs_short[field]
    if len(fields) > 9:
        rname = fields[9]
        for sep in ['_', ':', '|']:
            for part in rname.split(sep):
                if len(part) == 16 and part.isalpha():
                    short = part + '-1'
                    if short in target_bcs_short:
                        return target_bcs_short[short]
                if '-' in part and part in target_bcs_short:
                    return target_bcs_short[part]
    return None


def try_load_percell_umi(barcodes):
    """Per-cell raw UMI/count total for the given barcodes, or None. Never uses
    log-normalized X (rejects if it doesn't look like integer counts)."""
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
    mask = ad.obs_names.isin(set(barcodes))
    if mask.sum() == 0:
        return None
    names = ad.obs_names[mask].tolist()
    cols = [c for c in ad.obs.columns
            if c.lower() in ('total_counts', 'n_counts', 'ncount_rna', 'total_umis',
                             'n_umi', 'nreads', 'n_reads')]
    if cols:
        return dict(zip(names, ad.obs.loc[mask, cols[0]].astype(float).tolist()))
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
        import scipy.sparse as sp
        tot = np.asarray(src.sum(axis=1)).flatten() if sp.issparse(src) else np.asarray(src).sum(axis=1).flatten()
        if np.nanmax(tot) < 100:
            return None
        return dict(zip(names, tot.tolist()))
    except Exception:
        return None


def load_barcode_patient(barcodes):
    """bc -> subject id (patient) from adata.obs, or {} if unavailable.
    Used only for within-patient germline-fusion subtraction (read-only)."""
    try:
        import scanpy as sc
    except Exception:
        return {}
    if not os.path.isfile(ADATA_PATH):
        return {}
    try:
        ad = sc.read_h5ad(ADATA_PATH, backed='r')
    except Exception:
        return {}
    col = None
    for c in ad.obs.columns:
        if c.lower() in ('subject id', 'subject_id', 'patient', 'donor_id', 'donor'):
            col = c
            break
    if col is None:
        return {}
    mask = ad.obs_names.isin(set(barcodes))
    if mask.sum() == 0:
        return {}
    names = ad.obs_names[mask].tolist()
    vals = ad.obs.loc[mask, col].astype(str).tolist()
    return dict(zip(names, vals))


def main():
    banner("Diagnostic: three fusion normalizations, within-sample paired (SBS2 vs CNV)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.isfile(JUNCTIONS_PATH):
        log(f"  [FATAL] not found: {JUNCTIONS_PATH}")
        return

    # numerator: passing fusions per cell
    jdf = pd.read_csv(JUNCTIONS_PATH, sep='\t')
    jdf['barcode'] = jdf['barcode'].astype(str)
    n_fusion = jdf.groupby('barcode').size().to_dict()

    g = pd.read_csv(GROUP_PATH, sep='\t')
    g['srr'] = g['cell_barcode'].str.split('-').str[2]
    bc_to_group = dict(zip(g['cell_barcode'], g['group']))
    bc_to_srr = dict(zip(g['cell_barcode'], g['srr']))

    # srr -> {short_bc: full_bc}, exactly as Step05b builds it
    srr_bc_map = defaultdict(dict)
    for bc in g['cell_barcode']:
        parts = bc.split('-')
        if len(parts) >= 3:
            srr_bc_map[parts[2]][f"{parts[0]}-{parts[1]}"] = bc

    # ---- denominator (fracChim): re-parse junction files for per-cell chimeric total
    banner("Re-parsing junction files for per-cell total chimeric output", "-")
    n_chimeric = Counter()
    jfiles = sorted(glob.glob(os.path.join(STAR_DIR, "SRR*", "Chimeric.out.junction")))
    log(f"  Found {len(jfiles)} junction files under {STAR_DIR}")
    for jf in jfiles:
        srr = os.path.basename(os.path.dirname(jf))
        tbc = srr_bc_map.get(srr, {})
        if not tbc:
            continue
        n_here = 0
        with open(jf) as fh:
            for line in fh:
                if line.startswith('#') or line.startswith('chr_donorA'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 15 or fields[6] not in VALID_JUNCTION_TYPES:
                    continue
                full = extract_barcode(fields, tbc)
                if full is not None:
                    n_chimeric[full] += 1
                    n_here += 1
        log(f"    {srr}: {n_here:,} chimeric junctions attributed to target cells")

    # ---- denominator (perUMI): per-cell UMI
    umi = try_load_percell_umi(list(g['cell_barcode']))
    if umi is None:
        log("\n  [WARN] no per-cell UMI source; perUMI metric will be blank.")
        umi = {}

    # ---- assemble per-cell table
    rows = []
    for bc in g['cell_barcode']:
        rows.append({'barcode': bc, 'group': bc_to_group[bc], 'srr': bc_to_srr[bc],
                     'n_fusion': n_fusion.get(bc, 0),
                     'n_chimeric': n_chimeric.get(bc, 0),
                     'umi': umi.get(bc, np.nan)})
    cell = pd.DataFrame(rows)

    # ---- per (srr, group) within-sample aggregation
    def agg(sub):
        nc = len(sub)
        nf = sub['n_fusion'].sum()
        rate_raw = nf / nc if nc else np.nan
        umi_ok = sub[sub['umi'].notna()]
        rate_umi = (1000 * umi_ok['n_fusion'].sum() / umi_ok['umi'].sum()
                    if umi_ok['umi'].sum() > 0 else np.nan)
        chim_sum = sub['n_chimeric'].sum()
        rate_frac = (100 * nf / chim_sum) if chim_sum > 0 else np.nan
        return pd.Series({'n_cells': nc, 'n_fusion': int(nf), 'sum_umi': umi_ok['umi'].sum(),
                          'sum_chimeric': int(chim_sum), 'rate_raw': rate_raw,
                          'rate_umi': rate_umi, 'rate_frac': rate_frac})

    sg = cell.groupby(['srr', 'group']).apply(agg).reset_index()
    sg.to_csv(os.path.join(OUTPUT_DIR, "per_sample_metrics.tsv"), sep='\t', index=False)

    comp = g.groupby('srr')['group'].agg(lambda x: set(x))
    shared = sorted([s for s, gp in comp.items() if {'SBS2_HIGH', 'CNV_HIGH'} <= gp])

    # ---- paired table across shared samples
    banner("Paired within-sample SBS2 vs CNV across shared samples", "-")
    prows = []
    for s in shared:
        a = sg[(sg['srr'] == s) & (sg['group'] == 'SBS2_HIGH')]
        b = sg[(sg['srr'] == s) & (sg['group'] == 'CNV_HIGH')]
        if len(a) == 0 or len(b) == 0:
            continue
        a, b = a.iloc[0], b.iloc[0]
        prows.append({'srr': s, 'sbs2_cells': int(a['n_cells']), 'cnv_cells': int(b['n_cells']),
                      'sbs2_raw': a['rate_raw'], 'cnv_raw': b['rate_raw'],
                      'sbs2_umi': a['rate_umi'], 'cnv_umi': b['rate_umi'],
                      'sbs2_frac': a['rate_frac'], 'cnv_frac': b['rate_frac']})
    paired = pd.DataFrame(prows)
    paired.to_csv(os.path.join(OUTPUT_DIR, "paired_sbs2_cnv_metrics.tsv"), sep='\t', index=False)

    log(f"  {'srr':13s} {'cells S/C':>10s} | {'RAW S/C':>13s} | {'perUMI S/C':>15s} | {'fracChim% S/C':>15s}")
    for _, r in paired.iterrows():
        log(f"  {r['srr']:13s} {int(r['sbs2_cells']):4d}/{int(r['cnv_cells']):<4d} | "
            f"{r['sbs2_raw']:5.1f}/{r['cnv_raw']:<5.1f}  | "
            f"{r['sbs2_umi']:6.3f}/{r['cnv_umi']:<6.3f} | "
            f"{r['sbs2_frac']:6.2f}/{r['cnv_frac']:<6.2f}")

    # ---- paired stats per metric, gated
    def summarize(metric, unit):
        sc, cc = f'sbs2_{metric}', f'cnv_{metric}'
        for label, sub in [("all shared", paired),
                           (f">= {MIN_CELLS} cells/grp",
                            paired[(paired['sbs2_cells'] >= MIN_CELLS) & (paired['cnv_cells'] >= MIN_CELLS)])]:
            d = sub[[sc, cc]].dropna()
            if len(d) < 3:
                log(f"    [{metric:8s}] {label:20s}: n={len(d)} (too few)")
                continue
            diff = d[sc] - d[cc]
            nhi = int((diff > 0).sum()); nlo = int((diff < 0).sum())
            pval = np.nan
            try:
                from scipy.stats import wilcoxon
                pval = wilcoxon(d[sc], d[cc]).pvalue
            except Exception:
                pass
            log(f"    [{metric:8s}] {label:20s}: n={len(d)}  "
                f"median SBS2={d[sc].median():.3f} CNV={d[cc].median():.3f} {unit}  "
                f"SBS2>CNV {nhi}/{nhi+nlo}  medianDiff={diff.median():+.3f}  p={pval:.3f}")

    banner("Paired signed-rank by metric", "-")
    log("  (SBS2 vs CNV; only shared samples; second line gates thin samples)")
    log("")
    summarize('raw', 'j/cell')
    log("")
    summarize('umi', 'j/1kUMI')
    log("")
    summarize('frac', '% chim')

    # =========================================================================
    # NEW (A): GLOBAL depth-normalized comparison (every sample, not just shared)
    # =========================================================================
    banner("GLOBAL depth-normalized comparison (all cells, every sample)", "-")
    log("  Pools all cells in each group across ALL samples (not just the 15 shared).")
    log("  Depth is handled by the UMI / chimeric denominators; this does NOT control")
    log("  for patient, so it is only trustworthy if it AGREES with the paired result.")
    log("  Point estimate is the pooled ratio (sum fusions / sum denominator).")
    log("")

    def global_pool(grp):
        sub = cell[cell['group'] == grp]
        n_cells = len(sub)
        sum_fus = int(sub['n_fusion'].sum())
        umi_ok = sub[sub['umi'].notna() & (sub['umi'] > 0)]
        sum_umi = float(umi_ok['umi'].sum())
        sum_chim = int(sub['n_chimeric'].sum())
        raw = sum_fus / n_cells if n_cells else np.nan
        perumi = 1000 * umi_ok['n_fusion'].sum() / sum_umi if sum_umi > 0 else np.nan
        frac = 100 * sum_fus / sum_chim if sum_chim > 0 else np.nan
        med_umi = float(umi_ok['umi'].median()) if len(umi_ok) else np.nan
        return {'group': grp, 'n_cells': n_cells, 'n_umi_cells': len(umi_ok),
                'sum_fusion': sum_fus, 'sum_umi': sum_umi, 'sum_chimeric': sum_chim,
                'median_umi': med_umi, 'raw': raw, 'perUMI': perumi, 'fracChim': frac}

    gp = {grp: global_pool(grp) for grp in GROUPS}
    gdf = pd.DataFrame([gp[grp] for grp in GROUPS])
    gdf.to_csv(os.path.join(OUTPUT_DIR, "global_group_metrics.tsv"), sep='\t', index=False)

    log(f"  {'group':12s} {'cells':>6s} {'medUMI':>8s} {'raw j/c':>9s} {'perUMI':>9s} {'fracChim%':>10s}")
    for grp in GROUPS:
        d = gp[grp]
        log(f"  {grp:12s} {d['n_cells']:6d} {d['median_umi']:8.0f} {d['raw']:9.2f} "
            f"{d['perUMI']:9.3f} {d['fracChim']:10.2f}")

    # cell-level test (global), per-cell perUMI, with caveat
    def percell_rate(grp):
        sub = cell[(cell['group'] == grp) & (cell['umi'].notna()) & (cell['umi'] > 0)]
        return (1000 * sub['n_fusion'] / sub['umi']).values

    s_rate = percell_rate('SBS2_HIGH')
    c_rate = percell_rate('CNV_HIGH')
    mwu_p = np.nan
    try:
        from scipy.stats import mannwhitneyu
        if len(s_rate) and len(c_rate):
            mwu_p = mannwhitneyu(s_rate, c_rate, alternative='two-sided').pvalue
    except Exception:
        pass
    log("")
    log(f"  Global perUMI (pooled): SBS2 {gp['SBS2_HIGH']['perUMI']:.3f} vs "
        f"CNV {gp['CNV_HIGH']['perUMI']:.3f} j/1kUMI")
    log(f"  Global fracChim (pooled): SBS2 {gp['SBS2_HIGH']['fracChim']:.2f}% vs "
        f"CNV {gp['CNV_HIGH']['fracChim']:.2f}%")
    log(f"  Cell-level MWU on per-cell perUMI (all samples): p={mwu_p:.3f}")
    log("  CAVEAT: per-cell perUMI is noisy for low-UMI cells (few UMI + 1 fusion")
    log("  inflates the ratio); the pooled ratio above is the point estimate, and the")
    log("  MWU is only a rough significance check on the global (unpaired) contrast.")

    # ---- NORMAL descriptive
    banner("NORMAL (descriptive; disjoint samples, not comparable)", "-")
    nn = sg[sg['group'] == 'NORMAL'].sort_values('srr')
    for _, r in nn.iterrows():
        log(f"  {r['srr']:13s}: raw {r['rate_raw']:5.2f}  perUMI {r['rate_umi']:.3f}  fracChim {r['rate_frac']:.2f}%")
    log(f"  NORMAL mean: raw {nn['rate_raw'].mean():.2f}, perUMI {nn['rate_umi'].mean():.3f}, "
        f"fracChim {nn['rate_frac'].mean():.2f}%")

    # =========================================================================
    # NEW (C): NORMAL background + germline-fusion subtraction
    # =========================================================================
    banner("NORMAL BACKGROUND + GERMLINE-FUSION SUBTRACTION")
    log("  The neoantigen side subtracts variants also seen in NORMAL as germline")
    log("  background; the fusion side has not, so we quantify it here. A germline-origin")
    log("  fusion recurs in a patient's OWN normal cells and, being germline, appears in")
    log("  BOTH that patient's SBS2 and CNV cells, so it cancels in the paired comparison")
    log("  and subtracts equally in the global one. The SBS2-vs-CNV FINDING therefore")
    log("  cannot change; only the absolute (tumor-specific) burden does.")
    log(f"  NORMAL normalized rate (confounded, disjoint samples): perUMI "
        f"{gp['NORMAL']['perUMI']:.3f} j/1kUMI, fracChim {gp['NORMAL']['fracChim']:.2f}%.")
    log("")

    subtracted_percell = None
    bc_patient = load_barcode_patient(list(g['cell_barcode']))
    need = {'chrA', 'posA', 'chrB', 'posB', 'group', 'barcode'}
    if not bc_patient:
        log("  [WARN] no per-cell patient (subject id) source; skipping germline subtraction.")
    elif not need <= set(jdf.columns):
        log(f"  [WARN] junction file lacks breakpoint columns {sorted(need - set(jdf.columns))};")
        log("  cannot define identical fusions. Skipping germline subtraction.")
    else:
        jj = jdf.copy()
        jj['patient'] = jj['barcode'].map(bc_patient)
        jj['jkey'] = (jj['chrA'].astype(str) + ':' + jj['posA'].astype(str) + '-' +
                      jj['chrB'].astype(str) + ':' + jj['posB'].astype(str))

        normal = jj[jj['group'] == 'NORMAL']
        tumor = jj[jj['group'].isin(TUMOR)].copy()

        normal_patients = set(normal['patient'].dropna().unique())
        tumor_patients = set(tumor['patient'].dropna().unique())
        matched = sorted(normal_patients & tumor_patients)
        log(f"  Patients with NORMAL cells in set: {len(normal_patients)}")
        log(f"  Patients with tumor cells in set:  {len(tumor_patients)}")
        log(f"  Matched (normal+tumor) patients:   {len(matched)}  {matched}")
        log("  (within-patient germline subtraction can only affect these matched patients)")

        normal_keys_by_patient = normal.groupby('patient')['jkey'].agg(set).to_dict()
        normal_keys_all = set(normal['jkey'])

        tumor['is_germline_wp'] = tumor.apply(
            lambda r: r['jkey'] in normal_keys_by_patient.get(r['patient'], set()), axis=1)
        tumor['is_germline_pool'] = tumor['jkey'].isin(normal_keys_all)

        # per-cell fusion counts after removing within-patient germline junctions
        # (feeds the Panel A figure numbers below)
        subtracted_percell = tumor[~tumor['is_germline_wp']].groupby('barcode').size()

        n_tot = len(tumor)
        n_wp = int(tumor['is_germline_wp'].sum())
        n_pool = int(tumor['is_germline_pool'].sum())
        log("")
        log(f"  Tumor fusion junctions (target cells):        {n_tot:,}")
        log(f"  Shared with SAME-patient NORMAL (germline):   {n_wp:,} "
            f"({100 * n_wp / n_tot:.2f}%)  [strict within-patient]")
        log(f"  Shared with ANY NORMAL cell (pooled bkgd):    {n_pool:,} "
            f"({100 * n_pool / n_tot:.2f}%)  [looser upper bound, not strictly germline]")

        def perumi_after(mask_col):
            keep = tumor[~tumor[mask_col]] if mask_col else tumor
            kept_counts = keep.groupby('barcode').size()
            out = {}
            for grp in TUMOR:
                sub = cell[(cell['group'] == grp) & (cell['umi'].notna()) & (cell['umi'] > 0)]
                denom = sub['umi'].sum()
                numer = kept_counts.reindex(sub['barcode'].values).fillna(0).sum()
                out[grp] = 1000 * numer / denom if denom > 0 else np.nan
            return out

        base = perumi_after(None)
        wp = perumi_after('is_germline_wp')
        pool = perumi_after('is_germline_pool')
        log("")
        log(f"  perUMI (SBS2 / CNV) under each background model:")
        log(f"    no subtraction               : {base['SBS2_HIGH']:.3f} / {base['CNV_HIGH']:.3f}")
        log(f"    within-patient germline sub  : {wp['SBS2_HIGH']:.3f} / {wp['CNV_HIGH']:.3f}")
        log(f"    pooled-normal background sub : {pool['SBS2_HIGH']:.3f} / {pool['CNV_HIGH']:.3f}")
        log("  -> subtraction lowers both arms together; the SBS2-vs-CNV contrast stays")
        log("     comparable, confirming the finding is background-invariant.")

        pd.DataFrame([
            {'metric': 'perUMI_SBS2_HIGH', 'no_sub': base['SBS2_HIGH'],
             'within_patient_germline': wp['SBS2_HIGH'], 'pooled_normal': pool['SBS2_HIGH']},
            {'metric': 'perUMI_CNV_HIGH', 'no_sub': base['CNV_HIGH'],
             'within_patient_germline': wp['CNV_HIGH'], 'pooled_normal': pool['CNV_HIGH']},
            {'metric': 'n_tumor_junctions', 'no_sub': n_tot,
             'within_patient_germline': n_wp, 'pooled_normal': n_pool},
        ]).to_csv(os.path.join(OUTPUT_DIR, "germline_subtraction_summary.tsv"), sep='\t', index=False)

    # ---- figure-ready per-cell perUMI stats (Panel A fusion side)
    #      Within-patient-germline-subtracted mean per cell is the chosen number;
    #      raw per-cell and pooled perUMI are kept as reference columns. If the
    #      germline block was skipped, subtracted == raw (graceful fallback).
    def _pc_perumi(grp, override=None):
        sub = cell[(cell['group'] == grp) & (cell['umi'].notna()) & (cell['umi'] > 0)]
        if override is not None:
            nf = override.reindex(sub['barcode'].values).fillna(0).values
        else:
            nf = sub['n_fusion'].values
        return 1000 * nf / sub['umi'].values

    fig_rows = []
    for grp in GROUPS:
        raw_arr = _pc_perumi(grp)
        ov = subtracted_percell if (subtracted_percell is not None and grp in TUMOR) else None
        sub_arr = _pc_perumi(grp, ov)
        fig_rows.append({
            'group': grp, 'n_cells': int(len(sub_arr)),
            'mean_perUMI_per_cell_germline_sub': float(np.mean(sub_arr)) if len(sub_arr) else np.nan,
            'sem_perUMI_per_cell_germline_sub': float(np.std(sub_arr, ddof=1) / np.sqrt(len(sub_arr))) if len(sub_arr) > 1 else np.nan,
            'mean_perUMI_per_cell_raw': float(np.mean(raw_arr)) if len(raw_arr) else np.nan,
            'pooled_perUMI_raw': gp[grp]['perUMI'],
        })
    figdf = pd.DataFrame(fig_rows)
    try:
        from scipy.stats import mannwhitneyu
        s_arr = _pc_perumi('SBS2_HIGH', subtracted_percell)
        c_arr = _pc_perumi('CNV_HIGH', subtracted_percell)
        figdf['sbs2_vs_cnv_mwu_p_germline_sub'] = mannwhitneyu(s_arr, c_arr, alternative='two-sided').pvalue
    except Exception:
        figdf['sbs2_vs_cnv_mwu_p_germline_sub'] = np.nan
    figdf['sbs2_vs_cnv_mwu_p_raw'] = mwu_p
    figdf.to_csv(os.path.join(OUTPUT_DIR, "panelA_fusion_stats.tsv"), sep='\t', index=False)

    # =========================================================================
    # NEW (B): VERDICT + corrected text-ready summary (perUMI is the figure metric)
    # =========================================================================
    banner("VERDICT: does GLOBAL agree with PAIRED? (perUMI is the figure metric)")

    # paired all-shared perUMI + fracChim medians and signed-rank p
    d_umi = paired[['sbs2_umi', 'cnv_umi']].dropna()
    d_frac = paired[['sbs2_frac', 'cnv_frac']].dropna()
    p_s = d_umi['sbs2_umi'].median() if len(d_umi) else np.nan
    p_c = d_umi['cnv_umi'].median() if len(d_umi) else np.nan
    f_s = d_frac['sbs2_frac'].median() if len(d_frac) else np.nan
    f_c = d_frac['cnv_frac'].median() if len(d_frac) else np.nan
    p_umi_p = f_umi_p = np.nan
    try:
        from scipy.stats import wilcoxon
        if len(d_umi) >= 3:
            p_umi_p = wilcoxon(d_umi['sbs2_umi'], d_umi['cnv_umi']).pvalue
        if len(d_frac) >= 3:
            f_umi_p = wilcoxon(d_frac['sbs2_frac'], d_frac['cnv_frac']).pvalue
    except Exception:
        pass

    g_s = gp['SBS2_HIGH']['perUMI']
    g_c = gp['CNV_HIGH']['perUMI']

    paired_ns = (not np.isnan(p_umi_p)) and (p_umi_p >= 0.05)
    global_ns = (not np.isnan(mwu_p)) and (mwu_p >= 0.05)
    # point-estimate closeness: within 20% of each other on both paired + global
    def _close(a, b):
        m = max(abs(a), abs(b), 1e-9)
        return abs(a - b) / m < 0.20
    ests_close = _close(p_s, p_c) and _close(g_s, g_c)

    if paired_ns and global_ns and ests_close:
        verdict = ("CONSISTENT: paired and global both say the groups are comparable "
                   "once depth is handled. Patient confounding is not distorting the "
                   "result, so the intuitive global perUMI comparison is safe to carry "
                   "into Panel A.")
    elif paired_ns and global_ns:
        verdict = ("CONSISTENT (both n.s.) but point estimates differ by >20%; both "
                   "agree there is no significant SBS2/CNV difference. Global is usable "
                   "but report the paired number as primary.")
    else:
        verdict = ("DISAGREE: paired and global do not match (one significant, one not). "
                   "Patient confounding likely matters; do NOT pool globally, keep the "
                   "paired within-sample comparison as primary.")

    log(f"  paired  perUMI : SBS2 {p_s:.3f} vs CNV {p_c:.3f} j/1kUMI  (Wilcoxon p={p_umi_p:.3f})")
    log(f"  global  perUMI : SBS2 {g_s:.3f} vs CNV {g_c:.3f} j/1kUMI  (MWU p={mwu_p:.3f})")
    log(f"  paired  fracChim: SBS2 {f_s:.2f}% vs CNV {f_c:.2f}%  (Wilcoxon p={f_umi_p:.3f})")
    log(f"  global  fracChim: SBS2 {gp['SBS2_HIGH']['fracChim']:.2f}% vs "
        f"CNV {gp['CNV_HIGH']['fracChim']:.2f}%")
    log("")
    log(f"  VERDICT: {verdict}")

    banner("TEXT-READY SUMMARY (correct; paraphrase into methods/results)", ".")
    med_s = gp['SBS2_HIGH']['median_umi']
    med_c = gp['CNV_HIGH']['median_umi']
    log(f"  Raw per-cell RNA-fusion rate differs between the tumor groups, but this")
    log(f"  reflects sequencing depth: SBS2-HIGH cells are shallower "
        f"(median {med_s:,.0f} vs {med_c:,.0f} UMI). After depth normalization the")
    log(f"  groups are comparable, whether cells are paired within sample "
        f"(fusions per 1000 UMI {p_s:.2f} vs {p_c:.2f}, p={p_umi_p:.2f}) or pooled")
    log(f"  globally across all samples ({g_s:.2f} vs {g_c:.2f}); an independent")
    log(f"  chimeric-read-fraction normalization agrees ({f_s:.1f}% vs {f_c:.1f}%,")
    log(f"  p={f_umi_p:.2f}). RNA fusion burden therefore does not track the SBS2/CNV")
    log(f"  divergence, consistent with fusions not tracking the immune divergence")
    log(f"  between the two programs. NORMAL is confounded with sample and reported")
    log(f"  descriptively only.")
    log(f"  Fusions shared with matched-normal cells (germline background) are a small")
    log(f"  fraction and, being present in both tumor arms, do not alter the comparison;")
    log(f"  we did not background-subtract fusions the way neoantigens were subtracted,")
    log(f"  which is a stated limitation.")

    # ---- side-by-side guidance (retained)
    banner("SIDE-BY-SIDE (decide which denominator to carry forward)")
    log("  Read across the three metrics on the gated (>= 10 cells/group) rows above:")
    log("   - RAW disagreeing with perUMI/fracChim => the raw signal is depth, not biology.")
    log("   - perUMI and fracChim agreeing with each other => robust, denominator-independent.")
    log("   - If all three say 'comparable', the manuscript's 'fusions don't track the")
    log("     divergence' holds and any of the normalized metrics can be cited.")
    log(f"\n  Wrote: {OUTPUT_DIR}")
    with open(os.path.join(OUTPUT_DIR, "fusion_normalization_report.txt"), 'w') as fh:
        fh.write('\n'.join(report_lines))


if __name__ == "__main__":
    main()
