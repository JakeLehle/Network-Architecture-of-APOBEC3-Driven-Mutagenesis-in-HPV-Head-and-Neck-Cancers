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

READ-ONLY. Inputs:
  - all_filtered_junctions.tsv (Step05b passing set -> numerator per cell)
  - Chimeric.out.junction files (Step05a -> per-cell total chimeric, denominator)
  - three_group_assignments.tsv
  - adata_final.h5ad (per-cell UMI, optional but needed for perUMI)

Outputs to data/FIG_7/TROUBLESHOOTING/fusion_normalization_compare/:
  - per_sample_metrics.tsv        : srr, group, cells, fusions, umi, chimeric, three rates
  - paired_sbs2_cnv_metrics.tsv    : one row/shared sample, all three metrics both groups
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

    # ---- NORMAL descriptive
    banner("NORMAL (descriptive; disjoint samples, not comparable)", "-")
    nn = sg[sg['group'] == 'NORMAL'].sort_values('srr')
    for _, r in nn.iterrows():
        log(f"  {r['srr']:13s}: raw {r['rate_raw']:5.2f}  perUMI {r['rate_umi']:.3f}  fracChim {r['rate_frac']:.2f}%")
    log(f"  NORMAL mean: raw {nn['rate_raw'].mean():.2f}, perUMI {nn['rate_umi'].mean():.3f}, "
        f"fracChim {nn['rate_frac'].mean():.2f}%")

    # ---- side-by-side verdict
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
