#!/usr/bin/env python3
"""
Diagnostic_Fusion_Artifact_Audit.py
======================================
Figure 7 fusion: is the elevated NORMAL fusion count real, or artifact?

The per-cell fusion rate is only modestly higher in NORMAL (12.1 vs SBS2 9.4,
CNV 10.2), but recurrent and group-exclusive pairs blow out in NORMAL (exclusive
pairs SBS2 100 / CNV 129 / NORMAL 405). Inspection of the top NORMAL-exclusive
pairs points at two artifact-prone locus classes: hyperpolymorphic HLA loci
(HLA-DRB1, HLA-B, plus CD74) and paralogous/abundant ribosomal loci
(RPSA--RPSA2, RPL*). STAR chimeric detection over short reads is known to
produce recurrent false junctions at exactly these hard-to-align regions, the
same category the pipeline already blacklists (NRXN3, RASGRF2).

This diagnostic quantifies that. It is READ-ONLY: it reads the filtered junction
set produced by Step05b and characterizes it. It does not modify the pipeline
output. If the maintainer decides a filter is warranted, that decision and its
exact effect are recorded here so it can be described in the manuscript.

Two artifact signatures are measured:
  1. Promiscuous loci  -- a gene fusing to MANY distinct partners behaves like a
     mapping-noise hotspot (real recurrent fusions have one specific partner).
  2. Hyper-recurrent pairs -- a single pair recurring in an outlier number of
     cells, especially over HLA / ribosomal / paralog loci.

And one confound is controlled:
  3. Within-SRR comparison -- holding the sample (hence depth and batch) fixed,
     is NORMAL still higher than the tumor groups? If the gap survives within a
     sample, it is not a depth/composition artifact of which SRRs each group
     draws from.

Outputs (breadcrumbs for the text) to
  data/FIG_7/TROUBLESHOOTING/fusion_artifact_audit/:
    - dropped_loci_audit.tsv           : every flagged locus, class, partners, per-group junctions
    - per_group_before_after.tsv       : per-group totals/rates before vs after removal
    - promiscuity_sensitivity.tsv       : effect of the threshold choice (3/5/10)
    - hyper_recurrent_pairs.tsv         : single pairs recurring in >= outlier cells
    - within_srr_rates.tsv              : per-SRR per-group rate (depth/batch control)
    - fusion_artifact_audit_report.txt  : full log incl. a quotable summary paragraph

Env: NETWORK or NEOANTIGEN (pandas only). Read-only.
Belongs in scripts/NEOANTIGEN/TROUBLESHOOTING/ (excluded from walkthroughs).
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from collections import defaultdict
from datetime import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml")

# Resolve fusion dir + group file from config where possible, fall back to canonical paths.
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
STAR_DIR = os.path.join(FUSION_DIR, "star_chimeric")  # for optional Log.final.out depth read
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/TROUBLESHOOTING/fusion_artifact_audit")

GROUPS = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']

# Data-driven promiscuity: a locus is artifact-suspect if it fuses to >= this many
# DISTINCT partner genes. Reported at several thresholds; DEFAULT is the one used
# for the "corrected" numbers. Chosen from the distribution printed in STEP 2.
PROMISCUITY_THRESHOLDS = [3, 5, 10]
DEFAULT_THRESHOLD = 5

# A single pair seen in >= this many cells is called hyper-recurrent (outlier).
HYPER_RECURRENT_CELLS = 10

report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)


def locus_class(gene):
    """Descriptive class for a gene, for the audit breadcrumb (not used for filtering)."""
    g = str(gene).upper()
    if g.startswith("HLA-"):
        return "HLA"
    if g.startswith(("RPL", "RPS")) and any(c.isdigit() for c in g):
        return "ribosomal"
    if g.startswith("MRPL") or g.startswith("MRPS"):
        return "mito_ribosomal"
    if g.startswith("MT-"):
        return "mitochondrial"
    return "other"


def per_group_rates(jdf, ncells):
    out = {}
    for grp in GROUPS:
        n = int((jdf['group'] == grp).sum())
        out[grp] = (n, n / ncells[grp] if ncells[grp] else 0.0)
    return out


def main():
    banner("Diagnostic: fusion artifact audit (is the NORMAL excess real?)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    log(f"  Junctions: {JUNCTIONS_PATH}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.isfile(JUNCTIONS_PATH):
        log(f"  [FATAL] Not found: {JUNCTIONS_PATH}  (run Step05b first)")
        return

    jdf = pd.read_csv(JUNCTIONS_PATH, sep='\t')
    for c in ('geneA', 'geneB', 'group', 'barcode', 'srr_id'):
        jdf[c] = jdf[c].astype(str)
    log(f"  Loaded {len(jdf):,} filtered junctions")

    groups_df = pd.read_csv(GROUP_PATH, sep='\t')
    ncells = {g: int((groups_df['group'] == g).sum()) for g in GROUPS}
    log(f"  Group sizes: " + ", ".join(f"{g}={ncells[g]}" for g in GROUPS))

    # ------------------------------------------------------------------ STEP 1
    banner("STEP 1: Baseline per-group fusion counts", "-")
    base = per_group_rates(jdf, ncells)
    for g in GROUPS:
        log(f"  {g:12s}: {base[g][0]:5d} junctions  ({base[g][1]:.2f}/cell)")
    tumor_base = max(base['SBS2_HIGH'][1], base['CNV_HIGH'][1])
    log(f"  NORMAL excess over highest tumor group: {base['NORMAL'][1] - tumor_base:+.2f}/cell")

    # ------------------------------------------------------------------ STEP 2
    banner("STEP 2: Per-gene promiscuity (distinct fusion partners)", "-")
    partners = defaultdict(set)
    gene_junctions = defaultdict(int)
    gene_cells = defaultdict(set)
    for _, r in jdf.iterrows():
        a, b = r['geneA'], r['geneB']
        partners[a].add(b); partners[b].add(a)
        gene_junctions[a] += 1; gene_junctions[b] += 1
        gene_cells[a].add(r['barcode']); gene_cells[b].add(r['barcode'])

    prom = pd.DataFrame([
        {'gene': g, 'class': locus_class(g), 'n_distinct_partners': len(p),
         'n_junctions': gene_junctions[g], 'n_cells': len(gene_cells[g])}
        for g, p in partners.items()
    ]).sort_values('n_distinct_partners', ascending=False).reset_index(drop=True)

    # Distribution so the threshold choice is transparent (data-driven, not arbitrary)
    log("  Distribution of distinct-partner counts per gene:")
    for lo, hi in [(1, 1), (2, 2), (3, 4), (5, 9), (10, 10**9)]:
        n = ((prom['n_distinct_partners'] >= lo) & (prom['n_distinct_partners'] <= hi)).sum()
        label = f">={lo}" if hi >= 10**8 else (f"{lo}" if lo == hi else f"{lo}-{hi}")
        log(f"    partners {label:>5s}: {n} genes")
    log("\n  Top 15 most promiscuous loci:")
    log(f"    {'gene':14s} {'class':14s} {'partners':>8s} {'junctions':>9s} {'cells':>6s}")
    for _, r in prom.head(15).iterrows():
        log(f"    {r['gene']:14s} {r['class']:14s} {r['n_distinct_partners']:8d} "
            f"{r['n_junctions']:9d} {r['n_cells']:6d}")

    # ------------------------------------------------------------------ STEP 3
    banner("STEP 3: Threshold sensitivity (junctions removed per group)", "-")
    log(f"  {'thresh':>6s} {'n_loci':>6s} {'SBS2':>16s} {'CNV':>16s} {'NORMAL':>16s}")
    sens_rows = []
    for thr in PROMISCUITY_THRESHOLDS:
        flagged = set(prom.loc[prom['n_distinct_partners'] >= thr, 'gene'])
        keep = ~(jdf['geneA'].isin(flagged) | jdf['geneB'].isin(flagged))
        after = per_group_rates(jdf[keep], ncells)
        log(f"  {thr:6d} {len(flagged):6d} "
            + " ".join(f"{base[g][0]:5d}->{after[g][0]:5d}({after[g][1]:4.1f})" for g in GROUPS))
        sens_rows.append({'threshold': thr, 'n_loci_removed': len(flagged),
                          **{f'{g}_before': base[g][0] for g in GROUPS},
                          **{f'{g}_after': after[g][0] for g in GROUPS},
                          **{f'{g}_rate_after': round(after[g][1], 3) for g in GROUPS}})
    pd.DataFrame(sens_rows).to_csv(os.path.join(OUTPUT_DIR, "promiscuity_sensitivity.tsv"),
                                   sep='\t', index=False)

    # ------------------------------------------------------------------ STEP 4
    banner("STEP 4: Hyper-recurrent single pairs (outlier recurrence)", "-")
    jdf['_pair'] = jdf.apply(lambda r: tuple(sorted([r['geneA'], r['geneB']])), axis=1)
    pair_rows = []
    for grp in GROUPS:
        gj = jdf[jdf['group'] == grp]
        cells_per_pair = gj.groupby('_pair')['barcode'].nunique().sort_values(ascending=False)
        for pair, ncell in cells_per_pair[cells_per_pair >= HYPER_RECURRENT_CELLS].items():
            cls = "/".join(sorted({locus_class(pair[0]), locus_class(pair[1])} - {"other"})) or "other"
            pair_rows.append({'group': grp, 'geneA': pair[0], 'geneB': pair[1],
                              'n_cells': int(ncell), 'class': cls})
    pair_df = pd.DataFrame(pair_rows).sort_values(['group', 'n_cells'], ascending=[True, False]) \
        if pair_rows else pd.DataFrame(columns=['group', 'geneA', 'geneB', 'n_cells', 'class'])
    if len(pair_df):
        for grp in GROUPS:
            g = pair_df[pair_df['group'] == grp]
            log(f"  {grp}: {len(g)} pair(s) in >= {HYPER_RECURRENT_CELLS} cells")
            for _, r in g.iterrows():
                log(f"    {r['geneA']} -- {r['geneB']}  {r['n_cells']} cells  [{r['class']}]")
    else:
        log(f"  No single pair reaches {HYPER_RECURRENT_CELLS} cells in any group.")
    pair_df.to_csv(os.path.join(OUTPUT_DIR, "hyper_recurrent_pairs.tsv"), sep='\t', index=False)

    # ------------------------------------------------------------------ STEP 5
    banner("STEP 5: Within-SRR control (depth / batch confound)", "-")
    # target cells per (srr, group) from the assignments; barcode is BASE-1-SRR...
    groups_df = groups_df.copy()
    groups_df['srr'] = groups_df['cell_barcode'].str.split('-').str[2]
    cells_sg = groups_df.groupby(['srr', 'group']).size().rename('n_cells').reset_index()
    jun_sg = jdf.groupby(['srr_id', 'group']).size().rename('n_junctions').reset_index() \
        .rename(columns={'srr_id': 'srr'})
    wsr = cells_sg.merge(jun_sg, on=['srr', 'group'], how='left').fillna({'n_junctions': 0})
    wsr['rate'] = wsr['n_junctions'] / wsr['n_cells']
    wsr.to_csv(os.path.join(OUTPUT_DIR, "within_srr_rates.tsv"), sep='\t', index=False)

    # For SRRs containing NORMAL + at least one tumor group, does NORMAL exceed the tumor?
    wins = losses = 0
    for srr, sub in wsr.groupby('srr'):
        rates = dict(zip(sub['group'], sub['rate']))
        if 'NORMAL' not in rates:
            continue
        tumor = [rates[g] for g in ('SBS2_HIGH', 'CNV_HIGH') if g in rates]
        if not tumor:
            continue
        if rates['NORMAL'] > max(tumor):
            wins += 1
        else:
            losses += 1
    log(f"  SRRs with NORMAL + a tumor group: {wins + losses}")
    log(f"    NORMAL per-cell rate HIGHER than both tumor groups in-sample: {wins}")
    log(f"    NORMAL not higher (depth/batch would predict this): {losses}")
    if wins + losses:
        frac = wins / (wins + losses)
        if frac >= 0.6:
            log(f"  -> NORMAL stays higher within the same sample ({frac:.0%}); the excess is not")
            log(f"     merely which SRRs NORMAL draws from. Points to locus-class artifact, not depth.")
        else:
            log(f"  -> NORMAL is NOT consistently higher within-sample ({frac:.0%}); a depth/batch")
            log(f"     component is plausible and the promiscuity filter alone may not explain it.")

    # ------------------------------------------------------------------ STEP 6
    banner(f"STEP 6: Corrected numbers at default threshold (>= {DEFAULT_THRESHOLD} partners)", "-")
    flagged = set(prom.loc[prom['n_distinct_partners'] >= DEFAULT_THRESHOLD, 'gene'])
    cls_counts = prom.loc[prom['gene'].isin(flagged), 'class'].value_counts().to_dict()
    keep = ~(jdf['geneA'].isin(flagged) | jdf['geneB'].isin(flagged))
    after = per_group_rates(jdf[keep], ncells)

    # audit: per flagged locus, how many junctions it removes per group
    audit_rows = []
    for gene in sorted(flagged):
        involved = jdf[(jdf['geneA'] == gene) | (jdf['geneB'] == gene)]
        row = {'gene': gene, 'class': locus_class(gene),
               'n_distinct_partners': int(prom.loc[prom['gene'] == gene, 'n_distinct_partners'].iloc[0]),
               'n_junctions_total': len(involved)}
        for g in GROUPS:
            row[f'{g}_junctions'] = int((involved['group'] == g).sum())
        audit_rows.append(row)
    audit_df = pd.DataFrame(audit_rows).sort_values('n_junctions_total', ascending=False)
    audit_df.to_csv(os.path.join(OUTPUT_DIR, "dropped_loci_audit.tsv"), sep='\t', index=False)

    ba_rows = []
    for g in GROUPS:
        ba_rows.append({'group': g, 'n_cells': ncells[g],
                        'junctions_before': base[g][0], 'rate_before': round(base[g][1], 3),
                        'junctions_after': after[g][0], 'rate_after': round(after[g][1], 3),
                        'removed': base[g][0] - after[g][0],
                        'pct_removed': round(100 * (base[g][0] - after[g][0]) / base[g][0], 1) if base[g][0] else 0})
    pd.DataFrame(ba_rows).to_csv(os.path.join(OUTPUT_DIR, "per_group_before_after.tsv"),
                                 sep='\t', index=False)

    log(f"  Flagged loci: {len(flagged)}  (" + ", ".join(f"{k}:{v}" for k, v in cls_counts.items()) + ")")
    for g in GROUPS:
        log(f"  {g:12s}: {base[g][0]:5d} ({base[g][1]:4.1f}/cell) -> "
            f"{after[g][0]:5d} ({after[g][1]:4.1f}/cell)   "
            f"removed {base[g][0]-after[g][0]} ({ba_rows[GROUPS.index(g)]['pct_removed']}%)")

    excess_before = base['NORMAL'][1] - tumor_base
    tumor_after = max(after['SBS2_HIGH'][1], after['CNV_HIGH'][1])
    excess_after = after['NORMAL'][1] - tumor_after

    # ------------------------------------------------------------------ SUMMARY
    banner("TEXT-READY SUMMARY (paraphrase into methods/results)")
    hla = cls_counts.get('HLA', 0)
    ribo = cls_counts.get('ribosomal', 0) + cls_counts.get('mito_ribosomal', 0)
    verdict = ("largely attributable to recurrent chimeric artifacts over hyperpolymorphic HLA and "
               "paralogous/abundant ribosomal loci rather than a genuine fusion burden"
               if excess_after < 0.5 * max(excess_before, 1e-9)
               else "only partly explained by recurrent-artifact loci; a residual difference remains")
    summary = (
        f"After the standard filter cascade the per-cell fusion rate was highest in NORMAL "
        f"({base['NORMAL'][1]:.1f}/cell) versus SBS2-HIGH ({base['SBS2_HIGH'][1]:.1f}) and "
        f"CNV-HIGH ({base['CNV_HIGH'][1]:.1f}). Removing {len(flagged)} promiscuous loci "
        f"(genes fusing to >= {DEFAULT_THRESHOLD} distinct partners; {hla} HLA, {ribo} ribosomal, "
        f"{len(flagged)-hla-ribo} other) dropped {ba_rows[2]['pct_removed']:.0f}% of NORMAL junctions "
        f"versus {ba_rows[0]['pct_removed']:.0f}% (SBS2-HIGH) and {ba_rows[1]['pct_removed']:.0f}% "
        f"(CNV-HIGH). Corrected rates were SBS2-HIGH {after['SBS2_HIGH'][1]:.1f}, CNV-HIGH "
        f"{after['CNV_HIGH'][1]:.1f}, NORMAL {after['NORMAL'][1]:.1f} per cell, reducing the "
        f"NORMAL-over-tumor excess from {excess_before:+.1f} to {excess_after:+.1f} per cell. "
        f"The elevated NORMAL fusion signal is therefore {verdict}."
    )
    log("  " + summary.replace(". ", ".\n  "))

    with open(os.path.join(OUTPUT_DIR, "fusion_artifact_audit_report.txt"), 'w') as fh:
        fh.write('\n'.join(report_lines))
    log(f"\n  Wrote breadcrumbs to: {OUTPUT_DIR}")
    log("    dropped_loci_audit.tsv, per_group_before_after.tsv, promiscuity_sensitivity.tsv,")
    log("    hyper_recurrent_pairs.tsv, within_srr_rates.tsv, fusion_artifact_audit_report.txt")


if __name__ == "__main__":
    main()
