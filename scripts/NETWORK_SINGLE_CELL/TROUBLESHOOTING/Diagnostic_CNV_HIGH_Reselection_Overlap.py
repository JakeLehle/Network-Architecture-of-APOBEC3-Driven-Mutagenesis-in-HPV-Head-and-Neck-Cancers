#!/usr/bin/env python3
"""
Diagnostic_CNV_HIGH_Reselection_Overlap.py
===========================================

DIAGNOSTIC (read-only): compare candidate CNV-HIGH selection schemes against
the CNV-HIGH group used in the current results, and measure how far the
selection moves in cell-barcode space.

WHY
---
The original CNV-HIGH composite (Step00B) included an HPV16 late-gene (L1/L2)
term, which forward-references the virus analysis into group selection and was
meant to be diagnostic only. This script evaluates two late-gene-free
replacements before committing to a full network rerun:

  OPTION 1 (preferred; mirrors the SBS2-HIGH composite)
    CNV rank ............... 0.40   high CNV, the productive hallmark
    total A3 range match ... 0.20   A3A+A3B near the SBS2-HIGH mean
    stemness rank .......... 0.20   high CytoTRACE2; CNV-HIGH is more stem-like
                                    via p53 suppression and G2/M checkpoint loss
    A3B dominance rank ..... 0.20   A3B/(A3A+A3B); counterpart to the A3A
                                    dominance term in SBS2-HIGH

  OPTION 2 (equal weighting; sensitivity check)
    0.25 on each of the four terms above.

Both options DROP the old anti-correlated signature-profile term and the old
SBS2=0 scoring weight (SBS2=0 remains a hard pool filter), and neither uses the
late-gene term. SBS2-HIGH and NORMAL are unchanged and are not touched here.

WHAT IT DOES
------------
  1. Loads adata + signature weights and recomputes the per-cell metrics exactly
     as Step00B (adding A3B_fraction).
  2. Reads the OLD SBS2-HIGH and OLD CNV-HIGH barcodes from the existing
     three_group_assignments.tsv (the actual cells behind the current results).
  3. Rebuilds the CNV pool exactly as production (SBS2==0, cancer tissue, not
     already in SBS2-HIGH, has weights).
  4. Scores the pool under Option 1 and Option 2 and takes the top N cells,
     where N is the size of the old CNV-HIGH group (fair, equal-size overlap).
  5. Reports barcode overlap (OPT1 vs OLD, OPT2 vs OLD, OPT1 vs OPT2) as both
     percent of N and Jaccard, plus group-mean metrics for SBS2-HIGH, OLD,
     OPT1, and OPT2.

WHAT IT DOES NOT DO
-------------------
  - Does not modify any pipeline output or assignment file.
  - Does not run the networks. This is a selection-sensitivity check only.

OUTPUT (-> data/FIG_4/01_group_selection/CNV_HIGH_RESELECTION_DIAGNOSTIC/)
  CNV_HIGH_reselection_overlap_summary.tsv
  CNV_HIGH_reselection_group_means.tsv
  CNV_HIGH_reselection_membership.tsv

Usage:
  conda run -n NETWORK python Diagnostic_CNV_HIGH_Reselection_Overlap.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc

from network_config_SC import (
    DIR_01_GROUPS, ADATA_FINAL_PATH, WEIGHTS_PATH,
    TARGET_CELL_TYPE, banner, log, ensure_dir
)

# -----------------------------------------------------------------------------
# CONFIG (mirrors the Step00B locals so the pool matches production exactly)
# -----------------------------------------------------------------------------
PROJECT_ROOT  = "/master/jlehle/WORKING/2026_NMF_PAPER"
CANCER_SOURCE = 'head and neck squamous cell carcinoma'
HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_6", "03_hpv16_genome",
                             "per_cell_hpv16_gene_counts.tsv")
MASTER_PATH   = os.path.join(DIR_01_GROUPS, "three_group_assignments.tsv")
OUT_DIR       = os.path.join(DIR_01_GROUPS, "CNV_HIGH_RESELECTION_DIAGNOSTIC")

# Candidate weight schemes. Terms: cnv, a3range, stem, a3bdom (sum to 1.0 each).
OPTIONS = {
    "OPT1": {"cnv": 0.40, "a3range": 0.20, "stem": 0.20, "a3bdom": 0.20},
    "OPT2": {"cnv": 0.25, "a3range": 0.25, "stem": 0.25, "a3bdom": 0.25},
}

# late_fraction is carried only as a reference column (to show the old group
# leaned on late genes and the new ones do not); it never enters scoring here.
METRICS = ['cnv_score', 'A3_sum', 'A3A_fraction', 'A3B_fraction',
           'CytoTRACE2_Score', 'SBS2', 'late_fraction']


def load_basal():
    """Load adata, subset to basal, recompute per-cell metrics as in Step00B."""
    banner("LOAD ADATA + METRICS")
    log(f"  adata: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)

    if 'final_annotation' not in adata.obs.columns:
        log("  ERROR: 'final_annotation' not in adata.obs"); sys.exit(1)
    basal = adata[adata.obs['final_annotation'] == TARGET_CELL_TYPE].copy()
    if basal.n_obs == 0:
        log(f"  ERROR: no cells with final_annotation == '{TARGET_CELL_TYPE}'"); sys.exit(1)
    log(f"  basal cells: {basal.n_obs:,}")

    # SBS2 weights
    weights_df = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    if 'SBS2' not in weights_df.index:
        log("  ERROR: SBS2 not found in weights matrix"); sys.exit(1)
    sbs2_map = weights_df.loc['SBS2'].to_dict()
    basal.obs['SBS2'] = basal.obs_names.map(lambda x: sbs2_map.get(x, 0.0)).astype(float)
    basal.obs['has_weights'] = basal.obs_names.isin(weights_df.columns)

    # A3 expression
    gene_names = list(basal.var_names)
    for gene in ['APOBEC3A', 'APOBEC3B']:
        if gene in gene_names:
            vals = basal.X[:, gene_names.index(gene)]
            vals = vals.toarray().flatten() if hasattr(vals, 'toarray') else np.asarray(vals).flatten()
            basal.obs[gene] = vals
        else:
            log(f"  WARNING: {gene} not in var_names, set to 0"); basal.obs[gene] = 0.0
    basal.obs['A3_sum'] = basal.obs['APOBEC3A'] + basal.obs['APOBEC3B']
    denom = basal.obs['APOBEC3A'] + basal.obs['APOBEC3B'] + 0.01
    basal.obs['A3A_fraction'] = basal.obs['APOBEC3A'] / denom
    basal.obs['A3B_fraction'] = basal.obs['APOBEC3B'] / denom   # NEW: A3B dominance

    # CNV + stemness
    for col in ['cnv_score', 'CytoTRACE2_Score']:
        if col not in basal.obs.columns:
            log(f"  ERROR: required column '{col}' not in adata.obs"); sys.exit(1)
        basal.obs[col] = basal.obs[col].astype(float)

    # HPV late-gene fraction (reference / diagnostic only, not scored)
    if os.path.exists(HPV_GENE_PATH):
        hpv = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
        for col in ['L1', 'L2', 'total_hpv16_genome_reads']:
            if col in hpv.columns:
                basal.obs[col] = hpv[col].reindex(basal.obs_names).fillna(0).astype(float)
        if 'L1' in hpv.columns:
            l1 = basal.obs.get('L1', 0.0)
            l2 = basal.obs.get('L2', 0.0)
            tot = basal.obs.get('total_hpv16_genome_reads', 0.0)
            basal.obs['late_fraction'] = (l1 + l2) / (tot + 0.5)
            log(f"  HPV late-gene data present "
                f"({int((basal.obs['late_fraction'] > 0).sum()):,} cells with signal)")
        else:
            basal.obs['late_fraction'] = 0.0
    else:
        basal.obs['late_fraction'] = 0.0
        log("  HPV gene file not found; late_fraction set to 0 (reference only)")

    return basal


def load_old_groups():
    """Read OLD SBS2-HIGH and OLD CNV-HIGH barcodes from the master assignment."""
    banner("LOAD OLD GROUP ASSIGNMENTS")
    if not os.path.exists(MASTER_PATH):
        log(f"  ERROR: {MASTER_PATH} not found. Run Step00B first."); sys.exit(1)
    master = pd.read_csv(MASTER_PATH, sep='\t')
    old_high = master.loc[master['group'] == 'SBS2_HIGH', 'cell_barcode'].tolist()
    old_cnv  = master.loc[master['group'] == 'CNV_HIGH',  'cell_barcode'].tolist()
    log(f"  OLD SBS2-HIGH: {len(old_high):,} cells")
    log(f"  OLD CNV-HIGH:  {len(old_cnv):,} cells")
    if len(old_cnv) == 0 or len(old_high) == 0:
        log("  ERROR: empty SBS2-HIGH or CNV-HIGH group in master file"); sys.exit(1)
    return old_high, old_cnv


def build_pool(obs, old_high):
    """Rebuild the CNV pool exactly as Step00B step3_select_cnv."""
    banner("BUILD CNV POOL")
    if 'source_name' not in obs.columns:
        log("  ERROR: 'source_name' not in adata.obs"); sys.exit(1)
    cancer    = obs['source_name'] == CANCER_SOURCE
    sbs2_zero = obs['SBS2'] == 0
    not_high  = ~obs.index.isin(old_high)
    pool_mask = cancer & sbs2_zero & not_high & obs['has_weights']
    pool = obs.index[pool_mask].tolist()
    log(f"  pool (SBS2=0, cancer, not SBS2-HIGH, has weights): {len(pool):,} cells")
    return pool


def score_pool(obs, pool, old_high):
    """Per-cell term scores over the pool, then Option 1/2 composite scores."""
    banner("SCORE POOL UNDER OPTION 1 AND OPTION 2")

    # A3-range reference from the (unchanged) SBS2-HIGH group
    high_a3 = obs.loc[obs.index.isin(old_high), 'A3_sum']
    mu, sd = float(high_a3.mean()), float(high_a3.std())
    log(f"  SBS2-HIGH A3_sum reference: mean={mu:.4f}, std={sd:.4f}")

    df = pd.DataFrame(index=pd.Index(pool, name='cell_barcode'))
    # high CNV (percentile rank within pool)
    df['cnv'] = obs.loc[pool, 'cnv_score'].rank(pct=True).values
    # total A3 near the SBS2-HIGH mean (Gaussian proximity)
    a3 = obs.loc[pool, 'A3_sum']
    if sd > 0:
        z = (a3 - mu) / sd
        df['a3range'] = np.exp(-0.5 * z ** 2).values
    else:
        df['a3range'] = (a3 == mu).astype(float).values
    # high stemness (percentile rank)
    df['stem'] = obs.loc[pool, 'CytoTRACE2_Score'].rank(pct=True).values
    # high A3B dominance (percentile rank)
    df['a3bdom'] = obs.loc[pool, 'A3B_fraction'].rank(pct=True).values

    for name, w in OPTIONS.items():
        df[name] = (w['cnv'] * df['cnv'] + w['a3range'] * df['a3range'] +
                    w['stem'] * df['stem'] + w['a3bdom'] * df['a3bdom'])
        log(f"  {name} weights: {w}")
    return df


def top_n(df, score_col, n):
    return df.sort_values(score_col, ascending=False).head(n).index.tolist()


def overlap_stats(a, b):
    sa, sb = set(a), set(b)
    inter = len(sa & sb)
    union = len(sa | sb)
    pct = 100.0 * inter / len(sa) if len(sa) else 0.0
    jacc = inter / union if union else 0.0
    return inter, pct, jacc


def main():
    banner("CNV-HIGH RESELECTION OVERLAP DIAGNOSTIC")
    ensure_dir(OUT_DIR)

    basal = load_basal()
    obs = basal.obs
    old_high, old_cnv = load_old_groups()
    n_target = len(old_cnv)
    log(f"\n  Target group size (matched to OLD CNV-HIGH): {n_target}")

    pool = build_pool(obs, old_high)
    if len(pool) < n_target:
        log(f"  ERROR: pool ({len(pool)}) smaller than target ({n_target})"); sys.exit(1)

    df = score_pool(obs, pool, old_high)
    opt1 = top_n(df, 'OPT1', n_target)
    opt2 = top_n(df, 'OPT2', n_target)

    # ---- Overlap ----
    banner("BARCODE OVERLAP")
    comparisons = [
        ("OPT1 vs OLD",  opt1, old_cnv),
        ("OPT2 vs OLD",  opt2, old_cnv),
        ("OPT1 vs OPT2", opt1, opt2),
    ]
    rows = []
    log(f"  {'comparison':16s} {'shared':>8s} {'of':>6s} {'pct':>8s} {'jaccard':>9s}")
    for label, a, b in comparisons:
        inter, pct, jacc = overlap_stats(a, b)
        log(f"  {label:16s} {inter:8d} {n_target:6d} {pct:7.1f}% {jacc:9.3f}")
        rows.append({'comparison': label, 'shared': inter, 'group_size': n_target,
                     'pct_overlap': round(pct, 2), 'jaccard': round(jacc, 4)})
    pd.DataFrame(rows).to_csv(
        os.path.join(OUT_DIR, "CNV_HIGH_reselection_overlap_summary.tsv"),
        sep='\t', index=False)

    # ---- Group means ----
    banner("GROUP-MEAN METRICS")
    groups = [("SBS2_HIGH(ref)", old_high), ("OLD_CNV", old_cnv),
              ("OPT1_CNV", opt1), ("OPT2_CNV", opt2)]
    mean_rows = []
    header = f"  {'metric':20s}" + "".join(f"{g:>16s}" for g, _ in groups)
    log(header)
    for m in METRICS:
        line = f"  {m:20s}"
        rec = {'metric': m}
        for gname, cells in groups:
            val = obs.loc[obs.index.isin(cells), m].mean()
            line += f"{val:16.4f}"
            rec[gname] = round(float(val), 4)
        log(line)
        mean_rows.append(rec)
    pd.DataFrame(mean_rows).to_csv(
        os.path.join(OUT_DIR, "CNV_HIGH_reselection_group_means.tsv"),
        sep='\t', index=False)

    # ---- Per-cell membership (union of old / opt1 / opt2) ----
    union = sorted(set(old_cnv) | set(opt1) | set(opt2))
    mem = pd.DataFrame(index=pd.Index(union, name='cell_barcode'))
    mem['in_old']  = mem.index.isin(old_cnv)
    mem['in_opt1'] = mem.index.isin(opt1)
    mem['in_opt2'] = mem.index.isin(opt2)
    for col in ['OPT1', 'OPT2']:
        mem[f'{col}_score'] = df[col].reindex(mem.index)
    for m in METRICS:
        mem[m] = obs.loc[obs.index.isin(union), m].reindex(mem.index)
    mem.to_csv(os.path.join(OUT_DIR, "CNV_HIGH_reselection_membership.tsv"), sep='\t')

    banner("DONE")
    log(f"  Outputs in: {OUT_DIR}")
    log("    CNV_HIGH_reselection_overlap_summary.tsv")
    log("    CNV_HIGH_reselection_group_means.tsv")
    log("    CNV_HIGH_reselection_membership.tsv")
    log("\n  No pipeline files were modified; no networks were run.")


if __name__ == "__main__":
    main()
