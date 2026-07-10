#!/usr/bin/env python3
"""
Diagnostic_Selection_Robustness.py
==================================

Supplemental-methods diagnostic for the three-group selection in Step00B.

PURPOSE
  Show that the SBS2-HIGH and CNV-HIGH populations defined by the (admittedly
  complex) Step00B composite are largely recovered by far simpler selection
  rules, and that the A3A/A3B dominance contrast is an intrinsic property of
  those populations rather than something the selection criteria injected.

  Two simple methods are compared against the canonical Step00B groups:

    METHOD 2-TERM (50/50):
      SBS2-HIGH-simple = top 546 by  0.50 * SBS2_rank + 0.50 * (1 - CNV_rank)
      CNV-HIGH-simple  = top 546 by  0.50 * CNV_rank  + 0.50 * (1 - SBS2_rank)
      No A3 information of any kind enters selection.

    METHOD 3-TERM (33/33/33):
      SBS2-HIGH-simple = top 546 by  1/3 SBS2_rank + 1/3 (1 - CNV_rank)
                                     + 1/3 A3A_fraction_rank
      CNV-HIGH-simple  = top 546 by  1/3 CNV_rank  + 1/3 (1 - SBS2_rank)
                                     + 1/3 A3B_fraction_rank

  Both methods draw from a SINGLE shared pool: cancer-tissue basal cells that
  carry signature weights. Percentile ranks are computed once over that pool.
  Unlike Step00B, the CNV arm uses a SOFT inverse-SBS2 rank term rather than a
  hard SBS2 == 0 filter, so "simplest" stays honest.

  The canonical SBS2-HIGH / CNV-HIGH barcodes are READ FROM DISK
  (three_group_assignments.tsv) and never re-derived, so this script cannot
  drift from the populations used in the manuscript.

OUTPUT
  data/FIG_4/01_group_selection/selection_robustness_diagnostic/
    selection_overlap_summary.tsv     overlap (n, % , Jaccard) for every comparison
    simple_group_membership.tsv       per-cell membership flags (wide)
    simple_group_distributions.tsv    SBS2 / CNV / A3 / stemness profiles of each group

  Everything is also written to stdout for pasting back.

USAGE
  conda run -n NETWORK python Diagnostic_Selection_Robustness.py

NOTE
  Inputs (adata, weights, metric definitions) are loaded by the SAME code path
  as Step00B so the only difference between this diagnostic and the real
  selection is the scoring formula.

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime

# network_config_SC.py lives in the NETWORK_SINGLE_CELL directory, a parent of
# this TROUBLESHOOTING/ folder. Walk up from this script's own location until we
# find it, so the diagnostic runs regardless of the current working directory.
_here = os.path.dirname(os.path.abspath(__file__))
_search = _here
for _ in range(5):
    if os.path.exists(os.path.join(_search, "network_config_SC.py")):
        if _search not in sys.path:
            sys.path.insert(0, _search)
        break
    _parent = os.path.dirname(_search)
    if _parent == _search:
        break
    _search = _parent
else:
    raise ModuleNotFoundError(
        "network_config_SC.py not found in any parent of "
        f"{_here}. Run from within the NETWORK_SINGLE_CELL tree."
    )

from network_config_SC import (
    DIR_01_GROUPS,
    ADATA_FINAL_PATH, WEIGHTS_PATH,
    TARGET_CELL_TYPE,
    RANDOM_SEED, banner, log, ensure_dir
)

# =============================================================================
# CONFIGURATION
# =============================================================================

N_CELLS = 546
CANCER_SOURCE = 'head and neck squamous cell carcinoma'

# Canonical group labels as written by Step00B
CANON_FILE = os.path.join(DIR_01_GROUPS, "three_group_assignments.tsv")
LABEL_SBS2 = 'SBS2_HIGH'
LABEL_CNV = 'CNV_HIGH'

OUT_DIR = ensure_dir(os.path.join(DIR_01_GROUPS, "selection_robustness_diagnostic"))

np.random.seed(RANDOM_SEED)


# =============================================================================
# STEP 1: LOAD DATA  (mirrors Step00B.step1_load exactly)
# =============================================================================

def step1_load():
    """Load adata and weights; compute the identical per-cell metrics Step00B uses."""
    banner("STEP 1: Load Data (Step00B-identical metric computation)")

    log(f"  Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  Total cells: {adata.n_obs:,}, genes: {adata.n_vars:,}")

    basal_mask = adata.obs['final_annotation'] == TARGET_CELL_TYPE
    adata_basal = adata[basal_mask].copy()
    log(f"  Basal cells: {adata_basal.n_obs:,}")

    # SBS2 weights
    log(f"  Loading: {WEIGHTS_PATH}")
    weights_df = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    if 'SBS2' in weights_df.index:
        sbs2_map = weights_df.loc['SBS2'].to_dict()
        adata_basal.obs['SBS2'] = adata_basal.obs_names.map(
            lambda x: sbs2_map.get(x, 0.0)).astype(float)
    adata_basal.obs['has_weights'] = adata_basal.obs_names.isin(weights_df.columns)
    log(f"  Cells with weights: {adata_basal.obs['has_weights'].sum():,}")

    # A3 expression (identical to Step00B)
    gene_names = list(adata_basal.var_names)
    for gene in ['APOBEC3A', 'APOBEC3B']:
        if gene in gene_names:
            idx = gene_names.index(gene)
            vals = adata_basal.X[:, idx]
            if hasattr(vals, 'toarray'):
                vals = vals.toarray().flatten()
            else:
                vals = np.array(vals).flatten()
            adata_basal.obs[gene] = vals
        else:
            adata_basal.obs[gene] = 0.0
    adata_basal.obs['A3_sum'] = adata_basal.obs['APOBEC3A'] + adata_basal.obs['APOBEC3B']
    adata_basal.obs['A3A_fraction'] = (
        adata_basal.obs['APOBEC3A'] /
        (adata_basal.obs['APOBEC3A'] + adata_basal.obs['APOBEC3B'] + 0.01)
    )
    adata_basal.obs['A3B_fraction'] = (
        adata_basal.obs['APOBEC3B'] /
        (adata_basal.obs['APOBEC3A'] + adata_basal.obs['APOBEC3B'] + 0.01)
    )

    # CNV and CytoTRACE2
    for col in ['cnv_score', 'CytoTRACE2_Score']:
        if col in adata_basal.obs.columns:
            adata_basal.obs[col] = adata_basal.obs[col].astype(float)

    return adata_basal


# =============================================================================
# STEP 2: READ CANONICAL GROUPS FROM DISK
# =============================================================================

def step2_canonical():
    """Read the Step00B SBS2-HIGH / CNV-HIGH barcodes. Never re-derived."""
    banner("STEP 2: Read Canonical Groups (ground truth, on-disk)")

    if not os.path.exists(CANON_FILE):
        raise FileNotFoundError(
            f"Canonical assignments not found: {CANON_FILE}\n"
            f"  Run Step00B_Three_Group_Selection_and_Export.py first."
        )

    master = pd.read_csv(CANON_FILE, sep='\t')
    canon_sbs2 = set(master.loc[master['group'] == LABEL_SBS2, 'cell_barcode'])
    canon_cnv = set(master.loc[master['group'] == LABEL_CNV, 'cell_barcode'])

    log(f"  Source: {CANON_FILE}")
    log(f"  Canonical SBS2-HIGH: {len(canon_sbs2):,} cells")
    log(f"  Canonical CNV-HIGH:  {len(canon_cnv):,} cells")

    return canon_sbs2, canon_cnv


# =============================================================================
# STEP 3: BUILD SHARED POOL AND PERCENTILE RANKS
# =============================================================================

def step3_build_pool(adata_basal, canon_sbs2, canon_cnv):
    """Single shared pool: cancer-tissue basal cells with weights.

    Percentile ranks computed once over the whole pool and reused by both
    simple methods. Also confirms every canonical cell lives in this pool, so
    overlap denominators are clean.
    """
    banner("STEP 3: Build Shared Pool")

    obs = adata_basal.obs
    cancer_mask = obs['source_name'] == CANCER_SOURCE
    pool_mask = cancer_mask & obs['has_weights']
    pool = obs[pool_mask].copy()
    log(f"  Shared pool (cancer-tissue basal, has_weights): {len(pool):,} cells")

    # Sanity: canonical groups must be subsets of the pool universe
    canon_all = canon_sbs2 | canon_cnv
    in_pool = sum(1 for c in canon_all if c in pool.index)
    log(f"  Canonical cells found in pool: {in_pool:,} / {len(canon_all):,}")
    if in_pool != len(canon_all):
        log(f"  WARNING: {len(canon_all) - in_pool} canonical cells absent from pool "
            f"-- overlap denominators may be affected.")

    # Percentile ranks (Step00B uses rank(pct=True); ties averaged)
    pool['sbs2_pct'] = pool['SBS2'].rank(pct=True)
    pool['sbs2_pct_inv'] = 1.0 - pool['sbs2_pct']
    pool['cnv_pct'] = pool['cnv_score'].rank(pct=True)
    pool['cnv_pct_inv'] = 1.0 - pool['cnv_pct']
    pool['a3a_frac_pct'] = pool['A3A_fraction'].rank(pct=True)
    pool['a3b_frac_pct'] = pool['A3B_fraction'].rank(pct=True)

    return pool


# =============================================================================
# STEP 4: SIMPLE SELECTIONS
# =============================================================================

def _top_n(pool, score, n=N_CELLS):
    return set(score.sort_values(ascending=False).head(n).index)


def step4_simple_selections(pool):
    """Two-term (50/50) and three-term (33/33/33) selections."""
    banner("STEP 4: Simple Selections")

    # --- 2-TERM (50/50): no A3 information ---
    s2_sbs2_score = 0.50 * pool['sbs2_pct'] + 0.50 * pool['cnv_pct_inv']
    s2_cnv_score = 0.50 * pool['cnv_pct'] + 0.50 * pool['sbs2_pct_inv']
    s2_sbs2 = _top_n(pool, s2_sbs2_score)
    s2_cnv = _top_n(pool, s2_cnv_score)
    log(f"  2-TERM SBS2-HIGH-simple: {len(s2_sbs2):,} cells (SBS2 + low CNV)")
    log(f"  2-TERM CNV-HIGH-simple:  {len(s2_cnv):,} cells (CNV + low SBS2)")

    # --- 3-TERM (33/33/33): add A3 dominance ---
    third = 1.0 / 3.0
    s3_sbs2_score = (third * pool['sbs2_pct'] + third * pool['cnv_pct_inv']
                     + third * pool['a3a_frac_pct'])
    s3_cnv_score = (third * pool['cnv_pct'] + third * pool['sbs2_pct_inv']
                    + third * pool['a3b_frac_pct'])
    s3_sbs2 = _top_n(pool, s3_sbs2_score)
    s3_cnv = _top_n(pool, s3_cnv_score)
    log(f"  3-TERM SBS2-HIGH-simple: {len(s3_sbs2):,} cells (SBS2 + low CNV + A3A dom)")
    log(f"  3-TERM CNV-HIGH-simple:  {len(s3_cnv):,} cells (CNV + low SBS2 + A3B dom)")

    return {
        's2_sbs2': s2_sbs2, 's2_cnv': s2_cnv,
        's3_sbs2': s3_sbs2, 's3_cnv': s3_cnv,
    }


# =============================================================================
# STEP 5: OVERLAP METRICS
# =============================================================================

def _overlap(a, b):
    a, b = set(a), set(b)
    inter = len(a & b)
    union = len(a | b)
    pct = 100.0 * inter / len(a) if len(a) else 0.0   # both sets are 546 -> symmetric
    jac = inter / union if union else 0.0
    return inter, pct, jac


def step5_overlaps(simple, canon_sbs2, canon_cnv):
    """All overlap comparisons."""
    banner("STEP 5: Overlap with Canonical Groups")

    rows = []

    def add(method, comparison, a, b):
        inter, pct, jac = _overlap(a, b)
        rows.append({'method': method, 'comparison': comparison,
                     'n_intersect': inter, 'pct_overlap': round(pct, 1),
                     'jaccard': round(jac, 3)})
        log(f"  [{method:7s}] {comparison:38s} n={inter:4d}  "
            f"overlap={pct:5.1f}%  Jaccard={jac:.3f}")

    # Recovery of canonical groups
    add('2-term', 'SBS2-HIGH-simple vs canonical SBS2', simple['s2_sbs2'], canon_sbs2)
    add('2-term', 'CNV-HIGH-simple vs canonical CNV',   simple['s2_cnv'],  canon_cnv)
    add('3-term', 'SBS2-HIGH-simple vs canonical SBS2', simple['s3_sbs2'], canon_sbs2)
    add('3-term', 'CNV-HIGH-simple vs canonical CNV',   simple['s3_cnv'],  canon_cnv)

    log("")
    # Mutual overlap of the two simple arms (should be near zero)
    add('2-term', 'SBS2-HIGH-simple vs CNV-HIGH-simple', simple['s2_sbs2'], simple['s2_cnv'])
    add('3-term', 'SBS2-HIGH-simple vs CNV-HIGH-simple', simple['s3_sbs2'], simple['s3_cnv'])

    log("")
    # 2-term -> 3-term membership shift (does adding A3 reshuffle the group?)
    add('2->3',  'SBS2 arm: 2-term vs 3-term', simple['s2_sbs2'], simple['s3_sbs2'])
    add('2->3',  'CNV arm: 2-term vs 3-term',  simple['s2_cnv'],  simple['s3_cnv'])

    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(OUT_DIR, "selection_overlap_summary.tsv"),
                   sep='\t', index=False)
    log(f"\n  Wrote: selection_overlap_summary.tsv")
    return summary


# =============================================================================
# STEP 6: DISTRIBUTIONS (the intrinsic-dominance argument)
# =============================================================================

def step6_distributions(adata_basal, simple, canon_sbs2, canon_cnv):
    """Profile each group. Key point: the 2-term groups carry no A3 input, so
    their A3A/A3B dominance is intrinsic to the populations, not selected for.
    """
    banner("STEP 6: Group Distributions (intrinsic A3 dominance check)")

    obs = adata_basal.obs
    groups = {
        'canonical_SBS2': canon_sbs2,
        'canonical_CNV':  canon_cnv,
        '2term_SBS2':     simple['s2_sbs2'],
        '2term_CNV':      simple['s2_cnv'],
        '3term_SBS2':     simple['s3_sbs2'],
        '3term_CNV':      simple['s3_cnv'],
    }
    metrics = ['SBS2', 'cnv_score', 'A3A_fraction', 'A3B_fraction',
               'A3_sum', 'CytoTRACE2_Score']
    metrics = [m for m in metrics if m in obs.columns]

    rows = []
    for gname, cells in groups.items():
        sub = obs.loc[obs.index.isin(cells)]
        row = {'group': gname, 'n': len(sub)}
        for m in metrics:
            vals = sub[m].dropna()
            row[f'{m}_mean'] = round(float(vals.mean()), 4) if len(vals) else np.nan
            row[f'{m}_median'] = round(float(vals.median()), 4) if len(vals) else np.nan

        # Dominance call among A3-expressing cells only (clean A3A vs A3B split)
        a3pos = sub[sub['A3_sum'] > 0]
        n_a3pos = len(a3pos)
        if n_a3pos > 0:
            a3a_dom = (a3pos['A3A_fraction'] > 0.5).sum()
            row['n_A3pos'] = n_a3pos
            row['pct_A3A_dominant'] = round(100.0 * a3a_dom / n_a3pos, 1)
            row['pct_A3B_dominant'] = round(100.0 * (n_a3pos - a3a_dom) / n_a3pos, 1)
        else:
            row['n_A3pos'] = 0
            row['pct_A3A_dominant'] = np.nan
            row['pct_A3B_dominant'] = np.nan
        rows.append(row)

    dist = pd.DataFrame(rows)
    dist.to_csv(os.path.join(OUT_DIR, "simple_group_distributions.tsv"),
                sep='\t', index=False)

    # Console-friendly view of the headline columns
    show = ['group', 'n', 'SBS2_mean', 'cnv_score_mean',
            'A3A_fraction_mean', 'A3B_fraction_mean',
            'n_A3pos', 'pct_A3A_dominant', 'pct_A3B_dominant']
    show = [c for c in show if c in dist.columns]
    log("\n" + dist[show].to_string(index=False))
    log(f"\n  Wrote: simple_group_distributions.tsv")
    return dist


# =============================================================================
# STEP 7: MEMBERSHIP TABLE
# =============================================================================

def step7_membership(simple, canon_sbs2, canon_cnv):
    """Per-cell wide membership flags for every selected cell."""
    banner("STEP 7: Membership Table")

    all_cells = (set(canon_sbs2) | set(canon_cnv)
                 | simple['s2_sbs2'] | simple['s2_cnv']
                 | simple['s3_sbs2'] | simple['s3_cnv'])
    all_cells = sorted(all_cells)

    df = pd.DataFrame(index=all_cells)
    df.index.name = 'cell_barcode'
    df['canonical_SBS2'] = [c in canon_sbs2 for c in all_cells]
    df['canonical_CNV'] = [c in canon_cnv for c in all_cells]
    df['s2_SBS2'] = [c in simple['s2_sbs2'] for c in all_cells]
    df['s2_CNV'] = [c in simple['s2_cnv'] for c in all_cells]
    df['s3_SBS2'] = [c in simple['s3_sbs2'] for c in all_cells]
    df['s3_CNV'] = [c in simple['s3_cnv'] for c in all_cells]

    df.to_csv(os.path.join(OUT_DIR, "simple_group_membership.tsv"), sep='\t')
    log(f"  Union of all selected cells: {len(all_cells):,}")
    log(f"  Wrote: simple_group_membership.tsv")
    return df


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("Selection Robustness Diagnostic")
    log(f"  Start: {t0}")
    log(f"  Output: {OUT_DIR}")

    adata_basal = step1_load()
    canon_sbs2, canon_cnv = step2_canonical()
    pool = step3_build_pool(adata_basal, canon_sbs2, canon_cnv)
    simple = step4_simple_selections(pool)
    step5_overlaps(simple, canon_sbs2, canon_cnv)
    step6_distributions(adata_basal, simple, canon_sbs2, canon_cnv)
    step7_membership(simple, canon_sbs2, canon_cnv)

    banner("Diagnostic COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")


if __name__ == "__main__":
    main()
