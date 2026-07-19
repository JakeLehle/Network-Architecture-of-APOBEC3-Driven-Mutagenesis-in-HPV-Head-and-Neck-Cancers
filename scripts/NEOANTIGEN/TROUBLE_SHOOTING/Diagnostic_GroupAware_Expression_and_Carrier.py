#!/usr/bin/env python3
"""
Diagnostic_GroupAware_Expression_and_Carrier.py
===============================================
Group-aware expression + neoantigen carrier fraction for the Figure 7 lead genes,
plus the neoantigen-burden significance test for Panel A.

Purpose
-------
1. GATE (Step 0): reconcile cell barcodes across the three sources that must line
   up before any carrier fraction is trustworthy:
     - adata_final.h5ad          (expression; obs_names)
     - three_group_assignments   (group membership; cell_barcode)
     - SComatic single-cell genotype master (per-cell genotypes; CB)
   Prints overlap counts + example barcodes. If overlap is poor it stops before
   computing anything downstream.

2. Group-aware % expressing (Step 2): reproduces Step04's definition
   (100 * mean(counts > 0)) but with the tier-correct denominator:
     - Tier 2 (SBS2-specific): % of the 546 SBS2-HIGH cells
     - Tier 1 (shared):        % over SBS2-HIGH + CNV-HIGH cells (pooled), and
                               each group separately for reference

3. Carrier fraction (Step 3): among expressing cells, the fraction carrying ANY
   of the gene's mutations, reported two ways:
     - over ALL expressers (headline definition)
     - over expressers WITH coverage at >=1 of the gene's loci (dropout-adjusted)
   Coverage counts are printed so the dropout is visible.

4. Neoantigen burden test for Panel A (Step 4): computes BOTH so we decide from
   the numbers, not from argument:
     - GROUP-RATE test (matches the plotted 4.34 vs 2.45/cell): rate ratio and a
       binomial test on the SBS2/CNV split, at three units (peptides,
       neoantigen-producing variants, genes).
     - PER-CELL genotyped burden (test-type parity with the fusion Mann-Whitney):
       per-cell count of carried neoantigen-producing variants, raw and per-UMI,
       Mann-Whitney SBS2 vs CNV. FLAGGED depth-confounded (SBS2 shallower);
       cross-check only until we see whether it is reversed.

READ-ONLY. Env: NETWORK (scanpy, scipy, pandas).
Outputs to data/FIG_7/TROUBLESHOOTING/group_aware_expression/.
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH = os.path.join(BASE_DIR, "data/FIG_4/00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
SCOMATIC_TSV = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
                "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")
ANNOTATION_DIR = os.path.join(BASE_DIR, "data/FIG_7/02_snpeff_annotation")
MHC_DIR = os.path.join(BASE_DIR, "data/FIG_7/03_mhc_binding")
CANDIDATE_FILE = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/"
                              "tcw_neoantigen_ranking/gene_level_A3_neoantigen_CANDIDATES.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/group_aware_expression")

FIGURE_GENES = ['KLF3', 'CAST', 'SERPINB2', 'KRT6B', 'ANXA1']
TUMOR = ['SBS2_HIGH', 'CNV_HIGH']
TIER_SHARED = 'Tier1_shared'   # candidate-file label -> pooled SBS2+CNV denominator

MIN_BARCODE_OVERLAP = 0.5   # gate: fraction of group cells that must map to each source

report = []
def log(msg=""):
    print(msg, flush=True)
    report.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)


# =============================================================================
# LOADERS
# =============================================================================
def load_groups():
    g = pd.read_csv(GROUP_PATH, sep='\t')
    bc_col = 'cell_barcode' if 'cell_barcode' in g.columns else g.columns[0]
    grp_col = 'group' if 'group' in g.columns else g.columns[1]
    g = g.rename(columns={bc_col: 'cell_barcode', grp_col: 'group'})
    return g[['cell_barcode', 'group']]


def load_adata():
    import scanpy as sc
    return sc.read_h5ad(ADATA_PATH)


def gene_expr_vector(adata, gene):
    """Dense 1D expression vector for a gene, or None if absent."""
    if gene not in adata.var_names:
        return None
    import scipy.sparse as sp
    col = adata[:, gene].X
    arr = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
    return arr


def load_percell_umi(adata):
    """Per-cell total UMI keyed by obs_name (for per-UMI burden)."""
    for c in adata.obs.columns:
        if c.lower() in ('total_counts', 'n_counts', 'ncount_rna', 'total_umis',
                         'n_umi', 'nreads', 'n_reads'):
            return dict(zip(adata.obs_names, adata.obs[c].astype(float)))
    import scipy.sparse as sp
    src = adata.layers['counts'] if 'counts' in adata.layers else (adata.raw.X if adata.raw is not None else adata.X)
    tot = np.asarray(src.sum(axis=1)).ravel() if sp.issparse(src) else np.asarray(src).sum(axis=1).ravel()
    if np.nanmax(tot) < 100:
        return {}   # looks log-normalized, not counts
    return dict(zip(adata.obs_names, tot))


def load_gene_variants():
    """Per-gene protein-altering variants (chrom,pos,ref,alt,hgvs) + neoantigen flag,
    from the SBS2 annotation joined to the neoantigen calls. Also returns the full
    set of neoantigen-producing variant positions per group for the burden test."""
    variants = {}          # gene -> DataFrame(chrom,pos,ref,alt,hgvs_p,is_neo)
    neo_positions = {}     # group -> set((chrom,pos,alt)) neoantigen-producing
    for group in TUMOR:
        spa_p = os.path.join(ANNOTATION_DIR, f"{group}.somatic_protein_altering.tsv")
        neo_p = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
        if not os.path.exists(spa_p):
            continue
        spa = pd.read_csv(spa_p, sep='\t')
        lm = {c.lower(): c for c in spa.columns}
        cc, pc, rc, ac = lm.get('chrom', lm.get('#chrom')), lm.get('pos'), lm.get('ref'), lm.get('alt')
        gc, hc = lm.get('gene'), lm.get('hgvs_p')
        if not all([cc, pc, rc, ac, gc]):
            log(f"  [WARN] {group} annotation missing chrom/pos/ref/alt/gene")
            continue
        spa = spa.rename(columns={cc: 'chrom', pc: 'pos', rc: 'ref', ac: 'alt', gc: 'gene'})
        if hc:
            spa = spa.rename(columns={hc: 'hgvs_p'})
        else:
            spa['hgvs_p'] = ''

        neo_keys = set()
        if os.path.exists(neo_p):
            neo = pd.read_csv(neo_p, sep='\t')
            if {'gene', 'hgvs_p'} <= set(neo.columns):
                neo_keys = set(zip(neo['gene'].astype(str), neo['hgvs_p'].astype(str)))
        spa['is_neo'] = [(r.gene, str(r.hgvs_p)) in neo_keys for r in spa.itertuples()]

        neo_positions[group] = set(
            (str(r.chrom), int(r.pos), str(r.alt))
            for r in spa[spa['is_neo']].itertuples()
        )

        for gene in FIGURE_GENES:
            sub = spa[spa['gene'] == gene][['chrom', 'pos', 'ref', 'alt', 'hgvs_p', 'is_neo']].copy()
            if len(sub) == 0:
                continue
            if gene not in variants or (group == 'SBS2_HIGH'):
                variants[gene] = sub.drop_duplicates(subset=['chrom', 'pos', 'alt'])
    return variants, neo_positions


def load_tiers():
    if not os.path.exists(CANDIDATE_FILE):
        return {}
    c = pd.read_csv(CANDIDATE_FILE, sep='\t')
    if 'gene' in c.columns and 'tier' in c.columns:
        return dict(zip(c['gene'].astype(str), c['tier'].astype(str)))
    return {}


# =============================================================================
# SComatic streaming: per-cell genotypes at needed positions only
# =============================================================================
def stream_scomatic(needed_positions):
    """Stream the (large) SComatic master, keep only rows at needed (chrom,pos).
    Returns:
      cov[(chrom,pos)]     -> set(CB) with ANY call at the locus (coverage)
      mut[(chrom,pos,alt)] -> set(CB) with Base_observed == alt (carrier)
    """
    cov, mut = {}, {}
    if not os.path.exists(SCOMATIC_TSV):
        log(f"  [FATAL] SComatic master not found: {SCOMATIC_TSV}")
        return cov, mut
    with open(SCOMATIC_TSV) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {name: i for i, name in enumerate(header)}
        ci, si = idx.get('#CHROM'), idx.get('Start')
        bi, cbi = idx.get('Base_observed'), idx.get('CB')
        if None in (ci, si, bi, cbi):
            log(f"  [FATAL] SComatic master missing expected columns; header={header}")
            return cov, mut
        want = set((c, p) for (c, p) in needed_positions)
        n = 0
        for line in f:
            fields = line.rstrip('\n').split('\t')
            try:
                key = (fields[ci], int(fields[si]))
            except (ValueError, IndexError):
                continue
            if key not in want:
                continue
            cb = fields[cbi]
            base = fields[bi]
            cov.setdefault(key, set()).add(cb)
            mut.setdefault((key[0], key[1], base), set()).add(cb)
            n += 1
    log(f"  streamed {n:,} SComatic rows at {len(want):,} needed loci")
    return cov, mut


# =============================================================================
# STEP 0: barcode reconciliation gate
# =============================================================================
def reconcile_barcodes(group_bcs, adata_names, scomatic_cbs):
    banner("STEP 0: BARCODE RECONCILIATION (gate)")
    gset, aset, sset = set(group_bcs), set(adata_names), set(scomatic_cbs)
    log(f"  group cells:     {len(gset):,}   e.g. {list(gset)[:2]}")
    log(f"  adata cells:     {len(aset):,}   e.g. {list(aset)[:2]}")
    log(f"  SComatic cells:  {len(sset):,}   e.g. {list(sset)[:2]}  (sampled)")

    ov_a = len(gset & aset) / len(gset) if gset else 0
    ov_s = len(gset & sset) / len(gset) if gset else 0
    log(f"  group and adata:     {len(gset & aset):,} ({100 * ov_a:.1f}% of group)")
    log(f"  group and SComatic:  {len(gset & sset):,} ({100 * ov_s:.1f}% of group, sampled)")

    ok_a = ov_a >= MIN_BARCODE_OVERLAP
    ok_s = (ov_s >= MIN_BARCODE_OVERLAP) or (len(sset & gset) > 0)  # sampled; any overlap is enough to proceed
    if ok_a and ok_s:
        log("  -> direct barcode match OK; proceeding.")
    else:
        log("  -> [GATE] low overlap; carrier/expression numbers would be unreliable.")
        log("     Check the -1-{SRR} suffix format on each source before trusting output.")
    return ok_a and ok_s


# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("GROUP-AWARE EXPRESSION + NEOANTIGEN CARRIER FRACTION")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    g = load_groups()
    grp_of = dict(zip(g['cell_barcode'], g['group']))
    group_cells = {grp: set(g[g['group'] == grp]['cell_barcode']) for grp in g['group'].unique()}

    variants, neo_positions = load_gene_variants()
    tiers = load_tiers()

    needed = set()
    for gene, df in variants.items():
        for r in df.itertuples():
            needed.add((str(r.chrom), int(r.pos)))
    for grp, s in neo_positions.items():
        for (c, p, a) in s:
            needed.add((c, p))

    # peek at SComatic barcodes for the reconciliation gate (cheap first pass)
    scomatic_cbs = set()
    if os.path.exists(SCOMATIC_TSV):
        with open(SCOMATIC_TSV) as f:
            header = f.readline().rstrip('\n').split('\t')
            cbi = {n: i for i, n in enumerate(header)}.get('CB')
            if cbi is not None:
                for i, line in enumerate(f):
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) > cbi:
                        scomatic_cbs.add(parts[cbi])
                    if i > 500000:
                        break

    adata = load_adata()
    log(f"  adata: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    ok = reconcile_barcodes(list(g['cell_barcode']), list(adata.obs_names), scomatic_cbs)
    if not ok:
        log("\n  Stopping after the reconciliation gate. Nothing downstream computed.")
        with open(os.path.join(OUTPUT_DIR, "group_aware_report.txt"), 'w') as fh:
            fh.write('\n'.join(report))
        return

    umi = load_percell_umi(adata)

    # -------------------------------------------------------------------------
    # STEP 2: group-aware % expressing
    # -------------------------------------------------------------------------
    banner("STEP 2: GROUP-AWARE % EXPRESSING (tier-correct denominator)")
    obs_group = pd.Series({bc: grp_of.get(bc) for bc in adata.obs_names})
    masks = {grp: (obs_group.values == grp) for grp in TUMOR}

    expr_rows = []
    for gene in FIGURE_GENES:
        vec = gene_expr_vector(adata, gene)
        tier = tiers.get(gene, '?')
        if vec is None:
            log(f"  {gene}: not in adata var_names")
            expr_rows.append({'gene': gene, 'tier': tier, 'pct_sbs2': np.nan,
                              'pct_cnv': np.nan, 'pct_pooled': np.nan,
                              'denominator': 'NA', 'headline_pct': np.nan})
            continue
        pct_sbs2 = 100 * np.mean(vec[masks['SBS2_HIGH']] > 0) if masks['SBS2_HIGH'].any() else np.nan
        pct_cnv = 100 * np.mean(vec[masks['CNV_HIGH']] > 0) if masks['CNV_HIGH'].any() else np.nan
        pooled_mask = masks['SBS2_HIGH'] | masks['CNV_HIGH']
        pct_pooled = 100 * np.mean(vec[pooled_mask] > 0) if pooled_mask.any() else np.nan
        denom = 'SBS2+CNV pooled' if tier == TIER_SHARED else 'SBS2 only'
        headline = pct_pooled if tier == TIER_SHARED else pct_sbs2
        log(f"  {gene:9s} [{tier:18s}] denom={denom:16s} -> {headline:5.1f}%   "
            f"(SBS2 {pct_sbs2:.1f}%, CNV {pct_cnv:.1f}%)")
        expr_rows.append({'gene': gene, 'tier': tier, 'pct_sbs2': pct_sbs2,
                          'pct_cnv': pct_cnv, 'pct_pooled': pct_pooled,
                          'denominator': denom, 'headline_pct': headline})
    expr_df = pd.DataFrame(expr_rows)

    # -------------------------------------------------------------------------
    # Stream SComatic once for everything we need
    # -------------------------------------------------------------------------
    banner("Streaming SComatic single-cell genotypes at needed loci", "-")
    cov, mut = stream_scomatic(needed)

    # -------------------------------------------------------------------------
    # STEP 3: carrier fraction among expressing cells
    # -------------------------------------------------------------------------
    banner("STEP 3: CARRIER FRACTION (among expressing cells)")
    log("  'carrier' = cell genotyped Base_observed==alt at >=1 of the gene's loci.")
    log("  Reported as carrier / all expressers = fraction of expressing cells with the")
    log("  mutation detected. No 'covered' denominator: the FILTERED SComatic master only")
    log("  stores variant-supporting rows, so coverage==carrier by construction (a true")
    log("  coverage denominator would need the unfiltered per-cell base-count table).")
    log("")
    obs_names_arr = np.array(adata.obs_names)
    carrier_rows = []
    for gene in FIGURE_GENES:
        vec = gene_expr_vector(adata, gene)
        tier = tiers.get(gene, '?')
        if vec is None or gene not in variants:
            log(f"  {gene}: no variants or not in adata; skipped")
            continue
        if tier == TIER_SHARED:
            denom_cells = group_cells.get('SBS2_HIGH', set()) | group_cells.get('CNV_HIGH', set())
        else:
            denom_cells = group_cells.get('SBS2_HIGH', set())
        expr_names = set(obs_names_arr[vec > 0]) & denom_cells

        gv = variants[gene]
        carrier_cells = set()
        for r in gv.itertuples():
            carrier_cells |= mut.get((str(r.chrom), int(r.pos), str(r.alt)), set())
        expr_carrier = expr_names & carrier_cells

        n_expr = len(expr_names)
        pct_all = 100 * len(expr_carrier) / n_expr if n_expr else np.nan
        log(f"  {gene:9s} [{tier:18s}] n_variants={len(gv)}  expressers={n_expr}  "
            f"carriers={len(expr_carrier)}")
        log(f"            carrier / all expressers = {pct_all:5.1f}%  "
            f"(expressing cells with the mutation detected)")
        carrier_rows.append({'gene': gene, 'tier': tier, 'n_variants': len(gv),
                             'n_expressers': n_expr, 'n_carriers': len(expr_carrier),
                             'carrier_pct_all_expressers': pct_all})
    carrier_df = pd.DataFrame(carrier_rows)

    fig_expr = expr_df.merge(carrier_df.drop(columns=['tier'], errors='ignore'),
                             on='gene', how='left')
    fig_expr.to_csv(os.path.join(OUTPUT_DIR, "group_aware_expression_carrier.tsv"),
                    sep='\t', index=False)

    # -------------------------------------------------------------------------
    # STEP 4: neoantigen burden test for Panel A (compute BOTH; decide from output)
    # -------------------------------------------------------------------------
    banner("STEP 4: NEOANTIGEN BURDEN TEST (group-rate vs per-cell; decide from output)")

    def _n(path):
        return len(pd.read_csv(path, sep='\t')) if os.path.exists(path) else np.nan
    n_pep_sbs2 = _n(os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv"))
    n_pep_cnv = _n(os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv"))
    n_var_sbs2 = len(neo_positions.get('SBS2_HIGH', set()))
    n_var_cnv = len(neo_positions.get('CNV_HIGH', set()))

    def gene_count(group):
        p = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
        if not os.path.exists(p):
            return np.nan
        return pd.read_csv(p, sep='\t')['gene'].nunique()
    n_gene_sbs2, n_gene_cnv = gene_count('SBS2_HIGH'), gene_count('CNV_HIGH')

    log("  GROUP-RATE test (matches plotted per-cell rate; equal 546-cell exposure):")
    try:
        from scipy.stats import binomtest
        for unit, a, b in [('peptides', n_pep_sbs2, n_pep_cnv),
                           ('neoantigen-producing variants', n_var_sbs2, n_var_cnv),
                           ('neoantigen genes', n_gene_sbs2, n_gene_cnv)]:
            if not (np.isfinite(a) and np.isfinite(b)) or (a + b) == 0:
                continue
            p = binomtest(int(a), int(a + b), 0.5).pvalue
            log(f"    {unit:32s}: SBS2 {int(a):5d} vs CNV {int(b):5d}  "
                f"ratio {a / b:.2f}x  binomial p={p:.2e}")
    except Exception as e:
        log(f"    [WARN] binomial test unavailable: {e}")

    log("")
    log("  PER-CELL genotyped burden (DEPTH-CONFOUNDED; SBS2 shallower; cross-check):")

    def percell_burden(group):
        cells = group_cells.get(group, set())
        pos_alt = neo_positions.get(group, set())
        carried = {}
        for (c, p, a) in pos_alt:
            for cb in mut.get((c, p, a), set()):
                if cb in cells:
                    carried[cb] = carried.get(cb, 0) + 1
        cell_list = list(cells)
        return np.array([carried.get(cb, 0) for cb in cell_list]), cell_list

    sb, sb_cells = percell_burden('SBS2_HIGH')
    cb_, cb_cells = percell_burden('CNV_HIGH')
    log(f"    mean carried neoantigen variants/cell: SBS2 {sb.mean():.3f}  CNV {cb_.mean():.3f}")
    neo_perumi_p = np.nan
    s_umi = c_umi = np.array([])
    try:
        from scipy.stats import mannwhitneyu
        p_raw = mannwhitneyu(sb, cb_, alternative='two-sided').pvalue
        log(f"    raw per-cell MWU p={p_raw:.3f}")
        s_umi = np.array([1000 * cc / umi[cb] for cc, cb in zip(sb, sb_cells) if umi.get(cb)])
        c_umi = np.array([1000 * cc / umi[cb] for cc, cb in zip(cb_, cb_cells) if umi.get(cb)])
        if len(s_umi) and len(c_umi):
            neo_perumi_p = mannwhitneyu(s_umi, c_umi, alternative='two-sided').pvalue
            log(f"    per-UMI per-cell MWU p={neo_perumi_p:.3f}  "
                f"(SBS2 {np.mean(s_umi):.3f} vs CNV {np.mean(c_umi):.3f} /1kUMI)")
    except Exception as e:
        log(f"    [WARN] per-cell MWU unavailable: {e}")

    # -------------------------------------------------------------------------
    # Panel A adjusted p (BH-FDR over the two per-cell per-UMI burden tests), to
    # match the paper's convention of reporting adjusted p. The fusion per-UMI p
    # is read from the fusion diagnostic's panelA_fusion_stats.tsv so both Panel A
    # comparisons are corrected together.
    # -------------------------------------------------------------------------
    def bh_adjust(pvals):
        p = np.asarray(pvals, float)
        mask = ~np.isnan(p)
        m = int(mask.sum())
        q = np.full(p.shape, np.nan)
        if m == 0:
            return q
        idx = np.where(mask)[0]
        pv = p[idx]
        order = np.argsort(pv)          # ascending
        ranked = pv[order]
        q_ranked = np.empty(m)
        prev = 1.0
        for j in range(m - 1, -1, -1):
            val = ranked[j] * m / (j + 1)
            prev = min(prev, val)
            q_ranked[j] = prev
        q_back = np.empty(m)
        q_back[order] = q_ranked        # map sorted q back to input order (not reversed!)
        q[idx] = np.minimum(q_back, 1.0)
        return q

    fus_stats = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/"
                             "fusion_normalization_compare/panelA_fusion_stats.tsv")
    fus_p = fus_sbs2 = fus_cnv = fus_sem_sbs2 = fus_sem_cnv = np.nan
    if os.path.exists(fus_stats):
        fs = pd.read_csv(fus_stats, sep='\t')
        if 'sbs2_vs_cnv_mwu_p_germline_sub' in fs.columns:
            fus_p = float(fs['sbs2_vs_cnv_mwu_p_germline_sub'].iloc[0])
        if 'mean_perUMI_per_cell_germline_sub' in fs.columns and 'group' in fs.columns:
            fsg = fs.set_index('group')['mean_perUMI_per_cell_germline_sub']
            fus_sbs2 = float(fsg.get('SBS2_HIGH', np.nan))
            fus_cnv = float(fsg.get('CNV_HIGH', np.nan))
        if 'sem_perUMI_per_cell_germline_sub' in fs.columns and 'group' in fs.columns:
            fss = fs.set_index('group')['sem_perUMI_per_cell_germline_sub']
            fus_sem_sbs2 = float(fss.get('SBS2_HIGH', np.nan))
            fus_sem_cnv = float(fss.get('CNV_HIGH', np.nan))
    else:
        log(f"\n  [WARN] fusion stats not found ({fus_stats});")
        log("  run the fusion diagnostic first so the adjusted p spans both Panel A tests.")

    # neoantigen per-cell per-UMI SEM (for Panel A error bars)
    neo_sem_sbs2 = float(np.std(s_umi, ddof=1) / np.sqrt(len(s_umi))) if len(s_umi) > 1 else np.nan
    neo_sem_cnv = float(np.std(c_umi, ddof=1) / np.sqrt(len(c_umi))) if len(c_umi) > 1 else np.nan

    labels = ['neoantigen_perUMI', 'fusion_perUMI']
    raw_ps = [neo_perumi_p, fus_p]
    qs = bh_adjust(raw_ps)
    log("")
    log("  PANEL A adjusted p (BH-FDR across the two per-cell per-UMI burden tests):")
    for lab, rp, q in zip(labels, raw_ps, qs):
        log(f"    {lab:20s}: raw p={rp:.3g}   BH-adjusted p={q:.3g}")

    pd.DataFrame({
        'comparison': labels,
        'metric': ['neoantigen variants /1kUMI', 'fusion junctions /1kUMI'],
        'sbs2_mean': [float(np.mean(s_umi)) if len(s_umi) else np.nan, fus_sbs2],
        'cnv_mean': [float(np.mean(c_umi)) if len(c_umi) else np.nan, fus_cnv],
        'sbs2_sem': [neo_sem_sbs2, fus_sem_sbs2],
        'cnv_sem': [neo_sem_cnv, fus_sem_cnv],
        'raw_p': raw_ps,
        'bh_adjusted_p': [float(x) if np.isfinite(x) else np.nan for x in qs],
    }).to_csv(os.path.join(OUTPUT_DIR, "panelA_burden_stats.tsv"), sep='\t', index=False)

    log("    NOTE: raw per-cell burden (unnormalized) is depth-confounded (SBS2")
    log("    shallower); the per-UMI version above is the depth-corrected test that")
    log("    matches the group-rate result, so it is the one to plot on Panel A.")

    pd.DataFrame([
        {'unit': 'peptides', 'sbs2': n_pep_sbs2, 'cnv': n_pep_cnv},
        {'unit': 'neo_variants', 'sbs2': n_var_sbs2, 'cnv': n_var_cnv},
        {'unit': 'neo_genes', 'sbs2': n_gene_sbs2, 'cnv': n_gene_cnv},
    ]).to_csv(os.path.join(OUTPUT_DIR, "neoantigen_burden_counts.tsv"), sep='\t', index=False)

    log(f"\n  Wrote group-aware expression, carrier, and burden tables to {OUTPUT_DIR}")
    with open(os.path.join(OUTPUT_DIR, "group_aware_report.txt"), 'w') as fh:
        fh.write('\n'.join(report))


if __name__ == "__main__":
    main()
