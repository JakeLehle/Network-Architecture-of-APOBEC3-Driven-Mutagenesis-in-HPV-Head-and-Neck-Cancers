#!/usr/bin/env python3
"""
Diagnostic_TCW_Neoantigen_Ranking.py
=====================================
Figure 7 diagnostic: find the "definitive A3-induced" neoantigens.

Goal
----
ANXA1 tops the Step04 expression-weighted composite but none of its mutations
are in TCW (TpCpW APOBEC) context, so it cannot be tied mechanistically to
A3-driven deamination. This diagnostic re-ranks the neoantigen landscape on the
axis Step04 ignores: whether the mutation that flips a wild-type non-binder into
a mutant MHC-I binder is itself an A3 deamination product (C>T / C>G at TpCpW).

It restricts to Tier 1 (shared) and Tier 2 (SBS2-specific) genes, keeps only
DIFFERENTIAL binders (mutant IC50 < 500 nM AND wild-type IC50 > 500 nM, so the
mutation creates the epitope), classifies each binding-driving variant's
trinucleotide context, and produces two ranked tables:

  1. Variant-level  : one row per differential neoantigen variant, filtered to
                      TCW context + Tier 1/2, sorted by MHC binding gain.
  2. Gene-level     : per-gene rollup with percent-TCW of the binding-driving
                      variants + best binder, i.e. the ANXA1-replacement list.

Read-only. Reads existing pipeline outputs + the reference genome FASTA only.
Nothing on disk is modified.

Key design notes
----------------
- Tiers are derived from the neoantigen-gene Venn (Tier 1 = shared across
  SBS2/CNV, Tier 2 = SBS2-only), which reproduces the manuscript's
  105 / 276 / 135 split without needing Step06's tier file. Derived counts are
  reconciled against those anchors and a warning is printed on mismatch.
- Trinucleotide context is fetched from the reference genome FASTA (not the
  SComatic REF_TRI column). The coordinate convention of the annotation `pos`
  is NOT assumed: the script calibrates it against the genome by testing which
  0-based offset (pos-1 / pos / pos-2 / pos+1) reproduces the ref allele, picks
  the offset that agrees genome-wide, and HARD-STOPS if none clears 95% (a build
  mismatch rather than an offset). This replaced a v1 hardcoded window that was
  off relative to this annotation's convention (225/332 ref mismatches, 0 TCW).
- TCW is classified on the GENOMIC ref/alt in pyrimidine orientation (G>A / G>C
  are reverse-complemented to C>T / C>G), matching Step01. Minus-strand genes
  (e.g. ISG15, ACOT7: genomic G>A, coding C>T) are handled correctly.
- Two substitution scopes are carried per variant:
    is_tcw     : C>T and C>G at TpCpW  (SBS2 + SBS13, current manuscript text)
    is_tcw_ct  : C>T at TpCpW only     (SBS2 arm alone)
  Ranking uses is_tcw by default; is_tcw_ct is kept as a column to tighten later.
- Anchor annotation: a differential mutation at P2 or the C-terminus (PΩ) is the
  canonical binding-driving position. Reported as a column, not a hard filter.
- Expression (to argue "broadly expressed like ANXA1") is merged from Step04's
  {group}_expression_weighted_ranking.tsv IF present, otherwise skipped.

Run in the NEOANTIGEN conda env (needs pysam):
    conda run -n NEOANTIGEN python Diagnostic_TCW_Neoantigen_Ranking.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import math
import numpy as np
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIGURATION  (confirm these paths match your tree)
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG7         = os.path.join(PROJECT_ROOT, "data/FIG_7")
MHC_DIR      = os.path.join(FIG7, "03_mhc_binding")
ANN_DIR      = os.path.join(FIG7, "02_snpeff_annotation")

# Reference genome used for alignment (source of trinucleotide context).
# Requires a genome.fa.fai index alongside it (already present if STAR/CellRanger
# aligned against this FASTA).
REF_GENOME   = "/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"

# Diagnostic outputs (TROUBLESHOOTING -> excluded from walkthroughs by convention)
OUTPUT_DIR   = os.path.join(FIG7, "TROUBLESHOOTING", "tcw_neoantigen_ranking")

# Focus group for the A3 story; CNV is loaded only to derive shared vs SBS2-only.
FOCUS_GROUP  = "SBS2_HIGH"
OTHER_GROUP  = "CNV_HIGH"

# Manuscript tier anchors for reconciliation
MANUSCRIPT_SHARED   = 105   # Tier 1
MANUSCRIPT_SBS2ONLY = 276   # Tier 2
MANUSCRIPT_CNVONLY  = 135   # Tier 3

# MHCflurry thresholds (match pipeline)
BINDER_THRESH = 500.0       # nM
WT_NONBINDER  = 500.0       # nM  (differential requires wt > this)
WT_CAP        = 50000.0     # nM  (cap wt IC50 for the log gain; also guards the
                            #      999999 sentinel from blowing up the metric)

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)


# =============================================================================
# TCW CLASSIFICATION  (genomic ref/alt, pyrimidine-oriented; matches Step01)
# =============================================================================
def classify_tcw(ref, alt, tri):
    """
    Classify a genomic substitution's APOBEC context from the reference
    trinucleotide (as read on the + strand).

    Returns dict:
      valid_tri   : trinucleotide is a clean 3-mer of ACGT
      sbs_class   : 'C>T (SBS2)' | 'C>G (SBS13)' | 'C>A' | 'not_cytosine'
      oriented    : e.g. 'T[C>T]A' after pyrimidine orientation
      is_tcw      : C>T or C>G at TpCpW  (SBS2 + SBS13 scope)
      is_tcw_ct   : C>T at TpCpW only     (SBS2 arm)
    """
    ref = str(ref).upper()
    alt = str(alt).upper()
    tri = str(tri).upper()

    res = {'valid_tri': False, 'sbs_class': 'not_cytosine',
           'oriented': '', 'is_tcw': False, 'is_tcw_ct': False}

    if len(tri) != 3 or any(b not in 'ACGT' for b in tri):
        return res
    res['valid_tri'] = True

    # Orient to the pyrimidine (C). APOBEC deaminates C on either strand.
    if ref == 'C':
        pref, palt, ptri = ref, alt, tri
    elif ref == 'G':
        pref = 'C'
        palt = COMPLEMENT.get(alt, 'N')
        ptri = ''.join(COMPLEMENT[b] for b in reversed(tri))
    else:
        return res  # ref is A or T -> not a cytosine substitution

    up, down = ptri[0], ptri[2]
    in_motif = (up == 'T' and down in ('A', 'T'))
    res['oriented'] = f"{up}[{pref}>{palt}]{down}"

    if palt == 'T':
        res['sbs_class'] = 'C>T (SBS2)'
        res['is_tcw']    = in_motif
        res['is_tcw_ct'] = in_motif
    elif palt == 'G':
        res['sbs_class'] = 'C>G (SBS13)'
        res['is_tcw']    = in_motif      # counts under SBS2+SBS13 scope
        res['is_tcw_ct'] = False         # not C>T
    elif palt == 'A':
        res['sbs_class'] = 'C>A'         # not an APOBEC signature
    return res


def binding_gain(mut_ic50, wt_ic50):
    """log10 fold improvement of mutant over wild-type, wt capped at WT_CAP."""
    try:
        m = float(mut_ic50)
        w = min(float(wt_ic50), WT_CAP)
        if m <= 0:
            return float('nan')
        return math.log10(w) - math.log10(m)
    except Exception:
        return float('nan')


# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("Diagnostic: TCW-driven neoantigen ranking (Tier 1 + Tier 2)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    log(f"  Output dir: {OUTPUT_DIR}")

    # -- pysam (fail early with a helpful message) ---------------------------
    try:
        import pysam
    except ImportError:
        log("  [FATAL] pysam not importable. Run in the NEOANTIGEN env:")
        log("          conda run -n NEOANTIGEN python Diagnostic_TCW_Neoantigen_Ranking.py")
        sys.exit(1)

    if not os.path.exists(REF_GENOME):
        log(f"  [FATAL] Reference genome not found: {REF_GENOME}")
        sys.exit(1)
    if not os.path.exists(REF_GENOME + ".fai"):
        log(f"  [WARNING] No .fai index at {REF_GENOME}.fai")
        log( "            pysam will try to build one; if the ref dir is read-only")
        log( "            run:  conda run -n NEOANTIGEN samtools faidx " + REF_GENOME)
    fasta = pysam.FastaFile(REF_GENOME)
    fasta_chroms = set(fasta.references)

    # -- 1. Load neoantigen files, derive tiers ------------------------------
    banner("STEP 1: Load neoantigens and derive tiers from the Venn", char="-")
    sbs2_path = os.path.join(MHC_DIR, f"{FOCUS_GROUP}_neoantigens.tsv")
    cnv_path  = os.path.join(MHC_DIR, f"{OTHER_GROUP}_neoantigens.tsv")
    for p in (sbs2_path, cnv_path):
        if not os.path.exists(p):
            log(f"  [FATAL] Missing neoantigen file: {p}")
            sys.exit(1)

    sbs2_neo = pd.read_csv(sbs2_path, sep='\t')
    cnv_neo  = pd.read_csv(cnv_path, sep='\t')
    log(f"  {FOCUS_GROUP}: {len(sbs2_neo)} binder peptides, "
        f"{sbs2_neo['gene'].nunique()} genes")
    log(f"  {OTHER_GROUP}: {len(cnv_neo)} binder peptides, "
        f"{cnv_neo['gene'].nunique()} genes")

    sbs2_genes = set(sbs2_neo['gene'].dropna())
    cnv_genes  = set(cnv_neo['gene'].dropna())
    shared     = sbs2_genes & cnv_genes
    sbs2_only  = sbs2_genes - cnv_genes
    cnv_only   = cnv_genes - sbs2_genes

    def gene_tier(g):
        if g in shared:
            return "Tier1_shared"
        if g in sbs2_only:
            return "Tier2_sbs2_specific"
        return None  # cnv-only or absent -> excluded

    log(f"\n  Derived Venn:")
    log(f"    Tier 1 (shared)        : {len(shared):4d}   (manuscript {MANUSCRIPT_SHARED})")
    log(f"    Tier 2 (SBS2-specific) : {len(sbs2_only):4d}   (manuscript {MANUSCRIPT_SBS2ONLY})")
    log(f"    Tier 3 (CNV-specific)  : {len(cnv_only):4d}   (manuscript {MANUSCRIPT_CNVONLY})")
    for label, derived, anchor in [("Tier1", len(shared), MANUSCRIPT_SHARED),
                                    ("Tier2", len(sbs2_only), MANUSCRIPT_SBS2ONLY),
                                    ("Tier3", len(cnv_only), MANUSCRIPT_CNVONLY)]:
        if derived != anchor:
            log(f"  [WARNING] {label} derived {derived} != manuscript {anchor} "
                f"(diff {derived - anchor}). Neoantigen files may differ from the "
                f"submitted run; check before trusting downstream ranks.")

    # -- 2. Load annotation for genomic ref/alt ------------------------------
    banner("STEP 2: Load SnpEff annotation for genomic ref/alt", char="-")
    ann_path = os.path.join(ANN_DIR, f"{FOCUS_GROUP}.somatic_protein_altering.tsv")
    if not os.path.exists(ann_path):
        log(f"  [FATAL] Missing annotation file: {ann_path}")
        sys.exit(1)
    ann = pd.read_csv(ann_path, sep='\t')
    ann['location'] = ann['chrom'].astype(str) + ':' + ann['pos'].astype(str)
    # (location, hgvs_p) -> (chrom, pos, ref, alt); dedup keeps first
    ann_key = (ann.drop_duplicates(subset=['location', 'hgvs_p'])
                  .set_index(['location', 'hgvs_p'])[['chrom', 'pos', 'ref', 'alt']]
                  .to_dict('index'))
    log(f"  Loaded {len(ann)} protein-altering variants, "
        f"{len(ann_key)} unique (location, hgvs_p) keys")

    # -- 3. Differential peptides in the focus group -------------------------
    banner("STEP 3: Differential binders + context classification", char="-")
    diff = sbs2_neo[sbs2_neo['is_differential'] == True].copy()
    log(f"  Differential binder peptides ({FOCUS_GROUP}): {len(diff)}")
    log(f"  Distinct differential variants: "
        f"{diff.groupby(['location', 'hgvs_p']).ngroups}")

    diff['tier'] = diff['gene'].map(gene_tier)
    diff = diff[diff['tier'].notna()].copy()
    log(f"  Differential peptides in Tier 1/2 genes: {len(diff)} "
        f"({diff['gene'].nunique()} genes)")

    # Fetch trinucleotide + classify per UNIQUE variant.
    # The annotation 'pos' coordinate convention is auto-calibrated against the
    # genome rather than assumed: for each variant we test which 0-based genome
    # coordinate holds the ref allele, pick the offset that agrees genome-wide,
    # and require >=95% agreement before trusting any context call. This avoids
    # the silent off-by-one that produced 225/332 ref mismatches in v1.
    var_keys = diff[['location', 'hgvs_p']].drop_duplicates()
    n_no_annot  = 0
    n_bad_chrom = 0

    # candidate 0-based variant coordinate as a function of annotation pos
    CANDS = {
        'pos-1 (pos is 1-based VCF)': lambda p: p - 1,
        'pos   (pos is 0-based)':     lambda p: p,
        'pos-2 (double-shifted)':     lambda p: p - 2,
        'pos+1':                      lambda p: p + 1,
    }

    # Phase A: fetch a 6bp window per variant, tally which candidate coord == ref
    win = {}                                   # key -> (chrom,pos,ref,alt,seq,start)
    cand_hits = {name: 0 for name in CANDS}
    n_eval = 0
    for _, r in var_keys.iterrows():
        key = (r['location'], r['hgvs_p'])
        a = ann_key.get(key)
        if a is None:
            n_no_annot += 1
            win[key] = None
            continue
        chrom, pos = str(a['chrom']), int(a['pos'])
        ref, alt = str(a['ref']).upper(), str(a['alt']).upper()
        if chrom not in fasta_chroms:
            n_bad_chrom += 1
            win[key] = None
            continue
        s = max(0, pos - 3)                     # window covers 0-based pos-3..pos+2
        try:
            w = fasta.fetch(chrom, s, pos + 3).upper()
        except Exception:
            w = ''
        win[key] = (chrom, pos, ref, alt, w, s)
        if len(ref) == 1:                       # calibrate on SNVs only
            n_eval += 1
            for name, fn in CANDS.items():
                idx = fn(pos) - s
                if 0 <= idx < len(w) and w[idx] == ref:
                    cand_hits[name] += 1

    # Phase B: choose the offset that reproduces the ref allele genome-wide
    banner("Coordinate calibration (ref allele vs genome)", char=".")
    best_name, best_hits = None, -1
    for name, hits in cand_hits.items():
        frac = 100 * hits / n_eval if n_eval else 0
        log(f"    {name:28s}: {hits}/{n_eval} ref matches ({frac:5.1f}%)")
        if hits > best_hits:
            best_name, best_hits = name, hits
    best_frac = best_hits / n_eval if n_eval else 0
    log(f"  Chosen convention: {best_name}  ({100*best_frac:.1f}% ref match)")
    if best_frac < 0.95:
        log("  [FATAL] No coordinate offset reproduces the ref allele genome-wide.")
        log("          This is a genome-build / annotation mismatch, not an offset.")
        log("          Stopping before any context call is trusted.")
        sys.exit(1)
    coord_fn = CANDS[best_name]

    # Phase C: classify each variant's trinucleotide at the calibrated coordinate
    tri_cache = {}
    n_ref_mismatch = 0
    for key, v in win.items():
        if v is None:
            tri_cache[key] = None
            continue
        chrom, pos, ref, alt, w, s = v
        idx = coord_fn(pos) - s                 # window index of the variant base
        tri = w[idx - 1: idx + 2] if (idx - 1 >= 0 and idx + 2 <= len(w)) else ''
        if len(tri) == 3 and tri[1] != ref:
            n_ref_mismatch += 1
        c = classify_tcw(ref, alt, tri)
        c.update({'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                  'tri': tri, 'ref_matches': (len(tri) == 3 and tri[1] == ref)})
        tri_cache[key] = c

    log(f"\n  Unique variants classified: {sum(v is not None for v in tri_cache.values())}")
    log(f"    missing from annotation : {n_no_annot}")
    log(f"    chrom not in FASTA      : {n_bad_chrom}")
    log(f"    residual ref mismatches : {n_ref_mismatch} "
        f"(at calibrated offset; nonzero -> indels/edge cases, inspect individually)")

    # -- 4. Peptide-level enrichment -----------------------------------------
    def anchor_flag(row):
        L = int(row['peptide_length'])
        p = int(row['mut_position_in_peptide'])   # 0-indexed
        return (p == 1) or (p == L - 1)            # P2 or C-terminus (PΩ)

    recs = []
    for _, row in diff.iterrows():
        key = (row['location'], row['hgvs_p'])
        c = tri_cache.get(key)
        if c is None:
            continue
        recs.append({
            'gene': row['gene'], 'tier': row['tier'],
            'location': row['location'], 'hgvs_p': row['hgvs_p'],
            'chrom': c['chrom'], 'pos': c['pos'],
            'genomic_sub': f"{c['ref']}>{c['alt']}",
            'oriented': c['oriented'], 'sbs_class': c['sbs_class'],
            'tri': c['tri'], 'ref_matches': c['ref_matches'],
            'is_tcw': c['is_tcw'], 'is_tcw_ct': c['is_tcw_ct'],
            'mut_peptide': row['mut_peptide'], 'wt_peptide': row['wt_peptide'],
            'best_allele': row['best_allele'],
            'peptide_length': int(row['peptide_length']),
            'mut_pos_in_pep': int(row['mut_position_in_peptide']),
            'is_anchor': anchor_flag(row),
            'mut_ic50': float(row['mut_ic50']),
            'wt_ic50': float(row['wt_ic50']),
            'binding_gain': binding_gain(row['mut_ic50'], row['wt_ic50']),
        })
    pep = pd.DataFrame(recs)
    if pep.empty:
        log("  [FATAL] No differential peptides survived classification.")
        sys.exit(1)

    # -- 5. Variant-level table (best differential peptide per variant) ------
    banner("STEP 4: Variant-level ranked list", char="-")
    grp = pep.groupby(['gene', 'tier', 'location', 'hgvs_p', 'genomic_sub',
                       'oriented', 'sbs_class', 'tri', 'ref_matches',
                       'is_tcw', 'is_tcw_ct'], dropna=False)
    var_rows = []
    for gkey, sub in grp:
        best = sub.loc[sub['binding_gain'].idxmax()]
        var_rows.append({
            'gene': gkey[0], 'tier': gkey[1], 'location': gkey[2],
            'hgvs_p': gkey[3], 'genomic_sub': gkey[4], 'oriented': gkey[5],
            'sbs_class': gkey[6], 'tri': gkey[7], 'ref_matches': gkey[8],
            'is_tcw': gkey[9], 'is_tcw_ct': gkey[10],
            'n_diff_peptides': len(sub),
            'best_gain': round(float(best['binding_gain']), 3),
            'best_mut_ic50': round(float(sub['mut_ic50'].min()), 2),
            'gain_mut_ic50': round(float(best['mut_ic50']), 2),
            'gain_wt_ic50': round(float(best['wt_ic50']), 2),
            'best_allele': best['best_allele'],
            'best_peptide': best['mut_peptide'],
            'best_wt_peptide': best['wt_peptide'],
            'mut_pos_in_pep': int(best['mut_pos_in_pep']),
            'peptide_length': int(best['peptide_length']),
            'is_anchor': bool(sub['is_anchor'].any()),
        })
    variants = pd.DataFrame(var_rows)

    # Full differential variant table (all, tcw flag included) for context
    variants_all = variants.sort_values(
        ['is_tcw', 'best_gain'], ascending=[False, False]).reset_index(drop=True)
    variants_all.to_csv(
        os.path.join(OUTPUT_DIR, "differential_variants_tier12_ALL.tsv"),
        sep='\t', index=False)

    # The headline list: TCW-context differential neoantigens, ranked by gain
    tcw_ranked = (variants[variants['is_tcw']]
                  .sort_values(['best_gain', 'best_mut_ic50'],
                               ascending=[False, True])
                  .reset_index(drop=True))
    tcw_ranked.insert(0, 'rank', tcw_ranked.index + 1)
    tcw_ranked.to_csv(
        os.path.join(OUTPUT_DIR, "tcw_differential_variants_RANKED.tsv"),
        sep='\t', index=False)

    n_diff_var = len(variants)
    n_tcw_var  = int(variants['is_tcw'].sum())
    log(f"  Differential variants (Tier 1/2)      : {n_diff_var}")
    log(f"  ... in TCW context (SBS2+SBS13)        : {n_tcw_var} "
        f"({100*n_tcw_var/n_diff_var:.1f}%)")
    log(f"  ... in TCW context, C>T only (SBS2)    : {int(variants['is_tcw_ct'].sum())}")
    log(f"\n  TOP 25 TCW-driven differential neoantigens (by binding gain):")
    log(f"  {'rk':>3} {'gene':12} {'tier':9} {'hgvs_p':16} {'oriented':11} "
        f"{'gain':>5} {'mutIC50':>8} {'wtIC50':>9} {'anc':>3} {'allele':10}")
    for _, r in tcw_ranked.head(25).iterrows():
        log(f"  {r['rank']:>3} {str(r['gene'])[:12]:12} "
            f"{r['tier'].replace('Tier','T').replace('_shared','_sh').replace('_sbs2_specific','_s2'):9} "
            f"{str(r['hgvs_p'])[:16]:16} {r['oriented']:11} "
            f"{r['best_gain']:5.2f} {r['gain_mut_ic50']:8.1f} {r['gain_wt_ic50']:9.1f} "
            f"{'Y' if r['is_anchor'] else '.':>3} {str(r['best_allele'])[:10]:10}")

    # -- 6. Gene-level rollup (the ANXA1-replacement candidates) --------------
    banner("STEP 5: Gene-level rollup (percent-TCW of binding-driving variants)",
           char="-")
    gene_rows = []
    for (gene, tier), sub in variants.groupby(['gene', 'tier']):
        n_var  = len(sub)
        tcw    = sub[sub['is_tcw']]
        n_tcw  = len(tcw)
        gene_rows.append({
            'gene': gene, 'tier': tier,
            'n_diff_variants': n_var,
            'n_tcw_variants': n_tcw,
            'pct_tcw': round(100 * n_tcw / n_var, 1) if n_var else 0.0,
            'n_tcw_ct_variants': int(sub['is_tcw_ct'].sum()),
            'best_tcw_gain': round(float(tcw['best_gain'].max()), 3) if n_tcw else np.nan,
            'best_tcw_mut_ic50': round(float(tcw['best_mut_ic50'].min()), 2) if n_tcw else np.nan,
            'any_tcw_anchor': bool(tcw['is_anchor'].any()) if n_tcw else False,
        })
    genes = pd.DataFrame(gene_rows)

    # Optional expression merge (Step04 ranking, if it exists)
    expr_path = os.path.join(MHC_DIR, f"{FOCUS_GROUP}_expression_weighted_ranking.tsv")
    have_expr = os.path.exists(expr_path)
    if have_expr:
        er = pd.read_csv(expr_path, sep='\t')
        keep = {c: c for c in ['gene', 'mean_expr_group', 'pct_expressing_group',
                               'composite_score'] if c in er.columns}
        genes = genes.merge(er[list(keep)], on='gene', how='left')
        log(f"  Merged expression from {os.path.basename(expr_path)}")
    else:
        log(f"  [INFO] {os.path.basename(expr_path)} not found; expression columns "
            f"omitted. Run Step04, or point me at adata and I'll add the pull.")

    # Interpretive lead score (NOT a statistic): reward high TCW fraction AND a
    # strong TCW binder. Expression-weighted variant added when available.
    def lead_score(row):
        if not row['n_tcw_variants'] or pd.isna(row['best_tcw_mut_ic50']) \
                or row['best_tcw_mut_ic50'] <= 0:
            return 0.0
        s = (row['pct_tcw'] / 100.0) * (BINDER_THRESH / row['best_tcw_mut_ic50'])
        if have_expr and 'mean_expr_group' in row and pd.notna(row['mean_expr_group']):
            s *= float(row['mean_expr_group'])
        return round(s, 4)
    genes['lead_score'] = genes.apply(lead_score, axis=1)

    genes_all = genes.sort_values(
        ['lead_score', 'pct_tcw'], ascending=[False, False]).reset_index(drop=True)
    genes_all.to_csv(
        os.path.join(OUTPUT_DIR, "gene_level_tcw_summary_ALL.tsv"),
        sep='\t', index=False)

    candidates = (genes[genes['n_tcw_variants'] >= 1]
                  .sort_values(['lead_score', 'pct_tcw', 'best_tcw_gain'],
                               ascending=[False, False, False])
                  .reset_index(drop=True))
    candidates.insert(0, 'rank', candidates.index + 1)
    candidates.to_csv(
        os.path.join(OUTPUT_DIR, "gene_level_A3_neoantigen_CANDIDATES.tsv"),
        sep='\t', index=False)

    expr_hdr = f"{'meanExpr':>8} {'%expr':>6} " if have_expr else ""
    log(f"\n  TOP 25 candidate A3-induced neoantigen GENES "
        f"(>=1 TCW binding-driving variant):")
    log(f"  {'rk':>3} {'gene':12} {'tier':9} {'nTCW/nDiff':>10} {'%TCW':>5} "
        f"{'bestGain':>8} {'bestIC50':>8} {'anc':>3} {expr_hdr}{'lead':>7}")
    for _, r in candidates.head(25).iterrows():
        extra = ""
        if have_expr:
            me = r.get('mean_expr_group', np.nan)
            pe = r.get('pct_expressing_group', np.nan)
            extra = f"{(me if pd.notna(me) else 0):8.2f} {(pe if pd.notna(pe) else 0):6.1f} "
        log(f"  {r['rank']:>3} {str(r['gene'])[:12]:12} "
            f"{r['tier'].replace('Tier','T').replace('_shared','_sh').replace('_sbs2_specific','_s2'):9} "
            f"{str(r['n_tcw_variants'])+'/'+str(r['n_diff_variants']):>10} "
            f"{r['pct_tcw']:5.1f} {r['best_tcw_gain']:8.2f} {r['best_tcw_mut_ic50']:8.1f} "
            f"{'Y' if r['any_tcw_anchor'] else '.':>3} {extra}{r['lead_score']:7.3f}")

    # ANXA1 sanity check (should show 0 TCW variants)
    anx = genes_all[genes_all['gene'] == 'ANXA1']
    if len(anx):
        a = anx.iloc[0]
        log(f"\n  ANXA1 control: tier={a['tier']}, n_diff_variants={a['n_diff_variants']}, "
            f"n_tcw={a['n_tcw_variants']}, pct_tcw={a['pct_tcw']} "
            f"(expected 0 TCW -> confirms it is not an A3 substrate)")

    # -- 7. Write report -----------------------------------------------------
    banner("OUTPUTS", char="-")
    for f in ["tcw_differential_variants_RANKED.tsv",
              "differential_variants_tier12_ALL.tsv",
              "gene_level_A3_neoantigen_CANDIDATES.tsv",
              "gene_level_tcw_summary_ALL.tsv"]:
        log(f"  {os.path.join(OUTPUT_DIR, f)}")
    report_path = os.path.join(OUTPUT_DIR, "tcw_neoantigen_ranking_report.txt")
    with open(report_path, 'w') as fh:
        fh.write('\n'.join(report_lines))
    log(f"  {report_path}")
    log("\n  Done.")


if __name__ == "__main__":
    main()
