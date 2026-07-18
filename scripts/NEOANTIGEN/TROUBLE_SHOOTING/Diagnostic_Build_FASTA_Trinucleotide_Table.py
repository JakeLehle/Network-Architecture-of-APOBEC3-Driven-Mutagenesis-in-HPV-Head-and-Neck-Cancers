#!/usr/bin/env python3
"""
Diagnostic_Build_FASTA_Trinucleotide_Table.py
=============================================
Genome-correct trinucleotide context + TCW classification for EVERY
protein-altering variant, read from the reference FASTA (never SComatic REF_TRI).

Why
---
SComatic's REF_TRI column matches the genome only 23/355 (FASTA arbiter): its
centre base is right but the flanks are wrong. The ranking diagnostic already
proved that a calibrated FASTA read (ref base at 0-based pos-1) reproduces the
genome 355/355 on the differential set. This script generalises that read to all
protein-altering variants and writes one shared table so that:
  - the Figure 7 gene track can colour every lollipop by real TCW status, and
  - BEAT 3 / BEAT 5 / Step01's per-group APOBEC summary can be re-pointed off the
    broken REF_TRI onto a genome-correct source (Section 7 / F32 fix).

Coordinates
-----------
Annotation `pos` is 1-based; the reference base sits at 0-based pos-1 (arbiter-
validated). pysam fetch is 0-based half-open, so the trinucleotide centred on the
reference base is fetch(chrom, pos-2, pos+1) -> [5' neighbour, ref, 3' neighbour].
We assert the fetched centre equals REF for every variant; anything else means a
coordinate/contig mismatch and is flagged, never silently classified.

TCW classification (pyrimidine / SBS orientation)
-------------------------------------------------
APOBEC edits a C in a TpCpW context (5'=T, centre C, 3'=W in {A,T}). For a G-ref
call we reverse-complement the genome tri to reveal the C-centred context. Then:
  is_tcw    = TCW context AND substitution is C>T or C>G   (SBS2 + SBS13)
  is_tcw_ct = TCW context AND substitution is C>T          (SBS2 arm)

READ-ONLY inputs; single output table. Env: NEOANTIGEN (needs pysam + genome.fa).
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
ANNOTATION_DIR = os.path.join(BASE_DIR, "data/FIG_7/02_snpeff_annotation")
MHC_DIR = os.path.join(BASE_DIR, "data/FIG_7/03_mhc_binding")
GENOME_FA = "/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"
OUTPUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/fasta_context")
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "ref_tri_fasta.tsv")

GROUPS = ['SBS2_HIGH', 'CNV_HIGH']
FIGURE_GENES = ['KLF3', 'CAST', 'SERPINB2', 'KRT6B', 'ANXA1']
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

report = []
def log(msg=""):
    print(msg, flush=True)
    report.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)


# =============================================================================
# TCW classification
# =============================================================================
def classify_tcw(tri_genome, ref, alt):
    """Return (tri_pyr, sub_pyr, is_tcw, is_tcw_ct).
    tri_genome is the plus-strand trinucleotide centred on the reference base."""
    if tri_genome is None or len(tri_genome) != 3:
        return None, None, False, False
    ref, alt = ref.upper(), alt.upper()
    if ref == 'C':
        tri_p, alt_p = tri_genome, alt
    elif ref == 'G':
        tri_p = ''.join(COMP.get(b, 'N') for b in reversed(tri_genome))
        alt_p = COMP.get(alt, 'N')
    else:
        # A/T-ref call: not a C-centred (APOBEC-eligible) event
        return tri_genome, None, False, False
    context_tcw = (tri_p[0] == 'T' and tri_p[1] == 'C' and tri_p[2] in ('A', 'T'))
    sub_pyr = f"C>{alt_p}"
    is_tcw = context_tcw and alt_p in ('T', 'G')
    is_tcw_ct = context_tcw and alt_p == 'T'
    return tri_p, sub_pyr, is_tcw, is_tcw_ct


# =============================================================================
# LOAD variants (union of protein-altering across both groups)
# =============================================================================
def load_protein_altering():
    frames = []
    for group in GROUPS:
        p = os.path.join(ANNOTATION_DIR, f"{group}.somatic_protein_altering.tsv")
        if not os.path.exists(p):
            log(f"  [WARN] not found: {p}")
            continue
        df = pd.read_csv(p, sep='\t')
        lm = {c.lower(): c for c in df.columns}
        cc, pc, rc, ac = lm.get('chrom', lm.get('#chrom')), lm.get('pos'), lm.get('ref'), lm.get('alt')
        gc, hc = lm.get('gene'), lm.get('hgvs_p')
        if not all([cc, pc, rc, ac]):
            log(f"  [WARN] {group} missing chrom/pos/ref/alt (have {list(df.columns)})")
            continue
        ren = {cc: 'chrom', pc: 'pos', rc: 'ref', ac: 'alt'}
        if gc: ren[gc] = 'gene'
        if hc: ren[hc] = 'hgvs_p'
        df = df.rename(columns=ren)
        if 'gene' not in df.columns: df['gene'] = ''
        if 'hgvs_p' not in df.columns: df['hgvs_p'] = ''
        df['group'] = group
        frames.append(df[['chrom', 'pos', 'ref', 'alt', 'gene', 'hgvs_p', 'group']])
    if not frames:
        return pd.DataFrame()
    allv = pd.concat(frames, ignore_index=True)
    # neoantigen-producing flag by (gene, hgvs_p) present in the neoantigen calls
    neo_keys = set()
    for group in GROUPS:
        np_ = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
        if os.path.exists(np_):
            neo = pd.read_csv(np_, sep='\t')
            if {'gene', 'hgvs_p'} <= set(neo.columns):
                neo_keys |= set(zip(neo['gene'].astype(str), neo['hgvs_p'].astype(str)))
    allv['is_neo'] = [(g, str(h)) in neo_keys for g, h in zip(allv['gene'].astype(str), allv['hgvs_p'])]

    # collapse to unique variants, recording group membership
    grp = allv.groupby(['chrom', 'pos', 'ref', 'alt'], as_index=False).agg(
        gene=('gene', 'first'), hgvs_p=('hgvs_p', 'first'),
        is_neo=('is_neo', 'max'),
        in_sbs2=('group', lambda s: 'SBS2_HIGH' in set(s)),
        in_cnv=('group', lambda s: 'CNV_HIGH' in set(s)))
    return grp


# =============================================================================
# CONTIG resolution against the FASTA
# =============================================================================
def resolve_contig(chrom, refs):
    if chrom in refs:
        return chrom
    alt = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
    if alt in refs:
        return alt
    return None


def main():
    banner("BUILD SHARED FASTA TRINUCLEOTIDE TABLE")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    try:
        import pysam
    except Exception as e:
        log(f"  [FATAL] pysam not available ({e}); run in the NEOANTIGEN env.")
        return
    if not os.path.isfile(GENOME_FA):
        log(f"  [FATAL] genome FASTA not found: {GENOME_FA}")
        return

    v = load_protein_altering()
    if len(v) == 0:
        log("  [FATAL] no protein-altering variants loaded.")
        return
    log(f"  unique protein-altering variants: {len(v):,} "
        f"(SBS2 {int(v['in_sbs2'].sum())}, CNV {int(v['in_cnv'].sum())}, "
        f"neoantigen-producing {int(v['is_neo'].sum())})")

    fa = pysam.FastaFile(GENOME_FA)
    refs = set(fa.references)

    tri_g, tri_p, sub_p = [], [], []
    is_tcw, is_tcw_ct, refok, contig_used = [], [], [], []
    n_bad_contig = 0
    for r in v.itertuples():
        chrom = resolve_contig(str(r.chrom), refs)
        if chrom is None:
            n_bad_contig += 1
            tri_g.append(None); tri_p.append(None); sub_p.append(None)
            is_tcw.append(False); is_tcw_ct.append(False); refok.append(False)
            contig_used.append(None)
            continue
        pos = int(r.pos)
        try:
            tg = fa.fetch(chrom, pos - 2, pos + 1).upper()
        except Exception:
            tg = None
        centre = tg[1] if (tg and len(tg) == 3) else None
        ok = (centre == str(r.ref).upper())
        tp, sp, tcw, tcwct = classify_tcw(tg if ok else None, str(r.ref), str(r.alt))
        tri_g.append(tg); tri_p.append(tp); sub_p.append(sp)
        is_tcw.append(bool(tcw)); is_tcw_ct.append(bool(tcwct)); refok.append(bool(ok))
        contig_used.append(chrom)

    v = v.assign(contig=contig_used, tri_genome=tri_g, tri_pyr=tri_p, sub_pyr=sub_p,
                 is_tcw=is_tcw, is_tcw_ct=is_tcw_ct, ref_matches_genome=refok)

    def _category(ref, ok):
        if not ok:
            return 'genome_mismatch'
        return 'cytosine_ref' if str(ref).upper() in ('C', 'G') else 'non_cytosine_ref'
    v['category'] = [_category(r.ref, r.ref_matches_genome) for r in v.itertuples()]

    # -------------------------------------------------------------------------
    # Validation: every variant lands in exactly one category. The ONLY category
    # that is a problem is genome_mismatch; non_cytosine_ref (A/T reference) is a
    # correct, expected non-APOBEC classification, not a lookup failure.
    # -------------------------------------------------------------------------
    banner("VALIDATION", "-")
    n = len(v)
    n_ok = int(v['ref_matches_genome'].sum())
    log(f"  reference base matches genome at pos-1: {n_ok}/{n} ({100 * n_ok / n:.1f}%)")
    log("  category breakdown (mutually exclusive; every variant is accounted for):")
    for cat, cnt in v['category'].value_counts().items():
        log(f"    {cat:18s}: {cnt}")
    n_mismatch = int((v['category'] == 'genome_mismatch').sum())
    if n_bad_contig:
        log(f"  [WARN] {n_bad_contig} variants had an unresolved contig name")
    if n_mismatch:
        log(f"  [WARN][INVESTIGATE] {n_mismatch} genome mismatches:")
        for b in v[v['category'] == 'genome_mismatch'].head(10).itertuples():
            log(f"    {b.chrom}:{b.pos} ref={b.ref} genome_tri={b.tri_genome}")
    else:
        log("  -> 0 genome mismatches: every reference base is confirmed against GRCh38.")
    n_at = int((v['category'] == 'non_cytosine_ref').sum())
    log(f"  non_cytosine_ref (A/T reference; not APOBEC-eligible by definition -> tcw=False): {n_at}")
    log(f"  TCW (C>T + C>G, SBS2+SBS13): {int(v['is_tcw'].sum())}")
    log(f"  TCW clean C>T (SBS2 arm):    {int(v['is_tcw_ct'].sum())}")

    # -------------------------------------------------------------------------
    # Figure-gene sanity: label each variant by category (no false flags), then
    # confirm the handoff anchor context is PRESENT among the gene's clean-C>T set
    # (a gene may carry several TCW sites with different W bases).
    # -------------------------------------------------------------------------
    banner("FIGURE-GENE SANITY CHECK (categorised; no false flags)", "-")
    expected = {'KLF3': 'T[C>T]A', 'CAST': 'T[C>T]T'}   # main-figure anchor contexts
    for gene in FIGURE_GENES:
        sub = v[v['gene'] == gene]
        clean_ct_ctx = set()
        for b in sub.itertuples():
            if b.category == 'genome_mismatch':
                label = "GENOME MISMATCH -> investigate"
            elif b.category == 'non_cytosine_ref':
                label = f"{b.ref}-ref (non-cytosine); not APOBEC-eligible -> tcw=False (correct)"
            else:
                ctx = f"{b.tri_pyr[0]}[{b.sub_pyr}]{b.tri_pyr[2]}"
                label = f"{ctx}  tcw={b.is_tcw} ct={b.is_tcw_ct}"
                if b.is_tcw_ct:
                    clean_ct_ctx.add(ctx)
            log(f"  {gene:9s} {b.chrom}:{b.pos} {b.ref}>{b.alt}  {str(b.hgvs_p):14s}  {label}")
        if gene in expected:
            hit = expected[gene] in clean_ct_ctx
            log(f"    -> handoff anchor {expected[gene]} present among {gene} clean-C>T "
                f"variants: {'YES' if hit else 'NO -> investigate'}   "
                f"(clean-C>T set: {sorted(clean_ct_ctx)})")

    # -------------------------------------------------------------------------
    # Independent cross-validation against the arbiter-validated ranking file
    # (differential_variants_tier12_ALL.tsv). The FASTA arbiter confirmed that
    # file's `tri` matches the genome 355/355, so it is a known-good reference.
    # is_tcw / is_tcw_ct are orientation-independent, so they are the cleanest
    # thing to compare. 100% agreement on the overlap = high confidence.
    # -------------------------------------------------------------------------
    banner("CROSS-VALIDATION vs arbiter-validated ranking (differential_variants_tier12_ALL)")
    rank_path = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/"
                             "tcw_neoantigen_ranking/differential_variants_tier12_ALL.tsv")
    if not os.path.exists(rank_path):
        log(f"  [WARN] ranking file not found ({rank_path}); cross-validation skipped.")
    else:
        rk = pd.read_csv(rank_path, sep='\t')
        lm = {c.lower(): c for c in rk.columns}
        ren = {}
        for want, opts in [('r_chrom', ['chrom', '#chrom']), ('r_pos', ['pos', 'start']),
                           ('r_ref', ['ref']), ('r_alt', ['alt']), ('r_tri', ['tri']),
                           ('r_is_tcw', ['is_tcw']), ('r_is_tcw_ct', ['is_tcw_ct']),
                           ('r_gene', ['gene']), ('r_hgvs', ['hgvs_p'])]:
            for o in opts:
                if o in lm:
                    ren[lm[o]] = want
                    break
        rk = rk.rename(columns=ren)

        def _b(x):
            return str(x).strip().lower() in ('true', '1', '1.0', 'yes')

        have_bp = all(c in rk.columns for c in ['r_chrom', 'r_pos', 'r_ref', 'r_alt'])
        vv = v.copy()
        if have_bp:
            vv['_k'] = (vv['chrom'].astype(str) + ':' + vv['pos'].astype(int).astype(str) + ':' +
                        vv['ref'].astype(str).str.upper() + '>' + vv['alt'].astype(str).str.upper())
            rk['_k'] = (rk['r_chrom'].astype(str) + ':' + rk['r_pos'].astype(int).astype(str) + ':' +
                        rk['r_ref'].astype(str).str.upper() + '>' + rk['r_alt'].astype(str).str.upper())
            keyname = 'chrom:pos:ref:alt'
        elif 'r_gene' in rk.columns and 'r_hgvs' in rk.columns:
            vv['_k'] = vv['gene'].astype(str) + '|' + vv['hgvs_p'].astype(str)
            rk['_k'] = rk['r_gene'].astype(str) + '|' + rk['r_hgvs'].astype(str)
            keyname = 'gene|hgvs_p'
        else:
            log(f"  [WARN] ranking file lacks usable join keys (cols {list(rk.columns)}); skipped.")
            keyname = None

        if keyname:
            m = vv.merge(rk, on='_k', how='inner', suffixes=('', '_rk'))
            nc = len(m)
            log(f"  join key: {keyname}   overlap: {nc} variants "
                f"(ranking {len(rk)}, FASTA {len(vv)})")
            if nc:
                a_tcw = a_ct = a_tri = 0
                mism = []
                for r in m.itertuples():
                    r_tcw = _b(getattr(r, 'r_is_tcw', None)) if 'r_is_tcw' in m.columns else None
                    r_ct = _b(getattr(r, 'r_is_tcw_ct', None)) if 'r_is_tcw_ct' in m.columns else None
                    ok_tcw = (r_tcw is None) or (bool(r.is_tcw) == r_tcw)
                    ok_ct = (r_ct is None) or (bool(r.is_tcw_ct) == r_ct)
                    if r_tcw is not None and ok_tcw:
                        a_tcw += 1
                    if r_ct is not None and ok_ct:
                        a_ct += 1
                    if 'r_tri' in m.columns and str(getattr(r, 'r_tri', '')).upper() == str(r.tri_genome).upper():
                        a_tri += 1
                    if not (ok_tcw and ok_ct):
                        mism.append((r._k, r_tcw, r_ct, bool(r.is_tcw), bool(r.is_tcw_ct),
                                     str(getattr(r, 'r_tri', '')), str(r.tri_genome)))
                if 'r_is_tcw' in m.columns:
                    log(f"  is_tcw agreement:    {a_tcw}/{nc} ({100 * a_tcw / nc:.1f}%)")
                if 'r_is_tcw_ct' in m.columns:
                    log(f"  is_tcw_ct agreement: {a_ct}/{nc} ({100 * a_ct / nc:.1f}%)")
                if 'r_tri' in m.columns:
                    log(f"  tri (genome) agreement: {a_tri}/{nc} ({100 * a_tri / nc:.1f}%)")
                if mism:
                    log(f"  [INVESTIGATE] {len(mism)} classification mismatches:")
                    for x in mism[:15]:
                        log(f"    {x[0]}  ranking(tcw={x[1]},ct={x[2]},tri={x[5]}) "
                            f"vs FASTA(tcw={x[3]},ct={x[4]},tri={x[6]})")
                else:
                    log("  -> 100% agreement with the arbiter-validated ranking. High confidence.")
            else:
                log("  [WARN] no overlapping variants (check the join key columns).")


    v.to_csv(OUTPUT_TSV, sep='\t', index=False)
    log(f"\n  Wrote {len(v):,} rows -> {OUTPUT_TSV}")
    log("  Consumers: gene-track TCW colouring (join on chrom/pos/alt); and the")
    log("  BEAT 3 / BEAT 5 / Step01 re-point off SComatic REF_TRI (Section 7 / F32).")
    with open(os.path.join(OUTPUT_DIR, "build_fasta_trinucleotide_report.txt"), 'w') as fh:
        fh.write('\n'.join(report))


if __name__ == "__main__":
    main()
