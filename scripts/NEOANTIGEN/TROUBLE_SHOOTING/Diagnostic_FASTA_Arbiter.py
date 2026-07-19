#!/usr/bin/env python3
"""
Diagnostic_FASTA_Arbiter.py
===========================
THE arbiter. SComatic REF_TRI (used by Step01 APOBEC summary, BEAT 5 population
signal, BEAT 6's 36) and the ranking file's `tri` (its 13, built from a
calibrated FASTA read) DISAGREE on the trinucleotide flanks. Both are centered
on ref, so neither self-check catches it. This goes straight to the GRCh38 FASTA
-- the ground truth both were trying to represent -- for all 355 differential
variants and decides which source is correct.

For each variant:
  1. Find where the ref allele actually sits in the FASTA: test pos-1 (pos as
     1-based) and pos (pos as 0-based). Report the genome-wide winner.
  2. Extract the TRUE trinucleotide centered on that base.
  3. Compare TRUE tri to SComatic REF_TRI and to the ranking file's tri.
  4. Recount TCW under the TRUE FASTA tri (the definitive differential-TCW count).
  5. Report the four figure genes (+ ALKBH3, CENPC as the swap pair) under truth.

Read-only. Needs pysam + GRCh38 FASTA (NEOANTIGEN env).

Run:
    conda run -n NEOANTIGEN python Diagnostic_FASTA_Arbiter.py
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from datetime import datetime

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG7      = os.path.join(PROJECT_ROOT, "data/FIG_7")
RANK_DIR  = os.path.join(FIG7, "TROUBLESHOOTING/tcw_neoantigen_ranking")
DIFF_FILE = os.path.join(RANK_DIR, "differential_variants_tier12_ALL.tsv")
OUT_FILE  = os.path.join(RANK_DIR, "fasta_arbiter_reconciled.tsv")

REF_GENOME   = "/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"
SCOMATIC_TSV = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
                "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

FIGURE_GENES = ['KLF3', 'CAST', 'FABP5', 'ANXA1', 'ALKBH3', 'CENPC', 'CYTH1', 'KRT6B']
_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def log(m=""):
    print(m, flush=True)

def banner(t, c="="):
    log(""); log(c * 80); log(f"  {t}"); log(c * 80)


def classify(ref, alt, tri):
    """genomic ref/alt, pyrimidine-oriented (matches Step01 / ranking).
    returns (is_tcw_CTG, is_tcw_CTonly)."""
    if tri is None or len(str(tri)) != 3 or any(b not in 'ACGT' for b in str(tri)):
        return (False, False)
    tri = str(tri).upper()
    if ref == 'C':
        palt, ptri = alt, tri
    elif ref == 'G':
        palt = _COMP.get(alt, 'N')
        ptri = ''.join(_COMP[b] for b in reversed(tri))
    else:
        return (False, False)
    motif = (ptri[0] == 'T' and ptri[2] in ('A', 'T'))
    return (motif and palt in ('T', 'G'), motif and palt == 'T')


def parse_ref_alt(row):
    for col in ('genomic_sub', 'oriented', 'sbs_class'):
        v = row.get(col)
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        m = re.search(r'([ACGT])\s*[>/]\s*([ACGT])', str(v).upper())
        if m:
            return m.group(1), m.group(2)
    return None, None


def load_scomatic_tri(needed):
    idx = {}
    with open(SCOMATIC_TSV) as f:
        h = f.readline().rstrip('\n').split('\t')
        ci, si, ti = h.index('#CHROM'), h.index('Start'), h.index('REF_TRI')
        for line in f:
            fl = line.rstrip('\n').split('\t')
            try:
                key = (fl[ci], int(fl[si]))
            except (ValueError, IndexError):
                continue
            if key in needed and key not in idx:
                idx[key] = fl[ti].upper()
                if len(idx) == len(needed):
                    break
    return idx


def main():
    banner("FASTA ARBITER: raw GRCh38 vs SComatic REF_TRI vs ranking-file tri")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")

    try:
        import pysam
    except ImportError:
        log("  [FATAL] pysam missing. Run in NEOANTIGEN env:")
        log("          conda run -n NEOANTIGEN python Diagnostic_FASTA_Arbiter.py")
        sys.exit(1)
    if not os.path.exists(REF_GENOME):
        log(f"  [FATAL] FASTA not found: {REF_GENOME}"); sys.exit(1)
    fa = pysam.FastaFile(REF_GENOME)
    refs = set(fa.references)

    def chrom_name(c):
        if c in refs: return c
        if c.startswith('chr') and c[3:] in refs: return c[3:]
        if ('chr' + c) in refs: return 'chr' + c
        return None

    df = pd.read_csv(DIFF_FILE, sep="\t")
    df['REF'], df['ALT'] = zip(*df.apply(parse_ref_alt, axis=1))
    df['FILE_TRI'] = df['tri'].astype(str).str.upper()
    df['GENE'] = df['gene'].astype(str)
    p = df['location'].astype(str).str.split(':', expand=True)
    df['CHROM'] = p[0]
    df['POS'] = pd.to_numeric(p[1], errors='coerce').astype('Int64')

    sco = load_scomatic_tri(set((r['CHROM'], int(r['POS']))
                                for _, r in df.iterrows() if pd.notna(r['POS'])))

    # ---- per-variant resolution against the raw FASTA
    rows = []
    hit_posm1 = hit_pos = hit_neither = 0
    for _, r in df.iterrows():
        c = chrom_name(r['CHROM'])
        if c is None or pd.isna(r['POS']):
            rows.append({**r_base(r), 'true_tri': None, 'convention': 'no_chrom'})
            hit_neither += 1
            continue
        pos = int(r['POS'])
        b_m1 = fa.fetch(c, pos - 1, pos).upper()    # pos as 1-based -> ref at 0-based pos-1
        b_p0 = fa.fetch(c, pos, pos + 1).upper()    # pos as 0-based -> ref at 0-based pos
        if b_m1 == r['REF']:
            center = pos - 1; conv = 'pos-1 (1-based)'; hit_posm1 += 1
        elif b_p0 == r['REF']:
            center = pos; conv = 'pos (0-based)'; hit_pos += 1
        else:
            rows.append({**r_base(r), 'true_tri': None, 'convention': 'ref_not_found',
                         'fasta_posm1': b_m1, 'fasta_pos': b_p0})
            hit_neither += 1
            continue
        true_tri = fa.fetch(c, center - 1, center + 2).upper()
        rows.append({**r_base(r), 'true_tri': true_tri, 'convention': conv,
                     'sco_tri': sco.get((r['CHROM'], pos))})
    res = pd.DataFrame(rows)

    banner("WHERE DOES THE REF ALLELE SIT IN THE GENOME?", "-")
    n = len(df)
    log(f"  pos-1 (annotation pos is 1-based): {hit_posm1}/{n}")
    log(f"  pos   (annotation pos is 0-based): {hit_pos}/{n}")
    log(f"  ref not found at either / no chrom: {hit_neither}/{n}")
    conv = 'pos-1 (1-based)' if hit_posm1 >= hit_pos else 'pos (0-based)'
    log(f"  -> genome-wide convention: {conv}")

    ok = res[res['true_tri'].notna()].copy()
    ok['sco_tri'] = ok['sco_tri'].astype(str).str.upper()

    banner("WHICH SOURCE MATCHES THE RAW GENOME?")
    eq_sco  = int((ok['sco_tri'] == ok['true_tri']).sum())
    eq_file = int((ok['FILE_TRI'] == ok['true_tri']).sum())
    log(f"  variants resolved against FASTA: {len(ok)}/{n}")
    log(f"  SComatic REF_TRI == true FASTA tri : {eq_sco}/{len(ok)}   "
        f"(BEAT 5/6 + population signal ride on this)")
    log(f"  ranking file tri == true FASTA tri : {eq_file}/{len(ok)}   "
        f"(the 13 rides on this)")
    winner = 'SComatic REF_TRI (=> ~36 is right, KLF3 out, population signal safe)' \
             if eq_sco > eq_file else \
             'ranking-file tri (=> 13 is right, KLF3 valid, SComatic REF_TRI is the bug)'
    log(f"  -> source that matches the genome: {winner}")

    banner("DEFINITIVE TCW COUNT (under raw FASTA tri)")
    ctg = ok.apply(lambda r: classify(r['REF'], r['ALT'], r['true_tri'])[0], axis=1)
    cto = ok.apply(lambda r: classify(r['REF'], r['ALT'], r['true_tri'])[1], axis=1)
    log(f"  differential variants resolved : {len(ok)}")
    log(f"  TCW (C>T + C>G) under FASTA     : {int(ctg.sum())}   [file said 13, SComatic route 36]")
    log(f"  clean C>T under FASTA           : {int(cto.sum())}   [file said 9,  SComatic route 30]")
    ok['tcw_FASTA'] = ctg.values

    banner("FIGURE GENES UNDER RAW FASTA", "-")
    log(f"  {'gene':>9s}  {'sub':>4s}  {'FASTA':>6s}  {'SComa':>6s}  {'file':>6s}  {'TCW(FASTA)':>10s}")
    for g in FIGURE_GENES:
        for _, r in ok[ok['GENE'] == g].iterrows():
            t = classify(r['REF'], r['ALT'], r['true_tri'])[0]
            log(f"  {g:>9s}  {str(r['REF'])+'>'+str(r['ALT']):>4s}  "
                f"{str(r['true_tri']):>6s}  {str(r['sco_tri']):>6s}  "
                f"{str(r['FILE_TRI']):>6s}  {str(t):>10s}")

    res.to_csv(OUT_FILE, sep="\t", index=False)
    log(f"\n  Full reconciliation written: {OUT_FILE}")

    banner("VERDICT")
    log("  The 'matches the genome' line above is the whole answer:")
    log("   - if SComatic wins: true count ~36, KLF3 is not TCW, population signal safe.")
    log("   - if ranking file wins: true count 13, KLF3 valid, and SComatic REF_TRI")
    log("     (hence BEAT 5's population signal) is the thing that needs fixing.")


def r_base(r):
    return {'GENE': r['GENE'], 'chrom': r['CHROM'],
            'pos': int(r['POS']) if pd.notna(r['POS']) else None,
            'REF': r['REF'], 'ALT': r['ALT'], 'FILE_TRI': r['FILE_TRI'],
            'is_tcw_file': bool(r.get('is_tcw', False))}


if __name__ == "__main__":
    main()
