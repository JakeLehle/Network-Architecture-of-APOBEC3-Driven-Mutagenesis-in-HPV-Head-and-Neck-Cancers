#!/usr/bin/env python3
"""
Diagnostic_TCW_Tri_vs_SComatic.py
==================================
DECISIVE test. Two independent SComatic-based recomputes (BEAT 5 ~9.4%,
BEAT 6 10.1%) disagree with the ranking file's own `tri` column (3.7%, 13/355).
The v2 audit only proved the ranking file is SELF-consistent (classify on its
own tri == its own is_tcw); it never checked that tri against the genome.

For every differential variant, look up the canonical SComatic REF_TRI at the
confirmed no-offset position (chrom, pos) and:
  1. check sco_ref == genomic ref (position/strand sanity),
  2. check file_tri == sco_tri (agreement) and file_tri == reverse(sco_tri)
     (the flank-swap signature that would explain the whole discrepancy),
  3. recount TCW under BOTH tris (file vs SComatic),
  4. re-evaluate the four figure genes under the SComatic (truth) tri --
     KLF3's clean-A3 status depends on it (KLF3 = G>A, file tri TGA; reversed
     AGT is NOT TCW).

Read-only. SComatic REF_TRI is the ground truth (Step01/BEAT1/BEAT3/BEAT5 all
use it); the ranking file's `tri` column is the thing under test.

Run:
    conda activate NETWORK
    python Diagnostic_TCW_Tri_vs_SComatic.py
"""

import os
import re
import numpy as np
import pandas as pd
from datetime import datetime

BASE_DIR  = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG7_ROOT = os.path.join(BASE_DIR, "data/FIG_7")
RANK_DIR  = os.path.join(FIG7_ROOT, "TROUBLESHOOTING/tcw_neoantigen_ranking")
DIFF_FILE = os.path.join(RANK_DIR, "differential_variants_tier12_ALL.tsv")
OUT_FILE  = os.path.join(RANK_DIR, "tri_vs_scomatic_reconciled.tsv")

SCOMATIC_TSV = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
                "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

FIGURE_GENES = ['KLF3', 'CAST', 'FABP5', 'ANXA1', 'ALKBH3', 'CYTH1', 'KRT6B', 'SERPINB2']
_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def log(msg):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)

def banner(t, c="="):
    print(f"\n{c * 80}\n  {t}\n{c * 80}", flush=True)


def classify_tcw(ref, alt, tri):
    if tri is None or len(str(tri)) != 3:
        return False
    tri = str(tri).upper()
    if ref == 'C' and alt in ('T', 'G'):
        return tri[0] == 'T' and tri[2] in ('A', 'T')
    if ref == 'G' and alt in ('A', 'C'):
        rc = ''.join(_COMP.get(b, 'N') for b in reversed(tri))
        return rc[0] == 'T' and rc[2] in ('A', 'T')
    return False


def classify_tcw_ct(ref, alt, tri):
    if tri is None or len(str(tri)) != 3:
        return False
    tri = str(tri).upper()
    if ref == 'C' and alt == 'T':
        return tri[0] == 'T' and tri[2] in ('A', 'T')
    if ref == 'G' and alt == 'A':
        rc = ''.join(_COMP.get(b, 'N') for b in reversed(tri))
        return rc[0] == 'T' and rc[2] in ('A', 'T')
    return False


def parse_ref_alt(row):
    for col in ('genomic_sub', 'oriented', 'sbs_class'):
        v = row.get(col)
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        m = re.search(r'([ACGT])\s*[>/]\s*([ACGT])', str(v).upper())
        if m:
            return m.group(1), m.group(2)
    return None, None


def load_scomatic(needed):
    """(chrom, Start) -> (REF, REF_TRI) for the requested keys."""
    idx = {}
    with open(SCOMATIC_TSV) as f:
        h = f.readline().rstrip('\n').split('\t')
        ci, si, ri, ti = (h.index('#CHROM'), h.index('Start'),
                          h.index('REF'), h.index('REF_TRI'))
        for line in f:
            fl = line.rstrip('\n').split('\t')
            try:
                key = (fl[ci], int(fl[si]))
            except (ValueError, IndexError):
                continue
            if key in needed and key not in idx:
                idx[key] = (fl[ri].upper(), fl[ti].upper())
                if len(idx) == len(needed):
                    break
    return idx


def main():
    banner("DECISIVE: ranking-file tri vs canonical SComatic REF_TRI")

    df = pd.read_csv(DIFF_FILE, sep="\t")
    df['REF'], df['ALT'] = zip(*df.apply(parse_ref_alt, axis=1))
    df['FILE_TRI'] = df['tri'].astype(str).str.upper()
    df['GENE'] = df['gene'].astype(str)
    parts = df['location'].astype(str).str.split(':', expand=True)
    df['CHROM'] = parts[0]
    df['POS'] = pd.to_numeric(parts[1], errors='coerce').astype('Int64')

    needed = set((r['CHROM'], int(r['POS'])) for _, r in df.iterrows() if pd.notna(r['POS']))
    log(f"  variants: {len(df)}; resolving SComatic REF_TRI at (chrom, pos) [no offset] ...")
    sco = load_scomatic(needed)
    log(f"  resolved {len(sco)}/{len(needed)}")

    rows = []
    for _, r in df.iterrows():
        srec = sco.get((r['CHROM'], int(r['POS']))) if pd.notna(r['POS']) else None
        sco_ref, sco_tri = (srec if srec else (None, None))
        rel = 'unresolved'
        if sco_tri:
            if sco_tri == r['FILE_TRI']:
                rel = 'equal'
            elif sco_tri == r['FILE_TRI'][::-1]:
                rel = 'reversed'
            else:
                rel = 'other'
        rows.append({
            'gene': r['GENE'], 'chrom': r['CHROM'], 'pos': r['POS'],
            'ref': r['REF'], 'alt': r['ALT'],
            'file_tri': r['FILE_TRI'], 'sco_ref': sco_ref, 'sco_tri': sco_tri,
            'rel': rel,
            'file_is_tcw': bool(r.get('is_tcw', False)),
            'tcw_FILE': classify_tcw(r['REF'], r['ALT'], r['FILE_TRI']),
            'tcw_SCO':  classify_tcw(r['REF'], r['ALT'], sco_tri),
            'ct_SCO':   classify_tcw_ct(r['REF'], r['ALT'], sco_tri),
        })
    res = pd.DataFrame(rows)
    ok = res[res['sco_tri'].notna()]

    banner("POSITION / STRAND SANITY", "-")
    ref_match = int((ok['sco_ref'] == ok['ref']).sum())
    center_match = int(ok['sco_tri'].apply(lambda t: t[1]).eq(ok['ref']).sum())
    log(f"  resolved {len(ok)}/{len(res)}")
    log(f"  SComatic REF == genomic ref : {ref_match}/{len(ok)}")
    log(f"  SComatic tri center == ref  : {center_match}/{len(ok)}")

    banner("FILE tri vs SComatic tri", "-")
    for k in ('equal', 'reversed', 'other', 'unresolved'):
        log(f"  {k:11s}: {int((res['rel'] == k).sum())}")
    log("  -> a large 'reversed' bucket means the ranking diagnostic swapped the")
    log("     trinucleotide flanks (bug); 'equal' would vindicate the file.")

    banner("TCW RECOUNT: file tri vs SComatic tri (truth)")
    log(f"  TCW under FILE tri (C>T+C>G): {int(res['tcw_FILE'].sum())}   "
        f"[= file is_tcw 13, self-consistent]")
    log(f"  TCW under SCO  tri (C>T+C>G): {int(ok['tcw_SCO'].sum())}   "
        f"[BEAT 5/6 SComatic route ~= 36]")
    log(f"  clean C>T under SCO tri     : {int(ok['ct_SCO'].sum())}")
    agree = int((res['tcw_FILE'] == res['tcw_SCO']).sum())
    log(f"  rows where file and SComatic agree on TCW: {agree}/{len(res)}")

    banner("FIGURE-GENE STATUS UNDER TRUTH (SComatic tri)", "-")
    log(f"  {'gene':>9s}  {'sub':>4s}  {'file_tri':>8s}  {'sco_tri':>8s}  "
        f"{'rel':>8s}  {'tcw_file':>8s}  {'tcw_SCO':>7s}  {'cleanCT':>7s}")
    for g in FIGURE_GENES:
        for _, r in res[res['gene'] == g].iterrows():
            log(f"  {g:>9s}  {str(r['ref'])+'>'+str(r['alt']):>4s}  "
                f"{str(r['file_tri']):>8s}  {str(r['sco_tri']):>8s}  {r['rel']:>8s}  "
                f"{str(r['tcw_FILE']):>8s}  {str(r['tcw_SCO']):>7s}  {str(r['ct_SCO']):>7s}")

    banner("VARIANTS WHERE FILE AND TRUTH DISAGREE ON TCW", "-")
    flip = res[res['tcw_FILE'] != res['tcw_SCO']].sort_values('gene')
    log(f"  {len(flip)} disagreements (file_tri call != SComatic_tri call):")
    for _, r in flip.iterrows():
        log(f"    {r['gene']:>10s}  {r['chrom']}:{int(r['pos']) if pd.notna(r['pos']) else '?'}  "
            f"{r['ref']}>{r['alt']}  file={r['file_tri']}({r['tcw_FILE']})  "
            f"sco={r['sco_tri']}({r['tcw_SCO']})  {r['rel']}")

    res.to_csv(OUT_FILE, sep="\t", index=False)
    log(f"\n  Full reconciliation written: {OUT_FILE}")

    banner("VERDICT")
    log("  If 'reversed' dominates and TCW-under-SCO ~= 36, the ranking file's tri")
    log("  is flank-swapped: the true differential-TCW count is the SComatic one,")
    log("  BEAT 3's ANXA1 zero is correct, and the A3 lead must be re-picked from")
    log("  genes whose tcw_SCO is True (CAST survives; KLF3 likely does NOT).")
    log("  If 'equal' dominates, the file is right (13) and BEAT 6 has a separate bug.")


if __name__ == "__main__":
    main()
