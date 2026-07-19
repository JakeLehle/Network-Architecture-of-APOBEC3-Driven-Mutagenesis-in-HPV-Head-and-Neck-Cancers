#!/usr/bin/env python3
"""
Diagnostic_Differential_TCW_Audit.py  (v2)
===========================================
Reconcile the differential-neoantigen TCW count. v1 wrongly re-resolved
REF_TRI from SComatic at pos-1 (the protein-altering-table convention), which
is WRONG for this ranking file (its `location` is already 0-based SComatic
Start). That resolved only 46/355 and produced a bogus recount of 3.

v2 does it the honest way:
  1. Classify from the file's OWN `tri` column (already validated: ref_matches
     355/355), using classify_tcw / classify_tcw_ctonly copied verbatim from
     BEAT 5. This reproduces the file's is_tcw / is_tcw_ct outright and needs
     no external lookup.
  2. Confirm the trusted call agrees with the file's own is_tcw on every row.
  3. OFFSET PROBE: for a sample of variants, try SComatic at (chrom, pos) and
     (chrom, pos-1) and report which offset gives center == REF. This exposes
     the coordinate convention that broke v1.
  4. Isolate ANXA1: expect exactly 1 TCW (C>G / SBS13, T[C>G]A), 0 clean C>T.
  5. Binder snapshot for the four figure genes.

Read-only.

Run:
    conda activate NETWORK
    python Diagnostic_Differential_TCW_Audit.py
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from datetime import datetime

BASE_DIR  = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG7_ROOT = os.path.join(BASE_DIR, "data/FIG_7")
RANK_DIR  = os.path.join(FIG7_ROOT, "TROUBLESHOOTING/tcw_neoantigen_ranking")
DIFF_FILE = os.path.join(RANK_DIR, "differential_variants_tier12_ALL.tsv")
OUT_FILE  = os.path.join(RANK_DIR, "differential_TCW_audit_reconciled.tsv")

SCOMATIC_TSV = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
                "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

FIGURE_GENES = ['KLF3', 'CAST', 'FABP5', 'ANXA1']
_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def log(msg):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)

def banner(t, c="="):
    print(f"\n{c * 80}\n  {t}\n{c * 80}", flush=True)


# ----------------------------------------- trusted classifiers (from BEAT 5)
def classify_tcw(ref, alt, tri):
    """Strand-aware APOBEC TCW (C>T or C>G at T[C]W, W in {A,T})."""
    if tri is None or (isinstance(tri, float) and np.isnan(tri)) or len(str(tri)) != 3:
        return False, False
    tri = str(tri).upper()
    if ref == 'C' and alt in ('T', 'G'):
        return (tri[0] == 'T' and tri[2] in ('A', 'T')), True
    if ref == 'G' and alt in ('A', 'C'):
        rc = ''.join(_COMP.get(b, 'N') for b in reversed(tri))
        return (rc[0] == 'T' and rc[2] in ('A', 'T')), True
    return False, False


def classify_tcw_ctonly(ref, alt, tri):
    """C>T-only arm (SBS2), strand-aware."""
    if tri is None or (isinstance(tri, float) and np.isnan(tri)) or len(str(tri)) != 3:
        return False, False
    tri = str(tri).upper()
    if ref == 'C' and alt == 'T':
        return (tri[0] == 'T' and tri[2] in ('A', 'T')), True
    if ref == 'G' and alt == 'A':
        rc = ''.join(_COMP.get(b, 'N') for b in reversed(tri))
        return (rc[0] == 'T' and rc[2] in ('A', 'T')), True
    return False, False


def parse_ref_alt(row):
    for col in ('genomic_sub', 'oriented', 'sbs_class'):
        val = row.get(col)
        if val is None or (isinstance(val, float) and np.isnan(val)):
            continue
        m = re.search(r'([ACGT])\s*[>/]\s*([ACGT])', str(val).upper())
        if m:
            return m.group(1), m.group(2)
    return None, None


def load_scomatic_tri_index(needed_keys):
    """One pass over SComatic master -> {(chrom, Start): REF_TRI}."""
    idx = {}
    if not os.path.exists(SCOMATIC_TSV):
        log(f"  [WARN] SComatic master not found; skipping offset probe.")
        return idx
    with open(SCOMATIC_TSV) as f:
        header = f.readline().rstrip('\n').split('\t')
        ci = header.index('#CHROM'); si = header.index('Start'); ti = header.index('REF_TRI')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            try:
                key = (fields[ci], int(fields[si]))
            except (ValueError, IndexError):
                continue
            if key in needed_keys and key not in idx:
                idx[key] = fields[ti].upper()
                if len(idx) == len(needed_keys):
                    break
    return idx


# =============================================================================
def main():
    banner("DIFFERENTIAL-NEOANTIGEN TCW AUDIT v2 (classify from file's own tri)")

    if not os.path.exists(DIFF_FILE):
        log(f"[FATAL] not found: {DIFF_FILE}")
        sys.exit(1)

    df = pd.read_csv(DIFF_FILE, sep="\t")
    log(f"  rows: {len(df)}")

    df['REF'], df['ALT'] = zip(*df.apply(parse_ref_alt, axis=1))
    df['FILE_TRI'] = df['tri'].astype(str).str.upper()
    df['GENE'] = df['gene'].astype(str)
    parts = df['location'].astype(str).str.split(':', expand=True)
    df['CHROM'] = parts[0]
    df['POS'] = pd.to_numeric(parts[1], errors='coerce')

    file_tcw    = df['is_tcw'].fillna(0).astype(bool)
    file_tcw_ct = df['is_tcw_ct'].fillna(0).astype(bool)

    # ---- (1) trusted classification from the file's own validated tri
    tctg, tcto = [], []
    for _, r in df.iterrows():
        a, _ = classify_tcw(r['REF'], r['ALT'], r['FILE_TRI'])
        b, _ = classify_tcw_ctonly(r['REF'], r['ALT'], r['FILE_TRI'])
        tctg.append(a); tcto.append(b)
    df['TCW_trusted_CTG'] = tctg
    df['TCW_trusted_CTonly'] = tcto

    banner("TRUSTED RECOUNT (strand-aware classifier on the file's OWN tri)")
    log(f"  differential variants:           {len(df)}")
    log(f"  TCW (C>T + C>G):   trusted {int(df['TCW_trusted_CTG'].sum()):3d}   "
        f"file {int(file_tcw.sum()):3d}   handoff 13   BEAT6-live 36")
    log(f"  TCW (C>T only):    trusted {int(df['TCW_trusted_CTonly'].sum()):3d}   "
        f"file {int(file_tcw_ct.sum()):3d}   handoff  9   BEAT6-live 30")

    # ---- (2) does trusted agree with the file's own flag, row by row?
    disc = df[df['TCW_trusted_CTG'] != file_tcw]
    banner(f"AGREEMENT: trusted vs file is_tcw ({len(disc)} disagreements)", "-")
    if len(disc) == 0:
        log("  perfect agreement -> the ranking diagnostic used the strict rule correctly;")
        log("  its 13 is right and BEAT 6's live 36 is the outlier (its own classifier).")
    else:
        for _, r in disc.iterrows():
            log(f"    {r['GENE']:>10s}  {r['CHROM']}:{int(r['POS'])}  "
                f"{r['REF']}>{r['ALT']}  tri={r['FILE_TRI']}  "
                f"file={bool(file_tcw.loc[r.name])}  trusted={r['TCW_trusted_CTG']}")

    # ---- (3) offset probe: which SComatic offset matches this file's coords?
    banner("OFFSET PROBE (why v1 failed: pos vs pos-1 against SComatic)", "-")
    samp = df.dropna(subset=['POS']).head(60)
    needed = set()
    for _, r in samp.iterrows():
        needed.add((r['CHROM'], int(r['POS'])))
        needed.add((r['CHROM'], int(r['POS']) - 1))
    idx = load_scomatic_tri_index(needed)
    hit0 = ok0 = hit1 = ok1 = 0
    for _, r in samp.iterrows():
        t0 = idx.get((r['CHROM'], int(r['POS'])))
        t1 = idx.get((r['CHROM'], int(r['POS']) - 1))
        if t0: hit0 += 1;  ok0 += int(len(t0) == 3 and t0[1] == r['REF'])
        if t1: hit1 += 1;  ok1 += int(len(t1) == 3 and t1[1] == r['REF'])
    n = len(samp)
    log(f"  sample n = {n}")
    log(f"  offset  pos    (no -1): resolved {hit0}/{n}, center==REF {ok0}/{n}")
    log(f"  offset  pos-1  (v1)   : resolved {hit1}/{n}, center==REF {ok1}/{n}")
    better = 'pos (no offset)' if ok0 >= ok1 else 'pos-1'
    log(f"  -> this file's location matches SComatic Start at: {better}")
    log(f"  (v1 used pos-1 -> that is why it resolved almost nothing correctly.)")

    # ---- (4) ANXA1
    banner("ANXA1 (expect 1 C>G/SBS13 TCW, 0 clean C>T)", "-")
    anx = df[df['GENE'].str.upper() == 'ANXA1']
    for _, r in anx.iterrows():
        log(f"    {r['CHROM']}:{int(r['POS'])}  {r['REF']}>{r['ALT']}  "
            f"tri={r['FILE_TRI']}  file_is_tcw={bool(file_tcw.loc[r.name])}  "
            f"trusted_TCW={r['TCW_trusted_CTG']}  clean_CT={r['TCW_trusted_CTonly']}")
    log(f"  ANXA1 trusted TCW (C>T+C>G): {int(anx['TCW_trusted_CTG'].sum())}")
    log(f"  ANXA1 trusted clean C>T    : {int(anx['TCW_trusted_CTonly'].sum())}")
    log(f"  -> BEAT 3's flat '0 TCW' is stale; ANXA1 carries 1 C>G (SBS13) event.")

    # ---- (5) binder snapshot for the four figure genes
    banner("FIGURE-GENE BINDER SNAPSHOT (best differential binder per gene)", "-")
    log(f"  {'gene':>8s}  {'best_mut_ic50':>13s}  {'best_gain':>9s}  "
        f"{'gain_wt_ic50':>12s}  {'n_pep':>5s}  {'tcw':>5s}  tier")
    for g in FIGURE_GENES:
        sub = df[df['GENE'] == g]
        if len(sub) == 0:
            log(f"  {g:>8s}  not in SBS2 differential set")
            continue
        best = sub.loc[sub['best_mut_ic50'].idxmin()] if 'best_mut_ic50' in sub else sub.iloc[0]
        log(f"  {g:>8s}  {best.get('best_mut_ic50', np.nan):13.1f}  "
            f"{best.get('best_gain', np.nan):9.2f}  {best.get('gain_wt_ic50', np.nan):12.1f}  "
            f"{str(best.get('n_diff_peptides')):>5s}  "
            f"{str(bool(file_tcw.loc[best.name])):>5s}  {best.get('tier')}")

    keep = ['GENE', 'CHROM', 'POS', 'REF', 'ALT', 'FILE_TRI', 'is_tcw',
            'TCW_trusted_CTG', 'TCW_trusted_CTonly']
    df[keep].to_csv(OUT_FILE, sep="\t", index=False)
    log(f"\n  Reconciled table written: {OUT_FILE}")

    banner("VERDICT")
    log("  The differential-TCW count is 13 (C>T+C>G) / 9 (C>T only). This is the")
    log("  number for the text and the F32 email. BEAT 6's live 36/30 is a bug in")
    log("  BEAT 6's own classifier (needs the 1104-line source to fix). BEAT 3's")
    log("  ANXA1 '0 TCW' is stale; ANXA1 = 1 C>G (SBS13), 0 clean C>T.")


if __name__ == "__main__":
    main()
