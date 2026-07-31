#!/usr/bin/env python3
"""
Diagnostic_FullSpectrum_TCW_Enrichment.py   (read-only)
=======================================================
Does the TCW (TpCpW APOBEC) substitution fraction differ between SBS2-HIGH and
CNV-HIGH GENOME-WIDE, and does any genome-wide difference reach the
protein-altering compartment that actually produces neoantigens?

Why this exists
---------------
The Section 4.5 TCW-enrichment sentence (18.9% vs 9.0%, Fisher OR 2.36,
p=0.0037; and 17.8% vs 9.0% among neoantigen-forming variants) was computed on
the CODING subset using the SComatic REF_TRI field, which the FASTA arbiter
showed matches the genome only 23/355. Recomputed from genome-verified FASTA
context, that same coding subset is flat (11.7% vs 11.6% among C>T, OR 1.01,
p=1; see diagnostic_section7_numbers.py CHECK J). The open question is whether a
genome-wide TCW enrichment exists (it should, by the group definition) and
whether it propagates into the coding compartment.

What it does
------------
Computes BOTH compartments with ONE classifier so they are directly comparable:
  CODING       {group}.somatic_protein_altering.tsv  (already germline-subtracted)
  GENOME-WIDE  {group}.snpeff_all.tsv, deduped to unique (chrom,pos,ref,alt),
               minus NORMAL-background germline. Germline subtraction is
               replicated verbatim from Step02: a variant is germline if its
               (chrom, pos, alt) appears in NORMAL's annotated set.
Trinucleotide context is read from the GRCh38 FASTA via the classify_tcw used to
build ref_tri_fasta (copied verbatim). Fisher exact (the settled test) compares
SBS2-HIGH vs CNV-HIGH under three definitions:
  Def A  clean-TCW C>T among all C>T            (mirrors CHECK J Def A)
  Def B  TCW (C>T + C>G) among all SNVs         (mirrors CHECK J Def B)
  Def C  TCW (C>T + C>G) among C-class variants (the ORIGINAL sentence's denom)

READ-ONLY inputs. Env: NEOANTIGEN (needs pysam + genome FASTA), NOT NETWORK.
Run: conda run -n NEOANTIGEN python Diagnostic_FullSpectrum_TCW_Enrichment.py
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
ANNOTATION_DIR = os.path.join(BASE_DIR, "data/FIG_7/02_snpeff_annotation")
GENOME_FA = "/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"
OUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/fullspectrum_tcw")

GROUPS = ["SBS2_HIGH", "CNV_HIGH"]
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

_report = []
def log(m=""):
    print(m, flush=True)
    _report.append(str(m))
def banner(t, ch="="):
    log("")
    log(ch * 80)
    log(f"  {t}")
    log(ch * 80)

# =============================================================================
# TCW classification  (verbatim from Diagnostic_Build_FASTA_Trinucleotide_Table)
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

def resolve_contig(chrom, refs):
    if chrom in refs:
        return chrom
    alt = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
    return alt if alt in refs else None

# =============================================================================
# LOADING + CLASSIFICATION
# =============================================================================
def _first(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in cols:
            return c
        if c.lower() in low:
            return low[c.lower()]
    return None

def _load_variants(path):
    """Load a per-group table; return unique (chrom,pos,ref,alt), or None."""
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cc = _first(df.columns, ["chrom", "#chrom", "chr"])
    pc = _first(df.columns, ["pos", "start"])
    rc = _first(df.columns, ["ref"])
    ac = _first(df.columns, ["alt"])
    if not all([cc, pc, rc, ac]):
        return None
    df = df.rename(columns={cc: "chrom", pc: "pos", rc: "ref", ac: "alt"})
    df = df[["chrom", "pos", "ref", "alt"]].copy()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)
    return df.drop_duplicates()

def _snv_only(df):
    r = df["ref"].astype(str).str.upper()
    a = df["alt"].astype(str).str.upper()
    m = (r.str.len() == 1) & (a.str.len() == 1) & r.isin(list("ACGT")) & a.isin(list("ACGT"))
    return df[m].copy(), int((~m).sum())

def classify_frame(df, fa, refs):
    """Add sub_pyr / is_tcw / is_tcw_ct / ref_ok to a unique-variant frame."""
    sub, tcw, tcwct, okl = [], [], [], []
    for r in df.itertuples():
        chrom = resolve_contig(str(r.chrom), refs)
        if chrom is None:
            sub.append(None); tcw.append(False); tcwct.append(False); okl.append(False)
            continue
        try:
            tg = fa.fetch(chrom, int(r.pos) - 2, int(r.pos) + 1).upper()
        except Exception:
            tg = None
        centre = tg[1] if (tg and len(tg) == 3) else None
        ok = (centre == str(r.ref).upper())
        _, sp, t, tc = classify_tcw(tg if ok else None, str(r.ref), str(r.alt))
        sub.append(sp); tcw.append(bool(t)); tcwct.append(bool(tc)); okl.append(bool(ok))
    return df.assign(sub_pyr=sub, is_tcw=tcw, is_tcw_ct=tcwct, ref_ok=okl)

def counts(df):
    """Per-set counts feeding the three Fisher definitions (ref-matched only)."""
    d = df[df["ref_ok"]]
    return {
        "n_input": len(df),
        "n_ref_ok": len(d),
        "n_ct": int((d["sub_pyr"] == "C>T").sum()),
        "n_cclass": int(d["sub_pyr"].isin(["C>T", "C>G"]).sum()),
        "n_tcw": int(d["is_tcw"].sum()),
        "n_tcw_ct": int(d["is_tcw_ct"].sum()),
    }

# =============================================================================
# STATS + REPORTING
# =============================================================================
def fisher(a, b, cc, d):
    try:
        from scipy.stats import fisher_exact
    except Exception:
        return np.nan, np.nan
    if (a + b) == 0 or (cc + d) == 0:
        return np.nan, np.nan
    return fisher_exact([[a, b], [cc, d]])

def _pct(a, n):
    return f"{100 * a / n:.1f}%" if n else "n/a"

def run_compartment(name, cg):
    banner(f"COMPARTMENT: {name}", "-")
    log(f"  {'group':10s} {'SNV':>7s} {'C>T':>6s} {'Cclass':>7s} {'TCW':>5s} {'TCWct':>6s} "
        f"{'TCWct/C>T':>10s} {'TCW/SNV':>8s} {'TCW/Ccls':>9s}")
    for g in GROUPS:
        c = cg[g]
        fa_ = c['n_tcw_ct'] / c['n_ct'] if c['n_ct'] else float('nan')
        fb_ = c['n_tcw'] / c['n_ref_ok'] if c['n_ref_ok'] else float('nan')
        fc_ = c['n_tcw'] / c['n_cclass'] if c['n_cclass'] else float('nan')
        log(f"  {g:10s} {c['n_ref_ok']:7d} {c['n_ct']:6d} {c['n_cclass']:7d} "
            f"{c['n_tcw']:5d} {c['n_tcw_ct']:6d} {fa_:10.3f} {fb_:8.3f} {fc_:9.3f}")
    s, c = cg["SBS2_HIGH"], cg["CNV_HIGH"]
    # Def A: clean-TCW C>T among all C>T
    orv, p = fisher(s['n_tcw_ct'], s['n_ct'] - s['n_tcw_ct'], c['n_tcw_ct'], c['n_ct'] - c['n_tcw_ct'])
    log(f"\n  Def A  clean-TCW C>T among C>T:      SBS2 {s['n_tcw_ct']}/{s['n_ct']} "
        f"({_pct(s['n_tcw_ct'], s['n_ct'])}), CNV {c['n_tcw_ct']}/{c['n_ct']} "
        f"({_pct(c['n_tcw_ct'], c['n_ct'])})  Fisher OR={orv:.2f}, p={p:.3g}")
    # Def B: TCW among all SNVs
    orv, p = fisher(s['n_tcw'], s['n_ref_ok'] - s['n_tcw'], c['n_tcw'], c['n_ref_ok'] - c['n_tcw'])
    log(f"  Def B  TCW among all SNVs:           SBS2 {s['n_tcw']}/{s['n_ref_ok']} "
        f"({_pct(s['n_tcw'], s['n_ref_ok'])}), CNV {c['n_tcw']}/{c['n_ref_ok']} "
        f"({_pct(c['n_tcw'], c['n_ref_ok'])})  Fisher OR={orv:.2f}, p={p:.3g}")
    # Def C: TCW among C-class (original sentence denominator)
    orv, p = fisher(s['n_tcw'], s['n_cclass'] - s['n_tcw'], c['n_tcw'], c['n_cclass'] - c['n_tcw'])
    log(f"  Def C  TCW among C-class (C>T+C>G):  SBS2 {s['n_tcw']}/{s['n_cclass']} "
        f"({_pct(s['n_tcw'], s['n_cclass'])}), CNV {c['n_tcw']}/{c['n_cclass']} "
        f"({_pct(c['n_tcw'], c['n_cclass'])})  Fisher OR={orv:.2f}, p={p:.3g}")

def _defA_p(cg):
    s, c = cg["SBS2_HIGH"], cg["CNV_HIGH"]
    _, p = fisher(s['n_tcw_ct'], s['n_ct'] - s['n_tcw_ct'], c['n_tcw_ct'], c['n_ct'] - c['n_tcw_ct'])
    return p

def _defB_p(cg):
    s, c = cg["SBS2_HIGH"], cg["CNV_HIGH"]
    _, p = fisher(s['n_tcw'], s['n_ref_ok'] - s['n_tcw'], c['n_tcw'], c['n_ref_ok'] - c['n_tcw'])
    return p

# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("FULL-SPECTRUM TCW ENRICHMENT (SBS2-HIGH vs CNV-HIGH)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUT_DIR, exist_ok=True)

    try:
        import pysam
    except Exception as e:
        log(f"  [FATAL] pysam unavailable ({e}); run in the NEOANTIGEN env.")
        return
    if not os.path.isfile(GENOME_FA):
        log(f"  [FATAL] genome FASTA not found: {GENOME_FA}")
        return
    fa = pysam.FastaFile(GENOME_FA)
    refs = set(fa.references)

    # ---- CODING: somatic_protein_altering (already germline-subtracted) --------
    coding = {}
    for g in GROUPS:
        df = _load_variants(os.path.join(ANNOTATION_DIR, f"{g}.somatic_protein_altering.tsv"))
        if df is None:
            log(f"  [WARN] coding: missing {g}.somatic_protein_altering.tsv")
            coding = None
            break
        snv, n_non = _snv_only(df)
        coding[g] = counts(classify_frame(snv, fa, refs))
        log(f"  CODING {g}: {len(df)} unique variants, {len(snv)} SNVs "
            f"({n_non} non-SNV excluded), ref-match {coding[g]['n_ref_ok']}/{len(snv)}")

    # ---- GENOME-WIDE: snpeff_all minus NORMAL-background germline ---------------
    genome = None
    normal = _load_variants(os.path.join(ANNOTATION_DIR, "NORMAL.snpeff_all.tsv"))
    if normal is None:
        log("\n  [WARN] genome-wide: missing NORMAL.snpeff_all.tsv; cannot germline-subtract.")
    else:
        normal_keys = set(zip(normal["chrom"].astype(str),
                              normal["pos"].astype(int),
                              normal["alt"].astype(str)))
        log(f"\n  NORMAL germline keys (chrom,pos,alt): {len(normal_keys)}")
        genome = {}
        for g in GROUPS:
            df = _load_variants(os.path.join(ANNOTATION_DIR, f"{g}.snpeff_all.tsv"))
            if df is None:
                log(f"  [WARN] genome-wide: missing {g}.snpeff_all.tsv")
                genome = None
                break
            keys = list(zip(df["chrom"].astype(str), df["pos"].astype(int), df["alt"].astype(str)))
            is_germ = pd.Series([k in normal_keys for k in keys], index=df.index)
            som = df[~is_germ].copy()
            snv, n_non = _snv_only(som)
            genome[g] = counts(classify_frame(snv, fa, refs))
            log(f"  GENOME-WIDE {g}: {len(df)} annotated variants, {int(is_germ.sum())} germline removed, "
                f"{len(som)} somatic, {len(snv)} SNVs ({n_non} non-SNV excluded), "
                f"ref-match {genome[g]['n_ref_ok']}/{len(snv)}")

    # ---- report each compartment ----------------------------------------------
    if coding:
        run_compartment("CODING (protein-altering; produces neoantigens)", coding)
    if genome:
        run_compartment("GENOME-WIDE (all somatic SNVs)", genome)

    # ---- side-by-side read-out --------------------------------------------------
    if coding and genome:
        banner("READ-OUT: does the genome-wide signal reach the coding compartment?", "-")
        log(f"  Def A (clean-TCW C>T):  coding p={_defA_p(coding):.3g}   "
            f"genome-wide p={_defA_p(genome):.3g}")
        log(f"  Def B (TCW / all SNV):  coding p={_defB_p(coding):.3g}   "
            f"genome-wide p={_defB_p(genome):.3g}")
        log("")
        log("  Decide from the numbers above (do not assume the direction):")
        log("   - genome-wide SIGNIFICANT + coding FLAT -> APOBEC/TCW enrichment is a")
        log("     genome-wide property that does NOT concentrate in the protein-altering")
        log("     compartment, so it cannot be the mechanism behind the neoantigen excess.")
        log("     The excess is then a burden effect (more coding variants at equal TCW")
        log("     fraction). A rescued sentence would be a genome-wide / QC statement.")
        log("   - genome-wide FLAT as well -> no TCW enrichment to cite; drop the sentence.")
        log("   - coding SIGNIFICANT -> unexpected vs CHECK J; re-open the coding path.")

    # ---- persist ----------------------------------------------------------------
    rows = []
    for comp, cg in [("coding", coding), ("genome_wide", genome)]:
        if not cg:
            continue
        for g in GROUPS:
            row = {"compartment": comp, "group": g}
            row.update(cg[g])
            rows.append(row)
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "fullspectrum_tcw_counts.tsv"),
                                  sep="\t", index=False)
    with open(os.path.join(OUT_DIR, "fullspectrum_tcw_report.txt"), "w") as f:
        f.write("\n".join(_report))
    log(f"\n  Wrote counts + report to {OUT_DIR}")


if __name__ == "__main__":
    main()
