#!/usr/bin/env python3
"""
Diagnostic_BRD4_CCDS_Isoform_Identity.py
================================================================================
Resolve which CCDS-anchored BRD4 transcript group corresponds to BRD4-L and
which (if either) corresponds to the BRD4-S isoform assayed by Wu et al.
Mol Cell 2024 (PMID 38103559).

WHY THIS EXISTS
---------------
Diagnostic_BRD4_Isoform_Feasibility_Probe.py and
Diagnostic_BRD4_Isoform_Population_Ratio.py group BRD4 transcripts by CCDS:

    L    = CCDS12328  (ENST00000263377, ENST00000679869)   20 exons
    S(a) = CCDS46004  (ENST00000371835)                     12 exons
    S(b) = CCDS82307  (ENST00000360016)                     12 exons

Both short groups have 12 exons, so exon count cannot tell them apart. The
population-ratio result is reported entirely in terms of S(a). If S(a) is NOT
the isoform Wu et al. call BRD4-S, that result is mislabelled and cannot be
cited against their paper.

METHOD
------
Sum CDS feature lengths per transcript from the same GTF the pipeline used.
GENCODE CDS features EXCLUDE the stop codon, so:

    protein length (aa) = sum(CDS lengths) / 3

Compare against the published lengths for UniProt O60885:
    O60885-1  BRD4-L (long)   1362 aa
    O60885-2  BRD4-S (short)   722 aa

These reference lengths are hardcoded below as REFERENCE_AA and are the ONLY
external facts this script relies on. Verify them at
https://www.uniprot.org/uniprotkb/O60885/entry (Sequence & Isoforms) before
trusting the verdict; if they are wrong, change REFERENCE_AA and re-run.

An optional UniProt REST lookup is attempted at the end purely as an
independent confirmation. It is allowed to fail; the GTF result stands alone.

Env: NETWORK (no extra dependencies beyond the standard library)
Usage:
    python Diagnostic_BRD4_CCDS_Isoform_Identity.py
    python Diagnostic_BRD4_CCDS_Isoform_Identity.py --gtf /path/to/other.gtf

Author: Jake Lehle / Claude
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
DEFAULT_GTF = "/master/jlehle/WORKING/SC/ref/GRCh38/genes/genes_unzipped.gtf"
OUTPUT_DIR = ("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/"
              "DIAGNOSTIC_BRD4_ISOFORM")
GENE_NAME = "BRD4"

# CCDS groups as used by the two isoform scripts. Bare accession, no version.
CCDS_GROUPS = {
    "L":    "CCDS12328",
    "S(a)": "CCDS46004",
    "S(b)": "CCDS82307",
}

# UniProt O60885 isoform lengths. VERIFY THESE before trusting the verdict.
REFERENCE_AA = {
    "BRD4-L (O60885-1)": 1362,
    "BRD4-S (O60885-2)": 722,
}

# How close a computed length must be to a reference length to be called a match
AA_TOLERANCE = 3

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')

report = []


def log(msg=""):
    print(msg, flush=True)
    report.append(str(msg))


def banner(title, char="="):
    log()
    log(char * 78)
    log(f"  {title}")
    log(char * 78)


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", default=DEFAULT_GTF)
    ap.add_argument("--gene", default=GENE_NAME)
    ap.add_argument("--no-uniprot", action="store_true",
                    help="skip the optional UniProt REST confirmation")
    args = ap.parse_args()

    banner(f"STEP 0: Read GTF and collect {args.gene} features")
    if not os.path.exists(args.gtf):
        log(f"  ERROR: GTF not found: {args.gtf}")
        sys.exit(1)
    log(f"  GTF: {args.gtf}")

    # tx -> accumulated info
    cds_len = defaultdict(int)
    stop_len = defaultdict(int)
    n_cds = defaultdict(int)
    tx_name, tx_ccds, tx_biotype, tx_strand, tx_chrom = {}, {}, {}, {}, {}
    n_lines = 0

    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ("CDS", "stop_codon", "transcript"):
                continue
            if f"{args.gene}" not in f[8]:
                continue
            a = parse_attrs(f[8])
            if a.get("gene_name") != args.gene:
                continue
            tx = a.get("transcript_id")
            if tx is None:
                continue
            n_lines += 1
            length = int(f[4]) - int(f[3]) + 1
            if feat == "CDS":
                cds_len[tx] += length
                n_cds[tx] += 1
            elif feat == "stop_codon":
                stop_len[tx] += length
            else:  # transcript
                tx_name[tx] = a.get("transcript_name", "")
                tx_biotype[tx] = a.get("transcript_type",
                                       a.get("transcript_biotype", ""))
                tx_strand[tx] = f[6]
                tx_chrom[tx] = f[0]
            # ccdsid can appear on any feature line
            if "ccdsid" in a:
                tx_ccds[tx] = a["ccdsid"]

    log(f"  matched {n_lines} feature lines across "
        f"{len(set(list(cds_len) + list(tx_name)))} transcript(s)")
    if not cds_len:
        log("  ERROR: no CDS features found. Check --gene and the GTF build.")
        sys.exit(1)

    banner("STEP 1: Per-transcript CDS length and derived protein length")
    log(f"  {'transcript_id':<20s} {'name':<12s} {'CCDS':<14s} "
        f"{'nCDS':>4s} {'CDSbp':>7s} {'stop':>4s} {'aa':>6s}")
    log(f"  {'-'*20} {'-'*12} {'-'*14} {'-'*4} {'-'*7} {'-'*4} {'-'*6}")

    rows = []
    for tx in sorted(cds_len, key=lambda t: -cds_len[t]):
        bp = cds_len[tx]
        aa = bp / 3.0
        aa_disp = f"{aa:.1f}" if abs(aa - round(aa)) > 1e-6 else f"{int(round(aa))}"
        ccds = tx_ccds.get(tx, "-")
        rows.append({
            "transcript_id": tx,
            "transcript_name": tx_name.get(tx, ""),
            "ccdsid": ccds,
            "ccds_bare": ccds.split(".")[0] if ccds != "-" else "-",
            "n_cds_exons": n_cds[tx],
            "cds_bp": bp,
            "stop_bp": stop_len.get(tx, 0),
            "protein_aa": aa,
            "biotype": tx_biotype.get(tx, ""),
            "strand": tx_strand.get(tx, ""),
        })
        log(f"  {tx:<20s} {tx_name.get(tx,''):<12s} {ccds:<14s} "
            f"{n_cds[tx]:>4d} {bp:>7d} {stop_len.get(tx,0):>4d} {aa_disp:>6s}")

    for r in rows:
        if r["cds_bp"] % 3 != 0:
            log(f"  WARNING: {r['transcript_id']} CDS length {r['cds_bp']} bp is "
                f"not a multiple of 3 (incomplete CDS / cds_start_NF?). "
                f"Its aa length is unreliable.")

    banner("STEP 2: Collapse to the CCDS groups used by the isoform scripts")
    group_aa = {}
    for label, acc in CCDS_GROUPS.items():
        members = [r for r in rows if r["ccds_bare"] == acc]
        if not members:
            log(f"  {label:<5s} [{acc}]: NO TRANSCRIPT FOUND -- the isoform "
                f"scripts and this GTF disagree. Stop and reconcile.")
            group_aa[label] = None
            continue
        aas = sorted({round(m['protein_aa'], 1) for m in members})
        names = ", ".join(f"{m['transcript_id']}({m['transcript_name']})"
                          for m in members)
        log(f"  {label:<5s} [{acc}]: {names}")
        log(f"        CDS -> protein length(s): {aas} aa")
        if len(aas) > 1:
            log(f"        WARNING: members of one CCDS group disagree on length. "
                f"The group is not a single protein product.")
        group_aa[label] = aas[0]

    banner("STEP 3: VERDICT -- map CCDS groups onto BRD4-L / BRD4-S")
    log("  Reference lengths (VERIFY at uniprot.org/uniprotkb/O60885):")
    for k, v in REFERENCE_AA.items():
        log(f"    {k:<22s} {v} aa")
    log()

    assignment = {}
    for label, aa in group_aa.items():
        if aa is None:
            assignment[label] = "NOT FOUND"
            continue
        hit = None
        for name, ref in REFERENCE_AA.items():
            if abs(aa - ref) <= AA_TOLERANCE:
                hit = name
                break
        assignment[label] = hit if hit else f"no reference match ({aa:g} aa)"
        log(f"  {label:<5s} = {aa:>7.6g} aa  ->  {assignment[label]}")

    log()
    short_hits = [lab for lab, v in assignment.items()
                  if isinstance(v, str) and v.startswith("BRD4-S")]
    long_hits = [lab for lab, v in assignment.items()
                 if isinstance(v, str) and v.startswith("BRD4-L")]

    if len(long_hits) == 1 and len(short_hits) == 1:
        sh = short_hits[0]
        log(f"  RESOLVED: {long_hits[0]} is BRD4-L and {sh} is BRD4-S.")
        if sh == "S(a)":
            log("  -> The population-ratio result IS reported on the correct isoform.")
            log("     S(a) counts (SBS2 50 / CNV 78 / Normal 28) stand as BRD4-S.")
        else:
            log("  -> The population-ratio result is reported on the WRONG group.")
            log("     BRD4-S is S(b), whose counts are 0 / 2 / 0. There is")
            log("     effectively NO BRD4-S signal and the isoform claim cannot")
            log("     be made. Fall back to citation-only (option A).")
    elif not short_hits:
        log("  UNRESOLVED: neither short group matches the BRD4-S reference length.")
        log("  Either REFERENCE_AA is wrong, or this GTF does not annotate the")
        log("  Wu et al. BRD4-S isoform at all. Do not label either group")
        log("  'BRD4-S' in the manuscript until this is settled.")
    else:
        log("  AMBIGUOUS: more than one group matched a reference length.")
        log("  Inspect the per-transcript table above before proceeding.")

    # ---- optional independent confirmation --------------------------------
    if not args.no_uniprot:
        banner("STEP 4: Optional UniProt REST confirmation (allowed to fail)", char="-")
        try:
            import json
            import urllib.request
            url = "https://rest.uniprot.org/uniprotkb/O60885.json"
            with urllib.request.urlopen(url, timeout=20) as resp:
                data = json.load(resp)
            log("  UniProt O60885 CCDS cross-references:")
            found = False
            for xref in data.get("uniProtKBCrossReferences", []):
                if xref.get("database") != "CCDS":
                    continue
                found = True
                iso = ",".join(xref.get("isoformId", "")
                               if isinstance(xref.get("isoformId"), str)
                               else xref.get("isoformId", []) or ["(canonical)"])
                log(f"    {xref.get('id'):<14s} isoform: {iso or '(canonical)'}")
            if not found:
                log("    no CCDS cross-references returned")
        except Exception as e:
            log(f"  UniProt lookup unavailable ({type(e).__name__}: {e}).")
            log("  This is expected on a firewalled node. The STEP 3 verdict")
            log("  does not depend on it; confirm REFERENCE_AA manually instead.")

    # ---- write outputs ----------------------------------------------------
    banner("COMPLETE")
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        tsv = os.path.join(OUTPUT_DIR, "brd4_ccds_isoform_identity.tsv")
        cols = ["transcript_id", "transcript_name", "ccdsid", "ccds_bare",
                "biotype", "strand", "n_cds_exons", "cds_bp", "stop_bp",
                "protein_aa"]
        with open(tsv, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        txt = os.path.join(OUTPUT_DIR, "brd4_ccds_isoform_identity.txt")
        with open(txt, "w") as fh:
            fh.write("\n".join(report))
        log(f"  [SAVE] {tsv}")
        log(f"  [SAVE] {txt}")
    except Exception as e:
        log(f"  could not write outputs ({e}); console output above is the record")


if __name__ == "__main__":
    main()
