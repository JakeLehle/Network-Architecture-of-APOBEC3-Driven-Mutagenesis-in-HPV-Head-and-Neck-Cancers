#!/usr/bin/env python3
"""
resolve_KRT6B_codon.py  (read-only)
===================================
Resolve the full Glu342 codon of KRT6B p.Glu342Lys to distinguish GAA>AAA from
GAG>AAG. The trinucleotide in ref_tri_fasta covers positions 521-523, but the
third codon base is one position lower, so this reads one extra reference base.

Why the codon is fixed this way: Glu (GAA/GAG) and Lys (AAA/AAG) codons differ
only at position 1 (G vs A), so a single-base Glu>Lys change is a position-1 G>A
on the coding strand. Thus the mutated genomic base is codon position 1. KRT6B is
minus-strand, so the codon reads 5'->3' from HIGH to LOW coordinate: positions
pos, pos-1, pos-2. The script self-checks by confirming the reference codon is a
Glu codon and the mutant codon is a Lys codon.

Env: NEOANTIGEN (needs pysam + genome FASTA). Read-only.
Run: conda run -n NEOANTIGEN python resolve_KRT6B_codon.py
"""
import pandas as pd
import pysam

GENOME_FA = "/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"
FASTA_TSV = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_7/fasta_context/ref_tri_fasta.tsv"
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
AA = {"GAA": "Glu", "GAG": "Glu", "AAA": "Lys", "AAG": "Lys"}


def build_codon(plus, pos, alt):
    """Coding (minus-strand) codon 5'->3' = coords pos, pos-1, pos-2; the mutated
    base (pos) is codon position 1. Returns (ref_codon, mut_codon)."""
    ref_codon = COMP[plus[pos]] + COMP[plus[pos - 1]] + COMP[plus[pos - 2]]
    mut_codon = COMP[alt] + COMP[plus[pos - 1]] + COMP[plus[pos - 2]]
    return ref_codon, mut_codon


def main():
    df = pd.read_csv(FASTA_TSV, sep="\t")
    r = df[(df["gene"].astype(str) == "KRT6B") &
           (df["hgvs_p"].astype(str).str.contains("Glu342Lys"))].iloc[0]
    chrom = str(r["chrom"])
    pos = int(r["pos"])
    ref = str(r["ref"]).upper()
    alt = str(r["alt"]).upper()
    print(f"locus {chrom}:{pos}  ref {ref} -> alt {alt}   "
          f"tri_pyr={r['tri_pyr']}  sub_pyr={r['sub_pyr']}  is_tcw_ct={r['is_tcw_ct']}")

    fa = pysam.FastaFile(GENOME_FA)
    lo, hi = pos - 3, pos + 3                      # 1-based inclusive window
    win = fa.fetch(chrom, lo - 1, hi).upper()      # 0-based half-open
    plus = {lo + i: b for i, b in enumerate(win)}
    assert plus[pos] == ref, f"genome base {plus[pos]} != table ref {ref} at {pos}"

    print("\n  coord      plus   minus(coding)")
    for c in range(hi, lo - 1, -1):                # high->low = coding 5'->3'
        tag = "  <== mutated (codon position 1)" if c == pos else ""
        print(f"  {c}    {plus[c]}      {COMP[plus[c]]}{tag}")

    ref_codon, mut_codon = build_codon(plus, pos, alt)
    ok = AA.get(ref_codon) == "Glu" and AA.get(mut_codon) == "Lys"
    print(f"\n  Glu342 codon (coding 5'->3', coords {pos}, {pos - 1}, {pos - 2}):")
    print(f"    reference codon: {ref_codon}  ({AA.get(ref_codon, '??')})")
    print(f"    mutant   codon: {mut_codon}  ({AA.get(mut_codon, '??')})")
    print(f"\n  RESULT: {ref_codon} > {mut_codon}   "
          f"({AA.get(ref_codon, '?')}342{AA.get(mut_codon, '?')})")
    tp = str(r["tri_pyr"])
    print(f"  SBS2 context (pyrimidine strand): {tp} > {tp[0]}T{tp[2]}")
    print(f"  self-check (ref=Glu and mut=Lys): {'PASS' if ok else '**FAIL**'}")


if __name__ == "__main__":
    main()
