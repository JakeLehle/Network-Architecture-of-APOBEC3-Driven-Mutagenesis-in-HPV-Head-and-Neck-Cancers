#!/usr/bin/env python3
"""
resolve_KRT6B_HLA.py  (read-only)
=================================
Back out which HLA allele in the 10-allele panel presents the KRT6B Glu342Lys
neoantigen, so the co-culture cell line and PBMC donor can be HLA-matched.

Part 1 reads the winning allele straight from the per-peptide MHCflurry table
(the pipeline stores best_allele = the lowest-IC50 allele for each peptide, and
mut_ic50 = that best value). It lists every candidate peptide for the mutation,
sorted by predicted affinity, and names the strongest binder and its allele.

Part 2 re-predicts that top mutant peptide (and its wild-type) against ALL 10
panel alleles with the same MHCflurry predictor the pipeline uses, so you get the
full per-allele profile: which allele presents it, how strongly, and whether more
than one allele can present it (useful for finding a matching cell line / donor).

Env: NEOANTIGEN (pandas; mhcflurry for Part 2). Read-only.
Run: conda run -n NEOANTIGEN python resolve_KRT6B_HLA.py
"""
import os
import pandas as pd

MHC_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_7/03_mhc_binding"
GENE = "KRT6B"
HGVS = "Glu342Lys"
GROUPS = ["SBS2_HIGH", "CNV_HIGH"]   # KRT6B is SBS2-specific; CNV as fallback
PANEL = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*24:02",
         "HLA-B*07:02", "HLA-B*08:01", "HLA-B*35:01", "HLA-B*44:02",
         "HLA-C*04:01", "HLA-C*07:01"]


def _first(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in cols:
            return c
        if c.lower() in low:
            return low[c.lower()]
    return None


def part1():
    for group in GROUPS:
        p = os.path.join(MHC_DIR, f"{group}_all_peptide_results.tsv")
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p, sep="\t")
        g = _first(df.columns, ["gene", "gene_symbol"])
        h = _first(df.columns, ["hgvs_p"])
        if not (g and h):
            continue
        sub = df[(df[g].astype(str) == GENE) &
                 (df[h].astype(str).str.contains(HGVS))].copy()
        if not len(sub):
            continue
        m = _first(sub.columns, ["mut_ic50", "mut_IC50"])
        sub[m] = pd.to_numeric(sub[m], errors="coerce")
        sub = sub.sort_values(m)
        show = [c for c in ["mut_peptide", "wt_peptide", "best_allele", "mut_ic50",
                            "wt_ic50", "is_differential", "mut_position_in_peptide"]
                if c in sub.columns]
        print(f"=== {GENE} {HGVS}: candidate peptides ({group}_all_peptide_results.tsv) ===")
        print(sub[show].to_string(index=False))
        top = sub.iloc[0]
        print(f"\n-> strongest predicted binder: {top['mut_peptide']}")
        print(f"   presented by {top.get('best_allele')} at {top[m]:.1f} nM "
              f"(wild-type {top.get('wt_ic50')} nM)")
        return str(top["mut_peptide"]), str(top.get("wt_peptide", "")), str(top.get("best_allele"))
    print(f"[!] no {GENE} {HGVS} rows found in any group's all_peptide_results.tsv")
    return None


def part2(mut_pep, wt_pep):
    try:
        from mhcflurry import Class1AffinityPredictor
    except Exception as e:
        print(f"\n[Part 2 skipped: mhcflurry unavailable ({e})]")
        return
    pred = Class1AffinityPredictor.load()

    def ic50(pep, allele):
        if not pep or str(pep).lower() == "nan":
            return float("nan")
        try:
            d = pred.predict_to_dataframe(peptides=[pep], alleles=[allele])
            return float(d.iloc[0]["prediction"])
        except Exception:
            return float("nan")

    print(f"\n=== full per-allele binding across the 10-allele panel ===")
    print(f"    mutant peptide: {mut_pep}    wild-type peptide: {wt_pep}")
    rows = []
    for a in PANEL:
        mm = ic50(mut_pep, a)
        ww = ic50(wt_pep, a)
        rows.append({"allele": a,
                     "mut_ic50": round(mm, 1) if mm == mm else mm,
                     "wt_ic50": round(ww, 1) if ww == ww else ww,
                     "mut_binder_<500": (mm < 500) if mm == mm else False,
                     "mut_strong_<50": (mm < 50) if mm == mm else False})
    prof = pd.DataFrame(rows).sort_values("mut_ic50")
    print(prof.to_string(index=False))
    best = prof.iloc[0]
    binders = prof[prof["mut_ic50"] < 500]["allele"].tolist()
    print(f"\n-> presenting HLA (strongest): {best['allele']} at {best['mut_ic50']} nM")
    print(f"   panel alleles that present this peptide (<500 nM): "
          f"{', '.join(binders) if binders else 'none'}")
    print("   NOTE: MHCflurry is a prediction; the co-culture + TCR-seq is the validation.")


def main():
    r = part1()
    if r:
        part2(r[0], r[1])


if __name__ == "__main__":
    main()
