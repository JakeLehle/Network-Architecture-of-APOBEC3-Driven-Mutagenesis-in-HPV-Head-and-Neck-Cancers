#!/usr/bin/env python3
"""
Diagnostic_Section45_Verify_Numbers.py
======================================
Read-only audit that refreshes the numbers flagged with a dagger in the Section
4.5 results draft, from current pipeline outputs, so they can be quoted as final.

Checks
------
1. Binder / strong-binder / differential peptide counts per population, to
   reconfirm the 2,370 vs 1,339 binders and refresh strong binders (< 50 nM) and
   differential (mut < 500, wt > 500).
2. TCW-context enrichment, computed THREE explicit ways because the draft gave
   it three inconsistent ways. Pick the one that matches the claim wording.
3. Raw fusion junction totals per population (and normal-adjacent), plus per cell.
4. (optional) Raw per-cell protein-altering burden (the depth-confounded
   cross-check that was dropped from the results); reported for completeness.

Nothing here is a manuscript number until it prints cleanly; each value is
labeled with the exact definition used.

Run in NETWORK conda env (needs scipy for Fisher exact).
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MHC_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/03_mhc_binding")
ANNOT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/02_snpeff_annotation")
FASTA_TSV = os.path.join(PROJECT_ROOT, "data/FIG_7/fasta_context/ref_tri_fasta.tsv")
FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis")
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
GENO_PATH = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
             "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

GROUPS = ["SBS2_HIGH", "CNV_HIGH"]
STRONG = 50.0
BIND = 500.0

_report = []
def log(m=""):
    print(m, flush=True); _report.append(str(m))
def banner(t, ch="="):
    log(""); log(ch * 80); log(f"  {t}"); log(ch * 80)
def _first(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in cols: return c
        if c.lower() in low: return low[c.lower()]
    return None
def _read(p):
    return pd.read_csv(p, sep="\t")


# =============================================================================
# CHECK 1: BINDER / STRONG / DIFFERENTIAL PEPTIDE COUNTS
# =============================================================================
def check1():
    banner("CHECK 1: binder / strong-binder / differential peptide counts")
    for grp in GROUPS:
        neo_p = os.path.join(MHC_DIR, f"{grp}_neoantigens.tsv")
        all_p = os.path.join(MHC_DIR, f"{grp}_all_peptide_results.tsv")
        if not os.path.exists(neo_p) or not os.path.exists(all_p):
            log(f"  [WARN] missing peptide tables for {grp}"); continue
        neo = _read(neo_p)
        allp = _read(all_p)
        m = _first(allp.columns, ["mut_ic50", "mut_IC50"])
        w = _first(allp.columns, ["wt_ic50", "wt_IC50"])
        allp[m] = pd.to_numeric(allp[m], errors="coerce")
        allp[w] = pd.to_numeric(allp[w], errors="coerce")

        n_binder = len(neo)                                   # rows in neoantigens.tsv
        # strong binders
        if "is_strong_binder" in neo.columns:
            n_strong_neo = int(neo["is_strong_binder"].astype(str).str.lower().isin(
                ["true", "1", "1.0", "yes"]).sum())
        else:
            n_strong_neo = np.nan
        n_strong_all = int((allp[m] < STRONG).sum())
        # differential
        if "is_differential" in neo.columns:
            n_diff_neo = int(neo["is_differential"].astype(str).str.lower().isin(
                ["true", "1", "1.0", "yes"]).sum())
        else:
            n_diff_neo = np.nan
        n_diff_all = int(((allp[m] < BIND) & (allp[w] >= BIND)).sum())

        log(f"\n  {grp}:")
        log(f"    binders (rows in neoantigens.tsv):            {n_binder}")
        log(f"    strong binders (neoantigens.is_strong):       {n_strong_neo}")
        log(f"    strong binders (all_peptide mut_ic50 < 50):   {n_strong_all}")
        log(f"    differential (neoantigens.is_differential):   {n_diff_neo}")
        log(f"    differential (all_peptide mut<500 & wt>500):  {n_diff_all}")
    log("\n  -> use the neoantigens.tsv counts to match the 2,370 / 1,339 framing;")
    log("     the all_peptide counts are an independent cross-check.")


# =============================================================================
# CHECK 2: TCW-CONTEXT ENRICHMENT (three definitions)
# =============================================================================
def _fisher(a, b, c, d):
    """2x2 [[a,b],[c,d]] -> (OR, p) via scipy; safe if scipy missing."""
    try:
        from scipy.stats import fisher_exact
        return fisher_exact([[a, b], [c, d]])
    except Exception:
        return (np.nan, np.nan)

def check2():
    banner("CHECK 2: TCW-context enrichment, SBS2-HIGH vs CNV-HIGH (3 definitions)")
    if not os.path.exists(FASTA_TSV):
        log(f"  [WARN] ref_tri_fasta not found: {FASTA_TSV}"); return
    ref = _read(FASTA_TSV)
    rk = {c: _first(ref.columns, [c]) for c in
          ["chrom", "pos", "ref", "alt", "gene", "hgvs_p", "sub_pyr", "is_tcw", "is_tcw_ct", "is_neo"]}
    for col in ["is_tcw", "is_tcw_ct"]:
        if rk[col]:
            ref[rk[col]] = ref[rk[col]].astype(str).str.lower().isin(["true", "1", "1.0", "yes"])

    # per-group protein-altering variant sets from Step02
    grp_vars = {}
    for grp in GROUPS:
        p = os.path.join(ANNOT_DIR, f"{grp}.somatic_protein_altering.tsv")
        if not os.path.exists(p):
            log(f"  [WARN] missing somatic set for {grp}: {p}"); return
        df = _read(p)
        gc = _first(df.columns, ["chrom", "#chrom", "chr"])
        pc = _first(df.columns, ["pos", "gpos", "Start"])
        rc = _first(df.columns, ["ref"])
        ac = _first(df.columns, ["alt"])
        if None in (gc, pc, rc, ac):
            log(f"  [WARN] {grp} somatic set missing chrom/pos/ref/alt "
                f"(cols: {list(df.columns)})"); return
        key = df.rename(columns={gc: "chrom", pc: "pos", rc: "ref", ac: "alt"})
        grp_vars[grp] = key[["chrom", "pos", "ref", "alt"]].drop_duplicates()

    # join each group's variants to ref_tri on locus
    refkey = ref.rename(columns={rk["chrom"]: "chrom", rk["pos"]: "pos",
                                 rk["ref"]: "ref", rk["alt"]: "alt"})
    joined = {}
    for grp in GROUPS:
        gv = grp_vars[grp].astype({"pos": "int64"}, errors="ignore")
        rf = refkey.astype({"pos": "int64"}, errors="ignore")
        j = gv.merge(rf, on=["chrom", "pos", "ref", "alt"], how="left")
        n_join = j[rk["is_tcw"]].notna().sum() if rk["is_tcw"] else 0
        log(f"  {grp}: {len(gv)} protein-altering variants, {n_join} joined to ref_tri")
        joined[grp] = j

    def frac(df, mask_denom, mask_num):
        d = int(mask_denom.sum()); n = int((mask_denom & mask_num).sum())
        return n, d, (100.0 * n / d if d else np.nan)

    log("\n  Definition A: among C>T protein-altering variants, fraction in TCW (SBS2 signature)")
    A = {}
    for grp in GROUPS:
        j = joined[grp]
        sp = j[rk["sub_pyr"]].astype(str) if rk["sub_pyr"] else pd.Series([""] * len(j))
        ct = sp == "C>T"
        tcw = j[rk["is_tcw_ct"]].fillna(False) if rk["is_tcw_ct"] else pd.Series([False] * len(j))
        n, d, pct = frac(j, ct, tcw)
        A[grp] = (n, d)
        log(f"    {grp}: {n}/{d} = {pct:.1f}%")
    orv, pv = _fisher(A["SBS2_HIGH"][0], A["SBS2_HIGH"][1] - A["SBS2_HIGH"][0],
                      A["CNV_HIGH"][0], A["CNV_HIGH"][1] - A["CNV_HIGH"][0])
    log(f"    Fisher exact OR = {orv:.2f}, p = {pv:.4g}")

    log("\n  Definition B: among ALL protein-altering variants, fraction in TCW (any TCW)")
    B = {}
    for grp in GROUPS:
        j = joined[grp]
        allm = pd.Series([True] * len(j))
        tcw = j[rk["is_tcw"]].fillna(False) if rk["is_tcw"] else pd.Series([False] * len(j))
        n, d, pct = frac(j, allm, tcw)
        B[grp] = (n, d)
        log(f"    {grp}: {n}/{d} = {pct:.1f}%")
    orv, pv = _fisher(B["SBS2_HIGH"][0], B["SBS2_HIGH"][1] - B["SBS2_HIGH"][0],
                      B["CNV_HIGH"][0], B["CNV_HIGH"][1] - B["CNV_HIGH"][0])
    log(f"    Fisher exact OR = {orv:.2f}, p = {pv:.4g}")

    log("\n  Definition C: among NEOANTIGEN-forming variants, fraction in TCW (any TCW)")
    # neoantigen-forming = variant present in the group's neoantigens.tsv
    C = {}
    for grp in GROUPS:
        neo_p = os.path.join(MHC_DIR, f"{grp}_neoantigens.tsv")
        if not os.path.exists(neo_p):
            log(f"    [WARN] missing {grp} neoantigens"); continue
        neo = _read(neo_p)
        ng = _first(neo.columns, ["gene"]); nh = _first(neo.columns, ["hgvs_p"])
        neo_keys = set(zip(neo[ng].astype(str), neo[nh].astype(str)))
        j = joined[grp].copy()
        jg = j[rk["gene"]].astype(str) if rk["gene"] else pd.Series([""] * len(j))
        jh = j[rk["hgvs_p"]].astype(str) if rk["hgvs_p"] else pd.Series([""] * len(j))
        is_neo = pd.Series([(g, h) in neo_keys for g, h in zip(jg, jh)])
        tcw = j[rk["is_tcw"]].fillna(False) if rk["is_tcw"] else pd.Series([False] * len(j))
        n, d, pct = frac(j.reset_index(drop=True), is_neo, tcw.reset_index(drop=True))
        C[grp] = (n, d)
        log(f"    {grp}: {n}/{d} = {pct:.1f}%")
    if len(C) == 2:
        orv, pv = _fisher(C["SBS2_HIGH"][0], C["SBS2_HIGH"][1] - C["SBS2_HIGH"][0],
                          C["CNV_HIGH"][0], C["CNV_HIGH"][1] - C["CNV_HIGH"][0])
        log(f"    Fisher exact OR = {orv:.2f}, p = {pv:.4g}")
    log("\n  -> the draft's 'X% vs Y% of C>T variants' matches Definition A;")
    log("     the ANXA1-paragraph 'overall proportion of TCW' matches Definition B.")


# =============================================================================
# CHECK 3: RAW FUSION JUNCTION TOTALS
# =============================================================================
def check3():
    banner("CHECK 3: raw fusion junction totals per population")
    cand = ["per_group_junction_summary.tsv", "fusion_burden.tsv",
            "filtered_fusions.tsv", "all_fusions_filtered.tsv"]
    found = None
    for name in cand:
        p = os.path.join(FUSION_DIR, name)
        if os.path.exists(p):
            found = p; break
    if not found:
        log(f"  [WARN] no fusion summary found in {FUSION_DIR}. Files present:")
        if os.path.isdir(FUSION_DIR):
            for f in sorted(os.listdir(FUSION_DIR)):
                log(f"      {f}")
        log("  Point CHECK 3 at the right file and re-run.")
        return
    df = _read(found)
    log(f"  Using {os.path.basename(found)}: {df.shape}")
    log(f"  Columns: {list(df.columns)}")
    gcol = _first(df.columns, ["group", "population"])
    ncol = _first(df.columns, ["n_junctions", "junctions", "n_fusions", "count", "total"])
    if gcol and ncol:
        for _, r in df.iterrows():
            g = r[gcol]; n = r[ncol]
            per = n / 546.0
            log(f"    {g}: {n} junctions ({per:.2f} per cell)")
    else:
        log("  Could not auto-detect group/count columns; printing head so you can")
        log("  tell me which columns hold the per-group junction totals:")
        log(df.head(10).to_string())


# =============================================================================
# CHECK 4 (optional): RAW PER-CELL PROTEIN-ALTERING BURDEN (cross-check)
# =============================================================================
def check4():
    banner("CHECK 4 (optional): raw per-cell protein-altering burden [cross-check]")
    if not os.path.exists(GENO_PATH):
        log(f"  [WARN] genotype master not found; skipping."); return
    if not os.path.exists(FASTA_TSV):
        log(f"  [WARN] ref_tri_fasta not found; skipping."); return
    groups = _read(GROUP_PATH)
    bc = _first(groups.columns, ["cell_barcode", "CB", "barcode"])
    gp = _first(groups.columns, ["group"])
    cellsets = {g: set(groups.loc[groups[gp] == g, bc].astype(str)) for g in GROUPS}

    ref = _read(FASTA_TSV)
    rc = _first(ref.columns, ["chrom"]); rp = _first(ref.columns, ["pos"])
    pa_loci = set(zip(ref[rc].astype(str), pd.to_numeric(ref[rp], errors="coerce")))

    geno = _read(GENO_PATH)
    gch = _first(geno.columns, ["#CHROM", "chrom"]); gps = _first(geno.columns, ["Start", "pos"])
    gcb = _first(geno.columns, ["CB", "cell_barcode"])
    geno["_ch"] = geno[gch].astype(str); geno["_ps"] = pd.to_numeric(geno[gps], errors="coerce")

    # try coordinate offsets 0, +1, -1 and pick best locus-match rate
    best_off, best_rate = 0, -1
    for off in (0, 1, -1):
        hit = set(zip(geno["_ch"], geno["_ps"] + off)) & pa_loci
        rate = len(hit)
        if rate > best_rate:
            best_rate, best_off = rate, off
    log(f"  locus match offset = {best_off} (matched {best_rate} protein-altering loci)")
    geno["_key"] = list(zip(geno["_ch"], geno["_ps"] + best_off))
    geno_pa = geno[geno["_key"].isin(pa_loci)]

    for g in GROUPS:
        sub = geno_pa[geno_pa[gcb].astype(str).isin(cellsets[g])]
        n_events = len(sub)
        per_cell = n_events / max(len(cellsets[g]), 1)
        log(f"    {g}: {n_events} protein-altering carrier events / {len(cellsets[g])} cells "
            f"= {per_cell:.2f} per cell")
    log("  NOTE: this is the depth-confounded raw burden dropped from the results;")
    log("  the depth-corrected per-UMI value (0.343 vs 0.184) is the quoted one.")


def main():
    banner("SECTION 4.5 NUMBER AUDIT (refresh the daggered values)")
    check1()
    check2()
    check3()
    check4()
    banner("COMPLETE")
    out = os.path.join(PROJECT_ROOT, "data/FIG_7/05_summary/section45_number_audit.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write("\n".join(_report))
    log(f"  Report: {out}")


if __name__ == "__main__":
    main()
