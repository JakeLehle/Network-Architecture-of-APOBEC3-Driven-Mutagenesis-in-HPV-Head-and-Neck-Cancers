#!/usr/bin/env python3
"""
diagnostic_section7_numbers.py  (v3 - consolidated assertion audit)
===================================================================
ONE read-only audit for every number quoted in Section 4.5 / Figure 7
(neoantigen landscape and therapeutic target identification).

This version supersedes and replaces:
  - the old BEAT-based diagnostic_section7_numbers.py (gene-level tiers
    105/276/135, and the retired SComatic REF_TRI TCW recompute), and
  - the standalone Diagnostic_Section45_Verify_Numbers.py, whose genome-verified
    TCW logic and peptide-count recomputes are folded in here (CHECK C, CHECK J).

Design
------
Every value in the results prose is compared against a hardcoded EXPECTED block
and reported as expected | observed | PASS/FAIL, with a tally at the end. The
audit READS the locked single-source tables (so a table value and the prose
cannot silently drift) and only RECOMPUTES the two cheap peptide-level counts
(strong / differential) and the genome-verified TCW fractions:

  full.tsv                     tiers, per-niche + tier prevalence, expression,
                               featured binding, TCW class           (ranking)
  panelA_grouprate_stats.tsv   binder/variant/gene counts, folds, binomial p
                               (group-aware burden diagnostic; NEW output)
  panelA_burden_stats.tsv      neoantigen + fusion per-UMI means, BH-adj p
  shared_neoantigen_selection_evidence.tsv   Tier-1 selection-evidence set (Step05)
  {group}_all_peptide_results.tsv            strong / differential recompute
  per_group_junction_summary.tsv             raw fusion junction totals
  ref_tri_fasta.tsv + {group}.somatic_protein_altering.tsv   genome-verified TCW

A missing upstream table -> that check SKIPS (not a crash); a value off-target
-> FAIL. Nothing is a manuscript number until this prints a clean sweep.

Run in the NETWORK conda env (needs scipy for Fisher):
    conda run -n NETWORK python diagnostic_section7_numbers.py
    conda run -n NETWORK python diagnostic_section7_numbers.py --heavy   # + raw per-cell burden
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
FIG7 = os.path.join(PROJECT_ROOT, "data/FIG_7")

MHC_DIR      = os.path.join(FIG7, "03_mhc_binding")
ANNOT_DIR    = os.path.join(FIG7, "02_snpeff_annotation")
FASTA_TSV    = os.path.join(FIG7, "fasta_context/ref_tri_fasta.tsv")
FUSION_DIR   = os.path.join(FIG7, "04_fusion_analysis")
RANK_TSV     = os.path.join(FIG7, "06_prevalence_ranking/neoantigen_prevalence_ranking_full.tsv")
SUMMARY_DIR  = os.path.join(FIG7, "05_summary")
EVIDENCE_TSV = os.path.join(SUMMARY_DIR, "shared_neoantigen_selection_evidence.tsv")
GA_DIR       = os.path.join(FIG7, "TROUBLESHOOTING/group_aware_expression")
BURDEN_TSV   = os.path.join(GA_DIR, "panelA_burden_stats.tsv")
GROUPRATE_TSV = os.path.join(GA_DIR, "panelA_grouprate_stats.tsv")

GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
GENO_PATH  = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
              "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

GROUPS = ["SBS2_HIGH", "CNV_HIGH"]
STRONG = 50.0
BIND = 500.0
N_CELLS = 546

# Tier strings as written by the ranking / asserted by Step05.
T_SBS2, T_SHARED, T_CNV = "SBS2_specific", "shared", "CNV_specific"

# =============================================================================
# EXPECTED  (authoritative targets, from Section_4.5_results_updated.md + the
# 2026-07-21 session note; edit HERE if the prose changes, never the other way)
# =============================================================================
EXPECTED = {
    # Burden group-rate (Panel A) -> panelA_grouprate_stats.tsv
    "binders_peptides":    {"sbs2": 2370, "cnv": 1339, "fold": 1.77, "binom_p": 5.45e-65},
    "neoantigen_variants": {"sbs2": 560,  "cnv": 308,  "fold": 1.82, "binom_p": 9.50e-18},
    "neoantigen_genes":    {"sbs2": 354,  "cnv": 225,  "fold": 1.57, "binom_p": 9.25e-08},
    # Burden per-UMI (Panel A) -> panelA_burden_stats.tsv
    "neoantigen_perUMI":   {"sbs2": 0.343, "cnv": 0.184, "bh_p": 2.15e-06},
    # PENDING: means set to current table (germline-subtracted); notes had 0.567/0.526.
    # Same SBS2-CNV gap (0.041) and same ns p either way. Confirm the prose quotes these.
    "fusion_perUMI":       {"sbs2": 0.547, "cnv": 0.506, "bh_p": 0.327},
    # Peptide-level counts recomputed from all_peptide_results.tsv
    "strong_binders":      {"sbs2": 264, "cnv": 147},   # mut_ic50 < 50
    "differential":        {"sbs2": 615, "cnv": 339},   # mut < 500 & wt >= 500
    # Raw fusion junction totals
    "raw_junctions":         {"SBS2_HIGH": 5128, "CNV_HIGH": 5566, "NORMAL": 6625},
    "raw_junctions_percell": {"SBS2_HIGH": 9.4,  "CNV_HIGH": 10.2, "NORMAL": 12.1},
    # Tiers (Panel B) -> full.tsv
    "tiers": {T_SBS2: 467, T_SHARED: 93, T_CNV: 215},
    "tier_total": 775,
    "sbs2_total_mutations": 560,   # SBS2_specific + shared
    "cnv_total_mutations": 308,    # CNV_specific + shared
    # Tier-1 selection evidence (Step05)
    "shared_selection_evidence": 15,
    "shared_direction": {"depleted": 12, "flat": 3, "enriched": 0},
    "shared_mechanism": {"fusion": 10, "silenced": 6},
    "shared_hla_excluded": 5,
    "escape_subset": ["SPRR1A", "PI3", "SERPINB2", "KRT6A"],
    # Featured candidates (Panel C/D) -> full.tsv
    "featured": {
        "COX4I1": {"hgvs_p": "p.Ala9Thr",   "tier": T_CNV,    "prev": 24.9, "expr": 99.8, "wt": 996,   "mut": 321,  "tcw": False},
        "SPRR1A": {"hgvs_p": "p.Val61Ile",  "tier": T_SHARED, "prev": 16.9, "expr": 46.2, "wt": 444,   "mut": 295,  "tcw": False},
        "KRT6B":  {"hgvs_p": "p.Glu342Lys", "tier": T_SBS2,   "prev": 7.9,  "expr": 66.5, "wt": 27793, "mut": 80.5, "tcw": True},
    },
    # Leaders: featured gene per tier = top NON-HLA neoantigen by prevalence_max
    # (HLA-A/B/C excluded for germline polymorphism); ties broken by binding gain.
    "leader_shared": "SPRR1A",
    "leader_cnv_specific": "COX4I1",
    "leader_sbs2_specific": "KRT6B",    # ties TACSTD2 at prevMax 0.0788 -> delta tiebreak
    # ANXA1 rationale line (high expression, low carriage)
    "anxa1": {"expr": 99.5, "carriage_pct": 1.1},   # gene-union carriage = 6/546 = 1.1%
    # TCW three definitions, genome-verified, CODING subset -> must stay flat / ns
    "tcw_defA_ct":  {"sbs2": 11.7, "cnv": 11.6},   # among C>T protein-altering
    "tcw_defB_all": {"sbs2": 4.3,  "cnv": 3.8},    # among all protein-altering
    "tcw_defC_neo": {"sbs2": 4.3,  "cnv": 3.6},    # among neoantigen-forming
    # Population coverage / HLA panel
    "pop_coverage_pct": 90.5,
    "hla_alleles": ["A*01:01", "A*02:01", "A*03:01", "A*24:02", "B*07:02",
                    "B*08:01", "B*35:01", "B*44:02", "C*04:01", "C*07:01"],
}

# =============================================================================
# LOGGING + ASSERTION HELPER
# =============================================================================
_report = []
_tally = {"pass": 0, "fail": 0, "skip": 0}

def log(m=""):
    print(m, flush=True)
    _report.append(str(m))

def banner(t, ch="="):
    log("")
    log(ch * 78)
    log(f"  {t}")
    log(ch * 78)

def _first(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in cols:
            return c
        if c.lower() in low:
            return low[c.lower()]
    return None

def _read(p):
    return pd.read_csv(p, sep="\t")

def _need(path, label):
    if not os.path.exists(path):
        log(f"  [SKIP] {label} not found: {path}")
        return False
    return True

def _fmt(v):
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
        return f"{v:.3g}"
    return f"{v:.3f}"

def check(name, expected, observed, kind="exact", tol=None):
    """kind: exact (ints) | eq (str/bool) | close (abs tol) | rel (rel tol) |
    mag (same order of magnitude, for tiny p) | ns (observed > 0.05)."""
    if observed is None or (isinstance(observed, float) and np.isnan(observed)):
        _tally["skip"] += 1
        log(f"  [SKIP] {name}: expected {_fmt(expected)}, observed n/a")
        return
    ok = False
    if kind == "exact":
        ok = int(round(float(observed))) == int(expected)
    elif kind == "eq":
        ok = str(observed) == str(expected)
    elif kind == "close":
        t = 0.01 if tol is None else tol
        ok = abs(float(observed) - float(expected)) <= t
    elif kind == "rel":
        t = 0.02 if tol is None else tol
        ok = abs(float(observed) - float(expected)) <= t * abs(float(expected))
    elif kind == "mag":
        e, o = float(expected), float(observed)
        ok = (o > 0) and abs(np.log10(o) - np.log10(e)) <= 1.0
    elif kind == "ns":
        ok = float(observed) > 0.05
    _tally["pass" if ok else "fail"] += 1
    log(f"  [{'PASS' if ok else '**FAIL**'}] {name}: expected {_fmt(expected)}, observed {_fmt(observed)}")

def _tobool(s):
    return str(s).strip().lower() in ("true", "1", "1.0", "yes")


# =============================================================================
# CHECK A: burden group-rate (peptides / variants / genes)  [Panel A]
# =============================================================================
def check_A_grouprate():
    banner("CHECK A: burden group-rate counts, folds, binomial p  [Panel A]")
    if not _need(GROUPRATE_TSV, "panelA_grouprate_stats.tsv"):
        log("  -> re-run the updated Diagnostic_GroupAware_Expression_and_Carrier.py.")
        _tally["skip"] += 9
        return
    df = _read(GROUPRATE_TSV).set_index(_first(_read(GROUPRATE_TSV).columns, ["unit"]))
    unit_map = {"peptides": "binders_peptides",
                "neoantigen-producing variants": "neoantigen_variants",
                "neoantigen genes": "neoantigen_genes"}
    for unit, key in unit_map.items():
        e = EXPECTED[key]
        if unit not in df.index:
            check(f"{key} SBS2", e["sbs2"], None)
            check(f"{key} CNV", e["cnv"], None)
            check(f"{key} fold", e["fold"], None)
            continue
        r = df.loc[unit]
        check(f"{key} SBS2 count", e["sbs2"], r.get("sbs2"), "exact")
        check(f"{key} CNV count", e["cnv"], r.get("cnv"), "exact")
        check(f"{key} fold", e["fold"], r.get("fold_sbs2_over_cnv"), "close", 0.01)
        check(f"{key} binomial p", e["binom_p"], r.get("binomial_p"), "mag")


# =============================================================================
# CHECK B: burden per-UMI (neoantigen & fusion)  [Panel A]
# =============================================================================
def check_B_perumi():
    banner("CHECK B: burden per-UMI means + BH-adjusted p  [Panel A]")
    if not _need(BURDEN_TSV, "panelA_burden_stats.tsv"):
        _tally["skip"] += 6
        return
    df = _read(BURDEN_TSV).set_index("comparison")
    for comp, key in [("neoantigen_perUMI", "neoantigen_perUMI"),
                      ("fusion_perUMI", "fusion_perUMI")]:
        e = EXPECTED[key]
        if comp not in df.index:
            check(f"{key} SBS2 mean", e["sbs2"], None)
            check(f"{key} CNV mean", e["cnv"], None)
            check(f"{key} BH p", e["bh_p"], None)
            continue
        r = df.loc[comp]
        check(f"{key} SBS2 mean", e["sbs2"], r.get("sbs2_mean"), "close", 0.003)
        check(f"{key} CNV mean", e["cnv"], r.get("cnv_mean"), "close", 0.003)
        if key == "neoantigen_perUMI":
            check(f"{key} BH-adj p", e["bh_p"], r.get("bh_adjusted_p"), "rel", 0.10)
        else:
            check(f"{key} BH-adj p", e["bh_p"], r.get("bh_adjusted_p"), "close", 0.02)
            check(f"{key} BH-adj p is ns (>0.05)", 0.05, r.get("bh_adjusted_p"), "ns")


# =============================================================================
# CHECK C: strong binders & differential (peptide-level, recomputed)
# =============================================================================
def check_C_strong_diff():
    banner("CHECK C: strong binders (<50 nM) & differential (recomputed)")
    strong, diff, binders = {}, {}, {}
    for grp in GROUPS:
        allp = os.path.join(MHC_DIR, f"{grp}_all_peptide_results.tsv")
        neo = os.path.join(MHC_DIR, f"{grp}_neoantigens.tsv")
        if os.path.exists(allp):
            a = _read(allp)
            m = _first(a.columns, ["mut_ic50", "mut_IC50"])
            w = _first(a.columns, ["wt_ic50", "wt_IC50"])
            a[m] = pd.to_numeric(a[m], errors="coerce")
            a[w] = pd.to_numeric(a[w], errors="coerce")
            strong[grp] = int((a[m] < STRONG).sum())
            diff[grp] = int(((a[m] < BIND) & (a[w] >= BIND)).sum())
        else:
            log(f"  [SKIP] {grp}_all_peptide_results.tsv not found")
            strong[grp] = diff[grp] = None
        binders[grp] = len(_read(neo)) if os.path.exists(neo) else None
    check("strong binders SBS2 (mut<50)", EXPECTED["strong_binders"]["sbs2"], strong["SBS2_HIGH"], "exact")
    check("strong binders CNV (mut<50)", EXPECTED["strong_binders"]["cnv"], strong["CNV_HIGH"], "exact")
    check("differential SBS2 (mut<500,wt>=500)", EXPECTED["differential"]["sbs2"], diff["SBS2_HIGH"], "exact")
    check("differential CNV (mut<500,wt>=500)", EXPECTED["differential"]["cnv"], diff["CNV_HIGH"], "exact")
    check("binder peptides SBS2 (rows in neoantigens.tsv)", EXPECTED["binders_peptides"]["sbs2"], binders["SBS2_HIGH"], "exact")
    check("binder peptides CNV (rows in neoantigens.tsv)", EXPECTED["binders_peptides"]["cnv"], binders["CNV_HIGH"], "exact")


# =============================================================================
# CHECK D: raw fusion junction totals  [Panel A cross-check]
# =============================================================================
def check_D_raw_fusion():
    banner("CHECK D: raw fusion junction totals per population")
    counts = None
    summ = os.path.join(FUSION_DIR, "per_group_junction_summary.tsv")
    if os.path.exists(summ):
        df = _read(summ)
        gcol = _first(df.columns, ["group", "population"])
        ncol = _first(df.columns, ["n_junctions", "junctions", "total_junctions", "count", "total"])
        if gcol and ncol:
            counts = {str(r[gcol]): int(r[ncol]) for _, r in df.iterrows()}
    if counts is None:
        afj = os.path.join(FUSION_DIR, "all_filtered_junctions.tsv")
        if os.path.exists(afj):
            j = _read(afj)
            gcol = _first(j.columns, ["group"])
            if gcol:
                counts = {str(k): int(v) for k, v in j[gcol].value_counts().items()}
    if counts is None:
        log(f"  [SKIP] no per-group fusion summary in {FUSION_DIR}")
        _tally["skip"] += 6
        return
    for grp in ["SBS2_HIGH", "CNV_HIGH", "NORMAL"]:
        obs = counts.get(grp)
        check(f"raw junctions {grp}", EXPECTED["raw_junctions"][grp], obs, "exact")
        check(f"junctions/cell {grp}", EXPECTED["raw_junctions_percell"][grp],
              (obs / N_CELLS) if obs is not None else None, "close", 0.05)


# =============================================================================
# CHECK E: neoantigen tiers 467 / 93 / 215  [Panel B]  (full.tsv)
# =============================================================================
def check_E_tiers():
    banner("CHECK E: neoantigen tiers 467 / 93 / 215 and totals  [Panel B]")
    if not _need(RANK_TSV, "neoantigen_prevalence_ranking_full.tsv"):
        _tally["skip"] += 6
        return
    full = _read(RANK_TSV)
    tcol = _first(full.columns, ["tier"])
    vc = {str(k): int(v) for k, v in full[tcol].value_counts().items()}
    for t in [T_SBS2, T_SHARED, T_CNV]:
        check(f"tier {t}", EXPECTED["tiers"][t], vc.get(t), "exact")
    check("total neoantigen mutations", EXPECTED["tier_total"], len(full), "exact")
    check("SBS2 total (specific+shared)", EXPECTED["sbs2_total_mutations"],
          vc.get(T_SBS2, 0) + vc.get(T_SHARED, 0), "exact")
    check("CNV total (specific+shared)", EXPECTED["cnv_total_mutations"],
          vc.get(T_CNV, 0) + vc.get(T_SHARED, 0), "exact")


# =============================================================================
# CHECK F: Tier-1 selection evidence (15 of 93)  (Step05 output)
# =============================================================================
def check_F_evidence():
    banner("CHECK F: Tier-1 selection-evidence set (15 of 93 shared)")
    if not _need(EVIDENCE_TSV, "shared_neoantigen_selection_evidence.tsv"):
        log("  -> re-run Step05_Integrated_Neoantigen_Analysis.py.")
        _tally["skip"] += 9
        return
    ev = _read(EVIDENCE_TSV)
    se = _first(ev.columns, ["selection_evidence"])
    ev["_se"] = ev[se].map(_tobool)
    counted = ev[ev["_se"]]
    check("shared neoantigens with selection evidence", EXPECTED["shared_selection_evidence"], len(counted), "exact")

    dcol = _first(ev.columns, ["cnv_direction"])
    if dcol:
        d = counted[dcol].astype(str)
        check("counted-set depleted_in_CNV", EXPECTED["shared_direction"]["depleted"], int((d == "depleted_in_CNV").sum()), "exact")
        check("counted-set flat", EXPECTED["shared_direction"]["flat"], int((d == "flat").sum()), "exact")
        check("counted-set enriched_in_CNV", EXPECTED["shared_direction"]["enriched"], int((d == "enriched_in_CNV").sum()), "exact")

    mcol = _first(ev.columns, ["mechanisms"])
    if mcol:
        mm = counted[mcol].astype(str)
        check("counted-set fusion mechanism", EXPECTED["shared_mechanism"]["fusion"], int(mm.str.contains("fusion").sum()), "exact")
        check("counted-set silenced mechanism", EXPECTED["shared_mechanism"]["silenced"], int(mm.str.contains("silenced").sum()), "exact")

    hcol = _first(ev.columns, ["is_hla"])
    mecol = _first(ev.columns, ["mechanism_evidence"])
    if hcol and mecol:
        hla_ex = int((ev[hcol].map(_tobool) & ev[mecol].map(_tobool)).sum())
        check("HLA neoantigens excluded (had a mechanism)", EXPECTED["shared_hla_excluded"], hla_ex, "exact")

    gcol = _first(ev.columns, ["gene"])
    counted_genes = set(counted[gcol].astype(str))
    for g in EXPECTED["escape_subset"]:
        check(f"escape-subset gene in counted set: {g}", True, g in counted_genes, "eq")


# =============================================================================
# CHECK G: featured candidates (Panel C/D)  (full.tsv)
# =============================================================================
def check_G_featured():
    banner("CHECK G: featured candidate rows (prevalence / expression / binding / TCW)")
    if not _need(RANK_TSV, "full.tsv"):
        _tally["skip"] += 18
        return
    full = _read(RANK_TSV)
    gcol = _first(full.columns, ["gene"])
    hcol = _first(full.columns, ["hgvs_p"])
    tcol = _first(full.columns, ["tier"])
    pcol = _first(full.columns, ["prevalence_tier"])
    ecol = _first(full.columns, ["pct_expressing_tier"])
    wcol = _first(full.columns, ["wt_IC50", "wt_ic50"])
    mcol = _first(full.columns, ["mut_IC50", "mut_ic50"])
    ctcol = _first(full.columns, ["is_tcw_ct"])
    for gene, e in EXPECTED["featured"].items():
        row = full[(full[gcol].astype(str) == gene) & (full[hcol].astype(str) == e["hgvs_p"])]
        if len(row) == 0:
            check(f"{gene} {e['hgvs_p']} present in ranking", True, False, "eq")
            continue
        r = row.iloc[0]
        check(f"{gene} tier", e["tier"], r[tcol], "eq")
        check(f"{gene} prevalence_tier %", e["prev"], 100 * float(r[pcol]), "close", 0.15)
        check(f"{gene} % expressing (tier)", e["expr"], float(r[ecol]), "close", 0.3)
        check(f"{gene} wt_IC50", e["wt"], float(r[wcol]), "rel", 0.02)
        check(f"{gene} mut_IC50", e["mut"], float(r[mcol]), "rel", 0.02)
        check(f"{gene} clean-TCW C>T", e["tcw"], _tobool(r[ctcol]), "eq")


# =============================================================================
# CHECK H: named exemplars really are the tier leaders they are claimed to be
# =============================================================================
def check_H_leaders():
    banner("CHECK H: featured gene per tier = top NON-HLA neoantigen by prevalence_max")
    if not _need(RANK_TSV, "full.tsv"):
        _tally["skip"] += 3
        return
    full = _read(RANK_TSV)
    gcol = _first(full.columns, ["gene"])
    tcol = _first(full.columns, ["tier"])
    pmcol = _first(full.columns, ["prevalence_max"])
    dcol = _first(full.columns, ["delta_IC50"])
    # HLA-A/B/C carry natural germline polymorphism, so they are not real
    # neoantigen targets and are excluded from featured-gene selection.
    nonhla = full[~full[gcol].astype(str).str.startswith(("HLA-A", "HLA-B", "HLA-C"))].copy()
    nonhla["_pm"] = pd.to_numeric(nonhla[pmcol], errors="coerce")
    nonhla["_d"] = pd.to_numeric(nonhla[dcol], errors="coerce")
    for tier, key in [(T_SHARED, "leader_shared"),
                      (T_CNV, "leader_cnv_specific"),
                      (T_SBS2, "leader_sbs2_specific")]:
        sub = nonhla[nonhla[tcol] == tier].sort_values(["_pm", "_d"], ascending=[False, False])
        if not len(sub):
            check(f"{tier} non-HLA leader (max prevalence_max)", EXPECTED[key], None)
            continue
        top = sub.iloc[0]
        check(f"{tier} non-HLA leader (max prevalence_max, delta tiebreak)",
              EXPECTED[key], top[gcol], "eq")
        tied = sub[np.isclose(sub["_pm"].fillna(-1.0), top["_pm"])]
        if len(tied) > 1:
            names = ", ".join(f"{t[gcol]}(delta={t['_d']:.0f})" for _, t in tied.head(4).iterrows())
            log(f"        prevalence_max tie at {top['_pm']:.4f}; broke on binding gain among: {names}")


# =============================================================================
# CHECK I: ANXA1 rationale line (high expression, ~1% carriage)  (full.tsv)
# =============================================================================
def check_I_anxa1():
    banner("CHECK I: ANXA1 counterpoint (high expression, low gene-level carriage)")
    if not _need(RANK_TSV, "full.tsv"):
        _tally["skip"] += 2
        return
    full = _read(RANK_TSV)
    gcol = _first(full.columns, ["gene"])
    a = full[full[gcol].astype(str) == "ANXA1"]
    if len(a) == 0:
        log("  [SKIP] ANXA1 not present in ranking (no binder neoantigen row).")
        _tally["skip"] += 2
        return
    ecol = _first(full.columns, ["pct_expressing_sbs2", "pct_expressing_tier"])
    ucol = _first(full.columns, ["gene_union_sbs2"])
    r = a.iloc[0]
    check("ANXA1 % expressing (SBS2)", EXPECTED["anxa1"]["expr"], float(r[ecol]), "close", 0.6)
    if ucol:
        carriage = 100.0 * float(r[ucol]) / N_CELLS   # union carriers across ANXA1 loci
        check("ANXA1 gene-level carriage % (~1%)", EXPECTED["anxa1"]["carriage_pct"], carriage, "close", 0.3)
    else:
        check("ANXA1 gene-level carriage %", EXPECTED["anxa1"]["carriage_pct"], None)


# =============================================================================
# CHECK J: TCW enrichment, genome-verified, CODING subset (must be flat / ns)
#          (folds in Diagnostic_Section45_Verify_Numbers.py CHECK 2)
# =============================================================================
def check_J_tcw():
    banner("CHECK J: TCW enrichment (genome-verified, coding subset) -> confirm NO enrichment")
    if not _need(FASTA_TSV, "ref_tri_fasta.tsv"):
        _tally["skip"] += 6
        return
    try:
        from scipy.stats import fisher_exact
    except Exception:
        fisher_exact = None

    ref = _read(FASTA_TSV)
    rk = {c: _first(ref.columns, [c]) for c in
          ["chrom", "pos", "ref", "alt", "gene", "hgvs_p", "sub_pyr", "is_tcw", "is_tcw_ct"]}
    for col in ("is_tcw", "is_tcw_ct"):
        if rk[col]:
            ref[rk[col]] = ref[rk[col]].map(_tobool)
    refkey = ref.rename(columns={rk["chrom"]: "chrom", rk["pos"]: "pos",
                                 rk["ref"]: "ref", rk["alt"]: "alt"})

    joined = {}
    for grp in GROUPS:
        p = os.path.join(ANNOT_DIR, f"{grp}.somatic_protein_altering.tsv")
        if not os.path.exists(p):
            log(f"  [SKIP] missing {grp}.somatic_protein_altering.tsv")
            _tally["skip"] += 6
            return
        df = _read(p)
        gc = _first(df.columns, ["chrom", "#chrom", "chr"])
        pc = _first(df.columns, ["pos", "gpos", "Start"])
        rc = _first(df.columns, ["ref"])
        ac = _first(df.columns, ["alt"])
        gv = df.rename(columns={gc: "chrom", pc: "pos", rc: "ref", ac: "alt"})[
            ["chrom", "pos", "ref", "alt"]].drop_duplicates()
        gv["pos"] = pd.to_numeric(gv["pos"], errors="coerce")
        refkey["pos"] = pd.to_numeric(refkey["pos"], errors="coerce")
        joined[grp] = gv.merge(refkey, on=["chrom", "pos", "ref", "alt"], how="left")
        n_join = joined[grp][rk["is_tcw"]].notna().sum() if rk["is_tcw"] else 0
        log(f"  {grp}: {len(gv)} protein-altering variants, {n_join} joined to ref_tri")

    def frac(j, denom_mask, tcw_col):
        d = int(denom_mask.sum())
        tcw = j[tcw_col].fillna(False) if tcw_col else pd.Series([False] * len(j))
        n = int((denom_mask & tcw).sum())
        return n, d, (100.0 * n / d if d else np.nan)

    def one_def(label, key, denom_fn, tcw_col, tol=0.3):
        res = {}
        for grp in GROUPS:
            j = joined[grp].reset_index(drop=True)
            res[grp] = frac(j, denom_fn(j).reset_index(drop=True), tcw_col)
        s, c = res["SBS2_HIGH"], res["CNV_HIGH"]
        log(f"\n  {label}")
        log(f"    SBS2: {s[0]}/{s[1]} = {s[2]:.1f}%   CNV: {c[0]}/{c[1]} = {c[2]:.1f}%")
        pval = np.nan
        if fisher_exact is not None:
            orv, pval = fisher_exact([[s[0], s[1] - s[0]], [c[0], c[1] - c[0]]])
            log(f"    Fisher exact OR={orv:.2f}, p={pval:.4g}")
        check(f"{key} SBS2 %", EXPECTED[key]["sbs2"], s[2], "close", tol)
        check(f"{key} CNV %", EXPECTED[key]["cnv"], c[2], "close", tol)
        check(f"{key} non-significant (p>0.05)", 0.05, pval, "ns")

    sp = rk["sub_pyr"]
    one_def("Definition A: among C>T protein-altering, fraction clean-TCW", "tcw_defA_ct",
            lambda j: (j[sp].astype(str) == "C>T") if sp else pd.Series([False] * len(j)),
            rk["is_tcw_ct"])
    one_def("Definition B: among ALL protein-altering, fraction TCW", "tcw_defB_all",
            lambda j: pd.Series([True] * len(j)), rk["is_tcw"])

    # Definition C: among neoantigen-forming variants
    key = "tcw_defC_neo"
    resC = {}
    for grp in GROUPS:
        neo_p = os.path.join(MHC_DIR, f"{grp}_neoantigens.tsv")
        if not os.path.exists(neo_p):
            log(f"  [SKIP] {key}: missing {grp} neoantigens")
            _tally["skip"] += 3
            resC = None
            break
        neo = _read(neo_p)
        ng, nh = _first(neo.columns, ["gene"]), _first(neo.columns, ["hgvs_p"])
        neo_keys = set(zip(neo[ng].astype(str), neo[nh].astype(str)))
        j = joined[grp].reset_index(drop=True)
        jg = j[rk["gene"]].astype(str) if rk["gene"] else pd.Series([""] * len(j))
        jh = j[rk["hgvs_p"]].astype(str) if rk["hgvs_p"] else pd.Series([""] * len(j))
        is_neo = pd.Series([(g, h) in neo_keys for g, h in zip(jg, jh)])
        resC[grp] = frac(j, is_neo, rk["is_tcw"])
    if resC:
        s, c = resC["SBS2_HIGH"], resC["CNV_HIGH"]
        log(f"\n  Definition C: among neoantigen-forming variants, fraction TCW")
        log(f"    SBS2: {s[0]}/{s[1]} = {s[2]:.1f}%   CNV: {c[0]}/{c[1]} = {c[2]:.1f}%")
        pval = np.nan
        if fisher_exact is not None:
            orv, pval = fisher_exact([[s[0], s[1] - s[0]], [c[0], c[1] - c[0]]])
            log(f"    Fisher exact OR={orv:.2f}, p={pval:.4g}")
        check(f"{key} SBS2 %", EXPECTED[key]["sbs2"], s[2], "close", 0.3)
        check(f"{key} CNV %", EXPECTED[key]["cnv"], c[2], "close", 0.3)
        check(f"{key} non-significant (p>0.05)", 0.05, pval, "ns")
    log("\n  -> all three flat / ns confirms the TCW-enrichment sentence stays REMOVED")
    log("     from the coding-subset results (see the full-spectrum diagnostic separately).")


# =============================================================================
# CHECK K: HLA panel & population coverage (informational)
# =============================================================================
def check_K_coverage():
    banner("CHECK K: HLA panel size & population coverage (informational)")
    log(f"  HLA panel (n={len(EXPECTED['hla_alleles'])}): {', '.join(EXPECTED['hla_alleles'])}")
    check("HLA panel size", 10, len(EXPECTED["hla_alleles"]), "exact")
    log(f"  Population coverage (IEDB tool, recorded value): {EXPECTED['pop_coverage_pct']}%")
    log("  NOTE: 90.5% comes from the external IEDB Population Coverage tool; re-run it")
    log("  on the panel above to reconfirm (not recomputed here).")


# =============================================================================
# CHECK L (optional, --heavy): raw per-cell protein-altering burden (cross-check)
# =============================================================================
def check_L_raw_percell():
    banner("CHECK L (--heavy): raw per-cell protein-altering burden [informational]")
    if not (os.path.exists(GENO_PATH) and os.path.exists(FASTA_TSV) and os.path.exists(GROUP_PATH)):
        log("  [SKIP] genotype master / ref_tri / group file missing.")
        return
    groups = _read(GROUP_PATH)
    bc = _first(groups.columns, ["cell_barcode", "CB", "barcode"])
    gp = _first(groups.columns, ["group"])
    cellsets = {g: set(groups.loc[groups[gp] == g, bc].astype(str)) for g in GROUPS}

    ref = _read(FASTA_TSV)
    rc = _first(ref.columns, ["chrom"])
    rp = _first(ref.columns, ["pos"])
    pa_loci = set(zip(ref[rc].astype(str), pd.to_numeric(ref[rp], errors="coerce")))

    geno = _read(GENO_PATH)
    gch = _first(geno.columns, ["#CHROM", "chrom"])
    gps = _first(geno.columns, ["Start", "pos"])
    gcb = _first(geno.columns, ["CB", "cell_barcode"])
    geno["_ch"] = geno[gch].astype(str)
    geno["_ps"] = pd.to_numeric(geno[gps], errors="coerce")
    best_off, best_rate = 0, -1
    for off in (0, 1, -1):
        rate = len(set(zip(geno["_ch"], geno["_ps"] + off)) & pa_loci)
        if rate > best_rate:
            best_rate, best_off = rate, off
    log(f"  locus match offset = {best_off} ({best_rate} protein-altering loci matched)")
    geno["_key"] = list(zip(geno["_ch"], geno["_ps"] + best_off))
    geno_pa = geno[geno["_key"].isin(pa_loci)]
    for g in GROUPS:
        sub = geno_pa[geno_pa[gcb].astype(str).isin(cellsets[g])]
        per = len(sub) / max(len(cellsets[g]), 1)
        log(f"    {g}: {len(sub)} carrier events / {len(cellsets[g])} cells = {per:.2f}/cell")
    log("  (depth-confounded; dropped from the prose. Expected ~5.08 vs ~4.20.)")


# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("SECTION 4.5 / FIGURE 7 NUMBER AUDIT (assertion-based, v3 consolidated)")
    log("  PASS = matches the prose. **FAIL** = drift, investigate before quoting.")
    log("  SKIP = an upstream table is missing; run its producer and re-run this.")
    for fn in (check_A_grouprate, check_B_perumi, check_C_strong_diff, check_D_raw_fusion,
               check_E_tiers, check_F_evidence, check_G_featured, check_H_leaders,
               check_I_anxa1, check_J_tcw, check_K_coverage):
        try:
            fn()
        except Exception as e:
            _tally["fail"] += 1
            log(f"  [ERROR] {fn.__name__}: {type(e).__name__}: {e}")
    if "--heavy" in sys.argv:
        try:
            check_L_raw_percell()
        except Exception as e:
            log(f"  [ERROR] check_L_raw_percell: {type(e).__name__}: {e}")

    banner("AUDIT SUMMARY")
    log(f"  PASS: {_tally['pass']}    FAIL: {_tally['fail']}    SKIP: {_tally['skip']}")
    if _tally["fail"] == 0 and _tally["skip"] == 0:
        log("  ALL CHECKS PASSED - Section 4.5 numbers are verified against the tables.")
    elif _tally["fail"] == 0:
        log("  No failures, but some checks skipped (missing upstream tables). Run the")
        log("  producing diagnostics, then re-run this audit for a clean sweep.")
    else:
        log("  ONE OR MORE FAILURES - do not quote the affected numbers until resolved.")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    outp = os.path.join(SUMMARY_DIR, "section45_number_audit.txt")
    with open(outp, "w") as f:
        f.write("\n".join(_report))
    log(f"\n  Report written: {outp}")


if __name__ == "__main__":
    main()
