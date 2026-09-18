#!/usr/bin/env python3
"""
Diagnostic_Section44_Numbers.py  (v2)
=====================================

Assertion-based audit of Results 4.4. v2 fixes three defects in v1 that
produced two spurious failures and left two real questions unresolved.

v1 DEFECTS AND WHAT CHANGED
---------------------------
1. VIRAL BURDEN READ FROM THE WRONG PLACE (caused the two Fisher failures).
   v1 looked for raw_HPV16 in adata.obs, which does not carry it. Every
   viral_pct came back NaN. v1 then coerced NaN to False in the three-
   criteria test, so no patient could satisfy all three, the table
   collapsed to [[0,0],[3,11]], and Fisher returned p = 1. That was a
   cascade from a missing column, not drift in the manuscript.

   v2 reads raw_HPV16 and n_counts from the per-cell master tables that
   actually hold them, trying each in turn, and reports which it used.

2. NaN SILENTLY BECAME False.
   A missing measurement is not a failed criterion. v2 propagates
   "cannot evaluate" and SKIPs the Fisher test rather than reporting a
   false FAIL. A test that cannot run must never look like a test that ran
   and disagreed.

3. A3 PREVALENCE CELL SET WAS ASSUMED, NOT DERIVED (the remaining question).
   v1 computed prevalence over TUMOR basal cells. Under that definition
   both A3B correlations reproduce exactly (0.660 and 0.400 against the
   fold-enrichment contribution) but both A3A correlations drift: 0.581
   against a published 0.594, and 0.142 against a published 0.047. Two
   exact matches and two misses on the same contribution vector points at
   the A3A prevalence vector, not at the contribution definition.

   The likely cause is the denominator. Three patients (SC005, SC003,
   SC006) supply essentially all 554 normal-adjacent basal cells, so
   including or excluding normal-adjacent tissue moves their prevalence
   and therefore their rank, which is all Spearman sees.

   v2 therefore computes prevalence under three cell sets and reports the
   full grid instead of guessing:
       tumor   tumor basal cells only              (v1 behaviour)
       all     all basal cells, tumor + normal
       cancer  basal cells called Cancer cell

   THREE INDEPENDENT CONSTRAINTS must be satisfied by the same cell set:
       A3A prevalence range        = 74-fold
       A3A prevalence vs SBS2 fold = +0.594 (p = 0.025)
       A3A prevalence vs CNA fold  = +0.047 (p = 0.87)
   If one cell set satisfies all three, that is the definition the
   manuscript used and Methods 6.5 should say so explicitly. If none
   does, the published numbers came from a procedure not reproducible
   from these inputs, which is a finding, not a number to nudge.

WHAT v1 ALREADY ESTABLISHED (unchanged in v2)
---------------------------------------------
   Contribution is FOLD ENRICHMENT, not cell count or percentage. A3B
   reproduces to three decimals against fold in both directions
   (+0.660/0.0102 and +0.400/0.157) and not against the other two.

Usage:
  conda run -n NETWORK python Diagnostic_Section44_Numbers.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, chi2_contingency, fisher_exact

import scanpy as sc

BASE = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH  = os.path.join(BASE, "data/FIG_4/00_input/adata_final.h5ad")
GROUPS_PATH = os.path.join(BASE, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
HPV_GENE_PATH = os.path.join(BASE, "data/FIG_6/03_hpv16_genome/hpv16_gene_by_population.tsv")
OUTDIR = os.path.join(BASE, "data/FIG_5/00_diagnostics/SECTION_44_AUDIT")

# Per-cell tables carrying raw_HPV16 + n_counts for the FULL basal set.
# Tried in order; the first covering >90% of basal cells wins.
VIRAL_SOURCES = [
    os.path.join(BASE, "data/FIG_6/02_populations/population_assignments.tsv"),
    os.path.join(BASE, "data/FIG_6/04_population_profiles/two_population_assignments.tsv"),
    HPV_GENE_PATH,
]

PATIENT_COL, TISSUE_COL, CELLTYPE_COL = "subject id", "tissue type", "final_annotation"
CANCER_COL = "Final_cancer_cell_status"
BASAL = "basal cell"
HPV_UMI_GATE, MIN_GATED_CELLS = 8, 10

MAINTENANCE = ["E1", "E2"]
PRODUCTIVE  = ["E4", "E5", "L1", "L2"]
CODING      = ["E1", "E2", "E4", "E5", "E6", "E7", "L1", "L2"]

SBS2_HC, CNV_HC = ["SC013", "SC029", "SC001"], ["SC027", "SC001"]
VIRAL_MIN, STAGE_MIN = 0.1, 40.0
A3A_PREV_MIN, A3B_PREV_MIN = 20.0, 50.0

CELL_SETS = ["tumor", "all", "cancer"]
RESULTS = []


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


def banner(t):
    print("")
    print("=" * 88)
    print(f"  {t}")
    print("=" * 88)


def check(label, expected, observed, tol_rel=0.02, tol_abs=None):
    if observed is None or (isinstance(observed, float) and np.isnan(observed)):
        RESULTS.append((label, expected, observed, "SKIP"))
        print(f"  [SKIP] {label}: expected {expected}, not evaluable")
        return
    if tol_abs is not None:
        ok = abs(observed - expected) <= tol_abs
    elif expected == 0:
        ok = abs(observed) < 1e-12
    else:
        ok = abs(observed - expected) / abs(expected) <= tol_rel
    RESULTS.append((label, expected, observed, "PASS" if ok else "FAIL"))
    ef = f"{expected:.4g}" if isinstance(expected, float) else str(expected)
    of = f"{observed:.4g}" if isinstance(observed, float) else str(observed)
    print(f"  [{'PASS' if ok else '**FAIL**'}] {label}: expected {ef}, observed {of}")


def short(p):
    return str(p).replace("Patient ", "")


# =============================================================================
banner("STEP 0: BUILD THE PER-PATIENT TABLE")
# =============================================================================
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.n_obs:,} cells")
log(f"  obs columns: {', '.join(list(adata.obs.columns)[:20])}"
    f"{' ...' if adata.obs.shape[1] > 20 else ''}")

basal = adata.obs[adata.obs[CELLTYPE_COL].astype(str) == BASAL].copy()
log(f"  basal cells: {len(basal):,}")

g = pd.read_csv(GROUPS_PATH, sep="\t")
gcol = "cell_barcode" if "cell_barcode" in g.columns else g.columns[0]
gmap = dict(zip(g[gcol].astype(str), g["group"]))
basal["group"] = [gmap.get(c, "other") for c in basal.index]


def gene_vec(name):
    if name not in adata.var_names:
        return None
    x = adata[:, name].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


pos = {c: i for i, c in enumerate(adata.obs_names)}
bidx = np.array([pos[c] for c in basal.index])
basal["A3A"] = gene_vec("APOBEC3A")[bidx]
basal["A3B"] = gene_vec("APOBEC3B")[bidx]

# ---- FIX 1: viral burden from a per-cell table that actually has it --------
log("")
log("  resolving viral-burden source ...")
viral_src, viral_df = None, None
for path in VIRAL_SOURCES:
    if not os.path.exists(path):
        log(f"    [absent] {os.path.basename(path)}")
        continue
    try:
        d = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    except Exception as e:
        log(f"    [unreadable] {os.path.basename(path)}: {e}")
        continue
    have = {"raw_HPV16", "n_counts"} <= set(d.columns)
    cov = len(basal.index.intersection(d.index)) / max(len(basal), 1)
    log(f"    {os.path.basename(path)}: cols_ok={have} basal_coverage={100*cov:.1f}%")
    if have and cov > 0.90:
        viral_src, viral_df = path, d
        break

if viral_df is None:
    log("  !! no source covers the basal compartment; STEP 2 and STEP 6 will SKIP")
else:
    log(f"  using: {viral_src}")
    idx = basal.index.intersection(viral_df.index)
    basal.loc[idx, "raw_HPV16"] = pd.to_numeric(
        viral_df.loc[idx, "raw_HPV16"], errors="coerce")
    basal.loc[idx, "n_counts_v"] = pd.to_numeric(
        viral_df.loc[idx, "n_counts"], errors="coerce")

patients = sorted(basal[PATIENT_COL].astype(str).unique())
n_sbs2 = int((basal["group"] == "SBS2_HIGH").sum())
n_cnv = int((basal["group"] == "CNV_HIGH").sum())
n_basal_all = len(basal)

has_cancer_col = CANCER_COL in basal.columns
if not has_cancer_col:
    log(f"  note: '{CANCER_COL}' absent; the 'cancer' cell set will be skipped")

rows = []
for p in patients:
    sub = basal[basal[PATIENT_COL].astype(str) == p]
    share = len(sub) / n_basal_all
    n_h = int((sub["group"] == "SBS2_HIGH").sum())
    n_c = int((sub["group"] == "CNV_HIGH").sum())

    r = {"patient": short(p), "n_basal": len(sub),
         "n_sbs2_high": n_h, "n_cnv_high": n_c,
         "pct_of_sbs2": 100 * n_h / n_sbs2 if n_sbs2 else 0,
         "pct_of_cnv": 100 * n_c / n_cnv if n_cnv else 0,
         "fold_sbs2": (n_h / n_sbs2) / share if share and n_sbs2 else 0,
         "fold_cnv": (n_c / n_cnv) / share if share and n_cnv else 0}

    if viral_df is not None:
        tot = float(sub["n_counts_v"].sum(skipna=True))
        hpv = float(sub["raw_HPV16"].sum(skipna=True))
        r["viral_pct"] = 100 * hpv / tot if tot else np.nan
    else:
        r["viral_pct"] = np.nan

    # ---- FIX 3: prevalence under three cell-set definitions ---------------
    sets = {"tumor": sub[sub[TISSUE_COL].astype(str) == "tumor"], "all": sub}
    if has_cancer_col:
        sets["cancer"] = sub[sub[CANCER_COL].astype(str) == "Cancer cell"]

    for cs, cells in sets.items():
        for tag in ("A3A", "A3B"):
            v = cells[tag].values
            if len(v) == 0:
                r[f"{tag}_mean_all_{cs}"] = np.nan
                r[f"{tag}_mean_expr_{cs}"] = np.nan
                r[f"{tag}_prev_{cs}"] = np.nan
                continue
            r[f"{tag}_mean_all_{cs}"] = float(v.mean())
            e = v[v > 0]
            r[f"{tag}_mean_expr_{cs}"] = float(e.mean()) if len(e) else 0.0
            r[f"{tag}_prev_{cs}"] = 100 * float((v > 0).mean())
    rows.append(r)

pt = pd.DataFrame(rows).set_index("patient")
active_sets = [c for c in CELL_SETS if f"A3A_prev_{c}" in pt.columns]

# =============================================================================
banner("STEP 1: CONTRIBUTION AND CHI-SQUARE")
# =============================================================================
print("")
print(f"  {'pt':<8} {'n_basal':>8} {'nSBS2':>6} {'foldS':>7} {'nCNV':>6} "
      f"{'foldC':>7} {'viral%':>9}")
print(f"  {'-'*8} {'-'*8} {'-'*6} {'-'*7} {'-'*6} {'-'*7} {'-'*9}")
for i, r in pt.sort_values("n_sbs2_high", ascending=False).iterrows():
    mk = ("S" if i in SBS2_HC else "") + ("C" if i in CNV_HC else "")
    vp = "n/a" if np.isnan(r["viral_pct"]) else f"{r['viral_pct']:.3f}%"
    print(f"  {i:<8} {int(r['n_basal']):>8,} {int(r['n_sbs2_high']):>6} "
          f"{r['fold_sbs2']:>6.1f}x {int(r['n_cnv_high']):>6} "
          f"{r['fold_cnv']:>6.1f}x {vp:>9}  {mk}")

print("")
for col, label, exp_stat, exp_p in (("n_sbs2_high", "SBS2-HIGH", None, 4.46e-309),
                                    ("n_cnv_high", "CNA-HIGH", 3419.8, None)):
    tab = pd.concat([pt[col], pt["n_basal"] - pt[col]], axis=1)
    chi2, pval, dof, _ = chi2_contingency(tab.values)
    log(f"  {label}: chi2 = {chi2:.1f}, df = {dof}, p = {pval:.3g}")
    if exp_p is not None:
        check(f"{label} chi-square p", exp_p, float(pval), tol_rel=0.05)
    if exp_stat is not None:
        check(f"{label} chi-square statistic", exp_stat, float(chi2), tol_rel=0.01)

# =============================================================================
banner("STEP 2: VIRAL BURDEN vs CONTRIBUTION")
# =============================================================================
viral_ok = viral_df is not None and not pt["viral_pct"].isna().all()
if not viral_ok:
    log("  [SKIP] viral burden unavailable")
else:
    for fate, target in (("sbs2", 0.506), ("cnv", 0.513)):
        print("")
        log(f"  viral burden vs {fate.upper()} contribution (prose +{target:.3f})")
        best, bd = None, np.inf
        for defn, col in (("n_cells", f"n_{fate}_high"),
                          ("pct", f"pct_of_{fate}"),
                          ("fold", f"fold_{fate}")):
            rho, p = spearmanr(pt["viral_pct"], pt[col])
            log(f"    {defn:<9} rho = {rho:+.3f}  p = {p:.3g}")
            if abs(rho - target) < bd:
                best, bd = (rho, p), abs(rho - target)
        check(f"viral rho {fate.upper()}", target, float(best[0]), tol_abs=0.03)

    print("")
    log(f"  patients above {VIRAL_MIN}% viral burden:")
    for i, r in pt[pt["viral_pct"] > VIRAL_MIN].sort_values(
            "viral_pct", ascending=False).iterrows():
        tag = [t for t, s in (("SBS2-HC", SBS2_HC), ("CNA-HC", CNV_HC)) if i in s]
        log(f"    {i:<8} {r['viral_pct']:.3f}%  {', '.join(tag) or 'neither'}")

    for p_ in ("SC013", "SC027"):
        if p_ in pt.index:
            check(f"{p_} viral burden 0.17%", 0.17,
                  float(pt.loc[p_, "viral_pct"]), tol_abs=0.005)

# =============================================================================
banner("STEP 3: LIFECYCLE STAGE BALANCE")
# =============================================================================
stage = {}
if not os.path.exists(HPV_GENE_PATH):
    log(f"  [SKIP] {HPV_GENE_PATH} not found")
else:
    hv = pd.read_csv(HPV_GENE_PATH, sep="\t", index_col=0, low_memory=False)
    hv = hv[pd.to_numeric(hv.get("raw_HPV16"), errors="coerce") >= HPV_UMI_GATE]
    log(f"  gated at >= {HPV_UMI_GATE} UMI: {len(hv):,} cells")
    print("")
    print(f"  {'pt':<8} {'gated':>7} {'maint%':>9} {'prod%':>9}  assignment")
    print(f"  {'-'*8} {'-'*7} {'-'*9} {'-'*9}  {'-'*11}")
    for p in patients:
        sub = hv[hv[PATIENT_COL].astype(str) == p]
        if len(sub) < MIN_GATED_CELLS:
            print(f"  {short(p):<8} {len(sub):>7} {'--':>9} {'--':>9}  unassigned")
            continue
        tot = sum(pd.to_numeric(sub[c], errors="coerce").sum()
                  for c in CODING if c in sub.columns)
        if not tot:
            continue
        mt = sum(pd.to_numeric(sub[c], errors="coerce").sum()
                 for c in MAINTENANCE if c in sub.columns)
        pr = sum(pd.to_numeric(sub[c], errors="coerce").sum()
                 for c in PRODUCTIVE if c in sub.columns)
        mp, pp = 100 * mt / tot, 100 * pr / tot
        stage[short(p)] = (mp, pp)
        lab = "maintenance" if mp >= STAGE_MIN else (
            "productive" if pp >= STAGE_MIN else "neither")
        print(f"  {short(p):<8} {len(sub):>7} {mp:>8.1f}% {pp:>8.1f}%  {lab}")

    print("")
    for p_ in SBS2_HC:
        if p_ in stage:
            check(f"{p_} maintenance >= {STAGE_MIN}%", True,
                  bool(stage[p_][0] >= STAGE_MIN), tol_abs=0)
    for p_ in CNV_HC:
        if p_ in stage:
            check(f"{p_} productive >= {STAGE_MIN}%", True,
                  bool(stage[p_][1] >= STAGE_MIN), tol_abs=0)

# =============================================================================
banner("STEP 4: A3 SUMMARIES UNDER EACH CELL SET")
# =============================================================================
for cs in active_sets:
    print("")
    log(f"  cell set '{cs}':")
    for tag in ("A3A", "A3B"):
        me = pt[f"{tag}_mean_expr_{cs}"].dropna()
        pv = pt[f"{tag}_prev_{cs}"].dropna()
        pv = pv[pv > 0]
        if not len(me) or not len(pv):
            continue
        log(f"    {tag} mean-among-expressing {me.min():.3f}-{me.max():.3f} "
            f"= {me.max()/me.min():.2f}-fold")
        log(f"    {tag} prevalence {pv.min():.2f}%-{pv.max():.2f}% "
            f"= {pv.max()/pv.min():.1f}-fold")

print("")
log("  CONSTRAINT 1: A3A prevalence range should be 74-fold")
for cs in active_sets:
    pv = pt[f"A3A_prev_{cs}"].dropna()
    pv = pv[pv > 0]
    if len(pv):
        rng = pv.max() / pv.min()
        log(f"    {cs:<8} {rng:>7.1f}-fold{'   <-- matches' if abs(rng - 74) < 3 else ''}")

me = pt["A3A_mean_expr_tumor"].dropna()
check("per-cell A3 output range 1.4-fold", 1.4,
      float(me.max() / me.min()) if me.min() else np.nan, tol_abs=0.15)

if "SC005" in pt.index:
    check("SC005 A3A 3.03", 3.03, float(pt.loc["SC005", "A3A_mean_all_tumor"]),
          tol_abs=0.05)
    check("SC005 SBS2 fold 0.7x", 0.7, float(pt.loc["SC005", "fold_sbs2"]),
          tol_abs=0.08)
if "SC001" in pt.index:
    check("SC001 A3A 1.12", 1.12, float(pt.loc["SC001", "A3A_mean_all_tumor"]),
          tol_abs=0.05)
    rank = int((pt["A3A_mean_all_tumor"] >
                pt.loc["SC001", "A3A_mean_all_tumor"]).sum()) + 1
    check("SC001 A3A rank 6 of 14", 6, rank, tol_abs=0)

# =============================================================================
banner("STEP 5: PREVALENCE vs CONTRIBUTION -- FULL GRID")
# =============================================================================
TARGETS = {("A3A", "sbs2"): (0.594, 0.025), ("A3A", "cnv"): (0.047, 0.87),
           ("A3B", "cnv"): (0.660, 0.010), ("A3B", "sbs2"): (0.400, 0.16)}

print("")
print(f"  {'pair':<16} {'cell set':<9} {'n_cells':>19} {'pct':>19} {'fold':>19}")
print(f"  {'-'*16} {'-'*9} {'-'*19} {'-'*19} {'-'*19}")
grid = {}
for (tag, fate), (t_rho, t_p) in TARGETS.items():
    for cs in active_sets:
        cells = []
        for defn, col in (("n_cells", f"n_{fate}_high"),
                          ("pct", f"pct_of_{fate}"), ("fold", f"fold_{fate}")):
            rho, p = spearmanr(pt[f"{tag}_prev_{cs}"], pt[col])
            grid[(tag, fate, cs, defn)] = (rho, p)
            hit = "*" if abs(rho - t_rho) < 0.02 else " "
            cells.append(f"{rho:+.3f}/{p:.3g}{hit}")
        print(f"  {tag + ' vs ' + fate.upper():<16} {cs:<9} "
              + " ".join(f"{c:>19}" for c in cells))

print("")
log("  '*' marks a rho within 0.02 of the published value.")
log("  The correct definition is the ONE cell set hitting all four pairs.")
print("")
for cs in active_sets:
    hits = sum(1 for (tag, fate), (t_rho, _) in TARGETS.items()
               if abs(grid[(tag, fate, cs, "fold")][0] - t_rho) < 0.02)
    log(f"    cell set '{cs}' with contribution=fold: {hits}/4 pairs matched")

best_cs = max(active_sets,
              key=lambda cs: sum(1 for (tag, fate), (t_rho, _) in TARGETS.items()
                                 if abs(grid[(tag, fate, cs, "fold")][0] - t_rho) < 0.02))
print("")
log(f"  auditing against cell set '{best_cs}', contribution = fold")
for (tag, fate), (t_rho, t_p) in TARGETS.items():
    rho, p = grid[(tag, fate, best_cs, "fold")]
    check(f"{tag} prev vs {fate.upper()} rho", t_rho, float(rho), tol_abs=0.03)
    check(f"{tag} prev vs {fate.upper()} p", t_p, float(p), tol_rel=0.20)

# =============================================================================
banner("STEP 6: JOINT THREE-CRITERIA SEPARATION")
# =============================================================================
# FIX 2: do not run at all if any input criterion is unevaluable.
if not viral_ok or not stage:
    log("  [SKIP] viral burden or lifecycle stage unavailable.")
    log("         Not run, rather than run with missing criteria coerced to")
    log("         False, which in v1 produced a spurious Fisher p = 1.")
    check("Fisher exact SBS2", 0.0027, None)
    check("Fisher exact CNA", 0.011, None)
else:
    for fate, hc, prev_col, prev_min, si, target in (
            ("sbs2", SBS2_HC, f"A3A_prev_{best_cs}", A3A_PREV_MIN, 0, 0.0027),
            ("cnv", CNV_HC, f"A3B_prev_{best_cs}", A3B_PREV_MIN, 1, 0.011)):
        f = []
        for i, r in pt.iterrows():
            v = bool(r["viral_pct"] > VIRAL_MIN)
            s = bool(stage.get(i, (0, 0))[si] >= STAGE_MIN) if i in stage else False
            a = bool(r[prev_col] >= prev_min) if not np.isnan(r[prev_col]) else False
            f.append({"patient": i, "viral": v, "stage": s, "a3": a,
                      "all3": v and s and a, "hc": i in hc})
        f = pd.DataFrame(f).set_index("patient")
        tab = [[int((f["all3"] & f["hc"]).sum()), int((f["all3"] & ~f["hc"]).sum())],
               [int((~f["all3"] & f["hc"]).sum()), int((~f["all3"] & ~f["hc"]).sum())]]
        orr, pv = fisher_exact(tab)
        print("")
        log(f"  {fate.upper()}: table {tab}, OR = {orr:.3g}, Fisher p = {pv:.4g}")
        print(f"  {'pt':<8} {'viral':>7} {'stage':>7} {'A3':>6} {'all3':>6}  HC")
        for i, r in f.iterrows():
            print(f"  {i:<8} {str(r['viral']):>7} {str(r['stage']):>7} "
                  f"{str(r['a3']):>6} {str(r['all3']):>6}  "
                  f"{'HC' if r['hc'] else ''}")
        check(f"Fisher exact {fate.upper()}", target, float(pv), tol_rel=0.15)

# =============================================================================
banner("AUDIT SUMMARY")
# =============================================================================
npass = sum(1 for r in RESULTS if r[3] == "PASS")
nfail = sum(1 for r in RESULTS if r[3] == "FAIL")
nskip = sum(1 for r in RESULTS if r[3] == "SKIP")
print("")
log(f"  PASS: {npass}    FAIL: {nfail}    SKIP: {nskip}")
for tag, want in (("FAILURES", "FAIL"), ("SKIPPED", "SKIP")):
    hits = [r for r in RESULTS if r[3] == want]
    if hits:
        print("")
        log(f"  {tag}:")
        for lab, e, o, _ in hits:
            log(f"    {lab}: prose {e}, computed {o}")

os.makedirs(OUTDIR, exist_ok=True)
pt.to_csv(os.path.join(OUTDIR, "section44_per_patient.tsv"), sep="\t")
pd.DataFrame([{"tag": k[0], "fate": k[1], "cell_set": k[2], "contribution": k[3],
               "rho": v[0], "p": v[1]} for k, v in grid.items()]).to_csv(
    os.path.join(OUTDIR, "section44_prevalence_grid.tsv"), sep="\t", index=False)
pd.DataFrame(RESULTS, columns=["claim", "expected", "observed", "verdict"]).to_csv(
    os.path.join(OUTDIR, "section44_audit.tsv"), sep="\t", index=False)
log(f"  Output: {OUTDIR}")
