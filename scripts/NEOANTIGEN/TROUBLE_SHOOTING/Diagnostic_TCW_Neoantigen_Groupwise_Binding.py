#!/usr/bin/env python3
"""
Diagnostic_TCW_Neoantigen_Groupwise_Binding.py   (read-only)
============================================================
Does A3A-derived TCW beat SBS5-derived TCW at making MHC-I binders?

The hypothesis (Jake): CNV-HIGH is filtered to SBS2 = 0, so its TCW C>T
mutations come from broad background processes (SBS5 and friends), while
SBS2-HIGH's TCW C>T mutations are enriched for genuine A3A. If A3A deamination
is a "targeted" mutagen, then among TCW-context neoantigens the SBS2-tier ones
should bind MHC-I more strongly than the CNV-tier ones. This script tests that.

SCOPE / CAVEAT (printed in the output too)
------------------------------------------
This is CONDITIONAL on being a neoantigen: full.tsv holds binders only, so this
asks "given a TCW mutation makes a binder, is the SBS2-tier binder stronger than
the CNV-tier binder?" It does NOT address whether A3A-TCW become binders more
often (that is the unconditional follow-on, from all_peptide_results). Group
sizes are small, so this is a rigor-checked direction-finder, not a definitive
test.

Rigor for small N
-----------------
Primary metric: mut_IC50 (as featured); sensitivity: best_mut_IC50_any_allele.
Primary split: tier (SBS2_specific vs CNV_specific); sensitivity: carrier-
dominant group. For every comparison it reports, on the clean-TCW C>T set:
  - exact Mann-Whitney (asymptotic only if ties force it),
  - a label-permutation test on the median difference (PRIMARY inference),
  - Hodges-Lehmann shift + bootstrap 95% CI,
  - common-language effect size P(SBS2 stronger),
and it dumps every data point and checks the E>K confound.

READ-ONLY. Env: NETWORK (pandas + scipy). No pysam / no genome needed.
Run:            conda run -n NETWORK python Diagnostic_TCW_Neoantigen_Groupwise_Binding.py
Validate stats: conda run -n NETWORK python Diagnostic_TCW_Neoantigen_Groupwise_Binding.py --selftest
Custom input:   ... --input /path/to/full.tsv
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import re
import sys
from datetime import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
RANK_TSV = os.path.join(BASE_DIR, "data/FIG_7/06_prevalence_ranking/neoantigen_prevalence_ranking_full.tsv")
OUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/tcw_groupwise_binding")

SEED = 42
N_PERM = 20000
N_BOOT = 20000

T_SBS2, T_SHARED, T_CNV = "SBS2_specific", "shared", "CNV_specific"
AA3TO1 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
          'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
          'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
          'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}

_report = []
def log(m=""):
    print(m, flush=True)
    _report.append(str(m))
def banner(t, ch="="):
    log("")
    log(ch * 78)
    log(f"  {t}")
    log(ch * 78)

# =============================================================================
# STATS  (pure, unit-tested by --selftest)
# =============================================================================
def prob_less(x, y):
    """P(random x < random y), ties count 0.5. Lower IC50 = stronger binder."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    d = x[:, None] - y[None, :]
    return float((np.sum(d < 0) + 0.5 * np.sum(d == 0)) / (len(x) * len(y)))

def hodges_lehmann(x, y):
    """Median of all pairwise (x_i - y_j); robust location-shift estimate."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.median((x[:, None] - y[None, :]).ravel()))

def mannwhitney(x, y, alternative="two-sided"):
    from scipy.stats import mannwhitneyu
    x = np.asarray(x, float); y = np.asarray(y, float)
    pooled = np.concatenate([x, y])
    ties = len(np.unique(pooled)) < len(pooled)
    method = "asymptotic" if ties else "exact"
    U, p = mannwhitneyu(x, y, alternative=alternative, method=method)
    return float(U), float(p), method

def perm_test_median_diff(x, y, n_perm=N_PERM, seed=SEED):
    """Two-sided permutation p on the median difference (assumption-free)."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float); y = np.asarray(y, float)
    obs = float(np.median(x) - np.median(y))
    pool = np.concatenate([x, y]); n = len(x)
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(pool)
        d = np.median(perm[:n]) - np.median(perm[n:])
        if abs(d) >= abs(obs) - 1e-9:
            ge += 1
    return obs, (ge + 1) / (n_perm + 1)   # add-one smoothing

def bootstrap_ci(x, y, n_boot=N_BOOT, seed=SEED, ci=95):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float); y = np.asarray(y, float)
    stats = np.empty(n_boot)
    for i in range(n_boot):
        xb = rng.choice(x, len(x), replace=True)
        yb = rng.choice(y, len(y), replace=True)
        stats[i] = hodges_lehmann(xb, yb)
    lo = float(np.percentile(stats, (100 - ci) / 2))
    hi = float(np.percentile(stats, 100 - (100 - ci) / 2))
    return lo, hi

# =============================================================================
# SELF-TEST  (validate each statistic against a hand-checkable example)
# =============================================================================
def selftest():
    banner("SELF-TEST: validating statistics on controlled inputs")
    ok = True

    # clean separation: x all below y
    x, y = [1.0, 2.0, 3.0], [10.0, 20.0, 30.0]
    pl = prob_less(x, y)
    hl = hodges_lehmann(x, y)
    U, p, meth = mannwhitney(x, y, "two-sided")
    _, pp = perm_test_median_diff(x, y, n_perm=5000, seed=1)
    lo, hi = bootstrap_ci(x, y, n_boot=2000, seed=1)
    log(f"  [sep] prob_less={pl} (exp 1.0); HL={hl} (exp -18.0); "
        f"MW U={U} p={p:.3g} method={meth}; perm p={pp:.3g}; boot95=[{lo:.1f},{hi:.1f}]")
    ok &= abs(pl - 1.0) < 1e-9
    ok &= abs(hl - (-18.0)) < 1e-9
    ok &= abs(U - 0.0) < 1e-9            # no x>y pairs
    ok &= (meth == "exact")             # no ties
    ok &= abs(p - 0.1) < 1e-6           # exact two-sided p for n1=n2=3, U=0
    ok &= (pp < 0.2)
    ok &= (hi < 0)                      # CI entirely below 0

    # identical distributions: no effect
    x, y = [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]
    pl = prob_less(x, y); hl = hodges_lehmann(x, y)
    U, p, meth = mannwhitney(x, y, "two-sided")
    _, pp = perm_test_median_diff(x, y, n_perm=5000, seed=1)
    log(f"  [null] prob_less={pl} (exp 0.5); HL={hl} (exp 0.0); "
        f"MW U={U} p={p:.3g} method={meth}; perm p={pp:.3g}")
    ok &= abs(pl - 0.5) < 1e-9
    ok &= abs(hl - 0.0) < 1e-9
    ok &= abs(U - 8.0) < 1e-9           # n1*n2/2
    ok &= (meth == "asymptotic")        # ties present -> asymptotic
    ok &= (pp > 0.5)

    # directional one-sided: x < y should be significant one-sided at this sep
    x, y = [1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0, 10.0]
    _, p_less, _ = mannwhitney(x, y, "less")
    log(f"  [dir] one-sided p(x<y)={p_less:.4g} (should be small)")
    ok &= (p_less < 0.05)

    banner("SELF-TEST RESULT: " + ("PASS" if ok else "**FAIL**"), "-")
    return ok

# =============================================================================
# LOADING
# =============================================================================
def _first(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in cols:
            return c
        if c.lower() in low:
            return low[c.lower()]
    return None

def _tobool(s):
    return str(s).strip().lower() in ("true", "1", "1.0", "yes")

def parse_hgvsp(s):
    m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)', str(s))
    if not m:
        return None, None
    wt = AA3TO1.get(m.group(1))
    mut = AA3TO1.get(m.group(3), '*' if m.group(3) in ('Ter', '*') else None)
    return wt, mut

def load(input_path):
    df = pd.read_csv(input_path, sep="\t")
    C = {
        "gene": _first(df.columns, ["gene"]),
        "hgvs": _first(df.columns, ["hgvs_p"]),
        "tier": _first(df.columns, ["tier"]),
        "mut": _first(df.columns, ["mut_IC50", "mut_ic50"]),
        "wt": _first(df.columns, ["wt_IC50", "wt_ic50"]),
        "delta": _first(df.columns, ["delta_IC50", "delta_ic50"]),
        "best": _first(df.columns, ["best_mut_IC50_any_allele", "best_mut_ic50_any_allele"]),
        "tcw": _first(df.columns, ["is_tcw"]),
        "tcwct": _first(df.columns, ["is_tcw_ct"]),
        "cs": _first(df.columns, ["carriers_sbs2"]),
        "cc": _first(df.columns, ["carriers_cnv"]),
        "sub": _first(df.columns, ["sub_pyr"]),
    }
    df["_mut"] = pd.to_numeric(df[C["mut"]], errors="coerce")
    df["_best"] = pd.to_numeric(df[C["best"]], errors="coerce") if C["best"] else df["_mut"]
    df["_tcw"] = df[C["tcw"]].map(_tobool)
    df["_tcwct"] = df[C["tcwct"]].map(_tobool)
    df["_cs"] = pd.to_numeric(df[C["cs"]], errors="coerce") if C["cs"] else np.nan
    df["_cc"] = pd.to_numeric(df[C["cc"]], errors="coerce") if C["cc"] else np.nan
    df[["_wt1", "_mut1"]] = df[C["hgvs"]].apply(lambda s: pd.Series(parse_hgvsp(s)))
    df["_ek"] = (df["_wt1"] == "E") & (df["_mut1"] == "K")
    df["_carrier_grp"] = np.where(df["_cs"] > df["_cc"], "SBS2",
                                  np.where(df["_cc"] > df["_cs"], "CNV", "tie"))
    return df, C

# =============================================================================
# REPORTING
# =============================================================================
def compare(name, xdf, ydf, col, lx, ly):
    x = xdf[col].dropna().values
    y = ydf[col].dropna().values
    log(f"\n  [{name}]  {lx} (n={len(x)})  vs  {ly} (n={len(y)})   metric={col}")
    if len(x) < 2 or len(y) < 2:
        log(f"    [SKIP] a group has <2 values; report the points above and treat as descriptive.")
        return
    log(f"    {lx}: median={np.median(x):8.1f}  mean={np.mean(x):9.1f}  "
        f"range=[{np.min(x):.1f}, {np.max(x):.1f}]")
    log(f"    {ly}: median={np.median(y):8.1f}  mean={np.mean(y):9.1f}  "
        f"range=[{np.min(y):.1f}, {np.max(y):.1f}]")
    U2, p2, m = mannwhitney(x, y, "two-sided")
    log(f"    Mann-Whitney two-sided ({m}): U={U2:.0f}, p={p2:.4g}")
    Ul, pl, _ = mannwhitney(x, y, "less")   # H1: x (SBS2) < y (CNV) == stronger
    log(f"    directional H1 [{lx} < {ly}, i.e. stronger]: one-sided p={pl:.4g}, "
        f"P({lx} stronger)={prob_less(x, y):.2f}")
    hl = hodges_lehmann(x, y)
    lo, hi = bootstrap_ci(x, y)
    log(f"    Hodges-Lehmann shift ({lx}-{ly}) = {hl:.1f} nM   bootstrap95% = [{lo:.1f}, {hi:.1f}]")
    obs, pp = perm_test_median_diff(x, y)
    log(f"    PERMUTATION (median diff = {obs:.1f}, two-sided) p = {pp:.4g}   [primary inference]")

def dump_rows(d, C):
    log(f"\n  All {len(d)} clean-TCW C>T neoantigens (the entire evidence base):")
    log(f"    {'tier':14s} {'gene':10s} {'hgvs_p':14s} {'sub':5s} {'E>K':3s} "
        f"{'cS':>4s} {'cC':>4s} {'wt_IC50':>9s} {'mut_IC50':>9s} {'best':>8s}")
    dd = d.sort_values([C["tier"], "_mut"])
    for _, r in dd.iterrows():
        cs = "" if pd.isna(r["_cs"]) else str(int(r["_cs"]))
        cc = "" if pd.isna(r["_cc"]) else str(int(r["_cc"]))
        wt = float(r[C["wt"]]) if (C["wt"] and not pd.isna(r[C["wt"]])) else float("nan")
        log(f"    {str(r[C['tier']]):14s} {str(r[C['gene']]):10s} "
            f"{str(r[C['hgvs']]):14s} {str(r[C['sub']]):5s} "
            f"{'Y' if r['_ek'] else '-':3s} "
            f"{cs:>4} {cc:>4} {wt:9.1f} {r['_mut']:9.1f} {r['_best']:8.1f}")

# =============================================================================
# MAIN
# =============================================================================
def main(input_path):
    banner("TCW NEOANTIGEN BINDING BY GROUP OF ORIGIN (SBS2-tier vs CNV-tier)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    log(f"  seed={SEED}, permutations={N_PERM}, bootstraps={N_BOOT}")
    log(f"  input: {input_path}")
    log("  SCOPE: conditional on being a neoantigen (full.tsv = binders only);")
    log("  tests strength GIVEN a binder, not the odds of becoming one.")
    os.makedirs(OUT_DIR, exist_ok=True)
    if not os.path.exists(input_path):
        log(f"  [FATAL] input not found: {input_path}")
        return
    df, C = load(input_path)

    # PRIMARY set: clean-TCW C>T (SBS2 arm; excludes SBS13/C>G that is not established)
    ct = df[df["_tcwct"]].copy()
    sbs2 = ct[ct[C["tier"]] == T_SBS2]
    cnv = ct[ct[C["tier"]] == T_CNV]
    shared = ct[ct[C["tier"]] == T_SHARED]

    banner("GROUP SIZES (clean-TCW C>T)", "-")
    log(f"  SBS2_specific: {len(sbs2)}   CNV_specific: {len(cnv)}   shared: {len(shared)}   "
        f"total clean-TCW C>T: {len(ct)}")
    if min(len(sbs2), len(cnv)) < 5:
        log("  [WARNING] at least one group < 5. Treat p-values as a direction-finder;")
        log("  read the raw points and the permutation p, not the asymptotic p.")

    dump_rows(ct, C)

    banner("PRIMARY: binding strength, SBS2-tier vs CNV-tier (clean-TCW C>T)", "-")
    log("  Lower IC50 = stronger binder. Hypothesis: SBS2-tier (A3A) binds stronger.")
    compare("PRIMARY mut_IC50", sbs2, cnv, "_mut", "SBS2-tier", "CNV-tier")
    compare("SENSITIVITY best_mut_IC50_any_allele", sbs2, cnv, "_best", "SBS2-tier", "CNV-tier")

    # ---- confound 1: E>K composition differs by tier? ------------------------
    banner("CONFOUND CHECK: is any strength difference just E>K composition?", "-")
    log(f"  E>K fraction: SBS2-tier {int(sbs2['_ek'].sum())}/{len(sbs2)} "
        f"({100 * sbs2['_ek'].mean():.0f}%), CNV-tier {int(cnv['_ek'].sum())}/{len(cnv)} "
        f"({100 * cnv['_ek'].mean():.0f}%)")
    if len(sbs2) and len(cnv):
        from scipy.stats import fisher_exact
        a = int(sbs2['_ek'].sum()); b = len(sbs2) - a
        c = int(cnv['_ek'].sum()); d = len(cnv) - c
        orv, pf = fisher_exact([[a, b], [c, d]])
        log(f"  Fisher E>K x tier: [[{a},{b}],[{c},{d}]] OR={orv:.2f}, p={pf:.4g}")
    log("\n  Sensitivity: repeat the strength test with E>K removed (does it survive?):")
    compare("mut_IC50, E>K excluded", sbs2[~sbs2["_ek"]], cnv[~cnv["_ek"]], "_mut",
            "SBS2-tier", "CNV-tier")

    # ---- confound 2: carrier-dominant group instead of tier ------------------
    banner("ROBUSTNESS: split by carrier-dominant group instead of neoantigen tier", "-")
    cs = ct[ct["_carrier_grp"] == "SBS2"]
    cc = ct[ct["_carrier_grp"] == "CNV"]
    log(f"  carrier-dominant SBS2: {len(cs)}   CNV: {len(cc)}   tie: {int((ct['_carrier_grp']=='tie').sum())}")
    compare("mut_IC50 by carrier group", cs, cc, "_mut", "SBS2-carried", "CNV-carried")

    # ---- sensitivity: broaden to is_tcw (C>T + C>G) --------------------------
    banner("SENSITIVITY: broaden to TCW C>T + C>G (adds SBS13 arm, not established)", "-")
    tcw = df[df["_tcw"]].copy()
    s2 = tcw[tcw[C["tier"]] == T_SBS2]
    cv = tcw[tcw[C["tier"]] == T_CNV]
    log(f"  SBS2_specific: {len(s2)}   CNV_specific: {len(cv)}")
    compare("mut_IC50 (C>T + C>G)", s2, cv, "_mut", "SBS2-tier", "CNV-tier")

    # ---- persist -------------------------------------------------------------
    keep = [c for c in [C["gene"], C["hgvs"], C["tier"], C["sub"], C["wt"], C["mut"],
                        C["best"], C["cs"], C["cc"]] if c] + ["_ek", "_carrier_grp"]
    ct[keep].to_csv(os.path.join(OUT_DIR, "tcw_groupwise_binding_table.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT_DIR, "tcw_groupwise_binding_report.txt"), "w") as f:
        f.write("\n".join(_report))
    log(f"\n  Wrote table + report to {OUT_DIR}")
    banner("READ-OUT", "-")
    log("  If PRIMARY permutation p is small AND SBS2-tier median is lower AND the")
    log("  effect survives E>K exclusion, that supports: A3A-derived TCW mutations")
    log("  produce stronger MHC-I binders than SBS5/background TCW. If it only holds")
    log("  with E>K in, the effect is the charge-reversal composition, not the signature")
    log("  per se. If flat, the between-group functional difference is not established")
    log("  at this N and stays a future direction.")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    inp = RANK_TSV
    if "--input" in sys.argv:
        inp = sys.argv[sys.argv.index("--input") + 1]
    main(inp)
