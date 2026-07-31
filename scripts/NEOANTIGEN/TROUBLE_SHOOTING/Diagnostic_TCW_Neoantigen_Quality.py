#!/usr/bin/env python3
"""
Diagnostic_TCW_Neoantigen_Quality.py   (read-only)
==================================================
Within-mutation test: among the neoantigens these cells carry, do the APOBEC
(TCW-context) mutations preferentially become GOOD neoantigens?

Motivation
----------
The between-group TCW test is flat because SBS2-HIGH (A3A-active) and CNV-HIGH
(SBS2 = 0) have near-identical raw TpCpW-among-C>T fractions, both sitting at
the ~12.5% chance baseline (TpCpW is 2 of 16 NCN contexts). That comparison
cancels. The mechanistic question is not between groups but within the mutation
set, pooled across both groups for power: is an APOBEC-context substitution more
likely to create a de novo epitope, and does it produce a larger binding gain or
a stronger binder, than a non-APOBEC substitution?

What it does (reads data/FIG_7/06_prevalence_ranking/..._full.tsv)
------------------------------------------------------------------
Splits the 775 neoantigen mutations by TCW status and compares:
  1. binding gain   (delta_IC50 = wt_IC50 - mut_IC50)   Mann-Whitney
  2. binding strength (mut_IC50, lower = stronger)       Mann-Whitney
  run for is_tcw_ct (C>T only; the SBS2 arm; PRIMARY, since SBS13/C>G is not
  established in these tumors) and is_tcw (C>T + C>G; secondary).
Then the mechanistic "why": amino-acid change composition of TCW neoantigens
(APOBEC C>T at the codon produces E>K charge reversals), and whether E>K changes
carry the binding gain. Finally a Fisher test of whether TCW mutations are
over-represented among differential neoantigens (wt non-binder -> mut binder).

READ-ONLY. Env: NETWORK (pandas + scipy; no pysam / no genome needed).
Run: conda run -n NETWORK python Diagnostic_TCW_Neoantigen_Quality.py
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import re
from datetime import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
RANK_TSV = os.path.join(BASE_DIR, "data/FIG_7/06_prevalence_ranking/neoantigen_prevalence_ranking_full.tsv")
OUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/tcw_neoantigen_quality")

AA3TO1 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
          'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
          'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
          'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
ACIDIC = {'D', 'E'}
BASIC = {'K', 'R'}

_report = []
def log(m=""):
    print(m, flush=True)
    _report.append(str(m))
def banner(t, ch="="):
    log("")
    log(ch * 80)
    log(f"  {t}")
    log(ch * 80)

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
    """p.Glu342Lys -> ('E','K'); returns (wt1, mut1) or (None, None)."""
    m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)', str(s))
    if not m:
        return None, None
    wt = AA3TO1.get(m.group(1))
    mut = AA3TO1.get(m.group(3), '*' if m.group(3) in ('Ter', '*') else None)
    return wt, mut

# =============================================================================
# STATS
# =============================================================================
def mwu(x, y, lx, ly, note=""):
    """Two-sided Mann-Whitney with medians and common-language effect size."""
    from scipy.stats import mannwhitneyu
    x = np.asarray(pd.to_numeric(pd.Series(x), errors="coerce"), float); x = x[~np.isnan(x)]
    y = np.asarray(pd.to_numeric(pd.Series(y), errors="coerce"), float); y = y[~np.isnan(y)]
    if len(x) < 3 or len(y) < 3:
        log(f"    [SKIP] n too small ({lx}={len(x)}, {ly}={len(y)}) {note}")
        return
    U, p = mannwhitneyu(x, y, alternative="two-sided")
    cles = U / (len(x) * len(y))  # P(random x > random y)
    log(f"    {lx:14s} n={len(x):3d}  median={np.median(x):8.1f}  mean={np.mean(x):9.1f}")
    log(f"    {ly:14s} n={len(y):3d}  median={np.median(y):8.1f}  mean={np.mean(y):9.1f}")
    log(f"    Mann-Whitney U={U:.0f}, p={p:.4g}, P({lx} > {ly})={cles:.2f}   {note}")

def fisher2(a, b, c, d, label):
    from scipy.stats import fisher_exact
    orv, p = fisher_exact([[a, b], [c, d]])
    log(f"    {label}: [[{a},{b}],[{c},{d}]]  Fisher OR={orv:.2f}, p={p:.4g}")

# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("TCW NEOANTIGEN QUALITY (within-mutation, pooled across groups)")
    log(f"  {datetime.now().isoformat(timespec='seconds')}")
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(RANK_TSV):
        log(f"  [FATAL] ranking not found: {RANK_TSV}")
        return
    df = pd.read_csv(RANK_TSV, sep="\t")

    dcol = _first(df.columns, ["delta_IC50", "delta_ic50"])
    mcol = _first(df.columns, ["mut_IC50", "mut_ic50"])
    wcol = _first(df.columns, ["wt_IC50", "wt_ic50"])
    tcol = _first(df.columns, ["is_tcw"])
    ctcol = _first(df.columns, ["is_tcw_ct"])
    hcol = _first(df.columns, ["hgvs_p"])
    difcol = _first(df.columns, ["is_differential"])
    scol = _first(df.columns, ["sub_pyr"])

    df["_delta"] = pd.to_numeric(df[dcol], errors="coerce")
    df["_mut"] = pd.to_numeric(df[mcol], errors="coerce")
    df["_tcw"] = df[tcol].map(_tobool)
    df["_tcwct"] = df[ctcol].map(_tobool)
    df["_diff"] = df[difcol].map(_tobool) if difcol else False
    df[["_wt1", "_mut1"]] = df[hcol].apply(lambda s: pd.Series(parse_hgvsp(s)))
    df["_ek"] = (df["_wt1"] == "E") & (df["_mut1"] == "K")
    df["_chgrev"] = df.apply(
        lambda r: (r["_wt1"] in ACIDIC and r["_mut1"] in BASIC)
        or (r["_wt1"] in BASIC and r["_mut1"] in ACIDIC), axis=1)

    log(f"  neoantigen mutations: {len(df)}")
    log(f"  is_tcw (C>T + C>G): {int(df['_tcw'].sum())}   "
        f"is_tcw_ct (C>T only): {int(df['_tcwct'].sum())}   "
        f"differential: {int(df['_diff'].sum())}")

    # ---- 1-2. binding gain + strength, by TCW status --------------------------
    for split_col, split_name in [("_tcwct", "clean-TCW C>T (SBS2 arm) [PRIMARY]"),
                                   ("_tcw", "TCW C>T + C>G (SBS2 + SBS13)")]:
        banner(f"BINDING QUALITY by {split_name}", "-")
        tcw = df[df[split_col]]
        non = df[~df[split_col]]
        log("  Binding GAIN (delta_IC50 = wt_IC50 - mut_IC50; higher = bigger gain):")
        mwu(tcw["_delta"], non["_delta"], "TCW", "non-TCW",
            note="(does APOBEC create larger binding gains?)")
        log("\n  Binding STRENGTH (mut_IC50 nM; LOWER = stronger binder):")
        mwu(tcw["_mut"], non["_mut"], "TCW", "non-TCW",
            note="(P>0.5 here means TCW binds WEAKER)")

    # ---- 3. mechanistic why: amino-acid change composition --------------------
    banner("AMINO-ACID CHANGE COMPOSITION (why TCW mutations behave as they do)", "-")
    tcwct = df[df["_tcwct"]]
    log(f"  Among the {len(tcwct)} clean-TCW C>T neoantigens:")
    log(f"    E>K (Glu>Lys) charge reversals : {int(tcwct['_ek'].sum())} "
        f"({100 * tcwct['_ek'].mean():.0f}%)")
    log(f"    any acidic<->basic charge rev  : {int(tcwct['_chgrev'].sum())} "
        f"({100 * tcwct['_chgrev'].mean():.0f}%)")
    changes = (tcwct["_wt1"].fillna("?") + ">" + tcwct["_mut1"].fillna("?")).value_counts()
    log("    top substitutions:")
    for chg, n in changes.head(8).items():
        log(f"      {chg}: {n}")
    log(f"\n  For contrast, non-TCW E>K fraction: "
        f"{100 * df[~df['_tcwct']]['_ek'].mean():.1f}%  "
        f"(charge-rev {100 * df[~df['_tcwct']]['_chgrev'].mean():.1f}%)")

    log("\n  Is the binding gain really the E>K (charge-reversal) effect?")
    log("  Binding GAIN, E>K vs non-E>K (all neoantigens, pooled):")
    mwu(df[df["_ek"]]["_delta"], df[~df["_ek"]]["_delta"], "E>K", "non-E>K",
        note="(if E>K dominates the gain, TCW's effect is the charge reversal)")

    # ---- 4. do TCW mutations preferentially CREATE de novo epitopes? ----------
    banner("TCW vs DIFFERENTIAL (wt non-binder -> mut binder)", "-")
    if difcol:
        for split_col, name in [("_tcwct", "clean-TCW C>T"), ("_tcw", "TCW C>T + C>G")]:
            a = int((df[split_col] & df["_diff"]).sum())
            b = int((df[split_col] & ~df["_diff"]).sum())
            c = int((~df[split_col] & df["_diff"]).sum())
            d = int((~df[split_col] & ~df["_diff"]).sum())
            fisher2(a, b, c, d, f"{name} x differential")
        log("  (differential = the mutation creates the epitope, not just modifies a binder)")
    else:
        log("  [SKIP] is_differential column not found.")

    # ---- persist --------------------------------------------------------------
    keep = [c for c in [_first(df.columns, ["gene"]), hcol, tcol, ctcol, scol,
                        wcol, mcol, dcol, difcol] if c] + ["_wt1", "_mut1", "_ek", "_chgrev"]
    df[keep].to_csv(os.path.join(OUT_DIR, "tcw_neoantigen_quality_table.tsv"),
                    sep="\t", index=False)
    with open(os.path.join(OUT_DIR, "tcw_neoantigen_quality_report.txt"), "w") as f:
        f.write("\n".join(_report))
    log(f"\n  Wrote table + report to {OUT_DIR}")


if __name__ == "__main__":
    main()
