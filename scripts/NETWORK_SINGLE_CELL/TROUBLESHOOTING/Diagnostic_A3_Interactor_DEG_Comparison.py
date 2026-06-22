#!/usr/bin/env python3
"""
Diagnostic_A3_Interactor_DEG_Comparison.py
==========================================

DEG-level comparison of A3 interactors across the two tumor-vs-normal
networks, as a setup layer for the network section. The point: the
network does not produce de novo pathways, it shows which genes are
co-associated with a cell state. So DEG direction + network placement +
chain membership together argue functional role.

PARALOG ANNOTATION (transparency, not a specificity claim). Each interactor
is annotated with the A3 paralog bait(s) it was recovered with, read
directly from the A3_baits column of Harris_A3_interactors.txt (Jang 2024
AP-MS + McCann 2023). These baits are reported for transparency only. The
A3 paralogs are highly sequence-similar, and an overexpression affinity
pulldown is a different context from endogenous co-expression in tumor
cells, so we do NOT infer paralog selectivity in this system from the bait
assignment. (Concretely: RALY's documented interaction is with a non-A3A,
non-A3B paralog, yet it anchors the activator module adjacent to A3 here;
its relevance rests on co-expression and network placement, not on a
claimed A3A/A3B partnership.)

Builds a big-to-small funnel and shows findings at each level:
  Level 0  all Harris A3 interactors
  Level 1  interactors with DE data (per comparison)
  Level 2  significantly up in tumor (per comparison)
  Level 3  placed in the network partition (per comparison)
  Level 4  placed in an A3-containing community (per comparison)
  Level 5  placed in a concordant subnetwork (if the structure
           diagnostic's Harris index is on disk)

Direction buckets (significant direction in each comparison):
  up_both          up in tumor in BOTH  (activator overlap; RALY)
  down_both        up in normal in BOTH (candidate inhibitor lost in tumor)
  flip_up_to_down  up in SBS2, down in CNV (candidate brake)
  flip_down_to_up  down in SBS2, up in CNV

Ends on a RALY spotlight: its profile at every level.

Usage:
    conda run -n NETWORK python Diagnostic_A3_Interactor_DEG_Comparison.py
"""

import os
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4_ROOT = os.path.join(BASE_DIR, "data/FIG_4")
HARRIS_PATH = os.path.join(FIG4_ROOT, "00_input/Harris_A3_interactors.txt")
CONC_DIR = os.path.join(FIG4_ROOT, "DIAGNOSTIC_CONCORDANCE")
OUTPUT_DIR = os.path.join(FIG4_ROOT, "DIAGNOSTIC_DEG_COMPARISON")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COMPARISONS = [
    {"key": "SBS2", "name": "SBS2_VS_NORMAL", "label": "SBS2-HIGH vs NORMAL",
     "dir": os.path.join(FIG4_ROOT, "NETWORK_SBS2_VS_NORMAL")},
    {"key": "CNV", "name": "CNV_VS_NORMAL", "label": "CNV-HIGH vs NORMAL",
     "dir": os.path.join(FIG4_ROOT, "NETWORK_CNV_VS_NORMAL")},
]

A3_SYMBOLS = {"APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H"}
A3_ALIAS = {"APOBEC3A": "A3A", "APOBEC3B": "A3B"}
A3_TARGETS = ["APOBEC3A", "APOBEC3B"]

SIG_FDR = 0.05
FDR_CANDIDATES = ["FDR", "fdr", "pvals_adj", "padj", "pval_adj",
                  "p_adj", "qval", "qvalue", "adj_pval"]


def log(msg):
    print(msg, flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)

def section(title):
    print(f"\n  --- {title} ---", flush=True)


# =============================================================================
# LOAD
# =============================================================================

def resolve_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def parse_baits(s):
    if not isinstance(s, str):
        return set()
    return set(t.strip() for t in s.split(",") if t.strip())


def load_harris():
    """Return (all_genes, gene_to_baits) from the staged Harris file.

    A3_baits is column 2: the A3 paralog bait(s) each interactor was
    recovered with, comma-separated (read for transparency only).
    """
    h = pd.read_csv(HARRIS_PATH, sep="\t")
    genes = set(h["gene_symbol"].values)
    if "A3_baits" in h.columns:
        gene_to_baits = {r.gene_symbol: parse_baits(r.A3_baits)
                         for r in h.itertuples()}
    else:
        log("  [WARN] no A3_baits column; paralog annotation unavailable")
        gene_to_baits = {g: set() for g in genes}
    return genes, gene_to_baits


def load_comparison(cfg):
    out = {"key": cfg["key"], "label": cfg["label"]}

    de = pd.read_csv(os.path.join(cfg["dir"],
                     "02_differential_expression/SC_diffexpr_stats.csv"))
    fdr_col = resolve_col(de, FDR_CANDIDATES)
    out["has_fdr"] = fdr_col is not None
    out["fdr_col"] = fdr_col
    out["log2fc"] = dict(zip(de["gene"], de["log2FC"]))
    out["fdr"] = (dict(zip(de["gene"], de[fdr_col]))
                  if fdr_col is not None else {})

    part = pd.read_csv(os.path.join(cfg["dir"],
                       "04_communities/SC_best_partition.csv"))
    g2c = dict(zip(part["gene"], part["community"]))
    out["gene_to_comm"] = g2c

    comm_to_a3 = {}
    for a3g in A3_TARGETS:
        if a3g in g2c:
            comm_to_a3.setdefault(g2c[a3g], []).append(A3_ALIAS[a3g])
    out["a3_comms"] = set(comm_to_a3.keys())
    out["comm_to_a3"] = comm_to_a3

    idx_path = os.path.join(CONC_DIR,
                            f"{cfg['name']}_harris_subnetwork_index.tsv")
    gene_to_modes = {}
    if os.path.exists(idx_path):
        idx = pd.read_csv(idx_path, sep="\t")
        for _, r in idx.iterrows():
            gene_to_modes.setdefault(r["harris_gene"], set()).add(r["mode"])
    out["gene_to_modes"] = gene_to_modes
    out["has_index"] = bool(gene_to_modes)
    return out


def direction(gene, comp):
    fc = comp["log2fc"].get(gene)
    if fc is None or (isinstance(fc, float) and np.isnan(fc)):
        return "absent", None, None
    if comp["has_fdr"]:
        fdr = comp["fdr"].get(gene)
        if fdr is None or (isinstance(fdr, float) and np.isnan(fdr)) \
                or fdr >= SIG_FDR:
            return "NS", fc, fdr
        return ("up" if fc > 0 else "down"), fc, fdr
    return ("up" if fc > 0 else "down"), fc, None


def categorize(d_sbs2, d_cnv):
    if d_sbs2 == "up" and d_cnv == "up":
        return "up_both"
    if d_sbs2 == "down" and d_cnv == "down":
        return "down_both"
    if d_sbs2 == "up" and d_cnv == "down":
        return "flip_up_to_down"
    if d_sbs2 == "down" and d_cnv == "up":
        return "flip_down_to_up"
    return "other"


# =============================================================================
# BUILD MASTER TABLE
# =============================================================================

def build_master(harris_all, gene_to_baits, comps):
    sbs2, cnv = comps["SBS2"], comps["CNV"]
    genes = sorted(g for g in harris_all if g not in A3_SYMBOLS)

    rows = []
    for g in genes:
        d_s, fc_s, fdr_s = direction(g, sbs2)
        d_c, fc_c, fdr_c = direction(g, cnv)
        if d_s == "absent" and d_c == "absent":
            continue

        comm_s = sbs2["gene_to_comm"].get(g)
        comm_c = cnv["gene_to_comm"].get(g)
        a3comm_s = sbs2["comm_to_a3"].get(comm_s, []) if comm_s is not None else []
        a3comm_c = cnv["comm_to_a3"].get(comm_c, []) if comm_c is not None else []

        baits = gene_to_baits.get(g, set())
        rows.append({
            "gene": g,
            "A3_baits": ",".join(sorted(baits)),
            "binds_A3A": "A3A" in baits,
            "binds_A3B": "A3B" in baits,
            "log2FC_SBS2": round(fc_s, 4) if fc_s is not None else None,
            "FDR_SBS2": fdr_s,
            "dir_SBS2": d_s,
            "log2FC_CNV": round(fc_c, 4) if fc_c is not None else None,
            "FDR_CNV": fdr_c,
            "dir_CNV": d_c,
            "category": categorize(d_s, d_c),
            "comm_SBS2": comm_s,
            "in_A3comm_SBS2": ",".join(a3comm_s),
            "comm_CNV": comm_c,
            "in_A3comm_CNV": ",".join(a3comm_c),
            "subnet_SBS2": ",".join(sorted(sbs2["gene_to_modes"].get(g, []))),
            "subnet_CNV": ",".join(sorted(cnv["gene_to_modes"].get(g, []))),
        })
    return pd.DataFrame(rows)


# =============================================================================
# REPORTING
# =============================================================================

def funnel(df, comps, harris_all):
    sbs2, cnv = comps["SBS2"], comps["CNV"]
    banner("FUNNEL: A3 interactors, big scope to small")
    log(f"  {'level':45s}{'SBS2':>8s}{'CNV':>8s}")
    log(f"  {'-'*61}")

    n0 = len([g for g in harris_all if g not in A3_SYMBOLS])
    log(f"  {'0. Harris interactors (non-A3)':45s}{n0:>8d}{n0:>8d}")
    have = (df["dir_SBS2"] != "absent", df["dir_CNV"] != "absent")
    log(f"  {'1. with DE data':45s}{have[0].sum():>8d}{have[1].sum():>8d}")
    up = (df["dir_SBS2"] == "up", df["dir_CNV"] == "up")
    log(f"  {'2. significantly up in tumor':45s}{up[0].sum():>8d}{up[1].sum():>8d}")
    inpart = (df["comm_SBS2"].notna(), df["comm_CNV"].notna())
    log(f"  {'3. placed in network partition':45s}"
        f"{inpart[0].sum():>8d}{inpart[1].sum():>8d}")
    ina3 = (df["in_A3comm_SBS2"] != "", df["in_A3comm_CNV"] != "")
    log(f"  {'4. placed in an A3 community':45s}"
        f"{ina3[0].sum():>8d}{ina3[1].sum():>8d}")
    if sbs2["has_index"] or cnv["has_index"]:
        insub = (df["subnet_SBS2"] != "", df["subnet_CNV"] != "")
        log(f"  {'5. in a concordant subnetwork':45s}"
            f"{insub[0].sum():>8d}{insub[1].sum():>8d}")
        act = (df["subnet_SBS2"].str.contains("activator"),
               df["subnet_CNV"].str.contains("activator"))
        log(f"  {'5a. in an ACTIVATOR subnetwork':45s}"
            f"{act[0].sum():>8d}{act[1].sum():>8d}")
    else:
        log("  (level 5 skipped: run Diagnostic_A3_Network_Structure.py first)")


def cross_tab(df):
    section("Direction cross-tab (rows SBS2, cols CNV)")
    order = ["up", "down", "NS", "absent"]
    log("    " + "".join(f"{c:>8s}" for c in order))
    for r in order:
        cells = [f"{((df['dir_SBS2'] == r) & (df['dir_CNV'] == c)).sum():>8d}"
                 for c in order]
        log(f"  {r:>5s}" + "".join(cells))


def paralog_summary(sub):
    nA = int(sub["binds_A3A"].sum())
    nB = int(sub["binds_A3B"].sum())
    nOther = int(((~sub["binds_A3A"]) & (~sub["binds_A3B"])).sum())
    log(f"    -> baits (transparency, not specificity): binds A3A {nA}, "
        f"binds A3B {nB}, other paralog only {nOther}")


def report_categories(df, title):
    banner(title)
    cross_tab(df)
    labels = {
        "up_both": "UP in both (activator overlap)",
        "down_both": "DOWN in both / up in normal (candidate lost inhibitor)",
        "flip_up_to_down": "FLIP up(SBS2) -> down(CNV) (candidate brake)",
        "flip_down_to_up": "FLIP down(SBS2) -> up(CNV)",
    }
    for cat, lbl in labels.items():
        sub = df[df["category"] == cat].copy()
        section(f"{lbl}: {len(sub)}")
        if len(sub) == 0:
            log("    (none)")
            continue
        sub = sub.sort_values(["binds_A3A", "binds_A3B"], ascending=False)
        for _, r in sub.iterrows():
            loc = []
            if r["in_A3comm_SBS2"]:
                loc.append(f"SBS2:C{r['comm_SBS2']}({r['in_A3comm_SBS2']})")
            if r["in_A3comm_CNV"]:
                loc.append(f"CNV:C{r['comm_CNV']}({r['in_A3comm_CNV']})")
            loc_s = ("  [" + "; ".join(loc) + "]") if loc else ""
            baits = r["A3_baits"] or "none"
            log(f"    {r['gene']:<14s} "
                f"SBS2 log2FC={r['log2FC_SBS2']}, CNV log2FC={r['log2FC_CNV']}"
                f"  baits=[{baits}]{loc_s}")
        paralog_summary(sub)


def raly_spotlight(df, comps):
    banner("RALY SPOTLIGHT (the through-line)")
    r = df[df["gene"] == "RALY"]
    if len(r) == 0:
        log("  RALY not found in interactor DE tables."); return
    r = r.iloc[0]
    log(f"  A3 baits:   {r['A3_baits'] or 'none'}  "
        f"(binds A3A={r['binds_A3A']}, binds A3B={r['binds_A3B']})")
    log(f"  DEG:        SBS2 log2FC={r['log2FC_SBS2']} ({r['dir_SBS2']}), "
        f"CNV log2FC={r['log2FC_CNV']} ({r['dir_CNV']})  -> {r['category']}")
    log(f"  A3 comm:    SBS2 {r['in_A3comm_SBS2'] or 'no'}, "
        f"CNV {r['in_A3comm_CNV'] or 'no'}")
    log(f"  subnetwork: SBS2 {r['subnet_SBS2'] or 'n/a'}, "
        f"CNV {r['subnet_CNV'] or 'n/a'}")
    for c in comps.values():
        fc = c["log2fc"].get("APOBEC3A")
        log(f"  A3A in {c['label']}: log2FC={fc:.3f}"
            if fc is not None else f"  A3A in {c['label']}: absent")
    log("  Read: RALY is a documented A3-family interactor (non-A3A/A3B bait) "
        "that stays up and activator-anchored in both phases. Its relevance "
        "rests on co-expression and network placement, not a paralog claim; "
        "A3A itself falls in CNV so that route cannot drive mutations there.")


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("A3 INTERACTOR DEG-DIRECTION COMPARISON")

    harris_all, gene_to_baits = load_harris()
    n_a3a = sum(1 for b in gene_to_baits.values() if "A3A" in b)
    n_a3b = sum(1 for b in gene_to_baits.values() if "A3B" in b)
    log(f"Harris interactors: {len(harris_all)} total "
        f"(bait annotation: {n_a3a} bind A3A, {n_a3b} bind A3B)")

    comps = {}
    for cfg in COMPARISONS:
        log(f"Loading {cfg['name']}...")
        comps[cfg["key"]] = load_comparison(cfg)
        c = comps[cfg["key"]]
        log(f"  FDR column: {c['fdr_col'] or 'NONE (direction by sign only)'}; "
            f"A3 communities: {dict(c['comm_to_a3'])}; "
            f"subnet index: {'yes' if c['has_index'] else 'no'}")

    df = build_master(harris_all, gene_to_baits, comps)

    funnel(df, comps, harris_all)
    report_categories(df, "GLOBAL SCOPE: all A3 interactors")

    in_a3 = df[(df["in_A3comm_SBS2"] != "") | (df["in_A3comm_CNV"] != "")]
    report_categories(in_a3,
                      "A3-COMMUNITY SCOPE: interactors the network put near A3")

    if comps["SBS2"]["has_index"] or comps["CNV"]["has_index"]:
        banner("CHAIN SCOPE: interactors in concordant subnetworks")
        act = df[(df["subnet_SBS2"].str.contains("activator")) |
                 (df["subnet_CNV"].str.contains("activator"))]
        section(f"In an ACTIVATOR subnetwork: {len(act)}")
        for _, r in act.iterrows():
            log(f"    {r['gene']:<14s} baits=[{r['A3_baits'] or 'none'}], "
                f"cat={r['category']}, SBS2 subnet={r['subnet_SBS2'] or '-'}, "
                f"CNV subnet={r['subnet_CNV'] or '-'}")
        rep = df[(df["subnet_SBS2"].str.contains("repressor")) |
                 (df["subnet_CNV"].str.contains("repressor"))]
        section(f"In a REPRESSOR subnetwork: {len(rep)}")
        for _, r in rep.iterrows():
            log(f"    {r['gene']:<14s} baits=[{r['A3_baits'] or 'none'}], "
                f"cat={r['category']}")

    raly_spotlight(df, comps)

    df.sort_values(["category", "binds_A3A", "binds_A3B", "gene"],
                   ascending=[True, False, False, True], inplace=True)
    master_path = os.path.join(OUTPUT_DIR, "A3_interactor_DEG_comparison.tsv")
    df.to_csv(master_path, sep="\t", index=False)

    summ = (df.groupby("category")
              .agg(n=("gene", "size"),
                   n_bind_A3A=("binds_A3A", "sum"),
                   n_bind_A3B=("binds_A3B", "sum"))
              .reset_index())
    summ_path = os.path.join(OUTPUT_DIR, "A3_interactor_DEG_summary.tsv")
    summ.to_csv(summ_path, sep="\t", index=False)

    banner("DIAGNOSTIC COMPLETE")
    log(f"  [SAVE] {master_path}")
    log(f"  [SAVE] {summ_path}")


if __name__ == "__main__":
    main()
