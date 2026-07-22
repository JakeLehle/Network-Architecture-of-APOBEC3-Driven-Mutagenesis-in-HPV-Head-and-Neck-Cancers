#!/usr/bin/env python3
"""
Step05_Integrated_Neoantigen_Analysis.py  (rewritten: prevalence-first)
=======================================================================
Figure 7 integrated analysis, rebuilt around the prevalence-weighted ranking.

What this version does (and does not do)
----------------------------------------
- Reads neoantigen tier and per-niche prevalence from the single source,
  neoantigen_prevalence_ranking_full.tsv (per neoantigen mutation, not per gene).
- Confirms the neoantigen-level overlap (SBS2-specific / shared / CNV-specific)
  independently from the {group}_neoantigens.tsv files and asserts it matches
  the ranking's tiers (expected 467 / 93 / 215).
- Runs the selection-evidence analysis on the SHARED tier only: for each shared
  neoantigen, flags whether at least one active loss mechanism is present in
  CNV-HIGH (antigen loss, fusion disruption, or expression silencing), counts how
  many of the shared neoantigens carry evidence, and lists them.
- Folds in the per-gene aggregates that the retired Step04 used to provide
  (neoantigen and strong-binder counts per group), computed here from
  {group}_all_peptide_results.tsv.

What was removed (superseded by the prevalence ranking)
-------------------------------------------------------
- The composite vaccine_score and its breadth / pressure multipliers.
- The hardcoded key-target list (ANXA1, MDK, ...).
- The group-level neo:loss population statistic (antigen loss is now only a
  per-neoantigen mechanism flag).
- Dead imports and the unused shared_neoantigen_comparison load.

Mechanism definitions (gene-level signals; only the two strong, direct ones)
---------------------------------------------------------------------------
- fusion_disrupted: the gene appears in the cross-group fusion overlap as a
                   CNV fusion that removes an SBS2 neoantigen junction.
- silenced       : the gene's mean expression in CNV is < 0.5x its SBS2 mean
                   (>2-fold down) while still expressed in SBS2.
A shared neoantigen carries selection evidence if it shows at least one of these
two. The earlier antigen-loss flag (a gene-level, wt-binds/mut-destroys signal
on a DIFFERENT mutation than the neoantigen) was dropped as too weak and
indirect; it inflated the count and pulled in CNV-enriched, direction-wrong hits.

HLA-A/B/C neoantigens are retained in the binding analysis and the target table
(their MHC-I binding gains are real) but are EXCLUDED from the selection-evidence
count: those loci are hyperpolymorphic and mapping-artifact-prone, their
apparent gains sit at germline-like prevalence, and HLA loss is its own escape
story rather than a neoantigen being selected against. They are flagged, not
silently dropped, and the exclusion is documented in Methods.

Per-niche prevalence (prevalence_sbs2 vs prevalence_cnv) is reported alongside as
the descriptive "minority of CNV cells still carry it" observation, with a
depleted / flat / enriched breakdown, but it is NOT the evidence trigger.

Run in NETWORK conda env (needs anndata for the expression/silencing step).
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MHC_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/03_mhc_binding")
FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis")
RANK_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/06_prevalence_ranking")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/05_summary")

FULL_TSV = os.path.join(RANK_DIR, "neoantigen_prevalence_ranking_full.tsv")
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")

MHC_BIND = 500.0        # nM
SILENCE_RATIO = 0.5     # CNV mean < 0.5x SBS2 mean = >2-fold down
SILENCE_FLOOR = 1.0     # require SBS2 mean expression above this to call silencing
EXPECTED_TIERS = {"SBS2_specific": 467, "shared": 93, "CNV_specific": 215}
HLA_PREFIXES = ("HLA-A", "HLA-B", "HLA-C")   # excluded from selection-evidence count

# =============================================================================
# LOGGING
# =============================================================================
_report = []
def log(msg=""):
    print(msg, flush=True)
    _report.append(str(msg))
def banner(title, ch="="):
    log(""); log(ch * 80); log(f"  {title}"); log(ch * 80)
def _first(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in cols: return c
        if c.lower() in low: return low[c.lower()]
    return None


# =============================================================================
# BEAT 0: SCHEMA VERIFICATION
# =============================================================================
def beat0():
    banner("BEAT 0: SCHEMA + INPUT VERIFICATION")
    st = {}
    for label, path, required in [
        ("prevalence ranking", FULL_TSV, True),
        ("group assignments", GROUP_PATH, True),
        ("SBS2 peptides", os.path.join(MHC_DIR, "SBS2_HIGH_all_peptide_results.tsv"), True),
        ("CNV peptides", os.path.join(MHC_DIR, "CNV_HIGH_all_peptide_results.tsv"), True),
        ("SBS2 neoantigens", os.path.join(MHC_DIR, "SBS2_HIGH_neoantigens.tsv"), True),
        ("CNV neoantigens", os.path.join(MHC_DIR, "CNV_HIGH_neoantigens.tsv"), True),
    ]:
        if not os.path.exists(path):
            if required:
                log(f"  [FATAL] missing {label}: {path}"); sys.exit(1)
    log("  Required tables present.")

    full = pd.read_csv(FULL_TSV, sep="\t")
    need = {"gene", "hgvs_p", "tier", "peptide", "wt_IC50", "mut_IC50",
            "prevalence_sbs2", "prevalence_cnv", "prevalence_tier"}
    miss = need - set(full.columns)
    if miss:
        log(f"  [FATAL] full.tsv missing columns: {miss}"); sys.exit(1)
    log(f"  full.tsv: {full.shape}; tiers present: {sorted(full['tier'].unique())}")
    st["full"] = full

    groups = pd.read_csv(GROUP_PATH, sep="\t")
    bc = _first(groups.columns, ["cell_barcode", "CB", "barcode"])
    gp = _first(groups.columns, ["group", "assignment"])
    st["sbs2_cells"] = set(groups.loc[groups[gp] == "SBS2_HIGH", bc].astype(str))
    st["cnv_cells"] = set(groups.loc[groups[gp] == "CNV_HIGH", bc].astype(str))
    log(f"  cells: SBS2 {len(st['sbs2_cells'])}, CNV {len(st['cnv_cells'])}")

    # peptide tables (for antigen loss + per-gene aggregates)
    pep = {}
    for grp in ("SBS2_HIGH", "CNV_HIGH"):
        p = os.path.join(MHC_DIR, f"{grp}_all_peptide_results.tsv")
        df = pd.read_csv(p, sep="\t")
        g = _first(df.columns, ["gene"]); m = _first(df.columns, ["mut_ic50", "mut_IC50"])
        w = _first(df.columns, ["wt_ic50", "wt_IC50"])
        sb = _first(df.columns, ["is_strong_binder"])
        if None in (g, m, w):
            log(f"  [FATAL] {grp} peptides missing gene/ic50 columns"); sys.exit(1)
        df = df.rename(columns={g: "gene", m: "mut_ic50", w: "wt_ic50"})
        if sb: df = df.rename(columns={sb: "is_strong_binder"})
        pep[grp] = df
        log(f"  {grp} peptides: {df.shape}")
    st["pep"] = pep

    # neoantigen files (for independent neoantigen-level Venn)
    neo = {}
    for grp in ("SBS2_HIGH", "CNV_HIGH"):
        neo[grp] = pd.read_csv(os.path.join(MHC_DIR, f"{grp}_neoantigens.tsv"), sep="\t")
    st["neo"] = neo

    # fusion cross-reference (optional)
    st["fusion"] = None
    for name in ("cross_group_neoantigen_fusion_overlap.tsv", "neoantigen_fusion_crossref.tsv"):
        p = os.path.join(FUSION_DIR, name)
        if os.path.exists(p):
            st["fusion"] = pd.read_csv(p, sep="\t")
            st["fusion_file"] = name
            log(f"  fusion crossref: {name} {st['fusion'].shape}")
            break
    if st["fusion"] is None:
        log("  [WARN] no fusion cross-reference found; fusion_disrupted will be False for all.")
    log("  Beat 0 passed.")
    return st


# =============================================================================
# PART 1: NEOANTIGEN-LEVEL TIERS + INDEPENDENT VENN CONFIRMATION
# =============================================================================
def part1_tiers(st):
    banner("PART 1: NEOANTIGEN-LEVEL TIERS (from full.tsv) + Venn cross-check")
    full = st["full"]
    counts = full["tier"].value_counts().to_dict()
    log("  From full.tsv:")
    for t in ("SBS2_specific", "shared", "CNV_specific"):
        log(f"    {t:15s} {counts.get(t, 0)}   (expected {EXPECTED_TIERS[t]})")

    # independent neoantigen-level overlap from the neoantigen files
    def muts(grp):
        df = st["neo"][grp]
        g = _first(df.columns, ["gene"]); h = _first(df.columns, ["hgvs_p"])
        return set(zip(df[g].astype(str), df[h].astype(str)))
    s, c = muts("SBS2_HIGH"), muts("CNV_HIGH")
    venn = {"SBS2_specific": len(s - c), "shared": len(s & c), "CNV_specific": len(c - s)}
    log("  Independent Venn from neoantigen files (unique gene+hgvs_p):")
    for t in ("SBS2_specific", "shared", "CNV_specific"):
        agree = "ok" if venn[t] == counts.get(t, -1) else "[!] disagrees with full.tsv"
        log(f"    {t:15s} {venn[t]}   {agree}")
    if venn != {t: counts.get(t, -1) for t in venn}:
        log("  [WARN] ranking tiers and neoantigen-file Venn disagree. ref_tri_fasta "
            "in_sbs2/in_cnv may predate the current neoantigen files; reconcile before "
            "quoting the overlap.")
    st["venn"] = venn
    return st


# =============================================================================
# PART 2: PER-GENE AGGREGATES (folded in from retired Step04)
# =============================================================================
def part2_aggregates(st):
    banner("PART 2: PER-GENE AGGREGATES (neoantigen + strong-binder counts)")
    rows = []
    for grp in ("SBS2_HIGH", "CNV_HIGH"):
        df = st["pep"][grp]
        binder = df[df["mut_ic50"] < MHC_BIND]
        n_neo = binder.groupby("gene").size().rename("n_neoantigen_peptides")
        if "is_strong_binder" in df.columns:
            strong = df[(df["mut_ic50"] < 50)].groupby("gene").size().rename("n_strong_binders")
        else:
            strong = df[df["mut_ic50"] < 50].groupby("gene").size().rename("n_strong_binders")
        agg = pd.concat([n_neo, strong], axis=1).fillna(0).astype(int).reset_index()
        agg["group"] = grp
        rows.append(agg)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(os.path.join(OUTPUT_DIR, "per_gene_neoantigen_aggregates.tsv"), sep="\t", index=False)
    log(f"  Wrote per-gene aggregates: {out.shape[0]} gene-group rows")
    st["aggregates"] = out
    return st


# =============================================================================
# PART 3: GENE-LEVEL ESCAPE MECHANISMS
# =============================================================================
def part3_mechanisms(st):
    banner("PART 3: ESCAPE MECHANISMS (fusion + silencing)")

    # fusion disruption: CNV fusion removes an SBS2 neoantigen
    fusion_genes = set()
    if st["fusion"] is not None:
        f = st["fusion"]
        gcol = _first(f.columns, ["gene"])
        dcol = _first(f.columns, ["direction"])
        if gcol and dcol:
            sel = f[f[dcol].astype(str).str.contains("CNV", case=False, na=False)]
            fusion_genes = set(sel[gcol].astype(str))
        elif gcol:
            fusion_genes = set(f[gcol].astype(str))
    log(f"  fusion_disrupted genes: {len(fusion_genes)}")

    # silencing: mean expression CNV < 0.5x SBS2 (needs adata)
    silenced_genes = set()
    st["silence_ok"] = False
    if os.path.exists(ADATA_PATH):
        try:
            import anndata as ad
            adata = ad.read_h5ad(ADATA_PATH)
            obs = {str(b): i for i, b in enumerate(adata.obs_names)}
            var = {str(g): j for j, g in enumerate(adata.var_names)}
            neo_genes = set(st["full"]["gene"].astype(str))
            def mean_expr(cells):
                idx = [obs[b] for b in cells if b in obs]
                if not idx: return {}
                X = adata.X[idx, :]
                means = np.asarray(X.mean(axis=0)).ravel()
                return means
            m_sbs2 = mean_expr(st["sbs2_cells"])
            m_cnv = mean_expr(st["cnv_cells"])
            if len(m_sbs2) and len(m_cnv):
                for g in neo_genes:
                    j = var.get(g)
                    if j is None: continue
                    es, ec = float(m_sbs2[j]), float(m_cnv[j])
                    if es > SILENCE_FLOOR and ec < SILENCE_RATIO * es:
                        silenced_genes.add(g)
                st["silence_ok"] = True
                log(f"  silenced genes (CNV < 0.5x SBS2, SBS2 > {SILENCE_FLOOR}): {len(silenced_genes)}")
            else:
                log("  [WARN] no adata cells matched group barcodes; silencing skipped.")
        except Exception as e:
            log(f"  [WARN] adata/scanpy unavailable ({e}); silencing skipped.")
    else:
        log(f"  [WARN] adata not found ({ADATA_PATH}); silencing skipped.")

    st["fusion_genes"] = fusion_genes
    st["silenced_genes"] = silenced_genes
    return st


# =============================================================================
# PART 4: SHARED-TIER SELECTION-EVIDENCE ANALYSIS
# =============================================================================
def part4_shared_evidence(st):
    banner("PART 4: SHARED-TIER SELECTION EVIDENCE (the 93 shared neoantigens)")
    full = st["full"]
    shared = full[full["tier"] == "shared"].copy()
    log(f"  Shared neoantigens: {len(shared)} (expected {EXPECTED_TIERS['shared']})")

    g = shared["gene"].astype(str)
    shared["is_hla"] = g.str.startswith(HLA_PREFIXES)
    shared["fusion_disrupted"] = g.isin(st["fusion_genes"])
    shared["silenced"] = g.isin(st["silenced_genes"])
    shared["mechanism_evidence"] = shared["fusion_disrupted"] | shared["silenced"]
    # HLA is excluded from the count (artifact-prone); flagged, not dropped
    shared["selection_evidence"] = shared["mechanism_evidence"] & (~shared["is_hla"])

    def mechs(r):
        m = []
        if r["fusion_disrupted"]: m.append("fusion")
        if r["silenced"]: m.append("silenced")
        return "; ".join(m)
    shared["mechanisms"] = shared.apply(mechs, axis=1)

    # descriptive prevalence direction (NOT a gate)
    def direction(r):
        if r["prevalence_cnv"] < r["prevalence_sbs2"]: return "depleted_in_CNV"
        if r["prevalence_cnv"] > r["prevalence_sbs2"]: return "enriched_in_CNV"
        return "flat"
    shared["cnv_direction"] = shared.apply(direction, axis=1)

    cols = ["gene", "hgvs_p", "peptide", "is_hla", "prevalence_sbs2", "prevalence_cnv",
            "cnv_direction", "prevalence_tier", "fusion_disrupted", "silenced",
            "mechanism_evidence", "selection_evidence", "mechanisms"]
    cols = [c for c in cols if c in shared.columns]
    out = shared[cols].sort_values(
        ["selection_evidence", "prevalence_tier"], ascending=[False, False])
    out.to_csv(os.path.join(OUTPUT_DIR, "shared_neoantigen_selection_evidence.tsv"),
               sep="\t", index=False)

    n_counted = int(shared["selection_evidence"].sum())
    n_hla_ev = int((shared["mechanism_evidence"] & shared["is_hla"]).sum())
    log(f"  Shared neoantigens with selection evidence (fusion or silencing, "
        f"HLA excluded): {n_counted} / {len(shared)}")
    if not st["silence_ok"]:
        log("  NOTE: silencing was skipped (no adata), so this count is a floor.")
    log(f"    by mechanism (counted set): fusion "
        f"{int((shared['fusion_disrupted'] & ~shared['is_hla']).sum())}, "
        f"silenced {int((shared['silenced'] & ~shared['is_hla']).sum())}")
    log(f"    HLA neoantigens with a mechanism, EXCLUDED from the count: {n_hla_ev}")

    # direction breakdown of the counted set (keeps the 'lost in CNV' wording honest)
    counted = shared[shared["selection_evidence"]]
    dc = counted["cnv_direction"].value_counts().to_dict()
    log(f"    counted set direction: depleted {dc.get('depleted_in_CNV', 0)}, "
        f"flat {dc.get('flat', 0)}, enriched {dc.get('enriched_in_CNV', 0)}")

    log("\n  Shared neoantigens carrying selection evidence (HLA excluded):")
    for r in out[out["selection_evidence"]].itertuples():
        log(f"    {r.gene:12s} {r.hgvs_p:14s}  prevSBS2={r.prevalence_sbs2:.3f} "
            f"prevCNV={r.prevalence_cnv:.3f}  {r.cnv_direction:16s} [{r.mechanisms}]")

    if n_hla_ev:
        log("\n  HLA shared neoantigens with a mechanism (flagged, excluded from count):")
        for r in out[(out["is_hla"]) & (out["mechanism_evidence"])].itertuples():
            log(f"    {r.gene:12s} {r.hgvs_p:14s}  prevSBS2={r.prevalence_sbs2:.3f} "
                f"prevCNV={r.prevalence_cnv:.3f}  {r.cnv_direction:16s} [{r.mechanisms}]")

    # featured shared gene check
    feat = out[(out["gene"] == "SPRR1A") & (out["hgvs_p"] == "p.Val61Ile")]
    if len(feat):
        r = feat.iloc[0]
        log(f"\n  Featured shared candidate SPRR1A p.Val61Ile: "
            f"evidence={bool(r['selection_evidence'])} [{r['mechanisms'] or 'none'}]")
    else:
        log("\n  SPRR1A p.Val61Ile not in the shared set (check tier assignment).")
    st["shared_evidence"] = out
    return st


# =============================================================================
# PART 5: ANNOTATED TARGET TABLE (prevalence-ordered, escape flags joined)
# =============================================================================
def part5_targets(st):
    banner("PART 5: ANNOTATED NEOANTIGEN TARGET TABLE (prevalence-ordered)")
    full = st["full"].copy()
    g = full["gene"].astype(str)
    full["is_hla"] = g.str.startswith(HLA_PREFIXES)
    full["fusion_disrupted"] = g.isin(st["fusion_genes"])
    full["silenced"] = g.isin(st["silenced_genes"])
    full["mechanism_evidence"] = full["fusion_disrupted"] | full["silenced"]
    full["selection_evidence"] = full["mechanism_evidence"] & (~full["is_hla"])
    full = full.sort_values("prevalence_tier", ascending=False)
    full.to_csv(os.path.join(OUTPUT_DIR, "neoantigen_targets_annotated.tsv"), sep="\t", index=False)
    log(f"  Wrote annotated target table: {full.shape} (all neoantigens, prevalence-ordered)")
    log("  No composite score; ordering is prevalence_tier. Escape flags are annotation only.")
    return st


# =============================================================================
# MAIN
# =============================================================================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    banner("INTEGRATED NEOANTIGEN ANALYSIS (prevalence-first rewrite)")
    st = beat0()
    st = part1_tiers(st)
    st = part2_aggregates(st)
    st = part3_mechanisms(st)
    st = part4_shared_evidence(st)
    st = part5_targets(st)
    banner("COMPLETE")
    rp = os.path.join(OUTPUT_DIR, "step06_integrated_report.txt")
    with open(rp, "w") as f:
        f.write("\n".join(_report))
    log(f"  Report: {rp}")


if __name__ == "__main__":
    main()
