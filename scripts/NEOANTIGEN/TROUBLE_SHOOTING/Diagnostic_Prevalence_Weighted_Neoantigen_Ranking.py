#!/usr/bin/env python3
"""
Diagnostic_Prevalence_Weighted_Neoantigen_Ranking_v2.py
=======================================================
Read-only diagnostic. Ranks SBS2 + CNV differential neoantigens by CLONAL
PREVALENCE (primary) and MHC-I BINDING GAIN (secondary), with tier-aware
denominators so shared (Tier 1) hits are neither inflated by pool size nor
penalized for being spread across two niches.

What changed from v1 (this revision)
-------------------------------------
1. Candidate universe = BOTH SBS2_HIGH and CNV_HIGH neoantigen files, unioned
   and collapsed by (gene, hgvs_p). v1 loaded SBS2 only.
2. Tier per mutation from ref_tri_fasta.tsv (in_sbs2 / in_cnv), matching the
   manuscript Venn:
     shared (Tier 1) = in_sbs2 & in_cnv
     SBS2-specific   = in_sbs2 & ~in_cnv
     CNV-specific    = ~in_sbs2 & in_cnv
3. Carriers counted in BOTH groups (SBS2 546, CNV 546) for every mutation.
4. Denominators:
     prevalence_sbs2 = carriers_sbs2 / 546            (constant per-niche)
     prevalence_cnv  = carriers_cnv  / 546            (constant per-niche)
     prevalence_tier = tier-conditional HEADLINE:
                         SBS2-specific -> carriers_sbs2 / 546
                         CNV-specific  -> carriers_cnv  / 546
                         shared        -> (sbs2+cnv) / 1092
     prevalence_max  = max(prevalence_sbs2, prevalence_cnv)   (FAIRNESS column:
                         credits a shared hit for its strongest niche so it is
                         not diluted by a weak second niche; for two equal
                         groups mean/median == pooled, so max is the only
                         column that avoids the Tier-1 penalty)
5. Gene-union carrier counts per group + combined, to reconcile against the
   session-notes exemplars (KLF3 1, CAST 6, SERPINB2 42, KRT6B 51, ANXA1 6).
   The guardrail now compares gene-union to gene-union, not per-mutation.
6. NORMAL: nothing added. Step02 already drops any variant present in NORMAL
   before binding prediction, so every candidate here is already somatic
   (NORMAL-subtracted) by construction.

Carrier match: (chrom, pos_from_neo_location, Base_observed == alt_from_reftri).
Nucleotide alt comes from ref_tri via (gene, hgvs_p); the position used for the
genotype-master match is the neoantigen 'location' value (the convention that
matched 560/560 in v1), so no coordinate-system risk is introduced.

Binding gain: delta_IC50 = wt_IC50 - mut_IC50 (largest positive first). Each
mutation represented by its max-delta peptide; best_mut_IC50_any_allele carried
separately. Negatives retained and sorted to the bottom. No prevalence threshold.

Run in NETWORK conda env (needs scanpy/anndata for expression).
    conda run -n NETWORK python Diagnostic_Prevalence_Weighted_Neoantigen_Ranking_v2.py
"""

import os
import re
import sys
import pandas as pd
import numpy as np

# =============================================================================
# CONFIG
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"

NEO_SBS2_PATH = os.path.join(BASE_DIR, "data/FIG_7/03_mhc_binding/SBS2_HIGH_neoantigens.tsv")
NEO_CNV_PATH  = os.path.join(BASE_DIR, "data/FIG_7/03_mhc_binding/CNV_HIGH_neoantigens.tsv")
REFTRI_PATH   = os.path.join(BASE_DIR, "data/FIG_7/fasta_context/ref_tri_fasta.tsv")
GROUP_PATH    = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
ADATA_PATH    = os.path.join(BASE_DIR, "data/FIG_4/00_input/adata_final.h5ad")

# Corrected genotype-master location (from Jake).
GENO_PATH = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
             "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

OUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/06_prevalence_ranking")

SBS2_LABEL = "SBS2_HIGH"
CNV_LABEL  = "CNV_HIGH"
MHC_BIND   = 500.0   # nM; mutant binder threshold defining the candidate set

# Exemplar GENE-UNION carrier cross-check (session notes). Tier per notes.
EXEMPLAR_EXPECTED = {
    "KLF3": (1, "shared"), "CAST": (6, "SBS2_specific"),
    "SERPINB2": (42, "shared"), "KRT6B": (51, "SBS2_specific"),
    "ANXA1": (6, "SBS2_specific"),
}

# =============================================================================
# LOGGING
# =============================================================================
_report = []

def log(msg=""):
    print(msg, flush=True)
    _report.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)

def _first_present(cols, candidates):
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand in cols:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None

def _parse_location(series):
    chrom, pos = [], []
    for v in series.astype(str):
        m = re.match(r"^(chr[\w]+|[\dXYMT]+)[:_\-](\d+)", v)
        if m:
            chrom.append(m.group(1)); pos.append(int(m.group(2)))
        else:
            chrom.append(None); pos.append(np.nan)
    return pd.Series(chrom), pd.Series(pos)

def _as_bool(series):
    """Coerce a possibly-string/int column to boolean."""
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(
        ["true", "1", "yes", "t", "y"])


# =============================================================================
# BEAT 0: LOAD + SCHEMA VERIFICATION (both files, both groups)
# =============================================================================
def beat0_verify():
    banner("BEAT 0: SCHEMA + JOIN VERIFICATION")
    state = {}

    for label, path in [
        ("neo_SBS2", NEO_SBS2_PATH), ("neo_CNV", NEO_CNV_PATH),
        ("ref_tri_fasta", REFTRI_PATH), ("genotype_master", GENO_PATH),
        ("group_assignments", GROUP_PATH),
    ]:
        if not os.path.exists(path):
            log(f"  [FATAL] missing {label}: {path}")
            sys.exit(1)
    log("  All required tables present.")

    # --- Group assignments: SBS2 + CNV cell sets ---
    groups = pd.read_csv(GROUP_PATH, sep="\t")
    bc_col  = _first_present(groups.columns, ["cell_barcode", "CB", "barcode"])
    grp_col = _first_present(groups.columns, ["group", "assignment", "label"])
    if bc_col is None or grp_col is None:
        log(f"  [FATAL] group file columns unexpected: {list(groups.columns)}")
        sys.exit(1)
    sbs2_bc = set(groups.loc[groups[grp_col] == SBS2_LABEL, bc_col].astype(str))
    cnv_bc  = set(groups.loc[groups[grp_col] == CNV_LABEL,  bc_col].astype(str))
    log(f"  Group file: {groups.shape}")
    log(f"  {SBS2_LABEL}: {len(sbs2_bc)} cells   {CNV_LABEL}: {len(cnv_bc)} cells")
    if len(sbs2_bc) == 0 or len(cnv_bc) == 0:
        log("  [FATAL] one of the tumor groups is empty. Check group labels.")
        sys.exit(1)
    state.update(sbs2_bc=sbs2_bc, cnv_bc=cnv_bc,
                 n_sbs2=len(sbs2_bc), n_cnv=len(cnv_bc))

    # --- Genotype master ---
    geno = pd.read_csv(GENO_PATH, sep="\t")
    g_chrom = _first_present(geno.columns, ["#CHROM", "CHROM", "chrom", "chr"])
    g_pos   = _first_present(geno.columns, ["Start", "POS", "pos", "scstart"])
    g_alt   = _first_present(geno.columns, ["Base_observed", "ALT", "alt"])
    g_cb    = _first_present(geno.columns, ["CB", "cell_barcode", "barcode"])
    log(f"\n  Genotype master: {geno.shape}")
    log(f"  Mapped -> chrom='{g_chrom}', pos='{g_pos}', alt='{g_alt}', CB='{g_cb}'")
    if None in (g_chrom, g_pos, g_cb):
        log("  [FATAL] genotype master missing chrom/pos/CB.")
        sys.exit(1)
    geno_cb = set(geno[g_cb].astype(str))
    ov_sbs2 = len(sbs2_bc & geno_cb); ov_cnv = len(cnv_bc & geno_cb)
    log(f"  Barcode overlap: SBS2 {ov_sbs2}/{len(sbs2_bc)}, "
        f"CNV {ov_cnv}/{len(cnv_bc)}")
    if ov_sbs2 == 0 or ov_cnv == 0:
        log("  [FATAL] zero barcode overlap for a tumor group. Reconcile the "
            "'-1-{sample_id}' CB suffix before trusting carrier counts.")
        sys.exit(1)
    state["geno"] = geno
    state["g_cols"] = dict(chrom=g_chrom, pos=g_pos, alt=g_alt, cb=g_cb)

    # --- Neoantigen files (both) ---
    def load_neo(path, tag):
        df = pd.read_csv(path, sep="\t")
        cols = dict(
            gene=_first_present(df.columns, ["gene", "Gene", "symbol"]),
            hgvs=_first_present(df.columns, ["hgvs_p", "HGVSp"]),
            mut=_first_present(df.columns, ["mut_ic50", "mut_IC50"]),
            wt=_first_present(df.columns, ["wt_ic50", "wt_IC50"]),
            pep=_first_present(df.columns, ["mut_peptide", "peptide"]),
            allele=_first_present(df.columns, ["best_allele", "allele"]),
            diff=_first_present(df.columns, ["is_differential"]),
            loc=_first_present(df.columns, ["location", "loc"]),
        )
        for req in ("gene", "hgvs", "mut", "wt"):
            if cols[req] is None:
                log(f"  [FATAL] {tag} missing required column '{req}'")
                sys.exit(1)
        if cols["loc"] is None:
            log(f"  [FATAL] {tag} has no 'location' column for locus matching.")
            sys.exit(1)
        df = df.rename(columns={
            cols["gene"]: "gene", cols["hgvs"]: "hgvs_p",
            cols["mut"]: "mut_IC50", cols["wt"]: "wt_IC50",
            cols["loc"]: "location"})
        if cols["pep"]:    df = df.rename(columns={cols["pep"]: "peptide"})
        if cols["allele"]: df = df.rename(columns={cols["allele"]: "allele"})
        if cols["diff"]:   df = df.rename(columns={cols["diff"]: "is_differential"})
        df["_src"] = tag
        return df

    neo_sbs2 = load_neo(NEO_SBS2_PATH, "SBS2")
    neo_cnv  = load_neo(NEO_CNV_PATH,  "CNV")
    log(f"\n  Neoantigen rows: SBS2 {len(neo_sbs2)}, CNV {len(neo_cnv)}")
    neo = pd.concat([neo_sbs2, neo_cnv], ignore_index=True)
    log(f"  Combined peptide rows: {len(neo)}")
    log(f"  Unique mutations across both files: "
        f"{neo[['gene','hgvs_p']].drop_duplicates().shape[0]}")
    state["neo"] = neo

    # --- ref_tri_fasta: tier source, TCW, nucleotide alt, locus ---
    reftri = pd.read_csv(REFTRI_PATH, sep="\t")
    r = dict(
        gene=_first_present(reftri.columns, ["gene", "symbol"]),
        hgvs=_first_present(reftri.columns, ["hgvs_p", "HGVSp"]),
        in_sbs2=_first_present(reftri.columns, ["in_sbs2"]),
        in_cnv=_first_present(reftri.columns, ["in_cnv"]),
        tcw=_first_present(reftri.columns, ["is_tcw"]),
        tcwc=_first_present(reftri.columns, ["is_tcw_ct"]),
        sub=_first_present(reftri.columns, ["sub_pyr"]),
        alt=_first_present(reftri.columns, ["alt"]),
        chrom=_first_present(reftri.columns, ["chrom", "#CHROM", "chr"]),
        pos=_first_present(reftri.columns, ["pos", "Start", "POS"]),
    )
    log(f"\n  ref_tri_fasta: {reftri.shape}")
    log(f"  Mapped -> gene='{r['gene']}', hgvs_p='{r['hgvs']}', "
        f"in_sbs2='{r['in_sbs2']}', in_cnv='{r['in_cnv']}', alt='{r['alt']}'")
    for req in ("gene", "hgvs", "in_sbs2", "in_cnv", "alt"):
        if r[req] is None:
            log(f"  [FATAL] ref_tri_fasta missing '{req}' (needed for tier/alt).")
            sys.exit(1)
    state["reftri"] = reftri
    state["r_cols"] = r

    log("\n  Beat 0 verification passed.")
    return state


# =============================================================================
# BEAT 1: BINDING GAIN over the union
# =============================================================================
def beat1_binding(state):
    banner("BEAT 1: BINDING GAIN (delta_IC50 = wt_IC50 - mut_IC50), union")
    neo = state["neo"].copy()
    neo["mut_IC50"] = pd.to_numeric(neo["mut_IC50"], errors="coerce")
    neo["wt_IC50"]  = pd.to_numeric(neo["wt_IC50"], errors="coerce")
    neo = neo.dropna(subset=["mut_IC50", "wt_IC50"])
    n0 = len(neo)
    neo = neo[neo["mut_IC50"] < MHC_BIND].copy()
    neo["delta_IC50"] = neo["wt_IC50"] - neo["mut_IC50"]
    neo["fold_improvement"] = neo["wt_IC50"] / neo["mut_IC50"].replace(0, np.nan)
    if "is_differential" in neo.columns:
        neo["is_differential"] = _as_bool(neo["is_differential"])
    else:
        neo["is_differential"] = neo["wt_IC50"] > MHC_BIND
    log(f"  Peptide rows: {n0} -> mutant binders (<{MHC_BIND:.0f} nM): {len(neo)}")
    log(f"  Negative-delta rows (mut worse than wt): {(neo['delta_IC50'] < 0).sum()} "
        f"(retained, sink to bottom)")
    state["binders"] = neo
    return state


# =============================================================================
# BEAT 2: COLLAPSE to one row per (gene, hgvs_p) by MAX delta_IC50
# =============================================================================
def beat2_collapse(state):
    banner("BEAT 2: COLLAPSE TO ONE ROW PER MUTATION (max delta_IC50)")
    neo = state["binders"]
    key = ["gene", "hgvs_p"]
    best_abs = neo.groupby(key)["mut_IC50"].min().rename("best_mut_IC50_any_allele")
    rep = neo.loc[neo.groupby(key)["delta_IC50"].idxmax()].copy()
    rep = rep.merge(best_abs, left_on=key, right_index=True, how="left")
    log(f"  Unique mutations: {len(rep)}")
    log(f"  Differential among them: {int(rep['is_differential'].sum())}")
    state["rep"] = rep
    state["long"] = neo
    return state


# =============================================================================
# BEAT 3: JOIN ref_tri -> tier (Venn), TCW, nucleotide alt
# =============================================================================
def beat3_tier_tcw(state):
    banner("BEAT 3: TIER (Venn) + TCW + nucleotide alt from ref_tri")
    rep = state["rep"]; reftri = state["reftri"]; r = state["r_cols"]

    cols = [r["gene"], r["hgvs"], r["in_sbs2"], r["in_cnv"], r["alt"]]
    for k in ("tcw", "tcwc", "sub"):
        if r[k]: cols.append(r[k])
    sub = reftri[cols].drop_duplicates(subset=[r["gene"], r["hgvs"]])
    rep = rep.merge(sub, left_on=["gene", "hgvs_p"],
                    right_on=[r["gene"], r["hgvs"]], how="left", suffixes=("", "_rt"))

    ren = {r["in_sbs2"]: "in_sbs2", r["in_cnv"]: "in_cnv", r["alt"]: "alt_nt"}
    if r["tcw"]:  ren[r["tcw"]]  = "is_tcw"
    if r["tcwc"]: ren[r["tcwc"]] = "is_tcw_ct"
    if r["sub"]:  ren[r["sub"]]  = "sub_pyr"
    rep = rep.rename(columns=ren)

    for c in ["is_tcw", "is_tcw_ct"]:
        if c in rep.columns:
            rep[c] = _as_bool(rep[c])
        else:
            rep[c] = False
    rep["in_sbs2"] = _as_bool(rep["in_sbs2"])
    rep["in_cnv"]  = _as_bool(rep["in_cnv"])

    def tier(row):
        a, b = row["in_sbs2"], row["in_cnv"]
        if a and b:      return "shared"
        if a and not b:  return "SBS2_specific"
        if b and not a:  return "CNV_specific"
        return "unassigned"
    rep["tier"] = rep.apply(tier, axis=1)

    miss = rep["alt_nt"].isna().sum()
    log(f"  Mutations joined to ref_tri: {len(rep) - miss}/{len(rep)} "
        f"(unjoined: {miss})")
    log(f"  Tier counts: " + ", ".join(
        f"{t}={int((rep['tier'] == t).sum())}"
        for t in ["shared", "SBS2_specific", "CNV_specific", "unassigned"]))
    log(f"  TCW-positive (all TCW): {int(rep['is_tcw'].sum())}; "
        f"clean C>T: {int(rep['is_tcw_ct'].sum())}")
    state["rep"] = rep
    return state


# =============================================================================
# BEAT 4: CARRIERS in both groups, tier-conditional prevalence, max-niche,
#         gene-union, consistency flag
# =============================================================================
def beat4_prevalence(state):
    banner("BEAT 4: CARRIERS + TIER-CONDITIONAL PREVALENCE")
    rep = state["rep"]; geno = state["geno"]; g = state["g_cols"]
    n_sbs2, n_cnv = state["n_sbs2"], state["n_cnv"]

    geno = geno.copy()
    geno["_chrom"] = geno[g["chrom"]].astype(str)
    geno["_pos"] = pd.to_numeric(geno[g["pos"]], errors="coerce")
    geno["_alt"] = geno[g["alt"]].astype(str) if g["alt"] else ""
    geno["_cb"] = geno[g["cb"]].astype(str)

    def carrier_map_for(cellset):
        sub = geno[geno["_cb"].isin(cellset)]
        m = {}
        for (ch, ps, al), cb in sub.groupby(["_chrom", "_pos", "_alt"])["_cb"]:
            m[(ch, ps, al)] = set(cb)
        return m

    map_sbs2 = carrier_map_for(state["sbs2_bc"])
    map_cnv  = carrier_map_for(state["cnv_bc"])
    log(f"  Genotype rows -> SBS2 map keys: {len(map_sbs2)}, "
        f"CNV map keys: {len(map_cnv)}")

    ch_arr, ps_arr = _parse_location(rep["location"])
    rep = rep.assign(_chrom=ch_arr.values, _pos=ps_arr.values)

    def carriers(cmap, ch, ps, alt):
        if ch is None or (isinstance(ps, float) and np.isnan(ps)):
            return set()
        alt = str(alt)
        s = cmap.get((str(ch), float(ps), alt))
        if s is None and not (isinstance(ps, float) and np.isnan(ps)):
            s = cmap.get((str(ch), float(int(ps)), alt))
        return s or set()

    cs_sbs2, cs_cnv = [], []
    for _, row in rep.iterrows():
        cs_sbs2.append(carriers(map_sbs2, row["_chrom"], row["_pos"], row["alt_nt"]))
        cs_cnv.append(carriers(map_cnv,  row["_chrom"], row["_pos"], row["alt_nt"]))
    rep["_carr_sbs2_set"] = cs_sbs2
    rep["_carr_cnv_set"]  = cs_cnv
    rep["carriers_sbs2"] = [len(s) for s in cs_sbs2]
    rep["carriers_cnv"]  = [len(s) for s in cs_cnv]

    rep["prevalence_sbs2"] = rep["carriers_sbs2"] / n_sbs2
    rep["prevalence_cnv"]  = rep["carriers_cnv"] / n_cnv

    def tier_prev(row):
        if row["tier"] == "SBS2_specific": return row["carriers_sbs2"] / n_sbs2
        if row["tier"] == "CNV_specific":  return row["carriers_cnv"] / n_cnv
        if row["tier"] == "shared":
            return (row["carriers_sbs2"] + row["carriers_cnv"]) / (n_sbs2 + n_cnv)
        return np.nan
    rep["prevalence_tier"] = rep.apply(tier_prev, axis=1)
    rep["prevalence_max"] = rep[["prevalence_sbs2", "prevalence_cnv"]].max(axis=1)

    # Consistency: Venn-specific but carriers in the other group (eyeball flag).
    def consistency(row):
        if row["tier"] == "SBS2_specific" and row["carriers_cnv"] > 0:
            return f"SBS2-specific but {row['carriers_cnv']} CNV carriers"
        if row["tier"] == "CNV_specific" and row["carriers_sbs2"] > 0:
            return f"CNV-specific but {row['carriers_sbs2']} SBS2 carriers"
        return ""
    rep["tier_carrier_note"] = rep.apply(consistency, axis=1)
    n_flag = (rep["tier_carrier_note"] != "").sum()

    # Gene-union carriers (union of cells over the gene's candidate mutations).
    gu_sbs2, gu_cnv, gu_comb = {}, {}, {}
    for gene, grp in rep.groupby("gene"):
        s2 = set().union(*grp["_carr_sbs2_set"]) if len(grp) else set()
        cv = set().union(*grp["_carr_cnv_set"]) if len(grp) else set()
        gu_sbs2[gene] = len(s2); gu_cnv[gene] = len(cv)
        gu_comb[gene] = len(s2 | cv)
    rep["gene_union_sbs2"] = rep["gene"].map(gu_sbs2)
    rep["gene_union_cnv"]  = rep["gene"].map(gu_cnv)
    rep["gene_union_comb"] = rep["gene"].map(gu_comb)

    n_zero = int(((rep["carriers_sbs2"] + rep["carriers_cnv"]) == 0).sum())
    log(f"  Mutations with 0 carriers in either group: {n_zero}")
    log(f"  Tier/carrier consistency flags: {n_flag} "
        f"(informational: Venn membership is neoantigen-calling, carriers is "
        f"variant presence; some cross-group carriers are expected)")
    state["rep"] = rep
    return state


# =============================================================================
# BEAT 5: EXPRESSION per niche (both groups over their own 546)
# =============================================================================
def beat5_expression(state):
    banner("BEAT 5: EXPRESSION PER NICHE (pct expressing over each group's 546)")
    rep = state["rep"]
    if not os.path.exists(ADATA_PATH):
        log(f"  [WARN] adata missing; expression columns NaN.")
        rep["pct_expressing_sbs2"] = np.nan
        rep["pct_expressing_cnv"] = np.nan
        state["rep"] = rep; return state
    try:
        import scanpy as sc
        import scipy.sparse as sp
    except Exception as e:
        log(f"  [WARN] scanpy/scipy import failed ({e}); expression NaN. Use NETWORK env.")
        rep["pct_expressing_sbs2"] = np.nan
        rep["pct_expressing_cnv"] = np.nan
        state["rep"] = rep; return state

    adata = sc.read_h5ad(ADATA_PATH)
    obs = {str(b): i for i, b in enumerate(adata.obs_names)}
    var_index = {str(g): j for j, g in enumerate(adata.var_names)}
    X = adata.X
    is_sparse = sp.issparse(X)

    def pct_expr(cellset):
        idx = [obs[b] for b in cellset if b in obs]
        if not idx:
            return {}, 0
        sub = X[idx, :]
        out = {}
        return sub, len(idx)

    def compute(cellset, genes):
        idx = [obs[b] for b in cellset if b in obs]
        n = len(idx)
        if n == 0:
            return {gname: np.nan for gname in genes}, 0
        sub = X[idx, :]
        res = {}
        for gname in genes:
            j = var_index.get(gname)
            if j is None:
                res[gname] = np.nan; continue
            col = sub[:, j]
            nnz = col.getnnz() if is_sparse else int((np.asarray(col).ravel() > 0).sum())
            res[gname] = 100.0 * nnz / n
        return res, n

    genes = rep["gene"].astype(str).unique()
    e_sbs2, n1 = compute(state["sbs2_bc"], genes)
    e_cnv,  n2 = compute(state["cnv_bc"],  genes)
    rep["pct_expressing_sbs2"] = rep["gene"].astype(str).map(e_sbs2)
    rep["pct_expressing_cnv"]  = rep["gene"].astype(str).map(e_cnv)
    log(f"  Expression over SBS2 {n1} cells, CNV {n2} cells; "
        f"genes matched: {sum(v is not None for v in [var_index.get(g) for g in genes])}/{len(genes)}")
    state["rep"] = rep
    return state


# =============================================================================
# BEAT 6: RANK + ASSEMBLE + WRITE
# =============================================================================
def beat6_assemble(state):
    banner("BEAT 6: RANK + ASSEMBLE")
    rep = state["rep"].copy()

    gene_tcw = rep.groupby("gene")["is_tcw"].mean().mul(100).rename("gene_pct_TCW")
    rep = rep.merge(gene_tcw, on="gene", how="left")

    rep["prevalence_rank"] = rep["prevalence_tier"].rank(ascending=False, method="min")
    rep["delta_rank"] = rep["delta_IC50"].rank(ascending=False, method="min")
    rep["rank_product"] = rep["prevalence_rank"] * rep["delta_rank"]

    cols = ["gene", "hgvs_p", "tier", "peptide", "allele",
            "wt_IC50", "mut_IC50", "delta_IC50", "fold_improvement",
            "best_mut_IC50_any_allele",
            "is_tcw", "is_tcw_ct", "sub_pyr", "is_differential",
            "carriers_sbs2", "carriers_cnv",
            "prevalence_sbs2", "prevalence_cnv",
            "prevalence_tier", "prevalence_max",
            "gene_union_sbs2", "gene_union_cnv", "gene_union_comb",
            "pct_expressing_sbs2", "pct_expressing_cnv", "gene_pct_TCW",
            "tier_carrier_note",
            "prevalence_rank", "delta_rank", "rank_product"]
    cols = [c for c in cols if c in rep.columns]
    out = rep[cols].copy()
    out = out.sort_values(["prevalence_tier", "delta_IC50"],
                          ascending=[False, False]).reset_index(drop=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    full = os.path.join(OUT_DIR, "neoantigen_prevalence_ranking_full.tsv")
    lng  = os.path.join(OUT_DIR, "neoantigen_prevalence_ranking_long.tsv")
    out.to_csv(full, sep="\t", index=False)
    state["long"].to_csv(lng, sep="\t", index=False)
    log(f"  Wrote full ranking: {full}  ({out.shape})")
    log(f"  Wrote long companion: {lng}")
    state["out"] = out
    return state


# =============================================================================
# BEAT 7: REPORT + GUARDRAIL (gene-union reconciliation)
# =============================================================================
def beat7_report(state):
    banner("BEAT 7: RANKING REVIEW + GUARDRAIL")
    out = state["out"]

    def show(df, title, n=20):
        log(f"\n  {title} (top {n}):")
        h = ("  {:>3s}  {:10s}  {:13s}  {:13s}  {:>8s}  {:>9s}  {:>4s}  "
             "{:>4s}  {:>5s}  {:>8s}  {:>8s}")
        log(h.format("#", "gene", "hgvs_p", "tier", "mutIC50", "deltaIC50",
                     "cS", "cC", "TCW", "prevTier", "prevMax"))
        for i, row in df.head(n).iterrows():
            log("  {:>3d}  {:10s}  {:13s}  {:13s}  {:>8.1f}  {:>9.1f}  "
                "{:>4.0f}  {:>4.0f}  {:>5s}  {:>8.4f}  {:>8.4f}".format(
                    i + 1, str(row["gene"])[:10], str(row["hgvs_p"])[:13],
                    str(row["tier"])[:13], row.get("mut_IC50", np.nan),
                    row.get("delta_IC50", np.nan),
                    row.get("carriers_sbs2", 0), row.get("carriers_cnv", 0),
                    "Y" if bool(row.get("is_tcw", False)) else "-",
                    row.get("prevalence_tier", np.nan),
                    row.get("prevalence_max", np.nan)))

    show(out, "PRIMARY: prevalence_tier desc, then delta_IC50 desc")
    show(out[out["is_tcw"]].reset_index(drop=True),
         "TCW-ONLY (all TCW: SBS2 + SBS13), same sort")
    show(out.sort_values("prevalence_max", ascending=False).reset_index(drop=True),
         "FAIRNESS VIEW: sorted by prevalence_max (niche-dominance)")

    banner("GUARDRAIL: gene-union carrier reconciliation", char="-")
    log("  Comparing gene-union carriers to session-notes exemplars.")
    log("  (specific -> its group's union; shared -> combined union)\n")
    log("  {:10s}  {:14s}  {:>6s}  {:>6s}  {:>6s}  {:>6s}  {:>4s}".format(
        "gene", "tier(notes)", "expect", "uS", "uC", "comb", "ok?"))
    for gene, (exp, tier_note) in EXEMPLAR_EXPECTED.items():
        sub = out[out["gene"] == gene]
        if len(sub) == 0:
            log(f"  {gene:10s}  {tier_note:14s}  {exp:>6d}   NOT IN RANKED SET  [!]")
            continue
        us = int(sub["gene_union_sbs2"].max())
        uc = int(sub["gene_union_cnv"].max())
        cb = int(sub["gene_union_comb"].max())
        best = min([us, uc, cb], key=lambda v: abs(v - exp))
        ok = "ok" if best == exp else "[!]"
        log(f"  {gene:10s}  {tier_note:14s}  {exp:>6d}  {us:>6d}  {uc:>6d}  "
            f"{cb:>6d}  {ok:>4s}")

    top10 = set(out.head(10)["gene"])
    for gene in ["KRT6B", "SERPINB2", "SPRR1A"]:
        log(f"  {gene}: {'in top 10 by prevalence_tier' if gene in top10 else 'not in top 10'}")


# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("PREVALENCE-WEIGHTED NEOANTIGEN RANKING v2 (read-only)")
    log(f"  BASE_DIR: {BASE_DIR}")
    log(f"  Candidate universe: SBS2 + CNV neoantigen files")
    log(f"  Headline: tier-conditional prevalence; fairness: prevalence_max")

    state = beat0_verify()
    state = beat1_binding(state)
    state = beat2_collapse(state)
    state = beat3_tier_tcw(state)
    state = beat4_prevalence(state)
    state = beat5_expression(state)
    state = beat6_assemble(state)
    beat7_report(state)

    banner("DIAGNOSTIC COMPLETE")
    os.makedirs(OUT_DIR, exist_ok=True)
    rp = os.path.join(OUT_DIR, "prevalence_ranking_report.txt")
    with open(rp, "w") as f:
        f.write("\n".join(_report))
    log(f"  Report written: {rp}")


if __name__ == "__main__":
    main()
