#!/usr/bin/env python3
"""
Diagnostic_TCGA_Enrichment_in_SC.py
=====================================

Answers the question: Are any TCGA communities enriched in the SC network
when we ignore SC community boundaries?

Tests:
  1. Per-TCGA-community enrichment in the FULL SC LCC (1,202 genes)
  2. Per-TCGA-community enrichment in SC C0 specifically (299 genes)
  3. Per-TCGA-community enrichment in SC C0+C5 (top 2 SC communities by overlap)
  4. Aggregate view: which TCGA community contributes the most SC overlap genes?
  5. What if we treated ALL SC community genes as one set?

Uses hypergeometric test with shared gene universe.

Usage:
  conda run -n NETWORK python Diagnostic_TCGA_Enrichment_in_SC.py
"""

import os, json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from datetime import datetime

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
SC_PARTITION = os.path.join(BASE_DIR, "data/FIG_4/04_communities/SC_best_partition.csv")
SC_SELECTED  = os.path.join(BASE_DIR, "data/FIG_4/02_differential_expression/SC_selected_genes_filtered.csv")
TCGA_PARTITION = os.path.join(BASE_DIR, "data/FIG_2/05_communities/TCGA-HNSC/TCGA-HNSC_best_partition.csv")
ENSG_TO_SYMBOL = os.path.join(BASE_DIR, "data/FIG_2/01_cleaned_expression/ensg_to_symbol.json")
HARRIS_ALL = os.path.join(BASE_DIR, "data/FIG_4/00_input/Harris_A3_interactors.txt")

A3_ENSG = {"ENSG00000128383":"APOBEC3A","ENSG00000179750":"APOBEC3B","ENSG00000244509":"APOBEC3C",
           "ENSG00000243811":"APOBEC3D","ENSG00000128394":"APOBEC3F","ENSG00000239713":"APOBEC3G",
           "ENSG00000100298":"APOBEC3H"}

def banner(t): print(f"\n{'='*80}\n{t}\n{'='*80}", flush=True)
def log(m): print(m, flush=True)

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float); n = len(p)
    if n == 0: return np.array([])
    o = np.argsort(p); r = np.empty(n); r[o] = np.arange(1,n+1)
    a = p*n/r; a = np.minimum.accumulate(a[np.argsort(r)[::-1]])[::-1]
    return np.clip(a, 0, 1)

def main():
    start = datetime.now()
    banner("TCGA COMMUNITY ENRICHMENT IN SC NETWORK")

    # Load ENSG→symbol
    with open(ENSG_TO_SYMBOL) as f: e2s = json.load(f)

    # Load TCGA communities → symbols
    tcga_df = pd.read_csv(TCGA_PARTITION)
    gc = "gene" if "gene" in tcga_df.columns else tcga_df.columns[0]
    tcga_comms = {}
    tcga_all_symbols = set()
    for c in sorted(tcga_df["community"].unique()):
        ensg_set = set(tcga_df[tcga_df["community"]==c][gc])
        symbols = set()
        for g in ensg_set:
            s = e2s.get(g, A3_ENSG.get(g, g))
            symbols.add(s)
        tcga_comms[c] = symbols
        tcga_all_symbols |= symbols

    # Load SC communities
    sc_df = pd.read_csv(SC_PARTITION)
    sc_comms = {}
    for c in sorted(sc_df["community"].unique()):
        sc_comms[c] = set(sc_df[sc_df["community"]==c]["gene"])
    sc_all_comm_genes = set()
    for v in sc_comms.values(): sc_all_comm_genes |= v

    # Load SC selected genes (the full 9,244 that entered the network)
    sc_selected = set(pd.read_csv(SC_SELECTED)["gene"])

    # Load Harris
    harris = set()
    if os.path.exists(HARRIS_ALL):
        with open(HARRIS_ALL) as f:
            for line in f:
                g = line.strip().split("\t")[0].strip()
                if g and not g.startswith("#"): harris.add(g)

    # =========================================================================
    # UNIVERSE DEFINITION
    # =========================================================================
    banner("1. UNIVERSE AND OVERLAP BASICS")

    # The proper universe is genes that COULD appear in both analyses
    # = genes present in the SC selected set (9,244) that also have
    #   a symbol mapping from the TCGA pipeline
    # This is more conservative than union and gives proper enrichment stats
    sc_in_tcga_space = sc_selected & tcga_all_symbols
    universe = sc_selected  # All SC network-ready genes
    univ_size = len(universe)

    log(f"  SC community genes: {len(sc_all_comm_genes)}")
    log(f"  SC selected (network-ready): {len(sc_selected)}")
    log(f"  TCGA community genes (symbols): {len(tcga_all_symbols)}")
    log(f"  SC selected ∩ TCGA symbols: {len(sc_in_tcga_space)}")
    log(f"  Universe for enrichment: {univ_size} (all SC selected genes)")

    # How many TCGA community genes are in the SC selected gene pool?
    for c in sorted(tcga_comms.keys()):
        in_sc = tcga_comms[c] & sc_selected
        in_sc_comm = tcga_comms[c] & sc_all_comm_genes
        log(f"  TCGA C{c} ({len(tcga_comms[c])}): {len(in_sc)} in SC selected, {len(in_sc_comm)} in SC communities")

    # =========================================================================
    # TEST 1: Each TCGA community enriched in FULL SC LCC?
    # =========================================================================
    banner("2. TCGA COMMUNITY ENRICHMENT IN FULL SC LCC (1,202 genes)")

    results = []
    for tc in sorted(tcga_comms.keys()):
        tc_genes = tcga_comms[tc]
        # N = TCGA community genes present in universe
        tc_in_univ = tc_genes & universe
        N = len(tc_in_univ)
        # n = SC LCC size
        n = len(sc_all_comm_genes)
        # k = overlap
        overlap = tc_in_univ & sc_all_comm_genes
        k = len(overlap)
        # M = universe
        M = univ_size

        if k > 0 and N > 0:
            pval = hypergeom.sf(k-1, M, N, n)
        else:
            pval = 1.0

        results.append({
            "tcga_community": tc,
            "tcga_size": len(tc_genes),
            "tcga_in_universe": N,
            "overlap_with_SC_LCC": k,
            "expected": N * n / M if M > 0 else 0,
            "fold_enrichment": (k / (N * n / M)) if (N * n / M) > 0 else 0,
            "p_value": pval,
            "overlap_genes": sorted(overlap)
        })

    res_df = pd.DataFrame(results)
    res_df["fdr"] = bh_fdr(res_df["p_value"].values)
    res_df = res_df.sort_values("p_value")

    log(f"\n  {'TCGA':>6} {'Size':>5} {'InUniv':>7} {'Overlap':>8} {'Expected':>9} {'Fold':>6} {'p-value':>10} {'FDR':>10}")
    log(f"  {'-'*70}")
    for _, r in res_df.iterrows():
        sig = "***" if r["fdr"] < 0.05 else "**" if r["fdr"] < 0.1 else "*" if r["fdr"] < 0.25 else ""
        log(f"  C{int(r['tcga_community']):>4} {int(r['tcga_size']):>5} {int(r['tcga_in_universe']):>7} "
            f"{int(r['overlap_with_SC_LCC']):>8} {r['expected']:>9.1f} {r['fold_enrichment']:>6.2f} "
            f"{r['p_value']:>10.4e} {r['fdr']:>10.4e} {sig}")

    log(f"\n  Top overlapping TCGA communities with gene lists:")
    for _, r in res_df.head(5).iterrows():
        genes = r["overlap_genes"]
        if genes:
            log(f"  TCGA C{int(r['tcga_community'])}: {', '.join(genes[:20])}")

    # =========================================================================
    # TEST 2: Each TCGA community enriched in SC C0 specifically?
    # =========================================================================
    banner("3. TCGA COMMUNITY ENRICHMENT IN SC C0 (299 genes)")

    sc_c0 = sc_comms.get(0, set())
    results_c0 = []
    for tc in sorted(tcga_comms.keys()):
        tc_in_univ = tcga_comms[tc] & universe
        N = len(tc_in_univ)
        n = len(sc_c0)
        overlap = tc_in_univ & sc_c0
        k = len(overlap)
        M = univ_size
        pval = hypergeom.sf(k-1, M, N, n) if k > 0 and N > 0 else 1.0
        results_c0.append({
            "tcga_community": tc, "tcga_size": len(tcga_comms[tc]),
            "overlap": k, "expected": N*n/M if M>0 else 0,
            "fold": (k/(N*n/M)) if (N*n/M)>0 else 0,
            "p_value": pval, "genes": sorted(overlap)
        })

    rc0 = pd.DataFrame(results_c0)
    rc0["fdr"] = bh_fdr(rc0["p_value"].values)
    rc0 = rc0.sort_values("p_value")

    log(f"\n  {'TCGA':>6} {'Size':>5} {'Overlap':>8} {'Expected':>9} {'Fold':>6} {'p-value':>10} {'FDR':>10}")
    log(f"  {'-'*60}")
    for _, r in rc0.iterrows():
        if r["overlap"] > 0:
            log(f"  C{int(r['tcga_community']):>4} {int(r['tcga_size']):>5} {int(r['overlap']):>8} "
                f"{r['expected']:>9.1f} {r['fold']:>6.2f} {r['p_value']:>10.4e} {r['fdr']:>10.4e}")
            log(f"    Genes: {', '.join(r['genes'])}")

    # =========================================================================
    # TEST 3: What if we ignore SC community boundaries entirely?
    # =========================================================================
    banner("4. PERSPECTIVE: SC LCC AS ONE SET vs TCGA")

    log(f"  SC LCC contains {len(sc_all_comm_genes)} genes")
    log(f"  TCGA communities contain {len(tcga_all_symbols)} genes")
    total_overlap = sc_all_comm_genes & tcga_all_symbols
    log(f"  Total shared genes: {len(total_overlap)}")
    log(f"  Shared genes: {', '.join(sorted(total_overlap))}")

    # Which TCGA community contributes the most?
    log(f"\n  Breakdown by TCGA community:")
    tcga_contrib = []
    for tc in sorted(tcga_comms.keys()):
        contrib = tcga_comms[tc] & total_overlap
        tcga_contrib.append((tc, len(contrib), sorted(contrib)))
    tcga_contrib.sort(key=lambda x: -x[1])
    for tc, n, genes in tcga_contrib:
        if n > 0:
            log(f"    TCGA C{tc}: {n} genes — {', '.join(genes[:15])}")

    # =========================================================================
    # TEST 4: Known A3-Interactors status
    # =========================================================================
    banner("5. KNOWN A3-INTERACTORS — CROSS-REFERENCE")

    harris_in_sc_selected = harris & sc_selected
    harris_in_sc_comm = harris & sc_all_comm_genes
    harris_in_tcga = harris & tcga_all_symbols

    log(f"  Harris interactors: {len(harris)}")
    log(f"  In SC selected genes: {len(harris_in_sc_selected)} ({sorted(harris_in_sc_selected)})")
    log(f"  In SC communities: {len(harris_in_sc_comm)} ({sorted(harris_in_sc_comm)})")
    log(f"  In TCGA communities: {len(harris_in_tcga)}")
    if harris_in_tcga:
        log(f"    TCGA community assignments:")
        for g in sorted(harris_in_tcga):
            for tc, genes in tcga_comms.items():
                if g in genes:
                    log(f"      {g} → TCGA C{tc}")

    # Shared between SC and TCGA?
    harris_both = harris_in_sc_comm & harris_in_tcga
    log(f"  In BOTH SC and TCGA communities: {len(harris_both)} ({sorted(harris_both)})")

    # =========================================================================
    # TEST 5: What does the SC network look like if we DON'T split into communities?
    # =========================================================================
    banner("6. COMMUNITY RESOLUTION QUESTION")

    log(f"  Current SC communities: 14 (resolution 1.0)")
    log(f"  Largest: C0 with {len(sc_comms.get(0, set()))} genes")
    log(f"  If we merged all SC communities into 1:")
    log(f"    Would have {len(sc_all_comm_genes)} genes vs TCGA's {len(tcga_all_symbols)}")
    log(f"    Overlap: {len(total_overlap)} genes ({100*len(total_overlap)/len(sc_all_comm_genes):.1f}% of SC)")

    # Test at different merging levels
    # What if resolution 0.2 gives ~5 communities? Check the sweep files
    sweep_path = os.path.join(BASE_DIR, "data/FIG_4/04_communities/sweep")
    if os.path.exists(sweep_path):
        for res_file in sorted(os.listdir(sweep_path)):
            if not res_file.endswith(".csv"): continue
            res = res_file.replace("SC_partition_res","").replace(".csv","")
            rdf = pd.read_csv(os.path.join(sweep_path, res_file))
            n_comms = rdf["community"].nunique()
            sizes = rdf["community"].value_counts().sort_values(ascending=False)
            biggest = sizes.iloc[0] if len(sizes) > 0 else 0
            log(f"\n  Resolution {res}: {n_comms} communities, largest = {biggest} genes")

            # For lower resolutions, check if the biggest community has better TCGA overlap
            if float(res) <= 0.4:
                biggest_comm = sizes.index[0]
                biggest_genes = set(rdf[rdf["community"]==biggest_comm]["gene"])
                big_overlap = biggest_genes & tcga_all_symbols
                log(f"    Largest community overlap with TCGA: {len(big_overlap)} genes")
                if big_overlap:
                    # Which TCGA community?
                    for tc in sorted(tcga_comms.keys()):
                        tc_ov = tcga_comms[tc] & big_overlap
                        if tc_ov:
                            log(f"      TCGA C{tc}: {len(tc_ov)} genes — {', '.join(sorted(tc_ov)[:10])}")

    banner(f"COMPLETE | Elapsed: {datetime.now()-start}")

if __name__ == "__main__":
    main()
