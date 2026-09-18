#!/usr/bin/env python3
"""
Rescore_LOPO_Against_Current_Partition.py
=========================================

Re-scores the May 2026 leave-one-patient-out (LOPO) runs against the CURRENT
full-analysis partition and the CURRENT activating-chain gene set.

WHY THIS EXISTS
---------------
The three LOPO runs under data/FIG_5/03_sensitivity/ were produced on
2026-05-21. Two things have changed since:

  1. The full-analysis partition moved from 23 communities at Leiden
     resolution 0.70 to 25 communities at 0.80. This was NOT a data change:
     SC_diffexpr_stats.csv and SC_corr_DIFF.pkl are byte-identical between
     NETWORK_SBS2_VS_NORMAL_PRE_RESELECT/ and NETWORK_SBS2_VS_NORMAL/
     (verified by md5sum, 2026-09-15). The cause was the resolution sweep
     grid in network_config_SC.py changing from 0.1 steps (0.1 ... 0.8) to
     COMMUNITY_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0], which removed 0.70
     from the candidate set. 0.80 is the best score on the current grid.

     Consequence: community_ari_with_full in each LOPO summary was computed
     against the 23-community partition, which is no longer the reference.
     It must be recomputed. Jaccard and gene overlap are unaffected because
     the gene set did not change.

  2. ACTIVATING_CHAIN_GENES gained CHMP4B (9 genes -> 10). The recovery
     criterion is membership in the reconstructed network, confirmed by
     grepping LOPO_*_communities.tsv against the stored gene lists: the
     three stored lists match presence in the network exactly.

WHAT IT DOES NOT DO
-------------------
Does not rebuild any network. The LOPO reconstructions remain valid because
the underlying SBS2-HIGH and NORMAL cell sets are unchanged (n_high in each
summary equals 546 minus that patient's contribution, and normal_removed is
0 in all three).

SIDE EFFECT
-----------
harris_in_network / harris_total in the existing summaries read 0 / 175.
Both are wrong, and it is the known pre-June reader bug documented in the
patient_config.py docstring: the old loader read whole tab-delimited lines
including the header, so no gene symbol ever matched (hence 0) and the
header inflated the total to 175. No manuscript number depends on this
field (Methods 6.5 scores chain genes, ARI, Jaccard, and wall integrity),
but since the row is being rewritten anyway, this script recomputes it
properly. Set FIX_HARRIS = False to leave those two columns untouched.

Usage:
  conda run -n NETWORK python Rescore_LOPO_Against_Current_Partition.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import shutil
import ast
import pandas as pd
from datetime import datetime
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from patient_config import (
    COMMUNITIES_DIR, DIR_03_SENSITIVITY,
    ACTIVATING_CHAIN_GENES, INHIBITING_CHAIN_ANCHORS,
    HIGH_CONTRIBUTORS,
    load_harris_interactors,
    banner, log,
)

# =============================================================================
# SETTINGS
# =============================================================================
WRITE_IN_PLACE = True    # overwrite LOPO_*_summary.tsv (a .bak copy is kept)
FIX_HARRIS     = True    # also repair the broken harris_in_network / _total

FULL_PARTITION = os.path.join(COMMUNITIES_DIR, "SC_best_partition.csv")
FULL_PARAMS    = os.path.join(COMMUNITIES_DIR, "SC_selected_parameters.txt")

STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


# =============================================================================
# HELPERS
# =============================================================================

def load_full_partition():
    """Current full-analysis gene -> community map."""
    df = pd.read_csv(FULL_PARTITION)
    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    comm_col = "community" if "community" in df.columns else df.columns[1]
    mapping = dict(zip(df[gene_col].astype(str), df[comm_col]))
    log(f"  full partition: {len(mapping):,} genes, "
        f"{df[comm_col].nunique()} communities")
    return mapping


def load_full_params():
    params = {}
    if os.path.exists(FULL_PARAMS):
        with open(FULL_PARAMS) as f:
            for line in f:
                if "=" in line:
                    k, v = line.strip().split("=", 1)
                    params[k.strip()] = v.strip()
    return params


def load_lopo_partition(path):
    """LOPO gene -> community map from LOPO_*_communities.tsv."""
    df = pd.read_csv(path, sep="\t")
    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    comm_col = "community" if "community" in df.columns else df.columns[1]
    return dict(zip(df[gene_col].astype(str), df[comm_col]))


def score_run(patient, full_map, harris_all):
    """Recompute ARI, Jaccard, overlap, and chain recovery for one LOPO run."""
    run_dir = os.path.join(DIR_03_SENSITIVITY, f"LOPO_{patient}")
    comm_path = os.path.join(run_dir, f"LOPO_{patient}_communities.tsv")
    summ_path = os.path.join(run_dir, f"LOPO_{patient}_summary.tsv")

    for p in (comm_path, summ_path):
        if not os.path.exists(p):
            log(f"  [SKIP] {patient}: missing {os.path.basename(p)}")
            return None

    lopo_map = load_lopo_partition(comm_path)
    lopo_genes = set(lopo_map)
    full_genes = set(full_map)

    # ARI is only defined on genes present in BOTH partitions.
    shared = sorted(lopo_genes & full_genes)
    ari = adjusted_rand_score(
        [full_map[g] for g in shared],
        [lopo_map[g] for g in shared],
    )

    overlap = len(shared)
    union = len(lopo_genes | full_genes)
    jaccard = overlap / union if union else 0.0

    chain_present = [g for g in ACTIVATING_CHAIN_GENES if g in lopo_genes]
    chain_missing = [g for g in ACTIVATING_CHAIN_GENES if g not in lopo_genes]
    inhib_present = [g for g in INHIBITING_CHAIN_ANCHORS if g in lopo_genes]

    harris_n = len(harris_all & lopo_genes) if harris_all else None

    return {
        "patient": patient,
        "n_lopo_genes": len(lopo_genes),
        "n_shared": overlap,
        "ari": ari,
        "jaccard": jaccard,
        "chain_present": chain_present,
        "chain_missing": chain_missing,
        "chain_recovered": len(chain_present),
        "chain_total": len(ACTIVATING_CHAIN_GENES),
        "inhib_recovered": len(inhib_present),
        "inhib_total": len(INHIBITING_CHAIN_ANCHORS),
        "harris_n": harris_n,
        "harris_total": len(harris_all) if harris_all else None,
        "summ_path": summ_path,
    }


def rewrite_summary(res):
    """Update the four (or six) changed columns, preserving every other field."""
    path = res["summ_path"]
    df = pd.read_csv(path, sep="\t")

    before = {
        "chain": (df.at[0, "activating_chain_recovered"],
                  df.at[0, "activating_chain_total"]),
        "ari": df.at[0, "community_ari_with_full"],
        "jaccard": df.at[0, "jaccard_with_full"],
        "overlap": df.at[0, "gene_overlap_with_full"],
    }

    backup = f"{path}.bak_{STAMP}"
    shutil.copy2(path, backup)

    df.at[0, "activating_chain_recovered"] = res["chain_recovered"]
    df.at[0, "activating_chain_total"]     = res["chain_total"]
    df.at[0, "activating_chain_genes"]     = str(res["chain_present"])
    df.at[0, "inhibiting_anchors_recovered"] = res["inhib_recovered"]
    df.at[0, "community_ari_with_full"]    = res["ari"]
    df.at[0, "jaccard_with_full"]          = res["jaccard"]
    df.at[0, "gene_overlap_with_full"]     = res["n_shared"]

    if FIX_HARRIS and res["harris_n"] is not None:
        df.at[0, "harris_in_network"] = res["harris_n"]
        df.at[0, "harris_total"]      = res["harris_total"]

    df.to_csv(path, sep="\t", index=False)
    return before, backup


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("LOPO RE-SCORE AGAINST CURRENT PARTITION")

    params = load_full_params()
    log(f"  reference network: {params.get('N_GENES')} genes, "
        f"{params.get('N_COMMUNITIES')} communities, "
        f"threshold {params.get('DIFF_THRESHOLD')}, "
        f"resolution {params.get('LEIDEN_RESOLUTION')}")
    log(f"  activating chain set: {len(ACTIVATING_CHAIN_GENES)} genes")
    log(f"    {', '.join(ACTIVATING_CHAIN_GENES)}")

    full_map = load_full_partition()

    harris_all = set()
    if FIX_HARRIS:
        try:
            harris_all, _ = load_harris_interactors()
        except Exception as e:
            log(f"  [WARN] Harris load failed ({e}); leaving those columns alone")
            harris_all = set()

    patients = [p.replace("Patient ", "") for p in HIGH_CONTRIBUTORS]

    banner("PER-RUN RESULTS")
    results = []
    for patient in patients:
        res = score_run(patient, full_map, harris_all)
        if res is None:
            continue
        results.append(res)
        log("")
        log(f"  {patient}")
        log(f"    genes in LOPO network : {res['n_lopo_genes']:,}")
        log(f"    shared with full      : {res['n_shared']:,}")
        log(f"    ARI  (recomputed)     : {res['ari']:.4f}")
        log(f"    Jaccard               : {res['jaccard']:.4f}")
        log(f"    chain recovered       : {res['chain_recovered']}/{res['chain_total']}")
        log(f"      present : {', '.join(res['chain_present']) or 'none'}")
        log(f"      missing : {', '.join(res['chain_missing']) or 'none'}")
        log(f"    inhibiting anchors    : {res['inhib_recovered']}/{res['inhib_total']}")
        if res["harris_n"] is not None:
            log(f"    Harris interactors    : {res['harris_n']}/{res['harris_total']}")

    if not results:
        log("  nothing scored; check DIR_03_SENSITIVITY")
        return

    banner("SUMMARY TABLE (paste-ready)")
    log("")
    log(f"  {'Patient':<10} {'Chain':>8} {'Inhib':>7} {'ARI':>8} "
        f"{'Jaccard':>9} {'Overlap':>9} {'Nodes':>8}")
    log(f"  {'-'*10} {'-'*8} {'-'*7} {'-'*8} {'-'*9} {'-'*9} {'-'*8}")
    for r in results:
        log(f"  {r['patient']:<10} "
            f"{r['chain_recovered']}/{r['chain_total']:<6} "
            f"{r['inhib_recovered']}/{r['inhib_total']:<5} "
            f"{r['ari']:>8.3f} {r['jaccard']:>9.3f} "
            f"{r['n_shared']:>9,} {r['n_lopo_genes']:>8,}")

    # Genes lost in every run are the least robust members of the set.
    always_missing = set(ACTIVATING_CHAIN_GENES)
    for r in results:
        always_missing &= set(r["chain_missing"])
    if always_missing:
        log("")
        log(f"  absent from ALL reconstructions: {', '.join(sorted(always_missing))}")

    if WRITE_IN_PLACE:
        banner("REWRITING SUMMARY FILES")
        for r in results:
            before, backup = rewrite_summary(r)
            log("")
            log(f"  {r['patient']}")
            log(f"    chain   {before['chain'][0]}/{before['chain'][1]} "
                f"-> {r['chain_recovered']}/{r['chain_total']}")
            log(f"    ARI     {float(before['ari']):.4f} -> {r['ari']:.4f}")
            log(f"    Jaccard {float(before['jaccard']):.4f} -> {r['jaccard']:.4f}")
            log(f"    overlap {before['overlap']} -> {r['n_shared']}")
            log(f"    backup  {os.path.basename(backup)}")
        log("")
        log("  Re-run Generate_Supplemental_Patient_Effects.py to refresh Panel D.")
    else:
        banner("DRY RUN")
        log("  WRITE_IN_PLACE = False; no files modified.")

    banner("COMPLETE")


if __name__ == "__main__":
    main()
