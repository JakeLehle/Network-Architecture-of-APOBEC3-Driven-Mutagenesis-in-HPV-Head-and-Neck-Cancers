#!/usr/bin/env python3
"""
Diagnostic_DE_Threshold_Comparison.py
======================================

Compare gene counts at raw p < 0.05 vs BH-adjusted p < 0.05 across all four
network comparisons (1 TCGA + 3 single-cell).

IMPORTANT: The bh_fdr() function used in the pipeline had a bug where the
monotonicity enforcement step scrambled gene order. This script uses the
corrected implementation and recomputes FDR from raw p-values.

Usage:
    python Diagnostic_DE_Threshold_Comparison.py
"""

import os
import pandas as pd
import numpy as np

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG4 = os.path.join(BASE_DIR, "data/FIG_4")

DE_FILES = {
    "TCGA-HNSC": {
        "paths": [
            os.path.join(BASE_DIR, "data/FIG_2/03_differential_expression/TCGA-HNSC/TCGA-HNSC_diffexpr_stats.csv"),
        ],
        "gene_col": "gene",
        "symbol_col": "gene_symbol",
        "pval_col": "p_value",
        "padj_col": "p_adj",
        "is_tcga": True,
    },
    "SC: SBS2_VS_CNV": {
        "paths": [
            os.path.join(FIG4, "NETWORK_SBS2_VS_CNV/02_differential_expression/SC_diffexpr_stats.csv"),
            os.path.join(FIG4, "02_differential_expression/SC_diffexpr_stats.csv"),
        ],
        "gene_col": "gene",
        "symbol_col": None,
        "pval_col": "p_value",
        "padj_col": "fdr",
        "is_tcga": False,
    },
    "SC: SBS2_VS_NORMAL": {
        "paths": [
            os.path.join(FIG4, "NETWORK_SBS2_VS_NORMAL/02_differential_expression/SC_diffexpr_stats.csv"),
        ],
        "gene_col": "gene",
        "symbol_col": None,
        "pval_col": "p_value",
        "padj_col": "fdr",
        "is_tcga": False,
    },
    "SC: CNV_VS_NORMAL": {
        "paths": [
            os.path.join(FIG4, "NETWORK_CNV_VS_NORMAL/02_differential_expression/SC_diffexpr_stats.csv"),
        ],
        "gene_col": "gene",
        "symbol_col": None,
        "pval_col": "p_value",
        "padj_col": "fdr",
        "is_tcga": False,
    },
}

HARRIS_ALL_PATH = os.path.join(FIG4, "00_input/Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(FIG4, "00_input/Harris_A3_interactors_A3B_only.txt")

A3_ENSG = {
    "ENSG00000128383": "APOBEC3A", "ENSG00000179750": "APOBEC3B",
    "ENSG00000244509": "APOBEC3C", "ENSG00000243811": "APOBEC3D",
    "ENSG00000128394": "APOBEC3F", "ENSG00000239713": "APOBEC3G",
    "ENSG00000100298": "APOBEC3H",
}
A3_ALIASES = {
    "APOBEC3A": "A3A", "APOBEC3B": "A3B", "APOBEC3C": "A3C",
    "APOBEC3D": "A3D", "APOBEC3F": "A3F", "APOBEC3G": "A3G",
    "APOBEC3H": "A3H",
}

OUTPUT_DIR = os.path.join(BASE_DIR, "data/DIAGNOSTICS")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# CORRECTED BH-FDR
# =============================================================================
# The pipeline's bh_fdr() had a bug: the monotonicity step reversed the array
# instead of mapping back to original gene order. This version sorts, adjusts,
# enforces monotonicity on the SORTED array, then maps back correctly.
# =============================================================================

def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction (CORRECTED)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    # Sort p-values ascending
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    # BH formula: p * n / rank
    ranks = np.arange(1, n + 1)
    adjusted = sorted_p * n / ranks
    # Monotonicity: cumulative min from right to left (high rank to low)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    # Clip to [0, 1]
    adjusted = np.clip(adjusted, 0, 1)
    # Map back to original gene order
    result = np.empty(n)
    result[order] = adjusted
    return result


# =============================================================================
# HELPERS
# =============================================================================

def find_file(candidates):
    for p in candidates:
        if os.path.exists(p):
            return p
    return None


def load_harris():
    harris_all, harris_a3b = set(), set()
    for path, s in [(HARRIS_ALL_PATH, harris_all), (HARRIS_A3B_PATH, harris_a3b)]:
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t")
            col = next((c for c in ["gene_symbol", "gene"] if c in df.columns), df.columns[0])
            s.update(df[col].dropna().astype(str).tolist())
    return harris_all, harris_a3b


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 100)
    print("  DE THRESHOLD DIAGNOSTIC: raw p < 0.05 vs BH-adjusted p < 0.05")
    print("  (using CORRECTED BH-FDR recomputed from raw p-values)")
    print("=" * 100)

    harris_all, harris_a3b = load_harris()
    print(f"\nHarris interactors: {len(harris_all)} all, {len(harris_a3b)} A3B-specific")

    ensg_to_symbol = {}
    emap = os.path.join(BASE_DIR, "data/FIG_2/01_cleaned_expression/ensg_to_symbol.json")
    if os.path.exists(emap):
        import json
        with open(emap) as f:
            ensg_to_symbol = json.load(f)
    symbol_to_ensg = {v: k for k, v in ensg_to_symbol.items()}

    summary_rows = []

    for name, cfg in DE_FILES.items():
        print(f"\n{'=' * 100}")
        print(f"  {name}")
        print(f"{'=' * 100}")

        path = find_file(cfg["paths"])
        if path is None:
            print(f"  FILE NOT FOUND. Checked:")
            for p in cfg["paths"]:
                exists = "EXISTS" if os.path.exists(p) else "MISSING"
                print(f"    [{exists}] {p}")
            continue

        print(f"  File: {path}")
        df = pd.read_csv(path)
        print(f"  Genes tested: {len(df)}")

        gc = cfg["gene_col"]
        sc = cfg["symbol_col"]
        pc = cfg["pval_col"]
        ac = cfg["padj_col"]
        is_tcga = cfg["is_tcga"]

        # Verify columns
        for label, col in [("gene", gc), ("pval", pc), ("padj", ac)]:
            if col not in df.columns:
                print(f"  ERROR: column '{col}' not in {list(df.columns)}")

        # Check stored FDR for sanity
        raw_vals = df[pc].values.astype(float)
        stored_vals = df[ac].values.astype(float)
        n_bad = int(np.sum(stored_vals < raw_vals * 0.99))
        if n_bad > 0:
            print(f"  Stored FDR has {n_bad}/{len(df)} violations (adj_p < raw_p)")
            print(f"  This confirms the bh_fdr() bug in the pipeline.")

        # Always recompute with corrected function
        print(f"  Recomputing BH-FDR with corrected function...")
        df["_adj_p"] = bh_fdr(raw_vals)

        # Verify correction
        recomp = df["_adj_p"].values
        n_still_bad = int(np.sum(recomp < raw_vals * 0.99))
        if n_still_bad > 0:
            print(f"  ERROR: Recomputed FDR still has {n_still_bad} violations!")
        else:
            print(f"  Recomputed FDR: 0 violations (all adj_p >= raw_p)")

        # Symbol column
        if sc and sc in df.columns:
            df["_sym"] = df[sc].astype(str)
        elif is_tcga:
            df["_sym"] = df[gc].map(lambda x: ensg_to_symbol.get(str(x), str(x)))
        else:
            df["_sym"] = df[gc].astype(str)

        # A3 family
        print(f"\n  A3 family in DE results:")
        print(f"  {'Gene':>6s}  {'raw_p':>12s}  {'adj_p':>12s}  {'raw<0.05':>8s}  {'adj<0.05':>8s}")
        print(f"  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*8}")

        a3_info = {}
        for alias in ["A3A", "A3B", "A3C", "A3D", "A3F", "A3G", "A3H"]:
            full = f"APOBEC3{alias[2]}"
            ensg = symbol_to_ensg.get(full, "")
            match = df[(df[gc] == full) | (df[gc] == ensg)]
            if sc and sc in df.columns and len(match) == 0:
                match = df[df[sc] == full]

            if len(match) > 0:
                rp = float(match.iloc[0][pc])
                ap = float(match.iloc[0]["_adj_p"])
                rpass = "YES" if rp < 0.05 else "no"
                apass = "YES" if ap < 0.05 else "no"
                print(f"  {alias:>6s}  {rp:>12.2e}  {ap:>12.2e}  {rpass:>8s}  {apass:>8s}")
                a3_info[alias] = {"raw_p": rp, "adj_p": ap}
            else:
                print(f"  {alias:>6s}  {'not in DE':>12s}  {'not in DE':>12s}  {'--':>8s}  {'--':>8s}")

        # Threshold comparison
        print(f"\n  {'Threshold':<16s}  {'Total':>7s}  {'A3':>4s}  "
              f"{'Harris':>7s}  {'H-A3B':>6s}  {'A3 members passing'}")
        print(f"  {'-'*16}  {'-'*7}  {'-'*4}  {'-'*7}  {'-'*6}  {'-'*35}")

        for label, col in [("raw_p < 0.05", pc), ("adj_p < 0.05", "_adj_p")]:
            mask = df[col] < 0.05
            n_pass = int(mask.sum())

            a3_pass = [a for a, v in a3_info.items()
                       if (v["raw_p"] if col == pc else v["adj_p"]) < 0.05]

            syms = set(df.loc[mask, "_sym"].tolist())
            if is_tcga:
                ids = set(df.loc[mask, gc].astype(str).tolist())
                ha = sum(1 for g in harris_all if g in syms or symbol_to_ensg.get(g, "") in ids)
                hb = sum(1 for g in harris_a3b if g in syms or symbol_to_ensg.get(g, "") in ids)
            else:
                ha = len(harris_all & syms)
                hb = len(harris_a3b & syms)

            a3s = ", ".join(sorted(a3_pass)) if a3_pass else "none"
            print(f"  {label:<16s}  {n_pass:>7d}  {len(a3_pass):>4d}  "
                  f"{ha:>7d}  {hb:>6d}  {a3s}")

            summary_rows.append({
                "comparison": name, "threshold": label, "total_genes": n_pass,
                "a3_count": len(a3_pass), "a3_members": a3s,
                "harris_all": ha, "harris_a3b": hb,
            })

        # Harris detail at adj_p < 0.05
        adj_mask = df["_adj_p"] < 0.05
        adj_syms = set(df.loc[adj_mask, "_sym"].tolist())
        if is_tcga:
            adj_ids = set(df.loc[adj_mask, gc].astype(str).tolist())
            hin = [g for g in sorted(harris_all)
                   if g in adj_syms or symbol_to_ensg.get(g, "") in adj_ids]
        else:
            hin = sorted(harris_all & adj_syms)

        if hin:
            print(f"\n  Harris interactors at adj_p < 0.05 ({len(hin)}):")
            for i in range(0, len(hin), 6):
                print(f"    {', '.join(hin[i:i+6])}")
        else:
            print(f"\n  Harris interactors at adj_p < 0.05: NONE")

    # Summary
    print(f"\n\n{'=' * 100}")
    print(f"  SUMMARY TABLE")
    print(f"{'=' * 100}")

    sdf = pd.DataFrame(summary_rows)
    if len(sdf) > 0:
        for thr in ["raw_p < 0.05", "adj_p < 0.05"]:
            sub = sdf[sdf["threshold"] == thr]
            if len(sub) > 0:
                print(f"\n  {thr}:")
                print(f"  {'Comparison':<25s}  {'Genes':>7s}  {'A3':>4s}  "
                      f"{'Harris':>7s}  {'A3 members'}")
                print(f"  {'-'*25}  {'-'*7}  {'-'*4}  {'-'*7}  {'-'*30}")
                for _, r in sub.iterrows():
                    print(f"  {r['comparison']:<25s}  {r['total_genes']:>7d}  "
                          f"{r['a3_count']:>4d}  {r['harris_all']:>7d}  {r['a3_members']}")

        out = os.path.join(OUTPUT_DIR, "DE_threshold_comparison.csv")
        sdf.to_csv(out, index=False)
        print(f"\n  [SAVE] {out}")

    # Note about pipeline fix needed
    print(f"\n  NOTE: The bh_fdr() function in Step03 (TCGA) and Step01 (SC)")
    print(f"  has a gene-order scrambling bug. All stored FDR values in the")
    print(f"  pipeline DE output files are incorrect. The corrected function")
    print(f"  must be patched into both scripts and DE steps rerun.")

    print(f"\n{'=' * 100}")
    print(f"  DONE")
    print(f"{'=' * 100}")


if __name__ == "__main__":
    main()
