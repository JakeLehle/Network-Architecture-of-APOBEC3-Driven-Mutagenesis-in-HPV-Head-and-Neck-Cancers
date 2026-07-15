#!/usr/bin/env python3
"""
Diagnostic_DEG_KEGG_State.py
============================

READ-ONLY diagnostic. Confirms the on-disk state of the single-cell DEG stage
before any KEGG / volcano figure code is written. Modifies nothing.

Answers, from disk (not from assumptions):
  1. Where does each comparison's SC_diffexpr_stats.csv actually live?
     (Step01 writes fixed filenames to 02_differential_expression/, but the
      three-network runner moves outputs into data/FIG_4/NETWORK_*/ ; this
      script discovers the real layout by walking FIG_4.)
  2. What are the real column names + row counts in each stats table, and how
     many genes pass FDR < FDR_THRESHOLD?
  3. Are APOBEC3A, APOBEC3B, RALY, HNRNPA2B1 present, and what are their
     log2FC / p / FDR / selected status in each comparison?
  4. Does ANY KEGG / Enrichr output already exist anywhere under FIG_4?
  5. Flags the retired pre-reselection KEGG (Phase5A_v2 -> FIG_6/.../DEG) so it
     is not reused by accident.

Every file reported prints its mtime, because exact-integer matches and stale
artifacts after reselection are a known failure mode here.

Placement: scripts/NETWORK_SINGLE_CELL/TROUBLESHOOTING/  (excluded from walkthrough)
Run:       conda run -n NETWORK python Diagnostic_DEG_KEGG_State.py
"""

import os
import sys
import glob
from datetime import datetime

import pandas as pd

# ---------------------------------------------------------------------------
# Import the single-cell config (single source of truth for paths/thresholds).
# This diagnostic is meant to live in TROUBLESHOOTING/, one level below the
# config, so add the parent dir to the path before importing.
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

try:
    from network_config_SC import (
        BASE_DIR, FIG4_ROOT, DIR_02_DE,
        A3_GENES_SYMBOLS, A3_SYMBOL_TO_ALIAS,
        FDR_THRESHOLD,
    )
except Exception as e:  # pragma: no cover
    print(f"[FATAL] Could not import network_config_SC: {e}")
    print("        Run this from scripts/NETWORK_SINGLE_CELL/TROUBLESHOOTING/,")
    print("        or adjust sys.path to point at network_config_SC.py.")
    sys.exit(1)

# Cofactors we want to track alongside the A3 enzymes on the volcano.
COFACTORS = ["RALY", "HNRNPA2B1"]
TRACK_GENES = list(A3_GENES_SYMBOLS) + COFACTORS

# Retired pre-reselection KEGG location (do NOT reuse) -- flagged, not read.
RETIRED_KEGG_DIR = os.path.join(BASE_DIR, "data", "FIG_6",
                                "04_population_profiles_v2", "DEG")


def banner(title, char="=", width=88):
    print("\n" + char * width)
    print(f"  {title}")
    print(char * width)


def mtime(path):
    try:
        return datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")
    except OSError:
        return "??"


def infer_comparison(path):
    """Guess which comparison a stats file belongs to from its path."""
    p = path.upper()
    if "SBS2_VS_NORMAL" in p:
        return "SBS2_VS_NORMAL"
    if "CNV_VS_NORMAL" in p:
        return "CNV_VS_NORMAL"
    if "SBS2_VS_CNV" in p:
        return "SBS2_VS_CNV (dropped from paper)"
    return "UNKNOWN (shared/last-run dir?)"


def pick_col(cols, candidates):
    """Return the first candidate present in cols (case-insensitive), else None."""
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def summarize_stats_file(path):
    banner(f"DEG TABLE: {os.path.relpath(path, BASE_DIR)}", char="-")
    print(f"  comparison (inferred) : {infer_comparison(path)}")
    print(f"  mtime                 : {mtime(path)}")

    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"  [ERROR] could not read: {e}")
        return

    print(f"  shape                 : {df.shape[0]} rows x {df.shape[1]} cols")
    print(f"  columns               : {list(df.columns)}")

    gene_col = pick_col(df.columns, ["gene", "names", "gene_symbol", "symbol"])
    fdr_col = pick_col(df.columns, ["fdr", "pvals_adj", "padj", "adj_p", "qvalue"])
    p_col = pick_col(df.columns, ["p_value", "pvals", "pvalue", "pval"])
    fc_col = pick_col(df.columns, ["log2FC", "logfoldchanges", "log2fc", "lfc"])

    print(f"  resolved gene/p/fdr/fc: {gene_col} / {p_col} / {fdr_col} / {fc_col}")

    if fdr_col is not None:
        try:
            n_sig = int((pd.to_numeric(df[fdr_col], errors="coerce") < FDR_THRESHOLD).sum())
            print(f"  genes FDR < {FDR_THRESHOLD}       : {n_sig}")
            if fc_col is not None:
                fc = pd.to_numeric(df[fc_col], errors="coerce")
                fdr = pd.to_numeric(df[fdr_col], errors="coerce")
                for thr in (0.5, 1.0):
                    m = (fdr < FDR_THRESHOLD) & (fc.abs() >= thr)
                    print(f"    of those, |log2FC| >= {thr}: {int(m.sum())}")
        except Exception as e:
            print(f"  [WARN] could not count significant genes: {e}")

    # A3 + cofactor rows
    if gene_col is not None:
        print(f"\n  {'gene':12s} {'log2FC':>9s} {'p':>12s} {'FDR':>12s}  {'FDR<thr':>7s}")
        print(f"  {'-'*12} {'-'*9} {'-'*12} {'-'*12}  {'-'*7}")
        gset = df.set_index(df[gene_col].astype(str))
        for g in TRACK_GENES:
            alias = A3_SYMBOL_TO_ALIAS.get(g, g)
            label = f"{g}" if g == alias else f"{g}({alias})"
            if g in gset.index:
                row = gset.loc[g]
                if isinstance(row, pd.DataFrame):  # duplicate symbol guard
                    row = row.iloc[0]
                fc = row[fc_col] if fc_col else float("nan")
                pv = row[p_col] if p_col else float("nan")
                fd = row[fdr_col] if fdr_col else float("nan")
                try:
                    sig = "yes" if float(fd) < FDR_THRESHOLD else "no"
                except (TypeError, ValueError):
                    sig = "?"
                print(f"  {label:12s} {float(fc):9.3f} {float(pv):12.3e} "
                      f"{float(fd):12.3e}  {sig:>7s}")
            else:
                print(f"  {label:12s} {'--- not present in table ---':>44s}")


def main():
    banner("DEG / KEGG ON-DISK STATE  (read-only)")
    print(f"  BASE_DIR      : {BASE_DIR}")
    print(f"  FIG4_ROOT     : {FIG4_ROOT}")
    print(f"  DIR_02_DE     : {DIR_02_DE}")
    print(f"  FDR_THRESHOLD : {FDR_THRESHOLD}")
    print(f"  tracking      : {', '.join(TRACK_GENES)}")

    if not os.path.isdir(FIG4_ROOT):
        print(f"\n[FATAL] FIG4_ROOT does not exist: {FIG4_ROOT}")
        sys.exit(1)

    # -- 1. Locate every SC_diffexpr_stats.csv under FIG_4 --------------------
    banner("1. DEG STATS TABLES FOUND UNDER FIG_4")
    stats_files = sorted(glob.glob(os.path.join(FIG4_ROOT, "**", "SC_diffexpr_stats.csv"),
                                   recursive=True))
    if not stats_files:
        print("  NONE FOUND. Has Step01 been run on the current groups yet?")
    else:
        for f in stats_files:
            print(f"  [{mtime(f)}]  {os.path.relpath(f, BASE_DIR)}")

    # Also list the selected-gene tables (these are the network-entry lists)
    sel_files = sorted(glob.glob(os.path.join(FIG4_ROOT, "**", "SC_selected_genes*.csv"),
                                 recursive=True))
    if sel_files:
        print("\n  Selected-gene tables (network-entry lists):")
        for f in sel_files:
            try:
                n = sum(1 for _ in open(f)) - 1
            except OSError:
                n = "?"
            print(f"  [{mtime(f)}]  {os.path.relpath(f, BASE_DIR)}  ({n} genes)")

    # -- 2. Summarize each stats table ---------------------------------------
    for f in stats_files:
        summarize_stats_file(f)

    # -- 3. Any existing KEGG / Enrichr output under FIG_4? -------------------
    banner("3. EXISTING KEGG / ENRICHR OUTPUT UNDER FIG_4")
    kegg_hits = []
    for pat in ("*kegg*", "*KEGG*", "*enrichr*", "*Enrichr*", "*gsea*", "*GSEA*"):
        kegg_hits += glob.glob(os.path.join(FIG4_ROOT, "**", pat), recursive=True)
    kegg_hits = sorted(set(p for p in kegg_hits if os.path.isfile(p)))
    if not kegg_hits:
        print("  NONE. -> KEGG has not been computed on the current groups. "
              "This confirms it needs to be added.")
    else:
        for f in kegg_hits:
            print(f"  [{mtime(f)}]  {os.path.relpath(f, BASE_DIR)}")

    # -- 4. Flag the retired pre-reselection KEGG (do not reuse) --------------
    banner("4. RETIRED PRE-RESELECTION KEGG (FLAG ONLY -- DO NOT REUSE)")
    if os.path.isdir(RETIRED_KEGG_DIR):
        retired = sorted(glob.glob(os.path.join(RETIRED_KEGG_DIR, "*GSEA*")) +
                         glob.glob(os.path.join(RETIRED_KEGG_DIR, "*kegg*")) +
                         glob.glob(os.path.join(RETIRED_KEGG_DIR, "*KEGG*")))
        if retired:
            print("  Found retired-group KEGG/GSEA (Stealth_CNV / Normal_Control):")
            for f in retired:
                print(f"  [{mtime(f)}]  {os.path.relpath(f, BASE_DIR)}")
            print("  ^ pre-reselection group names. NOT for the current figures.")
        else:
            print("  Retired DEG dir exists but no KEGG/GSEA files in it. OK.")
    else:
        print(f"  Not present: {os.path.relpath(RETIRED_KEGG_DIR, BASE_DIR)} (fine).")

    banner("DIAGNOSTIC COMPLETE -- nothing was modified")


if __name__ == "__main__":
    main()
