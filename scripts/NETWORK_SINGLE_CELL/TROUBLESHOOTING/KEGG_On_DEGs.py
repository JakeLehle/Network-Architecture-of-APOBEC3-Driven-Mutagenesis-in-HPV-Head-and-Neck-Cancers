#!/usr/bin/env python3
"""
KEGG_On_DEGs.py
===============

Compute KEGG pathway enrichment on the single-cell DEGs for the two
tumor-vs-normal comparisons, using the shared Enrichr wrapper. This is the
compute half of the DEG->KEGG bridge (plots are a separate figure script).

Inputs (READ-ONLY):
  data/FIG_4/NETWORK_SBS2_VS_NORMAL/02_differential_expression/SC_diffexpr_stats.csv
  data/FIG_4/NETWORK_CNV_VS_NORMAL/02_differential_expression/SC_diffexpr_stats.csv
    columns: gene, p_value, fdr, log2FC, score, pct_test, pct_ref, selected

Decisions baked in (toggle at top to change):
  - Gene set: FDR < FDR_THRESHOLD, NO extra fold-change gate (matches the
    network-entry set), SPLIT by direction (up = log2FC>0, down = log2FC<0).
  - Background: the tested-gene set (all symbols in the DEG table) per
    comparison, so enrichment isn't saturated by huge lists vs a genome-wide
    background. Logs gseapy version + background mode so we can confirm on the
    first run that background actually took effect.

Outputs (per comparison, next to the DEG table):
  .../02_differential_expression/KEGG/KEGG_<COMP>_up.tsv
  .../02_differential_expression/KEGG/KEGG_<COMP>_down.tsv
  .../02_differential_expression/KEGG/KEGG_<COMP>_A3_cofactor_pathways.tsv
  data/FIG_4/02_differential_expression/KEGG_on_DEGs_report.txt   (top-level)

Placement: scripts/NETWORK_SINGLE_CELL/TROUBLESHOOTING/  (excluded from walkthrough)
Env:       conda run -n NETWORK python KEGG_On_DEGs.py
"""

import os
import sys
from datetime import datetime

import pandas as pd

# ---------------------------------------------------------------------------
# Config (single source of truth) + the shared wrapper, both one level up /
# alongside this script in TROUBLESHOOTING/.
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))
sys.path.insert(0, _HERE)  # so enrichr_kegg_wrapper.py (same dir) is importable

from network_config_SC import (
    BASE_DIR, FIG4_ROOT, DIR_02_DE,
    A3_SYMBOL_TO_ALIAS, FDR_THRESHOLD,
    ensure_dir,
)
from enrichr_kegg_wrapper import run_enrichr_kegg, KEGG_LIBRARY

# =============================================================================
# TOGGLES (the two locked decisions live here)
# =============================================================================
USE_TESTED_BACKGROUND = True    # decision 2: tested-gene set as Enrichr background
FC_GATE               = 0.0     # decision 1: extra |log2FC| gate on top of FDR (0 = none)
RUN_COMBINED          = False   # also run KEGG on all sig DEGs (up+down) per comparison

# The two canonical tumor-vs-normal comparisons (paper). SBS2_VS_CNV is dropped.
COMPARISONS = ["SBS2_VS_NORMAL", "CNV_VS_NORMAL"]

# Narrative genes to flag inside enriched pathways (for the volcano marking later).
A3_COFACTOR_TRACK = ["APOBEC3A", "APOBEC3B", "RALY", "HNRNPA2B1"]

# =============================================================================
# LOGGING
# =============================================================================
_report = []

def log(msg=""):
    print(msg, flush=True)
    _report.append(str(msg))

def banner(title, char="=", width=88):
    log("")
    log(char * width)
    log(f"  {title}")
    log(char * width)


# =============================================================================
# PATHS
# =============================================================================
def deg_path(comp):
    return os.path.join(FIG4_ROOT, f"NETWORK_{comp}",
                        "02_differential_expression", "SC_diffexpr_stats.csv")

def kegg_dir(comp):
    return ensure_dir(os.path.join(FIG4_ROOT, f"NETWORK_{comp}",
                                   "02_differential_expression", "KEGG"))


# =============================================================================
# HELPERS
# =============================================================================
def clean_symbols(series):
    """Symbols only; drop NaN and ENSG-looking (unannotated) IDs."""
    s = series.dropna().astype(str)
    return [g for g in s if not g.startswith("ENSG")]


def annotate_a3_pathways(kegg_df, direction):
    """Flag enriched pathways whose gene set contains an A3/cofactor gene.

    Enrichr's 'Genes' column is a ';'-separated list of hit symbols.
    """
    if kegg_df is None or len(kegg_df) == 0:
        return pd.DataFrame()
    track = set(A3_COFACTOR_TRACK)
    rows = []
    for _, r in kegg_df.iterrows():
        hits = {g.strip().upper() for g in str(r.get("Genes", "")).split(";") if g.strip()}
        matched = sorted(hits & track)
        if matched:
            rows.append({
                "direction": direction,
                "Term": r.get("Term", ""),
                "Adjusted P-value": r.get("Adjusted P-value", float("nan")),
                "Overlap": r.get("Overlap", ""),
                "A3_cofactor_hits": ", ".join(matched),
            })
    return pd.DataFrame(rows)


def show_top(kegg_df, n=8):
    if kegg_df is None or len(kegg_df) == 0:
        log("      (no results)")
        return
    sig = kegg_df[kegg_df["Adjusted P-value"] < 0.05]
    log(f"      {len(sig)} sig (FDR<0.05) of {len(kegg_df)} tested")
    for _, r in kegg_df.head(n).iterrows():
        star = "*" if r["Adjusted P-value"] < 0.05 else " "
        log(f"      {star} {str(r['Term'])[:52]:52s}  "
            f"adjP={r['Adjusted P-value']:.2e}  {r.get('Overlap','')}")


# =============================================================================
# PER-COMPARISON
# =============================================================================
def process_comparison(comp):
    banner(f"COMPARISON: {comp}")
    path = deg_path(comp)
    if not os.path.exists(path):
        log(f"  [SKIP] DEG table not found: {path}")
        return None

    df = pd.read_csv(path)
    log(f"  DEG table: {os.path.relpath(path, BASE_DIR)}")
    log(f"    {df.shape[0]} genes tested, columns: {list(df.columns)}")

    # Significance mask (explicit from fdr; cross-check the 'selected' column).
    mask = df["fdr"] < FDR_THRESHOLD
    if FC_GATE > 0:
        mask &= df["log2FC"].abs() >= FC_GATE
    if "selected" in df.columns and FC_GATE == 0:
        sel = df["selected"].astype(bool)
        if int(sel.sum()) != int(mask.sum()):
            log(f"    [WARN] 'selected' ({int(sel.sum())}) != fdr<{FDR_THRESHOLD} "
                f"({int(mask.sum())}); using explicit fdr mask.")

    up = df[mask & (df["log2FC"] > 0)]
    down = df[mask & (df["log2FC"] < 0)]
    log(f"    FDR<{FDR_THRESHOLD}"
        + (f" & |log2FC|>={FC_GATE}" if FC_GATE > 0 else "")
        + f": {int(mask.sum())}  (up={len(up)}, down={len(down)})")

    background = clean_symbols(df["gene"]) if USE_TESTED_BACKGROUND else None
    if background is not None:
        log(f"    Background: tested-gene set ({len(background)} symbols)")
    else:
        log(f"    Background: Enrichr default (genome-wide)")

    outdir = kegg_dir(comp)
    a3_frames = []

    to_run = [("up", up), ("down", down)]
    if RUN_COMBINED:
        to_run.append(("all", df[mask]))

    for direction, sub in to_run:
        genes = clean_symbols(sub["gene"])
        log(f"\n    --- KEGG: {comp} {direction} ({len(genes)} genes) ---")
        kegg = run_enrichr_kegg(genes, description=f"{comp}_{direction}",
                                background=background, logfn=log)
        show_top(kegg)
        if kegg is not None and len(kegg) > 0:
            out = os.path.join(outdir, f"KEGG_{comp}_{direction}.tsv")
            kegg.to_csv(out, sep="\t", index=False)
            log(f"      [SAVE] -> {os.path.relpath(out, BASE_DIR)}")
            a3_frames.append(annotate_a3_pathways(kegg, direction))

    # A3/cofactor-containing pathways across directions (volcano-marking seed).
    if a3_frames:
        a3_df = pd.concat([f for f in a3_frames if len(f) > 0], ignore_index=True) \
            if any(len(f) > 0 for f in a3_frames) else pd.DataFrame()
        if len(a3_df) > 0:
            a3_out = os.path.join(outdir, f"KEGG_{comp}_A3_cofactor_pathways.tsv")
            a3_df.sort_values("Adjusted P-value").to_csv(a3_out, sep="\t", index=False)
            log(f"\n    A3/cofactor-containing enriched pathways: {len(a3_df)}")
            for _, r in a3_df.sort_values("Adjusted P-value").head(10).iterrows():
                log(f"      [{r['direction']:4s}] {str(r['Term'])[:48]:48s}  "
                    f"adjP={r['Adjusted P-value']:.2e}  ({r['A3_cofactor_hits']})")
            log(f"      [SAVE] -> {os.path.relpath(a3_out, BASE_DIR)}")
        else:
            log(f"\n    No enriched pathway contains "
                f"{', '.join(A3_COFACTOR_TRACK)} (expected: A3A/A3B rarely in "
                f"canonical KEGG sets; cofactors may appear in RNA-processing terms).")

    return {"comp": comp, "n_sig": int(mask.sum()),
            "n_up": len(up), "n_down": len(down)}


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = datetime.now()
    banner("KEGG ON DEGs  (two tumor-vs-normal comparisons)")
    try:
        import gseapy as gp
        log(f"  gseapy version : {gp.__version__}")
    except Exception as e:
        log(f"  [FATAL] gseapy import failed: {e}")
        sys.exit(1)
    log(f"  KEGG library   : {KEGG_LIBRARY}")
    log(f"  FDR threshold  : {FDR_THRESHOLD}")
    log(f"  FC gate        : {FC_GATE} (0 = none)")
    log(f"  Background      : {'tested-gene set' if USE_TESTED_BACKGROUND else 'Enrichr default'}")
    log(f"  Comparisons    : {', '.join(COMPARISONS)}")

    summary = []
    for comp in COMPARISONS:
        try:
            r = process_comparison(comp)
            if r:
                summary.append(r)
        except Exception as e:
            log(f"\n  [ERROR] {comp}: {e}")
            import traceback
            log(traceback.format_exc())

    banner("SUMMARY")
    for r in summary:
        log(f"  {r['comp']:16s}  sig={r['n_sig']:5d}  "
            f"up={r['n_up']:5d}  down={r['n_down']:5d}")

    report_path = os.path.join(ensure_dir(DIR_02_DE), "KEGG_on_DEGs_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(_report))
    banner(f"COMPLETE | Elapsed: {datetime.now() - t0}")
    log(f"  Report -> {os.path.relpath(report_path, BASE_DIR)}")


if __name__ == "__main__":
    main()
