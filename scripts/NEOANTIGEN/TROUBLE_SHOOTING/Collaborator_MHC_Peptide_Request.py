#!/usr/bin/env python3
"""
Collaborator_MHC_Peptide_Request.py
===================================
One-off MHC-I binding prediction for an arbitrary list of peptides against an
arbitrary list of HLA class I alleles. Written to answer an external
collaborator request (G. Ippolito, Aug 2026) rather than as a pipeline stage.

This is deliberately NOT part of the Figure 7 neoantigen pipeline. It does no
proteome lookup, no variant parsing, and no peptide generation. It takes literal
peptide sequences in and writes a binding table out. It lives in
TROUBLE_SHOOTING/ and is therefore excluded from the Q5 walkthrough.

Reuses the predictor-loading and prediction pattern from:
    scripts/NEOANTIGEN/Step03_MHCflurry_Binding.py
    scripts/NEOANTIGEN/TROUBLE_SHOOTING/resolve_KRT6B_HLA.py

Both MHCflurry predictors are run:
    Class1AffinityPredictor      -> IC50 (nM), 5-95 percentile interval,
                                    affinity percentile rank
    Class1PresentationPredictor  -> antigen processing score, presentation
                                    score, presentation percentile

Every peptide x allele pair is reported. Nothing is filtered out. Binder flags
use the same thresholds as the manuscript pipeline (< 500 nM binder,
< 50 nM strong binder) but are reported as columns, not applied as filters.

Outputs (to data/FIG_7/COLLABORATOR_MHC_REQUEST/):
    {tag}_long.tsv              one row per peptide-allele pair, all metrics
    {tag}_ic50_matrix.tsv       wide matrix, peptides x alleles, IC50 nM
    {tag}_presentation_matrix.tsv  wide matrix, presentation score
    {tag}_table.txt             fixed-width, email-ready
    {tag}_heatmap.pdf / .png    two-panel heatmap, 300 DPI
    {tag}_run_report.txt        full console log

USAGE
-----
Full run (prediction + figure), NEOANTIGEN env:
    conda run -n NEOANTIGEN python Collaborator_MHC_Peptide_Request.py

Prediction only, no matplotlib needed:
    conda run -n NEOANTIGEN python Collaborator_MHC_Peptide_Request.py --table-only

Figure only, from a previously written long TSV (NETWORK env):
    conda run -n NETWORK python Collaborator_MHC_Peptide_Request.py \
        --figure-from-tsv <path to {tag}_long.tsv>

Override the peptide / allele lists without editing the file:
    ... --peptides SIINFEKL GILGFVFTL --alleles "A*02:01" "B*07:02"

Author: Jake Lehle / Claude
"""

import os
import sys
import argparse

import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/COLLABORATOR_MHC_REQUEST")

DEFAULT_TAG = "Ippolito_2026_08"

# Peptides as supplied by the collaborator, in the order supplied.
DEFAULT_PEPTIDES = [
    "NSSKVSQNY",
    "GSEELRSLY",
    "RDYVDRFFKTL",
    "YVDRFFKTL",
    "RYPLTFGWCF",
    "RYLKDQQLL",
]

# Alleles as supplied by the collaborator. Star notation is normalized below.
DEFAULT_ALLELES = [
    "A*01:01",
    "A*02:01",
    "A*24:02",
]

# Same thresholds the Figure 7 pipeline uses. Reported, not applied.
BIND_THRESH = 500.0
STRONG_THRESH = 50.0

# MHCflurry supported peptide length window for the affinity predictor.
MIN_PEP_LEN, MAX_PEP_LEN = 8, 15

# Probe peptide used to test whether an allele string is parseable.
PROBE_PEPTIDE = "GILGFVFTL"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# --- figure style -------------------------------------------------------------
DPI = 300
FS_TITLE = 34
FS_LABEL = 30
FS_TICK = 28
FS_ANNOT = 28

COLOR_STRONG = '#ed6a5a'   # coral,   low IC50 / high presentation
COLOR_MID = '#F6D155'      # mustard
COLOR_WEAK = '#EDEDED'     # light gray, high IC50 / low presentation
COLOR_OUTLINE = '#222222'

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []


def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))


def banner(title, char="="):
    log("")
    log(char * 78)
    log(f"  {title}")
    log(char * 78)


# =============================================================================
# INPUT VALIDATION
# =============================================================================
def validate_peptides(peptides):
    """Drop nothing silently. Report every rejection with a reason."""
    good, bad = [], []
    seen = set()
    for p in peptides:
        p = str(p).strip().upper()
        if not p:
            continue
        if p in seen:
            bad.append((p, "duplicate in input list"))
            continue
        seen.add(p)
        offenders = sorted(set(p) - VALID_AA)
        if offenders:
            bad.append((p, f"non-standard residue(s): {''.join(offenders)}"))
        elif not (MIN_PEP_LEN <= len(p) <= MAX_PEP_LEN):
            bad.append((p, f"length {len(p)} outside MHCflurry range "
                           f"{MIN_PEP_LEN}-{MAX_PEP_LEN}"))
        else:
            good.append(p)
    return good, bad


def normalize_allele(predictor, allele):
    """
    Return an allele string MHCflurry accepts, or None.

    Step03 writes alleles compressed ('HLA-A0201'); resolve_KRT6B_HLA.py writes
    them starred ('HLA-A*02:01'). MHCflurry takes either. Try the starred form
    first because that is how collaborators write them, then fall back.
    """
    raw = str(allele).strip()
    if not raw.upper().startswith("HLA-"):
        raw = "HLA-" + raw

    candidates = [raw, raw.replace("*", "").replace(":", "")]
    for cand in candidates:
        try:
            predictor.predict_to_dataframe(peptides=[PROBE_PEPTIDE],
                                           alleles=[cand])
            return cand
        except Exception:
            continue
    return None


def report_nested(peptides):
    """
    Flag peptides fully contained inside another peptide in the same list.
    Purely structural. Makes nested length variants obvious in the output
    without asserting anything about the origin of the sequences.
    """
    nested = []
    for short in peptides:
        for long_ in peptides:
            if short != long_ and short in long_:
                nested.append((short, long_))
    return nested


# =============================================================================
# PREDICTION
# =============================================================================
def run_affinity(peptides, alleles_raw):
    """Class1AffinityPredictor across every peptide x allele pair."""
    from mhcflurry import Class1AffinityPredictor

    log("  Loading Class1AffinityPredictor ...")
    predictor = Class1AffinityPredictor.load()
    log("  Loaded.")

    resolved, unsupported = {}, []
    for a in alleles_raw:
        norm = normalize_allele(predictor, a)
        if norm is None:
            unsupported.append(a)
            log(f"  WARNING: allele '{a}' not recognized by MHCflurry, skipping")
        else:
            resolved[a] = norm
            log(f"  allele '{a}' -> '{norm}'")

    if not resolved:
        log("  FATAL: no usable alleles.")
        sys.exit(1)

    rows = []
    for a_in, a_norm in resolved.items():
        preds = predictor.predict_to_dataframe(
            peptides=peptides,
            alleles=[a_norm] * len(peptides),
        )
        for i, pep in enumerate(peptides):
            r = preds.iloc[i]
            rows.append({
                "peptide": pep,
                "peptide_length": len(pep),
                "allele": a_in,
                "allele_mhcflurry": a_norm,
                "ic50_nM": float(r["prediction"]),
                "ic50_low_nM": float(r["prediction_low"])
                if "prediction_low" in r else np.nan,
                "ic50_high_nM": float(r["prediction_high"])
                if "prediction_high" in r else np.nan,
                "affinity_percentile": float(r["prediction_percentile"])
                if "prediction_percentile" in r
                and pd.notna(r["prediction_percentile"]) else np.nan,
            })
        log(f"  scored {len(peptides)} peptides on {a_norm}")

    return pd.DataFrame(rows), resolved, unsupported


def run_presentation(peptides, resolved):
    """
    Class1PresentationPredictor, one call per allele.

    predict() with a multi-allele genotype returns only the best allele per
    peptide. Calling it once per allele with a single-allele genotype forces a
    value for every peptide x allele pair. The dict form
    {name: [allele]} also stamps sample_name with the allele, which makes the
    merge unambiguous.

    Returns an empty frame if the presentation predictor is unavailable, so a
    failure here never costs us the affinity table.
    """
    try:
        from mhcflurry import Class1PresentationPredictor
    except Exception as e:
        log(f"  Class1PresentationPredictor unavailable ({e}); "
            f"skipping presentation metrics")
        return pd.DataFrame()

    try:
        log("  Loading Class1PresentationPredictor ...")
        pres = Class1PresentationPredictor.load()
        log("  Loaded.")
    except Exception as e:
        log(f"  Presentation predictor failed to load ({e}); skipping")
        return pd.DataFrame()

    rows = []
    for a_in, a_norm in resolved.items():
        try:
            df = pres.predict(
                peptides=peptides,
                alleles={a_norm: [a_norm]},
                include_affinity_percentile=True,
                verbose=0,
            )
        except Exception as e:
            log(f"  WARNING: presentation prediction failed for {a_norm}: "
                f"{str(e)[:120]}")
            continue

        for _, r in df.iterrows():
            rows.append({
                "peptide": str(r["peptide"]),
                "allele": a_in,
                "processing_score": float(r["processing_score"])
                if "processing_score" in r
                and pd.notna(r["processing_score"]) else np.nan,
                "presentation_score": float(r["presentation_score"])
                if "presentation_score" in r
                and pd.notna(r["presentation_score"]) else np.nan,
                "presentation_percentile": float(r["presentation_percentile"])
                if "presentation_percentile" in r
                and pd.notna(r["presentation_percentile"]) else np.nan,
            })
        log(f"  presentation scored {len(df)} peptides on {a_norm}")

    return pd.DataFrame(rows)


# =============================================================================
# TABLE ASSEMBLY
# =============================================================================
def assemble(aff, pres, peptides, allele_order):
    """Merge, flag, and order. Peptide and allele order follow the input."""
    df = aff.copy()
    if len(pres):
        df = df.merge(pres, on=["peptide", "allele"], how="left")
    else:
        for c in ("processing_score", "presentation_score",
                  "presentation_percentile"):
            df[c] = np.nan

    df["binder_lt500"] = df["ic50_nM"] < BIND_THRESH
    df["strong_binder_lt50"] = df["ic50_nM"] < STRONG_THRESH

    pep_rank = {p: i for i, p in enumerate(peptides)}
    all_rank = {a: i for i, a in enumerate(allele_order)}
    df["_p"] = df["peptide"].map(pep_rank)
    df["_a"] = df["allele"].map(all_rank)
    df = df.sort_values(["_p", "_a"]).drop(columns=["_p", "_a"])

    # Best allele per peptide, by lowest IC50.
    best = (df.loc[df.groupby("peptide")["ic50_nM"].idxmin(),
                   ["peptide", "allele", "ic50_nM"]]
              .rename(columns={"allele": "best_allele",
                               "ic50_nM": "best_ic50_nM"}))

    # Best allele per peptide, by lowest affinity percentile rank. IC50 is not
    # comparable across alleles, so these two criteria can disagree. Report both
    # rather than presenting one as if it were the answer.
    if df["affinity_percentile"].notna().any():
        bestp = (df.loc[df.groupby("peptide")["affinity_percentile"].idxmin(),
                        ["peptide", "allele", "affinity_percentile"]]
                   .rename(columns={"allele": "best_allele_percentile",
                                    "affinity_percentile": "best_percentile"}))
        best = best.merge(bestp, on="peptide", how="left")
    else:
        best["best_allele_percentile"] = np.nan
        best["best_percentile"] = np.nan
    best["criteria_agree"] = best["best_allele"] == best["best_allele_percentile"]

    # groupby sorts its keys alphabetically. Without this the summary block
    # prints in a different order from the matrix immediately above it.
    best["_p"] = best["peptide"].map(pep_rank)
    best = best.sort_values("_p").drop(columns=["_p"]).reset_index(drop=True)

    df = df.merge(best[["peptide", "best_allele"]], on="peptide", how="left")
    df["is_best_allele"] = df["allele"] == df["best_allele"]

    return df.reset_index(drop=True), best


def _clean(m):
    """Drop the 'peptide' / 'allele' axis names so to_string does not emit the
    two-row header pandas uses for named axes. Cosmetic, for the email block."""
    return m.rename_axis(index=None, columns=None)


def wide(df, value_col, peptides, allele_order):
    m = df.pivot(index="peptide", columns="allele", values=value_col)
    m = m.reindex(index=[p for p in peptides if p in m.index])
    m = m.reindex(columns=[a for a in allele_order if a in m.columns])
    return m


def email_table(df, ic50_m, pres_m, best, nested, unsupported):
    """Fixed-width text block, ready to paste into a reply."""
    L = []
    L.append("MHCflurry 2.x predicted MHC class I binding")
    L.append("")
    L.append("Predicted IC50 (nM), lower is stronger:")
    L.append("")
    L.append(_clean(ic50_m).round(1).to_string(
        float_format=lambda v: f"{v:,.1f}"))
    L.append("")
    has_pres = pres_m is not None and pres_m.notna().any().any()
    if has_pres:
        L.append("Presentation score (0-1, higher is more likely presented):")
        L.append("")
        L.append(_clean(pres_m).to_string(float_format=lambda v: f"{v:.3f}"))
        L.append("")
    L.append("Per-peptide detail:")
    L.append("")
    cols = ["peptide", "peptide_length", "allele", "ic50_nM",
            "affinity_percentile", "processing_score", "presentation_score",
            "binder_lt500", "strong_binder_lt50"]
    cols = [c for c in cols if c in df.columns]
    show = df[cols].copy()
    for c in ("ic50_nM",):
        show[c] = show[c].map(lambda v: f"{v:,.1f}")
    for c in ("affinity_percentile", "processing_score", "presentation_score"):
        if c in show.columns:
            show[c] = show[c].map(
                lambda v: "" if pd.isna(v) else f"{v:.3f}")
    L.append(show.to_string(index=False))
    L.append("")
    L.append("Strongest predicted allele per peptide (ranked by IC50):")
    disagree = []
    for _, r in best.iterrows():
        note = ""
        if not bool(r.get("criteria_agree", True)) and \
                pd.notna(r.get("best_allele_percentile", np.nan)):
            note = f"   [by percentile rank: {r['best_allele_percentile']}]"
            disagree.append(r["peptide"])
        L.append(f"  {r['peptide']:<12s} {r['best_allele']:<10s} "
                 f"{r['best_ic50_nM']:>10,.1f} nM{note}")
    L.append("")
    L.append("Notes:")
    L.append(f"  Thresholds shown are the conventional IC50 < {BIND_THRESH:,.0f} nM "
             f"(binder) and < {STRONG_THRESH:,.0f} nM (strong binder).")
    L.append("  IC50 is not directly comparable across alleles; the percentile "
             "rank is the")
    L.append("  allele-normalized value and is the better column for "
             "cross-allele comparison.")
    if disagree:
        L.append(f"  The two criteria disagree for: {', '.join(disagree)}. "
                 f"Both calls are shown above.")
    L.append("  All peptide-allele pairs are reported. Nothing was filtered.")

    # Calls sitting close enough to a threshold that the binary flag is fragile.
    borderline = df[df["ic50_nM"].between(0.8 * BIND_THRESH, 1.2 * BIND_THRESH) |
                    df["ic50_nM"].between(0.8 * STRONG_THRESH,
                                          1.2 * STRONG_THRESH)]
    if len(borderline):
        L.append("  Within 20% of a threshold, so the binary call is not "
                 "meaningful for:")
        for _, r in borderline.iterrows():
            L.append(f"    {r['peptide']:<12s} {r['allele']:<10s} "
                     f"{r['ic50_nM']:>10,.1f} nM")

    if has_pres:
        L.append("  The processing score is allele-independent by construction "
                 "(it models")
        L.append("  cleavage from sequence alone), so it is identical across "
                 "alleles for a given")
        L.append("  peptide and enters the allele-specific result only via the "
                 "presentation score.")
        L.append("  Peptides were submitted without flanking sequence, so the "
                 "processing and")
        L.append("  presentation columns are conservative relative to a run "
                 "with native flanks.")

    non9 = sorted({int(n) for n in df["peptide_length"].unique() if int(n) != 9})
    if non9:
        L.append(f"  MHCflurry training data are dominated by 9mers; "
                 f"{', '.join(f'{n}mer' for n in non9)} predictions carry wider "
                 f"uncertainty.")
    if nested:
        L.append("  Nested sequences in the submitted set:")
        for s, l in nested:
            L.append(f"    {s} is contained within {l}")
    if unsupported:
        L.append(f"  Alleles not recognized by MHCflurry and therefore not "
                 f"scored: {', '.join(unsupported)}")
    return "\n".join(L)


# =============================================================================
# FIGURE
# =============================================================================
def make_figure(df, peptides, allele_order, outdir, tag):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
    from matplotlib.patches import Rectangle

    plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42})

    ic50_m = wide(df, "ic50_nM", peptides, allele_order)
    pres_m = wide(df, "presentation_score", peptides, allele_order)
    has_pres = pres_m.notna().any().any()

    cmap_aff = LinearSegmentedColormap.from_list(
        "aff", [COLOR_STRONG, COLOR_MID, COLOR_WEAK])
    cmap_pres = LinearSegmentedColormap.from_list(
        "pres", [COLOR_WEAK, COLOR_MID, COLOR_STRONG])

    n_panels = 2 if has_pres else 1
    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(11 * n_panels + 3, 2.0 * len(ic50_m) + 5.5))
    axes = np.atleast_1d(axes)

    def draw(ax, mat, cmap, norm, title, fmt, outline_mask=None):
        im = ax.imshow(mat.values, cmap=cmap, norm=norm, aspect='auto')
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(mat.columns, fontsize=FS_TICK, rotation=30,
                           ha='right')
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index, fontsize=FS_TICK, family='monospace')
        ax.set_title(title, fontsize=FS_TITLE, pad=18, loc='left')
        ax.set_xticks(np.arange(-.5, mat.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-.5, mat.shape[0], 1), minor=True)
        ax.grid(which='minor', color='white', linewidth=2.5)
        ax.tick_params(which='minor', length=0)

        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.values[i, j]
                if pd.isna(v):
                    continue
                ax.text(j, i, fmt(v), ha='center', va='center',
                        fontsize=FS_ANNOT, color='#1A1A1A')
                if outline_mask is not None:
                    lw = outline_mask[i][j]
                    if lw:
                        ax.add_patch(Rectangle(
                            (j - .5, i - .5), 1, 1, fill=False,
                            edgecolor=COLOR_OUTLINE, linewidth=lw, zorder=5))

        cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
        cb.ax.tick_params(labelsize=FS_TICK)
        return cb

    # thicker outline for strong binders, thinner for binders
    outline = [[6 if v < STRONG_THRESH else (3 if v < BIND_THRESH else 0)
                for v in row] for row in ic50_m.values]

    cb = draw(axes[0], ic50_m, cmap_aff,
              LogNorm(vmin=1, vmax=50000),
              "a  predicted affinity (IC50, nM)",
              lambda v: f"{v:,.0f}", outline)
    cb.set_label("IC50 (nM)", fontsize=FS_LABEL, labelpad=14)

    if has_pres:
        # Fixed 0-1 scale. Autoscaling to the observed range would make a 0.86
        # and a 0.97 look identical and would mislabel the colorbar.
        cb2 = draw(axes[1], pres_m, cmap_pres, Normalize(vmin=0.0, vmax=1.0),
                   "b  presentation score",
                   lambda v: f"{v:.3f}")
        cb2.set_label("presentation score", fontsize=FS_LABEL, labelpad=14)

    fig.text(0.01, -0.02,
             f"Outlined cells: IC50 < {BIND_THRESH:,.0f} nM (binder); "
             f"heavy outline < {STRONG_THRESH:,.0f} nM (strong binder). "
             f"MHCflurry 2.x.",
             fontsize=FS_TICK, ha='left')

    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{tag}_heatmap.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
        log(f"  [SAVE] {os.path.basename(path)}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="MHC-I binding prediction for a literal peptide list.")
    ap.add_argument("--peptides", nargs="+", default=None,
                    help="peptide sequences (default: collaborator list)")
    ap.add_argument("--alleles", nargs="+", default=None,
                    help="HLA class I alleles (default: collaborator list)")
    ap.add_argument("--outdir", default=OUTPUT_DIR)
    ap.add_argument("--tag", default=DEFAULT_TAG,
                    help="output filename prefix")
    ap.add_argument("--table-only", action="store_true",
                    help="skip the figure (no matplotlib import)")
    ap.add_argument("--figure-from-tsv", default=None,
                    help="rebuild the figure from a previously written "
                         "{tag}_long.tsv; skips MHCflurry entirely")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- figure-only mode ---------------------------------------------------
    if args.figure_from_tsv:
        banner("FIGURE-ONLY MODE")
        df = pd.read_csv(args.figure_from_tsv, sep="\t")
        peptides = list(dict.fromkeys(df["peptide"].astype(str)))
        alleles = list(dict.fromkeys(df["allele"].astype(str)))
        log(f"  {len(df)} rows, {len(peptides)} peptides, "
            f"{len(alleles)} alleles")
        make_figure(df, peptides, alleles, args.outdir, args.tag)
        log("\n  Done.")
        return

    peptides_in = args.peptides if args.peptides else DEFAULT_PEPTIDES
    alleles_in = args.alleles if args.alleles else DEFAULT_ALLELES

    banner("INPUT")
    peptides, rejected = validate_peptides(peptides_in)
    log(f"  Submitted peptides: {len(peptides_in)}")
    for p in peptides:
        log(f"    {p:<14s} ({len(p)}mer)")
    if rejected:
        log("  Rejected:")
        for p, why in rejected:
            log(f"    {p:<14s} {why}")
    if not peptides:
        log("  FATAL: no usable peptides.")
        sys.exit(1)

    nested = report_nested(peptides)
    if nested:
        log("  Nested sequences detected:")
        for s, l in nested:
            log(f"    {s} is contained within {l}")

    log(f"  Submitted alleles: {', '.join(alleles_in)}")

    banner("AFFINITY PREDICTION (Class1AffinityPredictor)")
    aff, resolved, unsupported = run_affinity(peptides, alleles_in)

    banner("PRESENTATION PREDICTION (Class1PresentationPredictor)")
    pres = run_presentation(peptides, resolved)

    banner("ASSEMBLE")
    allele_order = [a for a in alleles_in if a in resolved]
    df, best = assemble(aff, pres, peptides, allele_order)
    ic50_m = wide(df, "ic50_nM", peptides, allele_order)
    pres_m = wide(df, "presentation_score", peptides, allele_order)

    long_path = os.path.join(args.outdir, f"{args.tag}_long.tsv")
    df.to_csv(long_path, sep="\t", index=False)
    log(f"  [SAVE] {os.path.basename(long_path)}  ({df.shape[0]} rows)")

    ic50_path = os.path.join(args.outdir, f"{args.tag}_ic50_matrix.tsv")
    ic50_m.round(2).to_csv(ic50_path, sep="\t")
    log(f"  [SAVE] {os.path.basename(ic50_path)}")

    if pres_m.notna().any().any():
        pres_path = os.path.join(args.outdir,
                                 f"{args.tag}_presentation_matrix.tsv")
        pres_m.round(4).to_csv(pres_path, sep="\t")
        log(f"  [SAVE] {os.path.basename(pres_path)}")

    txt = email_table(df, ic50_m, pres_m, best, nested, unsupported)
    txt_path = os.path.join(args.outdir, f"{args.tag}_table.txt")
    with open(txt_path, "w") as f:
        f.write(txt + "\n")
    log(f"  [SAVE] {os.path.basename(txt_path)}")

    banner("RESULT")
    log(txt)

    if not args.table_only:
        banner("FIGURE")
        try:
            make_figure(df, peptides, allele_order, args.outdir, args.tag)
        except ImportError as e:
            log(f"  matplotlib unavailable in this env ({e}).")
            log(f"  Re-run the figure in NETWORK env:")
            log(f"    conda run -n NETWORK python {os.path.basename(__file__)} "
                f"--figure-from-tsv {long_path} --tag {args.tag}")

    banner("COMPLETE")
    log(f"  Output directory: {args.outdir}")

    rep_path = os.path.join(args.outdir, f"{args.tag}_run_report.txt")
    with open(rep_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    print(f"\n  Report: {rep_path}")


if __name__ == "__main__":
    main()
