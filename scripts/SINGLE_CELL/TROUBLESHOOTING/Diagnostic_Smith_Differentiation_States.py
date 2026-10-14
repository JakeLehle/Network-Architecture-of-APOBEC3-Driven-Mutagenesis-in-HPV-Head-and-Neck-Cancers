#!/usr/bin/env python3
"""
Diagnostic_Smith_Differentiation_States.py
==========================================
2026_NMF_PAPER | Figure 3 diagnostic -- Smith et al. 2025 (EMBO J) keratinocyte
differentiation states scored across our basal epithelial compartment.

WHY THIS EXISTS
---------------
Tim Fenton's bioRxiv comment argues that the cells carrying our highest A3A and
SBS2 signal are what his group calls "differentiating keratinocytes", and that
our "basal cell" label understates them. His group defined four epithelial states
in healthy tonsil epithelium. This script scores those four states across our
basal compartment and asks two questions:

  Q1 (agreement, normal tissue)  In NORMAL-ADJACENT cells, do A3A-dominant cells
     sit in the differentiating / terminally differentiated states and A3B-dominant
     cells in the proliferating state, as Smith et al. report for healthy tonsil?

  Q2 (divergence, tumor tissue)  Does that same pairing hold in tumor cells, where
     HPV16 is present and E6/E7 have rewired the differentiation program?

TWO PANELS, BOTH RUN
--------------------
  CORE      the 12 marker genes named in their Figure 2 (3 per state). Clean,
            expression-based, defensible.
  COMBINED  CORE plus the transcription factors from their Figure 4.

  READ THE COMBINED PANEL WITH CARE. Their Figure 4 is a heatmap of SCENIC
  REGULON ACTIVITY, not marker expression, so averaging those genes' expression
  is NOT a reproduction of their analysis. Two specific hazards, both reported
  per-gene by this script so they can be checked rather than assumed:
    - the basal TF list carries EGR1, EGR3, ATF3, MYC, KLF4, CEBPD, which are
      classic tissue-dissociation stress genes;
    - the terminal TF list carries IRF7 and STAT2, which are interferon
      regulators, and our SBS2-HIGH population is defined in part by elevated
      type I interferon effectors. If the COMBINED terminal score separates our
      populations, check the per-gene table before believing it is differentiation.

CIRCULARITY NOTES (stated, not hidden)
--------------------------------------
  1. Tirosh phase overlap. MKI67 and TOP2A are in the Tirosh G2M gene set AND in
     Smith's proliferating panel. The proliferating score and the phase call are
     therefore NOT independent. Use one or the other as evidence, never both as
     mutual confirmation. The script prints the exact overlap.
  2. inferCNV reference. cnv_score used CytoTRACE2-defined normal cells as the
     reference, so normal-adjacent cells sit near baseline by construction. Low
     CNA in high-A3B normal cells is expected and is NOT independent evidence.
     The cell-cycle phase call does not use that reference and is clean.
  3. SPRR overlap. SPRR2D (terminal panel) sits in the CNA-HIGH cornification
     chain in Figure 5, wired to A3A. A leave-SPRR2D-out sensitivity score is
     computed alongside the full terminal score.

HARD STATE CALLS ARE MADE IN NORMAL-ADJACENT CELLS ONLY. Smith's states were
defined in healthy tonsil; in tumor the scores are shown as continuous overlays
because HPV rewires the program and an argmax label there would assert more than
the data supports.

PLUMBING matches Generate_Figure3_A3_Dominance_Slider.py exactly:
  basal = final_annotation == "basal cell"
  SBS2  = signature_weights_per_cell.txt row "SBS2" (0.0 where missing)
  CNV   = inferCNV cnv_score
  tissue split on source_name vs NORMAL_SOURCE (GEO typo intentional)

Run:
    conda run -n NETWORK python Diagnostic_Smith_Differentiation_States.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, chi2_contingency

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# ============================== CONFIG ==============================
PROJECT_ROOT     = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_FINAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH     = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")
FIG_DIR          = os.path.join(PROJECT_ROOT, "data/FIG_3/figures/TROUBLESHOOTING")
TSV_DIR          = os.path.join(PROJECT_ROOT, "data/FIG_3/00_diagnostics")

ANNOT_COL  = "final_annotation"
SOURCE_COL = "source_name"
CNV_COL    = "cnv_score"
CYTO_COL   = "CytoTRACE2_Score"
A3A_GENE   = "APOBEC3A"
A3B_GENE   = "APOBEC3B"
UMAP_KEY   = "X_umap"

TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE    = "normal tissue adjucent to head and neck squamous cell carcinoma"  # GEO typo intentional

# ---------------------------------------------------------------- gene panels
# CORE: named in Smith et al. Fig 2C/2D text.
CORE_PANELS = {
    "basal":          ["KRT14", "KRT15", "ITGB4"],
    "proliferating":  ["MKI67", "MCM7", "TOP2A"],
    "differentiating":["IVL", "KRT10", "KRT4"],
    "terminal":       ["SPRR2A", "S100P", "SPRR2D"],
}

# TF ADD-ON: Smith et al. Fig 4A, SCENIC regulon activity (NOT marker expression).
TF_PANELS = {
    "basal": ["NFIL3", "SOX6", "CEBPD", "IRF1", "EGR3", "SOX4", "KLF9", "ATF3",
              "MYC", "TEAD2", "NFIB", "KLF10", "ARID5B", "KLF4", "ERF", "RORB",
              "CREB5", "EGR1", "CHD2"],
    "proliferating": ["FOXM1"],
    "differentiating": ["GRHL3"],
    "terminal": ["LMO2", "ARNT", "ZNF101", "PPARD", "CEBPA", "ZNF860", "IRF7",
                 "NR1D1", "STAT2", "MXD1", "ELF3", "TEAD1", "NPAS2", "CD59",
                 "MAFF", "GRHL1", "PRDM1", "ELK3", "FOSL2", "PITX1", "EHF",
                 "FOXO3", "MEF2D"],
}

COMBINED_PANELS = {k: sorted(set(CORE_PANELS[k]) | set(TF_PANELS[k]))
                   for k in CORE_PANELS}

# genes flagged in the docstring, reported per-gene so they can be audited
FLAG_DISSOCIATION = {"EGR1", "EGR3", "ATF3", "MYC", "KLF4", "CEBPD", "FOSL2"}
FLAG_INTERFERON   = {"IRF7", "STAT2", "IRF1"}
FLAG_FIG5_CHAIN   = {"SPRR2D", "SPRR1A"}   # named in the CNA-HIGH cornification chain

STATE_ORDER  = ["basal", "proliferating", "differentiating", "terminal"]
STATE_LABELS = {"basal": "Basal", "proliferating": "Proliferating",
                "differentiating": "Differentiating", "terminal": "Term. diff."}
STATE_COLORS = {"basal": "#4c72b0", "proliferating": "#937860",
                "differentiating": "#55a868", "terminal": "#8172b3"}

# A3 dominance colors (match Diagnostic_Fig3_Dominance_Boxplots_SBS2_CNV.py)
COLOR_A3A = "#b5341f"
COLOR_A3B = "#3b6ea5"
# NMF palette (reserved: coral = SBS2/A3A pole, mustard = CNA/A3B pole)
COLOR_SBS2 = "#ed6a5a"
COLOR_CNV  = "#F6D155"

PHASE_ORDER  = ["G1", "S", "G2M"]
PHASE_COLORS = {"G1": "#9e9e9e", "S": "#7fa8c9", "G2M": "#c9737f"}

DOM_MARGIN = 0.0        # split A3A fraction at 0.5, matches the slider default
N_CTRL     = 50         # sc.tl.score_genes control bin size
SEED       = 0

# Boxplot overlay: boxes and the Mann-Whitney tests use ALL cells in a group.
# Overlaid points are a random subsample for legibility only (normal-adjacent
# groups are small enough that every cell is drawn). A diamond marks the group
# MEAN, computed on all cells.
N_JITTER = 1200
JIT_W    = 0.17
PT_ALPHA = 0.35
PT_SIZE  = 26

FS_TICK, FS_LAB, FS_TITLE, FS_ANNOT, FS_LEG = 28, 30, 32, 28, 26
FS_UMAP_STATE = 64      # state names over the UMAP grid (2x FS_TITLE, requested)
FS_PANEL_LTR  = 40      # a/b/c/d panel letters
DPI = 300

# UMAP rendering: matched to the published SBS2 / mutation / A3 panels in
# scripts/SINGLE_CELL/Step01_Generate_Figure3_Panels.py, which sets
# POINT_SIZE = 15 and CMAP_EXPRESSION = CMAP_MUTATIONS = 'plasma'. Keep these
# in step with that script so the diagnostic reads against the paper figures.
UMAP_POINT_SIZE = 15
UMAP_CMAP       = "plasma"

# Tirosh / Regev cell-cycle sets, as used by
# HPV_ANALYSIS/TROUBLESHOOTING/Diagnostic_Patient_CellCycle_by_Source.py.
# If patient_config.py carries edited lists, import those instead so the phase
# calls here match Supp. Fig. 6a exactly.
S_GENES = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
           "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS",
           "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2",
           "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6",
           "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
           "CHAF1B", "BRIP1", "E2F8"]
G2M_GENES = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
             "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
             "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
             "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
             "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
             "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
             "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
             "G2E3", "GAS2L3", "CBX5", "CENPA"]
# ====================================================================

RUN_TS = datetime.now().strftime("%Y%m%d_%H%M%S")
_loglines = []


def log(m=""):
    line = f"[{datetime.now():%H:%M:%S}] {m}"
    print(line, flush=True)
    _loglines.append(line)


def banner(t):
    log("")
    log("=" * 74)
    log(f"  {t}")
    log("=" * 74)


def save_fig(fig, name):
    for ext in ("pdf", "png"):
        p = os.path.join(FIG_DIR, f"{name}.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
    log(f"  [SAVE] {name}.pdf / .png")
    plt.close(fig)


def save_tsv(df, name):
    p = os.path.join(TSV_DIR, f"{name}.tsv")
    df.to_csv(p, sep="\t", index=False)
    log(f"  [SAVE] {name}.tsv  ({len(df)} rows)")


def get_expr(adata, gene):
    v = adata[:, gene].X
    return (v.toarray().flatten() if hasattr(v, "toarray")
            else np.asarray(v).flatten())


# =============================================================================
# BEAT 0: VERIFY ON-DISK STATE BEFORE COMPUTING ANYTHING
# =============================================================================
def beat0_verify():
    banner("BEAT 0: SCHEMA + PANEL VERIFICATION")

    for p in (ADATA_FINAL_PATH, WEIGHTS_PATH):
        if not os.path.exists(p):
            log(f"  ERROR: not found: {p}")
            sys.exit(1)
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(TSV_DIR, exist_ok=True)

    log(f"  Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for col in (ANNOT_COL, SOURCE_COL, CNV_COL, CYTO_COL):
        if col not in adata.obs.columns:
            numeric = [c for c in adata.obs.columns if adata.obs[c].dtype.kind in "fi"]
            log(f"  ERROR: '{col}' not in adata.obs. Numeric columns: {numeric}")
            sys.exit(1)
    log(f"  obs columns present: {ANNOT_COL}, {SOURCE_COL}, {CNV_COL}, {CYTO_COL}")

    if UMAP_KEY not in adata.obsm:
        log(f"  ERROR: '{UMAP_KEY}' not in adata.obsm. Available: {list(adata.obsm.keys())}")
        sys.exit(1)
    log(f"  obsm['{UMAP_KEY}'] present, shape {adata.obsm[UMAP_KEY].shape}")

    annots = sorted(adata.obs[ANNOT_COL].astype(str).unique())
    if TARGET_CELL_TYPE not in annots:
        log(f"  ERROR: '{TARGET_CELL_TYPE}' not in {ANNOT_COL}. Available: {annots}")
        sys.exit(1)
    if NORMAL_SOURCE not in set(adata.obs[SOURCE_COL].astype(str)):
        log(f"  ERROR: NORMAL_SOURCE not in {SOURCE_COL}.")
        log(f"    Available: {sorted(adata.obs[SOURCE_COL].astype(str).unique())}")
        sys.exit(1)

    # ---- panel coverage -----------------------------------------------------
    var = set(adata.var_names)
    rows = []
    log("")
    log("  Panel coverage (genes absent from var_names are DROPPED, not imputed):")
    for panel_name, panels in (("CORE", CORE_PANELS), ("COMBINED", COMBINED_PANELS)):
        for state in STATE_ORDER:
            genes = panels[state]
            present = [g for g in genes if g in var]
            missing = [g for g in genes if g not in var]
            log(f"    {panel_name:<9s} {state:<16s} {len(present):>2d}/{len(genes):<2d} present"
                + (f"   MISSING: {', '.join(missing)}" if missing else ""))
            for g in genes:
                rows.append({
                    "panel": panel_name, "state": state, "gene": g,
                    "present": g in var,
                    "source": "Fig2_core" if g in CORE_PANELS[state] else "Fig4_TF",
                    "flag_dissociation": g in FLAG_DISSOCIATION,
                    "flag_interferon": g in FLAG_INTERFERON,
                    "flag_fig5_chain": g in FLAG_FIG5_CHAIN,
                })
            if len(present) < 2:
                log(f"      WARNING: fewer than 2 genes present; this score will be unstable.")
    save_tsv(pd.DataFrame(rows), f"smith_panel_coverage_{RUN_TS}")

    # ---- circularity: Tirosh vs proliferating panel -------------------------
    tirosh = set(S_GENES) | set(G2M_GENES)
    for panel_name, panels in (("CORE", CORE_PANELS), ("COMBINED", COMBINED_PANELS)):
        ov = sorted(set(panels["proliferating"]) & tirosh)
        log("")
        log(f"  CIRCULARITY [{panel_name}] proliferating panel vs Tirosh phase genes: "
            f"{len(ov)} shared -> {', '.join(ov) if ov else 'none'}")
        if ov:
            log("    The proliferating score and the phase call share genes. Report one")
            log("    or the other as evidence, never both as mutual confirmation.")

    return adata


# =============================================================================
# BEAT 1: BUILD THE BASAL FRAME
# =============================================================================
def beat1_build(adata):
    banner("BEAT 1: BASAL SUBSET + PER-CELL METRICS")

    b = adata[adata.obs[ANNOT_COL].astype(str) == TARGET_CELL_TYPE].copy()
    log(f"  Basal cells: {b.n_obs:,}")

    if "SBS2" not in b.obs.columns:
        w = pd.read_csv(WEIGHTS_PATH, sep="\t", index_col=0)
        if "SBS2" not in w.index:
            log("  ERROR: SBS2 not a row in the weights matrix (needs .T on load?)")
            sys.exit(1)
        sbs2_map = w.loc["SBS2"].to_dict()
        b.obs["SBS2"] = b.obs_names.map(lambda x: sbs2_map.get(x, 0.0)).astype(float)
    b.obs["SBS2"] = b.obs["SBS2"].astype(float)

    for g in (A3A_GENE, A3B_GENE):
        if g not in b.var_names:
            log(f"  ERROR: {g} not in var_names")
            sys.exit(1)
        b.obs[g] = get_expr(b, g)
    b.obs[CNV_COL] = b.obs[CNV_COL].astype(float)
    b.obs[CYTO_COL] = b.obs[CYTO_COL].astype(float)

    b.obs["tissue_grp"] = np.where(
        b.obs[SOURCE_COL].astype(str) == NORMAL_SOURCE, "normal", "tumor")
    n_norm = int((b.obs["tissue_grp"] == "normal").sum())
    log(f"  Tumor: {b.n_obs - n_norm:,} | Normal-adjacent: {n_norm:,}")

    # A3 dominance, defined only where either enzyme is detected
    a3sum = b.obs[A3A_GENE] + b.obs[A3B_GENE]
    frac = np.where(a3sum > 0, b.obs[A3A_GENE] / a3sum.replace(0, np.nan), np.nan)
    b.obs["A3A_fraction"] = frac
    b.obs["A3_dominance"] = np.where(
        ~np.isfinite(frac), "none",
        np.where(frac > 0.5 + DOM_MARGIN, "A3A_dominant",
                 np.where(frac < 0.5 - DOM_MARGIN, "A3B_dominant", "tied")))
    for grp in ("A3A_dominant", "A3B_dominant", "tied", "none"):
        for t in ("tumor", "normal"):
            n = int(((b.obs["A3_dominance"] == grp) & (b.obs["tissue_grp"] == t)).sum())
            if n:
                log(f"    {grp:<14s} {t:<7s} n = {n:,}")

    # Tirosh phase (computed, not stored in obs; matches the Supp Fig 6a method)
    s_in = [g for g in S_GENES if g in b.var_names]
    g2_in = [g for g in G2M_GENES if g in b.var_names]
    log(f"  Tirosh sets present: S {len(s_in)}/{len(S_GENES)}, "
        f"G2M {len(g2_in)}/{len(G2M_GENES)}")
    sc.tl.score_genes_cell_cycle(b, s_genes=s_in, g2m_genes=g2_in)
    b.obs["phase"] = b.obs["phase"].astype(str)
    mix = b.obs["phase"].value_counts()
    log("  Basal phase mix: " + ", ".join(
        f"{p} {100 * mix.get(p, 0) / b.n_obs:.1f}%" for p in PHASE_ORDER))

    return b


# =============================================================================
# BEAT 2: SCORE THE PANELS
# =============================================================================
def beat2_score(b):
    banner("BEAT 2: DIFFERENTIATION STATE SCORES")

    var = set(b.var_names)
    score_cols = {}

    for panel_name, panels in (("CORE", CORE_PANELS), ("COMBINED", COMBINED_PANELS)):
        for state in STATE_ORDER:
            genes = [g for g in panels[state] if g in var]
            if not genes:
                log(f"  SKIP {panel_name}/{state}: no genes present")
                continue
            col = f"{panel_name}_{state}"
            sc.tl.score_genes(b, gene_list=genes, score_name=col,
                              ctrl_size=N_CTRL, random_state=SEED)
            # raw mean cross-check (no background subtraction)
            mat = b[:, genes].X
            mat = mat.toarray() if hasattr(mat, "toarray") else np.asarray(mat)
            b.obs[f"{col}_rawmean"] = mat.mean(axis=1)
            r = np.corrcoef(b.obs[col], b.obs[f"{col}_rawmean"])[0, 1]
            log(f"  {col:<26s} {len(genes):>2d} genes | score_genes vs raw mean r = {r:.3f}")
            score_cols.setdefault(panel_name, []).append(col)

    # sensitivity: terminal panels without SPRR2D (Fig 5 cornification chain)
    for panel_name, panels in (("CORE", CORE_PANELS), ("COMBINED", COMBINED_PANELS)):
        genes = [g for g in panels["terminal"] if g in var and g != "SPRR2D"]
        if genes:
            col = f"{panel_name}_terminal_noSPRR2D"
            sc.tl.score_genes(b, gene_list=genes, score_name=col,
                              ctrl_size=N_CTRL, random_state=SEED)
            r = np.corrcoef(b.obs[f"{panel_name}_terminal"], b.obs[col])[0, 1]
            log(f"  SENSITIVITY {col:<30s} vs full terminal r = {r:.3f}")
            if r < 0.90:
                log("    NOTE: dropping SPRR2D moves the score materially; report both.")

    # per-gene correlation with A3A / A3B so the flagged genes can be audited
    banner("BEAT 2b: PER-GENE AUDIT (flagged genes)")
    rows = []
    allg = sorted({g for p in COMBINED_PANELS.values() for g in p})
    for g in allg:
        if g not in var:
            continue
        x = get_expr(b, g)
        pct = 100.0 * (x > 0).mean()
        state = next(s for s in STATE_ORDER if g in COMBINED_PANELS[s])
        rows.append({
            "gene": g, "state": state,
            "source": "Fig2_core" if g in CORE_PANELS[state] else "Fig4_TF",
            "pct_cells_detected": round(pct, 2),
            "mean_expr": round(float(x.mean()), 4),
            "r_with_A3A": round(float(np.corrcoef(x, b.obs[A3A_GENE])[0, 1]), 4),
            "r_with_A3B": round(float(np.corrcoef(x, b.obs[A3B_GENE])[0, 1]), 4),
            "flag_dissociation": g in FLAG_DISSOCIATION,
            "flag_interferon": g in FLAG_INTERFERON,
            "flag_fig5_chain": g in FLAG_FIG5_CHAIN,
        })
    audit = pd.DataFrame(rows).sort_values(["state", "source", "gene"])
    save_tsv(audit, f"smith_pergene_audit_{RUN_TS}")

    thin = audit[audit["pct_cells_detected"] < 1.0]
    if len(thin):
        log(f"  {len(thin)} genes detected in under 1% of basal cells; these contribute")
        log(f"    almost nothing to a mean-based score: {', '.join(thin['gene'])}")
    for label, flag in (("dissociation", "flag_dissociation"),
                        ("interferon", "flag_interferon")):
        sub = audit[audit[flag]]
        if len(sub):
            log(f"  {label} flagged genes (check before trusting COMBINED):")
            for _, r in sub.iterrows():
                log(f"    {r['gene']:<9s} {r['state']:<16s} detected in "
                    f"{r['pct_cells_detected']:>5.1f}% of cells, "
                    f"r(A3A)={r['r_with_A3A']:+.3f} r(A3B)={r['r_with_A3B']:+.3f}")

    return b, score_cols


# =============================================================================
# BEAT 3: HARD STATE CALLS, NORMAL-ADJACENT ONLY
# =============================================================================
def beat3_assign(b):
    banner("BEAT 3: STATE ASSIGNMENT (NORMAL-ADJACENT ONLY)")

    for panel_name in ("CORE", "COMBINED"):
        cols = [f"{panel_name}_{s}" for s in STATE_ORDER if f"{panel_name}_{s}" in b.obs]
        if len(cols) < len(STATE_ORDER):
            log(f"  SKIP {panel_name}: not all four state scores available")
            continue
        # z-score each state score WITHIN normal-adjacent cells, then argmax.
        # Standardizing first stops a panel with more genes (and so a wider
        # score range) from winning the argmax on scale alone.
        nm = b.obs["tissue_grp"] == "normal"
        sub = b.obs.loc[nm, cols]
        z = (sub - sub.mean()) / sub.std(ddof=0)
        call = z.idxmax(axis=1).str.replace(f"{panel_name}_", "", regex=False)
        b.obs[f"state_{panel_name}"] = pd.Series("not_assigned", index=b.obs_names)
        b.obs.loc[nm, f"state_{panel_name}"] = call.values

        vc = call.value_counts()
        log(f"  {panel_name} state composition, normal-adjacent (n = {int(nm.sum()):,}):")
        for s in STATE_ORDER:
            n = int(vc.get(s, 0))
            log(f"    {STATE_LABELS[s]:<16s} {n:>4d}  ({100 * n / max(int(nm.sum()), 1):.1f}%)")
        log(f"    Smith et al. healthy tonsil for reference: Basal 1047, "
            f"Proliferating 353, Differentiating 1091, Term. diff. 158")
    return b


# =============================================================================
# BEAT 4: TESTS
# =============================================================================
def beat4_tests(b):
    banner("BEAT 4: A3 DOMINANCE vs STATE SCORES")

    rows = []
    for tissue in ("normal", "tumor"):
        tm = b.obs["tissue_grp"] == tissue
        a = b.obs.index[tm & (b.obs["A3_dominance"] == "A3A_dominant")]
        d = b.obs.index[tm & (b.obs["A3_dominance"] == "A3B_dominant")]
        log("")
        log(f"  {tissue.upper()}  A3A-dominant n = {len(a):,} | A3B-dominant n = {len(d):,}")
        if min(len(a), len(d)) < 20:
            log("    n is small; treat these as descriptive, not inferential.")
        for panel_name in ("CORE", "COMBINED"):
            for s in STATE_ORDER:
                col = f"{panel_name}_{s}"
                if col not in b.obs:
                    continue
                xa, xd = b.obs.loc[a, col].values, b.obs.loc[d, col].values
                if min(len(xa), len(xd)) < 3:
                    continue
                u, p = mannwhitneyu(xa, xd, alternative="two-sided")
                delta = float(np.median(xa) - np.median(xd))
                rows.append({"tissue": tissue, "panel": panel_name, "state": s,
                             "n_A3A_dom": len(xa), "n_A3B_dom": len(xd),
                             "median_A3A_dom": round(float(np.median(xa)), 4),
                             "median_A3B_dom": round(float(np.median(xd)), 4),
                             "median_delta": round(delta, 4),
                             "U": float(u), "p": float(p)})
                log(f"    {panel_name:<9s} {STATE_LABELS[s]:<16s} "
                    f"median delta (A3A minus A3B) {delta:+.4f}   p = {p:.2e}")
    tests = pd.DataFrame(rows)
    if len(tests):
        # BH within each tissue x panel family (four states per family)
        tests["q"] = np.nan
        for (t, pn), idx in tests.groupby(["tissue", "panel"]).groups.items():
            tests.loc[idx, "q"] = _bh(tests.loc[idx, "p"].values)
        save_tsv(tests, f"smith_state_vs_A3dominance_{RUN_TS}")

    # phase composition by dominance, normal-adjacent only
    banner("BEAT 4b: CELL-CYCLE PHASE BY A3 DOMINANCE, NORMAL-ADJACENT")
    nm = b.obs["tissue_grp"] == "normal"
    prows = []
    for grp in ("A3A_dominant", "A3B_dominant"):
        m = nm & (b.obs["A3_dominance"] == grp)
        n = int(m.sum())
        if n == 0:
            continue
        vc = b.obs.loc[m, "phase"].value_counts()
        pct = {p: 100.0 * vc.get(p, 0) / n for p in PHASE_ORDER}
        log(f"  {grp:<14s} n = {n:>4d}  " +
            "  ".join(f"{p} {pct[p]:.1f}%" for p in PHASE_ORDER))
        prows.append({"group": grp, "n": n, **{p: round(pct[p], 2) for p in PHASE_ORDER}})
    if len(prows) == 2:
        ct = np.array([[int(((nm) & (b.obs["A3_dominance"] == g) &
                             (b.obs["phase"] == p)).sum())
                        for p in PHASE_ORDER] for g in ("A3A_dominant", "A3B_dominant")])
        if ct.sum() > 0 and (ct.sum(axis=1) > 0).all():
            chi2, pv, dof, _ = chi2_contingency(ct)
            log(f"  chi-square {chi2:.2f}, dof {dof}, p = {pv:.2e}")
            log("  Quote percentage-point differences, not the p-value.")
    if prows:
        save_tsv(pd.DataFrame(prows), f"smith_phase_by_dominance_normal_{RUN_TS}")

    # CNA in normal-adjacent, reported with its circularity warning
    banner("BEAT 4c: CNA IN NORMAL-ADJACENT BY A3 DOMINANCE (circular, see header)")
    for grp in ("A3A_dominant", "A3B_dominant"):
        m = nm & (b.obs["A3_dominance"] == grp)
        if int(m.sum()):
            v = b.obs.loc[m, CNV_COL]
            log(f"  {grp:<14s} n = {int(m.sum()):>4d}  cnv_score median "
                f"{v.median():.4f}  mean {v.mean():.4f}")
    log("  inferCNV used CytoTRACE2-defined normal cells as reference, so these")
    log("  sit near baseline by construction. Descriptive only.")

    return b


# =============================================================================
# BEAT 5: FIGURES
# =============================================================================
def fig_umap_scores(b):
    banner("BEAT 5a: UMAP SCORE OVERLAYS")
    xy = b.obsm[UMAP_KEY]
    fig, axes = plt.subplots(2, 4, figsize=(34, 18.5))
    for ri, panel_name in enumerate(("CORE", "COMBINED")):
        for ci, s in enumerate(STATE_ORDER):
            ax = axes[ri, ci]
            col = f"{panel_name}_{s}"
            if col not in b.obs:
                ax.set_visible(False)
                continue
            v = b.obs[col].values
            lo, hi = np.percentile(v, [1, 99])
            order = np.argsort(v)
            sct = ax.scatter(xy[order, 0], xy[order, 1], c=v[order],
                             s=UMAP_POINT_SIZE, cmap=UMAP_CMAP,
                             vmin=lo, vmax=hi, linewidths=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if ri == 0:
                ax.set_title(STATE_LABELS[s], fontsize=FS_UMAP_STATE, pad=18)
            if ci == 0:
                ax.set_ylabel(panel_name, fontsize=FS_LAB)
            cb = fig.colorbar(sct, ax=ax, fraction=0.045, pad=0.02)
            cb.ax.tick_params(labelsize=FS_TICK - 6)
    fig.tight_layout()
    save_fig(fig, f"Smith_States_UMAP_Scores_{RUN_TS}")


def fig_umap_reference(b):
    banner("BEAT 5b: UMAP REFERENCE CHANNELS")
    xy = b.obsm[UMAP_KEY]
    channels = [A3A_GENE, A3B_GENE, "SBS2", CNV_COL]
    fig, axes = plt.subplots(1, 4, figsize=(34, 9))
    for ax, col in zip(axes, channels):
        v = b.obs[col].values.astype(float)
        lo, hi = np.percentile(v, [1, 99])
        if hi <= lo:
            hi = lo + 1e-6
        order = np.argsort(v)
        sct = ax.scatter(xy[order, 0], xy[order, 1], c=v[order],
                         s=UMAP_POINT_SIZE, cmap=UMAP_CMAP,
                         vmin=lo, vmax=hi, linewidths=0, rasterized=True)
        ax.set_title(col, fontsize=FS_TITLE, pad=12)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        cb = fig.colorbar(sct, ax=ax, fraction=0.045, pad=0.02)
        cb.ax.tick_params(labelsize=FS_TICK - 6)
    fig.tight_layout()
    save_fig(fig, f"Smith_States_UMAP_Reference_{RUN_TS}")


def fig_normal_states(b):
    banner("BEAT 5c: NORMAL-ADJACENT STATE COMPOSITION")
    nm = b.obs["tissue_grp"] == "normal"
    panels = [p for p in ("CORE", "COMBINED") if f"state_{p}" in b.obs]
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(13 * len(panels), 11))
    axes = np.atleast_1d(axes)
    for ax, panel_name in zip(axes, panels):
        col = f"state_{panel_name}"
        groups = ["A3A_dominant", "A3B_dominant"]
        labels, bottoms = [], []
        data = {}
        for g in groups:
            m = nm & (b.obs["A3_dominance"] == g)
            n = int(m.sum())
            vc = b.obs.loc[m, col].value_counts()
            data[g] = {s: 100.0 * vc.get(s, 0) / max(n, 1) for s in STATE_ORDER}
            labels.append(f"{'A3A' if g.startswith('A3A') else 'A3B'}-dominant\n(n = {n:,})")
        x = np.arange(len(groups))
        run = np.zeros(len(groups))
        for s in STATE_ORDER:
            vals = np.array([data[g][s] for g in groups])
            ax.bar(x, vals, bottom=run, color=STATE_COLORS[s], width=0.6,
                   edgecolor="white", linewidth=1.5)
            for xi, (v, r) in enumerate(zip(vals, run)):
                if v >= 7:
                    ax.text(xi, r + v / 2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=FS_ANNOT, fontweight="bold", color="white")
            run += vals
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=FS_TICK)
        ax.set_ylabel("Percent of cells", fontsize=FS_LAB)
        ax.set_title(panel_name, fontsize=FS_TITLE, pad=14)
        ax.tick_params(axis="y", labelsize=FS_TICK)
        ax.set_ylim(0, 100)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    handles = [Patch(facecolor=STATE_COLORS[s], label=STATE_LABELS[s])
               for s in STATE_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=FS_LEG,
               frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    save_fig(fig, f"Smith_States_Normal_Composition_{RUN_TS}")


def _bh(pvals):
    """
    Benjamini-Hochberg q-values for a 1-D array of p-values.

    The step-up enforcement (running minimum from the largest rank down) MUST be
    applied while the values are in rank order, and only then scattered back to
    the caller's ordering. Applying it to values already scattered back mixes up
    the ranks and collapses every q to the single smallest adjusted value, which
    is what made all four states report an identical q.
    """
    pv = np.asarray(pvals, dtype=float)
    n = len(pv)
    if n == 0:
        return pv
    order = np.argsort(pv)
    adj = pv[order] * n / (np.arange(n) + 1)     # rank-ordered adjustment
    adj = np.minimum.accumulate(adj[::-1])[::-1]  # step-up, still rank-ordered
    q = np.empty(n, dtype=float)
    q[order] = np.clip(adj, 0.0, 1.0)             # scatter back last
    return q


def _q_label(q):
    if not np.isfinite(q):
        return "n.s."
    if q >= 0.05:
        return f"n.s. (q = {q:.2f})"
    if q < 1e-4:
        return f"q = {q:.1e}"
    return f"q = {q:.4f}".rstrip("0")


def fig_score_boxes(b):
    """
    Panel a is the figure Jake sends: CORE scores in NORMAL-ADJACENT cells,
    A3A-dominant versus A3B-dominant. This is the direct agreement with
    Smith et al. Panels b-d carry the COMBINED panel and the tumor contrast.

    Boxes, mean diamonds and the Mann-Whitney tests use EVERY cell in a group.
    Overlaid points are a random subsample (N_JITTER) for legibility only; the
    normal-adjacent groups are smaller than that cap, so every normal cell is
    drawn. Mann-Whitney is two-sided, BH-corrected across the four states
    within each panel and tissue.
    """
    banner("BEAT 5d: STATE SCORES BY A3 DOMINANCE, TUMOR vs NORMAL")
    rng = np.random.default_rng(SEED)
    panels = ("CORE", "COMBINED")
    tissues = ("normal", "tumor")
    letters = ["a", "b", "c", "d"]

    fig, axes = plt.subplots(2, 2, figsize=(28, 22))
    ann_rows = []
    li = 0

    for ri, panel_name in enumerate(panels):
        # shared y-limits across the row so normal and tumor are comparable
        row_vals = []
        for tissue in tissues:
            tm = b.obs["tissue_grp"] == tissue
            dm = b.obs["A3_dominance"].isin(["A3A_dominant", "A3B_dominant"])
            for s in STATE_ORDER:
                col = f"{panel_name}_{s}"
                if col in b.obs:
                    row_vals.append(b.obs.loc[tm & dm, col].values)
        allv = np.concatenate(row_vals) if row_vals else np.array([0.0, 1.0])
        ylo, yhi = np.percentile(allv, [0.5, 99.5])
        span = max(yhi - ylo, 1e-6)
        ylo -= 0.08 * span
        ytop = yhi + 0.30 * span          # headroom for the test brackets

        for ci, tissue in enumerate(tissues):
            ax = axes[ri, ci]
            tm = b.obs["tissue_grp"] == tissue

            ticks, ticklabels, pending = [], [], []
            for si, s in enumerate(STATE_ORDER):
                col = f"{panel_name}_{s}"
                if col not in b.obs:
                    continue
                grp_vals = {}
                for oi, (grp, color) in enumerate((("A3A_dominant", COLOR_A3A),
                                                   ("A3B_dominant", COLOR_A3B))):
                    m = tm & (b.obs["A3_dominance"] == grp)
                    vals = b.obs.loc[m, col].values.astype(float)
                    grp_vals[grp] = vals
                    if len(vals) < 3:
                        continue
                    xpos = si * 3 + oi

                    bp = ax.boxplot([vals], positions=[xpos], widths=0.66,
                                    patch_artist=True, showfliers=False,
                                    zorder=3)
                    bp["boxes"][0].set_facecolor(color)
                    bp["boxes"][0].set_alpha(0.55)
                    bp["boxes"][0].set_edgecolor("#000000")
                    bp["boxes"][0].set_linewidth(2.0)
                    for el in ("whiskers", "caps"):
                        for ln in bp[el]:
                            ln.set_color("#000000")
                            ln.set_linewidth(2.0)
                    for ln in bp["medians"]:
                        ln.set_color("#000000")
                        ln.set_linewidth(3.0)

                    # jittered points (display subsample only)
                    if len(vals) > N_JITTER:
                        draw = vals[rng.choice(len(vals), N_JITTER, replace=False)]
                    else:
                        draw = vals
                    jx = xpos + rng.uniform(-JIT_W, JIT_W, size=len(draw))
                    ax.scatter(jx, draw, s=PT_SIZE, color=color, alpha=PT_ALPHA,
                               linewidths=0, zorder=2, rasterized=True)

                    # mean diamond, computed on ALL cells in the group
                    ax.scatter([xpos], [vals.mean()], marker="D", s=190,
                               facecolor="#ffffff", edgecolor="#000000",
                               linewidths=2.2, zorder=5)

                a = grp_vals.get("A3A_dominant", np.array([]))
                d = grp_vals.get("A3B_dominant", np.array([]))
                if len(a) >= 3 and len(d) >= 3:
                    u, p = mannwhitneyu(a, d, alternative="two-sided")
                    pending.append({
                        "panel": panel_name, "tissue": tissue, "state": s,
                        "x_left": si * 3, "x_right": si * 3 + 1,
                        "top": float(max(np.percentile(a, 97.5),
                                         np.percentile(d, 97.5))),
                        "n_A3A_dom": int(len(a)), "n_A3B_dom": int(len(d)),
                        "median_A3A_dom": float(np.median(a)),
                        "median_A3B_dom": float(np.median(d)),
                        "median_delta": float(np.median(a) - np.median(d)),
                        "U": float(u), "p": float(p),
                    })
                ticks.append(si * 3 + 0.5)
                ticklabels.append(STATE_LABELS[s])

            # BH across the four states in this panel x tissue
            if pending:
                qs = _bh([r["p"] for r in pending])
                for r, q in zip(pending, qs):
                    r["q"] = float(q)
                ann_rows.extend(pending)

                bar_y = yhi + 0.10 * span
                step = 0.055 * span
                for k, r in enumerate(pending):
                    y = bar_y + (k % 2) * step
                    ax.plot([r["x_left"], r["x_left"], r["x_right"], r["x_right"]],
                            [y - 0.018 * span, y, y, y - 0.018 * span],
                            color="#000000", lw=2.0, zorder=6,
                            solid_capstyle="butt")
                    ax.text((r["x_left"] + r["x_right"]) / 2, y + 0.012 * span,
                            _q_label(r["q"]), ha="center", va="bottom",
                            fontsize=FS_ANNOT - 6, zorder=6)

            ax.set_xticks(ticks)
            ax.set_xticklabels(ticklabels, fontsize=FS_TICK - 2,
                               rotation=18, ha="right")
            ax.tick_params(axis="y", labelsize=FS_TICK)
            ax.set_ylabel(f"{panel_name} score", fontsize=FS_LAB)
            ax.set_ylim(ylo, ytop)
            ax.set_xlim(-1.0, (len(STATE_ORDER) - 1) * 3 + 2.0)

            n_a = int((tm & (b.obs["A3_dominance"] == "A3A_dominant")).sum())
            n_b = int((tm & (b.obs["A3_dominance"] == "A3B_dominant")).sum())
            ax.set_title(
                f"{'Normal-adjacent' if tissue == 'normal' else 'Tumor'}   "
                f"A3A-dominant n = {n_a:,}, A3B-dominant n = {n_b:,}",
                fontsize=FS_TITLE - 6, pad=14)
            ax.text(-0.11, 1.06, letters[li], transform=ax.transAxes,
                    fontsize=FS_PANEL_LTR, fontweight="bold", va="top", ha="left")
            li += 1
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)

    handles = [Patch(facecolor=COLOR_A3A, alpha=0.55, edgecolor="#000000",
                     label="A3A-dominant"),
               Patch(facecolor=COLOR_A3B, alpha=0.55, edgecolor="#000000",
                     label="A3B-dominant")]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=FS_LEG,
               frameon=False, bbox_to_anchor=(0.5, -0.025))
    fig.tight_layout()
    save_fig(fig, f"Smith_States_Boxes_by_Dominance_{RUN_TS}")

    if ann_rows:
        df = pd.DataFrame(ann_rows)[
            ["panel", "tissue", "state", "n_A3A_dom", "n_A3B_dom",
             "median_A3A_dom", "median_A3B_dom", "median_delta", "U", "p", "q"]]
        # round the medians only; rounding p and q to 6 dp flattens anything
        # below 1e-6 to zero, which destroys the numbers the figure quotes
        for c in ("median_A3A_dom", "median_A3B_dom", "median_delta"):
            df[c] = df[c].round(4)
        save_tsv(df, f"smith_boxplot_tests_{RUN_TS}")
        log("")
        log("  Panel a source (CORE, normal-adjacent), Mann-Whitney two-sided, "
            "BH across the four states:")
        for _, r in df[(df["panel"] == "CORE") &
                       (df["tissue"] == "normal")].iterrows():
            log(f"    {STATE_LABELS[r['state']]:<16s} "
                f"A3A-dom median {r['median_A3A_dom']:+.3f} "
                f"(n = {int(r['n_A3A_dom'])})   "
                f"A3B-dom median {r['median_A3B_dom']:+.3f} "
                f"(n = {int(r['n_A3B_dom'])})   "
                f"delta {r['median_delta']:+.3f}   "
                f"p = {r['p']:.2e}   {_q_label(r['q'])}")
        log("  Each row is a separate two-sided Mann-Whitney U of A3A-dominant "
            "versus A3B-dominant")
        log("  cells within that one state. BH is applied across the four states "
            "only.")


# =============================================================================
def main():
    adata = beat0_verify()
    b = beat1_build(adata)
    del adata
    b, _ = beat2_score(b)
    b = beat3_assign(b)
    b = beat4_tests(b)

    fig_umap_scores(b)
    fig_umap_reference(b)
    fig_normal_states(b)
    fig_score_boxes(b)

    banner("BEAT 6: PER-CELL EXPORT")
    keep = ([ANNOT_COL, SOURCE_COL, "tissue_grp", A3A_GENE, A3B_GENE,
             "A3A_fraction", "A3_dominance", "SBS2", CNV_COL, CYTO_COL, "phase"]
            + [c for c in b.obs.columns
               if c.startswith(("CORE_", "COMBINED_", "state_"))])
    out = b.obs[[c for c in keep if c in b.obs.columns]].copy()
    out.insert(0, "cell_barcode", out.index)
    save_tsv(out, f"smith_state_scores_percell_{RUN_TS}")

    logp = os.path.join(TSV_DIR, f"smith_differentiation_states_{RUN_TS}.txt")
    with open(logp, "w") as fh:
        fh.write("\n".join(_loglines) + "\n")
    log(f"  [SAVE] {os.path.basename(logp)}")
    banner("DONE")


if __name__ == "__main__":
    main()
