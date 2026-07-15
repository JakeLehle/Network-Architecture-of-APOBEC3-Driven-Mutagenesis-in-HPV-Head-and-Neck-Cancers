#!/usr/bin/env python3
"""
Generate_DEG_KEGG_Supplement.py  (v7)
=====================================

Supplemental figure bridging single-cell DEG analysis to the differential
network analysis, tying DEGs to the A3 program.

Layout (2x2), one comparison per column:
                SBS2-HIGH vs NORMAL (coral)   |   CNV-HIGH vs NORMAL (mustard)
  top row:      diverging KEGG bar; UP pathways on top, DOWN pathways on bottom;
                pathway names as full y-tick labels
  bottom row:   volcano of ALL DEGs. Background points: light gray below
                FDR<0.05, darker gray above. Story-pathway significant genes
                colored (HPV coral, Cell cycle purple, Antigen presentation
                blue). Harris A3 interactors as gold dots; A3A/A3B as dark-green
                diamonds. Titles carry up/down counts. Labels stacked in edge
                columns with leader arrows. Per-panel legend on the right side.

v3 changes (from review):
  - KEGG order flipped (up on top).
  - Volcano story pathways = HPV (coral, SBS2 up), Cell cycle (purple, CNV up),
    Antigen presentation (blue, CNV down); HIV-1 removed.
  - Two-tone gray by FDR<0.05; threshold line drawn + labeled.
  - Interactors -> gold dots; A3A/A3B -> dark-green diamonds.
  - Up/down counts in each volcano title.
  - Legends moved to the side of each panel.

Reads (READ-ONLY):
  data/FIG_4/NETWORK_<COMP>/02_differential_expression/SC_diffexpr_stats.csv
  data/FIG_4/NETWORK_<COMP>/02_differential_expression/KEGG/KEGG_<COMP>_{up,down}.tsv
  Harris A3 interactor list (via network_config_SC.load_harris_interactors)

Writes:
  data/FIG_4/DIAGNOSTIC_DEG_KEGG_SUPPLEMENT/Supp_DEG_KEGG.{pdf,png}   (300 DPI)

Placement: scripts/NETWORK_SINGLE_CELL/TROUBLESHOOTING/  (excluded from walkthrough)
Env:       conda run -n NETWORK python Generate_DEG_KEGG_Supplement.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from network_config_SC import (
    BASE_DIR, FIG4_ROOT,
    A3_SYMBOL_TO_ALIAS, FDR_THRESHOLD,
    ensure_dir, load_harris_interactors,
)

# =============================================================================
# ADJUSTABLE KNOBS
# =============================================================================
OUTPUT_DIR   = os.path.join(FIG4_ROOT, "DIAGNOSTIC_DEG_KEGG_SUPPLEMENT")
OUTFILE_STEM = "Supp_DEG_KEGG"

N_PATHWAY_LABELS      = 25     # gene names labeled per volcano (pathway members, by FDR)
N_PATHWAYS_PER_DIR    = 8      # KEGG pathways shown per direction (up / down)
N_INTERACTOR_LABELS   = 15     # interactor names labeled per volcano (by FDR)
ALWAYS_LABEL          = ["APOBEC3A", "APOBEC3B"]     # enzymes always named
LABEL_ALL_INTERACTORS = False
CIRCLE_ONLY_SIG_INTERACTORS = False   # False = show every interactor present
EXPORT_SIG_ONLY = True    # DEG tables: True = significant only (FDR<thr); False = all tested genes
COLOR_FULL_PATHWAY = True   # color ALL members of a story pathway (up+down, sig+non-sig); False = enriched-direction hits only
RANK_PATHWAY_LABELS_BY = "pct_test"   # label ranking: "pct_test" (expression prevalence), "fdr", or "abs_log2FC"

# Per-pathway label ranking: (metric, direction). direction "desc" labels the
# highest-metric genes first, "asc" the lowest first. metric is a DEG column
# (pct_test = fraction of HIGH-group cells expressing the gene). This is what
# makes CNV differ: Cell cycle (up) labels its highest-expressed genes, while
# Antigen presentation (down) labels its lowest-expressed (most suppressed) ones,
# so the down story is not crowded out by the high-expression up genes.
PATHWAY_LABEL_RANK = {
    "Human papillomavirus infection":       ("pct_test", "desc"),
    "Cell cycle":                           ("pct_test", "desc"),
    "Antigen processing and presentation":  ("pct_test", "asc"),
}
# Interactor label ranking per comparison (CNV groups with Cell cycle: expression
# high-first; SBS2 keeps significance ordering).
INTERACTOR_LABEL_RANK = {
    "SBS2_VS_NORMAL": ("fdr", "asc"),
    "CNV_VS_NORMAL":  ("pct_test", "desc"),
}

# Fonts. Structural text stays in the 28-34 range; KEGG names and gene labels
# are smaller by necessity (many labels), both tunable here.
FONT_TITLE, FONT_AXIS, FONT_TICK, FONT_LEGEND = 34, 30, 28, 22
FONT_VOLC_TITLE = 30    # volcano title carries counts, kept slightly smaller to fit
FONT_BAR        = 20    # KEGG pathway y-tick labels
FONT_GENE_LABEL = 18    # volcano gene labels (deliberate deviation from 28-34)

FIGSIZE = (30, 26)
DPI = 300

HUE = {"SBS2_VS_NORMAL": "#ed6a5a", "CNV_VS_NORMAL": "#F6D155"}
LABELS = {"SBS2_VS_NORMAL": "SBS2-HIGH vs NORMAL",
          "CNV_VS_NORMAL":  "CNV-HIGH vs NORMAL"}
COMPARISONS = ["SBS2_VS_NORMAL", "CNV_VS_NORMAL"]

# Story pathways per column (HIV-1 removed).
INTERESTING_PATHWAYS = {
    "SBS2_VS_NORMAL": {"up": ["Human papillomavirus infection"], "down": []},
    "CNV_VS_NORMAL":  {"up": ["Cell cycle"],
                       "down": ["Antigen processing and presentation"]},
}

# Fixed color per story pathway.
PATHWAY_COLOR = {
    "Human papillomavirus infection":            "#ed6a5a",   # coral
    "Cell cycle":                                "#8d5a97",   # purple
    "Antigen processing and presentation":       "#2b6cb0",   # blue
}

# Point colors
GRAY_NS       = "#d9d9d9"   # below FDR<0.05
GRAY_SIG      = "#7a7f87"   # above FDR<0.05 (darker for emphasis)
INTERACTOR_COLOR = "#E8B923"   # gold
INTERACTOR_LABEL_COLOR = "#8a6d00"
A3_COLOR      = "#1b4332"   # dark green
A3_ENZYMES    = ["APOBEC3A", "APOBEC3B"]


def log(msg=""):
    print(msg, flush=True)


def mute(hexc, f=0.45):
    r, g, b = int(hexc[1:3], 16), int(hexc[3:5], 16), int(hexc[5:7], 16)
    r = int(r + (150 - r) * f); g = int(g + (150 - g) * f); b = int(b + (150 - b) * f)
    return f"#{r:02x}{g:02x}{b:02x}"


def deg_path(comp):
    return os.path.join(FIG4_ROOT, f"NETWORK_{comp}",
                        "02_differential_expression", "SC_diffexpr_stats.csv")


def kegg_path(comp, direction):
    return os.path.join(FIG4_ROOT, f"NETWORK_{comp}",
                        "02_differential_expression", "KEGG",
                        f"KEGG_{comp}_{direction}.tsv")


def load_kegg(comp, direction):
    p = kegg_path(comp, direction)
    if not os.path.exists(p):
        log(f"  [WARN] missing KEGG file: {p}")
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def pathway_hit_genes(kegg_df, pathway_names):
    out = {}
    if len(kegg_df) == 0:
        return out
    for _, r in kegg_df.iterrows():
        term = str(r.get("Term", ""))
        if term in pathway_names:
            genes = [g.strip().upper() for g in str(r.get("Genes", "")).split(";")
                     if g.strip()]
            out[term] = (float(r.get("Adjusted P-value", 1.0)), genes)
    return out


def get_pathway_members(comp, pathways):
    """Return ({pathway: set(member genes)}, source_tag).
    Primary (COLOR_FULL_PATHWAY): the FULL KEGG_2021_Human membership for each
    pathway (every annotated gene, both directions). Fallback / disabled: the
    union of significant up+down KEGG hits from the TSVs."""
    if COLOR_FULL_PATHWAY:
        try:
            import gseapy as gp
            lib = gp.get_library(name="KEGG_2021_Human")
            return ({p: {str(x).upper() for x in lib.get(p, [])} for p in pathways},
                    "full_KEGG_library")
        except Exception as e:
            log(f"  [WARN] KEGG library fetch failed ({e}); using significant hits only.")
    out = {p: set() for p in pathways}
    for direction in ("up", "down"):
        hits = pathway_hit_genes(load_kegg(comp, direction), set(pathways))
        for term, (_a, genes) in hits.items():
            out.setdefault(term, set()).update(genes)
    return out, "significant_hits"


# =============================================================================
# TOP ROW: diverging KEGG bar (UP on top, DOWN on bottom; full y-tick names)
# =============================================================================
def draw_kegg_bar(ax, comp):
    hue = HUE[comp]
    up, down = load_kegg(comp, "up"), load_kegg(comp, "down")

    def top(df, n):
        if len(df) == 0:
            return []
        d = df.sort_values("Adjusted P-value").head(n)
        return list(zip(d["Term"].astype(str), d["Adjusted P-value"].astype(float)))

    entries = []  # (term, signed_x, color, adjp) -- UP first so it lands on top
    for term, adjp in top(up, N_PATHWAYS_PER_DIR):
        entries.append((term, +(-np.log10(max(adjp, 1e-300))), hue, adjp))
    for term, adjp in top(down, N_PATHWAYS_PER_DIR):
        entries.append((term, -(-np.log10(max(adjp, 1e-300))), mute(hue), adjp))

    positions = list(range(len(entries) - 1, -1, -1))
    xs = [e[1] for e in entries] + [0.0]
    xmin, xmax = min(xs), max(xs)
    span = (xmax - xmin) or 1.0
    pad = 0.015 * span
    for pos, (term, x, color, adjp) in zip(positions, entries):
        ax.barh(pos, x, color=color, edgecolor="black", linewidth=0.6, height=0.72)
        # adjusted p-value at the bar tip, pointing outward
        tx = x + pad if x >= 0 else x - pad
        ha = "left" if x >= 0 else "right"
        ax.text(tx, pos, f"{adjp:.1e}", ha=ha, va="center",
                fontsize=FONT_BAR - 3, color="#333333")
    # leave room for the adjP text on both sides
    ax.set_xlim(xmin - 0.22 * span, xmax + 0.22 * span)

    ax.axvline(0, color="black", linewidth=1.2)
    ax.set_yticks(positions)
    ax.set_yticklabels([e[0] for e in entries], fontsize=FONT_BAR)
    ax.set_xlabel("signed  -log10(adj P)      (down  |  up)", fontsize=FONT_AXIS)
    ax.set_title(f"{LABELS[comp]}: KEGG enrichment",
                 fontsize=FONT_TITLE, color=hue, fontweight="bold", pad=14)
    ax.tick_params(axis="x", labelsize=FONT_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# =============================================================================
# Deterministic label placement: edge columns + leader arrows (no overlap)
# =============================================================================
def place_labels(ax, rows):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    xspan, yspan = xlim[1] - xlim[0], ylim[1] - ylim[0]
    ylo, yhi = ylim[0] + 0.04 * yspan, ylim[1] - 0.04 * yspan

    for items, xpos, ha in (
        ([r for r in rows if r["px"] < 0],  xlim[0] + 0.02 * xspan, "left"),
        ([r for r in rows if r["px"] >= 0], xlim[1] - 0.02 * xspan, "right"),
    ):
        items = sorted(items, key=lambda r: r["py"])
        n = len(items)
        if n == 0:
            continue
        if n == 1:
            ypositions = [(ylo + yhi) / 2]
        else:
            ypositions = [ylo + i * (yhi - ylo) / (n - 1) for i in range(n)]
        for r, ly in zip(items, ypositions):
            ax.annotate(r["text"], xy=(r["px"], r["py"]), xytext=(xpos, ly),
                        ha=ha, va="center", fontsize=r["size"],
                        fontweight=r["weight"], color=r["color"], zorder=9,
                        arrowprops=dict(arrowstyle="-", color="#b0b0b0",
                                        lw=0.6, shrinkA=2, shrinkB=3))


# =============================================================================
# BOTTOM ROW: volcano
# =============================================================================
def draw_volcano(ax, comp, interactors):
    hue = HUE[comp]
    df = pd.read_csv(deg_path(comp))
    df["gene"] = df["gene"].astype(str)
    df = df[~df["gene"].str.startswith("ENSG")].copy()
    df["nlp"] = -np.log10(df["p_value"].clip(lower=1e-300))
    sig = df["fdr"] < FDR_THRESHOLD
    n_up = int((sig & (df["log2FC"] > 0)).sum())
    n_down = int((sig & (df["log2FC"] < 0)).sum())

    # Pathway coloring (TRANSPARENCY): color the FULL membership of each story
    # pathway present in the data -- up and down, significant and not -- rather
    # than only the enriched-direction hits. The FDR line and volcano position
    # then show that the pathway skews up (more up than down), instead of us
    # pre-filtering to the up genes, which could read as cherry-picking.
    curated = INTERESTING_PATHWAYS[comp]["up"] + INTERESTING_PATHWAYS[comp]["down"]
    members, source = get_pathway_members(comp, curated)
    in_data = set(df["gene"])
    gene_color, gene_pathway, present_pathways = {}, {}, []
    for term in curated:
        genes = sorted(members.get(term, set()) & in_data)
        if not genes:
            continue
        present_pathways.append(term)
        for g in genes:
            gene_color.setdefault(g, PATHWAY_COLOR.get(term, "#444444"))
            gene_pathway.setdefault(g, term)
    df["pcolor"] = df["gene"].map(gene_color)   # NOT nulled by significance

    # background layers (two-tone gray by FDR<0.05)
    ns = df[~sig]
    ax.scatter(ns["log2FC"], ns["nlp"], s=10, c=GRAY_NS, alpha=0.35,
               linewidths=0, zorder=1)
    sig_plain = df[sig & df["pcolor"].isna()]
    ax.scatter(sig_plain["log2FC"], sig_plain["nlp"], s=16, c=GRAY_SIG,
               alpha=0.55, linewidths=0, zorder=2)

    # interactors as gold dots (below pathway genes)
    inter_df = df[df["gene"].isin(interactors)]
    if CIRCLE_ONLY_SIG_INTERACTORS:
        inter_df = inter_df[inter_df["fdr"] < FDR_THRESHOLD]
    ax.scatter(inter_df["log2FC"], inter_df["nlp"], s=95, c=INTERACTOR_COLOR,
               edgecolors="#8a6d00", linewidths=0.5, alpha=0.95, zorder=3)

    # story-pathway genes (FULL membership; drawn on top of interactors)
    path_genes = df[df["pcolor"].notna()]
    ax.scatter(path_genes["log2FC"], path_genes["nlp"], s=54,
               c=path_genes["pcolor"], alpha=0.95, linewidths=0.3,
               edgecolors="white", zorder=4)

    # provenance: exactly which genes got colored, with direction + significance
    prov = path_genes.copy()
    prov["pathway"] = prov["gene"].map(gene_pathway)
    prov["significant"] = prov["fdr"] < FDR_THRESHOLD
    prov = prov[["pathway", "gene", "log2FC", "p_value", "fdr", "significant"]] \
        .sort_values(["pathway", "fdr"])
    prov.to_csv(os.path.join(OUTPUT_DIR, f"Supp_pathway_membership_{comp}.tsv"),
                sep="\t", index=False)
    up_m = int((prov["log2FC"] > 0).sum()); dn_m = int((prov["log2FC"] < 0).sum())
    log(f"  [{comp}] pathway members colored (source={source}): "
        f"up={up_m} down={dn_m} sig={int(prov['significant'].sum())}")

    # A3 enzymes as dark-green diamonds
    a3_df = df[df["gene"].isin(A3_ENZYMES)]
    ax.scatter(a3_df["log2FC"], a3_df["nlp"], s=300, marker="D",
               facecolors=A3_COLOR, edgecolors="black", linewidths=2.0, zorder=7)

    # threshold lines
    if sig.any():
        yth = -np.log10(max(df.loc[sig, "p_value"].max(), 1e-300))
        ax.axhline(yth, ls="--", c="#555555", lw=1.5, zorder=0.5)
        ax.text(0.5, yth, "FDR < 0.05", transform=ax.get_yaxis_transform(),
                ha="center", va="bottom", fontsize=FONT_TICK - 4, color="#555555",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1),
                zorder=8)
    ax.axvline(0, ls="--", c="#888888", lw=1.0, zorder=0.5)

    # labels
    label_rows, labeled = [], set()

    def add(row, color="black", weight="normal", size=FONT_GENE_LABEL):
        alias = A3_SYMBOL_TO_ALIAS.get(row["gene"], row["gene"])
        label_rows.append(dict(text=alias, px=float(row["log2FC"]),
                               py=float(row["nlp"]), color=color,
                               weight=weight, size=size))
        labeled.add(row["gene"])

    def _sorted(frame, spec, default_metric=RANK_PATHWAY_LABELS_BY):
        metric, direction = spec
        asc = (direction == "asc")
        if metric == "abs_log2FC":
            f = frame.assign(_k=frame["log2FC"].abs())
            return f.sort_values("_k", ascending=asc)
        if metric in frame.columns:
            return frame.sort_values(metric, ascending=asc)
        return frame.sort_values(default_metric if default_metric in frame.columns
                                 else "fdr", ascending=asc)

    # pathway gene labels: ranked PER pathway by its configured metric/direction,
    # with the label budget split evenly across the pathways present.
    path_genes = path_genes.assign(pathway=path_genes["gene"].map(gene_pathway))
    n_pw = max(1, len(present_pathways))
    per_pw = -(-N_PATHWAY_LABELS // n_pw)   # ceil division
    for term in present_pathways:
        spec = PATHWAY_LABEL_RANK.get(term, (RANK_PATHWAY_LABELS_BY, "desc"))
        sub = _sorted(path_genes[path_genes["pathway"] == term], spec)
        for _, row in sub.head(per_pw).iterrows():
            if row["gene"] not in labeled:
                add(row, color=row["pcolor"])

    # interactor labels: ranked per comparison
    ispec = INTERACTOR_LABEL_RANK.get(comp, ("fdr", "asc"))
    ranked = _sorted(inter_df, ispec)
    ranked = ranked if LABEL_ALL_INTERACTORS else ranked.head(N_INTERACTOR_LABELS)
    for _, row in ranked.iterrows():
        if row["gene"] not in labeled:
            add(row, color=INTERACTOR_LABEL_COLOR)
    for _, row in a3_df.iterrows():
        if row["gene"] not in labeled:
            add(row, color=A3_COLOR, weight="bold", size=FONT_GENE_LABEL + 4)

    place_labels(ax, label_rows)

    ax.set_xlabel("log2 fold change", fontsize=FONT_AXIS)
    ax.set_ylabel("-log10(p)", fontsize=FONT_AXIS)
    ax.set_title(f"{LABELS[comp]}   (up {n_up} / down {n_down}, FDR<0.05)",
                 fontsize=FONT_VOLC_TITLE, color=hue, fontweight="bold", pad=14)
    ax.tick_params(labelsize=FONT_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # legend on the side (right of each panel)
    handles = [Line2D([0], [0], marker="o", color="none",
                      markerfacecolor=PATHWAY_COLOR.get(p, "#444444"),
                      markersize=16, label=p) for p in present_pathways]
    handles += [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=INTERACTOR_COLOR, markeredgecolor="#8a6d00",
               markersize=16, label="A3 interactor"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=A3_COLOR,
               markeredgecolor="black", markersize=16, label="A3 enzyme"),
    ]
    ax.legend(handles=handles, fontsize=FONT_LEGEND, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), frameon=False, handletextpad=0.4,
              borderaxespad=0.0)

    log(f"  [{comp}] up={n_up} down={n_down}  colored(sig)={int(path_genes.shape[0])}  "
        f"interactors={int(inter_df.shape[0])}  labels={len(label_rows)}")


# =============================================================================
# SUPPLEMENT DEG TABLES (csv + tsv)
# =============================================================================
def export_deg_tables(comp):
    """Write the DEG list for one comparison as CSV and TSV:
    columns gene, log2FC, p_value, p_adjusted. Significant only when
    EXPORT_SIG_ONLY, sorted most-significant first."""
    df = pd.read_csv(deg_path(comp))
    out = df[["gene", "log2FC", "p_value", "fdr"]].copy()
    out = out.rename(columns={"fdr": "p_adjusted"})
    if EXPORT_SIG_ONLY:
        out = out[out["p_adjusted"] < FDR_THRESHOLD]
    out = out.sort_values("p_adjusted").reset_index(drop=True)

    stem = os.path.join(OUTPUT_DIR, f"Supp_Table_DEG_{comp}")
    out.to_csv(stem + ".csv", index=False)
    out.to_csv(stem + ".tsv", sep="\t", index=False)
    scope = f"FDR<{FDR_THRESHOLD}" if EXPORT_SIG_ONLY else "all tested"
    log(f"  [{comp}] DEG table rows={len(out)} ({scope}) -> "
        f"{os.path.relpath(stem + '.csv', BASE_DIR)} (+ .tsv)")


# =============================================================================
# MAIN
# =============================================================================
def main():
    ensure_dir(OUTPUT_DIR)
    log("=" * 80)
    log("  DEG + KEGG SUPPLEMENT (v7)")
    log("=" * 80)

    interactors = set(load_harris_interactors()) - set(A3_ENZYMES)
    log(f"  Harris interactors loaded: {len(interactors)}")
    if not interactors:
        log("  [WARN] no interactors loaded; interactor dots skipped.")

    log("")
    log("  DEG supplement tables:")
    for comp in COMPARISONS:
        export_deg_tables(comp)

    fig = plt.figure(figsize=FIGSIZE)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.5],
                  hspace=0.30, wspace=0.48)

    for j, comp in enumerate(COMPARISONS):
        draw_kegg_bar(fig.add_subplot(gs[0, j]), comp)
    for j, comp in enumerate(COMPARISONS):
        draw_volcano(fig.add_subplot(gs[1, j]), comp, interactors)

    fig.suptitle("DEG programs and pathway enrichment across tumor populations",
                 fontsize=FONT_TITLE + 2, fontweight="bold", y=0.995)

    pdf = os.path.join(OUTPUT_DIR, OUTFILE_STEM + ".pdf")
    png = os.path.join(OUTPUT_DIR, OUTFILE_STEM + ".png")
    fig.savefig(pdf, dpi=DPI, bbox_inches="tight")
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    log(f"\n  [SAVE] {os.path.relpath(pdf, BASE_DIR)}")
    log(f"  [SAVE] {os.path.relpath(png, BASE_DIR)}")
    log("  DONE")


if __name__ == "__main__":
    main()
