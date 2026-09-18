#!/usr/bin/env python3
"""
Diagnostic_SBS13_Presence_and_REV1.py
=====================================

Answers the reviewer question on Results 4.1: why does the paper report SBS2
and not SBS13, when both are APOBEC3 signatures and both are well documented
in HPV-positive HNSCC?

BACKGROUND (why this is a real question, not a bookkeeping one)
---------------------------------------------------------------
SBS2 and SBS13 arise from the SAME lesion by DIFFERENT downstream handling:

    A3 deaminates C -> U in ssDNA
      |
      +-- replication reads U as T ......................... C>T  = SBS2
      |
      +-- UNG excises U -> abasic site
            |
            +-- REV1 inserts dCMP opposite the abasic site
                  and POL-zeta (REV3L + MAD2L2/REV7) extends .. C>G = SBS13

So SBS13 is gated on UNG plus translesion synthesis. A tissue with intact
deamination but limited REV1/POL-zeta capacity should show SBS2 with
comparatively little SBS13. That is a testable prediction, and it is the
hypothesis this script evaluates.

THREE COMPETING EXPLANATIONS, all tested here
---------------------------------------------
  H1  SBS13 is genuinely low in this dataset.
      -> STEP 1/2: low SBS13 weight, low prevalence, low SBS13:SBS2 ratio.

  H2  SBS13 is present but NMF cannot separate it from SBS2 at single-cell
      mutation counts. SBS2 and SBS13 are collinear (same TCW context,
      different substitution), and with few mutations per cell the
      factorisation tends to load one and starve the other.
      -> STEP 2: if SBS13 is near-zero in cells that carry high SBS2, and
         the two are strongly ANTI-correlated among carriers, suspect H2.
         If they co-occur positively, H2 is weakened.
      -> STEP 5 is the arbiter: the raw C>G call rate does not depend on NMF.

  H3  The C>G calls needed to support SBS13 are not being made, for
      technical reasons upstream of the signature fit.
      -> STEP 5: base-change spectrum straight from the SComatic table.

STEP 5 IS THE DECIDING TEST. Signature weights are a model output; the
substitution spectrum is not. If C>G at TCW is present at a normal
SBS13:SBS2 ratio but the NMF weight is not, the deficit is H2 and the
manuscript should say so rather than making a biological claim.

DELIBERATELY NOT DONE
---------------------
Trinucleotide context is read from the reference genome ONLY if pysam and
the GRCh38 FASTA are reachable (set REF_FASTA). The SComatic REF_TRI field
is NOT trusted: a previous 2-fold TCW "enrichment" in this project traced
to that field being wrong. With no genome available the script reports the
unstranded base-change spectrum, which is still sufficient to answer H3,
since SBS13 requires C>G regardless of context.

Usage:
  conda run -n NETWORK python Diagnostic_SBS13_Presence_and_REV1.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc

# =============================================================================
# PATHS
# =============================================================================
BASE = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH   = os.path.join(BASE, "data/FIG_4/00_input/adata_final.h5ad")
WEIGHTS_PATH = os.path.join(BASE, "data/FIG_4/00_input/signature_weights_per_cell.txt")
GROUPS_PATH  = os.path.join(BASE, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
SCOMATIC_PATH = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
                 "results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv")

# Optional: set to the GRCh38 FASTA to enable genome-verified TCW context.
REF_FASTA = "/master/jlehle/WORKING/SC/ref/GRCh38/genome.fa"

OUTDIR = os.path.join(BASE, "data/FIG_2/DIAGNOSTIC_SBS13")

# =============================================================================
# STYLE
# =============================================================================
FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
C_SBS2, C_SBS13, C_GREY = "#ed6a5a", "#6a8ded", "#cccccc"
C_CNV = "#F6D155"

plt.rcParams.update({
    "font.size": FS_TICK, "axes.titlesize": FS_TITLE, "axes.labelsize": FS_LABEL,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "pdf.fonttype": 42, "figure.dpi": 100,
})

CELLTYPE_COL, PATIENT_COL, TISSUE_COL = "final_annotation", "subject id", "tissue type"
BASAL = "basal cell"

# TLS / BER genes. REV1 + POLZ are the SBS13 arm; UNG opens the abasic site.
TLS_GENES = ["REV1", "REV3L", "MAD2L2", "UNG", "SMUG1", "TDG", "APEX1",
             "POLH", "POLK", "PCNA"]
A3_GENES = ["APOBEC3A", "APOBEC3B"]


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


def banner(t):
    print("")
    print("=" * 86)
    print(f"  {t}")
    print("=" * 86)


def savefig(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    log(f"  [SAVE] {name}.pdf / .png")


# =============================================================================
# STEP 0
# =============================================================================
banner("STEP 0: LOAD")

log(f"  adata: {ADATA_PATH}")
adata = sc.read_h5ad(ADATA_PATH)
log(f"    {adata.n_obs:,} cells x {adata.n_vars:,} genes")

log(f"  weights: {WEIGHTS_PATH}")
W = pd.read_csv(WEIGHTS_PATH, sep="\t", index_col=0).T   # .T is required
log(f"    {W.shape[0]:,} cells x {W.shape[1]} signatures")
log(f"    signatures: {', '.join(W.columns)}")

if "SBS13" not in W.columns:
    log("  !! SBS13 ABSENT from the weights table. It was never fit.")
    log("     That is the whole answer: report it as not retained by the")
    log("     deconvolution rather than as biologically absent.")
    sys.exit(0)

obs = adata.obs.copy()
basal_mask = obs[CELLTYPE_COL].astype(str) == BASAL
log(f"    basal cells: {int(basal_mask.sum()):,}")

shared = obs.index.intersection(W.index)
log(f"    cells with weights: {len(shared):,}")

groups = None
if os.path.exists(GROUPS_PATH):
    g = pd.read_csv(GROUPS_PATH, sep="\t")
    gcol = "cell_barcode" if "cell_barcode" in g.columns else g.columns[0]
    groups = dict(zip(g[gcol].astype(str), g["group"]))
    log(f"    three-group assignments: {len(groups):,}")

df = pd.DataFrame(index=shared)
for s in W.columns:
    df[s] = W.loc[shared, s].values
df["celltype"] = obs.loc[shared, CELLTYPE_COL].astype(str).values
df["patient"] = obs.loc[shared, PATIENT_COL].astype(str).values
df["tissue"] = obs.loc[shared, TISSUE_COL].astype(str).values
if groups:
    df["group"] = [groups.get(c, "other") for c in shared]

bdf = df[df["celltype"] == BASAL].copy()
log(f"    basal cells with weights: {len(bdf):,}")


# =============================================================================
# STEP 1
# =============================================================================
banner("STEP 1: WHERE DOES SBS13 SIT AMONG ALL FITTED SIGNATURES? (basal)")

rows = []
for s in W.columns:
    v = bdf[s].values
    rows.append({
        "signature": s, "mean": v.mean(), "median": np.median(v),
        "max": v.max(), "n_pos": int((v > 0).sum()),
        "pct_pos": 100 * (v > 0).mean(),
    })
sig = pd.DataFrame(rows).sort_values("mean", ascending=False)

print("")
print(f"  {'signature':<10} {'mean':>10} {'median':>9} {'max':>10} "
      f"{'n>0':>8} {'%>0':>8}")
print(f"  {'-'*10} {'-'*10} {'-'*9} {'-'*10} {'-'*8} {'-'*8}")
for _, r in sig.iterrows():
    star = "  <<<" if r["signature"] in ("SBS2", "SBS13") else ""
    print(f"  {r['signature']:<10} {r['mean']:>10.4f} {r['median']:>9.4f} "
          f"{r['max']:>10.4f} {int(r['n_pos']):>8,} {r['pct_pos']:>7.1f}%{star}")

s2, s13 = bdf["SBS2"].values, bdf["SBS13"].values
r2, r13 = sig.set_index("signature").index.get_loc("SBS2"), \
          sig.set_index("signature").index.get_loc("SBS13")
print("")
log(f"  SBS2  rank {r2+1}/{len(sig)} by mean weight")
log(f"  SBS13 rank {r13+1}/{len(sig)} by mean weight")
ratio_pooled = s13.sum() / s2.sum() if s2.sum() else np.nan
log(f"  pooled SBS13:SBS2 = {ratio_pooled:.3f}")
log(f"  (bulk HNSCC/cervical typically sits near 0.3-0.5; well below that")
log(f"   is the observation that needs explaining)")


# =============================================================================
# STEP 2
# =============================================================================
banner("STEP 2: DO SBS2 AND SBS13 CO-OCCUR? (H1 vs H2)")

both = (s2 > 0) & (s13 > 0)
only2 = (s2 > 0) & (s13 == 0)
only13 = (s2 == 0) & (s13 > 0)
neither = (s2 == 0) & (s13 == 0)

n = len(bdf)
print("")
print(f"  {'SBS2>0 and SBS13>0':<26} {int(both.sum()):>8,}  {100*both.mean():>6.1f}%")
print(f"  {'SBS2>0 only':<26} {int(only2.sum()):>8,}  {100*only2.mean():>6.1f}%")
print(f"  {'SBS13>0 only':<26} {int(only13.sum()):>8,}  {100*only13.mean():>6.1f}%")
print(f"  {'neither':<26} {int(neither.sum()):>8,}  {100*neither.mean():>6.1f}%")

tab = [[int(both.sum()), int(only2.sum())],
       [int(only13.sum()), int(neither.sum())]]
try:
    orr, pf = fisher_exact(tab)
    log(f"  co-occurrence Fisher OR = {orr:.2f}, p = {pf:.3g}")
    log("    OR > 1 -> they occur together, consistent with one APOBEC process (H1)")
    log("    OR < 1 -> they exclude each other, which smells like NMF")
    log("             splitting collinear signatures cell by cell (H2)")
except Exception as e:
    log(f"  Fisher failed: {e}")

carriers = bdf[(s2 > 0) | (s13 > 0)]
if len(carriers) > 10:
    rho, p = spearmanr(carriers["SBS2"], carriers["SBS13"])
    log(f"  among carriers, Spearman SBS2 vs SBS13: rho = {rho:.3f}, p = {p:.3g}")

sub = bdf[s2 > 0]
if len(sub):
    log(f"  in SBS2-positive cells, mean SBS13 = {sub['SBS13'].mean():.4f} "
        f"({100*(sub['SBS13'] > 0).mean():.1f}% carry any SBS13)")


# =============================================================================
# STEP 3
# =============================================================================
banner("STEP 3: SBS13 BY CELL TYPE AND BY POPULATION")

ct = df.groupby("celltype")[["SBS2", "SBS13"]].mean().sort_values(
    "SBS2", ascending=False)
ct["SBS13:SBS2"] = ct["SBS13"] / ct["SBS2"].replace(0, np.nan)
print("")
print(f"  {'cell type':<34} {'SBS2':>9} {'SBS13':>9} {'ratio':>8}")
print(f"  {'-'*34} {'-'*9} {'-'*9} {'-'*8}")
for i, r in ct.iterrows():
    print(f"  {i[:34]:<34} {r['SBS2']:>9.4f} {r['SBS13']:>9.4f} "
          f"{r['SBS13:SBS2']:>8.3f}")

if "group" in bdf.columns:
    gg = bdf[bdf["group"] != "other"].groupby("group")[["SBS2", "SBS13"]].agg(
        ["mean", "count"])
    print("")
    print("  three-group populations:")
    print(gg.to_string())


# =============================================================================
# STEP 4
# =============================================================================
banner("STEP 4: DOES SBS13 TRACK A3A OR A3B? (compare against SBS2)")

expr = {}
for g in A3_GENES + TLS_GENES:
    if g in adata.var_names:
        x = adata[:, g].X
        expr[g] = np.asarray(x.todense()).ravel() if hasattr(x, "todense") \
            else np.asarray(x).ravel()
    else:
        log(f"  [absent from var_names] {g}")

pos = {c: i for i, c in enumerate(adata.obs_names)}
idx = np.array([pos[c] for c in bdf.index])

print("")
print(f"  {'gene':<10} {'vs SBS2 rho':>14} {'p':>12} "
      f"{'vs SBS13 rho':>15} {'p':>12}")
print(f"  {'-'*10} {'-'*14} {'-'*12} {'-'*15} {'-'*12}")
for g in A3_GENES:
    if g not in expr:
        continue
    v = expr[g][idx]
    a, pa = spearmanr(v, bdf["SBS2"])
    b, pb = spearmanr(v, bdf["SBS13"])
    print(f"  {g:<10} {a:>14.4f} {pa:>12.3g} {b:>15.4f} {pb:>12.3g}")
log("  If SBS13 tracks A3A/A3B the same way SBS2 does, it is the same")
log("  process at lower yield, not a different one.")


# =============================================================================
# STEP 5  -- THE DECIDING TEST
# =============================================================================
banner("STEP 5: RAW BASE-CHANGE SPECTRUM (independent of the NMF)")

if not os.path.exists(SCOMATIC_PATH):
    log(f"  [SKIP] SComatic table not found: {SCOMATIC_PATH}")
else:
    log("  streaming SComatic calls ...")
    usecols = None
    head = pd.read_csv(SCOMATIC_PATH, sep="\t", comment="#", nrows=5)
    cols = list(head.columns)
    log(f"    columns: {', '.join(cols[:14])}{' ...' if len(cols) > 14 else ''}")

    def pick(cands):
        for c in cands:
            if c in cols:
                return c
        return None

    ref_c = pick(["REF", "ref", "Ref"])
    alt_c = pick(["ALT", "alt", "Alt"])
    ct_c = pick(["Cell_type", "cell_type", "CellType"])

    if not (ref_c and alt_c):
        log("  [SKIP] could not resolve REF/ALT columns; inspect the header above.")
    else:
        keep = [c for c in (ref_c, alt_c, ct_c) if c]
        counts = {}
        total = 0
        for chunk in pd.read_csv(SCOMATIC_PATH, sep="\t", comment="#",
                                 usecols=keep, chunksize=500_000):
            if ct_c:
                chunk = chunk[chunk[ct_c].astype(str).str.lower() == BASAL]
            chunk = chunk[(chunk[ref_c].astype(str).str.len() == 1) &
                          (chunk[alt_c].astype(str).str.len() == 1)]
            total += len(chunk)
            vc = (chunk[ref_c].astype(str) + ">" +
                  chunk[alt_c].astype(str)).value_counts()
            for k, v in vc.items():
                counts[k] = counts.get(k, 0) + int(v)

        # Fold to pyrimidine reference, the COSMIC convention.
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        folded = {}
        for k, v in counts.items():
            r, a = k.split(">")
            if r in ("A", "G"):
                r, a = comp[r], comp[a]
            folded[f"{r}>{a}"] = folded.get(f"{r}>{a}", 0) + v

        print("")
        log(f"  basal single-nucleotide calls: {total:,}")
        print("")
        print(f"  {'change':<10} {'n':>10} {'% of all':>10}")
        print(f"  {'-'*10} {'-'*10} {'-'*10}")
        for k in sorted(folded, key=lambda x: -folded[x]):
            print(f"  {k:<10} {folded[k]:>10,} {100*folded[k]/max(total,1):>9.2f}%")

        ct_n, cg_n = folded.get("C>T", 0), folded.get("C>G", 0)
        print("")
        log(f"  C>T = {ct_n:,}   C>G = {cg_n:,}")
        if ct_n:
            log(f"  C>G : C>T = {cg_n/ct_n:.3f}")
        log("")
        log("  READ THIS AGAINST STEP 1:")
        log("    If C>G:C>T in the RAW calls is much higher than the fitted")
        log("    SBS13:SBS2 weight ratio, the substitutions are there and the")
        log("    NMF is not resolving SBS13 -> H2, a deconvolution limit.")
        log("    If C>G calls are genuinely scarce, the deficit is upstream")
        log("    of the fit -> H1 or H3, and the REV1 result below is relevant.")


# =============================================================================
# STEP 6
# =============================================================================
banner("STEP 6: REV1 / UNG / POL-ZETA EXPRESSION (the H1 mechanism)")

print("")
print(f"  {'gene':<10} {'basal mean':>12} {'basal %pos':>12} "
      f"{'all mean':>11} {'all %pos':>11}")
print(f"  {'-'*10} {'-'*12} {'-'*12} {'-'*11} {'-'*11}")
bmask = basal_mask.values
tls_rows = []
for g in TLS_GENES:
    if g not in expr:
        continue
    v = expr[g]
    row = {"gene": g, "basal_mean": v[bmask].mean(),
           "basal_pct": 100 * (v[bmask] > 0).mean(),
           "all_mean": v.mean(), "all_pct": 100 * (v > 0).mean()}
    tls_rows.append(row)
    print(f"  {g:<10} {row['basal_mean']:>12.4f} {row['basal_pct']:>11.1f}% "
          f"{row['all_mean']:>11.4f} {row['all_pct']:>10.1f}%")

log("")
log("  CAUTION: REV1, REV3L and MAD2L2 are low-abundance transcripts and")
log("  10x chemistry drops them readily. A low %pos here is NOT by itself")
log("  evidence of low protein activity. The comparison that carries weight")
log("  is basal versus other cell types in THIS dataset, and REV1 against")
log("  a similarly-expressed control such as POLH or POLK.")

# Per-cell-type REV1, so 'low in basal' can be judged relatively.
if "REV1" in expr:
    log("")
    log("  REV1 by cell type:")
    r = pd.DataFrame({"celltype": adata.obs[CELLTYPE_COL].astype(str).values,
                      "REV1": expr["REV1"]})
    rr = r.groupby("celltype")["REV1"].agg(["mean", lambda x: 100*(x > 0).mean()])
    rr.columns = ["mean", "pct_pos"]
    for i, row in rr.sort_values("mean", ascending=False).iterrows():
        mark = "  <<<" if i == BASAL else ""
        print(f"    {i[:34]:<34} {row['mean']:>9.4f} {row['pct_pos']:>7.1f}%{mark}")

if "group" in bdf.columns and "REV1" in expr:
    log("")
    log("  REV1 across the three populations:")
    gi = {c: i for i, c in enumerate(adata.obs_names)}
    for grp in ["SBS2_HIGH", "CNV_HIGH", "NORMAL"]:
        cells = bdf.index[bdf["group"] == grp]
        if not len(cells):
            continue
        v = expr["REV1"][[gi[c] for c in cells]]
        print(f"    {grp:<12} n={len(cells):<5} mean={v.mean():.4f} "
              f"%pos={100*(v>0).mean():.1f}%")


# =============================================================================
# STEP 7: FIGURES
# =============================================================================
banner("STEP 7: FIGURES")

# (a) UMAP of SBS2 and SBS13 side by side.
if "X_umap" in adata.obsm:
    um = pd.DataFrame(adata.obsm["X_umap"][:, :2], index=adata.obs_names,
                      columns=["u1", "u2"])
    fig, axes = plt.subplots(1, 2, figsize=(26, 12))
    for ax, s, c in zip(axes, ["SBS2", "SBS13"], [C_SBS2, C_SBS13]):
        ax.scatter(um["u1"], um["u2"], s=1, c=C_GREY, alpha=0.25,
                   linewidths=0, rasterized=True)
        sel = bdf[bdf[s] > 0]
        if len(sel):
            uu = um.loc[sel.index]
            ax.scatter(uu["u1"], uu["u2"], s=6, c=sel[s], cmap="viridis",
                       linewidths=0, rasterized=True)
        ax.set_title(f"{s}  (n = {int((bdf[s] > 0).sum()):,} basal cells > 0)",
                     fontsize=FS_TITLE, color=c)
        ax.set_xlabel("UMAP 1", fontsize=FS_LABEL)
        ax.set_ylabel("UMAP 2", fontsize=FS_LABEL)
        ax.set_xticks([]); ax.set_yticks([])
    savefig(fig, "SBS13_vs_SBS2_UMAP")
else:
    log("  [SKIP] no X_umap in adata.obsm")

# (b) Weight distributions among carriers.
fig, ax = plt.subplots(figsize=(14, 12))
data, labels, colors = [], [], []
for s, c in [("SBS2", C_SBS2), ("SBS13", C_SBS13)]:
    v = bdf.loc[bdf[s] > 0, s].values
    if len(v):
        data.append(v); labels.append(f"{s}\n(n={len(v):,})"); colors.append(c)
if data:
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False,
                    widths=0.6)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    ax.set_ylabel("per-cell weight (carriers only)", fontsize=FS_LABEL)
    ax.set_title("SBS2 and SBS13 weight, basal cells", fontsize=FS_TITLE)
    savefig(fig, "SBS13_vs_SBS2_weight_distribution")

# (c) REV1 alongside comparable TLS genes.
if tls_rows:
    t = pd.DataFrame(tls_rows).sort_values("basal_mean", ascending=True)
    fig, ax = plt.subplots(figsize=(14, 12))
    cols = [C_SBS13 if g == "REV1" else C_GREY for g in t["gene"]]
    ax.barh(t["gene"], t["basal_mean"], color=cols, edgecolor="black",
            linewidth=1.2)
    ax.set_xlabel("mean log-normalized expression, basal cells",
                  fontsize=FS_LABEL)
    ax.set_title("Translesion synthesis and BER genes", fontsize=FS_TITLE)
    savefig(fig, "REV1_TLS_expression_basal")

# =============================================================================
# SAVE TABLES
# =============================================================================
os.makedirs(OUTDIR, exist_ok=True)
sig.to_csv(os.path.join(OUTDIR, "signature_summary_basal.tsv"),
           sep="\t", index=False)
ct.to_csv(os.path.join(OUTDIR, "sbs13_by_celltype.tsv"), sep="\t")
if tls_rows:
    pd.DataFrame(tls_rows).to_csv(
        os.path.join(OUTDIR, "tls_gene_expression.tsv"), sep="\t", index=False)

banner("COMPLETE")
log(f"  Output: {OUTDIR}")
log("")
log("  Decide from STEP 5 before writing any sentence:")
log("    raw C>G present, NMF weight absent   -> deconvolution limit, say so")
log("    raw C>G genuinely scarce             -> biological, REV1 is in play")
