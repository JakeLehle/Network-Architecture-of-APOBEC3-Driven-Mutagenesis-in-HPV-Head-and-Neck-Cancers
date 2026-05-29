#!/usr/bin/env python3
"""
Diagnostic_Regenerate_Panel_A_C.py
====================================

Standalone script to regenerate two supplemental panels:

  Panel A: Patient distribution of SBS2-HIGH cells
    1. High-contributor patients get bold coral labels
    2. Fold enrichment (obs/exp) shown alongside raw counts

  Panel C: PCA + correlation heatmap of SBS2-HIGH cells
    1. Point sizes 2x larger (HC: 400, non-HC: 160)
    2. Everything else identical

Usage:
    cd scripts/PATIENT_SPECIFIC_EFFECTS/
    conda run -n NETWORK python Diagnostic_Regenerate_Panel_A_C.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from patient_config import (
    DIR_00_DIAG, FIGURE_5_PANELS, DE_GENES_DIR,
    PATIENT_COL, CELLTYPE_COL,
    HIGH_CONTRIBUTORS, RANDOM_SEED,
    FONT_TITLE, FONT_LABEL, FONT_TICK,
    COLOR_SBS2_HIGH, COLOR_OTHER,
    banner, log, ensure_dir, load_adata, load_three_groups, load_groups,
)

# =============================================================================
# STYLE
# =============================================================================
DPI = 300
FONT_ANNOT = FONT_TICK - 2    # 26

COLOR_HC_LABEL = COLOR_SBS2_HIGH   # "#ed6a5a"
COLOR_HC_BAR   = COLOR_SBS2_HIGH   # coral for HC
COLOR_NON_HC   = "#f4a582"         # lighter salmon for non-HC HIGH cells
COLOR_BASAL    = "#d4d4d4"         # gray for other basal

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.size': FONT_TICK,
    'axes.titlesize': FONT_TITLE,
    'axes.labelsize': FONT_LABEL,
    'xtick.labelsize': FONT_TICK,
    'ytick.labelsize': FONT_TICK,
})

def short(p):
    return p.replace('Patient ', '')

# =============================================================================
# LOAD DATA
# =============================================================================
banner("PANEL A REGENERATION: PATIENT DISTRIBUTION + FOLD ENRICHMENT")

out_dir = ensure_dir(FIGURE_5_PANELS)

# Try loading precomputed enrichment table from the diagnostic
enrich_path = os.path.join(DIR_00_DIAG, "patient_enrichment_SBS2_HIGH_v3.tsv")
if os.path.exists(enrich_path):
    log(f"  Loading enrichment table: {enrich_path}")
    enrich_df = pd.read_csv(enrich_path, sep='\t')
    log(f"  {len(enrich_df)} patients loaded")
else:
    log(f"  v3 enrichment table not found, computing from scratch...")
    adata = load_adata()
    sbs2_high, cnv_high, normal = load_three_groups()

    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
    basal.obs['group'] = 'other'
    basal.obs.loc[basal.obs_names.isin(sbs2_high), 'group'] = 'SBS2_HIGH'

    patients = sorted(basal.obs[PATIENT_COL].unique())
    n_sbs2 = (basal.obs['group'] == 'SBS2_HIGH').sum()
    rows = []
    for patient in patients:
        pmask = basal.obs[PATIENT_COL] == patient
        n_basal = pmask.sum()
        n_high = (pmask & (basal.obs['group'] == 'SBS2_HIGH')).sum()
        pct = 100.0 * n_high / n_sbs2 if n_sbs2 > 0 else 0
        expected_frac = n_basal / basal.n_obs
        observed_frac = n_high / n_sbs2 if n_sbs2 > 0 else 0
        fold = observed_frac / expected_frac if expected_frac > 0 else 0
        rows.append({
            'patient': patient, 'n_basal': n_basal,
            'n_sbs2_high': n_high, 'pct_sbs2_high': pct,
            'fold_sbs2_high': fold,
        })
    enrich_df = pd.DataFrame(rows)

# Sort by SBS2-HIGH count (ascending for horizontal bars, top = highest)
enrich_df = enrich_df.sort_values('n_sbs2_high', ascending=True).reset_index(drop=True)

# =============================================================================
# PLOT
# =============================================================================
banner("PLOTTING")

n_patients = len(enrich_df)
fig, ax = plt.subplots(figsize=(16, 10))

y = np.arange(n_patients)
n_other = enrich_df['n_basal'].values - enrich_df['n_sbs2_high'].values
n_high = enrich_df['n_sbs2_high'].values
patients = enrich_df['patient'].values
folds = enrich_df['fold_sbs2_high'].values

# Determine bar color per patient: HC patients get strong coral, others get lighter
bar_colors = []
for p in patients:
    if p in HIGH_CONTRIBUTORS:
        bar_colors.append(COLOR_HC_BAR)
    else:
        bar_colors.append(COLOR_NON_HC)

# Draw other basal (gray)
ax.barh(y, n_other, color=COLOR_BASAL, edgecolor='black',
        linewidth=0.5, label='Other basal', zorder=2)

# Draw SBS2-HIGH bars individually for per-patient coloring
for i in range(n_patients):
    ax.barh(y[i], n_high[i], left=n_other[i], color=bar_colors[i],
            edgecolor='black', linewidth=0.5, zorder=2)

# Annotations: count + fold enrichment
x_max = enrich_df['n_basal'].max()
for i in range(n_patients):
    nh = n_high[i]
    fold = folds[i]
    total = enrich_df['n_basal'].values[i]

    if nh > 0:
        # Show count and fold
        fold_str = f'{fold:.1f}x'
        if fold >= 1.5:
            fold_color = COLOR_HC_LABEL
        elif fold >= 1.0:
            fold_color = '#333333'
        else:
            fold_color = '#888888'

        ax.text(total + x_max * 0.01, y[i],
                f'{nh}  ({fold_str})',
                va='center', ha='left',
                fontsize=FONT_ANNOT,
                fontweight='bold' if patients[i] in HIGH_CONTRIBUTORS else 'normal',
                color=fold_color)

# Y-axis labels: bold + coral for HC patients
ax.set_yticks(y)
labels = [short(p) for p in patients]
ax.set_yticklabels(labels)

for i, p in enumerate(patients):
    if p in HIGH_CONTRIBUTORS:
        ax.get_yticklabels()[i].set_color(COLOR_HC_LABEL)
        ax.get_yticklabels()[i].set_fontweight('bold')

ax.set_xlabel('Number of Basal Cells', fontsize=FONT_LABEL)
ax.set_title('Patient Distribution of SBS2-HIGH Basal Cells',
             fontsize=FONT_TITLE, pad=15)

# Legend
legend_handles = [
    Patch(facecolor=COLOR_BASAL, edgecolor='black', linewidth=0.5,
          label='Other basal'),
    Patch(facecolor=COLOR_HC_BAR, edgecolor='black', linewidth=0.5,
          label='SBS2-HIGH (enriched, fold > 2x)'),
    Patch(facecolor=COLOR_NON_HC, edgecolor='black', linewidth=0.5,
          label='SBS2-HIGH (not enriched)'),
]
ax.legend(handles=legend_handles, loc='lower right',
          fontsize=FONT_ANNOT, framealpha=0.9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

for ext in ['pdf', 'png']:
    path = os.path.join(out_dir, f"Supp_Panel_A_Patient_Distribution.{ext}")
    fig.savefig(path, dpi=DPI, bbox_inches='tight')
    log(f"  [SAVE] {path}")
plt.close()

# Also print the table to confirm
banner("SUMMARY TABLE")
enrich_sorted = enrich_df.sort_values('n_sbs2_high', ascending=False)
log(f"  {'Patient':<20s}  {'N_HIGH':>8s}  {'Fold':>6s}  {'HC?':>5s}")
log(f"  {'-'*20}  {'-'*8}  {'-'*6}  {'-'*5}")
for _, r in enrich_sorted.iterrows():
    hc = "YES" if r['patient'] in HIGH_CONTRIBUTORS else ""
    log(f"  {str(r['patient']):<20s}  {int(r['n_sbs2_high']):>8d}  "
        f"{r['fold_sbs2_high']:>5.1f}x  {hc:>5s}")

banner("PANEL A REGENERATION COMPLETE")


# #############################################################################
#
#   PANEL C: PCA + HIERARCHICAL CLUSTERING HEATMAP (2x point sizes)
#
# #############################################################################
banner("PANEL C: PCA + CORRELATION HEATMAP (2x point sizes)")

# Load adata (needed for Panel C even if Panel A used the enrichment table)
log("  Loading adata for Panel C...")
adata = load_adata()
sbs2_high_cells, normal_cells = load_groups()

high_in = sbs2_high_cells & set(adata.obs_names)
adata_high = adata[list(high_in)].copy()
pl = adata_high.obs[PATIENT_COL].values
up = sorted(set(pl))
log(f"  SBS2-HIGH cells in adata: {adata_high.n_obs:,}")
for p in up:
    log(f"    {short(p)}: {(pl == p).sum()}")

# Load DE genes
de_path = os.path.join(DE_GENES_DIR, "SC_selected_genes_filtered.csv")
if os.path.exists(de_path):
    de_genes = [g for g in pd.read_csv(de_path)['gene']
                if g in adata_high.var_names]
    log(f"  DE genes matched: {len(de_genes)}")
else:
    log(f"  DE file not found at {de_path}, using top 5000 var_names")
    de_genes = list(adata_high.var_names[:5000])

# PCA
X = adata_high[:, de_genes].X
if scipy.sparse.issparse(X):
    X = X.toarray()
X_s = StandardScaler().fit_transform(X)
n_comp = min(50, X_s.shape[1], X_s.shape[0] - 1)
pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
X_pca = pca.fit_transform(X_s)
var = pca.explained_variance_ratio_
log(f"  PC1: {var[0]*100:.1f}%, PC2: {var[1]*100:.1f}%, "
    f"PC1-10: {sum(var[:10])*100:.1f}%")

FONT_LEGEND = FONT_TICK - 4   # 24
FONT_ANNOT_C = FONT_TICK - 2  # 26

cm = plt.cm.tab20
pcols = {p: cm(i / max(len(up) - 1, 1)) for i, p in enumerate(sorted(up))}

fig, axes = plt.subplots(1, 2, figsize=(28, 12))

# --- PCA scatter (2x point sizes) ---
ax = axes[0]
for p in up:
    m = pl == p
    mk = '*' if p in HIGH_CONTRIBUTORS else 'o'
    # Original: 200 HC, 80 non-HC. Now 2x: 400 HC, 160 non-HC
    sz = 600 if p in HIGH_CONTRIBUTORS else 240
    ax.scatter(X_pca[m, 0], X_pca[m, 1], c=[pcols[p]], s=sz, marker=mk,
               alpha=0.6, edgecolors='black', linewidth=0.5, label=short(p))
ax.set_xlabel(f'PC1 ({var[0] * 100:.1f}%)')
ax.set_ylabel(f'PC2 ({var[1] * 100:.1f}%)')
ax.set_title('SBS2-HIGH Cells by Patient')
ax.legend(fontsize=FONT_LEGEND - 4, ncol=2, loc='best', markerscale=1.5)

# --- Correlation heatmap (unchanged) ---
ax = axes[1]
pts_c = [p for p in up if (pl == p).sum() >= 2]
profs = {p: X[pl == p].mean(axis=0) for p in pts_c}
corr = pd.DataFrame(profs).corr(method='spearman')
dm = np.maximum(1 - corr.values, 0)
np.fill_diagonal(dm, 0)
lo = leaves_list(linkage(squareform(dm), method='average'))
op = [corr.columns[i] for i in lo]
co = corr.loc[op, op]
labs = [short(p) for p in op]

im = ax.imshow(co.values, cmap='RdYlBu_r', vmin=0.5, vmax=1)
ax.set_xticks(range(len(labs)))
ax.set_xticklabels(labs, rotation=45, ha='right')
ax.set_yticks(range(len(labs)))
ax.set_yticklabels(labs)
for i, p in enumerate(op):
    if p in HIGH_CONTRIBUTORS:
        ax.get_xticklabels()[i].set_color(COLOR_SBS2_HIGH)
        ax.get_xticklabels()[i].set_fontweight('bold')
        ax.get_yticklabels()[i].set_color(COLOR_SBS2_HIGH)
        ax.get_yticklabels()[i].set_fontweight('bold')
for i in range(len(labs)):
    for j in range(len(labs)):
        ax.text(j, i, f'{co.values[i, j]:.2f}', ha='center', va='center',
                fontsize=FONT_ANNOT_C - 4,
                color='white' if co.values[i, j] < 0.7 else 'black')
plt.colorbar(im, ax=ax, label='Spearman rho', shrink=0.7)
ax.set_title('Mean HIGH-Cell Profile Correlation\n(hierarchically clustered)')

plt.tight_layout()
for ext in ['pdf', 'png']:
    path = os.path.join(out_dir, f"Supp_Panel_C_PCA_Correlation.{ext}")
    fig.savefig(path, dpi=DPI, bbox_inches='tight')
    log(f"  [SAVE] {path}")
plt.close()

banner("PANELS A + C REGENERATION COMPLETE")
