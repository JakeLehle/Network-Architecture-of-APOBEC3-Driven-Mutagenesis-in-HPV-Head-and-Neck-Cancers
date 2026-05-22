#!/usr/bin/env python3
"""
Generate_Supplemental_Patient_Effects_v4.py
=============================================

Updated from v3 for V4 pipeline compatibility:
  - Panel D redesigned: 4-metric LOPO sensitivity (chain recovery,
    community ARI, gene overlap Jaccard, A3 wall integrity)
  - Removed Tier3B tumor-only references
  - Uses updated patient_config gene sets (activating chain, anchors)
  - Font/color constants from config

Panels:
  A: Patient distribution of SBS2-HIGH cells (sorted by HIGH count)
  B: Per-patient A3A/A3B expression
  C: PCA + hierarchical clustering of HIGH cells
  D: LOPO sensitivity (V4 metrics)
  E: Variant sharing spectrum (5-tier stacked bar)
  F: Genome ideogram (sharing tiers)

Usage:
  conda run -n NETWORK python Generate_Supplemental_Patient_Effects_v4.py
"""

import os, sys, numpy as np, pandas as pd, scanpy as sc, scipy.sparse
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch

from patient_config import (
    DIR_00_DIAG, DIR_03_SENSITIVITY, FIGURE_5_PANELS, DE_GENES_DIR,
    PATIENT_COL, CELLTYPE_COL, TISSUE_COL, SCOMATIC_PATH, WEIGHTS_PATH,
    A3_GENES_SYMBOLS, A3_INTERACTOR_ANCHORS, ACTIVATING_CHAIN_GENES,
    HIGH_CONTRIBUTORS, RANDOM_SEED, MIN_COMMUNITY_SIZE,
    FONT_TITLE, FONT_LABEL, FONT_TICK,
    COLOR_SBS2_HIGH, COLOR_CNV_HIGH, COLOR_NORMAL,
    banner, log, ensure_dir, load_adata, load_groups, get_gene_expression,
)

# Derived font sizes for annotation and legend
FONT_ANNOT = FONT_TICK - 2    # 26
FONT_LEGEND = FONT_TICK - 4   # 24

plt.rcParams.update({
    'font.size': FONT_TICK, 'axes.titlesize': FONT_TITLE,
    'axes.labelsize': FONT_LABEL, 'xtick.labelsize': FONT_TICK,
    'ytick.labelsize': FONT_TICK, 'legend.fontsize': FONT_LEGEND,
})

def short(p):
    return p.replace('Patient ', '')

# Variant sharing tier colors
TIER_COLORS = {
    'Universal':          '#2d3436',
    'Broadly shared':     '#3cb44b',
    'Partially shared':   '#4363d8',
    'Patient-specific':   '#f58231',
    'HC-exclusive':       '#e6194b',
}
TIER_ORDER = ['Universal', 'Broadly shared', 'Partially shared',
              'HC-exclusive', 'Patient-specific']

CHR_LENGTHS = {
    'chr1': 248.956, 'chr2': 242.194, 'chr3': 198.296, 'chr4': 190.215,
    'chr5': 181.538, 'chr6': 170.805, 'chr7': 159.346, 'chr8': 145.138,
    'chr9': 138.395, 'chr10': 133.797, 'chr11': 135.087, 'chr12': 133.275,
    'chr13': 114.364, 'chr14': 107.044, 'chr15': 101.991, 'chr16': 90.339,
    'chr17': 83.257, 'chr18': 80.373, 'chr19': 58.617, 'chr20': 64.444,
    'chr21': 46.710, 'chr22': 50.818, 'chrX': 156.041, 'chrY': 57.228,
}
CHR_ORDER = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


def load_all_data():
    banner("LOADING ALL DATA")
    adata = load_adata()
    sbs2_high, normal = load_groups()

    enrich_path = os.path.join(DIR_00_DIAG, "patient_enrichment_SBS2_HIGH_v2.tsv")
    enrichment_df = pd.read_csv(enrich_path, sep='\t') if os.path.exists(enrich_path) else None

    log(f"  Loading SComatic: {SCOMATIC_PATH}")
    mut_df = pd.read_csv(SCOMATIC_PATH, sep='\t')
    mut_df = mut_df[mut_df['REF'] != mut_df['Base_observed']].copy()
    mut_df['variant_id'] = (mut_df['#CHROM'].astype(str) + ':' +
                            mut_df['Start'].astype(str) + ':' +
                            mut_df['REF'].astype(str) + '>' +
                            mut_df['Base_observed'].astype(str))
    bc_to_patient = adata.obs[PATIENT_COL].to_dict()
    mut_df['patient'] = mut_df['CB'].map(bc_to_patient)
    mut_df = mut_df[mut_df['patient'].notna()].copy()
    log(f"  Variants: {mut_df['variant_id'].nunique():,}")
    return adata, sbs2_high, enrichment_df, mut_df


# ── PANEL A ──────────────────────────────────────────────────────────────────
def plot_panel_a(adata, high_cells, enrichment_df, out_dir):
    banner("PANEL A: PATIENT DISTRIBUTION (sorted by HIGH cell count)")
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell']
    patients_all = basal.obs[PATIENT_COL].value_counts()
    high_in = high_cells & set(adata.obs_names)
    basal_obs = basal.obs.copy()
    basal_obs['is_HIGH'] = basal_obs.index.isin(high_in)
    patients_high = basal_obs[basal_obs['is_HIGH']].groupby(PATIENT_COL, observed=True).size()
    all_p = sorted(patients_all.index)
    order = sorted(all_p, key=lambda p: patients_high.get(p, 0))

    fig, ax = plt.subplots(figsize=(16, 10))
    y = range(len(order))
    n_other = [patients_all.get(p, 0) - patients_high.get(p, 0) for p in order]
    n_high = [patients_high.get(p, 0) for p in order]
    ax.barh(list(y), n_other, color='#d4d4d4', edgecolor='black',
            linewidth=0.5, label='Other basal')
    ax.barh(list(y), n_high, left=n_other, color=COLOR_SBS2_HIGH,
            edgecolor='black', linewidth=0.5, label='SBS2-HIGH')
    ax.set_yticks(list(y))
    ax.set_yticklabels([short(p) for p in order])
    ax.set_xlabel('Number of Basal Cells')
    ax.set_title('Patient Distribution of SBS2-HIGH Basal Cells')
    ax.legend(loc='lower right')
    for i, p in enumerate(order):
        nh = patients_high.get(p, 0)
        if nh > 0:
            ax.text(patients_all.get(p, 0) + 50, i, f'{nh}', va='center',
                    fontsize=FONT_ANNOT, color=COLOR_SBS2_HIGH)
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_A_Patient_Distribution.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved Panel A")


# ── PANEL B ──────────────────────────────────────────────────────────────────
def plot_panel_b(adata, out_dir):
    banner("PANEL B: A3 EXPRESSION (sorted by expression)")
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell']
    rows = []
    for p in sorted(basal.obs[PATIENT_COL].unique()):
        mask = basal.obs[PATIENT_COL] == p
        subset = basal[mask]
        row = {'patient': p}
        for gene in ['APOBEC3A', 'APOBEC3B']:
            expr = get_gene_expression(subset, gene)
            row[gene] = expr.mean() if expr is not None else 0.0
        rows.append(row)
    a3_df = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(28, 10))
    for ax, gene, title in zip(axes, ['APOBEC3A', 'APOBEC3B'], ['A3A', 'A3B']):
        sdf = a3_df.sort_values(gene, ascending=True)
        colors = [COLOR_SBS2_HIGH if p in HIGH_CONTRIBUTORS
                  else COLOR_NORMAL for p in sdf['patient']]
        ax.barh(range(len(sdf)), sdf[gene].values, color=colors,
                edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(sdf)))
        ax.set_yticklabels([short(p) for p in sdf['patient']])
        ax.set_xlabel(f'Mean {title} Expression (all basal cells)')
        ax.set_title(f'{title} Expression per Patient\n(sorted by expression)')
    axes[1].legend(
        handles=[Patch(fc=COLOR_SBS2_HIGH, ec='black', label='High contributors'),
                 Patch(fc=COLOR_NORMAL, ec='black', label='Other patients')],
        loc='lower right')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_B_A3_Expression.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved Panel B")


# ── PANEL C ──────────────────────────────────────────────────────────────────
def plot_panel_c(adata, high_cells, out_dir):
    banner("PANEL C: PCA + HIERARCHICAL CLUSTERING HEATMAP")
    high_in = high_cells & set(adata.obs_names)
    adata_high = adata[list(high_in)].copy()
    pl = adata_high.obs[PATIENT_COL].values
    up = sorted(set(pl))

    de_path = os.path.join(DE_GENES_DIR, "SC_selected_genes_filtered.csv")
    if os.path.exists(de_path):
        de_genes = [g for g in pd.read_csv(de_path)['gene']
                    if g in adata_high.var_names]
    else:
        de_genes = list(adata_high.var_names[:5000])
    log(f"  Using {len(de_genes)} genes")

    X = adata_high[:, de_genes].X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    X_s = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(50, X_s.shape[1], X_s.shape[0] - 1),
              random_state=RANDOM_SEED)
    X_pca = pca.fit_transform(X_s)
    var = pca.explained_variance_ratio_

    cm = plt.cm.tab20
    pcols = {p: cm(i / max(len(up) - 1, 1)) for i, p in enumerate(sorted(up))}

    fig, axes = plt.subplots(1, 2, figsize=(28, 12))

    # PCA scatter
    ax = axes[0]
    for p in up:
        m = pl == p
        mk = '*' if p in HIGH_CONTRIBUTORS else 'o'
        sz = 200 if p in HIGH_CONTRIBUTORS else 80
        ax.scatter(X_pca[m, 0], X_pca[m, 1], c=[pcols[p]], s=sz, marker=mk,
                   alpha=0.6, edgecolors='black', linewidth=0.5, label=short(p))
    ax.set_xlabel(f'PC1 ({var[0] * 100:.1f}%)')
    ax.set_ylabel(f'PC2 ({var[1] * 100:.1f}%)')
    ax.set_title('SBS2-HIGH Cells by Patient')
    ax.legend(fontsize=FONT_LEGEND - 4, ncol=2, loc='best', markerscale=1.5)

    # Correlation heatmap
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
                    fontsize=FONT_ANNOT - 4,
                    color='white' if co.values[i, j] < 0.7 else 'black')
    plt.colorbar(im, ax=ax, label='Spearman rho', shrink=0.7)
    ax.set_title('Mean HIGH-Cell Profile Correlation\n(hierarchically clustered)')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_C_PCA_Correlation.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved Panel C")


# ── PANEL D ──────────────────────────────────────────────────────────────────
def plot_panel_d(out_dir):
    """LOPO sensitivity: 4-metric comparison (V4 pipeline).

    Metrics:
      1. Activating chain recovery (out of 9 genes)
      2. Community structure ARI with full analysis
      3. Gene overlap Jaccard with full analysis
      4. A3 wall integrity (% negative DIFF edges for A3A)

    Each metric shown as grouped bars: Full baseline + 3 LOPO runs.
    """
    banner("PANEL D: LOPO SENSITIVITY (V4 metrics)")

    # Full analysis baselines (by definition)
    lopo_patients = ['SC029', 'SC013', 'SC001']
    configs = ['Full\n(n=546)'] + [f'\u2212{p}' for p in lopo_patients]

    # Initialize with full-analysis reference values
    chain_vals = [len(ACTIVATING_CHAIN_GENES)]  # 9/9
    ari_vals = [1.0]
    jaccard_vals = [1.0]
    wall_vals = [100.0]
    cells_removed = [0]

    # Load LOPO results
    lopo_found = 0
    for patient in lopo_patients:
        lp = os.path.join(DIR_03_SENSITIVITY, f"LOPO_{patient}",
                          f"LOPO_{patient}_summary.tsv")
        if os.path.exists(lp):
            r = pd.read_csv(lp, sep='\t').iloc[0]
            chain_vals.append(int(r.get('activating_chain_recovered', 0)))
            ari_vals.append(float(r.get('community_ari_with_full', 0)))
            jaccard_vals.append(float(r.get('jaccard_with_full', 0)))
            wall_pct = r.get('a3a_wall_pct_neg', None)
            wall_vals.append(float(wall_pct) if pd.notna(wall_pct) else 0.0)
            cells_removed.append(int(r.get('high_removed', 0)))
            lopo_found += 1
            log(f"  Loaded LOPO_{patient}: chain={chain_vals[-1]}/9, "
                f"ARI={ari_vals[-1]:.3f}, Jaccard={jaccard_vals[-1]:.3f}, "
                f"wall={wall_vals[-1]:.0f}%")
        else:
            log(f"  WARNING: LOPO_{patient} summary not found, using zeros")
            chain_vals.append(0)
            ari_vals.append(0)
            jaccard_vals.append(0)
            wall_vals.append(0)
            cells_removed.append(0)

    if lopo_found == 0:
        log("  No LOPO results found. Panel D will show baselines only.")

    # Update config labels with cell counts
    for i, patient in enumerate(lopo_patients):
        if cells_removed[i + 1] > 0:
            n_remaining = 546 - cells_removed[i + 1]
            configs[i + 1] = f'\u2212{patient}\n(n={n_remaining})'

    # Plot: 2x2 grid of bar charts
    fig, axes = plt.subplots(2, 2, figsize=(24, 18))
    x = np.arange(len(configs))
    bar_width = 0.55

    metrics = [
        (axes[0, 0], chain_vals, 'Activating Chain Recovery',
         'Genes Recovered (of 9)', len(ACTIVATING_CHAIN_GENES),
         '#f58231', True),
        (axes[0, 1], ari_vals, 'Community Structure Similarity',
         'Adjusted Rand Index', 1.0,
         '#4363d8', False),
        (axes[1, 0], jaccard_vals, 'Gene Overlap with Full Analysis',
         'Jaccard Index', 1.0,
         '#3cb44b', False),
        (axes[1, 1], wall_vals, 'A3A Wall Integrity',
         '% Negative DIFF Edges', 100.0,
         '#2d3436', False),
    ]

    for ax, vals, title, ylabel, baseline, color, is_integer in metrics:
        bars = ax.bar(x, vals, width=bar_width, color=color,
                      edgecolor='black', linewidth=1.5)

        # Color the full-analysis bar differently
        bars[0].set_facecolor('#d4d4d4')
        bars[0].set_edgecolor('black')

        ax.axhline(baseline, color=COLOR_SBS2_HIGH, linestyle='--',
                   linewidth=2, alpha=0.7, label='Full analysis')

        ax.set_xticks(x)
        ax.set_xticklabels(configs, fontsize=FONT_TICK - 4)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        # Value labels on bars
        for i, v in enumerate(vals):
            if is_integer:
                label = f'{int(v)}'
            elif v >= 10:
                label = f'{v:.0f}%'
            else:
                label = f'{v:.3f}'
            ax.text(i, v + (baseline * 0.02), label, ha='center',
                    va='bottom', fontsize=FONT_ANNOT, fontweight='bold')

        # Set y-axis limits with padding
        if is_integer:
            ax.set_ylim(0, baseline + 1.5)
        elif baseline == 100:
            ax.set_ylim(0, 115)
        else:
            ax.set_ylim(0, baseline + 0.15)

    plt.suptitle('Leave-One-Patient-Out Network Sensitivity',
                 fontsize=FONT_TITLE, fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_D_LOPO_Sensitivity.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved Panel D")


# ── PANELS E-F: 5-TIER VARIANT SHARING ───────────────────────────────────────
def variant_sharing_analysis(mut_df, out_dir):
    banner("PANELS E-F: VARIANT SHARING TIER ANALYSIS")
    patients = sorted(mut_df['patient'].unique())
    n_patients = len(patients)
    hc_set = set(HIGH_CONTRIBUTORS)
    var_patient = mut_df.groupby('variant_id')['patient'].apply(set).to_dict()
    all_variants = sorted(var_patient.keys())
    log(f"  {n_patients} patients, {len(all_variants):,} variants")

    # Assign tiers
    variant_tier = {}
    tier_counts = Counter()
    for v in all_variants:
        pts = var_patient[v]
        n = len(pts)
        if pts <= hc_set and n >= 2:
            tier = 'HC-exclusive'
        elif n >= 10:
            tier = 'Universal'
        elif n >= 5:
            tier = 'Broadly shared'
        elif n >= 2:
            tier = 'Partially shared'
        else:
            tier = 'Patient-specific'
        variant_tier[v] = tier
        tier_counts[tier] += 1

    log(f"\n  Tier distribution:")
    for t in TIER_ORDER:
        log(f"    {t}: {tier_counts[t]:,}")

    # Compute how many patients share variants in each tier
    tier_patient_overlap = {}
    for t in TIER_ORDER:
        tier_vars = [v for v, tt in variant_tier.items() if tt == t]
        patients_in_tier = set()
        for v in tier_vars:
            patients_in_tier |= var_patient[v]
        tier_patient_overlap[t] = len(patients_in_tier)

    tier_df = pd.DataFrame([
        {'variant_id': v, 'tier': t, 'n_patients': len(var_patient[v]),
         'patients': ','.join(sorted(short(p) for p in var_patient[v]))}
        for v, t in variant_tier.items()
    ])
    tier_df.to_csv(os.path.join(out_dir, "variant_sharing_tiers.tsv"),
                   sep='\t', index=False)

    hc_ex = tier_df[tier_df['tier'] == 'HC-exclusive']
    if len(hc_ex) > 0:
        log(f"\n  HC-exclusive variants ({len(hc_ex)}):")
        for _, r in hc_ex.iterrows():
            log(f"    {r['variant_id']} ({r['patients']})")

    # ── PANEL E ──────────────────────────────────────────────────────────────
    banner("  PANEL E: VARIANT SHARING STACKED BAR")
    patient_tier_counts = {p: Counter() for p in patients}
    for v, tier in variant_tier.items():
        for p in var_patient[v]:
            patient_tier_counts[p][tier] += 1

    count_matrix = np.zeros((n_patients, len(TIER_ORDER)))
    for i, p in enumerate(patients):
        for j, t in enumerate(TIER_ORDER):
            count_matrix[i, j] = patient_tier_counts[p].get(t, 0)
    row_sums = count_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    prop_matrix = count_matrix / row_sums

    # Sort by HIGH cell count
    enrich_path = os.path.join(DIR_00_DIAG, "patient_enrichment_SBS2_HIGH_v2.tsv")
    if os.path.exists(enrich_path):
        edf = pd.read_csv(enrich_path, sep='\t')
        hc_map = {r['patient']: r['n_high'] for _, r in edf.iterrows()}
    else:
        hc_map = {p: 0 for p in patients}
    sort_idx = sorted(range(n_patients),
                      key=lambda i: hc_map.get(patients[i], 0))
    sorted_patients = [patients[sort_idx[i]] for i in range(n_patients)]

    fig, ax = plt.subplots(figsize=(20, 12))
    bottom = np.zeros(n_patients)
    for j, t in enumerate(TIER_ORDER):
        color = TIER_COLORS[t]
        values = [prop_matrix[sort_idx[i], j] for i in range(n_patients)]
        label = (f'{t} ({tier_counts[t]:,} variants, '
                 f'{tier_patient_overlap[t]} patients)')
        ax.barh(range(n_patients), values, left=bottom.tolist(),
                color=color, edgecolor='black', linewidth=0.3, label=label)
        bottom += values
    ax.set_yticks(range(n_patients))
    ax.set_yticklabels([short(p) for p in sorted_patients])
    ax.set_xlabel('Proportion of Variants')
    ax.set_title('Variant Sharing Spectrum per Patient\n'
                 '(sorted by SBS2-HIGH cell count)')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
              title='Sharing Tier', fontsize=FONT_LEGEND - 4)
    for i, p in enumerate(sorted_patients):
        if p in HIGH_CONTRIBUTORS:
            ax.get_yticklabels()[i].set_color(COLOR_SBS2_HIGH)
            ax.get_yticklabels()[i].set_fontweight('bold')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_E_Variant_Sharing.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved Panel E")

    # Print table
    log(f"\n  Per-patient tier counts:")
    log(f"  {'Patient':<12s}  {'Univ':>8s}  {'Broad':>8s}  "
        f"{'Partial':>8s}  {'HC-ex':>8s}  {'Specif':>8s}  {'Total':>8s}")
    for i in reversed(range(n_patients)):
        p = sorted_patients[i]
        pi = sort_idx[i]
        c = [int(count_matrix[pi, j]) for j in range(len(TIER_ORDER))]
        flag = " <<<" if p in HIGH_CONTRIBUTORS else ""
        log(f"  {short(p):<12s}  {c[0]:>8d}  {c[1]:>8d}  {c[2]:>8d}  "
            f"{c[3]:>8d}  {c[4]:>8d}  {sum(c):>8d}{flag}")

    # ── PANEL F ──────────────────────────────────────────────────────────────
    banner("  PANEL F: GENOME IDEOGRAM (sharing tiers)")
    var_positions = []
    for v, t in variant_tier.items():
        parts = v.split(':')
        if len(parts) >= 2:
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if chrom in CHR_LENGTHS:
                var_positions.append({'chrom': chrom, 'pos': pos, 'tier': t})
    pos_df = pd.DataFrame(var_positions)
    log(f"  {len(pos_df):,} variants with positions")

    chroms = [c for c in CHR_ORDER if c in CHR_LENGTHS]
    n_chr = len(chroms)
    max_len = max(CHR_LENGTHS.values())

    fig, ax = plt.subplots(figsize=(28, 16))
    chr_height = 0.6
    chr_gap = 1.2

    for i, chrom in enumerate(chroms):
        yc = (n_chr - i - 1) * chr_gap
        xl = CHR_LENGTHS[chrom] / max_len * 100
        ax.add_patch(FancyBboxPatch(
            (0, yc - chr_height / 2), xl, chr_height,
            boxstyle="round,pad=0.1", facecolor='#f0f0f0',
            edgecolor='black', linewidth=1))
        ax.text(-2, yc, chrom.replace('chr', ''), ha='right', va='center',
                fontsize=FONT_TICK - 4, fontweight='bold')

    # Plot tiers: background first, important last
    plot_config = [
        ('Patient-specific',   0.15,  6, 1.0),
        ('Partially shared',   0.25,  6, 1.0),
        ('Broadly shared',     0.80, 10, 2.0),
        ('Universal',          0.90, 12, 2.5),
        ('HC-exclusive',       1.00, 16, 3.0),
    ]
    for tier, alpha, ms, mew in plot_config:
        tv = pos_df[pos_df['tier'] == tier]
        color = TIER_COLORS[tier]
        for ci, chrom in enumerate(chroms):
            yc = (n_chr - ci - 1) * chr_gap
            xl = CHR_LENGTHS[chrom] / max_len * 100
            cv = tv[tv['chrom'] == chrom]
            if len(cv) == 0:
                continue
            xs = cv['pos'].values / (CHR_LENGTHS[chrom] * 1e6) * xl
            ax.plot(xs, np.full(len(xs), yc), '|', color=color,
                    markersize=ms, markeredgewidth=mew, alpha=alpha)

    ax.set_xlim(-5, 105)
    ax.set_ylim(-chr_gap, n_chr * chr_gap)
    ax.set_xlabel('Genomic Position (normalized)')
    ax.set_title('Variant Distribution Across the Genome\n'
                 '(colored by sharing tier)')
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    legend_el = [
        Patch(fc=TIER_COLORS[t], ec='black',
              label=f'{t} ({tier_counts[t]:,} variants, '
                    f'{tier_patient_overlap[t]} patients)')
        for t in TIER_ORDER
    ]
    ax.legend(handles=legend_el, loc='lower right', title='Sharing Tier',
              fontsize=FONT_LEGEND - 2, title_fontsize=FONT_LEGEND, ncol=1)
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_F_Ideogram.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved Panel F")


# ── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    banner("SUPPLEMENTAL FIGURE: PATIENT-SPECIFIC EFFECTS (v4)")
    out_dir = ensure_dir(FIGURE_5_PANELS)
    adata, sbs2_high, enrichment_df, mut_df = load_all_data()
    plot_panel_a(adata, sbs2_high, enrichment_df, out_dir)
    plot_panel_b(adata, out_dir)
    plot_panel_c(adata, sbs2_high, out_dir)
    plot_panel_d(out_dir)
    variant_sharing_analysis(mut_df, out_dir)
    log(f"\nAll panels saved to: {out_dir}")


if __name__ == "__main__":
    main()
