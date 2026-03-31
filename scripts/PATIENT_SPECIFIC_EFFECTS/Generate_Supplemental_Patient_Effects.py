#!/usr/bin/env python3
"""
Generate_Supplemental_Patient_Effects_v3.py
=============================================

v3 changes (from v2):
  - Panel A: sorted by HIGH cell count (descending, top = most HIGH cells)
  - Panels E-F: replaced Leiden community detection with 5-tier sharing
    classification (universal, broadly shared, partially shared,
    patient-specific, high-contributor exclusive)
  - Panels B, C, D: unchanged from v2

Usage:
  conda run -n NETWORK python Generate_Supplemental_Patient_Effects_v3.py
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

from patient_config import *

# =============================================================================
# GLOBAL STYLE
# =============================================================================
FONT_TITLE  = 34
FONT_LABEL  = 30
FONT_TICK   = 28
FONT_ANNOT  = 26
FONT_LEGEND = 24

plt.rcParams.update({
    'font.size': FONT_TICK,
    'axes.titlesize': FONT_TITLE,
    'axes.labelsize': FONT_LABEL,
    'xtick.labelsize': FONT_TICK,
    'ytick.labelsize': FONT_TICK,
    'legend.fontsize': FONT_LEGEND,
})

def short(p): return p.replace('Patient ', '')

# 5-tier sharing colors
TIER_COLORS = {
    'Universal':          '#2d3436',   # dark gray — germline/hotspot
    'Broadly shared':     '#4363d8',   # blue
    'Partially shared':   '#3cb44b',   # green
    'Patient-specific':   '#f58231',   # orange
    'HC-exclusive':       '#e6194b',   # red — SC029+SC013+SC001 only
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
    high_cells, low_cells = load_groups()

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

    return adata, high_cells, low_cells, enrichment_df, mut_df


# =============================================================================
# PANEL A: Patient Distribution — SORTED BY HIGH CELL COUNT (v3 fix)
# =============================================================================
def plot_panel_a(adata, high_cells, enrichment_df, out_dir):
    banner("PANEL A: PATIENT DISTRIBUTION (sorted by HIGH cell count)")
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell']
    patients_all = basal.obs[PATIENT_COL].value_counts()
    high_in = high_cells & set(adata.obs_names)
    basal_obs = basal.obs.copy()
    basal_obs['is_HIGH'] = basal_obs.index.isin(high_in)
    patients_high = basal_obs[basal_obs['is_HIGH']].groupby(PATIENT_COL, observed=True).size()

    # Sort by HIGH cell count (ascending so highest is at top of barh)
    all_patients = sorted(patients_all.index)
    order = sorted(all_patients, key=lambda p: patients_high.get(p, 0))

    fig, ax = plt.subplots(figsize=(16, 10))
    y = range(len(order))
    n_other = [patients_all.get(p, 0) - patients_high.get(p, 0) for p in order]
    n_high = [patients_high.get(p, 0) for p in order]
    ax.barh(list(y), n_other, color='#d4d4d4', edgecolor='black', linewidth=0.5, label='Other basal')
    ax.barh(list(y), n_high, left=n_other, color='#ed6a5a', edgecolor='black', linewidth=0.5, label='SBS2-HIGH')
    ax.set_yticks(list(y))
    ax.set_yticklabels([short(p) for p in order])
    ax.set_xlabel('Number of Basal Cells')
    ax.set_title('Patient Distribution of SBS2-HIGH Basal Cells')
    ax.legend(loc='lower right')
    for i, p in enumerate(order):
        nh = patients_high.get(p, 0)
        if nh > 0:
            ax.text(patients_all.get(p, 0) + 50, i, f'{nh}', va='center',
                   fontsize=FONT_ANNOT, color='#ed6a5a')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_A_Patient_Distribution.{ext}"),
                   dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved Panel A")


# =============================================================================
# PANEL B: A3 Expression (unchanged from v2)
# =============================================================================
def plot_panel_b(adata, enrichment_df, out_dir):
    banner("PANEL B: A3 EXPRESSION (sorted by expression)")
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell']
    patients = sorted(basal.obs[PATIENT_COL].unique())
    rows = []
    for p in patients:
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
        sorted_df = a3_df.sort_values(gene, ascending=True)
        colors = ['#ed6a5a' if p in HIGH_CONTRIBUTORS else '#5e81ac' for p in sorted_df['patient']]
        y = range(len(sorted_df))
        ax.barh(list(y), sorted_df[gene].values, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_yticks(list(y))
        ax.set_yticklabels([short(p) for p in sorted_df['patient']])
        ax.set_xlabel(f'Mean {title} Expression (all basal cells)')
        ax.set_title(f'{title} Expression per Patient\n(sorted by expression)')
    legend_el = [Patch(fc='#ed6a5a', ec='black', label='High contributors (SC029, SC013, SC001)'),
                 Patch(fc='#5e81ac', ec='black', label='Other patients')]
    axes[1].legend(handles=legend_el, loc='lower right')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_B_A3_Expression.{ext}"),
                   dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved Panel B")


# =============================================================================
# PANEL C: PCA + Hierarchically Clustered Heatmap (unchanged from v2)
# =============================================================================
def plot_panel_c(adata, high_cells, out_dir):
    banner("PANEL C: PCA + HIERARCHICAL CLUSTERING HEATMAP")
    high_in = high_cells & set(adata.obs_names)
    adata_high = adata[list(high_in)].copy()
    patient_labels = adata_high.obs[PATIENT_COL].values
    unique_patients = sorted(set(patient_labels))

    de_path = os.path.join(DE_GENES_DIR, "SC_selected_genes_filtered.csv")
    if os.path.exists(de_path):
        de_df = pd.read_csv(de_path)
        de_genes = [g for g in de_df['gene'].tolist() if g in adata_high.var_names]
    else:
        de_genes = list(adata_high.var_names[:5000])
    log(f"  Using {len(de_genes)} genes")

    adata_de = adata_high[:, de_genes]
    X = adata_de.X
    if scipy.sparse.issparse(X): X = X.toarray()
    X_s = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(50, X_s.shape[1], X_s.shape[0]-1), random_state=RANDOM_SEED)
    X_pca = pca.fit_transform(X_s)
    var = pca.explained_variance_ratio_

    cmap = plt.cm.tab20
    pcols = {p: cmap(i / max(len(unique_patients)-1, 1)) for i, p in enumerate(sorted(unique_patients))}

    fig, axes = plt.subplots(1, 2, figsize=(28, 12))
    ax = axes[0]
    for p in unique_patients:
        m = patient_labels == p
        mk = '*' if p in HIGH_CONTRIBUTORS else 'o'
        sz = 200 if p in HIGH_CONTRIBUTORS else 80
        ax.scatter(X_pca[m, 0], X_pca[m, 1], c=[pcols[p]], s=sz, marker=mk,
                  alpha=0.6, edgecolors='black', linewidth=0.5, label=short(p))
    ax.set_xlabel(f'PC1 ({var[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({var[1]*100:.1f}%)')
    ax.set_title(f'SBS2-HIGH Cells by Patient\n(silhouette = 0.125)')
    ax.legend(fontsize=FONT_LEGEND-4, ncol=2, loc='best', markerscale=1.5)

    ax = axes[1]
    pts_corr = [p for p in unique_patients if (patient_labels == p).sum() >= 2]
    profiles = {p: X[patient_labels == p].mean(axis=0) for p in pts_corr}
    prof_df = pd.DataFrame(profiles)
    corr = prof_df.corr(method='spearman')
    dist_matrix = np.maximum(1 - corr.values, 0)
    np.fill_diagonal(dist_matrix, 0)
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')
    leaf_order = leaves_list(Z)
    ordered_patients = [corr.columns[i] for i in leaf_order]
    corr_ordered = corr.loc[ordered_patients, ordered_patients]

    labs = [short(p) for p in ordered_patients]
    im = ax.imshow(corr_ordered.values, cmap='RdYlBu_r', vmin=0.5, vmax=1)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=45, ha='right')
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs)
    for i, p in enumerate(ordered_patients):
        if p in HIGH_CONTRIBUTORS:
            ax.get_xticklabels()[i].set_color('#ed6a5a')
            ax.get_xticklabels()[i].set_fontweight('bold')
            ax.get_yticklabels()[i].set_color('#ed6a5a')
            ax.get_yticklabels()[i].set_fontweight('bold')
    for i in range(len(labs)):
        for j in range(len(labs)):
            ax.text(j, i, f'{corr_ordered.values[i,j]:.2f}', ha='center', va='center',
                   fontsize=FONT_ANNOT-4, color='white' if corr_ordered.values[i,j] < 0.7 else 'black')
    plt.colorbar(im, ax=ax, label='Spearman rho', shrink=0.7)
    ax.set_title('Mean HIGH-Cell Profile Correlation\n(hierarchically clustered)')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_C_PCA_Correlation.{ext}"),
                   dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved Panel C")


# =============================================================================
# PANEL D: LOPO Sensitivity (unchanged from v2)
# =============================================================================
def plot_panel_d(out_dir):
    banner("PANEL D: LOPO SENSITIVITY")
    results = [{'config': 'Full\n(n=546)', 'A3A_degree': 8, 'n_communities': 14, 'interactors': 4}]
    for patient, label in [('SC029', '−SC029\n(n=363)'), ('SC013', '−SC013\n(n=420)'),
                            ('SC001', '−SC001\n(n=486)')]:
        lopo_path = os.path.join(DIR_03_SENSITIVITY, f"LOPO_{patient}", f"LOPO_{patient}_summary.tsv")
        if os.path.exists(lopo_path):
            df = pd.read_csv(lopo_path, sep='\t'); row = df.iloc[0]
            n_int = len(eval(row['interactors_found'])) if isinstance(row['interactors_found'], str) else 0
            results.append({'config': label, 'A3A_degree': int(row['A3A_degree']),
                           'n_communities': int(row['n_communities']), 'interactors': n_int})
    tumor_path = os.path.join(DIR_03_SENSITIVITY, "tumor_only", "tumor_only_comparison.tsv")
    if os.path.exists(tumor_path):
        tdf = pd.read_csv(tumor_path, sep='\t')
        if 'approach' in tdf.columns:
            row_a = tdf[tdf['approach'] == 'A_original_LOW'].iloc[0]
        else:
            row_a = tdf.iloc[0]
        n_int = len(eval(row_a['interactors_found'])) if isinstance(row_a['interactors_found'], str) else 0
        results.append({'config': 'Tumor\nonly\n(n=521)', 'A3A_degree': int(row_a['A3A_degree']),
                       'n_communities': int(row_a['n_communities']), 'interactors': n_int})
    rdf = pd.DataFrame(results)
    configs = rdf['config'].tolist()

    fig, axes = plt.subplots(1, 3, figsize=(30, 10))
    metrics = [('A3A_degree', 'APOBEC3A Hub Degree', '#2d3436', 8),
               ('interactors', 'A3 Interactors in Network (of 4)', '#f58231', 4),
               ('n_communities', 'Leiden Communities', '#4363d8', 14)]
    for ax, (col, title, color, baseline) in zip(axes, metrics):
        ax.bar(range(len(configs)), rdf[col], color=color, edgecolor='black', linewidth=1.5, width=0.6)
        ax.set_xticks(range(len(configs))); ax.set_xticklabels(configs, fontsize=FONT_TICK-4)
        ax.set_ylabel(col.replace('_', ' ').title())
        ax.set_title(title)
        ax.axhline(baseline, color='#ed6a5a', linestyle='--', linewidth=2, label='Full analysis')
        ax.legend()
    plt.suptitle('Leave-One-Patient-Out Network Sensitivity', fontsize=FONT_TITLE, fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_D_LOPO_Sensitivity.{ext}"),
                   dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved Panel D")


# =============================================================================
# PANELS E-F: 5-TIER VARIANT SHARING (v3 — replaces Leiden)
# =============================================================================
def variant_sharing_analysis(mut_df, out_dir):
    banner("PANELS E-F: VARIANT SHARING TIER ANALYSIS")

    patients = sorted(mut_df['patient'].unique())
    n_patients = len(patients)
    hc_set = set(HIGH_CONTRIBUTORS)
    log(f"  {n_patients} patients, {mut_df['variant_id'].nunique():,} variants")

    # Build variant -> set of patients
    log("  Building variant-patient map...")
    var_patient = mut_df.groupby('variant_id')['patient'].apply(set).to_dict()
    all_variants = sorted(var_patient.keys())
    log(f"  {len(all_variants):,} unique variants")

    # Assign tiers
    log("  Assigning sharing tiers...")
    variant_tier = {}
    tier_counts = Counter()

    for v in all_variants:
        pts = var_patient[v]
        n = len(pts)

        # Check HC-exclusive first: found ONLY in high-contributors
        if pts <= hc_set and n >= 2:
            tier = 'HC-exclusive'
        elif n >= 10:  # >= ~70% of 14 patients
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
    for tier in TIER_ORDER:
        log(f"    {tier}: {tier_counts[tier]:,} variants")
    log(f"    Total: {sum(tier_counts.values()):,}")

    # Save tier assignments
    tier_df = pd.DataFrame([
        {'variant_id': v, 'tier': t, 'n_patients': len(var_patient[v]),
         'patients': ','.join(sorted(short(p) for p in var_patient[v]))}
        for v, t in variant_tier.items()
    ])
    tier_df.to_csv(os.path.join(out_dir, "variant_sharing_tiers.tsv"), sep='\t', index=False)
    log(f"  Saved tier assignments")

    # Detail the HC-exclusive variants
    hc_exclusive = tier_df[tier_df['tier'] == 'HC-exclusive']
    if len(hc_exclusive) > 0:
        log(f"\n  HC-exclusive variants ({len(hc_exclusive)}):")
        for _, row in hc_exclusive.head(20).iterrows():
            log(f"    {row['variant_id']} ({row['patients']})")

    # =========================================================================
    # PANEL E: Stacked bar — tier proportions per patient
    # =========================================================================
    banner("  PANEL E: VARIANT SHARING STACKED BAR")

    # Count tiers per patient
    patient_tier_counts = {p: Counter() for p in patients}
    for v, tier in variant_tier.items():
        for p in var_patient[v]:
            patient_tier_counts[p][tier] += 1

    # Build matrix
    count_matrix = np.zeros((n_patients, len(TIER_ORDER)))
    for i, p in enumerate(patients):
        for j, tier in enumerate(TIER_ORDER):
            count_matrix[i, j] = patient_tier_counts[p].get(tier, 0)

    # Normalize to proportions
    row_sums = count_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    prop_matrix = count_matrix / row_sums

    # Sort patients by HIGH cell count
    enrich_path = os.path.join(DIR_00_DIAG, "patient_enrichment_SBS2_HIGH_v2.tsv")
    if os.path.exists(enrich_path):
        edf = pd.read_csv(enrich_path, sep='\t')
        high_count_map = {r['patient']: r['n_high'] for _, r in edf.iterrows()}
    else:
        high_count_map = {p: 0 for p in patients}
    sort_idx = sorted(range(n_patients), key=lambda i: high_count_map.get(patients[i], 0))
    sorted_patients = [patients[sort_idx[i]] for i in range(n_patients)]

    fig, ax = plt.subplots(figsize=(20, 12))
    bottom = np.zeros(n_patients)

    for j, tier in enumerate(TIER_ORDER):
        color = TIER_COLORS[tier]
        values = [prop_matrix[sort_idx[i], j] for i in range(n_patients)]
        ax.barh(range(n_patients), values, left=bottom.tolist(),
               color=color, edgecolor='black', linewidth=0.3, label=f'{tier} ({tier_counts[tier]:,})')
        bottom += values

    ax.set_yticks(range(n_patients))
    ax.set_yticklabels([short(p) for p in sorted_patients])
    ax.set_xlabel('Proportion of Variants')
    ax.set_title('Variant Sharing Spectrum per Patient\n(sorted by SBS2-HIGH cell count)')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), title='Sharing Tier')

    for i, p in enumerate(sorted_patients):
        if p in HIGH_CONTRIBUTORS:
            ax.get_yticklabels()[i].set_color('#ed6a5a')
            ax.get_yticklabels()[i].set_fontweight('bold')

    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_E_Variant_Sharing.{ext}"),
                   dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved Panel E")

    # Also print raw counts per patient for the log
    log(f"\n  Per-patient tier counts:")
    log(f"  {'Patient':<15s}  {'Universal':>10s}  {'Broad':>10s}  {'Partial':>10s}  "
        f"{'HC-excl':>10s}  {'Specific':>10s}  {'Total':>10s}")
    for i in reversed(range(n_patients)):  # print high-to-low
        p = sorted_patients[i]
        pi = sort_idx[i]
        counts = [int(count_matrix[pi, j]) for j in range(len(TIER_ORDER))]
        total = int(sum(counts))
        flag = " <<<" if p in HIGH_CONTRIBUTORS else ""
        log(f"  {short(p):<15s}  {counts[0]:>10d}  {counts[1]:>10d}  {counts[2]:>10d}  "
            f"{counts[3]:>10d}  {counts[4]:>10d}  {total:>10d}{flag}")

    # =========================================================================
    # PANEL F: Genome ideogram colored by sharing tier
    # =========================================================================
    banner("  PANEL F: GENOME IDEOGRAM (sharing tiers)")

    var_positions = []
    for v, tier in variant_tier.items():
        parts = v.split(':')
        if len(parts) >= 2:
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if chrom in CHR_LENGTHS:
                var_positions.append({'chrom': chrom, 'pos': pos, 'tier': tier})
    pos_df = pd.DataFrame(var_positions)
    log(f"  {len(pos_df):,} variants with genomic positions")

    chroms_to_plot = [c for c in CHR_ORDER if c in CHR_LENGTHS]
    n_chr = len(chroms_to_plot)
    max_len = max(CHR_LENGTHS.values())

    fig, ax = plt.subplots(figsize=(28, 16))
    chr_height = 0.6
    chr_gap = 1.2

    # Draw chromosomes and variants
    # Plot tiers in order so rarer tiers render on top
    plot_order = ['Universal', 'Broadly shared', 'Partially shared',
                  'Patient-specific', 'HC-exclusive']

    for i, chrom in enumerate(chroms_to_plot):
        y_center = (n_chr - i - 1) * chr_gap
        x_len = CHR_LENGTHS[chrom] / max_len * 100
        rect = FancyBboxPatch((0, y_center - chr_height/2), x_len, chr_height,
                              boxstyle="round,pad=0.1", facecolor='#f0f0f0',
                              edgecolor='black', linewidth=1)
        ax.add_patch(rect)
        ax.text(-2, y_center, chrom.replace('chr', ''), ha='right', va='center',
                fontsize=FONT_TICK-4, fontweight='bold')

    # Plot variants tier by tier (common first, rare on top)
    for tier in plot_order:
        tier_vars = pos_df[pos_df['tier'] == tier]
        color = TIER_COLORS[tier]
        alpha = 0.4 if tier == 'Patient-specific' else 0.7
        ms = 12 if tier == 'HC-exclusive' else 8
        mew = 2.5 if tier == 'HC-exclusive' else 1.5

        for chrom_idx, chrom in enumerate(chroms_to_plot):
            y_center = (n_chr - chrom_idx - 1) * chr_gap
            x_len = CHR_LENGTHS[chrom] / max_len * 100
            chr_vars = tier_vars[tier_vars['chrom'] == chrom]
            if len(chr_vars) == 0:
                continue
            xs = chr_vars['pos'].values / (CHR_LENGTHS[chrom] * 1e6) * x_len
            ys = np.full(len(xs), y_center)
            ax.plot(xs, ys, '|', color=color, markersize=ms,
                   markeredgewidth=mew, alpha=alpha)

    ax.set_xlim(-5, 105)
    ax.set_ylim(-chr_gap, n_chr * chr_gap)
    ax.set_xlabel('Genomic Position (normalized)')
    ax.set_title('Variant Distribution Across the Genome\n(colored by sharing tier)')
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    legend_el = [Patch(fc=TIER_COLORS[t], ec='black', label=f'{t} ({tier_counts[t]:,})')
                 for t in TIER_ORDER]
    ax.legend(handles=legend_el, loc='upper right', title='Sharing Tier',
             fontsize=FONT_LEGEND-2, title_fontsize=FONT_LEGEND, ncol=1)

    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Supp_Panel_F_Ideogram.{ext}"),
                   dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved Panel F")

    return variant_tier


# =============================================================================
# MAIN
# =============================================================================
def main():
    banner("SUPPLEMENTAL FIGURE: PATIENT-SPECIFIC EFFECTS (v3)")
    out_dir = ensure_dir(FIGURE_5_PANELS)
    adata, high_cells, low_cells, enrichment_df, mut_df = load_all_data()

    plot_panel_a(adata, high_cells, enrichment_df, out_dir)
    plot_panel_b(adata, enrichment_df, out_dir)
    plot_panel_c(adata, high_cells, out_dir)
    plot_panel_d(out_dir)
    variant_sharing_analysis(mut_df, out_dir)

    log(f"\nAll panels saved to: {out_dir}")

if __name__ == "__main__":
    main()
