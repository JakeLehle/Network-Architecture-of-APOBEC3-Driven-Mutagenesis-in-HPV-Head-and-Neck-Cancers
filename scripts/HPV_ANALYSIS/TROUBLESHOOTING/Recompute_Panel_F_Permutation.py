#!/usr/bin/env python3
"""
Recompute_Panel_F_Permutation.py
==================================

Replaces Mann-Whitney U with a permutation test on the difference of means
for Panel F (HPV16 lifecycle fraction violins) in Figure 6.

Rationale: HPV16 gene fractions are zero-inflated (E6: 90%+ zeros).
Mann-Whitney detects distributional shifts (e.g., different zero fractions)
that don't correspond to biologically meaningful differences in expression.
A permutation test on means directly tests whether the observed mean
difference is larger than expected by chance.

Approach:
  1. For each gene x pair, compute observed difference of means
  2. Shuffle group labels N_PERM times, recompute mean difference each time
  3. p-value = proportion of permutations with |mean diff| >= |observed|
  4. Apply figure-wide BH correction (Panel B Mann-Whitney + Panel D
     Mann-Whitney + Panel F permutation = same total test count)

Outputs:
  - Panel_F_permutation_report.tsv (all raw + adjusted p-values)
  - Panel_6F_permutation_violins.pdf/.png (updated violin panel)
  - Comparison table: Mann-Whitney q vs Permutation q for Panel F

Usage:
    conda run -n NETWORK python Recompute_Panel_F_Permutation.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
HPV_GENE_PATH = os.path.join(BASE_DIR, "data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv")
GROUP_PATH = os.path.join(BASE_DIR, "data/FIG_4/01_group_selection/three_group_assignments.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "data/FIG_6/FIGURE_6_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

N_PERM = 10000
SEED = 42
MIN_CELLS = 10

POP_ORDER = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
POP_LABELS = {'SBS2_HIGH': 'SBS2-HIGH', 'CNV_HIGH': 'CNV-HIGH', 'NORMAL': 'Normal'}
POP_COLORS = {'SBS2_HIGH': '#ed6a5a', 'CNV_HIGH': '#F6D155', 'NORMAL': '#4682b4'}

HPV16_PHASES = OrderedDict([
    ('Maintenance',   ['E1', 'E2']),
    ('Amplification', ['E4', 'E5']),
    ('Oncogene',      ['E6', 'E7']),
    ('Capsid',        ['L1', 'L2']),
])

PHASE_COLORS = {
    'Maintenance':   '#2166ac',
    'Amplification': '#92c5de',
    'Oncogene':      '#f4a582',
    'Capsid':        '#b2182b',
}

# Font sizes
FONT_TITLE = 32
FONT_AXIS  = 30
FONT_TICK  = 26
FONT_PHASE = 28
FONT_PVAL  = 14
FONT_STARS = 20
FONT_LEGEND = 22
DPI = 300


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def banner(title, char="="):
    print(f"\n{char * 80}\n  {title}\n{char * 80}", flush=True)


# =============================================================================
# PERMUTATION TEST
# =============================================================================

def permutation_test_means(vals_a, vals_b, n_perm=N_PERM, seed=SEED):
    """
    Two-sided permutation test on the difference of means.

    Returns:
        obs_diff: observed difference (mean_a - mean_b)
        p_value: two-sided p-value
    """
    rng = np.random.RandomState(seed)
    obs_diff = np.mean(vals_a) - np.mean(vals_b)
    pooled = np.concatenate([vals_a, vals_b])
    n_a = len(vals_a)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(pooled)
        perm_diff = np.mean(pooled[:n_a]) - np.mean(pooled[n_a:])
        if abs(perm_diff) >= abs(obs_diff):
            count += 1
    p_value = (count + 1) / (n_perm + 1)  # +1 to include observed
    return obs_diff, p_value


def p_to_stars(p):
    if p < 1e-4: return '****'
    if p < 1e-3: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return 'ns'

def format_q(q):
    if q < 1e-4:
        return f'q={q:.0e}'
    elif q < 0.05:
        return f'q={q:.3f}'
    else:
        return f'q={q:.2f}'


# =============================================================================
# VIOLIN PANEL (adapted from figure generation script)
# =============================================================================

def make_violin_panel(ax, data_dict, pop_order, pop_colors, pop_labels,
                      ylabel='', adjusted_pvals=None):
    """Draw violin + box overlay with permutation-test brackets."""
    positions = list(range(len(pop_order)))
    plot_data = [data_dict[p] for p in pop_order]
    valid_groups = [len(data_dict[p]) >= MIN_CELLS for p in pop_order]

    if any(valid_groups):
        parts = ax.violinplot(
            [d if v else [0] for d, v in zip(plot_data, valid_groups)],
            positions=positions,
            showmeans=False, showmedians=False, showextrema=False)
        for idx_pc, (pc, pop) in enumerate(zip(parts['bodies'], pop_order)):
            if valid_groups[idx_pc]:
                pc.set_facecolor(pop_colors[pop])
                pc.set_edgecolor('#333333')
                pc.set_linewidth(1.2)
                pc.set_alpha(0.8)
            else:
                pc.set_visible(False)

    bp = ax.boxplot(
        [d if v else [0] for d, v in zip(plot_data, valid_groups)],
        positions=positions, widths=0.15, patch_artist=True, showfliers=False,
        medianprops=dict(color='black', linewidth=2),
        whiskerprops=dict(color='black', linewidth=1),
        capprops=dict(color='black', linewidth=1))
    for idx_bp, patch in enumerate(bp['boxes']):
        if valid_groups[idx_bp]:
            patch.set_facecolor('white')
            patch.set_edgecolor('#333333')
            patch.set_linewidth(1)
            patch.set_alpha(0.8)
        else:
            patch.set_visible(False)
            for element in ['whiskers', 'caps', 'medians']:
                n_per = 2 if element != 'medians' else 1
                start = idx_bp * n_per
                end = start + n_per
                for line in bp[element][start:end]:
                    line.set_visible(False)

    for i, pop in enumerate(pop_order):
        if valid_groups[i] and len(data_dict[pop]) > 0:
            mean_val = np.mean(data_dict[pop])
            ax.scatter(i, mean_val, color='black', s=40, zorder=5, marker='D')
        elif not valid_groups[i]:
            ax.text(i, 0, 'N.D.', ha='center', va='center',
                    fontsize=FONT_PVAL + 2, fontstyle='italic', color='#888888')

    # Brackets
    if adjusted_pvals is not None and len(pop_order) >= 3:
        pairs = [(0, 1), (1, 2), (0, 2)]
        valid_maxes = [np.percentile(data_dict[p], 99.5)
                       for p in pop_order if len(data_dict[p]) > 0]
        valid_mins = [np.min(data_dict[p])
                      for p in pop_order if len(data_dict[p]) > 0]
        y_data_max = max(valid_maxes) if valid_maxes else 1
        y_data_min = min(valid_mins) if valid_mins else 0
        y_range = max(y_data_max - y_data_min, 0.001)
        bracket_gap = 0.12 * y_range

        for bracket_idx, (i, j) in enumerate(pairs):
            q = adjusted_pvals[bracket_idx]
            y_bar = y_data_max + (bracket_idx + 0.5) * bracket_gap
            y_tip = y_bar + 0.015 * y_range

            if np.isnan(q):
                ax.plot([i, i, j, j], [y_bar, y_tip, y_tip, y_bar],
                        color='#cccccc', linewidth=1.0, clip_on=False)
                ax.text((i + j) / 2.0, y_tip + 0.005 * y_range,
                        'N.D.', ha='center', va='bottom',
                        fontsize=FONT_PVAL, fontstyle='italic', color='#888888')
                continue

            ax.plot([i, i, j, j], [y_bar, y_tip, y_tip, y_bar],
                    color='black', linewidth=1.2, clip_on=False)
            stars = p_to_stars(q)
            q_text = format_q(q)
            ax.text((i + j) / 2.0, y_tip + 0.005 * y_range,
                    f'{stars}  {q_text}',
                    ha='center', va='bottom',
                    fontsize=FONT_PVAL, fontweight='bold')

        ax.set_ylim(bottom=y_data_min - 0.05 * y_range,
                    top=y_data_max + (len(pairs) + 1.2) * bracket_gap)

    ax.set_xticks(positions)
    ax.set_xticklabels([pop_labels[p] for p in pop_order],
                       fontsize=FONT_TICK - 6, rotation=35, ha='right')
    ax.tick_params(axis='y', labelsize=FONT_TICK - 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=FONT_AXIS - 6)


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("PANEL F: PERMUTATION TEST ON MEANS")

    # Load data
    groups_df = pd.read_csv(GROUP_PATH, sep="\t")
    groups_map = dict(zip(groups_df.iloc[:, 0], groups_df.iloc[:, 1]))

    hpv_genes = pd.read_csv(HPV_GENE_PATH, sep="\t", index_col=0)
    hpv_genes['group'] = hpv_genes.index.map(lambda x: groups_map.get(x, 'other'))

    total_col = 'total_hpv16_genome_reads'
    hpv_pop = hpv_genes[hpv_genes['group'].isin(POP_ORDER)].copy()
    hpv_pop = hpv_pop[hpv_pop[total_col] > 0].copy()

    log(f"HPV+ cells in three populations: {len(hpv_pop)}")
    for pop in POP_ORDER:
        log(f"  {pop}: {(hpv_pop['group'] == pop).sum()}")

    # Compute fractions
    all_genes = []
    for phase_genes in HPV16_PHASES.values():
        all_genes.extend(phase_genes)

    for gene in all_genes:
        if gene in hpv_pop.columns:
            hpv_pop[f'{gene}_frac'] = hpv_pop[gene] / hpv_pop[total_col]

    # Compute permutation tests for all gene x pair combinations
    banner("Computing permutation tests (n_perm={})".format(N_PERM))

    pairs = [(0, 1), (1, 2), (0, 2)]
    pair_names = ['SBS2 vs CNV', 'CNV vs NORM', 'SBS2 vs NORM']

    all_perm_raw = []
    all_mw_raw = []
    all_labels = []
    all_obs_diffs = []
    all_data_dicts = {}

    for gene in all_genes:
        frac_col = f'{gene}_frac'
        if frac_col not in hpv_pop.columns:
            continue

        vdata = {}
        for pop in POP_ORDER:
            mask = hpv_pop['group'] == pop
            vdata[pop] = hpv_pop.loc[mask, frac_col].values.astype(float)
        all_data_dicts[gene] = vdata

        for pi, (i, j) in enumerate(pairs):
            pop_i = POP_ORDER[i]
            pop_j = POP_ORDER[j]
            vals_i = vdata[pop_i]
            vals_j = vdata[pop_j]

            if len(vals_i) >= MIN_CELLS and len(vals_j) >= MIN_CELLS:
                obs_diff, perm_p = permutation_test_means(vals_i, vals_j)
                _, mw_p = mannwhitneyu(vals_i, vals_j, alternative='two-sided')
                all_perm_raw.append(perm_p)
                all_mw_raw.append(mw_p)
                all_obs_diffs.append(obs_diff)
            else:
                all_perm_raw.append(np.nan)
                all_mw_raw.append(np.nan)
                all_obs_diffs.append(np.nan)
            all_labels.append(f"{gene} | {pair_names[pi]}")

    # BH correction (Panel F only)
    perm_arr = np.array(all_perm_raw, dtype=float)
    mw_arr = np.array(all_mw_raw, dtype=float)
    valid = ~np.isnan(perm_arr)

    adj_perm = np.full_like(perm_arr, np.nan)
    adj_mw = np.full_like(mw_arr, np.nan)
    if valid.sum() > 0:
        _, adj_p, _, _ = multipletests(perm_arr[valid], method='fdr_bh')
        adj_perm[valid] = adj_p
        _, adj_m, _, _ = multipletests(mw_arr[valid], method='fdr_bh')
        adj_mw[valid] = adj_m

    # Report
    banner("COMPARISON: Mann-Whitney vs Permutation (Panel F)")
    log(f"\n  {'Label':<25} {'Mean diff':>10} {'MW raw':>12} {'MW q':>12} "
        f"{'Perm raw':>12} {'Perm q':>12} {'MW':>6} {'Perm':>6}")
    log(f"  {'-'*95}")

    for k in range(len(all_labels)):
        if np.isnan(all_perm_raw[k]):
            log(f"  {all_labels[k]:<25} {'N.D.':>10} {'N.D.':>12} {'N.D.':>12} "
                f"{'N.D.':>12} {'N.D.':>12}")
            continue
        md = all_obs_diffs[k] * 100  # convert to percentage points
        mw_stars = p_to_stars(adj_mw[k])
        perm_stars = p_to_stars(adj_perm[k])
        log(f"  {all_labels[k]:<25} {md:>+9.3f}% {all_mw_raw[k]:>12.4e} "
            f"{adj_mw[k]:>12.4e} {all_perm_raw[k]:>12.4e} {adj_perm[k]:>12.4e} "
            f"{mw_stars:>6} {perm_stars:>6}")

    # Phase-level permutation tests
    banner("PHASE-LEVEL PERMUTATION TESTS")

    phase_perm_raw = []
    phase_labels = []
    phase_obs_diffs = []

    for phase_name, phase_genes in HPV16_PHASES.items():
        frac_cols = [f'{g}_frac' for g in phase_genes if f'{g}_frac' in hpv_pop.columns]
        if not frac_cols:
            continue

        for pi, (i, j) in enumerate(pairs):
            pop_i = POP_ORDER[i]
            pop_j = POP_ORDER[j]
            vals_i = hpv_pop.loc[hpv_pop['group'] == pop_i, frac_cols].sum(axis=1).values
            vals_j = hpv_pop.loc[hpv_pop['group'] == pop_j, frac_cols].sum(axis=1).values

            if len(vals_i) >= MIN_CELLS and len(vals_j) >= MIN_CELLS:
                obs_diff, perm_p = permutation_test_means(vals_i, vals_j)
                phase_perm_raw.append(perm_p)
                phase_obs_diffs.append(obs_diff)
            else:
                phase_perm_raw.append(np.nan)
                phase_obs_diffs.append(np.nan)
            phase_labels.append(f"{phase_name} | {pair_names[pi]}")

    phase_arr = np.array(phase_perm_raw, dtype=float)
    valid_ph = ~np.isnan(phase_arr)
    adj_phase = np.full_like(phase_arr, np.nan)
    if valid_ph.sum() > 0:
        _, adj_ph, _, _ = multipletests(phase_arr[valid_ph], method='fdr_bh')
        adj_phase[valid_ph] = adj_ph

    log(f"\n  {'Phase comparison':<35} {'Mean diff':>10} {'Perm raw':>12} {'Perm q':>12} {'Stars':>8}")
    log(f"  {'-'*80}")
    for k in range(len(phase_labels)):
        if np.isnan(phase_perm_raw[k]):
            log(f"  {phase_labels[k]:<35} {'N.D.':>10} {'N.D.':>12} {'N.D.':>12}")
            continue
        md = phase_obs_diffs[k] * 100
        stars = p_to_stars(adj_phase[k])
        log(f"  {phase_labels[k]:<35} {md:>+9.2f}% {phase_perm_raw[k]:>12.4e} "
            f"{adj_phase[k]:>12.4e} {stars:>8}")

    # Save report TSV
    report_rows = []
    for k in range(len(all_labels)):
        report_rows.append({
            'comparison': all_labels[k],
            'mean_diff_pct': all_obs_diffs[k] * 100 if not np.isnan(all_obs_diffs[k]) else np.nan,
            'mw_raw_p': all_mw_raw[k],
            'mw_bh_q': adj_mw[k],
            'perm_raw_p': all_perm_raw[k],
            'perm_bh_q': adj_perm[k],
        })
    report_df = pd.DataFrame(report_rows)
    report_path = os.path.join(OUTPUT_DIR, "Panel_F_permutation_report.tsv")
    report_df.to_csv(report_path, sep="\t", index=False)
    log(f"\n  Report saved: {report_path}")

    # =========================================================================
    # REGENERATE PANEL F VIOLINS WITH PERMUTATION Q-VALUES
    # =========================================================================
    banner("REGENERATING PANEL F VIOLINS (permutation q-values)")

    fig = plt.figure(figsize=(48, 12))
    gs = gridspec.GridSpec(1, 8, wspace=0.40)

    gene_idx = 0
    f_axes = []
    for phase_name, phase_genes in HPV16_PHASES.items():
        for gene in phase_genes:
            ax = fig.add_subplot(gs[0, gene_idx])
            f_axes.append(ax)

            vdata = all_data_dicts.get(gene, {pop: np.array([0.0]) for pop in POP_ORDER})

            # Get this gene's permutation q-values (3 pairs)
            base_idx = gene_idx * 3
            gene_qvals = [adj_perm[base_idx + pi] for pi in range(3)]

            make_violin_panel(ax, vdata, POP_ORDER, POP_COLORS, POP_LABELS,
                              ylabel='Fraction of HPV16 reads' if gene_idx == 0 else '',
                              adjusted_pvals=gene_qvals)

            ax.set_title(gene, fontsize=FONT_AXIS - 2, fontweight='bold', pad=15)
            gene_idx += 1

    # Phase headers
    fig.canvas.draw()
    phase_pairs_idx = [(0, 1), (2, 3), (4, 5), (6, 7)]
    for (left, right), (phase_name, _) in zip(phase_pairs_idx, HPV16_PHASES.items()):
        pos_l = f_axes[left].get_position()
        pos_r = f_axes[right].get_position()
        mid_x = (pos_l.x0 + pos_r.x1) / 2.0
        top_y = pos_l.y1 + 0.04
        fig.text(mid_x, top_y, phase_name,
                 fontsize=FONT_PHASE, fontweight='bold', ha='center', va='bottom',
                 color=PHASE_COLORS[phase_name])

    # Legend
    legend_elements = [
        Patch(facecolor=POP_COLORS[p], edgecolor='#333333', label=POP_LABELS[p])
        for p in POP_ORDER
    ]
    fig.legend(handles=legend_elements, loc='lower center',
               ncol=3, fontsize=FONT_LEGEND + 2, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.01))

    for ext in ['pdf', 'png']:
        out = os.path.join(OUTPUT_DIR, f"Panel_6F_permutation_violins.{ext}")
        fig.savefig(out, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    log(f"  [SAVE] Panel_6F_permutation_violins.pdf/.png")

    banner("COMPLETE")


if __name__ == "__main__":
    main()
