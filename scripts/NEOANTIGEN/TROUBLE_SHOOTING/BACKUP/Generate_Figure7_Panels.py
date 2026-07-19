#!/usr/bin/env python3
"""
Generate_Figure7_Panels.py
==========================
Figure 7 (MAIN): A3-driven neoantigen landscape of HPV16 lifecycle populations.

REDESIGN (replaces the ANXA1-centric v2):
  Panel A  Neoantigen burden (SBS2-HIGH vs CNV-HIGH, 1.77x) + depth-corrected
           RNA-fusion comparison (per-UMI, paired over shared samples -> n.s.)
  Panel B  Corrected neoantigen-gene Venn: 82 shared / 272 SBS2-specific /
           143 CNV-specific, labelled as the three tiers
  Panel C  KLF3 A3-lead card  (Tier 1, shared): TCW context + wt->mut binding + expression
  Panel D  CAST A3-lead card  (Tier 2, SBS2-specific): same layout

ANXA1 has moved to the supplemental figure as the off-signature counterpoint
(see Generate_Figure7_Supplemental.py).

Every number is read from disk via neoantigen_figure_utils (single source of
truth). Emits individual panels and a combined figure, PDF + PNG at 300 DPI.

Run in the NETWORK conda env.

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
import matplotlib.gridspec as gridspec
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

import neoantigen_figure_utils as U

plt.rcParams.update(U.BASE_RC)

OUTPUT_DIR = os.path.join(U.FIG7_ROOT, "FIGURE_7_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

LEAD_GENES = ['KLF3', 'CAST']

# Which color labels each burden series
COL_SBS2 = U.COL_SBS2
COL_CNV = U.COL_CNV


def _save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{ext}"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {name}.pdf / .png")


# =============================================================================
# PANEL A: burden (neoantigens) + depth-corrected fusion
# =============================================================================
def draw_panel_a_neo(ax, burden, letter='A'):
    labels = ['SBS2-HIGH', 'CNV-HIGH']
    totals = [burden['sbs2_total'], burden['cnv_total']]
    per_cell = [burden['sbs2_per_cell'], burden['cnv_per_cell']]
    strong = [burden['sbs2_strong'], burden['cnv_strong']]
    colors = [COL_SBS2, COL_CNV]

    x = np.arange(2)
    bars = ax.bar(x, totals, width=0.62, color=colors, edgecolor='white', linewidth=2)
    for i, bar in enumerate(bars):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + max(totals) * 0.015,
                f"{totals[i]:,}", ha='center', va='bottom',
                fontsize=U.FS_ANNOT, fontweight='bold', color=U.COL_INK)
        ax.text(bar.get_x() + bar.get_width() / 2, h * 0.5,
                f"{per_cell[i]:.2f}/cell\n{strong[i]} strong", ha='center', va='center',
                fontsize=U.FS_ANNOT - 2, color='white', fontweight='bold')

    # Fold-difference annotation
    ymax = max(totals)
    ax.annotate('', xy=(1, totals[1] + ymax * 0.10), xytext=(0, totals[0] + ymax * 0.10),
                arrowprops=dict(arrowstyle='<->', color='#555555', lw=2))
    ax.text(0.5, ymax + ymax * 0.14, f"{burden['fold']:.2f}x", ha='center', va='bottom',
            fontsize=U.FS_LABEL, fontstyle='italic', color='#555555')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=U.FS_TICK)
    ax.set_ylabel('predicted neoantigens', fontsize=U.FS_LABEL)
    ax.set_ylim(0, ymax * 1.30)
    ax.set_title(f"{letter}   Neoantigen burden", fontsize=U.FS_TITLE,
                 fontweight='bold', loc='left')


def draw_panel_a_fusion(ax, fusion):
    if fusion is None:
        ax.axis('off')
        ax.text(0.5, 0.5, "fusion metrics\nnot found", ha='center', va='center',
                fontsize=U.FS_LABEL, color='#999999', transform=ax.transAxes)
        return

    x = np.arange(2)
    meds = [fusion['sbs2_umi_med'], fusion['cnv_umi_med']]
    colors = [COL_SBS2, COL_CNV]
    ax.bar(x, meds, width=0.62, color=colors, edgecolor='white', linewidth=2,
           alpha=0.55, zorder=2)

    # Paired per-sample dots + connectors (the actual distribution behind the medians)
    s = np.asarray(fusion['sbs2_umi'], dtype=float)
    c = np.asarray(fusion['cnv_umi'], dtype=float)
    if len(s) and len(s) == len(c):
        jitter = (np.random.RandomState(0).rand(len(s)) - 0.5) * 0.12
        for i in range(len(s)):
            ax.plot([0 + jitter[i], 1 + jitter[i]], [s[i], c[i]],
                    color='#B0B0B0', lw=1.2, alpha=0.6, zorder=3)
        ax.scatter(np.zeros(len(s)) + jitter, s, s=90, color=COL_SBS2,
                   edgecolor=U.COL_INK, linewidth=1, zorder=4)
        ax.scatter(np.ones(len(c)) + jitter, c, s=90, color=COL_CNV,
                   edgecolor=U.COL_INK, linewidth=1, zorder=4)

    for i, mv in enumerate(meds):
        ax.text(i, mv + max(meds) * 0.04, f"{mv:.3f}", ha='center', va='bottom',
                fontsize=U.FS_ANNOT, fontweight='bold', color=U.COL_INK)

    ptxt = "n.s." if not np.isfinite(fusion['wilcoxon_p']) or fusion['wilcoxon_p'] >= 0.05 \
        else f"p={fusion['wilcoxon_p']:.3f}"
    ax.text(0.5, max(meds) * 1.28,
            f"comparable ({ptxt})\nfrac. chimeric {fusion['sbs2_frac_med']:.1f}% vs "
            f"{fusion['cnv_frac_med']:.1f}%",
            ha='center', va='bottom', fontsize=U.FS_ANNOT - 2, fontstyle='italic',
            color='#555555')

    ax.set_xticks(x)
    ax.set_xticklabels(['SBS2-HIGH', 'CNV-HIGH'], fontsize=U.FS_TICK)
    ax.set_ylabel('RNA fusions per UMI\n(depth-corrected)', fontsize=U.FS_LABEL)
    ax.set_ylim(0, max(meds) * 1.55)
    ax.set_title(f"RNA fusion burden (n={fusion['n_samples']} shared)",
                 fontsize=U.FS_LABEL, fontweight='bold', loc='left')


# =============================================================================
# PANEL B: corrected Venn (82 / 272 / 143), three tiers
# =============================================================================
def draw_panel_b(ax, venn, letter='B'):
    ax.set_xlim(-6.5, 6.5)
    ax.set_ylim(-5.5, 5.5)
    ax.set_aspect('equal')
    ax.axis('off')

    r, off = 3.0, 1.7
    ax.add_patch(Circle((-off, 0), r, facecolor=COL_SBS2, alpha=0.18,
                        edgecolor=COL_SBS2, linewidth=3))
    ax.add_patch(Circle((off, 0), r, facecolor=COL_CNV, alpha=0.30,
                        edgecolor='#C9A227', linewidth=3))

    # Counts
    ax.text(-off - 1.15, 0.4, str(venn['sbs2_only']), ha='center', va='center',
            fontsize=U.FS_TITLE + 6, fontweight='bold', color=COL_SBS2)
    ax.text(off + 1.15, 0.4, str(venn['cnv_only']), ha='center', va='center',
            fontsize=U.FS_TITLE + 6, fontweight='bold', color='#C9A227')
    ax.text(0, 0.4, str(venn['shared']), ha='center', va='center',
            fontsize=U.FS_TITLE + 6, fontweight='bold', color=U.COL_SHARED)

    # Tier labels under each count
    ax.text(-off - 1.15, -0.9, 'Tier 2\nSBS2-specific', ha='center', va='top',
            fontsize=U.FS_ANNOT, color=COL_SBS2)
    ax.text(off + 1.15, -0.9, 'Tier 3\nCNV-specific', ha='center', va='top',
            fontsize=U.FS_ANNOT, color='#C9A227')
    ax.text(0, -0.9, 'Tier 1\nshared', ha='center', va='top',
            fontsize=U.FS_ANNOT, color=U.COL_SHARED)

    # Group totals above circles
    ax.text(-off, r + 0.35, f"SBS2-HIGH ({venn['sbs2_total']})", ha='center', va='bottom',
            fontsize=U.FS_LABEL, fontweight='bold', color=COL_SBS2)
    ax.text(off, r + 0.35, f"CNV-HIGH ({venn['cnv_total']})", ha='center', va='bottom',
            fontsize=U.FS_LABEL, fontweight='bold', color='#C9A227')

    ax.set_title(f"{letter}   Neoantigen-gene overlap", fontsize=U.FS_TITLE,
                 fontweight='bold', loc='left')


# =============================================================================
# BUILD
# =============================================================================
def main():
    print("=" * 70)
    print("  FIGURE 7 MAIN  (KLF3 + CAST A3 leads)")
    print("=" * 70)
    burden = U.load_burden()
    venn = U.load_venn()
    fusion = U.load_fusion_depth_corrected()
    records = U.build_lead_records(LEAD_GENES)

    # ---- individual panels ----------------------------------------------
    fig, (axn, axf) = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={'width_ratios': [1, 1]})
    draw_panel_a_neo(axn, burden)
    draw_panel_a_fusion(axf, fusion)
    fig.tight_layout()
    _save(fig, "Figure7_PanelA_Burden")

    fig, ax = plt.subplots(figsize=(13, 11))
    draw_panel_b(ax, venn)
    fig.tight_layout()
    _save(fig, "Figure7_PanelB_Venn")

    letters = ['C', 'D']
    for rec, letter in zip(records, letters):
        fig = plt.figure(figsize=(26, 8))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.0, 1.7, 0.9], wspace=0.30)
        axc = fig.add_subplot(gs[0, 0])
        axb = fig.add_subplot(gs[0, 1])
        axe = fig.add_subplot(gs[0, 2])
        U.draw_lead_card(axc, axb, axe, rec, letter=letter)
        fig.subplots_adjust(top=0.78, bottom=0.22, left=0.04, right=0.98)
        _save(fig, f"Figure7_Panel{letter}_{rec['gene']}")

    # ---- combined figure -------------------------------------------------
    fig = plt.figure(figsize=(30, 30))
    outer = gridspec.GridSpec(3, 1, height_ratios=[1.05, 0.95, 0.95], hspace=0.55)

    top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0],
                                           width_ratios=[1.0, 0.85, 1.35], wspace=0.42)
    ax_neo = fig.add_subplot(top[0, 0])
    ax_fus = fig.add_subplot(top[0, 1])
    ax_venn = fig.add_subplot(top[0, 2])
    draw_panel_a_neo(ax_neo, burden)
    draw_panel_a_fusion(ax_fus, fusion)
    draw_panel_b(ax_venn, venn)

    for idx, (rec, letter) in enumerate(zip(records, letters)):
        row = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[idx + 1],
                                               width_ratios=[1.0, 1.7, 0.9], wspace=0.30)
        axc = fig.add_subplot(row[0, 0])
        axb = fig.add_subplot(row[0, 1])
        axe = fig.add_subplot(row[0, 2])
        U.draw_lead_card(axc, axb, axe, rec, letter=letter)

    _save(fig, "Figure7_Combined")

    print("\nAll main-figure panels written to:", OUTPUT_DIR)
    print("Done.")


if __name__ == "__main__":
    main()
