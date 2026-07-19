#!/usr/bin/env python3
"""
neoantigen_figure_utils.py
==========================
Shared loaders and renderers for the Figure 7 neoantigen landscape (main +
supplemental). Every number is read from a locked diagnostic output; nothing is
recomputed here.

Card layout (one gene): three coupled panels
  - MHC-binding bars  : one row per differential neoantigen mutation, best peptide,
                        wild-type -> mutant IC50 drop (dumbbell), sorted by mut IC50.
  - Gene track        : native-HGVS-frame backbone, domain boxes, disordered underlay,
                        lollipops colored by TCW status at their labeled positions.
  - Expression/carrier: tier-correct % expressing + fraction of expressing cells
                        carrying any of the gene's mutations.

Data sources (all under data/FIG_7/):
  03_mhc_binding/{group}_neoantigens.tsv                 (binding, differential)
  fasta_context/ref_tri_fasta.tsv                        (genome TCW per variant)
  protein_domains/protein_domains.tsv (+ {gene}_mutation_map.tsv)   (native tracks)
  TROUBLESHOOTING/group_aware_expression/
      group_aware_expression_carrier.tsv                 (expression + carrier)
      panelA_burden_stats.tsv                            (Panel A means/SEM/adj p)

Conventions: Agg backend, hex colors, fonts 28-34, PDF+PNG at 300 DPI.
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import re
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.gridspec import GridSpecFromSubplotSpec

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
MHC_DIR = os.path.join(BASE_DIR, "data/FIG_7/03_mhc_binding")
FASTA_TSV = os.path.join(BASE_DIR, "data/FIG_7/fasta_context/ref_tri_fasta.tsv")
DOMAINS_TSV = os.path.join(BASE_DIR, "data/FIG_7/protein_domains/protein_domains.tsv")
DOMAINS_DIR = os.path.join(BASE_DIR, "data/FIG_7/protein_domains")
GA_DIR = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/group_aware_expression")
EXPR_TSV = os.path.join(GA_DIR, "group_aware_expression_carrier.tsv")
PANELA_TSV = os.path.join(GA_DIR, "panelA_burden_stats.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/figures")

# =============================================================================
# STYLE
# =============================================================================
COLOR_SBS2 = '#ed6a5a'     # coral  (SBS2 / A3A)
COLOR_CNV = '#F6D155'      # mustard (CNV / A3B)
COLOR_SHARED = '#9B59B6'   # purple (shared tier)
TCW_CT = '#ed6a5a'         # clean C>T TCW  (SBS2 signature)
TCW_CG = '#E67E22'         # C>G TCW        (SBS13, off-signature)
NONTCW = '#9AA0A6'         # non-APOBEC     (gray)
DISORDER = '#B0BEC5'       # disordered regions
BACKBONE = '#D5DBDB'       # protein backbone bar
THRESH_IC50 = 500.0        # binder threshold (nM)

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT, FS_SMALL = 34, 30, 28, 28, 24

GENE_TIER = {'KLF3': 'Tier1_shared', 'CAST': 'Tier2_sbs2_specific',
             'SERPINB2': 'Tier1_shared', 'KRT6B': 'Tier2_sbs2_specific',
             'ANXA1': 'Tier2_sbs2_specific'}
GENE_GROUP = {g: 'SBS2_HIGH' for g in GENE_TIER}   # neoantigen source group

plt.rcParams.update({'font.size': FS_TICK, 'axes.linewidth': 1.4,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})


# =============================================================================
# HELPERS
# =============================================================================
def _pos(hgvs):
    m = re.match(r'p\.[A-Za-z]{3}(\d+)', str(hgvs))
    return int(m.group(1)) if m else None


def tcw_color(is_tcw, is_tcw_ct):
    if bool(is_tcw_ct):
        return TCW_CT
    if bool(is_tcw):
        return TCW_CG
    return NONTCW


def _read(path, what):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{what} not found: {path}")
    return pd.read_csv(path, sep='\t')


def _tobool(s):
    return str(s).strip().lower() in ('true', '1', '1.0', 'yes')


# =============================================================================
# LOADERS
# =============================================================================
def load_panelA():
    """means, SEM, and BH-adjusted p for the two per-cell per-UMI burden tests."""
    df = _read(PANELA_TSV, "Panel A burden stats").set_index('comparison')
    out = {}
    for key, label in [('neoantigen_perUMI', 'neoantigen'), ('fusion_perUMI', 'fusion')]:
        if key in df.index:
            r = df.loc[key]
            out[label] = {'sbs2': float(r['sbs2_mean']), 'cnv': float(r['cnv_mean']),
                          'sbs2_sem': float(r.get('sbs2_sem', np.nan)),
                          'cnv_sem': float(r.get('cnv_sem', np.nan)),
                          'adj_p': float(r['bh_adjusted_p']),
                          'metric': str(r.get('metric', ''))}
    return out


def load_venn():
    """Neoantigen gene sets -> (sbs2_only, shared, cnv_only)."""
    def genes(group):
        p = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
        return set(_read(p, f"{group} neoantigens")['gene'].dropna().astype(str))
    s, c = genes('SBS2_HIGH'), genes('CNV_HIGH')
    shared = s & c
    return len(s - c), len(shared), len(c - s)


def load_fasta_tcw():
    df = _read(FASTA_TSV, "FASTA trinucleotide table")
    df['_pos'] = df['hgvs_p'].map(_pos)
    return df


def load_gene_binding(gene, fasta=None):
    """Differential neoantigen mutations for a gene: best peptide per mutation,
    with TCW status joined from the FASTA table. Sorted by mut IC50 ascending."""
    group = GENE_GROUP.get(gene, 'SBS2_HIGH')
    neo = _read(os.path.join(MHC_DIR, f"{group}_neoantigens.tsv"), f"{group} neoantigens")
    sub = neo[neo['gene'].astype(str) == gene].copy()
    if 'is_differential' in sub.columns:
        sub = sub[sub['is_differential'].map(_tobool)]
    if len(sub) == 0:
        return pd.DataFrame(columns=['hgvs_p', 'pos', 'mut_peptide', 'wt_ic50',
                                     'mut_ic50', 'best_allele', 'is_tcw', 'is_tcw_ct'])
    sub = sub.sort_values('mut_ic50').drop_duplicates(subset=['hgvs_p'], keep='first').copy()
    sub['pos'] = sub['hgvs_p'].map(_pos)
    if fasta is None:
        fasta = load_fasta_tcw()
    key = fasta[fasta['gene'].astype(str) == gene][['hgvs_p', 'is_tcw', 'is_tcw_ct']].copy()
    key['is_tcw'] = key['is_tcw'].map(_tobool)
    key['is_tcw_ct'] = key['is_tcw_ct'].map(_tobool)
    sub = sub.merge(key, on='hgvs_p', how='left')
    sub['is_tcw'] = sub['is_tcw'].fillna(False)
    sub['is_tcw_ct'] = sub['is_tcw_ct'].fillna(False)
    return sub.sort_values('mut_ic50').reset_index(drop=True)


def load_gene_track(gene):
    """Native-frame domain/disordered segments + backbone length for a gene."""
    dm = _read(DOMAINS_TSV, "protein domains")
    g = dm[dm['gene'].astype(str) == gene].copy()
    native_len = int(g['native_length'].iloc[0]) if len(g) else 0
    acc = str(g['uniprot_acc'].iloc[0]) if len(g) else ''
    domains = g[g['category'] == 'domain'][['start', 'end', 'name', 'color']]
    disordered = g[g['category'] == 'disordered'][['start', 'end', 'name']]
    return {'native_len': native_len, 'acc': acc,
            'domains': domains.to_dict('records'),
            'disordered': disordered.to_dict('records')}


def load_expression_carrier(gene):
    df = _read(EXPR_TSV, "group-aware expression/carrier")
    r = df[df['gene'].astype(str) == gene]
    if len(r) == 0:
        return None
    r = r.iloc[0]
    return {'tier': str(r['tier']), 'headline_pct': float(r['headline_pct']),
            'pct_sbs2': float(r['pct_sbs2']), 'pct_cnv': float(r['pct_cnv']),
            'pct_pooled': float(r['pct_pooled']),
            'carrier_pct': float(r['carrier_pct_all_expressers']),
            'n_expressers': int(r['n_expressers']), 'n_carriers': int(r['n_carriers'])}


def load_gene_card(gene, fasta=None):
    b = load_gene_binding(gene, fasta=fasta)
    t = load_gene_track(gene)
    e = load_expression_carrier(gene)
    lolli = [{'pos': int(row.pos), 'mut_ic50': float(row.mut_ic50),
              'color': tcw_color(row.is_tcw, row.is_tcw_ct), 'hgvs': row.hgvs_p}
             for row in b.itertuples() if pd.notna(row.pos)]
    return {'gene': gene, 'tier': GENE_TIER.get(gene, '?'), 'binding': b,
            'track': t, 'expr': e, 'lollipops': lolli}


# =============================================================================
# RENDERERS
# =============================================================================
def draw_binding_bars(ax, gene, binding):
    ax.set_title(f"{gene}  MHC-I binding", fontsize=FS_LABEL, style='italic', pad=10)
    if len(binding) == 0:
        ax.text(0.5, 0.5, "no differential neoantigen", ha='center', va='center',
                fontsize=FS_SMALL, transform=ax.transAxes)
        ax.axis('off')
        return
    n = len(binding)
    ys = np.arange(n)[::-1]
    for y, row in zip(ys, binding.itertuples()):
        wt = min(float(row.wt_ic50), 50000.0)
        mut = float(row.mut_ic50)
        col = tcw_color(row.is_tcw, row.is_tcw_ct)
        ax.plot([mut, wt], [y, y], '-', color='#B0B0B0', lw=3, zorder=1)
        ax.scatter([wt], [y], s=180, facecolors='white', edgecolors='#555555',
                   linewidths=2.2, zorder=2)
        ax.scatter([mut], [y], s=260, facecolors=col, edgecolors='#333333',
                   linewidths=2.2, zorder=3)
        ax.text(mut, y + 0.30, f"{row.mut_peptide}", ha='center', va='bottom',
                fontsize=FS_SMALL - 4, family='monospace')
        ax.text(-0.02, y, f"{row.hgvs_p}", transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=FS_SMALL - 2, clip_on=False)
    ax.axvline(THRESH_IC50, color='#444444', ls='--', lw=2)
    ax.text(THRESH_IC50, ys.max() + 0.7, '500 nM', ha='center', va='bottom',
            fontsize=FS_SMALL - 2, color='#444444')
    ax.set_xscale('log')
    ax.set_xlim(1, 55000)
    ax.set_ylim(-0.7, n - 0.3 + 0.9)
    ax.set_yticks([])
    ax.set_xlabel('IC50 (nM, log)   wt (open) → mut (filled)', fontsize=FS_TICK - 2)
    ax.tick_params(labelsize=FS_TICK - 2)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)


def draw_gene_track(ax, gene, track, lollipops):
    L = max(track['native_len'], 1)
    ax.set_xlim(-L * 0.02, L * 1.02)
    ax.set_ylim(-1.2, 2.9)
    bar_y, bar_h = 0.0, 0.5
    for d in track['disordered']:
        ax.add_patch(Rectangle((d['start'], bar_y - 0.18), d['end'] - d['start'],
                               bar_h + 0.36, facecolor=DISORDER, edgecolor='none',
                               alpha=0.55, zorder=1))
    ax.add_patch(Rectangle((0, bar_y), L, bar_h, facecolor=BACKBONE,
                           edgecolor='#7F8C8D', linewidth=1.2, zorder=2))
    for dom in track['domains']:
        ax.add_patch(Rectangle((dom['start'], bar_y - 0.06), dom['end'] - dom['start'],
                               bar_h + 0.12, facecolor=dom['color'], edgecolor='#333333',
                               linewidth=1.4, zorder=3))
        ax.text((dom['start'] + dom['end']) / 2, bar_y + bar_h / 2, dom['name'],
                ha='center', va='center', fontsize=FS_SMALL - 8, zorder=4, clip_on=True)
    if lollipops:
        anchor = min(range(len(lollipops)), key=lambda i: lollipops[i]['mut_ic50'])
        for i, lp in enumerate(lollipops):
            x = lp['pos']
            ax.plot([x, x], [bar_y + bar_h, bar_y + bar_h + 1.3], '-',
                    color='#555555', lw=2, zorder=5)
            ax.scatter([x], [bar_y + bar_h + 1.4], s=520 if i == anchor else 320,
                       facecolors=lp['color'], edgecolors='#222222', linewidths=2.4,
                       zorder=6)
    ax.text(0, bar_y + bar_h + 2.15, f"{gene}  ({track['acc']}, {L} aa)",
            fontsize=FS_SMALL, style='italic', va='bottom')
    ax.set_yticks([])
    ax.set_xlabel('residue (native numbering)', fontsize=FS_TICK - 4)
    ax.tick_params(labelsize=FS_TICK - 4)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)


def draw_expression_carrier(ax, gene, expr):
    if expr is None:
        ax.axis('off')
        return
    shared = expr['tier'].startswith('Tier1')
    denom_lbl = 'SBS2+CNV' if shared else 'SBS2'
    expr_col = COLOR_SHARED if shared else COLOR_SBS2
    ax.barh([1.0], [expr['headline_pct']], height=0.55, color=expr_col,
            edgecolor='#333333', linewidth=1.4)
    ax.text(min(expr['headline_pct'] + 1.5, 101), 1.0, f"{expr['headline_pct']:.1f}%",
            va='center', fontsize=FS_SMALL)
    ax.text(0, 1.55, f"% expressing ({denom_lbl})", fontsize=FS_SMALL - 2, va='bottom')
    ax.barh([0.0], [expr['carrier_pct']], height=0.55, color='#34495E',
            edgecolor='#333333', linewidth=1.4)
    ax.text(expr['carrier_pct'] + 1.5, 0.0,
            f"{expr['carrier_pct']:.1f}%  ({expr['n_carriers']}/{expr['n_expressers']})",
            va='center', fontsize=FS_SMALL)
    ax.text(0, 0.55, "carrier / expressers", fontsize=FS_SMALL - 2, va='bottom')
    ax.set_xlim(0, 108)
    ax.set_ylim(-0.6, 2.1)
    ax.set_yticks([])
    ax.set_xlabel('percent', fontsize=FS_TICK - 4)
    ax.tick_params(labelsize=FS_TICK - 4)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)


def make_card_axes(fig, cell):
    gs = GridSpecFromSubplotSpec(2, 2, subplot_spec=cell, width_ratios=[0.44, 0.56],
                                 height_ratios=[0.62, 0.38], wspace=0.30, hspace=0.6)
    ax_bind = fig.add_subplot(gs[:, 0])
    ax_track = fig.add_subplot(gs[0, 1])
    ax_expr = fig.add_subplot(gs[1, 1])
    return ax_bind, ax_track, ax_expr


def draw_card(fig, cell, card):
    ax_bind, ax_track, ax_expr = make_card_axes(fig, cell)
    draw_binding_bars(ax_bind, card['gene'], card['binding'])
    draw_gene_track(ax_track, card['gene'], card['track'], card['lollipops'])
    draw_expression_carrier(ax_expr, card['gene'], card['expr'])


def draw_panelA(ax, stats):
    blocks = [('neoantigen', 'neoantigens'), ('fusion', 'RNA fusions')]
    xcenters = [0, 1.5]
    w = 0.5
    for (key, _), xc in zip(blocks, xcenters):
        d = stats.get(key)
        if not d:
            continue
        vals = [d['sbs2'], d['cnv']]
        sems = [d.get('sbs2_sem', np.nan), d.get('cnv_sem', np.nan)]
        xs = [xc - w / 2, xc + w / 2]
        ax.bar(xs, vals, width=w, color=[COLOR_SBS2, COLOR_CNV],
               edgecolor='#333333', linewidth=1.6,
               yerr=[[0, 0], [s if np.isfinite(s) else 0 for s in sems]],
               error_kw={'elinewidth': 2, 'capsize': 6, 'capthick': 2})
        top = max(v + (s if np.isfinite(s) else 0) for v, s in zip(vals, sems))
        br = top * 1.12
        ax.plot([xs[0], xs[0], xs[1], xs[1]], [top * 1.04, br, br, top * 1.04],
                color='#333333', lw=1.8)
        p = d['adj_p']
        ptxt = f"adj p = {p:.2g}" if p >= 1e-4 else f"adj p = {p:.1e}"
        ax.text(xc, br * 1.02, ptxt, ha='center', va='bottom', fontsize=FS_SMALL)
    ax.set_xticks(xcenters)
    ax.set_xticklabels(['neoantigens', 'RNA fusions'], fontsize=FS_TICK)
    ax.set_ylabel('events / 1000 UMI (per cell)', fontsize=FS_TICK)
    ax.tick_params(labelsize=FS_TICK - 2)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=COLOR_SBS2, edgecolor='#333', label='SBS2-HIGH'),
                       Patch(facecolor=COLOR_CNV, edgecolor='#333', label='CNV-HIGH')],
              fontsize=FS_SMALL, frameon=False, loc='upper right')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title('A  mutational burden per cell', fontsize=FS_LABEL, loc='left', pad=10)


def draw_panelB_venn(ax, counts):
    sbs2_only, shared, cnv_only = counts
    ax.set_xlim(0, 10); ax.set_ylim(0, 7); ax.axis('off')
    ax.add_patch(Circle((4.0, 3.4), 2.6, facecolor=COLOR_SBS2, edgecolor='#333',
                        linewidth=1.8, alpha=0.55))
    ax.add_patch(Circle((6.0, 3.4), 2.6, facecolor=COLOR_CNV, edgecolor='#333',
                        linewidth=1.8, alpha=0.55))
    ax.text(2.8, 3.4, f"{sbs2_only}", ha='center', va='center', fontsize=FS_TITLE)
    ax.text(7.2, 3.4, f"{cnv_only}", ha='center', va='center', fontsize=FS_TITLE)
    ax.text(5.0, 3.4, f"{shared}", ha='center', va='center', fontsize=FS_LABEL)
    ax.text(2.8, 0.5, 'SBS2-specific', ha='center', fontsize=FS_SMALL, color='#C0392B')
    ax.text(7.2, 0.5, 'CNV-specific', ha='center', fontsize=FS_SMALL, color='#B7950B')
    ax.text(5.0, 6.4, 'neoantigen genes', ha='center', fontsize=FS_SMALL)
    ax.set_title('B  neoantigen gene overlap', fontsize=FS_LABEL, loc='left')


def draw_carrier_comparison(ax, genes=('KLF3', 'CAST', 'SERPINB2', 'KRT6B', 'ANXA1')):
    """Carrier fraction across the five genes, highlighting the SERPINB2/KRT6B
    ~10x recurrence and flagging ANXA1 as the off-signature counterpoint."""
    data = []
    for g in genes:
        e = load_expression_carrier(g)
        if e:
            data.append((g, e['carrier_pct']))
    names = [x[0] for x in data]
    vals = [x[1] for x in data]
    colors = []
    for g, _ in data:
        if g in ('SERPINB2', 'KRT6B'):
            colors.append(COLOR_SBS2)          # the recurrent carriers
        elif g == 'ANXA1':
            colors.append(TCW_CG)              # off-signature (SBS13)
        else:
            colors.append('#BDC3C7')
    xs = np.arange(len(names))
    ax.bar(xs, vals, width=0.66, color=colors, edgecolor='#333333', linewidth=1.4)
    for x, v in zip(xs, vals):
        ax.text(x, v + 0.3, f"{v:.1f}%", ha='center', va='bottom', fontsize=FS_SMALL - 2)
    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=FS_TICK - 2, style='italic')
    ax.set_ylabel('carrier / expressers (%)', fontsize=FS_TICK - 2)
    ax.tick_params(labelsize=FS_TICK - 4)
    ax.set_ylim(0, max(vals) * 1.25 if vals else 1)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=COLOR_SBS2, edgecolor='#333', label='recurrent carrier'),
                       Patch(facecolor=TCW_CG, edgecolor='#333', label='off-signature (SBS13)'),
                       Patch(facecolor='#BDC3C7', edgecolor='#333', label='low carrier')],
              fontsize=FS_SMALL - 2, frameon=False, loc='upper left')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title('A  neoantigen carrier fraction among expressing cells',
                 fontsize=FS_LABEL, loc='left', pad=10)


def save_figure(fig, name):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{ext}"), dpi=300,
                    bbox_inches='tight', facecolor='white')
    return os.path.join(OUTPUT_DIR, f"{name}.pdf")
