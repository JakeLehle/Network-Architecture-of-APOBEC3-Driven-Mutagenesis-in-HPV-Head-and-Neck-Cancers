#!/usr/bin/env python3
"""
neoantigen_figure_utils.py  (v2, single-source)
===============================================
Shared loaders and renderers for the Figure 7 neoantigen landscape. Every
binding, prevalence, and expression number is read from ONE locked table,
neoantigen_prevalence_ranking_full.tsv, so a panel value and the manuscript
value cannot drift. Nothing is recomputed here.

Retargeted to three featured vaccine candidates, one per tier:
    COX4I1 p.Ala9Thr   CNV-specific (Tier 3)   prevalence 24.9%
    SPRR1A p.Val61Ile  shared       (Tier 1)   prevalence 16.9% (pooled /1092)
    KRT6B  p.Glu342Lys SBS2-specific(Tier 2)   prevalence  7.9%, TCW clean C>T

Two color languages, one per analytical axis:
    TIER palette  (prevalence / expression panels): CNV mustard, shared purple,
                  SBS2 coral. Encodes which niche(s) a hit belongs to.
    TCW palette   (binding / track panels): clean C>T coral, C>G orange,
                  non-APOBEC gray. Encodes mutational provenance.
Both legends are drawn so the dual encoding is unambiguous.

Building blocks:
  load_featured            the three featured rows from full.tsv (validated)
  load_gene_track          native-frame backbone + domains (protein_domains.tsv)
  load_panelA / load_venn  burden stats and gene-level Venn (own sources)
  draw_binding_featured    one panel, 3 horizontal stacked bars (wt + TCW gain)
  draw_expr_featured       one panel, 3 vertical stacked bars (expr + carriage)
  draw_gene_track          protein track with the featured lollipop
  draw_panelA / draw_panelB_venn   burden and Venn

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
from matplotlib.patches import Rectangle, Circle, Patch
from matplotlib.lines import Line2D
from matplotlib.offsetbox import TextArea, HPacker, AnnotationBbox

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
RANK_DIR = os.path.join(BASE_DIR, "data/FIG_7/06_prevalence_ranking")
FULL_TSV = os.path.join(RANK_DIR, "neoantigen_prevalence_ranking_full.tsv")   # SINGLE SOURCE
MHC_DIR = os.path.join(BASE_DIR, "data/FIG_7/03_mhc_binding")                 # Venn only
DOMAINS_TSV = os.path.join(BASE_DIR, "data/FIG_7/protein_domains/protein_domains.tsv")
GA_DIR = os.path.join(BASE_DIR, "data/FIG_7/TROUBLESHOOTING/group_aware_expression")
PANELA_TSV = os.path.join(GA_DIR, "panelA_burden_stats.tsv")                  # burden only
OUTPUT_DIR = os.path.join(BASE_DIR, "data/FIG_7/figures")

# =============================================================================
# FEATURED MUTATIONS  (order = prevalence_tier descending, used across bar panels)
# =============================================================================
FEATURED = {
    "COX4I1": "p.Ala9Thr",     # CNV-specific  (Tier 3)
    "SPRR1A": "p.Val61Ile",    # shared        (Tier 1)
    "KRT6B":  "p.Glu342Lys",   # SBS2-specific (Tier 2), TCW clean C>T
}
TRACK_LONG = "KRT6B"                       # full-width track (has real domains)
TRACK_SHORT = ["COX4I1", "SPRR1A"]         # paired short tracks

# =============================================================================
# STYLE
# =============================================================================
COLOR_SBS2 = '#ed6a5a'     # coral  (SBS2 tier / clean C>T TCW)
COLOR_CNV = '#F6D155'      # mustard (CNV tier)
COLOR_SHARED = '#9B59B6'   # purple (shared tier)
TCW_CT = '#ed6a5a'         # clean C>T TCW  (SBS2 signature)
TCW_CG = '#E67E22'         # C>G TCW        (SBS13, off-signature)
NONTCW = '#9AA0A6'         # non-APOBEC     (gray)
DISORDER = '#B0BEC5'
BACKBONE = '#D5DBDB'
WT_GRAY = '#D0D3D4'        # wild-type baseline / expressing (non-carrier) segment
MUT_RED = '#D1352B'        # mutated residue highlight in peptide sequences
THRESH_IC50 = 500.0
THRESH_SCORE = 1.0 - np.log10(THRESH_IC50) / np.log10(50000.0)   # ~0.426

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT, FS_SMALL = 34, 30, 28, 28, 24

TIER_COLOR = {'SBS2_specific': COLOR_SBS2, 'CNV_specific': COLOR_CNV,
              'shared': COLOR_SHARED}
TIER_LABEL = {'SBS2_specific': 'SBS2-specific', 'CNV_specific': 'CNV-specific',
              'shared': 'shared'}

plt.rcParams.update({'font.size': FS_TICK, 'axes.linewidth': 1.4,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})


# =============================================================================
# HELPERS
# =============================================================================
def _pos(hgvs):
    m = re.match(r'p\.[A-Za-z]{3}(\d+)', str(hgvs))
    return int(m.group(1)) if m else None


def _tobool(s):
    return str(s).strip().lower() in ('true', '1', '1.0', 'yes')


def tcw_color(is_tcw, is_tcw_ct):
    if _tobool(is_tcw_ct):
        return TCW_CT
    if _tobool(is_tcw):
        return TCW_CG
    return NONTCW


def _bind_score(ic50):
    """MHCflurry-style transformed affinity: 1 = strong (1 nM), 0 = weak (50000 nM)."""
    ic50 = float(min(max(float(ic50), 1.0), 50000.0))
    return float(np.clip(1.0 - np.log10(ic50) / np.log10(50000.0), 0.0, 1.0))


def _read(path, what):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{what} not found: {path}")
    return pd.read_csv(path, sep='\t')


def _mut_position_map():
    """(gene, hgvs_p, mut_peptide) -> 1-based index of the mutated residue in the
    peptide, read from the neoantigen files (mut_position_in_peptide). Position is
    a property of the peptide window, not the allele, so it is unambiguous. Used
    only to color the changed residue; it is not a manuscript number."""
    pos = {}
    for group in ('SBS2_HIGH', 'CNV_HIGH'):
        p = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p, sep='\t')
        if not {'gene', 'hgvs_p', 'mut_peptide', 'mut_position_in_peptide'} <= set(df.columns):
            continue
        d = df.dropna(subset=['mut_position_in_peptide'])
        for r in d.itertuples(index=False):
            pos[(str(r.gene), str(r.hgvs_p), str(r.mut_peptide))] = int(r.mut_position_in_peptide)
    return pos


def _colored_peptide(ax, x, y, peptide, mut_pos, suffix, fontsize):
    """Render a monospace peptide at data coords (x, y), left-anchored, with the
    residue at 1-based mut_pos drawn in red. `suffix` (e.g. IC50) is appended plain.
    Falls back to a single black label if mut_pos is missing/out of range."""
    tp = dict(family='monospace', size=fontsize)
    parts = []
    if mut_pos and 1 <= mut_pos <= len(peptide):
        pre, mid, post = peptide[:mut_pos - 1], peptide[mut_pos - 1], peptide[mut_pos:]
        if pre:
            parts.append(TextArea(pre, textprops={**tp, 'color': '#222222'}))
        parts.append(TextArea(mid, textprops={**tp, 'color': MUT_RED, 'weight': 'bold'}))
        if post:
            parts.append(TextArea(post, textprops={**tp, 'color': '#222222'}))
    else:
        parts.append(TextArea(peptide, textprops={**tp, 'color': '#222222'}))
    parts.append(TextArea(suffix, textprops={**tp, 'color': '#222222'}))
    box = HPacker(children=parts, align='baseline', pad=0, sep=0)
    ax.add_artist(AnnotationBbox(box, (x, y), xycoords=('data', 'data'),
                                 box_alignment=(0, 0.5), frameon=False, pad=0))


def binding_legend_handles():
    return [Patch(facecolor=WT_GRAY, edgecolor='#333', label='wild-type baseline'),
            Patch(facecolor=TCW_CT, edgecolor='#333', label='SBS2 gain (TCW C>T)'),
            Patch(facecolor=TCW_CG, edgecolor='#333', label='SBS13 gain (TCW C>G)'),
            Patch(facecolor=NONTCW, edgecolor='#333', label='non-APOBEC gain')]


def tier_legend_handles():
    return [Patch(facecolor=WT_GRAY, edgecolor='#333', label='expressing (no mutation)'),
            Patch(facecolor=COLOR_CNV, edgecolor='#333', label='carries mutation (CNV tier)'),
            Patch(facecolor=COLOR_SHARED, edgecolor='#333', label='carries mutation (shared tier)'),
            Patch(facecolor=COLOR_SBS2, edgecolor='#333', label='carries mutation (SBS2 tier)')]


# =============================================================================
# LOADERS  (single source: full.tsv)
# =============================================================================
def load_featured():
    """The three featured rows from full.tsv, validated, ordered by prevalence."""
    df = _read(FULL_TSV, "prevalence ranking full")
    rows = []
    for gene, hgvs in FEATURED.items():
        sub = df[(df['gene'].astype(str) == gene) & (df['hgvs_p'].astype(str) == hgvs)]
        if len(sub) == 0:
            near = df[df['gene'].astype(str) == gene]['hgvs_p'].tolist()[:6]
            raise ValueError(f"Featured mutation {gene} {hgvs} not in {FULL_TSV}. "
                             f"Re-run the ranking diagnostic (v3+). Gene rows present: {near}")
        rows.append(sub.iloc[0])
    out = pd.DataFrame(rows).reset_index(drop=True)
    out = out.sort_values('prevalence_tier', ascending=False).reset_index(drop=True)
    # attach the mutated-residue index within each featured peptide (for red highlight)
    pos = _mut_position_map()
    out['mut_pos_in_pep'] = [pos.get((str(r.gene), str(r.hgvs_p), str(r.peptide)))
                             for r in out.itertuples()]
    for r in out.itertuples():
        mp = r.mut_pos_in_pep
        pep = str(r.peptide)
        res = pep[mp - 1] if (mp and 1 <= mp <= len(pep)) else '?'
        print(f"  featured: {r.gene:8s} {r.hgvs_p:14s} {pep:12s} "
              f"mutated residue '{res}' at position {mp}")
    return out


def load_gene_track(gene):
    dm = _read(DOMAINS_TSV, "protein domains")
    g = dm[dm['gene'].astype(str) == gene].copy()
    if len(g) == 0:
        raise ValueError(
            f"{gene} is not in {DOMAINS_TSV}. Re-run Diagnostic_Fetch_Protein_Domains.py "
            f"with FIGURE_GENES = ['COX4I1', 'SPRR1A', 'KRT6B'] before building the figure.")
    native_len = int(g['native_length'].iloc[0])
    acc = str(g['uniprot_acc'].iloc[0])
    domains = g[g['category'] == 'domain'][['start', 'end', 'name', 'color']]
    disordered = g[g['category'] == 'disordered'][['start', 'end', 'name']]
    return {'native_len': native_len, 'acc': acc,
            'domains': domains.to_dict('records'),
            'disordered': disordered.to_dict('records')}


def featured_lollipop(feat_row):
    """Single featured lollipop for a gene's track, colored by TCW class."""
    p = _pos(feat_row['hgvs_p'])
    if p is None:
        return []
    return [{'pos': p, 'mut_ic50': float(feat_row['mut_IC50']),
             'color': tcw_color(feat_row['is_tcw'], feat_row['is_tcw_ct']),
             'hgvs': str(feat_row['hgvs_p'])}]


def load_panelA():
    df = _read(PANELA_TSV, "Panel A burden stats").set_index('comparison')
    out = {}
    for key, label in [('neoantigen_perUMI', 'neoantigen'), ('fusion_perUMI', 'fusion')]:
        if key in df.index:
            r = df.loc[key]
            out[label] = {'sbs2': float(r['sbs2_mean']), 'cnv': float(r['cnv_mean']),
                          'sbs2_sem': float(r.get('sbs2_sem', np.nan)),
                          'cnv_sem': float(r.get('cnv_sem', np.nan)),
                          'adj_p': float(r['bh_adjusted_p'])}
    return out


def load_venn():
    """Gene-level neoantigen overlap (should read the current 272/82/143)."""
    def genes(group):
        p = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
        return set(_read(p, f"{group} neoantigens")['gene'].dropna().astype(str))
    s, c = genes('SBS2_HIGH'), genes('CNV_HIGH')
    shared = s & c
    return len(s - c), len(shared), len(c - s)


# =============================================================================
# RENDERERS
# =============================================================================
def draw_binding_featured(ax, feat, panel=None):
    """One panel, one horizontal stacked bar per featured mutation: wild-type
    baseline binding score, then the mutation's binding gain stacked on top and
    colored by TCW provenance. Bars ordered by prevalence_tier (panel-consistent)."""
    ttl = "MHC-I binding gain of featured neoantigens"
    ax.set_title((f"{panel}  " if panel else "") + ttl, fontsize=FS_LABEL, loc='left', pad=10)
    n = len(feat)
    ys = np.arange(n)[::-1]
    for y, row in zip(ys, feat.itertuples()):
        s_wt = _bind_score(row.wt_IC50)
        s_mut = _bind_score(row.mut_IC50)
        col = tcw_color(row.is_tcw, row.is_tcw_ct)
        ax.barh(y, s_wt, height=0.58, color=WT_GRAY, edgecolor='#333333', linewidth=1.3, zorder=2)
        ax.barh(y, max(s_mut - s_wt, 0), left=s_wt, height=0.58, color=col,
                edgecolor='#333333', linewidth=1.3, zorder=3)
        ax.text(-0.02, y, f"{row.gene}  {row.hgvs_p}", transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=FS_SMALL - 2, clip_on=False)
        mut_pos = getattr(row, 'mut_pos_in_pep', None)
        _colored_peptide(ax, s_mut + 0.015, y, str(row.peptide), mut_pos,
                         f"  ({row.mut_IC50:.0f} nM)", FS_SMALL - 4)
    ax.axvline(THRESH_SCORE, color='#444444', ls='--', lw=2, zorder=4)
    ax.text(THRESH_SCORE, n - 0.35, 'binder\n(500 nM)', ha='center', va='bottom',
            fontsize=FS_SMALL - 6, color='#444444')
    ax.set_xlim(0, 1.32)
    ax.set_ylim(-0.7, n - 0.4 + 0.7)
    ax.set_yticks([])
    ax.set_xlabel('MHC-I binding score  (1 = strong)', fontsize=FS_TICK - 2)
    ax.tick_params(labelsize=FS_TICK - 4)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    # inline TCW legend (this panel and the tracks are colored by provenance)
    ax.legend(handles=binding_legend_handles(), fontsize=FS_SMALL - 6, frameon=False,
              loc='lower right', bbox_to_anchor=(1.0, -0.04), ncol=1)


def draw_expr_featured(ax, feat, panel=None):
    """One panel, one vertical stacked bar per featured mutation. Bar height =
    % of TIER cells expressing the gene (denominator 546 for specific tiers, 1092
    for shared). The colored base = prevalence_tier (the fraction that both
    express and carry the mutation), colored by tier. The colored base IS the
    headline prevalence number."""
    xs = np.arange(len(feat))
    for x, row in zip(xs, feat.itertuples()):
        total = float(row.pct_expressing_tier)
        carriage = float(row.prevalence_tier) * 100.0     # == coral segment height
        col = TIER_COLOR.get(row.tier, NONTCW)
        ax.bar(x, carriage, width=0.60, color=col, edgecolor='#333333', linewidth=1.4, zorder=3)
        ax.bar(x, max(total - carriage, 0), bottom=carriage, width=0.60, color=WT_GRAY,
               edgecolor='#333333', linewidth=1.4, zorder=2)
        ax.text(x, total + 1.6, f"{total:.0f}% expr", ha='center', va='bottom', fontsize=FS_SMALL - 2)
        ax.text(x, carriage + 0.6, f"{carriage:.1f}%", ha='center', va='bottom',
                fontsize=FS_SMALL - 2, color='#1a1a1a', fontweight='bold')
    ax.set_xticks(xs)
    ax.set_xticklabels([r.gene for r in feat.itertuples()], fontsize=FS_TICK - 2)
    # tier sublabel under each gene, colored by tier so the palette is self-documenting
    ymax = max(float(r.pct_expressing_tier) for r in feat.itertuples())
    ax.set_ylim(0, min(ymax * 1.25, 108))
    for x, row in zip(xs, feat.itertuples()):
        ax.annotate(TIER_LABEL.get(row.tier, row.tier), xy=(x, 0), xytext=(0, -34),
                    textcoords='offset points', ha='center', va='top',
                    fontsize=FS_SMALL - 6, color=TIER_COLOR.get(row.tier, '#333'),
                    fontweight='bold', annotation_clip=False)
    ax.set_ylabel('% of tier cells', fontsize=FS_TICK - 2)
    ax.tick_params(labelsize=FS_TICK - 4)
    ax.legend(handles=[Patch(facecolor=WT_GRAY, edgecolor='#333', label='expressing, no mutation'),
                       Patch(facecolor=NONTCW, edgecolor='#333',
                             label='carries mutation = prevalence (color = tier)')],
              fontsize=FS_SMALL - 6, frameon=False, loc='upper right')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ttl = "expression and neoantigen carriage"
    ax.set_title((f"{panel}  " if panel else "") + ttl, fontsize=FS_LABEL, loc='left', pad=10)


def draw_gene_track(ax, gene, track, lollipops, panel=None, legend='disorder'):
    """Native-frame backbone with functional domain boxes, muted disordered underlay,
    and the featured TCW-colored lollipop. Robust to a domain-empty track (short
    proteins render as backbone + one lollipop)."""
    L = max(track['native_len'], 1)
    ax.set_xlim(-L * 0.02, L * 1.02)
    ax.set_ylim(-2.0, 2.9)
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
                ha='center', va='center', fontsize=FS_SMALL - 6, zorder=4, clip_on=True)
    for lp in lollipops:
        x = lp['pos']
        ax.plot([x, x], [bar_y + bar_h, bar_y + bar_h + 1.3], '-', color='#555555', lw=2, zorder=5)
        ax.scatter([x], [bar_y + bar_h + 1.4], s=560, facecolors=lp['color'],
                   edgecolors='#222222', linewidths=2.4, zorder=6)
        ax.text(x, bar_y + bar_h + 1.75, lp['hgvs'], ha='center', va='bottom',
                fontsize=FS_SMALL - 6, zorder=7)
    label = f"{gene}  ({track['acc']}, {L} aa)"
    ax.text(0, bar_y + bar_h + 2.2, (f"{panel}  " if panel else "") + label,
            fontsize=FS_LABEL - 2, va='bottom')
    if legend:
        handles = []
        if track['disordered']:
            handles.append(Patch(facecolor=DISORDER, edgecolor='none', alpha=0.55,
                                 label='disordered region'))
        if legend == 'full':
            handles += [Line2D([0], [0], marker='o', linestyle='none', markersize=13,
                               markerfacecolor=c, markeredgecolor='#222', label=l)
                        for c, l in [(TCW_CT, 'SBS2 mutation (C>T)'),
                                     (TCW_CG, 'SBS13 mutation (C>G)'),
                                     (NONTCW, 'non-APOBEC mutation')]]
        if handles:
            ax.legend(handles=handles, loc='lower right', fontsize=FS_SMALL - 4, frameon=False)
    ax.set_yticks([])
    ax.set_xlabel('residue (native numbering)', fontsize=FS_TICK - 2)
    ax.tick_params(labelsize=FS_TICK - 4)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)


def draw_panelA(ax, stats, panel='A'):
    blocks = [('neoantigen', 'neoantigens'), ('fusion', 'RNA fusions')]
    xcenters, w = [0, 1.5], 0.5
    for (key, _), xc in zip(blocks, xcenters):
        d = stats.get(key)
        if not d:
            continue
        vals = [d['sbs2'], d['cnv']]
        sems = [d.get('sbs2_sem', np.nan), d.get('cnv_sem', np.nan)]
        xs = [xc - w / 2, xc + w / 2]
        ax.bar(xs, vals, width=w, color=[COLOR_SBS2, COLOR_CNV], edgecolor='#333333', linewidth=1.6,
               yerr=[[0, 0], [s if np.isfinite(s) else 0 for s in sems]],
               error_kw={'elinewidth': 2, 'capsize': 6, 'capthick': 2})
        top = max(v + (s if np.isfinite(s) else 0) for v, s in zip(vals, sems))
        br = top * 1.12
        ax.plot([xs[0], xs[0], xs[1], xs[1]], [top * 1.04, br, br, top * 1.04], color='#333333', lw=1.8)
        p = d['adj_p']
        ax.text(xc, br * 1.02, f"adj p = {p:.2g}" if p >= 1e-4 else f"adj p = {p:.1e}",
                ha='center', va='bottom', fontsize=FS_SMALL)
    ax.set_xticks(xcenters)
    ax.set_xticklabels(['neoantigens', 'RNA fusions'], fontsize=FS_TICK)
    ax.set_ylabel('events / 1000 UMI (per cell)', fontsize=FS_TICK - 2)
    ax.tick_params(labelsize=FS_TICK - 2)
    ax.legend(handles=[Patch(facecolor=COLOR_SBS2, edgecolor='#333', label='SBS2-HIGH'),
                       Patch(facecolor=COLOR_CNV, edgecolor='#333', label='CNV-HIGH')],
              fontsize=FS_SMALL, frameon=False, loc='upper left')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title((f"{panel}  " if panel else "") + 'mutational burden per cell',
                 fontsize=FS_LABEL, loc='left', pad=10)


def draw_panelB_venn(ax, counts, panel='B'):
    sbs2_only, shared, cnv_only = counts
    ax.set_xlim(0, 10); ax.set_ylim(0, 7); ax.axis('off')
    ax.add_patch(Circle((4.0, 3.4), 2.6, facecolor=COLOR_SBS2, edgecolor='#333', linewidth=1.8, alpha=0.55))
    ax.add_patch(Circle((6.0, 3.4), 2.6, facecolor=COLOR_CNV, edgecolor='#333', linewidth=1.8, alpha=0.55))
    ax.text(2.8, 3.4, f"{sbs2_only}", ha='center', va='center', fontsize=FS_TITLE)
    ax.text(7.2, 3.4, f"{cnv_only}", ha='center', va='center', fontsize=FS_TITLE)
    ax.text(5.0, 3.4, f"{shared}", ha='center', va='center', fontsize=FS_LABEL)
    ax.text(2.8, 0.5, 'SBS2-specific', ha='center', fontsize=FS_SMALL, color='#C0392B')
    ax.text(7.2, 0.5, 'CNV-specific', ha='center', fontsize=FS_SMALL, color='#B7950B')
    ax.text(5.0, 6.4, 'neoantigen genes', ha='center', fontsize=FS_SMALL)
    ax.set_title((f"{panel}  " if panel else "") + 'neoantigen gene overlap', fontsize=FS_LABEL, loc='left')


def save_figure(fig, name):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{ext}"), dpi=300,
                    bbox_inches='tight', facecolor='white')
    return os.path.join(OUTPUT_DIR, f"{name}.pdf")
