#!/usr/bin/env python3
"""
Generate_Figure7_Supplemental.py
================================
Supplemental Figure 7. Every element on its own full-width row so nothing is
compressed:
  Row A   expression + neoantigen carriage (vertical stacked bars, all five genes)
  then, for SERPINB2, KRT6B, ANXA1 (in order):
      binding gain (horizontal stacked bars)   -- its own row
      protein track (native frame)             -- its own full-width row

All numbers read from locked outputs via the utils module. Gene names plain.

Run in NETWORK env. Writes PDF + PNG (300 DPI) to data/FIG_7/figures/.
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import neoantigen_figure_utils as U

ALL_GENES = ['KLF3', 'CAST', 'SERPINB2', 'KRT6B', 'ANXA1']
SUPP_GENES = ['SERPINB2', 'KRT6B', 'ANXA1']


def main():
    fasta = U.load_fasta_tcw()
    cards = {g: U.load_gene_card(g, fasta=fasta) for g in SUPP_GENES}

    # row plan: Panel A, then (binding, track) per gene. Binding rows scale with
    # the number of mutations so ANXA1's nine rows get room and KRT6B's two don't
    # waste it; tracks get a fixed generous height.
    rows = [('exprA', None, 0.9)]
    for g in SUPP_GENES:
        n = max(len(cards[g]['binding']), 1)
        rows.append(('binding', g, 0.30 + 0.135 * n))
        rows.append(('track', g, 0.60))
    ratios = [r[2] for r in rows]

    fig = plt.figure(figsize=(22, sum(ratios) * 5.4))
    gs = GridSpec(len(rows), 1, height_ratios=ratios, hspace=0.75, figure=fig)

    letters = iter(['B', 'C', 'D'])
    for i, (kind, gene, _) in enumerate(rows):
        ax = fig.add_subplot(gs[i])
        if kind == 'exprA':
            U.draw_expr_carrier_stacked(ax, ALL_GENES, panel='A')
        elif kind == 'binding':
            U.draw_binding_bars(ax, gene, cards[gene]['binding'], panel=next(letters))
        else:
            U.draw_gene_track(ax, gene, cards[gene]['track'], cards[gene]['lollipops'],
                              legend='disorder')

    fig.legend(handles=U.binding_legend_handles(), fontsize=U.FS_SMALL - 2,
               frameon=False, loc='upper right', bbox_to_anchor=(0.995, 0.997))
    fig.suptitle("Figure 7 - supplement  Neoantigen carriage and the ANXA1 "
                 "off-signature counterpoint", fontsize=U.FS_TITLE, y=1.002)
    path = U.save_figure(fig, "Figure7_supplemental")
    plt.close(fig)
    print(f"Wrote {path} (+ .png)")


if __name__ == "__main__":
    main()
