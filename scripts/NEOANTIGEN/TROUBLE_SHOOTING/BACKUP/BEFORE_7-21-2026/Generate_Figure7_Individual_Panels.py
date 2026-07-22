#!/usr/bin/env python3
"""
Generate_Figure7_Individual_Panels.py
=====================================
Renders every Figure 7 panel as its own standalone PDF + PNG, in addition to the
combined main and supplemental figures. This makes it easy to drop individual
panels into Illustrator and rearrange them.

Outputs (data/FIG_7/figures/):
  Figure7_panel_burden.*                 mutational burden per cell
  Figure7_panel_venn.*                   neoantigen gene overlap
  Figure7_panel_expression_carriage.*    expression + carriage, all five genes
  Figure7_{GENE}_binding.*               binding gain, per gene (5)
  Figure7_{GENE}_track.*                 protein track, per gene (5)

Standalone binding panels carry the TCW color legend; standalone tracks carry the
full disordered + mutation legend, so each file is self-contained.

Run in NETWORK env. Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import neoantigen_figure_utils as U

ALL_GENES = ['KLF3', 'CAST', 'SERPINB2', 'KRT6B', 'ANXA1']


def _save_one(name, figsize, draw_fn, binding_legend=False):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    draw_fn(ax)
    if binding_legend:
        fig.legend(handles=U.binding_legend_handles(), fontsize=U.FS_SMALL - 2,
                   frameon=False, loc='upper right', bbox_to_anchor=(0.995, 0.995))
    U.save_figure(fig, name)
    plt.close(fig)
    print(f"  wrote {name}")


def main():
    fasta = U.load_fasta_tcw()
    panelA = U.load_panelA()
    venn = U.load_venn()
    cards = {g: U.load_gene_card(g, fasta=fasta) for g in ALL_GENES}

    # summary panels
    _save_one("Figure7_panel_burden", (11, 8),
              lambda ax: U.draw_panelA(ax, panelA, panel=None))
    _save_one("Figure7_panel_venn", (9, 7),
              lambda ax: U.draw_panelB_venn(ax, venn, panel=None))
    _save_one("Figure7_panel_expression_carriage", (16, 7),
              lambda ax: U.draw_expr_carrier_stacked(ax, ALL_GENES, panel=None))

    # per-gene binding + track
    for g in ALL_GENES:
        card = cards[g]
        n = max(len(card['binding']), 1)
        _save_one(f"Figure7_{g}_binding", (13, max(4.0, 1.15 * n + 2.2)),
                  lambda ax, c=card, gg=g: U.draw_binding_bars(ax, gg, c['binding'], panel=None),
                  binding_legend=True)
        _save_one(f"Figure7_{g}_track", (18, 5),
                  lambda ax, c=card, gg=g: U.draw_gene_track(
                      ax, gg, c['track'], c['lollipops'], panel=None, legend='full'))


if __name__ == "__main__":
    main()
