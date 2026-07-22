#!/usr/bin/env python3
"""
Generate_Figure7_Individual_Panels.py  (single-source, three featured candidates)
=================================================================================
Renders every Figure 7 panel as its own standalone PDF + PNG for Illustrator.
All numbers read from neoantigen_prevalence_ranking_full.tsv via the utils module.

Outputs (data/FIG_7/figures/):
  Figure7_panel_burden.*                mutational burden per cell
  Figure7_panel_venn.*                  neoantigen gene overlap
  Figure7_panel_expression_carriage.*   expression + carriage, three genes
  Figure7_panel_binding.*               MHC-I binding gain, three featured mutations
  Figure7_{GENE}_track.*                protein track, per featured gene (3)

Standalone binding panel carries the TCW legend; standalone tracks carry the
full disordered + mutation legend, so each file is self-contained.
Run in NETWORK env. Author: Jake Lehle / Claude (2026 NMF Paper)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import neoantigen_figure_utils as U


def _save_one(name, figsize, draw_fn):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    draw_fn(ax)                       # panels draw their own legends
    U.save_figure(fig, name)
    plt.close(fig)
    print(f"  wrote {name}")


def main():
    feat = U.load_featured()
    panelA = U.load_panelA()
    venn = U.load_venn()

    _save_one("Figure7_panel_burden", (11, 8),
              lambda ax: U.draw_panelA(ax, panelA, panel=None))
    _save_one("Figure7_panel_venn", (9, 7),
              lambda ax: U.draw_panelB_venn(ax, venn, panel=None))
    _save_one("Figure7_panel_expression_carriage", (12, 8),
              lambda ax: U.draw_expr_featured(ax, feat, panel=None))
    _save_one("Figure7_panel_binding", (14, 6),
              lambda ax: U.draw_binding_featured(ax, feat, panel=None))

    for _, row in feat.iterrows():
        gene = row['gene']
        track = U.load_gene_track(gene)
        lolli = U.featured_lollipop(row)
        width = 18 if gene == U.TRACK_LONG else 11
        _save_one(f"Figure7_{gene}_track", (width, 5),
                  lambda ax, g=gene, t=track, l=lolli: U.draw_gene_track(
                      ax, g, t, l, panel=None, legend='full'))


if __name__ == "__main__":
    main()
