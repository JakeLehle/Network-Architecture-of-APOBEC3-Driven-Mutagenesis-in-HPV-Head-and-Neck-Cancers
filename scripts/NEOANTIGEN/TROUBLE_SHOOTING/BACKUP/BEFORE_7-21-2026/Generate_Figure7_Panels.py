#!/usr/bin/env python3
"""
Generate_Figure7_Panels.py
==========================
Main Figure 7 (neoantigen landscape). Layout, top to bottom:
  Row 1  A KLF3 binding gain   |   B CAST binding gain   (horizontal stacked bars)
  Row 2  C KLF3 protein track  (full width, native frame)
  Row 3  D CAST protein track  (full width, native frame)
  Row 4  E burden per cell     |   F neoantigen gene overlap (Venn)

Binding is a stacked bar (wild-type baseline + mutation binding gain); the gene
tracks get their own full-width rows to spread out. Expression/carrier has moved
to the supplemental. All numbers read from locked outputs via the utils module.

Run in NETWORK env. Writes PDF + PNG (300 DPI) to data/FIG_7/figures/.
Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import neoantigen_figure_utils as U

MAIN_GENES = ['KLF3', 'CAST']


def main():
    fasta = U.load_fasta_tcw()
    panelA = U.load_panelA()
    venn = U.load_venn()
    cards = {g: U.load_gene_card(g, fasta=fasta) for g in MAIN_GENES}

    fig = plt.figure(figsize=(24, 26))
    gs = GridSpec(4, 1, height_ratios=[0.92, 0.58, 0.58, 0.78], hspace=0.62, figure=fig)

    # ---- Row 1: binding (KLF3 | CAST)
    r1 = gs[0].subgridspec(1, 2, wspace=0.34)
    ax_kb = fig.add_subplot(r1[0, 0])
    ax_cb = fig.add_subplot(r1[0, 1])
    U.draw_binding_bars(ax_kb, 'KLF3', cards['KLF3']['binding'], panel='A')
    U.draw_binding_bars(ax_cb, 'CAST', cards['CAST']['binding'], panel='B')

    # ---- Row 2 & 3: full-width tracks
    ax_kt = fig.add_subplot(gs[1])
    ax_ct = fig.add_subplot(gs[2])
    U.draw_gene_track(ax_kt, 'KLF3', cards['KLF3']['track'], cards['KLF3']['lollipops'], panel='C')
    U.draw_gene_track(ax_ct, 'CAST', cards['CAST']['track'], cards['CAST']['lollipops'], panel='D')

    # ---- Row 4: burden | Venn
    r4 = gs[3].subgridspec(1, 2, width_ratios=[0.55, 0.45], wspace=0.28)
    ax_e = fig.add_subplot(r4[0, 0])
    ax_f = fig.add_subplot(r4[0, 1])
    U.draw_panelA(ax_e, panelA, panel='E')
    U.draw_panelB_venn(ax_f, venn, panel='F')

    # ---- shared binding legend (one place, top-right of figure)
    fig.legend(handles=U.binding_legend_handles(), fontsize=U.FS_SMALL - 2,
               frameon=False, loc='upper right', bbox_to_anchor=(0.995, 0.995))

    fig.suptitle("Figure 7  Neoantigen landscape of APOBEC3-driven mutagenesis",
                 fontsize=U.FS_TITLE + 2, y=1.005)
    path = U.save_figure(fig, "Figure7_main")
    plt.close(fig)
    print(f"Wrote {path} (+ .png)")


if __name__ == "__main__":
    main()
