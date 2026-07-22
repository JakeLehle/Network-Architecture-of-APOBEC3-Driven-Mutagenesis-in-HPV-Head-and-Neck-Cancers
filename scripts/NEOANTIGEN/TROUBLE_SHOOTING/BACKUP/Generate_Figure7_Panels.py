#!/usr/bin/env python3
"""
Generate_Figure7_Panels.py  (single-source, three featured candidates)
======================================================================
Main Figure 7 (neoantigen landscape). Merged four-row layout:

  Row 1  A burden per cell | B neoantigen gene overlap (Venn) | C expression + carriage
  Row 2  D MHC-I binding gain (three featured mutations, horizontal stacked bars)
  Row 3  E KRT6B protein track (full width; the long protein with real domains)
  Row 4  F COX4I1 track | G SPRR1A track (paired short proteins)

All binding / prevalence / expression numbers are read from
neoantigen_prevalence_ranking_full.tsv via the utils module, so the panels and
the manuscript cannot disagree. Panel B (Venn) and Panel A (burden) read their
own locked sources. Two color languages: TIER palette for prevalence/expression,
TCW palette for binding/track provenance (see utils).

Run in NETWORK env. Writes PDF + PNG (300 DPI) to data/FIG_7/figures/.
Prerequisite: protein_domains.tsv must contain COX4I1, SPRR1A, KRT6B
(re-run Diagnostic_Fetch_Protein_Domains.py with the three-gene FIGURE_GENES).
Author: Jake Lehle / Claude (2026 NMF Paper)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import neoantigen_figure_utils as U


def main():
    feat = U.load_featured()                       # 3 rows, prevalence-ordered
    panelA = U.load_panelA()
    venn = U.load_venn()
    tracks = {g: U.load_gene_track(g) for g in [U.TRACK_LONG] + U.TRACK_SHORT}
    lolli = {row['gene']: U.featured_lollipop(row) for _, row in feat.iterrows()}

    fig = plt.figure(figsize=(24, 23))
    gs = GridSpec(4, 1, height_ratios=[1.02, 0.62, 0.52, 0.52], hspace=0.70, figure=fig)

    # ---- Row 1: burden | Venn | expression+carriage
    r1 = gs[0].subgridspec(1, 3, width_ratios=[0.36, 0.28, 0.36], wspace=0.34)
    U.draw_panelA(fig.add_subplot(r1[0, 0]), panelA, panel='A')
    U.draw_panelB_venn(fig.add_subplot(r1[0, 1]), venn, panel='B')
    U.draw_expr_featured(fig.add_subplot(r1[0, 2]), feat, panel='C')

    # ---- Row 2: binding (three featured mutations, full width)
    U.draw_binding_featured(fig.add_subplot(gs[1]), feat, panel='D')

    # ---- Row 3: KRT6B track, full width
    U.draw_gene_track(fig.add_subplot(gs[2]), U.TRACK_LONG, tracks[U.TRACK_LONG],
                      lolli[U.TRACK_LONG], panel='E', legend='full')

    # ---- Row 4: COX4I1 | SPRR1A tracks (paired short)
    r4 = gs[3].subgridspec(1, 2, wspace=0.16)
    for i, (gene, letter) in enumerate(zip(U.TRACK_SHORT, ['F', 'G'])):
        U.draw_gene_track(fig.add_subplot(r4[0, i]), gene, tracks[gene],
                          lolli[gene], panel=letter, legend='disorder')

    # TCW legend now lives inside the binding panel (D); tier palette is keyed by
    # the colored sublabels in the expression panel (C). No figure-level legend.
    fig.suptitle("Figure 7  Neoantigen landscape of APOBEC3-driven mutagenesis",
                 fontsize=U.FS_TITLE + 2, y=1.005)
    path = U.save_figure(fig, "Figure7_main")
    plt.close(fig)
    print(f"Wrote {path} (+ .png)")


if __name__ == "__main__":
    main()
