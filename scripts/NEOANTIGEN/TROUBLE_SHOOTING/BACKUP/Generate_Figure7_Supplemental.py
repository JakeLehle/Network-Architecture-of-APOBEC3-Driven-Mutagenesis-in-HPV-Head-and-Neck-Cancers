#!/usr/bin/env python3
"""
Generate_Figure7_Supplemental.py
================================
Figure 7 SUPPLEMENTAL: extended A3-lead cards.

Three stacked cards, same design as the main figure's KLF3 / CAST cards:
  a) SERPINB2 (Tier 1, shared)        clean-A3 shared mirror to KLF3
  b) KRT6B    (Tier 2, SBS2-specific) clean-A3 squamous-differentiation gene
  c) ANXA1    (Tier 2, SBS2-specific) OFF-SIGNATURE COUNTERPOINT

The ANXA1 card is deliberately not an A3 exemplar. Its single TpCpW event is a
C>G (SBS13), not a clean C>T (SBS2), so its context glyph is grayed and tagged
"SBS13, not SBS2". Its excellent immunogenicity (best binder ~15 nM) comes from
NON-APOBEC mutations, which the binding sub-panel marks separately so the card
reads as "highly expressed and immunogenic, but not A3-driven".

Every number is read from disk via neoantigen_figure_utils. PDF + PNG at 300 DPI.

Run in the NETWORK conda env.

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import warnings
warnings.filterwarnings('ignore')

import neoantigen_figure_utils as U

plt.rcParams.update(U.BASE_RC)

OUTPUT_DIR = os.path.join(U.FIG7_ROOT, "FIGURE_7_PANELS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

SUPP_GENES = ['SERPINB2', 'KRT6B', 'ANXA1']
OFFSIG = ('ANXA1',)


def _save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{ext}"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {name}.pdf / .png")


def main():
    print("=" * 70)
    print("  FIGURE 7 SUPPLEMENTAL  (SERPINB2 + KRT6B + ANXA1)")
    print("=" * 70)
    records = U.build_lead_records(SUPP_GENES, offsignature_genes=OFFSIG)

    letters = ['a', 'b', 'c']

    # ---- individual cards -----------------------------------------------
    for rec, letter in zip(records, letters):
        fig = plt.figure(figsize=(26, 8))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.0, 1.7, 0.9], wspace=0.30)
        axc = fig.add_subplot(gs[0, 0])
        axb = fig.add_subplot(gs[0, 1])
        axe = fig.add_subplot(gs[0, 2])
        U.draw_lead_card(axc, axb, axe, rec, letter=letter)
        fig.subplots_adjust(top=0.78, bottom=0.22, left=0.04, right=0.98)
        _save(fig, f"Figure7_Supp_{letter}_{rec['gene']}")

    # ---- combined supplemental (three stacked cards) --------------------
    fig = plt.figure(figsize=(28, 30))
    outer = gridspec.GridSpec(3, 1, hspace=0.55)
    for idx, (rec, letter) in enumerate(zip(records, letters)):
        row = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[idx],
                                               width_ratios=[1.0, 1.7, 0.9], wspace=0.30)
        axc = fig.add_subplot(row[0, 0])
        axb = fig.add_subplot(row[0, 1])
        axe = fig.add_subplot(row[0, 2])
        U.draw_lead_card(axc, axb, axe, rec, letter=letter)

    _save(fig, "Figure7_Supplemental_Combined")

    print("\nAll supplemental panels written to:", OUTPUT_DIR)
    print("Done.")


if __name__ == "__main__":
    main()
