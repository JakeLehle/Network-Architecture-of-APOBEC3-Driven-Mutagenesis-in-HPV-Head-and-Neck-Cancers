#!/usr/bin/env python3
"""
Tier2A_Expression_Haplotype_Proxies.py
=======================================

Expression-level haplotype proxy measurements per patient:
  1. A3A/A3B ratio (deletion polymorphism proxy)
  2. A3H expression (haplotype stability proxy)
  3. A3 interactor anchor expression (RALY, HNRNPA2B1)
  4. Activating chain gene expression (9-gene program from Figure 4)

Updated May 2026: gene sets reflect Figure 4 concordance analysis.
The activating chain (RALY/HNRNPA2B1 anchored) is the coordinately
upregulated program in SBS2-HIGH cells. Per-patient expression of
these genes tests whether high-contributor patients show elevated
baseline expression of the mutagenic program.

Usage: conda run -n NETWORK python Tier2A_Expression_Haplotype_Proxies.py
"""

import os, numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

from patient_config import (
    DIR_00_DIAG, DIR_02_SNP, PATIENT_COL, CELLTYPE_COL,
    A3_GENES_SYMBOLS, A3_INTERACTOR_ANCHORS, ACTIVATING_CHAIN_GENES,
    INHIBITING_CHAIN_ANCHORS, HIGH_CONTRIBUTORS,
    FONT_TITLE, FONT_LABEL, FONT_TICK,
    COLOR_SBS2_HIGH, COLOR_NORMAL,
    banner, log, ensure_dir, load_adata, get_gene_expression,
)


def main():
    banner("TIER 2A: EXPRESSION HAPLOTYPE PROXIES")
    out_dir = ensure_dir(DIR_02_SNP)
    adata = load_adata()

    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
    log(f"  Basal cells: {basal.n_obs:,}")

    # Load enrichment data
    enrich_path = os.path.join(DIR_00_DIAG, "patient_enrichment_SBS2_HIGH_v2.tsv")
    enrichment_df = pd.read_csv(enrich_path, sep='\t') if os.path.exists(enrich_path) else None

    patients = sorted(basal.obs[PATIENT_COL].unique())
    rows = []

    for patient in patients:
        mask = basal.obs[PATIENT_COL] == patient
        subset = basal[mask]
        row = {'patient': patient, 'n_basal': mask.sum()}

        # A3 family expression
        for gene in A3_GENES_SYMBOLS:
            expr = get_gene_expression(subset, gene)
            row[f'{gene}_mean'] = expr.mean() if expr is not None else 0.0

        # A3A/A3B ratio (deletion polymorphism proxy)
        a3a = row.get('APOBEC3A_mean', 0)
        a3b = row.get('APOBEC3B_mean', 0)
        row['A3A_A3B_ratio'] = a3a / a3b if a3b > 0 else np.nan

        # A3 interactor anchors (RALY, HNRNPA2B1)
        anchor_means = []
        for gene in A3_INTERACTOR_ANCHORS:
            expr = get_gene_expression(subset, gene)
            val = expr.mean() if expr is not None else 0.0
            row[f'{gene}_mean'] = val
            anchor_means.append(val)
        row['anchor_mean'] = np.mean(anchor_means)

        # Full activating chain (9 genes)
        chain_means = []
        chain_pct_nz = []
        for gene in ACTIVATING_CHAIN_GENES:
            expr = get_gene_expression(subset, gene)
            val = expr.mean() if expr is not None else 0.0
            pct = ((expr > 0).mean() * 100) if expr is not None else 0.0
            row[f'{gene}_mean'] = val
            row[f'{gene}_pct_nz'] = pct
            chain_means.append(val)
            chain_pct_nz.append(pct)
        row['activating_chain_mean'] = np.mean(chain_means)
        row['activating_chain_pct_detected'] = np.mean(chain_pct_nz)

        # Inhibiting chain anchors (for completeness)
        inh_means = []
        for gene in INHIBITING_CHAIN_ANCHORS:
            expr = get_gene_expression(subset, gene)
            val = expr.mean() if expr is not None else 0.0
            row[f'{gene}_mean'] = val
            inh_means.append(val)
        row['inhibiting_anchor_mean'] = np.mean(inh_means)

        rows.append(row)

    proxy_df = pd.DataFrame(rows)

    # Merge enrichment
    if enrichment_df is not None:
        proxy_df = proxy_df.merge(
            enrichment_df[['patient', 'fold_enrichment', 'n_high', 'high_fraction']],
            on='patient', how='left')
        proxy_df['fold_enrichment'] = proxy_df['fold_enrichment'].fillna(0)
        proxy_df['n_high'] = proxy_df['n_high'].fillna(0).astype(int)

    proxy_path = os.path.join(out_dir, "Tier2A_haplotype_proxies.tsv")
    proxy_df.to_csv(proxy_path, sep='\t', index=False)
    log(f"  Saved: {proxy_path}")

    # =========================================================================
    # SUMMARY TABLE
    # =========================================================================
    banner("PROXY SUMMARY")
    log(f"  {'Patient':<20s}  {'A3A':>6s}  {'A3B':>6s}  {'Ratio':>6s}  "
        f"{'A3H':>6s}  {'Anchor':>6s}  {'Chain':>6s}  {'Fold':>6s}")
    for _, r in proxy_df.sort_values('fold_enrichment', ascending=False).iterrows():
        flag = " <<<" if r['patient'] in HIGH_CONTRIBUTORS else ""
        log(f"  {r['patient']:<20s}  "
            f"{r['APOBEC3A_mean']:>6.3f}  {r['APOBEC3B_mean']:>6.3f}  "
            f"{r['A3A_A3B_ratio']:>6.2f}  {r['APOBEC3H_mean']:>6.3f}  "
            f"{r['anchor_mean']:>6.3f}  {r['activating_chain_mean']:>6.3f}  "
            f"{r.get('fold_enrichment', 0):>5.1f}x{flag}")

    # =========================================================================
    # CORRELATIONS WITH FOLD ENRICHMENT
    # =========================================================================
    if 'fold_enrichment' in proxy_df.columns:
        banner("CORRELATIONS WITH FOLD ENRICHMENT")
        test_cols = [
            ('A3A_A3B_ratio', 'A3A/A3B ratio'),
            ('APOBEC3A_mean', 'A3A expression'),
            ('APOBEC3B_mean', 'A3B expression'),
            ('APOBEC3H_mean', 'A3H expression'),
            ('anchor_mean', 'A3 interactor anchors (RALY+HNRNPA2B1)'),
            ('activating_chain_mean', 'Activating chain (9 genes)'),
            ('inhibiting_anchor_mean', 'Inhibiting chain anchors'),
        ]
        # Also test individual chain genes
        for gene in ACTIVATING_CHAIN_GENES:
            test_cols.append((f'{gene}_mean', f'{gene}'))

        for col, label in test_cols:
            if col not in proxy_df.columns:
                continue
            valid = proxy_df[[col, 'fold_enrichment']].dropna()
            if len(valid) >= 5:
                rho, pval = spearmanr(valid[col], valid['fold_enrichment'])
                sig = "*" if pval < 0.05 else ""
                log(f"  {label:<45s}  rho={rho:+.3f}, p={pval:.3f} {sig}")

    # =========================================================================
    # FIGURE: HAPLOTYPE PROXIES + CHAIN EXPRESSION
    # =========================================================================
    banner("PLOTTING")
    metrics = [
        ('A3A_A3B_ratio',
         'A3A/A3B Expression Ratio\n(deletion polymorphism proxy)'),
        ('APOBEC3H_mean',
         'Mean A3H Expression\n(haplotype stability proxy)'),
        ('anchor_mean',
         'A3 Interactor Anchor Expression\n(RALY + HNRNPA2B1)'),
        ('activating_chain_mean',
         'Activating Chain Expression\n(9-gene coordinated program)'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    for ax, (col, title) in zip(axes.flat, metrics):
        if col not in proxy_df.columns or 'fold_enrichment' not in proxy_df.columns:
            ax.set_visible(False)
            continue

        valid = proxy_df[[col, 'fold_enrichment', 'patient']].dropna()
        colors = [COLOR_SBS2_HIGH if p in HIGH_CONTRIBUTORS
                  else COLOR_NORMAL for p in valid['patient']]

        ax.scatter(valid[col], valid['fold_enrichment'], c=colors, s=120,
                   edgecolors='black', linewidth=0.8, zorder=3)

        for _, r in valid.iterrows():
            ax.annotate(r['patient'].replace('Patient ', ''),
                        (r[col], r['fold_enrichment']),
                        fontsize=8, ha='center', va='bottom',
                        xytext=(0, 6), textcoords='offset points')

        ax.axhline(1.0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)

        rho, pval = spearmanr(valid[col], valid['fold_enrichment'])
        short_title = title.split('\n')[0]
        ax.set_title(f'{short_title}\nrho={rho:+.3f}, p={pval:.3f}',
                     fontsize=FONT_TICK, fontweight='bold')
        ax.set_xlabel(title, fontsize=FONT_TICK)
        ax.set_ylabel('SBS2-HIGH Fold Enrichment', fontsize=FONT_TICK)
        ax.tick_params(labelsize=FONT_TICK - 4)

    plt.suptitle('Expression Haplotype Proxies vs SBS2-HIGH Enrichment',
                 fontsize=FONT_TITLE, fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(out_dir, f"Tier2A_haplotype_proxies.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved proxy scatter plots")

    # =========================================================================
    # FIGURE: PER-GENE CHAIN HEATMAP
    # =========================================================================
    banner("CHAIN GENE HEATMAP")
    chain_cols = [f'{g}_mean' for g in ACTIVATING_CHAIN_GENES
                  if f'{g}_mean' in proxy_df.columns]

    if chain_cols and 'fold_enrichment' in proxy_df.columns:
        plot_df = proxy_df.sort_values('fold_enrichment', ascending=False)
        mat = plot_df[chain_cols].values
        gene_labels = [c.replace('_mean', '') for c in chain_cols]
        patient_labels = [p.replace('Patient ', '') for p in plot_df['patient']]

        fig, ax = plt.subplots(figsize=(12, 8))
        im = ax.imshow(mat, cmap='YlOrRd', aspect='auto')

        ax.set_xticks(range(len(gene_labels)))
        ax.set_xticklabels(gene_labels, rotation=45, ha='right',
                           fontsize=FONT_TICK)
        ax.set_yticks(range(len(patient_labels)))
        ax.set_yticklabels(patient_labels, fontsize=FONT_TICK)

        # Highlight high contributors
        for i, p in enumerate(plot_df['patient']):
            if p in HIGH_CONTRIBUTORS:
                ax.get_yticklabels()[i].set_color(COLOR_SBS2_HIGH)
                ax.get_yticklabels()[i].set_fontweight('bold')

        # Highlight A3 interactor anchors
        for i, g in enumerate(gene_labels):
            if g in A3_INTERACTOR_ANCHORS:
                ax.get_xticklabels()[i].set_fontweight('bold')

        plt.colorbar(im, ax=ax, label='Mean Expression', shrink=0.7)
        ax.set_title('Activating Chain Gene Expression by Patient\n'
                     '(sorted by SBS2-HIGH fold enrichment)',
                     fontsize=FONT_LABEL, fontweight='bold')
        plt.tight_layout()
        for ext in ['pdf', 'png']:
            plt.savefig(os.path.join(out_dir,
                        f"Tier2A_activating_chain_heatmap.{ext}"),
                        dpi=300, bbox_inches='tight')
        plt.close()
        log("  Saved chain gene heatmap")

    log("\nTier 2A complete.")


if __name__ == "__main__":
    main()
