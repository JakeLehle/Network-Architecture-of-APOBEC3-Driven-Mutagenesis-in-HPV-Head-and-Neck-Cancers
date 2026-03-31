#!/usr/bin/env python3
"""
Tier2A_Expression_Haplotype_Proxies.py
=======================================

Expression-level haplotype proxy measurements per patient:
  1. A3A/A3B ratio (deletion polymorphism proxy)
  2. A3H expression (haplotype stability proxy)
  3. Mean C0 interactor expression

Usage: conda run -n NETWORK python Tier2A_Expression_Haplotype_Proxies.py
"""

import os, numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

from patient_config import *

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

        # A3 expression
        for gene in A3_GENES_SYMBOLS:
            expr = get_gene_expression(subset, gene)
            row[f'{gene}_mean'] = expr.mean() if expr is not None else 0.0

        # A3A/A3B ratio
        a3a = row.get('APOBEC3A_mean', 0)
        a3b = row.get('APOBEC3B_mean', 0)
        row['A3A_A3B_ratio'] = a3a / a3b if a3b > 0 else np.nan

        # C0 interactors (mean across the set)
        interactor_means = []
        for gene in C0_INTERACTORS:
            expr = get_gene_expression(subset, gene)
            val = expr.mean() if expr is not None else 0.0
            row[f'{gene}_mean'] = val
            interactor_means.append(val)
        row['C0_interactor_mean'] = np.mean(interactor_means)

        rows.append(row)

    proxy_df = pd.DataFrame(rows)

    # Merge enrichment
    if enrichment_df is not None:
        proxy_df = proxy_df.merge(enrichment_df[['patient','fold_enrichment','n_high','high_fraction']],
                                   on='patient', how='left')
        proxy_df['fold_enrichment'] = proxy_df['fold_enrichment'].fillna(0)
        proxy_df['n_high'] = proxy_df['n_high'].fillna(0).astype(int)

    proxy_path = os.path.join(out_dir, "Tier2A_haplotype_proxies.tsv")
    proxy_df.to_csv(proxy_path, sep='\t', index=False)
    log(f"  Saved: {proxy_path}")

    # Print summary
    banner("PROXY SUMMARY")
    log(f"  {'Patient':<20s}  {'A3A':>6s}  {'A3B':>6s}  {'Ratio':>6s}  {'A3H':>6s}  {'Intrc':>6s}  {'Fold':>6s}")
    for _, r in proxy_df.sort_values('fold_enrichment', ascending=False).iterrows():
        flag = " <<<" if r['patient'] in HIGH_CONTRIBUTORS else ""
        log(f"  {r['patient']:<20s}  {r['APOBEC3A_mean']:>6.3f}  {r['APOBEC3B_mean']:>6.3f}  "
            f"{r['A3A_A3B_ratio']:>6.2f}  {r['APOBEC3H_mean']:>6.3f}  "
            f"{r['C0_interactor_mean']:>6.3f}  {r.get('fold_enrichment',0):>5.1f}x{flag}")

    # Correlation tests
    if 'fold_enrichment' in proxy_df.columns:
        banner("CORRELATIONS WITH FOLD ENRICHMENT")
        for col in ['A3A_A3B_ratio', 'APOBEC3A_mean', 'APOBEC3B_mean',
                     'APOBEC3H_mean', 'C0_interactor_mean']:
            valid = proxy_df[[col, 'fold_enrichment']].dropna()
            if len(valid) >= 5:
                rho, pval = spearmanr(valid[col], valid['fold_enrichment'])
                sig = "*" if pval < 0.05 else ""
                log(f"  {col} vs fold: rho={rho:.3f}, p={pval:.3f} {sig}")

    # Visualization
    banner("PLOTTING")
    metrics = [
        ('A3A_A3B_ratio', 'A3A/A3B Expression Ratio\n(deletion polymorphism proxy)'),
        ('APOBEC3H_mean', 'Mean A3H Expression\n(haplotype stability proxy)'),
        ('C0_interactor_mean', 'Mean C0 Interactor Expression\n(HNRNPA2B1+HSPD1+RPL5+TIMM8B)'),
        ('APOBEC3A_mean', 'Mean A3A Expression'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    for ax, (col, title) in zip(axes.flat, metrics):
        if col not in proxy_df.columns or 'fold_enrichment' not in proxy_df.columns:
            ax.set_visible(False); continue
        valid = proxy_df[[col, 'fold_enrichment', 'patient']].dropna()
        colors = ['#ed6a5a' if p in HIGH_CONTRIBUTORS else '#5e81ac' for p in valid['patient']]
        ax.scatter(valid[col], valid['fold_enrichment'], c=colors, s=120,
                  edgecolors='black', linewidth=0.8, zorder=3)
        for _, r in valid.iterrows():
            ax.annotate(r['patient'].replace('Patient ',''), (r[col], r['fold_enrichment']),
                       fontsize=8, ha='center', va='bottom', xytext=(0,6), textcoords='offset points')
        ax.axhline(1.0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.set_xlabel(title, fontsize=11)
        ax.set_ylabel('SBS2-HIGH Fold Enrichment', fontsize=11)

        rho, pval = spearmanr(valid[col], valid['fold_enrichment'])
        ax.set_title(f'{title.split(chr(10))[0]}\nrho={rho:.3f}, p={pval:.3f}',
                    fontsize=12, fontweight='bold')

    plt.suptitle('Expression Haplotype Proxies vs SBS2-HIGH Enrichment',
                fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(os.path.join(out_dir, f"Tier2A_haplotype_proxies.{ext}"), dpi=300, bbox_inches='tight')
    plt.close()
    log("\nTier 2A complete.")

if __name__ == "__main__":
    main()
