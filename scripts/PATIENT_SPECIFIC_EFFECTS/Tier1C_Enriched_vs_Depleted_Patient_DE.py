#!/usr/bin/env python3
"""
Tier1C_Enriched_vs_Depleted_Patient_DE.py
==========================================

Tier 1C: DE + unbiased GSEA comparing ALL basal cells from high-contributor
patients (SC029, SC013, SC001) vs remaining 11 patients.

Asks: what is constitutively different about these patients' basal cells?

Usage: conda run -n NETWORK python Tier1C_Enriched_vs_Depleted_Patient_DE.py
"""

import os, sys, numpy as np, pandas as pd, scanpy as sc, scipy.sparse
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp

from patient_config import *

def main():
    banner("TIER 1C: HIGH-CONTRIBUTOR VS OTHER PATIENTS (DE + GSEA)")
    out_dir = ensure_dir(DIR_01_EXPRESSION)
    adata = load_adata()

    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
    log(f"  Basal cells: {basal.n_obs:,}")

    # Split by patient group
    in_high = basal.obs[PATIENT_COL].isin(HIGH_CONTRIBUTORS)
    log(f"  High-contributor basal cells: {in_high.sum():,} ({HIGH_CONTRIBUTORS})")
    log(f"  Other patients' basal cells: {(~in_high).sum():,}")

    # DE: Wilcoxon rank-sum per gene
    banner("DIFFERENTIAL EXPRESSION")
    X = basal.X
    if scipy.sparse.issparse(X):
        log("  Converting sparse -> dense..."); X = X.toarray()
    genes = basal.var_names.values
    mask_hc = in_high.values

    stats, pvals, mean_hc, mean_other = [], [], [], []
    for g in range(X.shape[1]):
        x_hc = X[mask_hc, g]
        x_ot = X[~mask_hc, g]
        mean_hc.append(x_hc.mean())
        mean_other.append(x_ot.mean())
        if x_hc.std() == 0 and x_ot.std() == 0:
            stats.append(0.0); pvals.append(1.0); continue
        s, p = ranksums(x_hc, x_ot)
        stats.append(s); pvals.append(p)

    de_df = pd.DataFrame({
        'gene': genes, 'stat': stats, 'pval': pvals,
        'mean_high_contrib': mean_hc, 'mean_other': mean_other,
    })
    de_df['log2FC'] = np.log2((de_df['mean_high_contrib'] + 1e-6) / (de_df['mean_other'] + 1e-6))
    _, de_df['padj'], _, _ = multipletests(de_df['pval'], method='fdr_bh')
    de_df = de_df.sort_values('stat', ascending=False)

    n_up = ((de_df['padj'] < 0.05) & (de_df['log2FC'] > 0.5)).sum()
    n_dn = ((de_df['padj'] < 0.05) & (de_df['log2FC'] < -0.5)).sum()
    log(f"  Significant DEGs (padj<0.05, |log2FC|>0.5): {n_up} up, {n_dn} down")

    de_path = os.path.join(out_dir, "Tier1C_high_contrib_vs_other_DE.tsv")
    de_df.to_csv(de_path, sep='\t', index=False)
    log(f"  Saved: {de_path}")

    # Volcano plot
    banner("VOLCANO PLOT")
    fig, ax = plt.subplots(figsize=(10, 8))
    sig_up = (de_df['padj'] < 0.05) & (de_df['log2FC'] > 0.5)
    sig_dn = (de_df['padj'] < 0.05) & (de_df['log2FC'] < -0.5)
    ns = ~sig_up & ~sig_dn
    ax.scatter(de_df.loc[ns, 'log2FC'], -np.log10(de_df.loc[ns, 'pval']+1e-300),
              c='#d4d4d4', s=5, alpha=0.3)
    ax.scatter(de_df.loc[sig_up, 'log2FC'], -np.log10(de_df.loc[sig_up, 'pval']+1e-300),
              c='#ed6a5a', s=10, alpha=0.5, label=f'Up in high-contrib ({sig_up.sum()})')
    ax.scatter(de_df.loc[sig_dn, 'log2FC'], -np.log10(de_df.loc[sig_dn, 'pval']+1e-300),
              c='#5e81ac', s=10, alpha=0.5, label=f'Down in high-contrib ({sig_dn.sum()})')

    # Label top genes
    top = pd.concat([de_df[sig_up].head(10), de_df[sig_dn].tail(10)])
    for _, r in top.iterrows():
        ax.annotate(r['gene'], (r['log2FC'], -np.log10(r['pval']+1e-300)),
                   fontsize=7, ha='center', va='bottom', xytext=(0,3), textcoords='offset points')

    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)
    ax.axvline(0.5, color='grey', linestyle='--', linewidth=0.5)
    ax.axvline(-0.5, color='grey', linestyle='--', linewidth=0.5)
    ax.set_xlabel('log2 Fold Change (high-contrib / other)', fontsize=12)
    ax.set_ylabel('-log10(p-value)', fontsize=12)
    ax.set_title(f'High-Contributor (SC029+SC013+SC001) vs Other Patients\nAll basal cells',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(os.path.join(out_dir, f"Tier1C_volcano_high_contrib_vs_other.{ext}"),
                   dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved volcano")

    # GSEA prerank
    banner("GSEA PRERANK (KEGG)")
    rnk = de_df.set_index('gene')['stat'].dropna()
    rnk = rnk[~rnk.index.duplicated(keep='first')]
    try:
        res = gp.prerank(rnk=rnk, gene_sets='KEGG_2021_Human',
                         min_size=GSEA_MIN_SIZE, max_size=GSEA_MAX_SIZE,
                         permutation_num=GSEA_PERMUTATIONS, seed=42,
                         no_plot=True, verbose=False)
        gsea_df = res.res2d
        gsea_path = os.path.join(out_dir, "Tier1C_high_contrib_KEGG_GSEA.tsv")
        gsea_df.to_csv(gsea_path, sep='\t', index=False)
        n_sig = (gsea_df['FDR q-val'].astype(float) < 0.25).sum()
        log(f"  {len(gsea_df)} pathways, {n_sig} FDR<0.25")
        log(f"  Saved: {gsea_path}")

        # Print top enriched and depleted
        gsea_df['NES'] = gsea_df['NES'].astype(float)
        gsea_df['FDR q-val'] = gsea_df['FDR q-val'].astype(float)
        top_up = gsea_df.sort_values('NES', ascending=False).head(10)
        top_dn = gsea_df.sort_values('NES', ascending=True).head(10)
        log(f"\n  Top 10 enriched in high-contributor patients:")
        for _, r in top_up.iterrows():
            fdr_str = f"FDR={r['FDR q-val']:.3f}"
            log(f"    NES={r['NES']:+.2f}  {fdr_str}  {r['Term']}")
        log(f"\n  Top 10 depleted in high-contributor patients:")
        for _, r in top_dn.iterrows():
            fdr_str = f"FDR={r['FDR q-val']:.3f}"
            log(f"    NES={r['NES']:+.2f}  {fdr_str}  {r['Term']}")

    except Exception as e:
        log(f"  ERROR: GSEA failed: {e}")

    log("\nTier 1C complete.")

if __name__ == "__main__":
    main()
