#!/usr/bin/env python3
"""
Tier1A_Patient_Expression_GSEA.py
==================================

Tier 1A: Unbiased per-patient expression analysis with KEGG GSEA.

For each of the 14 patients:
  1. Compare that patient's basal cells vs all other patients' basal cells
  2. Wilcoxon rank-sum per gene -> full ranked gene list
  3. GSEA prerank against KEGG pathways (unbiased)
  4. Track per-patient A3A/A3B expression as focused overlay

Output (-> data/FIG_5/01_patient_expression/):
  - patient_KEGG_NES_heatmap.pdf/.png
  - patient_A3_expression_summary.tsv
  - per_patient_GSEA_results/ (full tables)
  - patient_GSEA_top_pathways.tsv

Usage:
  conda run -n NETWORK python Tier1A_Patient_Expression_GSEA.py
"""

import os, sys, numpy as np, pandas as pd, scanpy as sc, scipy.sparse
from scipy.stats import ranksums
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import gseapy as gp

from patient_config import *

def compute_patient_rankings(adata_basal):
    """Wilcoxon rank-sum: each patient vs all others. Returns dict of DataFrames."""
    banner("COMPUTING PER-PATIENT GENE RANKINGS")
    X = adata_basal.X
    if scipy.sparse.issparse(X):
        log("  Converting sparse -> dense..."); X = X.toarray()
    gene_names = adata_basal.var_names.values
    labels = adata_basal.obs[PATIENT_COL].values
    patients = sorted(set(labels))
    rankings = {}
    for i, patient in enumerate(patients):
        log(f"  [{i+1}/{len(patients)}] {patient}...")
        m_in, m_out = labels == patient, labels != patient
        stats, pvals = [], []
        for g in range(X.shape[1]):
            if X[m_in, g].std() == 0 and X[m_out, g].std() == 0:
                stats.append(0.0); pvals.append(1.0); continue
            s, p = ranksums(X[m_in, g], X[m_out, g])
            stats.append(s); pvals.append(p)
        df = pd.DataFrame({'gene': gene_names, 'stat': stats, 'pval': pvals})
        df = df.sort_values('stat', ascending=False)
        rankings[patient] = df
        log(f"    Top: {df.iloc[0]['gene']} (stat={df.iloc[0]['stat']:.2f})")
    return rankings

def run_gsea_per_patient(rankings, out_dir):
    """GSEA prerank against KEGG for each patient."""
    banner("RUNNING GSEA PRERANK (KEGG)")
    gsea_dir = ensure_dir(os.path.join(out_dir, "per_patient_GSEA_results"))
    all_results = {}
    for patient, rank_df in rankings.items():
        short = patient.replace("Patient ", "")
        log(f"  {patient}...")
        rnk = rank_df.set_index('gene')['stat'].dropna()
        rnk = rnk[~rnk.index.duplicated(keep='first')]
        try:
            res = gp.prerank(rnk=rnk, gene_sets='KEGG_2021_Human',
                             min_size=GSEA_MIN_SIZE, max_size=GSEA_MAX_SIZE,
                             permutation_num=GSEA_PERMUTATIONS, seed=42,
                             no_plot=True, verbose=False)
            rdf = res.res2d; rdf['patient'] = patient
            rdf.to_csv(os.path.join(gsea_dir, f"{short}_KEGG_GSEA.tsv"), sep='\t', index=False)
            n_sig = (rdf['FDR q-val'].astype(float) < 0.25).sum()
            log(f"    {len(rdf)} pathways, {n_sig} FDR<0.25")
            all_results[patient] = rdf
        except Exception as e:
            log(f"    ERROR: {e}"); all_results[patient] = pd.DataFrame()
    return all_results

def compute_a3_expression(adata_basal):
    """Per-patient mean A3 and C0 interactor expression."""
    banner("PER-PATIENT A3 EXPRESSION")
    patients = sorted(adata_basal.obs[PATIENT_COL].unique())
    rows = []
    for patient in patients:
        mask = adata_basal.obs[PATIENT_COL] == patient
        subset = adata_basal[mask]
        row = {'patient': patient, 'n_basal': mask.sum()}
        for gene in A3_GENES_SYMBOLS + C0_INTERACTORS:
            expr = get_gene_expression(subset, gene)
            row[f'{gene}_mean'] = expr.mean() if expr is not None else 0.0
            row[f'{gene}_pct_nz'] = ((expr > 0).mean() * 100) if expr is not None else 0.0
        a3a, a3b = row.get('APOBEC3A_mean', 0), row.get('APOBEC3B_mean', 0)
        row['A3A_A3B_ratio'] = a3a / a3b if a3b > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)

def plot_gsea_heatmap(all_results, enrichment_df, out_dir):
    banner("GSEA NES HEATMAP")
    nes_dict = {}
    for patient, rdf in all_results.items():
        if len(rdf) == 0: continue
        for _, row in rdf.iterrows():
            nes_dict.setdefault(row['Term'], {})[patient] = float(row['NES'])
    if not nes_dict: log("  No results to plot"); return
    nes_df = pd.DataFrame(nes_dict).T.fillna(0)
    if enrichment_df is not None:
        order = enrichment_df.sort_values('fold_enrichment', ascending=False)['patient'].tolist()
        order = [p for p in order if p in nes_df.columns]
        for p in nes_df.columns:
            if p not in order: order.append(p)
        nes_df = nes_df[order]
    top30 = nes_df.abs().max(axis=1).sort_values(ascending=False).head(30).index
    nes_plot = nes_df.loc[top30]
    fig, ax = plt.subplots(figsize=(16, 12))
    vmax = max(abs(nes_plot.values.min()), abs(nes_plot.values.max()), 2.0)
    im = ax.imshow(nes_plot.values, cmap='RdBu_r', aspect='auto', vmin=-vmax, vmax=vmax)
    labs = [p.replace('Patient ', '') for p in nes_plot.columns]
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=45, ha='right', fontsize=11)
    pnames = [p[:47]+'...' if len(p)>50 else p.replace('_',' ') for p in nes_plot.index]
    ax.set_yticks(range(len(pnames))); ax.set_yticklabels(pnames, fontsize=9)
    for i, p in enumerate(nes_plot.columns):
        if p in HIGH_CONTRIBUTORS:
            ax.get_xticklabels()[i].set_color('#ed6a5a'); ax.get_xticklabels()[i].set_fontweight('bold')
    plt.colorbar(im, ax=ax, label='NES', shrink=0.6)
    ax.set_title('Per-Patient KEGG Pathway Enrichment (GSEA)\nSorted by SBS2-HIGH fold enrichment',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(os.path.join(out_dir, f"patient_KEGG_NES_heatmap.{ext}"), dpi=300, bbox_inches='tight')
    plt.close()
    log("  Saved heatmap")

def plot_a3_vs_enrichment(a3_df, enrichment_df, out_dir):
    banner("A3 EXPRESSION VS ENRICHMENT")
    if enrichment_df is None: log("  No enrichment data"); return
    m = a3_df.merge(enrichment_df[['patient','fold_enrichment','n_high']], on='patient', how='left')
    m['fold_enrichment'] = m['fold_enrichment'].fillna(0)
    cols = ['APOBEC3A_mean','APOBEC3B_mean','A3A_A3B_ratio','APOBEC3H_mean']
    titles = ['Mean A3A','Mean A3B','A3A/A3B Ratio','Mean A3H']
    fig, axes = plt.subplots(1, 4, figsize=(24, 6))
    for ax, col, title in zip(axes, cols, titles):
        if col not in m.columns: ax.set_visible(False); continue
        colors = ['#ed6a5a' if p in HIGH_CONTRIBUTORS else '#5e81ac' for p in m['patient']]
        ax.scatter(m[col], m['fold_enrichment'], c=colors, s=100, edgecolors='black', linewidth=0.8)
        for _, r in m.iterrows():
            ax.annotate(r['patient'].replace('Patient ',''), (r[col], r['fold_enrichment']),
                       fontsize=8, ha='center', va='bottom', xytext=(0,5), textcoords='offset points')
        ax.axhline(1.0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.set_xlabel(title, fontsize=12); ax.set_ylabel('Fold Enrichment', fontsize=12)
        ax.set_title(title, fontsize=13, fontweight='bold')
    plt.suptitle('Per-Patient A3 Expression vs SBS2-HIGH Enrichment', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(os.path.join(out_dir, f"patient_A3_vs_enrichment.{ext}"), dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved A3 scatter")

def main():
    banner("TIER 1A: PER-PATIENT EXPRESSION GSEA")
    out_dir = ensure_dir(DIR_01_EXPRESSION)
    adata = load_adata()
    basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
    log(f"  Basal cells: {basal.n_obs:,}")

    enrich_path = os.path.join(DIR_00_DIAG, "patient_enrichment_SBS2_HIGH_v2.tsv")
    enrichment_df = pd.read_csv(enrich_path, sep='\t') if os.path.exists(enrich_path) else None

    rankings = compute_patient_rankings(basal)
    all_results = run_gsea_per_patient(rankings, out_dir)

    combined = pd.concat([d for d in all_results.values() if len(d)>0], ignore_index=True)
    combined.to_csv(os.path.join(out_dir, "patient_GSEA_top_pathways.tsv"), sep='\t', index=False)

    a3_df = compute_a3_expression(basal)
    a3_df.to_csv(os.path.join(out_dir, "patient_A3_expression_summary.tsv"), sep='\t', index=False)

    plot_gsea_heatmap(all_results, enrichment_df, out_dir)
    plot_a3_vs_enrichment(a3_df, enrichment_df, out_dir)
    log("\nTier 1A complete.")

if __name__ == "__main__":
    main()
