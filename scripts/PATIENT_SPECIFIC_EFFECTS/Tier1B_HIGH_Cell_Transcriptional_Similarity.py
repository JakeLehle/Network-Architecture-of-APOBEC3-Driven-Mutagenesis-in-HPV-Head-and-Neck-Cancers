#!/usr/bin/env python3
"""
Tier1B_HIGH_Cell_Transcriptional_Similarity.py (v2 — fixed DE file detection)
===============================================================================

Among the 546 SBS2-HIGH basal cells:
  1. PCA on DE genes from Figure 4 Step 01 (SC_selected_genes_filtered.csv)
  2. Color by patient — do cells intermingle or segregate?
  3. Silhouette score (patient as label)
  4. Per-patient mean HIGH-cell profile correlations

Usage: conda run -n NETWORK python Tier1B_HIGH_Cell_Transcriptional_Similarity.py
"""

import os, sys, numpy as np, pandas as pd, scanpy as sc, scipy.sparse
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.preprocessing import StandardScaler
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

from patient_config import *

def load_de_genes():
    """Load DE genes from Figure 4 Step 01, trying multiple file patterns."""
    # Priority order: filtered CSV, selected CSV, stats CSV, then any TSV
    candidates = [
        ('SC_selected_genes_filtered.csv', ','),
        ('SC_selected_genes.csv', ','),
        ('SC_diffexpr_stats.csv', ','),
    ]
    for fname, sep in candidates:
        fpath = os.path.join(DE_GENES_DIR, fname)
        if not os.path.exists(fpath):
            continue
        log(f"  Found DE file: {fname}")
        df = pd.read_csv(fpath, sep=sep)
        log(f"    Columns: {list(df.columns)}, Rows: {len(df)}")
        # Find gene column
        for col in ['gene', 'gene_id', 'gene_ids', 'Gene']:
            if col in df.columns:
                genes = df[col].tolist()
                # For the full stats file, filter to significant
                if 'diffexpr_stats' in fname and 'pval' in df.columns:
                    genes = df.loc[df['pval'] < 0.05, col].tolist()
                    log(f"    Filtered stats to p<0.05: {len(genes)} genes")
                else:
                    log(f"    Using column '{col}': {len(genes)} genes")
                return genes
        # Try first column if it looks like gene names
        first_col = df.columns[0]
        samples = df[first_col].head(5).tolist()
        if all(isinstance(v, str) for v in samples):
            log(f"    Using first column '{first_col}': {len(df)} genes")
            return df[first_col].tolist()

    # Try TSV files
    for f in os.listdir(DE_GENES_DIR):
        if f.endswith('.tsv'):
            df = pd.read_csv(os.path.join(DE_GENES_DIR, f), sep='\t')
            log(f"  Fallback TSV: {f}, {len(df)} rows")
            return df.iloc[:, 0].tolist()

    log(f"  WARNING: No DE file found in {DE_GENES_DIR}")
    log(f"  Files present: {os.listdir(DE_GENES_DIR)}")
    return None

def main():
    banner("TIER 1B: HIGH-CELL TRANSCRIPTIONAL SIMILARITY (v2)")
    out_dir = ensure_dir(DIR_01_EXPRESSION)
    adata = load_adata()
    high_cells, low_cells = load_groups()
    high_in = high_cells & set(adata.obs_names)
    adata_high = adata[list(high_in)].copy()
    log(f"  HIGH cells: {adata_high.n_obs:,}")

    patient_labels = adata_high.obs[PATIENT_COL].values
    unique_patients = sorted(set(patient_labels))
    for p in unique_patients:
        log(f"    {p}: {(patient_labels==p).sum()} cells")

    # Load DE genes
    de_genes = load_de_genes()
    if de_genes is None:
        log("  Falling back to HVGs")
        sc.pp.highly_variable_genes(adata_high, n_top_genes=5000, flavor='seurat_v3')
        de_genes = adata_high.var_names[adata_high.var['highly_variable']].tolist()

    de_in = [g for g in de_genes if g in adata_high.var_names]
    log(f"  DE genes matched in adata: {len(de_in)}")
    adata_de = adata_high[:, de_in]
    X = adata_de.X
    if scipy.sparse.issparse(X): X = X.toarray()

    # PCA
    banner("PCA")
    X_s = StandardScaler().fit_transform(X)
    n_comp = min(50, X_s.shape[1], X_s.shape[0]-1)
    pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
    X_pca = pca.fit_transform(X_s)
    var = pca.explained_variance_ratio_
    log(f"  PC1: {var[0]*100:.1f}%, PC2: {var[1]*100:.1f}%, PC1-10: {sum(var[:10])*100:.1f}%")

    cmap = plt.cm.tab20
    pcols = {p: cmap(i/max(len(unique_patients)-1,1)) for i,p in enumerate(sorted(unique_patients))}
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    for ax_idx, pcy_i in enumerate([1, 2]):
        ax = axes[ax_idx]
        if pcy_i >= X_pca.shape[1]: pcy_i = 1
        for p in unique_patients:
            m = patient_labels == p
            mk = '*' if p in HIGH_CONTRIBUTORS else 'o'
            sz = 60 if p in HIGH_CONTRIBUTORS else 30
            ax.scatter(X_pca[m,0], X_pca[m,pcy_i], c=[pcols[p]], s=sz, marker=mk,
                      alpha=0.6, edgecolors='black', linewidth=0.3,
                      label=p.replace('Patient ',''))
        ax.set_xlabel(f'PC1 ({var[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC{pcy_i+1} ({var[pcy_i]*100:.1f}%)', fontsize=12)
        ax.set_title(f'HIGH Cells: PC1 vs PC{pcy_i+1}', fontsize=13, fontweight='bold')
        if ax_idx == 0: ax.legend(fontsize=8, ncol=2, loc='best')
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(os.path.join(out_dir, f"HIGH_cell_PCA_by_patient.{ext}"), dpi=300, bbox_inches='tight')
    plt.close(); log("  Saved PCA plots")

    # Silhouette
    banner("SILHOUETTE SCORE")
    pts_ok = [p for p in unique_patients if (patient_labels==p).sum() >= 3]
    if len(pts_ok) >= 2:
        m_ok = np.isin(patient_labels, pts_ok)
        sil = silhouette_score(X_pca[m_ok,:10], patient_labels[m_ok])
        sil_per = silhouette_samples(X_pca[m_ok,:10], patient_labels[m_ok])
        log(f"  Overall silhouette: {sil:.4f}")
        interp = ("VERY LOW - shared program" if sil < 0.1 else
                  "LOW - mostly shared" if sil < 0.25 else
                  "MODERATE - partial segregation" if sil < 0.5 else
                  "HIGH - patient-specific programs")
        log(f"  Interpretation: {interp}")
        for p in sorted(pts_ok):
            mp = patient_labels[m_ok] == p
            log(f"    {p}: {sil_per[mp].mean():.4f} (n={mp.sum()})")
        with open(os.path.join(out_dir, "HIGH_cell_silhouette_score.txt"), 'w') as f:
            f.write(f"Overall silhouette: {sil:.4f}\n{interp}\n")
            f.write(f"DE genes used: {len(de_in)}\n\n")
            for p in sorted(pts_ok):
                mp = patient_labels[m_ok] == p
                f.write(f"{p}: {sil_per[mp].mean():.4f} (n={mp.sum()})\n")

    # Per-patient mean profile correlations
    banner("MEAN PROFILE CORRELATIONS")
    pts_corr = [p for p in unique_patients if (patient_labels==p).sum() >= 2]
    profiles = {p: X[patient_labels==p].mean(axis=0) for p in pts_corr}
    prof_df = pd.DataFrame(profiles)
    corr = prof_df.corr(method='spearman')
    off_diag = corr.values[np.triu_indices_from(corr.values, k=1)]
    log(f"  Mean off-diag rho: {off_diag.mean():.4f}, min: {off_diag.min():.4f}, max: {off_diag.max():.4f}")
    prof_df.to_csv(os.path.join(out_dir, "HIGH_cell_patient_mean_profiles.tsv"), sep='\t')

    fig, ax = plt.subplots(figsize=(10, 8))
    labs = [p.replace('Patient ','') for p in corr.columns]
    im = ax.imshow(corr.values, cmap='RdYlBu_r', vmin=0, vmax=1)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=45, ha='right', fontsize=11)
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=11)
    for i, p in enumerate(corr.columns):
        if p in HIGH_CONTRIBUTORS:
            ax.get_xticklabels()[i].set_color('#ed6a5a'); ax.get_xticklabels()[i].set_fontweight('bold')
            ax.get_yticklabels()[i].set_color('#ed6a5a'); ax.get_yticklabels()[i].set_fontweight('bold')
    for i in range(len(labs)):
        for j in range(len(labs)):
            ax.text(j, i, f'{corr.values[i,j]:.2f}', ha='center', va='center', fontsize=8,
                   color='white' if corr.values[i,j]<0.5 else 'black')
    plt.colorbar(im, ax=ax, label='Spearman rho', shrink=0.8)
    ax.set_title('Mean HIGH-Cell Profile Correlation (per patient)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(os.path.join(out_dir, f"HIGH_cell_patient_correlation_heatmap.{ext}"), dpi=300, bbox_inches='tight')
    plt.close()
    log("\nTier 1B complete.")

if __name__ == "__main__":
    main()
