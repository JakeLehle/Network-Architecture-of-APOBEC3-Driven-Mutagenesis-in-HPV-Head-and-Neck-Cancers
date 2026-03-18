#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 12:44:15 2026

@author: jlehle
"""
#%% Read in the SRAscraper config.yaml to get the working dir
# This section is the only thing that need to be automated to bring it into the snakmake pipeline
import os
import yaml
# This is the one line you would have to change vvvvvvv
os.chdir('/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer')

with open('config_NMF_v0.1.1.yaml', 'r') as config:
    config_yaml = yaml.safe_load(config)


#%% Go load the meatdata of all the files downloaded
working_dir = config_yaml['output_dir']
figures_dir = os.path.join(working_dir, 'figures')
os.makedirs(figures_dir, exist_ok=True)
os.chdir(working_dir)
os.getcwd()
import pickle
with open(os.path.join(working_dir, 'metadata/dictionary_file_filtered.pkl'), 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')

#%%#%% load back in with all datasets
import scanpy as sc
import numpy as np
adata_pp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_pp_CD.h5ad')
)
#%%
import pandas as pd

# Run differential expression to find marker genes for each cell type
# Using Wilcoxon rank-sum test (robust for single-cell data)
sc.tl.rank_genes_groups(
    adata_pp,
    groupby='final_annotation',
    method='wilcoxon',
    n_genes=20
)

# Extract top 20 marker genes for each cell type into a DataFrame
cell_types = adata_pp.obs['final_annotation'].unique()
marker_dict = {}

for ct in cell_types:
    # Get ranked gene names for this cell type
    genes = sc.get.rank_genes_groups_df(adata_pp, group=ct)['names'].head(20).tolist()
    marker_dict[ct] = genes

# Create DataFrame: columns = cell types, rows = rank (1-20)
marker_df = pd.DataFrame(marker_dict)
marker_df.index = range(1, 21)
marker_df.index.name = 'rank'

print(marker_df)
#%%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# Define selected marker genes (common symbols only)
cell_type_markers = {
    'CD4-Positive alpha-beta T cell': ['IL7R', 'CD3E'],
    'B Cells': ['CD79A', 'MS4A1'],
    'fibroblasts': ['DCN', 'COL1A1'],
    'CD8-Positive, alpha-beta T cell': ['CD8A', 'CCL5'],
    'basal cell': ['KRT5', 'TACSTD2'],
    'macrophage': ['CD68', 'SPI1'],
    'endothelial cell': ['PECAM1', 'VWF'],
    'smooth muscle cell': ['TAGLN', 'RGS5'],
    'regulatory T cell': ['TIGIT', 'CTLA4'],
    'mast cell': ['TPSAB1', 'CPA3'],
    'plasmacytoid dendritic cell': ['IL3RA', 'LILRA4'],
    'myeloid dendritic cell': ['LAMP3', 'CCR7']
}

# --- Part 1: Create marker gene dataframe and save as TSV ---
marker_df = pd.DataFrame(cell_type_markers)
marker_df.index = ['Marker_1', 'Marker_2']
marker_df.index.name = 'marker_rank'

# Save to TSV
marker_df.to_csv('selected_marker_genes.tsv', sep='\t')
print("Saved marker gene table to 'selected_marker_genes.tsv'")
print(marker_df)

# --- Part 2: Create UMAP figure with custom colormap ---

# Custom colormap: slate gray to bright red
slate_gray = '#708090'
bright_red = '#FF2400'  # Scarlet red - vibrant
cmap_custom = LinearSegmentedColormap.from_list(
    'slate_to_red', 
    [slate_gray, '#A0522D', '#CD5C5C', '#FF4500', bright_red]
)

# Organize cell types into columns (2 per column, 6 columns total)
cell_type_order = [
    'CD4-Positive alpha-beta T cell', 'B Cells',                    # Column 0
    'CD8-Positive, alpha-beta T cell', 'regulatory T cell',         # Column 1
    'macrophage', 'myeloid dendritic cell',                         # Column 2
    'plasmacytoid dendritic cell', 'mast cell',                     # Column 3
    'fibroblasts', 'smooth muscle cell',                            # Column 4
    'basal cell', 'endothelial cell'                                # Column 5
]
#%%
# Create figure: 6 columns x 4 rows
fig, axes = plt.subplots(4, 6, figsize=(24, 16))
plt.subplots_adjust(wspace=0.1, hspace=0.35)  # Increased hspace for larger titles

# Plot each marker gene on UMAP
for col_idx in range(6):
    # Two cell types per column
    ct1 = cell_type_order[col_idx * 2]
    ct2 = cell_type_order[col_idx * 2 + 1]
    
    # Cell type 1: rows 0 and 1
    for gene_idx, gene in enumerate(cell_type_markers[ct1]):
        ax = axes[gene_idx, col_idx]
        sc.pl.umap(
            adata_pp,
            color=gene,
            cmap=cmap_custom,
            ax=ax,
            show=False,
            title=f'{ct1}\n{gene}',
            colorbar_loc='right' if col_idx == 5 else None,
            frameon=False,
            size=10,
            na_color=slate_gray
        )
        ax.set_title(f'{gene}', fontsize=38, fontweight='bold')  # Individual UMAP title
        ax.set_xlabel('')
        ax.set_ylabel('')
    
    # Cell type 2: rows 2 and 3
    for gene_idx, gene in enumerate(cell_type_markers[ct2]):
        ax = axes[gene_idx + 2, col_idx]
        sc.pl.umap(
            adata_pp,
            color=gene,
            cmap=cmap_custom,
            ax=ax,
            show=False,
            title=f'{ct2}\n{gene}',
            colorbar_loc='right' if col_idx == 5 else None,
            frameon=False,
            size=10,
            na_color=slate_gray
        )
        ax.set_title(f'{gene}', fontsize=38, fontweight='bold')  # Individual UMAP title
        ax.set_xlabel('')
        ax.set_ylabel('')

# Add overall title
fig.suptitle('Marker Gene Expression Validating popV Cell Type Annotations', 
             fontsize=26, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('supplemental_figure_marker_validation.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('supplemental_figure_marker_validation.pdf', bbox_inches='tight', facecolor='white')
plt.show()