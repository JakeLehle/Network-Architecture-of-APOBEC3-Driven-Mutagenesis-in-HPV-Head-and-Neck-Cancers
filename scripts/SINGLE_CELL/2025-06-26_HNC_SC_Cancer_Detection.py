#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 09:16:17 2025

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
     
#%% Convert dict to dataframe
import pandas as pd

x_data = pd.DataFrame()
for key in gse_dict.keys():
    gse_dict[key]['series_id'] = key
    x_data = pd.concat([x_data, gse_dict[key]])
x_data.reset_index(inplace=True)
x_data = x_data.drop(columns=['index'])
x_data['run_alias'] = x_data['run_alias'].str.split('_').str[0]
#%% display the automated cell types

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
import matplotlib.patheffects as pe
from adjustText import adjust_text
import matplotlib.pyplot as plt

def gen_mpl_labels(
    adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None, color_by_group=False
):
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    # Fill the text colors dictionary
    text_colors = {group: None for group in adata.obs[groupby].cat.categories}

    if color_by_group and groupby + "_colors" in adata.uns:
        for i, group in enumerate(adata.obs[groupby].cat.categories):
            if group in exclude:
                continue
            text_colors[group] = adata.uns[groupby + "_colors"][i]

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)


def nonoverlapping_UMAP(adata_obj, group_name):
    cmap = plt.get_cmap('turbo')
    # Generate evenly spaced values equal to the value of groups in the dataset
    adata_obj.obs[group_name] = adata_obj.obs[group_name].cat.remove_unused_categories()
    value_cat = pd.Categorical(adata_obj.obs[group_name])
    values = np.linspace(0, 1, len(value_cat.categories))
    # Get RGB values for each value in the colormap
    palette = [cmap(value) for value in values]
    
    # Combined path effects - white border first, then black border
    combined_effects = [
        pe.withStroke(linewidth=6, foreground="white"),  # Thick white border (outer)
        pe.withStroke(linewidth=1, foreground="black"),  # Thin black border (inner)
        pe.Normal()  # Original text color
    ]
    
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 150, "figure.frameon": False}):
        ax = sc.pl.umap(adata_obj, color=group_name, show=False, legend_loc=None, frameon=False, title='', size=5, palette=palette)
        gen_mpl_labels(
            adata_obj,
            group_name,
            exclude=("None",),
            ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(
                fontsize=26, 
                path_effects=combined_effects  # Use the combined effects here
            ),
            color_by_group=True
        )
        fig = ax.get_figure()
        fig.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'Final_Cell_Type_Annotation_UMAP.pdf'))
        plt.show()
        plt.close()
#%% load back in with all datasets
import scanpy as sc
adata = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata.h5ad')
)

adata_tmp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_tmp.h5ad')
)

adata_tmp_an = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_popv_an.h5ad')
)

adata_pp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_pp.h5ad')
)

#%% cytotrace2-py
# to install
# !git clone https://github.com/digitalcytometry/cytotrace2
# !cd cytotrace2/cytotrace2_python
# !pip install .
#%% Run CytoTRACE2-py
import pandas as pd
import scipy.sparse
import os
import time
import numpy as np
import shutil
from cytotrace2_py.cytotrace2_py import *

# Configuration - Adjust these parameters as needed
MAX_CELLS_PER_CHUNK = 200000  # Start with this, reduce if chunks fail
RANDOM_SEED = 42              # For reproducibility

# Initialize results tracking
all_cytotrace_results = []
project_ids = adata_tmp.obs['series_id'].unique()

# Define output directories
chunk_results_dir = os.path.join(working_dir, "cytotrace_chunk_results")
final_results_dir = os.path.join(working_dir, "final_results")
cytotrace2_results_dir = os.path.join(working_dir, "cytotrace2_results")

def safe_makedirs(dir_path):
    """Create directory, removing existing one if it exists"""
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)

# Create output directories with cleanup
safe_makedirs(chunk_results_dir)
safe_makedirs(final_results_dir)

def process_chunk(adata_chunk, chunk_name):
    """Process a single chunk of cells with strict error handling"""
    # Prepare cell annotations
    adata_chunk.obs['cell_barcodes'] = adata_chunk.obs.index
    adata_annotation_df = adata_chunk.obs[['cell_barcodes', 'final_annotation']]
    
    # Save annotations to temporary file
    annotation_path = os.path.join(working_dir, f"cell_annotations_{chunk_name}.txt")
    adata_annotation_df.to_csv(
        annotation_path,
        sep="\t",
        header=True,
        index=False,
    )
    
    # Check for duplicates
    assert adata_chunk.var_names.is_unique, "Duplicate gene names detected!"
    assert adata_chunk.obs_names.is_unique, "Duplicate cell barcodes detected!"
    
    # Prepare and save gene expression matrix
    adata_X_df_T = pd.DataFrame.sparse.from_spmatrix(
        scipy.sparse.csr_matrix(adata_chunk.X),
        index=adata_chunk.obs.index,
        columns=adata_chunk.var['gene_symbol'].index
    ).T
    
    expression_path = os.path.join(working_dir, f"gene_expression_matrix_{chunk_name}.txt")
    adata_X_df_T.to_csv(
        expression_path,
        sep="\t",
        header=True,
        index=True,
        chunksize=10000
    )
    
    # Make sure to remove any existitng cytotrace2 output dir
    if os.path.exists(cytotrace2_results_dir):
        shutil.rmtree(cytotrace2_results_dir)

    # Run CytoTRACE2
    results = cytotrace2(
        expression_path,
        annotation_path=annotation_path,
        species="human",
        max_cores=os.cpu_count(),
        seed=RANDOM_SEED
    )
    
    # Verify results file exists
    project_results_path = os.path.join(working_dir, "cytotrace2_results", "cytotrace2_results.txt")
    if not os.path.exists(project_results_path):
        raise FileNotFoundError(f"CytoTRACE2 output missing at {project_results_path}")
    
    # Read results
    cytotrace_txt = pd.read_csv(project_results_path, sep='\t')
    
    # Clean up temporary files
    if os.path.exists(annotation_path):
        os.remove(annotation_path)
    if os.path.exists(expression_path):
        os.remove(expression_path)    
    
    return cytotrace_txt

# Main processing loop with proper cell selection
try:
    for project_id in project_ids:
        print(f"\n{'='*50}")
        print(f"Processing project: {project_id}")
        project_start = time.time()
        
        # Subset the data for this project - creates a view
        project_mask = adata_tmp.obs['series_id'] == project_id
        adata_project = adata_tmp[project_mask].copy()
        n_cells = adata_project.n_obs
        print(f"Total cells in project: {n_cells}")
        
        # Get cell barcodes (index values) as an array
        cell_barcodes = adata_project.obs.index.values.copy()
        
        # Create a shuffled array of positions (not the barcodes themselves)
        shuffled_positions = np.arange(n_cells)
        np.random.shuffle(shuffled_positions)
        
        # Determine number of chunks needed
        n_chunks = max(1, (n_cells + MAX_CELLS_PER_CHUNK - 1) // MAX_CELLS_PER_CHUNK)
        print(f"Splitting into {n_chunks} chunks")

        # Process each chunk
        for chunk_idx in range(n_chunks):
            chunk_start = time.time()
            chunk_name = f"{project_id}_chunk{chunk_idx+1}"
            print(f"\nProcessing {chunk_name}...")
            
            # Get positions for this chunk
            start_idx = chunk_idx * MAX_CELLS_PER_CHUNK
            end_idx = min((chunk_idx + 1) * MAX_CELLS_PER_CHUNK, n_cells)
            chunk_positions = shuffled_positions[start_idx:end_idx]
            
            # Get the actual cell barcodes using the shuffled positions
            chunk_cell_barcodes = cell_barcodes[chunk_positions]
            
            # Subset the data using the actual barcodes
            adata_chunk = adata_project[adata_project.obs.index.isin(chunk_cell_barcodes)].copy()
            
            # Verify we got the correct number of cells
            assert adata_chunk.n_obs == len(chunk_cell_barcodes), "Cell count mismatch in chunk"
            
            # Process chunk
            chunk_results = process_chunk(adata_chunk, chunk_name)
            chunk_results = chunk_results.rename(columns={'Unnamed: 0': 'cell_barcodes'})
            
            # Verify results contain the same cells we processed
            if not all(cell in chunk_results['cell_barcodes'].values for cell in chunk_cell_barcodes):
                raise ValueError("Processed cells don't match input cells in chunk")
            
            # Add project and chunk info
            chunk_results['series_id'] = project_id
            chunk_results['chunk_id'] = chunk_name
            all_cytotrace_results.append(chunk_results)
            
            # Move results to chunk-specific directory
            chunk_output_dir = os.path.join(chunk_results_dir, chunk_name)
            safe_makedirs(chunk_output_dir)
            
            unique_dir_name = f"cytotrace2_results_{chunk_name}"
            os.rename(
                cytotrace2_results_dir,
                os.path.join(chunk_output_dir, unique_dir_name)
            )
            
            print(f"Completed {chunk_name} in {time.time()-chunk_start:.2f} seconds")

    # Combine all results
    final_cytotrace_results = pd.concat(all_cytotrace_results, axis=0)
    
    # Merge with original annotations
    Cell_type_anno_df = pd.DataFrame(adata_tmp.obs[['final_annotation']])
    Cell_type_anno_df['cell_barcodes'] = Cell_type_anno_df.index

    final_cytotrace_results = pd.merge(
        final_cytotrace_results,
        Cell_type_anno_df,
        on=['cell_barcodes'],
        how='outer'
    )

    # Save final results
    final_output_path = os.path.join(working_dir, "final_results", "combined_cytotrace_results.txt")
    final_cytotrace_results.to_csv(
        final_output_path,
        sep="\t",
        header=True,
        index=False
    )
    
    print(f"\nSUCCESS: Processed all {len(all_cytotrace_results)} chunks across {len(project_ids)} projects")
    print(f"Final results saved to: {final_output_path}")

except Exception as e:
    print(f"\nERROR: Processing failed at chunk {chunk_name if 'chunk_name' in locals() else 'unknown'}")
    print(f"Error details: {str(e)}")
    raise

#%% Re load back in the cytotrace dataset

final_cytotrace_results = pd.read_csv(os.path.join(working_dir, 'final_results/combined_cytotrace_results.txt'), sep='\t')

#%% Histogram of CytoTRACE2 scores broken down by cell type

import seaborn as sns
sns.set_context("paper", rc={'xtick.labelsize': 12, 'ytick.labelsize': 10, "axes.labelsize":14}) # Set specific sizes
ax = sns.histplot(data=final_cytotrace_results, x='CytoTRACE2_Score', bins=100, hue='final_annotation', palette='tab20', multiple='stack', edgecolor='darkgrey')
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1.05), ncol=1, title="Cell Types", title_fontsize=14, fontsize=11)
plt.savefig(os.path.join(figures_dir, 'cytotrace_score_distribution.pdf'), bbox_inches='tight')
plt.show()

#%% Smoothing histogram peaks from CytoTRACE
import numpy as np
from scipy.signal import convolve
from scipy.stats import norm
from scipy.signal import find_peaks

def convolve_hist_for_peaks(data, bins=100, density=False, sigma = 3):
    counts, edges = np.histogram(data, bins=bins, density=density)
    counts_data = np.asarray(counts)
    kernel_size = int(2 * np.ceil(3 * sigma) + .1)
    
    # Create Gaussian kernel
    x = np.linspace(-2.5*sigma, 2.5*sigma, kernel_size)
    kernel = norm.pdf(x, loc=0, scale=sigma)
    kernel = kernel / np.sum(kernel)
    convolved = convolve(counts_data, kernel, mode='same')
    peaks, _ = find_peaks(convolved)

    # This array will match the values form the cytotrace scores from 0 to 1
    arr = np.linspace(0, 1, 100)
    # Do some calculus to find the minimum point between the first two peaks
    valley_idx = np.argmin(convolved[peaks[0]:peaks[1]]) + peaks[0]
    threshold = arr[valley_idx]

    # --- Plotting ---
    plt.figure(figsize=(14, 4))
    
    # Original Data Plot
    plt.subplot(1, 2, 1)
    plt.stem(arr, counts_data, linefmt='C0-', markerfmt='C0o', basefmt='C0-')
    plt.ylabel('Cell count', fontsize=18)
    plt.xlabel('CytoTrace2 Score', fontsize=18)
    plt.title('Original Data', fontsize=20)
    
    # Smoothed Data Plot
    ax = plt.subplot(1, 2, 2)
    stem = plt.stem(arr, convolved, linefmt='C1-', markerfmt='C1o', basefmt='C1-')
    
    # Add background shading (must be done FIRST for proper z-ordering)
    ax.axvspan(0, threshold, facecolor='lightblue', alpha=0.3, zorder=0)
    ax.axvspan(threshold, 1, facecolor='mistyrose', alpha=0.3, zorder=0)
    
    # Add threshold line
    plt.axvline(x=threshold, color='black', linestyle='--', 
                linewidth=2, label=f'Threshold: {threshold:.3f}')
    
    # Add peak markers
    plt.scatter(arr[peaks], convolved[peaks], color='black', 
                s=100, zorder=3, label='Peaks')
    
    # Add text labels
    plt.text(threshold/2, 0.5*max(convolved), 'Normal', 
             ha='center', va='center', fontsize=18, weight='bold')
    plt.text((threshold+1)/2, 0.5*max(convolved), 'Cancer', 
             ha='center', va='center', fontsize=18, weight='bold')
    
    plt.ylabel('Cell count (smoothed)', fontsize=18)
    plt.xlabel('CytoTrace2 Score', fontsize=18)
    plt.title('Smoothed Data with Population Separation', fontsize=18)
    plt.legend(loc='upper right', fontsize=18)
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'cytotrace_score_threshold_detection.pdf'), bbox_inches='tight')
    plt.show()
 
    return counts, convolved, peaks, threshold

bin_counts, bin_convolved, bin_peaks, infer_cnv_threshold = convolve_hist_for_peaks(final_cytotrace_results['CytoTRACE2_Score'])

#%% Merge the cytotrace dataframe with the adata_pp.obs
import scanpy as sc
adata_pp.obs['cell_barcodes'] = adata_pp.obs.index

cytotrace_txt = final_cytotrace_results.drop(columns=['final_annotation', 'series_id'])
cytotrace_txt = cytotrace_txt.set_index('cell_barcodes')
cytotrace_txt = cytotrace_txt.reindex(adata_pp.obs['cell_barcodes'])



adata_pp.obs = adata_pp.obs.join(cytotrace_txt, on='cell_barcodes')

sc.pl.umap(adata_pp, color="final_annotation", show=False)
plt.savefig(os.path.join(figures_dir, 'UMAP_final_annotation.pdf'), bbox_inches='tight')
plt.show()
#%% Plot the CytoTRACE2 Score
sc.pl.umap(
    adata_pp, 
    color="CytoTRACE2_Score", 
    title="CytoTRACE2 Cancer Risk Score", 
    show=False,
    frameon=False,
    color_map='plasma',
    size=5
)
#%% Make a new function to plot the combined dataset from CytoTRACE2
import seaborn as sns
import matplotlib

def plot_combined_boxplot(adata, figures_dir):
    """Final version with all elements properly spaced"""
    # Calculate median scores for strict descending order
    pheno_median_score = (adata.obs.groupby("final_annotation", observed=True)['CytoTRACE2_Score']
                          .median()
                          .sort_values(ascending=False))
    
    # Create ordered lists
    pheno_order = list(pheno_median_score.index)
    median_values = list(pheno_median_score.values)
    
    # Create color palette
    cytotrace2_pal = [plt.get_cmap('Spectral', 51)(round(max(0, min(5.5 - 6 * score, 5)) * 10)) 
                      for score in median_values]
    
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'Helvetica'
    
    # Create figure with adjusted layout
    fig, ax = plt.subplots(figsize=(14, 7))  # Increased width for better spacing
    plt.subplots_adjust(right=0.5, top=0.9)  # More space for elements
    
    # Create boxplot in correct order
    sns.boxplot(
        data=adata.obs,
        x='final_annotation',
        y='CytoTRACE2_Score',
        order=pheno_order,
        palette=cytotrace2_pal,
        linewidth=1,
        showfliers=False,
        ax=ax
    )
    
    # Add stripplot
    sns.stripplot(
        data=adata.obs,
        x='final_annotation',
        y='CytoTRACE2_Score',
        order=pheno_order,
        alpha=0.25,
        jitter=0.25,
        size=1,
        color='k',
        ax=ax
    )
    
    # Format main axes
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Potency Score', fontsize=18)
    ax.set_ylim(0, 1.0)
    ax.tick_params(axis='y', labelsize=14)
    
    # Add potency category grid lines
    for i in range(6):
        ax.axhline((i + 1)/6, linestyle='-', color='grey', linewidth=0.15, zorder=0)
    
    # Add right axis with potency categories - moved higher
    ax_right = ax.twinx()
    ax_right.grid(False)
    ax_right.set_ylim([0,1])
    ax_right.set_yticks(
        np.linspace(1/12, 11/12, 6),
        ['Differentiated','Unipotent','Oligopotent','Multipotent','Pluripotent','Totipotent'],
        fontsize=16
    )
    ax_right.tick_params(axis='y', left=False, right=False, labelright=True, pad=5)
    ax_right.spines['right'].set_visible(False)
    
    # Create numbered legend
    legend_entries = []
    for idx, (cell_type, score) in enumerate(zip(pheno_order, median_values), 1):
        legend_entries.append(
            plt.Line2D([0], [0],
                      marker='s',
                      color='w',
                      label=f"{idx}. {cell_type} ({score:.3f})",
                      markerfacecolor=cytotrace2_pal[idx-1],
                      markersize=16)
        )
    
    # Add legend with proper spacing
    leg = ax.legend(
        handles=legend_entries,
        title='Cell Types (Median Score)',
        bbox_to_anchor=(1.4, 1),  # Adjusted position
        loc='upper left',
        frameon=True,
        title_fontsize=20,
        fontsize=15,
        handletextpad=0.5,
        labelspacing=0.8,
        borderpad=1
    )
    
    # Add column numbers below plot
    for idx in range(len(pheno_order)):
        ax.text(idx, -0.08, str(idx+1), 
                ha='center', va='top', 
                fontsize=14)
    
    # Final title
    ax.set_title(" ", 
                pad=20, fontsize=14)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save with proper bbox
    plt.savefig(
        os.path.join(figures_dir, 'final_potency_score_by_celltype.pdf'),
        bbox_extra_artists=(leg, ax_right),
        bbox_inches='tight',
        dpi=300
    )
    plt.show()

def plot_combined_umaps(adata, figures_dir):
    """Modified version of CytoTRACE2's UMAP plots for our combined data"""
    # Create potency categories based on your threshold
    adata.obs['Potency_Category'] = pd.cut(
        adata.obs['CytoTRACE2_Score'],
        bins=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
        labels=['Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent']
    )
    
    # Set up figure
    sc.set_figure_params(scanpy=True, fontsize=18, figsize=(4,4))
    plt.rc('legend', fontsize=14)
    
    # Plot 1: UMAP colored by cell type
    fig1 = sc.pl.umap(
        adata,
        color='final_annotation',
        title='Cell Types',
        palette='tab20',
        return_fig=True
    )
    ax1 = fig1.axes[0]
    ax1.spines[['right', 'top']].set_visible(False)
    ax1.set_aspect(1 / ax1.get_data_ratio())
    plt.savefig(os.path.join(figures_dir, 'combined_celltype_umap.pdf'), bbox_inches='tight')
    plt.show()
    
    # Plot 2: UMAP colored by potency category
    potency_palette = {
        'Differentiated': "#5E4FA2",
        'Unipotent': "#66C2A5",
        'Oligopotent': "#E6F598",
        'Multipotent': "#FEE08B",
        'Pluripotent': "#F46D43"
    }
    
    fig2 = sc.pl.umap(
        adata,
        color='Potency_Category',
        palette=potency_palette,
        title='Potency Categories',
        return_fig=True
    )
    ax2 = fig2.axes[0]
    ax2.legend_.set_title("Potency category", fontsize=16)
    ax2.spines[['right', 'top']].set_visible(False)
    ax2.set_aspect(1 / ax2.get_data_ratio())
    plt.savefig(os.path.join(figures_dir, 'combined_potency_category_umap.pdf'), bbox_inches='tight')
    plt.show()
    
    # Plot 3: UMAP colored by continuous potency score (spectral colors)
    adata.obs['CytoTRACE2_Score_Plot'] = 6 * adata.obs['CytoTRACE2_Score'] - 5.5
    colors = ["#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"]
    pad_labels = ['','Differentiated','Unipotent','Oligopotent','Multipotent','Pluripotent']
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "potency_cols", list(zip(np.linspace(0,1,6), colors))
    )
    
    fig3 = sc.pl.umap(
        adata,
        color='CytoTRACE2_Score_Plot',
        cmap=cmap,
        title='Potency Score',
        return_fig=True,
        vmin=-5,
        vmax=0,
        colorbar_loc=None
    )
    ax3 = plt.gca()
    clb = plt.colorbar(ax3.get_children()[0], ax=ax3, fraction=0.035)
    clb.ax.set_title('Potency score', pad=10, loc='left')
    clb.ax.set_yticks(list(np.linspace(-5,0,6)))
    clb.ax.set_yticklabels(pad_labels, va='top')
    ax3.spines[['right', 'top']].set_visible(False)
    
    # Offset labels
    offset = matplotlib.transforms.ScaledTranslation(0, -1/12, fig3.dpi_scale_trans)
    for label in clb.ax.yaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    ax3.set_aspect(1 / ax3.get_data_ratio())
    
    plt.savefig(os.path.join(figures_dir, 'combined_potency_score_umap.pdf'), bbox_inches='tight')
    plt.show()

#%% Generate the plots
plot_combined_boxplot(adata_pp, figures_dir)
plot_combined_umaps(adata_pp, figures_dir)
nonoverlapping_UMAP(adata_pp, 'final_annotation')

#%% Add a new metadata column where you set the Cancer status based on the cytotrace potency scores to be used to train the model for InferCNV

adata_pp.obs['infer_cnv_threshold'] = np.where(adata_pp.obs['CytoTRACE2_Score'] > infer_cnv_threshold, 'Cancer', 'Normal')

#%% Save your work
adata_pp.write("adata_pp_CT.h5ad")

#%% Load back in
#adata_pp = sc.read_h5ad(
#    filename=os.path.join(working_dir, 'adata_pp_CT.h5ad')
#)

#%% Add chromosomal information for each gene to adata.var layer before running InferCNV
import gffutils
import os

# Define database path
db_path = os.path.join(working_dir, "genes.db")
#gtf_path set at top of script!

# Check if database exists, create if not
if not os.path.exists(db_path):
    print("Creating new gene database...")
    db = gffutils.create_db(
        gtf_path,
        dbfn=db_path,
        force=True,
        keep_order=True
    )
    print(f"Database created at: {db_path}")
else:
    print(f"Using existing database at: {db_path}")

# Load database (whether newly created or existing)
db = gffutils.FeatureDB(db_path)

# Extract gene info
gene_data = []
for gene_id in adata_pp.var.gene_ids:
    try:
        gene = db[gene_id]
        gene_data.append({
            'chromosome': gene.chrom,
            'start': gene.start,
            'end': gene.end
        })
    except:
        gene_data.append({
            'chromosome': None,
            'start': None,
            'end': None
        })

tmp = adata_pp.var
tmp
# Add to adata.var
gene_df = pd.DataFrame(gene_data, index=adata_pp.var.index)

adata_pp.var = pd.concat([adata_pp.var, gene_df], axis=1)

# Double check everything got added correctly
tmp = adata_pp.var

# Make sure the index match the gene_ids
adata_pp.var = adata_pp.var.set_index('gene_ids')

#%% identify cnv
import infercnvpy as cnv

cnv.tl.infercnv(
    adata_pp,
    reference_key="infer_cnv_threshold",
    reference_cat=["Normal"],
    window_size=250,
)

cnv.pl.chromosome_heatmap(adata_pp, groupby="infer_cnv_threshold", show=False)
plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_infer_cnv_threshold.pdf'), bbox_inches='tight')
plt.show()

cnv.pl.chromosome_heatmap(adata_pp, groupby="final_annotation", show=False)
plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_final_annotation.pdf'), bbox_inches='tight')
plt.show()
#%% Do clustering based on copy number variation
cnv.tl.pca(adata_pp)
cnv.pp.neighbors(adata_pp)
cnv.tl.leiden(adata_pp)
#sc.tl.dendrogram(adata_pp, groupby="cnv_leiden")
cnv.pl.chromosome_heatmap(adata_pp, groupby="cnv_leiden", dendrogram=True, show=False)
plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_cnv_leiden.pdf'), bbox_inches='tight')
plt.show()
cnv.tl.umap(adata_pp)
cnv.tl.cnv_score(adata_pp)
#%% calculating the CNV scores on their own UMAP
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata_pp,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(adata_pp, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata_pp, color="final_annotation", ax=ax3, show=False)
plt.savefig(os.path.join(figures_dir, 'cnv_analysis_panel.pdf'), bbox_inches='tight')
plt.show()
#%% Mapping CNV scores back on to orginal UMAP plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
ax4.axis("off")
sc.pl.umap(adata_pp, color="cnv_leiden", ax=ax1, show=False)
sc.pl.umap(adata_pp, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata_pp, color="final_annotation", ax=ax3, show=False)
plt.savefig(os.path.join(figures_dir, 'umap_analysis_panel.pdf'), bbox_inches='tight')
plt.show()
#%% This is if you wanted to let infercnv decide your cancer cells using the clustering
#adata.obs["cnv_status"] = "normal"
#adata.obs.loc[adata.obs["cnv_leiden"].isin(["10", "16", "13", "8", "12", "17", "1", "14", "11"]), "cnv_status"] = (
#    "tumor"
#)
#%% Side by side of cancer vs normal cells
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"wspace": 0.5})
cnv.pl.umap(adata_pp, color="infer_cnv_threshold", ax=ax1, show=False)
sc.pl.umap(adata_pp, color="infer_cnv_threshold", ax=ax2, show=False)
plt.savefig(os.path.join(figures_dir, 'infer_cnv_threshold_comparison.pdf'), bbox_inches='tight')
plt.show()
#%% Side by side of Cancer risk score with two models
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
ax4.axis("off")
sc.pl.umap(adata_pp, color="CytoTRACE2_Score", ax=ax1, show=False)
sc.pl.umap(adata_pp, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata_pp, color="final_annotation", ax=ax3, show=False)
plt.savefig(os.path.join(figures_dir, 'score_comparison_panel.pdf'), bbox_inches='tight')
plt.show()

#%% Output cnv chromosomal map of Cancer 
cnv.pl.chromosome_heatmap(adata_pp[adata_pp.obs["infer_cnv_threshold"] == "Cancer", :], show=False)
plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_cancer_cells.pdf'), bbox_inches='tight')
plt.show()
#%% Output cnv chromosomal map of Normal 
cnv.pl.chromosome_heatmap(adata_pp[adata_pp.obs["infer_cnv_threshold"] == "Normal", :], show=False)
plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_normal_cells.pdf'), bbox_inches='tight')
plt.show()
#%% Initial Spearman and Pearson plot
from scipy.stats import spearmanr

rho, p_value = spearmanr(adata_pp.obs['CytoTRACE2_Score'], adata_pp.obs['cnv_score'])
print(f"Spearman rho: {rho:.3f}")
print(f"Spearman p-value: {p_value:.3e}")

from scipy.stats import pearsonr
r, p_value = pearsonr(adata_pp.obs['CytoTRACE2_Score'], adata_pp.obs['cnv_score'])
print(f"Pearson r: {r}")
print(f"Pearson p-value: {p_value:.3e}")

#%% weight the agrrement between the two model scores
array1 = adata_pp.obs['CytoTRACE2_Score']
array2 = adata_pp.obs['cnv_score']
rank1 = np.argsort(np.argsort(array1))  # Rank each cell (0 = lowest risk)
rank2 = np.argsort(np.argsort(array2))
normalized_rank1 = rank1 / len(array1)  # Scale to [0, 1]
normalized_rank2 = rank2 / len(array2)

# Normalize scores to [0, 1] (if not already)
norm_array1 = (array1 - np.min(array1)) / (np.max(array1) - np.min(array1))
norm_array2 = (array2 - np.min(array2)) / (np.max(array2) - np.min(array2))

alpha = 0.5  # Weight between ranks and values (0.5 = both models have equal importance in weight assignment)
rank_agreement = 1 - np.abs(normalized_rank1 - normalized_rank2)
value_agreement = 1 - np.abs(norm_array1 - norm_array2)
agreement_score = (alpha * rank_agreement) + ((1 - alpha) * value_agreement)

#%% Colored correlation plot between the models before filtereing
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr

# --- Compute correlations (same as before) ---
spearman_rho, spearman_p = spearmanr(array1, array2)
pearson_r, pearson_p = pearsonr(array1, array2)

plt.figure(figsize=(4.5, 6))

# --- Scatter plot (your original styling) ---
sc = plt.scatter(
    array1, 
    array2, 
    c=agreement_score,
    cmap="RdYlGn",
    alpha=0.5, 
    s=5
)

# --- Add Seaborn-style regression line ---
# Fit linear regression (same as sns.regplot's default)
slope, intercept = np.polyfit(array1, array2, deg=1)  # deg=1 for linear
x_vals = np.linspace(array1.min(), array1.max(), 100)
y_fit = slope * x_vals + intercept

plt.plot(
    x_vals, 
    y_fit, 
    color='royalblue',          # Matches sns.regplot's red dashed line
    linestyle='--',       # Dashed line
    linewidth=1.5

)

# --- Add correlation stats (your original text box) ---
plt.text(
    0.05, 0.95, 
    f"Spearman ρ = {spearman_rho:.2f} \n(p = {spearman_p:.2e})\n\nPearson r = {pearson_r:.2f} \n(p = {pearson_p:.2e})",
    transform=plt.gca().transAxes,
    fontsize=18,
    ha='left', 
    va='top',
    bbox=dict(facecolor='white', alpha=0, edgecolor='gray')
)

# --- Log-scale axes (your original setup) ---
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([0.01, 0.1, 1])
ax.set_yticks([0.01, 0.1])
ax.minorticks_off()
ax.tick_params(labelsize=14)
plt.xlabel("CytoTRACE2 Risk Score", fontsize=18)
plt.ylabel("InferCNV Risk Score", fontsize=18)

# --- Colorbar (your original styling) ---
cbar = plt.colorbar(sc)
cbar.set_label("Weighted Agreement Score", fontsize=16)
cbar.ax.tick_params(labelsize=12)

plt.legend(loc='lower right', fontsize=10) 
plt.savefig(os.path.join(figures_dir, 'score_correlation_plot.pdf'), bbox_inches='tight')
plt.show()


#%% Calculate quartiles and 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle

# Calculate quartiles and mean/median
quartiles = np.percentile(agreement_score, [25, 50, 75])
q1, q2, q3 = quartiles[0], quartiles[1], quartiles[2]
mean = np.mean(agreement_score)
median = np.median(agreement_score)

# Ensure agreement_score matches the full dataset length
assert len(agreement_score) == len(adata_pp.obs), "Length mismatch!"

# Convert columns to NumPy arrays
cyto_all = adata_pp.obs['CytoTRACE2_Score'].values
cnv_all = adata_pp.obs['cnv_score'].values

plt.figure(figsize=(3, 6))
for i, (q_range, q_name) in enumerate([((0, q1), "Q1"), ((q1, q2), "Q2"), ((q2, q3), "Q3"), ((q3, 1), "Q4")]):
    mask = (agreement_score >= q_range[0]) & (agreement_score < q_range[1])
    cyto_q = cyto_all[mask]
    cnv_q = cnv_all[mask]
    
    if len(cyto_q) > 1:
        rho, _ = spearmanr(cyto_q, cnv_q)
        print(f"{q_name} Spearman ρ: {rho:.3f} (n={len(cyto_q)})")
        plt.bar(i, rho, label=f"{q_name} (n={len(cyto_q)})")
    else:
        print(f"{q_name}: Not enough data")

plt.xticks(range(4), ['Q1\n (25%)', 'Q2\n (50%)', 'Q3\n (75%)', 'Q4\n (100%)'], fontsize=13)
plt.ylabel('Spearman Correlation (ρ)', fontsize=18)
plt.tick_params(axis='y', labelsize=14)
plt.title('')
plt.legend('')
plt.savefig(os.path.join(figures_dir, 'quartile_correlation_analysis.pdf'), bbox_inches='tight')
plt.show()


#%% Define quartile ranges (from lowest to highest agreement)
quartile_ranges = [(0, q1), (q1, q2), (q2, q3), (q3, 1)]  # Q1, Q2, Q3, Q4
quartile_names = ['Q1', 'Q2', 'Q3', 'Q4']

selected_threshold = None  # Initialize

# Check quartiles from lowest (Q1) to highest (Q4)
for q_range, q_name in zip(quartile_ranges, quartile_names):
    mask = (agreement_score >= q_range[0]) & (agreement_score < q_range[1])
    cyto_q = cyto_all[mask]
    cnv_q = cnv_all[mask]
    
    if len(cyto_q) > 1:
        rho, _ = spearmanr(cyto_q, cnv_q)
        print(f"{q_name} Spearman ρ: {rho:.3f} (range: {q_range[0]:.3f}-{q_range[1]:.3f})")
        
        if rho > 0.5:
            selected_threshold = q_range[1]  # Use UPPER bound of this quartile
            print(f"\nSELECTED: {q_name} (ρ = {rho:.3f}) with threshold > {q_range[0]:.3f}")
            break
    else:
        print(f"{q_name}: Not enough data")

# Fallback if no quartile meets criteria (use q1)
if selected_threshold is None:
    selected_threshold = q1
    print(f"\nNo quartile had ρ > 0.5. Using median threshold: {q2:.3f}")
    


#%% Plot quartiles as histogram
# Define colors for quartiles
colors = ['#ff0000', '#ff9900', '#ffff00', '#00cc00']

# Create figure and set x-axis range to start at 0.3
plt.figure(figsize=(6.5, 6))
ax = plt.gca()

# Plot histogram with specified x-range
n, bins, patches = plt.hist(agreement_score, bins=50, edgecolor='white', alpha=0.8, range=(0.3, 1))

# Color bins by quartile
for i in range(len(patches)):
    bin_center = (bins[i] + bins[i+1]) / 2
    if bin_center <= q1:
        color = colors[0]  # Q1: Red
    elif bin_center <= q2:
        color = colors[1]  # Q2: Orange
    elif bin_center <= q3:
        color = colors[2]  # Q3: Yellow
    else:
        color = colors[3]  # Q4: Green
    patches[i].set_facecolor(color)

y_max = plt.ylim()[1]

# Determine which quartile boundary is just below our selected threshold
if selected_threshold <= q1:
    filter_boundary = 0.3  # If threshold is in Q1, filter everything below 0.3
elif selected_threshold <= q2:
    filter_boundary = q1   # If threshold is in Q2, filter Q1
elif selected_threshold <= q3:
    filter_boundary = q2   # If threshold is in Q3, filter Q2
else:
    filter_boundary = q3   # If threshold is in Q4, filter Q3

# Add shaded region for filtered cells
ax.add_patch(Rectangle((0.3, 0), width=filter_boundary-0.3, height=y_max, 
              facecolor='#ffcccc', alpha=0.4, zorder=0))  # Light red with higher alpha

# Add text in the filtered region
filtered_text_x = 0.3 + (q1 - 0.3)/2  # Center of the filtered region
plt.text(filtered_text_x, y_max*0.4, "Filtered cells\nwith low score agreement", 
         fontsize=17, ha='center', va='center', color='darkred',
         bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

# Add quartile lines with staggered annotations

annotation_heights = [y_max * 0.5, y_max * 0.75, y_max * 0.95]  # Staggered heights

for q, q_label, height in zip(quartiles, ['Q1 (25%)', 'Q2 (50%)', 'Q3 (75%)'], annotation_heights):
    plt.axvline(q, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    plt.text(q, height, f'{q:.2f}', fontsize=18,
             ha='center', va='top', color='black', fontweight='bold',
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

# Custom legend (upper left)
legend_elements = [
    Patch(facecolor=colors[0], label='Q1 (0-25%): Lowest Agreement'),
    Patch(facecolor=colors[1], label='Q2 (25-50%)'),
    Patch(facecolor=colors[2], label='Q3 (50-75%)'),
    Patch(facecolor=colors[3], label='Q4 (75-100%): Highest Agreement')
]
legend = plt.legend(handles=legend_elements, loc='upper left', fontsize=13)

# Add stats box under the legend (left-aligned)
stats_text = (
    f"Mean: {mean:.2f} ± {np.std(agreement_score):.2f}\n"
    f"Median: {median:.2f}"
)
plt.text(0.02, 0.65, stats_text, ha='left', va='top', fontsize=14,
         transform=ax.transAxes, 
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

# Labels and title
plt.xlabel('Agreement Score \n(0 = No agreement, 1 = Perfect agreement)', fontsize=18)
plt.ylabel('Number of cells', fontsize=18)

plt.grid(axis='y', alpha=0.3)
ax.tick_params(labelsize=12)
ax.set_xlim(0.3, 1.0)  # Explicitly set x-axis range

plt.tight_layout()
plt.savefig(os.path.join(figures_dir, 'agreement_score_distribution.pdf'), bbox_inches='tight')
plt.show()


#%% Filtered scatter plot using the new threshold
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr

# Filter cells with agreement > LOWER bound of the selected quartile
filtered_idx = agreement_score > q_range[0]  # Changed to use lower bound
print(f"\nFinal filtering: agreement_score > {q_range[0]:.3f}")
print(f"Cells retained: {sum(filtered_idx)}/{len(agreement_score)}")

filtered_array1 = array1[filtered_idx]
filtered_array2 = array2[filtered_idx]
filtered_agreement = agreement_score[filtered_idx]

# Compute correlations on filtered data
spearman_rho, spearman_p = spearmanr(filtered_array1, filtered_array2)
pearson_r, pearson_p = pearsonr(filtered_array1, filtered_array2)

plt.figure(figsize=(5, 6))

# --- Scatter plot with filtered data ---
sc = plt.scatter(
    filtered_array1, 
    filtered_array2, 
    c=filtered_agreement,
    cmap="RdYlGn",
    alpha=0.5, 
    s=5,
    vmin=agreement_score.min(),  # Keep original color scale
    vmax=agreement_score.max()   # for comparison
)

# --- Add regression line for filtered data ---
slope, intercept = np.polyfit(filtered_array1, filtered_array2, deg=1)
x_vals = np.linspace(filtered_array1.min(), filtered_array1.max(), 100)
y_fit = slope * x_vals + intercept

plt.plot(
    x_vals, 
    y_fit, 
    color='royalblue',
    linestyle='--',
    linewidth=1.5
)

# --- Add correlation stats ---
plt.text(
    0.05, 0.95, 
    f"Agreement\n score > {q1:.2f}\n\n"
    f"Spearman\n ρ = {spearman_rho:.2f} \n(p = {spearman_p:.2e})\n\n"
    f"Pearson\n r = {pearson_r:.2f} \n(p = {pearson_p:.2e})\n\n"
    f"\n\n\nCells retained: \n{len(filtered_array1)}/{len(array1)}",
    transform=plt.gca().transAxes,
    fontsize=18,
    ha='left', 
    va='top',
    bbox=dict(facecolor='white', alpha=0, edgecolor='gray')
)

# --- Format axes ---
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([0.01, 0.1, 1])
ax.set_yticks([0.01, 0.1])
ax.minorticks_off()
ax.tick_params(labelsize=12)
plt.xlabel("CytoTRACE2 Risk Score", fontsize=18)
plt.ylabel("InferCNV Risk Score", fontsize=18)

# --- Colorbar ---
cbar = plt.colorbar(sc)
cbar.set_label("Weighted Agreement Score", fontsize=18)
cbar.ax.tick_params(labelsize=14)

plt.legend(loc='lower right', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, 'filtered_score_correlation_plot.pdf'), bbox_inches='tight')
plt.show()

#%% Filter the Non-cancer cells with low agreement before making the new UMAP

cancer_agreement = agreement_score[agreement_score > q_range[0]]

intersection = list(set(cancer_agreement.index).intersection(set(adata_pp.obs[adata_pp.obs['CytoTRACE2_Score'] > infer_cnv_threshold].index)))

adata_pp.obs['Final_cancer_cell_status'] = np.where(adata_pp.obs.index.isin(intersection), 'Cancer cell', 'Normal cell')
adata_pp.obs['Final_cancer_cell_status']

#%% Final umap of all aggreement between models
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.text import Text
import matplotlib.patheffects as pe
from adjustText import adjust_text  # Make sure to install: pip install adjust-text

# Set larger default font sizes
rcParams['font.size'] = 14
rcParams['axes.titlesize'] = 16
rcParams['axes.labelsize'] = 14
rcParams['legend.fontsize'] = 12
rcParams['figure.titlesize'] = 16

# Create figure with larger subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    2, 2, 
    figsize=(14, 13),
    gridspec_kw={"wspace": 0.2, "hspace": 0.2}
)

# Plot 1: CytoTRACE2 Score
sc.pl.umap(
    adata_pp, 
    color="CytoTRACE2_Score", 
    ax=ax1, 
    title="CytoTRACE2 Cancer Risk Score", 
    show=False,
    frameon=False,
    color_map='plasma',
    size=5
)
ax1.title.set_size(22)  

# Plot 2: CNV Score
sc.pl.umap(
    adata_pp, 
    color="cnv_score", 
    ax=ax2, 
    title="InferCNV Cancer Risk Score", 
    show=False,
    frameon=False,
    color_map='plasma',
    size=5
)
ax2.title.set_size(22) 

# Plot 3: Cancer Cell Status with custom colors
categories = adata_pp.obs['Final_cancer_cell_status'].cat.categories
color_map = {
    'Cancer cell': '#FF7F0E',  # Orange
    'Normal cell': '#1F77B4'   # Blue
}

combined_effects = [
    pe.withStroke(linewidth=6, foreground="white"),  # Thick white border (outer)
    pe.withStroke(linewidth=1, foreground="black"),  # Thin black border (inner)
    pe.Normal()  # Original text color
]

# Plot UMAP
sc.pl.umap(
    adata_pp, 
    color="Final_cancer_cell_status", 
    ax=ax3, 
    title="Cancer Cell Status", 
    show=False,
    frameon=False,
    palette=[color_map[c] for c in categories],
    size=5,
    legend_loc=None
)
ax3.title.set_size(22)

# Manually add non-overlapping labels using adjustText
if 'X_umap' in adata_pp.obsm:
    umap_coords = adata_pp.obsm['X_umap']
    texts = []
    for category in categories:
        # Get center point for each category
        mask = adata_pp.obs['Final_cancer_cell_status'] == category
        center = umap_coords[mask].mean(axis=0)
        texts.append(ax3.text(
            center[0], center[1], 
            category, 
            fontsize=18,
            ha='center',
            va='center',
            color=color_map[category],
            path_effects=combined_effects 
        ))
    
    # Adjust text positions to avoid overlap
    adjust_text(
        texts,
        ax=ax3,
        expand_points=(1.2, 1.2),
        expand_text=(1.2, 1.2),
        force_text=0.5
    )

# Plot 4: Final Annotation
sc.pl.umap(
    adata_pp, 
    color="final_annotation", 
    ax=ax4, 
    show=False,
    frameon=False,
    size=5
)
ax4.set_title("Cell Type Annotations", fontsize=22)

# Adjust scale bars and legends
for ax in [ax1, ax2, ax3, ax4]:
    # Increase scale bar text size
    for child in ax.get_children():
        if isinstance(child, Text):
            if 'μm' in child.get_text():
                child.set_size(14)
    
    # Adjust legend if present
    if ax.get_legend() is not None:
        legend = ax.get_legend()
        legend.set_title(legend.get_title().get_text(), prop={'size': 14})
        for text in legend.get_texts():
            text.set_fontsize(12)

plt.tight_layout()
plt.show()


#%% After defining the Cancer cells identifywhat cells are present with a plotly pie graph
# This section will throw an error if run in a script due to opening a browser.
# I am commenting this out for this reason.
# To uncomment thi quickly highlight the section and hit Ctlr + 1
import sys
import plotly.io as pio
import plotly.express as px
import pandas as pd
import numpy as np
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
pio.renderers.default = 'svg'
# Get all possible cell types and assign colors FIRST
all_celltypes = adata_pp.obs['final_annotation'].astype('category').cat.remove_unused_categories()
cmap = plt.get_cmap('turbo')
values = np.linspace(0, 1, len(all_celltypes.cat.categories))
color_map = dict(zip(all_celltypes.cat.categories, [to_hex(cmap(value)) for value in values]))

# Filter for cancer cells and calculate percentages
cancer_cells = adata_pp.obs[adata_pp.obs['Final_cancer_cell_status'] == 'Cancer cell']
celltype_counts = cancer_cells['final_annotation'].value_counts()
percentages = (celltype_counts / celltype_counts.sum() * 100).round(2)

# Filter out 0% categories while maintaining color mapping
nonzero_mask = percentages > 0
filtered_labels = percentages.index[nonzero_mask].tolist()
filtered_percentages = percentages[nonzero_mask].tolist()
filtered_colors = [color_map[label] for label in filtered_labels]  # Get colors for only non-zero types

# Create the pie chart with only non-zero categories
fig = {
    'data': [{
        'values': filtered_percentages,
        'labels': filtered_labels,
        'marker': {'colors': filtered_colors},
        'type': 'pie',
        'rotation': 140,
        'textinfo': 'percent+label',
        'insidetextorientation': 'radial',
        'hovertemplate': '<b>%{label}</b><br>%{value:.1f}%<extra></extra>',
        'showlegend': False
    }],
    'layout': {
        'font': {'size': 18},
        'margin': {'t': 120, 'b': 20, 'l': 20, 'r': 20},
        'uniformtext': {
            'minsize': 12,
            'mode': 'hide'
        }
    }
}
pio.write_image(fig, os.path.join(figures_dir, 'cancer_cell_type_distribution.pdf'), width=1200, height=900)
pio.show(fig)


#%% This ends the cancer detection part. We are going to save here and start up the next part with NMF

# Save your changes to the adata_pp object
adata_pp.write("adata_pp_CD.h5ad")

#%% End file
import sys

sys.exit()
