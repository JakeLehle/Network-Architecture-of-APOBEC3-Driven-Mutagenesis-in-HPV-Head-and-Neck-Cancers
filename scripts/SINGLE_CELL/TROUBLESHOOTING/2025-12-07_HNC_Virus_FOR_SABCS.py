#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 13:13:46 2025

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
        ax = sc.pl.umap(adata_obj, color=group_name, show=False, legend_loc=None, title='', frameon=False, size=5, palette=palette)
        gen_mpl_labels(
            adata_obj,
            group_name,
            exclude=("None",),
            ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(
                fontsize=36, 
                path_effects=combined_effects  # Use the combined effects here
            ),
            color_by_group=True
        )
        fig = ax.get_figure()
        fig.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'Final_Cell_Type_Annotation_UMAP.pdf'))
        plt.show()
        plt.close()

#%% Prep for reading in all the kraken2 files
# read in the hierarchy.txt file with all he viruses used form the kraken2 database. This is output in each folder and is the same for all of them so I just need to read in one.
with open(os.path.join(working_dir, 'fastq', x_data['series_id'][0], x_data['run_accession'][0], x_data['run_accession'][0]+'_S1_L001_/outs/kraken2_filtered_feature_bc_matrix/hierarchy.txt'), 'rt') as f:
    lines = f.readlines()
import re
column_4 = []
column_5 = []
for line in lines:
    # Assuming your columns are separated by spaces
    #row = re.split(r'\t+', line.rstrip('\n'))
    pattern = re.compile('[^\t]+')
    pattern_2 = re.compile(r'.*?([a-zA-Z].*)')
    row = pattern.findall(line.rstrip('\n'))
    # Replace 0 with the column index you want to extract (starting from 0)
    column_4.append(row[4])
    column_5.append(pattern_2.findall(row[5])[0])
     
v_hierarchy = {key: value for key, value in zip(column_4, column_5)}
#%% Make sure that all files have all viruses if they don't then you wont be able to concatenate them easily. 
counter = 0
import gzip
import csv

# read the tsv
for key in x_data['series_id'].unique():
    for accession in x_data['run_accession'][x_data['series_id'] == key]:
        try:
            genes_tsv_path = os.path.join(working_dir, 'fastq', key, accession, 
                                        accession+"_S1_L001_/outs/kraken2_filtered_feature_bc_matrix/genes.tsv.gz")
            matrix_mtx_path = os.path.join(working_dir, 'fastq', key, accession, 
                                         accession+"_S1_L001_/outs/kraken2_filtered_feature_bc_matrix/matrix.mtx.gz")
            
            if not os.path.exists(genes_tsv_path):
                raise FileNotFoundError(f"genes.tsv.gz not found for {accession}")
            if not os.path.exists(matrix_mtx_path):
                raise FileNotFoundError(f"matrix.mtx.gz not found for {accession}")

            accession_list = []
            with gzip.open(genes_tsv_path, 'rt') as f:
                tsv_reader = csv.reader(f, delimiter="\t")
                for row in tsv_reader:
                    accession_list.append(row[1])
                accession_list_final = [x for x in accession_list if x.startswith(tuple(v_hierarchy.values()))]
                
                for tax_id, virus in v_hierarchy.items():
                    if virus in accession_list_final:
                        print(f'{virus} reads found in {accession}')
                    else:
                        with gzip.open(genes_tsv_path, 'at') as f:
                            tsv_writer = csv.writer(f, delimiter="\t", lineterminator="\n")
                            tsv_writer.writerow([tax_id, virus])       
            
            with gzip.open(genes_tsv_path, 'rt') as f:
                feature_tsv = []
                tsv_reader = csv.reader(f, delimiter="\t")
                for row in tsv_reader:
                    feature_tsv.append(row[1])
                feature_tsv_len = len(feature_tsv)
                print(f'Finished writing {accession} : gene.tsv.gz')
            
            with gzip.open(matrix_mtx_path) as file:
                lines = file.readlines()
                tmp = str(feature_tsv_len) + ' ' + str(lines[3 - 1]).split(' ', 2)[1] + ' ' + str(len(lines) - 3) + '\n'
                lines[3 - 1] = bytes(tmp, 'utf-8')
            
            with gzip.open(matrix_mtx_path, "wb") as file:
                for line in lines:
                    file.write(line)
                print(f'Finished writing {accession} : matrix.mtx.gz')
            
            counter = counter + 1
            print(f'Counter: {counter}')
            
        except FileNotFoundError as e:
            print(f"ERROR: Skipping {accession} - {str(e)}")
            continue
        except Exception as e:
            print(f"ERROR: Unexpected error processing {accession}: {str(e)}")
            continue 
#%% Make the function that will create the adata_v object
import gzip
import csv
from tqdm import tqdm
import anndata as ad
import os

def create_adata_viral_sc():
    adata_v = None
    
    for i in tqdm(range(len(x_data['series_id'])), desc='Reading viral anndata'):
        try:
            # Build paths and check files
            data_file = os.path.join(working_dir, 'fastq', x_data['series_id'][i], 
                                   x_data['run_accession'][i], 
                                   x_data['run_accession'][i]+'_S1_L001_/outs/kraken2_filtered_feature_bc_matrix')
            
            # Check if directory exists
            if not os.path.exists(data_file):
                print(f"ERROR: Directory not found - {data_file}")
                continue
                
            file_list = ['barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz']
            
            # Check all required files exist
            missing_files = [f for f in file_list if not os.path.exists(os.path.join(data_file, f))]
            if missing_files:
                print(f"ERROR: Missing files in {data_file}: {', '.join(missing_files)}")
                continue
            
            # Process files
            for file in file_list:
                try:
                    with gzip.open(os.path.join(data_file, file), 'rt') as f_in, \
                         open(os.path.join(data_file, file[:-3]), 'wt') as f_out:
                        f_out.writelines(f_in)
                except Exception as e:
                    print(f"ERROR: Failed to decompress {file} in {data_file}: {str(e)}")
                    # Clean up partially extracted files
                    if os.path.exists(os.path.join(data_file, file[:-3])):
                        os.remove(os.path.join(data_file, file[:-3]))
                    continue
            
            # Read the data
            try:
                adata_v_tmp = sc.read_10x_mtx(data_file)
                
                # Add metadata
                for metadata in x_data.columns:
                    adata_v_tmp.obs[metadata] = [x_data[metadata][i]]*adata_v_tmp.n_obs
                
                # Process observation names
                result = list(map("-".join, zip(adata_v_tmp.obs_names.to_list(), 
                                             adata_v_tmp.obs.run_accession.to_list())))
                adata_v_tmp.obs_names = result    
                adata_v_tmp.var['gene_symbol'] = adata_v_tmp.var.index
                
                # Concatenate with main object
                if adata_v is None:
                    adata_v = adata_v_tmp
                else:
                    adata_v = ad.concat([adata_v, adata_v_tmp], join='outer', merge='same')
                
            except Exception as e:
                print(f"ERROR: Failed to process 10x data in {data_file}: {str(e)}")
                continue
            
            # Clean up extracted files
            for file in file_list:
                extracted_file = os.path.join(data_file, file[:-3])
                if os.path.exists(extracted_file):
                    os.remove(extracted_file)
                    
        except Exception as e:
            print(f"ERROR: Unexpected error processing sample {x_data['run_accession'][i]}: {str(e)}")
            continue
    
    if adata_v is None:
        print("ERROR: No valid data was processed - returning empty AnnData object")
        return ad.AnnData()
    
    return adata_v
#%% Create the adata_v object
adata_v = create_adata_viral_sc()

#%% Save your work
for target in list(adata_v.obs.columns):
    adata_v.obs[target] = [str(element) for element in adata_v.obs[target]]

adata_v.write("adata_v.h5ad")

#%% load back in with all datasets
import scanpy as sc
adata = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata.h5ad')
)

adata_v = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_v.h5ad')
)

adata_tmp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_tmp.h5ad')
)

adata_pp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_pp_CD.h5ad')
)
#%% Find all human viruses using the human kraken database
with open(os.path.join(os.environ['HOME'], 'WORKING/kraken2/human_viral/inspect.txt'), 'rt') as f:
    lines = f.readlines()

import re
v_hierarchy_human = {}
v_hierarchy_human_species_list = []

for line in lines:
    # Skip empty lines
    if not line.strip():
        continue
        
    # Split on whitespace (but keep the last columns together)
    parts = line.strip().split()
    
    # We need at least 6 columns to process
    if len(parts) < 6:
        continue
        
    # Extract columns (adjust indices based on actual structure)
    # The pattern appears to be: percentage, read_count, tax_count, rank_code, tax_id, name
    percentage = parts[0]
    read_count = parts[1]
    tax_count = parts[2]
    rank_code = parts[3]  # This is the S/D/K/etc. column
    tax_id = parts[4]
    # The name is everything after tax_id (columns 5+)
    name = ' '.join(parts[5:])
    
    # Store in dictionary with tax_id as key
    v_hierarchy_human[tax_id] = {
        'rank_code': rank_code,
        'name': name,
        'read_count': read_count,
        'tax_count': tax_count,
        'percentage': percentage
    }
    
    # Add to species list if it's a species (starts with S)
    if rank_code.startswith('S'):
        v_hierarchy_human_species_list.append({
            'tax_id': tax_id,
            'name': name,
            'read_count': read_count
        })

# Create the name lists you originally wanted
v_hierarchy_human_list = tuple(entry['name'] for entry in v_hierarchy_human.values())
v_hierarchy_human_species_names = tuple(entry['name'] for entry in v_hierarchy_human_species_list)

# I'll also make a list of all the gene names from adata_pp which will come in handy later
adata_pp_gene_names = list(adata_pp.var.index)

#%% Strip all viruses not found in humans away from the adata_v object
human_virus_names = tuple(name.strip() for name in v_hierarchy_human_species_names if name.strip())
#%%
# Create a mask for human viruses in adata_v
human_virus_mask = adata_v.var_names.isin(human_virus_names)
# Drop all the non-human viruses
adata_v = adata_v[:, human_virus_mask].copy()

#%% Normalize the viral reads
# There will be some cells that now have 0 reads and need to be dropped before normalizing
sc.pp.filter_cells(adata_v, min_counts=1)
# The viral data has no mitochonria data so I can't use the QC_on_data_normal() function on it
# Just gonna normalize it
sc.pp.normalize_total(adata_v, target_sum=1e6)
sc.pp.log1p(adata_v)

#%% filter adata_v by cell barcodes in adata_tmp
import anndata as ad
adata_v_filtered = adata_v[adata_v.obs.index.isin(list(adata_pp.obs_names))]
adata_filtered = adata[adata.obs.index.isin(list(adata_pp.obs_names))]

adata_filtered.varm = {}
adata_v_filtered.varm = {}
adata_v_filtered.var
adata_filtered.var

adata_joined = ad.concat(
    [adata_filtered, adata_v_filtered],
    join='outer',
    axis=1
)


adata_joined = adata_joined[adata_pp.obs_names, :].copy()


assert all(adata_joined.obs_names == adata_pp.obs_names), (
    f"Order mismatch!\n"
    f"Joined order: {adata_joined.obs_names[:3]}\n"
    f"Original order: {adata_pp.obs_names[:3]}"
)

adata_joined.obs = adata_pp.obs.copy()
adata_joined.var
assert list(adata_pp.obs.index) == list(adata_joined.obs.index), "Cell barcodes don't match"

#%% Make a post processing copy This is if you want to include virus
adata_pp = adata_joined.copy()
#%% Calculate the neighborhood plot and cluster
from os import cpu_count

NCPUS = cpu_count()
sc.pp.neighbors(adata_pp, n_pcs=NCPUS)
sc.tl.umap(adata_pp)
sc.tl.leiden(adata_pp, key_added='clusters', resolution=1, random_state=42)

def sort_by_substrings(list_of_strings, substring):
    with_substring = []
    without_substring = []
    for string in list_of_strings:
        if substring in string:
            with_substring.append(string)
        else:
            without_substring.append(string)

    return with_substring + without_substring

sorted_sources = sort_by_substrings(list(adata_pp.obs.source_name.unique()), 'norm')

sorted_sources
adata_pp.obs['source_name_sorted'] = pd.Categorical(
    values=adata_pp.obs.source_name, categories=sorted_sources, ordered=True
)
cmap = plt.get_cmap('turbo')

# Generate evenly spaced values equal to the value of groups in the dataset
values = np.linspace(0, 1, len(adata_pp.obs.run_accession.unique()))

# Get RGB values for each value in the colormap
color_list = [cmap(value) for value in values]

ncol = 3
nrow = 1
figsize = 4
wspace = 2
# Adapt figure size based on number of rows and columns and added space between them
# (e.g. wspace between columns)
fig, axs = plt.subplots(
    nrow, ncol, figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
)
plt.subplots_adjust(wspace=wspace)
sc.pl.umap(adata_pp, color=['source_name_sorted'], cmap=cmap, size=5, ax=axs[0], show=False)
sc.pl.umap(adata_pp, color=['clusters'], cmap=cmap, size=5, ax=axs[1], show=False)
sc.pl.umap(adata_pp, color=['run_accession'], palette=color_list, size=5, ax=axs[2], save='UMAP_before_batch_correction')

sc.external.pp.bbknn(adata_pp, batch_key='run_accession')  # running bbknn 1.3.6
sc.tl.umap(adata_pp)

fig, axs = plt.subplots(
    nrow, ncol, figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
)
plt.subplots_adjust(wspace=wspace)
sc.pl.umap(adata_pp, color=['source_name_sorted'], cmap=cmap, ax=axs[0], show=False)
sc.pl.umap(adata_pp, color=['clusters'], cmap=cmap, ax=axs[1],  show=False)
sc.pl.umap(adata_pp, color=['run_accession'], palette=color_list, ax=axs[2], save='UMAP_after_batch_correction')

#%% Rank the genes
sc.tl.rank_genes_groups(adata_pp, groupby='final_annotation', key_added='rank_genes', method='wilcoxon')

for target in list(adata_pp.obs.columns):
    adata_pp.obs[target] = [str(element) for element in adata_pp.obs[target]]

#del adata_pp.uns['markers']

adata_pp.var['gene_ids'] = [str(item) for item in adata_pp.var['gene_ids']]

adata_pp.write("adata_v_pp.h5ad")

#%%
adata_pp_v = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_v_pp.h5ad')
)
#%%
adata_pp
plt.figure(figsize=(4, 4))
sc.pl.umap(adata_pp, color=['Final_cancer_cell_status'], size=5, frameon=False)
sc.pl.umap(adata_pp, color=['final_annotation'], size=5 )
sc.pl.umap(adata_pp, color=['SBS2'], size=15, cmap = 'plasma', frameon=False, show=False)
ax = plt.gca()
ax.title.set_fontsize(26)
plt.show()
sc.pl.umap(adata_pp, color=['SBS1'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS2'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS3'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS4'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS5'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS7a'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS7b'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS7d'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS13'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS16'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS17a'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS17b'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS18'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS33'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['SBS40a'], size=15, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['APOBEC3B'], size=10, cmap = 'plasma', frameon=False, show=False)
ax = plt.gca()
ax.title.set_fontsize(36)
plt.show()
sc.pl.umap(adata_pp, color=['APOBEC3A'], size=10, cmap = 'plasma', frameon=False, show=False)
ax = plt.gca()
ax.title.set_fontsize(36)
plt.show()
#%%
sc.pl.umap(adata_pp, color=['Human papillomavirus 16'], size=5, cmap = 'plasma', frameon=False, show=False)
ax = plt.gca()
ax.title.set_fontsize(26)
plt.show()
adata_pp
#%%11-22-2025 THIS SECTION WAS FOR THE RFA FOR THE HVC 
# Task 1: Transfer HPV16 read counts from adata_pp_v to adata_pp
# Extract HPV16 counts from adata_pp_v
hpv16_gene = 'Human papillomavirus 16'
hpv16_counts = adata_pp_v[:, hpv16_gene].X.toarray().flatten()

# Create a Series with cell barcodes as index
hpv16_series = pd.Series(hpv16_counts, index=adata_pp_v.obs.index, name='HPV16_counts')

# Map to adata_pp, filling missing cells with NaN then converting to 0
adata_pp.obs['Human papillomavirus 16'] = adata_pp.obs.index.map(hpv16_series).fillna(0)
adata_pp.obs['Human papillomavirus 16'][adata_pp.obs['Human papillomavirus 16'] > 0].describe()
print(f"Added HPV16_counts to adata_pp.obs. Non-zero cells: {(adata_pp.obs['Human papillomavirus 16'] > 0).sum()}")
# Task 2: Add SBS COSMIC signature weights
import pandas as pd

# Load the signature weights file
sbs_weights = pd.read_csv('/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/signature_refitting_hnscc_all_cells/signature_weights_per_cell.txt', 
                          sep='\t', index_col=0)

# Each row is a signature, each column is a cell barcode
# Transpose so columns become rows (easier to iterate)
sbs_weights_t = sbs_weights.T

# Add each signature as a new column in adata_pp.obs
for signature in sbs_weights.index:
    # Map signature weights to cells, filling missing with 0
    adata_pp.obs[signature] = adata_pp.obs.index.map(sbs_weights.loc[signature]).fillna(0)
adata_pp.obs['SBS2'][adata_pp.obs['SBS2'] > 0].describe()    
print(f"Added {len(sbs_weights.index)} SBS signatures to adata_pp.obs")
print(f"Signature names: {list(sbs_weights.index)}")

#%% Recalculate the UMAP plot for the normal and cancerous basal cells
import anndata as ad
adata_pp
# Subset the data
basal_cell = adata_pp[adata_pp.obs['final_annotation'] == "basal cell"]
cancer_cells = basal_cell[
    (basal_cell.obs['source_name'] == 'head and neck squamous cell carcinoma') &
    (basal_cell.obs['Human papillomavirus 16'] > 0) &
    (basal_cell.obs['normalized_total_mutations'] > 1.899019e-05)
]
normal_cells = basal_cell[
    (basal_cell.obs['source_name'] == 'normal tissue adjucent to head and neck squamous cell carcinoma') &
    (basal_cell.obs['Human papillomavirus 16'] == 0) &
    (basal_cell.obs['normalized_total_mutations'] > 0)
]
cancer_cells.obs['SBS2'].describe()
cancer_cells.obs['Human papillomavirus 16'].describe()
custom_percentile = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
cancer_cells.obs['normalized_total_mutations'].describe(percentiles=custom_percentile)
normal_cells.obs['SBS2'].describe()

#%%
import matplotlib.pyplot as plt
import scanpy as sc

# Set font sizes (adjust these as needed)
AXIS_LABEL_SIZE = 18
TICK_LABEL_SIZE = 8

# Extract data
x = basal_cell.obs['SBS2']
y = basal_cell.X[:, basal_cell.var_names.get_loc('APOBEC3B')]

# If using sparse matrix, convert to dense
if hasattr(y, 'toarray'):
    y = y.toarray().flatten()

# # Filter out cells with zero values in either dimension
# mask = (x > 0) & (y > 0)
# x_filtered = x[mask]
# y_filtered = y[mask]

# # Create figure
# fig, ax = plt.subplots(figsize=(4.5, 5))

# # Create scatter plot
# ax.scatter(x_filtered, y_filtered, 
#           color='#3273a6', 
#           alpha=0.7,
#           s=10,  # adjust point size as needed
#           edgecolors='none')

# Create figure
fig, ax = plt.subplots(figsize=(4.5, 5))

# Create scatter plot
ax.scatter(x, y, 
          color='#3273a6', 
          alpha=0.7,
          s=10,  # adjust point size as needed
          edgecolors='none')

# Set labels with custom font sizes
ax.set_xlabel('SBS2 Weights', fontsize=AXIS_LABEL_SIZE, labelpad=20)
ax.set_ylabel('Log Normlaized A3B Expression \nfrom Basal Cells', fontsize=AXIS_LABEL_SIZE, labelpad=20)

# Set tick label sizes
ax.tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE)
ax.grid(False)
# Remove top and right spines for cleaner look (optional)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# Set remaining spines to black
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Set tick colors to black as well
ax.tick_params(axis='both', which='both', colors='black')

plt.tight_layout()
plt.savefig('apobec3b_sbs2_scatter.pdf', dpi=300, bbox_inches='tight')
plt.show()
#%%
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# Create custom colormap: light gray for zero, vibrant red for max
colors = ['#D3D3D3', '#FF0000']  # Light gray to bright red
n_bins = 100
cmap = LinearSegmentedColormap.from_list('gray_to_red', colors, N=n_bins)

# Create figure with 3 subplots side by side
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Plot 1: HPV16
sc.pl.umap(adata_pp, color=['Human papillomavirus 16'], 
           size=5, frameon=False, show=False, ax=axes[0], cmap=cmap)
axes[0].set_title('Human papillomavirus 16', fontsize=26, fontweight='bold')

# Plot 2: Normalized total mutations
sc.pl.umap(adata_pp, color=['normalized_total_mutations'], 
           size=20, frameon=False, show=False, ax=axes[1], cmap='plasma')
axes[1].set_title('Normalized Total Mutations', fontsize=26, fontweight='bold')

# Plot 3: SBS2
sc.pl.umap(adata_pp, color=['SBS2'], 
           size=20, frameon=False, show=False, ax=axes[2], cmap='plasma')
axes[2].set_title('SBS2', fontsize=26, fontweight='bold')

# Adjust layout
plt.tight_layout()
plt.show()
#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats

# Match sample sizes - take top 438 cancer cells with highest mutation counts
# Assuming you have a column with normalized mutation counts in cancer_cells.obs
# Adjust the column name if needed (e.g., 'n_counts', 'total_counts', or your mutation count column)
cancer_cells_matched = cancer_cells.obs.nlargest(438, 'normalized_total_mutations')  # Change 'n_counts' to your mutation count column name

# Create matched cancer AnnData subset
cancer_cells_subset = cancer_cells[cancer_cells_matched.index].copy()

# Combine the data for plotting
plot_data = pd.DataFrame({
    'SBS2': pd.concat([normal_cells.obs['SBS2'], cancer_cells_subset.obs['SBS2']]),
    'Cell_Type': ['Normal HPV16-'] * len(normal_cells) + ['Cancer HPV16+'] * len(cancer_cells_subset)
})

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Prepare data for violin plot
normal_sbs2 = normal_cells.obs['SBS2'].values
cancer_sbs2 = cancer_cells_subset.obs['SBS2'].values
data_for_violin = [normal_sbs2, cancer_sbs2]

# Create violin plot
parts = ax.violinplot(data_for_violin, 
                      positions=[1, 2],
                      widths=0.7,
                      showmeans=False,
                      showmedians=False,
                      showextrema=False)

# Color the violins
colors = ['#3498DB', '#E74C3C']  # Blue for normal, red for cancer
for pc, color in zip(parts['bodies'], colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.6)
    pc.set_edgecolor('black')
    pc.set_linewidth(2)

# Add jittered points
np.random.seed(42)  # For reproducibility
for i, (data, color) in enumerate(zip(data_for_violin, colors), 1):
    # Create jitter
    jitter = np.random.normal(0, 0.04, size=len(data))
    x = np.ones(len(data)) * i + jitter
    
    ax.scatter(x, data, alpha=0.3, s=30, color=color, edgecolors='black', linewidths=0.5)

# Calculate and add medians as horizontal lines
normal_median = np.median(normal_sbs2)
cancer_median = np.median(cancer_sbs2)

ax.hlines(normal_median, 0.7, 1.3, colors='black', linewidth=4, label='Median', zorder=10)
ax.hlines(cancer_median, 1.7, 2.3, colors='black', linewidth=4, zorder=10)

# Add median values as text
ax.text(1.35, normal_median, f'{normal_median:.3f}', 
        ha='left', va='center', fontsize=22, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', linewidth=2))
ax.text(2.35, cancer_median, f'{cancer_median:.3f}', 
        ha='left', va='center', fontsize=22, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', linewidth=2))

# Set x-axis labels
ax.set_xticks([1, 2])
ax.set_xticklabels(['Normal HPV16-\nBasal Cells\n(n = 438)', 
                    'Cancer HPV16+\nBasal Cells\n(n = 438)'])

# Labels and title
ax.set_ylabel('SBS2 Weight', fontsize=26, fontweight='bold')
ax.set_xlabel('Cell Population', fontsize=26, fontweight='bold')
ax.set_title('SBS2 Signature Weights:\nNormal vs Cancer Basal Cells (Matched n)', 
             fontsize=28, fontweight='bold', pad=20)

# Tick parameters
ax.tick_params(axis='both', which='major', labelsize=26, width=2, length=8)
ax.tick_params(axis='x', rotation=0)

# Grid for readability
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)
ax.set_axisbelow(True)

# Perform statistical tests
# 1. Shapiro-Wilk test for normality (for n=438, this is appropriate)
_, p_normal_shapiro = stats.shapiro(normal_sbs2)
_, p_cancer_shapiro = stats.shapiro(cancer_sbs2)

print(f"\nNormality Tests (Shapiro-Wilk):")
print(f"  Normal cells: p = {p_normal_shapiro:.4e}")
print(f"  Cancer cells: p = {p_cancer_shapiro:.4e}")

# 2. Choose appropriate test based on normality
if p_normal_shapiro > 0.05 and p_cancer_shapiro > 0.05:
    # Both normal - use t-test
    statistic, p_value = stats.ttest_ind(normal_sbs2, cancer_sbs2)
    test_used = "Two-sample t-test"
else:
    # Non-normal - use Mann-Whitney U test
    statistic, p_value = stats.mannwhitneyu(normal_sbs2, cancer_sbs2, alternative='two-sided')
    test_used = "Mann-Whitney U test"

# 3. Calculate effect size (Cohen's d)
pooled_std = np.sqrt((np.std(normal_sbs2, ddof=1)**2 + np.std(cancer_sbs2, ddof=1)**2) / 2)
cohens_d = (np.mean(cancer_sbs2) - np.mean(normal_sbs2)) / pooled_std

print(f"\n{test_used}:")
print(f"  Test statistic: {statistic:.2f}")
print(f"  P-value: {p_value:.2e}")
print(f"  Cohen's d (effect size): {cohens_d:.3f}")

# Add p-value to the plot
if p_value < 0.001:
    p_text = 'p < 0.001***'
elif p_value < 0.01:
    p_text = f'p = {p_value:.3f}**'
elif p_value < 0.05:
    p_text = f'p = {p_value:.3f}*'
else:
    p_text = f'p = {p_value:.3f} (n.s.)'

# Add significance bracket
y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
y_pos = ax.get_ylim()[1] - 0.05 * y_range
ax.plot([1, 1, 2, 2], [y_pos, y_pos + 0.02*y_range, y_pos + 0.02*y_range, y_pos], 
        'k-', linewidth=3)
ax.text(1.5, y_pos + 0.03*y_range, p_text, ha='center', fontsize=24, fontweight='bold')

# Adjust layout
plt.tight_layout()
plt.show()

# Print detailed statistics
print(f"\nNormal HPV16- Basal Cells (n=438):")
print(f"  Median SBS2: {normal_median:.4f}")
print(f"  Mean SBS2: {np.mean(normal_sbs2):.4f}")
print(f"  Std SBS2: {np.std(normal_sbs2, ddof=1):.4f}")
print(f"  IQR: {np.percentile(normal_sbs2, 75) - np.percentile(normal_sbs2, 25):.4f}")

print(f"\nCancer HPV16+ Basal Cells (n=438, top mutation burden):")
print(f"  Median SBS2: {cancer_median:.4f}")
print(f"  Mean SBS2: {np.mean(cancer_sbs2):.4f}")
print(f"  Std SBS2: {np.std(cancer_sbs2, ddof=1):.4f}")
print(f"  IQR: {np.percentile(cancer_sbs2, 75) - np.percentile(cancer_sbs2, 25):.4f}")
#%%
def plot_signature_cosmic_style(component_data, ax, title, bg_color='white'):
    """
    Plot mutation signature in COSMIC style on a given axis
    """
    # Normalize to percentages
    total = component_data.sum()
    percentages = (component_data / total * 100) if total > 0 else component_data
    
    # Set background color
    ax.set_facecolor(bg_color)
    
    # Extract contexts and mutation types
    contexts = component_data.index.tolist()
    
    # Parse mutation contexts correctly
    # Format should be like: A[C>T]G
    mutation_types = []
    left_bases = []
    right_bases = []
    
    for ctx in contexts:
        # Extract mutation type from brackets
        mut = ctx.split('[')[1].split(']')[0]  # Gets "C>T"
        mutation_types.append(mut)
        
        # Extract flanking bases
        left_base = ctx.split('[')[0]  # Gets "A"
        right_base = ctx.split(']')[1]  # Gets "G"
        left_bases.append(left_base)
        right_bases.append(right_base)
    
    # Create position array
    x = np.arange(len(contexts))
    bar_width = 0.8
    
    # Plot mutation spectrum
    bar_colors = [COSMIC_COLORS[mut] for mut in mutation_types]
    bars = ax.bar(x, percentages, width=bar_width, 
                 color=bar_colors, edgecolor='black', linewidth=0.5)
    
    # Add separation lines between mutation types
    for i in range(1, 6):
        ax.axvline(i*16 - 0.5, color='gray', linestyle='--', alpha=0.7)
    
    # Set y-axis limits
    y_max = max(percentages.max() * 1.2, 1)
    ax.set_ylim(0, y_max)
    
    # Add colored rectangles above each mutation type section
    mutation_categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    for i, mut in enumerate(mutation_categories):
        rect_x = i * 16
        rect = plt.Rectangle((rect_x - 0.5, y_max), 16, 0.03 * y_max,
                           facecolor=COSMIC_COLORS[mut], clip_on=False)
        ax.add_patch(rect)
        
        # Add mutation type label above rectangle
        ax.text(rect_x + 7.5, y_max * 1.05, mut,
               ha='center', va='bottom', fontsize=16, fontweight='bold')
    
    # Add trinucleotide context labels
    for i, (left, mut, right) in enumerate(zip(left_bases, mutation_types, right_bases)):
        center_base = mut[0]  # Reference base from mutation type (e.g., 'C' from 'C>T')
        
        y_offset = 0.03 * y_max
        y_pos = -0.03 * y_max
        
        # Left base (black) - top position
        ax.text(i, y_pos, left, 
               rotation=90, ha='center', va='top', 
               fontsize=7, color='black', fontfamily='monospace')
        
        # Center base (colored) - middle position
        ax.text(i, y_pos - y_offset, center_base, 
               rotation=90, ha='center', va='center', 
               fontsize=7, color=COSMIC_COLORS[mut], fontfamily='monospace')
        
        # Right base (black) - bottom position
        ax.text(i, y_pos - 2*y_offset, right, 
               rotation=90, ha='center', va='bottom', 
               fontsize=7, color='black', fontfamily='monospace')
    
    # Set title
    ax.set_title(title, fontsize=18, pad=45, fontweight='bold')
    ax.set_xlim(-0.5, len(contexts) - 0.5)
    ax.set_ylabel('Percentage', fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Remove x-axis ticks/labels
    ax.set_xticks([])
    ax.set_xticklabels([])
    
    # Add y-axis tick parameters
    ax.tick_params(axis='y', labelsize=12)
#%%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad

# Combine cancer and normal cells
combined_adata = ad.concat([cancer_cells_subset, normal_cells], join='outer', merge='same')

# Add a grouping column to identify cell type
combined_adata.obs['HPV_status'] = ['Cancer_HPV16+'] * len(cancer_cells_subset) + ['Normal_HPV16-'] * len(normal_cells)

# Calculate PCA, neighbors, UMAP (optional but recommended for QC)
from os import cpu_count
NCPUS = cpu_count()
sc.pp.pca(combined_adata, n_comps=min(NCPUS, 50))  # Cap at 50 components
sc.pp.neighbors(combined_adata, n_pcs=min(NCPUS, 50))
sc.tl.umap(combined_adata)

# Optional: visualize the combined data
sc.pl.umap(combined_adata, color=['HPV_status'], frameon=False)

# Perform differential expression analysis
sc.tl.rank_genes_groups(combined_adata, groupby='HPV_status', 
                        key_added='rank_genes', 
                        method='wilcoxon',
                        use_raw=False)

# Extract results for Cancer vs Normal
deg_results = sc.get.rank_genes_groups_df(combined_adata, 'Cancer_HPV16+', key='rank_genes')

# Add negative log10 p-value for volcano plot
deg_results['neg_log10_pvalue'] = -np.log10(deg_results['pvals_adj'])

# Replace any infinite values (from p-value = 0) with a high value
deg_results['neg_log10_pvalue'] = deg_results['neg_log10_pvalue'].replace([np.inf], 
                                                                           deg_results['neg_log10_pvalue'][deg_results['neg_log10_pvalue'] != np.inf].max() + 10)

# Print summary statistics
print(f"\nDifferential Expression Summary:")
print(f"Total genes tested: {len(deg_results)}")
print(f"Significant genes (FDR < 0.05): {(deg_results['pvals_adj'] < 0.05).sum()}")
print(f"Upregulated (logFC > 1.5, FDR < 0.05): {((deg_results['logfoldchanges'] > 1.5) & (deg_results['pvals_adj'] < 0.05)).sum()}")
print(f"Downregulated (logFC < -1.5, FDR < 0.05): {((deg_results['logfoldchanges'] < -1.5) & (deg_results['pvals_adj'] < 0.05)).sum()}")

# Display top upregulated genes
print("\nTop 10 upregulated genes:")
top_up = deg_results[(deg_results['logfoldchanges'] > 1.5) & (deg_results['pvals_adj'] < 0.05)].nlargest(10, 'logfoldchanges')
print(top_up[['names', 'logfoldchanges', 'pvals_adj']])

# Display top downregulated genes
print("\nTop 10 downregulated genes:")
top_down = deg_results[(deg_results['logfoldchanges'] < -1.5) & (deg_results['pvals_adj'] < 0.05)].nsmallest(10, 'logfoldchanges')
print(top_down[['names', 'logfoldchanges', 'pvals_adj']])

# Create volcano plot using your enhanced function
enhanced_volcano_plot(deg_results)
#%%
import gseapy as gp

# Get top upregulated genes for GSEA
top_up_genes = deg_results[(deg_results['logfoldchanges'] > 1.5) & 
                           (deg_results['pvals_adj'] < 0.05)].nlargest(100, 'neg_log10_pvalue')

# Run enrichment analysis
enr = gp.enrichr(gene_list=list(top_up_genes['names']),
                 gene_sets=['KEGG_2021_Human'],
                 organism='human',
                 outdir=None,
                )

# Plot results
from gseapy import dotplot
ax = dotplot(enr.results,
             column="Adjusted P-value",
             x='Gene_set',
             size=10,
             top_term=15,
             figsize=(8, 10),
             title="GSEA: Cancer HPV16+ vs Normal HPV16- Basal Cells",
             xticklabels_rot=45,
             show_ring=True,
             marker='o',
            )
plt.tight_layout()
plt.show()

# Display top enriched pathways
print("\nTop 10 enriched pathways:")
print(enr.results[['Term', 'Adjusted P-value', 'Overlap', 'Genes']].head(10))
#%%
combined_adata = ad.concat([cancer_cells, normal_cells], join='outer', merge='same')

# Calculate the neighborhood plot and cluster
from os import cpu_count

NCPUS = cpu_count()
sc.pp.pca(combined_adata, n_comps=NCPUS)
sc.pp.neighbors(combined_adata, n_pcs=NCPUS)
sc.tl.umap(combined_adata)
sc.tl.leiden(combined_adata, key_added='cancer_clusters', resolution=1, random_state=42)

#%% Visulization of the updated UMAP plot

sc.pp.normalize_total(combined_adata, target_sum=1e6)
sc.pp.log1p(combined_adata)

plt.figure(figsize=(4, 4))
sc.pl.umap(combined_adata, color=['Final_cancer_cell_status'], size=5, frameon=False)
sc.pl.umap(combined_adata, color=['cancer_clusters'], size=5, )
sc.pl.umap(combined_adata, color=['APOBEC3B'], size=10, cmap = 'plasma', frameon=False, show=False)
ax = plt.gca()
ax.title.set_fontsize(36)
plt.show()
sc.pl.umap(combined_adata, color=['APOBEC3A'], size=10, cmap = 'plasma', frameon=False, show=False)
ax = plt.gca()
ax.title.set_fontsize(36)
plt.show()
sc.pl.umap(combined_adata, color=['Human papillomavirus 16'], size=10, cmap = 'plasma', frameon=False, show=False)
ax = plt.gca()
ax.title.set_fontsize(26)
plt.show()
adata_pp_v = adata_pp.copy()
#%% Refine the UMAP for Cancer vs Normal basal cells

status_palette = {
    'Cancer cell': '#fda600',  # Orange
    'Normal cell': '#4781b5',  # Blue

}

from matplotlib.patches import Patch
with plt.rc_context({"figure.figsize": (7.5, 6), "figure.dpi": 150, "figure.frameon": False}):
    fig = sc.pl.umap(combined_adata, color=['Final_cancer_cell_status'], palette=status_palette, frameon=False, size=25, legend_loc=None, return_fig=True, title='')
    
    # Manually create legend with squares
    legend_elements = [
        Patch(facecolor='#fda600', edgecolor='k', label='Cancer cell'),
        Patch(facecolor='#4781b5', edgecolor='k', label='Normal cell')
    ]
    
    ax = fig.axes[0]
    
    # Place legend to the right of the plot
    legend = ax.legend(
        handles=legend_elements,
        fontsize=46,
        frameon=False,
        bbox_to_anchor=(0.6, 1),  # Coordinates for legend position (1.05 = right outside)
        loc='upper left',           # Anchor point for the legend
        borderaxespad=0            # Padding between legend and axes
    )
    
    # Adjust the figure to make room for the legend
    plt.tight_layout(rect=[0, 0, 1.2, 1])  # rect=[left, bottom, right, top] (0.85 leaves 15% space on right)
    plt.savefig(os.path.join(figures_dir, 'Cancer_Status_Basal_Cells.pdf'), bbox_inches='tight')
    plt.show()


#%% Patient specific analysis 
# Create the mapping dictionary
sample_to_patient = {
    'GSM5268284': 'SC003',
    'GSM5268298': 'SC003',
    'GSM5268285': 'SC005',
    'GSM5268299': 'SC005',
    'GSM5268286': 'SC006',
    'GSM5268300': 'SC006'
}

# Map the experiment_alias to patient IDs, filling with NaN for unmapped samples
combined_adata.obs['patient_mapping'] = combined_adata.obs['experiment_alias'].map(sample_to_patient)
#%% Plot UMAP for Patient Mapping with table
patient_palette = {
    'SC003': '#008f5a',  # Green
    'SC005': '#6d398b',  # Purple
    'SC006': '#f28e1c',  # Orange
    'NaN': '#f5f5f5',  # Gray
}

with plt.rc_context({"figure.figsize": (10, 8), "figure.dpi": 150, "figure.frameon": False}):
    fig = sc.pl.umap(combined_adata, color=['patient_mapping'], palette=patient_palette, frameon=False, size=25, legend_loc=None, return_fig=True, title='')
    
    # Manually create legend with squares
    legend_elements = [
        Patch(facecolor='#008f5a', edgecolor='k', label='SC003'),
        Patch(facecolor='#6d398b', edgecolor='k', label='SC005'),
        Patch(facecolor='#f28e1c', edgecolor='k', label='SC006'),
        Patch(facecolor='#f5f5f5', edgecolor='k', label='Non-matched')
    ]
    
    ax = fig.axes[0]
    
    # Place legend to the right of the plot
    legend = ax.legend(
        handles=legend_elements,
        fontsize=40,
        frameon=False,
        bbox_to_anchor=(1, 1),  # Coordinates for legend position (1.05 = right outside)
        loc='upper left',           # Anchor point for the legend
        borderaxespad=0            # Padding between legend and axes
    )
    
    # Adjust the figure to make room for the legend
    plt.tight_layout(rect=[0, 0, 1.2, 1])  # rect=[left, bottom, right, top] (0.85 leaves 15% space on right)
    plt.savefig(os.path.join(figures_dir, 'Patient_mapping_UMAP.pdf'), bbox_inches='tight')
    plt.show()

import pandas as pd

# Get counts (using crosstab for margins)
count_table = pd.crosstab(
    combined_adata.obs['patient_mapping'],
    combined_adata.obs['Final_cancer_cell_status'],
    margins=True,
    margins_name="Total"
)

print(count_table)  # Verify before plotting

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Create a figure with a table
fig, ax = plt.subplots(figsize=(8, 4))
ax.axis('off')  # Hide axes

# Plot table with styling
table = plt.table(
    cellText=count_table.values,
    colLabels=count_table.columns,
    rowLabels=count_table.index,
    loc='center',
    cellLoc='center',
    colColours=['#f0f0f0'] * len(count_table.columns),  # Header color
    rowColours=['#f0f0f0'] * len(count_table.index)     # Row label color
)

# Style adjustments
table.auto_set_font_size(False)
table.set_fontsize(20)
table.scale(1, 5)  # Expand cells

# Title
plt.title("Cell Counts by Patient and Cancer Status",pad=45, fontsize=20)
# Save to PDF
with PdfPages('cell_counts_by_patient.pdf') as pdf:
    pdf.savefig(fig, bbox_inches='tight')

plt.show()


#%%
adata_SC006 = combined_adata[combined_adata.obs['patient_mapping'] == 'SC006']

#%% Calculate the DEGs for these new subsetted basal cell groups
sc.tl.rank_genes_groups(adata_SC006, groupby='source_name', key_added='rank_genes', method='wilcoxon')
SC006_cell_cancer = sc.get.rank_genes_groups_df(adata_SC006, 'head and neck squamous cell carcinoma', key='rank_genes')
SC006_cell_cancer['neg_log10_pvalue'] = -np.log10(SC006_cell_cancer['pvals_adj'])

#%% Display the A3 specifically for any DEGs using a custom volcano plot
def A3_volcano(rg_df):
    A3 = rg_df[rg_df['names'].isin(['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H'])]
    plt.figure(figsize=(4.5, 4))
    plt.scatter(x=A3['logfoldchanges'],y=A3['neg_log10_pvalue'], s=35, color="black")
    down = A3[(A3['logfoldchanges']<=0)&(A3['pvals_adj']<=0.05)]
    up = A3[(A3['logfoldchanges']>=0)&(A3['pvals_adj']<=0.05)]
    plt.scatter(x=down['logfoldchanges'],y=down['neg_log10_pvalue'], s=35, color="blue")
    plt.scatter(x=up['logfoldchanges'],y=up['neg_log10_pvalue'], s=35, color="red")

    texts=[]
    for i,r in A3.iterrows():
        texts.append(plt.text(x=r['logfoldchanges'],y=r['neg_log10_pvalue'],s=r['names'], fontsize=16))
    adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=0.75))
    plt.axhline(1.3,color="grey",linestyle="--")
    plt.axvline(0,color="grey",linestyle="--")
    plt.xlabel("log"+'\u2082'+'FC', fontsize=18)
    plt.ylabel("-logFDR", fontsize=18)
    plt.show()
    plt.close()

#%%
A3_volcano(SC006_cell_cancer)

#%% Define a volcano plot that will indicate the top 100 genes based on significance but only plot them if they meet the logFC cutoff
def enhanced_volcano_plot(rg_df):
    # Create figure with larger size
    plt.figure(figsize=(18, 14))  # Increased figure size
    
    # Add all genes (grey)
    plt.scatter(x=rg_df['logfoldchanges'], 
                y=rg_df['neg_log10_pvalue'], 
                s=15, color="grey", alpha=0.5)
    
    # Add significant downregulated genes (blue)
    down = rg_df[(rg_df['logfoldchanges'] < -1.5) & (rg_df['pvals_adj'] < 0.05)]
    plt.scatter(x=down['logfoldchanges'], 
                y=down['neg_log10_pvalue'], 
                s=35, color="blue", alpha=0.7, label='Downregulated (logFC < -1.5)')
    
    # Add significant upregulated genes (red)
    up = rg_df[(rg_df['logfoldchanges'] > 1.5) & (rg_df['pvals_adj'] < 0.05)]
    plt.scatter(x=up['logfoldchanges'], 
                y=up['neg_log10_pvalue'], 
                s=35, color="red", alpha=0.7, label='Upregulated (logFC > 1.5)')
    
    # Add APOBEC3 family members (black outline)
    A3 = rg_df[rg_df['names'].isin(['APOBEC3B', 'APOBEC3C'])]
    plt.scatter(x=A3['logfoldchanges'], 
                y=A3['neg_log10_pvalue'], 
                s=150, facecolors='none', edgecolors='black', linewidths=2)
    
    # Get top genes with |logFC| > 1.5
    significant_genes = rg_df[(abs(rg_df['logfoldchanges']) > 1.5) & 
                             (rg_df['pvals_adj'] < 0.05)]
    top_25_genes = significant_genes.nlargest(25, 'neg_log10_pvalue')
    
    # Set axis limits with more buffer space for labels
    x_buffer = 0.2 * (rg_df['logfoldchanges'].max() - rg_df['logfoldchanges'].min())
    y_buffer = 0.2 * (rg_df['neg_log10_pvalue'].max() - rg_df['neg_log10_pvalue'].min())
    
    plt.xlim(rg_df['logfoldchanges'].min() - x_buffer, 
             rg_df['logfoldchanges'].max() + x_buffer)
    plt.ylim(0, rg_df['neg_log10_pvalue'].max() + y_buffer)
    
    xlim = plt.xlim()
    ylim = plt.ylim()
    
    # Parameters for label placement
    label_fontsize = 36  # Slightly smaller font
    line_alpha = 0.6
    line_width = 0.8
    
    # First place APOBEC3 labels with priority
    for _, row in A3.iterrows():
        # Calculate offset direction (away from center)
        x_offset = 0.075 if row['logfoldchanges'] > 0 else -0.075
        y_offset = 0.01
        
        # Initial position
        x_pos = row['logfoldchanges'] + x_offset * (xlim[1] - xlim[0])
        y_pos = row['neg_log10_pvalue'] + y_offset * (ylim[1] - ylim[0])
        
        # Create the label
        text = plt.text(x_pos, y_pos, row['names'], 
                       fontsize=label_fontsize,
                       fontweight='medium',
                       color='black')
        
        # Add connecting line
        plt.plot([row['logfoldchanges'], x_pos],
                 [row['neg_log10_pvalue'], y_pos],
                 color='gray', linestyle='-', linewidth=line_width, alpha=line_alpha)
    
    # Then place top DEG labels
    placed_labels = []
    for _, row in top_25_genes.iterrows():
        # Skip APOBEC genes since we already placed them
        if row['names'] in A3['names'].values:
            continue
            
        # Determine if gene is upregulated or downregulated
        is_upregulated = row['logfoldchanges'] > 0
        
        # Try multiple positions in an expanding grid pattern
        max_attempts = 40000
        placed = False
        
        for attempt in range(max_attempts):
            # Calculate grid positions
            grid_size = int(np.ceil(np.sqrt(max_attempts)))
            row_idx = attempt // grid_size
            col_idx = attempt % grid_size
            
            # Calculate position with directional bias
            if is_upregulated:
                # For upregulated genes, use right half (x > 0)
                x_frac = 0.5 + 0.4 * (col_idx / (grid_size-1))  # 0.5 to 0.9
                y_frac = 0.1 + 0.8 * (row_idx / (grid_size-1))  # 0.1 to 0.9
            else:
                # For downregulated genes, use left half (x < 0)
                x_frac = 0.1 + 0.4 * (col_idx / (grid_size-1))  # 0.1 to 0.5
                y_frac = 0.1 + 0.8 * (row_idx / (grid_size-1))  # 0.1 to 0.9
            
            # Calculate absolute coordinates
            x_pos = xlim[0] + x_frac * (xlim[1] - xlim[0])
            y_pos = ylim[0] + y_frac * (ylim[1] - ylim[0])
            
            # Check if this position is within plot bounds
            if (xlim[0] < x_pos < xlim[1] and ylim[0] < y_pos < ylim[1]):
                # Check for overlap with existing labels
                overlap = False
                for (existing_x, existing_y, existing_name) in placed_labels:
                    if ((abs(x_pos - existing_x) < 0.175 * (xlim[1] - xlim[0])) and 
                        (abs(y_pos - existing_y) < 0.09 * (ylim[1] - ylim[0]))):
                        overlap = True
                        break
                
                if not overlap:
                    # Place the label
                    text = plt.text(x_pos, y_pos, row['names'], 
                                   fontsize=label_fontsize,
                                   fontweight='medium',
                                   color='black')
                    
                    # Add connecting line
                    plt.plot([row['logfoldchanges'], x_pos],
                             [row['neg_log10_pvalue'], y_pos],
                             color='gray', linestyle='-', linewidth=line_width, alpha=line_alpha)
                    
                    # Record the placed label
                    placed_labels.append((x_pos, y_pos, row['names']))
                    placed = True
                    break
        
        if not placed:
            # Fallback - place at edge with arrow
            if is_upregulated:
                x_pos = xlim[1] - 0.05 * (xlim[1] - xlim[0])
            else:
                x_pos = xlim[0] + 0.05 * (xlim[1] - xlim[0])
            
            y_pos = row['neg_log10_pvalue']
            
            plt.annotate(row['names'], 
                         xy=(row['logfoldchanges'], row['neg_log10_pvalue']),
                         xytext=(x_pos, y_pos),
                         arrowprops=dict(arrowstyle="-", color='gray', 
                                        linestyle='-', linewidth=line_width, alpha=line_alpha),
                         fontsize=label_fontsize,
                         fontweight='medium',
                         color='black')
            placed_labels.append((x_pos, y_pos, row['names']))
    
    # Add significance threshold line
    plt.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=1.5)
    
    # Add logFC threshold lines
    plt.axvline(-1.5, color="blue", linestyle=":", linewidth=2, alpha=0.7)
    plt.axvline(1.5, color="red", linestyle=":", linewidth=2, alpha=0.7)
    plt.axvline(0, color="grey", linestyle="--", linewidth=1.5)
    
    # Labels and title
    plt.xlabel("log$_2$FC", fontsize=50)
    plt.ylabel("-log$_{10}$(FDR)", fontsize=50)
    plt.title("", fontsize=18, pad=20)
    
    # Legend
    plt.legend(fontsize=38, framealpha=1)
    
    # Adjust tick label sizes
    plt.xticks(fontsize=32)
    plt.yticks(fontsize=32)
    
    # Ensure all elements are visible
    plt.tight_layout()
    plt.show()

#%%
enhanced_volcano_plot(SC006_cell_cancer)

#%% Compared basal DEGs to those from network analysis.
top_SC006_cell_cancer = SC006_cell_cancer[
    (SC006_cell_cancer.pvals_adj < 0.05) & 
    (abs(SC006_cell_cancer.logfoldchanges) > 1.5)
]

network_genes = pd.read_csv(working_dir +'/top_genes.txt', sep='\t')
network_genes.index = network_genes['NETWORK_GENES']
network_genes
all_genes = combined_adata.var.loc[:, ['gene_ids', 'gene_symbol']]

all_genes_names = all_genes[all_genes['gene_ids'].isin(network_genes['NETWORK_GENES'])]

intersection = list(set(top_SC006_cell_cancer.names).intersection(set(list(all_genes_names['gene_symbol']))))
#%% Pick out genes with greater than a 0.8 cosine similarity to start 
top_gene_cor = combined_adata.uns['A3B_corr'][combined_adata.uns['A3B_corr'].Cosine_Similarity > 0.85]

list(top_gene_cor.Gene)
#%%


intersection = list(set(top_basal_cell_cancer.names).intersection(set(list(top_gene_cor.Gene))))


#%%


import gseapy as gp
names = gp.get_library_name(organism="Human")
names
# Run GSEA on KEGG pathways



# if you are only intrested in dataframe that enrichr returned, please set outdir=None
enr = gp.enrichr(gene_list=list(top_SC006_cell_cancer.names), # or "./tests/data/gene_list.txt",
                 #gene_sets=['WikiPathways_2024_Human', 'GO_Biological_Process_2025', 'GO_Molecular_Function_2025'],
                 gene_sets=['KEGG_2021_Human'] #'GO_Biological_Process_2025', 'GO_Molecular_Function_2025'],
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )

# simple plotting function
from gseapy import barplot, dotplot
# categorical scatterplot
ax = dotplot(enr.results,
              column="Adjusted P-value",
              x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
              size=10,
              top_term=15,
              figsize=(5,15),
              title = "Gene Set Enrichment Analysis",
              xticklabels_rot=45, # rotate xtick labels
              show_ring=True, # set to False to revmove outer ring
              marker='o',
             )

tmp = enr.results
    
#%%

#nonoverlapping_UMAP(adata_pp, 'final_annotation')

# #%% Strip all normal genes found away from the dataset to find the highest expression of viruses
# human_virus_names = tuple(name.strip() for name in v_hierarchy_human_species_names if name.strip())
# # Create a mask for the original genes in the adata_pp dataset
# human_virus_mask_updated = ~adata_pp.var_names.isin(adata_pp_gene_names)
# # Drop all the non-human viruses
# adata_v = adata_v[:, human_virus_mask].copy()


#%% Save this section it might be important later I will move it down to where it is needed 
# v_hierarchy_df = pd.DataFrame.from_dict(v_hierarchy, orient='index')

# virus_names = list(v_hierarchy_df[0])

# intersection = list(set(adata_pp.obs_names).intersection(set(adata_v.obs_names)))



# adata_joined = ad.concat([adata_pp, adata_v_filtered], join='outer', axis = 1)

# adata_joined


#%%
from tqdm import tqdm
results = adata_pp.uns['rank_genes']
len(results['names'])
remove_list = adata_pp_gene_names
out = np.array([[0,0,0,0,0]])
for group in tqdm(results['names'].dtype.names):
    out = np.vstack((out, np.vstack((results['names'][group],
                                     results['scores'][group],
                                     results['pvals_adj'][group],
                                     results['logfoldchanges'][group],
                                     np.array([group] * len(results['names'][group])).astype('object'))).T))
    markers = pd.DataFrame(out[1:], columns = ['virus', 'scores', 'pval_adj', 'lfc', 'cluster'])
    #This is the DEG object filtered for pval < 0.05 for ALL clusters
    markers = markers[(markers.scores > 0.05)]
    #markers = markers[(markers.pval_adj < 0.05)]
    mask = markers.virus.str.contains('|'.join(remove_list)) 
    # Create a mask to identify rows to remove
    markers = markers[~mask]

    markers['pval_adj'] = markers['pval_adj'].map('{:e}'.format)
    adata_pp.uns['markers'] = markers #save marker df to uns

#%% Filter out 

virus_markers = markers[markers['virus'].isin(human_virus_names)]

subset = virus_markers.iloc[:, [0,1]]

subset = subset.groupby('virus')['scores'].sum()
# List out the top 10 viruses found in the data
list(subset.sort_values(ascending=False).index)[:10]
subset.sort_values(ascending=False)

#%%
#This is just to set up how the figures will display
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, CSS4_COLORS
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3

sc.set_figure_params(scanpy=True, fontsize=25)

#Define a nice colour map for gene expression
colors_gray = plt.cm.Greys(np.linspace(0.4, 0.41, 30)) 
colors_red = plt.cm.Reds(np.linspace(0.3, 1, 226))    
colors_combined = np.vstack([colors_gray, colors_red])
mymap = LinearSegmentedColormap.from_list('my_improved_map', colors_combined)

sc.pl.umap(adata_pp, color=list(subset.sort_values(ascending=False).index)[:1], use_raw=False, size=5, color_map=mymap, frameon=False)
sc.pl.umap(adata_pp, color=['APOBEC3A'], use_raw=False, size=5, color_map=mymap, frameon=False)
sc.pl.matrixplot(adata_pp, list(['Human papillomavirus 16',
 'Molluscum contagiosum virus subtype 1',
 'Human papillomavirus 154',
 'Human betaherpesvirus 7',
 'NY_014 poxvirus',
 'Bourbon virus',
 'Human mastadenovirus B',
 'Songling virus',
 'Staphylococcus phage 6ec']), 'final_annotation', dendrogram=True, var_group_rotation=30, cmap = 'plasma', log=True)

plt.rcdefaults()


sc.pl.heatmap(
    adata_pp,
    list(subset.sort_values(ascending=False).index)[:],
    'final_annotation',
    cmap='plasma',
    log=True,
    show=False
)
ax = plt.gca()
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
plt.show()
#%%


del adata_pp.uns['markers']
adata_pp.var['gene_ids'] = [str(item) for item in adata_pp.var['gene_ids']]
adata_pp.var
adata_pp.write("adata_pp_HNC_V_HUMAN_VIRUS.h5ad")

adata_pp = sc.read_h5ad(
    filename="adata_pp_HNC_V.h5ad"
)


adata_pp

sc.pl.umap(adata_pp, color=['geo_accession'], palette=color_list)