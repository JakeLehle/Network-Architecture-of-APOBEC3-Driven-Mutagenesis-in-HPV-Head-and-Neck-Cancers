#%%
#!/usr/bin/env python3
"""
Calculate cosine, dot product, and inverse Eucledian distance similarity between single-cell samples and COSMIC SBS2/SBS13 signatures
Tailored for SigProfiler format mutation matrices
"""

import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def calculate_similarities(sample_matrix, cosmic_signatures, target_sigs=['SBS2', 'SBS13'], 
                         metrics=['cosine', 'dot_product', 'inverse_euclidean']):
    """
    Calculate multiple similarity metrics between each sample and target signatures
    
    Parameters:
    -----------
    sample_matrix : pd.DataFrame
        Mutation counts matrix (mutations x samples)
    cosmic_signatures : pd.DataFrame
        COSMIC signature matrix (mutations x signatures)
    target_sigs : list
        List of signature names to analyze
    metrics : list
        List of similarity metrics to calculate
    
    Returns:
    --------
    dict of pd.DataFrames with similarity scores for each metric
    """
    
    # Ensure mutation types are aligned
    common_mutations = sample_matrix.index.intersection(cosmic_signatures.index)
    print(f"\nFound {len(common_mutations)} common mutation types (out of 96 possible)")
    
    if len(common_mutations) < 96:
        print(f"Warning: Expected 96 mutation types but found {len(common_mutations)}")
    
    # Align both matrices to common mutations
    sample_matrix_aligned = sample_matrix.loc[common_mutations]
    cosmic_aligned = cosmic_signatures.loc[common_mutations]
    
    # Store results for each metric
    all_results = {}
    
    for metric in metrics:
        metric_results = {}
        
        for sig in target_sigs:
            if sig not in cosmic_aligned.columns:
                print(f"Error: {sig} not found in COSMIC signatures!")
                available_sbs = [col for col in cosmic_aligned.columns if col.startswith('SBS')]
                print(f"Available SBS signatures: {', '.join(available_sbs[:20])}...")
                continue
            
            # Get the signature vector (already normalized in COSMIC)
            sig_vector = cosmic_aligned[sig].values
            
            # Calculate similarity for each sample
            similarities = []
            
            for sample_id in sample_matrix_aligned.columns:
                sample_counts = sample_matrix_aligned[sample_id].values
                
                if metric == 'cosine':
                    # Cosine similarity (angle only, normalized vectors)
                    if sample_counts.sum() > 0:
                        sample_norm = sample_counts / sample_counts.sum()
                    else:
                        sample_norm = np.zeros_like(sample_counts)
                    
                    sim = cosine_similarity([sample_norm], [sig_vector])[0][0]
                
                elif metric == 'dot_product':
                    # Dot product (considers both angle and magnitude)
                    # Scale signature by total mutations in sample to make comparable
                    total_mutations = sample_counts.sum()
                    if total_mutations > 0:
                        # Scale signature to match sample's mutation burden
                        scaled_sig = sig_vector * total_mutations
                        sim = np.dot(sample_counts, scaled_sig)
                        # Normalize by the maximum possible dot product
                        max_possible = np.dot(scaled_sig, scaled_sig)
                        if max_possible > 0:
                            sim = sim / max_possible
                        else:
                            sim = 0
                    else:
                        sim = 0
                
                elif metric == 'inverse_euclidean':
                    # Inverse Euclidean distance
                    # Compare normalized counts to signature
                    if sample_counts.sum() > 0:
                        sample_norm = sample_counts / sample_counts.sum()
                    else:
                        sample_norm = np.zeros_like(sample_counts)
                    
                    # Calculate Euclidean distance
                    distance = np.sqrt(np.sum((sample_norm - sig_vector) ** 2))
                    
                    # Convert to similarity (inverse, bounded between 0 and 1)
                    # Maximum possible distance is sqrt(2) for normalized vectors
                    max_distance = np.sqrt(2)
                    sim = 1 - (distance / max_distance)
                
                similarities.append(sim)
            
            metric_results[sig] = pd.Series(similarities, index=sample_matrix_aligned.columns)
        
        all_results[metric] = pd.DataFrame(metric_results)
    
    return all_results

def create_visualizations(similarity_results, output_prefix='sbs_analysis'):
    """Create comprehensive visualizations for multiple similarity metrics"""
    
    # Set style
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")
    
    metrics = list(similarity_results.keys())
    
    # 1. Distribution plots for each metric
    fig, axes = plt.subplots(len(metrics), 2, figsize=(14, 6 * len(metrics)))
    if len(metrics) == 1:
        axes = axes.reshape(1, -1)
    
    for metric_idx, metric in enumerate(metrics):
        similarities_df = similarity_results[metric]
        
        for sig_idx, sig in enumerate(['SBS2', 'SBS13']):
            if sig in similarities_df.columns:
                ax = axes[metric_idx, sig_idx]
                
                # Histogram
                n, bins, patches = ax.hist(similarities_df[sig], bins=30, alpha=0.7, edgecolor='black')
                
                # Color bars by value
                for i, patch in enumerate(patches):
                    if bins[i] > 0.7:
                        patch.set_facecolor('darkgreen')
                    elif bins[i] > 0.5:
                        patch.set_facecolor('orange')
                    else:
                        patch.set_facecolor('lightblue')
                
                # Add statistics
                mean_val = similarities_df[sig].mean()
                median_val = similarities_df[sig].median()
                ax.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.3f}')
                ax.axvline(median_val, color='green', linestyle='--', label=f'Median: {median_val:.3f}')
                ax.axvline(0.7, color='black', linestyle=':', label='0.7 threshold')
                
                ax.set_title(f'{sig} {metric.replace("_", " ").title()} Distribution')
                ax.set_xlabel(f'{metric.replace("_", " ").title()} Similarity')
                ax.set_ylabel('Number of Cells')
                ax.legend()
                ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_distributions_all_metrics.png', dpi=300, bbox_inches='tight')
    print(f"Saved distribution plots as '{output_prefix}_distributions_all_metrics.png'")
    plt.close()
    
    # 2. Comparison scatter plots between metrics
    if len(metrics) > 1:
        fig, axes = plt.subplots(1, len(metrics)-1, figsize=(8*(len(metrics)-1), 8))
        if len(metrics) == 2:
            axes = [axes]
        
        # Compare each metric to cosine (or first metric)
        base_metric = 'cosine' if 'cosine' in metrics else metrics[0]
        other_metrics = [m for m in metrics if m != base_metric]
        
        for idx, metric in enumerate(other_metrics):
            ax = axes[idx]
            
            # Merge data for both signatures
            base_sbs2 = similarity_results[base_metric]['SBS2'] if 'SBS2' in similarity_results[base_metric].columns else pd.Series()
            metric_sbs2 = similarity_results[metric]['SBS2'] if 'SBS2' in similarity_results[metric].columns else pd.Series()
            
            # Plot SBS2
            if not base_sbs2.empty and not metric_sbs2.empty:
                ax.scatter(base_sbs2, metric_sbs2, alpha=0.6, label='SBS2', color='blue')
            
            # Plot SBS13
            if 'SBS13' in similarity_results[base_metric].columns and 'SBS13' in similarity_results[metric].columns:
                base_sbs13 = similarity_results[base_metric]['SBS13']
                metric_sbs13 = similarity_results[metric]['SBS13']
                ax.scatter(base_sbs13, metric_sbs13, alpha=0.6, label='SBS13', color='orange')
            
            # Add diagonal line
            lims = [0, 1]
            ax.plot(lims, lims, 'k--', alpha=0.3)
            
            ax.set_xlabel(f'{base_metric.replace("_", " ").title()} Similarity')
            ax.set_ylabel(f'{metric.replace("_", " ").title()} Similarity')
            ax.set_title(f'{base_metric.replace("_", " ").title()} vs {metric.replace("_", " ").title()}')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, 1.05)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_metric_comparisons.png', dpi=300, bbox_inches='tight')
        print(f"Saved metric comparison plots as '{output_prefix}_metric_comparisons.png'")
        plt.close()
    
    # 3. Heatmap of top cells across all metrics
    # Get top cells based on maximum similarity across all metrics and signatures
    all_max_sims = []
    for metric in metrics:
        df = similarity_results[metric]
        if 'SBS2' in df.columns and 'SBS13' in df.columns:
            max_sim = df[['SBS2', 'SBS13']].max(axis=1)
        elif 'SBS2' in df.columns:
            max_sim = df['SBS2']
        elif 'SBS13' in df.columns:
            max_sim = df['SBS13']
        else:
            continue
        all_max_sims.append(max_sim)
    
    if all_max_sims:
        # Get top 30 cells by average maximum similarity across metrics
        avg_max_sim = pd.concat(all_max_sims, axis=1).mean(axis=1)
        top_cells = avg_max_sim.nlargest(30).index
        
        # Create heatmap data
        heatmap_data = []
        row_labels = []
        
        for metric in metrics:
            df = similarity_results[metric]
            for sig in ['SBS2', 'SBS13']:
                if sig in df.columns:
                    heatmap_data.append(df.loc[top_cells, sig].values)
                    row_labels.append(f'{metric.replace("_", " ").title()} - {sig}')
        
        if heatmap_data:
            heatmap_df = pd.DataFrame(heatmap_data, index=row_labels, columns=top_cells)
            
            plt.figure(figsize=(16, 8))
            sns.heatmap(heatmap_df, cmap='YlOrRd', cbar_kws={'label': 'Similarity'}, 
                       fmt='.2f', linewidths=0.5)
            plt.title('Top 30 Cells - Similarity Across All Metrics', fontsize=16)
            plt.xlabel('Cell Barcode', fontsize=12)
            plt.ylabel('Metric - Signature', fontsize=12)
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.savefig(f'{output_prefix}_top_cells_heatmap.png', dpi=300, bbox_inches='tight')
            print(f"Saved top cells heatmap as '{output_prefix}_top_cells_heatmap.png'")
            plt.close()

def main():
    # File paths - update these to your actual paths
    matrix_file = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt"
    cosmic_file = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/COSMIC_v3.4_SBS_GRCh38.txt"
    output_prefix = "sbs_similarity_analysis"
    
    print("=== SBS2/SBS13 Multi-Metric Signature Analysis ===\n")
    
    # Load data
    print("Loading mutational matrix...")
    try:
        sample_matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)
        sample_matrix = sample_matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
        print(f"✓ Loaded matrix: {sample_matrix.shape[0]} mutations × {sample_matrix.shape[1]} cells")
    except Exception as e:
        print(f"✗ Error loading matrix: {e}")
        sys.exit(1)
    
    print("\nLoading COSMIC signatures...")
    try:
        cosmic_signatures = pd.read_csv(cosmic_file, sep='\t', index_col='Type')
        cosmic_signatures = cosmic_signatures.apply(pd.to_numeric, errors='coerce').fillna(0)
        print(f"✓ Loaded COSMIC: {cosmic_signatures.shape[0]} mutations × {cosmic_signatures.shape[1]} signatures")
    except Exception as e:
        print(f"✗ Error loading COSMIC: {e}")
        sys.exit(1)
    
    # Calculate similarities using multiple metrics
    print("\nCalculating similarities using multiple metrics...")
    print("  - Cosine similarity (angle only)")
    print("  - Dot product similarity (angle + magnitude)")
    print("  - Inverse Euclidean distance (overall difference)")
    
    similarity_results = calculate_similarities(
        sample_matrix, 
        cosmic_signatures, 
        ['SBS2', 'SBS13'],
        ['cosine', 'dot_product', 'inverse_euclidean']
    )
    
    # Process results for each metric
    for metric, similarities in similarity_results.items():
        if similarities.empty:
            print(f"✗ Error: No similarities calculated for {metric}")
            continue
        
        print(f"\n=== {metric.replace('_', ' ').title()} Results ===")
        
        # Add summary columns
        available_sigs = [sig for sig in ['SBS2', 'SBS13'] if sig in similarities.columns]
        if available_sigs:
            similarities['Max_Similarity'] = similarities[available_sigs].max(axis=1)
            similarities['Dominant_Signature'] = similarities[available_sigs].idxmax(axis=1)
            similarities['Similarity_Difference'] = abs(similarities['SBS2'] - similarities['SBS13'])
            
            # Sort by maximum similarity
            similarities_sorted = similarities.sort_values('Max_Similarity', ascending=False)
            
            # Save results
            output_file = f"{output_prefix}_{metric}_results.csv"
            similarities_sorted.to_csv(output_file, float_format='%.6f')
            print(f"✓ Results saved to: {output_file}")
            
            # Display summary statistics
            for sig in available_sigs:
                print(f"\n{sig}:")
                print(f"  Mean similarity: {similarities[sig].mean():.4f}")
                print(f"  Median similarity: {similarities[sig].median():.4f}")
                print(f"  Std deviation: {similarities[sig].std():.4f}")
                print(f"  Min: {similarities[sig].min():.4f}")
                print(f"  Max: {similarities[sig].max():.4f}")
                
                # Threshold analysis
                high_sim = (similarities[sig] > 0.7).sum()
                print(f"  Cells with similarity > 0.7: {high_sim} ({high_sim/len(similarities)*100:.1f}%)")
            
            # Store sorted results for visualization
            similarity_results[metric] = similarities_sorted
    
    # Create comprehensive output combining all metrics
    print("\n=== Creating Combined Analysis ===")
    
    # Create a summary DataFrame with best metric for each cell
    summary_data = []
    for cell in sample_matrix.columns:
        cell_data = {'Cell': cell}
        
        for metric in ['cosine', 'dot_product', 'inverse_euclidean']:
            if metric in similarity_results:
                df = similarity_results[metric]
                if cell in df.index:
                    for sig in ['SBS2', 'SBS13']:
                        if sig in df.columns:
                            cell_data[f'{metric}_{sig}'] = df.loc[cell, sig]
        
        summary_data.append(cell_data)
    
    summary_df = pd.DataFrame(summary_data).set_index('Cell')
    
    # Add consensus columns
    for sig in ['SBS2', 'SBS13']:
        metric_cols = [col for col in summary_df.columns if sig in col]
        if metric_cols:
            summary_df[f'{sig}_mean'] = summary_df[metric_cols].mean(axis=1)
            summary_df[f'{sig}_std'] = summary_df[metric_cols].std(axis=1)
    
    # Sort by average SBS2+SBS13 score
    if 'SBS2_mean' in summary_df.columns and 'SBS13_mean' in summary_df.columns:
        summary_df['Overall_mean'] = (summary_df['SBS2_mean'] + summary_df['SBS13_mean']) / 2
        summary_df = summary_df.sort_values('Overall_mean', ascending=False)
    
    summary_df.to_csv(f"{output_prefix}_all_metrics_summary.csv", float_format='%.6f')
    print(f"✓ Combined summary saved to: {output_prefix}_all_metrics_summary.csv")
    
    # Display top cells across all metrics
    print("\n=== Top 10 Cells by Average Score Across All Metrics ===")
    display_cols = [col for col in summary_df.columns if 'mean' in col or 'std' in col]
    print(summary_df[display_cols].head(10))
    
    # Identify consensus high-similarity cells
    consensus_threshold = 0.7
    consensus_cells = []
    
    for sig in ['SBS2', 'SBS13']:
        # Find cells where all metrics agree on high similarity
        high_sim_cells = set()
        first_metric = True
        
        for metric in ['cosine', 'dot_product', 'inverse_euclidean']:
            col_name = f'{metric}_{sig}'
            if col_name in summary_df.columns:
                metric_high = set(summary_df[summary_df[col_name] > consensus_threshold].index)
                if first_metric:
                    high_sim_cells = metric_high
                    first_metric = False
                else:
                    high_sim_cells = high_sim_cells.intersection(metric_high)
        
        if high_sim_cells:
            consensus_cells.extend([(cell, sig) for cell in high_sim_cells])
            print(f"\n{len(high_sim_cells)} cells show {sig} > {consensus_threshold} across ALL metrics")
            consensus_df = summary_df.loc[list(high_sim_cells)]
            consensus_df.to_csv(f"{output_prefix}_{sig}_consensus_high_similarity.csv", float_format='%.6f')
    
    # Create visualizations
    print("\nCreating visualizations...")
    create_visualizations(similarity_results, output_prefix)
    
    print("\n=== Analysis Complete ===")
    
    # Final insights
    print("\nKey Insights:")
    print("\n1. Metric Comparison:")
    print("   - Cosine: Measures pattern similarity (angle only)")
    print("   - Dot Product: Considers both pattern and mutation burden")
    print("   - Inverse Euclidean: Overall profile difference")
    
    if consensus_cells:
        print(f"\n2. Consensus high-similarity cells found: {len(set([c[0] for c in consensus_cells]))}")
        print("   These cells show consistent signatures across all metrics")
    else:
        print("\n2. No cells show consistent high similarity across all metrics")
        print("   This suggests weak or variable signature presence")
    
    print("\n3. Recommendations:")
    print("   - For sparse data, dot product may be most informative")
    print("   - Check cells with high dot product but low cosine similarity")
    print("   - These have the right mutation burden but different patterns")

if __name__ == "__main__":
    main()

#%%
#!/usr/bin/env python3
"""
Plot individual cell mutation profiles in COSMIC signature format
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def plot_cell_mutation_signature(mutation_matrix, cell_barcode, cosmic_signatures=None, 
                               reference_sig='SBS2', save_path=None):
    """
    Plot a single cell's mutation profile in COSMIC signature format
    
    Parameters:
    -----------
    mutation_matrix : pd.DataFrame
        Mutation counts (rows: mutation types, columns: cells)
    cell_barcode : str
        Cell barcode to plot
    cosmic_signatures : pd.DataFrame, optional
        COSMIC signatures for comparison
    reference_sig : str
        Which COSMIC signature to compare to
    save_path : str
        Path to save the figure
    """
    
    # Check if cell exists
    if cell_barcode not in mutation_matrix.columns:
        print(f"Error: Cell {cell_barcode} not found in mutation matrix")
        return None
    
    # Get cell data
    cell_data = mutation_matrix[cell_barcode]
    
    # Define mutation types and colors (COSMIC standard)
    mutation_types = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    colors = {
        'C>A': '#98D7EC',  # Light blue
        'C>G': '#212121',  # Black
        'C>T': '#FF003A',  # Red
        'T>A': '#A0A0A0',  # Gray
        'T>C': '#83A603',  # Green
        'T>G': '#F5ABCC'   # Pink
    }
    
    # Prepare figure
    n_plots = 2 if cosmic_signatures is not None else 1
    fig, axes = plt.subplots(n_plots, 1, figsize=(20, 6 * n_plots))
    if n_plots == 1:
        axes = [axes]
    
    # Plot cell profile
    ax = axes[0]
    
    # Organize data by mutation type
    x_labels = []
    y_values = []
    bar_colors = []
    
    for mut_type in mutation_types:
        # Get all contexts for this mutation type
        contexts = sorted([idx for idx in cell_data.index if f'[{mut_type}]' in idx])
        for context in contexts:
            x_labels.append(context)
            y_values.append(cell_data[context])
            bar_colors.append(colors[mut_type])
    
    # Normalize to frequencies
    total_mutations = sum(y_values)
    if total_mutations > 0:
        y_values_norm = [y / total_mutations for y in y_values]
    else:
        y_values_norm = y_values
    
    # Create bar plot
    x_pos = np.arange(len(x_labels))
    bars = ax.bar(x_pos, y_values_norm, color=bar_colors, width=1.0, edgecolor='none')
    
    # Formatting
    ax.set_xlim(-0.5, len(x_labels) - 0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([label.replace('[', '\n').replace(']', '') for label in x_labels], 
                       rotation=90, ha='center', fontsize=8)
    ax.set_ylabel('Relative Frequency', fontsize=12)
    ax.set_title(f'{cell_barcode} - Mutation Profile (n={total_mutations} mutations)', 
                fontsize=14, fontweight='bold')
    
    # Add mutation type labels
    current_x = 0
    for mut_type in mutation_types:
        n_contexts = sum(1 for label in x_labels if f'[{mut_type}]' in label)
        if n_contexts > 0:
            mid_x = current_x + n_contexts / 2 - 0.5
            ax.text(mid_x, ax.get_ylim()[1] * 0.95, mut_type, 
                   ha='center', va='top', fontsize=10, fontweight='bold')
            # Add background shading
            ax.add_patch(Rectangle((current_x - 0.5, 0), n_contexts, ax.get_ylim()[1], 
                                 facecolor=colors[mut_type], alpha=0.1, zorder=0))
            current_x += n_contexts
    
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_axisbelow(True)
    
    # Plot reference COSMIC signature if provided
    if cosmic_signatures is not None and reference_sig in cosmic_signatures.columns:
        ax2 = axes[1]
        
        cosmic_data = cosmic_signatures[reference_sig]
        
        # Use same order as cell plot
        cosmic_values = [cosmic_data[label] if label in cosmic_data.index else 0 
                        for label in x_labels]
        
        # Plot
        bars2 = ax2.bar(x_pos, cosmic_values, color=bar_colors, width=1.0, edgecolor='none')
        
        # Formatting
        ax2.set_xlim(-0.5, len(x_labels) - 0.5)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([label.replace('[', '\n').replace(']', '') for label in x_labels], 
                           rotation=90, ha='center', fontsize=8)
        ax2.set_ylabel('Probability', fontsize=12)
        ax2.set_title(f'COSMIC {reference_sig} Reference', fontsize=14, fontweight='bold')
        
        # Add mutation type labels
        current_x = 0
        for mut_type in mutation_types:
            n_contexts = sum(1 for label in x_labels if f'[{mut_type}]' in label)
            if n_contexts > 0:
                mid_x = current_x + n_contexts / 2 - 0.5
                ax2.text(mid_x, ax2.get_ylim()[1] * 0.95, mut_type, 
                       ha='center', va='top', fontsize=10, fontweight='bold')
                ax2.add_patch(Rectangle((current_x - 0.5, 0), n_contexts, ax2.get_ylim()[1], 
                                      facecolor=colors[mut_type], alpha=0.1, zorder=0))
                current_x += n_contexts
        
        ax2.grid(True, axis='y', alpha=0.3)
        ax2.set_axisbelow(True)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {save_path}")
    else:
        plt.show()
    
    return fig

def plot_multiple_cells(mutation_matrix, cell_list, cosmic_signatures=None, 
                       reference_sig='SBS2', output_prefix='cell_signatures'):
    """
    Plot mutation profiles for multiple cells
    """
    for i, cell in enumerate(cell_list):
        save_path = f"{output_prefix}_{i+1}_{cell}.png"
        plot_cell_mutation_signature(mutation_matrix, cell, cosmic_signatures, 
                                   reference_sig, save_path)

def extract_top_cells_from_results(results_file, signature='SBS2', n_top=5):
    """
    Extract top N cells from a similarity results file
    """
    df = pd.read_csv(results_file, index_col=0)
    if signature in df.columns:
        top_cells = df.nlargest(n_top, signature)
        print(f"\nTop {n_top} cells for {signature} from {results_file}:")
        for i, (cell, row) in enumerate(top_cells.iterrows()):
            print(f"  {i+1}. {cell}: {row[signature]:.4f}")
        return top_cells.index.tolist()
    else:
        print(f"Error: {signature} not found in {results_file}")
        return []

# Example usage function
def analyze_best_cells_for_each_metric(matrix_file, cosmic_file, signature='SBS2'):
    """
    Analyze and plot top cells from each metric
    """
    # Load data
    print("Loading data...")
    mutation_matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)
    cosmic_signatures = pd.read_csv(cosmic_file, sep='\t', index_col='Type')
    
    # Metrics to analyze
    metrics = ['cosine', 'dot_product', 'inverse_euclidean']
    
    all_top_cells = {}
    
    for metric in metrics:
        results_file = f"sbs_similarity_analysis_{metric}_results.csv"
        print(f"\nAnalyzing {metric} results...")
        
        # Get top cells
        top_cells = extract_top_cells_from_results(results_file, signature, n_top=5)
        all_top_cells[metric] = top_cells
        
        # Plot each cell
        for i, cell in enumerate(top_cells):
            save_path = f"{signature}_{metric}_top{i+1}_{cell}.png"
            plot_cell_mutation_signature(mutation_matrix, cell, cosmic_signatures, 
                                       signature, save_path)
    
    return all_top_cells

if __name__ == "__main__":
    # Example: Plot specific cells
    matrix_file = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt"
    cosmic_file = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/COSMIC_v3.4_SBS_GRCh38.txt"
    
    # Load data
    mutation_matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)
    cosmic_signatures = pd.read_csv(cosmic_file, sep='\t', index_col='Type')
    
    # Example 1: Plot a specific cell
    # cell_to_plot = "GTCACAAGTAAAGTCA-1-SRR14340896"  # Replace with actual cell barcode
    # plot_cell_mutation_signature(mutation_matrix, cell_to_plot, cosmic_signatures, 'SBS2')
    
    # Example 2: Analyze top cells from each metric
    analyze_best_cells_for_each_metric(matrix_file, cosmic_file, 'SBS2')
    
#%% 
#!/usr/bin/env python3
"""
Aggregate Mutation Data Across Similar Cells

This script reads similarity analysis results, filters cells that meet a threshold,
aggregates their mutation profiles, and generates COSMIC-style mutational signature plots.

Author: Created for mutational signature analysis
Date: 2025-11-06
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import sys


def load_similarity_results(similarity_file):
    """
    Load similarity results from CSV file.
    
    Parameters:
    -----------
    similarity_file : str or Path
        Path to the similarity results CSV file
        
    Returns:
    --------
    pd.DataFrame
        Dataframe containing similarity results
    """
    try:
        df = pd.read_csv(similarity_file)
        print(f"Loaded similarity results from: {similarity_file}")
        print(f"Shape: {df.shape}")
        print(f"Columns: {df.columns.tolist()}")
        return df
    except Exception as e:
        print(f"Error loading similarity file: {e}")
        sys.exit(1)


def filter_by_threshold(df, signature_column, threshold=0.7):
    """
    Filter cells that meet the similarity threshold for a specific signature.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Similarity results dataframe
    signature_column : str
        Name of the column containing similarity scores (e.g., 'SBS2', 'SBS13')
    threshold : float
        Minimum similarity threshold (default: 0.7)
        
    Returns:
    --------
    pd.DataFrame
        Filtered dataframe with cells meeting threshold
    list
        List of cell barcodes that meet the threshold
    """
    # Filter rows where the similarity score exceeds threshold
    filtered_df = df[df[signature_column] > threshold].copy()
    
    # Get cell barcodes (assuming first column contains cell IDs)
    cell_barcodes = filtered_df.iloc[:, 0].tolist()
    
    print(f"\nFiltering for {signature_column} similarity > {threshold}")
    print(f"Found {len(cell_barcodes)} cells meeting threshold")
    print(f"Score range: {filtered_df[signature_column].min():.3f} - {filtered_df[signature_column].max():.3f}")
    print(f"Mean score: {filtered_df[signature_column].mean():.3f}")
    
    return filtered_df, cell_barcodes


def load_mutation_matrix(matrix_file, cell_barcodes):
    """
    Load mutation matrix and filter for specific cell barcodes.
    
    Parameters:
    -----------
    matrix_file : str or Path
        Path to the mutation matrix file
    cell_barcodes : list
        List of cell barcodes to extract
        
    Returns:
    --------
    pd.DataFrame
        Filtered mutation matrix
    """
    try:
        # Load the full mutation matrix
        full_matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)
        print(f"\nLoaded mutation matrix from: {matrix_file}")
        print(f"Full matrix shape: {full_matrix.shape}")
        
        # Filter for the specified cell barcodes
        # Cell barcodes should be in the columns
        available_barcodes = [bc for bc in cell_barcodes if bc in full_matrix.columns]
        missing_barcodes = [bc for bc in cell_barcodes if bc not in full_matrix.columns]
        
        if missing_barcodes:
            print(f"Warning: {len(missing_barcodes)} cell barcodes not found in matrix")
            if len(missing_barcodes) <= 10:
                print(f"Missing barcodes: {missing_barcodes}")
        
        filtered_matrix = full_matrix[available_barcodes]
        print(f"Filtered matrix shape: {filtered_matrix.shape}")
        print(f"Successfully extracted {len(available_barcodes)} cells")
        
        return filtered_matrix
        
    except Exception as e:
        print(f"Error loading mutation matrix: {e}")
        sys.exit(1)


def aggregate_mutations(mutation_matrix):
    """
    Sum mutations across all cells and calculate the average profile.
    
    Parameters:
    -----------
    mutation_matrix : pd.DataFrame
        Filtered mutation matrix with cells as columns
        
    Returns:
    --------
    pd.Series
        Aggregated mutation profile (average across cells)
    pd.Series
        Total mutation counts (sum across cells)
    """
    # Sum across all cells (columns)
    total_counts = mutation_matrix.sum(axis=1)
    
    # Calculate average
    avg_counts = mutation_matrix.mean(axis=1)
    
    print(f"\nAggregation statistics:")
    print(f"Total mutations across all cells: {total_counts.sum():.0f}")
    print(f"Average mutations per context: {avg_counts.mean():.2f}")
    print(f"Mutations per cell (average): {mutation_matrix.sum(axis=0).mean():.2f}")
    
    return avg_counts, total_counts


def plot_cosmic_signature(mutation_profile, output_file, title="Aggregated Mutational Profile"):
    """
    Plot mutation profile in COSMIC mutational signature format.
    
    Parameters:
    -----------
    mutation_profile : pd.Series
        Mutation counts for 96 trinucleotide contexts
    output_file : str or Path
        Path to save the output plot
    title : str
        Title for the plot
    """
    # Define colors for each mutation type (matching COSMIC style)
    colors = {
        'C>A': '#03BCEE',  # Cyan
        'C>G': '#010101',  # Black
        'C>T': '#E32926',  # Red
        'T>A': '#CAC9C9',  # Gray
        'T>C': '#A1CE63',  # Green
        'T>G': '#EBC6C4'   # Pink
    }
    
    # Extract mutation types from the index
    # Format expected: MutationType (e.g., 'A[C>A]A')
    mutation_types = []
    contexts = []
    
    for idx in mutation_profile.index:
        # Parse the mutation context string
        # Format: X[Y>Z]W where Y>Z is the mutation type
        if '[' in idx and ']' in idx:
            start = idx.index('[')
            end = idx.index(']')
            mut_type = idx[start+1:end]
            mutation_types.append(mut_type)
            contexts.append(idx)
        else:
            # Fallback if format is different
            mutation_types.append('C>A')
            contexts.append(idx)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(20, 6))
    
    # Get mutation counts
    counts = mutation_profile.values
    
    # Normalize to show relative contributions
    if counts.sum() > 0:
        normalized_counts = counts / counts.sum()
    else:
        normalized_counts = counts
    
    # Create bar positions
    x_pos = np.arange(len(counts))
    
    # Plot bars colored by mutation type
    bar_colors = [colors.get(mt, '#808080') for mt in mutation_types]
    bars = ax.bar(x_pos, normalized_counts, color=bar_colors, width=0.8, edgecolor='none')
    
    # Customize plot
    ax.set_ylabel('Relative Contribution', fontsize=12, fontweight='bold')
    ax.set_xlabel('Trinucleotide Context', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Set x-axis labels (show every 4th context to avoid crowding)
    ax.set_xticks(x_pos[::4])
    ax.set_xticklabels([contexts[i] for i in range(0, len(contexts), 4)], 
                       rotation=90, fontsize=8, ha='center')
    
    # Add grid
    ax.yaxis.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Add mutation type labels with color boxes
    unique_types = []
    seen = set()
    for mt in mutation_types:
        if mt not in seen:
            unique_types.append(mt)
            seen.add(mt)
    
    # Create legend
    legend_elements = [plt.Rectangle((0,0),1,1, fc=colors.get(mt, '#808080'), 
                                    edgecolor='black', label=mt) 
                      for mt in unique_types]
    ax.legend(handles=legend_elements, loc='upper right', ncol=6, 
             frameon=True, fontsize=10)
    
    # Add dividing lines between mutation types
    type_boundaries = []
    current_type = mutation_types[0]
    for i, mt in enumerate(mutation_types[1:], 1):
        if mt != current_type:
            type_boundaries.append(i - 0.5)
            current_type = mt
    
    for boundary in type_boundaries:
        ax.axvline(x=boundary, color='black', linewidth=0.5, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved plot to: {output_file}")
    plt.close()


def save_cell_barcode_list(cell_barcodes, output_file, similarity_metric):
    """
    Save list of cell barcodes that met the filtering threshold.
    
    Parameters:
    -----------
    cell_barcodes : list
        List of cell barcodes
    output_file : str or Path
        Path to save the output file
    similarity_metric : str
        Name of the similarity metric used
    """
    df = pd.DataFrame({
        'cell_barcode': cell_barcodes,
        'similarity_metric': similarity_metric
    })
    
    df.to_csv(output_file, index=False)
    print(f"Saved cell barcode list to: {output_file}")


def aggregate_similar_cells(similarity_file, mutation_matrix_file, 
                           signature_column='SBS2', threshold=0.7, 
                           output_dir=None):
    """
    Main function to aggregate mutation data across cells that meet similarity threshold.
    
    Parameters:
    -----------
    similarity_file : str or Path
        Path to similarity results CSV file
    mutation_matrix_file : str or Path
        Path to the mutation matrix file
    signature_column : str
        Name of the signature column to filter on (default: 'SBS2')
    threshold : float
        Minimum similarity threshold (default: 0.7)
    output_dir : str or Path, optional
        Output directory for results (default: same as similarity file)
    
    Returns:
    --------
    dict
        Dictionary containing results including aggregated profile and cell list
    """
    # Setup output directory
    similarity_path = Path(similarity_file)
    if output_dir is None:
        output_dir = similarity_path.parent
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract similarity metric name from filename
    similarity_metric = similarity_path.stem.replace('sbs_similarity_analysis_', '').replace('_results', '')
    
    print(f"\n{'='*80}")
    print(f"AGGREGATING SIMILAR CELLS ANALYSIS")
    print(f"Similarity Metric: {similarity_metric}")
    print(f"Signature: {signature_column}")
    print(f"Threshold: {threshold}")
    print(f"{'='*80}")
    
    # Step 1: Load and filter similarity results
    similarity_df = load_similarity_results(similarity_file)
    filtered_df, cell_barcodes = filter_by_threshold(similarity_df, signature_column, threshold)
    
    if len(cell_barcodes) == 0:
        print(f"\nWarning: No cells met the threshold of {threshold} for {signature_column}")
        print("Try lowering the threshold or checking your similarity results.")
        return None
    
    # Step 2: Load and filter mutation matrix
    mutation_matrix = load_mutation_matrix(mutation_matrix_file, cell_barcodes)
    
    if mutation_matrix.shape[1] == 0:
        print("\nError: No cells found in mutation matrix")
        return None
    
    # Step 3: Aggregate mutations
    avg_profile, total_profile = aggregate_mutations(mutation_matrix)
    
    # Step 4: Generate plot
    plot_title = (f"Aggregated Mutational Profile\n"
                 f"{signature_column} - {similarity_metric} - "
                 f"Threshold: {threshold} - N={len(mutation_matrix.columns)} cells")
    
    plot_file = output_dir / f"aggregated_profile_{similarity_metric}_{signature_column}_threshold{threshold}.png"
    plot_cosmic_signature(avg_profile, plot_file, title=plot_title)
    
    # Step 5: Save cell barcode list
    barcode_file = output_dir / f"filtered_cells_{similarity_metric}_{signature_column}_threshold{threshold}.csv"
    save_cell_barcode_list(cell_barcodes, barcode_file, similarity_metric)
    
    # Step 6: Save aggregated profile as CSV
    profile_file = output_dir / f"aggregated_profile_{similarity_metric}_{signature_column}_threshold{threshold}.csv"
    profile_df = pd.DataFrame({
        'mutation_context': avg_profile.index,
        'average_count': avg_profile.values,
        'total_count': total_profile.values
    })
    profile_df.to_csv(profile_file, index=False)
    print(f"Saved aggregated profile to: {profile_file}")
    
    print(f"\n{'='*80}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nOutput files:")
    print(f"  - Plot: {plot_file.name}")
    print(f"  - Cell list: {barcode_file.name}")
    print(f"  - Profile data: {profile_file.name}")
    
    return {
        'filtered_cells': cell_barcodes,
        'avg_profile': avg_profile,
        'total_profile': total_profile,
        'filtered_df': filtered_df,
        'mutation_matrix': mutation_matrix,
        'output_files': {
            'plot': plot_file,
            'cell_list': barcode_file,
            'profile': profile_file
        }
    }


def process_all_similarity_files(results_dir, mutation_matrix_file, 
                                 signature='SBS2', threshold=0.7,
                                 output_dir=None):
    """
    Process all similarity result files in a directory.
    
    Parameters:
    -----------
    results_dir : str or Path
        Directory containing similarity results files
    mutation_matrix_file : str or Path
        Path to the mutation matrix file
    signature : str
        Signature column to analyze (default: 'SBS2')
    threshold : float
        Similarity threshold (default: 0.7)
    output_dir : str or Path, optional
        Output directory (default: same as results_dir)
    
    Returns:
    --------
    dict
        Dictionary with results for each similarity metric
    """
    results_dir = Path(results_dir)
    
    # Find all similarity result files
    similarity_files = [
        results_dir / 'sbs_similarity_analysis_cosine_results.csv',
        results_dir / 'sbs_similarity_analysis_dot_product_results.csv',
        results_dir / 'sbs_similarity_analysis_inverse_euclidean_results.csv'
    ]
    
    results = {}
    
    for sim_file in similarity_files:
        if sim_file.exists():
            print(f"\n\nProcessing: {sim_file.name}")
            result = aggregate_similar_cells(
                sim_file, 
                mutation_matrix_file,
                signature_column=signature,
                threshold=threshold,
                output_dir=output_dir
            )
            
            if result:
                metric_name = sim_file.stem.replace('sbs_similarity_analysis_', '').replace('_results', '')
                results[metric_name] = result
        else:
            print(f"\nWarning: File not found: {sim_file}")
    
    return results


def main():
    """Command-line interface for the script."""
    parser = argparse.ArgumentParser(
        description='Aggregate mutation data across cells meeting similarity threshold',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single similarity file
  python aggregate_similar_cells.py \\
    --similarity sbs_similarity_analysis_cosine_results.csv \\
    --matrix ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt \\
    --signature SBS2 \\
    --threshold 0.7

  # Process all similarity files in a directory
  python aggregate_similar_cells.py \\
    --results-dir /path/to/results_NMF_v0.1.1 \\
    --matrix ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt \\
    --signature SBS2 \\
    --threshold 0.7 \\
    --process-all

  # Use different signature and threshold
  python aggregate_similar_cells.py \\
    --similarity sbs_similarity_analysis_dot_product_results.csv \\
    --matrix ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt \\
    --signature SBS13 \\
    --threshold 0.6
        """
    )
    
    # Input arguments
    parser.add_argument('--similarity', type=str,
                       help='Path to similarity results CSV file')
    parser.add_argument('--matrix', type=str, required=True,
                       help='Path to mutation matrix file (required)')
    parser.add_argument('--results-dir', type=str,
                       help='Directory containing all similarity result files')
    
    # Analysis parameters
    parser.add_argument('--signature', type=str, default='SBS2',
                       help='Signature column to analyze (default: SBS2)')
    parser.add_argument('--threshold', type=float, default=0.7,
                       help='Minimum similarity threshold (default: 0.7)')
    
    # Output arguments
    parser.add_argument('--output-dir', type=str,
                       help='Output directory (default: same as input)')
    parser.add_argument('--process-all', action='store_true',
                       help='Process all similarity files in results-dir')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.process_all and not args.similarity:
        parser.error("Either --similarity or --process-all (with --results-dir) is required")
    
    if args.process_all and not args.results_dir:
        parser.error("--results-dir is required when using --process-all")
    
    # Process files
    if args.process_all:
        print("Processing all similarity files...")
        results = process_all_similarity_files(
            args.results_dir,
            args.matrix,
            signature=args.signature,
            threshold=args.threshold,
            output_dir=args.output_dir
        )
        
        print(f"\n\n{'='*80}")
        print("SUMMARY OF ALL ANALYSES")
        print(f"{'='*80}")
        for metric, result in results.items():
            if result:
                print(f"\n{metric.upper()}:")
                print(f"  Cells meeting threshold: {len(result['filtered_cells'])}")
                print(f"  Total mutations: {result['total_profile'].sum():.0f}")
                print(f"  Avg mutations/cell: {result['mutation_matrix'].sum(axis=0).mean():.2f}")
    
    else:
        print("Processing single similarity file...")
        result = aggregate_similar_cells(
            args.similarity,
            args.matrix,
            signature_column=args.signature,
            threshold=args.threshold,
            output_dir=args.output_dir
        )
        
        if not result:
            sys.exit(1)

#%%
# Example 1: Analyze a single similarity file
# ===========================================
print("Example 1: Single similarity file analysis")
print("=" * 80)

result = aggregate_similar_cells(
    similarity_file='/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/sbs_similarity_analysis_inverse_euclidean_results.csv',
    mutation_matrix_file='/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt',
    signature_column='SBS2',
    threshold=0.6,
    output_dir='/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/aggregated_analysis'
)

if result:
    print(f"\nAnalysis successful!")
    print(f"Number of cells: {len(result['filtered_cells'])}")
    print(f"Total mutations: {result['total_profile'].sum():.0f}")
    print(f"\nFirst 5 cell barcodes:")
    for bc in result['filtered_cells'][:5]:
        print(f"  - {bc}")
    print(f"\nOutput files saved to: {result['output_files']['plot'].parent}")

#%%
# ==========================================
print("\n\n" + "=" * 80)
print("Example 2: Process all similarity metrics")
print("=" * 80)

results = process_all_similarity_files(
    results_dir='/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1',
    mutation_matrix_file='/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt',
    signature='SBS2',
    threshold=0.6,
    output_dir='/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/aggregated_analysis'
)

# Compare results across metrics
print("\n\nComparison across metrics:")
print("-" * 80)
print(f"{'Metric':<20} {'N Cells':<10} {'Total Muts':<15} {'Avg Muts/Cell':<15}")
print("-" * 80)

for metric_name, result in results.items():
    if result:
        n_cells = len(result['filtered_cells'])
        total_muts = result['total_profile'].sum()
        avg_muts = result['mutation_matrix'].sum(axis=0).mean()
        print(f"{metric_name:<20} {n_cells:<10} {total_muts:<15.0f} {avg_muts:<15.2f}")

#%%
import pandas as pd

# Read barcodes
barcodes = pd.read_csv('/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/aggregated_analysis/filtered_cells_inverse_euclidean_SBS2_threshold0.6.csv')['cell_barcode'].tolist()

# Read and filter matrix
matrix = pd.read_csv('/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt', sep='\t', index_col=0)
filtered_matrix = matrix[barcodes]
#%%
# Save
filtered_matrix.to_csv('/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ED_0.6_Filtered_Basal_Cell_SNP_matrix_for_SigProfiler.txt', sep='\t')

#%%
inverse_ED_basal_cells = pd.read_csv('/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/sbs_similarity_analysis_inverse_euclidean_results.csv')
# Sort the dataframe
inverse_ED_basal_cells = inverse_ED_basal_cells.sort_values('SBS2', ascending=False)
top_df = inverse_ED_basal_cells.head(472)
bottom_df = inverse_ED_basal_cells.tail(472)
top_barcodes = top_df['Unnamed: 0'].tolist()  
bottom_barcodes = bottom_df['Unnamed: 0'].tolist()

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
import scipy
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

adata_pp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_pp_CD.h5ad')
)
#%% Do some QC to check how this looks on our single cell data
sc.pl.umap(adata_pp, color=['final_annotation'], size=5, frameon=False)
result = pd.read_csv(os.path.join(working_dir, "sigprofiler_mutation_counts.txt"), sep= '\t')
result
mutations_per_cell = result.sum(axis=0)

mutations_per_cell_df = mutations_per_cell.to_frame(name='Total_Mutations')
mutations_per_cell_df
# Make sure the cell barcodes match
mutations_per_cell_df = mutations_per_cell_df.reindex(adata_pp.obs.index.tolist())
mutations_per_cell_df['Total_Mutations'] = mutations_per_cell_df['Total_Mutations'].fillna(0)


# Add to AnnData obs
adata_pp.obs['total_mutations'] = mutations_per_cell_df['Total_Mutations']
adata_pp
sc.pl.umap(adata_pp, color=['CytoTRACE2_Score'], size=5, frameon=False, cmap='plasma')
sc.pl.umap(adata_pp, color=['cnv_score'], size=5, frameon=False, cmap = 'plasma')
sc.pl.umap(adata_pp, color=['Final_cancer_cell_status'], size=5, frameon=False)

#%% Load the cell annotations file
cell_annotations_path = os.path.join(working_dir, 'cell_annotations.txt')
cell_annotations = pd.read_csv(cell_annotations_path, sep='\t')

# Extract the complete list of cell barcodes
all_expected_barcodes = set(cell_annotations['cell_barcodes'])

print(f"Total expected cells: {len(all_expected_barcodes)}")
#%% Load the combined callable sites file
callable_sites_path = os.path.join(working_dir, 'CombinedCallableSites', 'complete_callable_sites.tsv')
callable_sites_df = pd.read_csv(callable_sites_path, sep='\t')

# Get the cell barcodes from callable sites
callable_barcodes = set(callable_sites_df['CB'])

print(f"Cells in callable sites: {len(callable_barcodes)}")

# Find missing cells
len(list(set(callable_barcodes).intersection(set(mutations_per_cell[1:].index))))
len(list(set(mutations_per_cell[1:].index) - set(callable_barcodes)))

missing_cells = all_expected_barcodes - callable_barcodes
print(f"Cells missing from callable sites: {len(missing_cells)}")
missing_cells
# Check if these cells are in the mutation data

missing_intersection = list(set(missing_cells).intersection(set(mutations_per_cell[1:].index)))
print(f"Missing cells also absent from mutations: {len(missing_intersection)}")
len(set(missing_intersection).intersection(list(set(mutations_per_cell[1:].index) - set(callable_barcodes))))

#%% Okay I found two ways to get out barcodes where we have mutations recorded 
# but the cells don't meet the cutoff that is used for the count base function. 
# So this will be an issue and I need to go and modify the snp call function to meet the same strengency that a cell has to have at least 5x coverage at the 
# Before I increase the stringency I want to see how the data looks for the SNP calls as is. In the fututre I'm going to drop any cells where mutations were found but those mutations sites weren't > 5X in coverage

cells_to_drop = list(set(mutations_per_cell[1:].index) - set(callable_barcodes))

adata_pp.obs.loc[adata_pp.obs.index.isin(cells_to_drop), 'total_mutations'] = 0

# Identify columns to zero out
cols_to_zero = result.columns[1:].intersection(cells_to_drop)
result.loc[:, cols_to_zero] = 0
result.index = result['Unnamed: 0']
result = result.drop(columns="Unnamed: 0")
result
result.index.names = ['']

result.to_csv(os.path.join(working_dir, 'SNP_matrix_for_SigProfiler.txt'), sep='\t', index=True, header=True)
#%%
# Check to see if you had any changes at this point you should have 0 if eveything is working
from pandas.testing import assert_frame_equal
# Read back in the file you just made 
result_2 = pd.read_csv(os.path.join(working_dir, "SNP_matrix_for_SigProfiler.txt"), sep= '\t')
result_2
result_2.index = result_2['Unnamed: 0']
result_2 = result_2.drop(columns="Unnamed: 0")
result_2
result_2.index.names = ['']
assert_frame_equal(result, result_2)
# If the assert function threw an error stop and try to figure out what is going on
# as long as there was no differences you can use the results object if needed you can use the results_2 object to keep pushing the pipeline
result.index
result.columns
callable_sites_df.index = callable_sites_df['CB']
callable_sites_df = callable_sites_df.reindex(adata_pp.obs.index.tolist())
callable_sites_df['CB'] = list(callable_sites_df.index)
callable_sites_df['SitesPerCell'] = callable_sites_df['SitesPerCell'].fillna(0)

callable_sites_df

#%% Plot some UMAPs 
# Go ahead and add the total mutations to the adata for both the raw and then the normalized mutation counts 
adata_pp.obs['normalized_total_mutations'] = adata_pp.obs['total_mutations'] / callable_sites_df['SitesPerCell']
adata_pp.obs['normalized_total_mutations'] = adata_pp.obs['normalized_total_mutations'].fillna(0)

sc.pl.umap(adata_pp, color=['total_mutations'], size=10, cmap = 'plasma', frameon=False)
sc.pl.umap(adata_pp, color=['normalized_total_mutations'], size=15, cmap = 'plasma', frameon=False)

#%% This was for modaheseh to run network analysis you can skip 
adata_pp.var.index = adata_pp.var['gene_ids']
adata_pp.var
# Subset the data to send over to mohadeseh for differential network analysis
basal_cell_high_SBS2_cor = adata_pp[adata_pp.obs.index.isin(top_barcodes)]
basal_cell_high_SBS2_cor.obs['final_annotation']

# Prepare and save gene expression matrix
adata_X_df_T_high_SBS2_corr = pd.DataFrame.sparse.from_spmatrix(
    scipy.sparse.csr_matrix(basal_cell_high_SBS2_cor.X),
    index=basal_cell_high_SBS2_cor.obs.index,
    columns=basal_cell_high_SBS2_cor.var['gene_symbol'].index
).T
    
expression_path = os.path.join(working_dir, f"Basal_Cells_High_SBS2_ED_Cor_0.6.tsv")
adata_X_df_T_high_SBS2_corr.to_csv(
    expression_path,
    sep="\t",
    header=True,
    index=True,
    chunksize=10000
)

basal_cell_low_SBS2_cor = adata_pp[adata_pp.obs.index.isin(bottom_barcodes)]
basal_cell_low_SBS2_cor.obs['final_annotation']

adata_X_df_T_low_SBS2_corr = pd.DataFrame.sparse.from_spmatrix(
    scipy.sparse.csr_matrix(basal_cell_low_SBS2_cor.X),
    index=basal_cell_low_SBS2_cor.obs.index,
    columns=basal_cell_low_SBS2_cor.var['gene_symbol'].index
).T
    
expression_path = os.path.join(working_dir, f"Basal_Cells_Low_SBS2_ED_Cor_0.6.tsv")
adata_X_df_T_high_SBS2_corr.to_csv(
    expression_path,
    sep="\t",
    header=True,
    index=True,
    chunksize=10000
)


#%%
import matplotlib.pyplot as plt
import numpy as np

# Subset the data
basal_cell = adata_pp[adata_pp.obs['final_annotation'] == "basal cell"]
cancer_cells = basal_cell[
    (basal_cell.obs['Final_cancer_cell_status'] == 'Cancer cell') & 
    (basal_cell.obs['source_name'] == 'head and neck squamous cell carcinoma')
]
normal_cells = basal_cell[
    (basal_cell.obs['Final_cancer_cell_status'] == 'Normal cell') & 
    (basal_cell.obs['source_name'] == 'normal tissue adjucent to head and neck squamous cell carcinoma')
]

import anndata as ad
combined_adata = ad.concat([cancer_cells, normal_cells], join='outer', merge='same')

#%% Calculate the neighborhood plot and cluster
from os import cpu_count

NCPUS = cpu_count()
sc.pp.pca(combined_adata, n_comps=NCPUS)
sc.pp.neighbors(combined_adata, n_pcs=NCPUS)
sc.tl.umap(combined_adata)
sc.tl.leiden(combined_adata, key_added='cancer_clusters', resolution=1, random_state=42)

#%% Visulization of the updated UMAP plot
plt.figure(figsize=(4, 4))
sc.pl.umap(combined_adata, color=['Final_cancer_cell_status'], size=5, frameon=False)

sc.pl.umap(adata_pp, color=['Final_cancer_cell_status'], size=5, frameon=False)
sc.pl.umap(adata_pp, color=['final_annotation'], size=5, frameon=False)
adata_pp
sc.pl.umap(adata_pp, color=['CytoTRACE2_Score'], size=5, frameon=False, cmap='plasma')
sc.pl.umap(adata_pp, color=['cnv_score'], size=5, frameon=False, cmap='plasma')
#%%
import pandas as pd
import numpy as np
import scanpy as sc

# Create the SBS2_corr column in adata_pp.obs
adata_pp.obs['SBS2_corr'] = np.nan  # Initialize with NA values

# Set 'High' for top barcodes
adata_pp.obs.loc[adata_pp.obs.index.isin(top_barcodes), 'SBS2_corr'] = 'High'

# Set 'Low' for bottom barcodes
adata_pp.obs.loc[adata_pp.obs.index.isin(bottom_barcodes), 'SBS2_corr'] = 'Low'

# Convert to categorical (optional but can be helpful)
adata_pp.obs['SBS2_corr'] = adata_pp.obs['SBS2_corr'].astype('category')

# Set up the color palette
custom_palette = {
    'High': 'red',
    'Low': 'green'
}

# Plot the UMAP with custom colors
sc.pl.umap(
    adata_pp, 
    color=['SBS2_corr'], 
    size=5, 
    frameon=False,
    palette=custom_palette
)
#%% That ended the preprocessing steps if we want to select cellsusing high Euclidian Distance similarity to SBS2 now we will work on if we want to use NMF to pull out signature
# Keep in mind for this section we would like to set up a method were we can subset cells above a cetrain threshold of mutations starting with > 20 
# After that we need to calculate some internal control metrics to make sure we are selecting the right number of componetes for the signals since we are forcing positivity in a signal in each componet.
# Before running NMF I need a PCA of the filtered dataset with the martix represented as bar plot that can have positive or negative values and the eigen vector displayed as a scree plot
# Durring the NMF I would like the signatures to be compared to COSMIC mutation signatures outputting the cosine, the dot produt, and the inverse Euclidian distance simialrity scores
# After running NMF I need to check the number of componets using a Frobenius graph and an average silhouette score

#%% Okay I spoke with Diako today and he wants to look at high mutation burden cells and run NMF so I'm going to look at the distribution of cells starting with 15 or more mutations and then try to pull the SBS2 mutation out 

#!/usr/bin/env python3
"""
Mutation Count Distribution Analysis for NMF Preprocessing

This script analyzes your sparse mutation matrix to help identify cells
with sufficient mutation counts for effective NMF signature extraction.

Key outputs:
1. Distribution statistics of mutations per cell
2. Breakdown by mutation count thresholds
3. Sparsity analysis at different thresholds
4. Visualizations to guide threshold selection
5. Filtered matrix export for NMF analysis

Run interactively in your Python IDE (Spyder, PyCharm, VSCode, etc.)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set better plotting defaults
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

#%% Configuration
# Update this path if needed
MATRIX_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt"
OUTPUT_DIR = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/"

#%% Load Data
print("="*80)
print("LOADING MUTATION MATRIX")
print("="*80)

# Load the mutation matrix
data = pd.read_csv(MATRIX_FILE, sep='\t', index_col=0)
print(f"\nOriginal matrix shape: {data.shape}")
print(f"Mutation contexts (rows): {data.shape[0]}")
print(f"Cells/samples (columns): {data.shape[1]}")

#%% Calculate Mutations Per Cell
print("\n" + "="*80)
print("CALCULATING MUTATIONS PER CELL")
print("="*80)

mutations_per_cell = data.sum(axis=0)

# print(f"\nTotal mutations across all cells: {mutations_per_cell.sum():,.0f}")
# print(f"Mean mutations per cell: {mutations_per_cell.mean():.2f}")
# print(f"Median mutations per cell: {mutations_per_cell.median():.0f}")
# print(f"Std deviation: {mutations_per_cell.std():.2f}")
# print(f"Min mutations: {mutations_per_cell.min():.0f}")
# print(f"Max mutations: {mutations_per_cell.max():.0f}")

# #%% Quartile Analysis
# quartiles = mutations_per_cell.quantile([0.25, 0.5, 0.75, 0.90, 0.95, 0.99])
# print("\nQuartile Analysis:")
# print(f"  25th percentile: {quartiles[0.25]:.0f} mutations")
# print(f"  50th percentile: {quartiles[0.50]:.0f} mutations")
# print(f"  75th percentile: {quartiles[0.75]:.0f} mutations")
# print(f"  90th percentile: {quartiles[0.90]:.0f} mutations")
# print(f"  95th percentile: {quartiles[0.95]:.0f} mutations")
# print(f"  99th percentile: {quartiles[0.99]:.0f} mutations")

# #%% Threshold Analysis - Focus on ≥10 mutations
# print("\n" + "="*80)
# print("THRESHOLD ANALYSIS")
# print("="*80)

# # Analyze different thresholds starting at 10
# thresholds = [10, 15, 20, 25, 30, 40, 50]

# results = []
# for thresh in thresholds:
#     cells_above = mutations_per_cell >= thresh
#     n_cells = cells_above.sum()
    
#     if n_cells == 0:
#         continue
    
#     # Get subset statistics
#     subset = mutations_per_cell[cells_above]
#     filtered_data = data.loc[:, cells_above]
    
#     # Calculate sparsity (percentage of zeros)
#     sparsity = (filtered_data == 0).sum().sum() / (filtered_data.shape[0] * filtered_data.shape[1])
    
#     # Average mutations per context (important for NMF)
#     avg_per_context = subset.mean() / data.shape[0]  # mutations / 96 contexts
    
#     # NMF feasibility assessment
#     if avg_per_context > 0.5:
#         feasibility = "Excellent"
#     elif avg_per_context > 0.3:
#         feasibility = "Good"
#     elif avg_per_context > 0.15:
#         feasibility = "Moderate"
#     else:
#         feasibility = "Challenging"
    
#     results.append({
#         'threshold': f'≥{thresh}',
#         'n_cells': n_cells,
#         'pct_cells': 100 * n_cells / len(mutations_per_cell),
#         'mean_muts': subset.mean(),
#         'median_muts': subset.median(),
#         'total_muts': subset.sum(),
#         'sparsity_pct': 100 * sparsity,
#         'muts_per_context': avg_per_context,
#         'nmf_feasibility': feasibility
#     })

# df_thresholds = pd.DataFrame(results)
# print("\nThreshold Comparison:")
# print(df_thresholds.to_string(index=False))

# # Save to file
# output_file = Path(OUTPUT_DIR) / "threshold_analysis.txt"
# df_thresholds.to_csv(output_file, sep='\t', index=False)
# print(f"\nSaved threshold analysis to: {output_file}")

# #%% Cells with ≥10 mutations - Detailed Analysis
# print("\n" + "="*80)
# print("DETAILED ANALYSIS: CELLS WITH ≥10 MUTATIONS")
# print("="*80)

# cells_10plus = mutations_per_cell >= 20
# n_cells_10plus = cells_10plus.sum()

# print(f"\nCells with ≥10 mutations: {n_cells_10plus:,} ({100*n_cells_10plus/len(mutations_per_cell):.2f}%)")

# if n_cells_10plus > 0:
#     subset_10plus = mutations_per_cell[cells_10plus]
    
#     print(f"\nStatistics for these {n_cells_10plus:,} cells:")
#     print(f"  Mean: {subset_10plus.mean():.2f} mutations")
#     print(f"  Median: {subset_10plus.median():.0f} mutations")
#     print(f"  Std Dev: {subset_10plus.std():.2f}")
#     print(f"  Min: {subset_10plus.min():.0f}")
#     print(f"  Max: {subset_10plus.max():.0f}")
#     print(f"  Total mutations: {subset_10plus.sum():,.0f}")
    
#     # Show distribution within this subset
#     print("\nDistribution within ≥10 mutation cells:")
#     bins = [10, 15, 20, 30, 50, 100, subset_10plus.max()+1]
#     for i in range(len(bins)-1):
#         count = ((subset_10plus >= bins[i]) & (subset_10plus < bins[i+1])).sum()
#         pct = 100 * count / n_cells_10plus
#         print(f"  {bins[i]}-{bins[i+1]-1} mutations: {count:,} cells ({pct:.1f}%)")

# #%% Visualization 1: Overall Distribution
# print("\n" + "="*80)
# print("CREATING VISUALIZATIONS")
# print("="*80)

# fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# # Plot 1: Histogram of all cells
# ax1 = axes[0, 0]
# ax1.hist(mutations_per_cell, bins=50, edgecolor='black', alpha=0.7)
# ax1.axvline(20, color='red', linestyle='--', linewidth=2, label='20 mutation threshold')
# ax1.set_xlabel('Mutations per Cell', fontsize=12)
# ax1.set_ylabel('Number of Cells', fontsize=12)
# ax1.set_title('Distribution of Mutations Across All Cells', fontsize=14, fontweight='bold')
# ax1.legend()
# ax1.grid(True, alpha=0.3)

# # Plot 2: Log scale histogram
# ax2 = axes[0, 1]
# ax2.hist(mutations_per_cell, bins=50, edgecolor='black', alpha=0.7)
# ax2.axvline(20, color='red', linestyle='--', linewidth=2, label='20 mutation threshold')
# ax2.set_xlabel('Mutations per Cell', fontsize=12)
# ax2.set_ylabel('Number of Cells (log scale)', fontsize=12)
# ax2.set_yscale('log')
# ax2.set_title('Distribution of Mutations (Log Scale)', fontsize=14, fontweight='bold')
# ax2.legend()
# ax2.grid(True, alpha=0.3)

# # Plot 3: Zoomed in on ≥10 mutations
# if n_cells_10plus > 0:
#     ax3 = axes[1, 0]
#     ax3.hist(subset_10plus, bins=30, edgecolor='black', alpha=0.7, color='green')
#     ax3.set_xlabel('Mutations per Cell', fontsize=12)
#     ax3.set_ylabel('Number of Cells', fontsize=12)
#     ax3.set_title(f'Distribution for Cells with ≥20 Mutations (n={n_cells_10plus:,})', 
#                   fontsize=14, fontweight='bold')
#     ax3.grid(True, alpha=0.3)

# # Plot 4: Threshold comparison
# ax4 = axes[1, 1]
# thresholds_plot = df_thresholds['threshold'].str.replace('≥', '').astype(int)
# ax4.plot(thresholds_plot, df_thresholds['n_cells'], 'o-', linewidth=2, markersize=8, label='Number of cells')
# ax4_twin = ax4.twinx()
# ax4_twin.plot(thresholds_plot, df_thresholds['muts_per_context'], 's-', 
#               color='orange', linewidth=2, markersize=8, label='Mutations/context')
# ax4.set_xlabel('Minimum Mutation Threshold', fontsize=12)
# ax4.set_ylabel('Number of Cells', fontsize=12, color='blue')
# ax4_twin.set_ylabel('Average Mutations per Context', fontsize=12, color='orange')
# ax4.set_title('Cells vs Threshold & NMF Feasibility', fontsize=14, fontweight='bold')
# ax4.tick_params(axis='y', labelcolor='blue')
# ax4_twin.tick_params(axis='y', labelcolor='orange')
# ax4.grid(True, alpha=0.3)
# ax4.legend(loc='upper left')
# ax4_twin.legend(loc='upper right')

# plt.tight_layout()
# fig_path = Path(OUTPUT_DIR) / "mutation_distribution_analysis.png"
# plt.savefig(fig_path, dpi=300, bbox_inches='tight')
# print(f"\nSaved figure to: {fig_path}")
# plt.show()

# #%% Sparsity Analysis Visualization
# fig2, ax = plt.subplots(figsize=(12, 6))

# thresholds_plot = df_thresholds['threshold'].str.replace('≥', '').astype(int)
# ax.bar(thresholds_plot, df_thresholds['sparsity_pct'], 
#        color='steelblue', edgecolor='black', alpha=0.7)
# ax.set_xlabel('Minimum Mutation Threshold', fontsize=12)
# ax.set_ylabel('Sparsity (%)', fontsize=12)
# ax.set_title('Matrix Sparsity by Mutation Threshold', fontsize=14, fontweight='bold')
# ax.grid(True, alpha=0.3, axis='y')

# # Add feasibility colors as background
# for i, row in df_thresholds.iterrows():
#     thresh_val = int(row['threshold'].replace('≥', ''))
#     color = {'Excellent': 'green', 'Good': 'yellow', 
#              'Moderate': 'orange', 'Challenging': 'red'}[row['nmf_feasibility']]
#     ax.axvspan(thresh_val-2.5, thresh_val+2.5, alpha=0.1, color=color)

# plt.tight_layout()
# fig2_path = Path(OUTPUT_DIR) / "sparsity_analysis.png"
# plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
# print(f"Saved figure to: {fig2_path}")
# plt.show()

#%%
#!/usr/bin/env python
"""
Comprehensive Single-Cell Mutation Signature Extraction Pipeline

This pipeline integrates all lessons learned from sparse single-cell mutation data:
1. Mutation count filtering (preprocessing)
2. PCA analysis for pattern detection
3. Enhanced NMF with multiple similarity metrics
4. Comprehensive quality assessment (Frobenius, Silhouette)
5. COSMIC signature matching

Author: Jake Lehle
Date: 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Core analysis libraries
from sklearn.decomposition import NMF, PCA
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cosine, euclidean
from scipy.optimize import linear_sum_assignment, nnls

# Set plotting style
plt.style.use('default')
sns.set_palette("husl")

# COSMIC signature colors (standard)
COSMIC_COLORS = {
    'C>A': '#1EBFF0',
    'C>G': '#050708',
    'C>T': '#E62725',
    'T>A': '#CBCACB',
    'T>C': '#A1CE63',
    'T>G': '#EDB6C2'
}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def log(msg, level="INFO"):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {level}: {msg}", flush=True)


def save_figure(fig, output_dir, filename, dpi=300):
    """Save figure to output directory"""
    filepath = Path(output_dir) / filename
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
    log(f"Saved: {filepath}")
    return filepath


# =============================================================================
# STEP 1: DATA LOADING AND FILTERING
# =============================================================================

def load_and_filter_mutation_matrix(matrix_file, mutation_threshold=20, 
                                   output_dir=None, save_filtered=True):
    """
    Load mutation matrix and filter cells by mutation count
    
    Parameters:
    -----------
    matrix_file : str
        Path to mutation matrix (96 contexts × cells)
    mutation_threshold : int
        Minimum mutations per cell to retain
    output_dir : str, optional
        Directory to save filtered matrix
    save_filtered : bool
        Whether to save filtered matrix
    
    Returns:
    --------
    dict with:
        - 'original_matrix': pd.DataFrame, original data
        - 'filtered_matrix': pd.DataFrame, filtered data
        - 'mutations_per_cell': pd.Series, mutation counts
        - 'stats': dict, filtering statistics
    """
    log("="*80)
    log("STEP 1: LOADING AND FILTERING MUTATION MATRIX")
    log("="*80)
    
    # Load data
    log(f"Loading: {matrix_file}")
    data_original = pd.read_csv(matrix_file, sep='\t', index_col=0)
    log(f"Original matrix: {data_original.shape[0]} contexts × {data_original.shape[1]} cells")
    
    # Calculate mutations per cell
    mutations_per_cell = data_original.sum(axis=0)
    
    # Original statistics
    log(f"\nOriginal dataset statistics:")
    log(f"  Total cells: {len(mutations_per_cell):,}")
    log(f"  Total mutations: {int(mutations_per_cell.sum()):,}")
    log(f"  Mean mutations/cell: {mutations_per_cell.mean():.2f}")
    log(f"  Median mutations/cell: {mutations_per_cell.median():.0f}")
    log(f"  Range: {mutations_per_cell.min():.0f} - {mutations_per_cell.max():.0f}")
    
    # Filter by threshold
    log(f"\nFiltering cells with ≥{mutation_threshold} mutations...")
    cells_pass_filter = mutations_per_cell >= mutation_threshold
    data_filtered = data_original.loc[:, cells_pass_filter]
    mutations_filtered = mutations_per_cell[cells_pass_filter]
    
    # Filtered statistics
    n_retained = data_filtered.shape[1]
    pct_retained = 100 * n_retained / data_original.shape[1]
    
    log(f"\nFiltered dataset statistics:")
    log(f"  Cells retained: {n_retained:,} ({pct_retained:.2f}%)")
    log(f"  Cells removed: {data_original.shape[1] - n_retained:,}")
    log(f"  Mean mutations/cell: {mutations_filtered.mean():.2f}")
    log(f"  Median mutations/cell: {mutations_filtered.median():.0f}")
    log(f"  Total mutations: {int(mutations_filtered.sum()):,}")
    
    # Calculate sparsity metrics
    sparsity = (data_filtered == 0).sum().sum() / (data_filtered.shape[0] * data_filtered.shape[1])
    avg_per_context = mutations_filtered.mean() / data_filtered.shape[0]
    
    log(f"\nData quality metrics:")
    log(f"  Sparsity: {100*sparsity:.2f}%")
    log(f"  Avg mutations/context: {avg_per_context:.3f}")
    
    # NMF feasibility assessment
    if avg_per_context > 0.5:
        feasibility = "EXCELLENT"
    elif avg_per_context > 0.3:
        feasibility = "GOOD"
    elif avg_per_context > 0.15:
        feasibility = "MODERATE"
    else:
        feasibility = "CHALLENGING"
    log(f"  NMF Feasibility: {feasibility}")
    
    # Save filtered matrix if requested
    if save_filtered and output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        filtered_file = output_path / f"filtered_matrix_min{mutation_threshold}muts.txt"
        data_filtered.to_csv(filtered_file, sep='\t', float_format='%.0f')
        log(f"\nSaved filtered matrix: {filtered_file}")
        
        # Save filtering stats
        stats_file = output_path / f"filtering_stats_min{mutation_threshold}muts.txt"
        with open(stats_file, 'w') as f:
            f.write(f"Mutation Matrix Filtering Statistics\n")
            f.write(f"{'='*60}\n\n")
            f.write(f"Threshold: ≥{mutation_threshold} mutations per cell\n\n")
            f.write(f"Original:\n")
            f.write(f"  Cells: {data_original.shape[1]:,}\n")
            f.write(f"  Total mutations: {int(mutations_per_cell.sum()):,}\n")
            f.write(f"  Mean/cell: {mutations_per_cell.mean():.2f}\n\n")
            f.write(f"Filtered:\n")
            f.write(f"  Cells: {n_retained:,} ({pct_retained:.2f}%)\n")
            f.write(f"  Total mutations: {int(mutations_filtered.sum()):,}\n")
            f.write(f"  Mean/cell: {mutations_filtered.mean():.2f}\n")
            f.write(f"  Sparsity: {100*sparsity:.2f}%\n")
            f.write(f"  Avg muts/context: {avg_per_context:.3f}\n")
            f.write(f"  NMF Feasibility: {feasibility}\n")
        log(f"Saved statistics: {stats_file}")
    
    # Prepare results
    stats = {
        'n_cells_original': data_original.shape[1],
        'n_cells_filtered': n_retained,
        'pct_retained': pct_retained,
        'mutations_original': int(mutations_per_cell.sum()),
        'mutations_filtered': int(mutations_filtered.sum()),
        'mean_muts_per_cell': mutations_filtered.mean(),
        'median_muts_per_cell': mutations_filtered.median(),
        'sparsity': sparsity,
        'avg_muts_per_context': avg_per_context,
        'feasibility': feasibility
    }
    
    results = {
        'original_matrix': data_original,
        'filtered_matrix': data_filtered,
        'mutations_per_cell': mutations_filtered,
        'stats': stats
    }
    
    log("\n" + "="*80)
    log("FILTERING COMPLETE")
    log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 2: PCA ANALYSIS FOR PATTERN DETECTION
# =============================================================================

def run_pca_analysis(data_matrix, output_dir=None, max_components=10):
    """
    Perform PCA analysis to check for patterns in sparse data
    
    PCA allows negative values and doesn't force positivity, making it useful
    for validating whether patterns exist before NMF analysis.
    
    Parameters:
    -----------
    data_matrix : pd.DataFrame
        Mutation matrix (96 contexts × cells)
    output_dir : str, optional
        Directory to save plots
    max_components : int
        Maximum number of PCs to compute (default: 10)
    
    Returns:
    --------
    dict with:
        - 'pca_model': fitted PCA object
        - 'principal_components': transformed data
        - 'explained_variance_ratio': variance explained by each PC
        - 'cumulative_variance': cumulative variance explained
    """
    log("="*80)
    log("STEP 2: PCA ANALYSIS FOR PATTERN DETECTION")
    log("="*80)
    
    # Prepare data (transpose: cells × contexts)
    X = data_matrix.T.values
    
    log(f"Data shape: {X.shape[0]} cells × {X.shape[1]} contexts")
    log(f"Running PCA with max {max_components} components...")
    
    # Limit max_components to data dimensions
    max_comp_possible = min(max_components, X.shape[0], X.shape[1])
    
    # Fit PCA
    pca = PCA(n_components=max_comp_possible)
    principal_components = pca.fit_transform(X)
    
    log(f"PCA complete. Computed {pca.n_components_} components")
    
    # Variance explained
    var_explained = pca.explained_variance_ratio_
    cumulative_var = np.cumsum(var_explained)
    
    log(f"\nVariance explained:")
    log(f"  PC1: {100*var_explained[0]:.2f}%")
    log(f"  PC1-2: {100*cumulative_var[1]:.2f}%")
    log(f"  PC1-3: {100*cumulative_var[2]:.2f}%")
    if len(cumulative_var) >= 5:
        log(f"  PC1-5: {100*cumulative_var[4]:.2f}%")
    
    # Create visualization
    fig = plt.figure(figsize=(18, 10))
    
    # -------------------------------------------------------------------------
    # Panel 1: Scree plot (variance explained)
    # -------------------------------------------------------------------------
    ax1 = plt.subplot(2, 3, 1)
    x_pos = np.arange(1, len(var_explained) + 1)
    ax1.bar(x_pos, var_explained * 100, color='steelblue', alpha=0.8)
    ax1.set_xlabel('Principal Component', fontsize=11)
    ax1.set_ylabel('Variance Explained (%)', fontsize=11)
    ax1.set_title('Scree Plot', fontsize=12, fontweight='bold')
    ax1.set_xticks(x_pos)
    ax1.grid(axis='y', alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 2: Cumulative variance explained
    # -------------------------------------------------------------------------
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(x_pos, cumulative_var * 100, 'o-', linewidth=2, markersize=6, color='darkgreen')
    ax2.axhline(80, color='red', linestyle='--', alpha=0.5, label='80% threshold')
    ax2.axhline(90, color='orange', linestyle='--', alpha=0.5, label='90% threshold')
    ax2.set_xlabel('Number of Components', fontsize=11)
    ax2.set_ylabel('Cumulative Variance (%)', fontsize=11)
    ax2.set_title('Cumulative Variance Explained', fontsize=12, fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 3: PC1 vs PC2 scatter (cells)
    # -------------------------------------------------------------------------
    ax3 = plt.subplot(2, 3, 3)
    ax3.scatter(principal_components[:, 0], principal_components[:, 1], 
                alpha=0.5, s=20, color='purple')
    ax3.set_xlabel(f'PC1 ({100*var_explained[0]:.1f}%)', fontsize=11)
    ax3.set_ylabel(f'PC2 ({100*var_explained[1]:.1f}%)', fontsize=11)
    ax3.set_title('PC1 vs PC2 (Cell Space)', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(0, color='gray', linewidth=0.5)
    ax3.axvline(0, color='gray', linewidth=0.5)
    
    # -------------------------------------------------------------------------
    # Panel 4-6: Top 3 PC loadings (trinucleotide context contributions)
    # -------------------------------------------------------------------------
    contexts = data_matrix.index.tolist()
    mutation_types = [ctx.split('[')[1].split(']')[0] for ctx in contexts]
    
    for pc_idx in range(min(3, pca.n_components_)):
        ax = plt.subplot(2, 3, 4 + pc_idx)
        
        loadings = pca.components_[pc_idx]
        x = np.arange(len(loadings))
        colors = [COSMIC_COLORS[mut] for mut in mutation_types]
        
        # Bar plot with COSMIC colors
        bars = ax.bar(x, loadings, color=colors, width=0.8, edgecolor='none')
        
        # Add separation lines
        for i in range(1, 6):
            ax.axvline(i*16 - 0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        
        ax.axhline(0, color='black', linewidth=1)
        ax.set_ylabel('Loading', fontsize=10)
        ax.set_title(f'PC{pc_idx+1} Loadings ({100*var_explained[pc_idx]:.1f}% var)', 
                     fontsize=11, fontweight='bold')
        
        # X-axis labels
        tick_positions = [8, 24, 40, 56, 72, 88]
        tick_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=9)
        ax.set_xlim(-0.5, 95.5)
        
        # Clean up
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.suptitle(f'PCA Analysis: {X.shape[0]} Cells with {X.shape[1]} Mutation Contexts',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    # Save figure
    if output_dir:
        save_figure(fig, output_dir, 'pca_analysis.png')
    
    plt.show()
    
    # Interpretation
    log("\nPCA Interpretation:")
    if var_explained[0] > 0.2:
        log("  ✓ Strong primary pattern detected (PC1 > 20%)")
    elif var_explained[0] > 0.1:
        log("  ○ Moderate primary pattern (PC1 = 10-20%)")
    else:
        log("  ⚠ Weak patterns (PC1 < 10%) - data may be very noisy")
    
    n_comps_80pct = np.where(cumulative_var >= 0.80)[0][0] + 1 if any(cumulative_var >= 0.80) else len(cumulative_var)
    log(f"  • {n_comps_80pct} components explain ≥80% variance")
    
    results = {
        'pca_model': pca,
        'principal_components': principal_components,
        'explained_variance_ratio': var_explained,
        'cumulative_variance': cumulative_var,
        'loadings': pca.components_
    }
    
    log("\n" + "="*80)
    log("PCA ANALYSIS COMPLETE")
    log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 3: SIMILARITY METRICS FOR COSMIC MATCHING
# =============================================================================

def calculate_cosine_similarity(sig1, sig2):
    """Cosine similarity (1 - cosine distance)"""
    return 1 - cosine(sig1, sig2)


def calculate_inverse_euclidean_similarity(sig1, sig2):
    """Inverse Euclidean distance similarity"""
    dist = euclidean(sig1, sig2)
    return 1 / (1 + dist)


def calculate_dot_product_similarity(sig1, sig2):
    """Normalized dot product similarity"""
    sig1_norm = sig1 / np.linalg.norm(sig1) if np.linalg.norm(sig1) > 0 else sig1
    sig2_norm = sig2 / np.linalg.norm(sig2) if np.linalg.norm(sig2) > 0 else sig2
    return np.dot(sig1_norm, sig2_norm)


def match_to_cosmic_all_metrics(extracted_signature, cosmic_sigs, top_n=5):
    """
    Match extracted signature to COSMIC using all three metrics
    
    Returns dict with results for each metric:
        {'cosine': [...], 'inverse_euclidean': [...], 'dot_product': [...]}
    """
    results = {}
    
    # Cosine similarity
    similarities = []
    for sig_name in cosmic_sigs.columns:
        cosmic_sig = cosmic_sigs[sig_name].values
        sim = calculate_cosine_similarity(extracted_signature, cosmic_sig)
        similarities.append((sig_name, sim))
    similarities.sort(key=lambda x: x[1], reverse=True)
    results['cosine'] = similarities[:top_n]
    
    # Inverse Euclidean
    similarities = []
    for sig_name in cosmic_sigs.columns:
        cosmic_sig = cosmic_sigs[sig_name].values
        sim = calculate_inverse_euclidean_similarity(extracted_signature, cosmic_sig)
        similarities.append((sig_name, sim))
    similarities.sort(key=lambda x: x[1], reverse=True)
    results['inverse_euclidean'] = similarities[:top_n]
    
    # Dot product
    similarities = []
    for sig_name in cosmic_sigs.columns:
        cosmic_sig = cosmic_sigs[sig_name].values
        sim = calculate_dot_product_similarity(extracted_signature, cosmic_sig)
        similarities.append((sig_name, sim))
    similarities.sort(key=lambda x: x[1], reverse=True)
    results['dot_product'] = similarities[:top_n]
    
    return results


# =============================================================================
# STEP 4: ENHANCED NMF WITH HUNGARIAN MATCHING
# =============================================================================

def run_nmf_replicate(matrix, n_components, max_iter=5000, tol=1e-6, random_state=None):
    """Run single NMF replicate"""
    model = NMF(
        n_components=n_components,
        init='random',
        random_state=random_state,
        max_iter=max_iter,
        tol=tol,
        solver='cd',
        beta_loss='frobenius',
        alpha_W=0.0,
        alpha_H=0.0
    )
    
    W = model.fit_transform(matrix.T)  # Transpose: samples × features
    H = model.components_  # Components × features
    
    # Reconstruction error
    reconstruction = W @ H
    error = np.linalg.norm(matrix.T - reconstruction, 'fro')
    
    return W, H, error, model.n_iter_ < max_iter


def match_signatures_hungarian(H_matrices):
    """Match signatures across replicates using Hungarian algorithm"""
    if len(H_matrices) < 2:
        return H_matrices
    
    reference = H_matrices[0]
    matched_list = [reference]
    
    for signatures in H_matrices[1:]:
        n_sigs = reference.shape[0]
        
        # Cosine similarity matrix
        sim_matrix = np.zeros((n_sigs, n_sigs))
        for i in range(n_sigs):
            for j in range(n_sigs):
                sim_matrix[i, j] = 1 - cosine(reference[i], signatures[j])
        
        # Hungarian algorithm for optimal matching
        row_ind, col_ind = linear_sum_assignment(-sim_matrix)
        matched_sigs = signatures[col_ind]
        matched_list.append(matched_sigs)
    
    return matched_list


def refine_exposures_nnls(signatures, data):
    """Refine exposures using Non-Negative Least Squares"""
    n_sigs = signatures.shape[0]
    n_samples = data.shape[1]
    exposures = np.zeros((n_sigs, n_samples))
    
    for i in range(n_samples):
        sample_mutations = data[:, i]
        try:
            result, residual = nnls(signatures.T, sample_mutations)
            exposures[:, i] = result
        except:
            exposures[:, i] = 0
    
    return exposures


def consensus_signatures(H_matrices, W_matrices, data_matrix, use_nnls=True):
    """
    Find consensus signatures across replicates
    
    Uses Hungarian algorithm to match signatures, then averages
    Optionally refines exposures with NNLS
    """
    # Match signatures across replicates
    matched_H = match_signatures_hungarian(H_matrices)
    
    # Compute consensus signatures
    consensus_H = np.mean(matched_H, axis=0)
    
    # Refine exposures
    if use_nnls:
        consensus_W = refine_exposures_nnls(consensus_H, data_matrix).T
    else:
        # Match and average exposures
        consensus_W = np.mean(W_matrices, axis=0)
    
    # Calculate reconstruction error
    reconstruction = consensus_W @ consensus_H
    error = np.linalg.norm(data_matrix.T - reconstruction, 'fro')
    
    return consensus_H, consensus_W, error


# =============================================================================
# STEP 5: COMPREHENSIVE NMF ANALYSIS
# =============================================================================

def run_comprehensive_nmf(data_matrix, cosmic_sigs, output_dir,
                         min_signatures=1, max_signatures=10,
                         n_replicates=10, max_iter=5000, tol=1e-6,
                         use_nnls=True):
    """
    Run comprehensive NMF analysis with quality metrics
    
    Parameters:
    -----------
    data_matrix : pd.DataFrame
        Filtered mutation matrix (96 contexts × cells)
    cosmic_sigs : pd.DataFrame
        COSMIC reference signatures
    output_dir : str
        Output directory
    min_signatures : int
        Minimum number of signatures to test
    max_signatures : int
        Maximum number of signatures to test
    n_replicates : int
        Number of NMF replicates per signature count
    max_iter : int
        Maximum NMF iterations
    tol : float
        Convergence tolerance
    use_nnls : bool
        Refine exposures with NNLS
    
    Returns:
    --------
    dict with comprehensive results
    """
    log("="*80)
    log("STEP 3: COMPREHENSIVE NMF ANALYSIS")
    log("="*80)
    
    log(f"\nConfiguration:")
    log(f"  Signature range: {min_signatures} to {max_signatures}")
    log(f"  Replicates: {n_replicates}")
    log(f"  Max iterations: {max_iter}")
    log(f"  Tolerance: {tol}")
    log(f"  NNLS refinement: {use_nnls}")
    
    results = {
        'n_components': [],
        'frobenius_error': [],
        'avg_silhouette': [],
        'convergence_rate': [],
        'best_H': [],
        'best_W': [],
        'cosmic_matches': []
    }
    
    for n_comp in range(min_signatures, max_signatures + 1):
        log(f"\n{'='*70}")
        log(f"Testing {n_comp} signature(s)")
        log(f"{'='*70}")
        
        H_replicates = []
        W_replicates = []
        errors = []
        converged_count = 0
        
        # Run replicates
        for rep in range(n_replicates):
            W, H, error, converged = run_nmf_replicate(
                data_matrix.values,
                n_components=n_comp,
                max_iter=max_iter,
                tol=tol,
                random_state=rep
            )
            
            H_replicates.append(H)
            W_replicates.append(W)
            errors.append(error)
            if converged:
                converged_count += 1
            
            if (rep + 1) % 5 == 0 or rep == 0:
                log(f"  Replicate {rep + 1}/{n_replicates}: error={error:.2f}, converged={converged}")
        
        # Consensus
        log(f"  Computing consensus with Hungarian matching...")
        best_H, best_W, best_error = consensus_signatures(
            H_replicates, W_replicates, data_matrix.values, use_nnls=use_nnls
        )
        
        # Silhouette score (signature stability)
        log(f"  Calculating silhouette score...")
        if n_comp > 1:
            # Stack all signatures from all replicates
            all_sigs = np.vstack(H_replicates)
            labels = np.repeat(range(n_comp), n_replicates)
            
            try:
                silhouette = silhouette_score(all_sigs, labels, metric='cosine')
            except:
                silhouette = np.nan
        else:
            silhouette = 1.0
        
        log(f"  Silhouette score: {silhouette:.4f}")
        
        # COSMIC matching
        log(f"  Matching to COSMIC signatures...")
        cosmic_matches = []
        for i in range(n_comp):
            matches = match_to_cosmic_all_metrics(best_H[i], cosmic_sigs, top_n=3)
            cosmic_matches.append(matches)
            
            log(f"    Signature {i+1}:")
            log(f"      Cosine:     {matches['cosine'][0][0]} ({matches['cosine'][0][1]:.4f})")
            log(f"      Inv.Euclid: {matches['inverse_euclidean'][0][0]} ({matches['inverse_euclidean'][0][1]:.4f})")
            log(f"      DotProduct: {matches['dot_product'][0][0]} ({matches['dot_product'][0][1]:.4f})")
        
        # Store results
        results['n_components'].append(n_comp)
        results['frobenius_error'].append(best_error)
        results['avg_silhouette'].append(silhouette)
        results['convergence_rate'].append(100 * converged_count / n_replicates)
        results['best_H'].append(best_H)
        results['best_W'].append(best_W)
        results['cosmic_matches'].append(cosmic_matches)
        
        log(f"\n  Summary:")
        log(f"    Frobenius error: {best_error:.2f}")
        log(f"    Silhouette: {silhouette:.4f}")
        log(f"    Convergence: {100*converged_count/n_replicates:.1f}%")
    
    log("\n" + "="*80)
    log("NMF ANALYSIS COMPLETE")
    log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 6: QUALITY METRICS VISUALIZATION
# =============================================================================

def plot_quality_metrics(results, output_dir):
    """
    Plot Frobenius error and Silhouette scores
    
    These help determine optimal number of signatures
    """
    log("="*80)
    log("STEP 4: QUALITY METRICS VISUALIZATION")
    log("="*80)
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
    
    n_comps = results['n_components']
    
    # -------------------------------------------------------------------------
    # Panel 1: Frobenius reconstruction error
    # -------------------------------------------------------------------------
    ax1 = axes[0]
    ax1.plot(n_comps, results['frobenius_error'], 'o-', linewidth=2, 
             markersize=8, color='steelblue')
    ax1.set_xlabel('Number of Signatures', fontsize=12)
    ax1.set_ylabel('Frobenius Reconstruction Error', fontsize=12)
    ax1.set_title('NMF Reconstruction Error (Lower is Better)', 
                  fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(n_comps)
    
    # Find elbow (if possible)
    errors = np.array(results['frobenius_error'])
    if len(errors) > 2:
        # Simple elbow detection: max curvature
        diffs = np.diff(errors)
        second_diffs = np.diff(diffs)
        if len(second_diffs) > 0:
            elbow_idx = np.argmax(np.abs(second_diffs)) + 1
            elbow_n = n_comps[elbow_idx]
            ax1.axvline(elbow_n, color='red', linestyle='--', alpha=0.6,
                       label=f'Potential elbow: {elbow_n}')
            ax1.legend()
    
    # -------------------------------------------------------------------------
    # Panel 2: Average Silhouette score
    # -------------------------------------------------------------------------
    ax2 = axes[1]
    ax2.plot(n_comps, results['avg_silhouette'], 'o-', linewidth=2,
             markersize=8, color='darkgreen')
    ax2.axhline(0.5, color='orange', linestyle='--', alpha=0.5, 
                label='Good threshold (0.5)')
    ax2.axhline(0.7, color='red', linestyle='--', alpha=0.5,
                label='Excellent threshold (0.7)')
    ax2.set_xlabel('Number of Signatures', fontsize=12)
    ax2.set_ylabel('Average Silhouette Score', fontsize=12)
    ax2.set_title('Signature Stability (Higher is Better)',
                  fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(n_comps)
    ax2.set_ylim([-1, 1])
    
    plt.tight_layout()
    save_figure(fig, output_dir, 'quality_metrics.png')
    plt.show()
    
    # Recommendation
    log("\nQuality Metric Recommendations:")
    
    # Best by Frobenius
    best_frob_idx = np.argmin(results['frobenius_error'])
    best_frob_n = n_comps[best_frob_idx]
    log(f"  Lowest Frobenius error: {best_frob_n} signatures (error={results['frobenius_error'][best_frob_idx]:.2f})")
    
    # Best by Silhouette (excluding n=1)
    silhouettes = np.array(results['avg_silhouette'])
    if len(silhouettes) > 1:
        best_sil_idx = np.argmax(silhouettes[1:]) + 1  # Skip first
        best_sil_n = n_comps[best_sil_idx]
        log(f"  Highest Silhouette (n>1): {best_sil_n} signatures (score={results['avg_silhouette'][best_sil_idx]:.4f})")
    
    log("\n" + "="*80)
    log("QUALITY METRICS COMPLETE")
    log("="*80 + "\n")


# =============================================================================
# STEP 7: VISUALIZATION OF EXTRACTED SIGNATURES
# =============================================================================

def plot_cosmic_signature(signature_vector, title, ax=None, normalize=True):
    """Plot signature in COSMIC style"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 4))
    
    # Normalize to percentages
    if normalize:
        total = signature_vector.sum()
        values = (signature_vector / total * 100) if total > 0 else signature_vector
    else:
        values = signature_vector
    
    # Get mutation types for coloring
    # Assumes 96 trinucleotide contexts in standard order
    mutation_types = []
    for i in range(6):
        mutation_types.extend([list(COSMIC_COLORS.keys())[i]] * 16)
    
    # Create bars
    x = np.arange(len(values))
    colors = [COSMIC_COLORS[mut] for mut in mutation_types]
    ax.bar(x, values, color=colors, width=0.8, edgecolor='none')
    
    # Separation lines
    for i in range(1, 6):
        ax.axvline(i*16 - 0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    # Labels
    ax.set_ylabel('% of mutations', fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # X-axis
    tick_positions = [8, 24, 40, 56, 72, 88]
    tick_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set_xlim(-0.5, 95.5)
    
    # Clean up
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return ax


def visualize_extracted_signatures(results, cosmic_sigs, output_dir, 
                                   signature_count=None):
    """
    Visualize extracted signatures with COSMIC comparisons
    
    If signature_count is None, visualizes the model with lowest Frobenius error
    """
    log("="*80)
    log("STEP 5: SIGNATURE VISUALIZATION")
    log("="*80)
    
    # Select which result to visualize
    if signature_count is None:
        # Use model with lowest Frobenius error
        best_idx = np.argmin(results['frobenius_error'])
        signature_count = results['n_components'][best_idx]
        log(f"Visualizing best model: {signature_count} signatures (lowest Frobenius error)")
    else:
        best_idx = results['n_components'].index(signature_count)
        log(f"Visualizing requested model: {signature_count} signatures")
    
    best_H = results['best_H'][best_idx]
    cosmic_matches = results['cosmic_matches'][best_idx]
    
    # Create figure with extracted + COSMIC comparisons
    n_rows = signature_count * 2  # Extracted + best COSMIC match for each
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 4 * n_rows))
    
    if signature_count == 1:
        axes = [axes] if n_rows == 1 else axes
    
    for i in range(signature_count):
        # Get best match (prioritize inverse Euclidean, which handles sparse data well)
        best_match_name = cosmic_matches[i]['inverse_euclidean'][0][0]
        best_match_score = cosmic_matches[i]['inverse_euclidean'][0][1]
        
        # Also get cosine and dot product
        cosine_match = cosmic_matches[i]['cosine'][0]
        dotprod_match = cosmic_matches[i]['dot_product'][0]
        
        # Plot extracted signature
        title_extracted = (f"Extracted Signature {i+1}\n"
                          f"Inv.Euc: {best_match_name} ({best_match_score:.3f}) | "
                          f"Cosine: {cosine_match[0]} ({cosine_match[1]:.3f}) | "
                          f"DotProd: {dotprod_match[0]} ({dotprod_match[1]:.3f})")
        
        plot_cosmic_signature(best_H[i], title_extracted, ax=axes[i*2])
        
        # Plot corresponding COSMIC signature
        cosmic_sig = cosmic_sigs[best_match_name].values
        title_cosmic = f"COSMIC {best_match_name} (Reference)"
        plot_cosmic_signature(cosmic_sig, title_cosmic, ax=axes[i*2 + 1])
    
    plt.tight_layout()
    save_figure(fig, output_dir, f'extracted_signatures_{signature_count}comp.png')
    plt.show()
    
    log("\n" + "="*80)
    log("SIGNATURE VISUALIZATION COMPLETE")
    log("="*80 + "\n")


# =============================================================================
# STEP 8: SAVE RESULTS
# =============================================================================

def save_nmf_results(results, data_matrix, output_dir, mutation_threshold):
    """Save all NMF results to files"""
    log("="*80)
    log("STEP 6: SAVING RESULTS")
    log("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Summary table
    summary_df = pd.DataFrame({
        'n_signatures': results['n_components'],
        'frobenius_error': results['frobenius_error'],
        'avg_silhouette': results['avg_silhouette'],
        'convergence_pct': results['convergence_rate']
    })
    
    summary_file = output_path / "nmf_summary.txt"
    summary_df.to_csv(summary_file, sep='\t', index=False, float_format='%.4f')
    log(f"Saved summary: {summary_file}")
    
    # Save best model (lowest Frobenius error)
    best_idx = np.argmin(results['frobenius_error'])
    best_n = results['n_components'][best_idx]
    best_H = results['best_H'][best_idx]
    best_W = results['best_W'][best_idx]
    
    # Signatures
    sig_df = pd.DataFrame(
        best_H.T,
        index=data_matrix.index,
        columns=[f"Signature_{i+1}" for i in range(best_n)]
    )
    sig_file = output_path / f"signatures_{best_n}comp.txt"
    sig_df.to_csv(sig_file, sep='\t', float_format='%.6f')
    log(f"Saved signatures: {sig_file}")
    
    # Exposures
    exp_df = pd.DataFrame(
        best_W,
        index=data_matrix.columns,
        columns=[f"Signature_{i+1}" for i in range(best_n)]
    )
    exp_file = output_path / f"exposures_{best_n}comp.txt"
    exp_df.to_csv(exp_file, sep='\t', float_format='%.6f')
    log(f"Saved exposures: {exp_file}")
    
    # COSMIC matches
    with open(output_path / "cosmic_matches_all_metrics.txt", 'w') as f:
        f.write("="*80 + "\n")
        f.write("COSMIC Signature Matches - All Metrics\n")
        f.write(f"Mutation threshold: ≥{mutation_threshold}\n")
        f.write("="*80 + "\n\n")
        
        for n_idx, n_comp in enumerate(results['n_components']):
            f.write(f"\n{n_comp} Signature(s):\n")
            f.write("-"*60 + "\n")
            
            for i, matches in enumerate(results['cosmic_matches'][n_idx]):
                f.write(f"\n  Signature {i+1}:\n")
                
                for metric in ['cosine', 'inverse_euclidean', 'dot_product']:
                    f.write(f"    {metric.upper().replace('_', ' ')}:\n")
                    for sig_name, similarity in matches[metric]:
                        f.write(f"      {sig_name:15s} {similarity:.4f}\n")
    
    log(f"Saved COSMIC matches: {output_path / 'cosmic_matches_all_metrics.txt'}")
    
    log("\n" + "="*80)
    log("RESULTS SAVED")
    log("="*80 + "\n")


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_complete_pipeline(matrix_file, cosmic_file, output_dir,
                          mutation_threshold=20,
                          min_signatures=1, max_signatures=10,
                          n_replicates=10, max_iter=5000, tol=1e-6,
                          use_nnls=True, run_pca=True):
    """
    Complete single-cell mutation signature extraction pipeline
    
    Parameters:
    -----------
    matrix_file : str
        Path to mutation matrix (96 contexts × cells)
    cosmic_file : str
        Path to COSMIC reference signatures
    output_dir : str
        Output directory for all results
    mutation_threshold : int
        Minimum mutations per cell
    min_signatures : int
        Minimum signatures to test
    max_signatures : int
        Maximum signatures to test
    n_replicates : int
        NMF replicates
    max_iter : int
        Max NMF iterations
    tol : float
        Convergence tolerance
    use_nnls : bool
        Refine exposures with NNLS
    run_pca : bool
        Run PCA preprocessing analysis
    
    Returns:
    --------
    dict with all analysis results
    """
    
    print("\n" + "="*80)
    print("COMPREHENSIVE SINGLE-CELL MUTATION SIGNATURE EXTRACTION PIPELINE")
    print("="*80)
    print(f"\nInput: {matrix_file}")
    print(f"COSMIC: {cosmic_file}")
    print(f"Output: {output_dir}")
    print(f"\nMutation threshold: ≥{mutation_threshold}")
    print(f"Signature range: {min_signatures}-{max_signatures}")
    print(f"Replicates: {n_replicates}")
    print("="*80 + "\n")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # -------------------------------------------------------------------------
    # STEP 1: Load and filter data
    # -------------------------------------------------------------------------
    filtering_results = load_and_filter_mutation_matrix(
        matrix_file, 
        mutation_threshold=mutation_threshold,
        output_dir=output_dir,
        save_filtered=True
    )
    
    filtered_matrix = filtering_results['filtered_matrix']
    
    # -------------------------------------------------------------------------
    # STEP 2: PCA analysis (optional)
    # -------------------------------------------------------------------------
    pca_results = None
    if run_pca:
        pca_results = run_pca_analysis(
            filtered_matrix,
            output_dir=output_dir,
            max_components=10
        )
    
    # -------------------------------------------------------------------------
    # STEP 3: Load COSMIC signatures
    # -------------------------------------------------------------------------
    log("="*80)
    log("LOADING COSMIC REFERENCE SIGNATURES")
    log("="*80)
    
    cosmic_sigs = pd.read_csv(cosmic_file, sep='\t', index_col=0)
    log(f"Loaded {cosmic_sigs.shape[1]} COSMIC signatures")
    
    # Ensure contexts match
    if not all(filtered_matrix.index == cosmic_sigs.index):
        log("Reordering COSMIC signatures to match matrix contexts...")
        cosmic_sigs = cosmic_sigs.loc[filtered_matrix.index]
    
    log("Contexts verified and aligned.\n")
    
    # -------------------------------------------------------------------------
    # STEP 4: Run comprehensive NMF
    # -------------------------------------------------------------------------
    nmf_results = run_comprehensive_nmf(
        filtered_matrix,
        cosmic_sigs,
        output_dir,
        min_signatures=min_signatures,
        max_signatures=max_signatures,
        n_replicates=n_replicates,
        max_iter=max_iter,
        tol=tol,
        use_nnls=use_nnls
    )
    
    # -------------------------------------------------------------------------
    # STEP 5: Plot quality metrics
    # -------------------------------------------------------------------------
    plot_quality_metrics(nmf_results, output_dir)
    
    # -------------------------------------------------------------------------
    # STEP 6: Visualize signatures
    # -------------------------------------------------------------------------
    visualize_extracted_signatures(
        nmf_results,
        cosmic_sigs,
        output_dir,
        signature_count=None  # Use best model
    )
    
    # -------------------------------------------------------------------------
    # STEP 7: Save all results
    # -------------------------------------------------------------------------
    save_nmf_results(
        nmf_results,
        filtered_matrix,
        output_dir,
        mutation_threshold
    )
    
    # -------------------------------------------------------------------------
    # Final summary
    # -------------------------------------------------------------------------
    print("\n" + "="*80)
    print("PIPELINE COMPLETE!")
    print("="*80)
    
    best_idx = np.argmin(nmf_results['frobenius_error'])
    best_n = nmf_results['n_components'][best_idx]
    
    print(f"\nBest model: {best_n} signatures")
    print(f"  Frobenius error: {nmf_results['frobenius_error'][best_idx]:.2f}")
    print(f"  Silhouette score: {nmf_results['avg_silhouette'][best_idx]:.4f}")
    
    print(f"\nTop signature matches:")
    for i, matches in enumerate(nmf_results['cosmic_matches'][best_idx]):
        print(f"  Signature {i+1}:")
        print(f"    Cosine:       {matches['cosine'][0][0]} ({matches['cosine'][0][1]:.4f})")
        print(f"    Inv.Euclidean: {matches['inverse_euclidean'][0][0]} ({matches['inverse_euclidean'][0][1]:.4f})")
        print(f"    Dot Product:   {matches['dot_product'][0][0]} ({matches['dot_product'][0][1]:.4f})")
    
    print(f"\nAll results saved to: {output_dir}")
    print("="*80 + "\n")
    
    # Return everything
    return {
        'filtering': filtering_results,
        'pca': pca_results,
        'nmf': nmf_results,
        'filtered_matrix': filtered_matrix,
        'cosmic_sigs': cosmic_sigs
    }

#%%
# =============================================================================
# EXAMPLE USAGE (FOR INTERACTIVE EXECUTION)
# =============================================================================

if __name__ == "__main__":
    
    # =========================================================================
    # CONFIGURATION - UPDATE THESE PATHS
    # =========================================================================
    
    MATRIX_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/ALL_Basal_Cell_SNP_matrix_for_SigProfiler.txt"
    COSMIC_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/COSMIC_v3.4_SBS_GRCh38.txt"
    OUTPUT_DIR = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/comprehensive_pipeline_v1"
    
    MUTATION_THRESHOLD = 30    # Minimum mutations per cell
    MIN_SIGNATURES = 1         # Start testing from 1 signature
    MAX_SIGNATURES = 8         # Test up to 8 signatures
    N_REPLICATES = 10          # NMF replicates per signature count
    MAX_ITER = 5000            # Maximum NMF iterations
    TOL = 1e-6                 # Convergence tolerance
    USE_NNLS = True            # Refine exposures with NNLS
    RUN_PCA = True             # Run PCA preprocessing
    
    # =========================================================================
    # RUN COMPLETE PIPELINE
    # =========================================================================
    
    results = run_complete_pipeline(
        matrix_file=MATRIX_FILE,
        cosmic_file=COSMIC_FILE,
        output_dir=OUTPUT_DIR,
        mutation_threshold=MUTATION_THRESHOLD,
        min_signatures=MIN_SIGNATURES,
        max_signatures=MAX_SIGNATURES,
        n_replicates=N_REPLICATES,
        max_iter=MAX_ITER,
        tol=TOL,
        use_nnls=USE_NNLS,
        run_pca=RUN_PCA
    )
    
    print("\n✓ Pipeline execution complete!")
    print(f"Check results in: {OUTPUT_DIR}")
#%%
###
### Jump to 3924
###
#!/usr/bin/env python
"""
Semi-Supervised Signature Refitting for Single-Cell Data

Instead of de novo signature extraction, this approach:
1. Uses KNOWN COSMIC signatures relevant to HNSCC
2. Solves for per-cell weights (exposures) using NNLS
3. Evaluates how well H × W approximates the original mutation matrix X

This version includes:
- Scree plot-based signature selection (elbow detection)
- Fixed SBS40 signature detection
- Proper evaluation using Eckart-Young theorem

Author: Jake Lehle
Date: 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
from scipy.optimize import nnls
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, r2_score
import warnings
warnings.filterwarnings('ignore')

# COSMIC signature colors
COSMIC_COLORS = {
    'C>A': '#1EBFF0',
    'C>G': '#050708',
    'C>T': '#E62725',
    'T>A': '#CBCACB',
    'T>C': '#A1CE63',
    'T>G': '#EDB6C2'
}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def log(msg, level="INFO"):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {level}: {msg}", flush=True)


def save_figure(fig, output_dir, filename, dpi=300):
    """Save figure to output directory"""
    filepath = Path(output_dir) / filename
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
    log(f"Saved: {filepath}")
    return filepath


# =============================================================================
# STEP 1: EXTRACT HNSCC-RELEVANT COSMIC SIGNATURES
# =============================================================================

def extract_hnscc_signatures(cosmic_file, output_dir=None):
    """
    Extract HNSCC-relevant COSMIC signatures
    
    HNSCC signatures based on literature:
    SBS1, SBS2, SBS3, SBS4, SBS5, SBS7a, SBS7b, SBS7d, SBS13, 
    SBS16, SBS17a, SBS17b, SBS18, SBS33, SBS40
    
    NOTE: SBS40 may appear as 'SBS40a', 'SBS40b', 'SBS40c' in COSMIC v3.4
    
    Parameters:
    -----------
    cosmic_file : str
        Path to COSMIC_v3.4_SBS_GRCh38.txt
    output_dir : str, optional
        Directory to save extracted signatures
    
    Returns:
    --------
    pd.DataFrame : HNSCC signature matrix (96 contexts × signatures)
    """
    log("="*80)
    log("EXTRACTING HNSCC-RELEVANT COSMIC SIGNATURES")
    log("="*80)
    
    # HNSCC-relevant signatures
    hnscc_signatures = [
        'SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5',
        'SBS7a', 'SBS7b', 'SBS7d', 'SBS13',
        'SBS16', 'SBS17a', 'SBS17b', 'SBS18',
        'SBS33', 'SBS40'
    ]
    
    log(f"\nLoading COSMIC signatures from: {cosmic_file}")
    cosmic_all = pd.read_csv(cosmic_file, sep='\t', index_col=0)
    log(f"Total COSMIC signatures available: {cosmic_all.shape[1]}")
    
    # Check for SBS40 variants (SBS40a, SBS40b, SBS40c)
    available_hnscc = []
    missing_hnscc = []
    
    for sig in hnscc_signatures:
        if sig in cosmic_all.columns:
            available_hnscc.append(sig)
        elif sig == 'SBS40':
            # Check for SBS40 variants
            variants = [col for col in cosmic_all.columns if col.startswith('SBS40')]
            if variants:
                log(f"\nℹ️  SBS40 variants found: {', '.join(variants)}", level="INFO")
                # Use SBS40a by default (most common)
                if 'SBS40a' in variants:
                    available_hnscc.append('SBS40a')
                    log(f"Using SBS40a as SBS40 representative")
                else:
                    available_hnscc.append(variants[0])
                    log(f"Using {variants[0]} as SBS40 representative")
            else:
                missing_hnscc.append(sig)
        else:
            missing_hnscc.append(sig)
    
    log(f"\nHNSCC signatures found: {len(available_hnscc)}/{len(hnscc_signatures)}")
    
    if available_hnscc:
        log(f"Available: {', '.join(available_hnscc)}")
    
    if missing_hnscc:
        log(f"⚠️  Missing: {', '.join(missing_hnscc)}", level="WARNING")
    
    # Extract HNSCC signatures
    hnscc_sigs = cosmic_all[available_hnscc].copy()
    
    log(f"\nExtracted signature matrix: {hnscc_sigs.shape[0]} contexts × {hnscc_sigs.shape[1]} signatures")
    
    # Verify normalization (should sum to 1)
    col_sums = hnscc_sigs.sum(axis=0)
    if not np.allclose(col_sums, 1.0, atol=1e-5):
        log("⚠️  Signatures not normalized, normalizing now...", level="WARNING")
        hnscc_sigs = hnscc_sigs / col_sums
    
    log("✓ Signatures are properly normalized (sum to 1.0)")
    
    # Save if output directory provided
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        hnscc_file = output_path / "hnscc_cosmic_signatures.txt"
        hnscc_sigs.to_csv(hnscc_file, sep='\t', float_format='%.6f')
        log(f"\nSaved HNSCC signatures: {hnscc_file}")
        
        # Save list of signatures
        sig_list_file = output_path / "hnscc_signature_list.txt"
        with open(sig_list_file, 'w') as f:
            f.write("HNSCC-Relevant COSMIC Signatures\n")
            f.write("="*50 + "\n\n")
            for i, sig in enumerate(available_hnscc, 1):
                f.write(f"{i:2d}. {sig}\n")
        log(f"Saved signature list: {sig_list_file}")
    
    log("\n" + "="*80)
    log("HNSCC SIGNATURE EXTRACTION COMPLETE")
    log("="*80 + "\n")
    
    return hnscc_sigs


# =============================================================================
# STEP 2: SIGNATURE SELECTION VIA SCREE PLOT ELBOW DETECTION
# =============================================================================

def select_signatures_via_scree_plot(mutation_matrix, hnscc_sigs, core_signatures,
                                     candidate_pool, output_dir, max_signatures=15,
                                     verbose=True):
    """
    Select optimal number of signatures using scree plot elbow detection
    
    Strategy:
    1. Start with core signatures (e.g., SBS2, SBS3, SBS5)
    2. For each additional signature, calculate explained variance
    3. Find elbow point where adding more signatures gives diminishing returns
    4. Select best candidates based on individual explanatory power
    
    The key insight: More signatures ALWAYS reduce error, but we want to stop
    when the improvement becomes marginal (elbow point).
    
    Parameters:
    -----------
    mutation_matrix : pd.DataFrame
        Original mutation matrix X
    hnscc_sigs : pd.DataFrame
        All available HNSCC signatures
    core_signatures : list
        Core signatures to always include
    candidate_pool : list
        Pool of candidate signatures to consider
    output_dir : str
        Directory to save results
    max_signatures : int
        Maximum signatures to test (default: 15)
    verbose : bool
        Print detailed progress
    
    Returns:
    --------
    dict with selected signatures and analysis results
    """
    if verbose:
        log("="*80)
        log("SIGNATURE SELECTION VIA SCREE PLOT ELBOW DETECTION")
        log("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Verify signatures available
    available_sigs = set(hnscc_sigs.columns)
    core_available = [sig for sig in core_signatures if sig in available_sigs]
    candidates_available = [sig for sig in candidate_pool if sig in available_sigs and sig not in core_available]
    
    if verbose:
        log(f"\nSignature pool:")
        log(f"  Core signatures: {core_available}")
        log(f"  Candidate pool: {len(candidates_available)} signatures")
        log(f"  Max signatures to test: {max_signatures}")
    
    # -------------------------------------------------------------------------
    # Step 1: Rank candidates by individual explanatory power
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 1: Ranking candidates by individual explanatory power")
        log(f"{'='*70}")
    
    # For each candidate, fit it alone and see how much variance it explains
    candidate_scores = []
    
    X = mutation_matrix.values
    X_norm = np.linalg.norm(X, 'fro')
    
    for candidate_sig in candidates_available:
        # Fit just this signature
        sig_matrix = hnscc_sigs[[candidate_sig]]
        fitting = fit_signatures_nnls(mutation_matrix, sig_matrix, verbose=False)
        
        # Calculate explained variance
        residual = X - fitting['reconstruction'].values
        residual_norm = np.linalg.norm(residual, 'fro')
        explained_var = 1 - (residual_norm / X_norm)**2
        
        candidate_scores.append({
            'signature': candidate_sig,
            'explained_variance': explained_var,
            'residual_norm': residual_norm
        })
    
    # Sort by explained variance
    candidate_scores = sorted(candidate_scores, key=lambda x: x['explained_variance'], reverse=True)
    
    if verbose:
        log(f"\nCandidate ranking (by individual explained variance):")
        for i, score in enumerate(candidate_scores[:10], 1):
            log(f"  {i:2d}. {score['signature']:<10} explains {100*score['explained_variance']:.2f}% variance")
    
    # Create ordered candidate list
    candidates_ordered = [s['signature'] for s in candidate_scores]
    
    # -------------------------------------------------------------------------
    # Step 2: Build scree plot - test increasing signature counts
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 2: Building scree plot (testing 1 to {max_signatures} signatures)")
        log(f"{'='*70}")
    
    scree_data = []
    
    # Always start with core
    current_signatures = core_available.copy()
    
    # Add candidates one by one in order of explanatory power
    all_signatures_ordered = core_available + candidates_ordered
    
    for n_sigs in range(len(core_available), min(max_signatures + 1, len(all_signatures_ordered) + 1)):
        # Get signatures to test
        test_signatures = all_signatures_ordered[:n_sigs]
        sig_matrix = hnscc_sigs[test_signatures]
        
        # Fit and evaluate
        fitting = fit_signatures_nnls(mutation_matrix, sig_matrix, verbose=False)
        evaluation = evaluate_reconstruction(mutation_matrix, fitting['reconstruction'], verbose=False)
        
        # Calculate explained variance
        explained_var = 1 - evaluation['relative_frobenius_error']**2
        
        scree_data.append({
            'n_signatures': n_sigs,
            'signatures': test_signatures.copy(),
            'frobenius_error': evaluation['frobenius_error'],
            'relative_error': evaluation['relative_frobenius_error'],
            'explained_variance': explained_var,
            'optimality_ratio': evaluation['optimality_ratio']
        })
        
        if verbose and n_sigs <= 10:
            log(f"  {n_sigs:2d} signatures: Error={evaluation['frobenius_error']:>8.2f}, "
                f"RelErr={100*evaluation['relative_frobenius_error']:>5.2f}%, "
                f"ExpVar={100*explained_var:>5.2f}%")
    
    # -------------------------------------------------------------------------
    # Step 3: Find elbow point using second derivative
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 3: Detecting elbow point in scree plot")
        log(f"{'='*70}")
    
    errors = np.array([s['frobenius_error'] for s in scree_data])
    n_sigs_array = np.array([s['n_signatures'] for s in scree_data])
    
    # Calculate first derivative (rate of error decrease)
    dy = np.gradient(errors)
    
    # Calculate second derivative (rate of change of decrease)
    d2y = np.gradient(dy)
    
    # Find elbow: where second derivative is most positive (curvature changes most)
    # We want the point where error reduction starts slowing down significantly
    elbow_idx = np.argmax(d2y)
    elbow_n_sigs = scree_data[elbow_idx]['n_signatures']
    
    if verbose:
        log(f"\nElbow detection:")
        log(f"  Elbow found at: {elbow_n_sigs} signatures")
        log(f"  Error at elbow: {scree_data[elbow_idx]['frobenius_error']:.2f}")
        log(f"  Relative error: {100*scree_data[elbow_idx]['relative_error']:.2f}%")
        log(f"  Explained variance: {100*scree_data[elbow_idx]['explained_variance']:.2f}%")
    
    # -------------------------------------------------------------------------
    # Step 4: Validate elbow with marginal improvement analysis
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 4: Validating elbow point (marginal improvement analysis)")
        log(f"{'='*70}")
    
    # Calculate marginal improvement for each step
    marginal_improvements = []
    for i in range(1, len(scree_data)):
        prev_error = scree_data[i-1]['frobenius_error']
        curr_error = scree_data[i]['frobenius_error']
        improvement = 100 * (prev_error - curr_error) / prev_error
        marginal_improvements.append(improvement)
    
    # Find where improvement drops below threshold (e.g., < 1%)
    threshold_pct = 1.0
    threshold_n_sigs = len(scree_data)  # default to all
    for i, imp in enumerate(marginal_improvements):
        if imp < threshold_pct:
            threshold_n_sigs = scree_data[i+1]['n_signatures']
            if verbose:
                log(f"\n  Marginal improvement drops below {threshold_pct}% at {threshold_n_sigs} signatures")
            break
    
    # Use the more conservative of elbow and threshold methods
    final_n_sigs = min(elbow_n_sigs, threshold_n_sigs)
    
    if verbose:
        log(f"\n  Final selection: {final_n_sigs} signatures")
        log(f"    (Elbow method: {elbow_n_sigs}, Threshold method: {threshold_n_sigs})")
    
    selected_signatures = scree_data[final_n_sigs - len(core_available)]['signatures']
    
    # -------------------------------------------------------------------------
    # Step 5: Create scree plot visualization
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"Creating scree plot visualization")
        log(f"{'='*70}")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Frobenius error vs. number of signatures
    ax1 = axes[0, 0]
    ax1.plot(n_sigs_array, errors, 'o-', linewidth=2, markersize=6, color='steelblue')
    ax1.axvline(final_n_sigs, color='red', linestyle='--', linewidth=2, 
                label=f'Selected: {final_n_sigs} sigs')
    ax1.set_xlabel('Number of Signatures', fontsize=11)
    ax1.set_ylabel('Frobenius Error', fontsize=11)
    ax1.set_title('Scree Plot: Reconstruction Error', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Explained variance vs. number of signatures
    ax2 = axes[0, 1]
    explained_vars = np.array([s['explained_variance'] for s in scree_data])
    ax2.plot(n_sigs_array, explained_vars * 100, 'o-', linewidth=2, markersize=6, color='darkgreen')
    ax2.axvline(final_n_sigs, color='red', linestyle='--', linewidth=2,
                label=f'Selected: {final_n_sigs} sigs')
    ax2.axhline(90, color='orange', linestyle=':', alpha=0.5, label='90%')
    ax2.set_xlabel('Number of Signatures', fontsize=11)
    ax2.set_ylabel('Explained Variance (%)', fontsize=11)
    ax2.set_title('Explained Variance', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 100])
    
    # Panel 3: Marginal improvement per signature added
    ax3 = axes[1, 0]
    ax3.bar(n_sigs_array[1:], marginal_improvements, color='coral', alpha=0.7, edgecolor='black')
    ax3.axhline(threshold_pct, color='red', linestyle='--', linewidth=2,
                label=f'{threshold_pct}% threshold')
    ax3.axvline(final_n_sigs, color='red', linestyle='--', linewidth=2, alpha=0.5)
    ax3.set_xlabel('Number of Signatures', fontsize=11)
    ax3.set_ylabel('Marginal Improvement (%)', fontsize=11)
    ax3.set_title('Marginal Error Reduction per Signature', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    # Panel 4: Second derivative (curvature) for elbow detection
    ax4 = axes[1, 1]
    ax4.plot(n_sigs_array, d2y, 'o-', linewidth=2, markersize=6, color='purple')
    ax4.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax4.axvline(elbow_n_sigs, color='red', linestyle='--', linewidth=2,
                label=f'Elbow: {elbow_n_sigs} sigs')
    ax4.scatter(elbow_n_sigs, d2y[elbow_idx], color='red', s=150, zorder=5,
                edgecolor='black', linewidth=2)
    ax4.set_xlabel('Number of Signatures', fontsize=11)
    ax4.set_ylabel('Second Derivative (Curvature)', fontsize=11)
    ax4.set_title('Elbow Detection (Maximum Curvature)', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle('Signature Selection: Scree Plot Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    scree_plot_file = output_path / "signature_selection_scree_plot.png"
    plt.savefig(scree_plot_file, dpi=300, bbox_inches='tight')
    if verbose:
        log(f"Saved scree plot: {scree_plot_file}")
    plt.show()
    
    # -------------------------------------------------------------------------
    # Step 6: Save detailed results
    # -------------------------------------------------------------------------
    results_file = output_path / "signature_selection_scree_analysis.txt"
    
    with open(results_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SIGNATURE SELECTION: SCREE PLOT ANALYSIS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Selection Method: Elbow detection in scree plot\n")
        f.write("Principle: Stop when adding more signatures gives diminishing returns\n\n")
        
        f.write("="*80 + "\n")
        f.write("CANDIDATE RANKING (by individual explanatory power)\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"{'Rank':<6} {'Signature':<12} {'Explained Variance':<20}\n")
        f.write("-"*40 + "\n")
        for i, score in enumerate(candidate_scores, 1):
            f.write(f"{i:<6} {score['signature']:<12} {100*score['explained_variance']:>18.2f}%\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("SCREE PLOT DATA\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"{'N':<4} {'Frobenius Error':<18} {'Rel.Error%':<12} {'Exp.Var%':<12} {'Marginal Imp%':<15}\n")
        f.write("-"*80 + "\n")
        
        for i, data in enumerate(scree_data):
            marg_imp = marginal_improvements[i-1] if i > 0 else 0
            marker = " ← SELECTED" if data['n_signatures'] == final_n_sigs else ""
            f.write(f"{data['n_signatures']:<4} {data['frobenius_error']:<18.2f} "
                   f"{100*data['relative_error']:<11.2f}% {100*data['explained_variance']:<11.2f}% "
                   f"{marg_imp:<14.2f}%{marker}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("ELBOW DETECTION RESULTS\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Elbow method: {elbow_n_sigs} signatures\n")
        f.write(f"Threshold method (<{threshold_pct}% improvement): {threshold_n_sigs} signatures\n")
        f.write(f"Final selection: {final_n_sigs} signatures (more conservative)\n\n")
        
        f.write("="*80 + "\n")
        f.write("SELECTED SIGNATURES\n")
        f.write("="*80 + "\n\n")
        
        for i, sig in enumerate(selected_signatures, 1):
            marker = "[CORE]" if sig in core_available else "[ADDED]"
            f.write(f"{i:2d}. {sig:<10} {marker}\n")
        
        f.write(f"\nTotal: {len(selected_signatures)} signatures\n")
        f.write(f"Core: {len([s for s in selected_signatures if s in core_available])}\n")
        f.write(f"Added: {len([s for s in selected_signatures if s not in core_available])}\n")
    
    if verbose:
        log(f"Saved analysis: {results_file}")
        
        log(f"\n{'='*70}")
        log(f"SIGNATURE SELECTION COMPLETE")
        log(f"{'='*70}")
        log(f"\nSelected {final_n_sigs} signatures:")
        for i, sig in enumerate(selected_signatures, 1):
            marker = "🔹" if sig in core_available else "✓"
            log(f"  {i:2d}. {marker} {sig}")
        log("")
    
    results = {
        'selected_signatures': selected_signatures,
        'signature_matrix': hnscc_sigs[selected_signatures],
        'n_signatures': final_n_sigs,
        'scree_data': scree_data,
        'candidate_ranking': candidate_scores,
        'elbow_n_sigs': elbow_n_sigs,
        'final_error': scree_data[final_n_sigs - len(core_available)]['frobenius_error'],
        'final_relative_error': scree_data[final_n_sigs - len(core_available)]['relative_error'],
        'final_explained_variance': scree_data[final_n_sigs - len(core_available)]['explained_variance']
    }
    
    return results


# =============================================================================
# STEP 3: SOLVE FOR WEIGHTS USING NNLS
# =============================================================================

def fit_signatures_nnls(mutation_matrix, signature_matrix, verbose=True):
    """
    Fit known signatures to mutation data using Non-Negative Least Squares
    
    Solves: X ≈ H × W
    
    For each cell (column in X), solves:
        min ||X_cell - H × W_cell||²  subject to W_cell ≥ 0
    
    Parameters:
    -----------
    mutation_matrix : pd.DataFrame
        Mutation count matrix (96 contexts × n_cells)
    signature_matrix : pd.DataFrame
        Known signature matrix (96 contexts × k_signatures)
    verbose : bool
        Print progress updates
    
    Returns:
    --------
    dict with:
        - 'weights': pd.DataFrame (k_signatures × n_cells)
        - 'residuals': np.array (residual per cell)
        - 'reconstruction': pd.DataFrame (96 contexts × n_cells)
    """
    if verbose:
        log("="*80)
        log("FITTING SIGNATURES USING NON-NEGATIVE LEAST SQUARES")
        log("="*80)
    
    # Ensure contexts match
    if not all(mutation_matrix.index == signature_matrix.index):
        log("⚠️  Reordering signatures to match mutation matrix contexts...", level="WARNING")
        signature_matrix = signature_matrix.loc[mutation_matrix.index]
    
    n_contexts = mutation_matrix.shape[0]
    n_cells = mutation_matrix.shape[1]
    n_sigs = signature_matrix.shape[1]
    
    if verbose:
        log(f"\nProblem dimensions:")
        log(f"  Mutation contexts: {n_contexts}")
        log(f"  Cells to fit: {n_cells:,}")
        log(f"  Signatures: {n_sigs}")
        log(f"\nSolving: X ({n_contexts}×{n_cells}) ≈ H ({n_contexts}×{n_sigs}) × W ({n_sigs}×{n_cells})")
    
    # Convert to numpy arrays
    X = mutation_matrix.values  # 96 × n_cells
    H = signature_matrix.values  # 96 × n_sigs
    
    # Initialize weight matrix
    W = np.zeros((n_sigs, n_cells))
    residuals = np.zeros(n_cells)
    
    if verbose:
        log(f"\nFitting signatures to {n_cells:,} cells...")
    
    # Solve for each cell
    progress_interval = max(1, n_cells // 20)  # Update every 5%
    
    for i in range(n_cells):
        # NNLS for this cell
        cell_mutations = X[:, i]
        
        try:
            weights, residual = nnls(H, cell_mutations)
            W[:, i] = weights
            residuals[i] = residual
        except Exception as e:
            if verbose and i == 0:
                log(f"⚠️  NNLS failed for cell {i}: {e}", level="WARNING")
            W[:, i] = 0
            residuals[i] = np.inf
        
        # Progress update
        if verbose and (i + 1) % progress_interval == 0:
            pct = 100 * (i + 1) / n_cells
            log(f"  Progress: {i+1:,}/{n_cells:,} cells ({pct:.1f}%)")
    
    if verbose:
        log(f"✓ Fitting complete for all {n_cells:,} cells")
    
    # Create weight DataFrame
    weights_df = pd.DataFrame(
        W,
        index=signature_matrix.columns,
        columns=mutation_matrix.columns
    )
    
    # Calculate reconstruction
    reconstruction = H @ W
    reconstruction_df = pd.DataFrame(
        reconstruction,
        index=mutation_matrix.index,
        columns=mutation_matrix.columns
    )
    
    results = {
        'weights': weights_df,
        'residuals': residuals,
        'reconstruction': reconstruction_df
    }
    
    if verbose:
        log("\n" + "="*80)
        log("SIGNATURE FITTING COMPLETE")
        log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 4: EVALUATE RECONSTRUCTION QUALITY
# =============================================================================

def evaluate_reconstruction(original_matrix, reconstructed_matrix, verbose=True):
    """
    Evaluate how well H × W approximates X using proper matrix norms
    
    Uses:
    - Frobenius norm for reconstruction error
    - Eckart-Young theorem for optimal low-rank approximation
    - Per-cell correlation metrics
    
    Metrics:
    - Frobenius norm ||X - (H×W)||_F (lower is better)
    - Relative Frobenius error ||X - (H×W)||_F / ||X||_F
    - Eckart-Young optimal rank-k approximation via SVD
    - Per-cell Pearson correlation
    - Per-cell Cosine similarity
    
    Parameters:
    -----------
    original_matrix : pd.DataFrame
        Original mutation matrix X
    reconstructed_matrix : pd.DataFrame
        Reconstructed matrix H × W
    verbose : bool
        Print detailed statistics
    
    Returns:
    --------
    dict with evaluation metrics
    """
    if verbose:
        log("="*80)
        log("EVALUATING RECONSTRUCTION QUALITY")
        log("="*80)
    
    X = original_matrix.values
    X_recon = reconstructed_matrix.values
    
    n_contexts = X.shape[0]
    n_cells = X.shape[1]
    
    if verbose:
        log(f"\nMatrix dimensions: {n_contexts} contexts × {n_cells} cells")
    
    # -------------------------------------------------------------------------
    # Calculate difference matrix (residual)
    # -------------------------------------------------------------------------
    if verbose:
        log("\nCalculating residual matrix (X - X_reconstructed)...")
    
    residual_matrix = X - X_recon
    
    # -------------------------------------------------------------------------
    # Frobenius Norm Analysis
    # -------------------------------------------------------------------------
    if verbose:
        log("\nCalculating Frobenius norms...")
    
    # ||X||_F
    frobenius_norm_original = np.linalg.norm(X, 'fro')
    
    # ||X_recon||_F
    frobenius_norm_reconstructed = np.linalg.norm(X_recon, 'fro')
    
    # ||X - X_recon||_F (reconstruction error)
    frobenius_error = np.linalg.norm(residual_matrix, 'fro')
    
    # Relative error: ||X - X_recon||_F / ||X||_F
    relative_frobenius_error = frobenius_error / frobenius_norm_original if frobenius_norm_original > 0 else 0
    
    if verbose:
        log(f"  ||X||_F (original):         {frobenius_norm_original:.2f}")
        log(f"  ||X_recon||_F:              {frobenius_norm_reconstructed:.2f}")
        log(f"  ||X - X_recon||_F (error):  {frobenius_error:.2f}")
        log(f"  Relative error:             {100*relative_frobenius_error:.2f}%")
    
    # -------------------------------------------------------------------------
    # Eckart-Young Theorem: Optimal Low-Rank Approximation
    # -------------------------------------------------------------------------
    if verbose:
        log("\nApplying Eckart-Young theorem (optimal rank-k approximation via SVD)...")
    
    # Perform SVD on original matrix: X = U Σ V^T
    U, singular_values, Vt = np.linalg.svd(X, full_matrices=False)
    
    # Rank of the matrix
    matrix_rank = np.linalg.matrix_rank(X)
    
    if verbose:
        log(f"  Matrix rank: {matrix_rank}")
        log(f"  Top 10 singular values: {singular_values[:10]}")
    
    # Calculate cumulative explained variance
    total_variance = np.sum(singular_values**2)
    cumulative_variance = np.cumsum(singular_values**2) / total_variance
    
    # Find rank needed for 90%, 95%, 99% reconstruction
    rank_90 = np.searchsorted(cumulative_variance, 0.90) + 1
    rank_95 = np.searchsorted(cumulative_variance, 0.95) + 1
    rank_99 = np.searchsorted(cumulative_variance, 0.99) + 1
    
    if verbose:
        log(f"\n  Eckart-Young optimal ranks:")
        log(f"    90% variance: rank-{rank_90}")
        log(f"    95% variance: rank-{rank_95}")
        log(f"    99% variance: rank-{rank_99}")
    
    # Calculate theoretical minimum error for rank-k approximation
    # More accurate: count non-zero singular values in reconstructed matrix
    _, sv_recon, _ = np.linalg.svd(X_recon, full_matrices=False)
    decomposition_rank_effective = np.sum(sv_recon > 1e-10)
    
    if decomposition_rank_effective < len(singular_values):
        # Theoretical minimum error for this rank
        theoretical_min_error = np.sqrt(np.sum(singular_values[decomposition_rank_effective:]**2))
        
        # How close are we to optimal?
        optimality_ratio = theoretical_min_error / frobenius_error if frobenius_error > 0 else 1.0
        
        if verbose:
            log(f"\n  Effective decomposition rank: {decomposition_rank_effective}")
            log(f"  Theoretical minimum error (Eckart-Young): {theoretical_min_error:.2f}")
            log(f"  Actual error: {frobenius_error:.2f}")
            log(f"  Optimality ratio: {optimality_ratio:.4f}")
            
            if optimality_ratio > 0.95:
                log(f"    ✓ Near-optimal reconstruction (>95% of theoretical best)")
            elif optimality_ratio > 0.80:
                log(f"    ○ Good reconstruction (>80% of theoretical best)")
            else:
                log(f"    ⚠ Suboptimal reconstruction (<80% of theoretical best)")
    else:
        theoretical_min_error = 0
        optimality_ratio = 1.0
        if verbose:
            log(f"\n  Effective decomposition rank: {decomposition_rank_effective} (full rank)")
    
    # -------------------------------------------------------------------------
    # Per-cell metrics
    # -------------------------------------------------------------------------
    if verbose:
        log("\nCalculating per-cell metrics...")
    
    pearson_corrs = []
    spearman_corrs = []
    cosine_sims = []
    cell_frobenius_errors = []
    
    for i in range(n_cells):
        x_cell = X[:, i]
        x_recon_cell = X_recon[:, i]
        
        # Frobenius error for this cell (L2 norm of difference)
        cell_error = np.linalg.norm(x_cell - x_recon_cell)
        cell_frobenius_errors.append(cell_error)
        
        # Pearson correlation
        if x_cell.sum() > 0 and x_recon_cell.sum() > 0:
            r, _ = pearsonr(x_cell, x_recon_cell)
            pearson_corrs.append(r)
            
            # Spearman correlation
            rho, _ = spearmanr(x_cell, x_recon_cell)
            spearman_corrs.append(rho)
            
            # Cosine similarity
            norm_x = np.linalg.norm(x_cell)
            norm_recon = np.linalg.norm(x_recon_cell)
            if norm_x > 0 and norm_recon > 0:
                cos_sim = np.dot(x_cell, x_recon_cell) / (norm_x * norm_recon)
                cosine_sims.append(cos_sim)
            else:
                cosine_sims.append(np.nan)
        else:
            pearson_corrs.append(np.nan)
            spearman_corrs.append(np.nan)
            cosine_sims.append(np.nan)
    
    pearson_corrs = np.array(pearson_corrs)
    spearman_corrs = np.array(spearman_corrs)
    cosine_sims = np.array(cosine_sims)
    cell_frobenius_errors = np.array(cell_frobenius_errors)
    
    # -------------------------------------------------------------------------
    # Additional metrics
    # -------------------------------------------------------------------------
    
    # Reconstruction rate (what % of mutations are explained)
    total_mutations_original = X.sum()
    total_mutations_reconstructed = X_recon.sum()
    reconstruction_rate = total_mutations_reconstructed / total_mutations_original if total_mutations_original > 0 else 0
    
    # Mean Absolute Error per element
    mae = np.mean(np.abs(residual_matrix))
    
    # Root Mean Squared Error per element
    rmse = np.sqrt(np.mean(residual_matrix**2))
    
    # -------------------------------------------------------------------------
    # Summary statistics
    # -------------------------------------------------------------------------
    if verbose:
        log("\n" + "="*60)
        log("RECONSTRUCTION QUALITY METRICS")
        log("="*60)
        
        log("\nMatrix Norm Metrics:")
        log(f"  Frobenius error:          {frobenius_error:.2f}")
        log(f"  Relative Frobenius error: {100*relative_frobenius_error:.2f}%")
        log(f"  MAE (per element):        {mae:.4f}")
        log(f"  RMSE (per element):       {rmse:.4f}")
        
        log("\nEckart-Young Optimal Approximation:")
        log(f"  Theoretical min error:    {theoretical_min_error:.2f}")
        log(f"  Optimality ratio:         {optimality_ratio:.4f}")
        log(f"  Effective rank:           {decomposition_rank_effective}")
        
        log("\nPer-Cell Metrics (averaged across cells):")
        log(f"  Mean cell error (L2):     {np.mean(cell_frobenius_errors):.2f}")
        log(f"  Median cell error:        {np.median(cell_frobenius_errors):.2f}")
        
        log(f"\n  Pearson correlation:")
        log(f"    Mean:   {np.nanmean(pearson_corrs):.4f}")
        log(f"    Median: {np.nanmedian(pearson_corrs):.4f}")
        log(f"    Std:    {np.nanstd(pearson_corrs):.4f}")
        
        log(f"\n  Spearman correlation:")
        log(f"    Mean:   {np.nanmean(spearman_corrs):.4f}")
        log(f"    Median: {np.nanmedian(spearman_corrs):.4f}")
        
        log(f"\n  Cosine similarity:")
        log(f"    Mean:   {np.nanmean(cosine_sims):.4f}")
        log(f"    Median: {np.nanmedian(cosine_sims):.4f}")
        
        log("\nMutation Counts:")
        log(f"  Original total:       {int(total_mutations_original):,}")
        log(f"  Reconstructed total:  {int(total_mutations_reconstructed):,}")
        log(f"  Difference:           {int(total_mutations_original - total_mutations_reconstructed):,}")
        log(f"  Reconstruction rate:  {100*reconstruction_rate:.2f}%")
    
    # Quality assessment based on Frobenius error
    if verbose:
        log("\nQuality Assessment:")
        if relative_frobenius_error < 0.10:
            quality = "EXCELLENT"
        elif relative_frobenius_error < 0.20:
            quality = "GOOD"
        elif relative_frobenius_error < 0.30:
            quality = "MODERATE"
        else:
            quality = "POOR"
        log(f"  Overall quality (Frobenius): {quality}")
        
        # Also assess by correlation
        mean_pearson = np.nanmean(pearson_corrs)
        if mean_pearson > 0.8:
            quality_corr = "EXCELLENT"
        elif mean_pearson > 0.6:
            quality_corr = "GOOD"
        elif mean_pearson > 0.4:
            quality_corr = "MODERATE"
        else:
            quality_corr = "POOR"
        log(f"  Overall quality (Correlation): {quality_corr}")
    
    results = {
        # Matrix difference
        'residual_matrix': residual_matrix,
        
        # Frobenius norms
        'frobenius_norm_original': frobenius_norm_original,
        'frobenius_norm_reconstructed': frobenius_norm_reconstructed,
        'frobenius_error': frobenius_error,
        'relative_frobenius_error': relative_frobenius_error,
        
        # Eckart-Young theorem results
        'singular_values': singular_values,
        'matrix_rank': matrix_rank,
        'decomposition_rank': decomposition_rank_effective,
        'theoretical_min_error': theoretical_min_error,
        'optimality_ratio': optimality_ratio,
        'rank_90pct': rank_90,
        'rank_95pct': rank_95,
        'rank_99pct': rank_99,
        'cumulative_variance': cumulative_variance,
        
        # Per-cell metrics
        'pearson_per_cell': pearson_corrs,
        'spearman_per_cell': spearman_corrs,
        'cosine_per_cell': cosine_sims,
        'cell_frobenius_errors': cell_frobenius_errors,
        
        # Summary statistics
        'mean_pearson': np.nanmean(pearson_corrs),
        'median_pearson': np.nanmedian(pearson_corrs),
        'std_pearson': np.nanstd(pearson_corrs),
        'mean_cosine': np.nanmean(cosine_sims),
        'median_cosine': np.nanmedian(cosine_sims),
        'mean_cell_error': np.mean(cell_frobenius_errors),
        
        # Other metrics
        'mae': mae,
        'rmse': rmse,
        'reconstruction_rate': reconstruction_rate,
        'total_mutations_original': total_mutations_original,
        'total_mutations_reconstructed': total_mutations_reconstructed,
        'quality': quality if verbose else None,
        'quality_correlation': quality_corr if verbose else None
    }
    
    if verbose:
        log("\n" + "="*80)
        log("EVALUATION COMPLETE")
        log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 5: VISUALIZATION
# =============================================================================

def plot_reconstruction_quality(evaluation_results, output_dir):
    """
    Visualize reconstruction quality metrics including SVD analysis
    """
    log("="*80)
    log("GENERATING RECONSTRUCTION QUALITY PLOTS")
    log("="*80)
    
    fig = plt.figure(figsize=(20, 14))
    
    pearson_corrs = evaluation_results['pearson_per_cell']
    cosine_sims = evaluation_results['cosine_per_cell']
    cell_errors = evaluation_results['cell_frobenius_errors']
    
    # Remove NaN values for plotting
    pearson_valid = pearson_corrs[~np.isnan(pearson_corrs)]
    cosine_valid = cosine_sims[~np.isnan(cosine_sims)]
    
    # -------------------------------------------------------------------------
    # Panel 1: Pearson correlation distribution
    # -------------------------------------------------------------------------
    ax1 = plt.subplot(3, 3, 1)
    ax1.hist(pearson_valid, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axvline(np.nanmean(pearson_corrs), color='red', linestyle='--', 
                linewidth=2, label=f'Mean: {np.nanmean(pearson_corrs):.3f}')
    ax1.axvline(np.nanmedian(pearson_corrs), color='orange', linestyle='--',
                linewidth=2, label=f'Median: {np.nanmedian(pearson_corrs):.3f}')
    ax1.set_xlabel('Pearson Correlation', fontsize=11)
    ax1.set_ylabel('Number of Cells', fontsize=11)
    ax1.set_title('Per-Cell Reconstruction Quality\n(Pearson Correlation)', 
                  fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 2: Cosine similarity distribution
    # -------------------------------------------------------------------------
    ax2 = plt.subplot(3, 3, 2)
    ax2.hist(cosine_valid, bins=50, color='darkgreen', alpha=0.7, edgecolor='black')
    ax2.axvline(np.nanmean(cosine_sims), color='red', linestyle='--',
                linewidth=2, label=f'Mean: {np.nanmean(cosine_sims):.3f}')
    ax2.axvline(np.nanmedian(cosine_sims), color='orange', linestyle='--',
                linewidth=2, label=f'Median: {np.nanmedian(cosine_sims):.3f}')
    ax2.set_xlabel('Cosine Similarity', fontsize=11)
    ax2.set_ylabel('Number of Cells', fontsize=11)
    ax2.set_title('Per-Cell Reconstruction Quality\n(Cosine Similarity)',
                  fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(axis='y', alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 3: Per-cell Frobenius error distribution
    # -------------------------------------------------------------------------
    ax3 = plt.subplot(3, 3, 3)
    ax3.hist(cell_errors, bins=50, color='coral', alpha=0.7, edgecolor='black')
    ax3.axvline(np.mean(cell_errors), color='red', linestyle='--',
                linewidth=2, label=f'Mean: {np.mean(cell_errors):.2f}')
    ax3.axvline(np.median(cell_errors), color='orange', linestyle='--',
                linewidth=2, label=f'Median: {np.median(cell_errors):.2f}')
    ax3.set_xlabel('Per-Cell Frobenius Error (L2 norm)', fontsize=11)
    ax3.set_ylabel('Number of Cells', fontsize=11)
    ax3.set_title('Per-Cell Reconstruction Error\n(Frobenius Norm)',
                  fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 4: Quality categories
    # -------------------------------------------------------------------------
    ax4 = plt.subplot(3, 3, 4)
    
    # Categorize cells
    excellent = np.sum(pearson_valid > 0.8)
    good = np.sum((pearson_valid > 0.6) & (pearson_valid <= 0.8))
    moderate = np.sum((pearson_valid > 0.4) & (pearson_valid <= 0.6))
    poor = np.sum(pearson_valid <= 0.4)
    
    categories = ['Excellent\n(>0.8)', 'Good\n(0.6-0.8)', 'Moderate\n(0.4-0.6)', 'Poor\n(<0.4)']
    counts = [excellent, good, moderate, poor]
    colors_cat = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c']
    
    bars = ax4.bar(categories, counts, color=colors_cat, alpha=0.7, edgecolor='black')
    ax4.set_ylabel('Number of Cells', fontsize=11)
    ax4.set_title('Reconstruction Quality Categories\n(Pearson Correlation)',
                  fontsize=12, fontweight='bold')
    ax4.grid(axis='y', alpha=0.3)
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}\n({100*count/len(pearson_valid):.1f}%)',
                ha='center', va='bottom', fontsize=9)
    
    # -------------------------------------------------------------------------
    # Panel 5: Pearson vs Cosine scatter
    # -------------------------------------------------------------------------
    ax5 = plt.subplot(3, 3, 5)
    
    # Match up valid values
    valid_mask = ~np.isnan(pearson_corrs) & ~np.isnan(cosine_sims)
    pearson_for_scatter = pearson_corrs[valid_mask]
    cosine_for_scatter = cosine_sims[valid_mask]
    
    ax5.scatter(pearson_for_scatter, cosine_for_scatter, alpha=0.3, s=10, color='purple')
    ax5.plot([0, 1], [0, 1], 'r--', linewidth=1, alpha=0.5, label='y=x')
    ax5.set_xlabel('Pearson Correlation', fontsize=11)
    ax5.set_ylabel('Cosine Similarity', fontsize=11)
    ax5.set_title('Pearson vs Cosine Similarity',
                  fontsize=12, fontweight='bold')
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim([0, 1])
    ax5.set_ylim([0, 1])
    
    # -------------------------------------------------------------------------
    # Panel 6: Cumulative distribution
    # -------------------------------------------------------------------------
    ax6 = plt.subplot(3, 3, 6)
    
    sorted_pearson = np.sort(pearson_valid)
    cumulative = np.arange(1, len(sorted_pearson) + 1) / len(sorted_pearson)
    
    ax6.plot(sorted_pearson, cumulative, linewidth=2, color='steelblue')
    ax6.axvline(0.5, color='orange', linestyle='--', alpha=0.5, label='r=0.5')
    ax6.axvline(0.7, color='red', linestyle='--', alpha=0.5, label='r=0.7')
    ax6.set_xlabel('Pearson Correlation', fontsize=11)
    ax6.set_ylabel('Cumulative Fraction of Cells', fontsize=11)
    ax6.set_title('Cumulative Distribution\n(Pearson Correlation)',
                  fontsize=12, fontweight='bold')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    ax6.set_xlim([0, 1])
    ax6.set_ylim([0, 1])
    
    # -------------------------------------------------------------------------
    # Panel 7: SVD Scree Plot (Singular Values)
    # -------------------------------------------------------------------------
    ax7 = plt.subplot(3, 3, 7)
    
    singular_values = evaluation_results['singular_values']
    n_show = min(50, len(singular_values))  # Show first 50
    
    ax7.plot(range(1, n_show + 1), singular_values[:n_show], 'o-', 
             linewidth=2, markersize=4, color='navy')
    ax7.set_xlabel('Rank', fontsize=11)
    ax7.set_ylabel('Singular Value', fontsize=11)
    ax7.set_title('SVD Scree Plot\n(Singular Values)',
                  fontsize=12, fontweight='bold')
    ax7.grid(True, alpha=0.3)
    ax7.set_yscale('log')
    
    # Mark decomposition rank
    decomp_rank = evaluation_results['decomposition_rank']
    if decomp_rank <= n_show:
        ax7.axvline(decomp_rank, color='red', linestyle='--', linewidth=2,
                   label=f'Decomp. rank: {decomp_rank}')
        ax7.legend(fontsize=9)
    
    # -------------------------------------------------------------------------
    # Panel 8: Cumulative Variance Explained (Eckart-Young)
    # -------------------------------------------------------------------------
    ax8 = plt.subplot(3, 3, 8)
    
    cumulative_variance = evaluation_results['cumulative_variance']
    
    ax8.plot(range(1, len(cumulative_variance) + 1), cumulative_variance * 100,
             linewidth=2, color='darkgreen')
    ax8.axhline(90, color='orange', linestyle='--', alpha=0.5, label='90%')
    ax8.axhline(95, color='red', linestyle='--', alpha=0.5, label='95%')
    ax8.axhline(99, color='purple', linestyle='--', alpha=0.5, label='99%')
    
    # Mark optimal ranks
    rank_90 = evaluation_results['rank_90pct']
    rank_95 = evaluation_results['rank_95pct']
    rank_99 = evaluation_results['rank_99pct']
    
    ax8.axvline(rank_90, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
    ax8.axvline(rank_95, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
    ax8.axvline(rank_99, color='purple', linestyle=':', linewidth=1.5, alpha=0.7)
    
    ax8.set_xlabel('Rank', fontsize=11)
    ax8.set_ylabel('Cumulative Variance Explained (%)', fontsize=11)
    ax8.set_title('Eckart-Young Optimal Rank\n(Cumulative Variance)',
                  fontsize=12, fontweight='bold')
    ax8.legend(fontsize=9, loc='lower right')
    ax8.grid(True, alpha=0.3)
    ax8.set_xlim([0, min(100, len(cumulative_variance))])
    ax8.set_ylim([0, 100])
    
    # -------------------------------------------------------------------------
    # Panel 9: Summary statistics
    # -------------------------------------------------------------------------
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    
    summary_text = f"""
    RECONSTRUCTION METRICS SUMMARY
    {'='*45}
    
    FROBENIUS NORMS:
      ||X||_F:           {evaluation_results['frobenius_norm_original']:>10.2f}
      ||X_recon||_F:     {evaluation_results['frobenius_norm_reconstructed']:>10.2f}
      ||Error||_F:       {evaluation_results['frobenius_error']:>10.2f}
      Relative error:    {100*evaluation_results['relative_frobenius_error']:>9.2f}%
    
    ECKART-YOUNG THEOREM:
      Matrix rank:       {evaluation_results['matrix_rank']:>10}
      Decomp. rank:      {evaluation_results['decomposition_rank']:>10}
      Theoretical min:   {evaluation_results['theoretical_min_error']:>10.2f}
      Optimality:        {evaluation_results['optimality_ratio']:>10.4f}
    
    CORRELATION METRICS:
      Mean Pearson:      {evaluation_results['mean_pearson']:>10.4f}
      Median Pearson:    {evaluation_results['median_pearson']:>10.4f}
      Mean Cosine:       {evaluation_results['mean_cosine']:>10.4f}
    
    MUTATION COUNTS:
      Original:    {int(evaluation_results['total_mutations_original']):>15,}
      Reconstructed:{int(evaluation_results['total_mutations_reconstructed']):>14,}
      Difference:  {int(evaluation_results['total_mutations_original'] - evaluation_results['total_mutations_reconstructed']):>15,}
    
    QUALITY ASSESSMENT:
      Frobenius:   {evaluation_results['quality']:>15s}
      Correlation: {evaluation_results['quality_correlation']:>15s}
    """
    
    ax9.text(0.05, 0.5, summary_text, fontsize=9, family='monospace',
             verticalalignment='center', transform=ax9.transAxes)
    
    plt.suptitle(f'Reconstruction Quality: H × W ≈ X ({len(pearson_valid):,} cells)',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    save_figure(fig, output_dir, 'reconstruction_quality.png')
    plt.show()
    
    log("\n" + "="*80)
    log("QUALITY PLOTS COMPLETE")
    log("="*80 + "\n")


def plot_signature_weights_summary(weights_df, output_dir, top_n=15):
    """
    Visualize distribution of signature weights across cells
    """
    log("="*80)
    log("GENERATING SIGNATURE WEIGHT PLOTS")
    log("="*80)
    
    n_sigs = weights_df.shape[0]
    n_cells = weights_df.shape[1]
    
    fig = plt.figure(figsize=(18, 10))
    
    # -------------------------------------------------------------------------
    # Panel 1: Mean weights per signature
    # -------------------------------------------------------------------------
    ax1 = plt.subplot(2, 2, 1)
    
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(mean_weights)))
    bars = ax1.barh(range(len(mean_weights)), mean_weights.values, color=colors)
    ax1.set_yticks(range(len(mean_weights)))
    ax1.set_yticklabels(mean_weights.index, fontsize=9)
    ax1.set_xlabel('Mean Weight Across Cells', fontsize=11)
    ax1.set_title(f'Average Signature Activity\n({n_cells:,} cells)',
                  fontsize=12, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    ax1.invert_yaxis()
    
    # -------------------------------------------------------------------------
    # Panel 2: Frequency of signature detection
    # -------------------------------------------------------------------------
    ax2 = plt.subplot(2, 2, 2)
    
    # Count cells with weight > 0.1 for each signature
    detection_threshold = 0.1
    detection_freq = (weights_df > detection_threshold).sum(axis=1) / n_cells * 100
    detection_freq = detection_freq.sort_values(ascending=False)
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(detection_freq)))
    bars = ax2.barh(range(len(detection_freq)), detection_freq.values, color=colors)
    ax2.set_yticks(range(len(detection_freq)))
    ax2.set_yticklabels(detection_freq.index, fontsize=9)
    ax2.set_xlabel('% of Cells with Weight > 0.1', fontsize=11)
    ax2.set_title(f'Signature Detection Frequency\n(threshold = {detection_threshold})',
                  fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    ax2.invert_yaxis()
    
    # -------------------------------------------------------------------------
    # Panel 3: Weight distribution heatmap (top signatures)
    # -------------------------------------------------------------------------
    ax3 = plt.subplot(2, 2, 3)
    
    # Select top N signatures by mean weight
    top_sigs = mean_weights.head(min(top_n, n_sigs)).index
    weights_top = weights_df.loc[top_sigs]
    
    # Sample cells if too many
    if n_cells > 1000:
        cell_sample = np.random.choice(weights_df.columns, 1000, replace=False)
        weights_plot = weights_top[cell_sample]
        title_suffix = f" (1000 random cells)"
    else:
        weights_plot = weights_top
        title_suffix = f" ({n_cells} cells)"
    
    im = ax3.imshow(weights_plot.values, aspect='auto', cmap='YlOrRd', 
                    interpolation='nearest')
    ax3.set_yticks(range(len(top_sigs)))
    ax3.set_yticklabels(top_sigs, fontsize=9)
    ax3.set_xlabel('Cells (sampled)', fontsize=11)
    ax3.set_title(f'Signature Weight Heatmap (Top {len(top_sigs)}){title_suffix}',
                  fontsize=12, fontweight='bold')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax3)
    cbar.set_label('Weight', fontsize=10)
    
    # -------------------------------------------------------------------------
    # Panel 4: Total mutations per cell
    # -------------------------------------------------------------------------
    ax4 = plt.subplot(2, 2, 4)
    
    total_mutations_per_cell = weights_df.sum(axis=0)
    
    ax4.hist(total_mutations_per_cell, bins=50, color='teal', 
             alpha=0.7, edgecolor='black')
    ax4.axvline(total_mutations_per_cell.mean(), color='red', 
                linestyle='--', linewidth=2, 
                label=f'Mean: {total_mutations_per_cell.mean():.1f}')
    ax4.axvline(total_mutations_per_cell.median(), color='orange',
                linestyle='--', linewidth=2,
                label=f'Median: {total_mutations_per_cell.median():.1f}')
    ax4.set_xlabel('Total Signature Weight (Sum Across Signatures)', fontsize=11)
    ax4.set_ylabel('Number of Cells', fontsize=11)
    ax4.set_title('Total Signature Activity Per Cell',
                  fontsize=12, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    save_figure(fig, output_dir, 'signature_weights_summary.png')
    plt.show()
    
    log("\n" + "="*80)
    log("WEIGHT PLOTS COMPLETE")
    log("="*80 + "\n")


# =============================================================================
# STEP 6: SAVE RESULTS
# =============================================================================

def save_refitting_results(weights_df, reconstruction_df, evaluation_results,
                           signature_matrix, mutation_matrix, output_dir):
    """
    Save all refitting results to files
    """
    log("="*80)
    log("SAVING REFITTING RESULTS")
    log("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # -------------------------------------------------------------------------
    # 1. Signature weights (W matrix)
    # -------------------------------------------------------------------------
    weights_file = output_path / "signature_weights_per_cell.txt"
    weights_df.to_csv(weights_file, sep='\t', float_format='%.6f')
    log(f"Saved weights: {weights_file}")
    
    # -------------------------------------------------------------------------
    # 2. Reconstructed mutation matrix (H × W)
    # -------------------------------------------------------------------------
    recon_file = output_path / "reconstructed_mutation_matrix.txt"
    reconstruction_df.to_csv(recon_file, sep='\t', float_format='%.2f')
    log(f"Saved reconstruction: {recon_file}")
    
    # -------------------------------------------------------------------------
    # 3. Residual matrix (X - X_reconstructed)
    # -------------------------------------------------------------------------
    residual_matrix = evaluation_results['residual_matrix']
    residual_df = pd.DataFrame(
        residual_matrix,
        index=mutation_matrix.index,
        columns=mutation_matrix.columns
    )
    residual_file = output_path / "residual_matrix_X_minus_Xrecon.txt"
    residual_df.to_csv(residual_file, sep='\t', float_format='%.2f')
    log(f"Saved residual matrix: {residual_file}")
    
    # -------------------------------------------------------------------------
    # 4. Evaluation metrics
    # -------------------------------------------------------------------------
    eval_file = output_path / "reconstruction_evaluation.txt"
    
    with open(eval_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SIGNATURE REFITTING EVALUATION\n")
        f.write("="*80 + "\n\n")
        
        f.write("Model:\n")
        f.write(f"  X ≈ H × W\n\n")
        f.write(f"Dimensions:\n")
        f.write(f"  X: {mutation_matrix.shape[0]} contexts × {mutation_matrix.shape[1]} cells\n")
        f.write(f"  H: {signature_matrix.shape[0]} contexts × {signature_matrix.shape[1]} signatures (FIXED)\n")
        f.write(f"  W: {weights_df.shape[0]} signatures × {weights_df.shape[1]} cells (FITTED)\n\n")
        
        f.write("="*80 + "\n")
        f.write("MATRIX NORM METRICS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Frobenius Norms:\n")
        f.write(f"  ||X||_F (original):           {evaluation_results['frobenius_norm_original']:.2f}\n")
        f.write(f"  ||X_recon||_F:                {evaluation_results['frobenius_norm_reconstructed']:.2f}\n")
        f.write(f"  ||X - X_recon||_F (error):    {evaluation_results['frobenius_error']:.2f}\n")
        f.write(f"  Relative error:               {100*evaluation_results['relative_frobenius_error']:.2f}%\n\n")
        
        f.write("Element-wise Metrics:\n")
        f.write(f"  MAE (per element):            {evaluation_results['mae']:.4f}\n")
        f.write(f"  RMSE (per element):           {evaluation_results['rmse']:.4f}\n\n")
        
        f.write("="*80 + "\n")
        f.write("ECKART-YOUNG THEOREM (OPTIMAL LOW-RANK APPROXIMATION)\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Matrix Properties:\n")
        f.write(f"  Full matrix rank:             {evaluation_results['matrix_rank']}\n")
        f.write(f"  Effective decomposition rank: {evaluation_results['decomposition_rank']}\n\n")
        
        f.write(f"Optimal Ranks for Variance Thresholds:\n")
        f.write(f"  90% variance: rank-{evaluation_results['rank_90pct']}\n")
        f.write(f"  95% variance: rank-{evaluation_results['rank_95pct']}\n")
        f.write(f"  99% variance: rank-{evaluation_results['rank_99pct']}\n\n")
        
        f.write(f"Reconstruction Optimality:\n")
        f.write(f"  Theoretical min error:        {evaluation_results['theoretical_min_error']:.2f}\n")
        f.write(f"  Actual error:                 {evaluation_results['frobenius_error']:.2f}\n")
        f.write(f"  Optimality ratio:             {evaluation_results['optimality_ratio']:.4f}\n")
        
        if evaluation_results['optimality_ratio'] > 0.95:
            f.write(f"  Assessment: Near-optimal (>95% of theoretical best)\n\n")
        elif evaluation_results['optimality_ratio'] > 0.80:
            f.write(f"  Assessment: Good (>80% of theoretical best)\n\n")
        else:
            f.write(f"  Assessment: Suboptimal (<80% of theoretical best)\n\n")
        
        f.write("="*80 + "\n")
        f.write("PER-CELL QUALITY METRICS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Per-Cell Metrics (averaged):\n")
        f.write(f"  Mean cell error (L2 norm):    {evaluation_results['mean_cell_error']:.2f}\n\n")
        
        f.write(f"  Mean Pearson correlation:     {evaluation_results['mean_pearson']:.4f}\n")
        f.write(f"  Median Pearson correlation:   {evaluation_results['median_pearson']:.4f}\n")
        f.write(f"  Std Pearson correlation:      {evaluation_results['std_pearson']:.4f}\n\n")
        
        f.write(f"  Mean Cosine similarity:       {evaluation_results['mean_cosine']:.4f}\n")
        f.write(f"  Median Cosine similarity:     {evaluation_results['median_cosine']:.4f}\n\n")
        
        f.write("="*80 + "\n")
        f.write("MUTATION COUNTS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Mutation Counts:\n")
        f.write(f"  Original total:       {int(evaluation_results['total_mutations_original']):>12,}\n")
        f.write(f"  Reconstructed total:  {int(evaluation_results['total_mutations_reconstructed']):>12,}\n")
        f.write(f"  Difference:           {int(evaluation_results['total_mutations_original'] - evaluation_results['total_mutations_reconstructed']):>12,}\n")
        f.write(f"  Reconstruction rate:  {100*evaluation_results['reconstruction_rate']:>11.2f}%\n\n")
        
        f.write("="*80 + "\n")
        f.write("OVERALL QUALITY ASSESSMENT\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Based on Frobenius norm:   {evaluation_results['quality']}\n")
        f.write(f"Based on correlation:      {evaluation_results['quality_correlation']}\n")
    
    log(f"Saved evaluation: {eval_file}")
    
    # -------------------------------------------------------------------------
    # 5. SVD analysis (singular values and variance explained)
    # -------------------------------------------------------------------------
    svd_file = output_path / "svd_analysis_singular_values.txt"
    
    singular_values = evaluation_results['singular_values']
    cumulative_variance = evaluation_results['cumulative_variance']
    
    svd_df = pd.DataFrame({
        'rank': np.arange(1, len(singular_values) + 1),
        'singular_value': singular_values,
        'variance_explained': (singular_values**2) / np.sum(singular_values**2),
        'cumulative_variance': cumulative_variance
    })
    
    svd_df.to_csv(svd_file, sep='\t', index=False, float_format='%.6f')
    log(f"Saved SVD analysis: {svd_file}")
    
    # -------------------------------------------------------------------------
    # 6. Per-cell quality metrics
    # -------------------------------------------------------------------------
    quality_df = pd.DataFrame({
        'cell_id': mutation_matrix.columns,
        'frobenius_error': evaluation_results['cell_frobenius_errors'],
        'pearson_correlation': evaluation_results['pearson_per_cell'],
        'cosine_similarity': evaluation_results['cosine_per_cell']
    })
    
    quality_file = output_path / "per_cell_quality_metrics.txt"
    quality_df.to_csv(quality_file, sep='\t', index=False, float_format='%.4f')
    log(f"Saved per-cell metrics: {quality_file}")
    
    # -------------------------------------------------------------------------
    # 7. Signature weight summary
    # -------------------------------------------------------------------------
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    median_weights = weights_df.median(axis=1).sort_values(ascending=False)
    max_weights = weights_df.max(axis=1).sort_values(ascending=False)
    detection_freq = (weights_df > 0.1).sum(axis=1) / weights_df.shape[1] * 100
    
    summary_df = pd.DataFrame({
        'signature': mean_weights.index,
        'mean_weight': mean_weights.values,
        'median_weight': median_weights.loc[mean_weights.index].values,
        'max_weight': max_weights.loc[mean_weights.index].values,
        'detection_frequency_pct': detection_freq.loc[mean_weights.index].values
    })
    
    summary_file = output_path / "signature_weight_summary.txt"
    summary_df.to_csv(summary_file, sep='\t', index=False, float_format='%.4f')
    log(f"Saved weight summary: {summary_file}")
    
    log("\n" + "="*80)
    log("ALL RESULTS SAVED")
    log("="*80 + "\n")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_signature_refitting_pipeline(mutation_matrix_file, cosmic_file, output_dir,
                                     mutation_threshold=0, use_scree_plot=True,
                                     core_signatures=None, candidate_order=None):
    """
    Complete semi-supervised signature refitting pipeline
    
    Steps:
    1. Extract HNSCC-relevant COSMIC signatures (H matrix)
    2. Load and optionally filter mutation matrix by mutation count
    3. Select optimal signatures using scree plot elbow detection
    4. Fit signatures to mutation data using NNLS (solve for W)
    5. Evaluate reconstruction quality (Frobenius norm, Eckart-Young)
    6. Visualize results
    7. Save everything
    
    Parameters:
    -----------
    mutation_matrix_file : str
        Path to mutation count matrix (96 contexts × cells)
    cosmic_file : str
        Path to COSMIC signature database
    output_dir : str
        Output directory for all results
    mutation_threshold : int
        Minimum mutations per cell to include (default: 0, no filtering)
        Set to >0 to filter out sparse cells
    use_scree_plot : bool
        If True, use scree plot for signature selection (default: True)
        If False, use all HNSCC signatures
    core_signatures : list, optional
        Core signatures to always include (default: ['SBS2', 'SBS3', 'SBS5'])
    candidate_order : list, optional
        Order to try adding candidate signatures (by mean weight)
    
    Returns:
    --------
    dict with all results
    """
    
    print("\n" + "="*80)
    print("SEMI-SUPERVISED SIGNATURE REFITTING PIPELINE")
    print("Using Known HNSCC COSMIC Signatures")
    print("="*80)
    print(f"\nMutation matrix: {mutation_matrix_file}")
    print(f"COSMIC database: {cosmic_file}")
    print(f"Output: {output_dir}")
    print(f"Mutation threshold: {'No filtering' if mutation_threshold == 0 else f'≥{mutation_threshold} mutations/cell'}")
    print(f"Scree plot selection: {'Enabled' if use_scree_plot else 'Disabled (using all signatures)'}")
    print("="*80 + "\n")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # -------------------------------------------------------------------------
    # STEP 1: Extract HNSCC signatures
    # -------------------------------------------------------------------------
    hnscc_sigs = extract_hnscc_signatures(cosmic_file, output_dir)
    
    # -------------------------------------------------------------------------
    # STEP 2: Load mutation matrix
    # -------------------------------------------------------------------------
    log("="*80)
    log("LOADING MUTATION MATRIX")
    log("="*80)
    
    log(f"\nLoading: {mutation_matrix_file}")
    mutation_matrix_full = pd.read_csv(mutation_matrix_file, sep='\t', index_col=0)
    
    log(f"Original matrix: {mutation_matrix_full.shape[0]} contexts × {mutation_matrix_full.shape[1]:,} cells")
    
    # -------------------------------------------------------------------------
    # STEP 2B: Filter by mutation count (if threshold > 0)
    # -------------------------------------------------------------------------
    if mutation_threshold > 0:
        log(f"\n{'='*70}")
        log(f"FILTERING CELLS BY MUTATION COUNT (threshold ≥{mutation_threshold})")
        log(f"{'='*70}")
        
        # Calculate mutations per cell
        mutations_per_cell = mutation_matrix_full.sum(axis=0)
        
        log(f"\nOriginal statistics:")
        log(f"  Total cells: {len(mutations_per_cell):,}")
        log(f"  Mean mutations/cell: {mutations_per_cell.mean():.2f}")
        log(f"  Median mutations/cell: {mutations_per_cell.median():.0f}")
        log(f"  Min: {mutations_per_cell.min():.0f}, Max: {mutations_per_cell.max():.0f}")
        
        # Filter
        cells_keep = mutations_per_cell >= mutation_threshold
        mutation_matrix = mutation_matrix_full.loc[:, cells_keep]
        mutations_filtered = mutations_per_cell[cells_keep]
        
        n_kept = mutation_matrix.shape[1]
        n_removed = mutation_matrix_full.shape[1] - n_kept
        pct_kept = 100 * n_kept / mutation_matrix_full.shape[1]
        
        log(f"\nFiltered statistics:")
        log(f"  Cells kept: {n_kept:,} ({pct_kept:.2f}%)")
        log(f"  Cells removed: {n_removed:,}")
        log(f"  Mean mutations/cell: {mutations_filtered.mean():.2f}")
        log(f"  Median mutations/cell: {mutations_filtered.median():.0f}")
        
        # Calculate sparsity
        sparsity = (mutation_matrix == 0).sum().sum() / (mutation_matrix.shape[0] * mutation_matrix.shape[1])
        log(f"  Sparsity: {100*sparsity:.2f}%")
        
        # Save filtered matrix
        filtered_matrix_file = output_path / f"filtered_mutation_matrix_min{mutation_threshold}muts.txt"
        mutation_matrix.to_csv(filtered_matrix_file, sep='\t', float_format='%.0f')
        log(f"\nSaved filtered matrix: {filtered_matrix_file}")
        
        log("\n" + "="*80)
        log("FILTERING COMPLETE")
        log("="*80 + "\n")
    else:
        log(f"\nNo filtering applied (using all {mutation_matrix_full.shape[1]:,} cells)")
        mutation_matrix = mutation_matrix_full
        mutations_per_cell = mutation_matrix.sum(axis=0)
        log(f"Mean mutations/cell: {mutations_per_cell.mean():.2f}")
        log(f"Median mutations/cell: {mutations_per_cell.median():.0f}")
    
    log(f"\nTotal mutations: {int(mutation_matrix.sum().sum()):,}")
    
    log("\n" + "="*80)
    log("MUTATION MATRIX LOADED")
    log("="*80 + "\n")
    
    # -------------------------------------------------------------------------
    # STEP 3: Signature selection via scree plot
    # -------------------------------------------------------------------------
    if use_scree_plot:
        # Set defaults if not provided
        if core_signatures is None:
            core_signatures = ['SBS2', 'SBS3', 'SBS5']
        
        if candidate_order is None:
            # Default order by mean weight
            candidate_order = ['SBS33', 'SBS7b', 'SBS1', 'SBS4', 'SBS17a', 'SBS17b', 
                             'SBS7d', 'SBS7a', 'SBS16', 'SBS18', 'SBS13']
        
        # Add SBS40a if it was used instead of SBS40
        if 'SBS40a' in hnscc_sigs.columns and 'SBS40a' not in candidate_order:
            candidate_order = ['SBS40a'] + candidate_order
        
        # Run scree plot selection
        selection_results = select_signatures_via_scree_plot(
            mutation_matrix,
            hnscc_sigs,
            core_signatures,
            candidate_order,
            output_dir,
            max_signatures=15,
            verbose=True
        )
        
        # Use selected signatures
        hnscc_sigs_final = selection_results['signature_matrix']
        
        log(f"Using {selection_results['n_signatures']} selected signatures via scree plot")
    else:
        # Use all HNSCC signatures
        hnscc_sigs_final = hnscc_sigs
        selection_results = None
        log(f"Using all {hnscc_sigs.shape[1]} HNSCC signatures (no scree plot selection)")
    
    # -------------------------------------------------------------------------
    # STEP 4: Fit signatures using NNLS
    # -------------------------------------------------------------------------
    fitting_results = fit_signatures_nnls(mutation_matrix, hnscc_sigs_final, verbose=True)
    
    weights_df = fitting_results['weights']
    reconstruction_df = fitting_results['reconstruction']
    
    # -------------------------------------------------------------------------
    # STEP 5: Evaluate reconstruction quality
    # -------------------------------------------------------------------------
    evaluation_results = evaluate_reconstruction(
        mutation_matrix, 
        reconstruction_df,
        verbose=True
    )
    
    # -------------------------------------------------------------------------
    # STEP 6: Visualize results
    # -------------------------------------------------------------------------
    plot_reconstruction_quality(evaluation_results, output_dir)
    plot_signature_weights_summary(weights_df, output_dir)
    
    # -------------------------------------------------------------------------
    # STEP 7: Save results
    # -------------------------------------------------------------------------
    save_refitting_results(
        weights_df,
        reconstruction_df,
        evaluation_results,
        hnscc_sigs_final,
        mutation_matrix,
        output_dir
    )
    
    # -------------------------------------------------------------------------
    # Final summary
    # -------------------------------------------------------------------------
    print("\n" + "="*80)
    print("PIPELINE COMPLETE!")
    print("="*80)
    
    print(f"\nSignatures fitted: {hnscc_sigs_final.shape[1]}")
    if use_scree_plot and selection_results:
        print(f"  Core signatures: {len(core_signatures)}")
        print(f"  Selected via scree plot: {', '.join(selection_results['selected_signatures'])}")
    
    print(f"\nCells analyzed: {mutation_matrix.shape[1]:,}")
    if mutation_threshold > 0:
        print(f"  Original cells: {mutation_matrix_full.shape[1]:,}")
        print(f"  Filtered (≥{mutation_threshold} muts): {mutation_matrix.shape[1]:,}")
    
    print(f"\nReconstruction Quality:")
    print(f"  Mean Pearson correlation: {evaluation_results['mean_pearson']:.4f}")
    print(f"  Frobenius error: {evaluation_results['frobenius_error']:.2f}")
    print(f"  Relative error: {100*evaluation_results['relative_frobenius_error']:.2f}%")
    print(f"  Optimality ratio: {evaluation_results['optimality_ratio']:.4f}")
    print(f"  Overall quality (Frobenius): {evaluation_results['quality']}")
    print(f"  Overall quality (Correlation): {evaluation_results['quality_correlation']}")
    
    print(f"\nTop 5 signatures by mean weight:")
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    for i, (sig, weight) in enumerate(mean_weights.head(5).items(), 1):
        print(f"  {i}. {sig}: {weight:.4f}")
    
    print(f"\nAll results saved to: {output_dir}")
    print("="*80 + "\n")
    
    return {
        'hnscc_signatures': hnscc_sigs_final,
        'mutation_matrix': mutation_matrix,
        'weights': weights_df,
        'reconstruction': reconstruction_df,
        'evaluation': evaluation_results,
        'selection': selection_results,
        'mutation_threshold': mutation_threshold
    }

#%%
# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    
    # =========================================================================
    # CONFIGURATION
    # =========================================================================
    
    MUTATION_MATRIX_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/SNP_matrix_for_SigProfiler.txt"
    COSMIC_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/COSMIC_v3.4_SBS_GRCh38.txt"
    OUTPUT_DIR = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/signature_refitting_hnscc_all_cells"
    
    # =========================================================================
    # PARAMETERS
    # =========================================================================
    
    # MUTATION FILTERING
    # Set to 0 to use all cells (no filtering)
    # Set to >0 to filter cells with fewer mutations
    MUTATION_THRESHOLD = 0  # Start with no filtering
    
    # SCREE PLOT SELECTION
    # If True: Use scree plot to find optimal number of signatures
    # If False: Use all 15 HNSCC signatures
    USE_SCREE_PLOT = False  # Recommended: True
    
    # Core signatures (always included)
    CORE_SIGNATURES = ['SBS2', 'SBS3', 'SBS5']
    
    # Candidate order (by mean weight in dataset)
    CANDIDATE_ORDER = ['SBS33', 'SBS7b', 'SBS1', 'SBS4', 'SBS17a', 'SBS17b',
                       'SBS7d', 'SBS7a', 'SBS16', 'SBS18', 'SBS13']
    
    # =========================================================================
    # RUN PIPELINE
    # =========================================================================
    
    results = run_signature_refitting_pipeline(
        mutation_matrix_file=MUTATION_MATRIX_FILE,
        cosmic_file=COSMIC_FILE,
        output_dir=OUTPUT_DIR,
        mutation_threshold=MUTATION_THRESHOLD,
        use_scree_plot=USE_SCREE_PLOT,
        core_signatures=CORE_SIGNATURES,
        candidate_order=CANDIDATE_ORDER
    )
    
    print("\n✓ Signature refitting complete!")
    print(f"Check results in: {OUTPUT_DIR}")
#%%
#!/usr/bin/env python
"""
Post-Processing for SBS2 Signature Analysis

This script performs advanced analysis on signature refitting results:
1. Improved elbow point detection in SBS2 weight distribution
2. Visualization of distribution and threshold
3. Aggregated signature profile for high SBS2 cells
4. Selection of matched SBS2-negative control cells
5. UMAP visualization of selected cells (with proper layering)
6. Export expression matrices for differential network analysis

Updates:
- Fixed UMAP plotting so colored points appear on top
- Improved elbow detection using angle-based method
- Export all inflection candidates with weights

Author: Jake Lehle
Date: 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# COSMIC signature colors
COSMIC_COLORS = {
    'C>A': '#1EBFF0',
    'C>G': '#050708',
    'C>T': '#E62725',
    'T>A': '#CBCACB',
    'T>C': '#A1CE63',
    'T>G': '#EDB6C2'
}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def log(msg, level="INFO"):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {level}: {msg}", flush=True)


def save_figure(fig, output_dir, filename, dpi=300):
    """Save figure to output directory"""
    filepath = Path(output_dir) / filename
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
    log(f"Saved: {filepath}")
    return filepath


def verify_context_alignment(mutation_matrix, cosmic_sigs):
    """
    Verify that mutation contexts are properly aligned between data and COSMIC
    
    This is critical for correct signature plotting!
    """
    log("="*80)
    log("VERIFYING CONTEXT ALIGNMENT")
    log("="*80)
    
    mut_contexts = mutation_matrix.index.tolist()
    cosmic_contexts = cosmic_sigs.index.tolist()
    
    log(f"\nMutation matrix contexts: {len(mut_contexts)}")
    log(f"COSMIC signature contexts: {len(cosmic_contexts)}")
    
    # Show first and last few contexts
    log(f"\nMutation matrix - First 5 contexts:")
    for i, ctx in enumerate(mut_contexts[:5]):
        log(f"  {i}: {ctx}")
    
    log(f"\nMutation matrix - Last 5 contexts:")
    for i, ctx in enumerate(mut_contexts[-5:], len(mut_contexts)-5):
        log(f"  {i}: {ctx}")
    
    log(f"\nCOSMIC - First 5 contexts:")
    for i, ctx in enumerate(cosmic_contexts[:5]):
        log(f"  {i}: {ctx}")
    
    log(f"\nCOSMIC - Last 5 contexts:")
    for i, ctx in enumerate(cosmic_contexts[-5:], len(cosmic_contexts)-5):
        log(f"  {i}: {ctx}")
    
    # Check if contexts match
    if len(mut_contexts) != len(cosmic_contexts):
        log(f"\n⚠️  WARNING: Different number of contexts!", level="WARNING")
        log(f"   Mutation matrix: {len(mut_contexts)} contexts")
        log(f"   COSMIC: {len(cosmic_contexts)} contexts")
        alignment_status = "mismatch"
    else:
        contexts_match = all(m == c for m, c in zip(mut_contexts, cosmic_contexts))
        
        if contexts_match:
            log("\n✓ Contexts are PERFECTLY aligned!", level="INFO")
            alignment_status = "perfect"
        elif set(mut_contexts) == set(cosmic_contexts):
            log("\n⚠️  Contexts are the SAME but in DIFFERENT ORDER!", level="WARNING")
            log("   COSMIC signatures will be reordered to match mutation matrix")
            alignment_status = "reorder_needed"
        else:
            log("\n⚠️  WARNING: Different contexts detected!", level="WARNING")
            
            mut_set = set(mut_contexts)
            cosmic_set = set(cosmic_contexts)
            
            missing_in_cosmic = mut_set - cosmic_set
            missing_in_mut = cosmic_set - mut_set
            
            if missing_in_cosmic:
                log(f"   Contexts in mutation matrix but not in COSMIC: {len(missing_in_cosmic)}")
                if len(missing_in_cosmic) <= 5:
                    for ctx in list(missing_in_cosmic)[:5]:
                        log(f"     - {ctx}")
            
            if missing_in_mut:
                log(f"   Contexts in COSMIC but not in mutation matrix: {len(missing_in_mut)}")
                if len(missing_in_mut) <= 5:
                    for ctx in list(missing_in_mut)[:5]:
                        log(f"     - {ctx}")
            
            alignment_status = "mismatch"
    
    # Verify SBS2 characteristic pattern
    log(f"\nVerifying COSMIC SBS2 characteristic pattern...")
    
    # Reorder if needed
    if alignment_status == "reorder_needed":
        cosmic_aligned = cosmic_sigs.loc[mut_contexts]
    else:
        cosmic_aligned = cosmic_sigs
    
    # Extract mutation types
    mutation_types = [ctx.split('[')[1].split(']')[0] for ctx in mut_contexts]
    
    # Calculate C>T percentage in SBS2
    ct_indices = [i for i, mt in enumerate(mutation_types) if mt == 'C>T']
    sbs2_values = cosmic_aligned['SBS2'].values
    ct_sum = np.sum(sbs2_values[ct_indices])
    
    log(f"  C>T percentage in COSMIC SBS2: {100*ct_sum:.2f}%")
    
    if ct_sum > 0.80:
        log(f"  ✓ CORRECT: SBS2 is dominated by C>T (APOBEC signature)")
    elif ct_sum > 0.60:
        log(f"  ⚠️  WARNING: C>T is only {100*ct_sum:.1f}% (expected >80%)", level="WARNING")
    else:
        log(f"  ✗ ERROR: C>T is only {100*ct_sum:.1f}% - CONTEXT MISMATCH!", level="ERROR")
    
    # Show top 3 contexts for SBS2
    sbs2_sorted = cosmic_aligned['SBS2'].sort_values(ascending=False)
    log(f"\n  Top 3 contexts in SBS2:")
    for i, (ctx, val) in enumerate(sbs2_sorted.head(3).items(), 1):
        log(f"    {i}. {ctx}: {100*val:.2f}%")
    
    log("\n" + "="*80)
    log("CONTEXT VERIFICATION COMPLETE")
    log("="*80 + "\n")
    
    return alignment_status


# =============================================================================
# STEP 1: IMPROVED ELBOW POINT DETECTION
# =============================================================================

def find_elbow_points_angle_method(weights, top_percentile=20, 
                                   smooth_window=51, smooth_poly=3):
    """
    Find elbow points in weight distribution using angle-based method
    
    Instead of looking for high second derivatives (which bias toward
    the steep left side), we look for points where the angle of the
    curve changes most dramatically - true "elbow" points.
    
    Method:
    1. Smooth the curve
    2. Calculate the angle between consecutive line segments
    3. Find points with maximum angle change (sharp turns)
    4. These are true elbow points regardless of position
    
    Parameters:
    -----------
    weights : np.array
        SBS2 weights sorted in descending order
    top_percentile : float
        Only consider top X% of cells
    smooth_window : int
        Window size for Savitzky-Golay smoothing (must be odd)
    smooth_poly : int
        Polynomial order for smoothing
    
    Returns:
    --------
    dict with elbow point analysis results
    """
    log("="*80)
    log("ELBOW POINT DETECTION USING ANGLE-BASED METHOD")
    log("="*80)
    
    # Sort weights in descending order
    weights_sorted = np.sort(weights)[::-1]
    
    # Focus on top percentile
    n_top = int(len(weights_sorted) * top_percentile / 100)
    weights_top = weights_sorted[:n_top]
    
    log(f"\nAnalyzing top {top_percentile}% of cells ({n_top:,} cells)")
    log(f"Weight range: {weights_top.min():.4f} to {weights_top.max():.4f}")
    
    # X-axis (cell rank)
    x = np.arange(len(weights_top))
    
    # Smooth the curve to reduce noise
    if len(weights_top) > smooth_window:
        weights_smooth = savgol_filter(weights_top, smooth_window, smooth_poly)
        log(f"Applied Savitzky-Golay smoothing (window={smooth_window}, poly={smooth_poly})")
    else:
        weights_smooth = gaussian_filter1d(weights_top, sigma=5)
        log(f"Applied Gaussian smoothing (sigma=5)")
    
    # Calculate first derivative (slope)
    dy = np.gradient(weights_smooth)
    
    # Calculate second derivative (curvature)
    d2y = np.gradient(dy)
    
    log("\nCalculated derivatives:")
    log(f"  First derivative (dy/dx) range: {dy.min():.6f} to {dy.max():.6f}")
    log(f"  Second derivative (d²y/dx²) range: {d2y.min():.6f} to {d2y.max():.6f}")
    
    # -------------------------------------------------------------------------
    # NEW METHOD: Calculate angles between line segments
    # -------------------------------------------------------------------------
    log("\nCalculating angles between line segments...")
    
    # For each point, calculate the angle between the line segment before
    # and after that point
    angles = []
    
    window = 10  # Look at segments of this length on each side
    
    for i in range(window, len(weights_smooth) - window):
        # Points before current point
        x1, y1 = i - window, weights_smooth[i - window]
        x2, y2 = i, weights_smooth[i]
        
        # Points after current point
        x3, y3 = i, weights_smooth[i]
        x4, y4 = i + window, weights_smooth[i + window]
        
        # Calculate slopes
        slope1 = (y2 - y1) / (x2 - x1) if (x2 - x1) != 0 else 0
        slope2 = (y4 - y3) / (x4 - x3) if (x4 - x3) != 0 else 0
        
        # Calculate angle between slopes (in degrees)
        # Using arctangent to get angle from slope
        angle1 = np.arctan(slope1)
        angle2 = np.arctan(slope2)
        
        # Angle change at this point
        angle_change = abs(angle2 - angle1) * 180 / np.pi  # Convert to degrees
        
        angles.append({
            'index': i,
            'angle_change': angle_change,
            'weight': weights_top[i],
            'slope_before': slope1,
            'slope_after': slope2,
            'second_derivative': d2y[i]
        })
    
    # Sort by angle change (largest = sharpest elbow)
    angles_sorted = sorted(angles, key=lambda x: x['angle_change'], reverse=True)
    
    log(f"\nFound {len(angles_sorted)} potential elbow points")
    log("\nTop 50 elbow candidates (by angle change):")
    log(f"{'Rank':<6} {'Index':<8} {'Weight':<10} {'Angle':<12} {'d²y/dx²':<12}")
    log("-" * 60)
    
    for i, candidate in enumerate(angles_sorted[:50], 1):
        log(f"{i:<6} {candidate['index']:<8} {candidate['weight']:<10.4f} "
            f"{candidate['angle_change']:<12.2f} {candidate['second_derivative']:<12.6f}")
    
    # -------------------------------------------------------------------------
    # Find specific weight targets (e.g., ~2.0, ~1.5)
    # -------------------------------------------------------------------------
    log(f"\n{'='*70}")
    log("FINDING ELBOW POINTS NEAR TARGET WEIGHTS")
    log(f"{'='*70}")
    
    target_weights = [2.0, 1.5, 1.0, 0.5]
    
    for target in target_weights:
        # Find candidates near this weight
        candidates_near_target = [c for c in angles_sorted 
                                 if abs(c['weight'] - target) < 0.3]
        
        if candidates_near_target:
            best = candidates_near_target[0]
            log(f"\nBest elbow near weight {target}:")
            log(f"  Rank: {angles_sorted.index(best) + 1}")
            log(f"  Index: {best['index']}")
            log(f"  Weight: {best['weight']:.4f}")
            log(f"  Angle change: {best['angle_change']:.2f}°")
        else:
            log(f"\nNo elbow found near weight {target}")
    
    results = {
        'weights_sorted': weights_sorted,
        'weights_top': weights_top,
        'weights_smooth': weights_smooth,
        'first_derivative': dy,
        'second_derivative': d2y,
        'elbow_candidates': angles_sorted,
        'n_top': n_top
    }
    
    log("\n" + "="*80)
    log("ELBOW POINT DETECTION COMPLETE")
    log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 2: VISUALIZE DISTRIBUTION AND THRESHOLD
# =============================================================================

def plot_sbs2_distribution_with_threshold(inflection_results, output_dir):
    """
    Create comprehensive visualization of SBS2 distribution with elbow points
    """
    log("="*80)
    log("VISUALIZING SBS2 DISTRIBUTION AND THRESHOLD")
    log("="*80)
    
    weights_top = inflection_results['weights_top']
    weights_smooth = inflection_results['weights_smooth']
    dy = inflection_results['first_derivative']
    d2y = inflection_results['second_derivative']
    candidates = inflection_results['elbow_candidates']
    
    fig = plt.figure(figsize=(18, 12))
    
    # -------------------------------------------------------------------------
    # Panel 1: Full distribution histogram
    # -------------------------------------------------------------------------
    ax1 = plt.subplot(2, 3, 1)
    
    ax1.hist(inflection_results['weights_sorted'], bins=100, 
             color='steelblue', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('SBS2 Weight', fontsize=11)
    ax1.set_ylabel('Number of Cells', fontsize=11)
    ax1.set_title('Full SBS2 Weight Distribution\n(All Cells)',
                  fontsize=12, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Mark top 20% boundary
    threshold_20pct = inflection_results['weights_sorted'][inflection_results['n_top']]
    ax1.axvline(threshold_20pct, color='red', linestyle='--', linewidth=2,
                label=f'Top 20% ({threshold_20pct:.4f})')
    ax1.legend(fontsize=10)
    
    # -------------------------------------------------------------------------
    # Panel 2: Top 20% distribution with elbow points
    # -------------------------------------------------------------------------
    ax2 = plt.subplot(2, 3, 2)
    
    ax2.hist(weights_top, bins=50, color='coral', alpha=0.7, edgecolor='black')
    ax2.set_xlabel('SBS2 Weight', fontsize=11)
    ax2.set_ylabel('Number of Cells', fontsize=11)
    ax2.set_title('Top 20% SBS2 Weight Distribution\nwith Elbow Points',
                  fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    
    # Mark top 10 elbow points
    colors_elbow = ['red', 'orange', 'gold', 'yellow', 'lime', 
                    'green', 'cyan', 'blue', 'purple', 'magenta']
    
    for i, candidate in enumerate(candidates[:10]):
        color = colors_elbow[i] if i < len(colors_elbow) else 'gray'
        linestyle = '--' if i < 5 else ':'
        linewidth = 2 if i < 5 else 1.5
        alpha = 0.8 if i < 5 else 0.6
        
        ax2.axvline(candidate['weight'], color=color, linestyle=linestyle, 
                   linewidth=linewidth, alpha=alpha,
                   label=f"#{i+1}: {candidate['weight']:.3f}")
    
    if len(candidates) > 0:
        ax2.legend(fontsize=7, ncol=2)
    
    # -------------------------------------------------------------------------
    # Panel 3: Smoothed curve with elbow points marked
    # -------------------------------------------------------------------------
    ax3 = plt.subplot(2, 3, 3)
    
    x = np.arange(len(weights_top))
    
    # Plot smoothed weights
    ax3.plot(x, weights_smooth, 'b-', linewidth=2, label='Smoothed weights')
    ax3.set_xlabel('Cell Rank (Descending)', fontsize=11)
    ax3.set_ylabel('SBS2 Weight', fontsize=11)
    ax3.set_title('Smoothed Weight Curve with Elbow Points\n(Top 20%)',
                  fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Mark top 10 elbow points
    for i, candidate in enumerate(candidates[:10]):
        color = colors_elbow[i] if i < len(colors_elbow) else 'gray'
        marker_size = 150 if i < 5 else 80
        edge_width = 2 if i < 5 else 1
        
        ax3.scatter(candidate['index'], candidate['weight'], 
                   color=color, s=marker_size, zorder=5, 
                   edgecolor='black', linewidth=edge_width,
                   label=f"#{i+1}" if i < 5 else None)
    
    ax3.legend(fontsize=10, loc='upper right')
    
    # -------------------------------------------------------------------------
    # Panel 4: Angle changes (NEW - this shows true elbow strength)
    # -------------------------------------------------------------------------
    ax4 = plt.subplot(2, 3, 4)
    
    # Plot angle changes for all candidates
    indices = [c['index'] for c in candidates]
    angle_changes = [c['angle_change'] for c in candidates]
    
    ax4.scatter(indices, angle_changes, alpha=0.5, s=20, color='purple')
    
    # Highlight top 10
    for i, candidate in enumerate(candidates[:10]):
        color = colors_elbow[i] if i < len(colors_elbow) else 'gray'
        marker_size = 150 if i < 5 else 80
        
        ax4.scatter(candidate['index'], candidate['angle_change'],
                   color=color, s=marker_size, zorder=5,
                   edgecolor='black', linewidth=2)
    
    ax4.set_xlabel('Cell Rank', fontsize=11)
    ax4.set_ylabel('Angle Change (degrees)', fontsize=11)
    ax4.set_title('Elbow Detection: Angle Changes\n(Larger = Sharper Elbow)',
                  fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 5: Second derivative (for comparison with angle method)
    # -------------------------------------------------------------------------
    ax5 = plt.subplot(2, 3, 5)
    
    ax5.plot(x, d2y, 'r-', linewidth=2)
    ax5.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax5.set_xlabel('Cell Rank', fontsize=11)
    ax5.set_ylabel('d²y/dx² (Curvature)', fontsize=11)
    ax5.set_title('Second Derivative\n(For Comparison)',
                  fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    
    # Mark top 10 elbow points
    for i, candidate in enumerate(candidates[:10]):
        color = colors_elbow[i] if i < len(colors_elbow) else 'gray'
        ax5.axvline(candidate['index'], color=color, linestyle=':', 
                   linewidth=1, alpha=0.5)
    
    # -------------------------------------------------------------------------
    # Panel 6: Summary statistics
    # -------------------------------------------------------------------------
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    summary_text = f"""
    SBS2 ELBOW POINT ANALYSIS
    {'='*50}
    
    Total cells analyzed:     {len(inflection_results['weights_sorted']):>10,}
    Top 20% cells:            {len(weights_top):>10,}
    
    Weight range (top 20%):
      Maximum:                {weights_top.max():>10.4f}
      Minimum:                {weights_top.min():>10.4f}
      Range:                  {weights_top.max() - weights_top.min():>10.4f}
    
    {'='*50}
    TOP 10 ELBOW CANDIDATES:
    {'='*50}
    """
    
    for i, candidate in enumerate(candidates[:10], 1):
        summary_text += f"\n    {i:2d}. Rank {candidate['index']:>5,} | "
        summary_text += f"Wt: {candidate['weight']:.4f}"
        summary_text += f"\n        Angle: {candidate['angle_change']:>6.2f}°\n"
    
    summary_text += f"""
    {'='*50}
    RECOMMENDED THRESHOLDS:
      Best elbow: {candidates[0]['weight']:.4f}
      Cells above: {np.sum(inflection_results['weights_sorted'] >= candidates[0]['weight']):,}
    """
    
    # Find elbows near specific weights
    for target in [2.0, 1.5]:
        near_target = [c for c in candidates[:50] 
                      if abs(c['weight'] - target) < 0.3]
        if near_target:
            best = near_target[0]
            rank = candidates.index(best) + 1
            summary_text += f"\n      Near {target}: #{rank} - {best['weight']:.4f}"
    
    ax6.text(0.05, 0.5, summary_text, fontsize=8, family='monospace',
             verticalalignment='center', transform=ax6.transAxes)
    
    plt.suptitle('SBS2 Weight Distribution: Elbow Point Analysis',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    save_figure(fig, output_dir, 'sbs2_elbow_point_analysis.png')
    plt.show()
    
    log("\n" + "="*80)
    log("VISUALIZATION COMPLETE")
    log("="*80 + "\n")
    
    # Save elbow candidates to file
    output_path = Path(output_dir)
    candidates_file = output_path / "sbs2_elbow_candidates_all.txt"
    
    with open(candidates_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SBS2 ELBOW POINT CANDIDATES (ALL)\n")
        f.write("="*80 + "\n\n")
        f.write(f"{'Rank':<6} {'Index':<8} {'Weight':<10} {'Angle°':<12} "
                f"{'d²y/dx²':<12} {'Slope Before':<15} {'Slope After':<15}\n")
        f.write("-" * 90 + "\n")
        
        for i, candidate in enumerate(candidates, 1):
            f.write(f"{i:<6} {candidate['index']:<8} {candidate['weight']:<10.4f} "
                   f"{candidate['angle_change']:<12.2f} "
                   f"{candidate['second_derivative']:<12.6f} "
                   f"{candidate['slope_before']:<15.6f} "
                   f"{candidate['slope_after']:<15.6f}\n")
    
    log(f"Saved all elbow candidates: {candidates_file}")


# =============================================================================
# STEP 3: SELECT HIGH SBS2 CELLS AND AGGREGATE SIGNATURE
# =============================================================================

def select_high_sbs2_cells(weights_df, mutation_matrix, threshold, 
                           cosmic_sigs, output_dir):
    """
    Select cells above SBS2 threshold and aggregate their mutation profile
    """
    log("="*80)
    log("SELECTING HIGH SBS2 CELLS")
    log("="*80)
    
    # Get SBS2 weights
    sbs2_weights = weights_df.loc['SBS2']
    
    # Select cells above threshold
    high_sbs2_cells = sbs2_weights[sbs2_weights >= threshold].index.tolist()
    
    log(f"\nThreshold: {threshold:.4f}")
    log(f"Cells above threshold: {len(high_sbs2_cells):,}")
    log(f"Percentage of total: {100 * len(high_sbs2_cells) / len(sbs2_weights):.2f}%")
    
    # Aggregate mutations from high SBS2 cells
    high_sbs2_mutations = mutation_matrix[high_sbs2_cells]
    aggregated_profile = high_sbs2_mutations.sum(axis=1)
    
    log(f"\nAggregated mutation profile:")
    log(f"  Total mutations: {int(aggregated_profile.sum()):,}")
    log(f"  Mean per context: {aggregated_profile.mean():.2f}")
    
    # Calculate statistics
    mean_weight = sbs2_weights[high_sbs2_cells].mean()
    median_weight = sbs2_weights[high_sbs2_cells].median()
    min_weight = sbs2_weights[high_sbs2_cells].min()
    max_weight = sbs2_weights[high_sbs2_cells].max()
    
    log(f"\nSBS2 weight statistics for selected cells:")
    log(f"  Mean:   {mean_weight:.4f}")
    log(f"  Median: {median_weight:.4f}")
    log(f"  Range:  {min_weight:.4f} - {max_weight:.4f}")
    
    # Normalize aggregated profile for COSMIC comparison
    aggregated_normalized = aggregated_profile / aggregated_profile.sum()
    
    # Compare to COSMIC SBS2
    cosmic_sbs2 = cosmic_sigs['SBS2'].values
    
    # Calculate similarity
    from scipy.spatial.distance import cosine, euclidean
    cosine_sim = 1 - cosine(aggregated_normalized, cosmic_sbs2)
    euclidean_dist = euclidean(aggregated_normalized, cosmic_sbs2)
    
    # Pearson correlation
    pearson_r, _ = pearsonr(aggregated_normalized, cosmic_sbs2)
    
    log(f"\nSimilarity to COSMIC SBS2:")
    log(f"  Cosine similarity: {cosine_sim:.4f}")
    log(f"  Pearson correlation: {pearson_r:.4f}")
    log(f"  Euclidean distance: {euclidean_dist:.4f}")
    
    results = {
        'high_sbs2_cells': high_sbs2_cells,
        'aggregated_profile': aggregated_profile,
        'aggregated_normalized': aggregated_normalized,
        'n_cells': len(high_sbs2_cells),
        'total_mutations': int(aggregated_profile.sum()),
        'cosine_similarity': cosine_sim,
        'pearson_correlation': pearson_r
    }
    
    log("\n" + "="*80)
    log("HIGH SBS2 CELL SELECTION COMPLETE")
    log("="*80 + "\n")
    
    return results


def plot_aggregated_signature(aggregated_profile, cosmic_sigs, output_dir,
                              title_prefix="High SBS2 Cells"):
    """
    Plot aggregated signature in COSMIC style with proper context alignment
    """
    log(f"Plotting aggregated signature profile for {title_prefix}...")
    
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    
    # Get contexts from aggregated profile
    contexts = aggregated_profile.index.tolist()
    
    # CRITICAL: Ensure COSMIC SBS2 is aligned to the same contexts
    contexts_match = all(c == a for c, a in zip(cosmic_sigs.index, aggregated_profile.index))
    
    if not contexts_match:
        log(f"  ⚠️  COSMIC contexts don't match mutation matrix, reordering...")
        cosmic_sigs_aligned = cosmic_sigs.loc[aggregated_profile.index]
    else:
        cosmic_sigs_aligned = cosmic_sigs
    
    # Extract mutation types from context strings
    mutation_types = [ctx.split('[')[1].split(']')[0] for ctx in contexts]
    
    # Normalize aggregated profile
    total = aggregated_profile.sum()
    values = (aggregated_profile.values / total * 100) if total > 0 else aggregated_profile.values
    
    # Get COSMIC SBS2 values
    cosmic_sbs2_values = cosmic_sigs_aligned['SBS2'].values * 100
    
    # Panel 1: Aggregated profile
    ax1 = axes[0]
    x = np.arange(len(contexts))
    colors = [COSMIC_COLORS[mut] for mut in mutation_types]
    ax1.bar(x, values, color=colors, width=0.8, edgecolor='none')
    
    # Separation lines
    for i in range(1, 6):
        ax1.axvline(i*16 - 0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    ax1.set_ylabel('% of mutations', fontsize=10)
    ax1.set_title(f'{title_prefix}: Aggregated Mutation Profile ({int(total):,} mutations)',
                  fontsize=12, fontweight='bold')
    
    # X-axis
    tick_positions = [8, 24, 40, 56, 72, 88]
    tick_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels, fontsize=10)
    ax1.set_xlim(-0.5, 95.5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Panel 2: COSMIC SBS2 reference
    ax2 = axes[1]
    ax2.bar(x, cosmic_sbs2_values, color=colors, width=0.8, edgecolor='none')
    
    # Separation lines
    for i in range(1, 6):
        ax2.axvline(i*16 - 0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    ax2.set_ylabel('% of mutations', fontsize=10)
    ax2.set_title(f'COSMIC SBS2 Reference (APOBEC)',
                  fontsize=12, fontweight='bold')
    
    # X-axis
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(tick_labels, fontsize=10)
    ax2.set_xlim(-0.5, 95.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    save_figure(fig, output_dir, f'{title_prefix.lower().replace(" ", "_")}_aggregated_signature.png')
    plt.show()
    
    log(f"  ✓ Plot saved successfully")


# =============================================================================
# STEP 4: SELECT MATCHED SBS2-NEGATIVE CONTROL CELLS
# =============================================================================

def select_matched_control_cells(weights_df, mutation_matrix, high_sbs2_cells,
                                 n_controls=None):
    """
    Select SBS2-negative control cells matched for:
    1. Similar total mutation count
    2. Similar average mutations per context
    3. Low or anti-correlated SBS2 signature
    4. Anti-correlated signature profiles (different patterns)
    """
    log("="*80)
    log("SELECTING MATCHED SBS2-NEGATIVE CONTROL CELLS")
    log("="*80)
    
    if n_controls is None:
        n_controls = len(high_sbs2_cells)
    
    # Calculate mutation statistics for high SBS2 cells
    high_sbs2_muts = mutation_matrix[high_sbs2_cells]
    high_sbs2_total = high_sbs2_muts.sum(axis=0).mean()
    high_sbs2_per_context = high_sbs2_total / mutation_matrix.shape[0]
    
    log(f"\nHigh SBS2 cell statistics:")
    log(f"  Number of cells: {len(high_sbs2_cells):,}")
    log(f"  Mean total mutations: {high_sbs2_total:.2f}")
    log(f"  Mean mutations per context: {high_sbs2_per_context:.4f}")
    
    # Get all other cells
    all_cells = set(mutation_matrix.columns)
    candidate_cells = list(all_cells - set(high_sbs2_cells))
    
    log(f"\nCandidate control cells: {len(candidate_cells):,}")
    
    # Calculate statistics for all candidate cells
    candidate_muts = mutation_matrix[candidate_cells]
    candidate_totals = candidate_muts.sum(axis=0)
    candidate_per_context = candidate_totals / mutation_matrix.shape[0]
    
    # Get SBS2 weights for candidates
    sbs2_weights = weights_df.loc['SBS2', candidate_cells]
    
    # Calculate signature profile correlations
    log("\nCalculating signature profile correlations...")
    high_sbs2_avg_profile = weights_df[high_sbs2_cells].mean(axis=1)
    
    profile_correlations = []
    for cell in candidate_cells:
        cell_profile = weights_df[cell]
        corr, _ = pearsonr(high_sbs2_avg_profile, cell_profile)
        profile_correlations.append(corr)
    
    profile_correlations = np.array(profile_correlations)
    
    # Scoring system
    log("\nScoring candidate control cells...")
    
    # 1. Similar mutation count
    mut_diff = np.abs(candidate_totals.values - high_sbs2_total)
    mut_score = 1 - (mut_diff / mut_diff.max())
    
    # 2. Low SBS2 weight
    sbs2_score = 1 - (sbs2_weights.values / sbs2_weights.max())
    
    # 3. Anti-correlated signature profile
    profile_score = (1 - profile_correlations) / 2
    
    # Combined score
    combined_score = (mut_score + sbs2_score + profile_score) / 3
    
    # Create DataFrame
    candidate_df = pd.DataFrame({
        'cell_id': candidate_cells,
        'total_mutations': candidate_totals.values,
        'mutations_per_context': candidate_per_context.values,
        'sbs2_weight': sbs2_weights.values,
        'profile_correlation': profile_correlations,
        'mutation_score': mut_score,
        'sbs2_score': sbs2_score,
        'profile_score': profile_score,
        'combined_score': combined_score
    })
    
    # Sort by combined score
    candidate_df = candidate_df.sort_values('combined_score', ascending=False)
    
    # Select top N controls
    control_cells = candidate_df.head(n_controls)['cell_id'].tolist()
    
    log(f"\nSelected {len(control_cells):,} control cells")
    log(f"\nControl cell statistics:")
    log(f"  Mean total mutations: {candidate_df.head(n_controls)['total_mutations'].mean():.2f}")
    log(f"  Mean mutations per context: {candidate_df.head(n_controls)['mutations_per_context'].mean():.4f}")
    log(f"  Mean SBS2 weight: {candidate_df.head(n_controls)['sbs2_weight'].mean():.4f}")
    log(f"  Mean profile correlation: {candidate_df.head(n_controls)['profile_correlation'].mean():.4f}")
    
    log("\n" + "="*80)
    log("CONTROL CELL SELECTION COMPLETE")
    log("="*80 + "\n")
    
    return control_cells, candidate_df


# =============================================================================
# STEP 5: UMAP VISUALIZATION (FIXED)
# =============================================================================

def plot_umap_with_sbs2_groups(adata, high_sbs2_cells, control_cells, 
                                output_dir, size=20):
    """
    Plot cells on UMAP colored by SBS2 group
    
    FIXED: Plots colored points on top of gray points
    """
    log("="*80)
    log("PLOTTING CELLS ON UMAP")
    log("="*80)
    
    try:
        import scanpy as sc
    except ImportError:
        log("⚠️  Scanpy not installed, skipping UMAP plot", level="WARNING")
        return
    
    # Create SBS2 group annotation
    adata.obs['SBS2_group'] = 'Other'
    adata.obs.loc[adata.obs.index.isin(high_sbs2_cells), 'SBS2_group'] = 'High SBS2'
    adata.obs.loc[adata.obs.index.isin(control_cells), 'SBS2_group'] = 'Low SBS2 (Control)'
    
    # Convert to categorical with specific order (Other first so it plots first)
    adata.obs['SBS2_group'] = pd.Categorical(
        adata.obs['SBS2_group'],
        categories=['Other', 'High SBS2', 'Low SBS2 (Control)'],
        ordered=True
    )
    
    log(f"\nCell group counts:")
    log(f"  High SBS2: {sum(adata.obs['SBS2_group'] == 'High SBS2'):,}")
    log(f"  Low SBS2 (Control): {sum(adata.obs['SBS2_group'] == 'Low SBS2 (Control)'):,}")
    log(f"  Other: {sum(adata.obs['SBS2_group'] == 'Other'):,}")
    
    # Set up color palette
    custom_palette = {
        'Other': 'lightgray',
        'High SBS2': 'red',
        'Low SBS2 (Control)': 'green'
    }
    
    # CRITICAL FIX: Plot manually to control layering
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    
    # Plot in correct order: Other (gray) first, then colored groups on top
    # This ensures colored points are visible
    
    # 1. Plot Other cells (background)
    other_mask = adata.obs['SBS2_group'] == 'Other'
    ax.scatter(umap_coords[other_mask, 0], umap_coords[other_mask, 1],
               c='lightgray', s=size, alpha=0.5, label='Other', rasterized=True)
    
    # 2. Plot control cells (on top)
    control_mask = adata.obs['SBS2_group'] == 'Low SBS2 (Control)'
    ax.scatter(umap_coords[control_mask, 0], umap_coords[control_mask, 1],
               c='green', s=size, alpha=0.8, label='Low SBS2 (Control)', 
               edgecolors='black', linewidths=0.5, rasterized=True)
    
    # 3. Plot high SBS2 cells (on top)
    high_mask = adata.obs['SBS2_group'] == 'High SBS2'
    ax.scatter(umap_coords[high_mask, 0], umap_coords[high_mask, 1],
               c='red', s=size, alpha=0.8, label='High SBS2',
               edgecolors='black', linewidths=0.5, rasterized=True)
    
    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title('UMAP: High SBS2 vs Control Cells',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=10, loc='best', frameon=True, fancybox=True, shadow=True)
    ax.set_aspect('equal')
    
    # Remove ticks
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.tight_layout()
    save_figure(fig, output_dir, 'umap_sbs2_groups.png')
    plt.show()
    
    log("\n" + "="*80)
    log("UMAP VISUALIZATION COMPLETE")
    log("="*80 + "\n")


# =============================================================================
# STEP 6: EXPORT EXPRESSION MATRICES
# =============================================================================

def export_expression_matrices(adata, high_sbs2_cells, control_cells, 
                                output_dir):
    """
    Export expression matrices for differential network analysis
    """
    log("="*80)
    log("EXPORTING EXPRESSION MATRICES")
    log("="*80)
    
    try:
        import scipy.sparse
    except ImportError:
        log("⚠️  Scipy not installed, cannot export matrices", level="ERROR")
        return
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Ensure gene IDs are set
    if 'gene_ids' in adata.var.columns:
        adata.var.index = adata.var['gene_ids']
    
    # -------------------------------------------------------------------------
    # High SBS2 cells
    # -------------------------------------------------------------------------
    log(f"\nExporting high SBS2 cells ({len(high_sbs2_cells):,})...")
    
    high_sbs2_adata = adata[adata.obs.index.isin(high_sbs2_cells)]
    
    # Get gene symbols
    if 'gene_ids' in high_sbs2_adata.var.columns:
        gene_names = high_sbs2_adata.var['gene_ids'].values
    else:
        gene_names = high_sbs2_adata.var.index
    
    # Convert to DataFrame
    if scipy.sparse.issparse(high_sbs2_adata.X):
        expression_df_high = pd.DataFrame.sparse.from_spmatrix(
            scipy.sparse.csr_matrix(high_sbs2_adata.X),
            index=high_sbs2_adata.obs.index,
            columns=gene_names
        ).T
    else:
        expression_df_high = pd.DataFrame(
            high_sbs2_adata.X.T,
            index=gene_names,
            columns=high_sbs2_adata.obs.index
        )
    
    # Save
    high_file = output_path / "High_SBS2_Cells_Expression.tsv"
    expression_df_high.to_csv(high_file, sep="\t", header=True, index=True, chunksize=10000)
    log(f"Saved: {high_file}")
    log(f"  Shape: {expression_df_high.shape[0]:,} genes × {expression_df_high.shape[1]:,} cells")
    
    # -------------------------------------------------------------------------
    # Control cells
    # -------------------------------------------------------------------------
    log(f"\nExporting control cells ({len(control_cells):,})...")
    
    control_adata = adata[adata.obs.index.isin(control_cells)]
    
    # Get gene symbols
    if 'gene_ids' in control_adata.var.columns:
        gene_names = control_adata.var['gene_ids'].values
    else:
        gene_names = control_adata.var.index
    
    # Convert to DataFrame
    if scipy.sparse.issparse(control_adata.X):
        expression_df_control = pd.DataFrame.sparse.from_spmatrix(
            scipy.sparse.csr_matrix(control_adata.X),
            index=control_adata.obs.index,
            columns=gene_names
        ).T
    else:
        expression_df_control = pd.DataFrame(
            control_adata.X.T,
            index=gene_names,
            columns=control_adata.obs.index
        )
    
    # Save
    control_file = output_path / "Control_Cells_Expression.tsv"
    expression_df_control.to_csv(control_file, sep="\t", header=True, index=True, chunksize=10000)
    log(f"Saved: {control_file}")
    log(f"  Shape: {expression_df_control.shape[0]:,} genes × {expression_df_control.shape[1]:,} cells")
    
    # -------------------------------------------------------------------------
    # Save cell lists
    # -------------------------------------------------------------------------
    cell_lists_file = output_path / "cell_lists.txt"
    with open(cell_lists_file, 'w') as f:
        f.write("HIGH SBS2 CELLS\n")
        f.write("="*50 + "\n")
        for cell in high_sbs2_cells:
            f.write(f"{cell}\n")
        
        f.write("\n\nCONTROL CELLS (SBS2-NEGATIVE)\n")
        f.write("="*50 + "\n")
        for cell in control_cells:
            f.write(f"{cell}\n")
    
    log(f"\nSaved cell lists: {cell_lists_file}")
    
    log("\n" + "="*80)
    log("EXPRESSION MATRICES EXPORTED")
    log("="*80 + "\n")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_sbs2_postprocessing_pipeline(weights_file, mutation_matrix_file,
                                     cosmic_file, output_dir, adata=None,
                                     top_percentile=20, elbow_rank=0):
    """
    Complete SBS2 post-processing pipeline
    
    Parameters:
    -----------
    weights_file : str
        Path to signature weights from refitting
    mutation_matrix_file : str
        Path to original mutation matrix
    cosmic_file : str
        Path to COSMIC signatures
    output_dir : str
        Output directory
    adata : AnnData, optional
        Single-cell expression data for UMAP plotting
    top_percentile : float
        Top X% to analyze for elbow points (default: 20)
    elbow_rank : int
        Which elbow point to use as threshold (0 = best, 1 = second best, etc.)
    
    Returns:
    --------
    dict with all results
    """
    
    print("\n" + "="*80)
    print("SBS2 POST-PROCESSING PIPELINE")
    print("Elbow Point Detection & Matched Control Selection")
    print("="*80)
    print(f"\nWeights: {weights_file}")
    print(f"Mutations: {mutation_matrix_file}")
    print(f"COSMIC: {cosmic_file}")
    print(f"Output: {output_dir}")
    print("="*80 + "\n")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------
    log("Loading data...")
    
    weights_df = pd.read_csv(weights_file, sep='\t', index_col=0)
    mutation_matrix = pd.read_csv(mutation_matrix_file, sep='\t', index_col=0)
    cosmic_sigs = pd.read_csv(cosmic_file, sep='\t', index_col=0)
    
    log(f"Weights: {weights_df.shape}")
    log(f"Mutations: {mutation_matrix.shape}")
    log(f"COSMIC: {cosmic_sigs.shape}")
    
    # -------------------------------------------------------------------------
    # Verify context alignment
    # -------------------------------------------------------------------------
    alignment_status = verify_context_alignment(mutation_matrix, cosmic_sigs)
    
    # Reorder COSMIC if needed
    if alignment_status == "reorder_needed":
        log("Reordering COSMIC signatures to match mutation matrix context order...")
        cosmic_sigs = cosmic_sigs.loc[mutation_matrix.index]
        log("✓ COSMIC signatures reordered")
    elif alignment_status == "mismatch":
        log("⚠️  WARNING: Context mismatch detected. Results may be incorrect!", level="WARNING")
    
    # -------------------------------------------------------------------------
    # Step 1: Find elbow points using angle method
    # -------------------------------------------------------------------------
    sbs2_weights = weights_df.loc['SBS2'].values
    
    elbow_results = find_elbow_points_angle_method(
        sbs2_weights,
        top_percentile=top_percentile
    )
    
    # -------------------------------------------------------------------------
    # Step 2: Visualize distribution
    # -------------------------------------------------------------------------
    plot_sbs2_distribution_with_threshold(elbow_results, output_dir)
    
    # -------------------------------------------------------------------------
    # Step 3: Select threshold and high SBS2 cells
    # -------------------------------------------------------------------------
    candidates = elbow_results['elbow_candidates']
    
    if len(candidates) == 0:
        log("⚠️  No elbow points found, using median of top 20%", level="WARNING")
        threshold = np.median(elbow_results['weights_top'])
    else:
        # Use specified elbow point
        elbow_idx = min(elbow_rank, len(candidates) - 1)
        threshold = candidates[elbow_idx]['weight']
        
        log(f"\nUsing elbow point #{elbow_idx + 1} as threshold")
        log(f"Threshold: {threshold:.4f}")
        log(f"Angle change: {candidates[elbow_idx]['angle_change']:.2f}°")
    
    high_sbs2_results = select_high_sbs2_cells(
        weights_df,
        mutation_matrix,
        threshold,
        cosmic_sigs,
        output_dir
    )
    
    # -------------------------------------------------------------------------
    # Step 4: Plot aggregated signature
    # -------------------------------------------------------------------------
    plot_aggregated_signature(
        high_sbs2_results['aggregated_profile'],
        cosmic_sigs,
        output_dir,
        title_prefix="High SBS2 Cells"
    )
    
    # -------------------------------------------------------------------------
    # Step 5: Select matched controls
    # -------------------------------------------------------------------------
    control_cells, candidate_df = select_matched_control_cells(
        weights_df,
        mutation_matrix,
        high_sbs2_results['high_sbs2_cells']
    )
    
    # Save candidate scoring
    candidate_file = output_path / "control_cell_candidates_scored.txt"
    candidate_df.to_csv(candidate_file, sep='\t', index=False, float_format='%.4f')
    log(f"Saved candidate scoring: {candidate_file}")
    
    # Plot control signature
    control_mutations = mutation_matrix[control_cells]
    control_profile = control_mutations.sum(axis=1)
    
    plot_aggregated_signature(
        control_profile,
        cosmic_sigs,
        output_dir,
        title_prefix="Control Cells"
    )
    
    # -------------------------------------------------------------------------
    # Step 6: UMAP visualization (if adata provided)
    # -------------------------------------------------------------------------
    if adata is not None:
        plot_umap_with_sbs2_groups(
            adata,
            high_sbs2_results['high_sbs2_cells'],
            control_cells,
            output_dir
        )
        
        # -------------------------------------------------------------------------
        # Step 7: Export expression matrices
        # -------------------------------------------------------------------------
        export_expression_matrices(
            adata,
            high_sbs2_results['high_sbs2_cells'],
            control_cells,
            output_dir
        )
    else:
        log("\n⚠️  No AnnData object provided, skipping UMAP and expression export",
            level="WARNING")
    
    # -------------------------------------------------------------------------
    # Save summary
    # -------------------------------------------------------------------------
    summary_file = output_path / "sbs2_analysis_summary.txt"
    
    with open(summary_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SBS2 POST-PROCESSING ANALYSIS SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        f.write("THRESHOLD SELECTION\n")
        f.write("-"*80 + "\n")
        f.write(f"Method: Angle-based elbow point detection\n")
        f.write(f"Top percentile analyzed: {top_percentile}%\n")
        f.write(f"Selected threshold: {threshold:.4f}\n")
        if len(candidates) > 0:
            f.write(f"Angle change: {candidates[elbow_idx]['angle_change']:.2f}°\n")
        f.write("\n")
        
        f.write("HIGH SBS2 CELLS\n")
        f.write("-"*80 + "\n")
        f.write(f"Number of cells: {high_sbs2_results['n_cells']:,}\n")
        f.write(f"Total mutations: {high_sbs2_results['total_mutations']:,}\n")
        f.write(f"Cosine similarity to COSMIC SBS2: {high_sbs2_results['cosine_similarity']:.4f}\n")
        f.write(f"Pearson correlation to COSMIC SBS2: {high_sbs2_results['pearson_correlation']:.4f}\n\n")
        
        f.write("CONTROL CELLS\n")
        f.write("-"*80 + "\n")
        f.write(f"Number of cells: {len(control_cells):,}\n")
        f.write(f"Selection criteria:\n")
        f.write(f"  - Similar mutation count\n")
        f.write(f"  - Low SBS2 weight\n")
        f.write(f"  - Anti-correlated signature profile\n\n")
        
        if adata is not None:
            f.write("EXPRESSION MATRICES\n")
            f.write("-"*80 + "\n")
            f.write(f"Exported for differential network analysis\n")
            f.write(f"  High SBS2: {high_sbs2_results['n_cells']:,} cells\n")
            f.write(f"  Controls: {len(control_cells):,} cells\n")
    
    log(f"Saved summary: {summary_file}")
    
    # -------------------------------------------------------------------------
    # Final summary
    # -------------------------------------------------------------------------
    print("\n" + "="*80)
    print("PIPELINE COMPLETE!")
    print("="*80)
    
    print(f"\nThreshold: {threshold:.4f}")
    print(f"High SBS2 cells: {high_sbs2_results['n_cells']:,}")
    print(f"Control cells: {len(control_cells):,}")
    print(f"Similarity to COSMIC SBS2: {high_sbs2_results['cosine_similarity']:.4f}")
    
    print(f"\nAll results saved to: {output_dir}")
    print(f"\nCheck 'sbs2_elbow_candidates_all.txt' for full list of elbow points")
    print("="*80 + "\n")
    
    return {
        'elbow_results': elbow_results,
        'threshold': threshold,
        'high_sbs2_cells': high_sbs2_results['high_sbs2_cells'],
        'control_cells': control_cells,
        'high_sbs2_profile': high_sbs2_results['aggregated_profile'],
        'control_profile': control_profile
    }
#%%

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    
    # =========================================================================
    # CONFIGURATION
    # =========================================================================
    
    # Paths to signature refitting results
    WEIGHTS_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/signature_refitting_hnscc_all_cells/signature_weights_per_cell.txt"
    MUTATION_MATRIX_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/SNP_matrix_for_SigProfiler.txt"
    COSMIC_FILE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/COSMIC_v3.4_SBS_GRCh38.txt"
    OUTPUT_DIR = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/sbs2_postprocessing_all_cells"
    
    # Optional: Load AnnData for UMAP and expression export
    ADATA_FILE = None  # Set to your .h5ad file path if available
    
    # Parameters
    TOP_PERCENTILE = 20  # Analyze top 20% for elbow points
    ELBOW_RANK = 270  # Use best elbow point (0=best, 1=second best, etc.)
    
    # =========================================================================
    # LOAD ADATA (OPTIONAL)
    # =========================================================================
    
    adata = None
    
    if ADATA_FILE is not None:
        try:
            import scanpy as sc
            log(f"Loading AnnData from: {ADATA_FILE}")
            adata = sc.read_h5ad(ADATA_FILE)
            log(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
        except Exception as e:
            log(f"⚠️  Could not load AnnData: {e}", level="WARNING")
            log("Continuing without UMAP visualization")
            adata = None
    
    # =========================================================================
    # RUN PIPELINE
    # =========================================================================
    
    results = run_sbs2_postprocessing_pipeline(
        weights_file=WEIGHTS_FILE,
        mutation_matrix_file=MUTATION_MATRIX_FILE,
        cosmic_file=COSMIC_FILE,
        output_dir=OUTPUT_DIR,
        adata=adata_pp,
        top_percentile=TOP_PERCENTILE,
        elbow_rank=ELBOW_RANK
    )
    
    print("\n✓ SBS2 post-processing complete!")
    print(f"Check results in: {OUTPUT_DIR}")
#%%
adata_pp
