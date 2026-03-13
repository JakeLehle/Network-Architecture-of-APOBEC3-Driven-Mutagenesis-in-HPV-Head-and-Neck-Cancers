#!/usr/bin/env Rscript
# ============================================================
# Differential Network Analysis: Identifying A3 Cofactors
# ============================================================
# 
# Rationale: Tumors with similar A3A+A3B expression show different
# SBS2 mutation levels, suggesting unknown cofactors regulate A3 activity.
# This analysis identifies genes with altered correlation patterns
# between high-SBS2 and low-SBS2 tumors (matched for A3 expression).
#
# ============================================================

# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
library(ggraph)
library(viridis)
library(parallel)

# Set working directory
setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

# ============================================================
# ANALYSIS PARAMETERS - Modify these as needed
# ============================================================
PARAMS <- list(
  # Sample selection
  a3_percentile_threshold = 0.50,    # Keep top 50% of A3 expressors
  sbs2_high_percentile = 0.80,       # Top 20% = high SBS2
  sbs2_low_percentile = 0.20,        # Bottom 20% = low SBS2
  
  # Gene filtering
  min_expression = 1,                 # Minimum mean FPKM-UQ across samples
  min_variance_percentile = 0.10,     # Remove bottom 10% low-variance genes
  max_zero_fraction = 0.50,           # Remove genes with >50% zeros
  protein_coding_only = TRUE,         # Keep only protein-coding genes
  

  # Network construction
  correlation_method = "spearman",    # "pearson" or "spearman"
  correlation_threshold = 0.75,       # |r| threshold for edge creation
  
  # Differential analysis
  n_top_differential_genes = 200,     # Number of top DC genes to highlight
  
  # Visualization
  network_layout = "fr",              # Fruchterman-Reingold layout
  show_full_network = TRUE,           # TRUE = show all genes, FALSE = subnetwork only
  subnetwork_depth = 1,               # How many hops from diff genes (only used if show_full_network = FALSE)
  
  # Visualization appearance
  edge_color = "gray30",              # Dark gray for edges
  edge_width = 0.3,                   # Edge line width (increase for thicker edges)
  edge_alpha_range = c(0.2, 0.6),     # Alpha range for edges based on weight
  nonsig_node_color = "#E0E0E0",      # Light gray for non-significant genes
  nonsig_node_alpha = 0.4,            # Alpha for non-significant nodes
  sig_node_color = "#E41A1C",         # Bright red for significant/differential genes
  sig_node_alpha = 0.9,               # Alpha for significant nodes
  label_all_sig_genes = TRUE,         # TRUE = label all differential genes, FALSE = top 20% only
  label_size = 2.5                    # Font size for gene labels
)

cat("============================================================\n")
cat("Differential Network Analysis: Identifying A3 Cofactors\n")
cat("============================================================\n\n")

cat("Analysis parameters:\n")
print(PARAMS)
cat("\n")

# ============================================================
# PHASE 1: Sample Selection & Stratification
# ============================================================
cat("============================================================\n")
cat("PHASE 1: Sample Selection & Stratification\n")
cat("============================================================\n\n")

# Step 1.1: Load and filter data
cat("Step 1.1: Loading matched expression and mutation data\n")
cat("------------------------------------------------------------\n")

expression_file <- "A3_Expression_FPKM_UQ_Matched.tsv"
mutation_file <- "Mutation_Signatures_Matched.tsv"

a3_expression <- fread(expression_file)
mutation_data <- fread(mutation_file)

cat("A3 expression file:", nrow(a3_expression), "samples\n")
cat("Mutation file:", nrow(mutation_data), "samples\n")

# Filter for HNSC tumors only
hnsc_expression <- a3_expression %>%
  filter((Project_ID == "TCGA-HNSC" | Project_ID == "TCGA-HNSCC") & 
         Tissue_Type == "Tumor")

cat("HNSC tumor samples in expression:", nrow(hnsc_expression), "\n")

# Get matching mutation data
hnsc_mutation <- mutation_data %>%
  filter(TCGA_Gene_Expression_Entity_ID %in% hnsc_expression$Entity_ID)

cat("HNSC tumor samples in mutation:", nrow(hnsc_mutation), "\n\n")

# Step 1.2: Calculate A3A + A3B sum
cat("Step 1.2: Calculating A3A + A3B expression sum\n")
cat("------------------------------------------------------------\n")

hnsc_expression <- hnsc_expression %>%
  mutate(A3_sum = APOBEC3A + APOBEC3B)

cat("A3_sum distribution:\n")
cat("  Min:", round(min(hnsc_expression$A3_sum), 4), "\n")
cat("  Q1:", round(quantile(hnsc_expression$A3_sum, 0.25), 4), "\n")
cat("  Median:", round(median(hnsc_expression$A3_sum), 4), "\n")
cat("  Q3:", round(quantile(hnsc_expression$A3_sum, 0.75), 4), "\n")
cat("  Max:", round(max(hnsc_expression$A3_sum), 4), "\n\n")

# Step 1.3: Apply A3 expression threshold (keep top 50%)
cat("Step 1.3: Filtering for high A3 expressors (top", 
    (1 - PARAMS$a3_percentile_threshold) * 100, "%)\n")
cat("------------------------------------------------------------\n")

a3_threshold <- quantile(hnsc_expression$A3_sum, PARAMS$a3_percentile_threshold)
cat("A3_sum threshold (", PARAMS$a3_percentile_threshold * 100, 
    "th percentile):", round(a3_threshold, 4), "\n")

high_a3_samples <- hnsc_expression %>%
  filter(A3_sum >= a3_threshold)

cat("Samples before A3 filter:", nrow(hnsc_expression), "\n")
cat("Samples after A3 filter:", nrow(high_a3_samples), "\n\n")

# Step 1.4: Stratify by SBS2
cat("Step 1.4: Stratifying by SBS2 levels\n")
cat("------------------------------------------------------------\n")

# Merge with mutation data to get SBS2 values
high_a3_with_sbs2 <- high_a3_samples %>%
  inner_join(
    hnsc_mutation %>% select(TCGA_Gene_Expression_Entity_ID, SBS2),
    by = c("Entity_ID" = "TCGA_Gene_Expression_Entity_ID")
  )

cat("Samples with both A3 expression and SBS2 data:", nrow(high_a3_with_sbs2), "\n")

# Calculate SBS2 thresholds
sbs2_high_threshold <- quantile(high_a3_with_sbs2$SBS2, PARAMS$sbs2_high_percentile)
sbs2_low_threshold <- quantile(high_a3_with_sbs2$SBS2, PARAMS$sbs2_low_percentile)

cat("SBS2 high threshold (", PARAMS$sbs2_high_percentile * 100, 
    "th percentile):", round(sbs2_high_threshold, 6), "\n")
cat("SBS2 low threshold (", PARAMS$sbs2_low_percentile * 100, 
    "th percentile):", round(sbs2_low_threshold, 6), "\n")

# Create high and low SBS2 groups
high_sbs2_samples <- high_a3_with_sbs2 %>%
  filter(SBS2 >= sbs2_high_threshold)

low_sbs2_samples <- high_a3_with_sbs2 %>%
  filter(SBS2 <= sbs2_low_threshold)

cat("\nHigh SBS2 group:", nrow(high_sbs2_samples), "samples\n")
cat("  SBS2 range:", round(min(high_sbs2_samples$SBS2), 6), "-", 
    round(max(high_sbs2_samples$SBS2), 6), "\n")
cat("  A3_sum range:", round(min(high_sbs2_samples$A3_sum), 4), "-", 
    round(max(high_sbs2_samples$A3_sum), 4), "\n")

cat("\nLow SBS2 group:", nrow(low_sbs2_samples), "samples\n")
cat("  SBS2 range:", round(min(low_sbs2_samples$SBS2), 6), "-", 
    round(max(low_sbs2_samples$SBS2), 6), "\n")
cat("  A3_sum range:", round(min(low_sbs2_samples$A3_sum), 4), "-", 
    round(max(low_sbs2_samples$A3_sum), 4), "\n")

# Verify A3 expression is similar between groups
cat("\nVerifying A3 expression similarity between groups:\n")
cat("  High SBS2 group - Mean A3_sum:", round(mean(high_sbs2_samples$A3_sum), 4), "\n")
cat("  Low SBS2 group - Mean A3_sum:", round(mean(low_sbs2_samples$A3_sum), 4), "\n")
t_test_result <- t.test(high_sbs2_samples$A3_sum, low_sbs2_samples$A3_sum)
cat("  T-test p-value:", format(t_test_result$p.value, scientific = TRUE, digits = 3), "\n")

if(t_test_result$p.value > 0.05) {
  cat("  >> GOOD: No significant difference in A3 expression between groups\n\n")
} else {
  cat("  >> WARNING: Significant A3 expression difference detected\n\n")
}

# Save sample lists
high_sbs2_ids <- high_sbs2_samples$Entity_ID
low_sbs2_ids <- low_sbs2_samples$Entity_ID

# Save sample stratification summary
stratification_summary <- data.frame(
  Group = c("High_SBS2", "Low_SBS2"),
  N_samples = c(length(high_sbs2_ids), length(low_sbs2_ids)),
  Mean_A3_sum = c(mean(high_sbs2_samples$A3_sum), mean(low_sbs2_samples$A3_sum)),
  Mean_SBS2 = c(mean(high_sbs2_samples$SBS2), mean(low_sbs2_samples$SBS2)),
  Min_SBS2 = c(min(high_sbs2_samples$SBS2), min(low_sbs2_samples$SBS2)),
  Max_SBS2 = c(max(high_sbs2_samples$SBS2), max(low_sbs2_samples$SBS2))
)

write.table(stratification_summary, "Differential_Network_Sample_Stratification.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved sample stratification to: Differential_Network_Sample_Stratification.tsv\n\n")

# ============================================================
# PHASE 2: Full Expression Matrix Retrieval
# ============================================================
cat("============================================================\n")
cat("PHASE 2: Full Expression Matrix Retrieval\n")
cat("============================================================\n\n")

# Step 2.1: Load full FPKM-UQ data
cat("Step 2.1: Loading full FPKM-UQ expression matrix\n")
cat("------------------------------------------------------------\n")

full_expr_file <- "TCGA_master_FPKM_UQ.tsv"
cat("Loading:", full_expr_file, "(this may take a moment...)\n")

# Read the file with its complex header structure
full_expr_raw <- fread(full_expr_file, header = FALSE, sep = "\t", fill = TRUE)
cat("Raw dimensions:", nrow(full_expr_raw), "rows x", ncol(full_expr_raw), "columns\n")

# Parse header structure
# Row 1: Column headers (metadata columns + gene_ids)
# Row 2: Gene names
# Row 3: Gene types
col_headers <- as.character(full_expr_raw[1, ])
gene_names <- as.character(full_expr_raw[2, ])
gene_types <- as.character(full_expr_raw[3, ])

# Find Entity_ID column
entity_id_col <- which(col_headers == "Entity_ID")
cat("Entity_ID column:", entity_id_col, "\n")

# Gene columns start after metadata (columns 6+)
gene_start_col <- 6
n_genes_total <- ncol(full_expr_raw) - gene_start_col + 1
cat("Total genes in file:", n_genes_total, "\n")

# Extract data rows (rows 4+)
data_rows <- full_expr_raw[4:nrow(full_expr_raw), ]

# Get Entity_IDs from the data
all_entity_ids <- as.character(data_rows[[entity_id_col]])

# Step 2.2: Filter to selected samples
cat("\nStep 2.2: Filtering to selected sample groups\n")
cat("------------------------------------------------------------\n")

# Find row indices for each group
high_sbs2_row_idx <- which(all_entity_ids %in% high_sbs2_ids)
low_sbs2_row_idx <- which(all_entity_ids %in% low_sbs2_ids)

cat("High SBS2 samples found in expression matrix:", length(high_sbs2_row_idx), "\n")
cat("Low SBS2 samples found in expression matrix:", length(low_sbs2_row_idx), "\n")

# Extract expression matrices for each group
gene_cols <- gene_start_col:ncol(full_expr_raw)

# High SBS2 expression matrix (genes as rows, samples as columns after transpose)
expr_high_raw <- data_rows[high_sbs2_row_idx, ..gene_cols]
expr_high_matrix <- as.matrix(expr_high_raw)
expr_high_matrix <- apply(expr_high_matrix, 2, as.numeric)
rownames(expr_high_matrix) <- all_entity_ids[high_sbs2_row_idx]

# Low SBS2 expression matrix
expr_low_raw <- data_rows[low_sbs2_row_idx, ..gene_cols]
expr_low_matrix <- as.matrix(expr_low_raw)
expr_low_matrix <- apply(expr_low_matrix, 2, as.numeric)
rownames(expr_low_matrix) <- all_entity_ids[low_sbs2_row_idx]

# Transpose so genes are rows, samples are columns (standard format)
expr_high_matrix <- t(expr_high_matrix)
expr_low_matrix <- t(expr_low_matrix)

# Set gene names as row names
gene_symbols <- gene_names[gene_cols]
rownames(expr_high_matrix) <- gene_symbols
rownames(expr_low_matrix) <- gene_symbols

cat("\nExpression matrix dimensions:\n")
cat("  High SBS2:", nrow(expr_high_matrix), "genes x", ncol(expr_high_matrix), "samples\n")
cat("  Low SBS2:", nrow(expr_low_matrix), "genes x", ncol(expr_low_matrix), "samples\n")

# Step 2.3: Gene filtering
cat("\nStep 2.3: Filtering genes\n")
cat("------------------------------------------------------------\n")

# Get gene types for filtering
gene_type_vector <- gene_types[gene_cols]

# Create gene info dataframe - handle duplicate gene symbols by making them unique
gene_info <- data.frame(
  gene_symbol = gene_symbols,
  gene_type = gene_type_vector,
  stringsAsFactors = FALSE
)

# Make gene symbols unique by adding suffix for duplicates
gene_info <- gene_info %>%
  group_by(gene_symbol) %>%
  mutate(
    n_duplicates = n(),
    dup_index = row_number(),
    gene_symbol_unique = if_else(n_duplicates > 1, 
                                  paste0(gene_symbol, "_", dup_index), 
                                  gene_symbol)
  ) %>%
  ungroup()

# Update row names in expression matrices to use unique gene symbols
rownames(expr_high_matrix) <- gene_info$gene_symbol_unique
rownames(expr_low_matrix) <- gene_info$gene_symbol_unique

# Calculate gene statistics across both groups combined
combined_expr <- cbind(expr_high_matrix, expr_low_matrix)

gene_stats <- data.frame(
  gene_symbol_unique = rownames(combined_expr),
  mean_expr = rowMeans(combined_expr, na.rm = TRUE),
  variance = apply(combined_expr, 1, var, na.rm = TRUE),
  zero_fraction = rowMeans(combined_expr == 0, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Join with gene info using unique symbol
gene_stats <- gene_stats %>%
  left_join(gene_info %>% select(gene_symbol_unique, gene_symbol, gene_type), 
            by = "gene_symbol_unique")

cat("Gene filtering criteria:\n")
cat("  - Minimum mean expression:", PARAMS$min_expression, "\n")
cat("  - Minimum variance percentile:", PARAMS$min_variance_percentile * 100, "%\n")
cat("  - Maximum zero fraction:", PARAMS$max_zero_fraction * 100, "%\n")
cat("  - Protein coding only:", PARAMS$protein_coding_only, "\n")

# Apply filters
variance_threshold <- quantile(gene_stats$variance, PARAMS$min_variance_percentile, na.rm = TRUE)

genes_to_keep <- gene_stats %>%
  filter(mean_expr >= PARAMS$min_expression) %>%
  filter(variance >= variance_threshold) %>%
  filter(zero_fraction <= PARAMS$max_zero_fraction)

if(PARAMS$protein_coding_only) {
  genes_to_keep <- genes_to_keep %>%
    filter(gene_type == "protein_coding")
}

cat("\nGene filtering results:\n")
cat("  Starting genes:", nrow(gene_stats), "\n")
cat("  After mean expression filter:", sum(gene_stats$mean_expr >= PARAMS$min_expression), "\n")
cat("  After variance filter:", sum(gene_stats$variance >= variance_threshold, na.rm = TRUE), "\n")
cat("  After zero fraction filter:", sum(gene_stats$zero_fraction <= PARAMS$max_zero_fraction), "\n")
if(PARAMS$protein_coding_only) {
  cat("  After protein coding filter:", nrow(genes_to_keep), "\n")
}
cat("  Final genes:", nrow(genes_to_keep), "\n")

# Filter expression matrices to keep only selected genes
keep_gene_symbols <- genes_to_keep$gene_symbol_unique

expr_high_filtered <- expr_high_matrix[rownames(expr_high_matrix) %in% keep_gene_symbols, ]
expr_low_filtered <- expr_low_matrix[rownames(expr_low_matrix) %in% keep_gene_symbols, ]

# Ensure same genes in same order
common_genes <- intersect(rownames(expr_high_filtered), rownames(expr_low_filtered))
expr_high_filtered <- expr_high_filtered[common_genes, ]
expr_low_filtered <- expr_low_filtered[common_genes, ]

cat("\nFinal filtered expression matrices:\n")
cat("  High SBS2:", nrow(expr_high_filtered), "genes x", ncol(expr_high_filtered), "samples\n")
cat("  Low SBS2:", nrow(expr_low_filtered), "genes x", ncol(expr_low_filtered), "samples\n\n")

# Clean up large objects to free memory
rm(full_expr_raw, data_rows, expr_high_raw, expr_low_raw, combined_expr)
gc()

# ============================================================
# PHASE 3: Correlation Network Construction
# ============================================================
cat("============================================================\n")
cat("PHASE 3: Correlation Network Construction\n")
cat("============================================================\n\n")

# Step 3.1: Calculate correlation matrices
cat("Step 3.1: Calculating correlation matrices\n")
cat("------------------------------------------------------------\n")
cat("Using", PARAMS$correlation_method, "correlation\n")
cat("This may take several minutes for", nrow(expr_high_filtered), "genes...\n\n")

# Calculate correlation matrices
cat("Calculating high SBS2 correlation matrix...\n")
cor_high <- cor(t(expr_high_filtered), method = PARAMS$correlation_method, 
                use = "pairwise.complete.obs")
cat("  Done. Dimensions:", dim(cor_high)[1], "x", dim(cor_high)[2], "\n")

cat("Calculating low SBS2 correlation matrix...\n")
cor_low <- cor(t(expr_low_filtered), method = PARAMS$correlation_method,
               use = "pairwise.complete.obs")
cat("  Done. Dimensions:", dim(cor_low)[1], "x", dim(cor_low)[2], "\n\n")

# Step 3.2: Apply correlation threshold to create adjacency matrices
cat("Step 3.2: Creating adjacency matrices (threshold = ", 
    PARAMS$correlation_threshold, ")\n")
cat("------------------------------------------------------------\n")

# Create adjacency matrices (1 if |correlation| > threshold, 0 otherwise)
adj_high <- abs(cor_high) >= PARAMS$correlation_threshold
adj_low <- abs(cor_low) >= PARAMS$correlation_threshold

# Remove self-loops (diagonal)
diag(adj_high) <- FALSE
diag(adj_low) <- FALSE

# Convert to numeric
adj_high <- adj_high * 1
adj_low <- adj_low * 1

cat("High SBS2 network edges:", sum(adj_high) / 2, "\n")
cat("Low SBS2 network edges:", sum(adj_low) / 2, "\n")

# Step 3.3: Build network graphs
cat("\nStep 3.3: Building network graphs\n")
cat("------------------------------------------------------------\n")

# Create igraph objects
graph_high <- graph_from_adjacency_matrix(adj_high, mode = "undirected", 
                                          weighted = NULL, diag = FALSE)
graph_low <- graph_from_adjacency_matrix(adj_low, mode = "undirected",
                                         weighted = NULL, diag = FALSE)

# Add ABSOLUTE correlation values as edge weights (must be positive for layout)
# For high SBS2 network
edges_high <- as_edgelist(graph_high)
if(nrow(edges_high) > 0) {
  edge_weights_high <- sapply(1:nrow(edges_high), function(i) {
    abs(cor_high[edges_high[i, 1], edges_high[i, 2]])  # Use absolute value
  })
  E(graph_high)$weight <- edge_weights_high
}

# For low SBS2 network
edges_low <- as_edgelist(graph_low)
if(nrow(edges_low) > 0) {
  edge_weights_low <- sapply(1:nrow(edges_low), function(i) {
    abs(cor_low[edges_low[i, 1], edges_low[i, 2]])  # Use absolute value
  })
  E(graph_low)$weight <- edge_weights_low
}

cat("High SBS2 network: ", vcount(graph_high), " nodes, ", ecount(graph_high), " edges\n")
cat("Low SBS2 network: ", vcount(graph_low), " nodes, ", ecount(graph_low), " edges\n")

# Calculate network statistics (using updated function name)
cat("\nNetwork statistics:\n")
cat("  High SBS2:\n")
cat("    - Density:", round(edge_density(graph_high), 6), "\n")
cat("    - Mean degree:", round(mean(degree(graph_high)), 2), "\n")
cat("    - Transitivity:", round(transitivity(graph_high, type = "global"), 4), "\n")

cat("  Low SBS2:\n")
cat("    - Density:", round(edge_density(graph_low), 6), "\n")
cat("    - Mean degree:", round(mean(degree(graph_low)), 2), "\n")
cat("    - Transitivity:", round(transitivity(graph_low, type = "global"), 4), "\n\n")

# ============================================================
# PHASE 4: Differential Network Analysis
# ============================================================
cat("============================================================\n")
cat("PHASE 4: Differential Network Analysis\n")
cat("============================================================\n\n")

# Step 4.1: Calculate correlation differences
cat("Step 4.1: Calculating correlation differences\n")
cat("------------------------------------------------------------\n")

# Ensure matrices are aligned
common_genes <- intersect(rownames(cor_high), rownames(cor_low))
cor_high_aligned <- cor_high[common_genes, common_genes]
cor_low_aligned <- cor_low[common_genes, common_genes]

# Calculate difference matrix
diff_matrix <- cor_high_aligned - cor_low_aligned

cat("Correlation difference matrix dimensions:", dim(diff_matrix)[1], "x", dim(diff_matrix)[2], "\n")
cat("Difference statistics:\n")
cat("  Min:", round(min(diff_matrix, na.rm = TRUE), 4), "\n")
cat("  Max:", round(max(diff_matrix, na.rm = TRUE), 4), "\n")
cat("  Mean:", round(mean(diff_matrix, na.rm = TRUE), 6), "\n")
cat("  SD:", round(sd(diff_matrix, na.rm = TRUE), 4), "\n\n")

# Step 4.2: Calculate Differential Connectivity (DC) scores
cat("Step 4.2: Calculating Differential Connectivity scores\n")
cat("------------------------------------------------------------\n")

# DC score = sum of absolute correlation differences for each gene
dc_scores <- data.frame(
  gene = common_genes,
  dc_score = rowSums(abs(diff_matrix), na.rm = TRUE),
  mean_cor_high = rowMeans(abs(cor_high_aligned), na.rm = TRUE),
  mean_cor_low = rowMeans(abs(cor_low_aligned), na.rm = TRUE),
  degree_high = degree(graph_high)[common_genes],
  degree_low = degree(graph_low)[common_genes],
  stringsAsFactors = FALSE
)

# Calculate degree change
dc_scores$degree_change <- dc_scores$degree_high - dc_scores$degree_low

# Rank by DC score
dc_scores <- dc_scores %>%
  arrange(desc(dc_score)) %>%
  mutate(rank = row_number())

cat("Top 20 differentially connected genes:\n")
print(head(dc_scores, 20))
cat("\n")

# Step 4.3: Select top differential genes
cat("Step 4.3: Selecting top", PARAMS$n_top_differential_genes, "differential genes\n")
cat("------------------------------------------------------------\n")

top_dc_genes <- dc_scores %>%
  head(PARAMS$n_top_differential_genes)

cat("Selected", nrow(top_dc_genes), "top differential genes\n")
cat("DC score range:", round(min(top_dc_genes$dc_score), 2), "-", 
    round(max(top_dc_genes$dc_score), 2), "\n\n")

# Save differential genes results
dc_output_file <- "Differential_Network_DC_Scores.tsv"
write.table(dc_scores, dc_output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved all DC scores to:", dc_output_file, "\n")

top_dc_output_file <- "Differential_Network_Top_DC_Genes.tsv"
write.table(top_dc_genes, top_dc_output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved top DC genes to:", top_dc_output_file, "\n\n")

# ============================================================
# PHASE 5: Visualization
# ============================================================
cat("============================================================\n")
cat("PHASE 5: Visualization\n")
cat("============================================================\n\n")

# Step 5.1: Prepare network for visualization
cat("Step 5.1: Preparing network for visualization\n")
cat("------------------------------------------------------------\n")

# Get the list of top differential genes
diff_gene_list <- top_dc_genes$gene

# Add node attributes to high SBS2 network
V(graph_high)$is_differential <- V(graph_high)$name %in% diff_gene_list
V(graph_high)$dc_score <- dc_scores$dc_score[match(V(graph_high)$name, dc_scores$gene)]

# Replace NA dc_scores with 0
V(graph_high)$dc_score[is.na(V(graph_high)$dc_score)] <- 0

# Color nodes: bright red for differential, light gray for others
V(graph_high)$color <- ifelse(V(graph_high)$is_differential, 
                               PARAMS$sig_node_color, 
                               PARAMS$nonsig_node_color)

# Node size based on degree
V(graph_high)$size <- scales::rescale(degree(graph_high), to = c(2, 15))

# For differential genes, make size based on DC score
diff_indices <- which(V(graph_high)$is_differential)
if(length(diff_indices) > 0) {
  dc_scores_for_diff <- V(graph_high)$dc_score[diff_indices]
  if(length(dc_scores_for_diff) > 1 && max(dc_scores_for_diff) > min(dc_scores_for_diff)) {
    V(graph_high)$size[diff_indices] <- scales::rescale(dc_scores_for_diff, to = c(5, 20))
  } else {
    V(graph_high)$size[diff_indices] <- 10
  }
}

cat("Network prepared with", sum(V(graph_high)$is_differential), 
    "differential genes highlighted\n\n")

# Step 5.2: Create subnetwork if requested
if(!PARAMS$show_full_network) {
  cat("Step 5.2: Extracting subnetwork around differential genes\n")
  cat("------------------------------------------------------------\n")
  
  # Get neighbors of differential genes (within specified depth)
  diff_vertices <- which(V(graph_high)$is_differential)
  
  # Get neighborhood
  neighborhood_vertices <- unique(unlist(
    ego(graph_high, order = PARAMS$subnetwork_depth, nodes = diff_vertices, mode = "all")
  ))
  
  # Create subgraph
  graph_viz <- induced_subgraph(graph_high, neighborhood_vertices)
  
  cat("Subnetwork: ", vcount(graph_viz), " nodes, ", ecount(graph_viz), " edges\n")
  cat("  - Differential genes:", sum(V(graph_viz)$is_differential), "\n")
  cat("  - Neighbor genes:", sum(!V(graph_viz)$is_differential), "\n\n")
} else {
  graph_viz <- graph_high
  cat("Using full network for visualization\n")
  cat("  - Total nodes:", vcount(graph_viz), "\n")
  cat("  - Total edges:", ecount(graph_viz), "\n")
  cat("  - Differential genes:", sum(V(graph_viz)$is_differential), "\n\n")
}

# Step 5.3: Generate network plot
cat("Step 5.3: Generating network visualization\n")
cat("------------------------------------------------------------\n")

# Calculate layout
cat("Calculating network layout...\n")
set.seed(42)  # For reproducibility

if(vcount(graph_viz) > 0 && ecount(graph_viz) > 0) {
  
  # Use layout without weights if network is very large, otherwise use weights
  if(vcount(graph_viz) > 3000) {
    cat("Large network detected (", vcount(graph_viz), " nodes), using optimized layout...\n")
    layout_coords <- layout_with_fr(graph_viz, niter = 500)
  } else {
    # Ensure all weights are positive before layout
    if(!is.null(E(graph_viz)$weight)) {
      E(graph_viz)$weight <- abs(E(graph_viz)$weight)
      # Replace any zeros or NAs with small positive value
      E(graph_viz)$weight[is.na(E(graph_viz)$weight) | E(graph_viz)$weight == 0] <- 0.01
    }
    layout_coords <- layout_with_fr(graph_viz, niter = 1000, weights = E(graph_viz)$weight)
  }
  
  cat("Layout complete.\n")
  cat("Building ggraph visualization...\n")
  
  # Create ggraph plot with improved visualization settings
  p_network <- ggraph(graph_viz, layout = layout_coords) +
    # Draw edges first (underneath nodes) - DARK GRAY with adjustable width
    geom_edge_link(
      color = PARAMS$edge_color,
      width = PARAMS$edge_width,
      alpha = 0.4
    ) +
    # Draw non-differential nodes - LIGHT GRAY
    geom_node_point(
      data = . %>% filter(!is_differential),
      aes(size = size),
      color = PARAMS$nonsig_node_color,
      alpha = PARAMS$nonsig_node_alpha
    ) +
    # Draw differential nodes on top - BRIGHT RED
    geom_node_point(
      data = . %>% filter(is_differential),
      aes(size = size),
      color = PARAMS$sig_node_color,
      alpha = PARAMS$sig_node_alpha
    ) +
    # Add labels for differential genes
    {
      if(PARAMS$label_all_sig_genes) {
        # Label ALL differential genes
        geom_node_text(
          data = . %>% filter(is_differential),
          aes(label = name),
          size = PARAMS$label_size,
          repel = TRUE,
          max.overlaps = Inf,  # Allow all labels
          segment.color = "gray50",
          segment.size = 0.2
        )
      } else {
        # Label only top 20% of differential genes by DC score
        geom_node_text(
          data = . %>% filter(is_differential & dc_score >= quantile(dc_score[is_differential], 0.8, na.rm = TRUE)),
          aes(label = name),
          size = PARAMS$label_size,
          repel = TRUE,
          max.overlaps = 20,
          segment.color = "gray50",
          segment.size = 0.2
        )
      }
    } +
    scale_size_continuous(range = c(1, 8), guide = "none") +
    theme_void() +
    theme(
      plot.margin = margin(20, 20, 20, 20),
      legend.position = "none"
    ) +
    # Add annotation for sample counts
    annotate("text", x = Inf, y = Inf,
             label = paste0("High SBS2: n = ", length(high_sbs2_ids), "\n",
                           "Total genes: ", vcount(graph_viz), "\n",
                           "Differential genes: ", sum(V(graph_viz)$is_differential), "\n",
                           "Edges: ", format(ecount(graph_viz), big.mark = ",")),
             hjust = 1.1, vjust = 1.1, size = 5)
  
  # Save network plot
  network_plot_file <- "Figure_Differential_Network_High_SBS2.pdf"
  cat("Saving PDF (this may take a while for large networks)...\n")
  ggsave(network_plot_file, p_network, width = 16, height = 14, units = "in")
  cat("Saved network plot to:", network_plot_file, "\n")
  
  network_plot_png <- "Figure_Differential_Network_High_SBS2.png"
  cat("Saving PNG...\n")
  ggsave(network_plot_png, p_network, width = 16, height = 14, units = "in", dpi = 300)
  cat("Saved network plot to:", network_plot_png, "\n\n")
  
} else {
  cat("WARNING: Network is empty, skipping visualization\n\n")
}

# Step 5.4: Create summary visualization of DC scores
cat("Step 5.4: Creating DC score distribution plot\n")
cat("------------------------------------------------------------\n")

# Histogram of DC scores with top genes highlighted
p_dc_hist <- ggplot(dc_scores, aes(x = dc_score)) +
  geom_histogram(bins = 50, fill = "gray70", color = "black", alpha = 0.7) +
  geom_vline(xintercept = min(top_dc_genes$dc_score), 
             color = PARAMS$sig_node_color, linetype = "dashed", linewidth = 1) +
  annotate("text", x = min(top_dc_genes$dc_score), y = Inf,
           label = paste0("Top ", PARAMS$n_top_differential_genes, " genes"),
           hjust = -0.1, vjust = 2, color = PARAMS$sig_node_color, size = 5) +
  labs(x = "Differential Connectivity Score",
       y = "Number of Genes") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )

dc_hist_file <- "Figure_DC_Score_Distribution.pdf"
ggsave(dc_hist_file, p_dc_hist, width = 10, height = 8, units = "in")
cat("Saved DC distribution plot to:", dc_hist_file, "\n")

# Bar plot of top 30 DC genes
top_30_dc <- dc_scores %>% head(30)

p_dc_bar <- ggplot(top_30_dc, aes(x = reorder(gene, dc_score), y = dc_score)) +
  geom_bar(stat = "identity", fill = PARAMS$sig_node_color, alpha = 0.8) +
  coord_flip() +
  labs(x = "Gene",
       y = "Differential Connectivity Score") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 16)
  )

dc_bar_file <- "Figure_Top_DC_Genes_Barplot.pdf"
ggsave(dc_bar_file, p_dc_bar, width = 10, height = 12, units = "in")
cat("Saved top DC genes barplot to:", dc_bar_file, "\n\n")

# ============================================================
# PHASE 6: Summary Report
# ============================================================
cat("============================================================\n")
cat("ANALYSIS COMPLETE - Summary Report\n")
cat("============================================================\n\n")

cat("SAMPLE SELECTION:\n")
cat("  - Starting HNSC tumors:", nrow(hnsc_expression), "\n")
cat("  - After A3 expression filter (top 50%):", nrow(high_a3_samples), "\n")
cat("  - High SBS2 group (top 20%):", length(high_sbs2_ids), "\n")
cat("  - Low SBS2 group (bottom 20%):", length(low_sbs2_ids), "\n")

cat("\nGENE FILTERING:\n")
cat("  - Starting genes:", n_genes_total, "\n")
cat("  - After filtering:", length(common_genes), "\n")

cat("\nNETWORK STATISTICS:\n")
cat("  High SBS2 network:\n")
cat("    - Nodes:", vcount(graph_high), "\n")
cat("    - Edges:", ecount(graph_high), "\n")
cat("  Low SBS2 network:\n")
cat("    - Nodes:", vcount(graph_low), "\n")
cat("    - Edges:", ecount(graph_low), "\n")

cat("\nDIFFERENTIAL ANALYSIS:\n")
cat("  - Top differential genes selected:", PARAMS$n_top_differential_genes, "\n")
cat("  - Highest DC score:", round(max(dc_scores$dc_score), 2), "\n")
cat("  - Top 10 differential genes:\n")
print(dc_scores %>% head(10) %>% select(gene, dc_score, degree_high, degree_low))

cat("\nVISUALIZATION SETTINGS:\n")
cat("  - Edge color:", PARAMS$edge_color, "\n")
cat("  - Edge width:", PARAMS$edge_width, "\n")
cat("  - Non-significant node color:", PARAMS$nonsig_node_color, "\n")
cat("  - Significant node color:", PARAMS$sig_node_color, "\n")
cat("  - Label all significant genes:", PARAMS$label_all_sig_genes, "\n")

cat("\nOUTPUT FILES:\n")
cat("  1. Differential_Network_Sample_Stratification.tsv\n")
cat("  2. Differential_Network_DC_Scores.tsv\n")
cat("  3. Differential_Network_Top_DC_Genes.tsv\n")
cat("  4. Figure_Differential_Network_High_SBS2.pdf/png\n")
cat("  5. Figure_DC_Score_Distribution.pdf\n")
cat("  6. Figure_Top_DC_Genes_Barplot.pdf\n")

cat("\n============================================================\n")
cat("Analysis complete!\n")
cat("============================================================\n")
