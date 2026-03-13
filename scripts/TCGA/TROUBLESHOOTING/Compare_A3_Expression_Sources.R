library(data.table)
library(dplyr)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("TCGA Expression Data Source Identification\n")
cat("============================================================\n\n")

# Define the A3 gene symbols we're looking for
a3_genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", 
              "APOBEC3F", "APOBEC3G", "APOBEC3H")

cat("Target A3 genes:", paste(a3_genes, collapse=", "), "\n\n")

# ============================================================
# PART 1: Read and filter the mystery file
# ============================================================
cat("============================================================\n")
cat("PART 1: Loading mystery file\n")
cat("============================================================\n")

mystery_file <- "A3_Table_All_Samples_TCGA_Mohadeseh.tsv"
mystery_data <- fread(mystery_file)

cat("Mystery file:", mystery_file, "\n")
cat("Dimensions:", nrow(mystery_data), "rows x", ncol(mystery_data), "columns\n")
cat("Columns:", paste(colnames(mystery_data), collapse=", "), "\n\n")

# Filter mystery file - keep Entity_ID and A3 genes, drop UUID column
mystery_filtered <- mystery_data %>%
  select(TCGA_Gene_Expression_Entity_ID, all_of(a3_genes)) %>%
  rename(Entity_ID = TCGA_Gene_Expression_Entity_ID)

# Check for duplicates in mystery file
cat("Checking for duplicate Entity_IDs in mystery file...\n")
dup_mystery <- sum(duplicated(mystery_filtered$Entity_ID))
cat("  Duplicate Entity_IDs:", dup_mystery, "\n")

if(dup_mystery > 0) {
  cat("  Removing duplicates (keeping first occurrence)...\n")
  mystery_filtered <- mystery_filtered %>% distinct(Entity_ID, .keep_all = TRUE)
}

cat("Mystery file filtered to", nrow(mystery_filtered), "rows x", ncol(mystery_filtered), "columns\n")
cat("Sample of mystery data (first 3 rows):\n")
print(head(mystery_filtered, 3))
cat("\n")

# Save filtered mystery file
mystery_filtered_file <- "A3_Table_Mystery_Filtered.tsv"
write.table(mystery_filtered, mystery_filtered_file, 
            sep="\t", row.names=FALSE, quote=FALSE)
cat("Saved filtered mystery file to:", mystery_filtered_file, "\n\n")

# ============================================================
# PART 2: Function to read and filter master files
# ============================================================

filter_master_file <- function(filepath, a3_gene_symbols) {
  cat("Processing:", filepath, "\n")
  
  # Read the file - it has a complex structure
  # Row 1: Column headers (Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID, gene_ids...)
  # Row 2: Gene names (NA, NA, NA, NA, NA, gene_names...)
  # Row 3: Gene types
  # Rows 4+: Data
  
  # Read all lines to parse structure
  all_lines <- fread(filepath, header=FALSE, sep="\t", fill=TRUE)
  
  cat("  File dimensions:", nrow(all_lines), "rows x", ncol(all_lines), "columns\n")
  
  # Row 1 is column headers (gene_ids in columns 6+)
  col_headers <- as.character(all_lines[1, ])
  
  # Row 2 is gene names (gene symbols in columns 6+)
  gene_names <- as.character(all_lines[2, ])
  
  # Find which columns have our A3 genes (by gene name/symbol)
  a3_col_indices <- c()
  a3_col_names <- c()
  
  for(gene in a3_gene_symbols) {
    idx <- which(gene_names == gene)
    if(length(idx) > 0) {
      a3_col_indices <- c(a3_col_indices, idx[1])
      a3_col_names <- c(a3_col_names, gene)
      cat("  Found", gene, "at column", idx[1], "(gene_id:", col_headers[idx[1]], ")\n")
    } else {
      cat("  WARNING:", gene, "not found in gene names!\n")
    }
  }
  
  # Entity_ID is column 5
  entity_col <- 5
  
  # Extract data rows (rows 4+)
  data_rows <- all_lines[4:nrow(all_lines), ]
  
  # Extract Entity_ID and A3 gene columns
  filtered_data <- data_rows[, c(entity_col, a3_col_indices), with=FALSE]
  
  # Set column names
  colnames(filtered_data) <- c("Entity_ID", a3_col_names)
  
  # Convert expression values to numeric
  for(gene in a3_col_names) {
    filtered_data[[gene]] <- as.numeric(filtered_data[[gene]])
  }
  
  # Check for duplicates
  dup_count <- sum(duplicated(filtered_data$Entity_ID))
  cat("  Duplicate Entity_IDs:", dup_count, "\n")
  
  if(dup_count > 0) {
    cat("  Removing duplicates (keeping first occurrence)...\n")
    filtered_data <- filtered_data %>% distinct(Entity_ID, .keep_all = TRUE)
  }
  
  cat("  Extracted", nrow(filtered_data), "samples with", ncol(filtered_data), "columns\n\n")
  
  return(filtered_data)
}

# ============================================================
# PART 3: Filter all three master files
# ============================================================
cat("============================================================\n")
cat("PART 2: Filtering master files\n")
cat("============================================================\n\n")

# Filter TPM file
tpm_filtered <- filter_master_file("TCGA_master_TPM.tsv", a3_genes)
tpm_filtered_file <- "A3_Table_TCGA_TPM_Filtered.tsv"
write.table(tpm_filtered, tpm_filtered_file, sep="\t", row.names=FALSE, quote=FALSE)
cat("Saved:", tpm_filtered_file, "\n\n")

# Filter FPKM file
fpkm_filtered <- filter_master_file("TCGA_master_FPKM.tsv", a3_genes)
fpkm_filtered_file <- "A3_Table_TCGA_FPKM_Filtered.tsv"
write.table(fpkm_filtered, fpkm_filtered_file, sep="\t", row.names=FALSE, quote=FALSE)
cat("Saved:", fpkm_filtered_file, "\n\n")

# Filter FPKM_UQ file
fpkm_uq_filtered <- filter_master_file("TCGA_master_FPKM_UQ.tsv", a3_genes)
fpkm_uq_filtered_file <- "A3_Table_TCGA_FPKM_UQ_Filtered.tsv"
write.table(fpkm_uq_filtered, fpkm_uq_filtered_file, sep="\t", row.names=FALSE, quote=FALSE)
cat("Saved:", fpkm_uq_filtered_file, "\n\n")

# ============================================================
# PART 4: Compare mystery file with each master file
# ============================================================
cat("============================================================\n")
cat("PART 3: Comparing mystery file with master files\n")
cat("============================================================\n\n")

compare_datasets <- function(mystery_df, master_df, master_name) {
  cat("Comparing mystery file with", master_name, "...\n")
  
  # Diagnostic info
  cat("  Mystery file rows:", nrow(mystery_df), "\n")
  cat("  Master file rows:", nrow(master_df), "\n")
  
  # Find common Entity_IDs
  common_ids <- intersect(mystery_df$Entity_ID, master_df$Entity_ID)
  cat("  Common Entity_IDs:", length(common_ids), "\n")
  
  if(length(common_ids) == 0) {
    cat("  No common Entity_IDs found!\n")
    cat("  Sample mystery Entity_IDs:", paste(head(mystery_df$Entity_ID, 3), collapse=", "), "\n")
    cat("  Sample master Entity_IDs:", paste(head(master_df$Entity_ID, 3), collapse=", "), "\n\n")
    return(list(
      name = master_name,
      common_samples = 0,
      correlation = NA,
      mean_diff = NA,
      exact_matches = 0,
      exact_match_pct = 0,
      close_matches = 0,
      close_match_pct = 0
    ))
  }
  
  # CRITICAL FIX: Filter BOTH dataframes to only common Entity_IDs
  # Then merge them to ensure perfect row alignment
  mystery_common <- mystery_df %>% 
    filter(Entity_ID %in% common_ids) %>%
    distinct(Entity_ID, .keep_all = TRUE) %>%
    arrange(Entity_ID)
  
  master_common <- master_df %>% 
    filter(Entity_ID %in% common_ids) %>%
    distinct(Entity_ID, .keep_all = TRUE) %>%
    arrange(Entity_ID)
  
  cat("  After filtering to common IDs:\n")
  cat("    Mystery rows:", nrow(mystery_common), "\n")
  cat("    Master rows:", nrow(master_common), "\n")
  
  # Verify the Entity_IDs are now aligned
  if(!all(mystery_common$Entity_ID == master_common$Entity_ID)) {
    cat("  WARNING: Entity_IDs not aligned after filtering, using merge...\n")
    
    # Use merge to ensure alignment
    merged_data <- merge(
      mystery_common, 
      master_common, 
      by = "Entity_ID", 
      suffixes = c("_mystery", "_master")
    )
    
    cat("    Merged rows:", nrow(merged_data), "\n")
    
    # Get gene columns
    gene_cols <- a3_genes
    
    # Extract values from merged data
    mystery_values <- c()
    master_values <- c()
    
    for(gene in gene_cols) {
      mystery_col <- paste0(gene, "_mystery")
      master_col <- paste0(gene, "_master")
      
      if(mystery_col %in% colnames(merged_data) && master_col %in% colnames(merged_data)) {
        mystery_values <- c(mystery_values, merged_data[[mystery_col]])
        master_values <- c(master_values, merged_data[[master_col]])
      }
    }
  } else {
    cat("    Entity_IDs aligned correctly\n")
    
    # Get gene columns (exclude Entity_ID)
    gene_cols <- setdiff(colnames(mystery_common), "Entity_ID")
    
    # Flatten all expression values for comparison
    mystery_values <- c()
    master_values <- c()
    
    for(gene in gene_cols) {
      mystery_values <- c(mystery_values, mystery_common[[gene]])
      master_values <- c(master_values, master_common[[gene]])
    }
  }
  
  # Final dimension check
  cat("  Values to compare: mystery =", length(mystery_values), ", master =", length(master_values), "\n")
  
  if(length(mystery_values) != length(master_values)) {
    cat("  ERROR: Dimension mismatch persists!\n\n")
    return(list(
      name = master_name,
      common_samples = length(common_ids),
      correlation = NA,
      mean_diff = NA,
      exact_matches = 0,
      exact_match_pct = 0,
      close_matches = 0,
      close_match_pct = 0
    ))
  }
  
  # Calculate correlation
  correlation <- cor(mystery_values, master_values, use = "complete.obs")
  
  # Calculate mean absolute difference
  mean_diff <- mean(abs(mystery_values - master_values), na.rm = TRUE)
  
  # Count exact matches (within floating point tolerance)
  exact_matches <- sum(abs(mystery_values - master_values) < 1e-6, na.rm = TRUE)
  
  # Count close matches (within 0.01)
  close_matches <- sum(abs(mystery_values - master_values) < 0.01, na.rm = TRUE)
  
  total_comparisons <- length(mystery_values)
  
  cat("  Total value comparisons:", total_comparisons, "\n")
  cat("  Correlation:", round(correlation, 6), "\n")
  cat("  Mean absolute difference:", round(mean_diff, 6), "\n")
  cat("  Exact matches (diff < 1e-6):", exact_matches, "(", 
      round(100 * exact_matches / total_comparisons, 2), "%)\n")
  cat("  Close matches (diff < 0.01):", close_matches, "(",
      round(100 * close_matches / total_comparisons, 2), "%)\n")
  
  # Show sample comparisons
  cat("\n  Sample value comparisons (first 5 common samples, APOBEC3A):\n")
  
  # Get aligned data for display
  display_mystery <- mystery_common %>% head(5)
  display_master <- master_common %>% head(5)
  
  comparison_sample <- data.frame(
    Entity_ID = display_mystery$Entity_ID,
    Mystery = display_mystery$APOBEC3A,
    Master = display_master$APOBEC3A
  )
  comparison_sample$Difference <- comparison_sample$Mystery - comparison_sample$Master
  print(comparison_sample)
  cat("\n")
  
  return(list(
    name = master_name,
    common_samples = length(common_ids),
    total_comparisons = total_comparisons,
    correlation = correlation,
    mean_diff = mean_diff,
    exact_matches = exact_matches,
    exact_match_pct = 100 * exact_matches / total_comparisons,
    close_matches = close_matches,
    close_match_pct = 100 * close_matches / total_comparisons
  ))
}

# Compare with each master file
results_tpm <- compare_datasets(mystery_filtered, tpm_filtered, "TPM")
results_fpkm <- compare_datasets(mystery_filtered, fpkm_filtered, "FPKM")
results_fpkm_uq <- compare_datasets(mystery_filtered, fpkm_uq_filtered, "FPKM_UQ")

# ============================================================
# PART 5: Generate comparison report
# ============================================================
cat("============================================================\n")
cat("PART 4: COMPARISON REPORT\n")
cat("============================================================\n\n")

# Create summary table
summary_df <- data.frame(
  Data_Type = c("TPM", "FPKM", "FPKM_UQ"),
  Common_Samples = c(results_tpm$common_samples, results_fpkm$common_samples, results_fpkm_uq$common_samples),
  Correlation = c(results_tpm$correlation, results_fpkm$correlation, results_fpkm_uq$correlation),
  Mean_Abs_Diff = c(results_tpm$mean_diff, results_fpkm$mean_diff, results_fpkm_uq$mean_diff),
  Exact_Match_Pct = c(results_tpm$exact_match_pct, results_fpkm$exact_match_pct, results_fpkm_uq$exact_match_pct),
  Close_Match_Pct = c(results_tpm$close_match_pct, results_fpkm$close_match_pct, results_fpkm_uq$close_match_pct)
)

cat("SUMMARY TABLE:\n")
print(summary_df)
cat("\n")

# Determine the best match
best_match_idx <- which.max(summary_df$Correlation)
best_match <- summary_df$Data_Type[best_match_idx]

# Additional criteria
highest_exact <- summary_df$Data_Type[which.max(summary_df$Exact_Match_Pct)]
lowest_diff <- summary_df$Data_Type[which.min(summary_df$Mean_Abs_Diff)]

cat("============================================================\n")
cat("CONCLUSION\n")
cat("============================================================\n")
cat("\n")
cat("Based on the comparison:\n")
cat("  - Highest correlation:", best_match, "(r =", round(max(summary_df$Correlation, na.rm = TRUE), 6), ")\n")
cat("  - Highest exact match %:", highest_exact, "(", round(max(summary_df$Exact_Match_Pct, na.rm = TRUE), 2), "%)\n")
cat("  - Lowest mean difference:", lowest_diff, "(", round(min(summary_df$Mean_Abs_Diff, na.rm = TRUE), 6), ")\n")
cat("\n")

if(max(summary_df$Correlation, na.rm = TRUE) > 0.999 && max(summary_df$Exact_Match_Pct, na.rm = TRUE) > 90) {
  cat(">>> STRONG MATCH: The mystery file 'A3_Table_All_Samples_TCGA_Mohadeseh.tsv'\n")
  cat("    appears to be ", best_match, " expression data.\n", sep = "")
} else if(max(summary_df$Correlation, na.rm = TRUE) > 0.99) {
  cat(">>> LIKELY MATCH: The mystery file appears to be ", best_match, " expression data,\n", sep = "")
  cat("    but with some differences (possibly different processing or version).\n")
} else {
  cat(">>> UNCLEAR MATCH: The mystery file does not closely match any of our\n")
  cat("    generated files. It may be from a different source or processing pipeline.\n")
}

cat("\n============================================================\n")

# Save the report
report_file <- "A3_Expression_Source_Identification_Report.txt"
sink(report_file)
cat("============================================================\n")
cat("TCGA A3 Expression Data Source Identification Report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n\n")

cat("Mystery file: A3_Table_All_Samples_TCGA_Mohadeseh.tsv\n")
cat("Mystery file samples:", nrow(mystery_filtered), "\n\n")

cat("Comparison files:\n")
cat("  - TCGA_master_TPM.tsv (", nrow(tpm_filtered), " samples)\n", sep = "")
cat("  - TCGA_master_FPKM.tsv (", nrow(fpkm_filtered), " samples)\n", sep = "")
cat("  - TCGA_master_FPKM_UQ.tsv (", nrow(fpkm_uq_filtered), " samples)\n\n", sep = "")

cat("A3 genes compared:", paste(a3_genes, collapse = ", "), "\n\n")

cat("SUMMARY TABLE:\n")
cat("-------------\n")
print(summary_df)
cat("\n\n")

cat("CONCLUSION:\n")
cat("-----------\n")
cat("Highest correlation:", best_match, "(r =", round(max(summary_df$Correlation, na.rm = TRUE), 6), ")\n")
cat("Highest exact match %:", highest_exact, "(", round(max(summary_df$Exact_Match_Pct, na.rm = TRUE), 2), "%)\n")
cat("Lowest mean difference:", lowest_diff, "(", round(min(summary_df$Mean_Abs_Diff, na.rm = TRUE), 6), ")\n\n")

if(max(summary_df$Correlation, na.rm = TRUE) > 0.999 && max(summary_df$Exact_Match_Pct, na.rm = TRUE) > 90) {
  cat("STRONG MATCH: The mystery file 'A3_Table_All_Samples_TCGA_Mohadeseh.tsv'\n")
  cat("appears to be ", best_match, " expression data.\n", sep = "")
} else if(max(summary_df$Correlation, na.rm = TRUE) > 0.99) {
  cat("LIKELY MATCH: The mystery file appears to be ", best_match, " expression data,\n", sep = "")
  cat("but with some differences (possibly different processing or version).\n")
} else {
  cat("UNCLEAR MATCH: The mystery file does not closely match any of our\n")
  cat("generated files. It may be from a different source or processing pipeline.\n")
}

cat("\n\nFILTERED FILES GENERATED:\n")
cat("-------------------------\n")
cat("1. A3_Table_Mystery_Filtered.tsv - Mystery file with Entity_ID + 7 A3 genes\n")
cat("2. A3_Table_TCGA_TPM_Filtered.tsv - TPM values for Entity_ID + 7 A3 genes\n")
cat("3. A3_Table_TCGA_FPKM_Filtered.tsv - FPKM values for Entity_ID + 7 A3 genes\n")
cat("4. A3_Table_TCGA_FPKM_UQ_Filtered.tsv - FPKM_UQ values for Entity_ID + 7 A3 genes\n")

sink()

cat("\nReport saved to:", report_file, "\n")
cat("\nFiltered files created:\n")
cat("  1.", mystery_filtered_file, "\n")
cat("  2.", tpm_filtered_file, "\n")
cat("  3.", fpkm_filtered_file, "\n")
cat("  4.", fpkm_uq_filtered_file, "\n")
