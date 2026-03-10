#!/usr/bin/env Rscript
# ============================================================
# Filter TCGA Master FPKM_UQ to HNSC Samples
# ============================================================
# 
# File structure:
#   Row 1: Column headers (Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID, gene_ids...)
#   Row 2: Gene symbols (NA, NA, NA, NA, NA, gene_symbols...)
#   Row 3: Gene types (NA, NA, NA, NA, NA, gene_types...)
#   Rows 4+: Sample data
#
# ============================================================

library(data.table)

# ============================================================
# Configuration
# ============================================================
input_path  <- "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1/TCGA_master_FPKM_UQ.tsv"
output_path <- "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1/TCGA_HNSC_master_FPKM_UQ.tsv"

cat("============================================================\n")
cat("Filtering TCGA Master FPKM_UQ to HNSC Samples\n")
cat("============================================================\n\n")

cat("Input file:", input_path, "\n")
cat("Output file:", output_path, "\n\n")

# ============================================================
# Step 1: Read the entire file preserving structure
# ============================================================
cat("Step 1: Reading input file...\n")

# Read all rows without headers (treat everything as data)
all_data <- fread(input_path, header = FALSE, sep = "\t", fill = TRUE)

cat("  Total rows read:", nrow(all_data), "\n")
cat("  Total columns:", ncol(all_data), "\n")

# ============================================================
# Step 2: Separate header rows from data rows
# ============================================================
cat("\nStep 2: Separating header rows from data...\n")

# Rows 1-3 are headers (column names, gene symbols, gene types)
header_rows <- all_data[1:3, ]

# Rows 4+ are sample data
data_rows <- all_data[4:nrow(all_data), ]

cat("  Header rows: 3\n")
cat("  Data rows:", nrow(data_rows), "\n")

# ============================================================
# Step 3: Identify Project_ID column and filter for HNSC
# ============================================================
cat("\nStep 3: Filtering for TCGA-HNSC samples...\n")

# Get column headers from row 1
col_headers <- as.character(header_rows[1, ])

# Find Project_ID column
project_id_col <- which(col_headers == "Project_ID")

if(length(project_id_col) == 0) {
  stop("ERROR: Could not find 'Project_ID' column!")
}

cat("  Project_ID column index:", project_id_col, "\n")

# Get Project_IDs from data rows
project_ids <- as.character(data_rows[[project_id_col]])

# Show distribution of cancer types
cat("\n  Cancer type distribution:\n")
cancer_counts <- table(project_ids)
print(head(sort(cancer_counts, decreasing = TRUE), 10))

# Filter for HNSC
hnsc_rows <- which(project_ids == "TCGA-HNSC")
cat("\n  TCGA-HNSC samples found:", length(hnsc_rows), "\n")

if(length(hnsc_rows) == 0) {
  stop("ERROR: No TCGA-HNSC samples found!")
}

# Filter data to HNSC only
data_filtered <- data_rows[hnsc_rows, ]

cat("  Filtered data rows:", nrow(data_filtered), "\n")

# ============================================================
# Step 4: Combine headers with filtered data
# ============================================================
cat("\nStep 4: Combining headers with filtered data...\n")

# Combine: header rows (1-3) + filtered data rows
output_data <- rbind(header_rows, data_filtered)

cat("  Output total rows:", nrow(output_data), "\n")
cat("    - Header rows: 3\n")
cat("    - Data rows:", nrow(data_filtered), "\n")

# ============================================================
# Step 5: Write output file
# ============================================================
cat("\nStep 5: Writing output file...\n")

fwrite(output_data, output_path, sep = "\t", col.names = FALSE, quote = FALSE)

cat("  File saved to:", output_path, "\n")

# ============================================================
# Step 6: Validation
# ============================================================
cat("\nStep 6: Validating output file...\n")

# Quick validation - read back and check
validation_data <- fread(output_path, header = FALSE, sep = "\t", nrows = 10)
validation_project_ids <- as.character(validation_data[4:nrow(validation_data), ..project_id_col][[1]])

if(all(validation_project_ids == "TCGA-HNSC")) {
  cat("  ✓ Validation passed: Output contains only TCGA-HNSC samples\n")
} else {
  cat("  ✗ Validation warning: Unexpected Project_IDs found\n")
}

# Check file size
file_size <- file.info(output_path)$size / (1024^2)  # Convert to MB
cat("  Output file size:", round(file_size, 2), "MB\n")

# ============================================================
# Summary
# ============================================================
cat("\n============================================================\n")
cat("COMPLETE\n")
cat("============================================================\n")
cat("Input samples:", nrow(data_rows), "\n")
cat("Output samples:", nrow(data_filtered), "(TCGA-HNSC only)\n")
cat("Output file:", output_path, "\n")
cat("============================================================\n")
