library(data.table)
library(dplyr)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("Preparing Matched Expression and Mutation Files\n")
cat("============================================================\n\n")

# Define the A3 gene symbols we want to keep (subset of 4)
a3_genes_subset <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3H")

cat("Target A3 genes:", paste(a3_genes_subset, collapse=", "), "\n\n")

# ============================================================
# STEP 1: Read and filter TCGA_master_FPKM_UQ.tsv
# ============================================================
cat("============================================================\n")
cat("STEP 1: Processing TCGA_master_FPKM_UQ.tsv\n")
cat("============================================================\n\n")

# Read the file
# Row 1: Column headers (Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID, gene_ids...)
# Row 2: Gene names (NA, NA, NA, NA, NA, gene_names...)
# Row 3: Gene types
# Rows 4+: Data

fpkm_uq_file <- "TCGA_master_FPKM_UQ.tsv"
cat("Reading:", fpkm_uq_file, "\n")

all_lines <- fread(fpkm_uq_file, header=FALSE, sep="\t", fill=TRUE)
cat("Original dimensions:", nrow(all_lines), "rows x", ncol(all_lines), "columns\n")

# Row 1 is column headers
col_headers <- as.character(all_lines[1, ])

# Row 2 is gene names (gene symbols)
gene_names <- as.character(all_lines[2, ])

# Identify column indices we need
# Metadata columns
project_id_col <- which(col_headers == "Project_ID")
tissue_type_col <- which(col_headers == "Tissue_Type")
entity_id_col <- which(col_headers == "Entity_ID")

cat("\nMetadata column indices:\n")
cat("  Project_ID: column", project_id_col, "\n")
cat("  Tissue_Type: column", tissue_type_col, "\n")
cat("  Entity_ID: column", entity_id_col, "\n")

# Find A3 gene columns by gene symbol (from row 2)
a3_col_indices <- c()
a3_col_names <- c()

cat("\nFinding A3 gene columns:\n")
for(gene in a3_genes_subset) {
  idx <- which(gene_names == gene)
  if(length(idx) > 0) {
    a3_col_indices <- c(a3_col_indices, idx[1])
    a3_col_names <- c(a3_col_names, gene)
    cat("  Found", gene, "at column", idx[1], "(gene_id:", col_headers[idx[1]], ")\n")
  } else {
    cat("  WARNING:", gene, "not found in gene names!\n")
  }
}

# Extract data rows (rows 4+)
data_rows <- all_lines[4:nrow(all_lines), ]

# Select the columns we want
cols_to_keep <- c(project_id_col, tissue_type_col, entity_id_col, a3_col_indices)
filtered_data <- data_rows[, ..cols_to_keep]

# Set new column names
new_col_names <- c("Project_ID", "Tissue_Type", "Entity_ID", a3_col_names)
colnames(filtered_data) <- new_col_names

# Convert expression values to numeric
for(gene in a3_col_names) {
  filtered_data[[gene]] <- as.numeric(filtered_data[[gene]])
}

cat("\nFiltered FPKM_UQ data:\n")
cat("  Dimensions:", nrow(filtered_data), "rows x", ncol(filtered_data), "columns\n")
cat("  Columns:", paste(colnames(filtered_data), collapse=", "), "\n")
cat("\nSample of filtered data (first 3 rows):\n")
print(head(filtered_data, 3))

# ============================================================
# STEP 2: Read Mutation_Table_Tumors_TCGA.tsv and get Entity_IDs
# ============================================================
cat("\n============================================================\n")
cat("STEP 2: Processing Mutation_Table_Tumors_TCGA.tsv\n")
cat("============================================================\n\n")

mutation_file <- "Mutation_Table_Tumors_TCGA.tsv"
cat("Reading:", mutation_file, "\n")

mutation_data <- fread(mutation_file)
cat("Dimensions:", nrow(mutation_data), "rows x", ncol(mutation_data), "columns\n")
cat("Columns:", paste(head(colnames(mutation_data), 10), collapse=", "), "...\n")

# Get the list of Entity_IDs from the mutation file
mutation_entity_ids <- mutation_data$TCGA_Gene_Expression_Entity_ID
cat("\nUnique Entity_IDs in mutation file:", length(unique(mutation_entity_ids)), "\n")
cat("Sample Entity_IDs:\n")
print(head(mutation_entity_ids, 5))

# ============================================================
# STEP 3: Filter FPKM_UQ data to match mutation file samples
# ============================================================
cat("\n============================================================\n")
cat("STEP 3: Filtering to matched samples\n")
cat("============================================================\n\n")

# Find common Entity_IDs
common_ids <- intersect(filtered_data$Entity_ID, mutation_entity_ids)
cat("Entity_IDs in FPKM_UQ file:", length(unique(filtered_data$Entity_ID)), "\n")
cat("Entity_IDs in mutation file:", length(unique(mutation_entity_ids)), "\n")
cat("Common Entity_IDs:", length(common_ids), "\n")

# Filter FPKM_UQ data to only include samples in mutation file
fpkm_uq_matched <- filtered_data %>%
  filter(Entity_ID %in% common_ids)

cat("\nFiltered FPKM_UQ data (matched samples):\n")
cat("  Dimensions:", nrow(fpkm_uq_matched), "rows x", ncol(fpkm_uq_matched), "columns\n")

# Also filter mutation data to only include samples in FPKM_UQ file (for consistency)
mutation_matched <- mutation_data %>%
  filter(TCGA_Gene_Expression_Entity_ID %in% common_ids)

cat("\nFiltered mutation data (matched samples):\n")
cat("  Dimensions:", nrow(mutation_matched), "rows x", ncol(mutation_matched), "columns\n")

# Verify both files have the same samples
fpkm_ids <- sort(unique(fpkm_uq_matched$Entity_ID))
mut_ids <- sort(unique(mutation_matched$TCGA_Gene_Expression_Entity_ID))

if(all(fpkm_ids == mut_ids)) {
  cat("\n>>> SUCCESS: Both files have identical Entity_IDs <<<\n")
} else {
  cat("\n>>> WARNING: Entity_IDs do not match perfectly <<<\n")
  cat("  In FPKM_UQ but not mutation:", length(setdiff(fpkm_ids, mut_ids)), "\n")
  cat("  In mutation but not FPKM_UQ:", length(setdiff(mut_ids, fpkm_ids)), "\n")
}

# ============================================================
# STEP 4: Save output files
# ============================================================
cat("\n============================================================\n")
cat("STEP 4: Saving output files\n")
cat("============================================================\n\n")

# Save filtered FPKM_UQ file
output_fpkm_file <- "A3_Expression_FPKM_UQ_Matched.tsv"
write.table(fpkm_uq_matched, output_fpkm_file, 
            sep="\t", row.names=FALSE, quote=FALSE)
cat("Saved:", output_fpkm_file, "\n")
cat("  Rows:", nrow(fpkm_uq_matched), "\n")
cat("  Columns:", paste(colnames(fpkm_uq_matched), collapse=", "), "\n\n")

# Save filtered mutation file
output_mutation_file <- "Mutation_Signatures_Matched.tsv"
write.table(mutation_matched, output_mutation_file,
            sep="\t", row.names=FALSE, quote=FALSE)
cat("Saved:", output_mutation_file, "\n")
cat("  Rows:", nrow(mutation_matched), "\n")
cat("  Columns:", ncol(mutation_matched), "total\n")

# ============================================================
# Summary
# ============================================================
cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")
cat("\nInput files:\n")
cat("  1. TCGA_master_FPKM_UQ.tsv\n")
cat("  2. Mutation_Table_Tumors_TCGA.tsv\n")

cat("\nOutput files (matched samples):\n")
cat("  1.", output_fpkm_file, "\n")
cat("     - Columns: Project_ID, Tissue_Type, Entity_ID, APOBEC3A, APOBEC3B, APOBEC3C, APOBEC3H\n")
cat("     - Rows:", nrow(fpkm_uq_matched), "samples\n")
cat("  2.", output_mutation_file, "\n")
cat("     - Columns: All original mutation signature columns\n")
cat("     - Rows:", nrow(mutation_matched), "samples\n")

cat("\nBoth files contain the same", length(common_ids), "samples based on Entity_ID matching.\n")
cat("============================================================\n")

# Show sample of final output
cat("\nSample of A3_Expression_FPKM_UQ_Matched.tsv (first 5 rows):\n")
print(head(fpkm_uq_matched, 5))

cat("\nSample of Mutation_Signatures_Matched.tsv (first 5 rows, first 8 columns):\n")
print(head(mutation_matched[, 1:8], 5))
