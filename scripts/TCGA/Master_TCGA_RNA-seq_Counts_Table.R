library(data.table)
library(dplyr)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

# Check if manifest file exists
manifest_file <- "TCGA_RNAseq_master_manifest.tsv"

if(!file.exists(manifest_file)) {
  stop(paste0(
    "\n\nSIGTERM: Manifest file not found!\n",
    "Expected file: ", manifest_file, "\n",
    "Working directory: ", getwd(), "\n\n",
    "Please run Organize_RNA-seq_Counts_TCGA.R first to generate the manifest.\n",
    "The manifest file is required to map UUIDs to Entity_IDs.\n"
  ))
}

# Read the manifest file to get mapping
cat("Reading master manifest...\n")
manifest <- fread(manifest_file)

cat("Manifest structure:\n")
cat("Columns:", paste(colnames(manifest), collapse=", "), "\n")
cat("Total rows:", nrow(manifest), "\n")

# Validate Entity_ID quality in manifest
cat("\nEntity_ID validation in manifest:\n")
cat("Total entries:", nrow(manifest), "\n")
cat("Complete 28-char barcodes:", sum(nchar(manifest$Entity_ID) == 28), "\n")
cat("Non-28-char barcodes:", sum(nchar(manifest$Entity_ID) != 28), "\n")

if(any(nchar(manifest$Entity_ID) != 28)) {
  cat("\nWARNING: Some Entity_IDs in manifest are not 28 characters!\n")
  bad_entries <- manifest[nchar(manifest$Entity_ID) != 28, ]
  print(head(bad_entries, 5))
  
  stop(paste0(
    "\n\nSIGTERM: Invalid Entity_IDs in manifest!\n",
    "Found ", sum(nchar(manifest$Entity_ID) != 28), " entries with non-28-char Entity_IDs.\n",
    "Please regenerate the manifest by deleting ", manifest_file, "\n",
    "and re-running Organize_RNA-seq_Counts_TCGA.R\n"
  ))
}

# Show sample Entity_IDs
cat("\nSample Entity_IDs from manifest (first 5):\n")
print(head(manifest$Entity_ID, 5))
cat("\n")

# Function to parse file path and name
# Filename format: Case_ID_UUID.rna_seq.augmented_star_gene_counts.tsv
# Example: TCGA-OR-A5J1_fe16b2d3-17b0-4e24-ab31-62d2e951b3a2.rna_seq.augmented_star_gene_counts.tsv
parse_file_info <- function(filepath, manifest_df) {
  # Extract project ID from path (e.g., "TCGA-BRCA")
  path_parts <- strsplit(filepath, "/")[[1]]
  project_idx <- grep("^TCGA-", path_parts)
  project_id <- path_parts[project_idx[1]]  # First TCGA- match is the project
  
  # Extract tissue type (Tumor/Normal/Other)
  tissue_type <- path_parts[project_idx[1] + 1]
  
  # Parse filename: TCGA-OR-A5J1_6bacf042...tsv
  filename <- basename(filepath)
  # Remove the .rna_seq.augmented_star_gene_counts.tsv suffix
  base_name <- sub("\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", "", filename)
  
  # Split on underscore: TCGA-OR-A5J1_uuid
  name_parts <- strsplit(base_name, "_", fixed = TRUE)[[1]]
  case_id <- name_parts[1]  # e.g., TCGA-OR-A5J1
  file_uuid <- name_parts[2]  # e.g., fe16b2d3-17b0-4e24-ab31-62d2e951b3a2
  
  # Look up Entity_ID from manifest using UUID (most reliable match)
  matched_rows <- manifest_df %>%
    filter(UUID == file_uuid)
  
  if(nrow(matched_rows) == 0) {
    # Fallback: try matching on Project + Tissue_Status + Case_ID
    matched_rows <- manifest_df %>%
      filter(
        Project == project_id,
        Tissue_Status == tissue_type,
        Case_ID == case_id
      )
  }
  
  if(nrow(matched_rows) == 0) {
    # Save diagnostic info before error
    error_info <- data.frame(
      Filepath = filepath,
      Project_ID = project_id,
      Tissue_Type = tissue_type,
      Case_ID = case_id,
      File_UUID = file_uuid,
      stringsAsFactors = FALSE
    )
    write.table(error_info, "FAILED_manifest_lookup.tsv",
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    stop(paste0(
      "\n\nSIGTERM: No manifest match found!\n",
      "File: ", filepath, "\n",
      "Project: ", project_id, "\n",
      "Tissue: ", tissue_type, "\n", 
      "Case_ID: ", case_id, "\n",
      "UUID: ", file_uuid, "\n\n",
      "Please check FAILED_manifest_lookup.tsv for details.\n"
    ))
  }
  
  if(nrow(matched_rows) > 1) {
    # Multiple matches - use first one but warn
    warning("Multiple manifest matches for UUID: ", file_uuid, " - using first match")
  }
  
  entity_id <- matched_rows$Entity_ID[1]
  
  # Validate Entity_ID is exactly 28 characters
  if(nchar(entity_id) != 28) {
    error_info <- data.frame(
      Filepath = filepath,
      Project_ID = project_id,
      Tissue_Type = tissue_type,
      Case_ID = case_id,
      File_UUID = file_uuid,
      Entity_ID = entity_id,
      Entity_ID_Length = nchar(entity_id),
      stringsAsFactors = FALSE
    )
    write.table(error_info, "FAILED_Entity_ID_validation.tsv",
                sep = "\t", row.names = FALSE, col.names = !file.exists("FAILED_Entity_ID_validation.tsv"),
                quote = FALSE, append = file.exists("FAILED_Entity_ID_validation.tsv"))
    
    stop(paste0(
      "\n\nSIGTERM: Entity_ID is not 28 characters!\n",
      "File: ", filepath, "\n",
      "Case_ID: ", case_id, "\n",
      "UUID: ", file_uuid, "\n",
      "Entity_ID: '", entity_id, "'\n",
      "Length: ", nchar(entity_id), " (expected 28)\n\n",
      "Please check FAILED_Entity_ID_validation.tsv for details.\n"
    ))
  }
  
  return(data.frame(
    Project_ID = project_id,
    Tissue_Type = tissue_type,
    Case_ID = case_id,
    File_ID = file_uuid,
    Entity_ID = entity_id,
    stringsAsFactors = FALSE
  ))
}

# Check if organized_counts directory exists
if(!dir.exists("organized_counts")) {
  stop(paste0(
    "\n\nSIGTERM: organized_counts directory not found!\n",
    "Working directory: ", getwd(), "\n\n",
    "Please run Organize_RNA-seq_Counts_TCGA.R first to organize the files.\n"
  ))
}

# Get all TSV files
cat("Finding all TSV files...\n")
all_files <- list.files("organized_counts",
                        pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
                        recursive = TRUE,
                        full.names = TRUE)

cat("Found", length(all_files), "files\n")

if(length(all_files) == 0) {
  stop(paste0(
    "\n\nSIGTERM: No TSV files found in organized_counts directory!\n",
    "Working directory: ", getwd(), "\n\n",
    "Please run Organize_RNA-seq_Counts_TCGA.R first to organize the files.\n"
  ))
}

# Read first file to establish gene order and metadata
cat("Reading first file to establish gene structure...\n")
first_file <- fread(all_files[1], skip = 1)  # Skip the comment line

# Filter out the mapping stats rows (N_unmapped, etc.)
gene_info <- first_file[!grepl("^N_", gene_id)]

# Extract gene metadata
gene_ids <- gene_info$gene_id
gene_names <- gene_info$gene_name
gene_types <- gene_info$gene_type

cat("Found", length(gene_ids), "genes\n")

# Initialize data structures for three expression types
n_samples <- length(all_files)
n_genes <- length(gene_ids)

cat("\nInitializing matrices for", n_samples, "samples and", n_genes, "genes\n")
cat("This may take a moment...\n\n")

# Create metadata data frame (5 columns - no duplicate UUID)
sample_metadata <- data.frame(
  Project_ID = character(n_samples),
  Tissue_Type = character(n_samples),
  Case_ID = character(n_samples),
  File_ID = character(n_samples),
  Entity_ID = character(n_samples),
  stringsAsFactors = FALSE
)

# Initialize expression matrices
tpm_matrix <- matrix(NA, nrow = n_samples, ncol = n_genes)
fpkm_matrix <- matrix(NA, nrow = n_samples, ncol = n_genes)
fpkm_uq_matrix <- matrix(NA, nrow = n_samples, ncol = n_genes)

# Process each file
cat("Processing files...\n")
pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for(i in 1:length(all_files)) {
  filepath <- all_files[i]
  
  # Parse file information with manifest lookup
  file_info <- parse_file_info(filepath, manifest)
  sample_metadata[i, ] <- file_info[1, 1:5]
  
  # Read expression data
  expr_data <- fread(filepath, skip = 1)
  
  # Filter out mapping stats
  expr_data <- expr_data[!grepl("^N_", gene_id)]
  
  # Verify gene order matches (safety check)
  if(!all(expr_data$gene_id == gene_ids)) {
    warning("Gene order mismatch in file: ", filepath)
    # Reorder to match if needed
    expr_data <- expr_data[match(gene_ids, expr_data$gene_id), ]
  }
  
  # Extract expression values
  tpm_matrix[i, ] <- expr_data$tpm_unstranded
  fpkm_matrix[i, ] <- expr_data$fpkm_unstranded
  fpkm_uq_matrix[i, ] <- expr_data$fpkm_uq_unstranded
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nCreating final data frames...\n")

# Function to create final output with gene metadata as first two rows
create_master_file <- function(expr_matrix, sample_meta, gene_ids, gene_names, gene_types, expr_type) {
  
  cat("Building", expr_type, "master file...\n")
  
  # Create column names (gene IDs)
  colnames(expr_matrix) <- gene_ids
  
  # Combine metadata and expression
  full_data <- cbind(sample_meta, expr_matrix)
  
  # Create gene_name row with NA for metadata columns (5 NAs)
  gene_name_row <- c(NA, NA, NA, NA, NA, gene_names)
  
  # Create gene_type row with NA for metadata columns (5 NAs)
  gene_type_row <- c(NA, NA, NA, NA, NA, gene_types)
  
  # Create column headers row
  col_headers <- c("Project_ID", "Tissue_Type", "Case_ID", "File_ID", "Entity_ID", gene_ids)
  
  # Write output file
  output_file <- paste0("TCGA_master_", expr_type, ".tsv")
  
  cat("Writing to", output_file, "...\n")
  
  # Write header rows first
  write.table(t(col_headers), output_file,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(t(gene_name_row), output_file,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  write.table(t(gene_type_row), output_file,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  
  # Write data
  write.table(full_data, output_file,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  
  cat("Finished writing", output_file, "\n")
  cat("File size:", file.size(output_file) / 1e9, "GB\n\n")
  
  return(output_file)
}

# Create all three master files
tpm_file <- create_master_file(tpm_matrix, sample_metadata, gene_ids, gene_names, gene_types, "TPM")
fpkm_file <- create_master_file(fpkm_matrix, sample_metadata, gene_ids, gene_names, gene_types, "FPKM")
fpkm_uq_file <- create_master_file(fpkm_uq_matrix, sample_metadata, gene_ids, gene_names, gene_types, "FPKM_UQ")

# Summary report
cat("\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("CONSOLIDATION COMPLETE\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Total samples processed:", n_samples, "\n")
cat("Total genes:", n_genes, "\n")
cat("\nOutput files:\n")
cat("  1.", tpm_file, "\n")
cat("  2.", fpkm_file, "\n")
cat("  3.", fpkm_uq_file, "\n")
cat("\nFile structure:\n")
cat("  - Row 1: Column headers (Project_ID, Tissue_Type, Case_ID, File_ID, Entity_ID, gene_ids...)\n")
cat("  - Row 2: Gene names (NA, NA, NA, NA, NA, gene_names...)\n")
cat("  - Row 3: Gene types (NA, NA, NA, NA, NA, gene_types...)\n")
cat("  - Rows 4+: Sample expression data\n")
cat("\nSample breakdown by project:\n")
print(table(sample_metadata$Project_ID))
cat("\nSample breakdown by tissue type:\n")
print(table(sample_metadata$Tissue_Type))

# Final Entity_ID validation
cat("\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("ENTITY_ID FINAL VALIDATION\n")
cat(paste(rep("=", 60), collapse=""), "\n")
entity_lengths <- nchar(sample_metadata$Entity_ID)
cat("All Entity_IDs are 28 characters:", all(entity_lengths == 28), "\n")
cat("Entity_ID length distribution:\n")
print(table(entity_lengths))
cat("\nExample Entity_IDs (first 5):\n")
print(head(sample_metadata$Entity_ID, 5))
cat("\nExample Case_IDs (first 5):\n")
print(head(sample_metadata$Case_ID, 5))
cat(paste(rep("=", 60), collapse=""), "\n")

# Create a sample metadata file for reference
write.table(sample_metadata, "TCGA_sample_metadata_final.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSample metadata saved to: TCGA_sample_metadata_final.tsv\n")

# Show example of final output format
cat("\nExample of first 3 samples (metadata columns only):\n")
print(head(sample_metadata, 3))
