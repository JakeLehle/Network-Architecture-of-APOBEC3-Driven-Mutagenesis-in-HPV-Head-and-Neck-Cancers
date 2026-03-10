library(TCGAbiolinks)
library(dplyr)
library(tidyr)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

# List of all TCGA cancer types
cancer_types <- c("BRCA", "LUAD", "LUSC", "PRAD", "COAD", "STAD",
                  "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML",
                  "PAAD", "ESCA", "PCPG", "READ", "TGCT", "THYM",
                  "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS",
                  "CHOL", "GBM", "HNSC", "KIRC", "LGG", "OV",
                  "SKCM", "THCA", "UCEC")

# Check if manifest exists
manifest_file <- "TCGA_RNAseq_master_manifest.tsv"

if(file.exists(manifest_file)) {
  cat("Manifest file already exists:", manifest_file, "\n")
  cat("Loading existing manifest...\n")
  master_manifest <- read.table(manifest_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  cat("Loaded", nrow(master_manifest), "samples from existing manifest\n\n")
  
  # Validate existing manifest
  cat("Validating existing manifest Entity_IDs...\n")
  entity_lengths <- nchar(master_manifest$Entity_ID)
  cat("28-char Entity_IDs:", sum(entity_lengths == 28), "/", nrow(master_manifest), "\n")
  
  if(all(entity_lengths == 28)) {
    cat("All Entity_IDs are valid. Skipping manifest regeneration.\n\n")
    skip_manifest_generation <- TRUE
  } else {
    cat("Some Entity_IDs are invalid. Regenerating manifest...\n\n")
    skip_manifest_generation <- FALSE
  }
} else {
  cat("Manifest file not found. Generating new manifest...\n\n")
  skip_manifest_generation <- FALSE
}

if(!skip_manifest_generation) {
  
  cat("Generating manifest with full TCGA aliquot barcodes...\n\n")
  
  # Verify the 'cases' column contains full barcodes in a test query
  cat("Verifying 'cases' column contains full 28-char barcodes...\n")
  test_query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  test_metadata <- getResults(test_query)
  
  if(!"cases" %in% colnames(test_metadata)) {
    stop("SIGTERM: 'cases' column not found in GDC metadata!")
  }
  
  # Check that cases column has 28-char barcodes
  test_cases <- as.character(test_metadata$cases[1:min(5, nrow(test_metadata))])
  cat("Sample 'cases' values:\n")
  for(v in test_cases) {
    cat("  '", v, "' (", nchar(v), " chars)\n", sep="")
  }
  
  if(!all(nchar(test_cases) == 28)) {
    cat("\nWARNING: Not all 'cases' values are 28 characters!\n")
    cat("Proceeding anyway but check results carefully.\n\n")
  } else {
    cat("\nConfirmed: 'cases' column contains 28-char aliquot barcodes.\n\n")
  }
  
  # Build comprehensive manifest
  all_metadata <- list()
  failed_samples <- list()
  
  for(cancer in cancer_types) {
    cat("\n", paste(rep("=", 50), collapse=""), "\n", sep="")
    cat("Processing metadata for TCGA-", cancer, "\n", sep="")
    cat(paste(rep("=", 50), collapse=""), "\n")
    
    # Query to get metadata
    query <- GDCquery(
      project = paste0("TCGA-", cancer),
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    # Extract metadata from query results
    metadata <- getResults(query)
    cat("Found", nrow(metadata), "samples in GDC query\n")
    
    # Find the data directory for this cancer type
    data_dir <- file.path("GDCdata", paste0("TCGA-", cancer),
                          "Transcriptome_Profiling",
                          "Gene_Expression_Quantification")
    
    if(!dir.exists(data_dir)) {
      cat("  WARNING: Data directory not found for TCGA-", cancer, "\n", sep="")
      cat("  Skipping this cancer type.\n")
      next
    }
    
    # Build metadata for each sample
    cancer_meta_list <- list()
    
    for(i in 1:nrow(metadata)) {
      file_uuid <- metadata$file_id[i]
      
      # Case_ID from cases.submitter_id (e.g., TCGA-BH-A18H)
      case_id <- metadata$cases.submitter_id[i]
      
      # Entity_ID from cases column (full 28-char barcode, e.g., TCGA-BH-A18H-01A-11R-A12D-07)
      entity_id <- as.character(metadata$cases[i])
      
      sample_type <- metadata$sample_type[i]
      
      # Validate Entity_ID is exactly 28 characters
      if(is.na(entity_id) || nchar(entity_id) != 28) {
        # Record the failure for error reporting
        failed_samples[[length(failed_samples) + 1]] <- data.frame(
          Cancer_Type = cancer,
          Case_ID = case_id,
          UUID = file_uuid,
          Attempted_Entity_ID = ifelse(is.na(entity_id), "NA", entity_id),
          Entity_ID_Length = ifelse(is.na(entity_id), 0, nchar(entity_id)),
          stringsAsFactors = FALSE
        )
        
        cat("\n")
        cat("!!! CRITICAL ERROR: Entity_ID is not 28 characters !!!\n")
        cat("Cancer Type: TCGA-", cancer, "\n", sep="")
        cat("Case_ID: ", case_id, "\n", sep="")
        cat("UUID: ", file_uuid, "\n", sep="")
        cat("Entity_ID from 'cases' column: ", ifelse(is.na(entity_id), "NA", entity_id), "\n", sep="")
        cat("Entity_ID Length: ", ifelse(is.na(entity_id), 0, nchar(entity_id)), " (expected 28)\n", sep="")
        cat("\n")
        
        # Save failure report before terminating
        if(length(failed_samples) > 0) {
          failure_df <- bind_rows(failed_samples)
          write.table(failure_df, "FAILED_Entity_ID_extraction.tsv",
                      sep = "\t", row.names = FALSE, quote = FALSE)
          cat("Failure report saved to: FAILED_Entity_ID_extraction.tsv\n")
        }
        
        # Also save partial manifest for debugging
        if(length(all_metadata) > 0) {
          partial_manifest <- bind_rows(all_metadata)
          write.table(partial_manifest, "PARTIAL_manifest_before_failure.tsv",
                      sep = "\t", row.names = FALSE, quote = FALSE)
          cat("Partial manifest saved to: PARTIAL_manifest_before_failure.tsv\n")
        }
        
        # SIGTERM - stop execution with error
        stop(paste0(
          "\n\nSIGTERM: Entity_ID extraction failed!\n",
          "Sample: TCGA-", cancer, " / ", case_id, " / ", file_uuid, "\n",
          "Entity_ID from 'cases' column: '", ifelse(is.na(entity_id), "NA", entity_id), "'\n",
          "Expected 28-character barcode, got ", ifelse(is.na(entity_id), 0, nchar(entity_id)), " chars\n\n",
          "Please check the FAILED_Entity_ID_extraction.tsv file for details.\n"
        ))
      }
      
      # Determine tissue status
      tissue_status <- case_when(
        grepl("Primary Tumor|Recurrent Tumor|Metastatic", sample_type, ignore.case = TRUE) ~ "Tumor",
        grepl("Solid Tissue Normal|Blood Derived Normal", sample_type, ignore.case = TRUE) ~ "Normal",
        TRUE ~ "Other"
      )
      
      cancer_meta_list[[length(cancer_meta_list) + 1]] <- data.frame(
        UUID = file_uuid,
        Project = paste0("TCGA-", cancer),
        Case_ID = case_id,
        Entity_ID = entity_id,
        Sample_Type = sample_type,
        Cancer_Type = cancer,
        Tissue_Status = tissue_status,
        stringsAsFactors = FALSE
      )
    }
    
    if(length(cancer_meta_list) > 0) {
      cancer_df <- bind_rows(cancer_meta_list)
      all_metadata[[cancer]] <- cancer_df
      cat("Successfully processed", nrow(cancer_df), "samples\n")
      cat("All Entity_IDs are 28 characters:", all(nchar(cancer_df$Entity_ID) == 28), "\n")
    }
  }
  
  # Combine all metadata
  master_manifest <- bind_rows(all_metadata)
  
  # Final validation
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("MANIFEST VALIDATION\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Total samples:", nrow(master_manifest), "\n")
  cat("All Entity_IDs are 28 characters:", all(nchar(master_manifest$Entity_ID) == 28), "\n")
  cat("\nEntity_ID length distribution:\n")
  print(table(nchar(master_manifest$Entity_ID)))
  cat("\nSample Entity_IDs (first 5):\n")
  print(head(master_manifest$Entity_ID, 5))
  cat("\nSample Case_IDs (first 5):\n")
  print(head(master_manifest$Case_ID, 5))
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  # Save the manifest
  write.table(master_manifest,
              manifest_file,
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\nManifest saved to:", manifest_file, "\n")
  cat("Tumor samples:", sum(master_manifest$Tissue_Status == "Tumor"), "\n")
  cat("Normal samples:", sum(master_manifest$Tissue_Status == "Normal"), "\n\n")
}

# Now organize files into directories
cat("Starting file organization...\n\n")
base_path <- "GDCdata"

# Create organized_counts directory if it doesn't exist
if(!dir.exists("organized_counts")) {
  dir.create("organized_counts", recursive = TRUE)
}

files_organized <- 0
files_not_found <- 0

for(i in 1:nrow(master_manifest)) {
  uuid <- master_manifest$UUID[i]
  cancer <- master_manifest$Cancer_Type[i]
  status <- master_manifest$Tissue_Status[i]
  case_id <- master_manifest$Case_ID[i]
  entity_id <- master_manifest$Entity_ID[i]
  
  # Find the actual file
  file_pattern <- file.path(base_path, paste0("TCGA-", cancer),
                            "Transcriptome_Profiling",
                            "Gene_Expression_Quantification",
                            uuid, "*.tsv")
  
  source_file <- Sys.glob(file_pattern)
  
  if(length(source_file) == 0) {
    cat("Warning: File not found for UUID", uuid, "\n")
    files_not_found <- files_not_found + 1
    next
  }
  
  # Create organized directory structure
  target_dir <- file.path("organized_counts",
                          paste0("TCGA-", cancer),
                          status)
  
  if(!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  # Create new filename with Case_ID and UUID
  # Entity_ID is stored in manifest for lookup
  new_filename <- paste0(case_id, "_", uuid, ".rna_seq.augmented_star_gene_counts.tsv")
  target_file <- file.path(target_dir, new_filename)
  
  # Copy file
  file.copy(source_file, target_file, overwrite = TRUE)
  files_organized <- files_organized + 1
  
  if(i %% 100 == 0) {
    cat("Processed", i, "of", nrow(master_manifest), "files\n")
  }
}

cat("\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("FILE ORGANIZATION COMPLETE\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Total samples in manifest:", nrow(master_manifest), "\n")
cat("Files successfully organized:", files_organized, "\n")
cat("Files not found:", files_not_found, "\n")
cat("\nFiles organized in: organized_counts/\n")
cat("Manifest saved as:", manifest_file, "\n")
cat(paste(rep("=", 60), collapse=""), "\n")
