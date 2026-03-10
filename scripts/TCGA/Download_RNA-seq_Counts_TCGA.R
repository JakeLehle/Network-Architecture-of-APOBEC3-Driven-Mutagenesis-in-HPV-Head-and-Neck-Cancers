library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

# List of all TCGA cancer types
cancer_types <- c("BRCA", "LUAD", "LUSC", "PRAD", "COAD", "STAD", 
                  "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML",
                  "PAAD", "ESCA", "PCPG", "READ", "TGCT", "THYM",
                  "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS",
                  "CHOL", "GBM", "HNSC", "KIRC", "LGG", "OV",
                  "SKCM", "THCA", "UCEC")

# Function to check if cancer type already processed successfully
is_processed <- function(cancer) {
  rds_file <- paste0("TCGA_", cancer, "_counts.rds")
  return(file.exists(rds_file))
}

# Function to check if data directory exists and has files
has_data_directory <- function(cancer) {
  data_dir <- file.path("GDCdata", paste0("TCGA-", cancer), 
                       "Transcriptome_Profiling", 
                       "Gene_Expression_Quantification")
  if(!dir.exists(data_dir)) {
    return(FALSE)
  }
  # Check if there are actual TSV files
  tsv_files <- list.files(data_dir, pattern = "\\.tsv$", recursive = TRUE)
  return(length(tsv_files) > 0)
}

# Clean up corrupted tar.gz files (these cause problems)
clean_corrupted_archives <- function() {
  tar_files <- list.files(pattern = "\\.tar\\.gz$", full.names = TRUE)
  if(length(tar_files) > 0) {
    cat("Removing", length(tar_files), "potentially corrupted tar.gz file(s)\n")
    file.remove(tar_files)
  }
}

# Download with retry logic
download_with_retry <- function(cancer, max_retries = 3) {
  for(attempt in 1:max_retries) {
    tryCatch({
      cat("\n=== Processing TCGA-", cancer, " (Attempt ", attempt, "/", max_retries, ") ===\n", sep="")
      
      # Query gene expression data
      query <- GDCquery(
        project = paste0("TCGA-", cancer),
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
      )
      
      n_samples <- nrow(getResults(query))
      cat("Found", n_samples, "samples for TCGA-", cancer, "\n", sep="")
      
      # Download
      cat("Downloading data for project TCGA-", cancer, "\n", sep="")
      GDCdownload(query, method = "api")
      
      # Clean up any tar.gz files immediately after download
      clean_corrupted_archives()
      
      # Verify download before preparing
      if(!has_data_directory(cancer)) {
        stop("Download appeared to complete but files not found")
      }
      
      # Prepare data matrix
      cat("Preparing data matrix for TCGA-", cancer, "\n", sep="")
      data <- GDCprepare(query)
      
      # Save
      rds_file <- paste0("TCGA_", cancer, "_counts.rds")
      saveRDS(data, file = rds_file)
      cat("Successfully saved:", rds_file, "\n")
      
      return(TRUE)
      
    }, error = function(e) {
      cat("ERROR on attempt", attempt, "for TCGA-", cancer, ":", e$message, "\n")
      
      # Clean up corrupted archives before retry
      clean_corrupted_archives()
      
      if(attempt < max_retries) {
        wait_time <- 30 * attempt  # Increasing wait time
        cat("  Waiting", wait_time, "seconds before retry...\n")
        Sys.sleep(wait_time)
      } else {
        cat("FAILED: TCGA-", cancer, " after", max_retries, "attempts\n", sep="")
        return(FALSE)
      }
    })
  }
  return(FALSE)
}

# Main processing loop
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Starting TCGA data download with resume capability\n")
cat("Total cancer types:", length(cancer_types), "\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Clean up any existing corrupted archives first
clean_corrupted_archives()

failed_types <- c()
skipped_types <- c()
recovered_types <- c()

for(cancer in cancer_types) {
  # Check if already processed successfully
  if(is_processed(cancer)) {
    cat("SKIPPING TCGA-", cancer, " - already processed (RDS file exists)\n", sep="")
    skipped_types <- c(skipped_types, cancer)
    next
  }
  
  # Check if data directory exists but RDS doesn't (recovery scenario)
  if(has_data_directory(cancer)) {
    cat("TCGA-", cancer, " - data directory exists, attempting to prepare without re-download...\n", sep="")
    
    tryCatch({
      query <- GDCquery(
        project = paste0("TCGA-", cancer),
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
      )
      
      data <- GDCprepare(query)
      rds_file <- paste0("TCGA_", cancer, "_counts.rds")
      saveRDS(data, file = rds_file)
      cat("Successfully recovered and saved:", rds_file, "\n")
      recovered_types <- c(recovered_types, cancer)
      next
      
    }, error = function(e) {
      cat("  Could not prepare existing data, will re-download. Error:", e$message, "\n")
    })
  }
  
  # Download with retry
  success <- download_with_retry(cancer)
  
  if(!success) {
    failed_types <- c(failed_types, cancer)
  }
}

# Summary
cat("\n\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("DOWNLOAD SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Total cancer types:       ", length(cancer_types), "\n")
cat("Skipped (already done):   ", length(skipped_types), "\n")
cat("Recovered (had data):     ", length(recovered_types), "\n")
cat("Newly downloaded:         ", length(cancer_types) - length(failed_types) - length(skipped_types) - length(recovered_types), "\n")
cat("Failed:                   ", length(failed_types), "\n")
cat(paste(rep("=", 60), collapse=""), "\n")

if(length(skipped_types) > 0) {
  cat("\nAlready completed:\n")
  cat(paste(skipped_types, collapse=", "), "\n")
}

if(length(recovered_types) > 0) {
  cat("\nRecovered from existing data:\n")
  cat(paste(recovered_types, collapse=", "), "\n")
}

if(length(failed_types) > 0) {
  cat("\nFailed cancer types:\n")
  cat(paste(failed_types, collapse=", "), "\n")
  cat("\nTo retry only failed types, run:\n")
  cat("cancer_types <- c('", paste(failed_types, collapse="', '"), "')\n", sep="")
  cat("# Then re-run the main loop\n")
}

# Save processing report
report <- data.frame(
  Cancer_Type = cancer_types,
  Status = sapply(cancer_types, function(c) {
    if(c %in% skipped_types) return("Already_Done")
    if(c %in% recovered_types) return("Recovered")
    if(c %in% failed_types) return("Failed")
    return("Newly_Downloaded")
  }),
  RDS_Exists = sapply(cancer_types, is_processed),
  Data_Dir_Exists = sapply(cancer_types, has_data_directory)
)

write.table(report, "download_status_report.tsv", 
            sep="\t", row.names=FALSE, quote=FALSE)
cat("\nStatus report saved to: download_status_report.tsv\n")
