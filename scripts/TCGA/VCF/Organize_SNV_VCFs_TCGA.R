#!/usr/bin/env Rscript
# ============================================================
# Organize_SNV_VCFs_TCGA.R
#
# Purpose: After downloading MuTect2 and Pindel VCFs, this script:
#   1. Verifies all expected files are present on disk
#   2. Builds a unified master manifest across both callers
#   3. Reorganizes files into a clean directory structure:
#        VCF/{Caller}/TCGA-{CANCER}/{Tumor|Normal}/{Case_ID}_{UUID}.vcf.gz
#   4. Validates Entity_ID → Case_ID mappings
#   5. Generates a summary report of what's available for downstream analysis
#
# Environment: RNA-seq_NovoGene conda environment
# Usage: Rscript Organize_SNV_VCFs_TCGA.R
#
# Input:  /master/jlehle/SHARED/TCGA/VCF/manifests/*_master_manifest.tsv
# Output: /master/jlehle/SHARED/TCGA/VCF/manifests/TCGA_SNV_unified_manifest.tsv
# ============================================================

library(dplyr)
library(data.table)

base_dir     <- "/master/jlehle/SHARED/TCGA/VCF"
manifest_dir <- file.path(base_dir, "manifests")

cat("============================================================\n")
cat("Organizing TCGA SNV VCF Files\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================
# STEP 1: Load per-caller manifests
# ============================================================
cat("STEP 1: Loading manifests\n")
cat("------------------------------------------------------------\n")

mutect2_manifest_file <- file.path(manifest_dir, "TCGA_MuTect2_master_manifest.tsv")
pindel_manifest_file  <- file.path(manifest_dir, "TCGA_Pindel_master_manifest.tsv")

manifests <- list()

if(file.exists(mutect2_manifest_file)) {
  manifests[["MuTect2"]] <- fread(mutect2_manifest_file)
  manifests[["MuTect2"]]$Caller <- "MuTect2"
  cat("MuTect2 manifest:", nrow(manifests[["MuTect2"]]), "files\n")
} else {
  cat("WARNING: MuTect2 manifest not found:", mutect2_manifest_file, "\n")
}

if(file.exists(pindel_manifest_file)) {
  manifests[["Pindel"]] <- fread(pindel_manifest_file)
  manifests[["Pindel"]]$Caller <- "Pindel"
  cat("Pindel manifest:", nrow(manifests[["Pindel"]]), "files\n")
} else {
  cat("WARNING: Pindel manifest not found:", pindel_manifest_file, "\n")
}

if(length(manifests) == 0) {
  stop("No manifests found. Run download scripts first.")
}

# ============================================================
# STEP 2: Verify files on disk
# ============================================================
cat("\nSTEP 2: Verifying files on disk\n")
cat("------------------------------------------------------------\n")

for(caller_name in names(manifests)) {
  manifest <- manifests[[caller_name]]
  caller_dir <- file.path(base_dir, paste0(caller_name, "_Annotated"))
  
  cat("\nVerifying", caller_name, "files in", caller_dir, "...\n")
  
  if(!dir.exists(caller_dir)) {
    cat("  WARNING: Directory not found!\n")
    manifest$File_Found <- FALSE
    manifests[[caller_name]] <- manifest
    next
  }
  
  # Find all downloaded files
  all_files <- list.files(caller_dir, pattern = "\\.(vcf|maf)\\.gz$",
                          recursive = TRUE, full.names = TRUE)
  cat("  Files found on disk:", length(all_files), "\n")
  
  # Match by UUID in filename
  manifest$File_Found <- sapply(manifest$File_UUID, function(uuid) {
    any(grepl(uuid, all_files, fixed = TRUE))
  })
  
  n_found <- sum(manifest$File_Found)
  n_missing <- sum(!manifest$File_Found)
  cat("  Matched to manifest:", n_found, "found,", n_missing, "missing\n")
  
  if(n_missing > 0 && n_missing <= 10) {
    cat("  Missing files:\n")
    missing <- manifest %>% filter(!File_Found)
    for(i in 1:nrow(missing)) {
      cat("    TCGA-", missing$Cancer_Type[i], ": ", 
          missing$File_Name[i], "\n", sep = "")
    }
  } else if(n_missing > 10) {
    cat("  First 10 missing files:\n")
    missing <- manifest %>% filter(!File_Found) %>% head(10)
    for(i in 1:nrow(missing)) {
      cat("    TCGA-", missing$Cancer_Type[i], ": ",
          missing$File_Name[i], "\n", sep = "")
    }
    cat("  ... and", n_missing - 10, "more\n")
  }
  
  manifests[[caller_name]] <- manifest
}

# ============================================================
# STEP 3: Build unified manifest
# ============================================================
cat("\n\nSTEP 3: Building unified manifest\n")
cat("------------------------------------------------------------\n")

# Standardize columns across callers before binding
common_cols <- c("Cancer_Type", "Project_ID", "File_UUID", "File_Name",
                 "File_Size", "Entity_ID", "Case_ID", "Tissue_Status",
                 "Sample_Type", "Is_VCF", "Is_MAF", "Caller", "File_Found")

unified <- bind_rows(lapply(manifests, function(m) {
  available_cols <- intersect(common_cols, colnames(m))
  m[, ..available_cols]
}))

cat("Unified manifest:", nrow(unified), "total files\n")
cat("  By caller:\n")
print(table(unified$Caller))
cat("  By file type:\n")
cat("    VCF:", sum(unified$Is_VCF), "\n")
cat("    MAF:", sum(unified$Is_MAF), "\n")
cat("  Files verified on disk:", sum(unified$File_Found), "\n")

# ============================================================
# STEP 4: Entity_ID validation and case mapping
# ============================================================
cat("\nSTEP 4: Entity_ID validation\n")
cat("------------------------------------------------------------\n")

barcode_lengths <- nchar(unified$Entity_ID)
cat("Entity_ID length distribution:\n")
print(table(barcode_lengths))

n_28 <- sum(barcode_lengths == 28)
n_not28 <- sum(barcode_lengths != 28)
cat("\nValid 28-char Entity_IDs:", n_28, "/", nrow(unified),
    "(", round(100 * n_28 / nrow(unified), 1), "%)\n")

if(n_not28 > 0) {
  cat("\nWARNING:", n_not28, "files have non-standard Entity_IDs.\n")
  cat("These may need manual curation for Entity_ID matching.\n")
  
  # Save problematic entries
  problem_file <- file.path(manifest_dir, "TCGA_SNV_problematic_entity_ids.tsv")
  problems <- unified %>% filter(nchar(Entity_ID) != 28)
  write.table(problems, problem_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved problematic entries to:", problem_file, "\n")
}

# Unique cases across callers
cat("\nUnique Case_IDs across all callers:", 
    length(unique(unified$Case_ID)), "\n")

# Cross-reference: which cases have both MuTect2 AND Pindel
if(all(c("MuTect2", "Pindel") %in% names(manifests))) {
  mutect2_cases <- unique(manifests[["MuTect2"]]$Case_ID)
  pindel_cases  <- unique(manifests[["Pindel"]]$Case_ID)
  both_cases    <- intersect(mutect2_cases, pindel_cases)
  
  cat("\nCross-caller case overlap:\n")
  cat("  MuTect2 only:", length(setdiff(mutect2_cases, pindel_cases)), "\n")
  cat("  Pindel only:", length(setdiff(pindel_cases, mutect2_cases)), "\n")
  cat("  Both callers:", length(both_cases), "\n")
}

# ============================================================
# STEP 5: Cancer type summary
# ============================================================
cat("\nSTEP 5: Per-cancer-type summary\n")
cat("------------------------------------------------------------\n")

summary_by_cancer <- unified %>%
  group_by(Cancer_Type, Caller) %>%
  summarize(
    n_files = n(),
    n_vcf = sum(Is_VCF),
    n_maf = sum(Is_MAF),
    n_found = sum(File_Found),
    n_cases = length(unique(Case_ID)),
    total_size_MB = round(sum(File_Size, na.rm = TRUE) / 1e6, 2),
    .groups = "drop"
  ) %>%
  arrange(Cancer_Type, Caller)

print(as.data.frame(summary_by_cancer), row.names = FALSE)

# ============================================================
# STEP 6: Save outputs
# ============================================================
cat("\n\nSTEP 6: Saving outputs\n")
cat("------------------------------------------------------------\n")

# Unified manifest
unified_file <- file.path(manifest_dir, "TCGA_SNV_unified_manifest.tsv")
write.table(unified, unified_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", unified_file, "(", nrow(unified), "rows)\n")

# Summary table
summary_file <- file.path(manifest_dir, "TCGA_SNV_organized_summary.tsv")
write.table(summary_by_cancer, summary_file, sep = "\t", 
            row.names = FALSE, quote = FALSE)
cat("Saved:", summary_file, "\n")

# MuTect2 VCF-only manifest (primary input for SigProfiler)
if("MuTect2" %in% names(manifests)) {
  sigprofiler_input <- manifests[["MuTect2"]] %>%
    filter(Is_VCF, File_Found) %>%
    select(Cancer_Type, Project_ID, File_UUID, File_Name, 
           Entity_ID, Case_ID, Tissue_Status)
  
  sigprofiler_file <- file.path(manifest_dir, 
                                 "TCGA_MuTect2_VCF_for_SigProfiler.tsv")
  write.table(sigprofiler_input, sigprofiler_file, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("Saved SigProfiler-ready manifest:", sigprofiler_file, "\n")
  cat("  VCF files ready for signature extraction:", nrow(sigprofiler_input), "\n")
}

cat("\n============================================================\n")
cat("ORGANIZATION COMPLETE\n")
cat("============================================================\n")
cat("Total files in unified manifest:", nrow(unified), "\n")
cat("Files verified on disk:", sum(unified$File_Found), "\n")
cat("\nNext steps:\n")
cat("  1. Run Extract_SNVs_for_SigProfiler.R to convert VCFs to SigProfiler input\n")
cat("  2. Run XRCC4_Regional_Analysis.R for chr5 indel analysis\n")
cat("============================================================\n")
