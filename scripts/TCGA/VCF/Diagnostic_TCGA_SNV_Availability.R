#!/usr/bin/env Rscript
# ============================================================
# Diagnostic_TCGA_SNV_Availability.R
#
# Purpose: Query GDC for all Simple Nucleotide Variation files
#          across all 33 TCGA cancer types. Reports available
#          variant callers, file counts, sizes, access levels,
#          and exact workflow.type strings needed for download.
#
# Environment: RNA-seq_NovoGene conda environment
# Usage: Rscript Diagnostic_TCGA_SNV_Availability.R
#
# Output:
#   TCGA_SNV_diagnostic_summary.tsv      — per cancer type x caller summary
#   TCGA_SNV_diagnostic_full_manifest.tsv — full file-level metadata
#   TCGA_SNV_diagnostic_report.txt       — human-readable report
# ============================================================

library(TCGAbiolinks)
library(dplyr)
library(data.table)

cat("============================================================\n")
cat("TCGA Simple Nucleotide Variation Diagnostic\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("TCGAbiolinks version:", as.character(packageVersion("TCGAbiolinks")), "\n\n")

# All 33 TCGA cancer types
cancer_types <- c("BRCA", "LUAD", "LUSC", "PRAD", "COAD", "STAD",
                  "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML",
                  "PAAD", "ESCA", "PCPG", "READ", "TGCT", "THYM",
                  "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS",
                  "CHOL", "GBM", "HNSC", "KIRC", "LGG", "OV",
                  "SKCM", "THCA", "UCEC")

# ============================================================
# PHASE 1: Discover available data.type and workflow.type values
# ============================================================
cat("============================================================\n")
cat("PHASE 1: Discovering available SNV data types and workflows\n")
cat("============================================================\n\n")

cat("Using TCGA-HNSC as test case to enumerate all available\n")
cat("data.type and workflow.type values under 'Simple Nucleotide Variation'...\n\n")

# First, query with just data.category to see all available data.types
tryCatch({
  test_query_broad <- GDCquery(
    project = "TCGA-HNSC",
    data.category = "Simple Nucleotide Variation"
  )
  test_meta_broad <- getResults(test_query_broad)
  
  cat("Available columns in SNV metadata:\n")
  for(i in seq_along(colnames(test_meta_broad))) {
    cat(sprintf("  %3d. %s\n", i, colnames(test_meta_broad)[i]))
  }
  
  cat("\n--- data_type values ---\n")
  print(table(test_meta_broad$data_type, useNA = "ifany"))
  
  cat("\n--- analysis_workflow_type values ---\n")
  if("analysis_workflow_type" %in% colnames(test_meta_broad)) {
    print(table(test_meta_broad$analysis_workflow_type, useNA = "ifany"))
  } else {
    cat("  Column 'analysis_workflow_type' not found.\n")
    cat("  Checking for workflow-related columns...\n")
    wf_cols <- grep("workflow|type|analysis", colnames(test_meta_broad), 
                    ignore.case = TRUE, value = TRUE)
    cat("  Workflow-related columns:", paste(wf_cols, collapse = ", "), "\n")
    for(col in wf_cols) {
      cat(sprintf("\n  Values in '%s':\n", col))
      print(table(test_meta_broad[[col]], useNA = "ifany"))
    }
  }
  
  cat("\n--- data_format values ---\n")
  print(table(test_meta_broad$data_format, useNA = "ifany"))
  
  cat("\n--- access values ---\n")
  print(table(test_meta_broad$access, useNA = "ifany"))
  
  cat("\n--- experimental_strategy values ---\n")
  if("experimental_strategy" %in% colnames(test_meta_broad)) {
    print(table(test_meta_broad$experimental_strategy, useNA = "ifany"))
  }
  
  # Cross-tabulate data_type x workflow
  cat("\n--- Cross-tabulation: data_type x workflow ---\n")
  if("analysis_workflow_type" %in% colnames(test_meta_broad)) {
    print(table(test_meta_broad$data_type, 
                test_meta_broad$analysis_workflow_type, useNA = "ifany"))
  }
  
  # Show sample filenames for each data_type to understand naming conventions
  cat("\n--- Sample filenames by data_type ---\n")
  for(dt in unique(test_meta_broad$data_type)) {
    cat(sprintf("\ndata_type = '%s':\n", dt))
    subset_files <- test_meta_broad %>% 
      filter(data_type == dt) %>%
      head(3)
    for(fn in subset_files$file_name) {
      cat("  ", fn, "\n")
    }
  }
  
  # Check for Entity_ID / barcode columns
  cat("\n--- Barcode / sample ID columns ---\n")
  id_cols <- grep("sample|aliquot|submitter|barcode|cases|entity",
                  colnames(test_meta_broad), ignore.case = TRUE, value = TRUE)
  cat("Potentially useful ID columns:", paste(id_cols, collapse = ", "), "\n")
  for(col in id_cols) {
    vals <- as.character(test_meta_broad[[col]][1:min(3, nrow(test_meta_broad))])
    cat(sprintf("  %s: '%s' (%d chars)\n", col, vals[1], nchar(vals[1])))
  }
  
}, error = function(e) {
  cat("ERROR in Phase 1:", e$message, "\n")
  cat("Continuing to Phase 2 with best-guess parameters...\n\n")
})

# ============================================================
# PHASE 2: Query all 33 cancer types for SNV data
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 2: Querying all 33 TCGA cancer types for SNV files\n")
cat("============================================================\n\n")

all_metadata <- list()
query_errors <- list()

for(cancer in cancer_types) {
  cat(sprintf("Querying TCGA-%s... ", cancer))
  
  tryCatch({
    query <- GDCquery(
      project = paste0("TCGA-", cancer),
      data.category = "Simple Nucleotide Variation"
    )
    
    metadata <- getResults(query)
    metadata$Cancer_Type <- cancer
    all_metadata[[cancer]] <- metadata
    
    cat(nrow(metadata), "files found\n")
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    query_errors[[cancer]] <<- e$message
  })
}

# Combine all metadata
cat("\nCombining metadata from all cancer types...\n")
full_manifest <- bind_rows(all_metadata)
cat("Total files across all cancer types:", nrow(full_manifest), "\n\n")

# ============================================================
# PHASE 3: Build summary tables
# ============================================================
cat("============================================================\n")
cat("PHASE 3: Summary statistics\n")
cat("============================================================\n\n")

# --- Overall summary ---
cat("--- Overall file counts by data_type ---\n")
print(table(full_manifest$data_type, useNA = "ifany"))

cat("\n--- Overall file counts by workflow ---\n")
if("analysis_workflow_type" %in% colnames(full_manifest)) {
  print(table(full_manifest$analysis_workflow_type, useNA = "ifany"))
  wf_col <- "analysis_workflow_type"
} else {
  # Try alternatives
  wf_col <- grep("workflow", colnames(full_manifest), ignore.case = TRUE, value = TRUE)[1]
  if(!is.null(wf_col) && length(wf_col) > 0) {
    cat("Using column:", wf_col, "\n")
    print(table(full_manifest[[wf_col]], useNA = "ifany"))
  } else {
    wf_col <- "data_type"
    cat("No workflow column found, using data_type as grouping\n")
  }
}

cat("\n--- Overall file counts by access level ---\n")
print(table(full_manifest$access, useNA = "ifany"))

cat("\n--- Cross-tabulation: data_type x access ---\n")
print(table(full_manifest$data_type, full_manifest$access, useNA = "ifany"))

# --- Per cancer type x caller summary ---
cat("\n--- Files per cancer type x workflow ---\n")
summary_table <- full_manifest %>%
  group_by(Cancer_Type, data_type) %>%
  summarize(
    n_files = n(),
    total_size_MB = round(sum(file_size, na.rm = TRUE) / 1e6, 2),
    access_types = paste(unique(access), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(Cancer_Type, data_type)

print(as.data.frame(summary_table), row.names = FALSE)

# --- MuTect2-specific summary (our target) ---
cat("\n\n--- MuTect2-specific files ---\n")
mutect2_files <- full_manifest %>%
  filter(grepl("mutect2|MuTect2", file_name, ignore.case = TRUE) |
         grepl("mutect2|MuTect2", get(wf_col), ignore.case = TRUE))

if(nrow(mutect2_files) > 0) {
  cat("Total MuTect2 files:", nrow(mutect2_files), "\n\n")
  
  mutect2_summary <- mutect2_files %>%
    group_by(Cancer_Type, data_type) %>%
    summarize(
      n_files = n(),
      total_size_MB = round(sum(file_size, na.rm = TRUE) / 1e6, 2),
      access = paste(unique(access), collapse = ", "),
      .groups = "drop"
    )
  
  cat("MuTect2 files per cancer type:\n")
  print(as.data.frame(mutect2_summary), row.names = FALSE)
  
  cat("\nMuTect2 file name patterns (first 5):\n")
  for(fn in head(unique(mutect2_files$file_name), 5)) {
    cat("  ", fn, "\n")
  }
  
  # Check for annotated vs raw
  cat("\nMuTect2 annotated vs raw:\n")
  mutect2_files$annotated <- grepl("annotation|annotated", mutect2_files$file_name, 
                                    ignore.case = TRUE)
  print(table(Annotated = mutect2_files$annotated))
  
  # Total download size
  cat("\nTotal MuTect2 download size:", 
      round(sum(mutect2_files$file_size, na.rm = TRUE) / 1e9, 2), "GB\n")
  
  # Annotated only size
  annotated_size <- sum(mutect2_files$file_size[mutect2_files$annotated], na.rm = TRUE)
  cat("Annotated-only download size:", round(annotated_size / 1e9, 2), "GB\n")
  
} else {
  cat("No MuTect2 files identified by filename pattern.\n")
  cat("Check workflow column values above for the correct filter string.\n")
}

# --- Entity_ID / barcode mapping diagnostic ---
cat("\n\n--- Sample-to-patient mapping columns ---\n")
id_cols <- grep("sample|aliquot|submitter|barcode|cases|entity",
                colnames(full_manifest), ignore.case = TRUE, value = TRUE)
cat("Available ID columns:", paste(id_cols, collapse = ", "), "\n")

for(col in id_cols) {
  vals <- as.character(full_manifest[[col]])
  vals <- vals[!is.na(vals) & vals != ""]
  if(length(vals) > 0) {
    cat(sprintf("\n  %s:\n", col))
    cat("    Example:", head(vals, 2), "\n")
    cat("    Char lengths:", paste(unique(head(nchar(vals), 10)), collapse = ", "), "\n")
    n_28 <- sum(nchar(vals) == 28)
    if(n_28 > 0) cat("    >>> 28-char barcodes:", n_28, "/", length(vals), "\n")
  }
}

# ============================================================
# PHASE 4: Pindel summary (for XRCC4 indel analysis)
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 4: Pindel files (for XRCC4 indel analysis)\n")
cat("============================================================\n\n")

pindel_files <- full_manifest %>%
  filter(grepl("pindel|Pindel", file_name, ignore.case = TRUE) |
         grepl("pindel|Pindel", get(wf_col), ignore.case = TRUE))

if(nrow(pindel_files) > 0) {
  cat("Total Pindel files:", nrow(pindel_files), "\n")
  
  pindel_summary <- pindel_files %>%
    group_by(Cancer_Type) %>%
    summarize(
      n_files = n(),
      total_size_MB = round(sum(file_size, na.rm = TRUE) / 1e6, 2),
      .groups = "drop"
    )
  
  cat("Pindel files per cancer type:\n")
  print(as.data.frame(pindel_summary), row.names = FALSE)
  cat("\nTotal Pindel download size:", 
      round(sum(pindel_files$file_size, na.rm = TRUE) / 1e9, 2), "GB\n")
} else {
  cat("No Pindel files identified.\n")
}

# ============================================================
# PHASE 5: Save outputs
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 5: Saving output files\n")
cat("============================================================\n\n")

# Full manifest
manifest_file <- "/master/jlehle/SHARED/TCGA/VCF/TCGA_SNV_diagnostic_full_manifest.tsv"
write.table(full_manifest, manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", manifest_file, "(", nrow(full_manifest), "rows )\n")

# Summary table
summary_file <- "/master/jlehle/SHARED/TCGA/VCF/TCGA_SNV_diagnostic_summary.tsv"
write.table(summary_table, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", summary_file, "\n")

# MuTect2-specific manifest if available
if(nrow(mutect2_files) > 0) {
  mutect2_file <- "/master/jlehle/SHARED/TCGA/VCF/TCGA_SNV_diagnostic_MuTect2_manifest.tsv"
  write.table(mutect2_files, mutect2_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved:", mutect2_file, "(", nrow(mutect2_files), "rows )\n")
}

# ============================================================
# PHASE 6: Report — exact GDCquery parameters for download
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 6: Recommended GDCquery parameters\n")
cat("============================================================\n\n")

cat("Based on the diagnostic results, use these parameters for the download script:\n\n")

if("analysis_workflow_type" %in% colnames(full_manifest)) {
  # Get unique workflow types that contain MuTect2
  mt2_workflows <- unique(full_manifest[[wf_col]][
    grepl("mutect2|MuTect2", full_manifest[[wf_col]], ignore.case = TRUE)
  ])
  
  if(length(mt2_workflows) > 0) {
    cat("MuTect2 annotated VCFs:\n")
    cat("  GDCquery(\n")
    cat("    project = 'TCGA-{CANCER}',\n")
    cat("    data.category = 'Simple Nucleotide Variation',\n")
    
    # Find the annotated data_type
    annotated_types <- unique(mutect2_files$data_type[mutect2_files$annotated])
    if(length(annotated_types) > 0) {
      cat(sprintf("    data.type = '%s',\n", annotated_types[1]))
    }
    
    cat(sprintf("    workflow.type = '%s'\n", mt2_workflows[1]))
    cat("  )\n\n")
  }
  
  # Pindel parameters
  pindel_workflows <- unique(full_manifest[[wf_col]][
    grepl("pindel|Pindel", full_manifest[[wf_col]], ignore.case = TRUE)
  ])
  
  if(length(pindel_workflows) > 0) {
    cat("Pindel VCFs (for XRCC4 indels):\n")
    cat("  GDCquery(\n")
    cat("    project = 'TCGA-{CANCER}',\n")
    cat("    data.category = 'Simple Nucleotide Variation',\n")
    
    annotated_pindel_types <- unique(pindel_files$data_type[
      grepl("annotation|annotated", pindel_files$file_name, ignore.case = TRUE)
    ])
    if(length(annotated_pindel_types) > 0) {
      cat(sprintf("    data.type = '%s',\n", annotated_pindel_types[1]))
    }
    
    cat(sprintf("    workflow.type = '%s'\n", pindel_workflows[1]))
    cat("  )\n")
  }
}

# ============================================================
# Summary report
# ============================================================
cat("\n\n============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("============================================================\n")
cat("Total cancer types queried:", length(all_metadata), "/", length(cancer_types), "\n")
cat("Failed queries:", length(query_errors), "\n")
if(length(query_errors) > 0) {
  for(cancer in names(query_errors)) {
    cat("  TCGA-", cancer, ":", query_errors[[cancer]], "\n", sep = "")
  }
}
cat("Total SNV files found:", nrow(full_manifest), "\n")
cat("MuTect2 files:", nrow(mutect2_files), "\n")
cat("Pindel files:", nrow(pindel_files), "\n")
cat("\nOutput files:\n")
cat("  1.", manifest_file, "\n")
cat("  2.", summary_file, "\n")
if(nrow(mutect2_files) > 0) cat("  3. TCGA_SNV_diagnostic_MuTect2_manifest.tsv\n")
cat("============================================================\n")
