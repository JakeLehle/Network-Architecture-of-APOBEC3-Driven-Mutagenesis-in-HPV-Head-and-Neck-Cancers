#!/usr/bin/env Rscript
# ============================================================
# Download_cBioPortal_Expression_TCGA.R
#
# Purpose: Download TCGA RNA-seq expression data from cBioPortal
#          via the public REST API for all 33 TCGA cancer types.
#
#          cBioPortal hosts the LEGACY GRCh37/hg19-aligned TCGA
#          RNA-seq data processed through the RSEM pipeline, which
#          is no longer available on the GDC portal (now GRCh38 only).
#          This is the most likely source of Diako's original
#          expression file (A3_Table_All_Samples_TCGA_Mohadeseh.tsv).
#
# Strategy: Uses the cBioPortal REST API (https://www.cbioportal.org/api/)
#           to query each TCGA study's RNA-seq V2 molecular profile
#           and fetch expression values for the 7 APOBEC3 genes.
#           Also fetches all-gene expression for a full comparison.
#
# Dependencies: httr, jsonlite, dplyr, data.table
#
# Output:
#   cBioPortal_A3_Expression_All_TCGA.tsv  — A3 genes only, all tumors
#   cBioPortal_Full_Expression_All_TCGA.tsv — all genes, all tumors (large)
#   cBioPortal_download_report.txt          — QC report
#
# ============================================================
# NORMALIZATION BACKGROUND
# ============================================================
#
# cBioPortal TCGA RNA-seq data uses RSEM (RNA-Seq by Expectation
# Maximization) normalization from the LEGACY Firehose pipeline,
# aligned to GRCh37/hg19. This is fundamentally different from the
# GRCh38/STAR-Counts pipeline on the current GDC portal.
#
# The key differences in normalization methods:
#
# 1. RSEM (cBioPortal / legacy GRCh37)
#    -----------------------------------------------------------------
#    RSEM models multi-mapped reads probabilistically using an EM
#    algorithm. The output "RSEM normalized count" (often called
#    "RNA-Seq V2 RSEM" on cBioPortal) is the estimated fraction of
#    reads from a gene, scaled by 10^6 (i.e., reads per million,
#    but accounting for multi-mapping uncertainty). In TCGA Firehose
#    data, these are further upper-quartile normalized across samples
#    within each cancer type. The values on cBioPortal are typically
#    the RSEM "scaled estimate" × 10^6 or the raw RSEM expected
#    counts, depending on the specific molecular profile.
#
#    cBioPortal exposes two RNA-seq profiles per TCGA study:
#      a) {study}_rna_seq_v2_mrna           — RSEM normalized (z-scores
#         or direct values; this is the log2-transformed or raw RSEM)
#      b) {study}_rna_seq_v2_mrna_median_Zscores — z-score normalized
#         relative to all samples in the study
#
#    We want the RAW expression values, not z-scores. The profile ID
#    pattern for raw RSEM is: {study_id}_rna_seq_v2_mrna
#
# 2. FPKM-UQ (GDC / current GRCh38)
#    -----------------------------------------------------------------
#    FPKM-UQ = Fragments Per Kilobase of transcript per Million mapped
#    reads, with the denominator replaced by the UPPER QUARTILE of
#    non-zero gene counts × 10^3. This replaces the total library size
#    with the 75th percentile gene count, making it more robust to
#    highly expressed outlier genes (e.g., ribosomal RNA, hemoglobin).
#
#    Formula:
#      FPKM-UQ = (raw_count × 10^9) / (gene_length × 75th_percentile_count)
#
#    This is computed from STAR 2-pass alignment to GRCh38.
#
# 3. TPM (also on GDC)
#    -----------------------------------------------------------------
#    TPM = Transcripts Per Million. Each gene's count is divided by its
#    effective length (from STAR), then the resulting rates are
#    normalized so they sum to 10^6 across all genes. TPM is comparable
#    between samples because the denominator is sample-specific.
#
# KEY COMPARISON IMPLICATIONS:
# -----------------------------------------------------------------
# - RSEM and FPKM-UQ are NOT directly comparable in absolute magnitude
#   because they use different normalization denominators and different
#   reference genomes (GRCh37 vs GRCh38 gene models).
# - However, RSEM and FPKM-UQ are highly RANK-correlated within a
#   sample (Spearman rho typically > 0.95) because both are proportional
#   to true expression.
# - If the mystery file came from cBioPortal, we expect:
#     * High correlation with RSEM values (r > 0.999)
#     * Moderate correlation with FPKM-UQ values (r ~ 0.85-0.95)
#       due to different reference genome + normalization
# - The correlation with FPKM-UQ should be LOWER than with RSEM
#   because FPKM-UQ and RSEM differ in both:
#     a) Reference genome (GRCh38 vs GRCh37 → different gene models,
#        coordinates, and effective lengths)
#     b) Normalization method (upper-quartile vs EM-based estimation)
#
# If the mystery file DOES match RSEM but NOT FPKM-UQ, this confirms
# it was sourced from cBioPortal or the legacy Firehose pipeline,
# and explains why it doesn't match our newly downloaded GRCh38 data.
# ============================================================

library(httr)
library(jsonlite)
library(dplyr)
library(data.table)

# ============================================================
# CONFIGURATION
# ============================================================

# cBioPortal API base URL
CBIO_API <- "https://www.cbioportal.org/api"

# A3 gene symbols (primary targets)
A3_GENES <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
              "APOBEC3F", "APOBEC3G", "APOBEC3H")

# Entrez gene IDs for A3 genes (needed for some API endpoints)
# These are stable identifiers that cBioPortal uses internally
A3_ENTREZ <- c(
  APOBEC3A = 200315,
  APOBEC3B = 9582,
  APOBEC3C = 27350,
  APOBEC3D = 140564,
  APOBEC3F = 200316,
  APOBEC3G = 60489,
  APOBEC3H = 164668
)

# All 32 Pan-Cancer Atlas 2018 TCGA study IDs on cBioPortal
# NOTE: COAD and READ do NOT have separate PCA studies — they are
#       combined as "coadread_tcga_pan_can_atlas_2018" (594 samples).
#       Confirmed via Diagnostic_cBioPortal_Study_IDs.R (2026-04-15).
#       The cancer_type column will report "COADREAD" for these samples;
#       downstream scripts can split them back using barcode matching
#       against our GDC metadata if needed.
TCGA_STUDIES <- c(
  "acc_tcga_pan_can_atlas_2018",
  "blca_tcga_pan_can_atlas_2018",
  "brca_tcga_pan_can_atlas_2018",
  "cesc_tcga_pan_can_atlas_2018",
  "chol_tcga_pan_can_atlas_2018",
  "coadread_tcga_pan_can_atlas_2018",  # COAD + READ combined (594 samples)
  "dlbc_tcga_pan_can_atlas_2018",
  "esca_tcga_pan_can_atlas_2018",
  "gbm_tcga_pan_can_atlas_2018",
  "hnsc_tcga_pan_can_atlas_2018",
  "kich_tcga_pan_can_atlas_2018",
  "kirc_tcga_pan_can_atlas_2018",
  "kirp_tcga_pan_can_atlas_2018",
  "laml_tcga_pan_can_atlas_2018",
  "lgg_tcga_pan_can_atlas_2018",
  "lihc_tcga_pan_can_atlas_2018",
  "luad_tcga_pan_can_atlas_2018",
  "lusc_tcga_pan_can_atlas_2018",
  "meso_tcga_pan_can_atlas_2018",
  "ov_tcga_pan_can_atlas_2018",
  "paad_tcga_pan_can_atlas_2018",
  "pcpg_tcga_pan_can_atlas_2018",
  "prad_tcga_pan_can_atlas_2018",
  "sarc_tcga_pan_can_atlas_2018",
  "skcm_tcga_pan_can_atlas_2018",
  "stad_tcga_pan_can_atlas_2018",
  "tgct_tcga_pan_can_atlas_2018",
  "thca_tcga_pan_can_atlas_2018",
  "thym_tcga_pan_can_atlas_2018",
  "ucec_tcga_pan_can_atlas_2018",
  "ucs_tcga_pan_can_atlas_2018",
  "uvm_tcga_pan_can_atlas_2018"
)

# Output directory (adjust for your HPC paths)
# OUTPUT_DIR <- "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1/TROUBLESHOOTING"
OUTPUT_DIR <- getwd()  # Default to working directory; change as needed

# Rate limiting
SLEEP_BETWEEN_REQUESTS <- 0.5  # seconds, be polite to cBioPortal

# ============================================================
# HELPER FUNCTIONS
# ============================================================

#' Query the cBioPortal REST API
#'
#' @param endpoint API endpoint (appended to CBIO_API)
#' @param method HTTP method ("GET" or "POST")
#' @param body POST body (list or JSON string)
#' @param max_retries Number of retry attempts on failure
#' @return Parsed JSON response (list or data.frame)
cbio_query <- function(endpoint, method = "GET", body = NULL, max_retries = 3) {
  url <- paste0(CBIO_API, endpoint)
  
  for(attempt in 1:max_retries) {
    tryCatch({
      if(method == "GET") {
        response <- GET(url, add_headers("Accept" = "application/json"))
      } else if(method == "POST") {
        response <- POST(
          url,
          add_headers(
            "Accept" = "application/json",
            "Content-Type" = "application/json"
          ),
          body = body,
          encode = "raw"
        )
      }
      
      if(status_code(response) == 200) {
        return(fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE))
      } else {
        cat("  WARNING: HTTP", status_code(response), "for", endpoint, "\n")
        if(attempt < max_retries) {
          Sys.sleep(2 * attempt)
        }
      }
    }, error = function(e) {
      cat("  ERROR on attempt", attempt, ":", e$message, "\n")
      if(attempt < max_retries) {
        Sys.sleep(5 * attempt)
      }
    })
  }
  
  return(NULL)
}


#' Discover available RNA-seq molecular profiles for a study
#'
#' @param study_id cBioPortal study ID
#' @return Data frame of molecular profiles, or NULL
get_rnaseq_profiles <- function(study_id) {
  profiles <- cbio_query(paste0("/studies/", study_id, "/molecular-profiles"))
  
  if(is.null(profiles) || nrow(profiles) == 0) return(NULL)
  
  # Filter to RNA-seq expression profiles (not z-scores, not microarray)
  rna_profiles <- profiles %>%
    filter(grepl("rna_seq", molecularProfileId, ignore.case = TRUE)) %>%
    filter(!grepl("[Zz]score|[Zz]_score|median_Zscores|_all_sample_Zscores", molecularProfileId))
  
  return(rna_profiles)
}


#' Get all sample IDs for a study
#'
#' @param study_id cBioPortal study ID
#' @return Character vector of sample IDs
get_study_samples <- function(study_id) {
  samples <- cbio_query(paste0("/studies/", study_id, "/samples"))
  
  if(is.null(samples) || nrow(samples) == 0) return(NULL)
  
  return(samples$sampleId)
}


#' Fetch expression data for specific genes from a molecular profile
#'
#' @param profile_id Molecular profile ID
#' @param sample_ids Character vector of sample IDs
#' @param entrez_ids Numeric vector of Entrez gene IDs
#' @return Data frame with columns: sampleId, entrezGeneId, value
fetch_expression <- function(profile_id, sample_ids, entrez_ids) {
  # cBioPortal API has a limit on the number of samples per request
  # Split into batches if needed
  batch_size <- 500
  all_results <- list()
  
  for(i in seq(1, length(sample_ids), batch_size)) {
    batch_end <- min(i + batch_size - 1, length(sample_ids))
    batch_samples <- sample_ids[i:batch_end]
    
    body <- toJSON(list(
      entrezGeneIds = as.integer(entrez_ids),
      sampleIds = batch_samples
    ), auto_unbox = FALSE)
    
    result <- cbio_query(
      paste0("/molecular-profiles/", profile_id, "/molecular-data/fetch"),
      method = "POST",
      body = body
    )
    
    if(!is.null(result) && nrow(result) > 0) {
      all_results[[length(all_results) + 1]] <- result
    }
    
    Sys.sleep(SLEEP_BETWEEN_REQUESTS)
  }
  
  if(length(all_results) == 0) return(NULL)
  
  return(bind_rows(all_results))
}

# ============================================================
# MAIN PIPELINE
# ============================================================

cat("============================================================\n")
cat("cBioPortal TCGA Expression Data Download\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("API endpoint:", CBIO_API, "\n")
cat("Target genes:", paste(A3_GENES, collapse = ", "), "\n")
cat("Number of TCGA studies:", length(TCGA_STUDIES), "(32 PCA studies; COAD+READ combined)\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# ============================================================
# PHASE 0: Verify API connectivity
# ============================================================
cat("============================================================\n")
cat("PHASE 0: Verifying API connectivity\n")
cat("============================================================\n")

server_info <- cbio_query("/info")
if(is.null(server_info)) {
  stop("Cannot connect to cBioPortal API. Check internet connection.")
}
cat("  Connected to cBioPortal\n")
cat("  Portal version:", server_info$portalVersion, "\n")
cat("  DB build:", server_info$dbVersion, "\n\n")

# ============================================================
# PHASE 1: Discover molecular profiles for each study
# ============================================================
cat("============================================================\n")
cat("PHASE 1: Discovering RNA-seq molecular profiles\n")
cat("============================================================\n\n")

study_profiles <- list()
study_sample_counts <- c()

for(study_id in TCGA_STUDIES) {
  cancer_type <- toupper(gsub("_tcga_pan_can_atlas_2018", "", study_id))
  cat("  ", cancer_type, ": ", sep = "")
  
  # Get profiles
  rna_profiles <- get_rnaseq_profiles(study_id)
  
  if(is.null(rna_profiles) || nrow(rna_profiles) == 0) {
    cat("NO RNA-seq profiles found\n")
    next
  }
  
  # Prefer the non-z-score RNA-seq V2 profile
  # Try exact match first, then fallback
  profile_id <- NULL
  
  # Priority order for profile selection:
  preferred_patterns <- c(
    paste0(study_id, "_rna_seq_v2_mrna$"),       # Pan-Cancer Atlas specific
    "_rna_seq_v2_mrna$",                          # Generic RNA-seq V2
    "_rna_seq_mrna$"                              # Fallback: RNA-seq V1
  )
  
  for(pattern in preferred_patterns) {
    match_idx <- grep(pattern, rna_profiles$molecularProfileId)
    if(length(match_idx) > 0) {
      profile_id <- rna_profiles$molecularProfileId[match_idx[1]]
      break
    }
  }
  
  if(is.null(profile_id)) {
    # Just take the first non-z-score profile
    profile_id <- rna_profiles$molecularProfileId[1]
  }
  
  # Get sample count
  samples <- get_study_samples(study_id)
  n_samples <- length(samples)
  
  study_profiles[[study_id]] <- list(
    cancer_type = cancer_type,
    profile_id = profile_id,
    samples = samples,
    n_samples = n_samples,
    all_profiles = rna_profiles$molecularProfileId
  )
  
  study_sample_counts[cancer_type] <- n_samples
  
  cat(profile_id, " (", n_samples, " samples)\n", sep = "")
  
  Sys.sleep(SLEEP_BETWEEN_REQUESTS)
}

cat("\nTotal studies with RNA-seq profiles:", length(study_profiles), "\n")
cat("Total samples across all studies:", sum(study_sample_counts), "\n\n")

# ============================================================
# PHASE 2: Download A3 gene expression for all studies
# ============================================================
cat("============================================================\n")
cat("PHASE 2: Downloading A3 gene expression\n")
cat("============================================================\n\n")

all_a3_data <- list()
download_log <- data.frame(
  Cancer_Type = character(),
  Study_ID = character(),
  Profile_ID = character(),
  Samples_Queried = integer(),
  Samples_Returned = integer(),
  Genes_Returned = integer(),
  Status = character(),
  stringsAsFactors = FALSE
)

for(study_id in names(study_profiles)) {
  info <- study_profiles[[study_id]]
  cat("  Fetching ", info$cancer_type, " (", info$n_samples, " samples)... ", sep = "")
  
  expr_data <- fetch_expression(
    profile_id = info$profile_id,
    sample_ids = info$samples,
    entrez_ids = A3_ENTREZ
  )
  
  if(!is.null(expr_data) && nrow(expr_data) > 0) {
    # Add cancer type annotation
    expr_data$cancer_type <- info$cancer_type
    expr_data$study_id <- study_id
    expr_data$profile_id <- info$profile_id
    
    # Map Entrez ID back to gene symbol
    entrez_to_symbol <- setNames(names(A3_ENTREZ), as.character(A3_ENTREZ))
    expr_data$gene_symbol <- entrez_to_symbol[as.character(expr_data$entrezGeneId)]
    
    all_a3_data[[study_id]] <- expr_data
    
    n_samples_returned <- length(unique(expr_data$sampleId))
    n_genes_returned <- length(unique(expr_data$entrezGeneId))
    cat(n_samples_returned, "samples, ", n_genes_returned, "genes\n")
    
    download_log <- rbind(download_log, data.frame(
      Cancer_Type = info$cancer_type,
      Study_ID = study_id,
      Profile_ID = info$profile_id,
      Samples_Queried = info$n_samples,
      Samples_Returned = n_samples_returned,
      Genes_Returned = n_genes_returned,
      Status = "SUCCESS"
    ))
  } else {
    cat("FAILED\n")
    download_log <- rbind(download_log, data.frame(
      Cancer_Type = info$cancer_type,
      Study_ID = study_id,
      Profile_ID = info$profile_id,
      Samples_Queried = info$n_samples,
      Samples_Returned = 0,
      Genes_Returned = 0,
      Status = "FAILED"
    ))
  }
  
  Sys.sleep(SLEEP_BETWEEN_REQUESTS)
}

# Combine all data
cat("\nCombining all results...\n")
combined_a3 <- bind_rows(all_a3_data)
cat("Total rows (long format):", nrow(combined_a3), "\n")
cat("Unique samples:", length(unique(combined_a3$sampleId)), "\n")
cat("Unique genes:", length(unique(combined_a3$gene_symbol)), "\n\n")

# ============================================================
# PHASE 3: Reshape to wide format (samples × genes)
# ============================================================
cat("============================================================\n")
cat("PHASE 3: Reshaping to wide format\n")
cat("============================================================\n\n")

# Pivot to wide format: one row per sample, one column per gene
# Keep the sampleId and cancer_type as metadata
a3_wide <- combined_a3 %>%
  select(sampleId, cancer_type, gene_symbol, value) %>%
  distinct(sampleId, gene_symbol, .keep_all = TRUE) %>%
  tidyr::pivot_wider(
    id_cols = c(sampleId, cancer_type),
    names_from = gene_symbol,
    values_from = value
  ) %>%
  arrange(cancer_type, sampleId)

cat("Wide format dimensions:", nrow(a3_wide), "samples x", ncol(a3_wide), "columns\n")
cat("Columns:", paste(colnames(a3_wide), collapse = ", "), "\n\n")

# Check for NAs
na_counts <- colSums(is.na(a3_wide[, A3_GENES[A3_GENES %in% colnames(a3_wide)]]))
if(any(na_counts > 0)) {
  cat("WARNING: Some genes have NA values:\n")
  print(na_counts[na_counts > 0])
  cat("\n")
}

# ============================================================
# PHASE 4: Create Entity_ID-compatible format
# ============================================================
cat("============================================================\n")
cat("PHASE 4: Formatting for comparison with mystery file\n")
cat("============================================================\n\n")

# cBioPortal sample IDs are typically the 15-character TCGA barcode
# (e.g., TCGA-A1-A0SB-01), which is the first 15 characters of
# the full 28-character Entity_ID. The mystery file uses the full
# 28-character barcode, so we'll need to match on a truncated key.
#
# For comparison purposes, we'll store the cBioPortal sample ID
# and also create a 12-char Case_ID and 15-char Sample_ID for matching.

a3_output <- a3_wide %>%
  mutate(
    # cBioPortal typically uses 15-char sample barcodes
    Sample_Barcode = sampleId,
    # Extract 12-char Case_ID (TCGA-XX-XXXX)
    Case_ID = substr(sampleId, 1, 12),
    # Extract sample type code (positions 14-15, e.g., "01" = tumor)
    Sample_Type_Code = substr(sampleId, 14, 15),
    # Classify as Tumor or Normal
    Tissue_Type = ifelse(
      as.numeric(Sample_Type_Code) < 10, "Tumor",
      ifelse(as.numeric(Sample_Type_Code) < 20, "Normal", "Other")
    )
  ) %>%
  select(Sample_Barcode, Case_ID, cancer_type, Tissue_Type, 
         all_of(A3_GENES[A3_GENES %in% colnames(a3_wide)]))

cat("Output table dimensions:", nrow(a3_output), "x", ncol(a3_output), "\n")

# Count tumors vs normals
tissue_counts <- table(a3_output$Tissue_Type)
cat("Tissue distribution:\n")
print(tissue_counts)
cat("\n")

# Filter to tumors only (for comparison with mystery file which is tumors)
a3_tumors <- a3_output %>% filter(Tissue_Type == "Tumor")
cat("Tumor samples:", nrow(a3_tumors), "\n\n")

# Per-cancer-type summary
cancer_summary <- a3_tumors %>%
  group_by(cancer_type) %>%
  summarise(
    n_samples = n(),
    mean_A3A = mean(APOBEC3A, na.rm = TRUE),
    mean_A3B = mean(APOBEC3B, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

cat("Per-cancer-type summary (tumors only):\n")
print(as.data.frame(cancer_summary), row.names = FALSE)
cat("\n")

# ============================================================
# PHASE 5: Quick self-check — value ranges
# ============================================================
cat("============================================================\n")
cat("PHASE 5: Value range diagnostics\n")
cat("============================================================\n\n")

cat("This helps determine which RSEM normalization cBioPortal is serving.\n")
cat("If values are in the thousands, it is likely RSEM 'scaled_estimate × 10^6'\n")
cat("or raw RSEM expected counts. If values are 0-20ish, it is log2-transformed.\n\n")

for(gene in A3_GENES) {
  if(gene %in% colnames(a3_tumors)) {
    vals <- a3_tumors[[gene]]
    vals <- vals[!is.na(vals)]
    cat(sprintf("  %-12s  min=%10.2f  median=%10.2f  mean=%10.2f  max=%10.2f  zeros=%d/%d\n",
                gene, min(vals), median(vals), mean(vals), max(vals),
                sum(vals == 0), length(vals)))
  }
}

cat("\n")

# Determine if values appear to be log2-transformed
sample_vals <- a3_tumors$APOBEC3B[!is.na(a3_tumors$APOBEC3B)]
if(max(sample_vals, na.rm = TRUE) < 25) {
  cat(">>> VALUES APPEAR LOG2-TRANSFORMED (max < 25)\n")
  cat("    cBioPortal often serves log2(RSEM + 1) for RNA-seq V2.\n")
  cat("    To compare with raw RSEM, use 2^value - 1.\n")
  cat("    To compare with FPKM-UQ, rank-based comparison is recommended.\n")
} else if(max(sample_vals, na.rm = TRUE) > 1000) {
  cat(">>> VALUES APPEAR TO BE RAW RSEM COUNTS (max > 1000)\n")
  cat("    These are likely RSEM estimated counts or scaled estimates.\n")
} else {
  cat(">>> VALUE RANGE IS AMBIGUOUS (max =", max(sample_vals, na.rm = TRUE), ")\n")
  cat("    May need further investigation to determine normalization.\n")
}
cat("\n")

# ============================================================
# PHASE 6: Save outputs
# ============================================================
cat("============================================================\n")
cat("PHASE 6: Saving outputs\n")
cat("============================================================\n\n")

# Save A3 expression (tumors only) — primary output for comparison
a3_output_file <- file.path(OUTPUT_DIR, "cBioPortal_A3_Expression_All_TCGA.tsv")
write.table(a3_tumors, a3_output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", a3_output_file, "\n")
cat("  Dimensions:", nrow(a3_tumors), "x", ncol(a3_tumors), "\n\n")

# Save full output including normals
a3_all_file <- file.path(OUTPUT_DIR, "cBioPortal_A3_Expression_All_TCGA_with_Normals.tsv")
write.table(a3_output, a3_all_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", a3_all_file, "\n")
cat("  Dimensions:", nrow(a3_output), "x", ncol(a3_output), "\n\n")

# Save download log
log_file <- file.path(OUTPUT_DIR, "cBioPortal_download_log.tsv")
write.table(download_log, log_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", log_file, "\n\n")

# ============================================================
# PHASE 7: Generate comparison report
# ============================================================
cat("============================================================\n")
cat("PHASE 7: Download report\n")
cat("============================================================\n\n")

report_lines <- c(
  "============================================================",
  "cBioPortal TCGA Expression Download Report",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "============================================================",
  "",
  "PURPOSE:",
  "  Download TCGA RNA-seq expression data from cBioPortal (GRCh37/RSEM)",
  "  for comparison with the mystery file A3_Table_All_Samples_TCGA_Mohadeseh.tsv",
  "  to determine if cBioPortal was the original data source.",
  "",
  "NORMALIZATION NOTE:",
  "  cBioPortal serves RSEM-normalized expression from the legacy GRCh37 pipeline.",
  "  This is DIFFERENT from the FPKM-UQ values we compute from GRCh38/STAR-Counts.",
  "  See script header for detailed normalization comparison.",
  "",
  paste("API endpoint:", CBIO_API),
  paste("Portal version:", server_info$portalVersion),
  paste("Studies queried:", length(TCGA_STUDIES)),
  paste("Studies with data:", sum(download_log$Status == "SUCCESS")),
  paste("Studies failed:", sum(download_log$Status == "FAILED")),
  paste("Total tumor samples:", nrow(a3_tumors)),
  paste("Total samples (incl. normals):", nrow(a3_output)),
  "",
  "PER-STUDY DOWNLOAD LOG:",
  "========================"
)

for(i in 1:nrow(download_log)) {
  report_lines <- c(report_lines, sprintf(
    "  %-6s  %-45s  %5d queried  %5d returned  %s",
    download_log$Cancer_Type[i],
    download_log$Profile_ID[i],
    download_log$Samples_Queried[i],
    download_log$Samples_Returned[i],
    download_log$Status[i]
  ))
}

report_lines <- c(report_lines, "",
  "VALUE RANGE SUMMARY:",
  "===================="
)

for(gene in A3_GENES) {
  if(gene %in% colnames(a3_tumors)) {
    vals <- a3_tumors[[gene]][!is.na(a3_tumors[[gene]])]
    report_lines <- c(report_lines, sprintf(
      "  %-12s  min=%.2f  median=%.2f  mean=%.2f  max=%.2f  zeros=%d/%d",
      gene, min(vals), median(vals), mean(vals), max(vals),
      sum(vals == 0), length(vals)
    ))
  }
}

report_lines <- c(report_lines, "",
  "OUTPUT FILES:",
  "=============",
  paste("  1.", a3_output_file, "  (", nrow(a3_tumors), " tumors × ", ncol(a3_tumors), " cols)"),
  paste("  2.", a3_all_file, "  (", nrow(a3_output), " samples × ", ncol(a3_output), " cols)"),
  paste("  3.", log_file),
  "",
  "NEXT STEPS:",
  "===========",
  "  1. Run the updated Compare_A3_Expression_Sources.R script to compare",
  "     cBioPortal RSEM values against the mystery file and GDC FPKM-UQ values.",
  "  2. If cBioPortal values match the mystery file (r > 0.999), this confirms",
  "     the original file was sourced from cBioPortal/Firehose GRCh37 data.",
  "  3. This would explain why the mystery file does not match GRCh38 FPKM-UQ."
)

report_file <- file.path(OUTPUT_DIR, "cBioPortal_download_report.txt")
writeLines(report_lines, report_file)
cat("Report saved to:", report_file, "\n\n")

cat("============================================================\n")
cat("DOWNLOAD COMPLETE\n")
cat("============================================================\n")
cat("Total tumor samples downloaded:", nrow(a3_tumors), "\n")
cat("Output file:", a3_output_file, "\n")
cat("\nTo compare with the mystery file, update Compare_A3_Expression_Sources.R\n")
cat("to include the cBioPortal data as a fourth comparison source.\n")
cat("============================================================\n")
