#!/usr/bin/env Rscript
# ============================================================
# Download_MuTect2_VCFs_TCGA.R
#
# Purpose: Download MuTect2 Annotated Somatic Mutation VCF files
#          from GDC for all 33 TCGA cancer types.
#
# Strategy: Uses TCGAbiolinks for metadata queries only, then
#           downloads each file individually via the GDC REST API
#           (curl). This bypasses TCGAbiolinks' GDCdownload which
#           fails with "gzip: stdin: not in gzip format" on the
#           bulk tar.gz endpoint.
#
# Environment: RNA-seq_NovoGene conda environment
# Usage: Rscript Download_MuTect2_VCFs_TCGA.R
# ============================================================

library(TCGAbiolinks)
library(dplyr)
library(data.table)

# ============================================================
# CONFIGURATION
# ============================================================

gdc_token_file <- "/master/jlehle/SHARED/TCGA/KEYS/gdc-user-token.2026-04-08T19_49_58.080Z.txt"

base_output_dir <- "/master/jlehle/SHARED/TCGA/VCF"
mutect2_dir     <- file.path(base_output_dir, "MuTect2_Annotated")
manifest_dir    <- file.path(base_output_dir, "manifests")
log_dir         <- file.path(base_output_dir, "logs")

max_retries  <- 5
base_wait    <- 10

cancer_types <- c("BRCA", "LUAD", "LUSC", "PRAD", "COAD", "STAD",
                  "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML",
                  "PAAD", "ESCA", "PCPG", "READ", "TGCT", "THYM",
                  "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS",
                  "CHOL", "GBM", "HNSC", "KIRC", "LGG", "OV",
                  "SKCM", "THCA", "UCEC")

# ============================================================
# SETUP
# ============================================================

cat("============================================================\n")
cat("TCGA MuTect2 Annotated VCF Download Pipeline\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Download method: Direct GDC API (curl per file)\n")
cat("Output directory:", base_output_dir, "\n\n")

if(!file.exists(gdc_token_file)) {
  stop("GDC token file not found: ", gdc_token_file)
}
gdc_token <- rawToChar(readBin(gdc_token_file, "raw", file.info(gdc_token_file)$size))
gdc_token <- gsub("[^[:print:]]", "", gdc_token)
gdc_token <- gsub("\\s", "", gdc_token)
cat("GDC token loaded. Length:", nchar(gdc_token), "characters\n\n")

for(d in c(base_output_dir, mutect2_dir, manifest_dir, log_dir)) {
  if(!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ============================================================
# FUNCTION: Download a single file from GDC via curl
# ============================================================
download_gdc_file <- function(file_uuid, output_path, token, expected_size = NULL,
                               max_retries = 5, base_wait = 10) {

  for(attempt in 1:max_retries) {

    # Skip if already downloaded with correct size
    if(file.exists(output_path)) {
      actual_size <- file.size(output_path)
      if(!is.null(expected_size) && actual_size == expected_size) {
        return(list(success = TRUE, status = "exists", size = actual_size))
      } else if(is.null(expected_size) && actual_size > 0) {
        return(list(success = TRUE, status = "exists", size = actual_size))
      }
      file.remove(output_path)
    }

    curl_cmd <- paste0(
      "curl -s -w '\\n%{http_code}' ",
      "-H 'X-Auth-Token: ", token, "' ",
      "'https://api.gdc.cancer.gov/data/", file_uuid, "' ",
      "-o '", output_path, "'"
    )

    result <- tryCatch({
      output <- system(curl_cmd, intern = TRUE)
      http_code <- as.integer(tail(output, 1))

      if(http_code == 200 && file.exists(output_path)) {
        actual_size <- file.size(output_path)
        if(!is.null(expected_size) && actual_size != expected_size) {
          file.remove(output_path)
          list(success = FALSE, status = paste0("size_mismatch:", actual_size, "vs", expected_size))
        } else {
          list(success = TRUE, status = "downloaded", size = actual_size)
        }
      } else {
        if(file.exists(output_path)) file.remove(output_path)
        list(success = FALSE, status = paste0("http_", http_code))
      }
    }, error = function(e) {
      list(success = FALSE, status = paste0("error:", e$message))
    })

    if(result$success) return(result)

    if(attempt < max_retries) Sys.sleep(base_wait * attempt)
  }

  return(list(success = FALSE, status = "all_retries_exhausted"))
}

# ============================================================
# PHASE 1: Query metadata for all cancer types
# ============================================================
cat("============================================================\n")
cat("PHASE 1: Querying GDC metadata\n")
cat("============================================================\n\n")

all_metadata <- list()

for(cancer in cancer_types) {
  cat(sprintf("Querying TCGA-%s... ", cancer))

  for(q_attempt in 1:3) {
    tryCatch({
      query <- GDCquery(
        project = paste0("TCGA-", cancer),
        data.category = "Simple Nucleotide Variation",
        data.type = "Annotated Somatic Mutation",
        workflow.type = "MuTect2 Annotation"
      )

      metadata <- getResults(query)

      cases_raw <- as.character(metadata$cases)
      parsed_barcodes <- sapply(cases_raw, function(x) {
        trimws(strsplit(x, ",")[[1]][1])
      }, USE.NAMES = FALSE)

      sample_codes <- substr(parsed_barcodes, 14, 15)
      tissue_status <- ifelse(as.integer(sample_codes) < 10, "Tumor",
                       ifelse(as.integer(sample_codes) < 20, "Normal", "Other"))

      cancer_metadata <- data.frame(
        Cancer_Type = cancer,
        Project_ID = paste0("TCGA-", cancer),
        File_UUID = metadata$file_id,
        File_Name = metadata$file_name,
        File_Size = metadata$file_size,
        Entity_ID = parsed_barcodes,
        Case_ID = substr(parsed_barcodes, 1, 12),
        Tissue_Status = tissue_status,
        Sample_Type = metadata$sample_type,
        Is_VCF = grepl("\\.vcf\\.gz$", metadata$file_name),
        Is_MAF = grepl("\\.maf\\.gz$", metadata$file_name),
        stringsAsFactors = FALSE
      )

      all_metadata[[cancer]] <- cancer_metadata
      cat(nrow(metadata), "files (", sum(cancer_metadata$Is_VCF), "VCF,",
          sum(cancer_metadata$Is_MAF), "MAF )\n")
      break

    }, error = function(e) {
      if(q_attempt < 3) {
        cat("(retry)... ")
        Sys.sleep(30 * q_attempt)
      } else {
        cat("FAILED:", e$message, "\n")
      }
    })
  }
}

master_manifest <- bind_rows(all_metadata)
cat("\nTotal files:", nrow(master_manifest), "\n")
cat("  VCF:", sum(master_manifest$Is_VCF), "  MAF:", sum(master_manifest$Is_MAF), "\n")
cat("  Valid Entity_IDs:", sum(nchar(master_manifest$Entity_ID) == 28), "/",
    nrow(master_manifest), "\n")

manifest_file <- file.path(manifest_dir, "TCGA_MuTect2_master_manifest.tsv")
write.table(master_manifest, manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", manifest_file, "\n")

# ============================================================
# PHASE 2: Download files via direct GDC API
# ============================================================
cat("\n============================================================\n")
cat("PHASE 2: Downloading via GDC REST API\n")
cat("============================================================\n\n")

master_manifest$Download_Status <- "pending"
master_manifest$Download_Path <- NA_character_

total_files <- nrow(master_manifest)
downloaded <- 0; skipped <- 0; failed <- 0
pipeline_start <- Sys.time()

for(cancer in cancer_types) {

  cancer_idx <- which(master_manifest$Cancer_Type == cancer)
  n_cancer <- length(cancer_idx)
  if(n_cancer == 0) next

  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat(sprintf("TCGA-%s: %d files (%.1f MB)\n", cancer, n_cancer,
              sum(master_manifest$File_Size[cancer_idx]) / 1e6))
  cat(paste(rep("=", 60), collapse = ""), "\n")

  cancer_dir <- file.path(mutect2_dir, paste0("TCGA-", cancer))
  if(!dir.exists(cancer_dir)) dir.create(cancer_dir, recursive = TRUE)

  cancer_start <- Sys.time()
  cancer_ok <- 0; cancer_skip <- 0; cancer_fail <- 0

  for(j in seq_along(cancer_idx)) {
    i <- cancer_idx[j]

    output_path <- file.path(cancer_dir, master_manifest$File_Name[i])
    master_manifest$Download_Path[i] <- output_path

    result <- download_gdc_file(
      file_uuid = master_manifest$File_UUID[i],
      output_path = output_path,
      token = gdc_token,
      expected_size = master_manifest$File_Size[i],
      max_retries = max_retries,
      base_wait = base_wait
    )

    master_manifest$Download_Status[i] <- result$status

    if(result$success) {
      if(result$status == "exists") { cancer_skip <- cancer_skip + 1; skipped <- skipped + 1 }
      else { cancer_ok <- cancer_ok + 1; downloaded <- downloaded + 1 }
    } else {
      cancer_fail <- cancer_fail + 1; failed <- failed + 1
    }

    if(j %% 50 == 0 || j == n_cancer) {
      elapsed <- as.numeric(difftime(Sys.time(), cancer_start, units = "mins"))
      total_done <- downloaded + skipped + failed
      cat(sprintf("  [%d/%d] ok=%d skip=%d fail=%d (%.1f min) | Total: %d/%d\n",
                  j, n_cancer, cancer_ok, cancer_skip, cancer_fail, elapsed,
                  total_done, total_files))
    }
  }

  # Checkpoint after each cancer type
  write.table(master_manifest, manifest_file, sep = "\t",
              row.names = FALSE, quote = FALSE)
}

# ============================================================
# PHASE 3: Final outputs
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 3: Final outputs\n")
cat("============================================================\n\n")

write.table(master_manifest, manifest_file, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Saved:", manifest_file, "\n")

vcf_manifest <- master_manifest %>%
  filter(Is_VCF, Download_Status %in% c("downloaded", "exists"))
vcf_file <- file.path(manifest_dir, "TCGA_MuTect2_VCF_for_SigProfiler.tsv")
write.table(vcf_manifest, vcf_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", vcf_file, "(", nrow(vcf_manifest), "VCFs )\n")

maf_manifest <- master_manifest %>%
  filter(Is_MAF, Download_Status %in% c("downloaded", "exists"))
maf_file <- file.path(manifest_dir, "TCGA_MuTect2_MAF_only_manifest.tsv")
write.table(maf_manifest, maf_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", maf_file, "(", nrow(maf_manifest), "MAFs )\n")

for(cancer in cancer_types) {
  cancer_meta <- master_manifest %>% filter(Cancer_Type == cancer)
  meta_file <- file.path(manifest_dir, paste0("TCGA-", cancer, "_MuTect2_metadata.tsv"))
  write.table(cancer_meta, meta_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

# ============================================================
# SUMMARY
# ============================================================
total_elapsed <- as.numeric(difftime(Sys.time(), pipeline_start, units = "hours"))

cat("\n============================================================\n")
cat("DOWNLOAD SUMMARY\n")
cat("============================================================\n")
cat("Total files:", total_files, "\n")
cat("  Downloaded:", downloaded, "\n")
cat("  Already existed:", skipped, "\n")
cat("  Failed:", failed, "\n")
cat("Total time:", round(total_elapsed, 2), "hours\n\n")
cat("Status breakdown:\n")
print(table(master_manifest$Download_Status))

if(failed > 0) {
  cat("\nFailed by cancer type:\n")
  fail_summary <- master_manifest %>%
    filter(!Download_Status %in% c("downloaded", "exists")) %>%
    group_by(Cancer_Type) %>%
    summarize(n_failed = n(), .groups = "drop")
  print(as.data.frame(fail_summary), row.names = FALSE)
  cat("\nRe-run this script to retry. Completed files are skipped (size-verified).\n")
}

cat("\nDisk verification:\n")
for(cancer in cancer_types) {
  cdir <- file.path(mutect2_dir, paste0("TCGA-", cancer))
  n_disk <- if(dir.exists(cdir)) length(list.files(cdir, pattern = "\\.(vcf|maf)\\.gz$")) else 0
  n_exp <- sum(master_manifest$Cancer_Type == cancer)
  cat(sprintf("  TCGA-%-4s: %4d/%4d [%s]\n", cancer, n_disk, n_exp,
              ifelse(n_disk == n_exp, "OK", "INCOMPLETE")))
}

cat("\n============================================================\n")
cat("Pipeline complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
