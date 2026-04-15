#!/usr/bin/env Rscript
# ============================================================
# Consolidate_MAFs_TCGA.R
#
# Step 4a: Consolidate all MuTect2 MAF files into:
#   1. A pan-cancer master MAF (maftools-ready)
#   2. A SNV-only table for SigProfiler input
#   3. Per-cancer-type MAF files for focused analysis
#
# The script first inspects one MAF to learn the structure,
# then efficiently reads and consolidates all 10,939 MAFs.
#
# Environment: RNA-seq_NovoGene conda environment
# Usage: Rscript Consolidate_MAFs_TCGA.R
#
# Input:  /master/jlehle/SHARED/TCGA/VCF/MuTect2_Annotated/TCGA-*/
#         /master/jlehle/SHARED/TCGA/VCF/manifests/TCGA_MuTect2_master_manifest.tsv
#
# Output: /master/jlehle/SHARED/TCGA/VCF/consolidated/
#           TCGA_pan_cancer_master.maf.tsv         — full master MAF (maftools-ready)
#           TCGA_pan_cancer_SNV_for_SigProfiler.tsv — SNV-only for SigProfiler
#           TCGA_pan_cancer_mutation_summary.tsv     — per-sample mutation counts
#           per_cancer/TCGA-{CANCER}_mutations.maf.tsv — per-cancer MAFs
# ============================================================

library(data.table)
library(dplyr)

# ============================================================
# CONFIGURATION
# ============================================================

base_dir       <- "/master/jlehle/SHARED/TCGA/VCF"
mutect2_dir    <- file.path(base_dir, "MuTect2_Annotated")
manifest_dir   <- file.path(base_dir, "manifests")
output_dir     <- file.path(base_dir, "consolidated")
per_cancer_dir <- file.path(output_dir, "per_cancer")

cancer_types <- c("BRCA", "LUAD", "LUSC", "PRAD", "COAD", "STAD",
                  "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML",
                  "PAAD", "ESCA", "PCPG", "READ", "TGCT", "THYM",
                  "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS",
                  "CHOL", "GBM", "HNSC", "KIRC", "LGG", "OV",
                  "SKCM", "THCA", "UCEC")

# Create output directories
for(d in c(output_dir, per_cancer_dir)) {
  if(!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("============================================================\n")
cat("TCGA MuTect2 MAF Consolidation Pipeline\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================
# PHASE 1: Inspect MAF structure from a sample file
# ============================================================
cat("============================================================\n")
cat("PHASE 1: Inspecting MAF file structure\n")
cat("============================================================\n\n")

# Find one MAF file to inspect
sample_maf <- list.files(mutect2_dir, pattern = "\\.maf\\.gz$",
                          recursive = TRUE, full.names = TRUE)[1]

if(is.na(sample_maf) || length(sample_maf) == 0) {
  stop("No MAF files found in ", mutect2_dir)
}

cat("Sample file:", basename(sample_maf), "\n")
cat("Size:", round(file.size(sample_maf) / 1024, 2), "KB\n\n")

# Read the file — GDC MAFs have comment lines starting with #
# First, count comment lines
con <- gzcon(file(sample_maf, "rb"))
comment_lines <- 0
header_line <- NULL
while(TRUE) {
  line <- readLines(con, n = 1)
  if(length(line) == 0) break
  if(startsWith(line, "#")) {
    comment_lines <- comment_lines + 1
    cat("Comment line", comment_lines, ":", substr(line, 1, 100), "\n")
  } else {
    header_line <- line
    break
  }
}
close(con)

cat("\nComment lines to skip:", comment_lines, "\n")

# Now read properly with data.table
sample_data <- fread(sample_maf, skip = comment_lines, header = TRUE,
                      sep = "\t", quote = "")

cat("Columns in MAF file:", ncol(sample_data), "\n")
cat("Rows (mutations):", nrow(sample_data), "\n\n")

cat("--- All column names ---\n")
for(i in seq_along(colnames(sample_data))) {
  cat(sprintf("  %3d. %s\n", i, colnames(sample_data)[i]))
}

# Key columns for maftools and SigProfiler
key_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
              "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
              "Variant_Classification", "Variant_Type",
              "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
              "FILTER", "t_depth", "t_ref_count", "t_alt_count",
              "n_depth", "n_ref_count", "n_alt_count",
              "Consequence", "IMPACT", "BIOTYPE")

cat("\n--- Key columns present ---\n")
for(col in key_cols) {
  present <- col %in% colnames(sample_data)
  cat(sprintf("  %-35s %s\n", col, ifelse(present, "YES", "*** MISSING ***")))
}

# Check Variant_Type values (SNP, INS, DEL, etc.)
cat("\n--- Variant_Type distribution ---\n")
if("Variant_Type" %in% colnames(sample_data)) {
  print(table(sample_data$Variant_Type, useNA = "ifany"))
}

# Check Variant_Classification values
cat("\n--- Variant_Classification distribution ---\n")
if("Variant_Classification" %in% colnames(sample_data)) {
  print(table(sample_data$Variant_Classification, useNA = "ifany"))
}

# Check FILTER values
cat("\n--- FILTER distribution ---\n")
if("FILTER" %in% colnames(sample_data)) {
  # FILTER can have complex values — show top 10
  filter_table <- sort(table(sample_data$FILTER), decreasing = TRUE)
  print(head(filter_table, 15))
}

# Check Tumor_Sample_Barcode format
cat("\n--- Tumor_Sample_Barcode examples ---\n")
if("Tumor_Sample_Barcode" %in% colnames(sample_data)) {
  barcodes <- unique(sample_data$Tumor_Sample_Barcode)
  cat("Unique barcodes:", length(barcodes), "\n")
  for(bc in head(barcodes, 5)) {
    cat("  '", bc, "' (", nchar(bc), " chars)\n", sep = "")
  }
}

# Show first few rows of key columns
cat("\n--- First 3 rows (key columns) ---\n")
display_cols <- intersect(c("Hugo_Symbol", "Chromosome", "Start_Position",
                             "Reference_Allele", "Tumor_Seq_Allele2",
                             "Variant_Classification", "Variant_Type",
                             "Tumor_Sample_Barcode", "FILTER"),
                           colnames(sample_data))
print(head(sample_data[, ..display_cols], 3))

# ============================================================
# PHASE 2: Load manifest for Entity_ID mapping
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 2: Loading manifest for metadata mapping\n")
cat("============================================================\n\n")

manifest_file <- file.path(manifest_dir, "TCGA_MuTect2_master_manifest.tsv")
manifest <- fread(manifest_file)

# Filter to MAF files only
maf_manifest <- manifest[Is_MAF == TRUE | Is_MAF == "TRUE"]
cat("MAF files in manifest:", nrow(maf_manifest), "\n")
cat("Cancer types:", length(unique(maf_manifest$Cancer_Type)), "\n")

# Build lookup: filename → metadata
maf_manifest$File_Base <- basename(maf_manifest$File_Name)
cat("Sample manifest entry:\n")
print(head(maf_manifest[, .(Cancer_Type, Entity_ID, Case_ID, File_Base)], 3))

# ============================================================
# PHASE 3: Read and consolidate all MAF files
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 3: Consolidating all MAF files\n")
cat("============================================================\n\n")

# Columns to keep — balance between completeness and memory
# These are the maftools-required + useful annotation columns
cols_to_keep <- intersect(
  c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
    "Chromosome", "Start_Position", "End_Position", "Strand",
    "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "dbSNP_RS", "dbSNP_Val_Status",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
    "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
    "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
    "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID",
    "Exon_Number", "FILTER",
    "t_depth", "t_ref_count", "t_alt_count",
    "n_depth", "n_ref_count", "n_alt_count",
    "Consequence", "IMPACT", "SYMBOL", "BIOTYPE",
    "SIFT", "PolyPhen", "EXON", "INTRON",
    "all_effects", "Existing_variation"),
  colnames(sample_data)
)

cat("Keeping", length(cols_to_keep), "columns per MAF\n\n")

# Process cancer types sequentially
all_mutations <- list()
total_mutations <- 0
total_files_read <- 0
failed_files <- character(0)

pipeline_start <- Sys.time()

for(cancer in cancer_types) {

  cancer_dir <- file.path(mutect2_dir, paste0("TCGA-", cancer))
  maf_files <- list.files(cancer_dir, pattern = "\\.maf\\.gz$", full.names = TRUE)

  if(length(maf_files) == 0) {
    cat(sprintf("TCGA-%-4s: No MAF files found, skipping\n", cancer))
    next
  }

  cat(sprintf("TCGA-%-4s: Reading %d MAF files... ", cancer, length(maf_files)))
  cancer_start <- Sys.time()

  cancer_mutations <- list()

  for(f in maf_files) {
    tryCatch({
      # Read MAF, skipping comment lines
      dt <- fread(f, skip = comment_lines, header = TRUE,
                  sep = "\t", quote = "", select = cols_to_keep)

      if(nrow(dt) > 0) {
        # Add cancer type metadata
        dt$Cancer_Type <- cancer
        dt$Project_ID <- paste0("TCGA-", cancer)
        dt$Source_File <- basename(f)

        cancer_mutations[[length(cancer_mutations) + 1]] <- dt
        total_files_read <- total_files_read + 1
      }

    }, error = function(e) {
      failed_files <<- c(failed_files, f)
    })
  }

  # Combine this cancer type
  if(length(cancer_mutations) > 0) {
    cancer_combined <- rbindlist(cancer_mutations, fill = TRUE)
    n_muts <- nrow(cancer_combined)
    total_mutations <- total_mutations + n_muts

    # Save per-cancer MAF
    cancer_out <- file.path(per_cancer_dir,
                            paste0("TCGA-", cancer, "_mutations.maf.tsv"))
    fwrite(cancer_combined, cancer_out, sep = "\t", quote = FALSE)

    all_mutations[[cancer]] <- cancer_combined

    elapsed <- as.numeric(difftime(Sys.time(), cancer_start, units = "secs"))
    n_samples <- length(unique(cancer_combined$Tumor_Sample_Barcode))
    cat(sprintf("%d mutations from %d samples (%.0fs)\n",
                n_muts, n_samples, elapsed))
  } else {
    cat("No mutations read\n")
  }
}

# ============================================================
# PHASE 4: Build pan-cancer master MAF
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 4: Building pan-cancer master MAF\n")
cat("============================================================\n\n")

master_maf <- rbindlist(all_mutations, fill = TRUE)

cat("Master MAF dimensions:", nrow(master_maf), "mutations x",
    ncol(master_maf), "columns\n")
cat("Unique samples:", length(unique(master_maf$Tumor_Sample_Barcode)), "\n")
cat("Cancer types:", length(unique(master_maf$Cancer_Type)), "\n\n")

# Summary statistics
cat("--- Variant_Type breakdown ---\n")
if("Variant_Type" %in% colnames(master_maf)) {
  vt <- sort(table(master_maf$Variant_Type), decreasing = TRUE)
  print(vt)
  cat("\n")
}

cat("--- Variant_Classification breakdown ---\n")
if("Variant_Classification" %in% colnames(master_maf)) {
  vc <- sort(table(master_maf$Variant_Classification), decreasing = TRUE)
  print(head(vc, 15))
  cat("\n")
}

cat("--- FILTER breakdown ---\n")
if("FILTER" %in% colnames(master_maf)) {
  ft <- sort(table(master_maf$FILTER), decreasing = TRUE)
  print(head(ft, 10))
  cat("\n")
}

cat("--- Mutations per cancer type ---\n")
cancer_summary <- master_maf %>%
  group_by(Cancer_Type) %>%
  summarize(
    n_mutations = n(),
    n_samples = length(unique(Tumor_Sample_Barcode)),
    n_SNP = sum(Variant_Type == "SNP", na.rm = TRUE),
    n_INS = sum(Variant_Type == "INS", na.rm = TRUE),
    n_DEL = sum(Variant_Type == "DEL", na.rm = TRUE),
    mutations_per_sample = round(n() / length(unique(Tumor_Sample_Barcode)), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_mutations))

print(as.data.frame(cancer_summary), row.names = FALSE)

# Save master MAF
master_maf_file <- file.path(output_dir, "TCGA_pan_cancer_master.maf.tsv")
cat("\nSaving master MAF:", master_maf_file, "\n")
fwrite(master_maf, master_maf_file, sep = "\t", quote = FALSE)
cat("File size:", round(file.size(master_maf_file) / 1e9, 2), "GB\n")

# Save cancer summary
summary_file <- file.path(output_dir, "TCGA_pan_cancer_mutation_summary.tsv")
write.table(cancer_summary, summary_file, sep = "\t",
            row.names = FALSE, quote = FALSE)

# ============================================================
# PHASE 5: Extract SNVs for SigProfiler
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 5: Extracting SNVs for SigProfiler\n")
cat("============================================================\n\n")

# SigProfiler needs SNPs only (single nucleotide variants)
# Filter criteria:
#   - Variant_Type == "SNP"
#   - Single base ref and alt
#   - Optionally filter by FILTER column (PASS only vs all)

# First, let's see what we have
cat("Total mutations in master MAF:", nrow(master_maf), "\n")

# Filter to SNPs
snv_data <- master_maf[Variant_Type == "SNP"]
cat("SNP mutations:", nrow(snv_data), "\n")

# Check single-base substitutions
snv_data <- snv_data[nchar(Reference_Allele) == 1 & nchar(Tumor_Seq_Allele2) == 1]
cat("Single-base SNPs:", nrow(snv_data), "\n")

# FILTER status
if("FILTER" %in% colnames(snv_data)) {
  cat("\nFILTER values for SNPs:\n")
  ft_snp <- sort(table(snv_data$FILTER), decreasing = TRUE)
  print(head(ft_snp, 10))

  # Note: GDC annotated MAFs include all calls, not just PASS
  # We'll create both PASS-only and all-calls versions
  n_pass <- sum(snv_data$FILTER == "PASS", na.rm = TRUE)
  cat("\nPASS SNPs:", n_pass, "\n")
  cat("Non-PASS SNPs:", nrow(snv_data) - n_pass, "\n")
}

# Create SigProfiler input format
# SigProfilerMatrixGenerator expects columns:
#   Sample  Genome  mut_type  chrom  pos_start  pos_end  ref  alt
#
# For the simple text input format:
#   Project  Sample  ID  Genome  mut_type  chrom  pos_start  pos_end  ref  alt  Type

# Map Tumor_Sample_Barcode to Entity_ID via manifest for consistency
# with the expression pipeline
barcode_to_entity <- maf_manifest %>%
  select(Entity_ID, Cancer_Type) %>%
  distinct()

# Build SigProfiler input — ALL SNPs (unfiltered)
sigprofiler_all <- data.frame(
  Project = snv_data$Project_ID,
  Sample = snv_data$Tumor_Sample_Barcode,
  ID = ".",
  Genome = "GRCh38",
  mut_type = "SNP",
  chrom = snv_data$Chromosome,
  pos_start = snv_data$Start_Position,
  pos_end = snv_data$End_Position,
  ref = snv_data$Reference_Allele,
  alt = snv_data$Tumor_Seq_Allele2,
  Type = "SOMATIC",
  stringsAsFactors = FALSE
)

sigprofiler_all_file <- file.path(output_dir,
                                   "TCGA_pan_cancer_SNV_for_SigProfiler_ALL.tsv")
fwrite(sigprofiler_all, sigprofiler_all_file, sep = "\t", quote = FALSE)
cat("\nSaved ALL SNVs for SigProfiler:", sigprofiler_all_file, "\n")
cat("  Rows:", nrow(sigprofiler_all), "\n")
cat("  Unique samples:", length(unique(sigprofiler_all$Sample)), "\n")

# Build SigProfiler input — PASS only
if("FILTER" %in% colnames(snv_data)) {
  snv_pass <- snv_data[FILTER == "PASS"]

  sigprofiler_pass <- data.frame(
    Project = snv_pass$Project_ID,
    Sample = snv_pass$Tumor_Sample_Barcode,
    ID = ".",
    Genome = "GRCh38",
    mut_type = "SNP",
    chrom = snv_pass$Chromosome,
    pos_start = snv_pass$Start_Position,
    pos_end = snv_pass$End_Position,
    ref = snv_pass$Reference_Allele,
    alt = snv_pass$Tumor_Seq_Allele2,
    Type = "SOMATIC",
    stringsAsFactors = FALSE
  )

  sigprofiler_pass_file <- file.path(output_dir,
                                      "TCGA_pan_cancer_SNV_for_SigProfiler_PASS.tsv")
  fwrite(sigprofiler_pass, sigprofiler_pass_file, sep = "\t", quote = FALSE)
  cat("\nSaved PASS SNVs for SigProfiler:", sigprofiler_pass_file, "\n")
  cat("  Rows:", nrow(sigprofiler_pass), "\n")
  cat("  Unique samples:", length(unique(sigprofiler_pass$Sample)), "\n")
}

# ============================================================
# PHASE 6: Per-sample mutation count table
# ============================================================
cat("\n\n============================================================\n")
cat("PHASE 6: Per-sample mutation count table\n")
cat("============================================================\n\n")

sample_counts <- master_maf %>%
  group_by(Cancer_Type, Tumor_Sample_Barcode) %>%
  summarize(
    total_mutations = n(),
    n_SNP = sum(Variant_Type == "SNP", na.rm = TRUE),
    n_INS = sum(Variant_Type == "INS", na.rm = TRUE),
    n_DEL = sum(Variant_Type == "DEL", na.rm = TRUE),
    n_DNP = sum(Variant_Type == "DNP", na.rm = TRUE),
    n_PASS = sum(FILTER == "PASS", na.rm = TRUE),
    n_SNP_PASS = sum(Variant_Type == "SNP" & FILTER == "PASS", na.rm = TRUE),
    .groups = "drop"
  )

counts_file <- file.path(output_dir, "TCGA_per_sample_mutation_counts.tsv")
fwrite(sample_counts, counts_file, sep = "\t", quote = FALSE)
cat("Saved per-sample counts:", counts_file, "\n")
cat("  Samples:", nrow(sample_counts), "\n")
cat("  Median mutations per sample:", median(sample_counts$total_mutations), "\n")
cat("  Median SNPs per sample:", median(sample_counts$n_SNP), "\n")
cat("  Median PASS SNPs per sample:", median(sample_counts$n_SNP_PASS), "\n")

# ============================================================
# SUMMARY
# ============================================================
total_elapsed <- as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))

cat("\n\n============================================================\n")
cat("CONSOLIDATION COMPLETE\n")
cat("============================================================\n")
cat("MAF files read:", total_files_read, "\n")
cat("Failed files:", length(failed_files), "\n")
cat("Total mutations:", total_mutations, "\n")
cat("Total time:", round(total_elapsed, 1), "minutes\n\n")

cat("Output files:\n")
cat("  1.", master_maf_file, "\n")
cat("     Pan-cancer master MAF (maftools-ready)\n")
cat("  2.", sigprofiler_all_file, "\n")
cat("     All SNVs for SigProfiler\n")
if(exists("sigprofiler_pass_file")) {
  cat("  3.", sigprofiler_pass_file, "\n")
  cat("     PASS-only SNVs for SigProfiler\n")
}
cat("  4.", counts_file, "\n")
cat("     Per-sample mutation counts\n")
cat("  5.", summary_file, "\n")
cat("     Per-cancer mutation summary\n")
cat("  6. per_cancer/ directory with 33 cancer-specific MAFs\n")

if(length(failed_files) > 0) {
  cat("\nFailed files:\n")
  for(f in failed_files) cat("  ", f, "\n")
}

cat("\nNext steps:\n")
cat("  1. Run SigProfiler (Python) on the SNV output\n")
cat("  2. Use maftools::read.maf() on the master MAF for visualization\n")
cat("  3. Quick maftools check:\n")
cat("     library(maftools)\n")
cat("     maf <- read.maf('", master_maf_file, "')\n", sep = "")
cat("     plotmafSummary(maf, rmOutlier = TRUE)\n")
cat("     oncoplot(maf, top = 20)\n")
cat("============================================================\n")
