#!/usr/bin/env Rscript
# Diagnostic script to identify the correct TCGAbiolinks column for Entity_ID
# Run this to see all available metadata columns and their values

library(TCGAbiolinks)
library(dplyr)

cat("============================================================\n")
cat("TCGAbiolinks Metadata Diagnostic Script\n")
cat("============================================================\n\n")

# Query TCGA-BRCA as test case
cat("Querying TCGA-BRCA for RNA-seq data...\n\n")

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Get the results/metadata
metadata <- getResults(query)

cat("============================================================\n")
cat("PART 1: ALL AVAILABLE COLUMNS\n")
cat("============================================================\n")
cat("Total columns:", ncol(metadata), "\n")
cat("Total rows (samples):", nrow(metadata), "\n\n")

cat("Column names:\n")
for(i in 1:length(colnames(metadata))) {
  cat(sprintf("  %3d. %s\n", i, colnames(metadata)[i]))
}

cat("\n============================================================\n")
cat("PART 2: COLUMNS CONTAINING 'sample', 'aliquot', 'submitter', or 'barcode'\n")
cat("============================================================\n")

interesting_cols <- grep("sample|aliquot|submitter|barcode|entity", 
                         colnames(metadata), 
                         ignore.case = TRUE, 
                         value = TRUE)

cat("Found", length(interesting_cols), "potentially relevant columns:\n\n")

for(col in interesting_cols) {
  cat("----------------------------------------\n")
  cat("Column:", col, "\n")
  cat("----------------------------------------\n")
  
  values <- as.character(metadata[[col]])
  values <- values[!is.na(values) & values != ""]
  
  if(length(values) > 0) {
    cat("  Sample values (first 5):\n")
    for(v in head(values, 5)) {
      cat("    '", v, "' (", nchar(v), " chars)\n", sep="")
    }
    
    # Check for 28-char values
    n_28_char <- sum(nchar(values) == 28)
    cat("  Values with exactly 28 chars:", n_28_char, "/", length(values), "\n")
    
    # Check if values match TCGA barcode pattern
    tcga_pattern <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[0-9]{2}$"
    n_match_pattern <- sum(grepl(tcga_pattern, values))
    cat("  Values matching full TCGA barcode pattern:", n_match_pattern, "/", length(values), "\n")
    
    if(n_match_pattern > 0) {
      cat("  >>> THIS COLUMN CONTAINS FULL TCGA ALIQUOT BARCODES! <<<\n")
    }
  } else {
    cat("  No non-empty values found\n")
  }
  cat("\n")
}

cat("\n============================================================\n")
cat("PART 3: EXAMINING ALL COLUMNS FOR 28-CHAR TCGA BARCODES\n")
cat("============================================================\n")

tcga_pattern <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[0-9]{2}$"

found_cols <- c()

for(col in colnames(metadata)) {
  values <- as.character(metadata[[col]])
  values <- values[!is.na(values) & values != ""]
  
  if(length(values) > 0) {
    n_match <- sum(grepl(tcga_pattern, values))
    if(n_match > 0) {
      cat("FOUND in column '", col, "': ", n_match, "/", length(values), " values match pattern\n", sep="")
      cat("  Example: '", values[grepl(tcga_pattern, values)][1], "'\n", sep="")
      found_cols <- c(found_cols, col)
    }
  }
}

if(length(found_cols) == 0) {
  cat("\nNo columns found with full 28-char TCGA aliquot barcodes.\n")
  cat("Will need to check if there's another way to get this data.\n")
}

cat("\n============================================================\n")
cat("PART 4: DETAILED VIEW OF FIRST ROW\n")
cat("============================================================\n")
cat("Showing all values for first sample:\n\n")

for(col in colnames(metadata)) {
  val <- as.character(metadata[[col]][1])
  if(!is.na(val) && val != "" && nchar(val) < 100) {
    cat(sprintf("%-60s = '%s'\n", col, val))
  } else if(!is.na(val) && val != "") {
    cat(sprintf("%-60s = '%s...'\n", col, substr(val, 1, 50)))
  }
}

cat("\n============================================================\n")
cat("PART 5: CHECK cases.samples NESTED STRUCTURE\n")
cat("============================================================\n")

# Sometimes the aliquot info is nested - let's check for any column with 'cases.samples'
cases_samples_cols <- grep("^cases\\.samples", colnames(metadata), value = TRUE)
cat("Found", length(cases_samples_cols), "columns starting with 'cases.samples':\n\n")

for(col in cases_samples_cols) {
  cat("Column:", col, "\n")
  val <- metadata[[col]][1]
  cat("  First value:", as.character(val), "\n")
  cat("  Length:", nchar(as.character(val)), "\n\n")
}

cat("\n============================================================\n")
cat("PART 6: TRYING ALTERNATIVE QUERY WITH LEGACY ARCHIVE\n")
cat("============================================================\n")

# Try to get sample sheet which might have more complete info
tryCatch({
  cat("Attempting to get sample sheet...\n")
  sample_sheet <- getSampleFilesSummary(query)
  cat("Sample sheet columns:\n")
  print(colnames(sample_sheet))
  cat("\nFirst few rows:\n")
  print(head(sample_sheet, 3))
}, error = function(e) {
  cat("Could not get sample sheet:", e$message, "\n")
})

cat("\n============================================================\n")
cat("PART 7: CHECK FOR BARCODE INFO IN FILE METADATA\n")
cat("============================================================\n")

# The barcode might be available through file annotation
tryCatch({
  # Get more detailed file info
  cat("Checking file_id and associated_entities columns...\n")
  
  if("associated_entities.entity_submitter_id" %in% colnames(metadata)) {
    cat("\nFound 'associated_entities.entity_submitter_id':\n")
    print(head(metadata$associated_entities.entity_submitter_id, 5))
  }
  
  if("cases.samples.portions.analytes.aliquots.submitter_id" %in% colnames(metadata)) {
    cat("\nFound 'cases.samples.portions.analytes.aliquots.submitter_id':\n")
    print(head(metadata$`cases.samples.portions.analytes.aliquots.submitter_id`, 5))
  }
  
}, error = function(e) {
  cat("Error checking columns:", e$message, "\n")
})

cat("\n============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("============================================================\n")
cat("\nPlease share the output above so we can identify the correct column\n")
cat("for the 28-character TCGA aliquot barcode (Entity_ID).\n")
