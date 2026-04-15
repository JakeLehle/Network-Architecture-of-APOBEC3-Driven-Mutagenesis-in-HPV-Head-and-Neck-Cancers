#!/usr/bin/env Rscript
# ============================================================
# Diagnostic_cBioPortal_Study_IDs.R
#
# Quick diagnostic to discover ALL TCGA study IDs on cBioPortal
# and identify the correct IDs for studies that returned 404
# (COAD, READ, and any others).
#
# Usage: Rscript Diagnostic_cBioPortal_Study_IDs.R
# ============================================================

library(httr)
library(jsonlite)

CBIO_API <- "https://www.cbioportal.org/api"

cat("============================================================\n")
cat("cBioPortal TCGA Study ID Diagnostic\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================
# STEP 1: Get ALL studies on cBioPortal
# ============================================================
cat("STEP 1: Fetching all studies from cBioPortal...\n")

response <- GET(
  paste0(CBIO_API, "/studies"),
  add_headers("Accept" = "application/json")
)

all_studies <- fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)
cat("  Total studies on cBioPortal:", nrow(all_studies), "\n\n")

# ============================================================
# STEP 2: Filter to TCGA studies only
# ============================================================
cat("STEP 2: Filtering to TCGA studies...\n")

tcga_studies <- all_studies[grepl("tcga", all_studies$studyId, ignore.case = TRUE), ]
cat("  TCGA studies found:", nrow(tcga_studies), "\n\n")

# Show all TCGA study IDs grouped by cancer type
cat("ALL TCGA study IDs:\n")
cat("-------------------\n")

# Extract cancer type from study ID
tcga_studies$cancer_abbrev <- toupper(gsub("_tcga.*", "", tcga_studies$studyId))

# Sort by cancer type
tcga_studies <- tcga_studies[order(tcga_studies$cancer_abbrev), ]

for(i in 1:nrow(tcga_studies)) {
  cat(sprintf("  %-50s  %-6s  %s\n",
              tcga_studies$studyId[i],
              tcga_studies$cancer_abbrev[i],
              substr(tcga_studies$name[i], 1, 60)))
}

cat("\n")

# ============================================================
# STEP 3: Focus on Pan-Cancer Atlas 2018 studies
# ============================================================
cat("STEP 3: Pan-Cancer Atlas 2018 studies...\n\n")

pca_studies <- tcga_studies[grepl("pan_can_atlas_2018", tcga_studies$studyId), ]
cat("  Pan-Cancer Atlas 2018 studies:", nrow(pca_studies), "\n")

cat("\n  Pan-Cancer Atlas study IDs:\n")
for(i in 1:nrow(pca_studies)) {
  cat(sprintf("    %-50s  %s\n", pca_studies$studyId[i], pca_studies$cancer_abbrev[i]))
}

# ============================================================
# STEP 4: Find COAD and READ specifically
# ============================================================
cat("\n\nSTEP 4: Finding COAD and READ studies...\n\n")

# COAD
coad_studies <- tcga_studies[grepl("coad", tcga_studies$studyId, ignore.case = TRUE), ]
cat("  COAD studies:\n")
if(nrow(coad_studies) > 0) {
  for(i in 1:nrow(coad_studies)) {
    cat(sprintf("    %s  →  %s\n", coad_studies$studyId[i], coad_studies$name[i]))
  }
} else {
  cat("    NONE found with 'coad' in study ID\n")
}

# READ
read_studies <- tcga_studies[grepl("read", tcga_studies$studyId, ignore.case = TRUE), ]
cat("\n  READ studies:\n")
if(nrow(read_studies) > 0) {
  for(i in 1:nrow(read_studies)) {
    cat(sprintf("    %s  →  %s\n", read_studies$studyId[i], read_studies$name[i]))
  }
} else {
  cat("    NONE found with 'read' in study ID\n")
}

# Also check for colorectal combined study
crc_studies <- tcga_studies[grepl("coadread|colorectal|crc", tcga_studies$studyId, ignore.case = TRUE), ]
cat("\n  Colorectal combined studies:\n")
if(nrow(crc_studies) > 0) {
  for(i in 1:nrow(crc_studies)) {
    cat(sprintf("    %s  →  %s\n", crc_studies$studyId[i], crc_studies$name[i]))
  }
} else {
  cat("    NONE found\n")
}

# ============================================================
# STEP 5: Check which of the 33 cancer types have PCA studies
# ============================================================
cat("\n\nSTEP 5: Coverage check — all 33 TCGA cancer types\n\n")

expected_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
                    "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                    "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
                    "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT",
                    "THCA", "THYM", "UCEC", "UCS", "UVM")

cat(sprintf("  %-6s  %-50s  %-10s  %s\n", "TYPE", "PAN-CANCER ATLAS STUDY ID", "STATUS", "ALT STUDY ID"))
cat("  ", paste(rep("-", 120), collapse = ""), "\n", sep = "")

for(ct in expected_types) {
  # Check for pan-cancer atlas study
  pca_id <- paste0(tolower(ct), "_tcga_pan_can_atlas_2018")
  pca_exists <- pca_id %in% pca_studies$studyId
  
  # Find any alternative TCGA study for this cancer type
  alt_studies <- tcga_studies[grepl(paste0("^", tolower(ct), "_tcga"), tcga_studies$studyId), ]
  alt_ids <- alt_studies$studyId[alt_studies$studyId != pca_id]
  alt_str <- ifelse(length(alt_ids) > 0, paste(alt_ids, collapse = ", "), "")
  
  status <- ifelse(pca_exists, "OK", "MISSING")
  
  cat(sprintf("  %-6s  %-50s  %-10s  %s\n", ct, pca_id, status, alt_str))
}

# ============================================================
# STEP 6: For missing PCA studies, check RNA-seq availability
#         in alternative studies
# ============================================================
cat("\n\nSTEP 6: Checking RNA-seq profiles for alternative studies...\n\n")

missing_types <- c()
for(ct in expected_types) {
  pca_id <- paste0(tolower(ct), "_tcga_pan_can_atlas_2018")
  if(!(pca_id %in% pca_studies$studyId)) {
    missing_types <- c(missing_types, ct)
  }
}

if(length(missing_types) == 0) {
  cat("  All 33 cancer types have Pan-Cancer Atlas studies!\n")
} else {
  cat("  Missing from Pan-Cancer Atlas:", paste(missing_types, collapse = ", "), "\n\n")
  
  for(ct in missing_types) {
    cat("  --- ", ct, " ---\n", sep = "")
    
    # Find all TCGA studies for this cancer type
    alt_studies <- tcga_studies[grepl(paste0("^", tolower(ct), "_tcga|coadread"),
                                      tcga_studies$studyId, ignore.case = TRUE), ]
    
    if(nrow(alt_studies) == 0) {
      # Broader search
      alt_studies <- tcga_studies[grepl(tolower(ct), tcga_studies$studyId, ignore.case = TRUE), ]
    }
    
    if(nrow(alt_studies) == 0) {
      cat("    No TCGA studies found for this cancer type!\n\n")
      next
    }
    
    for(j in 1:nrow(alt_studies)) {
      sid <- alt_studies$studyId[j]
      cat("    Checking:", sid, "\n")
      
      # Get molecular profiles
      Sys.sleep(0.5)
      resp <- GET(
        paste0(CBIO_API, "/studies/", sid, "/molecular-profiles"),
        add_headers("Accept" = "application/json")
      )
      
      if(status_code(resp) == 200) {
        profiles <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
        rna_profiles <- profiles[grepl("rna_seq", profiles$molecularProfileId, ignore.case = TRUE), ]
        rna_profiles <- rna_profiles[!grepl("[Zz]score", rna_profiles$molecularProfileId), ]
        
        if(nrow(rna_profiles) > 0) {
          for(k in 1:nrow(rna_profiles)) {
            # Get sample count
            Sys.sleep(0.3)
            samp_resp <- GET(
              paste0(CBIO_API, "/studies/", sid, "/samples"),
              add_headers("Accept" = "application/json")
            )
            n_samp <- 0
            if(status_code(samp_resp) == 200) {
              samps <- fromJSON(content(samp_resp, "text", encoding = "UTF-8"), flatten = TRUE)
              n_samp <- nrow(samps)
            }
            
            cat(sprintf("      FOUND: %s (%d samples)\n",
                        rna_profiles$molecularProfileId[k], n_samp))
          }
        } else {
          cat("      No RNA-seq profiles\n")
        }
      } else {
        cat("      HTTP", status_code(resp), "\n")
      }
    }
    cat("\n")
  }
}

cat("============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("============================================================\n")
cat("\nUse the output above to update TCGA_STUDIES in\n")
cat("Download_cBioPortal_Expression_TCGA.R\n")
