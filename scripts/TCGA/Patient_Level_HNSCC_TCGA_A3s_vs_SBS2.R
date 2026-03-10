library(data.table)
library(dplyr)
library(ggplot2)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("Patient-Level Nuclear A3 vs SBS2 Analysis in TCGA-HNSC\n")
cat("============================================================\n\n")

# Note: TCGA uses "HNSC" for Head and Neck Squamous Cell Carcinoma
# If your data uses "HNSCC", change the filter below accordingly

# ============================================================
# STEP 1: Load the matched files
# ============================================================
cat("STEP 1: Loading matched files\n")
cat("------------------------------------------------------------\n")

expression_file <- "A3_Expression_FPKM_UQ_Matched.tsv"
mutation_file <- "Mutation_Signatures_Matched.tsv"

expression_data <- fread(expression_file)
mutation_data <- fread(mutation_file)

cat("Expression file:", nrow(expression_data), "rows x", ncol(expression_data), "columns\n")
cat("Mutation file:", nrow(mutation_data), "rows x", ncol(mutation_data), "columns\n\n")

# Check what cancer types are available
cat("Available cancer types in expression data:\n")
print(sort(unique(expression_data$Project_ID)))
cat("\n")

# Check Tissue_Type distribution
cat("Tissue_Type distribution in expression data:\n")
print(table(expression_data$Tissue_Type))
cat("\n")

# ============================================================
# STEP 2: Filter for HNSC TUMOR samples only (explicitly exclude normals)
# ============================================================
cat("STEP 2: Filtering for TCGA-HNSC TUMOR samples only\n")
cat("------------------------------------------------------------\n")

# Filter expression data for HNSC AND Tumor samples only
# Try both "TCGA-HNSC" and "TCGA-HNSCC" in case of naming variation
hnsc_expression <- expression_data %>%
  filter((Project_ID == "TCGA-HNSC" | Project_ID == "TCGA-HNSCC") & 
         Tissue_Type == "Tumor")

if(nrow(hnsc_expression) == 0) {
  cat("WARNING: No tumor samples found with Project_ID 'TCGA-HNSC' or 'TCGA-HNSCC'\n")
  cat("Checking for partial matches...\n")
  hnsc_expression <- expression_data %>%
    filter(grepl("HNSC|HNSCC", Project_ID, ignore.case = TRUE) & 
           Tissue_Type == "Tumor")
}

cat("HNSC TUMOR samples in expression file:", nrow(hnsc_expression), "\n")

# Verify no normal samples
cat("Tissue types in filtered data:\n")
print(table(hnsc_expression$Tissue_Type))
cat("\n")

# Get Entity_IDs from HNSC tumor expression samples
hnsc_entity_ids <- hnsc_expression$Entity_ID
cat("HNSC tumor Entity_IDs extracted:", length(hnsc_entity_ids), "\n\n")

# ============================================================
# STEP 3: Harmonize data between both files for HNSC tumors only
# ============================================================
cat("STEP 3: Harmonizing HNSC tumor samples between files\n")
cat("------------------------------------------------------------\n")

# Filter mutation data to only HNSC tumor samples
hnsc_mutation <- mutation_data %>%
  filter(TCGA_Gene_Expression_Entity_ID %in% hnsc_entity_ids)

cat("HNSC tumor samples in mutation file:", nrow(hnsc_mutation), "\n")

# Find common Entity_IDs (samples present in both files)
common_ids <- intersect(hnsc_expression$Entity_ID, hnsc_mutation$TCGA_Gene_Expression_Entity_ID)
cat("Common HNSC tumor samples in both files:", length(common_ids), "\n\n")

# Filter both to only common samples
hnsc_expression_matched <- hnsc_expression %>%
  filter(Entity_ID %in% common_ids)

hnsc_mutation_matched <- hnsc_mutation %>%
  filter(TCGA_Gene_Expression_Entity_ID %in% common_ids)

cat("Final matched HNSC tumor expression samples:", nrow(hnsc_expression_matched), "\n")
cat("Final matched HNSC tumor mutation samples:", nrow(hnsc_mutation_matched), "\n\n")

# ============================================================
# STEP 4: Calculate Nuclear A3 sum for each patient
# ============================================================
cat("STEP 4: Calculating Nuclear A3 sum for each patient\n")
cat("------------------------------------------------------------\n")

hnsc_expression_matched <- hnsc_expression_matched %>%
  mutate(Nuclear_A3s_summed = APOBEC3A + APOBEC3B + APOBEC3C + APOBEC3H)

cat("Added Nuclear_A3s_summed column\n")
cat("Summary of Nuclear_A3s_summed:\n")
print(summary(hnsc_expression_matched$Nuclear_A3s_summed))
cat("\n")

# ============================================================
# STEP 5: Create patient-level dataframe for plotting
# ============================================================
cat("STEP 5: Creating patient-level dataframe\n")
cat("------------------------------------------------------------\n")

# Select relevant columns from expression data
expression_subset <- hnsc_expression_matched %>%
  select(Entity_ID, Nuclear_A3s_summed)

# Select relevant columns from mutation data
mutation_subset <- hnsc_mutation_matched %>%
  select(TCGA_Gene_Expression_Entity_ID, SBS2) %>%
  rename(Entity_ID = TCGA_Gene_Expression_Entity_ID)

# Merge to create final patient-level dataframe
patient_data <- expression_subset %>%
  inner_join(mutation_subset, by = "Entity_ID")

cat("Patient-level dataframe created:\n")
cat("  Rows (tumor patients):", nrow(patient_data), "\n")
cat("  Columns:", paste(colnames(patient_data), collapse = ", "), "\n\n")

cat("First 10 rows of patient data:\n")
print(head(patient_data, 10))
cat("\n")

cat("Summary statistics:\n")
cat("  Nuclear_A3s_summed - Mean:", round(mean(patient_data$Nuclear_A3s_summed), 4), 
    ", SD:", round(sd(patient_data$Nuclear_A3s_summed), 4), "\n")
cat("  SBS2 - Mean:", round(mean(patient_data$SBS2), 4), 
    ", SD:", round(sd(patient_data$SBS2), 4), "\n\n")

# Save the patient-level data
patient_data_file <- "HNSC_Patient_Level_A3_SBS2_Data.tsv"
write.table(patient_data, patient_data_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved patient-level data to:", patient_data_file, "\n\n")

# ============================================================
# STEP 6: Calculate correlation
# ============================================================
cat("STEP 6: Calculating correlation\n")
cat("------------------------------------------------------------\n")

correlation <- cor(patient_data$Nuclear_A3s_summed, patient_data$SBS2, method = "pearson")
cor_test <- cor.test(patient_data$Nuclear_A3s_summed, patient_data$SBS2, method = "pearson")

cat("Pearson correlation: r =", round(correlation, 4), "\n")
cat("P-value:", format(cor_test$p.value, scientific = TRUE, digits = 3), "\n\n")

# Also calculate Spearman correlation (robust to outliers)
spearman_cor <- cor(patient_data$Nuclear_A3s_summed, patient_data$SBS2, method = "spearman")
spearman_test <- cor.test(patient_data$Nuclear_A3s_summed, patient_data$SBS2, method = "spearman")

cat("Spearman correlation: rho =", round(spearman_cor, 4), "\n")
cat("P-value:", format(spearman_test$p.value, scientific = TRUE, digits = 3), "\n\n")

# ============================================================
# STEP 7: Create the patient-level plot
# ============================================================
cat("STEP 7: Creating patient-level association plot\n")
cat("------------------------------------------------------------\n")

# Get total number of tumors for annotation
n_tumors <- nrow(patient_data)

# Create the plot with outlined points and sample count annotation
p <- ggplot(patient_data, aes(x = Nuclear_A3s_summed, y = SBS2)) +
  geom_point(size = 5, alpha = 0.4, fill = "#E41A1C", color = "black",
             shape = 21, stroke = 0.8) +
  # Add sample count annotation in top right
  annotate("text", 
           x = Inf, y = Inf, 
           label = paste0("(n = ", n_tumors, ")"),
           hjust = 1.1, vjust = 1.5,
           size = 8) +
  labs(
    x = "Nuclear A3 expression (A3A/A3B/A3C/A3H) FPKM-UQ",
    y = "SBS2 Mutational Weight"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 26),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# Save the plot
plot_file <- "Figure_HNSC_Patient_Level_A3_vs_SBS2.pdf"
ggsave(plot_file, p, width = 12, height = 10, units = "in")
cat("Saved plot to:", plot_file, "\n")

# Also save as PNG for quick viewing
plot_file_png <- "Figure_HNSC_Patient_Level_A3_vs_SBS2.png"
ggsave(plot_file_png, p, width = 12, height = 10, units = "in", dpi = 300)
cat("Saved plot to:", plot_file_png, "\n\n")

# ============================================================
# Summary
# ============================================================
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nKey findings for TCGA-HNSC TUMORS:\n")
cat("  - Number of tumor patients:", nrow(patient_data), "\n")
cat("  - Pearson correlation (r):", round(correlation, 4), "\n")
cat("  - Pearson p-value:", format(cor_test$p.value, scientific = TRUE, digits = 3), "\n")
cat("  - Spearman correlation (rho):", round(spearman_cor, 4), "\n")
cat("  - Spearman p-value:", format(spearman_test$p.value, scientific = TRUE, digits = 3), "\n")

cat("\nOutput files:\n")
cat("  1.", patient_data_file, "- Patient-level data\n")
cat("  2.", plot_file, "- Publication-ready PDF figure\n")
cat("  3.", plot_file_png, "- PNG version for quick viewing\n")

cat("\nPlot specifications:\n")
cat("  - X-axis: Nuclear A3 expression (A3A/A3B/A3C/A3H) FPKM-UQ\n")
cat("  - Y-axis: SBS2 Mutational Weight\n")
cat("  - Points: Individual TUMOR patient samples, size 5, alpha 0.4, black outline\n")
cat("  - Sample count: (n =", n_tumors, ") displayed in top right\n")
cat("  - Axis text size: 26\n")
cat("  - No title (publication-ready)\n")
cat("  - NOTE: Normal samples explicitly excluded\n")
cat("============================================================\n")
