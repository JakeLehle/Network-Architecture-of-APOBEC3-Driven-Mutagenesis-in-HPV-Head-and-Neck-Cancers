library(data.table)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(viridis)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("3D Analysis: Individual A3 Expression vs SBS2 in TCGA-HNSC\n")
cat("============================================================\n\n")

# ============================================================
# PLOT SETTINGS - Easy to modify
# ============================================================
# Choose colormap: "viridis", "magma", "plasma", "inferno", or "cividis"
COLORMAP <- "viridis"  # Change to "magma" if preferred

cat("Plot settings:\n")
cat("  - Colormap:", COLORMAP, "\n")
cat("  - Point size: 3\n")
cat("  - Point opacity: 0.5\n")
cat("  - Line width: 0.5\n")
cat("  - Line opacity: 0.5\n\n")

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

# Check Tissue_Type distribution
cat("Tissue_Type distribution in expression data:\n")
print(table(expression_data$Tissue_Type))
cat("\n")

# ============================================================
# STEP 2: Filter for HNSC TUMORS only (explicitly exclude normals)
# ============================================================
cat("STEP 2: Filtering for TCGA-HNSC TUMOR samples only\n")
cat("------------------------------------------------------------\n")

# Filter expression data for HNSC AND Tumor samples only
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
# STEP 4: Create patient-level dataframe for 3D plotting
# ============================================================
cat("STEP 4: Creating patient-level dataframe for 3D plot\n")
cat("------------------------------------------------------------\n")

# Select relevant columns from expression data (individual A3s)
expression_subset <- hnsc_expression_matched %>%
  select(Entity_ID, APOBEC3A, APOBEC3B, APOBEC3C, APOBEC3H)

# Select SBS2 from mutation data
mutation_subset <- hnsc_mutation_matched %>%
  select(TCGA_Gene_Expression_Entity_ID, SBS2) %>%
  rename(Entity_ID = TCGA_Gene_Expression_Entity_ID)

# Merge to create final patient-level dataframe
patient_data <- expression_subset %>%
  inner_join(mutation_subset, by = "Entity_ID")

cat("Patient-level dataframe created:\n")
cat("  Rows (tumor patients):", nrow(patient_data), "\n")
cat("  Columns:", paste(colnames(patient_data), collapse = ", "), "\n\n")

# Summary of individual A3 expression
cat("Summary of individual A3 expression:\n")
cat("  APOBEC3A - Mean:", round(mean(patient_data$APOBEC3A), 4), 
    ", Range:", round(min(patient_data$APOBEC3A), 4), "-", round(max(patient_data$APOBEC3A), 4), "\n")
cat("  APOBEC3B - Mean:", round(mean(patient_data$APOBEC3B), 4), 
    ", Range:", round(min(patient_data$APOBEC3B), 4), "-", round(max(patient_data$APOBEC3B), 4), "\n")
cat("  APOBEC3C - Mean:", round(mean(patient_data$APOBEC3C), 4), 
    ", Range:", round(min(patient_data$APOBEC3C), 4), "-", round(max(patient_data$APOBEC3C), 4), "\n")
cat("  APOBEC3H - Mean:", round(mean(patient_data$APOBEC3H), 4), 
    ", Range:", round(min(patient_data$APOBEC3H), 4), "-", round(max(patient_data$APOBEC3H), 4), "\n")
cat("  SBS2 - Mean:", round(mean(patient_data$SBS2), 4), 
    ", Range:", round(min(patient_data$SBS2), 4), "-", round(max(patient_data$SBS2), 4), "\n\n")

# Save the patient-level data
patient_data_file <- "HNSC_Patient_Level_Individual_A3_SBS2_Data.tsv"
write.table(patient_data, patient_data_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved patient-level data to:", patient_data_file, "\n\n")

# ============================================================
# STEP 5: Calculate correlations for each A3 vs SBS2
# ============================================================
cat("STEP 5: Calculating individual A3 correlations with SBS2\n")
cat("------------------------------------------------------------\n")

a3_genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3H")

correlation_results <- data.frame(
  Gene = character(),
  Pearson_r = numeric(),
  Pearson_pvalue = numeric(),
  Spearman_rho = numeric(),
  Spearman_pvalue = numeric(),
  stringsAsFactors = FALSE
)

for(gene in a3_genes) {
  pearson_test <- cor.test(patient_data[[gene]], patient_data$SBS2, method = "pearson")
  spearman_test <- cor.test(patient_data[[gene]], patient_data$SBS2, method = "spearman")
  
  correlation_results <- rbind(correlation_results, data.frame(
    Gene = gene,
    Pearson_r = round(pearson_test$estimate, 4),
    Pearson_pvalue = pearson_test$p.value,
    Spearman_rho = round(spearman_test$estimate, 4),
    Spearman_pvalue = spearman_test$p.value,
    stringsAsFactors = FALSE
  ))
  
  cat(gene, "vs SBS2:\n")
  cat("  Pearson r =", round(pearson_test$estimate, 4), 
      ", p =", format(pearson_test$p.value, scientific = TRUE, digits = 3), "\n")
  cat("  Spearman rho =", round(spearman_test$estimate, 4), 
      ", p =", format(spearman_test$p.value, scientific = TRUE, digits = 3), "\n\n")
}

# Save correlation results
correlation_file <- "HNSC_Individual_A3_SBS2_Correlations.tsv"
write.table(correlation_results, correlation_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved correlation results to:", correlation_file, "\n\n")

# ============================================================
# STEP 6: Create the 3D plot
# ============================================================
cat("STEP 6: Creating 3D plot\n")
cat("------------------------------------------------------------\n")

# Get total number of tumors for annotation
n_tumors <- nrow(patient_data)

# Set colorscale based on COLORMAP setting
# Plotly accepts "Viridis", "Magma", "Plasma", "Inferno", "Cividis"
plotly_colorscale <- switch(COLORMAP,
  "viridis" = "Viridis",
  "magma" = "Magma",
  "plasma" = "Plasma",
  "inferno" = "Inferno",
  "cividis" = "Cividis",
  "Viridis"  # default
)

cat("Using plotly colorscale:", plotly_colorscale, "\n")

# Create 3D scatter plot with plotly
# Point size: 3, opacity: 0.5, line width: 0.5, line opacity: 0.5
fig <- plot_ly(
  data = patient_data,
  x = ~APOBEC3A,
  y = ~APOBEC3B,
  z = ~APOBEC3C,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 3,
    color = ~SBS2,
    colorscale = plotly_colorscale,
    colorbar = list(
      title = list(
        text = "SBS2 Weight",
        font = list(size = 20)
      ),
      tickfont = list(size = 16)
    ),
    opacity = 0.7,
    line = list(
      color = "rgba(0, 0, 0, 0.5)",
      width = 0.5
    )
  ),
  text = ~paste("Entity_ID:", Entity_ID,
                "<br>A3A:", round(APOBEC3A, 2),
                "<br>A3B:", round(APOBEC3B, 2),
                "<br>A3C:", round(APOBEC3C, 2),
                "<br>SBS2:", round(SBS2, 4)),
  hoverinfo = "text"
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = list(
          text = "APOBEC3A (FPKM-UQ)",
          font = list(size = 20)
        ),
        tickfont = list(size = 16)
      ),
      yaxis = list(
        title = list(
          text = "APOBEC3B (FPKM-UQ)",
          font = list(size = 20)
        ),
        tickfont = list(size = 16)
      ),
      zaxis = list(
        title = list(
          text = "APOBEC3C (FPKM-UQ)",
          font = list(size = 20)
        ),
        tickfont = list(size = 16)
      ),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      )
    ),
    annotations = list(
      list(
        x = 1,
        y = 1,
        text = paste0("(n = ", n_tumors, ")"),
        showarrow = FALSE,
        xref = "paper",
        yref = "paper",
        font = list(size = 20)
      )
    ),
    margin = list(l = 0, r = 0, b = 0, t = 50)
  )

# Save as interactive HTML
html_file <- "Figure_HNSC_3D_Individual_A3_vs_SBS2.html"
saveWidget(fig, html_file, selfcontained = TRUE)
cat("Saved interactive 3D plot to:", html_file, "\n")

# Also create a static image using orca if available, or provide instructions
cat("\nNote: To save as static PDF/PNG, you can:\n")
cat("  1. Open the HTML file in a browser\n")
cat("  2. Use the camera icon in the top right to download as PNG\n")
cat("  3. Or install orca for programmatic export\n\n")

# Try to save static image if orca is available
tryCatch({
  # Save as PNG
  png_file <- "Figure_HNSC_3D_Individual_A3_vs_SBS2.png"
  orca(fig, png_file, width = 1200, height = 1000)
  cat("Saved static PNG to:", png_file, "\n")
  
  # Save as PDF
  pdf_file <- "Figure_HNSC_3D_Individual_A3_vs_SBS2.pdf"
  orca(fig, pdf_file, width = 1200, height = 1000)
  cat("Saved static PDF to:", pdf_file, "\n")
}, error = function(e) {
  cat("Could not save static images (orca not available).\n")
  cat("Please use the interactive HTML file or export manually.\n")
})

# ============================================================
# STEP 7: Create alternative static 3D plot using scatterplot3d
# ============================================================
cat("\nSTEP 7: Creating static 3D plot using scatterplot3d\n")
cat("------------------------------------------------------------\n")

library(scatterplot3d)

# Create color gradient based on SBS2 values using viridis
sbs2_scaled <- (patient_data$SBS2 - min(patient_data$SBS2)) / 
               (max(patient_data$SBS2) - min(patient_data$SBS2))

# Use viridis package for colors
n_colors <- 100
if(COLORMAP == "viridis") {
  color_palette <- viridis(n_colors)
} else if(COLORMAP == "magma") {
  color_palette <- magma(n_colors)
} else if(COLORMAP == "plasma") {
  color_palette <- plasma(n_colors)
} else if(COLORMAP == "inferno") {
  color_palette <- inferno(n_colors)
} else if(COLORMAP == "cividis") {
  color_palette <- cividis(n_colors)
} else {
  color_palette <- viridis(n_colors)
}

point_colors <- color_palette[ceiling(sbs2_scaled * (n_colors - 1)) + 1]

# Save as PDF - using smaller point size and adjusted alpha
pdf_file_static <- "Figure_HNSC_3D_Individual_A3_vs_SBS2_static.pdf"
pdf(pdf_file_static, width = 12, height = 10)

# Create the 3D scatter plot with alpha = 0.5
s3d <- scatterplot3d(
  x = patient_data$APOBEC3A,
  y = patient_data$APOBEC3B,
  z = patient_data$APOBEC3C,
  color = adjustcolor(point_colors, alpha.f = 0.5),
  pch = 21,
  bg = adjustcolor(point_colors, alpha.f = 0.5),
  cex.symbols = 0.8,  # Smaller points to match size 3
  xlab = "APOBEC3A (FPKM-UQ)",
  ylab = "APOBEC3B (FPKM-UQ)",
  zlab = "APOBEC3C (FPKM-UQ)",
  cex.lab = 1.8,
  cex.axis = 1.5,
  angle = 45,
  box = TRUE,
  grid = TRUE
)

# Add sample count
legend("topright", 
       legend = paste0("(n = ", n_tumors, ")"),
       bty = "n",
       cex = 1.8)

# Add color legend using viridis colors
legend_colors <- color_palette[c(1, 25, 50, 75, 100)]
legend("right",
       legend = c(paste0("High SBS2 (", round(max(patient_data$SBS2), 3), ")"),
                  "", "", "",
                  paste0("Low SBS2 (", round(min(patient_data$SBS2), 3), ")")),
       fill = rev(legend_colors),
       border = "black",
       title = "SBS2 Weight",
       cex = 1.2,
       bty = "n")

dev.off()
cat("Saved static 3D PDF to:", pdf_file_static, "\n")

# Also save as PNG
png_file_static <- "Figure_HNSC_3D_Individual_A3_vs_SBS2_static.png"
png(png_file_static, width = 1200, height = 1000, res = 100)

s3d <- scatterplot3d(
  x = patient_data$APOBEC3A,
  y = patient_data$APOBEC3B,
  z = patient_data$APOBEC3C,
  color = adjustcolor(point_colors, alpha.f = 0.5),
  pch = 21,
  bg = adjustcolor(point_colors, alpha.f = 0.5),
  cex.symbols = 0.8,  # Smaller points to match size 3
  xlab = "APOBEC3A (FPKM-UQ)",
  ylab = "APOBEC3B (FPKM-UQ)",
  zlab = "APOBEC3C (FPKM-UQ)",
  cex.lab = 1.8,
  cex.axis = 1.5,
  angle = 45,
  box = TRUE,
  grid = TRUE
)

legend("topright", 
       legend = paste0("(n = ", n_tumors, ")"),
       bty = "n",
       cex = 1.8)

legend_colors <- color_palette[c(1, 25, 50, 75, 100)]
legend("right",
       legend = c(paste0("High SBS2 (", round(max(patient_data$SBS2), 3), ")"),
                  "", "", "",
                  paste0("Low SBS2 (", round(min(patient_data$SBS2), 3), ")")),
       fill = rev(legend_colors),
       border = "black",
       title = "SBS2 Weight",
       cex = 1.2,
       bty = "n")

dev.off()
cat("Saved static 3D PNG to:", png_file_static, "\n\n")

# ============================================================
# Summary
# ============================================================
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nKey findings for TCGA-HNSC tumors:\n")
cat("  - Number of tumor patients:", n_tumors, "\n")
cat("\nIndividual A3 correlations with SBS2:\n")
print(correlation_results)

cat("\nOutput files:\n")
cat("  1.", patient_data_file, "- Patient-level data\n")
cat("  2.", correlation_file, "- Correlation results\n")
cat("  3.", html_file, "- Interactive 3D plot (HTML)\n")
cat("  4.", pdf_file_static, "- Static 3D PDF figure\n")
cat("  5.", png_file_static, "- Static 3D PNG figure\n")

cat("\nPlot specifications:\n")
cat("  - X-axis: APOBEC3A expression (FPKM-UQ)\n")
cat("  - Y-axis: APOBEC3B expression (FPKM-UQ)\n")
cat("  - Z-axis: APOBEC3C expression (FPKM-UQ)\n")
cat("  - Colormap:", COLORMAP, "(low SBS2 -> high SBS2)\n")
cat("  - Point size: 3\n")
cat("  - Point opacity: 0.5\n")
cat("  - Line width: 0.5\n")
cat("  - Line opacity: 0.5\n")
cat("  - Sample count: (n =", n_tumors, ") displayed\n")
cat("  - APOBEC3H excluded due to low expression\n")
cat("============================================================\n")
