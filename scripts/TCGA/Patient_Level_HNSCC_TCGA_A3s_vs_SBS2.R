library(data.table)
library(dplyr)
library(ggplot2)
library(ggbreak)

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
# STEP 7: Compute data-driven region boundaries
# ============================================================
cat("STEP 7: Computing data-driven region boundaries\n")
cat("------------------------------------------------------------\n")

# The x-axis boundary between "has A3" and "no A3"
X_THRESHOLD <- 10

# Compute the slope of the diagonal from the origin.
# The diagonal must be steep enough that ALL points with x <= X_THRESHOLD
# fall BELOW the line (i.e., in the cream region, not the teal region).
# This guarantees the teal region contains zero points.
#
# For each point in the left zone (0 < x <= X_THRESHOLD), the slope
# needed to place it exactly on the line is y/x. We take the maximum
# and add a 5% buffer.

left_points <- patient_data %>%
  filter(Nuclear_A3s_summed > 0 & Nuclear_A3s_summed <= X_THRESHOLD)

if(nrow(left_points) > 0) {
  slopes <- left_points$SBS2 / left_points$Nuclear_A3s_summed
  max_slope <- max(slopes)
  safe_slope <- max_slope * 1.05  # 5% buffer above steepest point

  # Identify the point that defined the boundary
  steepest_idx <- which.max(slopes)
  cat("Steepest point in left zone:\n")
  cat("  Nuclear_A3s_summed =", left_points$Nuclear_A3s_summed[steepest_idx], "\n")
  cat("  SBS2 =", left_points$SBS2[steepest_idx], "\n")
  cat("  Slope (y/x) =", round(max_slope, 4), "\n")
} else {
  # No points with x > 0 and x <= threshold; use default
  safe_slope <- 1.0
  cat("No points in left zone with A3 > 0. Using default slope of 1.0\n")
}

cat("Max slope in left zone:", round(max_slope, 4), "\n")
cat("Safe slope (with 5% buffer):", round(safe_slope, 4), "\n")

# The diagonal intersects x = X_THRESHOLD at this y value
y_at_threshold <- safe_slope * X_THRESHOLD
cat("Diagonal intersects x =", X_THRESHOLD, "at y =", round(y_at_threshold, 4), "\n")

# The horizontal boundary extending rightward from (X_THRESHOLD, y_at_threshold)
# separates coral (above) from cream (below) for x > X_THRESHOLD
cat("Horizontal boundary for x >", X_THRESHOLD, ": y =", round(y_at_threshold, 4), "\n\n")

# ============================================================
# STEP 8: Assign region labels and count per region
# ============================================================
cat("STEP 8: Assigning region labels to each tumor\n")
cat("------------------------------------------------------------\n")

# Region assignment logic:
#   - Above the diagonal (y > safe_slope * x) AND x <= X_THRESHOLD  →  "no_A3" (teal)
#   - Above the horizontal (y > y_at_threshold) AND x > X_THRESHOLD →  "active" (coral)
#   - Everything else                                                →  "a3_low_sbs2" (cream)
#
# By construction, the "no_A3" region should have zero points.

patient_data <- patient_data %>%
  mutate(
    region = case_when(
      Nuclear_A3s_summed <= X_THRESHOLD & SBS2 > safe_slope * Nuclear_A3s_summed ~ "no_A3",
      Nuclear_A3s_summed >  X_THRESHOLD & SBS2 > y_at_threshold                  ~ "active",
      TRUE                                                                        ~ "a3_low_sbs2"
    )
  )

# Count per region
region_counts <- table(patient_data$region)
n_total    <- nrow(patient_data)
n_no_a3    <- sum(patient_data$region == "no_A3")
n_active   <- sum(patient_data$region == "active")
n_low_sbs2 <- sum(patient_data$region == "a3_low_sbs2")

cat("Region counts:\n")
cat("  Teal   (no A3 / no SBS2):", n_no_a3, "\n")
cat("  Coral  (A3 + active SBS2):", n_active, "\n")
cat("  Cream  (A3, low SBS2):", n_low_sbs2, "\n")
cat("  Total:", n_total, "\n\n")

# Sanity check: teal region must be empty
if(n_no_a3 > 0) {
  cat("WARNING: ", n_no_a3, " points fell in the teal (no_A3) zone!\n")
  cat("  This should not happen. Review the slope calculation.\n\n")
} else {
  cat("CONFIRMED: Teal region contains 0 points (A3 expression is necessary for SBS2)\n\n")
}

# ============================================================
# STEP 9: Build polygon data for background regions
# ============================================================
cat("STEP 9: Building background region polygons\n")
cat("------------------------------------------------------------\n")

# Determine plot extent (with padding for axis breaks)
x_max <- max(patient_data$Nuclear_A3s_summed) * 1.05
y_max <- max(patient_data$SBS2) * 1.05

cat("Plot extent: x_max =", round(x_max, 2), ", y_max =", round(y_max, 2), "\n")

# Teal region (#9bc1bc): "no A3" — above the diagonal, left of X_THRESHOLD
# Vertices: origin → top-left → top at x=threshold → diagonal intersection
poly_teal <- data.frame(
  x = c(0, 0,            X_THRESHOLD, X_THRESHOLD),
  y = c(0, y_max,        y_max,       y_at_threshold),
  region = "no_A3"
)

# Coral region (#ed6a5a): "A3 + active SBS2" — above horizontal, right of X_THRESHOLD
poly_coral <- data.frame(
  x = c(X_THRESHOLD, X_THRESHOLD, x_max, x_max),
  y = c(y_at_threshold, y_max,    y_max,  y_at_threshold),
  region = "active"
)

# Cream region (#f4f1bb): "A3, low SBS2" — below the diagonal + horizontal line
# This is one continuous L-shaped polygon
poly_cream <- data.frame(
  x = c(0, X_THRESHOLD,  X_THRESHOLD, x_max, x_max, 0),
  y = c(0, y_at_threshold, y_at_threshold, y_at_threshold, 0, 0),
  region = "a3_low_sbs2"
)

# Combine polygon data
poly_data <- rbind(poly_teal, poly_coral, poly_cream)

cat("Polygon vertices defined for 3 regions\n\n")

# ============================================================
# STEP 10: Create the figure
# ============================================================
cat("STEP 10: Creating figure with colored regions and broken x-axis\n")
cat("------------------------------------------------------------\n")

# Define color mapping
region_colors <- c(
  "no_A3"      = "#9bc1bc",
  "active"     = "#ed6a5a",
  "a3_low_sbs2" = "#f4f1bb"
)

# Build region count labels
label_teal  <- paste0("n = ", n_no_a3)
label_coral <- paste0("n = ", n_active)
label_cream <- paste0("n = ", n_low_sbs2)

# Position labels within each region
# Teal: upper-left area
label_teal_x <- X_THRESHOLD * 0.35
label_teal_y <- y_max * 0.85

# Coral: upper-right area
label_coral_x <- (X_THRESHOLD + x_max) / 2
label_coral_y <- y_max * 0.85

# Cream: lower-right area
label_cream_x <- (X_THRESHOLD + x_max) / 2
label_cream_y <- y_at_threshold * 0.4

p <- ggplot() +
  # Layer 1: Background polygons
  geom_polygon(data = poly_teal,
               aes(x = x, y = y),
               fill = "#9bc1bc", alpha = 0.75) +
  geom_polygon(data = poly_coral,
               aes(x = x, y = y),
               fill = "#ed6a5a", alpha = 0.75) +
  geom_polygon(data = poly_cream,
               aes(x = x, y = y),
               fill = "#f4f1bb", alpha = 0.75) +

  # Layer 2: Diagonal boundary line (origin to threshold intersection)
  geom_segment(aes(x = 0, y = 0,
                   xend = X_THRESHOLD, yend = y_at_threshold),
               linetype = "dashed", color = "grey30", linewidth = 0.6) +

  # Layer 3: Horizontal boundary line (threshold intersection to x_max)
  geom_segment(aes(x = X_THRESHOLD, y = y_at_threshold,
                   xend = x_max, yend = y_at_threshold),
               linetype = "dashed", color = "grey30", linewidth = 0.6) +

  # Layer 4: Vertical boundary line at X_THRESHOLD
  geom_segment(aes(x = X_THRESHOLD, y = 0,
                   xend = X_THRESHOLD, yend = y_max),
               linetype = "dashed", color = "grey30", linewidth = 0.6) +

  # Layer 5: Data points colored by region
  geom_point(data = patient_data,
             aes(x = Nuclear_A3s_summed, y = SBS2, fill = region),
             shape = 21, size = 5, alpha = 0.7, color = "black", stroke = 0.8) +
  scale_fill_manual(values = region_colors, guide = "none") +

  # Layer 6: Per-region count annotations
  annotate("text", x = label_teal_x, y = label_teal_y,
           label = label_teal, size = 7, fontface = "bold", color = "grey30") +
  annotate("text", x = label_coral_x, y = label_coral_y,
           label = label_coral, size = 7, fontface = "bold", color = "grey30") +
  annotate("text", x = label_cream_x, y = label_cream_y,
           label = label_cream, size = 7, fontface = "bold", color = "grey30") +

  # Layer 7: Total count annotation
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("(n = ", n_total, ")"),
           hjust = 1.1, vjust = 1.5,
           size = 8) +

  # Axis labels
  labs(
    x = "Nuclear A3 expression (A3A/A3B/A3C/A3H) FPKM-UQ",
    y = "SBS2 Mutational Weight"
  ) +

  # Theme
  theme_bw() +
  theme(
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 26),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +

  # Broken x-axis to accommodate outlier
  scale_x_break(c(400, 600), scales = 0.3)

# ============================================================
# STEP 11: Save the figure
# ============================================================
cat("STEP 11: Saving figure\n")
cat("------------------------------------------------------------\n")

plot_file <- "Figure_HNSC_Patient_Level_A3_vs_SBS2.pdf"
ggsave(plot_file, p, width = 14, height = 10, units = "in")
cat("Saved plot to:", plot_file, "\n")

plot_file_png <- "Figure_HNSC_Patient_Level_A3_vs_SBS2.png"
ggsave(plot_file_png, p, width = 14, height = 10, units = "in", dpi = 300)
cat("Saved plot to:", plot_file_png, "\n\n")

# ============================================================
# Summary
# ============================================================
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nKey findings for TCGA-HNSC TUMORS:\n")
cat("  - Total tumor patients:", n_total, "\n")
cat("  - Teal region (no A3, no SBS2):", n_no_a3, "tumors\n")
cat("  - Coral region (A3 + active SBS2):", n_active, "tumors\n")
cat("  - Cream region (A3, low SBS2):", n_low_sbs2, "tumors\n")
cat("  - Diagonal slope:", round(safe_slope, 4),
    "(boundary angle:", round(atan(safe_slope) * 180 / pi, 1), "degrees)\n")
cat("  - Pearson correlation (r):", round(correlation, 4), "\n")
cat("  - Pearson p-value:", format(cor_test$p.value, scientific = TRUE, digits = 3), "\n")
cat("  - Spearman correlation (rho):", round(spearman_cor, 4), "\n")
cat("  - Spearman p-value:", format(spearman_test$p.value, scientific = TRUE, digits = 3), "\n")

cat("\nRegion geometry:\n")
cat("  - X threshold:", X_THRESHOLD, "FPKM-UQ\n")
cat("  - Diagonal: y =", round(safe_slope, 4), "* x (from origin to x =", X_THRESHOLD, ")\n")
cat("  - Horizontal boundary: y =", round(y_at_threshold, 4), "(for x >", X_THRESHOLD, ")\n")
cat("  - X-axis break: 400–600 FPKM-UQ\n")

cat("\nOutput files:\n")
cat("  1.", patient_data_file, "- Patient-level data with region assignments\n")
cat("  2.", plot_file, "- Publication-ready PDF figure\n")
cat("  3.", plot_file_png, "- PNG version for quick viewing\n")
cat("============================================================\n")
