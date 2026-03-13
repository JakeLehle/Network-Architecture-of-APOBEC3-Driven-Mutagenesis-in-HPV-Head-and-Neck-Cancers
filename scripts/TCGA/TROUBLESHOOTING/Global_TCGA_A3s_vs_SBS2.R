library(data.table)
library(dplyr)
library(ggplot2)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("Nuclear A3 Expression vs SBS2 Mutation Signature Analysis\n")
cat("============================================================\n\n")

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

# ============================================================
# STEP 2: Filter to matched samples
# ============================================================
cat("STEP 2: Filtering to matched samples\n")
cat("------------------------------------------------------------\n")

# Get Entity_IDs from mutation file
mutation_ids <- mutation_data$TCGA_Gene_Expression_Entity_ID

# Filter expression data to only include samples in mutation file
expression_matched <- expression_data %>%
  filter(Entity_ID %in% mutation_ids)

# Also filter mutation data to only include samples in expression file
mutation_matched <- mutation_data %>%
  filter(TCGA_Gene_Expression_Entity_ID %in% expression_matched$Entity_ID)

cat("Matched samples in expression file:", nrow(expression_matched), "\n")
cat("Matched samples in mutation file:", nrow(mutation_matched), "\n\n")

# ============================================================
# STEP 3: Sum nuclear A3 expression for each sample
# ============================================================
cat("STEP 3: Calculating Nuclear A3 sum for each sample\n")
cat("------------------------------------------------------------\n")

expression_matched <- expression_matched %>%
  mutate(Nuclear_A3s_summed = APOBEC3A + APOBEC3B + APOBEC3C + APOBEC3H)

cat("Added Nuclear_A3s_summed column\n")
cat("Summary of Nuclear_A3s_summed:\n")
print(summary(expression_matched$Nuclear_A3s_summed))
cat("\n")

# ============================================================
# STEP 4: Create cancer type mapping for mutation data
# ============================================================
cat("STEP 4: Creating cancer type mapping\n")
cat("------------------------------------------------------------\n")

# Create a lookup table from expression data (Entity_ID -> Project_ID)
cancer_type_lookup <- expression_matched %>%
  select(Entity_ID, Project_ID) %>%
  distinct()

cat("Cancer type lookup table created with", nrow(cancer_type_lookup), "entries\n")
cat("Cancer types found:", length(unique(cancer_type_lookup$Project_ID)), "\n")
cat("Cancer types:", paste(sort(unique(cancer_type_lookup$Project_ID)), collapse=", "), "\n\n")

# Add cancer type to mutation data
mutation_matched <- mutation_matched %>%
  left_join(cancer_type_lookup, by = c("TCGA_Gene_Expression_Entity_ID" = "Entity_ID"))

cat("Added Project_ID to mutation data\n")
cat("Samples with assigned cancer type:", sum(!is.na(mutation_matched$Project_ID)), "\n\n")

# ============================================================
# STEP 5: Calculate averages by cancer type
# ============================================================
cat("STEP 5: Calculating averages by cancer type\n")
cat("------------------------------------------------------------\n")

# Average Nuclear_A3s_summed by cancer type
expression_avg <- expression_matched %>%
  group_by(Project_ID) %>%
  summarise(
    Avg_Nuclear_A3s = mean(Nuclear_A3s_summed, na.rm = TRUE),
    N_samples_expr = n(),
    .groups = "drop"
  )

cat("Expression averages by cancer type:\n")
print(expression_avg)
cat("\n")

# Average SBS2 by cancer type
mutation_avg <- mutation_matched %>%
  group_by(Project_ID) %>%
  summarise(
    Avg_SBS2 = mean(SBS2, na.rm = TRUE),
    N_samples_mut = n(),
    .groups = "drop"
  )

cat("SBS2 averages by cancer type:\n")
print(mutation_avg)
cat("\n")

# ============================================================
# STEP 6: Create combined dataframe for plotting
# ============================================================
cat("STEP 6: Creating combined dataframe\n")
cat("------------------------------------------------------------\n")

# Merge expression and mutation averages
plot_data <- expression_avg %>%
  inner_join(mutation_avg, by = "Project_ID") %>%
  mutate(Cancer_Type = gsub("TCGA-", "", Project_ID))  # Clean up cancer type names for labels

cat("Combined dataframe:\n")
print(plot_data)
cat("\n")

# Save the summary data
summary_file <- "A3_SBS2_Cancer_Type_Summary.tsv"
write.table(plot_data, summary_file, sep="\t", row.names=FALSE, quote=FALSE)
cat("Saved summary data to:", summary_file, "\n\n")

# ============================================================
# STEP 7: Calculate correlation
# ============================================================
cat("STEP 7: Calculating correlation\n")
cat("------------------------------------------------------------\n")

correlation <- cor(plot_data$Avg_Nuclear_A3s, plot_data$Avg_SBS2, method = "pearson")
cor_test <- cor.test(plot_data$Avg_Nuclear_A3s, plot_data$Avg_SBS2, method = "pearson")

cat("Pearson correlation: r =", round(correlation, 4), "\n")
cat("P-value:", format(cor_test$p.value, scientific = TRUE, digits = 3), "\n\n")

# ============================================================
# STEP 8: Create the plot
# ============================================================
cat("STEP 8: Creating association plot\n")
cat("------------------------------------------------------------\n")

# Create a color palette with enough distinct colors for all cancer types
n_cancers <- nrow(plot_data)
cancer_colors <- scales::hue_pal()(n_cancers)
names(cancer_colors) <- plot_data$Cancer_Type

# Create the plot
p <- ggplot(plot_data, aes(x = Avg_Nuclear_A3s, y = Avg_SBS2, color = Cancer_Type)) +
  geom_point(size = 8) +
  scale_color_manual(values = cancer_colors) +
  labs(
    x = "Average Nuclear A3 Expression (FPKM-UQ)",
    y = "Average SBS2 Mutational Weight",
    color = "Cancer Type"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 26),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.5, "lines"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

# Save the plot
plot_file <- "Figure_Nuclear_A3_vs_SBS2_by_Cancer_Type.pdf"
ggsave(plot_file, p, width = 14, height = 10, units = "in")
cat("Saved plot to:", plot_file, "\n")

# Also save as PNG for quick viewing
plot_file_png <- "Figure_Nuclear_A3_vs_SBS2_by_Cancer_Type.png"
ggsave(plot_file_png, p, width = 14, height = 10, units = "in", dpi = 300)
cat("Saved plot to:", plot_file_png, "\n\n")

# ============================================================
# Summary
# ============================================================
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nKey findings:\n")
cat("  - Number of cancer types:", n_cancers, "\n")
cat("  - Pearson correlation (r):", round(correlation, 4), "\n")
cat("  - P-value:", format(cor_test$p.value, scientific = TRUE, digits = 3), "\n")

cat("\nOutput files:\n")
cat("  1.", summary_file, "- Summary data by cancer type\n")
cat("  2.", plot_file, "- Publication-ready PDF figure\n")
cat("  3.", plot_file_png, "- PNG version for quick viewing\n")

cat("\nPlot specifications:\n")
cat("  - X-axis: Average Nuclear A3 Expression (sum of A3A, A3B, A3C, A3H)\n")
cat("  - Y-axis: Average SBS2 Mutational Weight\n")
cat("  - Points: Colored by cancer type, size 8\n")
cat("  - Axis text size: 26\n")
cat("  - No title (publication-ready)\n")
cat("============================================================\n")
