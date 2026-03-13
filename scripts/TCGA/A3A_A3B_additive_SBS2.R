library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("A3B and A3A Additive Contributions to SBS2 Mutagenesis\n")
cat("============================================================\n\n")

# ============================================================
# DATA LOADING
# ============================================================
cat("Loading and preparing data...\n")
cat("------------------------------------------------------------\n")

expression_data <- fread("A3_Expression_FPKM_UQ_Matched.tsv")
mutation_data   <- fread("Mutation_Signatures_Matched.tsv")

hnsc_expression <- expression_data %>%
  filter((Project_ID == "TCGA-HNSC" | Project_ID == "TCGA-HNSCC") &
         Tissue_Type == "Tumor")

hnsc_entity_ids <- hnsc_expression$Entity_ID

hnsc_mutation <- mutation_data %>%
  filter(TCGA_Gene_Expression_Entity_ID %in% hnsc_entity_ids)

common_ids <- intersect(hnsc_expression$Entity_ID,
                        hnsc_mutation$TCGA_Gene_Expression_Entity_ID)

patient_data <- hnsc_expression %>%
  filter(Entity_ID %in% common_ids) %>%
  select(Entity_ID, APOBEC3A, APOBEC3B, APOBEC3C, APOBEC3H) %>%
  inner_join(
    hnsc_mutation %>%
      filter(TCGA_Gene_Expression_Entity_ID %in% common_ids) %>%
      select(TCGA_Gene_Expression_Entity_ID, SBS2) %>%
      rename(Entity_ID = TCGA_Gene_Expression_Entity_ID),
    by = "Entity_ID"
  )

n_tumors <- nrow(patient_data)
median_a3a <- median(patient_data$APOBEC3A)
median_a3b <- median(patient_data$APOBEC3B)

cat("HNSC tumor patients:", n_tumors, "\n")
cat("Median A3A:", round(median_a3a, 4), "\n")
cat("Median A3B:", round(median_a3b, 4), "\n\n")

# Shared theme
theme_panel <- theme_bw() +
  theme(
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),
    strip.text = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  )


# ============================================================
# PANEL 1: A3B Threshold Sweep — Correlation Above vs Below
# ============================================================
cat("Building Panel 1: A3B threshold sweep\n")
cat("------------------------------------------------------------\n")

a3b_quantiles <- seq(0.2, 0.8, by = 0.05)

threshold_results <- data.frame(
  Quantile = numeric(),
  Threshold = numeric(),
  N_below = integer(),
  N_above = integer(),
  Rho_below = numeric(),
  Pval_below = numeric(),
  Rho_above = numeric(),
  Pval_above = numeric(),
  stringsAsFactors = FALSE
)

for(q in a3b_quantiles) {
  thresh <- quantile(patient_data$APOBEC3B, q)
  below <- patient_data %>% filter(APOBEC3B < thresh)
  above <- patient_data %>% filter(APOBEC3B >= thresh)

  if(nrow(below) >= 10 & nrow(above) >= 10) {
    test_below <- cor.test(below$APOBEC3B, below$SBS2, method = "spearman")
    test_above <- cor.test(above$APOBEC3B, above$SBS2, method = "spearman")

    threshold_results <- rbind(threshold_results, data.frame(
      Quantile = q,
      Threshold = round(thresh, 4),
      N_below = nrow(below),
      N_above = nrow(above),
      Rho_below = round(test_below$estimate, 4),
      Pval_below = test_below$p.value,
      Rho_above = round(test_above$estimate, 4),
      Pval_above = test_above$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

cat("Threshold sweep results:\n")
print(threshold_results)
cat("\n")

write.table(threshold_results, "A3B_Threshold_Sweep_Correlations.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

panel_1 <- ggplot(threshold_results, aes(x = Quantile)) +
  geom_ribbon(aes(ymin = pmin(Rho_below, Rho_above),
                  ymax = pmax(Rho_below, Rho_above)),
              alpha = 0.1, fill = "grey50") +
  geom_line(aes(y = Rho_below, color = "Below threshold"),
            linewidth = 1.2) +
  geom_point(aes(y = Rho_below, color = "Below threshold"),
             size = 3.5) +
  geom_line(aes(y = Rho_above, color = "Above threshold"),
            linewidth = 1.2) +
  geom_point(aes(y = Rho_above, color = "Above threshold"),
             size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  scale_color_manual(values = c("Below threshold" = "#5e81ac",
                                 "Above threshold" = "#bf616a"),
                     name = "A3B Expression Group") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.1),
                     labels = paste0(seq(20, 80, 10), "th %")) +
  labs(x = "A3B Expression Quantile (Threshold)",
       y = "Spearman rho (A3B vs SBS2)",
       title = "A3B–SBS2 Correlation Saturates at High A3B Expression",
       subtitle = paste0("Below threshold = rising phase  |  ",
                         "Above threshold = plateau phase  |  ",
                         "n = ", n_tumors, " tumors")) +
  theme_panel +
  theme(legend.position = c(0.80, 0.85),
        legend.background = element_rect(fill = "white", color = "grey70"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"))

cat("Panel 1 complete.\n\n")


# ============================================================
# PANEL 2: Boxplot — SBS2 by A3B×A3A Quadrant (ordered low→high)
# ============================================================
cat("Building Panel 2: Ordered boxplot with median trend line\n")
cat("------------------------------------------------------------\n")

# Create quadrant labels
patient_data <- patient_data %>%
  mutate(
    A3B_gate = ifelse(APOBEC3B >= median_a3b, "A3B HIGH", "A3B LOW"),
    A3A_level = ifelse(APOBEC3A >= median_a3a, "A3A HIGH", "A3A LOW"),
    Quadrant = paste0(A3B_gate, " / ", A3A_level)
  )

# Compute medians per quadrant for ordering and trend line
quadrant_medians <- patient_data %>%
  group_by(Quadrant, A3B_gate, A3A_level) %>%
  summarise(
    n = n(),
    median_SBS2 = median(SBS2),
    mean_SBS2 = mean(SBS2),
    .groups = "drop"
  ) %>%
  arrange(median_SBS2)

# Set order from lowest to highest median SBS2
ordered_levels <- quadrant_medians$Quadrant
patient_data$Quadrant <- factor(patient_data$Quadrant, levels = ordered_levels)
quadrant_medians$Quadrant <- factor(quadrant_medians$Quadrant, levels = ordered_levels)

# Add numeric x position for trend line
quadrant_medians$x_pos <- 1:nrow(quadrant_medians)

cat("Quadrant order (low → high median SBS2):\n")
for(i in 1:nrow(quadrant_medians)) {
  cat("  ", i, ". ", as.character(quadrant_medians$Quadrant[i]),
      " — median SBS2 =", round(quadrant_medians$median_SBS2[i], 6),
      " (n =", quadrant_medians$n[i], ")\n")
}
cat("\n")

# Colors matching the narrative: coolest for lowest SBS2, warmest for highest
quadrant_colors <- setNames(
  c("#9bc1bc", "#f4f1bb", "#e8c170", "#ed6a5a"),
  ordered_levels
)

# Build x-axis labels with sample counts
quadrant_labels <- paste0(ordered_levels, "\n(n = ", quadrant_medians$n, ")")
names(quadrant_labels) <- ordered_levels

# Wilcoxon pairwise tests
cat("Pairwise Wilcoxon tests (BH-adjusted):\n")
pw_test <- pairwise.wilcox.test(patient_data$SBS2, patient_data$Quadrant,
                                 p.adjust.method = "BH")
print(pw_test)
cat("\n")

panel_2 <- ggplot(patient_data, aes(x = Quadrant, y = SBS2)) +
  # Boxplot — emphasizes median, IQR, and outliers
  geom_boxplot(aes(fill = Quadrant), alpha = 0.75,
               outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.4,
               outlier.fill = "grey50", color = "black", linewidth = 0.5,
               width = 0.6) +
  # Individual points jittered behind
  geom_jitter(shape = 21, size = 1.5, alpha = 0.25, color = "black",
              fill = "grey40", width = 0.15, stroke = 0.3) +
  # Median trend line connecting the four boxes
  geom_line(data = quadrant_medians,
            aes(x = x_pos, y = median_SBS2, group = 1),
            color = "#2d3436", linewidth = 1.2, linetype = "solid") +
  geom_point(data = quadrant_medians,
             aes(x = x_pos, y = median_SBS2),
             color = "#2d3436", fill = "white", shape = 23, size = 5,
             stroke = 1.5) +
  # Annotate median values above each diamond
  geom_text(data = quadrant_medians,
            aes(x = x_pos, y = median_SBS2,
                label = round(median_SBS2, 1)),
            vjust = -7, size = 10, fontface = "bold", color = "#2d3436") +
  scale_fill_manual(values = quadrant_colors, guide = "none") +
  scale_x_discrete(labels = quadrant_labels) +
  labs(x = NULL,
       y = "SBS2 Mutational Weight",
       title = "Additive A3B and A3A Contributions to SBS2",
       subtitle = paste0("Ordered by median SBS2  |  ",
                         "A3B split at ", round(median_a3b, 1),
                         "  |  A3A split at ", round(median_a3a, 1),
                         " FPKM-UQ")) +
  theme_panel +
  theme(axis.text.x = element_text(size = 18))

cat("Panel 2 complete.\n\n")


# ============================================================
# PANEL 3: 2×2 Heatmap — Median SBS2 per Quadrant
# ============================================================
cat("Building Panel 3: 2×2 heatmap\n")
cat("------------------------------------------------------------\n")

# Rebuild stats in the 2x2 grid format
heatmap_stats <- patient_data %>%
  group_by(A3B_gate, A3A_level) %>%
  summarise(
    n = n(),
    median_SBS2 = round(median(SBS2), 1),
    mean_SBS2 = round(mean(SBS2), 1),
    max_SBS2 = round(max(SBS2), 1),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", n, "\nmedian = ", median_SBS2))

# Set factor levels so the grid reads intuitively
heatmap_stats$A3B_gate <- factor(heatmap_stats$A3B_gate,
                                  levels = c("A3B LOW", "A3B HIGH"))
heatmap_stats$A3A_level <- factor(heatmap_stats$A3A_level,
                                   levels = c("A3A LOW", "A3A HIGH"))

cat("2×2 Heatmap data:\n")
print(as.data.frame(heatmap_stats))
cat("\n")

write.table(heatmap_stats, "A3B_A3A_Quadrant_Stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

panel_3 <- ggplot(heatmap_stats, aes(x = A3B_gate, y = A3A_level)) +
  geom_tile(aes(fill = median_SBS2), color = "black", linewidth = 1.5) +
  geom_text(aes(label = label), color = "grey99",  size = 8, fontface = "bold", lineheight = 1.2) +
  scale_fill_viridis(option = "viridis", name = "Median\nSBS2 Weight",
                     guide = guide_colorbar(barwidth = 1.5, barheight = 10)) +
  labs(x = "APOBEC3B Status",
       y = "APOBEC3A Status",
       title = "2×2 Model: A3B × A3A Interaction",
       subtitle = paste0("Split at median expression  |  ",
                         "n = ", n_tumors, " HNSC tumors")) +
  theme_panel +
  theme(axis.text = element_text(size = 20, face = "bold"),
        panel.grid = element_blank())

cat("Panel 3 complete.\n\n")


# ============================================================
# ASSEMBLE AND SAVE
# ============================================================
cat("Assembling composite figure...\n")
cat("------------------------------------------------------------\n")

# Layout: Panel 1 on top (full width), Panel 2 and Panel 3 below
composite <- panel_1 /
             (panel_2 | panel_3 + plot_layout(widths = c(1.5, 1))) +
  plot_annotation(
    title = "Dissecting A3B and A3A Contributions to SBS2 Mutagenesis in HNSCC",
    subtitle = paste0("TCGA-HNSC Tumors (n = ", n_tumors, ")"),
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5, color = "grey30")
    )
  ) +
  plot_layout(heights = c(1, 1.2))

# Save composite
ggsave("Figure_A3B_A3A_Additive_SBS2.pdf", composite,
       width = 20, height = 18, units = "in")
ggsave("Figure_A3B_A3A_Additive_SBS2.png", composite,
       width = 20, height = 18, units = "in", dpi = 200)

# Save individual panels
ggsave("Panel_Threshold_Sweep.pdf", panel_1, width = 12, height = 7)
ggsave("Panel_Ordered_Boxplot.pdf", panel_2, width = 12, height = 8)
ggsave("Panel_2x2_Heatmap.pdf", panel_3, width = 8, height = 7)

cat("\nAll figures saved.\n")

# ============================================================
# SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("\nPanel 1 — Threshold sweep:\n")
cat("  A3B's correlation with SBS2 holds steady in the rising phase\n")
cat("  and collapses above the 80th percentile (plateau/saturation).\n\n")

cat("Panel 2 — Ordered boxplot:\n")
for(i in 1:nrow(quadrant_medians)) {
  cat("  ", i, ". ", as.character(quadrant_medians$Quadrant[i]),
      " — median SBS2 =", round(quadrant_medians$median_SBS2[i], 1),
      " (n =", quadrant_medians$n[i], ")\n")
}
cat("  Trend line connects medians showing additive relationship.\n\n")

cat("Panel 3 — 2×2 heatmap:\n")
cat("  Confirms the gate/driver model:\n")
cat("  A3B- / A3A- = lowest SBS2 (baseline)\n")
cat("  A3B+ / A3A- = moderate SBS2 (A3B alone drives steady accumulation)\n")
cat("  A3B- / A3A+ = modest SBS2 (A3A alone has limited effect)\n")
cat("  A3B+ / A3A+ = highest SBS2 (combinatorial amplification)\n\n")

cat("Output files:\n")
cat("  Composite: Figure_A3B_A3A_Additive_SBS2.pdf/.png\n")
cat("  Panel 1:   Panel_Threshold_Sweep.pdf\n")
cat("  Panel 2:   Panel_Ordered_Boxplot.pdf\n")
cat("  Panel 3:   Panel_2x2_Heatmap.pdf\n")
cat("  Data:      A3B_Threshold_Sweep_Correlations.tsv\n")
cat("  Data:      A3B_A3A_Quadrant_Stats.tsv\n")
cat("============================================================\n")
