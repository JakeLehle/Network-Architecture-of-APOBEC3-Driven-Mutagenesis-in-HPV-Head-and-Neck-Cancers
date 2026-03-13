library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("Exploration: Individual A3 Family Member Roles in SBS2\n")
cat("============================================================\n\n")

# ============================================================
# STEP 1: Load and prepare data (same as other scripts)
# ============================================================
cat("STEP 1: Loading and preparing data\n")
cat("------------------------------------------------------------\n")

expression_data <- fread("A3_Expression_FPKM_UQ_Matched.tsv")
mutation_data   <- fread("Mutation_Signatures_Matched.tsv")

# Filter HNSC tumors
hnsc_expression <- expression_data %>%
  filter((Project_ID == "TCGA-HNSC" | Project_ID == "TCGA-HNSCC") &
         Tissue_Type == "Tumor")

hnsc_entity_ids <- hnsc_expression$Entity_ID

hnsc_mutation <- mutation_data %>%
  filter(TCGA_Gene_Expression_Entity_ID %in% hnsc_entity_ids)

common_ids <- intersect(hnsc_expression$Entity_ID,
                        hnsc_mutation$TCGA_Gene_Expression_Entity_ID)

# Build patient-level dataframe
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
cat("HNSC tumor patients:", n_tumors, "\n\n")

# ============================================================
# STEP 2: Raw and partial correlations
# ============================================================
cat("STEP 2: Computing raw and partial correlations\n")
cat("------------------------------------------------------------\n")

a3_genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3H")
a3_labels <- c("A3A", "A3B", "A3C", "A3H")

# --- Raw Spearman correlations ---
raw_cors <- sapply(a3_genes, function(gene) {
  cor(patient_data[[gene]], patient_data$SBS2, method = "spearman")
})

raw_pvals <- sapply(a3_genes, function(gene) {
  cor.test(patient_data[[gene]], patient_data$SBS2, method = "spearman")$p.value
})

# --- Partial Spearman correlations ---
# For each A3, correlate its residuals (after regressing out the other 3 A3s)
# with SBS2 residuals (after regressing out the same 3 A3s).
# This isolates each A3's independent contribution to SBS2.

partial_cors <- numeric(length(a3_genes))
partial_pvals <- numeric(length(a3_genes))
names(partial_cors) <- a3_genes
names(partial_pvals) <- a3_genes

for(i in seq_along(a3_genes)) {
  target_gene <- a3_genes[i]
  control_genes <- setdiff(a3_genes, target_gene)

  # Regress target A3 on the other three A3s → residuals
  formula_gene <- as.formula(paste(target_gene, "~",
                                   paste(control_genes, collapse = " + ")))
  resid_gene <- residuals(lm(formula_gene, data = patient_data))

  # Regress SBS2 on the other three A3s → residuals
  formula_sbs2 <- as.formula(paste("SBS2 ~",
                                    paste(control_genes, collapse = " + ")))
  resid_sbs2 <- residuals(lm(formula_sbs2, data = patient_data))

  # Spearman correlation of the residuals
  test <- cor.test(resid_gene, resid_sbs2, method = "spearman")
  partial_cors[i] <- test$estimate
  partial_pvals[i] <- test$p.value
}

# Build summary table
cor_summary <- data.frame(
  Gene = a3_labels,
  Raw_Spearman = round(raw_cors, 4),
  Raw_pvalue = raw_pvals,
  Partial_Spearman = round(partial_cors, 4),
  Partial_pvalue = partial_pvals,
  stringsAsFactors = FALSE
)

cat("\nCorrelation summary:\n")
print(cor_summary)
cat("\n")

# Save
write.table(cor_summary, "HNSC_A3_Raw_vs_Partial_Correlations.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# STEP 3: Multivariate regression
# ============================================================
cat("STEP 3: Multivariate regression (SBS2 ~ A3A + A3B + A3C + A3H)\n")
cat("------------------------------------------------------------\n")

# Standardize predictors so coefficients are comparable
patient_scaled <- patient_data %>%
  mutate(across(all_of(a3_genes), ~ scale(.)[,1], .names = "{.col}_z"))

model <- lm(SBS2 ~ APOBEC3A_z + APOBEC3B_z + APOBEC3C_z + APOBEC3H_z,
            data = patient_scaled)

cat("\nRegression summary:\n")
print(summary(model))

# Extract coefficients for plotting
coef_df <- data.frame(
  Gene = a3_labels,
  Beta = coef(model)[2:5],
  SE = summary(model)$coefficients[2:5, "Std. Error"],
  pvalue = summary(model)$coefficients[2:5, "Pr(>|t|)"],
  stringsAsFactors = FALSE
)
coef_df$CI_lower <- coef_df$Beta - 1.96 * coef_df$SE
coef_df$CI_upper <- coef_df$Beta + 1.96 * coef_df$SE
coef_df$Significant <- ifelse(coef_df$pvalue < 0.05, "Yes", "No")

cat("\nStandardized coefficients:\n")
print(coef_df)
cat("\nModel R-squared:", round(summary(model)$r.squared, 4), "\n")
cat("Model adjusted R-squared:", round(summary(model)$adj.r.squared, 4), "\n\n")

# Save
write.table(coef_df, "HNSC_A3_Multivariate_Regression_Coefficients.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# STEP 4: Conditional analysis — split by A3B and A3A status
# ============================================================
cat("STEP 4: Conditional correlation analysis\n")
cat("------------------------------------------------------------\n")

# Split by A3B median
median_a3b <- median(patient_data$APOBEC3B)
patient_data <- patient_data %>%
  mutate(A3B_status = ifelse(APOBEC3B >= median_a3b, "A3B High", "A3B Low"))

# A3A vs SBS2 correlation within each A3B group
for(grp in c("A3B High", "A3B Low")) {
  subset <- patient_data %>% filter(A3B_status == grp)
  test <- cor.test(subset$APOBEC3A, subset$SBS2, method = "spearman")
  cat(grp, "group (n =", nrow(subset), "):\n")
  cat("  A3A vs SBS2 Spearman rho =", round(test$estimate, 4),
      ", p =", format(test$p.value, scientific = TRUE, digits = 3), "\n")
}

# Split by A3A median
median_a3a <- median(patient_data$APOBEC3A)
patient_data <- patient_data %>%
  mutate(A3A_status = ifelse(APOBEC3A >= median_a3a, "A3A High", "A3A Low"))

# A3B vs SBS2 correlation within each A3A group
cat("\n")
for(grp in c("A3A High", "A3A Low")) {
  subset <- patient_data %>% filter(A3A_status == grp)
  test <- cor.test(subset$APOBEC3B, subset$SBS2, method = "spearman")
  cat(grp, "group (n =", nrow(subset), "):\n")
  cat("  A3B vs SBS2 Spearman rho =", round(test$estimate, 4),
      ", p =", format(test$p.value, scientific = TRUE, digits = 3), "\n")
}
cat("\n")

# ============================================================
# STEP 5: Generate all panels
# ============================================================
cat("STEP 5: Generating figure panels\n")
cat("------------------------------------------------------------\n")

# --- Shared theme ---
theme_panel <- theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank()
  )

# =============================================
# PANEL A: Four individual A3 vs SBS2 scatters
# =============================================
cat("  Building Panel A: Individual A3 vs SBS2 scatters...\n")

# Reshape for faceting
plot_data_long <- patient_data %>%
  tidyr::pivot_longer(cols = all_of(a3_genes),
                      names_to = "A3_gene",
                      values_to = "Expression") %>%
  mutate(A3_label = case_when(
    A3_gene == "APOBEC3A" ~ "A3A",
    A3_gene == "APOBEC3B" ~ "A3B",
    A3_gene == "APOBEC3C" ~ "A3C",
    A3_gene == "APOBEC3H" ~ "A3H"
  ))

# Annotation dataframe with raw Spearman rho per gene
annot_raw <- data.frame(
  A3_label = a3_labels,
  rho = round(raw_cors, 3),
  pval = raw_pvals,
  stringsAsFactors = FALSE
) %>%
  mutate(label = paste0("rho == ", rho))

# Set factor order so A3B comes first (most important)
plot_data_long$A3_label <- factor(plot_data_long$A3_label,
                                   levels = c("A3A", "A3B", "A3C", "A3H"))
annot_raw$A3_label <- factor(annot_raw$A3_label,
                              levels = c("A3A", "A3B", "A3C", "A3H"))

panel_A <- ggplot(plot_data_long, aes(x = Expression, y = SBS2)) +
  geom_point(shape = 21, size = 3, alpha = 0.5,
             fill = "grey60", color = "black", stroke = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#E41A1C",
              linewidth = 1, linetype = "solid", alpha = 0.2) +
  facet_wrap(~ A3_label, scales = "free_x", nrow = 1) +
  geom_text(data = annot_raw,
            aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
            size = 6, parse = TRUE, color = "#E41A1C") +
  labs(x = "Expression (FPKM-UQ)",
       y = "SBS2 Weight",
       title = "Individual A3 Family Members vs SBS2",
       subtitle = paste0("HNSC tumors (n = ", n_tumors, ")")) +
  theme_panel

# =============================================
# PANEL B: Raw vs Partial correlation bar plot
# =============================================
cat("  Building Panel B: Raw vs Partial correlation comparison...\n")

cor_plot_data <- data.frame(
  Gene = rep(a3_labels, 2),
  Correlation = c(raw_cors, partial_cors),
  Type = rep(c("Raw Spearman", "Partial Spearman"), each = 4),
  pvalue = c(raw_pvals, partial_pvals),
  stringsAsFactors = FALSE
)
cor_plot_data$Gene <- factor(cor_plot_data$Gene,
                              levels = c("A3A", "A3B", "A3C", "A3H"))
cor_plot_data$Significant <- ifelse(cor_plot_data$pvalue < 0.05, "p < 0.05", "n.s.")

panel_B <- ggplot(cor_plot_data,
                  aes(x = Gene, y = Correlation, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = ifelse(Significant == "p < 0.05", "*", ""),
                y = Correlation + sign(Correlation) * 0.02),
            position = position_dodge(width = 0.7),
            size = 8, vjust = 0.5) +
  scale_fill_manual(values = c("Raw Spearman" = "#5e81ac",
                                "Partial Spearman" = "#bf616a"),
                    name = "Correlation Type") +
  labs(x = NULL,
       y = "Spearman Correlation with SBS2",
       title = "Raw vs Partial Correlation",
       subtitle = "Partial: controlling for other A3 family members") +
  theme_panel +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"))

# =============================================
# PANEL C: Multivariate regression coefficients
# =============================================
cat("  Building Panel C: Regression coefficient plot...\n")

coef_df$Gene <- factor(coef_df$Gene, levels = c("A3H", "A3C", "A3A", "A3B"))

panel_C <- ggplot(coef_df, aes(x = Beta, y = Gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.25, linewidth = 0.8, color = "grey30") +
  geom_point(aes(fill = Significant), shape = 21, size = 6,
             color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("Yes" = "#E41A1C", "No" = "grey70"),
                    name = "p < 0.05") +
  labs(x = "Standardized Beta Coefficient",
       y = NULL,
       title = "Multivariate Regression",
       subtitle = paste0("SBS2 ~ A3A + A3B + A3C + A3H  |  adj. R² = ",
                         round(summary(model)$adj.r.squared, 3))) +
  theme_panel +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"))

# =============================================
# PANEL D: A3A vs SBS2 faceted by A3B status
# =============================================
cat("  Building Panel D: A3A vs SBS2 by A3B status...\n")

# Compute per-facet correlations for annotation
annot_D <- patient_data %>%
  group_by(A3B_status) %>%
  summarise(
    rho = round(cor(APOBEC3A, SBS2, method = "spearman"), 3),
    pval = cor.test(APOBEC3A, SBS2, method = "spearman")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("rho == ", rho, "~~(n == ", n, ")"))

patient_data$A3B_status <- factor(patient_data$A3B_status,
                                   levels = c("A3B Low", "A3B High"))
annot_D$A3B_status <- factor(annot_D$A3B_status,
                              levels = c("A3B Low", "A3B High"))

panel_D <- ggplot(patient_data, aes(x = APOBEC3A, y = SBS2)) +
  geom_point(aes(fill = A3B_status), shape = 21, size = 3,
             alpha = 0.5, color = "black", stroke = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#E41A1C",
              linewidth = 1, alpha = 0.2) +
  facet_wrap(~ A3B_status) +
  scale_fill_manual(values = c("A3B Low" = "#f4f1bb", "A3B High" = "#ed6a5a"),
                    guide = "none") +
  geom_text(data = annot_D,
            aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
            size = 5.5, parse = TRUE, color = "#E41A1C") +
  labs(x = "APOBEC3A Expression (FPKM-UQ)",
       y = "SBS2 Weight",
       title = "A3A vs SBS2 — Conditional on A3B Expression",
       subtitle = paste0("Split at median A3B = ",
                         round(median_a3b, 2), " FPKM-UQ")) +
  theme_panel

# =============================================
# PANEL E: A3B vs SBS2 faceted by A3A status
# =============================================
cat("  Building Panel E: A3B vs SBS2 by A3A status...\n")

annot_E <- patient_data %>%
  group_by(A3A_status) %>%
  summarise(
    rho = round(cor(APOBEC3B, SBS2, method = "spearman"), 3),
    pval = cor.test(APOBEC3B, SBS2, method = "spearman")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("rho == ", rho, "~~(n == ", n, ")"))

patient_data$A3A_status <- factor(patient_data$A3A_status,
                                   levels = c("A3A Low", "A3A High"))
annot_E$A3A_status <- factor(annot_E$A3A_status,
                              levels = c("A3A Low", "A3A High"))

panel_E <- ggplot(patient_data, aes(x = APOBEC3B, y = SBS2)) +
  geom_point(aes(fill = A3A_status), shape = 21, size = 3,
             alpha = 0.5, color = "black", stroke = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#E41A1C",
              linewidth = 1, alpha = 0.2) +
  facet_wrap(~ A3A_status) +
  scale_fill_manual(values = c("A3A Low" = "#dcedc8", "A3A High" = "#5e81ac"),
                    guide = "none") +
  geom_text(data = annot_E,
            aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
            size = 5.5, parse = TRUE, color = "#E41A1C") +
  labs(x = "APOBEC3B Expression (FPKM-UQ)",
       y = "SBS2 Weight",
       title = "A3B vs SBS2 — Conditional on A3A Expression",
       subtitle = paste0("Split at median A3A = ",
                         round(median_a3a, 2), " FPKM-UQ")) +
  theme_panel

# =============================================
# PANEL F: A3A vs A3B colored by SBS2
# =============================================
cat("  Building Panel F: A3A vs A3B colored by SBS2...\n")

panel_F <- ggplot(patient_data, aes(x = APOBEC3B, y = APOBEC3A)) +
  geom_point(aes(fill = SBS2), shape = 21, size = 4,
             alpha = 0.7, color = "black", stroke = 0.5) +
  scale_fill_viridis(option = "viridis", name = "SBS2\nWeight",
                     guide = guide_colorbar(barwidth = 1.5, barheight = 10)) +
  geom_hline(yintercept = median_a3a, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = median_a3b, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("(n = ", n_tumors, ")"),
           hjust = 1.1, vjust = 1.5, size = 6) +
  labs(x = "APOBEC3B Expression (FPKM-UQ)",
       y = "APOBEC3A Expression (FPKM-UQ)",
       title = "A3B vs A3A — Colored by SBS2 Mutational Weight",
       subtitle = "Dashed lines = median expression") +
  theme_panel

# ============================================================
# STEP 6: Assemble and save
# ============================================================
cat("\nSTEP 6: Assembling composite figure\n")
cat("------------------------------------------------------------\n")

# Layout:
#   Row 1: Panel A (full width — 4 scatter panels)
#   Row 2: Panel B (left) + Panel C (right)
#   Row 3: Panel D (left) + Panel E (right)
#   Row 4: Panel F (full width)

composite <- (panel_A) /
             (panel_B | panel_C) /
             (panel_D | panel_E) /
             (panel_F) +
             plot_annotation(
               title = "Dissecting Individual A3 Family Member Contributions to SBS2 Mutagenesis",
               subtitle = paste0("TCGA-HNSC Tumors (n = ", n_tumors, ")"),
               theme = theme(
                 plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = 18, hjust = 0.5, color = "grey30")
               )
             ) +
             plot_layout(heights = c(1, 1, 1, 1.2))

# Save composite
composite_pdf <- "Exploration_A3_Individual_Contributions_to_SBS2.pdf"
ggsave(composite_pdf, composite, width = 20, height = 28, units = "in")
cat("Saved composite PDF:", composite_pdf, "\n")

composite_png <- "Exploration_A3_Individual_Contributions_to_SBS2.png"
ggsave(composite_png, composite, width = 20, height = 28, units = "in", dpi = 200)
cat("Saved composite PNG:", composite_png, "\n")

# Also save individual panels for flexibility
cat("\nSaving individual panels...\n")
ggsave("Panel_A_Individual_A3_Scatters.pdf", panel_A, width = 18, height = 5)
ggsave("Panel_B_Raw_vs_Partial_Correlation.pdf", panel_B, width = 8, height = 6)
ggsave("Panel_C_Regression_Coefficients.pdf", panel_C, width = 8, height = 6)
ggsave("Panel_D_A3A_vs_SBS2_by_A3B_status.pdf", panel_D, width = 14, height = 6)
ggsave("Panel_E_A3B_vs_SBS2_by_A3A_status.pdf", panel_E, width = 14, height = 6)
ggsave("Panel_F_A3A_vs_A3B_colored_SBS2.pdf", panel_F, width = 10, height = 8)
cat("All individual panels saved.\n")

# ============================================================
# Summary
# ============================================================
cat("\n============================================================\n")
cat("EXPLORATION COMPLETE\n")
cat("============================================================\n")

cat("\nKey findings:\n\n")

cat("POINT 1 — Ruling out A3C/A3H:\n")
cat("  A3C raw rho =", round(raw_cors["APOBEC3C"], 4),
    "→ partial rho =", round(partial_cors["APOBEC3C"], 4), "\n")
cat("  A3H raw rho =", round(raw_cors["APOBEC3H"], 4),
    "→ partial rho =", round(partial_cors["APOBEC3H"], 4), "\n")
cat("  Interpretation: If partial correlations collapse, A3C/A3H\n")
cat("  associations with SBS2 are driven by co-expression with A3A/A3B.\n\n")

cat("POINT 2 — A3B dominates:\n")
cat("  A3B standardized beta =", round(coef_df$Beta[coef_df$Gene == "A3B"], 4),
    "(p =", format(coef_df$pvalue[coef_df$Gene == "A3B"], scientific = TRUE, digits = 3), ")\n")
cat("  A3A standardized beta =", round(coef_df$Beta[coef_df$Gene == "A3A"], 4),
    "(p =", format(coef_df$pvalue[coef_df$Gene == "A3A"], scientific = TRUE, digits = 3), ")\n\n")

cat("POINT 3 — A3A conditional on A3B:\n")
for(grp in c("A3B Low", "A3B High")) {
  subset <- patient_data %>% filter(A3B_status == grp)
  rho <- round(cor(subset$APOBEC3A, subset$SBS2, method = "spearman"), 4)
  cat("  A3A vs SBS2 in", grp, "group: rho =", rho, "\n")
}
cat("  Compare to A3B vs SBS2 which should maintain correlation\n")
cat("  regardless of A3A status (asymmetric relationship).\n\n")

cat("POINT 4 — Synthesis:\n")
cat("  Panel F shows A3A vs A3B colored by SBS2.\n")
cat("  High SBS2 (bright) should cluster in the A3B-high region,\n")
cat("  with the highest values where both A3A and A3B are elevated.\n\n")

cat("OUTPUT FILES:\n")
cat("  Composite: ", composite_pdf, " / ", composite_png, "\n")
cat("  Individual panels: Panel_A through Panel_F PDFs\n")
cat("  Data: HNSC_A3_Raw_vs_Partial_Correlations.tsv\n")
cat("  Data: HNSC_A3_Multivariate_Regression_Coefficients.tsv\n")
cat("============================================================\n")
