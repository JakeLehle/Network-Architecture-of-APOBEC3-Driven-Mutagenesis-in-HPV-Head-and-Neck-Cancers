library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("Diagnostic: Non-Linear A3B and A3A Contributions to SBS2\n")
cat("============================================================\n\n")

# ============================================================
# DATA LOADING (identical to previous exploration script)
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
a3_genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3H")
a3_labels <- c("A3A", "A3B", "A3C", "A3H")

cat("HNSC tumor patients:", n_tumors, "\n\n")

# Shared theme
theme_panel <- theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey30"),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )

# Medians for splitting
median_a3a <- median(patient_data$APOBEC3A)
median_a3b <- median(patient_data$APOBEC3B)

cat("Median A3A:", round(median_a3a, 4), "\n")
cat("Median A3B:", round(median_a3b, 4), "\n\n")


# ############################################################
#
#  OPTION 1: Threshold-Based Spearman Correlation
#  (A3B rising phase vs plateau phase)
#
# ############################################################
cat("============================================================\n")
cat("OPTION 1: Threshold-Based Correlation Analysis\n")
cat("============================================================\n\n")

# Identify A3B's plateau region by looking at SBS2 as a function of A3B.
# We test multiple A3B thresholds and compute A3B-SBS2 correlation in the
# below-threshold (rising) and above-threshold (plateau) subsets.

a3b_quantiles <- quantile(patient_data$APOBEC3B,
                           probs = seq(0.2, 0.8, by = 0.1))

cat("Testing A3B thresholds at quantiles:\n")

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

for(q in seq(0.2, 0.8, by = 0.1)) {
  thresh <- quantile(patient_data$APOBEC3B, q)
  below <- patient_data %>% filter(APOBEC3B < thresh)
  above <- patient_data %>% filter(APOBEC3B >= thresh)

  # Need at least 10 points per group
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

cat("\nA3B threshold analysis results:\n")
print(threshold_results)
cat("\n")

write.table(threshold_results, "Diag_Option1_A3B_Threshold_Correlations.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Option 1 Plot ---
# Show A3B vs SBS2 scatter, colored by whether below/above the median split,
# with per-group correlation lines and rho values.
# Also include small inset-style panels at different thresholds.

# Use median as the primary visual threshold
patient_data <- patient_data %>%
  mutate(A3B_phase = ifelse(APOBEC3B < median_a3b, "Rising Phase", "Plateau Phase"))

patient_data$A3B_phase <- factor(patient_data$A3B_phase,
                                  levels = c("Rising Phase", "Plateau Phase"))

# Per-phase correlations for annotation
annot_opt1 <- patient_data %>%
  group_by(A3B_phase) %>%
  summarise(
    rho = round(cor(APOBEC3B, SBS2, method = "spearman"), 3),
    pval = cor.test(APOBEC3B, SBS2, method = "spearman")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("rho == ", rho, "~~(n == ", n, ")"))

opt1_plot <- ggplot(patient_data, aes(x = APOBEC3B, y = SBS2)) +
  geom_point(aes(fill = A3B_phase), shape = 21, size = 3,
             alpha = 0.5, color = "black", stroke = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#E41A1C",
              linewidth = 1, alpha = 0.2) +
  geom_vline(xintercept = median_a3b, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  facet_wrap(~ A3B_phase, scales = "free_x") +
  scale_fill_manual(values = c("Rising Phase" = "#5e81ac",
                                "Plateau Phase" = "#bf616a"),
                    guide = "none") +
  geom_text(data = annot_opt1,
            aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
            size = 5, parse = TRUE, color = "#E41A1C") +
  labs(x = "APOBEC3B Expression (FPKM-UQ)",
       y = "SBS2 Weight",
       title = "Option 1: A3B Rising Phase vs Plateau Phase",
       subtitle = paste0("Split at median A3B = ", round(median_a3b, 2),
                         " | Does correlation drop after plateau?")) +
  theme_panel

# Also make a summary plot showing rho at each threshold
opt1_threshold_plot <- ggplot(threshold_results, aes(x = Quantile)) +
  geom_line(aes(y = Rho_below, color = "Below threshold"), linewidth = 1) +
  geom_point(aes(y = Rho_below, color = "Below threshold"), size = 3) +
  geom_line(aes(y = Rho_above, color = "Above threshold"), linewidth = 1) +
  geom_point(aes(y = Rho_above, color = "Above threshold"), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Below threshold" = "#5e81ac",
                                 "Above threshold" = "#bf616a"),
                     name = "A3B Group") +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.1),
                     labels = paste0(seq(20, 80, 10), "th")) +
  labs(x = "A3B Quantile Threshold",
       y = "Spearman rho (A3B vs SBS2)",
       title = "A3B-SBS2 Correlation Across Thresholds",
       subtitle = "Rising phase should show stronger rho than plateau phase") +
  theme_panel +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12))

opt1_combined <- opt1_plot / opt1_threshold_plot +
  plot_annotation(title = "OPTION 1: Threshold-Based Correlation",
                  theme = theme(plot.title = element_text(size = 22, face = "bold",
                                                          hjust = 0.5)))

ggsave("Diag_Option1_Threshold_Correlation.pdf", opt1_combined,
       width = 16, height = 14)
ggsave("Diag_Option1_Threshold_Correlation.png", opt1_combined,
       width = 16, height = 14, dpi = 200)
cat("Option 1 figures saved.\n\n")


# ############################################################
#
#  OPTION 2: Log-Transformed Regression & Partial Correlations
#
# ############################################################
cat("============================================================\n")
cat("OPTION 2: Log-Transform A3B for Non-Linear Modeling\n")
cat("============================================================\n\n")

# Add log-transformed columns
patient_data <- patient_data %>%
  mutate(
    log_A3A = log1p(APOBEC3A),
    log_A3B = log1p(APOBEC3B),
    log_A3C = log1p(APOBEC3C),
    log_A3H = log1p(APOBEC3H)
  )

# --- 2A: Compare linear vs log-transformed regression ---
model_linear <- lm(SBS2 ~ APOBEC3A + APOBEC3B + APOBEC3C + APOBEC3H,
                   data = patient_data)
model_log    <- lm(SBS2 ~ APOBEC3A + log_A3B + APOBEC3C + APOBEC3H,
                   data = patient_data)
model_logall <- lm(SBS2 ~ log_A3A + log_A3B + log_A3C + log_A3H,
                   data = patient_data)

cat("Model comparison (adjusted R²):\n")
cat("  Linear (all raw):       ", round(summary(model_linear)$adj.r.squared, 4), "\n")
cat("  Log A3B only:           ", round(summary(model_log)$adj.r.squared, 4), "\n")
cat("  Log all A3s:            ", round(summary(model_logall)$adj.r.squared, 4), "\n\n")

# AIC comparison
cat("AIC comparison:\n")
cat("  Linear (all raw):       ", round(AIC(model_linear), 2), "\n")
cat("  Log A3B only:           ", round(AIC(model_log), 2), "\n")
cat("  Log all A3s:            ", round(AIC(model_logall), 2), "\n\n")

# --- 2B: Standardized coefficients from log-all model ---
patient_scaled_log <- patient_data %>%
  mutate(across(c(log_A3A, log_A3B, log_A3C, log_A3H),
                ~ scale(.)[,1], .names = "{.col}_z"))

model_logall_z <- lm(SBS2 ~ log_A3A_z + log_A3B_z + log_A3C_z + log_A3H_z,
                      data = patient_scaled_log)

cat("Log-transformed standardized regression:\n")
print(summary(model_logall_z))

coef_linear <- data.frame(
  Gene = a3_labels,
  Beta = coef(model_linear)[2:5],
  SE = summary(model_linear)$coefficients[2:5, "Std. Error"],
  pvalue = summary(model_linear)$coefficients[2:5, "Pr(>|t|)"],
  Model = "Linear",
  stringsAsFactors = FALSE
)

coef_log <- data.frame(
  Gene = a3_labels,
  Beta = coef(model_logall_z)[2:5],
  SE = summary(model_logall_z)$coefficients[2:5, "Std. Error"],
  pvalue = summary(model_logall_z)$coefficients[2:5, "Pr(>|t|)"],
  Model = "Log-transformed",
  stringsAsFactors = FALSE
)

coef_both <- rbind(coef_linear, coef_log)
coef_both$CI_lower <- coef_both$Beta - 1.96 * coef_both$SE
coef_both$CI_upper <- coef_both$Beta + 1.96 * coef_both$SE
coef_both$Significant <- ifelse(coef_both$pvalue < 0.05, "Yes", "No")
coef_both$Gene <- factor(coef_both$Gene, levels = c("A3H", "A3C", "A3A", "A3B"))

write.table(coef_both, "Diag_Option2_Linear_vs_Log_Coefficients.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- 2C: Partial correlations using log-transformed values ---
log_genes <- c("log_A3A", "log_A3B", "log_A3C", "log_A3H")

partial_cors_log <- numeric(4)
partial_pvals_log <- numeric(4)
names(partial_cors_log) <- a3_labels
names(partial_pvals_log) <- a3_labels

for(i in 1:4) {
  target <- log_genes[i]
  controls <- setdiff(log_genes, target)

  resid_gene <- residuals(lm(as.formula(paste(target, "~",
                              paste(controls, collapse = " + "))),
                              data = patient_data))
  resid_sbs2 <- residuals(lm(as.formula(paste("SBS2 ~",
                              paste(controls, collapse = " + "))),
                              data = patient_data))

  test <- cor.test(resid_gene, resid_sbs2, method = "spearman")
  partial_cors_log[i] <- test$estimate
  partial_pvals_log[i] <- test$p.value
}

cat("\nPartial correlations (log-transformed):\n")
for(i in 1:4) {
  cat("  ", a3_labels[i], ": rho =", round(partial_cors_log[i], 4),
      ", p =", format(partial_pvals_log[i], scientific = TRUE, digits = 3), "\n")
}

# Also compute raw partial correlations (linear) for comparison
raw_cors_vec <- sapply(a3_genes, function(gene) {
  cor(patient_data[[gene]], patient_data$SBS2, method = "spearman")
})

partial_cors_raw <- numeric(4)
partial_pvals_raw <- numeric(4)
names(partial_cors_raw) <- a3_labels
names(partial_pvals_raw) <- a3_labels

for(i in 1:4) {
  target <- a3_genes[i]
  controls <- setdiff(a3_genes, target)

  resid_gene <- residuals(lm(as.formula(paste(target, "~",
                              paste(controls, collapse = " + "))),
                              data = patient_data))
  resid_sbs2 <- residuals(lm(as.formula(paste("SBS2 ~",
                              paste(controls, collapse = " + "))),
                              data = patient_data))

  test <- cor.test(resid_gene, resid_sbs2, method = "spearman")
  partial_cors_raw[i] <- test$estimate
  partial_pvals_raw[i] <- test$p.value
}

# --- Option 2 Plots ---

# Plot 2A: Side-by-side coefficient plots (linear vs log)
opt2_coef_plot <- ggplot(coef_both, aes(x = Beta, y = Gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Model),
                 height = 0.2, linewidth = 0.8,
                 position = position_dodge(width = 0.5)) +
  geom_point(aes(fill = Significant, color = Model), shape = 21, size = 5,
             stroke = 0.8,
             position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("Yes" = "#E41A1C", "No" = "grey70"),
                    name = "p < 0.05") +
  scale_color_manual(values = c("Linear" = "#5e81ac", "Log-transformed" = "#bf616a")) +
  labs(x = "Standardized Beta Coefficient",
       y = NULL,
       title = "Linear vs Log-Transformed Regression Coefficients",
       subtitle = paste0("Linear adj.R² = ",
                         round(summary(model_linear)$adj.r.squared, 3),
                         "  |  Log adj.R² = ",
                         round(summary(model_logall_z)$adj.r.squared, 3))) +
  theme_panel +
  theme(legend.position = "bottom")

# Plot 2B: Partial correlations — raw vs log
partial_comparison <- data.frame(
  Gene = rep(a3_labels, 3),
  Correlation = c(raw_cors_vec, partial_cors_raw, partial_cors_log),
  Type = rep(c("Raw Spearman", "Partial (Linear)", "Partial (Log)"), each = 4),
  pvalue = c(rep(NA, 4), partial_pvals_raw, partial_pvals_log),
  stringsAsFactors = FALSE
)
partial_comparison$Gene <- factor(partial_comparison$Gene,
                                   levels = c("A3A", "A3B", "A3C", "A3H"))
partial_comparison$Type <- factor(partial_comparison$Type,
                                   levels = c("Raw Spearman",
                                              "Partial (Linear)",
                                              "Partial (Log)"))

opt2_partial_plot <- ggplot(partial_comparison,
                            aes(x = Gene, y = Correlation, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Raw Spearman" = "#a3be8c",
                                "Partial (Linear)" = "#5e81ac",
                                "Partial (Log)" = "#bf616a"),
                    name = "Correlation Type") +
  labs(x = NULL,
       y = "Spearman Correlation with SBS2",
       title = "Raw vs Partial Correlations: Linear vs Log-Transformed",
       subtitle = "Does log-transforming A3 expression recover A3B's contribution?") +
  theme_panel +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12))

# Plot 2C: Visual of log transform effect on A3B scatter
opt2_scatter <- ggplot(patient_data, aes(y = SBS2)) +
  geom_point(aes(x = APOBEC3B), shape = 21, size = 3,
             alpha = 0.5, fill = "#5e81ac", color = "black", stroke = 0.5) +
  geom_smooth(aes(x = APOBEC3B), method = "lm", se = TRUE,
              color = "#E41A1C", linewidth = 1, alpha = 0.2) +
  labs(x = "APOBEC3B Expression (FPKM-UQ)", y = "SBS2 Weight",
       title = "A3B vs SBS2 (Raw Scale)") +
  theme_panel

opt2_scatter_log <- ggplot(patient_data, aes(y = SBS2)) +
  geom_point(aes(x = log_A3B), shape = 21, size = 3,
             alpha = 0.5, fill = "#bf616a", color = "black", stroke = 0.5) +
  geom_smooth(aes(x = log_A3B), method = "lm", se = TRUE,
              color = "#E41A1C", linewidth = 1, alpha = 0.2) +
  labs(x = "log(APOBEC3B + 1)", y = "SBS2 Weight",
       title = "A3B vs SBS2 (Log Scale)") +
  theme_panel

# Annotate each with rho
rho_raw_a3b <- round(cor(patient_data$APOBEC3B, patient_data$SBS2, method = "spearman"), 3)
rho_log_a3b <- round(cor(patient_data$log_A3B, patient_data$SBS2, method = "spearman"), 3)

opt2_scatter <- opt2_scatter +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho == ", rho_raw_a3b), parse = TRUE,
           size = 6, color = "#E41A1C")

opt2_scatter_log <- opt2_scatter_log +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho == ", rho_log_a3b), parse = TRUE,
           size = 6, color = "#E41A1C")

opt2_combined <- (opt2_scatter | opt2_scatter_log) /
                 (opt2_partial_plot | opt2_coef_plot) +
  plot_annotation(title = "OPTION 2: Log-Transformed Non-Linear Modeling",
                  theme = theme(plot.title = element_text(size = 22, face = "bold",
                                                          hjust = 0.5)))

ggsave("Diag_Option2_Log_Transform.pdf", opt2_combined,
       width = 18, height = 16)
ggsave("Diag_Option2_Log_Transform.png", opt2_combined,
       width = 18, height = 16, dpi = 200)
cat("\nOption 2 figures saved.\n\n")


# ############################################################
#
#  OPTION 3: Sequential (Type I) ANOVA — Variance Explained
#
# ############################################################
cat("============================================================\n")
cat("OPTION 3: Sequential ANOVA — Order of Variance Explained\n")
cat("============================================================\n\n")

# The key question: does A3B explain baseline SBS2 variance, and does
# A3A explain additional variance on top of A3B?
# We test BOTH orderings to see the asymmetry.

# --- Order 1: A3B first, then A3A ---
model_b_first <- lm(SBS2 ~ APOBEC3B + APOBEC3A + APOBEC3C + APOBEC3H,
                     data = patient_data)
anova_b_first <- anova(model_b_first)

# --- Order 2: A3A first, then A3B ---
model_a_first <- lm(SBS2 ~ APOBEC3A + APOBEC3B + APOBEC3C + APOBEC3H,
                     data = patient_data)
anova_a_first <- anova(model_a_first)

cat("Order 1: A3B entered first\n")
print(anova_b_first)
cat("\nOrder 2: A3A entered first\n")
print(anova_a_first)

# Extract sum of squares for each gene in each ordering
total_ss <- sum(anova_b_first$`Sum Sq`)

anova_summary <- data.frame(
  Gene = c("A3B", "A3A", "A3C", "A3H", "Residual",
           "A3A", "A3B", "A3C", "A3H", "Residual"),
  SumSq = c(anova_b_first$`Sum Sq`, anova_a_first$`Sum Sq`),
  PctVariance = c(anova_b_first$`Sum Sq` / total_ss * 100,
                  anova_a_first$`Sum Sq` / total_ss * 100),
  Pvalue = c(anova_b_first$`Pr(>F)`, anova_a_first$`Pr(>F)`),
  Order = c(rep("A3B first", 5), rep("A3A first", 5)),
  Entry = c(1:5, 1:5),
  stringsAsFactors = FALSE
)

cat("\n\nVariance explained summary:\n")
print(anova_summary)

write.table(anova_summary, "Diag_Option3_Sequential_ANOVA.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Also do log-transformed sequential ANOVA ---
model_b_first_log <- lm(SBS2 ~ log_A3B + log_A3A + log_A3C + log_A3H,
                         data = patient_data)
anova_b_first_log <- anova(model_b_first_log)

model_a_first_log <- lm(SBS2 ~ log_A3A + log_A3B + log_A3C + log_A3H,
                         data = patient_data)
anova_a_first_log <- anova(model_a_first_log)

total_ss_log <- sum(anova_b_first_log$`Sum Sq`)

cat("\n\nLog-transformed sequential ANOVA:\n")
cat("Order 1 (A3B first):\n")
print(anova_b_first_log)
cat("\nOrder 2 (A3A first):\n")
print(anova_a_first_log)

# --- Option 3 Plots ---

# Build a focused comparison: % variance for A3A and A3B under both orderings
# Linear model
var_data_linear <- data.frame(
  Gene = c("A3B (1st)", "A3A (2nd)", "A3A (1st)", "A3B (2nd)"),
  PctVar = c(
    anova_b_first$`Sum Sq`[1] / total_ss * 100,    # A3B entered 1st
    anova_b_first$`Sum Sq`[2] / total_ss * 100,    # A3A entered 2nd
    anova_a_first$`Sum Sq`[1] / total_ss * 100,    # A3A entered 1st
    anova_a_first$`Sum Sq`[2] / total_ss * 100     # A3B entered 2nd
  ),
  Ordering = c("A3B → A3A", "A3B → A3A", "A3A → A3B", "A3A → A3B"),
  Role = c("A3B", "A3A", "A3A", "A3B"),
  Model = "Linear",
  stringsAsFactors = FALSE
)

# Log model
var_data_log <- data.frame(
  Gene = c("A3B (1st)", "A3A (2nd)", "A3A (1st)", "A3B (2nd)"),
  PctVar = c(
    anova_b_first_log$`Sum Sq`[1] / total_ss_log * 100,
    anova_b_first_log$`Sum Sq`[2] / total_ss_log * 100,
    anova_a_first_log$`Sum Sq`[1] / total_ss_log * 100,
    anova_a_first_log$`Sum Sq`[2] / total_ss_log * 100
  ),
  Ordering = c("A3B → A3A", "A3B → A3A", "A3A → A3B", "A3A → A3B"),
  Role = c("A3B", "A3A", "A3A", "A3B"),
  Model = "Log-transformed",
  stringsAsFactors = FALSE
)

var_data <- rbind(var_data_linear, var_data_log)
var_data$Gene <- factor(var_data$Gene,
                         levels = c("A3B (1st)", "A3A (2nd)",
                                    "A3A (1st)", "A3B (2nd)"))

opt3_plot <- ggplot(var_data, aes(x = Gene, y = PctVar, fill = Role)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.4, width = 0.7) +
  facet_grid(Model ~ Ordering, scales = "free_x") +
  scale_fill_manual(values = c("A3A" = "#5e81ac", "A3B" = "#bf616a"),
                    name = "Gene") +
  labs(x = NULL,
       y = "% Total Variance Explained",
       title = "Sequential ANOVA: Order of Entry Matters",
       subtitle = "Does A3B capture baseline SBS2 that A3A then amplifies?") +
  theme_panel +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

# Incremental R² plot
incremental_data <- data.frame(
  Model = c("A3B alone", "+ A3A", "+ A3C + A3H",
            "A3A alone", "+ A3B", "+ A3C + A3H",
            "log(A3B) alone", "+ log(A3A)", "+ log(A3C) + log(A3H)",
            "log(A3A) alone", "+ log(A3B)", "+ log(A3C) + log(A3H)"),
  R2 = c(
    summary(lm(SBS2 ~ APOBEC3B, data = patient_data))$r.squared,
    summary(lm(SBS2 ~ APOBEC3B + APOBEC3A, data = patient_data))$r.squared,
    summary(model_b_first)$r.squared,
    summary(lm(SBS2 ~ APOBEC3A, data = patient_data))$r.squared,
    summary(lm(SBS2 ~ APOBEC3A + APOBEC3B, data = patient_data))$r.squared,
    summary(model_a_first)$r.squared,
    summary(lm(SBS2 ~ log_A3B, data = patient_data))$r.squared,
    summary(lm(SBS2 ~ log_A3B + log_A3A, data = patient_data))$r.squared,
    summary(model_b_first_log)$r.squared,
    summary(lm(SBS2 ~ log_A3A, data = patient_data))$r.squared,
    summary(lm(SBS2 ~ log_A3A + log_A3B, data = patient_data))$r.squared,
    summary(model_a_first_log)$r.squared
  ),
  Ordering = rep(c("A3B first (linear)", "A3A first (linear)",
                    "A3B first (log)", "A3A first (log)"), each = 3),
  Step = rep(1:3, 4),
  stringsAsFactors = FALSE
)

incremental_data$Ordering <- factor(incremental_data$Ordering,
                                     levels = c("A3B first (linear)",
                                                "A3A first (linear)",
                                                "A3B first (log)",
                                                "A3A first (log)"))

opt3_incremental <- ggplot(incremental_data,
                           aes(x = Step, y = R2, color = Ordering, group = Ordering)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = round(R2, 3)), vjust = -1, size = 4) +
  scale_x_continuous(breaks = 1:3,
                     labels = c("Gene 1\nalone",
                                "+ Gene 2",
                                "+ A3C + A3H")) +
  scale_color_manual(values = c("A3B first (linear)" = "#bf616a",
                                 "A3A first (linear)" = "#5e81ac",
                                 "A3B first (log)" = "#d08770",
                                 "A3A first (log)" = "#81a1c1")) +
  labs(x = "Model Step",
       y = "R²",
       title = "Incremental R² — Building the Model Step by Step",
       subtitle = "How much does each gene add when entered sequentially?") +
  theme_panel +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 11))

opt3_combined <- opt3_plot / opt3_incremental +
  plot_annotation(title = "OPTION 3: Sequential Variance Decomposition",
                  theme = theme(plot.title = element_text(size = 22, face = "bold",
                                                          hjust = 0.5)))

ggsave("Diag_Option3_Sequential_ANOVA.pdf", opt3_combined,
       width = 16, height = 16)
ggsave("Diag_Option3_Sequential_ANOVA.png", opt3_combined,
       width = 16, height = 16, dpi = 200)
cat("\nOption 3 figures saved.\n\n")


# ############################################################
#
#  OPTION 4: A3B as Binary Gate, A3A as Continuous Driver
#
# ############################################################
cat("============================================================\n")
cat("OPTION 4: A3B as Gate / A3A as Driver Conceptual Model\n")
cat("============================================================\n\n")

# Classify tumors into a 2x2 framework
patient_data <- patient_data %>%
  mutate(
    A3B_gate = ifelse(APOBEC3B >= median_a3b, "A3B Present", "A3B Absent"),
    A3A_level = ifelse(APOBEC3A >= median_a3a, "A3A High", "A3A Low"),
    Quadrant = paste(A3B_gate, "/", A3A_level)
  )

patient_data$A3B_gate <- factor(patient_data$A3B_gate,
                                 levels = c("A3B Absent", "A3B Present"))
patient_data$A3A_level <- factor(patient_data$A3A_level,
                                  levels = c("A3A Low", "A3A High"))

# SBS2 statistics by quadrant
quadrant_stats <- patient_data %>%
  group_by(A3B_gate, A3A_level) %>%
  summarise(
    n = n(),
    mean_SBS2 = round(mean(SBS2), 6),
    median_SBS2 = round(median(SBS2), 6),
    max_SBS2 = round(max(SBS2), 6),
    pct_high_SBS2 = round(mean(SBS2 > median(patient_data$SBS2)) * 100, 1),
    .groups = "drop"
  )

cat("SBS2 by A3B/A3A quadrant:\n")
print(as.data.frame(quadrant_stats))
cat("\n")

write.table(quadrant_stats, "Diag_Option4_Quadrant_Stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Kruskal-Wallis and pairwise tests ---
cat("Kruskal-Wallis test across quadrants:\n")
kw_test <- kruskal.test(SBS2 ~ Quadrant, data = patient_data)
print(kw_test)

cat("\nPairwise Wilcoxon tests:\n")
pw_test <- pairwise.wilcox.test(patient_data$SBS2, patient_data$Quadrant,
                                 p.adjust.method = "BH")
print(pw_test)
cat("\n")

# --- Option 4 Plots ---

# Plot 4A: Box/violin plot of SBS2 by quadrant
quadrant_colors <- c(
  "A3B Absent / A3A Low"  = "#9bc1bc",
  "A3B Absent / A3A High" = "#5e81ac",
  "A3B Present / A3A Low" = "#f4f1bb",
  "A3B Present / A3A High" = "#ed6a5a"
)

opt4_boxplot <- ggplot(patient_data, aes(x = Quadrant, y = SBS2, fill = Quadrant)) +
  geom_violin(alpha = 0.5, color = "grey30", linewidth = 0.5) +
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = 21,
               outlier.size = 2, color = "black") +
  geom_jitter(shape = 21, size = 1.5, alpha = 0.3, color = "black",
              width = 0.1, stroke = 0.3) +
  scale_fill_manual(values = quadrant_colors, guide = "none") +
  labs(x = NULL,
       y = "SBS2 Weight",
       title = "SBS2 Distribution Across A3B/A3A Quadrants",
       subtitle = paste0("Split at median A3B = ", round(median_a3b, 2),
                         " and median A3A = ", round(median_a3a, 2))) +
  theme_panel +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12))

# Plot 4B: A3A vs SBS2 scatter, colored by A3B gate status
opt4_scatter <- ggplot(patient_data, aes(x = APOBEC3A, y = SBS2)) +
  geom_point(aes(fill = A3B_gate), shape = 21, size = 3.5,
             alpha = 0.6, color = "black", stroke = 0.5) +
  geom_smooth(aes(color = A3B_gate, linetype = A3B_gate),
              method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
  scale_fill_manual(values = c("A3B Absent" = "#9bc1bc",
                                "A3B Present" = "#ed6a5a"),
                    name = "A3B Status") +
  scale_color_manual(values = c("A3B Absent" = "#5a9b96",
                                 "A3B Present" = "#c0392b"),
                     name = "A3B Status") +
  scale_linetype_manual(values = c("A3B Absent" = "dashed",
                                    "A3B Present" = "solid"),
                        name = "A3B Status") +
  labs(x = "APOBEC3A Expression (FPKM-UQ)",
       y = "SBS2 Weight",
       title = "A3A as Continuous Driver — Gated by A3B Presence",
       subtitle = "A3A should only correlate with SBS2 when A3B is present") +
  theme_panel +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "grey70"))

# Plot 4C: Heatmap-style summary of the 2x2 model
heatmap_data <- quadrant_stats %>%
  mutate(label = paste0("n=", n, "\nmedian=", median_SBS2))

opt4_heatmap <- ggplot(heatmap_data, aes(x = A3B_gate, y = A3A_level)) +
  geom_tile(aes(fill = median_SBS2), color = "black", linewidth = 1) +
  geom_text(aes(label = label), size = 6, fontface = "bold") +
  scale_fill_viridis(option = "viridis", name = "Median\nSBS2") +
  labs(x = NULL, y = NULL,
       title = "2×2 Model: A3B Gate × A3A Driver",
       subtitle = "Median SBS2 weight per quadrant") +
  theme_panel +
  theme(axis.text = element_text(size = 16, face = "bold"),
        panel.grid = element_blank())

opt4_combined <- (opt4_boxplot) /
                 (opt4_scatter | opt4_heatmap) +
  plot_annotation(title = "OPTION 4: A3B as Binary Gate / A3A as Continuous Driver",
                  theme = theme(plot.title = element_text(size = 22, face = "bold",
                                                          hjust = 0.5)))

ggsave("Diag_Option4_Gate_Driver_Model.pdf", opt4_combined,
       width = 18, height = 16)
ggsave("Diag_Option4_Gate_Driver_Model.png", opt4_combined,
       width = 18, height = 16, dpi = 200)
cat("Option 4 figures saved.\n\n")


# ============================================================
# FINAL SUMMARY
# ============================================================
cat("============================================================\n")
cat("ALL DIAGNOSTICS COMPLETE\n")
cat("============================================================\n\n")

cat("Output files:\n")
cat("  Option 1: Diag_Option1_Threshold_Correlation.pdf/.png\n")
cat("            Diag_Option1_A3B_Threshold_Correlations.tsv\n")
cat("  Option 2: Diag_Option2_Log_Transform.pdf/.png\n")
cat("            Diag_Option2_Linear_vs_Log_Coefficients.tsv\n")
cat("  Option 3: Diag_Option3_Sequential_ANOVA.pdf/.png\n")
cat("            Diag_Option3_Sequential_ANOVA.tsv\n")
cat("  Option 4: Diag_Option4_Gate_Driver_Model.pdf/.png\n")
cat("            Diag_Option4_Quadrant_Stats.tsv\n")
cat("============================================================\n")
