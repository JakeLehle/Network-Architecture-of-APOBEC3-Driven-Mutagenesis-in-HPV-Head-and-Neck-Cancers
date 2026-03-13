library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(pROC)

setwd("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1")

cat("============================================================\n")
cat("Diagnostic: Are A3C and A3H Independent SBS2 Contributors\n")
cat("          or Bystanders Along for the Ride?\n")
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
median_sbs2 <- median(patient_data$SBS2)

a3_genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3H")
a3_labels <- c("A3A", "A3B", "A3C", "A3H")

# Define SBS2 groups
patient_data <- patient_data %>%
  mutate(SBS2_group = ifelse(SBS2 >= median_sbs2, "SBS2 High", "SBS2 Low"))

patient_data$SBS2_group <- factor(patient_data$SBS2_group,
                                   levels = c("SBS2 Low", "SBS2 High"))

n_high <- sum(patient_data$SBS2_group == "SBS2 High")
n_low  <- sum(patient_data$SBS2_group == "SBS2 Low")

cat("HNSC tumor patients:", n_tumors, "\n")
cat("Median SBS2:", round(median_sbs2, 6), "\n")
cat("SBS2 High (>= median):", n_high, "\n")
cat("SBS2 Low (< median):", n_low, "\n\n")

# Shared theme
theme_panel <- theme_bw() +
  theme(
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5, color = "grey30"),
    strip.text = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  )


# ############################################################
#
#  APPROACH A: "Who's Always at the Scene of the Crime?"
#  Fraction of high-SBS2 tumors with elevated expression
#  of each A3 family member
#
# ############################################################
cat("============================================================\n")
cat("APPROACH A: Presence at the Scene of the Crime\n")
cat("============================================================\n\n")

# For each A3, what fraction of high-SBS2 tumors have above-median
# expression of that A3? If A3B is necessary, it should be near 100%.
# If A3C is a bystander, it should be near 50% (coin flip).

scene_results <- data.frame(
  Gene = character(),
  Pct_elevated_in_SBS2_High = numeric(),
  Pct_elevated_in_SBS2_Low = numeric(),
  Enrichment_ratio = numeric(),
  Fisher_pvalue = numeric(),
  stringsAsFactors = FALSE
)

for(i in seq_along(a3_genes)) {
  gene <- a3_genes[i]
  label <- a3_labels[i]
  gene_median <- median(patient_data[[gene]])

  # Count elevated in each SBS2 group
  high_sbs2 <- patient_data %>% filter(SBS2_group == "SBS2 High")
  low_sbs2  <- patient_data %>% filter(SBS2_group == "SBS2 Low")

  n_elevated_high <- sum(high_sbs2[[gene]] >= gene_median)
  n_elevated_low  <- sum(low_sbs2[[gene]] >= gene_median)

  pct_high <- n_elevated_high / nrow(high_sbs2) * 100
  pct_low  <- n_elevated_low / nrow(low_sbs2) * 100

  # Fisher's exact test: is the A3 more likely to be elevated in
  # high-SBS2 tumors than low-SBS2 tumors?
  contingency <- matrix(c(
    n_elevated_high,
    nrow(high_sbs2) - n_elevated_high,
    n_elevated_low,
    nrow(low_sbs2) - n_elevated_low
  ), nrow = 2, byrow = TRUE)

  fisher <- fisher.test(contingency)

  enrichment <- (pct_high / 50) / (pct_low / 50)  # ratio vs expected 50%

  scene_results <- rbind(scene_results, data.frame(
    Gene = label,
    Pct_elevated_in_SBS2_High = round(pct_high, 1),
    Pct_elevated_in_SBS2_Low = round(pct_low, 1),
    Enrichment_ratio = round(enrichment, 2),
    Fisher_pvalue = fisher$p.value,
    stringsAsFactors = FALSE
  ))

  cat(label, ":\n")
  cat("  Elevated in SBS2-High:", round(pct_high, 1), "%\n")
  cat("  Elevated in SBS2-Low:", round(pct_low, 1), "%\n")
  cat("  Enrichment ratio:", round(enrichment, 2), "\n")
  cat("  Fisher p:", format(fisher$p.value, scientific = TRUE, digits = 3), "\n\n")
}

write.table(scene_results, "Diag_ApproachA_Scene_Presence.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Approach A Plot ---
# Grouped bar chart: % of high-SBS2 vs low-SBS2 tumors with elevated A3

scene_plot_data <- data.frame(
  Gene = rep(a3_labels, 2),
  SBS2_Group = rep(c("SBS2 High", "SBS2 Low"), each = 4),
  Pct_Elevated = c(scene_results$Pct_elevated_in_SBS2_High,
                   scene_results$Pct_elevated_in_SBS2_Low),
  stringsAsFactors = FALSE
)
scene_plot_data$Gene <- factor(scene_plot_data$Gene,
                                levels = c("A3B", "A3A", "A3C", "A3H"))
scene_plot_data$SBS2_Group <- factor(scene_plot_data$SBS2_Group,
                                      levels = c("SBS2 High", "SBS2 Low"))

# Add significance markers
sig_markers <- scene_results %>%
  mutate(sig = case_when(
    Fisher_pvalue < 0.001 ~ "***",
    Fisher_pvalue < 0.01  ~ "**",
    Fisher_pvalue < 0.05  ~ "*",
    TRUE                  ~ "n.s."
  ))

approach_A <- ggplot(scene_plot_data,
                     aes(x = Gene, y = Pct_Elevated, fill = SBS2_Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  annotate("text", x = 4.5, y = 52, label = "50% (random)",
           hjust = 0, size = 4, color = "grey40", fontface = "italic") +
  # Add significance stars above each gene pair
  geom_text(data = sig_markers,
            aes(x = Gene, y = pmax(Pct_elevated_in_SBS2_High,
                                    Pct_elevated_in_SBS2_Low) + 3,
                label = sig),
            inherit.aes = FALSE, size = 10, vjust = 0) +
  scale_fill_manual(values = c("SBS2 HIGH" = "#ed6a5a",
                                "SBS2 LOW" = "#9bc1bc"),
                    name = "SBS2 Group") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25)) +
  labs(x = NULL,
       y = "% of Tumors with Elevated Expression\n(above median)",
       title = "Approach A: Who Is Present When SBS2 Is High?",
       subtitle = paste0("If a gene is necessary for SBS2, it should be elevated ",
                         "in nearly all SBS2-High tumors")) +
  theme_panel +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"))


# ############################################################
#
#  APPROACH B: ROC/AUC — Predictive Power of Each A3
#
# ############################################################
cat("============================================================\n")
cat("APPROACH B: ROC/AUC — Predictive Power for SBS2 Status\n")
cat("============================================================\n\n")

# For each A3, use its expression as a single predictor to classify
# tumors as SBS2-High vs SBS2-Low. The AUC tells us how well knowing
# that A3's expression separates the two groups.

patient_data$SBS2_binary <- as.numeric(patient_data$SBS2_group == "SBS2 High")

roc_list <- list()
auc_results <- data.frame(
  Gene = character(),
  AUC = numeric(),
  AUC_CI_lower = numeric(),
  AUC_CI_upper = numeric(),
  stringsAsFactors = FALSE
)

for(i in seq_along(a3_genes)) {
  gene <- a3_genes[i]
  label <- a3_labels[i]

  roc_obj <- roc(patient_data$SBS2_binary, patient_data[[gene]],
                 quiet = TRUE, direction = "<")
  ci_obj <- ci.auc(roc_obj, conf.level = 0.95)

  roc_list[[label]] <- roc_obj

  auc_results <- rbind(auc_results, data.frame(
    Gene = label,
    AUC = round(auc(roc_obj), 4),
    AUC_CI_lower = round(ci_obj[1], 4),
    AUC_CI_upper = round(ci_obj[3], 4),
    stringsAsFactors = FALSE
  ))

  cat(label, ": AUC =", round(auc(roc_obj), 4),
      " (95% CI:", round(ci_obj[1], 4), "-", round(ci_obj[3], 4), ")\n")
}

# DeLong tests: compare each A3 AUC to A3B's AUC
cat("\nDeLong test comparing each A3's AUC to A3B:\n")
for(label in c("A3A", "A3C", "A3H")) {
  test <- roc.test(roc_list[["A3B"]], roc_list[[label]], method = "delong")
  cat("  A3B vs", label, ": p =",
      format(test$p.value, scientific = TRUE, digits = 3), "\n")
}
cat("\n")

write.table(auc_results, "Diag_ApproachB_AUC_Results.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Approach B Plot: Overlaid ROC curves ---

# Build ROC curve data for ggplot
roc_plot_data <- data.frame()

for(label in a3_labels) {
  roc_obj <- roc_list[[label]]
  roc_df <- data.frame(
    Sensitivity = roc_obj$sensitivities,
    Specificity = 1 - roc_obj$specificities,  # FPR = 1 - specificity
    Gene = label,
    stringsAsFactors = FALSE
  )
  roc_plot_data <- rbind(roc_plot_data, roc_df)
}

roc_plot_data$Gene <- factor(roc_plot_data$Gene,
                              levels = c("A3B", "A3A", "A3C", "A3H"))

# Legend labels with AUC values
auc_labels <- setNames(
  paste0(auc_results$Gene, " (AUC = ", round(auc_results$AUC, 3), ")"),
  auc_results$Gene
)

roc_colors <- c("A3B" = "#bf616a", "A3A" = "#5e81ac",
                "A3C" = "#a3be8c", "A3H" = "#d08770")

approach_B <- ggplot(roc_plot_data,
                     aes(x = Specificity, y = Sensitivity,
                         color = Gene, linewidth = Gene)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.5) +
  geom_line() +
  scale_color_manual(values = roc_colors, labels = auc_labels,
                     name = "A3 Family Member") +
  scale_linewidth_manual(values = c("A3B" = 1.5, "A3A" = 1.5,
                                     "A3C" = 1, "A3H" = 1),
                          guide = "none") +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       title = "Approach B: Which A3 Best Predicts SBS2 Status?",
       subtitle = paste0("ROC curves for classifying tumors as SBS2-High vs SBS2-Low  |  ",
                         "Diagonal = random (AUC = 0.5)")) +
  theme_panel +
  theme(legend.position = c(0.65, 0.25),
        legend.background = element_rect(fill = "white", color = "grey70"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        aspect.ratio = 1)

# --- Approach B supplement: AUC bar chart for quick comparison ---
auc_results$Gene <- factor(auc_results$Gene,
                            levels = c("A3B", "A3A", "A3C", "A3H"))

approach_B_bar <- ggplot(auc_results, aes(x = Gene, y = AUC, fill = Gene)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.4, width = 0.6) +
  geom_errorbar(aes(ymin = AUC_CI_lower, ymax = AUC_CI_upper),
                width = 0.2, linewidth = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  annotate("text", x = 4.5, y = 0.52, label = "Random",
           hjust = 0, size = 4, color = "grey40", fontface = "italic") +
  geom_text(aes(label = round(AUC, 3)), vjust = -0.8, size = 5,
            fontface = "bold") +
  scale_fill_manual(values = roc_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL,
       y = "AUC (Area Under ROC Curve)",
       title = "Predictive Power: AUC Summary",
       subtitle = "95% confidence intervals shown") +
  theme_panel


# ############################################################
#
#  APPROACH C: Expression Shift Between SBS2 Groups
#
# ############################################################
cat("============================================================\n")
cat("APPROACH C: Expression Shift Between SBS2 Groups\n")
cat("============================================================\n\n")

# For each A3, compare expression distributions between SBS2-High
# and SBS2-Low groups. A3C and A3H should show minimal shift if
# they are bystanders.

shift_results <- data.frame(
  Gene = character(),
  Median_SBS2_High = numeric(),
  Median_SBS2_Low = numeric(),
  Fold_change = numeric(),
  Wilcoxon_pvalue = numeric(),
  Effect_size_r = numeric(),
  stringsAsFactors = FALSE
)

for(i in seq_along(a3_genes)) {
  gene <- a3_genes[i]
  label <- a3_labels[i]

  high_expr <- patient_data %>%
    filter(SBS2_group == "SBS2 High") %>%
    pull(!!sym(gene))
  low_expr <- patient_data %>%
    filter(SBS2_group == "SBS2 Low") %>%
    pull(!!sym(gene))

  wtest <- wilcox.test(high_expr, low_expr)
  # Effect size r = Z / sqrt(N)
  z_val <- qnorm(wtest$p.value / 2)
  effect_r <- abs(z_val) / sqrt(n_tumors)

  shift_results <- rbind(shift_results, data.frame(
    Gene = label,
    Median_SBS2_High = round(median(high_expr), 4),
    Median_SBS2_Low = round(median(low_expr), 4),
    Fold_change = round(median(high_expr) / max(median(low_expr), 0.001), 2),
    Wilcoxon_pvalue = wtest$p.value,
    Effect_size_r = round(effect_r, 4),
    stringsAsFactors = FALSE
  ))

  cat(label, ":\n")
  cat("  Median in SBS2-High:", round(median(high_expr), 4), "\n")
  cat("  Median in SBS2-Low:", round(median(low_expr), 4), "\n")
  cat("  Fold change:", round(median(high_expr) / max(median(low_expr), 0.001), 2), "\n")
  cat("  Wilcoxon p:", format(wtest$p.value, scientific = TRUE, digits = 3), "\n")
  cat("  Effect size r:", round(effect_r, 4), "\n\n")
}

write.table(shift_results, "Diag_ApproachC_Expression_Shift.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Approach C Plot: Paired boxplots for each A3 ---

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

plot_data_long$A3_label <- factor(plot_data_long$A3_label,
                                   levels = c("A3B", "A3A", "A3C", "A3H"))

# Annotation with significance and effect size
annot_C <- shift_results %>%
  mutate(
    A3_label = factor(Gene, levels = c("A3B", "A3A", "A3C", "A3H")),
    sig = case_when(
      Wilcoxon_pvalue < 0.001 ~ "***",
      Wilcoxon_pvalue < 0.01  ~ "**",
      Wilcoxon_pvalue < 0.05  ~ "*",
      TRUE                    ~ "n.s."
    ),
    label = paste0(sig, "\nr = ", Effect_size_r)
  )

# Calculate y position for annotations (above the tallest whisker)
y_positions <- plot_data_long %>%
  group_by(A3_label) %>%
  summarise(y_max = quantile(Expression, 0.99), .groups = "drop")

annot_C <- annot_C %>%
  left_join(y_positions, by = "A3_label")

approach_C <- ggplot(plot_data_long,
                     aes(x = SBS2_group, y = Expression, fill = SBS2_group)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 21, outlier.size = 1.5,
               outlier.alpha = 0.3, color = "black", linewidth = 0.5,
               width = 0.6) +
  geom_jitter(shape = 21, size = 1, alpha = 0.15, color = "black",
              fill = "grey40", width = 0.12, stroke = 0.2) +
  facet_wrap(~ A3_label, scales = "free_y", nrow = 1) +
  # Add significance + effect size annotation
  geom_text(data = annot_C,
            aes(x = 1.5, y = y_max * 1.1, label = label),
            inherit.aes = FALSE, size = 5, fontface = "bold",
            color = "grey20", lineheight = 0.9) +
  scale_fill_manual(values = c("SBS2 High" = "#ed6a5a",
                                "SBS2 Low" = "#9bc1bc"),
                    name = "SBS2 Group") +
  labs(x = NULL,
       y = "Expression (FPKM-UQ)",
       title = "Approach C: Does A3 Expression Shift Between SBS2 Groups?",
       subtitle = paste0("Bystander genes should show minimal shift  |  ",
                         "Effect size r: small < 0.1, medium 0.1–0.3, large > 0.3")) +
  theme_panel +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"))

# --- Approach C supplement: Effect size bar chart ---
shift_results$Gene <- factor(shift_results$Gene,
                              levels = c("A3B", "A3A", "A3C", "A3H"))

# Color by effect size magnitude
shift_results <- shift_results %>%
  mutate(Effect_category = case_when(
    Effect_size_r >= 0.5  ~ "Medium",
    Effect_size_r >= 0.2  ~ "Small",
    TRUE                  ~ "None/Ignored"
  ))

approach_C_effect <- ggplot(shift_results,
                            aes(x = Gene, y = Effect_size_r,
                                fill = Effect_category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.4, width = 0.6) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey50") +
  annotate("text", x = 4.6, y = 0.105, label = "Small/Medium",
           hjust = 0, size = 3.5, color = "grey40", fontface = "italic") +
  annotate("text", x = 4.6, y = 0.305, label = "Medium/Large",
           hjust = 0, size = 3.5, color = "grey40", fontface = "italic") +
  geom_text(aes(label = round(Effect_size_r, 3)), vjust = -0.8,
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Large" = "#ed6a5a",
                                "Medium" = "#f4f1bb",
                                "Small/None" = "#9bc1bc"),
                    name = "Effect Size") +
  labs(x = NULL,
       y = "Effect Size (r)",
       title = "Effect Size: Expression Shift Between SBS2 Groups",
       subtitle = "Wilcoxon rank-sum effect size r = |Z| / sqrt(N)") +
  theme_panel +
  theme(legend.position = "bottom")


# ############################################################
#
#  BONUS: Correlation between A3C/A3H and A3B
#  (Are they co-expressed or independent?)
#
# ############################################################
cat("============================================================\n")
cat("BONUS: A3 Co-Expression Analysis\n")
cat("============================================================\n\n")

# If A3C and A3H are bystanders driven by immune infiltration,
# they should correlate with each other (shared immune source)
# but not necessarily with A3B (tumor-intrinsic).
# A3A and A3B may co-correlate (both tumor-intrinsic, linked regulation).

cat("Pairwise Spearman correlations between A3 family members:\n\n")

coexpr_matrix <- matrix(NA, nrow = 4, ncol = 4)
rownames(coexpr_matrix) <- a3_labels
colnames(coexpr_matrix) <- a3_labels

for(i in 1:4) {
  for(j in 1:4) {
    coexpr_matrix[i, j] <- cor(patient_data[[a3_genes[i]]],
                                patient_data[[a3_genes[j]]],
                                method = "spearman")
  }
}

print(round(coexpr_matrix, 3))
cat("\n")

# Reshape for heatmap
coexpr_df <- expand.grid(Gene1 = a3_labels, Gene2 = a3_labels,
                          stringsAsFactors = FALSE)
coexpr_df$Rho <- as.vector(coexpr_matrix)
coexpr_df$Gene1 <- factor(coexpr_df$Gene1, levels = a3_labels)
coexpr_df$Gene2 <- factor(coexpr_df$Gene2, levels = rev(a3_labels))

bonus_heatmap <- ggplot(coexpr_df, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = Rho), color = "black", linewidth = 1) +
  geom_text(aes(label = round(Rho, 2)), size = 7, fontface = "bold") +
  scale_fill_gradient2(low = "#5e81ac", mid = "white", high = "#bf616a",
                       midpoint = 0, limits = c(-0.5, 1),
                       name = "Spearman\nrho") +
  labs(x = NULL, y = NULL,
       title = "A3 Family Co-Expression",
       subtitle = "Do A3C/A3H cluster together (immune) vs A3A/A3B (tumor)?") +
  theme_panel +
  theme(axis.text = element_text(size = 16, face = "bold"),
        panel.grid = element_blank())


# ############################################################
#
#  ASSEMBLE AND SAVE
#
# ############################################################
cat("============================================================\n")
cat("Assembling composite figures...\n")
cat("============================================================\n\n")

# --- Composite 1: The three main approaches ---
# Row 1: Approach A (bar chart) + Approach B (ROC curves)
# Row 2: Approach C (paired boxplots, full width)

composite_main <- (approach_A | approach_B) /
                  approach_C +
  plot_annotation(
    title = "Are A3C and A3H Independent Contributors to SBS2 or Bystanders?",
    subtitle = paste0("TCGA-HNSC Tumors (n = ", n_tumors, ")  |  ",
                      "SBS2 split at median = ", round(median_sbs2, 4)),
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "grey30")
    )
  ) +
  plot_layout(heights = c(1, 0.8))

ggsave("Diag_Bridge_A3C_A3H_Bystander_Main.pdf", composite_main,
       width = 22, height = 18, units = "in")
ggsave("Diag_Bridge_A3C_A3H_Bystander_Main.png", composite_main,
       width = 22, height = 18, units = "in", dpi = 200)

# --- Composite 2: Supplementary panels ---
composite_supp <- (approach_B_bar | approach_C_effect | bonus_heatmap) +
  plot_annotation(
    title = "Supplementary: Quantitative Summary",
    theme = theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    )
  )

ggsave("Diag_Bridge_A3C_A3H_Bystander_Supp.pdf", composite_supp,
       width = 22, height = 8, units = "in")
ggsave("Diag_Bridge_A3C_A3H_Bystander_Supp.png", composite_supp,
       width = 22, height = 8, units = "in", dpi = 200)

# --- Individual panels ---
ggsave("Diag_ApproachA_Scene_Presence.pdf", approach_A, width = 10, height = 8)
ggsave("Diag_ApproachB_ROC_Curves.pdf", approach_B, width = 10, height = 9)
ggsave("Diag_ApproachB_AUC_Bar.pdf", approach_B_bar, width = 8, height = 7)
ggsave("Diag_ApproachC_Expression_Shift.pdf", approach_C, width = 18, height = 7)
ggsave("Diag_ApproachC_Effect_Size.pdf", approach_C_effect, width = 8, height = 7)
ggsave("Diag_Bonus_CoExpression_Heatmap.pdf", bonus_heatmap, width = 8, height = 7)

cat("All figures saved.\n\n")


# ============================================================
# SUMMARY
# ============================================================
cat("============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("============================================================\n\n")

cat("APPROACH A — Scene of the Crime:\n")
for(i in 1:nrow(scene_results)) {
  cat("  ", scene_results$Gene[i], ": ",
      scene_results$Pct_elevated_in_SBS2_High[i],
      "% elevated in SBS2-High vs ",
      scene_results$Pct_elevated_in_SBS2_Low[i],
      "% in SBS2-Low\n")
}

cat("\nAPPROACH B — ROC/AUC:\n")
for(i in 1:nrow(auc_results)) {
  cat("  ", as.character(auc_results$Gene[i]), ": AUC = ", auc_results$AUC[i],
      " (", auc_results$AUC_CI_lower[i], "-", auc_results$AUC_CI_upper[i], ")\n")
}

cat("\nAPPROACH C — Expression Shift:\n")
for(i in 1:nrow(shift_results)) {
  cat("  ", as.character(shift_results$Gene[i]),
      ": effect size r =", shift_results$Effect_size_r[i],
      " (", shift_results$Effect_category[i], ")\n")
}

cat("\nBONUS — Co-expression matrix:\n")
print(round(coexpr_matrix, 3))

cat("\nINTERPRETATION GUIDE:\n")
cat("  If A3C and A3H are bystanders, you should see:\n")
cat("  - Approach A: ~50% presence in both SBS2 groups (no enrichment)\n")
cat("  - Approach B: AUC near 0.5 (no predictive power)\n")
cat("  - Approach C: Small effect size, no meaningful expression shift\n")
cat("  - Bonus: A3C-A3H correlating with each other but not with A3B\n")
cat("\n  If A3A and A3B are true drivers, you should see:\n")
cat("  - Approach A: Strong enrichment in SBS2-High tumors\n")
cat("  - Approach B: AUC substantially above 0.5\n")
cat("  - Approach C: Large effect size, clear expression shift\n")
cat("  - Bonus: A3A-A3B may correlate (shared tumor regulation)\n")

cat("\nOUTPUT FILES:\n")
cat("  Main:  Diag_Bridge_A3C_A3H_Bystander_Main.pdf/.png\n")
cat("  Supp:  Diag_Bridge_A3C_A3H_Bystander_Supp.pdf/.png\n")
cat("  Individual panels: Diag_ApproachA/B/C PDFs\n")
cat("  Data:  Diag_ApproachA_Scene_Presence.tsv\n")
cat("         Diag_ApproachB_AUC_Results.tsv\n")
cat("         Diag_ApproachC_Expression_Shift.tsv\n")
cat("============================================================\n")
