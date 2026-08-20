#!/usr/bin/env Rscript
# Supplementary diagnostic figure relating expression characteristics (edgeR
# log2 CPM) to adjusted narrow-sense heritability and fine-mapped-variant count.
#
# Inputs (overridable with environment variables):
#   HERITABILITY_FILE
#   EXPRESSION_SUMMARY_FILE
# Output (PNG only): FIGURE_OUTPUT

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  heritability_file = Sys.getenv(
    "HERITABILITY_FILE",
    "results/summary/Heritability_and_r2/narrow_sense_h2_and_r2.tsv.gz"
  ),
  expression_summary_file = Sys.getenv(
    "EXPRESSION_SUMMARY_FILE",
    "results/full/GeneExpressionVariance/edgeR_log2_cpm_filtered_mean_sd.csv.gz"
  ),
  output_file = Sys.getenv("FIGURE_OUTPUT", "results/figures/supp/Fig_SX2_Diag.png")
)

for (path in unlist(cfg[1:2])) {
  if (!file.exists(path)) stop_with("Required input does not exist: %s", path)
}
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)

heritability <- fread(cfg$heritability_file)
expression_summary <- fread(cfg$expression_summary_file)
if ("Gene" %in% names(expression_summary)) setnames(expression_summary, "Gene", "gene")

required_heritability <- c("gene", "h2_pc_adjusted", "susie_n_variants")
required_expression <- c("gene", "Mean", "SD")
missing_heritability <- setdiff(required_heritability, names(heritability))
missing_expression <- setdiff(required_expression, names(expression_summary))
if (length(missing_heritability)) {
  stop_with("Heritability input is missing: %s", paste(missing_heritability, collapse = ", "))
}
if (length(missing_expression)) {
  stop_with("Expression-summary input is missing: %s", paste(missing_expression, collapse = ", "))
}
if (anyDuplicated(heritability$gene) || anyDuplicated(expression_summary$gene)) {
  stop_with("Input gene identifiers must be unique.")
}

plot_data <- merge(
  heritability[, .(gene, h2_pc_adjusted, susie_n_variants)],
  expression_summary[, .(gene, expression_mean_log2_cpm = Mean, expression_sd_log2_cpm = SD)],
  by = "gene",
  all = FALSE,
  sort = FALSE
)
plot_data <- plot_data[
  is.finite(h2_pc_adjusted) & is.finite(susie_n_variants) &
    is.finite(expression_mean_log2_cpm) & is.finite(expression_sd_log2_cpm)
]
if (!nrow(plot_data)) stop_with("No complete genes remain after joining the two diagnostic inputs.")
if (any(plot_data$susie_n_variants < 0L)) stop_with("susie_n_variants contains negative values.")

theme_set(theme_classic(base_size = 10))
panel_theme <- theme(
  axis.text = element_text(colour = "black"),
  plot.margin = margin(5.5, 5.5, 5.5, 5.5)
)

make_heritability_panel <- function(y_column, y_label) {
  correlation <- stats::cor.test(
    plot_data$h2_pc_adjusted,
    plot_data[[y_column]],
    method = "spearman",
    exact = FALSE
  )
  correlation_label <- sprintf(
    "Spearman rho = %.2f\nP %s",
    unname(correlation$estimate),
    format.pval(correlation$p.value, digits = 2, eps = 1e-300)
  )
  ggplot(plot_data, aes(x = h2_pc_adjusted, y = .data[[y_column]])) +
    geom_point(size = 0.45, alpha = 0.13, colour = "#0072B2") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8, colour = "#D55E00") +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = correlation_label,
      hjust = 1.05,
      vjust = 1.25,
      size = 3,
      colour = "black"
    ) +
    labs(x = expression("Adjusted " * italic(h)^2), y = y_label) +
    panel_theme
}

make_variant_count_panel <- function(y_column, y_label) {
  variant_levels <- sort(unique(plot_data$susie_n_variants))
  panel_data <- copy(plot_data)
  panel_data[, variant_group := factor(susie_n_variants, levels = variant_levels)]
  group_counts <- panel_data[, .(n_genes = .N), by = variant_group]
  variant_labels <- setNames(
    paste0(group_counts$variant_group, "\n", format(group_counts$n_genes, big.mark = ",")),
    group_counts$variant_group
  )
  ggplot(panel_data, aes(x = variant_group, y = .data[[y_column]])) +
    geom_boxplot(width = 0.7, outlier.shape = NA, fill = "#56B4E9", colour = "black", linewidth = 0.35) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.7, fill = "white", colour = "black") +
    scale_x_discrete(labels = variant_labels) +
    labs(x = "Number of fine-mapped variants\n(number of genes below each category)", y = y_label) +
    panel_theme +
    theme(axis.text.x = element_text(size = 7, lineheight = 0.9))
}

p_h2_mean <- make_heritability_panel("expression_mean_log2_cpm", "Mean expression (log2 CPM)")
p_h2_sd <- make_heritability_panel("expression_sd_log2_cpm", "Expression SD (log2 CPM)")
p_variant_mean <- make_variant_count_panel("expression_mean_log2_cpm", "Mean expression (log2 CPM)")
p_variant_sd <- make_variant_count_panel("expression_sd_log2_cpm", "Expression SD (log2 CPM)")

figure <- (p_h2_mean | p_h2_sd) / (p_variant_mean | p_variant_sd) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

print(figure)
ggsave(cfg$output_file, figure, width = 8.6, height = 6.5, units = "in", dpi = 400)
message("Wrote supplementary diagnostic figure to: ", cfg$output_file)
