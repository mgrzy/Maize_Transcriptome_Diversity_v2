#!/usr/bin/env Rscript
# Compare held-out prediction performance between SuSiE and OSCA across the
# number of variants in the SuSiE credible set.
#
# Delta R² is calculated as: SuSiE test R² - OSCA test R².
#
# Run from the project root:
#   Rscript scripts/local/figures/Fig_SX_Delta_R2_by_SuSiE_credible_set_size.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

config <- list(
  input = "results/summary/Heritability_and_r2/narrow_sense_h2_and_r2.tsv.gz",
  output = "results/figures/supp/Fig_Delta_R2_by_SuSiE_credible_set_size.png"
)

if (!file.exists(config$input)) stop("Input file does not exist: ", config$input)
dir.create(dirname(config$output), recursive = TRUE, showWarnings = FALSE)

required_columns <- c("susie_n_variants", "susie_r2_test", "osca_r2_test")
results <- fread(config$input, select = required_columns)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns)) {
  stop("Input is missing required column(s): ", paste(missing_columns, collapse = ", "))
}

plot_data <- results[
  is.finite(susie_n_variants) & is.finite(susie_r2_test) & is.finite(osca_r2_test),
  .(susie_n_variants, delta_r2 = susie_r2_test - osca_r2_test)
]
if (!nrow(plot_data)) stop("No rows with finite held-out R² values were found.")
if (any(plot_data$susie_n_variants < 0L)) stop("susie_n_variants contains negative values.")

group_counts <- plot_data[, .N, by = susie_n_variants][order(susie_n_variants)]
plot_data[, credible_set_size := factor(susie_n_variants, levels = group_counts$susie_n_variants)]
x_labels <- setNames(
  paste0(group_counts$susie_n_variants, "\n(n = ", format(group_counts$N, big.mark = ","), ")"),
  group_counts$susie_n_variants
)
# A small number of very negative test R² values would otherwise obscure the
# distribution around zero. The boxplot statistics still use every gene.
display_limits <- unname(stats::quantile(plot_data$delta_r2, probs = c(0.005, 0.995)))

figure <- ggplot(plot_data, aes(x = credible_set_size, y = delta_r2, fill = credible_set_size)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
  geom_boxplot(colour = "black", width = 0.7,
               outlier.shape = NA, linewidth = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 1.8,
               fill = "white", colour = "black") +
  scale_x_discrete(labels = x_labels) +
  scale_fill_viridis_d(name = "N causal sets", guide = "none") +
  labs(
    x = "Number of variants in SuSiE credible set\n(number of genes below each category)",
    y = expression(Delta * R^2 * " (SuSiE test " * R^2 * " - OSCA test " * R^2 * ")"),
    caption = "Display range is the central 99% of ΔR² values; boxplot statistics use all genes."
  ) +
  coord_cartesian(ylim = display_limits) +
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(lineheight = 0.9),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

ggsave(config$output, figure, width = 18, height = 12, units = "cm", dpi = 300, bg = "white")
message("Wrote figure to ", config$output)
