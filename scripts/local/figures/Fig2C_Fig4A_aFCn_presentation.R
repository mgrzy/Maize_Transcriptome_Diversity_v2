#!/usr/bin/env Rscript
# Run from the repository root:
# Rscript --vanilla scripts/local/figures/Fig2C_Fig4A_aFCn_presentation.R
# Combine the saved aFC-n histogram (Fig2C) and MAF boxplots (Fig4A).
# Preserve source data, histogram bins, reference lines, and MAF bins.
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

load_panel_definitions <- function(script, required_objects) {
  destination <- parent.frame()
  loaded <- character()
  for (source_expr in parse(script)) {
    if (!is.call(source_expr) || !identical(source_expr[[1L]], as.name("<-")) ||
        !is.symbol(source_expr[[2L]])) next
    object <- as.character(source_expr[[2L]])
    if (object %in% required_objects) {
      eval(source_expr, envir = destination)
      loaded <- c(loaded, object)
    }
  }
  if (length(setdiff(required_objects, loaded))) {
    stop("Required panel definitions are missing from ", script)
  }
}

theme_set(theme_classic(base_size = 10))
load_panel_definitions("scripts/local/figures/Fig2.R", c("myPal", "pip.af", "g3"))
load_panel_definitions("scripts/local/figures/Fig4.R", c("cfg", "myPal", "fm", "maf_plot"))
common_theme <- theme_classic(base_size = 10) +
  theme(axis.title = element_text(size = rel(1.15)),
        axis.text = element_text(colour = "black"))
# The histogram contains gene-variant estimates. For the requested variant
# annotation, count each variant once; a variant qualifies if any of its
# finite gene-specific effects lies outside the original +/- 2 lines.
finite_effects <- pip.af[is.finite(log2_aFC)]
n_effect_variants <- uniqueN(finite_effects$variant_id)
n_large_effect_variants <- uniqueN(finite_effects[abs(log2_aFC) > 2, variant_id])
percent_large_effect_variants <- 100 * n_large_effect_variants / n_effect_variants
effect_annotation <- paste0(
  "Unique variants\n|log₂(aFC-n)| > 2:\n",
  scales::comma(n_large_effect_variants),
  sprintf(" (%.1f%%)", percent_large_effect_variants)
)
histogram <- g3 + common_theme +
  annotate("text", x = -Inf, y = Inf, label = effect_annotation,
           hjust = -0.08, vjust = 1.15, size = 2.8)
maf_boxplots <- maf_plot + common_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
figure <- histogram | maf_boxplots
figure_file <- "results/figures/presentation/Fig2C_Fig4A_aFCn.png"
dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, figure, width = 18, height = 7, units = "cm", dpi = 300)
fwrite(data.table(
  abs_log2_aFCn_threshold = 2,
  unique_variants_with_finite_effect = n_effect_variants,
  unique_variants_above_threshold = n_large_effect_variants,
  percent_unique_variants_above_threshold = percent_large_effect_variants,
  gene_variant_estimates_above_threshold = sum(abs(finite_effects$log2_aFC) > 2),
  definition = "Variant has at least one gene-specific abs(log2_aFC) > 2"
), "results/figures/presentation/Fig2C_Fig4A_aFCn_summary.tsv", sep = "\t")
message("Saved presentation figure: ", figure_file)
message(effect_annotation)
