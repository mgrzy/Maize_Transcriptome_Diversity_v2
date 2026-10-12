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
histogram <- g3 + common_theme
maf_boxplots <- maf_plot + common_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
figure <- histogram | maf_boxplots
figure_file <- "results/figures/presentation/Fig2C_Fig4A_aFCn.png"
dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, figure, width = 18, height = 7, units = "cm", dpi = 300)
message("Saved presentation figure: ", figure_file)
