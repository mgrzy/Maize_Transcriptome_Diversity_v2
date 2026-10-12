#!/usr/bin/env Rscript
# Run from the repository root:
# Rscript --vanilla scripts/local/figures/Fig2_AB_presentation.R
# Reuse the panel A/B definitions from Fig2.R without loading other panels.
# Presentation figures omit panel letters.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(data.table)
})

source_expressions <- parse("scripts/local/figures/Fig2.R")
required_objects <- c("base_font_size", "myPal", "cs", "cs.var", "g1", "g2")
loaded_objects <- character()
for (expr in source_expressions) {
  if (!is.call(expr)) next
  operation <- as.character(expr[[1L]])
  if (operation == "<-" && is.symbol(expr[[2L]])) {
    object <- as.character(expr[[2L]])
    if (object %in% required_objects) {
      eval(expr)
      loaded_objects <- c(loaded_objects, object)
    }
  } else if (operation %in% c("theme_set", "theme_update")) {
    eval(expr)
  }
}
if (length(setdiff(required_objects, loaded_objects))) {
  stop("Required panel definitions are missing from Fig2.R.")
}

# Weight each credible-set count by the number of genes in that bin,
# including the zero-CS bin shown in the original panel.
median_cs <- median(rep(cs$cs, times = cs$n))
genes_multiple_cs <- sum(cs$n[cs$cs > 1])
annotation_text <- paste0(
  "Median CS per gene: ", median_cs,
  "\nGenes with >1 CS: ", comma(genes_multiple_cs)
)
g1 <- g1 +
  labs(x = "Number of credible causal sets") +
  annotate("text", x = Inf, y = Inf, label = annotation_text,
           hjust = 1.05, vjust = 1.3, size = 3)

figure_file <- "results/figures/presentation/Fig2_AB.png"
dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, plot = g1 | g2,
       width = 18, height = 6, units = "cm", dpi = 300)
message("Saved presentation figure: ", figure_file)
message(annotation_text)
