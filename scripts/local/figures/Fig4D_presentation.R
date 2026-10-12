#!/usr/bin/env Rscript
# Run from the repository root:
# Rscript --vanilla scripts/local/figures/Fig4D_presentation.R
# Reuse Figure 4D's deciles, gene filters, and colours without other panels.
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

required_objects <- c("cfg", "comparison_labels", "myPal", "comparison_colours",
                      "fm", "deg", "plot_data", "panel")
loaded_objects <- character()
for (expr in parse("scripts/local/figures/Fig4.R")) {
  if (!is.call(expr) || !identical(expr[[1L]], as.name("<-")) ||
      !is.symbol(expr[[2L]])) next
  object <- as.character(expr[[2L]])
  if (object %in% required_objects) {
    eval(expr)
    loaded_objects <- c(loaded_objects, object)
  }
  if (object == "panel") break
}
if (length(setdiff(required_objects, loaded_objects)) || !nrow(plot_data)) {
  stop("Required Figure 4D definitions or plotted data are missing.")
}
if (anyDuplicated(deg$gene_id)) stop("DEG summary has duplicate gene IDs.")

# Total DEGs counts all significant genes in the DEG summary, including genes
# without fine-mapped variants. Mean FST uses exactly the plotted finite
# gene-variant rows across all ten deciles, matching the boxplots' weighting.
total_degs <- imap_dfr(comparison_labels, function(label, comparison) {
  n_degs <- n_distinct(deg$gene_id[deg[[paste0(comparison, "_DEG")]] %in% TRUE])
  tibble(comparison = label, total_DEGs = n_degs)
})
annotations <- plot_data %>%
  group_by(comparison) %>%
  summarise(mean_Fst = mean(Fst), .groups = "drop") %>%
  left_join(total_degs, by = "comparison") %>%
  mutate(comparison = factor(comparison, levels = unname(comparison_labels)),
         label = sprintf("Total DEGs: %s", scales::comma(total_DEGs)))

figure <- panel +
  geom_text(data = annotations, aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 2.7) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
  theme(axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10), strip.text = element_text(size = 8))

figure_file <- "results/figures/presentation/Fig4D.png"
dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, figure, width = 30, height = 8, units = "cm", dpi = 300)
write_tsv(annotations %>% select(comparison, total_DEGs, mean_Fst),
          "results/figures/presentation/Fig4D_summary.tsv")
print(annotations %>% select(comparison, total_DEGs, mean_Fst))
message("Saved presentation figure: ", figure_file)
