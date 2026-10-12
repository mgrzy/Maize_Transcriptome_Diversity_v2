#!/usr/bin/env Rscript
# Run from the repository root:
# Rscript --vanilla scripts/local/figures/Fig2E_S8CD_presentation.R
# Preserve the source panels' separate gene filters; share styling and colours.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

heritability <- fread("results/summary/Heritability_and_r2/narrow_sense_h2_and_r2.tsv.gz")
leads <- fread("results/summary/AF/aFCn_MAF_se_filtered.tsv.gz", select = "gene_id")
expression_summary <- fread("results/full/GeneExpressionVariance/edgeR_log2_cpm_filtered_mean_sd.csv.gz")
if ("Gene" %in% names(expression_summary)) setnames(expression_summary, "Gene", "gene")
if (anyDuplicated(heritability$gene) || anyDuplicated(expression_summary$gene)) {
  stop("Input gene identifiers must be unique.")
}

# Retain Figure 2E's original positive-CS genes and add all zero-CS genes
# with finite adjusted heritability. Zero-CS genes cannot satisfy the
# original requirement for retained fine-mapped leads.
panel_e_data <- heritability[
  (gene %in% leads$gene_id & osca_r2_test > 0 & susie_r2_test > 0) |
    (susie_n_variants == 0 & is.finite(h2_pc_adjusted))
]
# S8 C/D use complete cases in the expression-summary/heritability join.
expression_data <- merge(
  heritability[, .(gene, h2_pc_adjusted, susie_n_variants)],
  expression_summary[, .(gene, expression_mean_log2_cpm = Mean,
                         expression_sd_log2_cpm = SD)],
  by = "gene", all = FALSE, sort = FALSE
)
expression_data <- expression_data[
  is.finite(h2_pc_adjusted) & is.finite(susie_n_variants) &
    is.finite(expression_mean_log2_cpm) & is.finite(expression_sd_log2_cpm)
]

# Match Figure 2E's viridis colours for positive CS counts. S8 also includes
# zero-CS genes, which receive grey without shifting the positive-count colours.
main_levels <- sort(unique(panel_e_data[susie_n_variants > 0, susie_n_variants]))
group_colours <- setNames(
  scales::viridis_pal()(length(main_levels)), as.character(main_levels)
)
extra_levels <- setdiff(sort(unique(c(panel_e_data$susie_n_variants,
                                     expression_data$susie_n_variants))), main_levels)
if (length(extra_levels)) {
  group_colours <- c(group_colours, setNames(rep("#BDBDBD", length(extra_levels)),
                                            as.character(extra_levels)))
}
group_levels <- sort(unique(c(main_levels, extra_levels)))
theme_set(theme_classic(base_size = 10))
theme_update(axis.title = element_text(size = rel(1.15)),
             axis.text = element_text(colour = "black"))

make_panel <- function(panel_data, y_column, y_label, expression_panel = FALSE) {
  panel_data <- copy(panel_data)
  panel_data[, cs_group := factor(susie_n_variants, levels = group_levels)]
  plot <- ggplot(panel_data, aes(cs_group, .data[[y_column]], fill = cs_group)) +
    geom_boxplot(width = 0.7, linewidth = 0.35,
                 outlier.shape = if (expression_panel) NA else 19) +
    scale_fill_manual(values = group_colours) +
    labs(x = "Number of credible sets", y = y_label) +
    theme(legend.position = "none")
  if (expression_panel) {
    plot <- plot + stat_summary(fun = median, geom = "point", shape = 21,
                                size = 1.7, fill = "white", colour = "black")
  } else {
    plot <- plot + stat_summary(fun = mean, geom = "point", shape = 21,
                                size = 1.7, fill = "white", colour = "black")
  }
  plot
}

figure <- make_panel(expression_data, "expression_mean_log2_cpm", "Mean expression (log2 CPM)", TRUE) |
  make_panel(expression_data, "expression_sd_log2_cpm", "Expression SD (log2 CPM)", TRUE) |
  make_panel(panel_e_data, "h2_pc_adjusted", "Narrow-sense heritability")
figure_file <- "results/figures/presentation/Fig2E_S8CD.png"
dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, figure, width = 27, height = 7, units = "cm", dpi = 300)
message("Saved presentation figure: ", figure_file)
