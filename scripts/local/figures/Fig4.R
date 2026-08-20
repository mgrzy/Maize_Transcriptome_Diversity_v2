#!/usr/bin/env Rscript
# Figure 4: Fst across DEG-significance deciles for fine-mapped variants linked
# to significantly differentially expressed genes. Each facet is one DEG/Fst
# comparison and uses the corresponding -log10(limma adjusted P value).

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  fine_mapped_file = "results/summary/AF/aFCn_MAF_se_filtered_DEG_Fst.tsv.gz",
  deg_file = "results/summary/DEG/DEG_limma_summary.csv",
  gerp_file = "results/summary/GERP/positive_GERP_mean_by_gene_and_promoter.csv",
  credible_set_file = "results/summary/SusieR/CS_Gene.csv",
  gerp_test_file = "results/summary/GERP/GERP_effect_size_and_credible_set_tests.csv",
  output_file = "results/figures/Fig4.png"
)
for (path in unlist(cfg[1:5])) {
  if (!file.exists(path)) stop_with("Required input does not exist: %s", path)
}
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)

comparison_labels <- c(
  "SS_vs_Tropical" = "SS vs Tropical",
  "NSS_vs_Tropical" = "NSS vs Tropical",
  "Iodent_vs_Tropical" = "Iodent vs Tropical",
  "Public_SS_vs_NSS" = "Public: SS vs NSS",
  "exPVP_SS_vs_NSS" = "exPVP: SS vs NSS"
)
myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")
comparison_colours <- setNames(myPal[seq_along(comparison_labels)], unname(comparison_labels))
gerp_colours <- c("Other genes" = myPal[1], "Top 10% GERP genes" = myPal[2])

fm <- read_tsv(cfg$fine_mapped_file, show_col_types = FALSE)
deg <- read_csv(cfg$deg_file, show_col_types = FALSE) %>%
  rename(gene_id = gene)
gerp <- read_csv(cfg$gerp_file, show_col_types = FALSE)
credible_sets <- read_csv(cfg$credible_set_file, show_col_types = FALSE) %>%
  rename(gene_id = gene, n_credible_sets = cs)
gerp_tests <- read_csv(cfg$gerp_test_file, show_col_types = FALSE)

required_fm_columns <- c("gene_id", "variant_id", paste0(names(comparison_labels), "_Fst"))
required_deg_columns <- c(
  "gene_id",
  paste0(names(comparison_labels), "_DEG"),
  paste0(names(comparison_labels), "_adj.P.Val")
)
missing_fm <- setdiff(required_fm_columns, names(fm))
missing_deg <- setdiff(required_deg_columns, names(deg))
if (length(missing_fm)) {
  stop_with("Fine-mapped DEG/Fst table is missing: %s", paste(missing_fm, collapse = ", "))
}
if (length(missing_deg)) {
  stop_with("DEG summary is missing: %s", paste(missing_deg, collapse = ", "))
}
required_gerp_columns <- c("gene_id", "gene_mean_positive_GERP")
required_cs_columns <- c("gene_id", "n_credible_sets")
if (length(missing_gerp <- setdiff(required_gerp_columns, names(gerp)))) {
  stop_with("GERP summary is missing: %s", paste(missing_gerp, collapse = ", "))
}
if (length(missing_cs <- setdiff(required_cs_columns, names(credible_sets)))) {
  stop_with("Credible-set summary is missing: %s", paste(missing_cs, collapse = ", "))
}
if (!all(c("test", "p_value") %in% names(gerp_tests))) {
  stop_with("GERP test summary must contain test and p_value columns.")
}
effect_p_value <- gerp_tests %>%
  filter(test == "Wilcoxon rank-sum") %>%
  pull(p_value)
cs_p_value <- gerp_tests %>%
  filter(test == "Quasi-Poisson GLM") %>%
  pull(p_value)
if (length(effect_p_value) != 1L || length(cs_p_value) != 1L ||
    !is.finite(effect_p_value) || !is.finite(cs_p_value)) {
  stop_with("GERP test summary must contain one finite P value for each plotted test.")
}
format_p_value <- function(p_value) formatC(p_value, format = "e", digits = 2)

# Top 90th percentile means the highest-scoring 10% of genes. GERP is
# classified at the gene level before joining to variant-level aFC estimates.
gerp_threshold <- quantile(gerp$gene_mean_positive_GERP, probs = 0.90, na.rm = TRUE)
gerp_groups <- gerp %>%
  filter(is.finite(gene_mean_positive_GERP)) %>%
  transmute(
    gene_id,
    gene_mean_positive_GERP,
    gerp_group = if_else(
      gene_mean_positive_GERP >= gerp_threshold,
      "Top 10% GERP genes",
      "Other genes"
    )
  ) %>%
  mutate(gerp_group = factor(gerp_group, levels = names(gerp_colours)))
if (anyDuplicated(gerp_groups$gene_id)) stop_with("GERP summary has duplicate gene IDs.")

# Each aFC row is a fine-mapped gene-variant pair. First calculate deciles once
# per unique DEG gene using -log10(adjusted P), separately for each comparison.
# Those gene-level deciles are then assigned to all fine-mapped variants linked
# to that gene, preventing genes with many fine-mapped variants from changing
# the decile cut points.
plot_data <- imap_dfr(comparison_labels, function(label, comparison) {
  deg_column <- paste0(comparison, "_DEG")
  fst_column <- paste0(comparison, "_Fst")
  adjusted_p_column <- paste0(comparison, "_adj.P.Val")

  gene_significance <- deg %>%
    transmute(
      gene_id,
      is_DEG = .data[[deg_column]],
      adjusted_p = .data[[adjusted_p_column]]
    ) %>%
    filter(is_DEG %in% TRUE, is.finite(adjusted_p), adjusted_p >= 0) %>%
    mutate(
      minus_log10_adjusted_p = -log10(pmax(adjusted_p, .Machine$double.xmin)),
      significance_decile = ntile(minus_log10_adjusted_p, 10L)
    )

  fm %>%
    filter(is.finite(.data[[fst_column]])) %>%
    transmute(
      gene_id,
      variant_id,
      Fst = .data[[fst_column]]
    ) %>%
    inner_join(gene_significance, by = "gene_id") %>%
    mutate(
      comparison = label
    )
})

if (!nrow(plot_data)) stop_with("No fine-mapped variants linked to significant DEGs had finite Fst values.")
if (any(table(plot_data$comparison) == 0L)) stop_with("At least one comparison has no plottable values.")

plot_data <- plot_data %>%
  mutate(
    comparison = factor(comparison, levels = unname(comparison_labels)),
    significance_decile = factor(
      significance_decile,
      levels = 1:10,
      labels = paste0("D", 1:10)
    )
  )

panel <- ggplot(plot_data, aes(x = significance_decile, y = Fst, fill = comparison)) +
  geom_boxplot(
    width = 0.7,
    colour = "#1D4F5A",
    linewidth = 0.25,
    outlier.size = 0.35,
    outlier.alpha = 0.65
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.1,
    stroke = 0.25,
    fill = "white",
    colour = "#1D4F5A"
  ) +
  scale_fill_manual(values = comparison_colours, guide = "none") +
  facet_wrap(~ comparison, nrow = 1) +
  labs(x = expression(-log[10]("adjusted P")~"decile (low to high)"), y = expression(F[ST])) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 9),
    strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.25),
    strip.text = element_text(size = 7),
    panel.spacing = unit(0.65, "lines")
  )

maf_plot <- fm %>%
  mutate(MAF_bin = cut(MAF, breaks = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))) %>%
  filter(!is.na(MAF_bin), is.finite(log2_aFC)) %>%
  ggplot(aes(x = MAF_bin, y = abs(log2_aFC), fill = MAF_bin)) +
  geom_boxplot(width = 0.7, linewidth = 0.25, outlier.size = 0.35, outlier.alpha = 0.65) +
  scale_fill_manual(values = myPal, guide = "none") +
  labs(x = "MAF", y = expression(abs(log[2](aFC)))) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.title = element_text(size = 9)
  )

gerp_effect_data <- fm %>%
  select(gene_id, variant_id, log2_aFC) %>%
  inner_join(gerp_groups, by = "gene_id") %>%
  filter(is.finite(log2_aFC))
if (!nrow(gerp_effect_data)) stop_with("No fine-mapped variants could be assigned to a GERP group.")

gerp_effect_plot <- ggplot(
  gerp_effect_data,
  aes(x = gerp_group, y = abs(log2_aFC), fill = gerp_group)
) +
  geom_boxplot(width = 0.7, linewidth = 0.25, outlier.size = 0.35, outlier.alpha = 0.65) +
  scale_fill_manual(values = gerp_colours, guide = "none") +
  annotate(
    "text", x = 1.5, y = Inf,
    label = paste0("Wilcoxon P = ", format_p_value(effect_p_value)),
    vjust = 1.5, size = 2.6
  ) +
  labs(x = NULL, y = expression(abs(log[2](aFC)))) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.margin = margin(t = 8, r = 5.5, b = 5.5, l = 5.5)
  )

credible_set_data <- credible_sets %>%
  inner_join(gerp_groups %>% select(gene_id, gerp_group), by = "gene_id") %>%
  filter(is.finite(n_credible_sets), n_credible_sets >= 0)
if (!nrow(credible_set_data)) stop_with("No credible-set counts could be assigned to a GERP group.")
credible_set_distribution <- credible_set_data %>%
  count(gerp_group, n_credible_sets, name = "n_genes") %>%
  group_by(gerp_group) %>%
  mutate(proportion_genes = n_genes / sum(n_genes)) %>%
  ungroup()

credible_set_plot <- ggplot(
  credible_set_distribution,
  aes(x = factor(n_credible_sets), y = proportion_genes, fill = gerp_group)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = gerp_colours, name = NULL) +
  annotate(
    "text", x = 5.5, y = Inf,
    label = paste0("Quasi-Poisson GLM P = ", format_p_value(cs_p_value)),
    vjust = 1.5, size = 2.6
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of credible sets", y = "Proportion of eGenes") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 9),
    legend.position = "top",
    legend.text = element_text(size = 7),
    plot.margin = margin(t = 8, r = 5.5, b = 5.5, l = 5.5)
  )

top_row <- maf_plot + gerp_effect_plot + credible_set_plot +
  plot_layout(widths = c(1, 1, 1.5))
figure <- top_row / panel + plot_layout(heights = c(0.9, 1))
ggsave(cfg$output_file, figure, width = 285, height = 145, units = "mm", dpi = 600)
message("Wrote Figure 4 with ", nrow(plot_data), " fine-mapped DEG-variant rows to ", cfg$output_file)
