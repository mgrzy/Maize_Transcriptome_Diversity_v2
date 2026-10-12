#!/usr/bin/env Rscript
# Run from the repository root:
# Rscript --vanilla scripts/local/figures/Fig4_EI_presentation.R
# Reuse Figure 4's complete E-I example, including the gene-track panel G.
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(data.table)
  library(bigmemory)
})

set.seed(42)
in_example <- FALSE
for (source_expr in parse("scripts/local/figures/Fig4.R")) {
  object <- if (is.call(source_expr) && identical(source_expr[[1L]], as.name("<-")) &&
                is.symbol(source_expr[[2L]])) as.character(source_expr[[2L]]) else ""
  if (object == "example_gene") in_example <- TRUE
  if (in_example || object %in% c("stop_with", "cfg", "myPal")) eval(source_expr)
  if (object == "example_panels") break
}
if (!exists("example_panels")) stop("Example panels are missing from Fig4.R.")

# Use R=reference (dosage 0), A=alternative (dosage 2), in CS1/CS2/CS3 order. Show the original
# four main combinations and pool every remaining combination, including
# heterozygous calls, as "Other" in both frequency and expression panels.
snp_order <- fine %>% arrange(cs) %>% pull(SNP)
if (length(snp_order) != 3L || anyDuplicated(snp_order)) {
  stop("Expected three distinct credible-set lead SNPs.")
}
all_haplotype_data <- geno_samples %>%
  filter(if_all(all_of(snp_order), ~ !is.na(.))) %>%
  mutate(
    genotype_combination = apply(select(., all_of(snp_order)), 1, paste, collapse = "/"),
    homozygous = if_all(all_of(snp_order), ~ . %in% c(0, 2)),
    haplotype = if_else(
      homozygous,
      apply(select(., all_of(snp_order)), 1,
            function(x) paste(ifelse(x == 0, "R", "A"), collapse = "/")),
      "Heterozygous"
    )
  )
if (any(!as.matrix(select(all_haplotype_data, all_of(snp_order))) %in% c(0, 1, 2))) {
  stop("Unexpected lead SNP genotype dosage.")
}
homozygous_data <- all_haplotype_data %>% filter(homozygous)
haplotype_levels <- c("R/A/R", "A/A/R", "R/R/R", "A/R/R", "Other")
haplotype_colours <- setNames(myPal[1:5], haplotype_levels)
all_haplotype_data <- all_haplotype_data %>%
  mutate(haplotype = if_else(haplotype %in% haplotype_levels, haplotype, "Other"))
expression_order <- all_haplotype_data %>%
  group_by(haplotype) %>%
  summarise(median_expression = median(expression_log2_cpm), .groups = "drop") %>%
  arrange(median_expression, haplotype) %>%
  pull(haplotype)
all_haplotype_data <- all_haplotype_data %>%
  mutate(haplotype = factor(haplotype, levels = expression_order))
haplotype_frequency_data <- all_haplotype_data %>%
  count(Group, haplotype, name = "n_lines", .drop = FALSE) %>%
  group_by(Group) %>%
  mutate(total_lines = sum(n_lines), frequency = n_lines / total_lines) %>%
  ungroup()
stopifnot(sum(haplotype_frequency_data$n_lines) == nrow(all_haplotype_data))
haplotype_frequency_plot <- ggplot(
  haplotype_frequency_data, aes(Group, frequency, fill = haplotype)
) +
  geom_col(width = .72, colour = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = haplotype_colours, drop = FALSE,
                    name = "CS1 / CS2 / CS3", guide = guide_legend(ncol = 3)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Haplotype frequency") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = grid::unit(3, "mm"))
haplotype_summary <- all_haplotype_data %>%
  group_by(haplotype) %>%
  summarise(n_lines = n(), n_tropical = sum(Group == "Tropical"),
            genotype_combinations = paste(sort(unique(genotype_combination)), collapse = ", "),
            mean_expression = mean(expression_log2_cpm),
            median_expression = median(expression_log2_cpm), .groups = "drop")
haplotype_effect_plot <- ggplot(
  all_haplotype_data, aes(haplotype, expression_log2_cpm, fill = haplotype)
) +
  geom_boxplot(width = .7, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = .14, height = 0, seed = 42),
              size = .7, alpha = .45, colour = "grey20") +
  geom_text(data = haplotype_summary,
             aes(x = haplotype, y = Inf, label = paste0("n=", n_lines)),
             inherit.aes = FALSE, vjust = 1.2, size = 2.3) +
  scale_fill_manual(values = haplotype_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(.04, .12))) +
  labs(x = "CS1 / CS2 / CS3\nR = reference; A = alternative",
       y = expression(log[2]~CPM)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_text(size = 8))
example_panels <- group_effect_plot + example_manhattan +
  haplotype_frequency_plot + haplotype_effect_plot +
  plot_layout(widths = c(1, 2, 1.35, 1.8))

figure_file <- "results/figures/presentation/Fig4_EI.png"
dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, example_panels, width = 34, height = 9, units = "cm",
       dpi = 300, device = ragg::agg_png)
write_tsv(haplotype_summary, "results/figures/presentation/Fig4_EI_haplotype_summary.tsv")
write_tsv(haplotype_frequency_data,
          "results/figures/presentation/Fig4_EI_simplified_haplotype_frequencies.tsv")

# Companion plot: compare CS3 dosages within each homozygous CS1/CS2 background.
cs3_data <- homozygous_data %>%
  mutate(background = paste(ifelse(.data[[snp_order[1]]] == 0, "R", "A"),
                            ifelse(.data[[snp_order[2]]] == 0, "R", "A"), sep = "/"),
         background = factor(background, levels = c("R/R", "R/A", "A/R", "A/A")),
         cs3_genotype = factor(if_else(.data[[snp_order[3]]] == 0, "R", "A"),
                               levels = c("R", "A")))
cs3_counts <- cs3_data %>% count(background, cs3_genotype, name = "n_lines")
cs3_plot <- ggplot(cs3_data, aes(cs3_genotype, expression_log2_cpm, fill = cs3_genotype)) +
  geom_boxplot(width = .65, outlier.shape = NA, linewidth = .35) +
  geom_jitter(position = position_jitter(width = .12, height = 0, seed = 42),
              size = .8, alpha = .4, colour = "grey20") +
  geom_text(data = cs3_counts,
             aes(x = cs3_genotype, y = Inf, label = paste0("n=", n_lines)),
             inherit.aes = FALSE, vjust = 1.2, size = 3) +
  scale_fill_manual(values = c("R" = myPal[1], "A" = myPal[5]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(.04, .12))) +
  facet_wrap(~ background, nrow = 1,
              labeller = labeller(background = function(x) paste0("CS1 / CS2: ", x))) +
  labs(x = "CS3 genotype\nR = reference; A = alternative", y = expression(log[2]~CPM)) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey95", colour = "grey70"))
ggsave("results/figures/presentation/Fig4_CS3_expression_by_background.png",
       cs3_plot, width = 22, height = 7, units = "cm", dpi = 300,
       device = ragg::agg_png)
message("Saved presentation figure: ", figure_file)
