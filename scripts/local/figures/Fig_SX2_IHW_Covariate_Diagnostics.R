## IHW covariate diagnostic plots
##
## RStudio use:
## 1. Set the working directory to the project root.
## 2. Run each section below, or source this entire file.
## 3. First run scripts/local/analysis/7.Select_P_adj.R to create the tables.

## =============================================================================
## 1. Packages and user settings
## =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

# Change this only if your RStudio working directory is not the project root.
diagnostic_dir <- "results/summary/ihw_covariate_diagnostics"

# Figure size in inches and raster-image resolution.
plot_width <- 13
density_plot_height <- 5.5
histogram_plot_height <- 10
plot_dpi <- 300
base_font_size <- 8

## =============================================================================
## 2. Load saved tables
## =============================================================================

files_needed <- file.path(diagnostic_dir, c(
  "distance_logp_density.tsv.gz",
  "maf_logp_density.tsv.gz",
  "distance_ihw_group_counts.tsv",
  "maf_ihw_group_counts.tsv",
  "distance_pvalue_histogram_counts.tsv.gz",
  "maf_pvalue_histogram_counts.tsv.gz"
))

if (any(!file.exists(files_needed))) {
  stop("Diagnostic tables are missing. Run 7.Select_P_adj.R before plotting.")
}

distance_density <- fread(files_needed[1])
maf_density <- fread(files_needed[2])
distance_groups <- fread(files_needed[3])
maf_groups <- fread(files_needed[4])
distance_histogram <- fread(files_needed[5])
maf_histogram <- fread(files_needed[6])

## =============================================================================
## 3. Density plots: distance/MAF versus association strength
## =============================================================================

maximum_logp <- ceiling(max(c(distance_density$logp, maf_density$logp)))
maximum_count <- max(c(distance_density$count, maf_density$count))

count_scale <- scale_fill_viridis_c(
  option = "C", trans = "log10", limits = c(1, maximum_count),
  labels = label_number(scale_cut = cut_short_scale()), name = "Number of tests"
)

distance_density_plot <- ggplot(distance_density, aes(x = abs_dist, y = logp, fill = count)) +
  geom_raster() +
  count_scale +
  scale_x_continuous(limits = c(0, max(distance_density$abs_dist)), breaks = pretty(c(0, max(distance_density$abs_dist))), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, maximum_logp), expand = expansion(mult = 0)) +
  labs(title = "Association strength by distance to TSS", x = "Absolute distance to TSS", y = expression(-log[10](p))) +
  theme_classic(base_size = base_font_size) +
  theme(plot.title.position = "plot")

maf_density_plot <- ggplot(maf_density, aes(x = maf, y = logp, fill = count)) +
  geom_raster() +
  count_scale +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, maximum_logp), expand = expansion(mult = 0)) +
  labs(title = "Association strength by minor allele frequency", x = "Minor allele frequency", y = expression(-log[10](p))) +
  theme_classic(base_size = base_font_size) +
  theme(plot.title.position = "plot")
density_plot <- distance_density_plot + maf_density_plot + plot_layout(guides = "collect") & theme(legend.position = "top")

## =============================================================================
## 4. Number of tests in each IHW covariate group
## =============================================================================

distance_group_plot <- ggplot(distance_groups, aes(x = factor(range, levels = range), y = n_tests, fill = group)) +
  geom_col() +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Absolute distance to TSS", y = "Number of tests") +
  theme_classic(base_size = base_font_size) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title.position = "plot")
maf_group_plot <- ggplot(maf_groups, aes(x = factor(range, levels = range), y = n_tests, fill = group)) +
  geom_col() +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Minor allele frequency", y = "Number of tests") +
  theme_classic(base_size = base_font_size) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title.position = "plot")

group_count_plot <- distance_group_plot | maf_group_plot


## =============================================================================
## 5. P-value histogram plots
## =============================================================================

distance_labels <- unique(distance_histogram[, .(group, group_range, group_total)])
distance_labels <- setNames(paste0(distance_labels$group_range, "\nN = ", label_number(scale_cut = cut_short_scale())(distance_labels$group_total)), distance_labels$group)
maf_labels <- unique(maf_histogram[, .(group, group_range, group_total)])
maf_labels <- setNames(paste0(maf_labels$group_range, "\nN = ", label_number(scale_cut = cut_short_scale())(maf_labels$group_total)), maf_labels$group)

distance_full_plot <- ggplot(distance_histogram, aes(p_mid, proportion_all, fill = group)) +
  geom_col(width = distance_histogram$p_right[1] - distance_histogram$p_left[1]) +
  facet_wrap(~group, ncol = 5, labeller = as_labeller(distance_labels)) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  scale_y_continuous(labels = label_percent()) +
  labs( x = "P-value", y = "Fraction of tests") +
  theme_classic(base_size = base_font_size) +
  theme(strip.background = element_blank(), strip.text = element_text(size = base_font_size), plot.title.position = "plot")
maf_full_plot <- ggplot(maf_histogram, aes(p_mid, proportion_all, fill = group)) +
  geom_col(width = maf_histogram$p_right[1] - maf_histogram$p_left[1]) +
  facet_wrap(~group, ncol = 5, labeller = as_labeller(maf_labels)) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  scale_y_continuous(labels = label_percent()) +
  labs(x = "P-value", y = "Fraction of tests") +
  theme_classic(base_size = base_font_size) +
  theme(strip.background = element_blank(), strip.text = element_text(size = base_font_size), plot.title.position = "plot")
full_histogram_plot <- distance_full_plot / maf_full_plot

distance_tail <- distance_histogram[p_left >= .05 & is.finite(proportion_tail)]
maf_tail <- maf_histogram[p_left >= .05 & is.finite(proportion_tail)]
distance_tail_plot <- ggplot(distance_tail, aes(p_mid, proportion_tail, fill = group)) +
  geom_col(width = distance_tail$p_right[1] - distance_tail$p_left[1]) +
  geom_hline(yintercept = 1 / length(unique(distance_tail$p_mid)), linetype = 2) +
  facet_wrap(~group, ncol = 5, labeller = as_labeller(distance_labels)) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(limits = c(.05, 1), breaks = c(.05, .5, 1)) +
  scale_y_continuous(labels = label_percent()) +
  labs(title = "Null-tail p-values by distance to TSS", x = "P-value", y = "Fraction within p-value tail") +
  theme_classic(base_size = base_font_size) +
  theme(strip.background = element_blank(), strip.text = element_text(size = base_font_size), plot.title.position = "plot")
maf_tail_plot <- ggplot(maf_tail, aes(p_mid, proportion_tail, fill = group)) +
  geom_col(width = maf_tail$p_right[1] - maf_tail$p_left[1]) +
  geom_hline(yintercept = 1 / length(unique(maf_tail$p_mid)), linetype = 2) +
  facet_wrap(~group, ncol = 5, labeller = as_labeller(maf_labels)) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(limits = c(.05, 1), breaks = c(.05, .5, 1)) +
  scale_y_continuous(labels = label_percent()) +
  labs(title = "Null-tail p-values by minor allele frequency", x = "P-value", y = "Fraction within p-value tail") +
  theme_classic(base_size = base_font_size) +
  theme(strip.background = element_blank(), strip.text = element_text(size = base_font_size), plot.title.position = "plot")
null_tail_histogram_plot <- distance_tail_plot / maf_tail_plot

## =============================================================================
## 6. Inspect plots interactively, then save PNG and PDF copies
## =============================================================================

print(density_plot)
print(group_count_plot)
print(full_histogram_plot)
print(null_tail_histogram_plot)

ggsave(file.path(diagnostic_dir, "covariate_logp_density.png"), density_plot, width = plot_width, height = density_plot_height, dpi = plot_dpi)
ggsave(file.path(diagnostic_dir, "covariate_logp_density.pdf"), density_plot, width = plot_width, height = density_plot_height)
ggsave(file.path(diagnostic_dir, "ihw_group_test_counts.png"), group_count_plot, width = plot_width, height = density_plot_height, dpi = plot_dpi)
ggsave(file.path(diagnostic_dir, "ihw_group_test_counts.pdf"), group_count_plot, width = plot_width, height = density_plot_height)
ggsave(file.path(diagnostic_dir, "ihw_pvalue_histograms_full.png"), full_histogram_plot, width = plot_width, height = histogram_plot_height, dpi = plot_dpi)
ggsave(file.path(diagnostic_dir, "ihw_pvalue_histograms_full.pdf"), full_histogram_plot, width = plot_width, height = histogram_plot_height)
ggsave(file.path(diagnostic_dir, "ihw_pvalue_histograms_null_tail.png"), null_tail_histogram_plot, width = plot_width, height = histogram_plot_height, dpi = plot_dpi)
ggsave(file.path(diagnostic_dir, "ihw_pvalue_histograms_null_tail.pdf"), null_tail_histogram_plot, width = plot_width, height = histogram_plot_height)

message("Plots saved in: ", normalizePath(diagnostic_dir))

###

density_plot / (distance_group_plot | distance_full_plot) + plot_annotation(tag_levels = "A")

ggsave("results/figures/supp/Fig_IHW.png", width = 12, height = 8)
