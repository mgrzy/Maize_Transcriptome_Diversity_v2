#!/usr/bin/env Rscript
# Plot per-TF matched-permutation enrichment results from
# 17.Fine_Mapped_TF_ChIP_Enrichment.R.
#
# Example:
#   TF_ENRICHMENT_FILE=results/summary/fine_mapped_TF_ChIP_enrichment/TF_overlap_matched_permutation.tsv \
#   OUTPUT_DIR=results/figures/TF_ChIP_enrichment \
#   Rscript scripts/local/analysis/18.Plot_TF_ChIP_Enrichment.R
#
# The main panel ranks TFs by effect size, not raw overlap count. With only 100
# permutations, many P values are tied at their minimum, so the script does not
# use significance stars. Re-run the enrichment with >=5,000 permutations for
# final inference.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

input_file <- Sys.getenv(
  "TF_ENRICHMENT_FILE",
  "results/summary/fine_mapped_TF_ChIP_enrichment/TF_overlap_matched_permutation.tsv"
)
output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "results/figures/TF_ChIP_enrichment"
)
top_n <- suppressWarnings(as.integer(Sys.getenv("TOP_N", "15")))
min_observed <- suppressWarnings(as.integer(Sys.getenv("MIN_OBSERVED", "200")))

if (!file.exists(input_file)) stop_with("Input file does not exist: %s", input_file)
if (is.na(top_n) || top_n < 1L) stop_with("TOP_N must be a positive integer.")
if (is.na(min_observed) || min_observed < 0L) {
  stop_with("MIN_OBSERVED must be a non-negative integer.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tf <- fread(input_file)
required <- c(
  "TF", "n_consensus_peaks", "consensus_covered_bp",
  "observed_unique_variants", "expected_unique_variants",
  "null_q025", "null_q975", "fold_enrichment", "log2_fold_enrichment",
  "p_enrichment", "p_enrichment_fdr_bh", "n_permutations"
)
missing_columns <- setdiff(required, names(tf))
if (length(missing_columns)) {
  stop_with("Missing required columns: %s", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(tf$TF)) stop_with("TF identifiers are not unique.")
if (tf[, any(!is.finite(fold_enrichment) | fold_enrichment <= 0)]) {
  stop_with("fold_enrichment must contain finite positive values.")
}

tf[, TF_family := sub("_TF.*$", "", TF)]
tf[, covered_Mb := consensus_covered_bp / 1e6]

# Convert the permutation null-count interval to an observed/expected-ratio
# envelope. This is a descriptive permutation envelope, not a parametric CI.
tf[, `:=`(
  fold_lower = fifelse(null_q975 > 0, observed_unique_variants / null_q975, NA_real_),
  fold_upper = fifelse(null_q025 > 0, observed_unique_variants / null_q025, NA_real_)
)]
tf[, `:=`(
  log2_fold_lower = log2(fold_lower),
  log2_fold_upper = log2(fold_upper)
)]

n_perm <- unique(tf$n_permutations)
if (length(n_perm) != 1L) warning("TFs have different permutation counts.")
if (min(n_perm) < 5000L) {
  warning(
    "Only ", min(n_perm), " permutations were used. Effect-size plots are valid, ",
    "but final empirical P values should use at least 5,000 permutations."
  )
}

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "right"
  )

# Main panel: select robustly estimated TFs by a prespecified minimum observed
# count, then show the largest fold enrichments. Selection is not based on P.
main_tf <- tf[observed_unique_variants >= min_observed]
setorder(main_tf, -log2_fold_enrichment, -observed_unique_variants, TF)
main_tf <- head(main_tf, top_n)
main_tf[, TF_plot := factor(TF, levels = rev(TF))]

p_main <- ggplot(main_tf, aes(x = log2_fold_enrichment, y = TF_plot)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = log2_fold_lower, xmax = log2_fold_upper),
    height = 0, linewidth = 0.45, colour = "grey45"
  ) +
  geom_point(
    aes(size = observed_unique_variants, colour = TF_family),
    alpha = 0.9
  ) +
  scale_size_continuous(range = c(2.5, 7), name = "Overlapping\nvariants") +
  labs(
    title = "Fine-mapped cis-eQTL enrichment in TF-binding sites",
    subtitle = sprintf(
      "Top %d TFs among those with at least %s observed variants",
      nrow(main_tf), format(min_observed, big.mark = ",")
    ),
    x = expression(log[2] * " fold enrichment over matched controls"),
    y = NULL,
    colour = "TF family"
  ) +
  base_theme

# Supplement: all TFs, ranked by effect size. Colour indicates the observed
# count, avoiding an unreadable 20+ category family legend.
all_tf <- copy(tf)
setorder(all_tf, log2_fold_enrichment, TF)
all_tf[, TF_plot := factor(TF, levels = TF)]
p_all <- ggplot(all_tf, aes(x = log2_fold_enrichment, y = TF_plot)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = log2_fold_lower, xmax = log2_fold_upper),
    height = 0, linewidth = 0.25, colour = "grey65"
  ) +
  geom_point(aes(colour = observed_unique_variants), size = 2) +
  scale_colour_viridis_c(option = "C", trans = "sqrt", name = "Observed\nvariants") +
  labs(
    title = "Enrichment across all 104 TF consensus peak sets",
    x = expression(log[2] * " fold enrichment over matched controls"),
    y = NULL
  ) +
  base_theme +
  theme(axis.text.y = element_text(size = 6.2))

# Diagnostic 1: raw overlap is expected to depend strongly on peak coverage.
p_coverage <- ggplot(
  tf,
  aes(x = covered_Mb, y = observed_unique_variants)
) +
  geom_point(aes(colour = fold_enrichment), alpha = 0.8, size = 2.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 0.6, colour = "grey30") +
  scale_colour_viridis_c(option = "C", name = "Fold\nenrichment") +
  labs(
    title = "Raw overlap depends on TF peak coverage",
    x = "Consensus peak coverage (Mb)",
    y = "Observed overlapping variants"
  ) +
  base_theme

# Diagnostic 2: observed versus matched-null expectation.
limits_oe <- range(c(tf$observed_unique_variants, tf$expected_unique_variants), finite = TRUE)
p_observed_expected <- ggplot(
  tf,
  aes(x = expected_unique_variants, y = observed_unique_variants)
) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey45") +
  geom_point(aes(colour = TF_family), alpha = 0.8, size = 2.4) +
  coord_equal(xlim = limits_oe, ylim = limits_oe) +
  labs(
    title = "Observed versus matched-null TF overlaps",
    x = "Expected overlapping variants",
    y = "Observed overlapping variants",
    colour = "TF family"
  ) +
  base_theme +
  theme(legend.position = "none")

# Family-level distributions are descriptive because TF peak sets overlap and
# TFs are not statistically independent. Only families represented by >=3 TFs
# are displayed.
family_counts <- tf[, .N, by = TF_family]
family_tf <- tf[family_counts[N >= 3L], on = "TF_family"]
family_order <- family_tf[, .(median_log2_enrichment = median(log2_fold_enrichment)),
                          by = TF_family][order(median_log2_enrichment)]
family_tf[, TF_family_plot := factor(TF_family, levels = family_order$TF_family)]
p_family <- ggplot(family_tf, aes(x = log2_fold_enrichment, y = TF_family_plot)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "grey90", colour = "grey35") +
  geom_jitter(height = 0.12, width = 0, size = 1.8, alpha = 0.75,
              aes(colour = TF_family)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  guides(colour = "none") +
  labs(
    title = "TF-family enrichment distributions",
    subtitle = "Descriptive only; families with at least three assayed TFs",
    x = expression(log[2] * " fold enrichment"),
    y = NULL
  ) +
  base_theme

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 400, bg = "white")
}

save_plot(p_main, "Fig3_TF_enrichment_main", 8.2, 5.8)
save_plot(p_all, "FigS_TF_enrichment_all_104", 7.5, 13.0)
save_plot(p_coverage, "FigS_TF_coverage_vs_raw_overlap", 6.4, 5.0)
save_plot(p_observed_expected, "FigS_TF_observed_vs_expected", 5.8, 5.2)
save_plot(p_family, "FigS_TF_family_enrichment", 7.0, 5.8)

fwrite(main_tf[, .(
  TF, TF_family, n_consensus_peaks, consensus_covered_bp,
  observed_unique_variants, expected_unique_variants,
  fold_enrichment, log2_fold_enrichment,
  p_enrichment, p_enrichment_fdr_bh, n_permutations
)], file.path(output_dir, "Fig3_TF_enrichment_main_data.tsv"), sep = "\t")

summary_table <- tf[, .(
  n_TFs = .N,
  n_enriched = sum(fold_enrichment > 1),
  minimum_fold = min(fold_enrichment),
  median_fold = median(fold_enrichment),
  maximum_fold = max(fold_enrichment),
  coverage_raw_overlap_Pearson_r = cor(consensus_covered_bp, observed_unique_variants),
  coverage_fold_enrichment_Pearson_r = cor(consensus_covered_bp, fold_enrichment),
  n_permutations = min(n_permutations)
)]
fwrite(summary_table, file.path(output_dir, "TF_enrichment_plot_summary.tsv"), sep = "\t")

message("Wrote main and supplementary TF-enrichment plots to: ", output_dir)
print(summary_table)
