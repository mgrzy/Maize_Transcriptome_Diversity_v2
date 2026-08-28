#!/usr/bin/env Rscript
# Manhattan plots of mean pairwise LD (mean_rsq) and LD score calculated by GCTA.
#
# Usage (from the repository root):
#   Rscript scripts/local/analysis/population_genetics/11.Plot_LD_Max_Rsq_Manhattan.R
#   Rscript scripts/local/analysis/population_genetics/11.Plot_LD_Max_Rsq_Manhattan.R \
#     results/full/LD/LDWiDiv.score.ld results/figures/supp
#   LD_METRICS=ldscore Rscript scripts/local/analysis/population_genetics/11.Plot_LD_Max_Rsq_Manhattan.R
#     # LD_METRICS can be "mean_rsq", "ldscore", or both (the default).

args <- commandArgs(trailingOnly = TRUE)
input_candidates <- c("results/full/LD/LDWiDiv.score.ld")
input_file <- if (length(args) >= 1L) args[[1L]] else input_candidates[file.exists(input_candidates)][1L]
output_dir <- if (length(args) >= 2L) args[[2L]] else "results/figures/supp"
metrics <- strsplit(Sys.getenv("LD_METRICS", "mean_rsq,ldscore"), ",", fixed = TRUE)[[1L]]
metrics <- trimws(metrics)
valid_metrics <- c("mean_rsq", "ldscore")
if (!length(metrics) || any(!metrics %in% valid_metrics)) {
  stop("LD_METRICS must be a comma-separated subset of: ", paste(valid_metrics, collapse = ", "))
}

if (is.na(input_file) || !file.exists(input_file)) {
  stop(
    "GCTA LD-score output not found. Expected one of: ",
    paste(input_candidates, collapse = ", "),
    ". Supply its path as the first argument."
  )
}

for (package in c("data.table", "ggplot2")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required; install it before running this script.")
  }
}

ld <- data.table::fread(input_file, select = c("chr", "bp", "mean_rsq", "ldscore"))
required_columns <- c("chr", "bp", "mean_rsq", "ldscore")
if (!all(required_columns %in% names(ld))) {
  stop(
    "Unexpected LD-score columns. Required columns: ",
    paste(required_columns, collapse = ", "),
    ". Found: ", paste(names(ld), collapse = ", ")
  )
}

ld[, chr := as.character(chr)]
ld[, bp := as.numeric(bp)]
ld[, mean_rsq := as.numeric(mean_rsq)]
ld[, ldscore := as.numeric(ldscore)]
ld <- ld[!is.na(chr) & nzchar(chr) & is.finite(bp) & is.finite(mean_rsq) & is.finite(ldscore) & bp >= 0 & mean_rsq >= 0 & ldscore >= 0]
if (!nrow(ld)) stop("No finite SNP position, mean_rsq, and ldscore values were found in ", input_file)

# Keep chromosome ordering numeric when possible (e.g. 1--10 for maize),
# while still supporting non-numeric chromosome labels.
chromosomes <- unique(ld$chr)
numeric_chromosomes <- suppressWarnings(as.numeric(chromosomes))
chromosomes <- chromosomes[order(is.na(numeric_chromosomes), numeric_chromosomes, chromosomes)]
ld[, chr := factor(chr, levels = chromosomes)]
data.table::setorder(ld, chr, bp)

chromosome_lengths <- ld[, .(chromosome_length = max(bp)), by = chr]
chromosome_lengths[, offset := c(0, head(cumsum(chromosome_length), -1L))]
chromosome_lengths[, centre := offset + chromosome_length / 2]
ld <- chromosome_lengths[ld, on = "chr"]
ld[, genomic_position := bp + offset]

make_manhattan_plot <- function(metric, y_label) {
  ggplot2::ggplot(ld, ggplot2::aes(x = genomic_position, y = .data[[metric]], colour = chr)) +
    ggplot2::geom_point(size = 0.12, alpha = 0.55, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = rep(c("#3B6EA5", "#E08E45"), length.out = nrow(chromosome_lengths))) +
    ggplot2::scale_x_continuous(
      breaks = chromosome_lengths$centre,
      labels = as.character(chromosome_lengths$chr),
      expand = ggplot2::expansion(mult = c(0.005, 0.005))
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.03))) +
    ggplot2::labs(x = "Chromosome", y = y_label) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(axis.line.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_files <- c(
  mean_rsq = file.path(output_dir, "LDWiDiv_mean_rsq_manhattan.png"),
  ldscore = file.path(output_dir, "LDWiDiv_ldscore_manhattan.png")
)
metric_labels <- list(
  mean_rsq = expression("Mean LD (" * r^2 * ")"),
  ldscore = "LD score"
)
for (metric in metrics) {
  ggplot2::ggsave(
    output_files[[metric]],
    make_manhattan_plot(metric, metric_labels[[metric]]),
    width = 16, height = 7, units = "in", dpi = 300, bg = "white"
  )
  gc()
}
message("Wrote Manhattan plot(s) to ", paste(output_files[metrics], collapse = " and "))
