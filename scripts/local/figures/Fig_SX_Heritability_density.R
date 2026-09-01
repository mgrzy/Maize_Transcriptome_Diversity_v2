#!/usr/bin/env Rscript
# Compare the distributions of raw and PC-covariate-adjusted narrow-sense
# heritability across genes.
#
# Run from the project root:
#   Rscript scripts/local/figures/Fig_SX_Heritability_density.R
#
# Optional paths:
#   Rscript scripts/local/figures/Fig_SX_Heritability_density.R \
#     --input path/to/narrow_sense_h2_and_r2.tsv.gz \
#     --output results/figures/supp/Fig_Heritability_density.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_file_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(script_file_arg)) {
  script_path <- normalizePath(script_file_arg[[1L]])
} else {
  stop("Run this script with Rscript so its project location can be determined.")
}
project_root <- normalizePath(file.path(dirname(script_path), "../../.."))

config <- list(
  input = file.path(project_root, "results/summary/Heritability_and_r2/narrow_sense_h2_and_r2.tsv.gz"),
  output = file.path(project_root, "results/figures/supp/Fig_Heritability_density.png")
)

args <- commandArgs(trailingOnly = TRUE)
# Some Rscript installations retain the internal --file argument here.
args <- args[!startsWith(args, "--file=")]
if (length(args) %% 2L != 0L || (length(args) && any(!args[seq.int(1L, length(args), by = 2L)] %in% c("--input", "--output")))) {
  stop("Usage: Rscript Fig_SX_Heritability_density.R [--input FILE] [--output FILE]")
}
if (length(args)) {
  for (i in seq.int(1L, length(args), by = 2L)) {
    config[[sub("^--", "", args[[i]])]] <- args[[i + 1L]]
  }
}

if (!file.exists(config$input)) stop("Heritability input file does not exist: ", config$input)
heritability <- fread(config$input)
required_columns <- c("h2_raw", "h2_pc_adjusted")
missing_columns <- setdiff(required_columns, names(heritability))
if (length(missing_columns)) {
  stop("Heritability input is missing required column(s): ", paste(missing_columns, collapse = ", "))
}

plot_data <- melt(
  heritability[, ..required_columns],
  measure.vars = required_columns,
  variable.name = "estimate",
  value.name = "heritability"
)[is.finite(heritability)]
plot_data[, estimate := factor(
  estimate,
  levels = required_columns,
  labels = c("Raw", "PC-covariate adjusted")
)]
if (!nrow(plot_data)) stop("No finite heritability values were found.")

medians <- plot_data[, .(median_h2 = median(heritability), n_genes = .N), by = estimate]
message(
  "Plotting ", paste(sprintf("%s: n = %d", medians$estimate, medians$n_genes), collapse = "; "),
  "."
)

figure <- ggplot(plot_data, aes(x = heritability, colour = estimate, fill = estimate)) +
  geom_density(linewidth = 0.8, alpha = 0.22, adjust = 1) +
  geom_vline(
    data = medians,
    aes(xintercept = median_h2, colour = estimate),
    linetype = "dashed",
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Raw" = "#4C78A8", "PC-covariate adjusted" = "#E45756")) +
  scale_fill_manual(values = c("Raw" = "#4C78A8", "PC-covariate adjusted" = "#E45756")) +
  labs(
    x = "Narrow-sense heritability (h²)",
    y = "Density",
    colour = NULL,
    fill = NULL
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top")

dir.create(dirname(config$output), recursive = TRUE, showWarnings = FALSE)
ggsave(config$output, figure, width = 6.5, height = 4, dpi = 300)
message("Wrote figure to ", config$output)
