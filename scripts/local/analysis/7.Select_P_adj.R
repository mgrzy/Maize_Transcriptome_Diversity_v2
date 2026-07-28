#!/usr/bin/env Rscript

# Combined OSCA covariate diagnostics for IHW
#
# Reads one compressed OSCA batch at a time and calculates:
#   1. distance-to-TSS versus -log10(p) density
#   2. MAF versus -log10(p) density
#   3. total test counts in each IHW covariate group
#   4. full p-value histograms by distance and MAF group
#   5. null-tail p-value histograms by distance and MAF group
#   6. total and per-file test counts
#
# It never combines all association rows into one R object.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =============================================================================
# 1. User settings
# =============================================================================

input_dir <- "results/full/osca"
output_dir <- "results/summary/ihw_covariate_diagnostics"

# Expected filenames: batch1.txt.gz, batch2.txt.gz, ..., batch10.txt.gz
file_pattern <- "^batch[0-9]+\\.txt\\.gz$"

setDTthreads(12)

max_dist <- 2e6

# Plotting-only cap. Original p-values are used for histograms and IHW counts.
max_logp <- 100

# Resolution of the two-dimensional density plots.
n_dist_bins <- 400L
n_maf_bins <- 200L
n_logp_bins <- 250L

# Fixed equal-width IHW groups.
n_ihw_dist_groups <- 10L
n_ihw_maf_groups <- 10L

# P-value histogram bins.
p_hist_breaks <- seq(0, 1, by = 0.025)
n_p_hist_bins <- length(p_hist_breaks) - 1L

# Beginning of the null-tail diagnostic.
null_tail_start <- 0.05

# Output dimensions.
plot_width <- 13
density_plot_height <- 5.5
histogram_plot_height <- 10
group_count_plot_height <- 5.5
plot_dpi <- 300

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. Find and order input files
# =============================================================================

lf <- list.files(
  path = input_dir,
  pattern = file_pattern,
  full.names = TRUE
)

if (length(lf) == 0L) {
  stop(
    "No files found in '",
    input_dir,
    "' matching ",
    file_pattern
  )
}

batch_number <- as.integer(
  sub(
    "^batch([0-9]+)\\.txt\\.gz$",
    "\\1",
    basename(lf)
  )
)

if (anyNA(batch_number)) {
  stop("Could not extract batch numbers from all filenames.")
}

lf <- lf[order(batch_number)]

message("Files to process:")
message(paste0("  ", basename(lf), collapse = "\n"))

# =============================================================================
# 3. Density-plot bins
# =============================================================================

dist_breaks <- seq(0, max_dist, length.out = n_dist_bins + 1L)
maf_breaks <- seq(0, 0.5, length.out = n_maf_bins + 1L)
logp_breaks <- seq(0, max_logp, length.out = n_logp_bins + 1L)

dist_counts <- numeric(n_dist_bins * n_logp_bins)
maf_counts <- numeric(n_maf_bins * n_logp_bins)

# =============================================================================
# 4. IHW groups
# =============================================================================

ihw_dist_breaks <- seq(
  0,
  max_dist,
  length.out = n_ihw_dist_groups + 1L
)

ihw_maf_breaks <- seq(
  0,
  0.5,
  length.out = n_ihw_maf_groups + 1L
)

dist_group_levels <- sprintf("D%02d", seq_len(n_ihw_dist_groups))
maf_group_levels <- sprintf("M%02d", seq_len(n_ihw_maf_groups))

dist_group_ranges <- sprintf(
  "%.1f-%.1f Mb",
  ihw_dist_breaks[-length(ihw_dist_breaks)] / 1e6,
  ihw_dist_breaks[-1L] / 1e6
)

maf_group_ranges <- sprintf(
  "%.2f-%.2f",
  ihw_maf_breaks[-length(ihw_maf_breaks)],
  ihw_maf_breaks[-1L]
)

# These are the m_groups vectors needed for incomplete-list IHW.
m_groups_dist <- numeric(n_ihw_dist_groups)
m_groups_maf <- numeric(n_ihw_maf_groups)

# Rows = covariate groups; columns = p-value histogram bins.
hist_dist <- matrix(
  0,
  nrow = n_ihw_dist_groups,
  ncol = n_p_hist_bins
)

hist_maf <- matrix(
  0,
  nrow = n_ihw_maf_groups,
  ncol = n_p_hist_bins
)

# =============================================================================
# 5. Running totals
# =============================================================================

total_rows_read <- 0
total_valid_p_tests <- 0
total_distance_tests <- 0
total_maf_tests <- 0
total_invalid_p <- 0
total_p_equal_zero <- 0
total_logp_capped_for_plot <- 0

file_summary <- vector("list", length(lf))

# =============================================================================
# 6. Process one batch at a time
# =============================================================================

for (j in seq_along(lf)) {
  
  file <- lf[j]
  
  message(
    "\nProcessing ",
    basename(file),
    " [",
    j,
    "/",
    length(lf),
    "]"
  )
  
  dt <- fread(
    file,
    select = c("BP", "Freq", "Probe_bp", "p"),
    showProgress = TRUE
  )
  
  n_rows_file <- nrow(dt)
  total_rows_read <- total_rows_read + n_rows_file
  
  valid_p <- (
    is.finite(dt$p) &
      dt$p >= 0 &
      dt$p <= 1
  )
  
  n_valid_p_file <- sum(valid_p)
  n_invalid_p_file <- n_rows_file - n_valid_p_file
  
  total_valid_p_tests <- total_valid_p_tests + n_valid_p_file
  total_invalid_p <- total_invalid_p + n_invalid_p_file
  
  if (n_valid_p_file == 0L) {
    warning("No valid p-values in ", basename(file))
    
    file_summary[[j]] <- data.table(
      file = basename(file),
      rows_read = n_rows_file,
      valid_p_tests = 0,
      invalid_p_rows = n_invalid_p_file,
      distance_tests = 0,
      maf_tests = 0,
      p_equal_zero = 0,
      logp_capped_for_plot = 0
    )
    
    rm(dt, valid_p)
    gc()
    next
  }
  
  pval <- dt$p[valid_p]
  bp <- dt$BP[valid_p]
  probe_bp <- dt$Probe_bp[valid_p]
  freq <- dt$Freq[valid_p]
  
  rm(dt, valid_p)
  gc()
  
  abs_dist <- abs(bp - probe_bp)
  maf <- pmin(freq, 1 - freq)
  
  # Plotting only. p-values themselves are not modified.
  logp <- -log10(pval)
  
  capped_for_plot <- (
    !is.finite(logp) |
      logp > max_logp
  )
  
  n_zero_file <- sum(pval == 0)
  n_capped_file <- sum(capped_for_plot)
  
  total_p_equal_zero <- total_p_equal_zero + n_zero_file
  total_logp_capped_for_plot <- (
    total_logp_capped_for_plot +
      n_capped_file
  )
  
  logp[capped_for_plot] <- max_logp
  
  p_hist_bin <- findInterval(
    pval,
    p_hist_breaks,
    rightmost.closed = TRUE,
    all.inside = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # 6a. Absolute distance to TSS
  # ---------------------------------------------------------------------------
  
  keep_dist <- (
    is.finite(abs_dist) &
      abs_dist >= 0 &
      abs_dist <= max_dist &
      is.finite(logp)
  )
  
  n_distance_file <- sum(keep_dist)
  total_distance_tests <- total_distance_tests + n_distance_file
  
  if (n_distance_file > 0L) {
    
    xbin_dist <- findInterval(
      abs_dist[keep_dist],
      dist_breaks,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
    
    ybin_dist <- findInterval(
      logp[keep_dist],
      logp_breaks,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
    
    density_cell_dist <- xbin_dist +
      (ybin_dist - 1L) * n_dist_bins
    
    dist_counts <- dist_counts + tabulate(
      density_cell_dist,
      nbins = length(dist_counts)
    )
    
    dist_group <- findInterval(
      abs_dist[keep_dist],
      ihw_dist_breaks,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
    
    dist_p_hist_bin <- p_hist_bin[keep_dist]
    
    m_groups_dist <- m_groups_dist + tabulate(
      dist_group,
      nbins = n_ihw_dist_groups
    )
    
    hist_cell_dist <- dist_group +
      (dist_p_hist_bin - 1L) * n_ihw_dist_groups
    
    hist_dist <- hist_dist + matrix(
      tabulate(
        hist_cell_dist,
        nbins = n_ihw_dist_groups * n_p_hist_bins
      ),
      nrow = n_ihw_dist_groups,
      ncol = n_p_hist_bins
    )
    
    rm(
      xbin_dist,
      ybin_dist,
      density_cell_dist,
      dist_group,
      dist_p_hist_bin,
      hist_cell_dist
    )
  }
  
  # ---------------------------------------------------------------------------
  # 6b. Minor allele frequency
  # ---------------------------------------------------------------------------
  
  keep_maf <- (
    is.finite(freq) &
      freq >= 0 &
      freq <= 1 &
      is.finite(maf) &
      maf >= 0 &
      maf <= 0.5 &
      is.finite(logp)
  )
  
  n_maf_file <- sum(keep_maf)
  total_maf_tests <- total_maf_tests + n_maf_file
  
  if (n_maf_file > 0L) {
    
    xbin_maf <- findInterval(
      maf[keep_maf],
      maf_breaks,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
    
    ybin_maf <- findInterval(
      logp[keep_maf],
      logp_breaks,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
    
    density_cell_maf <- xbin_maf +
      (ybin_maf - 1L) * n_maf_bins
    
    maf_counts <- maf_counts + tabulate(
      density_cell_maf,
      nbins = length(maf_counts)
    )
    
    maf_group <- findInterval(
      maf[keep_maf],
      ihw_maf_breaks,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
    
    maf_p_hist_bin <- p_hist_bin[keep_maf]
    
    m_groups_maf <- m_groups_maf + tabulate(
      maf_group,
      nbins = n_ihw_maf_groups
    )
    
    hist_cell_maf <- maf_group +
      (maf_p_hist_bin - 1L) * n_ihw_maf_groups
    
    hist_maf <- hist_maf + matrix(
      tabulate(
        hist_cell_maf,
        nbins = n_ihw_maf_groups * n_p_hist_bins
      ),
      nrow = n_ihw_maf_groups,
      ncol = n_p_hist_bins
    )
    
    rm(
      xbin_maf,
      ybin_maf,
      density_cell_maf,
      maf_group,
      maf_p_hist_bin,
      hist_cell_maf
    )
  }
  
  file_summary[[j]] <- data.table(
    file = basename(file),
    rows_read = n_rows_file,
    valid_p_tests = n_valid_p_file,
    invalid_p_rows = n_invalid_p_file,
    distance_tests = n_distance_file,
    maf_tests = n_maf_file,
    p_equal_zero = n_zero_file,
    logp_capped_for_plot = n_capped_file
  )
  
  rm(
    pval,
    bp,
    probe_bp,
    freq,
    abs_dist,
    maf,
    logp,
    capped_for_plot,
    p_hist_bin,
    keep_dist,
    keep_maf
  )
  
  gc()
}

# =============================================================================
# 7. Check count consistency
# =============================================================================

stopifnot(
  isTRUE(all.equal(sum(m_groups_dist), total_distance_tests))
)

stopifnot(
  isTRUE(all.equal(sum(m_groups_maf), total_maf_tests))
)

# =============================================================================
# 8. Save total counts and IHW group definitions
# =============================================================================

file_summary <- rbindlist(
  file_summary,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  file_summary,
  file.path(output_dir, "per_file_test_counts.tsv"),
  sep = "\t"
)

run_summary <- data.table(
  metric = c(
    "number_of_input_files",
    "total_rows_read",
    "total_valid_p_tests",
    "total_invalid_p_rows",
    "total_distance_tests",
    "total_maf_tests",
    "p_values_equal_zero",
    "values_capped_for_density_plot",
    "max_log10_p_for_density_plot",
    "number_of_distance_ihw_groups",
    "number_of_maf_ihw_groups",
    "number_of_p_histogram_bins"
  ),
  value = c(
    length(lf),
    total_rows_read,
    total_valid_p_tests,
    total_invalid_p,
    total_distance_tests,
    total_maf_tests,
    total_p_equal_zero,
    total_logp_capped_for_plot,
    max_logp,
    n_ihw_dist_groups,
    n_ihw_maf_groups,
    n_p_hist_bins
  )
)

fwrite(
  run_summary,
  file.path(output_dir, "run_summary.tsv"),
  sep = "\t"
)

writeLines(
  c(
    paste0("Total rows read: ", total_rows_read),
    paste0("Total valid p-value tests: ", total_valid_p_tests),
    paste0(
      "Total tests with valid distance covariate: ",
      total_distance_tests
    ),
    paste0(
      "Total tests with valid MAF covariate: ",
      total_maf_tests
    )
  ),
  file.path(output_dir, "total_number_of_tests.txt")
)

names(m_groups_dist) <- dist_group_levels
names(m_groups_maf) <- maf_group_levels

dist_group_counts <- data.table(
  group = dist_group_levels,
  lower_bp = ihw_dist_breaks[-length(ihw_dist_breaks)],
  upper_bp = ihw_dist_breaks[-1L],
  range = dist_group_ranges,
  n_tests = m_groups_dist,
  fraction_of_distance_tests = (
    m_groups_dist /
      sum(m_groups_dist)
  )
)

maf_group_counts <- data.table(
  group = maf_group_levels,
  lower_maf = ihw_maf_breaks[-length(ihw_maf_breaks)],
  upper_maf = ihw_maf_breaks[-1L],
  range = maf_group_ranges,
  n_tests = m_groups_maf,
  fraction_of_maf_tests = (
    m_groups_maf /
      sum(m_groups_maf)
  )
)

fwrite(
  dist_group_counts,
  file.path(output_dir, "distance_ihw_group_counts.tsv"),
  sep = "\t"
)

fwrite(
  maf_group_counts,
  file.path(output_dir, "maf_ihw_group_counts.tsv"),
  sep = "\t"
)

saveRDS(
  list(
    covariate = "absolute_distance_to_TSS",
    breaks = ihw_dist_breaks,
    levels = dist_group_levels,
    ranges = dist_group_ranges,
    m_groups = m_groups_dist,
    total_tests = sum(m_groups_dist)
  ),
  file.path(output_dir, "ihw_distance_group_definition.rds")
)

saveRDS(
  list(
    covariate = "minor_allele_frequency",
    breaks = ihw_maf_breaks,
    levels = maf_group_levels,
    ranges = maf_group_ranges,
    m_groups = m_groups_maf,
    total_tests = sum(m_groups_maf)
  ),
  file.path(output_dir, "ihw_maf_group_definition.rds")
)

# =============================================================================
# 9. Helper functions
# =============================================================================

counts_to_grid <- function(
    counts,
    n_x_bins,
    x_breaks,
    n_y_bins,
    y_breaks,
    x_name
) {
  
  cell <- which(counts > 0)
  
  xbin <- ((cell - 1L) %% n_x_bins) + 1L
  ybin <- ((cell - 1L) %/% n_x_bins) + 1L
  
  result <- data.table(
    x = (
      x_breaks[xbin] +
        x_breaks[xbin + 1L]
    ) / 2,
    logp = (
      y_breaks[ybin] +
        y_breaks[ybin + 1L]
    ) / 2,
    count = counts[cell]
  )
  
  setnames(result, "x", x_name)
  result
}

histogram_matrix_to_dt <- function(
    count_matrix,
    group_levels,
    group_ranges,
    m_groups
) {
  
  n_groups <- nrow(count_matrix)
  n_bins <- ncol(count_matrix)
  
  result <- data.table(
    group_index = rep(
      seq_len(n_groups),
      times = n_bins
    ),
    p_bin = rep(
      seq_len(n_bins),
      each = n_groups
    ),
    count = as.vector(count_matrix)
  )
  
  result[, `:=`(
    group = factor(
      group_levels[group_index],
      levels = group_levels
    ),
    group_range = group_ranges[group_index],
    group_total = m_groups[group_index],
    p_left = p_hist_breaks[p_bin],
    p_right = p_hist_breaks[p_bin + 1L],
    p_mid = (
      p_hist_breaks[p_bin] +
        p_hist_breaks[p_bin + 1L]
    ) / 2
  )]
  
  result[, proportion_all := fifelse(
    group_total > 0,
    count / group_total,
    NA_real_
  )]
  
  result[]
}

# =============================================================================
# 10. Build and save density tables
# =============================================================================

dist_grid <- counts_to_grid(
  counts = dist_counts,
  n_x_bins = n_dist_bins,
  x_breaks = dist_breaks,
  n_y_bins = n_logp_bins,
  y_breaks = logp_breaks,
  x_name = "abs_dist"
)

maf_grid <- counts_to_grid(
  counts = maf_counts,
  n_x_bins = n_maf_bins,
  x_breaks = maf_breaks,
  n_y_bins = n_logp_bins,
  y_breaks = logp_breaks,
  x_name = "maf"
)

fwrite(
  dist_grid,
  file.path(output_dir, "distance_logp_density.tsv.gz"),
  sep = "\t"
)

fwrite(
  maf_grid,
  file.path(output_dir, "maf_logp_density.tsv.gz"),
  sep = "\t"
)

# =============================================================================
# 11. Build and save histogram tables
# =============================================================================

dist_hist <- histogram_matrix_to_dt(
  count_matrix = hist_dist,
  group_levels = dist_group_levels,
  group_ranges = dist_group_ranges,
  m_groups = m_groups_dist
)

maf_hist <- histogram_matrix_to_dt(
  count_matrix = hist_maf,
  group_levels = maf_group_levels,
  group_ranges = maf_group_ranges,
  m_groups = m_groups_maf
)

dist_hist[
  p_left >= null_tail_start,
  proportion_tail := count / sum(count),
  by = group
]

maf_hist[
  p_left >= null_tail_start,
  proportion_tail := count / sum(count),
  by = group
]

fwrite(
  dist_hist,
  file.path(output_dir, "distance_pvalue_histogram_counts.tsv.gz"),
  sep = "\t"
)

fwrite(
  maf_hist,
  file.path(output_dir, "maf_pvalue_histogram_counts.tsv.gz"),
  sep = "\t"
)

# =============================================================================
# 12. Density plots
# =============================================================================

maximum_count <- max(
  dist_grid$count,
  maf_grid$count,
  na.rm = TRUE
)

count_scale <- function() {
  scale_fill_viridis_c(
    option = "C",
    trans = "log10",
    limits = c(1, maximum_count),
    breaks = breaks_log(n = 6),
    labels = label_number(
      scale_cut = cut_short_scale()
    ),
    oob = squish,
    name = "Number of tests"
  )
}

density_subtitle <- paste0(
  "All tests; values above -log10(p) = ",
  max_logp,
  " are capped for plotting"
)

p_dist_density <- ggplot(
  dist_grid,
  aes(
    x = abs_dist,
    y = logp,
    fill = count
  )
) +
  geom_raster() +
  count_scale() +
  scale_x_continuous(
    limits = c(0, max_dist),
    breaks = seq(0, max_dist, by = 0.5e6),
    labels = label_number(
      scale = 1e-6,
      suffix = " Mb"
    ),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, max_logp),
    expand = expansion(mult = 0)
  ) +
  labs(
    title = "Association strength by distance to TSS",
    subtitle = density_subtitle,
    x = "Absolute distance to TSS",
    y = expression(-log[10](p))
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title.position = "plot")

p_maf_density <- ggplot(
  maf_grid,
  aes(
    x = maf,
    y = logp,
    fill = count
  )
) +
  geom_raster() +
  count_scale() +
  scale_x_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, max_logp),
    expand = expansion(mult = 0)
  ) +
  labs(
    title = "Association strength by minor allele frequency",
    subtitle = density_subtitle,
    x = "Minor allele frequency",
    y = expression(-log[10](p))
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title.position = "plot")

p_density_combined <- (
  p_dist_density |
    p_maf_density
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# =============================================================================
# 13. IHW group-size plots
# =============================================================================

p_dist_group_counts <- ggplot(
  dist_group_counts,
  aes(
    x = factor(range, levels = range),
    y = n_tests,
    fill = group
  )
) +
  geom_col() +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_y_continuous(
    labels = label_number(
      scale_cut = cut_short_scale()
    )
  ) +
  labs(
    title = "Number of tests per distance group",
    subtitle = paste0(
      "Total tests = ",
      comma(sum(m_groups_dist))
    ),
    x = "Absolute distance to TSS",
    y = "Number of tests"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )

p_maf_group_counts <- ggplot(
  maf_group_counts,
  aes(
    x = factor(range, levels = range),
    y = n_tests,
    fill = group
  )
) +
  geom_col() +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_y_continuous(
    labels = label_number(
      scale_cut = cut_short_scale()
    )
  ) +
  labs(
    title = "Number of tests per MAF group",
    subtitle = paste0(
      "Total tests = ",
      comma(sum(m_groups_maf))
    ),
    x = "Minor allele frequency",
    y = "Number of tests"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )

p_group_counts_combined <- (
  p_dist_group_counts |
    p_maf_group_counts
)

# =============================================================================
# 14. Facet labels
# =============================================================================

format_n <- label_number(
  scale_cut = cut_short_scale(),
  accuracy = 0.1
)

dist_facet_labels <- setNames(
  paste0(
    dist_group_ranges,
    "\nN = ",
    format_n(m_groups_dist)
  ),
  dist_group_levels
)

maf_facet_labels <- setNames(
  paste0(
    maf_group_ranges,
    "\nN = ",
    format_n(m_groups_maf)
  ),
  maf_group_levels
)

# =============================================================================
# 15. Full p-value histograms
# =============================================================================

p_dist_hist_full <- ggplot(
  dist_hist,
  aes(
    x = p_mid,
    y = proportion_all,
    fill = group
  )
) +
  geom_col(width = diff(p_hist_breaks)[1]) +
  facet_wrap(
    ~group,
    ncol = 5,
    labeller = as_labeller(dist_facet_labels)
  ) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "P-value distributions by distance to TSS",
    subtitle = "Bar height is the fraction of all tests within each group",
    x = "P-value",
    y = "Fraction of tests"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    plot.title.position = "plot"
  )

p_maf_hist_full <- ggplot(
  maf_hist,
  aes(
    x = p_mid,
    y = proportion_all,
    fill = group
  )
) +
  geom_col(width = diff(p_hist_breaks)[1]) +
  facet_wrap(
    ~group,
    ncol = 5,
    labeller = as_labeller(maf_facet_labels)
  ) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "P-value distributions by minor allele frequency",
    subtitle = "Bar height is the fraction of all tests within each group",
    x = "P-value",
    y = "Fraction of tests"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    plot.title.position = "plot"
  )

p_hist_full_combined <- (
  p_dist_hist_full /
    p_maf_hist_full
)

# =============================================================================
# 16. Null-tail p-value histograms
# =============================================================================

dist_tail <- dist_hist[
  p_left >= null_tail_start &
    is.finite(proportion_tail)
]

maf_tail <- maf_hist[
  p_left >= null_tail_start &
    is.finite(proportion_tail)
]

n_tail_bins <- sum(
  p_hist_breaks[-length(p_hist_breaks)] >= null_tail_start
)

expected_tail_fraction <- 1 / n_tail_bins

p_dist_hist_tail <- ggplot(
  dist_tail,
  aes(
    x = p_mid,
    y = proportion_tail,
    fill = group
  )
) +
  geom_col(width = diff(p_hist_breaks)[1]) +
  geom_hline(
    yintercept = expected_tail_fraction,
    linetype = 2
  ) +
  facet_wrap(
    ~group,
    ncol = 5,
    labeller = as_labeller(dist_facet_labels)
  ) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(
    limits = c(null_tail_start, 1),
    breaks = c(null_tail_start, 0.5, 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "Null-tail p-value distributions by distance to TSS",
    subtitle = paste0(
      "P-values >= ",
      null_tail_start,
      "; dashed line indicates a flat tail"
    ),
    x = "P-value",
    y = "Fraction within the p-value tail"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    plot.title.position = "plot"
  )

p_maf_hist_tail <- ggplot(
  maf_tail,
  aes(
    x = p_mid,
    y = proportion_tail,
    fill = group
  )
) +
  geom_col(width = diff(p_hist_breaks)[1]) +
  geom_hline(
    yintercept = expected_tail_fraction,
    linetype = 2
  ) +
  facet_wrap(
    ~group,
    ncol = 5,
    labeller = as_labeller(maf_facet_labels)
  ) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_x_continuous(
    limits = c(null_tail_start, 1),
    breaks = c(null_tail_start, 0.5, 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "Null-tail p-value distributions by minor allele frequency",
    subtitle = paste0(
      "P-values >= ",
      null_tail_start,
      "; dashed line indicates a flat tail"
    ),
    x = "P-value",
    y = "Fraction within the p-value tail"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    plot.title.position = "plot"
  )

p_hist_tail_combined <- (
  p_dist_hist_tail /
    p_maf_hist_tail
)

# =============================================================================
# 17. Save plots
# =============================================================================

ggsave(
  file.path(output_dir, "covariate_logp_density.png"),
  p_density_combined,
  width = plot_width,
  height = density_plot_height,
  dpi = plot_dpi
)

ggsave(
  file.path(output_dir, "covariate_logp_density.pdf"),
  p_density_combined,
  width = plot_width,
  height = density_plot_height
)

ggsave(
  file.path(output_dir, "ihw_group_test_counts.png"),
  p_group_counts_combined,
  width = plot_width,
  height = group_count_plot_height,
  dpi = plot_dpi
)

ggsave(
  file.path(output_dir, "ihw_group_test_counts.pdf"),
  p_group_counts_combined,
  width = plot_width,
  height = group_count_plot_height
)

ggsave(
  file.path(output_dir, "ihw_pvalue_histograms_full.png"),
  p_hist_full_combined,
  width = plot_width,
  height = histogram_plot_height,
  dpi = plot_dpi
)

ggsave(
  file.path(output_dir, "ihw_pvalue_histograms_full.pdf"),
  p_hist_full_combined,
  width = plot_width,
  height = histogram_plot_height
)

ggsave(
  file.path(output_dir, "ihw_pvalue_histograms_null_tail.png"),
  p_hist_tail_combined,
  width = plot_width,
  height = histogram_plot_height,
  dpi = plot_dpi
)

ggsave(
  file.path(output_dir, "ihw_pvalue_histograms_null_tail.pdf"),
  p_hist_tail_combined,
  width = plot_width,
  height = histogram_plot_height
)

# =============================================================================
# 18. Print results
# =============================================================================

message("\nFinished.")
message("Output directory: ", normalizePath(output_dir))
message("Total rows read: ", comma(total_rows_read))
message("Total valid p-value tests: ", comma(total_valid_p_tests))
message("Distance IHW tests: ", comma(total_distance_tests))
message("MAF IHW tests: ", comma(total_maf_tests))
message("Exact zero p-values: ", comma(total_p_equal_zero))
message(
  "Values capped only for density plotting: ",
  comma(total_logp_capped_for_plot)
)

print(run_summary)
print(dist_group_counts)
print(maf_group_counts)

print(p_density_combined)
print(p_group_counts_combined)
print(p_hist_full_combined)
print(p_hist_tail_combined)
