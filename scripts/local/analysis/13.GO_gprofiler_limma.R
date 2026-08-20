#!/usr/bin/env Rscript
# Purpose: Run GO enrichment with g:Profiler for each Limma DEG contrast.
#
# The custom background is contrast-specific: every gene with a non-missing
# Limma adjusted P value in that contrast (i.e., every gene tested in it).
# GO is run for all DEGs together and separately for up- and downregulated
# DEGs, using adjusted P < 0.05 as the DEG threshold.
#
# Input:
#   - results/summary/DEG/DEG_limma_summary.csv
# Outputs:
#   - results/full/DEG/GO_gprofiler/<contrast>_{all,up,down}_GO.tsv
#   - results/full/DEG/GO_gprofiler/<contrast>_background_genes.tsv
#   - results/full/DEG/GO_gprofiler/GO_gprofiler_run_summary.tsv
#
# Usage:
#   Rscript scripts/local/analysis/13.GO_gprofiler_limma.R
#
# Optional settings:
#   GO_ADJ_P_THRESHOLD=0.05 GPROFILER_ORGANISM=zmays \
#     Rscript scripts/local/analysis/13.GO_gprofiler_limma.R

suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

get_env_number <- function(name, default, minimum = 0, maximum = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || !is.finite(value) || value < minimum || value > maximum) {
    stop_with("%s must be between %s and %s.", name, minimum, maximum)
  }
  value
}

cfg <- list(
  input_file = Sys.getenv("DEG_SUMMARY_FILE", "results/summary/DEG/DEG_limma_summary.csv"),
  output_dir = Sys.getenv("GO_OUTPUT_DIR", "results/full/DEG/GO_gprofiler"),
  organism = Sys.getenv("GPROFILER_ORGANISM", "zmays"),
  deg_threshold = get_env_number("GO_ADJ_P_THRESHOLD", 0.05, 0, 1),
  go_threshold = get_env_number("GPROFILER_P_THRESHOLD", 0.05, 0, 1)
)
if (!file.exists(cfg$input_file)) stop_with("DEG summary does not exist: %s", cfg$input_file)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

write_empty_go <- function(path) {
  fwrite(data.table(
    source = character(), term_id = character(), term_name = character(),
    term_size = integer(), query_size = integer(), intersection_size = integer(),
    effective_domain_size = integer(), intersection = character(), p_value = numeric()
  ), path, sep = "\t")
}

write_gene_list <- function(genes, path) {
  fwrite(data.table(gene = sort(unique(genes))), path, sep = "\t")
}

run_gprofiler <- function(query_genes, background_genes, contrast, direction) {
  output_file <- file.path(cfg$output_dir, paste0(contrast, "_", direction, "_GO.tsv"))
  query_genes <- sort(unique(query_genes))

  if (!length(query_genes)) {
    write_empty_go(output_file)
    return(list(n_query = 0L, n_terms = 0L, status = "no_query_genes", output_file = output_file))
  }

  result <- tryCatch(
    gprofiler2::gost(
      query = query_genes,
      organism = cfg$organism,
      ordered_query = FALSE,
      significant = TRUE,
      user_threshold = cfg$go_threshold,
      correction_method = "g_SCS",
      domain_scope = "custom",
      custom_bg = background_genes,
      sources = c("GO:BP", "GO:MF", "GO:CC"),
      evcodes = TRUE
    ),
    error = function(error) error
  )

  if (inherits(result, "error")) {
    warning("g:Profiler failed for ", contrast, " (", direction, "): ", conditionMessage(result))
    write_empty_go(output_file)
    return(list(
      n_query = length(query_genes), n_terms = 0L,
      status = paste0("gprofiler_error: ", conditionMessage(result)), output_file = output_file
    ))
  }

  terms <- as.data.table(result$result)
  if (!nrow(terms)) {
    write_empty_go(output_file)
    return(list(n_query = length(query_genes), n_terms = 0L, status = "no_significant_terms", output_file = output_file))
  }

  setorder(terms, p_value, source, term_id)
  fwrite(terms, output_file, sep = "\t")
  list(n_query = length(query_genes), n_terms = nrow(terms), status = "ok", output_file = output_file)
}

message("Reading Limma DEG summary: ", cfg$input_file)
deg <- fread(cfg$input_file)
if (!"gene" %in% names(deg)) stop_with("DEG summary is missing the gene column.")
deg[, gene := as.character(gene)]
if (deg[, anyDuplicated(gene)]) stop_with("DEG summary has duplicate gene IDs.")

adjusted_p_columns <- grep("_adj\\.P\\.Val$", names(deg), value = TRUE)
if (!length(adjusted_p_columns)) stop_with("No *_adj.P.Val columns were found in the DEG summary.")
contrasts <- sub("_adj\\.P\\.Val$", "", adjusted_p_columns)

required_columns <- c(
  "gene",
  paste0(contrasts, "_logFC"),
  paste0(contrasts, "_adj.P.Val")
)
missing_columns <- setdiff(required_columns, names(deg))
if (length(missing_columns)) {
  stop_with("DEG summary is missing: %s", paste(missing_columns, collapse = ", "))
}

run_summary <- rbindlist(lapply(contrasts, function(contrast) {
  adjusted_p_column <- paste0(contrast, "_adj.P.Val")
  logfc_column <- paste0(contrast, "_logFC")
  adjusted_p <- as.numeric(deg[[adjusted_p_column]])
  logfc <- as.numeric(deg[[logfc_column]])

  # The non-missing adjusted-P universe exactly mirrors the genes tested by
  # Limma for this contrast, and is supplied to g:Profiler as custom_bg.
  background <- unique(deg$gene[is.finite(adjusted_p)])
  if (!length(background)) stop_with("No tested genes found for contrast: %s", contrast)
  write_gene_list(background, file.path(cfg$output_dir, paste0(contrast, "_background_genes.tsv")))

  is_deg <- is.finite(adjusted_p) & adjusted_p < cfg$deg_threshold
  groups <- list(
    all = deg$gene[is_deg],
    up = deg$gene[is_deg & logfc > 0],
    down = deg$gene[is_deg & logfc < 0]
  )

  rbindlist(lapply(names(groups), function(direction) {
    go_result <- run_gprofiler(groups[[direction]], background, contrast, direction)
    data.table(
      contrast = contrast,
      direction = direction,
      n_background = length(background),
      n_DEG_query = go_result$n_query,
      n_significant_GO_terms = go_result$n_terms,
      status = go_result$status,
      output_file = go_result$output_file
    )
  }))
}), use.names = TRUE)

setorder(run_summary, contrast, direction)
summary_file <- file.path(cfg$output_dir, "GO_gprofiler_run_summary.tsv")
fwrite(run_summary, summary_file, sep = "\t")
message("Wrote g:Profiler GO results to ", cfg$output_dir)
message("Wrote run summary to ", summary_file)
