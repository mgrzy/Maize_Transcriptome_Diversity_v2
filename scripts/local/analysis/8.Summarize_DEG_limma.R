### Script: analysis/8.Summarize_DEG_limma.R
### Purpose: Combine all Limma DEG contrasts into one gene-level summary table.
### Inputs: results/full/DEG/*.csv from 7.DEG_limma.R
### Outputs:
###   - results/summary/DEG/DEG_limma_summary.csv
###   - results/summary/DEG/DEG_limma_contrast_summary.csv
### Usage: Rscript scripts/local/analysis/8.Summarize_DEG_limma.R

input_dir <- "results/full/DEG"
output_dir <- "results/summary/DEG"
output_file <- file.path(output_dir, "DEG_limma_summary.csv")
contrast_summary_file <- file.path(output_dir, "DEG_limma_contrast_summary.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result_files <- list.files(
  input_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)
result_files <- result_files[!grepl("_samples\\.csv$", result_files)]

if (length(result_files) == 0L) {
  stop("No Limma DEG result files found in ", input_dir)
}

summarise_contrast <- function(path) {
  contrast_name <- sub("\\.csv$", "", basename(path))
  result <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  missing_columns <- setdiff(required_columns, names(result))

  if (length(missing_columns) > 0L) {
    stop(
      "Missing expected Limma columns in ", basename(path), ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  result <- result[, required_columns, drop = FALSE]
  result$DEG <- ifelse(
    is.na(result$adj.P.Val),
    NA,
    result$adj.P.Val < 0.05
  )
  names(result)[-1L] <- paste(contrast_name, names(result)[-1L], sep = "_")
  result
}

contrast_tables <- lapply(result_files, summarise_contrast)
summary_table <- Reduce(
  function(x, y) merge(x, y, by = "gene", all = TRUE, sort = FALSE),
  contrast_tables
)
summary_table <- summary_table[order(summary_table$gene), , drop = FALSE]

write.csv(summary_table, output_file, row.names = FALSE, quote = FALSE)

contrast_summary <- do.call(
  rbind,
  lapply(result_files, function(path) {
    contrast_name <- sub("\\.csv$", "", basename(path))
    result <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

    data.frame(
      contrast = contrast_name,
      n_tested = sum(!is.na(result$adj.P.Val)),
      n_DEG = sum(result$adj.P.Val < 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
contrast_summary <- contrast_summary[order(contrast_summary$contrast), , drop = FALSE]
write.csv(contrast_summary, contrast_summary_file, row.names = FALSE, quote = FALSE)

message("Wrote ", nrow(summary_table), " genes to ", output_file)
message("Wrote contrast summary to ", contrast_summary_file)
