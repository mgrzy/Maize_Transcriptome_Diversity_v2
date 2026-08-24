#!/usr/bin/env Rscript
# Extract full OSCA association rows indexed by results/summary/osca/eqtl_padjusted.csv.gz.
# This script streams both inputs; it does not load all ~69 million hits in memory.

suppressPackageStartupMessages(library(data.table))

## RStudio settings: set the working directory to the project root, then source this file.
osca_dir <- "results/full/osca"
adjusted_file <- "results/summary/osca/eqtl_padjusted.csv.gz"
output_file <- "results/summary/osca/eqtl_padjusted_full_results.tsv.gz"
lines_per_chunk <- 500000L

if (!file.exists(adjusted_file)) stop("Adjusted-result index is missing: ", adjusted_file)
if (!dir.exists(osca_dir)) stop("OSCA directory is missing: ", osca_dir)
if (file.exists(output_file)) stop("Output already exists; rename or remove it first: ", output_file)

batch_files <- list.files(osca_dir, pattern = "^batch[0-9]+\\.txt\\.gz$", full.names = TRUE)
batch_number <- as.integer(sub("^batch([0-9]+)\\.txt\\.gz$", "\\1", basename(batch_files)))
batch_files <- batch_files[order(batch_number)]
if (!length(batch_files)) stop("No compressed OSCA batch files found.")

tmp_dir <- tempfile("osca-adjusted-index-")
dir.create(tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

## 1. Partition the adjusted index by OSCA batch.  Keeping only row number and
## adjusted p-value metadata makes each subsequent batch lookup memory-efficient.
index_connection <- gzfile(adjusted_file, open = "rt")
on.exit(close(index_connection), add = TRUE)
index_header <- readLines(index_connection, n = 1L)
index_columns <- strsplit(index_header, ",", fixed = TRUE)[[1L]]
needed_columns <- c("file_index", "row_in_file", "p", "abs_dist", "dist_group_index", "dist_group", "adj_p")
if (!all(needed_columns %in% index_columns)) stop("Adjusted-result file does not have the expected columns.")

repeat {
  lines <- readLines(index_connection, n = lines_per_chunk)
  if (!length(lines)) break
  # With `header = FALSE`, fread applies character `select` to its temporary
  # V1, V2, ... names rather than to `col.names`.  Read the small index rows
  # first, then subset by the names assigned through `col.names`.
  index_chunk <- fread(text = lines, header = FALSE, col.names = index_columns)
  index_chunk <- index_chunk[, ..needed_columns]
  for (current_index in unique(index_chunk$file_index)) {
    index_path <- file.path(tmp_dir, paste0("batch_", current_index, ".tsv"))
    fwrite(index_chunk[file_index == current_index], index_path,
           sep = "\t", append = file.exists(index_path), col.names = !file.exists(index_path))
  }
}
close(index_connection)

## 2. Stream every OSCA batch and write only indexed records.  row_in_file is
## the data-row number used by 8.IHW.R (the header is not counted).
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
output_connection <- gzfile(output_file, open = "wt")
on.exit(close(output_connection), add = TRUE)
header_written <- FALSE
total_written <- 0L

for (j in seq_along(batch_files)) {
  index_path <- file.path(tmp_dir, paste0("batch_", j, ".tsv"))
  if (!file.exists(index_path)) next
  message("Extracting ", basename(batch_files[j]), " (file index ", j, ")")
  wanted <- fread(index_path)
  setkey(wanted, row_in_file)
  input_connection <- gzfile(batch_files[j], open = "rt")
  batch_header <- readLines(input_connection, n = 1L)
  if (!header_written) {
    writeLines(paste0(batch_header, "\tfile_index\trow_in_file\tadj_p\tabs_dist\tdist_group"), output_connection)
    header_written <- TRUE
  }
  row_offset <- 0L
  repeat {
    lines <- readLines(input_connection, n = lines_per_chunk)
    if (!length(lines)) break
    row_ids <- row_offset + seq_along(lines)
    positions <- match(row_ids, wanted$row_in_file, nomatch = 0L)
    keep <- which(positions > 0L)
    if (length(keep)) {
      metrics <- wanted[positions[keep]]
      output_lines <- paste(lines[keep], j, row_ids[keep], metrics$adj_p, metrics$abs_dist, metrics$dist_group, sep = "\t")
      writeLines(output_lines, output_connection)
      total_written <- total_written + length(keep)
    }
    row_offset <- row_offset + length(lines)
  }
  close(input_connection)
}
close(output_connection)
message("Wrote ", format(total_written, big.mark = ","), " full adjusted OSCA records to ", output_file)
