#!/usr/bin/env Rscript
# Merge SuSiE batch results and retain the lead PIP variant from each credible set.
#
# Usage:
#   Rscript scripts/local/procesing/6.Merge_SusieR_BatchResults.R
#   Rscript scripts/local/procesing/6.Merge_SusieR_BatchResults.R --susie-dir results/full/SusieR
#
# Outputs:
#   <susie-dir>/PIP_all_batches.tsv.gz
#   <susie-dir>/CS_all_batches.tsv.gz
#   results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz

suppressPackageStartupMessages(library(data.table))

script_args <- commandArgs(trailingOnly = TRUE)
script_file_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
# `--file` is present with Rscript but not with source().  In the latter case,
# source() stores the file name in the calling frame.
if (length(script_file_arg)) {
  script_path <- normalizePath(script_file_arg[[1L]])
} else {
  source_files <- vapply(rev(seq_len(sys.nframe())), function(i) {
    ofile <- sys.frame(i)$ofile
    if (is.null(ofile)) "" else as.character(ofile)[[1L]]
  }, character(1L))
  source_files <- source_files[nzchar(source_files)]
  if (!length(source_files)) stop("Cannot determine the script location; use --susie-dir and --summary-dir.")
  script_path <- normalizePath(source_files[[1L]])
}
project_root <- normalizePath(file.path(dirname(script_path), "../../.."))

options <- list(
  susie_dir = file.path(project_root, "results/full/SusieR/"),
  summary_dir = file.path(project_root, "results/summary/SusieR")
)

usage <- function() cat(
  "Usage: Rscript 6.Merge_SusieR_BatchResults.R [options]\n\n",
  "Options:\n",
  "  --susie-dir DIR    Directory containing batchResults (default: results/full/SusieR)\n",
  "  --summary-dir DIR  Output directory for lead variants (default: results/summary/SusieR)\n",
  sep = ""
)

parse_args <- function(args, values) {
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) { usage(); quit(status = 0L) }
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- gsub("-", "_", substring(key, 3L))
    if (!name %in% names(values)) stop("Unknown option: ", key)
    if (i == length(args)) stop("Missing value for ", key)
    values[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values
}

options <- parse_args(script_args, options)
batch_dir <- file.path(options$susie_dir, "batchResults")
if (!dir.exists(batch_dir)) stop("Batch result directory does not exist: ", batch_dir)

find_batches <- function(prefix) {
  files <- list.files(batch_dir, pattern = paste0("^", prefix, ".*\\.(tsv|txt)(\\.gz)?$"),
                      full.names = TRUE, ignore.case = TRUE)
  # Do not accidentally re-read merged output if it was placed in batchResults.
  files[!grepl("all_batches", basename(files), fixed = TRUE)]
}

read_and_merge <- function(files, required_columns, label) {
  if (!length(files)) stop("No ", label, " files found in ", batch_dir)
  tables <- lapply(files, function(file) {
    x <- fread(file)
    missing <- setdiff(required_columns, names(x))
    if (length(missing)) stop("Missing column(s) in ", file, ": ", paste(missing, collapse = ", "))
    x[, source_batch_file := basename(file)]
    x
  })
  rbindlist(tables, use.names = TRUE, fill = TRUE)
}

pip_files <- find_batches("PIP")
cs_files <- find_batches("CS")
pip <- read_and_merge(pip_files, c("gene", "cs", "variable_prob", "SNP"), "PIP batch")
cs <- read_and_merge(cs_files, c("gene", "cs"), "credible-set batch")

setorderv(pip, c("gene", "cs", "variable_prob", "SNP"), c(1L, 1L, -1L, 1L), na.last = TRUE)
setorderv(cs, c("gene", "cs"))

dir.create(options$susie_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(options$summary_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(pip, file.path(options$susie_dir, "PIP_all_batches.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(cs, file.path(options$susie_dir, "CS_all_batches.tsv.gz"), sep = "\t", compress = "gzip")

# Retain a highest-PIP variant per credible set.  When the maximum PIP is tied,
# sample one of the tied SNPs at random.
lead_pip <- pip[, {
  max_pip <- max(variable_prob, na.rm = TRUE)
  .SD[sample(which(variable_prob == max_pip), size = 1L)]
}, by = .(gene, cs)]
lead_pip[, source_batch_file := NULL]
lead_pip <- merge(lead_pip, cs[, !"source_batch_file"], by = c("gene", "cs"), all.x = TRUE, sort = FALSE)
setorderv(lead_pip, c("gene", "cs"))
fwrite(lead_pip,
       file.path(options$summary_dir, "PIP_lead_variant_per_credible_set.tsv.gz"),
       sep = "\t", compress = "gzip")

message("Merged ", length(pip_files), " PIP and ", length(cs_files), " CS batch files.")
message("PIP rows: ", nrow(pip), "; CS rows: ", nrow(cs), "; lead variants: ", nrow(lead_pip))
