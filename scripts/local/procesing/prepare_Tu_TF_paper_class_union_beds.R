#!/usr/bin/env Rscript
# Create one union BED per paper_class from the Tu/Mejia-Guerra TF consensus
# peak sets.  Each output contains the union of all reproducible consensus
# peaks belonging to TFs assigned to that paper class.
#
# Inputs (overridable with environment variables):
#   TF_CONSENSUS_NAMED_BED
#     data/processed/chipseq/Tu_2020_TF_consensus_v5/all_TF_consensus_named.bed
#   TF_FUNCTIONAL_CLASSIFICATION
#     data/processed/chromatin/Tu_2020_TF_functional_classification.xlsx
#   TF_CLASSIFICATION_SHEET (default: 2)
#   OUTPUT_DIR
#     data/processed/chipseq/Tu_2020_TF_consensus_v5/by_paper_class
#
# Output files are named paper_class_<class>_union.bed; for example,
# paper_class_PS_union.bed and paper_class_metabolism_union.bed.  BED files
# contain three columns (chromosome, 0-based start, end), have no header, and
# merge both overlapping and directly adjacent intervals.
#
# Run:
#   Rscript scripts/local/procesing/prepare_Tu_TF_paper_class_union_beds.R

suppressPackageStartupMessages({
  library(data.table)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  named_bed = Sys.getenv(
    "TF_CONSENSUS_NAMED_BED",
    "data/processed/chipseq/Tu_2020_TF_consensus_v5/all_TF_consensus_named.bed"
  ),
  classification = Sys.getenv(
    "TF_FUNCTIONAL_CLASSIFICATION",
    "data/processed/chromatin/Tu_2020_TF_functional_classification.xlsx"
  ),
  classification_sheet = Sys.getenv("TF_CLASSIFICATION_SHEET", "2"),
  output_dir = Sys.getenv(
    "OUTPUT_DIR",
    "data/processed/chipseq/Tu_2020_TF_consensus_v5/by_paper_class"
  )
)

if (!file.exists(cfg$named_bed)) stop_with("Named TF consensus BED does not exist: %s", cfg$named_bed)
if (!file.exists(cfg$classification)) stop_with("TF functional-classification workbook does not exist: %s", cfg$classification)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop_with("Package 'readxl' is required to read the TF functional-classification workbook.")
}

sheet_number <- suppressWarnings(as.integer(cfg$classification_sheet))
if (is.na(sheet_number) || sheet_number < 1L) {
  stop_with("TF_CLASSIFICATION_SHEET must be a positive sheet number; received: %s", cfg$classification_sheet)
}

make_file_safe <- function(x) {
  x <- trimws(tolower(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("(^_+|_+$)", "", x)
  ifelse(nzchar(x), x, "unlabelled")
}

read_classification <- function(path, sheet) {
  classification <- as.data.table(readxl::read_excel(path, sheet = sheet))
  required_columns <- c("TF", "paper_class")
  missing_columns <- setdiff(required_columns, names(classification))
  if (length(missing_columns)) {
    stop_with(
      "Classification sheet is missing required column(s): %s",
      paste(missing_columns, collapse = ", ")
    )
  }

  classification <- classification[, ..required_columns]
  classification[, (required_columns) := lapply(.SD, function(x) trimws(as.character(x))), .SDcols = required_columns]
  classification <- classification[!is.na(TF) & nzchar(TF) & !is.na(paper_class) & nzchar(paper_class)]
  if (!nrow(classification)) stop_with("Classification sheet contains no TF/paper_class assignments.")
  if (anyDuplicated(classification$TF)) stop_with("Classification sheet contains duplicate TF identifiers.")
  classification
}

read_named_bed <- function(path) {
  peaks <- tryCatch(
    fread(path, header = FALSE, select = 1:4, colClasses = "character", showProgress = FALSE),
    error = function(error) stop_with("Could not read four BED columns from %s: %s", path, error$message)
  )
  if (ncol(peaks) < 4L) stop_with("Named TF consensus BED has fewer than four columns: %s", path)
  setnames(peaks, c("chromosome", "start", "end", "TF"))
  peaks[, `:=`(
    chromosome = trimws(chromosome),
    start = suppressWarnings(as.integer(start)),
    end = suppressWarnings(as.integer(end)),
    TF = trimws(TF)
  )]
  invalid <- is.na(peaks$chromosome) | !nzchar(peaks$chromosome) |
    is.na(peaks$start) | is.na(peaks$end) | peaks$start < 0L | peaks$end <= peaks$start |
    is.na(peaks$TF) | !nzchar(peaks$TF)
  if (any(invalid)) {
    stop_with("Named TF consensus BED contains %d invalid row(s): %s", sum(invalid), path)
  }
  peaks
}

merge_bed_intervals <- function(intervals) {
  # BED is 0-based and half-open.  start <= previous cumulative end therefore
  # combines both overlapping and book-ended intervals, as bedtools merge does.
  chromosome_number <- suppressWarnings(as.integer(sub("^chr", "", intervals$chromosome, ignore.case = TRUE)))
  intervals[, `:=`(
    chromosome_number = chromosome_number,
    chromosome_missing = is.na(chromosome_number),
    chromosome_label = tolower(chromosome)
  )]
  setorder(intervals, chromosome_missing, chromosome_number, chromosome_label, start, end)
  intervals[, group_id := cumsum(start > shift(cummax(end), fill = -1L)), by = chromosome]
  intervals[, .(start = min(start), end = max(end)), by = .(chromosome, group_id)][, group_id := NULL][]
}

message("Reading TF functional classifications from ", cfg$classification)
classification <- read_classification(cfg$classification, sheet_number)
message("Loaded ", nrow(classification), " TF/paper_class assignments across ", uniqueN(classification$paper_class), " paper classes.")

message("Reading TF-labelled consensus peaks from ", cfg$named_bed)
peaks <- read_named_bed(cfg$named_bed)
unclassified_tfs <- setdiff(sort(unique(peaks$TF)), classification$TF)
if (length(unclassified_tfs)) {
  stop_with(
    "%d TF(s) in the named consensus BED lack a paper_class assignment: %s",
    length(unclassified_tfs), paste(unclassified_tfs, collapse = ", ")
  )
}

peaks <- merge(peaks, classification, by = "TF", all.x = TRUE, sort = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

class_order <- sort(unique(classification$paper_class))
qc_rows <- vector("list", length(class_order))
for (i in seq_along(class_order)) {
  target_class <- class_order[[i]]
  class_tfs <- sort(classification[paper_class == target_class, TF])
  class_peaks <- peaks[paper_class == target_class, .(chromosome, start, end)]
  union_peaks <- if (nrow(class_peaks)) merge_bed_intervals(class_peaks) else class_peaks
  output_file <- file.path(
    cfg$output_dir,
    sprintf("paper_class_%s_union.bed", make_file_safe(target_class))
  )
  fwrite(union_peaks, output_file, sep = "\t", col.names = FALSE)

  qc_rows[[i]] <- data.table(
    paper_class = target_class,
    output_file = basename(output_file),
    n_TFs = length(class_tfs),
    TFs = paste(class_tfs, collapse = ";"),
    n_input_intervals = nrow(class_peaks),
    input_covered_bp = if (nrow(class_peaks)) sum(as.numeric(class_peaks$end) - as.numeric(class_peaks$start)) else 0,
    n_union_intervals = nrow(union_peaks),
    union_covered_bp = if (nrow(union_peaks)) sum(as.numeric(union_peaks$end) - as.numeric(union_peaks$start)) else 0
  )
  message(
    target_class, ": ", length(class_tfs), " TF(s), ", nrow(class_peaks),
    " input intervals -> ", nrow(union_peaks), " union intervals (", output_file, ")"
  )
}

qc <- rbindlist(qc_rows)
qc_file <- file.path(cfg$output_dir, "paper_class_union_QC.tsv")
fwrite(qc, qc_file, sep = "\t")
message("Wrote ", nrow(qc), " paper-class union BED file(s) and QC summary: ", qc_file)
