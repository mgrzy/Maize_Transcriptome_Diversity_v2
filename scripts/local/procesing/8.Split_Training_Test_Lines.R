#!/usr/bin/env Rscript
# Create stratified training and test line lists from OSCA PCs and line metadata.
#
# Usage:
#   Rscript scripts/local/procesing/8.Split_Training_Test_Lines.R
#   Rscript scripts/local/procesing/8.Split_Training_Test_Lines.R --test-fraction 0.2 --seed 1
#
# Outputs in data/processed/meta:
#   test_lines.txt, training_lines.txt, train_test_line_split.tsv.gz

suppressPackageStartupMessages(library(data.table))

script_args <- commandArgs(trailingOnly = TRUE)
script_file_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(script_file_arg)) {
  script_path <- normalizePath(script_file_arg[[1L]])
} else {
  source_files <- vapply(rev(seq_len(sys.nframe())), function(i) {
    ofile <- sys.frame(i)$ofile
    if (is.null(ofile)) "" else as.character(ofile)[[1L]]
  }, character(1L))
  source_files <- source_files[nzchar(source_files)]
  if (!length(source_files)) stop("Cannot determine the script location; run with Rscript or source a file.")
  script_path <- normalizePath(source_files[[1L]])
}
project_root <- normalizePath(file.path(dirname(script_path), "../../.."))

options <- list(
  pc_file = file.path(project_root, "data/processed/osca/WiDiv631_PC.txt"),
  metadata_file = file.path(project_root, "data/raw/meta/Lines_meta.csv"),
  output_dir = file.path(project_root, "data/processed/meta"),
  test_fraction = 0.20,
  seed = 1L
)

usage <- function() cat(
  "Usage: Rscript 8.Split_Training_Test_Lines.R [options]\n\n",
  "Options:\n",
  "  --pc-file FILE        PC file containing a Taxa column.\n",
  "  --metadata-file FILE  Line metadata containing Taxa and Group columns.\n",
  "  --output-dir DIR      Directory for split files.\n",
  "  --test-fraction NUM   Fraction assigned to test (default: 0.20).\n",
  "  --seed INT            Random seed (default: 1).\n",
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
options$test_fraction <- as.numeric(options$test_fraction)
options$seed <- as.integer(options$seed)
if (!is.finite(options$test_fraction) || options$test_fraction <= 0 || options$test_fraction >= 1) {
  stop("--test-fraction must be strictly between 0 and 1.")
}
if (is.na(options$seed)) stop("--seed must be an integer.")
if (!file.exists(options$pc_file)) stop("PC file does not exist: ", options$pc_file)
if (!file.exists(options$metadata_file)) stop("Metadata file does not exist: ", options$metadata_file)

pcs <- fread(options$pc_file)
metadata <- fread(options$metadata_file, na.strings = c("", "NA"))
required_pc <- "Taxa"
required_metadata <- c("Taxa", "Group")
if (!all(required_pc %in% names(pcs))) stop("PC file must contain a Taxa column.")
if (!all(required_metadata %in% names(metadata))) {
  stop("Metadata file must contain: ", paste(required_metadata, collapse = ", "))
}

pcs <- unique(pcs[, .(Taxa = as.character(Taxa))])
metadata[, Taxa := as.character(Taxa)]
if (anyDuplicated(metadata$Taxa)) stop("Metadata has duplicate Taxa values; resolve them before splitting.")

# The inner join ensures only lines available in the PC file are assigned.
merged <- merge(pcs, metadata, by = "Taxa", all = FALSE, sort = FALSE)
if (!nrow(merged)) stop("No shared Taxa values between the PC and metadata files.")
if (anyNA(merged$Group) || any(!nzchar(merged$Group))) stop("Every PC line must have a non-empty Group.")
if (nrow(merged) != nrow(pcs)) {
  warning(nrow(pcs) - nrow(merged), " PC lines lack metadata and were excluded.")
}

# Largest-remainder allocation preserves group proportions and makes the total
# test set exactly round(test_fraction * number_of_lines).
group_sizes <- merged[, .N, by = Group]
group_sizes[, target := N * options$test_fraction]
group_sizes[, n_test := floor(target)]
n_to_allocate <- as.integer(round(nrow(merged) * options$test_fraction) - sum(group_sizes$n_test))
if (n_to_allocate > 0L) {
  group_sizes[, remainder := target - n_test]
  setorder(group_sizes, -remainder, Group)
  group_sizes[seq_len(n_to_allocate), n_test := n_test + 1L]
}
setkey(group_sizes, Group)
setkey(merged, Group)
merged <- group_sizes[merged]

set.seed(options$seed)
merged[, split := {
  test_index <- sample.int(.N, size = n_test[[1L]], replace = FALSE)
  ifelse(seq_len(.N) %in% test_index, "test", "training")
}, by = Group]
setorder(merged, split, Group, Taxa)

dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(merged[split == "test", .(Taxa)], file.path(options$output_dir, "test_lines.txt"), col.names = FALSE)
fwrite(merged[split == "training", .(Taxa)], file.path(options$output_dir, "training_lines.txt"), col.names = FALSE)
fwrite(merged[, .(Taxa, Group, split)],
       file.path(options$output_dir, "train_test_line_split.tsv.gz"), sep = "\t", compress = "gzip")

split_counts <- merged[, .N, by = .(Group, split)]
print(dcast(split_counts, Group ~ split, value.var = "N", fill = 0L))
message("Saved ", sum(merged$split == "test"), " test and ", sum(merged$split == "training"), " training lines to ", options$output_dir)
