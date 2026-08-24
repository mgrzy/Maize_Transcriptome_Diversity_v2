#!/usr/bin/env Rscript
# Annotate fine-mapped aFCn variants with Tu and Mejia-Guerra (2020) ATAC-seq,
# replicate-supported TF consensus peaks, and the project's ATAC narrowPeak
# overlap.
#
# The sample description is read from the final column of
# Tu_Mejia-Guerra_2020_read_data.txt.  Samples containing "ATAC" contribute to
# in_ATAC_peak. TF annotations are read from the output of
# prepare_Tu_TF_consensus_peaks.sh: `in_TF_consensus` uses the union across
# reproducible TF peaks, while TF_peak_names and TF_peak_families use the
# TF-labelled consensus BED. These are replicate-supported consensus peaks,
# not formal IDR peaks. The project's narrowPeak files contribute their tissue
# name(s) (leaf and/or ear) to in_narrowPeak.
#
# Inputs (overridable with environment variables):
#   results/summary/SusieR/aFCn_MAF_se_filtered.tsv.gz
#   data/processed/chromatin/Tu_Mejia-Guerra_TF_Peaks/
#     Tu_Mejia-Guerra_2020_read_data.txt
#     <ATAC_run_accession>_region.bed[.gz]
#   data/processed/chipseq/Tu_2020_TF_consensus_v5/
#     all_TF_consensus_named.bed
#     any_TF_consensus_union.bed
#   data/processed/chromatin/Tu_2020_TF_functional_classification.xlsx
#     (sheet 2: TF_classification)
#
# Output:
#   results/summary/SusieR/aFCn_MAF_se_filtered_Tu_Mejia-Guerra_chromatin.tsv.gz
#
# Run:
#   Rscript scripts/local/analysis/17.Annotate_aFCn_Tu_Mejia-Guerra_Chromatin.R

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  afcn_path = Sys.getenv(
    "AFCN_PATH",
    "results/summary/SusieR/aFCn_MAF_se_filtered.tsv.gz"
  ),
  peak_dir = Sys.getenv(
    "TU_MEJIA_GUERRA_PEAK_DIR",
    "data/processed/chromatin/Tu_Mejia-Guerra_TF_Peaks"
  ),
  narrow_peak_dir = Sys.getenv("NARROW_PEAK_DIR", "data/processed/chromatin/Ricci"),
  metadata_path = Sys.getenv("TU_MEJIA_GUERRA_METADATA", ""),
  tf_consensus_dir = Sys.getenv(
    "TU_TF_CONSENSUS_DIR",
    "data/processed/chipseq/Tu_2020_TF_consensus_v5"
  ),
  tf_consensus_named_path = Sys.getenv("TU_TF_CONSENSUS_NAMED", ""),
  tf_consensus_union_path = Sys.getenv("TU_TF_CONSENSUS_UNION", ""),
  tf_functional_classification_path = Sys.getenv(
    "TU_TF_FUNCTIONAL_CLASSIFICATION",
    "data/processed/chromatin/Tu_2020_TF_functional_classification.xlsx"
  ),
  output_path = Sys.getenv(
    "OUTPUT_PATH",
    "results/summary/SusieR/aFCn_MAF_se_filtered_Tu_Mejia-Guerra_chromatin.tsv.gz"
  )
)
if (!nzchar(cfg$metadata_path)) {
  cfg$metadata_path <- file.path(cfg$peak_dir, "Tu_Mejia-Guerra_2020_read_data.txt")
}
if (!nzchar(cfg$tf_consensus_named_path)) {
  cfg$tf_consensus_named_path <- file.path(cfg$tf_consensus_dir, "all_TF_consensus_named.bed")
}
if (!nzchar(cfg$tf_consensus_union_path)) {
  cfg$tf_consensus_union_path <- file.path(cfg$tf_consensus_dir, "any_TF_consensus_union.bed")
}

for (path in unlist(cfg[c(
  "afcn_path", "metadata_path", "tf_consensus_named_path", "tf_consensus_union_path",
  "tf_functional_classification_path"
)])) {
  if (!file.exists(path)) stop_with("Required input file does not exist: %s", path)
}
if (!dir.exists(cfg$peak_dir)) stop_with("Peak directory does not exist: %s", cfg$peak_dir)
if (!dir.exists(cfg$narrow_peak_dir)) stop_with("narrowPeak directory does not exist: %s", cfg$narrow_peak_dir)

canonical_seqname <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chromosome[_.-]*", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  numeric_chromosome <- suppressWarnings(as.integer(x))
  out <- ifelse(
    !is.na(numeric_chromosome),
    paste0("chr", numeric_chromosome),
    paste0("chr", x)
  )
  out[is.na(x) | !nzchar(x)] <- NA_character_
  out
}

parse_variant_ranges <- function(variant_ids) {
  variant_ids <- as.character(variant_ids)
  matched <- regexec("^([^:_]+)[:_]([0-9]+)(?:[:_].*)?$", variant_ids)
  parts <- regmatches(variant_ids, matched)
  chromosome <- vapply(parts, function(x) if (length(x)) x[[2L]] else NA_character_, character(1L))
  position <- suppressWarnings(as.integer(vapply(
    parts, function(x) if (length(x)) x[[3L]] else NA_character_, character(1L)
  )))
  invalid <- is.na(chromosome) | is.na(position) | position < 1L
  if (any(invalid)) {
    stop_with(
      "Could not parse %d variant ID(s). Expected chr_position, chr:position, or an allele-suffixed version of either.",
      sum(invalid)
    )
  }
  GRanges(
    seqnames = canonical_seqname(chromosome),
    ranges = IRanges(start = position, end = position),
    strand = "*"
  )
}

first_bed_data_line <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  line_number <- 0L
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) break
    line_number <- line_number + 1L
    if (nzchar(line) && !startsWith(line, "#")) return(line_number - 1L)
  }
  stop_with("No BED data rows found in: %s", path)
}

read_peak_ranges <- function(path) {
  peaks <- fread(
    path,
    skip = first_bed_data_line(path),
    header = FALSE,
    select = 1:3,
    colClasses = "character",
    showProgress = FALSE
  )
  if (ncol(peaks) < 3L) stop_with("BED file has fewer than three columns: %s", path)
  starts <- suppressWarnings(as.integer(peaks[[2L]]))
  ends <- suppressWarnings(as.integer(peaks[[3L]]))
  chromosomes <- canonical_seqname(peaks[[1L]])
  keep <- !is.na(chromosomes) & !is.na(starts) & !is.na(ends) & starts >= 0L & ends > starts
  if (!any(keep)) return(GRanges())

  # BED intervals are 0-based, half-open; GRanges intervals are 1-based,
  # closed, hence start + 1 and unchanged end.
  GRanges(
    seqnames = chromosomes[keep],
    ranges = IRanges(start = starts[keep] + 1L, end = ends[keep]),
    strand = "*"
  )
}

annotate_binary_peak_class <- function(variant_ranges, peak_paths, label) {
  overlaps <- rep(FALSE, length(variant_ranges))
  for (i in seq_along(peak_paths)) {
    remaining <- which(!overlaps)
    if (!length(remaining)) break

    peak_ranges <- read_peak_ranges(peak_paths[[i]])
    if (length(peak_ranges)) {
      hits <- findOverlaps(variant_ranges[remaining], peak_ranges, ignore.strand = TRUE)
      if (length(hits)) overlaps[remaining[unique(queryHits(hits))]] <- TRUE
    }
    message(
      label, " peak file ", i, "/", length(peak_paths), ": ",
      sum(overlaps), "/", length(overlaps), " unique variants overlap."
    )
  }
  overlaps
}

annotate_named_peak_sets <- function(variant_ranges, peak_metadata, name_column, label) {
  feature_names <- vector("list", length(variant_ranges))
  peak_metadata <- copy(peak_metadata)
  setnames(peak_metadata, name_column, "feature_name")
  feature_groups <- split(peak_metadata, peak_metadata$feature_name)

  for (feature_index in seq_along(feature_groups)) {
    feature_name <- names(feature_groups)[[feature_index]]
    feature_paths <- feature_groups[[feature_index]]$peak_path
    feature_overlaps <- rep(FALSE, length(variant_ranges))

    # Replicates of a feature/tissue are treated as one peak set. Read each BED
    # sequentially, retaining only the variants not already hit by a replicate.
    for (path in feature_paths) {
      remaining <- which(!feature_overlaps)
      if (!length(remaining)) break
      peak_ranges <- read_peak_ranges(path)
      if (length(peak_ranges)) {
        hits <- findOverlaps(variant_ranges[remaining], peak_ranges, ignore.strand = TRUE)
        if (length(hits)) feature_overlaps[remaining[unique(queryHits(hits))]] <- TRUE
      }
    }

    hit_indices <- which(feature_overlaps)
    if (length(hit_indices)) {
      feature_names[hit_indices] <- lapply(feature_names[hit_indices], function(x) c(x, feature_name))
    }
    message(
      label, " ", feature_index, "/", length(feature_groups), " (", feature_name, "): ",
      length(hit_indices), "/", length(feature_names), " unique variants overlap."
    )
  }

  vapply(
    feature_names,
    function(x) if (length(x)) paste(sort(unique(x)), collapse = ";") else NA_character_,
    character(1L)
  )
}

read_named_consensus_ranges <- function(path) {
  peaks <- tryCatch(
    fread(
      path,
      skip = first_bed_data_line(path),
      header = FALSE,
      select = 1:4,
      colClasses = "character",
      showProgress = FALSE
    ),
    error = function(error) {
      stop_with("Could not read four BED columns from TF consensus file %s: %s", path, error$message)
    }
  )
  if (ncol(peaks) < 4L) stop_with("TF consensus BED has fewer than four columns: %s", path)

  starts <- suppressWarnings(as.integer(peaks[[2L]]))
  ends <- suppressWarnings(as.integer(peaks[[3L]]))
  chromosomes <- canonical_seqname(peaks[[1L]])
  tf_names <- trimws(as.character(peaks[[4L]]))
  valid <- !is.na(chromosomes) & !is.na(starts) & !is.na(ends) &
    starts >= 0L & ends > starts & !is.na(tf_names) & nzchar(tf_names)
  if (!all(valid)) {
    stop_with(
      "TF consensus BED contains %d invalid BED/name row(s): %s",
      sum(!valid), path
    )
  }

  # BED is 0-based and half-open; GRanges is 1-based and closed.
  ranges <- GRanges(
    seqnames = chromosomes,
    ranges = IRanges(start = starts + 1L, end = ends),
    strand = "*"
  )
  mcols(ranges)$tf_name <- tf_names
  ranges
}

annotate_named_consensus <- function(variant_ranges, consensus_path) {
  consensus_ranges <- read_named_consensus_ranges(consensus_path)
  if (!length(consensus_ranges)) return(rep(NA_character_, length(variant_ranges)))

  hits <- findOverlaps(variant_ranges, consensus_ranges, ignore.strand = TRUE)
  annotations <- rep(NA_character_, length(variant_ranges))
  if (length(hits)) {
    names_by_variant <- split(
      mcols(consensus_ranges)$tf_name[subjectHits(hits)],
      queryHits(hits)
    )
    hit_indices <- as.integer(names(names_by_variant))
    annotations[hit_indices] <- vapply(
      names_by_variant,
      function(tf_names) paste(sort(unique(tf_names)), collapse = ";"),
      character(1L)
    )
  }
  message(
    "TF consensus named BED: ", sum(!is.na(annotations)), "/", length(annotations),
    " unique variants overlap one or more reproducible TF peak(s)."
  )
  annotations
}

read_tf_functional_classification <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop_with(
      "Package 'readxl' is required to read the TF functional-classification workbook. "
    )
  }
  classification <- as.data.table(readxl::read_excel(path, sheet = 2))
  required_columns <- c("TF", "paper_class", "exploratory_broad_group")
  missing_columns <- setdiff(required_columns, names(classification))
  if (length(missing_columns)) {
    stop_with(
      "Sheet 2 of the TF functional-classification workbook is missing required column(s): %s",
      paste(missing_columns, collapse = ", ")
    )
  }

  classification <- classification[, ..required_columns]
  classification[, (required_columns) := lapply(.SD, function(x) trimws(as.character(x))), .SDcols = required_columns]
  classification <- classification[!is.na(TF) & nzchar(TF)]
  if (anyDuplicated(classification$TF)) {
    stop_with("Sheet 2 of the TF functional-classification workbook has duplicate TF identifiers.")
  }
  if (!nrow(classification)) stop_with("Sheet 2 of the TF functional-classification workbook has no TF records.")
  classification
}

annotate_tf_functional_field <- function(tf_peak_names, classification, field) {
  values_by_tf <- setNames(classification[[field]], classification$TF)
  annotations <- rep(NA_character_, length(tf_peak_names))
  hit_indices <- which(!is.na(tf_peak_names) & nzchar(tf_peak_names))
  if (!length(hit_indices)) return(annotations)

  annotations[hit_indices] <- vapply(
    strsplit(tf_peak_names[hit_indices], ";", fixed = TRUE),
    function(tf_names) {
      field_values <- values_by_tf[tf_names]
      field_values <- field_values[!is.na(field_values) & nzchar(field_values)]
      if (!length(field_values)) return(NA_character_)
      paste(sort(unique(field_values)), collapse = ";")
    },
    character(1L)
  )
  annotations
}

message("Reading aFCn variants from ", cfg$afcn_path)
afcn <- fread(cfg$afcn_path)
required_columns <- c("gene_id", "variant_id")
missing_columns <- setdiff(required_columns, names(afcn))
if (length(missing_columns)) {
  stop_with("aFCn input is missing required column(s): %s", paste(missing_columns, collapse = ", "))
}
if (anyNA(afcn$variant_id)) stop_with("aFCn input contains missing variant_id values.")

message(
  "Reading TF functional classifications from sheet 2 of ",
  cfg$tf_functional_classification_path
)
tf_functional_classification <- read_tf_functional_classification(
  cfg$tf_functional_classification_path
)
message("Loaded functional classifications for ", nrow(tf_functional_classification), " TF(s).")

message("Reading Tu and Mejia-Guerra sample metadata from ", cfg$metadata_path)
metadata <- fread(cfg$metadata_path)
required_metadata_columns <- c("run_accession")
missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
if (length(missing_metadata_columns)) {
  stop_with("Metadata is missing required column(s): %s", paste(missing_metadata_columns, collapse = ", "))
}
description <- as.character(metadata[[ncol(metadata)]])
metadata <- unique(metadata[grepl("ATAC", description, ignore.case = TRUE), .(run_accession)])
if (!nrow(metadata)) stop_with("No ATAC samples were identified from the metadata description column.")

metadata[, peak_path := file.path(cfg$peak_dir, paste0(run_accession, "_region.bed"))]
metadata[!file.exists(peak_path), peak_path := paste0(peak_path, ".gz")]
missing_peak_paths <- metadata[!file.exists(peak_path)]
if (nrow(missing_peak_paths)) {
  stop_with(
    "Missing BED file(s) for selected metadata run accession(s): %s",
    paste(missing_peak_paths$run_accession, collapse = ", ")
  )
}
message(
  "Selected ", nrow(metadata), " Tu and Mejia-Guerra ATAC peak file(s). ",
  "TF annotations will use consensus files: ", cfg$tf_consensus_union_path,
  " and ", cfg$tf_consensus_named_path, "."
)
narrow_peak_paths <- list.files(
  cfg$narrow_peak_dir,
  pattern = "\\.narrowPeak(\\.(gz|bgz))?$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (!length(narrow_peak_paths)) {
  stop_with("No narrowPeak files were found in: %s", cfg$narrow_peak_dir)
}
narrow_peak_tissues <- fifelse(
  grepl("(^|_)leaf(_|\\.)", basename(narrow_peak_paths), ignore.case = TRUE),
  "leaf",
  fifelse(
    grepl("(^|_)ear(_|\\.)", basename(narrow_peak_paths), ignore.case = TRUE),
    "ear",
    NA_character_
  )
)
if (anyNA(narrow_peak_tissues)) {
  stop_with(
    "Could not assign leaf/ear tissue from narrowPeak file name(s): %s",
    paste(basename(narrow_peak_paths[is.na(narrow_peak_tissues)]), collapse = ", ")
  )
}
narrow_peak_metadata <- data.table(peak_path = narrow_peak_paths, tissue = narrow_peak_tissues)
message(
  "Selected ", nrow(narrow_peak_metadata), " narrowPeak file(s) for: ",
  paste(sort(unique(narrow_peak_metadata$tissue)), collapse = ", "), "."
)

variants <- unique(afcn[, .(variant_id)])
message("Parsing ", nrow(variants), " unique aFCn variant position(s).")
variant_ranges <- parse_variant_ranges(variants$variant_id)

variants[, in_ATAC_peak := as.integer(annotate_binary_peak_class(
  variant_ranges,
  metadata$peak_path,
  "ATAC"
))]
variants[, in_narrowPeak := annotate_named_peak_sets(
  variant_ranges,
  narrow_peak_metadata,
  "tissue",
  "narrowPeak"
)]
variants[, in_TF_consensus := as.integer(annotate_binary_peak_class(
  variant_ranges,
  cfg$tf_consensus_union_path,
  "TF consensus union"
))]
variants[, TF_peak_names := annotate_named_consensus(
  variant_ranges,
  cfg$tf_consensus_named_path
)]
variants[, TF_peak_families := vapply(
  strsplit(TF_peak_names, ";", fixed = TRUE),
  function(tf_names) {
    if (!length(tf_names) || all(is.na(tf_names))) return(NA_character_)
    paste(sort(unique(sub("_TF.*$", "", tf_names))), collapse = ";")
  },
  character(1L)
)]
annotated_tfs <- sort(unique(unlist(strsplit(
  variants$TF_peak_names[!is.na(variants$TF_peak_names)], ";", fixed = TRUE
))))
unclassified_tfs <- setdiff(annotated_tfs, tf_functional_classification$TF)
if (length(unclassified_tfs)) {
  warning(
    length(unclassified_tfs), " consensus TF name(s) are absent from the functional-classification workbook: ",
    paste(unclassified_tfs, collapse = ", "), call. = FALSE
  )
}
variants[, TF_paper_classes := annotate_tf_functional_field(
  TF_peak_names,
  tf_functional_classification,
  "paper_class"
)]
variants[, TF_exploratory_broad_groups := annotate_tf_functional_field(
  TF_peak_names,
  tf_functional_classification,
  "exploratory_broad_group"
)]

variant_index <- match(afcn$variant_id, variants$variant_id)
afcn[, `:=`(
  in_ATAC_peak = variants$in_ATAC_peak[variant_index],
  in_narrowPeak = variants$in_narrowPeak[variant_index],
  in_TF_consensus = variants$in_TF_consensus[variant_index],
  TF_peak_names = variants$TF_peak_names[variant_index],
  TF_peak_families = variants$TF_peak_families[variant_index],
  TF_paper_classes = variants$TF_paper_classes[variant_index],
  TF_exploratory_broad_groups = variants$TF_exploratory_broad_groups[variant_index]
)]
if (anyNA(afcn$in_ATAC_peak)) {
  stop_with("Internal error: one or more aFCn rows could not be matched to a parsed variant.")
}

dir.create(dirname(cfg$output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(afcn, cfg$output_path, sep = "\t", compress = "gzip")
message(
  "Wrote ", nrow(afcn), " aFCn rows to ", cfg$output_path,
  " (ATAC overlap: ", sum(afcn$in_ATAC_peak == 1L),
  "; narrowPeak overlap: ", sum(!is.na(afcn$in_narrowPeak)),
  "; TF consensus union overlap: ", sum(afcn$in_TF_consensus == 1L),
  "; TF overlap: ", sum(!is.na(afcn$TF_peak_names)),
  "; TF paper-class annotation: ", sum(!is.na(afcn$TF_paper_classes)), ")."
)
