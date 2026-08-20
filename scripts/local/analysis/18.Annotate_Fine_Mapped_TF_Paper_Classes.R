#!/usr/bin/env Rscript
# Annotate fine-mapped variants with the Tu/Mejia-Guerra TF paper classes.
# Paper-class BEDs must first be created with:
#   bash scripts/local/procesing/prepare_Tu_TF_paper_class_union_beds.sh
#
# Inputs (overridable with environment variables):
#   FINE_MAPPED_FILE
#     results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz
#   PAPER_CLASS_BED_DIR
#     data/processed/chipseq/Tu_2020_TF_consensus_v5/by_paper_class
#   PAPER_CLASS_QC_FILE
#     <PAPER_CLASS_BED_DIR>/paper_class_union_QC.tsv
#   VARIANT_ID_COLUMN (default: SNP; CHROM and POS take precedence when present)
#   OUTPUT_FILE
#     results/summary/fine_mapped_TF_paper_class_annotations.tsv.gz
#   OVERLAP_OUTPUT_FILE
#     results/summary/fine_mapped_TF_paper_class_overlaps.tsv.gz
#
# The annotation output retains every input row and adds:
#   in_TF_paper_class_peak
#   n_TF_paper_class_peak_overlaps
#   TF_paper_classes (semicolon-delimited where a variant overlaps >1 class)
#
# Run:
#   Rscript scripts/local/analysis/18.Annotate_Fine_Mapped_TF_Paper_Classes.R

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  fine_mapped_file = Sys.getenv(
    "FINE_MAPPED_FILE",
    "results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz"
  ),
  paper_class_bed_dir = Sys.getenv(
    "PAPER_CLASS_BED_DIR",
    "data/processed/chipseq/Tu_2020_TF_consensus_v5/by_paper_class"
  ),
  paper_class_qc_file = Sys.getenv("PAPER_CLASS_QC_FILE", ""),
  variant_id_column = Sys.getenv("VARIANT_ID_COLUMN", "SNP"),
  output_file = Sys.getenv(
    "OUTPUT_FILE",
    "results/summary/fine_mapped_TF_paper_class_annotations.tsv.gz"
  ),
  overlap_output_file = Sys.getenv(
    "OVERLAP_OUTPUT_FILE",
    "results/summary/fine_mapped_TF_paper_class_overlaps.tsv.gz"
  )
)
if (!nzchar(cfg$paper_class_qc_file)) {
  cfg$paper_class_qc_file <- file.path(cfg$paper_class_bed_dir, "paper_class_union_QC.tsv")
}

canonical_seqname <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chromosome[_.-]*", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  chromosome_number <- suppressWarnings(as.integer(x))
  out <- ifelse(!is.na(chromosome_number), paste0("chr", chromosome_number), paste0("chr", x))
  out[is.na(x) | !nzchar(x)] <- NA_character_
  out
}

parse_variant_coordinates <- function(x, id_column) {
  if (all(c("CHROM", "POS") %in% names(x))) {
    chromosome <- canonical_seqname(x$CHROM)
    position <- suppressWarnings(as.integer(x$POS))
  } else {
    variant_ids <- as.character(x[[id_column]])
    matched <- regexec("^([^:_]+)[:_]([0-9]+)(?:[:_].*)?$", variant_ids)
    parts <- regmatches(variant_ids, matched)
    chromosome <- canonical_seqname(vapply(
      parts, function(part) if (length(part)) part[[2L]] else NA_character_, character(1L)
    ))
    position <- suppressWarnings(as.integer(vapply(
      parts, function(part) if (length(part)) part[[3L]] else NA_character_, character(1L)
    )))
  }
  invalid <- is.na(chromosome) | is.na(position) | position < 1L
  if (any(invalid)) {
    stop_with(
      "Could not parse %d fine-mapped variant coordinate(s). Supply CHROM/POS columns or IDs formatted as <chromosome>_<position>.",
      sum(invalid)
    )
  }
  data.table(variant_seqname = chromosome, variant_position = position)
}

read_paper_class_ranges <- function(qc_file, bed_dir) {
  if (!file.exists(qc_file)) stop_with("Paper-class BED QC file does not exist: %s", qc_file)
  qc <- fread(qc_file, select = c("paper_class", "output_file"), showProgress = FALSE)
  qc[, `:=`(
    paper_class = trimws(as.character(paper_class)),
    output_file = trimws(as.character(output_file))
  )]
  if (!nrow(qc) || qc[, any(is.na(paper_class) | !nzchar(paper_class) | is.na(output_file) | !nzchar(output_file))]) {
    stop_with("Paper-class BED QC file has missing paper_class or output_file values: %s", qc_file)
  }
  if (qc[, anyDuplicated(paper_class) | anyDuplicated(output_file)]) {
    stop_with("Paper-class BED QC file contains duplicated classes or output files: %s", qc_file)
  }

  range_parts <- lapply(seq_len(nrow(qc)), function(i) {
    path <- file.path(bed_dir, qc$output_file[[i]])
    if (!file.exists(path)) stop_with("Paper-class BED listed in QC file does not exist: %s", path)
    bed <- fread(path, header = FALSE, select = 1:3, colClasses = "character", showProgress = FALSE)
    if (ncol(bed) < 3L) stop_with("Paper-class BED has fewer than three columns: %s", path)
    starts <- suppressWarnings(as.integer(bed[[2L]]))
    ends <- suppressWarnings(as.integer(bed[[3L]]))
    chromosomes <- canonical_seqname(bed[[1L]])
    valid <- !is.na(chromosomes) & !is.na(starts) & !is.na(ends) & starts >= 0L & ends > starts
    if (!all(valid)) stop_with("Paper-class BED has %d invalid row(s): %s", sum(!valid), path)
    GRanges(
      seqnames = chromosomes,
      ranges = IRanges(start = starts + 1L, end = ends),
      strand = "*",
      paper_class = qc$paper_class[[i]],
      source_file = basename(path)
    )
  })
  do.call(c, range_parts)
}

if (!file.exists(cfg$fine_mapped_file)) stop_with("Fine-mapped input does not exist: %s", cfg$fine_mapped_file)
if (!dir.exists(cfg$paper_class_bed_dir)) stop_with("Paper-class BED directory does not exist: %s", cfg$paper_class_bed_dir)
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cfg$overlap_output_file), recursive = TRUE, showWarnings = FALSE)

message("Reading fine-mapped variants: ", cfg$fine_mapped_file)
fine_mapped <- fread(cfg$fine_mapped_file)
variant_id_column <- cfg$variant_id_column
if (!variant_id_column %in% names(fine_mapped)) {
  candidates <- c("SNP", "variant_id", "variant", "snp")
  variant_id_column <- candidates[candidates %in% names(fine_mapped)][1L]
}
if (is.na(variant_id_column) || !nzchar(variant_id_column)) {
  stop_with("Fine-mapped input has no variant ID column. Set VARIANT_ID_COLUMN.")
}
fine_mapped[, variant_id := as.character(get(variant_id_column))]
if (fine_mapped[, anyNA(variant_id) | any(!nzchar(variant_id))]) {
  stop_with("Fine-mapped input contains missing variant IDs in column %s.", variant_id_column)
}
coordinates <- parse_variant_coordinates(fine_mapped, variant_id_column)
fine_mapped[, `:=`(
  input_row = .I,
  variant_seqname = coordinates$variant_seqname,
  variant_position = coordinates$variant_position,
  variant_key = paste(variant_id, coordinates$variant_seqname, coordinates$variant_position, sep = "|")
)]

unique_variants <- unique(
  fine_mapped[, .(variant_key, variant_id, variant_seqname, variant_position)],
  by = "variant_key"
)
variant_ranges <- GRanges(
  seqnames = unique_variants$variant_seqname,
  ranges = IRanges(start = unique_variants$variant_position, width = 1L),
  strand = "*"
)
names(variant_ranges) <- unique_variants$variant_key

message("Reading paper-class union BEDs from ", cfg$paper_class_bed_dir)
paper_class_ranges <- read_paper_class_ranges(cfg$paper_class_qc_file, cfg$paper_class_bed_dir)
hits <- findOverlaps(variant_ranges, paper_class_ranges, ignore.strand = TRUE)
overlaps <- if (length(hits)) {
  unique(data.table(
    variant_key = names(variant_ranges)[queryHits(hits)],
    paper_class = as.character(mcols(paper_class_ranges)$paper_class[subjectHits(hits)]),
    source_file = as.character(mcols(paper_class_ranges)$source_file[subjectHits(hits)]),
    peak_start = start(paper_class_ranges)[subjectHits(hits)] - 1L,
    peak_end = end(paper_class_ranges)[subjectHits(hits)]
  ), by = c("variant_key", "paper_class", "peak_start", "peak_end"))
} else {
  data.table(
    variant_key = character(), paper_class = character(), source_file = character(),
    peak_start = integer(), peak_end = integer()
  )
}

annotation <- overlaps[, .(
  in_TF_paper_class_peak = TRUE,
  n_TF_paper_class_peak_overlaps = .N,
  TF_paper_classes = paste(sort(unique(paper_class)), collapse = ";")
), by = variant_key]
annotation <- merge(unique_variants[, .(variant_key)], annotation, by = "variant_key", all.x = TRUE, sort = FALSE)
annotation[is.na(in_TF_paper_class_peak), `:=`(
  in_TF_paper_class_peak = FALSE,
  n_TF_paper_class_peak_overlaps = 0L,
  TF_paper_classes = NA_character_
)]

fine_mapped <- merge(fine_mapped, annotation, by = "variant_key", all.x = TRUE, sort = FALSE)
setorder(fine_mapped, input_row)
fine_mapped[, c("input_row", "variant_key") := NULL]
fwrite(fine_mapped, cfg$output_file, sep = "\t")
fwrite(overlaps, cfg$overlap_output_file, sep = "\t")

message(
  "Annotated ", nrow(fine_mapped), " fine-mapped rows (", length(variant_ranges), " unique variants); ",
  sum(annotation$in_TF_paper_class_peak), " unique variants overlap one or more paper-class peak sets.\n",
  "Wrote annotation: ", cfg$output_file, "\n",
  "Wrote interval overlaps: ", cfg$overlap_output_file
)
