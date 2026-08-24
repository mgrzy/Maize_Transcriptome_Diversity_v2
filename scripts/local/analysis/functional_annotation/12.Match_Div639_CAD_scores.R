#!/usr/bin/env Rscript
# Purpose: Compare CAD ZeroShotScore values between fine-mapped variants and a
#          one-to-one, non-fine-mapped Div639 control population matched for
#          target-gene TSS distance and MAF.
#
# Each unique fine-mapped variant is assigned to its closest target-gene TSS
# (when it is linked to more than one gene), then receives one unique random
# control. Controls are drawn from Div639 variants with a ZeroShotScore, on the
# same chromosome and within the same MAF and absolute TSS-distance bins. If a
# strict bin has no usable control, the script expands by one then two adjacent
# bins and records this in `match_tier`.
#
# Inputs:
#   - results/summary/AF/aFCn_MAF_se_filtered.tsv.gz
#   - data/processed/geno/plink/plink.frq
#   - results/summary/Cad/Div639_zeroShotScore.tsv.gz
#   - ../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
# Output:
#   - results/summary/Cad/fine_mapped_vs_matched_random_zeroShotScore.tsv.gz
#
# Usage:
#   Rscript scripts/local/analysis/12.Match_Div639_CAD_scores.R
#
# Optional settings:
#   SEED=20260820 MAF_BIN_WIDTH=0.05 CIS_WINDOW_BP=2000000 \
#     Rscript scripts/local/analysis/12.Match_Div639_CAD_scores.R

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicFeatures)
  library(txdbmaker)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

get_env_number <- function(name, default, minimum = 0) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || !is.finite(value) || value < minimum) {
    stop_with("%s must be a number >= %s.", name, minimum)
  }
  value
}

cfg <- list(
  fine_mapped_file = Sys.getenv(
    "FINE_MAPPED_FILE",
    "results/summary/AF/aFCn_MAF_se_filtered.tsv.gz"
  ),
  plink_freq_file = Sys.getenv(
    "PLINK_FREQ_FILE",
    "data/processed/geno/plink/plink.frq"
  ),
  cad_score_file = Sys.getenv(
    "CAD_SCORE_FILE",
    "results/summary/Cad/Div639_zeroShotScore.tsv.gz"
  ),
  gff_file = Sys.getenv(
    "GFF_PATH",
    "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
  ),
  output_file = Sys.getenv(
    "OUTPUT_FILE",
    "results/summary/Cad/fine_mapped_vs_matched_random_zeroShotScore.tsv.gz"
  ),
  cis_window_bp = as.integer(get_env_number("CIS_WINDOW_BP", 2000000, 1)),
  maf_bin_width = get_env_number("MAF_BIN_WIDTH", 0.05, .Machine$double.eps),
  seed = as.integer(get_env_number("SEED", 20260820, 0))
)
if (cfg$maf_bin_width > 0.5) stop_with("MAF_BIN_WIDTH must be <= 0.5.")

for (path in unlist(cfg[c("fine_mapped_file", "plink_freq_file", "cad_score_file", "gff_file")])) {
  if (!file.exists(path)) stop_with("Required input does not exist: %s", path)
}
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)

canonical_seqname <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chromosome[_.-]*", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  numeric <- suppressWarnings(as.integer(x))
  out <- ifelse(!is.na(numeric), paste0("chr", numeric), paste0("chr", x))
  out[is.na(x) | !nzchar(x)] <- NA_character_
  out
}

parse_variant_coordinates <- function(ids) {
  ids <- as.character(ids)
  out <- data.table(
    variant_seqname = canonical_seqname(sub("_[^_]+$", "", ids)),
    variant_position = suppressWarnings(as.integer(sub("^.*_", "", ids)))
  )
  if (out[, anyNA(variant_seqname) || anyNA(variant_position) || any(variant_position < 1L)]) {
    stop_with("Variant IDs must have the format chr<chromosome>_<position>.")
  }
  out
}

read_gene_tss <- function(gff_file) {
  message("Reading gene TSS coordinates from ", gff_file)
  txdb <- txdbmaker::makeTxDbFromGFF(gff_file)
  gene_ranges <- GenomicFeatures::genes(txdb)
  gene_ids <- names(gene_ranges)
  if (is.null(gene_ids) || anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    stop_with("Could not recover gene IDs from the GFF-derived TxDb.")
  }
  tss <- data.table(
    gene_id = as.character(gene_ids),
    variant_seqname = canonical_seqname(as.character(GenomicRanges::seqnames(gene_ranges))),
    gene_tss = as.integer(ifelse(
      as.character(GenomicRanges::strand(gene_ranges)) == "-",
      GenomicRanges::end(gene_ranges),
      GenomicRanges::start(gene_ranges)
    ))
  )
  if (tss[, anyDuplicated(gene_id)]) stop_with("The GFF contains duplicate gene IDs.")
  tss
}

distance_breaks <- c(
  0L, 1000L, 2000L, 5000L, 10000L, 25000L,
  50000L, 100000L, 250000L, 500000L, 1000000L, Inf
)
maf_breaks <- unique(c(seq(0, 0.5, by = cfg$maf_bin_width), 0.5))

add_matching_bins <- function(x) {
  x[, maf_bin := as.integer(cut(MAF, maf_breaks, include.lowest = TRUE, labels = FALSE))]
  x[, distance_bin := as.integer(cut(
    abs_distance_TSS,
    distance_breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = FALSE
  ))]
  x
}

add_maf_bin <- function(x) {
  x[, maf_bin := as.integer(cut(MAF, maf_breaks, include.lowest = TRUE, labels = FALSE))]
  x
}

message("Reading CAD ZeroShotScore annotations...")
cad <- fread(cfg$cad_score_file, select = c("variant_id", "zeroShotScore"), showProgress = TRUE)
cad[, variant_id := as.character(variant_id)]
cad[, zeroShotScore := as.numeric(zeroShotScore)]
cad <- cad[!is.na(variant_id) & nzchar(variant_id) & is.finite(zeroShotScore)]
if (cad[, anyDuplicated(variant_id)]) stop_with("CAD score table has duplicate variant IDs.")

message("Reading unfiltered Div639 MAF values...")
plink <- fread(cfg$plink_freq_file, select = c("SNP", "MAF"), showProgress = TRUE)
setnames(plink, "SNP", "variant_id")
plink[, `:=`(variant_id = as.character(variant_id), MAF = as.numeric(MAF))]
plink <- plink[is.finite(MAF) & MAF > 0 & MAF <= 0.5]
if (plink[, anyDuplicated(variant_id)]) stop_with("PLINK frequency table has duplicate variant IDs.")

setkey(cad, variant_id)
setkey(plink, variant_id)
background <- cad[plink, on = "variant_id", nomatch = 0L,
                  .(variant_id = x.variant_id, zeroShotScore = x.zeroShotScore, MAF = i.MAF)]
rm(cad, plink)
gc(verbose = FALSE)
message("Background contains ", format(nrow(background), big.mark = ","), " Div639 variants with CAD scores.")

background_coordinates <- parse_variant_coordinates(background$variant_id)
background[, `:=`(
  variant_seqname = background_coordinates$variant_seqname,
  variant_position = background_coordinates$variant_position
)]
# The background is genome-wide and has no target-gene context yet, so only
# assign its MAF stratum here. TSS-distance strata are calculated per signal.
background <- add_maf_bin(background)

message("Reading fine-mapped variant-gene pairs...")
fine <- fread(cfg$fine_mapped_file, select = c("gene_id", "variant_id"), showProgress = TRUE)
fine[, `:=`(gene_id = as.character(gene_id), variant_id = as.character(variant_id))]
fine <- unique(fine[!is.na(gene_id) & nzchar(gene_id) & !is.na(variant_id) & nzchar(variant_id)])
fine_coordinates <- parse_variant_coordinates(fine$variant_id)
fine[, `:=`(
  variant_seqname = fine_coordinates$variant_seqname,
  variant_position = fine_coordinates$variant_position
)]

gene_tss <- read_gene_tss(cfg$gff_file)
fine <- merge(fine, gene_tss, by = c("gene_id", "variant_seqname"), all.x = TRUE, sort = FALSE)
fine[, abs_distance_TSS := abs(variant_position - gene_tss)]
fine <- fine[!is.na(gene_tss) & abs_distance_TSS <= cfg$cis_window_bp]

fine <- merge(
  fine,
  background[, .(variant_id, MAF, zeroShotScore)],
  by = "variant_id",
  all = FALSE,
  sort = FALSE
)
fine <- add_matching_bins(fine)

# A SNP can be linked to several eGenes. Retain its closest valid target TSS,
# so the matched random population has the same number of unique variants.
setorder(fine, variant_id, abs_distance_TSS, gene_id)
fine <- fine[!duplicated(variant_id)]
if (!nrow(fine)) stop_with("No fine-mapped variants have a valid TSS context and CAD score.")
fine[, signal_id := .I]
message("Matching ", format(nrow(fine), big.mark = ","), " unique fine-mapped variants.")

# Keep compact, sorted per-MAF/chromosome indexes. This avoids materialising
# every SNP in a broad TSS-distance window for every fine-mapped signal.
setorder(background, variant_seqname, maf_bin, variant_position)
background[, group_id__ := paste(variant_seqname, maf_bin, sep = "\r")]
group_indices <- split(seq_len(nrow(background)), background$group_id__)
group_positions <- lapply(group_indices, function(index) background$variant_position[index])

# Fine-mapped IDs and already chosen controls are excluded. The environment
# enables fast checks while preserving a random control population without
# repeated variants.
unavailable_ids <- new.env(hash = TRUE, parent = emptyenv())
for (id in fine$variant_id) assign(id, TRUE, envir = unavailable_ids)
is_available <- function(id) !exists(id, envir = unavailable_ids, inherits = FALSE)

make_candidate_pieces <- function(signal, maf_expand, distance_expand) {
  maf_bins <- seq.int(
    max(1L, signal$maf_bin - maf_expand),
    min(length(maf_breaks) - 1L, signal$maf_bin + maf_expand)
  )
  distance_lower <- distance_breaks[max(1L, signal$distance_bin - distance_expand)]
  distance_upper <- distance_breaks[min(
    length(distance_breaks), signal$distance_bin + distance_expand + 1L
  )] - 1L
  if (!is.finite(distance_upper)) distance_upper <- cfg$cis_window_bp
  distance_upper <- min(distance_upper, cfg$cis_window_bp)

  # Distance zero (the TSS itself) is represented in the downstream interval
  # only, preventing it from being sampled twice when lower == 0.
  upstream_start <- max(1L, signal$gene_tss - distance_upper)
  upstream_end <- signal$gene_tss - distance_lower
  if (distance_lower == 0L) upstream_end <- upstream_end - 1L
  intervals <- rbind(
    c(upstream_start, upstream_end),
    c(signal$gene_tss + distance_lower, signal$gene_tss + distance_upper)
  )
  pieces <- list()
  for (maf_bin in maf_bins) {
    group_id <- paste(signal$variant_seqname, maf_bin, sep = "\r")
    index <- group_indices[[group_id]]
    positions <- group_positions[[group_id]]
    if (is.null(index) || !length(index)) next
    for (interval_index in seq_len(nrow(intervals))) {
      start <- intervals[interval_index, 1L]
      end <- intervals[interval_index, 2L]
      first <- findInterval(start - 1L, positions) + 1L
      last <- findInterval(end, positions)
      if (start <= end && first <= last) {
        pieces[[length(pieces) + 1L]] <- list(index = index, first = first, last = last)
      }
    }
  }
  pieces
}

sample_control <- function(signal) {
  tiers <- list(
    list(label = "exact_bins", maf_expand = 0L, distance_expand = 0L),
    list(label = "adjacent_maf_distance_bins", maf_expand = 1L, distance_expand = 1L),
    list(label = "expanded_maf_distance_bins", maf_expand = 2L, distance_expand = 2L)
  )
  for (tier in tiers) {
    pieces <- make_candidate_pieces(signal, tier$maf_expand, tier$distance_expand)
    if (!length(pieces)) next
    sizes <- vapply(pieces, function(piece) piece$last - piece$first + 1L, integer(1))
    n_candidates <- sum(sizes)
    cumulative_sizes <- cumsum(sizes)

    # Rejection sampling is efficient because fine-mapped and previously
    # selected variants form only a small fraction of the Div639 background.
    for (attempt in seq_len(min(1000L, max(100L, n_candidates)))) {
      selected_position <- sample.int(n_candidates, 1L)
      piece_index <- which(cumulative_sizes >= selected_position)[1L]
      offset <- selected_position - if (piece_index == 1L) 0L else cumulative_sizes[piece_index - 1L]
      row_index <- pieces[[piece_index]]$index[pieces[[piece_index]]$first + offset - 1L]
      candidate_id <- background$variant_id[row_index]
      if (is_available(candidate_id)) {
        assign(candidate_id, TRUE, envir = unavailable_ids)
        return(list(
          row_index = row_index,
          match_tier = tier$label,
          n_candidate_positions = n_candidates
        ))
      }
    }

    # This exhaustive fallback is rarely needed, but guarantees an available
    # candidate is found even in a small, heavily used matching stratum.
    candidate_rows <- unlist(lapply(pieces, function(piece) {
      piece$index[seq.int(piece$first, piece$last)]
    }), use.names = FALSE)
    available_rows <- candidate_rows[vapply(
      background$variant_id[candidate_rows], is_available, logical(1)
    )]
    if (length(available_rows)) {
      row_index <- sample(available_rows, 1L)
      assign(background$variant_id[row_index], TRUE, envir = unavailable_ids)
      return(list(
        row_index = row_index,
        match_tier = tier$label,
        n_candidate_positions = n_candidates
      ))
    }
  }
  NULL
}

set.seed(cfg$seed)
controls <- vector("list", nrow(fine))
for (i in seq_len(nrow(fine))) {
  controls[[i]] <- sample_control(fine[i])
  if (i %% 1000L == 0L || i == nrow(fine)) {
    message("  matched ", i, "/", nrow(fine), " fine-mapped variants")
  }
}
matched <- !vapply(controls, is.null, logical(1))
if (!all(matched)) {
  warning(sum(!matched), " fine-mapped variants had no eligible unique random control and were excluded.")
  fine <- fine[matched]
  controls <- controls[matched]
}
if (!nrow(fine)) stop_with("No fine-mapped variants could be matched to a random control.")

control_rows <- vapply(controls, `[[`, integer(1), "row_index")
control_tier <- vapply(controls, `[[`, character(1), "match_tier")
control_n_candidates <- vapply(controls, `[[`, integer(1), "n_candidate_positions")
control <- background[control_rows]

fine_out <- fine[, .(
  signal_id,
  target_gene = gene_id,
  matched_fine_mapped_variant_id = variant_id,
  variant_id,
  variant_set = "fine_mapped",
  is_fine_mapped = TRUE,
  chr = variant_seqname,
  pos = variant_position,
  MAF,
  gene_tss,
  abs_distance_TSS,
  maf_bin,
  distance_bin,
  zeroShotScore,
  match_tier = "observed",
  n_candidate_positions = NA_integer_
)]
random_out <- data.table(
  signal_id = fine$signal_id,
  target_gene = fine$gene_id,
  matched_fine_mapped_variant_id = fine$variant_id,
  variant_id = control$variant_id,
  variant_set = "random",
  is_fine_mapped = FALSE,
  chr = control$variant_seqname,
  pos = control$variant_position,
  MAF = control$MAF,
  gene_tss = fine$gene_tss,
  abs_distance_TSS = abs(control$variant_position - fine$gene_tss),
  maf_bin = control$maf_bin,
  distance_bin = as.integer(cut(
    abs(control$variant_position - fine$gene_tss),
    distance_breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = FALSE
  )),
  zeroShotScore = control$zeroShotScore,
  match_tier = control_tier,
  n_candidate_positions = control_n_candidates
)

output <- rbindlist(list(fine_out, random_out), use.names = TRUE)
setorder(output, signal_id, -is_fine_mapped)
fwrite(output, cfg$output_file, sep = "\t", compress = "gzip")

message(
  "Wrote ", nrow(output), " rows (", nrow(fine_out), " fine-mapped and ",
  nrow(random_out), " matched random variants) to ", cfg$output_file
)
