#!/usr/bin/env Rscript
# Test enrichment of fine-mapped cis-eQTL variants in replicate-supported
# Tu/Mejia-Guerra 2020 TF ChIP-seq consensus peaks (B73 v5). Enrichment is
# evaluated separately for every TF relative to TSS-distance-, MAF-, and
# variant-type-matched SNPs. One control is sampled independently from every
# eligible fine-mapped SNP's candidate pool in each permutation.
#
# By default this uses one SuSiE lead variant per credible set. Override input
# paths or permutation settings with environment variables, for example:
#   TF_PEAKS_DIR=data/processed/chipseq/Tu_2020_TF_consensus_v5/per_tf \
#   N_PERMUTATIONS=5000 N_CORES=8 \
#     Rscript scripts/local/analysis/17.Fine_Mapped_TF_ChIP_Enrichment.R
#
# The TF_PEAKS_DIR must contain one <TF>_consensus.bed[.gz] file per TF. These
# are replicate-supported consensus peaks, not formal IDR peaks. BED starts are
# interpreted as zero-based and BED ends as exclusive. Do not point this script
# at the parent directory containing all_TF or any_TF aggregate BED files.
#
# Outputs:
#   results/summary/fine_mapped_TF_ChIP_enrichment/
#     fine_mapped_variant_TF_annotations.tsv.gz
#     fine_mapped_variant_TF_overlaps.tsv.gz
#     fine_mapped_variants_excluded_from_matched_controls.tsv.gz
#     matched_control_variants.tsv.gz
#     matched_control_variant_TF_annotations.tsv.gz
#     TF_overlap_matched_permutation.tsv
#     TF_overlap_matched_permutation_null_counts.rds

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(IRanges)
  library(rtracklayer)
  library(txdbmaker)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

get_env_integer <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop_with("%s must be an integer >= %d.", name, minimum)
  }
  value
}

get_env_logical <- function(name, default) {
  value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false")) stop_with("%s must be true or false.", name)
  identical(value, "true")
}

cfg <- list(
  chromatin_dir = Sys.getenv(
    "TF_PEAKS_DIR",
    "data/processed/chipseq/Tu_2020_TF_consensus_v5/per_tf"
  ),
  fine_mapped_file = Sys.getenv(
    "FINE_MAPPED_FILE",
    "results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz"
  ),
  plink_freq_file = Sys.getenv("PLINK_FREQ_FILE", "data/processed/geno/plink/plink.frq"),
  # Either provide a two-column chromosome-size file through GENOME_SIZES_PATH
  # (chromosome, length), or leave it empty to obtain chromosome lengths from
  # the GFF chromosome records below.
  genome_sizes_file = Sys.getenv("GENOME_SIZES_PATH", ""),
  gff_file = Sys.getenv(
    "GFF_PATH",
    "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
  ),
  output_dir = Sys.getenv(
    "OUTPUT_DIR",
    "results/summary/fine_mapped_TF_ChIP_enrichment"
  ),
  # These retain TF-specific defaults while allowing a wrapper to reuse the
  # matched-control implementation for another labelled BED feature set.
  feature_output_stem = Sys.getenv("FEATURE_OUTPUT_STEM", "TF"),
  feature_result_column = Sys.getenv("FEATURE_RESULT_COLUMN", "TF"),
  feature_label_map_file = Sys.getenv("FEATURE_LABEL_MAP_FILE", ""),
  variant_id_column = Sys.getenv("VARIANT_ID_COLUMN", "SNP"),
  target_gene_column = Sys.getenv("TARGET_GENE_COLUMN", "gene"),
  # Use +/- 2 Mb around target TSSs, matching the genomic-context enrichment
  # background and retaining all lead variants from the current SuSiE results.
  cis_window_bp = get_env_integer("CIS_WINDOW_BP", 2000000L),
  maf_bin_width = as.numeric(Sys.getenv("MAF_BIN_WIDTH", "0.05")),
  distance_breaks = c(0L, 1000L, 2000L, 5000L, 10000L, 25000L,
                      50000L, 100000L, 250000L, 500000L, 1000000L, Inf),
  max_controls_per_stratum = get_env_integer("MAX_CONTROLS_PER_STRATUM", 500L),
  n_permutations = get_env_integer("N_PERMUTATIONS", 5000L),
  seed = get_env_integer("SEED", 20260817L, minimum = 0L),
  use_multicore = get_env_logical("USE_MULTICORE", TRUE),
  n_cores = get_env_integer("N_CORES", 1L),
  save_null_counts = get_env_logical("SAVE_NULL_COUNTS", TRUE),
  expected_n_tfs = get_env_integer("EXPECTED_N_TFS", 104L),
  save_overlap_annotations = get_env_logical("SAVE_OVERLAP_ANNOTATIONS", TRUE)
)
if (!is.finite(cfg$maf_bin_width) || cfg$maf_bin_width <= 0 || cfg$maf_bin_width > 0.5) {
  stop_with("MAF_BIN_WIDTH must be > 0 and <= 0.5.")
}

# Map common chromosome spellings (1, Chr01, chromosome_1) to the convention
# used by fine-mapped SNP IDs (chr1). Non-numeric contigs are retained with a
# chr prefix, so an incompatible assembly fails explicitly below.
canonical_seqname <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chromosome[_.-]*", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  numeric <- suppressWarnings(as.integer(x))
  out <- ifelse(!is.na(numeric), paste0("chr", numeric), paste0("chr", x))
  out[!nzchar(x) | is.na(x)] <- NA_character_
  out
}

file_label <- function(path, root) {
  label <- sub(paste0("^", root, "/?"), "", path)
  label <- sub("\\.(gz|bgz)$", "", label, ignore.case = TRUE)
  label <- sub("\\.bed$", "", label, ignore.case = TRUE)
  label <- sub("^paper_class_", "", label, ignore.case = TRUE)
  sub("(_consensus|_union)$", "", label, ignore.case = TRUE)
}

parse_variant_coordinates <- function(fine_mapped, id_column) {
  if (all(c("CHROM", "POS") %in% names(fine_mapped))) {
    seqname <- canonical_seqname(fine_mapped$CHROM)
    position <- suppressWarnings(as.integer(fine_mapped$POS))
  } else {
    variant_id <- as.character(fine_mapped[[id_column]])
    seqname <- canonical_seqname(sub("_[^_]+$", "", variant_id))
    position <- suppressWarnings(as.integer(sub("^.*_", "", variant_id)))
  }
  invalid <- is.na(seqname) | is.na(position) | position < 1L
  if (any(invalid)) {
    stop_with(
      "Could not determine valid 1-based chromosome/position coordinates for %d fine-mapped rows. ",
      sum(invalid),
      "Supply CHROM and POS columns, or use SNP IDs formatted as <chromosome>_<position>."
    )
  }
  data.table(variant_seqname = seqname, variant_position = position)
}

read_bed_like <- function(path, label) {
  peak_table <- fread(
    path,
    header = FALSE,
    fill = TRUE,
    colClasses = "character",
    showProgress = FALSE
  )
  if (ncol(peak_table) < 3L) {
    stop_with("BED-like chromatin file has fewer than three columns: %s", path)
  }
  starts <- suppressWarnings(as.integer(peak_table[[2L]]))
  ends <- suppressWarnings(as.integer(peak_table[[3L]]))
  keep <- !is.na(starts) & !is.na(ends) & starts >= 0L & ends > starts
  if (!any(keep)) stop_with("No valid BED intervals found in: %s", path)
  GRanges(
    seqnames = canonical_seqname(peak_table[[1L]][keep]),
    ranges = IRanges(start = starts[keep] + 1L, end = ends[keep]),
    strand = "*",
    chromatin_feature = label,
    source_file = basename(path)
  )
}

read_chromatin_file <- function(path, label) {
  extension <- tolower(sub(".*\\.", "", sub("\\.(gz|bgz)$", "", path, ignore.case = TRUE)))
  if (extension %in% c("gff3", "gff", "gtf")) {
    peaks <- rtracklayer::import(path)
    peaks <- peaks[width(peaks) > 0L]
    seqnames(peaks) <- canonical_seqname(seqnames(peaks))
    mcols(peaks)$chromatin_feature <- label
    mcols(peaks)$source_file <- basename(path)
    return(peaks)
  }
  read_bed_like(path, label)
}

read_chromatin_features <- function(directory) {
  if (!dir.exists(directory)) stop_with("Chromatin directory does not exist: %s", directory)
  files <- list.files(
    directory,
    pattern = "(_consensus|_union)\\.bed(\\.gz)?$",
    recursive = FALSE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(files)) {
    stop_with("No *_consensus.bed[.gz] or *_union.bed[.gz] files were found directly under: %s", directory)
  }
  labels <- make.unique(vapply(files, file_label, character(1), root = directory))
  peaks <- Map(read_chromatin_file, files, labels)
  empty <- !vapply(peaks, length, integer(1))
  if (any(empty)) {
    warning("Skipping empty chromatin files: ", paste(basename(files[empty]), collapse = ", "))
    peaks <- peaks[!empty]
    labels <- labels[!empty]
  }
  names(peaks) <- labels
  if (anyDuplicated(names(peaks))) {
    stop_with("Chromatin-feature names are duplicated after removing BED filename suffixes.")
  }
  peaks
}

apply_feature_label_map <- function(features, path, feature_directory) {
  if (!nzchar(path)) return(features)
  if (!file.exists(path)) stop_with("FEATURE_LABEL_MAP_FILE does not exist: %s", path)
  label_map <- fread(path, select = c("paper_class", "output_file"), showProgress = FALSE)
  label_map[, `:=`(
    paper_class = trimws(as.character(paper_class)),
    output_file = trimws(as.character(output_file)),
    feature_label = file_label(output_file, feature_directory)
  )]
  if (
    !nrow(label_map) ||
    label_map[, any(is.na(paper_class) | !nzchar(paper_class) | is.na(feature_label) | !nzchar(feature_label))] ||
    label_map[, anyDuplicated(feature_label)] || label_map[, anyDuplicated(paper_class)]
  ) {
    stop_with("FEATURE_LABEL_MAP_FILE must contain unique, nonempty paper_class and output_file values: %s", path)
  }
  label_lookup <- setNames(label_map$paper_class, label_map$feature_label)
  missing_labels <- setdiff(names(features), names(label_lookup))
  if (length(missing_labels)) {
    stop_with(
      "FEATURE_LABEL_MAP_FILE lacks label(s) for BED feature(s): %s",
      paste(missing_labels, collapse = ", ")
    )
  }
  names(features) <- unname(label_lookup[names(features)])
  if (anyDuplicated(names(features))) stop_with("FEATURE_LABEL_MAP_FILE maps multiple BED features to the same label.")
  features
}

read_genome <- function() {
  if (nzchar(cfg$genome_sizes_file)) {
    if (!file.exists(cfg$genome_sizes_file)) {
      stop_with("GENOME_SIZES_PATH does not exist: %s", cfg$genome_sizes_file)
    }
    sizes <- fread(cfg$genome_sizes_file, header = FALSE, select = 1:2, colClasses = "character")
    setnames(sizes, c("seqname", "length"))
    sizes[, `:=`(
      seqname = canonical_seqname(seqname),
      length = suppressWarnings(as.integer(length))
    )]
  } else {
    if (!file.exists(cfg$gff_file)) {
      stop_with(
        "Neither GENOME_SIZES_PATH nor a readable GFF chromosome file is available. Set GENOME_SIZES_PATH, or provide GFF_PATH (missing: %s).",
        cfg$gff_file
      )
    }
    chromosome_ranges <- rtracklayer::import(
      cfg$gff_file,
      format = "gff3",
      feature.type = "chromosome"
    )
    sizes <- data.table(
      seqname = canonical_seqname(seqnames(chromosome_ranges)),
      length = as.integer(end(chromosome_ranges))
    )
  }
  sizes <- unique(sizes[!is.na(seqname) & !is.na(length) & length > 0L], by = "seqname")
  if (!nrow(sizes)) stop_with("No valid chromosome lengths were available for permutations.")
  GRanges(seqnames = sizes$seqname, ranges = IRanges(start = 1L, end = sizes$length), strand = "*")
}

read_gene_tss <- function(gff_file) {
  if (!file.exists(gff_file)) stop_with("GFF file does not exist: %s", gff_file)
  txdb <- txdbmaker::makeTxDbFromGFF(gff_file)
  gene_ranges <- GenomicFeatures::genes(txdb)
  gene_ids <- names(gene_ranges)
  if (is.null(gene_ids) || anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    stop_with("Could not recover gene IDs from the GFF-derived TxDb.")
  }
  gene_tss <- ifelse(as.character(strand(gene_ranges)) == "-", end(gene_ranges), start(gene_ranges))
  tss <- data.table(
    target_gene = as.character(gene_ids),
    variant_seqname = canonical_seqname(seqnames(gene_ranges)),
    gene_tss = as.integer(gene_tss)
  )
  if (tss[, anyDuplicated(target_gene)]) {
    stop_with("The GFF-derived TxDb has duplicate gene identifiers.")
  }
  tss
}

make_variant_type <- function(a1, a2) {
  fifelse(nchar(a1) == 1L & nchar(a2) == 1L, "SNP", "INDEL")
}

add_matching_bins <- function(x) {
  maf_breaks <- unique(c(seq(0, 0.5, by = cfg$maf_bin_width), 0.5))
  x[, maf_bin := cut(MAF, maf_breaks, include.lowest = TRUE, labels = FALSE)]
  x[, distance_bin := cut(
    abs_distance_TSS,
    cfg$distance_breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = FALSE
  )]
  if (x[, anyNA(maf_bin) || anyNA(distance_bin)]) {
    stop_with("MAF or target-TSS distance falls outside the configured matching bins.")
  }
  invisible(x)
}

# Build upstream and downstream intervals for each distance stratum. These are
# clipped both to chromosome ends and to the configured target-gene cis window.
make_stratum_intervals <- function(strata, chromosome_ends) {
  x <- merge(strata, chromosome_ends, by = "variant_seqname", all.x = TRUE, sort = FALSE)
  if (x[, anyNA(chromosome_end)]) {
    stop_with("Some target genes have no chromosome length for matched-control sampling.")
  }
  lower <- cfg$distance_breaks[x$distance_bin]
  upper <- pmin(cfg$distance_breaks[x$distance_bin + 1L] - 1L, cfg$cis_window_bp)
  upper[!is.finite(upper)] <- cfg$cis_window_bp
  downstream <- x[, .(
    stratum_id, variant_seqname, maf_bin, variant_type,
    interval_start = pmax(1L, gene_tss + lower),
    interval_end = pmin(chromosome_end, gene_tss + upper)
  )]
  upstream <- x[, .(
    stratum_id, variant_seqname, maf_bin, variant_type,
    interval_start = pmax(1L, gene_tss - upper),
    interval_end = pmin(chromosome_end, gene_tss - lower)
  )]
  unique(rbindlist(list(downstream, upstream))[interval_start <= interval_end])
}

read_plink_background <- function(path) {
  if (!file.exists(path)) stop_with("PLINK frequency table does not exist: %s", path)
  message("Reading PLINK frequency background (MAF and variant type): ", path)
  plink <- fread(path, select = c("SNP", "A1", "A2", "MAF"), showProgress = TRUE)
  setnames(plink, "SNP", "variant_id")
  plink[, variant_id := as.character(variant_id)]
  coordinates <- parse_variant_coordinates(plink, "variant_id")
  plink[, `:=`(
    variant_seqname = coordinates$variant_seqname,
    variant_position = coordinates$variant_position,
    MAF = as.numeric(MAF),
    variant_type = make_variant_type(as.character(A1), as.character(A2))
  )]
  plink <- plink[is.finite(MAF) & MAF > 0 & MAF <= 0.5]
  if (plink[, anyDuplicated(variant_id)]) {
    stop_with("PLINK frequency table has duplicate variant IDs.")
  }
  plink
}

sample_matched_controls <- function(fine_mapped, genome, target_gene_column) {
  message("Deriving target-gene TSS contexts for matched negative controls...")
  gene_tss <- read_gene_tss(cfg$gff_file)
  lead_context <- unique(fine_mapped[, .(
    variant_key, variant_id, target_gene = as.character(get(target_gene_column)),
    variant_seqname, variant_position
  )])
  lead_context <- merge(
    lead_context,
    gene_tss,
    by = c("target_gene", "variant_seqname"),
    all.x = TRUE,
    sort = FALSE
  )
  lead_context[, abs_distance_TSS := abs(variant_position - gene_tss)]
  lead_context[, exclusion_reason := fifelse(
    is.na(gene_tss),
    "no_matching_target_gene_TSS_on_variant_chromosome",
    fifelse(
      abs_distance_TSS > cfg$cis_window_bp,
      paste0("outside_", cfg$cis_window_bp, "bp_target_gene_cis_window"),
      NA_character_
    )
  )]
  valid_context <- lead_context[is.na(exclusion_reason)]
  excluded_variant_keys <- setdiff(unique(fine_mapped$variant_key), unique(valid_context$variant_key))
  excluded_context <- lead_context[variant_key %in% excluded_variant_keys]
  if (length(excluded_variant_keys)) {
    # A SNP may be linked to several target genes in the fine-mapping input.
    # Keep one deterministic diagnostic row only when none of those links
    # provides a valid same-chromosome cis context.
    setorder(excluded_context, variant_key, exclusion_reason, target_gene)
    excluded_context <- excluded_context[, .SD[1L], by = variant_key]
    warning(
      length(excluded_variant_keys),
      " unique fine-mapped SNPs lack a valid target-gene cis context and will be excluded from matched-control permutations. ",
      "They remain in the fine-mapped overlap annotation output."
    )
  }
  lead_context <- valid_context
  if (!nrow(lead_context)) stop_with("No fine-mapped variants have a matching target-gene TSS within the cis window.")

  # The same variant can lead credible sets for multiple genes. The enrichment
  # test uses unique SNP positions, so retain its closest target-gene context
  # deterministically to produce exactly one matched control per observed SNP.
  setorder(lead_context, variant_key, abs_distance_TSS, target_gene)
  lead_context <- lead_context[, .SD[1L], by = variant_key]
  plink <- read_plink_background(cfg$plink_freq_file)
  lead_context <- merge(
    lead_context,
    plink[, .(variant_id, MAF, variant_type)],
    by = "variant_id",
    all.x = TRUE,
    sort = FALSE
  )
  if (lead_context[, anyNA(MAF) | anyNA(variant_type)]) {
    stop_with("Some fine-mapped SNPs lack MAF or allele-type data in the PLINK frequency table.")
  }
  add_matching_bins(lead_context)
  lead_context[, signal_id := .I]

  chromosome_ends <- data.table(
    variant_seqname = as.character(seqnames(genome)),
    chromosome_end = as.integer(end(genome))
  )
  stratum_columns <- c(
    "target_gene", "variant_seqname", "gene_tss", "maf_bin", "distance_bin", "variant_type"
  )
  lead_context[, stratum_id := .GRP, by = stratum_columns]
  strata <- unique(lead_context[, c("stratum_id", stratum_columns), with = FALSE])
  intervals <- make_stratum_intervals(strata, chromosome_ends)

  maf_breaks <- unique(c(seq(0, 0.5, by = cfg$maf_bin_width), 0.5))
  plink[, maf_bin := cut(MAF, maf_breaks, include.lowest = TRUE, labels = FALSE)]
  setkey(plink, variant_seqname, maf_bin, variant_type, variant_position)
  set.seed(cfg$seed + 1L)
  batch_size <- 500L
  control_parts <- vector("list", ceiling(nrow(intervals) / batch_size))
  message("Building TSS-distance/MAF/type matched control pools...")
  for (batch_index in seq_along(control_parts)) {
    first <- (batch_index - 1L) * batch_size + 1L
    last <- min(nrow(intervals), batch_index * batch_size)
    interval_batch <- intervals[first:last]
    hits <- plink[interval_batch,
      on = .(variant_seqname, maf_bin, variant_type,
             variant_position >= interval_start, variant_position <= interval_end),
      nomatch = 0L,
      allow.cartesian = TRUE,
      .(stratum_id = i.stratum_id, control_variant_id = x.variant_id,
        control_seqname = x.variant_seqname, control_position = x.variant_position,
        control_MAF = x.MAF, control_variant_type = x.variant_type)
    ]
    if (nrow(hits)) {
      hits <- unique(hits, by = c("stratum_id", "control_variant_id"))
      hits <- hits[, .SD[sample.int(.N, min(.N, cfg$max_controls_per_stratum))], by = stratum_id]
    }
    control_parts[[batch_index]] <- hits
    if (batch_index %% 10L == 0L || batch_index == length(control_parts)) {
      message("  processed ", last, "/", nrow(intervals), " distance intervals")
    }
  }
  control_candidates <- rbindlist(control_parts)
  rm(control_parts, intervals)

  lead_summary <- lead_context[, .(
    signal_id, stratum_id, target_gene, lead_variant_id = variant_id,
    lead_seqname = variant_seqname, lead_position = variant_position,
    lead_MAF = MAF, abs_distance_TSS, maf_bin, distance_bin, variant_type
  )]
  matched_pool <- merge(
    lead_summary,
    control_candidates,
    by = "stratum_id",
    allow.cartesian = TRUE,
    sort = FALSE
  )
  observed_ids <- unique(fine_mapped$variant_id)
  matched_pool <- matched_pool[
    control_variant_id != lead_variant_id & !control_variant_id %chin% observed_ids
  ]
  matched_pool[, match_tier := "exact_bins"]

  # A small fraction of signals, especially rare variants in short TSS-distance
  # bins, have no non-fine-mapped SNP in the exact stratum. Expand matching
  # locally and record the tier used, rather than silently dropping signals.
  unmatched_signal_ids <- function(pool) {
    matched_ids <- unique(pool$signal_id)
    setdiff(lead_context$signal_id, matched_ids)
  }
  make_relaxed_intervals <- function(leads, maf_expand, distance_expand) {
    n_maf_bins <- length(maf_breaks) - 1L
    n_distance_bins <- length(cfg$distance_breaks) - 1L
    expanded <- leads[, {
      candidate_maf_bins <- seq.int(
        max(1L, maf_bin - maf_expand),
        min(n_maf_bins, maf_bin + maf_expand)
      )
      .(maf_bin = candidate_maf_bins)
    }, by = .(signal_id, variant_seqname, gene_tss, distance_bin, variant_type)]
    expanded <- chromosome_ends[expanded, on = "variant_seqname", nomatch = 0L]
    lower_index <- pmax(1L, expanded$distance_bin - distance_expand)
    upper_index <- pmin(n_distance_bins + 1L, expanded$distance_bin + distance_expand + 1L)
    lower <- cfg$distance_breaks[lower_index]
    upper <- pmin(cfg$distance_breaks[upper_index] - 1L, cfg$cis_window_bp)
    upper[!is.finite(upper)] <- cfg$cis_window_bp
    downstream <- expanded[, .(
      signal_id, variant_seqname, maf_bin, variant_type,
      interval_start = pmax(1L, gene_tss + lower),
      interval_end = pmin(chromosome_end, gene_tss + upper)
    )]
    upstream <- expanded[, .(
      signal_id, variant_seqname, maf_bin, variant_type,
      interval_start = pmax(1L, gene_tss - upper),
      interval_end = pmin(chromosome_end, gene_tss - lower)
    )]
    unique(rbindlist(list(downstream, upstream))[interval_start <= interval_end])
  }
  collect_relaxed_candidates <- function(intervals) {
    if (!nrow(intervals)) return(data.table())
    batch_size <- 500L
    parts <- vector("list", ceiling(nrow(intervals) / batch_size))
    for (batch_index in seq_along(parts)) {
      first <- (batch_index - 1L) * batch_size + 1L
      last <- min(nrow(intervals), batch_index * batch_size)
      hits <- plink[intervals[first:last],
        on = .(variant_seqname, maf_bin, variant_type,
               variant_position >= interval_start, variant_position <= interval_end),
        nomatch = 0L,
        allow.cartesian = TRUE,
        .(signal_id = i.signal_id, control_variant_id = x.variant_id,
          control_seqname = x.variant_seqname, control_position = x.variant_position,
          control_MAF = x.MAF, control_variant_type = x.variant_type)
      ]
      if (nrow(hits)) {
        hits <- unique(hits, by = c("signal_id", "control_variant_id"))
        hits <- hits[, .SD[sample.int(.N, min(.N, cfg$max_controls_per_stratum))], by = signal_id]
      }
      parts[[batch_index]] <- hits
    }
    candidates <- unique(rbindlist(parts), by = c("signal_id", "control_variant_id"))
    candidates[!control_variant_id %chin% observed_ids]
  }
  fallback_tiers <- list(
    list(label = "adjacent_maf_distance_bins", maf_expand = 1L, distance_expand = 1L),
    list(label = "expanded_maf_distance_bins", maf_expand = 2L, distance_expand = 2L)
  )
  for (tier in fallback_tiers) {
    missing_ids <- unmatched_signal_ids(matched_pool)
    if (!length(missing_ids)) break
    message("  relaxing matching for ", length(missing_ids), " signals: ", tier$label)
    fallback_intervals <- make_relaxed_intervals(
      lead_context[signal_id %in% missing_ids],
      tier$maf_expand,
      tier$distance_expand
    )
    fallback_candidates <- collect_relaxed_candidates(fallback_intervals)
    if (nrow(fallback_candidates)) {
      fallback_candidates[, match_tier := tier$label]
      fallback_pool <- merge(
        lead_summary[signal_id %in% missing_ids],
        fallback_candidates,
        by = "signal_id",
        allow.cartesian = TRUE,
        sort = FALSE
      )
      matched_pool <- rbindlist(list(matched_pool, fallback_pool), use.names = TRUE, fill = TRUE)
    }
  }

  # A candidate can be recovered from both upstream/downstream interval pieces
  # or from interval batches. It must occur only once within a signal so every
  # eligible control has equal sampling probability.
  matched_pool <- unique(matched_pool, by = c("signal_id", "control_variant_id"))

  pool_sizes <- matched_pool[, .(n_matched_candidates = uniqueN(control_variant_id)), by = signal_id]
  lead_context <- pool_sizes[lead_context, on = "signal_id"]
  lead_context[is.na(n_matched_candidates), n_matched_candidates := 0L]
  no_control_context <- lead_context[n_matched_candidates == 0L]
  if (nrow(no_control_context)) {
    no_control_diagnostics <- no_control_context[, .(
      variant_key, variant_id, target_gene, variant_seqname, variant_position,
      gene_tss, abs_distance_TSS, n_matched_candidates,
      exclusion_reason = "no_eligible_matched_control_after_local_MAF_TSS_bin_expansion"
    )]
    warning(
      nrow(no_control_diagnostics),
      " fine-mapped SNPs have no eligible matched control after local MAF/TSS-bin expansion and will be excluded from matched-control permutations."
    )
    excluded_context <- rbindlist(
      list(excluded_context, no_control_diagnostics),
      use.names = TRUE,
      fill = TRUE
    )
    matched_pool <- matched_pool[!signal_id %in% no_control_context$signal_id]
    lead_context <- lead_context[n_matched_candidates > 0L]
  }
  if (!nrow(lead_context)) stop_with("No fine-mapped SNPs have eligible matched controls.")
  matched_controls <- matched_pool[, .SD[sample.int(.N, 1L)], by = signal_id]
  matched_controls <- merge(
    lead_context[, .(signal_id, variant_key, n_matched_candidates)],
    matched_controls,
    by = "signal_id",
    all.x = TRUE,
    sort = FALSE
  )
  if (nrow(matched_controls) != nrow(lead_context)) {
    stop_with("Matched control sampling did not return one control for every eligible fine-mapped SNP.")
  }
  setorder(excluded_context, variant_key, exclusion_reason)
  list(
    matched_controls = matched_controls,
    matched_pool = matched_pool,
    eligible_leads = lead_context,
    excluded_context = unique(excluded_context, by = "variant_key")
  )
}

make_overlap_table <- function(query_ranges, query_table) {
  overlap_parts <- Map(
    function(label, peaks) {
      hits <- findOverlaps(query_ranges, peaks, ignore.strand = TRUE)
      if (!length(hits)) return(NULL)
      query_index <- queryHits(hits)
      subject_index <- subjectHits(hits)
      data.table(
        query_key = names(query_ranges)[query_index],
        chromatin_feature = label,
        source_file = mcols(peaks)$source_file[subject_index],
        peak_start = start(peaks)[subject_index],
        peak_end = end(peaks)[subject_index]
      )
    },
    names(chromatin_features),
    chromatin_features
  )
  overlaps <- rbindlist(overlap_parts, use.names = TRUE, fill = TRUE)
  if (!nrow(overlaps)) {
    overlaps <- data.table(
      query_key = character(), chromatin_feature = character(), source_file = character(),
      peak_start = integer(), peak_end = integer()
    )
  }
  merge(query_table, overlaps, by = "query_key", all = FALSE, allow.cartesian = TRUE, sort = FALSE)
}

add_chromatin_annotation <- function(query_table, overlaps) {
  annotation <- overlaps[, .(
    overlaps_chromatin = TRUE,
    n_chromatin_peak_overlaps = .N,
    chromatin_features = paste(sort(unique(chromatin_feature)), collapse = ";")
  ), by = query_key]
  annotated <- merge(query_table, annotation, by = "query_key", all.x = TRUE, sort = FALSE)
  annotated[is.na(overlaps_chromatin), `:=`(
    overlaps_chromatin = FALSE,
    n_chromatin_peak_overlaps = 0L,
    chromatin_features = NA_character_
  )]
  annotated
}

make_feature_membership <- function(ranges, feature_list) {
  query_ids <- names(ranges)
  if (is.null(query_ids) || anyNA(query_ids) || any(!nzchar(query_ids))) {
    stop_with("Feature-membership ranges must have unique, nonempty names.")
  }
  if (anyDuplicated(query_ids)) {
    stop_with("Feature-membership ranges have duplicated names.")
  }
  membership <- lapply(feature_list, function(peaks) {
    hits <- findOverlaps(ranges, peaks, ignore.strand = TRUE)
    if (!length(hits)) return(character())
    unique(query_ids[queryHits(hits)])
  })
  names(membership) <- names(feature_list)
  membership
}

combine_membership <- function(membership, labels) {
  labels <- intersect(labels, names(membership))
  if (!length(labels)) return(character())
  unique(unlist(membership[labels], use.names = FALSE))
}

summarise_matched_null <- function(label, observed, null_counts, n_signals) {
  expected <- mean(null_counts)
  p_enrichment <- (1 + sum(null_counts >= observed)) / (length(null_counts) + 1)
  p_depletion <- (1 + sum(null_counts <= observed)) / (length(null_counts) + 1)
  data.table(
    chromatin_feature = label,
    observed_unique_variants = as.integer(observed),
    observed_fraction = observed / n_signals,
    expected_unique_variants = expected,
    expected_fraction = expected / n_signals,
    null_sd = stats::sd(null_counts),
    null_q025 = as.numeric(stats::quantile(null_counts, 0.025, names = FALSE)),
    null_q975 = as.numeric(stats::quantile(null_counts, 0.975, names = FALSE)),
    fold_enrichment = if (expected > 0) observed / expected else NA_real_,
    log2_fold_enrichment = if (observed > 0 && expected > 0) log2(observed / expected) else NA_real_,
    p_enrichment = p_enrichment,
    p_depletion = p_depletion,
    p_two_sided = min(1, 2 * min(p_enrichment, p_depletion)),
    n_permutations = length(null_counts),
    n_variants = n_signals
  )
}
if (!file.exists(cfg$fine_mapped_file)) {
  stop_with("Fine-mapped input does not exist: %s", cfg$fine_mapped_file)
}
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

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
target_gene_column <- cfg$target_gene_column
if (!target_gene_column %in% names(fine_mapped)) {
  candidates <- c("gene", "target_gene", "gene_id")
  target_gene_column <- candidates[candidates %in% names(fine_mapped)][1L]
}
if (is.na(target_gene_column) || !nzchar(target_gene_column)) {
  stop_with("Fine-mapped input has no target-gene column. Set TARGET_GENE_COLUMN.")
}
fine_mapped[, variant_id := as.character(get(variant_id_column))]
if (anyNA(fine_mapped$variant_id) || any(!nzchar(fine_mapped$variant_id))) {
  stop_with("Fine-mapped input has missing variant IDs in column %s.", variant_id_column)
}
if (anyNA(fine_mapped[[target_gene_column]]) || any(!nzchar(fine_mapped[[target_gene_column]]))) {
  stop_with("Fine-mapped input has missing target genes in column %s.", target_gene_column)
}
coordinates <- parse_variant_coordinates(fine_mapped, variant_id_column)
fine_mapped[, `:=`(
  variant_seqname = coordinates$variant_seqname,
  variant_position = coordinates$variant_position,
  input_row = .I,
  variant_key = paste(variant_id, coordinates$variant_seqname, coordinates$variant_position, sep = "|")
)]

# The annotation output preserves every fine-mapped row. Permutations instead
# use distinct genomic variant positions so repeat variant-gene associations do
# not inflate the observed overlap count.
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

message("Reading chromatin intervals from: ", cfg$chromatin_dir)
chromatin_features <- read_chromatin_features(cfg$chromatin_dir)
chromatin_features <- apply_feature_label_map(
  chromatin_features, cfg$feature_label_map_file, cfg$chromatin_dir
)
if (length(chromatin_features) != cfg$expected_n_tfs) {
  warning(
    "Detected ", length(chromatin_features), " chromatin feature files; expected ",
    cfg$expected_n_tfs, ". All detected TFs will still be analysed."
  )
}
message("Detected ", length(chromatin_features), " chromatin feature peak sets.")
tf_peak_qc <- rbindlist(Map(
  function(tf, peaks) {
    reduced <- reduce(peaks, ignore.strand = TRUE)
    data.table(
      chromatin_feature = tf,
      n_consensus_peaks = length(peaks),
      consensus_covered_bp = sum(width(reduced))
    )
  },
  names(chromatin_features), chromatin_features
))
peak_seqlevels <- unique(unlist(lapply(chromatin_features, function(x) as.character(seqnames(x)))))
if (!all(as.character(seqnames(variant_ranges)) %in% peak_seqlevels)) {
  absent <- setdiff(unique(as.character(seqnames(variant_ranges))), peak_seqlevels)
  warning("Fine-mapped variants on chromosome(s) absent from chromatin files: ",
          paste(absent, collapse = ", "))
}

fine_query <- copy(unique_variants)
setnames(fine_query, "variant_key", "query_key")
fine_overlap_table <- make_overlap_table(variant_ranges, fine_query)
fine_annotation <- add_chromatin_annotation(fine_query, fine_overlap_table)
setnames(fine_annotation, "query_key", "variant_key")
fine_mapped <- merge(
  fine_mapped,
  fine_annotation[, .(
    variant_key, overlaps_chromatin, n_chromatin_peak_overlaps, chromatin_features
  )],
  by = "variant_key",
  all.x = TRUE,
  sort = FALSE
)
setorder(fine_mapped, input_row)

message("Reading chromosome lengths for matched-control sampling...")
genome <- read_genome()
missing_genome <- setdiff(unique(as.character(seqnames(variant_ranges))), as.character(seqnames(genome)))
if (length(missing_genome)) {
  stop_with("Fine-mapped variants use chromosome(s) missing from the permutation genome: %s",
            paste(missing_genome, collapse = ", "))
}

message("Sampling TSS-distance/MAF matched PLINK SNPs for the negative control...")
matched_control_result <- sample_matched_controls(fine_mapped, genome, target_gene_column)
matched_controls <- matched_control_result$matched_controls
matched_pool <- matched_control_result$matched_pool
eligible_leads <- matched_control_result$eligible_leads
excluded_context <- matched_control_result$excluded_context
permutation_variant_index <- match(matched_controls$variant_key, names(variant_ranges))
if (anyNA(permutation_variant_index) || anyDuplicated(permutation_variant_index)) {
  stop_with("Matched-control variants could not be mapped one-to-one to fine-mapped SNP positions.")
}
permutation_variant_ranges <- variant_ranges[permutation_variant_index]
control_ranges <- GRanges(
  seqnames = matched_controls$control_seqname,
  ranges = IRanges(start = matched_controls$control_position, width = 1L),
  strand = "*"
)
names(control_ranges) <- as.character(matched_controls$signal_id)
if (length(control_ranges) != length(permutation_variant_ranges)) {
  stop_with("The matched control set (%d) does not equal the eligible fine-mapped SNP set (%d).",
            length(control_ranges), length(permutation_variant_ranges))
}
if (nrow(excluded_context)) {
  message("Using ", length(permutation_variant_ranges), "/", length(variant_ranges),
          " unique fine-mapped SNPs in matched-control permutations.")
}
if (anyDuplicated(matched_controls$control_variant_id)) {
  warning(
    sum(duplicated(matched_controls$control_variant_id)),
    " matched controls repeat a SNP because their gene-specific matching strata overlap. ",
    "The negative-control set still has one sampled control per fine-mapped SNP."
  )
}
control_query <- copy(matched_controls)
control_query[, query_key := as.character(signal_id)]
control_query <- control_query[, .(
  query_key, control_variant_id, control_seqname, control_position,
  target_gene, lead_variant_id, lead_seqname, lead_position, lead_MAF,
  control_MAF, abs_distance_TSS, maf_bin, distance_bin, variant_type,
  match_tier, n_matched_candidates
)]
control_overlap_table <- make_overlap_table(control_ranges, control_query)
control_annotation <- add_chromatin_annotation(control_query, control_overlap_table)

# Annotate the complete matched-control candidate universe once. Permutations
# then resample one candidate per eligible fine-mapped SNP while preserving its
# target-gene/TSS-distance, MAF-bin, and variant-type matching pool.
message("Annotating the complete matched-control candidate universe...")
candidate_variants <- unique(
  matched_pool[, .(
    control_variant_id,
    control_seqname,
    control_position
  )],
  by = "control_variant_id"
)
if (candidate_variants[, anyDuplicated(control_variant_id)]) {
  stop_with("Matched candidate variant IDs are not unique after deduplication.")
}
candidate_ranges <- GRanges(
  seqnames = candidate_variants$control_seqname,
  ranges = IRanges(start = candidate_variants$control_position, width = 1L),
  strand = "*"
)
names(candidate_ranges) <- candidate_variants$control_variant_id

observed_membership <- make_feature_membership(
  permutation_variant_ranges,
  chromatin_features
)
candidate_membership <- make_feature_membership(
  candidate_ranges,
  chromatin_features
)
if (!identical(names(observed_membership), names(candidate_membership))) {
  stop_with("Observed and matched-candidate chromatin feature sets differ.")
}

# Order candidate pools contiguously by signal. This allows one vectorized
# uniform draw per signal in every permutation without repeatedly splitting a
# potentially very large candidate table.
setorder(matched_pool, signal_id, control_variant_id)
pool_index <- matched_pool[, .(
  first_row = .I[1L],
  n_candidates = .N
), by = signal_id]
setorder(pool_index, signal_id)
expected_signal_ids <- sort(unique(eligible_leads$signal_id))
if (!identical(pool_index$signal_id, expected_signal_ids)) {
  stop_with("Matched pools do not provide exactly one candidate group per eligible signal.")
}
n_signals <- nrow(pool_index)
if (n_signals != length(permutation_variant_ranges)) {
  stop_with(
    "Matched pool contains %d signals but the eligible observed set contains %d variants.",
    n_signals,
    length(permutation_variant_ranges)
  )
}

run_feature_null <- function(label) {
  # Resetting to the same seed for every feature ensures that permutation b
  # uses the same sampled control for all annotations without storing a huge
  # signal-by-permutation index matrix.
  set.seed(cfg$seed + 1000L)
  pool_overlap <- matched_pool$control_variant_id %chin% candidate_membership[[label]]
  null_counts <- integer(cfg$n_permutations)
  for (permutation_index in seq_len(cfg$n_permutations)) {
    sampled_rows <- pool_index$first_row + floor(
      stats::runif(n_signals) * pool_index$n_candidates
    )
    null_counts[permutation_index] <- sum(pool_overlap[sampled_rows])
  }
  observed <- sum(names(permutation_variant_ranges) %chin% observed_membership[[label]])
  list(
    summary = summarise_matched_null(label, observed, null_counts, n_signals),
    null_counts = null_counts
  )
}

feature_labels <- names(observed_membership)
message(
  "Running ", cfg$n_permutations,
  " matched-control permutations for ", length(feature_labels),
  " chromatin annotations..."
)
if (
  cfg$use_multicore && cfg$n_cores > 1L &&
  .Platform$OS.type != "windows" && length(feature_labels) > 1L
) {
  permutation_results <- parallel::mclapply(
    feature_labels,
    run_feature_null,
    mc.cores = min(cfg$n_cores, length(feature_labels)),
    mc.preschedule = TRUE
  )
} else {
  permutation_results <- lapply(feature_labels, run_feature_null)
}
names(permutation_results) <- feature_labels
enrichment_summary <- rbindlist(lapply(permutation_results, `[[`, "summary"))
enrichment_summary[, p_fdr_bh := stats::p.adjust(p_two_sided, method = "BH")]
enrichment_summary[, p_enrichment_fdr_bh := stats::p.adjust(p_enrichment, method = "BH")]
enrichment_summary <- merge(
  enrichment_summary,
  tf_peak_qc,
  by = "chromatin_feature",
  all.x = TRUE,
  sort = FALSE
)
setcolorder(enrichment_summary, c(
  "chromatin_feature", "n_consensus_peaks", "consensus_covered_bp",
  setdiff(names(enrichment_summary), c(
    "chromatin_feature", "n_consensus_peaks", "consensus_covered_bp"
  ))
))
enrichment_summary[, sufficiently_powered :=
  expected_unique_variants >= 10 & observed_unique_variants >= 10]
setnames(enrichment_summary, "chromatin_feature", cfg$feature_result_column)
setorderv(enrichment_summary, cfg$feature_result_column)
null_counts <- lapply(permutation_results, `[[`, "null_counts")

fwrite(
  fine_mapped[, c("input_row", "variant_key") := NULL],
  file.path(cfg$output_dir, paste0("fine_mapped_variant_", cfg$feature_output_stem, "_annotations.tsv.gz")),
  sep = "\t"
)
if (cfg$save_overlap_annotations) {
  fwrite(
    fine_overlap_table,
    file.path(cfg$output_dir, paste0("fine_mapped_variant_", cfg$feature_output_stem, "_overlaps.tsv.gz")),
    sep = "\t"
  )
}
fwrite(
  excluded_context,
  file.path(cfg$output_dir, "fine_mapped_variants_excluded_from_matched_controls.tsv.gz"),
  sep = "\t"
)
fwrite(
  control_annotation,
  file.path(cfg$output_dir, paste0("matched_control_variant_", cfg$feature_output_stem, "_annotations.tsv.gz")),
  sep = "\t"
)
if (cfg$save_overlap_annotations) {
  fwrite(
    control_overlap_table,
    file.path(cfg$output_dir, paste0("matched_control_variant_", cfg$feature_output_stem, "_overlaps.tsv.gz")),
    sep = "\t"
  )
}
fwrite(
  matched_controls,
  file.path(cfg$output_dir, "matched_control_variants.tsv.gz"),
  sep = "\t"
)
fwrite(
  enrichment_summary,
  file.path(cfg$output_dir, paste0(cfg$feature_output_stem, "_overlap_matched_permutation.tsv")),
  sep = "\t"
)
if (cfg$save_null_counts) {
  saveRDS(
    list(
      null_counts = null_counts,
      observed_membership = observed_membership,
      seed = cfg$seed + 1000L,
      n_permutations = cfg$n_permutations,
      n_signals = n_signals
    ),
    file.path(cfg$output_dir, paste0(cfg$feature_output_stem, "_overlap_matched_permutation_null_counts.rds"))
  )
}

message("Wrote chromatin-feature annotations and matched-permutation enrichment results to: ", cfg$output_dir)
print(enrichment_summary)
