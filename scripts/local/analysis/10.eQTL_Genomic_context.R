#!/usr/bin/env Rscript
# Exploratory genomic-context annotation of fine-mapped eQTL variants.
#
# The input contains one row per variant--target-eGene (and credible-set)
# association.  A variant can therefore occur in more than one row.  The
# genome-wide annotations are invariant for a SNP; target-gene annotations are
# intentionally evaluated for every variant--gene row.

library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)
library(S4Vectors)
library(txdbmaker)

## Configuration -------------------------------------------------------------

GFF_PATH <- "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
# Use one lead variant per SuSiE credible set, rather than the full PIP table.
FINE_MAPPED_PATH <- "results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz"
# Effect estimates are joined by target eGene and lead-variant ID, and are
# retained in the final annotated table (log2_aFC, its error, c0, and MAF).
AFCN_PATH <- "results/summary/SusieR/aFCn_MAF_se_filtered.tsv.gz"
OUTPUT_DIR <- "results/summary/eQTL_genomic_context"

# These boundaries are relative to each transcript TSS.  `promoters()` is
# strand-aware: 2 kb is upstream and 200 bp is downstream of the TSS.
PROMOTER_UPSTREAM_BP <- 2000L
PROMOTER_DOWNSTREAM_BP <- 200L

# Set these explicitly if a differently named fine-mapping file is used.
VARIANT_ID_COLUMN <- "SNP"
TARGET_GENE_COLUMN <- "gene"
COMPONENT_RANK_CANDIDATES <- c("component_rank", "cs_rank", "cs", "component")

# Permutation enrichment is restricted to the same cis search space used for
# eQTL mapping: the union of +/- 2 Mb windows around all gene TSSs.
RUN_ENRICHMENT <- TRUE
CIS_TSS_WINDOW_BP <- 2000000L
ENRICHMENT_PERMUTATIONS <- 1000L
ENRICHMENT_SEED <- 20260730L
# `regioneR` uses parallel::mclapply() for permutation iterations on Unix-like
# systems. Increase this only as memory permits; set FALSE/1L for serial runs.
ENRICHMENT_USE_MULTICORE <- TRUE
ENRICHMENT_CORES <- 24L

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Small helpers -------------------------------------------------------------

# TxDb columns use upper-case names (for example, GENEID), while metadata
# carried by imported objects may use lower-case names.  This extracts a
# scalar identifier without manufacturing one when the source does not have it.
metadata_id <- function(gr, candidates, fallback_to_names = FALSE) {
  available <- colnames(mcols(gr))
  column <- candidates[candidates %in% available][1L]

  if (!is.na(column)) {
    value <- mcols(gr)[[column]]
    if (is(value, "List")) {
      value <- vapply(
        value,
        function(x) if (length(x) == 1L) as.character(x) else NA_character_,
        character(1)
      )
    } else {
      value <- as.character(value)
    }
    return(value)
  }

  if (fallback_to_names && !is.null(names(gr))) return(as.character(names(gr)))
  rep(NA_character_, length(gr))
}

# Keep the TxDb-provided metadata intact, and add canonical convenience columns
# only when they can be recovered from the TxDb object itself.
add_canonical_ids <- function(gr, feature = c("gene", "transcript")) {
  feature <- match.arg(feature)
  if (feature == "gene") {
    gene_id <- metadata_id(gr, c("GENEID", "gene_id"), fallback_to_names = TRUE)
    if (any(!is.na(gene_id))) mcols(gr)$gene_id <- gene_id
  } else {
    transcript_id <- metadata_id(
      gr, c("TXNAME", "tx_name", "transcript_id"), fallback_to_names = TRUE
    )
    gene_id <- metadata_id(gr, c("GENEID", "gene_id"))
    if (any(!is.na(transcript_id))) mcols(gr)$transcript_id <- transcript_id
    if (any(!is.na(gene_id))) mcols(gr)$gene_id <- gene_id
  }
  gr
}

# Flatten a transcript-grouped GRangesList while retaining its transcript
# membership and any available transcript/gene metadata.  This is deliberately
# not reduced: a range retains the feature-level information needed below.
flatten_transcript_feature <- function(grl, transcript_ranges) {
  gr <- unlist(grl, use.names = FALSE)
  parent_transcript_id <- rep(names(grl), elementNROWS(grl))
  tx_index <- match(parent_transcript_id, mcols(transcript_ranges)$transcript_id)

  if (length(gr)) {
    # Add (rather than replace) transcript metadata.  Feature-level columns
    # such as EXONID remain intact, while TXID/TXTYPE/GENEID are available on
    # every flattened transcript-derived range when the TxDb supplied them.
    for (metadata_column in colnames(mcols(transcript_ranges))) {
      if (!metadata_column %in% colnames(mcols(gr))) {
        mcols(gr)[[metadata_column]] <- mcols(transcript_ranges)[[metadata_column]][tx_index]
      }
    }
    mcols(gr)$transcript_id <- parent_transcript_id
    if ("gene_id" %in% colnames(mcols(transcript_ranges))) {
      mcols(gr)$gene_id <- mcols(transcript_ranges)$gene_id[tx_index]
    }
  }
  gr
}

# The ranges in `subject` have already been reduced.  Consequently, no
# feature, transcript, or gene metadata are attached to the output of reduce().
overlaps_any <- function(query, subject) {
  answer <- rep(FALSE, length(query))
  if (!length(query) || !length(subject)) return(answer)
  hits <- findOverlaps(query, subject, ignore.strand = TRUE)
  answer[unique(queryHits(hits))] <- TRUE
  answer
}

id_strings_from_overlaps <- function(query, subject, ids) {
  answer <- rep(NA_character_, length(query))
  if (!length(query) || !length(subject) || !any(!is.na(ids))) return(answer)

  hits <- findOverlaps(query, subject, ignore.strand = TRUE)
  if (!length(hits)) return(answer)
  hit_table <- data.table(
    query = queryHits(hits),
    id = as.character(ids[subjectHits(hits)])
  )[!is.na(id) & nzchar(id)]
  if (!nrow(hit_table)) return(answer)

  collapsed <- hit_table[, .(ids = paste(sort(unique(id)), collapse = ";")), by = query]
  answer[collapsed$query] <- collapsed$ids
  answer
}

# This evaluates all target-feature matches in one findOverlaps() call per
# feature, rather than iterating through fine-mapped variants.
overlaps_target_feature <- function(query, target_gene, subject, subject_gene_id) {
  answer <- rep(FALSE, length(query))
  if (!length(query) || !length(subject) || !any(!is.na(subject_gene_id))) return(answer)

  hits <- findOverlaps(query, subject, ignore.strand = TRUE)
  if (!length(hits)) return(answer)
  matching_hits <- !is.na(subject_gene_id[subjectHits(hits)]) &
    target_gene[queryHits(hits)] == subject_gene_id[subjectHits(hits)]
  answer[unique(queryHits(hits)[matching_hits])] <- TRUE
  answer
}

## Build the TxDb feature sets ------------------------------------------------

message("Building TxDb from: ", GFF_PATH)
txdb <- makeTxDbFromGFF(GFF_PATH)
txdb_columns <- columns(txdb)
transcript_columns <- intersect(c("TXID", "TXNAME", "GENEID", "TXTYPE"), txdb_columns)
gene_columns <- intersect(c("GENEID"), txdb_columns)

gene_ranges <- add_canonical_ids(
  # Unlike transcripts(), genes() does not have a use.names argument.  Its
  # range names are the TxDb gene identifiers when available.
  genes(txdb, columns = gene_columns), "gene"
)
transcript_ranges <- add_canonical_ids(
  transcripts(txdb, columns = transcript_columns, use.names = TRUE), "transcript"
)

# Promoters are based on transcript, not gene, TSSs. trim() only affects a
# promoter at a sequence boundary and preserves the range metadata.
promoter_ranges <- trim(promoters(
  transcript_ranges,
  upstream = PROMOTER_UPSTREAM_BP,
  downstream = PROMOTER_DOWNSTREAM_BP
))

# Keep both the transcript-grouped source objects and flattened original ranges.
# The latter retain per-feature metadata plus transcript/gene identifiers where
# the TxDb exposed them, and are used for target-gene-aware overlaps.
transcript_grouped_features <- list(
  `5UTR` = fiveUTRsByTranscript(txdb, use.names = TRUE),
  `3UTR` = threeUTRsByTranscript(txdb, use.names = TRUE),
  CDS = cdsBy(txdb, by = "tx", use.names = TRUE),
  exon = exonsBy(txdb, by = "tx", use.names = TRUE),
  intron = intronsByTranscript(txdb, use.names = TRUE)
)

# Set 1: original, metadata-retaining feature ranges.  Do not replace these
# with their reduced counterparts when identifiers or transcript structure are
# required.
feature_ranges_original <- c(
  list(genes = gene_ranges, transcripts = transcript_ranges, promoters = promoter_ranges),
  lapply(transcript_grouped_features, flatten_transcript_feature, transcript_ranges = transcript_ranges)
)

# Set 2: metadata-free genomic unions for the simple genome-wide overlap flags.
# reduce() merges all overlapping feature ranges; transcript and gene IDs are
# intentionally not represented in these objects.
feature_ranges_reduced <- lapply(feature_ranges_original, reduce)

message("Feature ranges (original / reduced):")
print(vapply(feature_ranges_original, length, integer(1)))
print(vapply(feature_ranges_reduced, length, integer(1)))

## Read lead variants and retain aFCn effect sizes ----------------------------

fine_mapped <- fread(FINE_MAPPED_PATH)
required_columns <- c(VARIANT_ID_COLUMN, TARGET_GENE_COLUMN)
missing_columns <- setdiff(required_columns, names(fine_mapped))
if (length(missing_columns)) {
  stop("Fine-mapping input is missing: ", paste(missing_columns, collapse = ", "))
}

# This is deliberately a left join: a lead variant without an aFCn estimate is
# still a valid SuSiE lead variant and must remain in the genomic annotation.
afcn <- fread(AFCN_PATH)
afcn_key_columns <- c("gene_id", "variant_id")
missing_afcn_columns <- setdiff(afcn_key_columns, names(afcn))
if (length(missing_afcn_columns)) {
  stop("aFCn input is missing: ", paste(missing_afcn_columns, collapse = ", "))
}
if (afcn[, anyDuplicated(.SD), .SDcols = afcn_key_columns]) {
  stop("aFCn input has duplicate gene_id/variant_id keys; cannot safely retain one effect estimate per lead row.")
}

fine_mapped <- merge(
  fine_mapped,
  afcn,
  by.x = c(TARGET_GENE_COLUMN, VARIANT_ID_COLUMN),
  by.y = afcn_key_columns,
  all.x = TRUE,
  sort = FALSE
)
afcn_value_columns <- setdiff(names(afcn), afcn_key_columns)
missing_afcn_estimates <- is.na(fine_mapped[[afcn_value_columns[1L]]])
if (any(missing_afcn_estimates)) {
  warning(
    sum(missing_afcn_estimates),
    " lead-variant rows have no matching aFCn estimate; their effect-size columns are NA."
  )
}

# The merge retains every effect-size column from aFCn in `fine_mapped`; this
# includes log2_aFC, log2_aFC_error, log2_aFC_c0, and MAF in the current input.
setnames(fine_mapped, c(VARIANT_ID_COLUMN, TARGET_GENE_COLUMN), c("variant_id", "target_gene"))
fine_mapped[, `:=`(variant_id = as.character(variant_id), target_gene = as.character(target_gene))]

# SNP IDs are expected to end in _<1-based position>, e.g. chr1_248511.
variant_seqname <- sub("_[^_]+$", "", fine_mapped$variant_id)
variant_position <- suppressWarnings(as.integer(sub("^.*_", "", fine_mapped$variant_id)))
bad_variant_id <- is.na(variant_position) | !nzchar(variant_seqname) | variant_position < 1L
if (any(bad_variant_id)) {
  stop(
    "Could not parse chromosome and 1-based position from ", sum(bad_variant_id),
    " variant IDs. Expected <seqname>_<position>, e.g. chr1_248511."
  )
}

fine_mapped[, `:=`(seqname = variant_seqname, position = variant_position)]
on_txdb_seqlevel <- fine_mapped$seqname %in% seqlevels(txdb)
if (any(!on_txdb_seqlevel)) {
  warning(
    sum(!on_txdb_seqlevel),
    " fine-mapped rows use sequence names absent from the TxDb; their annotations are NA."
  )
}

variant_ranges <- GRanges(
  seqnames = fine_mapped$seqname[on_txdb_seqlevel],
  ranges = IRanges(start = fine_mapped$position[on_txdb_seqlevel], width = 1L),
  strand = "*"
)

## Genome-wide overlap annotations -------------------------------------------

# Variants are unstranded, so every overlap explicitly uses ignore.strand=TRUE.
logical_annotation_columns <- c(
  in_promoter = "promoters",
  in_5UTR = "5UTR",
  in_3UTR = "3UTR",
  in_CDS = "CDS",
  in_exon = "exon",
  in_intron = "intron",
  in_gene_body = "genes",
  in_transcript = "transcripts"
)

for (annotation_name in names(logical_annotation_columns)) {
  values <- rep(NA, nrow(fine_mapped))
  values[on_txdb_seqlevel] <- overlaps_any(
    variant_ranges,
    feature_ranges_reduced[[logical_annotation_columns[[annotation_name]]]]
  )
  fine_mapped[, (annotation_name) := values]
}
fine_mapped[, intergenic := !in_gene_body]

# This priority is only for plotting and summarization.  The separate logical
# overlap columns above retain all overlapping annotations.  Thus a promoter
# variant that also lies in a UTR is assigned to its UTR primary category.
fine_mapped[, primary_annotation := NA_character_]
fine_mapped[on_txdb_seqlevel, primary_annotation := "intergenic"]
fine_mapped[in_gene_body %in% TRUE, primary_annotation := "other_gene_body"]
fine_mapped[in_intron %in% TRUE, primary_annotation := "intron"]
fine_mapped[in_promoter %in% TRUE, primary_annotation := "promoter"]
fine_mapped[in_exon %in% TRUE, primary_annotation := "other_exon"]
fine_mapped[in_CDS %in% TRUE, primary_annotation := "CDS"]
fine_mapped[in_3UTR %in% TRUE, primary_annotation := "3UTR"]
fine_mapped[in_5UTR %in% TRUE, primary_annotation := "5UTR"]

## Target-eGene-aware annotations --------------------------------------------

gene_ids <- metadata_id(feature_ranges_original$genes, "gene_id", fallback_to_names = TRUE)
gene_ids_available <- any(!is.na(gene_ids))
target_ids_compatible <- gene_ids_available && any(fine_mapped$target_gene %in% gene_ids, na.rm = TRUE)

# GFF-derived distance to the target eGene TSS. The signed distance is in
# reference coordinates (variant position - TSS); the absolute counterpart is
# appropriate for cis-distance matching and summaries. These fields remain NA
# when the TxDb cannot map a target identifier to a gene on the same sequence.
target_gene_tss <- rep(NA_integer_, nrow(fine_mapped))
if (gene_ids_available) {
  target_gene_index <- match(fine_mapped$target_gene, gene_ids)
  target_gene_seqname <- as.character(seqnames(feature_ranges_original$genes))
  target_gene_tss_values <- ifelse(
    as.character(strand(feature_ranges_original$genes)) == "-",
    end(feature_ranges_original$genes),
    start(feature_ranges_original$genes)
  )
  same_sequence <- !is.na(target_gene_index) &
    fine_mapped$seqname == target_gene_seqname[target_gene_index]
  target_gene_tss[same_sequence] <- target_gene_tss_values[target_gene_index[same_sequence]]
}
target_tss_distance <- fine_mapped$position - target_gene_tss
fine_mapped[, `:=`(
  target_gene_tss = target_gene_tss,
  distance_to_target_TSS = target_tss_distance,
  abs_distance_TSS = abs(target_tss_distance)
)]

fine_mapped[, `:=`(
  overlaps_any_gene = in_gene_body,
  overlaps_target_gene = NA,
  overlaps_other_gene = NA,
  target_gene_promoter = NA,
  target_gene_5UTR = NA,
  target_gene_3UTR = NA,
  target_gene_CDS = NA,
  target_gene_exon = NA,
  target_gene_intron = NA,
  overlapping_gene_ids = NA_character_,
  overlapping_transcript_ids = NA_character_
)]

# Transcript identifiers can be useful even when gene IDs cannot be aligned to
# target eGenes.  They are reported independently of the target-gene checks.
transcript_ids <- metadata_id(feature_ranges_original$transcripts, "transcript_id", fallback_to_names = TRUE)
if (any(!is.na(transcript_ids))) {
  transcript_strings <- id_strings_from_overlaps(
    variant_ranges, feature_ranges_original$transcripts, transcript_ids
  )
  fine_mapped[on_txdb_seqlevel, overlapping_transcript_ids := transcript_strings]
}

if (!target_ids_compatible) {
  warning(
    "Compatible TxDb gene identifiers could not be recovered or do not overlap ",
    "the target eGene IDs. Genome-wide annotations were produced, but target-gene ",
    "annotations and overlapping_gene_ids remain NA."
  )
} else {
  gene_hits <- findOverlaps(
    variant_ranges, feature_ranges_original$genes, ignore.strand = TRUE
  )
  target_overlap <- rep(FALSE, length(variant_ranges))
  other_overlap <- rep(FALSE, length(variant_ranges))
  if (length(gene_hits)) {
    hit_query <- queryHits(gene_hits)
    hit_gene <- gene_ids[subjectHits(gene_hits)]
    target_overlap[unique(hit_query[!is.na(hit_gene) &
      fine_mapped$target_gene[on_txdb_seqlevel][hit_query] == hit_gene])] <- TRUE
    other_overlap[unique(hit_query[!is.na(hit_gene) &
      fine_mapped$target_gene[on_txdb_seqlevel][hit_query] != hit_gene])] <- TRUE
  }

  fine_mapped[on_txdb_seqlevel, `:=`(
    overlaps_target_gene = target_overlap,
    overlaps_other_gene = other_overlap,
    overlapping_gene_ids = id_strings_from_overlaps(
      variant_ranges, feature_ranges_original$genes, gene_ids
    )
  )]

  target_feature_map <- c(
    target_gene_promoter = "promoters",
    target_gene_5UTR = "5UTR",
    target_gene_3UTR = "3UTR",
    target_gene_CDS = "CDS",
    target_gene_exon = "exon",
    target_gene_intron = "intron"
  )
  for (annotation_name in names(target_feature_map)) {
    subject <- feature_ranges_original[[target_feature_map[[annotation_name]]]]
    subject_gene_ids <- metadata_id(subject, "gene_id")
    if (any(!is.na(subject_gene_ids))) {
      fine_mapped[on_txdb_seqlevel, (annotation_name) := overlaps_target_feature(
        variant_ranges,
        target_gene[on_txdb_seqlevel],
        subject,
        subject_gene_ids
      )]
    }
  }
}

## Exploratory summaries ------------------------------------------------------

# One row per SNP for summaries whose unit is a unique genomic variant.
unique_variants <- unique(fine_mapped, by = "variant_id")
unique_variant_denominator <- nrow(unique_variants)

unique_variant_logical_annotation_summary <- rbindlist(lapply(
  c(names(logical_annotation_columns), "intergenic"),
  function(annotation_name) data.table(
    annotation = annotation_name,
    n_unique_variants = sum(unique_variants[[annotation_name]] %in% TRUE),
    pct_unique_variants = 100 * sum(unique_variants[[annotation_name]] %in% TRUE) /
      unique_variant_denominator
  )
))

primary_annotation_summary <- fine_mapped[
  !is.na(primary_annotation),
  .(n_variant_gene_rows = .N),
  by = primary_annotation
][, pct_variant_gene_rows := 100 * n_variant_gene_rows / nrow(fine_mapped)]
primary_annotation_summary[, primary_annotation_rank := match(
  primary_annotation,
  c("5UTR", "3UTR", "CDS", "other_exon", "promoter", "intron", "other_gene_body", "intergenic")
)]
setorder(primary_annotation_summary, primary_annotation_rank)
primary_annotation_summary[, primary_annotation_rank := NULL]

component_rank_column <- COMPONENT_RANK_CANDIDATES[
  COMPONENT_RANK_CANDIDATES %in% names(fine_mapped)
][1L]
component_rank_annotation_summary <- NULL
signal_rank_annotation_summary <- NULL
if (!is.na(component_rank_column)) {
  component_rank <- suppressWarnings(as.integer(as.character(fine_mapped[[component_rank_column]])))
  if (all(is.na(component_rank) | (component_rank > 0L & component_rank == floor(component_rank)))) {
    fine_mapped[, component_rank := component_rank]
    component_rank_annotation_summary <- fine_mapped[
      !is.na(primary_annotation) & !is.na(component_rank),
      .(n_variant_gene_rows = .N),
      by = .(component_rank, primary_annotation)
    ][, pct_within_component := 100 * n_variant_gene_rows / sum(n_variant_gene_rows), by = component_rank]

    # SuSiE components are only labelled primary/secondary when the source
    # values are positive integer ranks starting at one; otherwise no such
    # summary is emitted rather than guessing a signal order.
    observed_ranks <- sort(unique(component_rank[!is.na(component_rank)]))
    if (length(observed_ranks) && observed_ranks[1L] == 1L) {
      fine_mapped[, signal_rank := fifelse(component_rank == 1L, "primary", "secondary")]
      signal_rank_annotation_summary <- fine_mapped[
        !is.na(primary_annotation) & !is.na(signal_rank),
        .(n_variant_gene_rows = .N),
        by = .(signal_rank, primary_annotation)
      ][, pct_within_signal_rank := 100 * n_variant_gene_rows / sum(n_variant_gene_rows), by = signal_rank]
    }
  } else {
    warning("Found ", component_rank_column, " but it is not a positive integer rank; rank summaries were skipped.")
  }
}

target_gene_overlap_summary <- data.table(
  target_gene_relationship = c("overlaps_target_gene", "overlaps_other_gene", "overlaps_no_gene"),
  # Target status is specific to a variant--target-gene association.  Deduplicate
  # repeated credible-set rows, but do not collapse the same SNP across eGenes.
  n_variant_target_gene_pairs = {
    unique_variant_target_pairs <- unique(fine_mapped, by = c("variant_id", "target_gene"))
    c(
      sum(unique_variant_target_pairs$overlaps_target_gene %in% TRUE),
      sum(unique_variant_target_pairs$overlaps_other_gene %in% TRUE),
      sum(unique_variant_target_pairs$overlaps_any_gene %in% FALSE)
    )
  }
)
target_gene_overlap_summary[, pct_variant_target_gene_pairs := 100 *
  n_variant_target_gene_pairs / uniqueN(fine_mapped[, .(variant_id, target_gene)])]

unique_target_genes_by_primary_annotation <- fine_mapped[
  !is.na(primary_annotation) & !is.na(target_gene),
  .(n_unique_target_genes = uniqueN(target_gene)),
  by = primary_annotation
]

print(unique_variant_logical_annotation_summary)
print(primary_annotation_summary)
print(component_rank_annotation_summary)
print(signal_rank_annotation_summary)
print(target_gene_overlap_summary)
print(unique_target_genes_by_primary_annotation)

## Cis-background permutation enrichment -------------------------------------

# Test each logical genomic annotation separately. These feature classes can
# overlap (for example, a promoter can be intergenic), so this analysis does
# not use `primary_annotation`, whose priority is only for plotting/summaries.
enrichment_summary <- NULL
enrichment_permutation_tests <- NULL
cis_tss_background <- NULL

if (RUN_ENRICHMENT) {
  if (!requireNamespace("regioneR", quietly = TRUE)) {
    stop("RUN_ENRICHMENT is TRUE but the regioneR package is not installed.")
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("RUN_ENRICHMENT is TRUE but the rtracklayer package is not installed.")
  }
  if (!is.logical(ENRICHMENT_USE_MULTICORE) || length(ENRICHMENT_USE_MULTICORE) != 1L) {
    stop("ENRICHMENT_USE_MULTICORE must be TRUE or FALSE.")
  }
  if (!is.numeric(ENRICHMENT_CORES) || length(ENRICHMENT_CORES) != 1L ||
    is.na(ENRICHMENT_CORES) || ENRICHMENT_CORES < 1L ||
    ENRICHMENT_CORES != floor(ENRICHMENT_CORES)) {
    stop("ENRICHMENT_CORES must be a positive integer.")
  }

  # TxDb does not retain chromosome lengths for this GFF, so retrieve the
  # chromosome records directly and clip cis windows at chromosome boundaries.
  chromosome_ranges <- rtracklayer::import(
    GFF_PATH,
    format = "gff3",
    feature.type = "chromosome"
  )
  chromosome_end <- end(chromosome_ranges)[match(
    as.character(seqnames(gene_ranges)),
    as.character(seqnames(chromosome_ranges))
  )]
  bounded_gene <- !is.na(chromosome_end)
  if (any(!bounded_gene)) {
    warning(
      sum(!bounded_gene),
      " genes are on sequences without a chromosome-length record and were excluded ",
      "from the cis enrichment background."
    )
  }

  gene_tss <- gene_ranges[bounded_gene]
  gene_tss_position <- ifelse(
    as.character(strand(gene_tss)) == "-",
    end(gene_tss),
    start(gene_tss)
  )
  cis_tss_background <- reduce(GRanges(
    seqnames = seqnames(gene_tss),
    ranges = IRanges(
      start = pmax(1L, gene_tss_position - CIS_TSS_WINDOW_BP),
      end = pmin(
        chromosome_end[bounded_gene],
        gene_tss_position + CIS_TSS_WINDOW_BP
      )
    ),
    strand = "*"
  ))

  # Use distinct genomic positions, not repeated variant--gene/credible-set
  # rows, as the observed region set for enrichment.
  lead_variant_ranges <- unique(variant_ranges)
  in_cis_background <- overlaps_any(lead_variant_ranges, cis_tss_background)
  lead_variant_ranges <- lead_variant_ranges[in_cis_background]
  if (!length(lead_variant_ranges)) {
    stop("No lead variants overlap the +/- 2 Mb gene-TSS cis background.")
  }

  restrict_to_cis_background <- function(gr) {
    reduce(intersect(gr, cis_tss_background, ignore.strand = TRUE))
  }
  enrichment_feature_ranges <- list(
    promoter = restrict_to_cis_background(feature_ranges_reduced$promoters),
    `5UTR` = restrict_to_cis_background(feature_ranges_reduced$`5UTR`),
    `3UTR` = restrict_to_cis_background(feature_ranges_reduced$`3UTR`),
    CDS = restrict_to_cis_background(feature_ranges_reduced$CDS),
    exon = restrict_to_cis_background(feature_ranges_reduced$exon),
    intron = restrict_to_cis_background(feature_ranges_reduced$intron),
    intragenic = restrict_to_cis_background(feature_ranges_reduced$genes),
    transcript = restrict_to_cis_background(feature_ranges_reduced$transcripts),
    # This follows the script's definition: no overlap with any gene body.
    intergenic = setdiff(cis_tss_background, feature_ranges_reduced$genes,
      ignore.strand = TRUE
    )
  )
  empty_categories <- names(enrichment_feature_ranges)[
    !vapply(enrichment_feature_ranges, length, integer(1))
  ]
  if (length(empty_categories)) {
    warning("Skipping empty enrichment categories: ", paste(empty_categories, collapse = ", "))
    enrichment_feature_ranges <- enrichment_feature_ranges[
      vapply(enrichment_feature_ranges, length, integer(1)) > 0L
    ]
  }

  # overlapPermTest randomizes the one-base variants within `genome`. Here the
  # genome is the cis-background GRanges rather than the whole genome. Keeping
  # per.chromosome=TRUE preserves the chromosome distribution of lead variants.
  # regioneR reads the number of workers from the mc.cores option.
  previous_mc_cores <- getOption("mc.cores")
  options(mc.cores = as.integer(ENRICHMENT_CORES))
  set.seed(ENRICHMENT_SEED)
  enrichment_permutation_tests <- regioneR::overlapPermTest(
    A = lead_variant_ranges,
    B = enrichment_feature_ranges,
    alternative = "greater",
    ntimes = ENRICHMENT_PERMUTATIONS,
    genome = cis_tss_background,
    per.chromosome = TRUE,
    allow.overlaps = TRUE,
    count.once = TRUE,
    force.parallel = ENRICHMENT_USE_MULTICORE,
    verbose = FALSE
  )
  if (is.null(previous_mc_cores)) {
    options(mc.cores = NULL)
  } else {
    options(mc.cores = previous_mc_cores)
  }

  enrichment_summary <- rbindlist(Map(
    function(category, test_result) {
      null_counts <- test_result$permuted
      expected <- mean(null_counts)
      data.table(
        category = category,
        observed_unique_lead_variants = test_result$observed,
        expected_unique_lead_variants = expected,
        null_sd = stats::sd(null_counts),
        fold_enrichment = if (expected > 0) test_result$observed / expected else NA_real_,
        z_score = test_result$zscore,
        p_value = test_result$pval,
        n_permutations = test_result$ntimes,
        n_unique_lead_variants_in_cis_background = length(lead_variant_ranges),
        cis_background_bp = sum(width(cis_tss_background))
      )
    },
    names(enrichment_permutation_tests),
    enrichment_permutation_tests
  ))
  enrichment_summary[, p_value_fdr_bh := stats::p.adjust(p_value, method = "BH")]
  setorder(enrichment_summary, p_value, category)
  print(enrichment_summary)
}

## Write tabular outputs ------------------------------------------------------

fwrite(fine_mapped, file.path(OUTPUT_DIR, "fine_mapped_variant_annotations.tsv.gz"), sep = "\t")
fwrite(unique_variant_logical_annotation_summary, file.path(OUTPUT_DIR, "unique_variant_logical_annotation_summary.tsv"), sep = "\t")
fwrite(primary_annotation_summary, file.path(OUTPUT_DIR, "primary_annotation_summary.tsv"), sep = "\t")
fwrite(target_gene_overlap_summary, file.path(OUTPUT_DIR, "target_gene_overlap_summary.tsv"), sep = "\t")
fwrite(unique_target_genes_by_primary_annotation, file.path(OUTPUT_DIR, "unique_target_genes_by_primary_annotation.tsv"), sep = "\t")
if (!is.null(component_rank_annotation_summary)) {
  fwrite(component_rank_annotation_summary, file.path(OUTPUT_DIR, "component_rank_annotation_summary.tsv"), sep = "\t")
}
if (!is.null(signal_rank_annotation_summary)) {
  fwrite(signal_rank_annotation_summary, file.path(OUTPUT_DIR, "signal_rank_annotation_summary.tsv"), sep = "\t")
}
if (!is.null(enrichment_summary)) {
  fwrite(enrichment_summary, file.path(OUTPUT_DIR, "cis_feature_enrichment_regioneR.tsv"), sep = "\t")
  fwrite(
    data.table(
      seqname = as.character(seqnames(cis_tss_background)),
      start = start(cis_tss_background),
      end = end(cis_tss_background)
    ),
    file.path(OUTPUT_DIR, "cis_tss_background_plus_minus_2Mb.bed"),
    sep = "\t",
    col.names = FALSE
  )
  saveRDS(
    enrichment_permutation_tests,
    file.path(OUTPUT_DIR, "cis_feature_enrichment_regioneR_permutation_tests.rds")
  )
}

message("Wrote annotations and summaries to: ", OUTPUT_DIR)
