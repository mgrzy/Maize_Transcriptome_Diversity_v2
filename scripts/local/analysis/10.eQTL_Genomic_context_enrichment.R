#!/usr/bin/env Rscript
# Matched functional-enrichment test for fine-mapped cis-eQTL lead variants.
#
# Lead annotations come from 10.eQTL_Genomic_context.R.  Controls are drawn
# from PLINK's complete variant-frequency table, within +/- 1 Mb of the same
# target gene TSS, and matched by MAF bin, TSS-distance bin, and SNP/indel.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(IRanges)
  library(txdbmaker)
})

cfg <- list(
  # Produced by 10.eQTL_Genomic_context.R; contains the lead-variant primary
  # annotation and the target eGene ID for every credible-set lead variant.
  lead_file = "results/summary/eQTL_genomic_context/fine_mapped_variant_annotations.tsv.gz",
  # This is the complete genotype/tested-variant universe and provides MAF and
  # alleles. Position is parsed from its chr_position SNP identifier.
  plink_freq_file = "data/processed/geno/plink/plink.frq",
  gff_file = "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz",
  output_dir = "results/summary/eQTL_genomic_context_enrichment",

  # Matches scripts/local/analysis/5.SusieR.R (window_bp = 1,000,000).
  cis_window_bp = 1000000L,
  maf_bin_width = 0.05,
  distance_breaks = c(
    0,
    1000,
    2000,
    5000,
    10000,
    25000,
    50000,
    100000,
    250000,
    500000,
    1000000,
    Inf
  ),
  minimum_controls = 5L,
  # Retaining a random subset avoids storing billions of valid controls in
  # broad distance bins while preserving a large candidate pool per signal.
  max_controls_per_stratum = 500L,
  n_permutations = 5000L,
  seed = 20260730L,
  # Permutations are independent. Keep this modest because each worker reads
  # the (potentially large) matched control pool.
  use_multicore = TRUE,
  n_cores = 10L,
  candidate_batch_size = 500L,
  save_null_counts = TRUE
)

annotation_levels <- c("5UTR", "3UTR", "CDS", "other_exon", "promoter",
                       "intron", "other_gene_body", "intergenic")

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

assert_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop_with("%s is missing: %s", label, paste(missing, collapse = ", "))
}

parse_variant_coordinates <- function(ids) {
  seqname <- sub("_[^_]+$", "", ids)
  position <- suppressWarnings(as.integer(sub("^.*_", "", ids)))
  if (anyNA(position) || any(position < 1L) || any(!nzchar(seqname))) {
    stop_with("SNP IDs must be formatted as <seqname>_<1-based position>.")
  }
  data.table(variant_seqname = seqname, variant_position = position)
}

make_variant_type <- function(a1, a2) {
  fifelse(nchar(a1) == 1L & nchar(a2) == 1L, "SNP", "INDEL")
}

add_matching_bins <- function(x) {
  maf_breaks <- unique(c(seq(0, 0.5, by = cfg$maf_bin_width), 0.5))
  x[, maf_bin := cut(MAF, maf_breaks, include.lowest = TRUE, labels = FALSE)]
  x[, distance_bin := cut(abs_distance_TSS, cfg$distance_breaks,
                          include.lowest = TRUE, right = FALSE, labels = FALSE)]
  if (x[, anyNA(maf_bin) || anyNA(distance_bin)]) {
    stop_with("MAF or TSS distance falls outside the configured matching bins.")
  }
  invisible(x)
}

make_primary_annotation <- function(flags) {
  primary <- rep("intergenic", nrow(flags))
  primary[flags$in_gene_body] <- "other_gene_body"
  primary[flags$in_intron] <- "intron"
  primary[flags$in_promoter] <- "promoter"
  primary[flags$in_exon] <- "other_exon"
  primary[flags$in_CDS] <- "CDS"
  primary[flags$in_3UTR] <- "3UTR"
  primary[flags$in_5UTR] <- "5UTR"
  factor(primary, levels = annotation_levels)
}

overlaps_any <- function(query, subject) {
  result <- rep(FALSE, length(query))
  if (!length(query) || !length(subject)) return(result)
  hits <- findOverlaps(query, subject, ignore.strand = TRUE)
  result[unique(queryHits(hits))] <- TRUE
  result
}

# Build the same genome-wide feature definitions used for the lead annotation.
build_gff_context <- function(gff_file) {
  txdb <- makeTxDbFromGFF(gff_file)
  genes_gr <- genes(txdb)
  tx_gr <- transcripts(txdb, use.names = TRUE)
  features <- list(
    genes = genes_gr,
    promoters = trim(promoters(tx_gr, upstream = 2000L, downstream = 200L)),
    `5UTR` = unlist(fiveUTRsByTranscript(txdb, use.names = TRUE), use.names = FALSE),
    `3UTR` = unlist(threeUTRsByTranscript(txdb, use.names = TRUE), use.names = FALSE),
    CDS = unlist(cdsBy(txdb, by = "tx", use.names = TRUE), use.names = FALSE),
    exon = unlist(exonsBy(txdb, by = "tx", use.names = TRUE), use.names = FALSE),
    intron = unlist(intronsByTranscript(txdb, use.names = TRUE), use.names = FALSE)
  )
  reduced_features <- lapply(features, reduce)

  gene_ids <- names(genes_gr)
  if (is.null(gene_ids) || any(!nzchar(gene_ids))) {
    stop_with("The GFF-derived TxDb does not expose usable gene identifiers.")
  }
  gene_tss <- ifelse(as.character(strand(genes_gr)) == "-", end(genes_gr), start(genes_gr))
  tss <- data.table(
    target_gene = as.character(gene_ids),
    variant_seqname = as.character(seqnames(genes_gr)),
    gene_tss = as.integer(gene_tss)
  )
  if (tss[, anyDuplicated(target_gene)]) stop_with("GFF has duplicate gene identifiers.")

  chromosome_ranges <- rtracklayer::import(gff_file, format = "gff3", feature.type = "chromosome")
  chromosome_ends <- data.table(
    variant_seqname = as.character(seqnames(chromosome_ranges)),
    chromosome_end = as.integer(end(chromosome_ranges))
  )
  list(tss = tss, chromosome_ends = chromosome_ends, features = reduced_features)
}

annotate_control_variants <- function(variants, reduced_features) {
  gr <- GRanges(
    seqnames = variants$variant_seqname,
    ranges = IRanges(start = variants$variant_position, width = 1L),
    strand = "*"
  )
  flags <- data.table(
    variant_id = variants$variant_id,
    in_promoter = overlaps_any(gr, reduced_features$promoters),
    in_5UTR = overlaps_any(gr, reduced_features$`5UTR`),
    in_3UTR = overlaps_any(gr, reduced_features$`3UTR`),
    in_CDS = overlaps_any(gr, reduced_features$CDS),
    in_exon = overlaps_any(gr, reduced_features$exon),
    in_intron = overlaps_any(gr, reduced_features$intron),
    in_gene_body = overlaps_any(gr, reduced_features$genes)
  )
  flags[, primary_annotation := make_primary_annotation(.SD),
        .SDcols = setdiff(names(flags), "variant_id")]
  flags[, .(variant_id, primary_annotation)]
}

# Every absolute-distance bin is represented by an upstream and downstream
# interval around the target gene's GFF-derived TSS.  Both are clipped to the
# chromosome bounds and to the +/- 1 Mb tested cis window.
make_stratum_intervals <- function(strata, chromosome_ends) {
  x <- merge(strata, chromosome_ends, by = "variant_seqname", all.x = TRUE, sort = FALSE)
  if (x[, anyNA(chromosome_end)]) {
    stop_with("Some target-gene sequence names have no chromosome bounds in the GFF.")
  }
  lower <- cfg$distance_breaks[x$distance_bin]
  upper <- pmin(cfg$distance_breaks[x$distance_bin + 1L] - 1L, cfg$cis_window_bp)
  if (any(!is.finite(upper))) upper[!is.finite(upper)] <- cfg$cis_window_bp

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

count_categories <- function(x, name) {
  observed <- x[, .N, by = primary_annotation]
  template <- data.table(primary_annotation = factor(annotation_levels, levels = annotation_levels))
  observed <- observed[template, on = "primary_annotation"]
  observed[is.na(N), N := 0L]
  setnames(observed, "N", name)
  observed
}

sample_one_per_signal <- function(pool) {
  indices <- pool[, .I[sample.int(.N, 1L)], by = signal_id]$V1
  pool[indices]
}

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading GFF and deriving gene TSS positions...")
gff_context <- build_gff_context(cfg$gff_file)

message("Reading annotated lead variants...")
lead <- fread(cfg$lead_file)
assert_columns(lead, c("target_gene", "variant_id", "cs", "primary_annotation"), "Annotated lead table")
lead <- unique(lead, by = c("target_gene", "cs"))
lead[, `:=`(target_gene = as.character(target_gene), variant_id = as.character(variant_id))]
lead_coordinates <- parse_variant_coordinates(lead$variant_id)
lead[, `:=`(variant_seqname = lead_coordinates$variant_seqname,
            variant_position = lead_coordinates$variant_position)]
lead[, primary_annotation := factor(as.character(primary_annotation), levels = annotation_levels)]
if (lead[, anyNA(primary_annotation)]) stop_with("Lead table has unrecognised primary annotations.")

message("Reading PLINK frequency universe (MAF, alleles, and SNP positions)...")
plink <- fread(cfg$plink_freq_file, select = c("CHR", "SNP", "A1", "A2", "MAF"))
setnames(plink, c("SNP", "A1", "A2"), c("variant_id", "A1", "A2"))
plink[, variant_id := as.character(variant_id)]
plink_coordinates <- parse_variant_coordinates(plink$variant_id)
plink[, `:=`(variant_seqname = plink_coordinates$variant_seqname,
             variant_position = plink_coordinates$variant_position,
             MAF = as.numeric(MAF),
             variant_type = make_variant_type(as.character(A1), as.character(A2)))]
plink <- plink[is.finite(MAF) & MAF > 0 & MAF <= 0.5]
if (plink[, anyDuplicated(variant_id)]) stop_with("PLINK frequency table has duplicate SNP IDs.")
maf_breaks <- unique(c(seq(0, 0.5, by = cfg$maf_bin_width), 0.5))
plink[, maf_bin := cut(MAF, maf_breaks, include.lowest = TRUE, labels = FALSE)]
if (plink[, anyNA(maf_bin)]) stop_with("PLINK MAF values are not covered by maf_bin_width.")

# Use PLINK MAF for every lead, including the few lead variants without aFCn.
# The annotated lead table already contains an aFCn-derived MAF column; remove
# it first so this merge has one explicit, complete PLINK MAF field.
if ("MAF" %in% names(lead)) lead[, MAF := NULL]
lead <- merge(lead, plink[, .(variant_id, MAF, variant_type)], by = "variant_id",
              all.x = TRUE, sort = FALSE)
lead <- merge(lead, gff_context$tss, by = c("target_gene", "variant_seqname"),
              all.x = TRUE, sort = FALSE)
if (lead[, anyNA(MAF) | anyNA(gene_tss)]) {
  stop_with("Some lead variants lack PLINK MAF or a matching GFF target-gene TSS.")
}
lead[, `:=`(
  # Positive values are downstream in reference coordinates; negative values
  # are upstream in reference coordinates. abs_distance_TSS is used for bins.
  distance_to_target_TSS = variant_position - gene_tss,
  abs_distance_TSS = abs(variant_position - gene_tss)
)]
lead[, within_gff_cis_window := abs_distance_TSS <= cfg$cis_window_bp]
outside_cis_window <- !lead$within_gff_cis_window
if (any(outside_cis_window)) {
  warning(
    sum(outside_cis_window),
    " lead signals are outside the GFF-derived +/- ", cfg$cis_window_bp,
    " bp window and will be excluded so controls use the same cis definition."
  )
  lead <- lead[!outside_cis_window]
}
if (!nrow(lead)) stop_with("No lead variants remain within the configured cis window.")
add_matching_bins(lead)
lead[, signal_id := .I]

# A stratum is shared only where every matching property, including target gene,
# is identical. Candidate variants are then expanded back to individual signals.
stratum_columns <- c("target_gene", "variant_seqname", "gene_tss", "maf_bin", "distance_bin", "variant_type")
lead[, stratum_id := .GRP, by = stratum_columns]
strata <- unique(lead[, c("stratum_id", stratum_columns), with = FALSE])
intervals <- make_stratum_intervals(strata, gff_context$chromosome_ends)

# Efficient interval joins against all PLINK variants. Each broad stratum is
# capped after a random draw so later permutations remain practical.
setkey(plink, variant_seqname, maf_bin, variant_type, variant_position)
set.seed(cfg$seed)
control_parts <- vector("list", ceiling(nrow(intervals) / cfg$candidate_batch_size))
message("Building GFF-distance/MAF/type matched control pools...")
for (batch_index in seq_along(control_parts)) {
  first <- (batch_index - 1L) * cfg$candidate_batch_size + 1L
  last <- min(nrow(intervals), batch_index * cfg$candidate_batch_size)
  interval_batch <- intervals[first:last]
  hits <- plink[interval_batch,
    on = .(variant_seqname, maf_bin, variant_type,
           variant_position >= interval_start, variant_position <= interval_end),
    nomatch = 0L, allow.cartesian = TRUE,
    .(stratum_id = i.stratum_id, variant_id = x.variant_id,
      variant_seqname = x.variant_seqname, variant_position = x.variant_position,
      MAF = x.MAF, variant_type = x.variant_type)
  ]
  if (nrow(hits)) {
    hits <- unique(hits, by = c("stratum_id", "variant_id"))
    hits <- hits[, .SD[sample.int(.N, min(.N, cfg$max_controls_per_stratum))], by = stratum_id]
  }
  control_parts[[batch_index]] <- hits
  if (batch_index %% 10L == 0L || batch_index == length(control_parts)) {
    message("  processed ", last, "/", nrow(intervals), " distance intervals")
  }
}
control_candidates <- rbindlist(control_parts)
rm(control_parts, intervals)

matched_pool <- merge(
  lead[, .(signal_id, stratum_id, lead_variant = variant_id)],
  control_candidates,
  by = "stratum_id", allow.cartesian = TRUE, sort = FALSE
)[variant_id != lead_variant]
rm(control_candidates)

pool_sizes <- matched_pool[, .(n_controls = uniqueN(variant_id)), by = signal_id]
lead <- pool_sizes[lead, on = "signal_id"]
lead[is.na(n_controls), n_controls := 0L]
lead[, eligible := n_controls >= cfg$minimum_controls]
matched_pool <- matched_pool[signal_id %in% lead[eligible == TRUE, signal_id]]
lead_matched <- lead[eligible == TRUE]
if (!nrow(lead_matched)) stop_with("No lead signals have the required number of controls.")

message("Annotating selected control variants from the GFF...")
control_annotation <- annotate_control_variants(
  unique(matched_pool[, .(variant_id, variant_seqname, variant_position)]),
  gff_context$features
)
matched_pool <- control_annotation[matched_pool, on = "variant_id"]
if (matched_pool[, anyNA(primary_annotation)]) stop_with("Some control variants could not be annotated.")

message("Retained ", nrow(lead_matched), "/", nrow(lead), " lead signals with >= ",
        cfg$minimum_controls, " controls.")
observed_counts <- count_categories(lead_matched, "observed")

set.seed(cfg$seed + 1L)
run_permutation <- function(iteration) {
  sampled <- sample_one_per_signal(matched_pool)
  count_categories(sampled, "null_count")[, permutation := iteration]
}

if (!is.logical(cfg$use_multicore) || length(cfg$use_multicore) != 1L ||
    !is.numeric(cfg$n_cores) || length(cfg$n_cores) != 1L ||
    is.na(cfg$n_cores) || cfg$n_cores < 1L || cfg$n_cores != floor(cfg$n_cores)) {
  stop_with("use_multicore must be logical and n_cores must be a positive integer.")
}
run_parallel <- isTRUE(cfg$use_multicore) && cfg$n_cores > 1L && .Platform$OS.type != "windows"
if (run_parallel) {
  message("Running ", cfg$n_permutations, " permutations on ", cfg$n_cores, " cores...")
  null_list <- parallel::mclapply(
    seq_len(cfg$n_permutations), run_permutation,
    mc.cores = as.integer(cfg$n_cores), mc.set.seed = TRUE
  )
} else {
  if (isTRUE(cfg$use_multicore) && .Platform$OS.type == "windows") {
    warning("Multicore permutations are unavailable on Windows; running serially.")
  }
  null_list <- vector("list", cfg$n_permutations)
  for (iteration in seq_len(cfg$n_permutations)) {
    null_list[[iteration]] <- run_permutation(iteration)
    if (iteration %% max(1L, floor(cfg$n_permutations / 10L)) == 0L) {
      message("  completed ", iteration, "/", cfg$n_permutations, " permutations")
    }
  }
}
null_counts <- rbindlist(null_list)

enrichment_summary <- merge(null_counts, observed_counts, by = "primary_annotation",
                            allow.cartesian = TRUE, sort = FALSE)[, {
  expected <- mean(null_count)
  .(observed = unique(observed), expected = expected, null_sd = sd(null_count),
    fold_enrichment = if (expected > 0) unique(observed) / expected else NA_real_,
    p_enrichment = (1 + sum(null_count >= unique(observed))) / (.N + 1L),
    p_depletion = (1 + sum(null_count <= unique(observed))) / (.N + 1L),
    n_permutations = .N)
}, by = primary_annotation]
enrichment_summary[, p_two_sided := pmin(1, 2 * pmin(p_enrichment, p_depletion))]
enrichment_summary[, p_fdr_bh := p.adjust(p_two_sided, method = "BH")]
enrichment_summary[, annotation_rank := match(primary_annotation, annotation_levels)]
setorder(enrichment_summary, annotation_rank)
enrichment_summary[, annotation_rank := NULL]

run_summary <- data.table(
  lead_rows = nrow(lead), eligible_lead_rows = nrow(lead_matched),
  plink_variants = nrow(plink), cis_window_bp = cfg$cis_window_bp,
  n_permutations = cfg$n_permutations, max_controls_per_stratum = cfg$max_controls_per_stratum,
  median_controls = median(lead_matched$n_controls), seed = cfg$seed
)
fwrite(enrichment_summary, file.path(cfg$output_dir, "matched_enrichment_summary.tsv"), sep = "\t")
fwrite(run_summary, file.path(cfg$output_dir, "run_summary.tsv"), sep = "\t")
fwrite(lead[, .(target_gene, cs, variant_id, MAF, distance_to_target_TSS,
                abs_distance_TSS, primary_annotation, n_controls, eligible)],
       file.path(cfg$output_dir, "signal_matching_diagnostics.tsv.gz"), sep = "\t")
if (cfg$save_null_counts) fwrite(null_counts,
  file.path(cfg$output_dir, "permutation_null_counts.tsv.gz"), sep = "\t")

message("\nMatched enrichment results:")
print(enrichment_summary)
message("Wrote results to: ", normalizePath(cfg$output_dir))
