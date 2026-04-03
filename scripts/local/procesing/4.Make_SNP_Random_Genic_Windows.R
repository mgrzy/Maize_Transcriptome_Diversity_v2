suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
})

get_script_path <- function() {
  script_path_arg <- "--file="
  script_path_match <- commandArgs()[grepl(paste0("^", script_path_arg), commandArgs())]

  if (length(script_path_match)) {
    return(normalizePath(sub(script_path_arg, "", script_path_match[1L]), winslash = "/", mustWork = TRUE))
  }

  source_refs <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  if (length(source_refs)) {
    return(normalizePath(source_refs[[length(source_refs)]], winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

resolve_project_root <- function() {
  script_path <- get_script_path()

  if (!is.na(script_path)) {
    return(normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE))
  }

  normalizePath(".", winslash = "/", mustWork = TRUE)
}

project_root <- resolve_project_root()

default_snp_csv <- file.path(project_root, "results", "full", "SusieR2", "PIP_Max.csv")
default_gff <- "/home/marcin/R/BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
default_out_dir <- file.path(project_root, "results", "summary", "bed_windows")

make_snp_random_genic_windows_config <- list(
  snp_csv = default_snp_csv,
  gff = default_gff,
  out_dir = default_out_dir,
  window_bp = 1000L,
  promoter_bp = 2000L,
  seed = 1L,
  chr_col = NA_character_,
  pos_col = "POS",
  id_col = "SNP",
  variant_col = "SNP"
)

find_first_column <- function(dt, candidates) {
  matches <- candidates[candidates %in% names(dt)]
  if (!length(matches)) {
    return(NA_character_)
  }

  matches[[1L]]
}

choose_column <- function(dt, explicit, candidates, label, required = TRUE) {
  if (!is.na(explicit)) {
    if (!explicit %in% names(dt)) {
      stop(
        "Requested ", label, " column '", explicit, "' was not found. Available columns: ",
        paste(names(dt), collapse = ", ")
      )
    }
    return(explicit)
  }

  detected <- find_first_column(dt, candidates)
  if (!is.na(detected) || !required) {
    return(detected)
  }

  stop(
    "Could not detect a ", label, " column. Available columns: ",
    paste(names(dt), collapse = ", ")
  )
}

normalize_seqname_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  tolower(x)
}

harmonize_seqnames <- function(x, target_seqlevels) {
  target_keys <- normalize_seqname_key(target_seqlevels)
  if (anyDuplicated(target_keys)) {
    stop("Could not build a unique chromosome mapping from the GFF seqlevels.")
  }

  idx <- match(normalize_seqname_key(x), target_keys)
  out <- target_seqlevels[idx]
  out[is.na(x) | !nzchar(trimws(as.character(x)))] <- NA_character_
  out
}

set_gr_seqlengths <- function(gr, genome_dt) {
  current_seqlevels <- seqlevels(gr)
  current_lengths <- genome_dt$chrom_length[match(current_seqlevels, genome_dt$chrom)]
  seqlengths(gr) <- current_lengths
  gr
}

parse_variant_column <- function(x) {
  x <- trimws(as.character(x))
  is_match <- grepl("^[^:_]+[:_][0-9]+(?:[:_].*)?$", x)

  chrom <- rep(NA_character_, length(x))
  pos <- rep(NA_real_, length(x))

  chrom[is_match] <- sub("^([^:_]+)[:_][0-9]+(?:[:_].*)?$", "\\1", x[is_match])
  pos[is_match] <- suppressWarnings(as.numeric(sub("^[^:_]+[:_]([0-9]+)(?:[:_].*)?$", "\\1", x[is_match])))

  list(chrom = chrom, pos = pos)
}

load_annotation <- function(gff_path) {
  message("Loading chromosome and gene features from: ", gff_path)

  gr <- import(gff_path, format = "gff3", feature.type = c("chromosome", "gene"))
  feature_type <- as.character(mcols(gr)$type)

  chrom_gr <- gr[feature_type == "chromosome"]
  gene_gr <- gr[feature_type == "gene"]

  if (!length(gene_gr)) {
    stop("No gene features were found in: ", gff_path)
  }

  if (length(chrom_gr)) {
    genome_dt <- unique(data.table(
      chrom = as.character(seqnames(chrom_gr)),
      chrom_length = as.integer(end(chrom_gr))
    ))
  } else {
    seq_lengths <- seqlengths(gr)
    seq_lengths <- seq_lengths[!is.na(seq_lengths)]
    if (!length(seq_lengths)) {
      stop("No chromosome lengths were available from the GFF.")
    }
    genome_dt <- data.table(
      chrom = names(seq_lengths),
      chrom_length = as.integer(seq_lengths)
    )
  }

  genome_dt <- genome_dt[!is.na(chrom_length) & chrom_length > 0L]
  if (!nrow(genome_dt)) {
    stop("No usable chromosome lengths were found in: ", gff_path)
  }

  genome_dt[, chrom := as.character(chrom)]
  genome_dt <- unique(genome_dt, by = "chrom")
  genome_dt[, chrom_order := seq_len(.N)]

  gene_gr <- keepSeqlevels(gene_gr, genome_dt$chrom, pruning.mode = "coarse")
  gene_gr <- set_gr_seqlengths(gene_gr, genome_dt)

  list(
    genome = genome_dt,
    genes = gene_gr
  )
}

read_snp_table <- function(path, genome_dt, chr_col = NA_character_, pos_col = NA_character_,
                           id_col = NA_character_, variant_col = NA_character_) {
  if (!file.exists(path)) {
    stop("SNP table not found: ", path)
  }

  message("Reading SNP table: ", path)
  dt <- fread(path, showProgress = FALSE)

  if (!nrow(dt)) {
    stop("SNP table is empty: ", path)
  }

  variant_col <- choose_column(
    dt,
    variant_col,
    c("SNP", "snp", "variant", "variant_id", "marker", "id", "ID"),
    "variant",
    required = FALSE
  )
  chr_col <- choose_column(
    dt,
    chr_col,
    c("chr", "chrom", "chromosome", "CHR", "CHROM", "Chromosome", "seqnames"),
    "chr",
    required = FALSE
  )
  pos_col <- choose_column(
    dt,
    pos_col,
    c("pos", "position", "bp", "BP", "POS", "Position", "start"),
    "pos",
    required = FALSE
  )

  if (is.na(chr_col) && is.na(variant_col)) {
    stop(
      "Could not detect a chromosome column. Available columns: ",
      paste(names(dt), collapse = ", "),
      ". Provide chr_col or variant_col."
    )
  }
  if (is.na(pos_col) && is.na(variant_col)) {
    stop(
      "Could not detect a position column. Available columns: ",
      paste(names(dt), collapse = ", "),
      ". Provide pos_col or variant_col."
    )
  }

  if (is.na(id_col)) {
    id_col <- choose_column(
      dt,
      id_col,
      c("SNP", "snp", "id", "ID", "rsid", "variant_id", "marker"),
      "id",
      required = FALSE
    )
  } else if (!id_col %in% names(dt)) {
    stop("Requested id column '", id_col, "' was not found.")
  }

  variant_info <- NULL
  if (!is.na(variant_col)) {
    variant_info <- parse_variant_column(dt[[variant_col]])
  }

  raw_chr <- if (!is.na(chr_col)) as.character(dt[[chr_col]]) else variant_info$chrom
  raw_pos <- if (!is.na(pos_col)) dt[[pos_col]] else variant_info$pos

  snp_dt <- data.table(
    raw_chr = raw_chr,
    raw_pos = raw_pos
  )

  if (!is.na(id_col)) {
    snp_dt[, snp_id := as.character(dt[[id_col]])]
  } else if (!is.na(variant_col)) {
    snp_dt[, snp_id := as.character(dt[[variant_col]])]
  } else {
    snp_dt[, snp_id := paste0("snp_", seq_len(.N))]
  }

  snp_dt[, pos := suppressWarnings(as.numeric(raw_pos))]
  invalid_numeric <- is.na(snp_dt$pos)
  if (any(invalid_numeric)) {
    stop("Found ", sum(invalid_numeric), " SNP rows with non-numeric positions.")
  }

  non_integer <- abs(snp_dt$pos - round(snp_dt$pos)) > .Machine$double.eps^0.5
  if (any(non_integer)) {
    stop("Found ", sum(non_integer), " SNP rows with non-integer positions.")
  }

  snp_dt[, chrom := harmonize_seqnames(raw_chr, genome_dt$chrom)]
  unmatched_chr <- is.na(snp_dt$chrom)
  if (any(unmatched_chr)) {
    bad_chr <- unique(snp_dt$raw_chr[unmatched_chr])
    bad_chr <- bad_chr[seq_len(min(length(bad_chr), 10L))]
    stop(
      "Found ", sum(unmatched_chr), " SNP rows with chromosomes not present in the GFF. Examples: ",
      paste(bad_chr, collapse = ", ")
    )
  }

  snp_dt[is.na(snp_id) | !nzchar(snp_id), snp_id := paste0("snp_", .I)]
  snp_dt[, pos := as.integer(round(pos))]

  chrom_lengths <- setNames(genome_dt$chrom_length, genome_dt$chrom)
  out_of_bounds <- snp_dt$pos < 1L | snp_dt$pos > chrom_lengths[snp_dt$chrom]
  if (any(out_of_bounds)) {
    stop("Found ", sum(out_of_bounds), " SNP rows with positions outside chromosome bounds.")
  }

  snp_dt <- unique(snp_dt[, .(chrom, pos, snp_id)])
  snp_dt[, chrom_order := genome_dt$chrom_order[match(chrom, genome_dt$chrom)]]
  setorder(snp_dt, chrom_order, pos, snp_id)
  snp_dt[, chrom_order := NULL]

  snp_dt
}

make_snp_windows <- function(snp_dt, genome_dt, window_bp) {
  left_flank <- as.integer(window_bp %/% 2L)
  right_flank <- as.integer(window_bp - left_flank - 1L)

  gr <- GRanges(
    seqnames = snp_dt$chrom,
    ranges = IRanges(
      start = snp_dt$pos - left_flank,
      end = snp_dt$pos + right_flank
    ),
    strand = "*"
  )

  gr <- set_gr_seqlengths(gr, genome_dt)
  gr <- trim(gr)
  mcols(gr)$name <- snp_dt$snp_id
  gr
}

make_random_windows <- function(n_windows, genome_dt, window_bp, seed) {
  width <- as.integer(window_bp)
  available <- pmax(genome_dt$chrom_length - width + 1L, 0L)

  if (!n_windows) {
    return(GRanges())
  }

  if (!any(available > 0L)) {
    stop("No chromosomes are long enough to place random windows of width ", width, ".")
  }

  if (n_windows > sum(available)) {
    stop(
      "Requested ", n_windows, " random windows, but only ",
      sum(available), " unique windows of width ", width, " are possible."
    )
  }

  set.seed(seed)

  valid_idx <- which(available > 0L)
  sampled_counts <- as.integer(rmultinom(1L, size = n_windows, prob = available[valid_idx]))

  chrom_vec <- character(0)
  start_vec <- integer(0)

  for (i in seq_along(valid_idx)) {
    chrom_idx <- valid_idx[[i]]
    n_chr <- sampled_counts[[i]]
    if (!n_chr) {
      next
    }

    chr_starts <- sample.int(available[[chrom_idx]], size = n_chr, replace = FALSE)
    chrom_vec <- c(chrom_vec, rep(genome_dt$chrom[[chrom_idx]], n_chr))
    start_vec <- c(start_vec, chr_starts)
  }

  random_dt <- data.table(
    chrom = chrom_vec,
    start = start_vec,
    end = start_vec + width - 1L
  )
  random_dt[, chrom_order := genome_dt$chrom_order[match(chrom, genome_dt$chrom)]]
  setorder(random_dt, chrom_order, start, end)
  random_dt[, chrom_order := NULL]
  random_dt[, name := sprintf("random_%06d", seq_len(.N))]

  gr <- GRanges(
    seqnames = random_dt$chrom,
    ranges = IRanges(start = random_dt$start, end = random_dt$end),
    strand = "*"
  )
  gr <- set_gr_seqlengths(gr, genome_dt)
  mcols(gr)$name <- random_dt$name
  gr
}

make_gene_outputs <- function(gene_gr, genome_dt, promoter_bp) {
  gene_ids <- as.character(mcols(gene_gr)$ID)
  missing_ids <- is.na(gene_ids) | !nzchar(gene_ids)
  gene_ids[missing_ids] <- sprintf("gene_%06d", which(missing_ids))

  mcols(gene_gr)$name <- gene_ids

  promoter_start <- ifelse(as.character(strand(gene_gr)) == "-", end(gene_gr) + 1L, start(gene_gr) - promoter_bp)
  promoter_end <- ifelse(as.character(strand(gene_gr)) == "-", end(gene_gr) + promoter_bp, start(gene_gr) - 1L)

  promoter_gr <- GRanges(
    seqnames = seqnames(gene_gr),
    ranges = IRanges(start = promoter_start, end = promoter_end),
    strand = strand(gene_gr)
  )
  promoter_gr <- set_gr_seqlengths(promoter_gr, genome_dt)
  promoter_gr <- trim(promoter_gr)
  keep_idx <- width(promoter_gr) > 0L
  promoter_gr <- promoter_gr[keep_idx]
  mcols(promoter_gr)$name <- mcols(gene_gr)$name[keep_idx]

  list(
    genes = gene_gr,
    promoters = promoter_gr
  )
}

gr_to_bed <- function(gr, genome_dt) {
  bed_dt <- data.table(
    chrom = as.character(seqnames(gr)),
    start = start(gr) - 1L,
    end = end(gr),
    name = if ("name" %in% colnames(mcols(gr))) as.character(mcols(gr)$name) else ".",
    score = ".",
    strand = as.character(strand(gr))
  )

  bed_dt[strand == "*" | is.na(strand), strand := "."]
  bed_dt[is.na(name) | !nzchar(name), name := "."]
  bed_dt[, chrom_order := genome_dt$chrom_order[match(chrom, genome_dt$chrom)]]
  setorder(bed_dt, chrom_order, start, end, name)
  bed_dt[, chrom_order := NULL]
  bed_dt
}

write_bed <- function(dt, path) {
  fwrite(dt, file = path, sep = "\t", col.names = FALSE, quote = FALSE)
}

run_make_snp_random_genic_windows <- function(config = make_snp_random_genic_windows_config) {
  required_names <- c("snp_csv", "gff", "out_dir", "window_bp", "promoter_bp", "seed")
  missing_names <- setdiff(required_names, names(config))
  if (length(missing_names)) {
    stop("Config is missing required entries: ", paste(missing_names, collapse = ", "))
  }

  config <- modifyList(make_snp_random_genic_windows_config, config)

  if (is.na(config$snp_csv) || !nzchar(config$snp_csv)) {
    stop("config$snp_csv must point to a SNP table.")
  }
  if (is.na(config$window_bp) || config$window_bp <= 0L) {
    stop("config$window_bp must be a positive integer.")
  }
  if (is.na(config$promoter_bp) || config$promoter_bp <= 0L) {
    stop("config$promoter_bp must be a positive integer.")
  }
  if (is.na(config$seed)) {
    stop("config$seed must be an integer.")
  }

  dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

  annotation <- load_annotation(config$gff)
  genome_dt <- annotation$genome
  gene_gr <- annotation$genes

  snp_dt <- read_snp_table(
    path = config$snp_csv,
    genome_dt = genome_dt,
    chr_col = config$chr_col,
    pos_col = config$pos_col,
    id_col = config$id_col,
    variant_col = config$variant_col
  )

  snp_window_gr <- make_snp_windows(snp_dt, genome_dt, as.integer(config$window_bp))
  random_window_gr <- make_random_windows(length(snp_window_gr), genome_dt, as.integer(config$window_bp), as.integer(config$seed))
  gene_outputs <- make_gene_outputs(gene_gr, genome_dt, as.integer(config$promoter_bp))

  snp_bed <- gr_to_bed(snp_window_gr, genome_dt)
  random_bed <- gr_to_bed(random_window_gr, genome_dt)
  genes_bed <- gr_to_bed(gene_outputs$genes, genome_dt)
  promoters_bed <- gr_to_bed(gene_outputs$promoters, genome_dt)

  snp_path <- file.path(config$out_dir, "snp_windows.bed")
  random_path <- file.path(config$out_dir, "random_windows.bed")
  genes_path <- file.path(config$out_dir, "genic_regions.bed")
  promoters_path <- file.path(config$out_dir, "promoter_regions.bed")

  write_bed(snp_bed, snp_path)
  write_bed(random_bed, random_path)
  write_bed(genes_bed, genes_path)
  write_bed(promoters_bed, promoters_path)

  message("Wrote SNP windows: ", nrow(snp_bed), " -> ", snp_path)
  message("Wrote random windows: ", nrow(random_bed), " -> ", random_path)
  message("Wrote genic regions: ", nrow(genes_bed), " -> ", genes_path)
  message("Wrote promoter regions: ", nrow(promoters_bed), " -> ", promoters_path)

  invisible(list(
    config = config,
    snp_table = snp_dt,
    snp_bed = snp_bed,
    random_bed = random_bed,
    genes_bed = genes_bed,
    promoters_bed = promoters_bed,
    output_paths = c(
      snp_windows = snp_path,
      random_windows = random_path,
      genic_regions = genes_path,
      promoter_regions = promoters_path
    )
  ))
}

# Source this file in an R session, then run:
# result <- run_make_snp_random_genic_windows()
