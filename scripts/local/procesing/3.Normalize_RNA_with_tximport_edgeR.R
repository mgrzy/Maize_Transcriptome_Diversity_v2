suppressPackageStartupMessages({
  library(edgeR)
  library(tximport)
})

counts_path <- "data/raw/rna/counts_NE2020_corrected_names.csv.gz"
geno_lines_path <- "data/processed/geno/rMVP/Div639.geno.ind"
output_dir <- "data/processed/rna"
output_path <- file.path(
  output_dir,
  "counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz"
)

read_counts_matrix <- function(path) {
  counts_df <- read.csv(path, check.names = FALSE)

  if (ncol(counts_df) < 2L) {
    stop("Expected a gene-by-sample matrix with at least two columns.")
  }

  gene_ids <- counts_df[[1]]
  count_mat <- as.matrix(counts_df[-1])
  rownames(count_mat) <- gene_ids
  storage.mode(count_mat) <- "double"
  count_mat
}

normalize_name <- function(x) {
  gsub("[^A-Za-z0-9]", "", toupper(x))
}

map_sample_to_line <- function(sample_name, geno_lines, geno_lines_normalized) {
  base_name <- sub("_R2$", "", sample_name)
  base_name_normalized <- normalize_name(base_name)
  trimmed_name_normalized <- sub("GOODMANBUCKLER$", "", base_name_normalized)

  candidates <- unique(c(
    base_name_normalized,
    trimmed_name_normalized,
    if (startsWith(base_name_normalized, "DK")) {
      substring(base_name_normalized, 3)
    },
    if (startsWith(trimmed_name_normalized, "DK")) {
      substring(trimmed_name_normalized, 3)
    }
  ))

  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- match(candidates, geno_lines_normalized)
  hit <- hit[!is.na(hit)]

  if (length(hit) == 0L) {
    return(NA_character_)
  }

  geno_lines[hit[1]]
}

write_tximport_inputs <- function(count_mat, tmp_dir) {
  sample_names <- colnames(count_mat)
  safe_names <- make.names(sample_names, unique = TRUE)
  files <- file.path(tmp_dir, paste0(safe_names, ".tsv"))
  names(files) <- sample_names

  for (i in seq_along(sample_names)) {
    sample_df <- data.frame(
      GeneID = rownames(count_mat),
      abundance = count_mat[, i],
      counts = count_mat[, i],
      length = 1,
      check.names = FALSE
    )

    write.table(
      sample_df,
      file = files[i],
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  files
}

aggregate_lines <- function(log2_mat, line_ids, line_order) {
  sample_by_gene <- t(log2_mat)
  line_sums <- rowsum(sample_by_gene, group = line_ids, reorder = FALSE)
  line_counts <- table(line_ids)

  # Average duplicated RNA profiles so the final matrix has one row per line.
  line_means <- sweep(
    line_sums,
    1,
    as.numeric(line_counts[rownames(line_sums)]),
    "/"
  )

  line_means[match(line_order, rownames(line_means)), , drop = FALSE]
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_mat <- read_counts_matrix(counts_path)
geno_lines <- scan(geno_lines_path, what = "character", quiet = TRUE)
geno_lines_normalized <- normalize_name(geno_lines)

sample_to_line <- vapply(
  colnames(count_mat),
  map_sample_to_line,
  character(1),
  geno_lines = geno_lines,
  geno_lines_normalized = geno_lines_normalized
)

keep_samples <- !is.na(sample_to_line)
count_mat <- count_mat[, keep_samples, drop = FALSE]
sample_to_line <- sample_to_line[keep_samples]

tmp_dir <- tempfile(pattern = "tximport_none_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

sample_files <- write_tximport_inputs(count_mat, tmp_dir)

txi <- tximport(
  files = sample_files,
  type = "none",
  txIn = TRUE,
  txOut = TRUE,
  txIdCol = "GeneID",
  abundanceCol = "abundance",
  countsCol = "counts",
  lengthCol = "length",
  importer = read.delim
)

design <- model.matrix(~1, data = data.frame(sample = colnames(txi$counts)))
y <- DGEList(counts = txi$counts)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

log2_normalized <- cpm(y, log = TRUE, prior.count = 2)
line_by_gene <- aggregate_lines(
  log2_mat = log2_normalized,
  line_ids = sample_to_line[colnames(log2_normalized)],
  line_order = geno_lines
)

line_by_gene_df <- data.frame(
  Line = rownames(line_by_gene),
  line_by_gene,
  check.names = FALSE
)

write.csv(
  line_by_gene_df,
  file = gzfile(output_path),
  row.names = FALSE,
  quote = FALSE
)

message("Saved filtered log2 CPM matrix to ", output_path)
