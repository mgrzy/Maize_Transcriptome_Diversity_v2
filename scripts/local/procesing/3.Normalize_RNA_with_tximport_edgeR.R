# Normalize the NE2020 RNA-seq gene-count matrix and summarize expression
# variance.
#
# Workflow:
#   1. Read the raw gene-by-sample count matrix and genotype line IDs.
#   2. Match RNA-seq samples to genotype lines and discard unmatched samples.
#   3. Pass the counts through tximport's "none" importer for a consistent
#      count object.
#   4. Apply edgeR::filterByExpr as the first gene-filtering step.
#   5. Normalize the retained genes with edgeR TMM and save log2-CPM values.
#   6. Normalize the same retained genes with DESeq2 variance stabilizing
#      transformation (VST) and save the VST values.
#   7. Aggregate replicate RNA profiles to line-level means.
#   8. Save gene-wise mean and standard deviation for raw, edgeR, and DESeq2
#      line-level expression in results/full/GeneExpressionVariance/.
#
# Raw expression summaries include all genes. The edgeR and DESeq2 summaries
# include only genes retained by filterByExpr. Missing genotype lines remain in
# the output matrices as rows with missing values.

suppressPackageStartupMessages({
  # edgeR: expression filtering, TMM normalization, and log2-CPM values.
  # DESeq2: variance stabilizing transformation.
  # tximport: standardized import of the temporary per-sample count files.
  library(edgeR)
  library(DESeq2)
  library(tximport)
})

# -----------------------------------------------------------------------------
# 1. Input and output paths
# -----------------------------------------------------------------------------
counts_path <- "data/raw/rna/counts_NE2020_corrected_names.csv.gz"
geno_lines_path <- "data/processed/geno/plink/Div639.geno.ind"
output_dir <- "data/processed/rna"
output_path <- file.path(
  output_dir,
  "counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz"
)
deseq2_vst_output_path <- file.path(
  output_dir,
  "counts_NE2020_corrected_names_deseq2_vst_filtered.csv.gz"
)
# Variance summaries are project-level results, separate from the normalized
# expression matrices written under data/processed/rna.
variance_dir <- "results/full/GeneExpressionVariance"

# -----------------------------------------------------------------------------
# 2. Helper functions for reading, matching, and preparing counts
# -----------------------------------------------------------------------------
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
  # Generate matching candidates for common sample-name suffix/prefix variants.
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
  # tximport expects one input file per sample. These temporary files preserve
  # the original counts, with a dummy transcript length of one.
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
  # Collapse technical/biological replicate profiles by taking their mean for
  # each genotype line, then restore the requested genotype-line order.
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

write_expression_matrix <- function(expression_mat, path) {
  # Save a line-by-gene matrix with the line identifier in the first column.
  expression_df <- data.frame(
    Line = rownames(expression_mat),
    expression_mat,
    check.names = FALSE
  )

  write.csv(
    expression_df,
    file = gzfile(path),
    row.names = FALSE,
    quote = FALSE
  )
}

write_gene_expression_variance <- function(expression_mat, path) {
  # Calculate one mean and SD per gene across genotype lines.
  gene_mean <- colMeans(expression_mat, na.rm = TRUE)
  gene_sd <- apply(expression_mat, 2L, sd, na.rm = TRUE)

  # Genes with no observed values should remain explicitly missing.
  no_values <- colSums(!is.na(expression_mat)) == 0L
  gene_mean[no_values] <- NA_real_
  gene_sd[no_values] <- NA_real_

  variance_df <- data.frame(
    Gene = colnames(expression_mat),
    Mean = gene_mean,
    SD = gene_sd,
    check.names = FALSE
  )

  write.csv(
    variance_df,
    file = gzfile(path),
    row.names = FALSE,
    quote = FALSE
  )
}

# -----------------------------------------------------------------------------
# 3. Create output directories and read the input data
# -----------------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(variance_dir, recursive = TRUE, showWarnings = FALSE)

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

# -----------------------------------------------------------------------------
# 4. Keep only samples that can be matched to genotype lines
# -----------------------------------------------------------------------------
keep_samples <- !is.na(sample_to_line)
count_mat <- count_mat[, keep_samples, drop = FALSE]
sample_to_line <- sample_to_line[keep_samples]

# -----------------------------------------------------------------------------
# 5. Create temporary tximport inputs and import the count matrix
# -----------------------------------------------------------------------------
tmp_dir <- tempfile(pattern = "tximport_none_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

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

# The script is commonly run with source(), so clean up only after tximport has
# finished reading the temporary files. A top-level on.exit() would remove the
# directory before the next source() expression is evaluated.
unlink(tmp_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 6. Filter genes with filterByExpr (first filtering step), then apply TMM
#    normalization with edgeR
# -----------------------------------------------------------------------------
design <- model.matrix(~1, data = data.frame(sample = colnames(txi$counts)))
y <- DGEList(counts = txi$counts)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# -----------------------------------------------------------------------------
# 7. Save the existing edgeR log2-CPM normalization at genotype-line level
# -----------------------------------------------------------------------------
log2_normalized <- cpm(y, log = TRUE, prior.count = 2)
line_by_gene <- aggregate_lines(
  log2_mat = log2_normalized,
  line_ids = sample_to_line[colnames(log2_normalized)],
  line_order = geno_lines
)

# Keep the existing edgeR output unchanged.
write_expression_matrix(line_by_gene, output_path)

# -----------------------------------------------------------------------------
# 8. Apply DESeq2 variance stabilizing transformation to the same filtered
#    genes and save the VST values at genotype-line level
# -----------------------------------------------------------------------------
# Counts are integer-valued for this tximport input, but round explicitly for
# DESeq2's count-matrix contract.
# y is an edgeR DGEList, so access its count matrix through $counts rather
# than the DESeq2 counts() accessor.
deseq2_counts <- round(y$counts)
deseq2_col_data <- data.frame(
  sample = colnames(deseq2_counts),
  row.names = colnames(deseq2_counts)
)
deseq2_dds <- DESeqDataSetFromMatrix(
  countData = deseq2_counts,
  colData = deseq2_col_data,
  design = ~1
)
deseq2_vst <- assay(vst(deseq2_dds, blind = TRUE))
deseq2_vst_line_by_gene <- aggregate_lines(
  log2_mat = deseq2_vst,
  line_ids = sample_to_line[colnames(deseq2_vst)],
  line_order = geno_lines
)
write_expression_matrix(deseq2_vst_line_by_gene, deseq2_vst_output_path)

# -----------------------------------------------------------------------------
# 9. Calculate and save gene-wise mean and SD for raw and normalized data
# -----------------------------------------------------------------------------
# Summaries are calculated on the same line-level matrices written above. Raw
# expression is retained without the filterByExpr gene selection; normalized
# summaries use the filtered genes.
raw_line_by_gene <- aggregate_lines(
  log2_mat = count_mat,
  line_ids = sample_to_line[colnames(count_mat)],
  line_order = geno_lines
)
write_gene_expression_variance(
  raw_line_by_gene,
  file.path(variance_dir, "raw_expression_mean_sd.csv.gz")
)
write_gene_expression_variance(
  line_by_gene,
  file.path(variance_dir, "edgeR_log2_cpm_filtered_mean_sd.csv.gz")
)
write_gene_expression_variance(
  deseq2_vst_line_by_gene,
  file.path(variance_dir, "DESeq2_vst_filtered_mean_sd.csv.gz")
)

message("Saved filtered log2 CPM matrix to ", output_path)
message("Saved DESeq2 VST matrix to ", deseq2_vst_output_path)
message("Saved mean/SD summaries to ", variance_dir)
