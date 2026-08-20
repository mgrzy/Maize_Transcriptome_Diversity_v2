### Script: scripts/local/analysis/13.Summarize_Downsampled_OSCA_Significant.R
### Purpose: Summarize significant downsampled and full-data OSCA associations
###          by counting significant genes and significant SNPs per gene.
### Inputs:
###  - results/full/osca/downsample/DownSampling_<n>/significant*
###  - results/full/osca/full_5Geno_5ExprPC/significant*
### Outputs:
###  - results/full/osca/significant_summary/significant_gene_counts.tsv
###    (includes the mean number of significant SNPs per significant gene)
###  - results/full/osca/significant_summary/significant_snps_per_gene.tsv
### Usage: Rscript scripts/local/analysis/13.Summarize_Downsampled_OSCA_Significant.R
### Dependencies: data.table

suppressPackageStartupMessages(library(data.table))

# Resolve project-relative paths even when the script is launched outside the
# repository root.
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
if (length(script_file) == 1L) {
  setwd(normalizePath(file.path(dirname(script_file), "../../..")))
}

results_root <- "results/full/osca"
output_dir <- file.path(results_root, "significant_summary")

significant_files <- list.files(
  results_root,
  pattern = "^significant",
  recursive = TRUE,
  full.names = TRUE
)

# OSCA writes queried BESD results as extensionless text files in this setup.
# Include downsampled and full-data exports, but exclude accompanying logs.
significant_files <- significant_files[
  (grepl("/downsample/DownSampling_[0-9]+/significant", significant_files) |
     grepl("/full_5Geno_5ExprPC/significant", significant_files)) &
    !grepl("\\.log$", significant_files)
]

if (length(significant_files) == 0L) {
  stop("No significant OSCA result files found under: ", results_root)
}

summaries <- vector("list", length(significant_files))
gene_summaries <- vector("list", length(significant_files))

for (i in seq_along(significant_files)) {
  result_file <- significant_files[[i]]
  downsample <- if (grepl("/downsample/DownSampling_[0-9]+/", result_file)) {
    sub(".*/(DownSampling_[0-9]+)/significant.*$", "\\1", result_file)
  } else {
    "Full_5Geno_5ExprPC"
  }
  sample_size <- if (downsample == "Full_5Geno_5ExprPC") {
    631L
  } else {
    as.integer(sub("DownSampling_", "", downsample))
  }
  header <- names(fread(result_file, nrows = 0L, check.names = FALSE))

  required_columns <- "Gene"
  missing_columns <- setdiff(required_columns, header)
  if (length(missing_columns) > 0L) {
    stop("Missing required column(s) in ", result_file, ": ",
         paste(missing_columns, collapse = ", "))
  }

  # OSCA has one result row for each tested SNP-gene pair. Read only Gene: the
  # queried tables can be several GB, while this is sufficient to count SNPs
  # per gene and keeps memory use substantially lower.
  results <- fread(result_file, select = "Gene", check.names = FALSE)
  results <- results[!is.na(Gene) & nzchar(Gene)]

  snps_per_gene <- results[,
    .(SignificantSNPs = .N),
    by = Gene
  ]
  setorder(snps_per_gene, Gene)
  snps_per_gene[, `:=`(
    Downsample = downsample,
    SampleSize = sample_size,
    ResultFile = basename(result_file)
  )]
  setcolorder(snps_per_gene, c("Downsample", "SampleSize", "ResultFile", "Gene", "SignificantSNPs"))

  summaries[[i]] <- data.table(
    Downsample = downsample,
    SampleSize = sample_size,
    ResultFile = basename(result_file),
    SignificantGenes = nrow(snps_per_gene),
    SignificantSNPs = sum(snps_per_gene$SignificantSNPs),
    MeanSignificantSNPsPerGene = mean(snps_per_gene$SignificantSNPs),
    SignificantAssociations = nrow(results)
  )
  gene_summaries[[i]] <- snps_per_gene
  rm(results)
  gc(verbose = FALSE)
}

gene_counts <- rbindlist(summaries)
setorder(gene_counts, SampleSize, ResultFile)

snps_per_gene <- rbindlist(gene_summaries)
setorderv(
  snps_per_gene,
  c("SampleSize", "SignificantSNPs", "Gene"),
  c(1L, -1L, 1L),
  na.last = TRUE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(gene_counts, file.path(output_dir, "significant_gene_counts.tsv"), sep = "\t")
fwrite(snps_per_gene, file.path(output_dir, "significant_snps_per_gene.tsv"), sep = "\t")

message("Summarized ", nrow(gene_counts), " significant OSCA file(s).")
