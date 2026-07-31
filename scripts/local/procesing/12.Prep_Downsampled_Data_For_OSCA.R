### Script: scripts/local/procesing/12.Prep_Downsampled_Data_For_OSCA.R
### Purpose: Prepare OSCA expression, probe, individual, and PC inputs for
###          each genotype downsample. Genes are never filtered: every gene
###          already present in the full OSCA expression file is retained.
### Inputs:
###  - data/processed/osca/WiDiv631.dat
###  - data/processed/osca/WiDiv631_OPI.txt
###  - data/processed/geno/plink/downsample/DownSampling_<n>.fam
###  - results/full/PCA/downsample/DownSampling_<n>.eigenvec
### Outputs (per n):
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>.dat
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>_OPI.txt
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>_OII.txt
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>_PC.txt
###    (five genetic PCs followed by five expression PCs)
### Usage:
###  Rscript scripts/local/procesing/12.Prep_Downsampled_Data_For_OSCA.R
###  Rscript scripts/local/procesing/12.Prep_Downsampled_Data_For_OSCA.R 100 200
### Dependencies: data.table

suppressPackageStartupMessages(library(data.table))

# Resolve project-relative paths even when the script is launched outside the
# repository root.
args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)])
if (length(script_file) == 1L) {
  setwd(normalizePath(file.path(dirname(script_file), "../../..")))
}

targets <- if (length(args) == 0L) {
  c(100L, 150L, 200L, 300L, 400L, 500L, 550L, 600L)
} else {
  suppressWarnings(as.integer(args))
}

if (anyNA(targets) || any(!targets %in% c(100L, 150L, 200L, 300L, 400L, 500L, 550L, 600L))) {
  stop("Optional arguments must be one or more supported downsample sizes: ",
       "100, 150, 200, 300, 400, 500, 550, 600.")
}

expression_file <- "data/processed/osca/WiDiv631.dat"
opi_file <- "data/processed/osca/WiDiv631_OPI.txt"
fam_dir <- "data/processed/geno/plink/downsample"
pca_dir <- "results/full/PCA/downsample"
output_root <- "data/processed/osca/downsample"

# These are the complete, pre-established OSCA expression genes.  Do not
# subset columns by expression level, annotation, or downsample membership.
expression <- fread(expression_file, data.table = FALSE, check.names = FALSE)
names(expression)[1] <- "Taxa"
expression$Taxa <- as.character(expression$Taxa)
gene_columns <- names(expression)[-1]

if (anyDuplicated(expression$Taxa)) {
  stop("The full OSCA expression file contains duplicate Taxa IDs.")
}

opi <- fread(opi_file, header = FALSE, data.table = FALSE)
if (ncol(opi) != 5L || anyDuplicated(opi[[2]])) {
  stop("OPI must contain five columns and unique probe IDs.")
}
if (!setequal(as.character(opi[[2]]), gene_columns)) {
  stop("The OPI probes do not exactly match the full OSCA expression gene ",
       "set; refusing to alter the gene set.")
}
# OSCA requires OPI rows to correspond to expression columns. This reorders
# metadata only; all expression genes remain present and in their source order.
opi <- opi[match(gene_columns, opi[[2]]), , drop = FALSE]

for (target in targets) {
  prefix <- sprintf("DownSampling_%d", target)
  fam_file <- file.path(fam_dir, paste0(prefix, ".fam"))
  pca_file <- file.path(pca_dir, paste0(prefix, ".eigenvec"))
  output_dir <- file.path(output_root, prefix)

  if (!file.exists(fam_file)) {
    stop("Missing downsample FAM file: ", fam_file)
  }
  if (!file.exists(pca_file)) {
    stop("Missing downsample PCA file: ", pca_file)
  }

  fam <- fread(fam_file, header = FALSE, data.table = FALSE)
  if (ncol(fam) < 2L || nrow(fam) != target) {
    stop("Expected ", target, " samples in FAM file: ", fam_file)
  }
  fid <- as.character(fam[[1]])
  iid <- as.character(fam[[2]])
  if (anyDuplicated(iid)) {
    stop("Duplicate IID values in FAM file: ", fam_file)
  }

  sample_index <- match(iid, expression$Taxa)
  if (anyNA(sample_index)) {
    stop("Expression is missing these downsampled samples: ",
         paste(iid[is.na(sample_index)], collapse = ", "))
  }

  # Keep all existing expression columns, in their original order; only rows
  # are selected and reordered to the genotype FAM order.
  expression_subset <- expression[sample_index, , drop = FALSE]
  expression_subset$Taxa <- iid

  pca <- fread(pca_file, header = FALSE, data.table = FALSE)
  if (ncol(pca) < 7L) {
    stop("Expected FID, IID, and at least five PCs in: ", pca_file)
  }
  names(pca)[1:7] <- c("FID", "IID", paste0("PC", 1:5))
  pca$IID <- as.character(pca$IID)
  if (anyDuplicated(pca$IID)) {
    stop("Duplicate IID values in PCA file: ", pca_file)
  }

  pca_index <- match(iid, pca$IID)
  if (anyNA(pca_index)) {
    stop("PCA is missing these downsampled samples: ",
         paste(iid[is.na(pca_index)], collapse = ", "))
  }

  # Calculate expression PCs from the complete, unchanged gene set. Samples are
  # rows and genes are columns, so each PC describes sample-level expression
  # variation within this downsample.
  expression_pca <- prcomp(
    expression_subset[, gene_columns, drop = FALSE],
    center = TRUE,
    scale. = FALSE,
    rank. = 5L
  )
  if (ncol(expression_pca$x) < 5L) {
    stop("Could not calculate five expression PCs for ", prefix)
  }
  expression_pcs <- as.data.frame(expression_pca$x[, 1:5, drop = FALSE])
  names(expression_pcs) <- paste0("ExprPC", 1:5)

  # OSCA qcovar input: two ID columns followed by five genetic PCs and five
  # expression PCs. Rows have the exact same order as expression and OII files.
  genetic_pcs <- pca[pca_index, paste0("PC", 1:5), drop = FALSE]
  names(genetic_pcs) <- paste0("GenoPC", 1:5)
  pc_subset <- data.frame(
    Taxa = fid,
    Taxa2 = iid,
    genetic_pcs,
    expression_pcs,
    check.names = FALSE
  )
  oii <- data.frame(Taxa = iid, a = 0, b = 0, z = NA)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(oii,
              file.path(output_dir, paste0(prefix, "_OII.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(opi,
              file.path(output_dir, paste0(prefix, "_OPI.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(expression_subset,
              file.path(output_dir, paste0(prefix, ".dat")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  write.table(pc_subset,
              file.path(output_dir, paste0(prefix, "_PC.txt")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

  message(prefix, ": wrote ", nrow(expression_subset), " samples, ",
          length(gene_columns), " unchanged genes, five genetic PCs, and ",
          "five expression PCs to ", output_dir)
}
