### Script: scripts/local/procesing/14.Prep_Full_OSCA_10PC_Covariates.R
### Purpose: Create full-data OSCA covariates containing only genetic PCs 1-5
###          and expression PCs 1-5, using the same expression-PCA procedure
###          as the downsampled OSCA preparation.
### Inputs:
###  - data/processed/osca/WiDiv631.oii
###  - data/processed/osca/WiDiv631.dat
###  - results/full/PCA/PCA_Geno.eigenvec
### Output:
###  - data/processed/osca/WiDiv631_PC_5Geno_5Expr.txt
### Usage: Rscript scripts/local/procesing/14.Prep_Full_OSCA_10PC_Covariates.R
### Dependencies: data.table

suppressPackageStartupMessages(library(data.table))

# Resolve project-relative paths even when the script is launched outside the
# repository root.
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
if (length(script_file) == 1L) {
  setwd(normalizePath(file.path(dirname(script_file), "../../..")))
}

oii_file <- "data/processed/osca/WiDiv631.oii"
expression_file <- "data/processed/osca/WiDiv631.dat"
genetic_pca_file <- "results/full/PCA/PCA_Geno.eigenvec"
output_file <- "data/processed/osca/WiDiv631_PC_5Geno_5Expr.txt"

# The BOD individual-information file defines the order OSCA will use.
oii <- fread(oii_file, header = FALSE, data.table = FALSE)
if (ncol(oii) < 2L || anyDuplicated(oii[[2]])) {
  stop("OSCA OII file must contain unique individual IDs in column 2: ", oii_file)
}
fid <- as.character(oii[[1]])
iid <- as.character(oii[[2]])

# Read the first five genetic PCs from PLINK and align them to OII order.
genetic_pca <- fread(genetic_pca_file, header = FALSE, data.table = FALSE)
if (ncol(genetic_pca) < 7L) {
  stop("Expected FID, IID, and at least five genetic PCs in: ", genetic_pca_file)
}
names(genetic_pca)[1:7] <- c("FID", "IID", paste0("GenoPC", 1:5))
genetic_pca$IID <- as.character(genetic_pca$IID)
if (anyDuplicated(genetic_pca$IID)) {
  stop("Genetic PCA file contains duplicate individual IDs: ", genetic_pca_file)
}
genetic_index <- match(iid, genetic_pca$IID)
if (anyNA(genetic_index)) {
  stop("Genetic PCA is missing OII samples: ",
       paste(iid[is.na(genetic_index)], collapse = ", "))
}
genetic_pcs <- genetic_pca[genetic_index, paste0("GenoPC", 1:5), drop = FALSE]

# Use exactly the full OSCA expression matrix, retaining all existing genes.
# This is the same sample-wise, centered, unscaled PCA used for downsampling.
expression <- fread(expression_file, data.table = FALSE, check.names = FALSE)
names(expression)[1] <- "Taxa"
expression$Taxa <- as.character(expression$Taxa)
if (anyDuplicated(expression$Taxa)) {
  stop("Expression file contains duplicate sample IDs: ", expression_file)
}
expression_index <- match(iid, expression$Taxa)
if (anyNA(expression_index)) {
  stop("Expression data are missing OII samples: ",
       paste(iid[is.na(expression_index)], collapse = ", "))
}
gene_columns <- names(expression)[-1]
expression_subset <- expression[expression_index, gene_columns, drop = FALSE]

expression_pca <- prcomp(
  expression_subset,
  center = TRUE,
  scale. = FALSE,
  rank. = 5L
)
if (ncol(expression_pca$x) < 5L) {
  stop("Could not calculate five expression PCs from: ", expression_file)
}
expression_pcs <- as.data.frame(expression_pca$x[, 1:5, drop = FALSE])
names(expression_pcs) <- paste0("ExprPC", 1:5)

covariates <- data.frame(
  Taxa = fid,
  Taxa2 = iid,
  genetic_pcs,
  expression_pcs,
  check.names = FALSE
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(
  covariates,
  output_file,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

message("Wrote ", nrow(covariates), " samples with 5 genetic PCs and 5 expression PCs: ",
        output_file)
