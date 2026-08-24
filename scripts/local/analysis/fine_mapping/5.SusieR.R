#!/usr/bin/env Rscript
# Fine-map one or many cis-eQTL regions with SuSiE.
#
# Batch example:
#   Rscript 5.SusieR.R --genes-file batch_0000.txt --batch-id 0000
# List eligible genes:
#   Rscript 5.SusieR.R --list-genes

suppressPackageStartupMessages({ library(bigmemory); library(data.table); library(susieR) })

script_path <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))
project_root <- normalizePath(file.path(dirname(script_path), "../../.."))
defaults <- list(
  gene = NULL, genes_file = NULL, batch_id = NULL, list_genes = FALSE,
  geno_desc = file.path(project_root, "data/processed/geno/rMVP/Div639.geno.desc"),
  geno_map = file.path(project_root, "data/processed/geno/rMVP/Div639.geno.map"),
  geno_ind = file.path(project_root, "data/processed/geno/rMVP/Div639.geno.ind"),
  expression = file.path(project_root, "data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz"),
  covariates = file.path(project_root, "data/processed/osca/WiDiv631_PC.txt"),
  tss = file.path(project_root, "data/processed/osca/WiDiv631_OPI.txt"),
  output_dir = file.path(project_root, "results/full/SusieR"),
  window_bp = 2000000L, L = 10L, max_iter = 500L, seed = 1L
)
usage <- function() cat(
  "Usage: Rscript 5.SusieR.R --genes-file FILE --batch-id ID [options]\n",
  "       Rscript 5.SusieR.R --gene GENE [options]\n",
  "       Rscript 5.SusieR.R --list-genes [options]\n\n",
  "Batch output: PIP_batch_<ID>.tsv.gz, CS_batch_<ID>.tsv.gz, errors_batch_<ID>.tsv.gz\n", sep = "")
parse_args <- function(args, values) {
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) { usage(); quit(status = 0L) }
    if (key == "--list-genes") { values$list_genes <- TRUE; i <- i + 1L; next }
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- gsub("-", "_", substring(key, 3L))
    if (!name %in% names(values)) stop("Unknown option: ", key)
    if (i == length(args)) stop("Missing value for ", key)
    values[[name]] <- args[[i + 1L]]; i <- i + 2L
  }
  values
}
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
for (name in c("window_bp", "L", "max_iter", "seed")) opt[[name]] <- as.integer(opt[[name]])
if (anyNA(unlist(opt[c("window_bp", "L", "max_iter", "seed")])) || opt$window_bp < 1L || opt$L < 1L || opt$max_iter < 1L) stop("--window-bp, --L and --max-iter must be positive integers.")
require_file <- function(path) if (!file.exists(path)) stop("Required input does not exist: ", path)
invisible(lapply(unlist(opt[c("geno_desc", "geno_map", "geno_ind", "expression", "covariates", "tss")]), require_file))
normalise_names <- function(x) gsub("[^a-z0-9]", "", tolower(x))
column_by_name <- function(x, choices, label) {
  hit <- match(normalise_names(choices), normalise_names(names(x)), nomatch = 0L)
  hit <- hit[!is.na(hit) & hit > 0L]
  if (!length(hit)) stop("Cannot find ", label, " column. Available columns: ", paste(names(x), collapse = ", "))
  names(x)[hit[1L]]
}

# Load shared inputs once per batch.
expression <- fread(opt$expression, data.table = FALSE, check.names = FALSE)
names(expression)[1L] <- "Taxa"
tss <- fread(opt$tss, header = FALSE, data.table = FALSE)[, 1:3, drop = FALSE]
names(tss) <- c("chr", "gene", "tss")
tss$chr <- as.character(tss$chr); tss$gene <- as.character(tss$gene); tss$tss <- suppressWarnings(as.numeric(tss$tss))
tss <- tss[!is.na(tss$tss) & !duplicated(tss$gene), , drop = FALSE]
available_genes <- intersect(names(expression)[-1L], tss$gene)
if (opt$list_genes) { writeLines(sort(available_genes)); quit(status = 0L) }
if (!is.null(opt$gene) && !is.null(opt$genes_file)) stop("Use either --gene or --genes-file, not both.")
if (!is.null(opt$genes_file)) { require_file(opt$genes_file); genes <- trimws(readLines(opt$genes_file, warn = FALSE)) } else if (!is.null(opt$gene)) { genes <- opt$gene } else { usage(); stop("--gene or --genes-file is required.") }
genes <- unique(genes[nzchar(genes)])
unknown_genes <- setdiff(genes, available_genes)
if (length(unknown_genes)) warning("Skipping ", length(unknown_genes), " genes absent from expression or TSS data.")
genes <- intersect(genes, available_genes)
if (!length(genes)) stop("No requested genes are present in both expression and TSS data.")

map_header <- fread(opt$geno_map, nrows = 0L, data.table = FALSE, check.names = FALSE)
snp_col <- column_by_name(map_header, c("SNP", "marker", "id"), "SNP")
chr_col <- column_by_name(map_header, c("CHROM", "CHR", "chromosome"), "chromosome")
pos_col <- column_by_name(map_header, c("POS", "BP", "position"), "position")
map <- fread(opt$geno_map, select = unique(c(snp_col, chr_col, pos_col)), data.table = FALSE, check.names = FALSE)
map[[chr_col]] <- as.character(map[[chr_col]]); map[[pos_col]] <- suppressWarnings(as.numeric(map[[pos_col]]))
individuals <- fread(opt$geno_ind, header = FALSE, data.table = FALSE); names(individuals)[1L] <- "Taxa"
covariates <- fread(opt$covariates, data.table = FALSE, check.names = FALSE); names(covariates)[1L] <- "Taxa"
covariate_columns <- setdiff(names(covariates)[vapply(covariates, is.numeric, logical(1L))], "Taxa")
cov_idx <- match(individuals$Taxa, covariates$Taxa)
Z <- if (length(covariate_columns)) as.matrix(covariates[cov_idx, covariate_columns, drop = FALSE]) else matrix(nrow = nrow(individuals), ncol = 0L)
expression_idx <- match(individuals$Taxa, expression$Taxa)
keep_covariates <- if (ncol(Z)) complete.cases(Z) else rep(TRUE, nrow(individuals))
snp <- attach.big.matrix(opt$geno_desc)
if (nrow(snp) != nrow(individuals)) stop("Genotype matrix and .ind file have different numbers of samples.")

run_gene <- function(gene, index) {
  gene_tss <- tss[tss$gene == gene, , drop = FALSE]
  variant_idx <- which(map[[chr_col]] == gene_tss$chr & abs(map[[pos_col]] - gene_tss$tss) <= opt$window_bp)
  if (length(variant_idx) < 2L) stop("Fewer than two variants in the cis window")
  y <- suppressWarnings(as.numeric(expression[[gene]][expression_idx]))
  keep <- !is.na(y) & keep_covariates
  if (sum(keep) < 10L) stop("Fewer than 10 samples with complete phenotype and covariates")
  X <- snp[keep, variant_idx, drop = FALSE]
  storage.mode(X) <- "double"; colnames(X) <- as.character(map[[snp_col]][variant_idx])
  if (anyNA(X)) stop("Missing genotype calls found; this dataset is expected to be complete")
  X <- X[, apply(X, 2L, function(x) var(x) > 0), drop = FALSE]
  if (ncol(X) < 2L) stop("Fewer than two variable variants in the cis window")
  Z_gene <- cbind(`(Intercept)` = 1, Z[keep, , drop = FALSE])
  qr_Z <- qr(Z_gene)
  y <- qr.resid(qr_Z, y[keep]); X <- qr.resid(qr_Z, X)
  if (!is.finite(var(y)) || var(y) == 0) stop("Residual expression has zero variance")
  set.seed(opt$seed + index)
  fit <- susie(X, y, L = min(opt$L, ncol(X)), max_iter = opt$max_iter, estimate_residual_variance = TRUE, coverage = 0.95)
  sets <- fit$sets$cs
  if (is.null(sets) || !length(sets)) return(list(pip = NULL, cs = NULL, n = ncol(X), ncs = 0L, converged = isTRUE(fit$converged)))
  members <- data.table(gene = gene, cs = rep(seq_along(sets), lengths(sets)), variable = unlist(sets, use.names = FALSE))
  members[, `:=`(variable_prob = fit$pip[variable], SNP = colnames(X)[variable])]
  purity <- fit$sets$purity
  if (is.null(purity)) purity <- matrix(NA_real_, nrow = length(sets), ncol = 3L)
  credible_sets <- data.table(gene = gene, cs = seq_along(sets), n_variants = lengths(sets), min_abs_corr = purity[, 1L], mean_abs_corr = purity[, 2L], median_abs_corr = purity[, 3L])
  list(pip = members, cs = credible_sets, n = ncol(X), ncs = length(sets), converged = isTRUE(fit$converged))
}

pip_results <- vector("list", length(genes)); cs_results <- vector("list", length(genes)); errors <- vector("list", length(genes))
for (i in seq_along(genes)) {
  gene <- genes[[i]]
  result <- tryCatch(run_gene(gene, i), error = function(e) e)
  if (inherits(result, "error")) {
    errors[[i]] <- data.table(gene = gene, error = conditionMessage(result)); message(gene, ": FAILED - ", conditionMessage(result))
  } else {
    pip_results[[i]] <- result$pip; cs_results[[i]] <- result$cs
    message(sprintf("%s: %d variants, %d credible sets%s", gene, result$n, result$ncs, if (result$converged) "" else " (not converged)"))
  }
}
empty_pip <- data.table(gene = character(), cs = integer(), variable = integer(), variable_prob = numeric(), SNP = character())
empty_cs <- data.table(gene = character(), cs = integer(), n_variants = integer(), min_abs_corr = numeric(), mean_abs_corr = numeric(), median_abs_corr = numeric())
empty_errors <- data.table(gene = character(), error = character())
pip_results <- Filter(Negate(is.null), pip_results); cs_results <- Filter(Negate(is.null), cs_results); errors <- Filter(Negate(is.null), errors)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- if (!is.null(opt$batch_id)) paste0("_batch_", opt$batch_id) else paste0("_", gsub("[^[:alnum:]_.-]", "_", genes[[1L]]))
fwrite(if (length(pip_results)) rbindlist(pip_results) else empty_pip, file.path(opt$output_dir, paste0("PIP", suffix, ".tsv.gz")), sep = "\t", compress = "gzip")
fwrite(if (length(cs_results)) rbindlist(cs_results) else empty_cs, file.path(opt$output_dir, paste0("CS", suffix, ".tsv.gz")), sep = "\t", compress = "gzip")
fwrite(if (length(errors)) rbindlist(errors) else empty_errors, file.path(opt$output_dir, paste0("errors", suffix, ".tsv.gz")), sep = "\t", compress = "gzip")
if (length(errors)) quit(status = 1L)
