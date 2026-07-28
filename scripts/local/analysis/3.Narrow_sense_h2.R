#!/usr/bin/env Rscript
# Per-gene narrow-sense heritability and eQTL/fine-mapped variance explained.
#
# Run: Rscript scripts/local/analysis/3.Narrow_sense_h2.R
# Test a small run: Rscript scripts/local/analysis/3.Narrow_sense_h2.R --max-genes 10 --cores 4

suppressPackageStartupMessages({
  library(bigmemory)
  library(data.table)
  library(parallel)
  library(sommer)
})

script_file_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(script_file_arg)) {
  script_path <- normalizePath(script_file_arg[[1L]])
} else {
  source_files <- vapply(rev(seq_len(sys.nframe())), function(i) {
    ofile <- sys.frame(i)$ofile
    if (is.null(ofile)) "" else as.character(ofile)[[1L]]
  }, character(1L))
  source_files <- source_files[nzchar(source_files)]
  if (!length(source_files)) stop("Cannot determine the script location; use Rscript or source a file.")
  script_path <- normalizePath(source_files[[1L]])
}
project_root <- normalizePath(file.path(dirname(script_path), "../../.."))
default_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

config <- list(
  ind_path = file.path(project_root, "data/processed/geno/rMVP/Div639.geno.ind"),
  kin_desc_path = file.path(project_root, "data/processed/geno/rMVP/Dib639_MAF005.kin.desc"),
  geno_desc_path = file.path(project_root, "data/processed/geno/rMVP/Div639.geno.desc"),
  map_path = file.path(project_root, "data/processed/geno/rMVP/Div639.geno.map"),
  pc_path = file.path(project_root, "data/processed/osca/WiDiv631_PC.txt"),
  rna_path = file.path(project_root, "data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz"),
  eqtl_path = file.path(project_root, "results/summary/osca/eQTL_P_min.txt.gz"),
  fm_path = file.path(project_root, "results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz"),
  training_lines_path = file.path(project_root, "data/processed/meta/training_lines.txt"),
  test_lines_path = file.path(project_root, "data/processed/meta/test_lines.txt"),
  output_path = file.path(project_root, "results/summary/Heritability_and_r2/"),
  cores = default_cores,
  max_genes = Inf
)

parse_args <- function(args, values) {
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      cat("Usage: Rscript 3.Narrow_sense_h2.R [--cores N] [--max-genes N]\n")
      quit(status = 0L)
    }
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- gsub("-", "_", substring(key, 3L))
    if (!name %in% names(values) || i == length(args)) stop("Invalid option: ", key)
    values[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values
}
config <- parse_args(commandArgs(trailingOnly = TRUE), config)
config$cores <- as.integer(config$cores); config$max_genes <- as.numeric(config$max_genes)
if (is.na(config$cores) || config$cores < 1L) stop("--cores must be a positive integer.")
if (is.na(config$max_genes) || config$max_genes < 1L) stop("--max-genes must be a positive number.")
required_paths <- unlist(config[c("ind_path", "kin_desc_path", "geno_desc_path", "map_path", "pc_path", "rna_path", "eqtl_path", "fm_path", "training_lines_path", "test_lines_path")])
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) stop("Required file(s) missing:\n", paste(missing_paths, collapse = "\n"))

safe_h2 <- function(y, ids, K) {
  keep <- is.finite(y) & !is.na(ids)
  y <- y[keep]; ids <- ids[keep]
  if (length(y) < 20L || var(y) == 0) return(c(h2 = NA_real_, genetic_var = NA_real_, residual_var = NA_real_))
  K_sub <- K[ids, ids, drop = FALSE]
  d <- data.frame(id = factor(ids, levels = ids), y = y)
  fit <- tryCatch(mmer(y ~ 1, random = ~vs(id, Gu = K_sub), rcov = ~units, data = d, verbose = FALSE), error = function(e) NULL)
  if (is.null(fit)) return(c(h2 = NA_real_, genetic_var = NA_real_, residual_var = NA_real_))
  genetic_var <- unname(fit$sigma[["u:id"]][1L]); residual_var <- unname(fit$sigma[["units"]][1L])
  total_var <- genetic_var + residual_var
  c(h2 = if (is.finite(total_var) && total_var > 0) genetic_var / total_var else NA_real_, genetic_var = genetic_var, residual_var = residual_var)
}

residualise <- function(y, covariates, train_index = seq_along(y), predict_index = seq_along(y)) {
  X_train <- cbind(`(Intercept)` = 1, covariates[train_index, , drop = FALSE])
  fit <- lm.fit(X_train, y[train_index])
  X_predict <- cbind(`(Intercept)` = 1, covariates[predict_index, , drop = FALSE])
  as.numeric(y[predict_index] - X_predict %*% fit$coefficients)
}

fit_r2 <- function(y_full, X_full, y_train, X_train, y_test, X_test) {
  result <- c(r2_full = NA_real_, r2_train = NA_real_, r2_test = NA_real_, rmse_train = NA_real_, n_variants = 0)
  if (!ncol(X_train)) return(result)
  variable <- vapply(seq_len(ncol(X_train)), function(j) var(X_train[, j], na.rm = TRUE) > 0, logical(1L))
  X_train <- X_train[, variable, drop = FALSE]; X_test <- X_test[, variable, drop = FALSE]
  X_full <- X_full[, variable, drop = FALSE]
  complete_full <- complete.cases(y_full, X_full); complete_train <- complete.cases(y_train, X_train); complete_test <- complete.cases(y_test, X_test)
  if (sum(complete_full) < ncol(X_train) + 3L || sum(complete_train) < ncol(X_train) + 3L || sum(complete_test) < 3L) return(result)
  X_full <- X_full[complete_full, , drop = FALSE]; y_full <- y_full[complete_full]
  X_train <- X_train[complete_train, , drop = FALSE]; y_train <- y_train[complete_train]
  X_test <- X_test[complete_test, , drop = FALSE]; y_test <- y_test[complete_test]
  fit_full <- lm.fit(cbind(`(Intercept)` = 1, X_full), y_full)
  pred_full <- as.numeric(cbind(`(Intercept)` = 1, X_full) %*% fit_full$coefficients)
  fit <- lm.fit(cbind(`(Intercept)` = 1, X_train), y_train)
  pred_train <- as.numeric(cbind(`(Intercept)` = 1, X_train) %*% fit$coefficients)
  pred_test <- as.numeric(cbind(`(Intercept)` = 1, X_test) %*% fit$coefficients)
  r2 <- function(y, prediction) {
    denominator <- sum((y - mean(y)) ^ 2)
    if (denominator > 0) 1 - sum((y - prediction) ^ 2) / denominator else NA_real_
  }
  result[c("r2_full", "r2_train", "r2_test")] <- c(r2(y_full, pred_full), r2(y_train, pred_train), r2(y_test, pred_test))
  result["rmse_train"] <- sqrt(mean((y_train - pred_train) ^ 2))
  result["n_variants"] <- ncol(X_train)
  result
}

individuals <- fread(config$ind_path, header = FALSE)[[1L]]
individuals <- as.character(individuals)
K <- as.matrix(attach.big.matrix(config$kin_desc_path))
if (nrow(K) != length(individuals)) stop("Kinship matrix and .ind file have different dimensions.")
rownames(K) <- colnames(K) <- individuals
geno <- attach.big.matrix(config$geno_desc_path)
if (nrow(geno) != length(individuals)) stop("Genotype matrix and .ind file have different dimensions.")

map <- fread(config$map_path, select = "SNP")
if (nrow(map) != ncol(geno)) stop("Genotype map and matrix have different numbers of variants.")
variant_index <- setNames(seq_len(nrow(map)), map$SNP)

expression <- fread(config$rna_path)
setnames(expression, 1L, "Taxa")
expression[, Taxa := as.character(Taxa)]
pcs <- fread(config$pc_path)
setnames(pcs, 1L, "Taxa")
pcs[, Taxa := as.character(Taxa)]
covariate_names <- setdiff(names(pcs)[vapply(pcs, is.numeric, logical(1L))], "Taxa")
if (!length(covariate_names)) stop("No numeric PC/covariate columns found in ", config$pc_path)

training_ids <- trimws(readLines(config$training_lines_path)); test_ids <- trimws(readLines(config$test_lines_path))
analysis_ids <- Reduce(intersect, list(individuals, expression$Taxa, pcs$Taxa, training_ids, test_ids))
if (length(analysis_ids)) stop("Training and test line lists overlap.")
analysis_ids <- intersect(individuals, intersect(expression$Taxa, pcs$Taxa))
train_ids <- intersect(analysis_ids, training_ids); test_ids <- intersect(analysis_ids, test_ids)
if (length(train_ids) < 20L || length(test_ids) < 3L) stop("Insufficient aligned training or test lines.")
sample_ids <- c(train_ids, test_ids)
sample_index <- match(sample_ids, individuals)
train_index <- seq_along(train_ids); test_index <- length(train_ids) + seq_along(test_ids)
pc_index <- match(sample_ids, pcs$Taxa)
covariates <- as.matrix(pcs[pc_index, ..covariate_names])
if (any(!complete.cases(covariates))) stop("PC covariates contain missing values for analysis lines.")

eqtl <- fread(config$eqtl_path, select = c("Gene", "SNP"))
eqtl[, gene := as.character(Gene)]; eqtl[, SNP := as.character(SNP)]
eqtl <- unique(eqtl[!is.na(gene) & !is.na(SNP), .(gene, SNP)])
fine_map <- fread(config$fm_path, select = c("gene", "SNP"))
fine_map[, `:=`(gene = as.character(gene), SNP = as.character(SNP))]
fine_map <- unique(fine_map[!is.na(gene) & !is.na(SNP)])

genes <- intersect(names(expression)[-1L], union(eqtl$gene, fine_map$gene))
genes <- head(sort(genes), config$max_genes)
if (!length(genes)) stop("No expression genes overlap the eQTL or fine-mapping results.")
expression_index <- match(sample_ids, expression$Taxa)

extract_genotypes <- function(snps) {
  cols <- unname(variant_index[unique(snps)])
  cols <- cols[!is.na(cols)]
  if (!length(cols)) return(matrix(numeric(), nrow = length(sample_index), ncol = 0L))
  X <- geno[sample_index, cols, drop = FALSE]
  storage.mode(X) <- "double"
  X
}

analyse_gene <- function(target_gene) {
  tryCatch({
    y <- as.numeric(expression[[target_gene]][expression_index])
    if (any(!is.finite(y))) stop("non-finite expression values")
    y_adjusted_full <- residualise(y, covariates)
    y_adjusted_train <- residualise(y, covariates, train_index, train_index)
    y_adjusted_test <- residualise(y, covariates, train_index, test_index)
    h_raw_full <- safe_h2(y, sample_ids, K)
    h_adjusted_full <- safe_h2(y_adjusted_full, sample_ids, K)
    fit_variants <- function(snps) {
      X <- extract_genotypes(snps)
      fit_r2(y_adjusted_full, X, y_adjusted_train, X[train_index, , drop = FALSE], y_adjusted_test, X[test_index, , drop = FALSE])
    }
    osca <- fit_variants(eqtl[eqtl$gene == target_gene, SNP])
    susie <- fit_variants(fine_map[fine_map$gene == target_gene, SNP])
    data.table(gene = target_gene, n_all = length(sample_ids), n_training = length(train_ids), n_test = length(test_ids),
      h2_raw = h_raw_full["h2"], h2_pc_adjusted = h_adjusted_full["h2"],
      genetic_var_raw = h_raw_full["genetic_var"], residual_var_raw = h_raw_full["residual_var"],
      genetic_var_pc_adjusted = h_adjusted_full["genetic_var"], residual_var_pc_adjusted = h_adjusted_full["residual_var"],
      osca_n_variants = osca["n_variants"], osca_r2_full = osca["r2_full"], osca_r2_training = osca["r2_train"], osca_rmse_training = osca["rmse_train"], osca_r2_test = osca["r2_test"],
      susie_n_variants = susie["n_variants"], susie_r2_full = susie["r2_full"], susie_r2_training = susie["r2_train"], susie_rmse_training = susie["rmse_train"], susie_r2_test = susie["r2_test"], error = NA_character_)
  }, error = function(e) data.table(gene = target_gene, error = conditionMessage(e)))
}

message("Analysing ", length(genes), " genes using ", config$cores, " core(s).")
results <- if (.Platform$OS.type == "unix" && config$cores > 1L) {
  mclapply(genes, analyse_gene, mc.cores = min(config$cores, length(genes)), mc.preschedule = FALSE)
} else lapply(genes, analyse_gene)
results <- rbindlist(results, fill = TRUE)
setorder(results, gene)
dir.create(config$output_path, recursive = TRUE, showWarnings = FALSE)
fwrite(results, file.path(config$output_path, "narrow_sense_h2_and_r2.tsv.gz"), sep = "\t", compress = "gzip")
message("Wrote ", nrow(results), " gene-level records to ", file.path(config$output_path, "narrow_sense_h2_and_r2.tsv.gz"))
