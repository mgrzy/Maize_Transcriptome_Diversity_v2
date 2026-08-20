#!/usr/bin/env Rscript
# Per-gene narrow-sense heritability and eQTL/fine-mapped variance explained.
#
# Run: Rscript scripts/local/analysis/3.Narrow_sense_h2.R
# Test a small run: Rscript scripts/local/analysis/3.Narrow_sense_h2.R --max-genes 10 --cores 4
# Split a full run into ten independent serial jobs:
#   Rscript scripts/local/analysis/3.Narrow_sense_h2.R --job-index 1 --job-count 10 --cores 1
# Merge completed job shards:
#   Rscript scripts/local/analysis/3.Narrow_sense_h2.R --merge-only --job-count 10
# `sommer::mmer()` can stall in a forked R worker.  Fits therefore run serially
# by default; override only after testing your R/sommer installation, e.g.
# MAX_CONCURRENT_MMER=2 Rscript ...

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
  max_genes = Inf,
  job_index = 1L,
  job_count = 1L,
  merge_only = FALSE,
  max_concurrent_mmer = as.integer(Sys.getenv("MAX_CONCURRENT_MMER", "1"))
)

parse_args <- function(args, values) {
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      cat(paste0(
        "Usage: Rscript 3.Narrow_sense_h2.R [--cores N] [--max-genes N] ",
        "[--job-index N] [--job-count N] [--merge-only]\n"
      ))
      quit(status = 0L)
    }
    if (key == "--merge-only") {
      values$merge_only <- TRUE
      i <- i + 1L
      next
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
config$cores <- as.integer(config$cores)
config$max_genes <- as.numeric(config$max_genes)
config$job_index <- as.integer(config$job_index)
config$job_count <- as.integer(config$job_count)
config$max_concurrent_mmer <- as.integer(config$max_concurrent_mmer)
if (is.na(config$cores) || config$cores < 1L) stop("--cores must be a positive integer.")
if (is.na(config$max_genes) || config$max_genes < 1L) stop("--max-genes must be a positive number.")
if (is.na(config$job_count) || config$job_count < 1L) stop("--job-count must be a positive integer.")
if (is.na(config$job_index) || config$job_index < 1L || config$job_index > config$job_count) {
  stop("--job-index must be between 1 and --job-count.")
}
if (is.na(config$max_concurrent_mmer) || config$max_concurrent_mmer < 1L) {
  stop("MAX_CONCURRENT_MMER must be a positive integer.")
}

# Progress messages are intentionally emitted from both the parent and forked
# workers.  This makes long mixed-model fits visible when stdout/stderr is
# redirected to a log with `tee`.
progress_message <- function(...) {
  message(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] [pid ", Sys.getpid(), "] ",
    paste0(..., collapse = "")
  )
}

job_output_path <- function(output_dir, job_index, job_count) {
  file.path(
    output_dir,
    sprintf("narrow_sense_h2_and_r2.job_%03d_of_%03d.tsv.gz", job_index, job_count)
  )
}

merge_job_outputs <- function(output_dir, job_count) {
  shard_paths <- vapply(
    seq_len(job_count),
    function(job_index) job_output_path(output_dir, job_index, job_count),
    character(1L)
  )
  missing_shards <- shard_paths[!file.exists(shard_paths)]
  if (length(missing_shards)) {
    stop("Cannot merge: missing job output(s):\n", paste(missing_shards, collapse = "\n"))
  }

  progress_message("Reading and merging ", job_count, " job shard(s).")
  merged <- rbindlist(lapply(shard_paths, fread), fill = TRUE, use.names = TRUE)
  if (anyDuplicated(merged$gene)) {
    duplicate_genes <- unique(merged$gene[duplicated(merged$gene)])
    stop("Cannot merge: duplicate gene(s) across shards: ", paste(head(duplicate_genes, 10L), collapse = ", "))
  }
  setorder(merged, gene)

  output_path <- file.path(output_dir, "narrow_sense_h2_and_r2.tsv.gz")
  fwrite(merged, output_path, sep = "\t", compress = "gzip")
  progress_message("Wrote ", nrow(merged), " merged gene-level records to ", output_path)
}

dir.create(config$output_path, recursive = TRUE, showWarnings = FALSE)
if (config$merge_only) {
  merge_job_outputs(config$output_path, config$job_count)
  quit(status = 0L)
}

progress_message("Validating input files.")
required_paths <- unlist(config[c("ind_path", "kin_desc_path", "geno_desc_path", "map_path", "pc_path", "rna_path", "eqtl_path", "fm_path", "training_lines_path", "test_lines_path")])
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) stop("Required file(s) missing:\n", paste(missing_paths, collapse = "\n"))

safe_h2 <- function(y, ids, K) {
  keep <- is.finite(y) & !is.na(ids)
  y <- y[keep]; ids <- ids[keep]
  if (length(y) < 20L || var(y) == 0) return(c(h2 = NA_real_, genetic_var = NA_real_, residual_var = NA_real_))
  K_sub <- K[ids, ids, drop = FALSE]
  d <- data.frame(id = factor(ids, levels = ids), y = y)
  # Only variance components are used below.  PEV calculation is unnecessary
  # here and is expensive when several per-gene fits run concurrently.
  fit <- tryCatch(
    mmer(
      y ~ 1,
      random = ~vs(id, Gu = K_sub),
      rcov = ~units,
      data = d,
      verbose = FALSE,
      getPEV = FALSE
    ),
    error = function(e) NULL
  )
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

progress_message("Loading genotype sample IDs and kinship matrix.")
individuals <- fread(config$ind_path, header = FALSE)[[1L]]
individuals <- as.character(individuals)
K <- as.matrix(attach.big.matrix(config$kin_desc_path))
if (nrow(K) != length(individuals)) stop("Kinship matrix and .ind file have different dimensions.")
rownames(K) <- colnames(K) <- individuals
progress_message("Loading genotype matrix and variant map.")
geno <- attach.big.matrix(config$geno_desc_path)
if (nrow(geno) != length(individuals)) stop("Genotype matrix and .ind file have different dimensions.")

map <- fread(config$map_path, select = "SNP")
if (nrow(map) != ncol(geno)) stop("Genotype map and matrix have different numbers of variants.")
variant_index <- setNames(seq_len(nrow(map)), map$SNP)

progress_message("Loading expression matrix and covariates.")
expression <- fread(config$rna_path)
setnames(expression, 1L, "Taxa")
expression[, Taxa := as.character(Taxa)]
pcs <- fread(config$pc_path)
setnames(pcs, 1L, "Taxa")
pcs[, Taxa := as.character(Taxa)]
covariate_names <- setdiff(names(pcs)[vapply(pcs, is.numeric, logical(1L))], "Taxa")
if (!length(covariate_names)) stop("No numeric PC/covariate columns found in ", config$pc_path)

progress_message("Aligning samples and loading eQTL/fine-mapping results.")
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
total_genes <- length(genes)
genes <- genes[((seq_along(genes) - 1L) %% config$job_count) == (config$job_index - 1L)]
if (!length(genes)) {
  stop("Job ", config$job_index, " of ", config$job_count, " was assigned no genes.")
}
expression_index <- match(sample_ids, expression$Taxa)
progress_message(
  "Job ", config$job_index, " of ", config$job_count, " received ", length(genes),
  " of ", total_genes, " eligible gene(s), with ", length(train_ids),
  " training and ", length(test_ids), " test samples."
)

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
    progress_message(target_gene, ": reading expression and fitting raw heritability.")
    y <- as.numeric(expression[[target_gene]][expression_index])
    if (any(!is.finite(y))) stop("non-finite expression values")
    y_adjusted_full <- residualise(y, covariates)
    y_adjusted_train <- residualise(y, covariates, train_index, train_index)
    y_adjusted_test <- residualise(y, covariates, train_index, test_index)
    raw_h2_start <- Sys.time()
    h_raw_full <- safe_h2(y, sample_ids, K)
    progress_message(
      target_gene, ": raw heritability complete in ",
      round(as.numeric(difftime(Sys.time(), raw_h2_start, units = "secs")), 1L), " s."
    )
    progress_message(target_gene, ": fitting PC-adjusted heritability.")
    adjusted_h2_start <- Sys.time()
    h_adjusted_full <- safe_h2(y_adjusted_full, sample_ids, K)
    progress_message(
      target_gene, ": PC-adjusted heritability complete in ",
      round(as.numeric(difftime(Sys.time(), adjusted_h2_start, units = "secs")), 1L), " s."
    )
    fit_variants <- function(snps) {
      X <- extract_genotypes(snps)
      fit_r2(y_adjusted_full, X, y_adjusted_train, X[train_index, , drop = FALSE], y_adjusted_test, X[test_index, , drop = FALSE])
    }
    progress_message(target_gene, ": fitting OSCA eQTL model.")
    osca <- fit_variants(eqtl[eqtl$gene == target_gene, SNP])
    progress_message(target_gene, ": fitting SuSiE fine-mapped model.")
    susie <- fit_variants(fine_map[fine_map$gene == target_gene, SNP])
    result <- data.table(gene = target_gene, n_all = length(sample_ids), n_training = length(train_ids), n_test = length(test_ids),
      h2_raw = h_raw_full["h2"], h2_pc_adjusted = h_adjusted_full["h2"],
      genetic_var_raw = h_raw_full["genetic_var"], residual_var_raw = h_raw_full["residual_var"],
      genetic_var_pc_adjusted = h_adjusted_full["genetic_var"], residual_var_pc_adjusted = h_adjusted_full["residual_var"],
      osca_n_variants = osca["n_variants"], osca_r2_full = osca["r2_full"], osca_r2_training = osca["r2_train"], osca_rmse_training = osca["rmse_train"], osca_r2_test = osca["r2_test"],
      susie_n_variants = susie["n_variants"], susie_r2_full = susie["r2_full"], susie_r2_training = susie["r2_train"], susie_rmse_training = susie["rmse_train"], susie_r2_test = susie["r2_test"], error = NA_character_)
    progress_message(target_gene, ": complete.")
    result
  }, error = function(e) {
    progress_message(target_gene, ": ERROR: ", conditionMessage(e))
    data.table(gene = target_gene, error = conditionMessage(e))
  })
}

worker_cores <- min(config$cores, config$max_concurrent_mmer, length(genes))
if (worker_cores < config$cores) {
  progress_message(
    "Limiting concurrent sommer fits to ", worker_cores,
    " (requested --cores ", config$cores,
    "; override with MAX_CONCURRENT_MMER if validated for this host)."
  )
}
progress_message("Analysing ", length(genes), " genes using ", worker_cores, " concurrent worker(s).")
results <- if (.Platform$OS.type == "unix" && worker_cores > 1L) {
  mclapply(genes, analyse_gene, mc.cores = worker_cores, mc.preschedule = FALSE, mc.silent = FALSE)
} else lapply(genes, analyse_gene)
results <- rbindlist(results, fill = TRUE)
setorder(results, gene)
output_path <- if (config$job_count > 1L) {
  job_output_path(config$output_path, config$job_index, config$job_count)
} else {
  file.path(config$output_path, "narrow_sense_h2_and_r2.tsv.gz")
}
progress_message("Writing ", nrow(results), " result row(s).")
fwrite(results, output_path, sep = "\t", compress = "gzip")
progress_message("Wrote ", nrow(results), " gene-level records to ", output_path)
