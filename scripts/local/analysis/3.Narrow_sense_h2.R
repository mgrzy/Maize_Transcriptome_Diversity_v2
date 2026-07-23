suppressPackageStartupMessages({
  library(bigmemory)
  library(data.table)
  library(parallel)
  library(sommer)
})

default_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

config <- list(
  ind_path = "data/processed/geno/rMVP/Div639.geno.ind",
  kin_desc_path = "data/processed/geno/rMVP/Dib639_MAF005.kin.desc",
  geno_desc_path = "data/processed/geno/rMVP/Div639.geno.desc",
  map_path = "data/processed/geno/rMVP/Div639.geno.map",
  pc_path = "results/full/covariateForEqlt/CovariateForEQTL.csv",
  rna_path = "../Transcriptome_diversity/data/RNA/work/Filtered_Noralised_Counts.csv",
  eqtl_path = "results/summary/eQTL_P_min.txt",
  fm_path = "../Transcriptome_diversity/results/SusieR2/PIP_Max.csv",
  output_path = "results/summary/Narrow_sense_h2_top_vs_fm.tsv",
  cores = default_cores,
  max_genes = Inf
)

parse_args <- function(args, config) {
  if (!length(args)) {
    return(config)
  }

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg == "--help") {
      cat(
        paste(
          "Usage:",
          "  Rscript scripts/local/analysis/3.Narrow_sense_h2.R [options]",
          "",
          "Options:",
          paste0("  --cores ", config$cores),
          paste0("  --max-genes ", config$max_genes),
          paste0("  --output-path ", config$output_path),
          paste0("  --eqtl-path ", config$eqtl_path),
          paste0("  --fm-path ", config$fm_path),
          sep = "\n"
        )
      )
      quit(save = "no", status = 0L)
    }

    if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=(.*)$", "\\1", arg)
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      if (i == length(args)) {
        stop("Missing value for argument: ", arg)
      }
      value <- args[[i + 1L]]
      i <- i + 1L
    } else {
      stop("Unrecognized argument: ", arg)
    }

    key <- gsub("-", "_", key)
    if (!key %in% names(config)) {
      stop("Unknown argument: --", gsub("_", "-", key))
    }

    if (key %in% c("cores", "max_genes")) {
      if (identical(value, "Inf")) {
        config[[key]] <- Inf
      } else {
        config[[key]] <- as.integer(value)
      }
    } else {
      config[[key]] <- value
    }

    i <- i + 1L
  }

  config
}

calc_h2 <- function(varcomp_df) {
  if (is.null(varcomp_df) || !nrow(varcomp_df)) {
    return(NA_real_)
  }

  line_idx <- grep("Line", rownames(varcomp_df), fixed = TRUE)
  if (!length(line_idx)) {
    return(NA_real_)
  }

  var_comp <- as.numeric(varcomp_df$VarComp)
  total_var <- sum(var_comp, na.rm = TRUE)

  if (!is.finite(total_var) || total_var <= 0) {
    return(NA_real_)
  }

  var_comp[[line_idx[[1L]]]] / total_var
}

calc_r2 <- function(y, X) {
  if (is.null(X)) {
    return(NA_real_)
  }

  X <- as.matrix(X)
  if (!ncol(X)) {
    return(NA_real_)
  }

  keep <- complete.cases(y, X)
  y <- y[keep]
  X <- X[keep, , drop = FALSE]

  if (length(y) < 3L || !nrow(X)) {
    return(NA_real_)
  }

  tss <- sum((y - mean(y))^2)
  if (!is.finite(tss) || tss <= 0) {
    return(NA_real_)
  }

  storage.mode(X) <- "double"
  fit <- lm.fit(x = cbind("(Intercept)" = 1, X), y = y)
  rss <- sum(fit$residuals^2)

  1 - (rss / tss)
}

lookup_named <- function(x, key, default = NA) {
  if (!key %in% names(x)) {
    return(default)
  }

  unname(x[[key]])
}

choose_first_column <- function(dt, candidates, label) {
  hit <- candidates[candidates %in% names(dt)]

  if (!length(hit)) {
    stop(
      "Could not detect a ", label, " column. Available columns: ",
      paste(names(dt), collapse = ", ")
    )
  }

  hit[[1L]]
}

config <- parse_args(commandArgs(trailingOnly = TRUE), config)
config$cores <- max(1L, as.integer(config$cores))

message("Loading genotype and covariate data.")

ind <- read.table(config$ind_path, stringsAsFactors = FALSE)
sample_ids <- ind$V1

kin <- attach.big.matrix(config$kin_desc_path)[]
geno <- attach.big.matrix(config$geno_desc_path)
map <- fread(config$map_path)
pc <- as.matrix(read.csv(config$pc_path, check.names = FALSE))
rna <- fread(config$rna_path)
eqtl <- fread(config$eqtl_path)
fm <- fread(config$fm_path)
rna_id_col <- choose_first_column(rna, c("Line", "Taxa"), "RNA sample ID")

if (nrow(pc) != length(sample_ids)) {
  stop("PC row count does not match the number of genotype samples.")
}

colnames(kin) <- sample_ids
rownames(kin) <- sample_ids
rownames(pc) <- sample_ids

# Numerical roundoff in the stored matrix can trigger symmetry warnings in sommer.
kin <- (kin + t(kin)) / 2

rna_idx <- match(sample_ids, rna[[rna_id_col]])
if (anyNA(rna_idx)) {
  stop("RNA matrix is missing some genotype samples.")
}
rna <- rna[rna_idx]

if (!identical(as.character(rna[[rna_id_col]]), as.character(sample_ids))) {
  stop("RNA samples could not be aligned to genotype sample order.")
}

setorder(eqtl, Gene, p)
eqtl_top <- eqtl[, .SD[1L], by = Gene]

snp_to_idx <- setNames(seq_len(nrow(map)), map$SNP)
top_snp_by_gene <- setNames(eqtl_top$SNP, eqtl_top$Gene)
top_p_by_gene <- setNames(eqtl_top$p, eqtl_top$Gene)
fm_snps_by_gene <- split(fm$SNP, fm$gene)

target_genes <- names(fm_snps_by_gene)
target_genes <- target_genes[target_genes %in% names(rna)]
fixed_formula <- as.formula(
  paste("expr ~", paste(colnames(pc), collapse = " + "))
)

if (is.finite(config$max_genes)) {
  target_genes <- head(target_genes, config$max_genes)
}

if (!length(target_genes)) {
  stop("No fine-mapped genes were found in both fm and RNA data.")
}

message("Genes to analyse: ", length(target_genes))
message("Parallel workers: ", config$cores)

analyze_gene <- function(gene_id) {
  result <- list(
    gene = gene_id,
    n_samples = NA_integer_,
    h2 = NA_real_,
    top_snp = NA_character_,
    top_p = NA_real_,
    top_r2 = NA_real_,
    n_fm_snps = NA_integer_,
    n_fm_snps_in_map = NA_integer_,
    fm_r2 = NA_real_,
    status = "ok",
    error_message = NA_character_
  )

  tryCatch({
    expr <- rna[[gene_id]]
    sample_idx <- which(is.finite(expr))

    if (length(sample_idx) < 3L) {
      stop("Too few non-missing expression values.")
    }

    result$n_samples <- length(sample_idx)

    pc_use <- pc[sample_idx, , drop = FALSE]
    pheno <- data.frame(
      Line = sample_ids[sample_idx],
      expr = expr[sample_idx],
      pc_use,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    kin_use <- kin[pheno$Line, pheno$Line, drop = FALSE]

    fit <- withCallingHandlers(
      sommer::mmer(
        fixed_formula,
        random = ~ vsr(Line, Gu = kin_use),
        rcov = ~ units,
        data = pheno,
        verbose = FALSE
      ),
      warning = function(w) {
        if (grepl("eig_sym\\(\\): given matrix is not symmetric", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )

    result$h2 <- calc_h2(summary(fit)$varcomp)
    residuals_vec <- as.numeric(fit$residuals)

    if (length(residuals_vec) != length(sample_idx)) {
      stop("Residual length does not match the phenotype length.")
    }

    result$top_snp <- lookup_named(top_snp_by_gene, gene_id, default = NA_character_)
    result$top_p <- lookup_named(top_p_by_gene, gene_id, default = NA_real_)

    if (!is.na(result$top_snp)) {
      top_idx <- unname(snp_to_idx[[result$top_snp]])
      if (!is.na(top_idx)) {
        top_geno <- geno[sample_idx, top_idx]
        result$top_r2 <- calc_r2(residuals_vec, top_geno)
      }
    }

    fm_snps <- unique(fm_snps_by_gene[[gene_id]])
    result$n_fm_snps <- length(fm_snps)

    fm_idx <- unname(snp_to_idx[fm_snps])
    fm_idx <- unique(fm_idx[!is.na(fm_idx)])
    result$n_fm_snps_in_map <- length(fm_idx)

    if (length(fm_idx)) {
      fm_geno <- geno[sample_idx, fm_idx, drop = FALSE]
      result$fm_r2 <- calc_r2(residuals_vec, fm_geno)
    }

    as.data.table(result)
  }, error = function(e) {
    result$status <- "error"
    result$error_message <- conditionMessage(e)
    as.data.table(result)
  })
}

results_list <- if (.Platform$OS.type == "unix" && config$cores > 1L) {
  parallel::mclapply(
    target_genes,
    analyze_gene,
    mc.cores = config$cores,
    mc.preschedule = FALSE
  )
} else {
  lapply(target_genes, analyze_gene)
}

results_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
setorder(results_dt, gene)

dir.create(dirname(config$output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(results_dt, config$output_path, sep = "\t", quote = FALSE, na = "NA")

message("Finished. Wrote ", nrow(results_dt), " rows to ", config$output_path)
