manPlotPrepareGG <- function(
    pvalues,
    chr,
    pos,
    snp = NULL,
    ismlog10 = FALSE,
    chrmargins = 5e6,
    thin_n = 300,
    y_bin_size = 0.01,
    seed = NULL
) {
  stopifnot(length(pvalues) == length(chr))
  stopifnot(length(pvalues) == length(pos))
  
  if (is.null(snp)) {
    snp <- names(pvalues)
  }
  if (is.null(snp)) {
    snp <- paste0("snp_", seq_along(pvalues))
  }
  stopifnot(length(pvalues) == length(snp))
  
  if (is.double(chr)) {
    chr <- as.integer(chr)
  }
  if (is.character(chr)) {
    chr <- factor(chr)
  }
  if (is.integer(chr)) {
    levels(chr) <- as.character(seq_len(tail(chr, 1)))
    class(chr) <- "factor"
  }
  
  pos <- as.numeric(pos)
  pvalues <- as.numeric(pvalues)
  
  poslist <- split(pos, chr, drop = FALSE)
  poslist[vapply(poslist, length, 0L) == 0L] <- list(0)
  chrmax <- vapply(poslist, max, 0) + 0
  names(chrmax) <- NULL
  
  offsets <- c(0, cumsum(chrmax)) + chrmargins
  names(offsets)[seq_along(poslist)] <- levels(chr)
  offsets <- offsets[seq_along(poslist)]
  
  x <- offsets[unclass(chr)] + pos
  if (ismlog10) {
    y <- pvalues
  } else {
    y <- -log10(pmax(pvalues, .Machine$double.xmin))
  }
  
  keep <- rep(TRUE, length(y))
  if (!is.null(thin_n) && is.finite(thin_n) && thin_n > 0L) {
    yfac <- as.integer(y / y_bin_size) + 1L
    ygroup <- split(seq_along(yfac), yfac)
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    keep <- rep(FALSE, length(y))
    for (idx in ygroup) {
      if (length(idx) > thin_n) {
        idx <- sample(idx, size = thin_n, replace = FALSE)
      }
      keep[idx] <- TRUE
    }
  }
  
  chr_index <- unclass(chr)
  chr_label <- gsub("chr", "", levels(chr))
  chr_table <- data.frame(
    chr = levels(chr),
    chr_label = chr_label,
    chr_index = seq_along(chr_label),
    chr_max = chrmax,
    offset = offsets,
    center = offsets + chrmax / 2,
    stringsAsFactors = FALSE
  )
  
  out <- data.frame(
    SNP = snp,
    CHR = as.character(chr),
    BP = pos,
    P = pvalues,
    x = x,
    y = y,
    chr_index = chr_index,
    keep = keep,
    stringsAsFactors = FALSE
  )
  
  attr(out, "chromosomes") <- chr_table
  attr(out, "chrmargins") <- chrmargins
  out
}

.resolve_manplot_colours <- function(colours, palette, chr_values) {
  if (!is.null(palette) && !is.null(colours)) {
    stop("Use either `palette` or `colours`, not both.")
  }
  
  palette_spec <- if (!is.null(palette)) palette else colours
  if (is.null(palette_spec)) {
    return(NULL)
  }
  
  chr_values <- unique(as.character(chr_values))
  n_chr <- length(chr_values)
  if (n_chr == 0L) {
    return(NULL)
  }
  
  resolved <- NULL
  if (is.function(palette_spec)) {
    resolved <- palette_spec(n_chr)
  } else if (is.character(palette_spec) && length(palette_spec) == 1L) {
    if (palette_spec %in% grDevices::hcl.pals()) {
      resolved <- grDevices::hcl.colors(n_chr, palette_spec)
    } else {
      stop(
        "Unknown palette name: ", palette_spec,
        ". Use a vector of colours, a palette function, or an HCL palette from grDevices::hcl.pals()."
      )
    }
  } else if (is.character(palette_spec) && length(palette_spec) >= 1L) {
    resolved <- if (length(palette_spec) >= n_chr) {
      palette_spec[seq_len(n_chr)]
    } else {
      grDevices::colorRampPalette(palette_spec)(n_chr)
    }
  } else {
    stop("`palette` must be NULL, a function, an HCL palette name, or a character vector of colours.")
  }
  
  stats::setNames(resolved, chr_values)
}

manPlotGG <- function(
    df,
    snp_col = "SNP",
    chr_col = "CHROM",
    pos_col = "POS",
    p_col = NULL,
    palette = NULL,
    colours = NULL,
    ismlog10 = FALSE,
    chrmargins = 5e6,
    thin_n = 300,
    y_bin_size = 0.01,
    seed = NULL,
    point_size = 0.7,
    point_alpha = 0.8
) {
  stopifnot(is.data.frame(df))
  stopifnot(snp_col %in% colnames(df))
  stopifnot(chr_col %in% colnames(df))
  stopifnot(pos_col %in% colnames(df))
  
  if (is.null(p_col)) {
    p_col <- tail(colnames(df), 1)
  }
  stopifnot(p_col %in% colnames(df))
  
  plot_df <- manPlotPrepareGG(
    pvalues = df[[p_col]],
    chr = df[[chr_col]],
    pos = df[[pos_col]],
    snp = df[[snp_col]],
    ismlog10 = ismlog10,
    chrmargins = chrmargins,
    thin_n = thin_n,
    y_bin_size = y_bin_size,
    seed = seed
  )
  
  chr_df <- attr(plot_df, "chromosomes")
  point_df <- plot_df[plot_df$keep, , drop = FALSE]
  point_colours <- .resolve_manplot_colours(
    colours = colours,
    palette = palette,
    chr_values = chr_df$chr
  )
  
  p <- ggplot2::ggplot(
    point_df,
    ggplot2::aes(x = x, y = y, colour = as.character(CHR))
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = expression(-log[10](P))) +
    ggplot2::scale_x_continuous(
      breaks = chr_df$center,
      labels = chr_df$chr_label
    )
  
  if (!is.null(point_colours)) {
    p <- p + ggplot2::scale_color_manual(values = point_colours)
  }
  
  list(
    data = plot_df,
    plot_data = point_df,
    chromosomes = chr_df,
    plot = p,
    p_col = p_col
  )
}

.susie_set_ld <- function(sets, ld) {
  set_names <- lapply(sets, names)
  
  if (!length(set_names)) {
    return(0)
  }
  
  if (length(set_names) == 1) {
    return(matrix(0, 1, 1, dimnames = list(names(set_names), names(set_names))))
  }
  
  set_ld <- matrix(
    0,
    nrow = length(set_names),
    ncol = length(set_names),
    dimnames = list(names(set_names), names(set_names))
  )
  
  for (i in 2:length(set_names)) {
    for (j in seq_len(i - 1)) {
      if (length(set_names[[i]]) && length(set_names[[j]])) {
        set_ld[i, j] <- max(ld[set_names[[i]], set_names[[j]]]^2, na.rm = TRUE)
      }
    }
  }
  
  set_ld
}

annotate_susie_fit <- function(res, snp = NULL, LD = NULL) {
  if (is.null(snp)) {
    snp <- colnames(res$alpha)
  }
  
  if (is.null(snp)) {
    snp <- paste0("snp_", seq_len(ncol(res$alpha)))
  }
  
  colnames(res$lbf_variable) <- c(snp, "null")[seq_len(ncol(res$lbf_variable))]
  colnames(res$alpha) <- c(snp, "null")[seq_len(ncol(res$alpha))]
  names(res$pip) <- snp
  
  if (length(res$sets$cs)) {
    res$sets$cs <- lapply(
      res$sets$cs,
      function(idx) {
        names(idx) <- snp[idx]
        idx
      }
    )
  }
  
  if (!is.null(LD)) {
    res$sld <- .susie_set_ld(res$sets$cs, LD)
  }
  
  res$pruned <- FALSE
  res
}

runsusie_raw <- function(
    X,
    y,
    suffix = 1,
    LD = NULL,
    maxit = 100,
    repeat_until_convergence = TRUE,
    s_init = NULL,
    retry_negative_residual = TRUE,
    ...
) {
  if (!(is.matrix(X) || inherits(X, "Matrix"))) {
    X <- as.matrix(X)
  }
  
  if (anyNA(X)) {
    stop("Input X must not contain missing values.")
  }
  
  if (anyNA(y)) {
    stop("Input y must not contain missing values.")
  }
  
  if (length(y) != nrow(X)) {
    stop("Input y must have the same length as the number of rows in X.")
  }
  
  snp <- colnames(X)
  if (is.null(snp)) {
    snp <- paste0("s", suffix, "_", seq_len(ncol(X)))
    colnames(X) <- snp
  }
  
  if (!is.null(LD)) {
    if (is.null(colnames(LD))) {
      colnames(LD) <- snp
    }
    if (is.null(rownames(LD))) {
      rownames(LD) <- snp
    }
    LD <- LD[snp, snp, drop = FALSE]
  }
  
  converged <- FALSE
  susie_args <- list(...)
  
  if ("max_iter" %in% names(susie_args)) {
    maxit <- susie_args$max_iter
    susie_args <- susie_args[setdiff(names(susie_args), "max_iter")]
  }
  
  if (!("s_init" %in% names(susie_args)) && !is.null(s_init)) {
    susie_args <- c(list(s_init = s_init), susie_args)
  }
  
  run_once <- function(max_iter, estimate_residual_variance) {
    do.call(
      susieR::susie,
      c(
        list(
          X = X,
          y = y,
          max_iter = max_iter,
          estimate_residual_variance = estimate_residual_variance
        ),
        susie_args
      )
    )
  }
  
  estimate_residual_variance <- TRUE
  if ("estimate_residual_variance" %in% names(susie_args)) {
    estimate_residual_variance <- isTRUE(susie_args$estimate_residual_variance)
    susie_args <- susie_args[setdiff(names(susie_args), "estimate_residual_variance")]
  }
  
  while (!converged) {
    message("running max iterations: ", maxit)
    
    res <- tryCatch(
      run_once(max_iter = maxit, estimate_residual_variance = estimate_residual_variance),
      error = function(e) {
        if (
          retry_negative_residual &&
          estimate_residual_variance &&
          grepl("negative", conditionMessage(e), ignore.case = TRUE)
        ) {
          message("retrying with estimate_residual_variance = FALSE")
          return(run_once(max_iter = maxit, estimate_residual_variance = FALSE))
        }
        
        stop(e)
      }
    )
    
    converged <- isTRUE(res$converged)
    message("\tconverged: ", converged)
    
    if (!converged && !repeat_until_convergence) {
      stop(
        "susie() did not converge in ",
        maxit,
        " iterations. Try running with repeat_until_convergence=TRUE"
      )
    }
    
    if (!converged) {
      maxit <- maxit * 100
    }
  }
  
  annotate_susie_fit(res = res, snp = snp, LD = LD)
}

run_coloc_susie_raw <- function(susie_fit1, susie_fit2, ...) {
  if (is.null(susie_fit1) || is.null(susie_fit2)) {
    stop("Both `susie_fit1` and `susie_fit2` are required.")
  }
  
  validate_fit <- function(fit, label) {
    required_slots <- c("alpha", "pip", "lbf_variable", "sets")
    missing_slots <- required_slots[!required_slots %in% names(fit)]
    if (length(missing_slots)) {
      stop(label, " is missing required SuSiE fields: ", paste(missing_slots, collapse = ", "))
    }
    if (!is.matrix(fit$lbf_variable)) {
      stop(label, "$lbf_variable must be a matrix.")
    }
  }
  
  validate_fit(susie_fit1, "susie_fit1")
  validate_fit(susie_fit2, "susie_fit2")
  
  # coloc::coloc.susie() assumes that coloc.bf_bf() always returns a
  # data.table. With some coloc/data.table versions, ret$summary is returned
  # as a data.frame and the package's subsequent := assignments fail. Keep
  # the same coloc logic here, while restoring the expected class explicitly.
  if ("susie" %in% class(susie_fit1)) {
    s1 <- susie_fit1
  } else {
    s1 <- runsusie(d = susie_fit1, suffix = 1)
  }
  if ("susie" %in% class(susie_fit2)) {
    s2 <- susie_fit2
  } else {
    s2 <- runsusie(d = susie_fit2, suffix = 2)
  }

  cs1 <- s1$sets
  cs2 <- s2$sets
  if (is.null(cs1$cs) || is.null(cs2$cs) ||
      length(cs1$cs) == 0 || length(cs2$cs) == 0 ||
      is.null(cs1$cs_index) || is.null(cs2$cs_index) ||
      length(cs1$cs_index) == 0 || length(cs2$cs_index) == 0) {
    return(list(
      summary = data.table::data.table(
        status = "no_credible_sets",
        nsnps = NA_integer_,
        idx1 = NA_integer_,
        idx2 = NA_integer_
      ),
      results = data.table::data.table(),
      priors = NULL
    ))
  }

  idx1 <- cs1$cs_index
  idx2 <- cs2$cs_index
  bf1 <- s1$lbf_variable[idx1, , drop = FALSE]
  bf2 <- s2$lbf_variable[idx2, , drop = FALSE]

  snp1 <- colnames(bf1)
  snp2 <- colnames(bf2)
  if (is.null(snp1) && is.null(snp2)) {
    if (ncol(bf1) != ncol(bf2)) {
      stop(
        "The SuSiE fits have no SNP names and contain different numbers of variants (",
        ncol(bf1), " versus ", ncol(bf2), ")."
      )
    }
    snp1 <- snp2 <- paste0("snp_", seq_len(ncol(bf1)))
    colnames(bf1) <- snp1
    colnames(bf2) <- snp2
    warning("The SuSiE fits have no SNP names; assuming identical variant order.")
  } else if (is.null(snp1) || is.null(snp2)) {
    stop("Both SuSiE fits must have SNP names, or neither can have SNP names.")
  } else {
    if (anyDuplicated(snp1) || anyDuplicated(snp2)) {
      stop("SNP names must be unique in both SuSiE fits.")
    }
    common_snp <- intersect(setdiff(snp1, "null"), setdiff(snp2, "null"))
    if (!length(common_snp)) {
      stop("The two SuSiE fits have no SNP names in common.")
    }
    keep1 <- c(common_snp, intersect("null", snp1))
    keep2 <- c(common_snp, intersect("null", snp2))
    bf1 <- bf1[, keep1, drop = FALSE]
    bf2 <- bf2[, keep2, drop = FALSE]
  }

  ret <- coloc::coloc.bf_bf(bf1, bf2, ...)
  if (!is.list(ret) || is.null(ret$summary)) {
    stop("coloc.bf_bf() did not return a coloc result; check SNP-name overlap.")
  }

  # coloc.bf_bf() may return a data.frame here, depending on package
  # versions. Rebuild the table before adding columns rather than using :=;
  # some data.table versions report an internal allocation error when the
  # returned table has been shallow-copied by coloc.
  summary_df <- as.data.frame(ret$summary)
  if (all(c("idx1", "idx2") %in% colnames(summary_df))) {
    summary_df[["idx1"]] <- idx1[summary_df[["idx1"]]]
    summary_df[["idx2"]] <- idx2[summary_df[["idx2"]]]
  }
  list(
    summary = data.table::as.data.table(summary_df),
    results = ret$results,
    priors = ret$priors
  )
}


###
remove.covariate.effects <- function (X, Z, y) {
  # include the intercept term
  if (any(Z[,1]!=1)) Z = cbind(1, Z)
  A   <- Matrix::forceSymmetric(crossprod(Z))
  SZy <- as.vector(solve(A,c(y %*% Z)))
  SZX <- as.matrix(solve(A,t(Z) %*% X))
  y <- y - c(Z %*% SZy)
  X <- X - Z %*% SZX
  return(list(X = X,y = y,SZy = SZy,SZX = SZX))
}

.fmc_require_cols <- function(x, cols, label) {
  missing_cols <- setdiff(cols, colnames(x))
  if (length(missing_cols)) {
    stop(label, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

.fmc_chr_to_int <- function(x) {
  x_chr <- as.character(x)
  x_chr <- sub("^chr", "", x_chr, ignore.case = TRUE)
  x_chr[!grepl("^[0-9]+$", x_chr)] <- NA_character_
  as.integer(x_chr)
}

.fmc_capture_fit <- function(expr_fun) {
  tryCatch(
    list(fit = expr_fun(), error = NA_character_),
    error = function(e) list(fit = NULL, error = conditionMessage(e))
  )
}

.fmc_sample_cis_snps <- function(snp, n = 3L, seed = NULL) {
  out <- rep(NA_character_, n)
  snp <- unique(as.character(snp))
  snp <- snp[!is.na(snp) & nzchar(snp)]
  
  if (!length(snp)) {
    return(out)
  }
  
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  picked <- sample(snp, size = min(n, length(snp)), replace = FALSE)
  out[seq_along(picked)] <- picked
  out
}

.fmc_prepare_condition_input <- function(condition, expr_tbl, expr_pc_tbl, taxa_tbl, geno_pc_tbl) {
  .fmc_require_cols(expr_tbl, "Taxa", paste0(condition, " expression table"))
  .fmc_require_cols(expr_pc_tbl, "Taxa", paste0(condition, " expression PC table"))
  .fmc_require_cols(taxa_tbl, "Taxa", "taxa table")
  .fmc_require_cols(geno_pc_tbl, "Taxa", "genotype PC table")
  
  taxa <- taxa_tbl$Taxa
  
  expr_idx <- match(taxa, expr_tbl$Taxa)
  expr_mat <- as.matrix(expr_tbl[expr_idx, setdiff(colnames(expr_tbl), "Taxa"), drop = FALSE])
  mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- taxa
  
  geno_pc_idx <- match(taxa, geno_pc_tbl$Taxa)
  expr_pc_idx <- match(taxa, expr_pc_tbl$Taxa)
  
  cov_mat <- cbind(
    as.matrix(geno_pc_tbl[geno_pc_idx, setdiff(colnames(geno_pc_tbl), "Taxa"), drop = FALSE]),
    as.matrix(expr_pc_tbl[expr_pc_idx, setdiff(colnames(expr_pc_tbl), "Taxa"), drop = FALSE])
  )
  mode(cov_mat) <- "numeric"
  rownames(cov_mat) <- taxa
  
  list(
    condition = condition,
    taxa = taxa,
    expr_mat = expr_mat,
    cov_mat = cov_mat,
    genes = colnames(expr_mat)
  )
}

.fmc_format_mvp_mlm <- function(gwas_fit, gene_id, condition) {
  gwas_dt <- data.table::as.data.table(data.frame(gwas_fit$map, gwas_fit$mlm.results, check.names = FALSE))
  p_col <- paste0(gene_id, ".MLM")
  
  if (!p_col %in% colnames(gwas_dt)) {
    stop("Expected MLM p-value column not found: ", p_col)
  }
  .fmc_require_cols(gwas_dt, c("SNP", "CHROM", "POS", "A1", "A2", "MAF", "Effect", "SE"), "MVP MLM results")
  
  gwas_dt[, `:=`(
    gene_id = gene_id,
    condition = condition,
    CHROM = .fmc_chr_to_int(CHROM),
    POS = as.integer(POS),
    MAF = as.numeric(MAF),
    Effect = as.numeric(Effect),
    SE = as.numeric(SE),
    p_value = as.numeric(get(p_col))
  )]
  
  gwas_dt[, .(
    gene_id,
    condition,
    SNP,
    CHROM,
    POS,
    A1,
    A2,
    MAF,
    Effect,
    SE,
    p_value
  )]
}

.fmc_extract_cs_membership <- function(fit, snp) {
  if (is.null(fit) || is.null(fit$sets$cs) || !length(fit$sets$cs)) {
    return(data.table::data.table(SNP = character(), cs_id = character()))
  }
  
  cs_list <- lapply(seq_along(fit$sets$cs), function(i) {
    idx <- fit$sets$cs[[i]]
    if (!length(idx)) {
      return(NULL)
    }
    data.table::data.table(
      SNP = snp[idx],
      cs_id = paste0("CS", i)
    )
  })
  
  data.table::rbindlist(cs_list, use.names = TRUE, fill = TRUE)
}

.fmc_extract_susie_tables <- function(gene_id, condition, fit, region_dt, n_samples, error_message) {
  if (is.null(fit)) {
    return(list(
      summary = data.table::data.table(
        gene_id = gene_id,
        condition = condition,
        status = "failed",
        error = error_message,
        n_samples = n_samples,
        n_snps = nrow(region_dt),
        n_credible_sets = NA_integer_,
        top_snp = NA_character_,
        top_pip = NA_real_
      ),
      pip = data.table::data.table(
        gene_id = character(),
        condition = character(),
        SNP = character(),
        CHROM = integer(),
        POS = integer(),
        pip = numeric(),
        cs_id = character(),
        in_credible_set = logical(),
        error = character()
      )
    ))
  }
  
  pip_dt <- data.table::copy(region_dt[, .(SNP, CHROM, POS)])
  pip_dt[, `:=`(
    gene_id = gene_id,
    condition = condition,
    pip = as.numeric(fit$pip),
    error = NA_character_
  )]
  
  cs_dt <- .fmc_extract_cs_membership(fit, pip_dt$SNP)
  if (nrow(cs_dt)) {
    pip_dt <- cs_dt[pip_dt, on = "SNP"]
  } else {
    pip_dt[, cs_id := NA_character_]
  }
  pip_dt[, in_credible_set := !is.na(cs_id)]
  pip_dt <- pip_dt[in_credible_set == TRUE & is.finite(pip)]
  data.table::setorder(pip_dt, cs_id, -pip, POS)
  data.table::setcolorder(
    pip_dt,
    c("gene_id", "condition", "SNP", "CHROM", "POS", "pip", "cs_id", "in_credible_set", "error")
  )
  
  top_idx <- which.max(fit$pip)
  summary_dt <- data.table::data.table(
    gene_id = gene_id,
    condition = condition,
    status = "ok",
    error = NA_character_,
    n_samples = n_samples,
    n_snps = length(fit$pip),
    n_credible_sets = if (is.null(fit$sets$cs)) 0L else length(fit$sets$cs),
    top_snp = pip_dt$SNP[top_idx],
    top_pip = max(fit$pip, na.rm = TRUE)
  )
  
  list(summary = summary_dt, pip = pip_dt)
}

.fmc_extract_coloc_tables <- function(gene_id, condition1, condition2, coloc_result) {
  pair_id <- paste(condition1, condition2, sep = "_vs_")
  
  if (is.null(coloc_result$fit)) {
    summary_dt <- data.table::data.table(
      gene_id = gene_id,
      pair_id = pair_id,
      condition1 = condition1,
      condition2 = condition2,
      status = "failed",
      error = coloc_result$error
    )
    return(list(summary = summary_dt, results = data.table::data.table()))
  }
  
  summary_dt <- data.table::as.data.table(coloc_result$fit$summary)
  summary_dt[, `:=`(
    gene_id = gene_id,
    pair_id = pair_id,
    condition1 = condition1,
    condition2 = condition2,
    status = "ok",
    error = NA_character_
  )]
  data.table::setcolorder(summary_dt, c("gene_id", "pair_id", "condition1", "condition2", "status", "error"))
  
  results_dt <- data.table::as.data.table(coloc_result$fit$results)
  if (nrow(results_dt)) {
    results_dt[, `:=`(
      gene_id = gene_id,
      pair_id = pair_id,
      condition1 = condition1,
      condition2 = condition2
    )]
    data.table::setcolorder(results_dt, c("gene_id", "pair_id", "condition1", "condition2"))
  }
  
  list(summary = summary_dt, results = results_dt)
}

.run_finemap_coloc_single_gene <- function(
    gene_id,
    gene_seed,
    tss_row,
    condition_inputs,
    geno,
    kin,
    map_dt,
    map_mvp,
    cis_window,
    gwas_threshold,
    gwas_ncpus,
    vc.method,
    susie_L,
    susie_max_iter,
    min_cis_snps,
    random_cis_n,
    verbose
) {
  chrom <- as.integer(tss_row$chromosome[[1]])
  tss <- as.integer(tss_row$TSS[[1]])
  cis_start <- max(0L, as.integer(tss - cis_window))
  cis_end <- as.integer(tss + cis_window)
  
  region_idx <- which(map_dt$CHROM == chrom & map_dt$POS >= cis_start & map_dt$POS <= cis_end)
  region_dt <- data.table::copy(map_dt[region_idx, .(SNP, CHROM, POS)])
  random_cis_snp <- .fmc_sample_cis_snps(region_dt$SNP, n = random_cis_n, seed = gene_seed)
  random_cis_snps <- data.table::data.table(
    gene_id = gene_id,
    chromosome = chrom,
    TSS = tss,
    cis_random_rank = seq_along(random_cis_snp),
    SNP = random_cis_snp
  )
  
  condition_summaries <- list()
  gwas_by_condition <- list()
  any_cis_signal <- FALSE
  
  for (condition in names(condition_inputs)) {
    condition_input <- condition_inputs[[condition]]
    
    if (!gene_id %in% condition_input$genes) {
      condition_summaries[[condition]] <- data.table::data.table(
        gene_id = gene_id,
        condition = condition,
        status = "gene_not_found",
        n_samples = 0L,
        n_cis_snps = nrow(region_dt),
        n_cis_significant = 0L,
        min_cis_p = NA_real_,
        lead_cis_snp = NA_character_,
        lead_cis_p = NA_real_,
        candidate_condition = FALSE
      )
      next
    }
    
    phenotype <- condition_input$expr_mat[, gene_id]
    phe_df <- data.frame(Taxa = condition_input$taxa, trait = phenotype, check.names = FALSE)
    colnames(phe_df)[2] <- gene_id
    
    gwas_fit <- rMVP::MVP(
      phe = phe_df,
      geno = geno,
      map = map_mvp,
      K = kin,
      CV.MLM = condition_input$cov_mat,
      ncpus = gwas_ncpus,
      vc.method = vc.method,
      method = "MLM",
      file.output = FALSE,
      verbose = verbose
    )
    
    gwas_dt <- .fmc_format_mvp_mlm(gwas_fit, gene_id = gene_id, condition = condition)
    gwas_dt[, is_cis := CHROM == chrom & POS >= cis_start & POS <= cis_end]
    
    cis_dt <- gwas_dt[is_cis == TRUE]
    n_cis_sig <- if (nrow(cis_dt)) sum(cis_dt$p_value <= gwas_threshold, na.rm = TRUE) else 0L
    lead_idx <- if (nrow(cis_dt)) which.min(cis_dt$p_value) else integer()
    
    condition_summaries[[condition]] <- data.table::data.table(
      gene_id = gene_id,
      condition = condition,
      status = "ok",
      n_samples = sum(stats::complete.cases(cbind(phenotype, condition_input$cov_mat))),
      n_cis_snps = nrow(cis_dt),
      n_cis_significant = n_cis_sig,
      min_cis_p = if (nrow(cis_dt)) min(cis_dt$p_value, na.rm = TRUE) else NA_real_,
      lead_cis_snp = if (length(lead_idx)) cis_dt$SNP[lead_idx] else NA_character_,
      lead_cis_p = if (length(lead_idx)) cis_dt$p_value[lead_idx] else NA_real_,
      candidate_condition = n_cis_sig > 0L
    )
    
    any_cis_signal <- any_cis_signal || (n_cis_sig > 0L)
    gwas_by_condition[[condition]] <- gwas_dt
  }
  
  gene_summary <- data.table::rbindlist(condition_summaries, use.names = TRUE, fill = TRUE)
  gene_summary <- gene_summary[, .(
    gene_id = gene_id[[1]],
    chromosome = chrom,
    TSS = tss,
    threshold = gwas_threshold,
    candidate_gene = any(candidate_condition),
    n_conditions_tested = sum(status == "ok"),
    n_conditions_with_cis_signal = sum(candidate_condition, na.rm = TRUE),
    min_cis_p_any = suppressWarnings(min(min_cis_p, na.rm = TRUE))
  )]
  random_snp_cols <- as.list(random_cis_snp)
  names(random_snp_cols) <- paste0("random_cis_snp_", seq_along(random_cis_snp))
  gene_summary[, (names(random_snp_cols)) := random_snp_cols]
  if (!is.finite(gene_summary$min_cis_p_any)) {
    gene_summary[, min_cis_p_any := NA_real_]
  }
  
  random_cis_snps <- data.table::rbindlist(
    lapply(names(condition_summaries), function(condition) {
      gwas_dt <- gwas_by_condition[[condition]]
      if (is.null(gwas_dt) || !nrow(random_cis_snps)) {
        out <- data.table::copy(random_cis_snps)
        out[, `:=`(
          condition = condition,
          A1 = NA_character_,
          A2 = NA_character_,
          MAF = NA_real_,
          Effect = NA_real_,
          SE = NA_real_,
          p_value = NA_real_
        )]
        return(out)
      }
      
      merge(
        data.table::copy(random_cis_snps),
        gwas_dt[, .(condition, SNP, A1, A2, MAF, Effect, SE, p_value)],
        by = "SNP",
        all.x = TRUE,
        sort = FALSE
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setcolorder(
    random_cis_snps,
    c("gene_id", "condition", "chromosome", "TSS", "cis_random_rank", "SNP", "A1", "A2", "MAF", "Effect", "SE", "p_value")
  )
  
  if (!any_cis_signal) {
    return(list(
      gene_summary = gene_summary,
      condition_summary = data.table::rbindlist(condition_summaries, use.names = TRUE, fill = TRUE),
      random_cis_snps = random_cis_snps,
      gwas_cis = data.table::data.table(),
      gwas_trans = data.table::data.table(),
      finemap_summary = data.table::data.table(),
      finemap_pip = data.table::data.table(),
      coloc_summary = data.table::data.table(),
      coloc_results = data.table::data.table()
    ))
  }
  
  cis_results <- data.table::rbindlist(
    lapply(gwas_by_condition, function(x) x[is_cis == TRUE]),
    use.names = TRUE,
    fill = TRUE
  )
  
  trans_results <- data.table::rbindlist(
    lapply(gwas_by_condition, function(x) x[is_cis == FALSE & p_value <= gwas_threshold]),
    use.names = TRUE,
    fill = TRUE
  )
  
  local_geno <- geno[, region_idx]
  if (length(region_idx) == 1L) {
    local_geno <- matrix(local_geno, ncol = 1L)
  }
  local_geno <- as.matrix(local_geno)
  rownames(local_geno) <- condition_inputs[[1]]$taxa
  colnames(local_geno) <- region_dt$SNP
  
  set.seed(gene_seed)
  
  finemap_summary <- list()
  finemap_pip <- list()
  susie_fits <- list()
  
  for (condition in names(condition_inputs)) {
    condition_input <- condition_inputs[[condition]]
    
    if (!gene_id %in% condition_input$genes || nrow(region_dt) < min_cis_snps) {
      fit_tables <- .fmc_extract_susie_tables(
        gene_id = gene_id,
        condition = condition,
        fit = NULL,
        region_dt = region_dt,
        n_samples = 0L,
        error_message = if (nrow(region_dt) < min_cis_snps) {
          paste0("Too few cis SNPs (", nrow(region_dt), ")")
        } else {
          "Gene not available in this condition"
        }
      )
      finemap_summary[[condition]] <- fit_tables$summary
      finemap_pip[[condition]] <- fit_tables$pip
      next
    }
    
    phenotype <- condition_input$expr_mat[, gene_id]
    keep_idx <- which(is.finite(phenotype) & stats::complete.cases(condition_input$cov_mat))
    region_geno <- local_geno[keep_idx, , drop = FALSE]
    cov_mat <- condition_input$cov_mat[keep_idx, , drop = FALSE]
    y <- phenotype[keep_idx]
    
    adjusted <- remove.covariate.effects(
      X = region_geno,
      Z = cov_mat,
      y = y
    )
    
    snp_var <- apply(adjusted$X, 2, stats::var)
    keep_snp <- is.finite(snp_var) & snp_var > 0
    
    if (sum(keep_snp) < min_cis_snps) {
      fit_tables <- .fmc_extract_susie_tables(
        gene_id = gene_id,
        condition = condition,
        fit = NULL,
        region_dt = region_dt[keep_snp],
        n_samples = length(keep_idx),
        error_message = paste0("Too few variable cis SNPs after adjustment (", sum(keep_snp), ")")
      )
      finemap_summary[[condition]] <- fit_tables$summary
      finemap_pip[[condition]] <- fit_tables$pip
      next
    }
    
    adjusted_x <- adjusted$X[, keep_snp, drop = FALSE]
    ld_mat <- stats::cor(adjusted_x)
    colnames(ld_mat) <- colnames(adjusted_x)
    rownames(ld_mat) <- colnames(adjusted_x)
    region_kept <- region_dt[keep_snp]
    
    susie_res <- .fmc_capture_fit(
      function() {
        runsusie_raw(
          X = adjusted_x,
          y = adjusted$y,
          LD = ld_mat,
          L = susie_L,
          standardize = TRUE,
          max_iter = susie_max_iter,
          verbose = FALSE
        )
      }
    )
    
    susie_fits[[condition]] <- susie_res$fit
    
    fit_tables <- .fmc_extract_susie_tables(
      gene_id = gene_id,
      condition = condition,
      fit = susie_res$fit,
      region_dt = region_kept,
      n_samples = length(keep_idx),
      error_message = susie_res$error
    )
    finemap_summary[[condition]] <- fit_tables$summary
    finemap_pip[[condition]] <- fit_tables$pip
  }
  
  coloc_pairs <- if (length(susie_fits) >= 2L) {
    utils::combn(names(susie_fits), 2, simplify = FALSE)
  } else {
    list()
  }
  coloc_summary <- list()
  coloc_results <- list()
  
  for (i in seq_along(coloc_pairs)) {
    pair <- coloc_pairs[[i]]
    condition1 <- pair[[1]]
    condition2 <- pair[[2]]
    
    if (is.null(susie_fits[[condition1]]) || is.null(susie_fits[[condition2]])) {
      coloc_res <- list(
        fit = NULL,
        error = paste("Missing fine-mapping fit for", paste(pair, collapse = " vs "))
      )
    } else {
      coloc_res <- .fmc_capture_fit(
        function() run_coloc_susie_raw(susie_fits[[condition1]], susie_fits[[condition2]])
      )
    }
    
    coloc_tables <- .fmc_extract_coloc_tables(
      gene_id = gene_id,
      condition1 = condition1,
      condition2 = condition2,
      coloc_result = coloc_res
    )
    coloc_summary[[i]] <- coloc_tables$summary
    coloc_results[[i]] <- coloc_tables$results
  }
  
  list(
    gene_summary = gene_summary,
    condition_summary = data.table::rbindlist(condition_summaries, use.names = TRUE, fill = TRUE),
    random_cis_snps = random_cis_snps,
    gwas_cis = cis_results,
    gwas_trans = trans_results,
    finemap_summary = data.table::rbindlist(finemap_summary, use.names = TRUE, fill = TRUE),
    finemap_pip = data.table::rbindlist(finemap_pip, use.names = TRUE, fill = TRUE),
    coloc_summary = if (length(coloc_summary)) data.table::rbindlist(coloc_summary, use.names = TRUE, fill = TRUE) else data.table::data.table(),
    coloc_results = if (length(coloc_results)) data.table::rbindlist(coloc_results, use.names = TRUE, fill = TRUE) else data.table::data.table()
  )
}

load_finemap_coloc_inputs <- function(
    taxa_path = "data/processed/rMVP/maize_262.geno.ind",
    geno_desc_path = "data/processed/rMVP/maize_262.geno.desc",
    kin_desc_path = "data/processed/rMVP/maize_262.kin.desc",
    map_path = "data/processed/rMVP/maize_262.geno.map",
    geno_pc_desc_path = "data/processed/rMVP/maize_262.pc.desc",
    control_expr_path = "data/processed/limma_normalized_counts/RNAseq_control_limma_voom_normalized_logCPM.csv",
    cold_expr_path = "data/processed/limma_normalized_counts/RNAseq_cold_limma_voom_normalized_logCPM.csv",
    chill_expr_path = "data/processed/limma_normalized_counts/RNAseq_chilling_limma_voom_normalized_logCPM.csv",
    control_pc_path = "results/full/PCAeQTL/PCA_control_scores_elbow.csv",
    cold_pc_path = "results/full/PCAeQTL/PCA_cold_scores_elbow.csv",
    chill_pc_path = "results/full/PCAeQTL/PCA_chilling_scores_elbow.csv",
    tss_path = "data/processed/meta/maize_v5_gene_tss.tsv"
) {
  taxa_tbl <- data.table::fread(taxa_path, header = FALSE, col.names = "Taxa", data.table = FALSE)
  geno <- bigmemory::attach.big.matrix(geno_desc_path)
  kin <- bigmemory::attach.big.matrix(kin_desc_path)
  map <- data.table::fread(map_path, data.table = FALSE)
  
  geno_pc_tbl <- as.data.frame(bigmemory::attach.big.matrix(geno_pc_desc_path)[])
  colnames(geno_pc_tbl) <- paste0("geno_pc_", seq_len(ncol(geno_pc_tbl)))
  geno_pc_tbl <- data.frame(Taxa = taxa_tbl$Taxa, geno_pc_tbl, check.names = FALSE)
  
  expr_tbls <- list(
    control = data.table::fread(control_expr_path, data.table = FALSE),
    cold = data.table::fread(cold_expr_path, data.table = FALSE),
    chill = data.table::fread(chill_expr_path, data.table = FALSE)
  )
  
  expr_pc_tbls <- list(
    control = data.table::fread(control_pc_path, data.table = FALSE),
    cold = data.table::fread(cold_pc_path, data.table = FALSE),
    chill = data.table::fread(chill_pc_path, data.table = FALSE)
  )
  
  tss_tbl <- data.table::fread(tss_path, data.table = FALSE)
  
  list(
    taxa_tbl = taxa_tbl,
    geno = geno,
    kin = kin,
    map = map,
    geno_pc_tbl = geno_pc_tbl,
    expr_tbls = expr_tbls,
    expr_pc_tbls = expr_pc_tbls,
    tss_tbl = tss_tbl
  )
}

run_finemap_coloc_batch <- function(
    taxa_tbl,
    geno,
    kin,
    map,
    geno_pc_tbl,
    expr_tbls,
    expr_pc_tbls,
    tss_tbl,
    genes = NULL,
    conditions = c("control", "cold", "chill"),
    cis_window = 1e6,
    gwas_threshold = NULL,
    gwas_ncpus = 16L,
    vc.method = "EMMA",
    susie_L = 10L,
    susie_max_iter = 500L,
    susie_seed = 1L,
    min_cis_snps = 10L,
    random_cis_n = 3L,
    out_dir = NULL,
    n_workers = 1L,
    verbose = TRUE
) {
  .fmc_require_cols(taxa_tbl, "Taxa", "taxa table")
  .fmc_require_cols(geno_pc_tbl, "Taxa", "genotype PC table")
  .fmc_require_cols(tss_tbl, c("gene_id", "chromosome", "TSS"), "TSS table")
  
  if (!all(conditions %in% names(expr_tbls))) {
    stop("Missing expression tables for conditions: ", paste(setdiff(conditions, names(expr_tbls)), collapse = ", "))
  }
  if (!all(conditions %in% names(expr_pc_tbls))) {
    stop("Missing expression PC tables for conditions: ", paste(setdiff(conditions, names(expr_pc_tbls)), collapse = ", "))
  }
  
  map_dt <- data.table::as.data.table(data.table::copy(map))
  .fmc_require_cols(map_dt, c("SNP", "CHROM", "POS", "A1", "A2"), "map table")
  map_dt[, `:=`(
    CHROM = .fmc_chr_to_int(CHROM),
    POS = as.integer(POS)
  )]
  map_mvp <- as.data.frame(map_dt[, .(SNP, CHROM, POS, A1, A2)])
  
  if (is.null(gwas_threshold)) {
    gwas_threshold <- 0.05 / nrow(map_dt)
  }
  
  tss_dt <- data.table::as.data.table(data.table::copy(tss_tbl))
  tss_dt[, chromosome := .fmc_chr_to_int(chromosome)]
  tss_dt[, TSS := as.integer(TSS)]
  tss_dt <- tss_dt[!is.na(gene_id) & !is.na(chromosome) & !is.na(TSS)]
  tss_dt <- unique(tss_dt, by = "gene_id")
  
  condition_inputs <- lapply(
    conditions,
    function(condition) {
      .fmc_prepare_condition_input(
        condition = condition,
        expr_tbl = expr_tbls[[condition]],
        expr_pc_tbl = expr_pc_tbls[[condition]],
        taxa_tbl = taxa_tbl,
        geno_pc_tbl = geno_pc_tbl
      )
    }
  )
  names(condition_inputs) <- conditions
  
  if (is.null(genes)) {
    available_genes <- unique(unlist(lapply(condition_inputs, `[[`, "genes"), use.names = FALSE))
    genes <- intersect(tss_dt$gene_id, available_genes)
  } else {
    genes <- intersect(as.character(genes), tss_dt$gene_id)
  }
  
  if (!length(genes)) {
    stop("No genes remain after intersecting the requested genes, TSS table, and expression matrices.")
  }
  
  if (verbose) {
    message("Testing ", length(genes), " genes with GWAS threshold ", signif(gwas_threshold, 4))
    if (n_workers > 1L && gwas_ncpus > 1L) {
      message("Requested parallel genes and threaded MVP calls; tune n_workers and gwas_ncpus together to avoid CPU oversubscription.")
    }
  }
  
  gene_rows <- split(tss_dt[match(genes, tss_dt$gene_id)], seq_len(length(genes)))
  
  worker_fun <- function(i) {
    .run_finemap_coloc_single_gene(
      gene_id = genes[[i]],
      gene_seed = susie_seed + i,
      tss_row = gene_rows[[i]],
      condition_inputs = condition_inputs,
      geno = geno,
      kin = kin,
      map_dt = map_dt,
      map_mvp = map_mvp,
      cis_window = cis_window,
      gwas_threshold = gwas_threshold,
      gwas_ncpus = gwas_ncpus,
      vc.method = vc.method,
      susie_L = susie_L,
      susie_max_iter = susie_max_iter,
      min_cis_snps = min_cis_snps,
      random_cis_n = random_cis_n,
      verbose = verbose
    )
  }
  
  if (.Platform$OS.type == "unix" && isTRUE(n_workers > 1L)) {
    results_list <- parallel::mclapply(
      seq_along(genes),
      worker_fun,
      mc.cores = n_workers,
      mc.preschedule = FALSE,
      mc.set.seed = TRUE
    )
  } else {
    results_list <- lapply(seq_along(genes), worker_fun)
  }
  
  combined <- list(
    gene_summary = data.table::rbindlist(lapply(results_list, `[[`, "gene_summary"), use.names = TRUE, fill = TRUE),
    condition_summary = data.table::rbindlist(lapply(results_list, `[[`, "condition_summary"), use.names = TRUE, fill = TRUE),
    random_cis_snps = data.table::rbindlist(lapply(results_list, `[[`, "random_cis_snps"), use.names = TRUE, fill = TRUE),
    gwas_cis = data.table::rbindlist(lapply(results_list, `[[`, "gwas_cis"), use.names = TRUE, fill = TRUE),
    gwas_trans = data.table::rbindlist(lapply(results_list, `[[`, "gwas_trans"), use.names = TRUE, fill = TRUE),
    finemap_summary = data.table::rbindlist(lapply(results_list, `[[`, "finemap_summary"), use.names = TRUE, fill = TRUE),
    finemap_pip = data.table::rbindlist(lapply(results_list, `[[`, "finemap_pip"), use.names = TRUE, fill = TRUE),
    coloc_summary = data.table::rbindlist(lapply(results_list, `[[`, "coloc_summary"), use.names = TRUE, fill = TRUE),
    coloc_results = data.table::rbindlist(lapply(results_list, `[[`, "coloc_results"), use.names = TRUE, fill = TRUE)
  )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    write_if_present <- function(dt, path) {
      if (!is.null(dt) && ncol(dt) > 0L) {
        data.table::fwrite(dt, path, sep = "\t")
      }
    }
    
    write_if_present(combined$gene_summary, file.path(out_dir, "gene_summary.tsv"))
    write_if_present(combined$condition_summary, file.path(out_dir, "condition_summary.tsv"))
    write_if_present(combined$random_cis_snps, file.path(out_dir, "random_cis_snps.tsv"))
    write_if_present(combined$gwas_cis, file.path(out_dir, "gwas_cis.tsv.gz"))
    write_if_present(combined$gwas_trans, file.path(out_dir, "gwas_trans_significant.tsv.gz"))
    write_if_present(combined$finemap_summary, file.path(out_dir, "finemap_summary.tsv"))
    write_if_present(combined$finemap_pip, file.path(out_dir, "finemap_pip_cs.tsv.gz"))
    write_if_present(combined$coloc_summary, file.path(out_dir, "coloc_summary.tsv"))
    write_if_present(combined$coloc_results, file.path(out_dir, "coloc_results.tsv.gz"))
  }
  
  combined
}
