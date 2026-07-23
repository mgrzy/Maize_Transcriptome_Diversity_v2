suppressPackageStartupMessages({
  library(data.table)
})

get_script_path <- function() {
  script_path_arg <- "--file="
  script_path_match <- commandArgs()[grepl(paste0("^", script_path_arg), commandArgs())]

  if (length(script_path_match)) {
    return(normalizePath(sub(script_path_arg, "", script_path_match[1L]), winslash = "/", mustWork = TRUE))
  }

  source_refs <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  if (length(source_refs)) {
    return(normalizePath(source_refs[[length(source_refs)]], winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

resolve_project_root <- function() {
  script_path <- get_script_path()

  if (!is.na(script_path)) {
    return(normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE))
  }

  normalizePath(".", winslash = "/", mustWork = TRUE)
}

find_first_existing <- function(paths) {
  paths <- paths[file.exists(paths)]

  if (!length(paths)) {
    return(NA_character_)
  }

  normalizePath(paths[[1L]], winslash = "/", mustWork = TRUE)
}

find_first_column <- function(dt, candidates) {
  matches <- candidates[candidates %in% names(dt)]

  if (!length(matches)) {
    return(NA_character_)
  }

  matches[[1L]]
}

choose_column <- function(dt, explicit, candidates, label) {
  if (!is.na(explicit)) {
    if (!explicit %in% names(dt)) {
      stop(
        "Requested ", label, " column '", explicit, "' was not found. Available columns: ",
        paste(names(dt), collapse = ", ")
      )
    }

    return(explicit)
  }

  detected <- find_first_column(dt, candidates)
  if (!is.na(detected)) {
    return(detected)
  }

  stop(
    "Could not detect a ", label, " column. Available columns: ",
    paste(names(dt), collapse = ", ")
  )
}

parse_args <- function(args, defaults) {
  config <- defaults

  if (!length(args)) {
    return(config)
  }

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-h", "--help")) {
      config$help <- TRUE
      i <- i + 1L
      next
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
      stop(
        "Unknown argument '--", gsub("_", "-", key), "'. Supported arguments: ",
        paste(paste0("--", gsub("_", "-", names(config))), collapse = ", ")
      )
    }

    config[[key]] <- value
    i <- i + 1L
  }

  config
}

print_usage <- function(defaults) {
  cat(
    paste(
      "Select the SNP with the minimum p-value for each gene.",
      "",
      "Usage:",
      "  Rscript scripts/local/analysis/X.Select_Min.P.R [options]",
      "",
      "Options:",
      paste0("  --input ", defaults$input),
      paste0("  --output ", defaults$output),
      "  --gene-col Gene",
      "  --p-col p",
      "  --help",
      "",
      "The output keeps one row per gene and preserves all original columns.",
      sep = "\n"
    )
  )
}

select_min_p_per_gene <- function(dt, gene_col, p_col) {
  dt <- copy(dt)
  dt[, (gene_col) := trimws(as.character(get(gene_col)))]
  dt[, (p_col) := suppressWarnings(as.numeric(get(p_col)))]

  dt <- dt[!is.na(get(gene_col)) & nzchar(get(gene_col))]
  dt <- dt[!is.na(get(p_col)) & is.finite(get(p_col))]

  if (!nrow(dt)) {
    stop("No rows with valid gene IDs and p-values were found after filtering.")
  }

  tie_break_cols <- intersect(
    c("Chr", "chr", "CHR", "CHROM", "chrom", "BP", "POS", "pos", "SNP", "snp"),
    names(dt)
  )
  tie_break_cols <- setdiff(tie_break_cols, c(gene_col, p_col))

  stats_dt <- dt[
    ,
    .(
      n_snps_tested_for_gene = .N,
      tie_min_p_count = sum(get(p_col) == min(get(p_col)))
    ),
    by = gene_col
  ]

  setorderv(dt, c(gene_col, p_col, tie_break_cols))
  selected_idx <- dt[, .I[1L], by = gene_col]$V1
  selected_dt <- merge(dt[selected_idx], stats_dt, by = gene_col, sort = FALSE)

  selected_dt
}

project_root <- resolve_project_root()

default_input <- find_first_existing(c(
  file.path(project_root, "..", "Transcriptome_diversity", "results", "eQTL", "eQTL_2Mb_DTScov.txt"),
  file.path(project_root, "..", "Transcriptome_diversity", "results", "eQTL", "eQTL_p001.txt"),
  file.path(project_root, "results", "full", "eQTL", "eQTL_2Mb_DTScov.txt"),
  file.path(project_root, "results", "full", "eQTL", "eQTL_p001.txt")
))

default_output <- file.path(project_root, "results", "full", "eQTL", "min_p_snp_per_gene.tsv")

defaults <- list(
  input = default_input,
  output = default_output,
  gene_col = NA_character_,
  p_col = NA_character_,
  help = FALSE
)

config <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

if (isTRUE(config$help)) {
  print_usage(defaults)
  quit(save = "no", status = 0L)
}

if (is.na(config$input) || !nzchar(config$input)) {
  stop(
    "No default input file was found. Provide one with --input. Looked for files under ../Transcriptome_diversity/results/eQTL and results/full/eQTL."
  )
}

input_path <- normalizePath(config$input, winslash = "/", mustWork = TRUE)
output_path <- normalizePath(config$output, winslash = "/", mustWork = FALSE)
output_dir <- dirname(output_path)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

message("Reading eQTL table: ", input_path)
eqtl_dt <- fread(input_path, showProgress = FALSE)

if (!nrow(eqtl_dt)) {
  stop("Input file is empty: ", input_path)
}

gene_col <- choose_column(eqtl_dt, config$gene_col, c("Gene", "gene", "Probe", "probe"), "gene")
p_col <- choose_column(eqtl_dt, config$p_col, c("p", "P", "pvalue", "p_value", "P_VALUE"), "p-value")

message("Using gene column: ", gene_col)
message("Using p-value column: ", p_col)

selected_dt <- select_min_p_per_gene(eqtl_dt, gene_col = gene_col, p_col = p_col)

fwrite(selected_dt, file = output_path, sep = "\t", quote = FALSE, na = "NA")

message("Selected ", nrow(selected_dt), " SNPs with minimum p-value across ", uniqueN(selected_dt[[gene_col]]), " genes.")
message("Saved output to: ", output_path)
