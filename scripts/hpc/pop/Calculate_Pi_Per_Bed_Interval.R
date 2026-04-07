suppressPackageStartupMessages({
  library(data.table)
})

slurm_threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (!is.na(slurm_threads) && slurm_threads > 0L) {
  setDTthreads(slurm_threads)
}

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

find_project_root <- function() {
  script_path <- get_script_path()
  base_dirs <- c(
    Sys.getenv("SLURM_SUBMIT_DIR", unset = NA_character_),
    getwd(),
    if (!is.na(script_path)) dirname(script_path) else NA_character_
  )

  for (base_dir in base_dirs) {
    if (is.na(base_dir) || !nzchar(base_dir) || !dir.exists(base_dir)) {
      next
    }

    candidate <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
    for (depth in 0:4) {
      if (
        dir.exists(file.path(candidate, "data")) &&
        dir.exists(file.path(candidate, "scripts"))
      ) {
        return(candidate)
      }

      parent <- dirname(candidate)
      if (identical(parent, candidate)) {
        break
      }
      candidate <- parent
    }
  }

  stop("Unable to determine project root from SLURM_SUBMIT_DIR, getwd(), or script path.")
}

resolve_path <- function(path, project_root) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(path)
  }

  if (grepl("^/", path)) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(project_root, path), winslash = "/", mustWork = FALSE)
  }
}

project_root <- find_project_root()

default_pestpg_root <- file.path(project_root, "results", "full", "pop", "angsd")
default_bed_dir <- file.path(project_root, "data", "processed", "bed_windows")

parse_cli_args <- function(args, project_root) {
  parsed <- list(
    group = NULL,
    pestpg_dir = NULL,
    bed_dir = default_bed_dir,
    out_dir = NULL,
    skip_empty_pestpg = TRUE
  )

  for (arg in args) {
    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/hpc/pop/Calculate_Pi_Per_Bed_Interval.R [--group=NAME] [--pestpg-dir=PATH]",
        "[--bed-dir=PATH] [--out-dir=PATH] [--keep-empty]\n",
        sep = " "
      )
      quit(save = "no", status = 0)
    } else if (grepl("^--group=", arg)) {
      parsed$group <- sub("^--group=", "", arg)
    } else if (grepl("^--pestpg-dir=", arg)) {
      parsed$pestpg_dir <- resolve_path(sub("^--pestpg-dir=", "", arg), project_root)
    } else if (grepl("^--bed-dir=", arg)) {
      parsed$bed_dir <- resolve_path(sub("^--bed-dir=", "", arg), project_root)
    } else if (grepl("^--out-dir=", arg)) {
      parsed$out_dir <- resolve_path(sub("^--out-dir=", "", arg), project_root)
    } else if (identical(arg, "--keep-empty")) {
      parsed$skip_empty_pestpg <- FALSE
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  if (is.null(parsed$pestpg_dir)) {
    if (!is.null(parsed$group) && nzchar(parsed$group)) {
      parsed$pestpg_dir <- file.path(default_pestpg_root, parsed$group)
    } else {
      parsed$pestpg_dir <- default_pestpg_root
    }
  }

  if (is.null(parsed$out_dir)) {
    parsed$out_dir <- file.path(parsed$pestpg_dir, "pi_by_interval")
  }

  parsed
}

find_first_column <- function(dt, candidates) {
  matches <- candidates[candidates %in% names(dt)]
  if (!length(matches)) {
    return(NA_character_)
  }

  matches[[1L]]
}

read_bed_intervals <- function(path) {
  if (!file.exists(path)) {
    stop("BED file not found: ", path)
  }

  dt <- fread(path, header = FALSE, sep = "\t", showProgress = FALSE)
  if (ncol(dt) < 3L) {
    stop("BED file must have at least 3 columns: ", path)
  }

  setnames(dt, seq_len(ncol(dt)), paste0("V", seq_len(ncol(dt))))

  bed_dt <- data.table(
    chrom = as.character(dt$V1),
    bed_start = as.integer(dt$V2),
    bed_end = as.integer(dt$V3)
  )

  if (ncol(dt) >= 4L) {
    bed_dt[, interval_id := as.character(dt$V4)]
  } else {
    bed_dt[, interval_id := paste0(chrom, ":", bed_start + 1L, "-", bed_end)]
  }

  bed_dt[is.na(interval_id) | !nzchar(interval_id), interval_id := paste0(chrom, ":", bed_start + 1L, "-", bed_end)]
  bed_dt[, interval_bp := bed_end - bed_start]
  bed_dt[, start := bed_start + 1L]
  bed_dt[, end := bed_end]
  bed_dt[, interval_index := .I]

  if (any(is.na(bed_dt$interval_bp) | bed_dt$interval_bp <= 0L)) {
    stop("BED file contains intervals with non-positive width: ", path)
  }

  setkey(bed_dt, chrom, start, end)
  bed_dt
}

read_pestpg <- function(path) {
  if (!file.exists(path)) {
    stop("pestPG file not found: ", path)
  }

  dt <- fread(path, header = TRUE, sep = "\t", showProgress = FALSE)
  if (!nrow(dt)) {
    stop("pestPG file is empty: ", path)
  }

  chrom_col <- find_first_column(dt, c("#Chromo", "Chromo", "Chr", "chr", "chrom"))
  pos_col <- find_first_column(dt, c("Pos", "WinCenter", "midPos", "Position", "pos"))
  pi_col <- find_first_column(dt, c("tP", "thetaPi", "pi", "Pairwise"))
  nsites_col <- find_first_column(dt, c("nSites", "Nsites"))

  if (is.na(chrom_col) || is.na(pos_col) || is.na(pi_col)) {
    stop(
      "Could not find required pestPG columns in ", path,
      ". Available columns: ", paste(names(dt), collapse = ", ")
    )
  }

  pi_values <- as.numeric(dt[[pi_col]])
  if (identical(pi_col, "Pairwise")) {
    pi_values <- exp(pi_values)
  }

  pest_dt <- data.table(
    chrom = as.character(dt[[chrom_col]]),
    pos = as.integer(dt[[pos_col]]),
    pi_site = pi_values
  )

  if (!is.na(nsites_col)) {
    pest_dt[, n_sites := as.numeric(dt[[nsites_col]])]
  } else {
    pest_dt[, n_sites := 1]
  }

  pest_dt <- pest_dt[!is.na(chrom) & !is.na(pos) & !is.na(pi_site)]
  pest_dt[, start := pos]
  pest_dt[, end := pos]
  setkey(pest_dt, chrom, start, end)
  pest_dt
}

calculate_pi_for_pair <- function(bed_path, pestpg_path, out_path = NULL) {
  bed_dt <- read_bed_intervals(bed_path)
  pest_dt <- read_pestpg(pestpg_path)

  joined_dt <- foverlaps(pest_dt, bed_dt, type = "within", nomatch = NA)
  joined_dt <- joined_dt[!is.na(interval_index)]

  if (nrow(joined_dt)) {
    summary_dt <- joined_dt[, .(
      n_pestpg_rows = .N,
      callable_sites = sum(n_sites, na.rm = TRUE),
      pi_total = sum(pi_site, na.rm = TRUE)
    ), by = .(interval_index)]
  } else {
    summary_dt <- data.table(
      interval_index = integer(0),
      n_pestpg_rows = integer(0),
      callable_sites = numeric(0),
      pi_total = numeric(0)
    )
  }

  result_dt <- merge(
    bed_dt[, .(interval_index, chrom, bed_start, bed_end, interval_id, interval_bp)],
    summary_dt,
    by = "interval_index",
    all.x = TRUE,
    sort = FALSE
  )

  result_dt[is.na(n_pestpg_rows), n_pestpg_rows := 0L]
  result_dt[is.na(callable_sites), callable_sites := 0]
  result_dt[is.na(pi_total), pi_total := 0]
  result_dt[, pi_per_bp := pi_total / interval_bp]
  result_dt[, pi_per_callable_site := fifelse(callable_sites > 0, pi_total / callable_sites, NA_real_)]
  result_dt[, interval_index := NULL]
  setcolorder(result_dt, c(
    "chrom", "bed_start", "bed_end", "interval_id", "interval_bp",
    "n_pestpg_rows", "callable_sites", "pi_total", "pi_per_bp", "pi_per_callable_site"
  ))

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    fwrite(result_dt, out_path, sep = "\t")
    message("Wrote interval-level pi: ", out_path)
  }

  result_dt[]
}

calculate_pi_for_matching_files <- function(pestpg_dir, bed_dir, out_dir, skip_empty_pestpg = TRUE) {
  if (!dir.exists(bed_dir)) {
    stop("BED directory not found: ", bed_dir)
  }
  if (!dir.exists(pestpg_dir)) {
    stop("pestPG directory not found: ", pestpg_dir)
  }

  bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
  pestpg_files <- list.files(pestpg_dir, pattern = "\\.pestPG$", full.names = TRUE)

  if (!length(bed_files)) {
    stop("No BED files found in: ", bed_dir)
  }
  if (!length(pestpg_files)) {
    stop("No pestPG files found in: ", pestpg_dir)
  }

  bed_map <- setNames(bed_files, sub("\\.bed$", "", basename(bed_files)))
  pestpg_map <- setNames(pestpg_files, sub("\\.pestPG$", "", basename(pestpg_files)))
  shared_names <- intersect(names(bed_map), names(pestpg_map))

  if (!length(shared_names)) {
    stop(
      "No matching BED / pestPG basenames found between ",
      bed_dir, " and ", pestpg_dir
    )
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- lapply(shared_names, function(name) {
    pestpg_path <- pestpg_map[[name]]
    out_path <- file.path(out_dir, paste0(name, ".interval_pi.tsv"))

    if (skip_empty_pestpg && isTRUE(file.info(pestpg_path)$size == 0)) {
      warning("Skipping empty pestPG file: ", pestpg_path, call. = FALSE, immediate. = TRUE)
      return(NULL)
    }

    message("Processing: ", name)
    calculate_pi_for_pair(
      bed_path = bed_map[[name]],
      pestpg_path = pestpg_path,
      out_path = out_path
    )
  })

  names(results) <- shared_names
  Filter(Negate(is.null), results)
}

run_calculate_pi_per_bed_interval <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- parse_cli_args(args, project_root)

  message("Using pestPG directory: ", parsed$pestpg_dir)
  message("Using BED directory: ", parsed$bed_dir)
  message("Writing results to: ", parsed$out_dir)

  invisible(calculate_pi_for_matching_files(
    pestpg_dir = parsed$pestpg_dir,
    bed_dir = parsed$bed_dir,
    out_dir = parsed$out_dir,
    skip_empty_pestpg = parsed$skip_empty_pestpg
  ))
}

should_autorun <- function() {
  autorun_option <- getOption("calculate_pi_per_bed_interval.autorun", default = NULL)
  if (!is.null(autorun_option)) {
    return(isTRUE(autorun_option))
  }

  file_args <- commandArgs()[grepl("^--file=", commandArgs())]
  if (length(file_args)) {
    script_arg <- normalizePath(sub("^--file=", "", file_args[1L]), winslash = "/", mustWork = FALSE)
    return(!is.na(get_script_path()) && identical(script_arg, get_script_path()))
  }

  interactive()
}

if (should_autorun()) {
  run_calculate_pi_per_bed_interval()
}
