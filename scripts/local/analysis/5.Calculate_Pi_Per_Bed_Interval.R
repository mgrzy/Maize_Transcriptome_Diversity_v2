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

project_root <- resolve_project_root()

default_pestpg_dir <- file.path(project_root, "results", "full", "pop", "angsd")
default_bed_dir <- file.path(project_root, "results", "summary", "bed_windows")
default_out_dir <- file.path(project_root, "results", "full", "pop", "angsd", "pi_by_interval")

calculate_pi_per_bed_interval_config <- list(
  pestpg_dir = default_pestpg_dir,
  bed_dir = default_bed_dir,
  out_dir = default_out_dir
)

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

  chrom_col <- find_first_column(dt, c("Chr", "chr", "chrom"))
  pos_col <- find_first_column(dt, c("WinCenter", "midPos", "Position", "pos"))
  pi_col <- find_first_column(dt, c("tP", "thetaPi", "pi"))
  nsites_col <- find_first_column(dt, c("nSites", "Nsites"))

  if (is.na(chrom_col) || is.na(pos_col) || is.na(pi_col)) {
    stop(
      "Could not find required pestPG columns in ", path,
      ". Available columns: ", paste(names(dt), collapse = ", ")
    )
  }

  pest_dt <- data.table(
    chrom = as.character(dt[[chrom_col]]),
    pos = as.integer(dt[[pos_col]]),
    pi_site = as.numeric(dt[[pi_col]])
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

calculate_pi_for_matching_files <- function(
  config = calculate_pi_per_bed_interval_config,
  skip_empty_pestpg = TRUE
) {
  config <- modifyList(calculate_pi_per_bed_interval_config, config)

  if (!dir.exists(config$bed_dir)) {
    stop("BED directory not found: ", config$bed_dir)
  }
  if (!dir.exists(config$pestpg_dir)) {
    stop("pestPG directory not found: ", config$pestpg_dir)
  }

  bed_files <- list.files(config$bed_dir, pattern = "\\.bed$", full.names = TRUE)
  pestpg_files <- list.files(config$pestpg_dir, pattern = "\\.pestPG$", full.names = TRUE)

  if (!length(bed_files)) {
    stop("No BED files found in: ", config$bed_dir)
  }
  if (!length(pestpg_files)) {
    stop("No pestPG files found in: ", config$pestpg_dir)
  }

  bed_map <- setNames(bed_files, sub("\\.bed$", "", basename(bed_files)))
  pestpg_map <- setNames(pestpg_files, sub("\\.pestPG$", "", basename(pestpg_files)))
  shared_names <- intersect(names(bed_map), names(pestpg_map))

  if (!length(shared_names)) {
    stop(
      "No matching BED / pestPG basenames found between ",
      config$bed_dir, " and ", config$pestpg_dir
    )
  }

  dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- lapply(shared_names, function(name) {
    pestpg_path <- pestpg_map[[name]]
    out_path <- file.path(config$out_dir, paste0(name, ".interval_pi.tsv"))

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
  results <- Filter(Negate(is.null), results)

  invisible(results)
}

run_calculate_pi_per_bed_interval <- function(
  config = calculate_pi_per_bed_interval_config,
  skip_empty_pestpg = TRUE,
  assign_global = interactive()
) {
  results <- calculate_pi_for_matching_files(
    config = config,
    skip_empty_pestpg = skip_empty_pestpg
  )

  if (assign_global) {
    assign("pi_by_interval_results", results, envir = .GlobalEnv)
    message("Saved results list to pi_by_interval_results")
  }

  invisible(results)
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
