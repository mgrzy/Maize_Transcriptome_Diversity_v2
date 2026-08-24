suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

get_script_path <- function() {
  script_arg <- "--file="
  script_match <- commandArgs()[grepl(paste0("^", script_arg), commandArgs())]

  if (length(script_match)) {
    return(normalizePath(sub(script_arg, "", script_match[1L]), winslash = "/", mustWork = TRUE))
  }

  source_refs <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  if (length(source_refs)) {
    return(normalizePath(source_refs[[length(source_refs)]], winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

find_project_root <- function() {
  script_path <- get_script_path()
  base_dirs <- c(getwd(), if (!is.na(script_path)) dirname(script_path) else NA_character_)

  for (base_dir in base_dirs) {
    if (is.na(base_dir) || !nzchar(base_dir) || !dir.exists(base_dir)) {
      next
    }

    candidate <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
    for (depth in 0:4) {
      if (dir.exists(file.path(candidate, "scripts")) && dir.exists(file.path(candidate, "results"))) {
        return(candidate)
      }

      parent <- dirname(candidate)
      if (identical(parent, candidate)) {
        break
      }
      candidate <- parent
    }
  }

  stop("Unable to determine project root.")
}

resolve_input_path <- function(path, project_root) {
  if (is.null(path)) {
    return(NULL)
  }

  if (grepl("^/", path)) {
    normalizePath(path, winslash = "/", mustWork = TRUE)
  } else {
    normalizePath(file.path(project_root, path), winslash = "/", mustWork = TRUE)
  }
}

resolve_output_path <- function(path, project_root) {
  if (is.null(path)) {
    return(NULL)
  }

  if (grepl("^/", path)) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(project_root, path), winslash = "/", mustWork = FALSE)
  }
}

load_pi_objects <- function(pi_dir = NULL, assign_global = FALSE) {
  project_root <- find_project_root()

  if (is.null(pi_dir)) {
    pi_dir <- file.path(project_root, "results", "full", "pop", "angsd", "pi")
  } else {
    pi_dir <- resolve_input_path(pi_dir, project_root)
  }

  file_paths <- list.files(
    path = pi_dir,
    pattern = "\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(file_paths)) {
    stop("No TSV files found in: ", pi_dir)
  }

  pi_raw_tbl <- map_dfr(file_paths, function(path) {
    read_tsv(path, show_col_types = FALSE) %>%
      mutate(
        region = str_remove(basename(path), "\\.interval_pi\\.tsv$"),
        subpopulation = str_remove(basename(dirname(path)), "^pi_")
      )
  })

  pi_summary_tbl <- pi_raw_tbl %>%
    group_by(subpopulation, region) %>%
    summarise(
      median_pi_per_callable_site = median(pi_per_callable_site, na.rm = TRUE),
      .groups = "drop"
    )

  pi_trop_tbl <- pi_summary_tbl %>%
    filter(str_to_lower(subpopulation) == "trop") %>%
    transmute(
      region,
      trop_median_pi_per_callable_site = median_pi_per_callable_site
    )

  pi_summary_tbl <- pi_summary_tbl %>%
    left_join(pi_trop_tbl, by = "region") %>%
    mutate(
      median_pi_per_callable_site_relative_to_trop =
        median_pi_per_callable_site / trop_median_pi_per_callable_site
    )

  pi_tbl <- pi_raw_tbl %>%
    left_join(pi_summary_tbl, by = c("subpopulation", "region")) %>%
    select(
      region,
      subpopulation,
      median_pi_per_callable_site,
      trop_median_pi_per_callable_site,
      median_pi_per_callable_site_relative_to_trop,
      everything()
    )

  pi_objects <- list(
    pi_raw_tbl = pi_raw_tbl,
    pi_summary_tbl = pi_summary_tbl,
    pi_trop_tbl = pi_trop_tbl,
    pi_tbl = pi_tbl
  )

  if (isTRUE(assign_global)) {
    list2env(pi_objects, envir = .GlobalEnv)
  }

  pi_objects
}

load_pi_tables <- function(pi_dir = NULL) {
  load_pi_objects(pi_dir = pi_dir, assign_global = FALSE)$pi_tbl
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  project_root <- find_project_root()
  pi_dir <- NULL
  out_path <- NULL

  for (arg in args) {
    if (grepl("^--pi-dir=", arg)) {
      pi_dir <- sub("^--pi-dir=", "", arg)
    } else if (grepl("^--out=", arg)) {
      out_path <- sub("^--out=", "", arg)
    } else if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/local/analysis/2.Load_Pi_Files.R [--pi-dir=PATH] [--out=PATH]\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  pi_objects <- load_pi_objects(pi_dir = pi_dir, assign_global = FALSE)
  pi_tbl <- pi_objects$pi_tbl

  if (!is.null(out_path)) {
    out_path <- resolve_output_path(out_path, project_root)
    write_tsv(pi_tbl, out_path)
    message("Wrote combined table to: ", out_path)
  } else {
    message(
      "Loaded ", nrow(pi_tbl), " rows from ",
      n_distinct(tibble(region = pi_tbl$region, subpopulation = pi_tbl$subpopulation)),
      " region/subpopulation combinations."
    )
    print(pi_tbl)
  }
} else if (interactive()) {
  load_pi_objects(assign_global = TRUE)
  message("Created objects in .GlobalEnv: pi_raw_tbl, pi_summary_tbl, pi_trop_tbl, pi_tbl")
}
