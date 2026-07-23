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

  stop("Could not resolve the script path.")
}

script_path <- get_script_path()
target_path <- file.path(dirname(script_path), "X.Select_Min.P.R")

source(target_path, local = FALSE)
