library(tidyverse)

bam_file <- "data/raw/meta/Temperate_BAM.txt"
meta_file <- "data/processed/meta/selected_lines_temperate.csv"
output_dir <- "data/processed/meta/temperate_bam_by_group"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bam_paths <- readLines(bam_file, warn = FALSE) %>%
  discard(~ .x == "")

bam_tbl <- tibble(
  bam_path = bam_paths,
  Taxa = tools::file_path_sans_ext(basename(bam_paths))
)

meta_tbl <- readr::read_csv(meta_file, show_col_types = FALSE) %>%
  select(Taxa, Group)

matched_tbl <- meta_tbl %>%
  inner_join(bam_tbl, by = "Taxa") %>%
  distinct(Group, Taxa, bam_path) %>%
  arrange(Group, Taxa)

unmatched_meta <- meta_tbl %>%
  anti_join(bam_tbl, by = "Taxa")

unmatched_bam <- bam_tbl %>%
  anti_join(meta_tbl, by = "Taxa")

if (nrow(unmatched_meta) > 0L) {
  warning(
    "Taxa present in metadata but missing BAM paths: ",
    paste(unmatched_meta$Taxa, collapse = ", ")
  )
}

if (nrow(unmatched_bam) > 0L) {
  warning(
    "BAM paths without a matching Taxa in metadata: ",
    paste(unmatched_bam$Taxa, collapse = ", ")
  )
}

sanitize_group_name <- function(group_name) {
  group_name %>%
    str_replace_all("[^A-Za-z0-9._-]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

matched_tbl %>%
  group_by(Group) %>%
  group_walk(~ {
    output_file <- file.path(
      output_dir,
      paste0(sanitize_group_name(.y$Group), "_BAM.txt")
    )
    writeLines(.x$bam_path, output_file)
  })

summary_tbl <- matched_tbl %>%
  count(Group, name = "bam_count") %>%
  arrange(Group)

readr::write_csv(
  summary_tbl,
  file.path(output_dir, "group_bam_counts.csv")
)
