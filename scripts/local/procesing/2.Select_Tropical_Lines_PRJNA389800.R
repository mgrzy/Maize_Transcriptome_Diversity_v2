library(tidyverse)

set.seed(123)

sample_size <- 15L

widiv_meta <- read.csv("data/raw/meta/WiDiv_meta.csv", check.names = FALSE)
filereport <- read.delim(
  "data/raw/meta/filereport_read_run_PRJNA389800.tsv",
  check.names = FALSE
)

normalize_name <- function(x) {
  gsub("[^A-Za-z0-9]", "", toupper(x))
}

tropical_lines <- widiv_meta %>%
  filter(
    Group == "Tropical",
    str_detect(Taxa, "^(CML|NC)")
  ) %>%
  distinct(Taxa) %>%
  arrange(Taxa)

selected_lines <- tropical_lines %>%
  slice_sample(n = min(sample_size, nrow(tropical_lines)))

filereport_normalized <- filereport %>%
  mutate(library_name_normalized = normalize_name(library_name))

matched_fastq <- selected_lines %>%
  mutate(Taxa_normalized = normalize_name(Taxa)) %>%
  mutate(
    matches = map2(
      Taxa,
      Taxa_normalized,
      function(taxa, taxa_normalized) {
        matched_rows <- filereport_normalized %>%
          filter(grepl(taxa_normalized, library_name_normalized, fixed = TRUE)) %>%
          transmute(fastq_ftp)

        if (nrow(matched_rows) == 0L) {
          return(tibble())
        }

        matched_rows
      }
    )
  ) %>%
  select(Taxa, matches) %>%
  unnest(matches) %>%
  separate_wider_delim(
    fastq_ftp,
    delim = ";",
    names = c("fastq_ftp_1", "fastq_ftp_2"),
    too_few = "align_start"
  ) %>%
  select(Taxa, fastq_ftp_1, fastq_ftp_2) %>%
  arrange(Taxa, fastq_ftp_1, fastq_ftp_2)

unmatched_lines <- selected_lines %>%
  anti_join(matched_fastq %>% distinct(Taxa), by = "Taxa")

if (nrow(unmatched_lines) > 0L) {
  warning(
    "No filereport match for: ",
    paste(unmatched_lines$Taxa, collapse = ", ")
  )
}

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  matched_fastq,
  "data/processed/tropical_lines_prjna389800_fastq.csv"
)
