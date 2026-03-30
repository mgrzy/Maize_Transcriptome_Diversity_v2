library(tidyverse)

set.seed(123)

lines_per_group <- 15L
balanced_groups <- c("SS", "NSS")
mandatory_taxa <- c("B73", "PH207")

meta <- readxl::read_xlsx("data/raw/meta/supplemetnal_tables_final.xlsx", skip = 1)
widiv_meta <- read.csv("data/raw/meta/WiDiv_meta.csv", check.names = FALSE)
rna_meta <- read.csv(
  "data/raw/meta/RNA_seq_sample_meta_data_selected.csv",
  check.names = FALSE
)
filereport <- read.delim(
  "data/raw/meta/filereport_read_run_PRJNA661271.tsv",
  check.names = FALSE
) %>%
  select(sample_accession, fastq_ftp) %>%
  separate_wider_delim(
    fastq_ftp,
    delim = ";",
    names = c("fastq_ftp_1", "fastq_ftp_2"),
    too_few = "align_start"
  )

colnames(meta)[1] <- "Taxa"

take_prioritized <- function(df, n_to_take) {
  if (n_to_take <= 0L || nrow(df) == 0L) {
    return(df[0, , drop = FALSE])
  }

  df %>%
    mutate(
      priority_bucket = case_when(
        is_mandatory ~ 3L,
        in_rna ~ 2L,
        TRUE ~ 1L
      ),
      priority_rand = runif(n())
    ) %>%
    arrange(desc(priority_bucket), priority_rand, Taxa) %>%
    slice_head(n = min(n_to_take, nrow(.))) %>%
    select(-priority_bucket, -priority_rand)
}

build_era_targets <- function(df, total_n) {
  era_counts <- df %>%
    filter(!is.na(Era)) %>%
    count(Era, name = "available")

  if (nrow(era_counts) == 0L) {
    return(tibble(Era = character(), target = integer()))
  }

  targets <- rep(total_n %/% nrow(era_counts), nrow(era_counts))
  remainder <- total_n %% nrow(era_counts)

  if (remainder > 0L) {
    extra_idx <- order(era_counts$available, era_counts$Era, decreasing = TRUE)[seq_len(remainder)]
    targets[extra_idx] <- targets[extra_idx] + 1L
  }

  tibble(Era = era_counts$Era, target = targets)
}

select_group_lines <- function(df) {
  total_n <- min(lines_per_group, nrow(df))

  selected <- take_prioritized(df %>% filter(is_mandatory), total_n)
  remaining <- anti_join(df, selected, by = "row_id")

  if (nrow(selected) < total_n && unique(df$Group) %in% balanced_groups) {
    era_targets <- build_era_targets(df, total_n) %>%
      left_join(selected %>% count(Era, name = "already_selected"), by = "Era") %>%
      mutate(
        already_selected = replace_na(already_selected, 0L),
        need = pmax(target - already_selected, 0L)
      )

    for (idx in seq_len(nrow(era_targets))) {
      if (nrow(selected) >= total_n) {
        break
      }

      era_name <- era_targets$Era[idx]
      need_n <- min(era_targets$need[idx], total_n - nrow(selected))

      era_pick <- take_prioritized(
        remaining %>% filter(Era == era_name),
        need_n
      )

      selected <- bind_rows(selected, era_pick)
      remaining <- anti_join(remaining, era_pick, by = "row_id")
    }
  }

  if (nrow(selected) < total_n) {
    fill_pick <- take_prioritized(remaining, total_n - nrow(selected))
    selected <- bind_rows(selected, fill_pick)
  }

  selected
}

rna_taxa <- rna_meta %>%
  pull(VCFname) %>%
  na.omit() %>%
  unique()

a <- meta %>%
  left_join(widiv_meta, by = "Taxa") %>%
  left_join(filereport, by = c("BioSample...6" = "sample_accession")) %>%
  mutate(
    row_id = row_number(),
    Group = str_squish(Group),
    Era = na_if(str_squish(Era), ""),
    in_rna = Taxa %in% rna_taxa,
    is_mandatory = Taxa %in% mandatory_taxa,
    has_fastq_match = if_any(
      starts_with("fastq_ftp_"),
      ~ !is.na(.x) & .x != ""
    )
  ) %>%
  filter(
    !is.na(Group),
    !str_to_lower(Group) %in% c("sweet corn", "popcorn", "landrace"),
    has_fastq_match | is_mandatory
  )

selected_lines <- a %>%
  group_by(Group) %>%
  group_split() %>%
  map_dfr(select_group_lines) %>%
  arrange(Group, desc(is_mandatory), desc(in_rna), Era, Taxa) %>%
  ungroup()

selection_summary <- selected_lines %>%
  count(Group, Era, name = "selected_n") %>%
  arrange(Group, Era)

selection_priority_summary <- selected_lines %>%
  summarise(
    total_selected = n(),
    selected_in_rna = sum(in_rna),
    mandatory_selected = sum(is_mandatory),
    fastq_matched = sum(has_fastq_match)
  )

readr::write_csv(selected_lines, "data/processed/selected_lines.csv")

