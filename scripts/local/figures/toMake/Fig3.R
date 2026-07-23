suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(data.table)
  library(bigmemory)
  library(ggrepel)
})

# Match Fig1 font sizing
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c(
  "#66C5CC", # teal
  "#F6CF71", # yellow
  "#DCB0F2", # lavender
  "#87C55F", # green
  "#9EB9F3", # light blue
  "#C3B091", # beige
  "#F49C9C", # soft coral
  "#FFA56A"  # warm apricot
)

hapPal <- myPal[1:4]

pip_df <- fread("results/SusieR2/PIP.csv") %>% as_tibble()
snp <- attach.big.matrix("data/geno/rMVP/Div639.geno.desc")
ind <- read.table("data/geno/rMVP/Div639.geno.ind", stringsAsFactors = FALSE) %>%
  as_tibble(.name_repair = "minimal") %>%
  setNames("Taxa")
map <- fread("data/geno/rMVP/Div639.geno.map") %>% as_tibble()
expr_counts <- fread("data/RNA/work/Filtered_Noralised_Counts.csv") %>% as_tibble()

gff_path <- "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
gff <- fread(
  cmd = paste("gzip -dc", shQuote(gff_path), "| grep -v '^#'"),
  sep = "\t",
  header = FALSE,
  col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
) %>%
  as_tibble()

# =========================
# Gene: Zm00001eb146510
# =========================

gene1_id <- "Zm00001eb146510"
gene1_hap_threshold <- 45

gene1_gwas_path <- file.path("results/eQTL", gene1_id)
gene1_gwas_df <- fread(gene1_gwas_path) %>% as_tibble()

gene1_fm_df <- pip_df %>%
  filter(gene == gene1_id) %>%
  filter(!is.na(POS), !is.na(variable_prob)) %>%
  mutate(variable_prob = pmin(variable_prob, 1))

gene1_cs_count <- n_distinct(gene1_fm_df$cs)
gene1_top_per_cs <- 1
gene1_fm_max <- gene1_fm_df %>%
  group_by(cs) %>%
  slice_max(variable_prob, n = gene1_top_per_cs, with_ties = FALSE) %>%
  arrange(desc(variable_prob)) %>%
  ungroup()

gene1_highlight_snps <- gene1_fm_max %>% pull(SNP)
gene1_highlight_levels <- c("Other", gene1_highlight_snps)

gene1_gwas_df <- gene1_gwas_df %>%
  mutate(
    logp = -log10(p),
    highlight = if_else(SNP %in% gene1_highlight_snps, SNP, "Other"),
    highlight = factor(highlight, levels = gene1_highlight_levels)
  )

gene1_fm_df <- gene1_fm_df %>%
  mutate(
    highlight = if_else(SNP %in% gene1_highlight_snps, SNP, "Other"),
    highlight = factor(highlight, levels = gene1_highlight_levels)
  )

gene1_bp_range <- range(gene1_gwas_df$BP, na.rm = TRUE)
gene1_window_padding <- 5e4
gene1_window_start <- max(gene1_bp_range[1] - gene1_window_padding, 0)
gene1_window_end <- gene1_bp_range[2] + gene1_window_padding

gene1_position_breaks <- pretty(c(gene1_window_start, gene1_window_end), n = 4)
gene1_position_labels <- label_number(scale = 1e-6, accuracy = 0.01)

gene1_chr_values <- gene1_gwas_df %>%
  pull(Chr) %>%
  na.omit() %>%
  unique()

gene1_target_seqids <- c(paste0("chr", gene1_chr_values), paste0("Chr", gene1_chr_values), as.character(gene1_chr_values))

gene1_gff_genes <- gff %>%
  filter(
    seqid %in% gene1_target_seqids,
    type == "gene",
    start <= gene1_window_end,
    end >= gene1_window_start
  ) %>%
  mutate(
    gene_id = str_extract(attributes, "(?<=ID=)[^;]+"),
    strand_label = factor(if_else(strand == "+", "+", "-"), levels = c("+", "-")),
    forward = strand == "+"
  ) %>%
  arrange(start)

gene1_selected_snps <- gene1_fm_max$SNP
gene1_selected_idx <- match(gene1_selected_snps, map$SNP)
gene1_selected_snps <- gene1_selected_snps[!is.na(gene1_selected_idx)]
gene1_selected_idx <- gene1_selected_idx[!is.na(gene1_selected_idx)]

gene1_expr_gene <- ind %>%
  left_join(expr_counts %>% select(Taxa, Expr = all_of(gene1_id)), by = "Taxa")

gene1_geno_df <- as.data.frame(snp[, gene1_selected_idx, drop = FALSE])
colnames(gene1_geno_df) <- gene1_selected_snps

gene1_expr_gene <- bind_cols(gene1_expr_gene, gene1_geno_df)

gene1_allele_lookup <- tibble(
  SNP = gene1_selected_snps,
  REF = map$REF[gene1_selected_idx],
  ALT = map$ALT[gene1_selected_idx]
)

gene1_hap_df <- gene1_expr_gene %>%
  select(all_of(gene1_selected_snps)) %>%
  as.matrix()

gene1_hap_code <- apply(gene1_hap_df, 1, function(row) {
  if (any(is.na(row))) {
    return(NA_character_)
  }
  paste0(row, collapse = "")
})

gene1_hap_nuc <- apply(gene1_hap_df, 1, function(row) {
  if (any(is.na(row))) {
    return(NA_character_)
  }
  alleles <- vapply(seq_along(row), function(idx) {
    dosage <- row[idx]
    ref_allele <- gene1_allele_lookup$REF[idx]
    alt_allele <- gene1_allele_lookup$ALT[idx]
    if (is.na(dosage)) {
      return(NA_character_)
    }
    if (dosage == 0) {
      return(ref_allele)
    }
    if (dosage == 2) {
      return(alt_allele)
    }
    return(NA_character_)
  }, character(1))
  if (any(is.na(alleles))) {
    return(NA_character_)
  }
  paste0(alleles, collapse = "")
})

gene1_combo_label <- apply(gene1_hap_df, 1, function(row) {
  if (any(is.na(row))) {
    return(NA_character_)
  }
  alleles <- vapply(seq_along(row), function(idx) {
    dosage <- row[idx]
    ref_allele <- gene1_allele_lookup$REF[idx]
    alt_allele <- gene1_allele_lookup$ALT[idx]
    if (is.na(dosage)) {
      return(NA_character_)
    }
    if (dosage == 0) {
      return(ref_allele)
    }
    if (dosage == 2) {
      return(alt_allele)
    }
    if (dosage == 1) {
      return(paste0(ref_allele, "/", alt_allele))
    }
    NA_character_
  }, character(1))
  if (any(is.na(alleles))) {
    return(NA_character_)
  }
  paste0(alleles, collapse = "-")
})

gene1_expr_gene <- gene1_expr_gene %>%
  mutate(
    Hap_code = gene1_hap_code,
    Hap = gene1_hap_nuc,
    Allele_combo = gene1_combo_label
  )

gene1_allele_combo_counts <- gene1_expr_gene %>%
  filter(!is.na(Allele_combo)) %>%
  count(Allele_combo, sort = TRUE)

gene1_hap_counts <- gene1_expr_gene %>%
  filter(!is.na(Hap)) %>%
  count(Hap, sort = TRUE)

gene1_hap_selected <- gene1_hap_counts %>%
  filter(n >= gene1_hap_threshold)

gene1_expr_hap45 <- gene1_expr_gene %>%
  filter(Hap %in% gene1_hap_selected$Hap)

gene1_hap_levels <- gene1_expr_hap45 %>%
  group_by(Hap) %>%
  summarise(median_expr = median(Expr, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_expr) %>%
  pull(Hap)

gene1_label_lookup <- tibble(
  Original = gene1_hap_levels,
  Short = paste0("Hap", seq_along(gene1_hap_levels))
)

gene1_expr_violin_df <- gene1_expr_hap45 %>%
  mutate(
    ViolinLabel = factor(
      gene1_label_lookup$Short[match(Hap, gene1_label_lookup$Original)],
      levels = gene1_label_lookup$Short
    )
  )

gene1_violin_lookup <- gene1_expr_violin_df %>%
  distinct(Allele_combo, ViolinLabel)

gene1_allele_combo_counts <- gene1_allele_combo_counts %>%
  left_join(gene1_violin_lookup, by = "Allele_combo") %>%
  mutate(
    HistGroup = if_else(is.na(ViolinLabel), "Other", as.character(ViolinLabel)),
    HistGroup = factor(HistGroup, levels = c(as.character(gene1_label_lookup$Short), "Other")),
    HapNumber = if_else(is.na(ViolinLabel), NA_character_, str_remove(as.character(ViolinLabel), "^Hap"))
  )

gene1_hap_counts_45 <- gene1_hap_selected %>%
  mutate(Display = gene1_label_lookup$Short[match(Hap, gene1_label_lookup$Original)])

gene1_hap_violin <- gene1_expr_violin_df %>%
  ggplot(aes(x = ViolinLabel, y = Expr, fill = ViolinLabel)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, colour = "black") +
  scale_fill_manual(values = rep(hapPal, length.out = length(gene1_label_lookup$Short))) +
  labs(
    x = "Haplotype",
    y = paste(gene1_id, "expression")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gene1_mapping_caption <- paste(gene1_label_lookup$Short, gene1_label_lookup$Original, sep = " = ", collapse = "; ")

gene1_hap_bar <- gene1_hap_counts_45 %>%
  mutate(Display = factor(Display, levels = gene1_label_lookup$Short)) %>%
  ggplot(aes(x = Display, y = n)) +
  geom_col(fill = myPal[1], width = 0.8) +
  geom_text(aes(label = n), vjust = -0.3, size = 2.1) +
  scale_y_continuous("Individuals", expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Haplotype",
    subtitle = gene1_mapping_caption
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gene1_hist_palette <- c(
  setNames(rep(hapPal, length.out = length(gene1_label_lookup$Short)), gene1_label_lookup$Short),
  Other = myPal[2]
)

gene1_allele_hist_data <- gene1_allele_combo_counts %>%
  mutate(Allele_combo = factor(Allele_combo, levels = Allele_combo))

gene1_allele_hist_plot <- ggplot(gene1_allele_hist_data, aes(x = Allele_combo, y = n, fill = HistGroup)) +
  geom_col(width = 0.75, colour = NA) +
  geom_text(aes(label = n), vjust = -0.3, size = 1.8, colour = "black") +
  geom_text(
    data = gene1_allele_hist_data %>% filter(!is.na(HapNumber)),
    aes(y = n / 2, label = HapNumber, colour = HistGroup),
    fontface = "bold",
    size = 2.0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = gene1_hist_palette, guide = "none") +
  scale_colour_manual(values = gene1_hist_palette, guide = "none") +
  scale_y_continuous("Individuals", expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Allele combination across fine-mapped SNPs",
    y = "Individuals"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gene1_gene_segments <- gene1_gff_genes %>%
  mutate(
    colour_group = if_else(gene_id == gene1_id, "target", "other"),
    y = case_when(
      colour_group == "target" ~ 1.15,
      forward ~ 1,
      TRUE ~ 0
    ),
    x = if_else(forward, start, end),
    xend = if_else(forward, end, start)
  ) %>%
  select(x, xend, y, colour_group)

gene1_gene_bg <- gene1_gene_segments %>% filter(colour_group == "other")
gene1_gene_target <- gene1_gene_segments %>% filter(colour_group == "target")

gene1_gene_track <- ggplot() +
  geom_hline(yintercept = c(0, 1), colour = "#E0E0E0", linewidth = 0.25) +
  geom_segment(
    data = gene1_gene_bg,
    aes(x = x, xend = xend, y = y, yend = y),
    colour = "#1A1A1A",
    linewidth = 0.45,
    lineend = "round",
    arrow = arrow(length = grid::unit(1.0, "mm"), type = "open", angle = 18)
  ) +
  geom_segment(
    data = gene1_gene_target,
    aes(x = x, xend = xend, y = y, yend = y),
    colour = "#D62728",
    linewidth = 0.85,
    lineend = "round",
    arrow = arrow(length = grid::unit(1.8, "mm"), type = "open", angle = 18)
  ) +
  coord_cartesian(xlim = c(gene1_window_start, gene1_window_end), ylim = c(-0.4, 1.4)) +
  scale_x_continuous(
    name = "Chromosome 3 (Mb)",
    breaks = gene1_position_breaks,
    labels = gene1_position_labels
  ) +
  scale_y_continuous(
    breaks = c(1, 0),
    labels = c("+", "-"),
    expand = expansion(add = 0)
  ) +
  labs(y = NULL) +
  theme(
    axis.text.y = element_text(size = base_font_size - 2),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

gene1_highlight_palette <- c("Other" = "#BDBDBD")
gene1_highlight_palette[gene1_highlight_snps] <- rep(myPal[5:8], length.out = length(gene1_highlight_snps))

gene1_gwas_top <- gene1_gwas_df %>% filter(highlight != "Other")
gene1_gwas_bg <- gene1_gwas_df %>% filter(highlight == "Other")

gene1_fm_top <- gene1_fm_df %>% filter(highlight != "Other")
gene1_fm_bg <- gene1_fm_df %>% filter(highlight == "Other")

gene1_gwas_plot <- ggplot(gene1_gwas_df, aes(x = BP, y = logp)) +
  geom_point(
    data = gene1_gwas_bg,
    colour = gene1_highlight_palette["Other"],
    size = 1,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  geom_point(
    data = gene1_gwas_top,
    aes(colour = highlight),
    size = 2.4,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = gene1_highlight_palette,
    breaks = gene1_highlight_snps,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(breaks = gene1_position_breaks, labels = gene1_position_labels) +
  coord_cartesian(xlim = c(gene1_window_start, gene1_window_end)) +
  labs(
    x = NULL,
    y = expression(-log[10](p))
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

gene1_fm_plot <- ggplot(gene1_fm_df, aes(x = POS, y = variable_prob)) +
  geom_point(
    data = gene1_fm_bg,
    colour = gene1_highlight_palette["Other"],
    size = 1.3,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  geom_point(
    data = gene1_fm_top,
    aes(colour = highlight),
    size = 2.6,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = gene1_fm_top,
    aes(label = SNP),
    colour = "black",
    size = 2.1,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0.05,
    seed = 123
  ) +
  scale_colour_manual(
    values = gene1_highlight_palette,
    breaks = gene1_highlight_snps,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(breaks = gene1_position_breaks, labels = gene1_position_labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(gene1_window_start, gene1_window_end)) +
  labs(
    x = NULL,
    y = "PIP"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

g1 <- (gene1_gwas_plot / gene1_fm_plot / gene1_gene_track) + plot_layout(heights = c(0.8, 0.8, 0.2)) | gene1_hap_violin

design <- "
AB
CB
DB
"


g1 <- wrap_plots(
  A = gene1_gwas_plot,
  C = gene1_fm_plot,
  D = gene1_gene_track,
  B = gene1_hap_violin,
  design = design
) +
  plot_layout(widths  = c(4, 1),
              heights = c(0.8, 0.8, 0.2))

# =========================
# Gene: Zm00001eb403750
# =========================

gene2_id <- "Zm00001eb403750"
gene2_hap_threshold <- 30

gene2_gwas_path <- file.path("results/eQTL", gene2_id)
gene2_gwas_df <- fread(gene2_gwas_path) %>% as_tibble()

gene2_fm_df <- pip_df %>%
  filter(gene == gene2_id) %>%
  filter(!is.na(POS), !is.na(variable_prob)) %>%
  mutate(variable_prob = pmin(variable_prob, 1))

gene2_cs_count <- n_distinct(gene2_fm_df$cs)
gene2_top_per_cs <- 1
gene2_fm_max <- gene2_fm_df %>%
  group_by(cs) %>%
  slice_max(variable_prob, n = gene2_top_per_cs, with_ties = FALSE) %>%
  arrange(desc(variable_prob)) %>%
  ungroup()

gene2_highlight_snps <- gene2_fm_max %>% pull(SNP)
gene2_highlight_levels <- c("Other", gene2_highlight_snps)

gene2_gwas_df <- gene2_gwas_df %>%
  mutate(
    logp = -log10(p),
    highlight = if_else(SNP %in% gene2_highlight_snps, SNP, "Other"),
    highlight = factor(highlight, levels = gene2_highlight_levels)
  )

gene2_fm_df <- gene2_fm_df %>%
  mutate(
    highlight = if_else(SNP %in% gene2_highlight_snps, SNP, "Other"),
    highlight = factor(highlight, levels = gene2_highlight_levels)
  )

gene2_bp_range <- range(gene2_gwas_df$BP, na.rm = TRUE)
gene2_window_padding <- 5e4
gene2_window_start <- max(gene2_bp_range[1] - gene2_window_padding, 0)
gene2_window_end <- gene2_bp_range[2] + gene2_window_padding

gene2_position_breaks <- pretty(c(gene2_window_start, gene2_window_end), n = 4)
gene2_position_labels <- label_number(scale = 1e-6, accuracy = 0.01)

gene2_chr_values <- gene2_gwas_df %>%
  pull(Chr) %>%
  na.omit() %>%
  unique()

gene2_target_seqids <- c(paste0("chr", gene2_chr_values), paste0("Chr", gene2_chr_values), as.character(gene2_chr_values))

gene2_gff_genes <- gff %>%
  filter(
    seqid %in% gene2_target_seqids,
    type == "gene",
    start <= gene2_window_end,
    end >= gene2_window_start
  ) %>%
  mutate(
    gene_id = str_extract(attributes, "(?<=ID=)[^;]+"),
    strand_label = factor(if_else(strand == "+", "+", "-"), levels = c("+", "-")),
    forward = strand == "+"
  ) %>%
  arrange(start)

gene2_selected_snps <- gene2_fm_max$SNP
gene2_selected_idx <- match(gene2_selected_snps, map$SNP)
gene2_selected_snps <- gene2_selected_snps[!is.na(gene2_selected_idx)]
gene2_selected_idx <- gene2_selected_idx[!is.na(gene2_selected_idx)]

gene2_expr_gene <- ind %>%
  left_join(expr_counts %>% select(Taxa, Expr = all_of(gene2_id)), by = "Taxa")

gene2_geno_df <- as.data.frame(snp[, gene2_selected_idx, drop = FALSE])
colnames(gene2_geno_df) <- gene2_selected_snps

gene2_expr_gene <- bind_cols(gene2_expr_gene, gene2_geno_df)

gene2_allele_lookup <- tibble(
  SNP = gene2_selected_snps,
  REF = map$REF[gene2_selected_idx],
  ALT = map$ALT[gene2_selected_idx]
)

gene2_hap_df <- gene2_expr_gene %>%
  select(all_of(gene2_selected_snps)) %>%
  as.matrix()

gene2_hap_code <- apply(gene2_hap_df, 1, function(row) {
  if (any(is.na(row))) {
    return(NA_character_)
  }
  paste0(row, collapse = "")
})

gene2_hap_nuc <- apply(gene2_hap_df, 1, function(row) {
  if (any(is.na(row))) {
    return(NA_character_)
  }
  alleles <- vapply(seq_along(row), function(idx) {
    dosage <- row[idx]
    ref_allele <- gene2_allele_lookup$REF[idx]
    alt_allele <- gene2_allele_lookup$ALT[idx]
    if (is.na(dosage)) {
      return(NA_character_)
    }
    if (dosage == 0) {
      return(ref_allele)
    }
    if (dosage == 2) {
      return(alt_allele)
    }
    return(NA_character_)
  }, character(1))
  if (any(is.na(alleles))) {
    return(NA_character_)
  }
  paste0(alleles, collapse = "")
})

gene2_combo_label <- apply(gene2_hap_df, 1, function(row) {
  if (any(is.na(row))) {
    return(NA_character_)
  }
  alleles <- vapply(seq_along(row), function(idx) {
    dosage <- row[idx]
    ref_allele <- gene2_allele_lookup$REF[idx]
    alt_allele <- gene2_allele_lookup$ALT[idx]
    if (is.na(dosage)) {
      return(NA_character_)
    }
    if (dosage == 0) {
      return(ref_allele)
    }
    if (dosage == 2) {
      return(alt_allele)
    }
    if (dosage == 1) {
      return(paste0(ref_allele, "/", alt_allele))
    }
    NA_character_
  }, character(1))
  if (any(is.na(alleles))) {
    return(NA_character_)
  }
  paste0(alleles, collapse = "-")
})

gene2_expr_gene <- gene2_expr_gene %>%
  mutate(
    Hap_code = gene2_hap_code,
    Hap = gene2_hap_nuc,
    Allele_combo = gene2_combo_label
  )

gene2_allele_combo_counts <- gene2_expr_gene %>%
  filter(!is.na(Allele_combo)) %>%
  count(Allele_combo, sort = TRUE)

gene2_hap_counts <- gene2_expr_gene %>%
  filter(!is.na(Hap)) %>%
  count(Hap, sort = TRUE)

gene2_combos_keep <- gene2_allele_combo_counts %>%
  filter(n >= gene2_hap_threshold) %>%
  pull(Allele_combo)

gene2_label_lookup <- tibble(
  Original = gene2_combos_keep,
  Short = paste0("Hap", seq_along(gene2_combos_keep))
)

gene2_expr_violin_df <- gene2_expr_gene %>%
  filter(!is.na(Allele_combo), Allele_combo %in% gene2_combos_keep) %>%
  mutate(
    ViolinLabel = factor(
      gene2_label_lookup$Short[match(Allele_combo, gene2_label_lookup$Original)],
      levels = gene2_label_lookup$Short
    )
  )

gene2_violin_lookup <- gene2_expr_violin_df %>%
  distinct(Allele_combo, ViolinLabel)

gene2_allele_combo_counts <- gene2_allele_combo_counts %>%
  left_join(gene2_violin_lookup, by = "Allele_combo") %>%
  mutate(
    HistGroup = if_else(is.na(ViolinLabel), "Other", as.character(ViolinLabel)),
    HistGroup = factor(HistGroup, levels = c(as.character(gene2_label_lookup$Short), "Other")),
    HapNumber = if_else(is.na(ViolinLabel), NA_character_, str_remove(as.character(ViolinLabel), "^Hap"))
  )

gene2_hap_counts_45 <- gene2_allele_combo_counts %>%
  filter(Allele_combo %in% gene2_combos_keep) %>%
  transmute(
    Hap = Allele_combo,
    Display = gene2_label_lookup$Short[match(Allele_combo, gene2_label_lookup$Original)],
    n = n
  )

gene2_hap_violin <- gene2_expr_violin_df %>%
  ggplot(aes(x = ViolinLabel, y = Expr, fill = ViolinLabel)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, colour = "black") +
  scale_fill_manual(values = rep(hapPal, length.out = length(gene2_label_lookup$Short))) +
  labs(
    x = "Allele combination",
    y = paste(gene2_id, "expression")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + 
  ylim(5.5, 9)

gene2_mapping_caption <- paste(gene2_label_lookup$Short, gene2_label_lookup$Original, sep = " = ", collapse = "; ")

gene2_hap_bar <- gene2_hap_counts_45 %>%
  mutate(Display = factor(Display, levels = gene2_label_lookup$Short)) %>%
  ggplot(aes(x = Display, y = n)) +
  geom_col(fill = myPal[1], width = 0.8) +
  geom_text(aes(label = n), vjust = -0.3, size = 2.1) +
  scale_y_continuous("Individuals", expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Allele combination",
    subtitle = gene2_mapping_caption
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gene2_hist_palette <- c(
  setNames(rep(hapPal, length.out = length(gene2_label_lookup$Short)), gene2_label_lookup$Short),
  Other = myPal[2]
)

gene2_allele_hist_data <- gene2_allele_combo_counts %>%
  mutate(Allele_combo = factor(Allele_combo, levels = Allele_combo))

gene2_allele_hist_plot <- ggplot(gene2_allele_hist_data, aes(x = Allele_combo, y = n, fill = HistGroup)) +
  geom_col(width = 0.75, colour = NA) +
  geom_text(aes(label = n), vjust = -0.3, size = 1.8, colour = "black") +
  geom_text(
    data = gene2_allele_hist_data %>% filter(!is.na(HapNumber)),
    aes(y = n / 2, label = HapNumber, colour = HistGroup),
    fontface = "bold",
    size = 2.0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = gene2_hist_palette, guide = "none") +
  scale_colour_manual(values = gene2_hist_palette, guide = "none") +
  scale_y_continuous("Individuals", expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Allele combination across fine-mapped SNPs",
    y = "Individuals"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gene2_gene_segments <- gene2_gff_genes %>%
  mutate(
    colour_group = if_else(gene_id == gene2_id, "target", "other"),
    y = case_when(
      colour_group == "target" ~ 1.15,
      forward ~ 1,
      TRUE ~ 0
    ),
    x = if_else(forward, start, end),
    xend = if_else(forward, end, start)
  ) %>%
  select(x, xend, y, colour_group)

gene2_gene_bg <- gene2_gene_segments %>% filter(colour_group == "other")
gene2_gene_target <- gene2_gene_segments %>% filter(colour_group == "target")

gene2_gene_track <- ggplot() +
  geom_hline(yintercept = c(0, 1), colour = "#E0E0E0", linewidth = 0.25) +
  geom_segment(
    data = gene2_gene_bg,
    aes(x = x, xend = xend, y = y, yend = y),
    colour = "#1A1A1A",
    linewidth = 0.45,
    lineend = "round",
    arrow = arrow(length = grid::unit(1.0, "mm"), type = "open", angle = 18)
  ) +
  geom_segment(
    data = gene2_gene_target,
    aes(x = x, xend = xend, y = y, yend = y),
    colour = "#D62728",
    linewidth = 0.85,
    lineend = "round",
    arrow = arrow(length = grid::unit(1.8, "mm"), type = "open", angle = 18)
  ) +
  coord_cartesian(xlim = c(gene2_window_start, gene2_window_end), ylim = c(-0.4, 1.4)) +
  scale_x_continuous(
    name = "Chromosome 9",
    breaks = gene2_position_breaks,
    labels = gene2_position_labels
  ) +
  scale_y_continuous(
    breaks = c(1, 0),
    labels = c("+", "-"),
    expand = expansion(add = 0)
  ) +
  labs(y = NULL) +
  theme(
    axis.text.y = element_text(size = base_font_size - 2),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

gene2_highlight_palette <- c("Other" = "#BDBDBD")
gene2_highlight_palette[gene2_highlight_snps] <- rep(myPal[5:8], length.out = length(gene2_highlight_snps))

gene2_gwas_top <- gene2_gwas_df %>% filter(highlight != "Other")
gene2_gwas_bg <- gene2_gwas_df %>% filter(highlight == "Other")

gene2_fm_top <- gene2_fm_df %>% filter(highlight != "Other")
gene2_fm_bg <- gene2_fm_df %>% filter(highlight == "Other")

gene2_gwas_plot <- ggplot(gene2_gwas_df, aes(x = BP, y = logp)) +
  geom_point(
    data = gene2_gwas_bg,
    colour = gene2_highlight_palette["Other"],
    size = 1,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  geom_point(
    data = gene2_gwas_top,
    aes(colour = highlight),
    size = 2.4,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = gene2_highlight_palette,
    breaks = gene2_highlight_snps,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(breaks = gene2_position_breaks, labels = gene2_position_labels) +
  coord_cartesian(xlim = c(gene2_window_start, gene2_window_end)) +
  labs(
    x = NULL,
    y = expression(-log[10](p))
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

gene2_fm_plot <- ggplot(gene2_fm_df, aes(x = POS, y = variable_prob)) +
  geom_point(
    data = gene2_fm_bg,
    colour = gene2_highlight_palette["Other"],
    size = 1.3,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  geom_point(
    data = gene2_fm_top,
    aes(colour = highlight),
    size = 2.6,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = gene2_fm_top,
    aes(label = SNP),
    colour = "black",
    size = 2.1,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0.05,
    seed = 123
  ) +
  scale_colour_manual(
    values = gene2_highlight_palette,
    breaks = gene2_highlight_snps,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(breaks = gene2_position_breaks, labels = gene2_position_labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(gene2_window_start, gene2_window_end)) +
  labs(
    x = NULL,
    y = "PIP"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# gene2_final_plot <- wrap_plots(
#   A = gene2_gwas_plot,
#   B = gene2_hap_violin,
#   C = gene2_fm_plot,
#   D = gene2_gene_track,
#   design = "
# AB
# CB
# DB
# "
# ) +
#   plot_annotation(tag_levels = "A", tag_suffix = ")")

# =========================
# Combined outputs
# =========================

# combined_locus_plot <- wrap_plots(
#   A = gene1_gwas_plot,
#   B = gene2_gwas_plot,
#   C = gene1_fm_plot,
#   D = gene2_fm_plot,
#   E = gene1_gene_track,
#   F = gene2_gene_track,
#   G = gene1_hap_violin,
#   H = gene2_hap_violin,
#   design = "
# ABG
# CDH
# EF.
# "
# ) +
#   plot_layout(
#     heights = c(1, 1, 0.08),
#     widths = c(1, 1, 0.8)
#   ) +
#   plot_annotation(tag_levels = "A", tag_suffix = ")")

g2 <- wrap_plots(
  A = gene2_gwas_plot,
  C = gene2_fm_plot,
  D = gene2_gene_track,
  B = gene2_hap_violin,
  design = design
) +
  plot_layout(widths  = c(4, 1),
              heights = c(0.8, 0.8, 0.2))
g1 / g2 + plot_annotation(tag_levels = "A")


output_dir <- file.path("results", "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, paste0("Fig3", ".png")),
  width = 17.5,
  height = 12,
  dpi = 300,
  units = "cm"
)

# Also save a vector PDF copy
ggsave(
  filename = file.path(output_dir, paste0("Fig3", ".pdf")),
  width = 17.5,
  height = 12,
  units = "cm"
)

gene1_allele_hist_plot + gene2_allele_hist_plot + plot_annotation(tag_levels = "A")

ggsave(
  filename = file.path("results/figures/Supp/", paste0("FigS_Hap_counts", ".png")),
  width = 8.75,
  height = 6,
  dpi = 300,
  units = "cm"
)

# Also save a vector PDF copy
ggsave(
  filename = file.path("results/figures/Supp/", paste0("FigS_Hap_counts", ".pdf")),
  width = 8.75,
  height = 6,
  units = "cm"
)
