#!/usr/bin/env Rscript
# Figure 4: Fst across DEG-significance deciles for fine-mapped variants linked
# to significantly differentially expressed genes. Each facet is one DEG/Fst
# comparison and uses the corresponding -log10(limma adjusted P value).

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(data.table)
  library(bigmemory)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  fine_mapped_file = "results/summary/AF/aFCn_MAF_se_filtered_DEG_Fst.tsv.gz",
  deg_file = "results/summary/DEG/DEG_limma_summary.csv",
  gerp_file = "results/summary/GERP/positive_GERP_mean_by_gene_and_promoter.csv",
  credible_set_file = "results/summary/SusieR/CS_Gene.csv",
  gerp_test_file = "results/summary/GERP/GERP_effect_size_and_credible_set_tests.csv",
  expression_file = "data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz",
  metadata_file = "data/raw/meta/WiDiv_meta.csv",
  eqtl_file = "results/summary/eqtl_example/Zm00001eb392880",
  fine_mapped_file_example = "results/summary/AF/aFCn_MAF_se_filtered_DEG_Fst.tsv.gz",
  fine_mapped_pip_file = "results/full/SusieR/PIP_Zm00001eb392880.tsv.gz",
  geno_descriptor = "data/processed/geno/rMVP/Div639.geno.desc",
  geno_map = "data/processed/geno/rMVP/Div639.geno.map",
  geno_individuals = "data/processed/geno/rMVP/Div639.geno.ind",
  gff_file = "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz",
  output_file = "results/figures/Fig4.png"
)
for (path in unlist(cfg[1:5])) {
  if (!file.exists(path)) stop_with("Required input does not exist: %s", path)
}
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)

comparison_labels <- c(
  "SS_vs_Tropical" = "SS vs Tropical",
  "NSS_vs_Tropical" = "NSS vs Tropical",
  "Iodent_vs_Tropical" = "Iodent vs Tropical",
  "Public_SS_vs_NSS" = "Public: SS vs NSS",
  "exPVP_SS_vs_NSS" = "exPVP: SS vs NSS"
)
myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")
comparison_colours <- setNames(myPal[seq_along(comparison_labels)], unname(comparison_labels))
gerp_colours <- c("Other genes" = myPal[1], "Top 10% GERP genes" = myPal[2])

fm <- read_tsv(cfg$fine_mapped_file, show_col_types = FALSE)
deg <- read_csv(cfg$deg_file, show_col_types = FALSE) %>%
  rename(gene_id = gene)
gerp <- read_csv(cfg$gerp_file, show_col_types = FALSE)
credible_sets <- read_csv(cfg$credible_set_file, show_col_types = FALSE) %>%
  rename(gene_id = gene, n_credible_sets = cs)
gerp_tests <- read_csv(cfg$gerp_test_file, show_col_types = FALSE)

required_fm_columns <- c("gene_id", "variant_id", paste0(names(comparison_labels), "_Fst"))
required_deg_columns <- c(
  "gene_id",
  paste0(names(comparison_labels), "_DEG"),
  paste0(names(comparison_labels), "_adj.P.Val")
)
missing_fm <- setdiff(required_fm_columns, names(fm))
missing_deg <- setdiff(required_deg_columns, names(deg))
if (length(missing_fm)) {
  stop_with("Fine-mapped DEG/Fst table is missing: %s", paste(missing_fm, collapse = ", "))
}
if (length(missing_deg)) {
  stop_with("DEG summary is missing: %s", paste(missing_deg, collapse = ", "))
}
required_gerp_columns <- c("gene_id", "gene_mean_positive_GERP")
required_cs_columns <- c("gene_id", "n_credible_sets")
if (length(missing_gerp <- setdiff(required_gerp_columns, names(gerp)))) {
  stop_with("GERP summary is missing: %s", paste(missing_gerp, collapse = ", "))
}
if (length(missing_cs <- setdiff(required_cs_columns, names(credible_sets)))) {
  stop_with("Credible-set summary is missing: %s", paste(missing_cs, collapse = ", "))
}
if (!all(c("test", "p_value") %in% names(gerp_tests))) {
  stop_with("GERP test summary must contain test and p_value columns.")
}
effect_p_value <- gerp_tests %>%
  filter(test == "Wilcoxon rank-sum") %>%
  pull(p_value)
cs_p_value <- gerp_tests %>%
  filter(test == "Quasi-Poisson GLM") %>%
  pull(p_value)
if (length(effect_p_value) != 1L || length(cs_p_value) != 1L ||
    !is.finite(effect_p_value) || !is.finite(cs_p_value)) {
  stop_with("GERP test summary must contain one finite P value for each plotted test.")
}
format_p_value <- function(p_value) formatC(p_value, format = "e", digits = 2)

# Top 90th percentile means the highest-scoring 10% of genes. GERP is
# classified at the gene level before joining to variant-level aFC estimates.
gerp_threshold <- quantile(gerp$gene_mean_positive_GERP, probs = 0.90, na.rm = TRUE)
gerp_groups <- gerp %>%
  filter(is.finite(gene_mean_positive_GERP)) %>%
  transmute(
    gene_id,
    gene_mean_positive_GERP,
    gerp_group = if_else(
      gene_mean_positive_GERP >= gerp_threshold,
      "Top 10% GERP genes",
      "Other genes"
    )
  ) %>%
  mutate(gerp_group = factor(gerp_group, levels = names(gerp_colours)))
if (anyDuplicated(gerp_groups$gene_id)) stop_with("GERP summary has duplicate gene IDs.")

# Each aFC row is a fine-mapped gene-variant pair. First calculate deciles once
# per unique DEG gene using -log10(adjusted P), separately for each comparison.
# Those gene-level deciles are then assigned to all fine-mapped variants linked
# to that gene, preventing genes with many fine-mapped variants from changing
# the decile cut points.
plot_data <- imap_dfr(comparison_labels, function(label, comparison) {
  deg_column <- paste0(comparison, "_DEG")
  fst_column <- paste0(comparison, "_Fst")
  adjusted_p_column <- paste0(comparison, "_adj.P.Val")

  gene_significance <- deg %>%
    transmute(
      gene_id,
      is_DEG = .data[[deg_column]],
      adjusted_p = .data[[adjusted_p_column]]
    ) %>%
    filter(is_DEG %in% TRUE, is.finite(adjusted_p), adjusted_p >= 0) %>%
    mutate(
      minus_log10_adjusted_p = -log10(pmax(adjusted_p, .Machine$double.xmin)),
      significance_decile = ntile(minus_log10_adjusted_p, 10L)
    )

  fm %>%
    filter(is.finite(.data[[fst_column]])) %>%
    transmute(
      gene_id,
      variant_id,
      Fst = .data[[fst_column]]
    ) %>%
    inner_join(gene_significance, by = "gene_id") %>%
    mutate(
      comparison = label
    )
})

if (!nrow(plot_data)) stop_with("No fine-mapped variants linked to significant DEGs had finite Fst values.")
if (any(table(plot_data$comparison) == 0L)) stop_with("At least one comparison has no plottable values.")

plot_data <- plot_data %>%
  mutate(
    comparison = factor(comparison, levels = unname(comparison_labels)),
    significance_decile = factor(
      significance_decile,
      levels = 1:10,
      labels = paste0("D", 1:10)
    )
  )

panel <- ggplot(plot_data, aes(x = significance_decile, y = Fst, fill = comparison)) +
  geom_boxplot(
    width = 0.7,
    colour = "#1D4F5A",
    linewidth = 0.25,
    outlier.size = 0.35,
    outlier.alpha = 0.65
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.1,
    stroke = 0.25,
    fill = "white",
    colour = "#1D4F5A"
  ) +
  scale_fill_manual(values = comparison_colours, guide = "none") +
  facet_wrap(~ comparison, nrow = 1) +
  labs(x = expression(-log[10]("adjusted P")~"decile (low to high)"), y = expression(F[ST])) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 9),
    strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.25),
    strip.text = element_text(size = 7),
    panel.spacing = unit(0.65, "lines")
  )

maf_plot <- fm %>%
  mutate(MAF_bin = cut(MAF, breaks = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))) %>%
  filter(!is.na(MAF_bin), is.finite(log2_aFC)) %>%
  ggplot(aes(x = MAF_bin, y = abs(log2_aFC), fill = MAF_bin)) +
  geom_boxplot(width = 0.7, linewidth = 0.25, outlier.size = 0.35, outlier.alpha = 0.65) +
  scale_fill_manual(values = myPal, guide = "none") +
  labs(x = "MAF", y = expression(abs(log[2](aFC)))) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.title = element_text(size = 9)
  )

gerp_effect_data <- fm %>%
  select(gene_id, variant_id, log2_aFC) %>%
  inner_join(gerp_groups, by = "gene_id") %>%
  filter(is.finite(log2_aFC))
if (!nrow(gerp_effect_data)) stop_with("No fine-mapped variants could be assigned to a GERP group.")

gerp_effect_plot <- ggplot(
  gerp_effect_data,
  aes(x = gerp_group, y = abs(log2_aFC), fill = gerp_group)
) +
  geom_boxplot(width = 0.7, linewidth = 0.25, outlier.size = 0.35, outlier.alpha = 0.65) +
  scale_fill_manual(values = gerp_colours, guide = "none") +
  annotate(
    "text", x = 1.5, y = Inf,
    label = paste0("Wilcoxon P = ", format_p_value(effect_p_value)),
    vjust = 1.5, size = 2.6
  ) +
  labs(x = NULL, y = expression(abs(log[2](aFC)))) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.margin = margin(t = 8, r = 5.5, b = 5.5, l = 5.5)
  )

credible_set_data <- credible_sets %>%
  inner_join(gerp_groups %>% select(gene_id, gerp_group), by = "gene_id") %>%
  filter(is.finite(n_credible_sets), n_credible_sets >= 0)
if (!nrow(credible_set_data)) stop_with("No credible-set counts could be assigned to a GERP group.")
credible_set_distribution <- credible_set_data %>%
  count(gerp_group, n_credible_sets, name = "n_genes") %>%
  group_by(gerp_group) %>%
  mutate(proportion_genes = n_genes / sum(n_genes)) %>%
  ungroup()

credible_set_plot <- ggplot(
  credible_set_distribution,
  aes(x = factor(n_credible_sets), y = proportion_genes, fill = gerp_group)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = gerp_colours, name = NULL) +
  annotate(
    "text", x = 5.5, y = Inf,
    label = paste0("Quasi-Poisson GLM P = ", format_p_value(cs_p_value)),
    vjust = 1.5, size = 2.6
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Number of credible sets", y = "Proportion of eGenes") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 9),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 7),
    plot.margin = margin(t = 8, r = 5.5, b = 5.5, l = 5.5)
  )

top_row <- maf_plot + gerp_effect_plot + credible_set_plot +
  plot_layout(widths = c(1, 1, 1.5))
####

# Figure 4 eQTL example: sourced by Fig4.R after its plotting setup.

example_gene <- "Zm00001eb392880"
groups <- c("Tropical", "NSS", "SS", "Iodent")
needed <- c("expression_file", "metadata_file", "eqtl_file", "geno_descriptor", "geno_map", "geno_individuals", "gff_file")
for (path in unlist(cfg[needed])) if (!file.exists(path)) stop_with("Required Figure 4 example input does not exist: %s", path)

fine_path <- c(cfg$fine_mapped_pip_file, cfg$fine_mapped_file_example)
fine_path <- fine_path[file.exists(fine_path)][1L]
if (!length(fine_path) || is.na(fine_path)) {
  stop_with(paste0("No fine-mapped-variant input found for %s. Expected %s or %s."),
            example_gene, cfg$fine_mapped_file_example, cfg$fine_mapped_pip_file, example_gene)
}

expr <- fread(cfg$expression_file, select = c("Line", example_gene), data.table = FALSE, check.names = FALSE)
names(expr) <- c("Line", "expression_log2_cpm")
meta <- fread(cfg$metadata_file, data.table = FALSE, check.names = FALSE)
if (!all(c("Taxa", "Group") %in% names(meta))) stop_with("WiDiv metadata must contain Taxa and Group.")
sample_data <- expr %>% mutate(Line = as.character(Line)) %>%
  inner_join(
    meta %>% transmute(
      Line = as.character(Taxa),
      Group = recode(as.character(Group), "IDT" = "Iodent")
    ),
    by = "Line"
  ) %>%
  filter(Group %in% groups, is.finite(expression_log2_cpm)) %>%
  mutate(Group = factor(Group, levels = groups))
if (!nrow(sample_data)) stop_with("No %s samples occur in the four requested groups.", example_gene)

group_effect_plot <- ggplot(sample_data, aes(Group, expression_log2_cpm, fill = Group)) +
  geom_boxplot(width = .7, outlier.size = .4) +
  geom_jitter(width = .14, size = .55, alpha = .4, colour = "grey20") +
  scale_fill_manual(values = setNames(myPal[1:4], groups), guide = "none") +
  labs(x = NULL, y = expression(log[2]~CPM)) +
  theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "italic", size = 9))

eqtl <- fread(cfg$eqtl_file, data.table = FALSE)
if (!all(c("SNP", "Chr", "BP", "p", "Probe_Chr", "Probe_bp") %in% names(eqtl))) stop_with("Unexpected eQTL columns for %s.", example_gene)
gene_chr <- as.character(eqtl$Probe_Chr[1]); tss <- as.numeric(eqtl$Probe_bp[1]); window <- 2e6
fine <- fread(fine_path, data.table = FALSE)
if ("gene_id" %in% names(fine)) fine <- fine[as.character(fine$gene_id) == example_gene, , drop = FALSE]
if ("gene" %in% names(fine)) fine <- fine[as.character(fine$gene) == example_gene, , drop = FALSE]
snp_col <- intersect(c("SNP", "variant_id", "variant"), names(fine))[1]
cs_col <- intersect(c("cs", "CS", "credible_set"), names(fine))[1]
if (is.na(snp_col)) stop_with("Fine-mapping file has no SNP or variant_id column: %s", fine_path)
cs_example_file <- file.path("results/eqtl_example", paste0(example_gene, "_fine_mapped_variants.tsv"))
cs_lookup <- if (file.exists(cs_example_file)) {
  fread(cs_example_file, data.table = FALSE) %>% transmute(SNP = as.character(SNP), cs_example = as.character(cs))
} else tibble(SNP = character(), cs_example = character())
fine <- tibble(SNP = as.character(fine[[snp_col]]), cs = if (is.na(cs_col)) "Fine-mapped" else as.character(fine[[cs_col]])) %>%
  left_join(cs_lookup, by = "SNP") %>%
  mutate(cs = if_else(!is.na(cs_example), cs_example, cs)) %>% select(-cs_example) %>%
  distinct(SNP, .keep_all = TRUE) %>% inner_join(eqtl %>% select(SNP, Chr, BP, p), by = "SNP") %>%
  mutate(cs = factor(paste0("CS ", cs)))
if (!nrow(fine)) stop_with("Fine-mapped variants for %s are not present in its association data.", example_gene)
xlim_mb <- (tss + c(-window, window)) / 1e6

manhattan_plot <- ggplot(eqtl, aes(BP / 1e6, -log10(pmax(p, .Machine$double.xmin)))) +
  geom_point(colour = "grey65", size = .55, alpha = .65) +
  geom_vline(xintercept = tss / 1e6, linetype = "dashed") +
  geom_point(
    data = fine,
    aes(x = BP / 1e6, y = -log10(pmax(p, .Machine$double.xmin)), colour = cs),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Fine-mapped") + coord_cartesian(xlim = xlim_mb, expand = FALSE) +
  labs(x = NULL, y = expression(-log[10](italic(P)))) + theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "inside", legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1), legend.text = element_text(size = 6), legend.title = element_text(size = 7)
  ) + 
  ylim(0, 25)

gff_cmd <- paste0("gzip -dc ", shQuote(cfg$gff_file), " | awk -F'\\t' '$1 == \"chr", gene_chr,
                  "\" && $3 == \"gene\" && $5 >= ", tss - window, " && $4 <= ", tss + window, " {print}'")
track <- fread(cmd = gff_cmd, sep = "\t", header = FALSE, col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
if (!nrow(track)) stop_with("No gene features were found in the %s eQTL window.", example_gene)
track[, gene_id := sub(".*ID=([^;]+).*", "\\1", attributes)]; setorder(track, start, end); track[, lane := 1L]
lane_ends <- numeric()
for (i in seq_len(nrow(track))) {
  available <- which(lane_ends < track$start[i])
  lane <- if (length(available)) available[1L] else length(lane_ends) + 1L
  track$lane[i] <- lane
  lane_ends[lane] <- track$end[i]
}
track[, `:=`(
  start_mb = start / 1e6,
  end_mb = end / 1e6,
  # Inset arrowheads from feature edges so they are not clipped at the plot boundary.
  arrow_inset_mb = pmin((end - start) / 1e6 * 0.15, 0.015),
  is_target = gene_id == example_gene
)]
track[, `:=`(
  arrow_x = ifelse(strand == "+", start_mb + arrow_inset_mb, end_mb - arrow_inset_mb),
  arrow_dx = ifelse(strand == "+", end_mb - start_mb - 2 * arrow_inset_mb, start_mb - end_mb + 2 * arrow_inset_mb),
  gene_y = ifelse(strand == "+", 1.1, 1.0)
)]

gene_track_plot <- ggplot(track) +
  metR::geom_arrow(
    aes(x = arrow_x, y = gene_y, dx = arrow_dx, dy = 0, colour = is_target),
    arrow.length = 0.35,
    arrow.ends = "last",
    arrow.type = "closed",
    linewidth = 0.45
  ) +
  scale_colour_manual(values = c("FALSE" = "grey45", "TRUE" = "firebrick"), guide = "none") +
  coord_cartesian(xlim = xlim_mb, ylim = c(0.9, 1.2), expand = FALSE) + scale_y_continuous(breaks = NULL) +
  labs(x = paste0("Chromosome ", gene_chr, " position (Mb)"), y = NULL) + theme_classic(base_size = 8) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

map <- fread(cfg$geno_map, data.table = FALSE); ind <- fread(cfg$geno_individuals, header = FALSE, data.table = FALSE); names(ind)[1] <- "Line"
idx <- match(fine$SNP, map$SNP); if (anyNA(idx)) stop_with("Fine-mapped variants are absent from Div639.geno.map.")
geno <- attach.big.matrix(cfg$geno_descriptor)
dosage <- as.data.frame(geno[, idx, drop = FALSE]); names(dosage) <- fine$SNP
geno_samples <- bind_cols(ind, dosage) %>% inner_join(sample_data, by = "Line")

hap_data <- geno_samples %>%
  filter(if_all(all_of(fine$SNP), ~ !is.na(.))) %>%
  mutate(haplotype = apply(select(., all_of(fine$SNP)), 1, paste, collapse = "/"))
top_haps <- hap_data %>% count(haplotype, sort = TRUE) %>% slice_head(n = 4) %>% pull(haplotype)
hap_order <- hap_data %>%
  filter(haplotype %in% top_haps) %>%
  group_by(haplotype) %>%
  summarise(median_expression = median(expression_log2_cpm, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_expression) %>%
  pull(haplotype)
hap_data <- hap_data %>% filter(haplotype %in% hap_order) %>%
  mutate(haplotype = factor(haplotype, levels = hap_order))

haplotype_frequency_data <- geno_samples %>%
  filter(if_all(all_of(fine$SNP), ~ !is.na(.))) %>%
  mutate(
    haplotype = apply(select(., all_of(fine$SNP)), 1, paste, collapse = "/"),
    haplotype = if_else(haplotype %in% hap_order, haplotype, "Other"),
    haplotype = factor(haplotype, levels = c(hap_order, "Other"))
  ) %>%
  count(Group, haplotype, name = "n_lines") %>%
  group_by(Group) %>%
  mutate(frequency = n_lines / sum(n_lines)) %>%
  ungroup() %>%
  mutate(Group = factor(Group, levels = groups))

haplotype_colours <- setNames(myPal[seq_along(c(hap_order, "Other"))], c(hap_order, "Other"))
haplotype_frequency_plot <- ggplot(haplotype_frequency_data, aes(Group, frequency, fill = haplotype)) +
  geom_col(width = .72, colour = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = haplotype_colours, name = "Three-SNP\nhaplotype") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Haplotype frequency") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
    legend.position = "none"
  )

haplotype_effect_plot <- ggplot(hap_data, aes(haplotype, expression_log2_cpm, fill = haplotype)) +
  geom_boxplot(width = .7, outlier.size = .4) +
  geom_jitter(width = .14, size = .55, alpha = .4, colour = "grey20") +
  scale_fill_manual(values = myPal, guide = "none") +
  labs(x = "Top four genotype haplotypes", y = expression(log[2]~CPM)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

example_manhattan <- manhattan_plot / gene_track_plot + plot_layout(heights = c(5, 1))
example_panels <- group_effect_plot + example_manhattan + haplotype_frequency_plot + haplotype_effect_plot +
  plot_layout(widths = c(1, 2, 1, 1))
###

figure <- (top_row / panel / example_panels) +
  plot_layout(heights = c(0.9, 1, 1.25)) +
  plot_annotation(tag_levels = "A")
ggsave(cfg$output_file, figure, width = 320, height = 220, units = "mm", dpi = 100, device = ragg::agg_png)
message("Wrote Figure 4 with ", nrow(plot_data), " fine-mapped DEG-variant rows and the eQTL example to ", cfg$output_file)
