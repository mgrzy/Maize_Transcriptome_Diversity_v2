# Figure 4 eQTL example: sourced by Fig4.R after its plotting setup.

example_gene <- "Zm00001eb146510"
groups <- c("Tropical", "NSS", "SS", "Iodent")
needed <- c("expression_file", "metadata_file", "eqtl_file", "geno_descriptor", "geno_map", "geno_individuals", "gff_file")
for (path in unlist(cfg[needed])) if (!file.exists(path)) stop_with("Required Figure 4 example input does not exist: %s", path)

fine_path <- c(cfg$fine_mapped_file_example, cfg$fine_mapped_pip_file)
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
  scale_fill_manual(values = setNames(myPal[1:4], groups), guide = "none") +
  labs(title = example_gene, x = NULL, y = expression(log[2]~CPM)) +
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
  )

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
track[, `:=`(start_mb = start / 1e6, end_mb = end / 1e6, is_target = gene_id == example_gene)]
gene_track_plot <- ggplot(track, aes(y = lane, yend = lane)) +
  geom_segment(aes(x = ifelse(strand == "+", start_mb, end_mb), xend = ifelse(strand == "+", end_mb, start_mb), colour = is_target),
    arrow = grid::arrow(length = grid::unit(1.2, "mm"), ends = "last", type = "closed")) +
  scale_colour_manual(values = c("FALSE" = "grey45", "TRUE" = "firebrick"), guide = "none") +
  coord_cartesian(xlim = xlim_mb, expand = FALSE) + scale_y_reverse(breaks = NULL) +
  labs(x = paste0("Chromosome ", gene_chr, " position (Mb)"), y = NULL) + theme_classic(base_size = 8) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

map <- fread(cfg$geno_map, data.table = FALSE); ind <- fread(cfg$geno_individuals, header = FALSE, data.table = FALSE); names(ind)[1] <- "Line"
idx <- match(fine$SNP, map$SNP); if (anyNA(idx)) stop_with("Fine-mapped variants are absent from Div639.geno.map.")
geno <- attach.big.matrix(cfg$geno_descriptor)
dosage <- as.data.frame(geno[, idx, drop = FALSE]); names(dosage) <- fine$SNP
geno_samples <- bind_cols(ind, dosage) %>% inner_join(sample_data, by = "Line")
allele_freq <- geno_samples %>% pivot_longer(all_of(fine$SNP), names_to = "SNP", values_to = "dosage") %>%
  group_by(Group, SNP) %>% summarise(alt_frequency = mean(as.numeric(dosage) / 2, na.rm = TRUE), .groups = "drop") %>%
  mutate(Group = factor(Group, levels = groups), SNP = factor(SNP, levels = fine$SNP))
major_minor_frequency <- allele_freq %>%
  transmute(Group, SNP, `Alternative allele` = alt_frequency, `Major allele` = 1 - alt_frequency) %>%
  pivot_longer(c(`Major allele`, `Alternative allele`), names_to = "allele", values_to = "frequency") %>%
  mutate(allele = factor(allele, levels = c("Major allele", "Alternative allele")))
stacked_allele_plot <- ggplot(major_minor_frequency, aes(Group, frequency, fill = allele)) +
  geom_col(width = .72) + facet_wrap(~SNP, nrow = 1) +
  scale_fill_manual(values = c("Major allele" = "grey75", "Alternative allele" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1)) + labs(x = NULL, y = "Allele frequency", fill = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6), legend.position = "bottom", legend.text = element_text(size = 6))

hap_data <- geno_samples %>% mutate(haplotype = apply(select(., all_of(fine$SNP)), 1, paste, collapse = "/"))
top_haps <- hap_data %>% count(haplotype, sort = TRUE) %>% slice_head(n = 4) %>% pull(haplotype)
hap_data <- hap_data %>% filter(haplotype %in% top_haps) %>%
  mutate(haplotype = factor(haplotype, levels = top_haps))
haplotype_effect_plot <- ggplot(hap_data, aes(haplotype, expression_log2_cpm, fill = haplotype)) +
  geom_boxplot(width = .7, outlier.size = .4) +
  scale_fill_manual(values = myPal, guide = "none") +
  labs(x = "Top four genotype haplotypes", y = expression(log[2]~CPM)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

example_manhattan <- manhattan_plot / gene_track_plot + plot_layout(heights = c(5, 1))
example_panels <- group_effect_plot | example_manhattan | stacked_allele_plot | haplotype_effect_plot +
  plot_layout(widths = c(0.9, 2.25, 1.15, 1.15))
