### Script: figures/Fig3.R
### Purpose: Generate Figure 3 panels: credible set counts per eGene, relationship between number of variants and credible sets, and distribution of aFC for eQTNs.
### Inputs:
###  - results/SusieR2/CS.csv, results/SusieR2/CS_VAR.csv
###  - results/SusieR2/PIP_dist.csv, results/SusieR2/PIP_Max.csv
###  - results/AF/Susie2_AF.txt
### Outputs: results/figures/Fig3.png
### Usage: Rscript figures/Fig3.R
### Dependencies: tidyverse, scales, ggpmisc, data.table, patchwork
### #####################
### Genetic diversity ###
#########################

library(tidyverse)
library(scales)
library(ggpmisc)
library(data.table)
library(patchwork)
library(ggbeeswarm)
library(ggsci)
library(ggrepel)

### ggplot2 settings (match Fig1 font sizing)
base_font_size <- 10
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")

cs <- read.csv("results/summary/SusieR/CS.csv")
cs.var <- read.csv("results/summary/SusieR/CS_VAR.csv")

pip.af <- fread("results/summary/AF/aFCn_MAF_se_filtered.tsv.gz")

h2r2 <- fread("results/summary/Heritability_and_r2/narrow_sense_h2_and_r2.tsv.gz")

h2r2 <- h2r2 %>%
  filter(gene %in% pip.af$gene_id & osca_r2_test > 0 & susie_r2_test > 0)

###
g1 <- ggplot(cs, aes(as.factor(cs), n)) + 
  geom_col(fill=myPal[1], colour="black") + 
  ylab("Number of eGenes") + 
  xlab("Number or credible causal sets") + 
  scale_y_continuous(labels = comma)

g2 <- ggplot(cs.var, aes(nVar, nCS)) + 
  geom_point(colour=myPal[2], size=2.5) + 
  geom_line() + 
  scale_x_log10(labels=comma) + 
  scale_y_continuous(labels = comma) + 
  annotation_logticks(sides = "b") + 
  xlab("Number of variants") +
  ylab("Number of credible sets")

g3 <- ggplot(pip.af, aes(x=log2_aFC)) +
  geom_histogram(position = "identity", alpha=.8, fill=myPal[3], colour="black", binwidth = 0.2) + 
  geom_vline(xintercept = c(-2, 2), linetype=2) + 
  ylab("Number of eQTNs") + 
  xlab(expression(italic(log)[2]~aFC))

g4 <- ggplot(h2r2, aes(as.factor(susie_n_variants), h2_pc_adjusted, fill=as.factor(susie_n_variants))) + 
  geom_boxplot() + 
  ylab("Narrow-sense heritability") +
  xlab("Number of credible sets") + 
  scale_fill_viridis_d(name="N causal sets") +
  theme(legend.position = "none")
  
g5 <- ggplot(h2r2, aes(osca_r2_test, susie_r2_test, colour=as.factor(susie_n_variants))) + 
  geom_point(alpha = 0.65) +
  geom_abline(slope = 1, linetype = 2, colour = "gray50") + 
  xlab(expression(Singe~variant~model~(italic(R^2)))) + 
  ylab(expression(Multiple~variant~model~(italic(R^2)))) + 
  scale_color_viridis_d(name="N causal sets") +
  guides(colour = guide_legend(ncol = 2)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 3)
  )

g6 <- ggplot(h2r2, aes(as.factor(susie_n_variants), susie_r2_test - osca_r2_test,  fill=as.factor(susie_n_variants))) + 
  geom_boxplot() + 
  ylab(expression(Delta~italic(R2))) +
  xlab("Number of credible sets") + 
  scale_fill_viridis_d(name="N causal sets") +
  theme(legend.position = "none")

fig2_summary <- g1 + g2 + g3 + g4 + g5 + g6

####################
### PSBS exmaple ###
####################
PSBS_GENE <- "Zm00001eb146510"
PSBS_TSS <- 179564407L
PSBS_CHR <- 3L
WINDOW_BP <- 2e6
window_start <- PSBS_TSS - WINDOW_BP
window_end <- PSBS_TSS + WINDOW_BP
window_mb <- c(window_start, window_end) / 1e6

psbs.gwas <- fread("results/full/osca/psbs")
psbs <- fread("results/eqtl_example/Zm00001eb146510_sample_expression_genotypes.tsv.gz")
psbs.fine_mapped <- fread("results/eqtl_example/Zm00001eb146510_fine_mapped_variants.tsv")

psbs.leads <- merge(
  pip.af[gene_id == PSBS_GENE],
  psbs.fine_mapped[, .(SNP, cs)],
  by.x = "variant_id",
  by.y = "SNP",
  all.x = TRUE,
  sort = FALSE
)
setorder(psbs.leads, cs)
if (nrow(psbs.leads) != 3L || anyNA(psbs.leads$cs)) {
  stop("Expected exactly three PSBS credible-set lead variants in pip.af.")
}
psbs.lead.gwas <- merge(
  psbs.gwas,
  psbs.leads[, .(SNP = variant_id, cs)],
  by = "SNP",
  all = FALSE,
  sort = FALSE
)
psbs.lead.gwas[, cs_label := factor(paste0("CS ", cs), levels = paste0("CS ", 1:3))]
lead_colours <- c("CS 1" = "#D55E00", "CS 2" = "#0072B2", "CS 3" = "#009E73")
manhattan_y_max <- ceiling(max(
  -log10(pmax(psbs.gwas$p, .Machine$double.xmin)),
  na.rm = TRUE
) + 2)

g_psbs_manhattan <- ggplot(
  psbs.gwas,
  aes(x = BP / 1e6, y = -log10(pmax(p, .Machine$double.xmin)))
) +
  geom_point(colour = "grey65", size = 1.2, alpha = 0.7) +
  geom_vline(xintercept = PSBS_TSS / 1e6, linetype = "dashed", colour = "grey30") +
  geom_point(
    data = psbs.lead.gwas,
    aes(colour = cs_label),
    size = 3.2
  ) +
  geom_text(
    data = psbs.lead.gwas,
    aes(label = SNP, colour = cs_label),
    size = 3,
    nudge_x = 0.05,
    nudge_y = 0.3,
    hjust = 0,
    vjust = 0,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = lead_colours, name = "Fine-mapped lead") +
  scale_y_continuous(
    limits = c(0, manhattan_y_max),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(xlim = window_mb, expand = FALSE) +
  labs(
    x = NULL,
    y = expression(-log[10](italic(P)))
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)
  )

# Stream only protein-coding gene features overlapping the PSBS +/- 2 Mb
# interval, then assign non-overlapping vertical lanes for the gene track.
gff_path <- "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
gff_cmd <- paste0(
  "gzip -dc ", shQuote(gff_path),
  " | awk -F'\\t' '$1 == \"chr", PSBS_CHR,
  "\" && $3 == \"gene\" && $5 >= ", window_start,
  " && $4 <= ", window_end, " { print }'"
)
gene_track <- fread(
  cmd = gff_cmd,
  sep = "\t",
  header = FALSE,
  col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
)
gene_track[, gene_id := sub(".*ID=([^;]+).*", "\\1", attributes)]
setorder(gene_track, start, end)
gene_track[, lane := NA_integer_]
lane_ends <- numeric()
for (i in seq_len(nrow(gene_track))) {
  available_lane <- which(lane_ends < gene_track$start[[i]])
  lane <- if (length(available_lane)) available_lane[[1L]] else length(lane_ends) + 1L
  gene_track$lane[[i]] <- lane
  lane_ends[[lane]] <- gene_track$end[[i]]
}
gene_track[, `:=`(
  start_mb = start / 1e6,
  end_mb = end / 1e6,
  midpoint_mb = (start + end) / 2e6,
  is_psbs = gene_id == PSBS_GENE
)]

g_psbs_genes <- ggplot() +
  geom_segment(
    data = gene_track[strand == "+" & is_psbs == FALSE],
    aes(x = start_mb, xend = end_mb, y = lane, yend = lane),
    colour = "grey45",
    arrow = grid::arrow(length = grid::unit(1.4, "mm"), ends = "last", type = "closed")
  ) +
  geom_segment(
    data = gene_track[strand == "-" & is_psbs == FALSE],
    aes(x = end_mb, xend = start_mb, y = lane, yend = lane),
    colour = "grey45",
    arrow = grid::arrow(length = grid::unit(1.4, "mm"), ends = "last", type = "closed")
  ) +
  geom_segment(
    data = gene_track[is_psbs == TRUE],
    aes(x = ifelse(strand == "+", start_mb, end_mb), xend = ifelse(strand == "+", end_mb, start_mb), y = lane, yend = lane),
    colour = "firebrick",
    linewidth = 1.1,
    arrow = grid::arrow(length = grid::unit(1.8, "mm"), ends = "last", type = "closed")
  ) +
  geom_text(
    data = gene_track[is_psbs == TRUE],
    aes(x = midpoint_mb, y = lane + 0.5, label = "PSBS"),
    colour = "firebrick",
    fontface = "bold",
    size = 2.8,
    vjust = 1.5
  ) +
  coord_cartesian(xlim = window_mb, expand = FALSE) +
  scale_y_reverse(
    limits = c(max(gene_track$lane) + 2.5, 0.5),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  labs(x = "Chromosome 3 position (Mb)", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 0, r = 5.5, b = 18, l = 5.5)
  )

lead_genotype_columns <- paste0("genotype_", psbs.leads$variant_id)
if (!all(lead_genotype_columns %in% names(psbs))) {
  stop("The PSBS example table does not contain all three lead-variant genotypes.")
}
psbs_haplotypes <- psbs[, c("expression_log2_cpm", lead_genotype_columns), with = FALSE]
psbs_haplotypes[, haplotype := do.call(paste, c(.SD, sep = "/")), .SDcols = lead_genotype_columns]
haplotype_counts <- psbs_haplotypes[, .N, by = haplotype]
retained_haplotypes <- haplotype_counts[N >= 30, haplotype]
psbs_haplotypes <- psbs_haplotypes[haplotype %chin% retained_haplotypes]
haplotype_order <- psbs_haplotypes[
  , .(median_expression = median(expression_log2_cpm, na.rm = TRUE)), by = haplotype
][order(median_expression), haplotype]
psbs_haplotypes[, haplotype := factor(haplotype, levels = haplotype_order)]

g_psbs_haplotypes <- ggplot(psbs_haplotypes, aes(x = haplotype, y = expression_log2_cpm, fill = haplotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.16, size = 1, alpha = 0.45) +
  scale_fill_manual(values = myPal) +
  labs(
    x = "Three-SNP genotype haplotype\n (CS 1 / CS 2 / CS 3)",
    y = expression(log[2]~CPM)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

fig2_psbs <- g_psbs_manhattan / g_psbs_genes +
  plot_layout(heights = c(5, 0.75))

fig2_bottom <- (fig2_psbs | g_psbs_haplotypes) +
  plot_layout(widths = c(2.3, 1))

fig2 <- (fig2_summary / fig2_bottom) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A")

# Save the summary and the PSBS example as separate figures.
ggsave(
  filename = "results/figures/Fig2.png",
  plot = fig2,
  width = 11.7,
  height = 10,
  dpi = 300
)
