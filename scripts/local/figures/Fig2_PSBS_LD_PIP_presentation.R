#!/usr/bin/env Rscript
# Rscript --vanilla scripts/local/figures/Fig2_PSBS_LD_PIP_presentation.R
# Add --keep-haplotypes to retain panel I between the locus and PIP panels.
# LD is raw dosage Pearson r^2 against the CS1 lead, using samples with
# complete PSBS expression/covariates. PIP values are the saved SuSiE
# credible-set members; absent variants are not assigned zero PIP.
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(bigmemory)
  library(patchwork)
  library(ggrepel)
})
set.seed(42)
in_example <- FALSE
for (source_expr in parse("scripts/local/figures/Fig2.R")) {
  object <- if (is.call(source_expr) && identical(source_expr[[1L]], as.name("<-")) &&
                is.symbol(source_expr[[2L]])) as.character(source_expr[[2L]]) else ""
  if (object == "PSBS_GENE") in_example <- TRUE
  operation <- if (is.call(source_expr)) as.character(source_expr[[1L]]) else ""
  if (in_example || object %in% c("base_font_size", "myPal", "pip.af") ||
      operation %in% c("theme_set", "theme_update")) eval(source_expr)
  if (object == "fig2_bottom") break
}

# Display +/- 1 Mb around the TSS while retaining the original analysis inputs.
display_window_mb <- (PSBS_TSS + c(-1e6, 1e6)) / 1e6

reference <- psbs.leads[cs == 1, variant_id]
if (length(reference) != 1L) stop("Expected exactly one CS1 lead.")
top_hit <- psbs.gwas[which.min(p), SNP]
if (psbs.leads[cs == 1, variant_id] != top_hit) {
  stop("CS1 lead is not the top association: ", top_hit)
}

# Stream the local genotype-map rows and their original matrix columns.
map_path <- "data/processed/geno/rMVP/Div639.geno.map"
map_command <- paste0(
  "awk 'BEGIN { OFS = \"\\t\" } FNR > 1 && $2 == ", PSBS_CHR,
  " && $3 >= ", window_start, " && $3 <= ", window_end,
  " { print FNR - 1, $1, $3 }' ", shQuote(map_path)
)
local_map <- fread(cmd = map_command,
                   col.names = c("matrix_column", "SNP", "BP"))
if (anyDuplicated(local_map$SNP)) stop("Duplicate local genotype SNP IDs.")
map_index <- match(psbs.gwas$SNP, local_map$SNP)
if (anyNA(map_index)) stop("Some association SNPs are absent from the genotype map.")
local_map <- local_map[map_index]

individuals <- fread("data/processed/geno/rMVP/Div639.geno.ind",
                     header = FALSE, col.names = "Line", colClasses = "character")
covariates <- fread("data/processed/osca/WiDiv631_PC.txt")
cov_index <- match(individuals$Line, as.character(covariates$Taxa))
expression_index <- match(individuals$Line, as.character(psbs$Line))
if (anyNA(cov_index) || anyNA(expression_index)) stop("Missing PSBS sample identifiers.")
cov_columns <- setdiff(names(covariates)[vapply(covariates, is.numeric, logical(1))], "Taxa")
keep_samples <- complete.cases(covariates[cov_index, ..cov_columns]) &
  is.finite(psbs$expression_log2_cpm[expression_index])
geno <- attach.big.matrix("data/processed/geno/rMVP/Div639.geno.desc")
if (nrow(geno) != nrow(individuals)) stop("Genotype sample count mismatch.")
dosages <- geno[which(keep_samples), local_map$matrix_column, drop = FALSE]
storage.mode(dosages) <- "double"
if (anyNA(dosages)) stop("Unexpected missing genotype dosages.")
reference_index <- match(reference, local_map$SNP)
if (is.na(reference_index) || var(dosages[, reference_index]) == 0) {
  stop("Reference SNP is absent or monomorphic.")
}
r_squared <- as.numeric(cor(dosages, dosages[, reference_index]))^2
r_squared[!is.finite(r_squared)] <- NA_real_
ld <- data.table(SNP = local_map$SNP, LD_r2 = pmin(1, pmax(0, r_squared)))
psbs.gwas[ld, LD_r2 := i.LD_r2, on = "SNP"]
psbs.gwas[, minus_log10_p := -log10(pmax(p, .Machine$double.xmin))]
setorder(psbs.gwas, LD_r2, na.last = FALSE)
lead_points <- merge(psbs.gwas, psbs.leads[, .(SNP = variant_id, cs)], by = "SNP")
setorder(lead_points, cs)
lead_points[, label := paste0("CS ", cs, "\n", SNP)]
lead_points[, annotation_colour := unname(lead_colours[paste0("CS ", cs)])]
fine_points <- psbs.gwas[SNP %in% psbs.fine_mapped$SNP]

ld_scale <- function() scale_colour_viridis_c(
  name = expression(LD~r^2~"with CS1"), limits = c(0, 1),
  na.value = "grey75", breaks = c(0, 0.5, 1),
  guide = guide_colourbar(direction = "horizontal",
                         barheight = grid::unit(3, "mm"),
                         barwidth = grid::unit(24, "mm"),
                         theme = theme(legend.title.position = "top"))
)
locus_theme <- theme_classic(base_size = 10) +
  theme(axis.title = element_text(size = rel(1.15)),
        plot.margin = margin(3, 3, 3, 3))

manhattan <- ggplot(psbs.gwas, aes(BP / 1e6, minus_log10_p)) +
  geom_vline(xintercept = PSBS_TSS / 1e6, linetype = "dashed", colour = "grey45") +
  geom_point(aes(colour = LD_r2), size = 1, alpha = 0.8) +
  geom_point(data = fine_points, aes(fill = LD_r2), shape = 21,
             colour = "black", size = 2.5, stroke = 0.4) +
  geom_point(data = lead_points, aes(fill = LD_r2), shape = 21,
             colour = "black", size = 3.2, stroke = 0.6) +
  geom_text_repel(data = lead_points, aes(label = label),
                  nudge_x = c(-0.45, 0.4, -0.78), nudge_y = c(2, 3, 4),
                  size = 3, seed = 42, min.segment.length = 0,
                  colour = lead_points$annotation_colour,
                  segment.colour = lead_points$annotation_colour, segment.size = 0.35,
                  arrow = grid::arrow(length = grid::unit(1.5, "mm")),
                  box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  ld_scale() +
  scale_fill_viridis_c(limits = c(0, 1), na.value = "grey75", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  coord_cartesian(xlim = display_window_mb, ylim = c(0, manhattan_y_max + 3), expand = FALSE) +
  labs(x = NULL, y = expression(-log[10](italic(P)))) + locus_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.background = element_rect(fill = "white"))

pip_points <- merge(psbs.fine_mapped, ld, by = "SNP", all.x = TRUE)
if (anyNA(pip_points$LD_r2) || any(!is.finite(pip_points$variable_prob)) ||
    any(pip_points$variable_prob < 0 | pip_points$variable_prob > 1)) {
  stop("Invalid saved PIP or missing LD for fine-mapped variants.")
}
pip_leads <- pip_points[SNP %in% lead_points$SNP]
setorder(pip_leads, cs)
pip_leads[, label := paste0("CS ", cs, "\n", SNP)]
pip_leads[, annotation_colour := unname(lead_colours[paste0("CS ", cs)])]
pip_plot <- ggplot(pip_points, aes(POS / 1e6, variable_prob)) +
  geom_vline(xintercept = PSBS_TSS / 1e6, linetype = "dashed", colour = "grey45") +
  geom_segment(aes(xend = POS / 1e6, yend = 0, colour = LD_r2), linewidth = 0.4) +
  geom_point(aes(fill = LD_r2), shape = 21, colour = "black", size = 2.5) +
  geom_text_repel(data = pip_leads, aes(label = label),
                  nudge_x = c(-0.55, 0.45, -0.45), nudge_y = c(0.08, 0.05, 0.27),
                  size = 3, seed = 42, min.segment.length = 0,
                  colour = pip_leads$annotation_colour,
                  segment.colour = pip_leads$annotation_colour, segment.size = 0.35,
                  arrow = grid::arrow(length = grid::unit(1.5, "mm")),
                  box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  scale_colour_viridis_c(limits = c(0, 1), guide = "none") +
  scale_fill_viridis_c(limits = c(0, 1), guide = "none") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(xlim = display_window_mb, ylim = c(0, 1.18), expand = FALSE) +
  labs(x = NULL, y = "Posterior inclusion probabilities") + locus_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

track_data <- copy(gene_track)
track_data <- track_data[
  strand %in% c("+", "-") & end_mb >= display_window_mb[1] &
    start_mb <= display_window_mb[2]
]
track_data[, `:=`(
  strand_row = ifelse(strand == "+", 2, 1),
  arrow_start = ifelse(strand == "+", start_mb, end_mb),
  arrow_end = ifelse(strand == "+", end_mb, start_mb)
)]
gene_track <- ggplot(track_data) +
  geom_hline(yintercept = c(1, 2), colour = "grey90", linewidth = 0.25) +
  geom_segment(
    data = track_data[is_psbs == FALSE],
    aes(x = arrow_start, xend = arrow_end, y = strand_row, yend = strand_row),
    colour = "grey45", linewidth = 0.5,
    arrow = grid::arrow(length = grid::unit(1.4, "mm"),
                        ends = "last", type = "closed")
  ) +
  # Draw PSBS last so overlapping genes cannot obscure its red arrow.
  geom_segment(
    data = track_data[is_psbs == TRUE],
    aes(x = arrow_start, xend = arrow_end, y = strand_row, yend = strand_row),
    colour = "firebrick", linewidth = 1.1,
    arrow = grid::arrow(length = grid::unit(1.4, "mm"),
                        ends = "last", type = "closed")
  ) +
  geom_text(
    data = track_data[is_psbs == TRUE],
    aes(x = midpoint_mb, y = strand_row + 0.35, label = "PSBS"),
    colour = "firebrick", fontface = "bold", size = 2.8
  ) +
  scale_y_continuous(breaks = c(1, 2), labels = c("−", "+")) +
  coord_cartesian(xlim = display_window_mb, ylim = c(0.4, 2.6), expand = FALSE) +
  labs(x = "Chromosome 3 position (Mb)", y = NULL) + locus_theme +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
left_panel <- (manhattan / gene_track) + plot_layout(heights = c(5, 0.8))
right_panel <- (pip_plot / gene_track) + plot_layout(heights = c(5, 0.8))

# Pairwise dosage r^2 among the three CS leads, using the same LD samples.
lead_snps <- psbs.leads[order(cs), variant_id]
lead_indices <- match(lead_snps, local_map$SNP)
if (anyNA(lead_indices)) stop("A lead SNP is absent from the LD matrix.")
lead_ld <- cor(dosages[, lead_indices, drop = FALSE])^2
cs_labels <- paste0("CS ", psbs.leads[order(cs), cs])
dimnames(lead_ld) <- list(cs_labels, cs_labels)
heatmap_data <- as.data.table(as.table(lead_ld))
setnames(heatmap_data, c("CS_row", "CS_column", "LD_r2"))
heatmap_data[, CS_row := factor(CS_row, levels = rev(cs_labels))]
heatmap_data[, CS_column := factor(CS_column, levels = cs_labels)]
ld_heatmap <- ggplot(heatmap_data, aes(CS_column, CS_row, fill = LD_r2)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", LD_r2),
                colour = LD_r2 > 0.5), size = 3.2) +
  scale_colour_manual(values = c("FALSE" = "white", "TRUE" = "black"), guide = "none") +
  scale_fill_viridis_c(limits = c(0, 1), guide = "none") +
  coord_fixed() +
  labs(x = NULL, y = NULL, title = expression("Lead SNP LD (" * r^2 * ")")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2))
# Show fully homozygous combinations (0 or 2 at all three SNPs) with at
# least 30 samples, matching the original panel I's minimum group size.
genotype_columns <- paste0("genotype_", lead_snps)
all_haplotypes <- psbs[
  expression_index[keep_samples],
  c("expression_log2_cpm", genotype_columns), with = FALSE
]
if (anyNA(all_haplotypes)) stop("Missing expression or lead genotypes.")
homozygous_samples <- all_haplotypes[
  , Reduce(`&`, lapply(.SD, function(dosage) dosage %in% c(0, 2))),
  .SDcols = genotype_columns
]
all_haplotypes <- all_haplotypes[homozygous_samples]
all_haplotypes[, haplotype := do.call(paste, c(.SD, sep = "/")),
               .SDcols = genotype_columns]
combination_counts <- all_haplotypes[
  , .(n_samples = .N), by = haplotype
]
retained_combinations <- combination_counts[
  n_samples >= 30, haplotype
]
all_haplotypes <- all_haplotypes[haplotype %in% retained_combinations]
all_haplotypes[, haplotype := factor(haplotype, levels = sort(unique(haplotype)))]
haplotype_plot <- ggplot(all_haplotypes,
                         aes(haplotype, expression_log2_cpm, fill = haplotype)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0, seed = 42),
              size = 0.7, alpha = 0.4, colour = "grey20") +
  scale_fill_viridis_d(guide = "none") +
  labs(x = "Three-SNP genotypes\n(CS1 / CS2 / CS3)",
       y = expression(log[2]~CPM)) + locus_theme +
  theme(axis.text.x = element_text(size = 8,
                                   angle = if (nlevels(all_haplotypes$haplotype) > 4) 60 else 0,
                                   hjust = if (nlevels(all_haplotypes$haplotype) > 4) 1 else 0.5),
        axis.title.x = element_text(size = 9))
small_heatmap <- (plot_spacer() | ld_heatmap | plot_spacer()) +
  plot_layout(widths = c(0.3, 1, 0.3))
compact_haplotypes <- (plot_spacer() | haplotype_plot | plot_spacer()) +
  plot_layout(widths = c(0.1, 1, 0.1))
heatmap_panel <- (small_heatmap / compact_haplotypes) +
  plot_layout(heights = c(1, 1.25))
keep_haplotypes <- "--keep-haplotypes" %in% commandArgs(trailingOnly = TRUE)
figure <- if (keep_haplotypes) {
  (left_panel | (g_psbs_haplotypes + locus_theme) | right_panel | heatmap_panel) +
    plot_layout(widths = c(1.5, 1, 1.5, 1.1))
} else (left_panel | right_panel | heatmap_panel) +
  plot_layout(widths = c(1, 1, 0.65))
output_dir <- "results/figures/presentation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_file <- file.path(output_dir, "Fig2_PSBS_LD_PIP.png")
ggsave(figure_file, figure, width = if (keep_haplotypes) 42 else 34,
       height = 10, units = "cm", dpi = 300)
fwrite(ld[, .(SNP, LD_r2, reference_SNP = reference,
              n_samples = sum(keep_samples))],
       file.path(output_dir, "Fig2_PSBS_LD.tsv"), sep = "\t")
fwrite(data.table(SNP = lead_snps, as.data.table(lead_ld)),
       file.path(output_dir, "Fig2_PSBS_lead_LD.tsv"), sep = "\t")
fwrite(all_haplotypes[, .(n_samples = .N,
                          mean_expression = mean(expression_log2_cpm),
                          median_expression = median(expression_log2_cpm)),
                      by = haplotype],
       file.path(output_dir, "Fig2_PSBS_haplotype_summary.tsv"), sep = "\t")
message("Saved: ", figure_file, "; LD reference: ", reference,
        "; samples: ", sum(keep_samples))
