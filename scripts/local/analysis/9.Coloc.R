library(rMVP)
library(susieR)
library(tidyverse)
library(data.table)
library(dtplyr)
library(patchwork)

geno <- attach.big.matrix("data/processed/geno/rMVP/Dib639_MAF005.geno.desc")
kin <- attach.big.matrix("data/processed/geno/rMVP/Dib639_MAF005.kin.desc")
pc <- attach.big.matrix("data/processed/geno/rMVP/Dib639_MAF005.pc.desc")[]
map <- fread("data/processed/geno/rMVP/Dib639_MAF005.geno.map")
ind <- fread("data/processed/geno/rMVP/Div639.geno.ind")
colnames(ind)[1] <- "Taxa"

blup <- fread("data/raw/meta/BLUP2020.csv")

rna <- fread("data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz")
colnames(rna)[1] <- "Taxa"

pc.rna <- fread("data/processed/osca/CovariateForEQTL.csv")
pc.rna <- as.matrix(pc.rna)

source("scripts/utilis/functions.R")


###
pc.rna[is.na(pc.rna[,"DaysToPollen"]),"DaysToPollen"] <- mean(pc.rna[,"DaysToPollen"], na.rm=T)

###
psbs.tss <- 179564407
mads.tss <- 161171989
chr <- 3
zcn8.gene <- "Zm00001eb353250"
zcn8.tss <- 126678681
zcn8.chr <- 8

###
blup <- plyr::join(ind, blup, by="Taxa")
rna <- plyr::join(ind, rna, by="Taxa")

mean_impute <- function(x, label) {
  missing <- !is.finite(x)
  if (all(missing)) {
    stop(label, " contains no finite values.")
  }
  if (any(missing)) {
    x[missing] <- mean(x[!missing])
    message("Mean-imputed ", sum(missing), " missing ", label, " value(s).")
  }
  x
}

blup$DaysToSilk <- mean_impute(blup$DaysToSilk, "DaysToSilk")
blup$NPQmaxme <- mean_impute(blup$NPQmaxme, "NPQmaxme")
blup$npqslop2me <- mean_impute(blup$npqslop2me, "npqslop2me")

local_window <- 2e6

subset_gwas_window <- function(gwas, chromosome, tss, window = local_window) {
  gwas[
    CHROM == chromosome &
      POS >= tss - window &
      POS <= tss + window
  ]
}

combine_mvp_results <- function(result) {
  cbind(as.data.table(result$map), as.data.table(result$mlm.results))
}

coloc_out_dir <- "results/full/coloc"
coloc_figure_dir <- file.path(coloc_out_dir, "figures")
dir.create(coloc_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(coloc_figure_dir, recursive = TRUE, showWarnings = FALSE)

save_full_mlm_manhattan <- function(
    gwas,
    p_col,
    trait,
    output_path,
    background_step = 50L,
    keep_p = 1e-5
) {
  p <- gwas[[p_col]]
  keep <- seq.int(1L, nrow(gwas), by = background_step)
  keep <- sort(unique(c(keep, which(is.finite(p) & p <= keep_p))))

  plot_data <- gwas[keep, .(
    CHROM,
    POS,
    P = pmax(get(p_col), .Machine$double.xmin)
  )]
  chr_info <- gwas[, .(chr_max = max(POS, na.rm = TRUE)), by = CHROM]
  setorder(chr_info, CHROM)
  chr_info[, offset := shift(cumsum(chr_max), fill = 0)]
  chr_info[, center := offset + chr_max / 2]
  plot_data[chr_info, offset := i.offset, on = "CHROM"]
  plot_data[, plot_pos := POS + offset]

  plot <- ggplot(
    plot_data,
    aes(
      x = plot_pos,
      y = -log10(P),
      colour = factor(as.integer(factor(CHROM)) %% 2L)
    )
  ) +
    geom_point(size = 0.45, alpha = 0.7) +
    geom_hline(
      yintercept = -log10(0.05 / nrow(gwas)),
      linetype = "dashed",
      colour = "firebrick"
    ) +
    scale_x_continuous(
      breaks = chr_info$center,
      labels = chr_info$CHROM,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_colour_manual(values = c("#3366AA", "#7A9CC6")) +
    labs(
      title = paste0(trait, " MLM GWAS"),
      x = "Chromosome",
      y = expression(-log[10](P))
    ) +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(output_path, plot, width = 12, height = 5.5, dpi = 300)
  invisible(NULL)
}

### GWAS DTS
gwas.dts <- MVP(
  phe = blup[,c("Taxa", "DaysToSilk")],
  geno = geno,
  map = map,
  K = kin,
  CV.MLM = pc,
  ncpus = 24,
  vc.method = "EMMA",
  method = "MLM",
  file.output = FALSE
)

gwas.dts <- combine_mvp_results(gwas.dts)
save_full_mlm_manhattan(
  gwas.dts,
  "DaysToSilk.MLM",
  "DaysToSilk",
  file.path(coloc_figure_dir, "DaysToSilk_MLM_full_manhattan.png")
)
gwas.dts.local <- subset_gwas_window(gwas.dts, chr, mads.tss)
gwas.dts.zcn8.local <- subset_gwas_window(gwas.dts, zcn8.chr, zcn8.tss)
rm(gwas.dts)

### GWAS NPQ
gwas.npq <- MVP(
  phe = blup[,c("Taxa", "NPQmaxme")],
  geno = geno,
  map = map,
  K = kin,
  CV.MLM = pc,
  ncpus = 24,
  vc.method = "EMMA",
  method = "MLM",
  file.output = FALSE
)

gwas.npq <- combine_mvp_results(gwas.npq)
save_full_mlm_manhattan(
  gwas.npq,
  "NPQmaxme.MLM",
  "NPQmaxme",
  file.path(coloc_figure_dir, "NPQmaxme_MLM_full_manhattan.png")
)
gwas.npq.local <- subset_gwas_window(gwas.npq, chr, psbs.tss)
rm(gwas.npq)

### GWAS NPQ slope 2
gwas.npqslope2 <- MVP(
  phe = blup[,c("Taxa", "npqslop2me")],
  geno = geno,
  map = map,
  K = kin,
  CV.MLM = pc,
  ncpus = 24,
  vc.method = "EMMA",
  method = "MLM",
  file.output = FALSE
)

gwas.npqslope2 <- combine_mvp_results(gwas.npqslope2)
save_full_mlm_manhattan(
  gwas.npqslope2,
  "npqslop2me.MLM",
  "NPQ slope 2",
  file.path(coloc_figure_dir, "NPQslope2me_MLM_full_manhattan.png")
)
gwas.npqslope2.local <- subset_gwas_window(gwas.npqslope2, chr, psbs.tss)
rm(gwas.npqslope2)

### GWAS MADS69
gwas.mads <- MVP(
  phe = rna[,c("Taxa", "Zm00001eb143080")],
  geno = geno,
  map = map,
  K = kin,
  CV.MLM = pc.rna[,-4],
  ncpus = 24,
  vc.method = "EMMA",
  method = "MLM",
  file.output = FALSE
)

gwas.mads <- combine_mvp_results(gwas.mads)
save_full_mlm_manhattan(
  gwas.mads,
  "Zm00001eb143080.MLM",
  "MADS69 expression",
  file.path(coloc_figure_dir, "MADS69_expression_MLM_full_manhattan.png")
)
gwas.mads.local <- subset_gwas_window(gwas.mads, chr, mads.tss)
rm(gwas.mads)

### GWAS PSBS
gwas.psbs <- MVP(
  phe = rna[,c("Taxa", "Zm00001eb146510")],
  geno = geno,
  map = map,
  K = kin,
  CV.MLM = pc.rna[,-4],
  ncpus = 24,
  vc.method = "EMMA",
  method = "MLM",
  file.output = FALSE
)

gwas.psbs <- combine_mvp_results(gwas.psbs)
save_full_mlm_manhattan(
  gwas.psbs,
  "Zm00001eb146510.MLM",
  "PSBS expression",
  file.path(coloc_figure_dir, "PSBS_expression_MLM_full_manhattan.png")
)
gwas.psbs.local <- subset_gwas_window(gwas.psbs, chr, psbs.tss)
rm(gwas.psbs)

### GWAS ZCN8 expression
gwas.zcn8 <- MVP(
  phe = as.data.frame(rna)[, c("Taxa", zcn8.gene), drop = FALSE],
  geno = geno,
  map = map,
  K = kin,
  CV.MLM = pc.rna[,-4],
  ncpus = 24,
  vc.method = "EMMA",
  method = "MLM",
  file.output = FALSE
)

gwas.zcn8 <- combine_mvp_results(gwas.zcn8)
save_full_mlm_manhattan(
  gwas.zcn8,
  paste0(zcn8.gene, ".MLM"),
  "ZCN8 expression",
  file.path(coloc_figure_dir, "ZCN8_expression_MLM_full_manhattan.png")
)
gwas.zcn8.local <- subset_gwas_window(gwas.zcn8, zcn8.chr, zcn8.tss)
rm(gwas.zcn8)

### Local genotypes and LD relative to the top expression SNP
variants.mads69 <- map %>%
  filter(
    CHROM == chr,
    POS >= mads.tss - local_window,
    POS <= mads.tss + local_window
  )
X.mads69 <- geno[, which(map$SNP %in% variants.mads69$SNP)]
colnames(X.mads69) <- variants.mads69$SNP

variants.psbs <- map %>%
  filter(
    CHROM == chr,
    POS >= psbs.tss - local_window,
    POS <= psbs.tss + local_window
  )
X.psbs <- geno[, which(map$SNP %in% variants.psbs$SNP)]
colnames(X.psbs) <- variants.psbs$SNP

variants.zcn8 <- map %>%
  filter(
    CHROM == zcn8.chr,
    POS >= zcn8.tss - local_window,
    POS <= zcn8.tss + local_window
  )
X.zcn8 <- geno[, which(map$SNP %in% variants.zcn8$SNP)]
colnames(X.zcn8) <- variants.zcn8$SNP

top_expression_snp <- function(gwas, p_col) {
  p <- gwas[[p_col]]
  p[!is.finite(p)] <- Inf
  if (all(is.infinite(p))) {
    stop("No finite expression GWAS p-values are available.")
  }
  gwas$SNP[which.min(p)]
}

ld_to_top_snp <- function(X, top_snp) {
  top_idx <- match(top_snp, colnames(X))
  if (is.na(top_idx)) {
    stop("Top expression SNP is absent from the local genotype matrix: ", top_snp)
  }

  r2 <- as.numeric(stats::cor(
    X,
    X[, top_idx],
    use = "pairwise.complete.obs"
  ))^2
  r2[!is.finite(r2)] <- NA_real_
  data.table(SNP = colnames(X), LD_r2 = r2)
}

top.mads69.snp <- top_expression_snp(gwas.mads.local, "Zm00001eb143080.MLM")
top.psbs.snp <- top_expression_snp(gwas.psbs.local, "Zm00001eb146510.MLM")
top.zcn8.snp <- top_expression_snp(gwas.zcn8.local, paste0(zcn8.gene, ".MLM"))

ld.mads69 <- ld_to_top_snp(X.mads69, top.mads69.snp)
ld.psbs <- ld_to_top_snp(X.psbs, top.psbs.snp)
ld.zcn8 <- ld_to_top_snp(X.zcn8, top.zcn8.snp)

gwas.dts.local[ld.mads69, LD_r2 := i.LD_r2, on = "SNP"]
gwas.mads.local[ld.mads69, LD_r2 := i.LD_r2, on = "SNP"]
gwas.npq.local[ld.psbs, LD_r2 := i.LD_r2, on = "SNP"]
gwas.npqslope2.local[ld.psbs, LD_r2 := i.LD_r2, on = "SNP"]
gwas.psbs.local[ld.psbs, LD_r2 := i.LD_r2, on = "SNP"]
gwas.dts.zcn8.local[ld.zcn8, LD_r2 := i.LD_r2, on = "SNP"]
gwas.zcn8.local[ld.zcn8, LD_r2 := i.LD_r2, on = "SNP"]

gwas.dts.local[, is_top_expression_snp := SNP == top.mads69.snp]
gwas.mads.local[, is_top_expression_snp := SNP == top.mads69.snp]
gwas.npq.local[, is_top_expression_snp := SNP == top.psbs.snp]
gwas.npqslope2.local[, is_top_expression_snp := SNP == top.psbs.snp]
gwas.psbs.local[, is_top_expression_snp := SNP == top.psbs.snp]
gwas.dts.zcn8.local[, is_top_expression_snp := SNP == top.zcn8.snp]
gwas.zcn8.local[, is_top_expression_snp := SNP == top.zcn8.snp]

### Local GWAS results and plots
local_manhattan_pair <- function(
    gwas1,
    p_col1,
    trait1,
    gwas2,
    p_col2,
    trait2,
    tss,
    gene,
    chromosome,
    ld_col = "LD_r2",
    top_col = "is_top_expression_snp"
) {
  plot_data <- rbindlist(list(
    gwas1[, .(
      SNP,
      POS,
      P = pmax(get(p_col1), .Machine$double.xmin),
      LD_r2 = get(ld_col),
      is_top_expression_snp = get(top_col),
      trait = trait1
    )],
    gwas2[, .(
      SNP,
      POS,
      P = pmax(get(p_col2), .Machine$double.xmin),
      LD_r2 = get(ld_col),
      is_top_expression_snp = get(top_col),
      trait = trait2
    )]
  ))
  plot_data[, trait := factor(trait, levels = c(trait1, trait2))]

  ggplot(plot_data, aes(x = POS / 1e6, y = -log10(P), colour = LD_r2)) +
    geom_point(size = 0.8, alpha = 0.8) +
    geom_point(
      data = plot_data[is_top_expression_snp == TRUE],
      aes(x = POS / 1e6, y = -log10(P)),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = 1,
      fill = NA,
      colour = "black"
    ) +
    geom_vline(
      xintercept = tss / 1e6,
      linetype = "dashed",
      colour = "firebrick"
    ) +
    annotate(
      "text",
      x = tss / 1e6,
      y = Inf,
      label = paste0(gene, " TSS"),
      colour = "firebrick",
      angle = 90,
      hjust = 1.1,
      vjust = -0.3,
      size = 3
    ) +
    facet_wrap(vars(trait), ncol = 1, scales = "free_y") +
    scale_colour_viridis_c(
      limits = c(0, 1),
      option = "D",
      name = expression(LD ~ r^2),
      na.value = "grey80"
    ) +
    labs(
      x = paste0("Chromosome ", chromosome, " position (Mb)"),
      y = expression(-log[10](P))
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
}

paired_gwas_scatter <- function(
    gwas1,
    p_col1,
    trait1,
    gwas2,
    p_col2,
    trait2,
    gene,
    ld_col = "LD_r2",
    top_col = "is_top_expression_snp"
) {
  paired <- merge(
    gwas1[, .(
      SNP,
      p1 = pmax(get(p_col1), .Machine$double.xmin),
      LD_r2 = get(ld_col),
      is_top_expression_snp = get(top_col)
    )],
    gwas2[, .(SNP, p2 = pmax(get(p_col2), .Machine$double.xmin))],
    by = "SNP"
  )

  plot <- ggplot(paired, aes(x = -log10(p1), y = -log10(p2))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(aes(colour = LD_r2), size = 0.9, alpha = 0.75) +
    geom_point(
      data = paired[is_top_expression_snp == TRUE],
      shape = 21,
      size = 3,
      stroke = 1,
      fill = NA,
      colour = "black"
    ) +
    scale_colour_viridis_c(
      limits = c(0, 1),
      option = "D",
      name = expression(LD ~ r^2),
      na.value = "grey80"
    ) +
    labs(
      title = paste0(gene, " local GWAS comparison"),
      x = bquote(-log[10](P) ~ .(trait1)),
      y = bquote(-log[10](P) ~ .(trait2))
    ) +
    theme_classic() +
    theme(legend.position = "bottom")

  list(data = paired, plot = plot)
}

g.local.mads.dts <- local_manhattan_pair(
  gwas.dts.local, "DaysToSilk.MLM", "DaysToSilk",
  gwas.mads.local, "Zm00001eb143080.MLM", "MADS69 expression",
  mads.tss, "MADS69", chr
)
g.local.psbs.npq <- local_manhattan_pair(
  gwas.npq.local, "NPQmaxme.MLM", "NPQmaxme",
  gwas.psbs.local, "Zm00001eb146510.MLM", "PSBS expression",
  psbs.tss, "PSBS", chr
)
g.local.psbs.npqslope2 <- local_manhattan_pair(
  gwas.npqslope2.local, "npqslop2me.MLM", "NPQ slope 2",
  gwas.psbs.local, "Zm00001eb146510.MLM", "PSBS expression",
  psbs.tss, "PSBS", chr
)
g.local.zcn8.dts <- local_manhattan_pair(
  gwas.dts.zcn8.local, "DaysToSilk.MLM", "DaysToSilk",
  gwas.zcn8.local, paste0(zcn8.gene, ".MLM"), "ZCN8 expression",
  zcn8.tss, "ZCN8", zcn8.chr
)

g.scatter.mads.dts <- paired_gwas_scatter(
  gwas.dts.local, "DaysToSilk.MLM", "DaysToSilk",
  gwas.mads.local, "Zm00001eb143080.MLM", "MADS69 expression",
  "MADS69"
)
g.scatter.psbs.npq <- paired_gwas_scatter(
  gwas.npq.local, "NPQmaxme.MLM", "NPQmaxme",
  gwas.psbs.local, "Zm00001eb146510.MLM", "PSBS expression",
  "PSBS"
)
g.scatter.psbs.npqslope2 <- paired_gwas_scatter(
  gwas.npqslope2.local, "npqslop2me.MLM", "NPQ slope 2",
  gwas.psbs.local, "Zm00001eb146510.MLM", "PSBS expression",
  "PSBS"
)
g.scatter.zcn8.dts <- paired_gwas_scatter(
  gwas.dts.zcn8.local, "DaysToSilk.MLM", "DaysToSilk",
  gwas.zcn8.local, paste0(zcn8.gene, ".MLM"), "ZCN8 expression",
  "ZCN8"
)

g.combined.mads.dts <- (
  g.local.mads.dts | g.scatter.mads.dts$plot
) +
  plot_layout(widths = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")

g.combined.psbs.npq <- (
  g.local.psbs.npq | g.scatter.psbs.npq$plot
) +
  plot_layout(widths = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")

g.combined.psbs.npqslope2 <- (
  g.local.psbs.npqslope2 | g.scatter.psbs.npqslope2$plot
) +
  plot_layout(widths = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")

g.combined.zcn8.dts <- (
  g.local.zcn8.dts | g.scatter.zcn8.dts$plot
) +
  plot_layout(widths = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")

fwrite(gwas.dts.local, file.path(coloc_out_dir, "gwas_DaysToSilk_MADS69_local.tsv.gz"), sep = "\t")
fwrite(gwas.dts.zcn8.local, file.path(coloc_out_dir, "gwas_DaysToSilk_ZCN8_local.tsv.gz"), sep = "\t")
fwrite(gwas.mads.local, file.path(coloc_out_dir, "gwas_MADS69_expression_local.tsv.gz"), sep = "\t")
fwrite(gwas.npq.local, file.path(coloc_out_dir, "gwas_NPQmaxme_PSBS_local.tsv.gz"), sep = "\t")
fwrite(gwas.npqslope2.local, file.path(coloc_out_dir, "gwas_NPQslope2me_PSBS_local.tsv.gz"), sep = "\t")
fwrite(gwas.psbs.local, file.path(coloc_out_dir, "gwas_PSBS_expression_local.tsv.gz"), sep = "\t")
fwrite(gwas.zcn8.local, file.path(coloc_out_dir, "gwas_ZCN8_expression_local.tsv.gz"), sep = "\t")
fwrite(
  data.table(
    gene = c("MADS69", "PSBS", "ZCN8"),
    top_expression_snp = c(top.mads69.snp, top.psbs.snp, top.zcn8.snp)
  ),
  file.path(coloc_out_dir, "top_expression_snps.tsv"),
  sep = "\t"
)

ggsave(
  file.path(coloc_figure_dir, "MADS69_DaysToSilk_local_combined.png"),
  g.combined.mads.dts,
  width = 13,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "PSBS_NPQmaxme_local_combined.png"),
  g.combined.psbs.npq,
  width = 13,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "PSBS_NPQslope2me_local_combined.png"),
  g.combined.psbs.npqslope2,
  width = 13,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "ZCN8_DaysToSilk_local_combined.png"),
  g.combined.zcn8.dts,
  width = 13,
  height = 7,
  dpi = 300
)

### remove covariance effects
fit_susie_complete <- function(X, Z, y, snp, L = 10) {
  missing_y <- !is.finite(y)
  if (all(missing_y)) {
    stop("The phenotype contains no finite values.")
  }
  if (any(missing_y)) {
    y[missing_y] <- mean(y[!missing_y])
    message("Mean-imputed ", sum(missing_y), " missing phenotype value(s).")
  }

  keep <- apply(Z, 1, function(z) all(is.finite(z)))
  if (sum(keep) < 2L) {
    stop("Too few complete samples for SuSiE.")
  }

  adjusted <- remove.covariate.effects(
    X = X[keep, , drop = FALSE],
    Z = Z[keep, , drop = FALSE],
    y = y[keep]
  )
  colnames(adjusted$X) <- snp

  fit <- susie(adjusted$X, adjusted$y, L = L)
  annotate_susie_fit(fit, snp = snp)
}

### run SuSiE
res.mads69 <- fit_susie_complete(
  X.mads69, pc.rna[,-4], rna$Zm00001eb143080, variants.mads69$SNP
)
res.psbs <- fit_susie_complete(
  X.psbs, pc.rna[,-4], rna$Zm00001eb146510, variants.psbs$SNP
)
res.zcn8 <- fit_susie_complete(
  X.zcn8, pc.rna[,-4], rna[[zcn8.gene]], variants.zcn8$SNP
)
res.dts <- fit_susie_complete(
  X.mads69, pc, blup$DaysToSilk, variants.mads69$SNP
)
res.dts.zcn8 <- fit_susie_complete(
  X.zcn8, pc, blup$DaysToSilk, variants.zcn8$SNP
)
res.npq <- fit_susie_complete(
  X.psbs, pc, blup$NPQmaxme, variants.psbs$SNP
)
res.npqslope2 <- fit_susie_complete(
  X.psbs, pc, blup$npqslop2me, variants.psbs$SNP
)

susie_pip_plot <- function(fit, variants, trait, gene, tss, chromosome) {
  pip_data <- data.table(
    SNP = names(fit$pip),
    PIP = as.numeric(fit$pip)
  )
  pip_data[variants[, .(SNP, CHROM, POS)], `:=`(
    CHROM = i.CHROM,
    POS = i.POS
  ), on = "SNP"]

  pip_data[, credible_set := "Not in CS"]
  if (length(fit$sets$cs)) {
    for (i in seq_along(fit$sets$cs)) {
      index <- fit$sets$cs[[i]]
      snp <- names(index)
      if (is.null(snp)) {
        snp <- names(fit$pip)[index]
      }
      pip_data[SNP %in% snp, credible_set := paste0("CS", i)]
    }
  }
  cs_levels <- c("Not in CS", paste0("CS", seq_along(fit$sets$cs)))
  pip_data[, credible_set := factor(credible_set, levels = cs_levels)]
  cs_palette <- if (length(fit$sets$cs)) {
    viridisLite::viridis(
      length(fit$sets$cs),
      option = "D",
      begin = 0.1,
      end = 0.9
    )
  } else {
    character()
  }
  cs_colours <- c(
    "Not in CS" = "grey70",
    stats::setNames(
      cs_palette,
      paste0("CS", seq_along(fit$sets$cs))
    )
  )

  top <- pip_data[which.max(PIP)]
  top_label <- paste0(top$SNP, "\nPIP = ", signif(top$PIP, 3))

  plot <- ggplot(pip_data, aes(x = POS / 1e6, y = PIP)) +
    geom_point(aes(colour = credible_set), size = 0.8, alpha = 0.75) +
    geom_vline(
      xintercept = tss / 1e6,
      linetype = "dashed",
      colour = "firebrick"
    ) +
    geom_point(
      data = top,
      shape = 21,
      size = 3.2,
      stroke = 1,
      fill = NA,
      colour = "black"
    ) +
    geom_text(
      data = top,
      aes(label = top_label),
      hjust = -0.05,
      vjust = -0.4,
      size = 3
    ) +
    scale_colour_manual(
      values = cs_colours,
      name = "95% credible set",
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.01, 0.15))
    ) +
    labs(
      title = paste0(trait, " SuSiE fine-mapping"),
      subtitle = paste0(gene, " ±2 Mb"),
      x = paste0("Chromosome ", chromosome, " position (Mb)"),
      y = "Posterior inclusion probability"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.margin = margin(5.5, 65, 5.5, 5.5)
    )

  list(data = pip_data, top = top, plot = plot)
}

g.susie.mads69 <- susie_pip_plot(
  res.mads69, variants.mads69, "MADS69 expression", "MADS69", mads.tss, chr
)
g.susie.psbs <- susie_pip_plot(
  res.psbs, variants.psbs, "PSBS expression", "PSBS", psbs.tss, chr
)
g.susie.zcn8 <- susie_pip_plot(
  res.zcn8, variants.zcn8, "ZCN8 expression", "ZCN8", zcn8.tss, zcn8.chr
)
g.susie.dts <- susie_pip_plot(
  res.dts, variants.mads69, "DaysToSilk", "MADS69", mads.tss, chr
)
g.susie.dts.zcn8 <- susie_pip_plot(
  res.dts.zcn8, variants.zcn8, "DaysToSilk", "ZCN8", zcn8.tss, zcn8.chr
)
g.susie.npq <- susie_pip_plot(
  res.npq, variants.psbs, "NPQmaxme", "PSBS", psbs.tss, chr
)
g.susie.npqslope2 <- susie_pip_plot(
  res.npqslope2, variants.psbs, "NPQ slope 2", "PSBS", psbs.tss, chr
)

ggsave(
  file.path(coloc_figure_dir, "MADS69_expression_SuSiE_PIP.png"),
  g.susie.mads69$plot,
  width = 8,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "PSBS_expression_SuSiE_PIP.png"),
  g.susie.psbs$plot,
  width = 8,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "ZCN8_expression_SuSiE_PIP.png"),
  g.susie.zcn8$plot,
  width = 8,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "DaysToSilk_MADS69_SuSiE_PIP.png"),
  g.susie.dts$plot,
  width = 8,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "NPQmaxme_PSBS_SuSiE_PIP.png"),
  g.susie.npq$plot,
  width = 8,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "NPQslope2me_PSBS_SuSiE_PIP.png"),
  g.susie.npqslope2$plot,
  width = 8,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(coloc_figure_dir, "DaysToSilk_ZCN8_SuSiE_PIP.png"),
  g.susie.dts.zcn8$plot,
  width = 8,
  height = 4.5,
  dpi = 300
)

fwrite(
  rbindlist(list(
    g.susie.mads69$top[, trait := "MADS69_expression"],
    g.susie.psbs$top[, trait := "PSBS_expression"],
    g.susie.zcn8$top[, trait := "ZCN8_expression"],
    g.susie.dts$top[, trait := "DaysToSilk"],
    g.susie.dts.zcn8$top[, trait := "DaysToSilk_ZCN8_locus"],
    g.susie.npq$top[, trait := "NPQmaxme"],
    g.susie.npqslope2$top[, trait := "NPQslope2me"]
  ), fill = TRUE),
  file.path(coloc_out_dir, "susie_top_pip_snps.tsv"),
  sep = "\t"
)

### run coloc

coloc.flower <- run_coloc_susie_raw(res.dts, res.mads69)

coloc.npq <- run_coloc_susie_raw(res.psbs, res.npq)
coloc.npqslope2 <- run_coloc_susie_raw(res.psbs, res.npqslope2)
coloc.zcn8 <- run_coloc_susie_raw(res.dts.zcn8, res.zcn8)

run_coloc_abf_gwas <- function(
    gwas1,
    beta_col1,
    se_col1,
    sdY1,
    N1,
    gwas2,
    beta_col2,
    se_col2,
    sdY2,
    N2
) {
  aligned <- merge(
    gwas1[, .(
      snp = SNP,
      position = POS,
      beta1 = get(beta_col1),
      varbeta1 = get(se_col1)^2
    )],
    gwas2[, .(
      snp = SNP,
      beta2 = get(beta_col2),
      varbeta2 = get(se_col2)^2
    )],
    by = "snp"
  )

  coloc::coloc.abf(
    dataset1 = list(
      beta = aligned$beta1,
      varbeta = aligned$varbeta1,
      snp = aligned$snp,
      position = aligned$position,
      type = "quant",
      N = N1,
      sdY = sdY1
    ),
    dataset2 = list(
      beta = aligned$beta2,
      varbeta = aligned$varbeta2,
      snp = aligned$snp,
      type = "quant",
      N = N2,
      sdY = sdY2
    )
  )
}

coloc.abf.flower.mlm <- run_coloc_abf_gwas(
  gwas.dts.local, "Effect", "SE", sd(blup$DaysToSilk), nrow(ind),
  gwas.mads.local, "Effect", "SE", sd(rna$Zm00001eb143080), nrow(ind)
)
coloc.abf.npq.mlm <- run_coloc_abf_gwas(
  gwas.npq.local, "Effect", "SE", sd(blup$NPQmaxme), nrow(ind),
  gwas.psbs.local, "Effect", "SE", sd(rna$Zm00001eb146510), nrow(ind)
)
coloc.abf.npqslope2.mlm <- run_coloc_abf_gwas(
  gwas.npqslope2.local, "Effect", "SE", sd(blup$npqslop2me), nrow(ind),
  gwas.psbs.local, "Effect", "SE", sd(rna$Zm00001eb146510), nrow(ind)
)
coloc.abf.zcn8.mlm <- run_coloc_abf_gwas(
  gwas.dts.zcn8.local, "Effect", "SE", sd(blup$DaysToSilk), nrow(ind),
  gwas.zcn8.local, "Effect", "SE", sd(rna[[zcn8.gene]]), nrow(ind)
)

print(coloc.flower$summary)
print(coloc.npq$summary)
print(coloc.npqslope2$summary)
print(coloc.zcn8$summary)

fwrite(
  coloc.flower$summary,
  file.path(coloc_out_dir, "coloc_DaysToSilk_MADS69_summary.tsv"),
  sep = "\t"
)
fwrite(
  coloc.flower$results,
  file.path(coloc_out_dir, "coloc_DaysToSilk_MADS69_snp_posterior.tsv.gz"),
  sep = "\t"
)
fwrite(
  coloc.npq$summary,
  file.path(coloc_out_dir, "coloc_PSBS_NPQmaxme_summary.tsv"),
  sep = "\t"
)
fwrite(
  coloc.npq$results,
  file.path(coloc_out_dir, "coloc_PSBS_NPQmaxme_snp_posterior.tsv.gz"),
  sep = "\t"
)
fwrite(
  coloc.npqslope2$summary,
  file.path(coloc_out_dir, "coloc_PSBS_NPQslope2me_summary.tsv"),
  sep = "\t"
)
fwrite(
  coloc.npqslope2$results,
  file.path(coloc_out_dir, "coloc_PSBS_NPQslope2me_snp_posterior.tsv.gz"),
  sep = "\t"
)
fwrite(
  coloc.zcn8$summary,
  file.path(coloc_out_dir, "coloc_DaysToSilk_ZCN8_summary.tsv"),
  sep = "\t"
)
fwrite(
  coloc.zcn8$results,
  file.path(coloc_out_dir, "coloc_DaysToSilk_ZCN8_snp_posterior.tsv.gz"),
  sep = "\t"
)
saveRDS(
  list(
    flower = coloc.flower,
    npq = coloc.npq,
    npqslope2 = coloc.npqslope2,
    zcn8 = coloc.zcn8
  ),
  file.path(coloc_out_dir, "coloc_results.rds")
)

format_abf_summary <- function(result, model, pair) {
  out <- as.data.table(as.list(result$summary))
  out[, `:=`(model = model, pair = pair)]
  setcolorder(out, c("model", "pair"))
  out
}

abf_summary <- rbindlist(list(
  format_abf_summary(coloc.abf.flower.mlm, "MLM", "DaysToSilk_MADS69"),
  format_abf_summary(coloc.abf.npq.mlm, "MLM", "PSBS_NPQmaxme"),
  format_abf_summary(coloc.abf.npqslope2.mlm, "MLM", "PSBS_NPQslope2me"),
  format_abf_summary(coloc.abf.zcn8.mlm, "MLM", "DaysToSilk_ZCN8")
), fill = TRUE)
fwrite(
  abf_summary,
  file.path(coloc_out_dir, "coloc_MLM_ABF_summary.tsv"),
  sep = "\t"
)
saveRDS(
  list(
    flower_mlm = coloc.abf.flower.mlm,
    npq_mlm = coloc.abf.npq.mlm,
    npqslope2_mlm = coloc.abf.npqslope2.mlm,
    zcn8_mlm = coloc.abf.zcn8.mlm
  ),
  file.path(coloc_out_dir, "coloc_MLM_ABF_results.rds")
)
