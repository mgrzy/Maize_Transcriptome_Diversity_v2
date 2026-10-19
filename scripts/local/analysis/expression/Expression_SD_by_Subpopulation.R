#!/usr/bin/env Rscript
# Compute per-gene sample SD by subpopulation using the saved VST expression
# and the flowering-time adjustment from 1.Transcriptome_diversity.R.
# Run from the repository root; no PCA or permutation analysis is rerun.
suppressPackageStartupMessages(library(data.table))

rna <- fread("data/processed/rna/counts_NE2020_corrected_names_deseq2_vst_filtered.csv.gz")
setnames(rna, names(rna)[1], "Taxa")
meta <- fread("data/raw/meta/Lines_meta.csv")
blup <- fread("data/raw/meta/BLUP2020.csv")
# Preserve the existing workflow's correction of the truncated sample name.
blup[Taxa == "Hy", Taxa := "ILLHy"]
mean_flowering_time <- mean(blup$DaysToPollen, na.rm = TRUE)
blup[is.na(DaysToPollen), DaysToPollen := mean_flowering_time]
if (anyDuplicated(rna$Taxa) || anyDuplicated(meta$Taxa) || anyDuplicated(blup$Taxa)) {
  stop("Input tables contain duplicated sample IDs.")
}
samples <- merge(meta, rna, by = "Taxa")
samples <- merge(blup[, .(Taxa, DaysToPollen)], samples, by = "Taxa")
gene_columns <- grep("^Zm", names(samples), value = TRUE)
if (!length(gene_columns) || anyNA(samples$Group) || anyNA(samples$DaysToPollen)) {
  stop("Missing genes, groups, or flowering-time covariate.")
}
expression <- as.matrix(samples[, ..gene_columns])
if (any(!is.finite(expression))) stop("Expression contains non-finite values.")
# Matrix QR residualization equals fitting expression ~ DaysToPollen per gene.
adjusted_expression <- qr.resid(qr(model.matrix(~ DaysToPollen, data = samples)), expression)
sd_rows <- list()
cv_check_rows <- list()
for (group in sort(unique(samples$Group))) {
  values <- adjusted_expression[samples$Group == group, , drop = FALSE]
  sd_rows[[group]] <- data.table(Group = group, gene = gene_columns,
                                 sd = apply(values, 2, sd), n_lines = nrow(values))
  absolute_values <- abs(values)
  cv_check_rows[[group]] <- data.table(
    Group = group, name = gene_columns,
    reconstructed_cv = apply(absolute_values, 2, var) / colMeans(absolute_values)
  )
}
# Verify that the reconstructed sample set/adjustment reproduces the old panel.
saved_cv <- fread("results/full/GeneExpressionVariance/All_CV.csv")
checked <- merge(saved_cv, rbindlist(cv_check_rows), by = c("Group", "name"))
if (nrow(checked) != nrow(saved_cv) ||
    any(abs(checked$cv - checked$reconstructed_cv) > 1e-8)) {
  stop("Reconstructed expression does not reproduce the saved variation summary.")
}
result <- rbindlist(sd_rows)
output <- "results/full/GeneExpressionVariance/All_SD.csv"
fwrite(result, output)
message("Saved ", nrow(result), " gene/subpopulation SD values to ", output)
