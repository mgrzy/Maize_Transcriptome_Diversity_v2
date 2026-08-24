#!/usr/bin/env Rscript
# Purpose: Test differences between genes in the top GERP decile and all other
#          genes for (1) fine-mapped aFC effect sizes and (2) credible-set
#          counts per eGene.
#
# Tests:
#   - Wilcoxon rank-sum test for absolute log2(aFC) of fine-mapped variants.
#   - Quasi-Poisson GLM for the number of credible sets per eGene. The
#     quasi-Poisson family retains the count-model mean structure while using
#     a dispersion estimate for potentially overdispersed credible-set counts.
#
# Inputs:
#   - results/summary/AF/aFCn_MAF_se_filtered.tsv.gz
#   - results/summary/GERP/positive_GERP_mean_by_gene_and_promoter.csv
#   - results/summary/SusieR/CS_Gene.csv
# Output:
#   - results/summary/GERP/GERP_effect_size_and_credible_set_tests.csv
#
# Usage:
#   Rscript scripts/local/analysis/15.Test_GERP_effect_sizes_and_CS.R

suppressPackageStartupMessages(library(data.table))

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  fine_mapped_file = Sys.getenv(
    "FINE_MAPPED_FILE",
    "results/summary/AF/aFCn_MAF_se_filtered.tsv.gz"
  ),
  gerp_file = Sys.getenv(
    "GERP_FILE",
    "results/summary/GERP/positive_GERP_mean_by_gene_and_promoter.csv"
  ),
  credible_set_file = Sys.getenv(
    "CREDIBLE_SET_FILE",
    "results/summary/SusieR/CS_Gene.csv"
  ),
  output_file = Sys.getenv(
    "OUTPUT_FILE",
    "results/summary/GERP/GERP_effect_size_and_credible_set_tests.csv"
  )
)
for (path in unlist(cfg[1:3])) {
  if (!file.exists(path)) stop_with("Required input does not exist: %s", path)
}
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)

fine_mapped <- fread(cfg$fine_mapped_file, select = c("gene_id", "variant_id", "log2_aFC"))
gerp <- fread(cfg$gerp_file, select = c("gene_id", "gene_mean_positive_GERP"))
credible_sets <- fread(cfg$credible_set_file, select = c("gene", "cs"))
setnames(credible_sets, c("gene", "cs"), c("gene_id", "n_credible_sets"))

gerp[, `:=`(
  gene_id = as.character(gene_id),
  gene_mean_positive_GERP = as.numeric(gene_mean_positive_GERP)
)]
gerp <- gerp[is.finite(gene_mean_positive_GERP)]
if (gerp[, anyDuplicated(gene_id)]) stop_with("GERP summary has duplicate gene IDs.")

gerp_threshold <- stats::quantile(gerp$gene_mean_positive_GERP, probs = 0.90, na.rm = TRUE)
gerp[, gerp_group := factor(
  ifelse(gene_mean_positive_GERP >= gerp_threshold, "Top 10% GERP genes", "Other genes"),
  levels = c("Other genes", "Top 10% GERP genes")
)]

effect_data <- merge(
  fine_mapped[, .(
    gene_id = as.character(gene_id),
    variant_id = as.character(variant_id),
    abs_log2_aFC = abs(as.numeric(log2_aFC))
  )],
  gerp[, .(gene_id, gerp_group)],
  by = "gene_id",
  all = FALSE,
  sort = FALSE
)
effect_data <- effect_data[is.finite(abs_log2_aFC)]
if (effect_data[, uniqueN(gerp_group)] != 2L) stop_with("Both GERP groups are required for the effect-size test.")

effect_wilcox <- stats::wilcox.test(
  abs_log2_aFC ~ gerp_group,
  data = effect_data,
  exact = FALSE,
  conf.int = TRUE
)

credible_set_data <- merge(
  credible_sets[, .(
    gene_id = as.character(gene_id),
    n_credible_sets = as.integer(n_credible_sets)
  )],
  gerp[, .(gene_id, gerp_group)],
  by = "gene_id",
  all = FALSE,
  sort = FALSE
)
credible_set_data <- credible_set_data[
  is.finite(n_credible_sets) & n_credible_sets >= 0L
]
if (credible_set_data[, uniqueN(gerp_group)] != 2L) stop_with("Both GERP groups are required for the credible-set GLM.")

cs_glm <- stats::glm(
  n_credible_sets ~ gerp_group,
  family = stats::quasipoisson(link = "log"),
  data = credible_set_data
)
glm_coefficients <- summary(cs_glm)$coefficients
group_term <- "gerp_groupTop 10% GERP genes"
if (!group_term %in% rownames(glm_coefficients)) stop_with("Could not find the top-GERP coefficient in the GLM.")
glm_estimate <- glm_coefficients[group_term, "Estimate"]
glm_se <- glm_coefficients[group_term, "Std. Error"]
glm_p_value <- glm_coefficients[group_term, grep("Pr", colnames(glm_coefficients), value = TRUE)]
glm_ci <- exp(glm_estimate + c(-1, 1) * stats::qnorm(0.975) * glm_se)

test_results <- rbindlist(list(
  data.table(
    test = "Wilcoxon rank-sum",
    outcome = "absolute_log2_aFC",
    comparison = "Top 10% GERP genes vs Other genes",
    gerp_90th_percentile = as.numeric(gerp_threshold),
    n_other_genes = uniqueN(effect_data[gerp_group == "Other genes", gene_id]),
    n_top_gerp_genes = uniqueN(effect_data[gerp_group == "Top 10% GERP genes", gene_id]),
    n_other_observations = nrow(effect_data[gerp_group == "Other genes"]),
    n_top_gerp_observations = nrow(effect_data[gerp_group == "Top 10% GERP genes"]),
    statistic = unname(effect_wilcox$statistic),
    estimate = NA_real_,
    estimate_label = NA_character_,
    conf_low = unname(effect_wilcox$conf.int[1L]),
    conf_high = unname(effect_wilcox$conf.int[2L]),
    p_value = effect_wilcox$p.value,
    dispersion = NA_real_
  ),
  data.table(
    test = "Quasi-Poisson GLM",
    outcome = "credible_set_count",
    comparison = "Top 10% GERP genes vs Other genes",
    gerp_90th_percentile = as.numeric(gerp_threshold),
    n_other_genes = uniqueN(credible_set_data[gerp_group == "Other genes", gene_id]),
    n_top_gerp_genes = uniqueN(credible_set_data[gerp_group == "Top 10% GERP genes", gene_id]),
    n_other_observations = nrow(credible_set_data[gerp_group == "Other genes"]),
    n_top_gerp_observations = nrow(credible_set_data[gerp_group == "Top 10% GERP genes"]),
    statistic = unname(glm_estimate / glm_se),
    estimate = exp(glm_estimate),
    estimate_label = "incidence_rate_ratio",
    conf_low = glm_ci[1L],
    conf_high = glm_ci[2L],
    p_value = unname(glm_p_value),
    dispersion = summary(cs_glm)$dispersion
  )
), use.names = TRUE)

fwrite(test_results, cfg$output_file)
message("Top-decile GERP threshold: ", signif(gerp_threshold, 6))
message("Wrote Wilcoxon and credible-set GLM results to ", cfg$output_file)
