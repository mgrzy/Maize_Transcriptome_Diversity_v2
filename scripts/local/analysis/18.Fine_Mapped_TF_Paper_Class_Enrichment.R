#!/usr/bin/env Rscript
# Run matched-control enrichment of fine-mapped variants in TF paper-class
# union peaks. This configures and reuses the tested implementation in
# 17.Fine_Mapped_TF_ChIP_Enrichment.R, with paper classes rather than
# individual TFs as the chromatin features.
#
# Before running, create the paper-class union BEDs:
#   bash scripts/local/procesing/prepare_Tu_TF_paper_class_union_beds.sh
#
# Defaults (all remain overridable with environment variables):
#   TF_PEAKS_DIR=data/processed/chipseq/Tu_2020_TF_consensus_v5/by_paper_class
#   OUTPUT_DIR=results/summary/fine_mapped_TF_paper_class_enrichment
#   EXPECTED_N_TFS=8
#   FEATURE_OUTPUT_STEM=paper_class
#   FEATURE_RESULT_COLUMN=paper_class
#
# Example:
#   N_PERMUTATIONS=5000 N_CORES=8 \
#     Rscript scripts/local/analysis/18.Fine_Mapped_TF_Paper_Class_Enrichment.R
#
# Outputs include paper_class_overlap_matched_permutation.tsv, annotated
# fine-mapped variants, matched controls, and (optionally) permutation nulls.

set_default_env <- function(name, value) {
  if (!nzchar(Sys.getenv(name))) do.call(Sys.setenv, setNames(list(value), name))
}

script_argument <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", script_argument, value = TRUE)
if (!length(script_flag)) stop("Could not determine this script's path.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_flag[[1L]]))
project_root <- normalizePath(file.path(dirname(script_path), "../../.."))

set_default_env(
  "TF_PEAKS_DIR",
  file.path(project_root, "data/processed/chipseq/Tu_2020_TF_consensus_v5/by_paper_class")
)
set_default_env(
  "OUTPUT_DIR",
  file.path(project_root, "results/summary/fine_mapped_TF_paper_class_enrichment")
)
set_default_env("EXPECTED_N_TFS", "8")
set_default_env("FEATURE_OUTPUT_STEM", "paper_class")
set_default_env("FEATURE_RESULT_COLUMN", "paper_class")
set_default_env(
  "FEATURE_LABEL_MAP_FILE",
  file.path(Sys.getenv("TF_PEAKS_DIR"), "paper_class_union_QC.tsv")
)

source(file.path(dirname(script_path), "17.Fine_Mapped_TF_ChIP_Enrichment.R"), local = FALSE)
