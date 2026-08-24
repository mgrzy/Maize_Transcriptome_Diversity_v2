#!/usr/bin/env Rscript
# Run leaf-ATAC-matched enrichment of fine-mapped variants in TF paper-class
# union peaks. This uses the matched-control design in
# 19.Fine_Mapped_TF_ChIP_Enrichment_ATAC_Matched.R: each sampled control has
# the same target-TSS-distance bin, MAF bin, variant type, and binary leaf-ATAC
# overlap status as its observed fine-mapped variant.
#
# Before running, create the paper-class union BEDs:
#   bash scripts/local/procesing/prepare_Tu_TF_paper_class_union_beds.sh
#
# LEAF_ATAC_FILE is required and should be in B73 v5 coordinates.
# Example:
#   LEAF_ATAC_FILE=data/processed/chromatin/Ricci/leaf_peaks.bed \
#   N_PERMUTATIONS=5000 N_CORES=8 \
#     Rscript scripts/local/analysis/20.Fine_Mapped_TF_Paper_Class_Enrichment_ATAC_Matched.R
#
# Output:
#   results/summary/fine_mapped_TF_paper_class_enrichment_leaf_ATAC_matched/
#     paper_class_overlap_matched_permutation.tsv
#     paper_class_overlap_matched_permutation_null_counts.rds
#     leaf_ATAC_status_matching_QC.tsv

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
  file.path(project_root, "results/summary/fine_mapped_TF_paper_class_enrichment_leaf_ATAC_matched")
)
set_default_env("EXPECTED_N_TFS", "8")
set_default_env("FEATURE_OUTPUT_STEM", "paper_class")
set_default_env("FEATURE_RESULT_COLUMN", "paper_class")
set_default_env(
  "FEATURE_LABEL_MAP_FILE",
  file.path(Sys.getenv("TF_PEAKS_DIR"), "paper_class_union_QC.tsv")
)

source(file.path(dirname(script_path), "19.Fine_Mapped_TF_ChIP_Enrichment_ATAC_Matched.R"), local = FALSE)
