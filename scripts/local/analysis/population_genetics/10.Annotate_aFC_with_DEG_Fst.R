### Script: analysis/10.Annotate_aFC_with_DEG_Fst.R
### Purpose: Annotate fine-mapped gene-variant pairs with Limma DEG status,
###          Limma adjusted P values, and GCTA Fst for every shared contrast.
### Inputs:
###   - results/summary/AF/aFCn_MAF_se_filtered.tsv.gz
###   - results/summary/DEG/DEG_limma_summary.csv
###   - results/full/Fst/*.fst
### Output:
###   - results/summary/AF/aFCn_MAF_se_filtered_DEG_Fst.tsv.gz
### Usage: Rscript scripts/local/analysis/10.Annotate_aFC_with_DEG_Fst.R
### Dependencies: data.table

suppressPackageStartupMessages(library(data.table))

afc_path <- "results/summary/AF/aFCn_MAF_se_filtered.tsv.gz"
deg_path <- "results/summary/DEG/DEG_limma_summary.csv"
fst_dir <- "results/full/Fst"
output_path <- "results/summary/AF/aFCn_MAF_se_filtered_DEG_Fst.tsv.gz"

afc <- fread(afc_path)
required_afc_columns <- c("gene_id", "variant_id")
if (!all(required_afc_columns %in% names(afc))) {
  stop("aFC input must contain: ", paste(required_afc_columns, collapse = ", "))
}
afc[, input_order__ := .I]

deg <- fread(deg_path)
deg_status_columns <- grep("_DEG$", names(deg), value = TRUE)
deg_padj_columns <- grep("_adj\\.P\\.Val$", names(deg), value = TRUE)
if (!"gene" %in% names(deg) || length(deg_status_columns) == 0L ||
    length(deg_padj_columns) == 0L) {
  stop(
    "DEG summary must contain a gene column plus *_DEG and *_adj.P.Val columns."
  )
}
# Keep the adjusted P value for each contrast together with its corresponding
# binary DEG annotation (DEG = adjusted P value < 0.05).
deg <- deg[, c("gene", deg_status_columns, deg_padj_columns), with = FALSE]
setkey(deg, gene)

# Retain every aFC gene-variant pair, including genes not tested by a contrast.
afc <- merge(
  afc,
  deg,
  by.x = "gene_id",
  by.y = "gene",
  all.x = TRUE,
  sort = FALSE
)

fst_files <- list.files(fst_dir, pattern = "\\.fst$", full.names = TRUE)
if (length(fst_files) == 0L) {
  stop("No GCTA .fst files found in ", fst_dir)
}

for (fst_file in sort(fst_files)) {
  contrast_name <- sub("\\.fst$", "", basename(fst_file))
  fst_column <- paste0(contrast_name, "_Fst")

  # Read only the columns needed for the variant-level annotation. Processing
  # one full-genome Fst result at a time keeps memory use bounded.
  fst <- fread(fst_file, select = c("SNP", "Fst"), showProgress = FALSE)
  setnames(fst, c("SNP", "Fst"), c("variant_id", fst_column))

  afc <- merge(
    afc,
    fst,
    by = "variant_id",
    all.x = TRUE,
    sort = FALSE
  )
  rm(fst)
  gc(verbose = FALSE)
  message("Annotated ", contrast_name)
}

setorder(afc, input_order__)
afc[, input_order__ := NULL]
fwrite(afc, output_path, sep = "\t", compress = "gzip")
message("Wrote ", nrow(afc), " annotated gene-variant pairs to ", output_path)
