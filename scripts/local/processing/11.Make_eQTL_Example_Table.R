#!/usr/bin/env Rscript
# Create an eQTL example dataset for Zm00001eb146510.
#
# Outputs:
#   results/eqtl_example/Zm00001eb146510_sample_expression_genotypes.tsv.gz
#   results/eqtl_example/Zm00001eb146510_fine_mapped_variants.tsv

library(bigmemory)
library(data.table)

TARGET_GENE <- "Zm00001eb146510"
PIP_PATH <- "results/full/SusieR/PIP_all_batches.tsv.gz"
GENO_DESC_PATH <- "data/processed/geno/rMVP/Div639.geno.desc"
GENO_IND_PATH <- "data/processed/geno/rMVP/Div639.geno.ind"
GENO_MAP_PATH <- "data/processed/geno/rMVP/Div639.geno.map"
EXPRESSION_PATH <- "data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz"
OUTPUT_DIR <- "results/eqtl_example"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# All rows in PIP_all_batches are members of SuSiE credible sets. Retain every
# distinct SNP for the target gene, rather than only each credible set's lead.
fine_mapped <- fread(PIP_PATH)[gene == TARGET_GENE]
if (!nrow(fine_mapped)) {
  stop("No fine-mapped variants found for ", TARGET_GENE)
}
setorder(fine_mapped, cs, -variable_prob)
fine_mapped <- unique(fine_mapped, by = "SNP")

# The rMVP map contains ~29 million rows. Stream it with awk to obtain only
# the selected SNPs and their 1-based big.matrix column indices.
snp_ids_path <- tempfile("eqtl_example_snps_", fileext = ".txt")
fwrite(fine_mapped[, .(SNP)], snp_ids_path, col.names = FALSE)
map_cmd <- paste(
  "awk 'BEGIN { OFS = \"\\t\" } NR == FNR { wanted[$1] = 1; next } FNR > 1 && ($1 in wanted) { print NR - 1, $1, $2, $3, $4, $5 }'",
  shQuote(snp_ids_path),
  shQuote(GENO_MAP_PATH)
)
map_subset <- fread(
  cmd = map_cmd,
  col.names = c("matrix_column", "SNP", "CHROM", "POS", "REF", "ALT")
)
unlink(snp_ids_path)

if (nrow(map_subset) != nrow(fine_mapped)) {
  missing_snps <- setdiff(fine_mapped$SNP, map_subset$SNP)
  stop("Fine-mapped SNPs absent from the rMVP map: ", paste(missing_snps, collapse = ", "))
}

variant_metadata <- map_subset[match(fine_mapped$SNP, SNP)]
variant_metadata[, `:=`(
  gene = fine_mapped$gene,
  cs = fine_mapped$cs,
  variable = fine_mapped$variable,
  variable_prob = fine_mapped$variable_prob,
  source_batch_file = fine_mapped$source_batch_file
)]
setcolorder(
  variant_metadata,
  c("gene", "cs", "SNP", "variable", "variable_prob", "CHROM", "POS", "REF", "ALT", "matrix_column", "source_batch_file")
)

individuals <- fread(GENO_IND_PATH, header = FALSE, col.names = "Line")
individuals[, Line := as.character(Line)]
expression <- fread(EXPRESSION_PATH, select = c("Line", TARGET_GENE))
expression[, Line := as.character(Line)]
expression_index <- match(individuals$Line, expression$Line)
if (anyNA(expression_index)) {
  stop("Some rMVP lines are absent from the expression matrix.")
}

geno <- attach.big.matrix(GENO_DESC_PATH)
if (nrow(geno) != nrow(individuals)) {
  stop("The rMVP genotype matrix and .ind file have different sample counts.")
}

genotypes <- geno[, variant_metadata$matrix_column, drop = FALSE]
storage.mode(genotypes) <- "double"
setnames(
  genotypes <- as.data.table(genotypes),
  paste0("genotype_", variant_metadata$SNP)
)

sample_table <- cbind(
  individuals,
  expression_log2_cpm = as.numeric(expression[[TARGET_GENE]][expression_index]),
  genotypes
)

sample_output <- file.path(
  OUTPUT_DIR,
  paste0(TARGET_GENE, "_sample_expression_genotypes.tsv.gz")
)
variant_output <- file.path(
  OUTPUT_DIR,
  paste0(TARGET_GENE, "_fine_mapped_variants.tsv")
)
fwrite(sample_table, sample_output, sep = "\t")
fwrite(variant_metadata, variant_output, sep = "\t")

message("Wrote ", nrow(sample_table), " samples and ", nrow(variant_metadata), " fine-mapped variants to ", OUTPUT_DIR)
