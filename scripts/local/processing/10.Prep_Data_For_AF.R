###
library(data.table)
###

taxa <- fread("data/processed/geno/vcf/Div639.nosex")

gen.exp <- fread("data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz")

taxa.names <- as.character(taxa[[2L]])

if (!identical(as.character(gen.exp$Line), taxa.names)) {
  stop("Expression samples do not match the Div639 VCF sample order.")
}

eqtl <- fread("results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz")

eGene <- fread("results/summary/osca/eGene.txt", col.names = "gene")

eqtl <- unique(eqtl[gene %in% eGene$gene], by = c("gene", "cs"))[
  , .(gene_id = gene, variant_id = SNP)
]

# transpose() returns a list; explicitly convert it to a data.table before
# adding gene and sample names.  The expression matrix has one non-expression
# column (Line), so retain every other column as a gene.
expression_values <- gen.exp[, setdiff(names(gen.exp), "Line"), with = FALSE]
g <- as.data.table(transpose(expression_values))
g[, Name := names(expression_values)]
setcolorder(g, "Name")
setnames(g, c("Name", taxa.names))
g <- g[Name %in% eqtl$gene_id]


dir.create("data/processed/AF", recursive = TRUE, showWarnings = FALSE)
fwrite(g, "data/processed/AF/Matrix_AF.csv.gz")

###

fwrite(eqtl, "data/processed/AF/eQTL_List.tsv", sep = "\t")
### Script: processing/9.Prep_Data_For_AF.R
### Purpose: Prepare inputs for allele frequency (AF) analysis of SusieR-lead SNPs/eQTLs across populations/eras.
### Inputs: lead SNP tables, population definitions, genotype data indices
### Outputs: results/AF/* (AF-ready inputs or merged AF tables)
### Usage: Rscript processing/9.Prep_Data_For_AF.R
### Dependencies: tidyverse, data.table (if used)
