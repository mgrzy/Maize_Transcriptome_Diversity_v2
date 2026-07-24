###
###
###
library(data.table)
library(bigsnpr)


### Find eGenes
eQTL <- fread("results/eQTL/testTable")

eGene <- unique(eQTL$Gene)

write.table(eGene, "results/eQTL/eGenes.txt", row.names = F, col.names = F, quote = F)
###
snp <- snp_readBed("data/geno/plink/Div639.bed")
### Script: processing/8.Prep_Data-For_SusiR.R
### Purpose: Prepare inputs for SusieR fine-mapping (e.g., expression, covariates, and variant/gene ranges).
### Inputs: normalized expression, covariates, genotype metadata
### Outputs: ../SusieR_Data/* (files consumed by HPCcode/SusieR.R)
### Usage: Rscript processing/8.Prep_Data-For_SusiR.R
### Dependencies: tidyverse and any file I/O helpers used
