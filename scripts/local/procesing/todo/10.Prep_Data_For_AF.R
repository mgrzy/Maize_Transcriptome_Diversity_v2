###
library(data.table)
library(dtplyr)
###

taxa <- fread("data/geno/vcf/taxa.txt")

gen.exp <- fread("data/RNA/work/Filtered_Noralised_Counts.csv")

taxa.names <- c("Name", gen.exp$Taxa)

gene.names <- colnames(gen.exp)[6:ncol(gen.exp)]

eqtl <- fread("results/SusieR2/PIP_Max.csv")

eqtl <- eqtl %>% 
  distinct(gene, cs, .keep_all = T) %>%
  select(gene, SNP) %>%
  rename_with(~c("gene_id", "variant_id"))

g <- gen.exp %>% 
  select(!c(Taxa,DaysToPollen, Taxa2, Group, Era)) %>%
  transpose() %>%
  mutate(Name=gene.names, .before = 1) %>%
  rename_with(~c(taxa.names)) %>%
  filter(Name %in% eqtl$gene_id)


fwrite(g, "data/RNA/work/Matrix_AF.csv.gz")

###

fwrite(eqtl, "results/SusieR2/eQTL_List.tsv", sep = "\t")
### Script: processing/9.Prep_Data_For_AF.R
### Purpose: Prepare inputs for allele frequency (AF) analysis of SusieR-lead SNPs/eQTLs across populations/eras.
### Inputs: lead SNP tables, population definitions, genotype data indices
### Outputs: results/AF/* (AF-ready inputs or merged AF tables)
### Usage: Rscript processing/9.Prep_Data_For_AF.R
### Dependencies: tidyverse, data.table (if used)
