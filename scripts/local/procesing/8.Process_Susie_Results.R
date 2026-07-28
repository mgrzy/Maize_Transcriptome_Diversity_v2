###
###
###
library(data.table)
library(tidyverse)

### Load data
pip <- fread("results/full/SusieR/PIP_all_batches.tsv.gz")

eGene <- read.table("results/eQTL/eGenes.txt", header = F)

genePos <- read.table("data/OSCA/rawData/WiDiv631_OPI.txt")
genePos <- genePos[,c(1:3, 5)]
colnames(genePos) <- c("chr", "gene", "TSS", "dir")

### Number of credible sets per gene
cs <- pip %>%
  group_by(gene) %>%
  summarise(cs=max(cs)) %>%
  count(cs) %>%
  add_row(.before = 1, cs=0, n=1123)

### Number of variants per credible sets
cs.var <- pip %>%
  group_by(cs) %>%
  count(gene, name = "nVar") %>%
  ungroup() %>%
  count(nVar, name = "nCS")

###
cs.gene <- pip %>%
  group_by(gene) %>%
  summarise(cs=max(cs))

###
pip$POS <- as.numeric(str_split_fixed(pip$SNP, "_", 2)[,2])
pip.pos <- merge(pip, genePos, by="gene")

pip.pos$dist <- ifelse(pip.pos$dir=="+", pip.pos$POS - pip.pos$TSS, pip.pos$TSS - pip.pos$POS)

pip.pos <- pip.pos[,c(1,5:8, 10)]

###
pip.max <- pip %>%
  group_by(gene, cs) %>%
  filter(variable_prob == max(variable_prob)) %>%
  arrange(gene, SNP, variable_prob, cs)

### Save results
write.csv(pip, "results/SusieR2/PIP.csv", row.names = F)
write.csv(cs, "results/SusieR2/CS.csv", row.names = F)
write.csv(cs.var, "results/SusieR2/CS_VAR.csv", row.names = F)
write.csv(pip.pos, "results/SusieR2/PIP_dist.csv", row.names = F)
write.csv(pip.max, "results/SusieR2/PIP_Max.csv", row.names = F)
write.csv(cs.gene, "results/SusieR2/CS_Gene.csv", row.names = F)
### Script: processing/8.Process_Susie_Results2.R
### Purpose: Alternate/extended post-processing of SusieR outputs; generates variant-level and eGene-level summaries.
### Inputs: ../results/SusieR/* (per-batch outputs)
### Outputs: results/SusieR2/* (supplementary CS/PIP summaries)
### Usage: Rscript processing/8.Process_Susie_Results2.R
### Dependencies: tidyverse, data.table (if used)
