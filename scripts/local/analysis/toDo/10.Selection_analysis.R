### Script: analysis/10.Selection_analysis.R
### Purpose: Analyze selection metrics (e.g., aFC, PIP, GERP) and
###          integrate with homologs and expression summaries.
### Inputs: results from AF/PIP pipelines, GERP, homolog mapping, RNA summaries
### Outputs: results/selection/* (selection metrics and merged outputs)
### Usage: Rscript analysis/10.Selection_analysis.R
### Dependencies: tidyverse, data.table, GenomicFeatures (and imported deps)

library(data.table)
library(GenomicFeatures)
library(tidyverse)


### Load data
# GERP scores
gerp <- fread("data/gerp/gerp_genes.bed")

# Gene information
txdb <- txdbmaker::makeTxDbFromGFF("../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz")
txdb <- keepStandardChromosomes(txdb)
gen <- genes(txdb)
seqlevelsStyle(gen) <- "NCBI"


### Load Allelic fold change results
af <- fread("results/AF/Susie2_AF.txt")

### Load SusieR results - number of credible sets
pip <- read.csv("results/SusieR2/PIP_Max.csv")


### Maize v4 v5 conversion 
v4v5 <- read.delim("../BigData/Maize/MetaData/B73v5_to_B73v4.tsv", header = F)
colnames(v4v5) <- c("v5", "v4")

### Homologs metadata
homologous <- fread("data/meta/sorghum3_intell_plusteff.csv")
homologous <- homologous[,c(1,4,5)]


# Select positive gerp sites 
gerp <- gerp %>%
  filter(V5>0)

gerp.db <- GenomicRanges::makeGRangesFromDataFrame(
  gerp,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3",
  keep.extra.columns = TRUE
)
seqlevelsStyle(gerp.db) <- "NCBI"

# Merge GERP with genes
gerp.gen.db <- GenomicRanges::mergeByOverlaps(gen, gerp.db)

# Gene and gerp data frame
data.gerp <- data.frame(gen=gerp.gen.db$gene_id, gerp=gerp.gen.db$V5)

# Create data frame with gene-level GERP summary
gerp.summary <- data.gerp %>%
  group_by(gen) %>%
  summarise(gen.mean=mean(gerp), n2=sum(gerp>2), n4=sum(gerp>4), gerp.sum=sum(gerp)) %>%
  dplyr::rename(gene=gen)

# Selection statistics
pip.res <- pip %>%
  group_by(gene) %>%
  filter(cs==max(cs)) %>%
  dplyr::select(gene, cs)

af.res <- af %>%
  dplyr::select(gene_id, log2_aFC) %>%
  dplyr::rename(gene=gene_id)

# Merge with GERP summaries and classify constraint
pip.res <- merge(pip.res, gerp.summary, by="gene")
pip.res$constrained <- ifelse(pip.res$gen.mean>2, "Top10", "Bottom90")

# Merge with GERP summaries and classify constraint
af.res <- merge(af.res, gerp.summary, by="gene")
af.res$constrained <- ifelse(af.res$gen.mean>2, "Top10", "Bottom90")


#
m12 <- homologous[,2:3]
m12$index <- 1:nrow(m12)

m1 <- m12 %>%
  dplyr::select(maize1_v4, index) %>%
  filter(maize1_v4 != "No Gene") %>%
  as.data.frame()
colnames(m1)[1] <- "v4"

m2 <- m12 %>%
  dplyr::select(maize2_v4, index) %>%
  filter(maize2_v4 != "No Gene") %>%
  as.data.frame()
colnames(m2)[1] <- "v4"

m1 <- left_join(m1, v4v5, by = "v4")

m1 <- m1[!duplicated(m1$index),]
m1 <- na.omit(m1)

m2 <- left_join(m2, v4v5, by = "v4")

m2 <- m2[!duplicated(m2$index),]
m2 <- na.omit(m2)

homolog <- merge(m1[,c(2,3)], m2[,c(2,3)], by="index")
colnames(homolog)[2:3] <- c("M1", "M2")

# RNA expression summary (mean per gene)
rna <- fread("data/RNA/work/Filtered_Noralised_Counts.csv")

rna.summary <- rna %>%
  dplyr::select(starts_with("Zm")) %>%
  pivot_longer(cols = starts_with("Zm")) %>%
  group_by(name) %>%
  summarise(meanExp=mean(value))
colnames(rna.summary)[1] <- "gene"

# Determine dominant/recessive expression within homologous pairs
dominant <- homolog %>% pivot_longer(2:3) %>%
  dplyr::rename(gene=value) %>%
  left_join(y = rna.summary, by="gene") %>%
  group_by(index) %>%
  mutate(highest = ifelse(meanExp == max(meanExp), "Dominant", "Recesive"))
dominant <- na.omit(dominant)


z <- plyr::join(dominant, gerp.summary)

z.af <- merge(z, af.res, by="gene")
z.pip <- merge(z, pip.res, by="gene")

# Write outputs
write.csv(af.res, "results/selection/AF_Selection.csv", row.names = F)
write.csv(pip.res, "results/selection/PIP_Selection.csv", row.names = F)

write.csv(z.af, "results/selection/WGS_AF.csv", row.names = F)
write.csv(z.pip, "results/selection/WGS_PIP.csv", row.names = F)
