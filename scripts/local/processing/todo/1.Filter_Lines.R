### Script: processing/1.Filter_Lines.R
### Purpose: Filter and harmonize lines/samples across metadata, RNA counts, and genotype availability; create downstream meta files and PLINK keep lists.
### Inputs:
###  - data/meta/WiDiv_meta.csv
###  - data/RNA/raw/counts_NE2020_corrected_names.csv
###  - ../BigData/Maize/v5/Genotypes/All/plink/All_maf001.fam
###  - data/meta/new1814lines.txt
### Outputs:
###  - data/meta/Lines_meta.csv, data/RNA/work/RNA_RawCounts_AllGenes.csv
###  - data/meta/LinesPlink.txt, data/meta/LinesPop.txt, data/meta/Pops.txt
###  - data/meta/LinesEra.txt, data/meta/Eras.txt
### Usage: Rscript processing/1.Filter_Lines.R
### Dependencies: tidyverse

#################################
### Filter lines for analysis ###
#################################
library(tidyverse)

# Load data
meta <- read_csv("data/meta/WiDiv_meta.csv")
rna <- read.csv("data/RNA/raw/counts_NE2020_corrected_names.csv")
fam <- read.table("../BigData/Maize/v5/Genotypes/All/plink/All_maf001.fam")

# Change line names which starts with number to correct one
colnames(rna)[2:17] <- gsub("X", "", colnames(rna)[2:17])

# Extract gene names
gene <- rna$GeneID

# Transpose data frame
rna <- t(rna[,-c(1)])

# Create data frame with lines in rows and gene in columns
rna <- data.frame(Taxa2=rownames(rna), rna)
colnames(rna)[2:ncol(rna)] <- gene
df <- merge(meta, rna, by="Taxa2")

# Select lines which belong to one of the six group
df <- df %>% 
  filter(Group %in% c("SS", "NSS", "IDT", "Tropical", "Mixed", "Broad origin-public"), Taxa %in% fam$V1)

meta <- meta[meta$Taxa %in% fam$V1,]

# Change Broad origin-public group name and Public sector to public
df$Group[df$Group=="Broad origin-public"] <- "Broad_origin-public"
df$Era[which(df$Era=="Public Sector")] <- "Public"
# Create file with lines in Plink format
plink <- data.frame(V1=df$Taxa, V2=df$Taxa, V3=0, V4=0, V5=0, V6=-9)

# Create files with line names for pop gene analisi
meta <- df[,2:4]
meta.era <- na.omit(meta)

meta.era <- meta.era %>%
  filter(Group %in% c("SS", "NSS", "Mixed"))

meta.era$Era <- paste0(meta.era$Group, "_", meta.era$Era)

###
new <- read.table("data/meta/new1814lines.txt")

meta <- meta[meta$Taxa %in% new$V1,]
meta.era <- meta.era[meta.era$Taxa %in% meta$Taxa,]
### Save data
write.csv(df[,c(1:4)], "data/meta/Lines_meta.csv", row.names = F)
write.csv(df[,c(1, 5:ncol(df))], "data/RNA/work/RNA_RawCounts_AllGenes.csv", row.names = F)
write.table(plink, "data/meta/LinesPlink.txt", row.names = F, col.names = F, quote = F)

write.table(meta[,1], "data/meta/LinesPop.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(meta[,1:2], "data/meta/Pops.txt", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(meta.era[,1], "data/meta/LinesEra.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(meta.era[,c(1,3)], "data/meta/Eras.txt", row.names = F, col.names = F, quote = F, sep = "\t")


