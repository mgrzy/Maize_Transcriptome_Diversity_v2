library(DESeq2)
library(mashr)

### Load data
rna <- read.csv("data/RNA/work/Era.csv")
meta <- read.csv("data/meta/Era_meta.csv")

###
colnames(rna)[1] <- "Taxa2"

rna <- merge(meta, rna, by="Taxa2")

meta <- rna[,c(1:5)]
rna <- rna[,-c(1:5)]

meta$GroupEra <- paste0(meta$Group, meta$Era)
### Normalised RNA counts 
dds <- DESeqDataSetFromMatrix(countData = round(t(rna)),
                              colData = meta,
                              design = ~ DaysToPollen + GroupEra)
  
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds) >= 10) >= 135
dds <- dds[keep,]
  
# run DEG analysis 
dds <- DESeq(dds, parallel = T)
  
# get raw estimates of effect sizes
res.Era.SS.raw <- results(dds, contrast=c("GroupEra","SSexPVP", "SSPublic"))
res.Era.NSS.raw <- results(dds, contrast=c("GroupEra","NSSexPVP", "NSSPublic"))

# DEG for Heterotic group
res.Era.Public.raw <- results(dds, contrast=c("GroupEra","SSPublic", "NSSPublic")) #3939
res.Era.exPVP.raw <- results(dds, contrast=c("GroupEra","SSexPVP", "NSSexPVP")) # 9615

# create mashr data object
beta <- cbind(res.Era.SS.raw$log2FoldChange, res.Era.NSS.raw$log2FoldChange)
colnames(beta) <- c("SS", "NSS")
rownames(beta) <- rownames(res.Era.SS.raw)
serror <- cbind(res.Era.SS.raw$lfcSE, res.Era.NSS.raw$lfcSE)
colnames(serror) <- c("SS", "NSS")
rownames(serror) <- rownames(res.Era.SS.raw)

data.mash <- mash_set_data(beta, serror)


dds2 <- DESeqDataSetFromMatrix(countData = round(t(rna)),
                              colData = meta,
                              design = ~ DaysToPollen + Group*Era)

dds2 <- DESeq(dds, parallel = T)

# set up the covariates matrix
U.c <- cov_canonical(data.mash)

# fit the model
m.c <- mash(data.mash, U.c)

# get the posterior mean effect sizes
beta.mash <- as.data.frame(m.c$result$PosteriorMean)
beta.mash$lfdr.SS <- m.c$result$lfsr[,1]
beta.mash$lfdr.NSS <- m.c$result$lfsr[,2]
beta.mash$n_sig <- get_n_significant_conditions(m.c)
beta.mash$SS_sd <- m.c$result$PosteriorSD[,1]
beta.mash$NSS_sd <- m.c$result$PosteriorSD[,2]
beta.mash <- beta.mash[order(beta.mash$n_sig, decreasing = F),]
beta.mash$gene <- rownames(beta.mash)

###
res.Era.SS.raw <- as.data.frame(res.Era.SS.raw)
res.Era.SS.raw$gene <- rownames(res.Era.SS.raw)
colnames(res.Era.SS.raw)[1:6] <- paste0(colnames(res.Era.SS.raw)[1:6], "_SS")

res.Era.NSS.raw <- as.data.frame(res.Era.NSS.raw)
res.Era.NSS.raw$gene <- rownames(res.Era.NSS.raw)
colnames(res.Era.NSS.raw)[1:6] <- paste0(colnames(res.Era.NSS.raw)[1:6], "_NSS")

beta.mash <- merge(beta.mash, res.Era.SS.raw, by="gene")

beta.mash <- merge(beta.mash, res.Era.NSS.raw, by="gene")

write.csv(beta.mash, "results/DEG/DEG_ERA.csv", row.names = F)
write.csv(res.Era.Public.raw, "results/DEG/DEG_Public_SS_NSS.csv")
write.csv(res.Era.exPVP.raw, "results/DEG/DEG_exPVP_SS_NSS.csv")
### Script: analysis/7.DEG.R
### Purpose: Differential expression analysis across specified contrasts/groups (details based on script contents) to support downstream figures.
### Inputs: RNA-seq normalized counts and metadata
### Outputs: results/DEG/* (DE results per contrast)
### Usage: Rscript analysis/7.DEG.R
### Dependencies: DESeq2, tidyverse, and any other packages loaded in the script
