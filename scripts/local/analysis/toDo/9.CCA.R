### Script: analysis/9.CCA.R
### Purpose: Compute canonical correlation analysis (CCA) between expression PCs and genotype PCs; outputs per-sample canonical scores.
### Inputs:
###  - results/PCA/PCAAll.csv (expression PCs)
###  - results/PCA/PCA_Geno.eigenvec (genotype PCs)
### Outputs:
###  - results/CCA/CCA.csv
### Usage: Rscript analysis/9.CCA.R
### Dependencies: CCA

library(CCA)

pcaAll <- read.csv("results/PCA/PCAAll.csv")

pca.gen <- read.table("results/PCA/PCA_Geno.eigenvec")
colnames(pca.gen)[1] <- "Taxa"
pca.gen <- pca.gen[,-c(2)]
colnames(pca.gen)[2:6] <- c("PC1", "PC2", "PC3", "PC4", "PC5")

X <- merge(pca.gen, pcaAll, by = "Taxa")

res <- cc(X[,9:ncol(X)], X[,2:6])

res <- as.data.frame(res$scores$xscores)

res <- cbind.data.frame(X[,c(1,7,8)], res)
colnames(res)[4:8] <- c("CCA1", "CCA2", "CCA3", "CCA4", "CCA5")

write.csv(res, "results/CCA/CCA.csv", row.names = F)
