### Script: analysis/1.Transcriptome_diversity.R
### Purpose: Normalize RNA-seq counts, remove flowering-time effects, compute per-gene expression variability (CV), fit LMMs across groups/eras, compute R2 for group/era effects, and run PCA on residual expression.
### Inputs:
###  - data/RNA/work/RNA_RawCounts_AllGenes.csv
###  - data/meta/Lines_meta.csv
###  - data/meta/BLUP2020.csv
### Outputs:
###  - results/CV/All_CV.csv, results/CV/Era_CV.csv
###  - results/R2/All_R2.csv, results/R2/Era_R2.csv
###  - results/PCA/PCA500.csv, results/PCA/PCA1000.csv, results/PCA/PCAAll.csv
###  - data/RNA/work/Filtered_Noralised_Counts.csv
### Usage: Rscript analysis/1.Transcriptome_diversity.R
### Dependencies: DESeq2, tidyverse, furrr, lme4, broom.mixed, FactoMineR
### Notes: Assumes sample/group columns as in inputs; count filter keeps genes with >=5 counts in >=160 samples. A variable 'residuals.mat.perm' is referenced but not defined before use in permutation section.

library(tidyverse)
library(furrr)
library(data.table)
library(FactoMineR)

### Load data
rna <- fread("data/processed/rna/counts_NE2020_corrected_names_deseq2_vst_filtered.csv.gz")
meta <- read.csv("data/raw/meta/Lines_meta.csv")
blup <- read.csv("data/raw/meta/BLUP2020.csv")

###
blup$Taxa[274] <- "ILLHy"

### Fill missing flowering time data with mean
blup$DaysToPollen[is.na(blup$DaysToPollen)] <- mean(blup$DaysToPollen, na.rm=T)

###
colnames(rna)[1] <- "Taxa"
norm.mat <- merge(meta, rna, by="Taxa")
norm.mat <- merge(blup[,c(1,2)], norm.mat, by.x="Taxa", by.y="Taxa")
rm(blup)

### Fit lm to normalised counts to remove flowering time effect 
residuals.mat <- norm.mat %>%
  select(Taxa, Group, DaysToPollen, starts_with("Zm")) %>%
  mutate(
    across(
      starts_with("Zm"),
      ~ residuals(lm(.x ~ DaysToPollen))
    )
  )


### Calculate coeffcient of variation for each gene in each group 
cv <- function(x){var(x, na.rm = T)/mean(x, na.rm=T)}

cv.exp <- residuals.mat %>% 
  pivot_longer(4:ncol(residuals.mat)) %>%
  group_by(Group, name) %>%
  summarise(cv = cv(abs(value)))

### PCA
genVar <- residuals.mat %>%
  select(starts_with("Zm")) %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "variance") %>%
  arrange(desc(variance))

resPCA500 <- PCA(residuals.mat[,genVar$Gene[1:500]], graph = F)
resPCA1000 <- PCA(residuals.mat[,genVar$Gene[1:1000]], graph = F)
resPCAAll <- PCA(residuals.mat[,4:ncol(residuals.mat)], graph = F)

resPCA500 <- resPCA500$ind$coord
resPCA1000 <- resPCA1000$ind$coord
resPCAAll <- resPCAAll$ind$coord

resPCA500 <- cbind.data.frame(residuals.mat[,1:3], resPCA500)
resPCA1000 <- cbind.data.frame(residuals.mat[,1:3], resPCA1000)
resPCAAll <- cbind.data.frame(residuals.mat[,1:3], resPCAAll)

### Calculate R2 of group effect for each gene
r2.df <- residuals.mat %>%
  select(Group, starts_with("Zm")) %>%
  summarise(
    future_map_dfc(
      select(., starts_with("Zm")),
      ~ {
        model <- lm(.x ~ residuals.mat$Group)
        summary(model)$r.squared  # Extract R-squared value
      },
      .names = "{col}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "gene",
    values_to = "r2"
  )

### Permutation 
set.seed(44)

res.r2.perm <- c()

residuals.mat.perm <- residuals.mat

for (i in 1:1000){

residuals.mat.perm$Group <- sample(residuals.mat.perm$Group)

r2.p <- residuals.mat.perm %>%
  select(Group, starts_with("Zm")) %>%
  summarise(
    future_map_dfc(
      select(., starts_with("Zm")),
      ~ {
        model <- lm(.x ~ residuals.mat.perm$Group)
        summary(model)$r.squared  # Extract R-squared value
      },
      .names = "{col}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "gene",
    values_to = "r2"
  )

r2.p$per <- i

res.r2.perm <- rbind.data.frame(res.r2.perm, r2.p)
print(i)
}

### Save results
write.csv(cv.exp, "results/full/GeneExpressionVariance/All_CV.csv", row.names = F)

write.csv(r2.df, "results/full/GeneExpressionVariance/All_R2.csv", row.names = F)

write.csv(resPCA500, "results/full/GeneExpressionVariance/PCA500.csv", row.names = F)
write.csv(resPCA1000, "results/full/GeneExpressionVariance/PCA1000.csv", row.names = F)
write.csv(resPCAAll, "results/full/GeneExpressionVariance/PCAAll.csv", row.names = F)
