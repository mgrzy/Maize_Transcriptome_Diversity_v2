### Script: analysis/3.Find_Covariate_for_eQTL.R
### Purpose: Identify covariates for eQTL mapping by selecting RNA PCs (via elbow and BE) and filtering known covariates correlated with RNA PCs.
### Inputs:
###  - data/RNA/work/Filtered_Noralised_Counts.csv
###  - results/PCA/PCA_Geno.eigenvec (genotype PCs)
### Outputs:
###  - data/meta/CovariateForEQTL.csv
### Usage: Rscript analysis/3.Find_Covariate_for_eQTL.R
### Dependencies: tidyverse, PCAForQTL

#########################################
### Find covariates for eQTL analysis ###
#########################################
library(tidyverse)
library(PCAForQTL)

base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))

### Load data
dat <- data.table::fread("data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz", data.table = F)

### 
genoPC <- bigmemory::attach.big.matrix("data/processed/geno/rMVP/Dib639_MAF005.pc.desc")[]
colnames(genoPC) <- c("PC_G1", "PC_G2", "PC_G3", "PC_G4", "PC_G5")

blup <- data.table::fread("data/raw/meta/BLUP2020.csv")

### 
know.cov <- data.frame(Taxa=dat$Line, genoPC)
rownames(know.cov) <- NULL

know.cov <- plyr::join(know.cov, blup[,c("Taxa", "DaysToPollen")], by="Taxa")

### Select column with expression
dat.exp <- dat %>% 
  select(starts_with("Zm"))

### Run PCA
resPCA <- prcomp(dat.exp, center = T, scale. = T)

### Run Elbow
resElbow <- runElbow(prcompResult = resPCA)

### Run BE
set.seed(44)
resBE <- runBE(dat.exp, B = 100, alpha = 0.05, mc.cores = 30)

### Selected number of PCs
K_elbow <- resElbow
K_BE <- resBE$numOfPCsChosen

### Scree plot
makeScreePlot(resPCA, labels = c("Elbow", "BE"), values = c(K_elbow, K_BE), maxNumOfPCsToPlot = 150)

ggsave(
  filename = "results/figures/supp/Fig_PCA_Var_exp.png",
  width = 12,
  height = 6,
)
