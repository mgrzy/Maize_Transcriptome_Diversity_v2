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

### Load data
dat <- read.csv("data/RNA/work/Filtered_Noralised_Counts.csv")

### 
genoPC <- read.table("results/PCA/PCA_Geno.eigenvec")
genoPC <- genoPC[,-2]
colnames(genoPC) <- c("Taxa", "PC_G1", "PC_G2", "PC_G3", "PC_G4", "PC_G5")

### 
know.cov <- left_join(dat[,c(1,2)], genoPC, by = "Taxa")
rownames(know.cov) <- NULL

### Select column with expression
dat.exp <- dat %>% 
  select(starts_with("Zm"))

### Run PCA
resPCA <- prcomp(dat.exp, center = T, scale. = T)

### Run Elbow
resElbow <- runElbow(prcompResult = resPCA)

### Run BE
set.seed(44)
resBE <- runBE(dat.exp, B = 20, alpha = 0.05, mc.cores = 30)

### Selected number of PCs
K_elbow <- resElbow
K_BE <- resBE$numOfPCsChosen

### Scree plot
makeScreePlot(resPCA, labels = c("Elbow", "BE"), values = c(K_elbow, K_BE), maxNumOfPCsToPlot = 150)

###  Select PCs based on BE
PCsTop <- resPCA$x[,1:K_BE]
PCsTop <- as.data.frame(PCsTop)

### Remove know covariats which correlates with PC based on RNA
knownCovariatesFiltered <- filterKnownCovariates(know.cov[,-c(1)], PCsTop, unadjustedR2_cutoff = 0.9)

### Create file with covariats for eQTL
PCsTop <- scale(PCsTop)
covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)

covariatesToUse <- covariatesToUse %>%
  mutate(Taxa = dat$Taxa, .before = DaysToPollen)

###
write.csv(covariatesToUse, "data/meta/CovariateForEQTL.csv", row.names = F)
