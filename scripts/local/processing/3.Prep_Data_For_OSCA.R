### Script: processing/3.Prep_Data_For_OSCA.R
### Purpose: Prepare OSCA input files from normalized expression, sample list, and covariates.
###
### What this creates (OSCA formats):
###  - OPI (probe info): per-gene metadata (chr, start, strand, geneID)
###  - OII (individual info): sample IDs with placeholders required by OSCA
###  - Expression .dat: tab-delimited expression matrix (Taxa + gene columns)
###  - Covariates .txt: sample-aligned covariate table (PCs, etc.)
###
### Inputs:
###  - data/geno/plink/Div639.fam                 (sample list from PLINK; used to define Taxa order)
###  - data/RNA/work/Filtered_Noralised_Counts.csv (normalized expression; must contain a Taxa column)
###  - data/meta/CovariateForEQTL.csv              (covariates; must contain a Taxa column)
###  - ../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz (gene annotation)
###
### Outputs:
###  - data/OSCA/rawData/WiDiv631_OPI.txt  (probe info)
###  - data/OSCA/rawData/WiDiv631_OII.txt  (individual info)
###  - data/OSCA/rawData/WiDiv631.dat      (expression matrix)
###  - data/OSCA/WiDiv631_PC.txt           (covariates aligned to OII order)
###
### Usage:
###   Rscript processing/3.Prep_Data_For_OSCA.R
###
### Dependencies:
###   - tidyverse (for dplyr verbs and pipes)
###   - data.table (fread) and stringr (used with explicit namespaces)

##########################
### Prep data for OSCA ###
##########################

library(tidyverse)

# 1) Define the set and order of individuals (Taxa) from PLINK .fam
#    We rename the first column to Taxa for consistency across inputs.
fam <- read.table("data/processed/geno/plink/Div639.fam")
colnames(fam)[1] <- "Taxa"

# 2) Load normalized expression matrix
#    Expect a column named Taxa followed by expression columns (gene IDs).
rna <- data.table::fread("data/processed/rna/counts_NE2020_corrected_names_log2_cpm_filtered.csv.gz", data.table = F)
colnames(rna)[1] <- "Taxa"
# 3) Load covariates (e.g., PCs) — must include a Taxa column
PC <- read.csv("results/full/covariateForEqlt/CovariateForEQTL.csv")
PC <- data.frame(Taxa=fam$Taxa, PC)

# 4) Build OPI (OSCA Probe Information) from GFF3
#    - Keep only rows annotated as "gene"
#    - Extract required fields: chr, start, end, strand (or), geneID
#    - Normalize chromosome labels and gene IDs
gene <- data.table::fread("../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz")
gene <- gene[gene$V3=="gene",]
gene <- gene[,c(1,4,5,7,9)]
colnames(gene) <- c("chr", "start", "end", "or","geneID")
gene$chr <- as.numeric(gsub("chr", "", gene$chr))
gene$geneID <- gsub("ID=", "", gene$geneID)
gene$geneID <- stringr::str_split_fixed(gene$geneID, ";", 2)[,1]
gene <- na.omit(gene)

opi <- data.frame(chr=gene$chr,
                  probeID=gene$geneID,
                  pos=gene$start, 
                  geneID=gene$geneID,
                  or=gene$or)

# Keep only genes present in the expression matrix (OPI rows must match data columns)
opi <- opi %>% filter(probeID %in% colnames(rna))

# 5) Create OII (OSCA Individual Information)
#    Format: Taxa, a, b, z
#    - a/b are typically allele codes for genotypes; for expression they can be placeholders
#    - z is a placeholder for phenotype; often NA for OSCA preprocessing
oii <- data.frame(Taxa=fam$Taxa, a=0, b=0, z =NA)

# 6) Create covariate table aligned to the OII sample order
#    - left_join with OII$Taxa preserves the OII order
#    - add Taxa2 duplicate after Taxa (some OSCA utilities expect two ID columns)
#    - drop the 3rd column if the raw covariate file has an extra non-PC column
PC <- left_join(oii[,"Taxa", drop=FALSE], PC, by="Taxa")
#PC <- PC[,-c(1)]
colnames(PC)[2] <- "Taxa2"
PC <- PC[,c(-18:-32)]
#PC <- PC %>% mutate(Taxa2 = PC$Taxa, .after = Taxa)

# 7) Reduce the expression matrix to the OSCA-required columns
#    - keep Taxa and gene columns that are present in OPI
genes <- colnames(rna)[colnames(rna) %in% opi$geneID]

rna2 <- rna %>%
  select(Taxa, genes)

# mat.exp.norm <- apply(mat.exp[,-c(1)], 2, function(x) {
#   x <- x + runif(length(x))
#   bestNormalize::boxcox(x)$x.t
# }
# )
# mat.exp.norm <- data.frame(Taxa=rownames(mat.exp.norm), mat.exp.norm)

# 8) Write OSCA input files
#    Note: downstream script processing/6.Make_Data_For_OSCA.sh uses these as inputs to OSCA
write.table(oii, file = "data/processed/osca/WiDiv631_OII.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(opi, file = "data/processed/osca/WiDiv631_OPI.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(rna2, file = "data/processed/osca/WiDiv631.dat", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(PC, file = "data/processed/osca/WiDiv631_PC.txt", quote = F, row.names = F, col.names = T, sep = "\t")

