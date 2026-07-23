library(data.table)
library(Hmisc)
library(IHW)
###
dt <- fread(cmd = paste("awk -F'\t' '{print $3, $9, $14}' results/eQTL/AlleQTL"), stringsAsFactors = T)

# Compute ranks
# dt[, rank := rank(p, ties.method = "first")]
# 
# # Compute FDR using BH method manually
# dt[, fdr := p * .N / rank]
# dt[, fdr := pmin(fdr, 1)]  
# 
# dt.sig <- dt[fdr<0.001]

#3.54368e-05
# n=1040807904

# Compute distance between SNPs (BP) and TSS (Probe_bp)
dt[, distance := BP - Probe_bp]
gc()

#number of bins
nbins <- 20

# Discretization of covariate (distance of SNPs to TSSs) with cut2 function from Hmisc
dt[, bin := cut2(distance, g = nbins)]
gc()

# Calulte number of tests per group 
m.groups <- table(dt$bin)

# Keep only low p-values
dt.small <- dt[p < 0.01] #105,442,590

rm(dt)
gc()

# Run IHW
res <- ihw(p ~ bin, dt.small, alpha=0.01, m_groups = m.groups)
### Script: analysis/15.Calulate_FDR_eQTL.R
### Purpose: Calculate FDR for eQTL signals using permutation/null distributions and consolidate results.
### Inputs: eQTL summary outputs and null/permuted results
### Outputs: results/FDR/* (FDR-adjusted results)
### Usage: Rscript analysis/15.Calulate_FDR_eQTL.R
### Dependencies: tidyverse and any necessary stats packages
