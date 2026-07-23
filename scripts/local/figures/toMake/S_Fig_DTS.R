### #####################
### Genetic diversity ###
#########################

library(tidyverse)
library(scales)
library(ggpmisc)

# Match main figure styling
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")

### Load and prep data

# PCA
pc <- read.table("results/PCA/PCA_Geno.eigenvec")
pc <- pc[,-1]
colnames(pc) <- c("Taxa", "PC1", "PC2", "PC3", "PC4", "PC5")

# Meta data 
meta <- read.csv("data/meta/WiDiv_meta.csv")
pc <- merge(meta, pc, by="Taxa")

#
blup <- read.csv("data/meta/BLUP2020.csv")

A <- merge(pc, blup, by="Taxa")
colnames(A)[3] <- "Subpopulation"
###
p <- ggplot(A, aes(Subpopulation, DaysToSilk, fill=Subpopulation)) +
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  stat_multcomp(label.type = "letters",
                adj.method.tag = F,
                geom = "label") + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = myPal) + 
  ylab("Days to silking")

# Export with main figure parameters (PNG + PDF, cm units)
ggsave(
  filename = "results/figures/Supp/FigS_DTS.png",
  plot = p,
  width = 8.75,
  height = 6,
  dpi = 300,
  units = "cm"
)
ggsave(
  filename = "results/figures/Supp/FigS_DTS.pdf",
  plot = p,
  width = 8.75,
  height = 6,
  units = "cm"
)
### Script: figures/S_Fig_DTS.R
### Purpose: Supplemental figure for DTS (days to pollen/flowering time) related analyses.
### Inputs: BLUP/phenotype summaries and related outputs
### Outputs: results/figures/S_Fig_DTS.png
### Usage: Rscript figures/S_Fig_DTS.R
### Dependencies: tidyverse, ggplot2
