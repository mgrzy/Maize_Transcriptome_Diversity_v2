### Script: figures/Fig3.R
### Purpose: Generate Figure 3 panels: credible set counts per eGene, relationship between number of variants and credible sets, and distribution of aFC for eQTNs.
### Inputs:
###  - results/SusieR2/CS.csv, results/SusieR2/CS_VAR.csv
###  - results/SusieR2/PIP_dist.csv, results/SusieR2/PIP_Max.csv
###  - results/AF/Susie2_AF.txt
### Outputs: results/figures/Fig3.png
### Usage: Rscript figures/Fig3.R
### Dependencies: tidyverse, scales, ggpmisc, data.table, patchwork
### #####################
### Genetic diversity ###
#########################

library(tidyverse)
library(scales)
library(ggpmisc)
library(data.table)
library(patchwork)

### ggplot2 settings (match Fig1 font sizing)
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")

cs <- read.csv("results/SusieR2/CS.csv")
cs.var <- read.csv("results/SusieR2/CS_VAR.csv")
snp.dist <- read.csv("results/SusieR2/PIP_dist.csv")
snp.lead <- read.csv("results/SusieR2/PIP_Max.csv")

snp.dist$Type <- ifelse(snp.dist$SNP %in% snp.lead$SNP, "Lead eQTL", "Rest")
snp.dist$Type <- factor(snp.dist$Type, levels = c("Rest", "Lead eQTL"))
snp.dist <- snp.dist[snp.dist$Type=="Lead eQTL",]

af <- fread("results/AF/Susie2_AF.txt")

###
g1.eqtl <- ggplot(cs, aes(as.factor(cs), n)) + 
  geom_col(fill=myPal[1], colour="black") + 
  ylab("Number of eGenes") + 
  xlab("Number or credible causal sets") + 
  scale_y_continuous(labels = comma)

g2.eqtl <- ggplot(cs.var, aes(nVar, nCS)) + 
  geom_point(colour=myPal[2], size=2.5) + 
  geom_line() + 
  scale_x_log10(labels=comma) + 
  scale_y_continuous(labels = comma) + 
  annotation_logticks(sides = "b") + 
  xlab("Number of variants") +
  ylab("Number of credible sets")

g3.eqtl <- ggplot(af, aes(x=log2_aFC)) +
  geom_histogram(position = "identity", alpha=.8, fill=myPal[3], colour="black", binwidth = 0.2) + 
  geom_vline(xintercept = c(-2, 2), linetype=2) + 
  ylab("Number of eQTNs") + 
  xlab(expression(italic(log)[2]~aFC))

g1.eqtl + g2.eqtl + g3.eqtl + plot_annotation(tag_levels = "A")

# Save with unified width for consistent perceived text size
ggsave(
  filename = "results/figures/Fig2.png",
  width = 17.5,
  height = 6,
  dpi = 300,
  units = "cm"
)

# Also save a vector PDF copy
ggsave(
  filename = "results/figures/Fig2.pdf",
  width = 17.5,
  height = 6,
  units = "cm"
)
