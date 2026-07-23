library(tidyverse)

# Match main figure styling
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")

###
pcaAll <- read.csv("results/RNA_PCA/PCAAll.csv")
pcaAll$Group <- factor(pcaAll$Group, levels = c("Tropical", "Broad_origin-public", "Mixed", "SS", "NSS", "IDT"))

###
p <- ggplot(pcaAll, aes(Dim.1, Dim.2, colour=Group)) + 
  geom_point() + 
  scale_color_manual(values = myPal, name="Subpopulation") +
  xlab("PC1") + ylab("PC2")

# Export with main figure parameters (PNG + PDF, cm units)
ggsave(
  filename = "results/figures/Supp/FigS_PCA_RNA_All.png",
  plot = p,
  width = 8.75,
  height = 6,
  dpi = 300,
  units = "cm"
)
ggsave(
  filename = "results/figures/Supp/FigS_PCA_RNA_All.pdf",
  plot = p,
  width = 8.75,
  height = 6,
  units = "cm"
)
### Script: figures/S_Fig_PCA_RNA_All.R
### Purpose: Supplemental figure for PCA on RNA data across all genes/samples.
### Inputs: PCA results from expression (e.g., results/PCA/PCAAll.csv)
### Outputs: results/figures/S_Fig_PCA_RNA_All.png
### Usage: Rscript figures/S_Fig_PCA_RNA_All.R
### Dependencies: tidyverse, ggplot2
