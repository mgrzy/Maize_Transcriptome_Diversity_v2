library(tidyverse)
library(scales)
library(patchwork)
library(data.table)
library(dtplyr)

# Match main figure styling
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")

eqtl <- fread("results/eQTL/Final_eQTL_IHW.csv.gz")

rl <- eqtl %>% 
  group_by(Probe) %>% summarise(L=abs(min(BP)-max(BP)))

eqtl$dist <- eqtl$BP - eqtl$Probe_bp

g1 <- ggplot(eqtl, aes(x = dist)) + 
  geom_histogram(binwidth = 50000, alpha=0.5, colour=myPal[1]) + 
  scale_x_continuous(labels = paste0(c(-2,-1,0,1,2), "Mb"),
                     breaks = 10^6 * c(-2,-1,0,1,2)) + 
  xlab(expression(italic(cis)-eQTN~distance~to~TSS))

g2 <- ggplot(eqtl, aes(x = Freq)) + 
  geom_histogram(binwidth = 0.025, alpha=0.5, colour=myPal[2]) + 
  xlab("Minor Allele Freqency")

g3 <- ggplot(rl, aes(x = L)) + 
  geom_histogram(binwidth = 25000, alpha=0.5, colour=myPal[3]) + 
  scale_x_continuous(labels = paste0(c(1,2,3,4), "Mb"),
                     breaks = 10^6 * c(1,2,3,4)) + 
  xlab(expression(italic(cis)-eQTL~lenght))

combined_plot <- g1 + g2 + g3 + plot_annotation(tag_levels = "A")

# Export with main figure parameters (PNG + PDF, cm units)
ggsave(
  filename = "results/figures/Supp/FigS_eQTL_length.png",
  plot = combined_plot,
  width = 17.5,
  height = 6,
  dpi = 300,
  units = "cm"
)
ggsave(
  filename = "results/figures/Supp/FigS_eQTL_length.pdf",
  plot = combined_plot,
  width = 17.5,
  height = 6,
  units = "cm"
)
### Script: figures/S_Fig_eQTL_Length.R
### Purpose: Supplemental figure showing distribution of distances/lengths for lead eQTL signals relative to genes.
### Inputs: results/SusieR2/* or eQTL summary with distances
### Outputs: results/figures/S_Fig_eQTL_Length.png
### Usage: Rscript figures/S_Fig_eQTL_Length.R
### Dependencies: tidyverse, ggplot2
