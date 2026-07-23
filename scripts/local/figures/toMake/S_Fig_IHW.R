library(tidyverse)

dat <- data.table::fread("results/eQTL/IHWcov.csv.gz")

dat$dist <- factor(dat$group, levels = c(dat$group)[1:20])

# Match main figure styling
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")

p <- ggplot(dat, aes(dist, weight, colour=as.character(fold), group=fold)) + 
  geom_point() + geom_line() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) + 
  xlab("Distance to TSS (bp)") + 
  scale_color_manual(values = myPal, name="Fold")

# Export with main figure parameters (PNG + PDF, cm units)
ggsave(
  filename = "results/figures/Supp/FigS_IHW.png",
  plot = p,
  width = 16,
  height = 8,
  dpi = 300,
  units = "cm"
)
ggsave(
  filename = "results/figures/Supp/FigS_IHW.pdf",
  plot = p,
  width = 17.5,
  height = 6,
  units = "cm"
)
### Script: figures/S_Fig_IHW.R
### Purpose: Supplemental figure summarizing IHW weights/adjustments and their impact on eQTL discovery.
### Inputs: IHW results
### Outputs: results/figures/S_Fig_IHW.png
### Usage: Rscript figures/S_Fig_IHW.R
### Dependencies: tidyverse, ggplot2
