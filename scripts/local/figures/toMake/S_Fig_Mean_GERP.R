library(tidyverse)
library(patchwork)

# Match main figure styling
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")
myPal2 <- c("#FC8D62", "#66C2A5")
###
af <- read.csv("results/selection/WGS_AF.csv")

af$highest2 <- af$highest

ann_text <- data.frame(gen.mean.x = c(3.5, 3.5), highest = c(1.5, 1.5), p = c(paste0(expression(italic(P) == 2.07~"x"~10^~-~6)), 
                                                                            paste0(expression(italic(P) == 2.53~"x"~10^~-15))),
                       name=c("M1", "M2"), highest2=c("Recesive", "Dominant"))

ggplot(af, aes(highest, gen.mean.x, fill=highest2)) + 
  geom_boxplot() + 
  geom_text(data = ann_text, aes(y=gen.mean.x, x=highest, label=p), parse = T, size = 2) + 
  scale_fill_manual(values = c("#66C5CC", "#FC8D62")) +
  ylab("Mean GERP score\n per gene") + 
  xlab("Homoeolog") + 
  theme(legend.position = "none") + 
  geom_segment(aes(x=1, xend=2, y=3.3, yend=3.3), linewidth=0.1) + 
  facet_wrap(name~.,)


ggsave(
  filename = "results/figures/Supp/FigS_GERP_mean.png",
  width = 8.75,
  height = 6,
  units = "cm",
  dpi = 300,
  bg = "white"
)
ggsave(
  filename = "results/figures/Supp/FigS_GERP_mean.pdf",
  width = 8.75,
  height = 6,
  units = "cm",
  bg = "white"
)
### Script: figures/S_Fig_Mean_GERP.R
### Purpose: Supplemental figure showing distribution/summary of mean GERP scores across sets.
### Inputs: GERP annotations merged with targets
### Outputs: results/figures/S_Fig_Mean_GERP.png
### Usage: Rscript figures/S_Fig_Mean_GERP.R
### Dependencies: tidyverse and plotting libs
