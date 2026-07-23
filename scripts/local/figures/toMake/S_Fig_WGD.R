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
pip <- read.csv("results/selection/WGS_PIP.csv")

af2 <- read.csv("results/selection/AF_Selection.csv")
pip2 <- read.csv("results/selection/PIP_Selection.csv")

pip.frac2 <- pip2 %>% 
  group_by(Sorghum_Ortholog, cs) %>%
  summarise(n=n()) %>% 
  mutate(freq=n/sum(n))

pip.frac <- pip %>% 
  group_by(name, highest, cs) %>%
  summarise(n=n()) %>% 
  mutate(freq=n/sum(n))

###
g1 <- ggplot(af2, aes(x=Sorghum_Ortholog, y=abs(log2_aFC), fill = Sorghum_Ortholog)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("#66C5CC", "#FC8D62"), name="Sorghum Ortholog") + 
  ylab(expression("|" ~ log[2](aFC) ~ "|")) + 
  annotate("text",x = 1.5, y=7.5, label = paste0(expression(italic(P) ==1.6~"x"~10^~-219)), parse=T) + 
  xlab("Sorgum ortholog") + 
  theme(legend.position = "none") + 
  geom_segment(aes(x=1, xend=2, y=7, yend=7), linewidth=0.1)
  
g2 <- ggplot(pip.frac2, aes(as.factor(cs), freq, fill=Sorghum_Ortholog)) + 
  geom_col(position = position_dodge()) + 
  scale_fill_manual(values = c("#66C5CC", "#FC8D62"), name="Sorghum Ortholog") + 
  theme(legend.position.inside = c(0.8, 0.8), legend.position = "inside") + 
  xlab("Number of credible sets") + 
  annotate("text",x = 8, y=0.125, label = paste0(expression(italic(P)==0.00695)), parse=T) + 
  ylab("Proportion of eGenes") 


###
af$highest2 <- af$highest
ann_text <- data.frame(log2_aFC = c(7.5, 7.5), highest = c(1.5, 1.5), p = c(paste0(expression(italic(P) > 0.05)), paste0(expression(italic(P) > 0.05))),
                                                                            name=c("M1", "M2"), highest2=c("Recesive", "Dominant"))

###
g3 <- ggplot(af, aes(x=highest, y=abs(log2_aFC), fill = highest2)) +
  geom_boxplot() + 
  geom_text(data = ann_text, aes(y=log2_aFC, x=highest, label=p), parse = T) + 
  scale_fill_manual(values = myPal[3:4]) + 
  ylab(expression("|" ~ log[2](aFC) ~ "|")) + 
  #annotate("text",x = 1.5, y=7.5, label = paste0(expression(italic(P) ==1.6~"x"~10^~-219)), parse=T) + 
  xlab("Homoeolog") + 
  theme(legend.position = "none") + 
  geom_segment(aes(x=1, xend=2, y=7, yend=7), linewidth=0.1) + 
  facet_wrap(name~.,)

g4 <- ggplot(pip.frac, aes(as.factor(cs), freq, fill=highest)) + 
  geom_col(position = position_dodge()) + 
  scale_fill_manual(values = myPal[3:4], name="Homoeolog") + 
  theme(legend.position.inside = c(0.8, 0.8), legend.position = "inside") + 
  xlab("Number of credible sets") + 
  annotate("text",x = 8, y=0.125, label = paste0(expression(italic(P) > 0.05)), parse=T) + 
  ylab("Proportion of eGenes") + facet_wrap(name~.,)

combined_plot <- g1 + g2 + g3 + g4 + plot_annotation(tag_levels = "A") + plot_layout(widths = c(0.5, 1))

# Export with main figure parameters (PNG + PDF, cm units)
ggsave(
  filename = "results/figures/Supp/FigS_WGD_effect.png",
  plot = combined_plot,
  width = 17.5,
  height = 12,
  dpi = 300,
  units = "cm"
)
ggsave(
  filename = "results/figures/Supp/FigS_WGD_effect.pdf",
  plot = combined_plot,
  width = 17.5,
  height = 12,
  units = "cm"
)
### Script: figures/S_Fig_WGD.R
### Purpose: Supplemental figure related to whole-genome duplication (WGD) analyses.
### Inputs: relevant WGD annotation and analysis outputs
### Outputs: results/figures/S_Fig_WGD.png
### Usage: Rscript figures/S_Fig_WGD.R
### Dependencies: tidyverse and plotting libs used
