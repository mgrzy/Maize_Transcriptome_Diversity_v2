library(tidyverse)
library(patchwork)

### (match Fig1 font sizing)
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
af <- read.csv("results/selection/AF_Selection.csv")
pip <- read.csv("results/selection/PIP_Selection.csv")

###

af$constrained <- factor(af$constrained, levels = c("Top10", "Bottom90"), labels = c("Top 10%", "Lower 90%"))

g1.effect <- ggplot(af, aes(constrained, abs(log2_aFC), fill=constrained)) +
  geom_boxplot(width=0.5) + 
  scale_fill_manual(values = c("#66C5CC", "#F6CF71")) + 
  ylab(expression("|" ~ log[2](aFC) ~ "|")) + 
  annotate("text",x = 1.5, y=7.5, label = paste0(expression(italic(P) ==8.8~"x"~10^~-13)), parse=T) + 
  xlab("GERP percentile") + 
  theme(legend.position = "none") + 
  geom_segment(aes(x=1, xend=2, y=7, yend=7), linewidth=0.1)

pip.frac <- pip %>% 
  group_by(constrained, cs) %>%
  summarise(n=n()) %>% 
  mutate(freq=n/sum(n))

pip.frac$constrained <- factor(pip.frac$constrained, levels = c("Top10", "Bottom90"), labels = c("Top 10%", "Lower 90%"))

g2.cs <- ggplot(pip.frac, aes(as.factor(cs), freq, fill=constrained)) + 
  geom_col(position = position_dodge(), colour="black") + 
  scale_fill_manual(values = c("#66C5CC", "#F6CF71"), name="GERP percentile") + 
  theme(legend.position.inside = c(0.8, 0.8), legend.position = "inside") + 
  xlab("Number of credible sets") + 
  annotate("text",x = 8, y=0.125, label = paste0(expression(italic(P) == 9.7~"x"~10^~-5)), parse=T) + 
  ylab("Proportion of eGenes")

g1.effect + g2.cs + plot_annotation(tag_levels = "A") + plot_layout(widths = c(0.25, 1))

# Save with unified export parameters for consistent perceived text size
ggsave(
  filename = "results/figures/Fig5.png",
  width = 17.5,
  height = 6,
  dpi = 300,
  units = "cm"
)

# Also save a vector PDF copy
ggsave(
  filename = "results/figures/Fig5.pdf",
  width = 17.5,
  height = 6,
  units = "cm"
)
### Script: figures/FigX.R
### Purpose: Prototype or auxiliary figure script; adapt to generate specific figure from processed results.
### Inputs: depends on the figure being prototyped
### Outputs: results/figures/FigX.png (or targeted filename)
### Usage: Rscript figures/FigX.R
### Dependencies: tidyverse and plotting libraries
