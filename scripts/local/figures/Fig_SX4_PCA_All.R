### #####################
###  ###
#########################

library(tidyverse)
library(scales)
library(patchwork)
library(ggsci)
library(ggnewscale)

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
dat <- read.csv("data/raw/meta/TableS1.csv")

pca.snp <- read.table("results/full/PCA/PCA_Geno.eigenvec")

dat2 <- dat[dat$ID %in% pca.snp$V1,]
dat <- dat[!dat$ID %in% pca.snp$V1,]

dat$Species <- ifelse(dat$Species=="Z. mays", "Z. mays", 
                      ifelse(dat$Species %in% c("Z. parviglumis", "Teosinte"), "Z. parviglumis", 
                             ifelse(dat$Species=="Z. mexicana","Z. mexicana", "Other maize relatives")))
dat[dat$PC1< -0.05 & dat$Species=="Z. mays",]$Species <- "Other maize relatives"

dat2$Het <- dat2$Group_As_On_Fig2
dat2$Het[dat2$Het=='Other'] <- "SS"
###
g1 <- ggplot() +
  geom_point(data=dat, aes(PC1, PC2, colour=Species), size=1.2,alpha=0.6) + 
  scale_color_d3() + 
  new_scale_color() +
  geom_point(data=dat2, aes(PC1, PC2), size=2, colour=pal_d3(palette = "category10")(6)[5], alpha=0.5)

g2 <- ggplot() +
  geom_point(data=dat, aes(PC2, PC3, colour=Species), size=1.2,alpha=0.6) + 
  scale_color_d3() + 
  new_scale_color() +
  geom_point(data=dat2, aes(PC2, PC3), size=2, colour=pal_d3(palette = "category10")(6)[5], alpha=0.5)

###
combined_plot <- g1 + g2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")

# Export with main figure parameters (PNG + PDF, cm units)
ggsave(
  filename = "results/figures/supp/Fig_PCA_All.png",
  plot = combined_plot,
  width = 18,
  height = 9,
  dpi = 300,
  units = "cm"
)
ggsave(
  filename = "results/figures/supp/Fig_PCA_All.pdf",
  plot = combined_plot,
  width = 18,
  height = 9,
  units = "cm"
)
