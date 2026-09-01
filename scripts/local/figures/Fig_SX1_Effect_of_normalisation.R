library(tidyverse)
library(data.table)
library(scales)
library(patchwork)

base_font_size <- 9
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.text = element_text(colour = "black", size = base_font_size, face = "plain"),
  axis.title = element_text(size = base_font_size, face = "plain"),
  axis.ticks = element_line(colour = "black"),
  legend.title = element_text(size = base_font_size, face = "plain"),
  legend.text = element_text(size = base_font_size, face = "plain"),
  plot.tag = element_text(size = base_font_size, face = "plain"),
  plot.title.position = "plot"
)


data.log2 <- fread("results/full/GeneExpressionVariance/edgeR_log2_cpm_filtered_mean_sd.csv.gz")
data.vst <- fread("results/full/GeneExpressionVariance/DESeq2_vst_filtered_mean_sd.csv.gz")
data.raw <- fread("results/full/GeneExpressionVariance/raw_expression_mean_sd.csv.gz")


###
data.raw <- data.raw[data.raw$Gene %in% data.log2$Gene,]

###
g1 <- ggplot(data.raw, aes(Mean, SD)) + 
  geom_point(size = 0.05) + 
  geom_smooth() + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides = "b") + 
  labs(x = expression(mu), y=expression(sigma))   

g2 <- ggplot(data.log2, aes(Mean, SD)) + 
  geom_point(size = 0.05) + 
  geom_smooth() + 
  labs(x = expression(mu), y=expression(sigma))   


g3 <- ggplot(data.vst, aes(Mean, SD)) + 
  geom_point(size = 0.05) + 
  geom_smooth() + 
  labs(x = expression(mu), y=expression(sigma))   

combined_plot <- g1 + g2 + g3 + plot_annotation(tag_levels = "A")

ggsave(
  filename = "results/figures/supp/Fig_GeneExpressionVaraianceNorm.png",
  plot = combined_plot,
  width = 21,
  height = 6,
  units = "cm",
  dpi = 300
)
