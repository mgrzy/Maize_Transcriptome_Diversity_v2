#!/usr/bin/env Rscript

# Combine Fig1 and Fig2 panels into one 5-panel layout (A-E)

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
})

# Consistent styling for the merged figure
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
# Harmonize global font sizing with Fig5.R: only set axis text colour
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#66C5CC", "#F6CF71", "#DCB0F2", "#87C55F", "#9EB9F3", "#C3B091")

# --- Fig1 panels ----------------------------------------------------------

# PCA
pc <- read.table("results/full/PCA/PCA_Geno.eigenvec") %>%
  select(-V1)
colnames(pc) <- c("Taxa", "PC1", "PC2", "PC3", "PC4", "PC5")

# Meta data
meta <- read.csv("data/raw/meta/WiDiv_meta.csv")
pc <- merge(meta, pc, by = "Taxa") %>%
  mutate(
    Group = str_replace_all(Group, "Broad[_ ]origin-public", "Broad origin"),
    Group = factor(Group, levels = c("Tropical", "Broad origin", "Mixed", "SS", "NSS", "IDT"))
  )

# CCA
cca <- read.csv("results/full/CCA/CCA.csv")

# R2
r2 <- read.csv("results/full/GeneExpressionVariance/All_R2.csv")

# PCA plot: PC1 vs PC2 with legend inside to save space
g1 <- ggplot(pc, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 2) +
  scale_color_manual(values = myPal, name = "Subpopulation") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(
    legend.position = c(0.26, 0.78),
    legend.background = element_rect(fill = alpha("white", 0.1), colour = NA),
    legend.text = element_text(size = base_font_size - 2),
    legend.title = element_text(size = base_font_size - 2),
    legend.key.size = grid::unit(0.35, "cm")
  )

# CCA plot
g2 <- ggplot(cca, aes(CCA1, CCA2, colour = Group)) +
  geom_point(size = 2) +
  scale_color_manual(values = myPal) +
  labs(x = "Coordinate 1", y = "Coordinate 2") +
  theme(legend.position = "none")

# Histogram of R^2 values
g3 <- ggplot(r2, aes(x = r2)) +
  geom_histogram(binwidth = 0.0025, alpha = 0.7, colour = myPal[2], fill = myPal[2]) +
  labs(x = "Proportion variance explained", y = "Count")

# --- Fig2 panels ----------------------------------------------------------

cv <- read_csv("results/full/GeneExpressionVariance/All_CV.csv") %>%
  mutate(
    Group = str_replace_all(Group, "Broad[_ ]origin-public", "Broad origin"),
    Group = factor(Group, levels = c("Tropical", "Broad origin", "Mixed", "SS", "NSS", "IDT"))
  )

n.div <- read_csv("results/full/Pi//Pi_all.csv")

n.div.pop <- n.div %>%
  filter(pop %in% c("Tropical", "Broad_origin-public", "Broad origin-public", "Mixed", "SS", "NSS", "IDT")) %>%
  mutate(
    pop = str_replace_all(pop, "Broad[_ ]origin-public", "Broad origin"),
    pop = factor(pop, levels = levels(cv$Group))
  )

# Nucleotide diversity by subpopulation
g.pi <- ggplot(n.div.pop, aes(pop, avg_pi, fill = pop)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 0.6) +
  scale_fill_manual(values = myPal) +
  scale_x_discrete(labels = c("Tropical", "Broad\norigin", "Mixed", "SS", "NSS", "Iodent")) +
  labs(x = "Subpopulation", y = expression(Nucleotide~diversity~(pi))) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.1))

# Transcriptomic CV by subpopulation
g.cv <- ggplot(cv, aes(Group, cv, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 0.6) +
  scale_fill_manual(values = myPal) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  scale_x_discrete(labels = c("Tropical", "Broad\norigin", "Mixed", "SS", "NSS", "Iodent")) +
  labs(x = "Subpopulation", y = "Coefficient of variation") +
  theme(legend.position = "none")

# --- Combined figure layout ----------------------------------------------

top_row <- wrap_plots(g1, g2, g3, nrow = 1)
bottom_row <- wrap_plots(g.pi, g.cv, nrow = 1)

combined_plot <- top_row / bottom_row +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = base_font_size, face = "plain"))

# Export figure
output_path <- "results/figures/Fig1.png"

ggsave(
  filename = output_path,
  plot = combined_plot,
  width = 17.5,
  height = 12,
  dpi = 300,
  units = "cm"
)

# Also save a vector PDF copy
ggsave(
  filename = sub("\\.png$", ".pdf", output_path),
  plot = combined_plot,
  width = 17.5,
  height = 12,
  units = "cm"
)

message("Combined figure saved to ", output_path)
