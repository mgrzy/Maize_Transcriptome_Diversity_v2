suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(data.table)
})

# Consistent styling for the merged figure
base_font_size <- 9
theme_set(theme_classic(base_size = base_font_size))
# Harmonize global font sizing with Fig5.R: only set axis text colour
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

res <- fread("results/full/osca/significant_summary/significant_gene_counts.tsv")

res$SignificantGenesProp <- res$SignificantGenes / 15587

g1 <- ggplot(res, aes(SampleSize, SignificantGenesProp)) + 
  geom_point() + 
  ylab("eQTL discoverability") + 
  xlab("Sample size") + 
  scale_x_continuous(labels = res$SampleSize, breaks = res$SampleSize)

g2 <- ggplot(res, aes(SampleSize, MeanSignificantSNPsPerGene)) + 
  geom_point() + 
  ylab("Mean Significant\n SNPs Per Gene") + 
  xlab("Sample size") + 
  scale_x_continuous(labels = res$SampleSize, breaks = res$SampleSize)

g1 + g2 + plot_annotation(tag_levels = "A")

ggsave(
  filename = "results/figures/supp/Fig_eQTL_discoverability.png",
  width = 18,
  height = 9,
  dpi = 300,
  units = "cm"
)
