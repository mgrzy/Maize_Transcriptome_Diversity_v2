suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(data.table)
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


res <- fread("results/summary/osca/eqtl_padjusted_full_results.tsv.gz")

###
snps_per_gene <- res[
  ,
  .(n_snps = .N),
  by = Gene
]

g1 <- ggplot(snps_per_gene, aes(x = n_snps)) +
  geom_histogram(
    binwidth = 1000,
    boundary = 0,
    fill = myPal[1], 
    colour="black"
  ) +
  labs(
    x = "Number of associated SNPs per eGene",
    y = "Number of eGenes"
  ) 

g2 <- ggplot(res, aes(x = Probe_bp - BP)) +
  geom_histogram(
    binwidth = 50000,
    boundary = 0,
    fill = myPal[2], 
    colour="black"
  ) +
  labs(
    x = "Distance to transcription start site",
  ) 

g3 <- ggplot(res, aes(x = Freq)) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0,
    fill = myPal[3], 
    colour="black"
  ) +
  labs(
    x = "Minor allele freqency",
  ) 

gall <- g1 + g2 + g3 + plot_annotation(tag_levels = "A")

ggsave(plot = gall,
  filename = "results/figures/supp/Fig_eGene_Prop.png", 
  width = 21,
  height = 6,
)

###
summary(as.numeric(table(res$Gene)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2    1362    3156    3967    5651   29932
### 


#Bonfferoni correction over total valid p-value tests: 917540789
bon <- 0.05 / 917540789

resBon <- res[res$p< bon,]

length(unique(resBon$Gene))
# 14690
###  save eGene list

eGenes <- unique(res$Gene)

write.table(eGenes, "results/summary/osca/eGene.txt", row.names = F, col.names = F, quote = F)
