library(tidyverse)
library(data.table)
library(car)
library(patchwork)

# Global plotting defaults (match Fig1) ------------------------------------
base_font_size <- 8
theme_set(theme_classic(base_size = base_font_size))
theme_update(
  axis.title = element_text(size = rel(1.15)),
  legend.title = element_text(size = rel(1.1)),
  strip.text = element_text(size = rel(1.05))
)

myPal <- c("#A8A8D1", "#E78AC3", "#80B1D3", "#F4A261", "#A6D854", "#C3B091")

# Data import --------------------------------------------------------------
deg <- read_csv("results/DEG/DEG_SS_vs_NSS_by_Era.csv", show_col_types = FALSE)
eqtl <- read_csv("results/SusieR2/PIP_Max.csv", show_col_types = FALSE)
af <- fread("results/AF/Susie2_AF.txt") %>%
  as_tibble()

pub_fst <- fread("results/Fst/SS_vs_NSS_Public_eqtl.fst") %>%
  as_tibble() %>%
  mutate(era = "Era I")
expvp_fst <- fread("results/Fst/SS_vs_NSS_exPVP_eqtl.fst") %>%
  as_tibble() %>%
  mutate(era = "Era II")

fst <- bind_rows(pub_fst, expvp_fst) %>%
  select(1:4, 7, 8)

# Harmonise era labels and merge datasets ----------------------------------
deg <- deg %>%
  mutate(era = dplyr::recode(era, "Public" = "Era I", "exPVP" = "Era II"))

###
res <- merge(deg, eqtl, by="gene")
res <- merge(res, fst, by=c("SNP", "era"))

### 
colnames(af)[1:2] <- c("gene", "SNP")

res <- deg %>%
  inner_join(eqtl, by = "gene") %>%
  inner_join(fst, by = c("SNP", "era"))

res2 <- res %>%
  inner_join(af, by = c("gene", "SNP")) %>%
  mutate(
    DEG = if_else(padj < 0.05, "DEG", "non-DEG"),
    DEG = factor(DEG, levels = c("non-DEG", "DEG"))
  )

# Fst enrichment across DEG expression bins --------------------------------
res_sig <- res %>%
  filter(padj < 0.05) %>%
  group_by(era) %>%
  mutate(q = ntile(dplyr::desc(padj), 10)) %>%
  ungroup()

g1 <- ggplot(res_sig, aes(as.factor(q), Fst, fill = era)) +
  geom_boxplot() +
  facet_wrap(~era, ncol = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values = myPal) +
  xlab("Differential expression -log10(p) decline") +
  ylab(expression(italic(F)[st] ~ "of lead eQTL"))

# Save plot 
ggsave(
  filename = "results/figures/Fig4.png",
  plot = g1,
  width = 8.75,
  height = 6,
  dpi = 300,
  units = "cm",
  bg = "white"
)

# Also save a vector PDF copy
ggsave(
  filename = "results/figures/Fig4.pdf",
  plot = g1,
  width = 8.75,
  height = 6,
  units = "cm",
  bg = "white"
)
