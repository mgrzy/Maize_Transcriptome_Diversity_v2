library(data.table)
library(ggplot2)

fm <- fread("results/summary/SusieR/aFCn_MAF_se_filtered.tsv.gz")

eqtl <- fread("results/summary/osca/eqtl_padjusted_full_results.tsv.gz")


colnames(fm)[7] <- "distance_TSS"
fm$distance_TSS <- abs(fm$distance_TSS)
colnames(eqtl)[18] <- "distance_TSS"


breaks <- c(
  0,
  1e3,
  2e3,
  5e3,
  10e3,
  50e3,
  100e3,
  250e3,
  500e3,
  1e6,
  2e6,
  Inf
)

labels <- c(
  "0–1 kb",
  "1–2 kb",
  "2–5 kb",
  "5–10 kb",
  "10–50 kb",
  "50–100 kb",
  "100–250 kb",
  "250–500 kb",
  "0.5–1 Mb",
  "1–2 Mb",
  ">2 Mb"
)

eqtl[, distance_bin := cut(
  distance_TSS,
  breaks = breaks,
  labels = labels,
  right = FALSE
)]

fm[, distance_bin := cut(
  distance_TSS,
  breaks = breaks,
  labels = labels,
  right = FALSE
)]

background_summary <- eqtl[
  ,
  .(n_background = .N),
  by = distance_bin
]

fine_summary <- fm[
  ,
  .(n_fm = .N),
  by = distance_bin
]

enrichment <- merge(
  background_summary,
  fine_summary,
  by = "distance_bin",
  all = TRUE
)

enrichment[is.na(n_background), n_background := 0L]
enrichment[is.na(n_fm), n_fm := 0L]

enrichment[
  ,
  background_fraction := n_background / sum(n_background)
]

enrichment[
  ,
  fm_fraction := n_fm / sum(n_fm)
]

enrichment[
  ,
  fold_enrichment := fm_fraction / background_fraction
]

enrichment[
  ,
  log2_fold_enrichment := log2(fold_enrichment)
]

enrichment

write.csv(enrichment, "results/summary/SusieR/enrichment_distance.csv", row.names = F)

all_bins <- levels(eqtl$distance_bin)

or_results <- rbindlist(lapply(all_bins, function(current_bin) {
  
  a <- fm[distance_bin == current_bin, .N]
  b <- fm[distance_bin != current_bin, .N]
  
  c <- eqtl[distance_bin == current_bin, .N]
  d <- eqtl[distance_bin != current_bin, .N]
  
  test <- fisher.test(
    matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  )
  
  data.table(
    distance_bin = current_bin,
    odds_ratio = unname(test$estimate),
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2],
    p_value = test$p.value
  )
}))

or_results[
  ,
  p_fdr := p.adjust(p_value, method = "BH")
]
