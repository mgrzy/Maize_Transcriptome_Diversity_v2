#!/usr/bin/env Rscript
# Summarise SuSiE results and join reliable aFCn estimates with PLINK MAF and
# target-gene TSS distances.
#
# The default aFCn standard-error cutoff is 0.5 on the log2(aFC) scale.
# Override it when running, for example:
#   MAX_LOG2_AFC_SE=1 Rscript scripts/local/procesing/8.Process_Susie_aFCn_Results.R

library(data.table)
library(tidyverse)

PIP_PATH <- "results/full/SusieR/PIP_all_batches.tsv.gz"
EGENE_PATH <- "results/summary/osca/eGene.txt"
AFCN_PATH <- "results/full/AF/aFCn.txt"
MAF_PATH <- "data/processed/geno/plink/plink.frq"
GFF_PATH <- "../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
SUSIE_OUTPUT_DIR <- "results/summary/SusieR"
AFCN_OUTPUT_DIR <- SUSIE_OUTPUT_DIR

MAX_LOG2_AFC_SE <- as.numeric(Sys.getenv("MAX_LOG2_AFC_SE", "0.5"))
if (!is.finite(MAX_LOG2_AFC_SE) || MAX_LOG2_AFC_SE < 0) {
  stop("MAX_LOG2_AFC_SE must be a non-negative number.")
}

dir.create(SUSIE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# aFCn variant IDs follow chr_position (for example, chr1_248511).  Keep the
# parser deliberately permissive so it also supports colon-delimited IDs and
# IDs with allele suffixes.
parse_variant_ids <- function(variant_ids) {
  variant_ids <- as.character(variant_ids)
  matches <- regexec("^([^:_]+)[:_]([0-9]+)(?:[:_].*)?$", variant_ids)
  parts <- regmatches(variant_ids, matches)

  data.table(
    variant_chromosome = vapply(parts, function(x) if (length(x)) x[[2L]] else NA_character_, character(1)),
    variant_position = suppressWarnings(as.integer(vapply(
      parts, function(x) if (length(x)) x[[3L]] else NA_character_, character(1)
    )))
  )
}

load_gene_tss <- function(gff_path) {
  if (!file.exists(gff_path)) {
    stop("GFF annotation file not found: ", gff_path)
  }

  annotation <- fread(
    gff_path,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    select = c(1, 3, 4, 5, 7, 9),
    col.names = c("chromosome", "feature", "start", "end", "strand", "attributes")
  )
  annotation <- annotation[feature == "gene"]
  annotation[, gene_id := sub("(^|.*;)ID=([^;]+)(;.*|$)", "\\2", attributes)]
  annotation[!grepl("(^|;)ID=", attributes), gene_id := NA_character_]
  annotation <- annotation[!is.na(gene_id) & strand %in% c("+", "-")]
  annotation[, tss := ifelse(strand == "-", end, start)]

  # A gene should have one annotation row. Retaining the first row makes the
  # join deterministic if a future annotation contains duplicated IDs.
  unique(annotation[, .(gene_id, gene_chromosome = chromosome, gene_tss = tss)], by = "gene_id")
}

pip <- fread(PIP_PATH)
eGene <- fread(EGENE_PATH, header = FALSE, col.names = "gene")
afcn <- fread(AFCN_PATH)

pip <- pip %>%
  filter(gene %in% eGene$gene)

eGene_ids <- unique(eGene$gene)
not_fine_mapped_genes <- setdiff(eGene_ids, unique(pip$gene))
n_not_fine_mapped <- length(not_fine_mapped_genes)

### Number of credible sets per gene
cs <- pip %>%
  group_by(gene) %>%
  summarise(cs = max(cs), .groups = "drop") %>%
  count(cs) %>%
  bind_rows(tibble(cs = 0L, n = n_not_fine_mapped)) %>%
  arrange(cs)

### Number of variants per credible set
cs.var <- pip %>%
  group_by(cs) %>%
  count(gene, name = "nVar") %>%
  ungroup() %>%
  count(nVar, name = "nCS")

cs.gene <- pip %>%
  group_by(gene) %>%
  summarise(cs = max(cs), .groups = "drop")

### aFCn results with MAF
# aFCn can contain identical gene-variant rows more than once. Keep one copy,
# discard non-finite or large standard errors, and then join one MAF per SNP.
afcn <- unique(afcn, by = c("gene_id", "variant_id"))
afcn <- afcn[
  is.finite(log2_aFC_error) & log2_aFC_error <= MAX_LOG2_AFC_SE
]

# plink.frq has ~29 million rows. Stream it once and retain only the variants
# that passed the aFCn standard-error filter instead of loading the full file
# into memory. Column 2 is SNP and column 5 is MAF in PLINK .frq output.
variant_ids_path <- tempfile("afcn_variant_ids_", fileext = ".txt")
fwrite(
  data.table(variant_id = unique(afcn$variant_id)),
  variant_ids_path,
  col.names = FALSE
)
maf_cmd <- paste(
  "awk 'NR == FNR { ids[$1] = 1; next } FNR > 1 && ($2 in ids) { print $2 \"\\t\" $5 }'",
  shQuote(variant_ids_path),
  shQuote(MAF_PATH)
)
maf <- fread(cmd = maf_cmd, col.names = c("variant_id", "MAF"))
unlink(variant_ids_path)

maf <- unique(
  maf,
  by = "variant_id"
)

afcn_maf <- merge(afcn, maf, by = "variant_id", all = FALSE, sort = FALSE)

# Add the signed distance in reference coordinates (variant position - TSS)
# and its absolute value. Both are NA when a variant ID cannot be parsed or
# its target gene is absent from the GFF annotation / lies on another sequence.
variant_coordinates <- parse_variant_ids(afcn_maf$variant_id)
gene_tss <- load_gene_tss(GFF_PATH)
target_gene_index <- match(afcn_maf$gene_id, gene_tss$gene_id)
same_chromosome <- !is.na(target_gene_index) &
  !is.na(variant_coordinates$variant_chromosome) &
  !is.na(variant_coordinates$variant_position) &
  variant_coordinates$variant_chromosome == gene_tss$gene_chromosome[target_gene_index]
distance_to_TSS <- rep(NA_integer_, nrow(afcn_maf))
distance_to_TSS[same_chromosome] <-
  variant_coordinates$variant_position[same_chromosome] - gene_tss$gene_tss[target_gene_index[same_chromosome]]
afcn_maf[, `:=`(
  distance_to_TSS = distance_to_TSS,
  abs_distance_to_TSS = abs(distance_to_TSS)
)]
setcolorder(
  afcn_maf,
  c("gene_id", "variant_id", "log2_aFC", "log2_aFC_error", "log2_aFC_c0", "MAF",
    "distance_to_TSS", "abs_distance_to_TSS")
)
setorder(afcn_maf, gene_id, variant_id)

afcn_maf_path <- file.path(AFCN_OUTPUT_DIR, "aFCn_MAF_se_filtered.tsv.gz")
fwrite(afcn_maf, afcn_maf_path, sep = "\t")
message(
  "Wrote ", nrow(afcn_maf), " gene-variant rows with log2_aFC_error <= ",
  MAX_LOG2_AFC_SE, " to ", afcn_maf_path
)
message(n_not_fine_mapped, " eGenes were not present in the SuSiE PIP results (cs = 0).")

### Save SuSiE summaries
write.csv(pip, file.path(SUSIE_OUTPUT_DIR, "PIP_eGenes.csv"), row.names = FALSE)
write.csv(cs, file.path(SUSIE_OUTPUT_DIR, "CS.csv"), row.names = FALSE)
write.csv(cs.var, file.path(SUSIE_OUTPUT_DIR, "CS_VAR.csv"), row.names = FALSE)
write.csv(cs.gene, file.path(SUSIE_OUTPUT_DIR, "CS_Gene.csv"), row.names = FALSE)
