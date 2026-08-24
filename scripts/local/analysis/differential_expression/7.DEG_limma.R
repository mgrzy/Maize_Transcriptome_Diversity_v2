### Script: analysis/7.DEG_limma.R
### Purpose: Run limma-voom differential-expression tests for heterotic-group
###          and breeding-era contrasts defined in WiDiv_meta.csv.
### Inputs:
###   - data/raw/rna/counts_NE2020_corrected_names.csv.gz
###   - data/raw/meta/WiDiv_meta.csv
### Outputs:
###   - results/full/DEG/*.csv
### Usage: Rscript scripts/local/analysis/7.DEG_limma.R
### Dependencies: edgeR, limma

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

counts_path <- "data/raw/rna/counts_NE2020_corrected_names.csv.gz"
metadata_path <- "data/raw/meta/WiDiv_meta.csv"
output_dir <- "results/full/DEG"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_counts <- function(path) {
  counts_df <- read.csv(path, check.names = FALSE)
  gene_ids <- counts_df[[1L]]
  counts <- as.matrix(counts_df[-1L])
  rownames(counts) <- gene_ids
  storage.mode(counts) <- "numeric"
  counts
}

run_limma_voom <- function(counts, sample_metadata, comparison_name,
                           numerator, denominator) {
  group <- factor(
    sample_metadata$comparison_group,
    levels = c(denominator, numerator)
  )
  group_sizes <- table(group)

  if (any(group_sizes < 2L)) {
    warning(
      comparison_name, ": at least one group has fewer than two lines (",
      paste(names(group_sizes), group_sizes, sep = "=", collapse = ", "),
      "). Interpret this contrast cautiously."
    )
  }

  y <- DGEList(
    counts = counts[, sample_metadata$Taxa, drop = FALSE],
    samples = sample_metadata
  )
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  keep_genes <- filterByExpr(y, design = design)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  voom_fit <- voom(y, design = design, plot = FALSE)
  contrast <- makeContrasts(
    contrasts = paste(numerator, denominator, sep = " - "),
    levels = design
  )
  fit <- contrasts.fit(lmFit(voom_fit, design), contrast)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  result <- topTable(
    fit,
    coef = 1L,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  result <- data.frame(gene = rownames(result), result, check.names = FALSE)

  write.csv(
    result,
    file.path(output_dir, paste0(comparison_name, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
  write.csv(
    data.frame(
      Taxa = sample_metadata$Taxa,
      Group = sample_metadata$Group,
      Era = sample_metadata$Era,
      comparison_group = group,
      check.names = FALSE
    ),
    file.path(output_dir, paste0(comparison_name, "_samples.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  message(
    "Wrote ", comparison_name, " (",
    paste(names(group_sizes), group_sizes, sep = "=", collapse = ", "), ")"
  )
}

# Read raw counts and map technical replicate names (e.g. *_R2) to WiDiv taxa.
counts <- read_counts(counts_path)
sample_taxa <- sub("_R2$", "", colnames(counts))
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata$Era[metadata$Era == "Public Sector"] <- "Public"
metadata$Era[metadata$Era == ""] <- NA_character_

metadata_index <- match(sample_taxa, metadata$Taxa)
matched <- !is.na(metadata_index)
if (!any(matched)) {
  stop("No count-matrix sample names matched WiDiv_meta.csv Taxa values.")
}

counts <- counts[, matched, drop = FALSE]
sample_metadata <- metadata[metadata_index[matched], c("Taxa", "Group", "Era")]

# Sum technical replicates to one count profile per genotype line.
line_by_gene <- rowsum(t(counts), group = sample_metadata$Taxa, reorder = FALSE)
counts <- t(line_by_gene)
line_metadata <- metadata[match(colnames(counts), metadata$Taxa), c("Taxa", "Group", "Era")]

if (anyNA(line_metadata$Taxa)) {
  stop("Some aggregated count profiles could not be matched back to metadata.")
}

comparison_specs <- list(
  SS_vs_Tropical = list(
    keep = line_metadata$Group %in% c("SS", "Tropical"),
    group = line_metadata$Group,
    numerator = "SS", denominator = "Tropical"
  ),
  NSS_vs_Tropical = list(
    keep = line_metadata$Group %in% c("NSS", "Tropical"),
    group = line_metadata$Group,
    numerator = "NSS", denominator = "Tropical"
  ),
  Iodent_vs_Tropical = list(
    keep = line_metadata$Group %in% c("IDT", "Tropical"),
    group = ifelse(line_metadata$Group == "IDT", "IDT", "Tropical"),
    numerator = "IDT", denominator = "Tropical"
  ),
  SS_exPVP_vs_Public = list(
    keep = line_metadata$Group == "SS" &
      line_metadata$Era %in% c("exPVP", "Public"),
    group = line_metadata$Era,
    numerator = "exPVP", denominator = "Public"
  ),
  NSS_exPVP_vs_Public = list(
    keep = line_metadata$Group == "NSS" &
      line_metadata$Era %in% c("exPVP", "Public"),
    group = line_metadata$Era,
    numerator = "exPVP", denominator = "Public"
  ),
  Public_SS_vs_NSS = list(
    keep = line_metadata$Group %in% c("SS", "NSS") &
      line_metadata$Era == "Public",
    group = line_metadata$Group,
    numerator = "SS", denominator = "NSS"
  ),
  exPVP_SS_vs_NSS = list(
    keep = line_metadata$Group %in% c("SS", "NSS") &
      line_metadata$Era == "exPVP",
    group = line_metadata$Group,
    numerator = "SS", denominator = "NSS"
  )
)

for (comparison_name in names(comparison_specs)) {
  spec <- comparison_specs[[comparison_name]]
  selected <- which(spec$keep %in% TRUE)
  comparison_metadata <- line_metadata[selected, , drop = FALSE]
  comparison_metadata$comparison_group <- spec$group[selected]

  if (!all(c(spec$numerator, spec$denominator) %in%
    unique(comparison_metadata$comparison_group))) {
    stop(comparison_name, ": one or both comparison groups are absent.")
  }

  run_limma_voom(
    counts = counts,
    sample_metadata = comparison_metadata,
    comparison_name = comparison_name,
    numerator = spec$numerator,
    denominator = spec$denominator
  )
}
