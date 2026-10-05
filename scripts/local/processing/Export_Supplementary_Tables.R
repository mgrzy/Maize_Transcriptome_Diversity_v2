#!/usr/bin/env Rscript
# Export supplementary results to one workbook and a full variant CSV.
# Run from the repository root:
# Rscript scripts/local/processing/Export_Supplementary_Tables.R [output.xlsx]
# Dependencies: install.packages(c("data.table", "openxlsx"))

suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop('Install the Excel writer with install.packages("openxlsx")')
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply at most one output path")
output <- if (length(args)) args[[1L]] else "results/tables/Supplementary_Tables_S1_S3.xlsx"
if (!grepl("\\.xlsx$", output, ignore.case = TRUE)) stop("Output must end in .xlsx")
csv_output <- sub("\\.xlsx$", "_Full_fine_mapped_variants.csv", output, ignore.case = TRUE)

# S1: saved multivariate scores with group, era and flowering-time metadata.
cca <- fread("results/full/CCA/CCA.csv")
expression <- fread("results/full/GeneExpressionVariance/PCAAll.csv")
genotype <- fread("results/full/PCA/PCA_Geno.eigenvec", header = FALSE,
                  colClasses = c(V1 = "character", V2 = "character"))
genotype[, V1 := NULL]
setnames(genotype, c("Taxa", paste0("Genotype_PC", seq_len(ncol(genotype) - 1L))))
for (x in list(cca, expression, genotype)) {
  if (anyNA(x$Taxa) || anyDuplicated(x$Taxa)) stop("Missing or duplicated sample IDs")
  if (!setequal(cca$Taxa, x$Taxa)) stop("CCA and PCA sample sets differ")
}
if (!identical(cca$Group, expression$Group[match(cca$Taxa, expression$Taxa)])) {
  stop("CCA and expression PCA group labels disagree")
}
expression[, Group := NULL]
pc_columns <- grep("^Dim[.]", names(expression), value = TRUE)
setnames(expression, pc_columns, sub("^Dim[.]", "Expression_PC", pc_columns))
s1 <- merge(cca, expression, by = "Taxa")
s1 <- merge(s1, genotype, by = "Taxa")

# Full CSV: every saved fine-mapped variant/credible-set membership.
# Keep gene and credible-set identifiers because a SNP can occur in several sets.
full_variants <- fread("results/full/SusieR/PIP_all_batches.tsv.gz",
            select = c("gene", "SNP", "cs", "variable_prob"))
setnames(full_variants, c("gene", "SNP", "cs", "variable_prob"),
         c("Gene", "Variant", "Credible_set", "PIP"))
setorder(full_variants, Gene, Credible_set, -PIP, Variant)

# S2: use the saved lead per credible set, preserving upstream PIP tie choices.
s2 <- fread("results/summary/SusieR/PIP_lead_variant_per_credible_set.tsv.gz",
            select = c("gene", "SNP", "cs", "variable_prob"))
setnames(s2, c("gene", "SNP", "cs", "variable_prob"),
         c("Gene", "Variant", "Credible_set", "PIP"))
if (anyDuplicated(s2, by = c("Gene", "Credible_set"))) {
  stop("More than one saved lead variant per credible set")
}
setorder(s2, Gene, Credible_set, Variant)

# S3: the aFCn results, with no extra flags or joins.
s3 <- fread("results/full/AF/aFCn.txt")
setorder(s3, gene_id, variant_id)

tables <- list(S1_Multivariate = s1, S2_Lead_variants = s2, S3_aFCn = s3)
# Round only continuous numeric results; preserve integer identifiers and counts.
# Source analysis files retain their original precision.
round_results <- function(x) {
  for (column in names(x)[vapply(x, is.double, logical(1))]) {
    set(x, j = column, value = round(x[[column]], digits = 5L))
  }
  x
}
tables <- lapply(tables, round_results)
full_variants <- round_results(full_variants)
# Split long tables only because Excel allows 1,048,576 rows including headers.
max_rows <- 1048575L
sheets <- list()
for (name in names(tables)) {
  x <- tables[[name]]
  if (!nrow(x)) stop("Empty result table: ", name)
  parts <- ceiling(nrow(x) / max_rows)
  for (part in seq_len(parts)) {
    sheet <- if (parts > 1L) paste0(name, "_", part) else name
    rows <- seq.int((part - 1L) * max_rows + 1L, min(part * max_rows, nrow(x)))
    sheets[[sheet]] <- as.data.frame(x[rows])
  }
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
workbook <- openxlsx::createWorkbook()
decimal_style <- openxlsx::createStyle(numFmt = "0.00000")
header_style <- openxlsx::createStyle(textDecoration = "bold")
for (name in names(sheets)) {
  x <- sheets[[name]]
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, x, headerStyle = header_style)
  decimal_columns <- which(vapply(x, is.double, logical(1)))
  if (length(decimal_columns)) {
    openxlsx::addStyle(workbook, name, decimal_style,
                       rows = seq_len(nrow(x)) + 1L, cols = decimal_columns,
                       gridExpand = TRUE)
  }
  # Allow space for headers, five-decimal numbers and full gene/variant IDs.
  widths <- vapply(seq_along(x), function(i) {
    values <- if (is.double(x[[i]])) sprintf("%.5f", x[[i]]) else as.character(x[[i]])
    max(20, nchar(names(x)[i]) + 2, nchar(values), na.rm = TRUE) + 2
  }, numeric(1))
  openxlsx::setColWidths(workbook, name, cols = seq_along(x), widths = widths)
}
openxlsx::saveWorkbook(workbook, file = output, overwrite = TRUE)
fwrite(full_variants, file = csv_output)
message("Created ", output)
for (name in names(sheets)) message(name, ": ", nrow(sheets[[name]]), " rows")
message("Created ", csv_output, ": ", nrow(full_variants), " rows")
