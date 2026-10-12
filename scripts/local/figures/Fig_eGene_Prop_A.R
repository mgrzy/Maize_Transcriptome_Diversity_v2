#!/usr/bin/env Rscript
# Run from the repository root:
# Rscript --vanilla scripts/local/figures/Fig_eGene_Prop_A.R
# Combine panel A of Fig_eGene_Prop.png with the eQTL discoverability
# plot from Fig S6 A on the right, using the same saved inputs.
# The upstream IHW workflow retains adj_p < 0.01; no new filtering is applied.
# "SNP" follows the original figure and includes all saved variant types.
# Percent of total genes uses the OSCA gene inventory, including failed probes
# with zero tests (not all genes in the reference annotation). The number of
# successfully tested genes and its corresponding percentage are also saved.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

association_file <- "results/summary/osca/eqtl_padjusted_full_results.tsv.gz"
tested_gene_file <- "results/summary/osca/cis_tested_variants_per_gene.tsv"
figure_file <- "results/figures/presentation/Fig_eGene_Prop_A.png"
summary_file <- "results/summary/osca/eGene_Prop_A_summary.tsv"
counts_file <- "results/summary/osca/eGene_Prop_A_snps_per_gene.tsv"
discoverability_file <- "results/full/osca/significant_summary/significant_gene_counts.tsv"

if (!file.exists(tested_gene_file)) {
  tested_gene_file <- paste0(tested_gene_file, ".gz")
}
for (path in c(association_file, tested_gene_file, discoverability_file)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

tested <- fread(tested_gene_file)
if (!all(c("gene", "n_cis_tested_variants") %in% names(tested))) {
  stop("Tested-gene table must contain gene and n_cis_tested_variants.")
}
if (anyNA(tested$gene) || any(!nzchar(tested$gene)) ||
    anyDuplicated(tested$gene) || anyNA(tested$n_cis_tested_variants) ||
    any(tested$n_cis_tested_variants < 0)) {
  stop("Invalid or duplicated identifiers/counts in tested-gene table.")
}
tested_genes <- tested[n_cis_tested_variants > 0, gene]
if (!length(tested_genes)) stop("No tested genes in denominator.")

# Read only the two necessary columns from the large association table.
message("Reading saved significant associations...")
associations <- fread(association_file, select = c("Gene", "SNP"),
                      colClasses = "character", nThread = 4L)
if (!nrow(associations) || anyNA(associations) ||
    any(!nzchar(associations$Gene)) || any(!nzchar(associations$SNP))) {
  stop("Association table is empty or contains missing gene/variant IDs.")
}
if (length(setdiff(unique(associations$Gene), tested_genes))) {
  stop("Some eGenes are absent from the tested-gene denominator.")
}

# Check uniqueness so row counts reproduce the original panel while also
# representing distinct SNPs per gene. Do not silently count duplicate hits.
message("Checking gene-SNP pairs and calculating statistics...")
if (anyDuplicated(associations, by = c("Gene", "SNP"))) {
  stop("Duplicate gene-SNP pairs found; resolve before reproducing panel A.")
}
snps_per_gene <- associations[, .(n_snps = .N), by = Gene]
setorder(snps_per_gene, Gene)
summary_stats <- data.table(
  n_eGenes = nrow(snps_per_gene),
  n_total_genes = nrow(tested),
  percent_total_genes = 100 * nrow(snps_per_gene) / nrow(tested),
  median_snps_per_eGene = median(snps_per_gene$n_snps),
  n_unique_snps = uniqueN(associations$SNP),
  n_total_tested_genes = length(tested_genes),
  percent_total_tested_genes = 100 * nrow(snps_per_gene) / length(tested_genes),
  n_gene_snp_associations = nrow(associations),
  denominator_definition = "All genes in OSCA inventory, including zero-test failed probes",
  association_input = association_file,
  denominator_input = tested_gene_file
)

panel_a <- ggplot(snps_per_gene, aes(x = n_snps)) +
  geom_histogram(binwidth = 1000, boundary = 0,
                 fill = "#66C5CC", colour = "black") +
  labs(x = "Number of associated SNPs per eGene", y = "Number of eGenes") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = rel(1.15)))

# Preserve the denominator and plotted values used in Fig S6 A.
discoverability <- fread(discoverability_file)
discoverability[, SignificantGenesProp := SignificantGenes / 15587]
panel_b <- ggplot(discoverability, aes(SampleSize, SignificantGenesProp)) +
  geom_point() +
  labs(x = "Sample size", y = "eQTL discoverability") +
  scale_x_continuous(breaks = discoverability$SampleSize,
                     labels = discoverability$SampleSize,
                     guide = guide_axis(check.overlap = TRUE)) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = rel(1.15)))
combined_figure <- panel_a | panel_b

dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, combined_figure, width = 18, height = 6, units = "cm", dpi = 300)
fwrite(snps_per_gene, counts_file, sep = "\t")
fwrite(summary_stats, summary_file, sep = "\t")
print(summary_stats[, 1:6])
message("Saved combined figure: ", figure_file,
        "\nSaved summary: ", summary_file,
        "\nSaved per-gene counts: ", counts_file)
