### Script: scripts/local/procesing/11.DowSampling.R
### Purpose: Create reproducible, approximately group-proportional subsets of
###          the PCA-matched lines for PLINK downsampling analyses.
### Inputs:
###  - data/raw/meta/TableS1.csv
###  - results/full/PCA/PCA_Geno.eigenvec
###  - data/processed/geno/plink/Div639.fam
### Outputs:
###  - data/processed/meta/downsample/DownSampling_<n>.fam
###  - data/processed/meta/downsample/DownSampling_group_counts.tsv
### Usage: Rscript scripts/local/procesing/11.DowSampling.R
### Dependencies: tidyverse

library(tidyverse)

# Resolve project-relative paths even when the script is run outside the
# repository root.
script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))
])
if (length(script_file) == 1L) {
  setwd(normalizePath(file.path(dirname(script_file), "../../..")))
}

set.seed(20260722)

dat <- read.csv("data/raw/meta/TableS1.csv", stringsAsFactors = FALSE)
pca.snp <- read.table("results/full/PCA/PCA_Geno.eigenvec", stringsAsFactors = FALSE)

# Keep the same PCA-matched lines used to generate the genetic PCA.
dat2 <- dat[as.character(dat$ID) %in% as.character(pca.snp$V1), , drop = FALSE]
dat2$ID <- as.character(dat2$ID)

fam <- read.table(
  "data/processed/geno/plink/Div639.fam",
  colClasses = "character",
  stringsAsFactors = FALSE
)

if (!all(dat2$ID %in% fam$V2)) {
  missing_ids <- setdiff(dat2$ID, fam$V2)
  stop("The following dat2 IDs are absent from Div639.fam: ",
       paste(missing_ids, collapse = ", "))
}

# Allocate each requested total to Groups using proportional allocation and
# largest-remainder rounding.  The minimum of one line keeps every Group
# represented, including the small Landrace group.
allocate_by_group <- function(n, group_counts) {
  expected <- n * group_counts / sum(group_counts)
  allocation <- pmax(1L, floor(expected))

  while (sum(allocation) < n) {
    candidates <- which(allocation < group_counts)
    score <- expected[candidates] - floor(expected[candidates])
    allocation[candidates[order(score, decreasing = TRUE)[1]]] <-
      allocation[candidates[order(score, decreasing = TRUE)[1]]] + 1L
  }

  while (sum(allocation) > n) {
    candidates <- which(allocation > 1L)
    score <- expected[candidates] - floor(expected[candidates])
    allocation[candidates[order(score, decreasing = FALSE)[1]]] <-
      allocation[candidates[order(score, decreasing = FALSE)[1]]] - 1L
  }

  allocation
}

group_counts <- table(dat2$Group)
targets <- c(100L, 150L, 200L, 300L, 400L, 500L, 550L, 600L)
out_dir <- "data/processed/meta/downsample"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

allocation_summary <- vector("list", length(targets))

for (target in targets) {
  allocation <- allocate_by_group(target, group_counts)
  names(allocation) <- names(group_counts)

  selected <- unlist(lapply(names(allocation), function(group) {
    group_ids <- dat2$ID[dat2$Group == group]
    sample(group_ids, size = allocation[[group]], replace = FALSE)
  }), use.names = FALSE)

  # Write the complete six-column PLINK .fam rows in the original genotype
  # order, rather than creating a two-column keep file with altered IDs.
  subset_fam <- fam[fam$V2 %in% selected, , drop = FALSE]
  if (nrow(subset_fam) != target) {
    stop("Expected ", target, " lines but selected ", nrow(subset_fam))
  }

  write.table(
    subset_fam,
    file = file.path(out_dir, paste0("DownSampling_", target, ".fam")),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  allocation_summary[[as.character(match(target, targets))]] <-
    data.frame(N = target, Group = names(allocation), Lines = as.integer(allocation))
}

write.table(
  bind_rows(allocation_summary),
  file.path(out_dir, "DownSampling_group_counts.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)
