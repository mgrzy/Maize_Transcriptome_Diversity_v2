#!/usr/bin/env Rscript
# Paper Figure 3: enrichment of fine-mapped cis-eQTL variants across genomic,
# chromatin-accessibility, and TF paper-class annotations. Interpretive detail
# belongs in the figure caption; consequently, the plot has no titles,
# subtitles, legends, point-value labels, or in-plot captions.
#
# Inputs (overridable with environment variables):
#   GENIC_ENRICHMENT_FILE
#   ATAC_ENRICHMENT_FILE
#   PAPER_CLASS_ENRICHMENT_FILE
#   PAPER_CLASS_ATAC_ENRICHMENT_FILE
#
# Outputs (default: results/figures/Fig3.{pdf,png}):
#   FIGURE_OUTPUT_DIR and FIGURE_STEM control the destination.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  genic_file = Sys.getenv(
    "GENIC_ENRICHMENT_FILE",
    "results/summary/eQTL_genomic_context_enrichment/matched_enrichment_summary.tsv"
  ),
  atac_file = Sys.getenv(
    "ATAC_ENRICHMENT_FILE",
    "results/summary/fine_mapped_chromatin_matched_test/chromatin_overlap_matched_permutation.tsv"
  ),
  paper_class_file = Sys.getenv(
    "PAPER_CLASS_ENRICHMENT_FILE",
    "results/summary/fine_mapped_TF_paper_class_enrichment/paper_class_overlap_matched_permutation.tsv"
  ),
  paper_class_atac_file = Sys.getenv(
    "PAPER_CLASS_ATAC_ENRICHMENT_FILE",
    "results/summary/fine_mapped_TF_paper_class_enrichment_leaf_ATAC_matched/paper_class_overlap_matched_permutation.tsv"
  ),
  genic_annotation_file = Sys.getenv(
    "GENIC_ANNOTATION_FILE",
    "results/summary/eQTL_genomic_context/fine_mapped_variant_annotations.tsv.gz"
  ),
  atac_annotation_file = Sys.getenv(
    "ATAC_ANNOTATION_FILE",
    "results/summary/fine_mapped_chromatin_matched_test/fine_mapped_variant_chromatin_annotations.tsv.gz"
  ),
  paper_class_annotation_file = Sys.getenv(
    "PAPER_CLASS_ANNOTATION_FILE",
    "results/summary/fine_mapped_TF_paper_class_enrichment/fine_mapped_variant_paper_class_annotations.tsv.gz"
  ),
  paper_class_atac_annotation_file = Sys.getenv(
    "PAPER_CLASS_ATAC_ANNOTATION_FILE",
    "results/summary/fine_mapped_TF_paper_class_enrichment_leaf_ATAC_matched/fine_mapped_variant_paper_class_annotations.tsv.gz"
  ),
  output_dir = Sys.getenv("FIGURE_OUTPUT_DIR", "results/figures"),
  figure_stem = Sys.getenv("FIGURE_STEM", "Fig3")
)

for (path in unlist(cfg[1:8])) {
  if (!file.exists(path)) stop_with("Required enrichment result does not exist: %s", path)
}
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

require_columns <- function(x, columns, label) {
  missing_columns <- setdiff(columns, names(x))
  if (length(missing_columns)) {
    stop_with("%s is missing required column(s): %s", label, paste(missing_columns, collapse = ", "))
  }
}

make_ordered_factor <- function(x, label_column, value_column) {
  levels <- x[order(get(value_column), get(label_column)), get(label_column)]
  factor(x[[label_column]], levels = levels)
}

annotation_variant_column <- function(x, label) {
  candidates <- c("variant_id", "SNP")
  column <- candidates[candidates %in% names(x)][1L]
  if (is.na(column) || !nzchar(column)) {
    stop_with("%s lacks a variant identifier column (expected variant_id or SNP).", label)
  }
  column
}

summarise_exclusive_proportions <- function(x, variant_column, category_column, categories, label) {
  x <- x[, .(variant_id = as.character(get(variant_column)), category = as.character(get(category_column)))]
  x <- x[!is.na(variant_id) & nzchar(variant_id)]
  if (x[, anyDuplicated(variant_id)]) {
    inconsistent <- x[, uniqueN(category) > 1L, by = variant_id][V1 == TRUE]
    if (nrow(inconsistent)) {
      stop_with("%s assigns %d fine-mapped variants to more than one exclusive category.", label, nrow(inconsistent))
    }
    x <- unique(x, by = "variant_id")
  }
  if (!setequal(unique(x$category), categories)) {
    stop_with("%s categories do not match the enrichment categories.", label)
  }
  counts <- x[, .(n_category_variants = .N), by = category]
  counts <- data.table(category = categories)[counts, on = "category"]
  counts[is.na(n_category_variants), n_category_variants := 0L]
  counts[, `:=`(
    n_unique_variants = uniqueN(x$variant_id),
    proportion = n_category_variants / uniqueN(x$variant_id)
  )]
  counts
}

summarise_overlapping_proportions <- function(x, variant_column, feature_column, categories, label) {
  x <- x[, .(
    variant_id = as.character(get(variant_column)),
    features = fifelse(is.na(get(feature_column)), "", as.character(get(feature_column)))
  )]
  x <- x[!is.na(variant_id) & nzchar(variant_id)]
  variants <- unique(x$variant_id)
  membership <- rbindlist(lapply(categories, function(category) {
    matching_variants <- unique(x[
      vapply(strsplit(features, ";", fixed = TRUE), function(values) category %chin% values, logical(1L)),
      variant_id
    ])
    data.table(category = category, n_category_variants = length(matching_variants))
  }))
  membership[, `:=`(
    n_unique_variants = length(variants),
    proportion = n_category_variants / length(variants)
  )]
  membership
}

align_proportions <- function(enrichment, proportions) {
  aligned <- merge(
    data.table(annotation = as.character(enrichment$annotation)),
    proportions,
    by.x = "annotation",
    by.y = "category",
    all.x = TRUE,
    sort = FALSE
  )
  if (aligned[, anyNA(proportion)]) stop_with("Missing fine-mapped-variant proportions for plotted category labels.")
  aligned[, annotation := factor(annotation, levels = levels(enrichment$annotation))]
  aligned
}

genic_labels <- c(
  "5UTR" = "5' UTR",
  "3UTR" = "3' UTR",
  "CDS" = "CDS",
  "promoter" = "Promoter",
  "intron" = "Intron",
  "intergenic" = "Intergenic"
)
atac_labels <- c(
  "ATAC_ear_specific" = "Ear-specific",
  "ATAC_leaf_specific" = "Leaf-specific",
  "ATAC_neither_leaf_nor_ear" = "Neither leaf nor ear",
  "ATAC_shared_leaf_ear" = "Shared leaf and ear"
)
paper_class_labels <- c(
  "PS" = "Photosynthesis",
  "cell wall" = "Cell wall",
  "hormone" = "Hormone",
  "metabolism" = "Metabolism",
  "other" = "Other",
  "redox" = "Redox",
  "signalling" = "Signalling",
  "transport" = "Transport"
)

genic <- fread(cfg$genic_file)
require_columns(genic, c("primary_annotation", "fold_enrichment"), "Genic enrichment input")
genic <- genic[
  primary_annotation %chin% names(genic_labels) & is.finite(fold_enrichment),
  .(
    annotation = unname(genic_labels[primary_annotation]),
    fold_enrichment = as.numeric(fold_enrichment)
  )
]
if (nrow(genic) != length(genic_labels)) {
  stop_with("Genic enrichment input does not contain all expected finite categories.")
}
genic[, annotation := make_ordered_factor(genic, "annotation", "fold_enrichment")]

atac <- fread(cfg$atac_file)
require_columns(atac, c("chromatin_feature", "fold_enrichment"), "ATAC enrichment input")
atac <- atac[
  chromatin_feature %chin% names(atac_labels) & is.finite(fold_enrichment),
  .(
    annotation = unname(atac_labels[chromatin_feature]),
    fold_enrichment = as.numeric(fold_enrichment)
  )
]
if (nrow(atac) != length(atac_labels)) {
  stop_with("ATAC enrichment input does not contain all four mutually exclusive ATAC categories.")
}
atac[, annotation := make_ordered_factor(atac, "annotation", "fold_enrichment")]

paper <- fread(cfg$paper_class_file)
paper_atac <- fread(cfg$paper_class_atac_file)
require_columns(paper, c("paper_class", "fold_enrichment"), "Paper-class enrichment input")
require_columns(paper_atac, c("paper_class", "fold_enrichment"), "ATAC-matched paper-class enrichment input")
if (anyDuplicated(paper$paper_class) || anyDuplicated(paper_atac$paper_class)) {
  stop_with("Paper-class enrichment inputs must contain one row per paper_class.")
}

tf_class <- merge(
  paper[, .(paper_class, fold_enrichment_standard = as.numeric(fold_enrichment))],
  paper_atac[, .(paper_class, fold_enrichment_atac_matched = as.numeric(fold_enrichment))],
  by = "paper_class",
  all = TRUE,
  sort = FALSE
)
if (
  nrow(tf_class) != length(paper_class_labels) ||
  anyNA(tf_class$fold_enrichment_standard) || anyNA(tf_class$fold_enrichment_atac_matched) ||
  !setequal(tf_class$paper_class, names(paper_class_labels))
) {
  stop_with("Paper-class enrichment results must contain the same eight complete paper classes.")
}
tf_class[, annotation := unname(paper_class_labels[paper_class])]
tf_class_order <- c(
  "Other",
  tf_class[annotation != "Other"][order(fold_enrichment_standard), annotation]
)
tf_class[, annotation := factor(annotation, levels = tf_class_order)]

# Fine-mapped-variant proportions -----------------------------------------

genic_annotation <- fread(cfg$genic_annotation_file)
genic_variant_column <- annotation_variant_column(genic_annotation, "Genic annotation input")
require_columns(genic_annotation, "primary_annotation", "Genic annotation input")
genic_proportions <- summarise_exclusive_proportions(
  genic_annotation,
  genic_variant_column,
  "primary_annotation",
  names(genic_labels),
  "Genic annotation input"
)
genic_proportions[, category := unname(genic_labels[category])]
genic_proportions <- align_proportions(genic, genic_proportions)

atac_annotation <- fread(cfg$atac_annotation_file)
atac_variant_column <- annotation_variant_column(atac_annotation, "ATAC annotation input")
require_columns(atac_annotation, "chromatin_features", "ATAC annotation input")
atac_variant_status <- atac_annotation[, .(
  has_leaf = any(grepl("leaf", chromatin_features, ignore.case = TRUE), na.rm = TRUE),
  has_ear = any(grepl("ear", chromatin_features, ignore.case = TRUE), na.rm = TRUE)
), by = .(variant_id = get(atac_variant_column))]
atac_variant_status[, category := fcase(
  has_leaf & has_ear, "ATAC_shared_leaf_ear",
  has_leaf, "ATAC_leaf_specific",
  has_ear, "ATAC_ear_specific",
  default = "ATAC_neither_leaf_nor_ear"
)]
atac_proportions <- summarise_exclusive_proportions(
  atac_variant_status,
  "variant_id",
  "category",
  names(atac_labels),
  "ATAC annotation input"
)
atac_proportions[, category := unname(atac_labels[category])]
atac_proportions <- align_proportions(atac, atac_proportions)

paper_annotation <- fread(cfg$paper_class_annotation_file)
paper_atac_annotation <- fread(cfg$paper_class_atac_annotation_file)
paper_variant_column <- annotation_variant_column(paper_annotation, "Paper-class annotation input")
paper_atac_variant_column <- annotation_variant_column(paper_atac_annotation, "ATAC-matched paper-class annotation input")
require_columns(paper_annotation, "chromatin_features", "Paper-class annotation input")
require_columns(paper_atac_annotation, "chromatin_features", "ATAC-matched paper-class annotation input")
paper_proportions <- summarise_overlapping_proportions(
  paper_annotation,
  paper_variant_column,
  "chromatin_features",
  names(paper_class_labels),
  "Paper-class annotation input"
)
paper_atac_proportions <- summarise_overlapping_proportions(
  paper_atac_annotation,
  paper_atac_variant_column,
  "chromatin_features",
  names(paper_class_labels),
  "ATAC-matched paper-class annotation input"
)
setorder(paper_proportions, category)
setorder(paper_atac_proportions, category)
if (!isTRUE(all.equal(paper_proportions, paper_atac_proportions, check.attributes = FALSE))) {
  warning("Standard and leaf-ATAC-matched paper-class annotation proportions differ; plotting the standard annotation proportions.")
}
paper_proportions[, category := unname(paper_class_labels[category])]
paper_proportions <- align_proportions(tf_class, paper_proportions)

proportion_output <- rbindlist(list(
  genic_proportions[, .(panel = "Genic", annotation, n_unique_variants, n_category_variants, proportion)],
  atac_proportions[, .(panel = "ATAC", annotation, n_unique_variants, n_category_variants, proportion)],
  paper_proportions[, .(panel = "TF paper class", annotation, n_unique_variants, n_category_variants, proportion)]
))
fwrite(
  proportion_output,
  file.path(cfg$output_dir, paste0(cfg$figure_stem, "_category_proportions.tsv")),
  sep = "\t"
)

panel_theme <- theme_classic(base_size = 9.5) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.text.y = element_text(size = 10.5, colour = "black"),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

make_lollipop_panel <- function(data, point_fill) {
  ggplot(data, aes(x = fold_enrichment, y = annotation)) +
    geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = "grey45") +
    geom_segment(aes(x = 1, xend = fold_enrichment, yend = annotation), linewidth = 0.55, colour = "grey70") +
    geom_point(shape = 21, size = 3, stroke = 0.55, fill = point_fill, colour = "black") +
    scale_x_continuous(
      limits = c(0.85, NA),
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(x = "Fold enrichment") +
    coord_cartesian(clip = "off") +
    panel_theme
}

make_proportion_panel <- function(data, bar_fill) {
  ggplot(data, aes(x = proportion, y = annotation)) +
    geom_col(width = 0.62, fill = bar_fill, colour = "black", linewidth = 0.25) +
    scale_x_continuous(
      limits = c(0, NA),
      labels = scales::label_percent(accuracy = 1),
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(x = "Proportion of\neQTLs") +
    panel_theme +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 8, lineheight = 0.9),
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )
}

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0.02, 0.98)
)

genic_plot <- make_lollipop_panel(genic, "#66C5CC") + labs(tag = "A") + tag_theme
atac_plot <- make_lollipop_panel(atac, "#F6CF71") + labs(tag = "B") + tag_theme

tf_class_plot <- ggplot(tf_class, aes(y = annotation)) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = "grey45") +
  geom_segment(
    aes(x = fold_enrichment_standard, xend = fold_enrichment_atac_matched, yend = annotation),
    linewidth = 0.65,
    colour = "grey65"
  ) +
  geom_point(
    aes(x = fold_enrichment_standard, fill = "Standard matched\ncontrols"),
    shape = 21, size = 3, stroke = 0.55, colour = "black"
  ) +
  geom_point(
    aes(x = fold_enrichment_atac_matched, fill = "Leaf-ATAC-matched\ncontrols"),
    shape = 21, size = 3, stroke = 0.55, colour = "black"
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Standard matched\ncontrols" = "white",
      "Leaf-ATAC-matched\ncontrols" = "#0072B2"
    )
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Fold enrichment") +
  coord_cartesian(clip = "off") +
  panel_theme +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0), colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.text = element_text(size = 7.5)
  ) +
  labs(tag = "C") +
  tag_theme

genic_proportion_plot <- make_proportion_panel(genic_proportions, "#66C5CC")
atac_proportion_plot <- make_proportion_panel(atac_proportions, "#F6CF71")
tf_class_proportion_plot <- make_proportion_panel(paper_proportions, "#0072B2")

genic_panel <- (genic_plot | genic_proportion_plot) + plot_layout(widths = c(2.5, 1))
atac_panel <- (atac_plot | atac_proportion_plot) + plot_layout(widths = c(2.5, 1))
tf_class_panel <- (tf_class_plot | tf_class_proportion_plot) + plot_layout(widths = c(2.5, 1))

figure <- wrap_plots(
  genic_plot,
  genic_proportion_plot,
  atac_plot,
  atac_proportion_plot,
  tf_class_plot,
  tf_class_proportion_plot,
  nrow = 1,
  widths = c(2.5, 1, 2.5, 1, 2.5, 1)
)

print(figure)

ggsave(
  file.path(cfg$output_dir, paste0(cfg$figure_stem, ".png")),
  figure,
  width = 16,
  height = 5.2,
  units = "in",
  dpi = 400
)

message("Wrote Figure 3 to: ", file.path(cfg$output_dir, paste0(cfg$figure_stem, ".png")))
