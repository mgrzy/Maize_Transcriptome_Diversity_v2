#!/usr/bin/env Rscript
# Supplementary figure: per-TF ChIP-seq enrichment for all 104 Tu/Mejia-Guerra
# TFs, grouped by the same eight functional paper classes used in Figure 3.
# Standard matched controls and leaf-ATAC-matched controls are shown together.
#
# Optional environment variables:
#   TF_ENRICHMENT_FILE
#   TF_ATAC_ENRICHMENT_FILE
#   TF_CLASSIFICATION_FILE
#   FIGURE_OUTPUT_DIR
#   FIGURE_STEM

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

stop_with <- function(...) stop(sprintf(...), call. = FALSE)

cfg <- list(
  standard_file = Sys.getenv(
    "TF_ENRICHMENT_FILE",
    "results/summary/fine_mapped_TF_ChIP_enrichment/TF_overlap_matched_permutation.tsv"
  ),
  atac_file = Sys.getenv(
    "TF_ATAC_ENRICHMENT_FILE",
    "results/summary/fine_mapped_TF_ChIP_enrichment_leaf_ATAC_matched_test/TF_overlap_matched_permutation.tsv"
  ),
  classification_file = Sys.getenv(
    "TF_CLASSIFICATION_FILE",
    "data/processed/chromatin/Tu_2020_TF_functional_classification.xlsx"
  ),
  output_dir = Sys.getenv("FIGURE_OUTPUT_DIR", "results/figures"),
  figure_stem = Sys.getenv("FIGURE_STEM", "FigS_TF_enrichment_all_104_by_functional_class")
)

for (path in unlist(cfg[1:3])) {
  if (!file.exists(path)) stop_with("Required input does not exist: %s", path)
}
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

read_enrichment <- function(path, label) {
  x <- fread(path)
  required <- c("TF", "fold_enrichment")
  missing_columns <- setdiff(required, names(x))
  if (length(missing_columns)) {
    stop_with("%s lacks required column(s): %s", label, paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(x$TF)) stop_with("%s has duplicated TF identifiers.", label)
  x <- x[, .(TF = as.character(TF), fold_enrichment = as.numeric(fold_enrichment))]
  if (x[, any(!is.finite(fold_enrichment) | fold_enrichment <= 0)]) {
    stop_with("%s has non-finite or non-positive fold-enrichment values.", label)
  }
  x
}

standard <- read_enrichment(cfg$standard_file, "Standard-control enrichment input")
atac <- read_enrichment(cfg$atac_file, "Leaf-ATAC-matched enrichment input")
if (!setequal(standard$TF, atac$TF)) {
  stop_with("The standard and leaf-ATAC-matched inputs do not contain the same TFs.")
}

classification <- as.data.table(read_excel(cfg$classification_file, sheet = "TF_classification"))
required_class_columns <- c("TF", "paper_class")
missing_class_columns <- setdiff(required_class_columns, names(classification))
if (length(missing_class_columns)) {
  stop_with("Classification workbook lacks required column(s): %s", paste(missing_class_columns, collapse = ", "))
}
classification <- classification[, .(
  TF = as.character(TF),
  paper_class = trimws(as.character(paper_class))
)]
classification <- classification[!is.na(TF) & nzchar(TF) & !is.na(paper_class) & nzchar(paper_class)]
if (anyDuplicated(classification$TF)) stop_with("Classification workbook has duplicated TF identifiers.")
if (!setequal(standard$TF, classification$TF)) {
  missing_class <- setdiff(standard$TF, classification$TF)
  extra_class <- setdiff(classification$TF, standard$TF)
  stop_with(
    "Classification and enrichment TF sets differ (missing classifications: %s; extra classifications: %s).",
    paste(missing_class, collapse = ", "), paste(extra_class, collapse = ", ")
  )
}

class_labels <- c(
  "PS" = "Photosynthesis",
  "cell wall" = "Cell wall",
  "hormone" = "Hormone",
  "metabolism" = "Metabolism",
  "other" = "Other",
  "redox" = "Redox",
  "signalling" = "Signalling",
  "transport" = "Transport"
)
if (!setequal(unique(classification$paper_class), names(class_labels))) {
  stop_with("Classification workbook does not contain exactly the eight Figure 3 paper classes.")
}

tf <- merge(
  merge(standard, atac, by = "TF", suffixes = c("_standard", "_atac_matched")),
  classification,
  by = "TF",
  all.x = TRUE,
  sort = FALSE
)
tf[, class_label := unname(class_labels[paper_class])]

# Class order and panel arrangement follow the Figure 3 labels, with the
# large signalling class paired with photosynthesis to keep the figure compact.
class_order <- c("PS", "signalling", "hormone", "cell wall", "metabolism", "transport", "redox", "other")
tf[, paper_class := factor(paper_class, levels = class_order)]
setorder(tf, paper_class, fold_enrichment_standard, TF)

long_tf <- melt(
  tf,
  id.vars = c("TF", "paper_class", "class_label"),
  measure.vars = c("fold_enrichment_standard", "fold_enrichment_atac_matched"),
  variable.name = "control_set",
  value.name = "fold_enrichment"
)
long_tf[, control_set := factor(
  control_set,
  levels = c("fold_enrichment_standard", "fold_enrichment_atac_matched"),
  labels = c("Standard matched controls", "Leaf-ATAC-matched controls")
)]

x_limit <- max(long_tf$fold_enrichment) * 1.08
base_font_size <- 9
panel_theme <- theme_classic(base_size = base_font_size) +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(colour = "black"),
    # The 104 TF labels are necessarily more compact than the shared 9 pt
    # figure text; 6.5 pt remains legible at the manuscript export size.
    axis.text.y = element_text(size = 6.5),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

make_class_panel <- function(class_name) {
  class_data <- long_tf[paper_class == class_name]
  tf_order <- tf[paper_class == class_name][order(fold_enrichment_standard, TF), TF]
  class_data[, TF_plot := factor(TF, levels = tf_order)]
  label <- unique(class_data$class_label)
  n_tfs <- uniqueN(class_data$TF)

  ggplot(class_data, aes(y = TF_plot)) +
    geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey45") +
    geom_segment(
      data = tf[paper_class == class_name],
      aes(x = fold_enrichment_standard, xend = fold_enrichment_atac_matched,
          y = factor(TF, levels = tf_order), yend = factor(TF, levels = tf_order)),
      inherit.aes = FALSE, linewidth = 0.45, colour = "grey65"
    ) +
    geom_point(aes(x = fold_enrichment, fill = control_set),
      shape = 21, size = 2.05, stroke = 0.38, colour = "black") +
    scale_fill_manual(values = c(
      "Standard matched controls" = "white",
      "Leaf-ATAC-matched controls" = "#0072B2"
    )) +
    scale_x_continuous(
      limits = c(0.75, x_limit),
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(title = sprintf("%s (n = %d)", label, n_tfs), x = "Fold enrichment", fill = NULL) +
    coord_cartesian(clip = "off") +
    panel_theme
}

plots <- lapply(class_order, make_class_panel)
figure <- wrap_plots(plots, ncol = 2, byrow = TRUE, guides = "collect") +
  plot_layout(heights = c(40, 18, 10, 4)) &
  theme(legend.position = "top", legend.text = element_text(size = base_font_size))

png_file <- file.path(cfg$output_dir, paste0(cfg$figure_stem, ".png"))
pdf_file <- file.path(cfg$output_dir, paste0(cfg$figure_stem, ".pdf"))
ggsave(png_file, figure, width = 12, height = 16, units = "in", dpi = 400, bg = "white")
ggsave(pdf_file, figure, width = 12, height = 16, units = "in", device = cairo_pdf, bg = "white")

output_table <- tf[, .(
  TF, paper_class, functional_class = class_label,
  fold_enrichment_standard, fold_enrichment_atac_matched
)]
fwrite(output_table,
  file.path(cfg$output_dir, paste0(cfg$figure_stem, "_data.tsv")), sep = "\t"
)

message("Wrote supplementary all-TF enrichment figure to: ", png_file)
