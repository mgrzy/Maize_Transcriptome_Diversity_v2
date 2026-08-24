### Script: analysis/9.Fst_gcta.R
### Purpose: Calculate per-SNP Fst with GCTA for the seven Limma DEG contrasts.
### Inputs:
###   - data/processed/geno/plink/Div639.{bed,bim,fam}
###   - data/raw/meta/WiDiv_meta.csv
### Outputs:
###   - results/full/Fst/<contrast>.fst
###   - results/full/Fst/<contrast>.{log,subpop}
###   - results/full/Fst/Fst_contrast_samples.csv
### Usage: Rscript scripts/local/analysis/9.Fst_gcta.R

bfile <- "data/processed/geno/plink/Div639"
metadata_path <- "data/raw/meta/WiDiv_meta.csv"
output_dir <- "results/full/Fst"
gcta_appimage <- Sys.getenv("GCTA_BIN", "../../soft/gcta/gcta")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ensure_gcta <- function(appimage) {
  appimage <- normalizePath(appimage, mustWork = TRUE)

  # Use GCTA directly when the AppImage can mount. In environments without
  # FUSE, extract its bundled executable into a temporary directory instead.
  direct_status <- suppressWarnings(
    system2(appimage, "--version", stdout = FALSE, stderr = FALSE)
  )
  if (identical(direct_status, 0L)) {
    return(list(executable = appimage, cleanup_dir = NULL))
  }

  extract_dir <- tempfile("gcta_appimage_")
  dir.create(extract_dir)
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  setwd(extract_dir)
  extract_status <- system2(appimage, "--appimage-extract")
  setwd(original_wd)

  executable <- file.path(extract_dir, "squashfs-root", "AppRun")
  if (!identical(extract_status, 0L) || !file.exists(executable)) {
    stop("Could not run or extract GCTA from ", appimage)
  }

  list(executable = executable, cleanup_dir = extract_dir)
}

run_fst <- function(gcta, fam_metadata, comparison_name, keep,
                    population, numerator, denominator) {
  selected <- which(keep %in% TRUE)
  subpop <- fam_metadata[selected, c("FID", "IID")]
  subpop$Population <- population[selected]
  group_sizes <- table(factor(subpop$Population, levels = c(denominator, numerator)))

  if (any(group_sizes < 2L)) {
    stop(
      comparison_name, ": each population needs at least two lines (",
      paste(names(group_sizes), group_sizes, sep = "=", collapse = ", "), ")"
    )
  }

  subpop_file <- file.path(output_dir, paste0(comparison_name, ".subpop"))
  output_prefix <- file.path(output_dir, comparison_name)
  write.table(
    subpop,
    subpop_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  gcta_args <- c(
    "--bfile", bfile,
    "--fst",
    "--sub-popu", subpop_file,
    "--out", output_prefix
  )
  status <- system2(gcta, gcta_args)
  if (!identical(status, 0L)) {
    stop("GCTA failed for ", comparison_name)
  }

  data.frame(
    contrast = comparison_name,
    population_1 = denominator,
    n_population_1 = unname(group_sizes[[denominator]]),
    population_2 = numerator,
    n_population_2 = unname(group_sizes[[numerator]]),
    fst_file = paste0(output_prefix, ".fst"),
    stringsAsFactors = FALSE
  )
}

required_bfile_files <- paste0(bfile, c(".bed", ".bim", ".fam"))
if (!all(file.exists(required_bfile_files))) {
  stop("Missing PLINK files: ", paste(required_bfile_files[!file.exists(required_bfile_files)], collapse = ", "))
}

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata$Era[metadata$Era == "Public Sector"] <- "Public"
metadata$Era[metadata$Era == ""] <- NA_character_

fam <- read.table(
  paste0(bfile, ".fam"),
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(fam)[1:2] <- c("FID", "IID")
metadata_index <- match(fam$IID, metadata$Taxa)
matched <- !is.na(metadata_index)
fam_metadata <- cbind(
  fam[matched, c("FID", "IID")],
  metadata[metadata_index[matched], c("Taxa", "Group", "Era")]
)

comparison_specs <- list(
  SS_vs_Tropical = list(
    keep = fam_metadata$Group %in% c("SS", "Tropical"),
    population = fam_metadata$Group,
    numerator = "SS", denominator = "Tropical"
  ),
  NSS_vs_Tropical = list(
    keep = fam_metadata$Group %in% c("NSS", "Tropical"),
    population = fam_metadata$Group,
    numerator = "NSS", denominator = "Tropical"
  ),
  Iodent_vs_Tropical = list(
    keep = fam_metadata$Group %in% c("IDT", "Tropical"),
    population = ifelse(fam_metadata$Group == "IDT", "IDT", "Tropical"),
    numerator = "IDT", denominator = "Tropical"
  ),
  SS_exPVP_vs_Public = list(
    keep = fam_metadata$Group == "SS" & fam_metadata$Era %in% c("exPVP", "Public"),
    population = fam_metadata$Era,
    numerator = "exPVP", denominator = "Public"
  ),
  NSS_exPVP_vs_Public = list(
    keep = fam_metadata$Group == "NSS" & fam_metadata$Era %in% c("exPVP", "Public"),
    population = fam_metadata$Era,
    numerator = "exPVP", denominator = "Public"
  ),
  Public_SS_vs_NSS = list(
    keep = fam_metadata$Group %in% c("SS", "NSS") & fam_metadata$Era == "Public",
    population = fam_metadata$Group,
    numerator = "SS", denominator = "NSS"
  ),
  exPVP_SS_vs_NSS = list(
    keep = fam_metadata$Group %in% c("SS", "NSS") & fam_metadata$Era == "exPVP",
    population = fam_metadata$Group,
    numerator = "SS", denominator = "NSS"
  )
)

gcta <- ensure_gcta(gcta_appimage)
if (!is.null(gcta$cleanup_dir)) {
  on.exit(unlink(gcta$cleanup_dir, recursive = TRUE), add = TRUE)
}

contrast_summary <- do.call(
  rbind,
  lapply(names(comparison_specs), function(comparison_name) {
    spec <- comparison_specs[[comparison_name]]
    run_fst(
      gcta = gcta$executable,
      fam_metadata = fam_metadata,
      comparison_name = comparison_name,
      keep = spec$keep,
      population = spec$population,
      numerator = spec$numerator,
      denominator = spec$denominator
    )
  })
)

write.csv(
  contrast_summary,
  file.path(output_dir, "Fst_contrast_samples.csv"),
  row.names = FALSE,
  quote = FALSE
)
message("Wrote Fst results to ", output_dir)
