# Maize Transcriptome Diversity

Code and lightweight documentation for the **Maize Transcriptome Diversity** project, focused on natural variation in gene expression and its genetic and evolutionary basis in the WiDiv maize diversity panel.

The repository contains analysis, figure-generation, preprocessing, and HPC workflow code. Large genomic, RNA-seq, and intermediate result files are intentionally kept outside Git.

## Project scope

The main analysis includes:

- transcriptomic diversity and population structure in the WiDiv panel;
- cis-eQTL mapping with OSCA mixed models;
- multiple-testing correction using IHW with distance to TSS;
- fine-mapping with SuSiE / SuSiE-RSS;
- allelic effect-size estimation (aFC-n);
- functional annotation and enrichment of fine-mapped regulatory variants;
- chromatin, ATAC-seq, and TF ChIP overlap analyses;
- population-genetic analyses including nucleotide diversity and FST;
- sensitivity analyses such as sample-size downsampling;
- manuscript and supplementary figure generation.

## Repository structure

```text
.
├── AGENTIC.md                         # guidance for AI/code agents working on the project
├── README.md                          # this file
├── Maize_Transcriptome_Diversity_v2.Rproj
├── data/
│   ├── raw/                           # local/raw inputs; large files ignored by Git
│   └── processed/                     # derived inputs; large files ignored by Git
├── notes/
│   └── Summary.Rmd                    # project notes / analysis summary
├── scripts/
│   ├── hpc/                           # HPC and SLURM workflows
│   │   ├── pop/                       # population-genetic workflows
│   │   └── rest/                      # download / auxiliary HPC workflows
│   ├── local/
│   │   ├── processing/                # preprocessing and dataset preparation
│   │   ├── analysis/
│   │   │   ├── expression/            # expression diversity, CCA and heritability
│   │   │   ├── eqtl/                  # OSCA mapping, significance and IHW
│   │   │   ├── fine_mapping/          # SuSiE, aFC and coloc
│   │   │   ├── differential_expression/ # limma DE and summaries
│   │   │   ├── functional_annotation/ # genomic context, GERP/CAD, ATAC/TF enrichment
│   │   │   ├── population_genetics/   # pi, genotype PCA, FST and DEG/eQTL integration
│   │   │   ├── downsampling/          # sample-size sensitivity analysis
│   │   │   └── todo/                  # historical/incomplete analysis scripts
│   │   ├── figures/                   # main and supplementary figure scripts
│   │   └── prepare_Tu_TF_consensus_peaks.sh
│   └── utils/
│       └── functions.R                # reusable R helper functions
└── sync.sh                            # conservative Git sync helper
```

Generated results are generally written under `results/` in the working copy. Large results and intermediate files should not be committed.

## Analysis environments

### Local workstation

Used for R/RStudio analyses, plotting, integration of summary files, and small-to-medium preprocessing jobs.

### HPC

Used for compute-intensive workflows such as genotype processing, OSCA/eQTL analyses, population-genetic calculations, and large batch operations. HPC scripts are under `scripts/hpc/`.

## Conventions

- Run scripts from the repository root unless a script explicitly documents otherwise.
- Prefer project-relative paths such as `data/...` and `results/...`.
- Do not commit FASTQ, BAM/CRAM, VCF/BCF, large BED files, RDS objects, or full analysis outputs.
- Keep reusable functions in `scripts/utils/functions.R` rather than copying helpers between scripts.
- New local data-preparation scripts belong in `scripts/local/processing/`.
- New analysis scripts should go into the relevant domain folder under `scripts/local/analysis/` rather than the analysis root.
- Manuscript figure scripts belong in `scripts/local/figures/`.
- HPC production jobs belong in `scripts/hpc/`.

See `scripts/local/analysis/README.md` for the analysis-folder map and execution-order notes.

To export supplementary Tables S1–S3 into one Excel workbook, install the R
packages `data.table` and `openxlsx`, then run from the repository root:

```bash
Rscript scripts/local/processing/Export_Supplementary_Tables.R
```

The output is `results/tables/Supplementary_Tables_S1_S3.xlsx`; an optional first
argument selects another output path. S1 contains saved CCA/PCA scores with
sample metadata. S2 lists one saved lead variant per credible set with its PIP.
S3 contains the saved aFCn results directly. The full list of fine-mapped
variant/credible-set memberships and PIPs is exported separately to
`results/tables/Supplementary_Tables_S1_S3_Full_fine_mapped_variants.csv`,
without eGene filtering. With a custom workbook path, the CSV uses the same
directory and filename stem plus `_Full_fine_mapped_variants.csv`.
The workbook contains only result sheets, with widened columns and continuous
numeric results displayed to five decimal places. Both exports round continuous
numeric results to five decimal places; source analysis files retain their
original precision. No analyses are refitted by this export.

## Reproducibility

The `.gitignore` is intentionally strict for large genomics files and generated results. If package reproducibility is managed with `renv`, commit `renv.lock`, `renv/activate.R`, and `renv/settings.json`, but never the local `renv/library/` cache.

For detailed instructions to coding agents and future automated maintenance, see [`AGENTIC.md`](AGENTIC.md).
