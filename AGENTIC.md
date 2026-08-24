# AGENTIC.md

## Purpose

This file defines the working rules for AI coding agents and automated assistants modifying the **Maize Transcriptome Diversity** repository.

The project analyzes transcriptomic and genetic diversity in maize, with particular emphasis on cis-eQTL discovery, fine-mapping, regulatory annotation, population genetics, and manuscript-quality figure generation. Changes should prioritize scientific reproducibility, transparent statistical assumptions, and preservation of established analysis outputs.

## Scientific context

Core project characteristics:

- species: maize (*Zea mays*);
- primary reference: B73 RefGen_v5 unless a script explicitly states otherwise;
- main panel: WiDiv diversity panel, approximately 631 inbred lines in the principal analysis;
- transcriptomic data: leaf RNA-seq;
- genotype scale: millions of variants;
- primary eQTL model: OSCA mixed linear model with a genomic relationship matrix;
- cis window: typically ±2 Mb around the gene unless explicitly overridden;
- multiple testing: IHW is used in the principal eQTL workflow, with distance to TSS as the covariate;
- fine-mapping: SuSiE / SuSiE-RSS;
- effect-size interpretation: aFC-n or related allelic fold-change summaries where available;
- downstream analyses include GERP, PlantCAD/functional scores, genomic annotations, ATAC-seq, TF ChIP, FST, nucleotide diversity, DEG/eQTL integration, and downsampling sensitivity analyses.

Do not silently change these assumptions. If a requested modification changes the biological model, statistical threshold, reference genome, sample definition, or variant filtering, make the change explicit in code comments and documentation.

## Repository layout

```text
scripts/
├── hpc/                         # production and SLURM/HPC workflows
├── local/
│   ├── processing/              # preparation of genotype/expression/covariate inputs
│   ├── analysis/
│   │   ├── expression/          # expression diversity, CCA, heritability
│   │   ├── eqtl/                # OSCA mapping and significance/IHW
│   │   ├── fine_mapping/        # SuSiE, aFC, coloc
│   │   ├── differential_expression/ # limma DE workflows
│   │   ├── functional_annotation/   # GERP/CAD, genomic context, ATAC/TF enrichment
│   │   ├── population_genetics/ # pi, PCA, FST, evolutionary integration
│   │   ├── downsampling/        # sample-size sensitivity analyses
│   │   └── todo/                # historical or incomplete analyses
│   └── figures/                 # manuscript and supplementary figures
└── utils/                       # reusable helper functions
```

Other important locations:

- `data/raw/`: raw or externally supplied inputs; large files are not tracked;
- `data/processed/`: derived inputs used by analyses; large files are not tracked;
- `notes/`: lightweight project notes and R Markdown summaries;
- `results/`: generated outputs in the working copy; generally not tracked;
- `scripts/utils/functions.R`: shared R helper functions;
- `scripts/local/analysis/README.md`: map of analysis domains;
- `sync.sh`: conservative helper for syncing code to GitHub.

## General rules

1. Preserve scientific reproducibility over cosmetic refactoring.
2. Never commit large genomic or sequencing datasets.
3. Do not delete or rename an analysis script unless its role and dependencies are understood.
4. Prefer small, reviewable changes over broad rewrites.
5. Keep existing output file names stable unless the task explicitly requires changing them.
6. Avoid changing statistical thresholds, covariates, sample filters, or model formulas without an explicit reason.
7. When fixing a bug, document whether previous outputs may need to be regenerated.
8. Do not fabricate missing metadata, paths, sample names, or biological annotations.
9. Assume scripts are run from the repository root unless documented otherwise.
10. Prefer relative repository paths over user-specific absolute paths.

## R conventions

- Prefer clear `tidyverse` or `data.table` code consistent with the surrounding script; do not convert an entire script between paradigms without need.
- Use `data.table::fread()`/`fwrite()` for very large tabular files where practical.
- Avoid unnecessary copies of large objects.
- Set explicit factor levels when figure ordering matters.
- Use deterministic seeds for resampling, matching, permutation, or downsampling analyses.
- When parallelizing, expose the number of workers/cores as a variable rather than hard-coding an excessive value.
- Check identifiers before joins and report unexpectedly duplicated keys when they could alter biological interpretation.
- For genomic coordinates, explicitly track whether intervals are 0-based or 1-based and whether endpoints are inclusive.
- Do not silently coerce chromosome naming conventions (`1` vs `chr1`).

## Shell and HPC conventions

- Use `#!/usr/bin/env bash` for new Bash scripts unless the HPC environment requires a specific shell.
- Prefer `set -Eeuo pipefail` for new robust shell workflows when compatible with the commands used.
- Quote paths and variables.
- For SLURM jobs, state CPUs, memory, wall time, and array dimensions explicitly.
- Do not increase CPU or memory requests without a computational reason.
- Do not run destructive commands on shared HPC paths without explicit safeguards.
- Keep temporary and intermediate outputs separate from final summaries.
- Check return codes for external tools where failure could otherwise generate apparently valid downstream files.

## eQTL-specific safeguards

Before modifying eQTL workflows, verify:

- the sample order matches between expression, genotype, GRM, and covariate files;
- gene/probe identifiers are consistent with the reference annotation;
- the intended cis window is preserved;
- MAF filtering is explicit;
- the genomic relationship matrix corresponds to the same sample set;
- covariate files contain finite, non-constant columns;
- multiple-testing correction is applied to the intended hypothesis family;
- primary and sensitivity analyses are not accidentally mixed.

The main analysis and sensitivity analyses may intentionally use different significance procedures. Do not standardize them automatically.

## Fine-mapping safeguards

For SuSiE / SuSiE-RSS analyses:

- preserve the mapping between variants, alleles, genomic positions, and LD matrices;
- verify allele orientation before combining effect sizes from different tools;
- do not interpret credible sets as independent causal variants without supporting evidence;
- keep posterior inclusion probabilities and credible-set membership distinct from ordinary eQTL significance;
- record filtering that removes variants before fine-mapping.

## Functional enrichment safeguards

When comparing fine-mapped variants with annotations such as promoters, UTRs, CDS, introns, ATAC peaks, conserved elements, TF ChIP peaks, GERP, or PlantCAD:

- distinguish raw overlap from enrichment;
- use matched or otherwise appropriate background variants when the analysis requires control for MAF, TSS distance, LD, or local variant density;
- avoid treating multiple variants from the same credible set as fully independent observations without justification;
- state whether analyses use lead variants, all credible-set variants, PIP-weighted variants, or unique variants;
- preserve the distinction between primary and secondary eQTL signals where relevant.

## Population-genetic safeguards

For FST, nucleotide diversity, or selection analyses:

- confirm population labels and sample membership before running contrasts;
- preserve the intended breeding-era or subpopulation definition;
- do not compare statistics produced with incompatible filters or callable regions without explicitly noting the difference;
- when integrating population-genetic statistics with eQTL/DEG results, keep gene-level and variant-level tests distinct.

## Figure scripts

Manuscript figure scripts should:

- read already generated summary/intermediate files rather than rerunning expensive upstream analyses;
- save figures to predictable paths under `results/figures/` or the established project location;
- keep dimensions and DPI explicit;
- use consistent subpopulation ordering and naming;
- avoid manual changes that cannot be regenerated from code;
- retain panel labels and manuscript figure numbering unless intentionally revised.

## Naming and file placement

For new files:

- processing/data-preparation code: `scripts/local/processing/`;
- expression diversity/heritability: `scripts/local/analysis/expression/`;
- eQTL mapping/significance: `scripts/local/analysis/eqtl/`;
- fine-mapping/aFC/coloc: `scripts/local/analysis/fine_mapping/`;
- differential expression: `scripts/local/analysis/differential_expression/`;
- regulatory/functional interpretation: `scripts/local/analysis/functional_annotation/`;
- population genetics: `scripts/local/analysis/population_genetics/`;
- downsampling/sensitivity: `scripts/local/analysis/downsampling/`;
- manuscript plotting: `scripts/local/figures/`;
- reusable helpers: `scripts/utils/`;
- production HPC workflows: `scripts/hpc/`.

Use descriptive names. Existing numeric prefixes may reflect historical workflow order, but new files should only receive a numeric prefix if their execution order is genuinely important. Do not infer global execution order from prefixes across different analysis subdirectories.

## Data and Git policy

Do not commit:

- FASTQ/FQ;
- BAM/SAM/CRAM and indexes;
- VCF/BCF and indexes;
- large BED/BigWig files;
- RDS/RData analysis objects;
- full OSCA output collections;
- full genotype matrices;
- large generated result directories;
- credentials, tokens, `.Renviron`, or local configuration containing secrets.

Small metadata files, lightweight summaries, scripts, configuration templates, and documentation are appropriate for Git when not sensitive.

## Validation before committing

For every code change, perform the relevant subset of the following checks:

1. inspect the script for hard-coded absolute paths;
2. verify input and output paths are consistent with repository structure;
3. run syntax/parsing checks when possible;
4. check that renamed files/directories are not referenced elsewhere;
5. ensure no large generated files were added accidentally;
6. verify sample/variant/gene identifiers are not unintentionally reordered;
7. confirm that statistical parameters remain unchanged unless intentionally modified;
8. summarize any outputs that must be regenerated because of the change.

## Preferred change workflow

1. Inspect the relevant script and immediate upstream/downstream dependencies.
2. Make the smallest coherent change.
3. Update comments or documentation where scientific behavior changes.
4. Validate syntax and paths.
5. Review the diff for accidental data/output additions.
6. Commit on a dedicated branch and use a pull request for structural changes.

## Current cleanup conventions

Canonical directory names are:

- `scripts/local/processing/` — not `procesing`;
- `scripts/utils/` — not `utilis`;
- `scripts/local/analysis/todo/` — not `toDo`.

Future analysis code should use the domain-specific subdirectories described above rather than accumulating scripts directly in `scripts/local/analysis/`.
