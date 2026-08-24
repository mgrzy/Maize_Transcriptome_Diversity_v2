# Analysis directory map

The local analysis scripts are grouped by biological/statistical role. Historical numeric prefixes are retained because they may encode execution history; the folder name is the primary organizational unit.

## `expression/`

Expression-level summaries and variance components:

- transcriptome diversity;
- CCA;
- narrow-sense expression heritability.

## `eqtl/`

Primary cis-eQTL mapping and significance processing:

- covariate diagnostics;
- OSCA mixed-model runs;
- adjusted/significant association selection;
- IHW multiple-testing workflow;
- full-analysis reruns and extraction.

## `fine_mapping/`

Post-eQTL causal-resolution and effect-size analyses:

- SuSiE fine-mapping;
- batch SuSiE execution;
- aFC calculations;
- coloc analyses.

## `differential_expression/`

Differential-expression analysis and summaries using limma-based workflows.

## `functional_annotation/`

Functional interpretation of eQTL/fine-mapped variants:

- genomic context and distance enrichment;
- GERP and CAD/zero-shot scores;
- GO enrichment;
- chromatin and ATAC overlap;
- TF ChIP and TF-class enrichment.

## `population_genetics/`

Population structure and evolutionary analyses:

- genotype PCA;
- nucleotide diversity;
- FST;
- integration of DEG/eQTL effects with population-genetic statistics.

## `downsampling/`

Sample-size sensitivity analysis for eQTL discovery.

## `todo/`

Historical, incomplete, or not-yet-integrated scripts. Do not treat these as part of the validated primary workflow without inspection.

## Execution notes

Scripts are generally expected to be run from the repository root, so moving them between these folders does not change project-relative `data/...` or `results/...` paths.

Do not infer execution order solely from numeric prefixes across folders. Within a domain, inspect input/output paths and comments before rerunning a workflow.
