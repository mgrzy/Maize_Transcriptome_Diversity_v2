#!/usr/bin/env bash
### Script: scripts/local/analysis/11.Run_Downsampled_OSCA.sh
### Purpose: Run one complete cis-eQTL OSCA GLM analysis for each genotype
###          downsample, without splitting analyses into task chunks.
### Inputs (per downsample):
###  - data/processed/geno/plink/downsample/DownSampling_<n> (PLINK files)
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n> (BOD files)
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>_PC.txt
###    (five genetic PCs and five expression PCs)
### Outputs:
###  - results/full/osca/downsample/DownSampling_<n>/cis2Mb.*
### Usage: bash scripts/local/analysis/11.Run_Downsampled_OSCA.sh
### Prerequisites: conda environment 'mkl'; OSCA at ../../soft/osca/osca

set -euo pipefail

# Work from the repository root so project-relative paths remain valid.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/../../.."

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mkl
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

OSCA="${OSCA:-../../soft/osca/osca}"
THREADS="${OSCA_THREADS:-48}"
TARGETS=(100 150 200 300 400 500 550 600)
GENO_ROOT="data/processed/geno/plink/downsample"
EXPR_ROOT="data/processed/osca/downsample"
RESULTS_ROOT="results/full/osca/downsample"

if [[ ! -x "${OSCA}" ]]; then
  echo "OSCA executable not found or not executable: ${OSCA}" >&2
  exit 1
fi

# Validate every dependency before starting, preventing a partially completed
# downsampling series when an upstream PCA or OSCA-preparation step is missing.
for n in "${TARGETS[@]}"; do
  prefix="DownSampling_${n}"
  geno_prefix="${GENO_ROOT}/${prefix}"
  expr_prefix="${EXPR_ROOT}/${prefix}/${prefix}"

  for file in \
    "${geno_prefix}.bed" \
    "${geno_prefix}.bim" \
    "${geno_prefix}.fam" \
    "${expr_prefix}.bod" \
    "${expr_prefix}.oii" \
    "${expr_prefix}.opi" \
    "${expr_prefix}_PC.txt"; do
    if [[ ! -f "${file}" ]]; then
      echo "Missing OSCA input: ${file}" >&2
      exit 1
    fi
  done
done

for n in "${TARGETS[@]}"; do
  prefix="DownSampling_${n}"
  geno_prefix="${GENO_ROOT}/${prefix}"
  expr_prefix="${EXPR_ROOT}/${prefix}/${prefix}"
  output_dir="${RESULTS_ROOT}/${prefix}"

  mkdir -p "${output_dir}"

  "${OSCA}" \
    --eqtl \
    --bfile "${geno_prefix}" \
    --befile "${expr_prefix}" \
    --qcovar "${expr_prefix}_PC.txt" \
    --thread-num "${THREADS}" \
    --cis \
    --cis-wind 2000 \
    --no-fid \
    --out "${output_dir}/cis2Mb"

  echo "Completed complete OSCA analysis for ${prefix}."
done
