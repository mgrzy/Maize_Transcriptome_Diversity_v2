#!/usr/bin/env bash
### Script: scripts/local/procesing/11.DowSampling.sh
### Purpose: Filter Div639 to each downsampled line set, then calculate
###          genotype PCs and a GCTA genomic relationship matrix (GRM).
### Inputs:
###  - data/processed/geno/plink/Div639 (PLINK bed/bim/fam)
###  - data/processed/meta/downsample/DownSampling_<n>.fam
### Outputs:
###  - data/processed/geno/plink/downsample/DownSampling_<n> (PLINK bed/bim/fam)
###  - results/full/PCA/downsample/DownSampling_<n>.eigenvec
###  - results/full/GRM/downsample/DownSampling_<n>.grm.{bin,N,id}
### Usage: bash scripts/local/procesing/11.DowSampling.sh

set -euo pipefail

# Work from the repository root so the project-relative paths below are valid.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/../../.."

# Initialise Conda for this non-interactive shell before activating the PLINK
# environment.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geno

PLINK="${PLINK:-plink}"
GCTA="${GCTA:-../../soft/gcta/gcta}"
GCTA_THREADS="${GCTA_THREADS:-24}"
BASE_BFILE="data/processed/geno/plink/Div639"
KEEP_DIR="data/processed/meta/downsample"
OUT_DIR="data/processed/geno/plink/downsample"
PCA_DIR="results/full/PCA/downsample"
GRM_DIR="results/full/GRM/downsample"

if [[ ! -x "${GCTA}" ]]; then
  echo "GCTA executable not found or not executable: ${GCTA}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}" "${PCA_DIR}" "${GRM_DIR}"

for n in 100 150 200 300 400 500 550 600; do
  keep_file="${KEEP_DIR}/DownSampling_${n}.fam"
  out_prefix="${OUT_DIR}/DownSampling_${n}"

  if [[ ! -f "${keep_file}" ]]; then
    echo "Missing keep file: ${keep_file}" >&2
    exit 1
  fi

  "${PLINK}" \
    --bfile "${BASE_BFILE}" \
    --keep "${keep_file}" \
    --maf 0.05 \
    --make-bed \
    --out "${out_prefix}"

  "${PLINK}" \
    --bfile "${out_prefix}" \
    --pca 5 \
    --out "${PCA_DIR}/DownSampling_${n}"

  "${GCTA}" \
    --bfile "${out_prefix}" \
    --make-grm \
    --thread-num "${GCTA_THREADS}" \
    --out "${GRM_DIR}/DownSampling_${n}"
done
