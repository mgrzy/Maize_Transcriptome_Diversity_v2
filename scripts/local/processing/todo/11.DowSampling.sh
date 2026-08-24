#!/usr/bin/env bash
### Script: processing/11.DowSampling.sh
### Purpose: Filter Div639 to each downsampled line set and apply MAF 0.05.
### Inputs:
###  - data/geno/plink/Div639 (PLINK bed/bim/fam)
###  - data/meta/downsample/DownSampling_<n>.fam
### Outputs:
###  - data/geno/plink/downsample/DownSampling_<n> (PLINK bed/bim/fam)
### Usage: bash processing/11.DowSampling.sh

set -euo pipefail

# Initialise Conda for this non-interactive shell before activating the PLINK
# environment.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geno

PLINK="${PLINK:-plink}"
BASE_BFILE="data/geno/plink/Div639"
KEEP_DIR="data/meta/downsample"
OUT_DIR="data/geno/plink/downsample"

mkdir -p "${OUT_DIR}"

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
done
