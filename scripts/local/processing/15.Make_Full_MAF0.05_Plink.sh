#!/usr/bin/env bash
### Script: scripts/local/procesing/15.Make_Full_MAF0.05_Plink.sh
### Purpose: Create the full-cohort PLINK genotype dataset filtered at MAF 0.05.
### Input:
###  - data/processed/geno/plink/Div639 (PLINK bed/bim/fam)
### Output:
###  - data/processed/geno/plink/Div639_maf0.05 (PLINK bed/bim/fam)
### Usage: bash scripts/local/procesing/15.Make_Full_MAF0.05_Plink.sh
### Prerequisites: conda environment 'geno' with PLINK

set -euo pipefail

# Work from the repository root so project-relative paths remain valid.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/../../.."

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geno

PLINK="${PLINK:-plink}"
INPUT_PREFIX="data/processed/geno/plink/Div639"
OUTPUT_PREFIX="data/processed/geno/plink/Div639_maf0.05"

for file in "${INPUT_PREFIX}.bed" "${INPUT_PREFIX}.bim" "${INPUT_PREFIX}.fam"; do
  if [[ ! -f "${file}" ]]; then
    echo "Missing PLINK input: ${file}" >&2
    exit 1
  fi
done

"${PLINK}" \
  --bfile "${INPUT_PREFIX}" \
  --maf 0.05 \
  --make-bed \
  --out "${OUTPUT_PREFIX}"
