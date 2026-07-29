#!/usr/bin/env bash
# Convert the filtered Div639 PLINK dataset to a bgzip-compressed VCF.
#
# Input:  data/processed/geno/plink/Div639.{bed,bim,fam}
# Output: data/processed/geno/vcf/Div639.vcf.gz
#
# Usage from any directory:
#   bash scripts/local/procesing/Convert_Div639_plink_to_vcf.sh
#
# Prerequisite: a Conda environment named `geno` containing PLINK.
#
# Set PLINK to use a non-default executable, e.g.:
#   PLINK=/path/to/plink bash scripts/local/procesing/Convert_Div639_plink_to_vcf.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

PLINK="${PLINK:-plink}"
INPUT_PREFIX="${PROJECT_ROOT}/data/processed/geno/plink/Div639"
OUTPUT_DIR="${PROJECT_ROOT}/data/processed/geno/vcf"
OUTPUT_PREFIX="${OUTPUT_DIR}/Div639"

if ! command -v conda >/dev/null 2>&1; then
  echo "Conda was not found in PATH." >&2
  exit 1
fi

# Make `conda activate` available in this non-interactive shell, then load PLINK.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geno

for extension in bed bim fam; do
  if [[ ! -f "${INPUT_PREFIX}.${extension}" ]]; then
    echo "Missing input file: ${INPUT_PREFIX}.${extension}" >&2
    exit 1
  fi
done

if ! command -v "${PLINK}" >/dev/null 2>&1; then
  echo "PLINK executable not found: ${PLINK}" >&2
  echo "Activate the genotype environment or set PLINK=/path/to/plink." >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

"${PLINK}" \
  --bfile "${INPUT_PREFIX}" \
  --recode vcf-iid bgz \
  --out "${OUTPUT_PREFIX}"

echo "Created ${OUTPUT_PREFIX}.vcf.gz"
