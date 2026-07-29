#!/usr/bin/env bash
# Run aFCn for Div639 lead eQTLs with the phased, sample-subsetted VCF.
#
# Inputs:
#   data/processed/geno/vcf/Div639.phased.vcf.gz
#   data/processed/AF/Matrix_AF.csv.gz
#   data/processed/AF/eQTL_List.tsv
# Output:
#   results/full/AF/aFCn.txt
#
# Usage:
#   bash scripts/local/analysis/8.aFC.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "Conda was not found in PATH." >&2
  exit 1
fi

# aFCn's Python dependencies are installed in the base Conda environment.
# Conda's Qt activation hook references an unset variable, so nounset must be
# disabled only while the environment is activated.
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${AFCN_CONDA_ENV:-base}"
set -u

AFCN="${AFCN:-${PROJECT_ROOT}/../../soft/aFCn/src/afcn.py}"
PYTHON="${PYTHON:-python}"
VCF="${PHASED_VCF:-${PROJECT_ROOT}/data/processed/geno/vcf/Div639.phased.vcf.gz}"
EXPR="${PROJECT_ROOT}/data/processed/AF/Matrix_AF.csv.gz"
EQTL="${PROJECT_ROOT}/data/processed/AF/eQTL_List.tsv"
OUTPUT_DIR="${PROJECT_ROOT}/results/full/AF"

for input_file in "${AFCN}" "${VCF}" "${EXPR}" "${EQTL}"; do
  if [[ ! -f "${input_file}" ]]; then
    echo "Missing input file: ${input_file}" >&2
    exit 1
  fi
done

if ! command -v "${PYTHON}" >/dev/null 2>&1; then
  echo "Python executable not found: ${PYTHON}" >&2
  exit 1
fi

if ! "${PYTHON}" -c 'import pandas, pysam, lmfit, sklearn' >/dev/null 2>&1; then
  echo "The selected Python lacks aFCn dependencies." >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

"${PYTHON}" "${AFCN}" \
  -j 36 \
  --vcf "${VCF}" \
  --expr "${EXPR}" \
  --eqtl "${EQTL}" \
  --output "${OUTPUT_DIR}/aFCn.txt"
