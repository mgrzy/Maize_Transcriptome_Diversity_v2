#!/usr/bin/env bash
# Compute allelic fold changes for Div639 lead eQTLs.
#
# Inputs:
#   data/processed/geno/vcf/Div639.phased.vcf.gz
#   data/processed/AF/Matrix_AF.csv.gz
#   data/processed/AF/eQTL_List.tsv
# Output:
#   results/full/AF/aFCn.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

AFCN="${AFCN:-${PROJECT_ROOT}/../../soft/aFCn/src/afcn.py}"
PYTHON="${PYTHON:-python}"
# aFC requires phased genotypes and chromosome names matching the eQTL IDs
# (e.g. chr1_248511).  Create this subset with
# scripts/local/procesing/Subset_Div639_phased_VCF.sh.
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
  echo "The selected Python lacks aFC dependencies." >&2
  echo "Activate the Conda base environment or set PYTHON to its Python executable." >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

"${PYTHON}" "${AFCN}" \
  -j 36 \
  --vcf "${VCF}" \
  --expr "${EXPR}" \
  --eqtl "${EQTL}" \
  --output "${OUTPUT_DIR}/aFCn.txt"
