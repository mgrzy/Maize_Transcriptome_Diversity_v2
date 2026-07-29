#!/usr/bin/env bash
# Subset the phased All_maf001 VCF to the Div639 variants and lines.
#
# Input:
#   ../BigData/Maize/v5/Genotypes/All/imputed/VCF/All_maf001.vcf.gz
#   data/processed/geno/plink/Div639.{bim,fam}
# Output:
#   data/processed/geno/vcf/Div639.phased.vcf.gz (+ .tbi index)
#
# Usage from the project root:
#   bash scripts/local/procesing/Subset_Div639_phased_VCF.sh
#
# Prerequisite: Conda environment `geno` with bcftools installed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

SOURCE_VCF="${SOURCE_VCF:-${PROJECT_ROOT}/../BigData/Maize/v5/Genotypes/All/imputed/VCF/All_maf001.vcf.gz}"
PLINK_PREFIX="${PLINK_PREFIX:-${PROJECT_ROOT}/data/processed/geno/plink/Div639}"
OUTPUT_VCF="${OUTPUT_VCF:-${PROJECT_ROOT}/data/processed/geno/vcf/Div639.phased.vcf.gz}"

if ! command -v conda >/dev/null 2>&1; then
  echo "Conda was not found in PATH." >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geno

BCFTOOLS="${BCFTOOLS:-bcftools}"

for input_file in "${SOURCE_VCF}" "${PLINK_PREFIX}.bim" "${PLINK_PREFIX}.fam"; do
  if [[ ! -f "${input_file}" ]]; then
    echo "Missing input file: ${input_file}" >&2
    exit 1
  fi
done

if ! command -v "${BCFTOOLS}" >/dev/null 2>&1; then
  echo "bcftools executable not found: ${BCFTOOLS}" >&2
  exit 1
fi

if [[ -e "${OUTPUT_VCF}" || -e "${OUTPUT_VCF}.tbi" ]]; then
  echo "Output already exists: ${OUTPUT_VCF}" >&2
  echo "Choose another OUTPUT_VCF path or remove the existing output first." >&2
  exit 1
fi

OUTPUT_DIR="$(dirname "${OUTPUT_VCF}")"
mkdir -p "${OUTPUT_DIR}"

TEMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TEMP_DIR}"' EXIT

SAMPLE_IDS="${TEMP_DIR}/Div639.samples.txt"
SNP_IDS="${TEMP_DIR}/Div639.snps.txt"

# .fam column 2 is the IID used in the imputed VCF; .bim column 2 is the
# chr<chromosome>_<position> variant ID used by All_maf001.
awk '{print $2}' "${PLINK_PREFIX}.fam" > "${SAMPLE_IDS}"
awk '{print $2}' "${PLINK_PREFIX}.bim" > "${SNP_IDS}"

"${BCFTOOLS}" view \
  --samples-file "${SAMPLE_IDS}" \
  --include "ID=@${SNP_IDS}" \
  --output-type z \
  --output-file "${OUTPUT_VCF}" \
  "${SOURCE_VCF}"

"${BCFTOOLS}" index --tbi --force "${OUTPUT_VCF}"

echo "Created ${OUTPUT_VCF}"
