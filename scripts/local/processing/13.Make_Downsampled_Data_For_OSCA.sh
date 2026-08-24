#!/usr/bin/env bash
### Script: scripts/local/procesing/13.Make_Downsampled_Data_For_OSCA.sh
### Purpose: Convert each prepared downsampled expression matrix to OSCA BOD
###          format and attach its matching probe information (OPI).
### Inputs (per downsample):
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>.dat
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>_OPI.txt
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>_OII.txt
###  - data/processed/osca/downsample/DownSampling_<n>/DownSampling_<n>_PC.txt
### Outputs (per downsample):
###  - DownSampling_<n>.bod, .oii, .opi, and related OSCA files
### Usage: bash scripts/local/procesing/13.Make_Downsampled_Data_For_OSCA.sh
### Prerequisites: OSCA at ../../soft/osca/osca

set -euo pipefail

# Work from the repository root so project-relative paths remain valid.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/../../.."

OSCA="${OSCA:-../../soft/osca/osca}"
INPUT_ROOT="data/processed/osca/downsample"
TARGETS=(100 150 200 300 400 500 550 600)

if [[ ! -x "${OSCA}" ]]; then
  echo "OSCA executable not found or not executable: ${OSCA}" >&2
  exit 1
fi

# Check every required input before writing anything, avoiding a partially
# converted set if the preceding R preparation has not completed.
for n in "${TARGETS[@]}"; do
  prefix="DownSampling_${n}"
  input_dir="${INPUT_ROOT}/${prefix}"

  for file in \
    "${input_dir}/${prefix}.dat" \
    "${input_dir}/${prefix}_OPI.txt" \
    "${input_dir}/${prefix}_OII.txt" \
    "${input_dir}/${prefix}_PC.txt"; do
    if [[ ! -f "${file}" ]]; then
      echo "Missing prepared OSCA input: ${file}" >&2
      exit 1
    fi
  done
done

for n in "${TARGETS[@]}"; do
  prefix="DownSampling_${n}"
  output_dir="${INPUT_ROOT}/${prefix}"
  output_prefix="${output_dir}/${prefix}"

  "${OSCA}" \
    --efile "${output_prefix}.dat" \
    --gene-expression \
    --make-bod \
    --no-fid \
    --out "${output_prefix}"

  "${OSCA}" \
    --befile "${output_prefix}" \
    --update-opi "${output_prefix}_OPI.txt"

  echo "Completed OSCA BOD preparation for ${prefix}."
done
