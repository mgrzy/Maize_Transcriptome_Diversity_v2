#!/usr/bin/env bash
### Script: scripts/local/analysis/12.Extract_Downsampled_OSCA_Significant.sh
### Purpose: Export significant cis-eQTL associations from each downsampled
###          OSCA BESD result using a p-value threshold of 3.221153e-09.
### Inputs (per downsample):
###  - results/full/osca/downsample/DownSampling_<n>/cis2Mb.{besd,esi,epi}
### Outputs (per downsample):
###  - results/full/osca/downsample/DownSampling_<n>/significant_p3.221153e-09
### Usage: bash scripts/local/analysis/12.Extract_Downsampled_OSCA_Significant.sh
### Prerequisites: OSCA at ../../soft/osca/osca

set -euo pipefail

# Work from the repository root so project-relative paths remain valid.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/../../.."

OSCA="${OSCA:-../../soft/osca/osca}"
P_THRESHOLD="${P_THRESHOLD:-3.221153e-09}"
TARGETS=(100 150 200 300 400 500 550 600)
RESULTS_ROOT="results/full/osca/downsample"

if [[ ! -x "${OSCA}" ]]; then
  echo "OSCA executable not found or not executable: ${OSCA}" >&2
  exit 1
fi

# Require the full series before exporting, avoiding a partial set of summary
# files when one of the downsampled OSCA analyses is still incomplete.
for n in "${TARGETS[@]}"; do
  result_prefix="${RESULTS_ROOT}/DownSampling_${n}/cis2Mb"
  for file in "${result_prefix}.besd" "${result_prefix}.esi" "${result_prefix}.epi"; do
    if [[ ! -f "${file}" ]]; then
      echo "Missing OSCA BESD result: ${file}" >&2
      exit 1
    fi
  done
done

for n in "${TARGETS[@]}"; do
  result_prefix="${RESULTS_ROOT}/DownSampling_${n}/cis2Mb"
  output_prefix="${RESULTS_ROOT}/DownSampling_${n}/significant_p${P_THRESHOLD}"

  "${OSCA}" \
    --beqtl-summary "${result_prefix}" \
    --query "${P_THRESHOLD}" \
    --out "${output_prefix}"

  echo "Extracted associations for DownSampling_${n} with p <= ${P_THRESHOLD}."
done
