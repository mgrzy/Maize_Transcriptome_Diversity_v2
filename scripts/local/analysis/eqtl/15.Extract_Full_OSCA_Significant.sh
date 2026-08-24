#!/usr/bin/env bash
### Script: scripts/local/analysis/15.Extract_Full_OSCA_Significant.sh
### Purpose: Export significant full-data OSCA cis-eQTL associations from the
###          BESD result using a p-value threshold of 3.221153e-09.
### Input:
###  - results/full/osca/full_5Geno_5ExprPC/cis2Mb.{besd,esi,epi}
### Output:
###  - results/full/osca/full_5Geno_5ExprPC/significant_p3.221153e-09
### Usage: bash scripts/local/analysis/15.Extract_Full_OSCA_Significant.sh
### Prerequisites: OSCA at ../../soft/osca/osca

set -euo pipefail

# Work from the repository root so project-relative paths remain valid.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/../../.."

OSCA="${OSCA:-../../soft/osca/osca}"
P_THRESHOLD="${P_THRESHOLD:-3.221153e-09}"
RESULT_PREFIX="results/full/osca/full_5Geno_5ExprPC/cis2Mb"
OUTPUT_PREFIX="results/full/osca/full_5Geno_5ExprPC/significant_p${P_THRESHOLD}"

if [[ ! -x "${OSCA}" ]]; then
  echo "OSCA executable not found or not executable: ${OSCA}" >&2
  exit 1
fi

for file in "${RESULT_PREFIX}.besd" "${RESULT_PREFIX}.esi" "${RESULT_PREFIX}.epi"; do
  if [[ ! -f "${file}" ]]; then
    echo "Missing OSCA BESD result: ${file}" >&2
    exit 1
  fi
done

"${OSCA}" \
  --beqtl-summary "${RESULT_PREFIX}" \
  --query "${P_THRESHOLD}" \
  --out "${OUTPUT_PREFIX}"

echo "Extracted full-data associations with p <= ${P_THRESHOLD}."
