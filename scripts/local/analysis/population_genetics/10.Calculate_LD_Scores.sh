#!/usr/bin/env bash
# Calculate LD scores for the WiDiv (Div639) PLINK dataset with GCTA.
#
# Usage (from the repository root):
#   bash scripts/local/analysis/population_genetics/10.Calculate_LD_Scores.sh
#   bash scripts/local/analysis/population_genetics/10.Calculate_LD_Scores.sh results/full/LD/LDWiDiv
#
# Optionally override the GCTA executable, e.g.:
#   GCTA_BIN=/path/to/gcta64 bash scripts/local/analysis/population_genetics/10.Calculate_LD_Scores.sh

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
bfile="${repo_root}/data/processed/geno/plink/Div639"
gcta_bin="${GCTA_BIN:-${repo_root}/../../soft/gcta/gcta}"
out_prefix="${1:-LDWiDiv}"

required_files=("${bfile}.bed" "${bfile}.bim" "${bfile}.fam")
for file in "${required_files[@]}"; do
  [[ -f "${file}" ]] || { echo "Missing required PLINK file: ${file}" >&2; exit 1; }
done

[[ -x "${gcta_bin}" ]] || { echo "GCTA executable not found or not executable: ${gcta_bin}" >&2; exit 1; }

mkdir -p "$(dirname "${out_prefix}")"

"${gcta_bin}" \
  --bfile "${bfile}" \
  --ld-score \
  --ld-wind 2000 \
  --ld-rsq-cutoff 0.01 \
  --thread-num 24 \
  --out "${out_prefix}"
