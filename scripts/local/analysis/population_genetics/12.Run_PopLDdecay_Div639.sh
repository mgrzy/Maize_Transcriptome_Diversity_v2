#!/usr/bin/env bash
# Calculate genome-wide LD decay in the Div639 VCF with PopLDdecay.
# The default common-variant filter and 250-kb maximum distance keep the
# pairwise calculation tractable for this dense sequencing dataset.
#
# Usage (from the repository root):
#   bash scripts/local/analysis/population_genetics/12.Run_PopLDdecay_Div639.sh
#   bash scripts/local/analysis/population_genetics/12.Run_PopLDdecay_Div639.sh results/full/LD/PopLDdecay_Div639
#
# Optionally override the executable, e.g.:
#   POPLDDECAY_BIN=/path/to/PopLDdecay bash scripts/local/analysis/population_genetics/12.Run_PopLDdecay_Div639.sh

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
vcf_file="${repo_root}/data/processed/geno/vcf/Div639.vcf.gz"
poplddecay_bin="${POPLDDECAY_BIN:-${repo_root}/../../soft/PopLDdecay/bin/PopLDdecay}"
out_prefix="${1:-results/full/LD/PopLDdecay_Div639_MAF0.05_MaxDist250}"

[[ -f "${vcf_file}" ]] || { echo "Input VCF not found: ${vcf_file}" >&2; exit 1; }
[[ -x "${poplddecay_bin}" ]] || { echo "PopLDdecay executable not found or not executable: ${poplddecay_bin}" >&2; exit 1; }

mkdir -p "$(dirname "${out_prefix}")"

"${poplddecay_bin}" \
  -InVCF "${vcf_file}" \
  -OutStat "${out_prefix}" \
  -MAF 0.05 \
  -MaxDist 250
