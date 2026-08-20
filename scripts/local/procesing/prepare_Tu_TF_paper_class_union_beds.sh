#!/usr/bin/env bash
# Bash entry point for prepare_Tu_TF_paper_class_union_beds.R.
#
# Run from anywhere:
#   bash scripts/local/procesing/prepare_Tu_TF_paper_class_union_beds.sh
#
# The environment-variable configuration documented in the R script is passed
# through unchanged, e.g. OUTPUT_DIR=/path/to/output bash <this script>.

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd -- "${script_dir}/../../.." && pwd)"

exec Rscript "${project_root}/scripts/local/procesing/prepare_Tu_TF_paper_class_union_beds.R" "$@"
