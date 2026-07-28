#!/usr/bin/env bash
# Run SuSiE in compressed batch outputs.
# Usage: bash Run_SusieR_All_Genes.sh [parallel_jobs=5] [genes_per_batch=1500] [extra R options]
# Example: bash scripts/local/analysis/Run_SusieR_All_Genes.sh 5 1500 --window-bp 500000

set -uo pipefail
jobs="${1:-5}"
batch_size="${2:-1500}"
if ! [[ "$jobs" =~ ^[1-9][0-9]*$ && "$batch_size" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Parallel jobs and genes per batch must be positive integers.\n' >&2; exit 2
fi
shift 2 || true
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
worker="$script_dir/5.SusieR.R"
project_root="$(cd -- "$script_dir/../../.." && pwd)"
output_dir="$project_root/results/full/SusieR"
batch_dir="$output_dir/batches"
log_dir="$output_dir/logs"
mkdir -p "$batch_dir" "$log_dir"

mapfile -t genes < <(Rscript "$worker" --list-genes "$@")
if ((${#genes[@]} == 0)); then printf 'No eligible genes found.\n' >&2; exit 1; fi
printf '%s\n' "${genes[@]}" | split -l "$batch_size" -d -a 4 --additional-suffix=.txt - "$batch_dir/genes_"

running=0; failed=0; pids=()
wait_for_one() {
  local pid status; pid="${pids[0]}"
  if wait "$pid"; then status=0; else status=$?; failed=1; fi
  pids=("${pids[@]:1}"); running=$((running - 1)); return "$status"
}
for batch_file in "$batch_dir"/genes_*.txt; do
  batch_id="${batch_file##*/genes_}"; batch_id="${batch_id%.txt}"
  pip_file="$output_dir/PIP_batch_${batch_id}.tsv.gz"
  cs_file="$output_dir/CS_batch_${batch_id}.tsv.gz"
  if [[ -s "$pip_file" && -s "$cs_file" ]]; then printf 'Skipping batch %s (outputs exist)\n' "$batch_id"; continue; fi
  while ((running >= jobs)); do wait_for_one || true; done
  printf 'Starting batch %s (%s genes)\n' "$batch_id" "$(wc -l < "$batch_file")"
  Rscript "$worker" --genes-file "$batch_file" --batch-id "$batch_id" "$@" >"$log_dir/batch_${batch_id}.log" 2>&1 &
  pids+=("$!"); running=$((running + 1))
done
while ((running > 0)); do wait_for_one || true; done
if ((failed)); then printf 'One or more batches had failed genes; inspect %s and errors_batch_*.tsv.gz.\n' "$log_dir" >&2; exit 1; fi
printf 'Completed all batches.\n'
