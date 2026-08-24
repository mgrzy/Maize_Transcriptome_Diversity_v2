#!/usr/bin/env bash
# Run narrow-sense heritability in independent, serial R jobs and merge them.
# This avoids forking sommer::mmer(), which can stall when multiple fits share
# a forked R parent.
#
# Usage:
#   bash scripts/local/analysis/3.Run_Narrow_sense_h2_Shards.sh [JOB_COUNT]
#
# Examples:
#   # Analyze every eligible gene using 10 independent R processes.
#   bash scripts/local/analysis/3.Run_Narrow_sense_h2_Shards.sh 10
#
#   # Smoke-test the first 10 eligible genes, split across 10 jobs.
#   MAX_GENES=10 bash scripts/local/analysis/3.Run_Narrow_sense_h2_Shards.sh 10
#
# Optional environment variables:
#   JOB_COUNT  Number of independent jobs (default: 10; positional arg wins).
#   MAX_GENES  Global cap on alphabetically sorted eligible genes (default: all).
#   R_SCRIPT   Path to 3.Narrow_sense_h2.R.
#   LOG_DIR    Directory for one log per job.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
R_SCRIPT="${R_SCRIPT:-${SCRIPT_DIR}/3.Narrow_sense_h2.R}"
JOB_COUNT="${1:-${JOB_COUNT:-10}}"
MAX_GENES="${MAX_GENES:-}"
LOG_DIR="${LOG_DIR:-${PROJECT_ROOT}/results/logs/narrow_sense_h2_jobs}"

if [[ "${JOB_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
  :
else
  echo "JOB_COUNT must be a positive integer; received: ${JOB_COUNT}" >&2
  exit 1
fi

if [[ -n "${MAX_GENES}" && ! "${MAX_GENES}" =~ ^[1-9][0-9]*$ ]]; then
  echo "MAX_GENES must be a positive integer when supplied; received: ${MAX_GENES}" >&2
  exit 1
fi

if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "R script not found: ${R_SCRIPT}" >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript was not found in PATH." >&2
  exit 1
fi

mkdir -p "${LOG_DIR}"

declare -a job_ids=()
declare -a job_pids=()

echo "Launching ${JOB_COUNT} independent narrow-sense h2 job(s)."
echo "Logs: ${LOG_DIR}"

for job_index in $(seq 1 "${JOB_COUNT}"); do
  job_log="${LOG_DIR}/job_$(printf '%03d' "${job_index}")_of_$(printf '%03d' "${JOB_COUNT}").log"
  r_args=(
    "${R_SCRIPT}"
    --job-index "${job_index}"
    --job-count "${JOB_COUNT}"
    --cores 1
  )
  if [[ -n "${MAX_GENES}" ]]; then
    r_args+=(--max-genes "${MAX_GENES}")
  fi

  echo "Starting job ${job_index}/${JOB_COUNT}: ${job_log}"
  Rscript "${r_args[@]}" >"${job_log}" 2>&1 &
  job_ids+=("${job_index}")
  job_pids+=("$!")
done

failed_jobs=()
for array_index in "${!job_pids[@]}"; do
  if wait "${job_pids[${array_index}]}"; then
    echo "Finished job ${job_ids[${array_index}]}/${JOB_COUNT}."
  else
    echo "FAILED job ${job_ids[${array_index}]}/${JOB_COUNT}; see its log." >&2
    failed_jobs+=("${job_ids[${array_index}]}")
  fi
done

if (( ${#failed_jobs[@]} )); then
  echo "Not merging because job(s) failed: ${failed_jobs[*]}" >&2
  exit 1
fi

echo "All jobs completed; merging shards."
Rscript "${R_SCRIPT}" --merge-only --job-count "${JOB_COUNT}"

