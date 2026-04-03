#!/bin/bash -l
#SBATCH --job-name=pi_bed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --account=g102-2488
#SBATCH --partition=okeanos
#SBATCH --array=1-6
#SBATCH --error=logs/pi_by_interval/%x.%A_%a.err
#SBATCH --output=logs/pi_by_interval/%x.%A_%a.out

module load common/anaconda/3.8
conda activate renv

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)

find_project_root() {
  local base level candidate

  for base in "${SLURM_SUBMIT_DIR:-}" "${PWD}" "${SCRIPT_DIR}"; do
    [[ -z "${base}" ]] && continue
    base=$(cd "${base}" 2>/dev/null && pwd) || continue

    for level in . .. ../.. ../../.. ../../../..; do
      candidate=$(cd "${base}/${level}" 2>/dev/null && pwd) || continue
      if [[ -d "${candidate}/data" && -d "${candidate}/scripts" ]]; then
        printf '%s\n' "${candidate}"
        return 0
      fi
    done
  done

  echo "Unable to determine project root from SLURM_SUBMIT_DIR, PWD, or script location" >&2
  exit 1
}

PROJECT_ROOT=$(find_project_root)
GROUP_DIR="${PROJECT_ROOT}/data/processed/meta/temperate_bam_by_group"
RSCRIPT_PATH="${PROJECT_ROOT}/scripts/hpc/pop/Calculate_Pi_Per_Bed_Interval.R"
LOG_DIR="${PROJECT_ROOT}/logs/pi_by_interval"

mkdir -p "${LOG_DIR}"
cd "${PROJECT_ROOT}"

if [[ ! -d "${GROUP_DIR}" ]]; then
  echo "Group BAM directory not found: ${GROUP_DIR}" >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is not available in the renv conda environment" >&2
  exit 1
fi

mapfile -t group_files < <(find "${GROUP_DIR}" -maxdepth 1 -type f -name '*_BAM.txt' | sort)

if [[ ${#group_files[@]} -ne 6 ]]; then
  echo "Expected 6 group BAM files in ${GROUP_DIR}, found ${#group_files[@]}" >&2
  exit 1
fi

task_id=${SLURM_ARRAY_TASK_ID:-1}
if (( task_id < 1 || task_id > ${#group_files[@]} )); then
  echo "Array task ${task_id} is out of range for ${#group_files[@]} groups" >&2
  exit 1
fi

group_file="${group_files[$((task_id - 1))]}"
group_name=$(basename "${group_file}")
group_name="${group_name%_BAM.txt}"

echo "Processing PI-by-interval for group: ${group_name}"

Rscript "${RSCRIPT_PATH}" --group="${group_name}"

conda deactivate
