#!/bin/sh
#SBATCH --job-name=downTemp        # Descriptive job name
#SBATCH --nodes=1                      # One node
#SBATCH --ntasks-per-node=1            # One task per node
#SBATCH --cpus-per-task=1              # CPU threads for bcftools
#SBATCH --mem=3G                       # Memory
#SBATCH --time=48:00:00                # Max runtime
#SBATCH --account=g102-2488            # Allocation/account
#SBATCH --partition=okeanos            # Partition/queue
#SBATCH --error=logs/download/dowloadTemp.%A_%a.err
#SBATCH --output=logs/download/dowloadTemp.%A_%a.out
#SBATCH --array=1-10                   # chr_1..chr_10

set -euo pipefail

FASTQ_LIST="data/raw/meta/fastq_links/Temperate.txt"
OUT_DIR="data/raw/fastq/temp"
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
TASK_COUNT="${SLURM_ARRAY_TASK_COUNT:-10}"

mkdir -p "$OUT_DIR" "logs/download"

mapfile -t FASTQ_URLS < <(grep -v '^[[:space:]]*$' "$FASTQ_LIST")

TOTAL_URLS="${#FASTQ_URLS[@]}"
if (( TOTAL_URLS == 0 )); then
  echo "No FASTQ URLs found in $FASTQ_LIST"
  exit 1
fi

echo "Task ${TASK_ID}/${TASK_COUNT}"
echo "Found ${TOTAL_URLS} FASTQ URLs in $FASTQ_LIST"
echo "Saving downloads to $OUT_DIR"

DOWNLOADS_FOR_TASK=0

for INDEX in "${!FASTQ_URLS[@]}"; do
  if (( INDEX % TASK_COUNT != TASK_ID - 1 )); then
    continue
  fi

  URL="${FASTQ_URLS[$INDEX]}"
  if [[ "$URL" != http://* && "$URL" != https://* && "$URL" != ftp://* ]]; then
    URL="https://${URL}"
  fi

  FILE_NAME="$(basename "$URL")"
  DOWNLOADS_FOR_TASK=$((DOWNLOADS_FOR_TASK + 1))

  echo "[$DOWNLOADS_FOR_TASK] Downloading $FILE_NAME"
  wget --continue --directory-prefix="$OUT_DIR" "$URL"
done

if (( DOWNLOADS_FOR_TASK == 0 )); then
  echo "No URLs assigned to task $TASK_ID"
fi
