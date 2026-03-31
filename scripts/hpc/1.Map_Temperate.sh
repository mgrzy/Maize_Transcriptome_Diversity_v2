#!/bin/bash -l
#SBATCH --job-name=bwa_temp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --account=g102-2488
#SBATCH --partition=okeanos
#SBATCH --error=logs/alignment/temp.%A_%a.err
#SBATCH --output=logs/alignment/temp.%A_%a.out
#SBATCH --array=1-30

module load common/anaconda/3.8
conda activate geno

set -euo pipefail

# Input fastq files and mapping list (Taxa <tab> SRR)
fastqPath="data/raw/fastq/temp"
mapFile="data/raw/meta/SrrToGeno/Temp.txt"

# Output paths
BAMpath="results/full/BAM/temperate"
repoPath="results/full/fastp"

# Reference genome
geno="../SNPcall/genome/maizeV5.fa.gz"
speedseq_bin="../soft/speedseq/bin/speedseq"

# Work directory
workDir="workdir"

mkdir -p "$BAMpath" "$repoPath" "$workDir"

wait_for_stable_outputs() {
  local sample_name="$1"
  shift
  local -a outputs=("$@")
  local attempt=0
  local max_attempts=60
  local sleep_seconds=10
  local previous_signature=""

  while (( attempt < max_attempts )); do
    local missing=0
    local current_signature=""

    for output in "${outputs[@]}"; do
      if [[ ! -s "$output" ]]; then
        missing=1
        break
      fi
      current_signature+=$(stat -c '%n:%s' "$output")
      current_signature+=$'\n'
    done

    if (( missing == 0 )); then
      if [[ "$current_signature" == "$previous_signature" ]]; then
        return 0
      fi
      previous_signature="$current_signature"
    fi

    sleep "$sleep_seconds"
    attempt=$((attempt + 1))
  done

  echo "Timed out waiting for stable speedseq outputs for $sample_name" >&2
  return 1
}

# Build sample list from mapping file
if [[ ! -f "$mapFile" ]]; then
  echo "Mapping file not found: $mapFile" >&2
  exit 1
fi

if [[ ! -d "$fastqPath" ]]; then
  echo "FASTQ directory not found: $fastqPath" >&2
  exit 1
fi

if [[ ! -f "$geno" ]]; then
  echo "Reference genome not found: $geno" >&2
  exit 1
fi

if [[ ! -x "$speedseq_bin" ]]; then
  echo "speedseq executable not found or not executable: $speedseq_bin" >&2
  exit 1
fi

# Read mapping into arrays preserving order
mapfile -t taxa_list < <(awk '{print $1}' "$mapFile" | sed '/^$/d')
mapfile -t srr_list < <(awk '{print $2}' "$mapFile" | sed '/^$/d')

if [[ ${#taxa_list[@]} -ne ${#srr_list[@]} ]]; then
  echo "Mapping file has inconsistent columns (Taxa vs SRR)" >&2
  exit 1
fi

total=${#taxa_list[@]}
if [[ $total -eq 0 ]]; then
  echo "No entries found in mapping file: $mapFile" >&2
  exit 1
fi

if [[ $(printf '%s\n' "${taxa_list[@]}" | sort | uniq -d | wc -l) -ne 0 ]]; then
  echo "Mapping file contains duplicate Taxa values; BAM outputs would collide" >&2
  exit 1
fi

if [[ $(printf '%s\n' "${srr_list[@]}" | sort | uniq -d | wc -l) -ne 0 ]]; then
  echo "Mapping file contains duplicate SRR values" >&2
  exit 1
fi

# Determine chunk size to process at least 2 lines per job
ARRAY_MAX=${SLURM_ARRAY_TASK_MAX:-30}
task_id=${SLURM_ARRAY_TASK_ID:-1}

# ceil(total / ARRAY_MAX)
chunk_size=$(( (total + ARRAY_MAX - 1) / ARRAY_MAX ))
if [[ $chunk_size -lt 2 ]]; then
  chunk_size=2
fi

start_idx=$(( (task_id - 1) * chunk_size ))
end_idx=$(( start_idx + chunk_size ))
if [[ $end_idx -gt $total ]]; then
  end_idx=$total
fi

if [[ $start_idx -ge $total ]]; then
  echo "No samples assigned to task $task_id/$ARRAY_MAX"
  conda deactivate
  exit 0
fi

echo "Processing samples $((start_idx+1)) to $end_idx of $total (task $task_id/$ARRAY_MAX; chunk_size=$chunk_size)"

for ((i=start_idx; i<end_idx; i++)); do
  name="${taxa_list[$i]}"
  srr="${srr_list[$i]}"
  r1="$fastqPath/${srr}_1.fastq.gz"
  r2="$fastqPath/${srr}_2.fastq.gz"
  sample_work_dir="$workDir/${name}"
  trimmed_r1="$sample_work_dir/${srr}_1.fastq.gz"
  trimmed_r2="$sample_work_dir/${srr}_2.fastq.gz"
  bam_prefix="$BAMpath/${name}"
  bam_file="${bam_prefix}.bam"
  bam_index="${bam_file}.bai"
  splitter_bam="${bam_prefix}.splitters.bam"
  splitter_index="${splitter_bam}.bai"
  discordant_bam="${bam_prefix}.discordants.bam"
  discordant_index="${discordant_bam}.bai"
  read_group=$(printf '@RG\tID:%s\tSM:%s\tLB:libX\tPL:ILLUMINA\tPU:unitX' "$name" "$name")

  if [[ ! -f "$r1" ]]; then
    echo "R1 missing for $name (SRR=$srr): $r1 not found, skipping" >&2
    continue
  fi
  if [[ ! -f "$r2" ]]; then
    echo "R2 missing for $name (SRR=$srr): $r2 not found, skipping" >&2
    continue
  fi

  echo "Running fastp for $name (SRR=$srr)"
  mkdir -p "$sample_work_dir"
  fastp \
    -i "$r1" -I "$r2" \
    -o "$trimmed_r1" -O "$trimmed_r2" \
    -w 12 -h "$repoPath/${name}.html" -j "$repoPath/${name}.json"

  if [[ ! -f "$trimmed_r1" || ! -f "$trimmed_r2" ]]; then
    echo "fastp did not produce expected outputs for $name (SRR=$srr)" >&2
    rm -rf "$sample_work_dir"
    continue
  fi

  echo "Aligning with speedseq for $name"
  "$speedseq_bin" align \
    -o "$bam_prefix" \
    -R "$read_group" \
    -t 12 \
    -T "$sample_work_dir" \
    "$geno" \
    "$trimmed_r1" \
    "$trimmed_r2"

  if ! wait_for_stable_outputs \
    "$name" \
    "$bam_file" \
    "$bam_index" \
    "$splitter_bam" \
    "$splitter_index" \
    "$discordant_bam" \
    "$discordant_index"; then
    echo "speedseq outputs incomplete for $name (SRR=$srr); keeping $sample_work_dir for inspection" >&2
    continue
  fi

  # Cleanup intermediates and extra outputs we don't keep
  rm -rf "$sample_work_dir"
  rm -f "$discordant_bam" "$discordant_index"
  rm -f "$splitter_bam" "$splitter_index"
done

conda deactivate
