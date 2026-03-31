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

# Work directory
workDir="workdir"

mkdir -p "$BAMpath" "$repoPath" "$workDir"

# Build sample list from mapping file
if [[ ! -f "$mapFile" ]]; then
  echo "Mapping file not found: $mapFile" >&2
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

echo "Processing samples $((start_idx+1)) to $end_idx of $total (task $task_id/$ARRAY_MAX; chunk_size=$chunk_size)"

for ((i=start_idx; i<end_idx; i++)); do
  name="${taxa_list[$i]}"
  srr="${srr_list[$i]}"
  r1="$fastqPath/${srr}_1.fastq.gz"
  r2="$fastqPath/${srr}_2.fastq.gz"
  sample_work_dir="$workDir/${name}"
  trimmed_r1="$sample_work_dir/${srr}_1.fastq.gz"
  trimmed_r2="$sample_work_dir/${srr}_2.fastq.gz"
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
  ../soft/speedseq/bin/speedseq align \
    -o "$BAMpath/${name}" \
    -R "$read_group" \
    -t 12 \
    -T "$sample_work_dir" \
    "$geno" \
    "$trimmed_r1" \
    "$trimmed_r2"

  # Cleanup intermediates and extra outputs we don't keep
  rm -rf "$sample_work_dir"
  rm -f "$BAMpath/${name}.discordants.bam" "$BAMpath/${name}.discordants.bam.bai"
  rm -f "$BAMpath/${name}.splitters.bam" "$BAMpath/${name}.splitters.bam.bai"
done

conda deactivate
