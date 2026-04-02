#!/bin/bash -l
#SBATCH --job-name=bwa_temp_fix
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --account=g102-2488
#SBATCH --partition=okeanos
#SBATCH --error=logs/alignment/temp_fix.%A_%a.err
#SBATCH --output=logs/alignment/temp_fix.%A_%a.out
#SBATCH --array=1-4

module load common/anaconda/3.8
conda activate geno

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
project_root=$(cd "$submit_dir" && pwd)

if [[ ! -f "$project_root/data/raw/meta/fastq_links/Temperate2.txt" && -f "$project_root/../../data/raw/meta/fastq_links/Temperate2.txt" ]]; then
  project_root=$(cd "$project_root/../.." && pwd)
fi

fastqPath="$project_root/data/raw/fastq/temp"
fastqList="$project_root/data/raw/meta/fastq_links/Temperate2.txt"
BAMpath="$project_root/results/full/BAM/temperate"
repoPath="$project_root/results/full/fastp"
geno="$project_root/../SNPcall/genome/maizeV5.fa.gz"
workDir="$project_root/workdir"
mapFile="$project_root/data/raw/meta/SrrToGeno/Temp.txt"
speedseqBin="$project_root/../soft/speedseq/bin/speedseq"
speedseqDir=$(dirname "$speedseqBin")
bwaBin="$speedseqDir/bwa"
samblasterBin="$speedseqDir/samblaster"
sambambaBin="$speedseqDir/sambamba"

mkdir -p "$BAMpath" "$repoPath" "$workDir" "$project_root/logs/alignment"
cd "$project_root"

if [[ ! -f "$fastqList" ]]; then
  echo "FASTQ list not found: $fastqList" >&2
  exit 1
fi

if [[ ! -f "$mapFile" ]]; then
  echo "Mapping file not found: $mapFile" >&2
  exit 1
fi

if [[ ! -x "$speedseqBin" ]]; then
  echo "speedseq executable not found: $speedseqBin" >&2
  exit 1
fi

if [[ ! -x "$bwaBin" ]]; then
  echo "bwa executable not found: $bwaBin" >&2
  exit 1
fi

if [[ ! -x "$samblasterBin" ]]; then
  echo "samblaster executable not found: $samblasterBin" >&2
  exit 1
fi

if [[ ! -x "$sambambaBin" ]]; then
  echo "sambamba executable not found: $sambambaBin" >&2
  exit 1
fi

if [[ ! -f "$geno" ]]; then
  echo "Reference genome not found: $geno" >&2
  exit 1
fi

manifest_file=$(mktemp "$workDir/temperate_redownload_manifest.XXXXXX.tsv")
trap 'rm -f "$manifest_file"; conda deactivate' EXIT

awk '
  NR == FNR {
    n = split($0, path_parts, "/")
    file_name = path_parts[n]
    sub(/_[12]\.fastq\.gz$/, "", file_name)
    if (file_name != "" && !seen[file_name]++) {
      order[++count] = file_name
    }
    next
  }
  FNR == 1 {
    next
  }
  ($2 in seen) {
    taxa[$2] = $1
  }
  END {
    missing = 0
    for (i = 1; i <= count; i++) {
      run = order[i]
      if (!(run in taxa)) {
        print "Run not found in mapping file: " run > "/dev/stderr"
        missing = 1
        continue
      }
      print taxa[run] "\t" run
    }
    if (missing) {
      exit 1
    }
  }
' "$fastqList" "$mapFile" > "$manifest_file"

mapfile -t selected_samples < "$manifest_file"
total_samples=${#selected_samples[@]}

if (( total_samples == 0 )); then
  echo "No temperate samples were derived from $fastqList" >&2
  exit 1
fi

process_sample() {
  local taxa="$1"
  local samp="$2"
  local sample_work_dir
  local full_tmp_dir
  local spl_tmp_dir
  local disc_tmp_dir
  local spl_pipe
  local disc_pipe
  local rgString
  local spl_pid
  local disc_pid

  echo "Processing $taxa ($samp)"

  if [[ ! -f "${fastqPath}/${samp}_1.fastq.gz" || ! -f "${fastqPath}/${samp}_2.fastq.gz" ]]; then
    echo "FASTQ files missing for $taxa ($samp) in $fastqPath" >&2
    exit 1
  fi

  fastp \
    -i "${fastqPath}/${samp}_1.fastq.gz" \
    -I "${fastqPath}/${samp}_2.fastq.gz" \
    -o "${workDir}/${taxa}_1.fastq.gz" \
    -O "${workDir}/${taxa}_2.fastq.gz" \
    -w 12 \
    -h "${repoPath}/${taxa}.html" \
    -j "${repoPath}/${taxa}.json"

  sample_work_dir="${workDir}/${taxa}"
  full_tmp_dir="${sample_work_dir}/full"
  spl_tmp_dir="${sample_work_dir}/spl"
  disc_tmp_dir="${sample_work_dir}/disc"
  spl_pipe="${sample_work_dir}/spl_pipe"
  disc_pipe="${sample_work_dir}/disc_pipe"
  rgString=$(printf '@RG\tID:%s\tSM:%s\tLB:libX\tPL:ILLUMINA\tPU:unitX' "$taxa" "$taxa")

  mkdir -p "$full_tmp_dir" "$spl_tmp_dir" "$disc_tmp_dir"
  rm -f "$spl_pipe" "$disc_pipe"
  mkfifo "$spl_pipe" "$disc_pipe"

  gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" "$spl_pipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 4 -m 1G --tmpdir="$spl_tmp_dir" -o "${BAMpath}/${taxa}.splitters.bam" /dev/stdin &
  spl_pid=$!

  gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" "$disc_pipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 4 -m 1G --tmpdir="$disc_tmp_dir" -o "${BAMpath}/${taxa}.discordants.bam" /dev/stdin &
  disc_pid=$!

  "$bwaBin" mem -t 12 -R "$rgString" "$geno" "${workDir}/${taxa}_1.fastq.gz" "${workDir}/${taxa}_2.fastq.gz" \
    | "$samblasterBin" --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
        --splitterFile "$spl_pipe" --discordantFile "$disc_pipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 12 -m 18G --tmpdir="$full_tmp_dir" -o "${BAMpath}/${taxa}.bam" /dev/stdin

  wait "$spl_pid"
  wait "$disc_pid"

  rm -f "$spl_pipe" "$disc_pipe"

  "$sambambaBin" index "${BAMpath}/${taxa}.bam"
  "$sambambaBin" index "${BAMpath}/${taxa}.discordants.bam"
  "$sambambaBin" index "${BAMpath}/${taxa}.splitters.bam"

  if [[ ! -f "${BAMpath}/${taxa}.bam" ]]; then
    echo "Alignment finished but no BAM was created: ${BAMpath}/${taxa}.bam" >&2
    exit 1
  fi

  rm -rf "${workDir}/${taxa}"
  rm -f "${BAMpath}/${taxa}.discordants.bam"
  rm -f "${BAMpath}/${taxa}.discordants.bam.bai"
  rm -f "${BAMpath}/${taxa}.splitters.bam"
  rm -f "${BAMpath}/${taxa}.splitters.bam.bai"
  rm -f "${workDir}/${taxa}_"*
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  task_id="$SLURM_ARRAY_TASK_ID"
  if (( task_id > total_samples )); then
    echo "No sample assigned to task $task_id; total selected samples: $total_samples"
    exit 0
  fi

  IFS=$'\t' read -r taxa samp <<< "${selected_samples[$((task_id - 1))]}"
  echo "Task $task_id/$total_samples"
  process_sample "$taxa" "$samp"
else
  echo "No SLURM array task detected; processing all $total_samples selected samples sequentially"
  for selected_sample in "${selected_samples[@]}"; do
    IFS=$'\t' read -r taxa samp <<< "$selected_sample"
    process_sample "$taxa" "$samp"
  done
fi
