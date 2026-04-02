#!/bin/bash -l
#SBATCH --job-name=bwa_trop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --account=g102-2488
#SBATCH --partition=okeanos
#SBATCH --error=logs/alignment/trop.%A_%a.err
#SBATCH --output=logs/alignment/trop.%A_%a.out
#SBATCH --array=1-15

module load common/anaconda/3.8
conda activate geno

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
project_root=$(cd "$submit_dir" && pwd)

if [[ ! -f "$project_root/data/raw/meta/fastq_links/Tropical.txt" && -f "$project_root/../../data/raw/meta/fastq_links/Tropical.txt" ]]; then
  project_root=$(cd "$project_root/../.." && pwd)
fi

fastqPath="$project_root/data/raw/fastq/tropical"
fastqList="$project_root/data/raw/meta/fastq_links/Tropical.txt"
filereport="$project_root/data/raw/meta/filereport_read_run_PRJNA389800.tsv"
BAMpath="$project_root/results/full/BAM/tropical"
repoPath="$project_root/results/full/fastp/tropical"
geno="$project_root/../SNPcall/genome/maizeV5.fa.gz"
toolsDir="$project_root/../soft/speedseq/bin"
bwaBin="$toolsDir/bwa"
samblasterBin="$toolsDir/samblaster"
sambambaBin="$toolsDir/sambamba"
workDir="$project_root/workdir/tropical"

mkdir -p "$BAMpath" "$repoPath" "$workDir" "$project_root/logs/alignment"
cd "$project_root"

resolve_executable() {
  local candidate="$1"
  local fallback="$2"

  if [[ -x "$candidate" ]]; then
    printf '%s\n' "$candidate"
    return 0
  fi

  if command -v "$fallback" >/dev/null 2>&1; then
    command -v "$fallback"
    return 0
  fi

  echo "Required executable not found: $candidate (or $fallback on PATH)" >&2
  return 1
}

derive_inbred_name() {
  local library_name="$1"
  local suffix=""

  while [[ "$library_name" == *_* ]]; do
    suffix="${library_name##*_}"
    if [[ "$suffix" =~ ^L[0-9]+$ || "$suffix" =~ ^[ACGTN]+$ || "$suffix" =~ ^H[A-Z0-9]+XX$ || "$suffix" =~ ^[A-Z0-9]{8,}$ ]]; then
      library_name="${library_name%_*}"
    else
      break
    fi
  done

  printf '%s\n' "$library_name"
}

run_alignment_pipeline() {
  local read_group="$1"
  local sample_work_dir="$2"
  local bam_file="$3"
  local splitter_bam="$4"
  local discordant_bam="$5"
  local trimmed_r1="$6"
  local trimmed_r2="$7"

  local full_tmp_dir="$sample_work_dir/full"
  local split_tmp_dir="$sample_work_dir/spl"
  local disc_tmp_dir="$sample_work_dir/disc"
  local splitter_pipe="$sample_work_dir/spl_pipe"
  local discordant_pipe="$sample_work_dir/disc_pipe"
  local split_pid
  local disc_pid
  local status=0

  mkdir -p "$full_tmp_dir" "$split_tmp_dir" "$disc_tmp_dir"
  rm -f "$splitter_pipe" "$discordant_pipe"
  mkfifo "$splitter_pipe" "$discordant_pipe"

  awk '{ if ($0 ~ "^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" "$splitter_pipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 4 -m 1G --tmpdir="$split_tmp_dir" -o "$splitter_bam" /dev/stdin &
  split_pid=$!

  awk '{ if ($0 ~ "^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" "$discordant_pipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 4 -m 1G --tmpdir="$disc_tmp_dir" -o "$discordant_bam" /dev/stdin &
  disc_pid=$!

  if ! "$bwaBin" mem -t 12 -R "$read_group" "$geno" "$trimmed_r1" "$trimmed_r2" \
    | "$samblasterBin" --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
        --splitterFile "$splitter_pipe" --discordantFile "$discordant_pipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 12 -m 18G --tmpdir="$full_tmp_dir" -o "$bam_file" /dev/stdin; then
    status=1
  fi

  if ! wait "$split_pid"; then
    status=1
  fi

  if ! wait "$disc_pid"; then
    status=1
  fi

  rm -f "$splitter_pipe" "$discordant_pipe"

  if (( status != 0 )); then
    return "$status"
  fi

  "$sambambaBin" index "$bam_file"
}

bwaBin=$(resolve_executable "$bwaBin" bwa)
samblasterBin=$(resolve_executable "$samblasterBin" samblaster)
sambambaBin=$(resolve_executable "$sambambaBin" sambamba)

if ! command -v fastp >/dev/null 2>&1; then
  echo "Required executable not found on PATH: fastp" >&2
  exit 1
fi

if [[ ! -f "$fastqList" ]]; then
  echo "FASTQ list not found: $fastqList" >&2
  exit 1
fi

if [[ ! -f "$filereport" ]]; then
  echo "Filereport not found: $filereport" >&2
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

manifest_file=$(mktemp "$workDir/tropical_manifest.XXXXXX.tsv")
groups_file=$(mktemp "$workDir/tropical_groups.XXXXXX.tsv")
trap 'rm -f "$manifest_file" "$groups_file"' EXIT

awk -F '\t' '
  NR == FNR {
    n = split($0, path_parts, "/")
    file_name = path_parts[n]
    sub(/_[12]\.fastq\.gz$/, "", file_name)
    if (file_name != "") {
      runs[file_name] = 1
    }
    next
  }
  FNR == 1 {
    next
  }
  ($1 in runs) {
    print $1 "\t" $3 "\t" $7
  }
' "$fastqList" "$filereport" | sort -u > "$manifest_file"

if [[ ! -s "$manifest_file" ]]; then
  echo "No tropical runs from $fastqList were matched in $filereport" >&2
  exit 1
fi

declare -A sample_inbred=()
declare -A sample_runs=()
declare -A sample_seen=()
declare -A inbred_seen=()
declare -a sample_order=()

while IFS=$'\t' read -r run_accession sample_accession library_name; do
  if [[ -z "$run_accession" || -z "$sample_accession" || -z "$library_name" ]]; then
    echo "Malformed manifest row: $run_accession $sample_accession $library_name" >&2
    exit 1
  fi

  inbred_name=$(derive_inbred_name "$library_name")
  if [[ -z "$inbred_name" ]]; then
    echo "Failed to derive inbred name from library name: $library_name" >&2
    exit 1
  fi

  if [[ -z "${sample_seen[$sample_accession]:-}" ]]; then
    sample_order+=("$sample_accession")
    sample_seen["$sample_accession"]=1
  fi

  if [[ -n "${sample_inbred[$sample_accession]:-}" && "${sample_inbred[$sample_accession]}" != "$inbred_name" ]]; then
    echo "Inconsistent inbred names for $sample_accession: ${sample_inbred[$sample_accession]} vs $inbred_name" >&2
    exit 1
  fi
  sample_inbred["$sample_accession"]="$inbred_name"
  sample_runs["$sample_accession"]+="${sample_runs[$sample_accession]:+ }$run_accession"
done < "$manifest_file"

for sample_accession in "${sample_order[@]}"; do
  inbred_name="${sample_inbred[$sample_accession]}"
  if [[ -n "${inbred_seen[$inbred_name]:-}" && "${inbred_seen[$inbred_name]}" != "$sample_accession" ]]; then
    echo "Duplicate inbred name resolved for multiple samples: $inbred_name" >&2
    exit 1
  fi
  inbred_seen["$inbred_name"]="$sample_accession"
  printf '%s\t%s\t%s\n' "$inbred_name" "$sample_accession" "${sample_runs[$sample_accession]}"
done | sort -k1,1 > "$groups_file"

mapfile -t group_lines < "$groups_file"
total_inbreds=${#group_lines[@]}
task_id=${SLURM_ARRAY_TASK_ID:-1}

if [[ $total_inbreds -eq 0 ]]; then
  echo "No tropical inbreds were derived from metadata" >&2
  exit 1
fi

if [[ $task_id -gt $total_inbreds ]]; then
  echo "No inbred assigned to task $task_id; total inbreds: $total_inbreds"
  conda deactivate
  exit 0
fi

IFS=$'\t' read -r inbred_name sample_accession run_string <<< "${group_lines[$((task_id - 1))]}"
read -r -a run_accessions <<< "$run_string"

if [[ ${#run_accessions[@]} -eq 0 ]]; then
  echo "No runs assigned to inbred $inbred_name ($sample_accession)" >&2
  exit 1
fi

echo "Task $task_id/$total_inbreds"
echo "Processing tropical inbred $inbred_name ($sample_accession) with ${#run_accessions[@]} libraries"

sample_work_dir="$workDir/$inbred_name"
trimmed_dir="$sample_work_dir/trimmed"
aligned_dir="$sample_work_dir/aligned"
mkdir -p "$trimmed_dir" "$aligned_dir"

run_bams=()

for run_accession in "${run_accessions[@]}"; do
  r1="$fastqPath/${run_accession}_1.fastq.gz"
  r2="$fastqPath/${run_accession}_2.fastq.gz"
  trimmed_r1="$trimmed_dir/${run_accession}_1.fastq.gz"
  trimmed_r2="$trimmed_dir/${run_accession}_2.fastq.gz"
  run_prefix="$aligned_dir/${run_accession}"
  run_bam="${run_prefix}.bam"
  run_splitters="${run_prefix}.splitters.bam"
  run_discordants="${run_prefix}.discordants.bam"
  read_group=$(printf '@RG\tID:%s\tSM:%s\tLB:%s\tPL:ILLUMINA\tPU:%s' "$run_accession" "$inbred_name" "$sample_accession" "$run_accession")

  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "FASTQ files missing for $run_accession in $fastqPath" >&2
    exit 1
  fi

  echo "Running fastp for $inbred_name library $run_accession"
  fastp \
    -i "$r1" -I "$r2" \
    -o "$trimmed_r1" -O "$trimmed_r2" \
    -w 12 \
    -h "$repoPath/${inbred_name}.${run_accession}.html" \
    -j "$repoPath/${inbred_name}.${run_accession}.json"

  if [[ ! -s "$trimmed_r1" || ! -s "$trimmed_r2" ]]; then
    echo "fastp did not create expected outputs for $run_accession" >&2
    exit 1
  fi

  echo "Aligning $run_accession into $run_bam"
  run_alignment_pipeline \
    "$read_group" \
    "$sample_work_dir/$run_accession" \
    "$run_bam" \
    "$run_splitters" \
    "$run_discordants" \
    "$trimmed_r1" \
    "$trimmed_r2"

  if [[ ! -s "$run_bam" ]]; then
    echo "Alignment finished but BAM was not created: $run_bam" >&2
    exit 1
  fi

  rm -rf "$sample_work_dir/$run_accession"
  rm -f "$run_splitters" "$run_splitters.bai"
  rm -f "$run_discordants" "$run_discordants.bai"
  rm -f "$trimmed_r1" "$trimmed_r2"

  run_bams+=("$run_bam")
done

final_bam="$BAMpath/${inbred_name}.bam"
final_bam_tmp="$aligned_dir/${inbred_name}.bam"
rm -f "$final_bam" "$final_bam.bai" "$final_bam_tmp" "$final_bam_tmp.bai"

if [[ ${#run_bams[@]} -eq 1 ]]; then
  mv "${run_bams[0]}" "$final_bam_tmp"
  if [[ -f "${run_bams[0]}.bai" ]]; then
    mv "${run_bams[0]}.bai" "${final_bam_tmp}.bai"
  fi
else
  echo "Merging ${#run_bams[@]} BAMs into $final_bam"
  "$sambambaBin" merge -t 12 "$final_bam_tmp" "${run_bams[@]}"
  "$sambambaBin" index "$final_bam_tmp"
  rm -f "${run_bams[@]}" "${run_bams[@]/%/.bai}"
fi

mv "$final_bam_tmp" "$final_bam"
if [[ -f "${final_bam_tmp}.bai" ]]; then
  mv "${final_bam_tmp}.bai" "${final_bam}.bai"
fi

if [[ ! -s "$final_bam" ]]; then
  echo "Final BAM missing after merge: $final_bam" >&2
  exit 1
fi

rm -rf "$sample_work_dir"

echo "Finished tropical inbred $inbred_name"

conda deactivate
