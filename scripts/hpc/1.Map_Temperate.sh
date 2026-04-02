#!/bin/bash -l
#SBATCH --job-name=bwa_temp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --account=g102-2488
#SBATCH --partition=okeanos
#SBATCH --error=logs/alignment/temp.%A_%a.err
#SBATCH --output=logs/alignment/temp.%A_%a.out
#SBATCH --array=1-26

module load common/anaconda/3.8
conda activate geno

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
project_root=$(cd "$submit_dir" && pwd)

if [[ ! -f "$project_root/data/raw/meta/SrrToGeno/Temp.txt" && -f "$project_root/../../data/raw/meta/SrrToGeno/Temp.txt" ]]; then
  project_root=$(cd "$project_root/../.." && pwd)
fi

# Paths
fastqPath="$project_root/data/raw/fastq/temp"
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

total_samples=$(awk 'NF > 0' "$mapFile" | wc -l)
samples_per_job=3

start_line=$(( (${SLURM_ARRAY_TASK_ID:-1} - 1) * samples_per_job + 1 ))
end_line=$(( start_line + samples_per_job - 1 ))

if [[ $start_line -gt $total_samples ]]; then
  echo "No samples assigned to task $SLURM_ARRAY_TASK_ID"
  conda deactivate
  exit 0
fi

if [[ $end_line -gt $total_samples ]]; then
  end_line=$total_samples
fi

echo "Processing samples $start_line to $end_line of $total_samples"

for line_no in $(seq "$start_line" "$end_line"); do
  samp=$(sed -n "${line_no}p" "$mapFile" | awk '{print $2}')
  Taxa=$(sed -n "${line_no}p" "$mapFile" | awk '{print $1}')

  if [[ -z "$samp" || -z "$Taxa" ]]; then
    echo "Malformed line $line_no in $mapFile" >&2
    conda deactivate
    exit 1
  fi

  echo "Processing sample $line_no of $total_samples: $Taxa ($samp)"

  if [[ ! -f "${fastqPath}/${samp}_1.fastq.gz" || ! -f "${fastqPath}/${samp}_2.fastq.gz" ]]; then
    echo "FASTQ files missing for $Taxa ($samp) in $fastqPath" >&2
    conda deactivate
    exit 1
  fi

  fastp \
    -i "${fastqPath}/${samp}_1.fastq.gz" \
    -I "${fastqPath}/${samp}_2.fastq.gz" \
    -o "${workDir}/${Taxa}_1.fastq.gz" \
    -O "${workDir}/${Taxa}_2.fastq.gz" \
    -w 12 \
    -h "${repoPath}/${Taxa}.html" \
    -j "${repoPath}/${Taxa}.json"

  sampleTmpDir="${workDir}/${Taxa}"
  fullTmpDir="${sampleTmpDir}/full"
  splTmpDir="${sampleTmpDir}/spl"
  discTmpDir="${sampleTmpDir}/disc"
  splPipe="${sampleTmpDir}/spl_pipe"
  discPipe="${sampleTmpDir}/disc_pipe"
  rgString=$(printf '@RG\tID:%s\tSM:%s\tLB:libX\tPL:ILLUMINA\tPU:unitX' "$Taxa" "$Taxa")

  mkdir -p "$fullTmpDir" "$splTmpDir" "$discTmpDir"
  rm -f "$splPipe" "$discPipe"
  mkfifo "$splPipe" "$discPipe"

  gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" "$splPipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 4 -m 1G --tmpdir="$splTmpDir" -o "${BAMpath}/${Taxa}.splitters.bam" /dev/stdin &
  spl_pid=$!

  gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" "$discPipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 4 -m 1G --tmpdir="$discTmpDir" -o "${BAMpath}/${Taxa}.discordants.bam" /dev/stdin &
  disc_pid=$!

  "$bwaBin" mem -t 12 -R "$rgString" "$geno" "${workDir}/${Taxa}_1.fastq.gz" "${workDir}/${Taxa}_2.fastq.gz" \
    | "$samblasterBin" --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
        --splitterFile "$splPipe" --discordantFile "$discPipe" \
    | "$sambambaBin" view -S -f bam -l 0 /dev/stdin \
    | "$sambambaBin" sort -t 12 -m 18G --tmpdir="$fullTmpDir" -o "${BAMpath}/${Taxa}.bam" /dev/stdin

  wait "$spl_pid"
  wait "$disc_pid"

  rm -f "$splPipe" "$discPipe"

  "$sambambaBin" index "${BAMpath}/${Taxa}.bam"
  "$sambambaBin" index "${BAMpath}/${Taxa}.discordants.bam"
  "$sambambaBin" index "${BAMpath}/${Taxa}.splitters.bam"

  if [[ ! -f "${BAMpath}/${Taxa}.bam" ]]; then
    echo "Alignment finished but no BAM was created: ${BAMpath}/${Taxa}.bam" >&2
    conda deactivate
    exit 1
  fi

  rm -rf "${workDir}/${Taxa}"
  rm -f "${BAMpath}/${Taxa}.discordants.bam"
  rm -f "${BAMpath}/${Taxa}.discordants.bam.bai"
  rm -f "${BAMpath}/${Taxa}.splitters.bam"
  rm -f "${BAMpath}/${Taxa}.splitters.bam.bai"
  rm -f "${workDir}/${Taxa}_"*
done

conda deactivate
