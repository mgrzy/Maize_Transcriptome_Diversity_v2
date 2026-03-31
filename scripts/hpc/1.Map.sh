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

# Paths
fastqPath="data/raw/fastq/temp"
BAMpath="results/full/BAM/temperate"
repoPath="results/full/fastp"
geno="../SNPcall/genome/maizeV5.fa.gz"
workDir="workdir"
mapFile="data/raw/meta/SrrToGeno/Temp.txt"
speedseqBin="../soft/speedseq/bin/speedseq"

samples_per_job=3

mkdir -p "$BAMpath" "$repoPath" "$workDir" "logs/alignment"

if [[ ! -f "$mapFile" ]]; then
  echo "Mapping file not found: $mapFile" >&2
  exit 1
fi

total_samples=$(awk 'NF > 0' "$mapFile" | wc -l)
start_line=$(( (SLURM_ARRAY_TASK_ID - 1) * samples_per_job + 1 ))
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
    echo "Skipping malformed line $line_no in $mapFile" >&2
    continue
  fi

  echo "$samp"

  fastp \
    -i "${fastqPath}/${samp}_1.fastq.gz" \
    -I "${fastqPath}/${samp}_2.fastq.gz" \
    -o "${workDir}/${Taxa}_1.fastq.gz" \
    -O "${workDir}/${Taxa}_2.fastq.gz" \
    -w 12 \
    -h "${repoPath}/${Taxa}.html" \
    -j "${repoPath}/${Taxa}.json"

  "$speedseqBin" align \
    -o "${BAMpath}/${Taxa}" \
    -R "$(printf '@RG\tID:%s\tSM:%s\tLB:libX\tPL:ILLUMINA\tPU:unitX' "$Taxa" "$Taxa")" \
    -t 12 \
    -T "${workDir}/${Taxa}" \
    "$geno" \
    "${workDir}/${Taxa}_1.fastq.gz" \
    "${workDir}/${Taxa}_2.fastq.gz"

  rm -rf "${workDir}/${Taxa}"
  rm -f "${BAMpath}/${Taxa}.discordants.bam"
  rm -f "${BAMpath}/${Taxa}.discordants.bam.bai"
  rm -f "${BAMpath}/${Taxa}.splitters.bam"
  rm -f "${BAMpath}/${Taxa}.splitters.bam.bai"
  rm -f "${workDir}/${Taxa}_"*
done

conda deactivate
