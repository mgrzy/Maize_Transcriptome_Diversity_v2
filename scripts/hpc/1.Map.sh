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
#SBATCH --array=1-2

module load common/anaconda/3.8
conda activate geno

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)

# Paths
fastqPath="$project_root/data/raw/fastq/temp"
BAMpath="$project_root/results/full/BAM/temperate"
repoPath="$project_root/results/full/fastp"
geno="$project_root/../SNPcall/genome/maizeV5.fa.gz"
workDir="$project_root/workdir"
mapFile="$project_root/data/raw/meta/SrrToGeno/Temp.txt"
speedseqBin="$project_root/../soft/speedseq/bin/speedseq"

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

if [[ ! -f "$geno" ]]; then
  echo "Reference genome not found: $geno" >&2
  exit 1
fi

total_samples=$(awk 'NF > 0' "$mapFile" | wc -l)

if [[ ${SLURM_ARRAY_TASK_ID:-1} -gt $total_samples ]]; then
  echo "No samples assigned to task $SLURM_ARRAY_TASK_ID"
  conda deactivate
  exit 0
fi

line_no=${SLURM_ARRAY_TASK_ID:-1}
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

"$speedseqBin" align \
  -o "${BAMpath}/${Taxa}" \
  -R "$(printf '@RG\tID:%s\tSM:%s\tLB:libX\tPL:ILLUMINA\tPU:unitX' "$Taxa" "$Taxa")" \
  -t 12 \
  -T "${workDir}/${Taxa}" \
  "$geno" \
  "${workDir}/${Taxa}_1.fastq.gz" \
  "${workDir}/${Taxa}_2.fastq.gz"

if [[ ! -f "${BAMpath}/${Taxa}.bam" ]]; then
  echo "speedseq finished but no BAM was created: ${BAMpath}/${Taxa}.bam" >&2
  conda deactivate
  exit 1
fi

rm -rf "${workDir}/${Taxa}"
rm -f "${BAMpath}/${Taxa}.discordants.bam"
rm -f "${BAMpath}/${Taxa}.discordants.bam.bai"
rm -f "${BAMpath}/${Taxa}.splitters.bam"
rm -f "${BAMpath}/${Taxa}.splitters.bam.bai"
rm -f "${workDir}/${Taxa}_"*

conda deactivate
