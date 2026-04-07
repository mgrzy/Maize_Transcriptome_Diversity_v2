#!/bin/bash -l
#SBATCH --job-name=angsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --account=g102-2488
#SBATCH --partition=okeanos
#SBATCH --error=logs/angsd/pi%x.%j.err
#SBATCH --output=logs/angsd/pi%x.%j.out

module load common/anaconda/3.8
conda activate geno

set -euo pipefail

# Slurm can execute the job from a spool directory, so we cannot rely on the
# script path alone to locate the repository root.
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

# Core tool and data locations for the ANGSD -> realSFS -> thetaStat workflow.
ANGSD_ROOT="${PROJECT_ROOT}/../soft/angsd"
ANGSD="${ANGSD_ROOT}/angsd"
REALSFS="${ANGSD_ROOT}/misc/realSFS"
THETASTAT="${ANGSD_ROOT}/misc/thetaStat"
BEDTOOLS="${BEDTOOLS:-$(command -v bedtools || true)}"
GENOME="${PROJECT_ROOT}/../SNPcall/genome/maizeV5.fa.gz"
REGION_DIR="${PROJECT_ROOT}/data/processed/bed_windows"
OUTBASE="${PROJECT_ROOT}/results/full/pop/angsd"

THREADS_DEFAULT="${SLURM_CPUS_PER_TASK:-12}"
THREADS_ANGSD="${THREADS_ANGSD:-${THREADS_DEFAULT}}"
THREADS_SFS="${THREADS_SFS:-${THREADS_DEFAULT}}"
THREADS_THETA="${THREADS_THETA:-${THREADS_DEFAULT}}"

usage() {
  echo "Usage: $(basename "$0") <bamlist>" >&2
  echo "Example: $(basename "$0") data/raw/meta/Tropical_BAM.txt" >&2
}

resolve_path() {
  local path="$1"

  if [[ "${path}" = /* ]]; then
    printf '%s\n' "${path}"
  else
    printf '%s\n' "${PROJECT_ROOT}/${path}"
  fi
}

require_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
}

require_exec() {
  local path="$1"
  if [[ ! -x "${path}" ]]; then
    echo "Missing required executable: ${path}" >&2
    exit 1
  fi
}

is_newer_than() {
  local target="$1"
  local dependency="$2"

  [[ ! -e "${target}" || "${dependency}" -nt "${target}" ]]
}

any_dependency_newer() {
  local target="$1"
  shift
  local dependency

  for dependency in "$@"; do
    if is_newer_than "${target}" "${dependency}"; then
      return 0
    fi
  done

  return 1
}

# A chromosome-level theta file is considered reusable only when it exists,
# is newer than the shared theta inputs, and matches the new streamed
# `thetaStat print` format. This prevents reusing partial files from older
# OOM-failed runs.
theta_output_is_current() {
  local theta_file="$1"
  local union_prefix="$2"

  if [[ ! -s "${theta_file}" ]]; then
    return 1
  fi

  if any_dependency_newer "${theta_file}" "${union_prefix}.thetas.idx" "${union_prefix}.thetas.gz"; then
    return 1
  fi

  if ! head -n 1 "${theta_file}" | grep -q $'^#Chromo\tPos\t'; then
    return 1
  fi

  return 0
}

# Emit the chromosome names present in a BED file, preserving first-seen order.
chrom_stream() {
  local bed_file="$1"

  awk '!seen[$1]++ { print $1 }' "${bed_file}"
}

# ANGSD requires an indexed sites file when `-sites` is used. Rebuild it only
# if the merged BED changed.
index_sites_file() {
  local bed_file="$1"

  if [[ -f "${bed_file}.idx" && -f "${bed_file}.bin" ]] \
    && [[ "${bed_file}" -ot "${bed_file}.idx" ]] \
    && [[ "${bed_file}" -ot "${bed_file}.bin" ]]; then
    echo "Reusing indexed sites file: ${bed_file}"
    return
  fi

  echo "Indexing sites file: ${bed_file}"
  "${ANGSD}" sites index "${bed_file}"
}

# Merge all BED region sets into one union BED. ANGSD is run once on the union,
# then the chromosome-level outputs are split back into the original BED sets.
make_union_bed() {
  local merged_bed="$1"
  shift
  local region_files=( "$@" )
  local tmp_bed

  tmp_bed=$(mktemp "${OUTDIR}/${POP_NAME}.all_region_sets.merged.XXXXXX.bed")

  cat "${region_files[@]}" \
    | "${BEDTOOLS}" sort -i - \
    | "${BEDTOOLS}" merge -i - > "${tmp_bed}"

  if [[ -f "${merged_bed}" ]] && cmp -s "${tmp_bed}" "${merged_bed}"; then
    rm -f "${tmp_bed}"
    echo "Reusing merged BED: ${merged_bed}"
    return
  fi

  mv "${tmp_bed}" "${merged_bed}"
  echo "Updated merged BED: ${merged_bed}"
}

# Run the expensive shared steps once per population:
# 1. `angsd -doSaf` on all sites in the merged BED
# 2. `realSFS` to estimate the folded SFS
# 3. `realSFS saf2theta` to generate theta files used downstream
# These outputs are reused on reruns unless any input changed.
run_union_pipeline() {
  local union_bed="$1"
  local union_prefix="$2"
  local union_sfs_file union_saf_idx union_saf_gz union_saf_pos union_thetas_idx union_thetas_gz

  union_sfs_file="${union_prefix}.sfs"
  union_saf_idx="${union_prefix}.saf.idx"
  union_saf_gz="${union_prefix}.saf.gz"
  union_saf_pos="${union_prefix}.saf.pos.gz"
  union_thetas_idx="${union_prefix}.thetas.idx"
  union_thetas_gz="${union_prefix}.thetas.gz"

  if ! any_dependency_newer \
    "${union_thetas_idx}" \
    "${union_bed}" "${union_bed}.idx" "${union_bed}.bin" "${BAMLIST}" "${GENOME}" \
    && [[ -s "${union_saf_idx}" && -s "${union_saf_gz}" && -s "${union_saf_pos}" && -s "${union_sfs_file}" && -s "${union_thetas_gz}" ]]; then
    echo "Reusing union ANGSD outputs: ${union_prefix}"
    return
  fi

  echo "Running ANGSD pipeline for ${POP_NAME}"

  "${ANGSD}" -bam "${BAMLIST}" -sites "${union_bed}" -doSaf 1 \
    -anc "${GENOME}" -GL 1 -P "${THREADS_ANGSD}" -out "${union_prefix}"

  "${REALSFS}" "${union_prefix}.saf.idx" -P "${THREADS_SFS}" -fold 1 > "${union_sfs_file}"
  "${REALSFS}" saf2theta "${union_prefix}.saf.idx" -outname "${union_prefix}" -P "${THREADS_THETA}" -fold 1 -sfs "${union_sfs_file}"
}

# Generate one chromosome-level theta table at a time. This is intentionally
# serial so reruns can resume from already completed chromosomes. We use
# `thetaStat print` instead of `do_stat -win 1 -step 1` to avoid OOM on large
# chromosomes.
run_theta_by_chromosome() {
  local union_bed="$1"
  local union_prefix="$2"
  local chrom theta_prefix

  while IFS= read -r chrom; do
    [[ -z "${chrom}" ]] && continue
    theta_prefix="${union_prefix}_${chrom}"
    if theta_output_is_current "${theta_prefix}.pestPG" "${union_prefix}"; then
      echo "Reusing thetaStat output: ${theta_prefix}.pestPG"
      continue
    fi
    # thetaStat do_stat with -win 1 -step 1 accumulates a chromosome-wide string in memory.
    # thetaStat print streams per-site theta values instead and avoids OOM on large chromosomes.
    "${THETASTAT}" print "${union_prefix}.thetas.idx" -r "${chrom}" > "${theta_prefix}.pestPG"
  done < <(chrom_stream "${union_bed}")
}

# Check whether a BED-specific split file can be reused. It must be present,
# newer than the BED and theta inputs, and newer than every chromosome file
# it depends on.
split_output_is_current() {
  local region_out="$1"
  local bed_file="$2"
  local union_prefix="$3"
  local chrom theta_pest

  if [[ ! -s "${region_out}" ]]; then
    return 1
  fi

  if any_dependency_newer "${region_out}" "${bed_file}" "${union_prefix}.thetas.idx"; then
    return 1
  fi

  while IFS= read -r chrom; do
    [[ -z "${chrom}" ]] && continue
    theta_pest="${union_prefix}_${chrom}.pestPG"
    if [[ ! -f "${theta_pest}" || "${region_out}" -ot "${theta_pest}" ]]; then
      return 1
    fi
  done < <(chrom_stream "${bed_file}")

  return 0
}

# Subset chromosome-level theta output back to one original BED file. The AWK
# block reads BED intervals first, then streams the chromosome theta table and
# keeps only rows whose positions fall inside those intervals.
split_region_set_output() {
  local bed_file="$1"
  local union_prefix="$2"
  local region_name region_out region_tmp chrom theta_pest header_written

  region_name=$(basename "${bed_file}" .bed)
  region_out="${OUTDIR}/${region_name}.pestPG"

  if [[ -f "${region_out}" ]] && split_output_is_current "${region_out}" "${bed_file}" "${union_prefix}"; then
    echo "Reusing split output: ${region_out}"
    return
  fi

  region_tmp=$(mktemp "${OUTDIR}/${region_name}.pestPG.tmp.XXXXXX")
  header_written=0

  while IFS= read -r chrom; do
    [[ -z "${chrom}" ]] && continue
    theta_pest="${union_prefix}_${chrom}.pestPG"
    if [[ ! -f "${theta_pest}" ]]; then
      continue
    fi

    awk -v chrom="${chrom}" -v print_header="${header_written}" '
      BEGIN {
        n = 0
        idx = 1
        pos_col = 0
      }
      FNR == NR {
        if ($1 != chrom) {
          next
        }
        n++
        start[n] = $2 + 1
        end[n] = $3
        next
      }
      FNR == 1 {
        if (print_header == 0) {
          print
        }
        for (i = 1; i <= NF; i++) {
          gsub(/^#/, "", $i)
          if ($i == "Pos" || $i == "Position" || $i == "midPos" || $i == "WinCenter" || $i == "pos") {
            pos_col = i
          }
        }
        if (pos_col == 0) {
          pos_col = 3
        }
        next
      }
      n == 0 {
        next
      }
      {
        pos = $(pos_col) + 0
        while (idx <= n && pos > end[idx]) {
          idx++
        }
        if (idx <= n && pos >= start[idx] && pos <= end[idx]) {
          print
        }
      }
    ' "${bed_file}" "${theta_pest}" >> "${region_tmp}"

    header_written=1
  done < <(chrom_stream "${bed_file}")

  mv "${region_tmp}" "${region_out}"
  echo "Wrote ${region_out}"
}

if [[ $# -ne 1 ]]; then
  usage
  exit 1
fi

BAMLIST=$(resolve_path "$1")
POP_NAME=$(basename "${BAMLIST}")
POP_NAME="${POP_NAME%.*}"
POP_NAME="${POP_NAME%_BAM}"
OUTDIR="${OUTBASE}/${POP_NAME}"

if [[ -z "${POP_NAME}" ]]; then
  echo "Failed to derive population name from BAM list: ${BAMLIST}" >&2
  exit 1
fi

require_exec "${ANGSD}"
require_exec "${REALSFS}"
require_exec "${THETASTAT}"

if [[ -z "${BEDTOOLS}" ]]; then
  echo "Required executable not found on PATH: bedtools" >&2
  exit 1
fi

require_file "${BAMLIST}"
require_file "${GENOME}"

mkdir -p "${PROJECT_ROOT}/logs/angsd" "${OUTDIR}"

shopt -s nullglob
region_files=( "${REGION_DIR}"/*.bed )
shopt -u nullglob

if [[ ${#region_files[@]} -eq 0 ]]; then
  echo "No BED region sets found in ${REGION_DIR}" >&2
  exit 1
fi

for bed_file in "${region_files[@]}"; do
  require_file "${bed_file}"
done

# Shared per-population filenames.
UNION_PREFIX="${OUTDIR}/${POP_NAME}.all_region_sets"
UNION_BED="${UNION_PREFIX}.merged.bed"

# Pipeline order:
# 1. build/reuse merged BED
# 2. build/reuse ANGSD + SFS + theta intermediates
# 3. build/reuse chromosome theta tables
# 4. split those chromosome tables into the original BED region sets
make_union_bed "${UNION_BED}" "${region_files[@]}"
if [[ ! -s "${UNION_BED}" ]]; then
  echo "Merged BED is empty: ${UNION_BED}" >&2
  exit 1
fi

index_sites_file "${UNION_BED}"
run_union_pipeline "${UNION_BED}" "${UNION_PREFIX}"
run_theta_by_chromosome "${UNION_BED}" "${UNION_PREFIX}"

for bed_file in "${region_files[@]}"; do
  split_region_set_output "${bed_file}" "${UNION_PREFIX}"
done

conda deactivate
