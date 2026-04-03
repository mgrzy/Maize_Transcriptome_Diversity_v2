#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Subset an ANGSD thetaStat .pestPG file by one or more BED files.

Usage:
  bash scripts/local/analysis/subset_pestPG_by_bed.sh \
    --pestpg results/full/pop/angsd/all_region_sets_chr9.pestPG \
    --bed_dir results/summary/bed_windows \
    --out_dir results/full/pop/angsd

Options:
  --pestpg FILE    Input .pestPG file from thetaStat.
  --bed FILE       BED file to subset against. Can be given multiple times.
  --bed_dir DIR    Directory of BED files. All *.bed files will be used.
  --out_dir DIR    Output directory. Defaults to the pestPG directory.
  --force          Overwrite outputs even if they are up to date.
  --help           Show this help.

Notes:
  - BED is assumed to be 0-based, half-open.
  - .pestPG is assumed to have chromosome in column 2 and position in column 3.
  - One output file per BED is written as <bed_basename>.pestPG.
EOF
}

require_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
}

subset_one() {
  local bed_file="$1"
  local pestpg_file="$2"
  local out_file="$3"

  awk '
    BEGIN {
      FS = OFS = "\t"
    }
    FNR == NR {
      if ($0 == "" || $1 ~ /^#/) {
        next
      }
      chrom = $1
      n[chrom]++
      start[chrom, n[chrom]] = $2 + 1
      end[chrom, n[chrom]] = $3
      next
    }
    FNR == 1 {
      print
      next
    }
    {
      chrom = $2
      pos = $3 + 0
      if (!(chrom in n)) {
        next
      }
      if (!(chrom in idx)) {
        idx[chrom] = 1
      }
      while (idx[chrom] <= n[chrom] && pos > end[chrom, idx[chrom]]) {
        idx[chrom]++
      }
      if (idx[chrom] <= n[chrom] && pos >= start[chrom, idx[chrom]] && pos <= end[chrom, idx[chrom]]) {
        print
      }
    }
  ' "${bed_file}" "${pestpg_file}" > "${out_file}"
}

pestpg_file=""
out_dir=""
force=0
declare -a bed_files=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pestpg)
      pestpg_file="$2"
      shift 2
      ;;
    --bed)
      bed_files+=( "$2" )
      shift 2
      ;;
    --bed_dir)
      while IFS= read -r bed_file; do
        bed_files+=( "${bed_file}" )
      done < <(find "$2" -maxdepth 1 -type f -name '*.bed' | sort)
      shift 2
      ;;
    --out_dir)
      out_dir="$2"
      shift 2
      ;;
    --force)
      force=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "${pestpg_file}" ]]; then
  echo "--pestpg is required" >&2
  usage >&2
  exit 1
fi

if [[ ${#bed_files[@]} -eq 0 ]]; then
  echo "Provide at least one --bed or --bed_dir" >&2
  usage >&2
  exit 1
fi

require_file "${pestpg_file}"

if [[ -z "${out_dir}" ]]; then
  out_dir=$(dirname "${pestpg_file}")
fi
mkdir -p "${out_dir}"

for bed_file in "${bed_files[@]}"; do
  require_file "${bed_file}"
done

for bed_file in "${bed_files[@]}"; do
  region_name=$(basename "${bed_file}" .bed)
  out_file="${out_dir}/${region_name}.pestPG"

  if [[ ${force} -eq 0 && -f "${out_file}" && "${out_file}" -nt "${bed_file}" && "${out_file}" -nt "${pestpg_file}" ]]; then
    echo "Reusing ${out_file}"
    continue
  fi

  echo "Subsetting $(basename "${pestpg_file}") with $(basename "${bed_file}")"
  subset_one "${bed_file}" "${pestpg_file}" "${out_file}"
  echo "Wrote ${out_file}"
done
