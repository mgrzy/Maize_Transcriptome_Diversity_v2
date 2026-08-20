#!/usr/bin/env bash
# Build replicate-supported consensus TF-binding BEDs from the Tu/Mejia-Guerra
# 2020 MaizeGDB B73 v5 PeakRanger region files. These are consensus peaks based
# on replicate support; they are not formal IDR peaks.
#
# Usage:
#   bash scripts/local/prepare_Tu_TF_consensus_peaks.sh \
#     --metadata path/to/Tu_Mejia-Guerra_2020_read_data.txt \
#     --peaks-dir path/to/maizegdb_region_peaks \
#     --out-dir data/processed/chipseq/Tu_2020_TF_consensus_v5 \
#     --support-mode at-least-two
#
# Options:
#   --metadata FILE
#   --peaks-dir DIR
#   --out-dir DIR
#   --support-mode at-least-two|majority|all   (default: at-least-two)
#   --force                                    replace prior generated outputs
#   --dry-run                                  inspect metadata and inputs only
#   --tf NAME                                  process one TF only (testing aid)
#   --help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

METADATA="${PROJECT_ROOT}/data/processed/chromatin/Tu_Mejia-Guerra_TF_Peaks/Tu_Mejia-Guerra_2020_read_data.txt"
PEAKS_DIR="${PROJECT_ROOT}/data/processed/chromatin/Tu_Mejia-Guerra_TF_Peaks"
OUT_DIR="${PROJECT_ROOT}/data/processed/chipseq/Tu_2020_TF_consensus_v5"
SUPPORT_MODE="at-least-two"
FORCE=0
DRY_RUN=0
TF_FILTER=""
TMP_DIR=""

usage() {
  sed -n '1,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

info() {
  printf '[%(%F %T)T] %s\n' -1 "$*"
}

cleanup() {
  if [[ -n "${TMP_DIR}" && -d "${TMP_DIR}" ]]; then
    rm -rf -- "${TMP_DIR}"
  fi
}
trap cleanup EXIT

while (( $# )); do
  case "$1" in
    --metadata)
      (( $# >= 2 )) || die "--metadata requires a file path"
      METADATA="$2"
      shift 2
      ;;
    --peaks-dir)
      (( $# >= 2 )) || die "--peaks-dir requires a directory path"
      PEAKS_DIR="$2"
      shift 2
      ;;
    --out-dir)
      (( $# >= 2 )) || die "--out-dir requires a directory path"
      OUT_DIR="$2"
      shift 2
      ;;
    --support-mode)
      (( $# >= 2 )) || die "--support-mode requires a value"
      SUPPORT_MODE="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --tf)
      (( $# >= 2 )) || die "--tf requires a TF name"
      TF_FILTER="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "Unknown option: $1"
      ;;
  esac
done

case "${SUPPORT_MODE}" in
  at-least-two|majority|all) ;;
  *) die "--support-mode must be at-least-two, majority, or all" ;;
esac

[[ -f "${METADATA}" ]] || die "Metadata file not found: ${METADATA}"
[[ -d "${PEAKS_DIR}" ]] || die "Peaks directory not found: ${PEAKS_DIR}"

required_commands=(awk sort gzip wc bedtools)
missing_commands=()
for command_name in "${required_commands[@]}"; do
  command -v "${command_name}" >/dev/null 2>&1 || missing_commands+=("${command_name}")
done
if (( ${#missing_commands[@]} )); then
  if (( DRY_RUN )); then
    info "DRY RUN: missing required workflow program(s): ${missing_commands[*]}"
  else
    die "Required program(s) not found in PATH: ${missing_commands[*]}"
  fi
fi

parse_metadata() {
  awk -F '\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "run_accession") run_column = i
        if ($i == "sample_title") title_column = i
      }
      if (!run_column || !title_column) {
        print "Metadata header must contain run_accession and sample_title" > "/dev/stderr"
        exit 2
      }
      next
    }
    {
      run = $run_column
      title = $title_column
      if (run == "" || title !~ /_rep[0-9]+$/) next
      tf = title
      sub(/_rep[0-9]+$/, "", tf)
      if (tf !~ /_TF/) next
      print tf "\t" run "\t" title
    }
  ' "${METADATA}"
}

resolve_peak_file() {
  local run_accession="$1"
  if [[ -f "${PEAKS_DIR}/${run_accession}_region.bed" ]]; then
    printf '%s\n' "${PEAKS_DIR}/${run_accession}_region.bed"
  elif [[ -f "${PEAKS_DIR}/${run_accession}_region.bed.gz" ]]; then
    printf '%s\n' "${PEAKS_DIR}/${run_accession}_region.bed.gz"
  else
    return 1
  fi
}

minimum_support() {
  local n_replicates="$1"
  case "${SUPPORT_MODE}" in
    at-least-two)
      printf '2\n'
      ;;
    majority)
      printf '%d\n' $(( n_replicates / 2 + 1 ))
      ;;
    all)
      printf '%d\n' "${n_replicates}"
      ;;
  esac
}

bed_sort() {
  LC_ALL=C sort -k1,1V -k2,2n -k3,3n
}

# Print valid chr1--chr10 BED3 rows. Malformed non-comment/nonblank lines fail
# the run rather than being allowed to become a header or invalid final record.
read_and_validate_bed3() {
  local input_file="$1"
  if [[ "${input_file}" == *.gz ]]; then
    gzip -cd -- "${input_file}"
  else
    cat -- "${input_file}"
  fi | awk -v source_file="${input_file}" '
    BEGIN { OFS = "\t"; invalid = 0 }
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    {
      if (NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || ($3 + 0) <= ($2 + 0)) {
        print "Invalid BED row in " source_file " at input line " NR ": " $0 > "/dev/stderr"
        invalid = 1
        next
      }
      if ($1 ~ /^chr([1-9]|10)$/) print $1, $2, $3
    }
    END { if (invalid) exit 2 }
  '
}

declare -A TF_RUNS=()
declare -A TF_TITLES=()
declare -A TF_METADATA_REPLICATES=()
declare -a TF_NAMES=()

metadata_records="$(parse_metadata)" || die "Could not parse TF metadata from: ${METADATA}"
[[ -n "${metadata_records}" ]] || die "No TF replicate samples were found in metadata: ${METADATA}"

while IFS=$'\t' read -r tf_name run_accession sample_title; do
  [[ -n "${tf_name}" && -n "${run_accession}" ]] || continue
  if [[ -z "${TF_METADATA_REPLICATES[${tf_name}]+x}" ]]; then
    TF_NAMES+=("${tf_name}")
    TF_METADATA_REPLICATES["${tf_name}"]=0
    TF_RUNS["${tf_name}"]=""
    TF_TITLES["${tf_name}"]=""
  fi
  TF_METADATA_REPLICATES["${tf_name}"]=$(( TF_METADATA_REPLICATES[${tf_name}] + 1 ))
  TF_RUNS["${tf_name}"]+="${TF_RUNS[${tf_name}]:+$'\n'}${run_accession}"
  TF_TITLES["${tf_name}"]+="${TF_TITLES[${tf_name}]:+$'\n'}${sample_title}"
done <<< "${metadata_records}"

mapfile -t TF_NAMES < <(printf '%s\n' "${TF_NAMES[@]}" | LC_ALL=C sort)

if [[ -n "${TF_FILTER}" ]]; then
  [[ -n "${TF_METADATA_REPLICATES[${TF_FILTER}]+x}" ]] || die "TF not found in parsed metadata: ${TF_FILTER}"
  TF_NAMES=("${TF_FILTER}")
fi

count_two=0
count_three=0
count_four=0
count_other=0
for tf_name in "${TF_NAMES[@]}"; do
  case "${TF_METADATA_REPLICATES[${tf_name}]}" in
    2) ((count_two += 1)) ;;
    3) ((count_three += 1)) ;;
    4) ((count_four += 1)) ;;
    *) ((count_other += 1)) ;;
  esac
done

info "Parsed ${#TF_NAMES[@]} TF(s): ${count_two} with 2 replicates, ${count_three} with 3, ${count_four} with 4, ${count_other} other."
if [[ -z "${TF_FILTER}" && ( ${#TF_NAMES[@]} -ne 104 || ${count_two} -ne 97 || ${count_three} -ne 5 || ${count_four} -ne 2 || ${count_other} -ne 0 ) ]]; then
  printf 'WARNING: expected 104 TFs (97 x 2 replicates, 5 x 3, 2 x 4), but observed %d TFs (%d x 2, %d x 3, %d x 4, %d other). No groups were dropped.\n' \
    "${#TF_NAMES[@]}" "${count_two}" "${count_three}" "${count_four}" "${count_other}" >&2
fi

missing_count=0
for tf_name in "${TF_NAMES[@]}"; do
  while IFS= read -r run_accession; do
    if ! resolve_peak_file "${run_accession}" >/dev/null; then
      ((missing_count += 1))
    fi
  done <<< "${TF_RUNS[${tf_name}]}"
done
info "Found ${missing_count} missing replicate peak file(s) among selected TF metadata rows."

if (( DRY_RUN )); then
  info "DRY RUN complete. No output directory or BED files were created."
  if (( ${#missing_commands[@]} )); then
    info "Install the missing program(s) before a real run: ${missing_commands[*]}"
  fi
  exit 0
fi

if [[ -d "${OUT_DIR}" ]] && [[ -n "$(find "${OUT_DIR}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  (( FORCE )) || die "Output directory is nonempty: ${OUT_DIR}. Use --force to replace generated outputs."
  info "--force supplied: removing only prior generated consensus outputs under ${OUT_DIR}."
  rm -rf -- "${OUT_DIR}/per_tf" "${OUT_DIR}/per_tf_named" "${OUT_DIR}/cleaned"
  rm -f -- "${OUT_DIR}/all_TF_consensus_named.bed" \
    "${OUT_DIR}/any_TF_consensus_union.bed" \
    "${OUT_DIR}/TF_consensus_QC.tsv" \
    "${OUT_DIR}/missing_inputs.tsv" \
    "${OUT_DIR}/run_to_TF_mapping.tsv"
fi

mkdir -p "${OUT_DIR}/per_tf" "${OUT_DIR}/per_tf_named" "${OUT_DIR}/cleaned"
TMP_DIR="$(mktemp -d "${OUT_DIR}/.prepare_Tu_TF_consensus.XXXXXX")"

qc_file="${OUT_DIR}/TF_consensus_QC.tsv"
missing_file="${OUT_DIR}/missing_inputs.tsv"
mapping_file="${OUT_DIR}/run_to_TF_mapping.tsv"
printf 'TF\tn_replicates\tn_metadata_replicates\trun_accessions\tsupport_mode\tminimum_replicate_support\tinput_peak_counts_by_replicate\tconsensus_peak_count\tconsensus_covered_bp\tfraction_of_each_replicate_peaks_overlapping_consensus\tstatus\n' > "${qc_file}"
printf 'TF\trun_accession\tsample_title\texpected_peak_file\treason\n' > "${missing_file}"
printf 'TF\trun_accession\tsample_title\tpeak_file\tstatus\n' > "${mapping_file}"

for tf_name in "${TF_NAMES[@]}"; do
  info "Processing ${tf_name}"
  tf_clean_dir="${OUT_DIR}/cleaned/${tf_name}"
  mkdir -p "${tf_clean_dir}"

  mapfile -t runs < <(printf '%s\n' "${TF_RUNS[${tf_name}]}")
  mapfile -t titles < <(printf '%s\n' "${TF_TITLES[${tf_name}]}")
  available_files=()
  available_runs=()
  input_counts=()
  fractions=()

  for run_index in "${!runs[@]}"; do
    run_accession="${runs[${run_index}]}"
    sample_title="${titles[${run_index}]}"
    if peak_file="$(resolve_peak_file "${run_accession}")"; then
      printf '%s\t%s\t%s\t%s\tfound\n' "${tf_name}" "${run_accession}" "${sample_title}" "${peak_file}" >> "${mapping_file}"
      raw_clean="${TMP_DIR}/${tf_name}.${run_accession}.valid.bed"
      clean_file="${tf_clean_dir}/${run_accession}.bed"
      read_and_validate_bed3 "${peak_file}" | bed_sort > "${raw_clean}"
      input_count="$(wc -l < "${raw_clean}")"
      bedtools merge -i "${raw_clean}" > "${clean_file}"
      available_files+=("${clean_file}")
      available_runs+=("${run_accession}")
      input_counts+=("${run_accession}=${input_count}")
    else
      expected_file="${PEAKS_DIR}/${run_accession}_region.bed[.gz]"
      printf '%s\t%s\t%s\t%s\tmissing_region_bed\n' "${tf_name}" "${run_accession}" "${sample_title}" "${expected_file}" >> "${mapping_file}"
      printf '%s\t%s\t%s\t%s\tmissing_region_bed\n' "${tf_name}" "${run_accession}" "${sample_title}" "${expected_file}" >> "${missing_file}"
    fi
  done

  n_available="${#available_files[@]}"
  min_support="$(minimum_support "${n_available}")"
  run_list="$(IFS=';'; printf '%s' "${runs[*]}")"
  input_count_list="$(IFS=';'; printf '%s' "${input_counts[*]:-}")"
  consensus_file="${OUT_DIR}/per_tf/${tf_name}_consensus.bed"
  named_file="${OUT_DIR}/per_tf_named/${tf_name}_consensus_named.bed"
  status="complete"

  if (( n_available < 2 )); then
    status="skipped_fewer_than_two_available_replicates"
    : > "${consensus_file}"
    : > "${named_file}"
    printf '%s\t%d\t%d\t%s\t%s\t%s\t%s\t0\t0\t\t%s\n' \
      "${tf_name}" "${n_available}" "${TF_METADATA_REPLICATES[${tf_name}]}" "${run_list}" \
      "${SUPPORT_MODE}" "${min_support}" "${input_count_list}" "${status}" >> "${qc_file}"
    info "  ${tf_name}: ${status}"
    continue
  fi

  if (( n_available < TF_METADATA_REPLICATES[${tf_name}] )); then
    status="complete_with_missing_replicates"
  fi

  multiinter_file="${TMP_DIR}/${tf_name}.multiinter.bed"
  bedtools multiinter -i "${available_files[@]}" > "${multiinter_file}"
  awk -v minimum_support="${min_support}" 'BEGIN { OFS = "\t" } $4 >= minimum_support { print $1, $2, $3 }' "${multiinter_file}" \
    | bed_sort \
    | bedtools merge -i - > "${consensus_file}"
  awk -v tf_name="${tf_name}" 'BEGIN { OFS = "\t" } { print $1, $2, $3, tf_name }' "${consensus_file}" > "${named_file}"

  consensus_count="$(wc -l < "${consensus_file}")"
  consensus_bp="$(awk '{ total += $3 - $2 } END { print total + 0 }' "${consensus_file}")"
  if (( consensus_count == 0 )); then
    status="${status}_no_consensus_segments"
  fi

  for replicate_index in "${!available_files[@]}"; do
    replicate_file="${available_files[${replicate_index}]}"
    replicate_count="$(wc -l < "${replicate_file}")"
    overlap_count=0
    if (( replicate_count > 0 && consensus_count > 0 )); then
      overlap_count="$(bedtools intersect -u -a "${replicate_file}" -b "${consensus_file}" | wc -l)"
    fi
    fraction="$(awk -v numerator="${overlap_count}" -v denominator="${replicate_count}" 'BEGIN { if (denominator) printf "%.6f", numerator / denominator; else printf "NA" }')"
    fractions+=("${available_runs[${replicate_index}]}=${fraction}")
  done
  fraction_list="$(IFS=';'; printf '%s' "${fractions[*]:-}")"

  printf '%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tf_name}" "${n_available}" "${TF_METADATA_REPLICATES[${tf_name}]}" "${run_list}" \
    "${SUPPORT_MODE}" "${min_support}" "${input_count_list}" "${consensus_count}" "${consensus_bp}" \
    "${fraction_list}" "${status}" >> "${qc_file}"
done

info "Writing combined TF-specific named BED. TF identity is retained; intervals are not merged across TFs."
find "${OUT_DIR}/per_tf_named" -type f -name '*_consensus_named.bed' -size +0c -print0 \
  | xargs -0 -r cat \
  | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -k4,4 \
  > "${OUT_DIR}/all_TF_consensus_named.bed"

info "Writing union of reproducible TF-binding regions across all TFs."
find "${OUT_DIR}/per_tf" -type f -name '*_consensus.bed' -size +0c -print0 \
  | xargs -0 -r cat \
  | bed_sort \
  | bedtools merge -i - \
  > "${OUT_DIR}/any_TF_consensus_union.bed"

info "Completed consensus construction. QC: ${qc_file}"
