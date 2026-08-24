#!/usr/bin/env bash
# Count cis variants tested for every gene in the full OSCA association output.
#
# OSCA writes one row per tested gene--variant pair in batch<N>.txt.gz.  This
# script streams the batches and counts those rows by Gene, so it is suitable
# for the full result set without loading it into memory.  Genes reported in
# OSCA *.failed.probe.list files are retained with zero tested variants.
#
# Usage:
#   bash scripts/local/procesing/16.Count_Cis_Tested_Variants_Per_Gene.sh
#   bash scripts/local/procesing/16.Count_Cis_Tested_Variants_Per_Gene.sh \
#     --input-dir results/full/osca \
#     --threads 10 \
#     --output results/summary/osca/cis_tested_variants_per_gene.tsv.gz
#
# Output columns:
#   gene                     OSCA Gene identifier
#   n_cis_tested_variants    Number of gene--variant association tests
#   osca_status              tested or failed_probe

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd -- "$script_dir/../../.." && pwd)"
input_dir="$project_root/results/full/osca"
output_file="$project_root/results/summary/osca/cis_tested_variants_per_gene.tsv.gz"
threads="${OSCA_COUNT_THREADS:-$(nproc)}"

usage() {
  cat <<'EOF'
Usage: 16.Count_Cis_Tested_Variants_Per_Gene.sh [options]

Options:
  --input-dir DIR  Directory containing OSCA batch<N>.txt.gz results.
  --threads N      Number of batches to process concurrently (default: CPU
                   cores, capped at the number of batches; can also be set
                   with OSCA_COUNT_THREADS).
  --output FILE    Output .tsv.gz file.
  -h, --help       Show this help message.
EOF
}

while (($#)); do
  case "$1" in
    --input-dir) input_dir="${2:?Missing directory after --input-dir}"; shift 2 ;;
    --threads) threads="${2:?Missing integer after --threads}"; shift 2 ;;
    --output) output_file="${2:?Missing file after --output}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -d "$input_dir" ]] || {
  printf 'Input directory does not exist: %s\n' "$input_dir" >&2
  exit 1
}
[[ "$output_file" == *.tsv.gz ]] || {
  printf 'Output file must end in .tsv.gz: %s\n' "$output_file" >&2
  exit 2
}
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || {
  printf '%s\n' "--threads must be a positive integer: $threads" >&2
  exit 2
}

mapfile -d '' batch_files < <(
  find "$input_dir" -maxdepth 1 -type f -regextype posix-extended \
    -regex '.*/batch[0-9]+\.txt\.gz' -print0 | sort -z -V
)
if ((${#batch_files[@]} == 0)); then
  printf 'No OSCA batch<N>.txt.gz files found in %s\n' "$input_dir" >&2
  exit 1
fi
if ((threads > ${#batch_files[@]})); then
  threads="${#batch_files[@]}"
fi

output_dir="$(dirname -- "$output_file")"
mkdir -p "$output_dir"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/osca-cis-test-counts.XXXXXX")"
trap 'rm -rf -- "$tmp_dir"' EXIT

tested_counts="$tmp_dir/tested_counts.tsv"
failed_probes="$tmp_dir/failed_probes.txt"
batch_count_dir="$tmp_dir/batch_counts"
mkdir -p "$batch_count_dir"

# Count every batch independently, allowing gzip decompression and parsing to
# run concurrently.  Only the small per-gene count tables are merged later.
count_batch() {
  local batch_file="$1"
  local batch_name
  local batch_output

  set -o pipefail
  batch_name="$(basename -- "$batch_file")"
  batch_output="$batch_count_dir/${batch_name%.txt.gz}.tsv"

  gzip -cd -- "$batch_file" | awk -F '\t' -v output_file="$batch_output" '
  function fail(message) {
    print "ERROR: " message > "/dev/stderr"
    exit 1
  }
  NR == 1 {
    header_seen = 1
    for (i = 1; i <= NF; i++) {
      if ($i == "Gene") gene_col = i
    }
    if (!gene_col) fail("OSCA header must contain a Gene column")
    next
  }
  {
    gene = $gene_col
    if (gene == "") fail("Encountered an association row without a Gene value")
    if (!(gene in n_tested)) n_genes++
    n_tested[gene]++
    total_tests++
  }
  END {
    if (!header_seen) fail("No readable OSCA header found")
    for (gene in n_tested) print gene "\t" n_tested[gene] > output_file
  }
'
}
export batch_count_dir
export -f count_batch

printf '%s\0' "${batch_files[@]}" |
  xargs -0 -r -n 1 -P "$threads" bash -c 'count_batch "$1"' _

awk -F '\t' '
  NF == 2 {
    n_tested[$1] += $2
    total_tests += $2
  }
  END {
    for (gene in n_tested) {
      print gene "\t" n_tested[gene]
      n_genes++
    }
    print "Counted " total_tests " cis association tests across " n_genes " tested genes." > "/dev/stderr"
  }
' "$batch_count_dir"/*.tsv > "$tested_counts"

# A failed probe has no association rows, but retaining it makes zero-test
# genes explicit in downstream joins and diagnostics.
find "$input_dir" -maxdepth 1 -type f -regextype posix-extended \
  -regex '.*/cis2Mb_[0-9]+_[0-9]+\.failed\.probe\.list' -print0 |
  sort -z -V |
  xargs -0r cat -- > "$failed_probes"

tmp_output="$tmp_dir/cis_tested_variants_per_gene.tsv.gz"
{
  printf 'gene\tn_cis_tested_variants\tosca_status\n'
  {
    awk -F '\t' 'NF == 2 { print $1 "\t" $2 "\ttested" }' "$tested_counts"
    awk 'NF { print $1 "\t0\tfailed_probe" }' "$failed_probes"
  } | LC_ALL=C sort -t $'\t' -k1,1 -k2,2nr |
    awk -F '\t' '!seen[$1]++'
} | gzip -c > "$tmp_output"
mv -f -- "$tmp_output" "$output_file"

n_tested_genes="$(wc -l < "$tested_counts")"
n_failed_probes="$(sort -u "$failed_probes" | wc -l)"
printf 'Processed %s OSCA batch file(s) using %s concurrent worker(s).\n' \
  "${#batch_files[@]}" "$threads"
printf 'Genes with cis tests: %s; failed probes retained with zero tests: %s.\n' \
  "$n_tested_genes" "$n_failed_probes"
printf 'Wrote %s\n' "$output_file"
