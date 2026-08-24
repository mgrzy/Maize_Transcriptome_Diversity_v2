#!/usr/bin/env bash
# Select one lowest-p-value eQTL variant per gene from compressed OSCA batches.
#
# Usage:
#   bash scripts/local/procesing/7.Select_OSCA_Min_P.sh
#   bash scripts/local/procesing/7.Select_OSCA_Min_P.sh --input-dir results/full/osca \
#     --output results/summary/eQTL_P_min.txt.gz
#
# The script streams the batches, so it never loads all association results into
# memory.  If equal minimum p-values occur for a gene, the first encountered
# variant is retained.

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd -- "$script_dir/../../.." && pwd)"
input_dir="$project_root/results/full/osca"
output_file="$project_root/results/summary/eQTL_P_min.txt.gz"

usage() {
  cat <<'EOF'
Usage: 7.Select_OSCA_Min_P.sh [options]

Options:
  --input-dir DIR  Directory containing batch<N>.txt.gz OSCA files.
  --output FILE    Compressed output file (must end in .gz).
  -h, --help       Show this help message.
EOF
}

while (($#)); do
  case "$1" in
    --input-dir) input_dir="${2:?Missing directory after --input-dir}"; shift 2 ;;
    --output) output_file="${2:?Missing file after --output}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -d "$input_dir" ]] || { printf 'Input directory does not exist: %s\n' "$input_dir" >&2; exit 1; }
[[ "$output_file" == *.gz ]] || { printf 'Output file must end in .gz: %s\n' "$output_file" >&2; exit 2; }

mapfile -d '' batch_files < <(find "$input_dir" -maxdepth 1 -type f -regextype posix-extended \
  -regex '.*/batch[0-9]+\.txt\.gz' -print0 | sort -z -V)
if ((${#batch_files[@]} == 0)); then
  printf 'No batch<N>.txt.gz files found in %s\n' "$input_dir" >&2
  exit 1
fi

output_dir="$(dirname -- "$output_file")"
mkdir -p "$output_dir"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/osca-min-p.XXXXXX")"
trap 'rm -rf -- "$tmp_dir"' EXIT
header_file="$tmp_dir/header.txt"
rows_file="$tmp_dir/rows.txt"

# OSCA association headers include Gene and p.  Their positions are detected
# from the header, keeping this valid if OSCA adds/reorders other columns.
for batch_file in "${batch_files[@]}"; do
  gzip -cd -- "$batch_file"
done | awk -F '\t' -v header_file="$header_file" -v rows_file="$rows_file" '
  function fail(message) { print "ERROR: " message > "/dev/stderr"; exit 1 }
  !header_seen {
    header_seen = 1
    header = $0
    for (i = 1; i <= NF; i++) {
      if ($i == "Gene") gene_col = i
      if ($i == "p") p_col = i
    }
    if (!gene_col || !p_col) fail("OSCA header must contain Gene and p columns")
    print header > header_file
    next
  }
  $0 == header { next }
  {
    gene = $gene_col
    p_text = $p_col
    if (gene == "" || p_text == "" || p_text == "NA") next
    p_value = p_text + 0
    if (!(gene in best_p) || p_value < best_p[gene]) {
      best_p[gene] = p_value
      best_row[gene] = $0
    }
  }
  END {
    if (!header_seen) fail("No readable OSCA header found")
    for (gene in best_row) print best_row[gene] > rows_file
  }
'

[[ -s "$header_file" && -s "$rows_file" ]] || { printf 'No valid eQTL rows were found.\n' >&2; exit 1; }
tmp_output="$tmp_dir/eQTL_P_min.txt.gz"
{
  cat "$header_file"
  # Gene is column 10 in the standard OSCA association output.
  LC_ALL=C sort -t $'\t' -k10,10 "$rows_file"
} | gzip -c > "$tmp_output"
mv -f -- "$tmp_output" "$output_file"

printf 'Selected %s genes from %s OSCA batch file(s).\n' "$(wc -l < "$rows_file")" "${#batch_files[@]}"
printf 'Wrote %s\n' "$output_file"
