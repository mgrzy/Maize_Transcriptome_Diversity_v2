#!/usr/bin/env bash
# Purpose: Extract ZeroShotScore annotations for variants present in the
#          unfiltered Div639 PLINK dataset.
# Inputs:
#   - data/processed/geno/plink/Div639.bim
#   - data/processed/Cad/chr*_schnable_scores.tsv.gz
# Output:
#   - results/summary/Cad/Div639_zeroShotScore.tsv.gz
# Usage: bash scripts/local/analysis/11.Extract_Div639_ZeroShotScore.sh

set -euo pipefail

bim_file="data/processed/geno/plink/Div639.bim"
cad_dir="data/processed/Cad"
output_dir="results/summary/Cad"
output_file="${output_dir}/Div639_zeroShotScore.tsv.gz"

[[ -s "$bim_file" ]] || {
  echo "Missing Div639 BIM file: $bim_file" >&2
  exit 1
}
[[ -d "$cad_dir" ]] || {
  echo "Missing CAD directory: $cad_dir" >&2
  exit 1
}

mkdir -p "$output_dir"
task_tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/div639_zeroshot.XXXXXX")
trap 'rm -rf -- "$task_tmp_dir"' EXIT

# Partition the Div639 identifiers by chromosome.  The score files identify a
# variant as chr<chromosome>_<position>, so retaining the BIM ID guarantees
# that the output IDs exactly match the PLINK data.
awk -v tmp_dir="$task_tmp_dir" '
  {
    chromosome = $1
    sub(/^chr/, "", chromosome)
    print $4 "\t" $2 > (tmp_dir "/chr" chromosome ".ids")
  }
' "$bim_file"

plain_output="$task_tmp_dir/Div639_zeroShotScore.tsv"
{
  printf 'variant_id\tchr\tpos\tref\talt\tzeroShotScore\n'

  for chromosome in {1..10}; do
    id_file="$task_tmp_dir/chr${chromosome}.ids"
    cad_file="$cad_dir/chr${chromosome}_schnable_scores.tsv.gz"

    [[ -s "$id_file" ]] || continue
    [[ -s "$cad_file" ]] || {
      echo "Missing CAD score file: $cad_file" >&2
      exit 1
    }

    # The CAD file has columns: chr, start, end, pos, ref, alt, sequences,
    # zeroShotScore.  Do not retain the large sequence field in the result.
    gzip -cd "$cad_file" |
      awk -F '\t' -v OFS='\t' -v ids="$id_file" '
        FILENAME == ids {
          variant_id[$1] = $2
          next
        }
        FNR == 1 { next }
        $4 in variant_id {
          print variant_id[$4], $1, $4, $5, $6, $8
        }
      ' "$id_file" -

    echo "Processed chromosome ${chromosome}" >&2
  done
} > "$plain_output"

gzip -n -c "$plain_output" > "${output_file}.tmp"
mv "${output_file}.tmp" "$output_file"

printf 'Wrote %s variants to %s\n' "$(awk 'END { print NR - 1 }' "$plain_output")" "$output_file"
