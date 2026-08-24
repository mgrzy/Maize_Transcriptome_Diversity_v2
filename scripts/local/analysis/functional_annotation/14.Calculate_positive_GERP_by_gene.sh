#!/usr/bin/env bash
# Purpose: Calculate mean positive GERP scores for every annotated gene body
#          and for its strand-aware promoter (-2,000 bp through the TSS).
#
# The GERP BED uses 0-based, half-open coordinates. Gene and promoter
# coordinates are derived from the v5 GFF and converted internally to the
# same convention. Only values > 0 in the final GERP BED column are included.
#
# Inputs:
#   - data/processed/gerp/gerp_genes.bed.gz
#   - ../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
# Output:
#   - results/summary/GERP/positive_GERP_mean_by_gene_and_promoter.csv
#
# Usage:
#   bash scripts/local/analysis/14.Calculate_positive_GERP_by_gene.sh

set -euo pipefail

gerp_file="${GERP_FILE:-data/processed/gerp/gerp_genes.bed.gz}"
gff_file="${GFF_PATH:-../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz}"
output_file="${OUTPUT_FILE:-results/summary/GERP/positive_GERP_mean_by_gene_and_promoter.csv}"

[[ -s "$gerp_file" ]] || {
  echo "Missing GERP BED file: $gerp_file" >&2
  exit 1
}
[[ -s "$gff_file" ]] || {
  echo "Missing GFF file: $gff_file" >&2
  exit 1
}

mkdir -p "$(dirname "$output_file")"
task_tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/positive_gerp_by_gene.XXXXXX")
trap 'rm -rf -- "$task_tmp_dir"' EXIT

features_file="$task_tmp_dir/gene_and_promoter_features.bed"
long_results_file="$task_tmp_dir/positive_gerp_feature_means.csv"

# Build gene-body and promoter features. Promoters are strand-aware:
#   + strand: [TSS - 2,000, TSS] inclusive
#   - strand: [TSS, TSS + 2,000] inclusive
# In BED coordinates this means a 2,001-bp interval that includes the TSS.
gzip -cd "$gff_file" |
  awk -F '\t' -v OFS='\t' '
    $3 == "gene" {
      gene_id = ""
      attributes_n = split($9, attributes, ";")
      for (attribute_index = 1; attribute_index <= attributes_n; attribute_index++) {
        if (attributes[attribute_index] ~ /^ID=/) {
          gene_id = attributes[attribute_index]
          sub(/^ID=/, "", gene_id)
          break
        }
      }
      if (gene_id == "") next

      gene_start = $4 - 1  # convert GFF 1-based start to BED 0-based start
      gene_end = $5         # BED end is exclusive, equal to GFF inclusive end
      strand = $7
      tss = (strand == "-") ? gene_end - 1 : gene_start

      print $1, gene_start, gene_end, gene_id, "gene", strand, gene_start, gene_end, tss

      if (strand == "-") {
        promoter_start = tss
        promoter_end = tss + 2001
      } else {
        promoter_start = tss - 2000
        if (promoter_start < 0) promoter_start = 0
        promoter_end = tss + 1
      }
      print $1, promoter_start, promoter_end, gene_id, "promoter", strand, gene_start, gene_end, tss
    }
  ' |
  LC_ALL=C sort -k1,1V -k2,2n > "$features_file"

[[ -s "$features_file" ]] || {
  echo "No gene features were read from: $gff_file" >&2
  exit 1
}

# Sweep the sorted, positive-score GERP sites across the sorted feature list.
# Only currently overlapping gene/promoter intervals are kept active, avoiding
# loading the 163-million-site GERP file into memory.
gzip -cd "$gerp_file" |
  awk -F '\t' -v OFS=',' -v feature_file="$features_file" '
    FILENAME == feature_file {
      chromosome = $1
      feature_index[chromosome]++
      key = chromosome SUBSEP feature_index[chromosome]
      feature_chromosome[key] = chromosome
      feature_start[key] = $2
      feature_end[key] = $3
      feature_gene[key] = $4
      feature_type[key] = $5
      feature_strand[key] = $6
      feature_gene_start[key] = $7
      feature_gene_end[key] = $8
      feature_tss[key] = $9
      feature_order[++n_features] = key
      next
    }

    {
      score = $NF + 0
      if (score <= 0) next

      chromosome = $1
      position = $2
      if (chromosome != active_chromosome) {
        for (key in active) delete active[key]
        active_chromosome = chromosome
      }

      next_index = next_feature_index[chromosome] + 0
      while (next_index < feature_index[chromosome]) {
        candidate_key = chromosome SUBSEP (next_index + 1)
        if (feature_start[candidate_key] > position) break
        active[candidate_key] = 1
        next_index++
      }
      next_feature_index[chromosome] = next_index

      for (key in active) {
        if (feature_end[key] <= position) {
          delete active[key]
        } else {
          score_sum[key] += score
          positive_site_count[key]++
        }
      }
    }

    END {
      print "chr", "start_0based", "end_0based", "gene_id", "feature", "strand", "gene_start_0based", "gene_end_0based", "tss_0based", "n_positive_GERP_sites", "mean_positive_GERP"
      for (feature_number = 1; feature_number <= n_features; feature_number++) {
        key = feature_order[feature_number]
        n = positive_site_count[key] + 0
        mean = (n > 0) ? score_sum[key] / n : "NA"
        print feature_chromosome[key], feature_start[key], feature_end[key], feature_gene[key], feature_type[key], feature_strand[key], feature_gene_start[key], feature_gene_end[key], feature_tss[key], n, mean
      }
    }
  ' "$features_file" - > "$long_results_file"

# Convert the two feature rows per gene into one analysis-ready row.
awk -F ',' -v OFS=',' '
  NR == 1 { next }
  {
    gene_id = $4
    feature = $5
    if (!(gene_id in seen)) {
      seen[gene_id] = 1
      gene_order[++n_genes] = gene_id
      chromosome[gene_id] = $1
      strand[gene_id] = $6
      gene_start[gene_id] = $7 + 1
      gene_end[gene_id] = $8
      tss[gene_id] = $9 + 1
    }
    if (feature == "gene") {
      gene_n[gene_id] = $10
      gene_mean[gene_id] = $11
    } else if (feature == "promoter") {
      promoter_start[gene_id] = $2 + 1
      promoter_end[gene_id] = $3
      promoter_n[gene_id] = $10
      promoter_mean[gene_id] = $11
    }
  }
  END {
    print "gene_id", "chr", "strand", "gene_start_1based", "gene_end_1based", "tss_1based", "promoter_start_1based", "promoter_end_1based", "gene_n_positive_GERP_sites", "gene_mean_positive_GERP", "promoter_n_positive_GERP_sites", "promoter_mean_positive_GERP"
    for (gene_number = 1; gene_number <= n_genes; gene_number++) {
      gene_id = gene_order[gene_number]
      print gene_id, chromosome[gene_id], strand[gene_id], gene_start[gene_id], gene_end[gene_id], tss[gene_id], promoter_start[gene_id], promoter_end[gene_id], gene_n[gene_id] + 0, gene_mean[gene_id], promoter_n[gene_id] + 0, promoter_mean[gene_id]
    }
  }
' "$long_results_file" > "$output_file"

printf 'Wrote positive-GERP means for %s genes to %s\n' "$(awk -F ',' 'NR > 1 { count++ } END { print count + 0 }' "$output_file")" "$output_file"
