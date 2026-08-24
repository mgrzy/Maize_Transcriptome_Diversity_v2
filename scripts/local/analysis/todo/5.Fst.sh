### Script: analysis/5.Fst.sh
### Purpose: Compute Fst for SS and NSS groups within each era (exPVP/Public) using GCTA, starting from the Div639 markers (both all variants and the eQTL marker subset).
### Inputs:
###  - data/geno/plink/Div639 (PLINK bed/bim/fam)
###  - data/meta/Eqtl_All.txt (marker list)
###  - data/meta/GCTA_exPVP.txt, data/meta/GCTA_Public.txt (era-specific population definitions)
### Outputs:
###  - results/Fst/SS_vs_NSS_<era>_all_variants.* (per-era Fst using all variants)
###  - results/Fst/SS_vs_NSS_<era>_eqtl.* (per-era Fst using eQTL subset)
###  - helper population files for each era in results/Fst
### Usage: bash analysis/5.Fst.sh (requires env with plink and gcta)
### Prerequisites: conda env 'geno' with plink; GCTA at ../../soft/gcta/gcta
conda activate geno

RESULTS_DIR=results/Fst
mkdir -p "${RESULTS_DIR}"

# Run SS vs NSS Fst separately for each era-specific population file
BASE_BFILE=data/geno/plink/Div639
EQTL_LIST=data/meta/Eqtl_All.txt

declare -A ERA_POP_FILES=(
  ["exPVP"]=data/meta/GCTA_exPVP.txt
  ["Public"]=data/meta/GCTA_Public.txt
)

eras=("exPVP" "Public")

for era in "${eras[@]}"; do
  pop_file="${ERA_POP_FILES[${era}]}"
  if [[ ! -f "${pop_file}" ]]; then
    echo "Skipping era ${era}: population file ${pop_file} not found." >&2
    continue
  fi
  era_prefix="${RESULTS_DIR}/GCTA_${era}_SS_NSS"
  keep_file="${era_prefix}.keep"
  subpop_file="${era_prefix}.subpop"

  # Prepare keep and subpopulation files (mapping SS/NSS labels)
  awk 'BEGIN {OFS="\t"} {print $1, $2}' "${pop_file}" > "${keep_file}"
  awk 'BEGIN {OFS="\t"}
  {
    split($3, parts, "_");
    group = parts[1];
    if (group == "SS" || group == "NSS") {
      print $1, $2, group
    }
  }' "${pop_file}" > "${subpop_file}"

  if [[ ! -s "${keep_file}" || ! -s "${subpop_file}" ]]; then
    echo "Skipping era ${era}: no SS/NSS entries found in ${pop_file}" >&2
    continue
  fi

  era_bfile_all="data/geno/plink/AllVariants_${era}_SS_NSS"
  era_bfile_eqtl="data/geno/plink/Eqtl_All_${era}_SS_NSS"

  # Subset Div639 markers to SS/NSS lines for the given era (all variants)
  plink --bfile "${BASE_BFILE}" \
    --keep "${keep_file}" \
    --maf 0.05 \
    --make-bed \
    --out "${era_bfile_all}"

  # Run Fst comparing SS vs NSS within the era using all variants
  ../../soft/gcta/gcta --bfile "${era_bfile_all}" \
    --fst \
    --sub-popu "${subpop_file}" \
    --out "${RESULTS_DIR}/SS_vs_NSS_${era}_all_variants"

  # Subset Div639 markers to SS/NSS lines for the given era and eQTL markers
  plink --bfile "${BASE_BFILE}" \
    --extract "${EQTL_LIST}" \
    --keep "${keep_file}" \
    --make-bed \
    --out "${era_bfile_eqtl}"

  # Run Fst comparing SS vs NSS within the era using only eQTL markers
  ../../soft/gcta/gcta --bfile "${era_bfile_eqtl}" \
    --fst \
    --sub-popu "${subpop_file}" \
    --out "${RESULTS_DIR}/SS_vs_NSS_${era}_eqtl"
done
