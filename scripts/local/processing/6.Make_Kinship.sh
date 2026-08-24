../../soft/gcta/gcta --bfile data/processed/geno/plink/Div639 --out data/processed/osca/WiDiv631_Kin --thread-num 24 --maf 0.1 --make-grm 
### Script: processing/7.Make_Kinship.sh
### Purpose: Compute kinship matrix (GRM) from genotype data for downstream mixed models.
### Inputs: PLINK bed/bim/fam for target samples
### Outputs: data/OSCA/*_Kin or similar GRM files
### Usage: bash processing/7.Make_Kinship.sh
### Prerequisites: Tool that computes GRM (e.g., OSCA, GCTA) installed and reachable
