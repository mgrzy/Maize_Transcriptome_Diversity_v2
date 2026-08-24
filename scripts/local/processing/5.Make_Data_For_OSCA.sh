### Make BOD file
../../soft/osca/osca --efile data/processed/osca/WiDiv631.dat --gene-expression \
--make-bod --out data/processed/osca/WiDiv631 --no-fid

../../soft/osca/osca --befile data/processed/osca/WiDiv631 --update-opi data/processed/osca/WiDiv631_OPI.txt
### Script: processing/6.Make_Data_For_OSCA.sh
### Purpose: Generate OSCA-compatible input files (binary expression, GRM, and covariates) from processed datasets.
### Inputs: PLINK genotypes, expression matrices, covariate files
### Outputs: data/OSCA/* (befile/GRM/PC/covariate files)
### Usage: bash processing/6.Make_Data_For_OSCA.sh
### Prerequisites: OSCA tools; PLINK; required data prepared
