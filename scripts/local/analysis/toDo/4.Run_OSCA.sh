### Script: analysis/4.Run_OSCA.sh
### Purpose: Run OSCA eQTL mapping with genotype PLINK files, binary expression, GRM, and PC covariates; cis window 2kb.
### Inputs:
###  - data/geno/plink/Div639 (PLINK bed/bim/fam)
###  - data/OSCA/WiDiv631 (binary expression befile)
###  - data/OSCA/WiDiv631_Kin (GRM)
###  - data/OSCA/WiDiv631_PC.txt (qcovar)
### Outputs:
###  - results/eQTL/test.*
### Usage: bash analysis/4.Run_OSCA.sh (from repo root; adjust paths to osca binary)
### Prerequisites: OSCA binary at ../../soft/osca/osca
### Run OSCA
../../soft/osca/osca --eqtl --bfile data/geno/plink/Div639 \
--befile data/OSCA/WiDiv631 \
--grm data/OSCA/WiDiv631_Kin \
--qcovar data/OSCA/WiDiv631_PC.txt \
--mlm --task-num 1 --task-id 1 --thread-num 8 --cis --cis-wind 2000 --no-fid --out results/eQTL/test
