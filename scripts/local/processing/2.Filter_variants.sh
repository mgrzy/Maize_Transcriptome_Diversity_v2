### Script: processing/2.Filter_variants.sh
### Purpose: Filter variants by MAF and sample list using PLINK; produce PLINK dataset for selected lines.
### Inputs:
###  - ../BigData/Maize/v5/Genotypes/All/plink/All_maf001 (PLINK bed/bim/fam)
###  - data/meta/LinesPlink.txt (keep list)
### Outputs:
###  - data/geno/plink/Div639 (filtered PLINK dataset)
### Usage: bash processing/2.Filter_variants.sh
### Prerequisites: conda env 'geno' with plink installed
conda activate geno

plink --bfile ../BigData/Maize/v5/Genotypes/All/plink/All_maf001 \
--maf 0.01 --keep data/meta/LinesPlink.txt \
--make-bed --out data/geno/plink/Div639

