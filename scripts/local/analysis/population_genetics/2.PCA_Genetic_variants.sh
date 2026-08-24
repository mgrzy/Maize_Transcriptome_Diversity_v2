### Calulate PCA based on genetic markers 

###
conda activate geno

###
plink --bfile data/processed/geno/plink/Div639 \
--maf 0.05 --pca 5 --out results/full/PCA/PCA_Geno
### Script: analysis/2.PCA_Genetic_variants.sh
### Purpose: Compute PCA for genetic variants using PLINK or related tools; prepares eigenvectors for downstream plots.
### Inputs: data/geno/plink/... (PLINK bed/bim/fam)
### Outputs: results/PCA/PCA_Geno.eigenvec and related files
### Usage: bash analysis/2.PCA_Genetic_variants.sh
### Prerequisites: PLINK installed and in PATH
