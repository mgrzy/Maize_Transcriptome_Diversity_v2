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

conda activate mkl

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

for task_id in {1..10}; do
  ../../soft/osca/osca --eqtl --bfile data/processed/geno/plink/Div639 \
  --befile data/processed/osca/WiDiv631 \
  --grm data/processed/osca/WiDiv631_Kin \
  --qcovar data/processed/osca/WiDiv631_PC.txt \
  --mlm --task-num 10 --task-id "$task_id" --thread-num 36 --cis --cis-wind 2000 --no-fid --out results/full/osca/cis2Mb
  echo "OSCA task ${task_id}/10 completed."
done
