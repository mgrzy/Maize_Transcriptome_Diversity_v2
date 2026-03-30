# Maize Cold Spring

## Overview
This repository contains the **code, configs, and documentation** for the Maize Transcriptome Diversity project.  
The project is executed in two environments:
- **Local (RStudio)** → exploratory analysis, visualization, small-scale testing  
- **HPC (SLURM)** → heavy pipelines, full-scale computations  

Raw and processed data are stored on HPC and synced via Google Drive.  
**GitHub contains only code, configs, and documentation.**

---

## Project structure

- `config/` → environment configs (paths_local.sh, paths_hpc.sh) and sample lists  
- `data/`  
  - `raw/` → raw input data (FASTQ/VCF, ignored in Git)  
  - `processed/` → processed intermediate files (ignored in Git)  
- `docs/` → project documentation and RMarkdown analysis notebooks  
- `logs/` → pipeline and SLURM logs (ignored in Git)  
- `results/`  
  - `full/` → full results (ignored in Git, synced via Drive)  
  - `summary/` → small CSVs/plots for reports (tracked in Git)  
- `scripts/`  
  - `hpc/` → production-level scripts and SLURM workflows  
  - `local/` → exploratory/local analysis scripts  
  - `utils/` → shared helper functions (R/Bash)  
- `README.md` → project overview  
- `.gitignore` - files ignored by git