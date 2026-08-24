###
###
###
library(data.table)
library(tidyverse)

### Load data
lf <- list.files("results/Pi/raw", pattern = "*pi.txt", full.names = T)
pip <- lapply(lf, fread)
pip <- rbindlist(pip)
pip <- na.omit(pip)
gc()

fwrite(pip, "results/Pi/Pi_all.csv")
### Script: processing/9.Process_Pixy_Results.R
### Purpose: Parse and summarize pixy outputs (pi/diversity etc.) into tidy tables per population and window.
### Inputs: pixy output files
### Outputs: results/Pi/* (summarized diversity metrics)
### Usage: Rscript processing/9.Process_Pixy_Results.R
### Dependencies: tidyverse, data.table
