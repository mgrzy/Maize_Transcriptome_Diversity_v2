# Load the required package
library(GenomicRanges)


# Example GRanges object
gr <- read.table("../Genomics/Maize_Genetic_Variants_v5/code/HCCscrpts/meta/ChrReg5Mb.txt")

head(gr)
gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = T, start.field = "V2", end.field = "V3", seqnames.field = "V1")

# Set the tile width to 250kb
tile_width <- 50000

# Split the ranges into 5MB tiles
tiled_gr <- tile(gr, width = tile_width)

# Unlist to flatten the GRanges object
tiled_gr <- unlist(tiled_gr)

# Add the original metadata to the tiled GRanges object
# Repeat the metadata to match the number of tiles
mcols(tiled_gr)$V4 <- rep(mcols(gr)$V4, elementNROWS(tile(gr, width = tile_width)))

# Check the resulting tiled GRanges
tiled_gr
gr <- as.data.frame(tiled_gr)

gr <- gr[,c(1,2,3,6)]

# Select first ranges for each DBs file

gr2 <- c()
for (i in unique(gr$V4)){
  new <- gr[gr$V4==i,][1:4,]
  gr2 <- rbind.data.frame(gr2, new)
}

write.table(gr2, "data/meta/PixyRanges.txt", row.names = F, quote = F, col.names = F, sep = "\t")
### Script: processing/3.Prep_Ranges_Pixy.R
### Purpose: Build genomic windows/regions and population files for pixy to compute diversity statistics.
### Inputs: variant metadata and population assignments
### Outputs: pixy input ranges/pop files in data/meta or results/
### Usage: Rscript processing/3.Prep_Ranges_Pixy.R
### Dependencies: tidyverse
