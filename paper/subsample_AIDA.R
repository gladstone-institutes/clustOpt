#!/usr/bin/env Rscript

###############################################################################
## Author: Natalie Elphick
##
## Script Goal: Create leverage score based subsamples of the 1M PBMC Asian
## Immune Diversity Atlas (AIDA)
##
## Usage example:
## Rscript subsample_AIDA.R \
## --input AIDA.h5ad \ 
## --output_dir /path/to/output/directory          
###############################################################################

# Get input arguments -----------------------------------------------------
library(optparse)
option_list <- list(
  make_option("--input",
    action = "store", default = NA, type = "character",
    help = "Seurat object with more than 1 sample  (required)"
  ),
  make_option("--output_dir",
    action = "store", default = NA, type = "character",
    help = "Output directory (required)"
  )
)
# Read in the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Check that all required arguments are present
if (is.na(opt$input) | is.na(opt$output_dir) ) {
  stop("***** ERROR: Missing required arguments! *****")
}


library(Seurat)
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
library(BPCells)
library(Azimuth)
library(tidyverse)

h_data <- open_matrix_anndata_hdf5(
  path = opt$input
)
count_dir <- paste0(opt$output_dir,"/raw_count_dir")
dir.create(count_dir,showWarnings = F, recursive = T)

# Write the matrix to a directory
write_matrix_dir(
  mat = h_data,
  overwrite = TRUE,
  dir = count_dir
)

h_data <- open_matrix_dir(dir = count_dir)
h_data <- Azimuth:::ConvertEnsembleToSymbol(mat = h_data, species = "human")

human_pbmc <- CreateSeuratObject(h_data,
  meta.data = LoadH5ADobs(path = opt$input)
)
rm(h_data)

human_pbmc <- NormalizeData(human_pbmc) |>
  FindVariableFeatures()

for (n in c(
  floor(ncol(human_pbmc) * 0.01), # 1%
  floor(ncol(human_pbmc) * 0.05), # 5%
  floor(ncol(human_pbmc) * 0.25), # 25%
  floor(ncol(human_pbmc) * 0.5),  # 50%
  floor(ncol(human_pbmc) * 0.75)  # 75%
)) {
  small_pbmc <- SketchData(human_pbmc, ncells = n, over.write = T)

  small_pbmc <- CreateSeuratObject(small_pbmc[["sketch"]]$counts,
    meta.data = small_pbmc@meta.data[colnames(small_pbmc[["sketch"]]$counts), ]
  )
  dir <- paste0(opt$output_dir,
                n,"_Cells_AIDA_counts")
  dir.create(dir,
    showWarnings = F,
    recursive = T
  )
  saveRDS(
    object = small_pbmc,
    file = paste0(n,"_cells_AIDA.Rds"),
    destdir = dir
  )
}
