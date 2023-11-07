library(Seurat)
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
library(BPCells)
library(Azimuth)
library(tidyverse)

set.seed(42)

# function to remove unused factor levels from the test data
refine_metadata_levels <- function(seurat_data){
  for (i in base::colnames(seurat_data@meta.data)){
    if (base::is.factor(seurat_data@meta.data[[i]])){
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
    }
  }
  return (seurat_data)
}

input <- "/Users/nelphick/Dropbox (Gladstone)/Work/Test_data/mapping_the_developing_human_immune_system_across_organs_HCA/9fcfb360-e8b2-47b9-a9ed-7af4d06b7f04/Visium10X_data_LI.h5ad"

h_data <- open_matrix_anndata_hdf5(
  path = input
)
count_dir <- paste0("/Users/nelphick/Dropbox (Gladstone)/Work/Test_data/mapping_the_developing_human_immune_system_across_organs_HCA/9fcfb360-e8b2-47b9-a9ed-7af4d06b7f04/","/raw_count_dir")
dir.create(count_dir,showWarnings = F, recursive = T)

# Write the matrix to a directory
write_matrix_dir(
  mat = h_data,
  overwrite = TRUE,
  dir = count_dir
)

h_data <- open_matrix_dir(dir = count_dir)
h_data <- Azimuth:::ConvertEnsembleToSymbol(mat = h_data, species = "human")


# Save the human dataset on its own
human <- CreateSeuratObject(counts = h_data, meta.data = LoadH5ADobs(input))
rm(data.list)
rm(metadata)


# Select only the control samples for each study
human <- subset(human, subset = disease == "normal")

# Get the 5 donors from each study with the highest cell counts
keep_donors <- human@meta.data |>
  as.data.frame() |>
  dplyr::select(donor_id, publication) |>
  group_by(donor_id, publication) |>
  summarise(n = n()) |>
  top_n(n = 5,wt = n) |>
  pull(donor_id) |>
  as.character()

Idents(human) <- "donor_id"
keep_cells <- WhichCells(human,idents = keep_donors)
human <- subset(human,cells = keep_cells)

human <- refine_metadata_levels(human)

human <- JoinLayers(human,assay = "RNA")

human <- NormalizeData(human) |>
  FindVariableFeatures()

human <- SketchData(human,ncells = 10000,verbose = T)




dir.create("~/Dropbox (Gladstone)/Work/Test_data/minimal_test_data_pbmc",
           recursive = T,
           showWarnings = F
)


saveRDS(
  object = human,
  file = "1M_human_PBMCs.Rds",
  destdir = "~/Dropbox (Gladstone)/Work/Test_data/1M_human_PBMC/human_data_all_genes"
)


