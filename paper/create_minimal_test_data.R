library(Seurat)
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
library(BPCells)

set.seed(42)

# Functions ---------------------------------------------------------------

# function to remove unused factor levels from the seurat metadata
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


# Read in input data ------------------------------------------------------

# https://data.humancellatlas.org/explore/projects/f2078d5f-2e7d-4844-8552-f7c41a231e52
input <- "~/Dropbox (Gladstone)/Work/Test_data/eQTLAutoimmune/eQTLAutoimmune.h5ad"

h_data <- open_matrix_anndata_hdf5(
  path = input
)
count_dir <- paste0("~/Dropbox (Gladstone)/Work/Test_data/eQTLAutoimmune/",
                    "/raw_count_dir")

# Write the matrix to a directory
write_matrix_dir(
  mat = h_data,
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


