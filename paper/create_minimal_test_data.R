library(Seurat)
options(Seurat.object.assay.version = "v5")
library(BPCells)
library(dplyr)
library(readr)
set.seed(42)

# Functions ---------------------------------------------------------------

# function to remove unused factor levels from the seurat metadata
refine_metadata_levels <- function(seurat_data) {
  for (i in base::colnames(seurat_data@meta.data)) {
    if (base::is.factor(seurat_data@meta.data[[i]])) {
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse = ", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse = ", ")
        )
      )
    }
  }
  return(seurat_data)
}


# Read in input data ------------------------------------------------------

# ROSMAP dataset
# https://www.synapse.org/#!Synapse:syn3219045
counts <- ReadMtx(
  mtx = "~/Dropbox (Gladstone)/Work/Test_data/MathysEtAl/filtered_count_matrix.mtx",
  cells = "~/Dropbox (Gladstone)/Work/Test_data/MathysEtAl/filtered_column_metadata.txt",
  features = "~/Dropbox (Gladstone)/Work/Test_data/MathysEtAl/filtered_gene_row_names.txt",
  feature.column = 1,
  skip.cell = 1
)

mdata <- read.delim("~/Dropbox (Gladstone)/Work/Test_data/MathysEtAl/filtered_column_metadata.txt")

rownames(mdata) <- mdata[, 1]
mdata <- mdata[, -1]

human <- CreateSeuratObject(counts = counts, meta.data = mdata)
rm(mdata)
rm(counts)
# Convert to on-disk counts
write_matrix_dir(
  mat = human[["RNA"]]$counts,
  dir = "~/Dropbox (Gladstone)/Work/Test_data/MathysEtAl_counts",
  overwrite = TRUE
)
h_data <- open_matrix_dir(dir = paste0("~/Dropbox (Gladstone)/Work/Test_data",
                                       "/MathysEtAl_counts"))

human[["RNA"]]$counts <- h_data
rm(h_data)


# Pre-sketch metadata
pre_sketch_mdata <- table(human@meta.data$broad.cell.type,
                          human@meta.data$projid) |>
  as.data.frame() |>
  rename(broad.cell.type = Var1, projid = Var2) 



human <- NormalizeData(human) |>
  FindVariableFeatures()
# Sketch 1000 cells
human <- SketchData(object = human, ncells = 1000)


DefaultAssay(human) <- "sketch"
human <- DietSeurat(human,
                    assays = "sketch",
                    layers = c("counts", "data"))


human <- RenameAssays(object = human, sketch = "RNA")

human <- refine_metadata_levels(human)

# Pre-sketch metadata
post_sketch_mdata <- table(human@meta.data$broad.cell.type,
                           human@meta.data$projid) |>
  as.data.frame() |>
  rename(broad.cell.type = Var1, projid = Var2)

dir.create("~/Dropbox (Gladstone)/Work/Test_data/minimal_test_data_rosmap",
  recursive = T,
  showWarnings = F
)

write_csv(x = pre_sketch_mdata,
          paste0("~/Dropbox (Gladstone)/Work/Test_data/",
                 "minimal_test_data_rosmap/pre-sketch_metadata.csv"))
write_csv(x = post_sketch_mdata,
          paste0("~/Dropbox (Gladstone)/Work/Test_data/",
                 "minimal_test_data_rosmap/post-sketch_metadata.csv"))
saveRDS(
  object = human,
  file = "~/Dropbox (Gladstone)/Work/Test_data/minimal_test_data_rosmap/1000_sketch_ROSMAP.Rds"
)
