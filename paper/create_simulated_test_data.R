library(Seurat)
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
library(BPCells)
library(tidyverse)
library(splatter)
library(SingleCellExperiment)

set.seed(1234)


params.batches <- newSplatPopParams(similarity.scale = 4,
                                    batchCells = c(20, 20),
                                    batch.size = 4,
                                    batch.facLoc = 0.1,
                                    batch.facScale = 0.2)

sim <- splatPopSimulate(params = params.batches, sparsify = TRUE)


sim <- CreateSeuratObject(sim@assays@data$counts,
  meta.data = as.data.frame(sim@colData)
)
sim <- NormalizeData(sim)
sim <- ScaleData(sim)
sim <- FindVariableFeatures(sim)
sim <- RunPCA(sim,layer = "data")


sim <- RunUMAP(sim, dims = 1:10)
