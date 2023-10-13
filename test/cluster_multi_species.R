library(Seurat)
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
library(BPCells)
library(Azimuth)
library(tidyverse)
library(biomaRt)

set.seed(42)



# Functions ---------------------------------------------------------------

plot_metadata_sketch <- function(data,output) {
  # get the metadata columns that not numeric
  meta_cols <- names(data@meta.data)[!sapply(data@meta.data, is.numeric)]
  # Exclude "orig.ident"
  meta_cols <- meta_cols[meta_cols != "orig.ident"]
  for (meta in meta_cols) {
    # Skip if there is only one value
    if (length(unique(data@meta.data[[meta]])) == 1) {
      next
    }
    print(paste0("***** Plotting ", meta, " *****"))
    png(file.path(paste0(output, "/", meta, "_sketch_UMAP.png")),
        width = 10,
        height = 10,
        units = "in",
        res = 300
    )
    print(DimPlot(data,
                  raster = FALSE,
                  pt.size = .75,
                  order = TRUE,
                  label = FALSE,
                  group.by = meta,
                  reduction = "umap"
    ))
    dev.off()
    
  }
}
plot_metadata_full <- function(data,output) {
  # get the metadata columns that not numeric
  meta_cols <- names(data@meta.data)[!sapply(data@meta.data, is.numeric)]
  # Exclude "orig.ident"
  meta_cols <- meta_cols[meta_cols != "orig.ident"]
  for (meta in meta_cols) {
    # Skip if there is only one value
    if (length(unique(data@meta.data[[meta]])) == 1) {
      next
    }
    print(paste0("***** Plotting ", meta, " *****"))
    png(file.path(paste0(output, "/", meta, "_full_UMAP.png")),
        width = 10,
        height = 10,
        units = "in",
        res = 300
    )
    print(DimPlot(data,
                  raster = FALSE,
                  alpha = 0.05,
                  pt.size = .75,
                  order = TRUE,
                  label = FALSE,
                  group.by = meta,
                  reduction = "full.umap"
    ))
    dev.off()
    
  }
}



# Read in and normalize data ----------------------------------------------

merged_seurat <- readRDS("~/Dropbox (Gladstone)/Work/Test_data/Merged_mouse_and_human_data/2_million_cells_mouse_neurons_human_pmbcs.Rds")
merged_seurat <- JoinLayers(merged_seurat)

# Filter out percent.mt > 15%
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 &
  percent.mt < 15)
ncol(merged_seurat)
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)

# Sketch data -------------------------------------------------------------

# Adds the "sketch" assay which is in-memory, we can then run the standard
# clustering workflow (including SCTransform which will run faster on the
# sketch)
merged_seurat <- SketchData(
  object = merged_seurat,
  ncells = 280000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(merged_seurat) <- "sketch"

# We redo the variable features since they might be slightly different for
# this set.
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, npcs = 30)
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30)

merged_seurat <- FindClusters(merged_seurat,resolution = 0.6)
# We have to save the model to project the full data later
merged_seurat <- RunUMAP(merged_seurat,
  reduction = "pca",
  dims = 1:30,
  return.model = TRUE
)


# Project onto the full dataset -------------------------------------------

merged_seurat <- ProjectData(
  object = merged_seurat,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  refdata = list(cluster_full = "seurat_clusters"),
  dims = 1:30
)

plot_metadata_sketch(merged_seurat, output = "~/tmp/")


# Now that we have projected the full dataset, switch back to analyzing all
# cells (on-disk) to save memory
DefaultAssay(merged_seurat) <- "RNA"

merged_seurat <- DietSeurat(merged_seurat,
                            assays = "RNA",
                            dimreducs = c("full.umap",
                                          "full.pca"))
gc()


plot_metadata_full(merged_seurat, output = "~/tmp/")

png(file.path(paste0("~/tmp", "/", "labeled_clusters_full_UMAP.png")),
    width = 10,
    height = 10,
    units = "in",
    res = 300
)
print(DimPlot(merged_seurat,
        raster = FALSE,
        alpha = 0.05,
        pt.size = .75,
        order = TRUE,
        label = TRUE,
        group.by = "cluster_full",
        reduction = "full.umap"))
dev.off()

png(file.path(paste0("~/tmp", "/", "cell_type_full_UMAP.png")),
    width = 20,
    height = 10,
    units = "in",
    res = 300
)
print(DimPlot(merged_seurat,
              raster = FALSE,
              alpha = 0.05,
              pt.size = .75,
              order = TRUE,
              label = FALSE,
              group.by = "cell_type",
              reduction = "full.umap"))
dev.off()

table(merged_seurat@meta.data$Species)

saveRDS(
  object = merged_seurat,
  file = "2_million_cells_mouse_neurons_human_pmbcs_clustered.Rds",
  destdir = "~/Dropbox (Gladstone)/Work/Test_data/Merged_mouse_and_human_data"
)

merged_seurat <- readRDS("~/Dropbox (Gladstone)/Work/Test_data/Merged_mouse_and_human_data/2_million_cells_mouse_neurons_human_pmbcs_clustered.Rds")

merged_seurat <- subset(merged_seurat,subset = cluster_full == 35)
table(merged_seurat@meta.data$cell_type) %>%
  as.data.frame() %>%
  write_csv(paste0("~/tmp/cluster_35_cell_counts.csv"))

png(file.path(paste0("~/tmp", "/", "35_cell_type_full_UMAP.png")),
    width = 15,
    height = 10,
    units = "in",
    res = 300
)
print(DimPlot(merged_seurat,
              raster = FALSE,
              pt.size = 1,
              order = TRUE,
              label = FALSE,
              group.by = "cell_type",
              reduction = "full.umap") +
        xlim(2,5) + 
        ylim(NA,-7))
dev.off()

merged_seurat <- subset(merged_seurat,subset = Species == "Human")
png(file.path(paste0("~/tmp", "/", "Human_only_cell_type.png")),
    width = 20,
    height = 10,
    units = "in",
    res = 300
)
print(DimPlot(merged_seurat,
              raster = FALSE,
              pt.size = .01,
              order = TRUE,
              label = FALSE,
              group.by = "cell_type",
              reduction = "full.umap"))
dev.off()