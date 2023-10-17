library(Seurat)
options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
library(BPCells)
library(Azimuth)
library(tidyverse)
library(biomaRt)


# For orthology mapping of mouse genes to human
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Read in mouse data ------------------------------------------------------

mouse_bm <- open_matrix_10x_hdf5(paste0(
  "/Users/nelphick/Dropbox (Gladstone)/Work/Test_data/",
  "mouse_bone_marrow_annotated_tabula_muris_cellxgene/",
  "local.h5ad"
))
write_matrix_dir(
  mat = mouse,
  dir = paste0(
    "~/Dropbox (Gladstone)/Work/Test_data/",
    "mouse_bone_marrow_annotated_tabula_muris_cellxgene/",
    "on_disk_mat"
  ),
  overwrite = T
)

mouse <- Azimuth:::ConvertEnsembleToSymbol(mat = mouse, species = "mouse")


# Have to use older archives here since the getLDS function fails for newer ones
human_mart <- useMart("ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "https://dec2021.archive.ensembl.org"
)
mouse_mart <- useMart("ensembl",
  dataset = "mmusculus_gene_ensembl",
  host = "https://dec2021.archive.ensembl.org"
)

gene_mapping <- getLDS(
  attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
  values = rownames(mouse), mart = human_mart,
  attributesL = c("mgi_symbol"), martL = mouse_mart,
  uniqueRows = T
)
rm(human_mart)
rm(mouse_mart)
duplicate_mappings <- gene_mapping %>%
  group_by(HGNC.symbol) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  pull(HGNC.symbol)

mouse_duplicate_mappings <- gene_mapping %>%
  group_by(MGI.symbol) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  pull(MGI.symbol)

# Remove duplicate mappings
gene_mapping <- gene_mapping %>%
  filter(!(HGNC.symbol %in% duplicate_mappings)) %>%
  filter(!(MGI.symbol %in% mouse_duplicate_mappings))


# Only ~50% of genes map uniquely to human ortholouges, the final dataset
# should keep only these from the human data to avoid a bias in the number of
# features.

mouse <- subset(mouse, rownames(mouse) %in% gene_mapping$MGI.symbol)

rownames(mouse) <- gene_mapping$HGNC.symbol[match(rownames(mouse), gene_mapping$MGI.symbol)]

mouse <- CreateSeuratObject(mouse)

dir.create("~/Dropbox (Gladstone)/Work/Test_data/1Million_mouse_brain_cells_10X/mouse_data",
  recursive = T,
  showWarnings = F
)

saveRDS(
  object = mouse,
  file = "1M_neurons_mouse_with_human_gene_names.Rds",
  destdir = "~/Dropbox (Gladstone)/Work/Test_data/1Million_mouse_brain_cells_10X/mouse_data"
)



# Read in human data ------------------------------------------------------

file.dir <- "~/Dropbox (Gladstone)/Work/Test_data/1M_human_PBMC/"
files.set <- c("ahern_pbmc.h5ad", "jin_pbmc.h5ad", "yoshida_pbmc.h5ad")
# Loop through h5ad files and output BPCells matrices on-disk
data.list <- c()
metadata.list <- c()

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i])
  data <- open_matrix_anndata_hdf5(path)
  write_matrix_dir(
    mat = data,
    dir = paste0(gsub(".h5ad", "", path), "_BP")
  )
  # Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(gsub(".h5ad", "", path), "_BP"))
  mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  # Get metadata
  metadata.list[[i]] <- LoadH5ADobs(path = path)
  data.list[[i]] <- mat
}
# Name layers
names(data.list) <- c("ahern", "jin", "yoshida")

# Add Metadata
for (i in 1:length(metadata.list)) {
  metadata.list[[i]]$publication <- names(data.list)[i]
}
metadata.list <- lapply(metadata.list, function(x) {
  x <- x[, c("publication", "sex", "cell_type", "donor_id", "disease")]
  return(x)
})
metadata <- Reduce(rbind, metadata.list)

# Save the human dataset on its own
human <- CreateSeuratObject(counts = data.list, meta.data = metadata)
rm(data.list)
rm(metadata)
dir.create("~/Dropbox (Gladstone)/Work/Test_data/1M_human_PBMC/human_data_all_genes",
  recursive = T,
  showWarnings = F
)

saveRDS(
  object = human,
  file = "1M_human_PBMCs.Rds",
  destdir = "~/Dropbox (Gladstone)/Work/Test_data/1M_human_PBMC/human_data_all_genes"
)

human <- human[rownames(human) %in% rownames(mouse), ]

human <- JoinLayers(human)

human <- AddMetaData(
  object = human,
  metadata = as.factor(rep(
    "Human",
    ncol(human)
  )),
  col.name = "Species"
)
mouse <- AddMetaData(
  object = mouse,
  metadata = as.factor(rep(
    "Mouse",
    ncol(mouse)
  )),
  col.name = "Species"
)

merged_seurat <- merge(human, y = mouse)
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat,
                                             pattern = "^MT-"
)

saveRDS(
  object = merged_seurat,
  file = "2_million_cells_mouse_neurons_human_pmbcs.Rds",
  destdir = "~/Dropbox (Gladstone)/Work/Test_data/Merged_mouse_and_human_data"
)


