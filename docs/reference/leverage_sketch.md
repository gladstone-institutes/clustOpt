# leverage_sketch

Uses leverage score-based sampling to reduce the size of large Seurat
objects by creating a representative sketch assay. This method preserves
the most informative cells while dramatically reducing computational
requirements. Supports both single-cell RNA-seq and CyTOF data.

## Usage

``` r
leverage_sketch(
  input,
  sketch_size,
  dtype = "scRNA",
  skip_norm = FALSE,
  on_disk = FALSE,
  output_dir = NULL,
  verbose = 1
)
```

## Arguments

- input:

  A Seurat object to be sketched

- sketch_size:

  Integer. Number of cells to include in the sketch assay. If NULL,
  defaults to 10% of total cells

- dtype:

  Character. Type of data: "scRNA" (default) for single-cell RNA-seq or
  "CyTOF" for mass cytometry. CyTOF data should be arcsinh normalized
  and stored in the counts slot

- skip_norm:

  Logical. Set to TRUE if scRNA-seq data has already been normalized
  with \`Seurat::NormalizeData()\` (default FALSE). CyTOF data
  normalization is always skipped

- on_disk:

  Logical. Whether to use BPCells on-disk count matrices to speed up
  sketching for very large datasets (default FALSE)

- output_dir:

  Character. Directory path for storing on-disk count matrices when
  \`on_disk = TRUE\`. If NULL, uses temporary directory

- verbose:

  Integer verbosity level: 0 = silent, 1 = milestones, 2 = detailed
  progress, 3 = includes Seurat function output (default 1)

## Value

A Seurat object containing only the sketch assay, renamed to "RNA" for
compatibility with downstream functions

## Details

Large datasets (\>200,000 cells) benefit from \`on_disk = TRUE\` to
reduce memory usage during sketching.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic sketching for scRNA-seq data (uses variable features)
sketched_obj <- leverage_sketch(seurat_obj,
  sketch_size = 5000,
  dtype = "scRNA"
)

# Sketching CyTOF data (uses ALL features from marker panel)
cytof_sketch <- leverage_sketch(cytof_obj,
  sketch_size = 2000,
  dtype = "CyTOF"
)

# Large dataset with on-disk matrices
large_sketch <- leverage_sketch(large_obj,
  sketch_size = 10000,
  on_disk = TRUE, verbose = 2
)
} # }
```
