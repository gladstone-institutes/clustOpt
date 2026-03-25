# convert_seurat_to_bpcells

Convert a Seurat object to BPCells on-disk format

Converts the counts matrix of specified assays in a Seurat object to
BPCells format, saving the matrices on disk and updating the Seurat
object to use these on-disk matrices. Only works for single count
layers.

## Usage

``` r
convert_seurat_to_bpcells(seurat_obj, output_dir = NULL, assays = "RNA")
```

## Arguments

- seurat_obj:

  A Seurat object to be converted.

- output_dir:

  Directory where the BPCells matrices will be saved. Defaults to a
  subdirectory in the system's temporary directory named after the
  Seurat object.

- assays:

  Character vector of assays to convert. Defaults to "RNA".

## Value

The updated Seurat object using on-disk matrices.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
seurat_obj <- readRDS("/path/to/your/seurat_object.rds")
seurat_obj <- convert_seurat_to_bpcells(seurat_obj)
} # }
```
