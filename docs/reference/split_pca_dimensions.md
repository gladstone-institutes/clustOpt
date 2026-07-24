# split_pca_dimensions

Takes a Seurat object with an existing PCA reduction and splits the
dimensions into 2 sets odd and even PCs.

## Usage

``` r
split_pca_dimensions(input, verbose = 0)
```

## Arguments

- input:

  Seurat object

- verbose:

  Integer verbosity level (0 = silent, 1 = milestones, 2 = detailed, 3 =
  includes Seurat output, 4 = includes output from other packages such
  as ranger)

## Value

Seurat object with new PCA reductions
