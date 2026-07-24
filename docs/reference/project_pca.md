# project_pca

Project Training and Test Seurat Objects onto Principal Components

This function projects both training and test Seurat objects onto a set
of principal components derived from the training data.

## Usage

``` r
project_pca(
  train_seurat,
  test_seurat,
  train_with_pcs,
  clust_pcs,
  dtype,
  verbose = 0,
  compute_train_eval = FALSE
)
```

## Arguments

- train_seurat:

  A Seurat object representing the training data set.

- test_seurat:

  A Seurat object representing the test data set.

- train_with_pcs:

  Which reduction should be used for training "odd_pca" or "even_pca"

- clust_pcs:

  Which reduction was used for clustering

- dtype:

  Type of data in the Seurat object "scRNA" or "CyTOF"

- verbose:

  Integer verbosity level (0 = silent, 1 = milestones, 2 = detailed, 3 =
  includes Seurat output, 4 = includes output from other packages such
  as ranger)

- compute_train_eval:

  Logical; if `TRUE`, also compute the training data projected onto
  clustering PCs (`train_proj_clust_pcs`). Default is `FALSE` because
  this projection is not used in the standard
  [`clust_opt`](https://gladstone-institutes.github.io/clustOpt/reference/clust_opt.md)
  pipeline.

## Value

A list containing projected data frames/matrices:

- train_proj_train_with_pcs:

  Training data projected onto training PCs (data.frame)

- test_proj_train_with_pcs:

  Test data projected onto training PCs (data.frame)

- test_proj_clust_pcs:

  Test data projected onto clustering PCs (matrix)

- train_proj_clust_pcs:

  Training data projected onto clustering PCs (data.frame); only present
  when `compute_train_eval = TRUE`

## Details

Identifies features that are common between the training and test data
sets. Extracts the PCA loadings from the training data for common
features. Both training and test data are projected onto these loadings.
Data is projected onto loadings for 2 PC sets (odd and even)
