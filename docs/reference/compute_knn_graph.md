# compute_knn_graph

Compute k-nearest neighbor graph from coordinate matrix

## Usage

``` r
compute_knn_graph(coords, k, mutual = FALSE, distance_metric = "euclidean")
```

## Arguments

- coords:

  Numeric matrix where rows are nodes and columns are coordinates

- k:

  Number of nearest neighbors

- mutual:

  Logical, whether to create mutual kNN graph (default: FALSE)

- distance_metric:

  Distance metric to use ("euclidean", "manhattan", "cosine")

## Value

Sparse adjacency matrix of the kNN graph
