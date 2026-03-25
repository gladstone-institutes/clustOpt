# Calculate modularity for a graph given coordinate matrix and cluster assignments

Calculate modularity for a graph given coordinate matrix and cluster
assignments

## Usage

``` r
calculate_modularity_from_coords(
  clusters,
  coords,
  k = 20,
  mutual = FALSE,
  distance_metric = "euclidean",
  directed = FALSE
)
```

## Arguments

- clusters:

  Vector of cluster assignments for each node

- coords:

  Numeric matrix where rows are nodes and columns are coordinates

- k:

  Number of nearest neighbors for kNN graph construction

- mutual:

  Logical, whether to create mutual kNN graph (default: FALSE)

- distance_metric:

  Distance metric to use ("euclidean", "manhattan", "cosine")

- directed:

  Logical, whether the graph is directed (default: FALSE)

## Value

List containing modularity value and the kNN adjacency matrix
