# calculate_modularity

Calculate modularity score for a graph given an adjacency matrix and
cluster assignments

## Usage

``` r
calculate_modularity(adj_matrix, clusters, directed = FALSE)
```

## Arguments

- adj_matrix:

  A sparse adjacency matrix (Matrix package format)

- clusters:

  Vector of cluster assignments for each node

- directed:

  Logical, whether the graph is directed (default: FALSE)

## Value

Modularity value (numeric)
