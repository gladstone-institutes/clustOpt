# calculate_silhouette_score

Calculate silhouette score

## Usage

``` r
calculate_silhouette_score(predicted, data_frame, precomputed_dist = NULL)
```

## Arguments

- predicted:

  Cluster assignments

- data_frame:

  Data frame containing the data

- precomputed_dist:

  Optional precomputed distance matrix (output of
  [`dist`](https://rdrr.io/r/stats/dist.html)). If `NULL` (default), the
  distance matrix is computed from `data_frame`.

## Value

A list containing the average silhouette score and the average
silhouette score for each cluster.
