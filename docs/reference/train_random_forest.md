# train_random_forest

Train the random forest and predict on the test sample

## Usage

``` r
train_random_forest(
  res,
  df_list,
  cluster_by_res,
  sam,
  num_trees,
  verbose = 0,
  snn_graph = NULL,
  precomputed_dist = NULL,
  rf_num_threads = 1
)
```

## Arguments

- res:

  Resolution to train on

- df_list:

  List containing training and test data

- cluster_by_res:

  Named list of cluster assignment vectors, keyed by resolution (as
  character). Pre-extracted from training metadata using exact column
  names.

- sam:

  Test sample

- num_trees:

  Number of trees for the random forest

- verbose:

  Integer verbosity level (0 = silent, 1 = milestones, 2 = detailed, 3 =
  includes Seurat output, 4 = includes output from other packages such
  as ranger)

- snn_graph:

  Precomputed SNN graph (sparse matrix). Required.

- precomputed_dist:

  Optional precomputed distance matrix (output of
  [`dist`](https://rdrr.io/r/stats/dist.html)). Passed through to
  [`calculate_silhouette_score`](https://gladstone-institutes.github.io/clustOpt/reference/calculate_silhouette_score.md).

- rf_num_threads:

  Number of threads for ranger (default 1).

## Value

A list containing the resolution, silhouette score, and number of
predicted clusters.
