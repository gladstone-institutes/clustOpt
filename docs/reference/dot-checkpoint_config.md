# Build the run configuration fingerprint stored in the manifest

Captures the run-defining arguments plus a cheap input signature (dims,
assay names, and a hash of the subject id column). The expression matrix
is deliberately not hashed.

## Usage

``` r
.checkpoint_config(
  input,
  ndim,
  dtype,
  sketch_size,
  skip_sketch,
  subject_ids,
  res_range,
  within_batch,
  num_trees,
  train_with,
  min_cells
)
```
