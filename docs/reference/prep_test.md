# prep_test

Prepare test data for random forest

## Usage

``` r
prep_test(
  input,
  subject_ids,
  dtype = "scRNA",
  test_id,
  verbose = 0,
  residual_features = NULL
)
```

## Arguments

- input:

  Seurat object

- subject_ids:

  Metadata field that identifies unique subjects.

- dtype:

  Type of data in the Seurat object "scRNA" or "CyTOF", default is
  "scRNA". CyTOF data is expected to be arcsinh normalized.

- test_id:

  subject_id for the test sample

- verbose:

  Integer verbosity level (0 = silent, 1 = milestones, 2 = detailed, 3 =
  includes Seurat output, 4 = includes output from other packages such
  as ranger)

- residual_features:

  Character vector of features to restrict the SCT residual computation
  to (passed to `SCTransform(residual.features=)`). Only used for scRNA
  data. Default `NULL` preserves prior behavior (residuals computed for
  all genes).

## Value

Training data formatted for sil_score, format depends on dtype
