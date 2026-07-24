# prep_train

Prepare training data for random forest

## Usage

``` r
prep_train(
  input,
  subject_ids,
  dtype = "scRNA",
  within_batch = NA,
  test_id,
  verbose = 0
)
```

## Arguments

- input:

  Seurat object

- subject_ids:

  Metadata field that identifies unique samples.

- dtype:

  Type of data in the Seurat object "scRNA" or "CyTOF", default is
  "scRNA". CyTOF data is expected to be arcsinh normalized.

- within_batch:

  Batch variable, for a given sample only those with the same value for
  the batch variable will be used for training.

- test_id:

  subject_id for the test sample

- verbose:

  Integer verbosity level (0 = silent, 1 = milestones, 2 = detailed, 3 =
  includes Seurat output, 4 = includes output from other packages such
  as ranger)

## Value

Training data formatted for sil_score, format depends on dtype
