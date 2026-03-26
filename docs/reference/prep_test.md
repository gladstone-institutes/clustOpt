# prep_test

Prepare test data for random forest

## Usage

``` r
prep_test(input, subject_ids, dtype = "scRNA", test_id, verbose = 0)
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
  includes Seurat output)

## Value

Training data formatted for sil_score, format depends on dtype
