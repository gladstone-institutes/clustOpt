# Canonical checkpoint path for one holdout subject

Index is zero-padded to the width of the subject count so files sort in
run order; the (sanitized) subject name is embedded only for
readability. Both index and name are deterministic given an identical
config and persisted sketch, so the exact path is reproducible on
resume.

## Usage

``` r
.subject_checkpoint_path(dir, idx, sam, n_total)
```
