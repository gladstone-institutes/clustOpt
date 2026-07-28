# Deterministic per-subject seed

Derives a stable integer seed from the run seed and the subject name, so
a subject reproduces exactly regardless of which subjects ran before it
in the process (the basis for bit-identical resume).

## Usage

``` r
.subject_seed(seed, sam)
```
