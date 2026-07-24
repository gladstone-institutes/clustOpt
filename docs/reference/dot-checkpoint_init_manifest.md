# Initialize or validate the checkpoint manifest

On first run, writes a manifest holding the config fingerprint and the
resolved seed. On resume, validates that the config matches (erroring on
drift) and returns the persisted seed so the sketch and per-subject
seeding are identical to the original run.

## Usage

``` r
.checkpoint_init_manifest(dir, config, seed)
```

## Arguments

- dir:

  Checkpoint directory.

- config:

  Config list from
  [`.checkpoint_config()`](https://gladstone-institutes.github.io/clustOpt/reference/dot-checkpoint_config.md).

- seed:

  Resolved integer seed for this run.

## Value

The seed to actually use (persisted seed on resume, else `seed`).
