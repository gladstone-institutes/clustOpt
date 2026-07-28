# Save an object to RDS atomically

Writes to a sibling `.tmp` file then renames, so a killed job never
leaves a half-written checkpoint. Only called for paths that do not yet
exist, so the rename never has to clobber a destination.

## Usage

``` r
.save_rds_atomic(obj, path)
```

## Arguments

- obj:

  Object to serialize.

- path:

  Destination path.

## Value

`path`, invisibly.
