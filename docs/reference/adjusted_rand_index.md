# Calculate Adjusted Rand Index

Computes the Adjusted Rand Index (ARI) to measure the similarity between
two clustering assignments. The ARI is a measure of agreement between
two partitions, adjusted for chance. Values range from 0 (random
partitioning) to 1 (perfect agreement), with negative values indicating
worse than random.

## Usage

``` r
adjusted_rand_index(seurat_obj, meta1, meta2)
```

## Arguments

- seurat_obj:

  A Seurat object containing the metadata with clustering assignments

- meta1:

  Character string specifying the first metadata column name containing
  cluster assignments

- meta2:

  Character string specifying the second metadata column name containing
  cluster assignments to compare against meta1

## Value

Numeric value representing the Adjusted Rand Index between the two
clustering assignments

## Details

The Adjusted Rand Index is calculated using the formula: ARI = (RI -
Expected_RI) / (max(RI) - Expected_RI)

Where: - RI is the Rand Index - Expected_RI is the expected value of RI
under random partitioning - max(RI) is the maximum possible value of RI

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare two clustering results in a Seurat object
ari_score <- adjusted_rand_index(seurat_obj, "seurat_clusters", "leiden_clusters")
} # }
```
