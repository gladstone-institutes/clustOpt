# calculate_kl_divergence

Calculate Kullback-Leibler divergence between two probability
distributions. D_KL(Q \|\| P) measures how much information is lost when
P is used to approximate Q. The result is always non-negative, with 0
indicating identical distributions.

## Usage

``` r
calculate_kl_divergence(q, p)
```

## Arguments

- q:

  Numeric vector representing the "true" probability distribution. Must
  be non-negative and sum to 1.

- p:

  Numeric vector representing the "approximating" probability
  distribution. Must be non-negative, sum to 1, and have the same length
  as q.

## Value

Numeric value representing D_KL(Q \|\| P), always \>= 0. Returns 0 when
distributions are identical.

## Details

The KL divergence is calculated as: D_KL(Q \|\| P) = sum(q \* log(q /
p))

Note that KL divergence is asymmetric: D_KL(Q \|\| P) != D_KL(P \|\| Q).
When q\[i\] \> 0 but p\[i\] = 0, the divergence is infinite. This
implementation requires all elements of p to be positive when
corresponding elements of q are positive.

## Examples

``` r
# Identical distributions
p <- c(0.25, 0.25, 0.25, 0.25)
calculate_kl_divergence(p, p)  # Returns 0
#> [1] 0

# Different distributions
q <- c(0.5, 0.3, 0.1, 0.1)
p <- c(0.25, 0.25, 0.25, 0.25)
calculate_kl_divergence(q, p)  # Returns positive value
#> [1] 0.2180119
```
