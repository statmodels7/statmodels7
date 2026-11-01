# The Weights of a Jeffreys Block Inside a Filter Along One Direction

Returns \\g\_{i,r} = \sum_q\mathcal I\_{qr}Z\_{q,i}\\ for every equation
\\r\\, the per-observation weights against which the derivatives of the
filtered predictors along the whitened direction \\d_i\\ are contracted.

## Usage

``` r
jfp_filter_weights(st, bs, i, Z = bs$Z)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- bs:

  A block's entry of the state, of kind `"filter"`.

- i:

  The direction.

- Z:

  The whitened rows to read, `bs$Z` by default.

## Value

A named list of numeric vectors, one per distribution parameter.
