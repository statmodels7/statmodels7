# The Per-Observation Information Blocks

Assembles \\\Omega_i\\, the \\K \times K\\ information of observation
\\i\\ in the link-scale predictors, weighted by the prior weight.

## Usage

``` r
info_blocks(spec, theta, expected = TRUE, approx = "bartlett")
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- theta:

  The per-observation parameters.

- expected:

  Whether the expected information is wanted.

- approx:

  The approximation, where the expected one is not closed.

## Value

An \\n \times K \times K\\ array.
