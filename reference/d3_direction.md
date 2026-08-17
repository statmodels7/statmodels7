# The Third Derivative Already Contracted in One Direction

\\\sum_k \ell\_{abk}(X_k v_k)\_i\\, the per-observation weight
[`contract3`](https://statmodels7.github.io/statmodels7/reference/contract3.md)
builds for the \\(a,b)\\ block.

## Usage

``` r
d3_direction(spec, d3, params, npar, a, b, tv)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- d3:

  The third derivatives on the link scale.

- params, npar:

  The parameter names and block sizes.

- a, b:

  The block indices.

- tv:

  The predictors of the direction.

## Value

A numeric vector as long as the sample.
