# The Log-Density's Curvature Used by the Refresh Corrections

Returns \\\ell\_{ab}\\ on the link scale, from whichever route was used
to assemble the criterion's \\K\\.

## Usage

``` r
refresh_hessian(spec, design, coef, expected = FALSE, approx = "opg")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- expected:

  Whether the information is the expected one.

- approx:

  The approximation for the expected information.

## Value

A named list of numeric vectors.
