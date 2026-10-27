# The Log-Density's Curvature Used by the Refresh Corrections

Returns \\\ell\_{ab}\\ on the link scale, from whichever route was used
to assemble the criterion's \\K\\.

## Usage

``` r
refresh_hessian(spec, design, coef, expected = FALSE, emethod = NULL)
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

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

## Value

A named list of numeric vectors.
