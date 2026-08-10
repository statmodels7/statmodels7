# The Information of the Weighted Log-Likelihood

The negative Hessian in the coefficients, assembled block by block from
the distribution's own link-scale second derivatives.

## Usage

``` r
statmod_information_at(
  spec,
  coef,
  design = statmod_design(spec),
  expected = TRUE,
  approx = "bartlett"
)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

- expected:

  Whether to use the expected information.

- approx:

  The approximation, when the expected one is not closed.

## Value

A square matrix over the stacked coefficients.

## Details

`expected = TRUE` gives the expected information, which is what Fisher
scoring inverts; `approx` is passed to distributions7 and read only
where the family has no closed expected information.
