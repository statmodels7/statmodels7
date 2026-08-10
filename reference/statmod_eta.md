# The Linear Predictors and the Parameters They Give

Turns a coefficient structure into the per-observation parameters the
distribution's generics take.

## Usage

``` r
statmod_eta(spec, design, coef)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, as
  [`statmod_design`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- coef:

  A named list of coefficient vectors, one per parameter.

## Value

A list with `eta` and `theta`, both named lists.
