# The Specification with the Criterion's Coefficients Held

Returns the specification in which the coefficients that the outer
criterion estimates beside the hyperparameters (the dispersion of the
marginal criteria) are held at their fitted values, together with the
description of those coordinates.

## Usage

``` r
criterion_held_spec(fit, design)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design.

## Value

A list with `spec` and `gam`.
