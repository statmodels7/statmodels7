# Reject a Fit That the HAC Variance Does Not Cover

Signals an error for a fit with a term whose block moves with its
coefficients, and for a model with no structural term and no `time`.

## Usage

``` r
hac_check(fit, design, sw)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- sw:

  The
  [HacSandwich](https://statmodels7.github.io/statmodels7/reference/HacSandwich.md)
  object.

## Value

`NULL`, invisibly.
