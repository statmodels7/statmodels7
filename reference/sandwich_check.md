# Reject a Fit That the Sandwich Variance Does Not Cover

Signals an error for a fit with a structural term, a term whose block
moves with its coefficients.

## Usage

``` r
sandwich_check(fit, design)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

## Value

`NULL`, invisibly.

## Details

Prior weights are accepted; see
[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md)
for their meaning.
