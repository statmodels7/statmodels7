# Resolve a Parameter Structure

Returns the fitted coefficients when `par` is `NULL`, and validates a
supplied structure against the design otherwise.

## Usage

``` r
par_at(fit, par, design)
```

## Arguments

- fit:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- par:

  A named list or `NULL`.

- design:

  The design to validate against.

## Value

A named list of coefficient vectors.
