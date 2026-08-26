# Resolve a Parameter Structure

Returns the fitted coefficients when `par` is `NULL`, and otherwise
checks a supplied structure against the design before returning it.

## Usage

``` r
par_at(fit, par, design)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md),
  read for its fitted coefficients.

- par:

  A named list of coefficient vectors, or `NULL`.

- design:

  The design to validate the lengths against.

## Value

A named list of coefficient vectors, one per distribution parameter in
the family's order.

## Details

Three things are checked: that `par` is a named list, that it names
every distribution parameter, and that each vector is as long as that
equation's design is wide. Each failure signals an error naming the
parameter, which is the alternative to R recycling a short vector into a
different model without a word.

## See also

[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md),
whose three methods all pass through this.
