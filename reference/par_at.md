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

  The design against which the lengths are validated.

## Value

A named list of coefficient vectors, one per distribution parameter in
the family's order.

## Details

Three things are checked: that `par` is a named list, that every name of
`par` is a distribution parameter, and that the vector given for each is
as long as that equation's design is wide. Each failure signals an
error. The message for a wrong name and the message for a wrong length
both give the name of the parameter, so that a short vector is not
recycled into a different model. A distribution parameter that `par`
omits keeps its fitted coefficients.

## See also

[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md),
whose three methods all pass through this.
