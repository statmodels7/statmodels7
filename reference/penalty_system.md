# The Penalized System Around One Term

Collects what the graphs of a penalized block need: the stacked
coefficients, the positions of the term's coefficients, the information
of the log-likelihood over all the coefficients, the Hessian of the
penalties of the other terms, and the coordinates that are held out of a
Newton step.

## Usage

``` r
penalty_system(fit, param, name, cols, information)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param, name, cols:

  The equation, the term and the columns of the term that are drawn
  (`NULL` for all).

- information:

  An object from
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  or
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html).

## Value

A list with `beta`, `idx` (positions of the drawn coefficients), `names`
(their names), `grad` (the gradient of the log-likelihood), `info` (the
information), `S_other`, `units` (the penalized units of the term),
`free` (the coordinates of a Newton step), `spec`, `design` and
`labels`.
