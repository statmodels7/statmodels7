# The Coefficients of a Random Term by Level

Returns the coefficients of a random term as a matrix with one row for
each level of the grouping factor and one column for each coefficient
within the level, with their positions in the stacked coefficients of
the model.

## Usage

``` r
random_block(fit, param, name, tm)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param, name:

  The equation and the name of the term.

- tm:

  The
  [modelterms7::RandomTerm](https://statmodels7.github.io/modelterms7/reference/RandomTerm.html).

## Value

A list with `levels`, `wnames`, `idx` (positions within the equation's
coefficients, group by group), `coef_names` (their names in
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md))
and `b`, the matrix.
