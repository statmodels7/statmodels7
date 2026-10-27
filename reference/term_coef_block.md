# The Coefficients and Variance of Part of a Term

Returns the coefficients of a term in its equation, or of the columns
`cols` within it, with their names in
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
and the block of the variance matrix that they span.

## Usage

``` r
term_coef_block(fit, param, name, cols, vcov_type)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param:

  The equation.

- name:

  The name of the term.

- cols:

  Positions within the term's coefficients, or `NULL`.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

## Value

A list with `coef` (a numeric vector) and `vcov` (a matrix).
