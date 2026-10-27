# The Contribution of Part of a Term and Its Standard Error

Evaluates \\X b\\ for the block of a term on `newdata` (or on the
fitting rows when `newdata` is `NULL`), restricted to the columns
`cols`, with the standard error \\\sqrt{x_i' V x_i}\\ of each row.

## Usage

``` r
term_contribution(fit, param, name, tm, cols, newdata, vcov_type)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param, name, cols:

  As in
  [`term_coef_block()`](https://statmodels7.github.io/statmodels7/reference/term_coef_block.md).

- tm:

  The term.

- newdata:

  A data frame, or `NULL`.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

## Value

A list with `fit` and `se`, numeric vectors.
