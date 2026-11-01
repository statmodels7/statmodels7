# A Block of the Jeffreys Penalty

Describes one block of coefficients on which the Jeffreys penalty
\\a\log\|\tilde I\|\\ is computed by
[`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md)
and its derivatives.

## Usage

``` r
jfp_block(
  param,
  cols,
  power = 0.5,
  const = FALSE,
  kind = c("static", "refresh", "filter"),
  term = NULL,
  slot = NULL,
  basis = NULL,
  C = NULL
)
```

## Arguments

- param:

  The name of the distribution parameter whose predictor the block
  enters.

- cols:

  Integer positions of the block's columns in that equation's design.

- power:

  The power \\a \> 0\\ of the determinant.

- const:

  Whether the constant direction is part of the block, as decided by
  [`jfp_constant()`](https://statmodels7.github.io/statmodels7/reference/jfp_constant.md).

- kind:

  One of `"static"`, `"refresh"` or `"filter"`.

- term, slot:

  For the kinds `"refresh"` and `"filter"`, the term and the parameter
  of the term whose subformula holds the block.

- basis, C:

  For the kind `"refresh"`: the columns of the equation's design from
  which the rows are combined, and the matrix that combines them, so
  that the rows are the refreshed design's columns `basis` times `C`.
  For the kind `"filter"`: positions among the filter's parameters and
  the matrix that combines the derivatives along them. The constant
  direction of the subformula's level, where it is included, is the
  first column of `C`.

## Value

A list of class `jfp_block` with the arguments as fields.

## Details

The block's rows \\\tilde V\_{q,t}\\ are the derivative of the predictor
of equation \\q\\ at observation \\t\\ with respect to the block's
coefficients, preceded by the constant direction of its formula level
when `const` is `TRUE`. For a block of `kind = "static"` the rows are
the block's columns of the equation's design, which do not depend on the
coefficients, and the constant direction is the column of ones.

The other kinds describe a block inside the subformula of a term:
`"refresh"` for a term whose block is recomputed at the coefficients
(`nl()`, a smoothed `seg()`), and `"filter"` for a parameter of `gas()`.
Their rows move with the coefficients. For the kind `"refresh"` they are
the refreshed design's columns `basis` times `C`; for the kind
`"filter"` they are the derivatives of the filtered predictors along the
columns of `C` placed at the positions `basis` among the filter's
parameters
([`jfp_filter_dirs()`](https://statmodels7.github.io/statmodels7/reference/jfp_filter_dirs.md)).

## See also

[`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md),
[`jfp_constant()`](https://statmodels7.github.io/statmodels7/reference/jfp_constant.md)
