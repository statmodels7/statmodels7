# The Criterion's Matrix With the Free Levels Appended

Assembles the matrix of
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
over the coefficients followed by one column of ones per free level
([`free_level_equations()`](https://statmodels7.github.io/statmodels7/reference/free_level_equations.md)),
factors it and returns the pieces that the criterion and its exact
gradient read.

## Usage

``` r
level_augmented(ctx, spec, design, coef, hyper, method, basis)
```

## Arguments

- ctx:

  An evaluation context from
  [`outer_context()`](https://statmodels7.github.io/statmodels7/reference/outer_context.md),
  or `NULL`.

- spec, design, coef, hyper, method, basis:

  As for
  [`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md).

## Value

`NULL` where no level is appended. Otherwise a list with `ok` (whether
the determinant exists), `logdet`, `q`, `K` (the pinned matrix,
coefficients first), `eq` (the equation of each appended column),
`npar`, `offs` and `total` (the coefficients' bookkeeping).

## Details

The matrix is assembled as the coefficients' own is, on the static
design: the information \\-X'WX\\ with the family's second derivatives
at the fitted predictors, plus the penalty's Hessian on the
coefficients. The level's column is a column of ones in its equation, so
its rows are the column sums of the weighted design blocks. The
coordinates that
[`laplace_pinned()`](https://statmodels7.github.io/statmodels7/reference/laplace_pinned.md)
names, and any with a non-finite curvature, are pinned as
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
pins them.

The result is cached in the evaluation context, so the criterion and its
gradient factor the matrix once at a point.

## See also

[`free_level_equations()`](https://statmodels7.github.io/statmodels7/reference/free_level_equations.md),
[`level_trace_pieces()`](https://statmodels7.github.io/statmodels7/reference/level_trace_pieces.md)
