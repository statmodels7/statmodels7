# The Trace Pieces of the Gradient on the Matrix With the Free Levels

Inverts the matrix of
[`level_augmented()`](https://statmodels7.github.io/statmodels7/reference/level_augmented.md)
and returns what
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
traces against: the inverse over the coefficients, and the
per-observation diagonals of the inverse on the design with the level
columns appended.

## Usage

``` r
level_trace_pieces(lev, spec, design, coef)
```

## Arguments

- lev:

  The result of
  [`level_augmented()`](https://statmodels7.github.io/statmodels7/reference/level_augmented.md).

- spec, design, coef:

  As for
  [`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md).

## Value

A list with `Mb` (the inverse over the coefficients), `M`, `G`, `npar`
and `offs` (the inverse, the diagonals and the bookkeeping on the
augmented design), or `NULL` where the matrix cannot be inverted.

## Details

The determinant moves with the mode through the family's third
derivative contracted with the per-observation diagonals \\x_i' M x_i\\
of the augmented design, which is
[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)'s
formula with that design; the level's column is constant, so it has no
derivative of its own. The penalty does not cover the level, so the
traces against the penalty's derivatives read the coefficients' block of
the inverse.

## See also

[`level_augmented()`](https://statmodels7.github.io/statmodels7/reference/level_augmented.md),
[`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
