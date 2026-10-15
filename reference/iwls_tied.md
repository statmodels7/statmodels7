# The Scoring Fit on the Subspace That a Kink Leaves Free

Runs
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
on \\\beta = \beta_0 + N w\\, with \\N\\ a basis of the directions that
the kink constraints and the held coordinates leave free
([`null_basis()`](https://statmodels7.github.io/statmodels7/reference/null_basis.md)),
and returns the point reached in the coefficients.

## Usage

``` r
iwls_tied(
  obj,
  beta,
  method,
  n,
  pieces_at,
  Ak,
  frozen,
  backup_at = NULL,
  newton_at = NULL,
  switch_at = NULL,
  switch_after = Inf,
  kinks_at = NULL,
  damp_on_reject = TRUE
)
```

## Arguments

- obj:

  The objective, with `fn` and `gr`.

- beta:

  The point \\\beta_0\\.

- method:

  The
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  method.

- n:

  The number of observations.

- pieces_at, backup_at, newton_at, switch_at, kinks_at:

  As
  [`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
  takes them, in the full coefficients.

- Ak:

  The rows of the kink constraints at `beta`.

- frozen:

  The positions held.

- switch_after, damp_on_reject:

  As
  [`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
  takes them.

## Value

`NULL` where the constraints leave no direction free, and otherwise a
list with `par`, the coefficients reached, `value`, `converged` and
`score`, the reduced fit's.

## Details

The reparametrization is linear, so the reduced problem's pieces are the
full problem's pieces carried through \\N\\: the square-root design and
the penalty's factor are multiplied by \\N\\ on the right, the assembled
matrix becomes \\N^\top A N\\, the gradient \\N^\top g\\ and the full
Hessian \\N^\top H N\\. A kink met in the reduced problem adds its row,
carried through \\N\\, in the same way.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
[`kink_constraints()`](https://statmodels7.github.io/statmodels7/reference/kink_constraints.md)
