# The Pieces One Scoring Step Needs

Assembles, at the current coefficients, whichever of the square-root
design, the penalty's factor and the penalized information the requested
decomposition will actually use.

## Usage

``` r
iwls_pieces(spec, design, coef, hyper, method)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, refreshed at `coef` if any term needs it.

- coef:

  A named list of coefficient vectors, one per distribution parameter.

- hyper:

  The hyperparameters, per penalized term.

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object, read for `decomposition`, `hessian` and `approx`.

## Value

A list of three, of which the ones the route does not need are `NULL`:

- `R`:

  the square-root design, as
  [`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
  returns it, or `NULL`.

- `C`:

  the penalty's factor, as
  [`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
  returns it, or `NULL`.

- `A`:

  the assembled penalized information, or `NULL`.

`R` is also `NULL` when the per-observation curvature had no Cholesky
factor, and
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
falls back to the assembled route there.

## Details

Only what the decomposition will read is built. `"qr"` and `"svd"` need
the square-root design `R` and the penalty's factor `C` and never form
the information; `"chol"` and `"chol_crossprod"` need the assembled `A`
and neither square root. Building all three unconditionally would pay
for the route not taken at every iteration of the loop.

## See also

[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md),
which consumes this,
[`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
and
[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
for the two factors.
