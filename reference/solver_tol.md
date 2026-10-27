# The Tolerances of a Solver

`solver_tol()` returns the `tol` of a
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
or
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
object and `NULL` for a
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
object, which has none. `solver_floor()` returns the `floor` of any of
the three.

## Usage

``` r
solver_tol(solver)

solver_floor(solver)
```

## Arguments

- solver:

  An object inheriting from
  [IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md).

## Value

A single number, or `NULL` for `solver_tol()` of a
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
object left at its default and of a
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
object.

## See also

[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
