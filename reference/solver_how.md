# The Name of a Solver

Internal generic that returns the word with which the numerical routines
of the step
([`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md),
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md))
identify a solver: `"qr"`, `"svd"` or `"chol"`. It dispatches on the
class of `solver`.

## Usage

``` r
solver_how(solver)
```

## Arguments

- solver:

  An object inheriting from
  [IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md).

## Value

A single string.

## See also

[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
