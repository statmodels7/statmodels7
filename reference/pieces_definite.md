# Whether Scoring Pieces Carry a Positive Definite Curvature

Returns `TRUE` where the pieces are square roots, which are positive
semidefinite by construction, or where the assembled matrix has its
smallest eigenvalue above a relative floor.
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
uses the result to check whether the observed pieces can take the step.

## Usage

``` r
pieces_definite(pc, rel = 1e-08)
```

## Arguments

- pc:

  The pieces, as
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  returns them.

- rel:

  The relative floor,
  [`pd_repair()`](https://statmodels7.github.io/statmodels7/reference/pd_repair.md)'s
  own.

## Value

A single logical. The result is `TRUE` for a matrix with a non-finite
entry; holding the boundary is left to
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md).

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
