# Whether Scoring Pieces Carry a Positive Definite Curvature

`TRUE` where the pieces are square roots, which are positive
semidefinite by construction, or where the assembled matrix has its
smallest eigenvalue above a relative floor.
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
reads it to decide whether the observed pieces can take the step.

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

A single logical. A matrix with a non-finite entry answers `TRUE`, the
boundary being
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)'s
to hold.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
