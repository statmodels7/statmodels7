# The Factor of a Diagonal Penalty

One row per coordinate the penalty reaches, carrying the square root of
that coordinate's entry, which is
[`penalty_sqrt`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)'s
answer where the matrix is diagonal.

## Usage

``` r
diagonal_sqrt(S, p)
```

## Arguments

- S:

  A diagonal penalty Hessian.

- p:

  Its dimension.

## Value

A matrix with `p` columns, or `NULL`.

## Details

The thresholds are the eigen route's, read on the diagonal, which for a
diagonal matrix IS its spectrum: a negative entry beyond the tolerance
makes the penalty indefinite and there is no factor to return, and an
entry at the tolerance is a null direction and contributes no row. The
class of the result mirrors the argument's rather than being chosen:
[`augmented_solve`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
routes on whether either of the two factors is sparse, so returning a
sparse factor for a dense design would send a dense fit through the
sparse route and a dense one through neither.
