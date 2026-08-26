# The Factor of a Diagonal Penalty

One row per coordinate the penalty reaches, carrying the square root of
that coordinate's entry, which is
[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)'s
answer where the matrix is diagonal.

## Usage

``` r
diagonal_sqrt(S, p)
```

## Arguments

- S:

  A diagonal penalty Hessian, `p x p`, dense or sparse. Only its
  diagonal is read.

- p:

  Its dimension. Passed in so that the caller's own count is used and no
  dimension is inferred here.

## Value

A matrix with `p` columns and one row per retained coordinate, each row
zero except for the square root of that coordinate's entry. Its class
mirrors `S`'s. `NULL` when any entry is negative beyond the tolerance.

## Details

The thresholds are the eigen route's, applied to the diagonal, which for
a diagonal matrix is its spectrum. A negative entry beyond the tolerance
makes the penalty indefinite and there is no factor to return; an entry
at or below the tolerance is a null direction and contributes no row.

The class of the result mirrors the argument's and is not chosen here.
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
routes on whether either of its two factors is sparse, so a sparse
factor returned for a dense design would send a dense fit down the
sparse route.

## See also

[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md),
which dispatches here for a diagonal penalty.
