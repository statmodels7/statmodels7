# The Factor of a Diagonal Penalty

Returns one row per coordinate that the penalty reaches, and each row
carries the square root of that coordinate's entry. This is the value of
[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
where the matrix is diagonal.

## Usage

``` r
diagonal_sqrt(S, p)
```

## Arguments

- S:

  A diagonal penalty Hessian, `p x p`, dense or sparse. Only its
  diagonal is read.

- p:

  Its dimension. The caller supplies it, so that the function uses the
  caller's own count and does not infer a dimension.

## Value

A matrix with `p` columns and one row per retained coordinate, in which
each row is zero except for the square root of that coordinate's entry.
Its class mirrors `S`'s. `NULL` when any entry is negative beyond the
tolerance.

## Details

The thresholds are the eigen route's, applied to the diagonal, which for
a diagonal matrix is its spectrum. A negative entry beyond the tolerance
makes the penalty indefinite and there is no factor to return; an entry
at or below the tolerance is a null direction and contributes no row.

The class of the result follows the argument's: a sparse `S` gives a
sparse factor and a dense `S` a dense one.
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
takes the sparse route whenever either of its two factors is sparse, so
the class is not changed here.

## See also

[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md),
which dispatches here for a diagonal penalty.
