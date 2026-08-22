# Pin the Coordinates a Boundary Has Frozen

Replaces the row and column of every coordinate whose curvature is not
finite by the identity's, leaving the rest of the matrix alone.

## Usage

``` r
pin_boundary(K)
```

## Arguments

- K:

  The penalized information.

## Value

`K` with the frozen coordinates pinned, or `K` unchanged where there are
none.

## Details

A parameter that has reached the clamp its link keeps it strictly inside
leaves the family's curvature there at `NaN`, and one such entry is
enough to deny the whole matrix a Cholesky factor. The coordinate is not
one the criterion integrates over – it is held at a boundary, where its
own score is exactly zero – so what the Laplace approximation wants from
it is nothing, which is what a unit diagonal contributes: \\\log\lvert
K\rvert\\ gains \\\log 1 = 0\\, and \\K^{-1}g\\ returns that
coordinate's own score, which is zero.

The shape is preserved deliberately. Dropping the coordinate would
change the dimension of `K` and of its inverse, which some twenty
consumers index into by position; pinning says the same thing and breaks
none of them.

## See also

[`ctx_penalized`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md),
[`iwls_solve`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
