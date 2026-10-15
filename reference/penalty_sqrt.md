# A Square-Root Factor of the Penalty

Returns \\C\\ with \\C'C = S(\theta)\\, the block-diagonal penalty
Hessian, omitting the rows that correspond to its null space.

## Usage

``` r
penalty_sqrt(S)
```

## Arguments

- S:

  The penalty Hessian, a `p x p` symmetric matrix, dense or sparse.
  [`Matrix::isDiagonal()`](https://rdrr.io/pkg/Matrix/man/isTriangular-methods.html)
  decides which route is taken, so a matrix stored as diagonal and one
  that only has zero off-diagonal entries are both treated as diagonal.

## Value

A matrix with `ncol(S)` columns and one row per retained eigendirection,
so at most `p` rows and fewer where the penalty has a null space. The
result is a Matrix sparse matrix only where `S` is a diagonal matrix
stored sparse; for any other `S`, including the zero matrix, it is a
dense base matrix. `NULL` when the penalty is indefinite beyond the
tolerance.

## Why an eigendecomposition

A penalty is positive semidefinite and is often rank deficient. A
spline's is deficient by exactly the dimension of its null space, the
set of coefficient directions that the penalty does not shrink. A
Cholesky factorization fails there, so the factor comes from an
eigendecomposition with the non-positive eigenvalues dropped, and the
rows that they would have contributed are absent. A zero matrix gives a
factor with no rows.

A non-convex penalty has an indefinite Hessian and no such factor.
`NULL` is returned and the caller falls back to the assembled route.

## The diagonal shortcut

A diagonal penalty is factored by taking the square root of its
diagonal. That is the same answer the eigendecomposition gives, since
the eigenvalues of a diagonal matrix are its diagonal entries, so the
two routes agree by construction and a test pins them together.

The case is common. A ridge penalty and the penalty of a random effect
with independent coefficients are diagonal, and so is the penalty of a
smoother in the Demmler-Reinsch basis (see
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html)),
which is zero on the null-space coordinates and a multiple of the
smoothing parameter elsewhere. Any block-diagonal assembly of such
penalties is diagonal as well. The factor is recomputed at every
iteration of the scoring loop, so the diagonal route avoids one dense
eigendecomposition per iteration.

## See also

[`diagonal_sqrt()`](https://statmodels7.github.io/statmodels7/reference/diagonal_sqrt.md)
for the diagonal route,
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
for the solve this feeds,
[`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
for the other half of the augmented matrix.
