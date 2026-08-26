# A Square-Root Factor of the Penalty

Returns \\C\\ with \\C'C = S(\theta)\\, the block-diagonal penalty
Hessian, dropping the rows a null space contributes nothing to.

## Usage

``` r
penalty_sqrt(S)
```

## Arguments

- S:

  The penalty Hessian, a `p x p` symmetric matrix, dense or sparse.
  [`Matrix::isDiagonal()`](https://rdrr.io/pkg/Matrix/man/isTriangular-methods.html)
  decides which route is taken, so a matrix stored as diagonal and one
  merely having zero off-diagonals are both caught.

## Value

A matrix with `ncol(S)` columns and one row per retained eigendirection,
so at most `p` and fewer where the penalty has a null space. Its class
mirrors `S`'s, dense for dense and sparse for sparse. `NULL` when the
penalty is indefinite beyond the tolerance.

## Why an eigendecomposition

A penalty is positive semidefinite and is often rank deficient. A
spline's is deficient by exactly the dimension of its null space, which
is the whole point of the construction: the null space is what the
penalty does not shrink. A Cholesky fails there, so the factor comes
from an eigendecomposition with the non-positive eigenvalues dropped,
and the rows they would have contributed are simply absent.

A non-convex penalty has an indefinite Hessian and no such factor at
all. `NULL` is returned and the caller falls back to the assembled
route.

## The diagonal shortcut

A diagonal penalty is factored by taking the square root of its
diagonal. That is the same answer the eigendecomposition gives, since
the eigenvalues of a diagonal matrix are its diagonal entries, so the
two routes agree by construction and a test pins them together.

The case is worth detecting for how often it arises. A ridge is
diagonal, a random effect is diagonal, and the Demmler-Reinsch penalty
of
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html)
is \\\mathrm{diag}(0, 1, \ldots, 1)\\ exactly. So is any block-diagonal
assembly of them. The factor is recomputed at every iteration of the
scoring loop, so its cost is multiplied by the iteration count: measured
on a random intercept over 1000 groups, one dense eigendecomposition of
the 1003 by 1003 penalty cost 0.63 s and was 83 per cent of the whole
fit.

## See also

[`diagonal_sqrt()`](https://statmodels7.github.io/statmodels7/reference/diagonal_sqrt.md)
for the diagonal route,
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
for the solve this feeds,
[`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
for the other half of the augmented matrix.
