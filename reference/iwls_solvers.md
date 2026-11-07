# Solvers Available to a Scoring Step

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
takes an argument `solver`, an object that selects the factorization
with which each scoring step is solved. Three objects are available,
built by
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md),
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md);
all three inherit from
[IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md).
A character string is not accepted. Each constructor has the numerical
tolerances of its factorization as arguments.

## Value

No value. The topic describes the three objects, of classes
[SolveQR](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md),
[SolveSVD](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md)
and
[SolveChol](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md),
that
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md),
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
return.

## Details

The step solves \\(X^\top W X + S)\\\delta = u\\, where \\X\\ is the
design, \\W\\ the diagonal matrix of the working weights (or, more
generally, the per-observation curvature), \\S\\ the Hessian of the
penalty, and \\u\\ the penalized score. The curvature is the matrix that
the argument `inner_information` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
selects (see
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)).
Let \\A = \[\sqrt{W}X;\\ \mathrm{chol}(S)\]\\, the design augmented by a
factor of the penalty, so that \\A^\top A = X^\top W X + S\\.

- [`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md):

  Factors \\A\\ by a QR decomposition, which does not form \\X^\top W
  X\\, so the condition number of the system is not squared. Dense
  matrices use a QR with column pivoting, which drops a column that is
  numerically dependent on the others (the coefficient is reported as
  aliased) and records the rank; the argument `tol` is the tolerance of
  that decision. Sparse matrices use a sparse QR, which declines a
  rank-deficient system, and the step is then solved by the dense route.
  This is the default.

- [`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md):

  Factors \\A\\ by a singular value decomposition and applies the
  pseudo-inverse. A singular value below the cut `tol` (by default
  \\\max(\mathrm{dim}(A))\\\epsilon\\\sigma\_{\max}\\, with \\\epsilon\\
  the machine precision) is discarded, so a deficient system gives a
  minimum-norm step and a numerical rank, and no coordinate is reported
  as aliased because a discarded direction is a combination of columns.
  It is the dearest of the three and has no sparse counterpart.

- [`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md):

  Factors the assembled information \\X^\top W X + S\\ by a Cholesky
  decomposition, after flooring the eigenvalues of a matrix that is not
  positive definite at the relative level `floor`. The matrix is formed
  explicitly, so the condition number is squared; the design of a
  break-point term is bounded only to \\\epsilon^{-1/2}\\, so the
  squaring can exhaust double precision.

[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
and
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
need a square root of the per-observation curvature and of the penalty.
Where either does not exist (an observed curvature far from the optimum,
a non-convex penalty), the step is solved from the assembled
information, as
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
does, and the eigenvalues of a matrix that is not positive definite are
floored at the `floor` of the solver. For this reason the three
constructors have the argument `floor`.

The tolerances and their defaults are:

|  |  |  |
|----|----|----|
| Argument | Meaning | Default |
| `tol` of [`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md) | pivoted QR rank tolerance | `1e-7` |
| `tol` of [`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md) | relative cut of the singular values | `NULL`, that is \\\max(\mathrm{dim}(A))\\\epsilon\\ |
| `floor` of the three | relative floor of the eigenvalues | `1e-8` |

The test of the sparse QR, which compares the equilibrated diagonal of
the triangular factor with \\\mathrm{ncol}(A)\\\epsilon\\, is not an
argument.

## See also

[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md),
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
for the constructors;
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the function that takes the object;
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
for the matrix that the solver factors.

## Examples

``` r
solve_qr()
#> QR solver (tol 1e-07, floor 1e-08)
solve_svd(tol = 1e-12)
#> SVD solver (tol 1e-12, floor 1e-08)
solve_chol(floor = 1e-10)
#> Cholesky solver (floor 1e-10)

# The three reach the same fit on a design of full rank.
set.seed(1)
dd <- data.frame(x = runif(120, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
G <- distributions7::gaussian1_distrib()
f_qr <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
                inner_optimizer = iwls(solver = solve_qr()),
                outer_criterion = NULL)
f_chol <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
                  inner_optimizer = iwls(solver = solve_chol()),
                  outer_criterion = NULL)
max(abs(unlist(coef(f_qr)) - unlist(coef(f_chol))))
#> [1] 5.551115e-17
```
