# The QR Solver

Selects a QR decomposition of the augmented design \\A = \[\sqrt{W}X;\\
\mathrm{chol}(S)\]\\ for the scoring step of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).
The cross-product \\X^\top W X\\ is never formed, so the condition
number of the system is not squared. A dense matrix is factored with
column pivoting: a column that is numerically dependent on the others is
dropped from the step, the coefficient is reported as aliased, and the
rank of the system is recorded. A sparse matrix is factored by a sparse
QR, which declines a rank-deficient system and leaves the step to the
dense route.

It is the default of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).
Choose it for ordinary designs and for a break-point term, whose design
is ill-conditioned.
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
compares it with the other two solvers.

## Usage

``` r
solve_qr(tol = 1e-07, floor = 1e-08)
```

## Arguments

- tol:

  The tolerance of the pivoted QR (the argument `tol` of
  [`base::qr()`](https://rdrr.io/r/base/qr.html)), a single positive
  number. A column is declared dependent on the previous ones when its
  norm after the projection is under `tol` times its original norm. The
  same value is the threshold of the equilibrated diagonal under which
  the threaded triangular factor, used for large designs with several
  threads, defers to the pivoted route. The default is `1e-7`.

- floor:

  The relative floor of the eigenvalues, a single positive number,
  applied when the step is solved from the assembled matrix because a
  square root of the curvature or of the penalty does not exist. The
  eigenvalues of a matrix that is not positive definite are raised to
  `floor` times the largest absolute eigenvalue before the
  factorization. The default is `1e-8`.

## Value

An S7 object of class `SolveQR`, which inherits from
[IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md),
with the properties `tol` and `floor`.

## See also

[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
for the other solvers;
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
for the comparison;
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the function that takes the object.

## Examples

``` r
solve_qr()
#> QR solver (tol 1e-07, floor 1e-08)
solve_qr(tol = 1e-9)
#> QR solver (tol 1e-09, floor 1e-08)

set.seed(1)
dd <- data.frame(x = runif(120, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
G <- distributions7::gaussian1_distrib()
fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
               inner_optimizer = iwls(solver = solve_qr(tol = 1e-9)))
fit@methods$smooth
#> iwls: qr, expected information (the observed after 10 iterations)
#>   maxit 1000, tol 1e-06
```
