# The SVD Solver

Selects a singular value decomposition of the augmented design \\A =
\[\sqrt{W}X;\\ \mathrm{chol}(S)\]\\ for the scoring step of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).
The step applies the pseudo-inverse of \\A^\top A\\: a singular value of
\\A\\ below the cut `tol` times the largest singular value is discarded.
A rank-deficient system gives a minimum-norm step and a numerical rank,
where a Cholesky decomposition would fail. No coordinate is reported as
aliased, because a discarded direction is a combination of columns. The
solver has no sparse counterpart, and it costs more than
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md).

[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
compares it with the other two solvers.

## Usage

``` r
solve_svd(tol = NULL, floor = 1e-08)
```

## Arguments

- tol:

  The relative cut of the singular values: `NULL` or a single positive
  number. With a number, a singular value \\\sigma\\ is discarded when
  \\\sigma \le \mathrm{tol}\\\sigma\_{\max}\\. With `NULL`, the default,
  the cut is \\\max(\mathrm{dim}(A))\\\epsilon\\, with \\\epsilon\\ the
  machine precision (`.Machine$double.eps`) and \\\mathrm{dim}(A)\\ the
  dimensions of the augmented design, which depend on the model.

- floor:

  The relative floor of the eigenvalues, as in
  [`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md).
  It is applied when the step is solved from the assembled matrix
  because a square root of the curvature or of the penalty does not
  exist. The default is `1e-8`.

## Value

An S7 object of class `SolveSVD`, which inherits from
[IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md),
with the properties `tol` and `floor`.

## See also

[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
for the other solvers;
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
for the comparison;
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the function that takes the object.

## Examples

``` r
solve_svd()
#> SVD solver (tol default, floor 1e-08)
solve_svd(tol = 1e-12)
#> SVD solver (tol 1e-12, floor 1e-08)

set.seed(1)
dd <- data.frame(x = runif(120, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
G <- distributions7::gaussian1_distrib()
fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
               inner_optimizer = iwls(solver = solve_svd(tol = 1e-12)))
fit@methods$smooth
#> iwls: svd, expected information (the observed after 10 iterations)
#>   maxit 1000, tol 1e-06
```
