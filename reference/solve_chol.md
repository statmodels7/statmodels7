# The Cholesky Solver

Selects a Cholesky decomposition of the assembled penalized information
\\X^\top W X + S\\ for the scoring step of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).
The information is formed explicitly, so the condition number of the
system is squared; a matrix that is not positive definite has its
eigenvalues floored, at the relative level `floor`, before it is
factored. The decomposition names no coordinate, so none is reported as
aliased. It is the route that
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
and
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
follow where the square root of the curvature does not exist.

[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
compares it with the other two solvers.

## Usage

``` r
solve_chol(floor = 1e-08)
```

## Arguments

- floor:

  The relative floor of the eigenvalues, a single positive number. A
  penalized matrix that is not positive definite is rebuilt from its
  eigendecomposition with every eigenvalue below `floor` times the
  largest absolute eigenvalue raised to that level. The default is
  `1e-8`.

## Value

An S7 object of class `SolveChol`, which inherits from
[IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md),
with the property `floor`.

## See also

[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
and
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
for the other solvers;
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
for the comparison;
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the function that takes the object.

## Examples

``` r
solve_chol()
#> Cholesky solver (floor 1e-08)
solve_chol(floor = 1e-10)
#> Cholesky solver (floor 1e-10)

set.seed(1)
dd <- data.frame(x = runif(120, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
G <- distributions7::gaussian1_distrib()
fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
               inner_optimizer = iwls(solver = solve_chol(floor = 1e-10)))
fit@methods$smooth
#> iwls: chol, expected information (the observed after 10 iterations)
#>   maxit 1000, tol 1e-06
```
