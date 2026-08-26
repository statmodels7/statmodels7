# The Size of a Penalty's Kink

Half the jump of \\\rho'\\ across the first point the penalty is not
differentiable at, which is the half-width of the subdifferential there.

## Usage

``` r
kink_scale(pen, theta, eps = 1e-04)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters, as a list or a named vector.

- eps:

  How far either side of the kink to read the derivative.

## Value

A single non-negative number, `0` where there is no kink.

## Details

A coefficient stays at the kink while the unpenalized score at that
point is inside the subdifferential, so this number is what a score has
to exceed for a coefficient to leave zero. It is \\\lambda\\ for the
lasso, SCAD and MCP, and \\\lambda\alpha\\ for the elastic net, but it
is read and never assumed: a penalty built from a density carries the
same information in a hyperparameter of its own, and a Laplace prior
written by its scale has a kink whose size falls as that scale grows.

The penalty is separable over the coefficients wherever it has a kink, a
kinked penalty under a general map being the generalized-lasso problem
and rejected upstream, so one coordinate answers for all of them.

The derivative just past the kink carries the distance it was read at,
MCP's being \\\lambda - \epsilon/\gamma\\, so the two readings are
extrapolated to the limit. Without that the shape parameters appear to
move the kink by a millionth of themselves, which is enough to be
selected for a path, and measures \\\epsilon\\ in place of the penalty.
The extrapolation is exact for a penalty that is affine on each side of
the kink, which the lasso, SCAD, MCP and the elastic net all are.

## See also

[`kink_hypers()`](https://statmodels7.github.io/statmodels7/reference/kink_hypers.md),
[`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
