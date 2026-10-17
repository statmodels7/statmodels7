# The Size of a Penalty's Kink

Returns half the jump of \\\rho'\\ across the first point at which the
penalty is not differentiable, which is the half-width of the
subdifferential there.

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

  The larger of the two distances from the kink at which the derivative
  is read.

## Value

A single non-negative number, `0` where there is no kink.

## Details

A coefficient stays at the kink while the unpenalized score at that
point is inside the subdifferential, so a score has to exceed this
number for a coefficient to leave zero. It is \\\lambda\\ for the lasso,
SCAD and MCP, and \\\lambda\alpha\\ for the elastic net, but it is read
and never assumed: a penalty built from a density carries the same
information in a hyperparameter of its own, and a Laplace prior written
by its scale has a kink whose size falls as that scale grows.

A kinked penalty is separable over the coefficients, and one under a
general map (the generalized-lasso problem) is rejected upstream. Under
a diagonal map, which `standardize = TRUE` writes, the kink differs by
coordinate. This function reads the first coordinate's kink, and
[`coord_kinks()`](https://statmodels7.github.io/statmodels7/reference/coord_kinks.md)
returns the kink of every coordinate.

The derivative just past the kink depends on the distance at which it is
read (MCP's, for example, is \\\lambda - \epsilon/\gamma\\), so the
one-sided derivatives at distances \\\epsilon/2\\ and \\\epsilon\\ are
combined by Richardson extrapolation to the limit. The extrapolation is
exact for a penalty that is affine on each side of the kink, as the
lasso, SCAD, MCP and the elastic net are.

## See also

[`kink_hypers()`](https://statmodels7.github.io/statmodels7/reference/kink_hypers.md),
[`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
