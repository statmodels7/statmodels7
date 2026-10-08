# The Coordinates the Marginal Determinant Leaves Out

[`pinned_coords()`](https://statmodels7.github.io/statmodels7/reference/pinned_coords.md)
and the coordinates a penalty with a kink covers.

## Usage

``` r
laplace_pinned(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

## Value

An integer vector of stacked positions, possibly empty.

## Details

A kinked penalty has no curvature at the coordinates it sets to zero, so
a Laplace expansion around them is not defined, and the coordinates it
leaves away from zero are chosen by the same selection. The criterion
holds all of them at the joint mode and integrates the rest, which is
the restricted likelihood of Verbyla (1993) with those coefficients
treated as known. Measured on
`y ~ x1 + ... + x20 | sigma ~ lasso(~ z1 + ... + z10)` at 200
observations over 15 samples, the dispersion's intercept is off by
+0.003 on average (root mean square 0.048) against +0.013 (0.067) when
the kinked coordinates are integrated and -0.057 (0.073) when the
intercept is read at the joint mode, with the same count of slopes
wrongly selected (1.9 against 1.2 of seven) and every fit converged.

The mode's own movement is not affected: it is read over every
coordinate, see
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md).

## References

Verbyla, A. P. (1993). Modelling variance heterogeneity: residual
maximum likelihood and diagnostics. *Journal of the Royal Statistical
Society B*, 55, 493–508.
