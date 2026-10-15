# The Coordinates That the Marginal Determinant Leaves Out

[`pinned_coords()`](https://statmodels7.github.io/statmodels7/reference/pinned_coords.md)
and the coordinates that a penalty with a kink covers.

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

A kinked penalty has no curvature at the coordinates that it sets to
zero, so a Laplace expansion around them is not defined, and the
coordinates it leaves away from zero are chosen by the same selection.
The criterion holds all of them at the joint mode and integrates the
rest, which is the restricted likelihood of Verbyla (1993) with those
coefficients treated as known.

The mode's own movement is not affected: it is read over every
coordinate, see
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md).

## References

Verbyla, A. P. (1993). Modelling variance heterogeneity: residual
maximum likelihood and diagnostics. *Journal of the Royal Statistical
Society B*, 55, 493–508.
