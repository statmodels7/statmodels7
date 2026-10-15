# How Far Above Its Mode a Filter Fit Stopped

The Newton decrement \\\tfrac12 g^\top H^{-1} g\\ of the joint objective
over the coefficients and the free parameters of a score-driven term, at
the point that an inner fit returned.

## Usage

``` r
filter_decrement(spec, design, res, hyper)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design, whose structural state holds the filter's parameters at
  the point.

- res:

  The inner fit's result, a list with `par` (the stacked coefficients)
  and `obj` (the objective).

- hyper:

  The hyperparameters.

## Value

A single non-negative number, or `Inf`.

## Details

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
reads it to decide whether the fit from the term's own start failed,
beside the optimizer's flag: where the filter explodes the objective is
flat to rounding, the flag at that point depends on the platform, and
the decrement is many orders of magnitude above
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).
The coefficients that the specification holds (those that the outer
criterion estimates) are left out, their score not being a residual of
the inner fit. A gradient or a curvature that is not finite, or a
curvature that cannot be solved, gives `Inf`.

## See also

[`statmod_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/statmod_joint_pieces.md),
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
