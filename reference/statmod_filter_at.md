# Run the Structural Terms at the Current Parameters

Returns the level that each filter adds to the equation in which it
sits, together with the derivative of that level in the term's own
parameters.

## Usage

``` r
statmod_filter_at(spec, design, eta_static, theta_static)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- eta_static:

  The static predictors, one per distribution parameter.

- theta_static:

  The parameters that they give.

## Value

A named list, one entry per structural term, carrying its predictor, the
derivative of that predictor in the term's own parameters, and the
curvature that the recursion read at each step.

## Details

The function keeps no memo of its own. The caller,
[`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md),
keeps one keyed on the coefficients and the term's parameters, because
the objective, its gradient and its curvature are evaluated at the same
point in turn and a filter is the expensive part of each.

## See also

[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
