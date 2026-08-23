# Run the Structural Terms at the Current Parameters

The level each filter adds to the equation it sits in, together with the
derivative of that level in the term's own parameters.

## Usage

``` r
statmod_filter_at(spec, design, eta_static, theta_static)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- eta_static:

  The static predictors, one per distribution parameter.

- theta_static:

  The parameters they give.

## Value

A named list, one entry per structural term, carrying its predictor, the
derivative of that predictor in the term's own parameters, and the
curvature the recursion read at each step.

## Details

The result is memoized on the coefficients and the term's parameters,
since the objective, its gradient and its curvature are asked for at the
same point in turn and a filter is the expensive part of each.

## See also

[`statmod_structural`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
