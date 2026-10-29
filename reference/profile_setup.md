# The Objective, Its Hessian and the Scale of the Held Coefficients at a Point

Prepares what the polynomial of the trace needs at the coefficients
`b0`: the objective of the fit with its gradient, the factorized Hessian
of the coefficients that are not held, and the standard error of each
held coefficient in the profile, which scales the stencil step.

## Usage

``` r
profile_setup(ctx, b0)
```

## Arguments

- ctx:

  The context of the call.

- b0:

  The stacked coefficients of the point.

## Value

A list.
