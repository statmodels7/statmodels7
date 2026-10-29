# The Exact Penalized Log-Likelihood at Predicted Coefficients

The Exact Penalized Log-Likelihood at Predicted Coefficients

## Usage

``` r
profile_exact(ctx, S, b, val)
```

## Arguments

- ctx:

  The context of the call.

- S:

  The list that
  [`profile_setup()`](https://statmodels7.github.io/statmodels7/reference/profile_setup.md)
  returns.

- b:

  The predicted stacked coefficients.

- val:

  The held values.

## Value

A point, as
[`profile_point()`](https://statmodels7.github.io/statmodels7/reference/profile_point.md)
returns it, with `converged` `NA`.
