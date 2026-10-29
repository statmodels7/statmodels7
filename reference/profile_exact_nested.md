# The Exact Quantities of a Point Predicted by the Nested Expansion

Computes the penalized log-likelihood, its two parts and, where the
profile reads it, the outer criterion, exactly at a point that
[`profile_nested_predict()`](https://statmodels7.github.io/statmodels7/reference/profile_nested_predict.md)
returns.

## Usage

``` r
profile_exact_nested(ctx, pr, val)
```

## Arguments

- ctx:

  The context of the call.

- pr:

  The predicted point.

- val:

  The held values.

## Value

A point, as
[`profile_point()`](https://statmodels7.github.io/statmodels7/reference/profile_point.md)
returns it, with `converged` `NA`.
