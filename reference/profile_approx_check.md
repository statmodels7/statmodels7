# Check That a Profile Can Be Approximated

Rejects, with a pointer to
[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md),
a call of an approximate method with `hyper = "refit"` or on a model
whose refitted coefficients are not smooth functions of the held ones.

## Usage

``` r
profile_approx_check(ctx, method)
```

## Arguments

- ctx:

  The context of the call.

- method:

  The name of the method, for the message.

## Value

`NULL`, invisibly.
