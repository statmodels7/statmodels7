# The Scale of Each Label of a Profile

The standard error of each label of `which`: for a coefficient the one
of the profile, from
[`profile_setup()`](https://statmodels7.github.io/statmodels7/reference/profile_setup.md),
and for a hyperparameter the one of the exact Hessian of the criterion
on the scale of its link, or 1 where it is not available. It scales the
step of the stencil.

## Usage

``` r
profile_scales(ctx)
```

## Arguments

- ctx:

  The context of the call.

## Value

A numeric vector.
