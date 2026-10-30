# The Interval of a Profile That Is Constant Between Observations

The Interval of a Profile That Is Constant Between Observations

## Usage

``` r
profile_step_limits(d, covariate, cut)
```

## Arguments

- d:

  The profile, as
  [`profile_curve()`](https://statmodels7.github.io/statmodels7/reference/profile_curve.md)
  returns it.

- covariate:

  The values of the covariate.

- cut:

  The quantile that the statistic must not exceed.

## Value

The two limits: the value of the covariate below the first interval in
which the statistic does not exceed `cut`, and the value above the last
one. A limit is `NA` where that interval is the first or the last of the
profile, the region possibly extending beyond it.
