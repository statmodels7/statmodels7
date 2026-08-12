# The Terms as the Fit Left Them

The specification with every refreshable term replaced by the one the
fit arrived at, so that a break-point, a nonlinear parameter and the
block they imply are read off the fitted object rather than off the
specification it started from.

## Usage

``` r
statmod_fitted_spec(spec, coef, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

## Value

A
[`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## See also

[`statmod_commit_refresh`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
