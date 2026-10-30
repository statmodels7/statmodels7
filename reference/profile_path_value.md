# The Criterion of a Path at One Point of a Profile

Scores the refit at one point of a profile with the criterion that chose
the hyperparameter of a kinked penalty along its path,
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
as the path scores each of its points: the curvature of a scaled SCAD or
MCP is read at the coefficients of the point, and the degrees of freedom
are counted on the active coordinates.

## Usage

``` r
profile_path_value(ctx, spec, design, par, cf, hyper)
```

## Arguments

- ctx:

  The context of the call.

- spec, design:

  The specification of the point and its design.

- par, cf:

  The stacked coefficients of the refit and their list.

- hyper:

  The hyperparameters of the point.

## Value

A single number, `NA` where the profile holds no hyperparameter of a
path, where the path's criterion is
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md), or
where the criterion cannot be read.
