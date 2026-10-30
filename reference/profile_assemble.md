# Assemble the Table of a Profile

Converts the points returned by
[`profile_compute()`](https://statmodels7.github.io/statmodels7/reference/profile_compute.md)
into the data frame of
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
and the row of the unrestricted fit.

## Usage

``` r
profile_assemble(ctx, grid, pts, method)
```

## Arguments

- ctx:

  The context of the call.

- grid:

  The data frame of held values.

- pts:

  The points.

- method:

  The profile method.

## Value

A list with `points`, the data frame, and `unrestricted_point`, the row
of the fit.

## Details

The objective of a point is the negative of the penalized
log-likelihood, \\-\ell + \tfrac12 P\_\lambda\\, so that the penalty
column is \\-(\mathrm{objective} + \ell)\\ and the penalized
log-likelihood is the negative of the objective.
