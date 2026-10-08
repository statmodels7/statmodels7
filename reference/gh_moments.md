# The First Two Moments of h^-1(eta) Under a Mixture of Gaussians

Each component is integrated by a 40-node Gauss-Hermite rule, or read at
its mean where its standard deviation is zero (a Monte Carlo draw).

## Usage

``` r
gh_moments(g, M, S, W)
```

## Arguments

- g:

  The link.

- M, S, W:

  As in
  [`group_sd_mixture()`](https://statmodels7.github.io/statmodels7/reference/group_sd_mixture.md).

## Value

A list with `m1` and `m2`, one value a row.
