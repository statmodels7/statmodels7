# Nodes and Weights for the Effects That a Marginal Prediction Integrates

A product Gauss-Hermite grid where every prior is Gaussian and the total
dimension is at most three; the prior itself, for
[`marginal_quad()`](https://statmodels7.github.io/statmodels7/reference/marginal_quad.md),
where the only prior is univariate, not Gaussian and has a parent
family; Monte Carlo draws otherwise.

## Usage

``` r
random_nodes(priors, nsim = 10000L)
```

## Arguments

- priors:

  A list of
  [`random_prior()`](https://statmodels7.github.io/statmodels7/reference/random_prior.md)
  results.

- nsim:

  The number of draws where the grid is not used.

## Value

A list with `b`, one matrix per prior with a row per node and a column
per coordinate, and `w`, the weights, summing to one. For one univariate
prior that is not Gaussian, a list with `quad`, that prior.

## Details

The grid has \\\min(40, \lfloor 8000^{1/D}\rfloor)\\ nodes per
dimension, so at most 8000 in all. Where the grid is not used, `nsim`
draws with equal weights are taken for each prior, a Gaussian prior
through its Cholesky factor and any other through
[`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html).
The draws are taken from a seed of their own and the caller's random
stream is put back afterwards, so the prediction is reproducible and
leaves the session's random stream unchanged.
