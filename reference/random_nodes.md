# Nodes and Weights for the Effects a Marginal Prediction Integrates

A product Gauss-Hermite grid where every prior is Gaussian and the total
dimension is at most three; the prior itself, for
[`marginal_quad()`](https://statmodels7.github.io/statmodels7/reference/marginal_quad.md),
where the only prior is univariate and not Gaussian; Monte Carlo draws
otherwise.

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
dimension, so at most 8000 in all. Measured on a random intercept, 20
nodes give the marginal mean exactly to 1e-15 under a log link and to
3e-10 and 9e-06 under a logit at prior standard deviations of 1 and 2,
where 10000 Monte Carlo draws stop between 1e-2 and 1e-3. The draws are
taken from a seed of their own and the caller's random stream is put
back afterwards, so the prediction is reproducible and moves nothing a
session draws next.
