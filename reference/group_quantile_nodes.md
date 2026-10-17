# Quantile Nodes for One Univariate Prior That Is Not Gaussian

Where one univariate prior of a family other than the Student t is the
only prior that is not Gaussian among those set aside for a new group's
interval, the mixture over its effect is built from the prior's own
quantiles, \\b_k = Q((k - 1/2)/K)\\ with \\K = 1000\\ and equal weights,
and each node is a Gaussian component that carries the estimation error
and the Gaussian effects. The result does not depend on the random seed,
whereas the Monte Carlo draws of
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md)
do.

## Usage

``` r
group_quantile_nodes(spec, eta, base, parts, K = 1000L)
```

## Arguments

- spec:

  The specification at the rows predicted.

- eta:

  A named list of the predictors' means.

- base:

  The covariance of the estimation error plus the Gaussian effects, an
  array `P x P x n`.

- parts:

  What
  [`prior_parts()`](https://statmodels7.github.io/statmodels7/reference/prior_parts.md)
  returns.

- K:

  The number of nodes.

## Value

A list of components in the shape that
[`predictive_response()`](https://statmodels7.github.io/statmodels7/reference/predictive_response.md)
reads, or `NULL` where the priors are of another kind.
