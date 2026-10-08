# Quantile Nodes for One Univariate Prior That Is Not Gaussian

Where the only prior a new group's interval sets aside that is not
Gaussian is one univariate prior of a family other than the Student t,
the mixture over its effect is built from the prior's own quantiles,
\\b_k = Q((k - 1/2)/K)\\ with \\K = 1000\\ and equal weights, each node
a Gaussian component carrying the estimation error and the Gaussian
effects. The result does not depend on the random seed, where the Monte
Carlo draws of
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md)
give the ends of a group interval with a standard deviation of 0.030 to
0.037 over 20 seeds on a logistic prior.

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

A list of components in the shape
[`predictive_response()`](https://statmodels7.github.io/statmodels7/reference/predictive_response.md)
reads, or `NULL` where the priors are of another kind.
