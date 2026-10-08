# The Predictive Mixture Over a New Group's Effects

The components
[`predictive_response()`](https://statmodels7.github.io/statmodels7/reference/predictive_response.md)
reads for one set of estimates: the predictors' mean `eta` and the
covariance `C0` of their estimation error, widened by the effects of a
new group for every term `random` sets aside.

## Usage

``` r
predictive_mixture(
  fit,
  spec,
  design,
  aside,
  eta,
  C0,
  weight = 1,
  n_draw = 5000L,
  n_nodes = 1000L
)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  whose hyperparameters give the priors.

- spec, design:

  The specification and design at the rows predicted.

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  that are not `"conditional"`.

- eta:

  A named list of the predictors' means.

- C0:

  The covariance of their estimation error, an array `P x P x n`, or
  `NULL` for none.

- weight:

  The total weight of the components returned.

- n_draw:

  The number of Monte Carlo draws, where they are needed.

- n_nodes:

  The number of quantile nodes for one univariate prior that is neither
  Gaussian nor a Student t
  ([`group_quantile_nodes()`](https://statmodels7.github.io/statmodels7/reference/group_quantile_nodes.md)).

## Value

A list of components, each a list with `eta`, `C` (an array or `NULL`)
and `weight`, with the attribute `infinite_variance`.

## Details

A Gaussian prior adds \\z^\top\Sigma_b z\\ to the covariance, so one
component carries it. A Student t prior, univariate on one coefficient a
group or multivariate, is a scale mixture of Gaussians, \$\$b \mid w
\sim \mathrm{N}(0, \Sigma/w), \qquad w \sim \mathrm{Gamma}(\nu/2,
\nu/2),\$\$ so it gives one component per node of a trapezoidal rule in
\\\log w\\
([`gamma_nodes()`](https://statmodels7.github.io/statmodels7/reference/gamma_nodes.md)),
each with the Gaussian covariance \\z^\top\Sigma z/w\\ and the node's
weight; two such terms give the product of the two rules. One univariate
prior of another family (a logistic, a Laplace, a Cauchy) gives one
component per quantile node, \\b_k = Q((k - 1/2)/K)\\ with \\K = 1000\\
([`group_quantile_nodes()`](https://statmodels7.github.io/statmodels7/reference/group_quantile_nodes.md)),
so the result does not depend on the random seed: measured on a logistic
prior the ends of the interval are 3.5e-6 from the exact ones, against
0.12 for 5000 draws, and on a Cauchy 1e-11 against 0.61. Any other prior
(several such terms, a prior over several coordinates, a prior shared by
a label that is not Gaussian), or more than two Student t terms, is
averaged by Monte Carlo: `n_draw` draws of a new group's effects, from
[`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html)
or from the scale mixture, each with a draw of the estimation error, so
the interval depends on the random seed. The quantiles of the mixture
exist whatever the prior; its standard deviation does not under a
Student t with \\\nu \le 2\\ or a univariate prior whose family reports
no finite variance (a Cauchy), and the result carries the attribute
`infinite_variance` to say so.

## See also

[`predictive_response()`](https://statmodels7.github.io/statmodels7/reference/predictive_response.md),
which reads the components,
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
the caller.

## Examples

``` r
set.seed(8)
gg <- data.frame(g = factor(rep(1:12, each = 5)), x = rnorm(60))
gg$y <- 1 + gg$x + 0.8 * rt(12, df = 3)[gg$g] + rnorm(60, sd = 0.4)
tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
fit <- statmod(y ~ x + random(~ 1 | g, distrib = tp),
               distributions7::gaussian1_distrib(), gg)
predict(fit, "response", data.frame(x = 0, g = "new"), random = "zero",
        interval = "prediction")
#>        fit       se     lower    upper
#> 1 1.435761 2.249694 -1.642905 4.514427
```
