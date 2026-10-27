# The Predictive Mixture by Monte Carlo

Returns `n_draw` point components. Each component is a draw of the
estimation error from \\\mathrm{N}(0, C)\\ at every row and a draw of a
new group's effects for every term in `heavy`. See
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md).

## Usage

``` r
prior_mc(spec, eta, base, heavy, weight, n_draw)
```

## Arguments

- spec:

  The specification at the rows predicted.

- eta:

  A named list of the predictors' means.

- base:

  The covariance of the estimation error plus the Gaussian effects, an
  array `P x P x n`.

- heavy:

  The Student t and other entries of
  [`prior_parts()`](https://statmodels7.github.io/statmodels7/reference/prior_parts.md).

- weight:

  The total weight.

- n_draw:

  The number of draws.

## Value

A list of `n_draw` components, each with `C = NULL`.

## Examples

``` r
set.seed(8)
gg <- data.frame(g = factor(rep(1:30, each = 5)), x = rnorm(150))
gg$y <- 1 + gg$x + (rexp(30) * sample(c(-1, 1), 30, TRUE))[gg$g] +
  rnorm(150, sd = 0.4)
lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
fit <- statmod(y ~ x + random(~ 1 | g, distrib = lp),
               distributions7::gaussian1_distrib(), gg)
predict(fit, what = "response", newdata = data.frame(x = 0, g = "new"),
        random = "zero", interval = prediction_interval())
#>        fit       se      lower    upper
#> 1 1.199763 1.009002 -0.8944612 3.293986
```
