# Refit One Bootstrap Replica

Re-estimates the coefficients of a fit on a simulated response, at its
hyperparameters or with them chosen again. See
[`predictive_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/predictive_bootstrap.md).

## Usage

``` r
bootstrap_refit(object, spec_b, refit, method, cfg, beta0)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec_b:

  Its specification with the simulated response.

- refit:

  `"coefficients"` or `"full"`.

- method:

  The fit's inner optimizer.

- cfg:

  Its
  [`inner_settings()`](https://statmodels7.github.io/statmodels7/reference/inner_settings.md).

- beta0:

  The fit's coefficients, the replica's start.

## Value

A list with `coefficients` and `hyper`.

## Examples

``` r
set.seed(4)
dd <- data.frame(x = runif(30))
dd$y <- 1 + dd$x + rnorm(30)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
sp <- S7::set_props(fit@spec, response = dd$y + rnorm(30, sd = 0.1))
statmodels7:::bootstrap_refit(fit, sp, "coefficients", fit@methods$smooth,
  statmodels7:::inner_settings(fit@methods$smooth),
  unlist(fit@coefficients, use.names = FALSE))$coefficients
#> $mu
#> [1] 1.6335080 0.3440787
#> 
#> $sigma
#> [1] -0.08989198
#> 
```
