# Graph of the Density of the Response under Each Regime

Selects the graph of the fitted density (the probability mass for a
discrete family) of the response under each regime of a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term, at the covariates of a reference observation, with one curve for
each regime.

## Usage

``` r
plot_regime_density(at = NULL, n_grid = 300)
```

## Arguments

- at:

  The row of the fitting data that supplies the covariates, or `NULL`
  for the observation with the median static predictor.

- n_grid:

  The number of points of the grid of the response of a continuous
  family.

## Value

A `PlotRegimeDensity` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

Under regime \\j\\ the predictor of the equation of the term is the
static predictor plus the level \\\mu_j\\ of the regime, and the other
parameters of the family are the fitted ones. The curves are the density
at those parameters for the reference observation `at`, which by default
is the one with the median static predictor. The horizontal range
extends from the 0.1 per cent quantile of the lowest regime to the 99.9
per cent quantile of the highest.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_states()`](https://statmodels7.github.io/statmodels7/reference/plot_states.md)

## Examples

``` r
set.seed(4)
n <- 80
s <- integer(n)
s[1] <- 1L
for (t in 2:n) s[t] <- if (runif(1) < 0.9) s[t - 1] else 3L - s[t - 1]
dd <- data.frame(t = seq_len(n), y = c(0, 3)[s] + rnorm(n, sd = 0.7))
fit <- statmod(y ~ regime(k = 2, time = t),
               distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, term = "regime(k = 2, time = t)",
                  plot_type = plot_regime_density())
```
