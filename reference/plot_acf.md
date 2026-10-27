# Autocorrelation of the Residuals

Selects the graph of the sample autocorrelation of the quantile
residuals, or of their squares, in the order of the observations, with
the limits \\\pm 1.96/\sqrt{n}\\.

## Usage

``` r
plot_acf(lag_max = 30, squared = FALSE, seed = NULL)
```

## Arguments

- lag_max:

  The largest lag.

- squared:

  Whether to use the squared residuals.

- seed:

  A single number that seeds the randomization, or `NULL`.

## Value

A `PlotAcf` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The autocorrelations are those of
[`stats::acf()`](https://rdrr.io/r/stats/acf.html) at the lags 1 to
`lag_max`. Under a correct model with independent observations the
quantile residuals are independent standard normal variables, the
autocorrelations have standard deviation near \\1/\sqrt{n}\\, and about
five percent of them fall outside the limits. The autocorrelation of the
squares (`squared = TRUE`) detects a variance that depends on the past.
The order is the row order of the data, so the graph is meaningful for a
series, or for data sorted within a panel; for a discrete family the
residuals are randomized and `seed` makes them reproducible.

## See also

[`plot_residuals()`](https://statmodels7.github.io/statmodels7/reference/plot_residuals.md),
[`diagnostics()`](https://statmodels7.github.io/statmodels7/reference/diagnostics.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:150)
dd$y <- as.numeric(arima.sim(list(ar = 0.5), n = 150))
fit <- statmod(y ~ 1, distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, type = plot_acf(lag_max = 20))
```
