# Graph of a Filtered Parameter

Selects the graph of the path that a score-driven term
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html))
gives to the parameter of its equation: the filtered parameter over the
order of the series, with a pointwise confidence band, and optionally
the observed response beside it. A model in which several equations
carry a `gas()` term joins them in one filter, and the graph then has a
panel for each of the filtered equations.

## Usage

``` r
plot_filter(
  level = 0.95,
  scale = "response",
  response = FALSE,
  groups = NULL,
  vcov_type = bayesian_vcov()
)
```

## Arguments

- level:

  The confidence level of the band.

- scale:

  `"response"` (the default), the parameter, or `"link"`, the predictor.

- response:

  Whether to draw the observed response. It is drawn on the scale of the
  parameter, so it requires `scale = "response"`, a response that is a
  vector, and a term in one equation.

- groups:

  `NULL` for every group, or the labels of the groups to draw (the
  levels of the `by` variable as text).

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance of the band.

## Value

A `PlotFilter` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The path is the fitted parameter \\\theta_t = g^{-1}(\eta_t)\\, where
\\\eta_t\\ is the static predictor plus the level of the filter, and the
band is the interval of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
with `se = TRUE`, which differentiates the filter in the coefficients
and in the parameters of the term (the standard error of a filtered
predictor). The horizontal axis is the variable named by `time` in the
term, or the row order when it has none.

For a panel (`gas(by = id)`) the graph has one panel for each group, at
most twelve; with more groups, `groups` names the ones to draw.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_news_impact()`](https://statmodels7.github.io/statmodels7/reference/plot_news_impact.md)

## Examples

``` r
set.seed(4)
n <- 100
x <- rnorm(n)
f <- 0
y <- integer(n)
for (t in seq_len(n)) {
  lam <- exp(0.7 + 0.3 * x[t] + f)
  y[t] <- rpois(1, lam)
  f <- 0.15 * (y[t] - lam) * lam^(-0.5) + 0.8 * f
}
dd <- data.frame(t = seq_len(n), x = x, y = y)
fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t, scaling = 0.5),
               distributions7::poisson_distrib(), dd)
ggplot2::autoplot(fit, term = "gas(p = 1, q = 1, time = t, scaling = 0.5)",
                  type = plot_filter(response = TRUE))
```
