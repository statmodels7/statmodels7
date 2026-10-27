# Graph of the Predicted Effects of a Random Effect

Selects the graph with one row for each level of the grouping factor of
a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term: the predicted effect with its confidence interval and, optionally,
the unshrunk estimate of the level as an open circle. A term with
several coefficients for each level, such as a random slope, gives one
panel for each coefficient.

## Usage

``` r
plot_caterpillar(
  raw = TRUE,
  sort = TRUE,
  level = 0.95,
  vcov_type = bayesian_vcov(),
  information = distributions7::expected_info()
)
```

## Arguments

- raw:

  Whether to draw the unshrunk estimates.

- sort:

  Whether to order the levels by the predicted effect of the first
  coefficient.

- level:

  The confidence level of the intervals.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance of the intervals.

- information:

  An object from
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  (the default) or
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
  the weight of the working residuals.

## Value

A `PlotCaterpillar` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The predicted effect is the coefficient of the level, and the interval
is \\\hat b \pm z\_{(1+\text{level})/2}\\\sqrt{V\_{jj}}\\, where \\V\\
is the variance matrix of the coefficients selected by `vcov_type`. The
unshrunk estimate of a level is the least squares estimate of the
coefficients of that level from the partial residuals \\Z_i \hat b_j +
(z_i - \eta_i)\\, over the observations of the level, where \\Z_i\\ is
the row of the design of the term, \\\hat b_j\\ the predicted effect of
the level, and \\z_i - \eta_i\\ the working residual of the equation
(see
[`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md)).
With one coefficient for each level it is the mean of the partial
residuals of the level, and for a gaussian random intercept it is the
mean of the residuals of the level plus its predicted effect. A level
with too few observations to estimate all its coefficients has no open
circle.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_prior()`](https://statmodels7.github.io/statmodels7/reference/plot_prior.md),
[`plot_random_cov()`](https://statmodels7.github.io/statmodels7/reference/plot_random_cov.md),
[`plot_shrinkage()`](https://statmodels7.github.io/statmodels7/reference/plot_shrinkage.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(120), id = factor(sample(1:12, 120, TRUE)))
dd$y <- 1 + dd$x + rnorm(12)[dd$id] + rnorm(120, sd = 0.5)
fit <- statmod(y ~ x + random(~ 1 | id), distributions7::gaussian1_distrib(),
               dd)
ggplot2::autoplot(fit, term = "random(~1|id)", type = plot_caterpillar())
```
