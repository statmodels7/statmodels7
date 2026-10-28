# Graph of the Covariance of a Random Effect

Selects the graph of the predicted effects of two coefficients of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term, such as an intercept and a slope, one point for each level, with
the ellipse of the fitted covariance of the effects.

## Usage

``` r
plot_random_cov(level = 0.95, which = c(1, 2))
```

## Arguments

- level:

  The probability content of the ellipse.

- which:

  The two coefficients within a level, by position.

## Value

A `PlotRandomCov` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The ellipse is the set \\\\b : b' \Sigma^{-1} b = c\\\\, centred at
zero, where \\\Sigma\\ is the fitted covariance matrix of the two
effects, built from their fitted standard deviation and correlation, and
\\c\\ is the quantile `level` of the chi-squared distribution with two
degrees of freedom. It holds the fraction `level` of the effects when
the prior is gaussian, which is the only prior for which the graph is
defined. The predicted effects are shrunk, so the points lie inside the
ellipse more often than that fraction.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_caterpillar()`](https://statmodels7.github.io/statmodels7/reference/plot_caterpillar.md),
[`plot_prior()`](https://statmodels7.github.io/statmodels7/reference/plot_prior.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(300), id = factor(sample(1:20, 300, TRUE)))
dd$y <- 1 + dd$x + rnorm(20)[dd$id] + 0.5 * rnorm(20)[dd$id] * dd$x +
  rnorm(300, sd = 0.5)
fit <- statmod(y ~ x + random(~ 1 + x | id), distributions7::gaussian1_distrib(),
               dd)
ggplot2::autoplot(fit, term = "random(~1+x|id)", plot_type = plot_random_cov())
```
