# Graph of the Prior of a Random Effect

Selects the graph that compares the predicted effects of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term with the prior density that the fit estimated for them: a histogram
of the predicted effects on the density scale, with the density on top.
A term with several coefficients for each level gives one panel for each
coefficient.

## Usage

``` r
plot_prior()
```

## Value

A `PlotPrior` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The density is the marginal density of the coefficient under the
distribution of the effects: the gaussian density with the fitted
standard deviation for the default prior, and the density of the
distribution passed to `random()` (a Student t, for example) with its
fitted hyperparameters for a univariate prior. For a multivariate
gaussian prior the marginal is the gaussian density with the fitted
standard deviation of the coefficient. The predicted effects are shrunk
towards zero, so their histogram is narrower than the prior.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_caterpillar()`](https://statmodels7.github.io/statmodels7/reference/plot_caterpillar.md),
[`plot_random_cov()`](https://statmodels7.github.io/statmodels7/reference/plot_random_cov.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200), id = factor(sample(1:25, 200, TRUE)))
dd$y <- 1 + dd$x + rnorm(25)[dd$id] + rnorm(200, sd = 0.5)
fit <- statmod(y ~ x + random(~ 1 | id), distributions7::gaussian1_distrib(),
               dd)
ggplot2::autoplot(fit, term = "random(~1|id)", plot_type = plot_prior())
```
