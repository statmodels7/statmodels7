# Graph of the Conditional Quantiles Against a Covariate

Selects the graph of the quantile curves of the fitted distribution as
functions of one covariate, with the observed points behind them. The
median is drawn thicker and the outermost curves dashed.

## Usage

``` r
plot_quantiles(
  over,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  at = NULL,
  points = TRUE,
  n_grid = 200,
  legend = TRUE
)
```

## Arguments

- over:

  The name of the covariate on the horizontal axis.

- probs:

  The probabilities of the curves, numbers in \\(0, 1)\\.

- at:

  A named list of values for the other covariates, or `NULL`.

- points:

  Whether to draw the observed values of the response against `over`.

- n_grid:

  The number of points of the grid of a numeric covariate.

- legend:

  Whether to draw the legend that names the elements of the graph.

## Value

A `PlotQuantiles` object, inheriting from
[model_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The curve of probability \\p\\ is \\q_p(x) = F^{-1}(p; \theta(x))\\,
where \\\theta(x)\\ is the vector of fitted parameters at the covariate
value \\x\\, from
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
with the random effects set to zero, and \\F^{-1}\\ is
[`distributions7::distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.html).
The other covariates are held at `at`, or at their median (numeric) or
first level (factor). For a discrete family the curves are step
functions, since the quantile takes integer values. The graph is not
available for a multivariate family.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_param()`](https://statmodels7.github.io/statmodels7/reference/plot_param.md),
[`plot_density()`](https://statmodels7.github.io/statmodels7/reference/plot_density.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- 1 + dd$x + rnorm(150, sd = exp(-1 + dd$x))
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_quantiles(over = "x"), data = dd)
```
