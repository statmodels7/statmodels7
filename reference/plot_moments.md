# Graph of the Moments of the Fitted Distribution Against a Covariate

Selects the graph of moments of the fitted distribution as functions of
one covariate, one panel for each moment with its own vertical scale.

## Usage

``` r
plot_moments(over, moments = c("mean", "variance"), at = NULL, n_grid = 200)
```

## Arguments

- over:

  The name of the covariate on the horizontal axis.

- moments:

  A character vector of moments.

- at:

  A named list of values for the other covariates, or `NULL`.

- n_grid:

  The number of points of the grid of a numeric covariate.

## Value

A `PlotMoments` object, inheriting from
[model_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The moments are those of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
with the random effects set to zero, which evaluates the closed forms of
the family or, for a family without them, a quadrature or a series. The
available moments are `"mean"`, `"variance"`, `"std_dev"`, `"skewness"`
and `"kurtosis"`. A moment that does not exist for some covariate value
(the mean of a Cauchy, the variance of a Student t with two degrees of
freedom or fewer) is returned as `NaN` or `Inf` by the family and is
left out of the graph. The other covariates are held at `at`, or at
their median (numeric) or first level (factor). The graph is not
available for a multivariate family.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_param()`](https://statmodels7.github.io/statmodels7/reference/plot_param.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- rgamma(150, shape = 3, rate = 3 / exp(1 + dd$x))
fit <- statmod(y ~ x, distributions7::gamma1_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_moments(over = "x",
                  moments = c("mean", "skewness")), data = dd)
```
