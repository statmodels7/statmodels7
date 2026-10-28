# Graph of a Fitted Parameter Against a Covariate

Selects the graph of a parameter of the fitted distribution as a
function of one covariate, on the scale of the parameter, with a
pointwise confidence band. For a numeric covariate the parameter is
drawn as a curve over the range of the covariate; for a factor, as one
estimate with its interval per level.

## Usage

``` r
plot_param(
  param,
  over,
  at = NULL,
  level = 0.95,
  vcov_type = bayesian_vcov(),
  n_grid = 200
)
```

## Arguments

- param:

  The name of a parameter of the family, a single string.

- over:

  The name of the covariate on the horizontal axis.

- at:

  A named list of values for the other covariates, or `NULL`.

- level:

  The confidence level of the band.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance of the band.

- n_grid:

  The number of points of the grid of a numeric covariate.

## Value

A `PlotParam` object, inheriting from
[model_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The values are those of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
with `what = param` and `se = TRUE` on a grid of `over`. The other
covariates are held at `at`, or at their median (numeric) or first level
(factor), and the random effects are set to zero. The band is the
interval of the parameter that
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
returns for the variance selected by `vcov_type`, that is, the interval
of the predictor mapped through the inverse link.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_quantiles()`](https://statmodels7.github.io/statmodels7/reference/plot_quantiles.md),
[`plot_moments()`](https://statmodels7.github.io/statmodels7/reference/plot_moments.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- 1 + dd$x + rnorm(150, sd = exp(-1 + dd$x))
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_param(param = "sigma", over = "x"),
                  data = dd)
```
