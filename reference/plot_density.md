# Graph of Fitted Densities at Chosen Covariate Values

Selects the graph of the density, or the mass function, that the fit
gives at each row of a data frame of covariate values. With `over`, each
density is drawn vertically at its position on the horizontal axis, over
the observed data; without it, the densities are overlaid and coloured
by row.

## Usage

``` r
plot_density(at, over = NULL, n_grid = 200, scale = 0.9)
```

## Arguments

- at:

  A data frame of covariate values, one fitted distribution per row.

- over:

  `NULL`, or the name of a numeric covariate of `at` that gives the
  position of each density on the horizontal axis.

- n_grid:

  The number of points of the response grid of a continuous family.

- scale:

  With `over`, the width of the widest density as a fraction of the
  smallest gap between positions (of 15 per cent of the range of the
  covariate for a single row).

## Value

A `PlotDensity` object, inheriting from
[model_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

Each row of `at` is completed with the other covariates of the model,
held at their median (numeric) or first level (factor), and the random
effects are set to zero. The response range of each row runs from the
quantile 0.001 to the quantile 0.999 of its fitted distribution. For a
continuous family the density is evaluated on `n_grid` points of that
range. For a discrete family the mass is drawn at the integers from the
smallest to the largest of those quantiles over the rows, as bars or,
beside the data, as horizontal segments. The graph is not available for
a multivariate family.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_quantiles()`](https://statmodels7.github.io/statmodels7/reference/plot_quantiles.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- 1 + dd$x + rnorm(150, sd = exp(-1 + dd$x))
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, type = plot_density(at = data.frame(x = c(0.1, 0.5,
                                                                  0.9)),
                                           over = "x"), data = dd)
```
