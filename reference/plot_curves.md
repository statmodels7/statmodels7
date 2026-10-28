# Graph of the Curves of a Random Smooth

Selects the graph of a random smooth, `random(~ s(x) | g)`: one thin
curve for each level of the grouping factor, over the range of the
variable within that level. The graph reads the covariates and requires
`data`.

## Usage

``` r
plot_curves(n_grid = 100, population = FALSE)
```

## Arguments

- n_grid:

  The number of points of each curve.

- population:

  Whether to add the population smooth of the variable.

## Value

A `PlotCurves` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The curve of a level is its contribution to the predictor of the
equation, \\g_j(x) = a\_{0j} + N(x) a_j + Z(x) \delta_j\\, the departure
of the level from the curve that the other terms of the equation
describe. It is centred around zero and does not include the population
curve. With `population = TRUE` the population curve, taken from the
single smooth of the same variable in the same equation, is added to
each curve and drawn as a thick line, so that the curves are the curves
of the levels.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_caterpillar()`](https://statmodels7.github.io/statmodels7/reference/plot_caterpillar.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(240), id = factor(rep(1:8, each = 30)))
dd$y <- sin(2 * pi * dd$x) + rnorm(8, sd = 0.4)[dd$id] * dd$x +
  rnorm(240, sd = 0.3)
fit <- statmod(y ~ s(x) + random(~ s(x) | id),
               distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, term = "random(~s(x)|id)", plot_type = plot_curves(),
                  data = dd)
```
