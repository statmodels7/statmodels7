# Graph of the Curves of a Random Smooth

Selects the graph of a random smooth, `random(~ s(x) | g)`: one thin
curve for each level of the grouping factor, over the range of the
variable within that level, and optionally the population curve as a
thick line. The graph reads the covariates and requires `data`.

## Usage

``` r
plot_curves(n_grid = 100, population = FALSE, legend = TRUE)
```

## Arguments

- n_grid:

  The number of points of each curve.

- population:

  Whether to draw the population curve \\s(x)\\ as well; it requires a
  smooth of the variable in the equation.

- legend:

  Whether to draw the legend that names the elements of the graph.

## Value

A `PlotCurves` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The contribution of a level to the predictor of the equation is \\g_j(x)
= a\_{0j} + N(x) a_j + Z(x) \delta_j\\, its departure from the
population curve \\s(x)\\, which is the single smooth of the same
variable in the same equation. The curve of a level is \\s(x) +
g_j(x)\\. When the equation has no smooth of the variable, the curve of
a level is \\g_j(x)\\ alone.

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
