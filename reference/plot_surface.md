# Graph of the Surface of a Smooth of Two Variables

Selects the graph of the effect of a smooth of two numeric variables (a
tensor product
[`modelterms7::te()`](https://statmodels7.github.io/modelterms7/reference/te.html)
or `s(x1, x2)`) as a raster with contour lines over the ranges of the
two variables. With `slices`, the graph instead shows curves of the
effect against the first variable at fixed values of the second, each
with a pointwise band.

## Usage

``` r
plot_surface(
  mask = TRUE,
  slices = NULL,
  level = 0.95,
  vcov_type = bayesian_vcov(),
  n_grid = 60
)
```

## Arguments

- mask:

  Whether to mark the cells in which the interval of the effect contains
  zero.

- slices:

  `NULL`, or a numeric vector of values of the second variable at which
  curves of the effect against the first variable are drawn.

- level:

  The confidence level.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance of the effect.

- n_grid:

  The number of points of the grid of each variable.

## Value

A `PlotSurface` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The effect is \\f(x_1, x_2) = X(x_1, x_2) b\\, evaluated on a regular
grid of `n_grid` points for each variable, with a standard error
\\\sqrt{X V X'}\\ from the block of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the coefficients of the term. The raster is colored on a scale that
is symmetric about zero. With `mask = TRUE`, a point is drawn over every
cell in which \\\|f\| \< z\_{(1+\text{level})/2}\\\mathrm{se}\\, that
is, where the pointwise interval of the effect contains zero. The slices
are the same function at the values `slices` of the second variable, and
the bands are \\f \pm z\\\mathrm{se}\\. A smooth with a `by` factor is
drawn in one panel for each level.

## See also

[`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md),
[`plot_derivative()`](https://statmodels7.github.io/statmodels7/reference/plot_derivative.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x1 = runif(200), x2 = runif(200))
dd$y <- sin(2 * pi * dd$x1) * dd$x2 + rnorm(200, sd = 0.2)
fit <- statmod(y ~ te(x1, x2), distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, term = "te(x1, x2)", type = plot_surface(n_grid = 30),
                  data = dd)

ggplot2::autoplot(fit, term = "te(x1, x2)",
                  type = plot_surface(slices = c(0.2, 0.8)), data = dd)
```
