# Histogram of the Probability Integral Transform

Selects the histogram of the values \\F(y_i; \hat\theta_i)\\ of the
distribution function at the observations, on the density scale, with
the uniform density and a pointwise binomial band for each bin.

## Usage

``` r
plot_pit(bins = 20, seed = NULL, level = 0.95, legend = TRUE)
```

## Arguments

- bins:

  The number of bins.

- seed:

  A single number that seeds the randomization, or `NULL`.

- level:

  The level of the band.

- legend:

  Whether to draw the legend that names the elements of the graph.

## Value

A `PlotPit` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

Under a correct model the values are uniform on \\(0, 1)\\. Each of the
`bins` classes of equal width receives a binomial number of the \\n\\
values, with probability \\1/\text{bins}\\, and the band of level
`level` is the interval of that binomial law divided by \\n\\ times the
bin width. For a discrete family the value is randomized between
\\F(y_i^-)\\ and \\F(y_i)\\ (the values are \\\Phi\\ of the quantile
residuals of
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)),
and `seed` makes it reproducible. A U shape indicates a response more
dispersed than the model, a hump one less dispersed, and a slope a bias
in the location.

## See also

[`plot_calibration()`](https://statmodels7.github.io/statmodels7/reference/plot_calibration.md),
[`plot_qq()`](https://statmodels7.github.io/statmodels7/reference/plot_qq.md),
[`diagnostics()`](https://statmodels7.github.io/statmodels7/reference/diagnostics.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- rpois(150, exp(1 + dd$x))
fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_pit(bins = 10, seed = 1))
```
