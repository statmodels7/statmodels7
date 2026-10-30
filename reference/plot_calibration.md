# Calibration of the Predicted Quantiles

Selects the graph of the observed coverage of the predicted quantiles
against their nominal level, with the diagonal and a pointwise band.
With `over`, one graph is drawn for each class of a covariate.

## Usage

``` r
plot_calibration(
  probs = seq(0.05, 0.95, by = 0.05),
  over = NULL,
  n_inter = 4,
  level = 0.95,
  seed = NULL,
  legend = TRUE
)
```

## Arguments

- probs:

  The nominal levels, numbers in \\(0, 1)\\.

- over:

  The name of a covariate, or `NULL`.

- n_inter:

  The number of classes of a numeric covariate.

- level:

  The level of the band.

- seed:

  A single number that seeds the randomization, or `NULL`.

- legend:

  Whether to draw the legend that names the elements of the graph.

## Value

A `PlotCalibration` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

For a probability \\p\\ the predicted quantile of observation \\i\\ is
\\F^{-1}(p; \hat\theta_i)\\, and the observed coverage is the fraction
of observations at or below it, \\\hat c(p) = n^{-1} \sum_i 1\\u_i \le
p\\\\, where \\u_i = F(y_i; \hat\theta_i)\\ is the probability integral
transform, randomized for a discrete family (the values are \\\Phi\\ of
the quantile residuals). Under a correct model \\\hat c(p)\\ has mean
\\p\\ and variance \\p(1 - p)/n\\, and the band is \\p \pm z\_{(1 +
\text{level})/2}\sqrt{p(1 - p)/n}\\, limited to \\\[0, 1\]\\. With
`over`, the classes are formed as in
[`plot_worm()`](https://statmodels7.github.io/statmodels7/reference/plot_worm.md)
and require `data`.

## See also

[`plot_pit()`](https://statmodels7.github.io/statmodels7/reference/plot_pit.md),
[`plot_worm()`](https://statmodels7.github.io/statmodels7/reference/plot_worm.md),
[`diagnostics()`](https://statmodels7.github.io/statmodels7/reference/diagnostics.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- 1 + dd$x + rnorm(150, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_calibration())

ggplot2::autoplot(fit, plot_type = plot_calibration(over = "x", n_inter = 3),
                  data = dd)
```
