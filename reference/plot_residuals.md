# Graph of the Residuals against a Quantity

Selects the graph of the quantile residuals against the predicted
conditional median (in ranks), a fitted parameter, the index of the
observation or a covariate, with a smoother and the horizontal lines at
0 and \\\pm 2\\.

## Usage

``` r
plot_residuals(against = "median", param = NULL, smooth = TRUE, seed = NULL)
```

## Arguments

- against:

  `"median"` (the default), `"fitted"`, `"index"`, or the name of a
  covariate.

- param:

  The parameter of `against = "fitted"`, otherwise `NULL`.

- smooth:

  Whether to draw the smoother.

- seed:

  A single number that seeds the randomization, or `NULL`.

## Value

A `PlotResiduals` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

With `against = "median"` the abscissa is the rank, among the
observations, of the conditional median \\F^{-1}(1/2; \hat\theta_i)\\;
the ranks spread the points evenly whatever the scale of the response,
and the median is defined for every family, including those without a
location parameter. With `against = "fitted"` the abscissa is the fitted
value of the parameter `param` on the parameter scale, and with
`"index"` the row number. Any other string names a covariate and
requires `data`. A numeric abscissa carries a loess smoother of the
residuals ([`loess`](https://rdrr.io/r/stats/loess.html)) when
`smooth = TRUE`, which should lie near zero under a correct model; a
factor is drawn with jittered points. For a discrete family the
residuals are randomized and `seed` makes them reproducible.

## See also

[`plot_qq()`](https://statmodels7.github.io/statmodels7/reference/plot_qq.md),
[`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md),
[`diagnostics()`](https://statmodels7.github.io/statmodels7/reference/diagnostics.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- 1 + dd$x + rnorm(150, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_residuals())

ggplot2::autoplot(fit, plot_type = plot_residuals(against = "x"), data = dd)
```
