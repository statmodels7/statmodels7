# Graph of the Scores of an Equation

Selects the graph of the contribution of each observation to the score
of one equation, on the scale of the link, against the index of the
observation or a covariate. The five largest absolute contributions are
highlighted and labeled with the row number.

## Usage

``` r
plot_score(param, against = "index", legend = TRUE)
```

## Arguments

- param:

  The equation, a parameter of the family.

- against:

  `"index"` (the default) or the name of a covariate.

- legend:

  Whether to draw the legend that names the elements of the graph.

## Value

A `PlotScore` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The contribution of observation \\i\\ to the score of the equation
`param` is \\s_i = \omega_i\\ \partial \ell_i / \partial \eta_i\\, with
\\\ell_i\\ the log-likelihood of the observation, \\\eta_i\\ the
predictor of the equation and \\\omega_i\\ the prior weight. The
contributions sum, against the design, to the gradient of the fit, so a
large \\\|s_i\|\\ marks an observation that strains that equation. They
are the numerator of the working residuals of
[`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md).
A covariate as abscissa requires `data`.

## See also

[`plot_influence()`](https://statmodels7.github.io/statmodels7/reference/plot_influence.md),
[`plot_residuals()`](https://statmodels7.github.io/statmodels7/reference/plot_residuals.md),
[`diagnostics()`](https://statmodels7.github.io/statmodels7/reference/diagnostics.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(120))
dd$y <- 1 + dd$x + rnorm(120, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_score(param = "mu"))

ggplot2::autoplot(fit, plot_type = plot_score(param = "sigma", against = "x"),
                  data = dd)
```
