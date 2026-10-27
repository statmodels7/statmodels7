# Graph of the Coefficient Paths of a Penalized Block

Selects the graph of the coefficients of a penalized block against the
logarithm of the penalty parameter \\\lambda\\, with the fitted value of
\\\lambda\\ marked by a dashed line. For a kinked penalty (a lasso, an
elastic net, a SCAD or an MCP) a second axis counts the nonzero
coefficients.

## Usage

``` r
plot_path(log_lambda = NULL, information = distributions7::expected_info())
```

## Arguments

- log_lambda:

  `NULL`, or an increasing numeric vector of values of \\\log\lambda\\.

- information:

  An object from
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  (the default) or
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
  the information of the log-likelihood in the closed form of a
  quadratic penalty.

## Value

A `PlotPath` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

For a kinked penalty the path is computed by refitting the model at each
value of \\\lambda\\, from the largest to the smallest, each fit
starting from the previous one. The other hyperparameters of the penalty
and all the coefficients of the other equations are held at their fitted
values, so that the path at the fitted \\\lambda\\ is the fit. By
default the grid runs from the smallest \\\lambda\\ that empties the
block, found by doubling the fitted value, down to a tenth of the fitted
value.

For a penalty that is quadratic in the coefficients, \\\lambda S\\ with
\\S\\ fixed, the path is the closed form \$\$b(\lambda) = (I + S_o +
\lambda S)^{-1} (I + S_o + \hat\lambda S)\\\hat b,\$\$ where \\\hat b\\
and \\\hat\lambda\\ are the fitted coefficients and penalty parameter,
\\I\\ is the information of the log-likelihood over all the coefficients
and \\S_o\\ the Hessian of the penalties of the other terms. It is the
exact mode for a gaussian mean with the identity link, with the other
equations held, and the first-order approximation obtained by replacing
the log-likelihood with its quadratic expansion at \\\hat b\\ otherwise.
By default \\\log\lambda\\ runs over five units on each side of the
fitted value.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_shrinkage()`](https://statmodels7.github.io/statmodels7/reference/plot_shrinkage.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
dd$y <- 1 + dd$x1 + 0.5 * dd$x2 + rnorm(100)
fit <- statmod(y ~ lasso(~ x1 + x2 + x3), distributions7::gaussian1_distrib(),
               dd)
ggplot2::autoplot(fit, term = "lasso(~x1 + x2 + x3)", type = plot_path())
```
