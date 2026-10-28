# Diagnostic Graphs of a Fitted Model

Draws a set of
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md)
graphs of a fit and returns them as a named list of ggplot2 objects, or
composed in one figure. The default set is the observed density of the
response with the densities of simulated responses
([`plot_ppc()`](https://statmodels7.github.io/statmodels7/reference/plot_ppc.md)),
the worm plot
([`plot_worm()`](https://statmodels7.github.io/statmodels7/reference/plot_worm.md)),
the calibration of the predicted quantiles
([`plot_calibration()`](https://statmodels7.github.io/statmodels7/reference/plot_calibration.md))
and the quantile residuals against the predicted median
([`plot_residuals()`](https://statmodels7.github.io/statmodels7/reference/plot_residuals.md)
with `against = "median"`). None of the four reads a particular
parameter, so the set applies to every family, including those without a
location parameter.

## Usage

``` r
diagnostics(fit, plot_types = NULL, data = NULL, compose = FALSE)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- plot_types:

  A list of objects inheriting from
  [diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md),
  or `NULL` for the default set.

- data:

  The data frame of the fit, or `NULL`; it is required by the graphs
  that read a covariate.

- compose:

  Whether to compose the graphs in one figure with patchwork.

## Value

A named list of ggplot2 objects, or a patchwork object when
`compose = TRUE`.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- 1 + dd$x + rnorm(150, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
d <- diagnostics(fit)
names(d)
#> [1] "ppc"         "worm"        "calibration" "residuals"  
diagnostics(fit, plot_types = list(plot_qq(), plot_pit()), compose = TRUE)
```
