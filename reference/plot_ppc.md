# Comparison of the Response with Simulated Responses

Selects the graph that compares the observed response with responses
simulated from the fitted distribution. Without `statistic`, a
continuous family gives the observed density with the densities of the
simulated responses, and a discrete family the observed counts of each
value with the mean and the 90 percent range of the simulated counts.
With a function `statistic`, the graph is the histogram of its value on
the simulated responses with the observed value as a line.

## Usage

``` r
plot_ppc(statistic = NULL, n_sim = 200, seed = NULL)
```

## Arguments

- statistic:

  `NULL`, or a function of a numeric vector that returns a number.

- n_sim:

  The number of simulated responses.

- seed:

  A single number that seeds the simulation, or `NULL`.

## Value

A `PlotPpc` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The responses are drawn by
[`simulate.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/simulate.StatmodFit.md)
at the fitted parameters, `n_sim` of them, so they carry the variation
of the response and not the uncertainty of the estimates. The densities
are computed with
[`stats::density()`](https://rdrr.io/r/stats/density.html) on a common
grid over the range of the observed response. The statistic is a
function of a numeric vector that returns one number, such as `max` or
`function(y) mean(y == 0)`. The graph is not defined for a multivariate
family.

## See also

[`simulate.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/simulate.StatmodFit.md),
[`plot_qq()`](https://statmodels7.github.io/statmodels7/reference/plot_qq.md),
[`diagnostics()`](https://statmodels7.github.io/statmodels7/reference/diagnostics.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(120))
dd$y <- rpois(120, exp(1 + dd$x))
fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
ggplot2::autoplot(fit, plot_type = plot_ppc(n_sim = 50, seed = 1))

ggplot2::autoplot(fit, plot_type = plot_ppc(statistic = max, n_sim = 50,
                                       seed = 1))
```
