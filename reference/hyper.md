# The Hyperparameters of a Fitted Model

One row per hyperparameter of every penalty the model carries, on the
scale the penalty declares them or on the free scale its links define,
with what put each value there.

## Usage

``` r
hyper(fit, scale = c("parameter", "link"))
```

## Arguments

- fit:

  A fitted model.

- scale:

  Which scale the values are reported on.

## Value

A data frame with `parameter`, `term`, `name`, `estimate`, `held` and
`source`, or a frame of no rows where the model carries no penalty.

## Details

A hyperparameter is not a coefficient and is not in
[`coef`](https://rdrr.io/r/stats/coef.html): it governs the coefficients
under it rather than sitting beside them, and the two are estimated by
different routes and reported with different qualifications. This is
where they are read.

The `parameter` scale is the one the penalty is written on, which is
what a reader wants: a smoothing parameter is a positive number and a
gaussian prior's `sigma` is a scale. The `link` scale is the free one
the outer search runs on, through each hyperparameter's own link, and is
what a caller comparing two fits' searches wants. Where a hyperparameter
carries no link the two coincide.

`source` says what put the value there, which `held` alone cannot: a
hyperparameter the term fixed reads `"fixed"`, one a marginal criterion
maximized reads that criterion's name, and one chosen along a path over
its own values reads the criterion that scored the path. A value chosen
along a path is the argument of a minimum over a grid rather than the
root of a derivative, so no standard error follows from it; one a
marginal criterion reached carries one, and
[`summary`](https://rdrr.io/r/base/summary.html) reports it.

## See also

[`coef.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md),
[`summary.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`statmod_held`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(80, 0, 1))
d$y <- sin(3 * d$x) + rnorm(80, 0, 0.3)
fit <- statmod(y ~ s(x, k = 6), distributions7::gaussian1_distrib(), d)
hyper(fit)
#>   parameter        term   name estimate  held source
#> 1        mu s(x, k = 6) lambda 29.67929 FALSE   reml
hyper(fit, scale = "link")
#>   parameter        term   name estimate  held source
#> 1        mu s(x, k = 6) lambda  3.39045 FALSE   reml
```
