# Bootstrap Intervals for the Coefficients of a Fit

Two intervals are read from the replicates that
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
stores in a fit: the percentile interval and the basic interval. Each is
an object built by `percentile_interval()` or `basic_interval()`, and
the object is passed as the argument `interval` of
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

`bootstrap_interval` is the abstract parent of the two classes. It
cannot be instantiated, and serves to test membership with
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html).

## Usage

``` r
bootstrap_interval(
  name = character(0),
  label = character(0),
  n = NULL,
  resample = NULL,
  refit = NULL
)

PercentileInterval(
  name = character(0),
  label = character(0),
  n = NULL,
  resample = NULL,
  refit = NULL
)

BasicInterval(
  name = character(0),
  label = character(0),
  n = NULL,
  resample = NULL,
  refit = NULL
)

percentile_interval(n = NULL, resample = NULL, refit = NULL)

basic_interval(n = NULL, resample = NULL, refit = NULL)
```

## Arguments

- name:

  A single string, the short name of the interval.

- label:

  A single string, the name under which an object prints.

- n:

  A single whole number, the total number of replicates, or `NULL`.

- resample:

  A
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  object, or `NULL`.

- refit:

  `"coefficients"`, `"full"` or `NULL`.

## Value

`percentile_interval()` and `basic_interval()` return an object that
inherits from `bootstrap_interval`, with the properties `name`, `label`,
`n`, `resample` and `refit`. `bootstrap_interval` returns nothing: its
constructor signals an error.

## Details

Write \\\hat\theta\\ for an estimate, \\\theta^\*\_{(\alpha)}\\ for the
\\\alpha\\ quantile of its bootstrap replicates and \\\alpha = 1 -\\
`level`. The percentile interval is \$\$\[\theta^\*\_{(\alpha/2)},\\
\theta^\*\_{(1 - \alpha/2)}\],\$\$ and the basic interval is
\$\$\[2\hat\theta - \theta^\*\_{(1 - \alpha/2)},\\ 2\hat\theta -
\theta^\*\_{(\alpha/2)}\],\$\$ the reflection of the percentile interval
about the estimate (Davison and Hinkley 1997, chapter 5). The quantiles
are those of
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) with its
default type.

The quantity is the one that
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
reports. With `readable = TRUE` the readable quantities are computed on
every replicate and the quantiles are taken afterwards, so the
percentile interval of a transformed quantity is the transformed
percentile interval. A coordinate that is missing in some replicate (the
individual random effects under
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md))
or in the fit (an aliased coefficient) has `NA` limits.

The replicates are those of the stored entry that the arguments select,
in the way
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
selects it: a `NULL` argument takes the value of the single entry that
matches the others, and the defaults when none matches. An entry with
fewer than `n` replicates is extended, and one with more is read up to
its first `n`. The default when nothing is stored is `n = 1000`, because
a tail of a 95 per cent interval holds only one replicate in forty. A
warning is signaled when a tail holds fewer than ten replicates, that
is, when \\n \alpha / 2 \< 10\\.

## References

Davison, A. C. and Hinkley, D. V. (1997). *Bootstrap Methods and Their
Application*. Cambridge University Press.

## See also

[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
fit <- add_bootstrap(fit, n = 100, resample = case_resample())
confint(fit, interval = percentile_interval(n = 100))
#> Warning: The 100 replicates leave 2.5 in each tail of a 95 per cent interval.
#>   The quantiles of a tail are unstable below 10 replicates: use n of at least 400.
#>                   parameter   term coefficient  estimate        se      lower
#> mu:(Intercept)           mu linpar (Intercept)  1.071069 0.1084330  0.8746006
#> mu:x                     mu linpar           x  1.956143 0.1567222  1.6758008
#> sigma:(Intercept)     sigma linpar (Intercept) -1.042678 0.1093931 -1.3106504
#>                        upper
#> mu:(Intercept)     1.2766207
#> mu:x               2.2323064
#> sigma:(Intercept) -0.8935903
confint(fit, interval = basic_interval(n = 100))
#> Warning: The 100 replicates leave 2.5 in each tail of a 95 per cent interval.
#>   The quantiles of a tail are unstable below 10 replicates: use n of at least 400.
#>                   parameter   term coefficient  estimate        se     lower
#> mu:(Intercept)           mu linpar (Intercept)  1.071069 0.1084330  0.865518
#> mu:x                     mu linpar           x  1.956143 0.1567222  1.679979
#> sigma:(Intercept)     sigma linpar (Intercept) -1.042678 0.1093931 -1.191766
#>                        upper
#> mu:(Intercept)     1.2675380
#> mu:x               2.2364845
#> sigma:(Intercept) -0.7747057
```
