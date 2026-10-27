# Bootstrap Intervals for the Coefficients of a Fit

Three intervals are read from the replicates that
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
stores in a fit: the percentile interval, the basic interval and the BCa
interval. Each is an object built by `percentile_interval()`,
`basic_interval()` or `bca_interval()`, and the object is passed as the
argument `interval` of
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

`bootstrap_interval` is the abstract parent of the three classes. It
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

BcaInterval(
  name = character(0),
  label = character(0),
  n = NULL,
  resample = NULL,
  refit = NULL
)

percentile_interval(n = NULL, resample = NULL, refit = NULL)

basic_interval(n = NULL, resample = NULL, refit = NULL)

bca_interval(n = NULL, resample = NULL, refit = NULL)
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

`percentile_interval()`, `basic_interval()` and `bca_interval()` return
an object that inherits from `bootstrap_interval`, with the properties
`name`, `label`, `n`, `resample` and `refit`. `bootstrap_interval`
returns nothing: its constructor signals an error.

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

The BCa interval (bias-corrected and accelerated) takes the quantiles of
the replicates at adjusted levels, \$\$\Phi\\\left(z_0 + \frac{z_0 +
z\_{\alpha'}}{1 - a\\(z_0 + z\_{\alpha'})} \right), \quad \alpha' \in
\\\alpha/2,\\ 1 - \alpha/2\\,\$\$ where \\\Phi\\ is the standard normal
distribution function and \\z\_{\alpha'}\\ its quantile. The bias
correction \\z_0 = \Phi^{-1}(\hat p)\\ uses the proportion \\\hat p\\ of
replicates below the estimate, with the replicates equal to it counted
by halves. The acceleration is \$\$a = \frac{\sum_i
(\bar\theta\_{(\cdot)} - \hat\theta\_{(-i)})^3} {6\\\[\sum_i
(\bar\theta\_{(\cdot)} - \hat\theta\_{(-i)})^2\]^{3/2}},\$\$ with
\\\hat\theta\_{(-i)}\\ the estimate of the fit refitted without
observation \\i\\ (the weight of the observation set to zero) and
\\\bar\theta\_{(\cdot)}\\ the mean of these estimates (Efron 1987). With
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
the unit deleted is a cluster. The refits are made once and stored in
the entry, one per observation or cluster, with the settings of the
entry (`refit`), so that the first call of
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
costs as many refits as there are units, and later calls read them. The
interval is available for
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md),
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
and
[`bayesian_resample()`](https://statmodels7.github.io/statmodels7/reference/bayesian_resample.md).
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md)
is rejected, because the acceleration needs the influence of the
observations of the data, and so is
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
beside a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term. A coordinate whose estimate lies outside the range of its
replicates, or whose adjusted level is undefined, has `NA` limits.

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

Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
American Statistical Association*, 82(397), 171–185.

## See also

[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + 2 * dd$x + rnorm(40, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
fit <- add_bootstrap(fit, n = 60, resample = case_resample())
confint(fit, interval = percentile_interval(n = 60))
#> Warning: The 60 replicates leave 1.5 in each tail of a 95 per cent interval.
#>   The quantiles of a tail are unstable below 10 replicates: use n of at least 400.
#>                   parameter   term coefficient  estimate         se      lower
#> mu:(Intercept)           mu linpar (Intercept)  1.120104 0.09288556  0.9708948
#> mu:x                     mu linpar           x  1.818008 0.15386898  1.5431546
#> sigma:(Intercept)     sigma linpar (Intercept) -1.097264 0.12477276 -1.4334584
#>                        upper
#> mu:(Intercept)     1.3282745
#> mu:x               2.0919233
#> sigma:(Intercept) -0.9032555
confint(fit, interval = basic_interval(n = 60))
#> Warning: The 60 replicates leave 1.5 in each tail of a 95 per cent interval.
#>   The quantiles of a tail are unstable below 10 replicates: use n of at least 400.
#>                   parameter   term coefficient  estimate         se     lower
#> mu:(Intercept)           mu linpar (Intercept)  1.120104 0.09288556  0.911934
#> mu:x                     mu linpar           x  1.818008 0.15386898  1.544093
#> sigma:(Intercept)     sigma linpar (Intercept) -1.097264 0.12477276 -1.291272
#>                        upper
#> mu:(Intercept)     1.2693137
#> mu:x               2.0928620
#> sigma:(Intercept) -0.7610695
confint(fit, interval = bca_interval(n = 60))
#> Warning: The 60 replicates leave 1.5 in each tail of a 95 per cent interval.
#>   The quantiles of a tail are unstable below 10 replicates: use n of at least 400.
#> bootstrap: computed 40 refits without one observation (case, refit coefficients).
#> Warning: The BCa levels of some coordinate fall in the first or last replicate, so an extreme order statistic is a limit.
#>   Use more replicates.
#>                   parameter   term coefficient  estimate         se      lower
#> mu:(Intercept)           mu linpar (Intercept)  1.120104 0.09288556  0.9807196
#> mu:x                     mu linpar           x  1.818008 0.15386898  1.5329340
#> sigma:(Intercept)     sigma linpar (Intercept) -1.097264 0.12477276 -1.2517273
#>                        upper
#> mu:(Intercept)     1.3363865
#> mu:x               2.0761335
#> sigma:(Intercept) -0.8839039
```
