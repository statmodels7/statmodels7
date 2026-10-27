# Bootstrap Variance of a Fit

Selects the variance that is the covariance of the coefficients of
refitted, resampled data sets. The matrix is not read from the curvature
at the mode: it is the empirical covariance of the replicates that
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
stores in the fit.

## Usage

``` r
bootstrap_vcov(n = NULL, resample = NULL, refit = NULL)
```

## Arguments

- n:

  A single whole number, the total number of replicates, or `NULL`.

- resample:

  A
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  object, or `NULL`.

- refit:

  `"coefficients"`, `"full"` or `NULL`.

## Value

A `BootstrapVcov` object, inheriting from
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
with the property `name` equal to `"bootstrap"` and the properties `n`,
`resample` and `refit`.

## Details

The arguments that are `NULL` are resolved against the replicates stored
in the fit. A `NULL` argument takes the value of the stored entry if
exactly one entry matches the arguments that were given, and the default
otherwise (`n = 200`,
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md),
`"coefficients"`). If several entries match, an error prints the table
of entries. A setting that matches no stored entry is computed at the
call and stored in the fit, with a message that says how many replicates
were computed. `n` is the total that is wanted: with an entry of \\m\\
replicates, \\n \> m\\ computes the \\n - m\\ missing ones, and \\n \le
m\\ uses the first \\n\\ without computing anything.
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
resolve the variance once, so the replicates are computed once.

The coordinates are those of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md):
the coefficients of every equation, with the aliased ones missing. A
coefficient under a kinked penalty
([`modelterms7::lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.html),
[`modelterms7::scad()`](https://statmodels7.github.io/modelterms7/reference/scad.html),
[`modelterms7::mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.html))
that the fit holds at zero has a variance here, which includes the mass
of the replicates at zero, where the other conventions report it as
missing. With `readable = TRUE` the readable quantities are computed for
each replicate and the covariance is taken afterwards.

The variance estimates the sampling variance of the penalized estimator.
It does not contain the smoothing bias that
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
carries as variance, so an interval around a penalized term built from
it covers less than the nominal level, as one built from
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
does. With `refit = "coefficients"` the hyperparameters are held at the
values of the fit, and the matrix is a relative of
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md).
With `refit = "full"` the hyperparameters are selected again on every
replicate, and the variability of the hyperparameters enters; it costs
about ten times more. The relative Monte Carlo error of a standard error
is about \\1/\sqrt{2n}\\. A replicate that fails (an error, a mode error
above
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
a coefficient that is not identified) is replaced; replacing conditions
the bootstrap on success, which understates the variability when
failures are frequent.

The convention has no `information` argument: the matrix is not built
from the information. The refits use the `inner_information` of the fit.

## See also

[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
which computes and stores the replicates,
[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md),
which reads them, and
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
bootstrap_vcov()
#> <vcov> bootstrap (n as stored, resample as stored, refit as stored)
fit <- add_bootstrap(fit, n = 40, resample = case_resample())
sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov())))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>         0.1263922         0.1736775         0.1165450 
sqrt(diag(vcov(fit)))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.09832324        0.17020018        0.09284767 
```
