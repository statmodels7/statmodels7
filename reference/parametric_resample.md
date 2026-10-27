# Parametric Resampling

Selects the parametric bootstrap: every replicate draws a new response
from the family at the fitted parameters and refits the model on it.

## Usage

``` r
parametric_resample()
```

## Value

A `ParametricResample` object, inheriting from
[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md).

## Details

The effects of every
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term are drawn afresh from the prior that the fit estimated, so the
groups of the data are new groups in every replicate. The replicate of a
random-effect coordinate is stored as \\\hat b^\* - b^\* + \hat b\\, the
estimate of the replicate minus the effect that was drawn, added to the
estimate of the fit. The covariance of these replicates is the variance
of the prediction error \\\mathrm{Var}(\hat b - b)\\, which is the
quantity that
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
reports for those coordinates. The covariance of \\\hat b^\*\\ alone
would be the variance of the shrunken estimate, which is smaller
wherever the shrinkage is strong.

The parametric bootstrap reproduces the variance that the model
attributes to the estimator. It does not detect a misspecified variance
function; the scheme of
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md)
does. A fit with a structural term
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html))
or a block that moves with its coefficients
([`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html))
is rejected.

## See also

[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
for the comparison of the three schemes,
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
and
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
for the functions that take it.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
fit <- add_bootstrap(fit, n = 30, resample = parametric_resample())
statmod_bootstrap(fit)
#>     resample        refit replicates replaced error mode_error unidentified
#> 1 parametric coefficients         30        0     0          0            0
```
