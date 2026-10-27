# Case Resampling

Selects the case bootstrap: every replicate draws the observations with
replacement and refits the model on them.

## Usage

``` r
case_resample()
```

## Value

A `CaseResample` object, inheriting from
[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md).

## Details

The draw is a vector of counts, how many times each observation is
drawn, which multiplies the weights of the fit. The design is not
rebuilt, so the knots of a basis stay those of the fit. A replicate in
which a coefficient is not identified on the observations of positive
weight (a level of a factor that the draw left out) fails and is
replaced.

The case bootstrap does not depend on the variance specification of the
model, and its variance converges to the sandwich variance. It breaks
the groups of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term, and a fit with one is rejected:
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
resamples the groups. A fit with a structural term or a block that moves
with its coefficients is rejected as well, and a
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
criterion is rejected with `refit = "full"`.

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
fit <- add_bootstrap(fit, n = 30, resample = case_resample())
sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov(resample = case_resample()))))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>         0.1217035         0.1698995         0.1178554 
```
