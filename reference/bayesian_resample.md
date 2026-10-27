# Bayesian Resampling

Selects the Bayesian bootstrap: every replicate gives each observation a
random positive weight and refits the model on the weighted data.

## Usage

``` r
bayesian_resample()
```

## Value

A `BayesianResample` object, inheriting from
[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md).

## Details

The weights of a replicate are \\n\\ times a draw from the Dirichlet
distribution with all parameters equal to one, which is a vector of
independent standard exponential variables divided by their mean. They
sum to \\n\\ and multiply the weights of the fit. The design is not
rebuilt, so the knots of a basis stay those of the fit. An observation
with prior weight zero keeps weight zero and is not drawn.

Every observation of a replicate has a positive weight, so no replicate
loses a level of a factor, and a rare level does not make a coefficient
unidentified or a logistic regression separable, as the case scheme does
when it leaves all the observations of the level out. The variance of
the replicates is that of the case scheme times a factor close to
\\n/(n + 1)\\, and it converges to the sandwich variance in the same
way. The scheme breaks the groups of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term, and a fit with one is rejected, like the case scheme. A fit with a
structural term is rejected as well, and a
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
criterion is rejected with `refit = "full"`.

## References

Rubin, D. B. (1981). The Bayesian bootstrap. *The Annals of Statistics*,
9(1), 130–134.

Newton, M. A. and Raftery, A. E. (1994). Approximate Bayesian inference
with the weighted likelihood bootstrap. *Journal of the Royal
Statistical Society: Series B*, 56(1), 3–48.

## See also

[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
for the comparison of the schemes,
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
fit <- add_bootstrap(fit, n = 30, resample = bayesian_resample())
sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov(
  resample = bayesian_resample()))))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.08195906        0.14916642        0.08869258 
```
