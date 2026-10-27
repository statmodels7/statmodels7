# Variance Conventions of a Fit

Three conventions set the variance matrix of the coefficients of a fit:
the Bayesian one, the frequentist one and the unconditional one. Each is
an object built by
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
or
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md),
and the object is passed as the argument `vcov_type` of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md)
and
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).
The default is
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md).

`vcov_method` is the abstract parent of the three classes. It cannot be
instantiated, and serves to test membership with
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html).
It has the property `name` only. The three conventions of this topic
compute the variance from the curvature of the likelihood and carry one
more property, `information`, which selects the matrix \\H\\ (see
below). A convention that does not read the curvature, such as
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md),
which is the covariance of refitted, resampled data sets, does not carry
it.

## Usage

``` r
vcov_method(name = character(0))
```

## Arguments

- name:

  A single string, the short name of the convention (`"bayesian"`,
  `"frequentist"` or `"unconditional"`), which is also the stem of the
  constructor
  ([`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
  and so on). A result prints it.

## Value

Nothing: the class is abstract and its constructor signals an error.
Every object that inherits from it has the property `name`; the three
conventions of this topic have the property `information` as well.

## Details

Write \\H\\ for the information of the log-likelihood and \\S\\ for the
second derivative of the penalty, both at the estimates. The argument
`information` of
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
selects \\H\\: `NULL` (the default), which is resolved for each fit by
the family,
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
the negative Hessian of the log-likelihood at the data, or
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html),
the Fisher information, with a `method` that selects how the expectation
is obtained. The default is the observed information where the
log-likelihood is twice differentiable in every parameter of the family
(every
[`distributions7::params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.html)
is at least 2), and the expected information, for the whole matrix,
otherwise: a family such as the Laplace has a singular observed
information in its location, and the closed form of the expected
information gives finite standard errors.
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
is not accepted, since the variance is defined on one named matrix.
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
states the matrix whenever it is not the observed one, and a printed
variance object states the information that it names. The effective
degrees of freedom that
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
prints are read on the same matrix. The choice is independent of the
information with which the fit was computed (`inner_information` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md))
and of the one of the outer criterion (the `information` of
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)).
A model that carries a structural term is inverted on the observed
information whatever the object says. The three matrices are \$\$V_b =
(H + S)^{-1}, \qquad V_f = (H+S)^{-1} H (H+S)^{-1}, \qquad V' = V_b + J
V\_\theta J'.\$\$

- [`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
  returns \\V_b\\, the posterior variance under the prior whose negative
  logarithm is the penalty, conditional on the hyperparameters. It
  carries the smoothing bias as variance, so an interval around a
  penalized term built from it covers at about its nominal rate.

- [`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
  returns \\V_f\\, the sandwich: the sampling variance of the penalized
  estimator at a fixed penalty. It is smaller than \\V_b\\ and covers
  less.

- [`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
  returns \\V'\\. The first two matrices read the hyperparameters
  (smoothing parameters, variances of random effects, scales of a prior)
  at the value where the outer search stopped, as though that value were
  known. \\V'\\ adds the variance that their estimation contributes,
  with \\V\_\theta\\ the variance of the estimated hyperparameters on
  their own free scale and \\J = -(H +
  S)^{-1}\\\partial^2\rho/\partial\beta\\\partial\theta\\ the movement
  of the penalized mode under them (Wood, Pya and Safken, 2016). The
  correction is positive semi-definite, so \\V'\\ is never narrower than
  \\V_b\\. It requires the curvature of the outer criterion, which the
  other two do not read, and it is the most expensive of the three.

With no penalty \\S = 0\\ and \\V_b\\ and \\V_f\\ are both \\H^{-1}\\
for a coefficient that the outer criterion does not estimate. Where the
observed information is singular at the estimates, as it is for a
coefficient on which the second derivative of the log-likelihood is zero
(the location of a Laplace response), the rows of the coefficients that
it does not identify are missing. The default `NULL` reads the expected
information for such a family, which is finite there. Where no
hyperparameter was estimated by a differentiable criterion the
correction of the unconditional matrix is zero and it equals \\V_b\\.
The coefficients that
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimate on the criterion (the unpenalized coefficients of every
parameter but the first, by default) take the curvature of that
criterion under the Bayesian and the unconditional conventions, which
therefore differ from the frequentist one for them even without a
penalty. The complete description of the matrices, of the coordinates
that are held and of the warnings is at
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## References

Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.

## See also

[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
for the three conventions,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the matrices,
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the tests that read a variance.

## Examples

``` r
set.seed(1)
ds <- data.frame(x = runif(200))
ds$y <- sin(2 * pi * ds$x) + rnorm(200, sd = 0.3)
fs <- statmod(y ~ s(x, bspline_smooth(k = 10)),
              distributions7::gaussian1_distrib(), ds)

# the standard error of one penalized coordinate under each convention
types <- list(frequentist_vcov(), bayesian_vcov(), unconditional_vcov())
vapply(types, function(ty) {
  sqrt(diag(vcov(fs, vcov_type = ty)))[["mu:s(x).z2"]]
}, numeric(1))
#> [1] 0.06053712 0.06083113 0.06134753

# the matrix H is chosen inside the convention
sqrt(diag(vcov(fs, vcov_type = bayesian_vcov(
  information = distributions7::expected_info()))))[["mu:s(x).z2"]]
#> [1] 0.06081721

# the class itself is abstract
try(vcov_method(name = "mine"))
#> Error in S7::new_object(S7::S7_object(), name = name) : 
#>   Can't construct an object from abstract class <vcov_method>
```
