# Unconditional Variance of a Fit

Selects the variance \\V' = V_b + J V\_\theta J'\\, which adds to the
Bayesian variance \\V_b\\ the variance that the estimation of the
hyperparameters (smoothing parameters, variances of random effects)
contributes. The Bayesian and the frequentist matrices read the
hyperparameters at the value where the outer search stopped, as though
it were known.

## Usage

``` r
unconditional_vcov(information = NULL)
```

## Arguments

- information:

  `NULL` (the default) or the matrix \\H\\, an object that inherits from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html):

  - `NULL`: the matrix is chosen for each fit by the family of the fit,
    as
    [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
    where the log-likelihood is twice differentiable in every parameter
    of the family (every
    [`distributions7::params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.html)
    is at least 2), and as
    [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
    with `method = NULL`, for the whole matrix, otherwise. The families
    with a kink in a parameter, such as the Laplace, have a singular
    observed information in that parameter and take the expected one;

  - [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html):
    the negative Hessian of the log-likelihood at the data;

  - [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html):
    the Fisher information, whose `method` is `NULL` (the closed form of
    the family where it has one,
    [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
    otherwise) or one of
    [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
    [`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
    and
    [`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
    (see
    [distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html)).

  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  is not accepted. For a model with a structural term, see
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).
  The choice does not depend on how the fit was computed
  (`inner_information` of
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)).
  [`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
  states the matrix when it is not the observed one, and a printed
  variance object states the information that it names.
  [information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
  compares the objects.

## Value

An `UnconditionalVcov` object, inheriting from
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
with the property `name` equal to `"unconditional"` and the property
`information`.

## Details

The correction is positive semi-definite and matters where a smoothing
parameter is poorly determined; it affects the penalized coordinates
much more than the unpenalized ones. It requires the curvature of the
outer criterion, so it costs more than the other two conventions. Where
no hyperparameter was estimated by a differentiable criterion the
correction is zero and the matrix equals the Bayesian one, and where the
curvature cannot be read the conditional matrix is returned with a
warning. It is the matrix that mgcv returns for `unconditional = TRUE`.
See
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the three conventions side by side and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the matrices in full.

## See also

[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
and
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
for the siblings,
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the comparison, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md)
for the functions that take it.

## Examples

``` r
set.seed(1)
ds <- data.frame(x = runif(200))
ds$y <- sin(2 * pi * ds$x) + rnorm(200, sd = 0.3)
fs <- statmod(y ~ s(x, bspline_smooth(k = 10)),
              distributions7::gaussian1_distrib(), ds)
unconditional_vcov()
#> <vcov> unconditional
sqrt(diag(vcov(fs, vcov_type = unconditional_vcov())))[["mu:s(x).z2"]]
#> [1] 0.06134753
unconditional_vcov(information = distributions7::expected_info())
#> <vcov> unconditional (expected information)
```
