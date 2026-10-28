# Bayesian Variance of a Fit

Selects the variance \\V_b = (H + S)^{-1}\\, the inverse of the
penalized information: the information \\H\\ of the log-likelihood plus
the second derivative \\S\\ of the penalty. It is the posterior variance
under the prior whose negative logarithm is the penalty, conditional on
the hyperparameters, and it is the default of every function that takes
a `vcov_type`.

## Usage

``` r
bayesian_vcov(information = NULL)
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

A `BayesianVcov` object, inheriting from
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
with the property `name` equal to `"bayesian"` and the property
`information`.

## Details

Choose it for an interval around a penalized term: it carries the
smoothing bias as variance, and such intervals cover at about their
nominal rate. It does not include the uncertainty of the
hyperparameters, for which
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
is the convention. See
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the three side by side and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the matrices in full.

## See also

[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
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
dd <- data.frame(x = runif(80))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
bayesian_vcov()
#> <vcov> bayesian
vcov(fit, vcov_type = bayesian_vcov())
#>                   mu:(Intercept)          mu:x sigma:(Intercept)
#> mu:(Intercept)      7.691473e-03 -1.161719e-02      1.055974e-18
#> mu:x               -1.161719e-02  2.209835e-02     -2.501978e-18
#> sigma:(Intercept)   1.055974e-18 -2.501978e-18      6.410256e-03
confint(fit, vcov_type = bayesian_vcov())
#>                   parameter   term coefficient  estimate         se      lower
#> mu:(Intercept)           mu linpar (Intercept)  1.043291 0.08770104  0.8714004
#> mu:x                     mu linpar           x  2.007856 0.14865515  1.7164973
#> sigma:(Intercept)     sigma linpar (Intercept) -1.032801 0.08006408 -1.1897235
#>                        upper
#> mu:(Intercept)     1.2151822
#> mu:x               2.2992148
#> sigma:(Intercept) -0.8758781

# the standard errors from the expected information instead
vcov(fit, vcov_type = bayesian_vcov(
  information = distributions7::expected_info()))
#>                   mu:(Intercept)        mu:x sigma:(Intercept)
#> mu:(Intercept)       0.007691473 -0.01161719       0.000000000
#> mu:x                -0.011617186  0.02209835       0.000000000
#> sigma:(Intercept)    0.000000000  0.00000000       0.006410256
```
