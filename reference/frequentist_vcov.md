# Frequentist Variance of a Fit

Selects the sandwich variance \\V_f = (H+S)^{-1} H (H+S)^{-1}\\, with
\\H\\ the information of the log-likelihood and \\S\\ the second
derivative of the penalty. It is the sampling variance of the penalized
estimator at a fixed penalty, conditional on the hyperparameters.

## Usage

``` r
frequentist_vcov(information = NULL, sandwich = NULL)
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
  is not accepted. A model with a structural term is inverted on the
  observed information whatever the object says. The choice does not
  depend on how the fit was computed (`inner_information` of
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)).
  [`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
  states the matrix when it is not the observed one, and a printed
  variance object states the information that it names.
  [information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
  compares the objects.

- sandwich:

  `NULL` (the default) or an object from
  [`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md)
  or
  [`hc1_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md).
  `NULL` keeps \\H\\ as the meat, the matrix above. Otherwise the meat
  is built from the scores of the observations and the matrix is the
  robust sandwich variance described at
  [`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md).

## Value

A `FrequentistVcov` object, inheriting from
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
with the property `name` equal to `"frequentist"` and the properties
`information` and `sandwich`.

## Details

The matrix is smaller than the Bayesian one wherever there is a penalty,
because it omits the smoothing bias, and an interval built from it
covers less. With no penalty it equals \\H^{-1}\\. See
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the three conventions side by side and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the matrices in full.

## See also

[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
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
frequentist_vcov()
#> <vcov> frequentist
vcov(fit, vcov_type = frequentist_vcov())
#>                   mu:(Intercept)          mu:x sigma:(Intercept)
#> mu:(Intercept)      7.691473e-03 -1.161719e-02      1.055974e-18
#> mu:x               -1.161719e-02  2.209835e-02     -2.501978e-18
#> sigma:(Intercept)   1.055974e-18 -2.501978e-18      6.410256e-03
summary(fit, vcov_type = frequentist_vcov())
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ x, distrib = distributions7::gaussian1_distrib(), 
#>             data = dd)
#> 
#> Distribution: gaussian1     Observations: 80
#> 
#> === mu   [identity link]
#> 
#> Parametric terms
#>                estimate     se     z       p  lower upper
#>   (Intercept)     1.043 0.0877  11.9 < 1e-16 0.8714 1.215
#>   x               2.008 0.1487 13.51 < 1e-16  1.716 2.299
#> 
#> === sigma   [log link]
#> 
#> Parametric terms
#>                estimate      se     z       p lower   upper
#>   (Intercept)    -1.033 0.08006 -12.9 < 1e-16 -1.19 -0.8759
#> 
#> 95% intervals, frequentist variance
#> conditional log-likelihood -29.891019    effective df 3.00
#> cAIC 65.782    cBIC 72.928
#> fitted in 103 ms
#> inner   max |grad|/se 8.8e-16   min eigen 0.11
#> outer   max |grad|/se 8.2e-15   min eigen 1
frequentist_vcov(information = distributions7::expected_info())
#> <vcov> frequentist (expected information)
frequentist_vcov(sandwich = hc1_sandwich())
#> <vcov> frequentist with HC1 sandwich
```
