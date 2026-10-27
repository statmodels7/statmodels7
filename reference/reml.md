# Estimate the Hyperparameters by a Marginal Likelihood

Build the criterion object that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
uses to choose a model's smooth hyperparameters by a marginal
likelihood: the coefficients are integrated out by a Laplace
approximation at the penalized mode, and the remainder is maximized in
the hyperparameters. Pass the result as `outer_criterion`.

`reml()` integrates every coefficient that `marginal` does not name.
`ml()` integrates only the penalized directions and profiles the rest.
`reml()` is the default of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md).

## Usage

``` r
reml(information = observed_info(), marginal = NULL)

ml(information = observed_info(), marginal = NULL)
```

## Arguments

- information:

  The information that enters the determinant, an object that inherits
  from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html):

  - [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
    (the default): the curvature of the log-likelihood at the data,
    which makes the criterion the Laplace approximation;

  - [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html):
    the Fisher information, a function of the parameters alone. Its
    `method` selects how the expectation is obtained: `NULL` (the closed
    form of the family where it has one, and
    [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
    otherwise),
    [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
    [`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
    or
    [`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
    (see
    [distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html));

  - [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md):
    settled against the family by
    [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
    through
    [`outer_resolve()`](https://statmodels7.github.io/statmodels7/reference/outer_resolve.md),
    as the expected information for a Student t, a skew t or a Cauchy
    response, whose observed information is not positive definite at an
    outlier, and as the observed information otherwise.

  The expected information is used where the caller names it, as for a
  Student t whose degrees of freedom are modeled. The fitted model
  records the information used.

  On the expected information, the `method` object defines the matrix
  and every derivative of it that the exact outer gradient and Hessian
  read, so that the derivatives are those of the criterion that is
  evaluated. With the closed form of the family, the derivative is the
  closed form of the derivative (every shipped family; the Hessian also
  needs the second derivative of the information, which every shipped
  family supplies except
  [`distributions7::enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.html)).
  With
  [`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
  or
  [`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html),
  it is the derivative of the expectation, which moves the density as
  well as the integrand. With
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
  the matrix is the outer product of the scores at the data and the
  derivative is its exact derivative at fixed data. A family without a
  closed form reads
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  for the matrix and for the derivatives when `method` is `NULL`. The
  outer product depends on the response, so the criterion built on it
  has exact derivatives and is not the Laplace approximation of an
  expectation. With `method = NULL` on a family without a closed form,
  the matrix of the criterion is therefore an estimate of the expected
  information, noisier than a computed expectation and about 25 times
  cheaper: on a density-only family (n = 300, ten basis functions) a fit
  took 0.8 s against 25.3 s for
  `method = expect_integrate(identity = "score")`, with smoothing
  parameters 12.57 and 13.00, against 13.02 for the observed
  information. The computed expectation is selected with
  `method = expect_integrate(identity = "score")`.
  [`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
  draws a new sample at each call of the matrix and of each derivative,
  so that its derivative is that of the expectation and not of the
  sampled matrix, and the outer search is noisy at the level of the
  Monte Carlo error.
  [information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
  compares the three objects.

  For a model with a structural term
  ([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
  [`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html),
  a marginal break-point term) the criterion reads the observed
  information over the coefficients and the term's own parameters:
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  is settled as
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
  and
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  is an error.

- marginal:

  Which distribution parameters have their unpenalized coefficients
  estimated by maximizing this criterion, instead of being read at the
  joint mode (`ml()`) or integrated (`reml()`). `NULL`, the default,
  names every parameter except the position, which is the family's first
  parameter. `"none"` names no parameter, `"all"` names every one, and a
  character vector names the parameters listed. With `"all"`, `reml()`
  and `ml()` are the same criterion. For a multivariate family the
  default names no parameter, and naming one explicitly is an error.
  Where the model carries a penalty with a kink (a lasso, an elastic
  net, a SCAD or an MCP, in any equation), `NULL` names no parameter:
  the path that chooses the kinked penalty scores the model at the joint
  mode, and the fit returned is the model that it scored. A parameter
  named explicitly is still estimated on the criterion, at the value
  that the path chose. Where the family lacks the two derivatives that
  the criterion reads in some parameter (a Laplace response) and the
  model carries no smooth hyperparameter, the default criterion is not
  run at all, and a `reml()` or `ml()` passed by name is rejected.

## Value

An
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
object of kind `"reml"` or `"ml"`, whose `hessian` holds the object
given as `information` and whose `marginal` is as supplied; its path
settings are unused.

## The criterion

At the penalized mode \\\hat\beta(\theta)\\, \$\$\log L(\theta) =
\ell(\hat\beta) - \rho(\hat\beta;\theta) + \frac{q}{2}\log 2\pi -
\frac12\log\|A'(H+S)A\|,\$\$ with \\H\\ the information of the
log-likelihood, \\S\\ the penalty's second derivative in the
coefficients, and \\A\\ an orthonormal basis of the subspace integrated
over. The term \\\rho\\ needs no correction: a penalties7 penalty keeps
its normalizing constant, so it is exactly minus a log prior density,
and for a quadratic penalty that constant carries the
\\-\frac{r}{2}\log\lambda\\ and the log pseudo-determinant that a
marginal criterion needs. Written out, the expression reproduces Wood's
(2011) REML criterion term for term.

## The effective degrees of freedom of the fit

The effective degrees of freedom stored in the fitted model (`fit@edf`),
the share of each term in \\\mathrm{tr}\[(H+S)^{-1}H\]\\, are read on
the `information` of the criterion:
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
by default, or the
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
object with its `method`. The count is not part of the criterion here,
but the standard errors of the same fit read the same default matrix, so
the count and the standard errors of a summary agree.
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
recomputes the count when its `vcov_type` names the other matrix. For a
positive semidefinite \\H\\ the total lies in \\\[0, p\]\\, with \\p\\
the number of coefficients; an indefinite observed information can give
a direction of negative curvature a negative share, and the summary then
adds a note.

## What each one integrates

Write \\\gamma\\ for the unpenalized coefficients of the parameters that
`marginal` names, and \\u\\ for every other coefficient. The criterion
is maximized in \\\gamma\\ together with the hyperparameters. At each
value of \\\gamma\\, the coefficients \\u\\ are set to the penalized
mode.

`reml()` integrates every coefficient in \\u\\. An unpenalized one is
integrated under a flat prior; the absence of a penalty amounts to a
flat prior. `ml()` integrates only the directions in the range space of
the penalty, so an unpenalized coefficient in \\u\\ is read at the mode.
An ordinary covariate is unpenalized. The linear component of a
Demmler-Reinsch smooth is unpenalized too, because the penalty of that
smooth does not act on it.

With `marginal = "all"`, \\\gamma\\ holds every unpenalized coefficient.
The two criteria then integrate the same directions, and they are the
same function.

Where the criterion estimates coefficients, a coefficient that the model
does not identify at the first fit is left out of the determinant for
the whole search and stays free in the inner fit. This corresponds to
the REML log-likelihood of
[`stats::lm()`](https://rdrr.io/r/stats/lm.html)
(`logLik(fit, REML = TRUE)`), which uses the rank of the design in place
of its number of columns.

## Where a dispersion and a shape are estimated

By default `marginal` names every parameter except the position, which
is the family's first parameter: the location, or the scale for a family
with no location, such as
[`distributions7::weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.html).
For a gaussian mixed model this is the REML of `lme` and `glmmTMB`. The
coefficients of the mean are integrated, and \\\sigma\\ is maximized on
the criterion. Without any penalty, `reml()` returns the \\\sigma\\ of
[`stats::lm()`](https://rdrr.io/r/stats/lm.html),
\\\sqrt{\mathrm{rss}/(n-p)}\\, and its criterion equals
`logLik(lm(...), REML = TRUE)`. `ml()` returns
\\\sqrt{\mathrm{rss}/n}\\. If `agec` is defined as `age - 11`, the model
`distance ~ agec + random(~ agec | Subject)` on
[`nlme::Orthodont`](https://rdrr.io/pkg/nlme/man/Orthodont.html) gives
the estimates and the criterion of `lme`.

`marginal = "none"` gives the convention of `gamlss` and of mgcv's
`gaulss`: every parameter of the distribution is read at the joint mode
of the penalized likelihood. That estimate ignores the degrees of
freedom of the coefficients estimated beside it, so a dispersion read
there is biased downwards. Where a dispersion carries a random effect of
its own, the bias at the joint mode is larger than the bias on the
criterion.

Where some coefficients of the location are integrated, the estimate
depends on the parametrization of the family. The dispersion of
[`distributions7::gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.html)
and the variance of
[`distributions7::gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.html),
\\\sigma^2 = \phi\mu^2\\, differ under `reml()`, because the second
contains the integrated mean. With `"all"` the criterion does not depend
on the parametrization.

Some configurations are not covered, and
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)
lists them: a structural term anywhere in the model, a block that is a
working linearization (a sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html))
anywhere in the model, a block that moves with its coefficients in the
parameter's own equation, and a quadratic penalty in that equation whose
null space is not spanned by coordinates. In each case a parameter that
the default names keeps the convention of `"none"`, and a parameter
named explicitly is rejected. Beside a penalty with a kink the default
names no parameter, and a parameter named explicitly is estimated on the
criterion, and the coordinates that the kinked penalty covers are
excluded. For a multivariate family the default names no parameter, and
an explicit `marginal` other than `"none"` is an error.

## Which hyperparameters

The criterion estimates the hyperparameters of the terms fitted in one
system, meaning the terms whose penalty is twice differentiable. A
lasso, a SCAD or an MCP has a kink, its coefficients are estimated by a
method of their own, and a Laplace approximation is undefined at a point
where the second derivative does not exist. Those hyperparameters are
not on the outer search.
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
chooses them by a path under `sparse_criterion`, which is
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) by
default and can also be
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md).

## The exact gradient

The criterion has an exact gradient on the observed and on the expected
information, for every penalty whose derivatives in the coefficients and
the hyperparameters are available, as those of
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html),
[`modelterms7::te()`](https://statmodels7.github.io/modelterms7/reference/te.html)
and any
[`penalties7::quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.html)
are; a model with a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term has none. The gradient, and the Hessian where it is exact, are
supplied to the search. The default optimizer is
[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
with a bounded step where the Hessian is exact and steers the search,
[`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
where only the gradient is exact,
[`optimizers7::brent()`](https://statmodels7.github.io/optimizers7/reference/brent.html)
for a single coordinate without a gradient, and
[`optimizers7::nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html)
otherwise. A simplex needs one vertex per dimension and a quasi-Newton
method does not, so an exact gradient saves criterion evaluations once
there are two or more hyperparameters.

## ML needs a null basis

`ml()` reads the null space of every penalty over which it integrates,
because the null space determines which directions are profiled.
[`penalties7::is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.html)
identifies a penalty with no null space, and
[`penalties7::penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.html)
returns the null space of a quadratic, structured, additive or separable
penalty. A penalty that provides neither is rejected by name, and the
criterion does not assume a subspace for it.

## References

Wood, S. N. (2011). Fast stable restricted maximum likelihood and
marginal likelihood estimation of semiparametric generalized linear
models. *Journal of the Royal Statistical Society, Series B*, 73(1),
3–36.

## See also

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
for the prediction-error criteria,
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
for the argument `information`,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for where these are passed,
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
for reading back what they chose.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(200, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)

fit <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
               outer_criterion = reml())

# The smoothing parameter was estimated, and hyper() reports it.
hyper(fit)
#>   parameter                         term   name estimate  held source   id
#> 1        mu s(x, bspline_smooth(k = 10)) lambda  3.21673 FALSE   reml <NA>

# ML profiles the unpenalized directions instead of integrating them, so
# it shrinks a little less. The gap is small here because only two of the
# ten coefficients are unpenalized; it widens with the fixed effects.
fml <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
               outer_criterion = ml())
c(reml = hyper(fit)$estimate, ml = hyper(fml)$estimate)
#>     reml       ml 
#> 3.216730 3.205381 
c(reml = sum(fit@edf$edf), ml = sum(fml@edf$edf))
#>     reml       ml 
#> 8.555773 8.567577 

# The marginal log-likelihood, which these criteria maximize, is only
# available where one of them ran.
logLik(fit, loglik_type = "marginal")
#> 'log Lik.' -53.66395 (df=4)
```
