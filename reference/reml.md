# Estimate the Hyperparameters by a Marginal Likelihood

Build the object that tells
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
to choose a model's smooth hyperparameters by a marginal likelihood: the
coefficients are integrated out by a Laplace approximation at the
penalized mode, and what is left is maximized in the hyperparameters.
Pass the result as `outer_criterion`.

`reml()` integrates **every** coefficient out. `ml()` integrates only
the penalized directions and profiles the rest. `reml()` is
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
default.

## Usage

``` r
reml(hessian = c("observed", "auto", "expected"), marginal = NULL)

ml(hessian = c("observed", "auto", "expected"), marginal = NULL)
```

## Arguments

- hessian:

  Which information enters the determinant: `"observed"` (the default),
  the curvature of the log-likelihood at the data, which makes the
  criterion the Laplace approximation, `"expected"`, the Fisher
  information, a function of the parameters alone, or `"auto"`, which
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  settles against the family through
  [`outer_resolve()`](https://statmodels7.github.io/statmodels7/reference/outer_resolve.md):
  the expected information for a Student t, a skew t or a Cauchy
  response, whose observed information is not positive definite at an
  outlier, and the observed information otherwise. The observed
  information is the default because, once a degrees-of-freedom start
  run to its limit is replaced (statmodels7 0.193.0), it reaches the
  same point as the expected one on every case measured and 2 to 20
  times faster; the case that needs the expected information
  ([`MASS::GAGurine`](https://rdrr.io/pkg/MASS/man/GAGurine.html) with
  `nu ~ Age`) is one where the caller names it (Giovanni, 2026-10-07).
  The fitted model records the information used. Matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html). Both carry an
  exact outer gradient and Hessian wherever the family writes its
  expected information out, which every shipped family does since
  distributions7 0.65.0. A family that does not would read the outer
  product of its scores at the data under the default
  `iwls(approx = "opg")`, which is not an expectation, and
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  rejects `"expected"` there: see
  [`assert_criterion_information()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_information.md).

- marginal:

  Which distribution parameters have their unpenalized coefficients
  estimated by maximizing this criterion, instead of being read at the
  joint mode (`ml()`) or integrated (`reml()`). `NULL`, the default,
  names every parameter except the position, which is the family's first
  parameter. `"none"` names no parameter, `"all"` names every one, and a
  character vector names the parameters listed. With `"all"`, `reml()`
  and `ml()` are the same criterion. Where the model carries a penalty
  with a kink (a lasso, an elastic net, a SCAD or an MCP, in any
  equation), `NULL` names no parameter: the path that chooses the kinked
  penalty scores the model at the joint mode, and the fit returned is
  the model it scored. A parameter named explicitly is still estimated
  on the criterion, at the value the path chose. Where the family lacks
  the two derivatives the criterion reads in some parameter (a Laplace
  response) and the model carries no smooth hyperparameter, the default
  criterion is not run at all, and a `reml()` or `ml()` passed by name
  is rejected.

## Value

An
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
object of kind `"reml"` or `"ml"`, with `hessian` and `marginal` as
supplied and the path settings unused.

## The criterion

At the penalized mode \\\hat\beta(\theta)\\, \$\$\log L(\theta) =
\ell(\hat\beta) - \rho(\hat\beta;\theta) + \frac{q}{2}\log 2\pi -
\frac12\log\|A'(H+S)A\|,\$\$ with \\H\\ the information of the
log-likelihood, \\S\\ the penalty's second derivative in the
coefficients, and \\A\\ an orthonormal basis of the subspace integrated
over. Nothing is added to \\\rho\\ to make this work: a penalties7
penalty keeps its normalizing constant, so it is exactly minus a log
prior density, and for a quadratic penalty that constant carries the
\\-\frac{r}{2}\log\lambda\\ and the log pseudo-determinant that a
marginal criterion needs. Written out, the expression reproduces Wood's
(2011) REML criterion term for term.

## What each one integrates

Write \\\gamma\\ for the unpenalized coefficients of the parameters that
`marginal` names, and \\u\\ for every other coefficient. The criterion
is maximized in \\\gamma\\ together with the hyperparameters. At each
value of \\\gamma\\, the coefficients \\u\\ are set to the penalized
mode.

`reml()` integrates every coefficient in \\u\\. An unpenalized one is
integrated under a flat prior, which is what the absence of a penalty
amounts to. `ml()` integrates only the directions in the range space of
the penalty, so an unpenalized coefficient in \\u\\ is read at the mode.
An ordinary covariate is unpenalized. The linear component of a
Demmler-Reinsch smooth is unpenalized too, because the penalty of that
smooth does not act on it.

With `marginal = "all"`, \\\gamma\\ holds every unpenalized coefficient.
The two criteria then integrate the same directions, and they are the
same function: measured on six families with a random intercept, the two
maxima and the two sets of coefficients agree exactly.

A coefficient that the model does not identify at the first fit is left
out of the determinant for the whole search, and it stays free in the
inner fit. This is how [`stats::lm()`](https://rdrr.io/r/stats/lm.html)
computes its REML criterion, with the rank of the design instead of its
number of columns.

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
\\\sqrt{\mathrm{rss}/n}\\. On
[`nlme::Orthodont`](https://rdrr.io/pkg/nlme/man/Orthodont.html),
`distance ~ agec + random(~ agec | Subject)` gives \\\sigma\\ =
1.310039, standard deviations 2.134332 and 0.226429, and a criterion of
-221.31834, which are the values of `lme` to six digits.

`marginal = "none"` gives the convention of `gamlss` and of mgcv's
`gaulss`: every parameter of the distribution is read at the joint mode
of the penalized likelihood. That estimate ignores the degrees of
freedom of the coefficients estimated beside it, so a dispersion read
there is biased downwards. Measured over 20 replicates of 40 groups of
8, a gamma dispersion with a random effect of its own is biased by -8.0
per cent at the joint mode and by -2.6 per cent on the criterion, and a
Weibull shape with a random intercept by +8.7 and +1.1 per cent.

Where some coefficients of the location are integrated, the estimate
depends on the parametrization of the family. The dispersion of
[`distributions7::gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.html)
and the variance of
[`distributions7::gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.html),
\\\sigma^2 = \phi\mu^2\\, differ by about one per cent under `reml()`,
because the second contains the integrated mean. With `"all"` the
criterion does not depend on the parametrization.

Five configurations are not covered yet, and
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)
lists them: a structural term, a block that moves with its coefficients,
a penalty with a kink, a block that is a working linearization (a sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)),
and a penalty whose null space is not spanned by coordinates. In each
case a parameter that the default names keeps the convention of
`"none"`, and a parameter named explicitly is refused.

## Which hyperparameters

Those of the terms fitted in one system, meaning those whose penalty is
twice differentiable. A lasso, a SCAD or an MCP has a kink, its
coefficients are estimated by a method of their own, and a Laplace
approximation at a point where the second derivative does not exist
would be arithmetic with no meaning. Those hyperparameters stay where
their term left them, and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) is
what
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
puts to them instead.

## The exact gradient

The criterion has one where the information is the observed one and
every penalty under estimation has a Hessian linear in its
hyperparameters, which covers
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html),
[`modelterms7::te()`](https://statmodels7.github.io/modelterms7/reference/te.html)
and any
[`penalties7::quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.html).
It is then supplied to the search and
[`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
becomes the default optimizer; otherwise the search compares values.

Measured in evaluations of the criterion, each a whole inner fit,
against
[`optimizers7::nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html):
40 against 32 with one smoothing parameter, 40 against 135 with two, 41
against 269 with three, and 12 against 283 with three and a modeled
scale. It does not pay in one dimension and pays from two on, a simplex
needing a vertex per dimension and a quasi-Newton method not.

## ML needs a null basis

For every penalty that has one, since that is what says which directions
are profiled.
[`penalties7::is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.html)
answers for a penalty with no null space at all, and
[`penalties7::penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.html)
for the quadratic and structured branches. A penalty offering neither is
rejected by name, instead of being integrated over a subspace guessed
at.

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

# The smoothing parameter was estimated, and hyper() says by what.
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
#> 8.562273 8.574067 

# The marginal log-likelihood is what these maximize, and it is only
# available where one of them ran.
logLik(fit, type = "marginal")
#> 'log Lik.' -53.66395 (df=3)
```
