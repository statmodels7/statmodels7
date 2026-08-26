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
reml(hessian = c("observed", "expected"))

ml(hessian = c("observed", "expected"))
```

## Arguments

- hessian:

  Which information enters the determinant: `"observed"` (the default)
  or `"expected"`. Matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html). The observed
  one is what the exact gradient needs.

## Value

An
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
object of kind `"reml"` or `"ml"`, with `hessian` as supplied and the
path settings unused.

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

`reml()` takes \\A = I\\: every coefficient is integrated, the
unpenalized ones under the flat prior their absence of a penalty amounts
to.

`ml()` takes \\A\\ spanning the range space of the penalty, so an
unpenalized coefficient is profiled. An ordinary covariate is one; so is
the linear component of a Demmler-Reinsch smooth, which that smooth's
penalty leaves alone.

This is the distinction between REML and ML for a variance component in
a mixed model, and `reml()` is the default for the same reason:
profiling a fixed effect leaves the variance estimate biased downwards.

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

fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
               outer_criterion = reml())

# The smoothing parameter was estimated, and hyper() says by what.
hyper(fit)
#>   parameter         term   name estimate  held source
#> 1        mu s(x, k = 10) lambda 3.245098 FALSE   reml

# ML profiles the unpenalized directions instead of integrating them, so
# it shrinks a little less. The gap is small here because only two of the
# ten coefficients are unpenalized; it widens with the fixed effects.
fml <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
               outer_criterion = ml())
c(reml = hyper(fit)$estimate, ml = hyper(fml)$estimate)
#>     reml       ml 
#> 3.245098 3.223994 
c(reml = sum(fit@edf$edf), ml = sum(fml@edf$edf))
#>     reml       ml 
#> 8.587662 8.593215 

# The marginal log-likelihood is what these maximize, and it is only
# available where one of them ran.
logLik(fit, type = "marginal")
#> 'log Lik.' -55.80844 (df=3)
```
