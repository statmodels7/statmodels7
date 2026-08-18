# Estimate the Hyperparameters by a Marginal Likelihood

`reml()` integrates every coefficient out of the likelihood before
maximizing in the hyperparameters; `ml()` integrates only the penalized
directions and profiles the rest.

## Usage

``` r
reml(hessian = c("observed", "expected"))

ml(hessian = c("observed", "expected"))
```

## Arguments

- hessian:

  Which information enters the determinant: `"expected"` or
  `"observed"`.

## Value

An
[`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Details

**The criterion.** At the penalized mode \\\hat\beta(\theta)\\, \$\$\log
L(\theta) = \ell(\hat\beta) - \rho(\hat\beta;\theta) + \frac{q}{2}\log
2\pi - \frac12\log\|A'(H+S)A\|,\$\$ with \\H\\ the information of the
log-likelihood, \\S\\ the penalty's second derivative in the
coefficients, and \\A\\ an orthonormal basis of the subspace integrated
over. Nothing is added to \\\rho\\ to make this work: a penalties7
penalty keeps its normalizing constant, so it is exactly minus a log
prior density, and for a quadratic penalty that constant carries the
\\-\frac{r}{2}\log\lambda\\ and the log pseudo-determinant that a
marginal criterion needs. Written out, the expression reproduces Wood's
(2011) REML criterion term for term.

**What each one integrates.** `reml()` takes \\A = I\\: every
coefficient is integrated, the unpenalized ones under the flat prior
their absence of a penalty amounts to. `ml()` takes \\A\\ spanning the
range space of the penalty, so a coefficient that is unpenalized – an
ordinary covariate, or the linear component of a Demmler-Reinsch smooth,
which its penalty leaves alone – is profiled rather than integrated.
This is the same distinction as between REML and ML for a variance
component in a mixed model, and `reml()` is the default for the same
reason: profiling a fixed effect leaves the estimate of the variance
biased downwards.

**Which hyperparameters.** Those of the terms fitted in one system,
which is to say those whose penalty is twice differentiable. A lasso, a
SCAD or an MCP has a kink, its coefficients are estimated by a method of
their own, and a Laplace approximation at a point where the second
derivative does not exist would be arithmetic without a meaning; those
hyperparameters stay where `hyper` put them.

**The criterion has an exact gradient** where the information is the
observed one and every penalty under estimation has a Hessian linear in
its hyperparameters, which covers `s()`, `te()` and any
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.html).
It is then supplied to the search and
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
becomes the default optimizer; otherwise the search compares values.
Measured, in evaluations of the criterion (each a whole inner fit)
against
[`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html):
40 against 32 with one smoothing parameter, 40 against 135 with two, 41
against 269 with three, and 12 against 283 with three and a modelled
scale. It does not pay in one dimension and pays from two on, a simplex
needing a vertex per dimension and a quasi-Newton method not. See
[`statmod_marginal_grad`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md).

**ML needs a null basis** for every penalty that has one, since that is
what says which directions are profiled.
[`is_proper`](https://statmodels7.github.io/penalties7/reference/is_proper.html)
answers for a penalty with no null space at all, and
[`penalty_null_basis`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.html)
for the quadratic and structured branches. A penalty that has neither is
rejected by name rather than integrated over a subspace guessed at.

## References

Wood, S. N. (2011). Fast stable restricted maximum likelihood and
marginal likelihood estimation of semiparametric generalized linear
models. *Journal of the Royal Statistical Society, Series B*, 73(1),
3–36.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(200, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
        outer_criterion = reml())
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ s(x, k = 10), distrib = distributions7::gaussian1_distrib(), 
#>             data = dd, outer_criterion = reml())
#> 
#> Distribution: gaussian1
#> Observations: 200
#> 
#>   mu         ~ s(x, k = 10)
#>                linpar           1 coef
#>                s(x, k = 10)     9 coef, edf 6.59
#>   sigma      ~ 1
#>                linpar           1 coef
#> 
#> log-likelihood -35.407348    objective 40.890371
#> REML -55.808442 over 5 hyperparameter evaluation(s)
#> fitted in 229 ms, converged
```
