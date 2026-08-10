# Fit a Model

Reads one formula carrying every parameter of a distribution, assembles
the terms it names into a penalized likelihood, and fits it.

## Usage

``` r
statmod(
  formula,
  distrib,
  data,
  weights = NULL,
  offsets = NULL,
  inner_method = iwls(),
  hyper = NULL,
  start = NULL,
  maxit = 50L,
  tol = 1e-08,
  verbose = 0
)
```

## Arguments

- formula:

  The model formula.

- distrib:

  A distributions7 distribution object.

- data:

  A data frame.

- weights:

  Optional prior weights, taken as given and not normalized.

- offsets:

  Optional named list of offsets, one per parameter.

- inner_method:

  How the smooth block is fitted:
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- hyper:

  Optional hyperparameters, a named list of named lists as
  `list(mu = list(lasso = c(lambda = 5)))`. They are held at these
  values.

- start:

  Optional starting coefficients, a named list.

- maxit:

  The alternation's iteration budget.

- tol:

  The alternation's tolerance, on the relative change in the objective.

- verbose:

  A level from 0 to 3, or a named logical vector.

## Value

An object of class
[`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Details

**The formula.** The equations of the distribution's parameters are
separated by `|`, the first carrying the response:

        y ~ x1 + ridge(R) + lasso(L)  |  sigma ~ z  |  nu ~ 1

A parameter with no equation gets an intercept. See
[`statmod_equations`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md),
whose recovery is not the obvious one.

**The fitting scheme.** The terms split in two by a property each one
already reports. Every term whose penalty is twice differentiable in its
coefficients – an unpenalized block, a ridge, a spline, a random effect
– is estimated in ONE system by `inner_method`, because their joint
curvature exists and using it is what makes a fit converge in a handful
of iterations. A term whose penalty has a kink – lasso, scad, mcp – is
estimated by a method of its own with everything else held fixed. The
fit alternates between the two until the objective and every block stop
moving.

**The objective is unaveraged**: minus the weighted log-likelihood plus
the penalties at full size, since a penalty is a negative log-prior and
a posterior adds the two at full size. What is scaled instead is the
stopping rule, so that a threshold means the same thing at \\n = 10\\
and at \\n = 10^7\\.

**The hyperparameters are held fixed.** Estimating a smoothing parameter
by an outer criterion is not written yet, so each one sits at the probe
value of its bounds unless `hyper` says otherwise. That value is a
placeholder rather than a choice, and it matters: a lasso at \\\lambda =
1\\ against an unaveraged log-likelihood of a few hundred observations
selects nothing at all.

**Verbosity** has three levels, naming the loops rather than counting
them: `1` the alternation, `2` the inner method's own iterations, `3`
the optimizers' traces as well. A named form is accepted too, as
`verbose = c(blocks = TRUE, inner = TRUE, optimizer = FALSE)`, since
watching the alternation while silencing a chatty inner optimizer is the
common case.

## See also

[`statmod_spec`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`loglik`](https://statmodels7.github.io/statmodels7/reference/loglik.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
fit
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ x, distrib = distributions7::gaussian1_distrib(), 
#>             data = dd)
#> 
#> Distribution: gaussian1
#> Observations: 60
#> 
#>   mu         ~ x
#>                linpar           2 coef
#>   sigma      ~ 1
#>                linpar           1 coef
#> 
#> log-likelihood -34.947195    objective 34.947195
#> fitted in 24 ms, converged

# every parameter can be modelled
statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ x | sigma ~ x, distrib = distributions7::gaussian1_distrib(), 
#>             data = dd)
#> 
#> Distribution: gaussian1
#> Observations: 60
#> 
#>   mu         ~ x
#>                linpar           2 coef
#>   sigma      ~ x
#>                linpar           2 coef
#> 
#> log-likelihood -33.446455    objective 33.446455
#> fitted in 36 ms, converged
```
