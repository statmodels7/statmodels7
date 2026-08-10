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
  outer_method = NULL,
  outer_optimizer = NULL,
  hyper = NULL,
  start = NULL,
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

- outer_method:

  How the hyperparameters are estimated:
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  or `NULL` to hold them.

- outer_optimizer:

  The optimizer that searches over them, or `NULL` to let the
  availability of the exact gradient decide.

- hyper:

  Optional hyperparameters, a named list of named lists as
  `list(mu = list(lasso = c(lambda = 5)))`. They are held at these
  values.

- start:

  Optional starting coefficients, a named list.

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

**The budget and the stopping rule belong to the method.** There is no
`maxit` and no `tol` here: they are set on `inner_method`, which is
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)`(maxit =, tol =)`
or an optimizer with its own `maxit` and `criterion`, and the
alternation reads them from there (see
[`method_budget`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)).
Carrying a second copy would let a caller set both and be obeyed by
neither.

**The hyperparameters.** With `outer_method = NULL`, the default, each
one sits where `hyper` put it, or at the probe value of its bounds
otherwise – a placeholder rather than a choice, and it matters, since a
lasso at \\\lambda = 1\\ against an unaveraged log-likelihood of a few
hundred observations selects nothing at all. With
`outer_method = `[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
they are estimated by a marginal criterion, `outer_optimizer` searching
over them and the coefficients being refitted at each. Only a twice
differentiable penalty takes part: a lasso, a SCAD or an MCP keeps the
value it was given.

**Verbosity** has three levels, naming the loops rather than counting
them: `1` the outer search and the alternation, `2` the inner method's
own iterations, `3` the optimizers' traces as well. A named form is
accepted too, as `verbose = c(outer = TRUE, blocks = FALSE)`, since
watching the hyperparameters move while silencing a chatty inner
optimizer is the common case.

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
#> fitted in 32 ms, converged

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
#> fitted in 38 ms, converged
```
