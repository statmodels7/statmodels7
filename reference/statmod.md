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
  inner_optimizer = iwls(),
  outer_criterion = reml(),
  sparse_criterion = bic(),
  outer_optimizer = NULL,
  start = NULL,
  linpar_control = linpar_options(),
  verbose = 0,
  ...
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

- inner_optimizer:

  How the smooth block is fitted:
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- outer_criterion:

  How the SMOOTH hyperparameters are estimated:
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  (the default),
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
  or `NULL` to hold them where they are.

- sparse_criterion:

  How the hyperparameter of a KINKED penalty – lasso, scad, mcp – is
  chosen:
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  (the default),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
  or `NULL` to hold it. A marginal criterion is rejected here, being
  read at a mode that sits on the kink.

- outer_optimizer:

  The optimizer that searches over them, or `NULL` to let the
  availability of the exact gradient decide.

- start:

  Where the fit begins: a named list of coefficients, a
  [`start_strategy`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
  such as
  [`start_search`](https://statmodels7.github.io/statmodels7/reference/start_search.md),
  or `NULL` for
  [`start_intercepts`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md).
  A strategy is asked once, before the alternation between the
  coefficients and the hyperparameters begins, which is why a global
  search belongs here and not in `inner_optimizer`: there it would rerun
  at every hyperparameter the outer search tried.

- linpar_control:

  How the unpenalized parametric block is built, as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns it: the storage and the contrasts for its factors. It reaches
  the IMPLICIT term, the one the bare covariates collapse into and which
  a caller never writes; a
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html)
  written out takes them directly. The argument and the function are
  named differently on purpose, as
  [`glm`](https://rdrr.io/r/stats/glm.html)'s `control` and
  [`glm.control`](https://rdrr.io/r/stats/glm.control.html) are.

- verbose:

  A level from 0 to 3, or a named logical vector.

- ...:

  Not used, and reported. `hyper` was removed from here: a
  hyperparameter is held in the term that carries the penalty, and a
  second place to say so would be read by nobody whenever the two
  disagreed.

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
– is estimated in ONE system by `inner_optimizer`, because their joint
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
`maxit` and no `tol` here: they are set on `inner_optimizer`, which is
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)`(maxit =, tol =)`
or an optimizer with its own `maxit` and `criterion`, and the
alternation reads them from there (see
[`method_budget`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)).
Carrying a second copy would let a caller set both and be obeyed by
neither.

**Every hyperparameter is ESTIMATED unless its own term holds it.**
Which ones are held is said by the TERM that carries the penalty –
`lasso(x, lambda = 3)`, `ridge(x, sigma = 0.5)`, `s(x, lambda = 2)`,
`enet(x, alpha = 0.5)` – and everything left `NULL`, which is the
default, is chosen from the data. The term is where the penalty is named
and so where that belongs; an argument here saying the same thing would
be read by nobody whenever the two disagreed.
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimates the smooth ones, with `outer_optimizer` searching over them
and the coefficients refitted at each.

A KINKED penalty is a different instrument and has its own argument.
`sparse_criterion`,
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) by
default, sweeps it along a PATH of its own values – from the kink that
empties the block down to `min_ratio` of it – because the penalized mode
is only piecewise smooth in that hyperparameter, turning a corner
whenever a coefficient joins the active set or leaves it, so a criterion
read there inherits the corners and a gradient search reads a slope
about to change. Where a model carries both kinds the path is outside
and the marginal criterion is estimated inside each of its points, so a
smoothing parameter can come from REML and a lasso's \\\lambda\\ from
BIC in the same fit.

The top of that path is DATA-DEPENDENT and depends on the rest of the
model: it is the kink that empties the block, found at the coefficients
in hand rather than at a refitted null, so the other terms' fits enter
it.

It comes into play IF AND ONLY IF the model carries a smooth penalty.
Where nothing is estimable – an ordinary `y ~ x`, or a model whose only
penalty is kinked – it is simply not run, and that is a property of the
model rather than of how the argument was written, so typing the default
changes nothing. `outer_criterion = NULL` holds every smooth
hyperparameter where its term left it.

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
#> fitted in 23 ms, converged

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
#> fitted in 28 ms, converged
```
