# Fit a Model

Fits a distributional regression. One formula carries an equation per
parameter of the response distribution, so a scale or a shape is modeled
as readily as a mean; the terms those equations name are assembled into
a penalized likelihood and fitted, and any smoothing parameters or prior
scales are estimated from the data.

Returns a
[`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md),
which
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
and the other accessors read.

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
  threads = numericals7::n_threads(),
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

  How the smooth hyperparameters are estimated:
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  (the default),
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
  or `NULL` to hold them where they are.

- sparse_criterion:

  How the hyperparameter of a kinked penalty – lasso, scad, mcp – is
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
  [`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
  such as
  [`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md),
  or `NULL` for
  [`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md).
  A strategy is asked once, before the alternation between the
  coefficients and the hyperparameters begins, which is why a global
  search belongs here instead of in `inner_optimizer`: there it would
  rerun at every hyperparameter the outer search tried.

- linpar_control:

  How the unpenalized parametric block is built, as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns it: the storage and the contrasts for its factors. It reaches
  the implicit term, the one the bare covariates collapse into and which
  a caller never writes; a
  [`modelterms7::linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html)
  written out takes them directly. The argument and the function are
  named differently on purpose, as
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html)'s `control` and
  [`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html)
  are.

- verbose:

  A level from 0 to 3, or a named logical vector.

- threads:

  How many threads the fit may use, as
  [`numericals7::n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.html)
  constructs it. The default, `n_threads(1)`, is sequential and takes
  exactly the sequential code path. A larger count reaches the compiled
  per-observation kernels of the family and the dense assembly products
  as an argument; below a kernel's measured internal threshold it stays
  sequential whatever the count says, and a sparse design keeps its
  Matrix route. The object's `workers` fans the units that are
  independent by construction out over separate R processes – the folds
  of a cross-validation, and the combinations of a kinked path's product
  grid, each of which restarts its warm chain from the sweep's own
  starting coefficients – each unit fitting sequentially, so the two
  levels never nest. The points within one chain stay sequential:
  measured, a point paid cold costs 2.2-3.2 times the warm chain, so
  splitting a chain would slow the single-process default or make the
  result depend on the count. The result does not depend on either
  count, bit for bit: every parallel region decomposes its work over the
  elements of its output and never splits a reduction, and the folds are
  seeded per fold and collected in fold order.

- ...:

  Not used, and reported. `hyper` was removed from here: a
  hyperparameter is held in the term that carries the penalty, and a
  second place to say so would be read by nobody whenever the two
  disagreed.

## Value

A
[`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
object.

## The formula

The equations are separated by `|`, and the first carries the response:

        y ~ x1 + ridge(R) + lasso(L)  |  sigma ~ z  |  nu ~ 1

A parameter with no equation of its own gets an intercept, so `y ~ x` on
a two-parameter family fits a constant scale.
[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
does the split; the recovery is not the obvious one, since R's
precedence makes the whole right-hand side of a three-equation formula
the last term alone.

The terms available are modelterms7's:
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html),
[`modelterms7::te()`](https://statmodels7.github.io/modelterms7/reference/te.html),
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html),
[`modelterms7::ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.html),
[`modelterms7::lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.html),
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
and the rest. Bare covariates collapse into one parametric block.

## The fitting scheme

The terms split in two by a property each one already reports.

A term whose penalty is **twice differentiable** in its coefficients, an
unpenalized block, a ridge, a spline or a random effect, is estimated in
one system by `inner_optimizer`. Their joint curvature exists, and using
it is what closes a fit in a handful of iterations.

A term whose penalty has a **kink**, a lasso, a SCAD or an MCP, is
estimated by a coordinate descent of its own with everything else held
fixed.

The fit alternates between the two until the objective and every block
stop moving.

## The objective is unaveraged

Minus the weighted log-likelihood, plus the penalties at full size. A
penalty is a negative log-prior, and a posterior adds a log-likelihood
and a log-prior at full size; averaging one of them would make a
hyperparameter mean something that depends on \\n\\. What is scaled
instead is the stopping rule, so a threshold means the same at \\n =
10\\ and at \\n = 10^7\\.

## The budget and the stopping rule belong to the method

There is no `maxit` and no `tol` here. Both are set on
`inner_optimizer`, which is
[`iwls(maxit =, tol =)`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
or an optimizers7 optimizer with its own `maxit` and `criterion`, and
the alternation reads them from there. A second copy here would let a
caller set both and be obeyed by neither.

## Every hyperparameter is estimated unless its term holds it

Which ones are held is said by the **term** that carries the penalty:
`lasso(x, lambda = 3)`, `ridge(x, sigma = 0.5)`, `s(x, lambda = 2)`,
`enet(x, alpha = 0.5)`. Everything left `NULL`, which is each term's
default, is chosen from the data. The term is where the penalty is named
and so is where that belongs; an argument here saying the same thing
would be read by nobody whenever the two disagreed.

The **smooth** hyperparameters go to `outer_criterion`,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
by default, with `outer_optimizer` searching over them and the
coefficients refitted at each point.

A **kinked** penalty is a different instrument with its own argument.
`sparse_criterion`,
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) by
default, sweeps it along a path of its own values, from the kink that
empties the block down to a fraction of it. The penalized mode is only
piecewise smooth in such a hyperparameter, turning a corner whenever a
coefficient joins the active set or leaves it, so a criterion read there
inherits the corners and a gradient search would read a slope about to
change.

The top of that path depends on the data **and on the rest of the
model**: it is the kink that empties the block, found at the
coefficients in hand never at a refitted null, so the other terms' fits
enter it.

Where a model carries both kinds the path is outside and the marginal
criterion is estimated inside each of its points, so a smoothing
parameter can come from REML and a lasso's \\\lambda\\ from BIC in one
fit.

`outer_criterion` runs if and only if the model carries a smooth
penalty. On an ordinary `y ~ x`, or on a model whose only penalty is
kinked, there is nothing for it to estimate and it is not run, so typing
the default changes nothing. `outer_criterion = NULL` holds every smooth
hyperparameter where its term left it.

## Verbosity

Three levels, naming the loops: `1` the outer search and the
alternation, `2` the inner method's own iterations, `3` the optimizers'
traces as well. A named form is accepted too,
`verbose = c(outer = TRUE, blocks = FALSE)`, for watching the
hyperparameters move while silencing a chatty inner optimizer.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for what to do with the result,
[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md)
to build the model without fitting it,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the inner method,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for the criteria,
[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
to simulate from a model.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(200, -2, 2))
dd$y <- 1 + 2 * dd$x + rnorm(200, sd = exp(0.3 * dd$x))

# An ordinary regression: the scale gets an intercept it was not given.
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
coef(fit)
#> $mu
#> (Intercept)           x 
#>   0.9956298   2.0078791 
#> 
#> $sigma
#> (Intercept) 
#>   0.0708999 
#> 

# The scale modeled too, which is what the framework is for. The data
# were drawn with log sigma = 0.3 x, and the interval covers it.
both <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
coef(both)$sigma
#>  (Intercept)            x 
#> -0.008201081  0.234254511 
confint(both)["sigma:x", c("estimate", "lower", "upper")]
#>          estimate     lower     upper
#> sigma:x 0.2342545 0.1429599 0.3255491

# It is the better model on this data.
c(one = AIC(fit), both = AIC(both))
#>      one     both 
#> 601.9354 578.9139 
```
