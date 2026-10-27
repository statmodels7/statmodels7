# Fit a Model

Fits a distributional regression. One formula carries an equation per
parameter of the response distribution, so a scale or a shape is modeled
as readily as a mean; the terms that those equations name are assembled
into a penalized likelihood and fitted, and any smoothing parameters or
prior scales are estimated from the data.

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
  inner_information = auto_info(),
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

  A data frame. Missing values in the response or in a variable that an
  equation uses are rejected, with the variable and the rows.

- weights:

  Optional prior weights, a numeric vector of length \\n\\ or the name
  of a column of `data`. The weights must be non-negative, finite and
  complete, and are taken as given and not normalized. They multiply the
  log-likelihood, so they act as frequency weights: a constant weight
  \\c\\ on every row is \\c\\ copies of the data, and the restricted
  (REML) criterion counts \\n c\\ observations against the coefficients.
  The REML estimate of a dispersion therefore changes under a constant
  rescaling of the weights, and the ML estimate does not.

- offsets:

  Optional named list of offsets, one entry per distribution parameter,
  each of length 1 or \\n\\. An offset written in the formula is added
  to the entry for the same parameter.

- inner_optimizer:

  How the smooth block is fitted:
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  (the default) or an optimizers7 optimizer. Beside a structural term of
  the filter shape, the coefficients and the term's parameters are
  fitted in one system, for which
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  stands for the system's own default optimizer and an optimizers7
  optimizer is used as given.

- inner_information:

  The matrix of second derivatives that the inner fit uses, an object
  that inherits from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html):
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
  (the observed Hessian),
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  (the Fisher information, whose `method` is `NULL` for the closed form
  of the family or
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
  or one of
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
  [`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
  and
  [`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html);
  see
  [distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html))
  or
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  (the default), which is settled against the family. With
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  the matrix is inverted at each step (Newton for the observed
  information, Fisher scoring for the expected one); with an optimizers7
  optimizer it is the Hessian that the optimizer receives, so that
  `inner_optimizer = optimizers7::newton()` with
  `inner_information = distributions7::expected_info()` takes the
  expected information. A named matrix is the only one that the inner
  fit uses.
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  chooses the expected information where the family writes it exactly
  and at an acceptable cost, and the observed information otherwise;
  with
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  it also lets the run change matrix and lets the expected information
  take a step that the observed one cannot (see
  [information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)).
  The fitted model records the settled object, and every later refit of
  the model (the folds of a cross-validation, the restricted fits of a
  test, the resampled fits of a bootstrap) reads the same one. Beside a
  score-driven filter
  ([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html))
  the coefficients and the filter's own parameters are fitted in one
  system, whose Hessian is the observed information, or the expected
  information conditional on the past where
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  is named here (Fisher scoring over that system);
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  leaves that system on the observed information. A model with a term of
  the likelihood shape
  ([`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html),
  a marginal break-point term) has no expected information, and
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  is an error there. The outer criterion has its own `information` (see
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md))
  and the standard errors have theirs (see
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md));
  the three choices are independent.

- outer_criterion:

  How the smooth hyperparameters are estimated:
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  (the default),
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  or
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  or `NULL` to hold them where their terms left them.
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  cross-validates the hyperparameter of a kinked penalty only: it is
  passed as `sparse_criterion`, and as `outer_criterion` it is an error.
  A prediction-error criterion beside a kinked penalty is rejected
  unless `sparse_criterion` is `NULL`.

- sparse_criterion:

  How the hyperparameter of a kinked penalty (lasso, elastic net, scad,
  mcp) is chosen:
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  (the default),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
  or `NULL`. With `NULL` the hyperparameter is held at its value, unless
  the model also carries a smooth penalty and `outer_criterion` is a
  prediction-error criterion, which then sweeps the path. A marginal
  criterion is rejected here, because it is read at a mode that sits on
  the kink.

- outer_optimizer:

  The optimizer that searches over the smooth hyperparameters, or `NULL`
  to choose one from the derivatives of the criterion that exist:
  [`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
  where the exact Hessian is usable,
  [`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
  with the exact gradient,
  [`optimizers7::brent()`](https://statmodels7.github.io/optimizers7/reference/brent.html)
  for one coordinate without a gradient, and
  [`optimizers7::nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html)
  otherwise.

- start:

  Where the fit begins: a named list of coefficients, a
  [`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
  such as
  [`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
  or
  [`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md),
  which takes the estimates of a model already fitted, or `NULL` for
  [`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md).
  A strategy runs once, before the alternation between the coefficients
  and the hyperparameters begins. A global search therefore belongs here
  and not in `inner_optimizer`, where it would run again at every
  hyperparameter value that the outer search tries.

- linpar_control:

  How the unpenalized parametric block is built, as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns it: the storage and the contrasts for its factors. It reaches
  the implicit term, the one into which the bare covariates collapse and
  which a caller never writes; a
  [`modelterms7::linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html)
  written out takes them directly. The argument and the function are
  named differently on purpose, as
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html)'s `control` and
  [`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html)
  are.

- verbose:

  A level from 0 to 3, or a named logical vector with any of `outer`,
  `blocks`, `inner` and `optimizer`.

- threads:

  How many threads the fit may use, as
  [`numericals7::n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.html)
  constructs it. The default, `n_threads(1)`, is sequential and takes
  the sequential code path. A larger count is passed to the compiled
  per-observation kernels of the family and to the dense assembly
  products; below a kernel's internal size threshold the computation
  stays sequential whatever the count, and a sparse design keeps its
  Matrix route. The object's `workers` fans the units that are
  independent by construction out over separate R processes: the folds
  of a cross-validation, and the combinations of a kinked path's product
  grid, each of which restarts its warm chain from the sweep's own
  starting coefficients. Each unit fits sequentially, so the two levels
  never nest. The points within one chain stay sequential, because a
  point fitted from a cold start costs more than one warm-started from
  the previous point, and splitting a chain would slow the
  single-process default or make the result depend on the count. For any
  two counts above one the result is identical bit for bit: every
  parallel region decomposes its work over the elements of its output
  and never splits a reduction, and the folds are seeded per fold and
  collected in fold order. Between the default and a count above one the
  route can change where the work is large (the BLAS crossproduct
  against the compiled one, and the pivoted
  [`qr()`](https://rdrr.io/r/base/qr.html) against the compiled
  triangular factor), and the results then agree to rounding.

- ...:

  Not used. A named argument signals an error. A hyperparameter is held
  in the term that carries the penalty.

## Value

A
[`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
object.

## The formula

The equations are separated by `|`, and the first carries the response:

        y ~ x1 + ridge(R) + lasso(L)  |  sigma ~ z  |  nu ~ 1

A parameter with no equation of its own is given an intercept, so
`y ~ x` on a two-parameter family fits a constant scale.
[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
does the split. R's operator precedence nests the equations on the left,
so the right-hand side of the whole formula holds only the last
equation, and the split walks down the left side instead.

The terms available are modelterms7's:
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html),
[`modelterms7::te()`](https://statmodels7.github.io/modelterms7/reference/te.html),
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html),
[`modelterms7::ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.html),
[`modelterms7::lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.html),
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
and the rest. Bare covariates collapse into one parametric block. A
formula can carry one structural term, which is
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
or the marginal break-point term, or one
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
term in each of several equations, which are filtered jointly as one
recursion.

A structural term whose level adds a constant to the predictor, the
\\\omega\\ of
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
or the first level of
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html),
is confounded with an intercept in the same equation. The intercept is
kept and the level held at zero, reported as `omega (held)` for
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
and as `level1 (held)` for
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html);
writing `0 + gas(...)` estimates the level instead. The two spellings
are the same model with the same log-likelihood, and the constants
differ because they measure different quantities: for a score-driven
term the intercept is the stationary level \\\omega / (1 - \sum_j
\xi_j)\\, which is \\\omega\\ divided by one minus the sum of the
persistences \\\xi_j\\. The page of
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
gives the identity.

## The fitting scheme

The terms split in two by a property that each one already reports.

A term whose penalty is twice differentiable in its coefficients, an
unpenalized block, a ridge, a spline or a random effect, is estimated in
one system by `inner_optimizer`. Their joint curvature exists, and using
it closes a fit in a handful of iterations.

A term whose penalty has a kink, a lasso, an elastic net, a SCAD or an
MCP, is estimated as a block of its own while the other blocks are held
fixed. The block is fitted by a compiled coordinate descent where the
penalty describes its proximal operator as a table at the step taken,
and otherwise by an accelerated proximal gradient. Where the equation
has a free unpenalized intercept, the descent centers the block and
moves that intercept jointly with it.

A structural term of the filter shape,
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
is fitted in one system with the coefficients. A structural term of the
likelihood shape,
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
or the marginal break-point term, is alternated with them.

The fit alternates between the blocks until the relative change of the
objective falls below the tolerance. Where a kinked block or a
structural term is present, the point reached must in addition satisfy
the optimality conditions of the smooth part and of every coefficient
held at a kink, and `converged` reports that reading.

## The objective is not averaged

The objective is minus the weighted log-likelihood plus the penalties at
full size. A penalty is a negative log-prior, and a posterior adds a
log-likelihood and a log-prior at full size; averaging one of them would
make a hyperparameter mean something that depends on \\n\\. The stopping
rule is scaled instead, so a threshold means the same at \\n = 10\\ and
at \\n = 10^7\\.

The weights multiply the log-likelihood, so they act as frequency
weights: a constant weight \\c\\ on every row is \\c\\ copies of the
data. The restricted (REML) criterion then counts \\n c\\ observations
against the \\p\\ coefficients, so the REML estimate of a dispersion
changes under a constant rescaling of the weights, and the ML estimate
does not.

## The budget and the stopping rule belong to the method

There is no `maxit` and no `tol` here. Both are set on
`inner_optimizer`, which is
[`iwls(maxit =, tol =)`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
or an optimizers7 optimizer with its own `maxit` and `criterion`. The
alternation reads them from there as its budget in passes and its
relative tolerance on the objective, so that one setting governs both
the inner fit and the alternation.

## Every hyperparameter is estimated unless its term holds it

The term that carries the penalty states which hyperparameters are held:
`lasso(x, lambda = 3)`, `ridge(x, lambda = 0.5)`,
`s(x, bspline_smooth(), hyper = c(lambda = 2))`, `enet(x, alpha = 0.5)`.
Every hyperparameter left `NULL`, which is each term's default, is
estimated from the data. A hyperparameter is held in the term and
nowhere else, so that the model has one statement of each held value.

The smooth hyperparameters go to `outer_criterion`,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
by default. `outer_optimizer` searches over them, and the coefficients
are refitted at each point.

A kinked penalty has its own argument. `sparse_criterion`,
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) by
default, sweeps its hyperparameters along a path of their own values,
from the kink that empties the block down to a fraction of it. The
penalized mode is only piecewise smooth in such a hyperparameter,
turning a corner whenever a coefficient joins the active set or leaves
it, so a gradient search on a criterion read there would meet those
corners.

The top of that path is the kink that empties the block, found at the
coefficients in hand and not at a refitted null. It therefore depends on
the data and on the rest of the model, since the fits of the other terms
enter it.

Where a model carries both kinds, the path is outside and the marginal
criterion is estimated inside each of its points, so a smoothing
parameter can come from REML and a lasso's \\\lambda\\ from BIC in one
fit.

With the default `marginal = NULL`,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
also estimate on the criterion the unpenalized coefficients of every
distribution parameter other than the first. `outer_criterion` therefore
runs when the model carries a smooth hyperparameter or such a
coefficient, and on an ordinary `y ~ x` with a two-parameter family it
estimates the dispersion by REML. Beside a kinked penalty, a structural
term or a block that moves with its coefficients, the default names no
coefficient, so a model whose only penalties are kinked does not run the
criterion. The default criterion is also not run when the family lacks
the two derivatives that it needs. With `outer_criterion = NULL`, every
smooth hyperparameter stays where its term left it and the coefficients
of all parameters are fitted at the joint mode, which gives the
maximum-likelihood estimate of a dispersion on an ordinary model.

## Verbosity

Three levels name the loops: `1` shows the outer search and the
alternation, `2` adds the inner method's own iterations, and `3` adds
the optimizers' traces. A named logical vector sets the four switches
`outer`, `blocks`, `inner` and `optimizer` independently, and a switch
that is left out is off. For example
`verbose = c(outer = TRUE, blocks = FALSE)` shows the hyperparameters
moving while a chatty inner optimizer stays silent.

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

# An ordinary regression: the scale receives an intercept that the
# formula does not specify.
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
coef(fit)
#> $mu
#> (Intercept)           x 
#>   0.9956298   2.0078791 
#> 
#> $sigma
#> (Intercept) 
#>  0.07592502 
#> 

# The scale is modeled too, which is the purpose of the framework. The data
# were drawn with log sigma = 0.3 x, and the interval covers it.
both <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
coef(both)$sigma
#>  (Intercept)            x 
#> -0.003093138  0.233104770 
confint(both)["sigma:x", c("estimate", "lower", "upper")]
#>          estimate     lower     upper
#> sigma:x 0.2331048 0.1427447 0.3234648

# It is the better model on this data.
c(one = AIC(fit), both = AIC(both))
#>      one     both 
#> 601.9454 578.9246 
```
