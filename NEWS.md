# statmodels7 0.9.0

* `cv()` chooses the hyperparameters of a penalty with a kink, which no
  criterion in the package could reach before. A marginal criterion
  approximates an integral at the penalized mode and asks for the
  penalty's second derivative there; the mode of a lasso, a SCAD or an
  MCP sits at the kink for every coefficient it sets to zero, which is
  where that derivative does not exist. Cross-validation asks about
  prediction instead and needs nothing from the penalty beyond a fit.
  `rule = "1se"` takes the largest kink whose criterion is within one
  standard error of the smallest.

* `aic()` and `bic()` reach those hyperparameters too, the effective
  degrees of freedom now being the trace over the coefficients that are
  away from the kink -- the submodel the fit selected, where the penalty
  is twice differentiable. For the lasso, which is linear there, that is
  exactly the number of surviving coefficients (Zou, Hastie and
  Tibshirani, 2007); the elastic net keeps its quadratic part and spends
  less, SCAD and MCP their own curvature and spend slightly more. The
  trace over the whole vector could not see the selection at all:
  measured on twenty noise columns it read 14 at every value of lambda,
  so a criterion built on it would have priced a model that selects
  nothing the same as one that selects everything.

* Neither is a gradient search. The penalized mode is piecewise smooth in
  the hyperparameter, differentiable while the active set holds and
  turning a corner whenever a coefficient joins it or leaves, so the
  value is swept over a grid. The grid is geometric in the size of the
  kink -- read from the penalty by probing its subdifferential rather
  than assumed to be lambda -- from the value that empties the block down
  to `min_ratio` of it, with every fit warm-started from the previous
  one. A choice at either end of the grid is reported: it is the grid's
  and not the criterion's.

  On 200 observations of twenty columns with three true predictors, every
  route keeps all three. `bic()` keeps four in all, `cv(rule = "1se")`
  six -- the same count as `cv.glmnet`'s own one-standard-error rule --
  and `cv(rule = "min")` sixteen against glmnet's nineteen.

* `predict()` and `loglik()` at new data reapply each term's recorded
  block instead of rebuilding it. A term records how its block was made,
  and `modelterms7::term_predict()` replays that record; rebuilding gives
  a block of the same shape multiplying the same coefficients that means
  something else. Measured on `y ~ s(x, k = 10)` at 200 observations,
  predicting on 40 of the rows the model was fitted to differed from the
  fitted values there by 0.237, and on the 51 rows with `|x| < 0.5`,
  where rebuilt knots move furthest, by 1.19. The whole data handed back
  agreed exactly, which is why nothing noticed, and cross-validation
  would have inherited it.

* The information matrix, its approximation and the budget are read off
  `inner_method` in one place, so a refit inside a path or a fold runs on
  the terms the caller asked for.

# statmodels7 0.8.1

* `scad()` and `mcp()` are fitted by the proximal scheme, as `lasso()`
  already was. The hyperparameters were passed to the penalty as a named
  numeric vector where the contract asks for a list, which `penalty_kinks()`
  accepts for some branches and rejects for these two; the failure was
  caught and read as "no kink", so both terms went into the jointly fitted
  system -- solved by the curvature of a function that has none. On a design
  of twenty noise columns the fit kept 19.00 effective degrees of freedom out
  of 20, which is no selection at all, and `y ~ scad(x)` stopped with
  "$ operator is invalid for atomic vectors" before reaching it.

  A penalty that stops when asked whether it has a kink is now reported,
  naming the term: the two schemes differ by exactly that property, so a
  penalty that does not answer cannot be assigned to one of them.

  The point the proximal scheme reaches satisfies the subdifferential
  conditions of the objective (stationarity 8e-8 where a coefficient is away
  from zero, and the unpenalized score inside the interval the kink opens
  where it is not), and the shrinkage answers the smoothing parameter: 17, 9
  and 0 of the twenty columns survive at lambda 1, 5 and 20. Against
  `ncvreg::ncvfit()` on the same objective the two land on different
  stationary points of a non-convex problem, ours the lower on both families
  (52.9966 against 53.5323 for SCAD, 52.9948 against 53.6638 for MCP), and
  starting the iteration from theirs moves 6e-2 back to ours.

# statmodels7 0.8.0

* Five defects a Student t fitted to `iris` exposed, each of which had
  been silent.

  `statmod_start()` read `distrib_start()`'s result by parameter name,
  where that result is a list of starts each keyed by parameter, so the
  name matched nothing and every fit began at zeros on the link scale --
  a location of 0 for a response centred at 5.84. The intercept of each
  equation now starts at the intercept-only maximum likelihood estimate,
  which for a gaussian is the sample mean and standard deviation
  exactly. That estimate draws its own starting values at random, so the
  stream is fixed for the length of the call and restored afterwards:
  without it the same call returned log-likelihoods of -103.49, -112.11
  and -111.83 on consecutive runs.

  `iwls()` stopped with "missing value where TRUE/FALSE needed" when the
  gradient at an accepted step was not finite, the line search having
  checked the objective alone, and `pd_repair()` raised from three
  frames further down on the same matrix. Such a point is now reported
  as unusable and the last good iterate is kept.

  The alternation set `converged` to `TRUE` as soon as there were no
  sparse blocks to alternate with, whatever the inner method had done,
  which is what let the three above go unnoticed. The verdict is the
  inner fit's.

  `vcov()` and `summary()` named two causes for a singular information,
  the fit not having reached a maximum and two columns carrying the same
  information, and on this fit neither applied: the design is full rank
  and the score is 4e-5. The message names the directions instead --
  which coefficients carry the flat eigenvector, or which rows are not
  finite.

* A term the fitting scheme does not cover is rejected when the
  specification is built, where it can be named, instead of reaching the
  design and raising on one of `term_matrix()` or `term_npar()`. Two
  shapes are outside the assembly, and both are read off the term rather
  than from a list of classes, so a term written later needs no edit
  here: a structural term, which rewrites the likelihood instead of
  contributing a predictor, and a term registering a `term_refresh()`
  method of its own, whose block moves with its coefficients.

  The second was the one that cost something. `seg()`, `jump()`,
  `jseg()` and `nl()` were assembled once and solved as though their
  block were a fixed design, so the break-point stayed at its starting
  value and `converged` was `TRUE`: measured against the same term
  iterated through `modelterms7::term_refresh()`, the break-point held
  at 2, 3, 5 and 8 from those four starts against 5.0681 from every one
  of them, at a residual sum of squares of 907.5, 478.1, 143.1 and 857.1
  against 143.1. Since the block of `seg()` is a Jacobian, the fit also
  used its break-point column as an ordinary regressor, giving a
  continuous construction a fitted mean with a step of -2.93.

# statmodels7 0.7.1

* The exact Hessian of a marginal or prediction-error criterion is now
  exact for a separable penalty too, penalties7 having stopped
  differencing its second-order pieces. A ridge and a random effect go
  through no difference at all.

* `outer_pieces()` takes an order and does not call the second-order
  generics at order 1, so a penalty that supplies only
  `penalty_dhessian()` gives an exact gradient rather than being
  rejected for a quantity the gradient does not use.

# statmodels7 0.7.0

* `aic()` and `bic()` choose the hyperparameters by an estimate of
  prediction error, `-2l + k*tau` with `tau = tr[(H+S)^-1 H]`, and both
  derivatives are exact. The envelope theorem does not apply to these:
  the log-likelihood alone is not stationary at the penalized mode, so
  its derivative carries `db/dt` from the first order, and what makes
  it computable is that at the mode `dl/db` equals `drho/db`.

# statmodels7 0.6.0

* The marginal criterion has an exact Hessian as well as an exact
  gradient, and both are modular: everything the penalty contributes is
  asked of the penalty through penalties7's generics, so a ridge, a
  random effect and a structured prior are covered by the same assembly
  as a spline. The route used to require a Hessian linear in the
  hyperparameters, measured rather than asked, which excluded every
  penalty built from a density.

* `outer_optimizer` defaults to `newton()` where the Hessian exists,
  `lbfgs()` where only the gradient does, and `nelder_mead()` otherwise.
  Measured in evaluations of the criterion: with three smoothing
  parameters, 31 against 41 against 269.

# statmodels7 0.5.0

* The marginal criterion has an exact gradient. The envelope theorem
  leaves only `-drho/dtheta` from the first two terms, and the
  determinant contributes `tr(M dS/dtheta)` and `u'v` with
  `v = db/dtheta` from the stationarity condition and
  `u_c = tr(M dK/db_c)` from the third derivative of the log-likelihood
  in the link-scale predictors. Checked against numDeriv at 1e-6 for
  one smoothing parameter, for several, with the scale modelled, and
  under `ml()`.

* `outer_optimizer` defaults to `lbfgs()` where the gradient exists and
  to `nelder_mead()` where it does not. Measured in evaluations of the
  criterion against the derivative-free search: 40 against 32 with one
  smoothing parameter, 40 against 135 with two, 41 against 269 with
  three.

# statmodels7 0.4.0

* The modelling layer. `statmod()` reads one formula carrying every
  parameter of a distribution, the equations separated by a bar, and
  fits it: the terms whose penalties are twice differentiable in one
  system by `iwls()` or any optimizer, each remaining block by a
  proximal method with the others held fixed, alternating until the
  objective stops moving.

* `reml()` and `ml()` estimate the hyperparameters by the Laplace
  approximation to the marginal likelihood, `outer_optimizer` searching
  over them. The criterion needs nothing added by hand: a penalties7
  penalty keeps its normalizing constant, so it is minus a log prior
  density, and written out the expression reproduces Wood's (2011) REML
  criterion term for term.

* `summary()` reads each distribution parameter as blocks -- the
  parametric terms together, then one per smooth, random effect or
  selection -- rather than as one list of coefficients. `vcov()`,
  `confint()`, `predict()` and `rstatmod()` come with it.

# statmodels7 0.3.0

* modelterms7 joins the toolkit as its eighth member: model terms as
  objects -- the parametric block, the penalized quartet, grouped
  random effects, the formula interpreter and the censored response.
  The member list is read from `Imports`, so nothing else changed.

# statmodels7 0.2.0

* numericals7 joins the toolkit as its sixth member and the root of its
  dependency graph: the numerical layer -- jets, and next the single stencil
  library and the parameter-vectorized quadrature -- written once where every
  package above can consume it. The member list is read from `Imports`, so
  nothing else changed.

# statmodels7 0.1.0

* First release. The package installs and attaches the member packages of the
  toolkit -- `linkfunctions7`, `distributions7`, `optimizers7`, `basis7` and
  `parameters7` -- so that installing it installs all of them and
  `library(statmodels7)` attaches all of them.

* `statmodels7_packages()` reads the members from this package's own
  `Imports`, keeping the names that end in `7`, so a member added there is a
  member everywhere and the list cannot disagree with itself.

* `statmodels7_versions()` reports the installed version of each member,
  `statmodels7_conflicts()` the exports that mask one another, and
  `statmodels7_update()` describes what is installed or installs it. The
  toolkit is not on CRAN, so the install path is GitHub through `pak`, which
  resolves the dependencies among the members itself.
