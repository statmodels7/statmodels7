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
