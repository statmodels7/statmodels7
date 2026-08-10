# Changelog

## statmodels7 0.6.0

- The marginal criterion has an exact Hessian as well as an exact
  gradient, and both are modular: everything the penalty contributes is
  asked of the penalty through penalties7’s generics, so a ridge, a
  random effect and a structured prior are covered by the same assembly
  as a spline. The route used to require a Hessian linear in the
  hyperparameters, measured rather than asked, which excluded every
  penalty built from a density.

- `outer_optimizer` defaults to `newton()` where the Hessian exists,
  `lbfgs()` where only the gradient does, and `nelder_mead()` otherwise.
  Measured in evaluations of the criterion: with three smoothing
  parameters, 31 against 41 against 269.

## statmodels7 0.5.0

- The marginal criterion has an exact gradient. The envelope theorem
  leaves only `-drho/dtheta` from the first two terms, and the
  determinant contributes `tr(M dS/dtheta)` and `u'v` with
  `v = db/dtheta` from the stationarity condition and
  `u_c = tr(M dK/db_c)` from the third derivative of the log-likelihood
  in the link-scale predictors. Checked against numDeriv at 1e-6 for one
  smoothing parameter, for several, with the scale modelled, and under
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).

- `outer_optimizer` defaults to `lbfgs()` where the gradient exists and
  to `nelder_mead()` where it does not. Measured in evaluations of the
  criterion against the derivative-free search: 40 against 32 with one
  smoothing parameter, 40 against 135 with two, 41 against 269 with
  three.

## statmodels7 0.4.0

- The modelling layer.
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  reads one formula carrying every parameter of a distribution, the
  equations separated by a bar, and fits it: the terms whose penalties
  are twice differentiable in one system by
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or any optimizer, each remaining block by a proximal method with the
  others held fixed, alternating until the objective stops moving.

- [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  and
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  estimate the hyperparameters by the Laplace approximation to the
  marginal likelihood, `outer_optimizer` searching over them. The
  criterion needs nothing added by hand: a penalties7 penalty keeps its
  normalizing constant, so it is minus a log prior density, and written
  out the expression reproduces Wood’s (2011) REML criterion term for
  term.

- [`summary()`](https://rdrr.io/r/base/summary.html) reads each
  distribution parameter as blocks – the parametric terms together, then
  one per smooth, random effect or selection – rather than as one list
  of coefficients. [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
  come with it.

## statmodels7 0.3.0

- modelterms7 joins the toolkit as its eighth member: model terms as
  objects – the parametric block, the penalized quartet, grouped random
  effects, the formula interpreter and the censored response. The member
  list is read from `Imports`, so nothing else changed.

## statmodels7 0.2.0

- numericals7 joins the toolkit as its sixth member and the root of its
  dependency graph: the numerical layer – jets, and next the single
  stencil library and the parameter-vectorized quadrature – written once
  where every package above can consume it. The member list is read from
  `Imports`, so nothing else changed.

## statmodels7 0.1.0

- First release. The package installs and attaches the member packages
  of the toolkit – `linkfunctions7`, `distributions7`, `optimizers7`,
  `basis7` and `parameters7` – so that installing it installs all of
  them and
  [`library(statmodels7)`](https://statmodels7.github.io/statmodels7/)
  attaches all of them.

- [`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
  reads the members from this package’s own `Imports`, keeping the names
  that end in `7`, so a member added there is a member everywhere and
  the list cannot disagree with itself.

- [`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
  reports the installed version of each member,
  [`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
  the exports that mask one another, and
  [`statmodels7_update()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)
  describes what is installed or installs it. The toolkit is not on
  CRAN, so the install path is GitHub through `pak`, which resolves the
  dependencies among the members itself.
