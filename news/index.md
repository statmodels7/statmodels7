# Changelog

## statmodels7 0.48.0

- A penalty may now answer with a sparse Hessian – a smooth repeated
  over the levels of a factor does – and the coercion is written ONCE,
  where the two kinds meet in
  [`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md),
  rather than at each of the eighteen places that read its result and
  the four that read `penalty_dhessian()`.

- The accumulator stays a base matrix, which is measured rather than
  timid. Making it the design’s own kind saved 0.8x end to end (17.0x
  against 16.2x), and it would have put twenty-two consumers on a
  contract only the two a sparse design happens to exercise are covered
  by. What the blocked penalty buys is its CONSTRUCTION and its own
  storage, and neither passes through this accumulator.

- Measured end to end on 120 levels and 3000 rows,
  `s(x, by = g, sparse = TRUE)` against the dense spelling: 26.17 s
  against 424.86 s and 0.55 MB against 20.44, with the log-likelihood
  identical to ten digits, the same lambda, the same effective degrees
  of freedom, [`vcov()`](https://rdrr.io/r/stats/vcov.html) agreeing
  exactly and [`predict()`](https://rdrr.io/r/stats/predict.html) to
  2.2e-16.

## statmodels7 0.47.0

- [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  takes `linpar_control`, as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns it: `sparse` for the storage of the unpenalized parametric
  block and `contrasts` for the coding of its factors. It reaches the
  IMPLICIT `linpar()` term, the one the bare covariates collapse into
  and which a caller never writes; a `linpar()` written out takes them
  directly.

- The specification CARRIES them, so a rebuild reproduces the storage: a
  fold of
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  that built a dense design where the fit built a sparse one would be
  paying for a storage the model did not ask for.

- Measured on 150 groups, the parametric block goes from 1.54 MB to 0.12
  MB with an identical log-likelihood and identical coefficients, and
  the fit is faster rather than slower.

- There is NO rescaling among the options, and that is measured rather
  than omitted. Scaling the columns and carrying the coefficients back
  is the remedy for a conditioning that squares, which is what forming
  `X'X` does;
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  fits through a QR and never forms it. On columns spanning fifteen
  decades the raw fit and the scaled one converge in the same number of
  iterations and both agree with
  [`lm()`](https://rdrr.io/r/stats/lm.html) to 1e-14. What moves is the
  SCORE the fit reports, 1.5e+02 against 9.2e-05, and the final verdict
  already arbitrates that on a dimensionless scale – so the argument
  would have changed a printed number and not an answer.

- The argument and the function are named differently on purpose. With
  one name for both, the argument’s default
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  resolves to the argument’s own promise and R reports “promise already
  under evaluation” – the shape section 7 records for an argument named
  after a class, met here for one named after a function.
  `glm(control = glm.control())` keeps them apart for the same reason.

## statmodels7 0.46.0

- A fold of
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  carries a term’s matrix input onto its own rows.
  `data.frame(X = X, y = y)` SPLITS a matrix into `X.x1 ... X.xp`,
  leaving no column `X`, so `lasso(X)` reaches past the data to the
  matrix in the calling environment – `interpret_formula()` evaluates
  the call as `eval(call, data, env)` and looks a name up in the data
  first, in the formula’s environment after. The fit was right, the
  matrix being captured once, and every fold then failed to rebuild
  because the name still resolved to all the rows.

- The matrix is already on the built term, and `term_build()` checked at
  the full fit that it has one row per observation, so a fold’s rows are
  the same rows by position.
  [`cv_bind_inputs()`](https://statmodels7.github.io/statmodels7/reference/cv_bind_inputs.md)
  binds the subset as a column of the fold, which builds for the fold
  the spelling the documentation asks the caller for. Measured, the two
  spellings now give identical coefficients, an identical criterion and
  an identical lambda on the same folds.

- Nothing is relearned that should be. A matrix carries no knots, no
  contrasts and no levels, so subsetting it and re-evaluating it give
  the same block; a FORMULA input is untouched and keeps being rebuilt
  on the fold’s own rows, which is what that rule exists for. A sparse
  input stays sparse – the column the fold is given is still a
  `dgCMatrix`, not its densification.

- The test fold is bound too: `term_predict()` evaluates a matrix
  input’s expression in the new data with
  [`baseenv()`](https://rdrr.io/r/base/environment.html) as its
  enclosure, so a name that is not a column there is not found at all.

## statmodels7 0.45.0

- `n_values` and `min_ratio` leave
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  and `OuterMethod`, as `search` did: what a PATH does is not the
  criterion’s business. The same criterion is put to every
  hyperparameter of a model, and a smooth one is read at the mode rather
  than swept, so three of the criterion’s arguments described something
  most of what it was asked about does not have.

- They live on the term’s own signature, where a reader can see the
  number: `lasso(x, n_lambda = 25, min_ratio = 1e-4)`,
  `enet(x, n_lambda = 25, n_alpha = 5)`.
  [`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
  is what remains here, and it is reached only by a term that declares a
  kinked penalty without offering an argument for the grid – `random()`
  under a Laplace prior is the case.

## statmodels7 0.44.0

- `search` leaves
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md).
  Whether a term’s own kinked hyperparameters are covered by a product
  or one at a time is the TERM’s, `enet(X, search = "cyclic")`, for the
  reason the whole enumeration is: a criterion is put to every
  hyperparameter of the model, and a smooth one is read at the mode
  rather than swept, so most of what it was asked about could not use
  the argument.
  [`statmod_search()`](https://statmodels7.github.io/statmodels7/reference/statmod_search.md)
  reads it, beside
  [`statmod_grid_size()`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md)
  and
  [`statmod_min_ratio()`](https://statmodels7.github.io/statmodels7/reference/statmod_min_ratio.md).

- Per term is also what keeps one term’s choice off another’s:
  `y ~ lasso(X) + enet(R, search = "cyclic")` sweeps the elastic net one
  coordinate at a time and leaves the lasso alone.

## statmodels7 0.43.0

- A hyperparameter the readable block does not DESCRIBE keeps its own
  row. The quantities of a multivariate Student t prior are the standard
  deviations and the correlations of its scale matrix, and its degrees
  of freedom are none of those, so replacing the coordinate rows
  wholesale dropped `nu` from the summary. The question is put to the
  Jacobian – a column that is zero throughout is a coordinate no
  quantity depends on – so a family that declares more later is covered
  without an edit.

- Measured across effect distributions that carry a shape: `skewnormal1`
  reports `sigma` and `alpha` with standard errors and intervals,
  `student_t1` reports `sigma` and `nu`, and a log-transformed gamma
  reports the variance of the gamma underneath at 0.2202 (se 0.0396,
  interval 0.155 to 0.313) against the 1/a = 0.25 the effects were drawn
  from. Where an interval is absent the cause is the point and not a
  missing derivative: a shape escaping towards a limit leaves an outer
  curvature of the wrong sign, and no interval follows from it.

## statmodels7 0.42.0

- [`summary()`](https://rdrr.io/r/base/summary.html) reports a
  correlated random effect by the quantities it is about. Where a
  penalty answers
  [`penalties7::penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.html)
  the coordinate rows are replaced by the standard deviations and
  correlations of the effects, with a standard error from the delta
  method – composing the penalty’s Jacobian with the link’s, the
  variance matrix being on the free scale the criterion was maximized on
  – and each interval built on the scale the quantity declares and
  mapped back: log for a standard deviation, Fisher’s z for a
  correlation. A standard deviation therefore cannot be given a negative
  lower end and a correlation cannot be given one that leaves (-1, 1),
  which is the rule every other interval in the toolkit follows. No test
  is printed, the null a z would report on being that the quantity is
  zero, which for a standard deviation is the edge of its range.

- Measured against `lme4` on random slopes over 40 groups, the effects’
  standard deviations and correlation are (1.2370, 0.4820, -0.1399)
  against (1.2395, 0.4898, -0.1373).

## statmodels7 0.41.0

- Nothing in the fitting layer changed: a random effect whose prior is
  read blockwise is one penalty entry like any other, and one per
  within-group column is the enumeration this layer has run on since
  0.12.0. What changed is what a fit reports – the hyperparameter of a
  gaussian random effect is now the standard deviation of the effects
  rather than a precision, so the tests that constructed or interpreted
  its value were turned round with it.

- Measured against `lme4` on a random intercept over 40 groups: the
  fixed effects agree to 1.3e-05 and the standard deviation of the
  effects is 1.14197 against 1.143856. With random slopes the fixed
  effects agree to 1.0e-03 and the effects’ standard deviations and
  correlation are (1.239, 0.568, 0.448) against (1.242, 0.575, 0.442).

## statmodels7 0.40.0

- [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  take `search`. A term carrying several hyperparameters with a kink has
  every combination of them visited under `"grid"`, the default, and one
  swept at a time under `"cyclic"`. Between two terms the search
  alternates either way, so `y ~ lasso(X) + enet(R)` costs the two
  blocks added and not multiplied, and two elastic nets in one formula
  are 100 + 100 points rather than 10^4.

- The product is the default because a cyclic sweep traverses a cross
  through the point in hand and can stop where each coordinate is
  separately best without being jointly so. It costs
  `n_lambda * n_alpha` fits against `n_lambda + n_alpha` per pass; with
  three or more estimated it grows exponentially, which is what
  `"cyclic"` is for.

- The grid is not a rectangle, and each axis is built at the settings of
  the axes outside it rather than fixed in advance. For the elastic net
  the kink is `lambda*alpha`, so every alpha carries its own lambda axis
  descending from its own `lambda_max = kink/alpha`; for SCAD and MCP
  the shape leaves the kink alone, so `lambda_max` is one number
  whatever the shape and those two axes really are a rectangle. Both
  follow from carrying the size of the kink back onto the hyperparameter
  where it stands, and a test pins the second, which is the controproof
  against writing the elastic net’s relation for a family that does not
  obey it.

- A hyperparameter a term wrote out as a vector is visited as it stands.
  The value emptying the block does not cap it and `min_ratio` does not
  extend it, both of those being ways to build a grid, and the order it
  is walked in is the penalty’s – from the emptiest fit towards the
  fullest – so the warm starts run as they do on a built grid.

- A shape parameter is swept above the smallest value at which the block
  can be FITTED. SCAD’s proximal operator needs `t < a - 1` and MCP’s
  `t < gamma`, tightened to `t*d^2` under standardization’s diagonal
  map, and `t = 1/sum(w x^2)` is a property of the data: measured, a
  standardized penalty on a column of spread 20 needs `a > 3` where
  SCAD’s own bound is 2, and a Poisson block whose fitted means are near
  1e-3 needs `a > 11`. The limit is asked of the penalty and bisected
  rather than written out, so neither constant appears in the path.

- The size of the kink is inverted in closed form. Measured over four
  decades, it is exactly a power of each hyperparameter – one for the
  lasso and for the elastic net in both of its own, minus one for a
  Laplace prior written by its scale, zero for the shapes of SCAD and
  MCP – with a spread of at most 5.6e-16, so the exponent is read once
  per grid and every value follows. The relation is measured and then
  checked at the values it produced, and a penalty that does not obey
  one falls back to bracketing. Measured, a bracketing solve cost 4.18
  ms against a fit’s 62.5 ms, so a path of twenty-five points spent 6.7
  per cent of itself locating the values it would visit.

- A pass that would visit the points just scored is not run. The top of
  the path is refreshed between passes because the rest of the model
  moves, and where it has not the grid is the one already in hand.

- `history$outer` gains `setting`, the rest of the combination each
  point belongs to. One row is still one point, and `name` and `value`
  still carry the axis the path descends.

## statmodels7 0.39.0

- The top of a path is read again at every sweep. It is the size of the
  kink that empties the block, taken from the score at the coefficients
  in hand, so it moves with the rest of the model: a smooth beside a
  kinked term has its own smoothing parameter estimated INSIDE each
  point of the path, and every such fit changes the score it is read
  from. Reading it once left the path anchored to the state the search
  began in. Measured on an elastic net, the largest lambda of the path
  moves from 37.0307 to 46.2884 between the first sweep and the second.

  It bites where there is more than one sweep, which today means more
  than one hyperparameter under estimation; a single kinked
  hyperparameter still makes one pass and reads the top once.

## statmodels7 0.38.0

- The path visits as many values as the TERM asked for, per
  hyperparameter: `lasso(~x, n_lambda = 50)`,
  `enet(~x, n_lambda = 40, n_alpha = 12)`. How finely a hyperparameter
  is swept is a property of the term for the same reason as whether it
  is swept at all, and the criterion applies to every term of the model
  at once and cannot know which it is looking at.
  [`statmod_grid_size()`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md)
  reads it from the penalty’s entry and falls back to the criterion’s
  `n_values` where the term said nothing.

- And as far DOWN as it asked: `lasso(~x, min_ratio = 1e-6)` sets the
  fraction of the emptying kink the path descends to. One number per
  term, since only the sweep by kink size uses it.
  [`statmod_min_ratio()`](https://statmodels7.github.io/statmodels7/reference/statmod_min_ratio.md)
  reads it and falls back to the criterion’s.

## statmodels7 0.36.0

- Which hyperparameters are estimated is said by the TERMS, and by
  nothing else. Every one is estimated unless the term that carries the
  penalty holds it: `lasso(x, lambda = 3)`, `ridge(x, sigma = 0.5)`,
  `enet(x, alpha = 0.5)`, `scad(x, a = 3.7)`, `s(x, lambda = 2)`,
  `te(x, z, lambda = c(1, 5))`,
  `random(~1 | g, hyper = c(sigma = 0.4))`. The term is where the
  penalty is named, so it is where that belongs.

- `statmod(hyper = )` is REMOVED, and so is `over` on
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md).
  Both said the same thing as the terms, and whichever of the two a
  reader believed, the other was read by nobody when they disagreed.
  Passing `hyper` now signals an error naming the spelling that works.

- A hyperparameter no path could reach is swept over a grid of its own.
  The elastic net’s `alpha` scales the kink and is bounded by one, so no
  admissible value of it empties the block; the shape of SCAD and MCP
  has no upper bound and does not move the kink at all. Neither is
  reachable by the geometric path over kink sizes, and leaving them out
  was how `alpha` stayed at 0.5 while the summary said a criterion had
  chosen it. A bounded one is swept over its interval and an unbounded
  shape over a geometric grid above its lower bound, which spans the
  values the literature uses (3.7 for the SCAD of Fan and Li, 3 for the
  MCP of Zhang).

- The summary reads what was held from the terms, so a hyperparameter is
  marked fixed if and only if a term fixed it.

## statmodels7 0.35.0

- A bounded hyperparameter is reported as held, not as chosen.
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) returns a result the
  length of its TEST, and the test was a scalar, so a penalty carrying
  two hyperparameters got one answer recycled over both: the elastic
  net’s `alpha`, which no path varies, was marked as chosen by the
  criterion that had chosen its `lambda`.

- [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  take `over`, as
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  already did, and a bounded hyperparameter named there is now really
  swept. It used to be accepted and ignored: a path walks the SIZE OF
  THE KINK, from the value that empties the block down, and no
  admissible `alpha` empties it at a given `lambda`, so every point was
  dropped and `alpha` came back at its default. A bounded hyperparameter
  is swept over its own interval instead, the endpoints excluded because
  the bounds are open – an elastic net at `alpha = 0` has no kink at all
  – and the sweeps being cyclic, the kink still moves through `lambda`.
  Measured on ten columns of which seven carried signal,
  `bic(over = c("lambda", "alpha"))` moves alpha off 0.5 to 0.615 and
  the criterion from 322.46 to 320.38.

  The default still holds it, as `glmnet` holds its `alpha` and `ncvreg`
  its `gamma`, and the end-of-path warning is not raised for such a
  sweep: the grid covers the whole interval, so an answer at either end
  is the interval’s and there is nothing to widen.

## statmodels7 0.34.0

- A hyperparameter chosen by `sparse_criterion` is reported as estimated
  by the criterion that chose it.
  [`summary()`](https://rdrr.io/r/base/summary.html) asked
  `methods$outer`, which carries the marginal criterion alone, so a
  lambda a path had selected was printed `(fixed)` beside a note saying
  it was held at the value it was given – of a number the caller had
  never seen. The fit now records which criterion swept the kinked
  penalties and exactly which hyperparameters it reached, and the cell
  says `(bic)`, `(cv)` or `(fixed)` accordingly.

- The hyperparameters come first in every penalized block. They govern
  every coefficient under them, and a table opening with a hundred
  selected coefficients buried the one number that produced the
  selection.

- A hyperparameter estimated by
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  or
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  carries a standard error and an interval. It maximizes a criterion
  that is twice differentiable in it, so its variance is the inverse of
  the negative of that criterion’s Hessian at the point reached, which
  [`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
  already computed exactly;
  [`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
  returns it. The interval is built on the free scale its link defines
  and mapped back, as every other interval in the toolkit is, so a
  positive hyperparameter keeps a positive lower end. No test is printed
  beside it: the null a `z` would report on is that the hyperparameter
  is zero, which is the edge of its range rather than an interior
  hypothesis. Validated against the criterion itself, refitted at each
  probe and differenced twice on the free scale.

  A hyperparameter chosen by a path –
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  over a kinked penalty – still carries none, and that is a refusal
  rather than an omission: it is the argument of a minimum over a grid,
  not the root of a derivative, so there is no curvature to read and its
  uncertainty is a resampling question.

- A path that cannot score a single point signals an error instead of
  keeping the hyperparameter it came in with.
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  refits on each fold, and a term built from a matrix that lives in the
  calling environment rather than in `data` cannot be re-evaluated on a
  subset of it, so every fold failed, every deviance came back `NA`,
  nothing was selected, and the fit returned the DEFAULT hyperparameter
  reporting success – a lasso at `lambda = 1` selecting every column.
  The error carries the message the refit failed with. A marginal
  criterion reaching a kinked row is unaffected: it scores nothing by
  construction, and the hyperparameter keeping its given value is the
  documented answer there.

## statmodels7 0.33.0

- A break-point term has a section of its own in
  [`summary()`](https://rdrr.io/r/base/summary.html), beside the
  smooths, the random effects and the penalized blocks, and it reports
  the quantities of the model rather than the coefficients of the
  working block: the linear effect, the changes, and the break-points
  under the name `psi`. The standard error of a position comes from
  [`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html)’s
  Jacobian by the delta method.

  A break-point gets an estimate, a standard error and an interval, and
  no test. The null a `z` of estimate over standard error reports on is
  that the position is zero, which is not a hypothesis anyone holds; the
  one a reader wants is that there is no break-point at all, and under
  it the position is a nuisance parameter that vanishes, so the
  classical p-value is wrong by a factor of three to five. `segmented`
  prints the estimate and the standard error alone for the same reason.

- The default start asks each term where its own block begins, through
  [`modelterms7::term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.html).
  It used to put every coefficient at zero except an equation’s
  intercept, and for a term that recomputes its block from its
  coefficients zero is a singular point rather than a neutral one:
  `y ~ jump(x, npsi = 3)` returned every coefficient exactly zero and
  reported failure, the three break-points having collapsed onto one
  clamped position. A term asking for zeros is left alone, so an
  ordinary block is unchanged, and
  [`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md)
  still means the origin.

## statmodels7 0.32.0

- The exact gradient of the marginal criterion reaches a penalty over a
  STRUCTURAL term’s own parameters.
  [`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
  no longer refuses `isTRUE(u$structural)`; it asks whether the term
  answers
  [`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html),
  read from the class the method is registered on, so a term written
  later is covered without a list of class names here.
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)/[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  on a score-driven panel therefore run `lbfgs()` on an exact derivative
  where they ran `nelder_mead()` on values.

  Where a filter is present the determinant spans the term’s parameters
  too, so `K` carries `-sum_t w_t l_p E_t` and differentiating it along
  the direction the mode moves in gives three contributions instead of
  one:
  [`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)’s
  formula with V in place of X (direction-free, computed once), the
  derivative of `V_p` itself, which is `E_t v`, and the derivative of
  `l_p E_t`, whose second half is the term’s third derivative. The full
  third derivative is never formed – it is contracted in one direction,
  so the cost is O(n m^2) per hyperparameter rather than O(n m^3) once.

  Validated against numDeriv on the criterion with the mode REFITTED at
  every hyperparameter: 2.7e-09 to 2.7e-07 relative under
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  3.3e-06 under
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  6.3e-06 with a smooth beside it, and 1.2e-07 at ten and at twenty
  groups. The control, an ordinary smooth, sits at 3.9e-07 on the same
  harness.

  Measured in evaluations of the criterion, each a whole inner fit,
  against `nelder_mead()`, all three reaching the same criterion and the
  same prior scale to the printed digit:

  ``` R
    3 groups x 40      6 against 59      4.1 s against 14.8 s
    8 groups x 80     10 against 27     18.7 s against 33.2 s
   20 groups x 60      9 against 41     88.0 s against 88.9 s
  ```

  The last row is worth reading as it is: the evaluation count is what
  the exact route buys and it is 4.6x there, while the WALL TIME is a
  wash, because each exact evaluation also pays for the gradient (one
  `term_third()` and one extra `term_curvature()` per hyperparameter)
  and one point on that path is still unavailable and costs a long inner
  fit that does not converge.

- An outer step the search TAKES BACK now leaves no trace, where until
  now only half of it did. The coefficients were already protected –
  `state$beta` is written after a point is known to be usable – but a
  structural term’s own parameters live in the design’s structural
  state, which is an environment the inner fit writes into as it goes,
  and
  [`statmod_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md)
  stores its optimizer’s last point whether that optimizer converged or
  not. An unavailable point therefore MOVED the filter’s parameters, and
  the next evaluation started from wherever the failure had left them.

  It ratchets, and the exact gradient is what made it reachable: with a
  derivative-free search the steps are small and the state never leaves
  the basin, while `lbfgs()` takes its first step as the full gradient,
  which on a panel of twenty groups is -12.3 on a log-scale
  hyperparameter (the criterion’s derivative grows with the number of
  penalized coordinates), so the first trial prior scale is 4.6e-06.
  Measured before the fix: the search then reported THIRTY consecutive
  unavailable points as the line search halved eta back to 1e-08,
  including points beside a start that had converged at the first
  evaluation, and gave up at 620 s with no fit. After it: nine recorded
  evaluations (ten attempts, the one unavailable point recovered from),
  and the same answer `nelder_mead()` reaches in 41 – prior scale
  0.3783702 against 0.37837, criterion -1696.89343 both.

  ⚠️ Two explanations were measured and refused before this one. The
  criterion is NOT noisy here: repeated at one hyperparameter from a
  running warm start, and arrived at from eight different previous
  hyperparameters as far as six units away, its spread is exactly zero.
  And the warm start is not poisoned by VISITING the degenerate point:
  stepping to eta = -12.3 and back returns the identical criterion. What
  ratchets is the structural state left behind by a failed inner search,
  which no snapshot protected.

- A filter beside a RANDOM EFFECT fits.
  [`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
  lays each equation’s design into its own columns of a row spanning
  every unknown, and an equation carrying a random effect has a sparse
  design: writing a `dgCMatrix` into a slice of a base matrix is a
  LENGTH ERROR rather than a conversion, so
  `y ~ x + random(~1|id) + gas(...)` stopped with “number of items to
  replace is not a multiple of replacement length”. It is the sixth
  round of the shape the sparse conversion keeps taking, and it was NOT
  reached by the exact gradient – the failure is in the joint inner
  step, which has called that assembly since 0.23.0, so the combination
  has never fitted. The block is densified at that one assignment and
  nowhere else: the matrix is one column per unknown with a single block
  filled, it is read one ROW at a time by the callback beneath it (the
  access a compressed-column matrix is worst at), and the fit, the
  penalties and the solve all go on seeing the sparse design.

- [`deriv4_key()`](https://statmodels7.github.io/statmodels7/reference/deriv4_key.md),
  [`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md),
  [`structural_joint_basis()`](https://statmodels7.github.io/statmodels7/reference/structural_joint_basis.md).
  The last is written once and used by both
  [`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
  and the gradient: two callers composing the integrated subspace
  separately would agree only by accident.

## statmodels7 0.31.0

- [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  and
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  default to the OBSERVED information (Giovanni). That is what turns on
  the exact outer gradient and Hessian, so the search becomes `newton()`
  instead of `nelder_mead()`. Measured in evaluations of the criterion,
  each a whole inner fit: 19 against 31 with one hyperparameter – where
  it does not pay, a simplex needing no derivative in one dimension –
  then 126 against 35 with a smooth and a random effect, 133 against 32
  with two smooths, and 166 against 6 with a modelled scale, the
  criterion agreeing to four decimals throughout.

  ⚠️ It has a cost, and it is on the criterion’s PLATEAU. Once the
  penalty dominates the information the criterion is flat – its kept
  normalizing constant cancels the Laplace determinant – and with the
  exact Hessian Newton walks the whole of it: on a response scaled by
  1e4 the smoothing parameter reaches 9.5e-20 where the simplex stopped
  at 3.2e-8, the fitted function identical to five decimals. The fit is
  unharmed, the reported hyperparameter is not meaningful there, and
  [`modelterms7::edf()`](https://statmodels7.github.io/modelterms7/reference/edf.html)’s
  plain [`solve()`](https://rdrr.io/r/base/solve.html) cannot invert
  `H + S` that far out. Three tests now ask for `reml("expected")` where
  what they check is their own reference or the inner rule rather than
  the outer surface.

  This does not reach a penalty over a STRUCTURAL term’s own parameters,
  where the exact-gradient route rejects for its own documented reason:
  the contraction assumes the predictor is `X beta` and a filter’s level
  is a recursion. Such a search is derivative-free whatever the
  information.

- The verbose trace is formatted. The method that will do the work is
  named on the right of every rule, blocks are separated by blank lines,
  and the hyperparameters’ names – which are long, a term’s key being
  the call that produced it – are said ONCE in a legend while the lines
  carry the values. A sparse block’s line names the route that actually
  ran, coordinate descent or proximal gradient, rather than one guessed
  at beforehand: whether the compiled descent applies depends on the
  operator being admissible at the step the fit takes.

- An outer point where the criterion is unavailable is now a FINITE
  barrier, strictly worse than everything seen, instead of an infinity.
  The difference matters to a search that differences its own gradient:
  the probe lands in the unavailable region, the difference comes back
  non-finite, the direction is meaningless and the line search stops.
  Measured on a panel with a ridge over a filter’s own parameters, the
  outer BFGS died at its 73rd evaluation with the criterion still
  improving; it now runs past it to 132, carrying the hyperparameters
  from (1, 1, 1) to (2.07, 1.44, 1.30) and the criterion from -1459.13
  to -1423.0.

- The trace names a term instead of repeating its whole specification. A
  key is the call that produced the term, and one line of an outer trace
  carried a deparsed
  `gas(p = 1, q = 1, time = t, by = ~ridge(~id), links = list(...))`
  three times over.
  [`short_keys()`](https://statmodels7.github.io/statmodels7/reference/short_keys.md)
  keeps the leading call’s first argument – so `s(x, k = 20)` and
  `s(z, k = 8)` stay apart – and everything after `::`, which is what
  distinguishes one entry of a term from another; where shortening would
  make two labels the same, none is shortened. The columns of
  `history$outer` keep the full key, where it has to stay unique and
  reconstructible.

- A fit records the optimizer that searched over its hyperparameters, in
  `methods$search`, the caller’s where one was given and otherwise the
  one chosen from what the criterion can supply. A fit says what fitted
  it rather than leaving a reader to reconstruct the default.

- `statmod(start =)` takes a STRATEGY as well as a list of values:
  [`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
  (the default, now sayable by name), `start_zeros()`,
  [`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md)
  and
  [`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md).
  A strategy answers the generic
  [`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)
  and is asked ONCE, before the alternation between the coefficients and
  the hyperparameters begins. That is the whole reason it is an argument
  of its own rather than an optimizer: handed to `inner_optimizer` as
  `chain(sa(), iwls())`, a search would rerun at every hyperparameter
  the outer criterion tried – 46 times on an ordinary fit, the folds
  times the path inside
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) –
  returning the same answer each time. `start = NULL` and a list of
  values are unchanged to the bit, which a test asserts.

- [`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
  runs a global search on the LIKELIHOOD ALONE, penalties off: what a
  starting value has to get right is the basin, which the fit will not
  correct by itself, while the penalties enter when their
  hyperparameters are estimated and at the probe values represent
  nobody’s choice. It searches only where the problem is not convex – a
  structural term’s own parameters, the blocks of terms that recompute
  their own design (`nl`, `seg` and the break-point terms), and each
  equation’s intercept – because a smooth or a random effect is a convex
  block the scoring step reaches from anywhere, and searching over a
  thousand random-effect coefficients would spend the budget where it
  buys nothing. Kinked penalties are never searched: their
  hyperparameter has a known upper end and is swept by a warm-started
  path, which a random jump would both fail to improve on and destroy. A
  structural term’s parameters do not live in the coefficient vector, so
  the search sets them into the design’s structural state and leaves its
  best there.

  ⚠️ Measured, on every model tried it changes NOTHING, and that is the
  honest report rather than an omission: `nl()` reaches
  [`nls()`](https://rdrr.io/r/stats/nls.html)’s answer from both a poor
  start and a good one, `jseg()` reaches the same likelihood from three
  starting positions, and a filter’s fit is unchanged. The toolkit’s
  structured initializations – the intercept-only MLE, `term_start()`,
  the break-point grid – already arrive. What the tests pin is therefore
  that it is SAFE (a convex fit is unchanged to 1e-10) and that it
  reaches the coordinates it claims to; it is an escape hatch for a
  model where the default fails, not a routine improvement.

## statmodels7 0.30.0

- [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  takes a `criterion`. A scoring step is not an optimizer and carries
  its own loop, but the rule that ends it is the caller’s to choose: any
  `optimizers7` criterion drives the loop, `NULL` keeps the built-in
  rule. The state the rule reads carries the score PER OBSERVATION as
  its gradient, so `crit_grad(t)` and `tol = t` are the same rule and a
  threshold means the same at n = 10 and at n = 1e7 – asserted by a test
  running both and comparing the iteration count and the coefficients.
  Passing `tol` beside a `criterion` is an error rather than a silent
  choice between them, and a rule needing a stationarity measure – which
  a scoring step does not compute – is rejected at construction, as
  [`optimizers7::check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.html)
  rejects it for a derivative-free method.

- [`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
  estimates the smallest eigenvalue from LAPACK’s condition estimator
  read on the Cholesky factor it needs anyway, instead of computing the
  whole spectrum: `rcond` is 1/(\|\|A\|\|\_1 \|\|A^-1\|\|\_1), so
  `rcond * ||A||_1` is 1/\|\|A^-1\|\|\_1, which for a symmetric matrix
  lies between lambda_min/sqrt(p) and lambda_min. The estimate therefore
  errs on the SMALL side – the verdict is conservative and can never
  accept a matrix the exact test would refuse – and by a bounded factor:
  measured over dimensions 5 to 300 and condition numbers 1e2 to 1e15
  the ratio to the true smallest eigenvalue stays between 0.29 and 0.75,
  where the two cases the test must keep apart are separated by some
  fifty orders of magnitude. Both sides are unchanged and still
  asserted: a smoothing parameter at 6e15 is accepted, two identical
  columns refused. Measured at p = 1022: **4.5x** on the operation (1.40
  s to 0.31 s) and a fit carrying a random effect over 1000 groups from
  4.00 s to 2.78 s. `src/lapack_rcond.cpp` is the wrapper and
  `src/Makevars` links R’s own LAPACK and BLAS.

- The penalty’s square-root factor is taken from the diagonal where the
  penalty is diagonal, instead of from a dense eigendecomposition. The
  eigenvalues of a diagonal matrix are its diagonal, so the two routes
  agree by construction and a test pins them together; what makes it
  worth having is how often the case arises, a ridge, a random effect
  and the Demmler-Reinsch penalty of `s()` all being diagonal, and that
  the factor is recomputed at every iteration of the scoring loop.
  Measured on a random intercept over 1000 groups at n = 20000: one
  factorization went from 0.63 s to under a millisecond and the fit from
  6.52 s to 1.97 s, **3.3x**, with the objective and every coefficient
  identical. With a smooth beside the random effect, 3.0x. A
  non-diagonal penalty – a tensor’s, an anisotropic sum – still goes
  through the eigen route, and a test asserts that it does.

## statmodels7 0.29.0

- [`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
  inverts through the eigendecomposition it already computes for its
  test, and the test takes an optional reference scale: a smoothing
  parameter a criterion legitimately sends to 6.4e15 separates the
  eigenvalues (min 29.7, max 6.4e15, the matrix strictly positive
  definite) without flattening a direction, and against `max(ev)` alone
  that read as singularity –
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and every `edf`
  refused a CORRECT fit on a poisson smooth over weak signal.
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)’s
  whole-model smoother pass the unpenalized information’s diagonal as
  the scale; a genuinely flat direction (two identical columns) is still
  refused, being small against that scale too, and a test asserts both
  sides. The eigen-inverse also removes the Cholesky that could fail
  outright at such conditioning, each eigenvalue inverted exactly and
  the shrunk directions simply reporting variances near zero.

- The inner stopping rule keeps its absolute form through the loop and
  adds a DIMENSIONLESS reading to the final verdict: max_j \|g_j\|/(n
  s_p) with s_p = sqrt(median H_jj / n) per equation. A location
  equation’s score carries the units 1/y, so on a response scaled by
  1e-3 the absolute rule’s floor sat at 1.37e-6 against the 1e-6
  threshold and a run stalled AT the optimum reported failure; the
  dimensionless reading relabels exactly that run and never stops one,
  so every trajectory is unchanged. Two stronger designs were tried and
  REFUSED by measurement: a per-coordinate normalization by (H+S)\_jj
  let penalized coordinates converge loosely at extreme shrinkage and
  moved every outer trajectory, and driving the loop with the
  per-equation form made the tolerance unreachable at y\*1e4 – the stall
  guard on the objective, whose magnitude grows with log y, fires before
  a rule 1/s_p stricter can, the inner reported failure across the whole
  corridor of smoothing parameters between the plateau and the optimum,
  and the outer search, reading those points as unavailable, never
  crossed it (1482 evaluations, cor 0.82 against the truth where the
  shipped design reaches 0.998 at every scale from 1e-3 to 1e4).

- An UNHELD structural level starts at the equation’s own data-based
  intercept rather than at zero – the mirror of “the intercept wins”:
  held, the intercept carries the response’s scale; unheld, which is the
  volatility spelling `sigma ~ gas(...) - 1`, the level itself must
  absorb it. Started at zero on a response of another scale, the filter
  reads a score of y^(2/sigma)2 at its first step and leaves the
  representable range before the search’s guard can step back: measured,
  `sigma = NaN` at y scaled by 100, where the intercept spelling of the
  same model fits. With the start from
  [`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md),
  the same volatility model now fits identically at y, 100 y and 10^4 y
  (persistence 0.813 at every scale). Only a fresh start is touched; a
  specification a fit has been through keeps the parameters it arrived
  at.

- An outer criterion that is unavailable at the STARTING hyperparameters
  is reported with its cause instead of dying inside the optimizer.
  `evaluate()`’s own comment already stated the intent – a failure at
  the start is the caller’s and is raised – but implemented it only for
  the branch that raises; the branch where the inner fit CONVERGES to a
  point whose penalized information has no Cholesky factor (a
  degenerated parameter, measured on a gamma model whose nl rate ran to
  exp(31) with `converged = TRUE`) fell through to a non-finite value
  and the optimizer stopped with “the objective is not finite at the
  starting value”, naming the point and not the reason. The error now
  says what is unavailable and why, and points at
  `outer_criterion = NULL` for inspection.

- The subformulas of a term’s own parameters (modelterms7 0.27.0) fit
  end to end, and the layer needed almost nothing: the penalty
  enumeration of 0.12.0 and the marginal route of 0.24.0 already reach a
  sub-term’s hyperparameter under the key `term::parameter::subterm`.
  Measured, not presumed: before modelterms7 retired the shorthand,
  `gas(omega ~ random(~1 | id), by = id)` was pinned against
  `deviations = "omega"` with a ridge EXACTLY – the REML hyperparameter
  to the printed digit (0.3668324), the coordinates one for one, the
  log-likelihood to 4 decimals – and `nl(a ~ ridge(~g))` and
  `seg(x, psi ~ id)` fit with their hyperparameters estimated and the
  truth recovered. All of it is pinned in `test-submodels.R`, the panel
  case against a per-group recursion written by hand.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) died on duplicate
  row names for ANY fit whose penalty sits on a structural term’s own
  parameters – the shorthand `gas(deviations =, penalty =)` of 0.24.0
  included, exposed only now because nothing called
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) on such a fit.
  [`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
  wrote a structural unit’s `cols`, which index the term’s parameter
  vector, into the per-column flags of an equation design that has no
  such columns; R grew the vector past the design and recycled the
  labels into duplicate rows. A structural unit is skipped there now,
  and both routes are asserted.

- The level’s confounding question generalizes from a name to a
  subspace. With `omega ~ Z gamma` the constant coordinates of the
  development are held exactly as the scalar level was
  (`term_level_param()` names them, and the intercept then carries the
  stationary level, with its standard error). For the rest, an
  UNPENALIZED coordinate whose column the equation’s design already
  spans is flagged with a warning through the new
  [`modelterms7::term_level_design()`](https://statmodels7.github.io/modelterms7/reference/term_level_design.html):
  holding it would change the model where the confounding is not exact,
  while a penalized coordinate is identified by its penalty, exactly as
  a deviation is.

- A structural term’s starting values come from the term itself,
  [`modelterms7::term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.html):
  zero on the unconstrained scale is no longer “the model without the
  term” for a loading on the log chart, whose natural point is a loading
  of ONE.

## statmodels7 0.28.0

- The compiled coordinate descent reads a SPARSE block, and the last
  densification in the chain is gone.

  [`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
  materialized the penalized block dense to hand the kernel an
  `arma::mat`. A coordinate descent reads one column at a time, so a
  compressed-column matrix is the storage the method wants rather than
  one it tolerates: column `j` of a `dgCMatrix` is a contiguous run of
  its own nonzeros, which is exactly the walk every step of this
  algorithm makes.

  The algorithm is written ONCE against a column accessor and
  instantiated twice, so the dense and sparse kernels cannot drift
  apart. They agree BIT FOR BIT rather than to a tolerance, and that is
  a property of the arithmetic and not luck: skipping a structural zero
  omits an addition of zero, which is exact. The tests assert
  [`identical()`](https://rdrr.io/r/base/identical.html) on the
  coefficients, the sweep count and the screening gradient, over both
  gradient routes (residual and covariance) and both densities.

  Measured at n = 20000 against the same design densified by hand, with
  the coefficients exactly identical:

  |                                | before | after      | vs dense   |
  |--------------------------------|--------|------------|------------|
  | lasso, p = 200, density 0.005  | 1.04 s | **0.23 s** | 14.4x      |
  | lasso, p = 1000, density 0.001 | 3.61 s | **0.25 s** | **233.7x** |

  The kinked branch is now faster than the smooth one on the same block
  (0.25 s against 2.37 s at p = 1000), which is the right way round: a
  lasso screens its coordinates where a ridge solves the whole system.

  The `dgCMatrix` is taken apart in R and its slots passed to the
  kernel, so the compiled code needs no dependency on the C API.

  ⚠️ The dense path’s residual is now accumulated column by column
  rather than row by row, a compressed-column matrix having no cheap row
  walk. The arithmetic is the same sum in a different order, so results
  move within rounding; every existing comparison is far above it (1e-10
  against the proximal route, 1e-13 against glmnet) and the suite is
  unchanged.

## statmodels7 0.27.0

- `methods` is declared in `Imports`.
  [`as_sparse()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  has called [`methods::as()`](https://rdrr.io/r/methods/as.html) since
  0.26.0 and nothing declared it, so `--as-cran` reported *“‘::’ or
  ‘:::’ import not declared from: ‘methods’”* – a WARNING, which the CI
  action treats as a failure. The local suite could not see it. Sixth
  instance of the rule that any package named with `::` goes in the
  dependencies in the same edit as the code that names it.

- A kinked penalty composed with a grouping indicator ABORTED THE
  PROCESS.
  [`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
  handed the compiled coordinate descent the penalized block sliced out
  of its equation’s design, and where the equation also carries a random
  effect that design is a `dgCMatrix`; the kernel takes an `arma::mat`,
  so the conversion threw at the `.Call` boundary and R died rather than
  raising something a `tryCatch` could see.
  `y ~ lasso(~x) + random(~1|g)` was unfittable. The block’s own columns
  are dense whatever the rest of the equation holds – the sparse columns
  belong to the other terms and are never touched – so they are
  densified at the slice. Sixth round of the sparse conversion, and the
  first where the failure was not an R error.

- [`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
  asked whether a penalty has a proximal table AT A STEP OF 1, which is
  a different question from the one it meant. SCAD and MCP have no table
  past their convex region, and under a diagonal map the condition is
  `t < (a-1)/d^2`, so a standardized block on a column of spread 20 was
  refused at the probe and sent to the general proximal route, which
  then raised on the same condition. The probe is taken at a step short
  enough to answer only the question about the family; the step that
  will be used is asked below it, once the working weights are known.

- Standardization needs nothing here, which was the point of putting it
  on the penalty: a standardized term fits what a hand-standardized
  design fits (coefficients `s*beta` to 8.9e-16, fitted values to
  3.6e-15), multiplying a column by a thousand leaves the standardized
  fit exactly where it was (3.6e-15) where it moves an unstandardized
  one, and a random effect’s block stays a `dgCMatrix` through the
  terms, the assembled design (density 0.064), the information and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) with a standardized
  lasso in the same equation.

## statmodels7 0.26.0

- A model carrying a grouping indicator is fitted SPARSE end to end, and
  the gain is an order of magnitude:

  | n    | groups | before  | after  | speedup   |
  |------|--------|---------|--------|-----------|
  | 2000 | 50     | 1.38 s  | 1.29 s | 1.1x      |
  | 4000 | 100    | 3.45 s  | 0.70 s | 4.9x      |
  | 6000 | 200    | 15.13 s | 0.86 s | 17.6x     |
  | 8000 | 400    | 92.81 s | 3.31 s | **28.0x** |

  The shape changes and not only the constant: 3.45, 15.13, 92.81 was
  accelerating, and 0.70, 0.86, 3.31 is nearly flat. The log-likelihood
  is identical to the last printed digit at every size and the fixed
  effects agree with `lmer` to 2.0e-05.

  It took two steps and the first alone was not enough. With the block
  sparse the fit was 0.9x to 1.4x: `crossprod` left the profile
  entirely, and `.Fortran` ROSE from 57.9 to 72.6 per cent, because
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  still ran a dense QR. **A profile share is not a speedup: removing the
  cheaper half leaves the dearer one.**
  [`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
  now keeps sparsity, where it assembled dense zero blocks per pair per
  iteration, and
  [`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
  reaches
  [`Matrix::qr()`](https://rdrr.io/pkg/Matrix/man/qr-methods.html): 695x
  at 100 groups and 75475x at 1000 on the augmented design alone.

  A sparse QR is a QR, so the conditioning property the augmented route
  exists for is kept exactly, and the choice against a Cholesky of the
  normal equations – which times the same and squares the conditioning –
  does not arise.

  ⚠️ `Matrix`’s generics must be IMPORTED: a package that only depends
  on it gets base [`t()`](https://rdrr.io/r/base/t.html), which reaches
  [`t.default()`](https://rdrr.io/r/base/t.html) and reports that its
  argument is not a matrix, from inside the assembly.
  [`bind_blocks()`](https://statmodels7.github.io/statmodels7/reference/bind_blocks.md),
  [`zero_information()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md),
  [`as_dense()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  and
  [`as_sparse()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  are where the two kinds meet; the full-inverse sites densify
  deliberately, the inverse of a sparse matrix being dense.

- `inner_method` is `inner_optimizer` and `outer_method` is
  `outer_criterion`. The first takes optimizers, like `outer_optimizer`;
  the second takes criteria.
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  stays a statmodels7 object rather than an optimizers7 one for the
  reason coordinate descent does: it needs the model, not `fn` and `gr`.

  The joint step passed `optimizer = NULL` unconditionally, so
  `inner_optimizer` was ACCEPTED AND IGNORED for exactly the models it
  matters most for. Measured now that it is honoured: `lbfgs` is 2x on a
  three-parameter filter (5.08 s against 9.95) and `newton` is 2.8x on a
  33-parameter panel (4.58 s against 12.97, 14 iterations against 130),
  the log-likelihood identical throughout. The default stays Newton,
  which wins where the parameters are many.

- The hyperparameters are ESTIMATED by default. `outer_criterion` is
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  and it applies to the SMOOTH penalties and comes into play if and only
  if the model carries one – a property of the model and not of how the
  argument was written, so typing the default changes nothing. A model
  with no penalty, or one whose only penalty is kinked, leaves it
  unused. On a smooth the difference is real: the effective degrees of
  freedom go from 9.01, which is no penalization at all, to 7.81.

  `hyper` still means HELD AT THESE VALUES and steps the default aside.

- `sparse_criterion`,
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  by default, chooses the hyperparameter of a kinked penalty along a
  path over its own values. Where a model carries both kinds the path is
  outside and the marginal criterion is estimated inside each of its
  points, so a smoothing parameter comes from REML and a lasso’s lambda
  from BIC in one fit – which one argument could not express.

  A prediction-error criterion for the smooth penalties nested inside
  that path is REJECTED: it scores the same quantity at two levels, and
  measured, every point of the path came back NA, so the path had
  nothing to choose between and the hyperparameter kept its starting
  value while the fit reported success.

- The top of a path now empties the block. The documentation described a
  check – “a top whose fit is not empty is doubled until it is” – that
  was not implemented, so the grid covered a nearly flat stretch of the
  criterion and the choice fell on its edge: on eight coefficients of
  which three carried real signal the top read 26.5 where the block
  first empties near

  500. `min_ratio` is 1e-4, which is glmnet’s own ratio now that the top
       is comparable to its `lambda.max`.

  ⚠️ The warning that reports a choice at an end fired on the INDEX
  alone while its text claimed the criterion was still falling. With the
  top now emptying the block the criterion is FLAT across that stretch,
  so index one is a legitimate minimum and the message was naming a
  cause that was not the real one. It compares against the neighbouring
  point now, and the two warnings the suite carried were both of that
  kind.

## statmodels7 0.25.0

- A structural term is reported under the names its literature uses, and
  under the quantities rather than the coordinates. For a score-driven
  term that is `omega`, `alpha1` and `beta1`: the last is the
  AUTOREGRESSIVE COEFFICIENT, taken through the Levinson-Durbin
  recursion, where the free coordinate is a partial autocorrelation and
  above `q = 1` a different number.

  The standard error is the delta method over the JOINT variance,
  `sqrt(J V J')` with `J` from
  [`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
  not one entry of a diagonal: a coefficient reads the whole chart, so
  the covariance between its coordinates enters. The interval is built
  on the scale that keeps each quantity in its own set and mapped back,
  as every interval in the toolkit is, and a quantity that reads a held
  parameter is reported without one rather than with the variance of the
  rest.

  Nothing else moved: for a term whose coordinates are already its
  quantities the base method reports them on the parameter scale with
  the diagonal Jacobian of their links, which is what this function did
  before.

## statmodels7 0.24.0

- [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  and
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  reach a penalty over a structural term’s own parameters, where they
  used to return a criterion that did not move.

  A marginal criterion integrates the quantities its penalty shrinks.
  Where that penalty is a term’s own – the deviations of a panel – those
  quantities are not coefficients, so a determinant taken over the
  coefficients alone does not depend on the hyperparameter at all and
  the criterion is the penalized likelihood, whose maximum in a
  shrinkage parameter is at no shrinkage.
  [`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
  spans the coefficients AND the term’s free parameters, which is the
  order
  [`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
  already carries them in, with the penalty’s own Hessian in the tail.
  Nothing was derived: both pieces existed.

  Measured on a panel of three groups, a ridge on every deviation, the
  mode refitted at each value: the criterion moves from -157.361 at a
  prior scale of 0.01 through a maximum of -157.347 at 0.05 to -189.023
  at 50, and the deviations follow, sd 0.0008 to 0.148. Against a
  Laplace formula assembled by `numDeriv` over the same unknowns it
  agrees to 1e-5.

  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  integrates the penalized coordinates alone, through the penalty’s own
  range basis, so a population value the penalty does not cover is
  profiled rather than integrated, exactly as an unpenalized coefficient
  is.

  The EXACT gradient rejects here and the search stays derivative-free.
  That route reads how the determinant moves with the mode through a
  contraction of the family’s third derivative against the design
  blocks, an assembly that assumes the predictor is `X beta`; a filter’s
  level is a recursion of the term’s parameters and the same contraction
  is not the derivative there.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) and the summary’s
  structural table report a term whose own parameters carry a penalty.

  Both inverted the unpenalized information, on the comment that a
  term’s parameters carry no penalty. That is not a conservative choice:
  the deviations of a panel are identified by their penalty and by
  nothing else, a constant added to a population value and taken off
  every deviation leaving the filter exactly unchanged, so the matrix is
  singular along that direction. The table reported a missing standard
  error for every parameter of the term.
  [`structural_penalty_block()`](https://statmodels7.github.io/statmodels7/reference/structural_penalty_block.md)
  places the penalty in the tail, written once because the fit, the
  variance matrix and the criterion all need it and a block each of them
  placed for itself would agree only by accident.

- A panel with deviations fits again.
  [`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
  reaches
  [`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html),
  which rejected deviations, so every such model – with or without a
  penalty – signalled an error rather than fitting. It needs
  `modelterms7` 0.21.0.

## statmodels7 0.23.0

- A structural term of the filter shape is fitted in the SAME system as
  the coefficients rather than alternated with them.

  The alternation was never a statement about the model: the exact
  gradient of both blocks and the exact observed information over both
  together were already available, the second as
  [`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md),
  which was built for [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  discarded for the fit. What it cost is filter runs. Each sweep handed
  the term’s parameters to an optimizer of their own – and always
  `lbfgs()`, whatever the caller asked for – whose every iteration
  re-ran the recursion and its adjoint with the coefficients held at a
  point that was about to move.

  [`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
  runs one Newton step over the stacked coefficients and the term’s free
  parameters. The unknowns are ordered as the information orders them,
  so no permutation is needed, and a level an intercept in the same
  equation carries is held and leaves the system exactly as it leaves
  the information. Measured on a panel of 25 groups and 750
  observations: 16.91 s to 3.21 s, a factor of 5.3, at a slightly higher
  maximum (-1248.7824347 against -1248.7825410). A term of the
  likelihood shape – `regime()` – keeps the alternation, its information
  being assembled by a different route.

- A penalty over a structural term’s own parameters is read from the
  term’s own vector.

  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  looked its positions up in the design, which a structural term does
  not have: `blocks[[term]]` was NULL and the positions then indexed the
  equation’s coefficients, 25 of them where the equation has one. The
  penalty was evaluated at `NA`, and everything built on it followed –
  an objective that was not finite, an outer criterion reported as
  unavailable, and an inner run that spent its whole budget.

  Such a unit is now marked as structural and evaluated from the term’s
  state, and
  [`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
  returns its derivative and its Hessian in the term’s own parameters.
  The gradient of the structural block adds that derivative, which its
  objective had been including all along: two functions differing by a
  penalty are not each other’s gradient, and optimizers7’s own check had
  been saying so. Verified against `numDeriv` at three points, agreeing
  to the reference’s own accuracy.

- The per-observation callbacks a filter runs reuse one buffer instead
  of rebuilding a list of scalars on every call. Worth 1.09x on its own;
  the number of filter runs, which the joint step addresses, is where
  the time actually was.

- [`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
  decides positive definiteness from the eigenvalues rather than from
  whether [`chol()`](https://rdrr.io/r/base/chol.html) raised.

## statmodels7 0.22.0

- A term’s effective degrees of freedom are its share of the WHOLE
  model’s smoother matrix, which is the definition, rather than a
  block-wise approximation to it.

  [`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
  computed `tr[(H_bb + S_b)^-1 H_bb]` from each term’s own block,
  dropping every coupling between that block and the rest of the model.
  The definition is `F = (H + S)^-1 H` over the stacked coefficients of
  every equation, with a term taking the trace of its diagonal block;
  the matrix is now formed once per call and every term reads its share
  off it.

  It was found verifying the counts on parameters other than the mean.
  Measured on a gaussian with a smooth in each equation, the two totals
  were 16.98939 and 16.98885, and the whole of the gap sat on the MU
  smooth while sigma’s agreed exactly: the Demmler-Reinsch basis is
  orthogonalized against the constant in the UNWEIGHTED metric, and the
  mean’s information carries weights `1/sigma^2` that vary as soon as
  the scale is modelled, so the orthogonality the construction arranged
  does not survive the weighting. The gap is small wherever the blocks
  are nearly orthogonal and is not bounded in general.

  Now exact against the full-model trace, per term and in total, on a
  gaussian with smooths in both equations, on a gumbel, whose location
  and scale are NOT information-orthogonal, and on a gamma with its
  dispersion modelled. A kinked penalty keeps its own count, the number
  of coefficients away from the kink, since the curvature the smoother
  matrix is built from does not exist at a coefficient sitting on one.

## statmodels7 0.21.0

- `summary(fit, correct = TRUE)` carries what estimating the
  hyperparameters cost into the degrees of freedom.

  The ordinary effective degrees of freedom read the smoothing
  parameters as though they were known, and they were chosen from the
  same data, so every criterion built on the count is too generous.
  [`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
  propagates their uncertainty into the coefficients by the implicit
  function theorem: with `J = -(H+S)^-1 d2rho/dbeta dtheta` from
  penalties7’s `penalty_cross()` and `V_theta` the inverse of the outer
  criterion’s own Hessian, the corrected count is
  `tr[(V_b + J V_theta J') H]`.

  It is the first of the two terms mgcv sums into `edf2`, and it
  reproduces mgcv’s to 3.0e-04, 4.7e-04 and 2.6e-05 at n = 200, 400 and
  2000 on a univariate smooth. mgcv’s second term corrects for the
  Gaussian scale, which it profiles out of the fit; here every
  distribution parameter carries its own equation and its own
  coefficients, so that uncertainty is already inside the information.

  Nothing in it is smooth-specific: it enumerates penalties, so a
  `random()` effect is corrected by the same code path with no branch.
  Where no hyperparameter was estimated by a marginal criterion – a
  kinked penalty, whose lambda
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
  skips by construction, or a held one – the correction is exactly zero
  and says so. That is a refusal rather than an approximation: the map
  from the hyperparameter to the penalized mode turns a corner whenever
  a coefficient joins or leaves the active set, and a delta method needs
  a derivative that does not exist there.

  It is off by default because it changes a number a reader may be
  comparing with an earlier fit.

## statmodels7 0.20.0

- A structural term’s own parameters are COUNTED and REPORTED.

  A term that rewrites the predictor contributes no design columns, and
  the degrees of freedom were read off the columns, so `gas(1, 1)` was
  counted as ZERO and a model carrying four parameters reported two.
  Every criterion built on the total was that much too generous. What
  such a term spends is its own parameters, one apiece, less any level
  an intercept in the same equation already carries – which is held
  rather than estimated and is not paid for twice.

  [`summary()`](https://rdrr.io/r/base/summary.html) now prints them
  with standard errors and intervals from the joint observed
  information, which spans the coefficients AND the term’s parameters
  and was already being computed for
  [`vcov()`](https://rdrr.io/r/stats/vcov.html). Each interval is built
  on the unconstrained scale and mapped back, so a persistence stays
  inside `(-1, 1)` and a gap stays positive; the standard error on the
  parameter scale is the delta method and is reported for reading. A
  held level is shown as held, with no standard error rather than a zero
  one.

- [`logLik()`](https://rdrr.io/r/stats/logLik.html) says which
  likelihood it is, and the other one is reachable.

  What it returns is the CONDITIONAL log-likelihood – the density at the
  fitted coefficients, a penalized coefficient among them – paired with
  the effective degrees of freedom. A criterion built on that pair is
  the conditional AIC, and
  [`summary()`](https://rdrr.io/r/base/summary.html) now labels it as
  such rather than printing a bare `AIC`. A mixed-model package reports
  the MARGINAL likelihood instead, integrating its random effects away
  and counting variance components, and its AIC is a different number
  answering a different question (Vaida and Blanchard, 2005). The two
  are not comparable and mixing their halves is what neither convention
  allows.

  `logLik(fit, type = "marginal")` returns the value
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  or
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  optimized while choosing the hyperparameters, with the estimated
  parameters as its degrees of freedom. On a random intercept at 3000
  observations it is -4372.79 on 4 degrees of freedom against lmer’s
  -4371.71 on 4. Where no marginal criterion ran it is an error rather
  than a number that would look like one.

## statmodels7 0.19.0

- A structural term’s callbacks go through
  [`distributions7::distrib_kernel()`](https://statmodels7.github.io/distributions7/reference/distrib_kernel.html)
  rather than through the derivative generics. A filter calls them once
  per time step, and each call was validating its arguments, aligning
  theta by name, dispatching, and assembling every component of its
  order to keep one of them.

  Measured on a gas(1,1) scale equation at 2000 observations, against
  rugarch on the same data: 123 seconds before this and the two changes
  below it, 17.7 after, with the log-likelihood identical to the digit
  at every size tried. `regime()` shares the callbacks and its standard
  errors are unchanged to the digit.

## statmodels7 0.18.0

- The callbacks a structural filter runs do their fixed work once
  instead of once per time step.
  [`structural_callbacks()`](https://statmodels7.github.io/statmodels7/reference/structural_callbacks.md)
  recycled every distribution parameter to the sample length inside the
  closure and rebuilt the second-derivative component key with
  [`match()`](https://rdrr.io/r/base/match.html) and
  [`paste()`](https://rdrr.io/r/base/paste.html), both on every step of
  every filter pass. Worth about 12 per cent of a gas fit.

  It is not a change of order:
  [`rep_len()`](https://rdrr.io/r/base/rep.html) on a vector already of
  the right length returns it untouched, so the old form was linear in
  the number of observations like this one. What remains is the two
  generic calls per step, which is where the rest of the time is.

## statmodels7 0.17.0

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) on a model carrying a
  `regime()` term reports the OBSERVED information. What it inverted
  before was the complete-data information, the ordinary one averaged
  over the smoothed states, which is the matrix a scoring step inverts
  and is not the information the fit carries: the two differ by the
  conditional variance of the complete-data score, so the complete-data
  one is the larger and every standard error read off it is too small.
  Measured on a two-regime gaussian at 300 observations, the matrices
  differ by 25 per cent and the worst standard error is understated by
  30 per cent where the regimes overlap; the gap closes as they separate
  and the states become known, which is what the missing-information
  principle says it should do.

- [`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md)
  supplies the model’s side of
  [`modelterms7::term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.html):
  how each equation’s unknowns reach each predictor, and the family’s
  first and second derivatives at the predictor each regime shifts to.
  Unlike a filter’s, those cannot be looked up from one evaluation, so
  the family is evaluated once per regime, vectorized over observations.

  As for a filter, the matrix spans the coefficients AND the term’s own
  parameters, and the coefficient block of the joint inverse is what is
  reported; a level an intercept already carries is dropped rather than
  estimated. Checked against `numDeriv` on the exact gradient away from
  the optimum, where a wrong Hessian and a right one differ.

## statmodels7 0.16.0

- [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  fits `regime()`, which is the last term modelterms7 defines. Every one
  of them is now fittable.

  A term of that shape does not report a predictor: its contribution is
  a likelihood mixed over latent states, so
  [`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
  takes the term’s own and not the density at any single point.
  Everything else follows from Fisher’s identity – the derivative of
  that likelihood in ANY predictor is the posterior-weighted derivative
  of the ordinary one – so the score differentiates the log-density once
  per state, vectorized, and weights by
  [`modelterms7::term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.html).
  Against `numDeriv`: 1.8e-07 in the coefficients on a scale of 246, and
  6.0e-09 in the term’s own parameters. On 600 observations simulated
  from the model the fit returns levels (-0.033, 3.057) against a truth
  of (0, 3) with the scale at 0.973.

  The matrix the scoring step inverts there is the COMPLETE-DATA
  information, the ordinary one averaged over the smoothed states, which
  is what an EM step inverts. It is positive definite and it is not the
  observed information of the mixture, which needs Louis’s identity;
  with an exact gradient a scoring matrix has only to be positive
  definite.

- Where a structural term’s level and a linear intercept are both
  present, THE LINEAR INTERCEPT WINS and the term’s level is held at
  zero. The two are exactly confounded – shifting the intercept by `c`
  and the level by `-c(1 - sum b)` leaves every predictor unchanged –
  and a fit reached that ridge reporting convergence, the score being
  small because the surface is flat. `y ~ x + gas(...)` and
  `y ~ x + regime(...)` are ordinary things to write again, and nothing
  is lost: what a constant cannot express is the dynamics, or the
  difference between one regime and another.

  Which parameter is the level is the term’s own answer, through
  [`modelterms7::term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.html);
  which one is dropped is this layer’s, since only it knows what else
  the equation carries. The question is asked of the SPAN of the
  equation’s design and not of a column named `"(Intercept)"`, a factor
  coded without one summing to the constant just as well. The held
  parameter leaves the information as well as the fit, its row and
  column dropped from the joint matrix.

  Measured, and this is what says the two parametrizations are one
  model: with the level held and with it free the fits reach the same
  maximum (-553.793723 against -553.7937276) and the fitted intercept,
  0.952936, is the other fit’s stationary level `omega/(1 - b)`,
  0.95295.

- `fit@structural` carries `held`, the parameters an intercept is
  carrying instead.

## statmodels7 0.15.0

- The observed information of a model carrying a filter is exact.
  [`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
  assembles

  ``` R
  -d2l/du2 = -sum_t w_t sum_{q,r} l_qr,t V_q,t' V_r,t
             - sum_t w_t l_p,t E_t
  ```

  over the coefficients of every equation AND the term’s own parameters,
  with `E`, the second derivative of the predictor through the
  recursion, from
  [`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html).
  Only the filter’s equation has a `V` that is not its own design; the
  third derivatives the second sum needs are distributions7’s, in closed
  form for every family.

  What it replaces is the naive `X'WX`, which the scoring step still
  uses and is right to – with an exact gradient any positive definite
  matrix converges to the right point – but which is not the
  information. Measured against numDeriv on the exact gradient, AWAY
  from the optimum, where a wrong and a right Hessian are told apart:
  4.3e-11 relative. The naive curvature differs from the exact one by
  17.5 per cent of the coefficient block on a one-equation model and
  35.3 per cent on a two-equation one, the second larger because the
  cross-equation feedback contributes there too.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) inverts the joint matrix
  and takes the coefficient block, rather than inverting the coefficient
  block alone. The two differ wherever the coefficients and the filter’s
  parameters are correlated, which is the situation the filter creates;
  the joint route is what accounts for having estimated the term’s
  parameters rather than known them.

- A filter’s level and a constant in its own equation are rejected. They
  are EXACTLY confounded: shifting an intercept by `c` and the level by
  `-c(1 - sum b)` leaves every predictor unchanged, the recursion being
  affine in the level given the score path and the score path depending
  on the predictor alone. The fit reached that ridge reporting
  convergence and only the information saw it, so
  `reject_confounded_level()` asks at design time, where the terms can
  be named.

  It asks whether the CONSTANT LIES IN THE SPAN of the equation’s design
  and not whether a column is called `"(Intercept)"`: a factor coded
  without one still sums to the constant, and the span test rejects it
  where a name test would not.

- No compiled code for any of it, and the measurement is the reason. One
  curvature costs 1.4 filter runs and is FLAT in the number of unknowns
  – five against fourteen changes nothing – so the matrix arithmetic is
  not where the time is; the family’s third derivatives cost 0.0000 s.
  What it is spent on is the per-observation callbacks, and there the
  win was in R: a curvature runs at a point already reached, so the
  score and the curvature of the density can be computed vectorized once
  and looked up, where the filter must evaluate them at a predictor it
  has just produced. That change took the curvature from 10.83 s to 3.95
  s at n = 8000, a factor of 2.7, and the answer is bit-identical.

## statmodels7 0.14.0

- [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  fits the terms whose block depends on their own coefficients: `seg()`,
  `jump()`, `jseg()` and `nl()`. They were rejected before, because
  assembling such a block once and solving for the coefficients
  estimates something else and reports convergence.

  The block is refreshed at the coefficients inside the objective, and
  the predictor takes the term’s CONTRIBUTION rather than the block
  times the coefficients. Writing the difference as a per-observation
  adjustment leaves every crossprod reading the block as it did, and
  makes the scoring step’s increment the Gauss-Newton one. Measured
  against modelterms7’s own iteration, which shares no code with the
  fitting layer: the break-point agrees to 0.01 and is the same from
  every start tried, and `nl()` agrees with
  [`nls()`](https://rdrr.io/r/stats/nls.html) to its own tolerance.

  The refresh is committed once per sweep and not once per objective
  evaluation, since the rescaling factor of a discontinuous term is a
  state of the iteration rather than a function of the point: a schedule
  advancing at the speed of a line search is not the one the
  construction was designed with, and one that never advanced would
  solve a permanently smoothed problem.

  The verdict on such a fit is the objective’s and the term’s own,
  through
  [`modelterms7::term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.html),
  not the inner score’s. Where the block is a working linearization with
  a frozen weight the score belongs to that working model: measured on
  `jump()`, the fit reaches the break-point and the jump size to three
  figures and the score stays at 0.176 forever.

- [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  fits a structural term of the FILTER shape, which is `gas()`. The term
  contributes no columns; its level is added to the predictor of the
  equation it sits in, and its own parameters are estimated on the
  unconstrained scale of `term_links()` by a general optimizer with an
  exact gradient, alternating with the coefficients.

  What makes the answer right is the gradient. The level was driven by
  scores read at predictors the coefficients also enter, so the
  derivative of the objective in a coefficient is not the block times
  the score;
  [`modelterms7::term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.html)
  supplies the correction and the layer chains it with its own mixed
  second derivatives, which reaches the coefficients of the OTHER
  equations too. Against `numDeriv` the whole gradient agrees to 2e-8,
  and the naive one is out by 9.7. On 2000 observations simulated from
  the model the fit returns (0.299, 0.440, 0.732) against a truth of
  (0.3, 0.4, 0.7), with the scale at 0.998 against 1.

  `fit@structural` carries each term’s estimated parameters, on their
  own scale and on the unconstrained one.

- `regime()` is still rejected, and the message now says what the layer
  is missing – the derivative of a likelihood mixed over latent states
  in the predictor it is handed – rather than naming an internal
  generic.

- A formula carrying more than one structural term is rejected, whatever
  equations they sit in. Each rewrites the predictor or the likelihood
  the others would be read at, so the pair is one recursion written as
  two and neither term implements it.

- An accepted step that moves the objective by less than its own
  rounding ends the scoring run. The sufficient-decrease test cannot
  tell such a step from no step – at a step of 1e-9 the decrease it asks
  for is of that order – so the run spent its budget standing still,
  which is what a term whose gradient belongs to a working model
  produces at its own fixed point.

## statmodels7 0.13.0

- A term carrying more than one penalty is fitted, counted and reported
  per penalty.
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  already enumerated them, so the fitting core was right; what still
  read one penalty per term were the three places that report a fit.

  - [`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
    hands
    [`modelterms7::edf()`](https://statmodels7.github.io/modelterms7/reference/edf.html)
    the hyperparameters of every penalty the term carries, keyed by the
    penalty names, which is the shape that function reads them in. It
    asked for `term_penalty()` before, so a term whose penalty covers
    part of its parameters was counted as unpenalized – its coefficient
    count – and a term carrying two had the hyperparameters of neither.
  - [`term_block_kind()`](https://statmodels7.github.io/statmodels7/reference/term_block_kind.md)
    reads the same enumeration, so a term penalized over part of itself
    is a penalized block, and a selection where any of its penalties has
    a kink, rather than a parametric one.
  - [`summary()`](https://rdrr.io/r/base/summary.html) shows one
    hyperparameter row per penalty. Where a term carries several the
    name carries the penalty as well (`delta.lambda`), two
    hyperparameters in one block not being the same number.

  The row of the degrees of freedom stays per TERM, which is the
  granularity a table of terms wants.

- [`statmod_entry_key()`](https://statmodels7.github.io/statmodels7/reference/statmod_entry_key.md)
  composes the key – the term’s name, with the entry’s after `::` where
  a term carries several – in one place. Two callers composing it apart
  would agree by accident.

## statmodels7 0.12.0

- [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  enumerates the penalized units of a model – one entry per penalty,
  whatever term it belongs to and whether or not that term has more than
  one – and the fitting core reads it instead of looping over terms
  itself. Twelve places ran the same loop, over the distribution
  parameters and over each one’s terms, and every one of them assumed a
  term carries at most one penalty. Three of the twelve are the fitting
  core and now read the enumeration: the block split, the penalty
  assembled at given coefficients, and the starting hyperparameters.

  The key is the term’s name in the formula, with the entry’s own name
  after `::` where a term carries several. A term with one penalty over
  the whole of itself keys exactly as it did, so nothing that reads a
  hyperparameter by term name changes, and two `ridge()` terms remain
  two terms with two hyperparameters – which they already were, at every
  layer.

  What this opens is a penalty over parameters that are not coefficients
  of a design block, which is what `gas()` on panel data, `nl()` and
  `seg()` need for a population value with shrunk deviations per group.

- The sweep is finished: the hyperparameter index, the integrated basis
  of
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  the exact gradient and Hessian of the marginal criterion, and the
  labels a summary reads all enumerate penalties rather than terms.
  [`statmod_unit()`](https://statmodels7.github.io/statmodels7/reference/statmod_unit.md)
  looks one up by parameter and key, which is the same answer where a
  term carries one penalty and the right answer where it carries
  several.

  What still reports per term rather than per penalty is `edf()`, which
  is the right granularity for a table of terms and would need
  [`modelterms7::edf()`](https://statmodels7.github.io/modelterms7/reference/edf.html)
  to accept more than one set of hyperparameters. A term with two
  penalties currently has its degrees of freedom reported as its
  coefficient count.

## statmodels7 0.11.0

- The path carries the size of the kink at the point just fitted onto
  the blocks, so the next point can screen against it. It travels there
  rather than through the argument list of every layer between the path
  and the descent: it is a property of the block, where its penalty was
  a moment ago, and the path rebuilds the blocks at each point anyway.

  **Measured, the rule still discards nothing**, and the reason is the
  grid rather than the data. On a geometric grid of ratio r the rule
  keeps a coordinate whose gradient exceeds `(2 - 1/r)` times the kink.
  At the default 25 points over `min_ratio` 1e-3 the ratio is 0.75 and
  the threshold two thirds of the kink, which almost every inactive
  coordinate clears; at `glmnet`’s 100 points over 1e-2 the ratio is
  0.955 and the threshold 0.95 of the kink, just under it, so most
  inactive coordinates fall below. The rule lives on fine grids and ours
  is deliberately coarse, each point being a whole fit. Measured over a
  path: 5 per cent of the columns discarded and 0.93 to 0.96 times the
  speed.

- A fit that did not converge now reports the range each distribution
  parameter reached. The case it exists for: a lasso at a fixed
  hyperparameter with a free scale, on a design the model can
  interpolate. Fitting the coefficients shrinks the residuals, which
  shrinks the scale, which raises the working weights, which makes the
  penalty count for relatively less, which lets more coefficients in. At
  200 observations and 400 columns the scale reached 3.8e-15 and 380 of
  the 400 coefficients survived, where the same block fitted at a held
  scale kept the five that were real.

  Nothing diagnoses that. The parameters it reached are a fact and a
  scale at 1e-15 says the rest, where naming a cause would mean choosing
  a threshold for running away – and the same fit at 100 columns
  converges to a scale of 0.77 that is nothing of the kind. What that
  fit is doing is the joint mode over coefficients and scale, which is
  degenerate wherever the likelihood can be driven up without bound; it
  is the reason the hyperparameters here are estimated by a marginal
  criterion or by prediction rather than jointly.

## statmodels7 0.10.1

- The coordinate descent kernel carries the two devices `glmnet` uses to
  go faster, and measured at a single hyperparameter **neither pays**.
  Both are path devices and the reason each one loses is visible:

  Covariance updating replaces an O(n) read of the gradient with an O(m)
  one, at the cost of building a column of X’WX the first time a
  coordinate moves off zero. That column costs O(nm), so it is worth it
  only while m is small: at 5000 observations with 200 columns and
  nothing screened away, taking that route cost 70 milliseconds against
  55 for the residual.
  [`coord_covariance()`](https://statmodels7.github.io/statmodels7/reference/coord_covariance.md)
  sets the threshold where the two cross rather than at a rule of thumb.

  The sequential strong rule discards a coordinate whose gradient at the
  previous point of a path is below `2*s_k - s_{k-1}`. With no previous
  point the reference is the kink that empties the block, and
  `2*s - max|g|` is negative at any value worth fitting at, so the rule
  discards nothing and costs one crossprod. It is therefore not
  attempted without a previous point, and passing one from the path is
  what is still owed.

  The rule assumes the gradient moves no faster than the threshold,
  which is not a theorem: it can discard a coordinate that belongs in
  the fit. What makes the answer exact is the check afterwards, reading
  the gradient over every column at the point reached and putting back
  whatever exceeds the kink. A test screens the block down to one column
  and to two and gets the same answer back; without the check, one
  column gives a different one.

## statmodels7 0.10.0

- A compiled coordinate descent fits a block whose penalty is separable,
  in `src/coord_descent.cpp`, the package’s first compiled code. It
  reads the block’s own columns and the running residual – the model
  rather than the objective – which is why it is not an optimizer and
  lives here. Measured against the proximal route on the same block,
  agreeing with it to 1e-10 on the coefficients and to 1e-9 on the
  objective:

  |             | n = 200, p = 20 | n = 1000, p = 100 | n = 5000, p = 200 |
  |-------------|-----------------|-------------------|-------------------|
  | lasso       | 26x             | 135x              | 27x               |
  | elastic net | 40x             | 9x                | 43x               |
  | SCAD        | 21x             | 154x              | 24x               |
  | MCP         | 353x            | 156x              | 7x                |

  The proximal route read the whole model at every step: one block fit
  at 200 observations made 88 evaluations of the objective, 75 of the
  gradient and 83 of the operator, each over every parameter of the
  distribution, and closed in 36 iterations where a coordinate descent
  closes in six sweeps.

- The penalty arrives as a piecewise linear table from
  [`penalties7::penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html),
  so the kernel names no family and a penalty that describes its
  operator gets the compiled route without an edit here. A penalty with
  no table – a heavy-tailed prior, whose operator is a root – keeps the
  proximal route.

- **Against the reference packages we are still slower**, and the lasso
  is the only comparison whose objectives match exactly, being
  homogeneous of degree one in lambda so that a rescaling of the loss is
  absorbed. There we agree with `glmnet` to 0, 1.4e-14 and 1.1e-13 at
  the three sizes and take 4.2, 10.6 and 12.6 times as long. For the
  elastic net, SCAD and MCP the shape parameters make the penalty
  non-homogeneous, so a lambda rescaling cannot match the two objectives
  and only the timing compares: 0.22x to 0.02x against `glmnet` and
  `ncvreg`.

  Of the remaining distance, at n = 5000 and p = 200 the kernel is 21 ms
  of a 60 ms block fit and the R around it is the other 36; what
  `glmnet` has beyond that is covariance updates and strong rules, which
  this kernel does not.

## statmodels7 0.9.0

- [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  chooses the hyperparameters of a penalty with a kink, which no
  criterion in the package could reach before. A marginal criterion
  approximates an integral at the penalized mode and asks for the
  penalty’s second derivative there; the mode of a lasso, a SCAD or an
  MCP sits at the kink for every coefficient it sets to zero, which is
  where that derivative does not exist. Cross-validation asks about
  prediction instead and needs nothing from the penalty beyond a fit.
  `rule = "1se"` takes the largest kink whose criterion is within one
  standard error of the smallest.

- [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  reach those hyperparameters too, the effective degrees of freedom now
  being the trace over the coefficients that are away from the kink –
  the submodel the fit selected, where the penalty is twice
  differentiable. For the lasso, which is linear there, that is exactly
  the number of surviving coefficients (Zou, Hastie and Tibshirani,
  2007); the elastic net keeps its quadratic part and spends less, SCAD
  and MCP their own curvature and spend slightly more. The trace over
  the whole vector could not see the selection at all: measured on
  twenty noise columns it read 14 at every value of lambda, so a
  criterion built on it would have priced a model that selects nothing
  the same as one that selects everything.

- Neither is a gradient search. The penalized mode is piecewise smooth
  in the hyperparameter, differentiable while the active set holds and
  turning a corner whenever a coefficient joins it or leaves, so the
  value is swept over a grid. The grid is geometric in the size of the
  kink – read from the penalty by probing its subdifferential rather
  than assumed to be lambda – from the value that empties the block down
  to `min_ratio` of it, with every fit warm-started from the previous
  one. A choice at either end of the grid is reported: it is the grid’s
  and not the criterion’s.

  On 200 observations of twenty columns with three true predictors,
  every route keeps all three.
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  keeps four in all, `cv(rule = "1se")` six – the same count as
  `cv.glmnet`’s own one-standard-error rule – and `cv(rule = "min")`
  sixteen against glmnet’s nineteen.

- [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  at new data reapply each term’s recorded block instead of rebuilding
  it. A term records how its block was made, and
  [`modelterms7::term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.html)
  replays that record; rebuilding gives a block of the same shape
  multiplying the same coefficients that means something else. Measured
  on `y ~ s(x, k = 10)` at 200 observations, predicting on 40 of the
  rows the model was fitted to differed from the fitted values there by
  0.237, and on the 51 rows with `|x| < 0.5`, where rebuilt knots move
  furthest, by 1.19. The whole data handed back agreed exactly, which is
  why nothing noticed, and cross-validation would have inherited it.

- The information matrix, its approximation and the budget are read off
  `inner_method` in one place, so a refit inside a path or a fold runs
  on the terms the caller asked for.

## statmodels7 0.8.1

- `scad()` and `mcp()` are fitted by the proximal scheme, as `lasso()`
  already was. The hyperparameters were passed to the penalty as a named
  numeric vector where the contract asks for a list, which
  `penalty_kinks()` accepts for some branches and rejects for these two;
  the failure was caught and read as “no kink”, so both terms went into
  the jointly fitted system – solved by the curvature of a function that
  has none. On a design of twenty noise columns the fit kept 19.00
  effective degrees of freedom out of 20, which is no selection at all,
  and `y ~ scad(x)` stopped with “\$ operator is invalid for atomic
  vectors” before reaching it.

  A penalty that stops when asked whether it has a kink is now reported,
  naming the term: the two schemes differ by exactly that property, so a
  penalty that does not answer cannot be assigned to one of them.

  The point the proximal scheme reaches satisfies the subdifferential
  conditions of the objective (stationarity 8e-8 where a coefficient is
  away from zero, and the unpenalized score inside the interval the kink
  opens where it is not), and the shrinkage answers the smoothing
  parameter: 17, 9 and 0 of the twenty columns survive at lambda 1, 5
  and 20. Against
  [`ncvreg::ncvfit()`](https://pbreheny.github.io/ncvreg/reference/ncvfit.html)
  on the same objective the two land on different stationary points of a
  non-convex problem, ours the lower on both families (52.9966 against
  53.5323 for SCAD, 52.9948 against 53.6638 for MCP), and starting the
  iteration from theirs moves 6e-2 back to ours.

## statmodels7 0.8.0

- Five defects a Student t fitted to `iris` exposed, each of which had
  been silent.

  [`statmod_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_start.md)
  read `distrib_start()`’s result by parameter name, where that result
  is a list of starts each keyed by parameter, so the name matched
  nothing and every fit began at zeros on the link scale – a location of
  0 for a response centred at 5.84. The intercept of each equation now
  starts at the intercept-only maximum likelihood estimate, which for a
  gaussian is the sample mean and standard deviation exactly. That
  estimate draws its own starting values at random, so the stream is
  fixed for the length of the call and restored afterwards: without it
  the same call returned log-likelihoods of -103.49, -112.11 and -111.83
  on consecutive runs.

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  stopped with “missing value where TRUE/FALSE needed” when the gradient
  at an accepted step was not finite, the line search having checked the
  objective alone, and
  [`pd_repair()`](https://statmodels7.github.io/statmodels7/reference/pd_repair.md)
  raised from three frames further down on the same matrix. Such a point
  is now reported as unusable and the last good iterate is kept.

  The alternation set `converged` to `TRUE` as soon as there were no
  sparse blocks to alternate with, whatever the inner method had done,
  which is what let the three above go unnoticed. The verdict is the
  inner fit’s.

  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) named two causes
  for a singular information, the fit not having reached a maximum and
  two columns carrying the same information, and on this fit neither
  applied: the design is full rank and the score is 4e-5. The message
  names the directions instead – which coefficients carry the flat
  eigenvector, or which rows are not finite.

- A term the fitting scheme does not cover is rejected when the
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
  iterated through
  [`modelterms7::term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.html),
  the break-point held at 2, 3, 5 and 8 from those four starts against
  5.0681 from every one of them, at a residual sum of squares of 907.5,
  478.1, 143.1 and 857.1 against 143.1. Since the block of `seg()` is a
  Jacobian, the fit also used its break-point column as an ordinary
  regressor, giving a continuous construction a fitted mean with a step
  of -2.93.

## statmodels7 0.7.1

- The exact Hessian of a marginal or prediction-error criterion is now
  exact for a separable penalty too, penalties7 having stopped
  differencing its second-order pieces. A ridge and a random effect go
  through no difference at all.

- [`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
  takes an order and does not call the second-order generics at order 1,
  so a penalty that supplies only `penalty_dhessian()` gives an exact
  gradient rather than being rejected for a quantity the gradient does
  not use.

## statmodels7 0.7.0

- [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  choose the hyperparameters by an estimate of prediction error,
  `-2l + k*tau` with `tau = tr[(H+S)^-1 H]`, and both derivatives are
  exact. The envelope theorem does not apply to these: the
  log-likelihood alone is not stationary at the penalized mode, so its
  derivative carries `db/dt` from the first order, and what makes it
  computable is that at the mode `dl/db` equals `drho/db`.

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
