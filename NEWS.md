# statmodels7 0.94.0

* `iwls()`'s sufficient-decrease condition had the WRONG SIGN, and is now
  `armijo_ok()`, named so that the sign carries a test of its own. The
  increment solves `(H + S) delta = -g`, so `g'delta` is negative -- measured,
  in 19 of 19 solves of one fit, over `[-6435, -0.0119]` -- and Armijo's
  bound `f(x) + c1 s g'delta` therefore sits BELOW `f(x)` and asks for a
  decrease. Written with that term SUBTRACTED it changes sign and permits an
  INCREASE of `c1 s |g'delta|` instead, up to 0.6435 a step on that fit,
  where the loop's own comment said sufficient decrease rather than mere
  non-increase.

  ⚠️ **It changed no fit measured.** On the model that motivated the
  investigation, 0 of 18 accepted steps had increased the objective and the
  trajectory after the fix is identical iteration by iteration; over a
  battery of seven models nothing moved. That is also why a test over a fit
  could not have caught it, and why the condition is named: it takes the
  form it replaced as its negative control, which accepts an increase where
  this one refuses.

* `iwls()` raises the Levenberg damping where the LINE SEARCH FINDS NO
  ACCEPTABLE STEP, which is the standard trust-region loop and the exit that
  a stalling fit actually takes. The loop broke there, so the damping the
  method already carries -- added for a Student t whose `nu` reached its
  clamp -- was unreachable on that path; it was raised only from the
  `stalled` branch, which such a run never sees.

  Measured over the same battery: six of the seven models are unchanged to
  six decimals in log-likelihood, in flag and in variance matrix. The
  seventh, a count model whose dispersion runs to its Poisson limit, goes
  from 18 iterations to 28, its objective from -5239.635133 to -5239.635085
  and its held coordinates from 6 to 5, at 3.4 s against 6.1 s.

  ⚠️ **It does NOT cure that model's drift**, and the drift is not a failure
  of the step: measured, every accepted step decreases the objective and the
  score stays at 9.0e-03 whatever the damping does. What the fit is doing is
  descending a flat tail -- from `alpha = 100` to infinity the likelihood
  moves by 0.12 -- and the coordinate is not identified there, which is what
  `vcov()`'s hold now reports rather than hides.

* Broadening the `stalled` trigger from `step_used < 1e-3` to `step_used < 1`
  was measured on the same seven models and changed NOTHING, so it was not
  taken: the depth of the shrink is not what separates the deadlock from a
  working linearization at its own fixed point, and that separation is
  already carried by the full step the second one takes.

* `vcov()` HOLDS a coordinate the penalized information carries nothing
  about and reports the rest, where it refused the whole matrix before. A
  combination `c'beta` has a finite variance exactly where `c` is orthogonal
  to the null space of `K`, so a coefficient the null space does not touch
  keeps its variance whatever happens to the others; the held ones come back
  `NA`, as `stats::vcov.lm(complete = TRUE)` returns for an aliased column.
  It runs ONLY where `solve_pd()` has already refused, so a fit whose matrix
  is invertible is untouched.

  `uninformative_coords()` is the new function and `held_condition()` the
  classed warning it raises. A candidate is a coordinate with no curvature of
  its own -- a diagonal that is not finite, which is what a parameter run out
  of its range leaves and is read by `boundary_coords()`' rule rather than by
  row, a coordinate at a bound making its whole row non-finite; a diagonal
  that is not positive; or, among the rest, one carrying the whole of a null
  direction of the equilibrated matrix.

  **A candidate is held only where holding it disturbs nothing**, and the
  test has two branches because it has two scales. Where the diagonal is
  positive the coordinate has a scale of its own and the question is the
  Schur correction `K_Aj K_jj^-1 K_jA` it removes, which in equilibrated
  units is `max_k (|K_jk|/sqrt(K_kk))^2 / K_jj` and is dimensionless: a
  coordinate scaled down by 1e-14 leaves that ratio at one and is not held.
  Where the diagonal is not positive the ratio is 0/0 and says nothing, so
  the row is read against the matrix's own scale instead. Measured on a
  Poisson-inverse gaussian regression whose dispersion left its range,
  `K_jj` is exactly zero, the row's largest entry is 1.2e-15 against a
  matrix whose largest is 2562, and `||K_j.||/||K||` is 4.4e-19.

  Where the flat direction is a COMBINATION -- two collinear columns, the
  null vector `(e_1 - e_2)/sqrt(2)`, each carrying half of it -- nothing is
  held and the refusal stands: dropping one would report the variance
  conditional on the other as if it were marginal.

  Measured, holding is the Moore-Penrose inverse restricted to what is kept,
  to 2.6e-18 and 1.7e-18 on the two shapes the tests pin. On the count model
  above, `vcov()` goes from an error to 51 usable coordinates of 57, with the
  45 coefficients of the mean carrying standard errors within 2.6 per cent of
  `gamlss`'; with the dispersion held inside its range one coordinate is held
  and they are within **0.17 per cent**.

* `inner_settings()` gives an \pkg{optimizers7} optimizer the OBSERVED
  information. It minimizes `-loglik + penalty` and `optimizers7::minimize()`
  documents `he` as that function's second derivative, which is what the
  observed information is; the expected one is a different matrix, so a
  method asked for a Newton step was performing Fisher scoring under another
  name. `approx` is `"opg"` there rather than `"bartlett"`, which is what
  `iwls()` itself defaults to; it is read only where `expected` is `TRUE`, so
  for an optimizer it records a default rather than a choice.

  What it cost before is not a constant factor. Where a family writes no
  closed expected information the `"bartlett"` route is the exact sum over
  the support: measured on a Poisson-inverse gaussian regression at
  n = 12096 with 57 coefficients, ONE evaluation of `he` took **211.409 s**
  against 0.032 s for the outer product and **0.041 s** for the observed
  information. `newton()` calls `he` once an iteration, so an
  `inner_optimizer = newton()` fit of that model had not finished after
  fifty minutes and was stopped. Six families reach that route: pig1, pig2,
  pseudohuber, skewnormal1, skewnormal2 and skewt.

  Measured before and after in one process, four families by three
  optimizers, every fit landing in the same place: the log-likelihood agrees
  to five decimals throughout; the coefficients agree to 5.3e-09
  (gaussian1), 0 (poisson), 7.4e-09 (gamma1) and 2.5e-07 (pig2) under
  `newton()`, and **exactly** under `bfgs()` and `lbfgs()`, which never call
  `he` and are the control. `newton()` on the pig2 model goes from
  **183.95 s to 0.05 s**, and on a gaussian from 1.07 s to 0.38 s.

* `fit_expected()` follows it: an optimizers7 optimizer no longer counts as
  having inverted the expected information, so `vcov()` reports the matrix
  the step used, which is the rule that function's page already stated.
  Measured on the three families that write a closed expected information --
  where the report therefore moves from the expected to the observed one --
  the standard errors change by at most **3.6e-11**: at the mode the score
  equations make the two sums coincide.

# statmodels7 0.93.0

* A family with no closed-form expected information no longer costs a sum
  over the support at every scoring iteration. `iwls(approx =)` defaults to
  `"opg"`, the outer product of the observed scores, and so do the twelve
  internal entry points that carried `"bartlett"` as a hardcoded default --
  which is what a caller who passes no `approx`, `vcov()` among them, was
  reading. Measured on `statmod(y ~ x, pig1_distrib())` at n = 500: **89.06 s
  to 0.64 s, 139x**, the same log-likelihood and coefficients agreeing to
  9.6e-07. The step is unaffected in substance, the score being exact and any
  positive definite matrix reaching the same stationary point.

* THE STEP AND THE REPORT NOW TAKE DIFFERENT MATRICES, deliberately.
  `fit_expected()` answers `TRUE` only where the fit inverted the expected
  information AND the family writes it out, so `vcov()` falls back to the
  observed Hessian where the step used an approximation. Measured on the same
  regression: standard errors 5.7 per cent from those of the exact
  expectation under the outer product, 0.6 per cent under the observed
  information.

* `vcov()` takes `approx` beside `expected`, and `summary()` takes both and
  passes them on, so which standard errors are computed is the caller's to
  say. `expected = NULL`, the default, is the policy above.

* `iwls()`'s page no longer claims that asking for `approx` where it would be
  ignored is an error. No such check exists, and it cannot live in the
  constructor: an `iwls()` object is built before it meets a distribution.

# statmodels7 0.92.0

* Terms carrying the same `id` label are fitted with ONE hyperparameter
  between them, on both machines: `outer_hyper_index()` gives such a group a
  single row, so a marginal criterion searches one coordinate for it, and
  `path_rows()` gives it a single axis, so a path sweeps it once. Measured on
  two smooths whose free smoothing parameters differ by a factor of 180000,
  sharing lands on one value and costs 28.79 of REML criterion, which is what
  a restriction must do.

  The value is written back under EVERY member's own key rather than kept in
  one home, which is what leaves the twenty-six readers of a hyperparameter
  untouched: each penalty goes on being a penalty with a hyperparameter of
  its own, `hyper()` reports one row per term with the same number in each,
  and the effective degrees of freedom are still counted term by term -- 10.745
  against 10.715 on the two smooths above, where a count taken per
  hyperparameter would have reported them equal.

  The group's exact gradient is the sum of its members'. Verified against two
  references sharing no arithmetic with it: the members' own gradients read
  from an index built without the label, agreeing to 0.0e+00 from values of
  -3.433 and +2.401, and a central difference of the criterion, 6.2e-07 at
  h = 3e-3 with clean second-order convergence.

* Two labelled groups are rejected rather than approximated. A label written
  on a smooth penalty and on a kinked one names two machines and cannot have
  one value; a label tying hyperparameters whose links differ asks one free
  value to lie in two domains.

* ⚠️ A shared group has no exact outer HESSIAN, so the search falls to
  `lbfgs()`. The criterion's second derivative is assembled from per-term
  tables keyed by the pair of hyperparameter names, and a row standing for
  members of two different terms would need the sum of two such tables under
  one key. The gradient is exact, which is what the search runs on; measured
  on the two smooths, `lbfgs()`, `bfgs()` and the default reach the same
  criterion to six decimals as `nelder_mead()`.

* `hyper()` gains an `id` column, so a reader seeing the same number twice is
  told that the agreement is the model rather than a coincidence.

# statmodels7 0.91.0

* The terms carrying the same covariance label and the same grouping are
  collected into one class, whose penalty covers their stacked coefficients.
  `y ~ random(~ 1 + x | b | id) | sigma ~ random(~ 1 + x | b | id)` is a
  random intercept and slope in the mean correlated with a random intercept
  and slope in the scale, one covariance block estimated across the two
  equations. Measured on 60 groups of 15 drawn at standard deviations 0.8 and
  0.5 and a correlation of 0.75, the fit returns 0.774, 0.502 and 0.759, the
  effects correlating 0.984 and 0.944 with the ones drawn, and the criterion
  is -810.59 against the -829.41 of two independent random effects.

  What moves is the penalty and not the coefficients. Each member's block
  stays in its own equation and `coef()`, `fitted()` and `predict()` are
  untouched; one penalty now reads columns from more than one of them.
  `statmod_classes()` collects the classes, `class_index()` interleaves the
  members' positions group by group -- which is the order a blockwise penalty
  reads, penalties7 reshaping its argument by row -- and `unit_beta()` is the
  one place a unit's coefficients are gathered, answering the same vector for
  an ordinary unit as the per-equation form it replaces.

  The joint prior belongs to the class. Where no member names a `distrib` the
  default is a centered multivariate gaussian on `parameters7::dr_prod()`,
  whose coordinates are the log standard deviations and the correlations'
  angles, so a printed hyperparameter is the quantity it names; at a total of
  one column there is no correlation and the default is the centered
  univariate gaussian a single unlabelled term builds, which makes a class of
  one member the same fit as `random(~ 1 | g)` to the last bit.

  Refused, at the specification and with the reason: a label used on two
  groupings, a `distrib` named on more than one member, and one whose
  dimension is not the class's total.

* A label written inside a subformula is refused, naming the parameter the
  labelled term develops. `term_tags_deep()` reaches it through
  `modelterms7::term_components()`, a labelled term inside
  `seg(x, psi ~ ...)` not being one of the equation's terms. The block would
  have to reach coefficients this layer addresses differently, which is work
  the shared penalty does not yet do.

* What a fit reports about a covariance class. `hyper()` names every equation
  the block spans, `mu, sigma` rather than the first of them: `param` on a
  class's unit is the first member's and is a convention the hyperparameter
  store is keyed by, which a table a reader reads must not repeat as though
  it were a fact. `summary()` reports the class's hyperparameters once, under
  its first member, and adds a note naming the label, the grouping and the
  terms that share the block; without it a reader sees a covariance of four
  coordinates under a term carrying two columns and nothing saying where the
  other two came from. A class of one member shares nothing and gets no note.

  A labelled random effect is also reported AS a random effect. The kind of a
  block was read after its penalties, and a labelled term declares none, so
  it came back `parametric`: forty grouping indicators printed one per line
  under a heading saying they are an unpenalized block, with the covariance
  that produced them nowhere on the page.

  The effective degrees of freedom needed no change and that was measured,
  not assumed: `statmod_edf()` reads one smoother matrix assembled over every
  equation at once, so the class's cross block is already in it. On a shared
  block over 40 groups the total agrees with the trace computed apart to
  1.4e-14, and dropping the cross block moves it by 0.2477, which is what
  says the coupling is really being counted.

* A label written inside the subformula of an additive term joins the class
  like any other. `y ~ seg(x, psi ~ random(~ 1 | u | id)) | sigma ~ random(~
  1 | u | id)` is a random break-point correlated with a random scale: the
  subjects whose break-point comes later are also the noisier ones, and one
  block estimates it. Measured on 40 subjects of 25 drawn at standard
  deviations 0.6 and 0.5 and a correlation of 0.6, the fit returns 0.538,
  0.532 and 0.612.

  It needed no new address space, which is what made it cheap. A labelled
  effect inside `seg()` or `nl()` occupies columns of that term's block --
  exactly the ones `modelterms7::term_components()` reports as the
  component's `sub_index` -- so `label_pieces()` records them as `within` and
  the design turns them into positions in the stacked vector as it does for
  any other term. The walk recurses, composing `within` on the way down.

  A label under a STRUCTURAL term is still refused, and the message now says
  why it is the one shape that cannot be reached: what such a term
  contributes is its own parameters rather than columns of the design, and
  only `statmod_marginal_full()` spans both -- for one filter and in one
  place.

  A parent keeps its own penalties beside a labelled sub-term: the class's
  entry is added to them and never replaces them.

  ⚠️ Measured, and worth knowing before reading a flag: a random break-point
  beside a modelled scale reports `converged = FALSE` whether or not the two
  share a block -- 48.9 s labelled and 6.4 s unlabelled, both not certified,
  against the same break-point without a modelled scale, which converges in
  2.0 s, and two labelled random effects without a break-point, which
  converge in 5.7 s. The label is not what stops it.

# statmodels7 0.90.0

* `ml()` fits a model carrying an anisotropic tensor smooth, where it
  refused. The criterion projects the determinant's matrix onto the range
  of each penalty, and `te()`'s is an `additive_penalty()`, which exposed
  no null basis to take the complement of; `penalties7` 0.19.0 supplies
  one, and it is a property of the components rather than of any setting
  of the hyperparameters, so there is nothing to guess at.

  Measured on the model the refusal was written for, 200 observations of a
  surface quadratic in one margin and linear in the other: `ml()` and
  `reml()` reach smoothing parameters agreeing to four significant figures
  (3.43e+07 and 0.438 either way), the same 8.50 effective degrees of
  freedom, and fitted values correlating at 1 with an identical rmse of
  0.0569 against the truth, where an intercept-only fit is at 0.638. The
  two criteria differ in value, as two criteria do.

  The guard itself is unchanged and still reachable: `scad_penalty()` and
  `mcp_penalty()` are improper by construction and are not quadratic, so
  they have no null basis, and `ml()` refuses by name rather than guess
  which directions are profiled.

# statmodels7 0.89.0

* `statmodels7_update(quiet =)` is removed. It was accepted and read by
  nobody: the report branch printed its whole table under it, and the install
  branch hands the work to `pak::pak()`, which takes no `quiet` argument of
  its own. `statmodels7_update(quiet = TRUE)` was as loud as the default and
  nothing said the setting had been dropped, which is the shape recorded twice
  for `fit_distrib(maxit =, tol =)` and once for `ridge(n_lambda =)`: an
  argument accepted and ignored is worse than one that errors.

  Nothing is lost. The report exists to be read, and its return value is
  `statmodels7_versions()`, so the versions without the console output are one
  call away and always were. The argument is now reported by name.

  A walk over every function definition in the nine packages, comparing each
  formal against the symbols in the body and in every default, put this at
  778 candidates over 3,280 functions, of which all but 57 are S7 or S3
  methods whose formals have to match the generic; of the 57 plain functions,
  this was the only exported one. So the toolkit has no second instance.

# statmodels7 0.88.1

* Documentation pass over all 329 help pages. Every page is now
  self-contained: a `@return` names the length, the class and the names of
  what comes back instead of a bare type, every argument states its accepted
  shape and its default, and the public pages carry examples that demonstrate
  a property rather than printing an object.

  Nine claims were corrected against measurement. `xtx()`'s page said a full
  p-squared kernel ships, where the kernel accumulates the upper triangle and
  mirrors it; `statmod_commit_refresh(which =)` defaults to `"all"` and not
  to `"jacobian"`; `statmod_fitted_spec()` commits only the Jacobian entries
  and also copies a structural term's own parameters; `start_strategy()` is
  abstract and its constructor signals an error; `penalty_theta_start()`
  starts a hyperparameter bounded below at `lower + 1` and not at the
  midpoint of the unit interval; and `iwls_info_diag()`,
  `statmod_structural_table()`, `rstatmod()` and `statmod_refresh_settled()`
  each described something the code does not do.

  Two rendering defects were fixed. A `\eqn{}` in a section heading made
  roxygen emit a random placeholder as the title, changing on every run, and
  a paragraph of `iwls_solve()`'s details sat between its `@param` and its
  `@return`, so it rendered inside the last argument's description.

  No behavior changed.

* The readable variance matrix is SYMMETRIZED rather than left symmetric by
  construction. The delta method's `J V J'` collects an entry and its
  transpose in a different order, and floating-point addition is not
  associative, so the two differed in the last bit -- on Linux and not on
  Windows or macOS, which turned a bit-for-bit assertion into a red job on
  one platform of five. A variance matrix is symmetric, so the halving is the
  contract and not a repair.

* `rstatmod(n_sim =)` draws several data sets in one call, and the truth is
  drawn ONCE and shared: what a study over replicates measures is the
  variability of an estimator at a set of parameters, so the replicates
  differ in what is random and not in what is being estimated. Varying the
  truth as well is a loop over calls, and reads differently. With
  `n_sim > 1` the per-replicate fields -- `data`, `theta`, `latent` -- come
  back as lists of that length while `par` and `structural` stay single; at
  one replicate the shape is what it was.

* `rstatmod(covariates =)` draws the covariates too, one function of the
  observation count per column, afresh at every replicate. The choice is not
  a detail: with `data` the covariates are held and what is measured is the
  estimator's behavior CONDITIONAL on that design, and with `covariates` it
  is measured over the design as well -- fixed-X against random-X, two
  studies that a report should tell apart. A design that changes shape
  between replicates, a factor that lost a level being the ordinary cause,
  is reported with the replicate rather than recycled against coefficients
  drawn for another shape.

* `rstatmod()` returns the truth BESIDE the data rather than attached to it,
  which is what a simulation study reads: a list of `data`, `par`, `theta`,
  `latent`, `structural` and `call`, with a `print` that says where the two
  halves are instead of dumping a parameter at every row. They were
  attributes of the data frame, and that was worse than it looks -- an
  attribute survives a row subset without being subset itself, so
  `sim[1:10, ]` kept a `theta` of the original length silently, while
  `subset()` and `merge()` dropped it altogether. Every parameter comes back
  at one value per observation, whether its equation varies it or not, so a
  study can bind them to the rows. Passing the whole result where a data
  frame is wanted is reported with the field to use.

* `rstatmod()` simulates from the model that was written, whatever terms it
  carries. Its predictor now comes from the assembly a fit reads --
  `statmod_design_at()` and the adjustment beside it -- so a term whose block
  moves with its coefficients contributes `term_value()` at them rather than
  its block times them; measured, the two differ by 3 on a `seg()`, 4.05 on
  an `nl()` and by a missing value on a `jseg()`, so the earlier version drew
  from a linearization about the term's starting position. A term with state
  goes through `modelterms7::term_simulate()`, so `gas()` and `regime()`,
  which used to fail naming an internal generic, are simulated by their own
  recursions; what each drew is returned in the `"latent"` attribute.

* `rstatmod()` takes `n` in place of `data` for a model with no covariates,
  `structural` for a structural term's own parameters (on the scale
  `term_params()` names, and worth naming: left out they take the term's
  deliberately weak starting values), and accepts a `par` entry that is one
  number, used for every coefficient, or a FUNCTION of the coefficient count
  -- which is how a structured truth is expressed without a vocabulary for
  it, `function(k) c(2, -1.5, rep(0, k - 2))` being a sparse vector for a
  lasso to find. A response that is not a symbol is now rejected rather than
  written to a column of the wrong name, and a censored one is rejected for
  the reason `statmod()` rejects one.

* `predict()` at new data for a model carrying a score-driven term reads the
  RESPONSE to decide what was asked. Rows carrying it are a re-reading, the
  filter run over them from the term's own seed, which is what a caller means
  by predicting a model on another series and is why
  `predict(fit, newdata = <the fitting data>)` returns the fitted values;
  rows without it are a continuation and must come after the observed series.
  A frame with the response on some rows and not others is rejected.

* The ordinary generics of `stats` answer, or say why they cannot.
  `nobs()`, `formula()`, `family()` (the distributions7 object itself, so
  everything the family can do is reachable from a fit), `weights()`,
  `df.residual()` (the count subtracted is the EFFECTIVE one, so it is not
  an integer), `sigma()` (a VECTOR, and the response's standard deviation
  under the fitted law rather than whichever parameter is spelled `sigma`),
  `model.matrix(fit, what)` and `simulate()`. Three signal an error naming
  what to ask instead: `terms()`, a fit having one set per distribution
  parameter and the formula's bars not being `stats`' syntax;
  `model.frame()`, a fit keeping each term's blueprint rather than the data;
  and `anova()`, a penalized fit whose hyperparameters were chosen from the
  same data having no null distribution to test against.

* `predict()` at new data works for a model carrying a score-driven term:
  the new rows continue the series rather than being read on their own, each
  group from its own state, through
  `modelterms7::term_continue()`. Before this the ordinary assembly returned
  the fitting data's values whatever `newdata` held -- a three-row `newdata`
  came back with three hundred numbers -- so the two paths are separated
  rather than merged. A forecast reports no standard error: what
  `se = TRUE` gives is the uncertainty of the parameters, and a forecast
  carries the uncertainty of the future scores as well, which is the larger
  part and is no delta method. A term reporting a likelihood mixed over
  latent states is rejected, its reading past the data being a predictive
  distribution and not a value.

* `predict()` takes `se = TRUE`, reporting the prediction's standard error
  and an interval built on the link scale and mapped back. Every kind of
  term is covered: a design's rows for the parametric and penalized blocks,
  a Jacobian for the terms whose block moves with the coefficients, and for
  a score-driven term the level's own derivative in its parameters together
  with the correction to the equation's design rows that
  `modelterms7::term_static_deriv()` supplies. Validated against a numerical
  derivative of the predictor -- 1.7e-12 to 8.9e-11 on the ordinary shapes,
  and 1.8e-09 or better on a filter alone, beside a covariate, in the scale's
  equation, over a panel and with its level developed over covariates. A
  discontinuous break-point reports `NA` rather than a number: its block is a
  frozen working linearization, whose curvature is not the likelihood's.

* `coef()`, `vcov()` and `confint()` report the quantities the model is
  written in, and take `readable = FALSE` for the coordinates it was
  estimated on. Most coefficients are the same either way; what moves is the
  two kinds of parameter reported under a name they are not carried under. A
  discontinuous break-point term carries `g` and reports
  \eqn{\psi = -g/\delta}, so the old vector held a number that is no quantity
  of the model and the position appeared nowhere; a score-driven term's
  persistence rides a partial autocorrelation and what the literature calls
  \eqn{\beta_j} is the autoregressive coefficient the whole chart produces.
  The variance follows by the delta method and each interval is built on the
  scale that keeps its quantity in its own set.

* The three views are built from ONE map, `readable_joint()`, so they cannot
  report a quantity under one name and index it under another.

* A structural term's parameters are in `coef()` under either reading. They
  were in neither: a model whose predictor is a score-driven filter answered
  `numeric(0)` while spending nine parameters, and they were reachable only
  through `fit@structural`.

* `vcov()` keeps the structural block of the joint inverse rather than
  computing it and dropping it, so a model with a filter has a variance to
  report at all -- and it is the same matrix `statmod_structural_table()`
  reads, which used to invert it a second time for itself.

* `vcov(parameter =)` selects one equation's submatrix, and the rows of both
  readings are grouped by equation. The joint coordinate order puts every
  design block first and the structural tail last, so a filter model listed
  the scale's intercept before the mean's own parameters and did not line up
  with `coef()`.

* `hyper()` reports every hyperparameter of every penalty -- the value, what
  put it there, and whether the term held it -- on the scale the penalty
  declares it or on the free scale its link defines. They are not
  coefficients and are not in `coef()`; before this they were on
  `fit@hyper` and nowhere else.

* `residuals()` gives the QUANTILE residual by default:
  \eqn{r_i = \Phi^{-1}(F(y_i; \hat\theta_i))}, which under a correct model is
  exactly standard normal whatever the family and whichever of its parameters
  are modeled. There is one residual per observation and not one per
  distribution parameter: a residual compares an observation with the whole
  law its row was given. Where the distribution function jumps -- every
  discrete family, and a mixed one at its atom -- the construction is
  randomized, which `seed` makes reproducible without disturbing the caller's
  stream. Pearson and response residuals are offered as options; both are
  defined against the mean, and for a skewed family the first is not standard
  normal even where the model is right.

* A term whose block is a working linearization with a frozen weight is held
  out of the variance matrix, and the rest of it is reported. Measured on a
  jump at 400 observations against a bootstrap of 200 resamples: the working
  information gives the change of level and the auxiliary coordinate a
  standard error of EXACTLY ZERO, against 0.063 and 0.540, and the position
  read off them 1.8e-05 against 0.090. A zero looks like a number and is
  worse than a gap. Before this the whole matrix was refused and a sharp
  `jump()` or `jseg()` had no summary at all; now the coefficients beside the
  term keep their inference, conditional on the break-points where they are.
  The question is asked of the term through
  `modelterms7::term_jacobian_block()`, so a continuous `seg()` and a
  smoothed discontinuous one, whose blocks ARE Jacobians, keep everything.

* `print()` heads a block by its term and nothing else, counts a structural
  term's own PARAMETERS where it used to report zero coefficients, pads the
  name column to the longest name in the table, and reports no likelihood:
  there are two of them, the conditional one the criteria are built on and
  the penalized one the inner fit minimizes, and `summary()` is where the
  pairing means something.

* A term that develops one of its own parameters over covariates reports
  that parameter as a compartment of its own. Its columns mean different
  things -- a break-point's population value and its per-group deviations
  are not comparable quantities -- and a table that stacks them reads as a
  list of numbers rather than as a model. Each compartment is indented under
  the term, headed by what develops the parameter, and each sub-term inside
  it is rendered the way a block of that kind is rendered at the top level.
  The division comes from `modelterms7::term_components()`; nothing here
  parses a coefficient name.

* A random development reports the scale of its effects with its interval
  and one line saying how many predictions there are and how far they
  spread. The predictions themselves are not printed: they are a column of
  numbers nobody reads to the end, and they are in `coef()`. This is what a
  top-level `random()` term has always done, applied where the effects sit
  inside another term.

* A compartment opens with its own hyperparameter, under a name that says
  what the hyperparameter is rather than which coordinate carries it. A
  gaussian prior's `sigma` is the scale of the effects it shrinks and a
  quadratic penalty's `lambda` is a precision; under a term of a model whose
  distribution has a `sigma` of its own, the coordinate's name is the name
  of a different quantity.

* A term that develops one of its parameters is printed with that parameter
  read at a glance before its tables: the population value of the
  development and what develops it. A parameter that is a number of its own
  is one row of the table below and gets no line above it; a developed one
  is spread over a compartment where its population value is labeled by the
  development's intercept, so this is the only place the parameter's own
  name appears beside a number.

* Each equation's heading names the link it is written on. Every coefficient
  in the blocks below it is a coefficient of the LINEAR PREDICTOR, so what
  it means for the parameter depends on the link, and a summary that did not
  say which one left the reader to remember the family's defaults.

* A structural term is headed by the call and nothing else. The word
  `Structural` said what the call already shows.

* The prefix a term repeats on every row is dropped for the printing, in its
  own table and inside each compartment: the heading has already said which
  term this is. `coef()` and the summary's own tables keep the names the fit
  was built with, which are the ones another call can be indexed by.

* A block of more than twelve rows shows the first ten coefficients and says
  how many it did not show. A hyperparameter is never among them: it governs
  every coefficient under it.

* The whole of a block is formatted in one pass, its own rows and every
  compartment's together, so the columns line up down the block rather than
  restarting at each section.

* A structural term is a block like any other. It contributes no design
  columns, so its block is built from what it REPORTS -- the quantities
  `term_readable()` names, with the variance of the joint information behind
  them -- and a developed parameter becomes a compartment there as it does in
  an additive term. `gas(omega ~ random(~1 | id))` used to print its
  population value and every one of its per-group deviations in a flat table
  below the equation's blocks; it prints the scale of the deviations and one
  line saying how many there are.

* And its hyperparameter sits with the coordinates it shrinks. It used to be
  reported in a block of its own headed `Penalized ... [0 coefficients]`,
  carrying nothing but that one number, while the coefficients it governs
  were printed further down under a different heading.

* `term_block_kind()` answers `"structural"`, asked before anything else. A
  structural term carrying no penalty came back `"parametric"`, which put a
  term with no columns among the terms whose columns are one unpenalized
  block.

* A column empty on every row of a block is not printed. A quantity reported
  without a test -- a break-point's position, a hyperparameter, everything a
  structural term reports -- left the two test columns blank throughout.

* The matrix `statmod_structural_table()` inverts carries the DESIGN's own
  penalties as well as the structural term's. A random effect is identified
  by its penalty and by nothing else, so the unpenalized information is
  singular along the direction that trades its population value against its
  deviations; with a penalized block beside a structural term that direction
  is in this matrix, and the solve failed. Measured on a break-point term
  with a random development beside a score-driven level: reciprocal condition
  number 1.2e-18, the smallest eigenvalue 2.7e-13 and its direction the six
  deviations at +0.378 against the intercept at -0.378. Every standard error
  of the structural term came back missing; they are there now, and they are
  read off the same penalized information `vcov()` uses for the coefficients,
  which the two did not share before.

# statmodels7 0.86.0

* The exact Hessian of a marginal criterion is complete where a design block
  moves with its coefficients. Three quantities were missing, and they are
  added together because adding any one alone is not an improvement --
  measured on the twice-contracted fourth derivative of a curved `nl`,
  2.62e-02 today, 2.88e-02 with one of them, 2.65e-03 with another, 5.19e-08
  with both.

* Two of the three reach `tr(M dK_m/dt_l)`. One is the block's SECOND
  derivative in the coefficients, which `modelterms7::term_block_deriv2()`
  now supplies; the other reads only the FIRST and had been missing since
  these corrections were written. The third-derivative contraction carries
  the direction's own predictor, and where the block is a Jacobian that
  predictor moves with the coefficients like everything else. On a bilinear
  `f`, where the block's second derivative is exactly zero, it is the whole
  of the gap: 2.32e-02 against an arbitrary trace matrix, 2.19e-08 with it.

* The third reaches the mode's SECOND movement. `b_ml` solves a system whose
  right-hand side carries the third derivative of the penalized OBJECTIVE,
  whose second derivative where a block moves is not `K` but `K + D` -- the
  two `statmod_marginal_grad()` has kept apart through `mode_curvature()`
  since the gradient was written, and which this equation did not. Measured
  against the mode refitted at four hyperparameter values and differenced
  twice, `b_m` was already right to 3.8e-08 and `b_ml` was wrong by 7.6 to
  9.0 per cent at every step tried; with `dD/dbeta` in, the gap is 1.8e-05,
  inside the spread between the reference's own consecutive steps.

* What it is worth, against a central difference of the exact gradient with
  the mode refitted warm, at an inner tolerance of 1e-10 and read against the
  same harness's own floor of 2.3e-07 on a fixed design:

  | model | before | after |
  |---|---|---|
  | `nl`, r*x_max = 1.50 | 5.48e-04 | 1.37e-07 |
  | `nl`, r*x_max = 1.75 | 2.28e-04 | 1.57e-07 |
  | `nl`, r*x_max = 7.00 | 1.25e-04 | 1.52e-07 |
  | `nl`, r*x_max = 2.40 | 2.15e-05 | 1.63e-07 |
  | `nl`, r*x_max = 3.00 | 1.09e-05 | 1.63e-07 |
  | `seg(gamma1 ~ 0 + ridge)` | 1.46e-07 | 1.62e-07 |
  | `s(x, k = 8)`, no moving block | 2.29e-07 | 2.29e-07 |

  Every cell now reads the floor whatever the nonlinearity, where before the
  error rose as the curve straightened and again as it sharpened. The
  standard error of the hyperparameter on the nearly straight model goes
  from 2.74e-04 wrong to 6.85e-08.

* Nothing else moves. A model with no block that moves is IDENTICAL TO THE
  BIT -- coefficients, log-likelihood, criterion, effective degrees of
  freedom, outer gradient and outer Hessian -- on a smooth, a smooth with a
  random effect and two smooths, because the corrections are not small there
  but exactly zero.

# statmodels7 0.85.0

* `statmod_certificate()` takes `edge`, the free value beyond which a
  hyperparameter whose gradient has already met `tol` is reported as sitting at
  a boundary. It was written inline as 8, and it belongs beside the tolerance
  it works with, where the signature shows both rules rather than one.

* The comment on that test now says why the threshold is safe rather than only
  that it is: a coordinate is called an edge only if it has ALREADY met `tol`,
  so removing it from the interior set cannot raise the maximum that decides
  the verdict, and moving the threshold changes only whether the state reads
  `converged` or `boundary`, both of which are certified. Delete the
  `abs(g) <= tol` conjunct and that stops being true, which is the thing a
  later reader needs told. The same explanation is on the page.

# statmodels7 0.84.0

* `statmod_certificate()` says WHICH coefficient is at a boundary where the
  outer gradient cannot be read there, instead of reporting only that it is
  not finite. `boundary_coords()` names the coordinates whose curvature is
  gone -- the diagonal and not the column, a frozen coordinate making its
  whole row non-finite -- and the equation they belong to goes into
  `boundary` and into the reason beside the mode error.

  On the reference battery's `fam-studentt` that turns

  ```
  the outer gradient is not finite at the reported point
  ```

  into a sentence naming `nu`, the bound its link keeps it inside, the fact
  that the family's third and fourth derivatives are not finite there, and
  the 1.3e-10 log-likelihood units the fit sits above its own mode.

* ⚠️ **The state stays `not converged` and the reason is a measurement, not
  a convention.** That fit reaches the best criterion of the three routes
  (gap 0.000e+00) and its mode is located to ten digits, so `boundary` is
  what it deserves -- and certifying it would claim the hyperparameters had
  been verified when the gradient that verifies them cannot be computed
  there. What blocks it is upstream: the Student t's THIRD and FOURTH
  derivatives are not finite at the `nu` its own chart can produce.
  Measured, the crossover is 1e150, which is `sqrt(double.xmax)` and so the
  signature of a product of two quantities of order `nu`:

  ```
  nu          <= 1e50   1e150   1e300   double.xmax
  d3, param        0       3       4        8
  d4, param        0       7      10       13
  d3, link         0       4       7       10
  d4, link         0      10      15       15
  ```

  distributions7 0.31.0 made the score and the observed Hessian finite to
  `double.xmax` and did not reach orders three and four, which is the shape
  this toolkit records as *when a defect is a shape of mistake, grep for the
  shape*. The outer gradient reads exactly those two orders. Closing them
  takes the battery from 24 of 29 to 25.

# statmodels7 0.83.0

* The scoring step no longer deadlocks where one coordinate's curvature is
  orders below its neighbors'. This is the second of the two regimes 0.82.0
  named: there the curvature was `NaN` and the solve died for everyone, here
  it is finite but negligible, so nothing is held and the step in that
  coordinate is astronomically long. Measured on a Student t whose `nu`
  passes through 9.1e+12, where the information is `3.5/nu^4` = 5e-52:

  ```
  the step        -1790 in nu     against at most 3.9 in every mean coordinate
  diag(K)        0.2357 in nu     against 2328
  the observed    -175.3          i.e. of the WRONG SIGN
  the search      1, 0.125, 1.5e-05, 1.5e-08
  the run         stalled at a score of 2.9e-02
  ```

  A scalar step length cannot treat one coordinate differently from the
  others, which is why the line search shrinking to 1.5e-08 kept `nu`
  admissible and stopped the mean moving at all. Levenberg's \eqn{\lambda}
  can: `iwls_solve()` gains a `damp` argument adding \eqn{\lambda I}, which
  shortens a coordinate in proportion to how little curvature it has --
  `(0.2357 + lambda)/0.2357` against `(2328 + lambda)/2328` -- and it is
  added to the augmented design as further rows, so the QR route keeps its
  conditioning. The identity and not `diag(K)`: a proportional damping
  shrinks every coordinate alike and would leave the disparity where it was.

  Measured on the case it exists for, the run goes from stalling at a score
  of 2.9e-02 to converging at **1.45e-07**: the damping fires at 6e-06, the
  objective moves again, `nu` reaches its clamp and 0.82.0's hold takes it
  from there. The two repairs compose.

* ⚠️ **It is raised ONLY where the step was shrunk to nothing**, and the
  first version was not and broke four tests. `stalled` is not always
  "stopped short": a term whose block is a working linearization stalls AT
  ITS OWN FIXED POINT with the full step taken, the gradient there belonging
  to the working model rather than to the objective. Escalating on every
  stall turned a converged break-point fit into a non-converged one and
  moved three exact-gradient checks from machine precision to 1e-3, the mode
  being left worse located. The deadlock's signature is the shrunken step,
  not the stall, and that is what the condition reads.

  The damping starts at zero, is raised only there, decays by a hundred on
  every step that moves and is floored back to zero, so a run that never
  stalls performs the plain scoring iteration and is unchanged. Eleven
  control fits are identical to the printed digit and the suite's two
  warnings are the two it had.

* `iwls_scale()` reads the curvature's largest diagonal off the pieces
  without assembling it -- column sums of squares on the augmented route,
  which keeps a sparse design sparse -- so `iwls_escalate()` carries no
  constant with units: the first \eqn{\lambda} is 1e-8 of that scale and
  each further try multiplies by a hundred, eight tries spanning sixteen
  orders.

* `iwls`'s inner history gains a `damp` column.

# statmodels7 0.82.0

* A coordinate at a boundary is one coordinate and no longer stops the
  others. Where a parameter reaches the clamp its link keeps it strictly
  inside, the family's curvature there is `NaN`, and one such entry was
  enough to deny the whole system a solve and the whole criterion a Cholesky
  factor. `iwls_solve()` holds those coordinates and solves the rest, which
  is the active set the boundary defines, and `pin_boundary()` puts a unit
  diagonal in their place wherever the penalized information is factorized.

  The Student t of the reference battery is what this is measured on. Its
  `nu` reaches `double.xmax`, and there:

  ```
                            before          after
  score/n at the stop       6.2e-01        6.1e-08
  score in sigma            -617.6           1.41
  iterations                     2              7   (all at full step)
  the fit                    ERROR      converged
  ```

  the run had been stopping with `sigma` -- an ordinary coordinate --
  utterly non-stationary, because the zero step the failed solve returned
  applied to every coordinate. `nu`'s own score there is exactly 0, which is
  what says the point is a stationary point of the constrained problem
  rather than a place the fit was stopped at.

  End to end, `fam-studentt` goes from raising on all three routes to
  fitting on all three, in 0.4 to 1.1 s. On the same data with the mean
  correctly specified -- `y ~ s(x) + s(z) + random(~1|g)` -- it reaches a
  log-likelihood of -869.40, against the -869.33 of a profile computed with
  the mean HELD AT THE TRUTH, which is the independent check that the answer
  is the right one.

  ⚠️ **The dimension of the Laplace approximation drops with the held
  coordinate.** A boundary coordinate is not one the criterion integrates
  over, so it contributes neither a determinant term (`log 1 = 0`) nor a
  `2*pi`. That is why a t at its clamp and a gaussian do NOT report the same
  criterion on the same data, and their hyperparameters differ; what they
  share is the predictor, which agrees to a correlation above 0.999.

  ⚠️ **By the DIAGONAL and not by the column**, which is measured rather
  than reasoned: a boundary coordinate makes its whole ROW non-finite, cross
  terms included, so a column test marks its neighbors too. The first
  version did, held `sigma` along with `nu`, and left the fit exactly where
  it had been. The test pins the distinction.

  The shape is preserved deliberately where the matrix is pinned rather than
  reduced: some twenty consumers index into `K` and its inverse by position.

* ⚠️ **What this does NOT repair, measured and stated.** There are two
  regimes and only the first is covered. AT the clamp the curvature is
  `NaN`, the solve dies for every coordinate, and that is what is fixed.
  APPROACHING the clamp without reaching it -- `nu` = 9.1e+12, where the
  information is `3.5/nu^4` = 5e-52 -- the curvature is finite, nothing is
  held, and the step in that coordinate is astronomically long: the line
  search shrinks the WHOLE step to keep it admissible (1, 0.125, 1.5e-05,
  1.5e-08) and the run stalls with a score of 2.9e-02. That is one coordinate
  holding the others hostage through a step length rather than through a
  failed solve, and it needs the cap `optimizers7` 0.5.0 put on `newton()`
  for the same reason, or a hold read off the curvature's scale rather than
  its finiteness.

* `iwls`'s inner history gains a `held` column, the number of coordinates
  dropped from that iteration's system.

# statmodels7 0.81.0

* Whether a hyperparameter is AVAILABLE to the outer search is decided by how
  far above its mode the inner fit stopped, and not by the inner optimizer's
  convergence flag.  The two answer different questions: the flag says whether
  a stopping rule fired, while availability asks whether the criterion -- a
  Laplace expansion AT THE MODE -- is valid at that point, which is a matter
  of distance and has a natural scale in log-likelihood units.

  Measured on `y ~ s(x, k = 15) | sigma ~ s(z, k = 10)`, a gaussian with a
  smooth in each equation and nothing else: of the 38 inner fits one search
  performs, **38 are at their mode** by the second reading -- between 1e-09
  and 3e-09 against `mode_error_limit()` of 1e-03, six orders of margin --
  and **four** report convergence.  The other 34 stopped on the
  objective-stall guard with the objective already fixed to twelve
  significant digits and the score oscillating between 2.5e-06 and 3.3e-06,
  just above its absolute tolerance of 1e-06.

  Read as unavailable, those points made the outer line search backtrack
  eleven times per iteration and accept a step of 0.0026 where the Newton
  step is 1.4, so the search moved 0.005 in eta over 38 evaluations and
  stopped 4.0 criterion units short of the optimum.  That the optimum was
  really there, and that the exact gradient was correctly pointing at it, is
  established three ways: a grid of the criterion with the hyperparameters
  held peaks at eta = (0, 2) with -1558.698 against the -1562.696 reached; a
  central difference of the criterion gives (-0.360, 2.721) at h = 0.2 and
  (-0.344, 2.723) at h = 0.05 against the exact (-0.343, 2.720); and all four
  routes -- default, `lbfgs()`, `newton()` and `nelder_mead()` -- stop at the
  same place, which a simplex could not do if the neighboring points were
  merely worse.

  End to end the model now reaches **-1558.352**, 4.33 better, and its
  certificate goes from `not converged` to `converged`.  The rule ADDS points
  and never removes one -- a run whose flag says converged stays usable
  whatever the mode error reads -- so no model that fitted before can stop
  fitting: a single smooth, a smooth with a random effect, two smooths with a
  random effect, and a smooth on pure noise whose smoothing parameter runs to
  its boundary are all **identical to the printed digit** in criterion,
  certificate, flag, effective degrees of freedom and hyperparameters.

  Making the flag STRICTER is the other direction, and `piano_stabilita.txt`
  13d measured it and withdrew it: it cost a false negative on a good fit.
  `inner_mode_error()` is the new reader, sharing the evaluation context's
  own factorization rather than assembling a second one.

# statmodels7 0.80.1

* The column norms of the equilibrated sparse rank test come from
  `colSums(A^2)` rather than `colSums(A * A)`.  The second is a binary
  operation between two sparse matrices and intersects their index sets;
  the first acts on the `x` slot with the pattern untouched, and the two
  return the same numbers to the bit.  The line arrived with the rank test
  in 0.79.0 and cost 20 to 31 ms at every inner iteration: measured on the
  three shapes a penalized fit meets here, 50.00 ms against 1.74, 26.56
  against 0.76 and 19.69 against 0.64, which is 29x to 35x, and on those
  shapes the norms were 70 to 74 per cent of the whole augmented solve
  against a QR of 21.25, 9.53 and 7.50 ms.  End to end, with the
  log-likelihood, the effective degrees of freedom and every coefficient
  identical: 1.29x on a gaussian smooth with a random effect over 500
  levels, 1.18x over 1000, 1.31x on the poisson, and 1.09x and 1.07x on the
  two shapes whose cost sits in the criterion rather than in the refits.

# statmodels7 0.80.0

* The leverage diagonal over the nonzeros of two sparse designs is a kernel.
  `leverage_pairs()` built seven vectors as long as the PAIR count -- about
  190 MB at the 3.38 million pairs a REML fit over 500 random-effect levels
  reaches -- for a quantity that is a short double loop per observation.
  Measured there, it was **31.6 per cent of the fit** at 262 ms a call. The
  kernel decomposes over the rows, so row `i` is written in full by one
  thread and the answer does not depend on the count. End to end, with the
  fit identical to the bit: **1.95x** over 500 levels (12.52 s to 6.42) and
  1.35x over 1000. Against the dense route the same function takes when the
  density gate refuses the sparse one -- a computation sharing no arithmetic
  with it -- the agreement is exact.

* `statmod_pe_derivs()` forms `P Am`, `P Am PH` and `P Bm` once per
  hyperparameter rather than inside the loop over PAIRS, where each was taken
  `nh + 1` times over. It is the hoist `statmod_marginal_hess()` received in
  0.49.0, which the prediction-error route never got, that release having
  deliberately left `aic()`, `bic()` and `cv()` alone. The saving grows as
  `p^3 nh^2`, so it is invisible where p is small and not where it is not:
  measured, 1.03x at p = 213 and **1.23x at p = 512** (28.94 s to 23.50),
  with the answer identical.

# statmodels7 0.79.0

* The rank of a penalized augmented matrix is read off the JACOBI-
  EQUILIBRATED diagonal of its triangular factor, which is the correction
  `solve_pd()` took in 0.70.0 and which had not been propagated to
  `sparse_augmented_solve()`. Since `R'R = A'A`, scaling A's columns by
  their norms scales that diagonal by the same factors, so `|R_jj| / ||a_j||`
  is the diagonal of a decomposition with unit column norms: per-direction
  scaling forgives separation from any source while an exact collinearity
  stays exactly singular. Without it a matrix was called deficient for
  having columns of different SIZE -- a large smoothing parameter makes the
  penalty rows of its own block enormous beside an unpenalized one -- and
  every such solve fell through to a dense QR of the whole thing. Measured
  on `y ~ s(x) + random(~1|g)` at n = 20000 over 200 levels, **87 of 127
  solves were rejected and the dense route found full rank in all 127**;
  equilibrated, the ratio runs from 0.445 to 1.000 where the raw one reached
  7.4e-30. The fit goes from 34.1 s to 3.9 s, and against the dense route
  forced on the same data it is 20.9x on a gaussian and 14.4x on a poisson,
  agreeing on the coefficients to 1.4e-06 and on the effective degrees of
  freedom to four decimals.

* A dense penalized system takes its triangular factor from a threaded
  kernel. `augmented_solve()` reads only R, the pivot and the rank -- Q is
  never accumulated, applied or returned -- so a kernel that produces R
  alone does the whole job, and the trailing columns of each Householder
  step are independent: column k is written in full by one thread in the
  order the sequential loop writes it, so the factor is bit-identical at any
  count BY CONSTRUCTION and no arithmetic is added over the sequential form.
  It is engaged above a measured gate and only where the matrix is
  comfortably of full rank on the equilibrated diagonal, at the tolerance
  `dqrdc2` itself uses, so anything the pivoted route would call deficient
  still goes there. Measured end to end at eight threads: `s(x, k = 50)`
  2.07x, `te(x, z)` 2.82x, `ridge(X)` at p = 300 **2.89x** (164 s to 57),
  with the log-likelihood identical to the bit in two of the three.

* The count reaches `distrib_pdf()`, which no call site had ever passed it.

# statmodels7 0.78.0

* A cross-validation fold fits with the thread count `statmod(threads =)`
  was given. A fold's specification is built by `statmod_spec()`, which
  makes a fresh one, so the count did not travel with it as it does through
  `statmod_respec()`, and every fold fell back to the class default of 1:
  measured on a lasso over a gamma response at 20000 observations, `bic()`
  gained 1.85x from eight threads while `cv()` gained 1.06x, and now gains
  1.91x at an identical answer. Where the folds themselves run in worker
  processes each one still fits on a single thread, the two levels of
  `numericals7::n_threads()` not nesting; the same rule is applied to the
  runs of a product grid, which reach `fit_at()` through a specification of
  their own.

* `sm7::par_for()` takes the shape `distributions7`'s `d7::par_for()`
  settled on: the worker's loop is noinline and the sequential branch runs
  through the worker, so both branches execute one compiled copy rather
  than two the compiler may optimize apart, and the worker installs the
  calling thread's floating-point environment. Nothing in these kernels
  calls an Rmath routine, so the second buys no correctness today; what
  both buy is that the bit-identity the twins assert holds by construction
  rather than by the optimizer's leave.

* `threads` is passed to `parallelFor()` instead of being left to
  `RCPP_PARALLEL_NUM_THREADS`, so a kernel called outside a fit honours the
  count it is given rather than using every core the machine has.

# statmodels7 0.77.0

* `statmod_certificate()` tells apart the TWO ways a fit can have no outer
  gradient, which it used to lump into one `"unknown"`.

  A model with NO PENALTY -- `linpar`, `nl`, `seg`, `jump`, `jseg` -- has no
  hyperparameter for a gradient to be about, so the only question left is
  whether the inner fit reached its mode, and the mode error answers it. Those
  now carry a state.

  A model whose only hyperparameters are KINKED -- `lasso`, `scad`, `mcp` --
  has one, but it is the argmin over a path rather than the root of a
  derivative, a Laplace approximation at a mode sitting on the kink having no
  meaning. It stays `"unknown"`, and the reason now says which of the two it
  is. ⚠️ The mode error is not a reading there either: at a coefficient the
  penalty has set to zero the score does not vanish but lies in the
  subdifferential. Measured on a lasso, its 4.7e-03 is carried by a coordinate
  whose coefficient is exactly 0 and whose score is -0.715.

* Measured on the reference battery, the certificate now reads 28 cases of 29
  where it read 22, and the share passing goes from **62.1% to 79.3%**. The
  disagreement between the convergence flag and the certificate falls from 11
  cases to 6.

* ⚠️ `jump` certifies as NOT CONVERGED on that battery and the reading is
  right. Three explanations were tested and two refuted: it is not that a
  frozen working block's score fails to vanish (`jseg` sharp is frozen too and
  reads 2.0e-06), and it is not arithmetic (equilibrating K leaves the reading
  identical to every printed digit). It is misspecification: `jump()` is a
  pure step model with no linear term and that truth carries a slope and a
  slope change, so its break-point iteration never settles, its annealing runs
  to the floor, and its block reaches 3.0e+13 on the information's diagonal
  against 5.4e+05. On data where the jump IS the truth the same term reads
  4.2e-04 and cond(K) falls from 6.25e+13 to 1.14e+06.

# statmodels7 0.76.0

* `statmod_certificate()` reads three things AT THE REPORTED POINT, so what a
  fit says about itself is a property of the point rather than of the search:
  the outer criterion's gradient, how far the coefficients sit above the
  penalized mode (`g'K^-1 g / 2`, in log-likelihood units), and which
  hyperparameters have run to a boundary. Four states: `converged`,
  `boundary`, `not converged`, and `unknown` where there is nothing to read.

* Why not the convergence flag. It says whether a search stopped on its own
  rule, and measured across shapes it does not order fits by quality: on one
  model the default reported success at a criterion of -1783.47 while the same
  data under `lbfgs()` reached -1664.43 and reported failure.

* ⚠️ The state comes from the gradient and the mode error is reported BESIDE
  it, not folded into it. Measured at the reported point over six shapes, the
  gradient separates by five orders -- 4.7e-07, 7.8e-07, 5.8e-05, 7.7e-05 and
  3.0e-04 on fits that are right against 28.8 on one that is not -- while the
  mode error does not: 1.8e-16 to 6.1e-12 on four of them, 22.8 on the failing
  one, and 0.114 on a random-changepoint `seg` whose answer is right to a
  correlation of 0.9932. A number that does not separate cannot decide a
  state, and a certificate that says how far from the mode is worth more than
  a boolean that hides it. That `seg` fit certifies as converged and carries
  the 0.114 as a line of its own.

* `tol` is 1e-2 rather than the geometric middle of the two groups because the
  two ways of being wrong are not symmetric: a certificate that says NOT
  CONVERGED at a good point is visible and checkable, one that certifies a bad
  point is the failure it exists to remove.

* It costs one outer gradient and one solve, once, and refits nothing. The
  criterion reconstructed from `fit@spec` equals the one the fit reports
  EXACTLY on every shape measured, so the reading is of the fitted model.
  A form with no exact outer gradient, or a fit with no marginal criterion,
  leaves the gradient `NA` and the state `unknown` rather than differencing it
  at 2p refits.

* `summary()` carries it and prints it, boundary coordinates by name and the
  reasons in full. A quantity no view displays is a quantity nothing checks.

# statmodels7 0.75.0

* `criterion_resolution()` refuses to report a resolution where the inner fit
  is not at a mode. The reading displaces the coefficients to where the inner
  score says the mode is and asks how far the criterion moved, which is a
  resolution while that displacement is a CORRECTION and something else once
  it is not. The test is the displacement's own predicted decrease,
  `g'K^-1 g / 2`, in log-likelihood units rather than the coefficients', so
  one limit (`mode_error_limit()`, 1e-3) serves every shape.

* ⚠️ It was not the estimate that was wrong. Measured on a hierarchical
  break-point model, the inner fit reports convergence with a score of
  **247.8** -- against 3.6e-04 on a smooth and 8.2e-06 on a random intercept
  -- and given a mode located that badly the criterion genuinely is uncertain
  by 28, which is more than its whole movement over the search (18.31) and
  four orders above the 1.6e-03 its reproducibility measures directly. The
  arithmetic was right and the point was not, so what is refused is the
  conversion of an unlocated mode into a stopping tolerance, not the reading.
  **The inner fit's convergence claim at a score of 247.8 is a separate and
  larger defect**, recorded in `piano_stabilita.txt`.

* Measured, and it does what it should and nothing else. A smooth, two smooths
  with a random effect and a random intercept are unchanged under both
  optimizers -- identical evaluation counts (4, 18, 7, 27, 9, 23) and criteria
  agreeing to 4e-07, which is their own reproducibility. On the break-point
  model `lbfgs` goes from **2 evaluations, -1764.283 and `converged = TRUE`**
  to **20 evaluations and -1664.434**, a gain of **99.85 REML units**;
  `newton()` is unchanged at -1815.372, stopping for another reason.

* The trace says when a resolution was refused and what the inner fit's excess
  was, once per fit rather than once per evaluation.

* ⚠️ The guard alone was not enough, and one test said so. The rule the
  criterion carries is a NUMBER settled before the run, so it can only be
  built from the reading at the starting point -- a cold start, the
  worst-located mode of the whole fit -- while the line search reads a CLOSURE
  over the running minimum. Nesting the line search inside the test for a
  usable first reading threw that away: `seg(x, psi ~ random(~1|id))`, whose
  first reading is refused because its cold mode sits 0.046 above its own
  minimum while the RUNNING MINIMUM is 3.1e-11, went from 5 evaluations to 18
  and from converged to not, at an identical answer -- cor 0.9932 and rmse
  0.0674 either way. The closure is now given whether or not the first reading
  was usable, and that fit is back to 6 evaluations and `TRUE`.

* `outer_default_optimizer()` is the policy `exact2 -> newton()`, extracted so
  that it can be read, pinned and swept. Behavior is unchanged.

# statmodels7 0.74.0

* The outer line search gets a backtracking budget of 12 rather than
  `optimizers7`'s default of 30, on the optimizer this package CHOOSES and on
  no other: an optimizer the caller named keeps its own budget as it keeps its
  own stopping rule. 30 suits an objective costing microseconds, and here every
  trial is a penalized refit.

* It costs nothing where nothing was wrong. Measured at 30, 12 and 8 backtracks
  on a smooth, on two smooths with a random effect and on a random intercept,
  all three are unchanged in evaluations, in criterion to six decimals, in
  effective degrees of freedom and in the convergence flag -- the resolution
  `criterion_resolution()` supplies ends the search before the budget is ever
  reached. A hierarchical break-point model goes from 31 evaluations and
  25.6 s to 13 and 20.3, with the criterion 1.3e-04 better and the same
  degrees of freedom.

* ⚠️ What it buys is smaller than the evaluation count implies, and the reason
  is that every trial warm-starts from the last ACCEPTED point: as the step
  shrinks the refit begins at nearly its own answer, so removing 22 evaluations
  of 38 removed 2.8 s of 30.8, which is 0.13 s each against an average
  evaluation's 0.81.

* 12 rather than 8 because 12 takes 5.3 s of the 6.6 there are to take, and
  because swept with the optimizer named -- where the resolution does not mask
  it -- the criterion given up against a budget of 30 is 4.8e-05 at 12,
  7.8e-04 at 8 and 3.3e-03 at 6, the last being past the 1.6e-03 that
  criterion can resolve.

# statmodels7 0.73.0

* `start_search()` leaves a PENALIZED coordinate out wherever it sits. The
  rule was applied per term -- a convex block is not searched -- and a
  penalty declared through a sub-formula covers columns of a refreshable
  term's own block, which the term-level rule takes whole. On
  `y ~ jseg(x, npsi = 2, by = ~random(~1 | id))` over thirty groups that was
  210 of 219 coordinates, and with no analytic gradient the optimizer
  differences centrally, so one iteration cost 2.3 s against the
  milliseconds of a scoring step (measured: 13.06 s against 0.55 s for
  `lbfgs(maxit = 5)`, 24x). The budget was the smaller half: the search runs
  on the likelihood with the penalties off, where a deviation is identified
  by nothing at all, so it fitted each group's own points and returned a
  start further from the penalized mode than the one it began with.

* `start_search(hyper =)` is removed. It was documented and never read by
  `start_at()`, and could not have been: the search objective is the
  likelihood with the penalties off, in which a hyperparameter does not
  appear, so no proposal could change one. A global search over the
  hyperparameters is a search over the outer criterion and is available as
  `statmod(outer_optimizer = optimizers7::sa())`; the page says so.

# statmodels7 0.72.1

* The thread twins separate the two claims they were conflating: the
  kernel's count-invariance is asserted `identical()` kernel against
  kernel, and the comparison against the BLAS expression a kernel replaces
  carries a tolerance -- OpenBLAS on the CI's Linux runners and Accelerate
  on macOS block their accumulations, so `identical()` there asserted a
  property of the reference BLAS and reddened four of the five platforms.
  The `wcrossprod()`/`xtv()`/`wxsq()`/`xtx()` pages state the same
  distinction.

# statmodels7 0.72.0

* The marginal break-point terms fit end to end:
  `y ~ jump(x, psi ~ random(~1 | id), marginal = TRUE)` -- and `seg`,
  `jseg`, several break-points, an explicit prior -- route through the
  existing likelihood-shape branch (the one `regime()` rides), so the
  coefficients are fitted by Fisher's identity with the term's posterior
  component weights, the term's own parameters by `lbfgs()` on the exact
  jacobian, and `vcov()`/`summary()` invert the exact joint observed
  Hessian from `term_hessian()`, the prior scales' intervals built on the
  log scale. Everything is plain maximum likelihood: the prior is part of
  the likelihood and no outer criterion is involved.

* What the branch needed was the levels generalized, not a route:
  `statmod_regime_at()` reads the mixture shifts off the new
  `modelterms7::term_levels()`, which may answer with a per-observation
  matrix (the quadrature nodes of a marginal `seg`/`jseg` shift each
  observation by its own hinge value), and the component loops of the
  score and the information read a column where a regime reads a number.
  A structural term's fresh start receives the equation's
  `predictor_target()`, off which the marginal term reads its exact
  two-stage profile -- the multimodality in the positions does not go
  away under the marginal, and a fit started conventionally was measured
  converging on a local optimum 33 log-likelihood units under the truth.

* New `statmod_latent(fit)`: the posterior mean and standard deviation of
  each group's latent break-points, from the same decomposition the
  likelihood is computed on.

* The step kind's marginal runs on the side chain's forward recursion
  (modelterms7 0.56.0), so several break-points stopped being exponential
  in the sample: measured through `statmod()` at identical estimates,
  K = 2 goes 13.2 s to 2.5 s, K = 3 goes 721 s to 8.8 s, and K = 5 --
  beyond the old cap of three -- fits in 94 s with every population
  position within 0.05 of the truth and `vcov()` finite from the
  propagated Hessian.

* Measured (the numbers are in `piano_marginal.txt`): on the 5b design the
  marginal jump reproduces the standalone reference to the third digit
  (m 4.898, tau 0.518, cor 0.877, rmse 0.292) and over five seeds
  converges 5/5 in 2-4 s where the production smoothed probit converges
  2/5, with tau always the closest to 0.5; a marginal `seg` recovers
  exactly what the native random-changepoint fit does (cor 0.984 both,
  rmse 0.088 against 0.086) at 47 times the cost, so the native route
  remains the recommended one there; a marginal `jseg` on a Poisson panel
  recovers the positions (cor 0.71, log-likelihood above the
  truth-started point) where the mode-based routes lose them; a t prior
  on contaminated positions beats the gaussian on the inliers (rmse 0.173
  against 0.241); and a lasso fits beside the marginal term with its
  selection intact, which is the composition `piano_marginal.txt` F4
  asked to be measured rather than promised.

# statmodels7 0.71.0

* A SMOOTHED break-point term (`seg`/`jump`/`jseg` with `smoothed =` a
  `penalties7` `abs_smoother`, modelterms7 0.55.0) is routed as what it is:
  a Jacobian block. `term_jacobian_block()` answers `TRUE` for it, so it
  takes the Gauss-Newton embedding of `seg()` -- no working-fit phase, no
  holding -- and a random or penalized development of its break-points fits
  through the machinery the random-changepoint model already uses. The
  layer needed no routing edit for that, which is what the predicate of
  0.70.0 was for.

* `summary()` reports a smoothed term's smoother and width -- the width of
  the transition, the bent-cable reading -- and, for random break-points
  under the probit smoother, the corrected scale beside the apparent one:
  the convolution identity `tau_apparent^2 = tau^2 + h^2` is the
  smoother's own declaration, read through its `tau_correction`.

# statmodels7 0.70.0

* `solve_pd()` tests and inverts on the JACOBI-EQUILIBRATED matrix
  (unit diagonal), which tells a flat direction from scale separation
  whatever produced the separation, and the caller-supplied reference
  scale is gone. The reference covered one source only: a smoothing
  parameter at 1e15 was forgiven, but a break-point term's committed
  working block carries auxiliary columns near 1/(2cd) -- 1e16 on the
  information's diagonal at the annealed floor -- so the DESIGN set the
  reference and ordinary curvature at 242 in the ordinary-scale
  directions read as flat: `summary()` on a converged three-break-point
  jseg refused its own variance matrix. Per-direction scaling forgives
  both sources; an exact collinearity stays exactly singular.

* An outer criterion, a path and a fold all HOLD a frozen break-point
  block at its committed positions (`statmod_alternate(hold_refresh =)`),
  and the positions are refined once, by a full alternation at the chosen
  hyperparameters before the restarts. Running the working phase inside
  every criterion evaluation had two costs measured together: a
  three-break-point fit beside an estimated smooth took 136 s and
  reported FALSE at the right answer -- the break-point moving between
  evaluations makes the criterion path-dependent, and the phase's cycling
  flags read as unavailable points to the search -- where the held route
  takes 1.1 s and reports TRUE. Measured across the combination battery:
  jseg beside a random intercept under REML goes 6.2 s FALSE to 2.9 s
  TRUE, beside a lasso swept by BIC 110 s to 2.9 s, with the break-points,
  the BLUPs and the selection unchanged; a per-group `by` development and
  a jump in the SCALE equation fit as they did.

* A break-point term whose block is a working linearization -- `jump()`
  and `jseg()`, read off `term_jacobian_block()` -- is fitted by its own
  construction's scheme (`fit_working()`): the smooth block is fitted
  exactly at the frozen block, the break-points are read off the fitted
  coefficients and committed, and the two alternate until the read-off
  settles or the working objective stalls, which is what `segmented`
  does. The previous embedding refreshed the read-off inside the inner
  optimizer's objective, whose line search then rejected the fixed-point
  iteration's own steps -- measured on a three-break-point jseg at
  n = 10000, the fit from the TRUE break-points ended at an rss worse
  than the mean-only fit and reported `DID NOT CONVERGE` in 3 passes,
  where the working-fit phase recovers psi = (-0.500, -0.000, 0.500)
  against a truth of (-0.5, 0, 0.5). Two bookkeeping defects went with
  it: a jseg's quadratic read-off is incremental in the committed
  position, so the pass-level commit took a hidden extra step that
  changed the objective at unchanged coefficients (measured at 0.71 per
  observation, exactly zero for `jump`), and the alternation's
  relative-change rule read that jump as its own progress measure.
  `statmod_commit_refresh()` now takes `which` and returns the committed
  coefficients, `seg()`'s Gauss-Newton embedding is unchanged.

* Bootstrap restarting (Wood 2001), the device `segmented` runs by
  default, with the observation that makes it cheap: the non-convexity of
  a break-point model lives entirely in the positions, so a proposal is a
  configuration of positions and two configurations are compared on the
  EXACT PROFILE -- least squares at fixed positions, one linear fit each
  -- rather than by refitting the model. Three proposal kinds: the
  deterministic sweep of `modelterms7::seg_polish()` first, each
  break-point descended over the profile with the others held, which
  walks straight to a feature the iteration pressed a break-point away
  from; then, alternating, the same sweep on a bootstrap resample's
  profile (the multinomial counts as weights) and from random positions
  over the confinement interval. Only a proposal the profile prefers
  earns a refit, accepted on the true objective; four consecutive dry
  proposals end the loop. Measured on three evident break-points at
  n = 10000: the default fit recovers the truth in 7.6 s where refitting
  every proposal took 945 s, a dry proposal costing half a second of
  linear fits instead of a 5 to 15 s refit. The random proposals exist
  because a resample's perturbation is of order 1/sqrt(n) and stops
  escaping a deep basin as the sample grows; the count is declared by the
  break-point terms (`n_boot`, default 10, 0 disables), and the loop runs
  once at the top level, never inside an outer criterion's search.

# statmodels7 0.69.0

* The COMBINATIONS of a kinked path's product grid run over the worker
  processes of `n_threads(workers =)`, through the same machinery as the
  cross-validation folds (`cv_fold_rows()` generalized to
  `worker_map()`). Each combination restarts its warm chain from the
  sweep's own starting coefficients and `cur` moves only after every one
  is scored, so the combinations are independent by construction and the
  result is identical at any worker count -- the same bodies in the same
  order. The points WITHIN one chain stay sequential, and the reason is
  measured rather than argued (voce 8 of piano_parallel.txt): a path
  point paid cold costs 2.2-3.2x the warm chain (lasso 400x30, lasso
  8000x200, scad 2000x100), so splitting a chain either slows the
  single-process default by that factor or makes the result depend on
  the count. The continuation worry -- that a cold SCAD point lands on a
  different local optimum -- was measured and refuted: 5 points of 128
  differ by at most 1.4e-4, convergence tolerance and not another
  optimum.

# statmodels7 0.68.0

* The exact gradient's reverse pass rides the fast route too:
  `term_adjoint()` receives the filter's fast context and the thread
  count, so at every gradient evaluation its forward pass runs without R
  callbacks where the registries cover the family and the link, and its
  reverse pass reads the curvature sequence that pass returns.

# statmodels7 0.67.0

* The three second-order structural call sites (the joint information,
  the exact outer gradient's shared parts, and the re-weighted pass of
  the chain term) hand `modelterms7::term_curvature()` the score and the
  curvature as LOOKUPS and the blocks callback's pieces as DATA
  (`structural_blocks_data()`), so an eligible subformula filter runs its
  second-order recursion compiled with the panel's groups over
  `spec@threads`. The callback stays beside the data for the cases the
  kernel declines (a developed autoregressive chart, the third order).

# statmodels7 0.66.0

* A score-driven filter runs on the C-callable fast route where the
  registries cover its family and its link (gaussian1 and gamma1, identity
  and log, for now): `structural_callbacks()` hands the filter a fast
  context beside the callbacks, and modelterms7's kernel reads the score
  and the curvature through the scalar entry points of distributions7 and
  linkfunctions7 instead of calling back into R at every step of the
  recursion. Bit-identical to the callback route by construction and by
  twin test; an uncovered pair leaves the context inert. With no R in the
  loop the panel's groups run over threads
  (`statmod(threads = n_threads(k))`), each group writing its own rows, so
  the result does not depend on the count, bit for bit.

# statmodels7 0.65.0

* A structural fit no longer recomputes the curvature recursion at a point
  it has already read. Measured on the gas panel at 60 groups, 62 of the
  154 recursions of one fit revisited a point -- the criterion, its exact
  gradient and the joint step's information all read the same mode, up to
  five times each, and the recursion is 35 per cent of the fit at 226 ms a
  call. `structural_memo()` is a depth-one exact memo on the design's
  structural state (the same environment that carries zeta), keyed with
  `identical()` on the full inputs and shared across the call sites, so a
  hit returns the previously computed object itself and is bit-identical
  by construction; coefficients, criterion and vcov are `identical()`
  across the change. It stands aside where the design carries refreshable
  terms, whose blocks advance a schedule the key cannot see. The
  recursions of the 60-group fit fall from 154 to 117; what remains of
  that fit's curvature cost is the R recursion itself, whose C++ port
  (licensed by its callbacks being lookups, with per-group contributions
  summed in group order for the thread stage) is the next piece of
  piano_parallel.txt's voce 5.

# statmodels7 0.64.0

* The folds of `cv()` run over worker processes:
  `statmod(threads = n_threads(workers = 4))`. Each fold is a complete,
  independent refit, so they go by PROCESSES with the safeguards
  `optimizers7::multistart` records (sequential under pkgload, sequential
  with a warning where a cluster cannot start or load the package), and a
  fit inside a worker takes a fresh specification and is sequential by
  construction -- the two levels of parallelism never nest. Each fold now
  draws from its own seed, taken once in the parent and applied whether
  the fold runs here or in a worker, and results are collected in fold
  order, so the answer does not depend on `workers`, bit for bit; the
  caller's own random stream is saved and restored around the folds, so
  what a session draws after a cross-validation no longer depends on
  whether a fold's fit happened to consume random numbers.

# statmodels7 0.63.0

* `statmod(threads = numericals7::n_threads())` accepts the toolkit's
  thread policy. The count is validated once, travels on the specification
  (`spec@threads`), and reaches the family's compiled per-observation
  kernels and the dense assembly products as an argument; the
  process-level RcppParallel setting is sized at the fit's entry and
  restored on exit. At the default the code takes exactly the sequential
  path, and a fold of `cv()` or a respec stays sequential by construction.
* The dense weighted cross product of the assembly, `X' diag(w) X`, runs
  through a threaded kernel (`wcrossprod()`) decomposed over the ELEMENTS
  of the output, each dot product accumulated in full by one thread in the
  sequential order: the result does not depend on the thread count, bit
  for bit, which a test asserts with `identical()` on a whole fit at 1 and
  2 threads. Only base dense blocks above a measured work threshold are
  eligible; a sparse design keeps its Matrix route, where a threaded dense
  kernel would do nothing.
* `.structural_blocks()` recycles the family's derivative components ONCE,
  outside its per-observation closure. The first version read them through
  `rep_len(x, n)[i]` -- an n-long allocation, and a sort-and-paste key
  rebuild, per observation per parameter pair, O(n^2) against a total of
  O(n), measured at 17.5% of a panel fit at 60 groups and growing with the
  groups. Same hoist `structural_callbacks()` received in 0.18.0; results
  are bit-identical (coefficients, criterion and vcov compared with
  `identical()` across the change). What remains on top of that fit's
  profile is the filter's per-step scalar callbacks (`distrib_kernel`),
  which is the C-callable item of piano_parallel.txt, not this one.

# statmodels7 0.62.0

* **The intercept-only fit is written to a PARAMETRIC intercept and to
  nothing else.** `start_intercepts()` recognized an equation's intercept by
  a coefficient name ending in `(Intercept)`, and `nl()` names the intercept
  of each of its own parameters the same way. Those live on those parameters'
  own charts rather than on the predictor's, so `y ~ 0 + nl(...)`, which puts
  one of them first, had `mean(y)` written into it: on a logistic growth curve
  whose asymptote rides a log link, `mean(y) = 23.94` started `phi` at
  `exp(23.94) = 2.5e10`, with an objective of 7.0e20 and a gradient of 1.4e21
  on data whose every scale is ordinary.

  Nothing reported it. The lasso path built at those coefficients spanned
  2.8e15 to 2.8e19 where the block empties at about 300, so all 25 of its
  points were the same empty fit, BIC reported that the criterion was still
  falling at the sparse end, and the model came back converged with every one
  of its 150 subject deviations estimated as exactly zero. `parametric_intercept()`
  asks which term is the parametric block instead of asking a name.

* **A term with parameters of its own is handed the response on the scale of
  the predictor**, `predictor_target()`, which `modelterms7::term_coef_start()`
  takes as its new `target`. It exists only where the response reads the
  parameter directly -- a mean or a location, which `params_interpretation`
  says -- and a term in a scale's equation is handed nothing rather than a
  quantity invented for it. The scale is not a detail: measured on a Poisson
  whose predictor is a logistic curve with `phi = 4`, a term handed the raw
  response estimated `phi` between 52.7 and 54.8 over five samples and one
  handed `g(y)` between 4.03 and 4.07. What does NOT matter is the other
  terms' contribution -- subtracting it moved the estimate from 39.68 to 39.87
  against a truth of 40 -- so nothing is residualized.

* Together these make a nonlinear mixed model fit as written. The model of the
  report, `y ~ 0 + nl(~ phi/(1 + exp(-(time - theta)/sigma)), phi ~ 1 +
  lasso(~id), ...)` over 50 individuals, went from a population-only fit with
  every deviation at zero to recovering them, at no cost the user has to know
  about.

# statmodels7 0.61.0

* **The outer Hessian reads a block that moves with its coefficients.**
  `nl()`, `seg()`, `jump()` and `jseg()` build their block as the Jacobian at
  the current coefficients, so `dX/dbeta` enters the criterion's second
  derivative in three places -- the matrix `dK/dt`, its trace against `M`, and
  the twice-contracted fourth derivative -- and none of them had it. The
  derivative is asked of the term through
  `modelterms7::term_block_deriv()` and never differenced here.

* **And the mode moves by the penalized likelihood's own curvature, not by
  `K`**, which is the same pair of matrices `statmod_marginal_grad()` already
  keeps apart through `mode_curvature()`. Differenced against two refits of the
  mode on a weakly identified `nl`, `b_m` read off `K` is `6.6e-01` wrong and
  `5.2e-05` read off `K + Dm`; everything downstream inherits it, which is why
  the three corrections above moved the Hessian by nothing until this was
  written with them.

* Measured against a central difference of the exact gradient with the mode
  refitted, over the dimensionless nonlinearity `r * x_max` of
  `a * exp(-r x)`:

  | `r * x_max` | before | after | standard error |
  |---|---|---|---|
  | 0.5 | 2.10e-01 | 3.54e-03 | 1.25e-01 -> 1.77e-03 |
  | 1.0 | 3.20e-02 | 9.51e-04 | 1.64e-02 -> 4.76e-04 |
  | 2.0 | 1.96e-03 | 1.42e-04 | 9.80e-04 -> 7.09e-05 |
  | 4.0 | 1.84e-05 | 1.49e-05 | |
  | 7.0 | 5.97e-05 | 1.56e-05 | |
  | 12.0 | 5.69e-04 | 3.27e-05 | |

  The error is largest where the model is NEARLY LINEAR and the parameters are
  weakly identified, which is the regime in which a standard error is read to
  decide whether a parameter is estimable at all. `seg()` goes from `1.59e-07`
  to `9.35e-08`.

* On a BILINEAR `f`, where the term's own third derivative is exactly zero and
  the one omitted piece therefore vanishes, the corrected Hessian reaches the
  reference's own floor: `1.33e-05` to `2.14e-07`, and the residual GROWS as
  the step shrinks (`2.14e-07`, `9.73e-07`, `4.88e-06` at `h` of `1e-3`,
  `1e-4`, `1e-5`), which is the mode's location and not a missing term.

* What is NOT computed is `d2X/dbeta2`, the term's own third derivative, for
  which no contract exists. It is what keeps the weakly identified cell at
  `3.5e-03` rather than at the floor.

* A model with no refreshable term is untouched bit for bit -- the criterion,
  the gradient, the Hessian, the effective degrees of freedom and the
  coefficients all `0.0e+00` against the stored twin on four shapes -- and the
  corrections there are the zero matrix and the number zero rather than small
  numbers.

* `statmod_hyper_vcov()` and `statmod_edf_correction()` read that Hessian, so
  both follow with no change of their own.

# statmodels7 0.60.0

* **A kinked penalty on a block that MOVES with its coefficients is fitted
  rather than refused.** `coord_fit()` read the block as it arrived and solved
  against the linear predictor `X beta`, and a term registering
  `term_refresh()` has neither property: its block is the Jacobian at the
  current coefficients, and what it contributes is `X beta + adj`. The sweep
  therefore updated the working residual of a different model.

* Measured on `nl(~ a * exp(-r * x), a ~ 0 + lasso(~grp))` over 300
  observations and ten groups, at a HELD lambda small enough that neither a
  lasso nor a ridge shrinks: the log-likelihood was `-339.74` against the ridge
  control's `+155.45`, with the rate at `0.22` against a truth of `0.70`, and
  the fit reported success. It is now `155.4618` against `155.4548`, with the
  rate at `0.7360` against `0.7377`.

* The block is read at the current coefficients and `adj` is subtracted from
  the working response, once per SWEEP and never per coordinate: the compiled
  descent exists because the design stands still while it walks the columns,
  and the refresh state advances where the alternation commits it. A model
  with no refreshable term is untouched, the question being asked only where
  something moves.

* The KKT conditions at the point reached, differentiated numerically so that
  the check shares no arithmetic with the proximal sweep: the gradient on the
  active coordinates is `2.6e-10` to `1.5e-13` across lambda, and where every
  coefficient is shrunk to zero they all sit inside the interval the kink
  opens.

* ⚠️ What remains, and it is not this route's: the alternation between a
  moving block and the free coordinates converges LINEARLY and can stop short.
  On this model the amplitudes and the rate are strongly coupled through
  `a * exp(-r x)`, so it takes 36 passes and stops with the rate `4e-4` from
  its optimum; at `iwls(tol = 1e-12)` the free gradient falls from `9.0e-01`
  to `8.4e-04`. The ridge control does it in one pass because it has no second
  block to alternate with.

# statmodels7 0.59.0

* **The criterion's resolution is COMPUTED at the fit instead of declared.**
  `criterion_resolution()` displaces the coefficients by the mode error the
  inner fit's own score implies, `K^-1 g`, and reads how far the criterion
  moves. It costs one assembly of the criterion at given coefficients and no
  refit.

* The quadratic form `0.5 * g' K^-1 g` alone is not enough, for a structural
  reason: the criterion carries `-0.5*log|K(beta)|`, which is NOT stationary in
  `beta`, so a mode error enters there at FIRST order. Measured, it is right to
  one per cent at an inner tolerance of `1e-4`, where the second-order term
  dominates, and undershoots by 50 to 1000 times at `1e-6` and below.

* Measured against the spread of the criterion at one hyperparameter reached
  from six different warm starts, over four shapes and five inner tolerances,
  the displaced reading tracks it across SIX ORDERS OF MAGNITUDE -- `8.6e-2`
  down to `2.2e-7` -- at a ratio between 0.05 and 0.99, and separates shapes a
  formula cannot: at one inner tolerance it reads `1.7e-4` on a random
  intercept over 500 levels against `7.9e-8` on a gaussian smooth, where the
  declared `|f| * tol` of 0.57.0 was 100 to 16000 times off.

* ⚠️ It is read off the criterion the SEARCH runs, passed in as a function
  rather than chosen inside. Reading the marginal criterion of a fit whose
  search is `aic()` answers for a quantity that search never sees, and the
  number that came back stopped two such fits short of their own optimum. A
  test pins it.

* The line search is given a TENTH of it. What is computed is the spread
  between evaluations reached from different warm starts, while a line search
  compares trials taken one after another from nearly the same state, where the
  criterion is far more reproducible; and it is the more aggressive consumer,
  ending a whole iteration rather than one comparison. Measured at the full
  number on a gaussian smooth the search stops after 3 evaluations against 29
  and gives up `3.2e-06` of criterion against a resolution of `9.7e-07`, so it
  stops just before it has to; at a tenth it reaches the same optimum as a
  search that was told nothing.

* It is RECOMPUTED at every usable point rather than once at the start, and
  reaches the line search as a closure over the running minimum. The first
  evaluation is a cold start, where the mode is least well located, so a
  reading taken there is the worst of the run and would govern every step after
  it: measured on a penalized smooth, the cold-start reading is 14 to 22 times
  the best later one. The minimum rather than the latest, because a resolution
  smaller than the truth leaves the search where it was while one larger stops
  a healthy search short.

* **The four reference shapes reach the same answer for a third of the
  evaluations.** Against the same model fitted with a fresh optimizer of the
  same class, which carries no resolution: 6 evaluations against 22 on a
  gaussian smooth, 6 against 19 on a Gamma one, 7 against 25 on three smooths
  with a random effect, and 31 against 31 on a random intercept over 500
  levels. The log-likelihood is identical in all four and the effective
  degrees of freedom agree to `1.2e-05` or better. A t-prior random effect goes
  from 43 evaluations reporting failure to 10 reporting success.

# statmodels7 0.57.0

* **The outer criterion has a RESOLUTION and the stopping rule is now told
  it.** Every evaluation refits the coefficients from the running warm start,
  so the value at one hyperparameter depends on the path taken to it, and a
  line search cannot verify a decrease smaller than that. Both halves of an
  optimizer's default rule are then unreachable: `crit_grad()` asks `1e-6` of
  a gradient that bottoms out where the decrease stops being verifiable, and
  `crit_rel_obj()` asks `1e-12` of a value of order `1e4`, which is an
  absolute `1e-8`. A run then reports failure at a point it does not leave.

* The resolution is DECLARED rather than discovered, as the criterion's own
  scale times the INNER tolerance, so the two move together and tightening
  `iwls(tol =)` buys both a better-located mode and a finer rule. Measured
  directly -- the spread of the criterion at one hyperparameter reached from
  six different warm starts, over four shapes and five inner tolerances from
  `1e-4` to `1e-8` -- it is an envelope on all twenty cells, with a margin of
  2.6 at the loosest tolerance and never below. The relation is not a clean
  power, the spread falling roughly as `tol^0.7` to `tol^0.8` and depending on
  the shape as well, so nothing tighter is claimed than the measurements
  support.

* It is ADDED to the chosen optimizer's own rule rather than replacing it, so
  a run can only stop earlier and never later, and only where this package
  CHOSE the optimizer: one given by name comes with its own rule and keeps it.

* **It is inert on a healthy fit and that is the point.** On the four
  reference shapes the log-likelihood, the effective degrees of freedom, the
  hyperparameters and the evaluation counts are identical to the previous
  rule, the gradient still firing first. Where the resolution is the binding
  constraint it is not: a t-prior random effect goes from 43 evaluations
  reporting failure to 9 reporting success, at the same answer.

* Known and not covered: a search that exhausts its line search instead of
  meeting a rule. `descent.cpp` asks the criterion there with `have_old =
  FALSE`, where a rule reading a CHANGE in the objective cannot fire by
  construction, so the backstop does not reach that path.

* ⚠️ `optimizers7` 0.4.0 gives its line searches a `resolution` for exactly
  that exit, and this package does NOT set it, which is a measurement and not
  an omission. The resolution was measured directly as the spread of the
  criterion at one hyperparameter reached from six different warm starts, and
  at ONE inner tolerance it ranges over 140 times across shapes after dividing
  by the criterion's own scale -- `6.2e-11` of it on a gaussian smooth against
  `8.6e-9` on a random intercept over 500 levels. A bound loose enough to
  cover the second is four orders too large for the first, where it fires
  while the outer gradient is still `1.8e-2` and a stationary point was within
  reach; six tests asserting stationarity said so. What the layer would need is
  the resolution of THIS fit, computable from the score the inner fit stopped
  at and the penalized information `outer_context()` already holds factorized.

# statmodels7 0.56.0

* `linpar_options(sparse = NULL)` is the default, so the IMPLICIT parametric
  block settles its own storage through `modelterms7::linpar()` instead of
  being built dense unless asked. A model carrying a factor of many levels is
  then fitted sparse without the caller knowing the argument exists: measured
  on `y ~ 0 + g + s(x)` at 20000 observations, 103.690 s against 2.370 s at
  four hundred levels and 5.160 s against 0.900 s at a hundred, with the
  log-likelihood identical.

* The storage is a storage and nothing else moves. On three shapes whose
  settled storage FLIPS under the new default -- a parametric block of
  indicators, a smooth under a factor `by`, and a penalized formula route --
  the fit against the same model built explicitly dense agrees on the
  log-likelihood exactly, on the coefficients to 1.2e-16, on the effective
  degrees of freedom to 3.3e-14, on the hyperparameters to 1.5e-13, on the
  standard errors to 7.3e-14 and on the predictions to 4.8e-16.

# statmodels7 0.55.0

* `reml("expected")` and `ml("expected")` carry an EXACT GRADIENT. The
  contraction is the one the observed route already performed -- one crossprod
  per distribution parameter against the same leverage diagonal -- with
  `distributions7::distrib_dexpected_hessian()` in place of the third
  derivative and its own key builder, the expected array being symmetric in
  its first two indices only. Measured against a central difference of the
  criterion with the mode refitted from a fixed start: 1.8e-08 on a gamma
  smooth, 4.3e-08 with a second penalized equation, and 1.2e-04 to 5.0e-04 at
  3000 observations across gamma, beta and a negative binomial with random
  effects in both equations.

* What it buys is the evaluation count of a search that had none, and it is
  concentrated where the hyperparameters are many or the family is dear.
  Against the derivative-free numbers on identical data: three smooths and a
  random effect go from 965 evaluations and 481.0 s to 121 and 18.0
  (**26.7x**), and a negative binomial with random effects in both equations,
  which had not finished after 94 minutes, fits in 88 evaluations and 34.3 s.
  On a single smoothing parameter it is a wash (0.8x to 1.3x), a simplex being
  efficient in one dimension.

* ⚠️ **The mode's movement is read off the penalized likelihood, not off the
  criterion.** `v = db/dt` solves `(H_obs + S) v = -d2rho/dbeta dt` whatever
  matrix the determinant is of, which is what this file's own derivation says;
  the code used the criterion's `K` for both, and the two coincide only on the
  observed route. It was therefore invisible until the expected route existed,
  and there it is a systematic error shrinking with n exactly as the observed
  information approaches the expected one: measured on a gamma smooth at 300,
  1000 and 3000 observations, 1.9e-03, 1.4e-03 and 1.1e-04, and FLAT in the
  inner tolerance while the mode's score fell two decades -- which is what said
  a reference was not the weak side. Corrected, the same comparison is
  1.8e-08. The observed route is bit-identical throughout.

* **The exact gradient is exact where the block MOVES with the coefficients.**
  `nl()`'s block is the Jacobian, so a penalty on a nonlinear parameter -- a
  ridge or a random effect over its groups -- used to be differentiated as
  though the design were fixed. Two pieces were missing and both are now asked
  of the TERM, through the new
  `modelterms7::term_block_contract()`: `dX/dbeta` enters `dK/dbeta`
  (`u_refresh()`), and the mode's own curvature, `l_a d2eta_a/dbeta2`, is what
  separates the true Hessian from the Gauss-Newton matrix the design gives
  (`mode_curvature()`). Against a finite difference of the criterion with the
  mode refitted from a fixed start: **4.6e-03 to 4.0e-09** on the observed
  route and **1.1e-09** on the expected one.

* ⚠️ **The derivative is asked of the term and never differenced here.** A
  break-point column is a step function in its break-point, so a difference
  quotient of it diverges as the step shrinks -- measured at h, h/4 and h/16:
  3.6e4, 1.4e5, 5.8e5, against `nl`'s 0.6038 throughout. A term that has not
  written the contraction inherits zeros, which is exactly right for a fixed
  design and leaves `seg()`, `jump()` and `jseg()` where they were.

* ⚠️ **A penalty with a kink on a block that moves is REJECTED** rather than
  fitted wrongly. `coord_fit()` reads the block as it arrives and solves
  against `X beta`; a refreshable term has neither property, its block being
  the Jacobian and its contribution `X beta + adj`. Measured on
  `nl(~ a * exp(-r * x), a ~ 0 + lasso(~grp))` at a HELD lambda of 0.01, where
  a lasso and a ridge must nearly agree because neither shrinks: the
  log-likelihood was -40.41 against the ridge's 226.88 and the rate came back
  0.885 against 0.712, with `converged = TRUE`. The same lasso on a linear
  model is exact, which is what says it is the moving block. The message names
  the term and the way out; `piano_kink_blocco_mobile.txt` is how the rejection
  is lifted.

* ⚠️ **The exact gradient reads the design block AT THE MODE.** It was handed
  the design as it arrived, which for a term registering
  `modelterms7::term_refresh()` -- `nl()`, `seg()`, `jump()`, `jseg()` -- is the
  block built at the coefficients the fit STARTED from, while
  `statmod_information_at()` refreshes internally so `K`, and therefore `M`,
  were assembled on the refreshed one. The two agree for every fixed design,
  which is why it surfaced only when a penalty was put INSIDE such a term.
  Measured on `nl(a ~ 0 + ridge(~grp))` the two blocks differ by 2.07 and
  `u = tr(M dK/dbeta)` was wrong by 110 -- including on the SIGMA rows, an
  equation carrying no refreshable term at all, because the leverage diagonal
  reads the mu block whatever row is being formed. Against a direct numerical
  `tr(M dK/dbeta)`, that row goes from -131.982189 to -21.902535 against
  -21.902535, and every row of an equation with no refreshable term is now
  exact. End to end the gradient on that model goes from 4.6e-03 to 1.4e-03.
  The same fix is applied in `statmod_marginal_hess()`.

* What remains there is the term's own SECOND derivative: with `X = X(beta)`,
  `dK/dbeta` gains everything coming from `dX/dbeta`, and the criterion's own
  second derivative asks for the third. The layer cannot difference the block
  to get them -- measured at h, h/4 and h/16, `nl`'s converges (0.6038
  throughout) and `seg`'s break-point column diverges as 1/h (3.6e4, 1.4e5,
  5.8e5), being a step function in psi -- so the term has to supply them.
  Designed in `piano_nl_derivate.txt`; the gradient is left admitted meanwhile,
  the search reaching the same hyperparameter to 0.01-0.14 per cent while
  refusing would cost 1.3x on those models and 3.6x to 7.4x beside a smooth.

* `ctx_penalized()` takes the criterion's own information rather than assuming
  the observed one, and caches the two separately. It was hard-coded, correctly,
  while the exact gradient ran on no other route. `ctx_trace_matrix()` is keyed
  the same way: the projection is of THAT matrix, so one cache entry cannot
  serve both. Nothing reaches it today, a search holding one `OuterMethod`
  throughout, but it is the twin of a defect that was unreachable in exactly
  the same way until it was not.

* `ml("expected")` holds as well as `reml("expected")`, which is where the two
  matrices are most easily confused: `ml()` projects the determinant's matrix
  onto the range basis while the MODE goes on moving in the full space.
  Measured against numDeriv, 1.0e-07 on one hyperparameter and 2.6e-08 /
  2.7e-07 on two, the same quality as `ml("observed")`; and the canonical-link
  identity under `ml()` is exact to 3.5e-14.

* The route is refused at order 2, the criterion's own second derivative
  wanting the next order of the same object, and where a penalty covers a
  STRUCTURAL term's own parameters, `statmod_marginal_full()` assembling the
  joint curvature from `term_curvature()` -- the observed one -- so that branch
  has no expected criterion for a gradient to be the derivative of.

* ⚠️ Consequence to know: a fit that reported `converged = TRUE` under the
  derivative-free search may now report FALSE, on three of six shapes measured.
  The fits are the same -- lambda to five significant figures, the correlation
  with the truth to five decimals, the criterion to the printed digit -- and
  the cause is the optimizer's ORDER: the default becomes `lbfgs()` where the
  observed route uses `newton()`, and `newton()` differencing this gradient
  converges on the same shapes. It is the flag item this project already has
  open, reached from a new direction rather than a new defect.

# statmodels7 0.54.0

* `offset()` written in the formula is an offset. It was SILENTLY DROPPED:
  `terms()` marks such a term in the `"offset"` attribute and
  `model.matrix()` excludes it from the design, so the term contributed no
  column, no offset and no message, and the model fitted was the one without
  it. On a count model over person-years -- 64800 observations, a negative
  binomial with a random effect in both equations -- that moved the intercept
  from -7.51 to -0.65, which is the difference between a log rate and a log
  count, and left the random effect absorbing the offset (sigma 1.21 against
  0.20).

  The offsets come out of each equation before the interpreter sees it. Any
  equation may carry one (`sigma ~ offset(s)`), several are summed as
  `glm()` sums them, and one given through the `offsets` argument as well is
  added to it.

* The offset survives prediction, which the `offsets` argument never did.
  The EXPRESSIONS are kept rather than the numbers, so `statmod_respec()`
  re-evaluates them against the new data; a vector supplied at fitting time
  has the length of the fitting data and cannot be reused, and prediction
  used to return the predictor of a model with no offset at all.

* An `offset()` buried inside another term's formula is REFUSED, naming the
  place it belongs. `ridge(~ z + offset(o))`, `random(~ 1 + offset(o) | g)`,
  `linpar(~ x + offset(o))` and a subformula's own `ridge(~ g + offset(o))`
  all used to fit, each dropping the offset in that term's own
  `model.matrix` -- the same defect one level down, and just as quiet. It is
  refused rather than routed up because the meaning differs by where it
  sits: in a penalized term's formula it would be a contribution to the
  equation's predictor, and in a subformula a contribution to that
  parameter's own chart.

* The starting values account for the offset. `statmod_intercepts()` fits the
  distribution to the response ALONE, so its answer is the predictor the
  model should have on average and the intercept has to carry that minus what
  the offset already contributes. Without the correction the run began
  wherever the offset happened to sit: measured on a count model whose offset
  averages 6.744, the starting mean was exp(7.64) = 2080 against a sample
  mean of 2.45, and the same fit took **4.9 s** corrected against more than
  25 minutes uncorrected, at 9 scoring iterations.

# statmodels7 0.53.0

* The penalty's Hessian is accumulated in the storage the DESIGN calls for, by
  `zero_information()`'s own rule, so it is stored the way the information it
  is added to already is. It was dense until now on the argument that eighteen
  places read it and only the two a sparse design exercises would be caught by
  the suite.

  What settled it is a number that argument did not have: `random(~1|g)` over
  1000 levels is 6.620 s against 2.867 s, **2.31x**, and over 500 1.24x, with
  the matrix identical to the last bit, `logLik` and `edf` identical, and 2.02
  MB becoming 0.010 MB. The dense shapes are unmoved (0.98x, 1.05x, 1.02x).

* The win is NOT the allocation, which is 0.229 ms of a 0.992 ms call and
  would have capped the whole thing at four per cent. It is that a sparse `S`
  stays sparse downstream: `penalty_sqrt()` returns a factor of 0.010 MB where
  it returned 2.0 MB, and `augmented_solve()` no longer converts it. Measuring
  the accumulator alone would have refused a change worth 2.31x.

* `zap_nonfinite()` replaces the seven copies of `S[!is.finite(S)] <- 0`.
  That expression is correct on a base matrix and a trap on a sparse one --
  the logical index is a dense p by p matrix, so the storage would have been
  thrown away at the first consumer. On a sparse matrix only the stored values
  can be non-finite, so the same answer costs O(nnz).

* Two sites assemble a matrix spanning the coefficients AND a structural
  term's own parameters, which is dense by construction, and now densify the
  penalty explicitly before writing it into a slice. One of them failed
  loudly in the suite -- a filter beside a random effect -- and the other,
  its twin in the marginal criterion, did not fire only because no test puts
  those two terms in one formula. Seventh round of this conversion, and the
  first where the sweep for the shape found the second site before it was
  reached.

# statmodels7 0.52.0

* The penalized information is factorized in the storage the MATRIX calls for.
  `H` is already sparse wherever the design is, and `ctx_penalized()` was
  densifying the sum only because the penalty's accumulator is a base matrix.
  Where the sum is large enough and sparse enough to be worth it
  (`worth_sparse()`, whose two thresholds are the measured crossover), it is
  kept sparse and factorized as such, and the log-determinant and the full
  inverse are read off that factor.

  Measured on the penalized information of a random intercept over 500 levels,
  p = 503 at a density of 0.014, each route timed with its own factorization:
  the factorization and its log-determinant cost 0.102 ms against 10.811 ms,
  and the full inverse 3.280 ms against 25.000 ms. End to end the fit goes
  1.25x at 500 levels and **2.01x at 1000**, the gap between the operation and
  the fit being the lesson this file records three times over -- removing the
  dearer half leaves the cheaper one.

  ⚠️ Those are the figures with `Matrix`'s factorization CACHE defeated.
  `Matrix::Cholesky` stores its result in the matrix's `factors` slot, so a
  benchmark refactorizing the same object measures a cache hit -- 0.004 ms
  rather than 0.102 -- and an earlier draft of this entry quoted 0.030 ms and
  a ratio of 414x on that basis. A fit never gets the hit, the penalized
  matrix being a new one at every point, so the end-to-end numbers were never
  affected.

* Below the crossover the sparse route LOSES, so the gate is what makes the
  change safe rather than an optimization with a tail of regressions: 0.13x on
  the inverse at p = 23, 0.33x at p = 53, and 0.01x on the fully dense
  penalized information of a single smooth. The three dense shapes of the
  battery are bit-identical to the previous release.

* `pd_factor()` is now the one place a penalized matrix is factorized, and
  `pd_logdet()` is it with the factor dropped, so the verdict on positive
  definiteness is written once for both storages. The sparse route carries its
  own condition estimate (`sparse_lmin()`, Higham's one-norm estimator against
  the factor's solves): `Matrix::rcond` costs 10.3 ms at p = 503 and 500 ms at
  p = 2003, more than the factorization it would be guarding, where the
  estimator is 0.58 ms to 0.80 ms and nearly flat. As on the dense route, the
  verdict never turns on whether a factorization raised.

* The route is a property of the matrix and not of the term that built it,
  which is checked rather than asserted. With every design built the same way
  it is worth 1.38x on `0 + g + s(x)` over 400 levels, 1.33x on
  `random(~1|g)` over 500 and 1.07x on `s(x, by = g)` over 60 -- an
  unpenalized indicator block, a random effect and a factor-`by` smooth,
  gaining together and in the order their sizes predict. It also survives the
  prior changing: with `random(~1|g, distrib = fixed(student_t1_distrib(),
  mu = 0))`, whose penalty is separable and moves with the coefficients at
  every inner iteration and whose `nu` is estimated, the same route is taken
  and is worth 1.24x against the Gaussian prior's 1.45x on the same data.

* Sharing that factorization with the criterion was MEASURED AND NOT TAKEN.
  Where the criterion's matrix is the one `ctx_penalized()` holds, the same
  matrix is factorized twice at one point -- 12.4 ms spent twice at p = 503 --
  and reading the determinant from the context removes one of them. End to end
  it is worth nothing: 1.01x at p = 503, and 0.92x to 1.04x over the four
  shapes. The criterion is evaluated at many points the gradient never reaches,
  every trial point of a line search among them, and there is nothing to share
  at those.

# statmodels7 0.51.0

* The leverage diagonal is read in whichever order gives the SMALLER
  intermediate, and the mirror block is not computed twice. `x_ai' M_ab x_bi`
  is a scalar, so the `(b, a)` block is the same number as the `(a, b)` one;
  and the dense route's intermediate is `n x p_b` written one way round and
  `n x p_a` the other, for the same answer.

  Both mattered on a mixed model. Fitting `y ~ x + random(~1|g)` over 500
  levels, the scale's equation carries ONE column and the mean's 502, so the
  `(sigma, mu)` block was materializing a dense 20000 x 502 product -- 80 MB
  -- to keep 20000 numbers out of it. `block_leverage()` was still 41 per cent
  of that fit after the sparse route of 0.50.0; it is 2.7 per cent now, and
  the fit goes from 3.700 s to 2.420 s.

* Where the four shapes stand against the release that began this work:

  | | before | now | | mgcv |
  |---|---|---|---|---|
  | one smooth, gaussian, n = 8000 | 0.603 s | 0.268 s | 2.25x | 4.2x faster than us |
  | one smooth, Gamma, n = 8000 | 1.670 s | 0.980 s | 1.70x | 4.9x faster than us |
  | three smooths + random(40) | 14.700 s | 2.910 s | 5.05x | 1.4x faster than us |
  | random(500), n = 20000 | 24.310 s | 2.420 s | **10.05x** | **we are 70.4x faster** |

  On the mixed model the fitted quantities agree with `lme4::lmer` as they
  did before: sd 0.558434 against 0.558939, fixed effects to five figures.

# statmodels7 0.50.0

* The inner optimizer is offered the objective's EXACT second derivative --
  the information plus the penalty's Hessian -- which it was not being given.
  `minimize()` was called with `fn` and `gr` alone, so `newton()` built a
  numerical Hessian by differencing the gradient once per coordinate. It is
  passed to every method, as the gradient is; whether it is read is the
  method's business and a closure costs nothing until it is called.

  The gradient count is the tell, and the answers are identical throughout:

  | | `fn` | `gr` | `he` | time | |
  |---|---|---|---|---|---|
  | one smooth, `newton()` before | 44 | 234 | 0 | 0.520 s | |
  | one smooth, `newton()` after | 44 | 14 | 10 | 0.210 s | 2.48x |
  | three smooths + random, before | 78 | 3484 | 0 | 9.360 s | |
  | three smooths + random, after | 76 | 34 | 25 | 1.690 s | **5.54x** |

  `iwls()`, `bfgs()` and `lbfgs()` are unchanged to the millisecond, none of
  them reading a Hessian, which is the control. With the exact one
  `newton()` becomes the fastest inner method: 0.210 s against `iwls()`'s
  0.280 on the first shape and 1.690 against 2.960 on the second.

  **The default was never affected**: `iwls()` fits through `iwls_pieces()`
  and has always used the exact expected or observed information. The gap was
  reachable only by naming an optimizer.

* The leverage diagonal `(X_a M_ab X_b')_ii`, which every trace against `M`
  reduces to, is computed over the NONZEROS of each row where the design is
  sparse. A grouping indicator puts one nonzero per block in a row, so the
  quantity is a quadratic form over a handful of entries where the dense
  route computed `p_a p_b` of them per observation to keep one.

  The gate is measured, not assumed: 14.2x on the operation at a combined
  density of 3.6e-05, 50x SLOWER at 0.18 and again on a dense block, R's
  per-element indexing being far dearer than a BLAS flop. Interpolating puts
  the crossover near 1.1e-03, which is the threshold, with an absolute cap on
  the number of pairs besides -- the density gate bounds the ratio of work and
  says nothing about its size.

  End to end at 20000 observations with a random intercept, the two routes
  forced either way: 1.17x at 200 levels, 1.89x at 500, 1.26x at 1000 and
  1.11x at 2000. The win peaks and falls back, which says where the
  bottleneck moves: past about a thousand levels the dense `O(p^3)` inverse
  takes over, and at 503 coefficients it had been 4.4 per cent of the fit.

* The exact gradient and the exact Hessian are computed WHEN THE SEARCH ASKS
  FOR THEM, not whenever the criterion could supply them. Which of the two a
  search reads varies -- `nelder_mead()` reads neither, `lbfgs()` and `bfgs()`
  read the gradient and never the Hessian, `newton()` reads both -- and a line
  search evaluates the objective at many trial points where it wants no
  derivative at all. Measured on the counts, the inner objective is asked for
  a value three to four times more often than for a gradient.

  This needs no predicate about the optimizer and follows what it does rather
  than what its class declares. Measured, with every fit landing on the same
  answer:

  | | before | after | |
  |---|---|---|---|
  | one smooth, `newton()` | 0.350 s | 0.260 s | 1.35x |
  | one smooth, `lbfgs()` | 0.230 s | 0.170 s | 1.35x |
  | one smooth, `nelder_mead()` | 0.720 s | 0.420 s | 1.71x |
  | three smooths + random, `newton()` | 4.790 s | 2.580 s | 1.86x |
  | three smooths + random, `lbfgs()` | 9.440 s | 3.290 s | 2.87x |
  | three smooths + random, `nelder_mead()` | 116.420 s | 38.140 s | 3.05x |

* Where the four shapes stand against the release that started this work, and
  against mgcv fitting the same model with the same basis and the same
  criterion:

  | | before | now | | mgcv |
  |---|---|---|---|---|
  | one smooth, gaussian, n = 8000 | 0.603 s | 0.280 s | 2.15x | 4.4x faster than us |
  | one smooth, Gamma, n = 8000 | 1.670 s | 0.995 s | 1.68x | 4.9x faster than us |
  | three smooths + random(40) | 14.700 s | 3.000 s | 4.90x | 1.4x faster than us |
  | random(500), n = 20000 | 24.310 s | 3.700 s | 6.57x | **we are 46.1x faster** |

* An INNER context was measured and NOT built. The objective's `fn`, `gr` and
  `he` are already separate closures, so an optimizer that reads no Hessian
  never pays for one -- measured, `he` is asked for zero times per fit under
  every inner method. What they do repeat at one point is the linear
  predictors, and that is 15 to 16 per cent of one `fn` plus `gr`; since a
  value is asked for three to four times more often than a gradient, most
  evaluations have no partner to share with, and the whole saving is about 7
  per cent of the inner objective. It does not pay for the machinery.

* The criterion, the gradient, the effective degrees of freedom and every
  fitted quantity are unchanged throughout: bit for bit on all four shapes at
  both a fitted and an off-optimum point, the Hessian to 2.2e-14 from the
  change of summation order in 0.49.0.

# statmodels7 0.49.0

* The marginal criterion, its exact gradient and its exact Hessian read ONE
  evaluation context per point instead of each rebuilding what it needs. The
  information used to be assembled three times at one point and the penalized
  matrix factorized twice, and the Hessian additionally called the gradient,
  which repeated the whole of it a fourth time -- measured by `Rprof`'s
  `by.total` on a random intercept over 500 levels, the gradient and the
  Hessian together accounted for 128 per cent of the fit, the overlap being
  exactly that repetition.

* A contraction of a third or fourth derivative is never assembled where it is
  only traced against `M`. `tr(M X'WX)` is a weighted sum of the
  per-observation diagonal of `X M X'`, which the gradient already computes,
  so `contract4` and the per-pair `contract3` leave the Hessian's pair loop
  entirely. Measured on the operation at 8000 observations and 69
  coefficients, forming the matrix and tracing costs 25.5 ms where the sum
  costs 0.031 ms; at 20000 and 503, 3480 ms against 0.066.

* Measured end to end, with the criterion, the gradient, the Hessian, the edf
  and every fitted quantity unchanged (the criterion, the gradient and the
  edf bit for bit, the Hessian to 2.2e-14 relative, the summation order having
  moved):

  | | before | after | |
  |---|---|---|---|
  | one smooth, gaussian, n = 8000 | 0.603 s | 0.380 s | 1.59x |
  | one smooth, Gamma, n = 8000 | 1.670 s | 1.080 s | 1.55x |
  | three smooths + random(40), n = 8000 | 14.700 s | 5.220 s | 2.82x |
  | random(500), gaussian, n = 20000 | 24.310 s | 10.320 s | 2.36x |

* Two changes were measured and NOT made. Carrying a block-supported quantity
  without its surrounding zeros is worth 0.023 ms a call at 69 coefficients,
  which is 5 ms of a 14.7 s fit, and at 503 the subsetting is 0.4x -- slower
  than the full product -- because the penalty covers 500 of the 503 columns.
  And computing `tr(M K_l M K_m)` through `Z = X U'` measured 0.8x and 0.5x:
  forming `Z` costs `O(np^2)` just as forming one `X'WX` does, and with one to
  ten pairs there is nothing to amortize it over.

* `aic()`, `bic()` and `cv()` are untouched, and deliberately: in
  `statmod_pe_derivs()` the same contractions enter genuine matrix products
  (`A_ml %*% PH`, `PA_m %*% (P %*% B_l)`) and not only traces, so the assembly
  is doing work there rather than being discarded.

# statmodels7 0.48.0

* A penalty may now answer with a sparse Hessian -- a smooth repeated over
  the levels of a factor does -- and the coercion is written ONCE, where the
  two kinds meet in `statmod_penalty_at()`, rather than at each of the
  eighteen places that read its result and the four that read
  `penalty_dhessian()`.

* The accumulator stays a base matrix, which is measured rather than timid.
  Making it the design's own kind saved 0.8x end to end (17.0x against
  16.2x), and it would have put twenty-two consumers on a contract only the
  two a sparse design happens to exercise are covered by. What the blocked
  penalty buys is its CONSTRUCTION and its own storage, and neither passes
  through this accumulator.

* Measured end to end on 120 levels and 3000 rows, `s(x, by = g,
  sparse = TRUE)` against the dense spelling: 26.17 s against 424.86 s and
  0.55 MB against 20.44, with the log-likelihood identical to ten digits, the
  same lambda, the same effective degrees of freedom, `vcov()` agreeing
  exactly and `predict()` to 2.2e-16.

# statmodels7 0.47.0

* `statmod()` takes `linpar_control`, as `linpar_options()` returns it:
  `sparse` for the storage of the unpenalized parametric block and
  `contrasts` for the coding of its factors. It reaches the IMPLICIT
  `linpar()` term, the one the bare covariates collapse into and which a
  caller never writes; a `linpar()` written out takes them directly.

* The specification CARRIES them, so a rebuild reproduces the storage: a fold
  of `cv()` that built a dense design where the fit built a sparse one would
  be paying for a storage the model did not ask for.

* Measured on 150 groups, the parametric block goes from 1.54 MB to 0.12 MB
  with an identical log-likelihood and identical coefficients, and the fit is
  faster rather than slower.

* There is NO rescaling among the options, and that is measured rather than
  omitted. Scaling the columns and carrying the coefficients back is the
  remedy for a conditioning that squares, which is what forming `X'X` does;
  `iwls()` fits through a QR and never forms it. On columns spanning fifteen
  decades the raw fit and the scaled one converge in the same number of
  iterations and both agree with `lm()` to 1e-14. What moves is the SCORE the
  fit reports, 1.5e+02 against 9.2e-05, and the final verdict already
  arbitrates that on a dimensionless scale -- so the argument would have
  changed a printed number and not an answer.

* The argument and the function are named differently on purpose. With one
  name for both, the argument's default `linpar_options()` resolves to the
  argument's own promise and R reports "promise already under evaluation" --
  the shape section 7 records for an argument named after a class, met here
  for one named after a function. `glm(control = glm.control())` keeps them
  apart for the same reason.

# statmodels7 0.46.0

* A fold of `cv()` carries a term's matrix input onto its own rows.
  `data.frame(X = X, y = y)` SPLITS a matrix into `X.x1 ... X.xp`, leaving no
  column `X`, so `lasso(X)` reaches past the data to the matrix in the calling
  environment -- `interpret_formula()` evaluates the call as
  `eval(call, data, env)` and looks a name up in the data first, in the
  formula's environment after. The fit was right, the matrix being captured
  once, and every fold then failed to rebuild because the name still resolved
  to all the rows.

* The matrix is already on the built term, and `term_build()` checked at the
  full fit that it has one row per observation, so a fold's rows are the same
  rows by position. `cv_bind_inputs()` binds the subset as a column of the
  fold, which builds for the fold the spelling the documentation asks the
  caller for. Measured, the two spellings now give identical coefficients, an
  identical criterion and an identical lambda on the same folds.

* Nothing is relearned that should be. A matrix carries no knots, no contrasts
  and no levels, so subsetting it and re-evaluating it give the same block; a
  FORMULA input is untouched and keeps being rebuilt on the fold's own rows,
  which is what that rule exists for. A sparse input stays sparse -- the
  column the fold is given is still a `dgCMatrix`, not its densification.

* The test fold is bound too: `term_predict()` evaluates a matrix input's
  expression in the new data with `baseenv()` as its enclosure, so a name that
  is not a column there is not found at all.

# statmodels7 0.45.0

* `n_values` and `min_ratio` leave `cv()` and `OuterMethod`, as `search` did:
  what a PATH does is not the criterion's business. The same criterion is put
  to every hyperparameter of a model, and a smooth one is read at the mode
  rather than swept, so three of the criterion's arguments described something
  most of what it was asked about does not have.

* They live on the term's own signature, where a reader can see the number:
  `lasso(x, n_lambda = 25, min_ratio = 1e-4)`,
  `enet(x, n_lambda = 25, n_alpha = 5)`. `path_fallbacks()` is what remains
  here, and it is reached only by a term that declares a kinked penalty
  without offering an argument for the grid -- `random()` under a Laplace
  prior is the case.

# statmodels7 0.44.0

* `search` leaves `aic()`, `bic()` and `cv()`. Whether a term's own kinked
  hyperparameters are covered by a product or one at a time is the TERM's,
  `enet(X, search = "cyclic")`, for the reason the whole enumeration is: a
  criterion is put to every hyperparameter of the model, and a smooth one is
  read at the mode rather than swept, so most of what it was asked about
  could not use the argument. `statmod_search()` reads it, beside
  `statmod_grid_size()` and `statmod_min_ratio()`.

* Per term is also what keeps one term's choice off another's:
  `y ~ lasso(X) + enet(R, search = "cyclic")` sweeps the elastic net one
  coordinate at a time and leaves the lasso alone.

# statmodels7 0.43.0

* A hyperparameter the readable block does not DESCRIBE keeps its own row.
  The quantities of a multivariate Student t prior are the standard
  deviations and the correlations of its scale matrix, and its degrees of
  freedom are none of those, so replacing the coordinate rows wholesale
  dropped `nu` from the summary. The question is put to the Jacobian -- a
  column that is zero throughout is a coordinate no quantity depends on --
  so a family that declares more later is covered without an edit.

* Measured across effect distributions that carry a shape: `skewnormal1`
  reports `sigma` and `alpha` with standard errors and intervals,
  `student_t1` reports `sigma` and `nu`, and a log-transformed gamma reports
  the variance of the gamma underneath at 0.2202 (se 0.0396, interval
  0.155 to 0.313) against the 1/a = 0.25 the effects were drawn from. Where
  an interval is absent the cause is the point and not a missing derivative:
  a shape escaping towards a limit leaves an outer curvature of the wrong
  sign, and no interval follows from it.

# statmodels7 0.42.0

* `summary()` reports a correlated random effect by the quantities it is
  about. Where a penalty answers `penalties7::penalty_readable()` the
  coordinate rows are replaced by the standard deviations and correlations of
  the effects, with a standard error from the delta method -- composing the
  penalty's Jacobian with the link's, the variance matrix being on the free
  scale the criterion was maximized on -- and each interval built on the scale
  the quantity declares and mapped back: log for a standard deviation,
  Fisher's z for a correlation. A standard deviation therefore cannot be given
  a negative lower end and a correlation cannot be given one that leaves
  (-1, 1), which is the rule every other interval in the toolkit follows.
  No test is printed, the null a z would report on being that the quantity is
  zero, which for a standard deviation is the edge of its range.

* Measured against `lme4` on random slopes over 40 groups, the effects'
  standard deviations and correlation are (1.2370, 0.4820, -0.1399) against
  (1.2395, 0.4898, -0.1373).

# statmodels7 0.41.0

* Nothing in the fitting layer changed: a random effect whose prior is read
  blockwise is one penalty entry like any other, and one per within-group
  column is the enumeration this layer has run on since 0.12.0. What changed
  is what a fit reports -- the hyperparameter of a gaussian random effect is
  now the standard deviation of the effects rather than a precision, so the
  tests that constructed or interpreted its value were turned round with it.

* Measured against `lme4` on a random intercept over 40 groups: the fixed
  effects agree to 1.3e-05 and the standard deviation of the effects is
  1.14197 against 1.143856. With random slopes the fixed effects agree to
  1.0e-03 and the effects' standard deviations and correlation are
  (1.239, 0.568, 0.448) against (1.242, 0.575, 0.442).

# statmodels7 0.40.0

* `aic()`, `bic()` and `cv()` take `search`. A term carrying several
  hyperparameters with a kink has every combination of them visited under
  `"grid"`, the default, and one swept at a time under `"cyclic"`. Between
  two terms the search alternates either way, so `y ~ lasso(X) + enet(R)`
  costs the two blocks added and not multiplied, and two elastic nets in one
  formula are 100 + 100 points rather than 10^4.

* The product is the default because a cyclic sweep traverses a cross through
  the point in hand and can stop where each coordinate is separately best
  without being jointly so. It costs `n_lambda * n_alpha` fits against
  `n_lambda + n_alpha` per pass; with three or more estimated it grows
  exponentially, which is what `"cyclic"` is for.

* The grid is not a rectangle, and each axis is built at the settings of the
  axes outside it rather than fixed in advance. For the elastic net the kink
  is `lambda*alpha`, so every alpha carries its own lambda axis descending
  from its own `lambda_max = kink/alpha`; for SCAD and MCP the shape leaves
  the kink alone, so `lambda_max` is one number whatever the shape and those
  two axes really are a rectangle. Both follow from carrying the size of the
  kink back onto the hyperparameter where it stands, and a test pins the
  second, which is the controproof against writing the elastic net's relation
  for a family that does not obey it.

* A hyperparameter a term wrote out as a vector is visited as it stands. The
  value emptying the block does not cap it and `min_ratio` does not extend
  it, both of those being ways to build a grid, and the order it is walked in
  is the penalty's -- from the emptiest fit towards the fullest -- so the warm
  starts run as they do on a built grid.

* A shape parameter is swept above the smallest value at which the block can
  be FITTED. SCAD's proximal operator needs `t < a - 1` and MCP's
  `t < gamma`, tightened to `t*d^2` under standardization's diagonal map, and
  `t = 1/sum(w x^2)` is a property of the data: measured, a standardized
  penalty on a column of spread 20 needs `a > 3` where SCAD's own bound is 2,
  and a Poisson block whose fitted means are near 1e-3 needs `a > 11`. The
  limit is asked of the penalty and bisected rather than written out, so
  neither constant appears in the path.

* The size of the kink is inverted in closed form. Measured over four
  decades, it is exactly a power of each hyperparameter -- one for the lasso
  and for the elastic net in both of its own, minus one for a Laplace prior
  written by its scale, zero for the shapes of SCAD and MCP -- with a spread
  of at most 5.6e-16, so the exponent is read once per grid and every value
  follows. The relation is measured and then checked at the values it
  produced, and a penalty that does not obey one falls back to bracketing.
  Measured, a bracketing solve cost 4.18 ms against a fit's 62.5 ms, so a
  path of twenty-five points spent 6.7 per cent of itself locating the values
  it would visit.

* A pass that would visit the points just scored is not run. The top of the
  path is refreshed between passes because the rest of the model moves, and
  where it has not the grid is the one already in hand.

* `history$outer` gains `setting`, the rest of the combination each point
  belongs to. One row is still one point, and `name` and `value` still carry
  the axis the path descends.

# statmodels7 0.39.0

* The top of a path is read again at every sweep. It is the size of the
  kink that empties the block, taken from the score at the coefficients in
  hand, so it moves with the rest of the model: a smooth beside a kinked
  term has its own smoothing parameter estimated INSIDE each point of the
  path, and every such fit changes the score it is read from. Reading it
  once left the path anchored to the state the search began in. Measured on
  an elastic net, the largest lambda of the path moves from 37.0307 to
  46.2884 between the first sweep and the second.

  It bites where there is more than one sweep, which today means more than
  one hyperparameter under estimation; a single kinked hyperparameter still
  makes one pass and reads the top once.

# statmodels7 0.38.0

* The path visits as many values as the TERM asked for, per
  hyperparameter: `lasso(~x, n_lambda = 50)`, `enet(~x, n_lambda = 40,
  n_alpha = 12)`. How finely a hyperparameter is swept is a property of the
  term for the same reason as whether it is swept at all, and the criterion
  applies to every term of the model at once and cannot know which it is
  looking at. `statmod_grid_size()` reads it from the penalty's entry and
  falls back to the criterion's `n_values` where the term said nothing.

* And as far DOWN as it asked: `lasso(~x, min_ratio = 1e-6)` sets the
  fraction of the emptying kink the path descends to. One number per term,
  since only the sweep by kink size uses it. `statmod_min_ratio()` reads it
  and falls back to the criterion's.

# statmodels7 0.36.0

* Which hyperparameters are estimated is said by the TERMS, and by nothing
  else. Every one is estimated unless the term that carries the penalty
  holds it: `lasso(x, lambda = 3)`, `ridge(x, sigma = 0.5)`,
  `enet(x, alpha = 0.5)`, `scad(x, a = 3.7)`, `s(x, lambda = 2)`,
  `te(x, z, lambda = c(1, 5))`, `random(~1 | g, hyper = c(sigma = 0.4))`.
  The term is where the penalty is named, so it is where that belongs.

* `statmod(hyper = )` is REMOVED, and so is `over` on `aic()`, `bic()` and
  `cv()`. Both said the same thing as the terms, and whichever of the two a
  reader believed, the other was read by nobody when they disagreed. Passing
  `hyper` now signals an error naming the spelling that works.

* A hyperparameter no path could reach is swept over a grid of its own. The
  elastic net's `alpha` scales the kink and is bounded by one, so no
  admissible value of it empties the block; the shape of SCAD and MCP has no
  upper bound and does not move the kink at all. Neither is reachable by the
  geometric path over kink sizes, and leaving them out was how `alpha`
  stayed at 0.5 while the summary said a criterion had chosen it. A bounded
  one is swept over its interval and an unbounded shape over a geometric
  grid above its lower bound, which spans the values the literature uses
  (3.7 for the SCAD of Fan and Li, 3 for the MCP of Zhang).

* The summary reads what was held from the terms, so a hyperparameter is
  marked fixed if and only if a term fixed it.

# statmodels7 0.35.0

* A bounded hyperparameter is reported as held, not as chosen. `ifelse()`
  returns a result the length of its TEST, and the test was a scalar, so a
  penalty carrying two hyperparameters got one answer recycled over both:
  the elastic net's `alpha`, which no path varies, was marked as chosen by
  the criterion that had chosen its `lambda`.

* `aic()` and `bic()` take `over`, as `cv()` already did, and a bounded
  hyperparameter named there is now really swept. It used to be accepted
  and ignored: a path walks the SIZE OF THE KINK, from the value that
  empties the block down, and no admissible `alpha` empties it at a given
  `lambda`, so every point was dropped and `alpha` came back at its default.
  A bounded hyperparameter is swept over its own interval instead, the
  endpoints excluded because the bounds are open -- an elastic net at
  `alpha = 0` has no kink at all -- and the sweeps being cyclic, the kink
  still moves through `lambda`. Measured on ten columns of which seven
  carried signal, `bic(over = c("lambda", "alpha"))` moves alpha off 0.5 to
  0.615 and the criterion from 322.46 to 320.38.

  The default still holds it, as `glmnet` holds its `alpha` and `ncvreg` its
  `gamma`, and the end-of-path warning is not raised for such a sweep: the
  grid covers the whole interval, so an answer at either end is the
  interval's and there is nothing to widen.

# statmodels7 0.34.0

* A hyperparameter chosen by `sparse_criterion` is reported as estimated by
  the criterion that chose it. `summary()` asked `methods$outer`, which
  carries the marginal criterion alone, so a lambda a path had selected was
  printed `(fixed)` beside a note saying it was held at the value it was
  given -- of a number the caller had never seen. The fit now records which
  criterion swept the kinked penalties and exactly which hyperparameters it
  reached, and the cell says `(bic)`, `(cv)` or `(fixed)` accordingly.

* The hyperparameters come first in every penalized block. They govern
  every coefficient under them, and a table opening with a hundred selected
  coefficients buried the one number that produced the selection.

* A hyperparameter estimated by `reml()` or `ml()` carries a standard error
  and an interval. It maximizes a criterion that is twice differentiable in
  it, so its variance is the inverse of the negative of that criterion's
  Hessian at the point reached, which `statmod_marginal_hess()` already
  computed exactly; `statmod_hyper_vcov()` returns it. The interval is
  built on the free scale its link defines and mapped back, as every other
  interval in the toolkit is, so a positive hyperparameter keeps a positive
  lower end. No test is printed beside it: the null a `z` would report on
  is that the hyperparameter is zero, which is the edge of its range rather
  than an interior hypothesis. Validated against the criterion itself,
  refitted at each probe and differenced twice on the free scale.

  A hyperparameter chosen by a path -- `aic()`, `bic()`, `cv()` over a
  kinked penalty -- still carries none, and that is a refusal rather than
  an omission: it is the argument of a minimum over a grid, not the root of
  a derivative, so there is no curvature to read and its uncertainty is a
  resampling question.

* A path that cannot score a single point signals an error instead of
  keeping the hyperparameter it came in with. `cv()` refits on each fold,
  and a term built from a matrix that lives in the calling environment
  rather than in `data` cannot be re-evaluated on a subset of it, so every
  fold failed, every deviance came back `NA`, nothing was selected, and the
  fit returned the DEFAULT hyperparameter reporting success -- a lasso at
  `lambda = 1` selecting every column. The error carries the message the
  refit failed with. A marginal criterion reaching a kinked row is
  unaffected: it scores nothing by construction, and the hyperparameter
  keeping its given value is the documented answer there.

# statmodels7 0.33.0

* A break-point term has a section of its own in `summary()`, beside the
  smooths, the random effects and the penalized blocks, and it reports the
  quantities of the model rather than the coefficients of the working
  block: the linear effect, the changes, and the break-points under the
  name `psi`. The standard error of a position comes from
  `modelterms7::term_readable()`'s Jacobian by the delta method.

  A break-point gets an estimate, a standard error and an interval, and no
  test. The null a `z` of estimate over standard error reports on is that
  the position is zero, which is not a hypothesis anyone holds; the one a
  reader wants is that there is no break-point at all, and under it the
  position is a nuisance parameter that vanishes, so the classical p-value
  is wrong by a factor of three to five. `segmented` prints the estimate
  and the standard error alone for the same reason.

* The default start asks each term where its own block begins, through
  `modelterms7::term_coef_start()`. It used to put every coefficient at
  zero except an equation's intercept, and for a term that recomputes its
  block from its coefficients zero is a singular point rather than a
  neutral one: `y ~ jump(x, npsi = 3)` returned every coefficient exactly
  zero and reported failure, the three break-points having collapsed onto
  one clamped position. A term asking for zeros is left alone, so an
  ordinary block is unchanged, and `start_origin()` still means the
  origin.

# statmodels7 0.32.0

* The exact gradient of the marginal criterion reaches a penalty over a
  STRUCTURAL term's own parameters. `outer_gradient_ok()` no longer refuses
  `isTRUE(u$structural)`; it asks whether the term answers
  `modelterms7::term_third()`, read from the class the method is registered
  on, so a term written later is covered without a list of class names here.
  `reml()`/`ml()` on a score-driven panel therefore run `lbfgs()` on an exact
  derivative where they ran `nelder_mead()` on values.

  Where a filter is present the determinant spans the term's parameters too,
  so `K` carries `-sum_t w_t l_p E_t` and differentiating it along the
  direction the mode moves in gives three contributions instead of one:
  `u_vector()`'s formula with V in place of X (direction-free, computed
  once), the derivative of `V_p` itself, which is `E_t v`, and the derivative
  of `l_p E_t`, whose second half is the term's third derivative. The full
  third derivative is never formed -- it is contracted in one direction, so
  the cost is O(n m^2) per hyperparameter rather than O(n m^3) once.

  Validated against numDeriv on the criterion with the mode REFITTED at
  every hyperparameter: 2.7e-09 to 2.7e-07 relative under `reml()`, 3.3e-06
  under `ml()`, 6.3e-06 with a smooth beside it, and 1.2e-07 at ten and at
  twenty groups. The control, an ordinary smooth, sits at 3.9e-07 on the
  same harness.

  Measured in evaluations of the criterion, each a whole inner fit, against
  `nelder_mead()`, all three reaching the same criterion and the same prior
  scale to the printed digit:

        3 groups x 40      6 against 59      4.1 s against 14.8 s
        8 groups x 80     10 against 27     18.7 s against 33.2 s
       20 groups x 60      9 against 41     88.0 s against 88.9 s

  The last row is worth reading as it is: the evaluation count is what the
  exact route buys and it is 4.6x there, while the WALL TIME is a wash,
  because each exact evaluation also pays for the gradient (one
  `term_third()` and one extra `term_curvature()` per hyperparameter) and
  one point on that path is still unavailable and costs a long inner fit
  that does not converge.

* An outer step the search TAKES BACK now leaves no trace, where until now
  only half of it did. The coefficients were already protected --
  `state$beta` is written after a point is known to be usable -- but a
  structural term's own parameters live in the design's structural state,
  which is an environment the inner fit writes into as it goes, and
  `statmod_fit_structural()` stores its optimizer's last point whether that
  optimizer converged or not. An unavailable point therefore MOVED the
  filter's parameters, and the next evaluation started from wherever the
  failure had left them.

  It ratchets, and the exact gradient is what made it reachable: with a
  derivative-free search the steps are small and the state never leaves the
  basin, while `lbfgs()` takes its first step as the full gradient, which on
  a panel of twenty groups is -12.3 on a log-scale hyperparameter (the
  criterion's derivative grows with the number of penalized coordinates), so
  the first trial prior scale is 4.6e-06. Measured before the fix: the
  search then reported THIRTY consecutive unavailable points as the line
  search halved eta back to 1e-08, including points beside a start that had
  converged at the first evaluation, and gave up at 620 s with no fit. After
  it: nine recorded evaluations (ten attempts, the one unavailable point
  recovered from), and the same answer `nelder_mead()` reaches in 41 --
  prior scale 0.3783702 against 0.37837, criterion -1696.89343 both.

  ⚠️ Two explanations were measured and refused before this one. The
  criterion is NOT noisy here: repeated at one hyperparameter from a running
  warm start, and arrived at from eight different previous hyperparameters
  as far as six units away, its spread is exactly zero. And the warm start
  is not poisoned by VISITING the degenerate point: stepping to eta = -12.3
  and back returns the identical criterion. What ratchets is the structural
  state left behind by a failed inner search, which no snapshot protected.

* A filter beside a RANDOM EFFECT fits. `statmod_full_information()` lays
  each equation's design into its own columns of a row spanning every
  unknown, and an equation carrying a random effect has a sparse design:
  writing a `dgCMatrix` into a slice of a base matrix is a LENGTH ERROR
  rather than a conversion, so `y ~ x + random(~1|id) + gas(...)` stopped
  with "number of items to replace is not a multiple of replacement length".
  It is the sixth round of the shape the sparse conversion keeps taking, and
  it was NOT reached by the exact gradient -- the failure is in the joint
  inner step, which has called that assembly since 0.23.0, so the
  combination has never fitted. The block is densified at that one
  assignment and nowhere else: the matrix is one column per unknown with a
  single block filled, it is read one ROW at a time by the callback beneath
  it (the access a compressed-column matrix is worst at), and the fit, the
  penalties and the solve all go on seeing the sparse design.

* `deriv4_key()`, `joint_design_rows()`, `structural_joint_basis()`. The last
  is written once and used by both `statmod_marginal_full()` and the
  gradient: two callers composing the integrated subspace separately would
  agree only by accident.

# statmodels7 0.31.0

* `reml()` and `ml()` default to the OBSERVED information (Giovanni). That is
  what turns on the exact outer gradient and Hessian, so the search becomes
  `newton()` instead of `nelder_mead()`. Measured in evaluations of the
  criterion, each a whole inner fit: 19 against 31 with one hyperparameter --
  where it does not pay, a simplex needing no derivative in one dimension --
  then 126 against 35 with a smooth and a random effect, 133 against 32 with
  two smooths, and 166 against 6 with a modeled scale, the criterion
  agreeing to four decimals throughout.

  ⚠️ It has a cost, and it is on the criterion's PLATEAU. Once the penalty
  dominates the information the criterion is flat -- its kept normalizing
  constant cancels the Laplace determinant -- and with the exact Hessian
  Newton walks the whole of it: on a response scaled by 1e4 the smoothing
  parameter reaches 9.5e-20 where the simplex stopped at 3.2e-8, the fitted
  function identical to five decimals. The fit is unharmed, the reported
  hyperparameter is not meaningful there, and `modelterms7::edf()`'s plain
  `solve()` cannot invert `H + S` that far out. Three tests now ask for
  `reml("expected")` where what they check is their own reference or the
  inner rule rather than the outer surface.

  This does not reach a penalty over a STRUCTURAL term's own parameters,
  where the exact-gradient route rejects for its own documented reason: the
  contraction assumes the predictor is `X beta` and a filter's level is a
  recursion. Such a search is derivative-free whatever the information.

* The verbose trace is formatted. The method that will do the work is named
  on the right of every rule, blocks are separated by blank lines, and the
  hyperparameters' names -- which are long, a term's key being the call that
  produced it -- are said ONCE in a legend while the lines carry the values.
  A sparse block's line names the route that actually ran, coordinate
  descent or proximal gradient, rather than one guessed at beforehand:
  whether the compiled descent applies depends on the operator being
  admissible at the step the fit takes.

* An outer point where the criterion is unavailable is now a FINITE barrier,
  strictly worse than everything seen, instead of an infinity. The
  difference matters to a search that differences its own gradient: the
  probe lands in the unavailable region, the difference comes back
  non-finite, the direction is meaningless and the line search stops.
  Measured on a panel with a ridge over a filter's own parameters, the outer
  BFGS died at its 73rd evaluation with the criterion still improving; it
  now runs past it to 132, carrying the hyperparameters from (1, 1, 1) to
  (2.07, 1.44, 1.30) and the criterion from -1459.13 to -1423.0.

* The trace names a term instead of repeating its whole specification.
  A key is the call that produced the term, and one line of an outer trace
  carried a deparsed `gas(p = 1, q = 1, time = t, by = ~ridge(~id), links =
  list(...))` three times over. `short_keys()` keeps the leading call's
  first argument -- so `s(x, k = 20)` and `s(z, k = 8)` stay apart -- and
  everything after `::`, which is what distinguishes one entry of a term
  from another; where shortening would make two labels the same, none is
  shortened. The columns of `history$outer` keep the full key, where it has
  to stay unique and reconstructible.

* A fit records the optimizer that searched over its hyperparameters, in
  `methods$search`, the caller's where one was given and otherwise the one
  chosen from what the criterion can supply. A fit says what fitted it
  rather than leaving a reader to reconstruct the default.

* `statmod(start =)` takes a STRATEGY as well as a list of values:
  `start_intercepts()` (the default, now sayable by name), `start_zeros()`,
  `start_random()` and `start_search()`. A strategy answers the generic
  `start_at()` and is asked ONCE, before the alternation between the
  coefficients and the hyperparameters begins. That is the whole reason it
  is an argument of its own rather than an optimizer: handed to
  `inner_optimizer` as `chain(sa(), iwls())`, a search would rerun at every
  hyperparameter the outer criterion tried -- 46 times on an ordinary fit,
  the folds times the path inside `cv()` -- returning the same answer each
  time. `start = NULL` and a list of values are unchanged to the bit, which
  a test asserts.

* `start_search()` runs a global search on the LIKELIHOOD ALONE, penalties
  off: what a starting value has to get right is the basin, which the fit
  will not correct by itself, while the penalties enter when their
  hyperparameters are estimated and at the probe values represent nobody's
  choice. It searches only where the problem is not convex -- a structural
  term's own parameters, the blocks of terms that recompute their own design
  (`nl`, `seg` and the break-point terms), and each equation's intercept --
  because a smooth or a random effect is a convex block the scoring step
  reaches from anywhere, and searching over a thousand random-effect
  coefficients would spend the budget where it buys nothing. Kinked
  penalties are never searched: their hyperparameter has a known upper end
  and is swept by a warm-started path, which a random jump would both fail
  to improve on and destroy. A structural term's parameters do not live in
  the coefficient vector, so the search sets them into the design's
  structural state and leaves its best there.

  ⚠️ Measured, on every model tried it changes NOTHING, and that is the
  honest report rather than an omission: `nl()` reaches `nls()`'s answer
  from both a poor start and a good one, `jseg()` reaches the same
  likelihood from three starting positions, and a filter's fit is
  unchanged. The toolkit's structured initializations -- the intercept-only
  MLE, `term_start()`, the break-point grid -- already arrive. What the
  tests pin is therefore that it is SAFE (a convex fit is unchanged to
  1e-10) and that it reaches the coordinates it claims to; it is an escape
  hatch for a model where the default fails, not a routine improvement.

# statmodels7 0.30.0

* `iwls()` takes a `criterion`. A scoring step is not an optimizer and
  carries its own loop, but the rule that ends it is the caller's to
  choose: any `optimizers7` criterion drives the loop, `NULL` keeps the
  built-in rule. The state the rule reads carries the score PER
  OBSERVATION as its gradient, so `crit_grad(t)` and `tol = t` are the
  same rule and a threshold means the same at n = 10 and at n = 1e7 --
  asserted by a test running both and comparing the iteration count and
  the coefficients. Passing `tol` beside a `criterion` is an error rather
  than a silent choice between them, and a rule needing a stationarity
  measure -- which a scoring step does not compute -- is rejected at
  construction, as `optimizers7::check_criterion()` rejects it for a
  derivative-free method.

* `solve_pd()` estimates the smallest eigenvalue from LAPACK's condition
  estimator read on the Cholesky factor it needs anyway, instead of
  computing the whole spectrum: `rcond` is
  1/(||A||_1 ||A^-1||_1), so `rcond * ||A||_1` is 1/||A^-1||_1, which for
  a symmetric matrix lies between lambda_min/sqrt(p) and lambda_min. The
  estimate therefore errs on the SMALL side -- the verdict is conservative
  and can never accept a matrix the exact test would refuse -- and by a
  bounded factor: measured over dimensions 5 to 300 and condition numbers
  1e2 to 1e15 the ratio to the true smallest eigenvalue stays between 0.29
  and 0.75, where the two cases the test must keep apart are separated by
  some fifty orders of magnitude. Both sides are unchanged and still
  asserted: a smoothing parameter at 6e15 is accepted, two identical
  columns refused. Measured at p = 1022: **4.5x** on the operation
  (1.40 s to 0.31 s) and a fit carrying a random effect over 1000 groups
  from 4.00 s to 2.78 s. `src/lapack_rcond.cpp` is the wrapper and
  `src/Makevars` links R's own LAPACK and BLAS.

* The penalty's square-root factor is taken from the diagonal where the
  penalty is diagonal, instead of from a dense eigendecomposition. The
  eigenvalues of a diagonal matrix are its diagonal, so the two routes
  agree by construction and a test pins them together; what makes it
  worth having is how often the case arises, a ridge, a random effect and
  the Demmler-Reinsch penalty of `s()` all being diagonal, and that the
  factor is recomputed at every iteration of the scoring loop. Measured
  on a random intercept over 1000 groups at n = 20000: one factorization
  went from 0.63 s to under a millisecond and the fit from 6.52 s to 1.97
  s, **3.3x**, with the objective and every coefficient identical. With a
  smooth beside the random effect, 3.0x. A non-diagonal penalty -- a
  tensor's, an anisotropic sum -- still goes through the eigen route, and
  a test asserts that it does.

# statmodels7 0.29.0

* `solve_pd()` inverts through the eigendecomposition it already computes
  for its test, and the test takes an optional reference scale: a
  smoothing parameter a criterion legitimately sends to 6.4e15 separates
  the eigenvalues (min 29.7, max 6.4e15, the matrix strictly positive
  definite) without flattening a direction, and against `max(ev)` alone
  that read as singularity -- `vcov()`, `summary()` and every `edf`
  refused a CORRECT fit on a poisson smooth over weak signal. `vcov()`
  and `statmod_edf()`'s whole-model smoother pass the unpenalized
  information's diagonal as the scale; a genuinely flat direction (two
  identical columns) is still refused, being small against that scale
  too, and a test asserts both sides. The eigen-inverse also removes the
  Cholesky that could fail outright at such conditioning, each eigenvalue
  inverted exactly and the shrunk directions simply reporting variances
  near zero.

* The inner stopping rule keeps its absolute form through the loop and
  adds a DIMENSIONLESS reading to the final verdict: max_j |g_j|/(n s_p)
  with s_p = sqrt(median H_jj / n) per equation. A location equation's
  score carries the units 1/y, so on a response scaled by 1e-3 the
  absolute rule's floor sat at 1.37e-6 against the 1e-6 threshold and a
  run stalled AT the optimum reported failure; the dimensionless reading
  relabels exactly that run and never stops one, so every trajectory is
  unchanged. Two stronger designs were tried and REFUSED by measurement:
  a per-coordinate normalization by (H+S)_jj let penalized coordinates
  converge loosely at extreme shrinkage and moved every outer trajectory,
  and driving the loop with the per-equation form made the tolerance
  unreachable at y*1e4 -- the stall guard on the objective, whose
  magnitude grows with log y, fires before a rule 1/s_p stricter can, the
  inner reported failure across the whole corridor of smoothing
  parameters between the plateau and the optimum, and the outer search,
  reading those points as unavailable, never crossed it (1482
  evaluations, cor 0.82 against the truth where the shipped design
  reaches 0.998 at every scale from 1e-3 to 1e4).

* An UNHELD structural level starts at the equation's own data-based
  intercept rather than at zero -- the mirror of "the intercept wins":
  held, the intercept carries the response's scale; unheld, which is the
  volatility spelling `sigma ~ gas(...) - 1`, the level itself must
  absorb it. Started at zero on a response of another scale, the filter
  reads a score of y^2/sigma^2 at its first step and leaves the
  representable range before the search's guard can step back: measured,
  `sigma = NaN` at y scaled by 100, where the intercept spelling of the
  same model fits. With the start from `statmod_intercepts()`, the same
  volatility model now fits identically at y, 100 y and 10^4 y
  (persistence 0.813 at every scale). Only a fresh start is touched; a
  specification a fit has been through keeps the parameters it arrived
  at.

* An outer criterion that is unavailable at the STARTING hyperparameters
  is reported with its cause instead of dying inside the optimizer.
  `evaluate()`'s own comment already stated the intent -- a failure at
  the start is the caller's and is raised -- but implemented it only for
  the branch that raises; the branch where the inner fit CONVERGES to a
  point whose penalized information has no Cholesky factor (a
  degenerated parameter, measured on a gamma model whose nl rate ran to
  exp(31) with `converged = TRUE`) fell through to a non-finite value
  and the optimizer stopped with "the objective is not finite at the
  starting value", naming the point and not the reason. The error now
  says what is unavailable and why, and points at
  `outer_criterion = NULL` for inspection.

* The subformulas of a term's own parameters (modelterms7 0.27.0) fit end
  to end, and the layer needed almost nothing: the penalty enumeration of
  0.12.0 and the marginal route of 0.24.0 already reach a sub-term's
  hyperparameter under the key `term::parameter::subterm`. Measured, not
  presumed: before modelterms7 retired the shorthand,
  `gas(omega ~ random(~1 | id), by = id)` was pinned against
  `deviations = "omega"` with a ridge EXACTLY -- the REML hyperparameter
  to the printed digit (0.3668324), the coordinates one for one, the
  log-likelihood to 4 decimals -- and `nl(a ~ ridge(~g))` and
  `seg(x, psi ~ id)` fit with their hyperparameters estimated and the
  truth recovered. All of it is pinned in `test-submodels.R`, the panel
  case against a per-group recursion written by hand.

* `vcov()` and `summary()` died on duplicate row names for ANY fit whose
  penalty sits on a structural term's own parameters -- the shorthand
  `gas(deviations =, penalty =)` of 0.24.0 included, exposed only now
  because nothing called `vcov()` on such a fit. `coef_labels()` wrote a
  structural unit's `cols`, which index the term's parameter vector, into
  the per-column flags of an equation design that has no such columns; R
  grew the vector past the design and recycled the labels into duplicate
  rows. A structural unit is skipped there now, and both routes are
  asserted.

* The level's confounding question generalizes from a name to a subspace.
  With `omega ~ Z gamma` the constant coordinates of the development are
  held exactly as the scalar level was (`term_level_param()` names them,
  and the intercept then carries the stationary level, with its standard
  error). For the rest, an UNPENALIZED coordinate whose column the
  equation's design already spans is flagged with a warning through the
  new `modelterms7::term_level_design()`: holding it would change the
  model where the confounding is not exact, while a penalized coordinate
  is identified by its penalty, exactly as a deviation is.

* A structural term's starting values come from the term itself,
  `modelterms7::term_start()`: zero on the unconstrained scale is no
  longer "the model without the term" for a loading on the log chart,
  whose natural point is a loading of ONE.

# statmodels7 0.28.0

* The compiled coordinate descent reads a SPARSE block, and the last
  densification in the chain is gone.

  `coord_fit()` materialized the penalized block dense to hand the kernel an
  `arma::mat`. A coordinate descent reads one column at a time, so a
  compressed-column matrix is the storage the method wants rather than one
  it tolerates: column `j` of a `dgCMatrix` is a contiguous run of its own
  nonzeros, which is exactly the walk every step of this algorithm makes.

  The algorithm is written ONCE against a column accessor and instantiated
  twice, so the dense and sparse kernels cannot drift apart. They agree BIT
  FOR BIT rather than to a tolerance, and that is a property of the
  arithmetic and not luck: skipping a structural zero omits an addition of
  zero, which is exact. The tests assert `identical()` on the coefficients,
  the sweep count and the screening gradient, over both gradient routes
  (residual and covariance) and both densities.

  Measured at n = 20000 against the same design densified by hand, with the
  coefficients exactly identical:

  | | before | after | vs dense |
  |---|---|---|---|
  | lasso, p = 200, density 0.005 | 1.04 s | **0.23 s** | 14.4x |
  | lasso, p = 1000, density 0.001 | 3.61 s | **0.25 s** | **233.7x** |

  The kinked branch is now faster than the smooth one on the same block
  (0.25 s against 2.37 s at p = 1000), which is the right way round: a lasso
  screens its coordinates where a ridge solves the whole system.

  The `dgCMatrix` is taken apart in R and its slots passed to the kernel, so
  the compiled code needs no dependency on the \pkg{Matrix} C API.

  ⚠️ The dense path's residual is now accumulated column by column rather
  than row by row, a compressed-column matrix having no cheap row walk. The
  arithmetic is the same sum in a different order, so results move within
  rounding; every existing comparison is far above it (1e-10 against the
  proximal route, 1e-13 against glmnet) and the suite is unchanged.

# statmodels7 0.27.0

* `methods` is declared in `Imports`. `as_sparse()` has called `methods::as()`
  since 0.26.0 and nothing declared it, so `--as-cran` reported *"'::' or
  ':::' import not declared from: 'methods'"* -- a WARNING, which the CI
  action treats as a failure. The local suite could not see it. Sixth
  instance of the rule that any package named with `::` goes in the
  dependencies in the same edit as the code that names it.

* A kinked penalty composed with a grouping indicator ABORTED THE PROCESS.
  `coord_fit()` handed the compiled coordinate descent the penalized block
  sliced out of its equation's design, and where the equation also carries a
  random effect that design is a `dgCMatrix`; the kernel takes an
  `arma::mat`, so the conversion threw at the `.Call` boundary and R died
  rather than raising something a `tryCatch` could see. `y ~ lasso(~x) +
  random(~1|g)` was unfittable. The block's own columns are dense whatever
  the rest of the equation holds -- the sparse columns belong to the other
  terms and are never touched -- so they are densified at the slice. Sixth
  round of the sparse conversion, and the first where the failure was not an
  R error.

* `coord_fit()` asked whether a penalty has a proximal table AT A STEP OF 1,
  which is a different question from the one it meant. SCAD and MCP have no
  table past their convex region, and under a diagonal map the condition is
  `t < (a-1)/d^2`, so a standardized block on a column of spread 20 was
  refused at the probe and sent to the general proximal route, which then
  raised on the same condition. The probe is taken at a step short enough to
  answer only the question about the family; the step that will be used is
  asked below it, once the working weights are known.

* Standardization needs nothing here, which was the point of putting it on
  the penalty: a standardized term fits what a hand-standardized design fits
  (coefficients `s*beta` to 8.9e-16, fitted values to 3.6e-15), multiplying a
  column by a thousand leaves the standardized fit exactly where it was
  (3.6e-15) where it moves an unstandardized one, and a random effect's block
  stays a `dgCMatrix` through the terms, the assembled design (density
  0.064), the information and `vcov()` with a standardized lasso in the same
  equation.

# statmodels7 0.26.0

* A model carrying a grouping indicator is fitted SPARSE end to end, and the
  gain is an order of magnitude:

  | n | groups | before | after | speedup |
  |---|---|---|---|---|
  | 2000 | 50 | 1.38 s | 1.29 s | 1.1x |
  | 4000 | 100 | 3.45 s | 0.70 s | 4.9x |
  | 6000 | 200 | 15.13 s | 0.86 s | 17.6x |
  | 8000 | 400 | 92.81 s | 3.31 s | **28.0x** |

  The shape changes and not only the constant: 3.45, 15.13, 92.81 was
  accelerating, and 0.70, 0.86, 3.31 is nearly flat. The log-likelihood is
  identical to the last printed digit at every size and the fixed effects
  agree with `lmer` to 2.0e-05.

  It took two steps and the first alone was not enough. With the block sparse
  the fit was 0.9x to 1.4x: `crossprod` left the profile entirely, and
  `.Fortran` ROSE from 57.9 to 72.6 per cent, because `iwls()` still ran a
  dense QR. **A profile share is not a speedup: removing the cheaper half
  leaves the dearer one.** `sqrt_design()` now keeps sparsity, where it
  assembled dense zero blocks per pair per iteration, and `augmented_solve()`
  reaches `Matrix::qr()`: 695x at 100 groups and 75475x at 1000 on the
  augmented design alone.

  A sparse QR is a QR, so the conditioning property the augmented route
  exists for is kept exactly, and the choice against a Cholesky of the normal
  equations -- which times the same and squares the conditioning -- does not
  arise.

  ⚠️ `Matrix`'s generics must be IMPORTED: a package that only depends on it
  gets base `t()`, which reaches `t.default()` and reports that its argument
  is not a matrix, from inside the assembly. `bind_blocks()`,
  `zero_information()`, `as_dense()` and `as_sparse()` are where the two
  kinds meet; the full-inverse sites densify deliberately, the inverse of a
  sparse matrix being dense.

* `inner_method` is `inner_optimizer` and `outer_method` is
  `outer_criterion`. The first takes optimizers, like `outer_optimizer`; the
  second takes criteria. `iwls()` stays a statmodels7 object rather than an
  optimizers7 one for the reason coordinate descent does: it needs the model,
  not `fn` and `gr`.

  The joint step passed `optimizer = NULL` unconditionally, so
  `inner_optimizer` was ACCEPTED AND IGNORED for exactly the models it
  matters most for. Measured now that it is honoured: `lbfgs` is 2x on a
  three-parameter filter (5.08 s against 9.95) and `newton` is 2.8x on a
  33-parameter panel (4.58 s against 12.97, 14 iterations against 130), the
  log-likelihood identical throughout. The default stays Newton, which wins
  where the parameters are many.

* The hyperparameters are ESTIMATED by default. `outer_criterion` is
  `reml()`, and it applies to the SMOOTH penalties and comes into play if and
  only if the model carries one -- a property of the model and not of how the
  argument was written, so typing the default changes nothing. A model with
  no penalty, or one whose only penalty is kinked, leaves it unused. On a
  smooth the difference is real: the effective degrees of freedom go from
  9.01, which is no penalization at all, to 7.81.

  `hyper` still means HELD AT THESE VALUES and steps the default aside.

* `sparse_criterion`, `bic()` by default, chooses the hyperparameter of a
  kinked penalty along a path over its own values. Where a model carries both
  kinds the path is outside and the marginal criterion is estimated inside
  each of its points, so a smoothing parameter comes from REML and a lasso's
  lambda from BIC in one fit -- which one argument could not express.

  A prediction-error criterion for the smooth penalties nested inside that
  path is REJECTED: it scores the same quantity at two levels, and measured,
  every point of the path came back NA, so the path had nothing to choose
  between and the hyperparameter kept its starting value while the fit
  reported success.

* The top of a path now empties the block. The documentation described a
  check -- "a top whose fit is not empty is doubled until it is" -- that was
  not implemented, so the grid covered a nearly flat stretch of the criterion
  and the choice fell on its edge: on eight coefficients of which three
  carried real signal the top read 26.5 where the block first empties near
  500. `min_ratio` is 1e-4, which is glmnet's own ratio now that the top is
  comparable to its `lambda.max`.

  ⚠️ The warning that reports a choice at an end fired on the INDEX alone
  while its text claimed the criterion was still falling. With the top now
  emptying the block the criterion is FLAT across that stretch, so index one
  is a legitimate minimum and the message was naming a cause that was not the
  real one. It compares against the neighboring point now, and the two
  warnings the suite carried were both of that kind.

# statmodels7 0.25.0

* A structural term is reported under the names its literature uses, and
  under the quantities rather than the coordinates. For a score-driven term
  that is `omega`, `alpha1` and `beta1`: the last is the AUTOREGRESSIVE
  COEFFICIENT, taken through the Levinson-Durbin recursion, where the free
  coordinate is a partial autocorrelation and above `q = 1` a different
  number.

  The standard error is the delta method over the JOINT variance,
  `sqrt(J V J')` with `J` from `modelterms7::term_readable()`, not one entry
  of a diagonal: a coefficient reads the whole chart, so the covariance
  between its coordinates enters. The interval is built on the scale that
  keeps each quantity in its own set and mapped back, as every interval in
  the toolkit is, and a quantity that reads a held parameter is reported
  without one rather than with the variance of the rest.

  Nothing else moved: for a term whose coordinates are already its
  quantities the base method reports them on the parameter scale with the
  diagonal Jacobian of their links, which is what this function did before.

# statmodels7 0.24.0

* `reml()` and `ml()` reach a penalty over a structural term's own
  parameters, where they used to return a criterion that did not move.

  A marginal criterion integrates the quantities its penalty shrinks. Where
  that penalty is a term's own -- the deviations of a panel -- those
  quantities are not coefficients, so a determinant taken over the
  coefficients alone does not depend on the hyperparameter at all and the
  criterion is the penalized likelihood, whose maximum in a shrinkage
  parameter is at no shrinkage. `statmod_marginal_full()` spans the
  coefficients AND the term's free parameters, which is the order
  `statmod_full_information()` already carries them in, with the penalty's
  own Hessian in the tail. Nothing was derived: both pieces existed.

  Measured on a panel of three groups, a ridge on every deviation, the mode
  refitted at each value: the criterion moves from -157.361 at a prior
  scale of 0.01 through a maximum of -157.347 at 0.05 to -189.023 at 50,
  and the deviations follow, sd 0.0008 to 0.148. Against a Laplace formula
  assembled by `numDeriv` over the same unknowns it agrees to 1e-5.

  `ml()` integrates the penalized coordinates alone, through the penalty's
  own range basis, so a population value the penalty does not cover is
  profiled rather than integrated, exactly as an unpenalized coefficient
  is.

  The EXACT gradient rejects here and the search stays derivative-free.
  That route reads how the determinant moves with the mode through a
  contraction of the family's third derivative against the design blocks,
  an assembly that assumes the predictor is `X beta`; a filter's level is a
  recursion of the term's parameters and the same contraction is not the
  derivative there.

* `vcov()` and the summary's structural table report a term whose own
  parameters carry a penalty.

  Both inverted the unpenalized information, on the comment that a term's
  parameters carry no penalty. That is not a conservative choice: the
  deviations of a panel are identified by their penalty and by nothing
  else, a constant added to a population value and taken off every
  deviation leaving the filter exactly unchanged, so the matrix is singular
  along that direction. The table reported a missing standard error for
  every parameter of the term. `structural_penalty_block()` places the
  penalty in the tail, written once because the fit, the variance matrix
  and the criterion all need it and a block each of them placed for itself
  would agree only by accident.

* A panel with deviations fits again. `statmod_fit_joint()` reaches
  `modelterms7::term_curvature()`, which rejected deviations, so every such
  model -- with or without a penalty -- signalled an error rather than
  fitting. It needs `modelterms7` 0.21.0.

# statmodels7 0.23.0

* A structural term of the filter shape is fitted in the SAME system as the
  coefficients rather than alternated with them.

  The alternation was never a statement about the model: the exact gradient
  of both blocks and the exact observed information over both together were
  already available, the second as `statmod_full_information()`, which was
  built for `vcov()` and discarded for the fit. What it cost is filter runs.
  Each sweep handed the term's parameters to an optimizer of their own -- and
  always `lbfgs()`, whatever the caller asked for -- whose every iteration
  re-ran the recursion and its adjoint with the coefficients held at a point
  that was about to move.

  `statmod_fit_joint()` runs one Newton step over the stacked coefficients
  and the term's free parameters. The unknowns are ordered as the information
  orders them, so no permutation is needed, and a level an intercept in the
  same equation carries is held and leaves the system exactly as it leaves
  the information. Measured on a panel of 25 groups and 750 observations:
  16.91 s to 3.21 s, a factor of 5.3, at a slightly higher maximum
  (-1248.7824347 against -1248.7825410). A term of the likelihood shape --
  `regime()` -- keeps the alternation, its information being assembled by a
  different route.

* A penalty over a structural term's own parameters is read from the term's
  own vector.

  `statmod_penalized()` looked its positions up in the design, which a
  structural term does not have: `blocks[[term]]` was NULL and the positions
  then indexed the equation's coefficients, 25 of them where the equation
  has one. The penalty was evaluated at `NA`, and everything built on it
  followed -- an objective that was not finite, an outer criterion reported
  as unavailable, and an inner run that spent its whole budget.

  Such a unit is now marked as structural and evaluated from the term's
  state, and `statmod_structural_penalty()` returns its derivative and its
  Hessian in the term's own parameters. The gradient of the structural block
  adds that derivative, which its objective had been including all along:
  two functions differing by a penalty are not each other's gradient, and
  optimizers7's own check had been saying so. Verified against `numDeriv` at
  three points, agreeing to the reference's own accuracy.

* The per-observation callbacks a filter runs reuse one buffer instead of
  rebuilding a list of scalars on every call. Worth 1.09x on its own; the
  number of filter runs, which the joint step addresses, is where the time
  actually was.

* `solve_pd()` decides positive definiteness from the eigenvalues rather than
  from whether `chol()` raised.

# statmodels7 0.22.0

* A term's effective degrees of freedom are its share of the WHOLE model's
  smoother matrix, which is the definition, rather than a block-wise
  approximation to it.

  `statmod_edf()` computed `tr[(H_bb + S_b)^-1 H_bb]` from each term's own
  block, dropping every coupling between that block and the rest of the
  model. The definition is `F = (H + S)^-1 H` over the stacked coefficients
  of every equation, with a term taking the trace of its diagonal block; the
  matrix is now formed once per call and every term reads its share off it.

  It was found verifying the counts on parameters other than the mean.
  Measured on a gaussian with a smooth in each equation, the two totals were
  16.98939 and 16.98885, and the whole of the gap sat on the MU smooth while
  sigma's agreed exactly: the Demmler-Reinsch basis is orthogonalized
  against the constant in the UNWEIGHTED metric, and the mean's information
  carries weights `1/sigma^2` that vary as soon as the scale is modeled, so
  the orthogonality the construction arranged does not survive the
  weighting. The gap is small wherever the blocks are nearly orthogonal and
  is not bounded in general.

  Now exact against the full-model trace, per term and in total, on a
  gaussian with smooths in both equations, on a gumbel, whose location and
  scale are NOT information-orthogonal, and on a gamma with its dispersion
  modeled. A kinked penalty keeps its own count, the number of coefficients
  away from the kink, since the curvature the smoother matrix is built from
  does not exist at a coefficient sitting on one.

# statmodels7 0.21.0

* `summary(fit, correct = TRUE)` carries what estimating the hyperparameters
  cost into the degrees of freedom.

  The ordinary effective degrees of freedom read the smoothing parameters as
  though they were known, and they were chosen from the same data, so every
  criterion built on the count is too generous. `statmod_edf_correction()`
  propagates their uncertainty into the coefficients by the implicit
  function theorem: with `J = -(H+S)^-1 d2rho/dbeta dtheta` from
  penalties7's `penalty_cross()` and `V_theta` the inverse of the outer
  criterion's own Hessian, the corrected count is `tr[(V_b + J V_theta J') H]`.

  It is the first of the two terms mgcv sums into `edf2`, and it reproduces
  mgcv's to 3.0e-04, 4.7e-04 and 2.6e-05 at n = 200, 400 and 2000 on a
  univariate smooth. mgcv's second term corrects for the Gaussian scale,
  which it profiles out of the fit; here every distribution parameter
  carries its own equation and its own coefficients, so that uncertainty is
  already inside the information.

  Nothing in it is smooth-specific: it enumerates penalties, so a
  `random()` effect is corrected by the same code path with no branch.
  Where no hyperparameter was estimated by a marginal criterion -- a kinked
  penalty, whose lambda `outer_hyper_index()` skips by construction, or a
  held one -- the correction is exactly zero and says so. That is a refusal
  rather than an approximation: the map from the hyperparameter to the
  penalized mode turns a corner whenever a coefficient joins or leaves the
  active set, and a delta method needs a derivative that does not exist
  there.

  It is off by default because it changes a number a reader may be
  comparing with an earlier fit.

# statmodels7 0.20.0

* A structural term's own parameters are COUNTED and REPORTED.

  A term that rewrites the predictor contributes no design columns, and the
  degrees of freedom were read off the columns, so `gas(1, 1)` was counted
  as ZERO and a model carrying four parameters reported two. Every criterion
  built on the total was that much too generous. What such a term spends is
  its own parameters, one apiece, less any level an intercept in the same
  equation already carries -- which is held rather than estimated and is not
  paid for twice.

  `summary()` now prints them with standard errors and intervals from the
  joint observed information, which spans the coefficients AND the term's
  parameters and was already being computed for `vcov()`. Each interval is
  built on the unconstrained scale and mapped back, so a persistence stays
  inside `(-1, 1)` and a gap stays positive; the standard error on the
  parameter scale is the delta method and is reported for reading. A held
  level is shown as held, with no standard error rather than a zero one.

* `logLik()` says which likelihood it is, and the other one is reachable.

  What it returns is the CONDITIONAL log-likelihood -- the density at the
  fitted coefficients, a penalized coefficient among them -- paired with the
  effective degrees of freedom. A criterion built on that pair is the
  conditional AIC, and `summary()` now labels it as such rather than
  printing a bare `AIC`. A mixed-model package reports the MARGINAL
  likelihood instead, integrating its random effects away and counting
  variance components, and its AIC is a different number answering a
  different question (Vaida and Blanchard, 2005). The two are not
  comparable and mixing their halves is what neither convention allows.

  `logLik(fit, type = "marginal")` returns the value `ml()` or `reml()`
  optimized while choosing the hyperparameters, with the estimated
  parameters as its degrees of freedom. On a random intercept at 3000
  observations it is -4372.79 on 4 degrees of freedom against lmer's
  -4371.71 on 4. Where no marginal criterion ran it is an error rather than
  a number that would look like one.

# statmodels7 0.19.0

* A structural term's callbacks go through `distributions7::distrib_kernel()`
  rather than through the derivative generics. A filter calls them once per
  time step, and each call was validating its arguments, aligning theta by
  name, dispatching, and assembling every component of its order to keep one
  of them.

  Measured on a gas(1,1) scale equation at 2000 observations, against
  rugarch on the same data: 123 seconds before this and the two changes
  below it, 17.7 after, with the log-likelihood identical to the digit at
  every size tried. `regime()` shares the callbacks and its standard errors
  are unchanged to the digit.

# statmodels7 0.18.0

* The callbacks a structural filter runs do their fixed work once instead of
  once per time step. `structural_callbacks()` recycled every distribution
  parameter to the sample length inside the closure and rebuilt the
  second-derivative component key with `match()` and `paste()`, both on
  every step of every filter pass. Worth about 12 per cent of a gas fit.

  It is not a change of order: `rep_len()` on a vector already of the right
  length returns it untouched, so the old form was linear in the number of
  observations like this one. What remains is the two generic calls per
  step, which is where the rest of the time is.

# statmodels7 0.17.0

* `vcov()` on a model carrying a `regime()` term reports the OBSERVED
  information. What it inverted before was the complete-data information,
  the ordinary one averaged over the smoothed states, which is the matrix a
  scoring step inverts and is not the information the fit carries: the two
  differ by the conditional variance of the complete-data score, so the
  complete-data one is the larger and every standard error read off it is
  too small. Measured on a two-regime gaussian at 300 observations, the
  matrices differ by 25 per cent and the worst standard error is understated
  by 30 per cent where the regimes overlap; the gap closes as they separate
  and the states become known, which is what the missing-information
  principle says it should do.

* `statmod_regime_information()` supplies the model's side of
  `modelterms7::term_hessian()`: how each equation's unknowns reach each
  predictor, and the family's first and second derivatives at the predictor
  each regime shifts to. Unlike a filter's, those cannot be looked up from
  one evaluation, so the family is evaluated once per regime, vectorized
  over observations.

  As for a filter, the matrix spans the coefficients AND the term's own
  parameters, and the coefficient block of the joint inverse is what is
  reported; a level an intercept already carries is dropped rather than
  estimated. Checked against `numDeriv` on the exact gradient away from the
  optimum, where a wrong Hessian and a right one differ.

# statmodels7 0.16.0

* `statmod()` fits `regime()`, which is the last term modelterms7 defines.
  Every one of them is now fittable.

  A term of that shape does not report a predictor: its contribution is a
  likelihood mixed over latent states, so `statmod_loglik_at()` takes the
  term's own and not the density at any single point. Everything else
  follows from Fisher's identity -- the derivative of that likelihood in ANY
  predictor is the posterior-weighted derivative of the ordinary one -- so
  the score differentiates the log-density once per state, vectorized, and
  weights by `modelterms7::term_posterior()`. Against `numDeriv`: 1.8e-07
  in the coefficients on a scale of 246, and 6.0e-09 in the term's own
  parameters. On 600 observations simulated from the model the fit returns
  levels (-0.033, 3.057) against a truth of (0, 3) with the scale at 0.973.

  The matrix the scoring step inverts there is the COMPLETE-DATA
  information, the ordinary one averaged over the smoothed states, which is
  what an EM step inverts. It is positive definite and it is not the
  observed information of the mixture, which needs Louis's identity; with an
  exact gradient a scoring matrix has only to be positive definite.

* Where a structural term's level and a linear intercept are both present,
  THE LINEAR INTERCEPT WINS and the term's level is held at zero. The two
  are exactly confounded -- shifting the intercept by `c` and the level by
  `-c(1 - sum b)` leaves every predictor unchanged -- and a fit reached that
  ridge reporting convergence, the score being small because the surface is
  flat. `y ~ x + gas(...)` and `y ~ x + regime(...)` are ordinary things to
  write again, and nothing is lost: what a constant cannot express is the
  dynamics, or the difference between one regime and another.

  Which parameter is the level is the term's own answer, through
  `modelterms7::term_level_param()`; which one is dropped is this layer's,
  since only it knows what else the equation carries. The question is asked
  of the SPAN of the equation's design and not of a column named
  `"(Intercept)"`, a factor coded without one summing to the constant just
  as well. The held parameter leaves the information as well as the fit,
  its row and column dropped from the joint matrix.

  Measured, and this is what says the two parametrizations are one model:
  with the level held and with it free the fits reach the same maximum
  (-553.793723 against -553.7937276) and the fitted intercept, 0.952936, is
  the other fit's stationary level `omega/(1 - b)`, 0.95295.

* `fit@structural` carries `held`, the parameters an intercept is carrying
  instead.

# statmodels7 0.15.0

* The observed information of a model carrying a filter is exact.
  `statmod_full_information()` assembles

      -d2l/du2 = -sum_t w_t sum_{q,r} l_qr,t V_q,t' V_r,t
                 - sum_t w_t l_p,t E_t

  over the coefficients of every equation AND the term's own parameters,
  with `E`, the second derivative of the predictor through the recursion,
  from `modelterms7::term_curvature()`. Only the filter's equation has a `V`
  that is not its own design; the third derivatives the second sum needs are
  distributions7's, in closed form for every family.

  What it replaces is the naive `X'WX`, which the scoring step still uses and
  is right to -- with an exact gradient any positive definite matrix
  converges to the right point -- but which is not the information. Measured
  against numDeriv on the exact gradient, AWAY from the optimum, where a
  wrong and a right Hessian are told apart: 4.3e-11 relative. The naive
  curvature differs from the exact one by 17.5 per cent of the coefficient
  block on a one-equation model and 35.3 per cent on a two-equation one, the
  second larger because the cross-equation feedback contributes there too.

* `vcov()` inverts the joint matrix and takes the coefficient block, rather
  than inverting the coefficient block alone. The two differ wherever the
  coefficients and the filter's parameters are correlated, which is the
  situation the filter creates; the joint route is what accounts for having
  estimated the term's parameters rather than known them.

* A filter's level and a constant in its own equation are rejected. They are
  EXACTLY confounded: shifting an intercept by `c` and the level by
  `-c(1 - sum b)` leaves every predictor unchanged, the recursion being
  affine in the level given the score path and the score path depending on
  the predictor alone. The fit reached that ridge reporting convergence and
  only the information saw it, so `reject_confounded_level()` asks at design
  time, where the terms can be named.

  It asks whether the CONSTANT LIES IN THE SPAN of the equation's design and
  not whether a column is called `"(Intercept)"`: a factor coded without one
  still sums to the constant, and the span test rejects it where a name test
  would not.

* No compiled code for any of it, and the measurement is the reason. One
  curvature costs 1.4 filter runs and is FLAT in the number of unknowns --
  five against fourteen changes nothing -- so the matrix arithmetic is not
  where the time is; the family's third derivatives cost 0.0000 s. What it
  is spent on is the per-observation callbacks, and there the win was in R:
  a curvature runs at a point already reached, so the score and the
  curvature of the density can be computed vectorized once and looked up,
  where the filter must evaluate them at a predictor it has just produced.
  That change took the curvature from 10.83 s to 3.95 s at n = 8000, a
  factor of 2.7, and the answer is bit-identical.

# statmodels7 0.14.0

* `statmod()` fits the terms whose block depends on their own coefficients:
  `seg()`, `jump()`, `jseg()` and `nl()`. They were rejected before, because
  assembling such a block once and solving for the coefficients estimates
  something else and reports convergence.

  The block is refreshed at the coefficients inside the objective, and the
  predictor takes the term's CONTRIBUTION rather than the block times the
  coefficients. Writing the difference as a per-observation adjustment leaves
  every crossprod reading the block as it did, and makes the scoring step's
  increment the Gauss-Newton one. Measured against modelterms7's own
  iteration, which shares no code with the fitting layer: the break-point
  agrees to 0.01 and is the same from every start tried, and `nl()` agrees
  with `nls()` to its own tolerance.

  The refresh is committed once per sweep and not once per objective
  evaluation, since the rescaling factor of a discontinuous term is a state
  of the iteration rather than a function of the point: a schedule advancing
  at the speed of a line search is not the one the construction was designed
  with, and one that never advanced would solve a permanently smoothed
  problem.

  The verdict on such a fit is the objective's and the term's own, through
  `modelterms7::term_converged()`, not the inner score's. Where the block is
  a working linearization with a frozen weight the score belongs to that
  working model: measured on `jump()`, the fit reaches the break-point and
  the jump size to three figures and the score stays at 0.176 forever.

* `statmod()` fits a structural term of the FILTER shape, which is `gas()`.
  The term contributes no columns; its level is added to the predictor of the
  equation it sits in, and its own parameters are estimated on the
  unconstrained scale of `term_links()` by a general optimizer with an exact
  gradient, alternating with the coefficients.

  What makes the answer right is the gradient. The level was driven by scores
  read at predictors the coefficients also enter, so the derivative of the
  objective in a coefficient is not the block times the score;
  `modelterms7::term_adjoint()` supplies the correction and the layer chains
  it with its own mixed second derivatives, which reaches the coefficients of
  the OTHER equations too. Against `numDeriv` the whole gradient agrees to
  2e-8, and the naive one is out by 9.7. On 2000 observations simulated from
  the model the fit returns (0.299, 0.440, 0.732) against a truth of
  (0.3, 0.4, 0.7), with the scale at 0.998 against 1.

  `fit@structural` carries each term's estimated parameters, on their own
  scale and on the unconstrained one.

* `regime()` is still rejected, and the message now says what the layer is
  missing -- the derivative of a likelihood mixed over latent states in the
  predictor it is handed -- rather than naming an internal generic.

* A formula carrying more than one structural term is rejected, whatever
  equations they sit in. Each rewrites the predictor or the likelihood the
  others would be read at, so the pair is one recursion written as two and
  neither term implements it.

* An accepted step that moves the objective by less than its own rounding
  ends the scoring run. The sufficient-decrease test cannot tell such a step
  from no step -- at a step of 1e-9 the decrease it asks for is of that order
  -- so the run spent its budget standing still, which is what a term whose
  gradient belongs to a working model produces at its own fixed point.

# statmodels7 0.13.0

* A term carrying more than one penalty is fitted, counted and reported per
  penalty. `statmod_penalized()` already enumerated them, so the fitting core
  was right; what still read one penalty per term were the three places that
  report a fit.

  - `statmod_edf()` hands `modelterms7::edf()` the hyperparameters of every
    penalty the term carries, keyed by the penalty names, which is the shape
    that function reads them in. It asked for `term_penalty()` before, so a
    term whose penalty covers part of its parameters was counted as
    unpenalized -- its coefficient count -- and a term carrying two had the
    hyperparameters of neither.
  - `term_block_kind()` reads the same enumeration, so a term penalized over
    part of itself is a penalized block, and a selection where any of its
    penalties has a kink, rather than a parametric one.
  - `summary()` shows one hyperparameter row per penalty. Where a term
    carries several the name carries the penalty as well (`delta.lambda`),
    two hyperparameters in one block not being the same number.

  The row of the degrees of freedom stays per TERM, which is the granularity
  a table of terms wants.

* `statmod_entry_key()` composes the key -- the term's name, with the entry's
  after `::` where a term carries several -- in one place. Two callers
  composing it apart would agree by accident.

# statmodels7 0.12.0

* `statmod_penalized()` enumerates the penalized units of a model -- one
  entry per penalty, whatever term it belongs to and whether or not that
  term has more than one -- and the fitting core reads it instead of
  looping over terms itself. Twelve places ran the same loop, over the
  distribution parameters and over each one's terms, and every one of them
  assumed a term carries at most one penalty. Three of the twelve are the
  fitting core and now read the enumeration: the block split, the penalty
  assembled at given coefficients, and the starting hyperparameters.

  The key is the term's name in the formula, with the entry's own name
  after `::` where a term carries several. A term with one penalty over the
  whole of itself keys exactly as it did, so nothing that reads a
  hyperparameter by term name changes, and two `ridge()` terms remain two
  terms with two hyperparameters -- which they already were, at every layer.

  What this opens is a penalty over parameters that are not coefficients of
  a design block, which is what `gas()` on panel data, `nl()` and `seg()`
  need for a population value with shrunk deviations per group.

* The sweep is finished: the hyperparameter index, the integrated basis of
  `ml()`, the exact gradient and Hessian of the marginal criterion, and the
  labels a summary reads all enumerate penalties rather than terms.
  `statmod_unit()` looks one up by parameter and key, which is the same
  answer where a term carries one penalty and the right answer where it
  carries several.

  What still reports per term rather than per penalty is `edf()`, which is
  the right granularity for a table of terms and would need
  `modelterms7::edf()` to accept more than one set of hyperparameters. A
  term with two penalties currently has its degrees of freedom reported as
  its coefficient count.

# statmodels7 0.11.0

* The path carries the size of the kink at the point just fitted onto the
  blocks, so the next point can screen against it. It travels there rather
  than through the argument list of every layer between the path and the
  descent: it is a property of the block, where its penalty was a moment
  ago, and the path rebuilds the blocks at each point anyway.

  **Measured, the rule still discards nothing**, and the reason is the grid
  rather than the data. On a geometric grid of ratio r the rule keeps a
  coordinate whose gradient exceeds `(2 - 1/r)` times the kink. At the
  default 25 points over `min_ratio` 1e-3 the ratio is 0.75 and the
  threshold two thirds of the kink, which almost every inactive coordinate
  clears; at `glmnet`'s 100 points over 1e-2 the ratio is 0.955 and the
  threshold 0.95 of the kink, just under it, so most inactive coordinates
  fall below. The rule lives on fine grids and ours is deliberately coarse,
  each point being a whole fit. Measured over a path: 5 per cent of the
  columns discarded and 0.93 to 0.96 times the speed.

* A fit that did not converge now reports the range each distribution
  parameter reached. The case it exists for: a lasso at a fixed
  hyperparameter with a free scale, on a design the model can interpolate.
  Fitting the coefficients shrinks the residuals, which shrinks the scale,
  which raises the working weights, which makes the penalty count for
  relatively less, which lets more coefficients in. At 200 observations and
  400 columns the scale reached 3.8e-15 and 380 of the 400 coefficients
  survived, where the same block fitted at a held scale kept the five that
  were real.

  Nothing diagnoses that. The parameters it reached are a fact and a scale
  at 1e-15 says the rest, where naming a cause would mean choosing a
  threshold for running away -- and the same fit at 100 columns converges
  to a scale of 0.77 that is nothing of the kind. What that fit is doing is
  the joint mode over coefficients and scale, which is degenerate wherever
  the likelihood can be driven up without bound; it is the reason the
  hyperparameters here are estimated by a marginal criterion or by
  prediction rather than jointly.

# statmodels7 0.10.1

* The coordinate descent kernel carries the two devices `glmnet` uses to
  go faster, and measured at a single hyperparameter **neither pays**.
  Both are path devices and the reason each one loses is visible:

  Covariance updating replaces an O(n) read of the gradient with an O(m)
  one, at the cost of building a column of X'WX the first time a
  coordinate moves off zero. That column costs O(nm), so it is worth it
  only while m is small: at 5000 observations with 200 columns and
  nothing screened away, taking that route cost 70 milliseconds against
  55 for the residual. `coord_covariance()` sets the threshold where the
  two cross rather than at a rule of thumb.

  The sequential strong rule discards a coordinate whose gradient at the
  previous point of a path is below `2*s_k - s_{k-1}`. With no previous
  point the reference is the kink that empties the block, and
  `2*s - max|g|` is negative at any value worth fitting at, so the rule
  discards nothing and costs one crossprod. It is therefore not attempted
  without a previous point, and passing one from the path is what is
  still owed.

  The rule assumes the gradient moves no faster than the threshold, which
  is not a theorem: it can discard a coordinate that belongs in the fit.
  What makes the answer exact is the check afterwards, reading the
  gradient over every column at the point reached and putting back
  whatever exceeds the kink. A test screens the block down to one column
  and to two and gets the same answer back; without the check, one column
  gives a different one.

# statmodels7 0.10.0

* A compiled coordinate descent fits a block whose penalty is separable,
  in `src/coord_descent.cpp`, the package's first compiled code. It reads
  the block's own columns and the running residual -- the model rather
  than the objective -- which is why it is not an optimizer and lives
  here. Measured against the proximal route on the same block, agreeing
  with it to 1e-10 on the coefficients and to 1e-9 on the objective:

  |  | n = 200, p = 20 | n = 1000, p = 100 | n = 5000, p = 200 |
  |---|---|---|---|
  | lasso | 26x | 135x | 27x |
  | elastic net | 40x | 9x | 43x |
  | SCAD | 21x | 154x | 24x |
  | MCP | 353x | 156x | 7x |

  The proximal route read the whole model at every step: one block fit at
  200 observations made 88 evaluations of the objective, 75 of the
  gradient and 83 of the operator, each over every parameter of the
  distribution, and closed in 36 iterations where a coordinate descent
  closes in six sweeps.

* The penalty arrives as a piecewise linear table from
  `penalties7::penalty_prox_spec()`, so the kernel names no family and a
  penalty that describes its operator gets the compiled route without an
  edit here. A penalty with no table -- a heavy-tailed prior, whose
  operator is a root -- keeps the proximal route.

* **Against the reference packages we are still slower**, and the lasso is
  the only comparison whose objectives match exactly, being homogeneous of
  degree one in lambda so that a rescaling of the loss is absorbed. There
  we agree with `glmnet` to 0, 1.4e-14 and 1.1e-13 at the three sizes and
  take 4.2, 10.6 and 12.6 times as long. For the elastic net, SCAD and MCP
  the shape parameters make the penalty non-homogeneous, so a lambda
  rescaling cannot match the two objectives and only the timing compares:
  0.22x to 0.02x against `glmnet` and `ncvreg`.

  Of the remaining distance, at n = 5000 and p = 200 the kernel is 21 ms
  of a 60 ms block fit and the R around it is the other 36; what
  `glmnet` has beyond that is covariance updates and strong rules, which
  this kernel does not.

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
  a location of 0 for a response centered at 5.84. The intercept of each
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
  one smoothing parameter, for several, with the scale modeled, and
  under `ml()`.

* `outer_optimizer` defaults to `lbfgs()` where the gradient exists and
  to `nelder_mead()` where it does not. Measured in evaluations of the
  criterion against the derivative-free search: 40 against 32 with one
  smoothing parameter, 40 against 135 with two, 41 against 269 with
  three.

# statmodels7 0.4.0

* The modeling layer. `statmod()` reads one formula carrying every
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
