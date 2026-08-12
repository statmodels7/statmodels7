# Changelog

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
