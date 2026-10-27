#' @include blocks.R
NULL

#' A Fitted Model
#'
#' @description
#' [statmod()] returns an object of this class. It holds three things: the
#' specification kept whole, so the model can be read back and reapplied to
#' new data; the
#' coefficients, hyperparameters and structural parameters that the fit
#' reached; and the record of the run that produced them.
#'
#' Reach the contents through the accessors rather than the properties where
#' one exists: [coef.StatmodFit()], [vcov.StatmodFit()],
#' [predict.StatmodFit()], [summary.StatmodFit()], [logLik.StatmodFit()].
#'
#' @param spec The [StatmodSpec()] that was fitted, with every refreshable
#'   term as the fit left it.
#' @param coefficients A named list, one numeric vector per distribution
#'   parameter.
#' @param structural A named list with one entry per structural term. Each
#'   entry is a list with the term's parameters on their own scale
#'   (`parameter`), on the unconstrained scale (`unconstrained`), the
#'   equation that carries the term (`equation`) and the names of the
#'   parameters that the term holds (`held`).
#' @param hyper The hyperparameters, keyed by distribution parameter and then
#'   by term.
#' @param loglik The weighted log-likelihood at the fitted coefficients, a
#'   single number. It is the conditional one.
#' @param objective The penalized objective there, not averaged.
#' @param edf A data frame of effective degrees of freedom with the columns
#'   `parameter`, `term`, `coefficients` and `edf`, one row per term of each
#'   equation. The count of a term is its share of
#'   \eqn{\mathrm{tr}[(H+S)^{-1}H]}, with \eqn{H} the information of the
#'   log-likelihood and \eqn{S} the second derivative of the penalty, and it
#'   is read at fit time on one information matrix, which `methods` records as
#'   `edf_information`. That matrix is the `information` of [aic()] or
#'   [bic()] where either prices the fit (the count is part of the
#'   criterion), otherwise the `information` of [reml()] or [ml()]
#'   ([distributions7::observed_info()] by default), and the information of
#'   the inner fit (`inner_information` of [statmod()], with [auto_info()]
#'   settled against the family) only where the model has no outer criterion.
#'   For a positive semidefinite \eqn{H} the total lies in \eqn{[0, p]}, with
#'   \eqn{p} the number of coefficients, and the count of a term lies between
#'   0 and its number of coefficients. An indefinite observed information can
#'   give a negative share to a direction of negative curvature, and then
#'   these bounds can fail. [summary.StatmodFit()] reports the count on the
#'   information of its `vcov_type` and leaves `fit@edf` unchanged.
#' @param fitted A named list of the fitted parameters, one vector per
#'   distribution parameter.
#' @param converged A single logical. Where the fit does not alternate between
#'   blocks, it is the verdict of the inner method. Where the model has a
#'   kinked block or a structural term, it is read at the point reached: the
#'   Newton decrement of the smooth part and the optimality decrement of every
#'   coefficient that a kinked penalty holds at zero are both within
#'   [mode_error_limit()]. A term that recomputes its own block, such as a
#'   break-point term, contributes its own verdict in place of that of the
#'   inner method.
#' @param aliased The labels, as [coef_labels()] names them, of the
#'   coefficients that the fitted point does not identify. A coefficient is
#'   listed when the pivoted decomposition of the inner fit drops it, when the
#'   penalized information at the fitted point is singular in its direction
#'   (this applies to every decomposition and to an \pkg{optimizers7} inner
#'   optimizer), or when its design column has vanished. The vector is empty
#'   where the design has full rank. The estimate, standard error and
#'   interval of a listed coefficient are reported as missing, and it is not
#'   counted in the degrees of freedom. The stored coefficient keeps the value
#'   that it had when the coordinate was left out of the increment, so every
#'   predictor reads the vector unchanged.
#' @param elapsed The elapsed time in seconds.
#' @param criterion The value of the criterion that chose the
#'   hyperparameters: [reml()], [ml()], [aic()], [bic()] or [cv()], including
#'   the criterion that swept the path of a kinked penalty. The value is `NA`
#'   when no criterion ran.
#' @param history A list of data frames recording the run: `outer`, `blocks`
#'   and `inner`. An entry is `NULL` where the run produced no such record.
#' @param methods A list of what was settled during the run, which can differ
#'   from what was supplied: `smooth`, the inner method after [iwls()] has
#'   settled its information (never `"auto"`); `outer`, the criterion after
#'   resolution, `NULL` when no marginal or prediction-error criterion ran;
#'   `search`, the optimizer that the outer search used; `pinned`, the
#'   coordinates that the criterion left out of its determinant; `kkt`, the
#'   edge violations that survived the restarts; `sparse`, the blocks of
#'   kinked penalties; `sparse_criterion`, the criterion that swept their
#'   paths; and `sparse_hyper`, the hyperparameters that it reached.
#' @param call The matched call.
#'
#' @return An object of class `StatmodFit` with one property per argument
#'   above.
#'
#' @seealso [statmod()], which builds one, [summary.StatmodFit()] for the
#'   printed report, [statmod_certificate()] for the verdict on convergence.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(30))
#' dd$y <- 1 + 2 * dd$x + rnorm(30, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' S7::S7_inherits(fit, StatmodFit)
#'
#' # The coefficients form one vector per parameter of the family.
#' fit@coefficients
#' fit@converged
#'
#' @name StatmodFit-class
#' @aliases StatmodFit
#' @keywords internal
#' @export
StatmodFit <- S7::new_class("StatmodFit",
  properties = list(
    spec = S7::class_any,
    coefficients = S7::class_list,
    structural = S7::class_list,
    hyper = S7::class_list,
    loglik = S7::class_numeric,
    objective = S7::class_numeric,
    edf = S7::class_any,
    fitted = S7::class_list,
    converged = S7::class_logical,
    aliased = S7::class_character,
    elapsed = S7::class_numeric,
    criterion = S7::class_numeric,
    history = S7::class_list,
    methods = S7::class_list,
    call = S7::class_any,
    # the bootstrap entries of add_bootstrap(), kept in an environment so that
    # a variance asked for later can store the replicates it computes; empty
    # after statmod(), and read only through statmod_bootstrap()
    bootstrap = S7::new_property(S7::class_environment,
                                 default = quote(new.env(parent = emptyenv())))
  )
)


#' Fit a Model
#'
#' @description
#' Fits a distributional regression. One formula carries an equation per
#' parameter of the response distribution, so a scale or a shape is modeled
#' as readily as a mean; the terms that those equations name are assembled into
#' a penalized likelihood and fitted, and any smoothing parameters or prior
#' scales are estimated from the data.
#'
#' Returns a [StatmodFit()], which [summary.StatmodFit()],
#' [predict.StatmodFit()], [coef.StatmodFit()] and the other accessors read.
#'
#' @details
#' # The formula
#'
#' The equations are separated by `|`, and the first carries the response:
#'
#' \preformatted{    y ~ x1 + ridge(R) + lasso(L)  |  sigma ~ z  |  nu ~ 1}
#'
#' A parameter with no equation of its own is given an intercept, so
#' `y ~ x` on a two-parameter family fits a constant scale.
#' [statmod_equations()] does the split. R's operator precedence nests the
#' equations on the left, so the right-hand side of the whole formula holds
#' only the last equation, and the split walks down the left side instead.
#'
#' The terms available are \pkg{modelterms7}'s: [modelterms7::s()],
#' [modelterms7::te()], [modelterms7::random()], [modelterms7::ridge()],
#' [modelterms7::lasso()], [modelterms7::seg()], [modelterms7::nl()],
#' [modelterms7::gas()] and the rest. Bare covariates collapse into one
#' parametric block. A formula can carry one structural term, which is
#' [modelterms7::gas()], [modelterms7::regime()] or the marginal break-point
#' term, or one [modelterms7::gas()] term in each of several equations,
#' which are filtered jointly as one recursion.
#'
#' A structural term whose level adds a constant to the predictor, the
#' \eqn{\omega} of [modelterms7::gas()] or the first level of
#' [modelterms7::regime()], is confounded with an intercept in the same
#' equation. The intercept is kept and the level held at zero, reported as
#' `omega (held)` for [modelterms7::gas()] and as `level1 (held)` for
#' [modelterms7::regime()]; writing `0 + gas(...)` estimates the level
#' instead. The two spellings are the same model with the same
#' log-likelihood, and the constants differ because they measure different
#' quantities: for a score-driven term the intercept is the stationary level
#' \eqn{\omega / (1 - \sum_j \xi_j)}, which is \eqn{\omega} divided by one
#' minus the sum of the persistences \eqn{\xi_j}. The page of
#' [modelterms7::gas()] gives the identity.
#'
#' # The fitting scheme
#'
#' The terms split in two by a property that each one already reports.
#'
#' A term whose penalty is twice differentiable in its coefficients, an
#' unpenalized block, a ridge, a spline or a random effect, is estimated in
#' one system by `inner_optimizer`. Their joint curvature exists, and using
#' it closes a fit in a handful of iterations.
#'
#' A term whose penalty has a kink, a lasso, an elastic net, a SCAD or an
#' MCP, is estimated as a block of its own while the other blocks are held
#' fixed.
#' The block is fitted by a compiled coordinate descent where the penalty
#' describes its proximal operator as a table at the step taken, and
#' otherwise by an accelerated proximal gradient. Where the equation has a
#' free unpenalized intercept, the descent centers the block and moves that
#' intercept jointly with it.
#'
#' A structural term of the filter shape, [modelterms7::gas()], is fitted in
#' one system with the coefficients. A structural term of the likelihood
#' shape, [modelterms7::regime()] or the marginal break-point term, is
#' alternated with them.
#'
#' The fit alternates between the blocks until the relative change of the
#' objective falls below the tolerance. Where a kinked block or a structural
#' term is present, the point reached must in addition satisfy the optimality
#' conditions of the smooth part and of every coefficient held at a kink, and
#' `converged` reports that reading.
#'
#' # The objective is not averaged
#'
#' The objective is minus the weighted log-likelihood plus the penalties at
#' full size. A penalty is a negative log-prior, and a posterior adds a
#' log-likelihood and a log-prior at full size; averaging one of them would
#' make a hyperparameter mean something that depends on \eqn{n}. The stopping
#' rule is scaled instead, so a threshold means the same at
#' \eqn{n = 10} and at \eqn{n = 10^7}.
#'
#' The weights multiply the log-likelihood, so they act as frequency
#' weights: a constant weight \eqn{c} on every row is \eqn{c} copies of the
#' data. The restricted (REML) criterion then counts \eqn{n c} observations
#' against the \eqn{p} coefficients, so the REML estimate of a dispersion
#' changes under a constant rescaling of the weights, and the ML estimate
#' does not.
#'
#' # The budget and the stopping rule belong to the method
#'
#' There is no `maxit` and no `tol` here. Both are set on `inner_optimizer`,
#' which is [`iwls(maxit =, tol =)`][iwls] or an \pkg{optimizers7} optimizer
#' with its own `maxit` and `criterion`. The alternation reads them from
#' there as its budget in passes and its relative tolerance on the
#' objective, so that one setting governs both the inner fit and the
#' alternation.
#'
#' # Every hyperparameter is estimated unless its term holds it
#'
#' The term that carries the penalty states which hyperparameters are held:
#' `lasso(x, lambda = 3)`, `ridge(x, lambda = 0.5)`,
#' `s(x, bspline_smooth(), hyper = c(lambda = 2))`,
#' `enet(x, alpha = 0.5)`. Every hyperparameter left `NULL`, which is each
#' term's default, is estimated from the data. A hyperparameter is held in
#' the term and nowhere else, so that the model has one statement of each
#' held value.
#'
#' The smooth hyperparameters go to `outer_criterion`, [reml()] by default.
#' `outer_optimizer` searches over them, and the coefficients are refitted at
#' each point.
#'
#' A kinked penalty has its own argument. `sparse_criterion`, [bic()] by
#' default, sweeps its hyperparameters along a path of their own values,
#' from the kink that empties the block down to a fraction of it. The
#' penalized mode is only piecewise smooth in such a hyperparameter, turning
#' a corner whenever a coefficient joins the active set or leaves it, so a
#' gradient search on a criterion read there would meet those corners.
#'
#' The top of that path is the kink that empties the block, found at the
#' coefficients in hand and not at a refitted null. It therefore depends on
#' the data and on the rest of the model, since the fits of the other terms
#' enter it.
#'
#' Where a model carries both kinds, the path is outside and the marginal
#' criterion is estimated inside each of its points, so a smoothing parameter
#' can come from REML and a lasso's \eqn{\lambda} from BIC in one fit.
#'
#' With the default `marginal = NULL`, [reml()] and [ml()] also estimate on
#' the criterion the unpenalized coefficients of every distribution parameter
#' other than the first. `outer_criterion` therefore runs when the model
#' carries a smooth hyperparameter or such a coefficient, and on an ordinary
#' `y ~ x` with a two-parameter family it estimates the dispersion by REML.
#' Beside a kinked penalty, a structural term or a block that moves with its
#' coefficients, the default names no coefficient, so a model whose only
#' penalties are kinked does not run the criterion. The default criterion is
#' also not run when the family lacks the two derivatives that it needs. With
#' `outer_criterion = NULL`, every smooth hyperparameter stays where its term
#' left it and the coefficients of all parameters are fitted at the joint
#' mode, which gives the maximum-likelihood estimate of a dispersion on an
#' ordinary model.
#'
#' # Verbosity
#'
#' Three levels name the loops: `1` shows the outer search and the
#' alternation, `2` adds the inner method's own iterations, and `3` adds the
#' optimizers' traces. A named logical vector sets the four switches
#' `outer`, `blocks`, `inner` and `optimizer` independently, and a switch
#' that is left out is off. For example `verbose = c(outer = TRUE, blocks =
#' FALSE)` shows the hyperparameters moving while a chatty inner optimizer
#' stays silent.
#'
#' @param formula The model formula.
#' @param distrib A \pkg{distributions7} distribution object.
#' @param data A data frame. Missing values in the response or in a variable
#'   that an equation uses are rejected, with the variable and the rows.
#' @param weights Optional prior weights, a numeric vector of length \eqn{n}
#'   or the name of a column of `data`. The weights must be non-negative,
#'   finite and complete, and are taken as given and not normalized. They
#'   multiply the log-likelihood, so they act as frequency weights: a constant
#'   weight \eqn{c} on every row is \eqn{c} copies of the data, and the
#'   restricted (REML) criterion counts \eqn{n c} observations against the
#'   coefficients. The REML estimate of a dispersion therefore changes under a
#'   constant rescaling of the weights, and the ML estimate does not.
#' @param offsets Optional named list of offsets, one entry per distribution
#'   parameter, each of length 1 or \eqn{n}. An offset written in the formula
#'   is added to the entry for the same parameter.
#' @param inner_optimizer How the smooth block is fitted: [iwls()] (the
#'   default) or an \pkg{optimizers7} optimizer. Beside a structural term of
#'   the filter shape, the coefficients and the term's parameters are fitted in
#'   one system, for which [iwls()] stands for the system's own default
#'   optimizer and an \pkg{optimizers7} optimizer is used as given.
#' @param inner_information The matrix of second derivatives that the inner fit
#'   uses, an object that inherits from [distributions7::Information]:
#'   [distributions7::observed_info()] (the observed Hessian),
#'   [distributions7::expected_info()] (the Fisher information, whose `method`
#'   is `NULL` for the closed form of the family or [distributions7::expect_opg()],
#'   or one of [distributions7::expect_opg()], [distributions7::expect_integrate()]
#'   and [distributions7::expect_mc()]; see
#'   [distributions7::expected_derivative_methods]) or [auto_info()] (the
#'   default), which is settled against the family. With [iwls()] the matrix
#'   is inverted at each step (Newton for the observed information, Fisher
#'   scoring for the expected one); with an \pkg{optimizers7} optimizer it is
#'   the Hessian that the optimizer receives, so that
#'   `inner_optimizer = optimizers7::newton()` with
#'   `inner_information = distributions7::expected_info()` takes the expected
#'   information. A named matrix is the only one that the inner fit uses.
#'   [auto_info()] chooses the expected information where the family writes it
#'   exactly and at an acceptable cost, and the observed information
#'   otherwise; with [iwls()] it also lets the run change matrix and lets the
#'   expected information take a step that the observed one cannot (see
#'   [information_objects]). The fitted model records the settled object, and
#'   every later refit of the model (the folds of a cross-validation, the
#'   restricted fits of a test, the resampled fits of a bootstrap) reads the
#'   same one. Beside a score-driven filter ([modelterms7::gas()]) the
#'   coefficients and the filter's own parameters are fitted in one system,
#'   whose Hessian is the observed information, or the expected information
#'   conditional on the past where [distributions7::expected_info()] is named
#'   here (Fisher scoring over that system); [auto_info()] leaves that system
#'   on the observed information. A model with a term of the likelihood shape
#'   ([modelterms7::regime()], a marginal break-point term) has no expected
#'   information, and [distributions7::expected_info()] is an error there. The
#'   outer criterion has its own
#'   `information` (see [reml()]) and the standard errors have theirs (see
#'   [vcov_method]); the three choices are independent.
#' @param outer_criterion How the smooth hyperparameters are estimated:
#'   [reml()] (the default), [ml()], [aic()] or [bic()], or `NULL` to hold
#'   them where their terms left them. [cv()] cross-validates the
#'   hyperparameter of a kinked penalty only: it is passed as
#'   `sparse_criterion`, and as `outer_criterion` it is an error. A
#'   prediction-error criterion beside a kinked penalty is rejected unless
#'   `sparse_criterion` is `NULL`.
#' @param sparse_criterion How the hyperparameter of a kinked penalty (lasso,
#'   elastic net, scad, mcp) is chosen: [bic()] (the default), [aic()],
#'   [cv()], or `NULL`. With `NULL` the hyperparameter is held at its value,
#'   unless the model also carries a smooth penalty and `outer_criterion` is a
#'   prediction-error criterion, which then sweeps the path.
#'   A marginal criterion is rejected here, because it is read at a mode that
#'   sits on the kink.
#' @param outer_optimizer The optimizer that searches over the smooth
#'   hyperparameters, or `NULL` to choose one from the derivatives of the
#'   criterion that exist: [optimizers7::newton()] where the exact Hessian is
#'   usable, [optimizers7::lbfgs()] with the exact gradient,
#'   [optimizers7::brent()] for one coordinate without a gradient, and
#'   [optimizers7::nelder_mead()] otherwise.
#' @param start Where the fit begins: a named list of coefficients, a
#'   [start_strategy()] such as [start_search()] or [start_from()], which
#'   takes the estimates of a model already fitted, or
#'   `NULL` for [start_intercepts()]. A strategy runs once, before the
#'   alternation between the coefficients and the hyperparameters begins. A
#'   global search therefore belongs here and not in `inner_optimizer`, where
#'   it would run again at every hyperparameter value that the outer search
#'   tries.
#' @param linpar_control How the unpenalized parametric block is built, as
#'   [linpar_options()] returns it: the storage and the contrasts
#'   for its factors. It reaches the implicit term, the one into which the bare
#'   covariates collapse and which a caller never writes; a
#'   [modelterms7::linpar()] written out takes them directly.
#'   The argument and the function are named differently on purpose, as
#'   [stats::glm()]'s `control` and
#'   [stats::glm.control()] are.
#' @param verbose A level from 0 to 3, or a named logical vector with any of
#'   `outer`, `blocks`, `inner` and `optimizer`.
#' @param threads How many threads the fit may use, as
#'   [numericals7::n_threads()] constructs it. The default,
#'   `n_threads(1)`, is sequential and takes the sequential code path. A
#'   larger count is passed to the compiled per-observation kernels of the
#'   family and to the dense assembly products; below a kernel's internal
#'   size threshold the computation stays sequential whatever the count,
#'   and a sparse design keeps its \pkg{Matrix} route. The object's `workers`
#'   fans the units that are independent by construction out over separate R
#'   processes: the folds of a cross-validation, and the combinations of a
#'   kinked path's product grid, each of which restarts its warm chain from
#'   the sweep's own starting coefficients. Each unit fits sequentially, so
#'   the two levels never nest. The points within one chain stay sequential,
#'   because a point fitted from a cold start costs more than one
#'   warm-started from the previous point, and splitting a chain would slow
#'   the single-process default or make the result depend on the count. For
#'   any two counts above one the result is identical bit for bit: every
#'   parallel region decomposes its work over the elements of its output and
#'   never splits a reduction, and the folds are seeded per fold and
#'   collected in fold order. Between the default and a count above one the
#'   route can change where the work is large (the BLAS crossproduct against
#'   the compiled one, and the pivoted `qr()` against the compiled
#'   triangular factor), and the results then agree to rounding.
#' @param ... Not used. A named argument signals an error. A hyperparameter
#'   is held in the term that carries the penalty.
#'
#' @return A [StatmodFit()] object.
#'
#' @seealso [summary.StatmodFit()] and [predict.StatmodFit()] for what to do
#'   with the result, [statmod_spec()] to build the model without fitting it,
#'   [iwls()] for the inner method, [reml()] and [bic()] for the criteria,
#'   [rstatmod()] to simulate from a model.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200, -2, 2))
#' dd$y <- 1 + 2 * dd$x + rnorm(200, sd = exp(0.3 * dd$x))
#'
#' # An ordinary regression: the scale receives an intercept that the
#' # formula does not specify.
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' coef(fit)
#'
#' # The scale is modeled too, which is the purpose of the framework. The data
#' # were drawn with log sigma = 0.3 x, and the interval covers it.
#' both <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' coef(both)$sigma
#' confint(both)["sigma:x", c("estimate", "lower", "upper")]
#'
#' # It is the better model on this data.
#' c(one = AIC(fit), both = AIC(both))
#'
#' @export
statmod <- function(formula, distrib, data, weights = NULL, offsets = NULL,
                    inner_optimizer = iwls(), inner_information = auto_info(),
                    outer_criterion = reml(),
                    sparse_criterion = bic(), outer_optimizer = NULL,
                    # the ARGUMENT and the FUNCTION are deliberately named
                    # differently, as `glm(control = glm.control())` is: with
                    # one name for both, the default resolves to the
                    # argument's own promise and R reports "promise already
                    # under evaluation"
                    start = NULL, linpar_control = linpar_options(),
                    verbose = 0, threads = numericals7::n_threads(), ...) {
  t0 <- proc.time()[["elapsed"]]
  cl <- match.call()
  # THE PRIOR WEIGHTS ARE LOOKED UP IN `data` FIRST, then where statmod() was
  # called, as glm() and lm() look them up: `weights = wt` naming a column
  # is the spelling a reader coming from those functions writes, and it used
  # to fail with "object 'wt' not found".
  weights <- eval(substitute(weights),
                  if (is.list(data)) data else NULL, parent.frame())
  # The count is validated once and travels DOWN on the specification; the
  # process-level RcppParallel setting is sized here and restored when this
  # frame exits, so a fit never leaves it moved for the code that runs
  # after it. At the default n_threads(1) neither call touches anything and
  # every kernel takes exactly the sequential path.
  n_thr <- numericals7::thread_count(threads)
  numericals7::local_threads(threads)
  asked <- !missing(outer_criterion)
  asked_sparse <- !missing(sparse_criterion)
  vb <- verbosity(verbose)
  budget <- method_budget(inner_optimizer)
  maxit <- budget$maxit
  tol <- budget$tol

  # WHICH hyperparameters are held is said by the terms, and by nothing
  # else. An argument here saying the same thing would be read by nobody
  # whenever the two disagreed.
  if (...length()) {
    nm <- ...names()
    stop(sprintf("unused argument (%s)", nm[[1L]]), call. = FALSE)
  }
  if (!is.list(linpar_control)) {
    stop("'linpar_control' must be the list that linpar_options() returns.",
         call. = FALSE)
  }
  spec <- statmod_spec(formula, distrib, data, weights, offsets,
                       linpar = linpar_control)
  spec@threads <- n_thr
  spec@workers <- numericals7::worker_count(threads)
  design <- statmod_design(spec)
  # a marginal criterion left at auto_info() is settled against the family
  # and the structural terms here, as an iwls() left at "auto" is below
  if (S7::S7_inherits(outer_criterion, OuterMethod)) {
    outer_criterion <- outer_resolve(outer_criterion, spec@distrib, design)
  }
  hyper <- statmod_hyper_start(spec, design)
  if (!is.null(outer_criterion)) {
    hyper <- statmod_hyper_scale_start(spec, design, hyper)
  }
  blocks <- statmod_blocks(spec, design)

  # the inner information is settled against the family once, here, and
  # travels with the method, so the objective, every refit and the fit's own
  # record read the same curvature
  check_information(inner_information, arg = "inner_information",
                    fn = "statmod")
  if (info_expected(inner_information)) {
    structural_expected_check(design, "'inner_information'")
  }
  inner_optimizer <- inner_bind(inner_optimizer, spec@distrib,
                                inner_information)
  cfg <- inner_settings(inner_optimizer)
  expected <- cfg$expected
  emethod <- cfg$emethod
  obj <- statmod_objective(spec, hyper, design, expected, emethod)
  beta <- statmod_start(spec, design, obj, start)

  if (vb$blocks || vb$outer) {
    vb_rule(sprintf("%s, %d observations", spec@distrib@distrib_name,
                    spec@n_obs), char = "=")
    vb_say("inner    %s", vb_name(inner_optimizer, "iwls"), indent = 5L)
    vb_say("outer    %s", if (is.null(outer_criterion)) "held" else
      paste0(toupper(outer_criterion@kind), " (",
             info_label(outer_criterion@hessian), " information)"), indent = 5L)
    if (!is.null(start)) {
      vb_say("start    %s",
             if (S7::S7_inherits(start, start_strategy_class())) start@label
             else "supplied values", indent = 5L)
    }
  }

  # A hyperparameter left where it started is a placeholder and not a choice,
  # so the criterion ESTIMATES by default: a model carrying a smooth, a ridge
  # or a random effect is one whose author wants those chosen from the data,
  # and mgcv and lme4 both estimate by default for the same reason.
  #
  # It applies to the SMOOTH penalties and to nothing else, and it comes into
  # play if and only if there is one. A kinked penalty -- lasso, scad, mcp --
  # is not estimated by a marginal criterion at all but by a path over its own
  # values, so it is absent from `outer_hyper_index()` and a model carrying
  # only kinked penalties, like a model carrying none, leaves the criterion
  # with nothing to do and it is not run. That is a property of the model
  # rather than of how the argument was written, so typing the default
  # changes nothing: there is one rule and no hidden second one.
  if (!is.null(outer_criterion) && !S7::S7_inherits(outer_criterion,
                                                   OuterMethod)) {
    stop("'outer_criterion' must be reml(), ml(), aic(), bic() or NULL.",
         call. = FALSE)
  }
  if (!is.null(sparse_criterion) && !S7::S7_inherits(sparse_criterion,
                                                    OuterMethod)) {
    stop("'sparse_criterion' must be bic(), aic(), cv() or NULL.",
         call. = FALSE)
  }
  # a kinked penalty is chosen by a path over its own values, which only a
  # prediction-error criterion scores: a marginal one is a Laplace expansion
  # at the mode and the mode sits AT the kink for every coefficient set to
  # zero, where the second derivative it asks for does not exist
  if (!is.null(sparse_criterion) && !outer_minimize(sparse_criterion)) {
    stop(paste0("'sparse_criterion' must be bic(), aic() or cv(): a kinked",
                " penalty is chosen by a path over its own values, and",
                " reml() and ml() are read at a mode that sits on the",
                " kink, where the curvature that they need does not exist."),
         call. = FALSE)
  }
  # cv() refits on folds, which only the path over a kinked penalty's own
  # values does, so it is a sparse criterion and never an outer one
  if (!is.null(outer_criterion) && identical(outer_criterion@kind, "cv")) {
    if (!nrow(path_rows(spec, blocks, hyper))) stop_cv_nothing()
    stop(paste0("'outer_criterion' is cv(), which chooses only the",
                " hyperparameter of a penalty\n  with a kink. Pass it as",
                " sparse_criterion = cv(); the smoothing parameters\n  are",
                " estimated by reml(), ml(), aic() or bic()."), call. = FALSE)
  }
  if (!is.null(sparse_criterion) &&
      !nrow(path_rows(spec, blocks, hyper))) {
    sparse_criterion <- NULL
  }
  structural_outer_check(sparse_criterion, design)
  # and it steps aside for a caller who set the hyperparameters by hand:
  # `hyper` says HELD AT THESE VALUES, so estimating them away would answer a
  # question nobody asked. Asking for a criterion explicitly overrides that,
  # in which case `hyper` is where the search starts.
  # "nothing to do" is asked OF THE CRITERION, since the two kinds reach
  # different penalties. A marginal one -- reml(), ml() -- reaches the smooth
  # penalties, which is what `outer_hyper_index()` enumerates. A
  # prediction-error one -- aic(), bic(), cv() -- walks a path over the
  # hyperparameter's own values and so reaches a KINKED penalty too, which is
  # how a lasso's lambda is chosen; for those, having nothing to do means
  # carrying no penalized term at all.
  # and it has something to do where a distribution parameter's unpenalized
  # coefficients are estimated on it, which reml() and ml() do by default for
  # every parameter but the position: a model with no penalty at all still
  # runs the criterion when it carries a dispersion or a shape. Asking here
  # also checks an explicit `marginal` against the family before any fit.
  # A FAMILY WITHOUT THE TWO DERIVATIVES the criterion reads keeps its
  # dispersion at the joint mode when the criterion would be there for that
  # alone: the default reml() carries no smooth hyperparameter to estimate,
  # and refusing `y ~ x` with a Laplace response for the sake of the
  # dispersion's intercept would refuse a model that needs no criterion at
  # all. A criterion the caller named is still refused below, by name.
  if (!is.null(outer_criterion) && !asked &&
      is.null(outer_criterion@marginal) &&
      !nrow(outer_hyper_index(spec, blocks)) &&
      !is.null(order_shortfall(spec@distrib, 2L, "the REML criterion"))) {
    outer_criterion <- NULL
  }
  # the default criterion, kept: a working break-point block refuses to put a
  # dispersion on it, and once statmod_settle_breakpoints() has held the
  # block the refusal no longer applies
  crit_default <- if (!asked) outer_criterion else NULL
  if (!is.null(outer_criterion)) {
    nothing <- !nrow(outer_hyper_index(spec, blocks)) &&
      !length(marginal_coords(spec, design, outer_criterion)$where)
    if (nothing) outer_criterion <- NULL
  }
  structural_outer_check(outer_criterion, design)

  # A criterion is a Laplace expansion at the penalized mode and reads the
  # curvature there, so it needs TWO derivatives of the log-density in every
  # parameter the model carries an equation for.  Where the family has fewer,
  # the expansion the criterion is built on does not exist, and no accuracy in
  # the search repairs it.  It is refused HERE, naming the family, the
  # parameter, the order and the remedy, rather than at the first evaluation,
  # where the same failure surfaces as an inner fit that did not converge --
  # which is true and is three layers from the cause.  See R/order.R.
  assert_criterion_order(spec@distrib, outer_criterion)
  assert_criterion_order(spec@distrib, sparse_criterion)
  # and a prediction-error criterion cannot select a penalty over a structural
  # term's own parameters, reading its degrees of freedom over the coefficients
  # alone, where that penalty covers nothing. Refused here, by name, before a
  # search that would stop at its first evaluation reporting success. An outer
  # criterion reaches the kinked penalties too when no sparse one is given, the
  # path then being swept by it. See R/order.R.
  assert_criterion_reach(spec, design, outer_criterion,
                         if (is.null(sparse_criterion)) "all" else "smooth")
  assert_criterion_reach(spec, design, sparse_criterion, "kinked")

  # A prediction-error criterion for the SMOOTH penalties cannot be nested
  # inside a path over the kinked ones: it scores the same quantity at two
  # levels, and measured, every point of the path came back NA -- the
  # penalized information at a smoothing parameter chosen by one criterion
  # and read by the other is not positive definite -- so the path had nothing
  # to choose between and the hyperparameter kept its starting value while
  # the fit reported success. It is rejected rather than left to do that.
  if (!is.null(sparse_criterion) && !is.null(outer_criterion) &&
      outer_minimize(outer_criterion)) {
    stop(paste0("'outer_criterion' is ", outer_criterion@kind, "() and the",
                " model also carries a penalty with a kink, whose",
                " hyperparameter is chosen by a path. A prediction-error",
                " criterion cannot be nested inside that path: it would",
                " score the same quantity at both levels. Use reml() or",
                " ml() for the smooth penalties and 'sparse_criterion'",
                " for the kinked ones, or hold the hyperparameters of one",
                " of the two in their terms, as in lasso(x, lambda = 3)."),
         call. = FALSE)
  }

  # ONE FIT, FROM ONE START. Written as a function so a structural term that
  # asks for several starts is fitted from each exactly as from the first:
  # the same branch, the same refinements, the same restarts.
  hyper0 <- hyper
  fit_once <- function(beta_start) {
    hyper <- hyper0
    spec_h <- spec
    # A criterion that cannot be read at its start because a break-point's
    # working iteration has not converged there is read again after the
    # positions are settled. Measured on segmented::globTempAnom with four
    # changes of slope, the iteration stopped on two observed years and the
    # REML criterion stopped the fit with an error; settled, the positions
    # reach the profile's minimum and the criterion is read there.
    deferred <- FALSE
    if (!is.null(outer_criterion) && length(attr(design, "refresh"))) {
      res <- tryCatch(
        statmod_select(spec, design, blocks, hyper, inner_optimizer,
                       outer_criterion, outer_optimizer, beta_start,
                       emethod, maxit, tol, vb, data,
                       sparse_criterion),
        error = function(e) e)
      if (inherits(res, "error")) {
        if (!grepl("unavailable at the starting hyperparameters",
                   conditionMessage(res), fixed = TRUE)) {
          stop(res)
        }
        deferred <- conditionMessage(res)
      }
    }
    if (!isFALSE(deferred) ||
        (is.null(outer_criterion) && is.null(sparse_criterion))) {
      res <- statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                               beta_start, expected, emethod, maxit, tol, vb)
      crit <- NA_real_
    } else {
      if (!length(attr(design, "refresh")) || is.null(outer_criterion)) {
        res <- statmod_select(spec, design, blocks, hyper, inner_optimizer,
                              outer_criterion, outer_optimizer, beta_start,
                              emethod, maxit, tol, vb, data,
                              sparse_criterion)
      }
      hyper <- res$hyper
      crit <- res$criterion
      # the coefficients the criterion estimated stay where it put them in
      # every refit below; the specification the fit returns holds nothing
      spec_h <- if (length(res$held_coef)) {
        S7::set_props(spec, held_coef = res$held_coef)
      } else spec
      # Everything inside the selection held the frozen break-point blocks at
      # their committed positions -- a break-point moving between criterion
      # evaluations makes the criterion path-dependent, and the phase's own
      # flags read as unavailable points to the search -- so the positions are
      # refined here, once, at the chosen hyperparameters.
      frozen <- any(vapply(attr(design, "refresh"),
                           function(r) isTRUE(r$frozen), logical(1)))
      if (frozen) {
        ro <- tryCatch(
          statmod_alternate(spec_h, design, blocks, hyper, inner_optimizer,
                            res$par, expected, emethod, maxit, tol, vb),
          error = function(e) NULL)
        if (!is.null(ro) && is.finite(ro$value)) {
          res[c("par", "value", "converged", "obj", "aliased",
                "hist_blocks", "hist_inner")] <-
            ro[c("par", "value", "converged", "obj", "aliased",
                 "hist_blocks", "hist_inner")]
        }
      }
    }

    # The bootstrap restarting the break-point terms declare, run ONCE at the
    # top level -- inside an outer search it would multiply by the number of
    # criterion evaluations -- and at the hyperparameters the fit ended at.
    nb <- seg_boot_total(spec)
    if (nb > 0L && length(attr(design, "refresh"))) {
      before <- res$par
      res <- statmod_boot_restart(spec_h, design, blocks, hyper, inner_optimizer,
                                  res, expected, emethod, maxit, tol, vb, nb)
      # A RESTART THAT MOVES THE FIT MOVES THE COEFFICIENTS THE CRITERION
      # ESTIMATES TOO. The search above read the criterion in the basin its
      # warm starts stayed in, and the restart may leave that basin; a
      # coefficient held at the old optimum is then the optimum of another
      # model. Measured on a smoothed jump, sigma stayed at 0.661 after the
      # restart had moved the break-point from 6.46 to 6.09, where the criterion
      # puts it at 0.561. The search is run once more from the restarted point.
      if (length(res$held_coef) && !identical(before, res$par)) {
        r2 <- tryCatch(
          statmod_select(spec, design, blocks, hyper, inner_optimizer,
                         outer_criterion, outer_optimizer, res$par, emethod,
                         maxit, tol, vb, data,
                         sparse_criterion),
          error = function(e) NULL)
        if (!is.null(r2)) {
          res <- r2
          hyper <- res$hyper
          crit <- res$criterion
        }
      }
    }
    # A SHARP BREAK-POINT IS SETTLED ON THE PROFILE AND HELD THERE. The
    # working construction stops at a fixed point of its own iteration, whose
    # coefficients are those of the working model and whose position may sit
    # one interval away from the profile's minimum; see
    # statmod_settle_breakpoints().
    for (round in seq_len(3L)) {
      # The settle and the criterion are ALTERNATED until the positions stop
      # moving, at most three rounds: the hyperparameters move the working
      # response the polish reads, and the polish moves the criterion. Measured
      # on MASS::mcycle with a jseg in sigma's equation, the first round left the
      # break-point at the confinement limit, where the smoothing parameter had
      # been chosen, and the second moved it to the profile's minimum.
      sg0 <- statmod_settle_seg(spec_h, design, blocks, hyper, inner_optimizer,
                                res, expected, emethod, maxit, tol, vb)
      res <- sg0$res
      st0 <- statmod_settle_breakpoints(spec_h, design, blocks, hyper,
                                        inner_optimizer, res, expected, emethod,
                                        maxit, tol, vb)
      res <- st0$res
      fresh <- !isFALSE(deferred)
      if (fresh) {
        # the criterion that could not be read at the start, read now
        r2 <- tryCatch(
          statmod_select(spec, design, blocks, hyper, inner_optimizer,
                         outer_criterion, outer_optimizer, res$par, emethod,
                         maxit, tol, vb, data,
                         sparse_criterion),
          error = function(e) NULL)
        if (is.null(r2)) stop(deferred, call. = FALSE)
        deferred <- FALSE
        res <- r2
        hyper <- res$hyper
        crit <- res$criterion
        spec_h <- if (length(res$held_coef)) {
          S7::set_props(spec, held_coef = res$held_coef)
        } else spec
        sg0$moved <- TRUE
      }
      if (round > 1L && !st0$moved && !sg0$moved) break
      # a held term is a Jacobian block from here on: its contribution is X beta
      # exactly, so it is fitted like any other block and the rules that apply
      # to a working linearization -- the working phase, and the refusal to put
      # a dispersion's coefficients on the criterion (marginal_coords()) -- no
      # longer apply to it
      if (st0$settled) {
        d2 <- design
        rf2 <- attr(d2, "refresh")
        stt <- attr(d2, "state")
        for (i in seq_along(rf2)) {
          tm <- stt$terms[[rf2[[i]]$param]][[rf2[[i]]$term]]
          if (isTRUE(tryCatch(tm@blueprint$held, error = function(e) FALSE))) {
            rf2[[i]]$frozen <- FALSE
          }
        }
        attr(d2, "refresh") <- rf2
        design <<- d2
      }
      # the criterion is re-read on the held block whether or not a position
      # moved: the search ran on the working block, whose determinant carries
      # the auxiliary column, so its hyperparameters are another criterion's
      if (st0$settled && is.null(outer_criterion) && is.null(sparse_criterion) &&
          !is.null(crit_default) &&
          length(marginal_coords(spec, design, crit_default)$where)) {
        outer_criterion <<- crit_default
      }
      if (!fresh && (st0$settled || sg0$moved) &&
          (!is.null(outer_criterion) || !is.null(sparse_criterion))) {
        r2 <- tryCatch(
          statmod_select(spec, design, blocks, hyper, inner_optimizer,
                         outer_criterion, outer_optimizer, res$par, emethod,
                         maxit, tol, vb, data,
                         sparse_criterion),
          error = function(e) NULL)
        if (!is.null(r2)) {
          res <- r2
          hyper <- res$hyper
          crit <- res$criterion
          spec_h <- if (length(res$held_coef)) {
            S7::set_props(spec, held_coef = res$held_coef)
          } else spec
        }
      }
      # with no criterion nothing changed the working response, so a second
      # round would read the same polish
      if ((!st0$settled && !sg0$moved) ||
          (is.null(outer_criterion) && is.null(sparse_criterion))) break
    }
    list(res = res, hyper = hyper, crit = crit, spec_h = spec_h)
  }
  # SCAD and MCP are scaled by the curvature of the likelihood in each of
  # their coordinates, read once at the starting coefficients, so that the
  # knee sits at a lambda / c_j. The blocks are rebuilt from the scaled
  # penalties: left as they were, the coordinate descent would minimize the
  # unscaled penalty while the objective reads the scaled one.
  sc <- statmod_curv(spec, design, beta, expected, emethod)
  if (!is.null(sc)) {
    spec <- sc
    blocks <- statmod_blocks(spec, design)
  }
  # A SECOND START FOR A SCORE LOADING. Where the fit from the term's own
  # start does not converge, it is fitted once more with each loading at
  # 0.1/I (structural_start_fixups()), and the better of the two is kept by
  # the rule the multistarts below use. A beta share in Seatbelts stops
  # unconverged at -3326.66 from the first and converges at 450.48 from the
  # second; the second fit is paid only on a failure.
  # ⚠️ A FAILURE IS READ ON THE MODE AND NOT ON THE FLAG ALONE: where the
  # filter explodes the likelihood is flat to rounding, and the optimizer's
  # flag at that point is noise. On ubuntu's R 4.5.3 the Seatbelts fit from
  # the first start stopped at -3299.1 with the flag TRUE and a Newton
  # decrement of 3e51, so the second start was never tried; on Windows the
  # same point carried the flag FALSE. filter_decrement() reads the decrement
  # over the coefficients and the filter's own parameters.
  sst0 <- statmod_structural_state(design)
  z_start <- if (!is.null(sst0)) sst0$zeta else NULL
  alt <- if (!is.null(sst0)) sst0$info_start else NULL
  if (length(alt)) {
    f1 <- tryCatch(fit_once(beta), error = function(e) e)
    if (inherits(f1, "error") || !isTRUE(f1$res$converged) ||
        !isTRUE(filter_decrement(f1$spec_h, design, f1$res, f1$hyper) <=
                  mode_error_limit())) {
      zinfo <- sst0$zeta
      for (tn in names(alt)) sst0$zeta[[tn]] <- alt[[tn]]
      f2 <- tryCatch(fit_once(beta), error = function(e) NULL)
      if (inherits(f1, "error")) {
        if (is.null(f2)) {
          sst0$zeta <- zinfo
          stop(f1)
        }
        f1 <- f2
      } else if (!is.null(f2) &&
                 (edge_restart_better(f2$res, f1$res, f1$crit,
                                      outer_criterion, sparse_criterion) ||
                  # a converged fit is kept unless the other is strictly
                  # better: on a tie the unconverged one is the same point
                  # with a flag that says less about it
                  (isTRUE(f2$res$converged) &&
                   (!edge_restart_better(f1$res, f2$res, f2$crit,
                                         outer_criterion, sparse_criterion) ||
                    retry_tie(f1, f2))))) {
        f1 <- f2
      } else {
        sst0$zeta <- zinfo
      }
    }
  } else {
    f1 <- fit_once(beta)
  }
  res <- f1$res
  hyper <- f1$hyper
  crit <- f1$crit
  spec_h <- f1$spec_h

  # TWO ROUTES FOR A TERM OF THE LIKELIHOOD SHAPE. The fit above solved the
  # coefficients and the term's own parameters in one system from the start;
  # the model is fitted once more from the same start by the alternation,
  # which ends in the joint system from the point it reached
  # (statmod_alternate()), and the better of the two is kept by the rule the
  # multistarts below use. Neither route dominates (measured on the battery:
  # regime-2 -1019.855 against -1019.368, the marginal jump() -859.796
  # against -859.898). The route that won stays recorded in the structural
  # state, so the edge restarts below take it too.
  su0 <- attr(design, "structural")
  two_routes <- !is.null(sst0) && length(blocks$smooth) &&
    any(vapply(su0, function(u) identical(u$kind, "loglik"), logical(1))) &&
    all(vapply(su0, function(u) u$kind %in% c("filter", "loglik"),
               logical(1)))
  if (two_routes) {
    zfirst <- sst0$zeta
    sst0$zeta <- z_start
    sst0$route <- "alternate"
    f2 <- tryCatch(fit_once(beta), error = function(e) NULL)
    if (!is.null(f2) && edge_restart_better(f2$res, res, crit,
                                            outer_criterion,
                                            sparse_criterion)) {
      res <- f2$res
      hyper <- f2$hyper
      crit <- f2$crit
      spec_h <- f2$spec_h
    } else {
      sst0$zeta <- zfirst
      sst0$route <- NULL
    }
  }

  # SEVERAL STARTS, WHERE A STRUCTURAL TERM ASKS FOR THEM. The likelihood of
  # a regime model has several maxima and one start reaches the basin it lies
  # in: measured on MASS::geyser with three regimes, eight starts reach five
  # maxima, from -1053.39 to -1210.49. modelterms7::term_starts() says where
  # to start, the first start being the one already fitted; each further one
  # is a whole fit from the fresh coefficient start, and the best is kept by
  # the rule the edge restarts below use. Run before the edge check, so the
  # point that check reads is the winner.
  ms <- structural_multistarts(spec, design)
  if (length(ms$zeta) > 1L) {
    sst <- statmod_structural_state(design)
    zbest <- sst$zeta
    su <- attr(design, "structural")
    for (s in seq.int(2L, length(ms$zeta))) {
      sst$zeta[[ms$term]] <- ms$zeta[[s]]
      structural_start_fixups(spec, sst, su,
                              vapply(su, function(u) u$term == ms$term,
                                     logical(1)))
      f2 <- tryCatch(fit_once(beta), error = function(e) NULL)
      if (!is.null(f2) && edge_restart_better(f2$res, res, crit,
                                              outer_criterion,
                                              sparse_criterion)) {
        res <- f2$res
        hyper <- f2$hyper
        crit <- f2$crit
        spec_h <- f2$spec_h
        zbest <- sst$zeta
      } else {
        sst$zeta <- zbest
      }
    }
    sst$zeta <- zbest
  }

  # A COORDINATE AT THE EDGE OF ITS CHART CAN LOOK STATIONARY WITHOUT BEING
  # AT A MAXIMUM: its score is the bounded-scale derivative times the chart's
  # own, and the second vanishes there. Measured on MASS::geyser,
  # `waiting ~ regime(k = 2)` stopped with a transition probability at
  # 2.8e-12, a mode error of 1.1e-04 and a log-likelihood of -1134.01, while
  # moving that probability to 3.4e-4 raises it by 0.06 and the fit restarted
  # from there reaches -1099.63. A violation moves the coordinate back to the
  # edge and refits, at most three times; one that survives is recorded and
  # statmod_certificate() reports it.
  kkt <- edge_violations(spec_h, design, res$obj$split(res$par))
  for (attempt in seq_len(3L)) {
    if (!nrow(kkt)) break
    sst <- statmod_structural_state(design)
    zsave <- if (is.null(sst)) NULL else sst$zeta
    cf <- res$obj$split(res$par)
    for (r in seq_len(nrow(kkt))) {
      if (kkt$kind[r] == "structural") {
        sst$zeta[[kkt$term[r]]][[kkt$name[r]]] <- kkt$target[r]
      } else {
        j <- match(kkt$name[r], design[[kkt$param[r]]]$coef_names)
        cf[[kkt$param[r]]][[j]] <- kkt$target[r]
      }
    }
    r2 <- tryCatch(
      if (is.null(outer_criterion) && is.null(sparse_criterion)) {
        statmod_alternate(spec_h, design, blocks, hyper, inner_optimizer,
                          res$obj$stack(cf), expected, emethod, maxit, tol, vb)
      } else {
        statmod_select(spec, design, blocks, hyper, inner_optimizer,
                       outer_criterion, outer_optimizer, res$obj$stack(cf),
                       emethod, maxit, tol, vb, data,
                       sparse_criterion)
      }, error = function(e) NULL)
    if (is.null(r2) || !edge_restart_better(r2, res, crit, outer_criterion,
                                            sparse_criterion)) {
      if (!is.null(sst)) sst$zeta <- zsave
      break
    }
    res <- r2
    if (!is.null(r2$hyper)) hyper <- r2$hyper
    if (!is.null(r2$criterion)) crit <- r2$criterion
    if (length(r2$held_coef)) {
      spec_h <- S7::set_props(spec, held_coef = r2$held_coef)
    }
    kkt <-edge_violations(spec_h, design, res$obj$split(res$par))
  }

  coef <- res$obj$split(res$par)
  # a scaled SCAD or MCP carries the curvature at the fitted coefficients, so
  # that what the fit reports -- the objective, the degrees of freedom, the
  # summary -- reads the scaling the alternation settled at
  sc <- statmod_curv(spec, design, coef, expected, emethod)
  if (!is.null(sc)) spec <- sc
  fitted <- statmod_eta(spec, design, coef)$theta
  # WHICH COLUMNS THE MODEL DOES NOT IDENTIFY, from the pivot that fitted it
  # where there was one and from the information at the mode otherwise. Only
  # iwls() on a pivoting decomposition reports a dropped column, so an
  # optimizers7 method and the chol and svd routes left a
  # rank-deficient design entirely unnamed and vcov() refused the whole
  # matrix. Measured, the after-the-fact test names the SAME coordinate the
  # pivot names wherever both speak, and costs 0.26 to 0.56 per cent of a
  # realistic fit -- 3.02 per cent, or 1.2 ms, of the cheapest one there is.
  alias <- res$aliased
  Kmode <- tryCatch(statmod_penalized_at(spec, coef, design, hyper,
                                         expected, emethod),
                    error = function(e) NULL)
  if (!length(alias) && !is.null(Kmode)) {
    alias <- tryCatch(deficient_coords(Kmode, zero_kinked(spec, design, coef)),
                      error = function(e) integer(0))
  }
  # and a coefficient whose own design column has vanished -- a break-point
  # run out of the data -- is named whatever the pivot said: the pivot names
  # one of two collinear columns, and here both are empty
  if (!is.null(Kmode)) {
    alias <- sort(unique(c(alias, tryCatch(
      vanished_coords(spec, coef, design, Kmode),
      error = function(e) integer(0)))))
  }
  # A READER WHO DOES NOT PRINT THE SUMMARY IS STILL TOLD. The `NA` that
  # follows is base R's convention and is right, but on its own it reports a
  # non-identified model by omission.
  # a coefficient a term holds is not estimated by the information and is
  # not a column the model fails to identify: its column is zero by design
  alias <- setdiff(alias, term_held_stack(spec, design))
  warn_aliased(aliased_labels(spec, design, alias))
  # the terms as the fit left them, so a break-point and a nonlinear
  # parameter are read off the fitted object and prediction reapplies them
  spec <- statmod_fitted_spec(spec, coef, design)
  # WHICH INFORMATION THE STORED COUNT READS. aic() and bic() price a fit as
  # -2 l + k tau, and tau is the trace read on the information their own
  # `information` names; the count stored beside that criterion is the same
  # quantity, so it reads the same matrix. It used to read the inner fit's
  # information whatever the criterion read, so under iwls()'s expected
  # information and aic()'s observed one the two disagreed: 6.754675 against
  # 6.754023 on one smooth, and the printed count no longer reproduced the
  # criterion printed above it. At most one criterion reading tau runs, a
  # prediction-error outer criterion being refused beside a path.
  # reml() and ml() are not priced by the count, but the standard errors of
  # the same fit read the information of the criterion that estimated the
  # hyperparameters (observed_info() by default), so the count reads it too.
  # Only a fit with no outer criterion falls back to the inner information.
  tau_reader <- Filter(function(m) !is.null(m) && m@kind %in% c("aic", "bic"),
                       list(sparse_criterion, outer_criterion))
  marg_reader <- Filter(function(m) !is.null(m) && m@kind %in% c("reml", "ml"),
                        list(outer_criterion))
  edf_reader <- c(tau_reader, marg_reader)
  edf_expected <- if (length(edf_reader)) {
    info_expected(edf_reader[[1L]]@hessian)
  } else expected
  edf_emethod <- if (length(edf_reader)) {
    select_emethod(edf_reader[[1L]], NULL, emethod)
  } else emethod
  # without a criterion, a model with a filter reads the count on the matrix
  # of its joint step, which is the observed information unless the caller
  # named the expected one
  if (!length(edf_reader) &&
      any(vapply(attr(design, "structural"),
                 function(u) identical(u$kind, "filter"), logical(1)))) {
    ji <- joint_inner_info(inner_optimizer)
    edf_expected <- info_expected(ji)
    edf_emethod <- info_emethod(ji)
  }
  edf_information <- if (edf_expected) {
    distributions7::expected_info(method = edf_emethod)
  } else distributions7::observed_info()
  StatmodFit(
    spec = spec, coefficients = coef,
    structural = statmod_structural_par(spec, design), hyper = hyper,
    loglik = statmod_loglik_at(spec, coef, design),
    objective = res$value,
    edf = statmod_edf(spec, coef, design, hyper, edf_expected, edf_emethod,
                      alias),
    fitted = fitted, converged = res$converged,
    aliased = aliased_labels(spec, design, alias),
    elapsed = proc.time()[["elapsed"]] - t0,
    criterion = crit,
    history = list(
      outer = res$hist_outer,
      blocks = if (length(res$hist_blocks)) do.call(rbind, res$hist_blocks)
        else NULL,
      inner = if (length(res$hist_inner)) do.call(rbind, res$hist_inner)
        else NULL
    ),
    methods = list(smooth = inner_optimizer,
                   inner_information = inner_information,
                   # the information that the stored edf was read on, and
                   # whether a criterion prices the fit with that count
                   edf_information = edf_information,
                   edf_priced = length(tau_reader) > 0L,
                   outer = outer_criterion,
                   search = res$optimizer,
                   # the coordinates the marginal criterion left out of its
                   # determinant, which statmod_certificate() and vcov() read
                   pinned = res$pinned,
                   # coordinates at the edge of their chart that still point
                   # inward after the restarts; see edge_violations()
                   kkt = kkt,
                   sparse = vapply(blocks$sparse, function(b)
                     paste(b$param, b$term, sep = "/"), character(1)),
                   # which criterion swept the kinked penalties, and exactly
                   # which of their hyperparameters it reached: a summary
                   # cannot otherwise tell a value the path chose from one
                   # the caller held, and reported both as held
                   sparse_criterion = sparse_criterion,
                   sparse_hyper = if (is.null(sparse_criterion)) character(0)
                     else with(path_rows(spec, blocks, hyper),
                               paste(parameter, term, name, sep = "\r"))),
    call = cl
  )
}


#' The Alternation Between the Smooth Block and the Rest
#'
#' @description
#' Fits the terms whose penalties are twice differentiable in one system and
#' each remaining block by a method of its own, alternating until the objective
#' stops moving.
#'
#' @details
#' The outer search calls this function once for each hyperparameter value
#' that it tries, warm-started from the previous coefficients, so a whole fit
#' at one point of the search is one call.
#'
#' One pass fits the smooth block (beside a structural term of the filter
#' shape, the joint system of the coefficients and the term's parameters),
#' then the structural block of a term of the likelihood shape, then each
#' kinked block while the other blocks are held, and then commits the
#' schedule of any refreshable term. The passes stop when the relative change of
#' the objective falls below `tol`. The verdict is the inner fit's flag where
#' the fit does not alternate between blocks. Otherwise it also requires
#' [alternation_settled()] (the Newton decrement of the smooth part and the KKT
#' decrement of each coefficient that a kinked penalty holds at zero, both
#' within [mode_error_limit()]), and [statmod_refresh_settled()] where the model
#' has refreshable terms.
#'
#' @param spec The specification.
#' @param design The design.
#' @param blocks The split of the terms into `smooth` and `sparse`, as
#'   [statmod_blocks()] returns it.
#' @param hyper The hyperparameters, held fixed for the whole call.
#' @param inner_optimizer How the smooth block is fitted: [iwls()] or an
#'   \pkg{optimizers7} optimizer.
#' @param beta The starting coefficients, stacked into one numeric vector.
#' @param expected `TRUE` for the expected information, `FALSE` for the
#'   observed one.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#' @param maxit,tol The budget in passes and the relative tolerance on the
#'   objective, as [method_budget()] read them off `inner_optimizer`.
#' @param vb The resolved verbosity, as [verbosity()] returns it.
#' @param working_budget How many working fits [fit_working()] may
#'   take. The bootstrap excursions of [statmod_boot_restart()]
#'   pass a short one: an excursion needs to travel, not to converge.
#' @param hold_refresh `TRUE` holds every frozen break-point block at
#'   its committed positions: the fit is then an ordinary smooth fit, and
#'   the schedule is neither read off nor advanced. The outer
#'   machinery passes it at every criterion evaluation, path point and fold. A
#'   working phase
#'   inside each evaluation would multiply the cost of the fit, and a
#'   break-point that moves between evaluations makes the criterion depend on
#'   the path of the search. The positions are refined once, by the full
#'   alternation that `statmod()` runs at the chosen hyperparameters before the
#'   restarts.
#'
#' @return A list of nine:
#'   \describe{
#'     \item{`par`}{the stacked coefficients reached.}
#'     \item{`value`}{the penalized objective there, not averaged.}
#'     \item{`converged`}{a single logical. Where a term recomputes its own
#'       block, the verdict is the objective's together with that term's own
#'       flags, and the flag of the smooth block is not part of it. Beside a
#'       kinked block or a structural term it also requires
#'       [alternation_settled()] at the point. Where the fit does not
#'       alternate between blocks, it is the flag of the inner method.}
#'     \item{`obj`}{the objective object that the last pass used.}
#'     \item{`note`}{the `note` element of the last smooth fit, as
#'       [fit_smooth()] returns it.}
#'     \item{`aliased`}{the sorted stacked positions that the pivot left
#'       out.}
#'     \item{`hist_blocks`}{a list of one-row data frames, one for each block
#'       fitted in each pass, with `pass`, `block`, `objective`, `change`,
#'       `iterations` and `converged`.}
#'     \item{`hist_inner`}{a list of data frames, one for each pass in which
#'       the smooth fit recorded steps, holding the inner history with the
#'       columns `pass` and `block` added.}
#'     \item{`spec`}{the specification after the last pass, with the SCAD and
#'       MCP curvatures that the fit used.}
#'   }
#'
#' @seealso [statmod()] for the fit that this performs,
#'   [fit_smooth()] and [sparse_fit()] for the two halves of a pass, and
#'   [statmod_boot_restart()], which calls this repeatedly.
#'
#' @keywords internal
statmod_alternate <- function(spec, design, blocks, hyper, inner_optimizer, beta,
                              expected, emethod, maxit, tol, vb,
                              working_budget = 500L, hold_refresh = FALSE) {
  obj <- statmod_objective(spec, hyper, design, expected, emethod)
  value <- obj$fn(beta)
  hist_blocks <- list()
  hist_inner <- list()
  # the coordinates the scoring step's pivot left out, accumulated over the
  # passes: a column aliased in one pass is aliased in all of them, the
  # deficiency being a property of the design
  aliased <- integer(0)
  converged <- FALSE
  smooth_ok <- TRUE
  smooth_note <- NULL
  has_refresh <- length(attr(design, "refresh")) > 0L
  has_frozen <- !isTRUE(hold_refresh) &&
    any(vapply(attr(design, "refresh"),
               function(r) isTRUE(r$frozen), logical(1)))
  has_structural <- length(attr(design, "structural")) > 0L
  terms_ok <- TRUE
  frozen_ok <- TRUE
  frozen_stalled <- FALSE
  struct_ok <- TRUE

  # A structural term is fitted in the same system as the coefficients rather
  # than alternated with them: the joint gradient and the joint observed
  # information both exist, and the alternation was paying for one optimizer
  # per pass whose every iteration re-ran the recursion. This holds for a term
  # of the FILTER shape and for one of the LIKELIHOOD shape (a regime, a
  # marginal break-point), whose observed information
  # statmod_full_information() assembles from the forward-backward pass or
  # the quadrature.
  #
  # A term of the likelihood shape is also fitted by the second route, the
  # alternation followed by the joint system from the point it reached, and
  # statmod() keeps the better of the two (`sst$route`, set there). Neither
  # route dominates. On MASS::geyser, `waiting ~ regime(k = 2)` alternated
  # stopped 1.3e-3 below the maximum (the intercept and the level gap sit in
  # different blocks and move together, the alternation zig-zagged between
  # them at a ratio of about 0.6 per pass, and the relative stopping rule
  # ended it with the joint Newton decrement at 1.275e-3); both routes reach
  # the maximum. On the battery's regime-2 the joint system from the start
  # reaches a stationary point at -1019.855, and the alternation with the
  # joint step -1019.368. On the battery's marginal jump() it is the other
  # way round: -859.796 in 1.8 s against -859.898 in 22.9 s.
  sst_route <- statmod_structural_state(design)
  alternate_route <- !is.null(sst_route) &&
    identical(sst_route$route, "alternate")
  joint <- has_structural && length(blocks$smooth) && !alternate_route &&
    all(vapply(attr(design, "structural"),
               function(u) u$kind %in% c("filter", "loglik"), logical(1)))

  # A SCAD or MCP scaled by the curvature of the likelihood is scaled
  # self-consistently: the curvature is read again at the start of every
  # pass, at the coefficients the previous pass reached, and the objective and
  # the blocks are rebuilt from it, so the point the alternation settles at is
  # one where the scaling is the curvature there. See statmod_curv().
  curv_scaled <- any(vapply(blocks$sparse, function(bl) {
    pen <- bl$penalty
    "curv" %in% S7::prop_names(pen) && length(pen@curv) > 0L
  }, logical(1)))

  for (pass in seq_len(as.integer(maxit))) {
    if (curv_scaled) {
      sc <- statmod_curv(spec, design, obj$split(beta), expected, emethod)
      if (!is.null(sc)) {
        spec <- sc
        blocks <- curv_blocks(blocks, spec, design)
        obj <- statmod_objective(spec, hyper, design, expected, emethod)
        value <- obj$fn(beta)
      }
    }
    before <- value

    if (joint) {
      if (vb$blocks) {
        vb_rule(sprintf("inner pass %d: joint system", pass),
                if (S7::S7_inherits(inner_optimizer, optimizers7::optimizer))
                  vb_name(inner_optimizer) else "newton", indent = 2L)
      }
      # the caller's optimizer, where they named one. `iwls()` is a scoring
      # iteration over a design block and has no meaning over a filter's own
      # parameters, so it is what asks for the joint step's own default
      # rather than a choice to be honored here; anything else is used as
      # given, and passing NULL unconditionally made `inner_optimizer`
      # accepted and ignored for exactly the models it matters most for
      res <- statmod_fit_joint(spec, design, obj, beta, hyper,
                               optimizer = if (S7::S7_inherits(
                                 inner_optimizer, optimizers7::optimizer))
                                 inner_optimizer else NULL,
                               verbose = vb$blocks,
                               information = joint_inner_info(inner_optimizer))
      beta <- res$par
      value <- res$value
      smooth_ok <- isTRUE(res$converged)
      struct_ok <- smooth_ok
      hist_blocks[[length(hist_blocks) + 1L]] <- data.frame(
        pass = pass, block = "joint", objective = value,
        change = before - value, iterations = res$iterations,
        converged = res$converged
      )
    }

    # the smooth block, all of it at once; where a term carries a FROZEN
    # working linearization the fit is the alternation of exact working
    # fits and read-offs instead, which is the construction's own scheme
    if (!joint && length(blocks$smooth)) {
      if (has_frozen) {
        if (vb$blocks) {
          vb_rule(sprintf(
            "inner pass %d: working fits over a frozen break-point block",
            pass), vb_name(inner_optimizer, "iwls"), indent = 2L)
        }
        res <- fit_working(obj, beta, blocks$smooth, spec, design, hyper,
                           inner_optimizer, vb, tol, budget = working_budget)
        frozen_ok <- isTRUE(res$converged)
        frozen_stalled <- isTRUE(res$stalled)
      } else {
        if (vb$blocks) {
          vb_rule(sprintf("inner pass %d: smooth block, %d coefficients",
                          pass, length(blocks$smooth)),
                  vb_name(inner_optimizer, "iwls"), indent = 2L)
        }
        res <- fit_smooth(obj, beta, blocks$smooth, spec, design, hyper,
                          inner_optimizer, vb)
      }
      beta <- res$par
      value <- res$value
      smooth_ok <- isTRUE(res$converged)
      smooth_note <- res$note
      aliased <- union(aliased, res$aliased)
      hist_blocks[[length(hist_blocks) + 1L]] <- data.frame(
        pass = pass, block = if (has_frozen) "working" else "smooth",
        objective = value,
        change = before - value, iterations = res$iterations,
        converged = res$converged
      )
      if (!is.null(res$history)) {
        res$history$pass <- pass
        res$history$block <- "smooth"
        hist_inner[[length(hist_inner) + 1L]] <- res$history
      }
    }

    # the structural terms' own parameters, the coefficients held where the
    # smooth block left them
    if (has_structural && !joint) {
      v0 <- value
      if (vb$blocks) {
        vb_rule(sprintf("inner pass %d: structural terms", pass),
                "lbfgs", indent = 2L)
      }
      res <- statmod_fit_structural(spec, design, obj, beta, hyper, NULL,
                                    verbose = vb$blocks)
      value <- res$value
      struct_ok <- isTRUE(res$converged)
      hist_blocks[[length(hist_blocks) + 1L]] <- data.frame(
        pass = pass, block = "structural", objective = value,
        change = v0 - value, iterations = res$iterations,
        converged = res$converged
      )
    }

    # each non-smooth block in turn, the others held fixed
    for (bl in blocks$sparse) {
      v0 <- value
      if (vb$blocks) {
        vb_rule(sprintf("inner pass %d: %s in %s, %d coefficients", pass,
                        short_keys(bl$term), bl$param, length(bl$index)),
                indent = 2L)
      }
      res <- sparse_fit(obj, beta, bl, hyper, verbose = vb$optimizer,
                        spec = spec, design = design, expected = expected,
                        emethod = emethod)
      # WHICH route ran is reported by the route itself rather than guessed
      # beforehand: whether the compiled coordinate descent applies depends on
      # the penalty's operator being admissible AT THE STEP the fit takes, and
      # a header that asked the question at a probe step of 1 would repeat the
      # defect statmodels7 0.27.0 fixed.
      if (vb$blocks) {
        vb_say("%s: %d iterations, converged %s",
               if (is.null(res$method)) "sparse fit" else res$method,
               as.integer(res$iterations), isTRUE(res$converged), indent = 5L)
      }
      beta <- res$par
      value <- res$value
      hist_blocks[[length(hist_blocks) + 1L]] <- data.frame(
        pass = pass, block = paste(bl$param, bl$term, sep = "/"),
        objective = value, change = v0 - value,
        iterations = res$iterations, converged = res$converged
      )
    }

    # A term whose block is a JACOBIAN of its contribution has its refresh
    # committed once here, not once per objective evaluation: what advances
    # is a state of the iteration, and a schedule advancing at the speed of
    # a line search is not the one the construction was designed with. The
    # frozen working blocks were committed by fit_working() at exactly
    # these coefficients, and a second commit is not free for them: a
    # jseg's incremental read-off would take a further step and its step
    # measure would read zero, which seg_converged() reads as settled.
    if (has_refresh) {
      cf <- statmod_commit_refresh(spec, obj$split(beta), design,
                                   which = "jacobian")
      beta <- obj$stack(cf)
      value <- obj$fn(beta)
      terms_ok <- statmod_refresh_settled(spec, design, which = "jacobian")
    }

    rel <- abs(before - value) / max(1, abs(value))
    if (vb$blocks) {
      vb_say("pass %d done: objective %.8f, relative change %.3e",
             pass, value, rel, indent = 2L)
    }
    # With nothing to alternate WITH, the pass is the inner fit and the
    # verdict is the inner fit's. This line used to set converged
    # unconditionally in that case, so every model without a kinked penalty
    # reported success whatever the inner method had done -- including a run
    # that stopped on a non-finite score.
    #
    # A JOINTLY fitted structural term is that same case and was excluded
    # from it: the coefficients and the term's own parameters are solved in
    # one system, so there is no second block for a second pass to move.
    # Every such fit therefore paid one extra pass to discover a relative
    # change of exactly zero -- visible in the trace as a second pass of one
    # iteration. The alternation is still entered where a sparse block is
    # present, or a term that recomputes its own design, because those really
    # do alternate.
    # a working phase that stopped for want of step control will stop the
    # same way at the next pass, the schedule staying at its floor, so the
    # alternation ends here too instead of repeating it up to maxit times
    if (frozen_stalled) {
      converged <- FALSE
      break
    }
    lone <- !length(blocks$sparse) && !has_refresh &&
      (!has_structural || isTRUE(joint))
    if (lone) {
      converged <- isTRUE(smooth_ok) && isTRUE(struct_ok)
      break
    }
    if (rel < tol) {
      # With a term that recomputes its own block the verdict is the
      # objective's and that term's own, not the inner score's: where the
      # block is a linearization rather than a Jacobian, the score belongs
      # to the working model and does not vanish at the answer. With a
      # structural term it is the objective's and the structural block's,
      # for the other half of the same reason: the smooth block is fitted
      # with the term's parameters held, so asking its score to reach the
      # tolerance at every pass asks each conditional optimum to be
      # located to a precision the joint answer does not need.
      converged <- TRUE
      if (has_refresh) converged <- isTRUE(terms_ok) && isTRUE(frozen_ok)
      else if (!length(blocks$sparse) && !has_structural) {
        converged <- isTRUE(smooth_ok)
      }
      # Where there is something to alternate between, the verdict is read at
      # the point and not off the inner optimizers' flags: a block already at
      # its mode takes one step that moves the objective by a relative 1e-15,
      # and its stall guard then reports it as not converged. Measured beside
      # gas(1, 1), a lasso fit at its KKT point to 7e-07 read FALSE for that
      # reason alone. See alternation_readings().
      if (length(blocks$sparse) || has_structural) {
        converged <- converged &&
          alternation_settled(spec, design, obj, beta, hyper, expected,
                              emethod, aliased)
      }
      break
    }
  }
  # THE SECOND ROUTE ENDS IN THE JOINT SYSTEM, from the point the alternation
  # reached. The step is taken where nothing else alternates (no kinked block,
  # no block that recomputes itself), and kept only where it does not raise
  # the objective; the flag is then the joint system's own, as on the first
  # route.
  if (alternate_route && has_structural && length(blocks$smooth) &&
      !length(blocks$sparse) && !has_refresh) {
    v0 <- value
    zsave <- sst_route$zeta
    rp <- tryCatch(statmod_fit_joint(spec, design, obj, beta, hyper,
                                     verbose = vb$blocks),
                   error = function(e) NULL)
    if (!is.null(rp) && is.finite(rp$value) && rp$value <= value) {
      beta <- rp$par
      value <- rp$value
      converged <- isTRUE(rp$converged)
      hist_blocks[[length(hist_blocks) + 1L]] <- data.frame(
        pass = pass, block = "joint", objective = value,
        change = v0 - value, iterations = rp$iterations,
        converged = rp$converged
      )
    } else {
      # statmod_fit_joint() writes the term's parameters into the state
      sst_route$zeta <- zsave
      sst_route$key <- NULL
      sst_route$value <- NULL
    }
  }
  list(par = beta, value = value, converged = converged, obj = obj,
       note = smooth_note, aliased = sort(as.integer(aliased)),
       hist_blocks = hist_blocks, hist_inner = hist_inner, spec = spec)
}


#' The Blocks With the Penalties That a Specification Now Carries
#'
#' @description
#' Replaces the penalty of every kinked block with the one that the
#' specification carries under the same key, leaving everything else that the
#' block records where it was.
#'
#' @details
#' [statmod_alternate()] rewrites a scaled SCAD or MCP's curvature at every
#' pass. Rebuilding the blocks with [statmod_blocks()] would lose what a path
#' has written onto them, including the previous point's kink, so only the
#' penalty is replaced.
#'
#' @param blocks The blocks, as [statmod_blocks()] returns them.
#' @param spec,design The specification carrying the new penalties, and its
#'   design.
#'
#' @return `blocks`, in which the `penalty` of each kinked entry is replaced.
#'
#' @seealso [statmod_alternate()], [statmod_curv()]
#'
#' @keywords internal
curv_blocks <- function(blocks, spec, design) {
  units <- statmod_penalized(spec, design)
  for (i in seq_along(blocks$sparse)) {
    bl <- blocks$sparse[[i]]
    for (u in units) {
      if (identical(u$param, bl$param) && identical(u$key, bl$term)) {
        blocks$sparse[[i]]$penalty <- u$penalty
        break
      }
    }
  }
  blocks
}


#' Whether the Alternation Has Reached Its Point
#'
#' @description
#' Reads, at the point that the alternation reached, how far the smooth part is
#' from its mode and how far each coefficient that a kinked penalty sets to
#' zero is from its optimality condition, both in log-likelihood units, and
#' compares them with [mode_error_limit()].
#'
#' @details
#' The smooth part is every coordinate that is not a zero of a kinked
#' penalty, not in a frozen working block, not aliased, not pinned by a
#' marginal criterion and not named by [deficient_coords()] on the
#' information over the others, together with a structural term's own free
#' parameters where the model carries one.
#' Its reading is the Newton decrement \eqn{\tfrac12 g^\top K^{-1} g}, with
#' \eqn{g} the gradient of the penalized objective and \eqn{K = H + S} the
#' penalized information over those coordinates (the joint one, through
#' [statmod_joint_pieces()], beside a structural term). A coefficient that a
#' kinked penalty holds at zero is optimal where the pull of the likelihood,
#' \eqn{s_j}, does not exceed the size \eqn{\kappa_j} of the kink, and its
#' reading is the KKT decrement
#' \deqn{\frac{\max(0,\ \lvert s_j\rvert - \kappa_j)^2}{2\,c_j},}
#' the gain that a Newton step on that coordinate alone would predict, with
#' \eqn{c_j} the diagonal of the unpenalized information. The largest such
#' decrement is reported.
#'
#' A reading of the smooth part that cannot be computed is `NA`, and a
#' reading of the kinked zeros that cannot be computed is `Inf`;
#' [alternation_settled()] reads both as not settled.
#'
#' @param spec,design,obj,beta,hyper The fit's specification, design,
#'   objective, stacked coefficients and hyperparameters.
#' @param expected,emethod Which information the fit uses.
#' @param aliased The positions that the scoring step's pivot left out.
#'
#' @return `alternation_readings()` returns a named vector of `mode` and
#'   `kkt`, where `kkt` is `NA` if no coefficient sits at a kink and `Inf` if
#'   one does and its reading cannot be computed; `alternation_settled()` a
#'   single logical.
#'
#' @seealso [statmod_alternate()], [zero_readings()], [inner_mode_error()]
#'
#' @keywords internal
alternation_readings <- function(spec, design, obj, beta, hyper, expected,
                                 emethod, aliased = integer(0)) {
  out <- c(mode = NA_real_, kkt = NA_real_)
  cf <- obj$split(beta)
  nb <- length(beta)
  lab <- tryCatch(coef_labels(spec, design), error = function(e) NULL)
  if (is.null(lab) || nrow(lab) != nb) return(out)
  zero <- lab$kinked & beta == 0
  drop <- zero | frozen_block(spec, lab) | seq_len(nb) %in% aliased |
    seq_len(nb) %in% pinned_coords(spec, design)
  su <- attr(design, "structural")
  pieces <- tryCatch({
    if (length(su)) {
      jp <- statmod_joint_pieces(spec, design, obj, hyper,
                                 kinds = c("filter", "loglik"))
      u <- c(beta, jp$zeta())
      g <- jp$gr(u)
      # the joint gradient returns zeros where it cannot be evaluated, which
      # would read as a stationary point
      if (!length(g) || all(g == 0)) stop("the gradient is unavailable")
      list(g = g, K = as_dense(jp$he(u)), keep = c(which(!drop), jp$ix))
    } else {
      list(g = obj$gr(beta),
           K = statmod_penalized_at(spec, cf, design, hyper, expected, emethod),
           keep = which(!drop))
    }
  }, error = function(e) NULL)
  if (is.null(pieces)) return(out)
  g <- pieces$g
  K <- pieces$K
  keep <- pieces$keep
  keep <- keep[is.finite(g[keep]) & is.finite(diag(K)[keep])]
  # a free coordinate the information does not identify is left out, as an
  # aliased one is: a lasso over two identical columns may share the effect
  # between them, both are then free and the information over them is
  # singular, and a reading that cannot be computed made every point of the
  # path below the empty fit read as unsettled, so the path chose the empty
  # fit (lambda 41.8 on an effect of 2)
  if (length(keep) > 1L) {
    A0 <- K[keep, keep, drop = FALSE]
    dep <- tryCatch(deficient_coords((A0 + t(A0)) / 2),
                    error = function(e) integer(0))
    if (length(dep)) keep <- keep[-dep]
  }
  if (!length(keep)) {
    out[["mode"]] <- 0
  } else {
    A <- K[keep, keep, drop = FALSE]
    A <- (A + t(A)) / 2
    sol <- tryCatch(solve(A, g[keep]), error = function(e) NULL)
    if (!is.null(sol) && all(is.finite(sol))) {
      out[["mode"]] <- 0.5 * sum(g[keep] * sol)
    }
  }
  z <- which(zero)
  if (length(z)) {
    kk <- tryCatch({
      flat <- function(c) unlist(statmod_penalty_at(spec, c, hyper, design,
                                                    "gradient"),
                                 use.names = FALSE)
      pull <- g[z] - flat(cf)[z]
      bp <- beta
      bp[z] <- sqrt(.Machine$double.eps) * .Machine$double.eps
      kink <- flat(obj$split(bp))[z]
      H <- statmod_information_at(spec, cf, design, expected, emethod)
      cj <- as.numeric(Matrix::diag(H))[z]
      dec <- pmax(0, abs(pull) - abs(kink))^2 / (2 * cj)
      dec <- dec[is.finite(dec)]
      if (length(dec)) max(dec) else NA_real_
    }, error = function(e) NA_real_)
    # a coefficient sits at a kink and its reading cannot be computed: Inf,
    # so that it is not taken for the absence of a kinked zero
    out[["kkt"]] <- if (is.na(kk)) Inf else kk
  }
  out
}

#' @rdname alternation_readings
#' @keywords internal
alternation_settled <- function(spec, design, obj, beta, hyper, expected,
                                emethod, aliased = integer(0)) {
  r <- alternation_readings(spec, design, obj, beta, hyper, expected, emethod,
                            aliased)
  lim <- mode_error_limit()
  isTRUE(r[["mode"]] <= lim) && (is.na(r[["kkt"]]) || r[["kkt"]] <= lim)
}


#' The Curvature by Which a SCAD or MCP Penalty Is Scaled
#'
#' @description
#' Writes into every SCAD and MCP penalty of the model the curvature of the
#' likelihood in each of its coordinates, read at the given coefficients.
#'
#' @details
#' The curvature of coordinate \eqn{j} is the one that the coordinate descent
#' uses in its steps. It is \eqn{c_j = \sum_i w_i (x_{ij} - \bar x_j)^2}
#' where the equation's intercept is profiled out of the descent
#' ([coord_centers()]), with \eqn{w_i} the working weights that
#' [coord_working()] returns for the coordinate's own equation and
#' \eqn{\bar x_j} the column's weighted mean, and it is
#' \eqn{c_j = \sum_i w_i x_{ij}^2} where the intercept is not. The uncentered
#' diagonal of the information is not used, because it differs from the
#' curvature of the centered step when the column means are large relative to
#' their standard deviations.
#'
#' It is called at the starting coefficients, at every pass of
#' [statmod_alternate()] and once more at the fitted coefficients, so that
#' the scaling is self-consistent: at the fit, \eqn{c_j} is the curvature at
#' the fit. The step condition of the descent is met because
#' [coord_table_penalty()] writes the current step's curvature into the
#' table. Under a diagonal map \eqn{u = D\beta} the curvature is divided by
#' \eqn{d_j^2}, because the penalty is written on \eqn{u}. A penalty under any
#' other map, or reached through a sub-term, is left as it is.
#'
#' @param spec,design The specification and its design.
#' @param coef The coefficients, split by parameter.
#' @param expected,emethod Which information.
#'
#' @return The specification in which the curvature is written into the terms,
#'   or `NULL` where no penalty was written, so that the caller rebuilds its
#'   blocks only when a penalty changed.
#'
#' @keywords internal
statmod_curv <- function(spec, design, coef, expected, emethod) {
  units <- statmod_penalized(spec, design)
  todo <- Filter(function(u) {
    pen <- u$penalty
    (S7::S7_inherits(pen, penalties7::ScadPenalty) ||
       S7::S7_inherits(pen, penalties7::McpPenalty)) &&
      is.null(u$class) && !isTRUE(u$structural) && !is.null(u$index) &&
      identical(u$key, u$term)
  }, units)
  if (!length(todo)) return(NULL)
  if (!is.list(coef)) {
    np <- vapply(design[spec@distrib@params], function(d) d$npar, integer(1))
    coef <- split(coef, rep(factor(names(np), levels = names(np)), np))
  }
  moves <- length(attr(design, "refresh")) > 0L
  ep <- statmod_eta(spec, design, coef)
  dat <- if (moves) statmod_design_at(spec, coef, design) else design
  terms <- spec@terms
  wrote <- FALSE
  for (u in todo) {
    pen <- u$penalty
    p <- u$param
    wq <- coord_working(spec, ep, coef, design, p, expected, emethod)
    if (is.null(wq) && !expected) {
      wq <- coord_working(spec, ep, coef, design, p, TRUE, emethod)
    }
    if (is.null(wq)) next
    X <- coord_block(dat[[p]]$X, u$cols)
    cj <- if (curv_centered(spec, design, p)) {
      coord_colsq(X, wq$w, seq_len(ncol(X)),
                  as.numeric(xtv(X, wq$w, spec@threads)) / sum(wq$w))
    } else {
      as.numeric(wxsq(X, wq$w, spec@threads))
    }
    if (!is.null(pen@map)) {
      if (!methods::is(pen@map, "diagonalMatrix")) next
      cj <- cj / as.numeric(Matrix::diag(pen@map))^2
    }
    if (length(cj) != pen@n_coef || any(!is.finite(cj)) || any(cj <= 0)) next
    tm <- terms[[u$param]][[u$term]]
    tm@penalty <- S7::set_props(pen, curv = cj)
    terms[[u$param]][[u$term]] <- tm
    wrote <- TRUE
  }
  if (!wrote) return(NULL)
  spec@terms <- terms
  spec
}


#' Whether a Kinked Block's Curvature Is Read Centered
#'
#' @description
#' Returns, for an equation, whether its parametric intercept is profiled out
#' of the coordinate descent, so that a column's curvature is its centered
#' weighted sum of squares, as [coord_centers()] does for a descent.
#'
#' @details
#' The intercept is profiled out where the equation has one and it is not
#' held at a value. [coord_centers()] derives this from the positions that a
#' hold reaches through the objective, and this function derives it from the
#' `held_coef` of the specification directly, because that slot is available
#' in [statmod_curv()].
#'
#' @param spec,design The specification and its design.
#' @param p The equation, a parameter name.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @seealso [coord_centers()], [statmod_curv()]
#'
#' @keywords internal
curv_centered <- function(spec, design, p) {
  j <- parametric_intercept(spec, design, p)
  if (is.na(j)) return(FALSE)
  hc <- spec@held_coef[[p]]
  !(design[[p]]$coef_names[j] %in% names(hc))
}


#' The Budget and the Stopping Rule of the Alternation
#'
#' @description
#' Reads the iteration budget and the tolerance off the method that fits the
#' smooth block, which is where they are set.
#'
#' @details
#' [statmod()] has no `maxit` and no `tol`. Both are read from
#' `inner_optimizer`, so that one setting governs the inner fit and the
#' alternation.
#'
#' [iwls()] carries both directly. An \pkg{optimizers7} optimizer carries
#' `maxit` and a `criterion`, and the tolerance is read from the criterion
#' tree by [criterion_tol()]: the largest tolerance of the parts of a rule
#' built by `crit_any()`, which stops when any part holds, and the smallest
#' of the parts of a rule built by `crit_all()`, which stops only when every
#' part holds. The alternation then requires no more precision than the loop
#' inside it delivers.
#'
#' A criterion carrying no tolerance at all takes the default of
#' [optimizers7::crit_grad()], read from that function's own formals and
#' not copied as a number, so that a change there reaches this function.
#'
#' @param method [iwls()] or an \pkg{optimizers7} optimizer.
#'
#' @return A list of two: `maxit`, the budget in passes, and `tol`, the
#'   relative tolerance on the objective.
#'
#' @seealso [statmod()], [iwls()], [criterion_tol()] for the tree walk.
#'
#' @keywords internal
method_budget <- function(method) {
  if (S7::S7_inherits(method, Iwls)) {
    return(list(maxit = as.integer(method@maxit), tol = method@tol))
  }
  if (!S7::S7_inherits(method, optimizers7::optimizer)) {
    stop(paste0("'inner_optimizer' must be iwls() or an optimizers7 optimizer.\n",
                "  The budget and the stopping rule are read from it."),
         call. = FALSE)
  }
  list(maxit = as.integer(method@maxit),
       tol = criterion_tol(method@criterion))
}


#' The Settings Read From the Inner Method
#'
#' @description
#' Reads the information matrix, its approximation and the budget from the
#' inner method in one place.
#'
#' @details
#' Every route that fits the coefficients reads these settings, and reading
#' them in one place keeps a refit inside a path or a fold on the same terms
#' as the fit of the caller.
#'
#' An \pkg{optimizers7} optimizer receives the observed information. It
#' minimizes \eqn{-\ell + \rho} and [optimizers7::minimize()] documents `he`
#' as that function's second derivative, which is the observed information;
#' the expected information is a different matrix, so a Newton step on it
#' would be Fisher scoring under another name. Where a family writes neither
#' out, the expected information is also far more expensive to form than the
#' observed one. `emethod` is `NULL` for such a method, the default of
#' [iwls()] itself; it is read only where `expected` is `TRUE`, so for an
#' optimizer it records a default and not a choice.
#'
#' @param method [iwls()] or an \pkg{optimizers7} optimizer.
#' @param distrib The distribution, read only to settle an [iwls()] left at
#'   [auto_info()] through [iwls_resolve()]; such a method with no
#'   distribution is an error rather than a guess.
#'
#' @return A list of four: `expected` (a logical), `emethod` (an object),
#'   `maxit` and `tol`.
#'
#' @seealso [method_budget()] for the last two,
#'   [iwls()] for where the first two are set.
#'
#' @keywords internal
inner_settings <- function(method, distrib = NULL) {
  b <- method_budget(method)
  info <- inner_info(method)
  if (info_auto(info)) {
    if (is.null(distrib)) {
      stop("auto_info() is settled by the family: pass 'distrib', ",
           "or settle the method with inner_bind() first.", call. = FALSE)
    }
    method <- inner_bind(method, distrib, info)
    info <- inner_info(method)
  }
  list(expected = info_expected(info), emethod = info_emethod(info),
       maxit = b$maxit, tol = b$tol)
}


#' Settle the Inner Information and Bind It to the Method
#'
#' @description
#' Settles `information` against `distrib` and stores the result in the inner
#' method, so that every function that receives the method reads the same
#' matrix. For [iwls()] the result is held in the properties `hessian`,
#' `backup`, `switching` and `repair_floor` (see [iwls_resolve()]). For an
#' \pkg{optimizers7} optimizer it is held in the attribute
#' `"inner_information"`.
#'
#' @details
#' A named [distributions7::observed_info()] or [distributions7::expected_info()]
#' is stored as it is. An [auto_info()] is settled by the rule of
#' [iwls_resolve()]: the expected information where the family writes it
#' exactly and at an acceptable cost, and the observed information otherwise.
#' An optimizer has no fallback and no switch, since both belong to the loop of
#' [iwls()], so the settled object is the whole of the choice there.
#'
#' @param method [iwls()] or an \pkg{optimizers7} optimizer.
#' @param distrib The distribution to which the model is fitted.
#' @param information An object that inherits from
#'   [distributions7::Information].
#'
#' @return `method`, with the settled information bound to it.
#'
#' @seealso [inner_info()], [inner_settings()], [iwls_resolve()]
#'
#' @keywords internal
inner_bind <- function(method, distrib, information) {
  if (S7::S7_inherits(method, Iwls)) {
    return(iwls_resolve(method, distrib, information))
  }
  if (info_auto(information)) {
    exact <- distributions7::expected_hessian_exact(distrib) &&
      !distributions7::expected_hessian_costly(distrib)
    information <- if (exact) {
      distributions7::expected_info(method = information@method)
    } else distributions7::observed_info()
    attr(method, "inner_auto") <- TRUE
  } else attr(method, "inner_auto") <- NULL
  attr(method, "inner_information") <- information
  method
}


#' The Information Bound to an Inner Method
#'
#' @description
#' Returns the information object that [inner_bind()] stored in `method`: the
#' property `hessian` of an [iwls()] object, and the attribute
#' `"inner_information"` of an \pkg{optimizers7} optimizer. An optimizer that
#' was never bound has the observed information.
#'
#' @param method [iwls()] or an \pkg{optimizers7} optimizer.
#'
#' @return An object that inherits from [distributions7::Information].
#'
#' @seealso [inner_bind()]
#'
#' @keywords internal
inner_info <- function(method) {
  if (S7::S7_inherits(method, Iwls)) return(method@hessian)
  info <- attr(method, "inner_information", exact = TRUE)
  if (is.null(info)) distributions7::observed_info() else info
}


#' The Information of the Joint Step Beside a Filter
#'
#' @description
#' Returns the information that [statmod_fit_joint()] uses for a model with a
#' score-driven filter: the [distributions7::expected_info()] object bound to
#' `method` where the caller named it as `inner_information`, and
#' [distributions7::observed_info()] otherwise.
#'
#' @details
#' The joint step runs over the coefficients and the filter's own parameters,
#' where the expected information is the one conditional on the past (see
#' [statmod_full_information()]). Fisher scoring on that matrix converges
#' linearly near the mode, so an [auto_info()] that the family settled on the
#' expected information leaves the joint step on the observed one, and only
#' an [distributions7::expected_info()] that the caller named is used. The
#' settlement is read from the method: an [iwls()] settled from
#' [auto_info()] on the expected information switches to the observed one, and
#' an \pkg{optimizers7} optimizer settled from it carries the attribute
#' `"inner_auto"`.
#'
#' @param method The inner method, bound by [inner_bind()].
#'
#' @return An object that inherits from [distributions7::Information].
#'
#' @seealso [inner_info()], [statmod_fit_joint()]
#'
#' @keywords internal
joint_inner_info <- function(method) {
  # measured on 2026-10-27 over 22 filter models (misure_info_strutt): Fisher
  # scoring on the conditional information took more iterations than Newton
  # on the observed one in 19 and did not converge in 2, so auto_info()
  # stays on the observed information here
  info <- inner_info(method)
  auto <- if (S7::S7_inherits(method, Iwls)) isTRUE(method@switching) else
    isTRUE(attr(method, "inner_auto", exact = TRUE))
  if (info_expected(info) && !auto) info else distributions7::observed_info()
}


#' The Tolerance of a Criterion
#'
#' @description
#' Returns the tolerance of a criterion, descending into a combined rule: the
#' largest tolerance of its parts for a rule built by `crit_any()` and the
#' smallest for a rule built by `crit_all()`.
#'
#' @details
#' A rule built by `crit_any()` stops when any part holds, so the loop that it
#' drives is no more precise than its loosest part, and the maximum is taken.
#' A rule built by `crit_all()` stops only when every part holds, so the loop
#' reaches the precision of its tightest part, and the minimum is taken. The
#' value is the precision that the loop delivers, and a caller of
#' [method_budget()] needs that value.
#'
#' @param crit An \pkg{optimizers7} criterion, possibly a combined one built
#'   by `crit_any()` or `crit_all()`.
#'
#' @return A single number. The value is [optimizers7::crit_grad()]'s own
#'   default where the criterion carries no tolerance at all, read from that
#'   function's formals.
#'
#' @seealso [method_budget()], its caller.
#'
#' @keywords internal
criterion_tol <- function(crit) {
  nms <- S7::prop_names(crit)
  if ("criteria" %in% nms) {
    parts <- vapply(crit@criteria, criterion_tol, numeric(1))
    # crit_any() stops on its loosest part, crit_all() only on its tightest
    return(if (identical(crit@how, "all")) min(parts) else max(parts))
  }
  if ("tol" %in% nms) return(as.numeric(crit@tol)[1L])
  eval(formals(optimizers7::crit_grad)$tol)
}


#' Fit the Smooth Block
#'
#' @description
#' Runs `inner_optimizer` on the jointly fitted coefficients, holding every
#' kinked block at its current values. This is the smooth half of one pass of
#' [statmod_alternate()].
#'
#' @details
#' The objective handed to the optimizer is the full one restricted to `idx`:
#' the coefficients outside the block enter it as constants, so the value that
#' it reports is the model's objective, never a block's own.
#'
#' With [iwls()] the step is solved by the scoring loop of [iwls_fit()], which
#' reads the expected or the observed information, as the settled method
#' specifies. With an \pkg{optimizers7} optimizer the objective's `fn`, `gr` and
#' `he` are all supplied, so `newton()` receives the second derivative of the
#' objective and does not difference the gradient.
#'
#' @param obj The objective, as [statmod_objective()] returns it.
#' @param beta The current stacked coefficients, all of them.
#' @param idx The smooth block's positions within `beta`, an integer vector.
#' @param spec The specification.
#' @param design The design.
#' @param hyper The hyperparameters, held fixed.
#' @param method [iwls()] or an \pkg{optimizers7} optimizer.
#' @param vb The resolved verbosity.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`par`}{the full stacked coefficient vector, in which the block's
#'       positions are replaced.}
#'     \item{`value`}{the objective there.}
#'     \item{`converged`}{a single logical, the method's own verdict.}
#'     \item{`iterations`}{how many the method took.}
#'     \item{`history`}{for [iwls()], a data frame with one row per accepted
#'       step, or `NULL` where no step was taken; `NULL` for an
#'       \pkg{optimizers7} method.}
#'     \item{`aliased`}{the stacked positions of the coefficients that the
#'       solve left out, `integer(0)` for an \pkg{optimizers7} method.}
#'     \item{`note`}{the note of [iwls_fit()], a string that states why the
#'       run stopped where it did not stop on its rule, or that it was read
#'       as converged by the dimensionless rule, and `NULL` otherwise; `NULL`
#'       for an \pkg{optimizers7} method.}
#'   }
#'
#' @seealso [statmod_alternate()], the caller,
#'   [sparse_fit()] for the kinked half of the same pass,
#'   [iwls()] for the default method.
#'
#' @keywords internal
fit_smooth <- function(obj, beta, idx, spec, design, hyper, method, vb) {
  whole <- length(idx) == length(beta)
  # A HELD COEFFICIENT IS WRITTEN IN HERE and not merely started from: the
  # constraint is a property of the fit and not of whoever supplied the
  # starting values, so it is enforced wherever the block is solved.
  hf <- held_positions(spec, design, obj, beta)
  if (length(hf$where)) {
    beta[hf$where] <- hf$value
    # an optimizers7 method sees only the block it is handed, so a held
    # coordinate is left out of that block; iwls() drops it from its own solve
    if (!S7::S7_inherits(method, Iwls)) {
      idx <- setdiff(idx, hf$where)
      whole <- FALSE
    }
  }
  fn <- function(b) {
    v <- beta
    v[idx] <- b
    obj$fn(v)
  }
  gr <- function(b) {
    v <- beta
    v[idx] <- b
    if (whole || is.null(obj$gr_sub)) obj$gr(v)[idx] else obj$gr_sub(v, idx)
  }
  # The objective carries an EXACT second derivative -- the information plus
  # the penalty's Hessian -- and it was not being offered, so an optimizer
  # that wanted one differenced the gradient instead: newton() without 'he'
  # builds a numerical Hessian at 2p gradient evaluations an iteration. It is
  # passed to every method, as the gradient is; whether it is read is the
  # method's business, and a closure costs nothing until it is called.
  he <- function(b) {
    v <- beta
    v[idx] <- b
    as_dense(obj$he(v))[idx, idx, drop = FALSE]
  }

  if (S7::S7_inherits(method, Iwls)) {
    method <- iwls_resolve(method, spec@distrib)
    sub <- list(fn = fn, gr = gr, npar = length(idx),
                split = obj$split, stack = obj$stack)
    pieces_at <- function(b) {
      v <- beta
      v[idx] <- b
      if (whole) return(iwls_pieces(spec, design, obj$split(v), hyper, method))
      subset_pieces(spec, design, obj$split(v), hyper, method, idx)
    }
    # THE EXPECTED INFORMATION STANDS IN where the observed step cannot be
    # used, on a method iwls_resolve() settled with the fallback: an observed
    # penalized information that is not positive definite, or an observed
    # step that finds no acceptable point
    backup_at <- NULL
    if (!is.null(method@backup)) {
      m_exp <- method
      m_exp@hessian <- method@backup
      m_exp@backup <- NULL
      backup_at <- function(b) {
        v <- beta
        v[idx] <- b
        if (whole) return(iwls_pieces(spec, design, obj$split(v), hyper, m_exp))
        subset_pieces(spec, design, obj$split(v), hyper, m_exp, idx)
      }
    }
    # AND THE OBSERVED INFORMATION TAKES OVER from the expected one after
    # switch_after scoring steps, on a method iwls_resolve() settled on
    # the expected information from "auto": see iwls_fit()
    switch_at <- NULL
    sw <- method@switch_after
    if (isTRUE(method@switching) && length(sw) && is.finite(sw)) {
      m_obs <- method
      m_obs@hessian <- observed_info()
      m_obs@switching <- FALSE
      switch_at <- function(b) {
        v <- beta
        v[idx] <- b
        if (whole) return(iwls_pieces(spec, design, obj$split(v), hyper, m_obs))
        subset_pieces(spec, design, obj$split(v), hyper, m_obs, idx)
      }
    }
    # WHERE A SHARP BREAK-POINT SITS ON AN OBSERVATION the objective has a
    # kink, and iwls_fit() solves on the subspace that keeps it there once a
    # line search rejects every step: see there. The constraints are read
    # over the block's coordinates, the others being held by the alternation.
    kinks_at <- NULL
    if (has_sharp_breakpoint(spec)) {
      sp <- obj$split(seq_along(beta))
      kinks_at <- function(b) {
        v <- beta
        v[idx] <- b
        A <- kink_constraints(spec, design, obj$split(v), sp, length(beta))
        A <- A[, idx, drop = FALSE]
        A[rowSums(abs(A)) > 0, , drop = FALSE]
      }
    }
    # THE FULL HESSIAN WHERE A BLOCK IS THE JACOBIAN OF ITS TERM. The scoring
    # pieces leave out the score times the block's own derivative, which
    # mode_curvature() builds; iwls_fit() reads the sum only where a scoring
    # step is rejected or shrunk below a tenth. A frozen working
    # linearization (a sharp jump() or jseg()) has no such derivative and is
    # left to the scoring step.
    newton_at <- NULL
    rf <- attr(design, "refresh")
    if (length(rf) && !any(vapply(rf, function(r) isTRUE(r$frozen),
                                  logical(1)))) {
      params <- spec@distrib@params
      npar <- vapply(params, function(p) design[[p]]$npar, integer(1))
      offs <- cumsum(npar) - npar
      newton_at <- function(b) {
        v <- beta
        v[idx] <- b
        D <- mode_curvature(spec, design, obj$split(v), params, npar, offs,
                            sum(npar))
        he(b) + D[idx, idx, drop = FALSE]
      }
    }
    # the equations' coordinate ranges, restated in the subset's own
    # numbering: the stopping rule's scale is per equation, and the
    # objective's split speaks the full vector's coordinates
    gfull <- obj$split(seq_along(beta))
    groups <- lapply(gfull, function(ix) match(intersect(ix, idx), idx))
    groups <- Filter(length, groups)
    frozen <- match(intersect(hf$where, idx), idx)
    res <- iwls_fit(sub, beta[idx], method, spec@n_obs, pieces_at,
                    verbose = vb$inner, groups = groups, frozen = frozen,
                    backup_at = backup_at,
                    damp_on_reject = !has_sharp_breakpoint(spec),
                    kinks_at = kinks_at, newton_at = newton_at,
                    switch_at = switch_at,
                    switch_after = if (is.null(switch_at)) Inf else sw)
    out <- beta
    out[idx] <- res$par
    # the block was solved in its own numbering, so the aliased coordinates
    # are carried back to the fit's before anything outside reads them
    return(list(par = out, value = obj$fn(out), converged = res$converged,
                iterations = res$iterations, history = res$history,
                aliased = idx[res$aliased], note = res$note))
  }

  # an optimizers7 method solves nothing by a pivot and names no coordinate
  res <- optimizers7::minimize(method, fn, beta[idx], gr = gr, he = he)
  out <- beta
  out[idx] <- res@par
  list(par = out, value = obj$fn(out), converged = res@converged,
       iterations = res@iterations, history = NULL, aliased = integer(0),
       note = NULL)
}


#' The Scoring Pieces of a Subset of the Coefficients
#'
#' @description
#' Returns the penalized information over the coordinates `idx`, in the
#' assembled form that [iwls_fit()] reads when a block is solved while the
#' others are held.
#'
#' @details
#' A subset has no square-root route of its own: that route needs the whole
#' design, and the solve on a subset reads the assembled matrix anyway. So
#' the information is formed over the subset's own columns through
#' [statmod_information_at()]'s `index`, and the penalty's Hessian is added
#' on the same rows and columns, which avoids forming the whole \eqn{X'WX} to
#' keep a submatrix.
#'
#' Where the whole route went through the square roots it read
#' \eqn{R'R + C'C}, which equals \eqn{X'WX + S} up to rounding wherever the
#' factors exist and was replaced by exactly this sum where they do not.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef A named list of coefficient vectors.
#' @param hyper The hyperparameters.
#' @param method An [iwls()] method, whose `hessian` is read.
#' @param idx Integer positions of the subset in the stacked coefficients.
#'
#' @return A list with `R = NULL`, `C = NULL` and `A`, the penalized
#'   information over `idx`.
#'
#' @seealso [fit_smooth()], [iwls_pieces()]
#'
#' @keywords internal
subset_pieces <- function(spec, design, coef, hyper, method, idx) {
  if (info_auto(method@hessian)) {
    method <- iwls_resolve(method, spec@distrib)
  }
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  H <- statmod_information_at(spec, coef, design,
                              info_expected(method@hessian),
                              info_emethod(method@hessian), index = idx)
  list(R = NULL, C = NULL, A = H + S[idx, idx, drop = FALSE])
}


#' Where a Specification's Held Coefficients Sit in the Stacked Vector
#'
#' @description
#' Translates `spec@held_coef`, which names coefficients, into positions in
#' the stacked coefficient vector on which the objective works, with the
#' values beside them.
#'
#' @details
#' A hold is written by name, because the column order of the design follows
#' the formula and a term added to an equation moves every coordinate after
#' it. The translation is the objective's own split of the stacked vector, so
#' the answer is in the numbering that every consumer of that vector uses.
#'
#' A name that matches no coefficient of its equation is an error, because a
#' silent `NA` would leave the hold unenforced.
#'
#' @param spec A [StatmodSpec()], read for `held_coef`.
#' @param design Its design, read for each equation's coefficient names.
#' @param obj The objective, read for its split of the stacked vector.
#' @param beta The stacked coefficients, read for their length and for the
#'   values of the coefficients that the terms hold.
#'
#' @return A list of two numeric vectors of equal length, `where` (integer
#'   positions) and `value`. Both are empty if the model holds no
#'   coefficient.
#'
#' @seealso [fit_smooth()], which enforces the hold.
#'
#' @keywords internal
held_positions <- function(spec, design, obj, beta) {
  hc <- spec@held_coef
  # the slots a term holds (a held break-point, see
  # modelterms7::seg_hold()) are held where they stand: their columns are
  # zero, so the value does not reach the objective, and left free the
  # pivot would drop them and report them as not identified
  th <- term_held_stack(spec, design)
  if (!length(hc)) {
    return(list(where = th, value = as.numeric(beta[th])))
  }
  pos <- obj$split(seq_along(beta))
  where <- integer(0)
  value <- numeric(0)
  for (p in names(hc)) {
    v <- hc[[p]]
    if (!length(v) || is.null(design[[p]])) next
    j <- match(names(v), design[[p]]$coef_names)
    if (anyNA(j)) {
      stop(sprintf("'%s' names no coefficient of the '%s' equation.",
                   paste(names(v)[is.na(j)], collapse = ", "), p),
           call. = FALSE)
    }
    where <- c(where, pos[[p]][j])
    value <- c(value, as.numeric(v))
  }
  extra <- setdiff(th, where)
  list(where = c(where, extra), value = c(value, as.numeric(beta[extra])))
}


#' The Stacked Positions That the Terms Hold
#'
#' @description
#' Returns the positions, in the stacked coefficient vector, of the
#' coefficients that the model's terms hold ([modelterms7::term_held()]),
#' namely the slot of a held break-point. The terms are read from the
#' design's state where the fit carries one, so a term held during the fit is
#' seen at once, and from the specification otherwise.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#'
#' @return An integer vector of positions, possibly empty.
#'
#' @keywords internal
term_held_stack <- function(spec, design) {
  st <- attr(design, "state")
  params <- spec@distrib@params
  npar <- vapply(params, function(p) {
    if (is.null(design[[p]])) 0L else as.integer(design[[p]]$npar)
  }, integer(1))
  offs <- cumsum(npar) - npar
  out <- integer(0)
  for (a in seq_along(params)) {
    p <- params[a]
    for (nm in names(spec@terms[[p]])) {
      tm <- if (!is.null(st) && !is.null(st$terms[[p]][[nm]]))
        st$terms[[p]][[nm]] else spec@terms[[p]][[nm]]
      h <- tryCatch(modelterms7::term_held(tm), error = function(e) integer(0))
      if (!length(h)) next
      cols <- design[[p]]$blocks[[nm]]
      out <- c(out, offs[a] + cols[h])
    }
  }
  sort(unique(as.integer(out)))
}


#' Fit the Smooth Block Around a Frozen Working Linearization
#'
#' @description
#' Runs the iteration of Fasola, Muggeo and Kuchenhoff (2018) for a term whose
#' block is a working linearization with a frozen weight
#' ([modelterms7::jump()] and [modelterms7::jseg()]). The smooth block is
#' fitted exactly at the committed block, the break-points are read off the
#' fitted coefficients and committed, and the two alternate until the
#' read-off settles or the working objective stops moving.
#'
#' @details
#' The two functions differ only in this sequencing, which is the sequencing
#' of the \pkg{segmented} package. The fixed-point
#' iteration to which these constructions belong is not a descent method on the
#' objective of the model: its early steps under a large scaling factor move
#' uphill on purpose, and in this way it leaves a spurious optimum. Embedding
#' the read-off inside the objective of the inner optimizer would put a
#' sufficient-decrease line search in its way and stall it. During the
#' working fit the frozen blocks contribute \eqn{X\beta} and nothing else
#' (`st$working`), which makes the inner fit the plain penalized working fit
#' of the papers; the commit then advances the read-off, the scaling
#' schedule and any relabeling of crossed break-point lineages, once per
#' working fit.
#'
#' Any inner method serves: the working fit goes through [fit_smooth()],
#' which takes [iwls()] or any \pkg{optimizers7} optimizer, and the read-off
#' never moves inside the objective of the inner method. The methods differ
#' in cost. Each working fit is solved afresh at a frozen block, so a method
#' that carries exact curvature closes it in a step or two, while a
#' quasi-Newton method rebuilds its own curvature from nothing every time.
#'
#' The exit is at a fixed point of the iteration or in the cycle into which it
#' settles, judged on the working objective. The iteration stops when the
#' read-off has settled ([modelterms7::term_converged()]) with the objective
#' stalled, when the objective has stalled three times in a row, or when the
#' objective equals its value two iterations back twice (the period-two
#' cycle of the break-point that Muggeo documents, which a
#' consecutive-change rule cannot detect). Only at a fixed point does the
#' working objective coincide with the model's, so the fit does not keep a
#' best-so-far iterate: mid-travel, the committed contribution of a good
#' working value can sit orders of magnitude off the data. If the budget runs
#' out, the verdict is `FALSE`.
#'
#' A phase that cannot settle is ended early and reports `FALSE` as well.
#' This happens where every break-point that has not settled sits at the
#' floor of its scaling factor ([modelterms7::term_stalled()]) for
#' `stall_limit` consecutive working fits: the factor can shrink no further,
#' and there is nothing left with which to bring the break-point to rest. The
#' profile
#' objective of a discontinuous term is constant between consecutive
#' observations, so such a break-point passes from one observation to the
#' next for as long as the budget lasts. The default `stall_limit` of 20 ends
#' such a phase, and the stall and cycle rules are read first at every
#' working fit.
#'
#' @param obj The objective.
#' @param beta The current stacked coefficients.
#' @param idx The smooth block's indices.
#' @param spec The specification.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param method [iwls()] or an optimizer.
#' @param vb The resolved verbosity.
#' @param tol The alternation's tolerance, read for the objective-stall rule.
#' @param budget How many working fits at most.
#' @param stall_limit How many consecutive working fits with every unsettled
#'   break-point at its scaling floor end the phase.
#'
#' @return A list with the elements of [fit_smooth()], where `iterations` is
#'   the total over the working fits and `history` is `NULL`, and with
#'   `fasola`, the number of working fits taken, and `stalled`, `TRUE` where
#'   the phase was ended by `stall_limit`. In the iteration proper the list
#'   also carries a `note` element, which is `NULL`.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. *Computational Statistics*, 33, 997--1015.
#'
#' @seealso [fit_smooth()], [statmod_commit_refresh()]
#'
#' @keywords internal
fit_working <- function(obj, beta, idx, spec, design, hyper, method, vb, tol,
                        budget = 500L, stall_limit = 20L) {
  st <- attr(design, "state")
  rf <- attr(design, "refresh")
  # The working fit covers the EQUATIONS that carry a frozen block and
  # nothing else, the other equations refitted after each commit on the
  # true objective. Fitting them jointly instead put a scale equation
  # inside the working fit's warm start: a commit moves the frozen columns,
  # the stale coefficients against the new block leave residuals orders of
  # magnitude off, and the JOINT scoring step's line search then halves the
  # exact mu-solve along with the diverging scale step until the stall
  # guard fires -- measured, the working objective jumped to 9.4e5 and
  # every later working fit bailed after one iteration, so the read-offs
  # were taken from unsolved fits. Alone, an identity-link working fit's
  # first full step IS the exact linear solve of the papers whatever the
  # warm start.
  fparams <- unique(vapply(Filter(function(r) isTRUE(r$frozen), rf),
                           function(r) r$param, character(1)))
  eq <- obj$split(seq_along(beta))
  idx_a <- intersect(idx, unlist(eq[fparams], use.names = FALSE))
  idx_b <- setdiff(idx, idx_a)
  if (!length(idx_a)) {
    # every column of the frozen equations sits outside the smooth block;
    # nothing to iterate over, so the ordinary fit is the whole answer
    res <- fit_smooth(obj, beta, idx, spec, design, hyper, method, vb)
    res$fasola <- 0L
    res$stalled <- FALSE
    return(res)
  }
  conv <- FALSE
  w_prev <- Inf
  w_prev2 <- Inf
  stall <- 0L
  cyc <- 0L
  floor_run <- 0L
  stalled <- FALSE
  it_total <- 0L
  it <- 0L
  aliased <- integer(0)
  for (it in seq_len(budget)) {
    st$working <- TRUE
    st$key <- NULL
    st$value <- NULL
    res <- fit_smooth(obj, beta, idx_a, spec, design, hyper, method, vb)
    st$working <- FALSE
    st$key <- NULL
    st$value <- NULL
    it_total <- it_total + as.integer(res$iterations)
    aliased <- res$aliased
    # the commit is where the break-points move: read off the fitted
    # coefficients, schedule advanced, crossed lineages relabeled -- and
    # the coefficients continue from what the terms stored, which is what
    # makes the relabeling invisible to the caller
    cf <- statmod_commit_refresh(spec, obj$split(res$par), design,
                                 which = "frozen")
    beta <- obj$stack(cf)
    # The exit is at a FIXED POINT of the iteration or in the cycle it
    # settles into, never mid-travel: only at a fixed point does the
    # working objective coincide with the model's, so an iterate kept
    # anywhere else -- a best-so-far, an early no-improvement stop -- is a
    # linearization whose committed contribution can sit orders of
    # magnitude off the data (measured: a working value of 1073 whose true
    # objective read 1.6e6). Three rules, each a measured failure of the
    # previous draft:
    # - the read-off settled AND the working objective stalled: the step
    #   rule alone fires during the annealing tail and handed the pass
    #   loop one Fasola step per pass, 427 passes;
    # - the objective stalled three times in a row: cycling in place;
    # - the objective equal to TWO ITERATIONS BACK, twice: the period-two
    #   cycle of the break-point Muggeo documents, which keeps a
    #   consecutive-change rule from ever firing -- measured, two passes
    #   of 500 working fits each.
    w <- res$value
    settled <- statmod_refresh_settled(spec, design, which = "frozen")
    near <- function(a, b) is.finite(a) && is.finite(b) &&
      abs(a - b) < tol * (abs(b) + 1)
    stall <- if (near(w, w_prev)) stall + 1L else 0L
    cyc <- if (near(w, w_prev2)) cyc + 1L else 0L
    if ((settled && stall >= 1L) || stall >= 3L || cyc >= 2L) {
      conv <- TRUE
      break
    }
    # a break-point still moving with its scaling factor at the floor is
    # under no step control, and one that has stayed there this long
    # settles only by chance: the profile is constant between observations
    floor_run <- if (statmod_refresh_stalled(spec, design, which = "frozen"))
      floor_run + 1L else 0L
    if (floor_run >= stall_limit) {
      stalled <- TRUE
      break
    }
    w_prev2 <- w_prev
    w_prev <- w
  }
  # the other equations, once, at the committed blocks and on the true
  # objective: inside the loop they were held, a scale riding along in the
  # working fit being what made a post-commit warm start explode
  if (length(idx_b)) {
    resb <- fit_smooth(obj, beta, idx_b, spec, design, hyper, method, vb)
    beta <- resb$par
    it_total <- it_total + as.integer(resb$iterations)
    aliased <- union(aliased, resb$aliased)
  }
  if (vb$blocks) {
    vb_say("%d working fits, %s", it,
           if (conv) "settled" else if (stalled)
             "stopped: the scaling factor is at its floor and a break-point is still moving"
           else "budget exhausted", indent = 5L)
  }
  list(par = beta, value = obj$fn(beta), converged = conv,
       iterations = it_total, history = NULL, aliased = aliased, note = res$note,
       fasola = it, stalled = stalled)
}


#' Starting Coefficients
#'
#' @description
#' Builds the starting coefficients of every equation: the intercepts from the
#' intercept-only fit of the distribution, a weighted least-squares start for
#' the parametric slopes of the location where it applies, and each term's
#' own start where the term declares one.
#'
#' @details
#' Each equation's intercept starts at the intercept-only maximum likelihood
#' estimate, which [distributions7::fit_distrib()] supplies: the model from
#' which every covariate is removed is a natural place for the model with
#' covariates to begin, and it costs one small fit. A block whose term declares
#' no start of its own starts at zero, where its penalty is smallest, and the
#' penalized blocks are of this kind. [start_intercepts()] describes the
#' strategy in full.
#'
#' [distributions7::distrib_start()] is the fallback for an intercept that the
#' fit cannot supply. It returns one list per start, and each list is keyed by
#' parameter, so the value for parameter `p` is `th[[1]][[p]]`.
#'
#' A start that cannot be obtained is not an error: the fit still runs, from
#' zero on the link scale. The two routes are tried in order.
#'
#' `start` is a named list of values, a [start_strategy()], or `NULL`, which
#' is [start_intercepts()]. A strategy is evaluated once, here, before the
#' alternation between the coefficients and the hyperparameters begins, and
#' it must return one numeric vector without missing values per distribution
#' parameter, of the length of that parameter's coefficient vector. A named
#' list replaces the coefficients of the parameters that it names, and may name
#' a subset; the other parameters start as under [start_intercepts()]. A list
#' without names, a name that is not a parameter, and a vector of the wrong
#' length are errors.
#'
#' @param spec The specification.
#' @param design The design.
#' @param obj The objective.
#' @param start Optional user starting values (a named list), or a
#'   [start_strategy()].
#'
#' @return A stacked numeric vector.
#'
#' @keywords internal
statmod_start <- function(spec, design, obj, start = NULL) {
  params <- spec@distrib@params
  if (S7::S7_inherits(start, start_strategy_class())) {
    out <- start_at(start, spec, design, obj)
    if (!is.list(out) || !setequal(names(out), params)) {
      stop("a start strategy must return one vector per distribution ",
           "parameter.", call. = FALSE)
    }
    for (p in params) {
      if (length(out[[p]]) != design[[p]]$npar) {
        stop(sprintf(paste0("the start strategy returned %d values for '%s'",
                            ", which has %d coefficients."),
                     length(out[[p]]), p, design[[p]]$npar), call. = FALSE)
      }
      if (!is.numeric(out[[p]]) || anyNA(out[[p]])) {
        stop(sprintf(paste0("the start strategy returned missing or",
                            " non-numeric values for '%s'."), p),
             call. = FALSE)
      }
    }
    return(obj$stack(out[params]))
  }

  out <- start_at(start_intercepts(), spec, design, obj)
  # Damping the intercepts and letting the objective pick among the results
  # was tried and measured: on the Student t of the report it changed nothing
  # for the better, the model's own flatness deciding where the run ends
  # whatever it starts from. It is left out rather than kept as machinery
  # that earns nothing.
  if (!is.null(start)) {
    if (is.null(names(start)) ||
        any(!nzchar(names(start)))) {
      stop(paste0("'start' must be a named list, one vector per distribution",
                  " parameter, or a start strategy."), call. = FALSE)
    }
    for (p in names(start)) {
      if (!p %in% params) {
        stop(sprintf("'start' names '%s', which is not a parameter.", p),
             call. = FALSE)
      }
      if (length(start[[p]]) != design[[p]]$npar) {
        stop(sprintf("'start$%s' has length %d but '%s' has %d coefficients.",
                     p, length(start[[p]]), p, design[[p]]$npar),
             call. = FALSE)
      }
      out[[p]] <- as.numeric(start[[p]])
    }
  }
  obj$stack(out)
}


#' The Intercept of Each Equation, on the Link Scale
#'
#' @description
#' Returns the intercept-only maximum likelihood estimate where it can be had,
#' and otherwise the data-based start of [distributions7::distrib_start()].
#'
#' @details
#' Two routes are tried in order. [distributions7::fit_distrib()] fits
#' the distribution to the response with no covariates, which is the same model
#' in which every slope is zero and therefore a natural place for the fit to
#' begin; its link-scale coefficients are the intercepts.
#' [distributions7::distrib_start()] is the fallback, and its result
#' is a list of starts, and each start is keyed by parameter, so a value is
#' reached at `[[1]][[p]]`.
#'
#' The random stream is pinned and restored. The intercept-only fit takes
#' several starts (`n_start = 10`) and keeps the best, and the starts of some
#' families are random draws. The seed is therefore fixed for the length of
#' this call and the caller's stream is put back afterwards, so that the same
#' call gives the same answer twice.
#'
#' An intercept that the fit ran to the limit of its chart is replaced by the
#' data-based start of [distributions7::distrib_start()]. On a log or logit
#' chart a value past \eqn{e^{\pm 16}} is a limit of the family and not an
#' estimate: a Student t fitted to a response no heavier-tailed than a
#' gaussian puts \eqn{\nu} at a value far beyond any useful estimate. A
#' marginal criterion that estimates such a coefficient is flat in it, and a
#' search that starts there does not leave it, so the data-based start is
#' used instead. An identity chart is left alone, because a location of any
#' size is an estimate.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A named list, one entry per distribution parameter, on the link
#'   scale; an entry is `NULL` where both routes failed to supply a value.
#'
#' @seealso [statmod_start()]
#'
#' @keywords internal
statmod_intercepts <- function(spec) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  none <- stats::setNames(vector("list", length(params)), params)

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(20260810L)

  fd <- tryCatch(
    suppressWarnings(distributions7::fit_distrib(spec@distrib,
                                                 spec@response,
                                                 n_start = 10L)),
    error = function(e) NULL)
  if (!is.null(fd)) {
    e <- tryCatch(stats::coef(fd, scale = "link"), error = function(e) NULL)
    if (!is.null(e) && length(e) == length(params) && all(is.finite(e))) {
      e <- stats::setNames(as.list(as.numeric(e)), params)
      run <- vapply(params, function(p)
        !identical(links[[p]]@link_name, "identity") &&
          abs(e[[p]]) > intercept_chart_limit(), logical(1))
      if (any(run)) {
        th <- tryCatch(
          distributions7::distrib_start(spec@distrib, spec@response, 1L),
          error = function(e) NULL)
        for (p in params[run]) {
          v <- if (length(th)) th[[1L]][[p]] else NULL
          if (is.null(v) || !is.finite(v[[1L]])) next
          e[[p]] <- linkfunctions7::linkfun(links[[p]], v[[1L]])
        }
      }
      return(intercept_start_from_data(spec, e))
    }
  }

  th <- tryCatch(
    distributions7::distrib_start(spec@distrib, spec@response, 1L),
    error = function(e) NULL)
  if (is.null(th) || !length(th)) return(none)
  th1 <- th[[1L]]
  intercept_start_from_data(spec, stats::setNames(lapply(params, function(p) {
    v <- th1[[p]]
    if (is.null(v) || !is.finite(v[[1L]])) return(NULL)
    linkfunctions7::linkfun(links[[p]], v[[1L]])
  }), params))
}

#' Limit of an Intercept on a Non-Identity Chart
#'
#' @description
#' Returns the absolute value, on the link scale, past which
#' [statmod_intercepts()] reads an intercept-only estimate as a limit to which
#' the fit ran. The value is 16, which corresponds to a parameter beyond
#' \eqn{e^{16} \approx 8.9 \times 10^6} or below its reciprocal on a log
#' chart.
#'
#' @return A single positive number.
#'
#' @keywords internal
intercept_chart_limit <- function() 16

#' Data-Based Intercepts of Families That Supply Them
#'
#' @description
#' Overwrites, in the list that [statmod_intercepts()] builds, the intercept of
#' each parameter for which [distributions7::distrib_intercept_start()]
#' returns a value, carried onto that parameter's own link.
#'
#' @details
#' The intercept-only fit is the right start for almost every family. It is
#' the wrong one for the mixing weight of [distributions7::zero_inflated()],
#' and for `skewnormal1_distrib()` and the pseudo-Huber families, for which
#' [distributions7::distrib_intercept_start()] has methods of their own.
#' For the mixing weight, the overdispersion of the parent can absorb the
#' excess zeros without covariates, the weight goes to the edge of its
#' domain, and the link is flat there, so that a model with covariates
#' started from it does not leave the edge. The family returns the
#' observed proportion of zeros, and the value is mapped through the
#' parameter's link, so a probit or a cloglog chart has its own linear
#' predictor and not a logit.
#'
#' @param spec A [StatmodSpec()].
#' @param eta The named list of link-scale intercepts.
#'
#' @return `eta`, in which the answered entries are replaced.
#'
#' @seealso [statmod_intercepts()]
#' @keywords internal
intercept_start_from_data <- function(spec, eta) {
  ov <- tryCatch(distributions7::distrib_intercept_start(spec@distrib,
                                                         spec@response),
                 error = function(e) list())
  links <- spec@distrib@link_params
  for (p in intersect(names(ov), names(eta))) {
    v <- tryCatch(linkfunctions7::linkfun(links[[p]], ov[[p]]),
                  error = function(e) NULL)
    if (!is.null(v) && length(v) == 1L && is.finite(v)) eta[[p]] <- v
  }
  eta
}


#' The Penalized Information at a Point
#'
#' @description
#' Returns \eqn{H + S} as a dense matrix and sets the penalty's non-finite
#' entries to zero.
#'
#' @details
#' This is the matrix from which a penalized fit's curvature is read: the
#' model's information from [statmod_information_at()] and the penalty's
#' Hessian from [statmod_penalty_at()]. It is written once because its
#' readers, the rank test of [deficient_coords()] in [statmod()], the readings
#' of the alternation ([alternation_readings()]) and the flat-direction test
#' of the outer search, must use the same matrix at the same point.
#'
#' A kinked penalty contributes no curvature away from its kink, and any
#' non-finite entry would be the kink itself reached by a hair, so the
#' penalty passes through [zap_nonfinite()] as it does everywhere else.
#'
#' @param spec The specification.
#' @param coef A named list of coefficients, one vector per distribution
#'   parameter.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param expected Whether the expected information is used.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#'
#' @return A dense symmetric matrix over the coefficients of every equation.
#'
#' @seealso [statmod_information_at()], [statmod_penalty_at()]
#'
#' @keywords internal
statmod_penalized_at <- function(spec, coef, design, hyper, expected = TRUE,
                                 emethod = NULL) {
  H <- statmod_information_at(spec, coef, design, expected, emethod)
  S <- zap_nonfinite(statmod_penalty_at(spec, coef, hyper, design, "hessian"))
  as_dense(H + S)
}


#' Effective Degrees of Freedom, Per Term
#'
#' @description
#' Computes the effective degrees of freedom of each term: from the diagonal
#' of the model's smoother where a penalized term is not kinked, from the
#' number of coefficients away from the kink for a kinked penalty, from the
#' term's free parameters for a structural term, and through
#' \pkg{modelterms7}'s `edf()` for the remaining terms.
#'
#' @details
#' A smooth penalized term takes its share from the diagonal of
#' \eqn{(H+S)^{-1}H} for the whole model, read over the coefficients of every
#' equation at once. The computation needs the unpenalized curvature of every
#' equation, and not only the coefficients and hyperparameters of the term.
#' That curvature is computed once for all terms and not per term.
#'
#' A term may carry more than one penalty, over different parameters of its
#' own, and each has hyperparameters of its own filed under a key of its own.
#' The result has one row per term, the granularity of a table of terms,
#' and the hyperparameters are handed over keyed by the penalty names that
#' [modelterms7::term_penalties()] returns, which is the shape in which
#' `edf()` reads them. Passing the hyperparameters of one penalty for
#' a term carrying two would count the whole block against it.
#'
#' The arguments of `edf()` are passed by name: its third argument is the
#' curvature and its fourth the hyperparameters.
#'
#' @param spec The specification.
#' @param coef The coefficients.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param expected Whether the curvature is the expected information.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#' @param aliased Stacked positions of the coefficients that the fit did not
#'   estimate; they are not counted.
#'
#' @return A data frame with one row per term and the columns `parameter`,
#'   `term`, `coefficients` and `edf`, or `NULL` where the model has no term.
#'   `edf` is `NA` for a term whose count could not be read.
#'
#' @keywords internal
statmod_edf <- function(spec, coef, design, hyper, expected = TRUE,
                        emethod = NULL, aliased = integer(0)) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  H <- NULL
  rows <- list()

  # The model's smoother matrix, ONCE, over the coefficients of every
  # equation together. A term's share of the degrees of freedom is the trace
  # of its own diagonal block of this, which is what the definition says and
  # what a block-wise tr[(H_bb + S_b)^-1 H_bb] only approximates: the latter
  # drops every coupling between a term's columns and the rest of the model.
  # Measured on a gaussian with a smooth in each equation the two totals were
  # 16.98939 and 16.98885, the whole of the gap sitting on the MU smooth
  # while sigma's agreed exactly -- because the Demmler-Reinsch basis is
  # orthogonalized against the constant in the UNWEIGHTED metric and mu's
  # information carries the weights 1/sigma^2, which vary here because sigma
  # is itself modeled. The gap is small wherever the blocks are nearly
  # orthogonal and is not bounded in general.
  #
  # A kinked penalty is NOT read from this matrix. There the count is the
  # number of coefficients away from the kink, after Zou, Hastie and
  # Tibshirani, and the curvature this matrix is built from does not exist at
  # a coefficient sitting on it.
  #
  # WHERE THE MODEL CARRIES A FILTER the same definition is read on the vector
  # the model actually estimates, which is the coefficients AND the term's own
  # free parameters: joint_smoother_diag() below. Its coefficient half stands
  # in for this one, so every branch underneath is unchanged.
  # A MODEL WITH A KINKED BLOCK reads the count the path that chose its
  # hyperparameters read: the trace over the coordinates away from a kink,
  # statmod_pe()'s tau, split by term. Until 0.183.0 a kinked block was
  # counted by its non-zero coefficients and the rest by the smoother over
  # every coordinate, the zero ones included, so the count a fit reported was
  # not the count its path had scored it with: a scaled MCP on MASS::UScrime
  # scored 9.72 at lambda = 37.45 where the fit there reported 8, a SCAD point
  # 70.3 against 34.4 in BIC, and the path chose a value whose fit was not the
  # best by the criterion the fit printed. The trace is Stein's count -- a
  # coefficient in MCP's concave zone counts 1/(1 - 1/gamma) -- and for a
  # lasso alone it is the number of non-zero coefficients, as before.
  kd <- kinked_edf_diag(spec, coef, design, hyper, expected, emethod, aliased)
  if (!is.null(kd)) {
    for (a in seq_along(params)) {
      p <- params[a]
      for (nm in names(spec@terms[[p]])) {
        cols <- design[[p]]$blocks[[nm]]
        rows[[length(rows) + 1L]] <- data.frame(
          parameter = p, term = nm, coefficients = length(cols),
          edf = sum(kd[offs[a] + cols]))
      }
    }
    if (!length(rows)) return(NULL)
    return(do.call(rbind, rows))
  }

  smoother <- NULL
  zsmooth <- NULL
  zfailed <- FALSE
  if (any(vapply(statmod_penalized(spec, design),
                 function(u) !penalty_has_kink(u$penalty, u$key),
                 logical(1)))) {
    js <- joint_smoother_diag(spec, coef, design, hyper, expected, emethod)
    if (!is.null(js) && !isTRUE(js$failed)) {
      smoother <- js$beta
      zsmooth <- js$zeta
    } else {
      zfailed <- isTRUE(js$failed)
      smoother <- tryCatch({
        H <- statmod_information_at(spec, coef, design, expected, emethod)
        S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
        S <- zap_nonfinite(S)
        # a coefficient a term holds (a held break-point's slot, a column of
        # zeros) counts one: the position it stands for was estimated, on
        # the profile rather than by the derivative, as segmented counts it
        hp <- term_held_stack(spec, design)
        if (length(hp)) {
          H <- as_dense(H)
          S <- as_dense(S)
          H[hp, ] <- 0
          H[, hp] <- 0
          H[cbind(hp, hp)] <- 1
          S[hp, ] <- 0
          S[, hp] <- 0
        }
        # through solve_pd, whose equilibrated test forgives scale
        # separation from any source: a smoothing parameter a criterion
        # sends to 1e15 separates the scales without flattening a
        # direction, and LAPACK's solve on the assembled system reported it
        # as "computationally singular" -- which left the fit standing with
        # every edf missing
        Hd <- as_dense(H)
        diag(solve_pd(as_dense(H + S), "the penalized information") %*% Hd)
      }, error = function(e) NULL)
    }
  }
  for (a in seq_along(params)) {
    p <- params[a]
    for (nm in names(spec@terms[[p]])) {
      cols <- design[[p]]$blocks[[nm]]
      ent <- modelterms7::term_penalties(spec@terms[[p]][[nm]])
      v <- if (S7::S7_inherits(spec@terms[[p]][[nm]],
                               modelterms7::structural_term)) {
        # A structural term contributes no columns, so counting its block
        # gave it ZERO and every criterion built on the total was that much
        # too generous: a gas(1,1) reported 2 degrees of freedom for a model
        # carrying four. What it spends is its own parameters, one apiece,
        # being estimated and unpenalized -- less any level an intercept in
        # the same equation already carries, which is held rather than
        # estimated and is not the model's to pay for twice.
        #
        # ONE APIECE IS THE SMOOTHER'S OWN ANSWER FOR AN UNPENALIZED
        # PARAMETER, and only for one. Since 0.24.0 a filter's own parameters
        # may carry a penalty -- a development of a loading over groups, a
        # covariance class -- and there the rule counts a coordinate the prior
        # has shrunk to nothing as a whole degree of freedom. Measured on the
        # joint smoother at the fitted point, the unpenalized coordinates come
        # back at 1.000000000000000 with a gap of 0.000e+00, so the two
        # readings are the same rule where the old one was right, while eight
        # shrunk deviations read between 0.080 and 0.247 against a count of
        # eight, and sixteen under a covariance class read between -0.0003 and
        # 0.829 against a count of sixteen.
        if (!is.null(zsmooth)) sum(zsmooth)
        else if (zfailed) {
          # THE MODEL HAS A JOINT VECTOR AND IT COULD NOT BE READ, so there is
          # no count here rather than the old one under a new name. Measured,
          # this fires where the joint matrix is indefinite -- a fit that has
          # already failed, reporting no criterion -- and the rule it replaces
          # returned the parameter count, 20.00 on a model whose penalized
          # count is 10.40, with nothing saying the definition had changed.
          warning(sprintf(paste0(
            "the degrees of freedom of '%s' in '%s' cannot be read: the ",
            "curvature over\n  the coefficients and the term's own parameters ",
            "is not positive definite at\n  this point, so the model's ",
            "smoother has no diagonal there. The fit itself\n  stands, but the ",
            "count is missing, and so is every criterion built on it."),
            nm, p), call. = FALSE)
          NA_real_
        } else {
          zn <- modelterms7::term_params(spec@terms[[p]][[nm]])
          st <- statmod_structural_state(design)
          as.numeric(length(setdiff(zn, st$held[[nm]])))
        }
      } else if (!length(ent) && is.null(smoother)) {
        # An ALIASED column is not a parameter the fit spent: the pivot left
        # it out and its coefficient is reported as missing, so counting it
        # here would put one degree of freedom too many into every criterion
        # built on the total. Measured on a design with one duplicated
        # column, the count read 5 where glm's rank plus the scale is 4.
        as.numeric(length(cols) -
                     sum((offs[a] + cols) %in% aliased))
      } else if (!length(ent)) {
        # An unpenalized block is not automatically worth one per column
        # either: coupled to a penalized one it takes whatever share of the
        # smoother matrix its own directions carry, and reading that off is
        # the same rule as everywhere else rather than an exception.
        sum(smoother[offs[a] + cols])
      } else if (!is.null(smoother) &&
                 !any(vapply(ent, function(e)
                   penalty_has_kink(e$penalty, nm), logical(1)))) {
        sum(smoother[offs[a] + cols])
      } else {
        if (is.null(H)) {
          H <- statmod_information_at(spec, coef, design, expected, emethod)
        }
        idx <- offs[a] + cols
        th <- if (length(ent) == 1L) {
          as.list(hyper[[p]][[statmod_entry_key(nm, ent, ent[[1L]])]])
        } else {
          stats::setNames(
            lapply(ent, function(e)
              as.list(hyper[[p]][[statmod_entry_key(nm, ent, e)]])),
            vapply(ent, function(e) e$name, character(1)))
        }
        tryCatch(
          modelterms7::edf(spec@terms[[p]][[nm]],
                           coef = coef[[p]][cols],
                           hessian = H[idx, idx, drop = FALSE],
                           theta = th),
          error = function(e) {
            # the count is not worth abandoning a fit for, but a message that
            # named none of this would leave the reader to guess why a column
            # went missing
            warning(sprintf("edf() failed for '%s' in '%s': %s", nm, p,
                            conditionMessage(e)), call. = FALSE)
            NA_real_
          })
      }
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = p, term = nm, coefficients = length(cols),
        edf = as.numeric(v)[1L]
      )
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}


#' The Count of Each Coordinate Where a Block Has a Kink
#'
#' @description
#' Returns the diagonal of \eqn{(H_{AA} + S_{AA})^{-1}H_{AA}} over the
#' coordinates \eqn{A} away from a kink, and zero elsewhere. This is the count
#' with which [statmod_pe()] prices a point of a path, one coordinate at a
#' time.
#'
#' @details
#' The total is [statmod_pe()]'s \eqn{\tau} at the same coefficients: where
#' the penalty's Hessian vanishes on \eqn{A} each active coordinate counts one,
#' which is the number of non-zero coefficients of a lasso, and otherwise the
#' trace is read through the same Cholesky factor. A coordinate in the
#' concave zone of SCAD or MCP has a negative penalty curvature and counts more
#' than one, which is Stein's count for those penalties.
#'
#' It returns a value only for a model with a kinked block and no structural
#' term, and returns `NULL` where the penalty has no curvature over the active
#' coordinates and the count of one per column cannot be certified
#' ([design_count_exact()] is `FALSE`, where [statmod_pe()] reads a rank in
#' place of a count), or where the penalized curvature of the active
#' coordinates is not positive definite, so that its Cholesky factor does not
#' exist; [statmod_edf()] then reads its other rules.
#'
#' @param spec,coef,design,hyper The specification, the coefficients, the
#'   design and the hyperparameters.
#' @param expected,emethod Which information.
#' @param aliased Stacked positions that the fit did not estimate.
#'
#' @return A numeric vector over the stacked coefficients, or `NULL`.
#'
#' @seealso [statmod_pe()], [statmod_edf()]
#'
#' @keywords internal
kinked_edf_diag <- function(spec, coef, design, hyper, expected, emethod,
                            aliased = integer(0)) {
  if (!is.null(statmod_structural_state(design))) return(NULL)
  blocks <- statmod_blocks(spec, design)
  if (!length(blocks$sparse)) return(NULL)
  params <- spec@distrib@params
  beta <- unlist(lapply(params, function(p) as.numeric(coef[[p]])),
                 use.names = FALSE)
  act <- statmod_active(spec, blocks, beta, hyper)
  act[aliased] <- FALSE
  ai <- which(act)
  d <- numeric(length(beta))
  if (!length(ai)) return(d)
  S <- zap_nonfinite(statmod_penalty_at(spec, coef, hyper, design, "hessian"))
  S_aa <- as_dense(S[ai, ai, drop = FALSE])
  if (!any(S_aa != 0)) {
    if (!design_count_exact(design)) return(NULL)
    d[ai] <- 1
    return(d)
  }
  H <- as_dense(statmod_information_at(spec, coef, design, expected, emethod,
                                       index = ai))
  P <- tryCatch(chol2inv(chol(H + S_aa)), error = function(e) NULL)
  if (is.null(P)) return(NULL)
  d[ai] <- rowSums(P * t(H))
  d
}


#' The Model's Smoother Over the Coefficients and a Filter's Parameters
#'
#' @description
#' Returns the diagonal of \eqn{F = (K+S)^{-1}K} on the vector that a model
#' carrying a score-driven filter actually estimates: the coefficients of
#' every equation followed by the term's own free parameters.
#'
#' @details
#' The effective degrees of freedom are the trace of the model's smoother and a
#' term's share is the trace of its own diagonal block. Where every unknown is
#' a coefficient that matrix is \eqn{(H+S)^{-1}H} and [statmod_edf()] reads it
#' directly. A filter's parameters are estimated beside the coefficients and
#' may carry a penalty of their own, so for such a model the same definition is
#' read on the joint vector: \eqn{K} is [statmod_full_information()] and
#' \eqn{K+S} the matrix whose determinant the criterion takes,
#' [statmod_marginal_full()].
#'
#' # The count of an unpenalized parameter
#'
#' Writing \eqn{P} for the coordinates that some penalty covers, \eqn{S} has a
#' zero row and a zero column outside \eqn{P}, so
#' \deqn{F = M^{-1}K = M^{-1}(M - S) = I - M^{-1}S,}
#' and for \eqn{j} outside \eqn{P} the whole column \eqn{S_{\cdot j}} vanishes,
#' which makes \eqn{F_{jj}} exactly one. The rule that counts one degree of
#' freedom per free parameter of a structural term, on the reading that those
#' parameters are estimated and unpenalized, therefore agrees with this matrix
#' wherever that rule is right. The count differs for a coordinate that
#' a prior shrinks, which a filter's own parameters may be (for example the
#' deviations of a developed loading, or the parameters under a covariance
#' class); such a coordinate counts less than one.
#'
#' A structural term that mixes over latent states rather than shifting the
#' predictor has no such joint determinant, because [statmod_marginal_full()]
#' returns `NULL` for it, so the count falls back to the coefficient-only
#' reading.
#'
#' The joint route reads the observed information by default, and the
#' expected information conditional on the past with `expected = TRUE` (see
#' [statmod_full_information()]); [statmod_marginal_full()] assembles from
#' [modelterms7::term_curvature()], and the criterion of such a model already
#' reads that same matrix. A model carrying no filter never reaches here and is
#' untouched.
#'
#' @param spec The specification.
#' @param coef The coefficients.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param expected,emethod Whether the matrix is the expected information
#'   conditional on the past, and its `method`.
#'
#' @return A list with `beta`, the diagonal over the stacked coefficients,
#'   `zeta`, the diagonal over the structural term's free parameters, named as
#'   the term names them, and `failed`. The result is `NULL` where the model
#'   carries no joint vector at all, and a list whose `failed` is `TRUE` where
#'   it carries one that could not be read. The two results differ: with the
#'   first, the caller may fall back on a rule that is right for such a model,
#'   and with the second, the caller has no reading, so that the count is
#'   reported missing rather than taken from a rule that means something else.
#'
#' @seealso [statmod_edf()] and [summary_class_blocks()], its callers.
#'
#' @keywords internal
joint_smoother_diag <- function(spec, coef, design, hyper, expected = FALSE,
                                emethod = NULL) {
  sst <- statmod_structural_state(design)
  if (is.null(sst)) return(NULL)
  su <- Filter(function(u) identical(u$kind, "filter"),
               attr(design, "structural"))
  if (!length(su)) return(NULL)
  key <- su[[1L]]$term
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(design[spec@distrib@params], function(d) d$npar,
                   integer(1)))
  d <- tryCatch({
    K <- statmod_full_information(spec, coef, design, expected, emethod)
    M <- statmod_marginal_full(spec, design, coef, hyper, expected = expected,
                               emethod = emethod)
    if (is.null(K) || is.null(M) || nrow(M) != nb + length(free) ||
        !identical(dim(K), dim(M))) {
      NULL
    } else {
      # solve_pd for the reason the coefficient route gives: a smoothing
      # parameter a criterion sends to 1e15 separates the scales without
      # flattening a direction, and its equilibrated test forgives that.
      # Where it refuses, the matrix is measured INDEFINITE rather than
      # ill-conditioned -- at one such point its smallest eigenvalue is
      # -7.1e-05, chol() refuses it too, and a plain solve returns a diagonal
      # with a negative entry, which is not a smoother's -- so the refusal is
      # right and there is nothing to relax.
      diag(solve_pd(as_dense(M), "the joint penalized information") %*%
             as_dense(K))
    }
  }, error = function(e) NULL)
  # A NEGATIVE COUNT IS NOT A COUNT, where M is accepted and K is indefinite.
  # Measured on a filter beside a random effect that stopped unconverged,
  # K's smallest eigenvalue was -3.1e-4 against 4.3e5, M passed, and the
  # filter's entries read -1.04 to -1.56, its count -3.21; it is reported
  # missing, as where M is refused. The rule is on the count and not on the
  # entries (Giovanni, 2026-10-05): entries of either sign near zero are
  # ordinary, -6e-8 for a shrunk deviation and -0.067 on a converged panel
  # whose K is indefinite as an observed information over penalized
  # deviations often is (count 3.73).
  if (!is.null(d) && sum(d[nb + seq_along(free)]) < 0) d <- NULL
  # THE TWO WAYS OF HAVING NOTHING ARE DIFFERENT and the caller must tell them
  # apart: a model with no joint vector falls back to a rule that is right for
  # it, while one whose joint matrix cannot be read has no reading at all, and
  # the old count would be a different quantity reported under the same name.
  if (is.null(d)) return(list(beta = NULL, zeta = NULL, failed = TRUE))
  z <- d[nb + seq_along(free)]
  names(z) <- free
  list(beta = d[seq_len(nb)], zeta = z, failed = FALSE)
}


#' Resolve the Verbosity Setting
#'
#' @description
#' Turns a level or a named logical vector into the four switches `outer`,
#' `blocks`, `inner` and `optimizer` that the fit reads.
#'
#' @details
#' Each level adds output: `1` shows the outer search and the alternation
#' (switches `outer` and `blocks`), `2` adds the inner method's own iterations
#' (`inner`), and `3` adds the traces of the optimizers (`optimizer`). The
#' named form sets the four switches independently, which allows the
#' alternation to be shown while a chatty inner optimizer is silent. A switch
#' that is left out is off.
#'
#' @param verbose A number from 0 to 3, or a named logical vector with any of
#'   `outer`, `blocks`, `inner`, `optimizer`.
#'
#' @return A list of four logicals named `outer`, `blocks`, `inner` and
#'   `optimizer`.
#'
#' @keywords internal
verbosity <- function(verbose) {
  switches <- c("outer", "blocks", "inner", "optimizer")
  if (is.logical(verbose) && !is.null(names(verbose))) {
    bad <- setdiff(names(verbose), switches)
    if (length(bad)) {
      stop(sprintf(paste0("'verbose' names '%s'. The switches are %s."),
                   bad[1L], paste0("'", switches, "'", collapse = ", ")),
           call. = FALSE)
    }
    # a switch left out is off, so that naming one of the four is a complete
    # request rather than an error
    on <- function(nm) nm %in% names(verbose) && isTRUE(verbose[[nm]])
    return(list(outer = on("outer"), blocks = on("blocks"),
                inner = on("inner"), optimizer = on("optimizer")))
  }
  if (!is.numeric(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a number from 0 to 3 or a named logical vector.",
         call. = FALSE)
  }
  lv <- as.integer(verbose)
  list(outer = lv >= 1L, blocks = lv >= 1L, inner = lv >= 2L,
       optimizer = lv >= 3L)
}


#' Whether a Model Carries a Sharp Break-Point Term
#'
#' @description
#' Returns `TRUE` where some equation carries [modelterms7::seg()],
#' [modelterms7::jump()] or [modelterms7::jseg()] with neither `smoothed` nor
#' `marginal = TRUE`, a term whose contribution has a kink in the break-point
#' at the observations.
#'
#' @details
#' [fit_smooth()] reads it to keep [iwls_fit()] from escalating the
#' Levenberg damping on a rejected step: at the minimum of an objective with
#' a kink every Gauss-Newton step is rejected whatever the damping, so the
#' escalation spends eight attempts at every inner fit and ends where a stop
#' would have.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A single logical.
#'
#' @seealso [iwls_fit()], [fit_smooth()]
#'
#' @keywords internal
has_sharp_breakpoint <- function(spec) {
  for (p in names(spec@terms)) {
    for (tm in spec@terms[[p]]) {
      if (S7::S7_inherits(tm, modelterms7::SegTerm) &&
          is.null(tm@spec$smoothed) && !isTRUE(tm@spec$marginal)) {
        return(TRUE)
      }
    }
  }
  FALSE
}
