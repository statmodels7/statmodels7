#' @include blocks.R
NULL

#' A Fitted Model
#'
#' @description
#' What [statmod()] returns. It holds three things: the specification kept
#' whole, so the model can be read back and reapplied to new data; the
#' coefficients, hyperparameters and structural parameters the fit reached;
#' and the record of how it got there.
#'
#' Reach the contents through the accessors rather than the properties where
#' one exists: [coef.StatmodFit()], [vcov.StatmodFit()],
#' [predict.StatmodFit()], [summary.StatmodFit()], [logLik.StatmodFit()].
#'
#' @param spec The [StatmodSpec()] that was fitted, with every refreshable
#'   term as the fit left it.
#' @param coefficients A named list, one numeric vector per distribution
#'   parameter.
#' @param structural A named list, one entry per structural term, holding
#'   that term's own parameters on the unconstrained scale.
#' @param hyper The hyperparameters, keyed by distribution parameter and then
#'   by term.
#' @param loglik The maximized weighted log-likelihood, a single number. The
#'   **conditional** one, read at the fitted coefficients.
#' @param objective The penalized objective there, unaveraged.
#' @param edf A data frame of effective degrees of freedom, one row per term.
#' @param fitted A named list of the fitted parameters, one vector per
#'   distribution parameter.
#' @param converged A single logical: whether every loop stopped on its own
#'   rule.
#' @param elapsed The elapsed time in seconds.
#' @param criterion The marginal criterion at the estimated hyperparameters,
#'   `NA` when no criterion ran.
#' @param history A list of data frames recording the run: `outer`, `blocks`
#'   and `inner`.
#' @param methods A list naming what fitted each block, and the optimizer the
#'   outer search used.
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
#' # One coefficient vector per parameter of the family.
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
    elapsed = S7::class_numeric,
    criterion = S7::class_numeric,
    history = S7::class_list,
    methods = S7::class_list,
    call = S7::class_any
  )
)


#' Fit a Model
#'
#' @description
#' Fits a distributional regression. One formula carries an equation per
#' parameter of the response distribution, so a scale or a shape is modeled
#' as readily as a mean; the terms those equations name are assembled into a
#' penalized likelihood and fitted, and any smoothing parameters or prior
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
#' A parameter with no equation of its own gets an intercept, so
#' `y ~ x` on a two-parameter family fits a constant scale.
#' [statmod_equations()] does the split; the recovery is not the obvious one,
#' since R's precedence makes the whole right-hand side of a three-equation
#' formula the last term alone.
#'
#' The terms available are \pkg{modelterms7}'s: [modelterms7::s()],
#' [modelterms7::te()], [modelterms7::random()], [modelterms7::ridge()],
#' [modelterms7::lasso()], [modelterms7::seg()], [modelterms7::nl()],
#' [modelterms7::gas()] and the rest. Bare covariates collapse into one
#' parametric block.
#'
#' # The fitting scheme
#'
#' The terms split in two by a property each one already reports.
#'
#' A term whose penalty is **twice differentiable** in its coefficients, an
#' unpenalized block, a ridge, a spline or a random effect, is estimated in
#' one system by `inner_optimizer`. Their joint curvature exists, and using
#' it is what closes a fit in a handful of iterations.
#'
#' A term whose penalty has a **kink**, a lasso, a SCAD or an MCP, is
#' estimated by a coordinate descent of its own with everything else held
#' fixed.
#'
#' The fit alternates between the two until the objective and every block
#' stop moving.
#'
#' # The objective is unaveraged
#'
#' Minus the weighted log-likelihood, plus the penalties at full size. A
#' penalty is a negative log-prior, and a posterior adds a log-likelihood and
#' a log-prior at full size; averaging one of them would make a
#' hyperparameter mean something that depends on \eqn{n}. What is scaled
#' instead is the stopping rule, so a threshold means the same at
#' \eqn{n = 10} and at \eqn{n = 10^7}.
#'
#' # The budget and the stopping rule belong to the method
#'
#' There is no `maxit` and no `tol` here. Both are set on `inner_optimizer`,
#' which is [`iwls(maxit =, tol =)`][iwls] or an \pkg{optimizers7} optimizer
#' with its own `maxit` and `criterion`, and the alternation reads them from
#' there. A second copy here would let a caller set both and be obeyed by
#' neither.
#'
#' # Every hyperparameter is estimated unless its term holds it
#'
#' Which ones are held is said by the **term** that carries the penalty:
#' `lasso(x, lambda = 3)`, `ridge(x, sigma = 0.5)`, `s(x, lambda = 2)`,
#' `enet(x, alpha = 0.5)`. Everything left `NULL`, which is each term's
#' default, is chosen from the data. The term is where the penalty is named
#' and so is where that belongs; an argument here saying the same thing would
#' be read by nobody whenever the two disagreed.
#'
#' The **smooth** hyperparameters go to `outer_criterion`, [reml()] by
#' default, with `outer_optimizer` searching over them and the coefficients
#' refitted at each point.
#'
#' A **kinked** penalty is a different instrument with its own argument.
#' `sparse_criterion`, [bic()] by default, sweeps it along a path of its own
#' values, from the kink that empties the block down to a fraction of it. The
#' penalized mode is only piecewise smooth in such a hyperparameter, turning
#' a corner whenever a coefficient joins the active set or leaves it, so a
#' criterion read there inherits the corners and a gradient search would read
#' a slope about to change.
#'
#' The top of that path depends on the data **and on the rest of the model**:
#' it is the kink that empties the block, found at the coefficients in hand
#' never at a refitted null, so the other terms' fits enter it.
#'
#' Where a model carries both kinds the path is outside and the marginal
#' criterion is estimated inside each of its points, so a smoothing parameter
#' can come from REML and a lasso's \eqn{\lambda} from BIC in one fit.
#'
#' `outer_criterion` runs if and only if the model carries a smooth penalty.
#' On an ordinary `y ~ x`, or on a model whose only penalty is kinked, there
#' is nothing for it to estimate and it is not run, so typing the default
#' changes nothing. `outer_criterion = NULL` holds every smooth
#' hyperparameter where its term left it.
#'
#' # Verbosity
#'
#' Three levels, naming the loops: `1` the outer search and the alternation,
#' `2` the inner method's own iterations, `3` the optimizers' traces as well.
#' A named form is accepted too, `verbose = c(outer = TRUE, blocks = FALSE)`,
#' for watching the hyperparameters move while silencing a chatty inner
#' optimizer.
#'
#' @param formula The model formula.
#' @param distrib A \pkg{distributions7} distribution object.
#' @param data A data frame.
#' @param weights Optional prior weights, taken as given and not normalized.
#' @param offsets Optional named list of offsets, one per parameter.
#' @param inner_optimizer How the smooth block is fitted: [iwls()] or
#'   an \pkg{optimizers7} optimizer.
#' @param outer_criterion How the smooth hyperparameters are estimated:
#'   [reml()] (the default), [ml()],
#'   [aic()], [bic()], [cv()], or
#'   `NULL` to hold them where they are.
#' @param sparse_criterion How the hyperparameter of a kinked penalty --
#'   lasso, scad, mcp -- is chosen: [bic()] (the default),
#'   [aic()], [cv()], or `NULL` to hold it.
#'   A marginal criterion is rejected here, being read at a mode that sits on
#'   the kink.
#' @param outer_optimizer The optimizer that searches over them, or
#'   `NULL` to let the availability of the exact gradient decide.
#' @param start Where the fit begins: a named list of coefficients, a
#'   [start_strategy()] such as [start_search()], or
#'   `NULL` for [start_intercepts()]. A strategy is asked once,
#'   before the alternation between the coefficients and the hyperparameters
#'   begins, which is why a global search belongs here instead of in
#'   `inner_optimizer`: there it would rerun at every hyperparameter the
#'   outer search tried.
#' @param linpar_control How the unpenalized parametric block is built, as
#'   [linpar_options()] returns it: the storage and the contrasts
#'   for its factors. It reaches the implicit term, the one the bare
#'   covariates collapse into and which a caller never writes; a
#'   [modelterms7::linpar()] written out takes them directly.
#'   The argument and the function are named differently on purpose, as
#'   [stats::glm()]'s `control` and
#'   [stats::glm.control()] are.
#' @param verbose A level from 0 to 3, or a named logical vector.
#' @param threads How many threads the fit may use, as
#'   [numericals7::n_threads()] constructs it. The default,
#'   `n_threads(1)`, is sequential and takes exactly the sequential
#'   code path. A larger count reaches the compiled per-observation kernels
#'   of the family and the dense assembly products as an argument; below a
#'   kernel's measured internal threshold it stays sequential whatever the
#'   count says, and a sparse design keeps its \pkg{Matrix} route. The
#'   object's `workers` fans the units that are independent by
#'   construction out over separate R processes -- the folds of a
#'   cross-validation, and the combinations of a kinked path's product
#'   grid, each of which restarts its warm chain from the sweep's own
#'   starting coefficients -- each unit fitting sequentially, so the two
#'   levels never nest. The points within one chain stay sequential:
#'   measured, a point paid cold costs 2.2-3.2 times the warm chain, so
#'   splitting a chain would slow the single-process default or make the
#'   result depend on the count. The result does
#'   not depend on either count, bit for bit: every parallel region
#'   decomposes its work over the elements of its output and never splits
#'   a reduction, and the folds are seeded per fold and collected in fold
#'   order.
#' @param ... Not used, and reported. `hyper` was removed from here: a
#'   hyperparameter is held in the term that carries the penalty, and a
#'   second place to say so would be read by nobody whenever the two
#'   disagreed.
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
#' # An ordinary regression: the scale gets an intercept it was not given.
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' coef(fit)
#'
#' # The scale modeled too, which is what the framework is for. The data
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
                    inner_optimizer = iwls(), outer_criterion = reml(),
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
    if (length(nm) && any(nm == "hyper")) {
      stop(paste0("'hyper' has been removed. A hyperparameter is held in",
                  " the term that\n  carries the penalty --",
                  " lasso(x, lambda = 3), ridge(x, sigma = 0.5),\n",
                  "  s(x, lambda = 2), enet(x, alpha = 0.5) -- and every",
                  " one left NULL,\n  which is the default, is",
                  " estimated."), call. = FALSE)
    }
    stop(sprintf("unused argument (%s)", nm[[1L]]), call. = FALSE)
  }
  if (!is.list(linpar_control)) {
    stop("'linpar_control' must be the list linpar_options() returns.",
         call. = FALSE)
  }
  spec <- statmod_spec(formula, distrib, data, weights, offsets,
                       linpar = linpar_control)
  spec@threads <- n_thr
  spec@workers <- numericals7::worker_count(threads)
  design <- statmod_design(spec)
  hyper <- statmod_hyper_start(spec, design)
  blocks <- statmod_blocks(spec, design)

  cfg <- inner_settings(inner_optimizer)
  expected <- cfg$expected
  approx <- cfg$approx
  obj <- statmod_objective(spec, hyper, design, expected, approx)
  beta <- statmod_start(spec, design, obj, start)

  if (vb$blocks || vb$outer) {
    vb_rule(sprintf("%s, %d observations", spec@distrib@distrib_name,
                    spec@n_obs), char = "=")
    vb_say("inner    %s", vb_name(inner_optimizer, "iwls"), indent = 5L)
    vb_say("outer    %s", if (is.null(outer_criterion)) "held" else
      paste0(toupper(outer_criterion@kind), " (",
             outer_criterion@hessian, " information)"), indent = 5L)
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
    stop("'outer_criterion' must be reml(), ml(), aic(), bic(), cv() or NULL.",
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
                "
  penalty is chosen by a path over its own values, and",
                " reml() and ml()
  are read at a mode that sits on the",
                " kink, where the curvature they
  need does not exist."),
         call. = FALSE)
  }
  if (!is.null(sparse_criterion) &&
      !nrow(path_rows(spec, blocks, hyper, sparse_criterion))) {
    sparse_criterion <- NULL
  }
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
  if (!is.null(outer_criterion)) {
    nothing <- !nrow(outer_hyper_index(spec, blocks))
    if (nothing) outer_criterion <- NULL
  }

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
                " model also carries a penalty
  with a kink, whose",
                " hyperparameter is chosen by a path.
  A prediction-error",
                " criterion cannot be nested inside that path: it
  would",
                " score the same quantity at both levels. Use reml() or",
                " ml() for
  the smooth penalties and 'sparse_criterion'",
                " for the kinked ones, or
  hold one of the two with",
                " 'hyper'."), call. = FALSE)
  }

  if (is.null(outer_criterion) && is.null(sparse_criterion)) {
    res <- statmod_alternate(spec, design, blocks, hyper, inner_optimizer, beta,
                             expected, approx, maxit, tol, vb)
    crit <- NA_real_
  } else {
    res <- statmod_select(spec, design, blocks, hyper, inner_optimizer,
                          outer_criterion, outer_optimizer, beta, approx, maxit,
                          tol, vb, data, weights, offsets, sparse_criterion)
    hyper <- res$hyper
    crit <- res$criterion
    # Everything inside the selection held the frozen break-point blocks at
    # their committed positions -- a break-point moving between criterion
    # evaluations makes the criterion path-dependent, and the phase's own
    # flags read as unavailable points to the search -- so the positions are
    # refined here, once, at the chosen hyperparameters.
    frozen <- any(vapply(attr(design, "refresh"),
                         function(r) isTRUE(r$frozen), logical(1)))
    if (frozen) {
      ro <- tryCatch(
        statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                          res$par, expected, approx, maxit, tol, vb),
        error = function(e) NULL)
      if (!is.null(ro) && is.finite(ro$value)) {
        res[c("par", "value", "converged", "obj",
              "hist_blocks", "hist_inner")] <-
          ro[c("par", "value", "converged", "obj",
               "hist_blocks", "hist_inner")]
      }
    }
  }

  # The bootstrap restarting the break-point terms declare, run ONCE at the
  # top level -- inside an outer search it would multiply by the number of
  # criterion evaluations -- and at the hyperparameters the fit ended at.
  nb <- seg_boot_total(spec)
  if (nb > 0L && length(attr(design, "refresh"))) {
    res <- statmod_boot_restart(spec, design, blocks, hyper, inner_optimizer,
                                res, expected, approx, maxit, tol, vb, nb)
  }

  coef <- res$obj$split(res$par)
  fitted <- statmod_eta(spec, design, coef)$theta
  # the terms as the fit left them, so a break-point and a nonlinear
  # parameter are read off the fitted object and prediction reapplies them
  spec <- statmod_fitted_spec(spec, coef, design)
  StatmodFit(
    spec = spec, coefficients = coef,
    structural = statmod_structural_par(spec, design), hyper = hyper,
    loglik = statmod_loglik_at(spec, coef, design),
    objective = res$value,
    edf = statmod_edf(spec, coef, design, hyper, expected, approx),
    fitted = fitted, converged = res$converged,
    elapsed = proc.time()[["elapsed"]] - t0,
    criterion = crit,
    history = list(
      outer = res$hist_outer,
      blocks = if (length(res$hist_blocks)) do.call(rbind, res$hist_blocks)
        else NULL,
      inner = if (length(res$hist_inner)) do.call(rbind, res$hist_inner)
        else NULL
    ),
    methods = list(smooth = inner_optimizer, outer = outer_criterion,
                   search = res$optimizer,
                   sparse = vapply(blocks$sparse, function(b)
                     paste(b$param, b$term, sep = "/"), character(1)),
                   # which criterion swept the kinked penalties, and exactly
                   # which of their hyperparameters it reached: a summary
                   # cannot otherwise tell a value the path chose from one
                   # the caller held, and reported both as held
                   sparse_criterion = sparse_criterion,
                   sparse_hyper = if (is.null(sparse_criterion)) character(0)
                     else with(path_rows(spec, blocks, hyper,
                                         sparse_criterion),
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
#' This is a function of its own because the outer search calls it once per
#' hyperparameter it tries, warm-started from the previous coefficients. A
#' whole fit at one point of the search is one call of this.
#'
#' One pass does the smooth block, then each kinked block with everything
#' else held, then commits any refreshable term's schedule. The passes repeat
#' until the objective's relative change falls below `tol`, every block
#' reports it has settled, and [statmod_refresh_settled()] agrees.
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
#' @param approx How the expected information is approximated for a family
#'   with no closed form.
#' @param maxit,tol The budget in passes and the relative tolerance on the
#'   objective, as [method_budget()] read them off `inner_optimizer`.
#' @param vb The resolved verbosity, as [verbosity()] returns it.
#' @param working_budget How many working fits [fit_working()] may
#'   take. The bootstrap excursions of [statmod_boot_restart()]
#'   pass a short one: an excursion needs to travel, not to converge.
#' @param hold_refresh `TRUE` holds every frozen break-point block at
#'   its committed positions: the fit is then an ordinary smooth fit, no
#'   read-off runs and no schedule advances. The outer machinery passes it
#'   at every criterion evaluation, path point and fold, for two reasons
#'   measured together: the working phase inside each of dozens of
#'   evaluations multiplied a fit's cost by twenty, and a break-point
#'   moving between evaluations makes the criterion path-dependent while
#'   its cycling flags read as unavailable points to the search. The
#'   positions are refined once, by the full alternation `statmod()`
#'   runs at the chosen hyperparameters before the restarts.
#'
#' @return A list of six:
#'   \describe{
#'     \item{`par`}{the stacked coefficients reached.}
#'     \item{`value`}{the penalized objective there, unaveraged.}
#'     \item{`converged`}{a single logical: the objective settled, every
#'       block settled, and every refreshable term reported it had.}
#'     \item{`obj`}{the objective object the last pass used.}
#'     \item{`hist_blocks`}{a data frame, one row per pass.}
#'     \item{`hist_inner`}{a data frame, one row per inner iteration.}
#'   }
#'
#' @seealso [statmod()] for the fit this performs,
#'   [fit_smooth()] and [sparse_fit()] for the two halves of a pass,
#'   [statmod_boot_restart()], which calls this repeatedly.
#'
#' @keywords internal
statmod_alternate <- function(spec, design, blocks, hyper, inner_optimizer, beta,
                              expected, approx, maxit, tol, vb,
                              working_budget = 500L, hold_refresh = FALSE) {
  obj <- statmod_objective(spec, hyper, design, expected, approx)
  value <- obj$fn(beta)
  hist_blocks <- list()
  hist_inner <- list()
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
  struct_ok <- TRUE

  # A structural term of the FILTER shape is fitted in the same system as the
  # coefficients rather than alternated with them: the joint gradient and the
  # joint observed information both exist, and the alternation was paying for
  # one optimizer per pass whose every iteration re-ran the recursion. A term
  # of the likelihood shape -- a regime -- keeps the alternation, its
  # information being assembled by a different route.
  joint <- has_structural && length(blocks$smooth) &&
    all(vapply(attr(design, "structural"),
               function(u) identical(u$kind, "filter"), logical(1)))

  for (pass in seq_len(as.integer(maxit))) {
    before <- value

    if (joint) {
      if (vb$blocks) {
        vb_rule(sprintf("inner pass %d: joint system", pass),
                vb_name(inner_optimizer, "newton"), indent = 2L)
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
                               verbose = vb$blocks)
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
                "newton", indent = 2L)
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
                        approx = approx)
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
    # iteration. The alternation is still entered where a term of the
    # LIKELIHOOD shape (a regime) is present, or a sparse block, or a term
    # that recomputes its own design, because those really do alternate.
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
      converged <- converged && isTRUE(struct_ok)
      break
    }
  }
  list(par = beta, value = value, converged = converged, obj = obj,
       note = smooth_note,
       hist_blocks = hist_blocks, hist_inner = hist_inner)
}


#' The Budget and the Stopping Rule of the Alternation
#'
#' @description
#' Reads the iteration budget and the tolerance off the method that fits the
#' smooth block, which is where they are set.
#'
#' @details
#' [statmod()] carries neither a `maxit` nor a `tol` of its own. An argument
#' accepted and ignored is worse than one that signals an error, and a second
#' copy would be exactly that: a caller setting `iwls(maxit = 20)` alongside
#' a loose `maxit = 100` would be obeyed in one and never told about the
#' other.
#' \pkg{distributions7}'s `fit_distrib()` shed the same pair for the same
#' reason.
#'
#' [iwls()] carries both directly. An \pkg{optimizers7} optimizer carries
#' `maxit` and a `criterion`, and the tolerance is taken as the **largest**
#' one in the criterion tree: a combined rule stops at whichever of its parts
#' fires first, so the alternation should not ask for more precision than the
#' loop inside it can deliver.
#'
#' A criterion carrying no tolerance at all leaves the default of
#' [optimizers7::crit_grad()], read from that function's own formals and
#' never copied as a number, so a change there reaches here.
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
                "  The budget and the stopping rule are read off it, so there",
                " is\n  nowhere else for them to come from."), call. = FALSE)
  }
  list(maxit = as.integer(method@maxit),
       tol = criterion_tol(method@criterion))
}


#' What the Inner Method Says About How to Fit
#'
#' @description
#' The information matrix, its approximation and the budget, read off the inner
#' method in one place.
#'
#' @details
#' Every route that fits the coefficients reads these, and reading them in
#' one place keeps a refit inside a path or a fold on the same terms as the
#' fit the caller asked for. Hard-coding them instead ran cross-validation at
#' a tolerance a hundred times tighter than [iwls()]'s own, which cost 26 per
#' cent in time and answered a question the caller had not asked.
#'
#' An \pkg{optimizers7} optimizer says nothing about which information to
#' use, having no such notion, so it gets the expected one and
#' `"bartlett"`.
#'
#' @param method [iwls()] or an \pkg{optimizers7} optimizer.
#'
#' @return A list of four: `expected` (a logical), `approx` (a string),
#'   `maxit` and `tol`.
#'
#' @seealso [method_budget()] for the last two,
#'   [iwls()] for where the first two are set.
#'
#' @keywords internal
inner_settings <- function(method) {
  b <- method_budget(method)
  list(expected = !S7::S7_inherits(method, Iwls) ||
         identical(method@hessian, "expected"),
       approx = if (S7::S7_inherits(method, Iwls)) method@approx else
         "bartlett",
       maxit = b$maxit, tol = b$tol)
}


#' The Tolerance a Criterion Asks For
#'
#' @description
#' Returns the largest tolerance a criterion asks for, descending into a
#' combined rule and taking the maximum over its parts.
#'
#' @details
#' The maximum is the right reduction because a combined rule stops when
#' **any** part fires, so the loop it drives is no more precise than its
#' loosest term. A caller of [method_budget()] wants the precision the loop
#' will actually deliver.
#'
#' @param crit An \pkg{optimizers7} criterion, possibly a combined one built
#'   by `crit_any()` or `crit_all()`.
#'
#' @return A single number. [optimizers7::crit_grad()]'s own default where
#'   the criterion carries no tolerance at all, read from that function's
#'   formals.
#'
#' @seealso [method_budget()], its caller.
#'
#' @keywords internal
criterion_tol <- function(crit) {
  nms <- S7::prop_names(crit)
  if ("criteria" %in% nms) {
    parts <- vapply(crit@criteria, criterion_tol, numeric(1))
    return(max(parts))
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
#' the coefficients outside the block enter it as constants, so the value it
#' reports is the model's objective, never a block's own.
#'
#' With [iwls()] the step is solved through [fit_smooth()]'s own scoring
#' loop, which reads the exact information. With an \pkg{optimizers7}
#' optimizer the objective's `fn`, `gr` and `he` are all supplied, so
#' `newton()` gets the exact second derivative instead of differencing the
#' gradient: measured, that is 2.5x on one smooth and 5.5x on three smooths
#' with a random effect, at identical answers.
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
#' @return A list of five:
#'   \describe{
#'     \item{`par`}{the full stacked coefficient vector, with the block's
#'       positions replaced.}
#'     \item{`value`}{the objective there.}
#'     \item{`converged`}{a single logical, the method's own verdict.}
#'     \item{`iterations`}{how many the method took.}
#'     \item{`history`}{a data frame, one row per iteration.}
#'   }
#'
#' @seealso [statmod_alternate()], the caller,
#'   [sparse_fit()] for the kinked half of the same pass,
#'   [iwls()] for the default method.
#'
#' @keywords internal
fit_smooth <- function(obj, beta, idx, spec, design, hyper, method, vb) {
  whole <- length(idx) == length(beta)
  fn <- function(b) {
    v <- beta
    v[idx] <- b
    obj$fn(v)
  }
  gr <- function(b) {
    v <- beta
    v[idx] <- b
    obj$gr(v)[idx]
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
    sub <- list(fn = fn, gr = gr, npar = length(idx),
                split = obj$split, stack = obj$stack)
    pieces_at <- function(b) {
      v <- beta
      v[idx] <- b
      p <- iwls_pieces(spec, design, obj$split(v), hyper, method)
      if (whole) return(p)
      # a subset of the coefficients takes the corresponding submatrix; the
      # square-root routes need the whole design, so the assembled one is
      # used. This crossprod is the one the profile of the dense lasso path
      # puts at 49 per cent of the whole fit (one X'X per scoring iteration
      # per alternation round per path point), which is why it reads the
      # thread count.
      A <- if (is.null(p$A)) xtx(p$R, spec@threads) + crossprod(p$C) else p$A
      list(R = NULL, C = NULL, A = A[idx, idx, drop = FALSE])
    }
    # the equations' coordinate ranges, restated in the subset's own
    # numbering: the stopping rule's scale is per equation, and the
    # objective's split speaks the full vector's coordinates
    gfull <- obj$split(seq_along(beta))
    groups <- lapply(gfull, function(ix) match(intersect(ix, idx), idx))
    groups <- Filter(length, groups)
    res <- iwls_fit(sub, beta[idx], method, spec@n_obs, pieces_at,
                    verbose = vb$inner, groups = groups)
    out <- beta
    out[idx] <- res$par
    return(list(par = out, value = obj$fn(out), converged = res$converged,
                iterations = res$iterations, history = res$history))
  }

  res <- optimizers7::minimize(method, fn, beta[idx], gr = gr, he = he)
  out <- beta
  out[idx] <- res@par
  list(par = out, value = obj$fn(out), converged = res@converged,
       iterations = res@iterations, history = NULL)
}


#' Fit the Smooth Block Around a Frozen Working Linearization
#'
#' @description
#' The iteration of \cite{fasola2018} for a term whose block is a working
#' linearization with a frozen weight -- [modelterms7::jump()] and
#' [modelterms7::jseg()]: the smooth block is fitted exactly at
#' the committed block, the break-points are read off the fitted
#' coefficients and committed, and the two alternate until the read-off
#' settles or the working objective stops moving.
#'
#' @details
#' The sequencing is the whole of the difference from
#' [fit_smooth()], and it is what `segmented` does. The
#' fixed-point iteration these constructions belong to is not a descent
#' method on the model's objective -- its early steps under a large scaling
#' factor move uphill on purpose, which is how it leaves a spurious optimum
#' -- so embedding the read-off inside the inner optimizer's objective put a
#' sufficient-decrease line search in its way and stalled it: measured on a
#' three-break-point jseg, the embedded route ended at an rss worse than the
#' mean-only fit from the TRUE break-points, while this iteration recovers
#' them from the same start. During the working fit the frozen blocks
#' contribute \eqn{X\beta} and nothing else (`st$working`), which makes
#' the inner fit the plain penalized working fit of the papers; the commit
#' then advances the read-off, the scaling schedule and any relabeling of
#' crossed break-point lineages, once per working fit.
#'
#' Any inner method serves: the working fit goes through
#' [fit_smooth()], which takes [iwls()] or any
#' \pkg{optimizers7} optimizer, and the read-off never moves inside
#' anyone's objective. What differs is the price. Each working fit is
#' solved afresh at a frozen block, so a method carrying exact curvature
#' closes it in a step or two while a quasi-Newton method rebuilds its own
#' from nothing every time: measured on three break-points at
#' \eqn{n = 10000}, the same answer to the digit costs 6.9 s under
#' `iwls()`, 9.2 s under `newton()` and 140 s under
#' `lbfgs()`.
#'
#' The exit is at a fixed point of the iteration or in the cycle it settles
#' into, judged on the working objective: the read-off settled
#' ([modelterms7::term_converged()]) with the objective stalled,
#' the objective stalled three times in a row, or the objective equal to
#' two iterations back twice -- the period-two cycle of the break-point
#' Muggeo documents, which a consecutive-change rule never sees. Only at a
#' fixed point does the working objective coincide with the model's, which
#' is why no best-so-far iterate is kept: mid-travel the committed
#' contribution of a good working value can sit orders of magnitude off
#' the data. Running out of the budget reports `FALSE`.
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
#' @param budget How many working fits at most. The default covers the
#'   measured runs (69 to 165 iterations on three break-points) with room.
#'
#' @return As [fit_smooth()], plus `fasola`, the number of
#'   working fits taken.
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
                        budget = 500L) {
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
    return(res)
  }
  conv <- FALSE
  w_prev <- Inf
  w_prev2 <- Inf
  stall <- 0L
  cyc <- 0L
  it_total <- 0L
  it <- 0L
  for (it in seq_len(budget)) {
    st$working <- TRUE
    st$key <- NULL
    st$value <- NULL
    res <- fit_smooth(obj, beta, idx_a, spec, design, hyper, method, vb)
    st$working <- FALSE
    st$key <- NULL
    st$value <- NULL
    it_total <- it_total + as.integer(res$iterations)
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
  }
  if (vb$blocks) {
    vb_say("%d working fits, %s", it,
           if (conv) "settled" else "budget exhausted", indent = 5L)
  }
  list(par = beta, value = obj$fn(beta), converged = conv,
       iterations = it_total, history = NULL, note = NULL, fasola = it)
}


#' Starting Coefficients
#'
#' @description
#' An intercept-only fit per distribution parameter, with every other
#' coefficient at zero.
#'
#' @details
#' Each equation's intercept starts at the intercept-only MLE, which
#' [distributions7::fit_distrib()] supplies: the model with every
#' covariate removed is the right place for the model with them to begin, and
#' it costs one small fit. The penalized blocks start at zero, where their
#' penalty is smallest.
#'
#' [distributions7::distrib_start()] is the fallback. It returns one
#' list per start, each keyed by parameter, so the value wanted is
#' `th[[1]][[p]]`; indexing the outer list by a parameter's name gives
#' `NULL`, and this function did that, so every start silently fell to
#' zero on the link scale. On a response centered at 5.84 that put the location
#' at 0 and sent the run traveling, which is how a Student t fitted to iris
#' reached a variance of \eqn{10^7}.
#'
#' A start that cannot be obtained is not an error: the fit still runs, from a
#' worse place. What would be an error is not noticing, which is why the two
#' routes are tried in order rather than one being assumed to work.
#'
#' `start` is either a named list of values, a
#' [start_strategy()] --- which is asked once, here, before the
#' alternation between the coefficients and the hyperparameters begins --- or
#' `NULL` for [start_intercepts()], which this function did before
#' strategies existed and still does.
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
      if (length(out[[p]]) != design[[p]]$npar || !is.numeric(out[[p]]) ||
          anyNA(out[[p]])) {
        stop(sprintf(paste0("the start strategy returned %d values for '%s'",
                            ", which has %d coefficients."),
                     length(out[[p]]), p, design[[p]]$npar), call. = FALSE)
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
#' The intercept-only maximum likelihood estimate, where it can be had, and a
#' draw from the parameter's domain otherwise.
#'
#' @details
#' Two routes, tried in order. [distributions7::fit_distrib()] fits
#' the distribution to the response with no covariates, which is the same model
#' with every slope set to zero and therefore exactly where the fit should
#' begin; its link-scale coefficients are the intercepts.
#' [distributions7::distrib_start()] is the fallback, and its result
#' is a list of starts, each keyed by parameter, so a value is reached at
#' `[[1]][[p]]`.
#'
#' **The random stream is pinned and restored.** That intercept-only fit
#' starts from draws over the parameters' domains, so it returns a different
#' answer on every call where a parameter is weakly identified -- fitted to
#' `iris`, a Student t's \eqn{\nu} came back at \eqn{e^{39}}, \eqn{e^{21}}
#' and \eqn{e^{17}} on three consecutive runs, and `statmod()` inherited
#' that: the same call gave log-likelihoods of -103.49, -112.11 and -111.83. A
#' fitting function has to give the same answer twice, so the seed is fixed for
#' the length of this call and the caller's stream is put back afterwards.
#'
#' Pinning makes it reproducible without making it good, one draw being
#' one draw; several are taken and the best kept. What would make it good is a
#' data-based `distrib_start` method on the univariate families, which is
#' the design \pkg{distributions7} already documents and which only its
#' multivariate gaussian implements.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A named list, one entry per distribution parameter, on the link
#'   scale; an entry is `NULL` where neither route answered.
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
      return(stats::setNames(as.list(as.numeric(e)), params))
    }
  }

  th <- tryCatch(
    distributions7::distrib_start(spec@distrib, spec@response, 1L),
    error = function(e) NULL)
  if (is.null(th) || !length(th)) return(none)
  th1 <- th[[1L]]
  stats::setNames(lapply(params, function(p) {
    v <- th1[[p]]
    if (is.null(v) || !is.finite(v[[1L]])) return(NULL)
    linkfunctions7::linkfun(links[[p]], v[[1L]])
  }), params)
}


#' Effective Degrees of Freedom, Per Term
#'
#' @description
#' Asks each term what it spends, through \pkg{modelterms7}'s `edf()`.
#'
#' @details
#' A smooth penalized term counts \eqn{\mathrm{tr}[(H+S)^{-1}H]} over its own
#' block, so it needs the unpenalized curvature there, never only its
#' coefficients and its hyperparameters. That block is cut out of the
#' likelihood's information, which is computed once for every term rather than
#' per term.
#'
#' A term may carry more than one penalty, over different parameters of its
#' own, and each has hyperparameters of its own filed under a key of its own.
#' The row stays per term, which is the granularity a table of terms wants,
#' and the hyperparameters are handed over keyed by the penalty names
#' [modelterms7::term_penalties()] gives, which is the shape
#' `edf()` reads them in. Passing the hyperparameters of one penalty for
#' a term carrying two would count the whole block against it.
#'
#' The arguments are passed by name. `edf()`'s third argument is the
#' curvature and its fourth the hyperparameters, and a positional call put the
#' hyperparameters where the curvature belongs: every smooth term then reported
#' `NA`, the total degrees of freedom counted the unpenalized terms alone,
#' and AIC and BIC were built on that count.
#'
#' @param spec The specification.
#' @param coef The coefficients.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param expected Whether the curvature is the expected information.
#' @param approx The approximation for the expected information.
#'
#' @return A data frame with one row per term.
#'
#' @keywords internal
statmod_edf <- function(spec, coef, design, hyper, expected = TRUE,
                        approx = "opg") {
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
  smoother <- NULL
  if (any(vapply(statmod_penalized(spec, design),
                 function(u) !penalty_has_kink(u$penalty, u$key),
                 logical(1)))) {
    smoother <- tryCatch({
      H <- statmod_information_at(spec, coef, design, expected, approx)
      S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
      S <- zap_nonfinite(S)
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
        zn <- modelterms7::term_params(spec@terms[[p]][[nm]])
        st <- statmod_structural_state(design)
        as.numeric(length(setdiff(zn, st$held[[nm]])))
      } else if (!length(ent) && is.null(smoother)) {
        as.numeric(length(cols))
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
          H <- statmod_information_at(spec, coef, design, expected, approx)
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


#' Resolve the Verbosity Setting
#'
#' @description
#' Turns a level or a named logical vector into the three switches the fit
#' reads.
#'
#' @details
#' The levels name the loops rather than counting them: `1` shows the
#' alternation, `2` the inner method's own iterations, `3` the
#' optimizers' traces. The named form exists because the three are genuinely
#' independent -- watching the alternation while silencing a chatty inner
#' optimizer is the common case.
#'
#' @param verbose A number from 0 to 3, or a named logical vector with any of
#'   `blocks`, `inner`, `optimizer`.
#'
#' @return A list of three logicals.
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
