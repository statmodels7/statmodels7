#' @include augmented.R solver.R
NULL

#' The Iterated Weighted Least Squares Method
#'
#' @description
#' `Iwls` is the S7 class of the scoring step that [statmod()] uses for the
#' jointly fitted smooth block. An object of this class carries the solver,
#' the budget, the stopping rule and the constants of the loop. Build one with
#' [iwls()], which validates every argument; this raw constructor does not.
#'
#' The matrix that the step inverts is not an argument of [iwls()]: it is the
#' argument `inner_information` of [statmod()]. [statmod()] stores the matrix
#' that it settled in the properties `hessian`, `backup`, `switching` and
#' `repair_floor`, which [iwls()] initializes to an unsettled state.
#'
#' @param solver How the step is solved, an object that inherits from
#'   [IwlsSolver]: [solve_qr()], [solve_svd()] or [solve_chol()].
#' @param maxit The iteration budget, a single positive whole number.
#' @param tol The stopping tolerance on the score per observation, a single
#'   positive number.
#' @param criterion An \pkg{optimizers7} `criterion` object driving the loop
#'   in place of `tol`, or `NULL` for the built-in rule.
#' @param step_halving The number of halvings tried after the full step before
#'   a step is abandoned, a single non-negative whole number. The line search
#'   tries the full step and then at most `step_halving` halved steps, so `0`
#'   tries only the full step.
#' @param decrement_limit The value under which the built-in rule reads the
#'   Newton decrement as small, a single positive number.
#' @param switch_after A single number: the iteration after which a run on
#'   the expected information continues on the observed one when the choice
#'   was settled from [auto_info()]. `Inf` never switches.
#' @param switch_rate A single number, or `NA`: the ratio of consecutive scores
#'   above which a run on the expected information counts as contracting
#'   slowly. `NA` never counts a run as slow.
#' @param hessian The information that the step inverts, once settled: an
#'   object that inherits from [distributions7::Information], set by
#'   [statmod()] from its argument `inner_information`. An [auto_info()] object
#'   stands for an unsettled choice, which `iwls_resolve()` settles against the
#'   family.
#' @param backup `NULL`, or the [distributions7::expected_info()] object that
#'   takes the step where the observed step finds no acceptable point. It is
#'   set only by settling an [auto_info()] on the observed information.
#' @param switching Whether the run continues on the observed information after
#'   `switch_after` iterations or a slow contraction. It is set only by
#'   settling an [auto_info()] on the expected information.
#' @param repair_floor The `repair_floor` of the [auto_info()] that was
#'   settled, or `NA` when the choice was named.
#'
#' @return An object of class `Iwls` with one property per argument above, and
#'   each property holds the value that was passed.
#'
#' @seealso [iwls()], the constructor to use.
#'
#' @examples
#' S7::S7_inherits(iwls(), Iwls)
#' iwls()@solver
#'
#' @name Iwls-class
#' @aliases Iwls
#' @keywords internal
#' @export
Iwls <- S7::new_class("Iwls",
  properties = list(
    solver = IwlsSolver,
    maxit = S7::class_numeric,
    tol = S7::class_numeric,
    criterion = S7::class_any,
    step_halving = S7::class_numeric,
    decrement_limit = S7::class_numeric,
    switch_after = S7::class_numeric,
    switch_rate = S7::class_numeric,
    hessian = distributions7::Information,
    backup = S7::class_any,
    switching = S7::class_logical,
    repair_floor = S7::class_numeric
  )
)


#' Iterated Weighted Least Squares
#'
#' @description
#' Configures the scoring step with which `statmod()` fits the smooth block,
#' and returns it as an object to pass as `inner_optimizer`. This is the default,
#' so `statmod(...)` and `statmod(..., inner_optimizer = iwls())` fit the
#' same model; call it by name to change the solver, the budget, the stopping
#' rule or the constants of the loop.
#'
#' The matrix that the step inverts is chosen by the argument
#' `inner_information` of [statmod()] (see [information_objects]), and not by
#' `iwls()`: the step is a Newton method on whatever curvature the objective
#' offers. Every argument of `iwls()` has a default, and the object carries no
#' data, so one `iwls()` can be passed to any number of fits.
#'
#' @details
#' # The name
#'
#' The name abbreviates iterated weighted least squares. The word iterated
#' already implies the re-weighting, so the abbreviation `irls` is not used.
#'
#' # The curvature
#'
#' The step inverts the matrix that `inner_information` of [statmod()]
#' selects. With [distributions7::expected_info()] it inverts the expected
#' information, which is Fisher scoring, and with
#' [distributions7::observed_info()] the observed Hessian, which is Newton.
#' The expected information is positive definite wherever the family is
#' regular, so Fisher scoring takes usable steps from a poor start where
#' Newton's matrix may be indefinite; Newton converges faster near the
#' optimum. An explicit choice is the only matrix that the step uses: the
#' observed information is never replaced by the expected one.
#'
#' The default `inner_information = auto_info()` is settled against the family
#' when the fit starts. The expected information is used where the family
#' writes it exactly and at an acceptable cost. Otherwise the observed
#' information is used, and the expected information takes the step where the
#' observed step, after the repair of an indefinite matrix, finds no
#' acceptable point. A trial point whose objective raises an error is read as
#' rejected. The fit records the settled choice, so `fit@methods$smooth` shows
#' which curvature ran.
#'
#' On a family that [auto_info()] settles on the expected information, a run
#' continues on the observed information as soon as Fisher scoring contracts
#' the score slowly, which means that each of two consecutive ratios of the
#' score is above `switch_rate`, and at the latest after `switch_after`
#' scoring steps. Fisher scoring converges only linearly near the mode
#' wherever the two informations differ there, which they do whenever the
#' dispersion has an equation of its own, whereas Newton's step converges
#' quadratically.
#'
#' Where the stopping rule is met, the run takes one more full step. The step
#' uses the observed information if the choice was settled on the expected
#' information from [auto_info()] and has not yet switched, and otherwise the
#' information that the run is using. The marginal criterion carries
#' \eqn{\log\lvert K(\beta)\rvert}, which is not stationary in \eqn{\beta},
#' so a mode left short by \eqn{\delta} enters it at first order; a
#' warm-started fit whose score is already under the rule would stop
#' without moving and leave the criterion out by an amount linear in the
#' outer step. The extra step takes the error to second order. It is kept
#' unless it raises the objective by more than a relative \eqn{10^{-8}}: a
#' sufficient-decrease test at a point already under the rule would accept or
#' reject by the last bits, and the criterion would jump between nearby
#' hyperparameters.
#'
#' The `method` carried by [distributions7::expected_info()] and [auto_info()]
#' is passed to \pkg{distributions7}. With `method = NULL`, the default, the
#' closed form of the family is used where the family has one
#' ([distributions7::expected_hessian_exact()] records which case a family is
#' in), and [distributions7::expect_opg()], the outer product of the observed
#' scores, which costs one gradient, otherwise. An object is used for every
#' family, also for one that has the closed form. [distributions7::expect_integrate()]
#' and [distributions7::expect_mc()] evaluate the expectation itself (a sum
#' over the support or a quadrature, or simulated draws) and cost much more
#' inside a loop that rebuilds the curvature at every iteration.
#' [distributions7::expected_derivative_methods] compares the three. The price
#' of the cheap matrix is paid in variance and not in the answer: the score is
#' exact, so a scoring step driven by any positive definite matrix reaches the
#' same stationary point, and the outer product is positive semidefinite by
#' construction, whereas the observed Hessian need not be. The information
#' behind a standard error is chosen separately; see [vcov_method].
#'
#' # The solver
#'
#' The step solves \eqn{(X'WX + S)\delta = X'g - S\beta}, and the argument
#' `solver` selects the method.
#'
#' \describe{
#'   \item{[solve_qr()]}{a QR of the augmented design
#'     \eqn{[\sqrt{W}X;\ \mathrm{chol}(S)]}, and the default. Forming
#'     \eqn{X'X} squares the condition number, and the design of a
#'     break-point term is bounded only to \eqn{\epsilon^{-1/2}} by
#'     construction, so that squaring it exhausts double precision.}
#'   \item{[solve_svd()]}{a singular value decomposition of the same augmented
#'     design. It reports a numerical rank where a deficient block would make
#'     a Cholesky fail.}
#'   \item{[solve_chol()]}{a Cholesky factor of the penalized information
#'     \eqn{X'WX + S}, assembled without the augmentation. The information is
#'     formed explicitly, so the condition number is squared.}
#' }
#'
#' [solve_qr()] and [solve_svd()] need a square root of the per-observation
#' curvature and of the penalty. Where either does not exist, the step uses
#' the assembled information, as [solve_chol()] does. [iwls_solvers] compares
#' the three and gives the tolerances that each one takes.
#'
#' # The stopping rule
#'
#' A scoring step is not an optimizer and carries its own loop, but the rule
#' that ends the loop belongs to the caller.
#'
#' With `criterion = NULL`, the default, the built-in rule applies. The rule
#' is met when the score per observation \eqn{\max_j \lvert g_j \rvert / n},
#' taken over the free coordinates, is below `tol` and the Newton decrement
#' \eqn{g^\top (H+S)^{-1} g / 2} is at most `decrement_limit`. The decrement is
#' in log-likelihood units and does not depend on the scale of the response.
#' When the loop ends without meeting the rule, the dimensionless reading of
#' [iwls_score()] can relabel the stopped run as converged. Any
#' \pkg{optimizers7} criterion may drive the loop instead. It replaces the
#' whole built-in rule, including the decrement condition and the
#' dimensionless reading, and `tol` is then not read, so passing both is an
#' error.
#'
#' The state that a criterion receives is the one that
#' `optimizers7::crit_met()` documents. In this loop the state has three
#' properties:
#'
#' - Its `gradient` is the score per observation, the quantity that the
#'   built-in rule compares with `tol`. The threshold of
#'   `optimizers7::crit_grad(t)` therefore applies to the same quantity as
#'   `tol = t` and has the same meaning at \eqn{n = 10} and at
#'   \eqn{n = 10^7}.
#' - Its `f_new` is the penalized objective, the negative penalized
#'   log-likelihood, not divided by the sample size, which is the scale on
#'   which the penalty is added. A rule that reads the absolute value of the
#'   objective therefore depends on the sample size; a rule that reads its
#'   relative change does not.
#' - On the first iteration `f_old` and `x_old` are `NULL`, because there is
#'   no previous point, so a rule reading a change returns `FALSE` there by
#'   construction.
#'
#' A criterion that needs a quantity that the step does not compute is
#' rejected at construction. A stationarity measure is an example: it belongs
#' to the derivative-free methods, and a rule that requires one would remain
#' in the loop without ever being met.
#'
#' # The constants of the loop
#'
#' Three constants of the loop are arguments. `decrement_limit` is the value
#' under which the built-in rule reads the Newton decrement as small
#' (default `1e-6`, three orders under the reading that the outer search uses
#' to decide whether a criterion is usable). `switch_after` is the iteration
#' after which a run settled on the expected information continues on the
#' observed one (default `10`; `Inf` never switches). `switch_rate` is the
#' ratio of consecutive scores above which two consecutive steps count as a
#' slow contraction (default `0.5`; `NA` never counts one). The last two are
#' read only when `inner_information` is settled from [auto_info()] on the
#' expected information.
#'
#' @param solver The factorization with which each step is solved, an object
#'   that inherits from [IwlsSolver]: [solve_qr()] (the default), [solve_svd()]
#'   or [solve_chol()]. [iwls_solvers] compares the three and gives their
#'   tolerances.
#' @param decrement_limit The value under which the Newton decrement is read as
#'   small by the built-in stopping rule, a single positive number. The
#'   default is `1e-6`.
#' @param switch_after The iteration after which a run settled on the expected
#'   information from [auto_info()] continues on the observed information, a
#'   single positive number; `Inf` never switches. The default is `10`.
#' @param switch_rate The ratio of consecutive scores above which two
#'   consecutive steps count as a slow contraction of Fisher scoring, a single
#'   number or `NA`, which never counts a run as slow. The default is `0.5`.
#' @inheritParams Iwls-class
#'
#' @return An [Iwls()] object holding the settings given and no data. Its
#'   `hessian` is an unsettled [auto_info()] until [statmod()] settles the
#'   argument `inner_information`.
#'
#' @seealso [statmod()] for the argument `inner_information` and for where the
#'   object is passed, [information_objects] and [iwls_solvers] for the matrix
#'   and the solver, [iwls_score()] for the built-in rule's final verdict,
#'   [optimizers7::crit_grad()] for the criteria that can replace it.
#'
#' @examples
#' iwls()
#'
#' # This call solves through an SVD, which reports a numerical rank where a
#' # deficient block would make a Cholesky fail.
#' iwls(solver = solve_svd())
#'
#' # This call uses an optimizers7 stopping rule in place of tol. The
#' # gradient in the state that the rule receives is the score per
#' # observation, as for tol.
#' iwls(criterion = optimizers7::crit_any(optimizers7::crit_grad(1e-8),
#'                                        optimizers7::crit_rel_obj(1e-12)))
#'
#' # The matrix is chosen in statmod(): this is Newton, on the observed
#' # information.
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(solver = solve_svd()),
#'                inner_information = distributions7::observed_info())
#' fit@methods$smooth
#'
#' # Both at once is an error: tol would not be read.
#' try(iwls(tol = 1e-8, criterion = optimizers7::crit_grad(1e-8)))
#'
#' @export
iwls <- function(solver = solve_qr(), maxit = 1000L, tol = 1e-6,
                 criterion = NULL, step_halving = 30L,
                 decrement_limit = 1e-6, switch_after = 10, switch_rate = 0.5) {
  check_solver(solver)
  if (!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) ||
      maxit < 1 || maxit != round(maxit)) {
    stop("'maxit' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(step_halving) || length(step_halving) != 1L ||
      !is.finite(step_halving) || step_halving < 0 ||
      step_halving != round(step_halving)) {
    stop("'step_halving' must be a single non-negative integer.",
         call. = FALSE)
  }
  check_positive_number(decrement_limit, "decrement_limit")
  if (!is.numeric(switch_after) || length(switch_after) != 1L ||
      is.na(switch_after) || switch_after <= 0) {
    stop("'switch_after' must be a single positive number or Inf.",
         call. = FALSE)
  }
  if (length(switch_rate) != 1L ||
      !(is.na(switch_rate) || (is.numeric(switch_rate) && is.finite(switch_rate)))) {
    stop("'switch_rate' must be a single number or NA.", call. = FALSE)
  }
  if (!is.null(criterion)) {
    if (!S7::S7_inherits(criterion, optimizers7::criterion)) {
      stop("'criterion' must be an optimizers7 criterion, or NULL.",
           call. = FALSE)
    }
    # the step computes a gradient and nothing else; a rule asking for more
    # would never fire, which is what optimizers7's own check_criterion()
    # refuses at construction for the same reason
    missing_ <- setdiff(optimizers7::crit_needs(criterion), "gradient")
    if (length(missing_)) {
      stop(sprintf(paste0("The stopping rule needs %s, which a scoring step ",
                          "does not compute.\n  Choose a rule that it can ",
                          "evaluate, or leave 'criterion' NULL."),
                   paste(missing_, collapse = ", ")), call. = FALSE)
    }
    if (!missing(tol)) {
      stop(paste0("'tol' and 'criterion' both set the stopping rule: pass one.\n",
                  "  'tol' belongs to the built-in rule, which a criterion ",
                  "replaces."), call. = FALSE)
    }
  }
  Iwls(solver = solver,
       maxit = as.numeric(maxit), tol = tol, criterion = criterion,
       step_halving = as.numeric(step_halving),
       decrement_limit = as.numeric(decrement_limit),
       switch_after = as.numeric(switch_after),
       switch_rate = as.numeric(switch_rate),
       hessian = auto_info(), backup = NULL, switching = FALSE,
       repair_floor = NA_real_)
}

#' @export
#' @param x An [Iwls()] object.
#' @param ... Unused.
#' @rdname iwls
print.Iwls <- function(x, ...) {
  settled <- !info_auto(x@hessian)
  cat(sprintf("iwls: %s%s\n", solver_how(x@solver),
              if (!settled) ""
              else sprintf(", %s information%s", info_label(x@hessian),
                           if (!is.null(x@backup))
                             " (the expected where it cannot step)"
                           else if (isTRUE(x@switching))
                             sprintf(" (the observed after %s iterations)",
                                     format(x@switch_after))
                           else "")))
  cat(sprintf("  maxit %d, %s\n", as.integer(x@maxit),
              if (is.null(x@criterion)) sprintf("tol %g", x@tol)
              else x@criterion@label))
  invisible(x)
}
S7::method(print, Iwls) <- print.Iwls


#' Settle the Curvature of an `iwls()` Against a Family
#'
#' @description
#' Turns the information of an [iwls()] into the matrix that a fit runs. With
#' `information` given, that object replaces the one that the method holds;
#' without it, a method that holds an unsettled [auto_info()] is settled and a
#' method whose curvature was named is returned as it is.
#'
#' @details
#' A named [distributions7::observed_info()] or [distributions7::expected_info()]
#' is the matrix of the step and nothing else: the method carries no
#' fallback and no switch.
#'
#' An [auto_info()] is settled by two questions.
#' [distributions7::expected_hessian_exact()] is the first. Where it returns
#' `FALSE` the expected information is an approximation, the outer product of
#' the scores by default, and the method is settled on the observed
#' information, with the expected information as the `backup` that takes the
#' step where the observed step finds no acceptable point. An observed
#' penalized information that is not positive definite is first floored by
#' [pd_repair()] and its step tried.
#'
#' [distributions7::expected_hessian_costly()] is the second question. A
#' family whose expected information is exact but costly is settled on the
#' observed information as well, because the expected information of such a
#' family requires a quadrature or a sum over the support at each
#' observation. The expected information is the `backup` there too.
#'
#' Only a family settled on the expected information has `switching` set.
#'
#' @param method An [Iwls()] object.
#' @param distrib The distribution to which the model is fitted.
#' @param information `NULL`, or an object that inherits from
#'   [distributions7::Information] that replaces the one of `method`.
#'
#' @return An [Iwls()] object whose `hessian` is an
#'   [distributions7::expected_info()] or an [distributions7::observed_info()]
#'   object, unless a named choice is returned as it is.
#'
#' @seealso [iwls()], [fit_smooth()], which reads `backup`.
#'
#' @keywords internal
iwls_resolve <- function(method, distrib, information = NULL) {
  if (is.null(information)) {
    if (!info_auto(method@hessian)) return(method)
    information <- method@hessian
  }
  method@backup <- NULL
  method@switching <- FALSE
  method@repair_floor <- NA_real_
  if (!info_auto(information)) {
    method@hessian <- information
    return(method)
  }
  exact <- distributions7::expected_hessian_exact(distrib) &&
    !distributions7::expected_hessian_costly(distrib)
  expd <- distributions7::expected_info(method = information@method)
  method@hessian <- if (exact) expd else distributions7::observed_info()
  method@backup <- if (exact) NULL else expd
  method@switching <- exact
  method@repair_floor <- information@repair_floor
  method
}


#' Whether Fisher Scoring Is Contracting Slowly
#'
#' @description
#' Returns `TRUE` where the last two ratios of consecutive scores both exceed
#' `rate`.
#'
#' @param scores The scores of the run so far, one per iteration.
#' @param rate The ratio above which a step counts as slow, or `NA`, for which
#'   the function never returns `TRUE`.
#'
#' @return A single logical.
#'
#' @seealso [iwls_fit()], [iwls()], whose argument `switch_rate` supplies
#'   `rate`.
#'
#' @keywords internal
iwls_slow <- function(scores, rate) {
  k <- length(scores)
  if (is.na(rate) || k < 3L) return(FALSE)
  r <- scores[(k - 1L):k] / scores[(k - 2L):(k - 1L)]
  all(is.finite(r)) && all(r > rate)
}


#' The Relative Level Under Which a Penalized Information Is Indefinite
#'
#' @description
#' Returns the relative floor with which [iwls_fit()] tests whether the
#' penalized information of the step is positive definite: the `repair_floor`
#' of the [auto_info()] that was settled, and the `floor` of the solver when the
#' curvature was named.
#'
#' @param method An [Iwls()] object.
#'
#' @return A single positive number.
#'
#' @seealso [pieces_definite()], [auto_info()]
#'
#' @keywords internal
iwls_definite_tol <- function(method) {
  rf <- method@repair_floor
  if (length(rf) && is.finite(rf)) rf else solver_floor(method@solver)
}


#' Whether an Inner Fit Takes One Step Past Its Stopping Rule
#'
#' @description
#' Returns `TRUE` where [iwls_fit()] takes a full step at the point that its
#' rule accepted, so a criterion read there is out by a second-order amount.
#'
#' @return A single logical.
#'
#' @seealso [iwls_fit()].
#'
#' @keywords internal
iwls_polish <- function() TRUE


#' Solve One Weighted Least Squares Step
#'
#' @description
#' Returns the increment \eqn{\delta} solving \eqn{(H + S)\delta = u} by the
#' requested decomposition.
#'
#' @details
#' # The two families of route
#'
#' `"qr"` and `"svd"` decompose the augmented matrix \eqn{[R;\, C]}, whose
#' cross-product is the penalized information, so \eqn{X'X} is never formed
#' and the condition number is never squared. That route needs the
#' per-observation curvature to be positive definite and the penalty positive
#' semidefinite. Where either fails, as with an observed Hessian far from the
#' optimum or a non-convex penalty, there is no square root to take: the step
#' falls back to the assembled matrix with its eigenvalues floored, and
#' records the fallback in `route` as `"qr->chol"` or `"svd->chol"`.
#'
#' `"chol"` factors the assembled information directly, after flooring the
#' eigenvalues of a matrix that is not positive definite.
#'
#' # Levenberg damping
#'
#' `damp` adds \eqn{\lambda I} to the system, which shortens a coordinate in
#' proportion to how little curvature it has. A coordinate whose diagonal is
#' \eqn{d} is shortened by \eqn{(d + \lambda)/d}, so that a coordinate with a
#' small diagonal is shortened more than one with a large diagonal: a
#' differential shrink, which a scalar step length cannot give. The damping is
#' zero unless a step has failed.
#'
#' On the augmented route the damping enters as further rows of the augmented
#' design, so the QR route keeps its conditioning; on the assembled route it is
#' added to the diagonal. The matrix added is the identity. A damping
#' proportional to \eqn{\mathrm{diag}(K)} would shrink every coordinate
#' alike and leave the disparity exactly where it was.
#'
#' # A coordinate at a boundary is held
#'
#' Where a parameter has reached the clamp with which its link keeps it
#' strictly inside its bounds, the family's curvature there is zero or `NaN`,
#' and that coordinate cannot be moved by any step. Those coordinates are
#' dropped from the system and the rest is solved, which is the active set
#' that the boundary defines; the step stays a descent step for the reduced
#' problem. Coordinates in `frozen` and coordinates with a non-finite
#' right-hand side are held in the same way.
#'
#' Holding them keeps a single boundary coordinate from preventing the others
#' from moving: a `NaN` in one entry of the curvature would otherwise make the
#' whole step zero.
#'
#' On the augmented route a coordinate is held when its column of the
#' square-root design or of the penalty factor has a non-finite entry. On the
#' assembled route the test is on the diagonal of the curvature. A boundary
#' coordinate makes its entire row non-finite, cross terms included, so a test
#' over whole columns would mark its neighbors too and hold coordinates that
#' can move. After the diagonal test, the columns of the reduced matrix are
#' tested again, and the test repeats until every column is finite.
#'
#' @param pieces A list with `R`, `C` and `A`, as [iwls_pieces()] builds it.
#'   Which of the three are present depends on `how`.
#' @param u The right-hand side \eqn{X'g - S\beta}, a numeric vector of
#'   length `p`.
#' @param how The decomposition: `"qr"`, `"svd"` or `"chol"`.
#' @param damp The Levenberg damping \eqn{\lambda}, a single non-negative
#'   number. `0`, the default, is the plain scoring step.
#' @param frozen Positions that the caller holds at their current values, an
#'   integer vector. They are dropped from the system exactly as a
#'   coordinate whose curvature is not finite is, so the step is that of the
#'   problem reduced to the rest. The default is empty.
#' @param tol The tolerance of the factorization, as the solver carries it:
#'   the pivoted QR tolerance for `how = "qr"`, the relative cut of the
#'   singular values for `how = "svd"` (`NULL` for the default cut), and
#'   unused for `"chol"`. The default is `NULL`.
#' @param floor The relative floor of the eigenvalues of an assembled matrix
#'   that is not positive definite, as [pd_repair()] takes it. The default is
#'   `1e-8`.
#'
#' @return A list of five:
#'   \describe{
#'     \item{`delta`}{the increment, a numeric vector of length `p`, zero in
#'       every held coordinate.}
#'     \item{`rank`}{the rank of the system solved, an integer.}
#'     \item{`route`}{which decomposition was used: the one requested, or
#'       `"qr->chol"` or `"svd->chol"` when the augmented route fell back to
#'       the assembled matrix with its eigenvalues floored.}
#'     \item{`held`}{the positions dropped from the system, an integer
#'       vector, empty when every coordinate was free.}
#'     \item{`dropped`}{the positions that the pivoted `qr()` left out as
#'       aliased, an integer vector, empty on every other route.}
#'   }
#'
#' @seealso [iwls_pieces()] for the input, [iwls()] for the `solver`
#'   argument that this function serves.
#'
#' @keywords internal
iwls_solve <- function(pieces, u, how, damp = 0, frozen = integer(0),
                       tol = NULL, floor = 1e-8) {
  p <- length(u)
  # A COORDINATE THE CALLER HOLDS is dropped from the solve exactly as a
  # coordinate whose curvature is gone: the step is that of the REDUCED
  # problem, which is what a constrained maximization means. Zeroing the
  # increment after an unconstrained solve would not do -- that gives the
  # free components of the full step, whose fixed point is not where the
  # free scores vanish.
  froz <- logical(p)
  if (length(frozen)) froz[frozen] <- TRUE
  # sqrt(lambda) I as further rows of the augmented design, which is how a
  # ridge is added there: [R; C; sqrt(lambda) I] has cross-product
  # K + lambda I and the route keeps its conditioning.
  damp_rows <- function(C, k) {
    if (!(damp > 0)) return(C)
    D <- if (isS4(C)) Matrix::Diagonal(k, sqrt(damp)) else diag(sqrt(damp), k)
    rbind(C, D)
  }
  # A COORDINATE WHOSE CURVATURE IS NOT FINITE IS ONE COORDINATE, NOT THE
  # WHOLE POINT. It happens where a parameter has reached the clamp its link
  # keeps it strictly inside: the family's derivatives there are 0 or NaN,
  # and the coordinate cannot move whatever is solved. Until 0.82.0 a single
  # such entry made the whole system unusable and the step came back zero for
  # EVERY coordinate, so a run stopped with the others far from stationary --
  # measured on a Student t whose nu reached double.xmax, the score in sigma
  # was -617.6 there while nu's was exactly 0. Those coordinates are held and
  # the rest is solved, which is the active set a boundary defines: the step
  # is a descent step for the reduced problem, and the held coordinate's own
  # score being zero is what says it is a stationary point of the constrained
  # problem rather than a place the fit was stopped at.
  # `colSums(abs(.))` and not `is.finite()` over the whole matrix: a column
  # sum keeps a sparse design sparse where an elementwise test would build a
  # dense logical of its size, and `abs` cannot overflow into a false
  # positive where a square could. A NaN or an infinity anywhere in the
  # column reaches the sum.
  bad_cols <- function(M) !is.finite(as.numeric(Matrix::colSums(abs(M))))
  held <- integer(0)
  if (how %in% c("qr", "svd") && !is.null(pieces$R) && !is.null(pieces$C)) {
    bad <- froz | !is.finite(u) | bad_cols(pieces$R) | bad_cols(pieces$C)
    if (any(bad)) {
      held <- which(bad)
      keep <- which(!bad)
      if (!length(keep)) {
        return(list(delta = numeric(p), rank = 0L, route = how, held = held,
                    dropped = integer(0)))
      }
      out <- augmented_lowrank(pieces$R[, keep, drop = FALSE],
                               damp_rows(pieces$C[, keep, drop = FALSE],
                                         length(keep)),
                               if (is.null(pieces$V)) NULL else
                                 pieces$V[keep, , drop = FALSE],
                               u[keep], how,
                               threads = if (is.null(pieces$threads)) 1L
                                         else pieces$threads, tol = tol)
      d <- numeric(p)
      d[keep] <- out$delta
      # the solve ran on the kept columns, so its coordinates are read back
      # through them before they mean anything to a caller
      return(list(delta = d, rank = out$rank, route = how, held = held,
                  dropped = keep[out$dropped]))
    }
    out <- augmented_lowrank(pieces$R, damp_rows(pieces$C, p), pieces$V, u,
                             how, threads = if (is.null(pieces$threads)) 1L
                                           else pieces$threads, tol = tol)
    out$route <- how
    out$held <- held
    return(out)
  }
  route <- if (how %in% c("qr", "svd")) paste0(how, "->chol") else how
  A0 <- pieces$A
  if (!is.null(A0)) {
    # BY THE DIAGONAL AND NOT BY THE COLUMN. A boundary coordinate makes its
    # whole ROW non-finite, cross terms included, so a column test marks its
    # neighbors too: measured on the Student t, testing columns held sigma
    # along with nu and left the fit exactly where it had been. The diagonal
    # names the coordinate whose own curvature is gone; dropping its row and
    # column then clears the cross terms, and the loop repeats in case the
    # reduced matrix still carries one from another source.
    bad <- froz | !is.finite(u) | !is.finite(as.numeric(Matrix::diag(A0)))
    repeat {
      keep <- which(!bad)
      if (!length(keep)) break
      more <- bad_cols(A0[keep, keep, drop = FALSE])
      if (!any(more)) break
      bad[keep[more]] <- TRUE
    }
    if (any(bad)) {
      held <- which(bad)
      keep <- which(!bad)
      if (!length(keep)) {
        return(list(delta = numeric(p), rank = 0L, route = route, held = held,
                    dropped = integer(0)))
      }
      Vk <- if (is.null(pieces$V)) NULL else
        pieces$V[keep, , drop = FALSE]
      lr <- lowrank_solve(A0[keep, keep, drop = FALSE], Vk, u[keep], damp)
      if (!is.null(lr)) {
        d <- numeric(p)
        d[keep] <- lr
        return(list(delta = d, rank = length(keep), route = route,
                    held = held, dropped = integer(0)))
      }
      A <- pd_repair(lowrank_dense(A0[keep, keep, drop = FALSE], Vk),
                     rel = floor)
      if (is.null(A)) {
        return(list(delta = numeric(p), rank = 0L, route = route, held = held,
                    dropped = integer(0)))
      }
      if (damp > 0) Matrix::diag(A) <- Matrix::diag(A) + damp
      ch <- chol(A)
      d <- numeric(p)
      d[keep] <- as.numeric(backsolve(ch, forwardsolve(t(ch), u[keep])))
      return(list(delta = d, rank = ncol(A), route = route, held = held,
                  dropped = integer(0)))
    }
  }
  if (!is.null(A0) && !is.null(pieces$V)) {
    lr <- lowrank_solve(A0, pieces$V, u, damp)
    if (!is.null(lr)) {
      return(list(delta = lr, rank = length(u), route = route, held = held,
                  dropped = integer(0)))
    }
    A0 <- lowrank_dense(A0, pieces$V)
  }
  A <- pd_repair(A0, rel = floor)
  # nothing can be solved at a point whose curvature is not finite ANYWHERE;
  # the caller reads a zero step as "this iterate is unusable" and keeps the
  # last good one
  if (is.null(A)) {
    return(list(delta = numeric(p), rank = 0L, route = route, held = held,
                dropped = integer(0)))
  }
  if (damp > 0) Matrix::diag(A) <- Matrix::diag(A) + damp
  ch <- chol(A)
  # a factorization either succeeds or is repaired: it names no coordinate,
  # so nothing on this route is reported aliased
  list(delta = as.numeric(backsolve(ch, forwardsolve(t(ch), u))),
       rank = ncol(A), route = route, held = held, dropped = integer(0))
}


#' The Pieces That One Scoring Step Needs
#'
#' @description
#' Assembles, at the current coefficients, whichever of the square-root
#' design, the penalty's factor and the penalized information that the
#' requested decomposition will actually use.
#'
#' @details
#' Only what the decomposition will read is built. `"qr"` and `"svd"` need
#' the square-root design `R` and the penalty's factor `C`, and build the
#' information only when one of the two factors does not exist; `"chol"` needs
#' the assembled `A` and neither square root. Building all three
#' unconditionally would pay for the route not taken at every iteration of the
#' loop.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, refreshed at `coef` if any term needs it.
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter.
#' @param hyper The hyperparameters, per penalized term.
#' @param method An [Iwls()] object, read for `solver` and `hessian`, whose
#'   `method` (see [distributions7::expected_info()]) gives the method of the
#'   expected information.
#'
#' @return A list with `R`, `C` and `A` and, on the square-root route,
#'   `threads`, the thread count. The matrices that the route does not need are
#'   `NULL`:
#'   \describe{
#'     \item{`R`}{the square-root design, as [sqrt_design()] returns it, or
#'       `NULL`.}
#'     \item{`C`}{the penalty's factor, as [penalty_sqrt()] returns it, or
#'       `NULL`.}
#'     \item{`A`}{the assembled penalized information, or `NULL`.}
#'   }
#'   `R` is also `NULL` when the per-observation curvature had no Cholesky
#'   factor, and [iwls_solve()] falls back to the assembled route there.
#'
#' @seealso [iwls_solve()], which consumes this,
#'   [sqrt_design()] and [penalty_sqrt()] for the two factors.
#'
#' @keywords internal
iwls_pieces <- function(spec, design, coef, hyper, method) {
  if (info_auto(method@hessian)) method <- iwls_resolve(method, spec@distrib)
  expected <- info_expected(method@hessian)
  emethod <- info_emethod(method@hessian)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian",
                          jeffreys = FALSE)
  # THE JEFFREYS CURVATURE IN TWO PARTS. The part with the sparsity of the
  # information is a per-observation weight on the same cross products, so it
  # is folded into the family's components and reaches both routes through
  # them; the rest is V V' with V of a few columns, which the square-root
  # route takes as further rows of the penalty's factor and the assembled
  # route by Woodbury's identity (iwls_solve()), so a sparse information is
  # never made dense by it.
  jp <- jfp_penalty(spec, design, coef, "parts")
  V <- NULL
  fold <- function(H) H
  if (!is.null(jp)) {
    if (!is.null(jp$refresh)) S <- S - jp$refresh
    V <- jp$lowrank
    fold <- function(H) {
      for (k in names(jp$weights)) {
        H[[k]] <- rep_len(H[[k]], spec@n_obs) - jp$weights[[k]]
      }
      H
    }
  }
  assembled <- function(H = NULL) {
    if (!is.null(jp)) {
      if (is.null(H)) {
        H <- statmod_family_hessian(spec, statmod_eta(spec, design,
                                                      coef)$theta,
                                    expected, emethod)
      }
      H <- fold(H)
    }
    list(R = NULL, C = NULL,
         A = statmod_information_at(spec, coef, design, expected,
                                    emethod, H = H) + S, V = V)
  }
  if (!(solver_how(method@solver) %in% c("qr", "svd"))) return(assembled())
  ev <- statmod_eta(spec, design, coef)
  # THE SQUARE-ROOT ROUTE IS ABANDONED AT EVERY POINT FOR SOME FAMILIES, and
  # the fallback below used to ask the family for the same components again.
  # Where the expected information is the outer product of scores the
  # per-observation block is -g g', of rank one, so it has no Cholesky factor
  # whenever the family carries more than one parameter -- measured on pig1
  # and pig2 before distributions7 0.65.0 made theirs exact, chol_blocks()
  # returned NULL for
  # every one of them at every point. Asking once and handing the answer to
  # both routes is what stops the abandoned one being paid for twice. A
  # mixture over regimes reads a different set per component, so it is left
  # to ask for its own.
  H <- if (length(ev$regimes)) NULL else
    statmod_family_hessian(spec, ev$theta, expected, emethod)
  L <- chol_blocks(info_blocks(spec, ev$theta, expected, emethod,
                               H = if (is.null(jp)) H else fold(H)))
  R <- sqrt_design(statmod_design_at(spec, coef, design), L)
  C <- penalty_sqrt(S)
  if (is.null(R) || is.null(C)) return(assembled(H))
  list(R = R, C = C, A = NULL, V = V, threads = spec@threads)
}


#' Fit the Smooth Block by Iterated Weighted Least Squares
#'
#' @description
#' Runs the scoring iteration on the objective of
#' [statmod_objective()] at fixed hyperparameters, returning the
#' coefficients and the history of the run.
#'
#' @details
#' # Repair of the curvature
#'
#' A curvature that is not positive definite is repaired by flooring its
#' eigenvalues, and the run continues. `solve()` would signal an error there
#' and abandon a start that is merely far from the optimum.
#'
#' # The stopping rule
#'
#' The built-in stopping rule compares a score divided by the sample size with
#' `tol`, so that a threshold has the same meaning at \eqn{n = 10} and at
#' \eqn{n = 10^7}, and it also requires the Newton decrement to be at most
#' \eqn{10^{-6}}. After the rule is met, one more full step is taken. The
#' objective itself is not averaged, because the penalty is added on that
#' scale.
#'
#' # The dimensionless final verdict
#'
#' The absolute score of a location equation carries the units \eqn{1/y}. On
#' a response three decades small its rounding floor sits above any fixed
#' threshold, and a run stalled exactly at the optimum reads as a failure.
#'
#' The verdict therefore includes a second, dimensionless reading,
#' \deqn{\max_j \lvert g_j \rvert / (n\,s_p), \qquad
#'       s_p = \sqrt{\mathrm{median}_j H_{jj} / n},}
#' where \eqn{s_p} is taken over the equation to which each coordinate belongs.
#'
#' It can relabel a run that has already stopped as converged, and it never
#' stops a run that is still moving. Driving the loop with it would make the
#' tolerance unreachable at the other end of the scale, where the magnitude of
#' the objective grows with \eqn{\log y} and the stall guard fires first.
#'
#' @param obj The objective, as [statmod_objective()] returns it.
#' @param start The starting coefficients, stacked into one numeric vector.
#' @param method An [Iwls()] object, carrying the curvature, the
#'   solver, the budget and the stopping rule.
#' @param n The number of observations, the divisor of the scaled stopping
#'   rule.
#' @param pieces_at A function of the stacked coefficients returning what
#'   [iwls_solve()] needs, as [iwls_pieces()] builds it.
#' @param verbose `TRUE` to print one line per iteration: the objective, the
#'   score per observation and the step length.
#' @param groups A list of integer index vectors, one per equation, in the
#'   numbering of the coefficients handed in. `NULL` treats them as one
#'   group, which is right where the whole vector is being fitted.
#' @param frozen Positions held at their starting values, an integer vector
#'   in the same numbering. They are dropped from every solve, so the run
#'   maximizes over the rest, and the stopping rule reads only the free
#'   coordinates, because the constrained optimum leaves in a held coordinate
#'   a score that does not vanish.
#' @param backup_at `NULL`, or a function, like `pieces_at`, that builds the
#'   pieces on the expected information. The expected information takes the
#'   step where the observed step finds no acceptable point. An observed
#'   penalized information that is not positive definite is first repaired by
#'   [pd_repair()] and its step tried. A trial point whose objective raises an
#'   error is read as rejected. [fit_smooth()] passes one for a method that
#'   [iwls_resolve()] settled with the fallback.
#' @param damp_on_reject Whether a line search that finds no acceptable step
#'   raises the Levenberg damping and retries (the default) or ends the run.
#'   [fit_smooth()] passes `FALSE` where the design carries a sharp
#'   break-point term, whose objective has a kink in the break-point at the
#'   observations: there every Gauss-Newton step from the minimum is
#'   rejected, and escalating the damping at every inner fit costs iterations
#'   without reaching anything that a stop would not.
#' @param kinks_at `NULL`, or a function of the coefficients returning the
#'   rows of the linear constraints that hold the objective on its kinks
#'   there, a matrix with one column per coefficient, as
#'   [kink_constraints()] finds them. Where a line search rejects every step
#'   length, the problem is solved on the subspace that those constraints and
#'   `frozen` leave free ([iwls_tied()]); from the point reached one step on
#'   every coordinate is tried, and where it is rejected too the run stops
#'   there, with the verdict of the constrained fit.
#' @param newton_at `NULL`, or a function of the coefficients returning the
#'   full Hessian of the objective, the scoring pieces' curvature plus the
#'   term that a block moving with its coefficients adds. Where the scoring step
#'   is rejected or shrunk below a tenth, the Newton step on it is tried by
#'   [iwls_newton_step()] and the one that decreases the objective more is
#'   taken.
#' @param switch_at `NULL`, or a function, like `pieces_at`, that builds the
#'   pieces on the observed information. After `switch_after` iterations without
#'   convergence, or earlier when [iwls_slow()] reports a slow contraction of
#'   the score, it takes the place of `pieces_at`, and `pieces_at` that of
#'   `backup_at`. [fit_smooth()] passes one for a method that [iwls_resolve()]
#'   settled on the expected information from [auto_info()].
#' @param switch_after The iteration after which `switch_at` replaces
#'   `pieces_at`.
#'
#' @return A list of ten: the six below; `note`, the reason for which a run
#'   stopped or `NULL`; `aliased`, the coordinates that the pivot left out;
#'   `switched`, the iteration at which the run moved to `switch_at`, or 0; and
#'   `fallback`, a named integer vector of four counts. The count `indefinite`
#'   is the number of iterations at which the observed penalized information was
#'   not positive definite, `search` the number at which the expected pieces
#'   stepped in for a step that found no acceptable point, `error` the number of
#'   trial points whose objective raised, and `newton` the number of iterations
#'   at which the Newton step on the full Hessian was taken in place of the
#'   scoring step. The six are: \describe{
#'     \item{`par`}{the stacked coefficients reached, a numeric vector.}
#'     \item{`value`}{the penalized objective there, not averaged over the
#'       observations.}
#'     \item{`converged`}{a single logical: whether a stopping rule fired,
#'       including the dimensionless verdict.}
#'     \item{`iterations`}{the index of the iteration at which the loop
#'       ended. `1` means the rule was met at the first check, as a warm start
#'       already at the mode reports; the full step taken after the rule is
#'       met is not counted.}
#'     \item{`score`}{the score per observation at the point reached.}
#'     \item{`history`}{a data frame with one row per accepted step, or `NULL`
#'       when no step was accepted.}
#'   }
#'
#' @seealso [iwls()] for the settings, [iwls_solve()] for one step,
#'   [iwls_score()] for the dimensionless reading.
#'
#' @keywords internal
iwls_fit <- function(obj, start, method, n, pieces_at, verbose = FALSE,
                     groups = NULL, frozen = integer(0), backup_at = NULL,
                     damp_on_reject = TRUE, kinks_at = NULL,
                     newton_at = NULL, switch_at = NULL,
                     switch_after = Inf) {
  beta <- start
  tol_def <- iwls_definite_tol(method)
  # the iteration at which the run moved to the observed information, 0 if
  # it never did
  switched <- 0L
  value <- obj$fn(beta)
  hist <- list()
  # a run that stops before its first solve has dropped nothing
  aliased <- integer(0)
  converged <- FALSE
  it <- 0L
  score <- Inf
  scores <- numeric(0)
  note <- NULL
  # NULL until a step has been taken: a rule reading a change in the
  # objective has nothing to read at the starting point, and says FALSE
  # there rather than comparing a number with itself
  f_old <- NULL
  x_old <- NULL
  # zero unless a step has failed, so the plain scoring iteration is what
  # every run that does not stall performs
  damp <- 0
  damp_tries <- 0L
  # how often the expected pieces stepped in, and for what
  fallback <- c(indefinite = 0L, search = 0L, error = 0L, newton = 0L)
  # the constrained fit on the subspace a kink leaves free, once one has run
  # from the current point; reset by every ordinary step that is accepted
  tied <- NULL
  # the equations' coordinate ranges, for the dimensionless reading of the
  # final verdict; the caller says them where the coefficients it hands in
  # are a subset, since the objective's own split maps the full vector
  if (is.null(groups)) groups <- list(seq_along(beta))

  if (verbose) {
    cat(sprintf("  %-5s %14s %12s %10s\n", "iter", "objective", "score/n",
                "step"))
  }
  for (it in seq_len(as.integer(method@maxit))) {
    # FISHER SCORING YIELDS TO NEWTON once it has taken `switch_after` steps
    # without converging. It converges only linearly near the mode wherever
    # the two informations differ there, which they do for every dispersion
    # with an equation of its own: measured on a gaussian with smooths in the
    # mean and in sigma, the first inner fit needed 235 scoring steps and the
    # REML criterion was unavailable at its start at the budget of 100. From
    # here the observed information takes the step and the expected one
    # stands in wherever the observed cannot, as on a family whose expected
    # information is not exact.
    g <- obj$gr(beta)
    # THE VERDICT READS THE FREE COORDINATES. A coordinate the caller holds
    # keeps whatever score the constrained optimum leaves in its direction --
    # that score is the quantity a Rao or a Terrell statistic is made of, so
    # it does not vanish and must not -- and reading it here would report a
    # perfectly converged constrained fit as still moving.
    gf <- g
    if (length(frozen)) gf[frozen] <- 0
    # a step whose OBJECTIVE is finite can still land where the score is not:
    # the line search below checks the value and nothing checked this, so a
    # non-finite gradient reached `if (score < tol)` and stopped the run with
    # "missing value where TRUE/FALSE needed", naming neither the iteration
    # nor the cause. The point is unusable and the last good one is kept.
    if (!all(is.finite(gf))) {
      note <- paste0("the score is not finite at iteration ", it,
                     ", so the run stopped at the last usable point")
      break
    }
    score <- max(abs(gf)) / n
    scores <- c(scores, score)
    # THE SWITCH, read off the run itself: Fisher scoring that contracts the
    # score slowly over two consecutive steps is in the linear regime where
    # Newton's step pays, whatever the count; the count is the ceiling.
    if (!switched && !is.null(switch_at) &&
        (it > switch_after || iwls_slow(scores, method@switch_rate))) {
      backup_at <- pieces_at
      pieces_at <- switch_at
      switched <- it
    }
    if (verbose) {
      cat(sprintf("  %-5d %14.6f %12.3e %10s\n", it - 1L, value, score,
                  if (it == 1L) "-" else fmt_step(step_used)))
    }
    met <- iwls_met(method, list(iter = it - 1L, f_new = value, f_old = f_old,
                                 x_new = beta, x_old = x_old, gradient = gf / n,
                                 stationarity = NULL))
    # THE MODE IS LOCATED IN LOG-LIKELIHOOD UNITS AS WELL. The score per
    # observation carries the units of the response: where the curvature
    # n / sigma^2 is small, a score under `tol` still leaves the mode far
    # away. Measured on a smooth with a random effect at a response of scale
    # 1000, the rule was met with the Newton decrement at 1.5e-02, and an
    # outer criterion read there differed between two searches by 0.017.
    # The built-in rule therefore also asks the decrement g'(H+S)^-1 g / 2,
    # which is scale-free, to be under the method's decrement_limit. A caller's
    # criterion is left as the caller wrote it.
    pc_met <- NULL
    if (met && is.null(method@criterion)) {
      pc_met <- pieces_at(beta)
      pcd <- if (!is.null(backup_at) && !pieces_definite(pc_met, tol_def))
        backup_at(beta) else pc_met
      dec <- iwls_decrement(pcd, g, method, damp, frozen)
      if (is.finite(dec) && dec > method@decrement_limit) met <- FALSE
    }
    if (met) {
      converged <- TRUE
      # ONE MORE STEP AT THE POINT THE RULE ACCEPTED. The marginal criterion
      # carries log|K(beta)|, which is not stationary in beta, so a mode left
      # short by delta enters it at FIRST order in delta. A warm-started fit
      # whose score is already under the rule stops without moving, and the
      # criterion is then read at the incumbent's coefficients: measured on a
      # gaussian location-scale smooth, the trial points of an outer line
      # search read the criterion off by an amount linear in the step,
      # -2.9e-05 halving with it, against a predicted gain of 1e-06, and the
      # search spent its twelve backtracks there. A full step on the observed
      # information takes the error to second order.
      if (iwls_polish()) {
        pp <- if (!is.null(switch_at) && !switched) switch_at(beta) else
          if (is.null(pc_met)) pieces_at(beta) else pc_met
        if (!pieces_definite(pp, tol_def)) pp <- NULL
        # The step is taken WHOLE and kept unless it makes the objective
        # clearly worse. At a point the rule accepted, the decrease a step
        # buys is at the objective's rounding, so a sufficient-decrease test
        # there accepts or rejects by the last bits, and the criterion read at
        # the mode then jumps between nearby hyperparameters: measured, a
        # central difference of aic() moved from 0.43769 to 0.43727 and
        # 0.43683 as the step halved, against 0.43769 at every step without it.
        if (!is.null(pp)) {
          sol <- tryCatch(iwls_solve(pp, -g, solver_how(method@solver), damp,
                                     frozen, tol = solver_tol(method@solver),
                                     floor = solver_floor(method@solver)),
                          error = function(e) NULL)
          if (!is.null(sol) && all(is.finite(sol$delta))) {
            cand <- beta + sol$delta
            vnew <- tryCatch(obj$fn(cand), error = function(e) Inf)
            if (is.finite(vnew) &&
                vnew <= value + 1e-8 * max(1, abs(value))) {
              beta <- cand
              value <- vnew
            }
          }
        }
      }
      break
    }
    pc <- if (is.null(pc_met)) pieces_at(beta) else pc_met
    # THE EXPECTED PIECES TAKE THE STEP where the observed ones cannot, on a
    # method iwls_resolve() settled with the fallback: an observed step that
    # finds no acceptable point. A trial point whose objective raises is then
    # read as rejected rather than as the end of the fit.
    #
    # AN OBSERVED PENALIZED INFORMATION THAT IS NOT POSITIVE DEFINITE IS
    # REPAIRED AND TRIED FIRST (iwls_solve() floors its eigenvalues with
    # pd_repair()); the expected pieces step in only if its search fails.
    # Until 0.197.0 they took every such step, and the expected information
    # was the slow route: on a Burr XII with every parameter over a factor
    # (rent, 1969 rows) the joint information is indefinite at 1207 steps of
    # 1323, and the fit took 221 s with the OPG standing in (a density-only
    # family) and 17.1 s with the closed forms, against 10.1 and 9.4 s now at
    # the same log-likelihood. On the seven families whose exact expected
    # information is costly, with a smooth on the location and the shape over
    # a covariate at n = 1000, the log-likelihood and the certificate's state
    # are the same to 10 digits at 1/4 to 1/6 of the time where the observed
    # information is indefinite (skew t 92.5 -> 15.4 s, pseudohuber2 254.4 ->
    # 57.2 s), and the repaired step's search did not fail once.
    if (!is.null(backup_at) && !pieces_definite(pc, tol_def)) {
      fallback[["indefinite"]] <- fallback[["indefinite"]] + 1L
    }
    st <- iwls_line_search(obj, beta, value, g, pc, method, damp, frozen,
                           guard = !is.null(backup_at))
    fallback[["error"]] <- fallback[["error"]] + st$errors
    if (!st$ok && !is.null(backup_at)) {
      pc <- backup_at(beta)
      fallback[["search"]] <- fallback[["search"]] + 1L
      st <- iwls_line_search(obj, beta, value, g, pc, method, damp, frozen,
                             guard = TRUE)
      fallback[["error"]] <- fallback[["error"]] + st$errors
    }
    # THE FULL CURVATURE WHERE A BLOCK MOVES WITH ITS COEFFICIENTS. The
    # scoring pieces carry the Gauss-Newton matrix, which leaves out the
    # score times the block's own derivative, and on a smoothed break-point
    # that term can be most of the curvature: measured on a jump() under the
    # quintic, the curvature in the break-point is 6.1 by Gauss-Newton and
    # 7853 in truth, so the scoring step in it is thousands of times too long,
    # the line search shrinks the whole step to 0.001, and the run crawls and
    # stalls 0.066 above the mode. Where the scoring step was shrunk that far
    # or rejected, the Newton step on the full Hessian is tried as well and
    # the better of the two is taken, so a run whose scoring steps are
    # accepted never sees it.
    if (!is.null(newton_at) && (!st$ok || st$step_used < 0.1)) {
      nt <- iwls_newton_step(obj, beta, value, g, newton_at(beta), frozen,
                             method@step_halving)
      if (!is.null(nt) && nt$ok && (!st$ok || nt$vnew < st$vnew)) {
        st[c("ok", "cand", "vnew", "step_used")] <-
          nt[c("ok", "cand", "vnew", "step_used")]
        fallback[["newton"]] <- fallback[["newton"]] + 1L
      }
    }
    sol <- st$sol
    # the coordinates the pivot left out, kept from the LAST step taken: a
    # design deficient at the fitted point was deficient at every point of
    # the run, the deficiency being a property of the columns and not of
    # where the coefficients are
    aliased <- sol$dropped
    ok <- st$ok
    cand <- st$cand
    vnew <- st$vnew
    step_used <- st$step_used
    # THE LINE SEARCH FOUND NO ACCEPTABLE STEP AT ANY LENGTH. That is the
    # deadlock in its purest form -- the direction cannot decrease the
    # objective however short it is made, so the quadratic it came out of is
    # wrong in some coordinate -- and it is what Levenberg's lambda is for:
    # the standard trust-region loop raises it on a rejected step and
    # retries. The loop used to break here, so the damping the method
    # already carries was never reached on this path, and a count model
    # whose dispersion runs to its Poisson limit exits exactly this way
    # after six halvings with the score still at 9.083e-03.
    if (!ok) {
      # ⚠️ NOT AT A KINK. A sharp break-point's objective is not
      # differentiable in the break-point at the observations, so at its
      # minimum every Gauss-Newton step is rejected whatever the damping:
      # measured on seg(x, psi ~ random(~1 | id)), the right and left slopes
      # there are +6.00 and -1.63 and the objective rises on both sides, and
      # escalating to a damping of 7.5e10 at every inner fit took the outer
      # search from 6 evaluations and 5.4 s of processor time to 28 and 25.1
      # s at an unchanged verdict. See piano_stabilita.txt 28d-28e.
      # A KINK, and the step that belongs there. Where the design carries a
      # sharp break-point sitting on an observation, the objective is not
      # differentiable in the coordinates that move it, the scoring direction
      # crosses the kink and no length of it decreases the objective -- while
      # every OTHER coordinate may still be far from its own optimum. Measured
      # on seg(x, psi ~ random(~1 | id)): the run stopped here with the other
      # groups' deviations at scores up to 1.6 and the objective 0.010 above
      # the mode, and holding the two coordinates that move the kinked
      # break-point reached the mode in four iterations. Damping does not
      # reach it: measured from the same point, it stops where it started.
      # The constraint is linear in the coefficients, so the problem is
      # solved on its null space: where the break-point's position is one
      # coefficient that coefficient is held, and where it is developed over
      # a design the coefficients move along the hyperplane that keeps the
      # break-point on the observation. Once the constrained fit has run
      # from here, a further rejection means the point is its optimum.
      if (!is.null(kinks_at)) {
        if (!is.null(tied)) {
          converged <- isTRUE(tied$converged)
          score <- tied$score
          note <- paste0("the objective has a kink at the point reached; ",
                         "the fit is the constrained one on the subspace ",
                         "that keeps it there")
          break
        }
        Ak <- kinks_at(beta)
        if (NROW(Ak)) {
          tied <- iwls_tied(obj, beta, method, n, pieces_at, Ak, frozen,
                            backup_at = backup_at, newton_at = newton_at,
                            switch_at = switch_at, switch_after = switch_after,
                            kinks_at = kinks_at,
                            damp_on_reject = damp_on_reject)
          if (!is.null(tied)) {
            if (tied$value <= value) {
              x_old <- beta
              f_old <- value
              beta <- tied$par
              value <- tied$value
            }
            next
          }
        }
      }
      if (damp_on_reject && damp_tries < 8L) {
        damp <- iwls_escalate(damp, pc)
        damp_tries <- damp_tries + 1L
        next
      }
      note <- paste0("the line search found no acceptable step at iteration ",
                     it, if (damp_on_reject) ", even with damping")
      break
    }
    hist[[length(hist) + 1L]] <- data.frame(
      iteration = it, objective = vnew, score = score, step = step_used,
      rank = sol$rank, route = sol$route, held = length(sol$held),
      damp = damp
    )
    # A step accepted for a decrease below the rounding of the objective has
    # not moved it, and the sufficient-decrease test cannot tell the two
    # apart: with a step of 1e-9 the decrease it asks for is of that order,
    # which any point satisfies. The run then spends its budget standing
    # still, which is what a term whose gradient belongs to a working model
    # rather than to the objective produces at its own fixed point.
    stalled <- value - vnew <= 1e-12 * max(1, abs(value))
    tied <- NULL
    x_old <- beta
    f_old <- value
    beta <- cand
    value <- vnew
    # A STALL WITH THE SCORE STILL LARGE IS ONE COORDINATE HOLDING THE OTHERS,
    # and a shorter step cannot cure it because the step length is scalar.
    # Where a coordinate's curvature is orders below its neighbors' -- a
    # shape approaching the clamp its link keeps it inside, whose information
    # falls as nu^-4 -- the scoring step in it is astronomically long, and the
    # line search shortens the WHOLE step to keep it admissible: measured on a
    # Student t at nu = 9.1e+12, the step was -1790 in nu against at most 3.9
    # in every mean coordinate, the search ran 1, 0.125, 1.5e-05, 1.5e-08 and
    # the run stalled at a score of 2.9e-02 with the mean nowhere near
    # stationary.
    #
    # Levenberg's lambda shortens that coordinate in proportion to how little
    # curvature it has and leaves the others where they were, which is the
    # differential shrink the step length cannot give. It is raised only
    # HERE, where the loop used to give up, so a run that never stalls never
    # sees it and is unchanged; it decays again on every step that moves.
    #
    # ⚠️ ONLY WHERE THE STEP WAS SHRUNK, which is the signature of the
    # deadlock and not of every stall. A term whose block is a working
    # linearization stalls AT ITS OWN FIXED POINT with the full step taken --
    # the gradient there belongs to the working model rather than to the
    # objective -- and damping past it costs iterations and ends a run that
    # had arrived: measured, escalating on every stall turned a converged
    # break-point fit into a non-converged one and moved three exact-gradient
    # checks from machine precision to 1e-3, the mode being left worse
    # located.
    # The deadlock is the case where the line search had to shrink the whole
    # step by orders to keep one coordinate admissible.
    #
    # ⚠️ Broadening this to `step_used < 1` was measured on seven models --
    # the count model above, a smooth with a random effect, two equations, a
    # Poisson smooth, a Gamma smooth, a Student t at its boundary and a
    # break-point term -- and changed NOTHING: the same log-likelihood to
    # six decimals, the same flags, the same variance matrices. The exit the
    # runaway takes is the one below, not this one.
    if (stalled) {
      if (step_used < 1e-3 && damp_tries < 8L) {
        damp <- iwls_escalate(damp, pc)
        damp_tries <- damp_tries + 1L
        next
      }
      note <- paste0("the objective stopped moving at iteration ", it)
      break
    }
    if (damp > 0) {
      damp <- damp / 100
      if (damp < iwls_scale(pc) * 1e-12) damp <- 0
    }
  }
  # The rule is asked once more at the point reached, so that a run which
  # could not move is still converged when the POINT says so. Two readings,
  # in order. The absolute score first, which is the rule that drove the
  # loop. Where it fails, the DIMENSIONLESS score arbitrates: the absolute
  # rule carries the units of the response (a location equation's score is
  # 1/y), so on a response three decades small its rounding floor sits
  # above the threshold and a run that stalled AT the optimum reads as a
  # failure -- measured, 1.37e-6 against 1e-6 at y ~ 1e-3, the objective no
  # longer moving and the fit right. The dimensionless reading, each
  # equation's score against its own information scale, survives the units.
  # It only ever RELABELS a run that has already stopped; it never stops
  # one, so every trajectory is the absolute rule's -- the lesson of
  # fit_distrib()'s dropped crit_rel_obj, and of a first version of this
  # rule that drove the loop dimensionless and made the tolerance
  # unreachable at y * 1e4, where the stall guard on the objective (whose
  # magnitude grows with log y) fires before a rule 1/s_p stricter can.
  # A caller's rule gets the same second reading and NOT the dimensionless
  # relabel, which is the built-in rule's own repair of its own units.
  if (!converged) {
    gfin <- obj$gr(beta)
    # the verdict reads the free coordinates, as the loop's own test does
    if (length(frozen)) gfin[frozen] <- 0
    if (all(is.finite(gfin))) {
      score <- max(abs(gfin)) / n
      converged <- iwls_met(method, list(iter = it, f_new = value,
                                         f_old = f_old, x_new = beta,
                                         x_old = x_old, gradient = gfin / n,
                                         stationarity = NULL))
      if (!converged && is.null(method@criterion)) {
        hj <- iwls_info_diag(pieces_at(beta))
        dimless <- iwls_score(gfin, hj, groups, n)
        if (dimless < method@tol) {
          converged <- TRUE
          note <- paste0("converged by the dimensionless rule (score ",
                         format(dimless, digits = 3), " against ",
                         format(score, digits = 3),
                         " absolute; the response's scale is the gap)")
        }
      }
    } else {
      score <- Inf
    }
  }
  if (verbose) {
    cat(sprintf("  %-5s %14.6f %12.3e %10s\n", "end", value, score,
                if (converged) "converged" else "stopped"))
    if (!is.null(note)) cat("  ", note, "\n", sep = "")
  }
  list(par = beta, value = value, converged = converged,
       iterations = it, score = score, note = note, aliased = aliased,
       fallback = fallback, switched = switched,
       history = if (length(hist)) do.call(rbind, hist) else NULL)
}


#' The Scoring Fit on the Subspace That a Kink Leaves Free
#'
#' @description
#' Runs [iwls_fit()] on \eqn{\beta = \beta_0 + N w}, with \eqn{N} a basis of
#' the directions that the kink constraints and the held coordinates leave
#' free ([null_basis()]), and returns the point reached in the coefficients.
#'
#' @details
#' The reparametrization is linear, so the reduced problem's pieces are the
#' full problem's pieces carried through \eqn{N}: the square-root design and
#' the penalty's factor are multiplied by \eqn{N} on the right, the assembled
#' matrix becomes \eqn{N^\top A N}, the gradient \eqn{N^\top g} and the full
#' Hessian \eqn{N^\top H N}. A kink met in the reduced problem adds its row,
#' carried through \eqn{N}, in the same way.
#'
#' @param obj The objective, with `fn` and `gr`.
#' @param beta The point \eqn{\beta_0}.
#' @param method The [Iwls()] method.
#' @param n The number of observations.
#' @param pieces_at,backup_at,newton_at,switch_at,kinks_at As [iwls_fit()]
#'   takes them, in the full coefficients.
#' @param Ak The rows of the kink constraints at `beta`.
#' @param frozen The positions held.
#' @param switch_after,damp_on_reject As [iwls_fit()] takes them.
#'
#' @return `NULL` where the constraints leave no direction free, and
#'   otherwise a list with `par`, the coefficients reached, `value`,
#'   `converged` and `score`, the reduced fit's.
#'
#' @seealso [iwls_fit()], [kink_constraints()]
#'
#' @keywords internal
iwls_tied <- function(obj, beta, method, n, pieces_at, Ak, frozen,
                      backup_at = NULL, newton_at = NULL, switch_at = NULL,
                      switch_after = Inf, kinks_at = NULL,
                      damp_on_reject = TRUE) {
  p <- length(beta)
  E <- matrix(0, length(frozen), p)
  if (length(frozen)) E[cbind(seq_along(frozen), frozen)] <- 1
  N <- null_basis(rbind(as.matrix(Ak), E), p)
  if (!ncol(N)) return(NULL)
  b0 <- beta
  at <- function(w) b0 + as.numeric(N %*% w)
  carry <- function(pc) {
    if (is.null(pc)) return(NULL)
    if (!is.null(pc$R)) pc$R <- pc$R %*% N
    if (!is.null(pc$C)) pc$C <- pc$C %*% N
    if (!is.null(pc$A)) pc$A <- crossprod(N, as_dense(pc$A) %*% N)
    if (!is.null(pc$V)) pc$V <- crossprod(N, pc$V)
    pc
  }
  wrap <- function(f) if (is.null(f)) NULL else function(w) carry(f(at(w)))
  sub <- list(fn = function(w) obj$fn(at(w)),
              gr = function(w) as.numeric(crossprod(N, obj$gr(at(w)))))
  nt <- if (is.null(newton_at)) NULL else function(w) {
    H <- newton_at(at(w))
    if (is.null(H)) NULL else crossprod(N, as_dense(H) %*% N)
  }
  kt <- if (is.null(kinks_at)) NULL else function(w) {
    A <- kinks_at(at(w))
    if (!NROW(A)) return(A)
    A <- as.matrix(A) %*% N
    A[rowSums(abs(A)) > 1e-12 * max(1, max(abs(A))), , drop = FALSE]
  }
  r <- iwls_fit(sub, numeric(ncol(N)), method, n, wrap(pieces_at),
                verbose = FALSE, groups = list(seq_len(ncol(N))),
                frozen = integer(0), backup_at = wrap(backup_at),
                damp_on_reject = damp_on_reject, kinks_at = kt,
                newton_at = nt, switch_at = wrap(switch_at),
                switch_after = switch_after)
  list(par = at(r$par), value = r$value, converged = r$converged,
       score = r$score)
}


#' One Scoring Step and Its Line Search
#'
#' @description
#' Solves for the scoring increment from the given pieces and halves the step
#' until Armijo's condition holds or the budget of halvings is spent. The full
#' step is tried first, so a budget of `step_halving` halvings allows
#' `step_halving + 1` trial points.
#'
#' @details
#' It is the body that [iwls_fit()] ran inline, moved into a function so that
#' the step can be tried twice at one iterate, first on the observed information
#' and then on the expected one, where [iwls_resolve()] settled the method with
#' the fallback. With `guard = FALSE` the operations are the inline ones, in the
#' same order. With `guard = TRUE` a trial point whose objective raises an error
#' is read as rejected rather than stopping the fit, and counted.
#'
#' @param obj The objective.
#' @param beta The current coefficients.
#' @param value The objective there.
#' @param g The gradient there.
#' @param pc The pieces that [iwls_solve()] reads.
#' @param method The [Iwls()] object.
#' @param damp The Levenberg damping.
#' @param frozen The held positions.
#' @param guard Whether an error at a trial point is read as a rejection.
#'
#' @return A list of `ok`, `cand`, `vnew`, `step_used`, `sol` and `errors`,
#'   the number of trial points whose objective raised.
#'
#' @seealso [iwls_fit()], its caller.
#'
#' @keywords internal
iwls_line_search <- function(obj, beta, value, g, pc, method, damp, frozen,
                             guard = FALSE) {
  sol <- iwls_solve(pc, -g, solver_how(method@solver), damp, frozen,
                    tol = solver_tol(method@solver),
                    floor = solver_floor(method@solver))
  delta <- sol$delta
  delta[!is.finite(delta)] <- 0
  # sufficient decrease, not mere non-increase: see armijo_ok()
  gd <- sum(g * delta)
  step_used <- 1
  ok <- FALSE
  cand <- beta
  vnew <- NA_real_
  errors <- 0L
  for (h in seq_len(as.integer(method@step_halving) + 1L)) {
    cand <- beta + step_used * delta
    vnew <- if (guard) {
      tryCatch(obj$fn(cand), error = function(e) {
        errors <<- errors + 1L
        Inf
      })
    } else {
      obj$fn(cand)
    }
    if (armijo_ok(vnew, value, step_used, gd)) {
      ok <- TRUE
      break
    }
    step_used <- step_used / 2
  }
  list(ok = ok, cand = cand, vnew = vnew, step_used = step_used, sol = sol,
       errors = errors)
}


#' Whether Scoring Pieces Carry a Positive Definite Curvature
#'
#' @description
#' Returns `TRUE` where the pieces are square roots, which are positive
#' semidefinite by construction, or where the assembled matrix has its
#' smallest eigenvalue above a relative floor. [iwls_fit()] uses the result to
#' check whether the observed pieces can take the step.
#'
#' @param pc The pieces, as [iwls_pieces()] returns them.
#' @param rel The relative floor, [pd_repair()]'s own by default; [iwls_fit()]
#'   passes [iwls_definite_tol()].
#'
#' @return A single logical. The result is `TRUE` for a matrix with a
#'   non-finite entry; holding the boundary is left to [iwls_solve()].
#'
#' @seealso [iwls_fit()]
#'
#' @keywords internal
pieces_definite <- function(pc, rel = 1e-8) {
  if (is.null(pc$A)) return(TRUE)
  A <- lowrank_dense(pc$A, pc$V)
  if (!all(is.finite(A))) return(TRUE)
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  min(ev) > rel * max(abs(ev))
}


#' Armijo's Sufficient-Decrease Condition
#'
#' @description
#' Tests whether a trial point may be accepted. The trial point is accepted
#' when the objective has fallen by at least a fixed share of what the model
#' predicted over that step.
#'
#' @details
#' \deqn{f(x + s d) \;\le\; f(x) + c_1 s\, g^\top d .}
#' The increment solves \eqn{(H + S)\delta = -g}, so
#' \eqn{g^\top\delta = -g^\top(H+S)^{-1}g} is negative and the bound lies
#' below \eqn{f(x)}, which requires a decrease. If the term were subtracted
#' instead, the bound would lie above \eqn{f(x)} and the test would permit an
#' increase of \eqn{c_1 s \lvert g^\top\delta\rvert}. The condition is a
#' separate function so that its sign can be tested directly.
#'
#' @param vnew The objective at the trial point.
#' @param value The objective at the current one.
#' @param step The step length.
#' @param gd \eqn{g^\top\delta}, negative for a descent direction.
#' @param c1 The share of the predicted decrease that the trial point must
#'   achieve.
#'
#' @return A single logical.
#'
#' @seealso [iwls_fit()], the loop that applies it.
#'
#' @keywords internal
armijo_ok <- function(vnew, value, step, gd, c1 = 1e-4) {
  is.finite(vnew) && vnew <= value + c1 * step * gd
}


#' One Newton Step on the Full Hessian and Its Line Search
#'
#' @description
#' Solves for the Newton increment on a given Hessian over the free
#' coordinates and halves the step until Armijo's condition holds or the
#' budget of halvings is spent.
#'
#' @details
#' The Hessian of an objective whose design moves with its coefficients need
#' not be positive definite away from the mode, so its eigenvalues are
#' replaced by their absolute values and floored at \eqn{10^{-8}} times the
#' largest, which keeps the increment a descent direction. The repair
#' differs from the one of [optimizers7::newton()], which raises each
#' eigenvalue below the floor to the floor without taking absolute values.
#'
#' @param obj The objective.
#' @param beta The current coefficients.
#' @param value The objective there.
#' @param g The gradient there.
#' @param H The full Hessian there.
#' @param frozen The held positions.
#' @param halvings The budget of step halvings.
#'
#' @return `NULL` where the Hessian is not usable, and otherwise a list of
#'   `ok`, `cand`, `vnew` and `step_used`.
#'
#' @seealso [iwls_fit()], its caller.
#'
#' @keywords internal
iwls_newton_step <- function(obj, beta, value, g, H, frozen, halvings) {
  H <- tryCatch(as_dense(H), error = function(e) NULL)
  if (is.null(H) || !all(is.finite(H))) return(NULL)
  free <- setdiff(seq_along(beta), frozen)
  if (!length(free)) return(NULL)
  e <- eigen((H[free, free, drop = FALSE] + t(H[free, free, drop = FALSE])) / 2,
             symmetric = TRUE)
  lam <- abs(e$values)
  if (!length(lam) || max(lam) <= 0) return(NULL)
  lam <- pmax(lam, 1e-8 * max(lam))
  delta <- numeric(length(beta))
  delta[free] <- -as.numeric(e$vectors %*% (crossprod(e$vectors, g[free]) / lam))
  gd <- sum(g * delta)
  step_used <- 1
  for (h in seq_len(as.integer(halvings) + 1L)) {
    cand <- beta + step_used * delta
    vnew <- tryCatch(obj$fn(cand), error = function(e) Inf)
    if (armijo_ok(vnew, value, step_used, gd)) {
      return(list(ok = TRUE, cand = cand, vnew = vnew, step_used = step_used))
    }
    step_used <- step_used / 2
  }
  list(ok = FALSE, cand = beta, vnew = value, step_used = step_used)
}


#' The Curvature's Own Scale at One Point
#'
#' @description
#' Returns the largest diagonal entry of the penalized information, read off
#' the pieces without assembling it.
#'
#' @details
#' The Levenberg damping is measured against this, so [iwls_escalate()]
#' carries no constant with units of its own.
#'
#' On the augmented route the diagonal of \eqn{K = R'R + C'C} is the column
#' sums of the squares of \eqn{R} and \eqn{C}, which reads a sparse design
#' without densifying it. On the assembled route it is the matrix's own
#' diagonal.
#'
#' @param pieces What [iwls_pieces()] built, carrying `R` and `C`, or `A`.
#'
#' @return A single positive number: the largest diagonal entry. The value is
#'   `1` where the pieces carry nothing usable, which makes the damping
#'   schedule absolute there instead of relative.
#'
#' @seealso [iwls_escalate()], which divides by this,
#'   [iwls_pieces()] for the input.
#'
#' @keywords internal
iwls_scale <- function(pieces) {
  v <- if (!is.null(pieces$R)) {
    as.numeric(Matrix::colSums(pieces$R^2)) +
      (if (is.null(pieces$C)) 0 else as.numeric(Matrix::colSums(pieces$C^2))) +
      (if (is.null(pieces$V)) 0 else rowSums(pieces$V^2))
  } else if (!is.null(pieces$A)) {
    as.numeric(Matrix::diag(pieces$A)) +
      if (is.null(pieces$V)) 0 else rowSums(pieces$V^2)
  } else {
    return(1)
  }
  v <- v[is.finite(v)]
  if (!length(v)) return(1)
  m <- max(v)
  if (is.finite(m) && m > 0) m else 1
}


#' Raise the Levenberg Damping
#'
#' @description
#' Returns the next \eqn{\lambda} to try after a step that failed or moved
#' nothing.
#'
#' @details
#' The first value is \eqn{10^{-8}} of the curvature's own scale, as
#' [iwls_scale()] measures it. That is negligible beside a well-curved
#' coordinate and already large beside one whose information has fallen by
#' eight orders, which is the disparity that the damping exists to correct.
#'
#' Each further escalation multiplies by a hundred, so that the eight
#' escalations that a run may use span fourteen orders of magnitude, from
#' \eqn{10^{-8}} to \eqn{10^{6}} times the largest diagonal.
#'
#' @param damp The current damping, a single non-negative number. `0` starts
#'   the schedule.
#' @param pieces What [iwls_pieces()] built, read only for the scale.
#'
#' @return The next damping to try, a single positive number.
#'
#' @seealso [iwls_scale()] for the scale, [iwls_solve()] for what the damping
#'   does to the system.
#'
#' @keywords internal
iwls_escalate <- function(damp, pieces) {
  if (damp > 0) return(damp * 100)
  iwls_scale(pieces) * 1e-8
}


#' The Newton Decrement of a Scoring Step
#'
#' @description
#' `iwls_decrement()` is \eqn{\tfrac{1}{2} g^\top (H+S)^{-1} g} over the free
#' coordinates, read from the step that [iwls_solve()] takes: the decrease in
#' the penalized log-likelihood that a full scoring step predicts.
#'
#' @details
#' It is in log-likelihood units, so it does not depend on the scale of the
#' response, where the score per observation does. [iwls_fit()]'s built-in
#' rule reads the mode as located when it is at most the `decrement_limit` of
#' [iwls()], whose default \eqn{10^{-6}} is three orders under
#' [mode_error_limit()], the reading used by the outer search to determine
#' whether a criterion is usable.
#'
#' @param pc Scoring pieces at the current coefficients.
#' @param g The gradient of the objective there.
#' @param method An [Iwls()] object.
#' @param damp The current Levenberg damping.
#' @param frozen Positions held out of the step.
#'
#' @return A number, `NA` where the solve fails.
#'
#' @keywords internal
iwls_decrement <- function(pc, g, method, damp, frozen) {
  sol <- tryCatch(iwls_solve(pc, -g, solver_how(method@solver), damp, frozen,
                             tol = solver_tol(method@solver),
                             floor = solver_floor(method@solver)),
                  error = function(e) NULL)
  if (is.null(sol)) return(NA_real_)
  delta <- sol$delta
  delta[!is.finite(delta)] <- 0
  gf <- g
  if (length(frozen)) gf[frozen] <- 0
  -0.5 * sum(gf * delta)
}




#' Has the Step's Stopping Rule Been Met?
#'
#' @description
#' Applies the stopping rule with which the step was configured: the built-in
#' one, the score per observation against `tol`, or the caller's
#' \pkg{optimizers7} criterion read on the same state.
#'
#' @details
#' Both routes are implemented here, although the loop evaluates the rule at
#' two places. This keeps the state that a criterion receives identical at
#' both places.
#'
#' In the state that a rule receives, `state$gradient` is already the score
#' per observation, so `crit_grad(t)` and `tol = t` apply the same threshold
#' to the score; the built-in rule also requires a small Newton decrement,
#' which a criterion does not receive. `state$f_new` is the negative
#' penalized log-likelihood, not divided by the sample size, which is the
#' scale on which the penalty is added.
#'
#' @param method An [Iwls()] object.
#' @param state The iteration state, as `optimizers7::crit_met`
#'   documents it.
#'
#' @return A single logical.
#'
#' @keywords internal
iwls_met <- function(method, state) {
  if (is.null(method@criterion)) {
    return(max(abs(state$gradient)) < method@tol)
  }
  isTRUE(optimizers7::crit_met(method@criterion, state))
}

#' The Dimensionless Reading of the Stopping Rule
#'
#' @description
#' Returns the absolute \eqn{\max\lvert g\rvert / n} read per equation
#' against that equation's own information scale,
#' \eqn{\max_j \lvert g_j\rvert / (n\,s_p)} with
#' \eqn{s_p = \sqrt{\mathrm{median}_j H_{jj} / n}} over the equation's
#' coordinates.
#'
#' @details
#' # Why the division is per equation
#'
#' The score of a location equation carries the units \eqn{1/y} and its
#' curvature \eqn{1/y^2}, so dividing by \eqn{s_p} survives any rescaling of
#' the response. Taking \eqn{s_p} over the whole equation changes nothing
#' within it: a stiff, heavily penalized coordinate is held to the same rule
#' as its neighbors, and the envelope identities on which the outer gradient
#' rests require this of the mode.
#'
#' This reading arbitrates the final verdict and never drives the loop.
#'
#' # Alternative readings that were rejected
#'
#' A per-coordinate normalization by \eqn{(H+S)_{jj}} is dimensionless too,
#' but at extreme shrinkage the penalized coordinates converge loosely and the
#' envelope identities cease to hold.
#'
#' When the per-equation form drives the loop, the rule works at small
#' \eqn{y} and cannot be met at large \eqn{y}. At \eqn{y \times 10^4} the
#' magnitude of the objective grows with \eqn{\log y}, so the stall guard on the
#' objective stops the run before a rule \eqn{1/s_p} stricter can be met.
#'
#' @param g The gradient at the point, a numeric vector over the stacked
#'   coefficients.
#' @param hj The information's diagonal at the same point, as
#'   [iwls_info_diag()] returns it, the same length as `g`.
#' @param groups A list of integer index vectors, one per equation, giving
#'   the coordinates belonging to that equation.
#' @param n The number of observations.
#'
#' @return A single non-negative number: the largest per-equation reading.
#'   An equation whose median curvature is zero or not finite is scaled by
#'   one, so that its score per observation enters the reading as it is.
#'
#' @seealso [fit_smooth()], whose verdict this arbitrates,
#'   [iwls()] for the absolute rule that drives the loop.
#'
#' @keywords internal
iwls_score <- function(g, hj, groups, n) {
  out <- 0
  for (idx in groups) {
    if (!length(idx)) next
    # a negative median curvature (a Student t's observed information at
    # its outliers) is tested before the root, which would warn
    m <- stats::median(hj[idx])
    s <- if (is.finite(m) && m > 0) sqrt(m / n) else 1
    out <- max(out, max(abs(g[idx])) / (n * s))
  }
  out
}

#' The Diagonal of the Information That a Step Uses
#'
#' @description
#' Returns the diagonal of the information by which [iwls_score()] normalizes,
#' one entry per stacked coefficient.
#'
#' @details
#' On the augmented route it is the column sums of the squares of the
#' square-root design, whose cross-product is the unpenalized information, so
#' the penalty does not enter.
#'
#' On the assembled route it is the diagonal of the matrix itself, which
#' folds the penalty in. That is a different normalizer. It is accepted
#' because the reading only arbitrates a verdict and does not drive the loop,
#' and because the assembled route is used either where the square
#' root does not exist or when `"chol"` is requested.
#'
#' @param pieces The pieces, as [iwls_pieces()] builds them, carrying `R` or
#'   `A`.
#'
#' @return An unnamed numeric vector, one entry per stacked coefficient.
#'   [iwls_pieces()] always builds one of the two matrices, so there is no
#'   third branch.
#'
#' @seealso [iwls_score()], which normalizes by it, and [iwls_fit()], the only
#'   caller.
#'
#' @keywords internal
iwls_info_diag <- function(pieces) {
  if (!is.null(pieces$R)) {
    return(as.numeric(Matrix::colSums(pieces$R^2)) +
             if (is.null(pieces$V)) 0 else rowSums(pieces$V^2))
  }
  as.numeric(diag(as_dense(pieces$A))) +
    if (is.null(pieces$V)) 0 else rowSums(pieces$V^2)
}


#' Floor the Eigenvalues of a Curvature Matrix
#'
#' @description
#' Returns the matrix unchanged when it is positive definite, and otherwise
#' returns the matrix rebuilt from its eigendecomposition after the
#' eigenvalues are floored.
#'
#' @details
#' A start at an indefinite curvature has to be abandoned if `solve()` is
#' used, and the repair avoids this. The floor is relative to the largest
#' absolute eigenvalue, since an absolute floor has no meaning across scales.
#' Positive definiteness is tested by a Cholesky factorization.
#'
#' @param A A symmetric matrix.
#' @param rel The floor, relative to the largest absolute eigenvalue.
#'
#' @return A symmetric positive definite matrix, or `NULL` where `A` has a
#'   non-finite entry. A matrix with no rows is returned as it is.
#'
#' @keywords internal
pd_repair <- function(A, rel = 1e-8) {
  if (nrow(A) == 0L) return(A)
  # a curvature with a non-finite entry is not a matrix to repair: eigen()
  # signals "infinite or missing values in 'x'" from three frames down, which
  # names the arithmetic and not the iterate that produced it
  if (!all(is.finite(A))) return(NULL)
  # the factorization is a probe: on a sparse matrix that is not positive
  # definite CHOLMOD warns before Matrix signals the error, and since 0.197.0
  # an indefinite observed information reaches this test inside a fit, where
  # the warning reached the user three times per fit
  ok <- tryCatch({
    suppressWarnings(chol(A))
    TRUE
  }, error = function(e) FALSE)
  if (ok) return(A)
  e <- eigen(A, symmetric = TRUE)
  floor_ <- rel * max(abs(e$values))
  if (!is.finite(floor_) || floor_ <= 0) floor_ <- rel
  e$vectors %*% (pmax(e$values, floor_) * t(e$vectors))
}


#' Format a Step Length
#'
#' @param x A number.
#'
#' @return A single string.
#'
#' @keywords internal
fmt_step <- function(x) {
  if (is.null(x) || !is.finite(x)) return("-")
  formatC(x, format = "g", digits = 3)
}


#' A Solve With a Matrix Plus a Product of Low Rank
#'
#' @description
#' `lowrank_solve()` solves \eqn{(A + VV^\top + \lambda I)\delta = u} by
#' Woodbury's identity, factoring \eqn{A + \lambda I} in its own storage and
#' a small matrix of the size of \eqn{V}'s columns. `lowrank_dense()` returns
#' \eqn{A + VV^\top} as a base matrix.
#'
#' @details
#' With \eqn{A} positive definite,
#' \deqn{(A + VV^\top)^{-1}u = A^{-1}u - A^{-1}V(I + V^\top A^{-1}V)^{-1}
#'   V^\top A^{-1}u,}
#' and the matrix in the middle is positive definite. The Jeffreys curvature
#' enters the assembled route in this form ([iwls_pieces()]), so a sparse
#' information keeps its sparse factorization. Where \eqn{A} has no Cholesky
#' factor the function returns `NULL`, and [iwls_solve()] assembles the sum
#' and repairs it as it repairs any indefinite matrix.
#'
#' @param A A symmetric matrix, dense or a \pkg{Matrix} object.
#' @param V A matrix with as many rows as `A`, or `NULL`.
#' @param u The right-hand side.
#' @param damp The Levenberg damping added to the diagonal.
#'
#' @return `lowrank_solve()`: a numeric vector, or `NULL`.
#'   `lowrank_dense()`: a base matrix.
#'
#' @seealso [iwls_solve()]
#'
#' @keywords internal
lowrank_solve <- function(A, V, u, damp = 0) {
  if (is.null(V) || !ncol(V)) return(NULL)
  if (!all(is.finite(V))) return(NULL)
  if (damp > 0) Matrix::diag(A) <- Matrix::diag(A) + damp
  if (isS4(A)) {
    ch <- tryCatch(Matrix::Cholesky(Matrix::forceSymmetric(
      methods::as(A, "CsparseMatrix")), LDL = FALSE, perm = TRUE),
      error = function(e) NULL)
    if (is.null(ch)) return(NULL)
    sol <- function(b) as.matrix(Matrix::solve(ch, b, system = "A"))
  } else {
    ch <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(ch)) return(NULL)
    sol <- function(b) backsolve(ch, forwardsolve(t(ch), b))
  }
  Au <- sol(u)
  AV <- sol(V)
  cap <- diag(ncol(V)) + crossprod(V, AV)
  cc <- tryCatch(chol((cap + t(cap)) / 2), error = function(e) NULL)
  if (is.null(cc)) return(NULL)
  w <- backsolve(cc, forwardsolve(t(cc), crossprod(V, Au)))
  out <- as.numeric(Au - AV %*% w)
  if (!all(is.finite(out))) return(NULL)
  out
}

#' @rdname lowrank_solve
#' @keywords internal
lowrank_dense <- function(A, V) {
  A <- as_dense(A)
  if (is.null(V) || !ncol(V)) return(A)
  A + tcrossprod(V)
}


#' The Augmented Solve With a Product of Low Rank
#'
#' @description
#' Solves \eqn{(R^\top R + C^\top C + VV^\top)\delta = u} from one
#' decomposition of \eqn{[R; C]}, by Woodbury's identity.
#'
#' @details
#' [augmented_solve()] solves the right-hand side and the columns of \eqn{V}
#' with one factorization, and the correction is a system of the size of
#' \eqn{V}'s columns. Appending \eqn{V^\top} as rows of the augmented
#' matrix would give the same solution, but each such row is dense, and a
#' sparse QR decomposition fills in completely beside it: measured on a
#' random intercept over 250 levels, the Jeffreys fit took 13.9 s against
#' 3.1 s without the block, almost all of it in that decomposition. A
#' coordinate that the pivot leaves out is read as zero in every column.
#'
#' @param R,C The square-root design and the penalty's factor.
#' @param V The low-rank factor, or `NULL` for the plain augmented solve.
#' @param u The right-hand side.
#' @param how,threads,tol As in [augmented_solve()].
#'
#' @return What [augmented_solve()] returns.
#'
#' @seealso [lowrank_solve()], the same identity on the assembled route.
#'
#' @keywords internal
augmented_lowrank <- function(R, C, V, u, how, threads = 1L, tol = NULL) {
  if (is.null(V) || !ncol(V)) {
    return(augmented_solve(R, C, u, how, threads = threads, tol = tol))
  }
  out <- augmented_solve(R, C, cbind(u, V), how, threads = threads, tol = tol)
  D <- as.matrix(out$delta)
  Au <- D[, 1L]
  AV <- D[, -1L, drop = FALSE]
  cap <- diag(ncol(V)) + crossprod(V, AV)
  w <- tryCatch(solve((cap + t(cap)) / 2, crossprod(V, Au)),
                error = function(e) NULL)
  if (is.null(w)) {
    out$delta <- Au
    return(out)
  }
  out$delta <- as.numeric(Au - AV %*% w)
  out
}
