# The restricted fit, and what the four likelihood statistics read off it.
#
# Three of the four statistics -- the likelihood ratio, Rao's score and
# Terrell's gradient -- are functions of the model refitted with one
# coefficient HELD at a value. Wald's is not: it reads the unrestricted fit
# alone, which is exactly why it is the one that is not invariant to how the
# parameter is written.


#' Refit With One Coefficient Held at a Value
#'
#' @description
#' Maximizes the penalized likelihood over every coefficient but one, which
#' is held at `value`. It returns the restricted maximum, the coefficients
#' that attain it and the score there.
#'
#' @details
#' The hold is written on the specification and enforced in the solve, so the
#' refit does not rewrite the formula and does not rebuild any term: a design
#' whose bases were rebuilt would carry different knots and would correspond
#' to a different model. See
#' [held_positions()] and [iwls_solve()] for the two halves.
#'
#' # The hyperparameters are held too
#'
#' The refit runs the alternation once at the hyperparameters of the
#' unrestricted fit, with no outer search. The statistics built on it are
#' therefore conditional on the smoothing that the data chose once. That is
#' the reading that a penalized fit supports, because the penalty is a
#' log-prior at the hyperparameters in force, and it keeps an interval by
#' inversion affordable: re-selecting the smoothing at every held value would
#' put a whole outer search inside every step of a root find. A profiled
#' version requires refitting by hand.
#'
#' # Coefficients that cannot be held
#'
#' A coefficient under a kinked penalty (a lasso, a SCAD, an MCP) cannot be
#' held, for two reasons. The block is fitted by a coordinate descent, which
#' updates every coordinate of the block and cannot leave one at a held
#' value, so that [coord_fit()] signals an error when a held coordinate lies
#' in its block. And at the kink the objective has no curvature for the
#' statistics to read, while the null distribution of a coefficient that a
#' selection kept is a different problem. An aliased coefficient has no
#' estimate to test. A model that carries a structural term cannot hold a
#' coefficient either, because the contribution of the term is the state of a
#' recursion and not \eqn{X\beta}, and the hold does not reach it.
#'
#' A coefficient under a penalty that is twice differentiable is held like
#' any other: a smooth's linear column and its rotated coordinates, a ridge,
#' a random effect. [statmod_stat_at()] states the meaning of such a fit.
#'
#' @param fit A [StatmodFit()].
#' @param param The distribution parameter whose equation carries the
#'   coefficient, a single string.
#' @param coefname The coefficient's name, a single string, as
#'   [coef.StatmodFit()] reports it.
#' @param value The value at which to hold it, a single number.
#'
#' @return A list:
#'   \describe{
#'     \item{`coefficients`}{the restricted estimates, a named list.}
#'     \item{`loglik`}{the log-likelihood there, unpenalized.}
#'     \item{`objective`}{the penalized objective there, the negative
#'       log-likelihood plus the penalty, which the fit minimizes.}
#'     \item{`score`}{the penalized score at the restricted point, a numeric
#'       vector over the stacked coefficients. Its free entries are zero to
#'       the tolerance; the statistics read the held one.}
#'     \item{`at`}{the held coordinate's position in that vector.}
#'     \item{`par`}{the restricted estimates, stacked.}
#'     \item{`information`}{a function of no arguments returning the
#'       penalized information at the restricted point. It is read at most
#'       once however often it is called, and it is not read if the function
#'       is never called, which is the ordinary case: the likelihood ratio and
#'       the gradient statistic need no matrix.}
#'     \item{`mode_error`}{a function of no arguments returning
#'       [restricted_mode_error()], how far above its own mode the refit
#'       stopped. It reads the information, so it costs one Hessian.}
#'     \item{`labels`}{the stacked coefficient labels.}
#'     \item{`converged`}{the inner optimizer's flag, a single logical. It
#'       records whether a stopping rule fired and does not measure the
#'       distance to the mode, which [restricted_mode_error()] measures.}
#'   }
#'
#' @seealso [iwls_solve()], which drops the held coordinate from the system.
#'
#' @keywords internal
statmod_restrict <- function(fit, param, coefname, value) {
  # a value of another length reaches the hold as a named vector of that
  # length and fails several frames down, inside the objective, on a names
  # assignment -- which names neither the argument nor the mistake
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("'value' must be a single finite number.", call. = FALSE)
  }
  spec0 <- fit@spec
  params <- spec0@distrib@params
  if (!param %in% params) {
    stop(sprintf("'%s' is not a parameter of %s. Its parameters are %s.",
                 param, spec0@distrib@distrib_name,
                 paste(params, collapse = ", ")), call. = FALSE)
  }
  design0 <- statmod_design(spec0)
  if (length(attr(design0, "structural"))) {
    stop("A model carrying a structural term cannot hold a coefficient: ",
         "its\n  contribution is a recursion's state and not a column times ",
         "a coefficient.", call. = FALSE)
  }
  nms <- design0[[param]]$coef_names
  j <- match(coefname, nms)
  if (is.na(j)) {
    stop(sprintf("'%s' names no coefficient of the '%s' equation.",
                 coefname, param), call. = FALSE)
  }
  npar <- vapply(design0, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  at <- offs[[match(param, params)]] + j
  if (aliased_labels(spec0, design0, at) %in% fit@aliased) {
    stop(sprintf(paste0("'%s' is aliased in the '%s' equation, so it has no ",
                        "estimate to test:\n  its column is a combination ",
                        "of the others and the fit left it out."),
                 coefname, param), call. = FALSE)
  }
  if (at %in% kinked_coords(spec0, design0)) {
    stop(sprintf(paste0("'%s' is under a kinked penalty (a lasso, a SCAD ",
                        "or an MCP) and\n  cannot be held. That block is ",
                        "fitted by a coordinate descent, which\n  updates ",
                        "every coordinate of the block and cannot leave one ",
                        "at a held value,\n  and at the kink the objective ",
                        "has no curvature for a statistic to read.\n  A ",
                        "penalty that is twice differentiable (a ridge, a ",
                        "random effect,\n  a smooth) is held like any other ",
                        "coefficient."), coefname), call. = FALSE)
  }

  hc <- spec0@held_coef
  hc[[param]] <- stats::setNames(as.numeric(value), coefname)
  spec <- S7::set_props(spec0, held_coef = hc)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$smooth
  cfg <- inner_settings(method)
  beta <- unlist(fit@coefficients[params], use.names = FALSE)
  beta[at] <- as.numeric(value)

  res <- statmod_alternate(spec, design, blocks, fit@hyper, method, beta,
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0))
  cf <- res$obj$split(res$par)
  # the objective is a MINIMAND -- the negative log-likelihood plus the
  # penalty -- so its gradient carries the opposite sign to the score
  sc <- -res$obj$gr(res$par)
  # the Hessian is read at most ONCE, and only where a caller asks: the
  # likelihood ratio and the gradient statistic never do. Measured on a
  # Poisson of twenty-six coefficients it is 63 per cent of the refit
  # (0.0044 s against 0.0070); at two coefficients, 2.6 per cent.
  hess <- NULL
  info <- function() {
    # the objective's Hessian IS the penalized information there, the
    # objective being the negative log-likelihood plus the penalty
    if (is.null(hess)) hess <<- as_dense(res$obj$he(res$par))
    hess
  }
  list(coefficients = cf,
       loglik = statmod_loglik_at(spec, cf, design),
       objective = res$value,
       score = sc,
       at = at,
       par = res$par,
       information = info,
       mode_error = function() restricted_mode_error(sc, info(), at),
       labels = rownames(coef_labels(spec, design)),
       converged = isTRUE(res$converged))
}


#' How Far Above Its Mode a Restricted Fit Stopped
#'
#' @description
#' Computes \eqn{\tfrac12 g'K^{-1}g} over the free coordinates of a
#' restricted fit, which is the decrease of the penalized objective that a
#' Newton step from the point returned by [statmod_restrict()] would predict,
#' in log-likelihood units.
#'
#' @details
#' [inner_mode_error()] computes the same quantity inside a search; here it
#' is applied to a restricted fit. The inner optimizer's flag records whether a
#' stopping rule fired, which is a boolean about a threshold on a score whose
#' size depends on the model, whereas whether a point is at its mode is a
#' matter of distance and has a natural scale. A stopping rule reads a change
#' in the objective and this quantity reads the point, so the two need not
#' agree, and the print method of [StatmodTest()] reads the mode error and
#' not the flag.
#'
#' The held coordinate is excluded because it is not free: its score is the
#' quantity that Rao's statistic reads, and it does not vanish under the
#' restriction.
#'
#' @param score The restricted fit's score, all coordinates.
#' @param information Its penalized information, all coordinates.
#' @param at The position of the held coordinate.
#'
#' @return A single number, or `NA` where the score on the free coordinates
#'   is not finite, where there is no free coordinate, or where the
#'   information on the free coordinates cannot be inverted.
#'
#' @seealso [mode_error_limit()], the limit against which it is read, and
#'   [inner_mode_error()], the same reading inside a search.
#'
#' @keywords internal
restricted_mode_error <- function(score, information, at) {
  free <- setdiff(seq_along(score), at)
  if (!length(free) || !all(is.finite(score[free]))) return(NA_real_)
  g <- score[free]
  v <- tryCatch(0.5 * sum(g * solve(information[free, free, drop = FALSE], g)),
                error = function(e) NA_real_)
  if (isTRUE(is.finite(v))) v else NA_real_
}


#' Which Stacked Coordinates a Kinked Penalty Covers
#'
#' @description
#' Returns the positions of the coefficients whose penalty is not
#' differentiable everywhere -- a lasso, a SCAD, an MCP -- as an integer
#' vector.
#'
#' @details
#' The function tests the coordinates in the way that [statmod_blocks()] tests
#' the blocks to decide which of the two fitting schemes is used for a block,
#' and it reads the same predicate ([penalty_has_kink()]), so that the two
#' cannot disagree. A penalty that is twice differentiable contributes
#' nothing, so that a ridge, a random effect and a smooth are testable.
#'
#' The result follows from the enumeration of each penalty, so a partial
#' penalty (one covering some of a term's coordinates and not the rest) names
#' only what it covers, and a term is not classified as a whole. The penalty of
#' a structural term covers positions among the term's own parameters and not
#' among the design's, and contributes nothing here.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#'
#' @return An integer vector of stacked positions, possibly empty.
#'
#' @seealso [statmod_penalized()], the enumeration that it reads, and
#'   [statmod_blocks()], which gives the same classification for a whole
#'   block.
#'
#' @keywords internal
kinked_coords <- function(spec, design) {
  out <- integer(0)
  for (u in statmod_penalized(spec, design)) {
    if (isTRUE(u$structural) || is.null(u$index)) next
    kink <- tryCatch(penalty_has_kink(u$penalty), error = function(e) TRUE)
    if (isTRUE(kink)) out <- c(out, as.integer(u$index))
  }
  sort(unique(out))
}


#' One Likelihood Statistic for One Coefficient
#'
#' @description
#' Computes the Wald, likelihood-ratio, Rao score or Terrell gradient
#' statistic for the hypothesis that one coefficient equals `value`, with its
#' p-value.
#'
#' @details
#' All four read the penalized objective, which is the function that the fit
#' maximizes. Writing \eqn{\ell_p(\beta) = \ell(\beta) - \rho(\beta;\theta)}
#' for it, \eqn{\hat\beta} for the unrestricted estimates, \eqn{\tilde\beta}
#' for the estimates in which \eqn{\beta_j} is held at \eqn{b},
#' \eqn{U = \partial\ell_p} for its score and \eqn{K} for the penalized
#' information, the four are
#' \deqn{W = (\hat\beta_j - b)^2 / \widehat{\mathrm{Var}}(\hat\beta_j),
#' \qquad LR = 2\{\ell_p(\hat\beta) - \ell_p(\tilde\beta)\},}
#' \deqn{S = U_j(\tilde\beta)^2 \, [K(\tilde\beta)^{-1}]_{jj}, \qquad
#' T = U_j(\tilde\beta)\,(\hat\beta_j - b),}
#' each compared with a \eqn{\chi^2_1}. Where the model carries no penalty
#' \eqn{\ell_p} is the log-likelihood and these are the classical four.
#'
#' The score statistic has the general form \eqn{U'K^{-1}U}. It reduces to
#' the product above because every other component of \eqn{U} vanishes at the
#' restricted maximum, which also indicates that the restricted fit has
#' converged.
#'
#' # The penalized objective and the likelihood
#'
#' The likelihood ratio is taken on \eqn{\ell_p} and not on \eqn{\ell}.
#' \eqn{\tilde\beta} maximizes \eqn{\ell_p} under the restriction and not
#' \eqn{\ell}, so \eqn{\ell(\tilde\beta)} can exceed
#' \eqn{\ell(\hat\beta)} and the difference of log-likelihoods can be
#' negative. The other three statistics read the penalized objective as well:
#' \eqn{U} and \eqn{K} are its gradient and its curvature, and the Wald
#' variance is \eqn{(H+S)^{-1}} under `vcov_type = bayesian_vcov()`.
#'
#' For a shrunk coordinate this is the reading that [vcov.StatmodFit()] gives
#' under `vcov_type = bayesian_vcov()`: the penalty is a log-prior, so \eqn{\ell_p} is
#' a log posterior at the hyperparameters in force, and the interval by
#' inversion is the set that the posterior does not reject, conditional on
#' them. The interpretation that the Wald interval on such a row carries
#' therefore extends to the other three statistics. On a smooth the inverted
#' interval and the Wald one agree to second order, and the inversion adds
#' the asymmetry.
#'
#' # Differences between the statistics
#'
#' The Wald statistic reads only the unrestricted fit and needs no refit,
#' which is why the summary reports it by default. It depends on how the
#' parameter is written: a nonlinear reparametrization of the coefficient
#' changes it, whereas the likelihood ratio, a difference of objective values
#' at the same two points, is not changed. The Wald statistic also loses power
#' against a distant alternative (the Hauck-Donner effect), because the
#' curvature is read at \eqn{\hat\beta} rather than under the null.
#'
#' The gradient statistic needs the restricted fit and no matrix inversion:
#' it is the score contracted with the distance between the two fits.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()].
#' @param value The value under the null, a single number. The default is
#'   `0`.
#' @param test A [test_statistic] object: [wald_test()], [lr_test()],
#'   [score_test()] or [gradient_test()].
#' @param vcov_type A [vcov_method] object, the variance that the Wald statistic
#'   reads, as [vcov.StatmodFit()] takes it.
#' @param mode_error Whether to read how far above its mode the restricted
#'   fit stopped ([restricted_mode_error()]). The default is `FALSE`, because
#'   reading it costs one Hessian where the likelihood ratio and the gradient
#'   statistic need none. [summary.StatmodFit()], which calls this once per
#'   row, and [statmod_invert()], which calls it once per evaluation of the
#'   statistic, do not read it.
#' @param objective0 The penalized objective at the unrestricted joint mode,
#'   as [joint_objective()] returns it, for a caller that reads the
#'   likelihood ratio many times on one fit. `NULL` computes it here.
#'
#' @return A list with `test` (the [test_statistic] object), `statistic`, `df`,
#'   `p.value`, `converged` --
#'   `NA` for the Wald statistic and the restricted fit's flag for the other
#'   three -- and `mode_error`, `NA` unless it is requested.
#'
#' @seealso [statmod_restrict()], which the three restricted statistics read.
#'
#' @keywords internal
statmod_stat_at <- function(fit, param, coefname, value = 0,
                         test = wald_test(), vcov_type = bayesian_vcov(),
                         mode_error = FALSE, objective0 = NULL) {
  test <- as_test_statistic(test)
  vcov_type <- as_vcov_method(vcov_type)
  tname <- test@name
  spec <- fit@spec
  design <- statmod_design(spec)
  params <- spec@distrib@params
  nms <- design[[param]]$coef_names
  j <- match(coefname, nms)
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  at <- if (is.na(j)) NA_integer_ else offs[[match(param, params)]] + j
  bhat <- if (is.na(j)) NA_real_ else fit@coefficients[[param]][[j]]

  out <- function(s, conv, me = NA_real_) {
    list(test = test, statistic = s, df = 1L,
         p.value = stats::pchisq(s, 1L, lower.tail = FALSE),
         converged = conv, mode_error = me)
  }

  if (!test_restricted(test)) {
    # the only one that needs no refit, and the only one a reparametrization
    # moves. Nothing was refitted, so there is no mode to be above
    V <- stats::vcov(fit, vcov_type = vcov_type)
    lab <- rownames(coef_labels(spec, design))[at]
    v <- if (is.na(at) || !lab %in% rownames(V)) NA_real_ else V[lab, lab]
    if (!isTRUE(is.finite(v)) || v <= 0) return(out(NA_real_, NA))
    return(out((bhat - value)^2 / v, NA))
  }

  r <- statmod_restrict(fit, param, coefname, value)
  me <- if (isTRUE(mode_error)) r$mode_error() else NA_real_
  if (identical(tname, "lr")) {
    # 2 times the rise in the PENALIZED OBJECTIVE, which is the function both
    # fits maximize, so the difference is non-negative by construction. Read
    # on the log-likelihood alone it is not: the restricted point maximizes
    # the penalized objective, not the likelihood, so it can sit ABOVE the
    # unrestricted point there -- measured at -0.0379 and -0.1075 on a
    # smooth's own coordinate. With no penalty the two are the same number.
    # The unrestricted side is the JOINT mode, which is what the restricted
    # refit reaches: where the criterion estimated some coefficients -- a
    # dispersion by REML, the default -- the fitted point is not that mode,
    # and differencing against it subtracted 2 (l_ML - l_REML), 0.118 on a
    # linear regression at n = 40, negative at the estimate itself.
    o0 <- if (is.null(objective0)) joint_objective(fit) else objective0
    return(out(2 * (r$objective - o0), r$converged, me))
  }
  u <- r$score[[r$at]]
  if (identical(tname, "gradient")) {
    return(out(u * (bhat - value), r$converged, me))
  }
  K <- r$information()
  V <- tryCatch(solve_pd(K, "the penalized information at the restricted fit",
                         r$labels),
                error = function(e) NULL)
  if (is.null(V)) return(out(NA_real_, r$converged, me))
  out(u^2 * V[r$at, r$at], r$converged, me)
}


#' The Penalized Objective at the Unrestricted Joint Mode
#'
#' @description
#' Computes the minimum of the penalized objective over every coefficient at
#' the fitted hyperparameters, which is the unrestricted side of the
#' likelihood ratio that [statmod_stat_at()] computes.
#'
#' @details
#' The restricted fit of [statmod_restrict()] holds one coefficient and takes
#' every other one to the joint mode of the penalized likelihood. The
#' unrestricted side has to be read at the same kind of point, or the
#' statistic compares two fits whose other coefficients were chosen by
#' different rules. Where the outer criterion estimated some coefficients
#' itself (by default [reml()] and [ml()] estimate on the criterion the
#' unpenalized coefficients of every parameter but the first, see
#' [marginal_coords()]), the point that the fit returned is not the joint
#' mode, and `fit@objective` is not its minimum. A likelihood ratio
#' differenced against `fit@objective` then contains the difference between
#' the objective at the fitted point and at the joint mode.
#'
#' The full model is therefore refitted from the fitted coefficients, which
#' are already close to the joint mode, at the fitted hyperparameters, and no
#' coefficient is held. Where the criterion estimated nothing, the fit is at the
#' joint mode and `fit@objective` is returned unchanged.
#'
#' @param fit A [StatmodFit()].
#'
#' @return A single number, on the scale of `fit@objective`.
#'
#' @seealso [statmod_stat_at()], [statmod_restrict()].
#'
#' @keywords internal
joint_objective <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  gam <- tryCatch(marginal_coords(spec, design, fit@methods$outer),
                  error = function(e) marginal_coords(spec, design, NULL))
  if (!length(gam$where)) return(fit@objective)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$smooth
  cfg <- inner_settings(method)
  beta <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
  res <- statmod_alternate(spec, design, blocks, fit@hyper, method, beta,
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0))
  res$value
}


#' A Test of One Coefficient Against One Value
#'
#' @description
#' Holds the result of testing that one coefficient of a fit equals one
#' value: the statistic, its distribution, the p-value and what was tested.
#'
#' @details
#' The four statistics are compared with the same \eqn{\chi^2_1}, so `df` is
#' always 1 and the object holds one kind of quantity whatever produced it.
#'
#' # Whether the refit reached its mode
#'
#' `mode_error` and `converged` both describe the restricted refit and both
#' are `NA` for the Wald statistic, which reads the unrestricted fit and
#' refits nothing. They measure different things, and only the first is
#' printed.
#'
#' `converged` is the inner optimizer's flag: whether a stopping rule fired,
#' a boolean about a threshold on a score whose size depends on the model.
#' `mode_error` is [restricted_mode_error()], the decrease of the penalized
#' objective that a Newton step from the point returned by the refit would
#' predict, in log-likelihood units, judged against [mode_error_limit()].
#' The two need not agree, because a stopping rule reads a change in the
#' objective whereas the mode error reads the distance to the mode. See
#' [restricted_mode_error()].
#'
#' A statistic computed from a refit that stopped short is the value at the
#' point where the refit stopped. `print()` reports this, reading
#' `mode_error`.
#'
#' The class declares its properties, as every object that this toolkit
#' returns does, whereas the `htest` returned by the tests of base R is a list
#' with a class attribute and no declared structure, so that it cannot be
#' validated.
#'
#' @param test The [test_statistic] object that was computed: [wald_test()],
#'   [lr_test()], [score_test()] or [gradient_test()]. Its `label` is the name
#'   under which `print()` shows the statistic.
#' @param statistic Its value, compared with a \eqn{\chi^2_1}.
#' @param df One, always.
#' @param p.value The upper tail of that \eqn{\chi^2_1}.
#' @param estimate The unrestricted estimate of the coefficient.
#' @param null_value The value against which the coefficient was tested.
#' @param parameter,coefficient Which coefficient of which equation.
#' @param mode_error How far above its own mode the restricted refit
#'   stopped, in log-likelihood units; `NA` for the Wald statistic.
#' @param converged The inner optimizer's flag at the restricted refit;
#'   `NA` for the Wald statistic.
#'
#' @return An object of class `StatmodTest`.
#'
#' @seealso [statmod_test()], which builds one, and
#'   [restricted_mode_error()], which `mode_error` reports.
#'
#' @examples
#' dd <- data.frame(x = seq(-1, 1, length.out = 40))
#' dd$y <- rpois(40, exp(0.3 + dd$x))
#' fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
#' S7::S7_inherits(statmod_test(fit, param = "mu", coefname = "x", value = 1),
#'                 StatmodTest)
#'
#' @name StatmodTest-class
#' @aliases StatmodTest
#' @keywords internal
#' @export
StatmodTest <- S7::new_class("StatmodTest",
  properties = list(
    # which of the four statistics, as the object `statmod_test()` received;
    # print() reads its label
    test = test_statistic,
    statistic = S7::class_numeric,
    df = S7::class_integer,
    p.value = S7::class_numeric,
    # the unrestricted estimate and the value it was tested against, so the
    # object says what the null WAS rather than leaving it to the call
    estimate = S7::class_numeric,
    null_value = S7::class_numeric,
    parameter = S7::class_character,
    coefficient = S7::class_character,
    # how far above its own mode the restricted refit stopped, and the raw
    # flag beside it. The first is what print() reads; the second is kept
    # because it is what the optimizer reported and a reader may want it
    mode_error = S7::class_numeric,
    converged = S7::class_logical
  )
)


#' Print a Test
#'
#' @description
#' Prints the statistic, its degrees of freedom and the p-value on one line,
#' under a title naming the test, in the layout that R uses for a hypothesis
#' test. A note follows when the restricted refit stopped above its mode, by
#' more than [mode_error_limit()], or when its mode error could not be read.
#'
#' @param x A [StatmodTest()].
#' @param digits Significant digits.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#'
#' @keywords internal
print.StatmodTest <- function(x, digits = 4L, ...) {
  cat("\n\t", x@test@label, " test on a statmod fit\n\n", sep = "")
  cat("coefficient:  ", x@parameter, ":", x@coefficient, "\n", sep = "")
  cat("X-squared = ", format(x@statistic, digits = digits),
      ",  df = ", x@df,
      ",  p-value ", format.pval(x@p.value, digits = digits, eps = 1e-16),
      "\n", sep = "")
  cat("alternative hypothesis: true ", x@coefficient, " is not equal to ",
      format(x@null_value, digits = digits), "\n", sep = "")
  cat("estimate: ", format(x@estimate, digits = digits), "\n", sep = "")
  # a statistic read off a refit that stopped short is the value where it
  # stopped, and saying so is the difference between a number and a result.
  # What is read is the MODE ERROR and not the optimizer's flag. Measured,
  # the flag rejects exactly the best-located points -- where the refit lands
  # on its mode in one step the objective does not move and the stall guard
  # fires -- so a warning printed off it fires loudest where the refit is
  # best. See restricted_mode_error().
  # Wald refits nothing, so there is no mode for it to be above, and
  # `converged` is NA exactly there
  refitted <- length(x@converged) == 1L && !is.na(x@converged)
  me <- if (length(x@mode_error) == 1L) x@mode_error else NA_real_
  if (refitted && is.na(me)) {
    cat("\nthe mode error of the restricted fit could not be read,\n",
        "  so the statistic is the value at wherever the refit stopped\n",
        sep = "")
  } else if (refitted && me > mode_error_limit()) {
    cat("\nthe restricted fit stopped above its own mode",
        sprintf(" (by %s, against %s)", format(me, digits = digits),
                format(mode_error_limit(), digits = digits)),
        ",\n  so the statistic is the value at wherever it stopped\n",
        sep = "")
  }
  invisible(x)
}
S7::method(print, StatmodTest) <- print.StatmodTest


#' Test One Coefficient Against One Value
#'
#' @description
#' Tests that one coefficient of a fit equals one value, by Wald's
#' statistic, the likelihood ratio, Rao's score or Terrell's gradient.
#'
#' @details
#' `summary(fit, test =)` reports one of the four statistics for the rows of
#' the coefficient tables that can be held (not a row under a kinked penalty,
#' not an aliased one, and none for a model with a structural term), always
#' against zero. This function tests one coefficient against any value, for
#' questions that a summary cannot answer, such as whether a slope is one,
#' whether an elasticity is unity, or whether a coefficient matches a value
#' fixed outside the data.
#'
#' The statistics, the quantities that they read and their cost are described
#' at [statmod_stat_at()]. The last three are built on the restricted fit of
#' [statmod_restrict()], and inverting one gives the interval of
#' [confint.StatmodFit()].
#'
#' Reading one row costs one restricted refit, whereas `summary(test =)`
#' costs one per row. The function also reads the mode error of the refit,
#' which costs one Hessian of the restricted fit; the Rao score statistic
#' reads that matrix anyway, so that it adds no cost there. The summary and
#' the interval do not read the mode error.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()]:
#'   the distribution parameter and the coefficient's name in its equation.
#' @param value The value under the null, a single number. The default is
#'   `0`, which is the value that a summary reports.
#' @param test A [test_statistic] object: [wald_test()] (the default),
#'   [lr_test()], [score_test()] or [gradient_test()]. The first reads the
#'   unrestricted fit only; the other three cost one restricted refit. The
#'   four are compared at [test_statistic]. A character string is rejected with
#'   an error.
#' @param vcov_type A [vcov_method] object, the variance that the Wald statistic
#'   reads, passed to [vcov.StatmodFit()]: [bayesian_vcov()] (the default),
#'   [frequentist_vcov()] or [unconditional_vcov()]. The conventions are
#'   compared at [vcov_method].
#'
#' @return A [StatmodTest()].
#'
#' @seealso [test_statistic] for the four statistics, [vcov_method] for the
#'   variance conventions, [summary.StatmodFit()] for the rows that can be
#'   held, against zero, [confint.StatmodFit()] for the interval that a test
#'   inverts.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200))
#' dd$y <- rpois(200, exp(0.3 + dd$x))
#' fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
#' ## the summary tests whether the slope is zero; this tests whether it is one
#' statmod_test(fit, param = "mu", coefname = "x", value = 1, test = lr_test())
#' statmod_test(fit, param = "mu", coefname = "x", value = 0, test = lr_test())
#'
#' @export
statmod_test <- function(fit, param, coefname, value = 0, test = wald_test(),
                         vcov_type = bayesian_vcov()) {
  test <- as_test_statistic(test)
  vcov_type <- as_vcov_method(vcov_type)
  # the mode error is asked for HERE and nowhere else: this is the one
  # surface a reader reads one row of, so one extra Hessian is affordable,
  # where the two loops that call the same helper in quantity are not
  s <- statmod_stat_at(fit, param, coefname, value, test, vcov_type,
                       mode_error = TRUE)
  design <- statmod_design(fit@spec)
  nms <- design[[param]]$coef_names
  j <- match(coefname, nms)
  est <- if (is.na(j)) NA_real_ else fit@coefficients[[param]][[j]]
  StatmodTest(test = s$test, statistic = s$statistic, df = s$df,
              p.value = s$p.value, estimate = est,
              null_value = as.numeric(value), parameter = param,
              coefficient = coefname, mode_error = s$mode_error,
              converged = s$converged)
}


#' A Confidence Interval by Inverting a Test
#'
#' @description
#' Returns the set of values that a test does not reject at `level`, found by
#' locating the two points where the statistic crosses its critical value.
#'
#' @details
#' An interval by inversion is \eqn{\{b : T(b) \le \chi^2_{1,\alpha}\}}, and
#' it is not symmetric about the estimate unless the statistic is a parabola
#' in \eqn{b}. The Wald statistic is a parabola by construction, because the
#' curvature is read once at \eqn{\hat\beta}, and the other three statistics
#' are not. The asymmetry is the reason to invert one of the other
#' three: the interval follows the shape of the likelihood and not a quadratic
#' fitted at its maximum.
#'
#' The search starts at the half-width of the Wald interval and widens the
#' distance from the estimate by a factor of 1.6, at most `maxit` times on
#' each side, until the statistic exceeds its critical value. It then locates
#' the crossing with [stats::uniroot()] (Brent's method) at tolerance 1e-7.
#' Each evaluation is one restricted refit, warm-started at the unrestricted
#' estimates and run at the fitted hyperparameters, so the interval requires
#' no outer search.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()].
#' @param level The confidence level, a single number in (0, 1).
#' @param test The [test_statistic] to invert: [lr_test()] (the default),
#'   [score_test()], [gradient_test()] or [wald_test()].
#' @param vcov_type The [vcov_method] that sets the starting bracket, and that the
#'   Wald statistic reads.
#' @param maxit How many widenings to try on each side before the search
#'   stops.
#'
#' @return A numeric vector of two, `lower` and `upper`. An end that the
#'   search could not bracket within `maxit` widenings, or at which the
#'   statistic could not be evaluated, is `NA`; both ends are `NA` when the
#'   standard error of the coefficient is not finite and positive.
#'
#' @seealso [statmod_stat_at()], the statistic that is inverted.
#'
#' @keywords internal
statmod_invert <- function(fit, param, coefname, level = 0.95,
                           test = lr_test(), vcov_type = bayesian_vcov(),
                           maxit = 40L) {
  test <- as_test_statistic(test)
  vcov_type <- as_vcov_method(vcov_type)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  q <- stats::qchisq(level, 1L)
  spec <- fit@spec
  design <- statmod_design(spec)
  params <- spec@distrib@params
  j <- match(coefname, design[[param]]$coef_names)
  if (is.na(j)) {
    stop(sprintf("'%s' names no coefficient of the '%s' equation.",
                 coefname, param), call. = FALSE)
  }
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  at <- offs[[match(param, params)]] + j
  bhat <- fit@coefficients[[param]][[j]]
  V <- stats::vcov(fit, vcov_type = vcov_type, readable = FALSE)
  lab <- rownames(coef_labels(spec, design))[at]
  se <- if (lab %in% rownames(V)) sqrt(V[lab, lab]) else NA_real_
  if (!isTRUE(is.finite(se)) || se <= 0) return(c(lower = NA_real_,
                                                  upper = NA_real_))

  o0 <- if (identical(test@name, "lr")) joint_objective(fit) else NULL
  gap <- function(b) {
    s <- statmod_stat_at(fit, param, coefname, b, test, vcov_type,
                         objective0 = o0)$statistic
    if (!isTRUE(is.finite(s))) return(NA_real_)
    s - q
  }
  side <- function(dir) {
    step <- sqrt(q) * se
    a <- bhat
    for (k in seq_len(as.integer(maxit))) {
      b <- bhat + dir * step
      v <- gap(b)
      if (is.na(v)) return(NA_real_)
      if (v > 0) {
        r <- tryCatch(stats::uniroot(gap, interval = sort(c(a, b)),
                                     tol = 1e-7)$root,
                      error = function(e) NA_real_)
        return(r)
      }
      a <- b
      step <- step * 1.6
    }
    NA_real_
  }
  c(lower = side(-1), upper = side(1))
}


#' Which Coefficients a Restricted Fit Can Hold
#'
#' @description
#' Returns the stacked positions of the coefficients for which the three
#' restricted statistics are defined: every coordinate less the ones under a
#' kinked penalty, less the aliased ones, and none at all where the model
#' carries a structural term.
#'
#' @details
#' It applies to the whole vector the enumeration that [statmod_restrict()]
#' applies to one coordinate at a time, signaling an error for a coordinate
#' that cannot be held. A caller reporting a table needs it in this form,
#' since each restricted statistic is a refit and an interval by inversion
#' requires many, so that a row that cannot be tested has to be recognized
#' before its refit is computed and not after.
#'
#' A coordinate under a penalty that is twice differentiable is included:
#' a smooth's linear column and its rotated coordinates, a ridge, a random
#' effect. The restricted fit holds it exactly as it holds an unpenalized
#' one, because [fit_smooth()] enforces the hold in the solve. Only a kinked
#' penalty is excluded, for the reasons that [statmod_restrict()] gives.
#'
#' @param fit A [StatmodFit()].
#' @param spec,design The specification and its design, where the caller has
#'   them already.
#'
#' @return An integer vector of stacked positions, possibly empty, named by
#'   the labels that [coef_labels()] gives them.
#'
#' @seealso [statmod_restrict()], which rejects the same coordinates one at a
#'   time, and [kinked_coords()].
#'
#' @keywords internal
testable_coords <- function(fit, spec = fit@spec,
                            design = statmod_design(spec)) {
  none <- stats::setNames(integer(0), character(0))
  if (length(attr(design, "structural"))) return(none)
  lab <- rownames(coef_labels(spec, design))
  free <- setdiff(seq_along(lab), kinked_coords(spec, design))
  if (!length(free)) return(none)
  keep <- free[!lab[free] %in% fit@aliased]
  if (!length(keep)) return(none)
  stats::setNames(as.integer(keep), lab[keep])
}
