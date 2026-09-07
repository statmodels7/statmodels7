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
#' is held at `value`. What comes back is the restricted maximum, the
#' coefficients that attain it and the score there.
#'
#' @details
#' The hold is written on the specification and enforced in the solve, so no
#' formula is rewritten and no term is rebuilt: a design whose bases were
#' rebuilt would carry different knots and answer for a different model. See
#' [held_positions()] and [iwls_solve()] for the two halves.
#'
#' # The hyperparameters are held too
#'
#' The refit runs the alternation ONCE at the hyperparameters the
#' unrestricted fit arrived at, with no outer search. That makes the
#' statistics built on it CONDITIONAL on the smoothing the data chose once,
#' which is the honest reading of a penalized fit and is also what makes an
#' interval by inversion affordable -- re-selecting the smoothing at every
#' held value would put a whole outer search inside every step of a root
#' find. A caller who wants the profiled version can refit by hand.
#'
#' # What is refused
#'
#' A coefficient under a KINKED penalty -- a lasso, a SCAD, an MCP -- for two
#' reasons that point the same way: that block is fitted by a coordinate
#' descent, and [coord_fit()] does not read `held_coef` at all, so a hold
#' there would be moved by the sweep without a word (measured: held at 3, the
#' refit reports 2.7348); and at the kink the objective has no curvature, so
#' the statistics have nothing to read and the null distribution of a
#' coefficient a selection kept is a different problem. An ALIASED
#' coefficient, which has no estimate to test. And a model carrying a
#' structural term, whose contribution is a recursion's state rather than
#' \eqn{X\beta} and which the hold does not reach.
#'
#' A coefficient under a penalty that IS twice differentiable is held like
#' any other: a smooth's linear column and its rotated coordinates, a ridge,
#' a random effect. What such a fit means is stated at [statmod_stat_at()].
#'
#' @param fit A [StatmodFit()].
#' @param param The distribution parameter whose equation carries the
#'   coefficient, a single string.
#' @param coefname The coefficient's name, a single string, as
#'   [coef.StatmodFit()] reports it.
#' @param value The value to hold it at, a single number.
#'
#' @return A list:
#'   \describe{
#'     \item{`coefficients`}{the restricted estimates, a named list.}
#'     \item{`loglik`}{the log-likelihood there, unpenalized.}
#'     \item{`objective`}{the penalized objective there.}
#'     \item{`score`}{the penalized score at the restricted point, a numeric
#'       vector over the stacked coefficients. Its free entries are zero to
#'       the tolerance; the held one is what the statistics read.}
#'     \item{`at`}{the held coordinate's position in that vector.}
#'     \item{`par`}{the restricted estimates, stacked.}
#'     \item{`information`}{a function of no arguments returning the
#'       penalized information at the restricted point. It is read at most
#'       once however often it is called, and not at all where nothing calls
#'       it, which is the ordinary case: the likelihood ratio and the
#'       gradient statistic need no matrix.}
#'     \item{`mode_error`}{a function of no arguments returning
#'       [restricted_mode_error()], how far above its own mode the refit
#'       stopped. It reads the information, so it costs one Hessian.}
#'     \item{`labels`}{the stacked coefficient labels.}
#'     \item{`converged`}{the inner optimizer's flag, a single logical.
#'       Measured, it is anti-correlated with how well the point is located
#'       -- see [restricted_mode_error()], which is the reading to prefer.}
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
    stop("'value' is the single finite number the coefficient is held at.",
         call. = FALSE)
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
    stop(sprintf(paste0("'%s' is under a kinked penalty -- a lasso, a SCAD ",
                        "or an MCP -- and\n  cannot be held. That block is ",
                        "fitted by a coordinate descent, which does\n  not ",
                        "read the hold, and at the kink the objective has ",
                        "no curvature for\n  a statistic to read. A penalty ",
                        "that is twice differentiable -- a ridge, a\n  ",
                        "random effect, a smooth -- is held like any other ",
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
#' \eqn{\tfrac12 g'K^{-1}g} over the FREE coordinates of a restricted fit:
#' how much of the penalized objective is still on the table at the point
#' [statmod_restrict()] returned, in log-likelihood units.
#'
#' @details
#' It is [inner_mode_error()]'s question asked of a restricted fit, and it
#' is asked for the reason statmodels7 0.81.0 established: the inner
#' optimizer's flag says whether a stopping rule fired, which is a boolean
#' about a threshold on a score whose size depends on the model, while
#' whether a point is at its mode is a matter of distance and has a natural
#' scale.
#'
#' The two come apart here as they do there, and they come apart BACKWARDS:
#' the flag rejects the points that are located best. Measured on `y ~ x`
#' with a Poisson response at \eqn{n = 200}, holding the slope at nineteen
#' values from 0.6 to 1.5, all nineteen are at their mode -- the worst at
#' 1.04e-11 against a limit of 1e-03 -- and fifteen report the flag. The four
#' it rejects are **the four best-located of the nineteen**, every one of
#' them near 1e-23 where the fifteen it accepts run out to 1.04e-11. On a
#' gaussian smooth holding `s(x).lin` at thirteen values it is the two
#' best-located of the thirteen, at 1.31e-12 against a worst of 6.67e-09.
#'
#' The mechanism is that a stopping rule reads a CHANGE and a distance reads
#' a POINT: where the refit lands on its mode in one step the objective does
#' not move between iterations, the objective-stall guard fires, and the run
#' is labelled stopped rather than converged. The better the point, the
#' likelier the flag denies it, which is why a warning printed off the flag
#' would fire loudest where there is least to warn about.
#'
#' The held coordinate is excluded because it is not free: its score is the
#' quantity Rao's statistic reads and does not vanish under the restriction.
#'
#' @param score The restricted fit's score, all coordinates.
#' @param information Its penalized information, all coordinates.
#' @param at The position of the held coordinate.
#'
#' @return A single number, or `NA` where the free block could not be
#'   inverted there -- which is itself a reason to doubt the point.
#'
#' @seealso [mode_error_limit()], the limit it is read against, and
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
#' The positions of the coefficients whose penalty is not differentiable
#' everywhere -- a lasso, a SCAD, an MCP -- as an integer vector.
#'
#' @details
#' It is the question [statmod_blocks()] asks to decide which of the two
#' fitting schemes a block goes to, asked of the coordinates instead of the
#' blocks, and it reads the same predicate ([penalty_has_kink()]) so the two
#' cannot disagree. A penalty that is twice differentiable contributes
#' nothing, which is what makes a ridge, a random effect and a smooth
#' testable.
#'
#' A penalty's own enumeration is what answers, so a partial penalty -- one
#' covering some of a term's coordinates and not the rest -- names only what
#' it covers, and a term is not classified as a whole. A structural term's
#' penalty covers positions among the TERM'S own parameters rather than the
#' design's, and contributes nothing here.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#'
#' @return An integer vector of stacked positions, possibly empty.
#'
#' @seealso [statmod_penalized()], the enumeration it reads, and
#'   [statmod_blocks()], which asks the same question of a whole block.
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
#' The Wald, likelihood-ratio, Rao score or Terrell gradient statistic for
#' the hypothesis that one coefficient equals `value`, with its p-value.
#'
#' @details
#' All four read the PENALIZED objective, which is the function the fit
#' maximizes. Writing \eqn{\ell_p(\beta) = \ell(\beta) - \rho(\beta;\theta)}
#' for it, \eqn{\hat\beta} for the unrestricted estimates, \eqn{\tilde\beta}
#' for those with \eqn{\beta_j} held at \eqn{b}, \eqn{U = \partial\ell_p} for
#' its score and \eqn{K} for the penalized information, the four are
#' \deqn{W = (\hat\beta_j - b)^2 / \widehat{\mathrm{Var}}(\hat\beta_j),
#' \qquad LR = 2\{\ell_p(\hat\beta) - \ell_p(\tilde\beta)\},}
#' \deqn{S = U_j(\tilde\beta)^2 \, [K(\tilde\beta)^{-1}]_{jj}, \qquad
#' T = U_j(\tilde\beta)\,(\hat\beta_j - b),}
#' each compared with a \eqn{\chi^2_1}. Where the model carries no penalty
#' \eqn{\ell_p} IS the log-likelihood and these are the classical four.
#'
#' The score statistic's general form is \eqn{U'K^{-1}U}, which reduces to
#' the product above because every other component of \eqn{U} vanishes at the
#' restricted maximum -- and that is worth knowing rather than assuming,
#' since it is what says the restricted fit converged.
#'
#' # Why the penalized objective and not the likelihood
#'
#' Reading the likelihood alone would make the four incoherent with one
#' another and the second of them not a ratio at all. \eqn{\tilde\beta}
#' maximizes \eqn{\ell_p} under the restriction and not \eqn{\ell}, so
#' \eqn{\ell(\tilde\beta)} can exceed \eqn{\ell(\hat\beta)} and the
#' difference come out NEGATIVE -- measured on a smooth, holding
#' `s(z).z3` at two values either side of its estimate gives -0.0379 and
#' -0.1075 where \eqn{\ell_p} gives +0.0054 and +0.1355. The other three
#' were already reading the penalized objective, \eqn{U} and \eqn{K} being
#' its gradient and its curvature and the Wald variance being
#' \eqn{(H+S)^{-1}}, so one table would have carried two different tests.
#'
#' What it means for a SHRUNK coordinate is the reading
#' [vcov.StatmodFit()] already gives it under `type = "bayesian"`: the
#' penalty is a log-prior, so \eqn{\ell_p} is a log posterior at the
#' hyperparameters in force and the interval by inversion is the set the
#' posterior does not reject, conditional on them. Nothing new is claimed;
#' the claim the Wald interval on such a row already carries is extended to
#' the other three. Measured, the two agree as they must: on a smooth at
#' \eqn{n = 300} the inverted interval and the Wald one differ by 8.7e-04
#' on `s(z).z1`, 1.1e-03 on an unpenalized covariate and 2.5e-03 on the
#' smooth's linear column, which is second-order agreement, and what the
#' inversion adds is the asymmetry.
#'
#' # What separates them
#'
#' Wald reads the unrestricted fit alone and needs no refit, which is why it
#' is the one every summary prints and also the one that is NOT invariant to
#' how the parameter is written: a reparametrization moves it and moves none
#' of the other three. It is also the one that loses power against a distant
#' alternative, the Hauck-Donner effect, the curvature being read at
#' \eqn{\hat\beta} rather than under the null.
#'
#' The gradient statistic needs the restricted fit but NO matrix inversion at
#' all, which is its whole point: it is the score contracted with the
#' distance between the two fits.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()].
#' @param value The value under the null, a single number. `0` by default.
#' @param test One of `"wald"`, `"lr"`, `"score"`, `"gradient"`.
#' @param type Which variance the Wald statistic reads, as
#'   [vcov.StatmodFit()] takes it.
#' @param mode_error Whether to read how far above its mode the restricted
#'   fit stopped ([restricted_mode_error()]). `FALSE` by default, because it
#'   costs one Hessian where the likelihood ratio and the gradient statistic
#'   need none, and the two loops that call this in quantity --
#'   [summary.StatmodFit()], one row at a time, and [statmod_invert()], five
#'   to seven times per interval -- do not read it.
#'
#' @return A list with `test`, `statistic`, `df`, `p.value`, `converged` --
#'   `NA` for the Wald statistic and the restricted fit's flag for the other
#'   three -- and `mode_error`, `NA` unless it was asked for.
#'
#' @seealso [statmod_restrict()], which the three restricted statistics read.
#'
#' @keywords internal
statmod_stat_at <- function(fit, param, coefname, value = 0, test = "wald",
                         type = "bayesian", mode_error = FALSE) {
  test <- match.arg(test, c("wald", "lr", "score", "gradient"))
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

  if (identical(test, "wald")) {
    # the only one that needs no refit, and the only one a reparametrization
    # moves. Nothing was refitted, so there is no mode to be above
    V <- stats::vcov(fit, type = type)
    lab <- rownames(coef_labels(spec, design))[at]
    v <- if (is.na(at) || !lab %in% rownames(V)) NA_real_ else V[lab, lab]
    if (!isTRUE(is.finite(v)) || v <= 0) return(out(NA_real_, NA))
    return(out((bhat - value)^2 / v, NA))
  }

  r <- statmod_restrict(fit, param, coefname, value)
  me <- if (isTRUE(mode_error)) r$mode_error() else NA_real_
  if (identical(test, "lr")) {
    # 2 times the rise in the PENALIZED OBJECTIVE, which is the function both
    # fits maximize, so the difference is non-negative by construction. Read
    # on the log-likelihood alone it is not: the restricted point maximizes
    # the penalized objective, not the likelihood, so it can sit ABOVE the
    # unrestricted point there -- measured at -0.0379 and -0.1075 on a
    # smooth's own coordinate. With no penalty the two are the same number.
    return(out(2 * (r$objective - fit@objective), r$converged, me))
  }
  u <- r$score[[r$at]]
  if (identical(test, "gradient")) {
    return(out(u * (bhat - value), r$converged, me))
  }
  K <- r$information()
  V <- tryCatch(solve_pd(K, "the penalized information at the restricted fit",
                         r$labels),
                error = function(e) NULL)
  if (is.null(V)) return(out(NA_real_, r$converged, me))
  out(u^2 * V[r$at, r$at], r$converged, me)
}


#' A Test of One Coefficient Against One Value
#'
#' @description
#' The result of testing that one coefficient of a fit equals one value,
#' carrying the statistic, its distribution, the p-value and what was
#' tested.
#'
#' @details
#' The four statistics are compared with the same \eqn{\chi^2_1}, so `df` is
#' always 1 and the object holds one kind of quantity whatever produced it.
#'
#' # Whether the refit reached its mode
#'
#' `mode_error` and `converged` both describe the restricted refit and both
#' are `NA` for the Wald statistic, which reads the unrestricted fit and
#' refits nothing. They answer different questions and only the first is
#' printed.
#'
#' `converged` is the inner optimizer's flag: whether a stopping rule fired,
#' a boolean about a threshold on a score whose size depends on the model.
#' `mode_error` is [restricted_mode_error()], how much of the penalized
#' objective is still on the table at the point the refit returned, in
#' log-likelihood units against [mode_error_limit()] -- the rule statmodels7
#' 0.81.0 established for exactly this question.
#'
#' The flag is not merely a false negative: it is anti-correlated with the
#' quality of the point, because a stopping rule reads a change and a
#' distance reads a point. Measured on `y ~ x` with a Poisson response at
#' \eqn{n = 200} over nineteen held values, all nineteen are at their mode
#' and the four the flag rejects are the four best-located of them. See
#' [restricted_mode_error()].
#'
#' A statistic read off a refit that stopped short is the value at wherever
#' it stopped, and `print()` says so -- reading `mode_error`, so that it says
#' it where it is true.
#'
#' The class exists rather than a bare list because everything a caller of
#' this toolkit receives is an object with declared properties. R's own
#' tests return an `htest`, which is a list with a class attribute and no
#' contract; nothing here can validate one, and a reader cannot ask it what
#' it carries.
#'
#' @param test Which of the four statistics, as [statmod_test()] was asked
#'   for it.
#' @param statistic Its value, compared with a \eqn{\chi^2_1}.
#' @param df One, always.
#' @param p.value The upper tail of that \eqn{\chi^2_1}.
#' @param estimate The unrestricted estimate of the coefficient.
#' @param null_value The value it was tested against.
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
#' S7::S7_inherits(statmod_test(fit, "mu", "x", 1), StatmodTest)
#'
#' @name StatmodTest-class
#' @aliases StatmodTest
#' @keywords internal
#' @export
StatmodTest <- S7::new_class("StatmodTest",
  properties = list(
    # which of the four statistics, as `statmod_test()` was asked for it
    test = S7::class_character,
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
#' The statistic, its degrees of freedom and the p-value on one line, under
#' a title naming the test, in the layout R uses for a hypothesis test.
#'
#' @param x A [StatmodTest()].
#' @param digits Significant digits.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#'
#' @keywords internal
print.StatmodTest <- function(x, digits = 4L, ...) {
  nm <- c(wald = "Wald", lr = "Likelihood-ratio", score = "Rao score",
          gradient = "Terrell gradient")[[x@test]]
  cat("\n\t", nm, " test on a statmod fit\n\n", sep = "")
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
  if (refitted && !isTRUE(me <= mode_error_limit())) {
    cat("\nthe restricted fit stopped above its own mode",
        if (is.na(me)) "" else
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
#' `summary(fit, test =)` reports the same four statistics for every row of
#' the coefficient tables, and always against zero. This tests ONE
#' coefficient against ANY value, which is the question a summary cannot be
#' asked: whether a slope is one, whether an elasticity is unity, whether a
#' coefficient matches a value fixed outside the data.
#'
#' What the statistics are, what they read and what they cost is at
#' [statmod_stat_at()]; the restricted fit the last three are built on is
#' [statmod_restrict()], and the interval obtained by inverting one is
#' [confint.StatmodFit()].
#'
#' Reading one row costs one restricted refit, where `summary(test =)` costs
#' one per row: measured on a Poisson fit of twenty-five covariates, one
#' statistic is 0.010 s against the summary's 0.220 s over twenty-six rows.
#' It also reads the refit's mode error, which is one Hessian -- 63 per cent
#' of the refit on that fit and 2.3 per cent on a penalized smooth, and none
#' at all for Rao's score, which reads the same matrix anyway. The summary
#' and the interval do not pay it.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()]:
#'   the distribution parameter and the coefficient's name in its equation.
#' @param value The value under the null, a single number. `0` by default,
#'   which is what a summary reports.
#' @param test One of `"wald"`, `"lr"`, `"score"`, `"gradient"`.
#' @param type Which variance the Wald statistic reads, passed to
#'   [vcov.StatmodFit()].
#'
#' @return A [StatmodTest()].
#'
#' @seealso [summary.StatmodFit()] for every row against zero,
#'   [confint.StatmodFit()] for the interval a test inverts.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200))
#' dd$y <- rpois(200, exp(0.3 + dd$x))
#' fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
#' ## the summary asks whether the slope is zero; this asks whether it is one
#' statmod_test(fit, "mu", "x", 1, "lr")
#' statmod_test(fit, "mu", "x", 0, "lr")
#'
#' @export
statmod_test <- function(fit, param, coefname, value = 0, test = "wald",
                         type = "bayesian") {
  # the mode error is asked for HERE and nowhere else: this is the one
  # surface a reader reads one row of, so one extra Hessian is affordable,
  # where the two loops that call the same helper in quantity are not
  s <- statmod_stat_at(fit, param, coefname, value, test, type,
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
#' The set of values a test does not reject at `level`, found by locating the
#' two points where the statistic crosses its critical value.
#'
#' @details
#' An interval by inversion is \eqn{\{b : T(b) \le \chi^2_{1,\alpha}\}}, and
#' it is not symmetric about the estimate unless the statistic is a parabola
#' in \eqn{b} -- which the Wald statistic is BY CONSTRUCTION, the curvature
#' being read once at \eqn{\hat\beta}, and which the other three are not.
#' That asymmetry is the whole reason to invert one of the others: it follows
#' the likelihood's own shape rather than a quadratic fitted at its top.
#'
#' The search starts from the Wald interval, widens by a factor of 1.6 until
#' the statistic exceeds its critical value, and then bisects. Each
#' evaluation is one restricted refit, warm-started at the unrestricted
#' estimates and run at the fitted hyperparameters, so an interval costs a
#' few dozen inner fits and no outer search at all.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()].
#' @param level The confidence level, a single number in (0, 1).
#' @param test One of `"lr"`, `"score"`, `"gradient"`, `"wald"`.
#' @param type Which variance sets the starting bracket, and which the Wald
#'   statistic reads.
#' @param maxit How many widenings to try on each side before giving up.
#'
#' @return A numeric vector of two, `lower` and `upper`. An end the search
#'   could not bracket is `NA`, which is what an unbounded interval looks
#'   like and is reported rather than guessed at.
#'
#' @seealso [statmod_stat_at()], the statistic being inverted.
#'
#' @keywords internal
statmod_invert <- function(fit, param, coefname, level = 0.95, test = "lr",
                           type = "bayesian", maxit = 40L) {
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
  V <- stats::vcov(fit, type = type, readable = FALSE)
  lab <- rownames(coef_labels(spec, design))[at]
  se <- if (lab %in% rownames(V)) sqrt(V[lab, lab]) else NA_real_
  if (!isTRUE(is.finite(se)) || se <= 0) return(c(lower = NA_real_,
                                                  upper = NA_real_))

  gap <- function(b) {
    s <- statmod_stat_at(fit, param, coefname, b, test, type)$statistic
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
#' The stacked positions of the coefficients the three restricted statistics
#' are defined for: every coordinate less the ones under a kinked penalty,
#' less the aliased ones, and none at all where the model carries a
#' structural term.
#'
#' @details
#' It is the enumeration [statmod_restrict()] refuses one coordinate at a
#' time, asked once for the whole vector instead. A caller reporting a table
#' needs it that way round, since each restricted statistic is a refit and an
#' interval by inversion is a few dozen: a row that cannot be tested has to
#' be recognized before it is paid for rather than after.
#'
#' A coordinate under a penalty that is twice differentiable IS included --
#' a smooth's linear column and its rotated coordinates, a ridge, a random
#' effect -- and the restricted fit holds it exactly as it holds an
#' unpenalized one, [fit_smooth()] enforcing the hold in the solve. Only a
#' kinked penalty is excluded, for the reasons [statmod_restrict()] gives.
#'
#' @param fit A [StatmodFit()].
#' @param spec,design The specification and its design, where the caller has
#'   them already.
#'
#' @return An integer vector of stacked positions, possibly empty, named by
#'   the labels [coef_labels()] gives them.
#'
#' @seealso [statmod_restrict()], which refuses the same coordinates one at a
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
