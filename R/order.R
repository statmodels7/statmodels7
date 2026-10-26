# WHAT EACH CONSUMER NEEDS, AND SAYING SO WHERE THE QUESTION IS POSED.
#
# distributions7 records how many derivatives a log-density has in each of its
# parameters (`params_order()`).  Every consumer here needs a different number
# of them, and the number is a property of the mathematics rather than a
# convention:
#
#   1  a score continuous in the coefficients          U = X' dl/deta
#   2  IWLS, Newton, vcov(), and the determinant       K = -grad^2 l + S
#      of a marginal criterion
#   3  the exact gradient of a marginal criterion      it reads dK/dbeta
#   4  its exact Hessian
#   5  the fourth derivative of a filtered predictor
#
# A family with an absolute value of the residual in its log-density has order
# 0 in its location, so no row of that scale holds: the curvature every one of
# them reads sits on a set of Lebesgue measure zero in the response.  Before
# this file the fit died at the first evaluation of the criterion with a
# message about the inner fit not having converged, which is true and is three
# layers from the cause.

#' Whether a Family Carries Enough Derivatives
#'
#' @description
#' Returns `TRUE` when the log-density has at least `need` continuous
#' derivatives in every parameter for which the model carries an equation.
#' The consumers that can fall back on something else use the result as a
#' gate, and the consumers that cannot use it through [order_shortfall()] as
#' an error message.
#'
#' @param distrib A distributions7 family.
#' @param need The order required, an integer from 1 to 5: 1 for a score
#'   continuous in the coefficients, 2 for the curvature, 3 and 4 for the exact
#'   gradient and Hessian of a marginal criterion, and 5 for the fourth
#'   derivative of a filtered predictor.
#'
#' @return A single logical. `NA` orders count as unavailable: a wrapper of a
#'   kinked family records the kink without declaring its order, and an order
#'   that has not been established is not evidence of smoothness.
#'
#' @examples
#' statmodels7:::order_available(distributions7::gaussian1_distrib(), 4L)
#' statmodels7:::order_available(distributions7::laplace_distrib(), 2L)
#'
#' @seealso [order_shortfall()] for the message, and
#'   [distributions7::params_order()] for the quantity.
#' @keywords internal
order_available <- function(distrib, need) {
  ord <- tryCatch(distributions7::params_order(distrib), error = function(e) NULL)
  if (is.null(ord)) return(TRUE)          # a family that cannot answer is not accused
  all(!is.na(ord) & ord >= need)
}

#' Why a Consumer Cannot Run on This Family
#'
#' @description
#' Builds the error message of a consumer that cannot run on the family. The
#' message names the family, the parameter, the order that the consumer needs
#' against the order that the family has, and what to write instead. The
#' result is `NULL` when the family has enough derivatives.
#'
#' @details
#' The remedy differs between two shapes of shortfall, and the message
#' distinguishes them. A family that DECLARES a kink has a known order and the
#' message names the composition that puts it there. A wrapper of such a family
#' declares none, so the order is not established rather than known to be low,
#' and the message reports that the order is not established and gives no
#' number.
#'
#' @param distrib A distributions7 family.
#' @param need The order that the consumer needs.
#' @param what A noun phrase naming the consumer, used as the subject.
#'
#' @return A single string, or `NULL`.
#'
#' @examples
#' statmodels7:::order_shortfall(distributions7::laplace_distrib(), 2L,
#'                               "the REML criterion")
#'
#' @seealso [order_available()], [distributions7::kink_decomposition()].
#' @keywords internal
order_shortfall <- function(distrib, need, what) {
  ord <- tryCatch(distributions7::params_order(distrib), error = function(e) NULL)
  if (is.null(ord)) return(NULL)
  bad <- names(ord)[is.na(ord) | ord < need]
  if (!length(bad)) return(NULL)

  fam <- distrib@distrib_name
  kd  <- tryCatch(distributions7::kink_decomposition(distrib), error = function(e) NULL)
  p1  <- bad[1L]
  have <- ord[[p1]]

  # what is wrong, in the terms the family itself declares
  cause <- if (is.na(have)) {
    sprintf(paste0("'%s' records '%s' as non-smooth and declares no",
                   " kink_decomposition(), so the order there is not established"),
            fam, p1)
  } else {
    sprintf("'%s' has %s of them in '%s'", fam, format(have), p1)
  }
  how <- if (!is.null(kd)) {
    sprintf(" -- its log-density carries %s of an argument that moves with '%s', so the curvature that %s reads sits on a set of measure zero in the response",
            switch(kd@phi,
                   abs = "an absolute value", hinge = "a hinge",
                   hinge2 = "a squared hinge", hinge3 = "a cubed hinge",
                   step = "a step", "a non-smooth composition"),
            p1, what)
  } else ""

  paste0(what, " needs ", format(need),
         " derivatives of the log-density in every parameter, and ", cause, how,
         ". Fit with outer_criterion = NULL, which holds each smoothing",
         " hyperparameter at the value written in its term or at its starting",
         " value, or model '", p1,
         "' with a family that is smooth in it.")
}

#' Derivative Orders Read by a Criterion
#'
#' @description
#' Raises [order_shortfall()]'s message when the family does not carry the two
#' derivatives that every criterion reads. [statmod()] calls the function once
#' the criterion is known to have something to do, so that a model whose
#' criterion is inert does not signal an error that cannot arise.
#'
#' @param distrib A distributions7 family.
#' @param method An [OuterMethod()], or `NULL`.
#'
#' @return `NULL`, invisibly. The function is called for the error.
#'
#' @examples
#' statmodels7:::assert_criterion_order(distributions7::gaussian1_distrib(), reml())
#'
#' @seealso [order_shortfall()]
#' @keywords internal
assert_criterion_order <- function(distrib, method) {
  if (is.null(method)) return(invisible(NULL))
  why <- order_shortfall(distrib, 2L,
                         sprintf("the %s criterion", toupper(method@kind)))
  if (!is.null(why)) stop(why, call. = FALSE)
  invisible(NULL)
}

#' Expected-Information Criteria and the Scores' Outer Product
#'
#' @description
#' Raises an error naming the family and the remedy when a criterion defined
#' on the expected information would read, for this family, the outer product
#' of the scores at the data instead of an expectation.
#'
#' @details
#' A criterion on the expected information replaces the observed curvature of
#' each observation, \eqn{-\partial^2\ell_i/\partial\eta\,\partial\eta^\top}
#' at \eqn{y_i}, with the Fisher information
#' \eqn{\mathcal{I}(\theta_i) = -\mathbb{E}[\partial^2\ell/\partial\eta\,
#' \partial\eta^\top]}, a function of the parameters alone. Where a family
#' writes that information out, [distributions7::expected_hessian_exact()]
#' returns `TRUE` and the check is skipped. Where a family does not write it
#' out, for example a family written outside the toolkit, the information
#' falls back on the strategy named by the `approx` of the criterion's
#' [expected_info()], and the default, [distributions7::expect_opg()],
#' returns \eqn{s_i s_i^\top} with \eqn{s_i} the score at \eqn{y_i}. That is an
#' unbiased estimate of \eqn{\mathcal{I}(\theta_i)} from one observation, not
#' the information itself: it depends on the response, as the observed curvature
#' does, and it is neither of them. The criterion built on it is neither the
#' Laplace approximation nor its Fisher variant, and it has no exact outer
#' gradient, so its search is derivative-free and [statmod_certificate()]
#' returns `unknown`.
#'
#' A criterion on the observed information has an exact outer gradient and
#' Hessian and is searched with Newton's method; the outer-product route has
#' neither and is searched without derivatives.
#'
#' The expectation can be requested with another `approx` in the criterion's
#' [expected_info()], and this check does not reject it:
#' [distributions7::expect_bartlett()] sums or integrates over the support and
#' gives a genuine expected information, at a much higher cost per evaluation.
#'
#' @param distrib A distributions7 family.
#' @param method An [OuterMethod()], or `NULL`.
#'
#' @return `NULL`, invisibly. The function is called for the error.
#'
#' @examples
#' statmodels7:::assert_criterion_information(
#'   distributions7::gaussian1_distrib(), reml(information = expected_info()))
#' # the Poisson-inverse Gaussian computes its expected information exactly
#' statmodels7:::assert_criterion_information(
#'   distributions7::pig1_distrib(), reml(information = expected_info()))
#'
#' @seealso [expected_is_opg()], the predicate, and [assert_criterion_order()].
#' @keywords internal
assert_criterion_information <- function(distrib, method) {
  if (is.null(method) || !info_expected(method@hessian)) {
    return(invisible(NULL))
  }
  if (!expected_is_opg(distrib, method@hessian@approx)) return(invisible(NULL))
  crit <- sprintf("%s(information = expected_info())", method@kind)
  stop(sprintf(paste0(
    "%s reads the expected information, which '%s' does not write out: with ",
    "approx = expect_opg() it would read the outer product of the scores at ",
    "the data, which depends on the response and is not an expectation, so ",
    "the criterion would be neither the Laplace approximation nor its Fisher ",
    "variant and would have no exact outer derivatives. Use ",
    "%s(information = observed_info()), whose outer gradient and Hessian are ",
    "exact on this family, or name the expectation with ",
    "%s(information = expected_info(approx = expect_bartlett())), which sums ",
    "or integrates over the support at a much higher cost per evaluation."),
    crit, distrib@distrib_name, method@kind, method@kind), call. = FALSE)
}

#' Does the Family's Expected Information Depend on the Response?
#'
#' @description
#' The result is `TRUE` where [distributions7::distrib_expected_hessian()],
#' called with `approx`, returns a matrix that changes with the response at
#' fixed parameters, which an expectation cannot do.
#'
#' @details
#' The answer comes from a probe and not from a list of families, so that a
#' family written later is covered. The probe draws two distinct responses at
#' one parameter value and compares the two answers. A family that
#' [distributions7::expected_hessian_exact()] reports as writing its
#' information out is not probed, and the probe is not run for any `approx`
#' other than [distributions7::expect_opg()], the only strategy that reads the
#' response. For a
#' truncated family, the predicate returns `FALSE` although the expected
#' information is a quadrature that does not depend on the response, so this
#' function also returns `FALSE`, and the answer is determined by the probe and
#' not by the predicate. The caller's random stream is restored.
#'
#' @param distrib A distributions7 family.
#' @param approx The approximation, an object that inherits from
#'   [distributions7::ExpectApprox].
#'
#' @return A single logical.
#'
#' @seealso [assert_criterion_information()]
#' @keywords internal
expected_is_opg <- function(distrib, approx) {
  if (!S7::S7_inherits(approx, distributions7::ExpectOPG)) return(FALSE)
  if (S7::S7_inherits(distrib, distributions7::multivariate_distrib)) return(FALSE)
  exact <- tryCatch(distributions7::expected_hessian_exact(distrib),
                    error = function(e) TRUE)
  if (isTRUE(exact)) return(FALSE)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(20260924L)
  th <- tryCatch(distributions7::generate_random_theta(distrib),
                 error = function(e) NULL)
  if (is.null(th)) return(FALSE)
  y <- tryCatch(unique(distributions7::distrib_rng(distrib, 20L, th)),
                error = function(e) NULL)
  y <- y[is.finite(y)]
  if (length(y) < 2L) return(FALSE)
  at <- function(v) {
    tryCatch(unlist(distributions7::distrib_expected_hessian(distrib, v, th,
                                                             approx = approx)),
             error = function(e) NULL)
  }
  e1 <- at(y[[1L]])
  e2 <- at(y[[2L]])
  if (is.null(e1) || is.null(e2) || length(e1) != length(e2)) return(FALSE)
  if (!all(is.finite(e1)) || !all(is.finite(e2))) return(FALSE)
  !isTRUE(all.equal(e1, e2, tolerance = 1e-8, check.attributes = FALSE))
}

#' Prediction-Error Criteria and a Structural Term's Own Parameters
#'
#' @description
#' Raises an error naming the penalty and the remedy when `aic()`, `bic()` or
#' `cv()` would have to select a hyperparameter of a penalty over the own
#' parameters of a structural term: the group departures of a score-driven
#' filter over a panel, or a covariance class that reaches into one.
#'
#' @details
#' A prediction-error criterion reads the effective degrees of freedom, and
#' their derivative, over the coefficients alone, where such a penalty covers
#' nothing. The exact gradient of `aic()` there is zero, so a search would
#' stop at its first evaluation and report convergence. `cv()` and the path
#' along which a kinked penalty is swept are built on the same blocks.
#' [reml()] and [ml()] span the term's parameters in their determinant and
#' estimate a twice differentiable penalty there. A penalty with a kink there
#' has no criterion that can select its hyperparameter. The hyperparameter can
#' be held by writing its value in the term, for example
#' `lasso(~id, lambda = 2)`, and the penalty is then fitted over the term's own
#' parameters by a proximal route, without a criterion.
#'
#' Which penalties a criterion reaches depends on its role, and only those are
#' checked: an `outer_criterion` reaches the smooth penalties, and the
#' kinked ones as well when no `sparse_criterion` is given, while a
#' `sparse_criterion` reaches only the kinked ones.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param method An [OuterMethod()], or `NULL`.
#' @param reach Which penalties the criterion selects: `"all"`, `"smooth"` or
#'   `"kinked"`.
#'
#' @return `NULL`, invisibly. The function is called for the error.
#'
#' @examples
#' dd <- data.frame(x = runif(60))
#' dd$y <- rnorm(60, dd$x)
#' spec <- statmodels7:::statmod_spec(y ~ x, distributions7::gaussian1_distrib(),
#'                                    dd)
#' statmodels7:::assert_criterion_reach(spec, statmodels7:::statmod_design(spec),
#'                                      aic())
#'
#' @seealso [assert_criterion_order()], [structural_penalized()]
#' @keywords internal
assert_criterion_reach <- function(spec, design, method,
                                   reach = c("all", "smooth", "kinked")) {
  reach <- match.arg(reach)
  if (is.null(method) || !method@kind %in% c("aic", "bic", "cv")) {
    return(invisible(NULL))
  }
  held_ids <- names(statmod_held_ids(statmod_penalty_keys(spec)))
  for (u in statmod_penalized(spec, design)) {
    if (!isTRUE(u$structural) && !isTRUE(u$mixed)) next
    kinked <- isTRUE(tryCatch(penalty_has_kink(u$penalty, u$key),
                              error = function(e) FALSE))
    if ((identical(reach, "smooth") && kinked) ||
        (identical(reach, "kinked") && !kinked)) next
    lab_of <- function(h) if (h %in% names(u$ids)) u$ids[[h]] else NA_character_
    free <- Filter(function(h) !(h %in% names(u$fixed)) && !(lab_of(h) %in% held_ids),
                   u$penalty@params)
    if (!length(free)) next
    what <- if (identical(method@kind, "cv")) {
      "cross-validation sweeps a path built over the coefficients alone"
    } else {
      paste0("a prediction-error criterion reads the degrees of freedom and ",
             "their gradient over the coefficients alone")
    }
    # ⚠️ THE MESSAGE NAMES NO REMEDY FOR A PENALTY WITH A KINK. Holding its
    # hyperparameter at a value written on the term does fit: sparse_fit()
    # sends the block to sparse_fit_structural(), a proximal route over the
    # term's own parameters, which is not a selection by a criterion.
    remedy <- if (kinked) {
      paste0("A penalty with a kink there has no criterion that can select ",
             "it: reml() and ml() read a curvature that does not exist at ",
             "the kink.")
    } else {
      "Use reml() or ml(), whose determinant spans the term's parameters as well."
    }
    stop(sprintf(paste0(
      "%s() cannot select the hyperparameter '%s' of '%s' in '%s': that ",
      "penalty covers the own parameters of a structural term, and %s, where ",
      "it covers nothing. %s"),
      method@kind, free[[1L]], short_keys(u$key), u$param, what, remedy),
      call. = FALSE)
  }
  invisible(NULL)
}
