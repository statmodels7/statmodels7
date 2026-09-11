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
#' `TRUE` when the log-density has at least `need` continuous derivatives in
#' every parameter the model carries an equation for. Used as a gate by the
#' consumers that can fall back on something else, and through
#' [order_shortfall()] as a refusal by the ones that cannot.
#'
#' @param distrib A distributions7 family.
#' @param need The order asked for; see the table in `R/order.R`.
#'
#' @return A single logical. `NA` orders count as unavailable: a wrapper of a
#'   kinked family records the kink without declaring its order, and an
#'   unestablished order is not evidence of smoothness.
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
#' The message a consumer refuses with: the family, the parameter, the order it
#' needs against the order there is, and what to write instead. `NULL` when
#' there is nothing to refuse.
#'
#' @details
#' Two shapes of shortfall are told apart, because the remedies differ. A
#' family that DECLARES a kink has a known order and the message says which
#' composition puts it there. A wrapper of such a family declares none, so the
#' order is unestablished rather than known to be low, and the message says
#' that instead of inventing a number.
#'
#' @param distrib A distributions7 family.
#' @param need The order the consumer needs.
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
    sprintf(" -- its log-density carries %s of an argument that moves with '%s', so the curvature %s reads sits on a set of measure zero in the response",
            switch(kd@phi,
                   abs = "an absolute value", hinge = "a hinge",
                   hinge2 = "a squared hinge", hinge3 = "a cubed hinge",
                   step = "a step", "a non-smooth composition"),
            p1, what)
  } else ""

  paste0(what, " needs ", format(need),
         " derivatives of the log-density in every parameter, and ", cause, how,
         ". Fit with outer_criterion = NULL, setting the hyperparameters",
         " yourself through hyper, or model '", p1,
         "' with a family that is smooth in it.")
}

#' Refuse Where a Criterion Cannot Be Computed
#'
#' @description
#' Raises [order_shortfall()]'s message when the family does not carry the two
#' derivatives every criterion reads. Called from [statmod()] once the
#' criterion is known to have something to do, so a model whose criterion is
#' inert is not refused for a reason that would never have arisen.
#'
#' @param distrib A distributions7 family.
#' @param method An [OuterMethod()], or `NULL`.
#'
#' @return `NULL`, invisibly. Called for the error.
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

#' Refuse a Prediction-Error Criterion on a Structural Term's Own Parameters
#'
#' @description
#' Raises an error naming the penalty and the remedy when `aic()`, `bic()` or
#' `cv()` would have to select a hyperparameter of a penalty over the own
#' parameters of a structural term: the deviations of a score-driven filter
#' over a panel, or a covariance class that reaches into one.
#'
#' @details
#' A prediction-error criterion reads the effective degrees of freedom, and
#' their derivative, over the coefficients alone, where such a penalty covers
#' nothing. Measured before this refusal existed, the exact gradient of
#' `aic()` there was exactly 0, the search stopped at its first evaluation
#' reporting convergence, and the traces it priced were 2 and 1 where the
#' traces over the joint vector are 7.433 and 22.410. `cv()` and the path a
#' kinked penalty is swept along are built on the same blocks, and for a lasso
#' over a filter's parameters the path died on "'from' must be a finite
#' number". [reml()] and [ml()] span the term's parameters in their
#' determinant and estimate a twice differentiable penalty there. A penalty
#' with a kink there has no criterion at all, and holding its hyperparameter
#' does not fit either: measured, a lasso over a filter's own parameters held
#' at `lambda = 2` dies inside the coordinate-descent route, in the published
#' release as well, so the message names no remedy for it.
#'
#' Which penalties a criterion reaches depends on its role, and only those are
#' asked about: an `outer_criterion` reaches the smooth penalties, and the
#' kinked ones as well when no `sparse_criterion` is given, while a
#' `sparse_criterion` reaches the kinked ones alone.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param method An [OuterMethod()], or `NULL`.
#' @param reach Which penalties the criterion selects: `"all"`, `"smooth"` or
#'   `"kinked"`.
#'
#' @return `NULL`, invisibly. Called for the error.
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
    # ⚠️ NO REMEDY IS NAMED FOR A PENALTY WITH A KINK, and a measurement is
    # why: holding its hyperparameter at a value written on the term does not
    # fit either. A lasso over a filter's own parameters held at lambda = 2
    # dies on "subscript out of bounds" in coord_block(), in the published
    # 0.124.0 as well, the coordinate-descent route reading the penalized
    # coordinates as columns of a design, which those parameters are not.
    remedy <- if (kinked) {
      paste0("A penalty with a kink there has no criterion that can select ",
             "it: reml() and ml() read a curvature the kink does not have.")
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
