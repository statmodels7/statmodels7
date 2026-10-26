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
