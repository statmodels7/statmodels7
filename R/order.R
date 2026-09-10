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
