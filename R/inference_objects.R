#' @include predict.R
NULL

# The two choices that every inference function of a fit shares, written as
# objects: the statistic with which a coefficient is tested, and the
# convention that sets the variance matrix. A character string used to select
# each one. Passing a string signals an error that names the object to use.


#' @title Statistics for Testing a Coefficient
#' @name test_statistic
#'
#' @description
#' Four statistics test the hypothesis that one coefficient of a fit equals a
#' value: Wald's, the likelihood ratio, Rao's score and Terrell's gradient.
#' Each is an object built by [wald_test()], [lr_test()], [score_test()] or
#' [gradient_test()], and the object is passed as the `test` argument of
#' [statmod_test()], [summary.StatmodFit()] and [confint.StatmodFit()]. The
#' default is [wald_test()].
#'
#' `test_statistic` is the abstract parent of the four classes. It cannot be
#' instantiated, and serves to test membership with `S7::S7_inherits()`.
#'
#' @details
#' Write \eqn{\ell_p} for the penalized log-likelihood (the log-likelihood
#' where the model carries no penalty), \eqn{\hat\beta} for the estimates of
#' the fit, \eqn{\tilde\beta} for the estimates in which the coefficient
#' \eqn{\beta_j} is held at the tested value \eqn{b}, \eqn{U = \partial\ell_p}
#' for the penalized score and \eqn{K} for the penalized information. The four
#' statistics are
#' \deqn{W = \frac{(\hat\beta_j - b)^2}{\widehat{\mathrm{Var}}(\hat\beta_j)},
#'   \qquad LR = 2\{\ell_p(\hat\beta) - \ell_p(\tilde\beta)\},}
#' \deqn{S = U_j(\tilde\beta)^2\,[K(\tilde\beta)^{-1}]_{jj}, \qquad
#'   T = U_j(\tilde\beta)\,(\hat\beta_j - b),}
#' and each is compared with a \eqn{\chi^2_1} distribution. In
#' [summary.StatmodFit()] the value of a table is the signed root of the
#' statistic, with the sign of the estimate, and the p-value is that of the
#' \eqn{\chi^2_1} value.
#'
#' | Constructor | Statistic | What it needs | Cost |
#' |:--|:--|:--|:--|
#' | [wald_test()] | \eqn{W} | the unrestricted fit and its variance matrix | no refit |
#' | [lr_test()] | \eqn{LR} | the restricted refit and the objective at the unrestricted joint mode | one refit per value \eqn{b} |
#' | [score_test()] | \eqn{S} | the restricted refit and its penalized information | one refit and one Hessian per value \eqn{b} |
#' | [gradient_test()] | \eqn{T} | the restricted refit and its score | one refit per value \eqn{b}, no matrix inversion |
#'
#' The three statistics that need the restricted refit are the restricted
#' tests, and [test_restricted()] returns `TRUE` for them. The refit maximizes
#' \eqn{\ell_p} over every coefficient but \eqn{\beta_j}, at the
#' hyperparameters of the fit, with the inner method that the fit used, and it
#' starts from the unrestricted estimates; the outer search over the
#' hyperparameters is not repeated. The statistics are therefore conditional
#' on the smoothing that the data chose once. Their cost in each consumer is
#' as follows.
#'
#' * [statmod_test()] runs one refit.
#' * [summary.StatmodFit()] runs one refit per coefficient row, against zero.
#' * [confint.StatmodFit()] inverts the statistic, which takes a few dozen
#'   refits per interval, and it requires `readable = FALSE` for a restricted
#'   test, because a readable quantity is a function of several coefficients
#'   at once and one coordinate cannot be held in its place.
#'
#' The restricted refit is defined for a coefficient that is not under a kinked
#' penalty (a lasso, a SCAD or an MCP), that is not aliased and that belongs to
#' a model without a structural term. For any other coefficient the restricted
#' statistics are `NA`.
#'
#' The Wald statistic depends on how the coefficient is written: a nonlinear
#' reparametrization changes it, whereas the likelihood ratio, a difference of
#' objective values at the same two points, is not changed. It also loses power
#' against a distant alternative (the Hauck-Donner effect), because the
#' curvature is read at \eqn{\hat\beta} and not under the null. The interval
#' that inverts the Wald statistic is symmetric about the estimate. The
#' intervals that invert the other three follow the shape of the likelihood
#' and are asymmetric where it is. The variance that the Wald statistic reads
#' is set by the argument `type`, described at [vcov_type].
#'
#' @param name A single string, the short name of the statistic (`"wald"`,
#'   `"lr"`, `"score"` or `"gradient"`), which also gives the name of the
#'   constructor (`wald_test()` and so on).
#' @param label A single string, the name under which a result prints the
#'   statistic (for example `"Likelihood-ratio"`).
#'
#' @return Nothing: the class is abstract and its constructor signals an
#'   error. Every object that inherits from it has the properties `name` and
#'   `label`.
#'
#' @seealso [wald_test()], [lr_test()], [score_test()] and [gradient_test()]
#'   for the four statistics, [test_restricted()] for the distinction between
#'   them, [statmod_test()], [summary.StatmodFit()] and
#'   [confint.StatmodFit()] for the functions that take one, and [vcov_type]
#'   for the variance that the Wald statistic reads.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # the four statistics of one hypothesis, that the slope is 1.5
#' tests <- list(wald_test(), lr_test(), score_test(), gradient_test())
#' vapply(tests, function(tt) {
#'   statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
#'                test = tt)@statistic
#' }, numeric(1))
#'
#' # which of them refit the model under the restriction
#' vapply(tests, test_restricted, logical(1))
#'
#' # the class itself is abstract
#' try(test_statistic(name = "mine", label = "Mine"))
#' @export
test_statistic <- S7::new_class(
  "test_statistic",
  properties = list(name = S7::class_character, label = S7::class_character),
  abstract = TRUE
)

# The concrete classes are not exported: the constructors are the interface.
WaldTest <- S7::new_class("WaldTest", parent = test_statistic)

RestrictedTest <- S7::new_class("RestrictedTest", parent = test_statistic,
                                abstract = TRUE)

LrTest <- S7::new_class("LrTest", parent = RestrictedTest)

ScoreTest <- S7::new_class("ScoreTest", parent = RestrictedTest)

GradientTest <- S7::new_class("GradientTest", parent = RestrictedTest)


#' @title Wald's Test of a Coefficient
#'
#' @description
#' Selects Wald's statistic for the hypothesis that one coefficient equals a
#' value \eqn{b}:
#' \deqn{W = \frac{(\hat\beta_j - b)^2}{\widehat{\mathrm{Var}}(\hat\beta_j)},}
#' compared with a \eqn{\chi^2_1} distribution. It reads the unrestricted fit
#' and the variance matrix of [vcov.StatmodFit()] and refits nothing, which is
#' why it is the default statistic of every function that takes a `test`.
#'
#' @details
#' Choose it when the cost of a refit per coefficient matters, or when the
#' interval symmetric about the estimate is what is wanted. The likelihood
#' ratio, the score and the gradient statistics follow the shape of the
#' likelihood, and the Wald statistic does not; see [test_statistic] for the
#' comparison. The variance on which it depends is selected by the argument
#' `type` of the consuming function (see [vcov_type]).
#'
#' The constructor has no arguments.
#'
#' @return A `WaldTest` object, inheriting from [test_statistic], with the
#'   properties `name` (`"wald"`) and `label` (`"Wald"`).
#'
#' @seealso [lr_test()], [score_test()] and [gradient_test()] for the
#'   siblings, [test_statistic] for the comparison of the four,
#'   [statmod_test()], [summary.StatmodFit()] and [confint.StatmodFit()] for
#'   the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' wald_test()
#' statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
#'              test = wald_test())
#' confint(fit, test = wald_test())
#' @export
wald_test <- function() WaldTest(name = "wald", label = "Wald")

#' @title Likelihood-Ratio Test of a Coefficient
#'
#' @description
#' Selects the likelihood-ratio statistic for the hypothesis that one
#' coefficient equals a value \eqn{b}:
#' \deqn{LR = 2\{\ell_p(\hat\beta) - \ell_p(\tilde\beta)\},}
#' where \eqn{\ell_p} is the penalized log-likelihood, \eqn{\hat\beta} the
#' estimates of the fit and \eqn{\tilde\beta} the estimates in which the
#' coefficient is held at \eqn{b}. It is compared with a \eqn{\chi^2_1}
#' distribution.
#'
#' @details
#' The statistic is a difference of the objective at two points, so a
#' reparametrization of the coefficient does not change it, and the interval
#' that inverts it follows the shape of the likelihood. It requires one
#' restricted refit per value \eqn{b} (a few dozen per interval), at the
#' hyperparameters of the fit; where the outer criterion estimated some of the
#' coefficients, the unrestricted side is the joint mode, obtained by one
#' refit of the full model from the fitted coefficients. See [test_statistic]
#' for the comparison with the other three.
#'
#' The constructor has no arguments.
#'
#' @return An `LrTest` object, inheriting from [test_statistic], with the
#'   properties `name` (`"lr"`) and `label` (`"Likelihood-ratio"`).
#'
#' @seealso [wald_test()], [score_test()] and [gradient_test()] for the
#'   siblings, [test_statistic] for the comparison of the four,
#'   [statmod_test()], [summary.StatmodFit()] and [confint.StatmodFit()] for
#'   the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' lr_test()
#' statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
#'              test = lr_test())
#' confint(fit, parm = "mu:x", test = lr_test(), readable = FALSE)
#' @export
lr_test <- function() LrTest(name = "lr", label = "Likelihood-ratio")

#' @title Rao's Score Test of a Coefficient
#'
#' @description
#' Selects Rao's score statistic for the hypothesis that one coefficient
#' equals a value \eqn{b}:
#' \deqn{S = U_j(\tilde\beta)^2\,[K(\tilde\beta)^{-1}]_{jj},}
#' where \eqn{\tilde\beta} are the estimates in which the coefficient is held
#' at \eqn{b}, \eqn{U_j} is the penalized score of the coefficient and \eqn{K}
#' the penalized information. It is compared with a \eqn{\chi^2_1}
#' distribution.
#'
#' @details
#' The statistic is evaluated at the restricted estimates, so the unrestricted
#' fit is not read. It requires one restricted refit per value \eqn{b}, at the
#' hyperparameters of the fit, and the inverse of the penalized information at
#' the restricted point, so each evaluation also costs one Hessian. The
#' general form \eqn{U'K^{-1}U} reduces to the expression above because every
#' other component of the score vanishes at the restricted maximum. See
#' [test_statistic] for the comparison with the other three.
#'
#' The constructor has no arguments.
#'
#' @return A `ScoreTest` object, inheriting from [test_statistic], with the
#'   properties `name` (`"score"`) and `label` (`"Rao score"`).
#'
#' @seealso [wald_test()], [lr_test()] and [gradient_test()] for the
#'   siblings, [test_statistic] for the comparison of the four,
#'   [statmod_test()], [summary.StatmodFit()] and [confint.StatmodFit()] for
#'   the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' score_test()
#' statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
#'              test = score_test())
#' confint(fit, parm = "mu:x", test = score_test(), readable = FALSE)
#' @export
score_test <- function() ScoreTest(name = "score", label = "Rao score")

#' @title Terrell's Gradient Test of a Coefficient
#'
#' @description
#' Selects Terrell's gradient statistic for the hypothesis that one
#' coefficient equals a value \eqn{b}:
#' \deqn{T = U_j(\tilde\beta)\,(\hat\beta_j - b),}
#' where \eqn{\hat\beta} are the estimates of the fit, \eqn{\tilde\beta} the
#' estimates in which the coefficient is held at \eqn{b} and \eqn{U_j} the
#' penalized score of the coefficient. It is compared with a \eqn{\chi^2_1}
#' distribution.
#'
#' @details
#' The statistic is the product of the score and the distance between the two
#' fits. It requires one restricted refit per value \eqn{b}, at the
#' hyperparameters of the fit. Unlike the score statistic it reads no Hessian
#' and inverts no matrix, so a refit is its whole cost. See [test_statistic]
#' for the comparison with the other three.
#'
#' The constructor has no arguments.
#'
#' @return A `GradientTest` object, inheriting from [test_statistic], with the
#'   properties `name` (`"gradient"`) and `label` (`"Terrell gradient"`).
#'
#' @seealso [wald_test()], [lr_test()] and [score_test()] for the siblings,
#'   [test_statistic] for the comparison of the four, [statmod_test()],
#'   [summary.StatmodFit()] and [confint.StatmodFit()] for the functions that
#'   take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' gradient_test()
#' statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
#'              test = gradient_test())
#' confint(fit, parm = "mu:x", test = gradient_test(), readable = FALSE)
#' @export
gradient_test <- function() {
  GradientTest(name = "gradient", label = "Terrell gradient")
}


#' @title Whether a Test Refits the Model Under the Restriction
#'
#' @description
#' Returns `TRUE` for a statistic that is evaluated at the estimates in which
#' the tested coefficient is held at a value ([lr_test()], [score_test()] and
#' [gradient_test()]), and `FALSE` for the statistic that reads the
#' unrestricted fit alone ([wald_test()]).
#'
#' @details
#' A restricted statistic is about one coefficient, since one coordinate is
#' held. [confint.StatmodFit()] reads this generic to require
#' `readable = FALSE` for such a statistic, and the summary reads it to decide
#' whether a row costs a refit.
#'
#' @param test A [test_statistic] object.
#' @param ... Unused.
#'
#' @return A single logical.
#'
#' @seealso [test_statistic] for the four statistics and what each needs.
#'
#' @examples
#' test_restricted(wald_test())
#' test_restricted(lr_test())
#' @export
test_restricted <- S7::new_generic("test_restricted", "test",
  function(test, ...) S7::S7_dispatch())

S7::method(test_restricted, WaldTest) <- function(test, ...) FALSE
S7::method(test_restricted, RestrictedTest) <- function(test, ...) TRUE


#' @title Print a Test Statistic
#' @name print.test_statistic
#' @description
#' Prints one line, the label of the statistic.
#' @param x A [test_statistic] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(lr_test())
#' @keywords internal
S7::method(print, test_statistic) <- function(x, ...) {
  cat("<test> ", x@label, "\n", sep = "")
  invisible(x)
}


#' Check the Argument That Selects a Test Statistic
#'
#' @description
#' Returns `test` if it is a [test_statistic] object. A character string
#' signals an error that names the constructor to use, and any other value
#' signals an error that lists the four constructors.
#'
#' @param test The value of a `test` argument.
#' @param arg The argument's name, for the message.
#'
#' @return `test`.
#'
#' @keywords internal
as_test_statistic <- function(test, arg = "test") {
  if (S7::S7_inherits(test, test_statistic)) return(test)
  ctor <- c(wald = "wald_test()", lr = "lr_test()", score = "score_test()",
            gradient = "gradient_test()")
  if (is.character(test) && length(test) >= 1L && !is.na(test[[1L]]) &&
      test[[1L]] %in% names(ctor)) {
    stop(sprintf(paste0("'%s' must be a test_statistic object such as ",
                        "wald_test(); the string \"%s\" is no longer ",
                        "accepted. Use %s."),
                 arg, test[[1L]], ctor[[test[[1L]]]]), call. = FALSE)
  }
  stop(sprintf(paste0("'%s' must be a test_statistic object: wald_test(), ",
                      "lr_test(), score_test() or gradient_test()."), arg),
       call. = FALSE)
}


#' @title Variance Conventions of a Fit
#' @name vcov_type
#'
#' @description
#' Three conventions set the variance matrix of the coefficients of a fit:
#' the Bayesian one, the frequentist one and the unconditional one. Each is
#' an object built by [bayesian_vcov()], [frequentist_vcov()] or
#' [unconditional_vcov()], and the object is passed as the argument `type` of
#' [vcov.StatmodFit()], [confint.StatmodFit()], [summary.StatmodFit()],
#' [statmod_test()] and [predict.StatmodFit()]. The default is
#' [bayesian_vcov()].
#'
#' `vcov_type` is the abstract parent of the three classes. It cannot be
#' instantiated, and serves to test membership with `S7::S7_inherits()`.
#'
#' @details
#' Write \eqn{H} for the information of the log-likelihood and \eqn{S} for the
#' second derivative of the penalty, both at the estimates. The three matrices
#' are
#' \deqn{V_b = (H + S)^{-1}, \qquad V_f = (H+S)^{-1} H (H+S)^{-1}, \qquad
#'   V' = V_b + J V_\theta J'.}
#'
#' * [bayesian_vcov()] returns \eqn{V_b}, the posterior variance under the
#'   prior whose negative logarithm is the penalty, conditional on the
#'   hyperparameters. It carries the smoothing bias as variance, so an
#'   interval around a penalized term built from it covers at about its
#'   nominal rate.
#' * [frequentist_vcov()] returns \eqn{V_f}, the sandwich: the sampling
#'   variance of the penalized estimator at a fixed penalty. It is smaller
#'   than \eqn{V_b} and covers less.
#' * [unconditional_vcov()] returns \eqn{V'}. The first two matrices read the
#'   hyperparameters (smoothing parameters, variances of random effects,
#'   scales of a prior) at the value where the outer search stopped, as though
#'   that value were known. \eqn{V'} adds the variance that their estimation
#'   contributes, with \eqn{V_\theta} the variance of the estimated
#'   hyperparameters on their own free scale and
#'   \eqn{J = -(H + S)^{-1}\,\partial^2\rho/\partial\beta\,\partial\theta} the
#'   movement of the penalized mode under them (Wood, Pya and Safken, 2016).
#'   The correction is positive semi-definite, so \eqn{V'} is never narrower
#'   than \eqn{V_b}. It requires the curvature of the outer criterion, which
#'   the other two do not read, and it is the most expensive of the three.
#'
#' With no penalty \eqn{S = 0} and \eqn{V_b} and \eqn{V_f} are both
#' \eqn{H^{-1}} for a coefficient that the outer criterion does not estimate.
#' Where no hyperparameter was estimated by a differentiable criterion the
#' correction of the unconditional matrix is zero and it equals \eqn{V_b}.
#' The coefficients that [reml()] and [ml()] estimate on the criterion (the
#' unpenalized coefficients of every parameter but the first, by default) take
#' the curvature of that criterion under the Bayesian and the unconditional
#' conventions, which therefore differ from the frequentist one for them even
#' without a penalty. The complete description of the matrices, of the
#' coordinates that are held and of the warnings is at [vcov.StatmodFit()].
#'
#' @param name A single string, the short name of the convention
#'   (`"bayesian"`, `"frequentist"` or `"unconditional"`), which is also the
#'   stem of the constructor (`bayesian_vcov()` and so on). A result prints
#'   it.
#'
#' @return Nothing: the class is abstract and its constructor signals an
#'   error. Every object that inherits from it has the property `name`.
#'
#' @references
#' Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and model
#' selection for general smooth models. *Journal of the American
#' Statistical Association*, 111(516), 1548--1563.
#'
#' @seealso [bayesian_vcov()], [frequentist_vcov()] and
#'   [unconditional_vcov()] for the three conventions, [vcov.StatmodFit()]
#'   for the matrices, [test_statistic] for the tests that read a variance.
#'
#' @examples
#' set.seed(1)
#' ds <- data.frame(x = runif(200))
#' ds$y <- sin(2 * pi * ds$x) + rnorm(200, sd = 0.3)
#' fs <- statmod(y ~ s(x, bspline_smooth(k = 10)),
#'               distributions7::gaussian1_distrib(), ds)
#'
#' # the standard error of one penalized coordinate under each convention
#' types <- list(frequentist_vcov(), bayesian_vcov(), unconditional_vcov())
#' vapply(types, function(ty) {
#'   sqrt(diag(vcov(fs, type = ty)))[["mu:s(x).z2"]]
#' }, numeric(1))
#'
#' # the class itself is abstract
#' try(vcov_type(name = "mine"))
#' @export
vcov_type <- S7::new_class(
  "vcov_type",
  properties = list(name = S7::class_character),
  abstract = TRUE
)

BayesianVcov <- S7::new_class("BayesianVcov", parent = vcov_type)

FrequentistVcov <- S7::new_class("FrequentistVcov", parent = vcov_type)

UnconditionalVcov <- S7::new_class("UnconditionalVcov", parent = vcov_type)


#' @title Bayesian Variance of a Fit
#'
#' @description
#' Selects the variance \eqn{V_b = (H + S)^{-1}}, the inverse of the
#' penalized information: the information \eqn{H} of the log-likelihood plus
#' the second derivative \eqn{S} of the penalty. It is the posterior variance
#' under the prior whose negative logarithm is the penalty, conditional on the
#' hyperparameters, and it is the default of every function that takes a
#' `type`.
#'
#' @details
#' Choose it for an interval around a penalized term: it carries the
#' smoothing bias as variance, and such intervals cover at about their nominal
#' rate. It does not include the uncertainty of the hyperparameters, for which
#' [unconditional_vcov()] is the convention. See [vcov_type] for the three side
#' by side and [vcov.StatmodFit()] for the matrices in full.
#'
#' The constructor has no arguments.
#'
#' @return A `BayesianVcov` object, inheriting from [vcov_type], with the
#'   property `name` equal to `"bayesian"`.
#'
#' @seealso [frequentist_vcov()] and [unconditional_vcov()] for the siblings,
#'   [vcov_type] for the comparison, and [vcov.StatmodFit()],
#'   [confint.StatmodFit()], [summary.StatmodFit()] and [statmod_test()] for
#'   the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' bayesian_vcov()
#' vcov(fit, type = bayesian_vcov())
#' confint(fit, type = bayesian_vcov())
#' @export
bayesian_vcov <- function() BayesianVcov(name = "bayesian")

#' @title Frequentist Variance of a Fit
#'
#' @description
#' Selects the sandwich variance
#' \eqn{V_f = (H+S)^{-1} H (H+S)^{-1}}, with \eqn{H} the information of the
#' log-likelihood and \eqn{S} the second derivative of the penalty. It is the
#' sampling variance of the penalized estimator at a fixed penalty,
#' conditional on the hyperparameters.
#'
#' @details
#' The matrix is smaller than the Bayesian one wherever there is a penalty,
#' because it omits the smoothing bias, and an interval built from it covers
#' less. With no penalty it equals \eqn{H^{-1}}. See [vcov_type] for the three
#' conventions side by side and [vcov.StatmodFit()] for the matrices in full.
#'
#' The constructor has no arguments.
#'
#' @return A `FrequentistVcov` object, inheriting from [vcov_type], with the
#'   property `name` equal to `"frequentist"`.
#'
#' @seealso [bayesian_vcov()] and [unconditional_vcov()] for the siblings,
#'   [vcov_type] for the comparison, and [vcov.StatmodFit()],
#'   [confint.StatmodFit()], [summary.StatmodFit()] and [statmod_test()] for
#'   the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' frequentist_vcov()
#' vcov(fit, type = frequentist_vcov())
#' summary(fit, type = frequentist_vcov())
#' @export
frequentist_vcov <- function() FrequentistVcov(name = "frequentist")

#' @title Unconditional Variance of a Fit
#'
#' @description
#' Selects the variance \eqn{V' = V_b + J V_\theta J'}, which adds to the
#' Bayesian variance \eqn{V_b} the variance that the estimation of the
#' hyperparameters (smoothing parameters, variances of random effects)
#' contributes. The Bayesian and the frequentist matrices read the
#' hyperparameters at the value where the outer search stopped, as though it
#' were known.
#'
#' @details
#' The correction is positive semi-definite and matters where a smoothing
#' parameter is poorly determined; it affects the penalized coordinates much
#' more than the unpenalized ones. It requires the curvature of the outer
#' criterion, so it costs more than the other two conventions. Where no
#' hyperparameter was estimated by a differentiable criterion the correction
#' is zero and the matrix equals the Bayesian one, and where the curvature
#' cannot be read the conditional matrix is returned with a warning. It is the
#' matrix that mgcv returns for `unconditional = TRUE`. See [vcov_type] for the
#' three conventions side by side and [vcov.StatmodFit()] for the matrices in
#' full.
#'
#' The constructor has no arguments.
#'
#' @return An `UnconditionalVcov` object, inheriting from [vcov_type], with
#'   the property `name` equal to `"unconditional"`.
#'
#' @seealso [bayesian_vcov()] and [frequentist_vcov()] for the siblings,
#'   [vcov_type] for the comparison, and [vcov.StatmodFit()],
#'   [confint.StatmodFit()], [summary.StatmodFit()] and [statmod_test()] for
#'   the functions that take it.
#'
#' @examples
#' set.seed(1)
#' ds <- data.frame(x = runif(200))
#' ds$y <- sin(2 * pi * ds$x) + rnorm(200, sd = 0.3)
#' fs <- statmod(y ~ s(x, bspline_smooth(k = 10)),
#'               distributions7::gaussian1_distrib(), ds)
#' unconditional_vcov()
#' sqrt(diag(vcov(fs, type = unconditional_vcov())))[["mu:s(x).z2"]]
#' @export
unconditional_vcov <- function() UnconditionalVcov(name = "unconditional")


#' @title Print a Variance Convention
#' @name print.vcov_type
#' @description
#' Prints one line, the name of the convention.
#' @param x A [vcov_type] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(unconditional_vcov())
#' @keywords internal
S7::method(print, vcov_type) <- function(x, ...) {
  cat("<vcov> ", x@name, "\n", sep = "")
  invisible(x)
}


#' Check the Argument That Selects a Variance Convention
#'
#' @description
#' Returns `type` if it is a [vcov_type] object. A character string signals an
#' error that names the constructor to use, and any other value signals an
#' error that lists the three constructors.
#'
#' @param type The value of a `type` argument.
#' @param arg The argument's name, for the message.
#'
#' @return `type`.
#'
#' @keywords internal
as_vcov_type <- function(type, arg = "type") {
  if (S7::S7_inherits(type, vcov_type)) return(type)
  ctor <- c(bayesian = "bayesian_vcov()", frequentist = "frequentist_vcov()",
            unconditional = "unconditional_vcov()")
  if (is.character(type) && length(type) >= 1L && !is.na(type[[1L]]) &&
      type[[1L]] %in% names(ctor)) {
    stop(sprintf(paste0("'%s' must be a vcov_type object such as ",
                        "bayesian_vcov(); the string \"%s\" is no longer ",
                        "accepted. Use %s."),
                 arg, type[[1L]], ctor[[type[[1L]]]]), call. = FALSE)
  }
  stop(sprintf(paste0("'%s' must be a vcov_type object: bayesian_vcov(), ",
                      "frequentist_vcov() or unconditional_vcov()."), arg),
       call. = FALSE)
}
