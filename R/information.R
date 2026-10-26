#' @include spec.R
NULL

#' @title S7 Classes for the Information Matrix of a Fit
#' @name Information
#' @aliases ObservedInfo ExpectedInfo AutoInfo
#'
#' @description
#' `Information` is the abstract parent of the objects that select which
#' matrix a fit inverts or takes the determinant of: the observed Hessian of
#' the log-likelihood, the Fisher information, or a choice between the two that
#' is settled once the family is known. Its three subclasses are
#' `ObservedInfo`, `ExpectedInfo` and `AutoInfo`. They are built by
#' [observed_info()], [expected_info()] and [auto_info()], which are the
#' functions to call; the classes are exported so that `S7::S7_inherits()` can
#' test an object. [information_objects] compares the three.
#'
#' `ObservedInfo` has no properties. `ExpectedInfo` and `AutoInfo` have one,
#' `approx`, an object that inherits from [distributions7::ExpectApprox] and
#' selects how the expected information is computed for a family that has no
#' closed form for it. The raw constructors do not validate `approx`; the
#' functions [expected_info()] and [auto_info()] do.
#'
#' @param approx An object inheriting from [distributions7::ExpectApprox], for
#'   `ExpectedInfo` and `AutoInfo`.
#'
#' @return An S7 object of the class named, which inherits from `Information`.
#'
#' @examples
#' a <- expected_info(approx = distributions7::expect_mc(nsim = 500))
#' S7::S7_inherits(a, Information)
#' S7::S7_inherits(a, ExpectedInfo)
#' a@approx
#'
#' @seealso [information_objects] for the comparison of the three
#'   constructors; [observed_info()], [expected_info()] and [auto_info()] for
#'   the constructors.
#' @export
Information <- S7::new_class("Information", abstract = TRUE)

#' @rdname Information
#' @export
ObservedInfo <- S7::new_class("ObservedInfo", parent = Information)

#' @rdname Information
#' @export
ExpectedInfo <- S7::new_class("ExpectedInfo", parent = Information,
  properties = list(approx = distributions7::ExpectApprox)
)

#' @rdname Information
#' @export
AutoInfo <- S7::new_class("AutoInfo", parent = Information,
  properties = list(approx = distributions7::ExpectApprox)
)


#' @title Information Matrices Available to a Fit
#' @name information_objects
#'
#' @description
#' [iwls()], [reml()], [ml()], [aic()] and [bic()] take an argument
#' `information`, an object that selects the matrix used for a Newton-type step
#' or for a Laplace determinant. Three objects are available, built by
#' [observed_info()], [expected_info()] and [auto_info()]; all three inherit
#' from [Information]. A character string is not accepted.
#'
#' @details
#' Let \eqn{\ell(\beta)} be the log-likelihood of the coefficients \eqn{\beta},
#' and \eqn{S} the Hessian of the penalty.
#'
#' \describe{
#'   \item{[observed_info()]}{The negative Hessian of the log-likelihood at the
#'     data, \eqn{-\partial^2 \ell / \partial\beta\,\partial\beta^\top}. In
#'     [iwls()] this is Newton's method, which converges quadratically near
#'     the mode. The matrix may be indefinite far from the mode, and for a
#'     family such as the Student t at an outlier. An observed penalized
#'     information that is not positive definite is repaired by flooring its
#'     eigenvalues before the step is taken. The exact outer gradient and
#'     Hessian of [reml()] and [ml()] are available on this matrix.}
#'   \item{[expected_info()]}{The Fisher information,
#'     \eqn{\mathbb{E}[-\partial^2 \ell / \partial\beta\,\partial\beta^\top]},
#'     a function of the parameters alone. In [iwls()] this is Fisher scoring.
#'     The matrix is positive definite wherever the family is regular, so the
#'     step is usable from a poor start. Where the family writes the
#'     expectation out, the closed form is used. Where it does not, the
#'     expectation is approximated by the object in `approx`, one of
#'     [distributions7::expect_opg()], [distributions7::expect_bartlett()],
#'     [distributions7::expect_integrate()] and
#'     [distributions7::expect_mc()]; see
#'     [distributions7::expected_derivative_methods].}
#'   \item{[auto_info()]}{A choice that [statmod()] makes once the family is
#'     known. In [iwls()], the expected information is chosen where the family
#'     writes it exactly and at an acceptable cost, and the observed
#'     information otherwise, with the expected one taking the step where the
#'     observed step finds no acceptable point. In [reml()] and [ml()], the
#'     expected information is chosen for a Student t, a skew t or a Cauchy
#'     response (wrapped or not), whose observed information is not positive
#'     definite at an outlier, and the observed information otherwise. The
#'     object carries an `approx`, which is the approximation used by the
#'     expected information whenever the choice or a fallback reaches it.
#'     [aic()] and [bic()] do not accept it.}
#' }
#'
#' # Which function uses which matrix
#'
#' \describe{
#'   \item{[iwls()]}{The matrix inverted at each step of the inner fit, with
#'     the penalty added. The default is [auto_info()].}
#'   \item{[reml()] and [ml()]}{The matrix whose determinant enters the Laplace
#'     approximation of the marginal likelihood, and whose derivatives enter
#'     the outer gradient. The default is [observed_info()].}
#'   \item{[aic()] and [bic()]}{The matrix whose trace gives the effective
#'     degrees of freedom of the fit. The default is [observed_info()]; only
#'     [observed_info()] and [expected_info()] are admissible.}
#' }
#'
#' # The approximation
#'
#' `approx` is read only when the family has no closed-form expected
#' information (see [distributions7::expected_hessian_exact()]). It defaults to
#' [distributions7::expect_opg()], the outer product of the scores, which costs
#' one gradient. In an [iwls()] object `approx` applies to the inner fit. In
#' [reml()], [ml()], [aic()] and [bic()] it applies to the criterion, so that
#' the two objects are independent. A criterion on the expected information of
#' a family with no closed form is rejected under the default
#' [distributions7::expect_opg()], because the outer product of the scores at
#' the data is not an expectation; the criterion needs
#' [distributions7::expect_bartlett()] or another approximation that takes an
#' expectation.
#'
#' @return No value. The topic describes the three objects, of classes
#'   [ObservedInfo][Information], [ExpectedInfo][Information] and
#'   [AutoInfo][Information], that [observed_info()], [expected_info()] and
#'   [auto_info()] return.
#'
#' @seealso [observed_info()], [expected_info()] and [auto_info()] for the
#'   constructors; [iwls()], [reml()], [ml()], [aic()] and [bic()] for the
#'   functions that take the object; [Information] for the classes.
#'
#' @examples
#' observed_info()
#' expected_info()
#' auto_info(approx = distributions7::expect_bartlett())
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(information = expected_info()),
#'                outer_criterion = reml(information = observed_info()))
#' fit@methods$smooth
NULL


#' @title The Observed Information
#' @name observed_info
#'
#' @description
#' Selects the observed information, the negative Hessian of the
#' log-likelihood at the data,
#' \deqn{-\frac{\partial^2 \ell}{\partial\beta\,\partial\beta^\top},}
#' as the matrix of a fit. In [iwls()] it gives Newton's method: the step
#' solves the system built on this matrix plus the penalty, and the iteration
#' converges quadratically near the mode. The matrix may fail to be positive
#' definite far from the mode; in that case its penalized eigenvalues are
#' floored before the step. In [reml()] and [ml()] the determinant of the
#' matrix enters the Laplace approximation, which is then the approximation of
#' the marginal likelihood itself, and the outer gradient and Hessian are exact
#' on it. In [aic()] and [bic()] its trace gives the effective degrees of
#' freedom.
#'
#' It is the default of [reml()], [ml()], [aic()] and [bic()].
#' [information_objects] compares it with the other two objects.
#'
#' @return An S7 object of class `ObservedInfo`, which inherits from
#'   [Information], with no properties.
#'
#' @examples
#' observed_info()
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(information = observed_info()),
#'                outer_criterion = reml(information = observed_info()))
#' hyper(fit)
#'
#' @seealso [expected_info()] and [auto_info()] for the other objects;
#'   [information_objects] for the comparison; [iwls()], [reml()], [ml()],
#'   [aic()] and [bic()] for the functions that take the object.
#' @export
observed_info <- function() ObservedInfo()


#' @title The Expected (Fisher) Information
#' @name expected_info
#'
#' @description
#' Selects the expected information, the Fisher information
#' \deqn{\mathbb{E}\left[-\frac{\partial^2 \ell}{\partial\beta\,\partial\beta^\top}
#' \right],}
#' as the matrix of a fit. In [iwls()] it gives Fisher scoring: the matrix is
#' positive definite wherever the family is regular, so the step is usable from
#' a poor start, and the iteration converges linearly near the mode wherever
#' the expected and the observed information differ there. In [reml()] and
#' [ml()] the determinant of the matrix enters the Laplace approximation, which
#' is then a Fisher variant of the marginal likelihood. In [aic()] and [bic()]
#' its trace gives the effective degrees of freedom.
#'
#' Where the family writes the expected information out, the closed form is
#' used and `approx` is not read. Where it does not (see
#' [distributions7::expected_hessian_exact()]), the expectation is approximated
#' by `approx`.
#'
#' @param approx An object that selects the approximation of the expectation:
#'   [distributions7::expect_opg()] (the outer product of the scores, the
#'   default), [distributions7::expect_bartlett()] (the Bartlett identity),
#'   [distributions7::expect_integrate()] (integration against the density) or
#'   [distributions7::expect_mc()] (Monte Carlo, with its own `nsim`). The
#'   comparison is in [distributions7::expected_derivative_methods]. A
#'   character string, or any object that does not inherit from
#'   [distributions7::ExpectApprox], signals an error.
#'
#' @return An S7 object of class `ExpectedInfo`, which inherits from
#'   [Information], with the property `approx`.
#'
#' @examples
#' expected_info()
#' expected_info(approx = distributions7::expect_mc(nsim = 2000))
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(information = expected_info()),
#'                outer_criterion = reml(information = expected_info()))
#' hyper(fit)
#'
#' @seealso [observed_info()] and [auto_info()] for the other objects;
#'   [information_objects] for the comparison; [iwls()], [reml()], [ml()],
#'   [aic()] and [bic()] for the functions that take the object.
#' @export
expected_info <- function(approx = distributions7::expect_opg()) {
  check_info_approx(approx)
  ExpectedInfo(approx = approx)
}


#' @title Information Chosen Against the Family
#' @name auto_info
#'
#' @description
#' Selects an information that [statmod()] settles once the family is known.
#' In [iwls()] the expected information is chosen where the family writes it
#' exactly and at an acceptable cost, and a run on it continues on the observed
#' information as soon as Fisher scoring contracts the score slowly. Otherwise
#' the observed information is chosen, and the expected one takes the step
#' where the observed step, after the repair of an indefinite matrix, finds no
#' acceptable point. In [reml()] and [ml()] the expected information is chosen
#' for a Student t, a skew t or a Cauchy response (wrapped or not), whose
#' observed information is not positive definite at an outlier, and the
#' observed information is chosen for every other family. The fitted model
#' records the settled choice.
#'
#' It is the default of [iwls()]. [aic()] and [bic()] do not accept it, since
#' the degrees of freedom they charge are defined on one named matrix.
#'
#' @param approx The approximation of the expected information, as in
#'   [expected_info()]. It is read only where the settled choice or the
#'   fallback uses the expected information of a family that has no closed form
#'   for it. The default is [distributions7::expect_opg()].
#'
#' @return An S7 object of class `AutoInfo`, which inherits from
#'   [Information], with the property `approx`.
#'
#' @examples
#' auto_info()
#' auto_info(approx = distributions7::expect_bartlett())
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(information = auto_info()),
#'                outer_criterion = reml(information = auto_info()))
#' fit@methods$smooth
#'
#' @seealso [observed_info()] and [expected_info()] for the objects that name
#'   a matrix; [information_objects] for the comparison; [iwls()], [reml()] and
#'   [ml()] for the functions that take the object.
#' @export
auto_info <- function(approx = distributions7::expect_opg()) {
  check_info_approx(approx)
  AutoInfo(approx = approx)
}


#' Check the Approximation Given to an Information Object
#'
#' @description
#' Returns `approx` invisibly when it inherits from
#' [distributions7::ExpectApprox] and signals an error otherwise. A character
#' string of the former interface produces a message that names the
#' constructor that replaces it.
#'
#' @param approx The value to check.
#'
#' @return `approx`, invisibly.
#'
#' @seealso [expected_info()], [auto_info()]
#' @keywords internal
check_info_approx <- function(approx) {
  if (S7::S7_inherits(approx, distributions7::ExpectApprox)) {
    return(invisible(approx))
  }
  repl <- c(opg = "expect_opg()", bartlett = "expect_bartlett()",
            integrate = "expect_integrate()", mc = "expect_mc(nsim = ...)")
  if (is.character(approx) && length(approx) == 1L && !is.na(approx) &&
      approx %in% names(repl)) {
    stop(sprintf(paste0(
      "'approx' must be an expectation object such as expect_opg(); the ",
      "string \"%s\" is no longer accepted. Use %s."),
      approx, repl[[approx]]), call. = FALSE)
  }
  stop(paste0("'approx' must be an expectation object: expect_opg(), ",
              "expect_bartlett(), expect_integrate() or expect_mc(nsim = ...)."),
       call. = FALSE)
}


#' Check That an Argument Is an Information Object
#'
#' @description
#' Returns `x` invisibly when it inherits from [Information] and signals an
#' error otherwise. A character string of the former interface
#' (`"observed"`, `"expected"` or `"auto"`) produces a message that names the
#' constructor that replaces it. With `auto = FALSE`, an [auto_info()] object is
#' rejected as well.
#'
#' @param x The value to check.
#' @param arg The argument name used in the message.
#' @param fn The function whose argument it is, used in the message.
#' @param auto Whether an [auto_info()] object is admissible.
#'
#' @return `x`, invisibly.
#'
#' @seealso [information_objects]
#' @keywords internal
check_information <- function(x, arg = "information", fn = "iwls",
                              auto = TRUE) {
  ok <- if (auto) "observed_info(), expected_info() or auto_info()"
        else "observed_info() or expected_info()"
  if (S7::S7_inherits(x, Information)) {
    if (!auto && S7::S7_inherits(x, AutoInfo)) {
      stop(sprintf(paste0("'%s' of %s() must name one matrix: %s. ",
                          "auto_info() is not admissible here."),
                   arg, fn, ok), call. = FALSE)
    }
    return(invisible(x))
  }
  repl <- c(observed = "observed_info()", expected = "expected_info()",
            auto = "auto_info()")
  if (is.character(x) && length(x) >= 1L && !anyNA(x) && x[[1L]] %in% names(repl)) {
    stop(sprintf(paste0("'%s' must be an information object such as ",
                        "observed_info(); the string \"%s\" is no longer ",
                        "accepted. Use %s."),
                 arg, x[[1L]],
                 if (!auto && identical(x[[1L]], "auto")) "observed_info() or expected_info()"
                 else repl[[x[[1L]]]]), call. = FALSE)
  }
  stop(sprintf("'%s' must be an information object: %s.", arg, ok),
       call. = FALSE)
}


#' Predicates on an Information Object
#'
#' @description
#' `info_expected()` returns `TRUE` for an [expected_info()] object and
#' `info_auto()` for an [auto_info()] object. `info_label()` returns the word
#' that names the matrix (`"observed"`, `"expected"` or `"auto"`).
#' `info_approx()` returns the approximation that the object carries, and
#' [distributions7::expect_opg()] for an [observed_info()] object.
#'
#' @param x An object inheriting from [Information].
#'
#' @return A logical, a string or an object inheriting from
#'   [distributions7::ExpectApprox].
#'
#' @seealso [information_objects]
#' @keywords internal
info_expected <- function(x) S7::S7_inherits(x, ExpectedInfo)

#' @rdname info_expected
#' @keywords internal
info_auto <- function(x) S7::S7_inherits(x, AutoInfo)

#' @rdname info_expected
#' @keywords internal
info_label <- function(x) {
  if (S7::S7_inherits(x, ExpectedInfo)) "expected"
  else if (S7::S7_inherits(x, AutoInfo)) "auto"
  else "observed"
}

#' @rdname info_expected
#' @keywords internal
info_approx <- function(x) {
  if (S7::S7_inherits(x, ExpectedInfo) || S7::S7_inherits(x, AutoInfo)) x@approx
  else distributions7::expect_opg()
}


#' The Approximation That a Criterion Reads
#'
#' @description
#' Returns the approximation of the expected information that the fits made
#' for a criterion use: the `approx` of the first of `sparse_method` and
#' `method` whose information is an [expected_info()] object, and `approx`,
#' that of the inner method, when neither is.
#'
#' @param method,sparse_method [OuterMethod()] objects, or `NULL`.
#' @param approx The approximation of the inner method.
#'
#' @return An object that inherits from [distributions7::ExpectApprox].
#'
#' @seealso [information_objects]
#' @keywords internal
select_approx <- function(method, sparse_method, approx) {
  for (m in list(sparse_method, method)) {
    if (!is.null(m) && info_expected(m@hessian)) return(m@hessian@approx)
  }
  approx
}


#' Print an Information Object
#'
#' @description
#' Writes one line naming the matrix and, for [expected_info()] and
#' [auto_info()], the approximation of the expectation.
#'
#' @param x An object inheriting from [Information].
#' @param ... Unused.
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' observed_info()
#' expected_info(approx = distributions7::expect_mc(nsim = 500))
#'
#' @seealso [information_objects]
#' @name print.Information
#' @keywords internal
print.Information <- function(x, ...) {
  cat(switch(info_label(x),
             observed = "Observed information",
             expected = "Expected information",
             auto = "Information chosen against the family"))
  if (!S7::S7_inherits(x, ObservedInfo)) {
    lab <- paste(utils::capture.output(print(x@approx)), collapse = " ")
    cat(sprintf(" (%s)", sub("^Expectation by ", "expectation by ", lab)))
  }
  cat("\n")
  invisible(x)
}
S7::method(print, Information) <- print.Information
