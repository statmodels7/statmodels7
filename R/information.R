#' @include spec.R
NULL

#' @title S7 Class for the Information Chosen Against the Family
#' @name AutoInfo
#'
#' @description
#' `AutoInfo` is the class of the object that [auto_info()] returns: an
#' information whose matrix, observed or expected, is settled once the family
#' is known. It inherits from [distributions7::Information], the abstract
#' parent of [distributions7::ObservedInfo] and [distributions7::ExpectedInfo].
#' Use [auto_info()], which validates its arguments; the class is exported so
#' that `S7::S7_inherits()` can test an object. [information_objects] compares
#' the three objects.
#'
#' @param approx An object that inherits from [distributions7::ExpectApprox],
#'   the approximation of the expected information for a family that has no
#'   closed form for it.
#' @param repair_floor A single positive number, the relative floor of the
#'   repair of an indefinite observed information.
#'
#' @return An S7 object of class `AutoInfo`, which inherits from
#'   [distributions7::Information], with the properties `approx` and
#'   `repair_floor`.
#'
#' @examples
#' a <- auto_info(approx = distributions7::expect_mc(nsim = 500))
#' S7::S7_inherits(a, distributions7::Information)
#' S7::S7_inherits(a, AutoInfo)
#' a@approx
#'
#' @seealso [auto_info()] for the constructor and [information_objects] for the
#'   comparison with [distributions7::observed_info()] and
#'   [distributions7::expected_info()].
#' @export
AutoInfo <- S7::new_class("AutoInfo", parent = distributions7::Information,
  properties = list(
    approx = distributions7::ExpectApprox,
    repair_floor = S7::class_numeric
  )
)


#' @title Information Matrices Available to a Fit
#' @name information_objects
#'
#' @description
#' An information object selects the matrix of second derivatives of the
#' log-likelihood that a computation inverts or takes the determinant of.
#' Three objects are available: [distributions7::observed_info()] and
#' [distributions7::expected_info()], which \pkg{distributions7} defines and
#' documents, and [auto_info()], which \pkg{statmodels7} defines. All three
#' inherit from [distributions7::Information]. A character string is not
#' accepted.
#'
#' Three arguments take them, and each one selects the matrix for its own
#' computation, so that the choices are independent:
#'
#' \describe{
#'   \item{`inner_information` of [statmod()]}{The matrix of the inner fit,
#'     which solves for the coefficients at fixed hyperparameters. With
#'     [iwls()] the matrix is inverted at each step; with an
#'     \pkg{optimizers7} optimizer it is the Hessian that the optimizer
#'     receives. The default is [auto_info()].}
#'   \item{`information` of [reml()], [ml()], [aic()] and [bic()]}{The matrix
#'     of the outer criterion: its determinant enters the Laplace
#'     approximation of [reml()] and [ml()], and its trace gives the effective
#'     degrees of freedom of [aic()] and [bic()]. The default is
#'     [distributions7::observed_info()].}
#'   \item{`information` of [bayesian_vcov()], [frequentist_vcov()] and
#'     [unconditional_vcov()]}{The matrix inverted for the standard errors.
#'     The default is [distributions7::observed_info()].}
#' }
#'
#' @details
#' Let \eqn{\ell(\beta)} be the log-likelihood of the coefficients \eqn{\beta}.
#'
#' \describe{
#'   \item{[distributions7::observed_info()]}{The negative Hessian of the
#'     log-likelihood at the data,
#'     \eqn{-\partial^2 \ell / \partial\beta\,\partial\beta^\top}. It is the
#'     curvature of the likelihood in the sample at hand, and every family has
#'     it. It may be indefinite far from the mode, and for a family such as
#'     the Student t at an outlier. It is always the matrix of the data, and
#'     the expected information is never substituted for it. A step that
#'     needs a positive definite matrix floors the eigenvalues of the
#'     penalized matrix (see [solve_chol()]). The observed information of a
#'     coefficient whose second derivative is zero at the estimates, as for
#'     the location of a Laplace response, is singular, and a variance that
#'     inverts it is not available there.}
#'   \item{[distributions7::expected_info()]}{The Fisher information,
#'     \eqn{\mathbb{E}[-\partial^2 \ell / \partial\beta\,\partial\beta^\top]},
#'     a function of the parameters alone. Where the family writes the
#'     expectation out, the closed form is used. Where it does not, the
#'     expectation is approximated by the `approx` of the object (see
#'     [distributions7::expected_derivative_methods]).}
#'   \item{[auto_info()]}{A choice that is settled once the family is known.
#'     The rule depends on where the object is passed and is given below.
#'     The fitted model records the settled choice.}
#' }
#'
#' # Which function accepts which object
#'
#' \tabular{lccc}{
#'   Argument \tab observed \tab expected \tab auto \cr
#'   `inner_information` of [statmod()] \tab yes \tab yes \tab yes (default) \cr
#'   `information` of [reml()] and [ml()] \tab yes (default) \tab yes \tab yes \cr
#'   `information` of [aic()] and [bic()] \tab yes (default) \tab yes \tab no \cr
#'   `information` of the three variance objects \tab yes (default) \tab yes
#'     \tab no \cr
#' }
#'
#' [aic()], [bic()] and the variance objects reject [auto_info()] because the
#' quantity that they compute is defined on one named matrix.
#'
#' # The rule of auto_info()
#'
#' In `inner_information`, the object is settled by the family and, with
#' [iwls()], by the run.
#'
#' - The expected information is chosen where the family writes it exactly
#'   (see [distributions7::expected_hessian_exact()]) and at an acceptable
#'   cost (see [distributions7::expected_hessian_costly()]). The observed
#'   information is chosen for every other family.
#' - With [iwls()] on the expected information, the run continues on the
#'   observed information as soon as Fisher scoring contracts the score
#'   slowly (see `switch_after` and `switch_rate` in [iwls()]).
#' - With [iwls()] on the observed information, the expected information takes
#'   the step where the observed step, after the repair of an indefinite
#'   matrix, finds no acceptable point. The `repair_floor` of [auto_info()]
#'   sets the relative eigenvalue under which the matrix counts as
#'   indefinite.
#' - With an \pkg{optimizers7} optimizer, the settled matrix is the Hessian
#'   that the optimizer receives (the argument `he` of
#'   `optimizers7::minimize()`). Neither the switch nor the fallback exists
#'   there, because both belong to the loop of [iwls()].
#'
#' In `information` of [reml()] and [ml()], the expected information is chosen
#' for a Student t, a skew t or a Cauchy response (wrapped or not), whose
#' observed information is not positive definite at an outlier, and the
#' observed information otherwise.
#'
#' # The approximation
#'
#' `approx` is read only when the family has no closed-form expected
#' information (see [distributions7::expected_hessian_exact()]). It defaults to
#' [distributions7::expect_opg()], the outer product of the scores, which
#' costs one gradient. The three arguments above are independent, and so are
#' their approximations. A criterion on the expected information of a family
#' with no closed form is rejected under the default
#' [distributions7::expect_opg()], because the outer product of the scores at
#' the data is not an expectation; the criterion needs
#' [distributions7::expect_bartlett()] or another approximation that takes an
#' expectation.
#'
#' A model that carries a structural term (a score-driven filter or a mixture
#' over states) is fitted and inverted on the observed information whatever
#' the object says, because those terms have no expected information.
#'
#' @return No value. The topic describes the objects returned by
#'   [distributions7::observed_info()], [distributions7::expected_info()] and
#'   [auto_info()].
#'
#' @seealso [distributions7::observed_info()] and
#'   [distributions7::expected_info()] for the two matrices, [auto_info()] for
#'   the third object, [AutoInfo] for its class, [statmod()], [iwls()],
#'   [reml()], [ml()], [aic()], [bic()] and [vcov_method] for the functions
#'   that take an object, and [iwls_solvers] for the factorization of the
#'   matrix.
#'
#' @examples
#' distributions7::observed_info()
#' distributions7::expected_info()
#' auto_info(approx = distributions7::expect_bartlett())
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_information = distributions7::expected_info(),
#'                outer_criterion = reml(
#'                  information = distributions7::observed_info()))
#' fit@methods$smooth
#'
#' # The standard errors are chosen after the fit, independently of both.
#' sqrt(diag(vcov(fit, vcov_type = bayesian_vcov(
#'   information = distributions7::expected_info()))))
NULL


#' @title Information Chosen Against the Family
#' @name auto_info
#'
#' @description
#' Selects an information that [statmod()] settles once the family is known.
#' It is the default of the argument `inner_information` of [statmod()], and
#' [reml()] and [ml()] accept it as well.
#'
#' In `inner_information`, the expected information is chosen where the family
#' writes it exactly and at an acceptable cost, and the observed information
#' otherwise. With [iwls()], a run on the expected information continues on
#' the observed information as soon as Fisher scoring contracts the score
#' slowly, and a run on the observed information takes the step on the
#' expected information where the repaired observed step finds no acceptable
#' point. With an \pkg{optimizers7} optimizer the settled matrix is the Hessian
#' that the optimizer receives. In [reml()] and [ml()] the expected information
#' is chosen for a Student t, a skew t or a Cauchy response (wrapped or not),
#' whose observed information is not positive definite at an outlier, and the
#' observed information is chosen for every other family. The fitted model
#' records the settled choice.
#'
#' [aic()], [bic()] and the variance objects [bayesian_vcov()],
#' [frequentist_vcov()] and [unconditional_vcov()] do not accept it, since the
#' quantity that each of them computes is defined on one named matrix.
#' [information_objects] compares the three information objects.
#'
#' @param approx The approximation of the expected information, as in
#'   [distributions7::expected_info()]. It is read only where the settled
#'   choice or the fallback uses the expected information of a family that
#'   has no closed form for it. The default is [distributions7::expect_opg()].
#' @param repair_floor A single positive number. A penalized observed
#'   information whose smallest eigenvalue is below `repair_floor` times the
#'   largest absolute eigenvalue is read as indefinite by [iwls()]. The step
#'   on it then floors the eigenvalues (the floor of the solver, see
#'   [solve_chol()]), the expected information takes the step only if the
#'   search on the repaired matrix fails, and a run on it does not take the
#'   extra step that follows convergence. The default is `1e-8`. The argument
#'   is read by [iwls()] only.
#'
#' @return An S7 object of class `AutoInfo`, which inherits from
#'   [distributions7::Information], with the properties `approx` and
#'   `repair_floor`.
#'
#' @examples
#' auto_info()
#' auto_info(approx = distributions7::expect_bartlett(), repair_floor = 1e-10)
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_information = auto_info(),
#'                outer_criterion = reml(information = auto_info()))
#' fit@methods$smooth
#'
#' @seealso [distributions7::observed_info()] and
#'   [distributions7::expected_info()] for the objects that name a matrix;
#'   [information_objects] for the comparison; [AutoInfo] for the class;
#'   [statmod()], [iwls()], [reml()] and [ml()] for the functions that take
#'   the object.
#' @export
auto_info <- function(approx = distributions7::expect_opg(),
                      repair_floor = 1e-8) {
  check_info_approx(approx)
  check_positive_number(repair_floor, "repair_floor")
  AutoInfo(approx = approx, repair_floor = as.numeric(repair_floor))
}


#' Check a Single Positive Number
#'
#' @description
#' Returns `x` invisibly when it is a single finite number above zero and
#' signals an error that names `arg` otherwise.
#'
#' @param x The value to check.
#' @param arg The argument name used in the message.
#'
#' @return `x`, invisibly.
#'
#' @seealso [auto_info()], [iwls()], [solve_qr()]
#' @keywords internal
check_positive_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be a single positive number.", arg), call. = FALSE)
  }
  invisible(x)
}


#' Check the Approximation Given to an Information Object
#'
#' @description
#' Returns `approx` invisibly when it inherits from
#' [distributions7::ExpectApprox] and signals an error otherwise.
#'
#' @param approx The value to check.
#'
#' @return `approx`, invisibly.
#'
#' @seealso [auto_info()]
#' @keywords internal
check_info_approx <- function(approx) {
  if (S7::S7_inherits(approx, distributions7::ExpectApprox)) {
    return(invisible(approx))
  }
  stop(paste0("'approx' must be an expectation object: expect_opg(), ",
              "expect_bartlett(), expect_integrate() or expect_mc(nsim = ...)."),
       call. = FALSE)
}


#' Check That an Argument Is an Information Object
#'
#' @description
#' Returns `x` invisibly when it inherits from
#' [distributions7::Information] and signals an error otherwise. With
#' `auto = FALSE`, an [auto_info()] object is rejected as well.
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
check_information <- function(x, arg = "information", fn = "reml",
                              auto = TRUE) {
  ok <- if (auto) "observed_info(), expected_info() or auto_info()"
        else "observed_info() or expected_info()"
  if (S7::S7_inherits(x, distributions7::Information)) {
    if (!auto && S7::S7_inherits(x, AutoInfo)) {
      stop(sprintf(paste0("'%s' of %s() must name one matrix: %s. ",
                          "auto_info() is not admissible here."),
                   arg, fn, ok), call. = FALSE)
    }
    return(invisible(x))
  }
  stop(sprintf("'%s' must be an information object: %s.", arg, ok),
       call. = FALSE)
}


#' Predicates on an Information Object
#'
#' @description
#' `info_expected()` returns `TRUE` for an [distributions7::expected_info()]
#' object and `info_auto()` for an [auto_info()] object. `info_label()`
#' returns the word that names the matrix (`"observed"`, `"expected"` or
#' `"auto"`). `info_approx()` returns the approximation that the object
#' carries, and [distributions7::expect_opg()] for an
#' [distributions7::observed_info()] object.
#'
#' @param x An object inheriting from [distributions7::Information].
#'
#' @return A logical, a string or an object inheriting from
#'   [distributions7::ExpectApprox].
#'
#' @seealso [information_objects]
#' @keywords internal
info_expected <- function(x) S7::S7_inherits(x, distributions7::ExpectedInfo)

#' @rdname info_expected
#' @keywords internal
info_auto <- function(x) S7::S7_inherits(x, AutoInfo)

#' @rdname info_expected
#' @keywords internal
info_label <- function(x) {
  if (S7::S7_inherits(x, distributions7::ExpectedInfo)) "expected"
  else if (S7::S7_inherits(x, AutoInfo)) "auto"
  else "observed"
}

#' @rdname info_expected
#' @keywords internal
info_approx <- function(x) {
  if (S7::S7_inherits(x, distributions7::ExpectedInfo) ||
      S7::S7_inherits(x, AutoInfo)) x@approx
  else distributions7::expect_opg()
}


#' The Approximation That a Criterion Reads
#'
#' @description
#' Returns the approximation of the expected information that the fits made
#' for a criterion use: the `approx` of the first of `sparse_method` and
#' `method` whose information is an [distributions7::expected_info()] object,
#' and `approx`, that of the inner method, when neither is.
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


#' Print an Information Object Chosen Against the Family
#'
#' @description
#' Writes one line naming the object and the approximation of the expectation.
#'
#' @param x An [auto_info()] object.
#' @param ... Unused.
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' auto_info()
#' auto_info(approx = distributions7::expect_mc(nsim = 500))
#'
#' @seealso [information_objects]
#' @name print.AutoInfo
#' @keywords internal
print.AutoInfo <- function(x, ...) {
  lab <- paste(utils::capture.output(print(x@approx)), collapse = " ")
  cat(sprintf("Information chosen against the family (%s)\n",
              sub("^Expectation by ", "expectation by ", lab)))
  invisible(x)
}
S7::method(print, AutoInfo) <- print.AutoInfo
