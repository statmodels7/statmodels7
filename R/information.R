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
#' @param method `NULL` or an object that inherits from
#'   [distributions7::ExpectMethod], the method that obtains the expected
#'   information wherever the settled choice or the fallback uses it. `NULL`
#'   selects the closed form of the family where it has one and
#'   [distributions7::expect_opg()] otherwise.
#' @param repair_floor A single positive number, the relative floor of the
#'   repair of an indefinite observed information.
#'
#' @return An S7 object of class `AutoInfo`, which inherits from
#'   [distributions7::Information], with the properties `method` and
#'   `repair_floor`.
#'
#' @examples
#' a <- auto_info(method = distributions7::expect_mc(nsim = 500))
#' S7::S7_inherits(a, distributions7::Information)
#' S7::S7_inherits(a, AutoInfo)
#' a@method
#'
#' @seealso [auto_info()] for the constructor and [information_objects] for the
#'   comparison with [distributions7::observed_info()] and
#'   [distributions7::expected_info()].
#' @export
AutoInfo <- S7::new_class("AutoInfo", parent = distributions7::Information,
  properties = list(
    method = S7::new_union(NULL, distributions7::ExpectMethod),
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
#' inherit from [distributions7::Information].
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
#'     degrees of freedom of [aic()] and [bic()]. The effective degrees of
#'     freedom stored in the fitted model (`fit@edf`) are read on this matrix
#'     whenever an outer criterion ran, and on the matrix of `inner_information`
#'     only for a model without one. The default is
#'     [distributions7::observed_info()].}
#'   \item{`information` of [bayesian_vcov()], [frequentist_vcov()] and
#'     [unconditional_vcov()]}{The matrix inverted for the standard errors.
#'     The default is `NULL`, which selects the matrix by the family of the
#'     fit: the observed information where the log-likelihood is twice
#'     differentiable in every parameter of the family, and
#'     [distributions7::expected_info()] otherwise (see "The default of the
#'     variance objects" below). [summary.StatmodFit()] reads the effective
#'     degrees of freedom it prints on this matrix too, so that they are shares
#'     of the same \eqn{H} as the standard errors beside them.}
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
#'     a function of the parameters alone. Its `method` selects how the
#'     expectation is obtained (see "The method of the expected information"
#'     below).}
#'   \item{[auto_info()]}{A choice that is settled once the family is known.
#'     The rule depends on where the object is passed and is given below.
#'     The fitted model records the settled choice.}
#' }
#'
#' # Which function accepts which object
#'
#' \tabular{lcccc}{
#'   Argument \tab observed \tab expected \tab auto \tab `NULL` \cr
#'   `inner_information` of [statmod()] \tab yes \tab yes \tab yes (default)
#'     \tab no \cr
#'   `information` of [reml()] and [ml()] \tab yes (default) \tab yes \tab yes
#'     \tab no \cr
#'   `information` of [aic()] and [bic()] \tab yes (default) \tab yes \tab no
#'     \tab no \cr
#'   `information` of the three variance objects \tab yes \tab yes \tab no
#'     \tab yes (default) \cr
#' }
#'
#' [aic()], [bic()] and the variance objects reject [auto_info()] because the
#' quantity that they compute is defined on one named matrix.
#'
#' # The method of the expected information
#'
#' The `method` of [distributions7::expected_info()] and of [auto_info()]
#' selects how the expected information, and the derivatives of it that an
#' exact outer gradient reads, are obtained. It is `NULL` or one of
#' [distributions7::expect_opg()], [distributions7::expect_integrate()] and
#' [distributions7::expect_mc()]; the three are compared in
#' [distributions7::expected_derivative_methods].
#'
#' - With `method = NULL`, the default, the closed form of the family is used
#'   where the family writes the expected information out (see
#'   [distributions7::expected_hessian_exact()]), and
#'   [distributions7::expect_opg()] otherwise. The values are those of the
#'   closed form for every shipped family, and the outer criterion under the
#'   default is the same criterion as before the argument existed.
#' - An object is always used, also for a family that has the closed form.
#' - The three arguments above are independent, and so are their methods.
#'
#' In an outer criterion ([reml()], [ml()], [aic()], [bic()]) on the expected
#' information, one method object defines the matrix and every derivative of
#' it that the exact outer gradient and Hessian read. The derivative is the
#' derivative of the matrix that the method defines:
#'
#' - for the closed form, the closed form of the derivative;
#' - for [distributions7::expect_integrate()] and [distributions7::expect_mc()],
#'   the derivative of the expectation, which moves the density as well as the
#'   integrand;
#' - for [distributions7::expect_opg()], the exact derivative of the outer
#'   product of the scores at the data, which treats the data as fixed.
#'
#' Where the method is `NULL` and the family has no closed form for the
#' derivative, the matrix and its derivatives are both
#' [distributions7::expect_opg()]. The outer product of the scores depends on
#' the response, so a criterion built on it is a smooth function of the
#' hyperparameters with exact derivatives, and it is not the Laplace
#' approximation of an expectation: the matrix is an estimate of the expected
#' information, noisier than a computed expectation and about 25 times
#' cheaper. Measured on a density-only family (n = 300, ten basis functions),
#' [reml()] took 0.8 s with `method = NULL` and 25.3 s with
#' `method = expect_integrate(identity = "score")`, and the smoothing
#' parameter was 12.57 and 13.00, against 13.02 on the observed information.
#' The computed expectation is selected with
#' `method = expect_integrate(identity = "score")`. [distributions7::expect_mc()] draws a new
#' sample at each call of the matrix and of each derivative, so that the
#' derivative is that of the expectation and not of the sampled matrix; the
#' outer search is then noisy at the level of the Monte Carlo error.
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
#' # The default of the variance objects
#'
#' The `information` of [bayesian_vcov()], [frequentist_vcov()] and
#' [unconditional_vcov()] is `NULL`. It is resolved for each fit by the family
#' of the fit, in the same way as the `information` of
#' [distributions7::vcov.distrib_fit()]: the matrix is
#' [distributions7::observed_info()] when [distributions7::params_order()] is
#' at least 2 for every parameter of the family, that is, when the
#' log-likelihood is twice differentiable in every parameter, and
#' [distributions7::expected_info()] with `method = NULL` otherwise, for the
#' whole matrix. A family with a kink in a parameter (the Laplace, the elastic
#' net) has a singular or undefined observed information in that parameter,
#' and the closed-form expected information gives finite standard errors.
#' [summary.StatmodFit()] states the matrix whenever it is not the observed
#' information, and a printed variance object states the information that it
#' names.
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
#'   the third object, [AutoInfo] for its class,
#'   [distributions7::expected_derivative_methods] for the methods of the
#'   expected information, [statmod()], [iwls()], [reml()], [ml()], [aic()],
#'   [bic()] and [vcov_method] for the functions that take an object, and
#'   [iwls_solvers] for the factorization of the matrix.
#'
#' @examples
#' distributions7::observed_info()
#' distributions7::expected_info()
#' auto_info(method = distributions7::expect_integrate(identity = "score"))
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
#' @param method `NULL` (the default) or an object that inherits from
#'   [distributions7::ExpectMethod]; it is used as the `method` of the
#'   [distributions7::expected_info()] that the settled choice or the fallback
#'   selects:
#'   - `NULL`: the closed form of the family where it has one, and
#'     [distributions7::expect_opg()] otherwise;
#'   - [distributions7::expect_opg()]: the outer product of the observed
#'     scores, which costs one gradient;
#'   - [distributions7::expect_integrate()]: quadrature (or summation over the
#'     support), with `identity = "hessian"` or `identity = "score"`;
#'   - [distributions7::expect_mc()]: Monte Carlo.
#'
#'   An object is used also for a family that has the closed form.
#'   [distributions7::expected_derivative_methods] compares the methods. Any
#'   other value signals an error.
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
#'   [distributions7::Information], with the properties `method` and
#'   `repair_floor`.
#'
#' @examples
#' auto_info()
#' auto_info(method = distributions7::expect_integrate(identity = "score"),
#'           repair_floor = 1e-10)
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
#'   [information_objects] for the comparison;
#'   [distributions7::expected_derivative_methods] for the methods; [AutoInfo]
#'   for the class; [statmod()], [iwls()], [reml()] and [ml()] for the
#'   functions that take the object.
#' @export
auto_info <- function(method = NULL, repair_floor = 1e-8) {
  check_info_method(method)
  check_positive_number(repair_floor, "repair_floor")
  AutoInfo(method = method, repair_floor = as.numeric(repair_floor))
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


#' Check the Method Given to an Information Object
#'
#' @description
#' Returns `method` invisibly when it is `NULL` or inherits from
#' [distributions7::ExpectMethod] and signals an error otherwise.
#'
#' @param method The value to check.
#'
#' @return `method`, invisibly.
#'
#' @seealso [auto_info()]
#' @keywords internal
check_info_method <- function(method) {
  if (is.null(method) ||
      S7::S7_inherits(method, distributions7::ExpectMethod)) {
    return(invisible(method))
  }
  stop(paste0("'method' must be NULL or an ExpectMethod object such as ",
              "expect_opg(), expect_integrate() or expect_mc()."),
       call. = FALSE)
}


#' Check That an Argument Is an Information Object
#'
#' @description
#' Returns `x` invisibly when it inherits from
#' [distributions7::Information] and signals an error otherwise. With
#' `auto = FALSE`, an [auto_info()] object is rejected as well, and with
#' `null = TRUE`, `NULL` is accepted.
#'
#' @param x The value to check.
#' @param arg The argument name used in the message.
#' @param fn The function whose argument it is, used in the message.
#' @param auto Whether an [auto_info()] object is admissible.
#' @param null Whether `NULL` is admissible.
#'
#' @return `x`, invisibly.
#'
#' @seealso [information_objects]
#' @keywords internal
check_information <- function(x, arg = "information", fn = "reml",
                              auto = TRUE, null = FALSE) {
  if (null && is.null(x)) return(invisible(x))
  ok <- if (auto) "observed_info(), expected_info() or auto_info()"
        else "observed_info() or expected_info()"
  if (null) ok <- paste0("NULL, ", ok)
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
#' `"auto"`). `info_emethod()` returns the method of the expected information
#' that the object carries, `NULL` for the closed form of the family or the
#' outer product of the scores, and `NULL` also for an
#' [distributions7::observed_info()] object.
#'
#' @param x An object inheriting from [distributions7::Information].
#'
#' @return A logical, a string, `NULL` or an object inheriting from
#'   [distributions7::ExpectMethod].
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
info_emethod <- function(x) {
  if (S7::S7_inherits(x, distributions7::ExpectedInfo) ||
      S7::S7_inherits(x, AutoInfo)) x@method
  else NULL
}


#' The Method That a Criterion Reads
#'
#' @description
#' Returns the method of the expected information that the fits made for a
#' criterion use: the `method` of the first of `sparse_method` and `method`
#' whose information is an [distributions7::expected_info()] object, and
#' `emethod`, that of the inner method, when neither is.
#'
#' @param method,sparse_method [OuterMethod()] objects, or `NULL`.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#'
#' @return `NULL` or an object that inherits from
#'   [distributions7::ExpectMethod].
#'
#' @seealso [information_objects]
#' @keywords internal
select_emethod <- function(method, sparse_method, emethod) {
  for (m in list(sparse_method, method)) {
    if (!is.null(m) && info_expected(m@hessian)) return(m@hessian@method)
  }
  emethod
}


# The one-slot memory of resolve_emethod(): the family last probed and the
# answer, compared by identical() (a pointer comparison for the same object).
.emethod_memory <- new.env(parent = emptyenv())

#' The Method That Obtains the Derivatives of the Expected Information
#'
#' @description
#' Returns the method object with which the derivatives of the expected
#' information (the functions
#' [distributions7::distrib_dexpected_hessian()] and
#' [distributions7::distrib_d2expected_hessian()]) are called, so that the
#' derivative is the derivative of the matrix that
#' [distributions7::distrib_expected_hessian()] returns for the same `emethod`.
#'
#' @details
#' An object is returned unchanged. For `NULL`, the matrix is the closed form
#' of the family where the family has one, and the derivatives of the family
#' are those of that closed form, so `NULL` is returned; where the family has
#' no closed form the matrix is [distributions7::expect_opg()] and the
#' derivatives signal an error with `NULL`, so [distributions7::expect_opg()]
#' is returned and used for both. A family is read as having the closed form of
#' its derivative when [distributions7::expected_hessian_exact()] is `TRUE`,
#' and otherwise (a truncated family, whose expected information is a
#' quadrature that the family registers with its derivative) when
#' [distributions7::distrib_dexpected_hessian()] returns a value at a random
#' point; the probe leaves the random stream of the caller as it was and the
#' answer for the last family probed is kept.
#'
#' @param distrib A [distributions7::distrib] object.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#'
#' @return `NULL` or an object that inherits from
#'   [distributions7::ExpectMethod].
#'
#' @seealso [information_objects], [expected_deriv_ok()]
#' @keywords internal
resolve_emethod <- function(distrib, emethod) {
  if (!is.null(emethod)) return(emethod)
  if (S7::S7_inherits(distrib, distributions7::multivariate_distrib)) {
    return(NULL)
  }
  if (isTRUE(tryCatch(distributions7::expected_hessian_exact(distrib),
                      error = function(e) FALSE))) {
    return(NULL)
  }
  if (!is.null(.emethod_memory$distrib) &&
      identical(.emethod_memory$distrib, distrib)) {
    return(.emethod_memory$value)
  }
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(20261026L)
  value <- if (isTRUE(expected_deriv_ok(distrib, NULL))) NULL else
    distributions7::expect_opg()
  .emethod_memory$distrib <- distrib
  .emethod_memory$value <- value
  value
}


#' Print an Information Object Chosen Against the Family
#'
#' @description
#' Writes one line naming the object and the method of the expectation.
#'
#' @param x An [auto_info()] object.
#' @param ... Unused.
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' auto_info()
#' auto_info(method = distributions7::expect_mc(nsim = 500))
#'
#' @seealso [information_objects]
#' @name print.AutoInfo
#' @keywords internal
print.AutoInfo <- function(x, ...) {
  lab <- if (is.null(x@method)) {
    "closed form where the family has one, otherwise the outer product of the scores"
  } else {
    sub("^Expectation by ", "expectation by ",
        paste(utils::capture.output(print(x@method)), collapse = " "))
  }
  cat(sprintf("Information chosen against the family (%s)\n", lab))
  invisible(x)
}
S7::method(print, AutoInfo) <- print.AutoInfo
