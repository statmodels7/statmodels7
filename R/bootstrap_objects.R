#' @include inference_objects.R
NULL

# The objects of the bootstrap variance: the way the data are resampled, the
# variance convention that reads the replicates, and the stored replicates.


#' @title Resampling Schemes for the Bootstrap of a Fit
#' @name resample_method
#'
#' @description
#' Three schemes generate the data sets from which [add_bootstrap()] refits a
#' model. [parametric_resample()] draws a new response from the fitted model,
#' [case_resample()] draws observations with replacement, and
#' [cluster_resample()] draws groups of observations with replacement. Each is
#' an object, and the object is passed as the argument `resample` of
#' [add_bootstrap()], [bootstrap_vcov()] and [statmod_bootstrap()].
#'
#' `resample_method` is the abstract parent of the three classes. It cannot be
#' instantiated, and serves to test membership with `S7::S7_inherits()`.
#'
#' @details
#' The schemes answer two different questions. The parametric scheme
#' simulates from the fitted model, so its replicates have the variance that
#' the model attributes to the estimator, and the bootstrap variance
#' converges to the model-based one (to [frequentist_vcov()] where the
#' penalties and hyperparameters are held). The case and cluster schemes
#' resample the data, so they do not depend on the variance specification of
#' the model: with a heteroscedastic response fitted by a homoscedastic model
#' the case bootstrap converges to the sandwich variance, which the parametric
#' bootstrap does not.
#'
#' A resampled data set is represented by frequency weights, which multiply
#' the weights of the fit. The design is not rebuilt, so the knots of a basis
#' and the levels of a factor stay those of the fit, and a refit on the
#' weights is the refit on the duplicated rows with the basis held.
#'
#' @param name A single string, the short name of the scheme.
#'
#' @return Nothing: the class is abstract and its constructor signals an
#'   error. Every object that inherits from it has the property `name`.
#'
#' @seealso [parametric_resample()], [case_resample()] and
#'   [cluster_resample()] for the three schemes, [add_bootstrap()] for the
#'   function that takes one.
#'
#' @examples
#' S7::S7_inherits(case_resample(), resample_method)
#' try(resample_method(name = "mine"))
#' @export
resample_method <- S7::new_class(
  "resample_method",
  properties = list(name = S7::class_character),
  abstract = TRUE
)

#' @title S7 Classes for the Shipped Resampling Schemes
#'
#' @description
#' These are the classes that [parametric_resample()], [case_resample()] and
#' [cluster_resample()] instantiate. Each inherits from [resample_method]. Use
#' the constructors, which set the properties.
#'
#' @param name A single string, the short name of the scheme.
#' @param by For a cluster scheme, the one-sided formula that names the
#'   variable that identifies the groups.
#'
#' @return An S7 object inheriting from [resample_method], with the property
#'   `name`, and `by` for a cluster scheme.
#'
#' @seealso [resample_method] for the comparison of the schemes.
#'
#' @examples
#' class(case_resample())
#' cluster_resample(~ g)@by
#'
#' @name ParametricResample
#' @aliases ParametricResample CaseResample ClusterResample
#' @keywords internal
ParametricResample <- S7::new_class("ParametricResample",
                                    parent = resample_method)

#' @rdname ParametricResample
#' @keywords internal
CaseResample <- S7::new_class("CaseResample", parent = resample_method)

#' @rdname ParametricResample
#' @keywords internal
ClusterResample <- S7::new_class("ClusterResample", parent = resample_method,
  properties = list(by = S7::class_any))


#' @title Parametric Resampling
#'
#' @description
#' Selects the parametric bootstrap: every replicate draws a new response from
#' the family at the fitted parameters and refits the model on it.
#'
#' @details
#' The effects of every [modelterms7::random()] term are drawn afresh from the
#' prior that the fit estimated, so the groups of the data are new groups in
#' every replicate. The replicate of a random-effect coordinate is stored as
#' \eqn{\hat b^* - b^* + \hat b}, the estimate of the replicate minus the
#' effect that was drawn, added to the estimate of the fit. The covariance of
#' these replicates is the variance of the prediction error
#' \eqn{\mathrm{Var}(\hat b - b)}, which is the quantity that
#' [bayesian_vcov()] reports for those coordinates. The covariance of
#' \eqn{\hat b^*} alone would be the variance of the shrunken estimate, which
#' is smaller wherever the shrinkage is strong.
#'
#' The parametric bootstrap reproduces the variance that the model attributes
#' to the estimator. It does not detect a misspecified variance function; the
#' scheme of [case_resample()] does. A fit with a structural term
#' ([modelterms7::gas()], [modelterms7::regime()]) or a block that moves with
#' its coefficients ([modelterms7::seg()], [modelterms7::jump()],
#' [modelterms7::jseg()], [modelterms7::nl()]) is rejected.
#'
#' @return A `ParametricResample` object, inheriting from [resample_method].
#'
#' @seealso [resample_method] for the comparison of the three schemes,
#'   [add_bootstrap()] and [bootstrap_vcov()] for the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' fit <- add_bootstrap(fit, n = 30, resample = parametric_resample())
#' statmod_bootstrap(fit)
#' @export
parametric_resample <- function() ParametricResample(name = "parametric")

#' @title Case Resampling
#'
#' @description
#' Selects the case bootstrap: every replicate draws the observations with
#' replacement and refits the model on them.
#'
#' @details
#' The draw is a vector of counts, how many times each observation is drawn,
#' which multiplies the weights of the fit. The design is not rebuilt, so the
#' knots of a basis stay those of the fit. A replicate in which a coefficient
#' is not identified on the observations of positive weight (a level of a
#' factor that the draw left out) fails and is replaced.
#'
#' The case bootstrap does not depend on the variance specification of the
#' model, and its variance converges to the sandwich variance. It breaks the
#' groups of a [modelterms7::random()] term, and a fit with one is rejected:
#' [cluster_resample()] resamples the groups. A fit with a structural term or
#' a block that moves with its coefficients is rejected as well, and a
#' [cv()] criterion is rejected with `refit = "full"`.
#'
#' @return A `CaseResample` object, inheriting from [resample_method].
#'
#' @seealso [resample_method] for the comparison of the three schemes,
#'   [add_bootstrap()] and [bootstrap_vcov()] for the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' fit <- add_bootstrap(fit, n = 30, resample = case_resample())
#' sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov(resample = case_resample()))))
#' @export
case_resample <- function() CaseResample(name = "case")

#' @title Cluster Resampling
#'
#' @description
#' Selects the cluster bootstrap: every replicate draws the groups of
#' observations with replacement and refits the model on them.
#'
#' @details
#' The groups are the values of the variable that the formula names. The
#' variable is evaluated in the data that the fit was called with (the `data`
#' argument of the call of [statmod()], found in the environment of the model
#' formula) and in the environment of the formula, and it must have one value
#' per observation. A group drawn \eqn{k} times gives all its observations the
#' count \eqn{k}. Without random effects the likelihood is a sum over
#' observations, so this is the fit to the duplicated groups with the design
#' held.
#'
#' In this version a fit with a [modelterms7::random()] term is rejected. The
#' coordinates of the individual random effects have no counterpart between
#' replicates, because a group can be absent or appear twice.
#'
#' @param by A one-sided formula with one variable, as in `~ id`.
#'
#' @return A `ClusterResample` object, inheriting from [resample_method], with
#'   the property `by`.
#'
#' @seealso [resample_method] for the comparison of the three schemes,
#'   [add_bootstrap()] and [bootstrap_vcov()] for the functions that take it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(g = rep(1:12, each = 5), x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(12)[dd$g] + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' fit <- add_bootstrap(fit, n = 30, resample = cluster_resample(~ g))
#' statmod_bootstrap(fit)
#' @export
cluster_resample <- function(by) {
  if (!inherits(by, "formula") || length(by) != 2L ||
      length(all.vars(by)) != 1L) {
    stop("'by' must be a one-sided formula with one variable, as in ~ id.",
         call. = FALSE)
  }
  environment(by) <- baseenv()
  ClusterResample(name = "cluster", by = by)
}

#' @title Print a Resampling Scheme
#' @name print.resample_method
#' @description Prints one line, the name of the scheme and, for a cluster
#'   scheme, its formula.
#' @param x A [resample_method] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(cluster_resample(~ id))
#' @keywords internal
S7::method(print, resample_method) <- function(x, ...) {
  cat("<resample> ", resample_label(x), "\n", sep = "")
  invisible(x)
}

#' The Words That Name a Resampling Scheme
#'
#' @description
#' Returns the name of the scheme, followed by the formula of a cluster
#' scheme. It is also the key by which stored bootstrap entries are matched.
#'
#' @param x A [resample_method] object.
#'
#' @return A single string.
#'
#' @keywords internal
resample_label <- function(x) {
  if (S7::S7_inherits(x, ClusterResample)) {
    paste0(x@name, " (", paste(deparse(x@by), collapse = ""), ")")
  } else x@name
}

#' Check the Argument That Selects a Resampling Scheme
#'
#' @description
#' Returns `x` if it is a [resample_method] object and signals an error that
#' lists the constructors otherwise.
#'
#' @param x The value of a `resample` argument.
#'
#' @return `x`.
#'
#' @keywords internal
as_resample_method <- function(x) {
  if (S7::S7_inherits(x, resample_method)) return(x)
  stop("'resample' must be a resample_method object: parametric_resample(), ",
       "case_resample() or cluster_resample().", call. = FALSE)
}


#' @rdname BayesianVcov
#' @keywords internal
BootstrapVcov <- S7::new_class("BootstrapVcov", parent = vcov_method,
  properties = list(n = S7::class_any, resample = S7::class_any,
                    refit = S7::class_any))

#' @title Bootstrap Variance of a Fit
#'
#' @description
#' Selects the variance that is the covariance of the coefficients of
#' refitted, resampled data sets. The matrix is not read from the curvature at
#' the mode: it is the empirical covariance of the replicates that
#' [add_bootstrap()] stores in the fit.
#'
#' @details
#' The arguments that are `NULL` are resolved against the replicates stored in
#' the fit. A `NULL` argument takes the value of the stored entry if exactly
#' one entry matches the arguments that were given, and the default otherwise
#' (`n = 200`, [parametric_resample()], `"coefficients"`). If several entries
#' match, an error prints the table of entries. A setting that matches no
#' stored entry is computed at the call and stored in the fit, with a message
#' that says how many replicates were computed. `n` is the total that is
#' wanted: with an entry of \eqn{m} replicates, \eqn{n > m} computes the
#' \eqn{n - m} missing ones, and \eqn{n \le m} uses the first \eqn{n} without
#' computing anything. [summary.StatmodFit()] and [confint.StatmodFit()]
#' resolve the variance once, so the replicates are computed once.
#'
#' The coordinates are those of [vcov.StatmodFit()]: the coefficients of every
#' equation, with the aliased ones missing. A coefficient under a kinked
#' penalty ([modelterms7::lasso()], [modelterms7::scad()],
#' [modelterms7::mcp()]) that the fit holds at zero has a variance here,
#' which includes the mass of the replicates at zero, where the other
#' conventions report it as missing. With `readable = TRUE` the readable
#' quantities are computed for each replicate and the covariance is taken
#' afterwards.
#'
#' The variance estimates the sampling variance of the penalized estimator. It
#' does not contain the smoothing bias that [bayesian_vcov()] carries as
#' variance, so an interval around a penalized term built from it covers less
#' than the nominal level, as one built from [frequentist_vcov()] does. With
#' `refit = "coefficients"` the hyperparameters are held at the values of the
#' fit, and the matrix is a relative of [frequentist_vcov()]. With
#' `refit = "full"` the hyperparameters are selected again on every replicate,
#' and the variability of the hyperparameters enters; it costs about ten times
#' more. The relative Monte Carlo error of a standard error is about
#' \eqn{1/\sqrt{2n}}. A replicate that fails (an error, a mode error above
#' [mode_error_limit()], a coefficient that is not identified) is replaced;
#' replacing conditions the bootstrap on success, which understates the
#' variability when failures are frequent.
#'
#' The convention has no `information` argument: the matrix is not built from
#' the information. The refits use the `inner_information` of the fit.
#'
#' @param n A single whole number, the total number of replicates, or `NULL`.
#' @param resample A [resample_method] object, or `NULL`.
#' @param refit `"coefficients"`, `"full"` or `NULL`.
#'
#' @return A `BootstrapVcov` object, inheriting from [vcov_method], with the
#'   property `name` equal to `"bootstrap"` and the properties `n`, `resample`
#'   and `refit`.
#'
#' @seealso [add_bootstrap()], which computes and stores the replicates,
#'   [statmod_bootstrap()], which reads them, and [vcov_method].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' bootstrap_vcov()
#' fit <- add_bootstrap(fit, n = 40, resample = case_resample())
#' sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov())))
#' sqrt(diag(vcov(fit)))
#' @export
bootstrap_vcov <- function(n = NULL, resample = NULL, refit = NULL) {
  if (!is.null(n)) n <- check_boot_n(n)
  if (!is.null(resample)) resample <- as_resample_method(resample)
  if (!is.null(refit)) refit <- check_boot_refit(refit)
  BootstrapVcov(name = "bootstrap", n = n, resample = resample, refit = refit)
}

#' Check the Arguments of the Bootstrap
#'
#' @description
#' `check_boot_n()` returns `n` as an integer if it is one positive whole
#' number, and `check_boot_refit()` returns `refit` if it is `"coefficients"`
#' or `"full"`. Both signal an error otherwise.
#'
#' @param n The number of replicates.
#' @param refit The refit.
#'
#' @return The checked value.
#'
#' @keywords internal
check_boot_n <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n != round(n)) {
    stop("'n' must be a positive whole number.", call. = FALSE)
  }
  as.integer(n)
}

#' @rdname check_boot_n
#' @keywords internal
check_boot_refit <- function(refit) {
  if (!is.character(refit) || length(refit) != 1L ||
      !refit %in% c("coefficients", "full")) {
    stop("'refit' must be \"coefficients\" or \"full\".", call. = FALSE)
  }
  refit
}

S7::method(print, BootstrapVcov) <- function(x, ...) {
  cat("<vcov> bootstrap (",
      if (is.null(x@n)) "n as stored" else paste(x@n, "replicates"), ", ",
      if (is.null(x@resample)) "resample as stored" else
        resample_label(x@resample), ", ",
      if (is.null(x@refit)) "refit as stored" else paste0("refit ", x@refit),
      ")\n", sep = "")
  invisible(x)
}


#' @title Stored Replicates of a Bootstrap
#' @name StatmodBootstrap
#'
#' @description
#' One entry of the bootstrap collection of a fit: the coefficients of the
#' refitted replicates and the record of how they were obtained.
#' [statmod_bootstrap()] returns it, and [coef()] on it returns the matrix of
#' replicates.
#'
#' @details
#' The properties are
#' \describe{
#'   \item{`coefficients`}{a matrix with one row per replicate and one column
#'     per coordinate, in the order of `coef(fit, readable = FALSE)`;}
#'   \item{`hyper`}{a matrix of the hyperparameters of each replicate, with
#'     `refit = "full"`, and `NULL` otherwise;}
#'   \item{`n`, `resample`, `refit`}{the number of replicates, the
#'     [resample_method] and the refit;}
#'   \item{`random`}{how the random effects were treated: `"redrawn"` or
#'     `"none"`;}
#'   \item{`status`}{an integer for every draw of the stream up to the last
#'     replicate kept: 0 for a replicate that was kept, 1 for an error, 2
#'     for a mode error above [mode_error_limit()] and 3 for a coefficient
#'     that is not identified;}
#'   \item{`seed`}{a list with the state of `.Random.seed` at the start of the
#'     stream and after the last draw that was kept, from which an extension
#'     continues;}
#'   \item{`fingerprint`}{the coefficients and hyperparameters of the fit that
#'     the replicates belong to.}
#' }
#'
#' @param x,object A `StatmodBootstrap` object.
#' @param ... Unused.
#'
#' @return `coef()` returns the matrix of replicates. The class itself has no
#'   constructor to call.
#'
#' @seealso [add_bootstrap()], [statmod_bootstrap()], [bootstrap_vcov()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' fit <- add_bootstrap(fit, n = 20, resample = case_resample())
#' entry <- statmod_bootstrap(fit, resample = case_resample())
#' dim(coef(entry))
#' @keywords internal
#' @export
StatmodBootstrap <- S7::new_class("StatmodBootstrap",
  properties = list(
    coefficients = S7::class_any,
    hyper = S7::class_any,
    n = S7::class_integer,
    resample = resample_method,
    refit = S7::class_character,
    random = S7::class_character,
    status = S7::class_integer,
    seed = S7::class_list,
    fingerprint = S7::class_numeric
  )
)

S7::method(coef, StatmodBootstrap) <- function(object, ...) {
  object@coefficients
}

S7::method(print, StatmodBootstrap) <- function(x, ...) {
  fl <- boot_failures(x)
  cat("<bootstrap> ", resample_label(x@resample), ", refit ", x@refit, ", ",
      x@n, " replicates (", sum(fl), " replaced: ", fl[["error"]],
      " errors, ", fl[["mode_error"]], " mode errors, ", fl[["unidentified"]],
      " unidentified)\n", sep = "")
  invisible(x)
}

#' The Failures Recorded in a Bootstrap Entry
#'
#' @description
#' Counts, by cause, the draws that failed before the last replicate kept.
#'
#' @param x A `StatmodBootstrap` entry.
#' @param n Count only the draws up to the `n`-th replicate kept; `NULL` counts
#'   all of them.
#'
#' @return An integer vector named `error`, `mode_error` and `unidentified`.
#'
#' @keywords internal
boot_failures <- function(x, n = NULL) {
  st <- x@status
  if (!is.null(n) && n < x@n) {
    ok <- which(st == 0L)
    st <- st[seq_len(ok[[n]])]
  }
  c(error = sum(st == 1L), mode_error = sum(st == 2L),
    unidentified = sum(st == 3L))
}
