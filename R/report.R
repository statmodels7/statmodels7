#' @include statmod.R
#' @importFrom stats logLik coef residuals
NULL

#' Format a Duration in a Readable Unit
#'
#' @description
#' Converts a number of seconds to a string in microseconds, milliseconds,
#' seconds, minutes, hours or days, whichever keeps the number readable.
#'
#' @details
#' Fit times range from a few hundred microseconds to several hours, so the
#' unit is chosen from the magnitude of the value.
#'
#' The thresholds are the usual ones: a value reads in microseconds below a
#' millisecond, in milliseconds below a second, in seconds below a minute,
#' and then in minutes, hours and days. The unit is chosen after the value
#' is rounded to `digits` significant digits, so that a value just below a
#' threshold reads in the larger unit (59.99 seconds reads `"1 min"`).
#'
#' A duration of zero is reported as `"0 s"`. On a coarse clock a fast fit
#' can measure exactly zero seconds, and suppressing that value would make
#' the printed output depend on the timer resolution of the platform.
#'
#' @param seconds A number of seconds. Vectorized: a vector in gives a vector
#'   out, and each element is in its own unit.
#' @param digits How many significant digits to keep, passed to [signif()].
#'   Defaults to 3.
#'
#' @return A character vector of the length of `seconds`. Each element is a
#'   number and a unit separated by a space, as `"340 us"` or `"2.6 h"`. A
#'   non-finite value gives `"NA"` and a negative value gives `"negative"`.
#'
#' @seealso [print.StatmodFit()], which reports the elapsed time of a fit
#'   with this function.
#'
#' @examples
#' format_duration(c(3.4e-4, 0.25, 90, 7200))
#'
#' # A duration of zero is reported in seconds, a missing value as "NA".
#' format_duration(c(0, NA))
#'
#' @export
format_duration <- function(seconds, digits = 3L) {
  vapply(seconds, function(s) {
    if (!is.finite(s)) return("NA")
    if (s < 0) return("negative")
    if (s == 0) return("0 s")
    # the first unit in which the ROUNDED value stays below the next unit's
    # threshold, so that 59.99 s reads "1 min" and not "60 s"
    scale <- c(us = 1e-6, ms = 1e-3, s = 1, min = 60, h = 3600, d = 86400)
    limit <- c(1000, 1000, 60, 60, 24, Inf)
    for (k in seq_along(scale)) {
      v <- signif(s / scale[[k]], digits)
      if (v < limit[[k]]) return(paste0(v, " ", names(scale)[[k]]))
    }
  }, character(1))
}




#' @title Print a Fitted Model
#' @name print.StatmodFit
#' @description
#' Prints the call, the distribution and the number of rows of a fitted
#' model. For each distribution parameter it prints the formula and one line
#' per term, giving the number of coefficients (or of structural parameters)
#' and, where it differs from that number, the effective degrees of freedom
#' of the term. Where an outer criterion ran, it prints the value of the
#' criterion with the number of hyperparameter evaluations. It then prints
#' the elapsed time and, where the fit made several passes over blocks, the
#' number of passes. It prints no log-likelihood and no convergence verdict.
#'
#' [summary.StatmodFit()] is the long form, with coefficients, standard
#' errors and the hyperparameters.
#' @param x A [StatmodFit()].
#' @param ... Unused.
#' @return `x`, invisibly, as a print method should.
#' @seealso [summary.StatmodFit()] for the full report,
#'   [statmod_certificate()] for the reading of the point that the fit
#'   reached.
#' @keywords internal
print.StatmodFit <- function(x, ...) {
  spec <- x@spec
  cat("A statmod fit\n\n")
  cat("Call:  ", paste(deparse(x@call), collapse = "\n        "), "\n\n",
      sep = "")
  cat("Distribution: ", spec@distrib@distrib_name, "\n", sep = "")
  cat("Observations: ", spec@n_obs, sep = "")
  sw <- sum(spec@weights)
  if (!isTRUE(all.equal(sw, spec@n_obs))) {
    cat(sprintf("  (prior weights summing to %s, not %d)",
                format(signif(sw, 6)), spec@n_obs))
  }
  cat("\n\n")

  w <- if (is.null(x@edf) || !nrow(x@edf)) 0L else max(nchar(x@edf$term))
  for (p in spec@distrib@params) {
    rhs <- paste(deparse(spec@equations[[p]][[2L]]), collapse = " ")
    cat(sprintf("  %-10s ~ %s\n", p, rhs))
    if (!is.null(x@edf)) {
      rows <- x@edf[x@edf$parameter == p, , drop = FALSE]
      for (i in seq_len(nrow(rows))) {
        tm <- spec@terms[[p]][[rows$term[i]]]
        st <- !is.null(tm) &&
          identical(term_block_kind(tm), "structural")
        k <- if (st) length(modelterms7::term_params(tm)) else
          rows$coefficients[i]
        cat(sprintf("  %-10s   %s %3d %s", "",
                    formatC(rows$term[i], width = -w, flag = " "), k,
                    if (st) "param" else "coef"))
        if (is.finite(rows$edf[i]) && !isTRUE(all.equal(rows$edf[i], k))) {
          cat(sprintf(", edf %.2f", rows$edf[i]))
        }
        cat("\n")
      }
    }
  }

  cat("\n")
  # NO LIKELIHOOD HERE. There are two of them and they answer different
  # questions -- the conditional one the criteria are built on, and the
  # penalized one the inner fit minimizes -- so a single line either carried
  # a number whose meaning changed with the model or invited the two views
  # to be read against each other. `summary()` reports the conditional one
  # beside the effective degrees of freedom and the criteria, which is where
  # the pairing means something; this view says what the model IS.
  if (!is.null(x@methods$outer) && length(x@criterion) &&
      is.finite(x@criterion)) {
    cat(sprintf("%s %.6f over %d hyperparameter evaluation(s)\n",
                toupper(x@methods$outer@kind), x@criterion,
                if (is.null(x@history$outer)) 0L else nrow(x@history$outer)))
  }
  # no verdict on the search here either (Giovanni, 2026-09-16): its flag
  # was wrong four times in five where it fired, and summary() reports the
  # point itself through statmod_certificate()'s readings
  cat(sprintf("fitted in %s\n", format_duration(x@elapsed)))
  if (!is.null(x@history$blocks) && nrow(x@history$blocks) > 1L) {
    cat(sprintf("%d pass(es) over %d block(s)\n",
                max(x@history$blocks$pass),
                length(unique(x@history$blocks$block))))
  }
  invisible(x)
}
S7::method(print, StatmodFit) <- print.StatmodFit


#' The Model as a Function of Parameters and Data
#'
#' @description
#' Turn a fitted model back into a function. `loglik()`, `gradient()` and
#' `hessian()` evaluate the model's log-likelihood and its first two
#' derivatives at coefficients and data of the caller's choosing, so a fit
#' can be profiled, differenced or handed to a user-written optimizer.
#'
#' If only the fit is given, these functions evaluate the model at the fitted
#' coefficients and on the fitting data.
#'
#' @details
#' # The name
#'
#' The generic is `loglik`, not `logLik`. R's [stats::logLik()] returns the
#' log-likelihood of a fitted model and carries `df` and `nobs` with it, and
#' overloading it would give one name two behaviors. Both functions are
#' available, and `loglik(fit)` with no further arguments equals
#' `logLik(fit)`, which checks that the callable route and the fitting route
#' evaluate the same model.
#'
#' # Generics instead of closures
#'
#' A closure captures its environment, which here means the data. A fit
#' carrying one would hold the frame twice and would keep a stale copy after
#' the data changed.
#'
#' The generics rebuild from the specification that the fit keeps whole,
#' running each term's `build_info` against `data` by the same path as
#' [predict.StatmodFit()]. A factor's levels, a spline's knots and a basis
#' reparametrization are therefore reapplied to the new rows and not
#' estimated again.
#'
#' # Likelihood without the penalty
#'
#' All three functions evaluate the likelihood without the penalty.
#' The penalized objective at given coefficients is a different quantity,
#' and the hyperparameters are not arguments of these functions.
#'
#' @param object A [StatmodFit()].
#' @param par A named list of coefficient vectors, each as long as the
#'   design of the equation that it names. A distribution parameter that the
#'   list omits keeps its fitted coefficients, so a list may change one equation
#'   and leave the others at the fit. Defaults to the fitted coefficients. A
#'   `par` that is not a named list, a name that is not a distribution
#'   parameter, or a vector of the wrong length is an error. The own
#'   parameters of a structural term are not among the arguments.
#' @param data A data frame with the columns that the model names. Defaults to
#'   the data used to fit the model. New data must carry the response, since a
#'   likelihood needs one.
#' @param information `hessian()` only: the matrix whose negative is returned,
#'   an object that inherits from [distributions7::Information].
#'   [distributions7::observed_info()] (the default) gives the Hessian of the
#'   log-likelihood at the data, and [distributions7::expected_info()] gives
#'   minus the Fisher information, whose `method` selects how the expectation
#'   is obtained. [auto_info()] is not accepted. The default is
#'   [distributions7::observed_info()], whatever the family, because the
#'   result is a Hessian and not a variance; the variance objects resolve a
#'   `NULL` information by the family (see [information_objects]).
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return `loglik()` gives a single number. `gradient()` gives a named list
#'   of numeric vectors, one per distribution parameter. `hessian()` gives a
#'   symmetric `p x p` matrix over the stacked regression coefficients, where
#'   `p` is their total count. The own parameters of a structural term are not
#'   among them and have no entry in the gradient or in the Hessian.
#'
#' @seealso [logLik.StatmodFit()] for the log-likelihood of the fit with its
#'   degrees of freedom, [vcov.StatmodFit()] for the variance matrix,
#'   [predict.StatmodFit()] for the same reapplication on new data.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # At the fitted coefficients it equals the conditional log-likelihood.
#' all.equal(loglik(fit), as.numeric(logLik(fit)))
#'
#' # Far from the fit it is lower, and the gradient gives the direction of
#' # increase.
#' loglik(fit, par = list(mu = c(0, 0), sigma = 0))
#' gradient(fit, par = list(mu = c(0, 0), sigma = 0))
#'
#' # The gradient with respect to the mu coefficients vanishes at the fit.
#' # The sigma coordinate is estimated on the REML criterion, so its
#' # component is minus the number of mu coefficients (here -2).
#' gradient(fit)
#'
#' # The observed information is positive definite at the fit.
#' eigen(-hessian(fit), only.values = TRUE)$values
#'
#' @export
loglik <- S7::new_generic("loglik", "object",
  function(object, par = NULL, data = NULL, ...) S7::S7_dispatch())

#' @rdname loglik
#' @export
gradient <- S7::new_generic("gradient", "object",
  function(object, par = NULL, data = NULL, ...) S7::S7_dispatch())

#' @rdname loglik
#' @export
hessian <- S7::new_generic("hessian", "object",
  function(object, par = NULL, data = NULL,
           information = distributions7::observed_info(), ...)
    S7::S7_dispatch())


#' Rebuild a Specification Against New Data
#'
#' @description
#' Returns the fit's own specification when `data` is `NULL`, and a
#' specification with the fitted terms reapplied to the new rows otherwise.
#' The two cases are distinguished here only, so that every route that
#' accepts a `data` argument treats `NULL` alike.
#'
#' @param fit A [StatmodFit()].
#' @param data A data frame, or `NULL` for the data used to fit the model.
#' @param need_response `TRUE` where the response must be present, as for a
#'   log-likelihood; `FALSE` for a prediction, because new data routinely has
#'   no response column.
#'
#' @return A [StatmodSpec()]: `fit@spec` itself when `data` is `NULL`, and
#'   [statmod_respec()]'s result otherwise.
#'
#' @keywords internal
spec_at <- function(fit, data, need_response = TRUE) {
  if (is.null(data)) return(fit@spec)
  statmod_respec(fit@spec, data, need_response = need_response)
}

#' Resolve a Parameter Structure
#'
#' @description
#' Returns the fitted coefficients when `par` is `NULL`, and otherwise checks
#' a supplied structure against the design before returning it.
#'
#' @details
#' Three things are checked: that `par` is a named list, that every name of
#' `par` is a distribution parameter, and that the vector given for each is
#' as long as that equation's design is wide. Each failure signals an error.
#' The message for a wrong name and the message for a wrong length both give
#' the name of the parameter, so that a short vector is not recycled into a
#' different model. A distribution
#' parameter that `par` omits keeps its fitted coefficients.
#'
#' @param fit A [StatmodFit()], read for its fitted coefficients.
#' @param par A named list of coefficient vectors, or `NULL`.
#' @param design The design against which the lengths are validated.
#'
#' @return A named list of coefficient vectors, one per distribution
#'   parameter in the family's order.
#'
#' @seealso [loglik()], whose three methods all pass through this.
#'
#' @keywords internal
par_at <- function(fit, par, design) {
  if (is.null(par)) return(fit@coefficients)
  params <- fit@spec@distrib@params
  if (!is.list(par) || is.null(names(par)) || anyNA(names(par)) ||
      !all(nzchar(names(par)))) {
    stop(sprintf(paste0("'par' must be a named list of coefficient vectors, ",
                        "one per distribution\n  parameter that it changes ",
                        "(%s)."), paste(params, collapse = ", ")),
         call. = FALSE)
  }
  bad <- setdiff(names(par), params)
  if (length(bad)) {
    stop(sprintf(paste0("'par' names '%s', which is not a parameter.\n",
                        "  They are: %s."),
                 bad[1L], paste(params, collapse = ", ")), call. = FALSE)
  }
  out <- fit@coefficients
  for (p in names(par)) {
    if (length(par[[p]]) != design[[p]]$npar) {
      stop(sprintf("'par$%s' has length %d but '%s' has %d coefficients.",
                   p, length(par[[p]]), p, design[[p]]$npar), call. = FALSE)
    }
    out[[p]] <- as.numeric(par[[p]])
  }
  out
}

#' @rdname loglik
#' @name loglik.StatmodFit
#' @keywords internal
S7::method(loglik, StatmodFit) <- function(object, par = NULL, data = NULL,
                                           ...) {
  check_dots("loglik", ...)
  spec <- spec_at(object, data)
  design <- statmod_design(spec)
  statmod_loglik_at(spec, par_at(object, par, design), design)
}

#' @rdname loglik
#' @name gradient.StatmodFit
#' @keywords internal
S7::method(gradient, StatmodFit) <- function(object, par = NULL, data = NULL,
                                             ...) {
  check_dots("gradient", ...)
  spec <- spec_at(object, data)
  design <- statmod_design(spec)
  statmod_score_at(spec, par_at(object, par, design), design)
}

#' @rdname loglik
#' @name hessian.StatmodFit
#' @keywords internal
S7::method(hessian, StatmodFit) <- function(object, par = NULL, data = NULL,
                                            information = distributions7::observed_info(),
                                            ...) {
  check_dots("hessian", ...)
  check_information(information, arg = "information", fn = "hessian",
                    auto = FALSE)
  spec <- spec_at(object, data)
  design <- statmod_design(spec)
  # the information is minus the Hessian, and what is asked for here is the
  # Hessian of the log-likelihood itself
  -statmod_information_at(spec, par_at(object, par, design), design,
                          info_expected(information), info_emethod(information))
}


#' @title The Log-Likelihood of a Fit
#' @name logLik.StatmodFit
#' @description
#' Returns the value that R's convention expects, with the degrees of freedom
#' and the number of observations as attributes. [loglik()] is a different
#' function: it evaluates the model at parameters and data chosen by the
#' caller.
#' @details
#' # Which likelihood, and why it matters
#'
#' The default is the conditional one: the log-density at the fitted
#' coefficients, including the penalized coefficients, paired with the
#' effective degrees of freedom \eqn{\mathrm{tr}[(H+S)^{-1}H]}. A criterion
#' built on that pair is the conditional AIC, which measures how well the
#' model describes the groups, curves and states actually observed.
#'
#' A mixed-model package reports the marginal likelihood instead: the
#' random effects are integrated out and the count is the number of estimated
#' parameters, including the variance components. That likelihood concerns the
#' population, and the AIC built on it is a different number, not comparable
#' with the conditional one.
#'
#' Each likelihood keeps its own count. A marginal likelihood against an
#' effective count, or a conditional one against a parameter count, mixes the
#' two conventions (Vaida and Blanchard, 2005).
#'
#' # When the marginal one is available
#'
#' `loglik_type = "marginal"` returns the value of the outer criterion at the fit.
#' Its degrees of freedom count the parameters that the marginal model
#' estimates: the coefficients that are not covered by any penalty, the
#' dimension of the null space of each penalty without a kink (the straight
#' line of a smooth), the hyperparameters that the outer search or a path
#' estimated, and the free parameters of a structural term. A held
#' hyperparameter is not
#' counted, and a hyperparameter shared by several terms is counted once. The
#' other coefficients of a penalized term are integrated out and are not
#' counted.
#'
#' The value exists only where the outer criterion is [ml()] or [reml()].
#' Where no outer criterion ran, or where it is [aic()] or [bic()], there is
#' no marginal likelihood, and requesting one signals an error.
#' @param object A [StatmodFit()].
#' @param loglik_type `"conditional"` (the default) or `"marginal"`.
#' @param ... Unused.
#' @return A `logLik` object: a single number with attributes `df`, the
#'   degrees of freedom, and `nobs`, the sum of the prior weights (the number
#'   of rows where the weights are all 1). [stats::AIC()] and [stats::BIC()]
#'   read it, so `AIC(fit)` is the conditional AIC.
#' @references
#' Vaida, F. and Blanchard, S. (2005). Conditional Akaike information for
#' mixed-effects models. *Biometrika*, 92(2), 351--370.
#' @seealso [loglik()]
#' @keywords internal
logLik.StatmodFit <- function(object,
                              loglik_type = c("conditional", "marginal"), ...) {
  check_dots("logLik", ...)
  loglik_type <- match.arg(loglik_type)
  if (identical(loglik_type, "marginal")) {
    outer <- object@methods$outer
    if (is.null(outer) || !outer@kind %in% c("ml", "reml")) {
      stop(paste0(
        "The marginal log-likelihood is the value that ml() or reml() ",
        "optimizes,\n  and ",
        if (is.null(outer)) "this fit ran no outer criterion" else
          sprintf("this fit used %s()", outer@kind),
        ". There is no marginal likelihood to\n  report; logLik(loglik_type = ",
        "\"conditional\") is available."), call. = FALSE)
    }
    if (!length(object@criterion) || !is.finite(object@criterion)) {
      stop("The outer criterion was not recorded for this fit.",
           call. = FALSE)
    }
    # The parameters of the MARGINAL model, which is a different count from
    # the conditional one and not an effective anything: the coefficients it
    # did not integrate away, one per hyperparameter it estimated, and a
    # structural term's own parameters, which are estimated rather than
    # integrated. A penalized coefficient is a random effect under another
    # name and is integrated out, so it does not appear.
    design <- statmod_design(object@spec)
    lab <- coef_labels(object@spec, design)
    # the hyperparameters ESTIMATED: one per row of the outer index (a held
    # one has none, a shared group has one) and the ones a path chose
    nh <- nrow(outer_hyper_index(object@spec,
                                 statmod_blocks(object@spec, design))) +
      length(object@methods$sparse_hyper)
    ns <- sum(vapply(statmod_structural_par(object@spec, design),
                     function(u) length(u$parameter) - length(u$held),
                     numeric(1)))
    # coef_labels() marks a whole penalized term, so the directions that its
    # penalty leaves free (a smooth's straight line) are added here: they are
    # estimated, not integrated, as lme4 counts its fixed effects
    nn <- sum(vapply(statmod_penalized(object@spec, design), function(u) {
      if (isTRUE(u$structural) ||
          isTRUE(tryCatch(penalty_has_kink(u$penalty),
                          error = function(e) TRUE))) return(0)
      N <- tryCatch(as.matrix(penalties7::penalty_null_basis(u$penalty)),
                    error = function(e) NULL)
      if (is.null(N)) 0 else ncol(N)
    }, numeric(1)))
    return(structure(as.numeric(object@criterion),
                     df = sum(!lab$penalized) + nn + nh + ns,
                     nobs = weight_total(object@spec), class = "logLik"))
  }
  df <- if (is.null(object@edf)) {
    sum(lengths(object@coefficients))
  } else {
    # a term whose count could not be obtained falls back to its number of
    # columns, which is the upper bound: dropping it would understate the
    # dimension of the model and flatter every criterion built on it
    #
    # A STRUCTURAL TERM HAS NO COLUMNS, so that bound is zero for it and the
    # rule above says the opposite of what it means: measured, a model whose
    # filter carries eighteen free parameters reported a df of 2. Its own
    # upper bound is the parameters it estimates, which is what the count
    # would have been had the smoother been readable.
    e <- object@edf
    bound <- e$coefficients
    if (anyNA(e$edf)) {
      sp <- tryCatch(statmod_structural_par(object@spec,
                                            statmod_design(object@spec)),
                     error = function(err) list())
      for (nm in names(sp)) {
        i <- which(e$term == nm)
        if (length(i))
          bound[i] <- length(sp[[nm]]$parameter) - length(sp[[nm]]$held)
      }
    }
    sum(ifelse(is.na(e$edf), bound, e$edf))
  }
  structure(object@loglik, df = df, nobs = weight_total(object@spec),
            class = "logLik")
}
S7::method(logLik, StatmodFit) <- logLik.StatmodFit

#' @title The Coefficients of a Fitted Model
#'
#' @description
#' Returns one named vector per distribution parameter: the quantities in
#' which the model is written, or the coordinates on which it was estimated.
#'
#' @details
#' # What `readable` moves, and what it does not
#'
#' Most coefficients are the same under either reading: the coefficient of a
#' linear predictor is reported as it is estimated. Three kinds of term
#' report quantities that differ from the coordinates on which they are
#' carried, and those are the parameters that this argument moves.
#'
#' A break-point term is fitted through a working block, and its
#' quantities are read from that block. A continuous term holds each quantity
#' as a coefficient, and the two readings coincide. A discontinuous term
#' carries \eqn{g} and \eqn{\delta}; where its position is not held, it
#' reports \eqn{\psi = -g/\delta}, so that at `readable = FALSE` the vector
#' holds \eqn{g}, which is no quantity of the model. Where the term holds its
#' position, as a sharp `jump()` or `jseg()` does after the fit has settled
#' it, the reported position is the held one.
#'
#' A score-driven term carries its level and its loadings
#' \eqn{\kappa_j} through their own links, the log for a loading, and its
#' persistence on partial autocorrelations, because the stationary region is
#' not a box. Each partial autocorrelation is carried on the unconstrained
#' scale of its link, the inverse hyperbolic tangent by default. The quantity
#' \eqn{\xi_j} of the literature is the autoregressive coefficient that the
#' whole chart produces. With one persistence, \eqn{\xi_1} is the partial
#' autocorrelation. With two, \eqn{\xi_1 = \mathrm{pacf}_1 (1 -
#' \mathrm{pacf}_2)} and \eqn{\xi_2 = \mathrm{pacf}_2}.
#'
#' A regime term is estimated on the additive log-ratios of each row of
#' its transition matrix, and it reports the transition probabilities
#' \eqn{p_{ij}} in their place. Its levels and gaps are reported through
#' their own links.
#'
#' Where a term declares no quantities of its own, the coordinates stand. So
#' do the coefficients of a parameter developed over covariates: a
#' development is a vector with no single value to report, so the term
#' declares nothing for it.
#'
#' # A structural term is present under both readings
#'
#' A structural term contributes no design columns, so its parameters are in
#' no block of the design. Both readings append them to the vector of the
#' equation, so that a model whose whole predictor is a score-driven filter
#' returns its filter parameters.
#'
#' # The coordinates of the fit
#'
#' `readable = FALSE` gives the coordinates on which the fit was estimated,
#' with a structural term's part on the unconstrained scale that its charts
#' define. For a model without a structural term, this vector has the order of
#' the rows of `vcov(object, readable = FALSE)`, and it can be passed as the
#' `par` argument of [loglik()] or to an optimizer. A structural level that
#' is held appears in the vector at its held value and has no row in that
#' matrix.
#'
#' Hyperparameters are not coefficients and are not here. [hyper()] reports
#' them.
#'
#' @param object A [StatmodFit()].
#' @param readable `TRUE`, the default, reports the quantities in which the
#'   model is written. `FALSE` reports the coordinates on which it was
#'   estimated.
#' @param ... Unused.
#'
#' @return A named list with one entry per distribution parameter, in the
#'   family's order, each a named numeric vector. The names are the
#'   coefficient labels, composed from each term's own label, so two terms of
#'   one kind in one formula stay apart.
#'
#' @seealso [hyper()] for the hyperparameters, [vcov.StatmodFit()] for the
#'   variance matrix, whose rows with `readable = FALSE` are the coordinates
#'   of `coef(object, readable = FALSE)` that are not held,
#'   [summary.StatmodFit()] for the two together,
#'   [modelterms7::term_readable()] for what a term declares.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(200, 0, 10)))
#' d$y <- 0.3 * d$x + 1.5 * pmax(d$x - 6, 0) + rnorm(200, 0, 0.4)
#' fit <- statmod(y ~ modelterms7::seg(x, psi = 4),
#'                distributions7::gaussian1_distrib(), d)
#'
#' # The break-point is a quantity of the model, near the true 6.
#' coef(fit)$mu
#'
#' # The coordinates on which the fit ran are the same numbers here, because
#' # seg() is continuous, and vcov() has one row for each of them.
#' coef(fit, readable = FALSE)$mu
#' rownames(vcov(fit, readable = FALSE))
#'
#' @keywords internal
coef.StatmodFit <- function(object, readable = TRUE, ...) {
  check_dots("coef", ...)
  spec <- object@spec
  params <- spec@distrib@params
  design <- statmod_design(spec)
  stats::setNames(lapply(params, function(p) {
    v <- stats::setNames(object@coefficients[[p]], design[[p]]$coef_names)
    # An aliased coordinate has no estimate: the design carries the same
    # information in another column and the pivot dropped this one, so
    # whatever the solve leaves there is not an estimate. `lm()` and `glm()`
    # report NA and so does this. The STORED value is left where it was --
    # the pivot drops the coordinate from the INCREMENT and not from the
    # parameter -- which keeps every predictor and every assembled quantity
    # reading the vector as before, and is the sharper reason for reporting
    # NA: measured, the same fit started with that column at 5 leaves 5
    # there and reaches the identical log-likelihood, so the number is one
    # point of a flat ridge rather than a quantity.
    if (length(object@aliased)) {
      na <- paste(p, names(v), sep = ":") %in% object@aliased
      if (any(na)) v[na] <- NA_real_
    }
    if (isTRUE(readable)) v <- coef_readable(spec, design, object, p, v)
    c(v, coef_structural(spec, object, p, readable))
  }), params)
}

#' The Quantities Reported by a Term in Place of Its Coordinates
#'
#' @description
#' Returns the coefficient vector of one equation, in which each term's
#' declared quantities replace the coordinates from which they are computed.
#'
#' @details
#' A term declares its quantities through [modelterms7::term_readable()],
#' which returns the quantities and their Jacobian with respect to the
#' coefficients from which they are computed. The columns that the Jacobian
#' touches are exactly the coordinates to replace. A coordinate that enters
#' no declared quantity stays in place.
#'
#' That rule keeps a developed parameter intact. A development is a vector of
#' coefficients over covariates with no single value to report, so the term
#' declares nothing for it and its coefficients stay.
#'
#' The names are composed from the term's own label, as the names of its
#' coefficients are, so that two terms of one kind in one formula stay apart.
#'
#' @param spec The fitted specification.
#' @param design The design.
#' @param fit The [StatmodFit()].
#' @param p The distribution parameter naming the equation, a string.
#' @param v The named coefficient vector of that equation.
#'
#' @return `v`, in which each term's declared quantities replace the
#'   coordinates from which they are computed. The vector `v` is returned
#'   unchanged where no term of the equation declares any.
#'
#' @seealso [coef.StatmodFit()], the caller,
#'   [modelterms7::term_readable()] for what a term declares.
#'
#' @keywords internal
coef_readable <- function(spec, design, fit, p, v) {
  keep <- rep(TRUE, length(v))
  ins <- list()
  for (nm in names(spec@terms[[p]])) {
    idx <- design[[p]]$blocks[[nm]]
    if (!length(idx)) next
    term <- spec@terms[[p]][[nm]]
    rd <- tryCatch(modelterms7::term_readable(term,
                                              fit@coefficients[[p]][idx]),
                   error = function(e) NULL)
    if (is.null(rd) || !length(rd$name)) next
    sup <- which(apply(rd$jacobian != 0, 2L, any))
    if (!length(sup)) next
    keep[idx[sup]] <- FALSE
    ins[[as.character(idx[sup][1L])]] <- stats::setNames(
      as.numeric(rd$value), paste(term@label, rd$name, sep = "."))
  }
  if (!length(ins)) return(v)
  out <- list()
  for (i in seq_along(v)) {
    k <- as.character(i)
    if (!is.null(ins[[k]])) out[[length(out) + 1L]] <- ins[[k]]
    if (keep[[i]]) out[[length(out) + 1L]] <- v[i]
  }
  unlist(out)
}

#' What a Structural Term Contributes to the Coefficients
#'
#' @description
#' Returns the own parameters of a structural term of one equation, as
#' quantities or as the coordinates on which they were estimated.
#'
#' @details
#' A structural term rewrites the likelihood and adds no columns to a design,
#' so its parameters sit in no block. [coef.StatmodFit()] appends them to the
#' vector of the equation under both readings, named from the term's label as
#' every other coefficient of that term is.
#'
#' A held parameter is one that an intercept in the same equation
#' carries, so it is not estimated. It is reported under either reading, at
#' the value at which it is held. Leaving it out would make the vector
#' shorter than the term's own parameter count, and the readable vector would
#' no longer correspond to the rows of `vcov(readable = TRUE)`, where the
#' held level has a row of missing values.
#'
#' @param spec The fitted specification.
#' @param fit The [StatmodFit()].
#' @param p The distribution parameter naming the equation, a string.
#' @param readable `TRUE` for the quantities that the term declares, `FALSE` for
#'   the coordinates on the unconstrained scale.
#'
#' @return A named numeric vector of that term's own parameters.
#'   The result is `numeric(0)` where the equation carries no structural term.
#'
#' @seealso [coef.StatmodFit()], the caller,
#'   [statmod_latent()] for a latent-state term's smoothed states.
#'
#' @keywords internal
coef_structural <- function(spec, fit, p, readable) {
  out <- numeric(0)
  sp <- fit@structural
  if (!length(sp)) return(out)
  for (nm in names(sp)) {
    eq <- sp[[nm]]$equation
    term <- spec@terms[[eq]][[nm]]
    if (!identical(eq, p) && !(is_joint_filter(term) && p %in% term@params)) {
      next
    }
    z <- sp[[nm]]$unconstrained
    if (isTRUE(readable)) {
      rd <- tryCatch(modelterms7::term_readable(term, z),
                     error = function(e) NULL)
      w <- if (is.null(rd) || !length(rd$name)) unlist(sp[[nm]]$parameter) else
        stats::setNames(as.numeric(rd$value), rd$name)
    } else {
      w <- unlist(z)
    }
    ds <- structural_display(term, eq, names(w))
    w <- w[ds$param == p]
    names(w) <- ds$name[ds$param == p]
    out <- c(out, w)
  }
  out
}
S7::method(coef, StatmodFit) <- coef.StatmodFit

#' @title The Residuals of a Fitted Model
#' @name residuals.StatmodFit
#' @description
#' Returns one residual per observation, which compares the observation with
#' the whole distribution that the model assigns to it and not with any single
#' parameter of that distribution.
#'
#' @details
#' # One residual per observation
#'
#' A residual measures whether an observation is consistent with the law that
#' was given to its row, and that law is one object that carries every
#' parameter at once. There is therefore one residual per observation,
#' however many parameters that the model develops over covariates.
#'
#' A per-parameter quantity is a different object: the contribution
#' \eqn{\partial \ell_i / \partial \eta_{ip}} identifies which equation an
#' observation strains, and the partial residuals of one equation are the
#' quantities against which a term's effect is drawn.
#'
#' # The quantile residual, and why it is the default
#'
#' \deqn{r_i = \Phi^{-1}(F(y_i; \hat\theta_i))}
#'
#' Under a correct model \eqn{F(y_i; \theta_i)} is uniform where the
#' distribution function is continuous, so \eqn{r_i} is standard normal, for
#' any such family and whichever of its parameters are modeled. The residual
#' does not single out any parameter, and its reference distribution needs no
#' asymptotics (Dunn and Smyth, 1996). Where the distribution function
#' jumps, the residual is randomized as described next.
#'
#' # Where the distribution function jumps
#'
#' The construction is randomized: \eqn{u_i} is drawn uniformly on
#' \eqn{(F(y_i^-), F(y_i))} and the residual is \eqn{\Phi^{-1}(u_i)}. The
#' residual is again exactly standard normal under a correct model, but it is
#' random, so that two calls give two answers.
#'
#' Where \eqn{F(y_i)} exceeds one half, the residual is read from the survival
#' function, \eqn{\Phi^{-1}(u_i) = -\Phi^{-1}(1 - u_i)}, and \eqn{1 - u_i} is
#' computed as \eqn{1 - F} directly, so a count far in the upper tail, where
#' \eqn{F} rounds to one, has a finite residual. A family whose
#' distribution function does not take `lower.tail` keeps the lower-tail
#' reading.
#'
#' The randomization applies to every discrete family and, for a mixed family
#' such as the zero-adjusted wrapper of a continuous parent, only at the
#' declared atoms. `seed` makes a call reproducible without disturbing the
#' stream of the caller. If `seed` is `NULL`, the ambient state is used and
#' the seed is not set.
#'
#' # Pearson and response residuals
#'
#' The Pearson residual is \eqn{(y_i - \mathbb{E}[Y_i]) / \mathrm{sd}(Y_i)},
#' and the response residual is its numerator. Both are defined against the
#' mean, and for a skewed family the Pearson residual is not standard normal
#' even where the model is right, so its quantile-quantile plot can mislead.
#' They are included for comparison with the residuals that other modeling
#' packages report.
#'
#' @param object A [StatmodFit()].
#' @param residual_type `"quantile"` (the default), `"pearson"` or
#'   `"response"`.
#'   The argument is matched with [match.arg()].
#' @param seed An integer used to seed the randomization, or `NULL` for the
#'   ambient state. Read only where the distribution function jumps, so it
#'   changes nothing for a continuous family.
#' @param ... Unused.
#'
#' @return A numeric vector with one entry per row of the data. The residuals
#'   are standard normal under a correct model for `"quantile"`;
#'   approximately standardized for `"pearson"`; on the response's own scale
#'   for `"response"`.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' *Journal of Computational and Graphical Statistics* 5(3), 236--244.
#'
#' @seealso [fitted.StatmodFit()] and [predict.StatmodFit()] for the fitted
#'   parameters against which these are read.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(200, -2, 2))
#' d$y <- 1 + 0.8 * d$x + rnorm(200, 0, 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
#'
#' # Under a correct model the quantile residuals are standard normal.
#' r <- residuals(fit)
#' c(mean = mean(r), sd = stats::sd(r))
#' stats::shapiro.test(r)$p.value
#'
#' # For a discrete family the same residual is randomized, so two calls
#' # differ unless a seed is given.
#' dp <- data.frame(x = runif(200, -1, 1))
#' dp$y <- rpois(200, exp(1 + dp$x))
#' fp <- statmod(y ~ x, distributions7::poisson_distrib(), dp)
#' identical(residuals(fp, seed = 1), residuals(fp, seed = 1))
#' identical(residuals(fp, seed = 1), residuals(fp, seed = 2))
#'
#' @keywords internal
residuals.StatmodFit <- function(object,
                                 residual_type = c("quantile", "pearson",
                                                   "response"),
                                 seed = NULL, ...) {
  check_dots("residuals", ...)
  residual_type <- match.arg(residual_type)
  spec <- object@spec
  d <- spec@distrib
  y <- spec@response
  th <- object@fitted
  if (!length(th)) {
    stop("The fit carries no fitted parameters against which to compare the response.",
         call. = FALSE)
  }
  if (identical(residual_type, "response") ||
      identical(residual_type, "pearson")) {
    m <- tryCatch(mean(d, th), error = function(e) {
      warning(sprintf(paste0(
        "The mean of '%s' could not be computed at the fitted parameters, ",
        "so the %s residuals are NA: %s It may not exist there."),
        d@distrib_name, residual_type, conditionMessage(e)), call. = FALSE)
      NA_real_
    })
    r <- as.numeric(y) - as.numeric(m)
    if (all(is.na(r))) return(r)
    if (identical(residual_type, "response")) return(r)
    # A family with no closed-form variance reaches it by quadrature, which
    # does not converge where the variance does not exist (a log-logistic of
    # shape below 2). That used to end in an error listing every observation;
    # the residual there has no value, and it is returned as NA.
    s <- tryCatch(as.numeric(distributions7::std_dev(d, th)),
                  error = function(e) {
                    warning(sprintf(paste0(
                      "The standard deviation of '%s' could not be computed ",
                      "at the fitted parameters, so the Pearson residuals are ",
                      "NA: %s It may not exist there."),
                      d@distrib_name, conditionMessage(e)), call. = FALSE)
                    rep(NA_real_, length(r))
                  })
    return(r / s)
  }
  fy <- as.numeric(distributions7::distrib_cdf(d, y, th))
  # WHERE THE DISTRIBUTION FUNCTION JUMPS, and by how much. A discrete family
  # jumps at every observation and its mass IS what distrib_pdf returns
  # there; a mixed one jumps at its declared atoms alone, where the same call
  # returns the atom's probability. Asking the density rather than
  # differencing the distribution function at the previous support point is
  # what keeps this right for a family whose support is not the integers.
  jump <- rep(FALSE, length(fy))
  if (S7::S7_inherits(d, distributions7::discrete_distrib)) {
    jump[] <- TRUE
  } else {
    at <- tryCatch(distributions7::distrib_atoms(d, th),
                   error = function(e) list(y = numeric(0)))
    if (length(at$y)) jump <- as.numeric(y) %in% at$y
  }
  u <- fy
  if (any(jump)) {
    py <- as.numeric(distributions7::distrib_pdf(d, y, th))
    lo <- pmax(fy - py, 0)
    if (!is.null(seed)) {
      # the caller's stream is left where it was: a residual is not a reason
      # to move it
      old <- if (exists(".Random.seed", envir = globalenv())) {
        get(".Random.seed", envir = globalenv())
      } else NULL
      set.seed(seed)
      on.exit({
        if (is.null(old)) {
          suppressWarnings(rm(".Random.seed", envir = globalenv()))
        } else {
          assign(".Random.seed", old, envir = globalenv())
        }
      }, add = TRUE)
    }
    v <- stats::runif(sum(jump))
    u[jump] <- lo[jump] + v * (fy[jump] - lo[jump])
  }
  r <- stats::qnorm(pmin(pmax(u, 0), 1))
  # IN THE UPPER TAIL THE RESIDUAL IS READ FROM THE SURVIVAL FUNCTION.
  # There F(y) rounds to 1 long before the residual is large: a Poisson
  # count of 67 at a mean of 15 has F(y - 1) = 1 in double precision and
  # came back Inf, where its residual is 9.9. A family whose distribution
  # function takes `lower.tail` returns S = 1 - F with its own digits; one
  # that ignores the argument returns F again, and the check S + F = 1 is
  # what tells the two apart, so such a family keeps the lower-tail reading.
  up <- which(fy > 0.5)
  if (length(up)) {
    sv <- tryCatch(as.numeric(distributions7::distrib_cdf(d, y, th,
                                                        lower.tail = FALSE)),
                   error = function(e) NULL)
    if (!is.null(sv) && length(sv) == length(fy)) {
      ok <- up[is.finite(sv[up]) & abs(sv[up] + fy[up] - 1) <= 1e-8]
      if (length(ok)) {
        w <- sv
        jj <- ok[jump[ok]]
        if (length(jj)) {
          vj <- v[match(jj, which(jump))]
          w[jj] <- sv[jj] + (1 - vj) * pmin(py[jj], fy[jj])
        }
        r[ok] <- stats::qnorm(pmin(pmax(w[ok], 0), 1), lower.tail = FALSE)
      }
    }
  }
  r
}
S7::method(residuals, StatmodFit) <- residuals.StatmodFit


#' @title The Hyperparameters of a Fitted Model
#'
#' @description
#' Reports every hyperparameter of every penalty that the model carries: one row
#' each, with the value, whether it was held, and what put it there. A
#' smoothing parameter, a prior scale, a lasso's \eqn{\lambda} and an elastic
#' net's \eqn{\alpha} all appear here.
#'
#' @details
#' # A hyperparameter is not a coefficient
#'
#' It governs the coefficients under it instead of sitting beside them, and
#' the two are estimated by different routes and reported with different
#' qualifications. [coef.StatmodFit()] reports the coefficients, and this
#' function reports the hyperparameters.
#'
#' # The two scales
#'
#' `scale = "parameter"`, the default, is the scale on which the penalty is
#' written: a smoothing parameter is a positive number, and the `sigma` of a
#' gaussian prior is a scale.
#'
#' `scale = "link"` is the free scale on which the outer search runs,
#' obtained through each hyperparameter's own link, and it is the scale on
#' which the searches of two fits are compared. Where a hyperparameter
#' carries no link, the two scales coincide.
#'
#' # The columns `held` and `source`
#'
#' `held` is a logical that records whether a term fixed the value, directly
#' or through the label that the hyperparameter shares. `source` records what
#' produced the value:
#'
#' - `"fixed"` for a hyperparameter that a term holds, as in
#'   `s(x, bspline_smooth(), hyper = c(lambda = 2))`.
#' - the name of the outer criterion (`"reml"`, `"ml"`, `"aic"` or `"bic"`)
#'   for a hyperparameter that the outer search chose.
#' - the name of the sparse criterion (`"aic"`, `"bic"` or `"cv"`) for a
#'   hyperparameter of a kinked penalty chosen along a path over its own
#'   values.
#' - `""` where no criterion ran, so that the value stayed at the starting
#'   value of its term.
#'
#' The distinction matters for inference. A value chosen by [aic()], [bic()]
#' or [cv()], whether by the outer search or along a path, is the argument of
#' a minimum of a prediction-error criterion and not the root of the
#' derivative of a likelihood, so the value does not yield a standard error.
#' A value that [reml()] or [ml()] estimated carries a standard error, computed
#' from the curvature of that criterion, and [summary.StatmodFit()] prints it.
#'
#' # Why two rows may carry the same number
#'
#' Terms whose [modelterms7::term_ids()] give the same label estimate one
#' value, and that value is written under each of their own keys. A shared
#' hyperparameter therefore appears once per term, with the same estimate and
#' the same `source`, and the `id` column marks the agreement as part of the
#' model and not as a coincidence. The penalties are not merged, so each term
#' still has a row of its own and the effective degrees of freedom are still
#' counted term by term.
#'
#' @param fit A [StatmodFit()].
#' @param scale `"parameter"` (the default) or `"link"`. The argument is
#'   matched with [match.arg()].
#'
#' @return A data frame with one row per hyperparameter and the columns:
#'   \describe{
#'     \item{`parameter`}{the distribution parameter whose equation contains
#'       the penalized term. For a covariance class that spans several
#'       equations, the names of all of them, separated by commas.}
#'     \item{`term`}{the term's key, its call as written.}
#'     \item{`name`}{the hyperparameter's own name, as the penalty names it.}
#'     \item{`estimate`}{its value, on the requested scale.}
#'     \item{`held`}{a logical: whether the term fixed it.}
#'     \item{`source`}{what produced the value, as above.}
#'     \item{`id`}{the sharing label, `NA` where the hyperparameter is not
#'       shared.}
#'   }
#'   The result is a data frame of no rows, with those columns, where the
#'   model carries no penalty at all.
#'
#' @seealso [coef.StatmodFit()] for the coefficients,
#'   [summary.StatmodFit()] for the two printed together with standard
#'   errors, [statmod_held()] for which are held.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(80, 0, 1))
#' d$y <- sin(3 * d$x) + rnorm(80, 0, 0.3)
#'
#' # The hyperparameter is estimated by REML, which the source column records.
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 6)),
#'                distributions7::gaussian1_distrib(), d)
#' hyper(fit)
#'
#' # The same value is shown on the scale on which the outer search ran.
#' hyper(fit, scale = "link")
#'
#' # Here the term holds the value instead, and it is reported as fixed.
#' held <- statmod(y ~ s(x, bspline_smooth(k = 6), hyper = c(lambda = 2)),
#'                 distributions7::gaussian1_distrib(), d)
#' hyper(held)[, c("name", "estimate", "held", "source")]
#'
#' @export
hyper <- function(fit, scale = c("parameter", "link")) {
  scale <- match.arg(scale)
  spec <- fit@spec
  held <- tryCatch(statmod_held(spec, statmod_design(spec)),
                   error = function(e) character(0))
  outer_kind <- if (is.null(fit@methods$outer)) NA_character_ else
    fit@methods$outer@kind
  spc <- fit@methods$sparse_criterion
  spc_kind <- if (is.null(spc)) NA_character_ else spc@kind
  rows <- list()
  units <- statmod_penalty_keys(spec)
  # a value one member of a sharing group holds is held for the group, so a
  # member carrying only the label is held too and says so
  held_ids <- names(tryCatch(statmod_held_ids(units), error = function(e) list()))
  for (u in units) {
    th <- fit@hyper[[u$param]][[u$key]]
    if (is.null(th) || !length(th)) next
    kink <- isTRUE(tryCatch(penalty_has_kink(u$penalty),
                            error = function(e) FALSE))
    for (h in names(th)) {
      v <- as.numeric(th[[h]])
      lk <- u$penalty@link_params[[h]]
      if (identical(scale, "link") && !is.null(lk)) {
        v <- linkfunctions7::linkfun(lk, v)
      }
      # HELD is the term's answer and nobody else's, and it is asked of the
      # same enumeration the outer index asks: a hyperparameter a term fixed
      # is fixed whatever criterion ran beside it.
      hk <- paste(u$param, u$key, h, sep = "\r")
      lab <- if (h %in% names(u$ids)) u$ids[[h]] else NA_character_
      is_held <- hk %in% held || h %in% names(u$fixed) ||
        (!is.na(lab) && lab %in% held_ids)
      src <- if (is_held) "fixed" else if (kink) spc_kind else outer_kind
      rows[[length(rows) + 1L]] <- data.frame(
        # A COVARIANCE CLASS has no single parameter. `param` on its unit is
        # the first member's, which is a convention the hyperparameter store
        # is keyed by; a table a reader reads says every equation the block
        # spans, since a shared covariance reported under one of them reads
        # as that equation's alone.
        parameter = if (is.null(u$params)) u$param else
          paste(unique(u$params), collapse = ", "),
        term = u$key, name = h, estimate = v,
        held = is_held, source = if (is.na(src)) "" else src,
        # WHY two rows carry the same number: they are shared, and the label
        # is what says so. Without it a reader sees a coincidence.
        id = lab,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) {
    return(data.frame(parameter = character(0), term = character(0),
                      name = character(0), estimate = numeric(0),
                      held = logical(0), source = character(0),
                      id = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' What Each Distribution Parameter Reached
#'
#' @description
#' Returns a single string that gives, for each parameter of the distribution,
#' the range of its fitted values. [summary.StatmodFit()] adds it as a note
#' where the fit did not converge.
#'
#' @details
#' The note is written for a lasso at a fixed hyperparameter, with a free
#' scale and a design that the model can interpolate. Fitting the
#' coefficients shrinks the residuals, which shrinks the scale and raises the
#' working weights, so that the penalty weighs relatively less and more
#' coefficients enter. A scale near zero in the note points to this case.
#'
#' The function reports where the parameters ended up and does not name a
#' cause, because naming one would require a threshold for what counts as
#' running away.
#'
#' @param x A [StatmodFit()].
#'
#' @return A single string, empty when the parameters cannot be read. It is
#'   a note of [summary.StatmodFit()] and not a line of
#'   [print.StatmodFit()], because it qualifies the fit instead of
#'   describing it.
#'
#' @seealso [statmod()]
#'
#' @keywords internal
fitted_ranges <- function(x) {
  th <- tryCatch(x@fitted, error = function(e) NULL)
  if (!length(th)) return("")
  rows <- vapply(names(th), function(p) {
    v <- as.numeric(th[[p]])
    v <- v[is.finite(v)]
    if (!length(v)) return(sprintf("  %-10s not finite", p))
    if (length(unique(v)) == 1L) {
      sprintf("  %-10s %s", p, format(signif(v[1L], 4)))
    } else {
      sprintf("  %-10s %s to %s", p, format(signif(min(v), 4)),
              format(signif(max(v), 4)))
    }
  }, character(1))
  paste0("The fit did not converge. The parameters that it reached: ",
         paste(trimws(gsub("[ ]{2,}", " ", rows)), collapse = "; "), ".")
}


#' A Term Key Shortened for Display
#'
#' @description
#' Returns the key in which the arguments of its leading call are reduced to
#' the first one, so that a trace line names the term instead of repeating its
#' whole specification.
#'
#' @details
#' A key is the call that produced the term, and for a penalty over a
#' sub-term it is that call followed by `::parameter::sub-term`. The call
#' grows when a term has many arguments, so that printing it in full
#' makes a line of an outer trace unreadable. Only the leading call is
#' shortened, and only past its first argument, so that
#' `s(x, bspline_smooth(k = 20))` and `s(z, bspline_smooth(k = 8))` stay
#' apart. Everything after `::` is kept whole, because that distinguishes
#' one entry of a term from another.
#'
#' Where shortening would make two labels the same, the full labels are
#' returned, all of them, so that every line of a trace follows one
#' convention.
#'
#' @param x A character vector of keys.
#'
#' @return A character vector the same length.
#'
#' @keywords internal
short_keys <- function(x) {
  one <- function(k) {
    parts <- strsplit(k, "::", fixed = TRUE)[[1L]]
    call <- parts[1L]
    op <- regexpr("(", call, fixed = TRUE)
    if (op < 1L || !endsWith(call, ")")) return(k)
    head_ <- substr(call, 1L, op)
    args <- substr(call, op + 1L, nchar(call) - 1L)
    if (!nzchar(args)) return(k)
    # the first argument, found by scanning for a comma at nesting depth zero
    depth <- 0L
    cut <- nchar(args) + 1L
    ch <- strsplit(args, "", fixed = TRUE)[[1L]]
    for (i in seq_along(ch)) {
      c_ <- ch[i]
      if (c_ %in% c("(", "[", "{")) depth <- depth + 1L
      else if (c_ %in% c(")", "]", "}")) depth <- depth - 1L
      else if (c_ == "," && depth == 0L) { cut <- i; break }
    }
    if (cut > nchar(args)) return(k)
    parts[1L] <- paste0(head_, substr(args, 1L, cut - 1L), ", ...)")
    paste(parts, collapse = "::")
  }
  out <- vapply(x, one, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(out) && !anyDuplicated(x)) x else out
}


#' A Titled Rule for a Verbose Trace
#'
#' @description
#' Prints one blank line, then a titled rule with the name of the method that
#' does the work on the right, so that a reader of a running fit can see where
#' one step ends and the next begins.
#'
#' @details
#' A verbose trace reports three nested levels: the outer step, the pass of
#' the alternation inside it, and the work of each block. The rule separates
#' these levels, and naming the method on every rule shows which method runs
#' at that step.
#'
#' @param title The step, e.g. `"outer 3"`.
#' @param method What runs it, or `NULL`.
#' @param indent How far in, in spaces.
#' @param char The rule's character.
#'
#' @return Invisibly `NULL`; prints.
#'
#' @keywords internal
vb_rule <- function(title, method = NULL, indent = 0L, char = "-") {
  pad <- strrep(" ", indent)
  left <- paste0(pad, strrep(char, 2L), " ", title, " ")
  right <- if (is.null(method) || !nzchar(method)) "" else
    paste0(" [", method, "]")
  fill <- max(0L, 72L - nchar(left) - nchar(right))
  cat("\n", left, strrep(char, fill), right, "\n", sep = "")
  invisible(NULL)
}

#' A Detail Line of a Verbose Trace
#'
#' @param ... Passed to `sprintf`.
#' @param indent How far in, in spaces.
#'
#' @return Invisibly `NULL`; prints.
#'
#' @keywords internal
vb_say <- function(..., indent = 5L) {
  cat(strrep(" ", indent), sprintf(...), "\n", sep = "")
  invisible(NULL)
}

#' The Name of Whatever Is About to Run
#'
#' @description
#' Returns an optimizer's own name, a criterion's kind, or a short description
#' of a step that is not an optimizer.
#'
#' @param x An optimizer, a criterion, or `NULL`.
#' @param default The string returned when `x` is `NULL` or carries neither a
#'   name nor a kind.
#'
#' @return A single string.
#'
#' @keywords internal
vb_name <- function(x, default = "") {
  if (is.null(x)) return(default)
  if (S7::S7_inherits(x) && "name" %in% names(S7::props(x))) return(x@name)
  if (S7::S7_inherits(x, Iwls)) return("iwls")
  if (S7::S7_inherits(x) && "kind" %in% names(S7::props(x))) {
    return(toupper(x@kind))
  }
  default
}
