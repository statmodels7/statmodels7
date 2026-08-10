#' @include predict.R
#' @importFrom stats vcov confint
NULL

# What a fitted model reports about its own uncertainty.
#
# The coefficients live on the LINK scale by construction -- they are the
# coefficients of a linear predictor, free on the whole line -- so a symmetric
# Wald interval is built where the quantity actually ranges and needs no
# mapping back. That is the one place this differs from
# distributions7::fit_distrib(), which estimates a bounded parameter and has to
# carry its interval through the link.

#' Where Each Stacked Coefficient Comes From
#'
#' @description
#' One row per stacked coefficient, naming the distribution parameter, the term
#' and the coefficient, and saying whether the term carries a penalty and
#' whether that penalty has a kink.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#'
#' @return A data frame with as many rows as there are coefficients.
#'
#' @keywords internal
coef_labels <- function(spec, design) {
  params <- spec@distrib@params
  rows <- list()
  for (p in params) {
    d <- design[[p]]
    if (d$npar == 0L) next
    term <- rep(NA_character_, d$npar)
    pen <- rep(FALSE, d$npar)
    kink <- rep(FALSE, d$npar)
    for (nm in names(d$blocks)) {
      cols <- d$blocks[[nm]]
      term[cols] <- nm
      tp <- modelterms7::term_penalty(spec@terms[[p]][[nm]])
      if (!is.null(tp)) {
        pen[cols] <- TRUE
        kink[cols] <- penalty_has_kink(tp)
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = p, term = term, coefficient = d$coef_names,
      penalized = pen, kinked = kink, stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    return(data.frame(parameter = character(0), term = character(0),
                      coefficient = character(0), penalized = logical(0),
                      kinked = logical(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- paste(out$parameter, out$coefficient, sep = ":")
  out
}


#' Which Information Matrix a Fit Used
#'
#' @description
#' \code{TRUE} when the fit inverted the expected information, which is what
#' \code{\link{iwls}()} does unless asked otherwise.
#'
#' @details
#' The default of \code{\link{vcov.StatmodFit}} follows this rather than
#' choosing for itself, so that a standard error comes from the same matrix the
#' fit did, and a caller who wants the other one asks for it.
#'
#' @param object A \code{\link{StatmodFit}}.
#'
#' @return A single logical.
#'
#' @keywords internal
fit_expected <- function(object) {
  m <- object@methods$smooth
  if (S7::S7_inherits(m, Iwls)) identical(m@hessian, "expected") else TRUE
}


#' @title The Variance Matrix of a Fit
#' @name vcov.StatmodFit
#' @description
#' The variance of the estimated coefficients, over every distribution
#' parameter's block at once.
#' @details
#' \strong{Two matrices, and they differ only when something is penalized.}
#' Writing \eqn{H} for the information of the log-likelihood and \eqn{S} for
#' the second derivative of the penalty,
#' \deqn{V_b = (H + S)^{-1}, \qquad V_f = (H+S)^{-1} H (H+S)^{-1}.}
#' The first is the posterior variance under the prior the penalty is the
#' negative logarithm of, and it is what an interval around a penalized term
#' should be built from: it carries the smoothing bias as though it were
#' variance, which is what makes such intervals cover at about their nominal
#' rate. The second is the sampling variance of the penalized estimator at a
#' fixed penalty, which is smaller and covers less. With no penalty \eqn{S = 0}
#' and both are \eqn{H^{-1}}.
#'
#' \strong{A coefficient a kinked penalty has set to zero has no row.} At zero
#' the penalty is not twice differentiable, so \eqn{S} does not exist there and
#' no curvature can be read; the entry is \code{NA}. The coefficients a lasso
#' or an MCP left non-zero do get a variance, and it is conditional on that
#' selection -- \code{\link{summary.StatmodFit}} says so in a note rather than
#' leaving the reader to assume otherwise.
#' @param object A \code{\link{StatmodFit}}.
#' @param type \code{"bayesian"} or \code{"frequentist"}.
#' @param expected Whether the expected information is used. Defaults to what
#'   the fit itself inverted.
#' @param ... Unused.
#' @return A square matrix over the stacked coefficients, with dimnames
#'   \code{parameter:coefficient}.
#' @seealso \code{\link{confint.StatmodFit}}, \code{\link{summary.StatmodFit}}
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' sqrt(diag(vcov(fit)))
#' @keywords internal
vcov.StatmodFit <- function(object, type = c("bayesian", "frequentist"),
                            expected = NULL, ...) {
  type <- match.arg(type)
  if (is.null(expected)) expected <- fit_expected(object)
  spec <- object@spec
  design <- statmod_design(spec)
  coef <- object@coefficients
  lab <- coef_labels(spec, design)
  total <- nrow(lab)
  nm <- rownames(lab)

  H <- statmod_information_at(spec, coef, design, expected)
  S <- statmod_penalty_at(spec, coef, object@hyper, design, "hessian")

  keep <- rep(TRUE, total)
  beta <- unlist(coef[spec@distrib@params], use.names = FALSE)
  keep[lab$kinked & beta == 0] <- FALSE
  # a kinked penalty contributes no curvature away from its kink either, and
  # any non-finite entry would be the kink itself reached by a hair
  S[!is.finite(S)] <- 0

  out <- matrix(NA_real_, total, total, dimnames = list(nm, nm))
  if (!any(keep)) return(out)
  A <- (H + S)[keep, keep, drop = FALSE]
  Vb <- solve_pd(A, "the penalized information")
  V <- if (type == "bayesian") Vb else
    Vb %*% H[keep, keep, drop = FALSE] %*% Vb
  out[keep, keep] <- V
  out
}
S7::method(vcov, StatmodFit) <- vcov.StatmodFit


#' Invert a Matrix That Ought to Be Positive Definite
#'
#' @description
#' A Cholesky inverse, signalling an error naming the matrix when the factor
#' does not exist.
#'
#' @details
#' A failure here is a statement about the fit rather than about the
#' arithmetic: at a maximum the penalized information is positive definite, so
#' a factor that does not exist says the run stopped somewhere that is not one,
#' or that two columns of the design carry the same information. Returning a
#' pseudo-inverse instead would give a standard error for a direction the data
#' does not identify.
#'
#' @param A A square matrix.
#' @param what What the matrix is, for the message.
#'
#' @return The inverse.
#'
#' @keywords internal
solve_pd <- function(A, what) {
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(R)) {
    stop(sprintf(paste0("%s is not positive definite, so there is no\n",
                        "  variance matrix at this point. Either the fit did",
                        " not reach a\n  maximum -- check its convergence --",
                        " or two columns of the design\n  carry the same",
                        " information."), what), call. = FALSE)
  }
  chol2inv(R)
}


#' @title Confidence Intervals for a Fit
#' @name confint.StatmodFit
#' @description
#' Wald intervals for the coefficients of every distribution parameter.
#' @details
#' The interval is symmetric about the estimate and needs no mapping back. A
#' coefficient of a linear predictor is unbounded whatever the distribution
#' parameter it belongs to, the link having already carried that parameter onto
#' the whole line, so the scale the interval is built on is the scale the
#' quantity lives on. What the interval does not do is respect a bound on the
#' parameter itself; for that, map an interval for the predictor through the
#' inverse link at the covariate values of interest.
#'
#' The variance comes from \code{\link{vcov.StatmodFit}}, so the same two
#' conventions apply, and a coefficient a kinked penalty set to zero has
#' \code{NA} rather than an interval.
#' @param object A \code{\link{StatmodFit}}.
#' @param parm Which coefficients: a distribution parameter's name, a vector of
#'   \code{parameter:coefficient} labels, or \code{NULL} for all of them.
#' @param level The confidence level.
#' @param type Passed to \code{\link{vcov.StatmodFit}}.
#' @param ... Passed to \code{\link{vcov.StatmodFit}}.
#' @return A data frame with the parameter, the term, the coefficient, the
#'   estimate, its standard error and the two limits.
#' @seealso \code{\link{vcov.StatmodFit}}, \code{\link{summary.StatmodFit}}
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' confint(fit)
#' confint(fit, "sigma")
#' @keywords internal
confint.StatmodFit <- function(object, parm = NULL, level = 0.95,
                               type = c("bayesian", "frequentist"), ...) {
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  spec <- object@spec
  design <- statmod_design(spec)
  lab <- coef_labels(spec, design)
  est <- unlist(object@coefficients[spec@distrib@params], use.names = FALSE)
  se <- sqrt(diag(vcov(object, type = type, ...)))
  z <- stats::qnorm((1 + level) / 2)

  out <- data.frame(lab[, c("parameter", "term", "coefficient")],
                    estimate = est, se = se,
                    lower = est - z * se, upper = est + z * se,
                    stringsAsFactors = FALSE)
  rownames(out) <- rownames(lab)
  if (is.null(parm)) return(out)

  keep <- rownames(out) %in% parm | out$parameter %in% parm |
    out$term %in% parm
  if (!any(keep)) {
    stop(sprintf(paste0("'parm' matched nothing. It takes a distribution",
                        " parameter\n  (%s), a term's name, or a label of",
                        " the form 'parameter:coefficient'."),
                 paste(spec@distrib@params, collapse = ", ")), call. = FALSE)
  }
  out[keep, , drop = FALSE]
}
S7::method(confint, StatmodFit) <- confint.StatmodFit


#' A Summary of a Fitted Model
#'
#' @description
#' What \code{\link{summary.StatmodFit}} returns: one coefficient table per
#' distribution parameter, the degrees of freedom, the information criteria and
#' whatever has to be said about how the numbers should be read.
#'
#' @param call The fit's call.
#' @param distrib_name The distribution's name.
#' @param n_obs The number of observations.
#' @param tables A named list of data frames, one per distribution parameter.
#' @param edf The per-term degrees of freedom.
#' @param loglik The maximized log-likelihood.
#' @param df The effective degrees of freedom in total.
#' @param aic,bic The information criteria.
#' @param converged Whether every loop stopped on its own rule.
#' @param elapsed The elapsed time, in seconds.
#' @param level The confidence level the intervals were built at.
#' @param type Which variance convention was used.
#' @param notes Character vector of things the reader has to know.
#'
#' @return An object of class \code{StatmodSummary}.
#'
#' @seealso \code{\link{summary.StatmodFit}}
#'
#' @examples
#' dd <- data.frame(y = rnorm(30), x = runif(30))
#' S7::S7_inherits(summary(statmod(y ~ x,
#'                                 distributions7::gaussian1_distrib(), dd)),
#'                 StatmodSummary)
#'
#' @name StatmodSummary-class
#' @aliases StatmodSummary
#' @keywords internal
#' @export
StatmodSummary <- S7::new_class("StatmodSummary",
  properties = list(
    call = S7::class_any,
    distrib_name = S7::class_character,
    n_obs = S7::class_numeric,
    tables = S7::class_list,
    edf = S7::class_any,
    loglik = S7::class_numeric,
    df = S7::class_numeric,
    aic = S7::class_numeric,
    bic = S7::class_numeric,
    converged = S7::class_logical,
    elapsed = S7::class_numeric,
    level = S7::class_numeric,
    type = S7::class_character,
    notes = S7::class_character
  )
)


#' @title Summarize a Fitted Model
#' @name summary.StatmodFit
#' @description
#' One coefficient table per distribution parameter -- estimate, standard
#' error, Wald statistic, p-value and interval -- with the degrees of freedom,
#' the information criteria and the qualifications the numbers carry.
#' @details
#' \strong{What a Wald p-value means here depends on the row}, and the summary
#' says which is which rather than printing one column and leaving it at that.
#' For an unpenalized coefficient it is the usual thing. For a coefficient in a
#' penalized block it is conditional on the smoothing parameter, which was not
#' estimated jointly with it, and it does not account for the shrinkage of the
#' estimate towards zero. For a block a kinked penalty selected -- a lasso, a
#' SCAD, an MCP -- the row exists only because that coefficient survived the
#' selection, and a naive interval there under-covers; the coefficients set
#' exactly to zero carry \code{NA}, since at the kink there is no curvature to
#' read.
#'
#' \strong{The degrees of freedom} are the effective ones, summed over the
#' terms, so that a penalized term counts what it spends rather than how many
#' columns it has. The information criteria are built on that count.
#' @param object A \code{\link{StatmodFit}}.
#' @param level The confidence level.
#' @param type Which variance matrix: passed to \code{\link{vcov.StatmodFit}}.
#' @param ... Passed to \code{\link{vcov.StatmodFit}}.
#' @return A \code{\link{StatmodSummary}}.
#' @seealso \code{\link{vcov.StatmodFit}}, \code{\link{confint.StatmodFit}}
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120))
#' dd$y <- 1 + 2 * dd$x + rnorm(120, sd = 0.4)
#' summary(statmod(y ~ x | sigma ~ x,
#'                 distributions7::gaussian1_distrib(), dd))
#' @keywords internal
summary.StatmodFit <- function(object, level = 0.95,
                               type = c("bayesian", "frequentist"), ...) {
  type <- match.arg(type)
  ci <- confint(object, level = level, type = type, ...)
  spec <- object@spec
  lab <- coef_labels(spec, statmod_design(spec))

  zval <- ci$estimate / ci$se
  ci$statistic <- zval
  ci$p_value <- 2 * stats::pnorm(-abs(zval))
  ci$penalized <- lab$penalized
  ci$selected <- lab$kinked

  tables <- lapply(spec@distrib@params, function(p) {
    r <- ci[ci$parameter == p, , drop = FALSE]
    r[, c("term", "coefficient", "estimate", "se", "statistic", "p_value",
          "lower", "upper", "penalized", "selected")]
  })
  names(tables) <- spec@distrib@params

  ll <- logLik.StatmodFit(object)
  df <- attr(ll, "df")
  notes <- character(0)
  if (any(lab$penalized & !lab$kinked)) {
    notes <- c(notes, paste0(
      "A penalized coefficient's interval is conditional on its smoothing ",
      "parameter,\n  which was held fixed and not estimated jointly with it."))
  }
  if (any(lab$kinked)) {
    nz <- sum(lab$kinked &
                unlist(object@coefficients[spec@distrib@params],
                       use.names = FALSE) != 0)
    notes <- c(notes, sprintf(paste0(
      "%d of %d coefficients under a kinked penalty are non-zero. Their ",
      "rows are\n  conditional on that selection, and the zeros carry NA: ",
      "at the kink there\n  is no curvature to read."),
      nz, sum(lab$kinked)))
  }
  if (!object@converged) {
    notes <- c(notes, paste0(
      "The fit did not converge, so everything below is read at a point that",
      "\n  is not a maximum."))
  }

  StatmodSummary(
    call = object@call, distrib_name = spec@distrib@distrib_name,
    n_obs = spec@n_obs, tables = tables, edf = object@edf,
    loglik = object@loglik, df = df,
    aic = -2 * object@loglik + 2 * df,
    bic = -2 * object@loglik + log(spec@n_obs) * df,
    converged = object@converged, elapsed = object@elapsed,
    level = level, type = type, notes = notes)
}
S7::method(summary, StatmodFit) <- summary.StatmodFit


#' @title Print a Model Summary
#' @name print.StatmodSummary
#' @description
#' The call, one coefficient table per distribution parameter, the degrees of
#' freedom, the criteria and the notes.
#' @param x A \code{\link{StatmodSummary}}.
#' @param digits Significant digits in the tables.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{summary.StatmodFit}}
#' @keywords internal
print.StatmodSummary <- function(x, digits = 4L, ...) {
  cat("A statmod fit\n\n")
  cat("Call:  ", paste(deparse(x@call), collapse = "\n        "), "\n\n",
      sep = "")
  cat("Distribution: ", x@distrib_name, "     Observations: ", x@n_obs,
      "\n", sep = "")

  for (p in names(x@tables)) {
    tb <- x@tables[[p]]
    cat("\n", p, "\n", sep = "")
    if (!nrow(tb)) {
      cat("  (no coefficients)\n")
      next
    }
    mark <- ifelse(tb$selected, "*", ifelse(tb$penalized, "+", " "))
    out <- data.frame(
      estimate = signif(tb$estimate, digits),
      se = signif(tb$se, digits),
      z = signif(tb$statistic, digits),
      p = format.pval(tb$p_value, digits = digits, eps = 1e-16),
      lower = signif(tb$lower, digits),
      upper = signif(tb$upper, digits),
      ` ` = mark, check.names = FALSE)
    rownames(out) <- paste(tb$term, tb$coefficient, sep = " / ")
    print(out)
  }

  cat(sprintf("\n%.0f%% intervals, %s variance\n", 100 * x@level, x@type))
  if (!is.null(x@edf) && any(x@edf$edf != x@edf$coefficients)) {
    cat("\nDegrees of freedom, per term\n")
    print(x@edf, row.names = FALSE)
  }
  cat(sprintf("\nlog-likelihood %.6f    df %.2f    AIC %.3f    BIC %.3f\n",
              x@loglik, x@df, x@aic, x@bic))
  cat(sprintf("fitted in %s, %s\n", format_duration(x@elapsed),
              if (x@converged) "converged" else "DID NOT CONVERGE"))
  if (length(x@notes)) {
    cat("\n")
    for (n in x@notes) cat("  ", n, "\n", sep = "")
  }
  invisible(x)
}
S7::method(print, StatmodSummary) <- print.StatmodSummary
