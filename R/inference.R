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
  Vb <- solve_pd(A, "the penalized information", nm[keep])
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
#' a factor that does not exist says something about where the run stopped.
#' Returning a pseudo-inverse instead would give a standard error for a
#' direction the data does not identify.
#'
#' The message names the directions rather than the causes. A first version
#' offered two -- the run not having reached a maximum, or two columns of the
#' design carrying the same information -- and on a Student t fitted to iris
#' NEITHER was right: the design was full rank and the score was 4e-5. What had
#' happened is the third and commonest case, a parameter drifting to where its
#' information vanishes, and no list of guesses would have said so. The
#' eigenvector of the smallest eigenvalue does: it is read off and the
#' coefficients that load on it are printed.
#'
#' @param A A square matrix.
#' @param what What the matrix is, for the message.
#' @param labels The names of the coefficients \code{A} is indexed by.
#'
#' @return The inverse.
#'
#' @keywords internal
solve_pd <- function(A, what, labels = NULL) {
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(R)) return(chol2inv(R))
  stop(sprintf(paste0("%s is not positive definite, so there is no variance",
                      "\n  matrix at this point.%s\n  A fit can reach such a",
                      " point without failing: a parameter that runs\n  to",
                      " the edge of its space leaves no information behind,",
                      " and the score\n  is then small because the surface is",
                      " flat, not because it is a maximum.\n  Compare the",
                      " fitted parameters against the data before reading",
                      " anything\n  else."),
               what, flat_directions(A, labels)), call. = FALSE)
}


#' Which Coefficients a Singular Curvature Is Flat In
#'
#' @description
#' The eigenvector of the smallest eigenvalue, reported as the coefficients
#' that load on it.
#'
#' @param A A square matrix.
#' @param labels Its coefficient names.
#'
#' @return A single string, empty when nothing can be said.
#'
#' @keywords internal
flat_directions <- function(A, labels) {
  if (is.null(labels) || !length(labels)) return("")
  # a parameter far enough out makes its own row non-finite rather than merely
  # small, and then there is no eigenvector to read: the rows themselves are
  # the answer, and they are the more direct one
  bad <- labels[apply(!is.finite(A), 1L, any)]
  if (length(bad)) {
    return(sprintf(paste0("\n  Its entries are not finite in the rows of: %s",
                          "\n  which is what a parameter that has run out of",
                          " its range looks like here."),
                   paste(bad, collapse = ", ")))
  }
  e <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) return("")
  k <- which.min(e$values)
  v <- abs(e$vectors[, k])
  hit <- labels[v > 0.2 * max(v)]
  if (!length(hit)) return("")
  sprintf(paste0("\n  It is flat along a direction carried by: %s\n  (its",
                 " smallest eigenvalue is %s, against %s at the largest)."),
          paste(hit, collapse = ", "), format(signif(e$values[k], 3)),
          format(signif(max(e$values), 3)))
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


#' What Kind of Block a Term Reports As
#'
#' @description
#' Which of the four readings a term gets in a summary: its coefficients, a
#' smooth's linear part and smoothing parameter, a random effect's variance
#' parameters, or a selection's survivors.
#'
#' @details
#' The classification is by the term's class and by its penalty, not by its
#' label, so a term given a name of its own is read the same way.
#'
#' @param term A built term.
#'
#' @return One of \code{"parametric"}, \code{"smooth"}, \code{"random"},
#'   \code{"selection"}, \code{"penalized"}.
#'
#' @keywords internal
term_block_kind <- function(term) {
  pen <- modelterms7::term_penalty(term)
  if (is.null(pen)) return("parametric")
  if (S7::S7_inherits(term, modelterms7::RandomTerm)) return("random")
  if (S7::S7_inherits(term, modelterms7::SmoothTerm)) return("smooth")
  if (penalty_has_kink(pen)) return("selection")
  "penalized"
}


#' Which Coefficients of a Smooth Are the Linear Part
#'
#' @description
#' \code{TRUE} for the columns a Demmler-Reinsch smooth carries its linear
#' effect in, which are the ones worth printing.
#'
#' @details
#' The rest of the block are coefficients of an orthonormal basis of the
#' wiggly part; individually they say nothing, and what they say jointly is the
#' effective degrees of freedom, which the block header reports instead.
#'
#' The question is asked of the term's own specification (\code{spec$linear})
#' rather than of a suffix in a coefficient's name, since a name is a label and
#' this is a fact about the construction.
#'
#' @param term A built smooth term.
#' @param k The number of columns in its block.
#'
#' @return A logical vector of length \code{k}.
#'
#' @keywords internal
smooth_linear_cols <- function(term, k) {
  out <- rep(FALSE, k)
  sp <- tryCatch(term@spec, error = function(e) NULL)
  if (is.list(sp) && isTRUE(sp$linear) && k > 0L) out[1L] <- TRUE
  out
}


#' A Summary of a Fitted Model
#'
#' @description
#' What \code{\link{summary.StatmodFit}} returns: the blocks of each
#' distribution parameter, the degrees of freedom, the information criteria and
#' whatever has to be said about how the numbers should be read.
#'
#' @param call The fit's call.
#' @param distrib_name The distribution's name.
#' @param n_obs The number of observations.
#' @param tables A named list, one entry per distribution parameter, each a
#'   list of block records.
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
#' \strong{Each distribution parameter is read as blocks, not as one list of
#' coefficients}, because most of a fitted model's coefficients are not
#' quantities anybody reads. The blocks are
#' \describe{
#'   \item{the parametric terms}{every unpenalized term together, one row per
#'     coefficient, which is the ordinary table.}
#'   \item{one block per smooth}{the linear component's coefficient where the
#'     construction carries one, the smoothing parameter, and the effective
#'     degrees of freedom. The coefficients of the wiggly part are not shown:
#'     they are coordinates in an orthonormal basis and say nothing one at a
#'     time, while what they say together is the edf.}
#'   \item{one block per random effect}{the parameters of the effects'
#'     distribution -- what is usually called the variance component -- and the
#'     edf. Not the effects themselves, of which there is one per level.}
#'   \item{one block per selection}{a lasso, a SCAD or an MCP: its
#'     hyperparameters, how many coefficients survived, and those coefficients.
#'     The ones set exactly to zero are counted, not listed.}
#'   \item{one block per other penalized term}{its coefficients, which stay
#'     interpretable under a ridge, together with its hyperparameters.}
#' }
#'
#' \strong{A hyperparameter carries no standard error yet.} It is held at the
#' value it was given rather than estimated, so the row reports the value and
#' marks it fixed; inventing an interval for a number nothing estimated would
#' be worse than the empty column. Estimating them by an outer criterion is
#' what fills those rows in.
#'
#' \strong{What a Wald p-value means here depends on the row}, and the summary
#' says which is which rather than printing one column and leaving it at that.
#' For an unpenalized coefficient it is the usual thing. For a coefficient in a
#' penalized block it is conditional on the smoothing parameter, which was not
#' estimated jointly with it, and it does not account for the shrinkage of the
#' estimate towards zero. For a coefficient a kinked penalty selected, the row
#' exists only because that coefficient survived the selection, and a naive
#' interval there under-covers.
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
  design <- statmod_design(spec)
  ci$statistic <- ci$estimate / ci$se
  ci$p_value <- 2 * stats::pnorm(-abs(ci$statistic))
  lab <- coef_labels(spec, design)

  tables <- lapply(spec@distrib@params, function(p)
    summary_blocks(object, spec, design, p, ci))
  names(tables) <- spec@distrib@params

  ll <- logLik.StatmodFit(object)
  df <- attr(ll, "df")
  notes <- character(0)
  if (any(lab$penalized)) {
    notes <- c(notes, if (is.null(object@methods$outer)) paste0(
      "A hyperparameter is held at the value it was given, not estimated, so ",
      "it has\n  no standard error and every interval beside it is ",
      "conditional on it.") else sprintf(paste0(
      "A hyperparameter marked estimated was found by %s and still carries no",
      "\n  standard error: its uncertainty is not this Hessian's to give, and",
      " every\n  interval beside it is conditional on the value reached."),
      toupper(object@methods$outer@kind)))
  }
  if (any(lab$kinked)) {
    nz <- sum(lab$kinked &
                unlist(object@coefficients[spec@distrib@params],
                       use.names = FALSE) != 0)
    notes <- c(notes, sprintf(paste0(
      "%d of %d coefficients under a kinked penalty survived. Their rows ",
      "are\n  conditional on that selection, and the ones at zero have no ",
      "variance at\n  all: at the kink there is no curvature to read."),
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


#' The Blocks of One Distribution Parameter
#'
#' @description
#' Groups a parameter's terms into the readings a summary prints: the
#' parametric terms together, and one block per penalized term.
#'
#' @param fit A \code{\link{StatmodFit}}.
#' @param spec The specification.
#' @param design The design.
#' @param p The distribution parameter.
#' @param ci The flat interval table, as \code{\link{confint.StatmodFit}}
#'   returns it with the statistic and the p-value added.
#'
#' @return A list of block records, each with \code{kind}, \code{label},
#'   \code{n_coef}, \code{edf}, \code{n_zero} and \code{table}.
#'
#' @keywords internal
summary_blocks <- function(fit, spec, design, p, ci) {
  rows <- ci[ci$parameter == p, , drop = FALSE]
  cols <- c("name", "estimate", "se", "statistic", "p_value", "lower",
            "upper", "role")
  empty <- stats::setNames(
    data.frame(character(0), numeric(0), numeric(0), numeric(0), numeric(0),
               numeric(0), numeric(0), character(0),
               stringsAsFactors = FALSE), cols)

  coef_rows <- function(nm) {
    r <- rows[rows$term == nm, , drop = FALSE]
    if (!nrow(r)) return(empty)
    out <- data.frame(name = r$coefficient, estimate = r$estimate, se = r$se,
                      statistic = r$statistic, p_value = r$p_value,
                      lower = r$lower, upper = r$upper, role = "coefficient",
                      stringsAsFactors = FALSE)
    stats::setNames(out, cols)
  }
  # a hyperparameter has an estimate and nothing else, whether an outer
  # criterion found it or the caller set it: the variance of a hyperparameter
  # estimated by a marginal criterion is not this Hessian's to give, and an
  # interval here would be invented rather than computed
  outer_ran <- !is.null(fit@methods$outer)
  hyper_rows <- function(nm) {
    th <- fit@hyper[[p]][[nm]]
    if (is.null(th) || !length(th)) return(empty)
    pen <- modelterms7::term_penalty(spec@terms[[p]][[nm]])
    role <- if (outer_ran && !penalty_has_kink(pen)) "estimated" else "fixed"
    out <- data.frame(name = names(th), estimate = as.numeric(th),
                      se = NA_real_, statistic = NA_real_, p_value = NA_real_,
                      lower = NA_real_, upper = NA_real_, role = role,
                      stringsAsFactors = FALSE)
    stats::setNames(out, cols)
  }
  term_edf <- function(nm) {
    if (is.null(fit@edf)) return(NA_real_)
    e <- fit@edf
    v <- e$edf[e$parameter == p & e$term == nm]
    if (length(v)) v[1L] else NA_real_
  }

  blocks <- list()
  para <- character(0)
  for (nm in names(spec@terms[[p]])) {
    if (identical(term_block_kind(spec@terms[[p]][[nm]]), "parametric")) {
      para <- c(para, nm)
    }
  }
  if (length(para)) {
    tb <- do.call(rbind, lapply(para, coef_rows))
    blocks[[length(blocks) + 1L]] <- list(
      kind = "parametric", label = "Parametric terms", term = NA_character_,
      n_coef = nrow(tb), edf = sum(vapply(para, term_edf, numeric(1))),
      n_zero = 0L, table = tb)
  }

  for (nm in names(spec@terms[[p]])) {
    term <- spec@terms[[p]][[nm]]
    kind <- term_block_kind(term)
    if (identical(kind, "parametric")) next
    k <- length(design[[p]]$blocks[[nm]])
    cr <- coef_rows(nm)
    keep <- switch(kind,
      smooth = smooth_linear_cols(term, nrow(cr)),
      # a coefficient a kinked penalty set to zero is counted, not listed
      selection = cr$estimate != 0,
      rep(TRUE, nrow(cr)))
    if (identical(kind, "random")) keep <- rep(FALSE, nrow(cr))
    tb <- rbind(cr[keep, , drop = FALSE], hyper_rows(nm))
    blocks[[length(blocks) + 1L]] <- list(
      kind = kind,
      label = switch(kind, smooth = "Smooth", random = "Random effect",
                     selection = "Selection", "Penalized"),
      term = nm, n_coef = k, edf = term_edf(nm),
      n_zero = if (identical(kind, "selection")) sum(cr$estimate == 0) else 0L,
      table = tb)
  }
  blocks
}


#' @title Print a Model Summary
#' @name print.StatmodSummary
#' @description
#' The call, then each distribution parameter's blocks, then the degrees of
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
    cat("\n", strrep("=", 3L), " ", p, "\n", sep = "")
    blocks <- x@tables[[p]]
    if (!length(blocks)) {
      cat("  (no coefficients)\n")
      next
    }
    for (b in blocks) print_block(b, digits)
  }

  cat(sprintf("\n%.0f%% intervals, %s variance\n", 100 * x@level, x@type))
  cat(sprintf("log-likelihood %.6f    df %.2f    AIC %.3f    BIC %.3f\n",
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


#' Print One Block of a Summary
#'
#' @description
#' A header saying what the block is and what it spends, then its rows.
#'
#' @details
#' A row whose quantity was held fixed prints its value and blanks the rest,
#' rather than showing \code{NA} four times over: what the columns say is that
#' nothing estimated it, and the mark in the header says so once.
#'
#' @param b A block record, as \code{\link{summary_blocks}} returns.
#' @param digits Significant digits.
#'
#' @return \code{NULL}, invisibly.
#'
#' @keywords internal
print_block <- function(b, digits = 4L) {
  head <- if (is.na(b$term)) b$label else sprintf("%s  %s", b$label, b$term)
  bits <- character(0)
  if (!identical(b$kind, "parametric")) {
    bits <- c(bits, sprintf("%d coefficients", b$n_coef))
    if (is.finite(b$edf)) bits <- c(bits, sprintf("edf %.2f", b$edf))
  }
  if (identical(b$kind, "selection")) {
    bits <- c(bits, sprintf("%d selected, %d at zero",
                            b$n_coef - b$n_zero, b$n_zero))
  }
  cat("\n", head, sep = "")
  if (length(bits)) cat("   [", paste(bits, collapse = ", "), "]", sep = "")
  cat("\n")

  tb <- b$table
  if (!nrow(tb)) {
    cat("  (nothing to report on its own)\n")
    return(invisible(NULL))
  }
  fixed <- tb$role %in% c("fixed", "estimated")
  num <- function(v) ifelse(is.na(v), "", format(signif(v, digits)))
  out <- data.frame(
    estimate = format(signif(tb$estimate, digits)),
    se = num(tb$se),
    z = num(tb$statistic),
    p = ifelse(is.na(tb$p_value), "",
               format.pval(tb$p_value, digits = digits, eps = 1e-16)),
    lower = num(tb$lower),
    upper = num(tb$upper),
    check.names = FALSE, stringsAsFactors = FALSE)
  # said once, in the column where a standard error would have been, rather
  # than four times across a row that has nothing else in it
  out$se[fixed] <- paste0("(", tb$role[fixed], ")")
  out[fixed, c("z", "p", "lower", "upper")] <- ""
  rownames(out) <- tb$name
  print(out)
  invisible(NULL)
}
