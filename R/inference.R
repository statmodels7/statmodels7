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
  units <- statmod_penalized(spec, design)
  rows <- list()
  for (p in params) {
    d <- design[[p]]
    if (d$npar == 0L) next
    term <- rep(NA_character_, d$npar)
    pen <- rep(FALSE, d$npar)
    kink <- rep(FALSE, d$npar)
    for (nm in names(d$blocks)) term[d$blocks[[nm]]] <- nm
    for (u in units) {
      if (!identical(u$param, p)) next
      # a structural unit's penalty sits on the term's OWN parameters,
      # which contribute no column: its cols index the term's parameter
      # vector, and writing them here grew pen past the design and
      # recycled the labels into duplicate rows
      if (isTRUE(u$structural)) next
      pen[u$cols] <- TRUE
      kink[u$cols] <- penalty_has_kink(u$penalty)
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

  # A model carrying a structural term is inverted over the coefficients AND
  # the term's own parameters together, and the coefficient block of that
  # inverse is taken. Inverting the coefficient block alone would report the
  # variance that would hold with the term's parameters known, which is not
  # what was estimated. The information there is the observed one: neither a
  # filter nor a mixture over states has an expected information to offer.
  fil <- length(attr(design, "structural")) > 0L
  H <- if (fil) statmod_full_information(spec, coef, design) else
    statmod_information_at(spec, coef, design, expected)
  S <- statmod_penalty_at(spec, coef, object@hyper, design, "hessian")
  nz <- nrow(H) - total
  if (nz > 0L) {
    S <- rbind(cbind(S, matrix(0, total, nz)),
               matrix(0, nz, total + nz))
    # A penalty over the term's OWN parameters belongs in that tail. Leaving
    # it out is not a conservative choice: the deviations of a panel are not
    # identified without it -- a constant added to a population value and
    # taken off every deviation leaves the filter unchanged -- so the matrix
    # is singular along exactly that direction and nothing can be reported.
    ps <- structural_penalty_block(spec, design, object@hyper, nz)
    if (!is.null(ps)) S[total + seq_len(nz), total + seq_len(nz)] <- ps
  }

  keep <- rep(TRUE, total)
  beta <- unlist(coef[spec@distrib@params], use.names = FALSE)
  keep[lab$kinked & beta == 0] <- FALSE
  # a kinked penalty contributes no curvature away from its kink either, and
  # any non-finite entry would be the kink itself reached by a hair
  S[!is.finite(S)] <- 0

  out <- matrix(NA_real_, total, total, dimnames = list(nm, nm))
  if (!any(keep)) return(out)
  keep_full <- c(keep, rep(TRUE, nz))
  A <- (H + S)[keep_full, keep_full, drop = FALSE]
  # the unpenalized information supplies the reference for what a SMALL
  # eigenvalue means: a smoothing parameter at 1e15 separates the scales
  # without flattening any direction, and against max(ev) alone that read
  # as singularity
  Hd <- diag(as_dense(H))[keep_full]
  Vb <- solve_pd(A, "the penalized information",
                 c(nm[keep], rep("", nz)),
                 scale = max(abs(Hd), na.rm = TRUE))
  V <- if (type == "bayesian") Vb else
    Vb %*% H[keep_full, keep_full, drop = FALSE] %*% Vb
  # the coefficient block of the joint inverse, which is not the inverse of
  # the coefficient block wherever the two are correlated
  V <- V[seq_len(sum(keep)), seq_len(sum(keep)), drop = FALSE]
  out[keep, keep] <- V
  out
}
S7::method(vcov, StatmodFit) <- vcov.StatmodFit


#' Invert a Matrix That Ought to Be Positive Definite
#'
#' @description
#' An inverse through the Cholesky factor, signalling an error naming the
#' matrix when a direction is flat.
#'
#' @details
#' A failure here is a statement about the fit rather than about the
#' arithmetic: at a maximum the penalized information is positive definite, so
#' a matrix that is not says something about where the run stopped. The test is
#' \code{lmin > tol * ref} on the smallest eigenvalue rather than whether
#' \code{chol()} raised, because on an exactly singular matrix the latter is
#' decided by rounding and differs between platforms; \code{ref} is the
#' matrix's own scale, or the scale of the unpenalized
#' information where the caller holds it, which is what tells a flat
#' direction from the scale separation a large smoothing parameter
#' legitimately produces. Returning a pseudo-inverse instead would give a
#' standard error for a direction the data does not identify.
#'
#' The smallest eigenvalue is ESTIMATED rather than computed, from LAPACK's
#' condition estimator (\code{dpocon}) read on the Cholesky factor the
#' inverse needs anyway: \code{rcond} is
#' \eqn{1/(\lVert A\rVert_1\lVert A^{-1}\rVert_1)}, so
#' \code{rcond * ||A||_1} is \eqn{1/\lVert A^{-1}\rVert_1}, which for a
#' symmetric matrix lies between \eqn{\lambda_{\min}/\sqrt{p}} and
#' \eqn{\lambda_{\min}}. The estimate therefore errs on the SMALL side and
#' the test is conservative by at most a factor \eqn{\sqrt{p}}, plus
#' whatever the estimator's own slack is; the two cases it has to keep
#' apart are separated by some fifty orders of magnitude, so neither
#' reaches the other. It replaced a full eigendecomposition, which answers
#' the same question exactly and costs O(p^3) with a large constant --
#' measured at p = 1022, 1.18 s against the Cholesky's 0.25.
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
#' @param scale An optional reference magnitude for the smallest
#'   eigenvalue, usually the largest diagonal entry of the UNPENALIZED
#'   information. Without it the reference is the matrix's own scale.
#'
#' @return The inverse.
#'
#' @keywords internal
solve_pd <- function(A, what, labels = NULL, scale = NULL) {
  # The verdict comes from the smallest eigenvalue and not from whether
  # chol() raised. On a matrix with an exactly zero eigenvalue -- two columns
  # of the design carrying the same information is the ordinary way to get
  # one -- the pivot that should be zero comes out positive or negative
  # according to rounding, so the same fit was refused on Windows and
  # accepted on Linux and macOS. Reading the condition estimate off the
  # factor is a statement about the matrix instead: where the factorization
  # succeeds on a singular matrix by luck, the estimate is at the rounding
  # scale and the answer is the same on every platform.
  #
  # The REFERENCE distinguishes two situations one ratio conflated. A
  # smoothing parameter a criterion legitimately sends to 1e15 puts
  # min(ev)/max(ev) at the rounding scale while the small eigenvalue is
  # ordinary curvature (measured: min 29.7 against max 6.4e15 on a poisson
  # smooth over weak signal, the matrix strictly positive definite and the
  # fit right); a FLAT direction is small against the unpenalized
  # information's own scale either way. A caller holding that information
  # passes its magnitude; the reference never exceeds the largest
  # eigenvalue, so the test only ever relaxes towards it.
  #
  # A non-finite entry is what a parameter run out of its range leaves behind,
  # and it is not a matrix to decompose: the factorization raises its own
  # error there, which would replace the message below with one naming
  # neither the fit nor the direction. flat_directions() reports those rows
  # instead.
  # The test and the message's flat direction both need a dense matrix, and
  # there is no cheap sparse counterpart to replace them with. Densifying
  # HERE is deliberate rather than a lapse: this runs once, at vcov(), on a
  # p by p matrix, where the sparsity that matters is in the n by p design
  # and the per-iteration products taken against it. The same judgement is
  # recorded for the observed Hessian of a regime mixture, which is also
  # computed once and left in R.
  A <- as_dense(A)
  if (ncol(A) > 0L && all(is.finite(A))) {
    ch <- tryCatch(chol(A), error = function(e) NULL)
    if (!is.null(ch)) {
      # the 1-norm, which is what the condition estimator is expressed in
      anorm <- max(colSums(abs(A)))
      rc <- chol_rcond_cpp(ch, anorm)
      lmin <- if (is.na(rc)) 0 else rc * anorm
      ref <- if (is.null(scale) || !is.finite(scale) || scale <= 0) anorm
        else min(scale, anorm)
      if (is.finite(lmin) && lmin > 1e-12 * ref) {
        # Inverted through the factor already in hand. A condition number of
        # 1e15 born of scale separation costs the well-determined directions
        # nothing here and the shrunk ones simply report variances near zero.
        out <- chol2inv(ch)
        return((out + t(out)) / 2)
      }
    }
  }
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


#' The Log-Determinant of a Penalized Information, Robustly
#'
#' @description
#' \eqn{\log|M|} for a matrix a Laplace approximation needs to be positive
#' definite, by the cheap route where that is safe and a costlier one where
#' it is not.
#'
#' @details
#' \strong{Why not \code{chol()} alone.} A marginal criterion read the
#' determinant off \code{chol(M)} and reported the criterion as NONEXISTENT
#' whenever the factorization raised. At a condition number near the rounding
#' floor whether it raises is decided by arithmetic and not by the matrix:
#' measured on a hierarchical score-driven panel, \eqn{K+S} had a smallest
#' eigenvalue of 4.3e-11 against a condition number of 8.0e15, and the outer
#' search then backtracked through a dozen points reported unavailable towards
#' one that had been available a moment earlier. The same doubt this package
#' already records for \code{\link{solve_pd}} and for basis7's rank tests.
#'
#' \strong{The three routes.} The factorization is tried first, being O(p^3/3)
#' and the common case. Where it succeeds, LAPACK's condition estimator reads
#' the smallest eigenvalue off the factor already in hand for O(p^2), and a
#' matrix comfortably away from the floor is accepted with the determinant the
#' factor gives. Only where that test is inconclusive -- or the factorization
#' raised at all -- is the eigendecomposition computed, which costs more and
#' answers about the MATRIX: a factorization that failed by rounding luck on a
#' matrix that is in fact positive definite is recovered there, and one that
#' is genuinely rank deficient is refused deterministically rather than
#' according to the platform.
#'
#' @param M A symmetric matrix.
#' @param scale A reference magnitude, as \code{\link{solve_pd}} takes: the
#'   unpenalized information's own scale, so that a hyperparameter legitimately
#'   sent to 1e15 is told apart from a flat direction.
#'
#' @return A list with \code{logdet} and \code{ok}, or \code{ok = FALSE} where
#'   the matrix is not positive definite.
#'
#' @seealso \code{\link{solve_pd}}, \code{\link{statmod_marginal}}
#'
#' @keywords internal
pd_logdet <- function(M, scale = NULL) {
  M <- as_dense(M)
  if (!ncol(M) || !all(is.finite(M))) return(list(logdet = NA_real_, ok = FALSE))
  anorm <- max(colSums(abs(M)))
  ref <- if (is.null(scale) || !is.finite(scale) || scale <= 0) anorm
         else min(scale, anorm)
  # ⚠️ The two thresholds are NOT solve_pd()'s, and the difference is the
  # quantity. Inverting a matrix at a condition number of 1e14 loses most of
  # the answer, so solve_pd() refuses there; a log-DETERMINANT is a sum of
  # logarithms and survives it -- measured, chol and eigen agree to six
  # significant figures at 1e14 (-644.725631 against -644.725128). What
  # matters here is only that every eigenvalue is genuinely POSITIVE, and an
  # eigenvalue below eps times the largest is not distinguishable from zero
  # in double precision whatever the factorization reports.
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(ch)) {
    rc <- tryCatch(chol_rcond_cpp(ch, anorm), error = function(e) NA_real_)
    lmin <- if (is.na(rc)) 0 else rc * anorm
    # comfortably clear of the floor: chol's answer is its own evidence, and
    # this is the common case, so the eigendecomposition is never computed
    if (is.finite(lmin) && lmin > sqrt(.Machine$double.eps) * ref) {
      return(list(logdet = 2 * sum(log(diag(ch))), ok = TRUE))
    }
  }
  # inconclusive or refused: ask the MATRIX rather than the arithmetic
  ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || !all(is.finite(ev))) {
    return(list(logdet = NA_real_, ok = FALSE))
  }
  if (min(ev) > .Machine$double.eps * max(ev)) {
    return(list(logdet = sum(log(ev)), ok = TRUE))
  }
  list(logdet = NA_real_, ok = FALSE, min_ev = min(ev), max_ev = max(ev))
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
#' The classification is by the term's class and by its penalties, not by its
#' label, so a term given a name of its own is read the same way. The
#' penalties are the ones the term declares through
#' \code{\link[modelterms7]{term_penalties}}, so a term penalized over part of
#' its parameters -- a segmented term's changes, a filter's deviations -- is
#' read as penalized rather than as parametric, and is a selection when any
#' of its penalties has a kink.
#'
#' @param term A built term.
#'
#' @return One of \code{"parametric"}, \code{"smooth"}, \code{"random"},
#'   \code{"selection"}, \code{"penalized"}.
#'
#' @keywords internal
term_block_kind <- function(term) {
  # a break-point term is asked about FIRST, before its penalties are
  # looked at: its block is a working linearization whose coefficients are
  # not the quantities of the model, so it wants a section of its own
  # whether or not a development of its coefficients carries a penalty
  if (S7::S7_inherits(term, modelterms7::SegTerm)) return("breakpoint")
  ent <- modelterms7::term_penalties(term)
  if (!length(ent)) return("parametric")
  if (S7::S7_inherits(term, modelterms7::RandomTerm)) return("random")
  if (S7::S7_inherits(term, modelterms7::SmoothTerm)) return("smooth")
  if (any(vapply(ent, function(e) penalty_has_kink(e$penalty), logical(1)))) {
    return("selection")
  }
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
    structural = S7::class_any,
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
#' \strong{A hyperparameter is the first row of its block}, since it governs
#' every coefficient under it, and the cell where its standard error would be
#' says what put the value there. One estimated by \code{\link{reml}()} or
#' \code{\link{ml}()} maximizes a twice differentiable criterion, so it carries
#' a standard error and an interval, both read on the free scale its link
#' defines and mapped back (\code{\link{statmod_hyper_vcov}}). One chosen by
#' \code{\link{aic}()}, \code{\link{bic}()} or \code{\link{cv}()} over a kinked
#' penalty is the argument of a minimum over a grid, so the row names the
#' criterion and leaves the remaining columns empty: there is no curvature at
#' such a point to read a standard error from. One the caller set is marked
#' fixed.
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
#' @param correct Whether the degrees of freedom carry what the estimation
#'   of the hyperparameters cost. The ordinary count reads them as known,
#'   and they were chosen from the same data, so a criterion built on it is
#'   too generous. See \code{\link{statmod_edf_correction}}. Defaults to
#'   \code{FALSE} because it changes a number a reader may be comparing with
#'   an earlier fit; it is zero where no hyperparameter was estimated.
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
                               type = c("bayesian", "frequentist"),
                               correct = FALSE, ...) {
  type <- match.arg(type)
  ci <- confint(object, level = level, type = type, ...)
  spec <- object@spec
  design <- statmod_design(spec)
  ci$statistic <- ci$estimate / ci$se
  ci$p_value <- 2 * stats::pnorm(-abs(ci$statistic))
  lab <- coef_labels(spec, design)

  # one variance matrix for every block that needs one: a term reported by
  # the quantities it is about carries them across by the delta method
  V <- tryCatch(vcov(object, type = type, ...), error = function(e) NULL)
  tables <- lapply(spec@distrib@params, function(p)
    summary_blocks(object, spec, design, p, ci, level, V))
  names(tables) <- spec@distrib@params

  ll <- logLik.StatmodFit(object)
  df <- attr(ll, "df")
  notes <- character(0)

  # The count above reads the hyperparameters as though they were known,
  # and they were estimated from the same data. Adding what that costs is
  # off by default because it changes a number a reader may be comparing
  # with an earlier fit.
  corr <- 0
  if (isTRUE(correct)) {
    cc <- tryCatch(statmod_edf_correction(spec, object@coefficients,
                                          object@hyper, design,
                                          object@methods$outer),
                   error = function(e) list(total = 0, per = numeric(0)))
    corr <- cc$total
    df <- df + corr
    notes <- c(notes, if (corr > 0) sprintf(paste0(
      "The degrees of freedom carry %.3f for the hyperparameters, which ",
      "were\n  estimated rather than given. Without it the criteria are ",
      "too generous."), corr) else paste0(
      "No hyperparameter here was estimated by a marginal criterion, so ",
      "there is\n  nothing for the correction to propagate and it is zero."))
  }
  if (any(lab$penalized)) {
    # WHICH criteria were in force, so a reader knows whether the number at
    # the head of a penalized block was chosen or given. The two kinds carry
    # different guarantees and are named separately.
    src <- unique(unlist(lapply(tables, function(bl)
      unlist(lapply(bl, function(b) b$table$source[b$table$role %in%
                                                     c("fixed", "estimated")]),
             use.names = FALSE)), use.names = FALSE))
    marg <- setdiff(intersect(src, c("reml", "ml")), NA)
    path <- setdiff(intersect(src, c("aic", "bic", "cv")), NA)
    if (length(marg)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s was estimated by that criterion, and its",
        " standard\n  error and interval are read on the free scale its link ",
        "defines and mapped\n  back. Every coefficient beside it is still ",
        "conditional on the value reached."),
        paste(toupper(marg), collapse = " or ")))
    }
    if (length(path)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s was chosen by a path over its own values,",
        " not held.\n  It is the argument of a minimum over a grid rather ",
        "than the root of a\n  derivative, so no standard error follows from ",
        "it; its uncertainty is a\n  resampling question."),
        paste(toupper(path), collapse = " or ")))
    }
    if ("fixed" %in% src) {
      notes <- c(notes, paste0(
        "A hyperparameter marked fixed is held at the value it was given, not",
        " estimated,\n  so it has no standard error and every interval beside",
        " it is conditional on it."))
    }
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

  # A structural term contributes no columns, so nothing above can report
  # it and its parameters were reachable only through fit@structural.
  strc <- tryCatch(statmod_structural_table(object, level),
                   error = function(e) NULL)
  if (!is.null(strc) && any(strc$held)) {
    notes <- c(notes, paste0(
      "A level marked held is carried by an intercept in the same equation ",
      "and is\n  not estimated: the two are exactly confounded, so only one ",
      "of them can be."))
  }

  StatmodSummary(
    call = object@call, distrib_name = spec@distrib@distrib_name,
    n_obs = spec@n_obs, tables = tables, edf = object@edf,
    structural = strc,
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
#' @param level The confidence level the intervals are built at.
#' @param V The variance matrix over the stacked coefficients, or
#'   \code{NULL}. It is needed only by a term reported through
#'   \code{\link[modelterms7]{term_readable}}, whose quantities are
#'   functions of several coefficients at once and whose standard errors
#'   are therefore the delta method rather than a diagonal entry.
#'
#' @return A list of block records, each with \code{kind}, \code{label},
#'   \code{n_coef}, \code{edf}, \code{n_zero} and \code{table}.
#'
#' @keywords internal
summary_blocks <- function(fit, spec, design, p, ci, level = 0.95,
                           V = NULL) {
  rows <- ci[ci$parameter == p, , drop = FALSE]
  cols <- c("name", "estimate", "se", "statistic", "p_value", "lower",
            "upper", "role")
  # `source` says what put the number there -- a criterion by name, or the
  # caller -- which `role` alone cannot, every criterion answering "estimated"
  all_cols <- c(cols, "source")
  empty <- stats::setNames(
    data.frame(character(0), numeric(0), numeric(0), numeric(0), numeric(0),
               numeric(0), numeric(0), character(0), character(0),
               stringsAsFactors = FALSE), all_cols)

  coef_rows <- function(nm) {
    r <- rows[rows$term == nm, , drop = FALSE]
    if (!nrow(r)) return(empty)
    out <- data.frame(name = r$coefficient, estimate = r$estimate, se = r$se,
                      statistic = r$statistic, p_value = r$p_value,
                      lower = r$lower, upper = r$upper, role = "coefficient",
                      source = "", stringsAsFactors = FALSE)
    stats::setNames(out, all_cols)
  }
  # WHICH hyperparameter was estimated and by what. A marginal criterion --
  # reml(), ml() -- maximizes a twice differentiable function of it, so it
  # carries a standard error and an interval, both read on the free scale its
  # link defines and mapped back (statmod_hyper_vcov). A path -- aic(), bic(),
  # cv() over a kinked penalty -- chooses the argument of a minimum over a
  # grid, so there is no curvature to read and the row reports the value and
  # the criterion that chose it. Only a value the caller set is held.
  outer_ran <- !is.null(fit@methods$outer)
  Vh <- if (outer_ran) tryCatch(
    statmod_hyper_vcov(spec, design, fit@coefficients, fit@hyper,
                       fit@methods$outer), error = function(e) NULL) else NULL
  spc <- fit@methods$sparse_criterion
  spc_keys <- fit@methods$sparse_hyper
  if (is.null(spc_keys)) spc_keys <- character(0)
  # A term may carry more than one penalty, each filed under a key of its
  # own, so the rows of a term are those of every key belonging to it. Where
  # there are several the hyperparameter is named for the penalty as well:
  # two lambdas in one block, one on the slope changes and one on the jumps,
  # are not the same number and cannot appear under the same name.
  hyper_rows <- function(nm) {
    ent <- modelterms7::term_penalties(spec@terms[[p]][[nm]])
    if (!length(ent)) return(empty)
    des <- statmod_design(spec)
    out <- lapply(ent, function(e) {
      key <- statmod_entry_key(nm, ent, e)
      th <- fit@hyper[[p]][[key]]
      if (is.null(th) || !length(th)) return(empty)
      u <- statmod_unit(spec, des, p, key)
      if (is.null(u)) return(empty)
      lab <- if (length(ent) > 1L && nzchar(e$name)) {
        paste(e$name, names(th), sep = ".")
      } else {
        names(th)
      }
      marginal <- outer_ran && !penalty_has_kink(u$penalty)
      hk <- paste(p, key, names(th), sep = "\r")
      # PER HYPERPARAMETER, and the test has to be the vector one: `ifelse`
      # returns a result the length of its TEST, so a scalar `marginal`
      # returned one answer and recycled it over a penalty carrying two
      # hyperparameters -- the elastic net's alpha, which no path varies,
      # was reported as chosen by the criterion that had chosen its lambda
      chosen <- rep(marginal, length(hk)) | hk %in% spc_keys
      role <- ifelse(chosen, "estimated", "fixed")
      by <- if (marginal) fit@methods$outer@kind else
        if (is.null(spc)) "estimated" else spc@kind
      src <- ifelse(chosen, by, "fixed")
      r <- data.frame(name = lab, estimate = as.numeric(th),
                      se = NA_real_, statistic = NA_real_, p_value = NA_real_,
                      lower = NA_real_, upper = NA_real_, role = role,
                      source = src, stringsAsFactors = FALSE)
      # the interval is built where the criterion was maximized and mapped
      # back, as every other interval in the toolkit is, so a positive
      # hyperparameter cannot be given a negative lower end; the standard
      # error printed beside it is the delta method onto its own scale. No
      # test accompanies it: the null a z of value/se reports on is that the
      # hyperparameter is zero, which for a smoothing parameter is the edge
      # of its range and not an interior hypothesis.
      if (!is.null(Vh) && marginal) {
        z <- stats::qnorm(1 - (1 - level) / 2)
        lk <- attr(Vh, "idx")
        for (i in seq_along(th)) {
          k <- paste(p, key, names(th)[[i]], sep = "\r")
          j <- match(k, rownames(Vh))
          if (is.na(j)) next
          link <- attr(lk, "links")[[j]]
          se_eta <- sqrt(Vh[j, j])
          eta <- linkfunctions7::linkfun(link, r$estimate[[i]])
          r$se[[i]] <- abs(linkfunctions7::dlinkinv(link, eta)) * se_eta
          ends <- sort(c(linkfunctions7::linkinv(link, eta - z * se_eta),
                         linkfunctions7::linkinv(link, eta + z * se_eta)))
          r$lower[[i]] <- ends[[1L]]
          r$upper[[i]] <- ends[[2L]]
        }
      }
      stats::setNames(r, c(cols, "source"))
    })
    do.call(rbind, out)
  }
  # A term whose block is a working linearization is reported by what it is
  # ABOUT and not by the coefficients it is fitted through: a break-point
  # term's auxiliary pair carries the position as -g/delta, which is a
  # number no reader wants to compute. term_readable() gives the quantities
  # and the Jacobian from the coefficients, so the variance comes across by
  # the delta method -- the same route segmented reports a break-point's
  # standard error by. A term that answers NULL, which is what a developed
  # one does, falls back to its coefficients.
  readable_rows <- function(nm) {
    term <- spec@terms[[p]][[nm]]
    idx <- design[[p]]$blocks[[nm]]
    rd <- tryCatch(modelterms7::term_readable(term, fit@coefficients[[p]][idx]),
                   error = function(e) NULL)
    if (is.null(rd) || !length(rd$name)) return(NULL)
    # V is indexed by NAME rather than by position: coef_labels() skips a
    # parameter with no coefficients, so a stacked offset computed here
    # would not be the one the matrix was built with
    key <- paste(p, design[[p]]$coef_names[idx], sep = ":")
    se <- rep(NA_real_, length(rd$name))
    if (!is.null(V) && all(key %in% rownames(V))) {
      Vb <- as.matrix(V[key, key, drop = FALSE])
      if (all(is.finite(Vb))) {
        se <- sqrt(pmax(diag(rd$jacobian %*% Vb %*% t(rd$jacobian)), 0))
      }
    }
    z <- stats::qnorm(1 - (1 - level) / 2)
    st <- rd$value / se
    # A break-point gets an estimate and an interval and NO test. The null
    # a z of value/se would report on is that the position is zero, which
    # is not a hypothesis anyone holds; the one a reader wants is that
    # there is no break-point at all, and under it the position is a
    # nuisance parameter that vanishes, so the classical p-value is wrong
    # by a factor of three to five (Davies' problem). segmented prints the
    # estimate and the standard error alone for the same reason.
    pos <- grepl("^psi[0-9]*$", rd$name)
    st[pos] <- NA_real_
    out <- data.frame(name = rd$name, estimate = rd$value, se = se,
                      statistic = st, p_value = 2 * stats::pnorm(-abs(st)),
                      lower = rd$value - z * se, upper = rd$value + z * se,
                      role = "coefficient", source = "",
                      stringsAsFactors = FALSE)
    stats::setNames(out, all_cols)
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
    body <- if (identical(kind, "breakpoint")) {
      rr <- readable_rows(nm)
      if (is.null(rr)) cr else rr
    } else {
      cr[keep, , drop = FALSE]
    }
    # the hyperparameters come FIRST in every penalized block: they govern
    # everything below them, and a table that opens with a hundred selected
    # coefficients buries the one number that produced that selection
    tb <- rbind(hyper_rows(nm), body)
    blocks[[length(blocks) + 1L]] <- list(
      kind = kind,
      label = switch(kind, smooth = "Smooth", random = "Random effect",
                     selection = "Selection", breakpoint = "Break-points",
                     "Penalized"),
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

    # a structural term of this equation: no columns, so no block above
    st <- x@structural
    if (!is.null(st) && any(st$parameter == p)) {
      for (tn in unique(st$term[st$parameter == p])) {
        r <- st[st$parameter == p & st$term == tn, , drop = FALSE]
        cat(sprintf("\n  %s   (structural: no design columns)\n", tn))
        tb <- data.frame(estimate = r$estimate, se = r$se,
                         lower = r$lower, upper = r$upper,
                         row.names = ifelse(r$held, paste0(r$name, " (held)"),
                                            r$name))
        print(format(tb, digits = digits))
      }
    }
  }

  cat(sprintf("\n%.0f%% intervals, %s variance\n", 100 * x@level, x@type))
  # WHICH log-likelihood and WHICH degrees of freedom, because the pairing
  # is what makes the criterion mean anything and the two conventions in
  # common use are not comparable. This one is conditional: the likelihood
  # is read at the fitted coefficients, a penalized coefficient among them,
  # and the count is the effective degrees of freedom. A mixed-model package
  # reporting a MARGINAL likelihood integrates its random effects out and
  # counts variance components instead, and its AIC is a different number
  # answering a different question (Vaida and Blanchard, 2005).
  cat(sprintf(paste0("conditional log-likelihood %.6f    effective df %.2f",
                     "\ncAIC %.3f    cBIC %.3f\n"),
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
  # a hyperparameter row prints numbers where there are any: one estimated by
  # a marginal criterion carries a standard error and an interval. Where there
  # is none the columns are blanked and the cell where the standard error
  # would have been says what put the value there instead
  hyp <- tb$role %in% c("fixed", "estimated")
  fixed <- hyp & !is.finite(tb$se)
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
  src <- if (is.null(tb$source)) tb$role else
    ifelse(nzchar(tb$source), tb$source, tb$role)
  out$se[fixed] <- paste0("(", src[fixed], ")")
  out[fixed, c("z", "p", "lower", "upper")] <- ""
  # an estimated one that does carry an interval still has no test: the null
  # a z would report on is that the hyperparameter is zero, the edge of its
  # range rather than an interior hypothesis
  out[hyp & !fixed, c("z", "p")] <- ""
  rownames(out) <- tb$name
  print(out)
  invisible(NULL)
}
