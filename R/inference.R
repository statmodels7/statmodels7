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
    # A structural term's own information is assembled dense, so the padded
    # penalty is dense too rather than a sum of two kinds whose result would
    # depend on which side of the bind it arrived on.
    S <- as_dense(S)
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
  S <- zap_nonfinite(S)

  out <- matrix(NA_real_, total, total, dimnames = list(nm, nm))
  if (!any(keep)) return(out)
  keep_full <- c(keep, rep(TRUE, nz))
  A <- (H + S)[keep_full, keep_full, drop = FALSE]
  # what a SMALL eigenvalue means is decided inside solve_pd(), on the
  # equilibrated matrix: a smoothing parameter at 1e15 and a break-point
  # term's annealed working columns both separate the scales without
  # flattening any direction, and per-direction scaling forgives both
  Vb <- solve_pd(A, "the penalized information",
                 c(nm[keep], rep("", nz)))
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
#'
#' @return The inverse.
#'
#' @keywords internal
solve_pd <- function(A, what, labels = NULL) {
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
  # The test runs on the JACOBI-EQUILIBRATED matrix, D^-1/2 A D^-1/2 with
  # D its diagonal, which is what tells a flat direction from scale
  # separation whatever produced the separation. A reference scale passed
  # by the caller (the unpenalized information's largest diagonal) used to
  # do that, and it covered one source only: a smoothing parameter at 1e15
  # separates the scales and the reference forgave it, but the committed
  # working block of a break-point term carries auxiliary columns near
  # 1/(2 c d) -- 1e8 at the annealed floor, 1e16 on the information's
  # diagonal -- so the DESIGN itself set the reference and ordinary
  # curvature at 242 in the ordinary-scale directions read as flat
  # (measured, on the flagship jseg's own summary). Equilibration makes
  # the test per-direction: unit diagonal, so the smallest eigenvalue is
  # small only where a direction is flat AGAINST ITS OWN SCALE, and an
  # exact collinearity stays exactly singular. The inverse is recovered
  # through the same scaling, D^-1/2 (D^-1/2 A D^-1/2)^-1 D^-1/2.
  #
  # A non-positive diagonal entry means that coordinate carries no
  # information at all, which is the flat case stated directly.
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
  if (ncol(A) > 0L && all(is.finite(A)) && all(diag(A) > 0)) {
    s <- 1 / sqrt(diag(A))
    Ae <- A * tcrossprod(s)
    ch <- tryCatch(chol(Ae), error = function(e) NULL)
    if (!is.null(ch)) {
      # the 1-norm, which is what the condition estimator is expressed in
      anorm <- max(colSums(abs(Ae)))
      rc <- chol_rcond_cpp(ch, anorm)
      lmin <- if (is.na(rc)) 0 else rc * anorm
      if (is.finite(lmin) && lmin > 1e-12 * anorm) {
        # Inverted through the factor already in hand. A condition number of
        # 1e15 born of scale separation costs the well-determined directions
        # nothing here and the shrunk ones simply report variances near zero.
        out <- chol2inv(ch) * tcrossprod(s)
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


#' Is a Matrix Worth Factorizing Sparsely?
#'
#' @description
#' Whether a symmetric matrix is large enough and sparse enough that a sparse
#' Cholesky beats a dense one, asked of the MATRIX and of nothing else.
#'
#' @details
#' The two conditions are the measured crossover and not a preference. On the
#' penalized information of a random intercept over \eqn{m} levels at 20000
#' observations, the sparse route against the dense one -- coercion,
#' factorization, log-determinant and full inverse, each timed with the
#' repetition loop sized by elapsed time:
#'
#' \tabular{rrrrr}{
#'   \strong{m} \tab \strong{p} \tab \strong{density} \tab \strong{whole route}
#'     \tab \strong{inverse} \cr
#'   20 \tab 23 \tab 0.282 \tab 1.08x \tab \strong{0.13x} \cr
#'   50 \tab 53 \tab 0.128 \tab 1.06x \tab \strong{0.33x} \cr
#'   100 \tab 103 \tab 0.067 \tab 1.18x \tab 1.14x \cr
#'   200 \tab 203 \tab 0.034 \tab 1.74x \tab 2.9x \cr
#'   500 \tab 503 \tab 0.014 \tab 5.28x \tab 7.6x \cr
#'   1000 \tab 1003 \tab 0.007 \tab 2.50x \tab 11.2x \cr
#'   2000 \tab 2003 \tab 0.003 \tab 4.20x \tab 7.0x
#' }
#'
#' The whole-route column builds the matrix afresh on each repetition and is
#' the one the thresholds are read from. The inverse column is the like-for-
#' like comparison, each route carrying its own factorization; an earlier
#' version of it timed the sparse solves against a factor built once outside
#' the loop, which flattered the sparse side without changing where it loses.
#'
#' Below about a hundred coefficients the sparse route LOSES, and loses badly:
#' its fixed cost is the coercion and the S4 dispatch around it, which does not
#' shrink with the matrix. On the fully dense penalized information of a single
#' smooth (p = 16, density 1) it measures 0.01x, a hundred times slower, which
#' is what the size condition is there to prevent.
#'
#' \strong{Both quantities are read off the matrix, and the first one is its
#' STORAGE.} A matrix held as a base matrix is refused whatever its zeros,
#' which reads like a test of the container rather than of the mathematics, so
#' it is worth saying why it is neither an oversight nor a term test.
#' \code{\link{statmod_information_at}} accumulates into the design's own kind,
#' so the penalized matrix is stored sparsely exactly when the design is, and
#' \pkg{modelterms7} builds a block sparse only when asked
#' (\code{sparse = TRUE}, whose default is \code{FALSE}). Measured on
#' \code{y ~ 0 + g + s(x)} over 400 levels at 20000 observations, whose
#' penalized matrix is 5 per cent nonzero either way: built dense the fit takes
#' 104.24 s and this factorization is \strong{0.16 per cent} of it, the time
#' being in the \eqn{O(np^2)} products against a dense design
#' (\code{statmod_information_at} 48.8 per cent, \code{crossprod} 57.2 per cent
#' of self time); built sparse the same fit takes 2.19 s. So where the storage
#' is dense the factorization is not what a fit is spending its time on, and
#' coercing a dense \eqn{p \times p} matrix here to save a share of that size
#' would cost more than it returns.
#'
#' \strong{The like-for-like comparison is the one that says this is not a term
#' test}, and it is the check \code{piano_lme4.txt} section 5 asks for. With
#' every design built the same way, this route is worth 1.38x on
#' \code{0 + g + s(x)} over 400 levels, 1.33x on \code{random(~1|g)} over 500
#' and 1.07x on \code{s(x, by = g)} over 60 -- an unpenalized indicator block,
#' a random effect and a factor-\code{by} smooth, gaining together and in the
#' order their sizes predict. Nothing here asks which term or which family
#' produced the matrix.
#'
#' @param M A symmetric matrix.
#' @param min_dim The smallest order worth the fixed cost.
#' @param max_density The largest fraction of nonzeros worth it.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{pd_factor}}
#'
#' @keywords internal
worth_sparse <- function(M, min_dim = 100L, max_density = 0.10) {
  if (!isS4(M)) return(FALSE)
  p <- ncol(M)
  if (!length(p) || p < min_dim) return(FALSE)
  nz <- tryCatch(length(M@x), error = function(e) NA_integer_)
  if (is.na(nz)) return(FALSE)
  nz / (as.numeric(p) * p) <= max_density
}


#' The Smallest Eigenvalue of a Sparse Factor's Matrix, Estimated
#'
#' @description
#' \eqn{1/\lVert A^{-1}\rVert_1} from a sparse Cholesky factor, which is the
#' quantity LAPACK's \code{dpocon} produces from a dense one.
#'
#' @details
#' The sparse route needs a condition estimate OF ITS OWN, and it cannot
#' borrow the dense one: \code{chol_rcond_cpp} reads a dense triangular
#' factor. \code{Matrix::rcond} is not the answer either -- measured, it costs
#' 10.3 ms at p = 503 and 500 ms at p = 2003, more than the dense
#' factorization the sparse route exists to replace. Higham's one-norm
#' estimator applied to the factor's own solves costs 0.58 ms at p = 53 and
#' 0.80 ms at p = 1003, nearly flat, because it is a handful of triangular
#' solves and an R loop around them.
#'
#' For a symmetric matrix \eqn{\lVert A^{-1}\rVert_1 \ge \lVert
#' A^{-1}\rVert_2 = 1/\lambda_{\min}}, so the quantity returned is at or below
#' the smallest eigenvalue, and the estimator's own error is a further
#' underestimate of the norm in the other direction. It is used exactly as the
#' dense estimate is: to separate a matrix comfortably positive definite from
#' one that is not, two situations that differ by some fifteen orders of
#' magnitude here (measured on a design with two identical columns, 1.4e-14
#' relative to the matrix's scale, against 1.5e3 for a hyperparameter driven
#' to 1e15). A factor of two either way cannot move that verdict, which is the
#' argument already recorded for the dense estimator.
#'
#' @param L A \code{CHMfactor}.
#' @param p The order of the matrix.
#'
#' @return A single number, or \code{NA_real_} where the estimate failed.
#'
#' @seealso \code{\link{pd_factor}}, \code{\link{pd_logdet}}
#'
#' @keywords internal
sparse_lmin <- function(L, p) {
  ax <- function(x) as.matrix(Matrix::solve(L, as.matrix(x)))
  est <- tryCatch(
    Matrix::onenormest(A.x = ax, At.x = ax, n = p, t = 1, silent = TRUE)$est,
    error = function(e) NA_real_)
  if (!is.finite(est) || est <= 0) return(NA_real_)
  1 / est
}


#' Factorize a Penalized Information Once
#'
#' @description
#' The Cholesky factor of a matrix a Laplace approximation needs positive
#' definite, together with its log-determinant, in whichever storage the
#' matrix itself calls for.
#'
#' @details
#' This is the one place the penalized matrix is factorized. The criterion
#' wants its log-determinant, the gradient wants the mode's movement and the
#' Hessian wants both plus the inverse; before this existed the criterion and
#' \code{\link{ctx_penalized}} each factorized the SAME matrix at the same
#' point, which at p = 503 was 12.4 ms spent twice.
#'
#' \strong{The verdict is unchanged and so is its property.} Whether the
#' matrix is accepted never turns on whether a factorization raised: where the
#' cheap test is inconclusive the eigendecomposition answers about the matrix.
#' The sparse route carries its own condition estimate
#' (\code{\link{sparse_lmin}}) rather than the dense one, and falls back to
#' the dense route where that estimate cannot be formed, so a refusal is
#' reached by the same reasoning on either storage.
#'
#' @param M A symmetric matrix, sparse or dense.
#' @param scale A reference magnitude, as \code{\link{pd_logdet}} takes.
#'
#' @return A list with \code{logdet}, \code{ok}, \code{factor} and
#'   \code{sparse}. The factor is \code{NULL} where the answer came from the
#'   eigendecomposition.
#'
#' @seealso \code{\link{pd_logdet}}, \code{\link{ctx_penalized}}
#'
#' @keywords internal
pd_factor <- function(M, scale = NULL) {
  no <- function() list(logdet = NA_real_, ok = FALSE, factor = NULL,
                        sparse = FALSE)
  if (!length(ncol(M)) || !ncol(M)) return(no())

  if (worth_sparse(M)) {
    p <- ncol(M)
    if (!all(is.finite(M@x))) return(no())
    anorm <- max(Matrix::colSums(abs(M)))
    ref <- if (is.null(scale) || !is.finite(scale) || scale <= 0) anorm
           else min(scale, anorm)
    Ms <- tryCatch(Matrix::forceSymmetric(M), error = function(e) NULL)
    # CHOLMOD reports a matrix it cannot factorize as a WARNING and not as an
    # error, where the dense route raises and is caught. A hyperparameter the
    # search should step away from is not something to warn a caller about --
    # the criterion returns NULL there and the search backtracks -- so the
    # warning is suppressed and the verdict is read off the return value, as
    # it is on the dense route.
    L <- if (is.null(Ms)) NULL else
      suppressWarnings(tryCatch(Matrix::Cholesky(Ms, LDL = FALSE, super = NA),
                                error = function(e) NULL))
    if (!is.null(L)) {
      lmin <- sparse_lmin(L, p)
      if (is.finite(lmin) && lmin > sqrt(.Machine$double.eps) * ref) {
        ld <- tryCatch(2 * as.numeric(Matrix::determinant(
          L, logarithm = TRUE, sqrt = TRUE)$modulus),
          error = function(e) NA_real_)
        if (is.finite(ld)) {
          return(list(logdet = ld, ok = TRUE, factor = L, sparse = TRUE))
        }
      }
    }
    # inconclusive: the dense route asks the matrix, which is what settles it
  }

  d <- pd_logdet_dense(as_dense(M), scale)
  c(d, list(sparse = FALSE))
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
  r <- pd_factor(M, scale)
  r$factor <- NULL
  r$sparse <- NULL
  r
}


#' The Dense Route of pd_logdet
#'
#' @description
#' The three routes described at \code{\link{pd_logdet}}, on a dense matrix.
#'
#' @details
#' Split out so that \code{\link{pd_factor}} can reach it as the fallback of
#' the sparse route without restating the verdict: there is one place that
#' decides whether a matrix is positive definite, and one set of thresholds.
#'
#' @param M A dense symmetric matrix.
#' @param scale A reference magnitude.
#'
#' @return A list with \code{logdet}, \code{ok} and, on a refusal reached
#'   through the eigendecomposition, \code{min_ev} and \code{max_ev}.
#'
#' @seealso \code{\link{pd_logdet}}
#'
#' @keywords internal
pd_logdet_dense <- function(M, scale = NULL) {
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
      # the factor travels with the answer so that a caller wanting the
      # inverse as well does not compute a second one of the same matrix
      return(list(logdet = 2 * sum(log(diag(ch))), ok = TRUE, factor = ch))
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
    notes = S7::class_character,
    # the reading of statmod_certificate(), or NULL where it could not be
    # taken. A list rather than columns: its four readings are of different
    # kinds and one of them is a set of names.
    certificate = S7::new_property(S7::class_any, default = NULL)
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

  # A smoothed break-point term: the smoother and its width are part of the
  # model -- the transition's width, the bent-cable reading -- so they are
  # reported; and where a break-point carries a random development, a
  # smoother declaring a scale correction (the probit's convolution
  # identity) gets the corrected scale printed beside the apparent one.
  notes <- c(notes, tryCatch(smoothed_notes(spec, object),
                             error = function(e) character(0)))

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
    level = level, type = type, notes = notes,
    # read HERE and not in statmod(): it costs one outer gradient and one
    # solve, which is nothing beside a summary that already inverts the
    # penalized information, and a fold of cv() or a path point that never
    # prints itself pays nothing at all.
    certificate = tryCatch(statmod_certificate(object),
                           error = function(e) NULL))
}
S7::method(summary, StatmodFit) <- summary.StatmodFit


#' The Quantities a Penalty's Hyperparameters Are About
#'
#' @description
#' Replaces the coordinate rows of a penalty whose hyperparameters are a chart
#' with the quantities it declares through
#' \code{\link[penalties7]{penalty_readable}}: the standard deviations and
#' correlations of a correlated random effect, rather than the logarithms of a
#' Cholesky diagonal and the entries below it.
#'
#' @details
#' The standard error is the delta method, and it composes two Jacobians: the
#' penalty's, which is in the parameter scale of its hyperparameters, and the
#' link's, the variance matrix being on the free scale the outer criterion was
#' maximized on. Each interval is built on the scale the quantity declares --
#' log for a standard deviation, Fisher's z for a correlation -- and mapped
#' back, so a standard deviation cannot be given a negative lower end and a
#' correlation cannot be given an interval that leaves \eqn{(-1, 1)}. That is
#' the rule every other interval in the toolkit follows.
#'
#' No test is printed, for the reason the coordinate rows print none: the null
#' a \eqn{z} of value over standard error reports on is that the quantity is
#' zero, which for a standard deviation is the edge of its range.
#'
#' @param rd The result of \code{\link[penalties7]{penalty_readable}}.
#' @param th The penalty's hyperparameters, as fitted.
#' @param Vh The hyperparameter variance matrix, or \code{NULL}.
#' @param p The parameter the term sits in.
#' @param key The penalty's key.
#' @param level The confidence level.
#' @param role,src What the coordinate rows reported.
#' @param cols The column names of a summary block.
#'
#' @return A data frame of rows, in the shape of a summary block.
#'
#' @keywords internal
readable_hyper_rows <- function(rd, th, Vh, p, key, level, role, src, cols) {
  nm <- names(th)
  k <- length(rd$value)
  se <- rep(NA_real_, k)
  if (!is.null(Vh)) {
    j <- match(paste(p, key, nm, sep = "\r"), rownames(Vh))
    if (!anyNA(j)) {
      vb <- as.matrix(Vh[j, j, drop = FALSE])
      lk <- attr(attr(Vh, "idx"), "links")[j]
      d <- vapply(seq_along(nm), function(i) {
        eta <- linkfunctions7::linkfun(lk[[i]], th[[i]])
        linkfunctions7::dlinkinv(lk[[i]], eta)
      }, 0)
      jac <- rd$jacobian * rep(d, each = k)
      if (all(is.finite(vb)) && all(is.finite(jac))) {
        se <- sqrt(pmax(diag(jac %*% vb %*% t(jac)), 0))
      }
    }
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  v <- as.numeric(rd$value)
  tr <- as.character(rd$transform)
  lo <- hi <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    if (!is.finite(se[i])) next
    ends <- switch(tr[i],
      log = exp(log(v[i]) + c(-1, 1) * z * se[i] / v[i]),
      atanh = tanh(atanh(v[i]) + c(-1, 1) * z * se[i] / (1 - v[i]^2)),
      v[i] + c(-1, 1) * z * se[i])
    ends <- sort(ends)
    lo[i] <- ends[[1L]]
    hi[i] <- ends[[2L]]
  }
  out <- data.frame(name = names(rd$value), estimate = v, se = se,
                    statistic = NA_real_, p_value = NA_real_,
                    lower = lo, upper = hi,
                    role = rep(role[[1L]], k), source = rep(src[[1L]], k),
                    stringsAsFactors = FALSE)
  stats::setNames(out, c(cols, "source"))
}

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
  # WHAT WAS HELD is the terms' answer and nobody else's: a hyperparameter
  # a term fixed is fixed whatever criterion ran beside it.
  held <- tryCatch(statmod_held(spec, design), error = function(e) character(0))
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
      chosen <- (rep(marginal, length(hk)) | hk %in% spc_keys) &
        !(hk %in% held)
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
      # WHERE THE HYPERPARAMETERS ARE A CHART and not the quantities, the
      # quantities are reported instead: a correlated random effect is about
      # the standard deviations and correlations of its effects, and the
      # logarithms of a Cholesky diagonal are the coordinates that produce
      # them. The penalty declares which case it is; every other branch
      # answers NULL and its rows stand.
      rd <- tryCatch(penalties7::penalty_readable(u$penalty, th),
                     error = function(e) NULL)
      if (!is.null(rd) && length(rd$value)) {
        rr <- readable_hyper_rows(rd, th, if (marginal) Vh else NULL, p, key,
                                  level, role, src, cols)
        # A hyperparameter the readable block does not DESCRIBE keeps its own
        # row: a multivariate Student t is about the standard deviations and
        # the correlations of its scale matrix, and its degrees of freedom are
        # none of those. The question is asked of the Jacobian -- a column that
        # is zero throughout is a coordinate no quantity depends on -- so a
        # family that declares more later is covered without an edit.
        keep <- apply(rd$jacobian, 2L, function(z) all(z == 0))
        if (any(keep)) {
          rr <- rbind(rr, stats::setNames(r[keep, , drop = FALSE],
                                          c(cols, "source")))
        }
        return(rr)
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


#' The Notes a Smoothed Break-Point Term Adds to a Summary
#'
#' @description
#' One note per smoothed term, naming the smoother and the width the build
#' resolved -- the width of the transition, which is part of the model and
#' not a detail -- and, where a break-point carries a random development
#' under a Gaussian precision and the smoother declares a scale correction,
#' the corrected scale beside the apparent one.
#'
#' @details
#' The correction is the smoother's own: the probit satisfies the exact
#' convolution identity \eqn{\tau^2_{\mathrm{apparent}} = \tau^2 + h^2}, so
#' the corrected scale is \eqn{\sqrt{\tau^2 - h^2}}; a smoother declaring
#' none (the hyperbolic, the quintic) gets the apparent scale alone, which
#' the random effect's own block already reports. The apparent scale is
#' read off the ridge precision as \eqn{1/\sqrt{\lambda}}, which is only a
#' scale where the penalty is the quadratic branch with that
#' hyperparameter; any other development is left without the note rather
#' than given a number of the wrong meaning.
#'
#' @param spec The fitted specification, whose terms are the ones the fit
#'   left.
#' @param object The fit.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso \code{\link[penalties7]{abs_smoother}}
#'
#' @keywords internal
smoothed_notes <- function(spec, object) {
  out <- character(0)
  for (p in names(spec@terms)) {
    for (nm in names(spec@terms[[p]])) {
      tm <- spec@terms[[p]][[nm]]
      if (!S7::S7_inherits(tm, modelterms7::SegTerm)) next
      bp <- tryCatch(tm@blueprint, error = function(e) NULL)
      smx <- if (is.list(bp)) bp$smooth else NULL
      if (is.null(smx)) next
      sm <- smx$sm
      out <- c(out, sprintf(paste0(
        "'%s' in '%s' is smoothed (%s, %s = %s%s): its break-points are ",
        "ordinary\n  parameters and that is the width of the transition, ",
        "so its rows are the\n  smoothed model's."),
        nm, p, sm@smoother_name, sm@width_name,
        format(smx$width, digits = 3),
        if (!is.null(smx$w_group)) ", per group" else ""))
      if (is.null(sm@tau_correction)) next
      ent <- tryCatch(modelterms7::term_penalties(tm),
                      error = function(e) list())
      for (e in ent) {
        if (!grepl("^psi", e$name)) next
        pen <- e$penalty
        key <- statmod_entry_key(nm, ent, e)
        th <- tryCatch(object@hyper[[p]][[key]], error = function(err) NULL)
        if (is.null(th)) next
        # The apparent scale of the break-point deviations. A gaussian
        # prior written by its scale carries it as sigma, which is what
        # random() declares; the quadratic ridge carries the precision, so
        # the scale is 1/sqrt(lambda). Any other prior is left without the
        # note rather than given a number of the wrong meaning: the
        # convolution identity composes GAUSSIAN variances.
        tau <- NULL
        if ("sigma" %in% pen@params &&
            grepl("gaussian", pen@penalty_name, fixed = TRUE)) {
          tau <- suppressWarnings(as.numeric(th[["sigma"]]))
        } else if (identical(pen@params, "lambda") &&
                   isTRUE(tryCatch(penalties7::is_quadratic(pen),
                                   error = function(err) FALSE))) {
          lam <- suppressWarnings(as.numeric(th[["lambda"]]))
          if (length(lam) == 1L && is.finite(lam) && lam > 0) {
            tau <- 1 / sqrt(lam)
          }
        }
        if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
            tau <= 0) next
        tauc <- sm@tau_correction(tau, smx$width)
        out <- c(out, sprintf(paste0(
          "The scale of the random break-points of '%s' composes with the ",
          "smoothing\n  width: apparent tau %s, corrected ",
          "sqrt(tau^2 - %s^2) = %s."),
          nm, format(tau, digits = 4), sm@width_name,
          format(tauc, digits = 4)))
      }
    }
  }
  out
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
  # THE CERTIFICATE, which is a property of the POINT where the line above is
  # a property of the search. The two disagree exactly where it matters: a
  # search can stop on its own rule far from an optimum, and a search that
  # kept going can end at one. See statmod_certificate().
  if (!is.null(x@certificate)) {
    ct <- x@certificate
    cat(sprintf("certificate: %s", toupper(ct$state)))
    if (is.finite(ct$gradient)) {
      cat(sprintf("   outer gradient %.3g", ct$gradient))
    }
    if (is.finite(ct$mode_error)) {
      cat(sprintf("   %.3g above the mode", ct$mode_error))
    }
    cat("\n")
    if (length(ct$boundary)) {
      cat("  at a boundary: ", paste(ct$boundary, collapse = ", "), "\n",
          sep = "")
    }
    for (r in ct$reason) cat("  ", r, "\n", sep = "")
  }
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


#' What the Fit Certifies About the Point It Reports
#'
#' @description
#' Three readings taken AT THE REPORTED POINT and independent of the path the
#' search took: the outer criterion's gradient, how far the coefficients sit
#' above the penalized mode, and which hyperparameters have run to a boundary.
#'
#' @details
#' \strong{Why a certificate rather than the optimizer's flag.} The flag says
#' whether a search stopped on its own rule, which is a statement about the
#' search. Measured across shapes, it does not order fits by quality: on one
#' model the default reported success at a criterion of -1783.47 while the same
#' data under \code{\link[optimizers7]{lbfgs}} reached -1664.43 and reported
#' failure. What a reader wants is a property of the point.
#'
#' \strong{The state comes from the gradient and the mode error is reported
#' beside it, not folded into it.} Measured at the reported point over six
#' shapes, the outer gradient separates by five orders -- 4.7e-07, 7.8e-07,
#' 5.8e-05, 7.7e-05 and 3.0e-04 on fits that are right, against 28.8 on one
#' that is not -- while the mode error does not: it reads 1.8e-16 to 6.1e-12 on
#' four of them, 22.8 on the failing one, and 0.114 on a random-changepoint
#' \code{seg} whose answer is right to a correlation of 0.9932. A number that
#' does not separate cannot decide a state, and a certificate that says how far
#' from the mode is worth more than a boolean that hides it.
#'
#' \code{tol} is 1e-2 rather than the geometric middle of the two groups: the
#' two ways of being wrong are not symmetric, and a certificate that says NOT
#' CONVERGED at a good point is visible and checkable where one that certifies
#' a bad point is the failure this exists to remove.
#'
#' \strong{What it costs} is one outer gradient and one solve, once, at a point
#' the fit already holds. Nothing is refitted: measured, the criterion
#' reconstructed from \code{fit@spec} equals the one the fit reports EXACTLY on
#' every shape, so the reading is of the fitted model and not of another one.
#'
#' \strong{Where there is no outer gradient there are two cases, and they get
#' different answers.} A model with NO PENALTY -- \code{linpar}, \code{nl},
#' \code{seg}, \code{jump}, \code{jseg} -- has no hyperparameter for a
#' gradient to be about, so the only question left is whether the inner fit
#' reached its mode, and the mode error answers it: measured over the
#' reference battery it reads 5.2e-11 to 7.9e-05 on fits that are right
#' against 1.215 on a \code{jump} fitted to data carrying a slope and a slope
#' change it has no term for. A model whose only hyperparameters are KINKED --
#' \code{lasso}, \code{scad}, \code{mcp}, swept along a path because a Laplace
#' approximation at a mode sitting on the kink has no meaning -- gets neither
#' reading and stays \code{"unknown"}: at a coefficient the penalty has set to
#' zero the score does not vanish but lies in the subdifferential, so the mode
#' error is not a statement about being at a mode. Measured on a lasso, its
#' 4.7e-03 is carried by a coordinate whose coefficient is exactly 0 and whose
#' score is -0.715.
#'
#' A form whose criterion has no EXACT gradient
#' (\code{\link{outer_gradient_ok}}) is also \code{"unknown"} rather than
#' approximated: 2p refits to difference it would cost more than the fit.
#'
#' \strong{The boundary label, and why its threshold needs no derivation.}
#' A hyperparameter may run to an edge and belong there: on a covariate that
#' is pure noise the smoothing parameter reaches 9.2e+08, the criterion is
#' genuinely flat, and calling that fit unconverged would be wrong. A
#' coordinate is reported as sitting at a boundary when its free value
#' exceeds \code{edge} AND its own gradient component has already met
#' \code{tol}. Because of that second condition the threshold cannot change
#' the verdict: a coordinate it moves out of the interior set had already
#' passed the test, so the maximum that decides the state is unaffected, and
#' both \code{"converged"} and \code{"boundary"} are certified. What
#' \code{edge} decides is how the point is described. The default separates
#' the measured cases with room on both sides: coordinates that ran to an
#' edge sit at 9.3, 10.5 and 20.6 on the free scale against 0.13, 0.30 and
#' 2.01 for the ones that did not.
#'
#' @param fit A \code{\link{StatmodFit}}.
#' @param tol The largest outer gradient a certified point may carry.
#' @param edge The free value beyond which a hyperparameter whose gradient
#'   has already met \code{tol} is reported as sitting at a boundary. It
#'   decides the label alone and never the verdict; see the details.
#'
#' @return A list with \code{state} (\code{"converged"}, \code{"boundary"},
#'   \code{"not converged"} or \code{"unknown"}), \code{gradient},
#'   \code{mode_error}, \code{boundary} and \code{reason}.
#'
#' @seealso \code{\link{statmod}}, \code{\link{mode_error_limit}},
#'   \code{\link{criterion_resolution}}
#'
#' @examples
#' dd <- data.frame(x = runif(120))
#' dd$y <- sin(4 * dd$x) + rnorm(120, 0, 0.3)
#' statmod_certificate(statmod(y ~ s(x, k = 8),
#'                             distributions7::gaussian1_distrib(), dd))$state
#'
#' @export
statmod_certificate <- function(fit, tol = 1e-2, edge = 8) {
  out <- list(state = "unknown", gradient = NA_real_, mode_error = NA_real_,
              boundary = character(0), reason = character(0))
  method <- fit@methods$outer
  spec <- fit@spec
  design <- tryCatch(statmod_design(spec), error = function(e) NULL)
  if (is.null(design)) {
    out$reason <- "the design could not be rebuilt from the fit"
    return(out)
  }
  cf <- fit@coefficients
  hy <- fit@hyper
  ctx <- tryCatch(outer_context(spec, design, cf, hy, "bartlett"),
                  error = function(e) NULL)

  # HOW FAR ABOVE THE MODE, in log-likelihood units: the decrease the mode's
  # own Newton correction predicts. Reported whether or not a criterion ran.
  if (!is.null(ctx)) {
    pen <- tryCatch(ctx_penalized(ctx, spec, design, cf, hy, FALSE),
                    error = function(e) NULL)
    if (!is.null(pen)) {
      obj <- tryCatch(statmod_objective(spec, hy, design, FALSE, "bartlett"),
                      error = function(e) NULL)
      sc <- if (is.null(obj)) NULL else
        tryCatch(obj$gr(obj$stack(cf)), error = function(e) NULL)
      if (!is.null(sc) && all(is.finite(sc))) {
        db <- tryCatch(as.numeric(as.matrix(pen$inv) %*% sc),
                       error = function(e) NULL)
        if (!is.null(db) && all(is.finite(db))) {
          out$mode_error <- 0.5 * sum(sc * db)
        }
      }
    }
  }

  blocks <- tryCatch(statmod_blocks(spec, design), error = function(e) NULL)
  idx <- if (is.null(blocks)) NULL else outer_hyper_index(spec, blocks)
  no_outer <- is.null(method) || !method@kind %in% c("ml", "reml") ||
    is.null(idx) || !nrow(idx)
  if (no_outer) {
    # ⚠️ TWO QUITE DIFFERENT WAYS TO HAVE NO OUTER GRADIENT, and they call for
    # different answers. Measured over the reference battery, seven cases of
    # twenty-nine land here:
    #
    #   NO PENALTY AT ALL -- linpar, nl, seg, jump, jseg sharp and smoothed.
    #     There is no hyperparameter to estimate, so there is nothing an outer
    #     gradient could say, and the only question left is whether the inner
    #     fit reached its mode. The mode error answers exactly that, and it
    #     answers well: 5.2e-11, 1.5e-10, 1.9e-10, 2.0e-06 and 7.9e-05 on fits
    #     that are right, against 1.215 on a `jump` fitted to data carrying a
    #     slope and a slope change it has no term for -- misspecified, so its
    #     break-point iteration never settles, its annealing runs to the floor
    #     and its block reaches 3.0e+13 on the information's diagonal. On data
    #     where the jump IS the truth the same term reads 4.2e-04.
    #
    #   A KINKED HYPERPARAMETER -- lasso, scad, mcp -- chosen by a PATH
    #     because a Laplace approximation at a mode sitting on the kink is
    #     arithmetic without a meaning. There is a hyperparameter, but it is
    #     the argmin over a grid rather than the root of a derivative. AND THE
    #     MODE ERROR IS NOT A READING HERE EITHER: at a coefficient the
    #     penalty has set to zero the score does not vanish, it lies in the
    #     subdifferential. Measured on a lasso, the mode error of 4.7e-03 is
    #     carried by a coordinate whose coefficient is exactly 0 and whose
    #     score is -0.715. So this stays unknown, and says why.
    pen_units <- tryCatch(statmod_penalized(spec, design),
                          error = function(e) list())
    if (length(pen_units)) {
      out$reason <- paste0(
        "the only hyperparameters here belong to a penalty with a kink, ",
        "chosen along a path rather than by a criterion with a derivative; ",
        "and at a coefficient the penalty has set to zero the score does not ",
        "vanish, so the mode error is not a reading either")
      return(out)
    }
    if (!is.finite(out$mode_error)) {
      out$reason <- "the model carries no penalty, and the mode error could not be read"
      return(out)
    }
    out$state <- if (out$mode_error <= mode_error_limit()) "converged" else
      "not converged"
    out$reason <- sprintf(paste0(
      "the model carries no penalty, so there is no outer gradient; the ",
      "reading is the inner fit's own, %.4g log-likelihood units above its ",
      "mode against %g"), out$mode_error, mode_error_limit())
    return(out)
  }
  if (!outer_gradient_ok(spec, design, idx, method, 1L)) {
    out$reason <- "this form has no exact outer gradient, and differencing it would cost more than the fit"
    return(out)
  }
  basis <- integrated_basis(spec, design, method@kind)
  g <- tryCatch(statmod_marginal_grad(spec, design, cf, hy, method, idx, basis,
                                      ctx = ctx),
                error = function(e) NULL)
  if (is.null(g) || !all(is.finite(g))) {
    # ⚠️ A NON-FINITE OUTER GRADIENT AT A BOUNDARY IS NOT THE SAME COMPLAINT,
    # and saying which it is costs one read. Where a coefficient sits at the
    # clamp its link keeps it strictly inside, the family's THIRD and FOURTH
    # derivatives there are not finite -- measured on the Student t, all ten
    # components of the third go from nu = 1e150, which is sqrt(double.xmax)
    # and so the signature of a product of two quantities of order nu -- and
    # the outer gradient reads exactly those. The fit itself may be perfectly
    # well located: on the reference battery's `fam-studentt` the mode error
    # is 1.3e-10 against a limit of 1e-3 and the criterion is the best of the
    # three routes, so what is missing is the reading and not the answer.
    #
    # The state stays `not converged` rather than `boundary`, deliberately:
    # certifying a point whose hyperparameters were never verified would
    # claim more than has been checked. What it can honestly do is name the
    # coordinate, so the reader is not left with "not finite".
    H0 <- tryCatch(ctx_information(ctx, spec, design, cf, hy, FALSE,
                                   "bartlett"), error = function(e) NULL)
    frozen <- if (is.null(H0)) integer(0) else boundary_coords(H0)
    if (length(frozen)) {
      npar <- vapply(design, function(d) d$npar, integer(1))
      ends <- cumsum(npar)
      who <- unique(names(npar)[vapply(frozen, function(j)
        which(j <= ends)[1], integer(1))])
      out$boundary <- who
      out$reason <- sprintf(paste0(
        "the outer gradient is not finite at the reported point, because a ",
        "coefficient of %s sits at the bound its link keeps it strictly ",
        "inside and the family's third and fourth derivatives are not finite ",
        "there. The fit itself is %.4g log-likelihood units above its mode"),
        paste(who, collapse = ", "), out$mode_error)
    } else {
      out$reason <- "the outer gradient is not finite at the reported point"
    }
    out$state <- "not converged"
    return(out)
  }
  # A COORDINATE AT A BOUNDARY is one whose criterion has stopped moving in it
  # while its value has run far from where it started. Both halves are needed:
  # a small gradient component alone is what convergence looks like, and a
  # large value alone is an ordinary answer on a wide scale.
  #
  # `edge` IS NOT DERIVED FROM ANYTHING, and it does not have to be, because
  # THE CONJUNCTION BELOW IS WHAT MAKES IT SAFE: a coordinate is called an
  # edge only if it has ALREADY met `tol`, so removing it from `interior`
  # cannot raise the maximum that decides the verdict. Move the threshold in
  # either direction and the only thing that changes is whether the state
  # reads `converged` or `boundary`, both of which are certified. Delete the
  # `abs(g) <= tol` conjunct, however, and the threshold starts excusing
  # coordinates from the gradient test, at which point this paragraph is
  # false and the number needs an argument of its own.
  #
  # What the default separates, measured: the coordinates that ran to an edge
  # sit at |eta| of 9.3, 10.5 and 20.6 -- a smoothing parameter of 9.2e+08 on
  # pure noise, prior scales of 9.2e-05 and 2.8e-05 -- against 0.13, 0.30 and
  # 2.01 for the ones that did not, so 8 sits in a wide gap rather than on a
  # boundary of its own.
  eta <- hyper_to_eta(hy, idx)
  at_edge <- which(abs(eta) > edge & abs(g) <= tol)
  if (length(at_edge)) {
    out$boundary <- vapply(at_edge, function(k)
      paste(idx$parameter[k], idx$term[k], idx$name[k], sep = "/"),
      character(1))
  }
  interior <- setdiff(seq_along(g), at_edge)
  out$gradient <- if (length(interior)) max(abs(g[interior])) else 0
  if (out$gradient <= tol) {
    out$state <- if (length(at_edge)) "boundary" else "converged"
  } else {
    out$state <- "not converged"
    out$reason <- sprintf(
      "the outer criterion's gradient is %.4g at the reported point, against %g",
      out$gradient, tol)
  }
  if (is.finite(out$mode_error) && out$mode_error > mode_error_limit()) {
    out$reason <- c(out$reason, sprintf(
      paste0("the coefficients sit %.4g log-likelihood units above the",
             " penalized mode, so everything read there -- the criterion,",
             " its gradient, vcov() -- carries that"),
      out$mode_error))
  }
  out
}
