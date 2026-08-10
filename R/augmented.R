#' @include objective.R
NULL

# The square-root design, so that the scoring step never forms X'X.
#
# With more than one modelled parameter the information is not X'WX for a
# diagonal W: it is Z' Omega Z, where Z is the block diagonal of the equations'
# designs and Omega is block diagonal BY OBSERVATION with a K by K block per
# observation, K being the number of modelled parameters. Those blocks are
# information matrices of one observation and are positive definite wherever
# the family is regular, so each has a Cholesky factor, and
#
#     Omega_i = L_i L_i'   =>   Z' Omega Z = R'R,  R = [L_1' Z_1; ...]
#
# which is the matrix to decompose. Building it costs one small factorization
# per observation, written out over the K by K indices and vectorized over the
# observations, so no loop runs over n.

#' Cholesky Factors of Small Blocks, Vectorized Over Observations
#'
#' @description
#' Returns the lower triangular \eqn{L_i} with \eqn{L_iL_i' = \Omega_i} for
#' every observation, as a \eqn{n \times K \times K} array.
#'
#' @details
#' The standard recursion is written over the \eqn{K} indices, which are few,
#' and evaluated over all observations at once, which are many.
#'
#' A block that is not positive definite has a non-positive pivot, and the
#' function returns \code{NULL} at that point rather than taking its square
#' root: the observed curvature far from the optimum is routinely indefinite,
#' so this is an ordinary outcome the caller answers by falling back to the
#' assembled route, and a warning about a \code{NaN} would report it as a
#' defect.
#'
#' @param Om An \eqn{n \times K \times K} array of symmetric blocks.
#'
#' @return An array of the same shape, lower triangular in its last two
#'   indices, or \code{NULL} when some block is not positive definite.
#'
#' @keywords internal
chol_blocks <- function(Om) {
  n <- dim(Om)[1L]
  K <- dim(Om)[2L]
  L <- array(0, dim = c(n, K, K))
  for (j in seq_len(K)) {
    s <- Om[, j, j]
    if (j > 1L) {
      for (m in seq_len(j - 1L)) s <- s - L[, j, m]^2
    }
    if (anyNA(s) || any(s <= 0)) return(NULL)
    L[, j, j] <- sqrt(s)
    if (j < K) {
      for (i in seq.int(j + 1L, K)) {
        v <- Om[, i, j]
        if (j > 1L) {
          for (m in seq_len(j - 1L)) v <- v - L[, i, m] * L[, j, m]
        }
        L[, i, j] <- v / L[, j, j]
      }
    }
  }
  L
}


#' The Per-Observation Information Blocks
#'
#' @description
#' Assembles \eqn{\Omega_i}, the \eqn{K \times K} information of observation
#' \eqn{i} in the link-scale predictors, weighted by the prior weight.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param theta The per-observation parameters.
#' @param expected Whether the expected information is wanted.
#' @param approx The approximation, where the expected one is not closed.
#'
#' @return An \eqn{n \times K \times K} array.
#'
#' @keywords internal
info_blocks <- function(spec, theta, expected = TRUE, approx = "bartlett") {
  params <- spec@distrib@params
  K <- length(params)
  n <- spec@n_obs
  H <- if (expected) {
    distributions7::distrib_expected_hessian(spec@distrib, spec@response,
                                             theta, scale = "link",
                                             approx = approx)
  } else {
    distributions7::distrib_hessian(spec@distrib, spec@response, theta,
                                    scale = "link")
  }
  Om <- array(0, dim = c(n, K, K))
  for (a in seq_len(K)) {
    for (b in seq.int(a, K)) {
      v <- -spec@weights * rep_len(H[[hess_key(params, a, b)]], n)
      Om[, a, b] <- v
      Om[, b, a] <- v
    }
  }
  Om
}


#' The Square-Root Design
#'
#' @description
#' Returns \eqn{R} with \eqn{R'R = Z'\Omega Z}, the matrix a scoring step
#' decomposes instead of the information itself.
#'
#' @details
#' Row block \eqn{a} of \eqn{R} carries, in the columns of parameter \eqn{b},
#' the entry \eqn{L_i[b, a]} times that parameter's design row; the factor
#' being lower triangular, the blocks with \eqn{b < a} are zero and are not
#' formed.
#'
#' @param design The design, as \code{\link{statmod_design}} returns it.
#' @param L The Cholesky factors, from \code{\link{chol_blocks}}.
#'
#' @return An \eqn{nK \times p} matrix, or \code{NULL} when a block was not
#'   positive definite.
#'
#' @keywords internal
sqrt_design <- function(design, L) {
  if (is.null(L) || anyNA(L) || any(!is.finite(L))) return(NULL)
  K <- dim(L)[2L]
  Xs <- lapply(design, function(d) d$X)
  npar <- vapply(design, function(d) d$npar, integer(1))
  n <- dim(L)[1L]
  rows <- vector("list", K)
  for (a in seq_len(K)) {
    blocks <- vector("list", K)
    for (b in seq_len(K)) {
      if (npar[b] == 0L) {
        blocks[[b]] <- matrix(0, n, 0L)
      } else if (b < a) {
        blocks[[b]] <- matrix(0, n, npar[b])
      } else {
        blocks[[b]] <- L[, b, a] * Xs[[b]]
      }
    }
    rows[[a]] <- do.call(cbind, blocks)
  }
  do.call(rbind, rows)
}


#' A Square-Root Factor of the Penalty
#'
#' @description
#' Returns \eqn{C} with \eqn{C'C = S(\theta)}, the block-diagonal penalty
#' Hessian, dropping the rows a null space contributes nothing to.
#'
#' @details
#' A penalty is positive semidefinite and may be rank deficient -- a spline's
#' is, by exactly the dimension of its null space -- so the factor comes from
#' an eigendecomposition with the non-positive eigenvalues dropped, rather
#' than from a Cholesky, which would fail there. A non-convex penalty, whose
#' Hessian is indefinite, has no such factor and the caller falls back.
#'
#' @param S The penalty Hessian.
#'
#' @return A matrix with \code{ncol(S)} columns, or \code{NULL}.
#'
#' @keywords internal
penalty_sqrt <- function(S) {
  p <- ncol(S)
  if (p == 0L) return(matrix(0, 0L, 0L))
  if (all(S == 0)) return(matrix(0, 0L, p))
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  tol <- p * .Machine$double.eps * max(abs(e$values))
  if (min(e$values) < -max(tol, 1e-10 * max(abs(e$values)))) return(NULL)
  keep <- e$values > tol
  if (!any(keep)) return(matrix(0, 0L, p))
  sqrt(e$values[keep]) * t(e$vectors[, keep, drop = FALSE])
}


#' Solve a Scoring Step From the Square-Root Design
#'
#' @description
#' Decomposes the augmented matrix \eqn{[R;\ C]}, whose cross-product is the
#' penalized information, and returns the increment solving
#' \eqn{(R'R + C'C)\delta = u}.
#'
#' @param R The square-root design.
#' @param C The penalty's factor.
#' @param u The right-hand side.
#' @param how Either \code{"qr"} or \code{"svd"}.
#'
#' @return A list with \code{delta} and \code{rank}.
#'
#' @keywords internal
augmented_solve <- function(R, C, u, how) {
  A <- rbind(R, C)
  if (how == "svd") {
    s <- svd(A, nu = 0L)
    tol <- max(dim(A)) * .Machine$double.eps * max(s$d)
    keep <- s$d > tol
    dinv <- numeric(length(s$d))
    dinv[keep] <- 1 / s$d[keep]^2
    return(list(delta = as.numeric(s$v %*% (dinv * crossprod(s$v, u))),
                rank = sum(keep)))
  }
  qrA <- qr(A)
  # (A'A)^-1 u from the R factor alone: A'A = R_qr' R_qr
  Rf <- qr.R(qrA)[seq_len(qrA$rank), seq_len(qrA$rank), drop = FALSE]
  piv <- qrA$pivot[seq_len(qrA$rank)]
  z <- forwardsolve(t(Rf), u[piv])
  d <- backsolve(Rf, z)
  delta <- numeric(ncol(A))
  delta[piv] <- d
  list(delta = delta, rank = qrA$rank)
}
