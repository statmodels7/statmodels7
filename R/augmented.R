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
  # The zero blocks are the larger part of this matrix -- one equation's rows
  # carry nothing of the equations before it -- and materializing them dense
  # was allocating n by npar of nothing per pair, per iteration. They follow
  # the design's own kind, so a sparse block stays sparse all the way into
  # the decomposition instead of pulling the whole augmented matrix dense.
  sp <- any(vapply(Xs, isS4, logical(1)))
  zero <- function(nc) {
    if (nc == 0L || !sp) return(matrix(0, n, nc))
    Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                         dims = c(n, nc))
  }
  rows <- vector("list", K)
  for (a in seq_len(K)) {
    blocks <- vector("list", K)
    for (b in seq_len(K)) {
      if (npar[b] == 0L) {
        blocks[[b]] <- zero(0L)
      } else if (b < a) {
        blocks[[b]] <- zero(npar[b])
      } else {
        blocks[[b]] <- L[, b, a] * Xs[[b]]
      }
    }
    rows[[a]] <- if (sp) Reduce(function(x, y) Matrix::cbind2(x, y), blocks)
      else do.call(cbind, blocks)
  }
  if (sp) Reduce(function(x, y) Matrix::rbind2(x, y), rows)
  else do.call(rbind, rows)
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
#' A DIAGONAL penalty is factored by taking the square root of its diagonal,
#' which is the same answer the eigendecomposition returns -- the eigenvalues
#' of a diagonal matrix are its diagonal, so the two routes agree by
#' construction and a test pins them together. It is not a special case worth
#' having for its own sake but for how often it is the one that arises: a
#' ridge, a random effect and the Demmler-Reinsch penalty of \code{s()},
#' which is \eqn{\mathrm{diag}(0, 1, \ldots, 1)} exactly, are all diagonal,
#' and so is any block-diagonal assembly of them. The factor is recomputed at
#' every iteration of the scoring loop, so the cost is the decomposition's
#' times the iteration count: measured on a random intercept over 1000 groups,
#' one dense eigendecomposition of the 1003 by 1003 penalty costs 0.63 s and
#' was 83 per cent of the whole fit.
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
  if (isTRUE(tryCatch(Matrix::isDiagonal(S), error = function(e) FALSE))) {
    return(diagonal_sqrt(S, p))
  }
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  tol <- p * .Machine$double.eps * max(abs(e$values))
  if (min(e$values) < -max(tol, 1e-10 * max(abs(e$values)))) return(NULL)
  keep <- e$values > tol
  if (!any(keep)) return(matrix(0, 0L, p))
  sqrt(e$values[keep]) * t(e$vectors[, keep, drop = FALSE])
}

#' The Factor of a Diagonal Penalty
#'
#' @description
#' One row per coordinate the penalty reaches, carrying the square root of
#' that coordinate's entry, which is \code{\link{penalty_sqrt}}'s answer
#' where the matrix is diagonal.
#'
#' @details
#' The thresholds are the eigen route's, read on the diagonal, which for a
#' diagonal matrix IS its spectrum: a negative entry beyond the tolerance
#' makes the penalty indefinite and there is no factor to return, and an
#' entry at the tolerance is a null direction and contributes no row. The
#' class of the result mirrors the argument's rather than being chosen:
#' \code{\link{augmented_solve}} routes on whether either of the two factors
#' is sparse, so returning a sparse factor for a dense design would send a
#' dense fit through the sparse route and a dense one through neither.
#'
#' @param S A diagonal penalty Hessian.
#' @param p Its dimension.
#'
#' @return A matrix with \code{p} columns, or \code{NULL}.
#'
#' @keywords internal
diagonal_sqrt <- function(S, p) {
  d <- as.numeric(Matrix::diag(S))
  mx <- max(abs(d))
  tol <- p * .Machine$double.eps * mx
  if (min(d) < -max(tol, 1e-10 * mx)) return(NULL)
  keep <- which(d > tol)
  if (!length(keep)) return(matrix(0, 0L, p))
  if (isS4(S)) {
    return(Matrix::sparseMatrix(i = seq_along(keep), j = keep,
                                x = sqrt(d[keep]),
                                dims = c(length(keep), p)))
  }
  C <- matrix(0, length(keep), p)
  C[cbind(seq_along(keep), keep)] <- sqrt(d[keep])
  C
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
  if (isS4(R) || isS4(C)) {
    out <- sparse_augmented_solve(R, C, u, how)
    if (!is.null(out)) return(out)
    # a rank-deficient or unfactorable augmented matrix falls through to the
    # dense route, which reports a rank and can drop columns
    R <- as_dense(R)
    C <- as_dense(C)
  }
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


#' Solve a Scoring Step From a Sparse Square-Root Design
#'
#' @description
#' The same increment \code{\link{augmented_solve}} returns, taken through a
#' sparse QR of \eqn{[R;\ C]}.
#'
#' @details
#' The augmented route exists so that \eqn{X'X} is never formed and the
#' conditioning is never squared, and a sparse QR is a QR: it keeps that
#' property exactly, which a sparse Cholesky of the normal equations would
#' not. Measured against the dense QR on the same augmented design of a
#' random-intercept model, it is 695 times faster at 100 groups and 75475 at
#' 1000, the dense factorization there costing 9.06 s against 0.00012.
#'
#' The factor is taken WITHOUT back-permuting. A sparse QR reorders the
#' columns to reduce fill, so \eqn{AP = QR} for the permutation \eqn{P} the
#' decomposition chose, and \eqn{(A'A)^{-1} = P(R'R)^{-1}P'}: the increment
#' is two triangular solves between a permutation and its inverse, which is
#' the same bookkeeping the dense route does with \code{qr}'s pivot.
#' \code{backPermute = TRUE} looks simpler and is a trap -- it returns a
#' factor that is no longer triangular, so its diagonal says nothing about
#' the rank and a solve against it is a general one rather than two
#' triangular ones.
#'
#' The rank is read off the diagonal of the triangular factor rather than
#' reported by the decomposition, which for a sparse QR does not give one.
#' Where the matrix is rank deficient there is no unique increment to return
#' and this route declines, leaving the caller its dense fallback, which can
#' drop columns and say how many it kept.
#'
#' @param R The square-root design.
#' @param C The penalty's factor.
#' @param u The right-hand side.
#' @param how The decomposition asked for; \code{"svd"} has no sparse
#'   counterpart and declines.
#'
#' @return A list with \code{delta} and \code{rank}, or \code{NULL}.
#'
#' @seealso \code{\link{augmented_solve}}
#'
#' @keywords internal
sparse_augmented_solve <- function(R, C, u, how) {
  if (!identical(how, "qr")) return(NULL)
  A <- tryCatch(Matrix::rbind2(as_sparse(R), as_sparse(C)),
                error = function(e) NULL)
  if (is.null(A) || nrow(A) < ncol(A)) return(NULL)
  qrA <- tryCatch(Matrix::qr(A), error = function(e) NULL)
  if (is.null(qrA)) return(NULL)
  Rf <- tryCatch(Matrix::qrR(qrA, backPermute = FALSE),
                 error = function(e) NULL)
  if (is.null(Rf) || nrow(Rf) != ncol(A)) return(NULL)
  dg <- abs(Matrix::diag(Rf))
  if (!length(dg) || !all(is.finite(dg)) || max(dg) == 0) return(NULL)
  if (min(dg) <= ncol(A) * .Machine$double.eps * max(dg)) return(NULL)
  # the columns in the order the decomposition put them
  q <- qrA@q
  ord <- if (length(q)) q + 1L else seq_len(ncol(A))
  d <- tryCatch({
    z <- Matrix::solve(Matrix::t(Rf), u[ord])
    as.numeric(Matrix::solve(Rf, z))
  }, error = function(e) NULL)
  if (is.null(d) || !all(is.finite(d))) return(NULL)
  delta <- numeric(ncol(A))
  delta[ord] <- d
  list(delta = delta, rank = ncol(A))
}
