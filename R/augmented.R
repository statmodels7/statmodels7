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
#' \eqn{K} is the number of distribution parameters, so it is 1, 2 or 3 for
#' most families and never large; \eqn{n} is the number of observations and
#' is. The standard Cholesky recursion is therefore written out over the
#' \eqn{K} indices and evaluated over all \eqn{n} observations at once, one
#' vectorized pass per entry of the factor, so no loop runs over
#' observations.
#'
#' A block that is not positive definite produces a non-positive pivot, and
#' the function returns `NULL` at that point instead of taking its square
#' root. That is an ordinary outcome and not a failure: the observed
#' curvature far from the optimum is routinely indefinite, and the caller
#' answers by falling back to the assembled route. A warning about a `NaN`
#' would report an expected branch as a defect.
#'
#' @param Om An \eqn{n \times K \times K} array of symmetric blocks, as
#'   [obs_information()] returns it. Only the lower triangle of each block is
#'   read.
#'
#' @return An \eqn{n \times K \times K} numeric array, lower triangular in
#'   its last two indices, with `Om[i, , ] == L[i, , ] %*% t(L[i, , ])` for
#'   every `i`. `NULL` as soon as any block fails to be positive definite,
#'   so a single bad observation declines for the whole sample.
#'
#' @seealso [obs_information()] for the input,
#'   [sqrt_design()] for what the factors are used to build.
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
#' \eqn{i} with respect to the \eqn{K} link-scale predictors, multiplied by
#' that observation's prior weight. This is the per-observation curvature a
#' scoring step is built from, before any design enters.
#'
#' @details
#' The derivatives come from the family on the link scale, so the chain rule
#' onto \eqn{\eta} has already been applied by \pkg{distributions7} and
#' nothing here multiplies by a link's derivative. With `expected = TRUE` the
#' expected information is taken and the blocks are positive definite
#' wherever the family is regular; with `expected = FALSE` the observed
#' Hessian is negated, which far from the optimum is routinely indefinite.
#' That is what makes [chol_blocks()]'s refusal an ordinary branch.
#'
#' The weights enter as given, without normalization, so a weight of two
#' counts an observation twice.
#'
#' @param spec A [StatmodSpec()], read for the distribution, the weights and
#'   the thread count.
#' @param theta The per-observation parameters on the parameter scale, a
#'   named list of \eqn{n}-vectors, one per distribution parameter.
#' @param expected `TRUE` for the expected information, `FALSE` for the
#'   negated observed Hessian.
#' @param approx How the expected information is approximated for a family
#'   that has no closed form, passed through to \pkg{distributions7}:
#'   `"bartlett"`, `"integrate"` or `"mc"`. Ignored when `expected` is
#'   `FALSE` or when the family computes its expected information exactly.
#'
#' @return An \eqn{n \times K \times K} numeric array, symmetric in its last
#'   two indices, with \eqn{K} the number of distribution parameters.
#'
#' @seealso [chol_blocks()], which factorizes these,
#'   [statmod_information_at()] for the assembled \eqn{Z'\Omega Z}.
#'
#' @keywords internal
info_blocks <- function(spec, theta, expected = TRUE, approx = "bartlett") {
  params <- spec@distrib@params
  K <- length(params)
  n <- spec@n_obs
  H <- if (expected) {
    distributions7::distrib_expected_hessian(spec@distrib, spec@response,
                                             theta, scale = "link",
                                             approx = approx, threads = spec@threads)
  } else {
    distributions7::distrib_hessian(spec@distrib, spec@response, theta,
                                    scale = "link", threads = spec@threads)
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
#' Builds \eqn{R} with \eqn{R'R = Z'\Omega Z}, the square root of the
#' assembled information. A scoring step decomposes this instead of the
#' information itself, so \eqn{Z'\Omega Z} is never formed and its condition
#' number is never squared.
#'
#' @details
#' Row block \eqn{a} of \eqn{R} carries, in the columns belonging to
#' parameter \eqn{b}, the entry \eqn{L_i[b, a]} times that parameter's design
#' row. The factor is lower triangular, so the blocks with \eqn{b < a} are
#' zero and are not formed at all.
#'
#' The result is sparse when any equation's design block is, which is what
#' sends the solve to [sparse_augmented_solve()].
#'
#' @param design The design, as [statmod_design()] returns it, holding one
#'   block of columns per distribution parameter.
#' @param L The per-observation Cholesky factors, as [chol_blocks()] returns
#'   them, or `NULL`.
#'
#' @return An \eqn{nK \times p} matrix, dense or a \pkg{Matrix} sparse
#'   matrix according to the design's own storage, where \eqn{p} is the total
#'   number of coefficients across the equations. `NULL` when `L` is `NULL`,
#'   so the refusal propagates from the factorization to the solve without a
#'   test at each step.
#'
#' @seealso [chol_blocks()] for `L`,
#'   [augmented_solve()] for the decomposition this feeds,
#'   [penalty_sqrt()] for the other half of the augmented matrix.
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
#' # Why an eigendecomposition and not a Cholesky
#'
#' A penalty is positive semidefinite and is often rank deficient. A spline's
#' is deficient by exactly the dimension of its null space, which is the
#' whole point of the construction: the null space is what the penalty does
#' not shrink. A Cholesky fails there, so the factor comes from an
#' eigendecomposition with the non-positive eigenvalues dropped, and the rows
#' they would have contributed are simply absent.
#'
#' A non-convex penalty has an indefinite Hessian and no such factor at all.
#' `NULL` is returned and the caller falls back to the assembled route.
#'
#' # The diagonal shortcut
#'
#' A diagonal penalty is factored by taking the square root of its diagonal.
#' That is the same answer the eigendecomposition gives, since the
#' eigenvalues of a diagonal matrix are its diagonal entries, so the two
#' routes agree by construction and a test pins them together.
#'
#' The case is worth detecting for how often it arises. A ridge is diagonal,
#' a random effect is diagonal, and the Demmler-Reinsch penalty of
#' [modelterms7::s()] is \eqn{\mathrm{diag}(0, 1, \ldots, 1)} exactly. So is
#' any block-diagonal assembly of them. The factor is recomputed at every
#' iteration of the scoring loop, so its cost is multiplied by the iteration
#' count: measured on a random intercept over 1000 groups, one dense
#' eigendecomposition of the 1003 by 1003 penalty cost 0.63 s and was 83 per
#' cent of the whole fit.
#'
#' @param S The penalty Hessian, a `p x p` symmetric matrix, dense or sparse.
#'   `Matrix::isDiagonal()` decides which route is taken, so a matrix stored
#'   as diagonal and one merely having zero off-diagonals are both caught.
#'
#' @return A matrix with `ncol(S)` columns and one row per retained
#'   eigendirection, so at most `p` and fewer where the penalty has a null
#'   space. Its class mirrors `S`'s, dense for dense and sparse for sparse.
#'   `NULL` when the penalty is indefinite beyond the tolerance.
#'
#' @seealso [penalty_sqrt_diag()] for the diagonal route,
#'   [augmented_solve()] for the solve this feeds,
#'   [sqrt_design()] for the other half of the augmented matrix.
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
#' that coordinate's entry, which is [penalty_sqrt()]'s answer
#' where the matrix is diagonal.
#'
#' @details
#' The thresholds are the eigen route's, applied to the diagonal, which for
#' a diagonal matrix is its spectrum. A negative entry beyond the tolerance
#' makes the penalty indefinite and there is no factor to return; an entry
#' at or below the tolerance is a null direction and contributes no row.
#'
#' The class of the result mirrors the argument's and is not chosen here.
#' [augmented_solve()] routes on whether either of its two factors is
#' sparse, so a sparse factor returned for a dense design would send a dense
#' fit down the sparse route.
#'
#' @param S A diagonal penalty Hessian, `p x p`, dense or sparse. Only its
#'   diagonal is read.
#' @param p Its dimension, passed rather than read so the caller's own count
#'   is used.
#'
#' @return A matrix with `p` columns and one row per retained coordinate,
#'   each row zero except for the square root of that coordinate's entry.
#'   Its class mirrors `S`'s. `NULL` when any entry is negative beyond the
#'   tolerance.
#'
#' @seealso [penalty_sqrt()], which dispatches here for a diagonal penalty.
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
#' Decomposes the augmented matrix \eqn{A = [R;\, C]}, whose cross-product
#' \eqn{A'A = R'R + C'C} is the penalized information, and returns the
#' increment solving
#'
#' \deqn{(R'R + C'C)\,\delta = u}
#'
#' The augmented form is the point: neither \eqn{R'R} nor the penalized
#' information is ever formed, so a design whose condition number is
#' \eqn{\kappa} is decomposed at \eqn{\kappa} and not at \eqn{\kappa^2}.
#'
#' @details
#' # The four routes
#'
#' Tried in this order.
#'
#' 1. **Sparse QR**, when either factor is a \pkg{Matrix} object. Handed to
#'    [sparse_augmented_solve()], which declines on a rank-deficient matrix.
#' 2. **SVD**, when `how = "svd"`. The singular values below
#'    \eqn{\max(\dim A)\,\epsilon\,\sigma_{\max}} are dropped and the
#'    increment is built from the retained right singular vectors, so a
#'    deficient system gets the minimum-norm answer.
#' 3. **The threaded triangular factor**, when `how = "qr"`, `threads > 1`
#'    and the work \eqn{n p^2} is at least `5e7`. Only the triangular factor
#'    is ever read, so a kernel that produces it and never accumulates
#'    \eqn{Q} does the whole job.
#' 4. **`qr()`**, the pivoted LINPACK route, which reports a rank and can
#'    drop columns. A sparse solve that declined falls through to here on
#'    densified factors.
#'
#' # Why the rank test is equilibrated
#'
#' The threaded route is taken only where the matrix is comfortably of full
#' rank, and the test reads the diagonal of the triangular factor **divided
#' by each column's norm**. Since \eqn{A'A = R'R}, scaling \eqn{A}'s columns
#' scales that diagonal by the same factors, so the ratio is the diagonal of
#' a decomposition with unit column norms.
#'
#' Reading the raw diagonal instead reports a matrix as near-singular
#' whenever its columns differ in size, which is what a large smoothing
#' parameter does to the penalty rows of its own block beside an unpenalized
#' one. Per-column scaling forgives separation from any source; an exact
#' collinearity stays exactly singular.
#'
#' The tolerance is `1e-7` on that ratio, which is what `dqrdc2` itself uses,
#' so anything the pivoted route would call deficient still goes there and is
#' reported with its rank.
#'
#' Measured against `qr()` at \eqn{n = 40000} with the column scales spread
#' over \eqn{10^8}: 2.6x at \eqn{p = 51} and 3.8x at \eqn{p = 145} to 600 on
#' eight threads, with the increment agreeing to \eqn{1.6 \times 10^{-16}}.
#'
#' @param R The square-root design, as [sqrt_design()] returns it,
#'   \eqn{nK \times p}.
#' @param C The penalty's factor, as [penalty_sqrt()] returns it, with `p`
#'   columns and at most `p` rows. May have zero rows for an unpenalized
#'   model.
#' @param u The right-hand side, a numeric vector of length `p`.
#' @param how `"qr"` or `"svd"`. `"svd"` has no sparse counterpart, so a
#'   sparse pair asked for it is densified first.
#' @param threads How many threads the triangular factor may use, a plain
#'   integer. `1L` takes `qr()` unconditionally.
#'
#' @return A list of two:
#'   \describe{
#'     \item{`delta`}{the increment, an unnamed numeric vector of length `p`.
#'       Zero in any coordinate the pivoted route dropped.}
#'     \item{`rank`}{the rank used, an integer. Equal to `p` on the threaded
#'       and sparse routes, which decline rather than report a deficiency.}
#'   }
#'
#' @seealso [sqrt_design()] and [penalty_sqrt()] for the two factors,
#'   [sparse_augmented_solve()] for the sparse route.
#'
#' @keywords internal
augmented_solve <- function(R, C, u, how, threads = 1L) {
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
  # THE THREADED FACTOR. Only the triangular factor is ever read here, so a
  # kernel that produces it and never accumulates Q does the whole job: the
  # trailing columns of each Householder step are independent and each is
  # written in full by one thread, so the factor is bit-identical at any
  # count. It is engaged above a measured gate and only where the matrix is
  # comfortably of full rank on the JACOBI-EQUILIBRATED diagonal, at the
  # tolerance dqrdc2 itself uses, so anything the pivoted route would call
  # deficient still goes there and is reported with its rank. Measured
  # against qr() at n = 40000 with column scales spread over 1e8: 2.6x at
  # p = 51 and 3.8x at p = 145 to 600 on eight threads, with the increment
  # agreeing to 1.6e-16.
  fast <- NULL
  if (threads > 1L && how == "qr" &&
      as.double(nrow(A)) * ncol(A) * ncol(A) >= 5e7) {
    cn <- sqrt(colSums(A * A))
    if (all(is.finite(cn)) && all(cn > 0)) {
      Rt <- qr_factor_cpp(A, threads)
      dg <- abs(diag(Rt))
      eq <- dg / cn
      if (all(is.finite(eq)) && max(eq) > 0 && min(eq) > 1e-7 * max(eq)) {
        d <- tryCatch(backsolve(Rt, forwardsolve(t(Rt), u)),
                      error = function(e) NULL)
        if (!is.null(d) && all(is.finite(d))) {
          fast <- list(delta = as.numeric(d), rank = ncol(A))
        }
      }
    }
  }
  if (!is.null(fast)) return(fast)

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
#' The same increment [augmented_solve()] returns, taken through a
#' sparse QR of \eqn{[R;\ C]}.
#'
#' @details
#' # Why a QR and not a Cholesky
#'
#' The augmented form exists so that \eqn{X'X} is never formed and the
#' conditioning is never squared. A sparse QR is a QR and keeps that
#' property exactly. A sparse Cholesky of the normal equations would be
#' faster still and would give it up.
#'
#' Measured against the dense QR on the same augmented design of a
#' random-intercept model: 695 times faster at 100 groups and 75475 times at
#' 1000, where the dense factorization costs 9.06 s against 0.00012 s.
#'
#' # The factor is not back-permuted
#'
#' A sparse QR reorders the columns to reduce fill, so \eqn{AP = QR} for the
#' permutation \eqn{P} the decomposition chose, and
#' \eqn{(A'A)^{-1} = P(R'R)^{-1}P'}. The increment is two triangular solves
#' between a permutation and its inverse, which is the bookkeeping the dense
#' route already does with `qr()`'s pivot.
#'
#' `qrR(backPermute = TRUE)` looks like the simpler choice and is a trap. The
#' back-permuted factor is no longer triangular, so its diagonal says nothing
#' about the rank, and a solve against it is a general solve instead of two
#' triangular ones.
#'
#' # The rank test
#'
#' Read off the diagonal of the triangular factor, since a sparse QR does not
#' report a rank. The test is on the **Jacobi-equilibrated** diagonal, each
#' entry divided by its column's norm, and declines when the smallest ratio
#' falls to `ncol(A) * .Machine$double.eps` of the largest.
#'
#' The equilibration is what makes the route usable. Reading the raw diagonal
#' rejected 87 of 127 solves on matrices the dense route finds at full rank,
#' the ratio there running down to \eqn{7 \times 10^{-30}} while the
#' equilibrated one stayed between 0.445 and 1. Those matrices are not
#' near-singular; their columns differ in size, which is what a large
#' smoothing parameter does. With the raw test a random-effect fit fell
#' through to a dense QR and cost 34.1 s where it now costs 3.9.
#'
#' The column norms come from `colSums(A^2)` and not `colSums(A * A)`. A
#' binary operation between two sparse matrices intersects their index sets,
#' which is 29 to 35 times slower here, and the norms were 70 to 74 per cent
#' of the whole solve before the change.
#'
#' Where the matrix really is rank deficient there is no unique increment,
#' and this route declines rather than choosing one. The caller falls back to
#' the dense route, which drops columns and says how many it kept.
#'
#' @param R The square-root design, \eqn{nK \times p}, sparse or dense.
#' @param C The penalty's factor, with `p` columns.
#' @param u The right-hand side, a numeric vector of length `p`.
#' @param how The decomposition asked for. Only `"qr"` is served; `"svd"` has
#'   no sparse counterpart and declines.
#'
#' @return A list with `delta` (numeric, length `p`) and `rank` (integer,
#'   equal to `p`), or `NULL` when the route declines: `how` is not `"qr"`,
#'   the factorization fails, or the equilibrated diagonal says the matrix is
#'   deficient.
#'
#' @seealso [augmented_solve()], the caller and the dense fallback.
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
  # the columns in the order the decomposition put them
  q <- qrA@q
  ord <- if (length(q)) q + 1L else seq_len(ncol(A))
  # The rank is read off the JACOBI-EQUILIBRATED diagonal, which is the same
  # correction solve_pd() took in 0.70.0 and which was not propagated here.
  # Since R'R = A'A, scaling A's columns by their norms scales that diagonal
  # by the same factors, so |R_jj| / ||a_j|| is the diagonal of the
  # decomposition of a matrix with unit column norms: per-direction scaling
  # forgives separation from any source, while an exact collinearity stays
  # exactly singular. Without it a penalized augmented matrix was rejected
  # for having columns of different SIZE -- a large smoothing parameter
  # makes the penalty rows of its own block enormous beside an unpenalized
  # one -- and every such solve fell through to a dense QR of the whole
  # thing. Measured on `s(x) + random(~1|g)`, 87 of 127 solves were rejected
  # and the dense route found full rank in all 127; equilibrated, the ratio
  # runs from 0.445 to 1.000 where the raw one reached 7.4e-30.
  # `A^2` and NOT `A * A`: the second is a binary operation between two
  # sparse matrices and intersects their index sets, where the first acts
  # on the `x` slot with the pattern untouched.  Measured on the three
  # shapes a penalized fit meets here, same values to the bit: 50.00 ms
  # against 1.74, 26.56 against 0.76, 19.69 against 0.64 -- 29x to 35x,
  # and the norms were 70 to 74 per cent of this solve.
  cn <- sqrt(Matrix::colSums(A^2))[ord]
  if (length(cn) != length(dg) || !all(is.finite(cn)) || any(cn == 0)) {
    return(NULL)
  }
  eq <- dg / cn
  if (!all(is.finite(eq)) || max(eq) == 0) return(NULL)
  if (min(eq) <= ncol(A) * .Machine$double.eps * max(eq)) return(NULL)
  d <- tryCatch({
    z <- Matrix::solve(Matrix::t(Rf), u[ord])
    as.numeric(Matrix::solve(Rf, z))
  }, error = function(e) NULL)
  if (is.null(d) || !all(is.finite(d))) return(NULL)
  delta <- numeric(ncol(A))
  delta[ord] <- d
  list(delta = delta, rank = ncol(A))
}
