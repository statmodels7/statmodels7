#' @include objective.R
NULL

# The square-root design, so that the scoring step never forms X'X.
#
# With more than one modeled parameter the information is not X'WX for a
# diagonal W: it is Z' Omega Z, where Z is the block diagonal of the equations'
# designs and Omega is block diagonal BY OBSERVATION with a K by K block per
# observation, K being the number of modeled parameters. Those blocks are
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
#' \eqn{K} is the number of distribution parameters, which is small for every
#' family, whereas \eqn{n}, the number of observations, can be large. The
#' standard Cholesky recursion is therefore written out over the
#' \eqn{K} indices and evaluated over all \eqn{n} observations at once, one
#' vectorized pass per entry of the factor, so no loop runs over
#' observations.
#'
#' A block that is not positive definite produces a non-positive pivot, and
#' the function returns `NULL` at that point instead of taking its square
#' root. This outcome is expected, because the observed curvature far from the
#' optimum is routinely indefinite, and the caller then falls back to the
#' assembled route.
#'
#' @param Om An \eqn{n \times K \times K} array of symmetric blocks, as
#'   [info_blocks()] returns it. Only the lower triangle of each block is
#'   read.
#'
#' @return An \eqn{n \times K \times K} numeric array, lower triangular in
#'   its last two indices, with `Om[i, , ] == L[i, , ] %*% t(L[i, , ])` for
#'   every `i`. `NULL` as soon as any block fails to be positive definite,
#'   so a single bad observation declines for the whole sample.
#'
#' @seealso [info_blocks()] for the input,
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
#' In that case [chol_blocks()] returns `NULL`.
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
#'   `"opg"` (the default), `"bartlett"`, `"integrate"` or `"mc"`. Ignored
#'   when `expected` is `FALSE` or when the family computes its expected
#'   information exactly.
#' @param H The family's link-scale second derivatives, as
#'   [statmod_family_hessian()] returns them, when the caller has already
#'   computed them at `theta`; `NULL` computes them here.
#'
#' @return An \eqn{n \times K \times K} numeric array, symmetric in its last
#'   two indices, with \eqn{K} the number of distribution parameters.
#'
#' @seealso [chol_blocks()], which factorizes these,
#'   [statmod_information_at()] for the assembled \eqn{Z'\Omega Z}.
#'
#' @keywords internal
info_blocks <- function(spec, theta, expected = TRUE, approx = "opg",
                        H = NULL) {
  params <- spec@distrib@params
  K <- length(params)
  n <- spec@n_obs
  # given by the caller where it has already asked for them at this point
  if (is.null(H)) H <- statmod_family_hessian(spec, theta, expected, approx)
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
#' The result is sparse when any equation's design block is, and the solve
#' then goes to [sparse_augmented_solve()].
#'
#' @param design The design, as [statmod_design()] returns it, holding one
#'   block of columns per distribution parameter.
#' @param L The per-observation Cholesky factors, as [chol_blocks()] returns
#'   them, or `NULL`.
#'
#' @return An \eqn{nK \times p} matrix, dense or a \pkg{Matrix} sparse
#'   matrix according to the design's own storage, where \eqn{p} is the total
#'   number of coefficients across the equations. `NULL` when `L` is `NULL`
#'   or has a missing or infinite entry, so the `NULL` propagates from the
#'   factorization to the solve without a test at each step.
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
#' Hessian, omitting the rows that correspond to its null space.
#'
#' @details
#' # Why an eigendecomposition
#'
#' A penalty is positive semidefinite and is often rank deficient. A spline's
#' is deficient by exactly the dimension of its null space, the set of
#' coefficient directions that the penalty does not shrink. A Cholesky
#' factorization fails there, so the factor comes from an eigendecomposition
#' with the non-positive eigenvalues dropped, and the rows that they would
#' have contributed are absent. A zero matrix gives a factor with no rows.
#'
#' A non-convex penalty has an indefinite Hessian and no such factor.
#' `NULL` is returned and the caller falls back to the assembled route.
#'
#' # The diagonal shortcut
#'
#' A diagonal penalty is factored by taking the square root of its diagonal.
#' That is the same answer the eigendecomposition gives, since the
#' eigenvalues of a diagonal matrix are its diagonal entries, so the two
#' routes agree by construction and a test pins them together.
#'
#' The case is common. A ridge penalty and the penalty of a random effect
#' with independent coefficients are diagonal, and so is the penalty of a
#' smoother in the Demmler-Reinsch basis (see [modelterms7::s()]), which is
#' zero on the null-space coordinates and a multiple of the smoothing
#' parameter elsewhere. Any block-diagonal assembly of such penalties is
#' diagonal as well. The factor is recomputed at every iteration of the
#' scoring loop, so the diagonal route avoids one dense eigendecomposition
#' per iteration.
#'
#' @param S The penalty Hessian, a `p x p` symmetric matrix, dense or sparse.
#'   `Matrix::isDiagonal()` decides which route is taken, so a matrix stored
#'   as diagonal and one that only has zero off-diagonal entries are both
#'   treated as diagonal.
#'
#' @return A matrix with `ncol(S)` columns and one row per retained
#'   eigendirection, so at most `p` rows and fewer where the penalty has a
#'   null space. The result is a \pkg{Matrix} sparse matrix only where `S` is
#'   a diagonal matrix stored sparse; for any other `S`, including the zero
#'   matrix, it is a dense base matrix. `NULL` when the penalty is indefinite
#'   beyond the tolerance.
#'
#' @seealso [diagonal_sqrt()] for the diagonal route,
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
#' One row per coordinate that the penalty reaches, carrying the square root
#' of that coordinate's entry, which is the value of [penalty_sqrt()] where
#' the matrix is diagonal.
#'
#' @details
#' The thresholds are the eigen route's, applied to the diagonal, which for
#' a diagonal matrix is its spectrum. A negative entry beyond the tolerance
#' makes the penalty indefinite and there is no factor to return; an entry
#' at or below the tolerance is a null direction and contributes no row.
#'
#' The class of the result follows the argument's: a sparse `S` gives a
#' sparse factor and a dense `S` a dense one. [augmented_solve()] takes the
#' sparse route whenever either of its two factors is sparse, so the class is
#' not changed here.
#'
#' @param S A diagonal penalty Hessian, `p x p`, dense or sparse. Only its
#'   diagonal is read.
#' @param p Its dimension. Passed in so that the caller's own count is used
#'   and no dimension is inferred here.
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
#' The augmented form never forms \eqn{R'R} or the penalized information, so
#' a design with condition number \eqn{\kappa} is decomposed at \eqn{\kappa}
#' and not at \eqn{\kappa^2}, which the normal equations would give.
#'
#' @details
#' # The four routes
#'
#' The routes are tried in this order.
#'
#' 1. **Sparse QR**, when either factor is a \pkg{Matrix} object. The pair is
#'    passed to [sparse_augmented_solve()], which declines on a
#'    rank-deficient matrix.
#' 2. **SVD**, when `how = "svd"`. The singular values below
#'    \eqn{\max(\dim A)\,\epsilon\,\sigma_{\max}} are dropped and the
#'    increment is built from the retained right singular vectors, so a
#'    deficient system gets the minimum-norm solution.
#' 3. **The threaded triangular factor**, when `how = "qr"`, `threads > 1`
#'    and \eqn{m p^2} is at least `5e7`, where \eqn{m} is the number of rows
#'    of \eqn{A}. Only the triangular factor is read, so a kernel that
#'    produces it without accumulating \eqn{Q} does the whole job.
#' 4. **`qr()`**, the pivoted LINPACK route, which returns a rank and can
#'    drop columns. A sparse solve that declined continues here on
#'    densified factors.
#'
#' # Why the rank test is equilibrated
#'
#' The threaded route is taken only where the matrix has full rank by a
#' comfortable margin. The test reads the diagonal of the triangular factor
#' **divided by each column's norm**. Writing \eqn{T} for the triangular
#' factor of \eqn{A}, so that \eqn{A'A = T'T}, scaling the columns of \eqn{A}
#' scales the diagonal of \eqn{T} by the same factors, so the ratio is the
#' diagonal of a decomposition with unit column norms.
#'
#' The raw diagonal would report a matrix as near-singular whenever its
#' columns differ in size, and a large smoothing parameter produces that
#' difference between the penalty rows of its own block and an unpenalized
#' block. Scaling each column to unit norm removes differences in column size
#' from any source and leaves an exact collinearity exactly singular.
#'
#' The route is taken when the smallest equilibrated entry exceeds `1e-7`
#' times the largest, the tolerance that `dqrdc2` uses, so a matrix that the
#' pivoted route would call deficient is solved there and its rank is
#' returned.
#'
#' @param R The square-root design, as [sqrt_design()] returns it,
#'   \eqn{nK \times p}.
#' @param C The penalty's factor, as [penalty_sqrt()] returns it, with `p`
#'   columns and at most `p` rows. May have zero rows for an unpenalized
#'   model.
#' @param u The right-hand side, a numeric vector of length `p`.
#' @param how `"qr"` or `"svd"`. `"svd"` has no sparse counterpart, so with
#'   `how = "svd"` a sparse pair is densified first.
#' @param threads How many threads the triangular factor may use, a plain
#'   integer. With `1L` the threaded route is not used.
#'
#' @return A list with the elements `delta`, `rank` and `dropped`:
#'   \describe{
#'     \item{`delta`}{the increment, an unnamed numeric vector of length `p`.
#'       Zero in every coordinate that the pivoted route dropped.}
#'     \item{`rank`}{the rank used, an integer. Always `p` on the threaded
#'       and sparse routes, which decline on a deficient matrix instead of
#'       reporting a smaller rank, and the number of retained singular
#'       values on the `"svd"` route.}
#'     \item{`dropped`}{the integer positions of the columns that the pivoted
#'       `qr()` route left out (the aliased coordinates); `integer(0)` on the
#'       sparse, threaded and `"svd"` routes.}
#'   }
#'
#' @seealso [sqrt_design()] and [penalty_sqrt()] for the two factors,
#'   [sparse_augmented_solve()] for the sparse route.
#'
#' @keywords internal
augmented_solve <- function(R, C, u, how, threads = 1L) {
  if (isS4(R) || isS4(C)) {
    out <- sparse_augmented_solve(R, C, u, how)
    # it declines rather than dropping anything, so a solve it returns has
    # nothing aliased
    if (!is.null(out)) return(c(out, list(dropped = integer(0))))
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
    # the svd removes singular DIRECTIONS, which are combinations of columns
    # and name no coordinate, so nothing is reported aliased on this route
    return(list(delta = as.numeric(s$v %*% (dinv * crossprod(s$v, u))),
                rank = sum(keep), dropped = integer(0)))
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
          fast <- list(delta = as.numeric(d), rank = ncol(A),
                       dropped = integer(0))
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
  # The coordinates the pivot left out are the ALIASED ones, and this is the
  # only place in the package where a column is dropped: the sparse route
  # declines on a rank-deficient matrix and falls through to here, the
  # threaded factor is engaged only above a full-rank test, and the svd route
  # removes singular DIRECTIONS rather than columns and so reports none.
  # Naming them is what lets a fit report NA where base R reports NA instead
  # of an estimate of zero, which reads as a coefficient that was estimated.
  list(delta = delta, rank = qrA$rank,
       dropped = setdiff(seq_len(ncol(A)), piv))
}


#' Solve a Scoring Step From a Sparse Square-Root Design
#'
#' @description
#' The same increment [augmented_solve()] returns, taken through a
#' sparse QR of \eqn{[R;\ C]}.
#'
#' @details
#' # Why a QR
#'
#' The augmented form exists so that \eqn{X'X} is never formed and the
#' conditioning is never squared. A sparse QR is a QR and keeps that
#' property. A sparse Cholesky factorization of the normal equations would
#' be faster and would give it up.
#'
#' # The factor is not back-permuted
#'
#' A sparse QR reorders the columns to reduce fill, so \eqn{AP = QR} for the
#' permutation \eqn{P} that the decomposition chose, and
#' \eqn{(A'A)^{-1} = P(R'R)^{-1}P'}. The increment is two triangular solves
#' between a permutation and its inverse, which is the bookkeeping that the
#' dense route does with the pivot of `qr()`.
#'
#' The function calls `qrR(backPermute = FALSE)`. The back-permuted factor
#' is not triangular, so its diagonal does not give the rank, and a solve
#' against it is a general solve instead of two triangular ones.
#'
#' # The rank test
#'
#' The rank is read off the diagonal of the triangular factor, since a sparse
#' QR does not report a rank. The test is on the **Jacobi-equilibrated**
#' diagonal, each entry divided by its column's norm, and the route declines
#' when the smallest ratio is at most `ncol(A) * .Machine$double.eps` times
#' the largest.
#'
#' Without the equilibration the raw diagonal would reject matrices that the
#' dense route finds at full rank, because their columns differ in size (a
#' large smoothing parameter makes them so), and the solve would fall through
#' to the dense QR.
#'
#' The column norms come from `colSums(A^2)`, which acts on the `x` slot of
#' the sparse matrix; `colSums(A * A)` would form a binary operation between
#' two sparse matrices, which intersects their index sets and is slower.
#'
#' Where the matrix is rank deficient there is no unique increment, and this
#' route declines instead of choosing one. The caller falls back to the dense
#' route, which drops columns and reports the rank and the dropped columns.
#'
#' @param R The square-root design, \eqn{nK \times p}, sparse or dense.
#' @param C The penalty's factor, with `p` columns.
#' @param u The right-hand side, a numeric vector of length `p`.
#' @param how The decomposition requested. Only `"qr"` is served; `"svd"` has
#'   no sparse counterpart and declines.
#'
#' @return A list with `delta` (numeric, length `p`) and `rank` (integer,
#'   equal to `p`), or `NULL` when the route declines: `how` is not `"qr"`,
#'   the stacked matrix has fewer rows than columns, a column has zero norm,
#'   the factorization or the triangular solves fail, or the equilibrated
#'   diagonal shows that the matrix is deficient. Unlike [augmented_solve()],
#'   the list has no `dropped` element.
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
