# The thread policy inside the fitting layer: the count travels on the
# specification (spec@threads) and reaches the compiled kernels as an
# argument. Nothing here reads global state.

#' The Weighted Cross Product of the Assembly
#'
#' @description
#' Computes \eqn{A^\top \mathrm{diag}(w) B}, the weighted cross product from
#' which the information matrix and the outer Hessian are assembled. Above a
#' work threshold, and only when more than one thread is requested, the work
#' goes to the package's own threaded kernel. Everywhere else the function
#' evaluates `crossprod(A * w, B)`, which is the expression that the kernel
#' replaces.
#'
#' @details
#' # Why the answer does not depend on the thread count
#'
#' The kernel splits the work over the elements of the output: one thread
#' owns one entry \eqn{(j,k)} and accumulates its whole dot product in
#' ascending row order, multiplying the weight onto `A`'s entry before
#' the product is formed. The kernel never splits an accumulation between
#' threads, so the sum is the same sequence of additions at any count and the
#' result is identical bit for bit at any thread count.
#'
#' That guarantee concerns only the kernel and does not cover the agreement
#' of the two routes with each other: engaging the kernel replaces a BLAS
#' call, and an optimized BLAS (OpenBLAS, Accelerate) blocks its
#' accumulations into a different order. The two routes therefore agree to
#' the last bit on R's reference BLAS and to the rounding of a single dot
#' product elsewhere.
#'
#' # When the kernel is used
#'
#' All of the following must hold: `threads > 1`; `A` and `B` are base dense
#' matrices with the same number of rows; `w` has one entry per row; and the
#' work \eqn{n \times p_A \times p_B} is at least `2e5` multiply-adds.
#' Otherwise `crossprod(A * w, B)` is evaluated, which signals its own error
#' for blocks of different heights. A sparse design keeps its \pkg{Matrix}
#' route, where a dense kernel would gain nothing: the cost there is set by
#' the number of stored nonzeros, and the threshold reads the dense shape.
#'
#' The threshold is a constant in the source, and the function has no argument
#' for it. It lies above the work at which opening a parallel region costs what
#' it saves.
#'
#' # Dimnames
#'
#' The kernel returns a bare matrix, so the column names of `A` and `B` are
#' put back as the row and column names of the result when either has them.
#' `crossprod()` does this itself, so both routes name their output alike.
#'
#' @param A,B Design blocks with the same number of rows, `n x pa` and
#'   `n x pb`. If either is not a base dense matrix, the sequential route is
#'   taken and the result is whatever `crossprod()` gives for that pair.
#' @param w Per-observation weights, a numeric vector of length `n`. A weight
#'   of any sign is accepted; the function does not assume that the working
#'   weights are positive.
#' @param threads The thread count, a plain integer as
#'   [numericals7::local_threads()] validates it. `1L`, the default, takes the
#'   sequential route unconditionally.
#'
#' @return A `pa x pb` numeric matrix, with the column names of `A` and `B`
#'   as its dimnames when either block carries them.
#'
#' @seealso [xtv()], [wxsq()] and [xtx()], the other threaded assembly
#'   kernels. They use the same constant, \eqn{2 \times 10^5}, on their own
#'   measures of work: \eqn{n p} for [xtv()] and [wxsq()], and
#'   \eqn{n p^2} for [xtx()].
#'
#' @keywords internal
wcrossprod <- function(A, w, B, threads = 1L) {
  if (threads > 1L && is.matrix(A) && is.matrix(B) &&
      length(w) == nrow(A) && nrow(B) == nrow(A) &&
      as.double(nrow(A)) * ncol(A) * ncol(B) >= 2e5) {
    out <- wcrossprod_cpp(A, w, B, threads)
    cn <- list(colnames(A), colnames(B))
    if (!is.null(cn[[1L]]) || !is.null(cn[[2L]])) dimnames(out) <- cn
    return(out)
  }
  crossprod(A * w, B)
}

#' The Two Vector Products of the Coordinate-Descent Loop
#'
#' @description
#' `xtv()` and `wxsq()` compute the two per-sweep reads of the compiled
#' coordinate descent; each returns one vector of length `p` over an `n x p`
#' design.
#'
#' `xtv()` computes \eqn{X^\top v}, one dot product per column. It is the
#' gradient read: `v` is the running residual, or the residual times the
#' working weights.
#'
#' `wxsq()` computes \eqn{\sum_i w_i x_{ij}^2} for each column \eqn{j}, the
#' column curvatures of a working model. The threaded kernel does not form the
#' `n x p` matrix of squares, whereas the sequential expression does.
#'
#' Both take the package's threaded kernel above a work threshold, the one
#' constant that [wcrossprod()] also uses, and otherwise evaluate the
#' expressions that they replace, `as.numeric(crossprod(X, v))` and
#' `as.numeric(crossprod(w, X^2))`.
#'
#' @details
#' # Why the answer does not depend on the thread count
#'
#' Each element of the output is one column's accumulation, owned in full by
#' one thread and summed in ascending row order. The kernel never splits an
#' accumulation between threads, so the result is identical bit for bit at
#' any count.
#'
#' The elementwise product inside each term is formed in the order in which
#' the replaced expression forms it: `v` arrives already multiplied by
#' whatever the caller wanted, and in `wxsq()` the entry is squared before
#' the weight multiplies it. Against an optimized BLAS, which blocks its own
#' accumulations, the two routes agree to the rounding of one dot product.
#'
#' # When the kernel is used
#'
#' The kernel is used when `threads > 1`, `X` is a base dense matrix, the
#' vector has length `nrow(X)`, and \eqn{n \times p} is at least `2e5`. A
#' sparse `X` takes the \pkg{Matrix} route.
#'
#' @param X A design block, `n x p`.
#' @param v The vector against which \eqn{X^\top v} is taken, length `n`.
#'   Only `xtv()` reads it.
#' @param w Per-observation weights, length `n`. Only `wxsq()` reads it, and
#'   it multiplies the squared entry.
#' @param threads The thread count, a plain integer as
#'   [numericals7::local_threads()] validates it. `1L` takes the sequential
#'   route unconditionally.
#'
#' @return An unnamed numeric vector of length `ncol(X)`: for `xtv()` the
#'   column dot products with `v`, for `wxsq()` the weighted sums of squares.
#'   Both functions omit the column names of `X`, because their results feed
#'   arithmetic indexed by position.
#'
#' @seealso [wcrossprod()] for the matrix-valued kernel and the work
#'   threshold, [coord_fit()] for the descent that calls these.
#'
#' @keywords internal
xtv <- function(X, v, threads = 1L) {
  if (threads > 1L && is.matrix(X) && length(v) == nrow(X) &&
      as.double(nrow(X)) * ncol(X) >= 2e5) {
    return(xtv_cpp(X, v, threads))
  }
  as.numeric(crossprod(X, v))
}

#' @rdname xtv
wxsq <- function(X, w, threads = 1L) {
  if (threads > 1L && is.matrix(X) && length(w) == nrow(X) &&
      as.double(nrow(X)) * ncol(X) >= 2e5) {
    return(wxsq_cpp(X, w, threads))
  }
  as.numeric(crossprod(w, X^2))
}

#' The Unweighted Cross Product of a Square-Root Design
#'
#' @description
#' Computes \eqn{A^\top A} for a square-root design, through the package's
#' threaded kernel above a work threshold and through
#' `crossprod(A)` everywhere else. The argument is a design that has already
#' been scaled by the square root of the working weights, so the function
#' takes no weight vector.
#'
#' @details
#' # Symmetry, and independence of the thread count
#'
#' Only the upper triangle is accumulated: one thread owns column \eqn{k} and
#' walks \eqn{j = 1, \ldots, k}, and each entry that it computes is written
#' to both \eqn{(j,k)} and \eqn{(k,j)}. The two halves are therefore the same
#' number, rounded once, so the result is exactly symmetric and needs no
#' symmetrization. Each pair belongs to exactly one thread, so the mirrored
#' write is disjoint. The kernel splits no accumulation, so the answer is
#' identical bit for bit at any thread count.
#'
#' The reference `dsyrk` behind `crossprod()` also halves the work in this
#' way and accumulates in the same ascending row order, so the two routes
#' agree to the last bit there. An optimized BLAS blocks its accumulations
#' and agrees to the rounding of one dot product.
#'
#' Only the upper triangle is computed, as `dsyrk` does, because computing all
#' \eqn{p^2} entries gave no gain.
#'
#' # When the kernel is used
#'
#' The kernel is used when `threads > 1`, `A` is a base dense matrix, and
#' \eqn{n \times p^2} is at least `2e5`. This kernel writes \eqn{p^2} outputs
#' where [xtv()] writes \eqn{p}, so a design wide enough to reach the
#' threshold here can be too narrow to reach it there.
#'
#' @param A A dense design block, `n x p`, already scaled by the square root
#'   of the working weights. A sparse `A` takes the \pkg{Matrix} route and
#'   returns whatever `crossprod()` gives for it.
#' @param threads The thread count, a plain integer as
#'   [numericals7::local_threads()] validates it. `1L` takes the sequential
#'   route unconditionally.
#'
#' @return A `p x p` symmetric numeric matrix, with the column names of `A`
#'   as both dimnames when `A` has them.
#'
#' @seealso [wcrossprod()] for the weighted two-block form,
#'   [xtv()] for the vector-valued kernels.
#'
#' @keywords internal
xtx <- function(A, threads = 1L) {
  if (threads > 1L && is.matrix(A) &&
      as.double(nrow(A)) * ncol(A) * ncol(A) >= 2e5) {
    out <- xtx_cpp(A, threads)
    cn <- colnames(A)
    if (!is.null(cn)) dimnames(out) <- list(cn, cn)
    return(out)
  }
  crossprod(A)
}
