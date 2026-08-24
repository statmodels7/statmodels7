# The thread policy inside the fitting layer: the count travels on the
# specification (spec@threads) and reaches the compiled kernels as an
# argument. Nothing here reads global state.

#' The Weighted Cross Product of the Assembly
#'
#' @description
#' Computes \eqn{A^\top \mathrm{diag}(w) B}, the weighted cross product the
#' information matrix and the outer Hessian are assembled from. Above a
#' measured work threshold, and only when the caller asked for more than one
#' thread, the work goes to the package's own threaded kernel. Everywhere
#' else it evaluates `crossprod(A * w, B)`, which is the expression the
#' kernel replaces.
#'
#' @details
#' # Why the answer does not depend on the thread count
#'
#' The kernel splits the work over the **elements of the output**: one thread
#' owns one entry \eqn{(j,k)} and accumulates its whole dot product, in
#' ascending row order, with the weight multiplied onto `A`'s entry before
#' the product is formed. No accumulation is ever split between threads, so
#' the sum is the same sequence of additions at any count and the result is
#' identical bit for bit at 1, 2 or 24 threads.
#'
#' That guarantee is about the kernel and not about the two routes agreeing
#' with each other. Engaging the kernel replaces a BLAS call, and an
#' optimized BLAS (OpenBLAS, Accelerate) blocks its accumulations into a
#' different order. The two routes therefore agree to the last bit on R's
#' reference BLAS and to the rounding of a single dot product elsewhere.
#'
#' # When the kernel is used
#'
#' All four conditions must hold: `threads > 1`, `A` and `B` are base dense
#' matrices, `w` has one entry per row, and the work
#' \eqn{n \times p_A \times p_B} is at least `2e5` multiply-adds. A sparse
#' design keeps its \pkg{Matrix} route, where a dense kernel would gain
#' nothing: the cost there is set by the number of stored nonzeros, and the
#' threshold above reads the dense shape.
#'
#' The threshold is a constant in the source and not an argument. It was
#' measured: the crossover where opening a parallel region costs what it
#' saves sits near \eqn{9 \times 10^4} multiply-adds on the development
#' machine (0.94x at \eqn{8 \times 10^4}, 1.7x at \eqn{10^5}, 10x at
#' \eqn{2 \times 10^6}, 19x to 20x at the shapes a real fit assembles), and
#' the gate is set above it, where a misjudgement in either direction costs
#' almost nothing.
#'
#' # Dimnames
#'
#' The kernel returns a bare matrix, so the column names of `A` and `B` are
#' put back as the row and column names of the result when either has them.
#' `crossprod()` does this itself, so both routes name their output alike.
#'
#' @param A,B Dense design blocks with the same number of rows, `n x pa` and
#'   `n x pb`. Either may be sparse, in which case the \pkg{Matrix} route is
#'   taken and the result is whatever `crossprod()` gives for that pair.
#' @param w Per-observation weights, a numeric vector of length `n`. A weight
#'   of any sign is accepted; nothing here assumes the working weights are
#'   positive.
#' @param threads The thread count, a plain integer as
#'   [numericals7::thread_count()] returns it. `1L`, the default, takes the
#'   sequential route unconditionally.
#'
#' @return A `pa x pb` numeric matrix, with the column names of `A` and `B`
#'   as its dimnames when either block carries them.
#'
#' @seealso [xtv()] and [xtx()], the other two threaded assembly kernels,
#'   which share this threshold.
#'
#' @keywords internal
wcrossprod <- function(A, w, B, threads = 1L) {
  if (threads > 1L && is.matrix(A) && is.matrix(B) &&
      length(w) == nrow(A) &&
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
#' The two per-sweep reads of the compiled coordinate descent, each one
#' vector of length `p` over an `n x p` design.
#'
#' `xtv()` computes \eqn{X^\top v}, one dot product per column. It is the
#' gradient read: `v` is the running residual, or the residual times the
#' working weights.
#'
#' `wxsq()` computes \eqn{\sum_i w_i x_{ij}^2} for each column \eqn{j}, the
#' column curvatures of a working model. The `n x p` matrix of squares is
#' never formed.
#'
#' Both take the package's threaded kernel above the same measured work
#' threshold [wcrossprod()] uses, and otherwise evaluate the expressions they
#' replace, `as.numeric(crossprod(X, v))` and
#' `as.numeric(crossprod(w, X^2))`.
#'
#' @details
#' # Why the answer does not depend on the thread count
#'
#' Each element of the output is one column's accumulation, owned in full by
#' one thread and summed in ascending row order. Nothing is split between
#' threads, so the result is identical bit for bit at any count.
#'
#' The elementwise product inside each term is formed in the order the
#' replaced expression forms it: `v` arrives already multiplied by whatever
#' the caller wanted, and in `wxsq()` the entry is squared before the weight
#' multiplies it. Against an optimized BLAS, which blocks its own
#' accumulations, the two routes agree to the rounding of one dot product.
#'
#' # When the kernel is used
#'
#' `threads > 1`, `X` a base dense matrix, the vector of length `nrow(X)`,
#' and \eqn{n \times p} at least `2e5`. A sparse `X` takes the
#' \pkg{Matrix} route.
#'
#' @param X A design block, `n x p`.
#' @param v The vector \eqn{X^\top v} is taken against, length `n`. Read by
#'   `xtv()` only.
#' @param w Per-observation weights, length `n`. Read by `wxsq()` only, and
#'   multiplied onto the squared entry.
#' @param threads The thread count, a plain integer as
#'   [numericals7::thread_count()] returns it. `1L` takes the sequential
#'   route unconditionally.
#'
#' @return An unnamed numeric vector of length `ncol(X)`: for `xtv()` the
#'   column dot products with `v`, for `wxsq()` the weighted sums of squares.
#'   Neither carries the column names of `X`, since both feed arithmetic
#'   indexed by position.
#'
#' @seealso [wcrossprod()] for the matrix-valued kernel and the measured
#'   threshold both share, [coord_fit()] for the descent that calls these.
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
#' threaded kernel above the measured work threshold and through
#' `crossprod(A)` everywhere else. This is what `fit_smooth()`'s subset route
#' assembles at every scoring iteration: the design has already been scaled
#' by the square root of the working weights, so no weight vector enters.
#'
#' @details
#' # Symmetry, and independence of the thread count
#'
#' Only the upper triangle is accumulated: one thread owns column \eqn{k} and
#' walks \eqn{j = 1, \ldots, k}, and each entry it computes is written to
#' both \eqn{(j,k)} and \eqn{(k,j)}. The two halves are therefore the same
#' number and not two roundings of one, so the result is exactly symmetric
#' with nothing to symmetrize away. Each pair belongs to exactly one thread,
#' so the mirrored write is disjoint, and no accumulation is split, so the
#' answer is identical bit for bit at any thread count.
#'
#' That halving is what the reference `dsyrk` behind `crossprod()` does too,
#' and it accumulates in the same ascending row order, so the two routes
#' agree to the last bit there. An optimized BLAS blocks its accumulations
#' and agrees to the rounding of one dot product.
#'
#' A full \eqn{p^2} version was measured and is not what ships: it was no
#' better than `dsyrk` at eight threads, so the second half of the work
#' bought nothing.
#'
#' # When the kernel is used
#'
#' `threads > 1`, `A` a base dense matrix, and \eqn{n \times p^2} at least
#' `2e5`. Note the \eqn{p^2}: this kernel writes \eqn{p^2} outputs where
#' [xtv()] writes \eqn{p}, so a design wide enough to reach the threshold
#' here can be too narrow to reach it there.
#'
#' @param A A dense design block, `n x p`, already scaled by the square root
#'   of the working weights. A sparse `A` takes the \pkg{Matrix} route and
#'   returns whatever `crossprod()` gives for it.
#' @param threads The thread count, a plain integer as
#'   [numericals7::thread_count()] returns it. `1L` takes the sequential
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
