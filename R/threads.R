# The thread policy inside the fitting layer: the count travels on the
# specification (spec@threads) and reaches the compiled kernels as an
# argument. Nothing here reads global state.

#' The Weighted Cross Product of the Assembly
#'
#' @description
#' `t(A) %*% (w * B)`, through the package's threaded kernel where the
#' thread count and the measured work threshold allow it, and through
#' `crossprod(A * w, B)` -- exactly the expression it replaces --
#' everywhere else.
#'
#' @details
#' The kernel decomposes over the ELEMENTS of the output, each dot product
#' accumulated in full by one thread in ascending row order with the weight
#' multiplied onto `A`'s entry first: the kernel's result does not
#' depend on the thread count, bit for bit. Engaging the kernel replaces the
#' BLAS expression, and an optimized BLAS (OpenBLAS, Accelerate) blocks its
#' accumulations, so the two implementations coincide to the last bit only
#' on the reference BLAS and to the rounding of one dot product elsewhere.
#' Only base dense matrices with a full-length weight are eligible; a sparse
#' design keeps its \pkg{Matrix} route, where a threaded dense kernel would
#' do nothing (piano_parallel.txt: the threshold is on the work, p and the
#' density, not on n alone).
#'
#' The threshold is internal and measured, not an argument: the crossover
#' where the cost of opening the region meets its gain sits near `9e4`
#' multiply-adds on this machine's grid (0.94x at 8e4, 1.7x at 1e5, 10x at
#' 2e6, 19-20x at the profile shapes), and the gate is set a factor of two
#' above it, where the asymmetry argument says an error either way costs
#' almost nothing.
#'
#' @param A,B Design blocks, `n x pa` and `n x pb`.
#' @param w The per-observation weights, length `n`.
#' @param threads The thread count, a plain integer.
#'
#' @return A `pa x pb` matrix.
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
#' `xtv()` is `as.numeric(crossprod(X, v))`, one dot product per
#' column; `wxsq()` is `as.numeric(crossprod(w, X^2))`, the column
#' curvatures of a working model, computed without materializing the
#' `n x p` square. Both run through the package's threaded kernels
#' above the same measured work gate as [wcrossprod()] and through
#' the exact expressions they replace everywhere else.
#'
#' @details
#' Each output element is accumulated in full by one thread in ascending row
#' order, with any elementwise product rounded exactly as the replaced
#' expression rounds it (`v` arrives precomputed; the square is taken
#' before the weight multiplies it), so the kernel's result does not depend
#' on the thread count, bit for bit; against an optimized BLAS's own
#' accumulation order it agrees to the rounding of one dot product. These
#' are the per-sweep reads the plan's
#' section 0quinquies classifies as reductions: they are parallel over the
#' OUTPUT elements, never over a split of one accumulation.
#'
#' @param X A design block, `n x p`.
#' @param v,w Per-observation vectors of length `n`.
#' @param threads The thread count, a plain integer.
#'
#' @return A numeric vector of length `ncol(X)`.
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
#' `crossprod(A)` through the threaded kernel where the count and the
#' measured work gate allow it, and through `crossprod(A)` itself
#' everywhere else. What `fit_smooth()`'s subset route assembles at
#' every scoring iteration.
#'
#' @details
#' Element \eqn{(j, k)} is one dot product accumulated in full by one
#' thread in ascending row order, and elements \eqn{(j, k)} and
#' \eqn{(k, j)} are the same products summed in the same order, so the
#' result is exactly symmetric and the kernel does not depend on the thread
#' count, bit for bit. The reference `dsyrk` behind `crossprod`
#' accumulates in the same order; an optimized BLAS does not, and agrees to
#' the rounding of one dot product.
#'
#' @param A A dense design block, `n x p`.
#' @param threads The thread count, a plain integer.
#'
#' @return A `p x p` matrix.
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
