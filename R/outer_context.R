#' One Evaluation Point, Shared
#'
#' @description
#' The quantities the marginal criterion, its gradient and its Hessian all
#' read at one \eqn{(\beta, \theta)}: the linear predictors, the information,
#' the penalty's Hessian, their sum and its factorization.
#'
#' @details
#' The three consumers each used to assemble these for themselves, so at one
#' point the information was built three times and the penalized matrix
#' factorized twice, and \code{\link{statmod_marginal_hess}} additionally
#' calls the gradient, which repeated the whole of it a fourth time. Measured
#' by \code{Rprof}'s \code{by.total} on a random intercept over 500 levels at
#' 20000 observations, the gradient and the Hessian together accounted for 128
#' per cent of the fit -- the overlap being exactly that repetition.
#'
#' The context is an environment, so an accessor fills it in place and later
#' readers find the quantity already there. Each is computed on first demand
#' and never speculatively: the criterion alone does not need an inverse, and
#' a search running without an exact gradient must not pay for one.
#'
#' Passing \code{NULL} wherever a context is accepted restores the old
#' behaviour exactly, which is what keeps every existing caller -- and a
#' caller's own code -- working unchanged.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param approx The approximation for the expected information.
#'
#' @return An environment carrying the point and, as they are asked for, the
#'   quantities derived from it.
#'
#' @seealso \code{\link{statmod_marginal}}, \code{\link{statmod_marginal_grad}},
#'   \code{\link{statmod_marginal_hess}}
#'
#' @keywords internal
outer_context <- function(spec, design, coef, hyper, approx = "bartlett") {
  e <- new.env(parent = emptyenv())
  e$spec <- spec
  e$design <- design
  e$coef <- coef
  e$hyper <- hyper
  e$approx <- approx
  e
}


#' Refuse a Context That Belongs Somewhere Else
#'
#' @description
#' Stops when a context is handed to a consumer reading a different point.
#'
#' @details
#' A stale cached information is a silently wrong gradient, which is the one
#' failure a shared cache can introduce and the one that would be hardest to
#' find. The comparison is over the coefficients and the hyperparameters,
#' which is linear in their length and so negligible beside the \eqn{O(np^2)}
#' assembly it guards.
#'
#' @param ctx A context, or \code{NULL}.
#' @param coef The coefficients the caller is reading.
#' @param hyper The hyperparameters the caller is reading.
#'
#' @return \code{TRUE} when the context is usable, \code{FALSE} when it is
#'   \code{NULL}; an error when it belongs to another point.
#'
#' @keywords internal
ctx_usable <- function(ctx, coef, hyper) {
  if (is.null(ctx)) return(FALSE)
  if (!identical(ctx$coef, coef) || !identical(ctx$hyper, hyper)) {
    stop(paste0("the evaluation context was built at a different point than\n",
                "  the one being read. A context belongs to one (coef,",
                " hyper);\n  build a new one rather than reusing it."),
         call. = FALSE)
  }
  TRUE
}


#' The Information at the Context's Point
#'
#' @param ctx A context, or \code{NULL}.
#' @param spec,design,coef The fallback arguments, used when \code{ctx} is
#'   \code{NULL}.
#' @param hyper The hyperparameters, for the context's own check.
#' @param expected Whether the expected information is wanted.
#' @param approx The approximation for the expected information.
#'
#' @return A square matrix.
#'
#' @keywords internal
ctx_information <- function(ctx, spec, design, coef, hyper, expected,
                            approx = "bartlett") {
  if (!ctx_usable(ctx, coef, hyper)) {
    return(statmod_information_at(spec, coef, design, expected, approx))
  }
  key <- if (expected) "H_expected" else "H_observed"
  v <- ctx[[key]]
  if (is.null(v)) {
    v <- statmod_information_at(ctx$spec, ctx$coef, ctx$design, expected,
                                ctx$approx)
    assign(key, v, envir = ctx)
  }
  v
}


#' The Penalty's Hessian at the Context's Point
#'
#' @details
#' The non-finite entries are zeroed here rather than by each caller: three
#' of them did it separately and a fourth would have had to remember.
#'
#' @param ctx A context, or \code{NULL}.
#' @param spec,design,coef,hyper The fallback arguments.
#'
#' @return A square matrix.
#'
#' @keywords internal
ctx_penalty <- function(ctx, spec, design, coef, hyper) {
  build <- function(sp, cf, hy, dz) {
    S <- statmod_penalty_at(sp, cf, hy, dz, "hessian")
    S <- zap_nonfinite(S)
    S
  }
  if (!ctx_usable(ctx, coef, hyper)) return(build(spec, coef, hyper, design))
  if (is.null(ctx$S)) {
    ctx$S <- build(ctx$spec, ctx$coef, ctx$hyper, ctx$design)
  }
  ctx$S
}


#' The Penalized Matrix and Its Inverse
#'
#' @description
#' \eqn{K = H + S} with the OBSERVED information, and its inverse, which the
#' gradient and the Hessian both read.
#'
#' @details
#' \strong{The storage is the matrix's own.} \eqn{H} is already sparse
#' wherever the design is -- a grouping indicator, a factor \code{by}, a
#' \code{linpar} over many levels -- and it was being densified here only
#' because the penalty's accumulator is a base matrix. Where the sum is large
#' enough and sparse enough to be worth it (\code{\link{worth_sparse}}) it is
#' kept sparse and factorized as such: measured on a random intercept over 500
#' levels, p = 503 at a density of 0.014, the factorization and its
#' log-determinant cost 0.102 ms against 10.811 ms dense, and the full inverse
#' 3.280 ms against 25.000 ms, each route timed WITH its own factorization.
#' End to end that is 1.25x at 500 levels and 2.01x at 1000, the difference
#' between the operation and the fit being the lesson this package records
#' three times over: removing the dearer half leaves the cheaper one. Nothing
#' here asks which term produced the matrix; both quantities are read off the
#' matrix.
#'
#' ⚠️ Those figures are the ones measured with \pkg{Matrix}'s factorization
#' CACHE defeated. \code{Matrix::Cholesky} stores its result in the matrix's
#' \code{factors} slot, so a benchmark that refactorizes the same object
#' measures a cache hit -- 0.004 ms rather than 0.102 -- and reports the
#' sparse route as three times better than it is. A fit never gets that hit,
#' the penalized matrix being a new one at every point.
#'
#' \strong{The inverse stays dense whatever the factor is.} Its readers take
#' full matrix products against it -- \code{\link{block_leverage}} and the
#' Hessian's pair loop -- so the inverse of a sparse matrix being dense costs
#' nothing sparsity could have kept. What the sparse factor buys is the cost
#' of producing it.
#'
#' \code{NULL} is returned where the matrix is not positive definite, which is
#' the answer both callers already gave there.
#'
#' @param ctx A context, or \code{NULL}.
#' @param spec,design,coef,hyper The fallback arguments.
#'
#' @return A list with \code{K}, \code{inv} and \code{logdet}, or \code{NULL}.
#'
#' @seealso \code{\link{pd_factor}}
#'
#' @keywords internal
ctx_penalized <- function(ctx, spec, design, coef, hyper) {
  build <- function() {
    H <- ctx_information(ctx, spec, design, coef, hyper, FALSE)
    S <- ctx_penalty(ctx, spec, design, coef, hyper)
    K <- H + S
    fac <- pd_factor(K)
    if (!isTRUE(fac$ok)) return(NULL)
    inv <- if (isTRUE(fac$sparse)) {
      as.matrix(Matrix::solve(fac$factor, Matrix::Diagonal(ncol(K))))
    } else if (!is.null(fac$factor)) {
      # the factor the verdict was read off, so the same matrix is not
      # factorized twice to invert it
      chol2inv(fac$factor)
    } else {
      # the dense route reached its answer through the eigendecomposition,
      # where there is no factor to invert from
      ch <- tryCatch(chol(as_dense(K)), error = function(e) NULL)
      if (is.null(ch)) return(NULL)
      chol2inv(ch)
    }
    if (is.null(inv)) return(NULL)
    list(K = K, inv = inv, logdet = fac$logdet)
  }
  if (!ctx_usable(ctx, coef, hyper)) return(build())
  if (is.null(ctx$penalized)) {
    # the miss is recorded so a second reader does not retry a factorization
    # that has already failed
    ctx$penalized <- list(value = build())
  }
  ctx$penalized$value
}


#' The Matrix the Traces Are Taken Against
#'
#' @description
#' \eqn{M}, which is \eqn{K^{-1}} when nothing is projected away and the
#' projected inverse otherwise.
#'
#' @param ctx A context, or \code{NULL}.
#' @param pen The result of \code{\link{ctx_penalized}}.
#' @param basis The integrated subspace, or \code{NULL}.
#'
#' @return A square matrix, or \code{NULL}.
#'
#' @keywords internal
ctx_trace_matrix <- function(ctx, pen, basis) {
  if (is.null(basis)) return(pen$inv)
  build <- function() {
    # K may now be sparse, and the projection onto a dense basis is dense
    # whatever it was; the result is read as a full matrix by
    # block_leverage() and by the Hessian's pair loop, so it is a base matrix
    # here as it was before the storage became the matrix's own choice
    inner <- tryCatch(chol2inv(chol(as_dense(crossprod(basis,
                                                      pen$K %*% basis)))),
                      error = function(e) NULL)
    if (is.null(inner)) return(NULL)
    basis %*% inner %*% t(basis)
  }
  if (is.null(ctx)) return(build())
  if (is.null(ctx$trace_matrix)) ctx$trace_matrix <- list(value = build())
  ctx$trace_matrix$value
}


#' The Linear Predictors' Parameters at the Context's Point
#'
#' @param ctx A context, or \code{NULL}.
#' @param spec,design,coef,hyper The fallback arguments.
#'
#' @return A named list, one entry per distribution parameter.
#'
#' @keywords internal
ctx_theta <- function(ctx, spec, design, coef, hyper) {
  if (!ctx_usable(ctx, coef, hyper)) {
    return(statmod_eta(spec, design, coef)$theta)
  }
  if (is.null(ctx$theta)) {
    ctx$theta <- statmod_eta(ctx$spec, ctx$design, ctx$coef)$theta
  }
  ctx$theta
}


#' The Per-Observation Diagonals at the Context's Point
#'
#' @description
#' \code{\link{block_leverage}}, computed once and read by the gradient's
#' contraction and by every pair of the Hessian.
#'
#' @param ctx A context, or \code{NULL}.
#' @param design The design.
#' @param M The matrix the traces are taken against.
#' @param params,npar,offs The block bookkeeping.
#'
#' @return A list of lists of vectors.
#'
#' @keywords internal
ctx_leverage <- function(ctx, design, M, params, npar, offs) {
  if (is.null(ctx)) return(block_leverage(design, M, params, npar, offs))
  if (is.null(ctx$leverage)) {
    ctx$leverage <- block_leverage(design, M, params, npar, offs)
  }
  ctx$leverage
}


#' A Higher Derivative of the Log-Density at the Context's Point
#'
#' @description
#' The third or fourth derivative on the link scale, which the gradient's
#' contraction and the Hessian's two contractions all read.
#'
#' @param ctx A context, or \code{NULL}.
#' @param spec,design,coef,hyper The fallback arguments.
#' @param order \code{3} or \code{4}.
#'
#' @return A named list of vectors.
#'
#' @keywords internal
ctx_deriv <- function(ctx, spec, design, coef, hyper, order) {
  th <- ctx_theta(ctx, spec, design, coef, hyper)
  build <- function() {
    if (order == 3L) {
      distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                     scale = "link")
    } else {
      distributions7::distrib_deriv4(spec@distrib, spec@response, th,
                                     scale = "link")
    }
  }
  if (!ctx_usable(ctx, coef, hyper)) return(build())
  key <- paste0("d", order)
  v <- ctx[[key]]
  if (is.null(v)) {
    v <- build()
    assign(key, v, envir = ctx)
  }
  v
}
