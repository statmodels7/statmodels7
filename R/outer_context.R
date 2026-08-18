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
#' \eqn{K = H + S} with the criterion's own information, and its inverse, which
#' the gradient and the Hessian both read.
#'
#' @param expected Whether \eqn{H} is the expected information, which is what
#'   the criterion carries under \code{reml(hessian = "expected")}. It used to
#'   be hard-coded to the observed one, correctly, because the exact gradient
#'   ran on no other route; admitting the expected route makes the criterion's
#'   determinant a different matrix, and reading the wrong one would be a
#'   gradient of the wrong function.
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
ctx_penalized <- function(ctx, spec, design, coef, hyper, expected = FALSE) {
  build <- function() {
    H <- ctx_information(ctx, spec, design, coef, hyper, expected,
                         ctx_approx(ctx))
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
  # the two informations give two different matrices, so they are two cache
  # entries: a criterion reading the expected one must not be handed the
  # observed one's factorization because a previous reader asked for it first
  slot <- if (expected) "penalized_expected" else "penalized_observed"
  if (is.null(ctx[[slot]])) {
    # the miss is recorded so a second reader does not retry a factorization
    # that has already failed
    assign(slot, list(value = build()), envir = ctx)
  }
  ctx[[slot]]$value
}


#' What the Marginal Criterion Can Resolve at This Fit
#'
#' @description
#' The size of the difference in the criterion that carries no information,
#' measured at the point rather than assumed.
#'
#' @details
#' The criterion is read at the penalized mode, and the inner fit stops short
#' of that mode by whatever its rule allows. Writing \eqn{g} for the score it
#' stopped at and \eqn{K = H + S} for the penalized information, the mode is
#' out by \eqn{\delta\beta = K^{-1} g}, and the criterion read at
#' \eqn{\hat\beta - \delta\beta} instead of at \eqn{\hat\beta} differs by an
#' amount that is exactly what an evaluation from a different warm start would
#' differ by. One assembly of the criterion at given coefficients answers it,
#' with no refit.
#'
#' \strong{The quadratic form alone is not enough}, and the reason is
#' structural rather than a matter of accuracy: the criterion carries
#' \eqn{-\tfrac{1}{2}\log\lvert K(\beta)\rvert}, which is NOT stationary in
#' \eqn{\beta}, so a mode error enters it at FIRST order. Measured,
#' \eqn{\tfrac{1}{2} g' K^{-1} g} is right to one per cent at an inner
#' tolerance of \code{1e-4}, where the second-order term dominates, and
#' undershoots by 50 to 1000 times at \code{1e-6} and below, where the first
#' order does.
#'
#' Measured against the spread of the criterion at one hyperparameter reached
#' from six different warm starts, over four shapes and five inner tolerances,
#' the displaced reading tracks it across six orders of magnitude -- from
#' \code{8.6e-2} down to \code{2.2e-7} -- at a ratio between 0.05 and 0.99, and
#' separates shapes that a formula cannot: at an inner tolerance of \code{1e-6}
#' it reads \code{1.7e-4} on a random intercept over 500 levels against
#' \code{7.9e-8} on a gaussian smooth.
#'
#' It reads LOW by two to three times, the spread being a range over six paths
#' where this is one deviation, and that is the side to err on. A resolution
#' smaller than the truth leaves the search where it was; one larger stops a
#' healthy search short.
#'
#' @param st An evaluation's stored state, carrying \code{par}, \code{split},
#'   \code{score}, \code{cf}, \code{hy}, \code{ctx} and \code{value}.
#' @param spec The specification.
#' @param design The design.
#' @param method An \code{\link{OuterMethod}}.
#' @param criterion_at A function of \code{(cf, hy, par, ctx)} returning what
#'   THIS search is running, which is a prediction-error criterion or a
#'   marginal one. It is passed in rather than chosen here: reading the
#'   marginal criterion of a fit whose search is \code{\link{aic}} answers for
#'   a quantity the search never sees, and the number that comes back stopped
#'   two such fits short of their own optimum.
#'
#' @return A single positive number, or \code{NA_real_} where the pieces are
#'   not available.
#'
#' @seealso \code{\link{outer_fit}}, \code{\link[optimizers7]{armijo}}
#'
#' @keywords internal
criterion_resolution <- function(st, spec, design, method, criterion_at) {
  if (!isTRUE(st$ok) || is.null(st$score) || is.null(st$split)) {
    return(NA_real_)
  }
  if (!all(is.finite(st$score))) return(NA_real_)
  expected <- identical(method@hessian, "expected")
  pen <- ctx_penalized(st$ctx, spec, design, st$cf, st$hy, expected)
  if (is.null(pen)) return(NA_real_)
  db <- tryCatch(as.numeric(as.matrix(pen$inv) %*% st$score),
                 error = function(e) NULL)
  if (is.null(db) || !all(is.finite(db))) return(NA_real_)
  par2 <- st$par - db
  # a FRESH context, the displaced coefficients not being the ones the cached
  # information, penalty and factorization were built at
  m2 <- tryCatch(criterion_at(st$split(par2), st$hy, par2, NULL),
                 error = function(e) NULL)
  if (is.null(m2) || !is.finite(m2$value)) return(NA_real_)
  out <- abs(m2$value - st$value)
  if (!is.finite(out) || out <= 0) NA_real_ else out
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
#' @param expected Which information \code{pen} was built with, which is the
#'   cache's key: the projection is of THAT matrix, so one entry cannot serve
#'   both. Nothing reaches it today, a search holding one
#'   \code{\link{OuterMethod}} throughout, but the slot is keyed rather than
#'   shared because the twin defect in \code{\link{ctx_penalized}} was
#'   unreachable in exactly the same way until it was not.
#'
#' @return A square matrix, or \code{NULL}.
#'
#' @keywords internal
ctx_trace_matrix <- function(ctx, pen, basis, expected = FALSE) {
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
  slot <- if (expected) "trace_expected" else "trace_observed"
  if (is.null(ctx[[slot]])) assign(slot, list(value = build()), envir = ctx)
  ctx[[slot]]$value
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


#' How the Penalized Matrix Moves With the Coefficients
#'
#' @description
#' The array \code{\link{u_vector}} contracts against the leverage diagonal,
#' together with the builder that reads it: the third derivative of the
#' log-density where the criterion uses the observed information, and the
#' derivative of the expected information where it uses the expected one.
#'
#' @details
#' \eqn{K} enters the criterion through its determinant, so the gradient needs
#' \eqn{\partial K/\partial\beta}. With the observed information that is
#' \eqn{-\ell'''}, which every family carries in closed form. With the expected
#' one it is \eqn{-\partial\,\mathbb{E}[\ell'']/\partial\eta}, which is NOT
#' \eqn{-\mathbb{E}[\ell''']}: differentiating an expectation moves the measure
#' as well as the integrand, and the missing piece
#' \eqn{\mathbb{E}[\ell_{ab}\ell_{c}]} is a mixed moment no Bartlett identity
#' isolates -- the third ties the symmetrized sum, not the single term.
#'
#' \pkg{distributions7} supplies it as
#' \code{\link[distributions7]{distrib_dexpected_hessian}}. The two arrays are
#' keyed differently, the observed one being symmetric in all three indices and
#' the expected one in its first two only, so the key builder travels with the
#' array rather than being assumed by the consumer.
#'
#' @param ctx A context, or \code{NULL}.
#' @param spec,design,coef,hyper The fallback arguments.
#' @param method An \code{\link{OuterMethod}}.
#'
#' @return A list with \code{deriv} and \code{key}.
#'
#' @seealso \code{\link{u_vector}}, \code{\link{outer_gradient_ok}}
#'
#' @keywords internal
ctx_kmove <- function(ctx, spec, design, coef, hyper, method) {
  params <- spec@distrib@params
  if (!identical(method@hessian, "expected")) {
    d3 <- ctx_deriv(ctx, spec, design, coef, hyper, 3L)
    keys <- names(d3)
    return(list(deriv = d3,
                key = function(a, b, k) d3_key(params, a, b, k, keys)))
  }
  build <- function() {
    th <- ctx_theta(ctx, spec, design, coef, hyper)
    distributions7::distrib_dexpected_hessian(spec@distrib, spec@response, th,
                                              scale = "link",
                                              approx = ctx_approx(ctx))
  }
  v <- if (ctx_usable(ctx, coef, hyper)) {
    if (is.null(ctx$dexp)) ctx$dexp <- build()
    ctx$dexp
  } else build()
  list(deriv = v,
       key = function(a, b, k) distributions7::dexpected_key(params, a, b, k))
}


#' The Approximation a Context Was Built With
#'
#' @param ctx A context, or \code{NULL}.
#'
#' @return A single string.
#'
#' @keywords internal
ctx_approx <- function(ctx) {
  if (is.null(ctx) || is.null(ctx$approx)) "bartlett" else ctx$approx
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
                                     scale = "link", threads = spec@threads)
    } else {
      distributions7::distrib_deriv4(spec@distrib, spec@response, th,
                                     scale = "link", threads = spec@threads)
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
