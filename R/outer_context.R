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
#' factorized twice, and [statmod_marginal_hess()] additionally
#' calls the gradient, which repeated the whole of it a fourth time. Measured
#' by `Rprof`'s `by.total` on a random intercept over 500 levels at
#' 20000 observations, the gradient and the Hessian together accounted for 128
#' per cent of the fit, the overlap being exactly that repetition.
#'
#' The context is an environment, so an accessor fills it in place and later
#' readers find the quantity already there. Each is computed on first demand
#' and never speculatively: the criterion alone does not need an inverse, and
#' a search running without an exact gradient must not pay for one.
#'
#' Passing `NULL` wherever a context is accepted restores the earlier
#' behavior exactly, so every existing caller goes on working unchanged, a
#' caller's own code included.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param approx The approximation for the expected information.
#'
#' @return An environment carrying the point and, as they are asked for, the
#'   quantities derived from it.
#'
#' @seealso [statmod_marginal()], [statmod_marginal_grad()],
#'   [statmod_marginal_hess()]
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
#' @param ctx A context, or `NULL`.
#' @param coef The coefficients the caller is reading.
#' @param hyper The hyperparameters the caller is reading.
#'
#' @return `TRUE` when the context is usable, `FALSE` when it is
#'   `NULL`; an error when it belongs to another point.
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
#' @param ctx A context, or `NULL`.
#' @param spec,design,coef The fallback arguments, used when `ctx` is
#'   `NULL`.
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
#' The non-finite entries are zeroed here, once. Three callers did it
#' separately before, and a fourth would have had to remember.
#'
#' @param ctx A context, or `NULL`.
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
#' @param expected Whether \eqn{H} is the expected information, as the
#'   criterion carries under `reml(hessian = "expected")`. It used to
#'   be hard-coded to the observed one, correctly, because the exact gradient
#'   ran on no other route; admitting the expected route makes the criterion's
#'   determinant a different matrix, and reading the wrong one would be a
#'   gradient of the wrong function.
#'
#' @details
#' **The storage is the matrix's own.** \eqn{H} is already sparse wherever
#' the design is, which covers a grouping indicator, a factor `by` and a
#' `linpar` over many levels, and it was being densified here only because
#' the penalty's accumulator is a base matrix. Where the sum is large
#' enough and sparse enough to be worth it ([worth_sparse()]) it is
#' kept sparse and factorized as such: measured on a random intercept over 500
#' levels, p = 503 at a density of 0.014, the factorization and its
#' log-determinant cost 0.102 ms against 10.811 ms dense, and the full inverse
#' 3.280 ms against 25.000 ms, each route timed with its own factorization.
#' End to end that is 1.25x at 500 levels and 2.01x at 1000, the difference
#' between the operation and the fit being the lesson this package records
#' three times over: removing the dearer half leaves the cheaper one. Nothing
#' here asks which term produced the matrix; both quantities are read off the
#' matrix.
#'
#' Those figures are the ones measured with \pkg{Matrix}'s factorization
#' cache defeated. `Matrix::Cholesky` stores its result in the matrix's
#' `factors` slot, so a benchmark that refactorizes the same object
#' measures a cache hit, 0.004 ms against 0.102, and reports the sparse
#' route as three times better than it is. A fit never gets that hit,
#' the penalized matrix being a new one at every point.
#'
#' **The inverse stays dense whatever the factor is.** Its readers take full
#' matrix products against it, [block_leverage()] and the Hessian's pair loop
#' among them, so a dense inverse of a sparse matrix costs nothing sparsity
#' could have kept. What the sparse factor buys is the cost
#' of producing it.
#'
#' `NULL` is returned where the matrix is not positive definite, which is
#' the answer both callers already gave there.
#'
#' @param ctx A context, or `NULL`.
#' @param spec,design,coef,hyper The fallback arguments.
#'
#' @return A list with `K`, `inv` and `logdet`, or `NULL`.
#'
#' @seealso [pd_factor()]
#'
#' @keywords internal
ctx_penalized <- function(ctx, spec, design, coef, hyper, expected = FALSE) {
  build <- function() {
    H <- ctx_information(ctx, spec, design, coef, hyper, expected,
                         ctx_approx(ctx))
    S <- ctx_penalty(ctx, spec, design, coef, hyper)
    K <- H + S
    K <- pin_boundary(K)
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


#' Pin the Coordinates a Boundary Has Frozen
#'
#' @description
#' Replaces the row and column of every coordinate whose curvature is not
#' finite by the identity's, leaving the rest of the matrix alone.
#'
#' @details
#' A parameter that has reached the clamp its link keeps it strictly inside
#' leaves the family's curvature there at `NaN`, and one such entry is
#' enough to deny the whole matrix a Cholesky factor. The coordinate is not
#' one the criterion integrates over, being held at a boundary where its own
#' score is exactly zero, so the Laplace approximation wants nothing from it.
#' A unit diagonal contributes exactly that:
#' \eqn{\log\lvert K\rvert} gains \eqn{\log 1 = 0}, and
#' \eqn{K^{-1}g} returns that coordinate's own score, which is zero.
#'
#' The shape is preserved deliberately. Dropping the coordinate would change
#' the dimension of `K` and of its inverse, which some twenty consumers
#' index into by position; pinning says the same thing and breaks none of
#' them.
#'
#' @param K The penalized information.
#'
#' @return `K` with the frozen coordinates pinned, or `K` unchanged
#'   where there are none.
#'
#' @seealso [ctx_penalized()], [iwls_solve()]
#'
#' @keywords internal
pin_boundary <- function(K) {
  j <- boundary_coords(K)
  if (!length(j)) {
    attr(K, "held") <- 0L
    return(K)
  }
  K[j, ] <- 0
  K[, j] <- 0
  K[cbind(j, j)] <- 1
  attr(K, "held") <- length(j)
  K
}


#' Which Coordinates a Boundary Has Frozen
#'
#' @description
#' The positions whose curvature is not finite, as a parameter
#' sitting at the clamp its link keeps it strictly inside leaves behind.
#'
#' @details
#' The test is on the **diagonal**. A frozen coordinate makes its whole row
#' non-finite, cross terms included, so a test over whole columns marks its
#' neighbors too. Measured on a Student t whose \eqn{\nu} had reached
#' `double.xmax`, a column test held \eqn{\sigma} along with \eqn{\nu} and
#' left the fit exactly where it had been.
#'
#' @param K A square matrix, the information or the penalized information.
#'
#' @return An integer vector of positions, empty where there are none.
#'
#' @seealso [pin_boundary()], [iwls_solve()]
#'
#' @keywords internal
boundary_coords <- function(K) {
  d <- tryCatch(as.numeric(Matrix::diag(K)), error = function(e) NULL)
  if (is.null(d)) return(integer(0))
  which(!is.finite(d))
}



#' What the Marginal Criterion Can Resolve at This Fit
#'
#' @description
#' The size of the difference in the criterion that carries no information,
#' measured at the point. Nothing about it is assumed.
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
#' **The quadratic form alone is not enough**, and the reason is structural
#' before it is a matter of accuracy: the criterion carries
#' \eqn{-\tfrac{1}{2}\log\lvert K(\beta)\rvert}, which is not stationary in
#' \eqn{\beta}, so a mode error enters it at **first** order. Measured,
#' \eqn{\tfrac{1}{2} g' K^{-1} g} is right to one per cent at an inner
#' tolerance of `1e-4`, where the second-order term dominates, and
#' undershoots by 50 to 1000 times at `1e-6` and below, where the first
#' order does.
#'
#' Measured against the spread of the criterion at one hyperparameter reached
#' from six different warm starts, over four shapes and five inner tolerances,
#' the displaced reading tracks it across six orders of magnitude, from
#' `8.6e-2` down to `2.2e-7`, at a ratio between 0.05 and 0.99, and
#' separates shapes that a formula cannot: at an inner tolerance of `1e-6`
#' it reads `1.7e-4` on a random intercept over 500 levels against
#' `7.9e-8` on a gaussian smooth.
#'
#' It reads low by two to three times, the spread being a range over six paths
#' where this is one deviation, and that is the side to err on. A resolution
#' smaller than the truth leaves the search where it was; one larger stops a
#' healthy search short.
#'
#' @param st An evaluation's stored state, carrying `par`, `split`,
#'   `score`, `cf`, `hy`, `ctx` and `value`.
#' @param spec The specification.
#' @param design The design.
#' @param method An [OuterMethod()].
#' @param criterion_at A function of `(cf, hy, par, ctx)` returning what
#'   **this** search is running, a prediction-error criterion or a marginal
#'   one. It is passed in and never chosen here: reading the marginal
#'   criterion of a fit whose search is [aic()] answers for a quantity the
#'   search never sees, and the number that came back stopped two such fits
#'   short of their own optimum.
#'
#' @return A single positive number, or `NA_real_` where the pieces are
#'   not available.
#'
#' @seealso [outer_fit()], [optimizers7::armijo()]
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
  # ⚠️ A RESOLUTION IS A CORRECTION'S WORTH OF CRITERION, AND ONLY WHERE THE
  # CORRECTION IS SMALL. What is computed below is how far the criterion moves
  # when the coefficients are displaced by the mode error the inner score
  # implies, and that reading is a resolution only while the inner fit is AT a
  # mode. Where it is not, the same arithmetic returns something much larger,
  # and handing THAT to `crit_abs_obj()` converts a badly located mode into a
  # declaration of convergence.
  #
  # Measured on a hierarchical break-point model: the inner fit reports
  # convergence with a score of 247.8 -- against 3.6e-04 on a smooth and
  # 8.2e-06 on a random intercept -- the implied displacement is 20 in
  # coefficient units, and the criterion duly moves by 28, which is MORE THAN
  # THE WHOLE MOVEMENT OF THE CRITERION over the search (18.31) and four
  # orders above the 1.6e-03 the criterion's reproducibility measures
  # directly. Told that, lbfgs stopped after TWO evaluations reporting success
  # a hundred REML units short of what the same search reaches when it is told
  # nothing.
  #
  # The test is the displacement's own predicted decrease, `g'K^-1 g / 2`,
  # which is the amount by which the inner objective sits above its minimum
  # and is therefore in log-likelihood units rather than in the coefficients'.
  # It separates the two situations by nine orders: measured over whole fits
  # it reaches 2.8e-08 on a smooth and 1.4e-10 on a random intercept against
  # 20.9 to 22.9 on the break-point model, so the limit's exact value is not
  # what decides anything -- 1e-03 has five orders of room on either side.
  #
  # Refusing is NOT masking the inner fit's failure: the fit still reports
  # what it reports, the criterion is still read there, and what stops is
  # only the conversion of an unlocated mode into a stopping tolerance. The
  # reason is carried on the returned NA so that the trace can say it.
  quad <- tryCatch(0.5 * sum(st$score * db), error = function(e) NA_real_)
  if (!is.finite(quad) || quad > mode_error_limit()) {
    return(structure(NA_real_, mode_error = quad))
  }
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
#' @param ctx A context, or `NULL`.
#' @param pen The result of [ctx_penalized()].
#' @param basis The integrated subspace, or `NULL`.
#' @param expected Which information `pen` was built with, which is the
#'   cache's key: the projection is of that matrix, so one entry cannot serve
#'   both. Nothing reaches it today, a search holding one [OuterMethod()]
#'   throughout. The slot is keyed all the same, the twin defect in
#'   [ctx_penalized()] having been unreachable in exactly the same way until
#'   it was not.
#'
#' @return A square matrix, or `NULL`.
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
#' @param ctx A context, or `NULL`.
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
#' [block_leverage()], computed once and read by the gradient's
#' contraction and by every pair of the Hessian.
#'
#' @param ctx A context, or `NULL`.
#' @param design The design.
#' @param M The matrix the traces are taken against.
#' @param params,npar,offs The block bookkeeping.
#' @param threads How many threads the sparse route's kernel may use.
#'
#' @return A list of lists of vectors.
#'
#' @keywords internal
ctx_leverage <- function(ctx, design, M, params, npar, offs, threads = 1L) {
  if (is.null(ctx)) {
    return(block_leverage(design, M, params, npar, offs, threads))
  }
  if (is.null(ctx$leverage)) {
    ctx$leverage <- block_leverage(design, M, params, npar, offs, threads)
  }
  ctx$leverage
}


#' How the Penalized Matrix Moves With the Coefficients
#'
#' @description
#' The array [u_vector()] contracts against the leverage diagonal,
#' together with the builder that reads it: the third derivative of the
#' log-density where the criterion uses the observed information, and the
#' derivative of the expected information where it uses the expected one.
#'
#' @details
#' \eqn{K} enters the criterion through its determinant, so the gradient needs
#' \eqn{\partial K/\partial\beta}. With the observed information that is
#' \eqn{-\ell'''}, which every family carries in closed form. With the expected
#' one it is \eqn{-\partial\,\mathbb{E}[\ell'']/\partial\eta}, which is not
#' \eqn{-\mathbb{E}[\ell''']}: differentiating an expectation moves the measure
#' as well as the integrand, and the missing piece
#' \eqn{\mathbb{E}[\ell_{ab}\ell_{c}]} is a mixed moment no Bartlett identity
#' isolates: the third ties the symmetrized sum, never the single term.
#'
#' \pkg{distributions7} supplies it as
#' [distributions7::distrib_dexpected_hessian()]. The two arrays are
#' keyed differently, the observed one being symmetric in all three indices and
#' the expected one in its first two only, so the key builder travels with
#' the array and no consumer assumes one.
#'
#' @param ctx A context, or `NULL`.
#' @param spec,design,coef,hyper The fallback arguments.
#' @param method An [OuterMethod()].
#'
#' @return A list with `deriv` and `key`.
#'
#' @seealso [u_vector()], [outer_gradient_ok()]
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
#' @param ctx A context, or `NULL`.
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
#' @param ctx A context, or `NULL`.
#' @param spec,design,coef,hyper The fallback arguments.
#' @param order `3` or `4`.
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


#' How Far Above Its Minimum an Inner Fit May Sit and Still Report a Resolution
#'
#' @description
#' The limit [criterion_resolution()] puts on \eqn{g'K^{-1}g/2}, the
#' decrease the mode's own Newton correction predicts, before it refuses to
#' report a resolution at all.
#'
#' @details
#' The resolution is read by displacing the coefficients to where the inner
#' score says the mode is and asking how far the criterion moved. That is a
#' resolution while the displacement is a correction and something else
#' entirely once it is not: an inner fit that stopped far from a mode produces
#' a large displacement, a large movement, and a number that
#' [optimizers7::crit_abs_obj()] would read as license to stop.
#'
#' The quantity tested is in log-likelihood units, not in the coefficients',
#' so one limit serves every shape. Measured over whole fits it
#' reaches 2.8e-08 on a smooth and 1.4e-10 on a random intercept, against 20.9
#' to 22.9 on a hierarchical break-point model whose inner fit reports
#' convergence at a score of 247.8. Nine orders separate them, so the value
#' below has five orders of room on either side and is not what decides
#' anything.
#'
#' @return A single number.
#'
#' @seealso [criterion_resolution()]
#'
#' @keywords internal
mode_error_limit <- function() 1e-3


#' How Far Above Its Mode the Inner Fit Stopped
#'
#' @description
#' \eqn{\tfrac12 g'K^{-1}g}, the decrease the penalized likelihood's own
#' Newton correction predicts at the point the inner fit returned: how much
#' log-likelihood is still on the table there.
#'
#' @details
#' It answers the question **availability** asks, which is a different
#' question from the one the inner optimizer's flag answers. The flag says
#' whether a stopping rule fired. Availability asks whether the criterion, a
#' Laplace expansion at the mode, is valid at this point.
#'
#' The second is a matter of distance and has a natural scale, log-likelihood
#' units. The first is a boolean about a threshold on a score whose size
#' depends on the model.
#'
#' The two come apart on nearly every point. Measured on
#' `y ~ s(x) | sigma ~ s(z)`, of 38 inner fits during one search, 38 are at
#' their mode by this reading, between 1e-09 and 3e-09 against a limit of
#' 1e-03, and **four** report convergence. The
#' other 34 stopped on the objective-stall guard with the objective already
#' fixed to twelve significant digits and a score oscillating between 2.5e-06
#' and 3.3e-06, just above the absolute tolerance of 1e-06. Read as
#' unavailable, they made the outer line search backtrack eleven times per
#' iteration and accept a step of 0.0026 where the Newton step is 1.4, so the
#' search moved 0.005 in eta over 38 evaluations and stopped 4.0 criterion
#' units below the optimum its own gradient was correctly pointing at.
#'
#' It is used to add points and never to remove one: a run whose flag says
#' converged stays usable whatever this reads, so no model that fitted before
#' can stop fitting. That is also why it is not folded into the flag itself,
#' which `piano_stabilita.txt` section 13 measured and withdrew -- there
#' the flag was made stricter, and it cost a false negative on a good fit.
#'
#' @param ctx The evaluation context, so the penalized factorization is the
#'   one the criterion will read rather than a second copy.
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients the inner fit returned.
#' @param hyper The hyperparameters.
#' @param score The inner objective's gradient at those coefficients.
#' @param expected Whether the penalized information is the expected one.
#'
#' @return A single number, or `NA` where the penalized system could not
#'   be read there -- which is itself a reason to call the point unavailable.
#'
#' @seealso [mode_error_limit()], [criterion_resolution()]
#'
#' @keywords internal
inner_mode_error <- function(ctx, spec, design, coef, hyper, score,
                             expected = FALSE) {
  if (is.null(score) || !length(score) || !all(is.finite(score))) {
    return(NA_real_)
  }
  pen <- tryCatch(ctx_penalized(ctx, spec, design, coef, hyper, expected),
                  error = function(e) NULL)
  if (is.null(pen)) return(NA_real_)
  db <- tryCatch(as.numeric(as.matrix(pen$inv) %*% score),
                 error = function(e) NULL)
  if (is.null(db) || !all(is.finite(db))) return(NA_real_)
  v <- 0.5 * sum(score * db)
  if (is.finite(v)) v else NA_real_
}
