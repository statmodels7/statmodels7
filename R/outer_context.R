#' One Evaluation Point, Shared
#'
#' @description
#' Creates the environment that holds the quantities that the marginal
#' criterion, its gradient and its Hessian all read at one
#' \eqn{(\beta, \theta)}: the linear predictors, the information, the
#' penalty's Hessian, their sum and its factorization.
#'
#' @details
#' The criterion, its gradient and its Hessian read the same information, the
#' same penalty Hessian and the same factorization at one point. The context
#' holds them, so that each is computed once.
#'
#' The context is an environment, so an accessor fills it in place and later
#' readers find the quantity already there. Each quantity is computed on first
#' demand and never speculatively: the criterion alone does not need an
#' inverse, and a search that runs without an exact gradient does not compute
#' one.
#'
#' Passing `NULL` wherever a context is accepted computes every quantity
#' afresh at the point given.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param approx The approximation for the expected information.
#'
#' @return An environment carrying the point and, as they are requested, the
#'   quantities derived from it.
#'
#' @seealso [statmod_marginal()], [statmod_marginal_grad()],
#'   [statmod_marginal_hess()]
#'
#' @keywords internal
outer_context <- function(spec, design, coef, hyper, approx = "opg") {
  e <- new.env(parent = emptyenv())
  e$spec <- spec
  e$design <- design
  e$coef <- coef
  e$hyper <- hyper
  e$approx <- approx
  e
}


#' Check That a Context Belongs to the Point Read
#'
#' @description
#' Signals an error when a context is handed to a consumer that reads a
#' different point.
#'
#' @details
#' A cached information that belongs to another point yields a wrong gradient
#' without any other sign of failure. The comparison is over the coefficients
#' and the hyperparameters, so its cost is linear in their length, which is
#' negligible beside the \eqn{O(np^2)} assembly that it guards.
#'
#' @param ctx A context, or `NULL`.
#' @param coef The coefficients that the caller is reading.
#' @param hyper The hyperparameters that the caller is reading.
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
#' @param expected Whether the expected information is read.
#' @param approx The approximation for the expected information.
#'
#' @return A square matrix.
#'
#' @keywords internal
ctx_information <- function(ctx, spec, design, coef, hyper, expected,
                            approx = "opg") {
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
#' The non-finite entries are set to zero here, once.
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
#' Builds \eqn{K = H + S} with the criterion's own information, and its
#' inverse, which the gradient and the Hessian both read.
#'
#' @details
#' The storage is the matrix's own. \eqn{H} is sparse wherever the design
#' is, which covers a grouping indicator, a factor `by` and a `linpar` over
#' many levels. Where the sum \eqn{H + S} is large enough and sparse enough
#' ([worth_sparse()]) it is kept sparse and factorized as such. The decision
#' depends only on the matrix and not on which term produced it. The
#' penalized matrix is a new object at every point, so a factorization stored
#' on an earlier matrix is never reused.
#'
#' The inverse stays dense whatever the factor is. Its readers take full
#' matrix products against it, [block_leverage()] and the Hessian's pair loop
#' among them, so a dense inverse of a sparse matrix costs nothing that
#' sparsity could have kept. The sparse factor saves the cost of producing the
#' inverse.
#'
#' `NULL` is returned where the matrix is not positive definite. For the
#' mode's matrix (`laplace = FALSE`) the directions that [deficient_coords()]
#' reports as unidentified are pinned and the factorization is attempted once
#' more before `NULL` is returned.
#'
#' Where a penalty has a kink, there are two matrices. The criterion's
#' determinant leaves out every coordinate that a kinked penalty covers
#' ([laplace_pinned()]), which is the matrix that `laplace = TRUE` selects. The
#' default is the mode's own curvature, which leaves out the coordinates that
#' the marginal criterion holds or does not identify ([pinned_coords()]), the
#' boundary coordinates ([pin_boundary()]) and the coordinates that a kinked
#' penalty holds at zero ([zero_kinked()]): those stay at zero when a
#' hyperparameter moves, and a Newton correction of the mode does not move
#' them either. Where no penalty has a kink the two are one matrix and one
#' cache entry.
#'
#' @param ctx A context, or `NULL`.
#' @param spec,design,coef,hyper The fallback arguments.
#' @param expected Whether \eqn{H} is the expected information, as the
#'   criterion carries under `reml(hessian = "expected")`, or the observed one.
#'   The two give different matrices and keep separate cache entries.
#' @param laplace Whether the matrix is the criterion's, which leaves out the
#'   kinked coordinates, or the mode's.
#'
#' @return A list with `K`, `inv` and `logdet`, or `NULL`.
#'
#' @seealso [pd_factor()]
#'
#' @keywords internal
ctx_penalized <- function(ctx, spec, design, coef, hyper, expected = FALSE,
                          laplace = FALSE) {
  kinked <- length(kinked_coords(spec, design)) > 0L
  if (laplace && !kinked) laplace <- FALSE
  build <- function() {
    H <- ctx_information(ctx, spec, design, coef, hyper, expected,
                         ctx_approx(ctx))
    S <- ctx_penalty(ctx, spec, design, coef, hyper)
    K <- H + S
    # a coefficient the marginal criterion holds is not integrated, and is
    # left out of K exactly as a boundary coordinate is; so, in the
    # criterion's matrix, is a coefficient a kinked penalty covers
    K <- pin_boundary(K, if (laplace) laplace_pinned(spec, design) else
                           c(pinned_coords(spec, design),
                             if (kinked) zero_kinked(spec, design, coef)))
    fac <- pd_factor(K)
    # THE MODE'S MATRIX LEAVES OUT THE DIRECTIONS THE MODEL DOES NOT
    # IDENTIFY, as the criterion's determinant leaves them out (Giovanni,
    # 2026-10-07). A lasso codes every level of the first factor of its
    # formula, and two of its coefficients away from zero span the intercept:
    # the matrix was singular, the exact outer gradient came back NULL, the
    # search read it as zero and stopped at its start. Measured on
    # gamlss.data::rent, a gamma with s(Fl) + s(A) + lasso(~ B + H + L + loc):
    # the smoothing parameters stayed at 1; with the directions held they are
    # 5446.2 and 121.92, a central difference of the criterion gives 5446.4
    # and 121.93, and the battery's verdicts do not move.
    if (!isTRUE(fac$ok) && !laplace) {
      fl <- deficient_coords(K, attr(K, "held_at"))
      if (length(fl)) {
        K <- pin_boundary(K, c(pinned_coords(spec, design),
                               if (kinked) zero_kinked(spec, design, coef), fl))
        fac <- pd_factor(K)
      }
    }
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
  slot <- paste0(if (expected) "penalized_expected" else "penalized_observed",
                 if (laplace) "_laplace" else "")
  if (is.null(ctx[[slot]])) {
    # the miss is recorded so a second reader does not retry a factorization
    # that has already failed
    assign(slot, list(value = build()), envir = ctx)
  }
  ctx[[slot]]$value
}


#' Pin the Coordinates That a Boundary Has Frozen
#'
#' @description
#' Replaces the row and column of every coordinate whose curvature is not
#' finite, and of every coordinate named in `hold`, by the identity's, leaving
#' the rest of the matrix alone.
#'
#' @details
#' A parameter that has reached the clamp of its link, which keeps it strictly
#' inside its bounds, leaves the family's curvature there at `NaN`, and one
#' such entry is enough to deny the whole matrix a Cholesky factor. The
#' criterion does not integrate over the coordinate, because it is held at a
#' boundary where its own score is exactly zero, so it contributes nothing to
#' the Laplace approximation. A unit diagonal contributes exactly that:
#' \eqn{\log\lvert K\rvert} gains \eqn{\log 1 = 0}, and
#' \eqn{K^{-1}g} returns that coordinate's own score, which is zero.
#'
#' The dimension of `K` is preserved, because its consumers index `K` and its
#' inverse by position. Dropping the coordinate would change both dimensions.
#'
#' @param K The penalized information.
#' @param hold Integer positions to pin in addition to the coordinates with a
#'   non-finite diagonal (the coefficients that the marginal criterion holds,
#'   and the coordinates that a kinked penalty covers). The default is
#'   `integer(0)`.
#'
#' @return `K`, in which the pinned coordinates are replaced by the identity's
#'   rows and columns, carrying the attributes `held` (their number), `held_at`
#'   (their positions), `boundary_at` (the positions with a non-finite diagonal)
#'   and `hold_at` (the positions pinned only through `hold`). Where there are
#'   no pinned coordinates the matrix is otherwise unchanged and the positions
#'   are empty.
#'
#' @seealso [ctx_penalized()], [iwls_solve()]
#'
#' @keywords internal
pin_boundary <- function(K, hold = integer(0)) {
  b <- boundary_coords(K)
  j <- sort(unique(c(b, as.integer(hold))))
  if (!length(j)) {
    attr(K, "held") <- 0L
    attr(K, "held_at") <- integer(0)
    attr(K, "boundary_at") <- integer(0)
    attr(K, "hold_at") <- integer(0)
    return(K)
  }
  K[j, ] <- 0
  K[, j] <- 0
  K[cbind(j, j)] <- 1
  attr(K, "held") <- length(j)
  # the positions as well as the count: the criterion's gradient must know
  # WHICH rows are constant, see statmod_marginal_grad()
  attr(K, "held_at") <- j
  # and which of them are a BOUNDARY rather than a coefficient the marginal
  # criterion holds. Both are left out of the determinant, but only the first
  # is a constant of the criterion: a held coefficient moves the information
  # through the fitted predictors, so the criterion's derivative in it reads
  # its row of the contraction. See statmod_marginal_grad().
  attr(K, "boundary_at") <- b
  attr(K, "hold_at") <- setdiff(j, b)
  K
}


#' Which Coordinates a Boundary Has Frozen
#'
#' @description
#' Returns the positions whose curvature is not finite, as left by a parameter
#' sitting at the clamp of its link, which keeps it strictly inside its bounds.
#'
#' @details
#' The test is on the diagonal. A frozen coordinate makes its whole row
#' non-finite, including the cross terms, so a test over whole columns marks
#' its neighbors too.
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



#' What the Outer Criterion Can Resolve at This Fit
#'
#' @description
#' Computes the size of a difference in the criterion that carries no
#' information, from the fit at the point.
#'
#' @details
#' The criterion is read at the penalized mode, and the inner fit stops short
#' of that mode by whatever its rule allows. Writing \eqn{g} for the score at
#' which the inner fit stopped and \eqn{K = H + S} for the penalized
#' information, the mode is out by \eqn{\delta\beta = K^{-1} g}, and the
#' criterion read at \eqn{\hat\beta - \delta\beta} instead of at
#' \eqn{\hat\beta} differs by an amount that equals the difference between
#' two evaluations from different warm starts. One assembly of the criterion
#' at given coefficients gives it, with no refit.
#'
#' The criterion carries \eqn{-\tfrac{1}{2}\log\lvert K(\beta)\rvert}, which
#' is not stationary in \eqn{\beta}, so a mode error enters it at first order
#' and the quadratic form \eqn{\tfrac{1}{2} g' K^{-1} g} does not bound
#' the change by itself. The criterion is therefore
#' evaluated at the displaced coefficients.
#'
#' The reading estimates the spread of the criterion across warm starts. It is
#' usually below that spread, because the spread is a range over several paths
#' and the reading is one deviation. A resolution smaller than the true one
#' leaves the search where it was, and a larger one stops a search before it
#' has converged.
#'
#' The reading is `NA` where \eqn{g' K^{-1} g / 2} exceeds
#' [mode_error_limit()], because a displacement computed from an inner fit that
#' is far from its mode is not a resolution.
#'
#' @param st An evaluation's stored state, carrying `ok`, `par`, `split`,
#'   `score`, `cf`, `hy`, `ctx` and `value`.
#' @param spec The specification.
#' @param design The design.
#' @param method An [OuterMethod()].
#' @param criterion_at A function of `(cf, hy, par, ctx)` returning what the
#'   search is running, a prediction-error criterion or a marginal one. It is
#'   passed in and not chosen here, because the marginal criterion of a fit
#'   whose search runs [aic()] is a quantity that the search does not evaluate.
#'
#' @return A single positive number, or `NA_real_` where the pieces are not
#'   available. Where the mode error exceeds [mode_error_limit()], the `NA`
#'   carries the attribute `mode_error`, which holds that error.
#'
#' @seealso [outer_fit()], [optimizers7::armijo()]
#'
#' @keywords internal
criterion_resolution <- function(st, spec, design, method, criterion_at) {
  if (!isTRUE(st$ok) || is.null(st$score) || is.null(st$split)) {
    return(NA_real_)
  }
  if (!all(is.finite(st$score))) return(NA_real_)
  st$score <- free_of_kinks(st$score, spec, design, st$cf)
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



#' The Resolution That the Search Is Given, Out of the Readings Taken So Far
#'
#' @description
#' Returns the largest of the usable readings that [criterion_resolution()]
#' has returned during this search.
#'
#' @details
#' A reading is taken at every usable evaluation, and each one measures how
#' far the criterion moves when the mode is displaced by the error at which
#' that evaluation's inner fit stopped. The displacement varies across
#' evaluations by far more than the quantity that it represents does, so the
#' answer depends on which reading is kept.
#'
#' The readings are already filtered to acceptable modes, because
#' [criterion_resolution()] returns `NA` where the mode error exceeds
#' [mode_error_limit()], so every reading that reaches this function comes
#' from a mode that the layer accepts. The largest of them is therefore the
#' one to keep. A resolution has to bound how far the criterion can move
#' between two such modes, and not how little it moved at the evaluation
#' whose inner score happened to stop lowest.
#'
#' The inner score at which a Newton step stops varies by orders of magnitude
#' across evaluations, and a minimum over such a sequence reads the luckiest
#' one, which understates the resolution. A resolution far below the true one
#' leaves the line search backtracking for improvements that the criterion
#' cannot resolve, and the run ends by exhausting [outer_backtracks()] and not
#' on a stopping rule.
#'
#' @param x Every usable reading so far.
#'
#' @return A single positive number, or `NA_real_` where there is none.
#'
#' @seealso [criterion_resolution()], [outer_fit()]
#'
#' @keywords internal
resolution_summary <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(NA_real_)
  max(x)
}


#' The Matrix Against Which the Traces Are Taken
#'
#' @description
#' Returns \eqn{M}, which is the projected inverse when a subspace is projected
#' away and \eqn{K^{-1}} otherwise.
#'
#' @param ctx A context, or `NULL`.
#' @param pen The result of [ctx_penalized()].
#' @param basis The integrated subspace, or `NULL`.
#' @param expected The information with which `pen` was built, which is the
#'   cache's key: the projection is of that matrix, so one entry cannot serve
#'   both information types.
#'
#' @return A square matrix in which the rows and columns of the pinned
#'   coordinates are zero, or `NULL`.
#'
#' @keywords internal
ctx_trace_matrix <- function(ctx, pen, basis, expected = FALSE) {
  # ⚠️ A COORDINATE pin_boundary() HELD IS A CONSTANT OF THE CRITERION, so
  # every derivative of log|K| reads M with its row and column zeroed: the
  # pinned block is the identity whatever the coefficients are, and
  # differentiates to zero. Read through the identity, the trace multiplies
  # the family's derivative at a parameter on its link's clamp -- NaN there --
  # by one. Zeroed here, where M is made, so the gradient and the Hessian,
  # which share the leverage built from it, cannot disagree about it.
  held <- attr(pen$K, "held_at")
  unpin <- function(M) {
    if (is.null(M) || !length(held)) return(M)
    M <- as_dense(M)
    M[held, ] <- 0
    M[, held] <- 0
    M
  }
  if (is.null(basis)) return(unpin(pen$inv))
  if (!ncol(basis)) return(matrix(0, nrow(basis), nrow(basis)))
  build <- function() {
    # K may now be sparse, and the projection onto a dense basis is dense
    # whatever it was; the result is read as a full matrix by
    # block_leverage() and by the Hessian's pair loop, so it is a base matrix
    # here as it was before the storage became the matrix's own choice
    inner <- tryCatch(chol2inv(chol(as_dense(crossprod(basis,
                                                      pen$K %*% basis)))),
                      error = function(e) NULL)
    if (is.null(inner)) return(NULL)
    unpin(basis %*% inner %*% t(basis))
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
#' Returns the result of [block_leverage()], computed once and read by the
#' gradient's contraction and by every pair of the Hessian.
#'
#' @param ctx A context, or `NULL`.
#' @param design The design.
#' @param M The matrix against which the traces are taken.
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
#' Returns the array that [u_vector()] contracts against the leverage
#' diagonal, together with the builder that reads it: the third derivative of
#' the log-density where the criterion uses the observed information, and the
#' derivative of the expected information where it uses the expected one.
#'
#' @details
#' \eqn{K} enters the criterion through its determinant, so the gradient needs
#' \eqn{\partial K/\partial\beta}. With the observed information, that
#' derivative is \eqn{-\ell'''}, which every family carries in closed form.
#' With the expected
#' one it is \eqn{-\partial\,\mathbb{E}[\ell'']/\partial\eta}, which is not
#' \eqn{-\mathbb{E}[\ell''']}: differentiating an expectation moves the measure
#' as well as the integrand, so that
#' \eqn{\partial\mathbb{E}[\ell_{ab}]/\partial\eta_c = \mathbb{E}[\ell_{abc}] +
#' \mathbb{E}[\ell_{ab}\ell_{c}]}. The mixed moment
#' \eqn{\mathbb{E}[\ell_{ab}\ell_{c}]} is not available from a Bartlett
#' identity on its own, which relates only symmetrized sums of such moments.
#'
#' \pkg{distributions7} supplies it as
#' [distributions7::distrib_dexpected_hessian()]. The two arrays are
#' keyed differently, because the observed one is symmetric in all three
#' indices and the expected one in its first two only, so the key builder
#' travels with the array and a consumer does not assume one key.
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


#' How the Expected Information Moves With the Coefficients, Twice
#'
#' @description
#' Returns \eqn{\partial^2\,\mathbb{E}[\ell'']/\partial\eta^2}, the array that
#' the expected route's outer Hessian traces, twice contracted against the
#' leverage diagonal, together with its key builder.
#'
#' @details
#' It is read only on the expected route, where [outer_gradient_ok()] has
#' already checked [expected_deriv2_ok()]. The components are symmetric
#' separately in the information's pair and in the pair with respect to which
#' it is differentiated, so they are read through
#' [distributions7::d2expected_key()] and never through the observed route's
#' sorted quadruple.
#'
#' @param ctx A context, or `NULL`.
#' @param spec,design,coef,hyper The fallback arguments.
#'
#' @return A list with `deriv` and `key`.
#'
#' @seealso [ctx_kmove()], [statmod_marginal_hess()]
#'
#' @keywords internal
ctx_kmove2 <- function(ctx, spec, design, coef, hyper) {
  params <- spec@distrib@params
  build <- function() {
    th <- ctx_theta(ctx, spec, design, coef, hyper)
    distributions7::distrib_d2expected_hessian(spec@distrib, spec@response, th,
                                               scale = "link",
                                               threads = spec@threads)
  }
  v <- if (ctx_usable(ctx, coef, hyper)) {
    if (is.null(ctx$d2exp)) ctx$d2exp <- build()
    ctx$d2exp
  } else build()
  list(deriv = v,
       key = function(a, b, k, q)
         distributions7::d2expected_key(params, a, b, k, q))
}


#' The Approximation With Which a Context Was Built
#'
#' @param ctx A context, or `NULL`.
#'
#' @return A single string: the approximation that the context carries, or
#'   `"opg"` where `ctx` is `NULL` or carries none.
#'
#' @keywords internal
ctx_approx <- function(ctx) {
  if (is.null(ctx) || is.null(ctx$approx)) "opg" else ctx$approx
}


#' A Higher Derivative of the Log-Density at the Context's Point
#'
#' @description
#' Returns the third or fourth derivative on the link scale, which the
#' gradient's contraction and the Hessian's two contractions all read.
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
#' Returns the limit that [criterion_resolution()] places on
#' \eqn{g'K^{-1}g/2}, the decrease that the mode's own Newton correction
#' predicts. Above this limit, [criterion_resolution()] reports no resolution.
#'
#' @details
#' The resolution is read by displacing the coefficients to the mode implied
#' by the inner score and measuring how far the criterion moves. That is a
#' resolution while the displacement is a correction and something else
#' entirely once it is not: an inner fit that stopped far from a mode produces
#' a large displacement, a large movement, and a number that
#' [optimizers7::crit_abs_obj()] would read as license to stop.
#'
#' The quantity tested is in log-likelihood units and not in the units of the
#' coefficients, so one limit serves every model.
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
#' Returns \eqn{\tfrac12 g'K^{-1}g}, the decrease that the penalized
#' likelihood's own Newton correction predicts at the point that the inner fit
#' returned, which is the amount of log-likelihood still to be gained there.
#'
#' @details
#' Its value determines whether a point is available. The inner optimizer's
#' flag concerns a different question: it reports whether a
#' stopping rule fired. Availability is whether the criterion, a Laplace
#' expansion at the mode, is valid at this point.
#'
#' Availability is a matter of distance and has a natural scale, log-likelihood
#' units. The flag is a boolean about a threshold on a score whose size depends
#' on the model. The two can differ: an inner fit may stop on the
#' objective-stall guard while the score is slightly above the absolute
#' tolerance and the penalized likelihood is already at its mode to many
#' digits. Such a
#' point is available by this reading and unavailable by the flag.
#'
#' It is used to add points and not to remove any: a run whose flag reports
#' convergence stays usable whatever this reads.
#'
#' @param ctx The evaluation context, so the penalized factorization is the
#'   one that the criterion will read rather than a second copy.
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients that the inner fit returned.
#' @param hyper The hyperparameters.
#' @param score The inner objective's gradient at those coefficients.
#' @param expected Whether the penalized information is the expected one.
#'
#' @return A single number, or `NA` where the penalized system could not be
#'   read there, which is itself a reason to call the point unavailable.
#'
#' @seealso [mode_error_limit()], [criterion_resolution()]
#'
#' @keywords internal
inner_mode_error <- function(ctx, spec, design, coef, hyper, score,
                             expected = FALSE) {
  if (is.null(score) || !length(score) || !all(is.finite(score))) {
    return(NA_real_)
  }
  score <- free_of_kinks(score, spec, design, coef)
  pen <- tryCatch(ctx_penalized(ctx, spec, design, coef, hyper, expected),
                  error = function(e) NULL)
  if (is.null(pen)) return(NA_real_)
  db <- tryCatch(as.numeric(as.matrix(pen$inv) %*% score),
                 error = function(e) NULL)
  if (is.null(db) || !all(is.finite(db))) return(NA_real_)
  v <- 0.5 * sum(score * db)
  if (is.finite(v)) v else NA_real_
}
