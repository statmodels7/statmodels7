#' @include path.R
NULL

#' Fit a Separable Block by Coordinate Descent
#'
#' @description
#' Estimates one penalized block by cycling over its coefficients on the
#' working quadratic of its own equation, the other blocks held fixed.
#'
#' @details
#' **Why not the proximal method.** A proximal gradient step reads the
#' whole model: measured on 200 observations and 20 columns, one block fit made
#' 88 evaluations of the objective, 75 of the gradient and 83 of the operator,
#' each over every parameter of the distribution, and closed in 36 iterations
#' at 0.17 seconds. A coordinate descent reads the block's own columns and the
#' running residual instead and closes in six sweeps.
#'
#' **The working quadratic.** With \eqn{\eta} the equation's linear
#' predictor, \eqn{s_i} the score in it and \eqn{h_i} the information,
#' \eqn{-\ell} is \eqn{\frac12\sum_i h_i(z_i - \eta_i)^2} up to a constant with
#' \eqn{z = \eta + s/h}, which is the weighted least squares problem of
#' [iwls()] restricted to one equation. The other columns of that
#' equation enter as an offset. For a Gaussian response with an identity link
#' the quadratic is exact and one pass is the answer; otherwise the weights are
#' rebuilt and the sweeps repeated.
#'
#' **The penalty arrives as a table.** The coordinate update is the
#' penalty's own proximal operator at the step \eqn{1/v_j}, with
#' \eqn{v_j = \sum_i w_i x_{ij}^2}, and \eqn{v_j} does not move while the
#' working weights are held. The whole table is therefore built once per
#' weighted least squares iteration by
#' [penalties7::penalty_prox_spec()] and the compiled sweeps read it,
#' so the kernel names no family and a penalty that describes its operator gets
#' the compiled route without an edit here.
#'
#' **Screening.** Passing from one point of a path to the next, a
#' coordinate can be discarded when the gradient it had at the previous point
#' is below \eqn{2s_k - s_{k-1}}, with \eqn{s} the size of the kink: the
#' sequential strong rule of Tibshirani and others (2012), which assumes the
#' gradient moves at most as fast as the threshold does. That assumption is not
#' a theorem, so the rule can discard a coordinate that belongs in the fit, and
#' what makes the answer exact is the check afterwards: the gradient is read
#' over every column at the point reached, any discarded coordinate whose
#' gradient exceeds the kink is put back, and the fit is repeated. Without
#' the check the route would be occasionally wrong instead of occasionally
#' slow.
#'
#' **Which update.** The gradient is kept either as a residual, at
#' \eqn{O(n)} a visit, or as itself through
#' \eqn{g_j = (X'Wz)_j - \sum_k (X'WX)_{jk}\beta_k}, at \eqn{O(m)} a change
#' with the Gram columns cached as coordinates come alive. The second wins when
#' \eqn{n} is large next to the number of live coordinates and pays in
#' memory, so [coord_covariance()] decides it from the two sizes.
#'
#' @param obj The full objective, as [statmod_objective()] returns it. Read
#'   for the value at the point reached, not for its gradient.
#' @param beta The current stacked coefficients, a named list with one vector
#'   per distribution parameter. Every block but this one is held at these.
#' @param block One entry of `statmod_blocks()$sparse`: the equation, the
#'   term, its column positions and its penalty.
#' @param hyper The hyperparameters, per penalized term, held fixed here.
#' @param spec A [StatmodSpec()].
#' @param design The design, refreshed at `beta` if any term needs it.
#' @param expected `TRUE` for the expected information in the working
#'   weights, `FALSE` for the observed one.
#' @param approx How the expected information is approximated for a family
#'   with no closed form.
#' @param maxit The budget in weighted least squares iterations, each of
#'   which rebuilds the weights and runs the compiled sweeps to convergence.
#' @param tol The stopping tolerance on the relative change in the block's
#'   coefficients.
#' @param prev_kink The size of the kink at the previous point of a path, a
#'   single number, or `NULL` to cycle over every coordinate. Only the strong
#'   rule reads it.
#'
#' @return A list shaped like [sparse_fit()]'s, with the block's fitted
#'   coefficients, the sweep count and the gradient the kernel ended at.
#'   `NULL` where the route does not apply: the penalty describes no
#'   proximal table at the steps this block's curvature produces, or the
#'   working weights are unusable.
#'
#' @references
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. *Journal of
#' Statistical Software* 33(1), 1--22.
#'
#' Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor, J.
#' and Tibshirani, R. J. (2012). Strong rules for discarding predictors in
#' lasso-type problems. *Journal of the Royal Statistical Society, Series
#' B* 74(2), 245--266.
#'
#' @seealso [sparse_fit()],
#'   [penalties7::penalty_prox_spec()]
#'
#' @keywords internal
coord_fit <- function(obj, beta, block, hyper, spec, design, expected, approx,
                      maxit = 100L, tol = 1e-8, prev_kink = NULL) {
  p <- block$param
  d <- design[[p]]
  th <- as.list(hyper[[p]][[block$term]])
  cols <- block$cols
  # A HOLD CANNOT BE HONOURED HERE, AND MUST NOT BE IGNORED. This loop
  # updates one coordinate at a time from the running residual and never
  # reads `held_coef`, so a coefficient held inside a kinked block would be
  # moved by the sweep without a word -- measured, held at 3 and returned at
  # 2.7348. `statmod_restrict()` refuses such a coordinate before it reaches
  # here; this is what stops that refusal from being the only guard, and it
  # costs one property read where nothing is held.
  hf <- held_positions(spec, design, obj, beta)
  if (length(hf$where) && length(intersect(hf$where, block$index))) {
    stop(sprintf(paste0("A coefficient of '%s' is held at a value, and a ",
                        "kinked penalty is\n  fitted by a coordinate ",
                        "descent, which updates every coordinate of its ",
                        "own\n  block and cannot leave one where it was."),
                 block$term), call. = FALSE)
  }
  # whether the block is solved with its equation's intercept profiled out,
  # which is what keeps the alternation with that intercept from zig-zagging
  # on columns that are not centered: see coord_centers(). Where it is, the
  # intercept is set to the value the profiling implies after every sweep,
  # so the step is a joint step in the block and the intercept and it lowers
  # the objective as a step in the block alone would.
  ctr <- coord_centers(spec, design, p, obj, beta)
  i0 <- if (ctr) {
    obj$split(seq_along(beta))[[p]][parametric_intercept(spec, design, p)]
  } else NA_integer_
  # The block is kept in whatever storage it arrived in. A coordinate
  # descent reads one column at a time, so a compressed-column matrix is the
  # storage the method wants rather than one it tolerates, and the kernel
  # walks the stored nonzeros; densifying here was the last densification in
  # the chain and it is gone.
  # the built block, which answers the two questions asked before the loop --
  # how many columns there are, and whether this penalty has a table at all.
  # Inside the loop it is read again at the current coefficients.
  X <- coord_block_at(design, p, d$X, cols)
  if (!ncol(X)) return(NULL)
  other <- setdiff(seq_len(d$npar), cols)
  # whether anything in this model recomputes its own block as the
  # coefficients move, which is what makes the loop below read the design
  # again rather than the block it was handed
  rf <- attr(design, "refresh")
  moves <- !is.null(rf) && length(rf) > 0L
  # Does this penalty have a table AT ALL? That is a question about the
  # family, and it is asked at a step short enough not to answer a different
  # one: SCAD and MCP have no table past their convex region, the condition
  # being t < (a-1)/d^2 and t < gamma/d^2 under a diagonal map, so a probe at
  # t = 1 rejects a standardized penalty whose real steps are 1/sum(w x^2)
  # and orders of magnitude shorter. The step that will actually be used is
  # asked below, once the working weights are known.
  step0 <- rep(1e-10, ncol(X))
  if (is.null(penalties7::penalty_prox_spec(block$penalty, th, step0))) {
    return(NULL)
  }

  n <- spec@n_obs
  cur <- beta
  sweeps <- 0L
  prev <- NULL
  # the curvature the previous table used, which a scaled SCAD or MCP damps
  # towards the current step's: see coord_table_penalty()
  pen0 <- block$penalty
  c_prev <- if ("curv" %in% S7::prop_names(pen0) && length(pen0@curv)) {
    pen0@curv
  }
  c_gap <- 0
  for (it in seq_len(maxit)) {
    coef <- obj$split(cur)
    # THE BLOCK AT THE CURRENT COEFFICIENTS, not the block as it was built.
    # A term registering term_refresh() has neither property this route was
    # written for: its block is the Jacobian at the coefficients, so it moves
    # as they do, and what it contributes is X beta + adj rather than X beta.
    # Reading the block as it arrived and dropping adj solves a different
    # model and converges to a point that is not the mode -- measured on
    # nl(~ a * exp(-r * x), a ~ 0 + lasso(~grp)) at a held lambda small enough
    # that neither a lasso nor a ridge shrinks, log-likelihood -339.74 against
    # the ridge control's +155.45 and a rate of 0.22 against a truth of 0.70.
    #
    # It is refreshed once per SWEEP of this loop and never per coordinate:
    # the compiled descent exists because the design stands still while it
    # walks the columns, and statmod_design_at() chains from the state the
    # alternation commits, so the rescaling schedule of a break-point term
    # advances at the speed of the fit rather than of this loop. The result is
    # memoized on the coefficients, so statmod_eta() below reuses it.
    # asked only where something moves: with no refreshable term
    # statmod_design_at() returns the design it was given, and re-subsetting
    # the same columns every sweep would cost a copy of the block for nothing
    dd <- if (moves) statmod_design_at(spec, coef, design)[[p]] else d
    if (moves) X <- coord_block(dd$X, cols)
    ep <- statmod_eta(spec, design, coef)
    wq <- coord_working(spec, ep, coef, design, p, expected, approx)
    # Where the OBSERVED curvature of this equation is not positive at some
    # observation there are no working weights to read, and abandoning the
    # descent sends the block to the proximal route. The expected information
    # stands in instead, as iwls(hessian = "auto") does for the smooth block:
    # the mode is the same, only the weights of the step change. Measured on
    # a negbin2 lasso over mu and theta at n = 1500, 164 of 198 theta-block
    # fits fell back and the fit took 486.9 s; it takes 20.5 s, and where the
    # observed curvature is positive nothing moves, identical() on five other
    # shapes. Always reading the expected one instead costs 1.6x to 1.8x on a
    # gaussian or gamma with a modelled dispersion, so it stays a fallback.
    if (is.null(wq) && !expected)
      wq <- coord_working(spec, ep, coef, design, p, TRUE, approx)
    if (is.null(wq)) return(NULL)
    # the weighted column means the kernel centers with; none where the block
    # is solved as it stands
    mw <- if (ctr) as.numeric(xtv(X, wq$w, spec@threads)) / sum(wq$w) else
      numeric(0)
    off <- if (length(other))
      as.numeric(dd$X[, other, drop = FALSE] %*% coef[[p]][other]) else
      rep(0, n)
    # everything the working response carries that these columns do not: the
    # other columns, the equation's offset, and what the term contributes
    # beyond its block, which statmod_eta() has already put into `z`
    adj <- if (is.null(dd$adj)) 0 else dd$adj
    z <- wq$z - off - coord_offset(spec, p, n) - adj
    # THE WORKING PROBLEM HAS NOT MOVED. Where the previous sweep's weights
    # are the same and its working response differs only by rounding -- a
    # gaussian mean, whose working response is y whatever the coefficients --
    # the previous call already returned this problem's solution, the strong
    # rule's discards checked against the kink. Calling again only confirmed
    # it: measured on a gaussian lasso path, 62 of 64 fits made that call,
    # it moved the coefficients by at most 1.7e-10 against a tolerance of
    # 1e-8, and it was half of the descent's time. So the loop ends, at a
    # point that differs from the confirmed one by that much. A block that
    # moves with its coefficients changes the problem through the design,
    # where the weights and the response cannot see it, and is never skipped.
    # A centered block does not see a constant added to the response, and the
    # intercept moved at the end of the previous sweep adds exactly that, so
    # the comparison is made on the centered response.
    zc <- if (ctr) z - sum(wq$w * z) / sum(wq$w) else z
    if (!moves && !is.null(prev) && identical(wq$w, prev$w) &&
        max(abs(zc - prev$z)) <= 64 * .Machine$double.eps * max(abs(z)) &&
        c_gap <= 1e-10) {
      break
    }
    prev <- list(w = wq$w, z = zc)
    # The column curvatures are read only on the columns the descent visits,
    # which the strong rule keeps to a few of the block. What the full vector
    # was also for is the check that every one is finite and positive, and
    # where the design stands still that is decided from its cached column
    # norms and the range of the weights, a sufficient condition; outside it
    # the whole vector is computed and checked as before.
    v <- rep(NA_real_, ncol(X))
    if (!coord_curv_bounded(X, wq$w, design, p, cols, moves)) {
      v <- wxsq(X, wq$w, spec@threads)
      if (any(!is.finite(v)) || any(v <= 0)) return(NULL)
    }
    # a centered coordinate's curvature is its centered sum of squares, read
    # on the columns the descent visits like the uncentered one
    if (ctr) v <- rep(NA_real_, ncol(X))
    b0 <- coef[[p]][cols]
    # The size of the kink COORDINATE BY COORDINATE. Under a diagonal map,
    # which is what `standardize` writes, coordinate j's kink is lambda
    # |d_j|, and reading the first coordinate's for all of them screened and
    # rechecked every other coordinate against a threshold that was not its
    # own. The previous point's kinks are this point's scaled by the ratio
    # the path recorded, every coordinate's kink moving with lambda alike.
    s_now <- coord_kinks(block$penalty, th)
    s_prev <- if (!is.null(prev_kink) && s_now[[1L]] > 0) {
      s_now * (prev_kink / s_now[[1L]])
    }
    keep <- coord_screen(X, wq$w, z, b0, s_now, s_prev, spec@threads,
                         center = ctr)

    repeat {
      need <- keep[is.na(v[keep])]
      if (length(need)) {
        v[need] <- if (ctr) coord_colsq(X, wq$w, need, mw) else
          coord_curv(X, wq$w, need, spec@threads)
      }
      # THE TABLE IS BUILT FOR THE KEPT COORDINATES, because the kernel reads
      # its row a for coordinate keep[a]. Built from the whole penalty with the
      # kept coordinates' steps, row a carried coordinate a's map entry and
      # curvature: under a strong rule that screened anything out, a
      # standardized lasso soft-thresholded coordinate keep[a] at another
      # coordinate's scale, and a SCAD or MCP read a curvature 13 to 28 per
      # cent off, enough to fail the step condition. See coord_table_penalty().
      tp <- coord_table_penalty(block$penalty, keep, v, c_prev)
      if (!is.null(c_prev) && length(tp@curv) == length(keep)) {
        # how far the damped curvature still is from the current step's, in
        # log units: where the working problem has not moved -- a gaussian
        # mean at a held scale -- the damping still has to reach it before
        # the loop may stop, or the fit is optimal for a curvature it does
        # not report (measured: KKT residuals of 1e-3 at the fit's curvature)
        tgt <- v[keep]
        if (!is.null(pen0@map)) {
          tgt <- tgt / as.numeric(Matrix::diag(pen0@map))[keep]^2
        }
        c_gap <- max(abs(log(tp@curv / tgt)))
        c_prev[keep] <- tp@curv
      }
      tab <- penalties7::penalty_prox_spec(tp, th, 1 / v[keep])
      if (is.null(tab)) return(NULL)
      out <- coord_call(X, z, wq$w, b0, tab, as.integer(keep - 1L), tol,
                        coord_covariance(n, length(keep)), means = mw)
      sweeps <- sweeps + as.integer(out$sweeps)
      # a strong rule is a heuristic: a coordinate it discarded whose gradient
      # exceeds the kink belongs in the fit, and only this comparison makes
      # the answer exact rather than usually right. With nothing discarded
      # there is nothing to check and the kernel does not compute it.
      if (length(keep) == ncol(X)) break
      back <- setdiff(which(abs(out$grad) > s_now * (1 + 1e-10)), keep)
      if (!length(back)) break
      keep <- sort(c(keep, back))
    }
    moved <- max(abs(out$beta - b0))
    cur[block$index] <- out$beta
    # The intercept the profiling implies: the weighted mean of the working
    # response net of the block, which `z` carries with the intercept's old
    # value already taken off.
    if (ctr) cur[i0] <- cur[i0] + sum(wq$w * z) / sum(wq$w) - sum(mw * out$beta)
    if (moved < tol && c_gap <= 1e-10) break
  }
  list(par = cur, value = obj$fn(cur), converged = TRUE,
       iterations = sweeps, method = "coordinate descent")
}


#' Which Coordinates a Path Point Has to Visit
#'
#' @description
#' The sequential strong rule: a coordinate whose gradient at the previous
#' point of the path is below \eqn{2s_k - s_{k-1}} is left out of the sweeps.
#'
#' @details
#' The rule rests on the gradient moving no faster than the threshold does.
#' That is an assumption and not a bound, so the rule screens without
#' proving, and the caller checks what it discarded. With no previous point
#' there is nothing to screen against and every coordinate is visited.
#'
#' @param X The block's own columns, `n x p`, dense or `dgCMatrix`.
#' @param w The working weights, length `n`.
#' @param z The working response, length `n`.
#' @param beta The block's coefficients at the previous point of the path,
#'   length `p`.
#' @param s_now The size of the kink at this point, one number per column
#'   (a single number is recycled).
#' @param s_prev The size of the kink at the previous point, one number per
#'   column, or `NULL` when there is no previous point.
#' @param threads The thread count the gradient read may use, a plain
#'   integer.
#'
#' @return An integer vector of one-based column indices to visit, in
#'   ascending order. Every column when `s_prev` is `NULL` or either kink
#'   size is not usable. A coordinate already away from zero is always kept,
#'   whatever its gradient, so the rule can only ever add coordinates to the
#'   active set. Never empty: where the test discards everything, the column
#'   with the largest gradient is kept.
#'
#' @seealso [coord_fit()]
#'
#' @keywords internal
coord_screen <- function(X, w, z, beta, s_now, s_prev, threads = 1L,
                         center = FALSE) {
  p <- ncol(X)
  # With no previous point there is nothing to screen against, and the rule in
  # its global form -- the reference being the kink that empties the block --
  # discards nothing at any smoothing parameter worth fitting at, since
  # 2s - max|g| is negative there. Measured at a single value it cost one
  # crossprod and saved no coordinate, so it is not attempted.
  if (is.null(s_prev) || any(!is.finite(s_prev)) || any(!is.finite(s_now)) ||
      any(s_now <= 0)) {
    return(seq_len(p))
  }
  r <- z - x_times_b(X, beta)
  # with the block centered the gradient is that of the centered columns,
  # which is X'W r once r is taken to weighted mean zero
  if (center) r <- r - sum(w * r) / sum(w)
  g <- abs(xtv(X, w * r, threads))
  keep <- which(g >= 2 * s_now - s_prev | beta != 0)
  if (!length(keep)) keep <- which.max(g)
  keep
}


#' The Penalized Block, in the Storage It Arrived In
#'
#' @description
#' Slices a penalized term's columns out of its equation's design without
#' densifying a sparse one, and normalizes a \pkg{Matrix} to the
#' compressed-column class the kernel reads.
#'
#' @details
#' A dense slice of a base matrix is returned as it is. Any \pkg{Matrix} is
#' carried to `dgCMatrix`: the general compressed-column form is the
#' one whose slots the kernel walks, and a symmetric or triangular
#' compression would describe the same entries differently. A dense
#' \pkg{Matrix} class is materialized as a base matrix instead, there being
#' nothing to save.
#'
#' @param X The equation's design, dense or any \pkg{Matrix} class.
#' @param cols The term's column positions within it, an integer vector.
#'
#' @return A base numeric matrix with `length(cols)` columns when `X` is
#'   dense or a dense \pkg{Matrix} class, and a `dgCMatrix` when `X` is
#'   sparse.
#'
#' @seealso [coord_call()], [coord_fit()]
#'
#' @keywords internal
coord_block <- function(X, cols) {
  B <- X[, cols, drop = FALSE]
  if (!isS4(B)) return(B)
  if (methods::is(B, "sparseMatrix")) {
    return(methods::as(methods::as(B, "generalMatrix"), "CsparseMatrix"))
  }
  as.matrix(B)
}

#' Run the Compiled Coordinate Descent on Either Storage
#'
#' @description
#' Sends the block to the dense kernel or to the sparse one, taking a
#' `dgCMatrix` apart into the slots the second reads.
#'
#' @details
#' The two kernels are one algorithm instantiated twice over a column
#' accessor, and they agree bit for bit. The arithmetic licenses that:
#' skipping a structural zero omits an addition of zero, which is exact. It
#' is the one place in this toolkit where an identity assertion over compiled
#' floating point is correct.
#'
#' The `dgCMatrix` is taken apart in R, so the compiled code needs no
#' dependency on the \pkg{Matrix} package's C API.
#'
#' @param X The block, `n x p`, dense or `dgCMatrix`.
#' @param z,w The working response and weights, each of length `n`.
#' @param b0 The starting coefficients, length `p`.
#' @param tab The piecewise linear proximal table, as
#'   [penalties7::penalty_prox_spec()] returns it.
#' @param screen The **zero-based** positions the strong rule kept, for the
#'   C++ indexing.
#' @param tol The stopping tolerance on the largest coefficient change of a
#'   sweep.
#' @param covariance `TRUE` to hold the gradient itself and cache Gram
#'   columns, `FALSE` to hold the running residual. [coord_covariance()]
#'   decides.
#'
#' @return The kernel's list of three: `beta` (the fitted coefficients,
#'   length `p`), `sweeps` (how many passes it took) and `grad` (the gradient
#'   at the point reached, length `p`).
#'
#' @seealso [coord_block()]
#'
#' @keywords internal
coord_call <- function(X, z, w, b0, tab, screen, tol, covariance,
                       means = numeric(0)) {
  if (isS4(X)) {
    return(coord_descent_sparse(X@i, X@p, X@x, nrow(X), ncol(X), z, w, b0,
                                tab$cut, tab$slope, tab$icept, screen, 500L,
                                tol, covariance, means))
  }
  coord_descent(X, z, w, b0, tab$cut, tab$slope, tab$icept, screen, 500L,
                tol, covariance, means)
}


#' The Size of the Kink in Each Coordinate
#'
#' @description
#' Returns, for every coordinate of a kinked penalty, the jump of its
#' derivative at the kink, which is the threshold a coordinate's gradient is
#' compared with to decide whether it can stay at zero.
#'
#' @details
#' [kink_scale()] reads the first coordinate only, which is what a path needs
#' to place its grid. A coordinate descent needs every coordinate's: under a
#' diagonal map \eqn{D}, which is what `standardize` writes, coordinate
#' \eqn{j}'s kink is \eqn{\lambda\lvert d_j\rvert}. The jump is measured as in
#' [kink_scale()], by a Richardson extrapolation of the one-sided derivatives,
#' on every coordinate at once. A scaled SCAD or MCP has its kink at
#' \eqn{\lambda\lvert d_j\rvert} whatever its curvature, the scaling leaving
#' the slope at the origin where it was.
#'
#' @param pen A kinked penalty.
#' @param theta Its hyperparameters.
#' @param eps The step of the one-sided derivatives.
#'
#' @return A numeric vector with one entry per coordinate of `pen`, all zero
#'   where the penalty reports no finite kink. Its first entry is
#'   [kink_scale()]'s answer.
#'
#' @seealso [kink_scale()], [coord_screen()], [coord_fit()]
#'
#' @keywords internal
coord_kinks <- function(pen, theta, eps = 1e-4) {
  th <- as.list(theta)
  p <- max(1L, as.integer(pen@n_coef))
  k <- penalties7::penalty_kinks(pen, th)
  k <- k[is.finite(k)]
  if (!length(k)) return(rep(0, p))
  at <- function(h) {
    up <- penalties7::penalty_gradient(pen, rep(k[[1L]] + h, p), th)
    dn <- penalties7::penalty_gradient(pen, rep(k[[1L]] - h, p), th)
    as.numeric(up - dn) / 2
  }
  2 * at(eps / 2) - at(eps)
}


#' The Penalty a Coordinate Descent Builds Its Table From
#'
#' @description
#' Restricts a kinked penalty to the coordinates the strong rule kept, and
#' writes into a scaled SCAD or MCP the curvature of the current step.
#'
#' @details
#' The compiled descent reads row \eqn{a} of the table for coordinate
#' \eqn{k_a}, the \eqn{a}-th kept one. A table built from the whole penalty
#' pairs row \eqn{a} with coordinate \eqn{a} instead, so the restriction is
#' what makes the two agree: every per-coordinate property is subset, the
#' map's diagonal and the curvature among them. It is done only where a
#' coordinate was screened out or a curvature is written, and a penalty under
#' a map that is not diagonal is returned as it stands, having no table.
#'
#' A SCAD or MCP carrying a curvature is scaled SELF-CONSISTENTLY: at the
#' point the descent settles at, its curvature is the one the step is taken
#' with, \eqn{c_j = v_j}, the (centered) weighted sum of squares of the column
#' at the current working weights, divided by \eqn{d_j^2} under a diagonal map,
#' the penalty being written on \eqn{u = D\beta}. On the way there it is
#' DAMPED, \eqn{c_j \leftarrow \sqrt{c_j^{\mathrm{prev}} v_j}}: taken
#' undamped, a coefficient between the two knees of a logistic SCAD
#' alternated between 1.009 and 1.496 with its curvature between 21.8 and
#' 30.7 and never settled, 3 fits of 16 not converging on a strong-effect
#' probe; damped, all 16 converge with the KKT conditions at the fit's
#' curvature met to \eqn{5 \times 10^{-8}}. The scaled step is
#' \eqn{\sqrt{c^{\mathrm{prev}}_j / v_j}}, one where the weights have not
#' moved, so the step condition \eqn{t c_j < a - 1} (SCAD) or
#' \eqn{t c_j < \gamma} (MCP) holds unless they moved by a factor of
#' \eqn{(a-1)^2}. [statmod_curv()] writes the undamped quantity into the
#' specification at every pass of the alternation, so the objective and the
#' degrees of freedom read it too.
#'
#' @param pen The block's penalty.
#' @param keep The kept coordinates, one-based.
#' @param v The column curvatures, one per coordinate of the block, with at
#'   least the kept ones filled in.
#' @param prev `NULL`, or the curvature the previous table used, on the scale
#'   of \eqn{u}, one per coordinate of the block; the new one is its geometric
#'   mean with the current step's.
#'
#' @return A penalty over `length(keep)` coordinates, or `pen` itself where
#'   nothing needs to change.
#'
#' @seealso [coord_fit()], [statmod_curv()]
#'
#' @keywords internal
coord_table_penalty <- function(pen, keep, v, prev = NULL) {
  if (!is.null(pen@map) && !methods::is(pen@map, "diagonalMatrix")) {
    return(pen)
  }
  q <- as.integer(pen@n_coef)
  sub <- !(length(keep) == q && all(keep == seq_len(q)))
  scaled <- "curv" %in% S7::prop_names(pen) && length(pen@curv) > 0L
  props <- list()
  if (sub) {
    props$n_coef <- length(keep)
    if (!is.null(pen@map)) props$map <- pen@map[keep, keep, drop = FALSE]
  }
  if (scaled) {
    cj <- v[keep]
    if (!is.null(pen@map)) {
      cj <- cj / as.numeric(Matrix::diag(pen@map))[keep]^2
    }
    if (!is.null(prev)) cj <- sqrt(prev[keep] * cj)
    if (all(is.finite(cj)) && all(cj > 0)) props$curv <- cj
    else if (sub) props$curv <- pen@curv[keep]
  }
  if (!length(props)) return(pen)
  do.call(S7::set_props, c(list(pen), props))
}


#' Record Where a Path Has Just Been
#'
#' @description
#' Writes the size of each kinked penalty's kink at the given hyperparameters
#' onto its block, so that the next point of a path can screen against it.
#'
#' @details
#' The previous point travels on the blocks themselves, sparing every layer
#' between the path and the descent an argument. Where a penalty was a moment
#' ago is a property of its block, and the path rebuilds the blocks at each
#' point in any case.
#'
#' @param blocks The blocks, as [statmod_blocks()] returns them, with a
#'   `sparse` list of the kinked entries.
#' @param hyper The hyperparameters at the point just fitted.
#'
#' @return `blocks`, with each entry of its `sparse` list carrying a
#'   `prev_kink` element: the size of that penalty's kink at `hyper`. The
#'   `smooth` list is untouched.
#'
#' @seealso [coord_screen()], [statmod_path()]
#'
#' @keywords internal
blocks_at_kink <- function(blocks, hyper) {
  for (i in seq_along(blocks$sparse)) {
    b <- blocks$sparse[[i]]
    th <- as.list(hyper[[b$param]][[b$term]])
    s <- tryCatch(kink_scale(b$penalty, th), error = function(e) NA_real_)
    blocks$sparse[[i]]$prev_kink <- if (is.finite(s) && s > 0) s else NULL
  }
  blocks
}


#' Which Way of Holding the Gradient Is Cheaper
#'
#' @description
#' Chooses how the compiled sweeps keep the gradient: `TRUE` for the
#' covariance form, which holds the gradient itself and caches columns of
#' \eqn{X'WX}, and `FALSE` for the running residual. The test is
#' `m <= 32 && n > 8 * m`.
#'
#' @details
#' The covariance form replaces an \eqn{O(n)} read of the gradient with an
#' \eqn{O(m)} one, and pays for it by building a column of \eqn{X'WX} at
#' \eqn{O(nm)} the first time a coordinate moves off zero. It is worth having
#' only while \eqn{m} is small next to \eqn{n}, and the two conditions say
#' exactly that.
#'
#' The measurement is unambiguous in the other direction: at 5000
#' observations with 200 columns and nothing screened away, the covariance
#' form cost 70 milliseconds against 55 for the residual, the Gram columns
#' being dearer than the residual passes they replaced.
#'
#' @param n The number of observations.
#' @param m How many coordinates the strong rule left to visit.
#'
#' @return A single logical.
#'
#' @seealso [coord_call()], which passes the answer to the kernel,
#'   [coord_screen()] for `m`.
#'
#' @keywords internal
coord_covariance <- function(n, m) {
  m <= 32L && n > 8L * m
}


#' The Working Response and Weights of One Equation
#'
#' @description
#' \eqn{h_i} and \eqn{z_i = \eta_i + s_i/h_i}, the weighted least squares
#' problem the log-likelihood is locally.
#'
#' @details
#' For a Gaussian response on the identity link the quadratic is exact and
#' one pass answers the problem. Elsewhere it is the local approximation a
#' scoring step works on, and the weights are rebuilt at each iteration.
#'
#' Where the design carries a structural term, the predictor of an equation
#' includes a part that belongs to that term, and no column of a block
#' carries it. The working response is then built on the static predictor,
#' with \eqn{s_i} the score of the model from [statmod_score_obs()], so the
#' quadratic has the gradient of the log-likelihood at the current
#' coefficients. The weights stay the family's curvature at the whole
#' predictor. Another positive weight would change the path of the descent
#' and would leave its fixed point where it is.
#'
#' @param spec A [StatmodSpec()].
#' @param ep The linear predictors and the parameters they imply, as
#'   [statmod_eta()] returns them.
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#' @param p Which distribution parameter's equation, a string.
#' @param expected `TRUE` for the expected information, `FALSE` for the
#'   observed one.
#' @param approx How the expected information is approximated for a family
#'   with no closed form.
#'
#' @return A list with `w` and `z`, each a numeric vector of length
#'   `spec@n_obs`. `NULL` where the curvature is not usable, which is any
#'   non-finite or non-positive \eqn{h_i}; the observed information can
#'   produce both far from the optimum.
#'
#' @keywords internal
coord_working <- function(spec, ep, coef, design, p, expected, approx) {
  n <- spec@n_obs
  params <- spec@distrib@params
  a <- match(p, params)
  # A STRUCTURAL TERM ADDS TO THE PREDICTOR A PART THAT NO COLUMN OF THE
  # BLOCK CARRIES: a filter's level, or the posterior mean of the shifts of a
  # likelihood mixed over latent states. The working response is then built
  # on the static predictor, from the score of the model rather than the
  # family's score at the whole predictor, and statmod_score_obs() supplies
  # that score through the filter's reverse recursion or Fisher's identity.
  # Built on the whole predictor, as it was until 0.154.0, it carried that
  # part and the descent solved a different problem: with a lasso held at a
  # fixed lambda, an active coefficient's score sat 38, 27 and 57 from its
  # KKT value beside gas(1, 1), regime(2) and a marginal jump(), and beside
  # regime(2) with a Poisson response the fit diverged to a log-likelihood
  # of -129223. The weights stay the family's curvature at the whole
  # predictor: any positive weight leaves the fixed point where the score
  # puts it.
  if (is.null(ep$eta_static)) {
    g <- distributions7::distrib_gradient(spec@distrib, spec@response,
                                          ep$theta, scale = "link",
                                          threads = spec@threads)
    s <- spec@weights * rep_len(g[[p]], n)
    e0 <- rep_len(ep$eta[[p]], n)
  } else {
    s <- statmod_score_obs(spec, coef, design)[[p]]
    e0 <- rep_len(ep$eta_static[[p]], n)
  }
  # the one diagonal entry, read where info_blocks() would have filed it in
  # an n x K x K array this function then took one slice of
  H <- statmod_family_hessian(spec, ep$theta, expected, approx)
  h <- -spec@weights * rep_len(H[[hess_key(params, a, a)]], n)
  if (any(!is.finite(h)) || any(h <= 0) || any(!is.finite(s))) return(NULL)
  list(w = h, z = e0 + s / h)
}


#' The Offset of One Equation
#'
#' @description
#' Returns the offset of one distribution parameter's equation, evaluated in
#' the fitting data and recycled to the sample size. An equation with no
#' offset gets zeros, so the caller adds the result unconditionally.
#'
#' @details
#' The offsets are stored on the specification as evaluated vectors, one per
#' parameter, with `NULL` for an equation that has none. This turns that
#' `NULL` into zeros at the point of use, so no caller has to test for it.
#'
#' An offset shorter than `n` is recycled with [rep_len()], so a single
#' number is a constant offset. That is the shape a caller writing
#' `offsets = list(mu = log(2))` gets.
#'
#' @param spec A [StatmodSpec()], read for its `offsets` list.
#' @param p Which distribution parameter, a string naming one of the
#'   family's.
#' @param n The number of observations to recycle to.
#'
#' @return A numeric vector of length `n`: the offset, or zeros where the
#'   equation has none.
#'
#' @seealso [eval_offsets()], which evaluates the expressions this reads,
#'   [statmod()] for the `offsets` argument.
#'
#' @keywords internal
coord_offset <- function(spec, p, n) {
  o <- spec@offsets[[p]]
  if (is.null(o) || !length(o)) return(rep(0, n))
  rep_len(as.numeric(o), n)
}


#' Whether Every Column Curvature Is Finite and Positive, Without Them
#'
#' @description
#' `coord_curv_bounded()` decides whether every column curvature
#' \eqn{v_j = \sum_i w_i x_{ij}^2} of a block would be finite and positive,
#' from the block's column norms and the range of the weights, so that
#' [coord_fit()] need not compute the curvatures of columns it never visits.
#' `coord_curv()` computes them on the columns asked for.
#'
#' @details
#' The condition is sufficient and not necessary. With the weights finite and
#' positive, which [coord_working()] has already checked, a column's
#' curvature is at least its largest term, which is at least
#' \eqn{\min_i w_i\, \lVert x_j\rVert^2 / n}, and at most
#' \eqn{n \max_i w_i \max_i x_{ij}^2}. So where every column norm is
#' positive and \eqn{\min w \min_j \lVert x_j\rVert^2 > 10^{-200}} and
#' \eqn{\max w \max_j \lVert x_j\rVert^2 < 10^{200}}, no curvature can be
#' zero, underflow or overflow. Where any of this fails, `FALSE` sends the
#' caller to the full computation and its own check, so the outcome is the
#' one the full vector would have given.
#'
#' The norms are cached in the design's `eta_memo` environment, keyed by the
#' equation and the block's columns. A design whose blocks move with the
#' coefficients carries no such environment, and there the answer is always
#' `FALSE`.
#'
#' `coord_curv()` takes the route [wxsq()] would take for the whole block:
#' each curvature is its own column's sum, so it is the same number either
#' way.
#'
#' @param X The block, dense or `dgCMatrix`.
#' @param w The working weights, finite and positive.
#' @param design The design the block was read from.
#' @param p The equation, a parameter name.
#' @param cols The block's column positions within the equation.
#' @param moves Whether any term recomputes its block with the coefficients.
#' @param k Integer column positions within `X`.
#' @param threads The thread count [wxsq()] would be given.
#'
#' @return `coord_curv_bounded()` gives a single logical. `coord_curv()`
#'   gives a numeric vector of `length(k)` curvatures.
#'
#' @seealso [coord_fit()], [wxsq()]
#'
#' @keywords internal
coord_curv_bounded <- function(X, w, design, p, cols, moves) {
  mm <- attr(design, "eta_memo")
  if (moves || is.null(mm)) return(FALSE)
  key <- paste0("colsq:", p)
  cq <- mm[[key]]
  if (is.null(cq) || !identical(cq$cols, cols)) {
    s <- if (isS4(X)) Matrix::colSums(X^2) else colSums(X^2)
    cq <- list(cols = cols, ok = length(s) > 0L && all(is.finite(s)) &&
                 all(s > 0), lo = min(s), hi = max(s))
    assign(key, cq, envir = mm)
  }
  if (!isTRUE(cq$ok)) return(FALSE)
  lo <- min(w) * cq$lo
  hi <- max(w) * cq$hi
  is.finite(lo) && is.finite(hi) && lo > 1e-200 && hi < 1e200
}

#' @rdname coord_curv_bounded
#' @keywords internal
coord_curv <- function(X, w, k, threads = 1L) {
  big <- threads > 1L && is.matrix(X) && length(w) == nrow(X) &&
    as.double(nrow(X)) * ncol(X) >= 2e5
  Xk <- X[, k, drop = FALSE]
  if (big) return(wxsq_cpp(Xk, w, threads))
  as.numeric(crossprod(w, Xk^2))
}


#' The Penalized Block, Kept on the Design
#'
#' @description
#' [coord_block()] of an equation's columns, kept in the design's `eta_memo`
#' environment so a descent called many times over one design copies the
#' block once.
#'
#' @details
#' Only a design whose blocks do not move with the coefficients carries that
#' environment; any other gets a fresh [coord_block()] at each call, as
#' before. The entry is keyed on the equation and the columns.
#'
#' @param design The design.
#' @param p The equation, a parameter name.
#' @param X The equation's design, `design[[p]]$X`.
#' @param cols The block's column positions within it.
#'
#' @return What [coord_block()] returns for `X` and `cols`.
#'
#' @seealso [coord_block()], [coord_fit()]
#'
#' @keywords internal
coord_block_at <- function(design, p, X, cols) {
  mm <- attr(design, "eta_memo")
  if (is.null(mm)) return(coord_block(X, cols))
  key <- paste0("block:", p)
  hit <- mm[[key]]
  if (!is.null(hit) && identical(hit$cols, cols)) return(hit$X)
  B <- coord_block(X, cols)
  assign(key, list(cols = cols, X = B), envir = mm)
  B
}


#' Whether a Penalized Block Is Solved with Its Equation's Intercept Profiled
#'
#' @description
#' `TRUE` where the equation carries an unpenalized intercept that is free to
#' move, in which case [coord_fit()] centers the block's columns with the
#' working weights and solves it with the intercept profiled out.
#'
#' @details
#' The block is fitted with the other columns of its equation held, and the
#' intercept is one of them. Where the columns are not centered, each change
#' of a coefficient moves the mean of the fit, which only the intercept can
#' take back, and the intercept is updated in another block. The alternation
#' between the two then converges at a rate set by how close each column is
#' to the constant. Measured on `MASS::UScrime`, whose columns have means up
#' to 33 times their spread, a lasso on the fifteen raw predictors of
#' `log(y)` stopped at a fixed \eqn{\lambda} of 17.4 with the objective
#' \eqn{-\ell + \rho} at -17.84 where its minimum is -21.93, six
#' coefficients against nine, reporting `converged = FALSE`, and the path
#' chose the empty model at \eqn{\lambda = 374.5}. With the columns centered
#' the path chooses nine of them at \eqn{\lambda = 17.4}, the point it
#' reaches with the predictors centered in the data, in 7.0 seconds of
#' processor time against 38.2.
#'
#' Solving with the intercept profiled out is the Frisch-Waugh-Lovell
#' reading of an unpenalized constant: the coefficients of the block are
#' those of the centered columns against the centered working response, and
#' they do not depend on the value the intercept holds. [coord_fit()] then
#' sets the intercept to the value that goes with them, the weighted mean of
#' the working response net of the block, so its step is a joint step in the
#' block and the intercept. A held intercept cannot take that value, so an
#' intercept named in `held_coef` leaves the block uncentered, and so does an
#' equation with no intercept at all.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param p The equation, a parameter name.
#' @param obj The objective, as [statmod_objective()] returns it.
#' @param beta The stacked coefficients.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @seealso [coord_fit()], [coord_colsq()], [parametric_intercept()]
#'
#' @keywords internal
coord_centers <- function(spec, design, p, obj, beta) {
  j <- parametric_intercept(spec, design, p)
  if (is.na(j)) return(FALSE)
  hf <- held_positions(spec, design, obj, beta)
  if (!length(hf$where)) return(TRUE)
  !(obj$split(seq_along(beta))[[p]][j] %in% hf$where)
}


#' Weighted Sums of Squares of Centered Columns
#'
#' @description
#' \eqn{\sum_i w_i (x_{ik} - m_k)^2} for the columns \eqn{k} asked for,
#' which is the curvature a coordinate of a centered block has.
#'
#' @details
#' It is computed as the sum it is, not as
#' \eqn{\sum_i w_i x_{ik}^2 - W m_k^2}, which is the same number in exact
#' arithmetic and loses as many digits as the mean exceeds the spread. On a
#' sparse block the sum runs over the stored entries and adds
#' \eqn{m_k^2} times the weight of the rows the column does not store,
#' which is what a stored zero would contribute, so the block is never
#' densified.
#'
#' A column that is constant over the rows is the intercept itself once
#' centered, and its sum is rounding. Where it falls below \eqn{10^{-10}}
#' of the uncentered sum the uncentered one is returned instead: the kernel
#' leaves such a coordinate where it is, and the value here only has to be
#' a usable step for the table.
#'
#' @param X The block, dense or `dgCMatrix`.
#' @param w The working weights, length `nrow(X)`.
#' @param k The columns, one-based.
#' @param m The weighted means of every column of the block.
#'
#' @return A numeric vector, one entry per column in `k`.
#'
#' @seealso [coord_centers()], [coord_fit()]
#'
#' @keywords internal
coord_colsq <- function(X, w, k, m) {
  Xk <- X[, k, drop = FALSE]
  mk <- m[k]
  if (isS4(Xk)) {
    Xk <- methods::as(Xk, "CsparseMatrix")
    cnt <- diff(Xk@p)
    ci <- rep(seq_along(k), cnt)
    ri <- Xk@i + 1L
    d <- Xk@x - mk[ci]
    s <- numeric(length(k))
    nzw <- numeric(length(k))
    raw <- numeric(length(k))
    if (length(ci)) {
      s[] <- vapply(split(w[ri] * d * d, factor(ci, levels = seq_along(k))),
                    sum, numeric(1))
      nzw[] <- vapply(split(w[ri], factor(ci, levels = seq_along(k))),
                      sum, numeric(1))
      raw[] <- vapply(split(w[ri] * Xk@x^2, factor(ci, levels = seq_along(k))),
                      sum, numeric(1))
    }
    s <- s + mk^2 * (sum(w) - nzw)
  } else {
    D <- sweep(Xk, 2L, mk)
    s <- as.numeric(crossprod(w, D^2))
    raw <- as.numeric(crossprod(w, Xk^2))
  }
  deg <- s <= 1e-10 * raw
  s[deg] <- raw[deg]
  s
}
