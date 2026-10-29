#' @include outer_criteria.R
NULL

#' The Size of a Penalty's Kink
#'
#' @description
#' Returns half the jump of \eqn{\rho'} across the first point at which the
#' penalty is not differentiable, which is the half-width of the
#' subdifferential there.
#'
#' @details
#' A coefficient stays at the kink while the unpenalized score at that point is
#' inside the subdifferential, so a score has to exceed this number for a
#' coefficient to leave zero. It is \eqn{\lambda} for the lasso, SCAD and
#' MCP, and \eqn{\lambda\alpha} for the elastic net, but it is read and never
#' assumed: a penalty built from a density carries the same information in a
#' hyperparameter of its own, and a Laplace prior written by its scale has a
#' kink whose size falls as that scale grows.
#'
#' A kinked penalty is separable over the coefficients, and one under a general
#' map (the generalized-lasso problem) is rejected upstream. Under a diagonal
#' map, which `standardize = TRUE` writes, the kink differs by coordinate. This
#' function reads the first coordinate's kink, and [coord_kinks()] returns the
#' kink of every coordinate.
#'
#' The derivative just past the kink depends on the distance at which it is
#' read (MCP's, for example, is \eqn{\lambda - \epsilon/\gamma}), so the
#' one-sided derivatives at distances \eqn{\epsilon/2} and \eqn{\epsilon} are
#' combined by Richardson extrapolation to the limit. The extrapolation is
#' exact for a penalty that is affine on each side of the kink, as the lasso,
#' SCAD, MCP and the elastic net are.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters, as a list or a named vector.
#' @param eps The larger of the two distances from the kink at which the
#'   derivative is read.
#'
#' @return A single non-negative number, `0` where there is no kink.
#'
#' @seealso [kink_hypers()], [kink_solve()]
#'
#' @keywords internal
kink_scale <- function(pen, theta, eps = 1e-4) {
  th <- as.list(theta)
  k <- penalties7::penalty_kinks(pen, th)
  k <- k[is.finite(k)]
  if (!length(k)) return(0)
  p <- max(1L, as.integer(pen@n_coef))
  at <- function(h) {
    up <- penalties7::penalty_gradient(pen, rep(k[[1L]] + h, p), th)
    dn <- penalties7::penalty_gradient(pen, rep(k[[1L]] - h, p), th)
    as.numeric(up[[1L]] - dn[[1L]]) / 2
  }
  2 * at(eps / 2) - at(eps)
}


#' Which Hyperparameters Set the Size of the Kink
#'
#' @description
#' Returns the names of the hyperparameters whose value changes
#' [kink_scale()].
#'
#' @details
#' Each hyperparameter is moved to another admissible value and the size of the
#' kink is read again, so that a list of families is not consulted. The shape
#' parameters of SCAD and MCP leave the subdifferential at zero unchanged and
#' govern how fast the penalty flattens further out, while \eqn{\lambda} and
#' the elastic net's \eqn{\alpha} both scale it.
#'
#' `unbounded` restricts the answer to the hyperparameters with no upper
#' bound. [path_by_kink()] calls this function with `unbounded = FALSE` and
#' excludes the bounded hyperparameters itself; those are swept over their
#' interval by [path_grid()].
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param unbounded Whether to keep only those with an infinite upper bound.
#'   The default is `TRUE`.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso [kink_scale()]
#'
#' @keywords internal
kink_hypers <- function(pen, theta, unbounded = TRUE) {
  th <- as.list(theta)
  s0 <- kink_scale(pen, th)
  out <- character(0)
  for (h in pen@params) {
    b <- pen@params_bounds[[h]]
    if (unbounded && is.finite(b[2L])) next
    v <- th[[h]]
    th2 <- th
    th2[[h]] <- bounded_bump(v, b)
    s1 <- tryCatch(kink_scale(pen, th2), error = function(e) NA_real_)
    if (is.finite(s1) && abs(s1 - s0) > 1e-10 * max(1, abs(s0))) {
      out <- c(out, h)
    }
  }
  out
}


#' Move a Hyperparameter Without Leaving Its Interval
#'
#' @param v The current value.
#' @param b Its bounds.
#'
#' @return A different admissible value.
#'
#' @keywords internal
bounded_bump <- function(v, b) {
  up <- v * 1.5 + 0.5
  if (up < b[2L]) return(up)
  dn <- (v + b[1L]) / 2
  if (dn > b[1L]) return(dn)
  v
}


#' The Hyperparameter That Gives the Kink a Chosen Size
#'
#' @description
#' Solves `kink_scale(pen, theta) == target` in one named hyperparameter.
#'
#' @details
#' The size of the kink is monotone in such a hyperparameter but not
#' necessarily increasing: a Laplace prior written by its scale has a kink of
#' \eqn{1/\sigma}, which narrows as the hyperparameter grows. The direction of
#' the walk is therefore measured before the root is bracketed, by comparing
#' the size at the current value with the size at twice it. The root is then
#' bracketed by doubling or halving and found with [stats::uniroot()].
#'
#' Where the size is a power of the hyperparameter, as it is for the shipped
#' kinked penalties because the hyperparameter enters a separable penalty as a
#' scale, the inversion is closed and the search is not run. [kink_power()]
#' measures the exponent and [kink_by_power()] checks the answer against the
#' size before it is returned, so that a penalty that does not obey a power
#' law falls back to the search.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force; the value of `name` is where
#'   the search starts.
#' @param name Which one to solve for.
#' @param target The size that the kink should have.
#'
#' @return A single value, or `NA` where the target is out of reach.
#'
#' @seealso [kink_power()], [path_values()]
#'
#' @keywords internal
kink_solve <- function(pen, theta, name, target) {
  th <- as.list(theta)
  b <- pen@params_bounds[[name]]
  size <- function(v) {
    th[[name]] <- v
    kink_scale(pen, th)
  }
  v <- kink_by_power(pen, th, name, target, kink_power(pen, th, name))
  if (!is.na(v)) return(v)
  at <- function(v) size(v) - target
  lo <- max(b[1L], 0) + .Machine$double.eps
  v0 <- as.numeric(th[[name]])[1L]
  s0 <- size(v0)
  if (!is.finite(s0)) return(NA_real_)
  if (abs(s0 - target) <= 1e-12 * max(1, target)) return(v0)
  s1 <- size(min(v0 * 2, b[2L] / 2))
  rising <- is.finite(s1) && s1 > s0
  # walk towards the target, which is up when the size rises with the value
  # and the target is above the current size, and up as well when it falls and
  # the target is below
  up <- xor(target < s0, rising)

  v <- v0
  f <- at(v)
  for (i in seq_len(500L)) {
    nxt <- if (up) v * 2 else v / 2
    if (!is.finite(nxt) || nxt <= lo || nxt >= b[2L]) break
    f_nxt <- at(nxt)
    if (!is.finite(f_nxt)) break
    if (f * f_nxt <= 0) {
      return(tryCatch(stats::uniroot(at, sort(c(v, nxt)), tol = 1e-14)$root,
                      error = function(e) NA_real_))
    }
    v <- nxt
    f <- f_nxt
  }
  NA_real_
}


#' How the Size of the Kink Scales With a Hyperparameter
#'
#' @description
#' Returns the exponent \eqn{k} in \eqn{s(v) = c\,v^k}, read from the size of
#' the kink at the current value and at twice it (or at half it, where twice
#' it is outside the hyperparameter's interval), or `NULL` where there is no
#' such exponent to read.
#'
#' @details
#' A separable penalty is minus a log density and its hyperparameter enters as
#' a scale, so the width of the subdifferential at the kink is a power of it.
#' The exponent is one for the lasso, for SCAD and MCP in \eqn{\lambda} and for
#' the elastic net in both \eqn{\lambda} and \eqn{\alpha}. For the shapes of
#' SCAD and MCP, which do not move the kink, there is no exponent and the
#' result is `NULL`. A Laplace prior written by its scale has exponent
#' \eqn{-1}, and its kink narrows as the hyperparameter grows.
#'
#' The exponent is measured and not assumed, and [kink_by_power()] checks the
#' inversion, so that a penalty whose kink is not a power of its
#' hyperparameter costs two evaluations of the size and then takes the
#' search of [kink_solve()].
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one.
#'
#' @return A list of the exponent `k`, the value `v0` at which it was read
#'   and the size `s0` there; `NULL` where the size does not move or where
#'   either size or the value is not positive and finite.
#'
#' @seealso [kink_solve()], [path_values()]
#'
#' @keywords internal
kink_power <- function(pen, theta, name) {
  th <- as.list(theta)
  b <- pen@params_bounds[[name]]
  v0 <- as.numeric(th[[name]])[1L]
  if (!is.finite(v0) || v0 <= 0) return(NULL)
  s0 <- kink_scale(pen, th)
  if (!is.finite(s0) || s0 <= 0) return(NULL)
  v1 <- if (v0 * 2 < b[2L]) v0 * 2 else v0 / 2
  if (!is.finite(v1) || v1 <= b[1L] || v1 >= b[2L]) return(NULL)
  th[[name]] <- v1
  s1 <- kink_scale(pen, th)
  if (!is.finite(s1) || s1 <= 0) return(NULL)
  k <- log(s1 / s0) / log(v1 / v0)
  if (!is.finite(k) || abs(k) < 1e-8) return(NULL)
  list(k = k, v0 = v0, s0 = s0)
}


#' Invert the Size of the Kink Through a Power Law
#'
#' @description
#' Returns, for each target size, the value that gives it, checked against the
#' size that the value produces.
#'
#' @details
#' The exponent comes from two points, and the relation is checked at the
#' first and the last admissible target before any answer is used. If either
#' check misses by more than a rounding, every target comes back `NA`; a
#' target whose value lands outside the hyperparameter's own interval comes
#' back `NA` on its own. The caller then falls back to bracketing.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one to solve for.
#' @param target The sizes that the kink should have.
#' @param pw What [kink_power()] returned, or `NULL`.
#'
#' @return A numeric vector as long as `target`, `NA` where the
#'   power law gives no value.
#'
#' @seealso [kink_power()], [kink_solve()]
#'
#' @keywords internal
kink_by_power <- function(pen, theta, name, target, pw) {
  out <- rep(NA_real_, length(target))
  if (is.null(pw)) return(out)
  b <- pen@params_bounds[[name]]
  v <- pw$v0 * (target / pw$s0)^(1 / pw$k)
  ok <- is.finite(v) & v > b[1L] & v < b[2L] & target > 0
  if (!any(ok)) return(out)
  # asserted at the extremes of what was asked for: the relation is a power
  # law or it is not, and where it holds at both ends it holds between them
  th <- as.list(theta)
  for (j in unique(c(which(ok)[[1L]], rev(which(ok))[[1L]]))) {
    th[[name]] <- v[[j]]
    s <- tryCatch(kink_scale(pen, th), error = function(e) NA_real_)
    if (!is.finite(s) ||
        abs(s - target[[j]]) > 1e-8 * max(1, abs(target[[j]]))) {
      return(out)
    }
  }
  out[ok] <- v[ok]
  out
}


#' The Largest Score That a Kinked Block Has to Beat
#'
#' @description
#' Returns \eqn{\max_j |\partial(-\ell)/\partial\beta_j|} over the block's own
#' coefficients, evaluated while the block is held at the kink.
#'
#' @details
#' A coefficient leaves the kink when the unpenalized score there exceeds the
#' half-width of the subdifferential, so a kink at least this wide leaves the
#' whole block at zero. That is where a path starts: at the smallest
#' hyperparameter for which the term contributes nothing, which is
#' \pkg{glmnet}'s `lambda.max` written for any separable penalty.
#'
#' The other coefficients are held where the caller left them and are not
#' refitted, so the number is a starting point and not a boundary. The path
#' checks it: a top whose fit is not empty is doubled until it is.
#'
#' Where the kinks of the coordinates differ, as under `standardize`, each
#' score is read against its own coordinate's kink, and the result is
#' expressed in the units of the first coordinate's kink, which is how
#' [kink_solve()] reads it.
#'
#' @param obj The stacked objective.
#' @param beta The current coefficients.
#' @param block One entry of `statmod_blocks()$sparse`.
#' @param hyper The hyperparameters.
#' @param also Other blocks to hold at their kinks while the score is read.
#'   Each block is a list with the stacked positions `index` and the value
#'   `at`, and these blocks are the members of a shared hyperparameter other
#'   than `block`.
#'
#' @return A single number.
#'
#' @keywords internal
path_null_score <- function(obj, beta, block, hyper, also = NULL) {
  th <- as.list(hyper[[block$param]][[block$term]])
  k <- penalties7::penalty_kinks(block$penalty, th)
  k <- k[is.finite(k)]
  at <- if (length(k)) k[[1L]] else 0
  v <- beta
  v[block$index] <- at
  # the other members of a shared axis are held at their own kinks too: the
  # top is the value that empties all of them at once, so none of them may
  # carry a coefficient the score of this one is read against
  for (o in also) v[o$index] <- o$at
  g <- obj$gr(v)[block$index] -
    penalties7::penalty_gradient(block$penalty,
                                 rep(at, length(block$index)), th)
  # EVERY COORDINATE AGAINST ITS OWN KINK. Under a diagonal map, which is what
  # `standardize` writes, coordinate j's kink is the first coordinate's times
  # |d_j|/|d_1|, so the block empties when the largest |g_j|/w_j is at most
  # one. The size returned is the first coordinate's, which is the size
  # kink_scale() and kink_solve() speak. Reading max|g| against the first
  # coordinate's kink alone made the top depend on which covariate was
  # written first: measured on MASS::UScrime, a lasso over fifteen
  # standardized covariates started at 374.5 and the same model written as
  # two blocks sharing an id at 158.4, so the two built different grids and
  # chose different values of lambda.
  w <- tryCatch(coord_kinks(block$penalty, th), error = function(e) NULL)
  if (is.null(w) || length(w) != length(g) || !all(is.finite(w)) ||
      any(w <= 0)) {
    return(max(abs(g)))
  }
  w[[1L]] * max(abs(g) / w)
}


#' The Values That a Path Visits
#'
#' @description
#' Returns a geometric grid of kink sizes from the one that empties the block
#' down to `min_ratio` of it, carried back onto the hyperparameter.
#'
#' @details
#' The grid is geometric in the size of the kink, never in the
#' hyperparameter, so that a penalty whose kink narrows as its hyperparameter
#' grows is swept in the same order as one whose kink widens: from the empty
#' model towards the full one. Values that the penalty cannot reach are
#' dropped.
#'
#' The exponent relating the two is read once by [kink_power()] and every
#' target is inverted through it, in place of bracketing each target on its
#' own. [kink_by_power()] checks the relation before the values are used and
#' returns `NA` where it does not hold, and those targets fall back to
#' [kink_solve()] one at a time.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one the path varies.
#' @param s_max The size of the kink at the top of the path.
#' @param n_values How many points.
#' @param min_ratio The smallest kink size, as a fraction of `s_max`.
#'
#' @return A numeric vector of values for `name`, from the emptiest fit to
#'   the fullest.
#'
#' @seealso [kink_power()], [path_forced()]
#'
#' @keywords internal
path_values <- function(pen, theta, name, s_max, n_values = 40L,
                        min_ratio = 1e-4) {
  s <- exp(seq(log(s_max), log(s_max * min_ratio), length.out = n_values))
  v <- kink_by_power(pen, theta, name, s, kink_power(pen, theta, name))
  todo <- which(is.na(v))
  if (length(todo)) {
    v[todo] <- vapply(s[todo],
                      function(target) kink_solve(pen, theta, name, target),
                      numeric(1))
  }
  v[is.finite(v) & v > 0]
}


#' The Values That a Caller Wrote Out
#'
#' @description
#' Returns the grid that the term named, ordered from the emptiest fit towards
#' the fullest so that the warm starts run in the same direction as in every
#' other path here.
#'
#' @details
#' Which end is the sparse one is a property of the penalty and not of the
#' numbers: the kink of a lasso widens with \eqn{\lambda} and that of a
#' Laplace prior written by its scale narrows with \eqn{\sigma}. The order is
#' therefore taken from the sign of the exponent that [kink_power()] reads,
#' descending where the kink widens with the hyperparameter and ascending
#' otherwise. The function applies no other rule: the value that empties the
#' block does not cap the grid, `min_ratio` does not extend it, and both of
#' those are ways to build a grid.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one the path varies.
#' @param values What the term wrote out.
#'
#' @return The distinct values sorted, from the emptiest fit to the fullest.
#'
#' @seealso [modelterms7::term_values()], [path_values()]
#'
#' @keywords internal
path_forced <- function(pen, theta, name, values) {
  v <- sort(unique(as.numeric(values)))
  pw <- kink_power(pen, theta, name)
  if (!is.null(pw) && pw$k > 0) rev(v) else v
}


#' The Values That a Path Visits Over a Hyperparameter Not Swept by Kink Size
#'
#' @description
#' Returns the grid of a hyperparameter that [path_values()] does not reach:
#' equally spaced points strictly inside its interval where the interval is
#' bounded, and a geometric grid above [shape_floor()] where it has no upper
#' bound.
#'
#' @details
#' [path_values()] walks the size of the kink, from the value that
#' empties the block down, and a bounded hyperparameter cannot reach that end:
#' the elastic net's kink is \eqn{\lambda\alpha}, so at a given \eqn{\lambda}
#' the block stays nonempty for every admissible \eqn{\alpha} and every point
#' of such a path is dropped. A bounded hyperparameter is therefore swept over
#' its interval, and \eqn{\alpha} runs between the ridge and the lasso. The kink
#' still moves through \eqn{\lambda}, which is swept on its own axis.
#'
#' The endpoints are excluded because the bounds are open: the elastic net at
#' \eqn{\alpha = 0} has no kink at all, and the path would be scoring a
#' penalty of another kind.
#'
#' A hyperparameter with no upper bound that does not move the kink is the
#' shape of SCAD or MCP. Its grid is the floor plus
#' `exp(seq(log(0.25), log(25), length.out = n_values))`, which covers the
#' values in common use (3.7 for the SCAD of Fan and Li, 2001, and 3 for the
#' MCP of Zhang, 2010) and does not reach the floor itself, where the operator
#' is set-valued or the penalty is not defined.
#'
#' The floor is the smallest value at which the block can be fitted, which
#' [shape_floor()] derives from the proximal condition at the steps that the
#' block's coordinate descent takes, and not the constant over which the
#' penalty is defined. The two coincide on an ordinary well-conditioned block
#' and differ where the steps are long.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which hyperparameter the path varies.
#' @param n_values How many points.
#' @param steps What [path_steps()] returned, or `NULL`.
#'
#' @return A numeric vector of values for `name`.
#'
#' @seealso [path_values()], [path_bounded()],
#'   [shape_floor()]
#'
#' @keywords internal
path_grid <- function(pen, theta, name, n_values = 25L, steps = NULL) {
  b <- pen@params_bounds[[name]]
  if (is.null(b) || b[2L] <= b[1L]) return(numeric(0))
  n <- max(2L, as.integer(n_values))
  if (is.finite(b[2L])) {
    return(seq(b[1L], b[2L], length.out = n + 2L)[seq_len(n) + 1L])
  }
  # No upper bound and no effect on the kink: this is the SHAPE of SCAD and
  # MCP, which govern how fast the penalty flattens beyond the kink. Neither
  # route above reaches it -- the kink-size path cannot move a value the kink
  # does not depend on, and there is no interval to span -- so the grid is
  # geometric above the floor. It spans the values the literature uses, 3.7
  # for the SCAD of Fan and Li and 3 for the MCP of Zhang, and cannot reach
  # the floor itself, where the operator is set-valued or the penalty is not
  # defined.
  shape_floor(pen, theta, name, steps) +
    exp(seq(log(0.25), log(25), length.out = n))
}


#' The Step That a Coordinate Descent Would Take on a Block
#'
#' @description
#' Returns \eqn{1/v_j}, with \eqn{v_j = \sum_i w_i x_{ij}^2}, one per column,
#' at the coefficients in hand.
#'
#' @details
#' This is the step that the sweeps of [coord_fit()] use, and it determines
#' whether a shape parameter is admissible: SCAD's proximal operator requires
#' the step to be below \eqn{a - 1} and MCP's below \eqn{\gamma}, where the
#' step is multiplied by the squared entries of the diagonal map that
#' standardization applies. The lower limit of a shape is therefore a property
#' of the data, and not only the constant 2 or 1 over which the penalty is
#' defined.
#'
#' The working weights are those of [coord_working()] for the observed
#' information at the coefficients in hand. The result is `NULL` for a
#' structural block, for a block without columns, and wherever the working
#' weights or the squared column sums are not finite and positive; the caller
#' reads `NULL` as the penalty's own bound and nothing more. The function also
#' returns `NULL` when any error occurs inside the computation, since the grid
#' built from the steps is a starting grid and may be approximate.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param block One entry of `statmod_blocks()$sparse`.
#' @param beta The current coefficients.
#' @param split The objective's own splitter, which puts a stacked
#'   coefficient vector back into one piece per distribution parameter.
#'
#' @return A numeric vector, one step per column, or `NULL`.
#'
#' @seealso [shape_floor()], [coord_fit()]
#'
#' @keywords internal
path_steps <- function(spec, design, block, beta, split) {
  if (isTRUE(block$structural) || is.null(block$cols)) return(NULL)
  out <- tryCatch({
    d <- design[[block$param]]
    X <- coord_block(d$X, block$cols)
    coef <- split(beta)
    ep <- statmod_eta(spec, design, coef)
    wq <- coord_working(spec, ep, coef, design, block$param, FALSE,
                          NULL)
    if (is.null(wq)) return(NULL)
    v <- wxsq(X, wq$w, spec@threads)
    if (any(!is.finite(v)) || any(v <= 0)) return(NULL)
    # the step coord_fit() will take: on a block solved with its intercept
    # profiled out that is the centered column's, see coord_centers()
    if (coord_centers(spec, design, block$param, list(split = split), beta)) {
      m <- as.numeric(xtv(X, wq$w, spec@threads)) / sum(wq$w)
      v <- coord_colsq(X, wq$w, seq_len(ncol(X)), m)
    }
    1 / v
  }, error = function(e) NULL)
  out
}


#' The Smallest Admissible Value of a Shape Parameter
#'
#' @description
#' Returns the lower end of the range over which a shape may be swept: the
#' penalty's own bound, raised to wherever its proximal operator starts to
#' exist at the steps that the block's coordinate descent will take.
#'
#' @details
#' The limit is found by bisection on whether
#' [penalties7::penalty_prox_spec()] returns a table for the penalty at the
#' block's steps, so that a family added later is covered, and the
#' \eqn{a - 1} of SCAD and the \eqn{\gamma} of MCP do not appear here. A grid
#' starting just above the constant over which the penalty is defined would
#' otherwise contain points at which the block cannot be fitted with those data
#' by the only route that a kinked penalty has.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one is the shape.
#' @param steps What [path_steps()] returned, or `NULL`.
#'
#' @return A single number: the penalty's lower bound where `steps` is `NULL`
#'   or not finite, or where the search does not find an admissible value
#'   below the upper bound; otherwise the smallest value, to bisection
#'   accuracy, at which the proximal table exists.
#'
#' @seealso [path_steps()], [path_grid()]
#'
#' @keywords internal
shape_floor <- function(pen, theta, name, steps = NULL) {
  b <- pen@params_bounds[[name]]
  lo <- b[1L]
  if (is.null(steps) || !length(steps) || any(!is.finite(steps))) return(lo)
  th <- as.list(theta)
  admits <- function(v) {
    th[[name]] <- v
    !is.null(tryCatch(penalties7::penalty_prox_spec(pen, th, steps),
                      error = function(e) NULL))
  }
  hi <- lo + max(1, max(steps))
  for (i in seq_len(60L)) {
    if (!is.finite(hi) || hi >= b[2L]) return(lo)
    if (admits(hi)) break
    hi <- lo + (hi - lo) * 2
  }
  if (!admits(hi)) return(lo)
  for (i in seq_len(200L)) {
    mid <- (lo + hi) / 2
    if (mid <= lo || mid >= hi) break
    if (admits(mid)) hi <- mid else lo <- mid
  }
  hi
}


#' Default Grid Settings of a Path
#'
#' @description
#' Returns the number of values and the depth used for a kinked hyperparameter
#' whose term named neither.
#'
#' @details
#' The constructors [modelterms7::lasso()], [modelterms7::enet()],
#' [modelterms7::scad()] and [modelterms7::mcp()] carry these settings on
#' their own signatures, for example `lasso(x, n_lambda = 25, min_ratio =
#' 1e-4)` and `enet(x, n_lambda = 25, n_alpha = 5)`, so that these
#' constructors do not use the settings here. The settings here are used by a
#' term that declares a kinked penalty without an argument for the grid, such
#' as [modelterms7::random()] under a Laplace prior, whose hyperparameters are
#' those of the effects' distribution.
#'
#' `kink` is the number of points on the axis swept by the size of the kink,
#' which runs geometrically over a factor of `1/min_ratio` (four decades at
#' the default). `other` is the number of points on an axis that is swept over
#' an interval instead, such as \eqn{\alpha} between the ridge and the lasso
#' or the shape of SCAD or MCP over its usual range; with a product of axes,
#' every point added there multiplies the number of fits.
#'
#' @return A named list of `kink`, `other` and `min_ratio`.
#'
#' @seealso [path_grid()], [statmod_grid_size()],
#'   [modelterms7::term_grid()]
#'
#' @keywords internal
path_fallbacks <- function() {
  list(kink = 25L, other = 5L, min_ratio = 1e-4)
}


#' Is a Hyperparameter Bounded Above?
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param name Which hyperparameter.
#'
#' @return A single logical.
#'
#' @keywords internal
path_bounded <- function(pen, name) {
  b <- pen@params_bounds[[name]]
  !is.null(b) && is.finite(b[2L])
}


#' Is a Hyperparameter Swept by the Size of Its Kink?
#'
#' @description
#' Returns TRUE where the geometric path of path_values() reaches it: the
#' hyperparameter scales the kink and has no upper bound, so the value that
#' empties the block is admissible.
#'
#' @details
#' The two conditions fail in different ways and both have to hold. The
#' elastic net's alpha scales the kink and is bounded by one, so an admissible
#' value of it cannot empty the block at a given lambda; the shape of
#' SCAD and MCP has no upper bound and does not move the kink at all, so the
#' solve has nothing to solve. Either way the sweep is
#' [path_grid()].
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which hyperparameter.
#'
#' @return A single logical.
#'
#' @seealso [path_values()], [path_grid()]
#'
#' @keywords internal
path_by_kink <- function(pen, theta, name) {
  !path_bounded(pen, name) &&
    name %in% kink_hypers(pen, theta, unbounded = FALSE)
}


#' Which Hyperparameters a Path Has to Select
#'
#' @description
#' Returns the rows of an index like [outer_hyper_index()]'s, for the
#' hyperparameters whose penalty has a kink.
#'
#' @details
#' A gradient search is not the instrument for these. The penalized mode is a
#' piecewise smooth function of the hyperparameter, differentiable while the
#' active set holds and turning a corner every time a coefficient joins it or
#' leaves, so a criterion read at that mode inherits the corners and a
#' quasi-Newton step is reading a slope that is about to change. A grid is not
#' affected by the corners, and warm starts make it cheap.
#'
#' Every hyperparameter of a kinked penalty is varied, except those that
#' its term holds and those of an `id` group that one member holds.
#'
#' @param spec A [StatmodSpec()].
#' @param blocks The blocks.
#' @param hyper The hyperparameters.
#'
#' @return A data frame with one row per axis, with columns `parameter`,
#'   `term` and `name`, and a list column `members` that holds, for each axis,
#'   the data frame of the hyperparameters represented by the axis. It has
#'   zero rows if there is no hyperparameter to estimate.
#'
#' @keywords internal
path_rows <- function(spec, blocks, hyper) {
  rows <- list()
  held <- statmod_held(spec)
  # a label ONE member holds is held for the group, which is the same rule
  # outer_hyper_index() applies on the other machine
  held_ids <- names(statmod_held_ids(statmod_penalty_keys(spec)))
  for (b in blocks$sparse) {
    th <- as.list(hyper[[b$param]][[b$term]])
    # EVERY hyperparameter of the penalty, less the ones the term holds.
    # What is estimated is what the caller did not fix, and a shape that no
    # longer scales the kink is swept over a grid of its own rather than
    # left out for want of one.
    # matched LITERALLY: a term's name is its call deparsed, so it carries
    # parentheses and dots, and a regex over it silently matches nothing --
    # which is a held hyperparameter swept anyway, the defect this exists to
    # prevent
    parts <- strsplit(held, "\r", fixed = TRUE)
    mine <- vapply(parts, function(q)
      length(q) == 3L && identical(q[[1L]], b$param) &&
        identical(q[[2L]], b$term), logical(1))
    want <- setdiff(b$penalty@params,
                    vapply(parts[mine], `[[`, character(1), 3L))
    ids <- if (is.null(b$ids)) character(0) else b$ids
    for (h in want) {
      if (h %in% names(ids) && ids[[h]] %in% held_ids) next
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = b$param, term = b$term, name = h,
        id = if (h %in% names(ids)) ids[[h]] else NA_character_,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) {
    out <- data.frame(parameter = character(0), term = character(0),
                      name = character(0), stringsAsFactors = FALSE)
    out$members <- list()
    return(out)
  }
  path_group(do.call(rbind, rows))
}


#' Collapse the Shared Hyperparameters of a Path
#'
#' @description
#' Turns one line per estimated hyperparameter into one line per axis: those
#' carrying the same label become a single axis, swept once, whose value is
#' written into each of them.
#'
#' @details
#' It plays for a path the role that [index_group()] plays for the outer
#' search. The axis is identified by its first member, so the code that reads
#' `parameter`, `term` and `name` (to find the block to which the axis
#' belongs, to find whether it scales the kink, or to read the history)
#' continues to read the same fields, and the `members` list column lists the
#' hyperparameters that the axis represents.
#'
#' If a line shares its label with no other line, it is its own axis and
#' carries itself as its single member.
#'
#' @param rows One line per estimated hyperparameter, with an `id` column.
#'
#' @return One line per axis, with a `members` list column and no `id`.
#'
#' @keywords internal
path_group <- function(rows) {
  key <- ifelse(is.na(rows$id), paste0("\r", seq_len(nrow(rows))), rows$id)
  first <- !duplicated(key)
  grp <- match(key, key[first])
  out <- rows[first, setdiff(names(rows), "id"), drop = FALSE]
  rownames(out) <- NULL
  out$members <- lapply(seq_len(nrow(out)), function(g)
    rows[grp == g, c("parameter", "term", "name"), drop = FALSE])
  out
}


#' Assign the Observations to Folds
#'
#' @description
#' Returns a fold number per observation, either the assignment that the
#' method carries or a permutation drawn with a fixed seed.
#'
#' @details
#' The permutation is drawn after `set.seed(20260811L)`, and the caller's
#' random-number state is restored afterwards, so the partition is the same
#' whatever the caller's seed. Passing `folds` explicitly makes two criteria
#' comparable on a partition of the caller's own. The values of
#' `folds` are relabeled to the integers `1, 2, ...` in the order of their
#' factor levels, and a length other than `n` signals an error.
#'
#' @param n The number of observations.
#' @param k How many folds, at most `n`.
#' @param folds The method's own assignment, or `numeric(0)`.
#'
#' @return An integer vector of length `n`.
#'
#' @keywords internal
cv_folds <- function(n, k, folds) {
  if (length(folds)) {
    if (length(folds) != n) {
      stop(sprintf("'folds' has %d entries but there are %d observations.",
                   length(folds), n), call. = FALSE)
    }
    return(as.integer(as.factor(folds)))
  }
  k <- min(as.integer(k), n)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(20260811L)
  sample(rep_len(seq_len(k), n))
}


#' The Hyperparameters as Plain Named Vectors
#'
#' @param hyper The nested list of hyperparameters, one list per distribution
#'   parameter and one entry per penalized term.
#'
#' @return The same list, in which each entry is converted to a named numeric
#'   vector and the distribution parameters that carry no penalized term are
#'   dropped.
#'
#' @keywords internal
hyper_plain <- function(hyper) {
  out <- lapply(hyper, function(per) {
    lapply(per, function(v) stats::setNames(as.numeric(v), names(v)))
  })
  out[lengths(out) > 0L]
}


#' Set One Hyperparameter
#'
#' @param hyper The hyperparameters.
#' @param row One row of [path_rows()]'s index.
#' @param value The value to write.
#'
#' @return The hyperparameters, in which that one is changed.
#'
#' @keywords internal
hyper_set <- function(hyper, row, value) {
  # ONE value written under EVERY member's key: a shared axis is one number
  # standing for several penalties, which stay several penalties. With nothing
  # shared the axis carries itself alone and this is the single assignment
  # that was here before.
  mem <- row$members
  if (is.null(mem)) {
    hyper[[row$parameter]][[row$term]][[row$name]] <- value
    return(hyper)
  }
  mem <- mem[[1L]]
  for (i in seq_len(nrow(mem))) {
    hyper[[mem$parameter[i]]][[mem$term[i]]][[mem$name[i]]] <- value
  }
  hyper
}


#' Carry a Term's Matrix Input Onto a Subset of the Rows
#'
#' @description
#' Adds, as a column of the subset, each matrix that a term was given, so that
#' rebuilding the model on those rows finds it there.
#'
#' @details
#' `interpret_formula()` evaluates a term's call as
#' `eval(call, data, env)`, so a name is looked up in `data` first
#' and in the formula's environment after. `data.frame(X = X, y = y)`
#' splits a matrix into `X.x1 ... X.xp`, leaving no column `X`, so
#' `lasso(X)` reaches past the data to the matrix in the calling
#' environment. The fit itself is correct (the matrix is captured once, and
#' the coefficients are identical to those of the other spelling), but the
#' fold cannot evaluate the term on its rows: the name still resolves to all
#' the rows.
#'
#' The matrix is already on the built term, and `term_build()` checked at
#' the full fit that it has one row per observation, so the rows of a fold are
#' the same rows by position. Binding the subset here gives the fold the
#' spelling that the documentation of the matrix input recommends, with the
#' matrix as a column of the data.
#'
#' A matrix carries no knots, contrasts or levels, so subsetting it and
#' evaluating it again give the same block. A formula input is left
#' untouched.
#'
#' It applies where `input_expr` is a plain symbol, which is the case for
#' which the name can be bound. A call such as `lasso(scale(X))` keeps only its
#' own value on the term, never the `X` that its re-evaluation would need, so
#' it is left to the error that names it.
#'
#' @param spec A [StatmodSpec()] whose terms are built.
#' @param sub The subset of the data, already taken.
#' @param i The rows used to take it.
#' @param n The number of rows on which the fit was built.
#'
#' @return `sub`, which has a column per matrix input that the terms carry.
#'
#' @seealso [cv_curve()]
#'
#' @keywords internal
cv_bind_inputs <- function(spec, sub, i, n) {
  for (p in spec@distrib@params) {
    for (tm in spec@terms[[p]]) {
      pn <- tryCatch(S7::prop_names(tm), error = function(e) character(0))
      if (!all(c("input", "input_expr") %in% pn)) next
      ex <- tm@input_expr
      if (!is.symbol(ex)) next
      nm <- as.character(ex)
      # already a column of the data: the caller wrote it the other way and
      # the subset carries it
      if (nm %in% names(sub)) next
      x <- tm@input
      if (inherits(x, "formula") || is.null(dim(x)) || nrow(x) != n) next
      ok <- tryCatch({
        sub[[nm]] <- x[i, , drop = FALSE]
        TRUE
      }, error = function(e) FALSE)
      # a storage a data frame will not hold as one column is left to the
      # error, which says what to write; it is not worth densifying here,
      # which is what the sparse input was passed to avoid
      if (!ok) next
    }
  }
  sub
}


#' The Held-Out Criterion of Every Point of a Path
#'
#' @description
#' Refits the model on each training fold along the whole path and scores it on
#' the fold left out. Returns, at each point, the mean over folds of minus twice
#' the weighted log-likelihood of the held-out rows divided by the sum of their
#' weights, and its standard error across folds.
#'
#' @details
#' The path is run fold by fold and not point by point, so that each fit
#' starts from the previous point's coefficients. Each fold reapplies to its
#' rows the terms built on all the rows, as [predict.StatmodFit()] reapplies
#' them to new data: the knots of a basis, the levels of a factor and the
#' scales of a standardized block come from the covariates of every row and
#' never from the response. A basis rebuilt on the training rows alone would
#' not cover a held-out row beyond their range, and the fold would fail.
#'
#' The training rows and the held-out rows carry the weights and offsets of
#' `spec` at those rows, including the offsets given through the `offsets`
#' argument of [statmod()]. Where the fit of a fold at a point does not
#' converge, or its criterion is not finite, the fold's value at that point is
#' missing and is counted in `n_fail`.
#'
#' @param spec A [StatmodSpec()], holding the weights and offsets of the fit.
#' @param data The data with which the fit was called.
#' @param inner_optimizer The inner method.
#' @param hypers A list of hyperparameter settings, one per path point.
#' @param folds A fold number per observation.
#' @param run The combination of the outer axes to which each point belongs.
#'   The warm start begins again at the head of each combination, where the
#'   kink jumps back up. `NULL` treats the whole list as one run.
#'
#' @return A list with `cvm`, the mean over folds at each point; `cvse`, its
#'   standard error across the folds with a finite value; `n_fail`, the number
#'   of folds without a value at each point; and `error`, the message of the
#'   first fold that could not rebuild the model, or `NULL`. `cvm` is missing
#'   at a point where every fold lacks a value.
#'
#' @seealso [cv()]
#'
#' @keywords internal
cv_curve <- function(spec, data, inner_optimizer, hypers,
                     folds, run = NULL) {
  nf <- max(folds)
  m <- length(hypers)
  # One seed per fold, drawn HERE: a fold's fit is then a deterministic
  # function of the fold whether it runs in this process or in a worker,
  # which is what lets the result be independent of `cores` bit for bit
  # even where a family's starting values draw at random. The caller's own
  # stream is put back afterwards, so a fold cannot shift what the session
  # draws next.
  seeds <- sample.int(.Machine$integer.max, nf)
  one_fold <- function(f) {
    if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", globalenv())
      on.exit(assign(".Random.seed", old_seed, globalenv()), add = TRUE)
    }
    set.seed(seeds[f])
    dev_f <- rep(NA_real_, m)
    keep <- folds != f
    # the TEST fold needs them too: term_predict() evaluates a matrix input's
    # expression in the new data with baseenv() as its enclosure, so a name
    # that is not a column there is not found at all
    train <- cv_bind_inputs(spec, data[keep, , drop = FALSE], keep, nrow(data))
    test <- cv_bind_inputs(spec, data[!keep, , drop = FALSE], !keep, nrow(data))
    # the design depends on the fold and not on the path, so it is built once
    # here rather than once per point. Measured, that is worth about 4 per cent
    # and not the most of it, which was the guess: at 200 observations and 20
    # columns a cross-validated path costs 0.88 seconds a fit and almost all of
    # it is the proximal iteration.
    # THE FOLD REAPPLIES THE TERMS BUILT ON ALL THE ROWS (Giovanni,
    # 2026-10-07), as predict() does: rebuilt on the training rows, a basis
    # did not reach a held-out row past their range and statmod_respec() of
    # the test rows stopped the whole fit. Only covariates enter a term's
    # construction, so no response of the held-out rows reaches the training
    # fit. The weights and the offsets are the fit's own, at these rows.
    at_rows <- function(s0, rows) {
      S7::set_props(
        s0,
        weights = if (length(spec@weights) == length(rows)) spec@weights[rows]
                  else rep(1, sum(rows)),
        offsets = lapply(spec@offsets, function(o)
          if (length(o) == length(rows)) o[rows] else o))
    }
    ts <- tryCatch(at_rows(statmod_respec(spec, train), keep),
                   error = function(e) conditionMessage(e))
    # statmod_spec() builds a FRESH specification, so the thread count does
    # not travel with it as it does through statmod_respec(), which starts
    # from the one it is given: until 2026-08-21 a fold fell back to the
    # class default of 1 and a cross-validated path was single-threaded
    # whatever statmod(threads =) asked for. Measured on a lasso over a
    # gamma response, bic() gained 1.85x from eight threads where cv()
    # gained 1.06x. The two levels still do not nest -- where the folds
    # themselves run in separate PROCESSES each one fits on one thread, which
    # is what local_threads() documents -- so the count is passed
    # on only when this process is running them.
    if (!is.character(ts)) {
      ts@threads <- if (spec@cores > 1L) 1L else spec@threads
    }
    if (is.character(ts)) {
      # a fold that cannot be rebuilt is a configuration error and not a
      # numerical one -- a term whose input is not in `data` cannot be
      # re-evaluated on a subset of it -- so the message is carried out
      # rather than dropped: every fold fails for the same reason, and the
      # path would otherwise score NA everywhere and keep its starting value
      return(list(dev = dev_f, err = ts))
    }
    # the terms are reapplied when the design is built, so an input the fold
    # cannot evaluate fails here, and is carried out as above
    td <- tryCatch(statmod_design(ts), error = function(e) conditionMessage(e))
    if (is.character(td)) return(list(dev = dev_f, err = td))
    cfgs <- inner_settings(inner_optimizer, ts@distrib)
    obj <- statmod_objective(ts, hypers[[1L]], td, cfgs$expected, cfgs$emethod)
    warm <- statmod_start(ts, td, obj, NULL)
    # a fold scales its SCAD and MCP penalties as statmod() scales the whole
    # data's, at its own starting coefficients: left unscaled, the path would
    # score a different penalty from the one the fit it selects for
    tc <- statmod_curv(ts, td, warm, cfgs$expected, cfgs$emethod)
    if (!is.null(tc)) ts <- tc
    tb <- statmod_blocks(ts, td)
    # the held-out rows carry their own weights and offsets as well, so that
    # the deviance they score is the fit's weighted log-likelihood
    hs <- at_rows(statmod_respec(ts, test), !keep)
    hd <- statmod_design(hs)
    tbj <- tb
    warm0 <- warm
    for (j in seq_len(m)) {
      if (!is.null(run) && j > 1L && run[[j]] != run[[j - 1L]]) {
        warm <- warm0
        tbj <- tb
      }
      r <- tryCatch(statmod_alternate(ts, td, tbj, hypers[[j]], inner_optimizer,
                                      warm, cfgs$expected, cfgs$emethod,
                                      cfgs$maxit, cfgs$tol, verbosity(0),
                                      hold_refresh = TRUE),
                    error = function(e) NULL)
      tbj <- blocks_at_kink(tb, hypers[[j]])
      if (is.null(r) || !isTRUE(r$converged)) next
      warm <- r$par
      cf <- r$obj$split(r$par)
      d <- tryCatch(-2 * statmod_loglik_at(hs, cf, hd),
                    error = function(e) NA_real_)
      if (is.finite(d) && sum(hs@weights) > 0) dev_f[j] <- d / sum(hs@weights)
    }
    list(dev = dev_f, err = NULL)
  }
  rows <- worker_map(spec@cores, nf, one_fold)
  dev <- do.call(rbind, lapply(rows, `[[`, "dev"))
  errs <- unlist(lapply(rows, `[[`, "err"))
  err <- if (length(errs)) errs[[1L]] else NULL
  ok <- colSums(is.finite(dev))
  cvm <- apply(dev, 2L, mean, na.rm = TRUE)
  cvse <- apply(dev, 2L, stats::sd, na.rm = TRUE) / sqrt(pmax(ok, 1))
  cvm[ok == 0L] <- NA_real_
  list(cvm = cvm, cvse = cvse, n_fail = nf - ok, error = err)
}


#' Run Independent Units, in This Process or Over Workers
#'
#' @description
#' Applies a body to each of `n` independent units (a cross-validation fold,
#' a combination of a path's product grid, a bootstrap replicate) over
#' `cores` separate R processes, and in this process when `cores` is 1.
#' Results come back in unit order whatever the number of processes, which
#' makes the answer independent of the count: the units share nothing, so the
#' same bodies run either way.
#'
#' @details
#' The units are independent by construction (a fold is a complete refit
#' on its own rows, a path combination restarts its warm chain from the
#' sweep's own starting coefficients), so they go by processes, with the
#' safeguards that `optimizers7::multistart` records: under `pkgload`
#' the run stays sequential, because a worker loads the installed copy and
#' S7 objects built in the development namespace do not dispatch correctly
#' against it; a cluster that cannot start, or workers that cannot load
#' the package, fall back to sequential with a warning rather than fail
#' the fit. A fit inside a worker takes a fresh specification and is
#' therefore sequential by construction: the two levels of parallelism do
#' not nest.
#'
#' @param cores The number of separate R processes, a single whole number
#'   (see [check_cores()]). It is distinct from the `threads` of a fit, which
#'   counts threads inside one process. At most `n` processes are opened.
#' @param n How many units.
#' @param body The unit's body, a function of the unit index.
#' @param what The unit's name, for the warnings.
#'
#' @return A list of the bodies' results, in unit order.
#'
#' @seealso [cv_curve()], [statmod_path()]
#'
#' @keywords internal
worker_map <- function(cores, n, body, what = "folds") {
  workers <- min(check_cores(cores), n)
  sequential <- function() lapply(seq_len(n), body)
  if (workers <= 1L) return(sequential())
  if (isNamespaceLoaded("statmodels7") &&
      exists(".__DEVTOOLS__", asNamespace("statmodels7"), inherits = FALSE)) {
    return(sequential())
  }
  cl <- try(parallel::makePSOCKcluster(workers), silent = TRUE)
  if (inherits(cl, "try-error")) {
    warning("Could not start ", workers, " worker processes; running the ",
            what, " sequentially.", call. = FALSE)
    return(sequential())
  }
  on.exit(parallel::stopCluster(cl), add = TRUE)
  have <- try(parallel::clusterEvalQ(
    cl, requireNamespace("statmodels7", quietly = TRUE)), silent = TRUE)
  if (inherits(have, "try-error") || !all(unlist(have))) {
    warning("The worker processes could not load statmodels7, so the ", what,
            " were run\n  sequentially. They are separate R sessions and an ",
            "uninstalled or\n  differently-located copy is invisible to ",
            "them.", call. = FALSE)
    return(sequential())
  }
  parallel::parLapply(cl, seq_len(n), body)
}


#' Check the Number of Processes
#'
#' @description
#' Returns `cores` as an integer if it is a single whole number of at least 1,
#' and signals an error otherwise. `cores` is the number of separate R
#' processes over which a function distributes independent fits; it is not the
#' `threads` of a fit, which counts threads that share memory inside one
#' process.
#'
#' @param cores The value of a `cores` argument.
#'
#' @return A single integer, at least 1.
#'
#' @examples
#' statmodels7:::check_cores(2)
#' try(statmodels7:::check_cores(0))
#'
#' @keywords internal
check_cores <- function(cores) {
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != round(cores)) {
    stop("'cores' must be a single whole number, at least 1.", call. = FALSE)
  }
  as.integer(cores)
}


#' Choose a Point of the Path
#'
#' @description
#' Returns the index of the smallest criterion, or of the largest kink whose
#' criterion is within one standard error of it.
#'
#' @details
#' The one-standard-error rule takes the sparsest fit that is not measurably
#' worse than the best one, which is Breiman's rule as \pkg{glmnet} applies it.
#' The path runs from the emptiest fit to the fullest, so the largest kink
#' among the admissible points is the first of them. The rule falls back to
#' the smallest criterion where `se` is `NULL` or not finite there.
#'
#' @param value The criterion at each point.
#' @param se Its standard error, or `NULL`.
#' @param rule `"min"` or `"1se"`.
#'
#' @return A single index, or `NA` if every point was unusable.
#'
#' @references
#' Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
#' *Classification and Regression Trees*. Wadsworth.
#'
#' @keywords internal
path_pick <- function(value, se = NULL, rule = "min") {
  if (!any(is.finite(value))) return(NA_integer_)
  j <- which.min(value)
  if (!identical(rule, "1se") || is.null(se) || !is.finite(se[j])) return(j)
  ok <- which(is.finite(value) & value <= value[j] + se[j])
  if (!length(ok)) return(j)
  min(ok)
}


#' Select the Hyperparameters of a Kinked Penalty Along a Path
#'
#' @description
#' Sweeps the hyperparameters of the kinked penalties over grids of kink sizes
#' and keeps the setting at which the criterion is best, by its minimum or by
#' the one-standard-error rule.
#'
#' @details
#' The grid runs from the kink that empties the block down to `min_ratio`
#' of it, so the sweep goes from the sparsest fit towards the fullest and every
#' fit starts from the previous one's coefficients. The top of the grid is
#' computed at the coefficients in hand and not at a refitted null; where the
#' fit at the top does not empty the block, the top is doubled until it does.
#'
#' For a term that carries several of these hyperparameters, every combination
#' is visited where the term sets `search = "grid"` and one coordinate at a
#' time where it sets `"cyclic"`. Between terms the alternation is cyclic either
#' way, so that the cost is the product within a term and the sum across terms.
#' Each axis is built at the settings of the axes outside it, which makes the
#' elastic net's grid a family of \eqn{\lambda} axes and not a single one, and
#' the axis swept by kink size is put innermost so that the warm starts walk
#' along it. The runs of a product grid are independent and are distributed over
#' worker processes where the specification requests them.
#'
#' A pass that would visit the points scored by the previous one is not run.
#' The top of the path is computed again at every pass, because the rest of
#' the model moves, and where it has not moved the grid is the one already in
#' hand. The function makes at most `sweeps` passes, and only where more than
#' one hyperparameter is selected.
#'
#' Where the model also carries hyperparameters that are twice differentiable,
#' those are estimated by [outer_fit()] inside each point of the path, by
#' `nested_method`, so that the two kinds are not mixed into one search. The
#' unpenalized coefficients of a dispersion are then held at the joint mode
#' along the path when `nested_method` is a marginal criterion, and are
#' estimated on that criterion once, at the setting chosen.
#'
#' An error is signaled where `method` is [aic()], [bic()] or [cv()] and the
#' path scores no point. A warning is signaled where the rule is `"min"`,
#' the axis is swept by kink size, and the minimum lies at the dense end of the
#' path and the criterion is still falling there.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param blocks The blocks.
#' @param hyper The hyperparameters.
#' @param inner_optimizer The inner method.
#' @param method The [OuterMethod()] that scores the points of the path.
#' @param optimizer The optimizer for the differentiable hyperparameters.
#' @param beta The starting coefficients.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#' @param maxit,tol The inner budget.
#' @param vb The verbosity.
#' @param data The data with which the fit was called, which cross-validation
#'   needs in order to refit on part of it; the weights and offsets are those
#'   held by `spec`.
#' @param rows The axes to select, as [path_rows()] returns them.
#' @param sweeps The largest number of cyclic passes, `2` by default.
#' @param nested_method The [OuterMethod()] that estimates the differentiable
#'   hyperparameters inside each point of the path, or `NULL`, in which case
#'   `method` is used.
#'
#' @return A list with `par`, `hyper`, `value`, `criterion`, `converged`,
#'   `obj`, `hist_blocks`, `hist_inner`, `hist_outer`, `iterations`,
#'   `evaluations`, `exact_gradient` and `exact_hessian`: the fields of
#'   [outer_fit()]'s result that a path produces. `criterion` is the value of
#'   `method` at the setting chosen, and `hist_outer` has one row for each
#'   point visited.
#'
#' @seealso [cv()], [path_rows()]
#'
#' @keywords internal
statmod_path <- function(spec, design, blocks, hyper, inner_optimizer, method,
                         optimizer, beta, emethod, maxit, tol, vb, data,
                         rows, sweeps = 2L,
                         nested_method = NULL) {
  expected <- info_expected(method@hessian)
  # the information statmod() reads the curvature of a scaled SCAD or MCP
  # with at the fitted coefficients, which a point is scored with as well
  curv_expected <- inner_settings(inner_optimizer)$expected
  smooth_idx <- outer_hyper_index(spec, blocks)
  is_cv <- identical(method@kind, "cv")
  # the smooth hyperparameters are estimated inside each point of the path by
  # THEIR OWN criterion, which is not the path's: a smoothing parameter is
  # read at the mode by reml() while a lasso's lambda is swept by bic()
  inner_crit <- if (is.null(nested_method)) method else nested_method
  # A DISPERSION IS READ AT THE JOINT MODE ALONG THE PATH, and on the
  # criterion once, at the value the path chose. A marginal criterion also
  # estimates the unpenalized coefficients of a dispersion or a shape, and
  # searching them at every point failed where the active set jumps: on
  # UScrime two SCAD points of 25 did not converge and were scored NA, and
  # under scad(a = 30) the point skipped was the best, BIC 30.40 against the
  # 37.25 chosen. The smooth hyperparameters, where there are any, are
  # still estimated inside each point.
  path_crit <- inner_crit
  if (!is.null(path_crit) && !outer_minimize(path_crit)) {
    path_crit@marginal <- "none"
  }
  gam_final <- !is.null(inner_crit) && !outer_minimize(inner_crit) &&
    !identical(inner_crit@kind, "cv") &&
    length(marginal_coords(spec, design, inner_crit)$where) > 0L
  nested <- nrow(smooth_idx) > 0L && !is.null(inner_crit) &&
    !identical(inner_crit@kind, "cv")
  obj0 <- statmod_objective(spec, hyper, design, expected, emethod)

  # one fit at one setting, warm-started, with the differentiable
  # hyperparameters estimated inside it where there are any
  # `sp` is the specification a fit runs under, and it is an argument
  # because one caller needs a different one: a run dispatched to a worker
  # PROCESS must fit on a single thread, the two levels of
  # threads and cores not nesting, while every other call here runs
  # in this process and uses the count the caller asked for.
  fit_at <- function(hy, warm, bk = blocks, sp = spec, final = FALSE) {
    if (nested || (final && gam_final)) {
      crit <- if (final) inner_crit else path_crit
      r <- tryCatch(outer_fit(sp, design, bk, hy, inner_optimizer,
                              crit, optimizer, warm, emethod, maxit, tol,
                              vb_inner(vb)), error = function(e) NULL)
      if (!is.null(r)) return(r)
    }
    r <- statmod_alternate(sp, design, bk, hy, inner_optimizer, warm,
                           expected, emethod, maxit, tol, vb_inner(vb),
                           hold_refresh = TRUE)
    r$hyper <- hy
    r$criterion <- NA_real_
    r
  }

  # The criterion is read at the hyperparameters the fit ACTUALLY used. Where
  # the smooth ones were estimated inside this point, `r$hyper` carries what
  # they were estimated to and `hy` still carries what the path came in with:
  # scoring at `hy` then reads a penalty at one smoothing parameter against
  # coefficients fitted at another, and the penalized information that
  # mismatch produces need not be positive definite. Measured, it was not:
  # every point of the path scored NA, path_pick() had nothing to choose
  # between, and the hyperparameter kept its starting value while the sweep
  # looked like it had run.
  score_at <- function(r, hy) {
    cf <- r$obj$split(r$par)
    at <- if (!is.null(r$hyper)) r$hyper else hy
    # the specification the point was fitted under: a scaled SCAD or MCP
    # carries the curvature the alternation settled at there, and the degrees
    # of freedom have to read that one and not the curvature at the start
    sp <- if (!is.null(r$spec)) r$spec else spec
    # and read once more at the point's own coefficients, which is what the
    # fit returned at this value does before it reports anything: the
    # curvature is a function of the coefficients, so the two then score the
    # same penalty. Read at the curvature of the alternation's last pass, an
    # MCP point on MASS::UScrime scored a trace of 9.0 where the fit at the
    # same value reported 10.96
    sc <- tryCatch(statmod_curv(sp, design, cf, curv_expected, emethod),
                   error = function(e) NULL)
    if (!is.null(sc)) sp <- sc
    act <- statmod_active(sp, blocks, r$par, at)
    m <- statmod_pe(sp, design, cf, at, method, emethod, act)
    if (is.null(m)) NA_real_ else m$value
  }

  # the top of the path: the kink that leaves every coefficient of the block
  # at zero
  # THE TOP OF THE PATH IS NOT A CONSTANT. It is the size of the kink that
  # empties the block, read from the score at the coefficients in hand, so
  # it moves with the rest of the model: a smooth beside this term has its
  # own smoothing parameter estimated INSIDE each point of the path, and
  # every such fit changes the score this is read from. Reading it once at
  # the start left the path anchored to the state the search began in.
  path_top <- function(hy, at) {
    tp <- stats::setNames(numeric(nrow(rows)), rownames(rows))
    for (i in seq_len(nrow(rows))) {
      row <- rows[i, ]
      b <- path_block(blocks, row)
      ob <- statmod_objective(spec, hy, design, expected, emethod)
      tp[[i]] <- path_top_shared(ob, at, blocks, row, hy, b)
      # a hyperparameter that does not scale the kink has no top of this
      # kind: no value of it empties the block
      if (!path_by_kink(b$penalty, hy[[row$parameter]][[row$term]],
                        row$name)) next
      # path_null_score() reads the score with the OTHER coefficients where
      # the fit left them rather than at their optimum with this block
      # zeroed, so it is a starting point and not a boundary -- and
      # measured, it can be an order of magnitude short: on eight
      # coefficients of which three carried real signal it returned 26.5
      # where the block first empties near 500, so the path covered a nearly
      # flat stretch of the criterion and every fit chose its sparse end.
      # Doubling until the fit really is empty is the check the
      # documentation already described, and it costs a handful of fits.
      for (step in seq_len(24L)) {
        v <- kink_solve(b$penalty, hy[[row$parameter]][[row$term]],
                        row$name, tp[[i]])
        if (!is.finite(v) || v <= 0) break
        r <- tryCatch(fit_at(hyper_set(hy, row, v), at),
                      error = function(e) NULL)
        if (is.null(r)) break
        # every MEMBER's block must be empty, not only the first's: one value
        # multiplies each penalty as it is, so the value that empties one may
        # leave another full and the sparse end of the path would not be
        # sparse there
        if (max(abs(r$par[path_member_index(blocks, row)])) <= 1e-8) break
        tp[[i]] <- tp[[i]] * 2
      }
    }
    tp
  }
  top <- path_top(hyper, beta)

  cur <- hyper
  hist <- list()
  best <- list(value = Inf, hyper = cur, fit = NULL)

  # A PRODUCT WITHIN A TERM AND AN ALTERNATION BETWEEN TERMS. What a term's
  # own hyperparameters do is not separable -- the elastic net's kink is
  # lambda*alpha, so the criterion has a valley along the hyperbolas
  # lambda*alpha = constant and a sweep moving one coordinate at a time
  # descends it in steps of its own grid -- while two terms in one formula
  # are two blocks, and taking their product would make two elastic nets
  # 10^4 points where alternating makes them 100 + 100.
  keys <- paste(rows$parameter, rows$term, sep = "\r")
  groups <- unname(split(seq_len(nrow(rows)),
                         factor(keys, levels = unique(keys))))

  # Whether a hyperparameter SCALES the kink, which is what decides where its
  # axis sits in the product: that axis is the one running from the emptiest
  # fit towards the fullest, so it goes innermost, where the warm starts walk
  # along it and where the history reports it. Asked without building the
  # grid, the ordering asking it of every axis.
  #
  # It is asked of the hyperparameter and NOT of how its values were arrived
  # at. A written-out lambda is still the axis the path descends -- what the
  # caller fixed is which values, not what they mean -- and reading the two
  # as one question put alpha in the history where lambda belonged.
  scales_kink <- function(i, th) {
    row <- rows[i, ]
    path_by_kink(path_block(blocks, row)$penalty,
                 th[[row$parameter]][[row$term]], row$name)
  }

  # THE VALUES ONE AXIS TAKES AT A GIVEN SETTING OF THE OTHERS, which is what
  # makes the grid a product of one axis by a FAMILY of axes rather than of
  # two axes fixed in advance. A hyperparameter that scales the kink descends
  # from the value emptying the block AT THOSE OTHER SETTINGS: the elastic
  # net's lambda_max is kink/alpha and moves by a factor of nine across
  # alpha, while a SCAD's is one number whatever its shape. Neither is
  # written down; both follow from carrying the size of the kink back onto
  # the hyperparameter where it stands.
  #
  # The steps a shape's floor is read at come from `beta`, the coefficients
  # the sweep starts from, and not from the point of the path in hand: the
  # floor says which shapes the block can be FITTED at, and a floor that
  # moved from point to point would make the axis a different axis at each
  # of them.
  axis <- function(i, th) {
    row <- rows[i, ]
    b <- path_block(blocks, row)
    thb <- th[[row$parameter]][[row$term]]
    # the caller wrote the values out, so there is nothing to build: the
    # emptying value does not cap them and min_ratio does not extend them
    forced <- statmod_values(spec, row)
    if (!is.null(forced)) {
      return(list(vals = path_forced(b$penalty, thb, row$name, forced),
                  kink = FALSE))
    }
    # HOW FINELY is the TERM's answer, and it is on the term's own
    # signature: a block of four columns and one of four hundred want
    # different grids, and a criterion is put to every hyperparameter of the
    # model, the smooth ones included, so it cannot know which it is looking
    # at. `path_fallbacks()` is reached only by a term that declares a
    # kinked penalty without offering an argument for the grid.
    fb <- path_fallbacks()
    if (!path_by_kink(b$penalty, thb, row$name)) {
      nv <- statmod_grid_size(spec, row, fb$other)
      return(list(vals = path_grid(b$penalty, thb, row$name, nv,
                                   path_steps(spec, design, b, beta,
                                              obj0$split)),
                  kink = FALSE))
    }
    nv <- statmod_grid_size(spec, row, fb$kink)
    mr <- statmod_min_ratio(spec, row, fb$min_ratio)
    list(vals = path_values(b$penalty, thb, row$name, top[[i]], nv, mr),
         kink = TRUE)
  }

  # the product over a set of axes, built one axis at a time so that each is
  # computed at the settings of those already fixed
  expand <- function(idx, th) {
    if (!length(idx)) return(list(th))
    v <- axis(idx[[1L]], th)$vals
    if (!length(v)) return(expand(idx[-1L], th))
    out <- list()
    for (val in v) {
      out <- c(out, expand(idx[-1L], hyper_set(th, rows[idx[[1L]], ], val)))
    }
    out
  }

  # The settings one group visits. The axis swept by KINK SIZE goes last, so
  # it is the one rebuilt at every combination of the others and the one the
  # warm starts walk along, from the emptiest fit towards the fullest. `run`
  # says which combination each point belongs to, which is what lets the
  # warm start restart at the head of a row and the end-of-path warning ask
  # about that row rather than about the whole product.
  plan_of <- function(idx) {
    ord <- idx[order(vapply(idx, scales_kink, logical(1), th = cur))]
    inner <- ord[[length(ord)]]
    combos <- expand(ord[-length(ord)], cur)
    hys <- list()
    run <- integer(0)
    for (ci in seq_along(combos)) {
      v <- axis(inner, combos[[ci]])$vals
      for (val in v) {
        hys[[length(hys) + 1L]] <- hyper_set(combos[[ci]], rows[inner, ], val)
        run <- c(run, ci)
      }
    }
    # the ends of a WRITTEN-OUT grid are the caller's, so there is nothing to
    # widen and nothing to warn about
    list(idx = ord, inner = inner, hys = hys, run = run,
         kink = scales_kink(inner, cur) &&
           is.null(statmod_values(spec, rows[inner, ])))
  }

  seen <- vector("list", length(groups))
  for (s in seq_len(if (nrow(rows) > 1L) sweeps else 1L)) {
    moved <- FALSE
    # ... and again here, at the coefficients and hyperparameters this sweep
    # starts from, which the previous sweep moved
    if (s > 1L) top <- path_top(cur, beta)
    for (gi in seq_along(groups)) {
      idx <- groups[[gi]]
      # ASKED OF THE TERM, and a group is one term. It is not the criterion's
      # to answer: the same criterion is put to the smooth hyperparameters of
      # the model, which are read at the mode and not swept at all, so most
      # of what it is asked about could not use the answer. Per term is also
      # what keeps one term's choice off another's.
      product <- identical(statmod_search(spec, rows[idx[[1L]], ]), "grid")
      # cyclic is one plan per axis, each picked on its own; the product is
      # one plan over all of them, picked jointly
      specs <- if (product) list(idx) else as.list(idx)
      sig <- vector("list", length(specs))

      for (pn in seq_along(specs)) {
        # BUILT HERE and not before the loop, so that a cyclic sweep reads
        # the value the previous axis was just set to: an axis's grid is a
        # function of the others, which is the whole reason the two searches
        # differ.
        pl <- plan_of(specs[[pn]])
        hys <- pl$hys
        if (!length(hys)) next
        inner <- pl$inner
        row <- rows[inner, ]
        sig[[pn]] <- vapply(hys, function(h)
          as.numeric(h[[row$parameter]][[row$term]][[row$name]]), numeric(1))
        # A PASS THAT WOULD VISIT THE SAME POINTS AGAIN IS NOT RUN. The top
        # of the path is refreshed between passes because the rest of the
        # model moves, and where it has not moved the grid is the one just
        # scored: re-fitting it is a whole product of fits for an answer
        # already in hand.
        if (!is.null(seen[[gi]]) &&
            isTRUE(all.equal(seen[[gi]][[pn]], sig[[pn]]))) next

        cv_err <- NULL
        if (is_cv) {
          folds <- cv_folds(spec@n_obs, method@nfolds, method@folds)
          cc <- cv_curve(spec, data, inner_optimizer, hys,
                         folds, pl$run)
          value <- cc$cvm
          se <- cc$cvse
          cv_err <- cc$error
        } else {
          # THE COMBINATIONS ARE INDEPENDENT BY CONSTRUCTION: each run of
          # the product restarts its warm chain from `beta` with the full
          # blocks, and `cur` moves only after every run is scored, so the
          # runs share nothing and go through worker_map() -- the same
          # bodies in the same order whatever the count. What is NOT
          # parallelized is the points WITHIN a run: measured
          # (piano_parallel.txt, voce 8), a point paid cold costs 2.2-3.2x
          # the warm chain, so splitting a chain either slows the
          # single-process default by that factor or makes the result
          # depend on the worker count.
          #
          # Within a run, the point just fitted is what the next one
          # screens against: the grid runs from the emptiest fit towards
          # the fullest, so the kink shrinks and the strong rule has a
          # previous size to compare with. At the head of a combination
          # the kink jumps back up, so the screening starts again from the
          # whole block.
          runs <- split(seq_along(hys), factor(pl$run, levels = unique(pl$run)))
          # where the runs go to worker processes each one fits on a single
          # thread; where worker_map() keeps them here they use the count
          # the caller asked for
          run_spec <- spec
          if (spec@cores > 1L) run_spec@threads <- 1L
          one_run <- function(ri) {
            idxs <- runs[[ri]]
            warm <- beta
            bk <- blocks
            out <- rep(NA_real_, length(idxs))
            for (k in seq_along(idxs)) {
              j <- idxs[[k]]
              r <- fit_at(hys[[j]], warm, bk, sp = run_spec)
              bk <- blocks_at_kink(blocks, hys[[j]])
              # a nested search whose stopping rule did not fire is scored
              # where its decrement says it stands at the optimum, see
              # outer_fit()
              if (!isTRUE(r$converged) && !isTRUE(r$settled)) next
              warm <- r$par
              out[[k]] <- score_at(r, r$hyper)
            }
            out
          }
          got <- worker_map(spec@cores, length(runs), one_run,
                            what = "path combinations")
          value <- rep(NA_real_, length(hys))
          for (ri in seq_along(runs)) value[runs[[ri]]] <- got[[ri]]
          se <- NULL
        }

        at <- function(j, i) as.numeric(
          hys[[j]][[rows$parameter[[i]]]][[rows$term[[i]]]][[rows$name[[i]]]])
        vals <- vapply(seq_along(hys), at, numeric(1), i = inner)
        # the whole combination, for a reader: `name` and `value` carry the
        # axis the path descends and this carries the rest, so one row of the
        # history is one point whatever the search visited
        other <- setdiff(pl$idx, inner)
        setting <- if (!length(other)) rep("", length(hys)) else
          vapply(seq_along(hys), function(j) paste(vapply(other, function(i)
            sprintf("%s=%s", rows$name[[i]], format(signif(at(j, i), 4))),
            character(1)), collapse = ", "), character(1))

        hist[[length(hist) + 1L]] <- data.frame(
          sweep = s, parameter = row$parameter, term = row$term,
          name = row$name, value = vals, setting = setting, criterion = value,
          se = if (is.null(se)) NA_real_ else se, stringsAsFactors = FALSE)

        # NOT A POINT OF THE PATH WAS SCORED. Skipping quietly leaves the
        # hyperparameter at whatever it came in with -- its default, where the
        # caller asked for it to be chosen -- and the fit then reports success
        # at a value nothing selected, which is the failure this refuses to
        # produce. Where the criterion said why, it says so here.
        # ... and only where the criterion is one that selects. A marginal
        # criterion reaching a kinked row scores nothing BY CONSTRUCTION --
        # the Laplace expansion it is has no second derivative at the kink --
        # and the hyperparameter keeping the value it was given is the
        # documented answer there, not a failure.
        if (method@kind %in% c("aic", "bic", "cv") && !any(is.finite(value))) {
          stop(sprintf(paste0(
            "%s() could not score a single point of the path for '%s' in",
            " '%s',\n  so %s was not chosen and would have kept the value",
            " with which it entered.%s"),
            method@kind, row$term, row$parameter, row$name,
            if (!is.null(cv_err))
              paste0("\n  Every fold failed to rebuild the model: ", cv_err,
                     "\n  A term's input must be a column of 'data', or it",
                     " cannot be re-evaluated on\n  a subset of it.")
            else paste0("\n  No fit along it converged, or the criterion was",
                        " not defined at any of them.")),
            call. = FALSE)
        }
        j <- path_pick(value, se, method@rule)
        if (is.na(j)) next
        # A choice at either end MAY be the grid's rather than the criterion's,
        # and the two are told apart by asking whether the criterion is still
        # improving there. It is not always: the top of the path now empties
        # the block by construction, so the criterion is FLAT across the
        # emptied stretch and index 1 is a legitimate minimum. Warning on the
        # index alone said "the criterion was still falling there" of a
        # criterion that was doing nothing of the kind -- a message naming a
        # cause that is not the real one, which section 7 records as worse than
        # no message.
        # The question is asked WITHIN the chosen point's own combination: in
        # a product each combination carries a path of its own, and the ends
        # that could be widened are that path's.
        mine <- which(pl$run == pl$run[[j]])
        k <- match(j, mine)
        falling <- function() {
          o <- if (k == 1L) mine[[2L]] else mine[[length(mine) - 1L]]
          v <- value[c(j, o)]
          if (!all(is.finite(v))) return(FALSE)
          v[1L] < v[2L] - 1e-8 * max(1, abs(v[2L]))
        }
        # ... and only for a hyperparameter swept by kink size. A bounded one
        # is swept over the whole of its own interval, so an answer at either
        # end is the interval's and there is nothing to widen. Of the two ends
        # of a kink-size path only the dense one can be widened: the sparse
        # end empties the block, so an answer there is the empty block and
        # min_ratio, which moves the dense end, does not reach it.
        if (pl$kink && identical(method@rule, "min") &&
            length(mine) > 1L && k == length(mine) && falling()) {
          warning(sprintf(paste0("The path for '%s' in '%s' stopped at its ",
                                 "dense end (%s = %s).\n  The criterion was ",
                                 "still falling there, so widen the path ",
                                 "with min_ratio\n  or set the value ",
                                 "yourself."),
                          row$term, row$parameter, row$name,
                          format(signif(vals[[j]], 4))), call. = FALSE)
        }
        for (i in pl$idx) {
          r2 <- rows[i, ]
          v <- hys[[j]][[r2$parameter]][[r2$term]][[r2$name]]
          if (!identical(cur[[r2$parameter]][[r2$term]][[r2$name]], v)) {
            moved <- TRUE
          }
          cur <- hyper_set(cur, r2, v)
        }
        best$value <- value[[j]]
      }
      seen[[gi]] <- sig
    }
    if (!moved) break
  }

  final <- fit_at(cur, beta, final = TRUE)
  cf <- final$obj$split(final$par)
  act <- statmod_active(spec, blocks, final$par, final$hyper)
  crit <- if (is_cv) best$value else {
    m <- statmod_pe(spec, design, cf, final$hyper, method, emethod, act)
    if (is.null(m)) NA_real_ else m$value
  }
  list(par = final$par, hyper = final$hyper, value = final$value,
       criterion = crit,
       converged = isTRUE(final$converged) || isTRUE(final$settled),
       obj = final$obj, hist_blocks = final$hist_blocks,
       hist_inner = final$hist_inner,
       hist_outer = do.call(rbind, hist),
       iterations = length(hist), evaluations = nrow(do.call(rbind, hist)),
       exact_gradient = FALSE, exact_hessian = FALSE)
}


#' The Block to Which a Path Row Belongs
#'
#' @param blocks The blocks.
#' @param row One row of [path_rows()]'s index.
#'
#' @return One entry of `blocks$sparse`.
#'
#' @keywords internal
path_block <- function(blocks, row) {
  for (b in blocks$sparse) {
    if (identical(b$param, row$parameter) && identical(b$term, row$term)) {
      return(b)
    }
  }
  stop(sprintf("There is no penalized block for '%s' in '%s'.", row$term,
               row$parameter), call. = FALSE)
}


#' The Coefficients That an Axis Covers, Across Its Members
#'
#' @description
#' Returns the union of the blocks that the axis's members penalize, as
#' positions in the stacked coefficient vector.
#'
#' @details
#' A shared axis is one value multiplying several penalties, so whether the
#' block is empty has to be decided over all of them at once. If the axis has a
#' single member, the result is that block's own index.
#'
#' @param blocks The block split.
#' @param row One row of [path_rows()]'s index.
#'
#' @return An integer vector of positions.
#'
#' @keywords internal
path_member_index <- function(blocks, row) {
  mem <- row$members
  if (is.null(mem)) return(path_block(blocks, row)$index)
  mem <- mem[[1L]]
  sort(unique(unlist(lapply(seq_len(nrow(mem)), function(i)
    path_block(blocks, mem[i, , drop = FALSE])$index), use.names = FALSE)))
}


#' The Top of a Path Over Every Member of an Axis
#'
#' @description
#' Returns the size of the kink, in the axis's own block's units, at which
#' every member of the axis has all its coefficients at the kink.
#'
#' @details
#' A shared axis is one value multiplying several penalties, which may sit in
#' different equations (a lasso on the mean and a lasso on the scale with
#' one `id`). Each member's top is read from its own score, while every member
#' is held at its kink, and carried onto the value of the hyperparameter through
#' its own penalty; the axis starts at the largest of those values, written
#' back as a size of the first member's kink, which is the size in which the
#' rest of the path is expressed. Two blocks of one equation sharing an `id`
#' therefore start where a single block over both starts, because the score of
#' a column does not depend on how the columns were grouped. If the axis has a
#' single member, this is [path_null_score()] of that block.
#'
#' @param obj The stacked objective.
#' @param beta The coefficients at which the scores are read.
#' @param blocks The block split.
#' @param row One row of [path_rows()]'s index.
#' @param hyper The hyperparameters.
#' @param b The axis's own block, `path_block(blocks, row)`.
#'
#' @return A single number, a size of `b`'s kink.
#'
#' @seealso [path_null_score()], [path_member_index()]
#'
#' @keywords internal
path_top_shared <- function(obj, beta, blocks, row, hyper, b) {
  mem <- row$members
  if (is.null(mem) || nrow(mem[[1L]]) < 2L) {
    return(path_null_score(obj, beta, b, hyper))
  }
  mem <- mem[[1L]]
  mb <- lapply(seq_len(nrow(mem)), function(i)
    path_block(blocks, mem[i, , drop = FALSE]))
  kink_at <- function(bl) {
    k <- penalties7::penalty_kinks(bl$penalty,
                                   as.list(hyper[[bl$param]][[bl$term]]))
    k <- k[is.finite(k)]
    if (length(k)) k[[1L]] else 0
  }
  held <- lapply(mb, function(bl) list(index = bl$index, at = kink_at(bl)))
  vals <- vapply(seq_along(mb), function(i) {
    bl <- mb[[i]]
    th <- as.list(hyper[[bl$param]][[bl$term]])
    s <- path_null_score(obj, beta, bl, hyper, also = held[-i])
    tryCatch(kink_solve(bl$penalty, th, mem$name[[i]], s),
             error = function(e) NA_real_)
  }, numeric(1))
  vals <- vals[is.finite(vals) & vals > 0]
  if (!length(vals)) return(path_null_score(obj, beta, b, hyper))
  th <- as.list(hyper[[b$param]][[b$term]])
  th[[row$name]] <- max(vals)
  kink_scale(b$penalty, th)
}


#' Choose the Hyperparameters by Cross-Validation
#'
#' @description
#' Builds the object with which [statmod()] chooses the hyperparameters of
#' the penalties with a kink by cross-validation: the points of a path of
#' values are scored by the log-likelihood that the fit assigns to
#' observations on which it was not fitted, and the best is kept. The result
#' is passed as `sparse_criterion`. It is not an outer criterion:
#' `outer_criterion = cv()` is an error, because the smoothing parameters are
#' estimated by [reml()], [ml()], [aic()] or [bic()].
#'
#' @details
#' # Cross-validation for a penalty with a kink
#'
#' A marginal criterion, [reml()] or [ml()], approximates an integral by a
#' Laplace expansion at the penalized mode, which requires the second
#' derivative of the penalty there. The mode of a lasso, a SCAD or an MCP
#' sits at the kink for every coefficient that it sets to zero, and the
#' derivative does not exist at the kink, so the criterion is undefined at the
#' point at which it would be read.
#'
#' Cross-validation scores prediction and not a posterior, and it needs
#' nothing from the penalty beyond a fit.
#'
#' Passed as `sparse_criterion` to a model with no hyperparameter of a kinked
#' penalty left to estimate, `cv()` has no effect.
#'
#' # The criterion
#'
#' The criterion is the mean over folds of \eqn{-2\ell_f/w_f}, where
#' \eqn{\ell_f} is the weighted log-likelihood of the fold left out at the
#' coefficients fitted on the other folds, and \eqn{w_f} is the sum of the
#' weights of its rows. Each
#' training fit reapplies to its rows the terms built on all the rows (knots,
#' factor levels and scales of standardization), as [predict.StatmodFit()]
#' does for new data, and carries the weights and offsets of the fit.
#'
#' `rule = "1se"` takes the largest kink whose criterion is within one
#' standard error of the smallest, which is the sparsest fit that is not
#' measurably worse.
#'
#' # The folds
#'
#' Unless `folds` is given, the observations are assigned to folds by a
#' permutation drawn with a fixed seed, and the caller's random-number state
#' is restored afterwards, so that the partition does not depend on the
#' caller's seed.
#'
#' # The path
#'
#' The penalized mode is a piecewise smooth function of the hyperparameter:
#' differentiable while the active set holds, turning a corner whenever a
#' coefficient joins it or leaves. A criterion read there inherits the
#' corners, so a gradient search would read a slope about to change.
#'
#' The hyperparameter is therefore swept over a grid, geometric in the size
#' of the kink, from the value that leaves the block empty down to
#' `min_ratio` of it. Every fit starts from the previous one's coefficients.
#'
#' # Which hyperparameters
#'
#' The path covers every hyperparameter of the term's penalty that the
#' constructor did not hold at a number.
#'
#' # A product within a term, an alternation between terms
#'
#' For a term with several kinked hyperparameters, every combination is
#' visited under `search = "grid"` and one hyperparameter is swept at a time
#' under `"cyclic"`.
#'
#' The choice belongs to the term, through `enet(X, search =)` and
#' [modelterms7::term_search()], and not to this criterion. The criterion is
#' also applied to the model's smooth hyperparameters, which are read at the
#' mode and not swept, so that an argument of this kind would have no meaning
#' for them.
#'
#' Between two terms the sweep alternates whichever search each one names, so
#' that `y ~ lasso(X) + enet(R)` costs the two blocks added and not
#' multiplied, and a term that sets a product does not make the other term
#' pay for it.
#'
#' A term that names neither is swept over the product, because the cyclic sweep
#' traverses a cross through the point in hand and can stop where each
#' coordinate is separately best without being jointly so. Its cost is the
#' product of the term's grids, whereas the cost of the cyclic sweep is their
#' sum. At two hyperparameters this is `n_lambda * n_alpha` fits against
#' `n_lambda + n_alpha` per pass. If three or more hyperparameters are
#' estimated, the cost of the product grows exponentially, and `"cyclic"` is
#' there for that case.
#'
#' # How long the path is, and how far down it reaches
#'
#' The length of the path and its depth are settings of the term, written on
#' its own signature: `lasso(x, n_lambda = 25, min_ratio = 1e-4)` and
#' `enet(x, n_lambda = 25, n_alpha = 5)`.
#'
#' The two numbers differ because the axes do. \eqn{\lambda} descends the
#' size of the kink over four decades, while \eqn{\alpha} spans one bounded
#' interval, so that the default product is 25 by 5.
#'
#' # The shape of the grid
#'
#' The grid is not a rectangle for the elastic net and is a rectangle for
#' SCAD and MCP. The elastic net's kink is \eqn{\lambda\alpha}, so the value
#' emptying the block is \eqn{\lambda_{\max} = \kappa/\alpha} and every
#' \eqn{\alpha} carries its own \eqn{\lambda} axis, descending from its own
#' top. The shapes of SCAD and MCP do not change the kink, so that
#' \eqn{\lambda_{\max}} is one number whatever the shape and the two axes form
#' a rectangle. Each axis is built at the settings of the axes outside it,
#' which covers both cases without a rule about families.
#'
#' # The cost
#'
#' Each point of the path costs `nfolds` fits. The number of points is the
#' term's `n_lambda` for one kinked hyperparameter, and the product of its
#' axes for several. The fits along a path are warm-started within each fold.
#'
#' The folds run over separate R processes when
#' `statmod(cores = 4)` requests them. One seed is
#' drawn per fold in the parent, so that the answer is identical at any
#' number of processes.
#'
#' @param nfolds How many folds, a single number of at least 2. The default
#'   is `10`. It is ignored when `folds` is given.
#' @param folds A fold number per observation, an integer vector as long as
#'   the data, for a user-defined partition or to compare two criteria on the
#'   same one. `NULL`, the default, draws a partition with a fixed seed.
#' @param rule `"min"` for the best criterion, `"1se"` for the sparsest fit
#'   within one standard error of it. It is matched with [match.arg()].
#'
#' @return An [OuterMethod()] object of kind `"cv"`, carrying `nfolds`,
#'   `rule` and `folds`. Its `hessian` is `"observed"` and its `k` is
#'   `NA_real_`; this criterion reads neither.
#'
#' @references
#' Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
#' *Classification and Regression Trees*. Wadsworth.
#'
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. *Journal of
#' Statistical Software* 33(1), 1--22.
#'
#' @seealso [aic()] and [bic()] for the other prediction-error criteria,
#'   [reml()] for the marginal ones, [statmod()] for where this is passed,
#'   [hyper()] for reading back what it chose.
#'
#' @examples
#' # Two of eight columns carry signal; the rest are noise for the path to
#' # shrink away.
#' set.seed(1)
#' X <- matrix(rnorm(120 * 8), 120, 8)
#' dd <- data.frame(y = 1 + X %*% c(2, -1.5, rep(0, 6)) + rnorm(120, sd = 0.5))
#' dd$x <- X
#'
#' fit <- statmod(y ~ lasso(x, n_lambda = 12),
#'                distributions7::gaussian1_distrib(), dd,
#'                sparse_criterion = cv(nfolds = 5))
#' hyper(fit)
#'
#' # The two columns that carry signal are recovered.
#' round(coef(fit)$mu[2:3], 3)
#'
#' # "1se" takes the largest kink within one standard error of the best, so
#' # it never keeps more coefficients than "min". On this data the two agree,
#' # because the only sparser point of the path has a criterion more than one
#' # standard error above the minimum.
#' one_se <- statmod(y ~ lasso(x, n_lambda = 12),
#'                   distributions7::gaussian1_distrib(), dd,
#'                   sparse_criterion = cv(nfolds = 5, rule = "1se"))
#' c(min = hyper(fit)$estimate, one_se = hyper(one_se)$estimate)
#' c(min = sum(coef(fit)$mu[-1] != 0), one_se = sum(coef(one_se)$mu[-1] != 0))
#'
#' @export
cv <- function(nfolds = 10, folds = NULL, rule = c("min", "1se")) {
  OuterMethod(kind = "cv", hessian = observed_info(), k = NA_real_,
              nfolds = as.numeric(nfolds), rule = match.arg(rule),
              folds = if (is.null(folds)) numeric(0) else as.numeric(folds))
}


#' The Error for a Cross-Validation With Nothing to Select
#'
#' @description
#' Signals the error that [statmod()] gives when [cv()] is passed as
#' `outer_criterion` to a model in which there is no penalty with a kink that
#' has a hyperparameter left to estimate.
#'
#' @return Does not return.
#'
#' @keywords internal
stop_cv_nothing <- function() {
  stop(paste0("cv() has nothing to select: there is no penalty with a",
              "\n  kink here that has a hyperparameter left to estimate. A",
              " smoothing parameter that is twice\n  differentiable is estimated by reml(),",
              " ml(), aic() or bic(), which read\n  a criterion rather than",
              " refitting on folds."), call. = FALSE)
}


#' Estimate the Hyperparameters by the Route That Each One Admits
#'
#' @description
#' Routes the twice differentiable hyperparameters to [outer_fit()]
#' and the rest to [statmod_path()].
#'
#' @details
#' The split is the same one that determines how the coefficients are fitted,
#' so for a term whose penalty has a kink, both the coefficients and the
#' hyperparameters are handled by methods that do not require a curvature that
#' the penalty does not have.
#'
#' @inheritParams statmod_path
#' @param method The [OuterMethod()] that estimates the twice differentiable
#'   hyperparameters, or `NULL`.
#' @param sparse_method The [OuterMethod()] that scores the path, or `NULL`,
#'   in which case `method` is used.
#'
#' @return The list that [outer_fit()] or [statmod_path()] returns, with the
#'   element `criterion_by`, the [OuterMethod()] whose value is `criterion`;
#'   the list that [statmod_alternate()] returns when the fit runs without a
#'   criterion.
#'
#' @seealso [statmod()]
#'
#' @keywords internal
statmod_select <- function(spec, design, blocks, hyper, inner_optimizer, method,
                           optimizer, beta, emethod, maxit, tol, vb, data,
                           sparse_method = NULL) {
  # Two families of penalty, two criteria. A kinked one is swept over a path
  # of its own values by `sparse_method`; a twice differentiable one is read
  # by `method` at the mode. Where both are present the path is outside and
  # the marginal criterion is estimated inside each of its points, so a model
  # can have its smoothing parameter by REML and its lasso by BIC at once --
  # which one argument could not express.
  sm <- if (is.null(sparse_method)) method else sparse_method
  # a criterion on the expected information reads the approximation of its
  # own object, in the fits it makes as well as in the determinant
  emethod <- select_emethod(method, sparse_method, emethod)
  rows <- path_rows(spec, blocks, hyper)
  if (!nrow(rows)) {
    if (is.null(method)) {
      # held here too: statmod() runs the full alternation, phase on, at
      # whatever the selection ends with, so every select route refines the
      # break-points exactly once
      return(statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                               beta, info_expected(sm@hessian),
                               emethod, maxit, tol, vb, hold_refresh = TRUE))
    }
    if (identical(method@kind, "cv")) stop_cv_nothing()
    res <- outer_fit(spec, design, blocks, hyper, inner_optimizer, method,
                     optimizer, beta, emethod, maxit, tol, vb)
    res$criterion_by <- method
    return(res)
  }
  # where a path runs, the criterion it returns is the path's own
  res <- statmod_path(spec, design, blocks, hyper, inner_optimizer, sm,
                      optimizer, beta, emethod, maxit, tol, vb, data, rows,
                      nested_method = method)
  res$criterion_by <- sm
  res
}
