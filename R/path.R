#' @include outer_criteria.R
NULL

#' The Size of a Penalty's Kink
#'
#' @description
#' Half the jump of \eqn{\rho'} across the first point the penalty is not
#' differentiable at, which is the half-width of the subdifferential there.
#'
#' @details
#' A coefficient stays at the kink while the unpenalized score at that point is
#' inside the subdifferential, so this number is what a score has to exceed for
#' a coefficient to leave zero. It is \eqn{\lambda} for the lasso, SCAD and
#' MCP, and \eqn{\lambda\alpha} for the elastic net, but it is read rather than
#' assumed: a penalty built from a density carries the same information in a
#' hyperparameter of its own, and a Laplace prior written by its scale has a
#' kink whose size falls as that scale grows.
#'
#' The penalty is separable over the coefficients wherever it has a kink, a
#' kinked penalty under a general map being the generalized-lasso problem and
#' rejected upstream, so one coordinate answers for all of them.
#'
#' The derivative just past the kink carries the distance it was read at --
#' MCP's is \eqn{\lambda - \epsilon/\gamma} -- so the two readings are
#' extrapolated to the limit. Without that the shape parameters appear to move
#' the kink by a millionth of themselves, which is enough to be selected for a
#' path and is a measurement of \eqn{\epsilon} rather than of the penalty. The
#' extrapolation is exact for a penalty that is affine on each side of the
#' kink, which the lasso, SCAD, MCP and the elastic net all are.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters, as a list or a named vector.
#' @param eps How far either side of the kink to read the derivative.
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
#' The names whose value moves [kink_scale()], which are the ones a
#' path over the penalty has to vary.
#'
#' @details
#' The question is put to the penalty rather than answered from a list of
#' families: the shape parameters of SCAD and MCP leave the subdifferential at
#' zero unchanged and govern how fast the penalty flattens further out, while
#' \eqn{\lambda} and the elastic net's \eqn{\alpha} both scale it.
#'
#' `unbounded` restricts the answer to the hyperparameters with no upper
#' bound, which is the default choice of what to select. The reference
#' implementations do the same by convention -- \pkg{glmnet} holds
#' \eqn{\alpha} fixed and \pkg{ncvreg} holds \eqn{\gamma} -- and a bounded
#' shape is better swept by hand over a few values than searched.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param unbounded Whether to keep only those with an infinite upper bound.
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
#' \eqn{1/\sigma}, which narrows as the hyperparameter grows. Which way to walk
#' is therefore measured before walking, by comparing the size at the current
#' value with the size at twice it, and only then is the root bracketed by
#' doubling and found with [stats::uniroot()]. A version that assumed
#' the size increases returned `NA` for every target on the Laplace,
#' having walked away from the answer.
#'
#' Where the size is a POWER of the hyperparameter -- which every kinked
#' penalty here turns out to be, the hyperparameter entering a separable
#' penalty as a scale -- the inversion is closed and the search is not run at
#' all. [kink_power()] measures the exponent and the answer is
#' checked against the size before it is returned, so a penalty that does not
#' obey a power law falls back to the search rather than being assumed into
#' one.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The other hyperparameters.
#' @param name Which one to solve for.
#' @param target The size the kink should have.
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
#' The exponent \eqn{k} in \eqn{s(v) = c\,v^k}, read from the size of the kink
#' at the current value and at twice it, or `NA` where there is no such
#' exponent to read.
#'
#' @details
#' A separable penalty is minus a log density and its hyperparameter enters as
#' a scale, so the width of the subdifferential at the kink is a power of it.
#' Measured over four decades, the exponent is exactly one for the lasso, for
#' SCAD and MCP in \eqn{\lambda} and for the elastic net in both \eqn{\lambda}
#' and \eqn{\alpha}, and exactly zero for the shapes of SCAD and MCP, which do
#' not move the kink at all; the largest spread across a decade sweep is
#' \eqn{5.6\times 10^{-16}}. A Laplace prior written by its scale has exponent
#' \eqn{-1}, its kink narrowing as the hyperparameter grows.
#'
#' The exponent is MEASURED rather than assumed, and whoever inverts it checks
#' the answer, so a penalty whose kink is not a power of its hyperparameter
#' costs two evaluations and then takes the search.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one.
#'
#' @return A list of the exponent `k`, the value `v0` it was read at
#'   and the size `s0` there; `NULL` where the size does not move.
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
#' The value giving each target size, checked against the size it produces.
#'
#' @details
#' The check is what licenses the closed route: the exponent came from two
#' points and the relation is asserted at a third before the answer is used.
#' A target that misses by more than a rounding, or that lands outside the
#' hyperparameter's own interval, comes back `NA` and the caller falls
#' back to bracketing.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one to solve for.
#' @param target The sizes the kink should have.
#' @param pw What [kink_power()] returned, or `NULL`.
#'
#' @return A numeric vector as long as `target`, `NA` where the
#'   power law did not answer.
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


#' The Largest Score a Kinked Block Has to Beat
#'
#' @description
#' \eqn{\max_j |\partial(-\ell)/\partial\beta_j|} over the block's own
#' coefficients, with the block held at the kink.
#'
#' @details
#' A coefficient leaves the kink when the unpenalized score there exceeds the
#' half-width of the subdifferential, so a kink at least this wide leaves the
#' whole block at zero. That is where a path starts: at the smallest
#' hyperparameter for which the term contributes nothing, which is
#' \pkg{glmnet}'s `lambda.max` written for any separable penalty.
#'
#' The other coefficients are held where the caller left them rather than
#' refitted, so the number is a starting point and not a boundary. The path
#' checks it: a top whose fit is not empty is doubled until it is.
#'
#' @param obj The stacked objective.
#' @param beta The current coefficients.
#' @param block One entry of `statmod_blocks()$sparse`.
#' @param hyper The hyperparameters.
#'
#' @return A single number.
#'
#' @keywords internal
path_null_score <- function(obj, beta, block, hyper) {
  th <- as.list(hyper[[block$param]][[block$term]])
  k <- penalties7::penalty_kinks(block$penalty, th)
  k <- k[is.finite(k)]
  at <- if (length(k)) k[[1L]] else 0
  v <- beta
  v[block$index] <- at
  g <- obj$gr(v)[block$index] -
    penalties7::penalty_gradient(block$penalty,
                                 rep(at, length(block$index)), th)
  max(abs(g))
}


#' The Values a Path Visits
#'
#' @description
#' A geometric grid of kink sizes from the one that empties the block down to
#' `min_ratio` of it, carried back onto the hyperparameter.
#'
#' @details
#' The grid is geometric in the size of the kink rather than in the
#' hyperparameter, so that a penalty whose kink narrows as its hyperparameter
#' grows is swept in the same order as one whose kink widens: from the empty
#' model towards the full one. Values the penalty cannot reach are dropped.
#'
#' The exponent relating the two is read ONCE and every target inverted
#' through it, rather than each target bracketed on its own. Measured, a
#' bracketing solve costs 4.18 ms against a fit's 62.5 ms, so a path of
#' twenty-five points spent 6.7 per cent of itself locating the values it
#' would visit; through the exponent the whole grid costs four evaluations of
#' the size. [kink_by_power()] checks the relation before the values
#' are used and returns `NA` where it does not hold, and those fall back
#' to [kink_solve()] one at a time.
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


#' The Values a Caller Wrote Out
#'
#' @description
#' The grid the term named, ordered from the emptiest fit towards the fullest
#' so that the warm starts run the way every other path here runs.
#'
#' @details
#' Which end is the sparse one is a property of the penalty and not of the
#' numbers: the kink of a lasso widens with \eqn{\lambda} and that of a
#' Laplace prior written by its scale narrows with \eqn{\sigma}, so the order
#' is settled by asking the penalty which way its kink moves rather than by
#' sorting downwards. Nothing else is applied -- the value that empties the
#' block does not cap the grid and `min_ratio` does not extend it, both
#' of those being ways to build one.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one the path varies.
#' @param values What the term wrote out.
#'
#' @return A numeric vector, from the emptiest fit to the fullest.
#'
#' @seealso [modelterms7::term_values()], [path_values()]
#'
#' @keywords internal
path_forced <- function(pen, theta, name, values) {
  v <- sort(unique(as.numeric(values)))
  pw <- kink_power(pen, theta, name)
  if (!is.null(pw) && pw$k > 0) rev(v) else v
}


#' The Values a Path Visits Over a Bounded Hyperparameter
#'
#' @description
#' An equally spaced grid strictly inside the hyperparameter's own interval.
#'
#' @details
#' [path_values()] walks the SIZE OF THE KINK, from the value that
#' empties the block down, and a bounded hyperparameter cannot reach that end:
#' the elastic net's kink is \eqn{\lambda\alpha}, so at a given \eqn{\lambda}
#' no admissible \eqn{\alpha} empties the block and every point of such a path
#' is dropped. The interval itself is what a bounded shape is swept over
#' instead -- \eqn{\alpha} between the ridge and the lasso, SCAD's and MCP's
#' shape between their limits -- and the sweeps being cyclic, one coordinate
#' at a time, the kink still moves through \eqn{\lambda}.
#'
#' The endpoints are excluded because the bounds are open: the elastic net at
#' \eqn{\alpha = 0} has no kink at all, and the path would be scoring a
#' penalty of another kind.
#'
#' A shape parameter is swept above the smallest value at which the block can
#' be FITTED, which [shape_floor()] derives from the proximal
#' condition at the steps the block's coordinate descent will take, rather
#' than above the constant the penalty is defined over. The two coincide on an
#' ordinary well-conditioned block and differ where the steps are long: with a
#' standardized penalty on a column of spread 20 the floor is 3 where SCAD's
#' own bound is 2, so a quarter of the old grid named shapes no fit could
#' reach.
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


#' The Step a Coordinate Descent Would Take on a Block
#'
#' @description
#' \eqn{1/v_j} with \eqn{v_j = \sum_i w_i x_{ij}^2}, one per column, at the
#' coefficients in hand.
#'
#' @details
#' This is the step the sweeps of [coord_fit()] use, and it is what
#' decides whether a shape parameter is admissible: SCAD's proximal operator
#' needs \eqn{t < a - 1} and MCP's \eqn{t < \gamma}, tightened to
#' \eqn{t d^2} under the diagonal map standardization applies. So the useful
#' lower limit of a shape is a property of the DATA, not the constant 2 or 1
#' the penalty is defined above, and it binds at ordinary settings: measured,
#' a standardized penalty on a column of spread 20 at \eqn{n = 200} needs
#' \eqn{a > 3}, and a Poisson block whose fitted means are near
#' \eqn{10^{-3}} needs \eqn{a > 11}.
#'
#' Everything is guarded, and `NULL` -- which the caller reads as "the
#' penalty's own bound and nothing more" -- is the answer wherever the
#' working weights are not usable. A starting grid may be approximate; what
#' it may not do is fail.
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
    wq <- coord_working(spec, ep, coef, design, block$param, FALSE, "bartlett")
    if (is.null(wq)) return(NULL)
    v <- wxsq(X, wq$w, spec@threads)
    if (any(!is.finite(v)) || any(v <= 0)) return(NULL)
    1 / v
  }, error = function(e) NULL)
  out
}


#' The Smallest Admissible Value of a Shape Parameter
#'
#' @description
#' The lower end of the range a shape may be swept over: the penalty's own
#' bound, raised to wherever its proximal operator starts to exist at the
#' steps the block's coordinate descent will take.
#'
#' @details
#' The limit is DERIVED from the condition and not written down. The question
#' is put to the penalty -- does it produce a table at this step? -- and
#' bisected, so a family added later is covered and neither the \eqn{a - 1} of
#' SCAD nor the \eqn{\gamma} of MCP appears here. A grid starting just above
#' the constant the penalty is defined over would otherwise contain points at
#' which THAT block, with THOSE data, cannot be fitted by the only route a
#' kinked penalty has.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one is the shape.
#' @param steps What [path_steps()] returned, or `NULL`.
#'
#' @return A single number strictly above the penalty's lower bound.
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


#' What a Path Does Where the Term Says Nothing
#'
#' @description
#' The number of values and the depth used for a kinked hyperparameter whose
#' term named neither.
#'
#' @details
#' The five penalized constructors carry these on their own signatures, where
#' a reader can see them -- `lasso(x, n_lambda = 25, min_ratio = 1e-4)`,
#' `enet(x, n_lambda = 25, n_alpha = 5)` -- so nothing here is reached
#' for them. What reaches it is a term that declares a kinked penalty without
#' offering an argument for the grid: [modelterms7::random()] under
#' a Laplace prior is the case, its hyperparameters being whatever the
#' effects' distribution happens to carry.
#'
#' `kink` is the length of the path over the SIZE OF THE KINK, which
#' runs geometrically over `1/min_ratio` -- four decades -- and wants
#' that many points to be smooth in. `other` serves an axis that spans
#' one bounded interval instead, \eqn{\alpha} between the ridge and the lasso
#' or a shape over its useful range, and needs fewer; with a product every
#' extra point there multiplies the fits.
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
#' TRUE where the geometric path of path_values() reaches it: the
#' hyperparameter scales the kink and has no upper bound, so the value that
#' empties the block is admissible.
#'
#' @details
#' The two conditions fail in different ways and both have to hold. The
#' elastic net's alpha scales the kink and is bounded by one, so no
#' admissible value of it empties the block at a given lambda; the shape of
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
#' The rows of an index like [outer_hyper_index()]'s, for the
#' hyperparameters whose penalty has a kink.
#'
#' @details
#' A gradient search is not the instrument for these. The penalized mode is a
#' piecewise smooth function of the hyperparameter, differentiable while the
#' active set holds and turning a corner every time a coefficient joins it or
#' leaves, so a criterion read at that mode inherits the corners and a
#' quasi-Newton step is reading a slope that is about to change. A grid does
#' not care, and warm starts make it cheap.
#'
#' Which hyperparameters are varied is read from the penalty by
#' [kink_hypers()] unless the method names them.
#'
#' @param spec A [StatmodSpec()].
#' @param blocks The blocks.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#'
#' @return A data frame of `parameter`, `term` and `name`.
#'
#' @keywords internal
path_rows <- function(spec, blocks, hyper, method) {
  rows <- list()
  held <- statmod_held(spec)
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
    for (h in want) {
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = b$param, term = b$term, name = h,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) {
    return(data.frame(parameter = character(0), term = character(0),
                      name = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}


#' Assign the Observations to Folds
#'
#' @description
#' A fold number per observation, either the one the method carries or a fresh
#' permutation.
#'
#' @details
#' The permutation is drawn from the caller's stream and put back, so a fit is
#' not silently reproducible only when the caller happens to have seeded.
#' Passing `folds` explicitly is what makes two criteria comparable on the
#' same partition.
#'
#' @param n The number of observations.
#' @param k How many folds.
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


#' The Hyperparameters in the Shape statmod() Accepts
#'
#' @param hyper The nested list the machinery carries.
#'
#' @return The same, with the parameters that carry no penalized term dropped.
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
#' @return The hyperparameters, with that one changed.
#'
#' @keywords internal
hyper_set <- function(hyper, row, value) {
  hyper[[row$parameter]][[row$term]][[row$name]] <- value
  hyper
}


#' Carry a Term's Matrix Input Onto a Subset of the Rows
#'
#' @description
#' Adds each matrix a term was given as a column of the subset, so that
#' rebuilding the model on those rows finds it there.
#'
#' @details
#' `interpret_formula()` evaluates a term's call as
#' `eval(call, data, env)`, so a name is looked up in `data` first
#' and in the formula's environment after. `data.frame(X = X, y = y)`
#' SPLITS a matrix into `X.x1 ... X.xp`, leaving no column `X`, so
#' `lasso(X)` reaches past the data to the matrix in the calling
#' environment. The fit is right -- the matrix is captured once and the
#' coefficients are identical to the other spelling -- but the fold cannot
#' rebuild: the name still resolves to all the rows.
#'
#' The matrix is already on the built term, and `term_build()` checked at
#' the full fit that it has one row per observation, so the rows of a fold are
#' the same rows by position. Binding the subset here builds, for the fold,
#' the spelling the documentation asks the caller for.
#'
#' Nothing is relearned that should be: a matrix carries no knots, no
#' contrasts and no levels, so subsetting it and re-evaluating it give the
#' same block. A FORMULA input is untouched and keeps being rebuilt on the
#' fold's own rows, which is what that rule exists for.
#'
#' It applies where `input_expr` is a plain symbol, which is the case the
#' name can be bound for. A call -- `lasso(scale(X))` -- keeps only its
#' own value on the term and not the `X` its re-evaluation would need, so
#' it is left to the error that names it.
#'
#' @param spec A [StatmodSpec()] whose terms are built.
#' @param sub The subset of the data, already taken.
#' @param i The rows it was taken with.
#' @param n The number of rows the fit was built on.
#'
#' @return `sub`, with a column per matrix input the terms carry.
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


#' The Held-Out Deviance of Every Point of a Path
#'
#' @description
#' Refits the model on each training fold along the whole path and scores it on
#' the fold left out, returning the mean deviance per observation and its
#' standard error across folds.
#'
#' @details
#' The path is run fold by fold rather than point by point so that each fit
#' starts from the previous point's coefficients, which is what makes a path
#' cheaper than its length suggests. Each training fit rebuilds the design on
#' its own rows: a term is re-evaluated in the data it is fitted to, so a basis
#' or a set of contrasts is not carried over from rows the fit did not see.
#'
#' @param spec A [StatmodSpec()].
#' @param data The data the fit was called on.
#' @param weights,offsets As [statmod()] received them.
#' @param inner_optimizer The inner method.
#' @param hypers A list of hyperparameter settings, one per path point.
#' @param folds A fold number per observation.
#' @param run Which combination of the outer axes each point belongs to. The
#'   warm start begins again at the head of each, the kink jumping back up
#'   there. `NULL` treats the whole list as one run.
#'
#' @return A list with `cvm`, `cvse` and `n_fail`.
#'
#' @seealso [cv()]
#'
#' @keywords internal
cv_curve <- function(spec, data, weights, offsets, inner_optimizer, hypers,
                     folds, run = NULL) {
  nf <- max(folds)
  m <- length(hypers)
  # One seed per fold, drawn HERE: a fold's fit is then a deterministic
  # function of the fold whether it runs in this process or in a worker,
  # which is what lets the result be independent of `workers` bit for bit
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
    w <- if (is.null(weights)) NULL else weights[keep]
    off <- if (is.null(offsets)) NULL else lapply(offsets, function(o) o[keep])
    # the design depends on the fold and not on the path, so it is built once
    # here rather than once per point. Measured, that is worth about 4 per cent
    # and not the most of it, which was the guess: at 200 observations and 20
    # columns a cross-validated path costs 0.88 seconds a fit and almost all of
    # it is the proximal iteration.
    # the fold is rebuilt with the SAME linpar options: one that built a
    # dense design where the fit built a sparse one would be paying for a
    # storage the model did not ask for
    ts <- tryCatch(statmod_spec(spec@formula, spec@distrib, train, w, off,
                                linpar = spec@linpar),
                   error = function(e) conditionMessage(e))
    # statmod_spec() builds a FRESH specification, so the thread count does
    # not travel with it as it does through statmod_respec(), which starts
    # from the one it is given: until 2026-08-21 a fold fell back to the
    # class default of 1 and a cross-validated path was single-threaded
    # whatever statmod(threads =) asked for. Measured on a lasso over a
    # gamma response, bic() gained 1.85x from eight threads where cv()
    # gained 1.06x. The two levels still do not nest -- where the folds
    # themselves run in worker PROCESSES each one fits on one thread, which
    # is what numericals7::n_threads() documents -- so the count is passed
    # on only when this process is running them.
    if (!is.character(ts)) {
      ts@threads <- if (spec@workers > 1L) 1L else spec@threads
    }
    if (is.character(ts)) {
      # a fold that cannot be rebuilt is a configuration error and not a
      # numerical one -- a term whose input is not in `data` cannot be
      # re-evaluated on a subset of it -- so the message is carried out
      # rather than dropped: every fold fails for the same reason, and the
      # path would otherwise score NA everywhere and keep its starting value
      return(list(dev = dev_f, err = ts))
    }
    td <- statmod_design(ts)
    tb <- statmod_blocks(ts, td)
    hs <- statmod_respec(ts, test)
    hd <- statmod_design(hs)
    cfgs <- inner_settings(inner_optimizer)
    obj <- statmod_objective(ts, hypers[[1L]], td, cfgs$expected, cfgs$approx)
    warm <- statmod_start(ts, td, obj, NULL)
    tbj <- tb
    warm0 <- warm
    for (j in seq_len(m)) {
      if (!is.null(run) && j > 1L && run[[j]] != run[[j - 1L]]) {
        warm <- warm0
        tbj <- tb
      }
      r <- tryCatch(statmod_alternate(ts, td, tbj, hypers[[j]], inner_optimizer,
                                      warm, cfgs$expected, cfgs$approx,
                                      cfgs$maxit, cfgs$tol, verbosity(0),
                                      hold_refresh = TRUE),
                    error = function(e) NULL)
      tbj <- blocks_at_kink(tb, hypers[[j]])
      if (is.null(r) || !isTRUE(r$converged)) next
      warm <- r$par
      cf <- r$obj$split(r$par)
      d <- tryCatch(-2 * statmod_loglik_at(hs, cf, hd),
                    error = function(e) NA_real_)
      if (is.finite(d)) dev_f[j] <- d / nrow(test)
    }
    list(dev = dev_f, err = NULL)
  }
  rows <- worker_map(spec, nf, one_fold)
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
#' Applies a body to each of `n` independent units -- a
#' cross-validation fold, a combination of a path's product grid -- over
#' the worker processes the specification asks for (`spec@workers`,
#' from [`n_threads(workers =)`][numericals7::n_threads]) and in this
#' process otherwise. Results come back in unit order whatever the number
#' of workers, which is what makes the answer independent of the count:
#' the units share nothing, so the same bodies run either way.
#'
#' @details
#' The units are independent BY CONSTRUCTION -- a fold is a complete refit
#' on its own rows, a path combination restarts its warm chain from the
#' sweep's own starting coefficients -- so they go by PROCESSES, with the
#' safeguards `optimizers7::multistart` records: under `pkgload`
#' the run stays sequential, because a worker loads the installed copy and
#' S7 objects built in the development namespace do not dispatch correctly
#' against it; a cluster that cannot start, or workers that cannot load
#' the package, fall back to sequential with a warning rather than fail
#' the fit. A fit inside a worker takes a fresh specification and is
#' therefore sequential by construction: the two levels of parallelism do
#' not nest.
#'
#' @param spec A [StatmodSpec()].
#' @param n How many units.
#' @param body The unit's body, a function of the unit index.
#' @param what The unit's name, for the warnings.
#'
#' @return A list of the bodies' results, in unit order.
#'
#' @seealso [cv_curve()], [statmod_path()]
#'
#' @keywords internal
worker_map <- function(spec, n, body, what = "folds") {
  workers <- min(spec@workers, n)
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


#' Choose a Point of the Path
#'
#' @description
#' The index of the smallest criterion, or of the largest kink whose criterion
#' is within one standard error of it.
#'
#' @details
#' The one-standard-error rule takes the sparsest fit that is not measurably
#' worse than the best one, which is Breiman's rule as \pkg{glmnet} applies it.
#' The path runs from the emptiest fit to the fullest, so the largest kink
#' among the admissible points is the first of them.
#'
#' @param value The criterion at each point.
#' @param se Its standard error, or `NULL`.
#' @param rule `"min"` or `"1se"`.
#'
#' @return A single index, or `NA` where no point was usable.
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
#' Sweeps them over grids of kink sizes and keeps the setting the criterion
#' prefers.
#'
#' @details
#' The grid runs from the kink that empties the block down to `min_ratio`
#' of it, so the sweep goes from the sparsest fit towards the fullest and every
#' fit starts from the previous one's coefficients. Where the top of the grid
#' does not empty the block it is doubled until it does, the starting value
#' being computed at the coefficients in hand rather than at a refitted null.
#'
#' A term carrying several of them has every combination visited where the
#' term asks for `search = "grid"` and one coordinate at a time where it
#' asks for `"cyclic"`. Between terms the alternation is cyclic either
#' way, so the cost is the product WITHIN a term and the sum ACROSS them. Each
#' axis is built at the settings of the axes outside it, which is what makes
#' the elastic net's grid a family of \eqn{\lambda} axes rather than one, and
#' the axis swept by kink size is put innermost so that the warm starts walk
#' along it.
#'
#' A pass that would visit the points the last one scored is not run. The top
#' of the path is read again at every pass because the rest of the model moves,
#' and where it has not the grid is the one already in hand.
#'
#' Where the model also carries hyperparameters that are twice differentiable,
#' those are estimated by [outer_fit()] inside each point of the
#' path, so the two kinds are not mixed into one search.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param blocks The blocks.
#' @param hyper The hyperparameters.
#' @param inner_optimizer The inner method.
#' @param method An [OuterMethod()].
#' @param optimizer The optimizer for the differentiable hyperparameters.
#' @param beta The starting coefficients.
#' @param approx The approximation for the expected information.
#' @param maxit,tol The inner budget.
#' @param vb The verbosity.
#' @param data,weights,offsets What the fit was called on, which
#'   cross-validation needs in order to refit on part of it.
#' @param rows Which hyperparameters to select.
#' @param sweeps How many cyclic passes.
#'
#' @return The same list [outer_fit()] returns.
#'
#' @seealso [cv()], [path_rows()]
#'
#' @keywords internal
statmod_path <- function(spec, design, blocks, hyper, inner_optimizer, method,
                         optimizer, beta, approx, maxit, tol, vb, data,
                         weights, offsets, rows, sweeps = 2L,
                         nested_method = NULL) {
  expected <- identical(method@hessian, "expected")
  smooth_idx <- outer_hyper_index(spec, blocks)
  is_cv <- identical(method@kind, "cv")
  # the smooth hyperparameters are estimated inside each point of the path by
  # THEIR OWN criterion, which is not the path's: a smoothing parameter is
  # read at the mode by reml() while a lasso's lambda is swept by bic()
  inner_crit <- if (is.null(nested_method)) method else nested_method
  nested <- nrow(smooth_idx) > 0L && !is.null(inner_crit) &&
    !identical(inner_crit@kind, "cv")
  obj0 <- statmod_objective(spec, hyper, design, expected, approx)

  # one fit at one setting, warm-started, with the differentiable
  # hyperparameters estimated inside it where there are any
  # `sp` is the specification a fit runs under, and it is an argument
  # because one caller needs a different one: a run dispatched to a worker
  # PROCESS must fit on a single thread, the two levels of
  # numericals7::n_threads() not nesting, while every other call here runs
  # in this process and uses the count the caller asked for.
  fit_at <- function(hy, warm, bk = blocks, sp = spec) {
    if (nested) {
      r <- tryCatch(outer_fit(sp, design, bk, hy, inner_optimizer,
                              inner_crit, optimizer, warm, approx, maxit, tol,
                              vb_inner(vb)), error = function(e) NULL)
      if (!is.null(r)) return(r)
    }
    r <- statmod_alternate(sp, design, bk, hy, inner_optimizer, warm,
                           expected, approx, maxit, tol, vb_inner(vb),
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
    act <- statmod_active(spec, blocks, r$par, at)
    m <- statmod_pe(spec, design, cf, at, method, approx, act)
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
      ob <- statmod_objective(spec, hy, design, expected, approx)
      tp[[i]] <- path_null_score(ob, at, b, hy)
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
        if (max(abs(r$par[b$index])) <= 1e-8) break
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
          cc <- cv_curve(spec, data, weights, offsets, inner_optimizer, hys,
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
          if (spec@workers > 1L) run_spec@threads <- 1L
          one_run <- function(ri) {
            idxs <- runs[[ri]]
            warm <- beta
            bk <- blocks
            out <- rep(NA_real_, length(idxs))
            for (k in seq_along(idxs)) {
              j <- idxs[[k]]
              r <- fit_at(hys[[j]], warm, bk, sp = run_spec)
              bk <- blocks_at_kink(blocks, hys[[j]])
              if (!isTRUE(r$converged)) next
              warm <- r$par
              out[[k]] <- score_at(r, r$hyper)
            }
            out
          }
          got <- worker_map(spec, length(runs), one_run,
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
            " '%s',\n  so %s was not chosen and would have kept the value it",
            " came in with.%s"),
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
        # end is the interval's and there is nothing to widen.
        if (pl$kink && identical(method@rule, "min") &&
            length(mine) > 1L && k %in% c(1L, length(mine)) && falling()) {
          warning(sprintf(paste0("The path for '%s' in '%s' stopped at its %s ",
                                 "end (%s = %s).\n  The criterion was still ",
                                 "falling there, so widen the path with ",
                                 "min_ratio\n  or set the value yourself."),
                          row$term, row$parameter,
                          if (k == 1L) "sparse" else "dense", row$name,
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

  final <- fit_at(cur, beta)
  cf <- final$obj$split(final$par)
  act <- statmod_active(spec, blocks, final$par, final$hyper)
  crit <- if (is_cv) best$value else {
    m <- statmod_pe(spec, design, cf, final$hyper, method, approx, act)
    if (is.null(m)) NA_real_ else m$value
  }
  list(par = final$par, hyper = final$hyper, value = final$value,
       criterion = crit, converged = isTRUE(final$converged),
       obj = final$obj, hist_blocks = final$hist_blocks,
       hist_inner = final$hist_inner,
       hist_outer = do.call(rbind, hist),
       iterations = length(hist), evaluations = nrow(do.call(rbind, hist)),
       exact_gradient = FALSE, exact_hessian = FALSE)
}


#' The Block a Path Row Belongs To
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
  stop(sprintf("No penalized block for '%s' in '%s'.", row$term,
               row$parameter), call. = FALSE)
}


#' Choose the Hyperparameters by Cross-Validation
#'
#' @description
#' Builds the object that tells [statmod()] to choose a model's
#' hyperparameters by cross-validation: a path of values is scored by the
#' log-likelihood the fit assigns to observations it was not fitted on, and
#' the best is kept. Pass the result as `sparse_criterion`, or as
#' `outer_criterion`, or as both.
#'
#' @details
#' # Why a penalty with a kink needs this
#'
#' A marginal criterion, [reml()] or [ml()], approximates an integral by a
#' Laplace expansion at the penalized mode, which asks for the second
#' derivative of the penalty there. The mode of a lasso, a SCAD or an MCP
#' sits at the kink for every coefficient it sets to zero, and that is where
#' the derivative does not exist, so the criterion is undefined at the point
#' it would be read at.
#'
#' Cross-validation asks a different question, about prediction and not about
#' a posterior, and asking it needs nothing from the penalty beyond a fit.
#'
#' # The criterion
#'
#' The mean over folds of \eqn{-2\ell/n_f} on the fold left out, each
#' training fit rebuilding the design on its own rows.
#'
#' `rule = "1se"` takes the largest kink whose criterion is within one
#' standard error of the smallest, which is the sparsest fit that is not
#' measurably worse.
#'
#' # A path, not a search
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
#' The term's answer: every one of its penalty's that the constructor did not
#' hold at a number. What the criterion decides is how they are covered.
#'
#' # A product within a term, an alternation between terms
#'
#' A term carrying several kinked hyperparameters has every combination
#' visited under `search = "grid"` and one swept at a time under
#' `"cyclic"`.
#'
#' The choice belongs to the term, through `enet(X, search =)` and
#' [modelterms7::term_search()], and not to this criterion. The same
#' criterion is put to the model's smooth hyperparameters as well, and those
#' are read at the mode instead of being swept, so most of what it is asked
#' about could not use such an argument.
#'
#' Between two terms the sweep alternates whichever each one named, so
#' `y ~ lasso(X) + enet(R)` costs the two blocks added and not multiplied,
#' and one term asking for a product does not make the other pay for it.
#'
#' A term that names neither gets the product, because the cyclic sweep
#' traverses a cross through the point in hand and can stop where each
#' coordinate is separately best without being jointly so. Its cost is the
#' product of the term's grids where the cyclic sweep's is their sum, which at
#' two hyperparameters is `n_lambda * n_alpha` fits against
#' `n_lambda + n_alpha` per pass; with three or more estimated it grows
#' exponentially, and `"cyclic"` is there for that.
#'
#' # How long the path is, and how far down it reaches
#'
#' The term's answers again, written on its own signature where a reader can
#' see them: `lasso(x, n_lambda = 25, min_ratio = 1e-4)`,
#' `enet(x, n_lambda = 25, n_alpha = 5)`.
#'
#' The two numbers differ because the axes do. \eqn{\lambda} descends the
#' size of the kink over four decades, while \eqn{\alpha} spans one bounded
#' interval, so the shipped product is 25 by 5.
#'
#' # The grid is not a rectangle
#'
#' For two different reasons. The
#' elastic net's kink is \eqn{\lambda\alpha}, so the value emptying the block
#' is \eqn{\lambda_{\max} = \kappa/\alpha} and every \eqn{\alpha} carries its
#' own \eqn{\lambda} axis, descending from its own top. The shapes of SCAD and
#' MCP leave the kink alone, so there \eqn{\lambda_{\max}} is one number
#' whatever the shape and the two axes really are a rectangle. Each axis is
#' built at the settings of the axes outside it, which covers both cases with
#' no rule about families written down anywhere.
#'

#' # The cost
#'
#' `nfolds` fits per point of the path, and how many points there are is the
#' term's `n_lambda` for one kinked hyperparameter, or the product of its
#' axes for several.
#'
#' The warm starts are worth 1.8 times, and building each fold's design once
#' instead of once per point another 4 per cent. What remains is the proximal
#' iteration: measured at 200 observations and 20 columns, 0.88 seconds a
#' fit, against `cv.glmnet`'s 0.03 seconds for its whole path of 100 values
#' on five folds. That distance is why `n_lambda` is 25 here and 100 there.
#'
#' The folds themselves run over separate R processes when
#' `statmod(threads = n_threads(workers = 4))` asks for it, with one seed
#' drawn per fold in the parent, so the answer is identical at any worker
#' count. Measured, a ten-fold lasso cross-validation goes from 31.45 s to
#' 14.18 s at four workers.
#'
#' @param nfolds How many folds, a single number of at least 2. `10` by
#'   default. Ignored when `folds` is given.
#' @param folds A fold number per observation, an integer vector as long as
#'   the data, for a partition of your own or to compare two criteria on the
#'   same one. `integer(0)`, the default, draws them at fit time.
#' @param rule `"min"` for the best criterion, `"1se"` for the sparsest fit
#'   within one standard error of it. Matched with [match.arg()].

#' @return An [OuterMethod()] object of kind `"cv"`, carrying `nfolds`,
#'   `rule` and `folds`. Its `hessian` is `"observed"` and its `k` is
#'   `NA_real_`, neither being read by this criterion.
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
#'                outer_criterion = cv(nfolds = 5))
#' hyper(fit)
#'
#' # The two columns that carry signal are recovered.
#' round(coef(fit)$mu[2:3], 3)
#'
#' # "1se" takes the largest kink within one standard error of the best, so
#' # it never keeps more coefficients than "min". On this data the two agree,
#' # the criterion being flat enough that the best point is already within
#' # one standard error of itself.
#' one_se <- statmod(y ~ lasso(x, n_lambda = 12),
#'                   distributions7::gaussian1_distrib(), dd,
#'                   outer_criterion = cv(nfolds = 5, rule = "1se"))
#' c(min = hyper(fit)$estimate, one_se = hyper(one_se)$estimate)
#' c(min = sum(coef(fit)$mu[-1] != 0), one_se = sum(coef(one_se)$mu[-1] != 0))
#'
#' @export
cv <- function(nfolds = 10, folds = NULL, rule = c("min", "1se")) {
  OuterMethod(kind = "cv", hessian = "observed", k = NA_real_,
              nfolds = as.numeric(nfolds), rule = match.arg(rule),
              folds = if (is.null(folds)) numeric(0) else as.numeric(folds))
}


#' Estimate the Hyperparameters, by Whichever Route Each One Admits
#'
#' @description
#' Routes the twice differentiable hyperparameters to [outer_fit()]
#' and the rest to [statmod_path()].
#'
#' @details
#' The split is the same one that decides how the coefficients are fitted, so a
#' term whose penalty has a kink has both its coefficients and its
#' hyperparameters handled by methods that do not ask for a curvature it does
#' not have.
#'
#' @inheritParams statmod_path
#'
#' @return The list [outer_fit()] returns.
#'
#' @seealso [statmod()]
#'
#' @keywords internal
statmod_select <- function(spec, design, blocks, hyper, inner_optimizer, method,
                           optimizer, beta, approx, maxit, tol, vb, data,
                           weights, offsets, sparse_method = NULL) {
  # Two families of penalty, two criteria. A kinked one is swept over a path
  # of its own values by `sparse_method`; a twice differentiable one is read
  # by `method` at the mode. Where both are present the path is outside and
  # the marginal criterion is estimated inside each of its points, so a model
  # can have its smoothing parameter by REML and its lasso by BIC at once --
  # which one argument could not express.
  sm <- if (is.null(sparse_method)) method else sparse_method
  rows <- path_rows(spec, blocks, hyper, sm)
  if (!nrow(rows)) {
    if (is.null(method)) {
      # held here too: statmod() runs the full alternation, phase on, at
      # whatever the selection ends with, so every select route refines the
      # break-points exactly once
      return(statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                               beta, identical(sm@hessian, "expected"),
                               approx, maxit, tol, vb, hold_refresh = TRUE))
    }
    if (identical(method@kind, "cv")) {
      stop(paste0("cv() has nothing to select: no term here carries a penalty",
                  "\n  with a kink. A smoothing parameter that is twice",
                  " differentiable is\n  estimated by reml(), ml(), aic() or",
                  " bic(), which read a criterion\n  rather than refitting on",
                  " folds."), call. = FALSE)
    }
    return(outer_fit(spec, design, blocks, hyper, inner_optimizer, method,
                     optimizer, beta, approx, maxit, tol, vb))
  }
  statmod_path(spec, design, blocks, hyper, inner_optimizer, sm, optimizer,
               beta, approx, maxit, tol, vb, data, weights, offsets, rows,
               nested_method = method)
}
