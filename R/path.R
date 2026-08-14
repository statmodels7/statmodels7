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
#' @return A single non-negative number, \code{0} where there is no kink.
#'
#' @seealso \code{\link{kink_hypers}}, \code{\link{kink_solve}}
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
#' The names whose value moves \code{\link{kink_scale}}, which are the ones a
#' path over the penalty has to vary.
#'
#' @details
#' The question is put to the penalty rather than answered from a list of
#' families: the shape parameters of SCAD and MCP leave the subdifferential at
#' zero unchanged and govern how fast the penalty flattens further out, while
#' \eqn{\lambda} and the elastic net's \eqn{\alpha} both scale it.
#'
#' \code{unbounded} restricts the answer to the hyperparameters with no upper
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
#' @seealso \code{\link{kink_scale}}
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
#' Solves \code{kink_scale(pen, theta) == target} in one named hyperparameter.
#'
#' @details
#' The size of the kink is monotone in such a hyperparameter but not
#' necessarily increasing: a Laplace prior written by its scale has a kink of
#' \eqn{1/\sigma}, which narrows as the hyperparameter grows. Which way to walk
#' is therefore measured before walking, by comparing the size at the current
#' value with the size at twice it, and only then is the root bracketed by
#' doubling and found with \code{\link[stats]{uniroot}}. A version that assumed
#' the size increases returned \code{NA} for every target on the Laplace,
#' having walked away from the answer.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The other hyperparameters.
#' @param name Which one to solve for.
#' @param target The size the kink should have.
#'
#' @return A single value, or \code{NA} where the target is out of reach.
#'
#' @keywords internal
kink_solve <- function(pen, theta, name, target) {
  th <- as.list(theta)
  b <- pen@params_bounds[[name]]
  size <- function(v) {
    th[[name]] <- v
    kink_scale(pen, th)
  }
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
#' \pkg{glmnet}'s \code{lambda.max} written for any separable penalty.
#'
#' The other coefficients are held where the caller left them rather than
#' refitted, so the number is a starting point and not a boundary. The path
#' checks it: a top whose fit is not empty is doubled until it is.
#'
#' @param obj The stacked objective.
#' @param beta The current coefficients.
#' @param block One entry of \code{statmod_blocks()$sparse}.
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
#' \code{min_ratio} of it, carried back onto the hyperparameter.
#'
#' @details
#' The grid is geometric in the size of the kink rather than in the
#' hyperparameter, so that a penalty whose kink narrows as its hyperparameter
#' grows is swept in the same order as one whose kink widens: from the empty
#' model towards the full one. Values the penalty cannot reach are dropped.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which one the path varies.
#' @param s_max The size of the kink at the top of the path.
#' @param n_values How many points.
#' @param min_ratio The smallest kink size, as a fraction of \code{s_max}.
#'
#' @return A numeric vector of values for \code{name}, from the emptiest fit to
#'   the fullest.
#'
#' @keywords internal
path_values <- function(pen, theta, name, s_max, n_values = 40L,
                        min_ratio = 1e-4) {
  s <- exp(seq(log(s_max), log(s_max * min_ratio), length.out = n_values))
  v <- vapply(s, function(target) kink_solve(pen, theta, name, target),
              numeric(1))
  v[is.finite(v) & v > 0]
}


#' The Values a Path Visits Over a Bounded Hyperparameter
#'
#' @description
#' An equally spaced grid strictly inside the hyperparameter's own interval.
#'
#' @details
#' \code{\link{path_values}} walks the SIZE OF THE KINK, from the value that
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
#' @param pen A \pkg{penalties7} penalty.
#' @param name Which hyperparameter the path varies.
#' @param n_values How many points.
#'
#' @return A numeric vector of values for \code{name}.
#'
#' @seealso \code{\link{path_values}}, \code{\link{path_bounded}}
#'
#' @keywords internal
path_grid <- function(pen, name, n_values = 25L) {
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
  # geometric above the lower bound. It spans the values the literature uses,
  # 3.7 for the SCAD of Fan and Li and 3 for the MCP of Zhang, and cannot
  # reach the bound itself, where neither penalty is defined.
  b[1L] + exp(seq(log(0.25), log(25), length.out = n))
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
#' \code{\link{path_grid}}.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param theta The hyperparameters in force.
#' @param name Which hyperparameter.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{path_values}}, \code{\link{path_grid}}
#'
#' @keywords internal
path_by_kink <- function(pen, theta, name) {
  !path_bounded(pen, name) &&
    name %in% kink_hypers(pen, theta, unbounded = FALSE)
}


#' Which Hyperparameters a Path Has to Select
#'
#' @description
#' The rows of an index like \code{\link{outer_hyper_index}}'s, for the
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
#' \code{\link{kink_hypers}} unless the method names them.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param blocks The blocks.
#' @param hyper The hyperparameters.
#' @param method An \code{\link{OuterMethod}}.
#'
#' @return A data frame of \code{parameter}, \code{term} and \code{name}.
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
#' Passing \code{folds} explicitly is what makes two criteria comparable on the
#' same partition.
#'
#' @param n The number of observations.
#' @param k How many folds.
#' @param folds The method's own assignment, or \code{numeric(0)}.
#'
#' @return An integer vector of length \code{n}.
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
#' @param row One row of \code{\link{path_rows}}'s index.
#' @param value The value to write.
#'
#' @return The hyperparameters, with that one changed.
#'
#' @keywords internal
hyper_set <- function(hyper, row, value) {
  hyper[[row$parameter]][[row$term]][[row$name]] <- value
  hyper
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
#' @param spec A \code{\link{StatmodSpec}}.
#' @param data The data the fit was called on.
#' @param weights,offsets As \code{\link{statmod}} received them.
#' @param inner_optimizer The inner method.
#' @param hypers A list of hyperparameter settings, one per path point.
#' @param folds A fold number per observation.
#'
#' @return A list with \code{cvm}, \code{cvse} and \code{n_fail}.
#'
#' @seealso \code{\link{cv}}
#'
#' @keywords internal
cv_curve <- function(spec, data, weights, offsets, inner_optimizer, hypers,
                     folds) {
  nf <- max(folds)
  m <- length(hypers)
  dev <- matrix(NA_real_, nf, m)
  err <- NULL
  for (f in seq_len(nf)) {
    keep <- folds != f
    train <- data[keep, , drop = FALSE]
    test <- data[!keep, , drop = FALSE]
    w <- if (is.null(weights)) NULL else weights[keep]
    off <- if (is.null(offsets)) NULL else lapply(offsets, function(o) o[keep])
    # the design depends on the fold and not on the path, so it is built once
    # here rather than once per point. Measured, that is worth about 4 per cent
    # and not the most of it, which was the guess: at 200 observations and 20
    # columns a cross-validated path costs 0.88 seconds a fit and almost all of
    # it is the proximal iteration.
    ts <- tryCatch(statmod_spec(spec@formula, spec@distrib, train, w, off),
                   error = function(e) conditionMessage(e))
    if (is.character(ts)) {
      # a fold that cannot be rebuilt is a configuration error and not a
      # numerical one -- a term whose input is not in `data` cannot be
      # re-evaluated on a subset of it -- so the message is carried out
      # rather than dropped: every fold fails for the same reason, and the
      # path would otherwise score NA everywhere and keep its starting value
      err <- ts
      next
    }
    td <- statmod_design(ts)
    tb <- statmod_blocks(ts, td)
    hs <- statmod_respec(ts, test)
    hd <- statmod_design(hs)
    cfgs <- inner_settings(inner_optimizer)
    obj <- statmod_objective(ts, hypers[[1L]], td, cfgs$expected, cfgs$approx)
    warm <- statmod_start(ts, td, obj, NULL)
    tbj <- tb
    for (j in seq_len(m)) {
      r <- tryCatch(statmod_alternate(ts, td, tbj, hypers[[j]], inner_optimizer,
                                      warm, cfgs$expected, cfgs$approx,
                                      cfgs$maxit, cfgs$tol, verbosity(0)),
                    error = function(e) NULL)
      tbj <- blocks_at_kink(tb, hypers[[j]])
      if (is.null(r) || !isTRUE(r$converged)) next
      warm <- r$par
      cf <- r$obj$split(r$par)
      d <- tryCatch(-2 * statmod_loglik_at(hs, cf, hd),
                    error = function(e) NA_real_)
      if (is.finite(d)) dev[f, j] <- d / nrow(test)
    }
  }
  ok <- colSums(is.finite(dev))
  cvm <- apply(dev, 2L, mean, na.rm = TRUE)
  cvse <- apply(dev, 2L, stats::sd, na.rm = TRUE) / sqrt(pmax(ok, 1))
  cvm[ok == 0L] <- NA_real_
  list(cvm = cvm, cvse = cvse, n_fail = nf - ok, error = err)
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
#' @param se Its standard error, or \code{NULL}.
#' @param rule \code{"min"} or \code{"1se"}.
#'
#' @return A single index, or \code{NA} where no point was usable.
#'
#' @references
#' Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
#' \emph{Classification and Regression Trees}. Wadsworth.
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
#' Sweeps each of them over a grid of kink sizes, holding the others, and keeps
#' the setting the criterion prefers.
#'
#' @details
#' The grid runs from the kink that empties the block down to \code{min_ratio}
#' of it, so the sweep goes from the sparsest fit towards the fullest and every
#' fit starts from the previous one's coefficients. Where the top of the grid
#' does not empty the block it is doubled until it does, the starting value
#' being computed at the coefficients in hand rather than at a refitted null.
#'
#' With several such hyperparameters the sweeps are cyclic, one coordinate at a
#' time, which is what keeps the cost linear in their number where a full grid
#' would be exponential in it.
#'
#' Where the model also carries hyperparameters that are twice differentiable,
#' those are estimated by \code{\link{outer_fit}} inside each point of the
#' path, so the two kinds are not mixed into one search.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param blocks The blocks.
#' @param hyper The hyperparameters.
#' @param inner_optimizer The inner method.
#' @param method An \code{\link{OuterMethod}}.
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
#' @return The same list \code{\link{outer_fit}} returns.
#'
#' @seealso \code{\link{cv}}, \code{\link{path_rows}}
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
  fit_at <- function(hy, warm, bk = blocks) {
    if (nested) {
      r <- tryCatch(outer_fit(spec, design, bk, hy, inner_optimizer,
                              inner_crit, optimizer, warm, approx, maxit, tol,
                              vb_inner(vb)), error = function(e) NULL)
      if (!is.null(r)) return(r)
    }
    r <- statmod_alternate(spec, design, bk, hy, inner_optimizer, warm,
                           expected, approx, maxit, tol, vb_inner(vb))
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
  top <- stats::setNames(numeric(nrow(rows)), rownames(rows))
  for (i in seq_len(nrow(rows))) {
    b <- path_block(blocks, rows[i, ])
    top[[i]] <- path_null_score(obj0, beta, b, hyper)
  }

  # path_null_score() reads the score with the OTHER coefficients where the
  # caller left them rather than at their optimum with this block zeroed, so
  # it is a starting point and not a boundary -- and measured, it can be an
  # order of magnitude short: on eight coefficients of which three carried
  # real signal it returned 26.5 where the block first empties near 500, so
  # the path covered a nearly flat stretch of the criterion and every fit
  # chose its sparse end. Doubling until the fit really is empty is the check
  # the documentation already described, and it costs a handful of fits.
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    b <- path_block(blocks, row)
    # a hyperparameter that does not scale the kink has no top of this kind:
    # no value of it empties the block
    if (!path_by_kink(b$penalty, hyper[[row$parameter]][[row$term]],
                      row$name)) next
    for (step in seq_len(24L)) {
      v <- kink_solve(b$penalty, hyper[[row$parameter]][[row$term]],
                      row$name, top[[i]])
      if (!is.finite(v) || v <= 0) break
      r <- tryCatch(fit_at(hyper_set(hyper, row, v), beta),
                    error = function(e) NULL)
      if (is.null(r)) break
      if (max(abs(r$par[b$index])) <= 1e-8) break
      top[[i]] <- top[[i]] * 2
    }
  }

  cur <- hyper
  hist <- list()
  best <- list(value = Inf, hyper = cur, fit = NULL)
  for (s in seq_len(if (nrow(rows) > 1L) sweeps else 1L)) {
    moved <- FALSE
    for (i in seq_len(nrow(rows))) {
      row <- rows[i, ]
      b <- path_block(blocks, row)
      # the kink-size path is for a hyperparameter that SCALES the kink and
      # has no upper bound; everything else is swept over a grid of its own
      bounded <- !path_by_kink(b$penalty, cur[[row$parameter]][[row$term]],
                               row$name)
      # HOW FINELY is the term's answer where it gave one: a block of four
      # columns and one of four hundred want different grids, and the
      # criterion applies to every term at once and cannot know which it is
      # looking at. Where the term says nothing the criterion's default
      # stands.
      nv <- statmod_grid_size(spec, row, as.integer(method@n_values))
      mr <- statmod_min_ratio(spec, row, method@min_ratio)
      vals <- if (bounded) {
        path_grid(b$penalty, row$name, nv)
      } else {
        path_values(b$penalty, cur[[row$parameter]][[row$term]],
                    row$name, top[[i]], nv, mr)
      }
      if (!length(vals)) next
      hys <- lapply(vals, function(v) hyper_set(cur, row, v))

      cv_err <- NULL
      if (is_cv) {
        folds <- cv_folds(spec@n_obs, method@nfolds, method@folds)
        cc <- cv_curve(spec, data, weights, offsets, inner_optimizer, hys, folds)
        value <- cc$cvm
        se <- cc$cvse
        cv_err <- cc$error
      } else {
        warm <- beta
        value <- rep(NA_real_, length(hys))
        # the point just fitted is what the next one screens against: the grid
        # runs from the emptiest fit towards the fullest, so the kink shrinks
        # and the strong rule has a previous size to compare with
        bk <- blocks
        for (j in seq_along(hys)) {
          r <- fit_at(hys[[j]], warm, bk)
          bk <- blocks_at_kink(blocks, hys[[j]])
          if (!isTRUE(r$converged)) next
          warm <- r$par
          value[[j]] <- score_at(r, r$hyper)
        }
        se <- NULL
      }

      hist[[length(hist) + 1L]] <- data.frame(
        sweep = s, parameter = row$parameter, term = row$term,
        name = row$name, value = vals, criterion = value,
        se = if (is.null(se)) NA_real_ else se, stringsAsFactors = FALSE)

      # NOT A POINT OF THE PATH WAS SCORED. Skipping quietly leaves the
      # hyperparameter at whatever it came in with -- its default, where the
      # caller asked for it to be chosen -- and the fit then reports success
      # at a value nothing selected, which is the failure this refuses to
      # produce. Where the criterion said why, it says so here.
      # ... and only where the criterion is one that selects. A marginal
      # criterion reaching a kinked row scores nothing BY CONSTRUCTION -- the
      # Laplace expansion it is has no second derivative at the kink -- and
      # the hyperparameter keeping the value it was given is the documented
      # answer there, not a failure.
      if (method@kind %in% c("aic", "bic", "cv") && !any(is.finite(value))) {
        stop(sprintf(paste0(
          "%s() could not score a single point of the path for '%s' in '%s',",
          " so\n  %s was not chosen and would have kept the value it came in",
          " with.%s"),
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
      falling <- function(j) {
        k <- if (j == 1L) 2L else length(value) - 1L
        if (k < 1L || k > length(value)) return(FALSE)
        v <- value[c(j, k)]
        if (!all(is.finite(v))) return(FALSE)
        v[1L] < v[2L] - 1e-8 * max(1, abs(v[2L]))
      }
      # ... and only for a hyperparameter swept by kink size. A bounded one
      # is swept over the whole of its own interval, so an answer at either
      # end is the interval's and there is nothing to widen.
      if (!bounded && identical(method@rule, "min") &&
          j %in% c(1L, length(value)) && falling(j)) {
        warning(sprintf(paste0("The path for '%s' in '%s' stopped at its %s ",
                               "end (%s = %s).\n  The criterion was still ",
                               "falling there, so widen the path with ",
                               "min_ratio\n  or set the value yourself."),
                        row$term, row$parameter,
                        if (j == 1L) "sparse" else "dense", row$name,
                        format(signif(vals[[j]], 4))), call. = FALSE)
      }
      if (!identical(cur[[row$parameter]][[row$term]][[row$name]], vals[[j]])) {
        moved <- TRUE
      }
      cur <- hyper_set(cur, row, vals[[j]])
      best$value <- value[[j]]
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
#' @param row One row of \code{\link{path_rows}}'s index.
#'
#' @return One entry of \code{blocks$sparse}.
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
#' \code{cv()} scores a path of hyperparameter values by the log-likelihood the
#' fit assigns to observations it was not fitted on.
#'
#' @details
#' \strong{Why a penalty with a kink needs this.} A marginal criterion --
#' \code{\link{reml}()} or \code{\link{ml}()} -- approximates an integral by a
#' Laplace expansion at the penalized mode, which asks for the second
#' derivative of the penalty there. The mode of a lasso, a SCAD or an MCP sits
#' at the kink for every coefficient it sets to zero, which is where that
#' derivative does not exist, so the criterion is not defined at the point it
#' would be read at. Cross-validation asks a different question, about
#' prediction rather than about a posterior, and asking it needs nothing from
#' the penalty beyond a fit.
#'
#' \strong{The criterion} is the mean over folds of
#' \eqn{-2\ell/n_f} on the fold left out, each training fit rebuilding the
#' design on its own rows. \code{rule = "1se"} takes the largest kink whose
#' criterion is within one standard error of the smallest, which is the
#' sparsest fit that is not measurably worse.
#'
#' \strong{The path, not a search.} The penalized mode is a piecewise smooth
#' function of the hyperparameter: differentiable while the active set holds,
#' turning a corner whenever a coefficient joins it or leaves. A criterion read
#' there inherits the corners, so the hyperparameter is swept over a grid
#' rather than searched by slope. The grid is geometric in the size of the
#' kink, from the value that leaves the block empty down to \code{min_ratio} of
#' it, and every fit starts from the previous one's coefficients.
#'
#' \strong{Which hyperparameters.} Those whose value sets the size of the kink,
#' read from the penalty by probing the subdifferential rather than taken from
#' a list of families: \eqn{\lambda} for the lasso, SCAD and MCP, and the
#' elastic net's \eqn{\lambda} with \eqn{\alpha} held, as \pkg{glmnet} holds it
#' and \pkg{ncvreg} holds \eqn{\gamma}. Name others in \code{over} to sweep
#' them too.
#'
#' \strong{The cost} is \code{nfolds * n_values} fits per hyperparameter. The
#' warm starts are worth 1.8 times, and building each fold's design once rather
#' than once per point another 4 per cent, but what remains is the proximal
#' iteration: measured at 200 observations and 20 columns, 0.88 seconds a fit,
#' against \code{cv.glmnet}'s 0.03 seconds for its whole path of 100 values on
#' five folds. That distance is the reason \code{n_values} is 25 here and 100
#' there, and closing it needs the compiled coordinate descent that a separable
#' penalty on a linear predictor admits.
#'
#' @param nfolds How many folds. Ignored when \code{folds} is given.
#' @param folds A fold number per observation, for a partition of your own or
#'   to compare two criteria on the same one.
#' @param rule \code{"min"} for the best criterion, \code{"1se"} for the
#'   sparsest fit within one standard error of it.
#' @param n_values How many points the path visits.
#' @param min_ratio The smallest kink the path reaches, as a fraction of the
#'   one that empties the block.
#'
#' @return An \code{\link{OuterMethod}}.
#'
#' @references
#' Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
#' \emph{Classification and Regression Trees}. Wadsworth.
#'
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. \emph{Journal of
#' Statistical Software} 33(1), 1--22.
#'
#' @seealso \code{\link{reml}}, \code{\link{aic}}, \code{\link{statmod}}
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(y = rnorm(60))
#' dd$x <- matrix(rnorm(60 * 5), 60, 5)
#' statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dd,
#'         outer_criterion = cv(nfolds = 3, n_values = 6))
#'
#' @export
cv <- function(nfolds = 10, folds = NULL, rule = c("min", "1se"),
               n_values = 25, min_ratio = 1e-4) {
  OuterMethod(kind = "cv", hessian = "observed", k = NA_real_,
              n_values = as.numeric(n_values),
              min_ratio = as.numeric(min_ratio),
              nfolds = as.numeric(nfolds), rule = match.arg(rule),
              folds = if (is.null(folds)) numeric(0) else as.numeric(folds))
}


#' Estimate the Hyperparameters, by Whichever Route Each One Admits
#'
#' @description
#' Routes the twice differentiable hyperparameters to \code{\link{outer_fit}}
#' and the rest to \code{\link{statmod_path}}.
#'
#' @details
#' The split is the same one that decides how the coefficients are fitted, so a
#' term whose penalty has a kink has both its coefficients and its
#' hyperparameters handled by methods that do not ask for a curvature it does
#' not have.
#'
#' @inheritParams statmod_path
#'
#' @return The list \code{\link{outer_fit}} returns.
#'
#' @seealso \code{\link{statmod}}
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
      return(statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                               beta, identical(sm@hessian, "expected"),
                               approx, maxit, tol, vb))
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
