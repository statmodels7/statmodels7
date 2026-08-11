#' @include path.R
NULL

#' Fit a Separable Block by Coordinate Descent
#'
#' @description
#' Estimates one penalized block by cycling over its coefficients on the
#' working quadratic of its own equation, the other blocks held fixed.
#'
#' @details
#' \strong{Why not the proximal method.} A proximal gradient step reads the
#' whole model: measured on 200 observations and 20 columns, one block fit made
#' 88 evaluations of the objective, 75 of the gradient and 83 of the operator,
#' each over every parameter of the distribution, and closed in 36 iterations
#' at 0.17 seconds. A coordinate descent reads the block's own columns and the
#' running residual instead and closes in six sweeps.
#'
#' \strong{The working quadratic.} With \eqn{\eta} the equation's linear
#' predictor, \eqn{s_i} the score in it and \eqn{h_i} the information,
#' \eqn{-\ell} is \eqn{\frac12\sum_i h_i(z_i - \eta_i)^2} up to a constant with
#' \eqn{z = \eta + s/h}, which is the weighted least squares problem of
#' \code{\link{iwls}} restricted to one equation. The other columns of that
#' equation enter as an offset. For a Gaussian response with an identity link
#' the quadratic is exact and one pass is the answer; otherwise the weights are
#' rebuilt and the sweeps repeated.
#'
#' \strong{The penalty arrives as a table.} The coordinate update is the
#' penalty's own proximal operator at the step \eqn{1/v_j}, with
#' \eqn{v_j = \sum_i w_i x_{ij}^2}, and \eqn{v_j} does not move while the
#' working weights are held. The whole table is therefore built once per
#' weighted least squares iteration by
#' \code{\link[penalties7]{penalty_prox_spec}} and the compiled sweeps read it,
#' so the kernel names no family and a penalty that describes its operator gets
#' the compiled route without an edit here.
#'
#' \strong{Screening.} Passing from one point of a path to the next, a
#' coordinate can be discarded when the gradient it had at the previous point
#' is below \eqn{2s_k - s_{k-1}}, with \eqn{s} the size of the kink: the
#' sequential strong rule of Tibshirani and others (2012), which assumes the
#' gradient moves at most as fast as the threshold does. That assumption is not
#' a theorem, so the rule can discard a coordinate that belongs in the fit, and
#' what makes the answer exact is the check afterwards: the gradient is read
#' over every column at the point reached, any discarded coordinate whose
#' gradient exceeds the kink is put back, and the fit is repeated. Without the
#' check the route would be wrong now and then rather than slow.
#'
#' \strong{Which update.} The gradient is kept either as a residual, at
#' \eqn{O(n)} a visit, or as itself through
#' \eqn{g_j = (X'Wz)_j - \sum_k (X'WX)_{jk}\beta_k}, at \eqn{O(m)} a change
#' with the Gram columns cached as coordinates come alive. The second wins when
#' \eqn{n} is large next to the number of live coordinates and pays in memory,
#' so the choice is made by size rather than declared.
#'
#' @param obj The full objective.
#' @param beta The current stacked coefficients.
#' @param block One entry of \code{statmod_blocks()$sparse}.
#' @param hyper The hyperparameters.
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param expected Whether the information is the expected one.
#' @param approx How the expected information is approximated.
#' @param maxit The iteration budget.
#' @param tol The stopping tolerance.
#' @param prev_kink The size of the kink at the previous point of a path, or
#'   \code{NULL} to cycle over every coordinate.
#'
#' @return A list shaped like \code{\link{sparse_fit}}'s, or \code{NULL} where
#'   the route does not apply.
#'
#' @references
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. \emph{Journal of
#' Statistical Software} 33(1), 1--22.
#'
#' Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor, J.
#' and Tibshirani, R. J. (2012). Strong rules for discarding predictors in
#' lasso-type problems. \emph{Journal of the Royal Statistical Society, Series
#' B} 74(2), 245--266.
#'
#' @seealso \code{\link{sparse_fit}},
#'   \code{\link[penalties7]{penalty_prox_spec}}
#'
#' @keywords internal
coord_fit <- function(obj, beta, block, hyper, spec, design, expected, approx,
                      maxit = 100L, tol = 1e-8, prev_kink = NULL) {
  p <- block$param
  d <- design[[p]]
  th <- as.list(hyper[[p]][[block$term]])
  cols <- block$cols
  X <- d$X[, cols, drop = FALSE]
  if (!ncol(X)) return(NULL)
  other <- setdiff(seq_len(d$npar), cols)
  step0 <- rep(1, ncol(X))
  if (is.null(penalties7::penalty_prox_spec(block$penalty, th, step0))) {
    return(NULL)
  }

  n <- spec@n_obs
  cur <- beta
  sweeps <- 0L
  for (it in seq_len(maxit)) {
    coef <- obj$split(cur)
    ep <- statmod_eta(spec, design, coef)
    wq <- coord_working(spec, ep, coef, design, p, expected, approx)
    if (is.null(wq)) return(NULL)
    off <- if (length(other))
      as.numeric(d$X[, other, drop = FALSE] %*% coef[[p]][other]) else
      rep(0, n)
    z <- wq$z - off - coord_offset(spec, p, n)
    v <- as.numeric(crossprod(wq$w, X^2))
    if (any(!is.finite(v)) || any(v <= 0)) return(NULL)
    b0 <- coef[[p]][cols]
    s_now <- kink_scale(block$penalty, th)
    keep <- coord_screen(X, wq$w, z, b0, s_now, prev_kink)

    repeat {
      tab <- penalties7::penalty_prox_spec(block$penalty, th, 1 / v[keep])
      if (is.null(tab)) return(NULL)
      out <- coord_descent(X, z, wq$w, b0, tab$cut, tab$slope, tab$icept,
                           as.integer(keep - 1L), 500L, tol,
                           coord_covariance(n, length(keep)))
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
    if (moved < tol) break
  }
  list(par = cur, value = obj$fn(cur), converged = TRUE, iterations = sweeps)
}


#' Which Coordinates a Path Point Has to Visit
#'
#' @description
#' The sequential strong rule: a coordinate whose gradient at the previous
#' point of the path is below \eqn{2s_k - s_{k-1}} is left out of the sweeps.
#'
#' @details
#' The rule rests on the gradient moving no faster than the threshold, which is
#' an assumption rather than a bound, so it screens rather than proves and the
#' caller checks what it discarded. With no previous point there is nothing to
#' screen against and every coordinate is visited.
#'
#' @param X The block's columns.
#' @param w The working weights.
#' @param z The working response.
#' @param beta The coefficients at the previous point.
#' @param s_now The size of the kink here.
#' @param s_prev The size of the kink at the previous point, or \code{NULL}.
#'
#' @return An integer vector of column indices, never empty.
#'
#' @seealso \code{\link{coord_fit}}
#'
#' @keywords internal
coord_screen <- function(X, w, z, beta, s_now, s_prev) {
  p <- ncol(X)
  # With no previous point there is nothing to screen against, and the rule in
  # its global form -- the reference being the kink that empties the block --
  # discards nothing at any smoothing parameter worth fitting at, since
  # 2s - max|g| is negative there. Measured at a single value it cost one
  # crossprod and saved no coordinate, so it is not attempted.
  if (is.null(s_prev) || !is.finite(s_prev) || !is.finite(s_now) ||
      s_now <= 0) {
    return(seq_len(p))
  }
  r <- z - as.numeric(X %*% beta)
  g <- abs(as.numeric(crossprod(X, w * r)))
  keep <- which(g >= 2 * s_now - s_prev | beta != 0)
  if (!length(keep)) keep <- which.max(g)
  keep
}


#' Which Way of Holding the Gradient Is Cheaper
#'
#' @description
#' \code{TRUE} for the covariance form, \code{FALSE} for the residual.
#'
#' @details
#' The covariance form replaces an \eqn{O(n)} read with an \eqn{O(m)} one, and
#' pays for it by building a column of \eqn{X'WX} at \eqn{O(nm)} the first time
#' a coordinate moves off zero. It is therefore worth it only when \eqn{m} is
#' small next to \eqn{n}, and the measurement is unambiguous: at 5000
#' observations with 200 columns screened to 200, taking the covariance form
#' cost 70 milliseconds against 55 for the residual, the Gram columns being
#' dearer than the residual passes they replaced. The threshold is set where
#' the two cross rather than at a rule of thumb.
#'
#' @param n The number of observations.
#' @param m How many coordinates are visited.
#'
#' @return A single logical.
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
#' @param spec A \code{\link{StatmodSpec}}.
#' @param ep The linear predictors and parameters, from
#'   \code{\link{statmod_eta}}.
#' @param coef The coefficients.
#' @param design The design.
#' @param p Which distribution parameter.
#' @param expected Whether the information is the expected one.
#' @param approx How it is approximated.
#'
#' @return A list with \code{w} and \code{z}, or \code{NULL} where the
#'   curvature is not usable.
#'
#' @keywords internal
coord_working <- function(spec, ep, coef, design, p, expected, approx) {
  n <- spec@n_obs
  params <- spec@distrib@params
  a <- match(p, params)
  g <- distributions7::distrib_gradient(spec@distrib, spec@response, ep$theta,
                                        scale = "link")
  s <- spec@weights * rep_len(g[[p]], n)
  Om <- info_blocks(spec, ep$theta, expected, approx)
  h <- Om[, a, a]
  if (any(!is.finite(h)) || any(h <= 0) || any(!is.finite(s))) return(NULL)
  list(w = h, z = rep_len(ep$eta[[p]], n) + s / h)
}


#' The Offset of One Equation
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param p Which distribution parameter.
#' @param n The number of observations.
#'
#' @return A numeric vector of length \code{n}.
#'
#' @keywords internal
coord_offset <- function(spec, p, n) {
  o <- spec@offsets[[p]]
  if (is.null(o) || !length(o)) return(rep(0, n))
  rep_len(as.numeric(o), n)
}
