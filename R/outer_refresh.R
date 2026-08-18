#' @include outer_gradient.R
NULL

# How a block that moves with its coefficients enters the outer derivatives.
#
# With H_{(a,j),(b,k)} = -sum_i w_i l_ab,i X_a[i,j] X_b[i,k] and a block whose
# X is the Jacobian of the term's own contribution, differentiating in beta
# gives three contributions rather than one: the log-likelihood's own third
# derivative with X held, and two from dX/dbeta that are transposes of each
# other under the trace. The gradient consumes them through `u_refresh()`; the
# HESSIAN needs the same derivative in three further places, and the three are
# written here because correcting one of them alone is worse than correcting
# none -- measured, a weakly identified nl goes from 2.02e-01 to 6.35e+02 with
# a standard error of NaN when only dK/dt is corrected.
#
# What is NOT here is d2X/dbeta2, the term's own third derivative. A bilinear f
# is quadratic jointly in its parameters, so that quantity is exactly zero
# there while dX/dbeta is not, and the residual gap on such a model is what
# says the first derivative is the piece that was missing.


#' The Terms Whose Block Moves With Their Coefficients
#'
#' @description
#' One entry per refreshable term, carrying everything its derivative needs:
#' the distribution parameter it sits in, its columns in the stacked
#' coefficient vector, its own coefficients and the term rebuilt at them.
#'
#' @details
#' Resolved once and passed to every consumer, so the gradient's contraction
#' and the Hessian's three corrections cannot disagree about which terms move
#' or where their columns are.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design, already refreshed at \code{coef}.
#' @param coef The coefficients.
#' @param params,npar,offs The block bookkeeping.
#'
#' @return A list, empty where no block moves.
#'
#' @seealso \code{\link{u_refresh}}, \code{\link{contract3_refresh}}
#'
#' @keywords internal
refresh_units <- function(spec, design, coef, params, npar, offs) {
  rf <- attr(design, "refresh")
  st <- attr(design, "state")
  if (is.null(rf) || !length(rf) || is.null(st)) return(list())
  out <- list()
  for (r in rf) {
    p <- r$param
    a <- match(p, params)
    if (is.na(a) || npar[a] == 0L) next
    cols <- design[[p]]$blocks[[r$term]]
    if (!length(cols)) next
    tm <- st$terms[[p]][[r$term]]
    if (is.null(tm)) next
    bt <- coef[[p]][cols]
    out[[length(out) + 1L]] <- list(
      a = a, param = p, term = r$term, cols = cols,
      ra = offs[a] + cols, bt = bt,
      built = modelterms7::term_refresh(tm, bt))
  }
  out
}


#' The Log-Density's Curvature the Refresh Corrections Read
#'
#' @description
#' \eqn{\ell_{ab}} on the link scale, from whichever route the criterion's
#' \eqn{K} was assembled on.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients.
#' @param expected Whether the information is the expected one.
#' @param approx The approximation for the expected information.
#'
#' @return A named list of numeric vectors.
#'
#' @keywords internal
refresh_hessian <- function(spec, design, coef, expected = FALSE,
                            approx = "bartlett") {
  th <- statmod_eta(spec, design, coef)$theta
  if (expected) {
    distributions7::distrib_expected_hessian(spec@distrib, spec@response, th,
                                             scale = "link", approx = approx, threads = spec@threads)
  } else {
    distributions7::distrib_hessian(spec@distrib, spec@response, th,
                                    scale = "link", threads = spec@threads)
  }
}


#' The Weight a Refreshable Block Is Contracted Against
#'
#' @description
#' \eqn{A_a[i,j] = w_i\sum_b c_{ab}(i)\,(X_b M_{ab}^\top)[i,j]}, the matrix a
#' term's own derivative is paired with.
#'
#' @details
#' The same shape serves three quantities and differs only in \eqn{c_{ab}}: the
#' curvature \eqn{\ell_{ab}} for the gradient's contraction, and the third
#' derivative already contracted in a direction for each of the two mixed terms
#' of the Hessian.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param M The matrix the trace is taken against.
#' @param params,npar,offs The block bookkeeping.
#' @param ra The term's rows in the stacked coefficient vector.
#' @param cw One length-\eqn{n} vector per distribution parameter, or
#'   \code{NULL} where that parameter carries no coefficient.
#'
#' @return A numeric matrix, one row per observation and one column per
#'   coefficient of the term.
#'
#' @keywords internal
refresh_amat <- function(spec, design, M, params, npar, offs, ra, cw) {
  n <- spec@n_obs
  A <- matrix(0, n, length(ra))
  for (b in seq_along(params)) {
    if (npar[b] == 0L) next
    rb <- offs[b] + seq_len(npar[b])
    Mab <- as_dense(M[ra, rb, drop = FALSE])
    A <- A + rep_len(cw[[b]], n) * as_dense(design[[params[b]]]$X %*% t(Mab))
  }
  A * spec@weights
}


#' A Refreshable Block's Derivative Along One Direction
#'
#' @description
#' \code{\link[modelterms7]{term_block_deriv}} on one unit, with the shape it
#' returns checked rather than assumed.
#'
#' @param un One entry of \code{\link{refresh_units}}.
#' @param v The direction, as long as the term's coefficients.
#' @param n The number of observations.
#'
#' @return A numeric matrix, \code{n} by the term's coefficient count.
#'
#' @keywords internal
refresh_dblock <- function(un, v, n) {
  D <- as_dense(modelterms7::term_block_deriv(un$built, coef = un$bt, v = v))
  if (!is.matrix(D) || nrow(D) != n || ncol(D) != length(un$cols)) {
    stop(sprintf(paste("term_block_deriv() returned a %s block where %d by %d",
                       "was wanted\n  in '%s'."),
                 paste(dim(as.matrix(D)), collapse = " by "), n,
                 length(un$cols), un$term), call. = FALSE)
  }
  D
}


#' The Third Derivative Already Contracted in One Direction
#'
#' @description
#' \eqn{\sum_k \ell_{abk}(X_k v_k)_i}, the per-observation weight
#' \code{\link{contract3}} builds for the \eqn{(a,b)} block.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param d3 The third derivatives on the link scale.
#' @param params,npar The parameter names and block sizes.
#' @param a,b The block indices.
#' @param tv The predictors of the direction.
#'
#' @return A numeric vector as long as the sample.
#'
#' @keywords internal
d3_direction <- function(spec, d3, params, npar, a, b, tv) {
  n <- spec@n_obs
  keys <- names(d3)
  w <- numeric(n)
  for (k in seq_along(params)) {
    if (npar[k] == 0L) next
    w <- w + rep_len(d3[[d3_key(params, a, b, k, keys)]], n) * tv[[k]]
  }
  w
}


#' How a Moving Block Enters the Contracted Third Derivative
#'
#' @description
#' The part of \eqn{T[v] = (\partial K/\partial\beta)\cdot v} that
#' \code{\link{contract3}} does not compute, as a matrix over the stacked
#' coefficients.
#'
#' @details
#' Differentiating \eqn{H} in \eqn{\beta} and contracting with \eqn{v} leaves,
#' beside the third derivative \code{\link{contract3}} carries, two terms in
#' \eqn{\partial X/\partial\beta}: with \eqn{D_a = (\partial X_a/\partial\beta)v},
#' \deqn{R[(a,j),(b,k)] = -\sum_i w_i\,\ell_{ab,i}\,D_a[i,j]\,X_b[i,k],}
#' and the other is its transpose, so the correction is \eqn{R + R^\top}. The
#' derivative is asked of the TERM through
#' \code{\link[modelterms7]{term_block_deriv}} and never differenced here, for
#' the reason \code{\link{u_refresh}} records: a break-point column is a step
#' function in its break-point.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param params,npar,offs,total The block bookkeeping.
#' @param dir One entry per unit, from \code{\link{refresh_direction}}.
#' @param Hl The link-scale curvature, from \code{\link{refresh_hessian}}.
#' @param units The refreshable terms, from \code{\link{refresh_units}}.
#'
#' @return A square matrix, whose transpose completes the correction.
#'
#' @seealso \code{\link{contract3}}, \code{\link{u_refresh}}
#'
#' @keywords internal
contract3_refresh <- function(spec, design, params, npar, offs, total, dir, Hl,
                              units) {
  R <- zero_information(design, total)
  if (!length(units)) return(R)
  n <- spec@n_obs
  w <- spec@weights
  for (s in seq_along(units)) {
    un <- units[[s]]
    D <- dir[[s]]$D
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      rb <- offs[b] + seq_len(npar[b])
      lab <- rep_len(Hl[[hess_key(params, un$a, b)]], n)
      R[un$ra, rb] <- R[un$ra, rb] -
        as_dense(crossprod(D * (w * lab), design[[params[b]]]$X))
    }
  }
  R
}


#' What One Direction Costs a Moving Block
#'
#' @description
#' The two quantities every consumer of \eqn{\partial X/\partial\beta} reads
#' for a given direction: the block's own derivative along it, and the third
#' derivative already contracted in it and carried onto the term's columns.
#'
#' @details
#' Both depend on ONE direction and not on a pair, so they are built once per
#' hyperparameter and combined in the Hessian's pair loop: computing them
#' inside that loop would repeat an \eqn{O(np^2)} product for every pair.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param M The matrix the trace is taken against.
#' @param params,npar,offs The block bookkeeping.
#' @param d3 The third derivatives on the link scale.
#' @param tv The predictors of the direction.
#' @param v The direction, over the stacked coefficients.
#' @param units The refreshable terms, from \code{\link{refresh_units}}.
#'
#' @return A list, one entry per unit, with \code{D} and \code{A}.
#'
#' @seealso \code{\link{contract3_refresh}}, \code{\link{trace_refresh4}}
#'
#' @keywords internal
refresh_direction <- function(spec, design, M, params, npar, offs, d3, tv, v,
                              units) {
  lapply(units, function(un) {
    cw <- vector("list", length(params))
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      cw[[b]] <- d3_direction(spec, d3, params, npar, un$a, b, tv)
    }
    list(D = refresh_dblock(un, v[un$ra], spec@n_obs),
         A = refresh_amat(spec, design, M, params, npar, offs, un$ra, cw))
  })
}


#' How a Moving Block Enters the Twice-Contracted Fourth Derivative
#'
#' @description
#' The part of \eqn{\mathrm{tr}(M\,U[v,u])} that \code{\link{trace_design_form}}
#' does not compute, where \eqn{U} is the second derivative of \eqn{K} in the
#' coefficients contracted in two directions.
#'
#' @details
#' Differentiating the three contributions of \eqn{\partial K/\partial\beta}
#' once more leaves three further terms that read \eqn{\partial X/\partial\beta}
#' and one that would read \eqn{\partial^2X/\partial\beta^2}. Writing
#' \eqn{D_a = (\partial X_a/\partial\beta)v} and
#' \eqn{E_a = (\partial X_a/\partial\beta)u},
#' \deqn{N = N_1 + N_2 + N_3, \qquad
#'   N_1[(a,j),(b,k)] = -\sum_i w_i\Big(\sum_m \ell_{abm}(X_mv_m)_i\Big)
#'     E_a[i,j]X_b[i,k],}
#' \eqn{N_2} the same with the two directions exchanged, and
#' \eqn{N_3[(a,j),(b,k)] = -\sum_i w_i\ell_{ab,i}D_a[i,j]E_b[i,k]}, which is
#' supported where BOTH blocks move. The correction is
#' \eqn{\mathrm{tr}(M(N + N^\top)) = 2\sum M\odot N}, and the first two pieces
#' are read off \code{\link{refresh_amat}} with the third derivative in place of
#' the curvature, so neither \eqn{N} nor any contraction is assembled.
#'
#' \strong{The omitted term} carries \eqn{\partial^2X/\partial\beta^2}, the
#' term's own third derivative, which no contract supplies. Measured on a
#' bilinear \eqn{f}, where it is exactly zero, the corrected Hessian reaches the
#' reference's own floor.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param M The matrix the trace is taken against.
#' @param params The parameter names.
#' @param Hl The link-scale curvature, from \code{\link{refresh_hessian}}.
#' @param dv,du The two directions, from \code{\link{refresh_direction}}.
#' @param units The refreshable terms, from \code{\link{refresh_units}}.
#'
#' @return A single number.
#'
#' @seealso \code{\link{trace_design_form}}, \code{\link{contract3_refresh}}
#'
#' @keywords internal
trace_refresh4 <- function(spec, M, params, Hl, dv, du, units) {
  if (!length(units)) return(0)
  n <- spec@n_obs
  w <- spec@weights
  out <- 0
  for (s in seq_along(units)) {
    out <- out - sum(du[[s]]$D * dv[[s]]$A) - sum(dv[[s]]$D * du[[s]]$A)
    for (t2 in seq_along(units)) {
      Mab <- as_dense(M[units[[s]]$ra, units[[t2]]$ra, drop = FALSE])
      lab <- rep_len(Hl[[hess_key(params, units[[s]]$a, units[[t2]]$a)]], n)
      out <- out - sum(w * lab * rowSums((dv[[s]]$D %*% Mab) * du[[t2]]$D))
    }
  }
  2 * out
}
