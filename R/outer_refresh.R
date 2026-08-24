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
# The HESSIAN's twice-contracted term needs two further pieces, and they are
# written together for the same reason: measured against a mixed second
# difference of H, adding either one alone leaves a curved nl WORSE or barely
# better (2.62e-02 today, 2.88e-02 with the first alone, 2.65e-03 with the
# second alone, 5.19e-08 with both).
#
# The first is d2X/dbeta2, the term's own third derivative, which
# term_block_deriv2() supplies. The second reads only dX/dbeta and had been
# missing since these corrections were written: the third-derivative
# contraction carries the DIRECTION'S PREDICTOR (X_m v_m), and where X_m is a
# Jacobian that predictor moves with beta too. On a bilinear f, where
# d2X/dbeta2 is exactly zero, it is the whole of the gap -- 2.32e-02 against
# an arbitrary trace matrix, 2.19e-08 once it is in.


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
#' @param spec A [StatmodSpec()].
#' @param design The design, already refreshed at `coef`.
#' @param coef The coefficients.
#' @param params,npar,offs The block bookkeeping.
#'
#' @return A list, empty where no block moves.
#'
#' @seealso [u_refresh()], [contract3_refresh()]
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
#' @param spec A [StatmodSpec()].
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
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param M The matrix the trace is taken against.
#' @param params,npar,offs The block bookkeeping.
#' @param ra The term's rows in the stacked coefficient vector.
#' @param cw One length-\eqn{n} vector per distribution parameter, or
#'   `NULL` where that parameter carries no coefficient.
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
#' [modelterms7::term_block_deriv()] on one unit, with the shape it
#' returns checked rather than assumed.
#'
#' @param un One entry of [refresh_units()].
#' @param v The direction, as long as the term's coefficients.
#' @param n The number of observations.
#'
#' @return A numeric matrix, `n` by the term's coefficient count.
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


#' A Refreshable Block's Second Derivative Along Two Directions
#'
#' @description
#' [modelterms7::term_block_deriv2()] on one unit, with the shape it
#' returns checked rather than assumed.
#'
#' @details
#' A term that has not written the method inherits the base one and gets
#' zeros, which is exact for a design that does not move and a deliberate
#' refusal for a block that is a working linearization rather than a Jacobian.
#' Nothing is differenced here in either case.
#'
#' @param un One entry of [refresh_units()].
#' @param v,u The two directions, each as long as the term's coefficients.
#' @param n The number of observations.
#'
#' @return A numeric matrix, `n` by the term's coefficient count.
#'
#' @seealso [refresh_dblock()], [trace_refresh4()]
#'
#' @keywords internal
refresh_dblock2 <- function(un, v, u, n) {
  D <- as_dense(modelterms7::term_block_deriv2(un$built, coef = un$bt, v = v,
                                               u = u))
  if (!is.matrix(D) || nrow(D) != n || ncol(D) != length(un$cols)) {
    stop(sprintf(paste("term_block_deriv2() returned a %s block where %d by %d",
                       "was wanted\n  in '%s'."),
                 paste(dim(as.matrix(D)), collapse = " by "), n,
                 length(un$cols), un$term), call. = FALSE)
  }
  D
}


#' The Weight the Second Derivative of a Block Is Paired With
#'
#' @description
#' [refresh_amat()] with the curvature, one matrix per unit.
#'
#' @details
#' It is the same weight [u_refresh()] builds for the gradient, and
#' it depends on neither direction, so it is built once per Hessian rather
#' than once per pair of hyperparameters.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param M The matrix the trace is taken against.
#' @param params,npar,offs The block bookkeeping.
#' @param units The refreshable terms, from [refresh_units()].
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#'
#' @return A list of matrices, one per unit.
#'
#' @seealso [trace_refresh4()], [u_refresh()]
#'
#' @keywords internal
refresh_curv_amat <- function(spec, design, M, params, npar, offs, units, Hl) {
  lapply(units, function(un) {
    cw <- vector("list", length(params))
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      cw[[b]] <- Hl[[hess_key(params, un$a, b)]]
    }
    refresh_amat(spec, design, M, params, npar, offs, un$ra, cw)
  })
}


#' The Third Derivative Already Contracted in One Direction
#'
#' @description
#' \eqn{\sum_k \ell_{abk}(X_k v_k)_i}, the per-observation weight
#' [contract3()] builds for the \eqn{(a,b)} block.
#'
#' @param spec A [StatmodSpec()].
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
#' [contract3()] does not compute, as a matrix over the stacked
#' coefficients.
#'
#' @details
#' Differentiating \eqn{H} in \eqn{\beta} and contracting with \eqn{v} leaves,
#' beside the third derivative [contract3()] carries, two terms in
#' \eqn{\partial X/\partial\beta}: with \eqn{D_a = (\partial X_a/\partial\beta)v},
#' \deqn{R[(a,j),(b,k)] = -\sum_i w_i\,\ell_{ab,i}\,D_a[i,j]\,X_b[i,k],}
#' and the other is its transpose, so the correction is \eqn{R + R^\top}. The
#' derivative is asked of the TERM through
#' [modelterms7::term_block_deriv()] and never differenced here, for
#' the reason [u_refresh()] records: a break-point column is a step
#' function in its break-point.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param params,npar,offs,total The block bookkeeping.
#' @param dir One entry per unit, from [refresh_direction()].
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#' @param units The refreshable terms, from [refresh_units()].
#'
#' @return A square matrix, whose transpose completes the correction.
#'
#' @seealso [contract3()], [u_refresh()]
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
#' Both depend on one direction alone, never on a pair, so they are built
#' once per
#' hyperparameter and combined in the Hessian's pair loop: computing them
#' inside that loop would repeat an \eqn{O(np^2)} product for every pair.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param M The matrix the trace is taken against.
#' @param params,npar,offs The block bookkeeping.
#' @param d3 The third derivatives on the link scale.
#' @param tv The predictors of the direction.
#' @param v The direction, over the stacked coefficients.
#' @param units The refreshable terms, from [refresh_units()].
#'
#' @return A list, one entry per unit, with `D` and `A`.
#'
#' @seealso [contract3_refresh()], [trace_refresh4()]
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
#' The part of \eqn{\mathrm{tr}(M\,U[v,u])} that [trace_design_form()]
#' does not compute, where \eqn{U} is the second derivative of \eqn{K} in the
#' coefficients contracted in two directions.
#'
#' @details
#' Differentiating the three contributions of \eqn{\partial K/\partial\beta}
#' once more leaves FIVE terms. Writing \eqn{D_a = (\partial X_a/\partial\beta)v}
#' and \eqn{E_a = (\partial X_a/\partial\beta)u}, three of them read the first
#' derivative alone,
#' \deqn{N_1[(a,j),(b,k)] = -\sum_i w_i\Big(\sum_m \ell_{abm}(X_mv_m)_i\Big)
#'     E_a[i,j]X_b[i,k],}
#' \eqn{N_2} the same with the two directions exchanged, and
#' \eqn{N_3[(a,j),(b,k)] = -\sum_i w_i\ell_{ab,i}D_a[i,j]E_b[i,k]}, which is
#' supported where BOTH blocks move. A fourth reads the SECOND,
#' \eqn{N_4[(a,j),(b,k)] = -\sum_i w_i\ell_{ab,i}F_a[i,j]X_b[i,k]} with
#' \eqn{F_a = (\partial^2X_a/\partial\beta^2)[v,u]} from
#' [modelterms7::term_block_deriv2()]. Those four enter as
#' \eqn{\mathrm{tr}(M(N + N^\top)) = 2\sum M\odot N}, and each is read off
#' [refresh_amat()], with the third derivative for the first two and the
#' curvature for the fourth, so neither \eqn{N} nor any contraction is
#' assembled.
#'
#' The fifth is of another shape and is added separately. The
#' third-derivative contraction carries the DIRECTION'S PREDICTOR
#' \eqn{(X_mv_m)_i}, which where \eqn{X_m} is a Jacobian moves with
#' \eqn{\beta} as everything else does, and differentiating it leaves
#' \deqn{-\sum_i w_i\Big(\sum_m \ell_{abm}\,\dot v_m(i)\Big)X_a[i,j]X_b[i,k],
#'   \qquad \dot v_m = (\partial X_m/\partial\beta[u])\,v_m.}
#' That is the shape [trace_design_form()] already computes, so the
#' piece is one further call of it with \eqn{\dot v} in place of the
#' predictors. \eqn{\dot v} is the predictor's own second derivative
#' contracted in both directions and is therefore symmetric in them, which is
#' a free check on the assembly.
#'
#' **Measured**, against a mixed second difference of \eqn{H} that shares
#' no arithmetic with any of this: on a bilinear \eqn{f}, where the fourth
#' term is exactly zero, the fifth is the whole of the gap and takes
#' \eqn{2.32\times10^{-2}} to \eqn{2.19\times10^{-8}}; on a curved one neither
#' alone suffices: \eqn{2.62\times10^{-2}} today, \eqn{2.88\times10^{-2}}
#' with the fifth alone, \eqn{2.65\times10^{-3}} with the fourth alone,
#' \eqn{5.19\times10^{-8}} with both.
#'
#' @param spec A [StatmodSpec()].
#' @param M The matrix the trace is taken against.
#' @param params,npar The parameter names and block sizes.
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#' @param dv,du The two directions, from [refresh_direction()].
#' @param units The refreshable terms, from [refresh_units()].
#' @param G The leverage diagonal, from [block_leverage()].
#' @param d3 The third derivatives on the link scale.
#' @param v The first direction over the stacked coefficients.
#' @param f2 The block's second derivative in the two directions, one matrix
#'   per unit, from [refresh_dblock2()].
#' @param acurv The curvature weights, from [refresh_curv_amat()].
#'
#' @return A single number.
#'
#' @seealso [trace_design_form()], [contract3_refresh()],
#'   [refresh_dblock2()]
#'
#' @keywords internal
trace_refresh4 <- function(spec, M, params, npar, Hl, dv, du, units, G, d3,
                           v, f2, acurv) {
  if (!length(units)) return(0)
  n <- spec@n_obs
  w <- spec@weights
  out <- 0
  vdot <- lapply(seq_along(params), function(m) numeric(n))
  for (s in seq_along(units)) {
    un <- units[[s]]
    out <- out - sum(du[[s]]$D * dv[[s]]$A) - sum(dv[[s]]$D * du[[s]]$A)
    for (t2 in seq_along(units)) {
      Mab <- as_dense(M[un$ra, units[[t2]]$ra, drop = FALSE])
      lab <- rep_len(Hl[[hess_key(params, un$a, units[[t2]]$a)]], n)
      out <- out - sum(w * lab * rowSums((dv[[s]]$D %*% Mab) * du[[t2]]$D))
    }
    out <- out - sum(f2[[s]] * acurv[[s]])
    vdot[[un$a]] <- vdot[[un$a]] + as.numeric(du[[s]]$D %*% v[un$ra])
  }
  2 * out + trace_design_form(spec, G, d3, params, npar, vdot)
}


#' How a Moving Block Enters the Mode's Second Movement
#'
#' @description
#' The part of \eqn{\partial^3L/\partial\beta^3[v,u]} that
#' [contract3()] and [contract3_refresh()] do not carry,
#' as a vector over the stacked coefficients.
#'
#' @details
#' \eqn{b_{ml}} solves \eqn{J b_{ml} = -[(S_l + T[b_l])b_m + S_m b_l + c_{ml}]}
#' with \eqn{T} the THIRD derivative of the penalized objective in \eqn{\beta}.
#' Where a block moves, the objective's second derivative is not \eqn{K} but
#' \eqn{K + D}, with \eqn{D} the term [mode_curvature()] builds, so
#' its third derivative is not \eqn{\partial K/\partial\beta} either. What is
#' missing is \eqn{\partial D/\partial\beta}, and differentiating
#' \eqn{D_{cd} = -\sum_i w_i\ell_a\,\Xi_a[i,c,d]} once and contracting gives two
#' pieces,
#' \deqn{-\sum_i w_i\Big(\sum_k \ell_{ak}(X_kv_k)_i\Big)\,
#'   \big[(\partial X_a/\partial\beta)u\big][i,c]
#'   \;-\;\sum_i w_i\,\ell_a(i)\,F_a[i,c],}
#' the first reading the block's first derivative and the second its SECOND,
#' \eqn{F_a = (\partial^2X_a/\partial\beta^2)[v,u]}. Both are one
#' `crossprod` against a per-observation weight, so neither the
#' third-derivative array of the predictor nor any contraction of it is formed.
#'
#' **Measured** against the mode refitted at four hyperparameter values
#' and differenced twice, on `nl(a ~ 0 + ridge(~grp))`: \eqn{b_m} is right
#' to 3.8e-08 without this and \eqn{b_{ml}} is wrong by 7.6 to 9.0 per cent at
#' every step tried, which is a systematic error and no reference's noise.
#' With it
#' the gap falls to 1.8e-05 at the reference's best step, which is inside the
#' spread between the reference's own consecutive steps.
#'
#' @param spec A [StatmodSpec()].
#' @param params,npar,total The block bookkeeping.
#' @param units The refreshable terms, from [refresh_units()].
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#' @param gl The link-scale score.
#' @param tv The predictors of the first direction.
#' @param du The second direction, from [refresh_direction()].
#' @param f2 The block's second derivative, one matrix per unit.
#'
#' @return A numeric vector as long as the stacked coefficients.
#'
#' @seealso [mode_curvature()], [trace_refresh4()]
#'
#' @keywords internal
refresh_mode_third <- function(spec, params, npar, units, Hl, gl, tv, du, f2,
                               total) {
  out <- numeric(total)
  if (!length(units)) return(out)
  n <- spec@n_obs
  w <- spec@weights
  for (s in seq_along(units)) {
    un <- units[[s]]
    hk <- numeric(n)
    for (k in seq_along(params)) {
      if (npar[k] == 0L) next
      hk <- hk + rep_len(Hl[[hess_key(params, un$a, k)]], n) * tv[[k]]
    }
    out[un$ra] <- out[un$ra] -
      as.numeric(crossprod(du[[s]]$D, w * hk)) -
      as.numeric(crossprod(f2[[s]], w * rep_len(gl[[un$param]], n)))
  }
  out
}
