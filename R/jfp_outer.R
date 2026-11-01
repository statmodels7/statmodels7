#' @include jfp_engine.R jfp_move.R
NULL

# The third and fourth derivatives of log|I_tilde| in the joint vector, in a
# model without a filter. Every quantity is written in
# the whitened coordinates of each block, I_tilde = U'U and Z_q = V_q U^-1,
# where the derivatives of the block's matrix are
#
#   A_{k x}  = sum_t w_t sum_qr d^{1+|x|} I_qr / d eta_c d x  X_{c,t,k} Z_q'Z_r
#   A_{kl x} = sum_t w_t sum_qr d^{2+|x|} I_qr / d eta_c d eta_d d x
#              X_{c,t,k} X_{d,t,l} Z_q'Z_r
#   A_x      = sum_t w_t sum_qr d^{|x|} I_qr / d x  Z_q'Z_r
#
# for a set x of directions in u, and I_qr = -E_qr is the family's expected
# information on the link scale. Differentiating log|I_tilde| gives
#
#   D3(v)_kl   = tr A_klv - tr(A_v A_kl) - tr(A_kv A_l) - tr(A_k A_lv)
#                + tr(A_v A_k A_l) + tr(A_k A_v A_l)
#
# and D4(v, w) its derivative along w, in which every A_x gains the direction
# w and every inverse of I_tilde inserts -A_w. A term with k and l inside one
# A is a weighted cross product of the designs, with the sparsity of the
# information; a term with k and l in two different A's is B_x K B_y', B_x
# holding vec(A_{k x}) in row k, with K a Kronecker product of p_tilde x
# p_tilde matrices, by tr(X P Y Q) = vec(X)'(Q' (x) P) vec(Y) for symmetric X.
# Where a predictor is not linear in u, or the block's rows move, each A
# gains the shares that R/jfp_move.R computes, and a term with k and l inside
# one A gains dense entries over the coordinates of the moving terms.


#' The Third and Fourth Derivatives of the Jeffreys Penalty
#'
#' @description
#' `jfp_deriv3_dir()` returns the derivative of the Hessian of
#' \eqn{\sum_b a_b\log|\tilde I_b|} along a direction \eqn{v} of the joint
#' vector, and `jfp_deriv4_dir()` its second derivative along two directions
#' \eqn{v} and \eqn{w}, as matrices over the joint vector. `jfp_deriv3()` and
#' `jfp_deriv4()` return their traces against a symmetric matrix \eqn{M},
#' without assembling them.
#'
#' @details
#' The derivatives are those of the determinant with the derivatives of the
#' family's expected information to the fourth order, from
#' [distributions7::distrib_dexpected_hessian()] and its siblings of orders
#' two to four. A term with both indices inside one derivative of the block's
#' matrix is a weighted cross product of the designs and has the sparsity of
#' the information; the others are products of low rank through matrices of
#' \eqn{\tilde p^2} columns. Where a term recomputes its block at the
#' coefficients ([modelterms7::nl()], a smoothed [modelterms7::seg()]), the
#' derivatives of the predictors to the fourth order and, for a block in the
#' subformula of such a term, the movement of the block's rows enter as well
#' ([jfp_mv_ctx()]); the derivatives of the term's block come from
#' [modelterms7::term_block_deriv()] and its siblings up to
#' [modelterms7::term_block_deriv4()]. In a model with a filter the
#' predictors are those of the joint vector ([jfp_joint_rows()]): their
#' derivatives come from the filter's kernels along the directions, and for a
#' block in the subformula of a filter's parameter the rows' derivatives are
#' read along the whitened directions of the block, one kernel call per
#' direction ([modelterms7::term_third()], [modelterms7::term_fourth()],
#' [modelterms7::term_fifth()]).
#'
#' The outer criterion reads them: the derivative of its penalized matrix
#' along the movement of the mode carries \eqn{-a\,D_3(v)}, and its second
#' derivative \eqn{-a\,D_4(v, w)} ([statmod_penalty_dbeta()],
#' [statmod_penalty_second()]).
#'
#' @param st A state, from [jfp_state()].
#' @param v,w Directions over the joint vector.
#' @param M A symmetric matrix over the joint vector.
#'
#' @return `jfp_deriv3_dir()` and `jfp_deriv4_dir()` a base matrix over the
#'   joint vector; `jfp_deriv3()` and `jfp_deriv4()` a single number.
#'
#' @seealso [jfp_hessian()]
#'
#' @keywords internal
jfp_deriv3_dir <- function(st, v) {
  jfp_terms_matrix(st, jfp_terms3(st, v))
}

#' @rdname jfp_deriv3_dir
#' @keywords internal
jfp_deriv3 <- function(st, M, v) {
  jfp_terms_trace(st, jfp_terms3(st, v), M)
}

#' @rdname jfp_deriv3_dir
#' @keywords internal
jfp_deriv4_dir <- function(st, v, w) {
  jfp_terms_matrix(st, jfp_terms4(st, v, w))
}

#' @rdname jfp_deriv3_dir
#' @keywords internal
jfp_deriv4 <- function(st, M, v, w) {
  jfp_terms_trace(st, jfp_terms4(st, v, w), M)
}


#' The Terms of the Third and Fourth Derivatives
#'
#' @description
#' Returns the terms of \eqn{D_3(v)} (`jfp_terms3()`) and \eqn{D_4(v, w)}
#' (`jfp_terms4()`), summed over the blocks with their powers, in the form
#' that [jfp_terms_matrix()] assembles and [jfp_terms_trace()] traces.
#'
#' @param st A state, from [jfp_state()].
#' @param v,w Directions over the joint vector.
#'
#' @return A list with `sparse`, a list of per-observation weights keyed as
#'   [statmod_information_at()] reads its `H` (the matrix is
#'   \eqn{\sum_{c,d}X_c^\top\mathrm{diag}(w\,W_{cd})X_d}), `lr`, a list of
#'   triples `X`, `K`, `Y` (the matrix is \eqn{XKY^\top}), and `dense`, a
#'   list of entries `idx`, `X` (the matrix places `X` in the rows `idx` and
#'   its transpose in the columns `idx`).
#'
#' @keywords internal
jfp_terms3 <- function(st, v) {
  jfp_check_outer(st)
  ds <- list(v = jfp_dir(st, v))
  out <- list(sparse = list(), lr = list(), dense = list())
  e <- character(0)
  for (bs in st$blocks) {
    cx <- jfp_mv_ctx(st, bs, ds)
    a <- bs$block$power
    p <- cx$p
    Ip <- diag(p)
    Av <- jfp_mv_white(cx, "v")
    B <- jfp_mv_B(cx, e)
    Bv <- jfp_mv_B(cx, "v")
    out <- jfp_add_pair(out, jfp_mv_pair(cx, Ip, "v"), a)
    out <- jfp_add_pair(out, jfp_mv_pair(cx, Av, e), -a)
    out$lr <- c(out$lr, list(
      list(X = Bv, K = -a * diag(p * p), Y = B),
      list(X = B, K = -a * diag(p * p), Y = Bv),
      list(X = B, K = a * (kronecker(Av, Ip) + kronecker(Ip, Av)), Y = B)))
  }
  out
}

#' @rdname jfp_terms3
#' @keywords internal
jfp_terms4 <- function(st, v, w) {
  jfp_check_outer(st)
  ds <- list(v = jfp_dir(st, v), w = jfp_dir(st, w))
  out <- list(sparse = list(), lr = list(), dense = list())
  e <- character(0)
  vw <- c("v", "w")
  for (bs in st$blocks) {
    cx <- jfp_mv_ctx(st, bs, ds)
    a <- bs$block$power
    p <- cx$p
    Ip <- diag(p)
    I2 <- diag(p * p)
    Av <- jfp_mv_white(cx, "v")
    Aw <- jfp_mv_white(cx, "w")
    Avw <- jfp_mv_white(cx, vw)
    B <- jfp_mv_B(cx, e)
    Bv <- jfp_mv_B(cx, "v")
    Bw <- jfp_mv_B(cx, "w")
    Bvw <- jfp_mv_B(cx, vw)
    # the two terms without a direction share one trace, being linear in Y
    out <- jfp_add_pair(out, jfp_mv_pair(cx, Ip, vw), a)
    out <- jfp_add_pair(out, jfp_mv_pair(cx, Aw, "v"), -a)
    out <- jfp_add_pair(out, jfp_mv_pair(cx, Av, "w"), -a)
    out <- jfp_add_pair(out, jfp_mv_pair(cx, Aw %*% Av + Av %*% Aw - Avw, e),
                        a)
    kr <- kronecker
    lr <- list(
      # minus the derivative of tr(A_kv A_l)
      list(Bvw, -I2, B), list(Bv, -I2, Bw),
      list(Bv, kr(Aw, Ip), B), list(Bv, kr(Ip, Aw), B),
      # minus the derivative of tr(A_k A_lv)
      list(Bw, -I2, Bv), list(B, -I2, Bvw),
      list(B, kr(Aw, Ip), Bv), list(B, kr(Ip, Aw), Bv),
      # the derivative of tr(A_v A_k A_l)
      list(B, kr(Avw, Ip), B), list(Bw, kr(Av, Ip), B),
      list(B, kr(Av, Ip), Bw), list(B, -kr(Av %*% Aw, Ip), B),
      list(B, -kr(Aw %*% Av, Ip), B), list(B, -kr(Av, Aw), B),
      # the derivative of tr(A_k A_v A_l)
      list(Bw, kr(Ip, Av), B), list(B, kr(Ip, Avw), B),
      list(B, kr(Ip, Av), Bw), list(B, -kr(Aw, Av), B),
      list(B, -kr(Ip, Aw %*% Av), B), list(B, -kr(Ip, Av %*% Aw), B))
    out$lr <- c(out$lr, lapply(lr, function(z)
      list(X = z[[1L]], K = a * z[[2L]], Y = z[[3L]])))
  }
  out
}


#' Assemble or Trace a List of Derivative Terms
#'
#' @param st A state, from [jfp_state()].
#' @param terms What [jfp_terms3()] or [jfp_terms4()] returns.
#' @param M A symmetric matrix over the joint vector.
#'
#' @return `jfp_terms_matrix()` a symmetric base matrix over the joint
#'   vector; `jfp_terms_trace()` the trace of that matrix against `M`, a
#'   single number.
#'
#' @keywords internal
jfp_terms_matrix <- function(st, terms) {
  if (!is.null(st$joint)) {
    out <- jfp_joint_sparse(st, terms$sparse)
  } else out <- as_dense(statmod_information_at(
    st$spec, st$coef, st$design, expected = TRUE, emethod = st$emethod,
    H = lapply(terms$sparse, function(x) -x)))
  for (z in terms$lr) out <- out + z$X %*% z$K %*% t(z$Y)
  for (z in terms$dense) {
    out[z$idx, ] <- out[z$idx, ] + z$X
    out[, z$idx] <- out[, z$idx] + t(z$X)
  }
  (out + t(out)) / 2
}

#' @rdname jfp_terms_matrix
#' @keywords internal
jfp_terms_trace <- function(st, terms, M) {
  M <- as_dense(M)
  params <- st$params
  np <- length(params)
  n <- st$spec@n_obs
  w <- st$spec@weights
  if (!is.null(st$joint)) {
    # the weighted cross products of the joint rows, traced densely
    s <- 0
    V <- st$joint$V
    for (cc in seq_len(np)) {
      VM <- V[[cc]] %*% M
      for (dd in seq_len(np)) {
        Wcd <- terms$sparse[[hess_key(params, cc, dd)]]
        if (is.null(Wcd)) next
        s <- s + sum(w * rep_len(Wcd, n) * rowSums(VM * V[[dd]]))
      }
    }
  } else {
  # the outer derivatives trace against one M several times at a point
  if (is.null(st$lev_M) || !identical(st$lev_M, M)) {
    st$lev_M <- M
    st$lev_G <- block_leverage(st$design, M, params, st$npar, st$offs,
                               st$spec@threads)
  }
  G <- st$lev_G
  s <- 0
  for (cc in seq_len(np)) {
    if (st$npar[cc] == 0L) next
    for (dd in seq_len(np)) {
      if (st$npar[dd] == 0L) next
      Wcd <- terms$sparse[[hess_key(params, cc, dd)]]
      if (is.null(Wcd)) next
      s <- s + sum(w * rep_len(Wcd, n) * G[[cc]][[dd]])
    }
  }
  }
  for (z in terms$lr) {
    s <- s + sum((crossprod(z$X, M %*% z$Y)) * z$K)
  }
  for (z in terms$dense) s <- s + 2 * sum(M[z$idx, , drop = FALSE] * z$X)
  s
}


#' The Movement of the Predictors and the Information Along Directions
#'
#' @description
#' `jfp_eta_dir()` returns the movement of each equation's predictor along a
#' direction, the Jacobian times the direction; `jfp_dI()` a derivative of
#' the family's expected information, with some predictor indices fixed and
#' the others contracted against per-equation movements.
#'
#' @param st A state, from [jfp_state()].
#' @param v A direction over the joint vector.
#' @param q,r The pair of the information.
#' @param fixed The equations of the derivatives that are not contracted.
#' @param dirs A list of movements, each a list of numeric vectors, one per
#'   distribution parameter.
#'
#' @return `jfp_eta_dir()` a list of numeric vectors, one per distribution
#'   parameter; `jfp_dI()` a numeric vector over the observations.
#'
#' @keywords internal
jfp_eta_dir <- function(st, v) {
  if (!is.null(st$joint)) {
    return(lapply(st$joint$V, function(X) as.numeric(X %*% v)))
  }
  lapply(seq_along(st$params), function(cc) {
    if (st$npar[cc] == 0L) return(numeric(st$spec@n_obs))
    as.numeric(st$design[[st$params[cc]]]$X %*%
                 v[st$offs[cc] + seq_len(st$npar[cc])])
  })
}

#' @rdname jfp_eta_dir
#' @keywords internal
jfp_dI <- function(st, q, r, fixed, dirs) {
  params <- st$params
  np <- length(params)
  n <- st$spec@n_obs
  k <- length(fixed) + length(dirs)
  if (k == 0L) return(-rep_len(st$E[[hess_key(params, q, r)]], n))
  D <- jfp_darray(st, k)
  key <- switch(k,
    function(i) distributions7::dexpected_key(params, q, r, i[1L]),
    function(i) distributions7::d2expected_key(params, q, r, i[1L], i[2L]),
    function(i) distributions7::d3expected_key(params, q, r, i[1L], i[2L],
                                               i[3L]),
    function(i) distributions7::d4expected_key(params, q, r, i[1L], i[2L],
                                               i[3L], i[4L]))
  out <- numeric(n)
  combos <- if (length(dirs)) {
    as.matrix(expand.grid(rep(list(seq_len(np)), length(dirs))))
  } else matrix(integer(0), 1L, 0L)
  for (j in seq_len(nrow(combos))) {
    e <- combos[j, ]
    m <- rep(1, n)
    for (i in seq_along(dirs)) m <- m * dirs[[i]][[e[i]]]
    if (!any(m != 0)) next
    out <- out - rep_len(D[[key(c(fixed, e))]], n) * m
  }
  out
}


#' Small Helpers of the Third and Fourth Derivatives
#'
#' @description
#' `jfp_live()` returns the equations in which a block has rows;
#' `jfp_add_weights()` sums weight lists with coefficients, and
#' `jfp_add_pair()` adds a trace term of [jfp_mv_pair()] with a coefficient;
#' `jfp_darray()` returns the family's derivative of the expected
#' information of order `k`, kept in the state; `jfp_check_outer()` signals
#' an error where a block's information is not positive definite.
#'
#' @param bs A block's entry of the state.
#' @param acc,lst,coef An accumulated weight list, the lists to add and
#'   their coefficients.
#' @param pair What [jfp_mv_pair()] returns.
#' @param st A state, from [jfp_state()].
#' @param k The order, 1 to 4.
#'
#' @return As described.
#'
#' @keywords internal
jfp_live <- function(bs) which(!vapply(bs$rows, is.null, logical(1)))

#' @rdname jfp_live
#' @keywords internal
jfp_add_weights <- function(acc, lst, coef) {
  for (i in seq_along(lst)) {
    for (k in names(lst[[i]])) {
      x <- coef[i] * lst[[i]][[k]]
      acc[[k]] <- if (is.null(acc[[k]])) x else acc[[k]] + x
    }
  }
  acc
}

#' @rdname jfp_live
#' @keywords internal
jfp_add_pair <- function(acc, pair, coef) {
  acc$sparse <- jfp_add_weights(acc$sparse, list(pair$sparse), coef)
  for (z in pair$dense) {
    acc$dense[[length(acc$dense) + 1L]] <- list(idx = z$idx, X = coef * z$X)
  }
  acc
}

#' @rdname jfp_live
#' @keywords internal
jfp_darray <- function(st, k) {
  slot <- paste0("dE", k)
  if (is.null(st[[slot]])) {
    st[[slot]] <- switch(k,
      jfp_dexpected(st),
      jfp_d2expected(st),
      jfp_scale_apply(st, 3L),
      jfp_scale_apply(st, 4L))
  }
  st[[slot]]
}

#' @rdname jfp_live
#' @keywords internal
jfp_check_outer <- function(st) {
  jfp_check_ok(st)
  invisible(NULL)
}


#' The Weighted Cross Product of the Joint Rows
#'
#' @description
#' Assembles \eqn{\sum_{c,d}V_c^\top\mathrm{diag}(w\,W_{cd})V_d} over the
#' joint vector of a model with a filter, where \eqn{V_c} are the rows of
#' every equation in the joint vector ([jfp_joint_rows()]) and the weights
#' are keyed as [jfp_terms3()] returns them.
#'
#' @param st A state, from [jfp_state()], of a model with a filter.
#' @param sparse The list of per-observation weights.
#'
#' @return A base matrix over the joint vector.
#'
#' @keywords internal
jfp_joint_sparse <- function(st, sparse) {
  params <- st$params
  n <- st$spec@n_obs
  w <- st$spec@weights
  V <- st$joint$V
  out <- matrix(0, st$total, st$total)
  for (cc in seq_along(params)) {
    for (dd in seq_along(params)) {
      Wcd <- sparse[[hess_key(params, cc, dd)]]
      if (is.null(Wcd) || !any(Wcd != 0)) next
      out <- out + crossprod(V[[cc]], (w * rep_len(Wcd, n)) * V[[dd]])
    }
  }
  out
}
