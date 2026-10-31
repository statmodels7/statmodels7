#' @include jfp_engine.R
NULL

# The third and fourth derivatives of log|I_tilde| in the joint vector, for
# blocks whose rows do not move, in a model without a filter and without a
# block that is recomputed at the coefficients. Every quantity is written in
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
#' \eqn{\tilde p^2} columns. They are written for blocks whose rows do not
#' move, in a model without a filter and without a term whose block is
#' recomputed at the coefficients, and they signal an error elsewhere.
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
#'   \eqn{\sum_{c,d}X_c^\top\mathrm{diag}(w\,W_{cd})X_d}), and `lr`, a list of
#'   triples `X`, `K`, `Y` (the matrix is \eqn{XKY^\top}).
#'
#' @keywords internal
jfp_terms3 <- function(st, v) {
  jfp_check_outer(st)
  ev <- jfp_eta_dir(st, v)
  out <- list(sparse = list(), lr = list())
  for (bs in st$blocks) {
    a <- bs$block$power
    p <- ncol(bs$Z[[jfp_live(bs)[1L]]])
    Ip <- diag(p)
    Av <- jfp_white(st, bs, list(ev))
    B <- jfp_rows_B(st, bs, list())
    Bv <- jfp_rows_B(st, bs, list(ev))
    out$sparse <- jfp_add_weights(out$sparse, list(
      jfp_pair_weights(st, bs, Ip, list(ev)),
      jfp_pair_weights(st, bs, Av, list())), c(a, -a))
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
  ev <- jfp_eta_dir(st, v)
  ew <- jfp_eta_dir(st, w)
  out <- list(sparse = list(), lr = list())
  for (bs in st$blocks) {
    a <- bs$block$power
    p <- ncol(bs$Z[[jfp_live(bs)[1L]]])
    Ip <- diag(p)
    I2 <- diag(p * p)
    Av <- jfp_white(st, bs, list(ev))
    Aw <- jfp_white(st, bs, list(ew))
    Avw <- jfp_white(st, bs, list(ev, ew))
    B <- jfp_rows_B(st, bs, list())
    Bv <- jfp_rows_B(st, bs, list(ev))
    Bw <- jfp_rows_B(st, bs, list(ew))
    Bvw <- jfp_rows_B(st, bs, list(ev, ew))
    kr <- kronecker
    out$sparse <- jfp_add_weights(out$sparse, list(
      jfp_pair_weights(st, bs, Ip, list(ev, ew)),
      jfp_pair_weights(st, bs, Aw, list(ev)),
      jfp_pair_weights(st, bs, Avw, list()),
      jfp_pair_weights(st, bs, Av, list(ew)),
      jfp_pair_weights(st, bs, Aw %*% Av + Av %*% Aw, list())),
      a * c(1, -1, -1, -1, 1))
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
  out <- as_dense(statmod_information_at(
    st$spec, st$coef, st$design, expected = TRUE, emethod = st$emethod,
    H = lapply(terms$sparse, function(x) -x)))
  for (z in terms$lr) out <- out + z$X %*% z$K %*% t(z$Y)
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
  for (z in terms$lr) {
    s <- s + sum((crossprod(z$X, M %*% z$Y)) * z$K)
  }
  s
}


#' The Pieces From Which the Terms Are Built
#'
#' @description
#' `jfp_eta_dir()` returns the movement of each equation's predictor along a
#' direction; `jfp_dI()` a derivative of the family's expected information,
#' contracted along directions; `jfp_white()`, `jfp_rows_B()` and
#' `jfp_pair_weights()` the three shapes in which a derivative of the
#' block's matrix enters (see the comment at the top of `R/jfp_outer.R`).
#'
#' @param st A state, from [jfp_state()].
#' @param bs A block's entry of the state.
#' @param v A direction over the joint vector.
#' @param q,r The pair of the information.
#' @param fixed The equations of the derivatives that are not contracted.
#' @param dirs A list of directions, each as `jfp_eta_dir()` returns it.
#' @param Y A symmetric \eqn{\tilde p\times\tilde p} matrix.
#'
#' @return `jfp_eta_dir()` a list of numeric vectors, one per distribution
#'   parameter; `jfp_dI()` a numeric vector over the observations;
#'   `jfp_white()` a \eqn{\tilde p\times\tilde p} matrix; `jfp_rows_B()` a
#'   matrix over the joint vector and the \eqn{\tilde p^2} entries;
#'   `jfp_pair_weights()` a list keyed as [statmod_information_at()] reads
#'   its `H`.
#'
#' @keywords internal
jfp_eta_dir <- function(st, v) {
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

#' @rdname jfp_eta_dir
#' @keywords internal
jfp_white <- function(st, bs, dirs) {
  live <- jfp_live(bs)
  w <- st$spec@weights
  p <- ncol(bs$Z[[live[1L]]])
  out <- matrix(0, p, p)
  for (q in live) {
    for (r in live) {
      out <- out + crossprod(bs$Z[[q]] * (w * jfp_dI(st, q, r, integer(0),
                                                    dirs)), bs$Z[[r]])
    }
  }
  (out + t(out)) / 2
}

#' @rdname jfp_eta_dir
#' @keywords internal
jfp_rows_B <- function(st, bs, dirs) {
  live <- jfp_live(bs)
  w <- st$spec@weights
  n <- st$spec@n_obs
  p <- ncol(bs$Z[[live[1L]]])
  ii <- rep(seq_len(p), times = p)
  jj <- rep(seq_len(p), each = p)
  B <- matrix(0, st$total, p * p)
  for (cc in seq_along(st$params)) {
    if (st$npar[cc] == 0L) next
    Fc <- matrix(0, n, p * p)
    for (q in live) {
      for (r in live) {
        s <- jfp_dI(st, q, r, cc, dirs)
        if (!any(s != 0)) next
        Fc <- Fc + s * bs$Z[[q]][, ii, drop = FALSE] *
          bs$Z[[r]][, jj, drop = FALSE]
      }
    }
    B[st$offs[cc] + seq_len(st$npar[cc]), ] <-
      as.matrix(crossprod(st$design[[st$params[cc]]]$X, w * Fc))
  }
  B
}

#' @rdname jfp_eta_dir
#' @keywords internal
jfp_pair_weights <- function(st, bs, Y, dirs) {
  live <- jfp_live(bs)
  params <- st$params
  np <- length(params)
  y <- list()
  for (q in live) {
    for (r in live) {
      y[[paste(q, r)]] <- rowSums((bs$Z[[q]] %*% Y) * bs$Z[[r]])
    }
  }
  out <- list()
  for (cc in seq_len(np)) {
    for (dd in cc:np) {
      v <- numeric(st$spec@n_obs)
      for (q in live) {
        for (r in live) {
          v <- v + y[[paste(q, r)]] * jfp_dI(st, q, r, c(cc, dd), dirs)
        }
      }
      out[[hess_key(params, cc, dd)]] <- v
    }
  }
  out
}


#' Small Helpers of the Third and Fourth Derivatives
#'
#' @description
#' `jfp_live()` returns the equations in which a block has rows;
#' `jfp_add_weights()` sums weight lists with coefficients;
#' `jfp_darray()` returns the family's derivative of the expected
#' information of order `k`, kept in the state; `jfp_check_outer()` signals
#' an error where the third and fourth derivatives are not written.
#'
#' @param bs A block's entry of the state.
#' @param acc,lst,coef An accumulated weight list, the lists to add and
#'   their coefficients.
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
jfp_darray <- function(st, k) {
  slot <- paste0("dE", k)
  if (is.null(st[[slot]])) {
    d <- st$spec@distrib
    y <- st$spec@response
    th <- st$ev$theta
    st[[slot]] <- switch(k,
      jfp_dexpected(st),
      jfp_d2expected(st),
      distributions7::distrib_d3expected_hessian(d, y, th, scale = "link",
                                                 method = st$emethod),
      distributions7::distrib_d4expected_hessian(d, y, th, scale = "link",
                                                 method = st$emethod))
  }
  st[[slot]]
}

#' @rdname jfp_live
#' @keywords internal
jfp_check_outer <- function(st) {
  jfp_check_ok(st)
  if (!is.null(st$joint)) {
    stop("The third derivative of a Jeffreys penalty beside a filter is not ",
         "implemented yet.", call. = FALSE)
  }
  rf <- attr(st$design, "refresh")
  if (!is.null(rf) && length(rf)) {
    stop(paste0("The third derivative of a Jeffreys penalty beside a block ",
                "recomputed at the\n  coefficients is not implemented yet."),
         call. = FALSE)
  }
  invisible(NULL)
}
