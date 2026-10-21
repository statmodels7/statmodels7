#' @include joint_filter.R
NULL

# ------------------------------------------------- the scaled driving vector
#
# Matrices over the filtered parameters are held as arrays n x K x K, one
# matrix per observation, so that every step below is one vectorized
# operation over the observations.

.jm_mult <- function(A, B) {
  K <- dim(A)[2L]
  out <- array(0, dim(A))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    s <- 0
    for (k in seq_len(K)) s <- s + A[, i, k] * B[, k, j]
    out[, i, j] <- s
  }
  out
}

.jm_t <- function(A) aperm(A, c(1L, 3L, 2L))

#' The Derivatives of the Scaling Matrix in the Predictors
#'
#' @description
#' Returns a function of a multiset of predictor indices that gives the
#' derivative of \eqn{S = \mathcal{I}_{GG}^{-d}} in those predictors, one
#' matrix per observation, for \eqn{d} equal to 1/2 or 1.
#'
#' @details
#' With \eqn{M = \mathcal{I}} for \eqn{d = 1} and \eqn{M = \mathcal{I}^{1/2}}
#' for \eqn{d = 1/2}, \eqn{S = M^{-1}}, and Leibniz's rule over the positions
#' of the multiset \eqn{B} on \eqn{MS = 1} gives
#' \deqn{S_B = -S \sum_{\emptyset \ne A \subseteq B} M_A S_{B \setminus A}.}
#' For the symmetric root the same rule on \eqn{XX = \mathcal{I}} gives
#' \deqn{X X_B + X_B X = \mathcal{I}_B
#'   - \sum_{\emptyset \ne A \subsetneq B} X_A X_{B \setminus A},}
#' a Sylvester equation, which is solved in the eigenbasis of
#' \eqn{\mathcal{I}} by dividing by \eqn{\sqrt{\lambda_i} + \sqrt{\lambda_j}}.
#' The derivatives of \eqn{\mathcal{I}} are the family's
#' ([distributions7::distrib_dexpected_hessian()] and its siblings, on the
#' link scale). Every derivative is memoized by its multiset.
#'
#' @param spec A [StatmodSpec()].
#' @param theta The per-observation parameters at the fitted predictors.
#' @param G The positions of the filtered parameters.
#' @param d The power, 1/2 or 1.
#' @param order The highest order of derivative of \eqn{\mathcal{I}} read.
#'
#' @return A function of an integer vector of predictor positions.
#'
#' @keywords internal
joint_scale_derivs <- function(spec, theta, G, d, order) {
  dist <- spec@distrib
  y <- spec@response
  params <- dist@params
  n <- spec@n_obs
  K <- length(G)
  rn <- function(x) rep_len(as.numeric(x), n)
  E0 <- distributions7::distrib_expected_hessian(dist, y, theta, scale = "link")
  Es <- list()
  if (order >= 1L) Es[[1L]] <- distributions7::distrib_dexpected_hessian(
    dist, y, theta, scale = "link")
  if (order >= 2L) Es[[2L]] <- distributions7::distrib_d2expected_hessian(
    dist, y, theta, scale = "link")
  if (order >= 3L) Es[[3L]] <- distributions7::distrib_d3expected_hessian(
    dist, y, theta, scale = "link")
  if (order >= 4L) Es[[4L]] <- distributions7::distrib_d4expected_hessian(
    dist, y, theta, scale = "link")
  keyf <- list(distributions7::dexpected_key, distributions7::d2expected_key,
               distributions7::d3expected_key, distributions7::d4expected_key)
  # the information's derivative in a multiset of predictors
  I_of <- function(B) {
    out <- array(0, c(n, K, K))
    for (a in seq_len(K)) for (b in seq_len(K)) {
      out[, a, b] <- if (!length(B)) {
        -rn(E0[[hess_key(params, G[a], G[b])]])
      } else {
        -rn(Es[[length(B)]][[do.call(keyf[[length(B)]],
                                     c(list(params, G[a], G[b]),
                                       as.list(B)))]])
      }
    }
    out
  }
  I0 <- I_of(integer(0))
  # the eigenbasis of the information, observation by observation
  Q <- array(0, c(n, K, K))
  r <- matrix(0, n, K)
  for (i in seq_len(n)) {
    ev <- eigen(matrix(I0[i, , ], K, K), symmetric = TRUE)
    if (any(ev$values <= 0)) {
      stop("the expected information of the filtered parameters is not positive definite.",
           call. = FALSE)
    }
    Q[i, , ] <- ev$vectors
    r[i, ] <- if (d == 0.5) sqrt(ev$values) else ev$values
  }
  diag_of <- function(v) {
    out <- array(0, c(n, K, K))
    for (k in seq_len(K)) out[, k, k] <- v[, k]
    out
  }
  M0 <- .jm_mult(.jm_mult(Q, diag_of(r)), .jm_t(Q))
  S0 <- .jm_mult(.jm_mult(Q, diag_of(1 / r)), .jm_t(Q))
  den <- array(0, c(n, K, K))
  for (a in seq_len(K)) for (b in seq_len(K)) den[, a, b] <- r[, a] + r[, b]
  key <- function(B) paste(sort(B), collapse = ".")
  subsets <- function(m) {
    if (!m) return(list())
    lapply(seq_len(2^m - 1L), function(mask)
      which(bitwAnd(mask, 2^(seq_len(m) - 1L)) > 0))
  }
  memo_M <- new.env()
  memo_S <- new.env()
  M_of <- function(B) {
    if (!length(B)) return(M0)
    k <- key(B)
    if (!is.null(memo_M[[k]])) return(memo_M[[k]])
    v <- if (d == 1) I_of(B) else {
      R <- I_of(B)
      for (A in subsets(length(B))) {
        if (length(A) == length(B)) next
        R <- R - .jm_mult(M_of(B[A]), M_of(B[-A]))
      }
      .jm_mult(.jm_mult(Q, .jm_mult(.jm_mult(.jm_t(Q), R), Q) / den),
               .jm_t(Q))
    }
    memo_M[[k]] <- v
    v
  }
  S_of <- function(B) {
    if (!length(B)) return(S0)
    k <- key(B)
    if (!is.null(memo_S[[k]])) return(memo_S[[k]])
    acc <- array(0, c(n, K, K))
    for (A in subsets(length(B))) {
      acc <- acc + .jm_mult(M_of(B[A]),
                            if (length(A) == length(B)) S0 else S_of(B[-A]))
    }
    v <- -.jm_mult(S0, acc)
    memo_S[[k]] <- v
    v
  }
  S_of
}


#' The Scaled Driving Vector's Derivatives in Every Predictor
#'
#' @description
#' Computes the entries of [joint_drive_data()] for a joint filter whose
#' `scaling` is 1/2 or 1.
#'
#' @details
#' With \eqn{u_l = \sum_c S_{lc}\, g_c} and \eqn{g_c} the link-scale
#' gradient, Leibniz's rule over the positions of a multiset \eqn{T} of
#' predictors gives
#' \deqn{U_{l;T} = \sum_{A \subseteq T} \sum_c S_{lc;T \setminus A}\,
#'   g_{c;A},}
#' where \eqn{g_{c;A}} is the family's derivative of order \eqn{1 + |A|}
#' with indices \eqn{c} and \eqn{A}, and the derivatives of \eqn{S} come
#' from [joint_scale_derivs()]. The entries are symmetric in \eqn{T}, so each
#' multiset is computed once and copied to its permutations.
#'
#' @param spec A [StatmodSpec()].
#' @param tm The joint filter.
#' @param theta The per-observation parameters at the fitted predictors.
#' @param gl,H,D3,D4,D5 The family's link-scale derivatives.
#' @param order As for [joint_drive_data()].
#'
#' @return The list that [joint_drive_data()] returns.
#'
#' @keywords internal
joint_drive_scaled <- function(spec, tm, theta, gl, H, D3, D4, D5, order) {
  params <- spec@distrib@params
  P <- length(params)
  G <- match(tm@params, params)
  K <- length(G)
  n <- spec@n_obs
  d <- tm@terms[[1L]]@scaling
  rn <- function(x) rep_len(as.numeric(x), n)
  S_of <- joint_scale_derivs(spec, theta, G, d, order)
  # the family's derivative with index c and the multiset A
  g_of <- function(c, A) {
    m <- length(A)
    idx <- c(G[c], A)
    rn(switch(m + 1L,
              gl[[params[G[c]]]],
              H[[hess_key(params, idx[1L], idx[2L])]],
              D3[[do.call(deriv3_key, c(list(params), as.list(idx)))]],
              D4[[do.call(deriv4_key, c(list(params), as.list(idx)))]],
              D5[[do.call(deriv5_key, c(list(params), as.list(idx)))]]))
  }
  U_of <- function(T_) {
    m <- length(T_)
    out <- matrix(0, n, K)
    for (mask in 0:(2^m - 1L)) {
      inA <- bitwAnd(mask, 2^(seq_len(m) - 1L)) > 0
      S <- S_of(T_[!inA])
      for (c in seq_len(K)) {
        gc <- g_of(c, T_[inA])
        for (l in seq_len(K)) out[, l] <- out[, l] + S[, l, c] * gc
      }
    }
    out
  }
  u <- U_of(integer(0))
  U <- list()
  for (r in seq_len(order)) {
    A <- array(0, c(n, K, rep(P, r)))
    tup <- as.matrix(do.call(expand.grid, rep(list(seq_len(P)), r)))
    done <- new.env()
    for (k in seq_len(nrow(tup))) {
      T_ <- as.integer(tup[k, ])
      kk <- paste(sort(T_), collapse = ".")
      v <- done[[kk]]
      if (is.null(v)) {
        v <- U_of(sort(T_))
        done[[kk]] <- v
      }
      for (l in seq_len(K)) {
        A[cbind(seq_len(n), l, matrix(T_, n, r, byrow = TRUE))] <- v[, l]
      }
    }
    U[[r]] <- A
  }
  list(filtered = G, u = u, U1 = U[[1L]],
       U2 = if (order >= 2L) U[[2L]], U3 = if (order >= 3L) U[[3L]],
       U4 = if (order >= 4L) U[[4L]])
}
