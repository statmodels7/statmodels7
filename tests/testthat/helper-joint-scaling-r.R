# The R implementation of joint_scale_derivs() and joint_drive_scaled() as
# they stood before the compiled route (statmodels7 0.241.0), kept as the
# reference that test-joint-scaling.R compares the compiled route against.
# S_B is built by Leibniz's rule over the positions of the multiset with one
# eigen() per observation; see R/joint_scaling.R for the formulas.

.jm_mult_r <- function(A, B) {
  K <- dim(A)[2L]
  out <- array(0, dim(A))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    s <- 0
    for (k in seq_len(K)) s <- s + A[, i, k] * B[, k, j]
    out[, i, j] <- s
  }
  out
}

.jm_t_r <- function(A) aperm(A, c(1L, 3L, 2L))

joint_scale_derivs_r <- function(spec, theta, G, d, order) {
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
  if (order >= 5L) Es[[5L]] <- distributions7::distrib_d5expected_hessian(
    dist, y, theta, scale = "link")
  keyf <- list(distributions7::dexpected_key, distributions7::d2expected_key,
               distributions7::d3expected_key, distributions7::d4expected_key,
               distributions7::d5expected_key)
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
  M0 <- .jm_mult_r(.jm_mult_r(Q, diag_of(r)), .jm_t_r(Q))
  S0 <- .jm_mult_r(.jm_mult_r(Q, diag_of(1 / r)), .jm_t_r(Q))
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
        R <- R - .jm_mult_r(M_of(B[A]), M_of(B[-A]))
      }
      .jm_mult_r(.jm_mult_r(Q, .jm_mult_r(.jm_mult_r(.jm_t_r(Q), R), Q) / den),
               .jm_t_r(Q))
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
      acc <- acc + .jm_mult_r(M_of(B[A]),
                            if (length(A) == length(B)) S0 else S_of(B[-A]))
    }
    v <- -.jm_mult_r(S0, acc)
    memo_S[[k]] <- v
    v
  }
  S_of
}


joint_drive_scaled_r <- function(spec, tm, theta, gl, H, D3, D4, D5, order,
                               D6 = NULL) {
  params <- spec@distrib@params
  P <- length(params)
  G <- match(tm@params, params)
  K <- length(G)
  n <- spec@n_obs
  d <- tm@terms[[1L]]@scaling
  rn <- function(x) rep_len(as.numeric(x), n)
  S_of <- joint_scale_derivs_r(spec, theta, G, d, order)
  # the family's derivative with index c and the multiset A
  g_of <- function(c, A) {
    m <- length(A)
    idx <- c(G[c], A)
    rn(switch(m + 1L,
              gl[[params[G[c]]]],
              H[[hess_key(params, idx[1L], idx[2L])]],
              D3[[do.call(deriv3_key, c(list(params), as.list(idx)))]],
              D4[[do.call(deriv4_key, c(list(params), as.list(idx)))]],
              D5[[do.call(deriv5_key, c(list(params), as.list(idx)))]],
              D6[[do.call(deriv6_key, c(list(params), as.list(idx)))]]))
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
       U4 = if (order >= 4L) U[[4L]], U5 = if (order >= 5L) U[[5L]])
}
