#' @include joint_filter.R
NULL

# ------------------------------------------------- the scaled driving vector
#
# Matrices over the filtered parameters are held per order k as arrays
# n x K x K x C_k, one matrix per observation and per multiset of k
# predictor indices; the compiled functions of src/joint_scaling.cpp, one per
# order, compute every multiset of their order.

.jsc_cache <- new.env(parent = emptyenv())

#' The Multisets of Predictor Indices and Their Sub-Multisets
#'
#' @description
#' Returns, for \eqn{P} predictors and the orders 0 to `order`, the multisets
#' of each size, the table of sub-multisets that the compiled functions of
#' [joint_scale_derivs()] read, and the position of each ordered tuple's
#' multiset.
#'
#' @details
#' The multisets of size \eqn{k} are the sorted tuples over \eqn{1, \dots,
#' P}, in lexicographic order. The table of order \eqn{m} has one row per
#' multiset \eqn{B} and one column per subset of its positions, a bit mask
#' with position \eqn{j} as bit \eqn{j - 1}; the entry is the zero-based
#' index of the sub-multiset among the multisets of its size. The map of
#' order \eqn{m} gives, for each tuple in the order of [expand.grid()] (the
#' first index running fastest), the one-based index of its multiset. The
#' result is cached by \eqn{P} and `order`.
#'
#' @param P The number of predictors.
#' @param order The highest size.
#'
#' @return A list with `ms` (the multisets, one matrix per size from 0),
#'   `tab` (one integer matrix per order from 1) and `map` (one integer
#'   vector per order from 1).
#'
#' @keywords internal
joint_scale_tables <- function(P, order) {
  ck <- paste(P, order, sep = ".")
  hit <- .jsc_cache[[ck]]
  if (!is.null(hit)) return(hit)
  key <- function(x) paste(x, collapse = ".")
  ms <- lapply(0:order, function(k) {
    if (!k) return(matrix(integer(0), 1L, 0L))
    rows <- as.matrix(do.call(expand.grid, rep(list(seq_len(P)), k)))
    rows <- rows[apply(rows, 1L, function(r) !is.unsorted(r)), , drop = FALSE]
    rows <- rows[do.call(base::order, unname(as.data.frame(rows))), , drop = FALSE]
    unname(rows)
  })
  keys <- lapply(ms, function(x) apply(x, 1L, key))
  keys[[1L]] <- ""
  tab <- lapply(seq_len(order), function(m) {
    B <- ms[[m + 1L]]
    out <- matrix(0L, nrow(B), 2^m)
    for (b in seq_len(nrow(B))) {
      for (mask in 0:(2^m - 1L)) {
        inA <- bitwAnd(mask, 2^(seq_len(m) - 1L)) > 0
        k <- sum(inA)
        out[b, mask + 1L] <- if (!k) 0L else
          match(key(sort(B[b, inA])), keys[[k + 1L]]) - 1L
      }
    }
    out
  })
  map <- lapply(seq_len(order), function(m) {
    tup <- as.matrix(do.call(expand.grid, rep(list(seq_len(P)), m)))
    match(apply(tup, 1L, function(r) key(sort(r))), keys[[m + 1L]])
  })
  out <- list(ms = ms, tab = tab, map = map)
  .jsc_cache[[ck]] <- out
  out
}


#' The Derivatives of the Scaling Matrix in the Predictors
#'
#' @description
#' Returns the derivatives of \eqn{S = \mathcal{I}_{GG}^{-d}} in every
#' multiset of predictor indices up to the order given, one matrix per
#' observation, for \eqn{d} equal to 1/2 or 1.
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
#' ([distributions7::distrib_dexpected_hessian()] and its siblings to
#' [distributions7::distrib_d5expected_hessian()], on the link scale). The
#' eigendecomposition (LAPACK's `dsyev`, one per observation) and each order
#' are compiled, one function per order with its terms written out.
#'
#' @param spec A [StatmodSpec()].
#' @param theta The per-observation parameters at the fitted predictors.
#' @param G The positions of the filtered parameters.
#' @param d The power, 1/2 or 1.
#' @param order The highest order of derivative of \eqn{\mathcal{I}} read.
#'
#' @return A list with `S`, the derivatives of orders 0 to `order` (arrays
#'   \eqn{n \times K \times K \times C_k} with \eqn{C_k} the number of
#'   multisets of size \eqn{k}, in the order of [joint_scale_tables()]), and
#'   `tables`, what [joint_scale_tables()] returns.
#'
#' @keywords internal
joint_scale_derivs <- function(spec, theta, G, d, order) {
  dist <- spec@distrib
  y <- spec@response
  params <- dist@params
  P <- length(params)
  n <- spec@n_obs
  K <- length(G)
  rn <- function(x) rep_len(as.numeric(x), n)
  tables <- joint_scale_tables(P, order)
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
  # the information's derivatives of order k, one matrix per multiset
  I_of <- function(k) {
    B <- tables$ms[[k + 1L]]
    out <- array(0, c(n, K, K, nrow(B)))
    for (b in seq_len(nrow(B))) {
      for (a in seq_len(K)) for (c in a:K) {
        v <- if (!k) {
          -rn(E0[[hess_key(params, G[a], G[c])]])
        } else {
          -rn(Es[[k]][[do.call(keyf[[k]], c(list(params, G[a], G[c]),
                                             as.list(B[b, ])))]])
        }
        out[, a, c, b] <- v
        out[, c, a, b] <- v
      }
    }
    out
  }
  s0 <- jsc_s0_cpp(I_of(0L), n, K, d)
  if (!s0$ok) {
    stop("the expected information of the filtered parameters is not positive definite.",
         call. = FALSE)
  }
  S <- list(s0$S)
  M <- list()
  for (k in seq_len(order)) {
    Ik <- I_of(k)
    tab <- tables$tab[[k]]
    res <- switch(k,
      jsc_s1_cpp(Ik, M, S, tab, s0$Q, s0$r, n, K, d),
      jsc_s2_cpp(Ik, M, S, tab, s0$Q, s0$r, n, K, d),
      jsc_s3_cpp(Ik, M, S, tab, s0$Q, s0$r, n, K, d),
      jsc_s4_cpp(Ik, M, S, tab, s0$Q, s0$r, n, K, d),
      jsc_s5_cpp(Ik, M, S, tab, s0$Q, s0$r, n, K, d))
    M[[k]] <- res$M
    S[[k + 1L]] <- res$S
  }
  list(S = S, tables = tables)
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
#' from [joint_scale_derivs()]. Each order is compiled, one function per
#' order; the entries are symmetric in \eqn{T}, so each multiset is computed
#' once and copied to its permutations.
#'
#' @param spec A [StatmodSpec()].
#' @param tm The joint filter.
#' @param theta The per-observation parameters at the fitted predictors.
#' @param gl,H,D3,D4,D5 The family's link-scale derivatives.
#' @param order As for [joint_drive_data()].
#' @param D6 The family's sixth link-scale derivative, read at order 5.
#'
#' @return The list that [joint_drive_data()] returns.
#'
#' @keywords internal
joint_drive_scaled <- function(spec, tm, theta, gl, H, D3, D4, D5, order,
                               D6 = NULL) {
  params <- spec@distrib@params
  P <- length(params)
  G <- match(tm@params, params)
  K <- length(G)
  n <- spec@n_obs
  d <- tm@terms[[1L]]@scaling
  rn <- function(x) rep_len(as.numeric(x), n)
  sc <- joint_scale_derivs(spec, theta, G, d, order)
  tables <- sc$tables
  # the family's derivatives with index c and a multiset of size k
  fam <- list(H, D3, D4, D5, D6)
  keyf <- list(hess_key, deriv3_key, deriv4_key, deriv5_key, deriv6_key)
  g_of <- function(k) {
    B <- tables$ms[[k + 1L]]
    out <- array(0, c(n, K, nrow(B)))
    for (b in seq_len(nrow(B))) {
      for (c in seq_len(K)) {
        out[, c, b] <- if (!k) rn(gl[[params[G[c]]]]) else
          rn(fam[[k]][[do.call(keyf[[k]], c(list(params, G[c]),
                                              as.list(B[b, ])))]])
      }
    }
    out
  }
  g <- lapply(0:order, g_of)
  u <- matrix(jsc_u0_cpp(sc$S[[1L]], g[[1L]], n, K), n, K)
  U <- list()
  for (r in seq_len(order)) {
    tab <- tables$tab[[r]]
    Sr <- sc$S[seq_len(r + 1L)]
    gr <- g[seq_len(r + 1L)]
    Uc <- switch(r,
      jsc_u1_cpp(Sr, gr, tab, n, K), jsc_u2_cpp(Sr, gr, tab, n, K),
      jsc_u3_cpp(Sr, gr, tab, n, K), jsc_u4_cpp(Sr, gr, tab, n, K),
      jsc_u5_cpp(Sr, gr, tab, n, K))
    Uc <- array(Uc, c(n, K, nrow(tab)))[, , tables$map[[r]], drop = FALSE]
    U[[r]] <- array(Uc, c(n, K, rep(P, r)))
  }
  list(filtered = G, u = u, U1 = U[[1L]],
       U2 = if (order >= 2L) U[[2L]], U3 = if (order >= 3L) U[[3L]],
       U4 = if (order >= 4L) U[[4L]], U5 = if (order >= 5L) U[[5L]])
}
