# The quantity that drives a score-driven filter, and its derivatives.
#
# gas(scaling = d) drives its recursion by u = s * I^(-d), where s = l_p is the
# score in the filter's own predictor and I = -E[l_pp] the expected information
# of that predictor, both on the link scale. The recursions of modelterms7 read
# the driving quantity and its derivatives in every predictor through the same
# arrays that carry the family's derivatives (gl, H, D3, D4, D5), and only the
# entries whose multiset contains the filter's parameter: H[p, q] as the
# derivative of the driving quantity in eta_q, D3[p, r, r2] as its second
# derivative, and so on. filter_driving() therefore returns copies of those
# arrays with exactly those entries replaced by the derivatives of u; every
# other entry, and every array the caller uses for the information itself or
# for the weights of a contraction, stays the log-likelihood's.
#
# At d = 0 the arrays are returned as they came, the same objects, so the
# default route is the one the package ran before the scaling existed.

#' The Scaling Exponent of a Filter
#'
#' @description
#' The `scaling` of a [modelterms7::gas()] term, and zero for a structural
#' term that has none.
#'
#' @param tm A structural term.
#'
#' @return A single number.
#'
#' @keywords internal
filter_scaling <- function(tm) {
  d <- tryCatch(tm@scaling, error = function(e) NULL)
  if (is.null(d) || !length(d)) 0 else d
}


#' Check That a Family Supports a Scaled Score
#'
#' @description
#' Signals an error unless the family is univariate and carries an analytic
#' second derivative of its expected information, which the derivatives of a
#' filter driven by a scaled score reach (the fourth from a stencil on it).
#'
#' @param d A distribution object.
#' @param y The response.
#' @param theta The parameters, one value per observation or one in all.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
check_filter_scaling <- function(d, y, theta) {
  if (S7::S7_inherits(d, distributions7::multivariate_distrib)) {
    stop("a scaled score is available for univariate families only.",
         call. = FALSE)
  }
  th1 <- lapply(theta, function(v) v[1L])
  ok <- tryCatch({
    distributions7::distrib_d2expected_hessian(d, y[1L], th1, scale = "link")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    stop(sprintf(paste0(
      "gas(scaling = ) needs the expected information of '%s' and its\n",
      "  analytic second derivative, which the family does not carry."),
      d@distrib_name), call. = FALSE)
  }
  invisible(NULL)
}


#' The Driving Quantity of a Filter and Its Derivatives
#'
#' @description
#' Replaces, in the family's link-scale derivative arrays, the entries that a
#' filter's recursion reads as derivatives of its driving quantity with the
#' derivatives of the scaled score \eqn{u = s\,\mathcal{I}^{-d}}.
#'
#' @details
#' With \eqn{w = \mathcal{I}^{-d}}, Leibniz's rule over the positions of an
#' index tuple \eqn{T} gives
#' \deqn{u_T = \sum_{A \subseteq T} s_{A}\, w_{T \setminus A},}
#' where \eqn{s_A} is the family's derivative of order \eqn{1 + |A|} with the
#' filter's parameter as one index, and Faa di Bruno's formula over the set
#' partitions \eqn{\pi} of \eqn{B} gives
#' \deqn{w_B = \sum_{\pi} \phi^{(|\pi|)}(\mathcal{I})
#'   \prod_{\beta \in \pi} \mathcal{I}_\beta, \qquad
#'   \phi^{(k)}(x) = (-d)(-d-1)\cdots(-d-k+1)\, x^{-d-k}.}
#' The derivatives \eqn{\mathcal{I}_\beta} are minus those of
#' [distributions7::distrib_dexpected_hessian()] and its higher-order
#' siblings, all on the link scale. The highest array given sets the order:
#' `D5` needs the fourth derivative of the information.
#'
#' @param spec A [StatmodSpec()].
#' @param theta The per-observation parameters.
#' @param ap The index of the parameter the filter sits in.
#' @param scaling The exponent \eqn{d}.
#' @param gl,H,D3,D4,D5 The family's link-scale derivatives; the higher ones
#'   may be `NULL`.
#'
#' @return A list with `gl`, `H`, `D3`, `D4` and `D5`, the same objects where
#'   `scaling` is zero.
#'
#' @keywords internal
filter_driving <- function(spec, theta, ap, scaling, gl, H, D3 = NULL,
                           D4 = NULL, D5 = NULL) {
  if (scaling == 0) return(list(gl = gl, H = H, D3 = D3, D4 = D4, D5 = D5))
  d <- spec@distrib
  y <- spec@response
  n <- spec@n_obs
  params <- d@params
  np <- length(params)
  order <- if (!is.null(D5)) 4L else if (!is.null(D4)) 3L else
    if (!is.null(D3)) 2L else if (!is.null(H)) 1L else 0L
  rn <- function(x) rep_len(as.numeric(x), n)
  tkey <- function(idx) paste(sort(idx), collapse = ".")

  info <- -rn(distributions7::distrib_expected_hessian(
    d, y, theta, scale = "link")[[hess_key(params, ap, ap)]])
  Ider <- list()
  if (order >= 1L) {
    E1 <- distributions7::distrib_dexpected_hessian(d, y, theta, scale = "link")
    for (c1 in seq_len(np)) {
      Ider[[tkey(c1)]] <- -rn(E1[[distributions7::dexpected_key(params, ap, ap, c1)]])
    }
  }
  # the non-decreasing index tuples of length m, in the order the family's
  # own component names enumerate them
  tuples <- function(m) {
    if (m == 0L) return(list(integer(0)))
    idx <- as.matrix(do.call(expand.grid, rep(list(seq_len(np)), m)))
    idx <- idx[, rev(seq_len(m)), drop = FALSE]
    idx <- idx[apply(idx, 1L, function(r) all(diff(r) >= 0)), , drop = FALSE]
    lapply(seq_len(nrow(idx)), function(k) as.integer(idx[k, ]))
  }
  if (order >= 2L) {
    E2 <- distributions7::distrib_d2expected_hessian(d, y, theta, scale = "link")
    for (tu in tuples(2L)) {
      Ider[[tkey(tu)]] <- -rn(E2[[distributions7::d2expected_key(
        params, ap, ap, tu[1L], tu[2L])]])
    }
  }
  if (order >= 3L) {
    E3 <- distributions7::distrib_d3expected_hessian(d, y, theta, scale = "link")
    for (tu in tuples(3L)) {
      Ider[[tkey(tu)]] <- -rn(E3[[distributions7::d3expected_key(
        params, ap, ap, tu[1L], tu[2L], tu[3L])]])
    }
  }
  if (order >= 4L) {
    E4 <- distributions7::distrib_d4expected_hessian(d, y, theta, scale = "link")
    for (tu in tuples(4L)) {
      Ider[[tkey(tu)]] <- -rn(E4[[distributions7::d4expected_key(
        params, ap, ap, tu[1L], tu[2L], tu[3L], tu[4L])]])
    }
  }

  # phi^(k)(I) = (-d)(-d-1)...(-d-k+1) I^(-d-k)
  phik <- function(k) {
    co <- if (k == 0L) 1 else prod(-scaling - seq.int(0L, k - 1L))
    co * info^(-scaling - k)
  }
  w_of <- function(B) {
    if (!length(B)) return(phik(0L))
    out <- numeric(n)
    for (part in numericals7::set_partitions(length(B))) {
      term <- phik(length(part))
      for (beta in part) term <- term * Ider[[tkey(B[beta])]]
      out <- out + term
    }
    out
  }
  # the score's derivatives are read from the arrays as they came in, which
  # the assignments below replace one order at a time
  s0 <- list(gl = gl, H = H, D3 = D3, D4 = D4, D5 = D5)
  s_of <- function(A) {
    switch(length(A) + 1L,
           rn(s0$gl[[params[ap]]]),
           rn(s0$H[[hess_key(params, ap, A[1L])]]),
           rn(s0$D3[[deriv3_key(params, ap, A[1L], A[2L])]]),
           rn(s0$D4[[deriv4_key(params, ap, A[1L], A[2L], A[3L])]]),
           rn(s0$D5[[deriv5_key(params, ap, A[1L], A[2L], A[3L], A[4L])]]))
  }
  u_of <- function(Tt) {
    m <- length(Tt)
    out <- numeric(n)
    for (mask in 0:(2^m - 1L)) {
      inA <- bitwAnd(mask, 2^(seq_len(m) - 1L)) > 0
      out <- out + s_of(Tt[inA]) * w_of(Tt[!inA])
    }
    out
  }

  gl[[params[ap]]] <- u_of(integer(0))
  if (order >= 1L) {
    for (q in seq_len(np)) H[[hess_key(params, ap, q)]] <- u_of(q)
  }
  if (order >= 2L) {
    for (tu in tuples(2L)) D3[[deriv3_key(params, ap, tu[1L], tu[2L])]] <- u_of(tu)
  }
  if (order >= 3L) {
    for (tu in tuples(3L)) {
      D4[[deriv4_key(params, ap, tu[1L], tu[2L], tu[3L])]] <- u_of(tu)
    }
  }
  if (order >= 4L) {
    for (tu in tuples(4L)) {
      D5[[deriv5_key(params, ap, tu[1L], tu[2L], tu[3L], tu[4L])]] <- u_of(tu)
    }
  }
  list(gl = gl, H = H, D3 = D3, D4 = D4, D5 = D5)
}
