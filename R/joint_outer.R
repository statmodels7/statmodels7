#' @include joint_predict.R
NULL

# The exact outer derivatives of a model carrying a joint filter. The
# marginal criterion's matrix is
#   K = -sum_t w_t sum_{a,b} l_ab V_a' V_b - sum_t w_t sum_k l_k Phi_k,
# with V_k the joint filter's Jacobian for each filtered equation and Phi_k
# the second derivative of its level, so every piece of the single-filter
# gradient becomes a sum over the filtered parameters, and the recursions are
# the joint ones, read with weight matrices with one column per filtered
# parameter.

#' The Direction-Free Quantities of the Joint Chain Term, Joint Filter
#'
#' @description
#' Computes, for a model carrying a joint filter, what
#' [structural_grad_parts()] computes for one filter: the family's
#' derivatives, the rows of every equation with the joint Jacobian in place,
#' their products with \eqn{M}, the per-observation diagonals of \eqn{M} and
#' the contraction \eqn{u}.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param jd The joint rows, from [joint_design_rows()].
#' @param M The matrix against which the trace is taken.
#' @param order The highest order of the driving vector's derivatives read:
#'   3 for the gradient, 4 for the Hessian.
#'
#' @return A list with `u`, `V`, `Vk`, `VM`, `G`, `H`, `D3`, `D4`, `D5`,
#'   `gl`, `gw` (the weights of the recursion, one column per filtered
#'   parameter), `bd` (the driving data), `GF` (the positions of the
#'   filtered parameters) and `w`.
#'
#' @keywords internal
joint_grad_parts <- function(spec, design, coef, jd, M, order = 3L) {
  params <- jd$params
  n <- jd$n
  w <- spec@weights
  keep <- jd$keep
  ev <- jd$ev
  f <- jd$f
  d <- spec@distrib
  y <- spec@response
  th <- ev$theta
  gl <- distributions7::distrib_gradient(d, y, th, scale = "link",
                                         threads = spec@threads)
  H <- distributions7::distrib_hessian(d, y, th, scale = "link",
                                       threads = spec@threads)
  D3 <- distributions7::distrib_deriv3(d, y, th, scale = "link",
                                       threads = spec@threads)
  D4 <- distributions7::distrib_deriv4(d, y, th, scale = "link",
                                       threads = spec@threads)
  D5 <- if (order >= 4L) distributions7::distrib_deriv5(
    d, y, th, scale = "link", threads = spec@threads)
  GF <- match(f$params, params)
  gw <- matrix(vapply(f$params, function(p) w * rep_len(as.numeric(gl[[p]]), n),
                      numeric(n)), n, length(GF))
  bd <- joint_drive_data(spec, f$tm, gl, H, D3, D4, D5, order = order,
                         theta = th)
  cv <- joint_filter_curvature(spec, design, f, jd$V, gl, H, D3, th)
  V <- jd$V
  for (k in seq_along(GF)) V[[GF[k]]] <- cv$jacobian[[k]]
  Vk <- lapply(V, function(x) x[, keep, drop = FALSE])
  VM <- lapply(Vk, function(x) x %*% M)
  G <- vector("list", length(params))
  for (a in seq_along(params)) {
    G[[a]] <- lapply(seq_along(params), function(b) rowSums(VM[[a]] * Vk[[b]]))
  }
  keys3 <- names(D3)
  u <- numeric(length(keep))
  for (k in seq_along(params)) {
    s <- numeric(n)
    for (a in seq_along(params)) {
      for (b in seq_along(params)) {
        s <- s + rep_len(D3[[d3_key(params, a, b, k, keys3)]], n) * G[[a]][[b]]
      }
    }
    u <- u - as.numeric(crossprod(Vk[[k]], w * s))
  }
  list(u = u, V = V, Vk = Vk, VM = VM, G = G, H = H, D3 = D3, D4 = D4,
       D5 = D5, gl = gl, gw = gw, bd = bd, GF = GF, w = w)
}


#' The Chain Term of the Joint Gradient Along a Direction, Joint Filter
#'
#' @description
#' Computes, for a model carrying a joint filter, the three contributions
#' that [structural_chain_extra()] computes for one filter: the third
#' derivative of the levels along \eqn{v} contracted against the weighted
#' score, the curvature of the levels re-weighted by the derivative of the
#' score along \eqn{v}, and the movement of the joint Jacobian through every
#' cross derivative of a filtered equation.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows.
#' @param M The matrix against which the trace is taken.
#' @param st The shared quantities, from [joint_grad_parts()].
#' @param v The direction, over the estimated coordinates.
#'
#' @return A single number.
#'
#' @keywords internal
joint_chain_extra <- function(spec, design, jd, M, st, v) {
  params <- jd$params
  n <- jd$n
  f <- jd$f
  w <- st$w
  keep <- jd$keep
  GF <- st$GF
  vfull <- numeric(jd$mf)
  vfull[keep] <- v
  seed <- stats::setNames(jd$V, params)
  dV <- lapply(st$Vk, function(x) as.numeric(x %*% v))

  cv3 <- modelterms7::term_third(f$tm, f$eta_static, NULL, NULL, NULL, f$psi, st$gw,
                                 seed, st$bd, vfull)
  W3 <- cv3$curvature[keep, keep, drop = FALSE]

  # the weight of each level's curvature moves with the score of its own
  # parameter along v
  kappa <- matrix(0, n, length(GF))
  for (k in seq_along(GF)) {
    for (b in seq_along(params)) {
      kappa[, k] <- kappa[, k] +
        rep_len(st$H[[hess_key(params, GF[k], b)]], n) * dV[[b]]
    }
  }
  cvk <- modelterms7::term_curvature(f$tm, f$eta_static, NULL, NULL, NULL, f$psi,
                                     w * kappa, seed, st$bd)
  Wk <- cvk$curvature[keep, keep, drop = FALSE]

  # each filtered row moves by Phi_k v and enters through every l_kb; the two
  # orderings of a pair are transposes and contribute equally
  b2 <- numeric(n)
  for (k in seq_along(GF)) {
    dphi <- cv3$dphi[[k]][[1L]][, keep, drop = FALSE]
    for (b in seq_along(params)) {
      b2 <- b2 + rep_len(st$H[[hess_key(params, GF[k], b)]], n) *
        rowSums(st$VM[[b]] * dphi)
    }
  }
  -2 * sum(w * b2) - sum(M * Wk) - sum(M * W3)
}
