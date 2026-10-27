#' @include joint_outer.R
NULL

#' The Derivative of the Joint Matrix Along a Direction, Joint Filter
#'
#' @description
#' Computes, for a model carrying a joint filter, the matrix that
#' [structural_dk_matrix()] computes for one filter: the part of
#' \eqn{\partial K/\partial u[v]} that comes from the family and the
#' recursion.
#'
#' @details
#' Each term of the single-filter matrix that reads the filter's equation is
#' a sum over the filtered parameters here. The third derivative of the levels
#' and their re-weighted curvature come from one call each to the joint
#' recursion, with one column of weights per filtered parameter, and the
#' movement of every filtered row enters through the cross derivatives of its
#' own parameter.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows.
#' @param st The shared quantities, from [joint_grad_parts()].
#' @param v The direction, over the estimated coordinates.
#'
#' @return A square symmetric matrix over the estimated coordinates.
#'
#' @keywords internal
joint_dk_matrix <- function(spec, design, jd, st, v) {
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
  out <- -cv3$curvature[keep, keep, drop = FALSE]
  kappa <- matrix(0, n, length(GF))
  for (k in seq_along(GF)) {
    for (b in seq_along(params)) {
      kappa[, k] <- kappa[, k] +
        rep_len(st$H[[hess_key(params, GF[k], b)]], n) * dV[[b]]
    }
  }
  cvk <- modelterms7::term_curvature(f$tm, f$eta_static, NULL, NULL, NULL, f$psi,
                                     w * kappa, seed, st$bd)
  out <- out - cvk$curvature[keep, keep, drop = FALSE]

  keys3 <- names(st$D3)
  for (a in seq_along(params)) {
    for (b in seq_along(params)) {
      cab <- numeric(n)
      for (k in seq_along(params)) {
        cab <- cab + rep_len(st$D3[[d3_key(params, a, b, k, keys3)]], n) *
          dV[[k]]
      }
      out <- out - crossprod(st$Vk[[a]] * (w * cab), st$Vk[[b]])
    }
  }
  for (k in seq_along(GF)) {
    dphi <- cv3$dphi[[k]][[1L]][, keep, drop = FALSE]
    for (b in seq_along(params)) {
      hb <- rep_len(st$H[[hess_key(params, GF[k], b)]], n)
      Z <- crossprod(dphi * (w * hb), st$Vk[[b]])
      out <- out - Z - t(Z)
    }
  }
  (out + t(out)) / 2
}


#' The Joint Second Derivative of the Penalized Information, Joint Filter
#'
#' @description
#' Returns, for a model carrying a joint filter, the quantity that
#' [structural_chain_extra2()] returns for one filter:
#' \eqn{\mathrm{tr}(M\,\partial^2 K/\partial u^2[v, w])} less its design
#' part.
#'
#' @details
#' The terms are those of the single-filter function, with every quantity
#' that reads the filter's equation summed over the filtered parameters. The
#' second derivative of the score of parameter \eqn{k} along the two
#' directions is
#' \deqn{\sum_{a,b}\ell_{kab}(V_a\cdot v)(V_b\cdot w)
#'   + \sum_j \ell_{kj}\, v^\top\Phi_j w,}
#' in which only the filtered equations have a second derivative
#' \eqn{\Phi_j}. The movements of the filtered rows along the two directions
#' pair every two filtered parameters through \eqn{\ell_{jj'}}.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows.
#' @param M The matrix against which the trace is taken.
#' @param st The shared quantities, from [joint_grad_parts()] at order 4.
#' @param v,w The two directions, over the estimated coordinates.
#'
#' @return A single number.
#'
#' @keywords internal
joint_chain_extra2 <- function(spec, design, jd, M, st, v, w) {
  params <- jd$params
  n <- jd$n
  f <- jd$f
  wt <- st$w
  keep <- jd$keep
  GF <- st$GF
  K <- length(GF)
  vfull <- numeric(jd$mf)
  vfull[keep] <- v
  wfull <- numeric(jd$mf)
  wfull[keep] <- w
  seed <- stats::setNames(jd$V, params)
  dVv <- lapply(st$Vk, function(x) as.numeric(x %*% v))
  dVw <- lapply(st$Vk, function(x) as.numeric(x %*% w))
  H <- st$H
  D3 <- st$D3
  D4 <- st$D4
  keys3 <- names(D3)
  hk <- function(a, b) rep_len(H[[hess_key(params, a, b)]], n)

  cv4 <- modelterms7::term_fourth(f$tm, f$eta_static, NULL, NULL, NULL, f$psi, st$gw,
                                  seed, st$bd, list(vfull, wfull))
  dphi_v <- lapply(seq_len(K), function(k)
    cv4$dphi[[k]][[1L]][, keep, drop = FALSE])
  dphi_w <- lapply(seq_len(K), function(k)
    cv4$dphi[[k]][[2L]][, keep, drop = FALSE])
  dpsi <- lapply(seq_len(K), function(k) cv4$dpsi[[k]][, keep, drop = FALSE])
  phi_vw <- lapply(seq_len(K), function(k) as.numeric(dphi_w[[k]] %*% v))

  # the levels' term at the fourth order
  tot <- sum(M * cv4$curvature[keep, keep, drop = FALSE])

  # the levels' term re-weighted by the second derivative of each score
  gk <- matrix(0, n, K)
  for (k in seq_len(K)) {
    for (a in seq_along(params)) {
      for (b in seq_along(params)) {
        gk[, k] <- gk[, k] +
          rep_len(D3[[d3_key(params, GF[k], a, b, keys3)]], n) *
          dVv[[a]] * dVw[[b]]
      }
    }
    for (j in seq_len(K)) gk[, k] <- gk[, k] + hk(GF[k], GF[j]) * phi_vw[[j]]
  }
  cv0 <- modelterms7::term_curvature(f$tm, f$eta_static, NULL, NULL, NULL, f$psi,
                                     wt * gk, seed, st$bd)
  tot <- tot + sum(M * cv0$curvature[keep, keep, drop = FALSE])

  # the levels' third derivative along one direction, weighted by how each
  # score moves along the other
  kap <- function(dv) {
    out <- matrix(0, n, K)
    for (k in seq_len(K)) {
      for (b in seq_along(params)) out[, k] <- out[, k] + hk(GF[k], b) * dv[[b]]
    }
    out
  }
  for (pair in list(list(kap(dVv), wfull), list(kap(dVw), vfull))) {
    cv3 <- modelterms7::term_third(f$tm, f$eta_static, NULL, NULL, NULL, f$psi,
                                   wt * pair[[1L]], seed, st$bd, pair[[2L]])
    tot <- tot + sum(M * cv3$curvature[keep, keep, drop = FALSE])
  }

  # the family's fourth and third derivatives against the leverage diagonal
  for (a in seq_along(params)) {
    for (b in seq_along(params)) {
      co <- numeric(n)
      for (k in seq_along(params)) {
        for (k2 in seq_along(params)) {
          co <- co + rep_len(D4[[deriv4_key(params, a, b, k, k2)]], n) *
            dVv[[k]] * dVw[[k2]]
        }
      }
      for (j in seq_len(K)) {
        co <- co + rep_len(D3[[d3_key(params, a, b, GF[j], keys3)]], n) *
          phi_vw[[j]]
      }
      tot <- tot + sum(wt * co * st$G[[a]][[b]])
    }
  }

  # a filtered row moves along one direction inside a third derivative read
  # along the other
  for (dd in list(list(dVv, dphi_w), list(dVw, dphi_v))) {
    for (j in seq_len(K)) {
      for (b in seq_along(params)) {
        co <- numeric(n)
        for (k in seq_along(params)) {
          co <- co + rep_len(D3[[d3_key(params, GF[j], b, k, keys3)]], n) *
            dd[[1L]][[k]]
        }
        tot <- tot + 2 * sum(wt * co * rowSums(st$VM[[b]] * dd[[2L]][[j]]))
      }
    }
  }

  # the second movement of a filtered row, which is the third derivative of
  # its predictor along both directions
  for (j in seq_len(K)) {
    for (b in seq_along(params)) {
      tot <- tot + 2 * sum(wt * hk(GF[j], b) *
                             rowSums(st$VM[[b]] * dpsi[[j]]))
    }
  }

  # two filtered rows moving at once, one along each direction
  for (j in seq_len(K)) {
    Mv <- dphi_v[[j]] %*% M
    for (j2 in seq_len(K)) {
      tot <- tot + 2 * sum(wt * hk(GF[j], GF[j2]) *
                             rowSums(Mv * dphi_w[[j2]]))
    }
  }
  -tot
}
