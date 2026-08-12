#' @include outer_gradient.R
NULL

# The exact Hessian of the marginal criterion.
#
# Differentiating dV/dt_m = -rho_m - (1/2) tr(M K_m) once more:
#
#   d2V/dt_m dt_l = -rho_ml + b_m' J b_l
#                   + (1/2) tr(M K_l M K_m) - (1/2) tr(M dK_m/dt_l)
#
# The first two terms come from the envelope form of the gradient: rho_m is
# read at the mode, so its total derivative picks up c_m' b_l, and c_m = -J b_m
# turns that into b_m' J b_l.
#
# The determinant contributes twice. M itself moves, dM/dt_l = -M K_l M, which
# is the same expression whether M is K^-1 or A(A'KA)^-1 A'; and K_m moves,
#
#   dK_m/dt_l = S_ml + U[b_l, b_m] + T[b_ml],
#
# with T[v] the third derivative of the penalized objective in beta contracted
# once, U[v, u] the fourth contracted twice, and
#
#   b_ml = -J^-1 ( (S_l + T[b_l]) b_m + S_m b_l + c_ml ).
#
# Every term of dK_m/dt_l that would carry a third or fourth derivative of the
# PENALTY in beta is absent because the route requires beta_quadratic(); what
# is left of T and U is the log-likelihood's own third and fourth derivatives
# in the link-scale predictors, which distributions7 carries in closed form for
# every family.

#' The Exact Hessian of the Marginal Criterion
#'
#' @description
#' \eqn{\partial^2 V/\partial\eta^2} at the penalized mode, over the free scale
#' of the hyperparameters under estimation.
#'
#' @details
#' The pieces are those of the gradient differentiated once more: the
#' hyperparameter Hessian of the penalty, the movement of the mode read
#' through the penalized curvature, and two contributions from the
#' determinant, one from \eqn{M} moving and one from \eqn{K_m} moving.
#'
#' \strong{Everything the penalty contributes is asked of the penalty}
#' (\code{\link[penalties7]{penalty_hess_theta}},
#' \code{\link[penalties7]{penalty_dhessian}},
#' \code{\link[penalties7]{penalty_d2hessian}},
#' \code{\link[penalties7]{penalty_dcross}}), so a penalty that is not
#' quadratic in its hyperparameters -- a ridge, a random effect, a structured
#' prior -- is covered by the same assembly with no branch here.
#'
#' \strong{Onto the free scale} the chain rule is second order and diagonal,
#' each hyperparameter having its own link: with \eqn{\theta = h(\eta)},
#' \deqn{\partial^2 V/\partial\eta_m\partial\eta_l =
#'   h_m'h_l'\,\partial^2 V/\partial\theta_m\partial\theta_l
#'   + \delta_{ml}\,h_m''\,\partial V/\partial\theta_m.}
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An \code{\link{OuterMethod}}.
#' @param idx The outer index.
#' @param basis The integrated subspace, or \code{NULL}.
#'
#' @return A square matrix, one row per row of \code{idx}, or \code{NULL} where
#'   the determinant does not exist.
#'
#' @seealso \code{\link{statmod_marginal_grad}}, \code{\link{reml}}
#'
#' @keywords internal
statmod_marginal_hess <- function(spec, design, coef, hyper, method, idx,
                                  basis = NULL) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  nh <- nrow(idx)

  H <- statmod_information_at(spec, coef, design, expected = FALSE)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S[!is.finite(S)] <- 0
  K <- H + S
  Kfac <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(Kfac)) return(NULL)
  Kinv <- chol2inv(Kfac)
  M <- if (is.null(basis)) Kinv else {
    inner <- tryCatch(chol2inv(chol(crossprod(basis, K %*% basis))),
                      error = function(e) NULL)
    if (is.null(inner)) return(NULL)
    basis %*% inner %*% t(basis)
  }

  th <- statmod_eta(spec, design, coef)$theta
  d3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                       scale = "link")
  d4 <- distributions7::distrib_deriv4(spec@distrib, spec@response, th,
                                       scale = "link")

  # the pieces each hyperparameter owns, in the stacked coefficient space
  pieces <- outer_pieces(spec, design, coef, hyper, idx, offs, total)
  bhat <- lapply(seq_len(nh), function(m) -as.numeric(Kinv %*% pieces$c[[m]]))
  tv <- lapply(bhat, function(v) block_predictors(design, params, npar, offs,
                                                 v))
  Tm <- lapply(tv, function(t) contract3(spec, design, d3, params, npar, offs,
                                         total, t))
  Km <- lapply(seq_len(nh), function(m) pieces$S[[m]] + Tm[[m]])

  out <- matrix(0, nh, nh)
  for (m in seq_len(nh)) {
    for (l in m:nh) {
      key <- pieces$pair[m, l]
      Sml <- pieces$S2[[key]]
      cml <- pieces$c2[[key]]
      rhs <- (pieces$S[[l]] + Tm[[l]]) %*% bhat[[m]] +
        pieces$S[[m]] %*% bhat[[l]] + cml
      bml <- -as.numeric(Kinv %*% rhs)
      dKm <- Sml +
        contract4(spec, design, d4, params, npar, offs, total, tv[[l]],
                  tv[[m]]) +
        contract3(spec, design, d3, params, npar, offs, total,
                  block_predictors(design, params, npar, offs, bml))
      v <- -pieces$rho2[m, l] +
        sum(bhat[[m]] * as.numeric(K %*% bhat[[l]])) +
        sum((M %*% Km[[l]]) * t(M %*% Km[[m]])) / 2 -
        sum(M * dKm) / 2
      out[m, l] <- v
      out[l, m] <- v
    }
  }

  # and onto the free scale the search runs on
  links <- attr(idx, "links")
  g <- statmod_marginal_grad(spec, design, coef, hyper, method, idx, basis,
                             free = FALSE)
  h1 <- numeric(nh)
  h2 <- numeric(nh)
  for (r in seq_len(nh)) {
    v <- hyper[[idx$parameter[r]]][[idx$term[r]]][[idx$name[r]]]
    e <- linkfunctions7::linkfun(links[[r]], v)
    h1[r] <- linkfunctions7::dlinkinv(links[[r]], e)
    h2[r] <- linkfunctions7::d2linkinv(links[[r]], e)
  }
  out <- out * outer(h1, h1)
  diag(out) <- diag(out) + h2 * g
  out
}


#' The Per-Hyperparameter Pieces of the Outer Derivatives
#'
#' @description
#' The penalty's contributions to the criterion's gradient and Hessian, placed
#' in the stacked coefficient space.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param idx The outer index.
#' @param offs,total The block offsets and the total width.
#' @param order \code{1} for the first-order pieces alone, \code{2} for the
#'   second-order ones as well.
#'
#' @details
#' At order 1 the second-order generics are not called at all. That is what
#' lets a penalty supplying only \code{penalty_dhessian()} give an exact
#' gradient: asking it for a derivative the gradient does not use would have
#' rejected it for a quantity nobody wanted.
#'
#' @return A list with \code{S} (one matrix per hyperparameter) and \code{c}
#'   (one vector), and at order 2 also \code{S2} and \code{c2} (one per pair,
#'   keyed), \code{rho2} (the hyperparameter Hessian) and \code{pair} (the key
#'   of each pair).
#'
#' @keywords internal
outer_pieces <- function(spec, design, coef, hyper, idx, offs, total,
                         order = 2L) {
  params <- spec@distrib@params
  nh <- nrow(idx)
  Sm <- vector("list", nh)
  cm <- vector("list", nh)
  S2 <- list()
  c2 <- list()
  rho2 <- matrix(0, nh, nh)
  pair <- matrix("", nh, nh)

  terms <- unique(paste(idx$parameter, idx$term, sep = "\r"))
  for (s in terms) {
    bits <- strsplit(s, "\r", fixed = TRUE)[[1L]]
    p <- bits[1L]
    nm <- bits[2L]
    a <- match(p, params)
    un <- statmod_unit(spec, design, p, nm)
    cols <- un$cols
    pos <- un$index
    pen <- un$penalty
    bt <- coef[[p]][cols]
    th <- as.list(hyper[[p]][[nm]])
    dS <- penalties7::penalty_dhessian(pen, bt, th)
    cr <- penalties7::penalty_cross(pen, bt, th)
    rows <- which(idx$parameter == p & idx$term == nm)
    for (r in rows) {
      Sm[[r]] <- matrix(0, total, total)
      Sm[[r]][pos, pos] <- dS[[idx$name[r]]]
      cm[[r]] <- numeric(total)
      cm[[r]][pos] <- as.numeric(cr[[idx$name[r]]])
    }
    if (order < 2L) next
    d2S <- penalties7::penalty_d2hessian(pen, bt, th)
    dcr <- penalties7::penalty_dcross(pen, bt, th)
    ht <- penalties7::penalty_hess_theta(pen, bt, th)
    for (r in rows) {
      for (q in rows) {
        key <- pair_key(idx$name[r], idx$name[q], names(ht))
        pair[r, q] <- key
        rho2[r, q] <- as.numeric(ht[[key]])[1L]
        if (is.null(S2[[key]])) {
          S2[[key]] <- matrix(0, total, total)
          S2[[key]][pos, pos] <- d2S[[key]]
          v <- numeric(total)
          v[pos] <- as.numeric(dcr[[key]])
          c2[[key]] <- v
        }
      }
    }
  }
  if (order < 2L) return(list(S = Sm, c = cm))
  # two hyperparameters of DIFFERENT terms are independent: the penalty is a
  # sum, so no second derivative mixes them
  zeroS <- matrix(0, total, total)
  zeroc <- numeric(total)
  for (m in seq_len(nh)) {
    for (l in seq_len(nh)) {
      if (nzchar(pair[m, l])) next
      key <- paste0("\r", m, "_", l)
      pair[m, l] <- key
      S2[[key]] <- zeroS
      c2[[key]] <- zeroc
    }
  }
  list(S = Sm, c = cm, S2 = S2, c2 = c2, rho2 = rho2, pair = pair)
}


#' The Key of a Hyperparameter Pair
#'
#' @description
#' Locates the entry of a penalty's hyperparameter Hessian, which is keyed by
#' name and not by position.
#'
#' @param a,b The two hyperparameter names.
#' @param keys The names the penalty actually returned.
#'
#' @return A single string.
#'
#' @keywords internal
pair_key <- function(a, b, keys) {
  k <- paste(a, b, sep = "_")
  if (k %in% keys) return(k)
  k <- paste(b, a, sep = "_")
  if (k %in% keys) return(k)
  stop(sprintf("No hyperparameter component for '%s' and '%s'.", a, b),
       call. = FALSE)
}


#' The Predictors a Coefficient Direction Induces
#'
#' @description
#' \eqn{(X_k v_k)_i} for each distribution parameter, which is how a movement
#' of the coefficients is felt by the log-density.
#'
#' @param design The design.
#' @param params The parameter names.
#' @param npar,offs The block sizes and offsets.
#' @param v A stacked coefficient vector.
#'
#' @return A list of numeric vectors, one per parameter.
#'
#' @keywords internal
block_predictors <- function(design, params, npar, offs, v) {
  lapply(seq_along(params), function(k) {
    if (npar[k] == 0L) return(NULL)
    as.numeric(design[[params[k]]]$X %*% v[offs[k] + seq_len(npar[k])])
  })
}


#' The Third Derivative of the Objective Contracted Once
#'
#' @description
#' \eqn{T[v] = (\partial K/\partial\beta)\cdot v}, a matrix over the stacked
#' coefficients.
#'
#' @details
#' Each block is a weighted crossproduct, the weight being
#' \eqn{w_i\sum_k \ell'''_{abk}(X_k v_k)_i}: the third derivative never
#' appears as an array.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param d3 The third derivatives in the link scale.
#' @param params,npar,offs,total The block bookkeeping.
#' @param tv The predictors of the direction.
#'
#' @return A square matrix.
#'
#' @keywords internal
contract3 <- function(spec, design, d3, params, npar, offs, total, tv) {
  keys <- names(d3)
  n <- spec@n_obs
  out <- matrix(0, total, total)
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in a:length(params)) {
      if (npar[b] == 0L) next
      w <- numeric(n)
      for (k in seq_along(params)) {
        if (npar[k] == 0L) next
        w <- w + rep_len(d3[[d3_key(params, a, b, k, keys)]], n) * tv[[k]]
      }
      blk <- -crossprod(design[[params[a]]]$X * (spec@weights * w),
                        design[[params[b]]]$X)
      ra <- offs[a] + seq_len(npar[a])
      rb <- offs[b] + seq_len(npar[b])
      out[ra, rb] <- blk
      if (a != b) out[rb, ra] <- t(blk)
    }
  }
  out
}


#' The Fourth Derivative of the Objective Contracted Twice
#'
#' @description
#' \eqn{U[v, u] = (\partial^2 K/\partial\beta^2)\cdot(v, u)}, a matrix over the
#' stacked coefficients.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param d4 The fourth derivatives in the link scale.
#' @param params,npar,offs,total The block bookkeeping.
#' @param tv,tu The predictors of the two directions.
#'
#' @return A square matrix.
#'
#' @keywords internal
contract4 <- function(spec, design, d4, params, npar, offs, total, tv, tu) {
  keys <- names(d4)
  n <- spec@n_obs
  out <- matrix(0, total, total)
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in a:length(params)) {
      if (npar[b] == 0L) next
      w <- numeric(n)
      for (k in seq_along(params)) {
        if (npar[k] == 0L) next
        for (q in seq_along(params)) {
          if (npar[q] == 0L) next
          w <- w + rep_len(d4[[d4_key(params, a, b, k, q, keys)]], n) *
            tv[[k]] * tu[[q]]
        }
      }
      blk <- -crossprod(design[[params[a]]]$X * (spec@weights * w),
                        design[[params[b]]]$X)
      ra <- offs[a] + seq_len(npar[a])
      rb <- offs[b] + seq_len(npar[b])
      out[ra, rb] <- blk
      if (a != b) out[rb, ra] <- t(blk)
    }
  }
  out
}


#' The Name of a Fourth-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, k, q)} entry, by a name BUILT from the parameter
#' names in the family's own order.
#'
#' @param params The parameter names.
#' @param a,b,k,q Indices into \code{params}.
#' @param keys The names the derivative returned.
#'
#' @return A single string.
#'
#' @keywords internal
d4_key <- function(params, a, b, k, q, keys) {
  want <- paste(params[sort(c(a, b, k, q))], collapse = "_")
  if (!want %in% keys) {
    stop(sprintf("No fourth-derivative component '%s'.", want), call. = FALSE)
  }
  want
}


#' The Smoothing-Parameter Correction to the Effective Degrees of Freedom
#'
#' @description
#' The amount by which \eqn{\mathrm{tr}[(H+S)^{-1}H]} understates the
#' complexity of a fit whose hyperparameters were themselves estimated.
#'
#' @details
#' The ordinary effective degrees of freedom read the smoothing parameters
#' as though they were known, and they were not: they were chosen from the
#' same data. Propagating their uncertainty into the coefficients gives the
#' corrected Bayesian covariance
#'
#' \deqn{V' = V_\beta + J V_\theta J^\top, \qquad
#'   J = \partial\hat\beta/\partial\theta,}
#'
#' and the corrected count is \eqn{\mathrm{tr}(V' H)}. \eqn{J} comes from
#' the implicit function theorem at the penalized mode: differentiating
#' \eqn{\partial(-\ell + \rho)/\partial\beta = 0} in the hyperparameter
#' gives \eqn{(H+S)J_k = -\partial^2\rho/\partial\beta\,\partial\theta_k},
#' whose right-hand side is \pkg{penalties7}'s \code{penalty_cross()}.
#' \eqn{V_\theta} is the inverse of the outer criterion's own Hessian,
#' which \code{\link{statmod_marginal_hess}} returns, negated because that
#' criterion is a maximand.
#'
#' Everything is on the hyperparameter's LINK scale, which is where the
#' outer criterion optimizes and therefore the only scale on which its
#' Hessian is a variance.
#'
#' \strong{Against mgcv.} This is the first of the two terms mgcv sums into
#' \code{edf2}, and it agrees with mgcv's to about 3e-4 on a univariate
#' smooth once the difference between the two bases is allowed for. mgcv
#' adds a second term for the Gaussian scale, which it profiles out of the
#' fit and whose uncertainty it must therefore add back; here every
#' distribution parameter carries its own equation and its own coefficients,
#' so that uncertainty is already inside \eqn{H}. The residual difference is
#' measured at 0.16, 0.10 and 0.05 effective parameters at n = 200, 400 and
#' 2000, falling with the sample size.
#'
#' \strong{Where it does not apply.} A kinked penalty has no hyperparameter
#' the outer criterion estimates -- \code{\link{outer_hyper_index}} skips it
#' -- so there is nothing to propagate and the correction is zero. That is
#' not an approximation: the map from the hyperparameter to the penalized
#' mode turns a corner whenever a coefficient joins or leaves the active
#' set, and a delta method needs a derivative that does not exist there.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param design The design.
#' @param method The outer method that estimated them, or \code{NULL}.
#' @param expected Whether the information is the expected one.
#' @param approx The approximation for the expected information.
#'
#' @return A list with \code{total}, the scalar correction, and \code{per},
#'   one entry per penalty key. Zero throughout where no hyperparameter was
#'   estimated.
#'
#' @references
#' Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and model
#' selection for general smooth models. \emph{Journal of the American
#' Statistical Association}, 111(516), 1548--1563.
#'
#' @seealso \code{\link{statmod_marginal_hess}},
#'   \code{\link[penalties7]{penalty_cross}}
#'
#' @keywords internal
statmod_edf_correction <- function(spec, coef, hyper, design, method,
                                   expected = TRUE, approx = "bartlett") {
  zero <- list(total = 0, per = numeric(0))
  if (is.null(method) || !method@kind %in% c("ml", "reml")) return(zero)

  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  if (!nrow(idx)) return(zero)

  H <- statmod_information_at(spec, coef, design, expected, approx)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S[!is.finite(S)] <- 0
  Vb <- tryCatch(solve(H + S), error = function(e) NULL)
  if (is.null(Vb)) return(zero)

  # J = -Vb %*% d2rho/dbeta dtheta, one column per estimated hyperparameter
  J <- matrix(0, nrow(H), nrow(idx))
  for (un in statmod_penalized(spec, design)) {
    cr <- tryCatch(
      penalties7::penalty_cross(un$penalty, coef[[un$param]][un$cols],
                                as.list(hyper[[un$param]][[un$key]]),
                                scale = "link"),
      error = function(e) NULL)
    if (is.null(cr)) next
    for (h in names(cr)) {
      k <- which(idx$parameter == un$param & idx$term == un$key &
                   idx$name == h)
      if (!length(k)) next
      J[un$index, k] <- as.numeric(cr[[h]])
    }
  }
  J <- -Vb %*% J

  Ho <- tryCatch(statmod_marginal_hess(spec, design, coef, hyper, method,
                                       idx, NULL),
                 error = function(e) NULL)
  if (is.null(Ho)) return(zero)
  # the criterion is a maximand, so its Hessian is negative definite at the
  # optimum and the variance is the inverse of its negative
  Vth <- tryCatch(solve(-as.matrix(Ho)), error = function(e) NULL)
  if (is.null(Vth)) return(zero)

  # per hyperparameter, so a summary can say which penalty the extra
  # complexity belongs to, and in total
  contrib <- vapply(seq_len(nrow(idx)), function(k) {
    v <- J[, k, drop = FALSE] %*% Vth[k, k, drop = FALSE] %*% t(J[, k, drop = FALSE])
    sum(v * t(H))
  }, numeric(1))
  per <- tapply(contrib, paste(idx$parameter, idx$term, sep = "\r"), sum)
  total <- sum((J %*% Vth %*% t(J)) * t(H))
  list(total = total, per = per)
}
