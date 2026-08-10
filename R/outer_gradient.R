#' @include outer.R
NULL

# The exact gradient of the marginal criterion.
#
# Writing V(t) = l(b) - rho(b; t) + (q/2) log 2pi - (1/2) log|A'KA| at the
# penalized mode b(t), with K = H + S:
#
#   dV/dt_m = -drho/dt_m - (1/2) [ tr(M dS/dt_m) + u' v_m ]
#
# where M = A(A'KA)^-1 A' (which is K^-1 for REML), v_m = db/dt_m and
# u_c = tr(M dK/db_c).
#
# The first two terms of V contribute only -drho/dt_m, by the envelope
# theorem: their derivative through b carries the factor dl/db - drho/db,
# which is the score of the penalized objective and vanishes at the mode.
# Nothing is left of db/dt there, and it survives only inside the determinant.
#
# v_m comes from differentiating the stationarity condition:
#   (H_obs + S) db/dt_m + d2rho/db dt_m = 0.
#
# u is the contraction of the third derivative of the log-likelihood in the
# link-scale predictors with the diagonal of the "hat" matrix of M, one
# crossprod per distribution parameter, so nothing of order p^3 beyond the
# factorization of K that the criterion has already paid for.

#' Can the Exact Gradient Be Computed Here?
#'
#' @description
#' \code{TRUE} when every hyperparameter under estimation belongs to a penalty
#' whose second derivative in the coefficients is linear in the
#' hyperparameters and free of the coefficients, and the criterion uses the
#' observed information.
#'
#' @details
#' \strong{Why the observed information.} \eqn{K} enters the criterion through
#' its determinant, so the gradient needs \eqn{\partial K/\partial\beta}. With
#' the observed information that is the third derivative of the log-likelihood
#' in the link-scale predictors, which every family of \pkg{distributions7}
#' carries in closed form. With the expected information it would be the
#' derivative in \eqn{\beta} of \eqn{-E[\ell'']}, which is not
#' \eqn{-E[\ell''']} and is not one of that package's generics. So
#' \code{reml(hessian = "observed")} is what the exact route asks for, and
#' \code{"expected"} keeps the derivative-free search.
#'
#' \strong{Why linear.} \eqn{\partial S/\partial\theta} is not a generic of
#' \pkg{penalties7}: what that package exposes is the penalty, its gradient,
#' its Hessian and the mixed block, not the third derivative
#' \eqn{\partial^3\rho/\partial\beta^2\partial\theta}. For a penalty whose
#' Hessian is linear in the hyperparameters the derivative is recoverable from
#' the Hessian itself and nothing has to be differentiated; for one that is
#' not -- a ridge, whose Hessian carries \eqn{1/\sigma^2}, or any penalty built
#' from a density -- it is not, and the search stays derivative-free. Closing
#' that case is a generic in \pkg{penalties7}, and for the density branch it
#' would need \eqn{\partial^3\ell/\partial y^2\partial\theta} from
#' \pkg{distributions7}, which does not exist either.
#'
#' The linearity is \strong{checked and not assumed}: a penalty is asked for
#' its Hessian at \eqn{\theta} and at \eqn{2\theta}, and at two coefficient
#' vectors, and admitted only if the first doubles and the second does not
#' move.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param idx The outer index, from \code{\link{outer_hyper_index}}.
#' @param method An \code{\link{OuterMethod}}.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @seealso \code{\link{statmod_marginal_grad}}
#'
#' @keywords internal
outer_gradient_ok <- function(spec, design, idx, method) {
  if (!identical(method@hessian, "observed")) return(FALSE)
  if (!nrow(idx)) return(FALSE)
  seen <- unique(paste(idx$parameter, idx$term, sep = "\r"))
  for (s in seen) {
    bits <- strsplit(s, "\r", fixed = TRUE)[[1L]]
    pen <- modelterms7::term_penalty(spec@terms[[bits[1L]]][[bits[2L]]])
    if (!penalty_hessian_linear(pen)) return(FALSE)
  }
  TRUE
}


#' Is a Penalty's Hessian Linear in the Hyperparameters?
#'
#' @description
#' Asks the penalty rather than assuming: its Hessian must double when the
#' hyperparameters double, and must not move when the coefficients do.
#'
#' @param pen A \pkg{penalties7} penalty.
#'
#' @return A single logical.
#'
#' @keywords internal
penalty_hessian_linear <- function(pen) {
  th <- penalty_theta_start(pen)
  k <- as.integer(pen@n_coef)
  if (!length(th) || is.na(k) || k < 1L) return(FALSE)
  b1 <- rep(0.37, k)
  b2 <- rep(-1.21, k)
  ok <- tryCatch({
    S1 <- penalties7::penalty_hessian(pen, b1, as.list(th))
    S2 <- penalties7::penalty_hessian(pen, b1, as.list(th * 2))
    S3 <- penalties7::penalty_hessian(pen, b2, as.list(th))
    isTRUE(all.equal(S2, 2 * S1, tolerance = 1e-10)) &&
      isTRUE(all.equal(S3, S1, tolerance = 1e-10)) &&
      any(abs(S1) > 0)
  }, error = function(e) FALSE)
  isTRUE(ok)
}


#' The Derivative of a Penalty's Hessian in Its Hyperparameters
#'
#' @description
#' One matrix per hyperparameter, for a penalty whose Hessian is linear in
#' them.
#'
#' @details
#' Under that linearity \eqn{S(\theta + h e_m) = S(\theta) + h\,\partial
#' S/\partial\theta_m} holds for every \eqn{h}, so the difference below is
#' exact arithmetic and not an approximation: there is no truncation error to
#' shrink and no step to choose. It is written as a difference rather than by
#' evaluating at a unit vector because a hyperparameter of zero is outside the
#' bounds a penalty validates against, while \eqn{2\theta_m} is not.
#'
#' The linearity is verified by \code{\link{penalty_hessian_linear}} before
#' this is called.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param beta The term's coefficients.
#' @param theta The term's hyperparameters.
#'
#' @return A named list of matrices.
#'
#' @keywords internal
penalty_dhessian <- function(pen, beta, theta) {
  th <- as.list(theta)
  S0 <- penalties7::penalty_hessian(pen, beta, th)
  stats::setNames(lapply(names(th), function(m) {
    h <- as.numeric(th[[m]])
    tp <- th
    tp[[m]] <- h * 2
    (penalties7::penalty_hessian(pen, beta, tp) - S0) / h
  }), names(th))
}


#' The Exact Gradient of the Marginal Criterion
#'
#' @description
#' \eqn{\partial V/\partial\eta} at the penalized mode, over the free scale of
#' the hyperparameters under estimation.
#'
#' @details
#' The three pieces are the envelope term \eqn{-\partial\rho/\partial\theta},
#' the explicit derivative of the determinant \eqn{\mathrm{tr}(M\,\partial
#' S/\partial\theta)}, and the implicit one \eqn{u'v}, where
#' \eqn{v = -(H+S)^{-1}\partial^2\rho/\partial\beta\partial\theta} is how the
#' mode moves and \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)} is how
#' the determinant reads that movement.
#'
#' \eqn{u} is assembled without forming any third-derivative array. Writing
#' \eqn{G_{ab,i} = x_{ia}'M_{[a][b]}x_{ib}} for the per-observation diagonal of
#' the block of \eqn{M},
#' \deqn{u_{k} = -X_k'\Big(w \sum_{a,b} \ell'''_{abk}\, G_{ab}\Big),}
#' one crossprod per distribution parameter. The component
#' \eqn{\ell'''_{abk}} is looked up by a name BUILT from the parameter names in
#' the family's own order, never parsed out of one.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An \code{\link{OuterMethod}}.
#' @param idx The outer index.
#' @param basis The integrated subspace, or \code{NULL}.
#'
#' @return A numeric vector, one entry per row of \code{idx}, or \code{NULL}
#'   where the determinant does not exist.
#'
#' @seealso \code{\link{statmod_marginal}}, \code{\link{reml}}
#'
#' @keywords internal
statmod_marginal_grad <- function(spec, design, coef, hyper, method, idx,
                                  basis = NULL) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)

  H <- statmod_information_at(spec, coef, design, expected = FALSE)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S[!is.finite(S)] <- 0
  K <- H + S
  Kfac <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(Kfac)) return(NULL)
  Kinv <- chol2inv(Kfac)
  M <- if (is.null(basis)) Kinv else {
    A <- basis
    inner <- tryCatch(chol2inv(chol(crossprod(A, K %*% A))),
                      error = function(e) NULL)
    if (is.null(inner)) return(NULL)
    A %*% inner %*% t(A)
  }

  u <- u_vector(spec, design, coef, M, params, npar, offs, total)

  out <- numeric(nrow(idx))
  links <- attr(idx, "links")
  for (a in seq_along(params)) {
    p <- params[a]
    rows <- which(idx$parameter == p)
    if (!length(rows)) next
    for (nm in unique(idx$term[rows])) {
      cols <- design[[p]]$blocks[[nm]]
      pos <- offs[a] + cols
      pen <- modelterms7::term_penalty(spec@terms[[p]][[nm]])
      bt <- coef[[p]][cols]
      th <- as.list(hyper[[p]][[nm]])
      gt <- penalties7::penalty_grad_theta(pen, bt, th)
      cr <- penalties7::penalty_cross(pen, bt, th)
      dS <- penalty_dhessian(pen, bt, th)
      for (r in rows[idx$term[rows] == nm]) {
        h <- idx$name[r]
        # the mode moves by -(H+S)^-1 d2rho/dbeta dtheta
        c_m <- numeric(total)
        c_m[pos] <- as.numeric(cr[[h]])
        v <- -as.numeric(Kinv %*% c_m)
        dS_m <- matrix(0, total, total)
        dS_m[pos, pos] <- dS[[h]]
        dtheta <- -as.numeric(gt[[h]]) -
          (sum(M * dS_m) + sum(u * v)) / 2
        # and onto the free scale the search runs on
        eta <- linkfunctions7::linkfun(links[[r]], hyper[[p]][[nm]][[h]])
        out[r] <- dtheta * linkfunctions7::dlinkinv(links[[r]], eta)
      }
    }
  }
  out
}


#' The Trace of the Determinant's Movement With the Coefficients
#'
#' @description
#' \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)}, assembled one
#' crossprod per distribution parameter.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients.
#' @param M The matrix the trace is taken against.
#' @param params The distribution's parameter names.
#' @param npar,offs,total The block sizes, their offsets and the total.
#'
#' @return A numeric vector as long as the stacked coefficients.
#'
#' @keywords internal
u_vector <- function(spec, design, coef, M, params, npar, offs, total) {
  n <- spec@n_obs
  th <- statmod_eta(spec, design, coef)$theta
  d3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                       scale = "link")
  keys <- names(d3)

  # the per-observation diagonal of each block of M
  G <- vector("list", length(params))
  for (a in seq_along(params)) {
    G[[a]] <- vector("list", length(params))
    if (npar[a] == 0L) next
    Xa <- design[[params[a]]]$X
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      Mab <- M[offs[a] + seq_len(npar[a]), offs[b] + seq_len(npar[b]),
               drop = FALSE]
      G[[a]][[b]] <- rowSums((Xa %*% Mab) * design[[params[b]]]$X)
    }
  }

  out <- numeric(total)
  for (k in seq_along(params)) {
    if (npar[k] == 0L) next
    s <- numeric(n)
    for (a in seq_along(params)) {
      if (npar[a] == 0L) next
      for (b in seq_along(params)) {
        if (npar[b] == 0L) next
        s <- s + rep_len(d3[[d3_key(params, a, b, k, keys)]], n) * G[[a]][[b]]
      }
    }
    out[offs[k] + seq_len(npar[k])] <-
      -as.numeric(crossprod(design[[params[k]]]$X, spec@weights * s))
  }
  out
}


#' The Name of a Third-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, k)} entry of a distribution's third derivative,
#' which is keyed by name and not by position.
#'
#' @details
#' The name is BUILT by putting the three parameter names in the family's own
#' order and joining them, the direction \pkg{distributions7} sanctions, and
#' then checked against the enumeration rather than trusted.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,k Indices into \code{params}.
#' @param keys The names the derivative actually returned.
#'
#' @return A single string.
#'
#' @keywords internal
d3_key <- function(params, a, b, k, keys) {
  want <- paste(params[sort(c(a, b, k))], collapse = "_")
  if (!want %in% keys) {
    stop(sprintf("No third-derivative component '%s'.", want), call. = FALSE)
  }
  want
}
