#' @include outer_hessian.R
NULL

# The prediction-error criteria, and their exact derivatives.
#
# A marginal criterion integrates the coefficients out; a prediction-error one
# estimates how well the fit would predict, and every member is a function of
# two quantities:
#
#   the log-likelihood at the penalized mode,   l(b(t))
#   the effective degrees of freedom,           tau(t) = tr(J^-1 H)
#
# with J = H + S the penalized information and H the likelihood's own. An
# unpenalized model has J = H and tau = p, and every degree of freedom a
# penalty takes away shows up there.
#
# Both are differentiated by the implicit function theorem, exactly as the
# marginal criterion is, and out of the same pieces: b_m = -J^-1 c_m, the
# contractions T and U of the third and fourth derivatives of the objective in
# the coefficients, and the penalty's own theta-derivatives. What differs is
# ONE thing, and it is the thing that makes these criteria different objects:
# the envelope theorem does not apply. l alone is not stationary at the mode --
# the penalized objective is -- so dl/dt carries db/dt from the start, and at
# the mode dl/db equals drho/db, which is what makes it computable.

#' Prediction-Error Criteria for the Hyperparameters
#'
#' @description
#' `aic()` and `bic()` choose the hyperparameters by an estimate of
#' prediction error rather than by a marginal likelihood.
#'
#' @details
#' **The criterion** is
#' \deqn{C(\theta) = -2\ell(\hat\beta(\theta)) + \kappa\,\tau(\theta),}
#' with \eqn{\tau = \mathrm{tr}[(H+S)^{-1}H]} the effective degrees of freedom
#' and \eqn{\kappa} the price of one of them: \eqn{2} for `aic()` and
#' \eqn{\log n} for `bic()`, which is `aic()` with that
#' \eqn{\kappa} and is offered separately because it is what anybody would
#' look for.
#'
#' **Both derivatives are exact**, by the implicit function theorem, and
#' come from the same pieces the marginal criterion uses. One thing differs and
#' it is what separates the two families: the envelope theorem does not apply
#' here. \eqn{\ell} alone is not stationary at the penalized mode, so its
#' derivative carries \eqn{d\hat\beta/d\theta} from the first order, and what
#' makes it computable is that at the mode \eqn{\partial\ell/\partial\beta} is
#' \eqn{\partial\rho/\partial\beta}.
#'
#' \strong{Which \eqn{\tau}.} The trace is taken over the whole coefficient
#' vector, so a term's contribution reads the full penalized information and
#' not only its own block. [statmod_edf()]'s per-term numbers invert
#' the block instead, which is what a per-term reading has to do, and the two
#' agree when the blocks are orthogonal and differ when they are not.
#'
#' **GCV is not among these**, and the reason is the framework rather
#' than the work: classical GCV divides a residual sum of squares by
#' \eqn{(n-\tau)^2}, which estimates an unknown scale. Here every distribution
#' parameter has its own equation and nothing is profiled, so there is no
#' unknown scale for that ratio to estimate, and the criterion it degenerates
#' to is `aic()` -- the substitution Wood (2008) makes in the other
#' direction when the scale is known. A GCV on the squared error of the fitted
#' mean is a different and well-defined object; it needs the derivative of that
#' mean in the parameters, which is not one of \pkg{distributions7}'s generics.
#'
#' **Over a penalty with a kink** these criteria sweep a path rather than
#' differentiate. How that path covers a term carrying more than one such
#' hyperparameter is the TERM's own, [modelterms7::term_search()],
#' since the same criterion is asked of the smooth hyperparameters of the
#' model as well and those are not swept at all.
#'
#' @param k The price of one degree of freedom. Defaults to 2; `bic()`
#'   uses \eqn{\log n}, resolved when the model is fitted.
#' @param hessian Which information is used, `"expected"` or
#'   `"observed"`. The exact derivatives need the observed one.
#'
#' @return An [OuterMethod()].
#'
#' @references
#' Wood, S. N. (2008). Fast stable direct fitting and smoothness selection for
#' generalized additive models. *Journal of the Royal Statistical Society,
#' Series B*, 70(3), 495--518.
#'
#' @seealso [reml()], [statmod()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
#' statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
#'         outer_criterion = aic())
#'
#' @export
aic <- function(k = 2, hessian = c("observed", "expected")) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 0) {
    stop("'k' must be a single non-negative number.", call. = FALSE)
  }
  do.call(OuterMethod, c(list(kind = "aic", hessian = match.arg(hessian),
                             k = as.numeric(k)), outer_path_defaults()))
}

#' @rdname aic
#' @export
bic <- function(hessian = c("observed", "expected")) {
  do.call(OuterMethod, c(list(kind = "bic", hessian = match.arg(hessian),
                             k = NA_real_), outer_path_defaults()))
}


#' Is a Criterion Minimized?
#'
#' @description
#' `TRUE` for a prediction-error criterion, `FALSE` for a marginal
#' likelihood, which is maximized.
#'
#' @details
#' The search always minimizes; this is what says whether the criterion goes in
#' as it is or with its sign turned. It is asked of the method rather than
#' written into the search, so a criterion added later declares its own
#' orientation.
#'
#' @param method An [OuterMethod()].
#'
#' @return A single logical.
#'
#' @keywords internal
outer_minimize <- function(method) {
  method@kind %in% c("aic", "bic", "cv")
}


#' The Price of One Degree of Freedom
#'
#' @description
#' \eqn{\kappa}, resolved against the sample size where the method left it
#' open.
#'
#' @param method An [OuterMethod()].
#' @param n The number of observations.
#'
#' @return A single number.
#'
#' @keywords internal
outer_k <- function(method, n) {
  if (identical(method@kind, "bic")) return(log(n))
  method@k
}


#' The Effective Degrees of Freedom of a Whole Fit
#'
#' @description
#' \eqn{\tau = \mathrm{tr}[(H+S)^{-1}H]}, the trace of the influence matrix.
#'
#' A penalty with a kink is not twice differentiable everywhere, and the mode
#' of such a fit sits at the kink for every coefficient it sets to zero. It is
#' twice differentiable away from the kink, so the trace is taken over the
#' active coordinates alone, which is the submodel the fit has selected. For a
#' penalty that is linear there -- the lasso -- \eqn{S_{AA}} vanishes and
#' \eqn{\tau} is the number of coefficients that survived, which is the
#' unbiased count of Zou, Hastie and Tibshirani (2007). The elastic net keeps
#' its quadratic part, and SCAD and MCP their own curvature.
#'
#' @param J The penalized information.
#' @param H The likelihood's information.
#' @param active A logical vector over the stacked coefficients, or
#'   `NULL` for all of them.
#'
#' @return A single number.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of freedom of
#' the lasso. *The Annals of Statistics* 35(5), 2173--2192.
#'
#' @keywords internal
outer_tau <- function(J, H, active = NULL) {
  J <- as_dense(J)
  H <- as_dense(H)
  if (!is.null(active)) {
    if (!any(active)) return(0)
    J <- J[active, active, drop = FALSE]
    H <- H[active, active, drop = FALSE]
  }
  P <- tryCatch(chol2inv(chol(J)), error = function(e) NULL)
  if (is.null(P)) return(NA_real_)
  sum(P * t(H))
}


#' A Prediction-Error Criterion at Given Coefficients and Hyperparameters
#'
#' @description
#' \eqn{-2\ell + \kappa\tau} at the penalized mode.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param approx The approximation for the expected information.
#' @param active Which coefficients are away from a kink, as
#'   [statmod_active()] reports them, or `NULL` where no
#'   penalty has one.
#'
#' @return A list with `value`, `loglik`, `penalty` and
#'   `edf`, or `NULL` where the information is not invertible.
#'
#' @seealso [aic()]
#'
#' @keywords internal
statmod_pe <- function(spec, design, coef, hyper, method,
                       approx = "bartlett", active = NULL) {
  expected <- identical(method@hessian, "expected")
  ll <- statmod_loglik_at(spec, coef, design)
  H <- statmod_information_at(spec, coef, design, expected, approx)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S <- zap_nonfinite(S)
  tau <- outer_tau(H + S, H, active)
  if (!is.finite(tau)) return(NULL)
  k <- outer_k(method, spec@n_obs)
  list(value = -2 * ll + k * tau, loglik = ll,
       penalty = statmod_penalty_at(spec, coef, hyper, design, "value"),
       edf = tau)
}


#' The Exact Derivatives of a Prediction-Error Criterion
#'
#' @description
#' The gradient and, when asked, the Hessian of \eqn{-2\ell + \kappa\tau} over
#' the free scale of the hyperparameters.
#'
#' @details
#' Writing \eqn{P = J^{-1}}, \eqn{A_m = \partial J/\partial\theta_m = S_m +
#' T[\hat\beta_m]} and \eqn{B_m = \partial H/\partial\theta_m = T[\hat\beta_m]},
#' \deqn{\tau_m = \mathrm{tr}(PB_m) - \mathrm{tr}(PA_mPH),}
#' and the log-likelihood contributes \eqn{(\partial\rho/\partial\beta)'
#' \hat\beta_m}. Differentiating once more brings in \eqn{A_{ml}}, \eqn{B_{ml}}
#' and \eqn{\hat\beta_{ml}}, which are the quantities
#' [statmod_marginal_hess()] already assembles.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param order `1` for the gradient alone, `2` for both.
#'
#' @return A list with `grad` and, at order 2, `hess`; or
#'   `NULL` where the information is not invertible.
#'
#' @seealso [aic()], [statmod_marginal_grad()]
#'
#' @keywords internal
statmod_pe_derivs <- function(spec, design, coef, hyper, method, idx,
                              order = 1L) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  nh <- nrow(idx)
  n <- spec@n_obs
  kap <- outer_k(method, n)

  H <- statmod_information_at(spec, coef, design, expected = FALSE)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S <- zap_nonfinite(S)
  # The routes below form a FULL inverse, and the inverse of a sparse matrix
  # is dense, so densifying here gives up nothing sparsity could have kept.
  J <- as_dense(H + S)
  Jfac <- tryCatch(chol(J), error = function(e) NULL)
  if (is.null(Jfac)) return(NULL)
  P <- chol2inv(Jfac)

  th <- statmod_eta(spec, design, coef)$theta
  d3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                       scale = "link", threads = spec@threads)
  d4 <- if (order >= 2L) {
    distributions7::distrib_deriv4(spec@distrib, spec@response, th,
                                   scale = "link", threads = spec@threads)
  } else NULL

  pieces <- outer_pieces(spec, design, coef, hyper, idx, offs, total, order)
  grho <- statmod_penalty_at(spec, coef, hyper, design, "gradient")
  grho <- unlist(grho[params], use.names = FALSE)

  bhat <- lapply(seq_len(nh), function(m) -as.numeric(P %*% pieces$c[[m]]))
  tv <- lapply(bhat, function(v) block_predictors(design, params, npar, offs,
                                                  v))
  Bm <- lapply(tv, function(t) contract3(spec, design, d3, params, npar, offs,
                                         total, t))
  Am <- lapply(seq_len(nh), function(m) pieces$S[[m]] + Bm[[m]])
  PH <- P %*% H

  # one per hyperparameter, and the gradient and the Hessian both read them
  PA <- lapply(Am, function(A) P %*% A)

  g <- numeric(nh)
  for (m in seq_len(nh)) {
    tau_m <- sum(P * t(Bm[[m]])) - sum(PA[[m]] * t(PH))
    g[m] <- -2 * sum(grho * bhat[[m]]) + kap * tau_m
  }
  if (order < 2L) return(list(grad = free_scale(g, hyper, idx)))

  # These depend on ONE hyperparameter and were formed inside the PAIR loop,
  # so each dense p by p product was taken nh + 1 times over: measured, 6
  # against 2 at two hyperparameters and 12 against 3 at three. It is the
  # hoist statmod_marginal_hess received in 0.49.0 -- "the two expensive
  # quantities each depend on ONE direction and not on a pair, so they are
  # built once per hyperparameter" -- which the prediction-error route never
  # got, that release having deliberately left aic(), bic() and cv() alone.
  # Nothing about the arithmetic changes: the same product is taken once
  # instead of many times, so the answer is identical to the bit.
  PAPH <- lapply(PA, function(M) M %*% PH)
  PB <- lapply(Bm, function(B) P %*% B)
  Hb <- lapply(bhat, function(v) as.numeric(H %*% v))

  Hm <- matrix(0, nh, nh)
  for (m in seq_len(nh)) {
    for (l in m:nh) {
      key <- pieces$pair[m, l]
      Bml <- contract4(spec, design, d4, params, npar, offs, total, tv[[l]],
                       tv[[m]])
      rhs <- Am[[l]] %*% bhat[[m]] + pieces$S[[m]] %*% bhat[[l]] +
        pieces$c2[[key]]
      bml <- -as.numeric(P %*% rhs)
      Tml <- contract3(spec, design, d3, params, npar, offs, total,
                       block_predictors(design, params, npar, offs, bml))
      Bml <- Bml + Tml
      Aml <- pieces$S2[[key]] + Bml

      # -2 l(b(t)) : the second derivative carries the curvature of the
      # log-likelihood along the two directions and the mode's own second
      # movement, the score at the mode being drho/dbeta
      ll_ml <- -sum(bhat[[l]] * Hb[[m]]) + sum(grho * bml)

      tau_ml <- sum(PA[[l]] * t(PAPH[[m]])) - sum(P * t(Aml %*% PH)) +
        sum(PA[[m]] * t(PAPH[[l]])) - sum(PA[[m]] * t(PB[[l]])) -
        sum(PA[[l]] * t(PB[[m]])) + sum(P * t(Bml))

      v <- -2 * ll_ml + kap * tau_ml
      Hm[m, l] <- v
      Hm[l, m] <- v
    }
  }
  list(grad = free_scale(g, hyper, idx),
       hess = free_scale2(Hm, g, hyper, idx))
}


#' Carry a Hyperparameter Derivative onto the Free Scale
#'
#' @description
#' The diagonal chain rule each hyperparameter's own link induces, at first
#' order for a gradient and at second order for a Hessian.
#'
#' @param g A gradient on the parameter scale.
#' @param M A Hessian on the parameter scale.
#' @param hyper The hyperparameters.
#' @param idx The outer index.
#'
#' @return A vector, or a matrix.
#'
#' @keywords internal
free_scale <- function(g, hyper, idx) {
  g * link_slopes(hyper, idx)$h1
}

#' @rdname free_scale
#' @keywords internal
free_scale2 <- function(M, g, hyper, idx) {
  s <- link_slopes(hyper, idx)
  out <- M * outer(s$h1, s$h1)
  diag(out) <- diag(out) + s$h2 * g
  out
}

#' @rdname free_scale
#' @keywords internal
link_slopes <- function(hyper, idx) {
  links <- attr(idx, "links")
  nh <- nrow(idx)
  h1 <- numeric(nh)
  h2 <- numeric(nh)
  for (r in seq_len(nh)) {
    v <- hyper[[idx$parameter[r]]][[idx$term[r]]][[idx$name[r]]]
    e <- linkfunctions7::linkfun(links[[r]], v)
    h1[r] <- linkfunctions7::dlinkinv(links[[r]], e)
    h2[r] <- linkfunctions7::d2linkinv(links[[r]], e)
  }
  list(h1 = h1, h2 = h2)
}
