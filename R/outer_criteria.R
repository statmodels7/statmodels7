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
#' Build the object that tells [statmod()] to choose a model's
#' hyperparameters by an estimate of prediction error. Pass the result as
#' `outer_criterion` for the smooth hyperparameters, as `sparse_criterion`
#' for a lasso's or a SCAD's, or as both.
#'
#' The criterion is
#' \deqn{C(\theta) = -2\ell(\hat\beta(\theta)) + \kappa\,\tau(\theta),}
#' with \eqn{\tau = \mathrm{tr}[(H+S)^{-1}H]} the effective degrees of freedom
#' and \eqn{\kappa} the price of one of them. `aic()` charges 2 by default and
#' `bic()` charges \eqn{\log n}. Both are minimized.
#'
#' @details
#' # The alternative, and when to prefer this
#'
#' [reml()] and [ml()] integrate the coefficients out and maximize what is
#' left, which needs a penalty that is twice differentiable at the mode.
#' A kinked penalty puts the mode on the kink for every coefficient it sets
#' to zero, so a marginal criterion cannot be read there at all. That is the
#' case these two exist for, and `bic()` is `statmod()`'s default
#' `sparse_criterion`.
#'
#' For a smooth penalty both families apply and the choice is the usual one:
#' a marginal criterion asks which hyperparameter makes the data most likely
#' with the coefficients integrated out, and these ask which one predicts
#' best.
#'
#' # Both derivatives are exact
#'
#' By the implicit function theorem, from the same pieces the marginal
#' criteria use. One thing differs and it separates the two families: the
#' envelope theorem does not apply here. \eqn{\ell} alone is not stationary
#' at the penalized mode, so its derivative carries \eqn{d\hat\beta/d\theta}
#' from the first order. What makes that computable is the mode's own
#' condition, \eqn{\partial\ell/\partial\beta = \partial\rho/\partial\beta}.
#'
#' The exact derivatives need `hessian = "observed"`, which is the default.
#'
#' # Which effective degrees of freedom
#'
#' The trace runs over the whole coefficient vector, so a term's contribution
#' reads the full penalized information and not only its own block.
#' [statmod_edf()]'s per-term numbers invert the block instead, which is what
#' a per-term reading has to do. The two agree when the blocks are orthogonal
#' and differ when they are not.
#'
#' Where a penalty has a kink the trace is restricted to the active
#' coordinates, so for a lasso \eqn{\tau} is the number of surviving
#' coefficients, the unbiased count of Zou, Hastie and Tibshirani (2007).
#'
#' # GCV is not offered
#'
#' Classical GCV divides a residual sum of squares by \eqn{(n-\tau)^2}, which
#' estimates an unknown scale. Here every distribution parameter carries its
#' own equation and nothing is profiled, so there is no unknown scale for
#' that ratio to estimate and the criterion degenerates to `aic()`. That is
#' the substitution Wood (2008) makes in the other direction when the scale
#' is known.
#'
#' A GCV on the squared error of the fitted mean is a different and
#' well-defined object. It would need the derivative of that mean in the
#' parameters, which is not one of \pkg{distributions7}'s generics.
#'
#' # Over a kinked penalty these sweep a path
#'
#' The penalized mode is only piecewise smooth in such a hyperparameter,
#' turning a corner whenever a coefficient joins the active set or leaves it,
#' so the criterion inherits the corners and is swept rather than
#' differentiated. How the path covers a term carrying more than one kinked
#' hyperparameter belongs to the term, through
#' [modelterms7::term_search()]: the same criterion object is asked about the
#' model's smooth hyperparameters as well, and those are not swept at all.
#'
#' @param k The price of one degree of freedom, a single non-negative finite
#'   number. Defaults to 2, which is Akaike's. `k = log(n)` by hand is `bic()`.
#'   Anything else is an error. `bic()` does not take it: its \eqn{\kappa} is
#'   `NA` on the object and resolved to \eqn{\log n} when the model is fitted,
#'   the sample size being unknown until then.
#' @param hessian Which information the criterion is built on, `"observed"`
#'   (the default) or `"expected"`. Matched with [match.arg()]. The exact
#'   gradient and Hessian of the criterion need the observed one; with
#'   `"expected"` the outer search is derivative-free.
#'
#' @return An [OuterMethod()] object with properties `kind` (`"aic"` or
#'   `"bic"`), `hessian` as supplied, `k` (the number given, or `NA_real_`
#'   for `bic()`), and the path settings `nfolds`, `rule` and `folds`, which
#'   only [cv()] reads.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of freedom of
#' the lasso. *The Annals of Statistics* **35**(5), 2173--2192.
#'
#' Wood, S. N. (2008). Fast stable direct fitting and smoothness selection for
#' generalized additive models. *Journal of the Royal Statistical Society,
#' Series B* **70**(3), 495--518.
#'
#' @seealso [reml()] and [ml()] for the marginal criteria, [cv()] for
#'   cross-validation, [statmod()] for where these are passed.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
#'
#' fa <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
#'               outer_criterion = aic())
#' fb <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
#'               outer_criterion = bic())
#'
#' # BIC charges log(n) = 5.3 per degree of freedom against AIC's 2, so it
#' # buys a smoother fit: a larger smoothing parameter and fewer edf.
#' c(aic = unlist(fa@hyper), bic = unlist(fb@hyper))
#' c(aic = sum(fa@edf$edf), bic = sum(fb@edf$edf))
#'
#' # The object is a specification and carries no data.
#' aic()
#' aic(k = 4)@k
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
#' Reports whether a criterion is to be made small: `TRUE` for [aic()],
#' [bic()] and [cv()], which estimate prediction error, and `FALSE` for
#' [reml()] and [ml()], which are log-likelihoods and are maximized.
#'
#' @details
#' The outer search always minimizes, so this decides whether the criterion
#' goes in as it stands or with its sign turned. Asking the method, rather
#' than writing the list of kinds into the search, means a criterion added
#' later declares its own orientation in one place.
#'
#' @param method An [OuterMethod()], as [reml()], [ml()], [aic()], [bic()] or
#'   [cv()] returns it. Only its `kind` is read.
#'
#' @return A single logical.
#'
#' @seealso [outer_k()] for the other property of a criterion the search
#'   reads, [aic()] for the criteria this answers `TRUE` for.
#'
#' @keywords internal
outer_minimize <- function(method) {
  method@kind %in% c("aic", "bic", "cv")
}


#' The Price of One Degree of Freedom
#'
#' @description
#' Returns \eqn{\kappa}, the price a prediction-error criterion charges for
#' one effective degree of freedom, resolved against the sample size where
#' the method left it open. [bic()] stores `NA` and gets \eqn{\log n} here;
#' every other criterion returns the `k` it was constructed with.
#'
#' @details
#' The resolution happens at fit time and not at construction because
#' [bic()] takes no data and cannot know \eqn{n}. A criterion object is a
#' specification and is reusable across models of different sizes.
#'
#' @param method An [OuterMethod()]. Its `kind` and `k` are read.
#' @param n The number of observations, a single positive number.
#'
#' @return A single number: \eqn{\log n} for [bic()], and `method@k` for
#'   anything else, which is `NA_real_` for [reml()], [ml()] and [cv()],
#'   none of which charges per degree of freedom.
#'
#' @seealso [aic()] for where `k` is set, [outer_minimize()] for the other
#'   property read at the same point.
#'
#' @keywords internal
outer_k <- function(method, n) {
  if (identical(method@kind, "bic")) return(log(n))
  method@k
}


#' The Effective Degrees of Freedom of a Whole Fit
#'
#' @description
#' Computes \eqn{\tau = \mathrm{tr}[(H+S)^{-1}H]}, the trace of the influence
#' matrix, which is the effective degrees of freedom a prediction-error
#' criterion charges for. An unpenalized fit has \eqn{S = 0} and
#' \eqn{\tau = p}; shrinkage lowers it, and \eqn{\tau} falls to the dimension
#' of the penalty's null space as the hyperparameter grows.
#'
#' @details
#' # Restricting to the active set
#'
#' A penalty with a kink is not twice differentiable everywhere, and the mode
#' of such a fit sits on the kink for every coefficient it sets to zero. It
#' is twice differentiable away from the kink, so the trace is taken over the
#' active coordinates alone, which is the submodel the fit has selected.
#'
#' For a penalty that is linear there, the lasso, \eqn{S_{AA}} vanishes and
#' \eqn{\tau} is exactly the number of coefficients that survived, the
#' unbiased count of Zou, Hastie and Tibshirani (2007). The elastic net keeps
#' its quadratic part and spends slightly less; SCAD and MCP contribute their
#' own curvature and spend slightly more.
#'
#' Taking the trace over the whole vector instead cannot see the selection at
#' all: measured on twenty noise columns it read 14 at every value of
#' \eqn{\lambda}, so a criterion built on it prices a model that selects
#' nothing the same as one that selects everything.
#'
#' # The inverse
#'
#' Through a Cholesky of \eqn{J}, and `NA_real_` when that fails, which is
#' what a caller reads as an unusable point. Both matrices are densified
#' first: the trace needs the full inverse, which is dense whatever \eqn{J}
#' was.
#'
#' @param J The penalized information \eqn{H + S}, a `p x p` symmetric
#'   matrix, dense or sparse.
#' @param H The likelihood's information alone, the same shape as `J`.
#' @param active A logical vector over the stacked coefficients marking the
#'   ones away from a kink, or `NULL` (the default) for all of them. Both
#'   matrices are subset to it before the trace.
#'
#' @return A single number: the trace, `0` when `active` selects nothing, and
#'   `NA_real_` when `J` has no Cholesky factor.
#'
#' @seealso [aic()] for the criterion this feeds,
#'   [statmod_edf()] for the per-term reading, which inverts each block
#'   instead.
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
#' Evaluates \eqn{-2\ell + \kappa\tau} at one set of coefficients and
#' hyperparameters. The coefficients are taken as given and no fitting
#' happens here, so a caller who wants the criterion at the penalized mode
#' has to pass the mode.
#'
#' @details
#' The pieces are assembled once each: the weighted log-likelihood, the
#' penalized information \eqn{H + S} and the information \eqn{H} alone, then
#' [outer_tau()] for the trace and [outer_k()] for \eqn{\kappa}.
#'
#' The penalty's own value is reported alongside but does not enter the
#' criterion. A prediction-error criterion prices the fit's complexity
#' through \eqn{\tau}, not through the size of the penalty.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, refreshed at `coef` if any term needs it.
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter.
#' @param hyper The hyperparameters, per penalized term.
#' @param method An [OuterMethod()] of kind `"aic"` or `"bic"`, read for
#'   `hessian` and for \eqn{\kappa}.
#' @param approx How the expected information is approximated for a family
#'   with no closed form. Read only when `method@hessian` is `"expected"`.
#' @param active Which coefficients are away from a kink, as
#'   [statmod_active()] reports them, or `NULL` where no penalty has one.
#'
#' @return A list of four, or `NULL` when the penalized information has no
#'   Cholesky factor:
#'   \describe{
#'     \item{`value`}{the criterion, \eqn{-2\ell + \kappa\tau}.}
#'     \item{`loglik`}{the weighted log-likelihood at `coef`.}
#'     \item{`penalty`}{the penalties' total value there, reported and not
#'       used in `value`.}
#'     \item{`edf`}{the trace \eqn{\tau}.}
#'   }
#'
#' @seealso [aic()] for the criterion, [outer_tau()] for the trace,
#'   [statmod_marginal()] for the marginal criteria's counterpart.
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
