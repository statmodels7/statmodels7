#' @include spec.R
NULL

# The objective, its score and its information.
#
# eta_k = X_k beta_k + offset_k is the linear predictor of the k-th
# distribution parameter, on the link scale; theta_k = g_k^-1(eta_k) is what
# the distribution's generics take. Those generics return derivatives with
# respect to the LINK scale when asked, so the chain rule through the link is
# already done and what is left is one crossprod per block:
#
#   dl/dbeta_k          = X_k' (w * g_k)
#   d2l/dbeta_k dbeta_l = X_k' diag(w * W_kl) X_l
#
# Measured on 20000 x 221, one iteration of this spends 98.1% of its time
# inside BLAS, so there is nothing here to compile.

#' The Linear Predictors and the Parameters They Give
#'
#' @description
#' Turns a coefficient structure into the per-observation parameters every
#' \pkg{distributions7} generic takes. Each equation's design, coefficients,
#' offset and adjustment give a linear predictor \eqn{\eta_k}, and the
#' parameter's own link inverts it to \eqn{\theta_k}.
#'
#' @details
#' \deqn{\eta_k = X_k \beta_k + o_k + a_k, \qquad \theta_k = g_k^{-1}(\eta_k)}
#'
#' with \eqn{o_k} the equation's offset and \eqn{a_k} the adjustment a
#' refreshable term contributes, zero for an ordinary design.
#'
#' The inverse link clamps \eqn{\theta_k} strictly inside its own open
#' bounds, so a predictor far out on the link scale gives a parameter a
#' density can be evaluated at. A structural term of the filter shape adds
#' its level to \eqn{\eta_k} before the link is inverted.
#'
#' @param spec A [StatmodSpec()], read for the distribution, its links, the
#'   offsets and the sample size.
#' @param design The design, as [statmod_design()] returns it, refreshed at
#'   `coef` if any term needs it.
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter, each as long as its equation's design is wide.
#'
#' @return A list of two named lists, each with one entry per distribution
#'   parameter in the family's order:
#'   \describe{
#'     \item{`eta`}{the linear predictors, numeric vectors of length
#'       `spec@n_obs`.}
#'     \item{`theta`}{the parameters they imply, the same shape, each
#'       strictly inside its parameter's bounds.}
#'   }
#'
#' @seealso [statmod_loglik_at()] and [statmod_score_at()], which read this,
#'   [statmod_design_at()] for the refresh.
#'
#' @keywords internal
statmod_eta <- function(spec, design, coef) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  sst <- statmod_structural_state(design)
  if (!is.null(sst)) {
    # a filter is the expensive part of every evaluation and the objective,
    # its gradient and its curvature are asked for at the same point in turn
    key <- list(unlist(coef, use.names = FALSE), sst$zeta)
    if (!is.null(sst$key) && identical(sst$key, key)) return(sst$value)
  }
  design <- statmod_design_at(spec, coef, design)
  eta <- stats::setNames(vector("list", length(params)), params)
  theta <- eta
  for (p in params) {
    d <- design[[p]]
    e <- if (d$npar == 0L) rep(0, spec@n_obs) else
      as.numeric(d$X %*% coef[[p]])
    # what a term contributes is not always its block times its coefficients:
    # where the block is a Jacobian the two differ, and the difference is
    # carried here so that every crossprod elsewhere still reads the block
    if (!is.null(d$adj)) e <- e + d$adj
    off <- spec@offsets[[p]]
    if (!is.null(off)) e <- e + off
    eta[[p]] <- e
    theta[[p]] <- linkfunctions7::linkinv(links[[p]], e)
  }
  # a structural term rewrites the predictor of the equation it sits in: the
  # static part above is what its filter is handed, and what comes back is
  # the predictor with the term's own level in it
  filters <- list()
  regimes <- list()
  if (!is.null(sst)) {
    eta_static <- eta
    filters <- statmod_filter_at(spec, design, eta, theta)
    for (f in filters) {
      eta[[f$param]] <- f$eta
      theta[[f$param]] <- linkfunctions7::linkinv(links[[f$param]], f$eta)
    }
    # a term whose contribution is a likelihood mixed over states reports no
    # predictor of its own; what is recorded here is the posterior-weighted
    # one, which is what a fitted value means for it
    regimes <- statmod_regime_at(spec, design, eta_static, theta)
    for (r in regimes) {
      eta[[r$param]] <- r$eta_bar
      theta[[r$param]] <- linkfunctions7::linkinv(links[[r$param]], r$eta_bar)
    }
    out <- list(eta = eta, theta = theta, filters = filters,
                regimes = regimes, eta_static = eta_static)
    sst$key <- list(unlist(coef, use.names = FALSE), sst$zeta)
    sst$value <- out
    return(out)
  }
  list(eta = eta, theta = theta, filters = filters, regimes = regimes)
}


#' The Weighted Log-Likelihood of a Specification at Given Coefficients
#'
#' @description
#' Computes \deqn{\ell(\beta) = \sum_i w_i \log f(y_i; \theta_i),} the
#' weighted log-likelihood of the whole model at one set of coefficients. The
#' weights enter as supplied and are never normalized, so a weight of two
#' counts an observation twice.
#'
#' @details
#' No penalty enters. This is the likelihood alone; [statmod_penalty_at()]
#' gives the other half of the objective, and [statmod_objective()] combines
#' them.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter.
#' @param design The design. Rebuilt from `spec` when absent, which is
#'   convenient at the console and wasteful in a loop.
#'
#' @return A single number. `-Inf` where the density is zero at some
#'   observation, and `NaN` where the parameters are outside the family's
#'   support, both of which a search reads as an unusable point.
#'
#' @seealso [statmod_score_at()] for its derivative,
#'   [statmod_objective()] for the penalized objective it enters.
#'
#' @keywords internal
statmod_loglik_at <- function(spec, coef, design = statmod_design(spec)) {
  ev <- statmod_eta(spec, design, coef)
  # a latent Markov term replaces the likelihood rather than shifting a
  # predictor, so the contribution is the term's own and not the density at
  # any single point
  if (length(ev$regimes)) {
    r <- ev$regimes[[1L]]
    ll <- modelterms7::term_loglik(r$tm, r$eta_static, spec@response,
                                   r$cb$logdens, r$cb$score, r$psi)$loglik
    return(sum(spec@weights * ll))
  }
  ll <- distributions7::distrib_pdf(spec@distrib, spec@response, ev$theta,
                                    log = TRUE, threads = spec@threads)
  sum(spec@weights * ll)
}


#' The Score of the Weighted Log-Likelihood
#'
#' @description
#' Computes the score of the weighted log-likelihood in the coefficients, one
#' block per distribution parameter:
#'
#' \deqn{\frac{\partial \ell}{\partial \beta_k} = X_k'\,(w \odot g_k)}
#'
#' with \eqn{g_k} the per-observation derivative of the log-density in the
#' link-scale predictor \eqn{\eta_k}, which \pkg{distributions7} supplies
#' already chained onto that scale.
#'
#' @details
#' This is the score of the likelihood alone, with no penalty in it, and it
#' is the gradient of a model whose design is fixed.
#'
#' A model carrying a structural term needs a correction that is not here:
#' the term's level is driven by scores read at earlier predictors, so a
#' coefficient reaches it through the recursion. [statmod_structural_score()]
#' and `modelterms7::term_adjoint()` supply that.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design, refreshed at `coef` if any term needs it.
#'
#' @return A named list of numeric vectors, one per distribution parameter in
#'   the family's order, each as long as that equation's design is wide.
#'
#' @seealso [statmod_loglik_at()] for the value,
#'   [statmod_information_at()] for the curvature,
#'   [statmod_structural_score()] for a filter's own parameters.
#'
#' @keywords internal
statmod_score_at <- function(spec, coef, design = statmod_design(spec)) {
  params <- spec@distrib@params
  ev <- statmod_eta(spec, design, coef)
  design <- statmod_design_at(spec, coef, design)
  n <- spec@n_obs

  # Fisher's identity: the derivative of a likelihood mixed over states is
  # the posterior-weighted derivative of the ordinary one, in EVERY
  # predictor and not only the one the regimes shift. K vectorized passes
  # replace the per-observation callback a filter needs.
  if (length(ev$regimes)) {
    r <- ev$regimes[[1L]]
    gv <- stats::setNames(lapply(params, function(p) numeric(n)), params)
    for (k in seq_len(component_count(r$mu))) {
      thk <- statmod_theta_shifted(spec, ev$eta_static, r$param,
                                   component_shift(r$mu, k))
      gk <- distributions7::distrib_gradient(spec@distrib, spec@response, thk,
                                             scale = "link", threads = spec@threads)
      for (p in params) {
        gv[[p]] <- gv[[p]] + r$gamma[, k] * spec@weights * rep_len(gk[[p]], n)
      }
    }
    return(stats::setNames(lapply(params, function(p) {
      d <- design[[p]]
      if (d$npar == 0L) return(numeric(0))
      as.numeric(crossprod(d$X, gv[[p]]))
    }), params))
  }

  g <- distributions7::distrib_gradient(spec@distrib, spec@response, ev$theta,
                                        scale = "link", threads = spec@threads)
  gv <- stats::setNames(lapply(params, function(p)
    spec@weights * rep_len(g[[p]], n)), params)

  # Where a filter drives one equation, the derivative of the objective in a
  # coefficient of ANY equation carries a term the block does not: the level
  # was driven by scores read at predictors those coefficients enter. The
  # reverse recursion returns the derivative in that score sequence, and the
  # derivative of the score in another equation's predictor is the mixed
  # second derivative of the log-density. For the filter's own equation this
  # reproduces its `deta`, the curvature being the second derivative there.
  for (f in ev$filters) {
    ad <- modelterms7::term_adjoint(f$tm, f$eta_static, spec@response,
                                    f$cb$score, f$cb$curvature, f$psi,
                                    g = gv[[f$param]],
                                    fast = f$cb$fast, threads = spec@threads)
    H <- distributions7::distrib_hessian(spec@distrib, spec@response,
                                         ev$theta, scale = "link", threads = spec@threads)
    a <- match(f$param, params)
    add <- stats::setNames(lapply(params, function(q) {
      ad$dscore * rep_len(H[[hess_key(params, a, match(q, params))]], n)
    }), params)
    for (q in params) gv[[q]] <- gv[[q]] + add[[q]]
  }

  stats::setNames(lapply(params, function(p) {
    d <- design[[p]]
    if (d$npar == 0L) return(numeric(0))
    as.numeric(crossprod(d$X, gv[[p]]))
  }), params)
}


#' The Score in a Structural Term's Own Parameters
#'
#' @description
#' The derivative of the weighted log-likelihood in the unconstrained
#' parameters of each structural term, one entry per term.
#'
#' @details
#' The filter's Jacobian is already the **total** derivative of the predictor
#' in the term's own parameters, with the recursion propagated inside it. So
#' a chain rule over the observations is the whole of the work here and no
#' reverse pass is needed.
#'
#' The reverse pass answers a different question: the derivative in the
#' coefficients of the equations, where a coefficient reaches the level only
#' through the scores at earlier times. That is `modelterms7::term_adjoint()`.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design, whose structural state holds each term's current
#'   parameters.
#'
#' @return A named list of numeric vectors, one entry per structural term,
#'   keyed by the term's key in the specification, each as long as that term
#'   has parameters. An empty list when the model carries no structural term.
#'
#' @seealso [statmod_filter_at()] for the filter run this reads,
#'   [statmod_score_at()] for the coefficients' own score.
#'
#' @keywords internal
statmod_structural_score <- function(spec, coef,
                                     design = statmod_design(spec)) {
  ev <- statmod_eta(spec, design, coef)
  sst <- statmod_structural_state(design)
  chain_of <- function(tm, key) {
    links <- modelterms7::term_links(tm)
    z <- sst$zeta[[key]]
    vapply(names(z), function(j)
      linkfunctions7::dlinkinv(links[[j]], z[[j]]), numeric(1))
  }
  if (length(ev$regimes)) {
    out <- list()
    for (r in ev$regimes) {
      ll <- modelterms7::term_loglik(r$tm, r$eta_static, spec@response,
                                     r$cb$logdens, r$cb$score, r$psi)
      dpsi <- as.numeric(crossprod(ll$jacobian, spec@weights))
      out[[r$term]] <- stats::setNames(dpsi * chain_of(r$tm, r$term),
                                       names(sst$zeta[[r$term]]))
    }
    return(out)
  }
  if (!length(ev$filters)) return(list())
  g <- distributions7::distrib_gradient(spec@distrib, spec@response, ev$theta,
                                        scale = "link", threads = spec@threads)
  n <- spec@n_obs
  out <- list()
  for (f in ev$filters) {
    gp <- spec@weights * rep_len(g[[f$param]], n)
    dpsi <- as.numeric(crossprod(f$jacobian, gp))
    links <- modelterms7::term_links(f$tm)
    z <- sst$zeta[[f$term]]
    chain <- vapply(names(z), function(j)
      linkfunctions7::dlinkinv(links[[j]], z[[j]]), numeric(1))
    out[[f$term]] <- stats::setNames(dpsi * chain, names(z))
  }
  out
}


#' The Information of the Weighted Log-Likelihood
#'
#' @description
#' The negative Hessian in the coefficients, assembled block by block from the
#' distribution's own link-scale second derivatives.
#'
#' @details
#' Block \eqn{(a, b)} is \eqn{X_a'\,\mathrm{diag}(w\,h_{ab})\,X_b} with
#' \eqn{h_{ab}} the per-observation second derivative of the log-density in
#' \eqn{(\eta_a, \eta_b)}, negated. The blocks are assembled in the storage
#' the designs call for, so a model with a sparse equation gets a sparse
#' information.
#'
#' `expected = TRUE` gives the expected information, which Fisher scoring
#' inverts and which is positive definite wherever the family is regular.
#' `expected = FALSE` gives the negated observed Hessian, which far from the
#' optimum may be indefinite.
#'
#' `approx` reaches \pkg{distributions7} and is read only where the family
#' has no closed expected information.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design, refreshed at `coef` if any term needs it.
#' @param expected `TRUE` for the expected information, `FALSE` for the
#'   negated observed Hessian.
#' @param approx How the expected information is approximated for a family
#'   with no closed form: `"bartlett"`, `"integrate"` or `"mc"`.
#'
#' @return A symmetric `p x p` matrix over the stacked coefficients, `p`
#'   being their total count across the equations. A \pkg{Matrix} object when
#'   any equation's design is sparse, a base matrix otherwise.
#'
#' @seealso [statmod_score_at()] for the first derivative,
#'   [info_blocks()] for the per-observation blocks a square-root route uses
#'   instead.
#'
#' @keywords internal
statmod_information_at <- function(spec, coef, design = statmod_design(spec),
                                   expected = TRUE, approx = "bartlett") {
  params <- spec@distrib@params
  ev <- statmod_eta(spec, design, coef)
  th <- ev$theta
  design <- statmod_design_at(spec, coef, design)

  # For a likelihood mixed over states the matrix assembled here is the
  # COMPLETE-DATA information, the ordinary one averaged over the smoothed
  # states. It is the matrix an EM step inverts and it is positive definite,
  # which is all a scoring matrix has to be once the gradient is exact. It
  # is not the observed information of the mixture: by the
  # missing-information principle the two differ by the conditional variance
  # of the complete-data score, so this one is the LARGER and a standard
  # error read off it is too small. What vcov() reports comes from
  # statmod_regime_information(), which differentiates the forward recursion
  # twice.
  if (length(ev$regimes)) {
    r <- ev$regimes[[1L]]
    npar <- vapply(design, function(d) d$npar, integer(1))
    offs <- cumsum(npar) - npar
    total <- sum(npar)
    out <- zero_information(design, total)
    for (k in seq_len(component_count(r$mu))) {
      thk <- statmod_theta_shifted(spec, ev$eta_static, r$param,
                                   component_shift(r$mu, k))
      Hk <- if (expected) {
        distributions7::distrib_expected_hessian(spec@distrib, spec@response,
                                                 thk, scale = "link",
                                                 approx = approx, threads = spec@threads)
      } else {
        distributions7::distrib_hessian(spec@distrib, spec@response, thk,
                                        scale = "link", threads = spec@threads)
      }
      for (a in seq_along(params)) {
        if (npar[a] == 0L) next
        for (b in seq_along(params)) {
          if (npar[b] == 0L || b < a) next
          wv <- spec@weights * r$gamma[, k] *
            rep_len(Hk[[hess_key(params, a, b)]], spec@n_obs)
          blk <- wcrossprod(design[[params[a]]]$X, wv,
                        design[[params[b]]]$X, spec@threads)
          ra <- offs[a] + seq_len(npar[a])
          rb <- offs[b] + seq_len(npar[b])
          out[ra, rb] <- out[ra, rb] + blk
          if (a != b) out[rb, ra] <- out[rb, ra] + t(blk)
        }
      }
    }
    return(-out)
  }

  H <- if (expected) {
    distributions7::distrib_expected_hessian(spec@distrib, spec@response, th,
                                             scale = "link", approx = approx, threads = spec@threads)
  } else {
    distributions7::distrib_hessian(spec@distrib, spec@response, th,
                                    scale = "link", threads = spec@threads)
  }
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  out <- zero_information(design, total)
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in seq_along(params)) {
      if (npar[b] == 0L || b < a) next
      key <- hess_key(params, a, b)
      wv <- spec@weights * rep_len(H[[key]], spec@n_obs)
      blk <- wcrossprod(design[[params[a]]]$X, wv,
                        design[[params[b]]]$X, spec@threads)
      ra <- offs[a] + seq_len(npar[a])
      rb <- offs[b] + seq_len(npar[b])
      out[ra, rb] <- blk
      if (a != b) out[rb, ra] <- t(blk)
    }
  }
  # the information is minus the Hessian
  -out
}


#' The Name of a Second-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b)} entry of a distribution's Hessian list, which is
#' keyed by [distributions7::hess_names()] and never by position.
#'
#' @details
#' The key is built by pasting the two parameter names in the order
#' [distributions7::hess_names()] uses, and is never parsed back out of a
#' name. That is the discipline \pkg{distributions7} records for a family
#' whose parameter name contains the separator: splitting `"log_scale_mu"` on
#' the underscore gives the wrong pair, while generating the name from the
#' same enumeration that generated the list cannot.
#'
#' @param params The parameter names, in the family's order, as
#'   `distrib@params` gives them.
#' @param a,b Indices into `params`, in either order: the Hessian is
#'   symmetric and this returns the key the list actually holds.
#'
#' @return A single string, one of the names of the list
#'   [distributions7::distrib_hessian()] returns.
#'
#' @seealso [statmod_information_at()], its caller.
#'
#' @keywords internal
hess_key <- function(params, a, b) {
  nm <- distributions7::hess_names(params)
  want <- paste(params[sort(c(a, b))], collapse = "_")
  if (want %in% nm) return(want)
  alt <- paste(params[rev(sort(c(a, b)))], collapse = "_")
  if (alt %in% nm) return(alt)
  stop(sprintf("No Hessian component for '%s' and '%s'.",
               params[a], params[b]), call. = FALSE)
}


#' The Penalty of a Specification at Given Coefficients
#'
#' @description
#' The sum of every penalized term's `penalty_value()`, with its gradient
#' and Hessian placed in the columns that term owns.
#'
#' @details
#' Each penalized term is asked for the quantity at its own coefficients and
#' its own hyperparameters, and the answer is placed in the columns that term
#' owns. A term declaring several penalties, each over a subset of its
#' coefficients, contributes each one at its own coordinates.
#'
#' The accumulator is a base matrix even where a penalty answers with a
#' sparse Hessian, and the coercion is written once at that boundary. Making
#' it the design's own kind was measured at 0.8x end to end against a blast
#' radius of twenty-two consumers.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param hyper A named list with one entry per distribution parameter, each
#'   a named list of hyperparameter vectors keyed by term.
#' @param design The design.
#' @param what Which quantity: `"value"`, `"gradient"` or `"hessian"`.
#'
#' @return Depends on `what`:
#'   \describe{
#'     \item{`"value"`}{a single number, the penalties summed. `0` for an
#'       unpenalized model.}
#'     \item{`"gradient"`}{a named list of numeric vectors, one per
#'       distribution parameter, zero in every unpenalized coordinate.}
#'     \item{`"hessian"`}{a symmetric `p x p` base matrix over the stacked
#'       coefficients, zero outside the penalized blocks.}
#'   }
#'
#' @seealso [statmod_loglik_at()] for the other half of the objective,
#'   [statmod_structural_penalty()] for penalties over a structural term's
#'   own parameters, which this does not cover.
#'
#' @keywords internal
statmod_penalty_at <- function(spec, coef, hyper,
                               design = statmod_design(spec),
                               what = c("value", "gradient", "hessian")) {
  what <- match.arg(what)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)

  value <- 0
  grad <- stats::setNames(lapply(npar, numeric), params)
  # The accumulator follows the DESIGN, by zero_information()'s own rule, so
  # the penalty's Hessian is stored the way the information it is added to
  # already is.
  #
  # It was dense until 0.53.0, on the argument that eighteen places read this
  # result and only the two a sparse design exercises would be caught by the
  # suite. What settled it is a number the argument did not have: a random
  # intercept over 1000 levels is 4.030 s against 2.943 s, 1.37x, and over 500
  # 1.11x, with the matrix identical to the last bit and 2.02 MB becoming
  # 0.010 MB. And the win is not this allocation -- that is 0.229 ms of a
  # 0.992 ms call -- it is that a sparse S stays sparse through
  # penalty_sqrt(), whose factor goes from 2.0 MB to 0.010 MB, and through
  # augmented_solve(), which no longer converts it.
  #
  # What the conversion needed at the consumers is one function:
  # zap_nonfinite(), since `S[!is.finite(S)] <- 0` builds a dense logical
  # index and would have thrown the storage away at the first reader.
  hess <- zero_information(design, total)

  # A penalty over a structural term's own parameters is evaluated from the
  # term's state, not from the coefficients: it enters the VALUE, which is one
  # number and belongs to the whole objective, while its derivatives are in a
  # different vector and are returned by statmod_structural_penalty().
  sst <- statmod_structural_state(design)
  for (u in statmod_penalized(spec, design)) {
    p <- u$param
    if (isTRUE(u$structural)) {
      if (what != "value" || is.null(sst)) next
      z <- sst$zeta[[u$term]][u$cols]
      value <- value + penalties7::penalty_value(u$penalty, as.numeric(z),
                                                 as.list(hyper[[p]][[u$key]]))
      next
    }
    b <- coef[[p]][u$cols]
    th <- as.list(hyper[[p]][[u$key]])
    if (what == "value") {
      value <- value + penalties7::penalty_value(u$penalty, b, th)
    } else if (what == "gradient") {
      grad[[p]][u$cols] <- grad[[p]][u$cols] +
        penalties7::penalty_gradient(u$penalty, b, th)
    } else {
      # the one point the two kinds meet. A penalty that avoids assembling
      # its own matrix -- a smooth repeated over the levels of a factor --
      # answers sparse, and a dense accumulator cannot take that block, so
      # the coercion is written once here rather than at every reader.
      blk <- penalties7::penalty_hessian(u$penalty, b, th)
      hess[u$index, u$index] <- hess[u$index, u$index] +
        (if (isS4(hess)) blk else as_dense(blk))
    }
  }
  switch(what, value = value, gradient = grad, hessian = hess)
}


#' The Penalty Over a Structural Term's Own Parameters
#'
#' @description
#' The derivative of the penalties a structural term declares, in the term's
#' own parameters, in place of the coefficients the other penalty function
#' reads.
#'
#' @details
#' The objective of the structural block includes these penalties through
#' [statmod_penalty_at()], so its gradient must include their
#' derivative: without it the two describe different functions, an optimizer
#' walks until its budget runs out, and `optimizers7`'s own check reports
#' that the objective changes at a rate the gradient does not predict.
#'
#' The penalty is read on the **unconstrained scale**, which is where the
#' term carries its parameters and where a deviation from a population value
#' is defined. For a deviation, whose link is the identity, the two scales
#' coincide.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, whose structural state holds each term's current
#'   parameters.
#' @param hyper The hyperparameters, per penalized term.
#' @param what Which quantity: `"value"`, `"gradient"` or `"hessian"`.
#'
#' @return A named list, one entry per structural term, each a numeric vector
#'   or matrix over that term's parameters in their own order; empty when no
#'   structural term carries a penalty.
#'
#' @seealso [statmod_penalty_at()]
#'
#' @keywords internal
statmod_structural_penalty <- function(spec, design, hyper,
                                       what = c("value", "gradient",
                                                "hessian")) {
  what <- match.arg(what)
  sst <- statmod_structural_state(design)
  if (is.null(sst)) return(list())
  out <- list()
  for (u in statmod_penalized(spec, design)) {
    if (!isTRUE(u$structural)) next
    z <- sst$zeta[[u$term]]
    b <- as.numeric(z[u$cols])
    th <- as.list(hyper[[u$param]][[u$key]])
    if (is.null(out[[u$term]])) {
      out[[u$term]] <- switch(
        what,
        value = 0,
        gradient = stats::setNames(numeric(length(z)), names(z)),
        hessian = matrix(0, length(z), length(z),
                         dimnames = list(names(z), names(z))))
    }
    if (what == "value") {
      out[[u$term]] <- out[[u$term]] +
        penalties7::penalty_value(u$penalty, b, th)
    } else if (what == "gradient") {
      out[[u$term]][u$cols] <- out[[u$term]][u$cols] +
        penalties7::penalty_gradient(u$penalty, b, th)
    } else {
      out[[u$term]][u$cols, u$cols] <- out[[u$term]][u$cols, u$cols] +
        penalties7::penalty_hessian(u$penalty, b, th)
    }
  }
  out
}


#' The Hyperparameters a Specification Starts From
#'
#' @description
#' Builds the hyperparameter structure a fit begins from: one entry per
#' penalized term, each at the midpoint of its own bounds. That is the probe
#' rule \pkg{modelterms7} already uses when it reads a penalty's kinks, so
#' the two layers start a penalty at the same place.
#'
#' @details
#' A hyperparameter the term holds keeps the held value instead of the probe.
#' Which ones those are is said by the term and by nothing else, through
#' `modelterms7::term_hyper()`.
#'
#' The probe is a placeholder. Every hyperparameter left
#' `NULL` on its term is estimated afterwards, by a marginal criterion or
#' along a path.
#'
#' @param spec A [StatmodSpec()], whose terms are walked equation by
#'   equation.
#'
#' @return A named list with one entry per distribution parameter, each a
#'   named list of numeric vectors keyed by term key, each vector named by
#'   that penalty's own hyperparameters. An equation with no penalized term
#'   holds an empty list.
#'
#' @seealso [penalty_theta_start()] for one penalty's midpoints,
#'   [statmod_hyper_merge()] for a caller's overrides,
#'   [statmod_held()] for which are held.
#'
#' @keywords internal
statmod_hyper_start <- function(spec, design = NULL) {
  params <- spec@distrib@params
  if (is.null(design)) design <- statmod_design(spec)
  out <- stats::setNames(lapply(params, function(p) list()), params)
  for (u in statmod_penalized(spec, design)) {
    th <- penalty_theta_start(u$penalty)
    # A HELD value is not a start, it is the answer: the term said so, and
    # nothing below estimates it away. A hyperparameter the term did not
    # name keeps the probe value, which is a placeholder until a criterion
    # reaches it.
    for (h in names(u$fixed)) th[[h]] <- as.numeric(u$fixed[[h]])
    out[[u$param]][[u$key]] <- th
  }
  out
}


#' Override the Starting Hyperparameters
#'
#' @description
#' Merges a caller's hyperparameters into the ones
#' [statmod_hyper_start()] computed, by parameter and by term.
#'
#' @details
#' Until an outer criterion estimates it, a hyperparameter sits at the probe
#' value, which is a placeholder: a lasso at \eqn{\lambda = 1} against an
#' unaveraged log-likelihood over a few hundred observations selects
#' everything. This is the route by which a caller sets it instead.
#'
#' A vector is matched by name against the penalty's own hyperparameters, so
#' `c(lambda = 5)` sets that one and leaves the rest. An unnamed vector of
#' the full length replaces them all.
#'
#' A term may be named either by the key the specification holds it under,
#' which is its call deparsed, or by its `label`: `"lasso"` in place of
#' `"lasso(~noise1 + noise2)"`. Where two terms share a label the request is
#' ambiguous and the keys are asked for.
#'
#' @param spec A [StatmodSpec()].
#' @param start The hyperparameters as [statmod_hyper_start()] computed them.
#' @param user A named list of named lists, keyed by distribution parameter
#'   and then by term, or `NULL` for no overrides.
#'
#' @return The same structure as `start`, with the caller's values merged in.
#'   `start` unchanged when `user` is `NULL`.
#'
#' @seealso [statmod_hyper_start()] for the structure,
#'   [hyper_key()] for the name matching.
#'
#' @keywords internal
statmod_hyper_merge <- function(spec, start, user) {
  if (is.null(user)) return(start)
  if (!is.list(user) || is.null(names(user))) {
    stop("'hyper' must be a named list, one entry per distribution parameter.",
         call. = FALSE)
  }
  bad <- setdiff(names(user), names(start))
  if (length(bad)) {
    stop(sprintf(paste0("'hyper' names '%s', which is not a parameter.\n",
                        "  They are: %s."),
                 bad[1L], paste(names(start), collapse = ", ")), call. = FALSE)
  }
  for (p in names(user)) {
    u <- user[[p]]
    if (!is.list(u) || is.null(names(u))) {
      stop(sprintf(paste0("'hyper$%s' must be a named list, one entry per\n",
                          "  penalized term. This one has: %s."),
                   p, paste(names(start[[p]]), collapse = ", ")),
           call. = FALSE)
    }
    keys <- vapply(names(u), function(k) hyper_key(spec, start, p, k),
                   character(1))
    for (i in seq_along(u)) {
      nm <- keys[[i]]
      cur <- start[[p]][[nm]]
      v <- u[[i]]
      if (!is.null(names(v))) {
        miss <- setdiff(names(v), names(cur))
        if (length(miss)) {
          stop(sprintf(paste0("'hyper$%s$%s' names '%s'. That penalty's\n",
                              "  hyperparameters are: %s."),
                       p, nm, miss[1L], paste(names(cur), collapse = ", ")),
               call. = FALSE)
        }
        cur[names(v)] <- as.numeric(v)
      } else {
        if (length(v) != length(cur)) {
          stop(sprintf(paste0("'hyper$%s$%s' has length %d but that penalty",
                              " has %d\n  hyperparameters: %s."),
                       p, nm, length(v), length(cur),
                       paste(names(cur), collapse = ", ")), call. = FALSE)
        }
        cur[] <- as.numeric(v)
      }
      start[[p]][[nm]] <- cur
    }
  }
  start
}


#' Resolve a Term's Name Against a Specification
#'
#' @description
#' Turns the name a caller used into the key the specification holds the term
#' under, accepting either that key or the term's `label`.
#'
#' @details
#' A specification keys its terms by the call as written, so a lasso is
#' `"lasso(~noise1 + noise2)"`. The call is what distinguishes two lassos on
#' different covariates, and it is also nothing anybody wants to type. The
#' `label` the term constructor carries is the short form, `"lasso"`, and
#' both are accepted here. Where two terms share a label the request is
#' ambiguous and the keys are asked for, guessing having a fair chance of
#' setting a hyperparameter on the wrong block.
#'
#' @param spec A [StatmodSpec()], read for the terms' labels.
#' @param start The hyperparameter structure, whose names are the keys.
#' @param p The distribution parameter whose equation to look in, a string.
#' @param name What the caller wrote: a key or a label.
#'
#' @return A single string, the key `start[[p]]` holds the term under.
#'   Signals an error when `name` matches nothing, and when it matches the
#'   labels of two terms at once.
#'
#' @seealso [statmod_hyper_merge()], its caller.
#'
#' @keywords internal
hyper_key <- function(spec, start, p, name) {
  keys <- names(start[[p]])
  if (name %in% keys) return(name)
  labels <- vapply(keys, function(k) {
    lb <- tryCatch(spec@terms[[p]][[k]]@label, error = function(e) "")
    if (is.character(lb) && length(lb) == 1L) lb else ""
  }, character(1))
  hit <- keys[labels == name]
  if (length(hit) == 1L) return(hit)
  if (length(hit) > 1L) {
    stop(sprintf(paste0("'hyper$%s$%s' is ambiguous: %d terms of '%s' carry",
                        " that label.\n  Name one of them: %s."),
                 p, name, length(hit), p, paste(hit, collapse = ", ")),
         call. = FALSE)
  }
  stop(sprintf(paste0("'hyper$%s$%s' names no penalized term of '%s'.\n",
                      "  The penalized terms there are: %s."),
               p, name, p,
               if (length(keys)) paste(keys, collapse = ", ") else "none"),
       call. = FALSE)
}


#' A Penalty's Starting Hyperparameters
#'
#' @description
#' Picks a starting value for each of a penalty's hyperparameters from its
#' `params_bounds`, one unit inside whichever ends are infinite:
#'
#' - both ends finite: their midpoint;
#' - bounded below only, the common case: `lower + 1`, so a hyperparameter on
#'   \eqn{[0, \infty)} starts at 1;
#' - bounded above only: `upper - 1`;
#' - unbounded: 1.
#'
#' @details
#' One is the scale a smoothing parameter lives on before anything is known
#' about it. The value matters only as somewhere to begin: it is a probe, and
#' the criterion moves it at the first opportunity.
#'
#' @param pen A \pkg{penalties7} penalty object, read for its
#'   `params_bounds` property alone.
#'
#' @return A named numeric vector, one entry per hyperparameter of `pen`,
#'   named as the penalty names them. `numeric(0)` for a penalty with none,
#'   which a fixed prior is.
#'
#' @seealso [statmod_hyper_start()], its caller.
#'
#' @keywords internal
penalty_theta_start <- function(pen) {
  bounds <- pen@params_bounds
  vapply(bounds, function(b) {
    lo <- b[1L]
    hi <- b[2L]
    if (is.finite(lo) && is.finite(hi)) (lo + hi) / 2
    else if (is.finite(lo)) lo + 1
    else if (is.finite(hi)) hi - 1
    else 1
  }, numeric(1))
}


#' The Objective, Its Gradient and Its Hessian, Stacked
#'
#' @description
#' \eqn{F(\beta) = -\ell(\beta) + \sum_t \rho_t}, unaveraged, over the
#' coefficients of every parameter stacked into one vector.
#'
#' @details
#' # The objective is not divided by the sample size
#'
#' A penalty is a negative log-prior at full size, and a posterior adds a
#' log-likelihood and a log-prior at full size. Averaging the likelihood
#' alone would make a hyperparameter mean something that depends on \eqn{n}.
#'
#' What is scaled instead is the stopping rule, in the one place it is read.
#' See [iwls()] and [iwls_score()].
#'
#' # The five closures share one design
#'
#' All five close over `design`, `hyper` and the two information settings, so
#' a caller who changes a hyperparameter builds a new objective. That is what
#' the outer search does at every point it visits.
#'
#' @param spec A [StatmodSpec()].
#' @param hyper The hyperparameters, per penalized term, held for the life of
#'   the objective.
#' @param design The design. Rebuilt from `spec` when absent.
#' @param expected `TRUE` for the expected information in `he`, `FALSE` for
#'   the observed one.
#' @param approx How the expected information is approximated for a family
#'   with no closed form.
#'
#' @return A list of five functions:
#'   \describe{
#'     \item{`fn(b)`}{the objective at the stacked coefficients `b`, a single
#'       number.}
#'     \item{`gr(b)`}{its gradient, a numeric vector as long as `b`.}
#'     \item{`he(b)`}{its Hessian, a `p x p` matrix.}
#'     \item{`split(b)`}{the stacked vector as a named list, one entry per
#'       distribution parameter.}
#'     \item{`stack(l)`}{the inverse of `split`.}
#'   }
#'
#' @seealso [iwls()] and [fit_smooth()], which consume this,
#'   [statmod_loglik_at()] and [statmod_penalty_at()] for its two halves.
#'
#' @keywords internal
statmod_objective <- function(spec, hyper, design = statmod_design(spec),
                              expected = TRUE, approx = "bartlett") {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar

  split <- function(v) {
    stats::setNames(lapply(seq_along(params), function(a)
      v[offs[a] + seq_len(npar[a])]), params)
  }
  stack <- function(l) unlist(l[params], use.names = FALSE)

  list(
    npar = npar,
    split = split,
    stack = stack,
    fn = function(v) {
      cf <- split(v)
      -statmod_loglik_at(spec, cf, design) +
        statmod_penalty_at(spec, cf, hyper, design, "value")
    },
    gr = function(v) {
      cf <- split(v)
      s <- statmod_score_at(spec, cf, design)
      pg <- statmod_penalty_at(spec, cf, hyper, design, "gradient")
      stack(Map(function(a, b) -a + b, s, pg))
    },
    he = function(v) {
      cf <- split(v)
      statmod_information_at(spec, cf, design, expected, approx) +
        statmod_penalty_at(spec, cf, hyper, design, "hessian")
    }
  )
}
