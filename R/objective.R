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
#' Turns a coefficient structure into the per-observation parameters the
#' distribution's generics take.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design, as \code{\link{statmod_design}} returns it.
#' @param coef A named list of coefficient vectors, one per parameter.
#'
#' @return A list with \code{eta} and \code{theta}, both named lists.
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
#' \eqn{\sum_i w_i \log f(y_i; \theta_i)}, the weights entering as given.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design; recomputed when absent.
#'
#' @return A single number.
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
                                    log = TRUE)
  sum(spec@weights * ll)
}


#' The Score of the Weighted Log-Likelihood
#'
#' @description
#' One block per distribution parameter, \eqn{X_k'(w\,g_k)} with \eqn{g_k} the
#' per-observation derivative of the log-density in the link-scale predictor.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return A named list of gradient vectors, one per parameter.
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
    for (k in seq_along(r$mu)) {
      thk <- statmod_theta_shifted(spec, ev$eta_static, r$param, r$mu[[k]])
      gk <- distributions7::distrib_gradient(spec@distrib, spec@response, thk,
                                             scale = "link")
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
                                        scale = "link")
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
                                    g = gv[[f$param]])
    H <- distributions7::distrib_hessian(spec@distrib, spec@response,
                                         ev$theta, scale = "link")
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
#' The filter's jacobian is the total derivative of the predictor in the
#' term's parameters, the recursion included, so the chain rule over the
#' observations is all that is needed here and no reverse pass is: what the
#' reverse pass answers is the other question, the derivative in the
#' coefficients of the equations.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return A named list of numeric vectors, one per structural term.
#'
#' @seealso \code{\link{statmod_filter_at}}
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
                                        scale = "link")
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
#' `expected = TRUE` gives the expected information, which is what Fisher
#' scoring inverts; `approx` is passed to \pkg{distributions7} and read only
#' where the family has no closed expected information.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#' @param expected Whether to use the expected information.
#' @param approx The approximation, when the expected one is not closed.
#'
#' @return A square matrix over the stacked coefficients.
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
    out <- matrix(0, total, total)
    for (k in seq_along(r$mu)) {
      thk <- statmod_theta_shifted(spec, ev$eta_static, r$param, r$mu[[k]])
      Hk <- if (expected) {
        distributions7::distrib_expected_hessian(spec@distrib, spec@response,
                                                 thk, scale = "link",
                                                 approx = approx)
      } else {
        distributions7::distrib_hessian(spec@distrib, spec@response, thk,
                                        scale = "link")
      }
      for (a in seq_along(params)) {
        if (npar[a] == 0L) next
        for (b in seq_along(params)) {
          if (npar[b] == 0L || b < a) next
          wv <- spec@weights * r$gamma[, k] *
            rep_len(Hk[[hess_key(params, a, b)]], spec@n_obs)
          blk <- crossprod(design[[params[a]]]$X * wv, design[[params[b]]]$X)
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
                                             scale = "link", approx = approx)
  } else {
    distributions7::distrib_hessian(spec@distrib, spec@response, th,
                                    scale = "link")
  }
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  out <- matrix(0, total, total)
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in seq_along(params)) {
      if (npar[b] == 0L || b < a) next
      key <- hess_key(params, a, b)
      wv <- spec@weights * rep_len(H[[key]], spec@n_obs)
      blk <- crossprod(design[[params[a]]]$X * wv, design[[params[b]]]$X)
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
#' keyed by \code{\link[distributions7]{hess_names}} and not by position.
#'
#' @details
#' The name is BUILT from the parameter names rather than parsed out of one,
#' the discipline \pkg{distributions7} records for a parameter whose own name
#' contains the separator.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b Indices into \code{params}.
#'
#' @return A single string.
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
#' The sum of every penalized term's \code{penalty_value()}, with its gradient
#' and Hessian placed in the columns that term owns.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param hyper A named list, one entry per parameter, each a named list of
#'   hyperparameter vectors per penalized term.
#' @param design The design.
#' @param what One of \code{"value"}, \code{"gradient"}, \code{"hessian"}.
#'
#' @return A number, a named list of vectors, or a square matrix.
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
  hess <- matrix(0, total, total)

  for (u in statmod_penalized(spec, design)) {
    p <- u$param
    b <- coef[[p]][u$cols]
    th <- as.list(hyper[[p]][[u$key]])
    if (what == "value") {
      value <- value + penalties7::penalty_value(u$penalty, b, th)
    } else if (what == "gradient") {
      grad[[p]][u$cols] <- grad[[p]][u$cols] +
        penalties7::penalty_gradient(u$penalty, b, th)
    } else {
      hess[u$index, u$index] <- hess[u$index, u$index] +
        penalties7::penalty_hessian(u$penalty, b, th)
    }
  }
  switch(what, value = value, gradient = grad, hessian = hess)
}


#' The Hyperparameters a Specification Starts From
#'
#' @description
#' One entry per penalized term, at the midpoint of its bounds -- the probe
#' rule \pkg{modelterms7} already uses when it reads a penalty's kinks.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#'
#' @return A named list, one entry per parameter.
#'
#' @keywords internal
statmod_hyper_start <- function(spec, design = NULL) {
  params <- spec@distrib@params
  if (is.null(design)) design <- statmod_design(spec)
  out <- stats::setNames(lapply(params, function(p) list()), params)
  for (u in statmod_penalized(spec, design)) {
    out[[u$param]][[u$key]] <- penalty_theta_start(u$penalty)
  }
  out
}


#' Override the Starting Hyperparameters
#'
#' @description
#' Merges a caller's hyperparameters into the ones
#' \code{\link{statmod_hyper_start}} computed, by parameter and by term.
#'
#' @details
#' Until a hyperparameter is estimated by an outer criterion it is held at the
#' probe value, which is a placeholder and not a choice: a lasso at
#' \eqn{\lambda = 1} against an unaveraged log-likelihood of a few hundred
#' observations selects nothing. This is what lets a caller set it.
#'
#' A vector is matched by name against the penalty's own hyperparameters, so
#' \code{c(lambda = 5)} sets that one and leaves the rest where they were; an
#' unnamed vector of the full length replaces them all.
#'
#' A term is named either by the key the specification holds it under, which is
#' its call deparsed, or by its \code{label} -- \code{"lasso"} rather than
#' \code{"lasso(~noise1 + noise2)"}. Two terms sharing a label are ambiguous
#' and the keys are asked for instead.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param start The hyperparameters as computed.
#' @param user A named list of named lists, or \code{NULL}.
#'
#' @return The merged structure.
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
#' under, accepting either that key or the term's \code{label}.
#'
#' @details
#' A specification keys its terms by the call as written, so a lasso is
#' \code{"lasso(~noise1 + noise2)"}. That is what makes two lassos on different
#' covariates distinct, and it is not what anybody wants to type; the label the
#' term constructor carries is. Where two terms share a label the request is
#' ambiguous and the keys are asked for, since guessing would set a
#' hyperparameter on the wrong block.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param start The hyperparameter structure, whose names are the keys.
#' @param p The distribution parameter.
#' @param name What the caller wrote.
#'
#' @return A single key.
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
#' The midpoint of each hyperparameter's bounds, finite ends taken as they
#' are and an infinite one replaced by one, which is the scale a smoothing
#' parameter lives on before anything is known about it.
#'
#' @param pen A \pkg{penalties7} penalty.
#'
#' @return A named numeric vector.
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
#' The objective is NOT divided by the sample size. A penalty is a negative
#' log-prior at full size, and a posterior is a log-likelihood plus a log-prior
#' at full size; averaging the first would make a hyperparameter mean something
#' that depends on \eqn{n}. What is scaled instead is the stopping rule, in the
#' one place it is read.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param hyper The hyperparameters.
#' @param design The design.
#' @param expected Whether the information is the expected one.
#' @param approx The approximation for the expected information.
#'
#' @return A list of functions \code{fn}, \code{gr} and \code{he} of the
#'   stacked coefficient vector, plus \code{split} and \code{stack} to move
#'   between that vector and the per-parameter list.
#'
#' @seealso \code{\link{iwls}}
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
