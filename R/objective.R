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
  eta <- stats::setNames(vector("list", length(params)), params)
  theta <- eta
  for (p in params) {
    d <- design[[p]]
    e <- if (d$npar == 0L) rep(0, spec@n_obs) else
      as.numeric(d$X %*% coef[[p]])
    off <- spec@offsets[[p]]
    if (!is.null(off)) e <- e + off
    eta[[p]] <- e
    theta[[p]] <- linkfunctions7::linkinv(links[[p]], e)
  }
  list(eta = eta, theta = theta)
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
  th <- statmod_eta(spec, design, coef)$theta
  ll <- distributions7::distrib_pdf(spec@distrib, spec@response, th,
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
  th <- statmod_eta(spec, design, coef)$theta
  g <- distributions7::distrib_gradient(spec@distrib, spec@response, th,
                                        scale = "link")
  stats::setNames(lapply(params, function(p) {
    d <- design[[p]]
    if (d$npar == 0L) return(numeric(0))
    as.numeric(crossprod(d$X, spec@weights * rep_len(g[[p]], spec@n_obs)))
  }), params)
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
  th <- statmod_eta(spec, design, coef)$theta
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

  for (a in seq_along(params)) {
    p <- params[a]
    for (nm in names(spec@terms[[p]])) {
      pen <- modelterms7::term_penalty(spec@terms[[p]][[nm]])
      if (is.null(pen)) next
      cols <- design[[p]]$blocks[[nm]]
      b <- coef[[p]][cols]
      th <- hyper[[p]][[nm]]
      if (what == "value") {
        value <- value + penalties7::penalty_value(pen, b, th)
      } else if (what == "gradient") {
        grad[[p]][cols] <- grad[[p]][cols] +
          penalties7::penalty_gradient(pen, b, th)
      } else {
        idx <- offs[a] + cols
        hess[idx, idx] <- hess[idx, idx] +
          penalties7::penalty_hessian(pen, b, th)
      }
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
statmod_hyper_start <- function(spec) {
  params <- spec@distrib@params
  stats::setNames(lapply(params, function(p) {
    tms <- spec@terms[[p]]
    out <- list()
    for (nm in names(tms)) {
      pen <- modelterms7::term_penalty(tms[[nm]])
      if (is.null(pen)) next
      out[[nm]] <- penalty_theta_start(pen)
    }
    out
  }), params)
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
