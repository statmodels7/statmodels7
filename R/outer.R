#' @include inference.R
NULL

# The hyperparameters, estimated by a marginal criterion.
#
# The penalized objective is F(b; t) = -l(b) + sum_j rho_j(b; t), and every
# penalties7 penalty keeps its normalizing constant, so rho IS minus a log
# prior density and nothing has to be added by hand. The marginal likelihood
# is then
#
#   L(t) = integral exp(-F(b; t)) db
#
# and its Laplace approximation at the penalized mode b(t) is
#
#   log L(t) = l(b) - rho(b; t) + (q/2) log(2 pi) - (1/2) log|A'(H + S)A|
#
# with H the information of the likelihood, S the penalty's Hessian, and A an
# orthonormal basis of whichever subspace is integrated over. Writing it out
# with the constants in place reproduces Wood's (2011) REML criterion term for
# term, the -(M_p/2) log(2 pi) of that paper appearing here as the difference
# between the rank of the penalty and the dimension of the coefficient vector.

#' How the Hyperparameters Are Estimated
#'
#' @description
#' The criterion \code{\link{reml}()} and \code{\link{ml}()} build: which
#' subspace of the coefficients is integrated over, and which information
#' matrix enters the determinant.
#'
#' @param kind \code{"reml"} or \code{"ml"}.
#' @param hessian \code{"expected"} or \code{"observed"}.
#'
#' @return An object of class \code{OuterMethod}.
#'
#' @seealso \code{\link{reml}}, \code{\link{ml}}, \code{\link{statmod}}
#'
#' @examples
#' reml()
#' ml(hessian = "observed")
#'
#' @name OuterMethod-class
#' @aliases OuterMethod
#' @keywords internal
#' @export
OuterMethod <- S7::new_class("OuterMethod",
  properties = list(
    kind = S7::class_character,
    hessian = S7::class_character
  ),
  validator = function(self) {
    if (!identical(length(self@kind), 1L) ||
        !self@kind %in% c("reml", "ml")) {
      return("Property 'kind' must be \"reml\" or \"ml\".")
    }
    if (!identical(length(self@hessian), 1L) ||
        !self@hessian %in% c("expected", "observed")) {
      return("Property 'hessian' must be \"expected\" or \"observed\".")
    }
    NULL
  }
)


#' Estimate the Hyperparameters by a Marginal Likelihood
#'
#' @description
#' \code{reml()} integrates every coefficient out of the likelihood before
#' maximizing in the hyperparameters; \code{ml()} integrates only the
#' penalized directions and profiles the rest.
#'
#' @details
#' \strong{The criterion.} At the penalized mode \eqn{\hat\beta(\theta)},
#' \deqn{\log L(\theta) = \ell(\hat\beta) - \rho(\hat\beta;\theta)
#'   + \frac{q}{2}\log 2\pi - \frac12\log|A'(H+S)A|,}
#' with \eqn{H} the information of the log-likelihood, \eqn{S} the penalty's
#' second derivative in the coefficients, and \eqn{A} an orthonormal basis of
#' the subspace integrated over. Nothing is added to \eqn{\rho} to make this
#' work: a \pkg{penalties7} penalty keeps its normalizing constant, so it is
#' exactly minus a log prior density, and for a quadratic penalty that constant
#' carries the \eqn{-\frac{r}{2}\log\lambda} and the log pseudo-determinant
#' that a marginal criterion needs. Written out, the expression reproduces
#' Wood's (2011) REML criterion term for term.
#'
#' \strong{What each one integrates.} \code{reml()} takes \eqn{A = I}: every
#' coefficient is integrated, the unpenalized ones under the flat prior their
#' absence of a penalty amounts to. \code{ml()} takes \eqn{A} spanning the
#' range space of the penalty, so a coefficient that is unpenalized -- an
#' ordinary covariate, or the linear component of a Demmler-Reinsch smooth,
#' which its penalty leaves alone -- is profiled rather than integrated. This
#' is the same distinction as between REML and ML for a variance component in a
#' mixed model, and \code{reml()} is the default for the same reason: profiling
#' a fixed effect leaves the estimate of the variance biased downwards.
#'
#' \strong{Which hyperparameters.} Those of the terms fitted in one system,
#' which is to say those whose penalty is twice differentiable. A lasso, a SCAD
#' or an MCP has a kink, its coefficients are estimated by a method of their
#' own, and a Laplace approximation at a point where the second derivative does
#' not exist would be arithmetic without a meaning; those hyperparameters stay
#' where \code{hyper} put them.
#'
#' \strong{The criterion has an exact gradient} where the information is the
#' observed one and every penalty under estimation has a Hessian linear in its
#' hyperparameters, which covers \code{s()}, \code{te()} and any
#' \code{\link[penalties7]{quadratic_penalty}}. It is then supplied to the
#' search and \code{\link[optimizers7]{lbfgs}} becomes the default optimizer;
#' otherwise the search compares values. Measured, in evaluations of the
#' criterion (each a whole inner fit) against
#' \code{\link[optimizers7]{nelder_mead}}: 40 against 32 with one smoothing
#' parameter, 40 against 135 with two, 41 against 269 with three, and 12
#' against 283 with three and a modelled scale. It does not pay in one
#' dimension and pays from two on, a simplex needing a vertex per dimension
#' and a quasi-Newton method not. See \code{\link{statmod_marginal_grad}}.
#'
#' \strong{ML needs a null basis} for every penalty that has one, since that is
#' what says which directions are profiled. \code{\link[penalties7]{is_proper}}
#' answers for a penalty with no null space at all, and
#' \code{\link[penalties7]{penalty_null_basis}} for the quadratic and
#' structured branches. A penalty that has neither is rejected by name rather
#' than integrated over a subspace guessed at.
#'
#' @param hessian Which information enters the determinant: \code{"expected"}
#'   or \code{"observed"}.
#'
#' @return An \code{\link{OuterMethod}}.
#'
#' @references
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal
#' likelihood estimation of semiparametric generalized linear models.
#' \emph{Journal of the Royal Statistical Society, Series B}, 73(1), 3--36.
#'
#' @seealso \code{\link{statmod}}, \code{\link{iwls}}
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
#' statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
#'         outer_method = reml())
#'
#' @export
reml <- function(hessian = c("expected", "observed")) {
  OuterMethod(kind = "reml", hessian = match.arg(hessian))
}

#' @rdname reml
#' @export
ml <- function(hessian = c("expected", "observed")) {
  OuterMethod(kind = "ml", hessian = match.arg(hessian))
}


#' @title Print an Outer Method
#' @name print.OuterMethod
#' @description Which criterion it is and which information it uses.
#' @param x An \code{\link{OuterMethod}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{reml}}
#' @keywords internal
print.OuterMethod <- function(x, ...) {
  cat(sprintf("<%s>  %s information\n", toupper(x@kind), x@hessian))
  invisible(x)
}
S7::method(print, OuterMethod) <- print.OuterMethod


#' The Hyperparameters an Outer Method Estimates
#'
#' @description
#' One row per hyperparameter of a penalized term that is fitted in the joint
#' system, with the link that carries it onto the whole line.
#'
#' @details
#' A term whose penalty has a kink is left out. Its coefficients are estimated
#' by a proximal method with everything else held fixed, and the criterion is a
#' Laplace approximation, which asks for a second derivative that does not
#' exist at a kink.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param blocks The block split, as \code{\link{statmod_blocks}} returns it.
#'
#' @return A data frame with \code{parameter}, \code{term} and \code{name}, and
#'   a list column-free \code{link} carried as an attribute list.
#'
#' @keywords internal
outer_hyper_index <- function(spec, blocks) {
  kinked <- vapply(blocks$sparse, function(b)
    paste(b$param, b$term, sep = "\r"), character(1))
  rows <- list()
  links <- list()
  for (p in spec@distrib@params) {
    for (nm in names(spec@terms[[p]])) {
      pen <- modelterms7::term_penalty(spec@terms[[p]][[nm]])
      if (is.null(pen)) next
      if (paste(p, nm, sep = "\r") %in% kinked) next
      for (h in pen@params) {
        rows[[length(rows) + 1L]] <- data.frame(
          parameter = p, term = nm, name = h, stringsAsFactors = FALSE)
        links[[length(links) + 1L]] <- pen@link_params[[h]]
      }
    }
  }
  if (!length(rows)) {
    return(structure(data.frame(parameter = character(0), term = character(0),
                                name = character(0),
                                stringsAsFactors = FALSE),
                     links = list()))
  }
  structure(do.call(rbind, rows), links = links)
}


#' Move Between the Hyperparameters and the Free Vector
#'
#' @description
#' \code{hyper_to_eta()} carries the estimated hyperparameters onto the whole
#' line through their links; \code{eta_to_hyper()} puts a free vector back.
#'
#' @details
#' The outer search runs on the free scale for the reason every other search in
#' the toolkit does: a smoothing parameter is positive, and an optimizer that
#' does not know that will step outside its domain.
#'
#' @param hyper The hyperparameter structure.
#' @param idx The index, as \code{\link{outer_hyper_index}} returns it.
#' @param eta A free vector.
#'
#' @return A numeric vector, or the hyperparameter structure.
#'
#' @keywords internal
hyper_to_eta <- function(hyper, idx) {
  links <- attr(idx, "links")
  vapply(seq_len(nrow(idx)), function(i) {
    v <- hyper[[idx$parameter[i]]][[idx$term[i]]][[idx$name[i]]]
    linkfunctions7::linkfun(links[[i]], v)
  }, numeric(1))
}

#' @rdname hyper_to_eta
#' @keywords internal
eta_to_hyper <- function(eta, idx, hyper) {
  links <- attr(idx, "links")
  for (i in seq_len(nrow(idx))) {
    hyper[[idx$parameter[i]]][[idx$term[i]]][[idx$name[i]]] <-
      linkfunctions7::linkinv(links[[i]], eta[[i]])
  }
  hyper
}


#' The Subspace a Marginal Criterion Integrates Over
#'
#' @description
#' \code{NULL} for REML, which integrates everything, and an orthonormal basis
#' of the penalty's range space for ML.
#'
#' @details
#' The total penalty is block diagonal over the terms, each owning its own
#' columns, so its null space is the direct sum of the unpenalized columns and
#' each penalty's own null space. It is assembled that way and never read off
#' the eigenvalues of the assembled matrix, which measure the arithmetic rather
#' than the family: with two smoothing parameters ten orders of magnitude
#' apart, an eigenvalue count of a sum falls while the null space does not
#' move.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param kind \code{"reml"} or \code{"ml"}.
#'
#' @return \code{NULL}, or a matrix with one column per integrated direction.
#'
#' @keywords internal
integrated_basis <- function(spec, design, kind) {
  if (identical(kind, "reml")) return(NULL)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  cols <- list()
  for (a in seq_along(params)) {
    p <- params[a]
    for (nm in names(spec@terms[[p]])) {
      pen <- modelterms7::term_penalty(spec@terms[[p]][[nm]])
      if (is.null(pen)) next
      blk <- design[[p]]$blocks[[nm]]
      k <- length(blk)
      R <- penalty_range_basis(pen, k, p, nm)
      if (!ncol(R)) next
      M <- matrix(0, total, ncol(R))
      M[offs[a] + blk, ] <- R
      cols[[length(cols) + 1L]] <- M
    }
  }
  if (!length(cols)) {
    stop(paste0("ml() found no penalized direction to integrate over.\n",
                "  Every penalty here is entirely null, which is not a model",
                " ml() can\n  read; reml() integrates the whole coefficient",
                " vector instead."), call. = FALSE)
  }
  do.call(cbind, cols)
}


#' An Orthonormal Basis of a Penalty's Range Space
#'
#' @description
#' The directions a penalty constrains, which are the ones \code{\link{ml}()}
#' integrates over.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param k The number of coefficients in the term's block.
#' @param p The distribution parameter, for the message.
#' @param nm The term's name, for the message.
#'
#' @return A \code{k} by \code{r} matrix.
#'
#' @keywords internal
penalty_range_basis <- function(pen, k, p, nm) {
  if (isTRUE(penalties7::is_proper(pen))) return(diag(1, k, k))
  N <- tryCatch(penalties7::penalty_null_basis(pen), error = function(e) NULL)
  if (is.null(N)) {
    stop(sprintf(paste0("ml() cannot read the null space of the penalty on",
                        " '%s' in '%s'.\n  Its class does not expose one, so",
                        " which directions are profiled is\n  not something",
                        " to guess at. Use reml(), which integrates every",
                        "\n  coefficient and needs no such basis."), nm, p),
         call. = FALSE)
  }
  N <- as.matrix(N)
  if (!ncol(N)) return(diag(1, k, k))
  # the range is the orthogonal complement of the null space inside the block
  Q <- qr.Q(qr(N), complete = TRUE)
  Q[, seq.int(ncol(N) + 1L, k), drop = FALSE]
}


#' The Marginal Criterion at Given Coefficients and Hyperparameters
#'
#' @description
#' The Laplace approximation to the log marginal likelihood, evaluated at the
#' penalized mode.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients, at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An \code{\link{OuterMethod}}.
#' @param approx The approximation for the expected information.
#' @param basis The integrated subspace, from \code{\link{integrated_basis}}.
#'
#' @return A list with \code{value}, \code{loglik}, \code{penalty},
#'   \code{logdet} and \code{q}, or \code{NULL} where the determinant does not
#'   exist.
#'
#' @keywords internal
statmod_marginal <- function(spec, design, coef, hyper, method,
                             approx = "bartlett", basis = NULL) {
  expected <- identical(method@hessian, "expected")
  ll <- statmod_loglik_at(spec, coef, design)
  rho <- statmod_penalty_at(spec, coef, hyper, design, "value")
  H <- statmod_information_at(spec, coef, design, expected, approx)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S[!is.finite(S)] <- 0
  M <- H + S
  if (!is.null(basis)) M <- crossprod(basis, M %*% basis)
  R <- tryCatch(chol(M), error = function(e) NULL)
  # a determinant that does not exist is a hyperparameter the search should
  # step away from, not an error: at a far-out value the penalized information
  # can lose definiteness while the fit itself is sound
  if (is.null(R)) return(NULL)
  q <- nrow(M)
  logdet <- 2 * sum(log(diag(R)))
  list(value = ll - rho + q / 2 * log(2 * pi) - logdet / 2,
       loglik = ll, penalty = rho, logdet = logdet, q = q)
}


#' Estimate the Hyperparameters
#'
#' @description
#' Runs \code{outer_optimizer} on the marginal criterion, refitting the
#' coefficients at every hyperparameter it tries.
#'
#' @details
#' Each evaluation is a whole inner fit, warm-started from the previous one,
#' which is what makes the search affordable: after the first few
#' hyperparameters the coefficients move very little and the inner loop
#' converges in two or three iterations.
#'
#' \strong{The optimizer is chosen by whether the gradient exists.} Where
#' \code{\link{statmod_marginal_grad}} applies -- the observed information, and
#' penalties whose Hessian is linear in their hyperparameters -- the criterion
#' is handed its exact derivative and \code{\link[optimizers7]{lbfgs}} is the
#' default; otherwise the search compares values and
#' \code{\link[optimizers7]{nelder_mead}} is. An optimizer given explicitly is
#' used as given, and one that needs a gradient it cannot be given will say so
#' itself.
#'
#' @param spec The specification.
#' @param design The design.
#' @param blocks The block split.
#' @param hyper The starting hyperparameters.
#' @param inner_method How the smooth block is fitted.
#' @param method An \code{\link{OuterMethod}}.
#' @param optimizer An \pkg{optimizers7} optimizer, or \code{NULL} to let the
#'   availability of the exact gradient decide.
#' @param beta The starting coefficients, stacked.
#' @param approx The approximation for the expected information.
#' @param maxit,tol The alternation's budget and tolerance.
#' @param vb The resolved verbosity.
#'
#' @return A list with \code{par}, \code{hyper}, \code{value}, \code{criterion},
#'   \code{converged}, \code{history} and the inner results.
#'
#' @seealso \code{\link{reml}}, \code{\link{statmod}}
#'
#' @keywords internal
outer_fit <- function(spec, design, blocks, hyper, inner_method, method,
                      optimizer, beta, approx, maxit, tol, vb) {
  idx <- outer_hyper_index(spec, blocks)
  if (!nrow(idx)) {
    stop(paste0("outer_method was given but there is no hyperparameter to\n",
                "  estimate: no term here carries a penalty that is twice\n",
                "  differentiable. A lasso, a scad or an mcp keeps the value",
                "\n  'hyper' gave it, its penalty having a kink."),
         call. = FALSE)
  }
  expected <- identical(method@hessian, "expected")
  basis <- integrated_basis(spec, design, method@kind)
  labels <- paste(idx$parameter, idx$term, idx$name, sep = "/")
  exact <- outer_gradient_ok(spec, design, idx, method)
  if (is.null(optimizer)) {
    optimizer <- if (exact) optimizers7::lbfgs() else
      optimizers7::nelder_mead()
  }

  state <- new.env(parent = emptyenv())
  state$beta <- beta
  state$inner <- NULL
  state$rows <- list()
  state$evals <- 0L
  state$key <- NULL

  # one inner fit serves the value and the gradient: an optimizer asks for
  # both at the same point, and refitting for the second would double the cost
  # of every step
  evaluate <- function(eta) {
    key <- paste(format(eta, digits = 17), collapse = ",")
    if (identical(state$key, key)) return(state$last)
    hy <- eta_to_hyper(eta, idx, hyper)
    res <- statmod_alternate(spec, design, blocks, hy, inner_method,
                             state$beta, expected, approx, maxit, tol,
                             vb_inner(vb))
    cf <- res$obj$split(res$par)
    m <- statmod_marginal(spec, design, cf, hy, method, approx, basis)
    state$evals <- state$evals + 1L
    if (is.null(m)) {
      if (vb$outer) {
        cat(sprintf("[outer %3d] criterion unavailable at %s\n", state$evals,
                    paste(signif(eta, 4), collapse = ", ")))
      }
      out <- list(value = -Inf, grad = rep(0, nrow(idx)))
      state$key <- key
      state$last <- out
      return(out)
    }
    state$beta <- res$par
    state$inner <- res
    state$hyper <- hy
    g <- if (exact) statmod_marginal_grad(spec, design, cf, hy, method, idx,
                                          basis) else NULL
    row <- data.frame(evaluation = state$evals, criterion = m$value,
                      loglik = m$loglik, penalty = m$penalty)
    for (j in seq_along(labels)) {
      row[[labels[j]]] <- hyper_value(hy, idx, j)
    }
    state$rows[[length(state$rows) + 1L]] <- row
    if (vb$outer) {
      cat(sprintf("[outer %3d] %s = %s    %s %.6f\n", state$evals,
                  paste(labels, collapse = ", "),
                  paste(signif(vapply(seq_along(labels), function(j)
                    hyper_value(hy, idx, j), numeric(1)), 5), collapse = ", "),
                  toupper(method@kind), m$value))
    }
    out <- list(value = m$value,
                grad = if (is.null(g)) rep(0, nrow(idx)) else g)
    state$key <- key
    state$last <- out
    out
  }

  fn <- function(eta) -evaluate(eta)$value
  gr <- function(eta) -evaluate(eta)$grad

  eta0 <- hyper_to_eta(hyper, idx)
  res <- if (exact) optimizers7::minimize(optimizer, fn, eta0, gr = gr) else
    optimizers7::minimize(optimizer, fn, eta0)

  # the last evaluation is not necessarily the optimum, so the fit is taken at
  # the reported point rather than at whatever was tried last
  hy <- eta_to_hyper(res@par, idx, hyper)
  inner <- statmod_alternate(spec, design, blocks, hy, inner_method,
                             state$beta, expected, approx, maxit, tol, vb)
  m <- statmod_marginal(spec, design, inner$obj$split(inner$par), hy, method,
                        approx, basis)
  list(par = inner$par, hyper = hy, value = inner$value,
       criterion = if (is.null(m)) NA_real_ else m$value,
       converged = res@converged && inner$converged,
       obj = inner$obj, hist_blocks = inner$hist_blocks,
       hist_inner = inner$hist_inner,
       hist_outer = if (length(state$rows)) do.call(rbind, state$rows) else
         NULL,
       iterations = res@iterations, evaluations = state$evals,
       exact_gradient = exact)
}


#' One Hyperparameter of an Index
#'
#' @description
#' The value the \code{j}-th row of an outer index points at.
#'
#' @param hyper The hyperparameter structure.
#' @param idx The index.
#' @param j A row.
#'
#' @return A single number.
#'
#' @keywords internal
hyper_value <- function(hyper, idx, j) {
  as.numeric(hyper[[idx$parameter[j]]][[idx$term[j]]][[idx$name[j]]])[1L]
}


#' The Verbosity of an Inner Fit Inside the Outer Search
#'
#' @description
#' The same switches with the block trace off, since one line per sweep per
#' hyperparameter tried is not a trace anybody reads.
#'
#' @param vb The resolved verbosity.
#'
#' @return A list of switches.
#'
#' @keywords internal
vb_inner <- function(vb) {
  vb$blocks <- vb$blocks && vb$inner
  vb
}
