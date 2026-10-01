#' @include statmod.R
NULL

#' The Random-Effect Terms a Prediction Can Set Aside
#'
#' @description
#' Lists the [modelterms7::random()] terms written in the model's equations,
#' with the mode `predict()` reads each at.
#'
#' @details
#' `random` is a single string, applied to every such term, or a character
#' vector named by the terms' keys, one mode per term named, the rest staying
#' conditional. A key is matched with its white space removed, so
#' `"random(~1|g)"` finds the term written `random(~ 1 | g)`. A name that
#' matches no term signals an error listing those there are.
#'
#' Only a term written in an equation is listed. One written inside a
#' subformula develops another term's parameter and is predicted as that
#' term predicts it.
#'
#' @param spec The fit's specification.
#' @param random The mode or modes, as [predict.StatmodFit()] takes them.
#'
#' @return A data frame with columns `param`, `key` and `mode`, one row per
#'   random-effect term, or zero rows where the model has none.
#'
#' @keywords internal
random_modes <- function(spec, random) {
  modes <- c("conditional", "zero", "marginal")
  rows <- list()
  for (p in names(spec@terms)) {
    for (k in names(spec@terms[[p]])) {
      if (S7::S7_inherits(spec@terms[[p]][[k]], modelterms7::RandomTerm)) {
        rows[[length(rows) + 1L]] <- data.frame(param = p, key = k,
                                                stringsAsFactors = FALSE)
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(param = character(0), key = character(0),
               stringsAsFactors = FALSE)
  if (!is.character(random) || !length(random) || anyNA(random) ||
      !all(random %in% modes)) {
    stop("'random' must be \"conditional\", \"zero\" or \"marginal\", or a ",
         "character vector of\n  those named by the random-effect terms.",
         call. = FALSE)
  }
  if (is.null(names(random))) {
    if (length(random) != 1L) {
      stop("An unnamed 'random' must be a single mode.", call. = FALSE)
    }
    out$mode <- rep(random, nrow(out))
    return(out)
  }
  squash <- function(x) gsub("[[:space:]]", "", x)
  out$mode <- rep("conditional", nrow(out))
  for (nm in names(random)) {
    hit <- squash(out$key) == squash(nm)
    if (!any(hit)) {
      stop(sprintf(paste0("'random' names '%s', which is not a random-effect ",
                          "term of the model.\n  They are: %s."), nm,
                   if (nrow(out)) paste(unique(out$key), collapse = ", ")
                   else "none"), call. = FALSE)
    }
    out$mode[hit] <- random[[nm]]
  }
  out
}


#' The Prior of One Random-Effect Term, for Integrating a New Group
#'
#' @description
#' Nodes and weights, or draws, for one group's effect under the prior the
#' fit estimated for it.
#'
#' @details
#' A term whose penalty is quadratic in the coefficients carries a Gaussian
#' prior, read off the penalty's Hessian: its first \eqn{d\times d} block is
#' one group's precision \eqn{\Omega}, the block being
#' \eqn{I_m \otimes \Omega}. Any other prior, a Student t for one, is
#' sampled with [penalties7::penalty_draw()], one draw of the block giving
#' \eqn{m} independent groups.
#'
#' @param spec,design The specification and its design.
#' @param fit The fitted model, for the hyperparameters.
#' @param param,key The term.
#'
#' @return A list with `gaussian` (logical), `dim` and, for a Gaussian
#'   prior, `chol` (the lower Cholesky factor of the covariance); for any
#'   other, `pen`, `theta` and the group count `m`.
#'
#' @keywords internal
random_prior <- function(spec, design, fit, param, key) {
  tm <- spec@terms[[param]][[key]]
  if (!is.na(modelterms7::term_tag(tm))) {
    stop(sprintf(paste0("'%s' shares a covariance with other terms, so its ",
                        "prior is not its own\n  and random = \"marginal\" ",
                        "cannot integrate it alone. random = \"zero\" can."),
                 key), call. = FALSE)
  }
  units <- Filter(function(u) identical(u$param, param) &&
                    identical(u$term, key), statmod_penalized(spec, design))
  if (length(units) != 1L) {
    stop(sprintf("'%s' has no single prior to integrate over.", key),
         call. = FALSE)
  }
  u <- units[[1L]]
  pen <- u$penalty
  th <- as.list(fit@hyper[[param]][[u$key]])
  gr <- modelterms7::term_group(tm)
  d <- gr$dim
  m <- length(gr$levels)
  if (isTRUE(penalties7::beta_quadratic(pen, th))) {
    H <- as.matrix(penalties7::penalty_hessian(pen, rep(0, pen@n_coef), th))
    Om <- H[seq_len(d), seq_len(d), drop = FALSE]
    L <- t(chol(solve(Om)))
    return(list(gaussian = TRUE, dim = d, chol = L))
  }
  # A HEAVY-TAILED PRIOR MAY HAVE NO AVERAGE TO TAKE. A Student t has no
  # moment generating function, so under a log link E[exp(b)] is infinite,
  # and a Monte Carlo mean returns whatever its largest draw makes it:
  # measured on a Poisson random intercept at nu = 3.04, 10000 draws gave
  # 31.5, 2.9, 28.7, 57.9 and 331 on five seeds, against 1.81 at a zero
  # effect. Where the inverse link is bounded the average always exists, and
  # the same draws gave 0.628 to 0.631 under a logit.
  lb <- spec@distrib@link_params[[param]]@link_bounds
  if (!all(is.finite(lb))) {
    stop(sprintf(paste0("The prior of '%s' is not Gaussian and the link of ",
                        "'%s' is not bounded,\n  so the average over new ",
                        "groups may not exist (under a log link it does\n",
                        "  not for any Student t). random = \"zero\" gives ",
                        "the typical group."), key, param), call. = FALSE)
  }
  list(gaussian = FALSE, dim = d, pen = pen, theta = th, m = m)
}


#' Nodes and Weights for the Effects a Marginal Prediction Integrates
#'
#' @description
#' A product Gauss-Hermite grid where every prior is Gaussian and the total
#' dimension is at most three; the prior itself, for [marginal_quad()], where
#' the only prior is univariate and not Gaussian; Monte Carlo draws otherwise.
#'
#' @details
#' The grid has \eqn{\min(40, \lfloor 8000^{1/D}\rfloor)} nodes per
#' dimension, so at most 8000 in all. Measured on a random intercept, 20
#' nodes give the marginal mean exactly to 1e-15 under a log link and to
#' 3e-10 and 9e-06 under a logit at prior standard deviations of 1 and 2,
#' where 10000 Monte Carlo draws stop between 1e-2 and 1e-3. The draws are
#' taken from a seed of their own and the caller's random stream is put back
#' afterwards, so the prediction is reproducible and moves nothing a session
#' draws next.
#'
#' @param priors A list of [random_prior()] results.
#' @param nsim The number of draws where the grid is not used.
#'
#' @return A list with `b`, one matrix per prior with a row per node and a
#'   column per coordinate, and `w`, the weights, summing to one.
#'   For one univariate prior that is not Gaussian, a list with `quad`, that
#'   prior.
#'
#' @keywords internal
random_nodes <- function(priors, nsim = 10000L) {
  # ONE UNIVARIATE PRIOR THAT IS NOT GAUSSIAN is integrated exactly, by
  # adaptive quadrature against its own density, rather than by draws: on a
  # logit with a Student t prior 10000 draws were 3e-3 out
  if (length(priors) == 1L && !priors[[1L]]$gaussian &&
      priors[[1L]]$dim == 1L &&
      "parent" %in% S7::prop_names(priors[[1L]]$pen)) {
    return(list(quad = priors[[1L]]))
  }
  D <- sum(vapply(priors, `[[`, integer(1), "dim"))
  gauss <- all(vapply(priors, `[[`, logical(1), "gaussian"))
  if (gauss && D <= 3L) {
    k <- min(40L, floor(8000^(1 / D)))
    q <- gauss_hermite(k)
    grid <- as.matrix(expand.grid(rep(list(seq_len(k)), D)))
    x <- matrix(q$x[grid], nrow(grid), D) * sqrt(2)
    w <- apply(matrix(q$w[grid], nrow(grid), D), 1L, prod) / pi^(D / 2)
    b <- list()
    at <- 0L
    for (pr in priors) {
      cols <- at + seq_len(pr$dim)
      b[[length(b) + 1L]] <- x[, cols, drop = FALSE] %*% t(pr$chol)
      at <- at + pr$dim
    }
    return(list(b = b, w = w))
  }
  if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", old, globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(20260928L)
  b <- lapply(priors, function(pr) {
    if (pr$gaussian) {
      return(matrix(stats::rnorm(nsim * pr$dim), nsim, pr$dim) %*%
               t(pr$chol))
    }
    out <- NULL
    while (is.null(out) || nrow(out) < nsim) {
      v <- penalties7::penalty_draw(pr$pen, pr$theta)
      out <- rbind(out, matrix(v, ncol = pr$dim, byrow = TRUE))
    }
    out[seq_len(nsim), , drop = FALSE]
  })
  list(b = b, w = rep(1 / nsim, nsim))
}


#' Gauss-Hermite Nodes and Weights
#'
#' @description
#' For \eqn{\int f(x) e^{-x^2} dx}, by the eigendecomposition of the Jacobi
#' matrix of the Hermite polynomials (Golub and Welsch 1969).
#'
#' @param k The number of nodes.
#'
#' @return A list with `x` and `w`.
#'
#' @keywords internal
gauss_hermite <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1L, ]^2)
}


#' A Marginal Prediction: the Population Average Over New Groups
#'
#' @description
#' Integrates the quantity asked for over the prior of every random-effect
#' term read at `"marginal"`, at the fixed part of each equation's predictor.
#'
#' @details
#' At each node the predictor of equation \eqn{p} is
#' \eqn{\eta_{0p} + \sum_t Z_t b_t}, with \eqn{Z_t} the term's within-group
#' design ([modelterms7::term_within()]). A parameter is averaged,
#' \eqn{E_b[h^{-1}(\eta_0 + Zb)]}; the mean is averaged too, by the law of
#' total expectation; the variance adds the variance of the conditional mean,
#' \eqn{E_b[\mathrm{Var}(Y\mid b)] + \mathrm{Var}_b(E[Y\mid b])}, which is
#' where a random effect on the scale enters. A predictor, `"link"` or
#' `"link:p"`, is \eqn{\eta_0} itself, the effects having mean zero.
#'
#' @param spec,design The specification at the new rows and its design, with
#'   the marginal terms' columns set to zero.
#' @param ep The predictors at the fixed part, as [statmod_eta()] returns
#'   them.
#' @param mt The marginal terms, rows of [random_modes()].
#' @param nodes What [random_nodes()] returns, in the order of `mt`.
#' @param what What was asked for.
#' @param eta_at A named list of predictors to integrate at in place of
#'   `ep$eta`, for the ends of an interval; `NULL` for `ep$eta`.
#' @param deriv `TRUE` to average the derivative of the inverse link instead
#'   of the inverse link, which is what the delta method of a marginal
#'   parameter multiplies the predictor's standard error by. Read only for a
#'   parameter.
#'
#' @return A numeric vector, or a named list of them for `"parameter"`.
#'
#' @keywords internal
random_marginal <- function(spec, design, ep, mt, nodes, what, eta_at = NULL,
                            deriv = FALSE) {
  inv <- if (deriv) linkfunctions7::dlinkinv else linkfunctions7::linkinv
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  n <- spec@n_obs
  eta0 <- if (is.null(eta_at)) ep$eta else eta_at
  eta0 <- lapply(eta0, function(e) rep_len(as.numeric(e), n))
  Zs <- lapply(seq_len(nrow(mt)), function(i) {
    tm <- spec@terms[[mt$param[i]]][[mt$key[i]]]
    if (!is.null(spec@newdata)) {
      return(modelterms7::term_within(tm, spec@newdata))
    }
    # at the fitting rows the block itself says it: a row is non-zero in its
    # own group's d columns alone, so summing the groups' columns returns
    # the within-group row
    gr <- modelterms7::term_group(tm)
    as.matrix(modelterms7::term_matrix(tm) %*%
                kronecker(rep(1, length(gr$levels)), diag(gr$dim)))
  })
  mom <- predict_moments()
  kind <- if (what %in% params) "param" else if (identical(what, "parameter"))
    "all" else if (what %in% c("mean", "variance", "std_dev")) "moment" else
      stop(sprintf(paste0("'%s' has no marginal prediction here. Ask for a ",
                          "parameter (%s),\n  \"parameter\", \"mean\", ",
                          "\"variance\", \"std_dev\" or a predictor."),
                   what, paste(params, collapse = ", ")), call. = FALSE)
  if (!is.null(nodes$quad)) {
    return(marginal_quad(spec, eta0, Zs[[1L]], mt$param[1L], nodes$quad,
                         what, kind, inv))
  }
  acc <- stats::setNames(lapply(params, function(p) numeric(n)), params)
  m1 <- m2 <- v1 <- numeric(n)
  for (k in seq_along(nodes$w)) {
    eta <- eta0
    for (i in seq_len(nrow(mt))) {
      p <- mt$param[i]
      eta[[p]] <- eta[[p]] + as.numeric(Zs[[i]] %*% nodes$b[[i]][k, ])
    }
    th <- stats::setNames(lapply(params, function(p)
      as.numeric(inv(links[[p]], eta[[p]]))), params)
    w <- nodes$w[k]
    if (kind %in% c("param", "all")) {
      for (p in params) acc[[p]] <- acc[[p]] + w * th[[p]]
    } else {
      mu <- rep_len(mom[["mean"]](spec@distrib, th), n)
      m1 <- m1 + w * mu
      m2 <- m2 + w * mu^2
      if (what != "mean") {
        v1 <- v1 + w * rep_len(mom[["variance"]](spec@distrib, th), n)
      }
    }
  }
  # A PREDICTOR WHOSE DOMAIN IS NOT THE WHOLE LINE (square root, inverse)
  # leaves it with positive probability under a Gaussian effect, and there the
  # average over new groups does not exist: E[1/(eta0 + u)] is infinite for a
  # Gaussian u. Where that probability is negligible the average over the
  # nodes is the meaningful number; above 1e-8 it is reported as NA.
  checked <- switch(kind, param = intersect(what, mt$param),
                    unique(mt$param))
  for (p in checked) {
    bad <- marginal_out_of_domain(links[[p]], eta0[[p]], Zs, nodes,
                                  which(mt$param == p))
    if (!any(bad)) next
    if (kind == "moment") {
      m1[bad] <- m2[bad] <- v1[bad] <- NA_real_
    } else {
      acc[[p]][bad] <- NA_real_
    }
  }
  switch(kind,
         param = acc[[what]],
         all = acc,
         moment = switch(what,
                         mean = m1,
                         variance = v1 + m2 - m1^2,
                         std_dev = sqrt(v1 + m2 - m1^2)))
}


#' Rows Where a New Group's Predictor Leaves the Domain of Its Link
#'
#' @description
#' A predictor whose domain is not the whole line (the square root and
#' inverse links) leaves it with positive probability under a Gaussian
#' effect, and there the average over new groups does not exist:
#' \eqn{E[1/(\eta_0 + u)]} is infinite for a Gaussian \eqn{u}. Where that
#' probability is below 1e-8 the average over the nodes is the meaningful
#' number; above it the row is reported as `NA`, with a warning.
#'
#' @details
#' The variance of the effects' contribution at each row is read off the
#' nodes, \eqn{\sum_k w_k (z^\top b_k)^2}, which a Gauss-Hermite grid gives
#' exactly, the effects having mean zero.
#'
#' @param link The parameter's link.
#' @param eta0 The predictor at the fixed part, one value a row.
#' @param Zs The within-group rows of the marginal terms.
#' @param nodes What [random_nodes()] returns.
#' @param idx The marginal terms in this parameter's equation.
#'
#' @return A logical vector, `TRUE` at a row to report as `NA`.
#'
#' @keywords internal
marginal_out_of_domain <- function(link, eta0, Zs, nodes, idx) {
  n <- length(eta0)
  eb <- linkfunctions7::eta_bounds(link)
  if (!any(is.finite(eb)) || is.null(nodes$w)) return(rep(FALSE, n))
  v <- numeric(n)
  for (k in seq_along(nodes$w)) {
    add <- numeric(n)
    for (i in idx) add <- add + as.numeric(Zs[[i]] %*% nodes$b[[i]][k, ])
    v <- v + nodes$w[k] * add^2
  }
  s <- sqrt(v)
  pr <- ifelse(s > 0,
               stats::pnorm((eb[1L] - eta0) / s) +
                 stats::pnorm((eta0 - eb[2L]) / s),
               as.numeric(eta0 <= eb[1L] | eta0 >= eb[2L]))
  bad <- pr > 1e-8
  if (any(bad)) {
    warning(sprintf(paste0("A new group's predictor leaves the domain of the ",
                           "%s link with probability up\n  to %.2g, so its ",
                           "average over new groups does not exist there; ",
                           "those rows are NA."), link@link_name,
                    max(pr[bad])), call. = FALSE)
  }
  bad
}


#' A Marginal Prediction by Adaptive Quadrature Over One Univariate Prior
#'
#' @description
#' The average over a new group's effect when the only effect set aside is
#' one coordinate with a prior that is not Gaussian: each quantity is
#' \eqn{\int q(\eta_0 + z u)\, f_b(u)\, du}, with \eqn{f_b} the prior's own
#' density, computed by [numericals7::quad_vec()] over the rows at once.
#'
#' @details
#' Measured on a Bernoulli with a logit link and a Student t prior, 10000
#' Monte Carlo draws were 3e-3 out against `integrate()`; this is within the
#' default tolerances of `quad_vec()`.
#'
#' @param spec The specification at the rows predicted.
#' @param eta0 The predictors at the fixed part, a named list.
#' @param Z The within-group column of the term, one value a row.
#' @param param The parameter whose equation holds the term.
#' @param prior What [random_prior()] returns for the term.
#' @param what,kind What was asked for and its kind, as in
#'   [random_marginal()].
#' @param inv The inverse link or its derivative.
#'
#' @return As [random_marginal()].
#'
#' @keywords internal
marginal_quad <- function(spec, eta0, Z, param, prior, what, kind, inv) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  n <- spec@n_obs
  Z <- rep_len(as.numeric(Z), n)
  par <- prior$pen@parent
  th0 <- prior$theta
  mom <- predict_moments()
  theta_at <- function(x, i) {
    stats::setNames(lapply(params, function(p) {
      e <- eta0[[p]][i]
      if (identical(p, param)) e <- e + Z[i] * x
      as.numeric(inv(links[[p]], e))
    }), params)
  }
  integral <- function(q) {
    f <- function(x, i) {
      xv <- as.vector(x)
      ii <- rep(i, length.out = length(xv))
      dens <- as.numeric(distributions7::distrib_pdf(par, xv, th0))
      dens * q(theta_at(xv, ii))
    }
    numericals7::quad_vec(f, -Inf, rep(Inf, n))
  }
  if (kind == "param") return(integral(function(th) th[[what]]))
  if (kind == "all") {
    return(stats::setNames(lapply(params, function(p)
      integral(function(th) th[[p]])), params))
  }
  m1 <- integral(function(th) mom[["mean"]](spec@distrib, th))
  if (what == "mean") return(m1)
  m2 <- integral(function(th) mom[["mean"]](spec@distrib, th)^2)
  v1 <- integral(function(th) mom[["variance"]](spec@distrib, th))
  switch(what, variance = v1 + m2 - m1^2, std_dev = sqrt(v1 + m2 - m1^2))
}
