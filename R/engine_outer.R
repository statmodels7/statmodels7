#' The Outer Criterion on the Compiled Route
#'
#' @description
#' `engine_outer_at()` returns the compiled context of an outer point, from
#' which [statmod_marginal()], [statmod_marginal_grad()] and
#' [statmod_marginal_hess()] read the criterion, its gradient and its Hessian
#' on the models that the compiled route covers, or `NULL` where the R route
#' computes them.
#'
#' @details
#' A point is covered when `spec@engine` is `"auto"`, the method is
#' [reml()] or [ml()] on the observed information, the model is one that
#' [engine_model()] describes, the family carries its derivatives to the
#' fourth order in the registry, and every penalized unit is a quadratic
#' penalty or a gaussian prior with a fixed mean and no map. The compiled
#' context computes, at the mode, the derivatives of the log-density to the
#' fourth order on the link scale, the observed information, the penalized
#' information pinned at [pinned_coords()] with its inverse, the matrix of
#' the traces and the leverages, as the evaluation context of
#' [outer_context()] holds them. Where the R route would take a branch that
#' the compiled context does not reproduce (a curvature that is not finite,
#' a Cholesky factor that fails or whose condition is within a factor of 100
#' of the acceptance threshold of [pd_factor()], and therefore the
#' eigendecomposition and the retry over [deficient_coords()]), the context
#' reports it and the R route computes the point. The context is stored in
#' `ctx`, keyed by the integrated basis, so that the criterion, the gradient
#' and the Hessian at one point read one context.
#'
#' `engine_outer_members()` describes the outer index to the compiled
#' context: for each member line of [index_members()] its row, the position
#' of its penalized unit in [statmod_penalized()], and the derivative of its
#' hyperparameter's inverse link at its value; and, for each row, the first
#' and second derivatives of the inverse link at the first member's value,
#' which the Hessian reads.
#'
#' @param ctx An evaluation context from [outer_context()], or `NULL`.
#' @param spec,design,coef,hyper,method,basis As for [statmod_marginal()].
#' @param idx The outer index.
#'
#' @return `engine_outer_at()`: an external pointer, or `NULL`.
#'   `engine_outer_members()`: a list with `mrow`, `munit`, `mchain`, `h1`
#'   and `h2`, or `NULL` where a member is not covered.
#'
#' @seealso [engine_inner_fit()], [statmod_marginal()]
#'
#' @keywords internal
engine_outer_at <- function(ctx, spec, design, coef, hyper, method,
                            basis = NULL) {
  if (!identical(spec@engine, "auto")) return(NULL)
  if (!method@kind %in% c("reml", "ml")) return(NULL)
  if (info_expected(method@hessian)) return(NULL)
  if (!is.null(ctx) && exists("engine_outer", envir = ctx, inherits = FALSE)) {
    e <- get("engine_outer", envir = ctx)
    if (identical(e$basis, basis)) return(e$handle)
  }
  h <- engine_outer_build(spec, design, coef, hyper, basis)
  if (!is.null(ctx)) {
    assign("engine_outer", list(basis = basis, handle = h), envir = ctx)
  }
  h
}

#' @rdname engine_outer_at
#' @keywords internal
engine_outer_build <- function(spec, design, coef, hyper, basis) {
  kc <- kinked_coords(spec, design)
  if (length(kc) && !is.null(basis)) return(NULL)
  em <- tryCatch(engine_model(spec, design, coef, hyper, "chol"),
                 error = function(e) NULL)
  if (is.null(em) || em$covered < 4L) return(NULL)
  units <- statmod_penalized(spec, design)
  for (k in seq_along(em$pen$units)) {
    u <- em$pen$units[[k]]
    # a t prior (kind 9): the positions of its hyperparameters in (sigma, nu),
    # and whether its Hessian is read as fixed at this point
    # a multivariate t prior (kind 10): its derivatives in the
    # hyperparameters at the mode, the dS blocks group by group, the chart's
    # derivatives of Sinv and of N
    # (given at the first gradient or Hessian, engine_outer_complete())
    if (identical(u$kind, 10L)) next
    if (identical(u$kind, 9L)) {
      pen <- units[[k]]$penalty
      th <- as.list(hyper[[units[[k]]$param]][[units[[k]]$key]])
      tf <- match(pen@params, c("sigma", "nu"))
      if (anyNA(tf)) return(NULL)
      em$pen$units[[k]]$tfree <- as.integer(tf)
      em$pen$units[[k]]$bq <- isTRUE(penalties7::beta_quadratic(pen, th))
      next
    }
    if (identical(u$kind, 1L) || u$kind >= 5L) next
    if (identical(u$kind, 3L)) next
    if (identical(u$kind, 2L)) {
      pen <- units[[k]]$penalty
      m <- if (length(pen@blocks)) as.integer(pen@blocks[[1L]]) else 1L
      em$pen$units[[k]]$add <- list(
        m = m, rank = as.numeric(pen@rank), tol = as.numeric(pen@tol),
        block_mats = if (m > 1L) lapply(pen@block_mats, as_dense) else list())
      next
    }
    if (!identical(u$kind, 4L)) return(NULL)
    ch <- engine_outer_chart(units[[k]], hyper)
    if (is.null(ch)) return(NULL)
    em$pen$units[[k]]$chart <- ch
  }
  beta <- unlist(coef[spec@distrib@params], use.names = FALSE)
  pc <- as.integer(pinned_coords(spec, design))
  # beside a kinked penalty, the criterion's matrix is pinned at every kinked
  # coordinate (laplace_pinned()) and the mode's at the kinked ones at zero
  hold <- if (length(kc)) as.integer(laplace_pinned(spec, design)) else pc
  hm <- if (length(kc)) as.integer(c(pc, zero_kinked(spec, design, coef))) else NULL
  B <- if (is.null(basis)) NULL else as_dense(basis)
  r <- engine_outer_new(em$model, em$pen, as.numeric(beta), hold, B, 2L, hm)
  if (!identical(r$status, 0L)) return(NULL)
  r$ctx
}

#' @rdname engine_outer_at
#' @keywords internal
engine_outer_members <- function(spec, design, idx, hyper) {
  mem <- index_members(idx)
  links <- attr(idx, "links")
  units <- statmod_penalized(spec, design)
  keys <- vapply(units, function(u) paste(u$param, u$key, sep = "\r"),
                 character(1))
  nm <- nrow(mem)
  munit <- integer(nm)
  mname <- integer(nm)
  mchain <- numeric(nm)
  for (i in seq_len(nm)) {
    k <- match(paste(mem$parameter[i], mem$term[i], sep = "\r"), keys)
    if (is.na(k)) return(NULL)
    pen <- units[[k]]$penalty
    h <- match(mem$name[i], pen@params)
    if (is.na(h)) return(NULL)
    mname[i] <- h
    munit[i] <- k
    r <- mem$row[i]
    v <- hyper[[mem$parameter[i]]][[mem$term[i]]][[mem$name[i]]]
    e <- linkfunctions7::linkfun(links[[r]], v)
    mchain[i] <- linkfunctions7::dlinkinv(links[[r]], e)
  }
  nh <- nrow(idx)
  h1 <- numeric(nh)
  h2 <- numeric(nh)
  for (r in seq_len(nh)) {
    v <- hyper[[idx$parameter[r]]][[idx$term[r]]][[idx$name[r]]]
    e <- linkfunctions7::linkfun(links[[r]], v)
    h1[r] <- linkfunctions7::dlinkinv(links[[r]], e)
    h2[r] <- linkfunctions7::d2linkinv(links[[r]], e)
  }
  list(mrow = as.integer(mem$row), munit = munit, mname = mname,
       mchain = mchain, h1 = h1, h2 = h2, nh = nh)
}

#' @rdname engine_outer_at
#' @param x A vector of the length of the stacked coefficients.
#' @details
#' `engine_ctx_solve()` returns \eqn{K^{-1}x} with \eqn{K} the mode's pinned
#' penalized information, from the compiled context that `ctx` already holds,
#' or `NULL` where `ctx` holds none or the context does not cover the point,
#' in which case the caller reads the inverse of [ctx_penalized()].
#' [inner_mode_error()] and [criterion_resolution()] read it. The switch
#' `engine_flags$mode_inverse` turns it off.
#' @keywords internal
engine_ctx_solve <- function(ctx, x) {
  if (!isTRUE(engine_flags$mode_inverse)) return(NULL)
  if (is.null(ctx) || !exists("engine_outer", envir = ctx, inherits = FALSE)) {
    return(NULL)
  }
  h <- get("engine_outer", envir = ctx)$handle
  if (is.null(h)) return(NULL)
  r <- engine_outer_kinv(h, as.numeric(x))
  if (!identical(r$status, 0L)) return(NULL)
  r$x
}

#' @rdname engine_outer_at
#' @format `engine_flags`: an environment of switches of the compiled route.
#' @keywords internal
engine_flags <- new.env(parent = emptyenv())
engine_flags$mode_inverse <- TRUE

#' @rdname engine_outer_at
#' @param eo An outer context of [engine_outer_at()].
#' @details
#' `engine_outer_complete()` gives the context the derivatives of each
#' multivariate t prior in its hyperparameters ([engine_outer_mvt()]), which
#' the criterion does not read, before the first gradient or Hessian. It
#' returns the context, or `NULL` where a prior is not covered.
#' @keywords internal
engine_outer_complete <- function(eo, spec, design, coef, hyper) {
  need <- engine_outer_mvt_need(eo)
  if (!length(need)) return(eo)
  units <- statmod_penalized(spec, design)
  for (k in need) {
    ob <- engine_outer_mvt(units[[k]], NULL, coef, spec, hyper)
    if (is.null(ob)) return(NULL)
    engine_outer_mvt_set(eo, k, ob)
  }
  eo
}

#' @rdname engine_outer_at
#' @param un The unit's description by [engine_penalty_units()].
#' @details
#' `engine_outer_mvt()` describes a multivariate Student t prior to the
#' compiled context at the mode: the derivatives of the penalty in its
#' hyperparameters ([penalties7::penalty_grad_theta()],
#' [penalties7::penalty_hess_theta()], [penalties7::penalty_cross()],
#' [penalties7::penalty_dcross()]), the blocks of
#' [penalties7::penalty_dhessian()] and [penalties7::penalty_d2hessian()]
#' group by group, whether [penalties7::beta_quadratic()] holds, and per
#' hyperparameter the derivative of the inverse scale matrix (zero for `nu`)
#' and of `nu + p` (one for `nu`), from which the compiled context forms the
#' movement of the penalty's Hessian with the coefficients. It returns
#' `NULL` where the chart's free values do not match the hyperparameters.
#' @keywords internal
engine_outer_mvt <- function(u, un, coef, spec, hyper) {
  pen <- u$penalty
  th <- as.list(hyper[[u$param]][[u$key]])
  b <- unit_beta(u, coef, spec@distrib@params)
  par <- pen@params
  H <- length(par)
  p <- as.integer(pen@block)
  G <- length(b) %/% p
  full <- distributions7::fixed_full_theta(pen@parent, th)
  pc <- distributions7:::mvt_pieces(pen@parent@parent_distrib, full, derivs = TRUE)
  chart <- setdiff(par, "nu")
  if (length(pc$ai) != length(chart)) return(NULL)
  key <- function(a, b) {
    if (match(a, par) > match(b, par)) { t <- a; a <- b; b <- t }
    paste(a, b, sep = "_")
  }
  pairs <- expand.grid(k = seq_len(H), h = seq_len(H))
  gt <- penalties7::penalty_grad_theta(pen, b, th)
  ht <- penalties7::penalty_hess_theta(pen, b, th)
  cr <- penalties7::penalty_cross(pen, b, th)
  dc <- penalties7::penalty_dcross(pen, b, th)
  ds <- penalties7::penalty_dhessian(pen, b, th)
  d2 <- penalties7::penalty_d2hessian(pen, b, th)
  blocks <- function(M) {
    M <- as_dense(M)
    unlist(lapply(seq_len(G), function(g) {
      ix <- (g - 1L) * p + seq_len(p)
      M[ix, ix]
    }))
  }
  pk <- vapply(seq_len(nrow(pairs)), function(i) key(par[pairs$h[i]], par[pairs$k[i]]),
               character(1))
  list(
    H = H, bq = isTRUE(penalties7::beta_quadratic(pen, th)),
    gt = as.numeric(unlist(gt[par])),
    rho2 = as.numeric(unlist(ht[pk])),
    c = as.numeric(unlist(cr[par])),
    c2 = as.numeric(unlist(dc[pk])),
    dS = unlist(lapply(ds[par], blocks)),
    d2S = unlist(lapply(d2[pk], blocks)),
    ai = lapply(par, function(nm) {
      if (nm == "nu") matrix(0, p, p) else unname(as.matrix(pc$ai[[match(nm, chart)]]))
    }),
    Na = as.numeric(par == "nu"))
}

#' @rdname engine_outer_at
#' @param u A penalized unit of [statmod_penalized()] whose penalty is a
#'   multivariate gaussian prior.
#' @details
#' `engine_outer_chart()` describes the covariance chart of a multivariate
#' gaussian prior to the compiled context: a code (1 for
#' [parameters7::log_cholesky()], of the covariance or of the precision, 2
#' for [parameters7::compound_symmetry()], 3 for [parameters7::ar1()]), the
#' free values of the chart in its own order, which is the order of the
#' penalty's hyperparameters, and, for the last two, the inverse links of
#' the scale and of the correlation with their first two derivatives at the
#' free values. It returns `NULL` for any other chart.
#' @keywords internal
engine_outer_chart <- function(u, hyper) {
  pen <- u$penalty
  par <- pen@parent
  th <- as.list(hyper[[u$param]][[u$key]])
  full <- distributions7::fixed_full_theta(par, th)
  pc <- distributions7::mvg_pieces(par@parent_distrib, full)
  cls <- attr(S7::S7_class(pc$s), "name")
  eta <- as.numeric(pc$eta)
  if (length(eta) != length(pen@params)) return(NULL)
  code <- switch(cls, LogCholeskyParam = 1L, CompoundSymmetryParam = 2L,
                 Ar1Param = 3L, 9L)
  if (code == 9L) {
    # a chart that the compiled context does not write: the derivatives of
    # Sigma^-1 and of log|Sigma| in its free values, from mv_matrix_pieces(),
    # the second ones as full K x K arrays (t * K + s)
    m <- distributions7:::mv_matrix_pieces(pc$s, eta, isTRUE(par@parent_distrib@inverted),
                                           TRUE, TRUE)
    nm <- names(m$ai)
    K <- length(nm)
    if (K != length(eta) || is.null(m$ai2) || is.null(m$d2logdet)) return(NULL)
    pair <- function(a, b) {
      k <- paste(a, b, sep = ":")
      if (!k %in% names(m$ai2)) k <- paste(b, a, sep = ":")
      k
    }
    keys <- unlist(lapply(nm, function(a) vapply(nm, function(b) pair(a, b), "")))
    if (!all(keys %in% names(m$ai2)) || !all(keys %in% names(m$d2logdet))) return(NULL)
    return(list(code = 9L, eta = eta, h = numeric(0),
                ai = lapply(m$ai, function(x) unname(as.matrix(x))),
                ai2 = lapply(m$ai2[keys], function(x) unname(as.matrix(x))),
                dld = as.numeric(m$dlogdet),
                d2ld = as.numeric(m$d2logdet[keys])))
  }
  h <- numeric(0)
  if (code > 1L) {
    lk <- pc$s@param_params
    hs <- function(l, e) c(linkfunctions7::linkinv(l, e),
                           linkfunctions7::dlinkinv(l, e),
                           linkfunctions7::d2linkinv(l, e))
    h <- c(hs(lk$link_scale, eta[1L]), hs(lk$link_rho, eta[2L]))
  }
  list(code = code, eta = eta, h = h)
}
