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
  em <- engine_outer_pen(spec, design, coef, hyper)
  if (is.null(em)) return(NULL)
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
#' @details
#' `engine_outer_pen()` returns the model and the penalty as the outer
#' context reads them, `list(model, pen)`, or `NULL` where the point is not
#' covered: [engine_model()] on the Cholesky route, with the pieces that the
#' context reads beside the penalty's units (the positions of a t prior's
#' hyperparameters, an additive penalty's blocks, a covariance chart).
#' @keywords internal
engine_outer_pen <- function(spec, design, coef, hyper) {
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
  em
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
engine_flags$outer_search <- TRUE
engine_flags$outer_point <- TRUE
engine_flags$outer_pen <- TRUE

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

#' @rdname engine_outer_at
#' @param optimizer An optimizer of optimizers7.
#' @param nv The length of the outer vector.
#' @details
#' `outer_search_settings()` returns the settings of a descent method for
#' the compiled outer search ([optimizers7::descent_settings()]), or `NULL`
#' where the search runs in R: an optimizer that does not run the descent
#' loop, and a stopping rule without a compiled form. The switch
#' `engine_flags$outer_search` turns the compiled search on, `engine_flags$outer_point`
#' the point computed in compiled code within it, and `engine_flags$outer_pen`
#' the penalty built in compiled code at each point; all three are on.
#' @keywords internal
outer_search_settings <- function(optimizer, nv) {
  s <- optimizers7::descent_settings(optimizer, nv)
  if (is.null(s) || is.null(s$rule_spec)) return(NULL)
  s
}

#' @rdname engine_outer_at
#' @param beta The stacked coefficients at the start.
#' @param how The solver of the inner fit, as for [engine_model()].
#' @details
#' `engine_pen_spec()` describes the penalty to the compiled construction of
#' the outer search (`engine_pen_build()`): per penalized unit its kind, its
#' positions, its description by [engine_penalty_units()] for the inner fit
#' and for the outer context, its hyperparameters in the order of the
#' penalty's parameters, and the static pieces of its kind; and the members of
#' the outer index with the links of its rows. It returns `NULL` where a unit
#' is not one that the compiled construction writes (kinds 1 to 4, a kind 4
#' unit on a log-Cholesky, compound symmetry or AR(1) chart, no map on kinds 3
#' and 4), or where a link has no compiled route.
#' @keywords internal
engine_pen_spec <- function(spec, design, idx, hyper, beta, how) {
  cf <- NULL
  em <- tryCatch({
    obj <- statmod_objective(spec, hyper, design, FALSE, NULL)
    cf <- obj$split(beta)
    engine_model(spec, design, cf, hyper, how)
  }, error = function(e) NULL)
  eo <- if (is.null(cf)) NULL else engine_outer_pen(spec, design, cf, hyper)
  if (is.null(em) || is.null(eo)) return(NULL)
  units <- statmod_penalized(spec, design)
  link_ref <- function(l) {
    r <- linkfunctions7::link_scalar_route(l)
    if (is.null(r)) return(NULL)
    list(name = r$name, par = as.numeric(r$par), bounds = as.numeric(l@link_bounds))
  }
  out <- vector("list", length(units))
  for (k in seq_along(units)) {
    u <- units[[k]]
    pen <- u$penalty
    ul <- em$pen$units[[k]]
    ulo <- eo$pen$units[[k]]
    th <- as.list(hyper[[u$param]][[u$key]])
    if (!identical(names(th), pen@params) ||
        !all(vapply(th, function(x) is.numeric(x) && length(x) == 1L, logical(1)))) {
      return(NULL)
    }
    s <- list(kind = ul$kind, idx = as.integer(ul$idx), ul = ul, ulo = ulo,
              theta = as.numeric(unlist(th)))
    if (identical(ul$kind, 3L)) {
      if (!is.null(ul$map)) return(NULL)
      s$si <- match("sigma", pen@params)
      if (is.na(s$si) || !identical(ul$sigma, th$sigma)) return(NULL)
    } else if (identical(ul$kind, 10L)) {
      if (!is.null(ul$map)) return(NULL)
      ms <- engine_mvt_static(u, ul)
      if (is.null(ms)) return(NULL)
      s$mvt <- ms
    } else if (identical(ul$kind, 4L) && identical(ulo$chart$code, 9L)) {
      # dr_prod with correlation_matrix (D49)
      if (!is.null(ul$map)) return(NULL)
      full <- distributions7::fixed_full_theta(pen@parent, th)
      pc <- distributions7::mvg_pieces(pen@parent@parent_distrib, full)
      cs <- engine_chart_static(pc$s, isTRUE(ul$inverted))
      if (is.null(cs) || !identical(cs$code, 4L) ||
          !identical(as.numeric(ulo$chart$eta), s$theta)) {
        return(NULL)
      }
      s$chart7 <- cs
    } else if (identical(ul$kind, 4L)) {
      ch <- ulo$chart
      if (!is.null(ul$map) || is.null(ch) || !ch$code %in% 1:3 ||
          !identical(as.numeric(ch$eta), s$theta)) {
        return(NULL)
      }
      full <- distributions7::fixed_full_theta(pen@parent, th)
      pc <- distributions7::mvg_pieces(pen@parent@parent_distrib, full)
      s$code <- ch$code
      s$pd <- as.integer(pen@block)
      s$inverted <- isTRUE(ul$inverted)
      if (ch$code == 1L) {
        pos <- pc$s@param_params$positions
        s$lc_row <- as.integer(pos$row)
        s$lc_col <- as.integer(pos$col)
        s$lc_diag <- as.logical(pos$on_diagonal)
      } else {
        lk <- pc$s@param_params
        s$link_scale <- link_ref(lk$link_scale)
        s$link_rho <- link_ref(lk$link_rho)
        if (is.null(s$link_scale) || is.null(s$link_rho)) return(NULL)
      }
    } else if (ul$kind %in% 5:8) {
      second <- switch(as.character(ul$kind), "5" = NA_character_, "6" = "alpha",
                       "7" = "a", "8" = "gamma")
      s$tl <- match("lambda", pen@params, nomatch = 0L)
      s$tp2 <- if (is.na(second)) 0L else match(second, pen@params, nomatch = 0L)
      s$p20 <- if (identical(ul$kind, 6L)) {
        as.numeric(distributions7::fixed_full_theta(pen@parent, th)[[3]])
      } else if (ul$kind >= 7L) as.numeric(ul$par) else 0
      if (!s$tl || !identical(as.numeric(ul$lambda), th$lambda)) return(NULL)
    } else if (identical(ul$kind, 9L)) {
      s$tsig <- match("sigma", pen@params, nomatch = 0L)
      s$tnu <- match("nu", pen@params, nomatch = 0L)
      if ((s$tsig && !identical(ul$sigma, th$sigma)) ||
          (s$tnu && !identical(ul$nu, th$nu))) {
        return(NULL)
      }
    } else if (!ul$kind %in% 1:2) {
      return(NULL)
    }
    out[[k]] <- s
  }
  mb <- engine_outer_members(spec, design, idx, hyper)
  if (is.null(mb)) return(NULL)
  if (!length(mb$mrow)) mb$mrow <- mb$munit <- mb$mname <- integer(0)
  links <- lapply(attr(idx, "links"), link_ref)
  if (any(vapply(links, is.null, logical(1)))) return(NULL)
  list(p = length(beta), sparse = isTRUE(em$model$sparse),
       need_C = how %in% c("qr", "svd"), units = out,
       members = list(mrow = as.integer(mb$mrow), munit = as.integer(mb$munit),
                      mname = as.integer(mb$mname), links = links))
}

#' @rdname engine_outer_at
#' @param s A matrix parameter of parameters7 (a covariance chart).
#' @param inverted Whether the chart is of the precision.
#' @details
#' `engine_chart_static()` describes a covariance chart to the compiled
#' construction (`engine_charts.h`): its code (1 log-Cholesky, 2 compound
#' symmetry, 3 AR(1), 4 `dr_prod()` over a `correlation_matrix()`), its
#' dimension and number of free values, the orientation, the positions of the
#' log-Cholesky factor or of the angles (0-based), and the links of the
#' scales and of the correlation. It returns `NULL` for any other chart, and
#' for a `dr_prod()` whose links are not the log and the doubly bounded link.
#' `engine_mvt_static()` describes a multivariate t prior: its chart, its
#' location, the role of each hyperparameter (the chart's free value, 0-based,
#' or -1 for `nu`), and for an odd dimension the coefficients of the series of
#' `mvt_cp_d()`.
#' @keywords internal
engine_chart_static <- function(s, inverted = FALSE) {
  cls <- attr(S7::S7_class(s), "name")
  pp <- s@param_params
  out <- list(code = switch(cls, LogCholeskyParam = 1L, CompoundSymmetryParam = 2L,
                            Ar1Param = 3L, DrProdParam = 4L, NA_integer_),
              p = as.integer(s@dimension), K = as.integer(s@n_free),
              inverted = isTRUE(inverted))
  if (is.na(out$code)) return(NULL)
  lspec <- function(lk, pre) {
    r <- linkfunctions7::link_scalar_route(lk)
    if (is.null(r)) return(NULL)
    stats::setNames(list(r$name, as.numeric(r$par), as.numeric(lk@link_bounds[1]),
                         as.numeric(lk@link_bounds[2])),
                    paste0(pre, c("_name", "_par", "_lwr", "_upr")))
  }
  if (cls == "LogCholeskyParam") {
    pos <- pp$positions
    out$row <- as.integer(pos$row - 1L)
    out$col <- as.integer(pos$col - 1L)
    out$ondiag <- as.logical(pos$on_diagonal)
  } else if (cls %in% c("CompoundSymmetryParam", "Ar1Param")) {
    a <- lspec(pp$link_scale, "ls")
    b <- lspec(pp$link_rho, "lr")
    if (is.null(a) || is.null(b)) return(NULL)
    out <- c(out, a, b)
  } else {
    cor <- pp$cor
    if (!identical(attr(S7::S7_class(cor), "name"), "CorrelationParam") || out$p < 2L) {
      return(NULL)
    }
    a <- lspec(pp$link, "ls")
    b <- lspec(cor@param_params$link, "lr")
    if (is.null(a) || is.null(b) || !identical(a$ls_name, "LogLink") ||
        !identical(b$lr_name, "DoublyBoundedLink")) {
      return(NULL)
    }
    out <- c(out, a, b)
    out$arow <- as.integer(cor@param_params$row - 1L)
    out$acol <- as.integer(cor@param_params$col - 1L)
  }
  out
}

#' @rdname engine_outer_at
#' @param ul The unit's description by [engine_penalty_units()].
#' @keywords internal
engine_mvt_static <- function(u, ul) {
  pen <- u$penalty
  parent <- pen@parent@parent_distrib
  p <- as.integer(pen@block)
  if (p < 2L) return(NULL)
  ch <- engine_chart_static(parent@param, FALSE)
  if (is.null(ch)) return(NULL)
  par <- pen@params
  chart <- setdiff(par, "nu")
  role <- vapply(par, function(nm) if (nm == "nu") -1L else match(nm, chart) - 1L, 0L)
  if (length(chart) != ch$K) return(NULL)
  cf <- list(numeric(0), numeric(0), numeric(0))
  if (p %% 2L == 1L) {
    g <- distributions7:::mvt_nu_coef(p)
    n <- seq_len(length(g) - 1L)
    cf[[1]] <- -g[n + 1L] * factorial(n - 1L)
    for (k in 1:2) cf[[k + 1]] <- (-1)^(k - 1L) * g[n + 1L] * factorial(n + k - 1L)
  }
  list(chart = ch, mu = as.numeric(ul$mu), role = unname(role),
       cf0 = cf[[1]], cf1 = cf[[2]], cf2 = cf[[3]])
}
