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
  if (length(kinked_coords(spec, design))) return(NULL)
  em <- tryCatch(engine_model(spec, design, coef, hyper, "chol"),
                 error = function(e) NULL)
  if (is.null(em) || em$covered < 4L) return(NULL)
  units <- statmod_penalized(spec, design)
  for (k in seq_along(em$pen$units)) {
    u <- em$pen$units[[k]]
    if (identical(u$kind, 1L)) next
    if (identical(u$kind, 3L) && is.null(u$map)) next
    if (identical(u$kind, 2L)) {
      pen <- units[[k]]$penalty
      m <- if (length(pen@blocks)) as.integer(pen@blocks[[1L]]) else 1L
      em$pen$units[[k]]$add <- list(
        m = m, rank = as.numeric(pen@rank), tol = as.numeric(pen@tol),
        block_mats = if (m > 1L) lapply(pen@block_mats, as_dense) else list())
      next
    }
    if (!identical(u$kind, 4L) || !is.null(u$map)) return(NULL)
    ch <- engine_outer_chart(units[[k]], hyper)
    if (is.null(ch)) return(NULL)
    em$pen$units[[k]]$chart <- ch
  }
  beta <- unlist(coef[spec@distrib@params], use.names = FALSE)
  hold <- as.integer(pinned_coords(spec, design))
  B <- if (is.null(basis)) NULL else as_dense(basis)
  r <- engine_outer_new(em$model, em$pen, as.numeric(beta), hold, B, 2L)
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
                 Ar1Param = 3L, 0L)
  if (code == 0L) return(NULL)
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
