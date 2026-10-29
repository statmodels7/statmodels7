# The approximate methods of statmod_profile(): a continuation along the
# profile and a local expansion of the trace of the refitted coefficients.
# Both follow the coefficients that are not held as a function of the held
# ones, and evaluate the log-likelihood and the penalty exactly at the
# coefficients they predict, so that the penalized log-likelihood of a row is
# a lower bound of the profile.


#' @rdname profile_continuation
#' @aliases ProfileContinuation
ProfileContinuation <- S7::new_class(
  "ProfileContinuation",
  parent = profile_method,
  properties = list(order = S7::class_integer, maxit = S7::class_numeric),
  validator = function(self) {
    if (length(self@order) != 1L || !self@order %in% 1:2) {
      return("'order' must be 1 or 2.")
    }
    if (length(self@maxit) != 1L || is.na(self@maxit) || self@maxit < 1) {
      return("'maxit' must be a single number of at least 1, or Inf.")
    }
    NULL
  }
)


#' @rdname profile_taylor
#' @aliases ProfileTaylor
ProfileTaylor <- S7::new_class(
  "ProfileTaylor",
  parent = profile_method,
  properties = list(order = S7::class_integer),
  validator = function(self) {
    if (length(self@order) != 1L || !self@order %in% 2:5) {
      return("'order' must be an integer from 2 to 5.")
    }
    NULL
  }
)


#' @title Profile by Continuation
#' @name profile_continuation
#'
#' @description
#' Selects a method of [statmod_profile()] that follows the refitted
#' coefficients from one value of the grid to the next. At each value a
#' predictor extrapolates the coefficients that are not held from the previous
#' point, and a Newton corrector brings them to the optimum with the held
#' coefficients fixed.
#'
#' @details
#' The predictor of order 1 is the tangent of the trace of the refitted
#' coefficients, \eqn{d\hat\eta/d\theta = -H_{\eta\eta}^{-1} H_{\eta\theta}},
#' where \eqn{\theta} are the held coefficients, \eqn{\eta} the others and
#' \eqn{H} the Hessian of the objective of the fit (the negative penalized
#' log-likelihood) at the previous point. Order 2 adds the curvature of the
#' trace, obtained from the second derivative along the step of the gradient
#' of the objective with the trace developed to order 1, by one central
#' stencil. The corrector is the refit of [profile_grid()] started from the
#' prediction. With `maxit = Inf` it runs to convergence, and the rows equal
#' those of [profile_grid()] within the tolerance of the fit. A finite `maxit`
#' limits the passes of the corrector, and the rows are then approximations
#' whose `converged` is `FALSE` where the limit stopped the corrector.
#'
#' If the corrector signals an error, returns a non-finite objective or, with
#' `maxit = Inf`, does not converge, the step is halved: a point halfway
#' between the previous point and the requested one is computed first, and the
#' requested point is attempted again from it. The intermediate points are not
#' part of the output.
#'
#' Coordinates that the outer criterion estimates (the unpenalized
#' coefficients of the parameters other than the first under [reml()] or
#' [ml()]) are re-estimated by the criterion at every point of the exact
#' method, and with `hyper = "refit"` so are the hyperparameters that `which`
#' does not name. They are followed through the nested problem: the held
#' quantities (the coefficients and the hyperparameters of `which`) and the
#' ones that the criterion estimates are criterion coordinates together, the
#' exact Hessian of the criterion over them gives the derivative of the
#' estimated quantities with respect to the held ones, and the other
#' coefficients are the penalized mode at the predicted values. The first step
#' starts from the refit with the held quantities at the estimate, which
#' differs from the fit where the criterion does not integrate the held
#' coefficients.
#'
#' A hyperparameter in `which` is held on the grid of its link scale, and the
#' profile is then the column `outer_criterion`, evaluated exactly at the
#' predicted point. Where nothing else is estimated by the criterion, the
#' predicted point is the penalized mode at the held values.
#'
#' The method requires a model whose coefficients are smooth functions of the
#' held quantities: the terms [modelterms7::linpar()], [modelterms7::s()],
#' [modelterms7::te()], [modelterms7::random()], [modelterms7::ridge()],
#' [modelterms7::nl()] and the smoothed forms of [modelterms7::seg()],
#' [modelterms7::jump()] and [modelterms7::jseg()]. A model with a
#' break-point term that is not smoothed, a kinked penalty or a structural
#' term is rejected with a pointer to [profile_grid()].
#'
#' @param order The order of the predictor: 1 for the tangent, 2 for the
#'   tangent and the curvature.
#' @param maxit The number of passes the corrector may take at each point, or
#'   `Inf` for convergence.
#' @param name The short name of the method, a single string. The class
#'   constructor `ProfileContinuation()` takes it, and `profile_continuation()`
#'   sets it to `"continuation"`.
#'
#' @return A `ProfileContinuation` object, inheriting from [profile_method],
#'   with the properties `name`, `order` and `maxit`.
#'
#' @seealso [statmod_profile()], [profile_grid()], [profile_taylor()].
#'
#' @examples
#' profile_continuation()
#' profile_continuation(order = 2, maxit = 1)
#' @export
profile_continuation <- function(order = 1, maxit = Inf) {
  if (!is.numeric(order) || length(order) != 1L || is.na(order)) {
    stop("'order' must be 1 or 2.", call. = FALSE)
  }
  ProfileContinuation(name = "continuation", order = as.integer(order),
                      maxit = as.numeric(maxit))
}


#' @title Profile by Local Expansion
#' @name profile_taylor
#'
#' @description
#' Selects a method of [statmod_profile()] that does not refit. The
#' coefficients that are not held are developed as a polynomial in the held
#' ones around the estimate, and the log-likelihood and the penalty are
#' computed exactly at the coefficients that the polynomial gives.
#'
#' @details
#' With \eqn{\theta} the held coefficients and \eqn{\eta} the others, the
#' refitted \eqn{\hat\eta(\theta)} solves \eqn{\nabla_\eta F(\eta, \theta) = 0}
#' for the objective \eqn{F} of the fit. Along the segment from the estimate
#' to a point of the grid, the derivatives of \eqn{\hat\eta} are obtained one
#' order after another: the first is
#' \eqn{-H_{\eta\eta}^{-1} H_{\eta\theta}\,\delta}, and each following one is
#' minus \eqn{H_{\eta\eta}^{-1}} times the derivative of that order, along the
#' segment, of the gradient evaluated at the polynomial of the lower orders.
#' That derivative is taken by one central stencil of numericals7 on the
#' analytic gradient. The polynomial has degree `order - 1`, so that `order`
#' is the highest order of the derivatives of the objective that the
#' expansion uses. The Hessian is the observed information of the fit, with,
#' for a model with nonlinear or smoothed break-point terms, the term of the
#' block that moves with its coefficients.
#' With the log-likelihood and the penalty computed exactly at the predicted
#' coefficients, the penalized log-likelihood of every row is a lower bound
#' of the profile. The column `converged` is `NA`.
#'
#' Quantities that the outer criterion estimates (coordinates, and
#' hyperparameters with `hyper = "refit"`) are developed as in
#' [profile_continuation()], from the derivatives of the criterion in those
#' quantities, and the coefficients that are neither held nor estimated by the
#' criterion are then the penalized mode at the predicted values (one inner fit
#' for each point). The lower bound does not hold in that case, because the
#' exact method maximizes the criterion and not the penalized log-likelihood
#' over those quantities. A hyperparameter in `which` is held at the values of
#' its grid on the scale of its link, and the column `outer_criterion` is the
#' profile.
#'
#' The requirements on the model are those of [profile_continuation()].
#'
#' @param order The highest order of the derivatives of the objective used, an
#'   integer from 2 to 5.
#' @param name The short name of the method, a single string. The class
#'   constructor `ProfileTaylor()` takes it, and `profile_taylor()` sets it to
#'   `"taylor"`.
#'
#' @return A `ProfileTaylor` object, inheriting from [profile_method], with the
#'   properties `name` and `order`.
#'
#' @seealso [statmod_profile()], [profile_grid()], [profile_continuation()].
#'
#' @examples
#' profile_taylor()
#' profile_taylor(order = 4)
#' @export
profile_taylor <- function(order = 2) {
  if (!is.numeric(order) || length(order) != 1L || is.na(order)) {
    stop("'order' must be an integer from 2 to 5.", call. = FALSE)
  }
  ProfileTaylor(name = "taylor", order = as.integer(order))
}


#' Check That a Profile Can Be Approximated
#'
#' @description
#' Rejects, with a pointer to [profile_grid()], a call of an approximate
#' method with `hyper = "refit"` or on a model whose refitted coefficients are
#' not smooth functions of the held ones.
#'
#' @param ctx The context of the call.
#' @param method The name of the method, for the message.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
profile_approx_check <- function(ctx, method) {
  hint <- "Use profile_grid()."
  spec <- ctx$spec
  design <- statmod_design(spec)
  pth <- ctx$hyp$key[ctx$hyp$path]
  if (length(pth)) {
    stop(sprintf(paste0("%s cannot follow %s. The penalty of that ",
                        "hyperparameter has a kink, so the penalized mode ",
                        "is a piecewise smooth function of it and the ",
                        "expansion has no derivative at the corners. %s"),
                 method, paste(sprintf("'%s'", pth), collapse = ", "), hint),
         call. = FALSE)
  }
  if (length(kinked_coords(spec, design))) {
    stop(sprintf(paste0("%s needs coefficients that are smooth functions of ",
                        "the held ones, and the model has a kinked penalty. ",
                        "%s"), method, hint), call. = FALSE)
  }
  st <- attr(design, "state")
  for (p in ctx$params) {
    for (nm in names(st$terms[[p]])) {
      tm <- st$terms[[p]][[nm]]
      if (S7::S7_inherits(tm, modelterms7::SegTerm) &&
          is.null(tm@build_info$smooth)) {
        stop(sprintf(paste0("%s needs coefficients that are smooth functions ",
                            "of the held ones, and the model has the ",
                            "break-point term '%s' without smoothing. %s"),
                     method, nm, hint), call. = FALSE)
      }
    }
  }
  invisible(NULL)
}


#' The Objective, Its Hessian and the Scale of the Held Coefficients at a Point
#'
#' @description
#' Prepares what the polynomial of the trace needs at the coefficients `b0`:
#' the objective of the fit with its gradient, the factorized Hessian of the
#' coefficients that are not held, and the standard error of each held
#' coefficient in the profile, which scales the stencil step.
#'
#' @param ctx The context of the call.
#' @param b0 The stacked coefficients of the point.
#'
#' @return A list.
#'
#' @keywords internal
profile_setup <- function(ctx, b0) {
  spec <- ctx$spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, ctx$hyper0, design, FALSE, NULL)
  lab <- coef_labels(spec, design)
  p <- length(b0)
  alias <- which(rownames(lab) %in% ctx$fit@aliased)
  free <- setdiff(seq_len(p), c(ctx$at, alias))
  H <- as.matrix(obj$he(b0))
  # a block that moves with its coefficients adds the score times the
  # derivative of the block to the scoring pieces
  rf <- attr(design, "refresh")
  if (length(rf)) {
    if (any(vapply(rf, function(r) isTRUE(r$frozen), logical(1)))) {
      stop("A working linearization without a derivative (a sharp jump() or ",
           "jseg()) has no Hessian here. Use profile_grid().", call. = FALSE)
    }
    npar <- vapply(ctx$params, function(q) design[[q]]$npar, integer(1))
    H <- H + mode_curvature(spec, design, obj$split(b0), ctx$params, npar,
                            cumsum(npar) - npar, sum(npar))
  }
  H <- (H + t(H)) / 2
  fac <- tryCatch(chol(H[free, free, drop = FALSE]), error = function(e) NULL)
  if (is.null(fac)) {
    stop("The Hessian of the coefficients that are not held is not positive ",
         "definite at this point. Use profile_grid().", call. = FALSE)
  }
  Hft <- H[free, ctx$at, drop = FALSE]
  schur <- H[ctx$at, ctx$at, drop = FALSE] - crossprod(Hft, chol2inv(fac) %*% Hft)
  se <- sqrt(diag(solve(schur)))
  list(obj = obj, design = design, b0 = b0, at = ctx$at, free = free, H = H,
       fac = fac, Hft = Hft, se = se)
}


#' The Coefficients That a Polynomial of the Trace Predicts at a Point
#'
#' @description
#' Develops the coefficients that are not held as a polynomial in the step
#' `delta` of the held ones, to degree `order - 1`, and returns its value at
#' the end of the step.
#'
#' @details
#' The coefficient of degree 1 is \eqn{-H_{ff}^{-1} H_{f\theta}\delta}. For
#' each higher degree \eqn{j}, the gradient of the free coefficients is
#' evaluated along the step at the polynomial of degree \eqn{j - 1}, its
#' \eqn{j}-th derivative at the origin is taken by one central stencil of
#' numericals7 with a step that moves each held coefficient by a tenth of its
#' standard error, and the coefficient of degree \eqn{j} is minus
#' \eqn{H_{ff}^{-1}} times that derivative.
#'
#' @param S The list that [profile_setup()] returns.
#' @param delta The step of the held coefficients, one number for each.
#' @param order One more than the degree of the polynomial.
#'
#' @return The stacked coefficients at the end of the step.
#'
#' @keywords internal
profile_ray <- function(S, delta, order) {
  p <- length(S$b0)
  nrm <- max(abs(delta) / S$se)
  if (!is.finite(nrm) || nrm <= 0) return(S$b0)
  a1 <- numeric(p)
  a1[S$at] <- delta
  a1[S$free] <- -backsolve(S$fac, forwardsolve(t(S$fac), S$Hft %*% delta))
  coefs <- list(a1)
  poly <- function(s) {
    b <- S$b0
    fact <- 1
    for (j in seq_along(coefs)) {
      fact <- fact * j
      b <- b + coefs[[j]] * s^j / fact
    }
    b
  }
  h <- 0.1 / nrm
  for (j in seq_len(max(order - 1L, 1L))[-1L]) {
    m <- ceiling(j / 2) + 2L
    off <- -m:m
    w <- numericals7::fd_weights(off, j)
    d <- 0
    for (k in seq_along(off)) {
      d <- d + w[k] * S$obj$gr(poly(off[k] * h))[S$free]
    }
    d <- d / h^j
    aj <- numeric(p)
    aj[S$free] <- -backsolve(S$fac, forwardsolve(t(S$fac), d))
    coefs[[j]] <- aj
  }
  poly(1)
}


#' A Subtable of the Hyperparameters, With Its Attributes
#'
#' @param htab The table of [profile_hyper_table()], with the column `row`.
#' @param rows Row numbers of `htab`.
#'
#' @return The rows, with the attributes `links` and `members` of `htab`.
#'
#' @keywords internal
profile_hyper_rows <- function(htab, rows) {
  out <- htab[rows, , drop = FALSE]
  attr(out, "links") <- attr(htab, "links")
  attr(out, "members") <- attr(htab, "members")
  out
}


#' The Nested Problem at a Point
#'
#' @description
#' Prepares what the expansion of the quantities that the outer criterion
#' estimates needs at a point of the profile. The held quantities are the
#' coefficients and the hyperparameters of `which`; the unknown ones are the
#' coefficients that the criterion estimates and, with `hyper = "refit"`, the
#' other hyperparameters. Both sets are coordinates of the criterion, and the
#' function holds the exact Hessian of the criterion over them, with the
#' functions that solve the inner problem and read the exact gradient of the
#' criterion at given values.
#'
#' @details
#' With the held quantities \eqn{h} and the unknown ones \eqn{u}, the
#' criterion \eqn{V(u, h)} is evaluated at the penalized mode of the other
#' coefficients, and none of the coefficient coordinates enters its
#' determinant. The result of [statmod_marginal_hess()] over all of them holds
#' the blocks \eqn{V_{uu}} and \eqn{V_{uh}}, and the derivative of the unknown
#' quantities is \eqn{du/dh = -V_{uu}^{-1} V_{uh}}. The hyperparameters are on
#' the scale of their links, the scale of the grid.
#'
#' @param ctx The context of the call.
#' @param pt A point, as [profile_point()] returns it, with the held values
#'   in `val`.
#'
#' @return A list, or `NULL` where the Hessian of the criterion is not
#'   available.
#'
#' @keywords internal
profile_nested_state <- function(ctx, pt) {
  htab <- ctx$htab
  links <- attr(htab, "links")
  spec_t <- profile_hold(ctx$spec, ctx$param, ctx$coefname,
                         pt$val[ctx$cidx])
  if (identical(ctx$hyper, "held")) {
    oth <- setdiff(which(!htab$path), ctx$hyp$row)
    if (length(oth)) {
      eta0 <- vapply(oth, function(r) {
        linkfunctions7::linkfun(links[[r]],
          pt$hyper[[htab$parameter[r]]][[htab$term[r]]][[htab$name[r]]])
      }, numeric(1))
      spec_t <- profile_hold_hyper(spec_t, profile_hyper_rows(htab, oth), eta0)
    }
  }
  design <- statmod_design(spec_t)
  blocks <- statmod_blocks(spec_t, design)
  meth <- ctx$outer
  gam0 <- if (outer_minimize(meth)) marginal_coords(spec_t, design, NULL) else
    marginal_coords(spec_t, design, meth)
  gam1 <- list(where = c(gam0$where, ctx$at),
               param = c(gam0$param, ctx$param),
               name = c(gam0$name, ctx$coefname))
  idx <- outer_hyper_index(spec_t, blocks)
  nh <- nrow(idx)
  ng <- length(gam0$where)
  k <- length(ctx$at)
  # which row of the criterion each held hyperparameter is
  hrow <- match(paste(idx$parameter, idx$term, idx$name, sep = "/"),
                htab$key)
  held_row <- match(ctx$hyp$row, hrow)
  unk_h <- setdiff(seq_len(nh), held_row)
  upos <- c(unk_h, nh + seq_len(ng))
  # the position in the criterion of each label of `which`
  kpos <- integer(length(ctx$which))
  kpos[ctx$hidx] <- held_row
  kpos[ctx$cidx] <- nh + ng + seq_len(k)
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  split <- function(b) {
    stats::setNames(lapply(seq_along(npar), function(a)
      b[offs[a] + seq_len(npar[a])]), ctx$params)
  }
  hyper_of <- function(uh, val) {
    eta <- numeric(nh)
    eta[unk_h] <- uh
    eta[held_row] <- val[ctx$hidx]
    eta_to_hyper(eta, idx, pt$hyper)
  }
  at_values <- function(g, val, start) {
    sp <- profile_hold(spec_t, ctx$param, ctx$coefname, val[ctx$cidx])
    sp <- statmod_hold(sp, gam0, g)
    b0 <- start
    b0[gam0$where] <- g
    b0[ctx$at] <- val[ctx$cidx]
    list(spec = sp, start = b0)
  }
  # the penalized mode of the other coefficients with the held ones fixed
  mode <- function(g, val, hy, start) {
    h <- at_values(g, val, start)
    d <- statmod_design(h$spec)
    res <- statmod_alternate(h$spec, d, statmod_blocks(h$spec, d), hy,
                             ctx$inner, h$start, ctx$cfg$expected,
                             ctx$cfg$emethod, ctx$cfg$maxit, ctx$cfg$tol,
                             verbosity(0))
    list(par = res$par, spec = h$spec, design = d)
  }
  # the exact gradient and Hessian of the criterion at a mode
  derivs <- function(m, hy, hessian) {
    basis <- integrated_basis(m$spec, m$design, meth@kind, gamma = ng > 0L)
    cf <- split(m$par)
    oc <- outer_context(m$spec, m$design, cf, hy, NULL)
    if (hessian) {
      statmod_marginal_hess(m$spec, m$design, cf, hy, meth, idx, basis,
                            ctx = oc, gam = gam1)
    } else {
      statmod_marginal_grad(m$spec, m$design, cf, hy, meth, idx, basis,
                            ctx = oc, gam = gam1)
    }
  }
  g0 <- pt$par[gam0$where]
  eta_h <- hyper_to_eta(pt$hyper, idx)
  m0 <- list(par = pt$par, spec = at_values(g0, pt$val, pt$par)$spec,
             design = design)
  V <- derivs(m0, pt$hyper, TRUE)
  if (is.null(V) || !all(is.finite(V))) return(NULL)
  list(ctx = ctx, nu = length(upos), upos = upos, kpos = kpos, unk_h = unk_h,
       gam = gam0$where, par = pt$par, val = pt$val, hyper = pt$hyper,
       u0 = c(eta_h[unk_h], g0), Vuu = V[upos, upos, drop = FALSE],
       Vuk = V[upos, kpos, drop = FALSE], mode = mode, derivs = derivs,
       hyper_of = hyper_of, nh = nh, ng = ng, meth = meth)
}


#' The Coefficients That the Nested Expansion Predicts at a Step
#'
#' @description
#' Develops the unknown quantities of the nested problem as a polynomial in
#' the step of the held ones, and returns the point at its end: the
#' coefficients at their penalized mode, the hyperparameters, and the
#' specification that holds the coefficients that the criterion does not
#' integrate.
#'
#' @details
#' The coefficient of degree 1 is \eqn{-V_{uu}^{-1} V_{uh}\delta}. For a
#' higher degree \eqn{j}, the exact gradient of the criterion in the unknown
#' quantities is evaluated along the step, at the penalized mode of the other
#' coefficients, with the unknown quantities at the polynomial of lower
#' degree. Its \eqn{j}-th derivative at the origin is taken by one central
#' stencil of numericals7, and the coefficient of degree \eqn{j} is minus
#' \eqn{V_{uu}^{-1}} times that derivative. Where nothing is unknown the
#' point is the penalized mode at the held values.
#'
#' @param N The list that [profile_nested_state()] returns.
#' @param delta The step of the held values, one number for each label of
#'   `which`.
#' @param order One more than the degree of the polynomial.
#' @param se The scale of each label of `which`.
#'
#' @return A list with `par`, `hyper`, `spec` and `design`.
#'
#' @keywords internal
profile_nested_predict <- function(N, delta, order, se) {
  val0 <- N$val
  u0 <- N$u0
  nu <- N$nu
  nrm <- max(abs(delta) / se)
  ng <- N$ng
  gam_of <- function(u) u[length(N$unk_h) + seq_len(ng)]
  hyp_of <- function(u) u[seq_along(N$unk_h)]
  coefs <- list()
  if (nu && is.finite(nrm) && nrm > 0) {
    coefs[[1L]] <- -solve(N$Vuu, N$Vuk %*% delta)
  }
  poly <- function(s) {
    u <- u0
    fact <- 1
    for (j in seq_along(coefs)) {
      fact <- fact * j
      u <- u + as.numeric(coefs[[j]]) * s^j / fact
    }
    u
  }
  start <- N$par
  if (nu && length(coefs)) {
    h <- 0.1 / nrm
    for (j in seq_len(max(order - 1L, 1L))[-1L]) {
      m <- ceiling(j / 2) + 2L
      off <- -m:m
      w <- numericals7::fd_weights(off, j)
      d <- 0
      for (i in seq_along(off)) {
        s <- off[i] * h
        u <- poly(s)
        val <- val0 + s * delta
        hy <- N$hyper_of(hyp_of(u), val)
        mo <- N$mode(gam_of(u), val, hy, start)
        d <- d + w[i] * N$derivs(mo, hy, FALSE)[N$upos]
      }
      coefs[[j]] <- -solve(N$Vuu, d / h^j)
    }
  }
  u1 <- poly(1)
  val1 <- val0 + delta
  hy1 <- N$hyper_of(hyp_of(u1), val1)
  mo <- N$mode(gam_of(u1), val1, hy1, start)
  c(mo, list(hyper = hy1))
}


#' The Exact Quantities of a Point Predicted by the Nested Expansion
#'
#' @description
#' Computes the penalized log-likelihood, its two parts and, where the profile
#' reads it, the outer criterion, exactly at a point that
#' [profile_nested_predict()] returns.
#'
#' @param ctx The context of the call.
#' @param pr The predicted point.
#' @param val The held values.
#'
#' @return A point, as [profile_point()] returns it, with `converged` `NA`.
#'
#' @keywords internal
profile_exact_nested <- function(ctx, pr, val) {
  design <- pr$design
  cf <- stats::setNames(
    lapply(seq_along(ctx$params), function(a) {
      npar <- vapply(design, function(d) d$npar, integer(1))
      pr$par[(cumsum(npar) - npar)[a] + seq_len(npar[a])]
    }), ctx$params)
  ll <- statmod_loglik_at(pr$spec, cf, design)
  pen <- statmod_penalty_at(pr$spec, cf, pr$hyper, design, "value")
  obj <- -ll + pen
  if (!is.finite(obj)) {
    return(profile_failed(ctx, "the predicted point is infeasible"))
  }
  m <- statmod_marginal(pr$spec, design, cf, pr$hyper, ctx$outer,
                        basis = integrated_basis(pr$spec, design,
                                                 ctx$outer@kind,
                                                 gamma = ctx$ng > 0L))
  list(par = pr$par, hyper = pr$hyper, objective = obj, loglik = ll,
       criterion = if (is.null(m)) NA_real_ else m$value, converged = NA)
}


#' The Scale of Each Label of a Profile
#'
#' @description
#' The standard error of each label of `which`: for a coefficient the one of
#' the profile, from [profile_setup()], and for a hyperparameter the one of
#' the exact Hessian of the criterion on the scale of its link, or 1 where it
#' is not available. It scales the step of the stencil.
#'
#' @param ctx The context of the call.
#'
#' @return A numeric vector.
#'
#' @keywords internal
profile_scales <- function(ctx) {
  se <- rep(1, length(ctx$which))
  if (length(ctx$cidx)) {
    se[ctx$cidx] <- profile_setup(ctx, ctx$beta_fit)$se
  }
  if (length(ctx$hidx)) {
    sh <- profile_hyper_se(ctx)
    se[ctx$hidx] <- ifelse(is.finite(sh), sh, 1)
  }
  se
}


#' The Exact Penalized Log-Likelihood at Predicted Coefficients
#'
#' @param ctx The context of the call.
#' @param S The list that [profile_setup()] returns.
#' @param b The predicted stacked coefficients.
#' @param val The held values.
#'
#' @return A point, as [profile_point()] returns it, with `converged` `NA`.
#'
#' @keywords internal
profile_exact <- function(ctx, S, b, val) {
  b[ctx$at] <- val
  obj <- S$obj$fn(b)
  if (!is.finite(obj)) return(profile_failed(ctx, "the predicted point is infeasible"))
  list(par = b, hyper = ctx$hyper0, objective = obj,
       loglik = statmod_loglik_at(ctx$spec, S$obj$split(b), S$design),
       criterion = NA_real_, converged = NA)
}


#' Whether a Profile Follows the Nested Problem
#'
#' @description
#' The approximate methods follow the nested problem where the outer criterion
#' estimates anything at a point (coordinates, or hyperparameters with
#' `hyper = "refit"`) and where a hyperparameter is held.
#'
#' @param ctx The context of the call.
#'
#' @return A single logical.
#'
#' @keywords internal
profile_nested_needed <- function(ctx) {
  isTRUE(ctx$use_outer) || length(ctx$sidx) > 0L
}


#' @rdname profile_compute
#' @name profile_compute.ProfileTaylor
#' @keywords internal
S7::method(profile_compute, ProfileTaylor) <- function(method, ctx, grid, cores) {
  profile_approx_check(ctx, "profile_taylor()")
  nested <- profile_nested_needed(ctx)
  vals <- as.matrix(grid)
  if (!nested) {
    S <- profile_setup(ctx, ctx$beta_fit)
    return(lapply(seq_len(nrow(vals)), function(r) {
      tryCatch(profile_exact(ctx, S,
                             profile_ray(S, vals[r, ] - ctx$estimate,
                                         method@order), vals[r, ]),
               error = function(e) profile_failed(ctx, conditionMessage(e)))
    }))
  }
  # the expansion of a nested problem starts at the refit with the held
  # quantities at the estimate: with the criterion estimating the others,
  # that refit is not the fit
  base <- profile_point(ctx, ctx$estimate)
  base$val <- ctx$estimate
  N <- profile_nested_state(ctx, base)
  if (is.null(N)) {
    stop("The Hessian of the criterion is not available at the estimate. ",
         "Use profile_grid().", call. = FALSE)
  }
  se <- profile_scales(ctx)
  lapply(seq_len(nrow(vals)), function(r) {
    tryCatch({
      pr <- profile_nested_predict(N, vals[r, ] - ctx$estimate, method@order,
                                   se)
      profile_exact_nested(ctx, pr, vals[r, ])
    }, error = function(e) profile_failed(ctx, conditionMessage(e)))
  })
}


#' One Step of a Continuation
#'
#' @description
#' Predicts the coefficients at `val` from the point `prev` and corrects them.
#' A failed corrector halves the step through a point that is not returned.
#'
#' @param ctx The context of the call.
#' @param method The [profile_continuation()] object.
#' @param prev The previous point, with the held values in `val`.
#' @param val The held values to reach.
#' @param nested Whether the profile follows the nested problem.
#' @param se The scale of each label of `which`, for the nested problem.
#' @param depth The number of halvings so far.
#'
#' @return The point at `val`, with `val` in it.
#'
#' @keywords internal
profile_cont_step <- function(ctx, method, prev, val, nested, se, depth = 0L) {
  delta <- val - prev$val
  pred <- tryCatch({
    if (nested) {
      N <- profile_nested_state(ctx, prev)
      if (is.null(N)) {
        list(par = prev$par, hyper = prev$hyper)
      } else {
        profile_nested_predict(N, delta, method@order + 1L, se)
      }
    } else {
      list(par = profile_ray(profile_setup(ctx, prev$par), delta,
                             method@order + 1L), hyper = prev$hyper)
    }
  }, error = function(e) list(par = prev$par, hyper = prev$hyper))
  start <- list(par = pred$par, hyper = pred$hyper)
  pt <- tryCatch(
    profile_point(ctx, val, start, maxit = method@maxit),
    error = function(e) profile_failed(ctx, conditionMessage(e)))
  bad <- !is.null(pt$error) || !is.finite(pt$objective) ||
    (is.infinite(method@maxit) && !isTRUE(pt$converged))
  if (bad && depth < 8L) {
    mid <- prev$val + delta / 2
    pm <- profile_cont_step(ctx, method, prev, mid, nested, se, depth + 1L)
    if (!is.null(pm$error)) return(pm)
    return(profile_cont_step(ctx, method, pm, val, nested, se, depth + 1L))
  }
  pt$val <- val
  pt
}


#' @rdname profile_compute
#' @name profile_compute.ProfileContinuation
#' @keywords internal
S7::method(profile_compute, ProfileContinuation) <- function(method, ctx, grid, cores) {
  profile_approx_check(ctx, "profile_continuation()")
  nested <- profile_nested_needed(ctx)
  base <- list(par = ctx$beta_fit, hyper = ctx$hyper0, val = ctx$estimate)
  se <- NULL
  if (nested) {
    base <- profile_point(ctx, ctx$estimate)
    base$val <- ctx$estimate
    se <- profile_scales(ctx)
  }
  chains <- profile_chains(ctx, grid, "neighbor")
  if (cores > 1L) ctx$spec <- S7::set_props(ctx$spec, threads = 1L)
  vals <- as.matrix(grid)
  body <- function(i) {
    ch <- chains[[i]]
    prev <- base
    out <- vector("list", length(ch$rows))
    for (m in seq_along(ch$rows)) {
      pt <- profile_cont_step(ctx, method, prev, vals[ch$rows[m], ], nested, se)
      out[[m]] <- pt
      if (!is.null(pt$error)) break
      prev <- pt
    }
    out
  }
  res <- profile_map(ctx, cores, length(chains), body)
  pts <- vector("list", nrow(grid))
  for (i in seq_along(chains)) {
    for (m in seq_along(chains[[i]]$rows)) {
      pt <- res[[i]][[m]]
      pts[[chains[[i]]$rows[m]]] <- if (is.null(pt)) {
        profile_failed(ctx, "the continuation stopped at an earlier point")
      } else pt
    }
  }
  pts
}
