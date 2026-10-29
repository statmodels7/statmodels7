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
#' method. Their tangent is computed from the nested problem: the refit of
#' [profile_grid()] at two values on either side of the estimate gives, by a
#' central difference, the derivative of all coefficients with respect to each
#' held one, and this tangent is applied from every previous point. The first
#' step starts from the refit with the held coefficients at the estimate,
#' which differs from the fit where the criterion does not integrate the held
#' coefficients. In that case only `order = 1` is available.
#'
#' The method requires `hyper = "held"` and a model whose coefficients are
#' smooth functions of the held ones: the terms [modelterms7::linpar()],
#' [modelterms7::s()], [modelterms7::te()], [modelterms7::random()],
#' [modelterms7::ridge()], [modelterms7::nl()] and the smoothed forms of
#' [modelterms7::seg()], [modelterms7::jump()] and [modelterms7::jseg()].
#' A model with a break-point term that is not smoothed, a kinked penalty or
#' a structural term is rejected with a pointer to [profile_grid()].
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
#' expansion uses. The Hessian is the analytic one, or, for a model with
#' nonlinear or smoothed break-point terms, the difference of the gradient.
#' With the log-likelihood and the penalty computed exactly at the predicted
#' coefficients, the penalized log-likelihood of every row is a lower bound
#' of the profile. The column `converged` is `NA`.
#'
#' Coordinates that the outer criterion estimates are handled as in
#' [profile_continuation()], by a tangent from the nested problem, and the
#' order is then limited to 2. The lower bound does not hold in that case,
#' because the exact method maximizes the criterion and not the penalized
#' log-likelihood over those coordinates.
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
  if (!identical(ctx$hyper, "held")) {
    stop(sprintf("%s requires hyper = \"held\". %s", method, hint),
         call. = FALSE)
  }
  spec <- ctx$spec
  design <- statmod_design(spec)
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


#' The Hessian of an Objective by Differences of Its Gradient
#'
#' @param gr The gradient, a function of the stacked coefficients.
#' @param b The coefficients.
#'
#' @return A symmetric matrix.
#'
#' @keywords internal
profile_fd_hessian <- function(gr, b) {
  p <- length(b)
  H <- matrix(0, p, p)
  for (j in seq_len(p)) {
    h <- 1e-5 * max(1, abs(b[j]))
    e <- numeric(p)
    e[j] <- h
    H[, j] <- (gr(b + e) - gr(b - e)) / (2 * h)
  }
  (H + t(H)) / 2
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
  refresh <- length(attr(design, "refresh")) > 0L
  H <- if (refresh) profile_fd_hessian(obj$gr, b0) else
    (function(h) (h + t(h)) / 2)(as.matrix(obj$he(b0)))
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


#' The Jacobian of the Coefficients With Respect to the Held Ones, From the
#' Nested Problem
#'
#' @description
#' Differentiates by a central difference the refits of [profile_point()] at
#' two values of each held coefficient, a tenth of a standard error on either
#' side of the estimate. The refits re-estimate the coordinates that the outer
#' criterion estimates, so the result is the tangent of the nested problem.
#'
#' @param ctx The context of the call.
#' @param S The list that [profile_setup()] returns at the estimate.
#'
#' @return A matrix with one row per coefficient and one column per held one.
#'
#' @keywords internal
profile_nested_jacobian <- function(ctx, S) {
  k <- length(ctx$at)
  J <- matrix(0, length(S$b0), k)
  for (m in seq_len(k)) {
    eps <- 0.1 * S$se[m]
    up <- dn <- ctx$estimate
    up[m] <- up[m] + eps
    dn[m] <- dn[m] - eps
    zp <- profile_point(ctx, up)
    zm <- profile_point(ctx, dn)
    J[, m] <- (zp$par - zm$par) / (2 * eps)
  }
  J
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


#' @rdname profile_compute
#' @name profile_compute.ProfileTaylor
#' @keywords internal
S7::method(profile_compute, ProfileTaylor) <- function(method, ctx, grid, cores) {
  profile_approx_check(ctx, "profile_taylor()")
  nested <- ctx$ng > 0L
  if (nested && method@order > 2L) {
    stop("profile_taylor() with order above 2 does not handle coordinates ",
         "that the outer criterion estimates. Use order = 2 or profile_grid().",
         call. = FALSE)
  }
  S <- profile_setup(ctx, ctx$beta_fit)
  # the expansion of a nested problem starts at the refit with the held
  # coefficients at the estimate: with the criterion estimating the other
  # coordinates, that refit is not the fit
  base <- if (nested) profile_point(ctx, ctx$estimate)
  J <- if (nested) profile_nested_jacobian(ctx, S)
  vals <- as.matrix(grid)
  lapply(seq_len(nrow(vals)), function(r) {
    tryCatch({
      delta <- vals[r, ] - ctx$estimate
      b <- if (nested) base$par + as.numeric(J %*% delta) else
        profile_ray(S, delta, method@order)
      profile_exact(ctx, S, b, vals[r, ])
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
#' @param J The nested Jacobian, or `NULL`.
#' @param depth The number of halvings so far.
#'
#' @return The point at `val`, with `val` in it.
#'
#' @keywords internal
profile_cont_step <- function(ctx, method, prev, val, J, depth = 0L) {
  delta <- val - prev$val
  pred <- if (is.null(J)) {
    profile_ray(profile_setup(ctx, prev$par), delta, method@order + 1L)
  } else {
    prev$par + as.numeric(J %*% delta)
  }
  pred[ctx$at] <- val
  pt <- tryCatch(
    profile_point(ctx, val, list(par = pred, hyper = prev$hyper),
                  maxit = method@maxit),
    error = function(e) profile_failed(ctx, conditionMessage(e)))
  bad <- !is.null(pt$error) || !is.finite(pt$objective) ||
    (is.infinite(method@maxit) && !isTRUE(pt$converged))
  if (bad && depth < 8L) {
    mid <- prev$val + delta / 2
    pm <- profile_cont_step(ctx, method, prev, mid, J, depth + 1L)
    if (!is.null(pm$error)) return(pm)
    return(profile_cont_step(ctx, method, pm, val, J, depth + 1L))
  }
  pt$val <- val
  pt
}


#' @rdname profile_compute
#' @name profile_compute.ProfileContinuation
#' @keywords internal
S7::method(profile_compute, ProfileContinuation) <- function(method, ctx, grid, cores) {
  profile_approx_check(ctx, "profile_continuation()")
  nested <- ctx$ng > 0L
  if (nested && method@order > 1L) {
    stop("profile_continuation() with order = 2 does not handle coordinates ",
         "that the outer criterion estimates. Use order = 1 or profile_grid().",
         call. = FALSE)
  }
  base <- list(par = ctx$beta_fit, hyper = ctx$hyper0, val = ctx$estimate)
  J <- NULL
  if (nested) {
    base <- profile_point(ctx, ctx$estimate)
    base$val <- ctx$estimate
    J <- profile_nested_jacobian(ctx, profile_setup(ctx, ctx$beta_fit))
  }
  chains <- profile_chains(ctx, grid, "neighbor")
  if (cores > 1L) ctx$spec <- S7::set_props(ctx$spec, threads = 1L)
  vals <- as.matrix(grid)
  body <- function(i) {
    ch <- chains[[i]]
    prev <- base
    out <- vector("list", length(ch$rows))
    for (m in seq_along(ch$rows)) {
      pt <- profile_cont_step(ctx, method, prev, vals[ch$rows[m], ], J)
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
