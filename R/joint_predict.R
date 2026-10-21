#' @include joint_scaling.R
NULL

#' The Derivative Rows of Every Equation Under a Joint Filter
#'
#' @description
#' Returns the derivative of each equation's predictor in the coefficients of
#' every equation and the joint filter's free parameters, at the fitted
#' point, together with the first derivative of the driving vector in every
#' predictor.
#'
#' @details
#' The rows of a filtered equation are the forward Jacobian of the joint
#' recursion from [joint_filter_curvature()]; those of the other equations
#' are their static designs. The columns of a held level are left out, as
#' the variance matrix leaves them out.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#'
#' @return `NULL` where the model carries no joint filter, otherwise a list
#'   with `V` (one matrix per distribution parameter), `U1` (the derivative
#'   of each driving score in every predictor), `G` (the positions of the
#'   filtered parameters), `nb`, `key` (the names of the columns in the
#'   variance matrix), `free` (the filter's free parameters), `f` and `ev`.
#'
#' @seealso [structural_se_columns()], [joint_continued_deriv_inputs()]
#'
#' @keywords internal
joint_filter_jacobian <- function(spec, design, coef) {
  jd <- joint_design_rows(spec, design, coef)
  if (is.null(jd) || !isTRUE(jd$f$joint)) return(NULL)
  d <- spec@distrib
  th <- jd$ev$theta
  gl <- distributions7::distrib_gradient(d, spec@response, th, scale = "link",
                                         threads = spec@threads)
  H <- distributions7::distrib_hessian(d, spec@response, th, scale = "link",
                                       threads = spec@threads)
  D3 <- distributions7::distrib_deriv3(d, spec@response, th, scale = "link",
                                       threads = spec@threads)
  f <- jd$f
  params <- jd$params
  G <- match(f$params, params)
  cv <- joint_filter_curvature(spec, design, f, jd$V, gl, H, D3, th)
  V <- jd$V
  for (k in seq_along(G)) V[[G[k]]] <- cv$jacobian[[k]]
  V <- lapply(V, function(x) x[, jd$keep, drop = FALSE])
  names(V) <- params
  U1 <- joint_drive_data(spec, f$tm, gl, H, order = 1L, theta = th)$U1
  ckey <- unlist(lapply(params, function(q) {
    if (design[[q]]$npar) paste(q, design[[q]]$coef_names, sep = ":") else
      character(0)
  }), use.names = FALSE)
  ds <- structural_display(f$tm, f$param, jd$zn)
  tkey <- paste(ds$param, ds$name, sep = ":")
  list(V = V, U1 = U1, G = G, nb = jd$nb, key = c(ckey, tkey)[jd$keep],
       free = jd$zn[jd$keep[jd$keep > jd$nb] - jd$nb], f = f, ev = jd$ev,
       keep = jd$keep)
}


#' The Inputs to the Derivative of a Continued Joint Filter
#'
#' @description
#' Returns the derivatives of the levels and of the driving scores at the
#' observed rows, and of the filter's parameters, in the coordinates of the
#' variance matrix, as [modelterms7::term_continue()] takes them for a joint
#' filter.
#'
#' @details
#' The level of a filtered equation is its predictor less the static part,
#' so its derivative row is the equation's row less its design. The driving
#' score of parameter \eqn{l} depends on the predictors of every equation, so
#' its derivative row is \eqn{\sum_b U_{l;b} V_b}.
#'
#' @param ospec,odesign The fit's specification and design.
#' @param coef The coefficients.
#' @param ost The structural state.
#'
#' @return A list with `deriv`, `cols` (for each filtered parameter, the
#'   columns of its equation's coefficients) and `key`.
#'
#' @seealso [statmod_eta_continued()]
#'
#' @keywords internal
joint_continued_deriv_inputs <- function(ospec, odesign, coef, ost) {
  fj <- joint_filter_jacobian(ospec, odesign, coef)
  params <- ospec@distrib@params
  n <- ospec@n_obs
  npar <- vapply(odesign, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  # the coefficient columns survive every held level, which sits in the tail
  m <- length(fj$key)
  cols <- lapply(fj$G, function(a) offs[a] + seq_len(npar[a]))
  df_past <- lapply(seq_along(fj$G), function(k) {
    a <- fj$G[k]
    S <- matrix(0, n, m)
    if (npar[a]) S[, cols[[k]]] <- as_dense(odesign[[params[a]]]$X)
    fj$V[[a]] - S
  })
  ds_past <- lapply(seq_along(fj$G), function(l) {
    out <- matrix(0, n, m)
    for (b in seq_along(params)) out <- out + fj$U1[, l, b] * fj$V[[b]]
    out
  })
  z <- ost$zeta[[fj$f$term]]
  nm <- names(z)
  links <- modelterms7::term_links(fj$f$tm)
  dpsi <- matrix(0, length(nm), m)
  for (i in seq_along(fj$free)) {
    j <- fj$free[i]
    dpsi[match(j, nm), fj$nb + i] <-
      linkfunctions7::dlinkinv(links[[j]], z[[j]])
  }
  list(deriv = list(df_past = df_past, ds_past = ds_past, dpsi = dpsi),
       cols = cols, key = fj$key)
}


#' Simulate a Model Carrying a Joint Filter
#'
#' @description
#' Runs the joint recursion forward drawing the response at each step, the
#' branch of [rstatmod_eta()] for a joint filter.
#'
#' @details
#' At each row the response is drawn at the filtered predictors that the
#' recursion has reached there and the static predictors of the other
#' equations, and the driving vector is read at the drawn value, the same
#' vector, scaled as in the fit, that drove the fitted filter.
#'
#' @param spec A [StatmodSpec()].
#' @param tm The joint filter.
#' @param eta,theta The static predictors and parameters, named lists.
#' @param psi The filter's parameters on the parameter scale.
#'
#' @return The list that [rstatmod_eta()] returns.
#'
#' @keywords internal
rstatmod_eta_joint <- function(spec, tm, eta, theta, psi) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  n <- spec@n_obs
  G <- match(tm@params, params)
  theta_n <- lapply(params, function(q) rep_len(theta[[q]], n))
  names(theta_n) <- params
  draw <- function(e, i) {
    th <- lapply(theta_n, function(v) v[[i]])
    for (k in seq_along(G)) {
      th[[params[G[k]]]] <- linkfunctions7::linkinv(links[[params[G[k]]]],
                                                    e[k])
    }
    as.numeric(distributions7::distrib_rng(spec@distrib, 1L, th))
  }
  cb <- joint_callbacks(spec, tm, eta)
  es <- vapply(tm@params, function(p) rep_len(as.numeric(eta[[p]]), n),
               numeric(n))
  es <- matrix(es, n, length(G))
  res <- modelterms7::term_simulate(tm, psi, es, draw, score = cb$score_y,
                                    curvature = cb$curvature_y)
  for (k in seq_along(G)) {
    pk <- params[G[k]]
    eta[[pk]] <- as.numeric(res$eta[, k])
    theta[[pk]] <- linkfunctions7::linkinv(links[[pk]], eta[[pk]])
  }
  list(eta = eta, theta = theta, y = res$y, latent = res$latent,
       structural = as.list(psi))
}
