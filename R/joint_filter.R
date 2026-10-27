#' @include structural.R filter_driving.R
NULL

# Several gas() terms in the equations of different parameters are one
# recursion: the score of each filtered parameter is read at the levels of
# all the others. statmod_terms() joins them into one modelterms7::GasJoint
# placed where the first of them sat, so that the rest of the layer sees one
# structural term with one vector of own parameters; what differs is that its
# predictor is a matrix, one column per filtered equation, and its driving
# score a vector.

#' Join the Score-Driven Terms of Several Equations
#'
#' @description
#' Replaces the [modelterms7::gas()] terms of a model, when there is more
#' than one, with a single [modelterms7::GasJoint()] placed in the equation
#' of the first of them under the key `"gas_joint"`.
#'
#' @details
#' [reject_incompatible()] has already checked that each equation carries at
#' most one such term and that no other structural term sits beside them.
#' [modelterms7::gas_joint()] checks that the terms walk the same series and
#' share one `scaling`, and fills the cross loadings that the `p_cross` of
#' each term declares. A parameter that a `p_cross` names and that carries no
#' gas() term joins the filter as a driver
#' ([modelterms7::GasDriver()]): its predictor stays static and its score
#' enters the driving vector. A single score-driven term with no driver is
#' left as it is, after the same check of its `p_cross`.
#'
#' @param terms The built terms, a named list of named lists, one per
#'   distribution parameter.
#'
#' @return `terms`, with the score-driven terms joined where there are two or
#'   more.
#'
#' @keywords internal
join_gas_terms <- function(terms) {
  found <- list()
  for (p in names(terms)) {
    for (nm in names(terms[[p]])) {
      if (S7::S7_inherits(terms[[p]][[nm]], modelterms7::GasTerm)) {
        found[[length(found) + 1L]] <- list(param = p, term = nm)
      }
    }
  }
  if (!length(found)) return(terms)
  parts <- stats::setNames(lapply(found, function(f)
    terms[[f$param]][[f$term]]), vapply(found, `[[`, "", "param"))
  # a parameter that a p_cross names with a positive count and that carries
  # no gas() term drives the filter with its score: it must be a parameter
  # of the family, whose equations are the names of `terms`
  named <- unique(unlist(lapply(parts, function(tm)
    names(tm@p_cross)[tm@p_cross > 0L])))
  bad <- setdiff(named, names(terms))
  if (length(bad)) {
    stop(sprintf(paste0(
      "a gas() term loads the score of '%s' (p_cross), which is not a ",
      "parameter of the family;\n  the parameters are %s."), bad[1L],
      paste(names(terms), collapse = ", ")), call. = FALSE)
  }
  drivers <- setdiff(named, names(parts))
  # one term alone, with no driver, is filtered as it is, and gas_joint()
  # only checks its p_cross
  if (length(found) < 2L && !length(drivers)) {
    modelterms7::gas_joint(parts)
    return(terms)
  }
  jt <- modelterms7::gas_joint(parts, drivers = drivers)
  jt@keys <- stats::setNames(vapply(found, `[[`, "", "term"), names(parts))
  for (f in found) terms[[f$param]][[f$term]] <- NULL
  terms[[found[[1L]]$param]][["gas_joint"]] <- jt
  terms
}


#' Is a Structural Term a Joint Filter?
#'
#' @param tm A term.
#'
#' @return A single logical.
#'
#' @keywords internal
is_joint_filter <- function(tm) S7::S7_inherits(tm, modelterms7::GasJoint)


#' The Family at One Observation, in Every Predictor
#'
#' @description
#' Returns a function of the row and of the predictors of every equation at
#' that row, which gives the link-scale gradient and Hessian of the
#' log-density there and, where `info` is `TRUE`, the expected information
#' and its derivative in each predictor.
#'
#' @details
#' A joint filter reads its driving vector at the predictor that it has just
#' produced, one observation at a time, so the family is evaluated at one
#' point per call. The S7 methods of the family and of its links, and the
#' names of the components, are resolved once, and the chain rule onto the
#' link scale is applied here: with \eqn{h_a} the inverse link of parameter
#' \eqn{a}, the gradient is \eqn{\ell_a h_a'} and the Hessian
#' \eqn{\ell_{ab} h_a' h_b' + \delta_{ab}\,\ell_a h_a''}; the expected
#' information \eqn{\mathcal{I}_{ab} = -\mathbb{E}[\ell_{ab}]\,h_a' h_b'}
#' carries no \eqn{h''} term, because the score has mean zero, and its
#' derivative in the predictor of parameter \eqn{c} is
#' \deqn{-\frac{\partial \mathbb{E}[\ell_{ab}]}{\partial\theta_c}\,
#'   h_c' h_a' h_b' - \mathbb{E}[\ell_{ab}]\,(\delta_{ac} h_a'' h_b'
#'   + \delta_{bc} h_a' h_b''),}
#' as [distributions7::distrib_dexpected_hessian()] composes it. Where the
#' family's expected information is not written out
#' ([distributions7::expected_hessian_exact()] is `FALSE`) or its derivative
#' is the default difference, the two generics are called on the link scale
#' instead, so that their own approximation is the one read.
#'
#' @param spec A [StatmodSpec()].
#' @param info Whether the expected information and its derivative are read.
#'
#' @return A function `function(i, e, yi)` of the row, of the vector of all
#'   the predictors at that row and of the response there (by default the
#'   specification's own), returning a list with `g`, `H` and, where
#'   `info` is `TRUE`, `I` and `dI` (an array whose slice `[, , c]` is the
#'   derivative of `I` in the predictor of parameter `c`).
#'
#' @keywords internal
joint_row_family <- function(spec, info = FALSE) {
  d <- spec@distrib
  y <- spec@response
  params <- d@params
  P <- length(params)
  cls <- S7::S7_class(d)
  m_grad <- S7::method(distributions7::distrib_gradient, cls)
  m_hess <- S7::method(distributions7::distrib_hessian, cls)
  # the expected information and its derivative are read on the parameter
  # scale only where both are written out: the default derivative differences
  # on the scale it is asked for, so it is left to the generic
  exact <- info && distributions7::expected_hessian_exact(d)
  if (exact) {
    m_E <- S7::method(distributions7::distrib_expected_hessian, cls)
    m_dE <- S7::method(distributions7::distrib_dexpected_hessian, cls)
    exact <- !identical(m_dE, S7::method(
      distributions7::distrib_dexpected_hessian, distributions7::distrib))
  }
  lks <- d@link_params[params]
  m_inv <- lapply(lks, function(lk)
    S7::method(linkfunctions7::linkinv, S7::S7_class(lk)))
  m_h1 <- lapply(lks, function(lk)
    S7::method(linkfunctions7::dlinkinv, S7::S7_class(lk)))
  m_h2 <- lapply(lks, function(lk)
    S7::method(linkfunctions7::d2linkinv, S7::S7_class(lk)))
  bounds <- lapply(lks, function(lk) lk@link_bounds)
  # the keys of the components, built once: hkeys[a, b] indexes hk, and
  # dkeys[a, b, c] indexes dk
  hkeys <- matrix("", P, P)
  for (a in seq_len(P)) for (b in seq_len(P)) hkeys[a, b] <- hess_key(params, a, b)
  hk <- unique(as.vector(hkeys))
  hidx <- match(hkeys, hk)
  dk <- character(0)
  didx <- integer(0)
  if (info) {
    dkeys <- array("", c(P, P, P))
    for (a in seq_len(P)) for (b in seq_len(P)) for (c in seq_len(P)) {
      dkeys[a, b, c] <- distributions7::dexpected_key(params, a, b, c)
    }
    dk <- unique(as.vector(dkeys))
    didx <- match(dkeys, dk)
  }
  pick <- function(x, keys) unlist(x[keys], use.names = FALSE)
  function(i, e, yi = y[i]) {
    th <- vector("list", P)
    h1 <- numeric(P)
    h2 <- numeric(P)
    for (a in seq_len(P)) {
      th[[a]] <- linkfunctions7::link_bounds_clamp(
        m_inv[[a]](lks[[a]], e[a]), bounds[[a]])
      h1[a] <- m_h1[[a]](lks[[a]], e[a])
      h2[a] <- m_h2[[a]](lks[[a]], e[a])
    }
    names(th) <- params
    gr <- pick(m_grad(d, yi, th, scale = "parameter"), params)
    hh <- h1 %o% h1
    # (l_ab h_a') h_b', in the order of the Hessian's chain rule
    H <- matrix(pick(m_hess(d, yi, th, scale = "parameter"), hk)[hidx], P, P) *
      h1 * rep(h1, each = P)
    diag(H) <- diag(H) + gr * h2
    out <- list(g = gr * h1, H = H)
    if (!info) return(out)
    if (exact) {
      E <- matrix(pick(m_E(d, yi, th, scale = "parameter"), hk)[hidx], P, P)
      dE <- array(pick(m_dE(d, yi, th, scale = "parameter"), dk)[didx],
                  c(P, P, P))
      dI <- array(0, c(P, P, P))
      for (c in seq_len(P)) {
        x <- -dE[, , c] * h1[c] * hh
        x[c, ] <- x[c, ] - E[c, ] * h2[c] * h1
        x[, c] <- x[, c] - E[, c] * h1 * h2[c]
        dI[, , c] <- x
      }
      out$I <- -E * hh
      out$dI <- dI
      return(out)
    }
    E <- distributions7::distrib_expected_hessian(d, yi, th, scale = "link")
    dE <- distributions7::distrib_dexpected_hessian(d, yi, th, scale = "link")
    out$I <- -matrix(pick(E, hk)[hidx], P, P)
    out$dI <- -array(pick(dE, dk)[didx], c(P, P, P))
    out
  }
}


#' The Scaling Matrix of a Joint Filter and Its First Derivative
#'
#' @description
#' Returns \eqn{S = \mathcal{I}^{-d}} for a symmetric positive definite
#' block, with \eqn{d} in 0, 1/2 or 1, and its derivative along each of a
#' set of derivatives of the block.
#'
#' @details
#' For \eqn{d = 1}, \eqn{\dot S = -S\dot{\mathcal{I}}S}. For \eqn{d = 1/2},
#' \eqn{X = \mathcal{I}^{1/2}} is the symmetric root, whose derivative solves
#' the Sylvester equation \eqn{X\dot X + \dot X X = \dot{\mathcal{I}}}; in
#' the eigenbasis of \eqn{\mathcal{I}} that is a division by
#' \eqn{\sqrt{\lambda_i} + \sqrt{\lambda_j}}, and then
#' \eqn{\dot S = -S\dot X S}.
#'
#' @param I The block.
#' @param dI A list of its derivatives.
#' @param d The power.
#'
#' @return A list with `S` and `dS`, a list like `dI`.
#'
#' @keywords internal
joint_scale_matrix <- function(I, dI, d) {
  K <- nrow(I)
  if (d == 0) {
    return(list(S = diag(K), dS = lapply(dI, function(x) matrix(0, K, K))))
  }
  if (d == 1) {
    S <- solve(I)
    return(list(S = S, dS = lapply(dI, function(x) -S %*% x %*% S)))
  }
  ev <- eigen(I, symmetric = TRUE)
  if (any(ev$values <= 0)) {
    stop("the expected information of the filtered parameters is not positive definite.",
         call. = FALSE)
  }
  Q <- ev$vectors
  r <- sqrt(ev$values)
  S <- Q %*% (t(Q) / r)
  den <- outer(r, r, `+`)
  dS <- lapply(dI, function(x) {
    dX <- Q %*% ((t(Q) %*% x %*% Q) / den) %*% t(Q)
    -S %*% dX %*% S
  })
  list(S = S, dS = dS)
}


#' The Driving Vector of a Joint Filter, One Observation at a Time
#'
#' @description
#' Returns the callbacks that [modelterms7::term_filter()] reads for a
#' [modelterms7::GasJoint()]: the driving vector
#' \eqn{u = \mathcal{I}_{GG}^{-d}\,\partial\ell/\partial\eta_G} at the
#' filtered predictors of one row, the other predictors held at their static
#' values, and its derivative in the filtered predictors.
#'
#' @details
#' With \eqn{g} the link-scale gradient, \eqn{H} the Hessian and
#' \eqn{S = \mathcal{I}_{GG}^{-d}},
#' \deqn{\frac{\partial u_l}{\partial \eta_k} = \sum_c S_{lc} H_{ck}
#'   + \sum_c \frac{\partial S_{lc}}{\partial\eta_k}\, g_c .}
#' The two callbacks are called in turn at the same point, so the family is
#' evaluated once per point.
#'
#' Where the response is a vector and the C registries of
#' \pkg{distributions7} and \pkg{linkfunctions7} know the family and every
#' link, the result also carries `fast`, the context with which
#' [modelterms7::term_filter()] runs the recursion compiled and reads the
#' family through the joint registry; a family that the joint registry does
#' not cover leaves that context inert, and the callbacks run.
#'
#' @param spec A [StatmodSpec()].
#' @param tm The joint filter.
#' @param eta The static predictors, a named list over the distribution
#'   parameters.
#'
#' @return A list with `score` and `curvature`, functions `(e, i)`,
#'   `score_y` and `curvature_y`, the same at a response `y` given as their
#'   first argument, and `fast`, the compiled route's context or `NULL`.
#'
#' @keywords internal
joint_callbacks <- function(spec, tm, eta) {
  params <- spec@distrib@params
  G <- match(tm@params, params)
  K <- length(G)
  n <- spec@n_obs
  d <- tm@terms[[1L]]@scaling
  fam <- joint_row_family(spec, info = d != 0)
  E <- vapply(params, function(p) rep_len(as.numeric(eta[[p]]), n),
              numeric(n))
  E <- matrix(E, n, length(params))
  last <- NULL
  # the response at the row may be the caller's, which is how a simulation
  # reads the value it has just drawn
  at <- function(e, i, yi = NULL) {
    if (is.null(yi) && !is.null(last) && last$i == i &&
        identical(last$e, e)) return(last)
    full <- E[i, ]
    full[G] <- e
    r <- if (is.null(yi)) fam(i, full) else fam(i, full, yi)
    gG <- r$g[G]
    if (d == 0) {
      u <- gG
      C <- r$H[G, G, drop = FALSE]
    } else {
      sm <- joint_scale_matrix(r$I[G, G, drop = FALSE],
                               lapply(G, function(k) r$dI[G, G, k]), d)
      u <- as.numeric(sm$S %*% gG)
      C <- sm$S %*% r$H[G, G, drop = FALSE]
      for (k in seq_len(K)) C[, k] <- C[, k] + as.numeric(sm$dS[[k]] %*% gG)
    }
    out <- list(i = i, e = e, u = u, C = C)
    if (is.null(yi)) last <<- out
    out
  }
  list(score = function(e, i) at(e, i)$u,
       curvature = function(e, i) at(e, i)$C,
       score_y = function(y, e, i) at(e, i, y)$u,
       curvature_y = function(y, e, i) at(e, i, y)$C,
       fast = joint_fast_context(spec, G, d, E))
}


#' The Compiled Route's Context of a Joint Filter
#'
#' @description
#' Returns the list that [modelterms7::term_filter()] reads to run a joint
#' filter compiled, or `NULL` where the response is a matrix or the scalar
#' registries do not know the family or one of the links.
#'
#' @details
#' [distributions7::distrib_scalar_route()] gives the name under which the
#' registry knows the family and its constants, and
#' [linkfunctions7::link_scalar_route()] the name and the parameters of each
#' link. Whether the family's joint entries are written out is read by the
#' compiled code, which leaves the context inert where they are not.
#'
#' @param spec A [StatmodSpec()].
#' @param G The positions of the filtered parameters.
#' @param d The scaling.
#' @param E The static predictors of every parameter, one column each.
#'
#' @return A list, or `NULL`.
#'
#' @keywords internal
joint_fast_context <- function(spec, G, d, E) {
  dist <- spec@distrib
  y <- spec@response
  if (is.matrix(y)) return(NULL)
  route <- distributions7::distrib_scalar_route(dist)
  if (is.null(route)) return(NULL)
  lks <- dist@link_params[dist@params]
  lroutes <- lapply(lks, linkfunctions7::link_scalar_route)
  if (any(vapply(lroutes, is.null, logical(1)))) return(NULL)
  n <- spec@n_obs
  list(family = route$name,
       constants = lapply(route$constants, function(v) {
         v <- as.numeric(v)
         if (length(v) == 1L) v else rep_len(v, n)
       }),
       links = vapply(lroutes, function(r) r$name, character(1)),
       link_par = lapply(lroutes, function(r) r$par),
       bounds = lapply(lks, function(lk) as.numeric(lk@link_bounds)),
       y = as.numeric(y), eta = E, filtered = as.integer(G), scaling = d)
}


#' The Driving Vector's Derivatives in Every Predictor, at the Fit
#'
#' @description
#' Returns the data that [modelterms7::term_curvature()] and its siblings
#' read for a joint filter: the driving vector at the fitted predictors and
#' its derivatives in the predictors of every equation, to the order that
#' `order` gives.
#'
#' @details
#' With `scaling` zero the driving vector is the link-scale gradient of the
#' filtered parameters, and its derivatives are the family's own: the entry
#' \eqn{U_{l;a_1\dots a_r}} is the derivative of order \eqn{r + 1} with
#' indices \eqn{l, a_1, \dots, a_r}. Otherwise [joint_drive_scaled()]
#' composes them with the derivatives of the scaling matrix.
#'
#' @param spec A [StatmodSpec()].
#' @param tm The joint filter.
#' @param gl,H,D3,D4,D5 The family's link-scale derivatives at the fitted
#'   predictors; those above the order that `order` gives may be `NULL`.
#' @param order The highest order of the predictor's derivative to be
#'   propagated: 1 for the reverse pass, 2 for the curvature, 3 for the third
#'   derivative, 4 for the fourth.
#' @param theta The per-observation parameters at the fitted predictors,
#'   from which the expected information is read where `scaling` is not
#'   zero.
#'
#' @return A list with `filtered`, `u`, `U1`, `U2` and, to the order given,
#'   `U3` and `U4`.
#'
#' @keywords internal
joint_drive_data <- function(spec, tm, gl, H, D3 = NULL, D4 = NULL, D5 = NULL,
                             order = 2L, theta = NULL) {
  params <- spec@distrib@params
  P <- length(params)
  G <- match(tm@params, params)
  K <- length(G)
  n <- spec@n_obs
  if (tm@terms[[1L]]@scaling != 0) {
    return(joint_drive_scaled(spec, tm, theta, gl, H, D3, D4, D5, order))
  }
  rn <- function(x) rep_len(as.numeric(x), n)
  u <- vapply(G, function(l) rn(gl[[params[l]]]), numeric(n))
  U1 <- array(0, c(n, K, P))
  U2 <- if (order >= 2L) array(0, c(n, K, P, P))
  U3 <- if (order >= 3L) array(0, c(n, K, P, P, P))
  U4 <- if (order >= 4L) array(0, c(n, K, P, P, P, P))
  for (l in seq_len(K)) {
    L <- G[l]
    for (a in seq_len(P)) {
      U1[, l, a] <- rn(H[[hess_key(params, L, a)]])
      if (order < 2L) next
      for (b in seq_len(P)) {
        U2[, l, a, b] <- rn(D3[[deriv3_key(params, L, a, b)]])
        if (order < 3L) next
        for (c in seq_len(P)) {
          U3[, l, a, b, c] <- rn(D4[[deriv4_key(params, L, a, b, c)]])
          if (order < 4L) next
          for (e in seq_len(P)) {
            U4[, l, a, b, c, e] <- rn(D5[[deriv5_key(params, L, a, b, c, e)]])
          }
        }
      }
    }
  }
  list(filtered = G, u = matrix(u, n, K), U1 = U1, U2 = U2, U3 = U3,
       U4 = U4)
}


#' Run a Joint Filter at the Current Parameters
#'
#' @description
#' Returns the entry of [statmod_filter_at()] for a joint filter: the
#' filtered predictors of every filtered equation, their Jacobian in the
#' filter's parameters, and the callbacks that produced them.
#'
#' @param spec A [StatmodSpec()].
#' @param u The structural unit.
#' @param tm The joint filter.
#' @param psi Its parameters on the parameter scale, as [structural_psi()]
#'   returns them.
#' @param eta_static The static predictors, a named list.
#'
#' @return A list with `param` (the first filtered parameter), `params`,
#'   `joint` (`TRUE`), `term`, `tm`, `eta` (a matrix with one column per
#'   filtered parameter), `jacobian` (a named list of matrices), `curv`,
#'   `psi`, `cb` and `eta_static` (a matrix like `eta`).
#'
#' @keywords internal
joint_filter_entry <- function(spec, u, tm, psi, eta_static) {
  n <- spec@n_obs
  es <- vapply(tm@params, function(p) rep_len(as.numeric(eta_static[[p]]), n),
               numeric(n))
  es <- matrix(es, n, length(tm@params), dimnames = list(NULL, tm@params))
  cb <- joint_callbacks(spec, tm, eta_static)
  f <- modelterms7::term_filter(tm, es, spec@response, cb$score,
                                cb$curvature, psi, fast = cb$fast,
                                threads = spec@threads)
  list(param = tm@params[1L], params = tm@params, joint = TRUE,
       term = u$term, tm = tm, eta = f$eta, jacobian = f$jacobian,
       curv = f$curv, psi = psi, cb = cb, eta_static = es)
}


#' A Joint Filter's Second-Order Recursion at the Fitted Point
#'
#' @description
#' Runs [modelterms7::term_curvature()] for a joint filter at the parameters
#' at which it was fitted and returns the Jacobian of every filtered
#' predictor beside the second derivative contracted against the weighted
#' score.
#'
#' @details
#' The weights are the link-scale gradient of each filtered parameter times
#' the observation weights, so the contraction is the term of the observed
#' information that the curvature of the levels contributes; the driving
#' vector's derivatives are read at the fitted predictors.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param f The filter's entry of [statmod_eta()]'s `filters`.
#' @param Vs The static rows, one matrix per distribution parameter over the
#'   coefficients followed by the filter's parameters.
#' @param gl,H,D3 The family's first three link-scale derivatives at the
#'   fitted predictors.
#' @param theta The per-observation parameters at the fitted predictors.
#'
#' @return The list that [modelterms7::term_curvature()] returns.
#'
#' @keywords internal
joint_filter_curvature <- function(spec, design, f, Vs, gl, H, D3, theta) {
  params <- spec@distrib@params
  n <- spec@n_obs
  w <- spec@weights
  gw <- vapply(f$params, function(p) w * rep_len(as.numeric(gl[[p]]), n),
               numeric(n))
  gw <- matrix(gw, n, length(f$params))
  bd <- joint_drive_data(spec, f$tm, gl, H, D3, order = 2L, theta = theta)
  seed <- stats::setNames(Vs, params)
  structural_memo(design, "curv", list(zeta = f$psi, eta = f$eta_static,
                                       g = gw), function() {
    modelterms7::term_curvature(f$tm, f$eta_static, NULL, NULL, NULL, f$psi, gw,
                                seed, bd)
  })
}


#' The Equation and the Reported Name of a Structural Term's Parameters
#'
#' @description
#' Returns, for names of a structural term's parameters or of its reported
#' quantities, the equation to which each belongs and the name under which
#' it is reported: the term's label, a dot, and the name.
#'
#' @details
#' A term of one equation belongs to that equation. A joint filter spans
#' several, and its names carry the parameter of their term as a prefix
#' (`sigma.kappa1`): the prefix gives the equation and is replaced by the
#' label of that term, so the loading of the scale's filter is reported as
#' `gas.kappa1` among the coefficients of `sigma`, as it would be without
#' the joint filter.
#'
#' @param tm The structural term.
#' @param param The equation in which the term is kept.
#' @param nm The names.
#'
#' @return A list with `param` and `name`, character vectors like `nm`.
#'
#' @keywords internal
structural_display <- function(tm, param, nm) {
  if (is_joint_filter(tm)) {
    eq <- character(length(nm))
    out <- character(length(nm))
    for (k in tm@params) {
      hit <- startsWith(nm, paste0(k, "."))
      eq[hit] <- k
      lb <- tm@terms[[k]]@label
      rest <- substring(nm[hit], nchar(k) + 2L)
      out[hit] <- if (nzchar(lb)) paste(lb, rest, sep = ".") else rest
    }
    return(list(param = eq, name = out))
  }
  lb <- tryCatch(tm@label, error = function(e) "")
  list(param = rep(param, length(nm)),
       name = if (length(lb) == 1L && nzchar(lb)) paste(lb, nm, sep = ".")
              else nm)
}


#' The Joint Filters That Span an Equation
#'
#' @description
#' Returns the joint filters of a specification that filter the predictor of
#' a parameter, wherever the filter is kept, with the key under which each is
#' kept.
#'
#' @param spec A [StatmodSpec()].
#' @param p The distribution parameter.
#'
#' @return A list of entries with `key`, `param` (the equation in which the
#'   filter is kept) and `term`.
#'
#' @keywords internal
joint_filters_of <- function(spec, p) {
  out <- list()
  for (q in names(spec@terms)) {
    for (nm in names(spec@terms[[q]])) {
      tm <- spec@terms[[q]][[nm]]
      if (is_joint_filter(tm) && p %in% tm@params) {
        out[[length(out) + 1L]] <- list(key = nm, param = q, term = tm)
      }
    }
  }
  out
}


#' A Term of an Equation, a Joint Filter Included
#'
#' @description
#' Returns the term kept under a key in the equation of a parameter, or the
#' joint filter of that key that spans the parameter where the filter is
#' kept in another equation.
#'
#' @param spec A [StatmodSpec()].
#' @param p The distribution parameter.
#' @param nm The key.
#'
#' @return The term, or `NULL`.
#'
#' @keywords internal
term_of_equation <- function(spec, p, nm) {
  tm <- spec@terms[[p]][[nm]]
  if (!is.null(tm)) return(tm)
  for (jf in joint_filters_of(spec, p)) if (identical(jf$key, nm)) return(jf$term)
  NULL
}
