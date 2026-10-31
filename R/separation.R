#' Directions Along Which the Estimate Does Not Exist
#'
#' @description
#' Searches the coefficients of a fit for a direction of recession, along
#' which the penalized log-likelihood keeps increasing, or stays level, as the
#' linear predictors move without bound. Such a direction exists when the
#' maximum likelihood estimate does not exist (separation in a binary or
#' binomial response, a level of a factor with all its counts at zero, a group
#' without zeros in the equation of a zero-inflation probability), and the fit
#' then stops at a finite point that depends only on its stopping rule.
#'
#' @details
#' The objective is \eqn{F} of [statmod_objective()], at the fitted
#' hyperparameters, with gradient \eqn{g} and Hessian \eqn{K}. The
#' coefficients are split into those that carry no smooth penalty (\eqn{u})
#' and those that do (\eqn{v}), and the latter are profiled out through the
#' Schur complement \eqn{S = K_{uu} - K_{uv} K_{vv}^{-1} K_{vu}}, so that the
#' eigenvalue problem has the size of the unpenalized coefficients whatever
#' the number of random effects. A penalized coefficient moves along each
#' candidate as the profile requires, \eqn{v = -K_{vv}^{-1} K_{vu} u}, and
#' only the null space of a penalty can carry a direction of recession.
#'
#' The candidates are the Newton step at the fitted point,
#' \eqn{\delta = -K^{-1} g}, with each eigenvalue of \eqn{S} floored at
#' \eqn{10^{-14}} times the largest, and the eigenvectors of \eqn{S} against
#' the metric of the linear predictors, \eqn{M = X_u^\top X_u} over the
#' equations, whose eigenvalues fall below `flat` times the largest one. At a
#' point where the fit has drifted along a direction of recession the gradient
#' and the curvature along it vanish together, so the Newton step remains of
#' order one in the linear predictor while the decrease that it predicts is
#' below the stopping rule.
#'
#' Each candidate \eqn{d} is scaled so that the largest change it makes in a
#' linear predictor is one, and \eqn{F} is evaluated at
#' \eqn{\hat\beta + t d} for every `t` in `steps`. A candidate is a direction
#' of recession when \eqn{F} does not rise by more than `tol` times
#' \eqn{\max(1, |F(\hat\beta)|)} at any of those points, nor from one point
#' to the next, and the whole decrease stays within [mode_error_limit()]. It
#' is a direction of descent when \eqn{F} falls by more than
#' [mode_error_limit()] at one of them, which means that the point is not a
#' minimum at all (a scale whose group has a constant response, where the
#' log-likelihood is unbounded, is an example). At an ordinary optimum a
#' change of one in a linear predictor already raises \eqn{F}.
#'
#' Coefficients held at the kink of a penalty, those of a block that is a
#' working linearization ([frozen_block()]) and aliased coefficients are held
#' fixed. Coefficients estimated on the outer criterion move with the others,
#' since a direction of recession of \eqn{F} is one of that criterion too
#' when no random effect integrates over them, but their gradient in \eqn{F}
#' is set to zero, because they sit at the root of the outer criterion and not
#' of \eqn{F}. For the same reason a direction of descent is read on the change
#' of \eqn{F} after removing the linear term of those coordinates,
#' \eqn{t\,g_m^\top d}, where \eqn{g_m} is their gradient in \eqn{F}: under
#' [reml()] the scale of a gaussian is not at the minimum of \eqn{F}, and
#' without the correction the difference between the two criteria reads as a
#' direction of descent. A direction of recession is read on \eqn{F} itself.
#'
#' @param fit A [StatmodFit()].
#' @param steps The multiples of the scaled direction at which the objective is
#'   evaluated, in units of the linear predictor.
#' @param tol The largest rise of the objective, relative to
#'   \eqn{\max(1, |F(\hat\beta)|)}, that is still read as level.
#' @param flat The largest ratio of an eigenvalue to the largest one at which
#'   an eigenvector of \eqn{S} against \eqn{M} is probed beside the Newton
#'   step.
#' @param max_eigen The largest number of eigenvectors that are probed.
#'
#' @return A list with `found` (logical, a direction of recession), `descent`
#'   (logical, a direction of descent), `directions` (one entry per direction
#'   found, each with `kind`, `"recession"` or `"descent"`; `coefficients`,
#'   the named components of the direction over the coefficients that it
#'   moves, scaled so that the largest change in a linear predictor is one;
#'   `free`, the names of those coefficients that carry no smooth penalty;
#'   `equations`; `n_moved`, the number of observations whose linear predictor
#'   moves by more than 1e-3 in each equation; `change`, the change of the
#'   objective at each step; `source`, `"newton"` or `"eigen"`; and
#'   `curvature`, the curvature along the direction against \eqn{M} as a ratio
#'   to the largest eigenvalue), `newton_eta`, the largest change that the
#'   Newton step makes in a linear predictor, `min_eigen`, the smallest ratio
#'   of an eigenvalue of \eqn{S} to the largest one, `n_eval`, the number of
#'   evaluations of the objective, and `reason`, a character string where the
#'   check could not be made.
#'
#' @keywords internal
separation_diagnosis <- function(fit, steps = 2^(0:6), tol = 1e-8,
                                 flat = 1e-6, max_eigen = 4L) {
  out <- list(found = FALSE, descent = FALSE, directions = list(),
              newton_eta = NA_real_, min_eigen = NA_real_, n_eval = 0L,
              reason = character(0))
  spec <- fit@spec
  design <- tryCatch(statmod_design(spec), error = function(e) NULL)
  if (is.null(design)) {
    out$reason <- "the design could not be rebuilt from the fit"
    return(out)
  }
  if (length(fit@methods$pinned)) attr(design, "pinned") <- fit@methods$pinned
  params <- spec@distrib@params
  obj <- tryCatch(statmod_objective(spec, fit@hyper, design, FALSE, NULL),
                  error = function(e) NULL)
  if (is.null(obj)) {
    out$reason <- "the objective could not be built"
    return(out)
  }
  b <- obj$stack(fit@coefficients)
  p <- length(b)
  F0 <- obj$fn(b)
  g <- tryCatch(obj$gr(b), error = function(e) NULL)
  K <- tryCatch(obj$he(b), error = function(e) NULL)
  if (!is.finite(F0) || is.null(g) || is.null(K) || !all(is.finite(g))) {
    out$reason <- "the objective, its gradient or its Hessian is not finite"
    return(out)
  }

  # the coordinates that move, and which of them carry a smooth penalty
  lab <- coef_labels(spec, design)
  if (nrow(lab) != p) {
    out$reason <- "the coefficients could not be labeled"
    return(out)
  }
  drop <- (lab$kinked & b == 0) | frozen_block(spec, lab) |
    rownames(lab) %in% fit@aliased
  iu <- which(!drop & !(lab$penalized & !lab$kinked))
  iv <- which(!drop & lab$penalized & !lab$kinked)
  if (!length(iu)) {
    out$reason <- "every coefficient that moves carries a smooth penalty"
    return(out)
  }
  # a coefficient estimated on the outer criterion sits at the root of that
  # criterion and not of F, so its gradient in F is not a residual of the
  # fit: it is set to zero for the Newton step and the slopes, as the
  # certificate sets it to zero for the mode error
  gam <- tryCatch(marginal_coords(spec, design, fit@methods$outer),
                  error = function(e) NULL)
  g_marg <- numeric(p)
  if (!is.null(gam) && length(gam$where)) {
    g_marg[gam$where] <- g[gam$where]
    g[gam$where] <- 0
  }

  npar <- obj$npar
  offs <- cumsum(npar) - npar
  Xs <- lapply(params, function(a) design[[a]]$X)
  delta_eta <- function(d) {
    lapply(seq_along(params), function(a) {
      if (!npar[a]) return(numeric(0))
      as.numeric(as.matrix(Xs[[a]] %*% d[offs[a] + seq_len(npar[a])]))
    })
  }
  eq_of <- rep(params, npar)

  Kuu <- as_dense(K[iu, iu, drop = FALSE])
  if (length(iv)) {
    Kvv <- K[iv, iv, drop = FALSE]
    Kvu <- as_dense(K[iv, iu, drop = FALSE])
    if (!all(is.finite(as_dense(Kvv))) || !all(is.finite(Kvu))) {
      out$reason <- "the Hessian is not finite"
      return(out)
    }
    fac <- tryCatch(if (isS4(Kvv)) Matrix::Cholesky(Matrix::forceSymmetric(Kvv))
                    else chol(Kvv), error = function(e) NULL)
    if (is.null(fac)) {
      out$reason <- "the penalized block of the Hessian could not be factorized"
      return(out)
    }
    solve_v <- function(r) {
      as.matrix(if (isS4(fac)) Matrix::solve(fac, r)
                else backsolve(fac, backsolve(fac, r, transpose = TRUE)))
    }
    KvvKvu <- solve_v(Kvu)
    S <- Kuu - crossprod(Kvu, KvvKvu)
  } else {
    S <- Kuu
  }
  if (!all(is.finite(S))) {
    out$reason <- "the Hessian is not finite"
    return(out)
  }
  S <- (S + t(S)) / 2
  # a direction over u, completed by the profile over v
  complete <- function(du) {
    d <- numeric(p)
    d[iu] <- du
    if (length(iv)) d[iv] <- -as.numeric(KvvKvu %*% du)
    d
  }

  # the metric of the linear predictors over u
  M <- matrix(0, length(iu), length(iu))
  for (a in seq_along(params)) {
    w <- which(eq_of[iu] == params[a])
    if (!length(w)) next
    Xa <- as_dense(Xs[[a]][, iu[w] - offs[a], drop = FALSE])
    M[w, w] <- crossprod(Xa)
  }
  em <- eigen(M, symmetric = TRUE)
  ok <- em$values > max(em$values) * 1e-10
  if (!any(ok)) {
    out$reason <- "the linear predictors do not move with the coefficients"
    return(out)
  }
  W <- em$vectors[, ok, drop = FALSE] %*%
    diag(1 / sqrt(em$values[ok]), sum(ok))
  C <- crossprod(W, S %*% W)
  ec <- eigen((C + t(C)) / 2, symmetric = TRUE)
  lmax <- max(abs(ec$values))
  rel <- ec$values / lmax
  out$min_eigen <- min(rel)

  scale_to_eta <- function(d) {
    m <- max(abs(unlist(delta_eta(d))))
    if (!is.finite(m) || m == 0) return(NULL)
    d / m
  }

  # the Newton step on the profile, with the curvature floored at rounding of
  # the largest eigenvalue so that a direction whose curvature is lost to
  # rounding is still followed
  gu <- g[iu]
  if (length(iv)) gu <- gu - as.numeric(crossprod(KvvKvu, g[iv]))
  lam <- pmax(ec$values, lmax * 1e-14)
  du <- -as.numeric(W %*% (ec$vectors %*%
                              (crossprod(ec$vectors, crossprod(W, gu)) / lam)))
  newton <- numeric(p)
  newton[iu] <- du
  if (length(iv)) {
    newton[iv] <- -as.numeric(solve_v(g[iv] + as.numeric(Kvu %*% du)))
  }
  out$newton_eta <- max(abs(unlist(delta_eta(newton))))
  cand <- list()
  dn <- scale_to_eta(newton)
  if (!is.null(dn)) cand[[1]] <- list(d = dn, source = "newton")
  small <- order(rel)[seq_len(min(max_eigen, length(rel)))]
  small <- small[rel[small] <= flat]
  for (k in small) {
    v <- scale_to_eta(complete(as.numeric(W %*% ec$vectors[, k])))
    if (is.null(v)) next
    slope <- sum(g * v)
    signs <- if (abs(slope) > 0) -sign(slope) else c(1, -1)
    for (s in signs) cand[[length(cand) + 1L]] <-
      list(d = s * v, source = "eigen")
  }

  # the probe
  lev <- tol * max(1, abs(F0))
  fall <- mode_error_limit()
  n_eval <- 0L
  found <- list()
  for (cd in cand) {
    d <- cd$d
    dup <- any(vapply(found, function(f) {
      abs(sum(f$d * d)) / sqrt(sum(f$d^2) * sum(d^2)) > 0.999
    }, logical(1)))
    if (dup) next
    change <- corrected <- rep(NA_real_, length(steps))
    for (i in seq_along(steps)) {
      Ft <- obj$fn(b + steps[i] * d)
      n_eval <- n_eval + 1L
      change[i] <- Ft - F0
      corrected[i] <- change[i] - steps[i] * sum(g_marg * d)
      if (is.nan(change[i]) || change[i] > lev) break
    }
    seen <- change[!is.na(change)]
    if (!length(seen) || any(is.nan(seen))) next
    kind <- if (min(corrected, na.rm = TRUE) < -fall) {
      "descent"
    } else if (length(seen) == length(steps) && all(seen <= lev) &&
               all(diff(seen) <= lev)) {
      "recession"
    } else next
    du <- d[iu]
    curv <- as.numeric(crossprod(du, S %*% du)) /
      as.numeric(crossprod(du, M %*% du)) / lmax
    found[[length(found) + 1L]] <- list(d = d, change = change, kind = kind,
                                        source = cd$source, curvature = curv)
  }
  out$n_eval <- n_eval
  if (!length(found)) return(out)

  colscale <- numeric(p)
  for (a in seq_along(params)) {
    if (!npar[a]) next
    colscale[offs[a] + seq_len(npar[a])] <-
      sqrt(Matrix::colSums(Xs[[a]]^2))
  }
  out$found <- any(vapply(found, function(f) f$kind == "recession", TRUE))
  out$descent <- any(vapply(found, function(f) f$kind == "descent", TRUE))
  out$directions <- lapply(found, function(f) {
    contrib <- abs(f$d) * colscale
    who <- which(contrib > 1e-3 * max(contrib))
    de <- delta_eta(f$d)
    moved <- vapply(de, function(z) sum(abs(z) > 1e-3), integer(1))
    names(moved) <- params
    list(kind = f$kind,
         coefficients = stats::setNames(f$d[who], rownames(lab)[who]),
         free = rownames(lab)[intersect(who, iu)],
         equations = unique(eq_of[who]), n_moved = moved[moved > 0],
         change = f$change, source = f$source, curvature = f$curvature)
  })
  out
}


#' A Sentence for Each Direction That separation_diagnosis() Found
#'
#' @description
#' Describes each direction of recession and each direction of descent in
#' the result of [separation_diagnosis()] in one sentence, naming the
#' coefficients that move and the number of observations whose linear
#' predictor moves with them.
#'
#' @param sep The list that [separation_diagnosis()] returns.
#'
#' @return A character vector with one sentence per direction, empty where
#'   there is none.
#'
#' @keywords internal
separation_message <- function(sep) {
  if (!length(sep$directions)) return(character(0))
  vapply(sep$directions, function(d) {
    cf <- d$coefficients
    where <- if (length(cf) == 1L) {
      sprintf("as %s moves toward %s infinity", names(cf),
              if (cf > 0) "plus" else "minus")
    } else {
      sprintf("along the direction %s",
              paste(sprintf("%.3g %s", cf, names(cf)), collapse = " + "))
    }
    obs <- paste(sprintf("%s at %d observations", names(d$n_moved),
                         d$n_moved), collapse = " and ")
    if (identical(d$kind, "recession")) {
      sprintf(paste0(
        "the log-likelihood increases without reaching a maximum %s, which ",
        "moves the linear predictor of %s, so the estimate does not exist ",
        "and the value reported depends on where the fit stopped"),
        where, obs)
    } else {
      sprintf(paste0(
        "the point is not a minimum of the objective, which decreases by ",
        "%.4g %s; this moves the linear predictor of %s"),
        -min(d$change, na.rm = TRUE), where, obs)
    }
  }, character(1))
}
