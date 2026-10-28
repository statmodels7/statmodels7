# The restricted fit, and what the four likelihood statistics read off it.
#
# Three of the four statistics -- the likelihood ratio, Rao's score and
# Terrell's gradient -- are functions of the model refitted with one
# coefficient HELD at a value. Wald's is not: it reads the unrestricted fit
# alone, which is exactly why it is the one that is not invariant to how the
# parameter is written.


#' Refit With One Coefficient Held at a Value
#'
#' @description
#' Maximizes the penalized likelihood over every coefficient but one, which
#' is held at `value`. It returns the restricted maximum, the coefficients
#' that attain it and the score there.
#'
#' @details
#' The hold is written on the specification and enforced in the solve, so the
#' refit does not rewrite the formula and does not rebuild any term: a design
#' whose bases were rebuilt would carry different knots and would correspond
#' to a different model. See
#' [held_positions()] and [iwls_solve()] for the two halves.
#'
#' # The hyperparameters are held too
#'
#' The refit runs the alternation once at the hyperparameters of the
#' unrestricted fit, with no outer search. The statistics built on it are
#' therefore conditional on the smoothing that the data chose once. That is
#' the reading that a penalized fit supports, because the penalty is a
#' log-prior at the hyperparameters in force, and it keeps an interval by
#' inversion affordable: re-selecting the smoothing at every held value would
#' put a whole outer search inside every step of a root find. A profiled
#' version requires refitting by hand.
#'
#' # Coefficients that cannot be held
#'
#' A coefficient under a kinked penalty (a lasso, a SCAD, an MCP) cannot be
#' held, for two reasons. The block is fitted by a coordinate descent, which
#' updates every coordinate of the block and cannot leave one at a held
#' value, so that [coord_fit()] signals an error when a held coordinate lies
#' in its block. And at the kink the objective has no curvature for the
#' statistics to read, while the null distribution of a coefficient that a
#' selection kept is a different problem. An aliased coefficient has no
#' estimate to test. A model that carries a structural term cannot hold a
#' coefficient either, because the contribution of the term is the state of a
#' recursion and not \eqn{X\beta}, and the hold does not reach it.
#'
#' A coefficient under a penalty that is twice differentiable is held like
#' any other: a smooth's linear column and its rotated coordinates, a ridge,
#' a random effect. [statmod_stat_at()] states the meaning of such a fit.
#'
#' # Break-points
#'
#' A coefficient of a continuous [modelterms7::seg()] term with two or more
#' break-points is held with the other break-points free, and the fit is the
#' one of a break-point at the held value. A free break-point cannot reach the
#' held one: it stays on one side of it with at least two observations
#' between the two ([modelterms7::seg_pin()]), because two break-points at
#' one position make two columns equal and the changes of slope grow without
#' bound in opposite directions. The held position stands for each of the
#' break-points of the term in turn, taking the change of slope of that
#' break-point while the others keep their positions and changes of slope, so
#' that the free break-point starts on one side of the held position or on the
#' other, and the refit with the smaller objective is kept. The coefficients
#' of the result are in the order of the positions, and `[held]` follows the
#' names of the held break-point's coefficients ([pin_held_breakpoints()]).
#'
#' @param fit A [StatmodFit()].
#' @param param The distribution parameter whose equation carries the
#'   coefficient, a single string.
#' @param coefname The coefficient's name, a single string, as
#'   [coef.StatmodFit()] reports it.
#' @param value The value at which to hold it, a single number.
#'
#' @return A list:
#'   \describe{
#'     \item{`coefficients`}{the restricted estimates, a named list.}
#'     \item{`loglik`}{the log-likelihood there, unpenalized.}
#'     \item{`objective`}{the penalized objective there, the negative
#'       log-likelihood plus the penalty, which the fit minimizes.}
#'     \item{`score`}{the penalized score at the restricted point, a numeric
#'       vector over the stacked coefficients. Its free entries are zero to
#'       the tolerance; the statistics read the held one.}
#'     \item{`at`}{the held coordinate's position in that vector.}
#'     \item{`par`}{the restricted estimates, stacked.}
#'     \item{`information`}{a function of no arguments returning the
#'       penalized information at the restricted point. It is read at most
#'       once however often it is called, and it is not read if the function
#'       is never called, which is the ordinary case: the likelihood ratio and
#'       the gradient statistic need no matrix.}
#'     \item{`mode_error`}{a function of no arguments returning
#'       [restricted_mode_error()], how far above its own mode the refit
#'       stopped. It reads the information, so it costs one Hessian.}
#'     \item{`labels`}{the stacked coefficient labels.}
#'     \item{`converged`}{the inner optimizer's flag, a single logical. It
#'       records whether a stopping rule fired and does not measure the
#'       distance to the mode, which [restricted_mode_error()] measures.}
#'   }
#'
#' @seealso [iwls_solve()], which drops the held coordinate from the system.
#'
#' @keywords internal
statmod_restrict <- function(fit, param, coefname, value) {
  # a value of another length reaches the hold as a named vector of that
  # length and fails several frames down, inside the objective, on a names
  # assignment -- which names neither the argument nor the mistake
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("'value' must be a single finite number.", call. = FALSE)
  }
  spec0 <- fit@spec
  params <- spec0@distrib@params
  if (!param %in% params) {
    stop(sprintf("'%s' is not a parameter of %s. Its parameters are %s.",
                 param, spec0@distrib@distrib_name,
                 paste(params, collapse = ", ")), call. = FALSE)
  }
  design0 <- statmod_design(spec0)
  if (length(attr(design0, "structural"))) {
    stop("A model carrying a structural term cannot hold a coefficient: ",
         "its\n  contribution is a recursion's state and not a column times ",
         "a coefficient.", call. = FALSE)
  }
  nms <- design0[[param]]$coef_names
  j <- match(coefname, nms)
  if (is.na(j)) {
    stop(sprintf("'%s' names no coefficient of the '%s' equation.",
                 coefname, param), call. = FALSE)
  }
  npar <- vapply(design0, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  at <- offs[[match(param, params)]] + j
  if (aliased_labels(spec0, design0, at) %in% fit@aliased) {
    stop(sprintf(paste0("'%s' is aliased in the '%s' equation, so it has no ",
                        "estimate to test:\n  its column is a combination ",
                        "of the others and the fit left it out."),
                 coefname, param), call. = FALSE)
  }
  if (at %in% kinked_coords(spec0, design0)) {
    stop(sprintf(paste0("'%s' is under a kinked penalty (a lasso, a SCAD ",
                        "or an MCP) and\n  cannot be held. That block is ",
                        "fitted by a coordinate descent, which\n  updates ",
                        "every coordinate of the block and cannot leave one ",
                        "at a held value,\n  and at the kink the objective ",
                        "has no curvature for a statistic to read.\n  A ",
                        "penalty that is twice differentiable (a ridge, a ",
                        "random effect,\n  a smooth) is held like any other ",
                        "coefficient."), coefname), call. = FALSE)
  }

  hc <- spec0@held_coef
  hc[[param]] <- stats::setNames(as.numeric(value), coefname)
  spec <- S7::set_props(spec0, held_coef = hc)
  method <- fit@methods$smooth
  cfg <- inner_settings(method)
  beta_fit <- unlist(fit@coefficients[params], use.names = FALSE)
  beta0 <- beta_fit
  beta0[at] <- as.numeric(value)
  # ONE FIT FOR EACH BREAK-POINT THE HELD ONE CAN STAND FOR. The held
  # position takes the change of slope of one of the unrestricted
  # break-points and the free ones keep the rest, so that the free
  # break-point starts on one side of the held one or on the other. The fit
  # that ends with the smaller objective is kept.
  run <- function(variant) {
    design <- statmod_design(spec)
    blocks <- statmod_blocks(spec, design)
    pins <- pin_held_breakpoints(spec, design, beta0, from = beta_fit,
                                 variant = variant)
    res <- statmod_alternate(spec, design, blocks, fit@hyper, method,
                             pins$beta, cfg$expected, cfg$emethod, cfg$maxit,
                             cfg$tol, verbosity(0))
    list(design = design, pins = pins, res = res)
  }
  best <- run(1L)
  for (v in seq_len(best$pins$variants)[-1L]) {
    other <- run(v)
    if (is.finite(other$res$value) && other$res$value < best$res$value) {
      best <- other
    }
  }
  design <- best$design
  pins <- best$pins
  res <- best$res
  cf <- res$obj$split(res$par)
  # the objective is a MINIMAND -- the negative log-likelihood plus the
  # penalty -- so its gradient carries the opposite sign to the score
  sc <- -res$obj$gr(res$par)
  # the Hessian is read at most ONCE, and only where a caller asks: the
  # likelihood ratio and the gradient statistic never do. Measured on a
  # Poisson of twenty-six coefficients it is 63 per cent of the refit
  # (0.0044 s against 0.0070); at two coefficients, 2.6 per cent.
  hess <- NULL
  info <- function() {
    # the objective's Hessian IS the penalized information there, the
    # objective being the negative log-likelihood plus the penalty
    if (is.null(hess)) hess <<- as_dense(res$obj$he(res$par))
    hess
  }
  list(coefficients = pins$present(cf),
       loglik = statmod_loglik_at(spec, cf, design),
       objective = res$value,
       score = sc,
       at = at,
       par = res$par,
       information = info,
       mode_error = function() restricted_mode_error(sc, info(), at),
       labels = rownames(coef_labels(spec, design)),
       converged = isTRUE(res$converged))
}


#' Keep the Free Break-Points Away From a Held One
#'
#' @description
#' Prepares the break-point terms of a restricted fit in which a coefficient
#' of a continuous [modelterms7::seg()] term with two or more break-points is
#' held. The free break-points are confined to the side of the held position
#' on which they start, with at least `min_obs` observations between them and
#' it ([modelterms7::seg_pin()]), and the starting values of their positions
#' are brought inside the interval. The position of a free break-point is kept
#' inside its interval by an infinite objective outside it ([pin_infeasible()]).
#'
#' @details
#' Two break-points at one position make two columns of the block equal. A
#' free break-point that reaches the held one then has coefficients of the
#' changes of slope that grow without bound in opposite directions, and the
#' relabeling of crossed break-points of [modelterms7::term_refresh()] writes
#' the held value onto the coefficients of the free one. A restricted fit of
#' this kind was measured at 4006 scoring iterations, without converging, for
#' a statistic that was not the likelihood ratio.
#'
#' The coefficients of such a term are returned in the order of the positions
#' of the break-points, by `present()`, and the coefficient labeled as held
#' carries `[held]` after its name.
#'
#' @param spec A [StatmodSpec()] with the coefficient held.
#' @param design Its design, whose state holds the terms that a fit
#'   refreshes. The terms in the state are replaced.
#' @param beta The starting coefficients, stacked, with the held value in
#'   place.
#' @param from The coefficients of the unrestricted fit, stacked, from which
#'   the starting pairs of position and change of slope are taken.
#' @param min_obs The smallest number of observations between a free and a
#'   held break-point.
#' @param variant Which unrestricted break-point the held position stands
#'   for, counted from the nearest to it: the held position takes the change
#'   of slope of that break-point and the free ones keep the other pairs.
#'
#' @return A list with `beta`, the starting coefficients with the free
#'   positions inside their intervals, `present`, a function of the
#'   coefficient list of the restricted fit that returns it with the
#'   break-points in order, and `variants`, the number of break-points that
#'   the held position can stand for: the number of break-points of the term
#'   when one position is held, and one otherwise. A model without such a
#'   term is left unchanged.
#'
#' @seealso [statmod_restrict()], which calls it.
#'
#' @keywords internal
pin_held_breakpoints <- function(spec, design, beta, from = beta,
                                 min_obs = 2L, variant = 1L) {
  out <- list(beta = beta, present = function(cf) cf, variants = 1L)
  st <- attr(design, "state")
  hc <- spec@held_coef
  if (is.null(st) || !length(hc)) return(out)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  entries <- list()
  for (p in names(hc)) {
    a <- match(p, params)
    for (nm in names(st$terms[[p]])) {
      tm <- st$terms[[p]][[nm]]
      if (!S7::S7_inherits(tm, modelterms7::SegTerm) ||
          !identical(tm@kind, "seg") || length(tm@subformulas) ||
          tm@npsi < 2L) {
        next
      }
      cols <- design[[p]]$blocks[[nm]]
      cn <- design[[p]]$coef_names[cols]
      hv <- hc[[p]][names(hc[[p]]) %in% cn]
      if (!length(hv)) next
      K <- tm@npsi
      bp <- tm@blueprint
      own <- tm@coef_names[match(names(hv), cn)]
      at_of <- function(prefix) {
        vapply(paste0(prefix, seq_len(K)), function(q) bp$index[[q]][1L],
               integer(1))
      }
      psi_cols <- cols[at_of("psi")]
      gam_cols <- cols[at_of("gamma")]
      ix <- offs[[a]] + psi_cols
      held_psi <- match(own, tm@coef_names[at_of("psi")])
      held_gam <- match(own, tm@coef_names[at_of("gamma")])
      held_values <- stats::setNames(as.numeric(hv), own)
      # THE STARTING PAIRS. The held position takes the change of slope of
      # one of the unrestricted break-points, the `variant`-th nearest to it,
      # and the free positions keep the pairs that remain. Left as they were,
      # the held position would carry the change of slope of another feature
      # of the data, and the fit would start from a bad local optimum
      # (measured, the objective at the end of the fit was 178 against 97 with
      # the pairs matched).
      if (all(is.na(held_gam))) {
        gix <- offs[[a]] + gam_cols
        lineage <- cbind(psi = from[ix], gamma = from[gix])
        hslot <- held_psi[!is.na(held_psi)]
        hpos <- as.numeric(hv)[!is.na(held_psi)]
        free_k <- setdiff(seq_len(K), hslot)
        left <- seq_len(K)
        take <- integer(length(hslot))
        for (j in seq_along(hslot)) {
          ranked <- left[order(abs(lineage[left, "psi"] - hpos[j]))]
          take[j] <- ranked[min(if (j == 1L) variant else 1L, length(ranked))]
          left <- setdiff(left, take[j])
        }
        if (length(hslot) == 1L) out$variants <- K
        out$beta[ix[hslot]] <- hpos
        out$beta[gix[hslot]] <- lineage[take, "gamma"]
        out$beta[ix[free_k]] <- lineage[left, "psi"]
        out$beta[gix[free_k]] <- lineage[left, "gamma"]
      }
      start <- out$beta[ix]
      pinned <- modelterms7::seg_pin(tm, held = held_values, start = start,
                                     min_obs = min_obs)
      free_k <- setdiff(seq_len(K), held_psi[!is.na(held_psi)])
      st$terms[[p]][[nm]] <- pinned
      lm <- modelterms7::seg_limits(pinned)
      out$beta[ix] <- start
      narrow <- integer(0)
      for (k in free_k) {
        if (lm[k, 1L] > bp$lim[1L] || lm[k, 2L] < bp$lim[2L]) {
          m <- min(1e-3 * (bp$hi - bp$lo), 0.25 * (lm[k, 2L] - lm[k, 1L]))
          out$beta[ix[k]] <- min(max(start[k], lm[k, 1L] + m), lm[k, 2L] - m)
          narrow <- c(narrow, k)
        }
      }
      # the position stays inside its interval during the fit: a step that
      # takes the coefficient out of it is rejected by the objective
      if (length(narrow)) {
        st$pins <- c(st$pins, list(list(p = p, cols = psi_cols[narrow],
                                        lo = lm[narrow, 1L],
                                        hi = lm[narrow, 2L])))
      }
      held_slot <- seq_len(K) %in% c(held_psi, held_gam)
      entries[[length(entries) + 1L]] <- list(
        p = p, psi = psi_cols, gam = gam_cols, held = held_slot, lim = lm,
        names = design[[p]]$coef_names)
    }
  }
  if (!length(entries)) return(out)
  out$present <- function(cf) {
    for (e in entries) {
      v <- cf[[e$p]]
      # the positions that the fit used: a coefficient beyond its interval is
      # read at the end of the interval
      pos <- pmin(pmax(v[e$psi], e$lim[, 1L]), e$lim[, 2L])
      o <- order(pos)
      w <- v
      w[e$gam] <- v[e$gam[o]]
      w[e$psi] <- pos[o]
      nms <- if (is.null(names(v))) e$names else names(v)
      tag <- ifelse(e$held[o], " [held]", "")
      nms[e$gam] <- paste0(e$names[e$gam], tag)
      nms[e$psi] <- paste0(e$names[e$psi], tag)
      names(w) <- nms
      cf[[e$p]] <- w
    }
    cf
  }
  out
}


#' Whether a Free Break-Point Has Left Its Interval
#'
#' @description
#' Reads the intervals that [pin_held_breakpoints()] recorded in the state of
#' a design and returns `TRUE` when the position of a free break-point lies
#' outside its interval at the coefficients `cf`. The objective is infinite
#' there, so that a line search backtracks to the inside.
#'
#' @param design A design.
#' @param cf The coefficient list, one vector for each parameter.
#'
#' @return A single logical, `FALSE` for a design with no recorded interval.
#'
#' @keywords internal
pin_infeasible <- function(design, cf) {
  st <- attr(design, "state")
  pins <- if (is.null(st)) NULL else st$pins
  if (is.null(pins)) return(FALSE)
  for (e in pins) {
    v <- cf[[e$p]][e$cols]
    if (any(v < e$lo | v > e$hi)) return(TRUE)
  }
  FALSE
}


#' How Far Above Its Mode a Restricted Fit Stopped
#'
#' @description
#' Computes \eqn{\tfrac12 g'K^{-1}g} over the free coordinates of a
#' restricted fit, which is the decrease of the penalized objective that a
#' Newton step from the point returned by [statmod_restrict()] would predict,
#' in log-likelihood units.
#'
#' @details
#' [inner_mode_error()] computes the same quantity inside a search; here it
#' is applied to a restricted fit. The inner optimizer's flag records whether a
#' stopping rule fired, which is a boolean about a threshold on a score whose
#' size depends on the model, whereas whether a point is at its mode is a
#' matter of distance and has a natural scale. A stopping rule reads a change
#' in the objective and this quantity reads the point, so the two need not
#' agree, and the print method of [StatmodTest()] reads the mode error and
#' not the flag.
#'
#' The held coordinate is excluded because it is not free: its score is the
#' quantity that Rao's statistic reads, and it does not vanish under the
#' restriction.
#'
#' @param score The restricted fit's score, all coordinates.
#' @param information Its penalized information, all coordinates.
#' @param at The position of the held coordinate.
#'
#' @return A single number, or `NA` where the score on the free coordinates
#'   is not finite, where there is no free coordinate, or where the
#'   information on the free coordinates cannot be inverted.
#'
#' @seealso [mode_error_limit()], the limit against which it is read, and
#'   [inner_mode_error()], the same reading inside a search.
#'
#' @keywords internal
restricted_mode_error <- function(score, information, at) {
  free <- setdiff(seq_along(score), at)
  if (!length(free) || !all(is.finite(score[free]))) return(NA_real_)
  g <- score[free]
  v <- tryCatch(0.5 * sum(g * solve(information[free, free, drop = FALSE], g)),
                error = function(e) NA_real_)
  if (isTRUE(is.finite(v))) v else NA_real_
}


#' Which Stacked Coordinates a Kinked Penalty Covers
#'
#' @description
#' Returns the positions of the coefficients whose penalty is not
#' differentiable everywhere -- a lasso, a SCAD, an MCP -- as an integer
#' vector.
#'
#' @details
#' The function tests the coordinates in the way that [statmod_blocks()] tests
#' the blocks to decide which of the two fitting schemes is used for a block,
#' and it reads the same predicate ([penalty_has_kink()]), so that the two
#' cannot disagree. A penalty that is twice differentiable contributes
#' nothing, so that a ridge, a random effect and a smooth are testable.
#'
#' The result follows from the enumeration of each penalty, so a partial
#' penalty (one covering some of a term's coordinates and not the rest) names
#' only what it covers, and a term is not classified as a whole. The penalty of
#' a structural term covers positions among the term's own parameters and not
#' among the design's, and contributes nothing here.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#'
#' @return An integer vector of stacked positions, possibly empty.
#'
#' @seealso [statmod_penalized()], the enumeration that it reads, and
#'   [statmod_blocks()], which gives the same classification for a whole
#'   block.
#'
#' @keywords internal
kinked_coords <- function(spec, design) {
  out <- integer(0)
  for (u in statmod_penalized(spec, design)) {
    if (isTRUE(u$structural) || is.null(u$index)) next
    kink <- tryCatch(penalty_has_kink(u$penalty), error = function(e) TRUE)
    if (isTRUE(kink)) out <- c(out, as.integer(u$index))
  }
  sort(unique(out))
}


#' One Likelihood Statistic for One Coefficient
#'
#' @description
#' Computes the Wald, likelihood-ratio, Rao score or Terrell gradient
#' statistic for the hypothesis that one coefficient equals `value`, with its
#' p-value.
#'
#' @details
#' All four read the penalized objective, which is the function that the fit
#' maximizes. Writing \eqn{\ell_p(\beta) = \ell(\beta) - \rho(\beta;\theta)}
#' for it, \eqn{\hat\beta} for the unrestricted estimates, \eqn{\tilde\beta}
#' for the estimates in which \eqn{\beta_j} is held at \eqn{b},
#' \eqn{U = \partial\ell_p} for its score and \eqn{K} for the penalized
#' information, the four are
#' \deqn{W = (\hat\beta_j - b)^2 / \widehat{\mathrm{Var}}(\hat\beta_j),
#' \qquad LR = 2\{\ell_p(\hat\beta) - \ell_p(\tilde\beta)\},}
#' \deqn{S = U_j(\tilde\beta)^2 \, [K(\tilde\beta)^{-1}]_{jj}, \qquad
#' T = U_j(\tilde\beta)\,(\hat\beta_j - b),}
#' each compared with a \eqn{\chi^2_1}. Where the model carries no penalty
#' \eqn{\ell_p} is the log-likelihood and these are the classical four.
#'
#' The score statistic has the general form \eqn{U'K^{-1}U}. It reduces to
#' the product above because every other component of \eqn{U} vanishes at the
#' restricted maximum, which also indicates that the restricted fit has
#' converged.
#'
#' # The penalized objective and the likelihood
#'
#' The likelihood ratio is taken on \eqn{\ell_p} and not on \eqn{\ell}.
#' \eqn{\tilde\beta} maximizes \eqn{\ell_p} under the restriction and not
#' \eqn{\ell}, so \eqn{\ell(\tilde\beta)} can exceed
#' \eqn{\ell(\hat\beta)} and the difference of log-likelihoods can be
#' negative. The other three statistics read the penalized objective as well:
#' \eqn{U} and \eqn{K} are its gradient and its curvature, and the Wald
#' variance is \eqn{(H+S)^{-1}} under `vcov_type = bayesian_vcov()`.
#'
#' For a shrunk coordinate this is the reading that [vcov.StatmodFit()] gives
#' under `vcov_type = bayesian_vcov()`: the penalty is a log-prior, so \eqn{\ell_p} is
#' a log posterior at the hyperparameters in force, and the interval by
#' inversion is the set that the posterior does not reject, conditional on
#' them. The interpretation that the Wald interval on such a row carries
#' therefore extends to the other three statistics. On a smooth the inverted
#' interval and the Wald one agree to second order, and the inversion adds
#' the asymmetry.
#'
#' # Differences between the statistics
#'
#' The Wald statistic reads only the unrestricted fit and needs no refit,
#' which is why the summary reports it by default. It depends on how the
#' parameter is written: a nonlinear reparametrization of the coefficient
#' changes it, whereas the likelihood ratio, a difference of objective values
#' at the same two points, is not changed. The Wald statistic also loses power
#' against a distant alternative (the Hauck-Donner effect), because the
#' curvature is read at \eqn{\hat\beta} rather than under the null.
#'
#' The gradient statistic needs the restricted fit and no matrix inversion:
#' it is the score contracted with the distance between the two fits.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()].
#' @param value The value under the null, a single number. The default is
#'   `0`.
#' @param test A [test_statistic] object: [wald_test()], [lr_test()],
#'   [score_test()] or [gradient_test()].
#' @param vcov_type A [vcov_method] object, the variance that the Wald statistic
#'   reads, as [vcov.StatmodFit()] takes it.
#' @param mode_error Whether to read how far above its mode the restricted
#'   fit stopped ([restricted_mode_error()]). The default is `FALSE`, because
#'   reading it costs one Hessian where the likelihood ratio and the gradient
#'   statistic need none. [summary.StatmodFit()], which calls this once per
#'   row, and [statmod_invert()], which calls it once per evaluation of the
#'   statistic, do not read it.
#' @param objective0 The penalized objective at the unrestricted joint mode,
#'   as [joint_objective()] returns it, for a caller that reads the
#'   likelihood ratio many times on one fit. `NULL` computes it here.
#'
#' @return A list with `test` (the [test_statistic] object), `statistic`, `df`,
#'   `p.value`, `converged` --
#'   `NA` for the Wald statistic and the restricted fit's flag for the other
#'   three -- and `mode_error`, `NA` unless it is requested.
#'
#' @seealso [statmod_restrict()], which the three restricted statistics read.
#'
#' @keywords internal
statmod_stat_at <- function(fit, param, coefname, value = 0,
                         test = wald_test(), vcov_type = bayesian_vcov(),
                         mode_error = FALSE, objective0 = NULL) {
  test <- as_test_statistic(test)
  vcov_type <- as_vcov_method(vcov_type)
  tname <- test@name
  spec <- fit@spec
  design <- statmod_design(spec)
  params <- spec@distrib@params
  nms <- design[[param]]$coef_names
  j <- match(coefname, nms)
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  at <- if (is.na(j)) NA_integer_ else offs[[match(param, params)]] + j
  bhat <- if (is.na(j)) NA_real_ else fit@coefficients[[param]][[j]]

  out <- function(s, conv, me = NA_real_) {
    list(test = test, statistic = s, df = 1L,
         p.value = stats::pchisq(s, 1L, lower.tail = FALSE),
         converged = conv, mode_error = me, df_denom = Inf)
  }

  if (!test_restricted(test)) {
    # the only one that needs no refit, and the only one a reparametrization
    # moves. Nothing was refitted, so there is no mode to be above
    V <- stats::vcov(fit, vcov_type = vcov_type)
    lab <- rownames(coef_labels(spec, design))[at]
    v <- if (is.na(at) || !lab %in% rownames(V)) NA_real_ else V[lab, lab]
    if (!isTRUE(is.finite(v)) || v <= 0) return(out(NA_real_, NA))
    res <- out((bhat - value)^2 / v, NA)
    # CR2: the t reference with Satterthwaite degrees of freedom
    dfs <- sandwich_wald_df(fit, vcov_type, lab)
    if (is.finite(dfs)) {
      res$p.value <- wald_pvalue(sqrt(res$statistic), dfs)
      res$df_denom <- dfs
    }
    return(res)
  }

  r <- statmod_restrict(fit, param, coefname, value)
  me <- if (isTRUE(mode_error)) r$mode_error() else NA_real_
  if (identical(tname, "lr")) {
    # 2 times the rise in the PENALIZED OBJECTIVE, which is the function both
    # fits maximize, so the difference is non-negative by construction. Read
    # on the log-likelihood alone it is not: the restricted point maximizes
    # the penalized objective, not the likelihood, so it can sit ABOVE the
    # unrestricted point there -- measured at -0.0379 and -0.1075 on a
    # smooth's own coordinate. With no penalty the two are the same number.
    # The unrestricted side is the JOINT mode, which is what the restricted
    # refit reaches: where the criterion estimated some coefficients -- a
    # dispersion by REML, the default -- the fitted point is not that mode,
    # and differencing against it subtracted 2 (l_ML - l_REML), 0.118 on a
    # linear regression at n = 40, negative at the estimate itself.
    o0 <- if (is.null(objective0)) joint_objective(fit) else objective0
    return(out(2 * (r$objective - o0), r$converged, me))
  }
  u <- r$score[[r$at]]
  if (identical(tname, "gradient")) {
    return(out(u * (bhat - value), r$converged, me))
  }
  K <- r$information()
  V <- tryCatch(solve_pd(K, "the penalized information at the restricted fit",
                         r$labels),
                error = function(e) NULL)
  if (is.null(V)) return(out(NA_real_, r$converged, me))
  out(u^2 * V[r$at, r$at], r$converged, me)
}


#' The Penalized Objective at the Unrestricted Joint Mode
#'
#' @description
#' Computes the minimum of the penalized objective over every coefficient at
#' the fitted hyperparameters, which is the unrestricted side of the
#' likelihood ratio that [statmod_stat_at()] computes.
#'
#' @details
#' The restricted fit of [statmod_restrict()] holds one coefficient and takes
#' every other one to the joint mode of the penalized likelihood. The
#' unrestricted side has to be read at the same kind of point, or the
#' statistic compares two fits whose other coefficients were chosen by
#' different rules. Where the outer criterion estimated some coefficients
#' itself (by default [reml()] and [ml()] estimate on the criterion the
#' unpenalized coefficients of every parameter but the first, see
#' [marginal_coords()]), the point that the fit returned is not the joint
#' mode, and `fit@objective` is not its minimum. A likelihood ratio
#' differenced against `fit@objective` then contains the difference between
#' the objective at the fitted point and at the joint mode.
#'
#' The full model is therefore refitted from the fitted coefficients, which
#' are already close to the joint mode, at the fitted hyperparameters, and no
#' coefficient is held. Where the criterion estimated nothing, the fit is at the
#' joint mode and `fit@objective` is returned unchanged.
#'
#' @param fit A [StatmodFit()].
#'
#' @return A single number, on the scale of `fit@objective`.
#'
#' @seealso [statmod_stat_at()], [statmod_restrict()].
#'
#' @keywords internal
joint_objective <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  gam <- tryCatch(marginal_coords(spec, design, fit@methods$outer),
                  error = function(e) marginal_coords(spec, design, NULL))
  if (!length(gam$where)) return(fit@objective)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$smooth
  cfg <- inner_settings(method)
  beta <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
  res <- statmod_alternate(spec, design, blocks, fit@hyper, method, beta,
                           cfg$expected, cfg$emethod, cfg$maxit, cfg$tol,
                           verbosity(0))
  res$value
}


#' A Test of One Coefficient Against One Value
#'
#' @description
#' Holds the result of testing that one coefficient of a fit equals one
#' value: the statistic, its distribution, the p-value and what was tested.
#'
#' @details
#' The four statistics are compared with the same \eqn{\chi^2_1}, so `df` is
#' always 1 and the object holds one kind of quantity whatever produced it.
#'
#' # Whether the refit reached its mode
#'
#' `mode_error` and `converged` both describe the restricted refit and both
#' are `NA` for the Wald statistic, which reads the unrestricted fit and
#' refits nothing. They measure different things, and only the first is
#' printed.
#'
#' `converged` is the inner optimizer's flag: whether a stopping rule fired,
#' a boolean about a threshold on a score whose size depends on the model.
#' `mode_error` is [restricted_mode_error()], the decrease of the penalized
#' objective that a Newton step from the point returned by the refit would
#' predict, in log-likelihood units, judged against [mode_error_limit()].
#' The two need not agree, because a stopping rule reads a change in the
#' objective whereas the mode error reads the distance to the mode. See
#' [restricted_mode_error()].
#'
#' A statistic computed from a refit that stopped short is the value at the
#' point where the refit stopped. `print()` reports this, reading
#' `mode_error`.
#'
#' The class declares its properties, as every object that this toolkit
#' returns does, whereas the `htest` returned by the tests of base R is a list
#' with a class attribute and no declared structure, so that it cannot be
#' validated.
#'
#' @param test The [test_statistic] object that was computed: [wald_test()],
#'   [lr_test()], [score_test()] or [gradient_test()]. Its `label` is the name
#'   under which `print()` shows the statistic.
#' @param statistic Its value, compared with a \eqn{\chi^2_1}.
#' @param df One, always.
#' @param p.value The upper tail of that \eqn{\chi^2_1}.
#' @param estimate The unrestricted estimate of the coefficient.
#' @param null_value The value against which the coefficient was tested.
#' @param param,coefficient Which coefficient of which equation.
#' @param df_denom The denominator degrees of freedom when the statistic is
#'   referred to an \eqn{F(1, \nu)} distribution, which is the case of a Wald
#'   statistic under the CR2 variance with a \eqn{t} reference; `Inf`
#'   otherwise, when the reference is the chi-squared distribution.
#' @param mode_error How far above its own mode the restricted refit
#'   stopped, in log-likelihood units; `NA` for the Wald statistic.
#' @param converged The inner optimizer's flag at the restricted refit;
#'   `NA` for the Wald statistic.
#'
#' @return An object of class `StatmodTest`.
#'
#' @seealso [statmod_test()], which builds one, and
#'   [restricted_mode_error()], which `mode_error` reports.
#'
#' @examples
#' dd <- data.frame(x = seq(-1, 1, length.out = 40))
#' dd$y <- rpois(40, exp(0.3 + dd$x))
#' fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
#' S7::S7_inherits(statmod_test(fit, param = "mu", coefname = "x", value = 1),
#'                 StatmodTest)
#'
#' @name StatmodTest-class
#' @aliases StatmodTest
#' @keywords internal
#' @export
StatmodTest <- S7::new_class("StatmodTest",
  properties = list(
    # which of the four statistics, as the object `statmod_test()` received;
    # print() reads its label
    test = test_statistic,
    statistic = S7::class_numeric,
    df = S7::class_integer,
    p.value = S7::class_numeric,
    # the unrestricted estimate and the value it was tested against, so the
    # object says what the null WAS rather than leaving it to the call
    estimate = S7::class_numeric,
    null_value = S7::class_numeric,
    param = S7::class_character,
    coefficient = S7::class_character,
    # how far above its own mode the restricted refit stopped, and the raw
    # flag beside it. The first is what print() reads; the second is kept
    # because it is what the optimizer reported and a reader may want it
    mode_error = S7::class_numeric,
    converged = S7::class_logical,
    # the denominator degrees of freedom of an F statistic: finite when the
    # Wald statistic is referred to a t distribution (CR2), infinite otherwise
    df_denom = S7::new_property(S7::class_numeric, default = Inf)
  )
)


#' Print a Test
#'
#' @description
#' Prints the statistic, its degrees of freedom and the p-value on one line,
#' under a title naming the test, in the layout that R uses for a hypothesis
#' test. A note follows when the restricted refit stopped above its mode, by
#' more than [mode_error_limit()], or when its mode error could not be read.
#'
#' @param x A [StatmodTest()].
#' @param digits Significant digits.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#'
#' @keywords internal
print.StatmodTest <- function(x, digits = 4L, ...) {
  cat("\n\t", x@test@label, " test on a statmod fit\n\n", sep = "")
  cat("coefficient:  ", x@param, ":", x@coefficient, "\n", sep = "")
  dd <- if (length(x@df_denom) == 1L) x@df_denom else Inf
  if (is.finite(dd)) {
    cat("F(", x@df, ", ", format(dd, digits = digits), ") = ",
        format(x@statistic, digits = digits), ",  p-value ",
        format.pval(x@p.value, digits = digits, eps = 1e-16), "\n", sep = "")
  } else cat("X-squared = ", format(x@statistic, digits = digits),
      ",  df = ", x@df,
      ",  p-value ", format.pval(x@p.value, digits = digits, eps = 1e-16),
      "\n", sep = "")
  cat("alternative hypothesis: true ", x@coefficient, " is not equal to ",
      format(x@null_value, digits = digits), "\n", sep = "")
  cat("estimate: ", format(x@estimate, digits = digits), "\n", sep = "")
  # a statistic read off a refit that stopped short is the value where it
  # stopped, and saying so is the difference between a number and a result.
  # What is read is the MODE ERROR and not the optimizer's flag. Measured,
  # the flag rejects exactly the best-located points -- where the refit lands
  # on its mode in one step the objective does not move and the stall guard
  # fires -- so a warning printed off it fires loudest where the refit is
  # best. See restricted_mode_error().
  # Wald refits nothing, so there is no mode for it to be above, and
  # `converged` is NA exactly there
  refitted <- length(x@converged) == 1L && !is.na(x@converged)
  me <- if (length(x@mode_error) == 1L) x@mode_error else NA_real_
  if (refitted && is.na(me)) {
    cat("\nthe mode error of the restricted fit could not be read,\n",
        "  so the statistic is the value at wherever the refit stopped\n",
        sep = "")
  } else if (refitted && me > mode_error_limit()) {
    cat("\nthe restricted fit stopped above its own mode",
        sprintf(" (by %s, against %s)", format(me, digits = digits),
                format(mode_error_limit(), digits = digits)),
        ",\n  so the statistic is the value at wherever it stopped\n",
        sep = "")
  }
  invisible(x)
}
S7::method(print, StatmodTest) <- print.StatmodTest


#' Test One Coefficient Against One Value
#'
#' @description
#' Tests that one coefficient of a fit equals one value, by Wald's
#' statistic, the likelihood ratio, Rao's score or Terrell's gradient.
#'
#' @details
#' `summary(fit, test =)` reports one of the four statistics for the rows of
#' the coefficient tables that can be held (not a row under a kinked penalty,
#' not an aliased one, and none for a model with a structural term), always
#' against zero. This function tests one coefficient against any value, for
#' questions that a summary cannot answer, such as whether a slope is one,
#' whether an elasticity is unity, or whether a coefficient matches a value
#' fixed outside the data.
#'
#' The statistics, the quantities that they read and their cost are described
#' at [statmod_stat_at()]. The last three are built on the restricted fit of
#' [statmod_restrict()], and inverting one gives the interval of
#' [confint.StatmodFit()].
#'
#' Reading one row costs one restricted refit, whereas `summary(test =)`
#' costs one per row. The function also reads the mode error of the refit,
#' which costs one Hessian of the restricted fit; the Rao score statistic
#' reads that matrix anyway, so that it adds no cost there. The summary and
#' the interval do not read the mode error.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()]:
#'   the distribution parameter and the coefficient's name in its equation.
#' @param value The value under the null, a single number. The default is
#'   `0`, which is the value that a summary reports.
#' @param test A [test_statistic] object: [wald_test()] (the default),
#'   [lr_test()], [score_test()] or [gradient_test()]. The first reads the
#'   unrestricted fit only; the other three cost one restricted refit. The
#'   four are compared at [test_statistic]. A character string is rejected with
#'   an error.
#' @param vcov_type A [vcov_method] object, the variance that the Wald statistic
#'   reads, passed to [vcov.StatmodFit()]: [bayesian_vcov()] (the default),
#'   [frequentist_vcov()] or [unconditional_vcov()]. The conventions are
#'   compared at [vcov_method].
#'
#' @return A [StatmodTest()].
#'
#' @seealso [test_statistic] for the four statistics, [vcov_method] for the
#'   variance conventions, [summary.StatmodFit()] for the rows that can be
#'   held, against zero, [confint.StatmodFit()] for the interval that a test
#'   inverts.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200))
#' dd$y <- rpois(200, exp(0.3 + dd$x))
#' fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
#' ## the summary tests whether the slope is zero; this tests whether it is one
#' statmod_test(fit, param = "mu", coefname = "x", value = 1, test = lr_test())
#' statmod_test(fit, param = "mu", coefname = "x", value = 0, test = lr_test())
#'
#' @export
statmod_test <- function(fit, param, coefname, value = 0, test = wald_test(),
                         vcov_type = bayesian_vcov()) {
  test <- as_test_statistic(test)
  vcov_type <- as_vcov_method(vcov_type)
  # the mode error is asked for HERE and nowhere else: this is the one
  # surface a reader reads one row of, so one extra Hessian is affordable,
  # where the two loops that call the same helper in quantity are not
  s <- statmod_stat_at(fit, param, coefname, value, test, vcov_type,
                       mode_error = TRUE)
  design <- statmod_design(fit@spec)
  nms <- design[[param]]$coef_names
  j <- match(coefname, nms)
  est <- if (is.na(j)) NA_real_ else fit@coefficients[[param]][[j]]
  StatmodTest(test = s$test, statistic = s$statistic, df = s$df,
              p.value = s$p.value, estimate = est,
              null_value = as.numeric(value), param = param,
              coefficient = coefname, mode_error = s$mode_error,
              converged = s$converged, df_denom = s$df_denom)
}


#' A Confidence Interval by Inverting a Test
#'
#' @description
#' Returns the set of values that a test does not reject at `level`, found by
#' locating the two points where the statistic crosses its critical value.
#'
#' @details
#' An interval by inversion is \eqn{\{b : T(b) \le \chi^2_{1,\alpha}\}}, and
#' it is not symmetric about the estimate unless the statistic is a parabola
#' in \eqn{b}. The Wald statistic is a parabola by construction, because the
#' curvature is read once at \eqn{\hat\beta}, and the other three statistics
#' are not. The asymmetry is the reason to invert one of the other
#' three: the interval follows the shape of the likelihood and not a quadratic
#' fitted at its maximum.
#'
#' The search starts at the half-width of the Wald interval and widens the
#' distance from the estimate by a factor of 1.6, at most `maxit` times on
#' each side, until the statistic exceeds its critical value. It then locates
#' the crossing with [stats::uniroot()] (Brent's method) at tolerance 1e-7.
#' Each evaluation is one restricted refit, warm-started at the unrestricted
#' estimates and run at the fitted hyperparameters, so the interval requires
#' no outer search.
#'
#' The position of a break-point has a range: a value at or beyond the
#' smallest or the largest observation of the variable leaves the term with no
#' effect, so that every such value gives the same refit. [break_edges()]
#' returns that range, and the search never leaves it. On a side where the
#' statistic is still under its critical value when the search reaches the end
#' of the range, the data do not bound the position, and the limit is
#' `-Inf` or `Inf`.
#'
#' A term with two or more break-points restricts the search further. The
#' search for a break-point stays between the fitted positions of the
#' break-points next to it ([neighbor_positions()]), because a refit that
#' holds one break-point is a fit with a break-point at that value and the
#' other free, and beyond a neighbor the interval would become that of the
#' neighbor. The statistic is not monotone next to a neighbor, so that side is
#' evaluated at the points of the widening and at points whose distances from
#' the neighbor grow by a factor of 1.3, and on a side where the statistic is
#' still under its critical value at the neighbor's position, the limit is that
#' position.
#'
#' @param fit A [StatmodFit()].
#' @param param,coefname Which coefficient, as in [statmod_restrict()].
#' @param level The confidence level, a single number in (0, 1).
#' @param test The [test_statistic] to invert: [lr_test()] (the default),
#'   [score_test()], [gradient_test()] or [wald_test()].
#' @param vcov_type The [vcov_method] that sets the starting bracket, and that the
#'   Wald statistic reads.
#' @param maxit How many widenings to try on each side before the search
#'   stops.
#'
#' @return A numeric vector of two, `lower` and `upper`. A side of the
#'   position of a break-point that the data do not bound is `-Inf` or `Inf`,
#'   and a side that stays under the critical value as far as the fitted
#'   position of the next break-point is that position.
#'   An end that the search could not bracket within `maxit` widenings, or at
#'   which the statistic could not be evaluated, is `NA`; both ends are `NA`
#'   when the standard error of the coefficient is not finite and positive.
#'
#' @seealso [statmod_stat_at()], the statistic that is inverted.
#'
#' @keywords internal
statmod_invert <- function(fit, param, coefname, level = 0.95,
                           test = lr_test(), vcov_type = bayesian_vcov(),
                           maxit = 40L) {
  test <- as_test_statistic(test)
  vcov_type <- as_vcov_method(vcov_type)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  q <- stats::qchisq(level, 1L)
  spec <- fit@spec
  design <- statmod_design(spec)
  params <- spec@distrib@params
  j <- match(coefname, design[[param]]$coef_names)
  if (is.na(j)) {
    stop(sprintf("'%s' names no coefficient of the '%s' equation.",
                 coefname, param), call. = FALSE)
  }
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  at <- offs[[match(param, params)]] + j
  bhat <- fit@coefficients[[param]][[j]]
  V <- stats::vcov(fit, vcov_type = vcov_type, readable = FALSE)
  lab <- rownames(coef_labels(spec, design))[at]
  se <- if (lab %in% rownames(V)) sqrt(V[lab, lab]) else NA_real_
  if (!isTRUE(is.finite(se)) || se <= 0) return(c(lower = NA_real_,
                                                  upper = NA_real_))

  o0 <- if (identical(test@name, "lr")) joint_objective(fit) else NULL
  gap <- function(b) {
    s <- statmod_stat_at(fit, param, coefname, b, test, vcov_type,
                         objective0 = o0)$statistic
    if (!isTRUE(is.finite(s))) return(NA_real_)
    s - q
  }
  edges <- break_edges(spec, design, param, j)
  edges_range <- edges
  nb <- neighbor_positions(fit, spec, design, param, j)
  stopped <- c(FALSE, FALSE)
  if (is.finite(nb[[1L]]) && nb[[1L]] > edges[[1L]]) {
    edges[[1L]] <- nb[[1L]]
    stopped[[1L]] <- TRUE
  }
  if (is.finite(nb[[2L]]) && nb[[2L]] < edges[[2L]]) {
    edges[[2L]] <- nb[[2L]]
    stopped[[2L]] <- TRUE
  }
  # THE STATISTIC NEXT TO A NEIGHBOR IS NOT MONOTONE: a break-point held a
  # few observations from another one fits worse than either of the two fits
  # that hold one at the position of the other, so an interval above the
  # critical value can lie between the estimate and the neighbor, and a
  # bracket that widens by a fixed factor can jump over it. A side that ends
  # at a neighbor is evaluated at the points of that widening and at points
  # that grow their distance from the neighbor by a factor of 1.3, from the
  # estimate outward, and the first crossing is located between the two points
  # around it.
  bounded <- function(dir, end) {
    span <- abs(end - bhat)
    d0 <- 0.005 * diff(edges_range)
    far <- end - dir * d0 * 1.3^(0:60)
    far <- far[dir * (far - bhat) > 0]
    wide <- bhat + dir * sqrt(q) * se * 1.6^(0:60)
    wide <- wide[dir * (wide - end) < 0]
    pts <- sort(unique(c(wide, far, end)), decreasing = dir < 0)
    pts <- pts[abs(pts - bhat) <= span]
    a <- bhat
    for (b in pts) {
      v <- gap(b)
      if (is.na(v)) return(NA_real_)
      if (v > 0) {
        return(tryCatch(stats::uniroot(gap, interval = sort(c(a, b)),
                                       tol = 1e-7)$root,
                        error = function(e) NA_real_))
      }
      a <- b
    }
    end
  }
  side <- function(dir) {
    end <- if (dir < 0) edges[[1L]] else edges[[2L]]
    if (stopped[[if (dir < 0) 1L else 2L]]) return(bounded(dir, end))
    step <- sqrt(q) * se
    a <- bhat
    for (k in seq_len(as.integer(maxit))) {
      b <- bhat + dir * step
      # the refit at the end of the range is the last one: beyond it the
      # statistic does not change
      last <- dir * (b - end) >= 0
      if (last) b <- end
      v <- gap(b)
      if (is.na(v)) return(NA_real_)
      if (v > 0) {
        r <- tryCatch(stats::uniroot(gap, interval = sort(c(a, b)),
                                     tol = 1e-7)$root,
                      error = function(e) NA_real_)
        return(r)
      }
      if (last) return(dir * Inf)
      a <- b
      step <- step * 1.6
    }
    NA_real_
  }
  c(lower = side(-1), upper = side(1))
}


#' The Range of a Break-Point Position
#'
#' @description
#' Returns the smallest and the largest value that the position of a
#' break-point can be held at before the term stops depending on it: the range
#' of the variable the break-point is on.
#'
#' @details
#' A position at or beyond the largest observation makes every hinge
#' \eqn{(x - \psi)_+} and every step \eqn{1(x > \psi)} of the term zero, and a
#' position at or below the smallest makes the hinge a line that the term's
#' linear part absorbs. Either way the restricted refit is the same for every
#' such value, so that an interval by inversion can stop at the end of the
#' range. A coefficient that is not the position of a break-point of a
#' [modelterms7::seg()] term, or a term whose positions are developed over
#' covariates, has no range.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param param The distribution parameter of the coefficient.
#' @param j The position of the coefficient in the parameter's equation.
#'
#' @return A numeric vector of two, `-Inf` and `Inf` where there is no range.
#'
#' @seealso [statmod_invert()], which stops its search there.
#'
#' @keywords internal
break_edges <- function(spec, design, param, j) {
  none <- c(-Inf, Inf)
  blocks <- design[[param]]$blocks
  nm <- names(blocks)[vapply(blocks, function(ix) j %in% ix, logical(1))]
  if (length(nm) != 1L) return(none)
  term <- spec@terms[[param]][[nm]]
  if (!S7::S7_inherits(term, modelterms7::SegTerm) ||
      !identical(term@kind, "seg") || length(term@subformulas)) {
    return(none)
  }
  if (!grepl("psi[0-9]+$", design[[param]]$coef_names[[j]])) return(none)
  xv <- term@blueprint$xv
  if (!is.numeric(xv) || !length(xv) || !all(is.finite(xv))) return(none)
  range(xv)
}


#' Fitted Positions of the Break-Points Next to One
#'
#' @description
#' Returns the fitted positions of the break-points of the same term that lie
#' immediately below and above the break-point whose position is the
#' coefficient `j`, which bound the search for an interval by inversion.
#'
#' @details
#' Only a position of a [modelterms7::seg()] term with two or more
#' break-points and positions that carry no development has neighbors. A
#' neighbor is a position strictly below or above the estimate, read at the
#' confinement of the term.
#'
#' @param fit A [StatmodFit()].
#' @param spec,design The specification and its design.
#' @param param The distribution parameter of the coefficient.
#' @param j The position of the coefficient in the parameter's equation.
#'
#' @return A numeric vector of two, the nearest position below and the
#'   nearest above, `NA` where there is none.
#'
#' @seealso [statmod_invert()], which stops its search there, and
#'   [break_edges()].
#'
#' @keywords internal
neighbor_positions <- function(fit, spec, design, param, j) {
  none <- c(NA_real_, NA_real_)
  blocks <- design[[param]]$blocks
  nm <- names(blocks)[vapply(blocks, function(ix) j %in% ix, logical(1))]
  if (length(nm) != 1L) return(none)
  term <- spec@terms[[param]][[nm]]
  if (!S7::S7_inherits(term, modelterms7::SegTerm) ||
      !identical(term@kind, "seg") || length(term@subformulas) ||
      term@npsi < 2L) {
    return(none)
  }
  bp <- term@blueprint
  cols <- blocks[[nm]]
  at <- cols[vapply(paste0("psi", seq_len(term@npsi)),
                    function(q) bp$index[[q]][1L], integer(1))]
  if (!j %in% at) return(none)
  pos <- pmin(pmax(fit@coefficients[[param]][at], bp$lim[1L]), bp$lim[2L])
  b <- pos[match(j, at)]
  others <- pos[at != j]
  c(if (any(others < b)) max(others[others < b]) else NA_real_,
    if (any(others > b)) min(others[others > b]) else NA_real_)
}


#' Which Coefficients a Restricted Fit Can Hold
#'
#' @description
#' Returns the stacked positions of the coefficients for which the three
#' restricted statistics are defined: every coordinate less the ones under a
#' kinked penalty, less the aliased ones, and none at all where the model
#' carries a structural term.
#'
#' @details
#' It applies to the whole vector the enumeration that [statmod_restrict()]
#' applies to one coordinate at a time, signaling an error for a coordinate
#' that cannot be held. A caller reporting a table needs it in this form,
#' since each restricted statistic is a refit and an interval by inversion
#' requires many, so that a row that cannot be tested has to be recognized
#' before its refit is computed and not after.
#'
#' A coordinate under a penalty that is twice differentiable is included:
#' a smooth's linear column and its rotated coordinates, a ridge, a random
#' effect. The restricted fit holds it exactly as it holds an unpenalized
#' one, because [fit_smooth()] enforces the hold in the solve. Only a kinked
#' penalty is excluded, for the reasons that [statmod_restrict()] gives.
#'
#' @param fit A [StatmodFit()].
#' @param spec,design The specification and its design, where the caller has
#'   them already.
#'
#' @return An integer vector of stacked positions, possibly empty, named by
#'   the labels that [coef_labels()] gives them.
#'
#' @seealso [statmod_restrict()], which rejects the same coordinates one at a
#'   time, and [kinked_coords()].
#'
#' @keywords internal
testable_coords <- function(fit, spec = fit@spec,
                            design = statmod_design(spec)) {
  none <- stats::setNames(integer(0), character(0))
  if (length(attr(design, "structural"))) return(none)
  lab <- rownames(coef_labels(spec, design))
  free <- setdiff(seq_along(lab), kinked_coords(spec, design))
  if (!length(free)) return(none)
  keep <- free[!lab[free] %in% fit@aliased]
  if (!length(keep)) return(none)
  stats::setNames(as.integer(keep), lab[keep])
}
