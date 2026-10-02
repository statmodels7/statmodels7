#' @include outer.R
NULL

#' The Distribution Parameters a Marginal Criterion Estimates
#'
#' @description
#' Resolves the `marginal` argument of [reml()] and [ml()] against the
#' family: the parameters whose unpenalized coefficients are estimated by
#' maximizing the marginal criterion.
#'
#' @details
#' `NULL` names every parameter except the position, which is the family's
#' first parameter. In every shipped univariate family the first parameter is
#' the location or, where the family has none, the scale (the `mu` of
#' [distributions7::weibull1_distrib()], the `sigma` of
#' [distributions7::gpd_distrib()]).
#'
#' A multivariate family is not covered. An explicit request is refused and
#' the default names nothing, which is the convention of the criterion
#' without the argument.
#'
#' @param method An [OuterMethod()], or `NULL`.
#' @param distrib The family.
#'
#' @return A list with `named`, a character vector of parameter names in the
#'   family's order, and `explicit`, a single logical saying whether the
#'   caller named them (`TRUE`) or the default did (`FALSE`).
#'
#' @seealso [marginal_coords()], which turns the names into coefficients.
#'
#' @keywords internal
marginal_params <- function(method, distrib) {
  none <- list(named = character(0), explicit = FALSE)
  if (is.null(method) || !method@kind %in% c("reml", "ml")) return(none)
  params <- distrib@params
  m <- method@marginal
  if (is.null(m)) {
    if (S7::S7_inherits(distrib, distributions7::multivariate_distrib)) {
      return(none)
    }
    return(list(named = params[-1L], explicit = FALSE))
  }
  if (identical(m, "none")) return(none)
  named <- if (identical(m, "all")) params else m
  bad <- setdiff(named, params)
  if (length(bad)) {
    stop(sprintf(paste0("'marginal' names %s, which %s not a parameter of",
                        " the family.\n  Its parameters are: %s."),
                 paste0("'", bad, "'", collapse = ", "),
                 if (length(bad) == 1L) "is" else "are",
                 paste(params, collapse = ", ")), call. = FALSE)
  }
  if (S7::S7_inherits(distrib, distributions7::multivariate_distrib)) {
    stop(paste0("'marginal' does not cover a multivariate family yet.",
                " Use\n  marginal = \"none\"."), call. = FALSE)
  }
  list(named = params[params %in% named], explicit = TRUE)
}


#' The Coefficients a Marginal Criterion Estimates
#'
#' @description
#' The positions, in the stacked coefficient vector, of the coordinates that
#' [reml()] and [ml()] estimate by maximizing the criterion: in the equation
#' of every parameter [marginal_params()] names, the coefficients no penalty
#' shrinks.
#'
#' @details
#' A coefficient qualifies when no penalty covers it, or when it lies in the
#' null space of the quadratic penalty covering it and that null space is
#' spanned by coordinates. The free columns of [modelterms7::s()] are the
#' second case: the Demmler-Reinsch penalty is \eqn{\mathrm{diag}(0, 1,
#' \ldots, 1)}, so its null space is the linear column.
#'
#' A coefficient a penalty with a kink covers (lasso, SCAD, MCP) is not
#' one of them: it stays at the joint mode, and the criterion leaves it out
#' of its determinant ([laplace_pinned()]). The unpenalized coefficients
#' beside it are estimated on the criterion as anywhere else, which corrects
#' a dispersion's intercept for the columns of the mean it is fitted with.
#'
#' Four configurations are not covered yet, and for each one a parameter
#' the caller named is refused while a parameter the default named keeps the
#' convention of the criterion without the argument:
#' \itemize{
#'   \item a model carrying a structural term ([modelterms7::gas()],
#'     [modelterms7::regime()]), whose joint fit does not hold a coefficient;
#'   \item an equation carrying a block that moves with its coefficients
#'     ([modelterms7::nl()], [modelterms7::seg()] and the other break-point
#'     terms);
#'   \item a model carrying a block that is a working linearization rather
#'     than a Jacobian, which is the case of a sharp [modelterms7::jump()]
#'     and [modelterms7::jseg()], in any equation: the determinant over its
#'     columns reads no curvature;
#'   \item an equation carrying a penalty whose null space is not spanned by
#'     coordinates, which is the case of [modelterms7::te()].
#' }
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param method An [OuterMethod()], or `NULL`.
#'
#' @return A list with `where` (integer positions in the stacked vector),
#'   `param` and `name` (the equation and the coefficient name of each
#'   position), and `skipped`, a named character vector giving, for every
#'   parameter the default named and this function left out, the reason.
#'   Every field is empty where nothing is estimated this way.
#'
#' @seealso [marginal_params()], [statmod_hold()], which holds these
#'   coefficients in the inner fit.
#'
#' @keywords internal
marginal_coords <- function(spec, design, method) {
  empty <- list(where = integer(0), param = character(0),
                name = character(0), skipped = character(0))
  mp <- marginal_params(method, spec@distrib)
  if (!length(mp$named)) return(empty)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  units <- statmod_penalized(spec, design)
  # BESIDE A PENALTY WITH A KINK THE DEFAULT NAMES NOTHING. Its coordinates
  # stay at the joint mode, and a coefficient moved on the criterion moves
  # that mode with it: a lasso's lambda is chosen by a path that scores the
  # model with every other coefficient at the joint mode, and reading them on
  # the criterion afterwards returned another model -- on UScrime a SCAD
  # scored with 7 covariates came back with 6 and not converged, the
  # criterion jumping where the active set changes. A parameter the caller
  # names is still estimated on the criterion.
  if (!mp$explicit) {
    kinked <- any(vapply(units, function(u)
      isTRUE(tryCatch(penalty_has_kink(u$penalty), error = function(e) FALSE)),
      logical(1)))
    if (kinked) {
      empty$skipped <- stats::setNames(
        rep("the model carries a penalty with a kink", length(mp$named)),
        mp$named)
      return(empty)
    }
  }
  structural <- length(attr(design, "structural")) > 0L
  refresh <- vapply(attr(design, "refresh"), function(r) r$param,
                    character(1))
  # a block that is a working linearization rather than a Jacobian has no
  # curvature for a determinant to read: at one break-point the criterion
  # read -132.9 or -109.8 according to where the iteration began
  frozen <- any(vapply(attr(design, "refresh"), function(r) isTRUE(r$frozen),
                       logical(1)))

  out <- empty
  for (p in mp$named) {
    a <- match(p, params)
    if (is.null(design[[p]]) || npar[a] == 0L) next
    pos <- offs[a] + seq_len(npar[a])
    reason <- NULL
    if (structural) {
      reason <- "the model carries a structural term"
    } else if (frozen) {
      reason <- paste0("the model carries a block that is a working",
                       " linearization rather than a Jacobian")
    } else if (p %in% refresh) {
      reason <- "its equation carries a block that moves with its coefficients"
    }
    free <- rep(TRUE, length(pos))
    if (is.null(reason)) {
      for (u in units) {
        if (isTRUE(u$structural)) next
        ix <- if (isTRUE(u$mixed)) u$beta_index else u$index
        if (is.null(ix)) next
        inside <- ix %in% pos
        if (!any(inside)) next
        pen <- u$penalty
        kink <- tryCatch(penalty_has_kink(pen), error = function(e) TRUE)
        if (isTRUE(kink)) {
          free[match(ix[inside], pos)] <- FALSE
          next
        }
        nul <- if (isTRUE(penalties7::is_proper(pen))) integer(0) else
          null_coordinates(pen, length(ix))
        if (is.null(nul)) {
          reason <- paste0("its equation carries a penalty whose null space",
                           " is not spanned by coordinates")
          break
        }
        covered <- setdiff(ix[inside], ix[nul])
        free[match(covered, pos)] <- FALSE
      }
    }
    if (!is.null(reason)) {
      if (mp$explicit) {
        stop(sprintf(paste0("'marginal' names '%s', and %s. That case is not",
                            " covered yet;\n  leave '%s' out of 'marginal'."),
                     p, reason, p), call. = FALSE)
      }
      out$skipped[[p]] <- reason
      next
    }
    w <- pos[free]
    out$where <- c(out$where, w)
    out$param <- c(out$param, rep(p, length(w)))
    out$name <- c(out$name, design[[p]]$coef_names[w - offs[a]])
  }
  out
}


#' The Null Coordinates of a Penalty
#'
#' @description
#' The positions, within a penalty's block, that span its null space, or
#' `NULL` where the null space is not spanned by coordinates.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param k The block's width.
#'
#' @return An integer vector of positions within the block, empty for a
#'   penalty of full rank, or `NULL`.
#'
#' @keywords internal
null_coordinates <- function(pen, k) {
  N <- tryCatch(as.matrix(penalties7::penalty_null_basis(pen)),
                error = function(e) NULL)
  if (is.null(N)) return(NULL)
  if (!ncol(N)) return(integer(0))
  sup <- apply(N, 2L, function(v) {
    s <- which(abs(v) > 1e-8 * max(abs(v)))
    if (length(s) == 1L) s else NA_integer_
  })
  if (anyNA(sup) || anyDuplicated(sup)) return(NULL)
  sort(as.integer(sup))
}


#' The Coordinates Left Out of the Marginal Determinant
#'
#' @description
#' The positions [pin_boundary()] pins besides a boundary: the coefficients
#' the criterion estimates, which the specification holds, and the
#' coordinates the model does not identify, which [outer_fit()] records on
#' the design at its first fit.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#'
#' @return An integer vector of stacked positions, possibly empty.
#'
#' @keywords internal
pinned_coords <- function(spec, design) {
  sort(unique(c(held_stack(spec, design),
                as.integer(attr(design, "pinned")))))
}


#' The Coordinates the Marginal Determinant Leaves Out
#'
#' @description
#' [pinned_coords()] and the coordinates a penalty with a kink covers.
#'
#' @details
#' A kinked penalty has no curvature at the coordinates it sets to zero, so
#' a Laplace expansion around them is not defined, and the coordinates it
#' leaves away from zero are chosen by the same selection. The criterion
#' holds all of them at the joint mode and integrates the rest, which is the
#' restricted likelihood of Verbyla (1993) with those coefficients treated
#' as known. Measured on `y ~ x1 + ... + x20 | sigma ~ lasso(~ z1 + ... +
#' z10)` at 200 observations over 15 samples, the dispersion's intercept is
#' off by +0.003 on average (root mean square 0.048) against +0.013 (0.067)
#' when the kinked coordinates are integrated and -0.057 (0.073) when the
#' intercept is read at the joint mode, with the same count of slopes
#' wrongly selected (1.9 against 1.2 of seven) and every fit converged.
#'
#' The mode's own movement is not affected: it is read over every
#' coordinate, see [ctx_penalized()].
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#'
#' @return An integer vector of stacked positions, possibly empty.
#'
#' @references Verbyla, A. P. (1993). Modelling variance heterogeneity:
#'   residual maximum likelihood and diagnostics. \emph{Journal of the Royal
#'   Statistical Society B}, 55, 493--508.
#'
#' @keywords internal
laplace_pinned <- function(spec, design) {
  sort(unique(c(pinned_coords(spec, design), kinked_coords(spec, design))))
}


#' The Coordinates a Kinked Penalty Holds at Zero
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param coef The coefficients, a named list.
#'
#' @return An integer vector of stacked positions, possibly empty.
#'
#' @keywords internal
zero_kinked <- function(spec, design, coef) {
  kc <- kinked_coords(spec, design)
  if (!length(kc)) return(integer(0))
  b <- unlist(coef[spec@distrib@params], use.names = FALSE)
  kc[b[kc] == 0]
}


#' A Score Over the Coordinates a Kink Leaves Free
#'
#' @description
#' Zeroes the entries of a score at the coordinates [zero_kinked()] returns.
#'
#' @details
#' At a coordinate a kinked penalty holds at zero the objective has no
#' derivative, and the entry the smooth part reports is the log-likelihood's
#' score, which the kink's subgradient interval contains rather than
#' cancels. Read as a residual it inflated the mode error of a fit sitting at
#' its mode: on a lasso over a dispersion's ten covariates at 1000
#' observations, seven of them at zero, it read 2.47 log-likelihood units
#' where the free coordinates' own reading is 4.5e-08, and a marginal search
#' whose every point was refused a resolution for that reason ran out of
#' backtracks and reported failure. Whether such a coordinate should leave
#' zero is the kinked block's own condition, read by
#' [alternation_readings()].
#'
#' @param score A score over the stacked coefficients.
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param coef The coefficients, a named list.
#'
#' @return `score`, with those entries set to zero.
#'
#' @keywords internal
free_of_kinks <- function(score, spec, design, coef) {
  zk <- zero_kinked(spec, design, coef)
  if (length(zk)) score[zk] <- 0
  score
}


#' A Specification With Coefficients Held
#'
#' @description
#' Writes values for the coefficients [marginal_coords()] returned into the
#' specification's `held_coef`, which the inner fit enforces.
#'
#' @param spec A [StatmodSpec()].
#' @param coords The result of [marginal_coords()].
#' @param values The values, one per position in `coords$where`.
#'
#' @return `spec` with `held_coef` set. The holds already present are kept
#'   for coefficients `coords` does not name.
#'
#' @keywords internal
statmod_hold <- function(spec, coords, values) {
  if (!length(coords$where)) return(spec)
  hc <- spec@held_coef
  for (p in unique(coords$param)) {
    j <- coords$param == p
    v <- stats::setNames(as.numeric(values[j]), coords$name[j])
    old <- hc[[p]]
    if (length(old)) old <- old[setdiff(names(old), names(v))]
    hc[[p]] <- c(old, v)
  }
  S7::set_props(spec, held_coef = hc)
}


#' The Stacked Positions a Specification Holds
#'
#' @description
#' Translates `spec@held_coef` into positions in the stacked coefficient
#' vector, from the design's own coefficient names.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#'
#' @return An integer vector of positions, empty where nothing is held.
#'
#' @keywords internal
held_stack <- function(spec, design) {
  hc <- spec@held_coef
  th <- term_held_stack(spec, design)
  if (!length(hc)) return(th)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  out <- integer(0)
  for (p in names(hc)) {
    a <- match(p, params)
    if (is.na(a) || is.null(design[[p]])) next
    j <- match(names(hc[[p]]), design[[p]]$coef_names)
    out <- c(out, offs[a] + j[!is.na(j)])
  }
  sort(unique(c(out, th)))
}


#' Name Stacked Positions
#'
#' @description
#' The equation and the coefficient name of each position in the stacked
#' coefficient vector, in the shape [marginal_coords()] returns.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param pos Integer positions in the stacked vector.
#'
#' @return A list with `where`, `param` and `name`.
#'
#' @keywords internal
stack_coords <- function(spec, design, pos) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  a <- findInterval(pos, offs + 1L)
  list(where = as.integer(pos), param = params[a],
       name = vapply(seq_along(pos), function(i)
         design[[params[a[i]]]]$coef_names[pos[i] - offs[a[i]]], character(1)))
}


#' The Variance Where the Criterion Estimates Coefficients
#'
#' @description
#' The variance of the coefficients of a fit whose marginal criterion
#' estimated some of them ([marginal_coords()]), in the block form of the
#' penalized information with the curvature of the criterion in place of the
#' Schur complement of those coefficients.
#'
#' @details
#' Write \eqn{\gamma} for the estimated coefficients, \eqn{u} for the rest and
#' \eqn{K} for the penalized information. The coefficients \eqn{u} are the
#' mode of the penalized likelihood at \eqn{\gamma}, and
#' \eqn{b_j = e_j - K_{uu}^{-1}K_{u\gamma}e_j} is how that mode moves with
#' \eqn{\gamma_j}. With \eqn{v} the outer vector,
#' \deqn{\mathrm{Var}(\hat\beta) = K_{uu}^{-1} + B\,V_v B^\top,\qquad
#'   V_v = \big[-\nabla^2_{vv}\ell_M\big]^{-1},}
#' where \eqn{K_{uu}^{-1}} is padded with zeros on \eqn{\gamma} and \eqn{B}
#' holds one direction per coordinate of \eqn{v}. For `type = "bayesian"`
#' \eqn{v = \gamma} and the hyperparameters are held; for
#' `type = "unconditional"` \eqn{v = (\eta, \gamma)}, the hyperparameters'
#' directions being \eqn{-K_{uu}^{-1}\partial^2\rho/\partial u\,\partial\eta}.
#'
#' Replacing \eqn{-\nabla^2_{\gamma\gamma}\ell_M} by the Schur complement
#' \eqn{K_{\gamma\gamma} - K_{\gamma u}K_{uu}^{-1}K_{u\gamma}} gives back
#' \eqn{K^{-1}}, which is the convention of a fit that reads every
#' coefficient at the joint mode.
#'
#' @param object A [StatmodFit()].
#' @param design Its design.
#' @param A The penalized information over the kept coordinates.
#' @param keep The kept coordinates, a logical vector over the stacked
#'   coefficients.
#' @param type `"bayesian"` or `"unconditional"`.
#'
#' @return The variance over the kept coordinates, or `NULL` where the fit
#'   estimated no coefficient on its criterion, or where a piece cannot be
#'   read.
#'
#' @keywords internal
marginal_vcov <- function(object, design, A, keep, type) {
  method <- object@methods$outer
  if (is.null(method) || !method@kind %in% c("ml", "reml")) return(NULL)
  if (length(attr(design, "structural"))) return(NULL)
  mc <- fit_marginal_context(object@spec, design, object@coefficients,
                             method, object@methods$pinned)
  if (!length(mc$gam$where)) return(NULL)
  mb <- marginal_blocks(mc, object@coefficients, object@hyper, method, A,
                        keep, full = !identical(type, "bayesian"))
  if (is.null(mb)) return(NULL)
  if (identical(type, "bayesian") || !mb$nh) {
    sel <- mb$nh + seq_len(mb$ng)
    Vv <- hyper_variance(-mb$Hv[sel, sel, drop = FALSE])
    B <- mb$Bg
  } else {
    Vv <- hyper_variance(-mb$Hv)
    B <- cbind(mb$Be, mb$Bg)
  }
  if (is.null(Vv)) return(NULL)
  Vv <- zero_unread(Vv)
  V <- matrix(0, mb$p, mb$p)
  V[mb$ui, mb$ui] <- as_dense(mb$Pu)
  V <- V + B %*% Vv %*% t(B)
  (V + t(V)) / 2
}


#' The Pieces of the Block Variance Over the Estimated Coefficients
#'
#' @description
#' The quantities [marginal_vcov()] and [marginal_edf_correction()] share:
#' the inverse of the penalized information over the integrated
#' coefficients, the movement of those coefficients with the estimated ones
#' and with the hyperparameters, and the Hessian of the criterion over
#' \eqn{(\eta, \gamma)}.
#'
#' @details
#' Write \eqn{u} for the coefficients the criterion integrates and
#' \eqn{\gamma} for the ones it estimates. At the mode,
#' \eqn{\partial u/\partial\gamma = -A_{uu}^{-1}A_{u\gamma}}, and a
#' hyperparameter moves \eqn{u} by
#' \eqn{-A_{uu}^{-1}\,\partial^2\rho/\partial u\,\partial\eta}. The
#' columns of `Bg` and `Be` are these movements written over the kept
#' coordinates, with a unit entry at \eqn{\gamma}'s own position.
#'
#' @param mc The list [fit_marginal_context()] returns.
#' @param coef,hyper The fitted coefficients and hyperparameters.
#' @param method The [OuterMethod()] the fit ran.
#' @param A The penalized information over the kept coordinates.
#' @param keep A logical vector over every coefficient, `TRUE` where it is
#'   kept.
#' @param full Whether to build `Be`, the movement with the
#'   hyperparameters.
#'
#' @return A list with `Pu`, `ui`, `p`, `Bg`, `Be`, `Hv`, `nh`, `ng` and
#'   `idx`, or `NULL` where a piece cannot be computed.
#'
#' @keywords internal
marginal_blocks <- function(mc, coef, hyper, method, A, keep, full = TRUE) {
  spec <- mc$spec
  design <- mc$design
  gam <- mc$gam
  kf <- which(keep)
  gi <- match(gam$where, kf)
  if (anyNA(gi)) return(NULL)
  ui <- setdiff(seq_along(kf), gi)
  A <- as_dense(A)
  Pu <- tryCatch(solve_pd(A[ui, ui, drop = FALSE], "the penalized information"),
                 error = function(e) NULL)
  if (is.null(Pu)) return(NULL)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  nh <- nrow(idx)
  ng <- length(gi)
  basis <- integrated_basis(spec, design, method@kind, gamma = TRUE)
  Hv <- tryCatch(statmod_marginal_hess(spec, design, coef, hyper, method, idx,
                                       basis, gam = gam),
                 error = function(e) NULL)
  if (is.null(Hv) || !all(is.finite(Hv))) return(NULL)
  p <- length(kf)
  Bg <- matrix(0, p, ng)
  for (j in seq_len(ng)) {
    Bg[gi[j], j] <- 1
    Bg[ui, j] <- -as.numeric(Pu %*% A[ui, gi[j]])
  }
  Be <- matrix(0, p, 0)
  if (full && nh) {
    X <- hyper_mode_cross(spec, design, coef, hyper, idx, length(keep))$cross
    X <- X[keep, , drop = FALSE]
    Be <- matrix(0, p, nh)
    Be[ui, ] <- -as.matrix(Pu %*% X[ui, , drop = FALSE])
  }
  list(Pu = Pu, ui = ui, p = p, Bg = Bg, Be = Be, Hv = as.matrix(Hv),
       nh = nh, ng = ng, idx = idx)
}


#' A Variance Matrix With Its Unreadable Coordinates Set to Zero
#'
#' @description
#' Sets to zero the rows and columns of a variance matrix whose diagonal is
#' not finite, which is how [hyper_variance()] marks a coordinate whose
#' curvature it could not read.
#'
#' @param V A variance matrix.
#'
#' @return `V` with those rows and columns set to zero.
#'
#' @keywords internal
zero_unread <- function(V) {
  ok <- is.finite(diag(V))
  if (!all(ok)) {
    V[!ok, ] <- 0
    V[, !ok] <- 0
  }
  V
}


#' The Correction for the Estimated Hyperparameters, With Estimated
#' Coefficients
#'
#' @description
#' The branch of [statmod_edf_correction()] for a fit whose criterion also
#' estimated coefficients. The correction is the trace of the information
#' against the part of the unconditional variance that the hyperparameters
#' add to the block variance of [marginal_vcov()].
#'
#' @details
#' With \eqn{B = (B_\eta, B_\gamma)} the movement of the integrated
#' coefficients, \eqn{V_v} the inverse of the negative Hessian of the
#' criterion over \eqn{(\eta, \gamma)} and \eqn{V_\gamma} the one over
#' \eqn{\gamma} alone, the correction is
#' \deqn{\mathrm{tr}\{(B V_v B^\top - B_\gamma V_\gamma B_\gamma^\top) H\}.}
#' The part for one hyperparameter is
#' \eqn{\mathrm{tr}(b_k v_{kk} b_k^\top H)}, with \eqn{b_k} its column of
#' \eqn{B_\eta}.
#'
#' @param mc The list [fit_marginal_context()] returns.
#' @param coef,hyper The fitted coefficients and hyperparameters.
#' @param method The [OuterMethod()] the fit ran.
#' @param idx The hyperparameters' index.
#' @param expected,approx Passed to [statmod_information_at()].
#' @param zero The answer where the correction cannot be computed.
#'
#' @return A list with `total`, `per` and `n_hyper`, as
#'   [statmod_edf_correction()] returns it.
#'
#' @keywords internal
marginal_edf_correction <- function(mc, coef, hyper, method, idx, expected,
                                    approx, zero) {
  spec <- mc$spec
  design <- mc$design
  H <- as_dense(statmod_information_at(spec, coef, design, expected, approx))
  S <- zap_nonfinite(statmod_penalty_at(spec, coef, hyper, design, "hessian"))
  A <- as_dense(H + S)
  keep <- rep(TRUE, nrow(A))
  # the flat coordinates only: the estimated ones stay, being what B moves
  keep[as.integer(attr(design, "pinned"))] <- FALSE
  mb <- marginal_blocks(mc, coef, hyper, method, A[keep, keep, drop = FALSE],
                        keep, full = TRUE)
  if (is.null(mb) || !mb$nh) return(zero)
  Vv <- hyper_variance(-mb$Hv)
  sel <- mb$nh + seq_len(mb$ng)
  Vg <- hyper_variance(-mb$Hv[sel, sel, drop = FALSE])
  if (is.null(Vv) || is.null(Vg)) return(zero)
  Vv <- zero_unread(Vv)
  Vg <- zero_unread(Vg)
  Hk <- H[keep, keep, drop = FALSE]
  B <- cbind(mb$Be, mb$Bg)
  D <- B %*% Vv %*% t(B) - mb$Bg %*% Vg %*% t(mb$Bg)
  total <- sum(D * t(Hk))
  contrib <- vapply(seq_len(mb$nh), function(k) {
    b <- mb$Be[, k, drop = FALSE]
    sum((b %*% Vv[k, k, drop = FALSE] %*% t(b)) * t(Hk))
  }, numeric(1))
  per <- tapply(contrib, paste(idx$parameter, idx$term, sep = "\r"), sum)
  list(total = total, per = per, n_hyper = nrow(idx))
}


#' The Criterion a Fit Maximized, Rebuilt From the Fit
#'
#' @description
#' The specification and the design on which a fit's marginal criterion was
#' maximized: the coefficients [marginal_coords()] names are held at their
#' fitted values, and the coordinates left out of the determinant are
#' recorded on the design.
#'
#' @details
#' A consumer that reads the criterion's derivatives at a fitted object has
#' to read them on this criterion. On the fit's own specification, which
#' holds nothing, the same functions return the derivatives of a criterion
#' that integrates the estimated coefficients, which is another function.
#'
#' @param spec The fit's specification.
#' @param design Its design.
#' @param coef The fitted coefficients, a named list.
#' @param method The [OuterMethod()] the fit ran, or `NULL`.
#' @param pinned The coordinates left out of the determinant, as the fit
#'   records them in `methods$pinned`.
#'
#' @return A list with `spec`, `design` and `gam`, the last as
#'   [marginal_coords()] returns it.
#'
#' @keywords internal
fit_marginal_context <- function(spec, design, coef, method, pinned = NULL) {
  gam <- tryCatch(marginal_coords(spec, design, method),
                  error = function(e) marginal_coords(spec, design, NULL))
  if (length(gam$where)) {
    beta <- unlist(coef[spec@distrib@params], use.names = FALSE)
    spec <- statmod_hold(spec, gam, beta[gam$where])
  }
  if (length(pinned)) attr(design, "pinned") <- as.integer(pinned)
  list(spec = spec, design = design, gam = gam)
}
