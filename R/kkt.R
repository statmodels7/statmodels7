#' Coordinates at the Edge of Their Chart That Point Inward
#'
#' @description
#' Checks the first-order condition for a maximum at every free coordinate
#' that has run past the edge of a chart mapping onto a bounded set, reading
#' it on the bounded scale rather than on the free one. Returns the
#' coordinates where moving back inside raises the log-likelihood, so the
#' point the fit stopped at is not a maximum.
#'
#' @details
#' At a coordinate \eqn{\eta} whose chart \eqn{\theta = h(\eta)} saturates,
#' \eqn{\partial\ell/\partial\eta = (\partial\ell/\partial\theta)\,h'(\eta)},
#' and \eqn{h'} tends to zero at the edge whatever
#' \eqn{\partial\ell/\partial\theta} is. A fit can therefore report a
#' vanishing score, and a certificate a vanishing mode error, at a point
#' where the log-likelihood still rises towards the interior. The condition
#' for a maximum on the boundary of a bounded parameter (Karush, Kuhn and
#' Tucker) is that the derivative on the bounded scale points outward, and
#' that is what is checked.
#'
#' The derivative is read as a one-sided difference: each coordinate beyond
#' `edge` in absolute value is moved to \eqn{\pm}`edge`, the rest held, and
#' the log-likelihood is evaluated there. A coordinate is reported where the
#' gain exceeds `tol`, in log-likelihood units. The coordinates checked are
#' those [modelterms7::term_charted()] names for a structural term and for
#' [modelterms7::nl()], and the intercept of an equation that has no other
#' column and carries a link that is not the identity. None of them carries
#' a penalty, so the log-likelihood difference is the objective's.
#'
#' @param spec,design,coef The specification, its design and the
#'   coefficients by parameter, as [statmod()] holds them at the end of a fit.
#' @param edge The free value past which a coordinate is at the edge of its
#'   chart, the same default as [statmod_certificate()]'s.
#' @param tol The gain, in log-likelihood units, above which the point is
#'   not a maximum: [mode_error_limit()].
#'
#' @return A data frame with one row per violating coordinate and columns
#'   `kind` (`"structural"` or `"coefficient"`), `param`, `term`, `name`,
#'   `eta`, `target` and `gain`; zero rows where every coordinate at an edge
#'   points outward. The design's structural state is left as it was found.
#'
#' @seealso [statmod_certificate()], which reports a violation.
#' @keywords internal
edge_violations <- function(spec, design, coef, edge = 8,
                            tol = mode_error_limit()) {
  out <- data.frame(kind = character(0), param = character(0),
                    term = character(0), name = character(0),
                    eta = numeric(0), target = numeric(0), gain = numeric(0))
  base <- tryCatch(statmod_loglik_at(spec, coef, design),
                   error = function(e) NA_real_)
  if (!is.finite(base)) return(out)
  gain_at <- function(cf) {
    l <- tryCatch(statmod_loglik_at(spec, cf, design),
                  error = function(e) NA_real_)
    if (is.finite(l)) l - base else NA_real_
  }
  add <- function(kind, param, term, name, eta, g) {
    if (is.finite(g) && g > tol) {
      out[nrow(out) + 1L, ] <<- list(kind, param, term, name, eta,
                                     sign(eta) * edge, g)
    }
  }

  # a structural term's own parameters live in the design's state
  sst <- statmod_structural_state(design)
  for (u in attr(design, "structural")) {
    tm <- spec@terms[[u$param]][[u$term]]
    ch <- setdiff(modelterms7::term_charted(tm), sst$held[[u$term]])
    z0 <- sst$zeta[[u$term]]
    for (j in intersect(ch, names(z0))) {
      if (!(abs(z0[[j]]) > edge)) next
      z <- z0
      z[[j]] <- sign(z0[[j]]) * edge
      sst$zeta[[u$term]] <- z
      g <- gain_at(coef)
      sst$zeta[[u$term]] <- z0
      add("structural", u$param, u$term, j, z0[[j]], g)
    }
  }

  # coefficients that are the free value of a chart
  for (p in spec@distrib@params) {
    d <- design[[p]]
    pos <- integer(0)
    if (identical(d$npar, 1L) && identical(d$coef_names, "(Intercept)") &&
        !identical(spec@distrib@link_params[[p]]@link_name, "identity")) {
      pos <- stats::setNames(1L, "(Intercept)")
    }
    for (k in names(spec@terms[[p]])) {
      tm <- spec@terms[[p]][[k]]
      if (!S7::S7_inherits(tm, modelterms7::NlTerm)) next
      ch <- modelterms7::term_charted(tm)
      if (length(ch)) pos <- c(pos, stats::setNames(d$blocks[[k]][ch],
                                                    d$coef_names[d$blocks[[k]][ch]]))
    }
    for (i in seq_along(pos)) {
      b0 <- coef[[p]][[pos[[i]]]]
      if (!(abs(b0) > edge)) next
      cf <- coef
      cf[[p]][[pos[[i]]]] <- sign(b0) * edge
      add("coefficient", p, "", names(pos)[i], b0, gain_at(cf))
    }
  }
  out
}

#' Whether a Restart From the Edge Improved the Fit
#'
#' @description
#' Compares a fit restarted by [statmod()] after [edge_violations()] with the
#' one it replaces: on the criterion that chose the hyperparameters where one
#' ran, in its own direction, and on the penalized objective otherwise.
#'
#' @param new,old The two results of [statmod_alternate()] or
#'   [statmod_select()].
#' @param crit The criterion `old` reached, `NA` where none ran.
#' @param outer_criterion,sparse_criterion As passed to [statmod()].
#'
#' @return `TRUE` where `new` is strictly better.
#' @keywords internal
edge_restart_better <- function(new, old, crit, outer_criterion,
                                sparse_criterion) {
  cr <- if (!is.null(outer_criterion)) outer_criterion else sparse_criterion
  if (!is.null(cr) && is.finite(crit) && !is.null(new$criterion) && is.finite(new$criterion)) {
    if (outer_minimize(cr)) new$criterion < crit else new$criterion > crit
  } else {
    is.finite(new$value) && new$value < old$value
  }
}
