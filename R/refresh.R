#' @include spec.R
NULL

# A term whose design block depends on its own coefficients is carried here.
# The block is a Jacobian and the contribution is not the block times the
# coefficients, so two things follow: the predictor takes term_value() and not
# X beta, and the scoring step's increment is the Gauss-Newton one. Writing
# the difference as a per-observation adjustment,
#
#   eta = X(beta) beta + adj(beta),   adj = term_value(beta) - X(beta) beta,
#
# leaves every crossprod in objective.R reading X as it did, since the
# adjustment does not depend on beta to first order and the derivative of the
# predictor in the coefficients is the block itself.

#' Which Terms Recompute Their Own Block
#'
#' @description
#' The parameter and name of every term whose design block is a function of
#' its own coefficients, in the order the design holds them.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A list of entries with `param` and `term`, possibly
#'   empty.
#'
#' @seealso [statmod_design_at()], [refreshes_own_block()]
#'
#' @keywords internal
statmod_refreshable <- function(spec) {
  out <- list()
  for (p in names(spec@terms)) {
    for (nm in names(spec@terms[[p]])) {
      if (refreshes_own_block(spec@terms[[p]][[nm]])) {
        out[[length(out) + 1L]] <- list(param = p, term = nm)
      }
    }
  }
  out
}


#' The Design at Given Coefficients
#'
#' @description
#' The design with every refreshable term's block recomputed at the
#' coefficients it currently holds, and the difference between its
#' contribution and its linearization carried as an adjustment to the
#' predictor.
#'
#' @details
#' A design with no refreshable term is returned unchanged, so a model of
#' ordinary terms pays nothing and reaches exactly the same arithmetic as
#' before.
#'
#' The refresh is CHAINED from the term the state holds rather than from the
#' specification, because the rescaling factor of a discontinuous break-point
#' term is a state of the iteration and not a function of the coefficients:
#' it halves when the break-point reverses direction, which is a fact about
#' the path and not about the point. The state advances only when
#' [statmod_commit_refresh()] is called, so the trial points of a
#' line search all see the same schedule and the schedule advances once per
#' sweep.
#'
#' The result is memoized on the coefficients, since the objective, its
#' gradient and its curvature are asked for at the same point in turn.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design, as [statmod_design()] returns it.
#'
#' @return A design.
#'
#' @seealso [statmod_commit_refresh()]
#'
#' @keywords internal
statmod_design_at <- function(spec, coef, design) {
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(design)
  st <- attr(design, "state")
  key <- c(unlist(coef, use.names = FALSE),
           as.numeric(isTRUE(st$working)))
  if (!is.null(st$key) && identical(st$key, key)) return(st$value)

  out <- design
  for (p in unique(vapply(rf, function(r) r$param, character(1)))) {
    out[[p]]$adj <- numeric(spec@n_obs)
  }
  # On other rows nothing is refreshed: a term's break-points and nonlinear
  # parameters are what the fit left them, and reading them off the new rows
  # would be the rebuild term_predict() exists to avoid. What the predictor
  # still needs there is the contribution, which the block times the
  # coefficients is not.
  nd <- spec@newdata
  for (r in rf) {
    p <- r$param
    cols <- design[[p]]$blocks[[r$term]]
    b <- coef[[p]][cols]
    tm <- st$terms[[p]][[r$term]]
    if (isTRUE(r$frozen)) {
      # A FROZEN working block is read as the fit last committed it and is
      # never refreshed at a trial point. Reading it chained instead made
      # the objective a moving target with a step at every data point the
      # implied break-point crosses: the line search inside the inner
      # optimizer then rejected the fixed-point iteration's own steps
      # (measured, a jseg fitted from the TRUE break-points landed at an
      # rss worse than the mean-only fit), and for a jseg the incremental
      # quadratic read-off took a further hidden step at every commit, so
      # the objective changed at unchanged coefficients. During the
      # working fit of fit_working() the block IS the model (eta = X beta,
      # no adjustment, the plain linear fit of Fasola et al.); everywhere
      # else the term contributes its committed value, which is what
      # "held fixed while another block moves" means for it.
      if (is.null(nd)) {
        Xt <- modelterms7::term_matrix(tm)
        out[[p]]$X[, cols] <- Xt
        val <- as.numeric(modelterms7::term_value(tm))
      } else {
        Xt <- design[[p]]$X[, cols, drop = FALSE]
        val <- as.numeric(modelterms7::term_value(tm, coef = b, newdata = nd))
      }
      if (!isTRUE(st$working)) {
        # anchored at the coefficients the term was COMMITTED at, so eta
        # stays linear in the current ones through the frozen block: at the
        # committed point the two coincide and eta is the contribution
        # exactly, while a b-dependent adjustment would make eta constant
        # in these columns and flat to anything that moves them
        bc <- tryCatch(tm@blueprint$coef, error = function(e) NULL)
        if (is.null(bc) || length(bc) != length(b)) bc <- b
        out[[p]]$adj <- out[[p]]$adj + (val - as.numeric(Xt %*% bc))
      }
      next
    }
    if (is.null(nd)) {
      tm <- modelterms7::term_refresh(tm, b)
      Xt <- modelterms7::term_matrix(tm)
      out[[p]]$X[, cols] <- Xt
      val <- as.numeric(modelterms7::term_value(tm))
    } else {
      Xt <- design[[p]]$X[, cols, drop = FALSE]
      val <- as.numeric(modelterms7::term_value(tm, coef = b, newdata = nd))
    }
    out[[p]]$adj <- out[[p]]$adj + (val - as.numeric(Xt %*% b))
  }
  st$key <- key
  st$value <- out
  out
}


#' Advance the Refresh State
#'
#' @description
#' Replaces the terms the design refreshes from by the ones they refresh to
#' at the given coefficients, so that whatever a term carries as a state of
#' its iteration moves on one step.
#'
#' @details
#' What moves is the rescaling factor of a discontinuous break-point term and
#' the direction it last travelled in, which \pkg{modelterms7} halves on a
#' reversal: a schedule that advanced once per objective evaluation would
#' anneal at the speed of the line search rather than at the speed of the
#' fit, and one that never advanced would solve a permanently smoothed
#' problem, whose fixed point is not the model's.
#'
#' For a term whose block is a Jacobian the value it reports is unchanged by
#' the schedule -- a break-point is read off the coefficients and the
#' rescaling reaches only the columns -- so committing does not move the
#' objective at the same coefficients. For a FROZEN working block that
#' sentence is false in two ways, which is why those terms are committed by
#' [fit_working()] and skipped here: a jseg's quadratic read-off
#' is incremental in the committed position, so a second commit at the same
#' coefficients takes a second step, and a refresh may relabel crossed
#' break-point lineages, after which the caller's coefficients are stale.
#' The relabeling is why the COMMITTED coefficients are returned: a caller
#' continues from what the terms stored, not from what it passed in.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#' @param which Which refresh entries to commit: `"all"`,
#'   `"jacobian"` (the default at the alternation's pass level, where
#'   the frozen ones are already committed by their own phase) or
#'   `"frozen"`.
#'
#' @return The coefficient list, with each committed term's stretch replaced
#'   by the coefficients the term stored, invisibly.
#'
#' @seealso [statmod_design_at()], [fit_working()]
#'
#' @keywords internal
statmod_commit_refresh <- function(spec, coef, design, which = "all") {
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(invisible(coef))
  st <- attr(design, "state")
  for (r in rf) {
    if (identical(which, "jacobian") && isTRUE(r$frozen)) next
    if (identical(which, "frozen") && !isTRUE(r$frozen)) next
    cols <- design[[r$param]]$blocks[[r$term]]
    tm <- modelterms7::term_refresh(st$terms[[r$param]][[r$term]],
                                    coef[[r$param]][cols])
    st$terms[[r$param]][[r$term]] <- tm
    bc <- tryCatch(tm@blueprint$coef, error = function(e) NULL)
    if (!is.null(bc) && length(bc) == length(cols)) {
      coef[[r$param]][cols] <- as.numeric(bc)
    }
  }
  st$key <- NULL
  st$value <- NULL
  invisible(coef)
}


#' Have the Refreshable Terms Settled?
#'
#' @description
#' `TRUE` when every term that recomputes its own block reports that its
#' own iteration has nothing further to say.
#'
#' @details
#' The question cannot always be answered by the score. Where a term's block
#' is the Jacobian of its contribution, the gradient of the model's objective
#' is the model's and its vanishing is the test; where the block is a working
#' linearization with a frozen weight, as in a discontinuous break-point
#' term, it is not, and the profile objective there is a step function in the
#' break-point with no gradient to vanish. Measured on
#' `jump()`: the fit reaches the break-point and the jump size to three
#' figures, the objective stops moving at the twelfth digit, and the score of
#' the working model stays at 0.176 forever.
#' [modelterms7::term_converged()] is what each construction
#' answers instead.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param which Which refresh entries to ask: `"all"`,
#'   `"jacobian"` or `"frozen"`.
#'
#' @return A single logical; `TRUE` when there is nothing to ask.
#'
#' @seealso [statmod_design_at()]
#'
#' @keywords internal
statmod_refresh_settled <- function(spec, design, which = "all") {
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(TRUE)
  st <- attr(design, "state")
  for (r in rf) {
    if (identical(which, "jacobian") && isTRUE(r$frozen)) next
    if (identical(which, "frozen") && !isTRUE(r$frozen)) next
    if (!isTRUE(modelterms7::term_converged(st$terms[[r$param]][[r$term]]))) {
      return(FALSE)
    }
  }
  TRUE
}


#' The Terms as the Fit Left Them
#'
#' @description
#' The specification with every refreshable term replaced by the one the fit
#' arrived at, so that a break-point, a nonlinear parameter and the block
#' they imply are read off the fitted object rather than off the
#' specification it started from.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return A [StatmodSpec()].
#'
#' @seealso [statmod_commit_refresh()]
#'
#' @keywords internal
statmod_fitted_spec <- function(spec, coef, design) {
  sst <- statmod_structural_state(design)
  if (!is.null(sst)) spec@structural <- sst$zeta
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(spec)
  # jacobian entries only: a frozen term was committed by its own phase at
  # exactly these coefficients, and committing a jseg again at the same
  # point takes a further step of its incremental read-off, so the reported
  # break-points would not be the fitted ones
  statmod_commit_refresh(spec, coef, design, which = "jacobian")
  st <- attr(design, "state")
  tms <- spec@terms
  for (r in rf) {
    tms[[r$param]][[r$term]] <- st$terms[[r$param]][[r$term]]
  }
  spec@terms <- tms
  spec
}
