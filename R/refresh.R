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
#' @param spec A \code{\link{StatmodSpec}}.
#'
#' @return A list of entries with \code{param} and \code{term}, possibly
#'   empty.
#'
#' @seealso \code{\link{statmod_design_at}}, \code{\link{refreshes_own_block}}
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
#' \code{\link{statmod_commit_refresh}} is called, so the trial points of a
#' line search all see the same schedule and the schedule advances once per
#' sweep.
#'
#' The result is memoized on the coefficients, since the objective, its
#' gradient and its curvature are asked for at the same point in turn.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design, as \code{\link{statmod_design}} returns it.
#'
#' @return A design.
#'
#' @seealso \code{\link{statmod_commit_refresh}}
#'
#' @keywords internal
statmod_design_at <- function(spec, coef, design) {
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(design)
  st <- attr(design, "state")
  key <- unlist(coef, use.names = FALSE)
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
#' The value a term reports is unchanged by the schedule -- a break-point is
#' read off the coefficients and the rescaling reaches only the columns -- so
#' committing does not move the objective at the same coefficients.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return \code{TRUE} when there was something to commit, invisibly.
#'
#' @seealso \code{\link{statmod_design_at}}
#'
#' @keywords internal
statmod_commit_refresh <- function(spec, coef, design) {
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(invisible(FALSE))
  st <- attr(design, "state")
  for (r in rf) {
    cols <- design[[r$param]]$blocks[[r$term]]
    st$terms[[r$param]][[r$term]] <-
      modelterms7::term_refresh(st$terms[[r$param]][[r$term]],
                                coef[[r$param]][cols])
  }
  st$key <- NULL
  st$value <- NULL
  invisible(TRUE)
}


#' Have the Refreshable Terms Settled?
#'
#' @description
#' \code{TRUE} when every term that recomputes its own block reports that its
#' own iteration has nothing further to say.
#'
#' @details
#' The question cannot always be answered by the score. Where a term's block
#' is the Jacobian of its contribution, the gradient of the model's objective
#' is the model's and its vanishing is the test; where the block is a working
#' linearization with a frozen weight, as in a discontinuous break-point
#' term, it is not, and the profile objective there is a step function in the
#' break-point with no gradient to vanish. Measured on
#' \code{jump()}: the fit reaches the break-point and the jump size to three
#' figures, the objective stops moving at the twelfth digit, and the score of
#' the working model stays at 0.176 forever.
#' \code{\link[modelterms7]{term_converged}} is what each construction
#' answers instead.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#'
#' @return A single logical; \code{TRUE} when there is nothing to ask.
#'
#' @seealso \code{\link{statmod_design_at}}
#'
#' @keywords internal
statmod_refresh_settled <- function(spec, design) {
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(TRUE)
  st <- attr(design, "state")
  for (r in rf) {
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
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return A \code{\link{StatmodSpec}}.
#'
#' @seealso \code{\link{statmod_commit_refresh}}
#'
#' @keywords internal
statmod_fitted_spec <- function(spec, coef, design) {
  sst <- statmod_structural_state(design)
  if (!is.null(sst)) spec@structural <- sst$zeta
  rf <- attr(design, "refresh")
  if (is.null(rf) || !length(rf)) return(spec)
  statmod_commit_refresh(spec, coef, design)
  st <- attr(design, "state")
  tms <- spec@terms
  for (r in rf) {
    tms[[r$param]][[r$term]] <- st$terms[[r$param]][[r$term]]
  }
  spec@terms <- tms
  spec
}
