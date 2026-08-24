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
#' Locates every term of the specification whose design block is a function
#' of its own coefficients, and reports where each one sits. These are the
#' terms whose block has to be rebuilt whenever the coefficients move:
#' [modelterms7::seg()], [modelterms7::jump()], [modelterms7::jseg()] and
#' [modelterms7::nl()]. Every other term's block is fixed once at build
#' time.
#'
#' @details
#' Membership is decided by asking whether the term registers a
#' `term_refresh()` method of its own, which is [refreshes_own_block()]. The
#' base method on `model_term` is the identity, so a term written later is
#' covered without an edit here.
#'
#' This is the list every other function in the file walks. Its emptiness is
#' what makes a model of ordinary terms pay nothing for the refresh
#' machinery.
#'
#' @param spec A [StatmodSpec()], whose `terms` are walked equation by
#'   equation.
#'
#' @return A list with one element per refreshable term, in the order the
#'   design holds them: equations in the family's parameter order, and terms
#'   within an equation in the order they were written. Each element is a
#'   list of two:
#'   \describe{
#'     \item{`param`}{the distribution parameter whose equation the term sits
#'       in, a string.}
#'     \item{`term`}{the term's index within that equation, an integer.}
#'   }
#'   An empty list when no term refreshes, which is the common case.
#'
#' @seealso [statmod_design_at()] for the rebuild,
#'   [refreshes_own_block()] for the predicate,
#'   [statmod_commit_refresh()] for advancing the state.
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
#' Recomputes every refreshable term's design block at the coefficients
#' currently in hand, and carries the difference between what the term
#' contributes and what its linearization contributes as a per-observation
#' adjustment to the predictor. The adjustment is what lets every
#' cross-product in the objective read the block as an ordinary design while
#' the predictor stays exact:
#'
#' \deqn{\mathrm{adj} = \mathrm{term\_value}(\beta) - X(\beta)\,\beta}
#'
#' For a term whose block is a Jacobian, [modelterms7::seg()] and
#' [modelterms7::nl()], the adjustment is non-zero and the resulting scoring
#' step is the Gauss-Newton one. For [modelterms7::jump()] the columns
#' satisfy \eqn{X\beta = \mathrm{value}} exactly, so the adjustment is zero
#' and the same step is Fasola's fixed-point iteration.
#'
#' @details
#' # A model with no refreshable term pays nothing
#'
#' [refresh_units()] is empty there, the design is returned as it arrived,
#' and the arithmetic downstream is untouched.
#'
#' # Chained from the term, not from the specification
#'
#' The refresh reads the term the design state currently holds, and not the
#' one the specification was built with. The reason is that a discontinuous
#' break-point term carries a rescaling factor that is a state of the
#' iteration and not a function of the coefficients: \pkg{modelterms7} halves
#' it whenever the break-point reverses direction, which is a fact about the
#' path taken and not about the point reached. Refreshing from the
#' specification each time would reset that factor to its starting value and
#' solve a permanently smoothed problem, whose fixed point is not the
#' model's.
#'
#' The state advances only when [statmod_commit_refresh()] is called, so
#' every trial point of a line search sees one schedule and the schedule
#' advances once per sweep.
#'
#' # Memoization
#'
#' The result is cached on the coefficients, because the objective, its
#' gradient and its curvature are all asked for at the same point in turn
#' and each would otherwise rebuild the same blocks.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter, each as long as its equation's design is wide.
#' @param design The design to refresh, as [statmod_design()] returns it.
#'
#' @return A design of the same shape as `design`, with the refreshable
#'   terms' column blocks replaced and each equation's `adj` set to the
#'   per-observation adjustment above. `design` itself when nothing
#'   refreshes.
#'
#' @seealso [statmod_commit_refresh()] to advance the state afterwards,
#'   [refresh_units()] for the terms this walks,
#'   [statmod_refresh_settled()] for the verdict.
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
#' # What moves
#'
#' The rescaling factor of a discontinuous break-point term, and the
#' direction it last travelled in, which \pkg{modelterms7} halves on a
#' reversal. Advancing that once per objective evaluation would anneal at
#' the speed of the line search instead of the speed of the fit; never
#' advancing it would solve a permanently smoothed problem, whose fixed
#' point is not the model's. Once per sweep is what this call is for.
#'
#' # Why a frozen block is committed elsewhere
#'
#' Where a term's block is a Jacobian, committing does not move the
#' objective at the same coefficients: the break-point is read off the
#' coefficients and the rescaling reaches only the columns.
#'
#' Where the block is a frozen working linearization, that is false twice
#' over, which is why those terms are committed by [fit_working()] and
#' skipped here at the default. A [modelterms7::jseg()] reads its position
#' from a quadratic that is incremental in the position already committed,
#' so a second commit at the same coefficients takes a further hidden step,
#' measured at up to 0.71 per observation on the contribution. And a refresh
#' may relabel crossed break-point lineages, after which the caller's own
#' copy of the coefficients names them in the old order.
#'
#' The relabeling is why the committed coefficients are returned. A caller
#' continues from what the terms stored, not from what it passed in.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter.
#' @param design The design, whose refresh state is what this advances.
#' @param which Which entries to commit: `"all"` (the default), `"jacobian"`
#'   or `"frozen"`. The alternation's pass level and [statmod_fitted_spec()]
#'   both pass `"jacobian"`, the frozen terms having been committed already
#'   by their own phase.
#'
#' @return The coefficient list, invisibly, with each committed term's
#'   stretch replaced by the coefficients that term stored. Identical to
#'   `coef` when nothing was committed or when no term relabeled.
#'
#' @seealso [statmod_design_at()] for the refresh this advances,
#'   [fit_working()] for the phase that commits the frozen terms,
#'   [statmod_refresh_settled()] for the verdict.
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
#' Asks every term that recomputes its own block whether its own iteration
#' has anything further to say, through
#' [modelterms7::term_converged()], and returns `TRUE` when none of them
#' does. This is the verdict for the refreshable half of a fit.
#'
#' @details
#' # Why the score cannot answer this
#'
#' Where a term's block is the Jacobian of its contribution, the gradient of
#' the model's objective is the model's own gradient and its vanishing is
#' the test. Where the block is a working linearization with a frozen
#' weight, as in a discontinuous break-point term, the gradient belongs to
#' the working model and not to the objective, and the profile objective is
#' a step function in the break-point with no gradient to vanish at all.
#'
#' Measured on `y ~ jump(x)` at \eqn{n = 400}, a Gaussian response with a
#' step of 2 at \eqn{x = 6}: the fit recovers the position at 6.004 and
#' reports `converged = TRUE`, while the score of the working model at that
#' point is \eqn{7.8 \times 10^{7}}. The size is the annealed rescaling
#' factor, whose auxiliary columns grow as the schedule tightens. A rule
#' reading that score would never stop.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, whose refresh state holds the terms asked.
#' @param which Which entries to ask: `"all"` (the default), `"jacobian"` or
#'   `"frozen"`.
#'
#' @return A single logical. `TRUE` when every term asked reports it has
#'   settled, and `TRUE` when there is nothing to ask, so a model with no
#'   refreshable term never blocks a fit's verdict on this.
#'
#' @seealso [modelterms7::term_converged()] for what each construction
#'   answers, [statmod_design_at()] for the refresh,
#'   [statmod_commit_refresh()] for the state it reads.
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
#' Returns the specification with every refreshable term replaced by the one
#' the fit arrived at. A break-point, a nonlinear parameter and the design
#' block they imply are then read off the fitted object, and the
#' specification a caller passed to [statmod()] no longer decides what
#' `summary()`, `predict()` or [modelterms7::seg_psi()] report.
#'
#' @details
#' This is what a [StatmodFit()] stores in its `spec` property. The commit is
#' `which = "jacobian"`: a frozen block was committed by its own phase at
#' exactly these coefficients, and committing a [modelterms7::jseg()] again
#' at the same point would take a further step of its incremental read-off,
#' so the break-points reported would not be the fitted ones. Both kinds of
#' term are then copied across from the design's state, committed or not.
#'
#' A structural term's own parameters are copied too, from the design's
#' structural state onto `spec@structural`, so a filter's persistence and
#' loadings are read off the fit as well.
#'
#' @param spec A [StatmodSpec()], the one the fit started from.
#' @param coef A named list of coefficient vectors, the ones the fit reached.
#' @param design The design at those coefficients, whose refresh state holds
#'   the terms to copy across.
#'
#' @return A [StatmodSpec()] identical to `spec` except that each
#'   refreshable term is the object the fit left behind.
#'
#' @seealso [statmod_commit_refresh()] for the commit this performs,
#'   [refresh_units()] for the terms replaced.
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
