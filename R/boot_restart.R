#' @include statmod.R
NULL

#' How Many Restarts the Terms Ask For
#'
#' @description
#' Returns the largest `n_boot` any break-point term of the specification
#' declares, and zero when the model carries no such term. This is the
#' budget [statmod_boot_restart()] spends: how many proposals it may try
#' before giving up on improving the fit.
#'
#' @details
#' The number is declared on the term. [modelterms7::seg()],
#' [modelterms7::jump()] and [modelterms7::jseg()] each take `n_boot`, with
#' a default of 10, because those are the terms whose objective has the
#' spurious local optima the device exists for. Running the restarts belongs
#' here, in the layer that can refit the model, and the split follows the
#' one a penalty already uses: the term says what it needs, this package
#' does it.
#'
#' A model with two break-point terms asking for 10 and 25 gets 25. The
#' restart loop works on all of them together, since a proposal moves every
#' break-point in the model at once, so the budget is one number and the
#' largest request is the one honored.
#'
#' @param spec A [StatmodSpec()], whose terms are walked in every equation.
#'
#' @return A single non-negative integer. Zero when no term declares
#'   `n_boot`, and a budget of zero turns the restart loop off.
#'
#' @seealso [statmod_boot_restart()], which spends this budget,
#'   [modelterms7::seg()] for where the number is set.
#'
#' @keywords internal
seg_boot_total <- function(spec) {
  nb <- 0L
  for (p in names(spec@terms)) {
    for (tm in spec@terms[[p]]) {
      if (!S7::S7_inherits(tm, modelterms7::SegTerm)) next
      v <- tryCatch(tm@spec$n_boot, error = function(e) NULL)
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
        nb <- max(nb, as.integer(v))
      }
    }
  }
  nb
}


#' Restarting Around a Fitted Model, Screened on the Exact Profile
#'
#' @description
#' Improves a fitted break-point model by bootstrap restarting (Wood 2001),
#' screening each proposal on the exact profile before paying for a refit.
#' The objective of a model with break-points has spurious local minima, and
#' an ordinary fit converges into whichever one its starting positions sit
#' in; this searches for a better one.
#'
#' @details
#' # What makes the screen cheap
#'
#' The non-convexity of a break-point model lives entirely in the positions.
#' Hold them and everything left is convex, so the exact profile at a fixed
#' configuration of positions is one linear fit. A proposal is therefore a
#' configuration of positions and nothing else, and two proposals are ranked
#' by their profiles at a cost of one linear fit each. Only a proposal the
#' profile prefers earns a refit of the whole model, and the refit's answer
#' is accepted or rejected on the true objective.
#'
#' Measured at \eqn{n = 10^4}: a proposal that goes nowhere costs about half
#' a second, against 5 to 15 seconds for the refit it would otherwise have
#' triggered. The design this replaced refitted every proposal and spent
#' 945 seconds re-verifying an optimum the sweep had already found.
#'
#' # The three proposal kinds
#'
#' Tried in this order.
#'
#' 1. **The deterministic sweep.** Each break-point in turn is swept over a
#'    grid on the profile with the others held, which is
#'    [modelterms7::seg_polish()]. This walks straight to a feature the
#'    fitting iteration pressed a break-point away from.
#' 2. **The bootstrap sweep.** The same descent on the profile of a
#'    resample, the multinomial counts entering as weights, which moves the
#'    profile's optima the way refitting the resample would.
#' 3. **The random sweep.** The same descent from positions drawn uniformly
#'    over the confinement interval.
#'
#' The two stochastic kinds alternate. Four consecutive proposals that fail
#' the screen end the loop, whatever budget is left.
#'
#' # The profile is exact for an identity link and a proposal elsewhere
#'
#' It reads the response net of the other contributions in the term's
#' equation, on the predictor scale. For a Gaussian response and an identity
#' link that is the model's own least-squares objective. For anything else
#' it is an approximation used to rank candidates, and the true objective
#' decides the acceptance, so a poor ranking costs time and never
#' correctness.
#' [modelterms7::seg_start()] makes the same argument for the same reason.
#'
#' # Reproducibility
#'
#' The draws come from the session's generator, so a fit with restarts
#' repeats under [set.seed()]. The refreshable and structural state of the
#' design is snapshotted at the incumbent and restored whenever a candidate
#' loses, so a rejected proposal leaves nothing behind.
#'
#' @param spec The [StatmodSpec()] being fitted.
#' @param design The assembled design, as [statmod_design()] returns it.
#' @param blocks The split of the terms into the jointly fitted smooth block
#'   and the kinked ones, as [statmod_blocks()] returns it.
#' @param hyper The hyperparameters the fit ended at, held fixed throughout:
#'   the restarts search over positions, not over hyperparameters.
#' @param inner_optimizer How the smooth block is fitted, [iwls()] or an
#'   \pkg{optimizers7} optimizer.
#' @param res The fitted result to improve, as [statmod_alternate()] returns
#'   it. Returned unchanged when nothing better is found.
#' @param expected,approx,maxit,tol Passed to [statmod_alternate()] for each
#'   refit, with the same meanings they have there.
#' @param vb The resolved verbosity, as [verbosity()] returns it.
#' @param nb The budget: at most how many proposals to try, as
#'   [seg_boot_total()] reports it. A budget of zero returns `res`
#'   untouched.
#'
#' @return `res`, with `par`, `value`, `converged`, `obj` and the block
#'   histories replaced when a restart improved the objective. Everything
#'   else is carried over untouched, an outer search's history and its
#'   optimizer among them, so the result is the same shape either way.
#'
#' @references
#' Wood, S. N. (2001). Minimizing model fitting objectives that contain
#' spurious local minima by bootstrap restarting. *Biometrics*, 57(1),
#' 240--244.
#'
#' @seealso [seg_boot_total()], [statmod_alternate()]
#'
#' @keywords internal
statmod_boot_restart <- function(spec, design, blocks, hyper, inner_optimizer,
                                 res, expected, approx, maxit, tol, vb, nb) {
  st <- attr(design, "state")
  if (is.null(st)) return(res)
  y0 <- spec@response
  if (!is.numeric(y0) || is.matrix(y0) || length(y0) != spec@n_obs ||
      anyNA(y0)) {
    return(res)
  }
  sst <- attr(design, "structure")
  n <- spec@n_obs
  vbq <- lapply(vb, function(x) FALSE)
  obj_split <- res$obj$split
  obj_stack <- res$obj$stack
  snap <- function() {
    list(terms = st$terms, zeta = if (!is.null(sst)) sst$zeta)
  }
  restore <- function(s) {
    st$terms <- s$terms
    st$key <- NULL
    st$value <- NULL
    if (!is.null(sst)) {
      sst$zeta <- s$zeta
      sst$key <- NULL
      sst$value <- NULL
    }
  }
  best <- res
  keep <- snap()
  fields <- c("par", "value", "converged", "obj", "hist_blocks", "hist_inner")

  # the response net of what the rest of the term's equation contributes,
  # on the predictor scale
  net_y <- function(r, tm) {
    d1 <- statmod_design_at(spec, obj_split(best$par), design)
    p <- r$param
    cf <- obj_split(best$par)[[p]]
    eta <- as.numeric(d1[[p]]$X %*% cf)
    if (!is.null(d1[[p]]$adj)) eta <- eta + d1[[p]]$adj
    y0 - (eta - as.numeric(modelterms7::term_value(tm)))
  }

  # One proposal: new positions for every break-point term, each kept only
  # when its own profile improves. NULL when nothing improved anywhere,
  # which is what makes a dry round cost linear fits and no refit.
  propose <- function(kind) {
    par <- best$par
    gain <- FALSE
    for (r in attr(design, "refresh")) {
      tm <- st$terms[[r$param]][[r$term]]
      if (!S7::S7_inherits(tm, modelterms7::SegTerm)) next
      # per term, so a term the profile machinery rejects -- a developed
      # break-point has one position per observation and no single profile
      # -- leaves the other terms' proposals standing
      cand <- tryCatch({
        yn <- net_y(r, tm)
        base <- modelterms7::seg_profile_rss(tm, yn)
        switch(kind,
          sweep = modelterms7::seg_polish(tm, yn),
          boot = {
            w <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
            modelterms7::seg_polish(tm, yn, weights = spec@weights * w)
          },
          random = {
            lim <- tm@blueprint$lim
            modelterms7::seg_polish(
              modelterms7::seg_relocate(tm, stats::runif(tm@npsi, lim[1L],
                                                         lim[2L])), yn)
          })
      }, error = function(e) NULL)
      if (is.null(cand)) next
      v <- modelterms7::seg_profile_rss(cand, yn)
      if (v < base - 1e-6 * (base + 1)) {
        gain <- TRUE
        st$terms[[r$param]][[r$term]] <- cand
        cols <- design[[r$param]]$blocks[[r$term]]
        cf <- obj_split(par)
        cf[[r$param]][cols] <- as.numeric(cand@blueprint$coef)
        par <- obj_stack(cf)
      }
    }
    st$key <- NULL
    st$value <- NULL
    if (gain) par else NULL
  }

  dry <- 0L
  kinds <- c("sweep", rep(c("boot", "random"), length.out = max(0L, nb - 1L)))
  for (b in seq_along(kinds)) {
    # four consecutive proposals the screen turns away end the loop: after
    # the sweep has landed the optimum, the remaining draws keep polishing
    # back onto it
    if (dry >= 4L) break
    restore(keep)
    start <- tryCatch(propose(kinds[b]), error = function(e) NULL)
    if (is.null(start)) {
      dry <- dry + 1L
      next
    }
    ro <- tryCatch(
      statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                        start, expected, approx, maxit, tol, vbq,
                        working_budget = 200L),
      error = function(e) NULL)
    if (!is.null(ro) && is.finite(ro$value) &&
        ro$value < best$value - 1e-8 * (abs(best$value) + 1)) {
      best[fields] <- ro[fields]
      keep <- snap()
      dry <- 0L
      if (vb$outer || vb$blocks) {
        vb_say("restart %d (%s) improved the objective to %.6f",
               b, kinds[b], ro$value)
      }
    } else {
      dry <- dry + 1L
    }
  }
  restore(keep)
  best
}
