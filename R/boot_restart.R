#' @include statmod.R
NULL

#' How Many Restarts the Terms Ask For
#'
#' @description
#' The largest `n_boot` any break-point term of the specification
#' declares, or zero. The value is declared on [modelterms7::seg()],
#' [modelterms7::jump()] and [modelterms7::jseg()], the
#' terms whose objective has the spurious local optima the device exists
#' for; running the restarts is this layer's, the way a penalty's
#' hyperparameters are declared by a term and estimated here.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A single non-negative integer.
#'
#' @seealso [statmod_boot_restart()]
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
#' Bootstrap restarting (Wood 2001) with the observation that makes it
#' cheap: the non-convexity of a break-point model lives entirely in the
#' positions, everything else being convex once they are held, so a
#' proposal is a CONFIGURATION OF POSITIONS and two configurations are
#' compared on the exact profile -- least squares at fixed positions, one
#' linear fit each -- rather than by refitting the model. Only a proposal
#' the profile prefers earns a refit, and the refit's answer is accepted on
#' the true objective.
#'
#' @details
#' Three proposal kinds, in order. The first is deterministic: each
#' break-point swept over the profile with the others held
#' ([modelterms7::seg_polish()]), which walks straight to a
#' feature the iteration pressed a break-point away from. The stochastic
#' ones alternate a BOOTSTRAP sweep -- the same descent on the profile of a
#' resample, the multinomial counts entering as weights, which moves the
#' profile's optima the way refitting the resample would -- with a sweep
#' from RANDOM positions drawn over the confinement interval. The screen
#' compares every proposal's unweighted profile against the incumbent's,
#' so a proposal that lands back on the incumbent costs a few grid-many
#' linear fits and no refit: measured at n = 10000, a dry proposal fell
#' from a whole refit (5 to 15 s) to about half a second, and the earlier
#' design that refitted every proposal had spent 945 s re-verifying an
#' optimum the sweep had already found. Four consecutive dry proposals end
#' the loop.
#'
#' The profile reads the response net of the other contributions in the
#' term's equation, on the predictor scale: exact for an identity link, a
#' proposal elsewhere -- the argument [modelterms7::seg_start()]
#' already makes. The refreshable and structural state of the design is
#' snapshotted at the incumbent and restored whenever a candidate loses,
#' and the draws come from the session's generator, so a fit with restarts
#' is reproducible under `set.seed()`.
#'
#' @param spec The specification.
#' @param design The design.
#' @param blocks The block split.
#' @param hyper The hyperparameters the fit ended at.
#' @param inner_optimizer How the smooth block is fitted.
#' @param res The fitted result the restarts try to improve.
#' @param expected,approx,maxit,tol As in [statmod_alternate()].
#' @param vb The resolved verbosity.
#' @param nb At most how many proposals.
#'
#' @return `res`, with `par`, `value`, `converged`,
#'   `obj` and the block histories replaced when a restart improved
#'   the fit; any other field -- an outer search's history, its optimizer
#'   -- is kept.
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
