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

  # the working response of the term's own equation net of the rest of that
  # equation, with its working weights: for a gaussian mean the response
  # itself, elsewhere the quadratic model of the objective a scoring step
  # reads. It was the RESPONSE net of the equation whatever the equation, so
  # a break-point in sigma's equation was proposed on a quantity that has
  # nothing to do with sigma: measured on MASS::mcycle with a jseg there,
  # every proposal passed the screen and was refitted, and the fit took
  # 2961 s.
  net_y <- function(r, tm) {
    cfl <- obj_split(best$par)
    ep <- statmod_eta(spec, design, cfl)
    wk <- coord_working(spec, ep, cfl, design, r$param, TRUE, approx)
    if (is.null(wk)) stop("no working response")
    eta <- rep_len(if (is.null(ep$eta_static)) ep$eta[[r$param]] else
      ep$eta_static[[r$param]], n)
    list(y = wk$z - (eta - as.numeric(modelterms7::term_value(tm))),
         w = wk$w)
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
        nw <- net_y(r, tm)
        yn <- nw$y
        base <- modelterms7::seg_profile_rss(tm, yn, weights = nw$w)
        switch(kind,
          sweep = modelterms7::seg_polish(tm, yn, weights = nw$w),
          boot = {
            w <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
            modelterms7::seg_polish(tm, yn, weights = nw$w * w)
          },
          random = {
            lim <- tm@blueprint$lim
            modelterms7::seg_polish(
              modelterms7::seg_relocate(tm, stats::runif(tm@npsi, lim[1L],
                                                         lim[2L])), yn,
              weights = nw$w)
          })
      }, error = function(e) NULL)
      if (is.null(cand)) next
      v <- modelterms7::seg_profile_rss(cand, yn, weights = nw$w)
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


#' Settle the Sharp Break-Points on the Profile and Hold Them
#'
#' @description
#' After the working phase and the restarts, every sharp [modelterms7::jump()]
#' or [modelterms7::jseg()] term is moved to the minimum of its exact
#' least-squares profile ([modelterms7::seg_polish_exact()]), held there
#' ([modelterms7::seg_hold()]), and the coefficients of the model are refitted
#' with the positions held.
#'
#' @details
#' The working construction of Fasola, Muggeo and Kuchenhoff stops at a fixed
#' point of its own iteration. Its coefficients there are those of the working
#' model, whose block carries a weight frozen at the previous iterate, and not
#' the coefficients that maximize the likelihood at the positions reached; and
#' the profile of a discontinuous term is constant between consecutive
#' observations, so the iteration can stop one interval away from the
#' minimum. Measured on 16 samples of 200 observations (a jseg and a jump, 8
#' seeds each), the fit reached the interval of the profile's minimum in 11,
#' its coefficients sat up to 0.04 log-likelihood units below least squares
#' at the position reached, and `statmod_certificate()` reported
#' `not converged` in 12, reading the mode along the working block's auxiliary
#' column, which is not a direction of the model.
#'
#' Held, a break-point term contributes a linear function of its remaining
#' coefficients, so the refit is an ordinary fit, and the coordinate that
#' carried the position is held in the solve and left out of the information
#' ([term_held_stack()]). This is how `segmented` and `stepmented` report
#' their coefficients: at the estimated positions, conditional on them.
#'
#' The profile is weighted least squares of the working response of the
#' term's own equation, net of the rest of that equation, on the term's own
#' columns, with the working weights of a scoring step: exact for a gaussian
#' mean and the quadratic model of the objective for any other equation. A polished position is kept
#' only where the refit at it reaches a better objective than the refit at the
#' position the iteration reached. A term whose per-break-point coefficients
#' carry a development is held without polishing, its positions being one per
#' observation.
#'
#' @param spec The specification, holds included.
#' @param design The design, carrying the terms in its state.
#' @param blocks,hyper,inner_optimizer,expected,approx,maxit,tol,vb As
#'   [statmod_alternate()].
#' @param res The fit to settle, as [statmod_alternate()] returns it.
#'
#' @return A list: `res`, the settled fit, `moved`, `TRUE` where a
#'   polished position replaced the one the iteration reached, and `settled`,
#'   `TRUE` where a term was held.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change models.
#' *Computational Statistics*, 33, 997--1015.
#'
#' @seealso [statmod_boot_restart()], which runs before it.
#'
#' @keywords internal
statmod_settle_breakpoints <- function(spec, design, blocks, hyper,
                                       inner_optimizer, res, expected, approx,
                                       maxit, tol, vb, rounds = 5L) {
  st <- attr(design, "state")
  rf <- attr(design, "refresh")
  none <- list(res = res, moved = FALSE, settled = FALSE)
  if (is.null(st) || !length(rf)) return(none)
  sharp <- Filter(function(r) {
    tm <- st$terms[[r$param]][[r$term]]
    S7::S7_inherits(tm, modelterms7::SegTerm) &&
      tm@kind %in% c("jump", "jseg") && is.null(tm@spec$smoothed)
  }, rf)
  if (!length(sharp)) return(none)
  obj_split <- res$obj$split
  obj_stack <- res$obj$stack
  vbq <- lapply(vb, function(x) FALSE)
  y0 <- spec@response
  can_polish <- is.numeric(y0) && !is.matrix(y0) &&
    length(y0) == spec@n_obs && !anyNA(y0)
  reset <- function() {
    st$key <- NULL
    st$value <- NULL
  }
  # the terms with the given positions written in and every sharp term held,
  # and the coefficients that go with them
  place <- function(terms, par) {
    cf <- obj_split(par)
    for (r in sharp) {
      tm <- modelterms7::seg_hold(terms[[r$param]][[r$term]])
      st$terms[[r$param]][[r$term]] <- tm
      cols <- design[[r$param]]$blocks[[r$term]]
      cf[[r$param]][cols] <- as.numeric(tm@blueprint$coef)
    }
    reset()
    obj_stack(cf)
  }
  refit <- function(par) {
    tryCatch(statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                               par, expected, approx, maxit, tol, vbq),
             error = function(e) NULL)
  }
  ok_fit <- function(r) !is.null(r) && is.finite(r$value)
  keep <- st$terms

  # FIRST the hold where the iteration stopped and the exact refit there:
  # the polish reads the working response of the term's equation, and at the
  # working fixed point the coefficients beside the break-point are those of
  # the working model. Measured on MASS::mcycle with a jseg in sigma's
  # equation, the polish read at the fixed point put the minimum at the
  # confinement limit where the iteration had stopped, and read at the exact
  # fit there it put it at 16.3, beside the profile's minimum at 14.7.
  best <- refit(place(keep, res$par))
  if (!ok_fit(best)) {
    st$terms <- keep
    reset()
    return(none)
  }
  best_terms <- st$terms
  moved <- FALSE
  # THEN the polish, repeated while it improves the objective: the working
  # response moves with the coefficients the refit changes, so one reading
  # is a step and not the answer
  for (it in seq_len(if (can_polish) as.integer(rounds) else 0L)) {
    cfl <- obj_split(best$par)
    ep <- tryCatch(statmod_eta(spec, design, cfl), error = function(e) NULL)
    if (is.null(ep)) break
    cand <- best_terms
    any_move <- FALSE
    for (r in sharp) {
      tm <- best_terms[[r$param]][[r$term]]
      p <- r$param
      # the working response of the term's own equation, net of the rest of
      # that equation, with the working weights: for a gaussian mean it is the
      # response itself, and in any other equation it is the quadratic model
      # of the objective a scoring step reads
      wk <- tryCatch(coord_working(spec, ep, cfl, design, p, TRUE, approx),
                     error = function(e) NULL)
      if (is.null(wk)) next
      eta <- rep_len(if (is.null(ep$eta_static)) ep$eta[[p]] else
        ep$eta_static[[p]], spec@n_obs)
      yn <- wk$z - (eta - as.numeric(modelterms7::term_value(tm)))
      pt <- tryCatch(modelterms7::seg_polish_exact(tm, yn, weights = wk$w),
                     error = function(e) NULL)
      if (is.null(pt)) next
      if (!isTRUE(all.equal(modelterms7::seg_psi(pt), modelterms7::seg_psi(tm),
                            tolerance = 0))) {
        cand[[p]][[r$term]] <- pt
        any_move <- TRUE
      }
    }
    if (!any_move) break
    st$terms <- best_terms
    r1 <- refit(place(cand, best$par))
    if (ok_fit(r1) && r1$value < best$value - 1e-8 * (abs(best$value) + 1)) {
      best <- r1
      best_terms <- st$terms
      moved <- TRUE
    } else {
      break
    }
  }
  # AND WHERE THE PROFILE IS ONLY A MODEL OF THE OBJECTIVE, a local search
  # on the objective itself. The working response is exact for a gaussian
  # mean and a quadratic model anywhere else, so there the polish stops near
  # the minimum and not at it: measured on MASS::mcycle with a jseg in sigma's
  # equation, it stopped at 16.3, a local minimum of the objective six
  # intervals from the global one at 14.7 (656.08 against 657.40 at the fit's
  # smoothing parameter). The five deepest local minima of the working
  # profile are refitted with the positions held and compared on the
  # objective, and then the neighbouring intervals of the best, the window
  # recentred while the best sits at its edge.
  for (r in sharp) {
    if (!can_polish || working_exact(spec, r$param)) next
    tm0 <- best_terms[[r$param]][[r$term]]
    if (tm0@blueprint$developed) next
    xv <- tm0@blueprint$xv
    lim <- tm0@blueprint$lim
    u <- sort(unique(xv))
    mid <- (u[-1L] + u[-length(u)]) / 2
    mid <- mid[mid > lim[1L] & mid < lim[2L]]
    if (length(mid) < 2L) next
    for (k in seq_len(tm0@npsi)) {
      # the candidates: the five deepest local minima of the working profile
      # over every interval, read at the exact fit reached so far, each
      # refitted and compared on the objective
      cfl <- obj_split(best$par)
      ep <- tryCatch(statmod_eta(spec, design, cfl), error = function(e) NULL)
      wk <- if (is.null(ep)) NULL else
        tryCatch(coord_working(spec, ep, cfl, design, r$param, TRUE, approx),
                 error = function(e) NULL)
      if (!is.null(wk)) {
        tmb <- best_terms[[r$param]][[r$term]]
        eta <- rep_len(if (is.null(ep$eta_static)) ep$eta[[r$param]] else
          ep$eta_static[[r$param]], spec@n_obs)
        yn <- wk$z - (eta - as.numeric(modelterms7::term_value(tmb)))
        pr <- tryCatch(modelterms7::seg_profile_intervals(tmb, yn, k,
                                                          weights = wk$w),
                       error = function(e) NULL)
        if (!is.null(pr) && nrow(pr) >= 2L) {
          v <- pr$rss
          nv <- length(v)
          lo <- c(Inf, v[-nv])
          hi <- c(v[-1L], Inf)
          locmin <- which(is.finite(v) & v <= lo & v <= hi)
          cands <- pr$psi[locmin[order(v[locmin])]][seq_len(min(5L,
                                                               length(locmin)))]
          psi_b <- as.numeric(modelterms7::seg_psi(tmb))
          for (q in cands) {
            if (isTRUE(abs(q - psi_b[k]) < 1e-12)) next
            pj <- psi_b
            pj[k] <- q
            cj <- best_terms
            cj[[r$param]][[r$term]] <- tryCatch(
              modelterms7::seg_relocate(tmb, pj), error = function(e) NULL)
            if (is.null(cj[[r$param]][[r$term]])) next
            st$terms <- best_terms
            r1 <- refit(place(cj, best$par))
            if (ok_fit(r1) &&
                r1$value < best$value - 1e-8 * (abs(best$value) + 1)) {
              best <- r1
              best_terms <- st$terms
              moved <- TRUE
            }
          }
          st$terms <- best_terms
          reset()
        }
      }
      for (step in seq_len(10L)) {
        tmb <- best_terms[[r$param]][[r$term]]
        psi <- as.numeric(modelterms7::seg_psi(tmb))
        i0 <- which.min(abs(mid - psi[k]))
        win <- setdiff(max(1L, i0 - 3L):min(length(mid), i0 + 3L), i0)
        found <- 0L
        for (j in win) {
          pj <- psi
          pj[k] <- mid[j]
          cj <- best_terms
          cj[[r$param]][[r$term]] <- tryCatch(
            modelterms7::seg_relocate(tmb, pj), error = function(e) NULL)
          if (is.null(cj[[r$param]][[r$term]])) next
          st$terms <- best_terms
          r1 <- refit(place(cj, best$par))
          if (ok_fit(r1) &&
              r1$value < best$value - 1e-8 * (abs(best$value) + 1)) {
            best <- r1
            best_terms <- st$terms
            moved <- TRUE
            found <- j
          }
        }
        st$terms <- best_terms
        reset()
        # recentre only where the best sits at the edge of the window
        if (!found || abs(found - i0) < 3L) break
      }
    }
  }
  st$terms <- best_terms
  reset()
  fields <- c("par", "value", "converged", "obj", "aliased", "hist_blocks",
              "hist_inner")
  res[fields] <- best[fields]
  list(res = res, moved = moved, settled = TRUE)
}


#' Settle the Changes of Slope on the Profile
#'
#' @description
#' After the working phase and the restarts, every sharp [modelterms7::seg()]
#' term is moved to the minimum of its exact least-squares profile over every
#' interval of the covariate ([modelterms7::seg_polish_exact()]), and the
#' model is refitted from there with the positions free.
#'
#' @details
#' The profile of a change of slope is smooth inside each interval between
#' consecutive observed values and has a kink at each of them, so the
#' iteration of Muggeo (2003) can stop on a kink that is not a minimum:
#' measured on `segmented::globTempAnom` with four break-points, it stopped
#' with two positions on observed years and the inner fit 0.106
#' log-likelihood units above its mode, and the REML criterion could not be
#' evaluated at its start. Inside an interval the profile has one stationary
#' point, so its minimum over every interval is exact
#' ([modelterms7::seg_polish_exact()]). The profile is read on the working
#' response of the term's equation, and a polished position is kept only
#' where the refit from it reaches a better objective. In an equation whose
#' working profile only approximates the objective, the five deepest local
#' minima of the profile are refitted as well.
#'
#' @inheritParams statmod_settle_breakpoints
#'
#' @return A list: `res`, the fit, and `moved`, `TRUE` where a polished
#'   position replaced the one the iteration reached.
#'
#' @references
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. *Statistics in Medicine*, 22, 3055--3071.
#'
#' @keywords internal
statmod_settle_seg <- function(spec, design, blocks, hyper, inner_optimizer,
                               res, expected, approx, maxit, tol, vb,
                               rounds = 5L) {
  st <- attr(design, "state")
  rf <- attr(design, "refresh")
  none <- list(res = res, moved = FALSE)
  if (is.null(st) || !length(rf)) return(none)
  y0 <- spec@response
  if (!is.numeric(y0) || is.matrix(y0) || length(y0) != spec@n_obs ||
      anyNA(y0)) {
    return(none)
  }
  segs <- Filter(function(r) {
    tm <- st$terms[[r$param]][[r$term]]
    S7::S7_inherits(tm, modelterms7::SegTerm) && identical(tm@kind, "seg") &&
      is.null(tm@spec$smoothed) && !isTRUE(tm@blueprint$developed) &&
      !isTRUE(tm@spec$marginal)
  }, rf)
  if (!length(segs)) return(none)
  obj_split <- res$obj$split
  obj_stack <- res$obj$stack
  vbq <- lapply(vb, function(x) FALSE)
  reset <- function() {
    st$key <- NULL
    st$value <- NULL
  }
  # the terms at the given positions and the coefficients that go with them
  place <- function(terms, par) {
    cf <- obj_split(par)
    for (r in segs) {
      tm <- terms[[r$param]][[r$term]]
      st$terms[[r$param]][[r$term]] <- tm
      cols <- design[[r$param]]$blocks[[r$term]]
      cf[[r$param]][cols] <- as.numeric(tm@blueprint$coef)
    }
    reset()
    obj_stack(cf)
  }
  refit <- function(par) {
    tryCatch(statmod_alternate(spec, design, blocks, hyper, inner_optimizer,
                               par, expected, approx, maxit, tol, vbq),
             error = function(e) NULL)
  }
  better <- function(r1, r0) {
    !is.null(r1) && is.finite(r1$value) &&
      (!is.finite(r0$value) || r1$value < r0$value - 1e-8 * (abs(r0$value) + 1))
  }
  # the working response of an equation net of one term, and its weights
  working_of <- function(par, r, tm) {
    cfl <- obj_split(par)
    ep <- tryCatch(statmod_eta(spec, design, cfl), error = function(e) NULL)
    if (is.null(ep)) return(NULL)
    wk <- tryCatch(coord_working(spec, ep, cfl, design, r$param, TRUE, approx),
                   error = function(e) NULL)
    if (is.null(wk)) return(NULL)
    eta <- rep_len(if (is.null(ep$eta_static)) ep$eta[[r$param]] else
      ep$eta_static[[r$param]], spec@n_obs)
    list(y = wk$z - (eta - as.numeric(modelterms7::term_value(tm))), w = wk$w)
  }
  best <- res
  best_terms <- st$terms
  moved <- FALSE
  for (it in seq_len(as.integer(rounds))) {
    cand <- best_terms
    any_move <- FALSE
    for (r in segs) {
      tm <- best_terms[[r$param]][[r$term]]
      wk <- working_of(best$par, r, tm)
      if (is.null(wk)) next
      pt <- tryCatch(modelterms7::seg_polish_exact(tm, wk$y, weights = wk$w),
                     error = function(e) NULL)
      if (is.null(pt)) next
      if (!isTRUE(all.equal(modelterms7::seg_psi(pt), modelterms7::seg_psi(tm),
                            tolerance = 1e-10))) {
        cand[[r$param]][[r$term]] <- pt
        any_move <- TRUE
      }
    }
    if (!any_move) break
    st$terms <- best_terms
    r1 <- refit(place(cand, best$par))
    if (better(r1, best)) {
      best <- r1
      best_terms <- st$terms
      moved <- TRUE
    } else {
      break
    }
  }
  # where the working profile is a model of the objective, the deepest local
  # minima of the profile are refitted and compared on the objective
  for (r in segs) {
    if (working_exact(spec, r$param)) next
    tm0 <- best_terms[[r$param]][[r$term]]
    for (k in seq_len(tm0@npsi)) {
      tmb <- best_terms[[r$param]][[r$term]]
      wk <- working_of(best$par, r, tmb)
      if (is.null(wk)) next
      pr <- tryCatch(modelterms7::seg_profile_intervals(tmb, wk$y, k,
                                                        weights = wk$w),
                     error = function(e) NULL)
      if (is.null(pr) || nrow(pr) < 2L) next
      v <- pr$rss
      nv <- length(v)
      locmin <- which(is.finite(v) & v <= c(Inf, v[-nv]) & v <= c(v[-1L], Inf))
      cands <- pr$psi[locmin[order(v[locmin])]][seq_len(min(5L,
                                                           length(locmin)))]
      psi_b <- as.numeric(modelterms7::seg_psi(tmb))
      for (q in cands) {
        if (isTRUE(abs(q - psi_b[k]) < 1e-10)) next
        pj <- psi_b
        pj[k] <- q
        cj <- best_terms
        cj[[r$param]][[r$term]] <- tryCatch(
          modelterms7::seg_relocate(tmb, pj), error = function(e) NULL)
        if (is.null(cj[[r$param]][[r$term]])) next
        st$terms <- best_terms
        r1 <- refit(place(cj, best$par))
        if (better(r1, best)) {
          best <- r1
          best_terms <- st$terms
          moved <- TRUE
        }
      }
    }
  }
  st$terms <- best_terms
  reset()
  fields <- c("par", "value", "converged", "obj", "aliased", "hist_blocks",
              "hist_inner")
  res[intersect(fields, names(best))] <- best[intersect(fields, names(best))]
  list(res = res, moved = moved)
}


#' Is the Working Profile of an Equation Its Objective?
#'
#' @description
#' `TRUE` where the working response of a distribution parameter's equation
#' gives the objective exactly: the location of a gaussian family on the
#' identity link, whose log-likelihood is quadratic in that predictor. Every
#' other equation reads a quadratic model of its objective.
#'
#' @param spec A [StatmodSpec()].
#' @param p The distribution parameter, a string.
#'
#' @return A single logical.
#'
#' @keywords internal
working_exact <- function(spec, p) {
  d <- spec@distrib
  isTRUE(identical(p, d@params[1L]) &&
           grepl("^gaussian", d@distrib_name) &&
           identical(d@link_params[[1L]]@link_name, "identity"))
}
