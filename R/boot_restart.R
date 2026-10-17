#' @include statmod.R
NULL

#' Restart Budget of the Break-Point Terms
#'
#' @description
#' Returns the largest `n_boot` set in a break-point term that is written
#' directly in an equation of the specification, and zero when there is no
#' such term. A break-point term nested inside a parameter of
#' [modelterms7::nl()] is not counted. [statmod_boot_restart()] takes this
#' number as the maximum number of proposals that it tries.
#'
#' @details
#' The number is set in the term. [modelterms7::seg()], [modelterms7::jump()]
#' and [modelterms7::jseg()] each take `n_boot`, with a default of 10,
#' because the objective of these terms has spurious local optima, and the
#' restarts search for a better one. The restarts are run by this package, which
#' can refit the model.
#'
#' A model with two break-point terms that set `n_boot` to 10 and to 25 has a
#' budget of 25. The restart loop works on all of the terms together, since a
#' proposal is a configuration of positions for every break-point term in the
#' model, so the budget is a single number and the largest of the set values
#' is used.
#'
#' @param spec A [StatmodSpec()]. The terms of every equation are examined.
#'
#' @return A single non-negative integer. It is zero when no term sets
#'   `n_boot`, and [statmod()] does not run the restart loop for a budget of
#'   zero.
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


#' Bootstrap Restarting of a Break-Point Fit
#'
#' @description
#' Improves a fitted break-point model by bootstrap restarting (Wood 2001),
#' screening each proposal on the profile of the term's equation before a
#' refit of the model. The objective of a model with break-points has
#' spurious local minima, and an ordinary fit converges into the one in which
#' its starting positions lie; this search looks for a better one.
#'
#' @details
#' # The screen on the profile
#'
#' The non-convexity of a break-point model is in the positions. Given the
#' positions, the coefficients of a term enter linearly, so the profile at a
#' fixed configuration of positions is one weighted linear fit of the working
#' response of the term's equation. A proposal is therefore a configuration of
#' positions, and two proposals are ranked by their profiles at the cost of
#' one linear fit each. Only a proposal that ranks higher on the profile leads
#' to a refit of the whole model, and the refit is accepted or rejected on the
#' true objective. A proposal that the profile rejects costs profile evaluations
#' only, with no refit of the model.
#'
#' # The three proposal kinds
#'
#' The kinds are tried in this order.
#'
#' 1. The deterministic sweep. Each break-point in turn is swept over a
#'    grid on the profile while the others are held, as
#'    [modelterms7::seg_polish()] does. This reaches a feature of the profile
#'    from which the fitting iteration had moved a break-point away.
#' 2. The bootstrap sweep. The same descent on the profile of a
#'    resample, in which the multinomial counts enter as weights. This moves
#'    the profile's optima in the way that refitting the resample would.
#' 3. The random sweep. The same descent from positions drawn uniformly
#'    over the confinement interval.
#'
#' The two stochastic kinds alternate. Four consecutive proposals that do not
#' improve the objective, either because the screen rejects them or because
#' the refit does not decrease the objective, end the loop, whatever budget
#' is left.
#'
#' # The profile for a gaussian location and for other equations
#'
#' The profile is the weighted least-squares fit of the working response of
#' the term's equation, net of the other contributions to that equation, with
#' the working weights of a scoring step. For the location parameter of a
#' gaussian family on the identity link, this is the least-squares objective
#' of the model. In every other equation, it is a quadratic model of the
#' objective that ranks the candidates, and each candidate is accepted or
#' rejected on the true objective. [modelterms7::seg_start()] uses the same
#' profile to place starting values.
#'
#' # Reproducibility
#'
#' The draws come from the session's generator, so a fit with restarts
#' repeats under [set.seed()]. The refreshable and structural state of the
#' design is saved at the incumbent and restored whenever a candidate
#' is rejected, so that a rejected proposal leaves the state unchanged.
#'
#' @param spec The [StatmodSpec()] being fitted.
#' @param design The assembled design, as [statmod_design()] returns it.
#' @param blocks The split of the terms into the jointly fitted smooth block
#'   and the kinked ones, as [statmod_blocks()] returns it.
#' @param hyper The hyperparameters at which the fit ended, held fixed
#'   throughout. The restarts search over positions only.
#' @param inner_optimizer How the smooth block is fitted, [iwls()] or an
#'   \pkg{optimizers7} optimizer.
#' @param res The fitted result to improve, as [statmod_alternate()] returns
#'   it. It is returned unchanged when no proposal improves the objective, and
#'   also when the design carries no state, or when the response is a matrix,
#'   is not numeric or has missing values.
#' @param expected,approx,maxit,tol Passed to [statmod_alternate()] for each
#'   refit, with the meanings that they have there.
#' @param vb The resolved verbosity, as [verbosity()] returns it.
#' @param nb The budget: the maximum number of proposals to try, as
#'   [seg_boot_total()] reports it. The first proposal is the deterministic
#'   sweep. A budget below one returns `res` unchanged.
#'
#' @return `res`, in which `par`, `value`, `converged`, `obj` and the block
#'   histories are replaced when a restart improves the objective. Every other
#'   element is carried over unchanged, including the history of an outer
#'   search and its optimizer, so the result has the same structure either
#'   way.
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
    if (is.null(wk)) stop("there is no working response")
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

  if (nb < 1L) return(res)
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
#' or [modelterms7::jseg()] term is moved to the minimum of the least-squares
#' profile of its working response ([modelterms7::seg_polish_exact()]), which
#' is exact only for the location of a gaussian family, held there
#' ([modelterms7::seg_hold()]), and the coefficients of the model are refitted
#' while the positions are held.
#'
#' @details
#' The working construction of Fasola, Muggeo and Kuchenhoff stops at a fixed
#' point of its own iteration. Its coefficients there are those of the working
#' model, whose block carries a weight frozen at the previous iterate, and not
#' the coefficients that maximize the likelihood at the positions reached; and
#' the profile of a discontinuous term is constant between consecutive
#' observations, so the iteration can stop one interval away from the
#' minimum. The certificate of such a fit can then read the mode along the
#' auxiliary column of the working block, which is not a direction of the
#' model.
#'
#' Held, a break-point term contributes a linear function of its remaining
#' coefficients, so the refit is an ordinary fit, and the coordinate that
#' carried the position is held in the solve and left out of the information
#' ([term_held_stack()]). `segmented` and `stepmented` report their
#' coefficients in this way, at the estimated positions and conditional on
#' them.
#'
#' The profile is weighted least squares of the working response of the
#' term's own equation, net of the rest of that equation, on the term's own
#' columns, with the working weights of a scoring step. It is exact for the
#' location of a gaussian family on the identity link and is a quadratic model
#' of the objective in any other equation. A polished position is kept only
#' where the refit at it reaches a better objective than the refit at the
#' position that the iteration reached. The polish is repeated, up to
#' `rounds` times, while the refit improves the objective.
#'
#' A term whose break-point or per-break-point coefficients carry a development
#' is held and polished like any other term, with the exact line searches of
#' [modelterms7::seg_polish_exact()] in the coefficients of the development.
#' In an equation whose working profile is only a model of the objective, the
#' five deepest local minima of the profile are also refitted while the
#' positions are held, and so are the neighboring intervals of the best
#' position; these two searches skip a developed term, whose positions are one
#' per observation.
#'
#' @param spec The specification, including the holds.
#' @param design The design, carrying the terms in its state.
#' @param blocks,hyper,inner_optimizer,expected,approx,maxit,tol,vb As
#'   [statmod_alternate()].
#' @param res The fit to settle, as [statmod_alternate()] returns it.
#' @param rounds The maximum number of times that the polished positions are
#'   read again after a refit, 5 by default.
#'
#' @return A list: `res`, the settled fit, `moved`, `TRUE` where a
#'   polished position replaced the one that the iteration reached, and
#'   `settled`, `TRUE` where a term was held.
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
#' After the working phase and the restarts, every [modelterms7::seg()] term
#' that is not smoothed, not developed and not marginal is moved to the
#' minimum of the least-squares profile of its working response over every
#' interval of the covariate ([modelterms7::seg_polish_exact()]), and the
#' model is refitted from there while the positions are free.
#'
#' @details
#' The profile of a change of slope is smooth inside each interval between
#' consecutive observed values and has a kink at each of them, so the
#' iteration of Muggeo (2003) can stop on a kink that is not a minimum.
#' Inside an interval the profile has one stationary
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
#'   position replaced the one that the iteration reached.
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
#' identity link, whose log-likelihood is quadratic in that predictor. In
#' every other equation the profile is a quadratic model of the objective.
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
