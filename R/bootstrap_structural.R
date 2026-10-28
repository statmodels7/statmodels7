#' @include bootstrap_engine.R
NULL

# The bootstrap of a model with a structural term (gas(), regime()) or with a
# term whose block moves with its coefficients (seg(), jump(), jseg(), nl()).
# The engine of bootstrap_engine.R calls these functions at four points: the
# plan of a replicate, its refit and validation, the names of the coordinates,
# and the covariance of the stored replicates.


#' The Kind of Model That the Bootstrap Treats Apart
#'
#' @description
#' Returns `"structural"` for a model with a [modelterms7::gas()] or
#' [modelterms7::regime()] term, `"refresh"` for a model with a term whose
#' block moves with its coefficients, and `NULL` for every other model.
#'
#' @param design0 The design of the fit.
#'
#' @return A string or `NULL`.
#'
#' @keywords internal
boot_mode <- function(design0) {
  if (length(attr(design0, "structural"))) return("structural")
  if (length(attr(design0, "refresh"))) return("refresh")
  NULL
}


#' The Draw of a Replicate of a Model With a Structural Term
#'
#' @description
#' Returns a function of no arguments that draws a new path of the series from
#' the fitted model.
#'
#' @details
#' The path comes from [rstatmod_eta()], that is from
#' [modelterms7::term_simulate()] run at the fitted coefficients and the
#' fitted parameters of the term. A score-driven term draws the response as
#' it runs, and the level at one time is driven by the score of the value just
#' drawn, so the replicate is a new realization of the process and not the
#' fitted path with new noise. A latent chain draws its states. The effects of
#' a [modelterms7::random()] term keep their estimates.
#'
#' Only the parametric scheme applies: a resample of rows or of groups breaks
#' the recursion of the filter at every join.
#'
#' @param object The [StatmodFit()].
#' @param resample The [resample_method].
#' @param spec0,design0 Its specification and design.
#'
#' @return A function that returns `list(y = )`, with the attributes `kind`,
#'   `targets`, `random` and `mode` that [bootstrap_engine()] reads.
#'
#' @keywords internal
boot_plan_structural <- function(object, resample, spec0, design0) {
  cf <- object@coefficients
  out <- function() {
    ep <- rstatmod_eta(spec0, design0, cf)
    y <- ep$y
    if (is.null(y)) {
      y <- distributions7::distrib_rng(spec0@distrib, spec0@n_obs, ep$theta)
    }
    list(y = as.numeric(y))
  }
  attr(out, "kind") <- "parametric"
  attr(out, "targets") <- list()
  attr(out, "random") <- "none"
  attr(out, "mode") <- "structural"
  out
}


#' Add the Refit Seed to the Draw of a Break-Point Model
#'
#' @description
#' Wraps a draw so that every replicate carries an integer seed, drawn from
#' the stream of the session, and tags the draw with its mode.
#'
#' @details
#' The restart cycle of a break-point fit ([statmod_boot_restart()]) draws
#' from the generator. The refit of a replicate may run in a worker process,
#' so it is seeded from the replicate's own draw and not from the state of
#' whichever process runs it. The result then depends on the stream of the
#' session and not on the number of workers.
#'
#' @param draw The function from [bootstrap_plan()].
#' @param design0 The design of the fit.
#'
#' @return `draw`, wrapped for a break-point model and unchanged otherwise.
#'
#' @keywords internal
boot_plan_wrap <- function(draw, design0) {
  if (!identical(boot_mode(design0), "refresh")) return(draw)
  out <- function() {
    pl <- draw()
    pl$seed <- sample.int(.Machine$integer.max, 1L)
    pl
  }
  for (a in c("kind", "targets", "random")) attr(out, a) <- attr(draw, a)
  attr(out, "mode") <- "refresh"
  out
}


#' The Names of the Coordinates of the Stored Replicates
#'
#' @description
#' Returns the names of the coefficients of every equation followed by the
#' free parameters of a structural term, the joint order of
#' [vcov.StatmodFit()].
#'
#' @param spec0,design0 The specification and design of the fit.
#'
#' @return A character vector.
#'
#' @keywords internal
boot_coord_names <- function(spec0, design0) {
  nm <- rownames(coef_labels(spec0, design0))
  c(nm, unname(structural_tail_names(spec0, design0)))
}


#' The Free Parameters of a Structural Term in a Refit
#'
#' @description
#' Returns the unconstrained values of the free parameters of every structural
#' term, in the order of the joint vector.
#'
#' @param spec,design The specification and the design that the refit ran on.
#'
#' @return A numeric vector.
#'
#' @keywords internal
boot_tail_values <- function(spec, design) {
  su <- attr(design, "structural")
  sst <- statmod_structural_state(design)
  out <- numeric(0)
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    free <- setdiff(modelterms7::term_params(tm), sst$held[[u$term]])
    out <- c(out, unlist(sst$zeta[[u$term]])[free])
  }
  unname(out)
}


#' Refit One Replicate of a Model With a Break-Point Term
#'
#' @description
#' Refits the model on the replicate from the break-points and coefficients of
#' the fit, and then runs the restart cycle of [statmod_boot_restart()] and
#' the settling of the sharp break-points, as [statmod()] does at the end of
#' a fit.
#'
#' @details
#' The refit from the fit's positions is [bootstrap_refit()]. A start from the
#' positions of the fit alone would find the local optimum nearest to them,
#' which is another estimator than the one that produced the fit; the restart
#' cycle searches for a better optimum of the replicate with the budget `n_boot`
#' that the terms set, and the settling puts a sharp break-point at the minimum
#' of its profile. The hyperparameters are held at the fit's values with
#' `refit = "coefficients"` and selected again with `"full"`, before the
#' restart.
#'
#' @param object The [StatmodFit()].
#' @param spec_b The specification of the replicate.
#' @param refit,method,cfg,beta0 As in [bootstrap_refit()].
#'
#' @return A list with `coefficients`, `hyper` and `design`.
#'
#' @keywords internal
boot_refit_breakpoint <- function(object, spec_b, refit, method, cfg, beta0) {
  r <- bootstrap_refit(object, spec_b, refit, method, cfg, beta0)
  design <- r$design
  blocks <- statmod_blocks(spec_b, design)
  res <- r$res
  hyper <- r$hyper
  vb <- verbosity(0)
  spec_h <- if (length(res$held_coef)) {
    S7::set_props(spec_b, held_coef = res$held_coef)
  } else spec_b
  ex <- cfg$expected
  em <- cfg$emethod
  nb <- seg_boot_total(spec_b)
  if (nb > 0L) {
    res <- statmod_boot_restart(spec_h, design, blocks, hyper, method, res, ex,
                                em, cfg$maxit, cfg$tol, vb, nb)
  }
  for (round in seq_len(2L)) {
    sg0 <- statmod_settle_seg(spec_h, design, blocks, hyper, method, res, ex,
                              em, cfg$maxit, cfg$tol, vb)
    res <- sg0$res
    st0 <- statmod_settle_breakpoints(spec_h, design, blocks, hyper, method,
                                      res, ex, em, cfg$maxit, cfg$tol, vb)
    res <- st0$res
    if (st0$settled) {
      rf2 <- attr(design, "refresh")
      stt <- attr(design, "state")
      for (i in seq_along(rf2)) {
        tm <- stt$terms[[rf2[[i]]$param]][[rf2[[i]]$term]]
        if (isTRUE(tryCatch(tm@build_info$held, error = function(e) FALSE))) {
          rf2[[i]]$frozen <- FALSE
        }
      }
      attr(design, "refresh") <- rf2
    }
    if (!st0$moved && !sg0$moved) break
  }
  cf <- boot_record_positions(spec_b, design, res$obj$split(res$par))
  list(coefficients = cf, hyper = hyper, design = design)
}


#' Store the Position of a Held Break-Point in Its Slot
#'
#' @description
#' A sharp [modelterms7::jump()] or [modelterms7::jseg()] is settled on the
#' profile and held: its slot in the coefficient vector is a working column
#' whose value is the same in every replicate, and the position of the
#' break-point is kept in the term. The function writes the position of each
#' held break-point into its slot, so that the replicates carry it.
#'
#' @param spec_b The specification of the replicate.
#' @param design The design that the refit ran on, whose state holds the terms
#'   as the refit left them.
#' @param cf The coefficients of the refit.
#'
#' @return The coefficients, with the held slots replaced.
#'
#' @keywords internal
boot_record_positions <- function(spec_b, design, cf) {
  st <- attr(design, "state")
  if (is.null(st)) return(cf)
  for (p in spec_b@distrib@params) {
    for (nm in names(spec_b@terms[[p]])) {
      tm <- st$terms[[p]][[nm]]
      if (is.null(tm)) next
      h <- tryCatch(modelterms7::term_held(tm), error = function(e) integer(0))
      if (!length(h)) next
      cols <- design[[p]]$blocks[[nm]]
      rd <- tryCatch(modelterms7::term_readable(tm, cf[[p]][cols]),
                     error = function(e) NULL)
      if (is.null(rd)) next
      pos <- rd$value[startsWith(rd$name, "psi")]
      if (length(pos) == length(h)) cf[[p]][cols[h]] <- as.numeric(pos)
    }
  }
  cf
}


#' The Status of a Replicate of a Special Model
#'
#' @description
#' Returns 0 for a replicate that is kept and 2 for a replicate that did not
#' reach a mode.
#'
#' @details
#' A model with a score-driven term is checked on the Newton decrement over
#' the coefficients and the term's parameters, [filter_decrement()], against
#' [mode_error_limit()]. A model with a regime term, or with a term whose block
#' moves with its coefficients, has no information that reads the mode (the
#' block is a working linearization), and its replicate is kept when every
#' coefficient and the objective are finite.
#'
#' @param mode `"structural"` or `"refresh"`.
#' @param spec_b The specification of the replicate.
#' @param r The refit.
#'
#' @return A single integer.
#'
#' @keywords internal
boot_special_status <- function(mode, spec_b, r) {
  cf <- unlist(r$coefficients)
  if (!all(is.finite(cf))) return(2L)
  if (identical(mode, "structural")) {
    if (!all(is.finite(boot_tail_values(spec_b, r$design)))) return(2L)
    su <- attr(r$design, "structural")
    if (identical(su[[1L]]$kind, "filter")) {
      d <- filter_decrement(spec_b, r$design, r$res, r$hyper)
      if (!is.finite(d) || d > mode_error_limit()) return(2L)
    }
  }
  0L
}


#' One Replicate of a Special Model
#'
#' @description
#' The body of the refit of a replicate for a model with a structural term or a
#' break-point term, as [boot_body()] runs it.
#'
#' @param object,spec_b,pl The fit, the specification and the draw of the
#'   replicate.
#' @param mode `"structural"` or `"refresh"`.
#' @param refit,method,cfg,beta0 As in [bootstrap_refit()].
#' @param keep_fit Whether to return the refit.
#'
#' @return A list as [boot_body()] returns.
#'
#' @keywords internal
boot_special_replicate <- function(object, spec_b, pl, mode, refit, method,
                                   cfg, beta0, keep_fit) {
  params <- spec_b@distrib@params
  if (!is.null(pl$seed)) set.seed(pl$seed)
  r <- tryCatch(withCallingHandlers(
    if (identical(mode, "refresh")) {
      boot_refit_breakpoint(object, spec_b, refit, method, cfg, beta0)
    } else {
      bootstrap_refit(object, spec_b, refit, method, cfg, beta0)
    },
    warning = function(w) {
      if (grepl("CHOLMOD", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }), error = function(e) NULL)
  if (is.null(r)) return(list(status = 1L))
  st <- tryCatch(boot_special_status(mode, spec_b, r),
                 error = function(e) 2L)
  if (st != 0L) return(list(status = st))
  beta <- unlist(r$coefficients[params], use.names = FALSE)
  if (identical(mode, "structural")) {
    beta <- c(beta, boot_tail_values(spec_b, r$design))
  }
  list(status = 0L, beta = beta,
       hyper = if (identical(refit, "full")) unlist(r$hyper) else NULL,
       fit = if (keep_fit) list(coefficients = r$coefficients,
                                hyper = r$hyper))
}


#' The Covariance of the Replicates of a Model With a Structural Term
#'
#' @description
#' Returns the covariance that [boot_vcov_matrix()] reports for a model whose
#' replicates carry the free parameters of a structural term after the
#' coefficients.
#'
#' @details
#' The coordinates are the joint ones of [vcov.StatmodFit()]. With
#' `readable = TRUE` the reported quantities (a persistence, a transition
#' probability) are evaluated at every replicate, from the replicate's
#' coefficients and parameters, and their covariance is taken. With
#' `readable = FALSE` the matrix is ordered by equation as [vcov.StatmodFit()]
#' orders it.
#'
#' @param object The [StatmodFit()].
#' @param B The matrix of replicates, over the joint coordinates.
#' @param readable Whether to return the reported quantities.
#' @param spec,design The specification and design of the fit.
#' @param nm The names of the coefficient coordinates.
#'
#' @return A square matrix.
#'
#' @keywords internal
boot_vcov_structural <- function(object, B, readable, spec, design, nm) {
  params <- spec@distrib@params
  lab <- coef_labels(spec, design)
  total <- nrow(lab)
  nz <- ncol(B) - total
  tl <- structural_tail_names(spec, design)
  ali <- c(nm %in% object@aliased, rep(FALSE, nz))
  B[, ali] <- NA_real_
  if (!isTRUE(readable)) {
    out <- stats::cov(B)
    rn <- c(nm, unname(tl))
    dimnames(out) <- list(rn, rn)
    who <- c(lab$parameter, names(tl))
    ord <- order(match(who, params))
    return(out[ord, ord, drop = FALSE])
  }
  rj <- readable_joint(spec, design, object)
  J <- rj$jacobian
  if (nrow(J) == 0L || ncol(J) != ncol(B)) {
    out <- stats::cov(B)
    dimnames(out) <- list(c(nm, unname(tl)), c(nm, unname(tl)))
    return(out)
  }
  sst <- statmod_structural_state(design)
  su <- attr(design, "structural")
  lens <- vapply(object@coefficients[params], length, integer(1))
  grp <- factor(rep(params, lens), levels = params)
  Q <- t(vapply(seq_len(nrow(B)), function(b) {
    cf <- split(B[b, seq_len(total)], grp)
    cf <- Map(function(v, p) stats::setNames(v, names(object@coefficients[[p]])),
              cf, params)
    at <- total
    for (u in su) {
      tm <- spec@terms[[u$param]][[u$term]]
      free <- setdiff(modelterms7::term_params(tm), sst$held[[u$term]])
      sst$zeta[[u$term]][free] <- B[b, at + seq_along(free)]
      at <- at + length(free)
    }
    readable_joint(spec, design,
                   S7::set_props(object, coefficients = cf))$value
  }, numeric(nrow(J))))
  out <- stats::cov(Q)
  ok <- !apply(B, 2L, anyNA)
  bad <- rj$held | apply(J[, !ok, drop = FALSE] != 0, 1L, any)
  out[bad, ] <- NA_real_
  out[, bad] <- NA_real_
  nmq <- paste(rj$parameter, rj$name, sep = ":")
  dimnames(out) <- list(nmq, nmq)
  out
}


#' The Replicates of a Model With a Structural Term
#'
#' @description
#' Returns the matrix of bootstrap replicates over the joint coordinates of
#' [vcov.StatmodFit()], as [boot_vcov_matrix()] reads it, with the aliased
#' coordinates missing.
#'
#' @param object The [StatmodFit()].
#' @param vt A `BootstrapVcov` object.
#'
#' @return A matrix with one row for each replicate used.
#'
#' @keywords internal
boot_replicate_matrix <- function(object, vt) {
  vt <- boot_resolve(object, vt)
  e <- boot_entries(object)
  e <- boot_match(e, vt@resample, vt@refit)[[1L]]
  B <- e@coefficients[seq_len(vt@n), , drop = FALSE]
  spec <- object@spec
  design <- statmod_design(spec)
  nm <- rownames(coef_labels(spec, design))
  nz <- ncol(B) - length(nm)
  B[, c(nm %in% object@aliased, rep(FALSE, nz))] <- NA_real_
  B
}


#' The Quantities of a Structural Term at Each Replicate
#'
#' @description
#' Evaluates, for every bootstrap replicate, the quantities that a structural
#' term reports ([modelterms7::term_readable()]) from the replicate's free
#' parameters of the term.
#'
#' @param object The [StatmodFit()].
#' @param B The matrix of replicates of [boot_replicate_matrix()].
#' @param spec,design The specification and design of the fit.
#'
#' @return A named list, one entry for each structural term, holding a matrix
#'   with one row for each replicate and one column for each quantity of the
#'   term, in the order of [modelterms7::term_readable()]. A replicate with a
#'   missing parameter of the term gives a row of missing values.
#'
#' @keywords internal
boot_term_values <- function(object, B, spec, design) {
  su <- attr(design, "structural")
  sst <- statmod_structural_state(design)
  total <- nrow(coef_labels(spec, design))
  at <- total
  out <- list()
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    nmp <- modelterms7::term_params(tm)
    free <- setdiff(nmp, sst$held[[u$term]])
    cols <- at + seq_along(free)
    at <- at + length(free)
    nq <- length(modelterms7::term_readable(tm, sst$zeta[[u$term]])$name)
    Q <- matrix(NA_real_, nrow(B), nq)
    for (b in seq_len(nrow(B))) {
      if (anyNA(B[b, cols])) next
      z <- sst$zeta[[u$term]]
      z[free] <- B[b, cols]
      Q[b, ] <- modelterms7::term_readable(tm, z)$value
    }
    out[[u$term]] <- Q
  }
  out
}
