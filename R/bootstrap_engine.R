#' @include bootstrap_objects.R
NULL

#' The Engine of the Bootstrap
#'
#' @description
#' Draws resampled data sets from a fit, refits the model on each and returns
#' the coefficients of the replicates that were kept. It is the engine that
#' [add_bootstrap()] and the predictive bootstrap of [pred_bootstrap()] share.
#'
#' @details
#' The draws are made in the calling process, one replicate after another from
#' the random number stream, and the refits go to [worker_map()] with the
#' number of workers that the fit carries. A replicate is therefore a function
#' of its position in the stream and not of the number of workers, and the
#' first `n` replicates of a stream are those of a run of `n`. The state of
#' `.Random.seed` after the last replicate that was kept is returned, so that a
#' later call continues the stream: 500 replicates followed by 300 more are the
#' 800 replicates of one run from the same seed. A refit does not draw from the
#' stream of the calling process: the seed is saved before it and restored
#' after it.
#'
#' A replicate fails when the refit signals an error (status 1), when the mode
#' error of the refit exceeds [mode_error_limit()] (status 2), or when a
#' coefficient is unidentified in the replicate and not in the fit (status 3),
#' which is the case of a level of a factor whose observations have all a zero
#' weight. The information of the mode error is that of the inner fit. A
#' failed replicate is replaced by the next draw of the stream with
#' `replace = TRUE`, and the engine stops with an error when the failures
#' outnumber the successes after at least ten draws.
#'
#' @param object The [StatmodFit()].
#' @param n The number of replicates to keep.
#' @param resample A [resample_method] object.
#' @param refit `"coefficients"` or `"full"`.
#' @param seed `NULL`, or the value of `.Random.seed` from which the stream
#'   continues.
#' @param aside The rows of [random_modes()] whose effects are redrawn by a
#'   parametric replicate; `NULL` for all of them.
#' @param replace `TRUE` replaces a failed replicate by the next draw; `FALSE`
#'   makes exactly `n` draws and leaves the failures out.
#' @param validate `TRUE` applies the mode error and the identification
#'   criteria; `FALSE` fails a replicate on an error only.
#' @param hook `NULL`, or a function of a refit and its position that is
#'   called in the calling process for every replicate kept.
#' @param batch The largest number of draws made before their refits run;
#'   `NULL` draws all that are missing.
#' @param stream `TRUE` manages the state of the random number generator of the
#'   session as described above.
#' @param adjust `TRUE` stores a random-effect coordinate of a parametric
#'   replicate as \eqn{\hat b^* - b^* + \hat b}.
#' @param check `TRUE` rejects the fits that the bootstrap does not handle.
#'
#' @return A list with `coefficients` (a matrix with a row for each replicate
#'   kept), `hyper` (a matrix, or `NULL`), `status` (an integer for every
#'   draw), `seed_start`, `seed_end`, `random`, `raw` (the matrix of the
#'   random-effect coefficients of the replicates as refitted, without the
#'   adjustment of the random effects, for a parametric scheme with
#'   `adjust = TRUE` beside random effects, and `NULL` otherwise) and
#'   `raw_cols` (the positions of the columns of `raw` among the columns of
#'   `coefficients`).
#'
#' @seealso [add_bootstrap()], [predictive_bootstrap()], [worker_map()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' eng <- statmodels7:::bootstrap_engine(fit, 5, case_resample(),
#'                                       "coefficients")
#' dim(eng$coefficients)
#'
#' @keywords internal
bootstrap_engine <- function(object, n, resample, refit, seed = NULL,
                             aside = NULL, replace = TRUE, validate = TRUE,
                             hook = NULL, batch = NULL, stream = TRUE,
                             adjust = TRUE, check = TRUE) {
  spec0 <- object@spec
  design0 <- statmod_design(spec0)
  if (check) bootstrap_check(object, resample, refit, spec0, design0)
  object <- boot_light(object)
  if (is.null(aside)) aside <- random_modes(spec0, "zero")
  sm <- object@methods
  params <- spec0@distrib@params
  method <- sm$smooth
  cfg <- inner_settings(method)
  beta0 <- unlist(object@coefficients[params], use.names = FALSE)
  draw <- boot_plan_wrap(
    bootstrap_plan(object, resample, spec0, design0, aside, adjust), design0)
  ref <- list(def = integer(0), singular = FALSE)
  if (validate && is.null(boot_mode(design0))) {
    ref <- boot_read_mode(object, spec0, design0, object@coefficients,
                          object@hyper, cfg)
  }
  if (stream) {
    if (!is.null(seed)) {
      assign(".Random.seed", seed, envir = globalenv())
    } else if (!exists(".Random.seed", envir = globalenv(),
                       inherits = FALSE)) {
      stats::runif(1)
    }
    seed_start <- get(".Random.seed", envir = globalenv())
  } else seed_start <- NULL
  seed_end <- seed_start
  rows <- list()
  raws <- list()
  hypers <- list()
  kept <- list()
  status <- integer(0)
  accepted <- 0L
  repeat {
    need <- if (replace) n - accepted else n - length(status)
    if (need <= 0L) break
    k <- if (is.null(batch)) need else min(as.integer(batch), need)
    plans <- vector("list", k)
    seeds <- vector("list", k)
    for (j in seq_len(k)) {
      plans[[j]] <- draw()
      if (stream) seeds[[j]] <- get(".Random.seed", envir = globalenv())
    }
    body <- boot_body(object, spec0, draw, refit, method, cfg, beta0, plans,
                      ref, validate, isolate = stream,
                      keep_fit = !is.null(hook))
    res <- worker_map(spec0, k, body, what = "bootstrap replicas")
    for (j in seq_len(k)) {
      st <- res[[j]]$status
      status <- c(status, st)
      if (stream) seed_end <- seeds[[j]]
      if (st != 0L) next
      accepted <- accepted + 1L
      rows[[accepted]] <- res[[j]]$beta
      raws[[accepted]] <- res[[j]]$raw
      hypers[[accepted]] <- res[[j]]$hyper
      if (!is.null(hook)) hook(res[[j]]$fit, accepted)
      if (replace && accepted >= n) break
    }
    if (replace && accepted < n && length(status) >= 10L &&
        sum(status > 0L) > sum(status == 0L)) {
      stop(sprintf(paste0("More bootstrap replicas failed than succeeded ",
                          "(%d failed, %d kept): %d with an error, %d with a ",
                          "mode error above the limit and %d with an ",
                          "unidentified coefficient."),
                   sum(status > 0L), accepted, sum(status == 1L),
                   sum(status == 2L), sum(status == 3L)), call. = FALSE)
    }
  }
  if (stream) assign(".Random.seed", seed_end, envir = globalenv())
  nm <- boot_coord_names(spec0, design0)
  coefficients <- if (length(rows)) {
    m <- do.call(rbind, rows)
    dimnames(m) <- list(NULL, nm)
    m
  } else matrix(numeric(0), 0L, length(nm), dimnames = list(NULL, nm))
  hyper <- if (identical(refit, "full") && length(hypers)) {
    do.call(rbind, hypers)
  } else NULL
  raw_cols <- boot_raw_cols(attr(draw, "targets"), object@coefficients,
                            spec0@distrib@params)
  raw <- if (length(rows) && length(raws) == length(rows) &&
             !any(vapply(raws, is.null, TRUE))) {
    m <- do.call(rbind, raws)
    dimnames(m) <- list(NULL, nm[raw_cols])
    m
  } else NULL
  list(coefficients = coefficients, hyper = hyper, status = status,
       seed_start = seed_start, seed_end = seed_end,
       random = attr(draw, "random"), raw = raw,
       raw_cols = if (is.null(raw)) NULL else raw_cols)
}


#' The Fit as the Workers Receive It
#'
#' @description
#' Returns the fit without the optimizer of its outer search.
#'
#' @details
#' The closures of the engine hold the fit, and [worker_map()] serializes them
#' for every worker. The optimizer that the outer search stores in
#' `object@methods$search` keeps the environment of that search, which holds the
#' design and the inner fit, and it is not read by a refit. Dropping it leaves
#' the properties that [bootstrap_refit()], [boot_read_mode()] and the
#' draws read.
#'
#' @param object The [StatmodFit()].
#'
#' @return The [StatmodFit()] with `methods$search` set to `NULL`.
#'
#' @keywords internal
boot_light <- function(object) {
  m <- object@methods
  m["search"] <- list(NULL)
  S7::set_props(object, methods = m)
}


#' Reject a Fit That the Bootstrap Does Not Handle
#'
#' @description
#' Signals an error for a model that the bootstrap of [add_bootstrap()] does
#' not cover.
#'
#' @details
#' A response that is not one number per observation is rejected. So are
#' `refit = "full"` with a [cv()] criterion, which needs the folds of the
#' data, the case scheme beside a [modelterms7::random()] term, and the
#' cluster scheme beside a random-effect term that
#' [boot_cluster_random_check()] cannot regroup. A model with a structural
#' term ([modelterms7::gas()], [modelterms7::regime()]) takes
#' [parametric_resample()] only; a model with a block that moves with its
#' coefficients ([modelterms7::seg()], [modelterms7::jump()],
#' [modelterms7::jseg()], [modelterms7::nl()]) takes any scheme and is
#' refitted by [boot_refit_breakpoint()].
#'
#' @param object The [StatmodFit()].
#' @param resample The [resample_method].
#' @param refit `"coefficients"` or `"full"`.
#' @param spec0,design0 Its specification and design.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
bootstrap_check <- function(object, resample, refit, spec0, design0) {
  if (length(attr(design0, "structural")) &&
      !S7::S7_inherits(resample, ParametricResample)) {
    stop(paste0("With a structural term (gas(), regime()) a bootstrap replicate ",
                "needs a new path of the series, because a resample of rows ",
                "or of groups breaks the recursion of the filter. Use ",
                "parametric_resample()."), call. = FALSE)
  }
  y0 <- spec0@response
  if (!is.numeric(y0) || !is.null(dim(y0)) || length(y0) != spec0@n_obs) {
    stop(paste0("The bootstrap resamples the observations, and this ",
                "response is not one number per observation."), call. = FALSE)
  }
  sm <- object@methods
  is_cv <- function(m) !is.null(m) && identical(m@kind, "cv")
  if (identical(refit, "full") &&
      (is_cv(sm$outer) || is_cv(sm$sparse_criterion))) {
    stop(paste0("With refit = \"full\" the hyperparameters are chosen again ",
                "on every replicate, and\n  a cross-validated criterion ",
                "needs the data's folds, which a fit does not keep.\n  ",
                "refit = \"coefficients\" holds them at the fit's values."),
         call. = FALSE)
  }
  has_random <- nrow(random_modes(spec0, "zero")) > 0L ||
    nrow(attr(random_modes(spec0, "zero"), "nested")) > 0L
  if (has_random && (S7::S7_inherits(resample, CaseResample) ||
                     S7::S7_inherits(resample, BayesianResample))) {
    stop(paste0(resample@name, "_resample() breaks the groups of a random() ",
                "term. cluster_resample() resamples the groups."),
         call. = FALSE)
  }
  if (has_random && S7::S7_inherits(resample, ClusterResample)) {
    boot_cluster_random_check(spec0)
  }
  invisible(NULL)
}


#' The Draw of One Bootstrap Replicate
#'
#' @description
#' Returns a function of no arguments that draws the description of one
#' replicate: the simulated response of a parametric replicate, or the
#' frequency weights of a case or a cluster replicate.
#'
#' @details
#' The weights multiply the weights of the fit. A cluster replicate gives the
#' observations of a group the number of times that the group was drawn. The
#' attribute `kind` of the result is `"parametric"` or `"weights"`; `targets`
#' lists the random-effect coordinates of a parametric draw, and `random` says
#' whether the effects are redrawn. A cluster scheme beside random effects
#' returns the plan of [boot_cluster_random_plan()], of kind `"regrouped"`.
#'
#' @param object The [StatmodFit()].
#' @param resample The [resample_method].
#' @param spec0,design0 Its specification and design.
#' @param aside The rows of [random_modes()] that a parametric draw redraws.
#' @param adjust Whether a parametric plan carries the effects drawn.
#'
#' @return A function.
#'
#' @keywords internal
bootstrap_plan <- function(object, resample, spec0, design0, aside, adjust) {
  if (identical(boot_mode(design0), "structural")) {
    return(boot_plan_structural(object, resample, spec0, design0))
  }
  if (S7::S7_inherits(resample, ParametricResample)) {
    sim <- bootstrap_simulator(object, spec0, design0, aside, keep = adjust)
    out <- if (adjust) {
      function() sim()
    } else {
      function() list(y = sim())
    }
    attr(out, "kind") <- "parametric"
    attr(out, "targets") <- if (adjust) attr(sim, "targets") else list()
    attr(out, "random") <- if (nrow(aside)) "redrawn" else "none"
    return(out)
  }
  n <- spec0@n_obs
  w0 <- spec0@weights
  if (S7::S7_inherits(resample, ClusterResample) &&
      nrow(random_modes(spec0, "zero")) > 0L) {
    return(boot_cluster_random_plan(object, resample, spec0, design0))
  }
  if (S7::S7_inherits(resample, CaseResample)) {
    out <- function() list(w = tabulate(sample.int(n, n, TRUE), n) * w0)
  } else if (S7::S7_inherits(resample, BayesianResample)) {
    pos <- w0 > 0
    np <- sum(pos)
    out <- function() {
      e <- stats::rexp(np)
      w <- numeric(n)
      w[pos] <- w0[pos] * e * (np / sum(e))
      list(w = w)
    }
  } else {
    id <- boot_cluster_ids(object, resample, spec0)
    gi <- match(id, unique(id))
    ng <- max(gi)
    out <- function() {
      list(w = tabulate(sample.int(ng, ng, TRUE), ng)[gi] * w0)
    }
  }
  attr(out, "kind") <- "weights"
  attr(out, "targets") <- list()
  attr(out, "random") <- "none"
  out
}


#' The Groups of a Cluster Bootstrap
#'
#' @description
#' Evaluates the variable of [cluster_resample()] on the data of the fit.
#'
#' @details
#' The data frame is the `data` of the scheme if given, and otherwise the
#' `data` of the call that built the fit, evaluated in the environment of the
#' model formula. The fit does not keep its data, so a `data` argument that no
#' longer evaluates is an error that says to pass `data`.
#'
#' @param object The [StatmodFit()].
#' @param resample A [ClusterResample].
#' @param spec0 Its specification.
#'
#' @return A vector with one value per observation.
#'
#' @keywords internal
boot_cluster_ids <- function(object, resample, spec0) {
  data <- boot_cluster_data(object, resample, spec0)
  v <- all.vars(resample@cluster)
  if (!v %in% names(data)) {
    stop(sprintf("The cluster variable '%s' is not a column of the data.", v),
         call. = FALSE)
  }
  id <- data[[v]]
  if (anyNA(id)) {
    stop(sprintf("The cluster variable '%s' has missing values.", v),
         call. = FALSE)
  }
  id
}

#' @rdname boot_cluster_ids
#' @keywords internal
boot_cluster_data <- function(object, resample, spec0) {
  n <- spec0@n_obs
  data <- resample@data
  if (is.null(data)) {
    data <- tryCatch({
      env <- environment(spec0@formula)
      if (is.null(env)) env <- globalenv()
      eval(object@call$data, env)
    }, error = function(e) NULL)
  }
  if (!is.data.frame(data) || nrow(data) != n) {
    stop(paste0("The data frame with the cluster variable could not be ",
                "recovered from the call of the fit. Pass it as ",
                "cluster_resample(cluster = , data = ), with one row for ",
                "each of the ", n, " observations."), call. = FALSE)
  }
  data
}


#' The Function That Refits the Replicates of a Round
#'
#' @description
#' Returns the function of a replicate's position that [worker_map()] runs. It
#' is built here so that its environment holds only what the refit needs, and
#' that is what a worker process receives.
#'
#' @param object,spec0 The fit and its specification.
#' @param draw The function from [bootstrap_plan()], whose attributes give the
#'   kind of plan and the targets.
#' @param refit,method,cfg,beta0 As in [bootstrap_refit()].
#' @param plans The draws of the round, or a function of the position that
#'   returns the draw.
#' @param ref The reading of the fit itself, from [boot_read_mode()].
#' @param validate,isolate,keep_fit As in [bootstrap_engine()], where `isolate`
#'   is `stream`.
#'
#' @return A function of the position of a replicate in `plans`, returning a
#'   list with `status`, `beta`, `hyper` and `fit`.
#'
#' @keywords internal
boot_body <- function(object, spec0, draw, refit, method, cfg, beta0, plans,
                      ref, validate, isolate, keep_fit) {
  kind <- attr(draw, "kind")
  mode <- attr(draw, "mode")
  targets <- attr(draw, "targets")
  build <- attr(draw, "build")
  params <- spec0@distrib@params
  cf_fit <- object@coefficients
  raw_cols <- boot_raw_cols(targets, cf_fit, params)
  force(plans)
  function(i) {
    pl <- if (is.function(plans)) plans(i) else plans[[i]]
    if (isolate && exists(".Random.seed", envir = globalenv(),
                          inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    }
    cb <- NULL
    start <- beta0
    spec_b <- if (identical(kind, "parametric")) {
      S7::set_props(spec0, response = pl$y)
    } else if (identical(kind, "regrouped")) {
      cb <- tryCatch(build(pl), error = function(e) NULL)
      if (is.null(cb)) return(list(status = 1L))
      start <- cb$start
      cb$spec
    } else S7::set_props(spec0, weights = pl$w)
    if (!is.null(mode)) {
      return(boot_special_replicate(object, spec_b, pl, mode, refit, method,
                                    cfg, beta0, keep_fit))
    }
    # a trial point of a refit may meet a penalized information the sparse
    # factorization reports as not positive definite; the refit steps away
    # from it, and the warning says nothing about the replica
    r <- tryCatch(withCallingHandlers(
      bootstrap_refit(object, spec_b, refit, method, cfg, start),
      warning = function(w) {
        if (grepl("CHOLMOD", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }), error = function(e) NULL)
    if (is.null(r)) return(list(status = 1L))
    status <- 0L
    if (validate) {
      status <- tryCatch(
        if (is.null(cb)) boot_validate(object, spec_b, r, cfg, ref)
        else boot_cluster_validate(object, spec_b, r, cfg, ref, cb),
        error = function(e) 2L)
    }
    if (status != 0L) return(list(status = status))
    cf <- r$coefficients
    # the replicate of a random effect is b_hat* - b* + b_hat, whose
    # covariance is that of the prediction error
    if (length(targets) && !is.null(pl$cf)) {
      for (tg in targets) {
        cf[[tg$param]][tg$pos] <- cf[[tg$param]][tg$pos] -
          pl$cf[[tg$param]][tg$pos] + cf_fit[[tg$param]][tg$pos]
      }
    }
    if (!is.null(cb)) {
      return(list(status = 0L, beta = cb$expand(cf),
                  hyper = if (identical(refit, "full")) unlist(r$hyper) else NULL,
                  fit = if (keep_fit) list(coefficients = r$coefficients,
                                           hyper = r$hyper)))
    }
    list(status = 0L, beta = unlist(cf[params], use.names = FALSE),
         raw = if (length(targets) && !is.null(pl$cf)) {
           unlist(r$coefficients[params], use.names = FALSE)[raw_cols]
         } else NULL,
         hyper = if (identical(refit, "full")) unlist(r$hyper) else NULL,
         fit = if (keep_fit) list(coefficients = r$coefficients,
                                  hyper = r$hyper))
  }
}


#' The Columns of the Refitted Random Effects
#'
#' @description
#' Returns the positions, among the stacked coefficients, of the
#' random-effect coordinates that a parametric replicate redraws.
#'
#' @param targets The coordinates of the draw, from `attr(draw, "targets")`.
#' @param cf The coefficients of the fit, a list by equation.
#' @param params The names of the equations, in the order of the stack.
#'
#' @return A sorted integer vector, empty when there are no targets.
#'
#' @keywords internal
boot_raw_cols <- function(targets, cf, params) {
  if (!length(targets)) return(integer(0))
  lens <- vapply(cf[params], length, integer(1))
  offs <- cumsum(lens) - lens
  sort(unique(unlist(lapply(targets, function(tg)
    offs[[match(tg$param, params)]] + tg$pos))))
}


#' The Status of One Refitted Replicate
#'
#' @description
#' Returns 0 for a replicate that is kept, 2 when its mode error exceeds
#' [mode_error_limit()] or cannot be read, and 3 when it leaves a coordinate
#' unidentified that the fit identifies, or has an information that cannot be
#' factorized where the fit's can.
#'
#' @param object The [StatmodFit()].
#' @param spec_b The specification of the replicate.
#' @param r The refit, as [bootstrap_refit()] returns it.
#' @param cfg The inner settings of the fit.
#' @param ref The reading of the fit itself, from [boot_read_mode()].
#'
#' @return A single integer.
#'
#' @keywords internal
boot_validate <- function(object, spec_b, r, cfg, ref) {
  rd <- boot_read_mode(object, spec_b, r$design, r$coefficients, r$hyper, cfg)
  # an information that cannot be factorized where the fit's can is a
  # coefficient that the replicate does not identify
  if (isTRUE(rd$singular)) return(if (isTRUE(ref$singular)) 0L else 3L)
  if (length(setdiff(rd$def, ref$def))) return(3L)
  if (!is.finite(rd$me) || rd$me > mode_error_limit()) return(2L)
  0L
}


#' How Far From Its Mode a Fit Is, and What It Does Not Identify
#'
#' @description
#' Returns \eqn{\tfrac12 g'K^{-1}g} at the coefficients, as the certificate
#' reads it, and the coordinates that the penalized information does not
#' identify there.
#'
#' @details
#' The coefficients that the marginal criterion estimates are held, as the
#' certificate holds them, and the coordinates of a kinked penalty that are
#' zero are left out of the score. The information is that of the inner fit.
#' The unidentified coordinates are those that [ctx_penalized()] pins, apart
#' from the ones pinned for another reason: a coefficient whose observations
#' all have a zero weight is among them, although the inner fit does not report
#' it as aliased.
#'
#' @param object The [StatmodFit()], for its criterion and pinned coordinates.
#' @param spec,design The specification and design where the point is read.
#' @param cf,hy The coefficients and hyperparameters.
#' @param cfg The inner settings of the fit.
#' @param pinned The pinned coordinates at the point read, in its own stacked
#'   positions; those of the fit by default.
#'
#' @return A list with `me` (a number, `NA` if the point could not be read),
#'   `def` (stacked positions) and `singular` (whether the penalized
#'   information could not be factorized).
#'
#' @keywords internal
boot_read_mode <- function(object, spec, design, cf, hy, cfg,
                           pinned = object@methods$pinned) {
  method <- object@methods$outer
  gam <- tryCatch(marginal_coords(spec, design, method),
                  error = function(e) marginal_coords(spec, design, NULL))
  ng <- length(gam$where)
  beta_all <- unlist(cf[spec@distrib@params], use.names = FALSE)
  sp <- if (ng) statmod_hold(spec, gam, beta_all[gam$where]) else spec
  if (length(pinned)) attr(design, "pinned") <- pinned
  bad <- list(me = NA_real_, def = integer(0), singular = FALSE)
  ctx <- tryCatch(outer_context(sp, design, cf, hy, NULL),
                  error = function(e) NULL)
  if (is.null(ctx)) return(bad)
  pen <- tryCatch(ctx_penalized(ctx, sp, design, cf, hy, cfg$expected),
                  error = function(e) NULL)
  if (is.null(pen)) return(list(me = NA_real_, def = integer(0),
                                singular = TRUE))
  obj <- statmod_objective(sp, hy, design, cfg$expected, cfg$emethod)
  sc <- obj$gr(obj$stack(cf))
  if (ng) sc[gam$where] <- 0
  sc <- free_of_kinks(sc, sp, design, cf)
  if (!all(is.finite(sc))) return(bad)
  db <- as.numeric(as.matrix(pen$inv) %*% sc)
  known <- c(pinned_coords(sp, design), zero_kinked(sp, design, cf),
             attr(pen$K, "boundary_at"))
  list(me = 0.5 * sum(sc * db),
       def = setdiff(attr(pen$K, "held_at"), known), singular = FALSE)
}
