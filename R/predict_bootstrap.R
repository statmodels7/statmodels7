#' @include predict_interval.R
NULL

#' The Parametric Bootstrap of a Prediction Interval
#'
#' @description
#' Returns the components of the predictive mixture that
#' [predictive_response()] reads for [pred_bootstrap()]. There is one
#' component per replica, which is the fit refitted on a response simulated
#' from the fitted model and read at the rows predicted with no estimation
#' uncertainty of its own.
#'
#' @details
#' A replica draws the response from the family at the fitted model. The
#' effects of a [modelterms7::random()] term that `random` sets aside are
#' drawn afresh from the prior that the fit estimated, so the groups of the
#' data are new groups in every replica, as they are in the sampling
#' distribution of the population-level estimates on which a new group's
#' prediction rests. The
#' effects of a term read `"conditional"` keep their estimates: the
#' prediction is about those groups, and their effects are the quantity
#' predicted rather than a source of variation. Every other coefficient keeps
#' its estimate. The replica is then refitted in one of two ways:
#' \describe{
#'   \item{`"coefficients"`}{every coefficient is re-estimated at the
#'     hyperparameters that the fit reached. A coefficient that the fit
#'     estimated on its criterion (a dispersion under [reml()]) is estimated
#'     there again while the hyperparameters are held, so the replica's
#'     estimator is the fit's.}
#'   \item{`"full"`}{the hyperparameters are chosen again by the fit's own
#'     criteria, so their uncertainty enters the interval too, and the path
#'     of a kinked penalty is walked again, warm started from the fit's
#'     coefficients. A cross-validated criterion needs the data's folds and
#'     is rejected.}
#' }
#' Every refit starts at the optimum that the fit found: its coefficients,
#' including those that the fit estimated on its criterion, and its
#' hyperparameters.
#' The bootstrap is rejected for a model with a block that moves with its
#' coefficients ([modelterms7::seg()], [modelterms7::jump()],
#' [modelterms7::jseg()], [modelterms7::nl()]) and for a response that is not
#' one number per observation.
#'
#' Each replica gives the component \eqn{(\hat\eta^*_b, C^*_b)}: its
#' predictors at the rows predicted and the covariance that the effects set
#' aside add, at its own hyperparameters. The interval is the mixture's, which
#' is the parametric bootstrap predictive distribution of Harris (1989); the
#' default [pred_averaged()] is its normal approximation, in which the
#' sampling distribution of the estimates is replaced by \eqn{\mathrm{N}(\hat
#' \beta, \widehat{V})}. A replica whose refit fails is left out, and the
#' warning reports how many were.
#'
#' @param object The [StatmodFit()].
#' @param spec,design The specification and design at the rows predicted,
#'   in which the columns of the effects set aside are zero.
#' @param aside The rows of [random_modes()] that are not `"conditional"`.
#' @param n_boot The number of replicas, from [pred_bootstrap()].
#' @param refit `"coefficients"` or `"full"`.
#'
#' @return A list of components, each a list with `eta`, `C` and `weight`.
#'   The list has the attribute `infinite_variance`, which is `TRUE` where a
#'   prior has no variance. The weights of the components add up to one.
#'
#' @references Harris, I. R. (1989). Predictive fit for natural exponential
#'   families. *Biometrika*, 76, 675--684.
#'
#' @seealso [predict.StatmodFit()], the caller, [predictive_response()] for
#'   the mixture.
#'
#' @examples
#' set.seed(2)
#' dd <- data.frame(x = runif(30))
#' dd$y <- 1 + dd$x + rnorm(30, sd = 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' predict(fit, what = "response", newdata = data.frame(x = 0.5),
#'         interval = prediction_interval(
#'           predictive = pred_bootstrap(n = 20)))
#'
#' @keywords internal
predictive_bootstrap <- function(object, spec, design, aside, n_boot, refit) {
  spec0 <- object@spec
  design0 <- statmod_design(spec0)
  if (length(attr(design0, "refresh"))) {
    stop(paste0("pred_bootstrap() is not available for a term ",
                "whose block moves with\n  its coefficients (seg(), jump(), ",
                "jseg(), nl()): a replica would refit\n  its break-points ",
                "or its nonlinear parameters from where the fit left them,\n",
                "  which is a different estimator. pred_averaged() ",
                "is available."), call. = FALSE)
  }
  y0 <- spec0@response
  if (!is.numeric(y0) || !is.null(dim(y0)) ||
      length(y0) != spec0@n_obs) {
    stop(paste0("pred_bootstrap() simulates the response, and ",
                "this response is not\n  one number per observation. ",
                "pred_averaged() is available."), call. = FALSE)
  }
  sm <- object@methods
  is_cv <- function(m) !is.null(m) && identical(m@kind, "cv")
  if (identical(refit, "full") &&
      (is_cv(sm$outer) || is_cv(sm$sparse_criterion))) {
    stop(paste0("With refit = \"full\" the hyperparameters are chosen ",
                "again on every replica, and\n  a cross-validated criterion ",
                "needs the data's folds, which a fit does not\n  keep. ",
                "refit = \"coefficients\" holds them at the fit's ",
                "values."), call. = FALSE)
  }
  draw <- bootstrap_simulator(object, spec0, design0, aside)
  params <- spec0@distrib@params
  method <- sm$smooth
  cfg <- inner_settings(method)
  # every replica starts at the fit's optimum: its coefficients here, its
  # hyperparameters through object@hyper in bootstrap_refit()
  beta0 <- unlist(object@coefficients[params], use.names = FALSE)
  comps <- list()
  inf_var <- FALSE
  failed <- 0L
  for (b in seq_len(n_boot)) {
    y <- draw()
    spec_b <- S7::set_props(spec0, response = y)
    # a trial point of a refit may meet a penalized information the sparse
    # factorization reports as not positive definite; the refit steps away
    # from it, and the warning says nothing about the replica
    r <- tryCatch(withCallingHandlers(
      bootstrap_refit(object, spec_b, refit, method, cfg, beta0),
      warning = function(w) {
        if (grepl("CHOLMOD", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }), error = function(e) NULL)
    if (is.null(r)) {
      failed <- failed + 1L
      next
    }
    fit_b <- S7::set_props(object, coefficients = r$coefficients,
                           hyper = r$hyper)
    cb <- predictive_mixture(fit_b, spec, design, aside,
                             statmod_eta(spec, design, r$coefficients)$eta,
                             NULL, weight = 1,
                             n_draw = max(1L, ceiling(5000 / n_boot)))
    inf_var <- inf_var || isTRUE(attr(cb, "infinite_variance"))
    comps <- c(comps, cb)
  }
  if (!length(comps)) {
    stop("Every bootstrap replica failed to refit.", call. = FALSE)
  }
  if (failed) {
    warning(sprintf(paste0("%d of %d bootstrap replicas could not be ",
                           "refitted and are left out."), failed, n_boot),
            call. = FALSE)
  }
  # every replica carries weight one among its own components, so the
  # replicas count alike however many components each has
  nb <- n_boot - failed
  comps <- lapply(comps, function(cc) {
    cc$weight <- cc$weight / nb
    cc
  })
  attr(comps, "infinite_variance") <- inf_var
  comps
}


#' A Simulator of the Response From a Fitted Model
#'
#' @description
#' Returns a function of no arguments that draws one response at the fitting
#' rows. The function draws the effects of the [modelterms7::random()] terms
#' in `aside` from the prior that the fit estimated, keeps every other
#' coefficient at its estimate, and draws the response from the family at the
#' parameters that these coefficients give.
#'
#' @details
#' The effects are drawn by reading the penalty as a prior at the fit's
#' hyperparameters with [penalties7::penalty_draw()], and are placed by the
#' map that [rstatmod()] uses, so a structured covariance (an AR(1), a
#' compound symmetry), a covariance that a label shares between terms and a
#' prior that is not Gaussian (a multivariate Student t, a Laplace) are drawn
#' alike. A smooth or a ridge is not a random effect and keeps its estimate.
#'
#' @param object The [StatmodFit()].
#' @param spec0,design0 Its specification and design at the fitting rows.
#' @param aside The rows of [random_modes()] whose effects are redrawn; the
#'   other effects keep their estimates.
#'
#' @return A function returning a numeric vector of `spec0@n_obs` values.
#'
#' @examples
#' set.seed(3)
#' gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
#' gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40)
#' fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
#'                gg)
#' rm <- statmodels7:::random_modes(fit@spec, "zero")
#' draw <- statmodels7:::bootstrap_simulator(fit, fit@spec,
#'                                           statmod_design(fit@spec), rm)
#' length(draw())
#'
#' @keywords internal
bootstrap_simulator <- function(object, spec0, design0, aside) {
  params <- spec0@distrib@params
  set_aside <- paste(aside$param, aside$key)
  is_random <- function(u) {
    if (!is.null(u$class)) {
      return(any(vapply(u$class$pieces, function(pc)
        paste(pc$param, pc$term) %in% set_aside, logical(1))))
    }
    tm <- spec0@terms[[u$param]][[u$term]]
    !is.null(tm) && S7::S7_inherits(tm, modelterms7::RandomTerm) &&
      paste(u$param, u$term) %in% set_aside
  }
  units <- Filter(function(u) !isTRUE(u$structural) && is_random(u),
                  statmod_penalized(spec0, design0))
  draws <- lapply(units, function(u) {
    tg <- unit_draw_targets(u, design0, params)
    if (is.null(tg) || tg$space != "beta") {
      stop(sprintf(paste0("The effects of '%s' cannot be drawn for a ",
                          "bootstrap replica."), u$key), call. = FALSE)
    }
    list(pen = u$penalty, th = as.list(object@hyper[[u$param]][[u$key]]),
         tg = tg)
  })
  cf0 <- object@coefficients
  function() {
    cf <- cf0
    for (d in draws) {
      v <- penalties7::penalty_draw(d$pen, d$th)
      if (length(v) != length(d$tg$pos) || anyNA(v)) {
        stop("A random effect's prior did not give a draw for every effect.",
             call. = FALSE)
      }
      for (p in unique(d$tg$param)) {
        j <- d$tg$param == p
        cf[[p]][d$tg$pos[j]] <- v[j]
      }
    }
    th <- statmod_eta(spec0, design0, cf)$theta
    distributions7::distrib_rng(spec0@distrib, spec0@n_obs, th)
  }
}


#' Refit One Bootstrap Replica
#'
#' @description
#' Re-estimates the coefficients of a fit on a simulated response, either at
#' the fit's hyperparameters or after choosing the hyperparameters again. See
#' [predictive_bootstrap()].
#'
#' @param object The [StatmodFit()].
#' @param spec_b Its specification with the simulated response.
#' @param refit `"coefficients"` or `"full"`.
#' @param method The fit's inner optimizer.
#' @param cfg Its [inner_settings()].
#' @param beta0 The fit's coefficients, the replica's start.
#'
#' @return A list with `coefficients` and `hyper`.
#'
#' @examples
#' set.seed(4)
#' dd <- data.frame(x = runif(30))
#' dd$y <- 1 + dd$x + rnorm(30)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' sp <- S7::set_props(fit@spec, response = dd$y + rnorm(30, sd = 0.1))
#' statmodels7:::bootstrap_refit(fit, sp, "coefficients", fit@methods$smooth,
#'   statmodels7:::inner_settings(fit@methods$smooth),
#'   unlist(fit@coefficients, use.names = FALSE))$coefficients
#'
#' @keywords internal
bootstrap_refit <- function(object, spec_b, refit, method, cfg, beta0) {
  design <- statmod_design(spec_b)
  blocks <- statmod_blocks(spec_b, design)
  sm <- object@methods
  hyper <- object@hyper
  vb <- verbosity(0)
  crit <- sm$outer
  marg <- !is.null(crit) &&
    length(marginal_coords(spec_b, design,
                           if (outer_minimize(crit)) NULL else crit)$where) > 0L
  if (identical(refit, "full") &&
      (!is.null(crit) || !is.null(sm$sparse_criterion))) {
    res <- statmod_select(spec_b, design, blocks, hyper, method, crit, NULL,
                          beta0, cfg$approx, cfg$maxit, cfg$tol, vb, NULL,
                          sm$sparse_criterion)
    hyper <- res$hyper
  } else if (marg) {
    res <- outer_fit(spec_b, design, blocks, hyper, method, crit, NULL, beta0,
                     cfg$approx, cfg$maxit, cfg$tol, vb, hold_hyper = TRUE)
  } else {
    res <- statmod_alternate(spec_b, design, blocks, hyper, method, beta0,
                             cfg$expected, cfg$approx, cfg$maxit, cfg$tol, vb)
  }
  list(coefficients = res$obj$split(res$par), hyper = hyper)
}
