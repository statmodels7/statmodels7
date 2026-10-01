#' @include predict_interval.R
NULL

#' The Parametric Bootstrap of a Prediction Interval
#'
#' @description
#' The components of the predictive mixture [predictive_response()] reads
#' for `predictive = "bootstrap"`: one per replica, each the fit refitted on a
#' response simulated from the fitted model, read at the rows predicted with
#' no estimation uncertainty of its own.
#'
#' @details
#' A replica draws the effects of every [modelterms7::random()] term afresh
#' from the Gaussian prior the fit estimated, keeps every other coefficient at
#' its estimate, and draws the response from the family at the parameters
#' this gives, so the groups of the data are new groups in every replica, as
#' they are in the sampling distribution the bootstrap stands in for. The
#' replica is then refitted in one of two ways:
#' \describe{
#'   \item{`"coefficients"`}{at the hyperparameters the fit reached, every
#'     coefficient re-estimated. A coefficient the fit estimated on its
#'     criterion (a dispersion under [reml()]) is estimated there again, the
#'     hyperparameters held, so the replica's estimator is the fit's.}
#'   \item{`"full"`}{with the hyperparameters chosen again by the fit's own
#'     criteria, so their uncertainty enters the interval too. A
#'     cross-validated criterion needs the data's folds and is rejected.}
#' }
#' Every refit starts at the optimum the fit found: its coefficients, the
#' criterion's own coefficients among them, and its hyperparameters. A
#' kinked penalty's path is the one search that walks its own grid, warm
#' started from the fit's coefficients.
#' Each replica gives the component \eqn{(\hat\eta^*_b, C^*_b)}: its
#' predictors at the rows predicted and the covariance the effects set aside
#' add, at its own hyperparameters. The interval is the mixture's, which is
#' the parametric bootstrap predictive distribution of Harris (1989); the
#' default `predictive = "averaged"` is its normal approximation, the
#' sampling distribution of the estimates replaced by \eqn{\mathrm{N}(\hat
#' \beta, \widehat{V})}. A replica whose refit fails is left out, and a
#' warning says how many were.
#'
#' @param object The [StatmodFit()].
#' @param spec,design The specification and design at the rows predicted,
#'   with the columns of the effects set aside at zero.
#' @param aside The rows of [random_modes()] that are not `"conditional"`.
#' @param n_boot The number of replicas.
#' @param refit `"coefficients"` or `"full"`.
#'
#' @return A list of components, each a list with `eta` and `C`.
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
#' predict(fit, "response", data.frame(x = 0.5), interval = "prediction",
#'         predictive = "bootstrap", n_boot = 20)
#'
#' @keywords internal
predictive_bootstrap <- function(object, spec, design, aside, n_boot, refit) {
  spec0 <- object@spec
  design0 <- statmod_design(spec0)
  if (length(attr(design0, "refresh"))) {
    stop(paste0("predictive = \"bootstrap\" is not available for a term ",
                "whose block moves with\n  its coefficients (seg(), jump(), ",
                "jseg(), nl()): a replica would refit\n  its break-points ",
                "or its nonlinear parameters from where the fit left them,\n",
                "  which is a different estimator. predictive = \"averaged\" ",
                "is available."), call. = FALSE)
  }
  y0 <- spec0@response
  if (!is.numeric(y0) || !is.null(dim(y0)) ||
      length(y0) != spec0@n_obs) {
    stop(paste0("predictive = \"bootstrap\" simulates the response, and ",
                "this response is not\n  one number per observation. ",
                "predictive = \"averaged\" is available."), call. = FALSE)
  }
  sm <- object@methods
  is_cv <- function(m) !is.null(m) && identical(m@kind, "cv")
  if (identical(refit, "full") &&
      (is_cv(sm$outer) || is_cv(sm$sparse_criterion))) {
    stop(paste0("boot_refit = \"full\" chooses the hyperparameters again ",
                "on every replica, and\n  a cross-validated criterion needs ",
                "the data's folds, which a fit does not\n  keep. ",
                "boot_refit = \"coefficients\" holds them at the fit's ",
                "values."), call. = FALSE)
  }
  draw <- bootstrap_simulator(object, spec0, design0)
  params <- spec0@distrib@params
  method <- sm$smooth
  cfg <- inner_settings(method)
  # every replica starts at the fit's optimum: its coefficients here, its
  # hyperparameters through object@hyper in bootstrap_refit()
  beta0 <- unlist(object@coefficients[params], use.names = FALSE)
  comps <- list()
  failed <- 0L
  for (b in seq_len(n_boot)) {
    y <- draw()
    spec_b <- S7::set_props(spec0, response = y)
    r <- tryCatch(bootstrap_refit(object, spec_b, refit, method, cfg, beta0),
                  error = function(e) NULL)
    if (is.null(r)) {
      failed <- failed + 1L
      next
    }
    fit_b <- S7::set_props(object, coefficients = r$coefficients,
                           hyper = r$hyper)
    blocks <- random_blocks(spec, design, fit_b, aside)
    C <- predictive_cov(fit_b, spec, design, blocks, fixed = FALSE)$random
    comps[[length(comps) + 1L]] <- list(
      eta = statmod_eta(spec, design, r$coefficients)$eta, C = C)
  }
  if (!length(comps)) {
    stop("No bootstrap replica could be refitted.", call. = FALSE)
  }
  if (failed) {
    warning(sprintf(paste0("%d of %d bootstrap replicas could not be ",
                           "refitted and are left out."), failed, n_boot),
            call. = FALSE)
  }
  comps
}


#' A Simulator of the Response From a Fitted Model
#'
#' @description
#' Returns a function of no arguments that draws one response at the fitting
#' rows: the effects of every [modelterms7::random()] term from the Gaussian
#' prior the fit estimated, every other coefficient at its estimate, and the
#' response from the family at the parameters that gives.
#'
#' @param object The [StatmodFit()].
#' @param spec0,design0 Its specification and design at the fitting rows.
#'
#' @return A function returning a numeric vector of `spec0@n_obs` values.
#'
#' @examples
#' set.seed(3)
#' gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
#' gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40)
#' fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
#'                gg)
#' draw <- statmodels7:::bootstrap_simulator(fit, fit@spec,
#'                                           statmod_design(fit@spec))
#' length(draw())
#'
#' @keywords internal
bootstrap_simulator <- function(object, spec0, design0) {
  rm <- random_modes(spec0, "zero")
  rm <- rm[rm$mode != "conditional", , drop = FALSE]
  bl <- if (nrow(rm)) random_blocks(spec0, design0, object, rm) else list()
  # the square root of each block's covariance, read once: an eigenvalue
  # floor rather than chol(), a covariance at the edge of its chart being
  # singular to the last digit
  roots <- lapply(bl, function(b) {
    ev <- eigen((b$Sigma + t(b$Sigma)) / 2, symmetric = TRUE)
    ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow(b$Sigma))
  })
  cf0 <- object@coefficients
  function() {
    cf <- cf0
    for (j in seq_along(bl)) {
      mem <- bl[[j]]$members
      D <- sum(mem$dim)
      idx1 <- design0[[mem$param[1L]]]$blocks[[mem$key[1L]]]
      m <- length(idx1) %/% mem$dim[1L]
      # one row per group, the members' coordinates side by side
      B <- t(roots[[j]] %*% matrix(stats::rnorm(D * m), D, m))
      at <- cumsum(c(0L, mem$dim))
      for (k in seq_len(nrow(mem))) {
        idx <- design0[[mem$param[k]]]$blocks[[mem$key[k]]]
        # the coefficients are group-major: a group's d coordinates together
        cf[[mem$param[k]]][idx] <- as.vector(
          t(B[, at[k] + seq_len(mem$dim[k]), drop = FALSE]))
      }
    }
    th <- statmod_eta(spec0, design0, cf)$theta
    distributions7::distrib_rng(spec0@distrib, spec0@n_obs, th)
  }
}


#' Refit One Bootstrap Replica
#'
#' @description
#' Re-estimates the coefficients of a fit on a simulated response, at its
#' hyperparameters or with them chosen again. See [predictive_bootstrap()].
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
                          NULL, NULL, sm$sparse_criterion)
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
